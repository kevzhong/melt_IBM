! Calculations of ambient fluid velocity "seen" by the particle
! cf. Kidanemariam et al. (2013)
! Use to estimate the slip (relative) velocity between the particle and ambient fluid

! Consistent with Kidanemariam et al. (2013), we take the velocity "seen" by the particle to be
! the velocity sampled on a spherical shell of radiius R_shell, centered at the COM of the particle

! Slip velocity is just the relative motion of this velocity wrt. the particle centroid
! u_slip = u_S - U_centroid,    where u_S is the velocity sampled on the spherical shell 

! For melting objects, since we don't necessarily have an object radius R_sphere, we can take the equivalent radius instead


! Also a primitive measure included: the relative motion in the entire computaitonal domain (box-averaged)
! Given by calcRelBoxVel()

      subroutine calcRelBoxVel
        use mpih
        use param
        use local_arrays,only: vx,vy,vz
        use stat_arrays
        use mls_param
        use mpi_param, only: kstart,kend
  
        implicit none
        integer :: ic,jc,kc
        integer :: counter
        real :: Urel_box, Vrel_box, Wrel_box
        character(70) namfile
        

        counter = 0
        Urel_box = 0.0
        Vrel_box = 0.0
        Wrel_box = 0.0

  
         do kc=kstart,kend
         do jc=1,n2m
         do ic=1,n1m  

            if ( solid_mask(ic,jc,kc) .eqv. .false. ) then
            ! Box-averaged relative velocity
                counter = counter + 1

                Urel_box = Urel_box + ( vx(ic,jc,kc) - vel_CM(1,1) )
                Vrel_box = Vrel_box + ( vy(ic,jc,kc) - vel_CM(2,1) )
                Wrel_box = Wrel_box + ( vz(ic,jc,kc) - vel_CM(3,1) )
            endif

         end do
         end do
         end do
             
        call MpiAllSumRealScalar(Urel_box)
        call MpiAllSumRealScalar(Vrel_box)
        call MpiAllSumRealScalar(Wrel_box)
        call MPI_ALLREDUCE(MPI_IN_PLACE,counter,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        Urel_box = Urel_box / dble(counter)
        Vrel_box = Vrel_box / dble(counter)
        Wrel_box = Wrel_box / dble(counter)

        ! kenerg =0.5d0* kenerg / (dble(n1m*n2m*n3m))
        ! nu=1.0d0/ren
        ! diss_volAvg = 2.0d0 * nu*diss_volAvg / (dble(n1m*n2m*n3m))  
        ! urms = sqrt( 2.0 / 3.0 * kenerg )
        ! eta = (nu**3/diss_volAvg)**(0.25)
        ! kmax_eta = pi * dx1 * eta ! kmax_eta where kmax = pi /  dx is the Nyquist limit
        ! Re_L = kenerg**2 / (nu * diss_volAvg)
        ! lambda_t = sqrt( 15.0 * nu * urms**2 / diss_volAvg )
        ! re_lam = urms * lambda_t / nu
        ! L_int = sqrt(kenerg**3) / diss_volAvg ! Large-eddy length scale

        if(ismaster) then
        namfile='stringdata/rel_velBox.txt'
        open(unit=92,file=namfile, Access='append', Status='unknown')
        write(92,'(100E15.7)') time, Urel_box, Vrel_box, Wrel_box
        close(92)
        end if
     
        return
    end subroutine calcRelBoxVel


    subroutine calcRelShellVel
        use mpih
        use param
        use local_arrays,only: vx,vy,vz
        use stat_arrays
        use mls_param
        use mpi_param, only: kstart,kend
  
        implicit none
        integer :: i,j,k,ii,jj,kk
        real :: R_shell, R_eqv
        real :: rr_x, rr_y, rr_z, phi
        integer, dimension(3,2) :: bbox_inds
        real, dimension(3,2) :: lim
        real, dimension(3) :: x_GC
        character(70) namfile

        integer :: nshell
        real, dimension(6) :: Rshell_on_Reqv
        real, dimension(6) :: sumPhix, sumPhiy, sumPhiz, Ushell, Vshell, Wshell

        ! Varying choices of shell radius to test
        ! Up to 6 equivalent radii (= 3 equivalent diameters, cf. Kidanemariam et al. 2013)
        Rshell_on_Reqv = [6.0, 5.0, 4.0, 3.0, 2.0, 1.5]

        

        ! Running averages of the shell kernel functions
        sumPhix = 0.0
        sumPhiy = 0.0
        sumPhiz = 0.0

        ! Shell-averaged velocities
        Ushell = 0.0
        Vshell = 0.0
        Wshell = 0.0

        ! First calculate equivalent radius based on current object volume
        R_eqv =  (3.0 * Volume(1) / (4.0 * pi))**(1.0/3.0)

        !----------------- Shell radius  ---------------------------------------
        ! Compute the largest shell radius: use to construct the bounding box to loop over
        R_shell = Rshell_on_Reqv(1) * R_eqv

        !-----------------------------------------------------------------------

          ! get bounding box
        do i = 1,3
            lim(i,1) =  pos_CM(i,1) - R_shell ! min
            lim(i,2) =  pos_CM(i,1) + R_shell ! max
        end do
        bbox_inds = floor(lim*dx1) + 1 ! compute indices cell centered

        do i = bbox_inds(1,1),bbox_inds(1,2)
            do j = bbox_inds(2,1),bbox_inds(2,2)
              do kk = bbox_inds(3,1),bbox_inds(3,2)

                k = kk
                call get_periodic_indices(k,x_GC)

                if (k.ge.kstart.and.k.le.kend) then

                    ii = modulo(i-1,n1m) + 1
                    jj = modulo(j-1,n2m) + 1

                    ! Abs distances to centroid from the (i,j,k) cell
                    rr_x =  norm2 (  [ xc(i), ym(j), zm(kk) ]  - pos_CM(:,1)  ) 
                    rr_y =  norm2 (  [ xm(i), yc(j), zm(kk) ]  - pos_CM(:,1)  ) 
                    rr_z =  norm2 (  [ xm(i), ym(j), zc(kk) ]  - pos_CM(:,1)  )

                    do nshell = 1,6 ! loop over each shell

                      R_shell = Rshell_on_Reqv(nshell) * R_eqv

                      ! x component
                      phi = 1.0 / ( cosh( (rr_x - R_shell) * dx1 ) )
                      sumPhix(nshell) = sumPhix(nshell) + phi
                      Ushell(nshell) = Ushell(nshell) + phi * vx(ii,jj,k)

                      ! y component
                      phi = 1.0 / ( cosh( (rr_y - R_shell) * dx1 ) )
                      sumPhiy(nshell) = sumPhiy(nshell) + phi
                      Vshell(nshell) = Vshell(nshell) + phi * vy(ii,jj,k)

                      ! z component
                      phi = 1.0 / ( cosh( (rr_z - R_shell) * dx1 ) )
                      sumPhiz(nshell) = sumPhiz(nshell) + phi
                      Wshell(nshell) = Wshell(nshell) + phi * vz(ii,jj,k)

                    enddo

                endif
              enddo
            enddo
        enddo
             
        !call MpiAllSumRealScalar(Ushell)
        !call MpiAllSumRealScalar(Vshell)
        !call MpiAllSumRealScalar(Wshell)
        !call MpiAllSumRealScalar(sumPhix)
        !call MpiAllSumRealScalar(sumPhiy)
        !call MpiAllSumRealScalar(sumPhiz)

        call MpiSumReal1D(Ushell, 6)
        call MpiSumReal1D(Vshell, 6)
        call MpiSumReal1D(Wshell, 6)
        call MpiSumReal1D(sumPhix, 6)
        call MpiSumReal1D(sumPhiy, 6)
        call MpiSumReal1D(sumPhiz, 6)

        do i = 1,6
          Ushell(i) = Ushell(i) / sumPhix(i)
          Vshell(i) = Vshell(i) / sumPhiy(i)
          Wshell(i) = Wshell(i) / sumPhiz(i)
        enddo

        if(ismaster) then

        ! x-component
          namfile='stringdata/rel_UShell.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, Ushell
          close(92)

          ! y-component
          namfile='stringdata/rel_VShell.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, Vshell
          close(92)

          ! z-component
          namfile='stringdata/rel_WShell.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, Wshell
          close(92)

        end if
     
        return
    end subroutine calcRelShellVel

  
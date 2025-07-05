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

      subroutine calcFluidVelAvgs
        ! Instantaneous fluid-domain-avergaed velocity components
        ! For expressions like (3), (4) in Bellani & Variano (2012)
        use mpih
        use param
        use local_arrays,only: vx,vy,vz,temp
        use stat_arrays
        use mls_param
        use mpi_param, only: kstart,kend
  
        implicit none
        integer :: ic,jc,kc
        integer :: counter
        real :: ufluid, vfluid, wfluid
        real :: Tfluid, DT_fluid, DT2_fluid
        real :: ufluid_var, vfluid_var, wfluid_var
        real :: Urel_box, Vrel_box, Wrel_box
        real :: uu_rel, vv_rel, ww_rel
        integer :: ip, im, jp, jm, kp, km

        character(70) namfile
        

        counter = 0

        ufluid = 0.0
        vfluid = 0.0
        wfluid = 0.0
        
        Tfluid = 0.0
        DT_fluid = 0.0
        DT2_fluid = 0.0


        ufluid_var = 0.0
        vfluid_var = 0.0
        wfluid_var = 0.0

        Urel_box = 0.0
        Vrel_box = 0.0
        Wrel_box = 0.0

        uu_rel = 0.0
        vv_rel = 0.0
        ww_rel = 0.0

        ! For scalar dissipation rate
        call update_both_ghosts(n1,n2,temp,kstart,kend)



         ! First, calculate instantaneous box-averages of instantaneous velocity

         do kc=kstart,kend
          kp=kc+1
          km=kc-1
         do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
         do ic=1,n1m  
          im= imv(ic)
          ip= ipv(ic)

            ! ignoring grid-stagger, should be mostly negligible since no. of interface cells is not too significant
            if ( VOFp(ic,jc,kc) .eq. 1.0 ) then

                counter = counter + 1

                ! Fluid-domain averaged velocities
                ufluid = ufluid +  vx(ic,jc,kc) 
                vfluid = vfluid +  vy(ic,jc,kc) 
                wfluid = wfluid +  vz(ic,jc,kc) 

                ! First part of (3) in B&V(2012)
                ! Component-wise fluid velocity variance
                ufluid_var = ufluid_var + ( vx(ic,jc,kc) )**2
                vfluid_var = vfluid_var + ( vy(ic,jc,kc) )**2
                wfluid_var = wfluid_var + ( vz(ic,jc,kc) )**2

                ! Relative velocities
                Urel_box = Urel_box + ( vx(ic,jc,kc) - vel_CM(1,1) )
                Vrel_box = Vrel_box + ( vy(ic,jc,kc) - vel_CM(2,1) )
                Wrel_box = Wrel_box + ( vz(ic,jc,kc) - vel_CM(3,1) )


                ! Variance of relative velocity
                uu_rel = uu_rel + ( vx(ic,jc,kc) - vel_CM(1,1) )**2
                vv_rel = vv_rel + ( vy(ic,jc,kc) - vel_CM(2,1) )**2
                ww_rel = ww_rel + ( vz(ic,jc,kc) - vel_CM(3,1) )**2


                ! Temperature
                Tfluid = Tfluid + temp(ic,jc,kc)
                DT_fluid = DT_fluid + ( temp(ic,jc,kc) - Tliq )
                DT2_fluid = DT2_fluid + ( temp(ic,jc,kc) - Tliq )**2

            endif

         end do
         end do
         end do
             
        ! Number of fluid cells
        call MPI_ALLREDUCE(MPI_IN_PLACE,counter,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)


        call MpiAllSumRealScalar(ufluid)
        call MpiAllSumRealScalar(vfluid)
        call MpiAllSumRealScalar(wfluid)
        ufluid = ufluid / dble(counter)
        vfluid = vfluid / dble(counter)
        wfluid = wfluid / dble(counter)

        call MpiAllSumRealScalar(ufluid_var)
        call MpiAllSumRealScalar(vfluid_var)
        call MpiAllSumRealScalar(wfluid_var)    
        ufluid_var = ufluid_var / dble(counter)
        vfluid_var = vfluid_var / dble(counter)
        wfluid_var = wfluid_var / dble(counter)       

        call MpiAllSumRealScalar(Urel_box)
        call MpiAllSumRealScalar(Vrel_box)
        call MpiAllSumRealScalar(Wrel_box)
        Urel_box = Urel_box / dble(counter)
        Vrel_box = Vrel_box / dble(counter)
        Wrel_box = Wrel_box / dble(counter)

        call MpiAllSumRealScalar(uu_rel)
        call MpiAllSumRealScalar(vv_rel)
        call MpiAllSumRealScalar(ww_rel)
        uu_rel = uu_rel / dble(counter)
        vv_rel = vv_rel / dble(counter)
        ww_rel = ww_rel / dble(counter)

        call MpiAllSumRealScalar(Tfluid)
        call MpiAllSumRealScalar(DT_fluid)
        call MpiAllSumRealScalar(DT2_fluid)
        Tfluid = Tfluid / dble(counter)
        DT_fluid = DT_fluid / dble(counter)
        DT2_fluid = DT2_fluid / dble(counter)


        if(ismaster) then
          namfile='stringdata/fluidAvgVels.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(E15.7, I12, 6E15.7)') time, counter, ufluid, vfluid, wfluid, ufluid_var, vfluid_var, wfluid_var
          close(92)
        end if

        if(ismaster) then
          namfile='stringdata/rel_velBox.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(E15.7, I12, 6E15.7)') time, counter, Urel_box, Vrel_box, Wrel_box, uu_rel, vv_rel, ww_rel
          close(92)
        end if

        if(ismaster) then
          namfile='stringdata/fluidAvgTemps.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(E15.7, I12, 5E15.7)') time, counter, Tfluid, DT_fluid, DT2_fluid
          close(92)
        end if

        return
    end subroutine calcFluidVelAvgs

  subroutine calcLocalShellFlow
        use mpih
        use param
        use local_arrays,only: vx,vy,vz,temp
        use local_aux,only: vorx, vory, vorz, diss, tke, chi
        use stat_arrays
        use mls_param
        use mpi_param, only: kstart,kend
  
        implicit none
        integer :: i,j,k,ii,jj,kk
        integer :: iip, jjp
        real :: R_eqv
        real :: rr_t
        integer, dimension(3,2) :: bbox_inds
        real, dimension(3,2) :: lim
        real, dimension(3) :: x_GC
        character(70) namfile
        real :: Rshell_min, Rshell_max
        real :: u_interp, v_interp, w_interp, urel, vrel, wrel

        real :: sumVOFmin, sumVOFmax
        real :: sumU_min    , sumU_max    
        real :: sumV_min    , sumV_max    
        real :: sumW_min    , sumW_max    
        real :: sumUrel_min , sumUrel_max 
        real :: sumVrel_min , sumVrel_max 
        real :: sumWrel_min , sumWrel_max 
        real :: sumTKE_min  , sumTKE_max  
        real :: sumDISS_min , sumDISS_max 
        real :: sumTemp_min , sumTemp_max 
        real :: sumCHImin   , sumCHImax   
        real :: sumUTmin    , sumUTmax    
        real :: sumVTmin    , sumVTmax    
        real :: sumWTmin    , sumWTmax    
        real :: sumUTrelmin , sumUTrelmax 
        real :: sumVTrelmin , sumVTrelmax 
        real :: sumWTrelmin , sumWTrelmax 


        ! Quantities to sample:
        ! 1) <u - Uc>
        ! 2) <v - Vc>
        ! 3) <w - Wc>
        ! 4) tke
        ! 5) dissipation
        ! 6) Temperature
        ! 7) chi
        ! 8) < (u-Uc)T >
        ! 9) < (v - Vc) > T
        ! 10) < (w - Wc) > T

        ! Boundary layer sample: <= 1.5 R_eqv
        ! Local sample: 1.5 to 3R_eqv

        ! Running averages of the shell kernel functions
        sumVOFmin = 0.0
        sumVOFmax = 0.0
 
        ! Shell-averaged quantities
        sumU_min    = 0.0   ; sumU_max    = 0.0
        sumV_min    = 0.0   ; sumV_max    = 0.0
        sumW_min    = 0.0   ; sumW_max    = 0.0
        sumUrel_min = 0.0   ; sumUrel_max = 0.0
        sumVrel_min = 0.0   ; sumVrel_max = 0.0
        sumWrel_min = 0.0   ; sumWrel_max = 0.0
        sumTKE_min  = 0.0   ; sumTKE_max  = 0.0
        sumDISS_min = 0.0   ; sumDISS_max = 0.0
        sumTemp_min = 0.0   ; sumTemp_max = 0.0
        sumCHImin   = 0.0   ; sumCHImax   = 0.0
        sumUTmin    = 0.0   ; sumUTmax    = 0.0
        sumVTmin    = 0.0   ; sumVTmax    = 0.0
        sumWTmin    = 0.0   ; sumWTmax    = 0.0
        sumUTrelmin = 0.0   ; sumUTrelmax = 0.0
        sumVTrelmin = 0.0   ; sumVTrelmax = 0.0
        sumWTrelmin = 0.0   ; sumWTrelmax = 0.0

        ! First calculate equivalent radius based on current object volume
        R_eqv =  (3.0 * Volume(1) / (4.0 * pi))**(1.0/3.0)

        !----------------- Shell radii  ---------------------------------------
        ! 
        Rshell_min = 1.5 * R_eqv
        Rshell_max = 3.0 * R_eqv

        !-----------------------------------------------------------------------

        ! ----------- GET BOUNDING BOX AND ITS INDICES ------------------------
        do i = 1,3
            lim(i,1) =  pos_CM(i,1) - Rshell_max ! min
            lim(i,2) =  pos_CM(i,1) + Rshell_max ! max
        end do
        bbox_inds = floor(lim*dx1) + 1 ! compute indices cell centered
        ! ----------------------------------------------------------------------

        do i = bbox_inds(1,1),bbox_inds(1,2)
            do j = bbox_inds(2,1),bbox_inds(2,2)
              do kk = bbox_inds(3,1),bbox_inds(3,2)

                k = kk
                call get_periodic_indices(k,x_GC)

                if (k.ge.kstart.and.k.le.kend) then

                    ii = modulo(i-1,n1m) + 1
                    jj = modulo(j-1,n2m) + 1

                    iip = modulo(i,n1m) + 1
                    jjp = modulo(j,n2m) + 1

                    ! Abs distances to centroid from the (i,j,k) cell
                    rr_t =  norm2 (  [ xm(i), ym(j), zm(kk) ]  - pos_CM(:,1)  )

                    if (VOFp(ii,jj,k) .eq. 1.0 ) then ! Only accumulate avg. in fluid domain

                      u_interp = 0.5 * ( vx(ii,jj,k) + vx(iip,jj ,k) )
                      v_interp = 0.5 * ( vy(ii,jj,k) + vx(ii ,jjp,k) )
                      w_interp = 0.5 * ( vz(ii,jj,k) + vz(ii ,jj ,k+1) )
                      urel = u_interp - vel_CM(1,1)
                      vrel = v_interp - vel_CM(2,1)
                      wrel = w_interp - vel_CM(3,1)

                      ! SMALL SHELL: <= 1.5R_EQV
                      if (rr_t .le. Rshell_min) then ! r_ijk < 1.5 R_Eqv
                        sumVOFmin = sumVOFmin + VOFp(ii,jj,k)

                        ! Velocities
                        sumU_min = sumU_min + u_interp
                        sumV_min = sumV_min + v_interp
                        sumW_min = sumW_min + w_interp
                        ! Relative
                        sumUrel_min = sumUrel_min + urel
                        sumVrel_min = sumVrel_min + vrel
                        sumWrel_min = sumWrel_min + wrel
                        ! Turbulence quantities
                        sumTKE_min = sumTKE_min + tke(ii,jj,k)
                        sumDISS_min = sumDISS_min + diss(ii,jj,k)
                        ! Temperature
                        sumTemp_min = sumTemp_min + temp(ii,jj,k)
                        sumCHImin = sumCHImin + chi(ii,jj,k)
                        ! Temperature transport fluxes
                        sumUTmin = sumUTmin + u_interp * (temp(ii,jj,k) - Tliq)
                        sumVTmin = sumVTmin + v_interp * (temp(ii,jj,k) - Tliq)
                        sumWTmin = sumWTmin + w_interp * (temp(ii,jj,k) - Tliq)

                        ! Relative
                        sumUTrelmin = sumUTrelmin + urel * (temp(ii,jj,k) - Tliq)
                        sumVTrelmin = sumVTrelmin + vrel * (temp(ii,jj,k) - Tliq)
                        sumWTrelmin = sumWTrelmin + wrel * (temp(ii,jj,k) - Tliq)
                      endif

                      if ( (rr_t .gt. Rshell_min) .and. (rr_t .le. Rshell_max) ) then
                        sumVOFmax = sumVOFmax + VOFp(ii,jj,k)

                        ! Velocities
                        sumU_max = sumU_max + u_interp
                        sumV_max = sumV_max + v_interp
                        sumW_max = sumW_max + w_interp
                        ! Relative
                        sumUrel_max = sumUrel_max + urel
                        sumVrel_max = sumVrel_max + vrel
                        sumWrel_max = sumWrel_max + wrel
                        ! Turbulence quantities
                        sumTKE_max = sumTKE_max + tke(ii,jj,k)
                        sumDISS_max = sumDISS_max + diss(ii,jj,k)
                        ! Temperature
                        sumTemp_max = sumTemp_max + temp(ii,jj,k)
                        sumCHImax = sumCHImax + chi(ii,jj,k)
                        ! Temperature transport fluxes
                        sumUTmax = sumUTmax + u_interp * (temp(ii,jj,k) - Tliq)
                        sumVTmax = sumVTmax + v_interp * (temp(ii,jj,k) - Tliq)
                        sumWTmax = sumWTmax + w_interp * (temp(ii,jj,k) - Tliq)

                        ! Relative
                        sumUTrelmax = sumUTrelmax + urel * (temp(ii,jj,k) - Tliq)
                        sumVTrelmax = sumVTrelmax + vrel * (temp(ii,jj,k) - Tliq)
                        sumWTrelmax = sumWTrelmax + wrel * (temp(ii,jj,k) - Tliq)

                      endif

                    endif

                endif
              enddo
            enddo
        enddo

        call MpiAllSumRealScalar(sumVOFmin)
        call MpiAllSumRealScalar(sumU_min)
        call MpiAllSumRealScalar(sumV_min)
        call MpiAllSumRealScalar(sumW_min)
        call MpiAllSumRealScalar(sumUrel_min)
        call MpiAllSumRealScalar(sumVrel_min)
        call MpiAllSumRealScalar(sumWrel_min)
        call MpiAllSumRealScalar(sumTKE_min)
        call MpiAllSumRealScalar(sumDISS_min)
        call MpiAllSumRealScalar(sumTemp_min)
        call MpiAllSumRealScalar(sumCHImin)
        call MpiAllSumRealScalar(sumUTmin)
        call MpiAllSumRealScalar(sumVTmin)
        call MpiAllSumRealScalar(sumWTmin)
        call MpiAllSumRealScalar(sumUTrelmin)
        call MpiAllSumRealScalar(sumVTrelmin)
        call MpiAllSumRealScalar(sumWTrelmin)
        sumU_min  =   sumU_min    / sumVOFmin
        sumV_min  =   sumV_min    / sumVOFmin
        sumW_min  =   sumW_min    / sumVOFmin
        sumUrel_min = sumUrel_min / sumVOFmin
        sumVrel_min = sumVrel_min / sumVOFmin
        sumWrel_min = sumWrel_min / sumVOFmin
        sumTKE_min  = sumTKE_min  / sumVOFmin
        sumDISS_min = sumDISS_min / sumVOFmin
        sumTemp_min = sumTemp_min / sumVOFmin
        sumCHImin   = sumCHImin   / sumVOFmin
        sumUTmin    = sumUTmin    / sumVOFmin
        sumVTmin    = sumVTmin    / sumVOFmin
        sumWTmin    = sumWTmin    / sumVOFmin
        sumUTrelmin = sumUTrelmin / sumVOFmin
        sumVTrelmin = sumVTrelmin / sumVOFmin
        sumWTrelmin = sumWTrelmin / sumVOFmin

        
        call MpiAllSumRealScalar(sumVOFmax)
        call MpiAllSumRealScalar(sumU_max)
        call MpiAllSumRealScalar(sumV_max)
        call MpiAllSumRealScalar(sumW_max)
        call MpiAllSumRealScalar(sumUrel_max)
        call MpiAllSumRealScalar(sumVrel_max)
        call MpiAllSumRealScalar(sumWrel_max)
        call MpiAllSumRealScalar(sumTKE_max)
        call MpiAllSumRealScalar(sumDISS_max)
        call MpiAllSumRealScalar(sumTemp_max)
        call MpiAllSumRealScalar(sumCHImax)
        call MpiAllSumRealScalar(sumUTmax)
        call MpiAllSumRealScalar(sumVTmax)
        call MpiAllSumRealScalar(sumWTmax)
        call MpiAllSumRealScalar(sumUTrelmax)
        call MpiAllSumRealScalar(sumVTrelmax)
        call MpiAllSumRealScalar(sumWTrelmax)
        sumU_max  =   sumU_max    / sumVOFmax
        sumV_max  =   sumV_max    / sumVOFmax
        sumW_max  =   sumW_max    / sumVOFmax
        sumUrel_max = sumUrel_max / sumVOFmax
        sumVrel_max = sumVrel_max / sumVOFmax
        sumWrel_max = sumWrel_max / sumVOFmax
        sumTKE_max  = sumTKE_max  / sumVOFmax
        sumDISS_max = sumDISS_max / sumVOFmax
        sumTemp_max = sumTemp_max / sumVOFmax
        sumCHImax   = sumCHImax   / sumVOFmax
        sumUTmax    = sumUTmax    / sumVOFmax
        sumVTmax    = sumVTmax    / sumVOFmax
        sumWTmax    = sumWTmax    / sumVOFmax
        sumUTrelmax = sumUTrelmax / sumVOFmax
        sumVTrelmax = sumVTrelmax / sumVOFmax
        sumWTrelmax = sumWTrelmax / sumVOFmax

        if(ismaster) then
          !------------------- MINIMUM SHELL -----------------------------------------------
          ! Momentum field
          namfile='stringdata/minShell_mmtm.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumU_min, sumV_min, sumW_min, sumTKE_min, sumDISS_min
          close(92)

          ! Momentum field with relative motion
          namfile='stringdata/minShell_Relmmtm.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumUrel_min, sumVrel_min, sumWrel_min
          close(92)

          ! Temperature field
          namfile='stringdata/minShell_temp.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumTemp_min, sumCHImin, sumUTmin, sumVTmin, sumWTmin
          close(92)

          ! Relative temperature field
          namfile='stringdata/minShell_Reltemp.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumUTrelmin, sumVTrelmin, sumWTrelmin
          close(92)

          !------------------- MAXIMUM SHELL -----------------------------------------------
          ! Momentum field
          namfile='stringdata/maxShell_mmtm.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumU_max, sumV_max, sumW_max, sumTKE_max, sumDISS_max
          close(92)

          ! Momentum field with relative motion
          namfile='stringdata/maxShell_Relmmtm.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumUrel_max, sumVrel_max, sumWrel_max
          close(92)

          ! Temperature field
          namfile='stringdata/maxShell_temp.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumTemp_max, sumCHImax, sumUTmax, sumVTmax, sumWTmax
          close(92)

          ! Relative temperature field
          namfile='stringdata/maxShell_Reltemp.txt'
          open(unit=92,file=namfile, Access='append', Status='unknown')
          write(92,'(100E15.7)') time, sumUTrelmax, sumVTrelmax, sumWTrelmax
          close(92)

        end if
     
        return
    end subroutine calcLocalShellFlow

    ! subroutine calcRelShellVel
    !     use mpih
    !     use param
    !     use local_arrays,only: vx,vy,vz,temp
    !     use stat_arrays
    !     use mls_param
    !     use mpi_param, only: kstart,kend
  
    !     implicit none
    !     integer :: i,j,k,ii,jj,kk
    !     integer :: iip, jjp
    !     real :: R_shell, R_eqv
    !     real :: rr_x, rr_y, rr_z, rr_t, phi
    !     integer, dimension(3,2) :: bbox_inds
    !     real, dimension(3,2) :: lim
    !     real, dimension(3) :: x_GC
    !     character(70) namfile

    !     integer :: nshell
    !     real, dimension(6) :: Rshell_on_Reqv
    !     real, dimension(6) :: sumPhix, sumPhiy, sumPhiz, sumPhiT, Ushell, Vshell, Wshell, Tshell
    !     real, dimension(6) :: uTshell, vTshell, wTshell
    !     real :: u_interp, v_interp, w_interp

    !     ! Varying choices of shell radius to test
    !     ! Up to 6 equivalent radii (= 3 equivalent diameters, cf. Kidanemariam et al. 2013)
    !     Rshell_on_Reqv = [6.0, 5.0, 4.0, 3.0, 2.0, 1.5]

        

    !     ! Running averages of the shell kernel functions
    !     sumPhix = 0.0
    !     sumPhiy = 0.0
    !     sumPhiz = 0.0
    !     sumPhiT = 0.0

    !     ! Shell-averaged velocities
    !     Ushell = 0.0
    !     Vshell = 0.0
    !     Wshell = 0.0
    !     Tshell = 0.0

    !     ! First calculate equivalent radius based on current object volume
    !     R_eqv =  (3.0 * Volume(1) / (4.0 * pi))**(1.0/3.0)

    !     !----------------- Shell radius  ---------------------------------------
    !     ! Compute the largest shell radius: use to construct the bounding box to loop over
    !     R_shell = Rshell_on_Reqv(1) * R_eqv

    !     !-----------------------------------------------------------------------

    !       ! get bounding box
    !     do i = 1,3
    !         lim(i,1) =  pos_CM(i,1) - R_shell ! min
    !         lim(i,2) =  pos_CM(i,1) + R_shell ! max
    !     end do
    !     bbox_inds = floor(lim*dx1) + 1 ! compute indices cell centered

    !     do i = bbox_inds(1,1),bbox_inds(1,2)
    !         do j = bbox_inds(2,1),bbox_inds(2,2)
    !           do kk = bbox_inds(3,1),bbox_inds(3,2)

    !             k = kk
    !             call get_periodic_indices(k,x_GC)

    !             if (k.ge.kstart.and.k.le.kend) then

    !                 ii = modulo(i-1,n1m) + 1
    !                 jj = modulo(j-1,n2m) + 1

    !                 iip = modulo(i,n1m) + 1
    !                 jjp = modulo(j,n2m) + 1

    !                 ! Abs distances to centroid from the (i,j,k) cell
    !                 rr_x =  norm2 (  [ xc(i), ym(j), zm(kk) ]  - pos_CM(:,1)  ) 
    !                 rr_y =  norm2 (  [ xm(i), yc(j), zm(kk) ]  - pos_CM(:,1)  ) 
    !                 rr_z =  norm2 (  [ xm(i), ym(j), zc(kk) ]  - pos_CM(:,1)  )
    !                 rr_t =  norm2 (  [ xm(i), ym(j), zm(kk) ]  - pos_CM(:,1)  )

    !                 do nshell = 1,6 ! loop over each shell

    !                   R_shell = Rshell_on_Reqv(nshell) * R_eqv

    !                   ! x component
    !                   phi = 1.0 / ( cosh( 0.5*(rr_x - R_shell) * dx1 ) )**2
    !                   sumPhix(nshell) = sumPhix(nshell) + phi
    !                   Ushell(nshell) = Ushell(nshell) + phi * vx(ii,jj,k)

    !                   ! y component
    !                   phi = 1.0 / ( cosh( 0.5*(rr_y - R_shell) * dx1 ) )**2
    !                   sumPhiy(nshell) = sumPhiy(nshell) + phi
    !                   Vshell(nshell) = Vshell(nshell) + phi * vy(ii,jj,k)

    !                   ! z component
    !                   phi = 1.0 / ( cosh( 0.5*(rr_z - R_shell) * dx1 ) )**2
    !                   sumPhiz(nshell) = sumPhiz(nshell) + phi
    !                   Wshell(nshell) = Wshell(nshell) + phi * vz(ii,jj,k)


    !                   ! temperature component
    !                   phi = 1.0 / ( cosh( 0.5*(rr_t - R_shell) * dx1 ) )**2
    !                   sumPhiT(nshell) = sumPhiT(nshell) + phi
    !                   Tshell(nshell) = Tshell(nshell) + phi * temp(ii,jj,k)


    !                   ! Temperature transport fluxes
    !                   ! phi from cell-centered temperature is re-used
    !                   ! temperature relative to Tamb
    !                   u_interp = 0.5 * ( vx(ii,jj,k) + vx(iip,jj ,k) )
    !                   v_interp = 0.5 * ( vy(ii,jj,k) + vx(ii ,jjp,k) )
    !                   w_interp = 0.5 * ( vz(ii,jj,k) + vz(ii ,jj ,k+1) )

    !                   uTshell(nshell) = uTshell(nshell) + phi * u_interp * (temp(ii,jj,k) - Tliq)
    !                   vTshell(nshell) = vTshell(nshell) + phi * v_interp * (temp(ii,jj,k) - Tliq)
    !                   wTshell(nshell) = wTshell(nshell) + phi * w_interp * (temp(ii,jj,k) - Tliq)
    !                 enddo

    !             endif
    !           enddo
    !         enddo
    !     enddo
             
    !     !call MpiAllSumRealScalar(Ushell)
    !     !call MpiAllSumRealScalar(Vshell)
    !     !call MpiAllSumRealScalar(Wshell)
    !     !call MpiAllSumRealScalar(sumPhix)
    !     !call MpiAllSumRealScalar(sumPhiy)
    !     !call MpiAllSumRealScalar(sumPhiz)

    !     call MpiSumReal1D(Ushell, 6)
    !     call MpiSumReal1D(Vshell, 6)
    !     call MpiSumReal1D(Wshell, 6)
    !     call MpiSumReal1D(Tshell, 6)
    !     call MpiSumReal1D(uTshell, 6)
    !     call MpiSumReal1D(vTshell, 6)
    !     call MpiSumReal1D(wTshell, 6)

    !     call MpiSumReal1D(sumPhix, 6)
    !     call MpiSumReal1D(sumPhiy, 6)
    !     call MpiSumReal1D(sumPhiz, 6)
    !     call MpiSumReal1D(sumPhiT, 6)

    !     do i = 1,6
    !       Ushell(i) = Ushell(i) / sumPhix(i)
    !       Vshell(i) = Vshell(i) / sumPhiy(i)
    !       Wshell(i) = Wshell(i) / sumPhiz(i)
    !       Tshell(i) = Tshell(i) / sumPhiT(i)
    !       uTshell(i) = uTshell(i) / sumPhiT(i)
    !       vTshell(i) = vTshell(i) / sumPhiT(i)
    !       wTshell(i) = wTshell(i) / sumPhiT(i)
    !     enddo

    !     if(ismaster) then

    !     ! x-component
    !       namfile='stringdata/rel_UShell.txt'
    !       open(unit=92,file=namfile, Access='append', Status='unknown')
    !       write(92,'(100E15.7)') time, Ushell
    !       close(92)

    !       ! y-component
    !       namfile='stringdata/rel_VShell.txt'
    !       open(unit=92,file=namfile, Access='append', Status='unknown')
    !       write(92,'(100E15.7)') time, Vshell
    !       close(92)

    !       ! z-component
    !       namfile='stringdata/rel_WShell.txt'
    !       open(unit=92,file=namfile, Access='append', Status='unknown')
    !       write(92,'(100E15.7)') time, Wshell
    !       close(92)

    !       ! temperature-component
    !       namfile='stringdata/rel_TShell.txt'
    !       open(unit=92,file=namfile, Access='append', Status='unknown')
    !       write(92,'(100E15.7)') time, Tshell
    !       close(92)

    !       ! temperature fluxes
    !       namfile='stringdata/rel_uTShell.txt'
    !       open(unit=92,file=namfile, Access='append', Status='unknown')
    !       write(92,'(100E15.7)') time, uTshell
    !       close(92)

    !       namfile='stringdata/rel_vTShell.txt'
    !       open(unit=92,file=namfile, Access='append', Status='unknown')
    !       write(92,'(100E15.7)') time, vTshell
    !       close(92)

    !       namfile='stringdata/rel_wTShell.txt'
    !       open(unit=92,file=namfile, Access='append', Status='unknown')
    !       write(92,'(100E15.7)') time, wTshell
    !       close(92)
    !     end if
     
    !     return
    ! end subroutine calcRelShellVel

  
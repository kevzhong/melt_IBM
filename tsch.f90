      subroutine tschem 
      use param
      use local_arrays
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_aux
      use stat_arrays, only: vxvyvz_rms_vol
      implicit none
      real,dimension(3,3)     :: AA, AAT
      real,dimension(3,2)     :: bbox_inds
      real,dimension(3,Nparticle) :: vel_m1,pos_m1,pos_k,om_m1
      real :: tstart, tend

      integer :: ns, inp, ntr, nsub
      integer :: i,j,k
      

      beta=dt/ren*0.5d0

      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        wtime_vof = 0.
        !eul_solve_wtime = 0.
        !mls_wtime = 0.
        !pressure_wtime = 0.
        !hit_wtime = 0.

        tstart = MPI_WTIME()
        
        VOFx(:,:,:) = 1.
        VOFy(:,:,:) = 1.
        VOFz(:,:,:) = 1.
        VOFp(:,:,:) = 1.
        if(imlsfor.eq.1)then
          do inp=1,Nparticle
            !call calc_rot_matrix(quat(:,inp),AA)

            call get_bbox_inds(bbox_inds,inp)

            call convex_hull_q12(bbox_inds,inp)
            call convex_hull_q22(bbox_inds,inp)
            call convex_hull_q32(bbox_inds,inp)
            call convex_hull_qc2(bbox_inds,inp)
        enddo
        endif
        
        tend = MPI_WTIME()
        wtime_vof = wtime_vof + (tend - tstart)

        tstart = MPI_WTIME()
        call hdnl1
        call hdnl2
        call hdnl3
        call hdnlte


        call invtr1 
        call invtr2      
        call invtr3
        call invtrte
        tend = MPI_WTIME()
        eul_solve_wtime = eul_solve_wtime + (tend - tstart)


        VOFx(:,:,:) = 1.
        VOFy(:,:,:) = 1.
        VOFz(:,:,:) = 1.
        VOFp(:,:,:) = 1.

        tstart = MPI_WTIME()
        call particle
        tend = MPI_WTIME()
        mls_wtime = mls_wtime + (tend - tstart)

        tstart = MPI_WTIME()

        !!------------ KZ: update VOF, remove later since called by Newton--Euler -----
        !if ( (imlsfor.eq.1) .and. (imlsstr .eq. 0 ) ) then
        !do inp=1,Nparticle
          !call calc_rot_matrix(quat(:,inp),AA)
        !  call get_bbox_inds(bbox_inds,inp)
        !  call convex_hull_q12(bbox_inds,inp)
        !  call convex_hull_q22(bbox_inds,inp)
        !  call convex_hull_q32(bbox_inds,inp)
        !  call convex_hull_qc2(bbox_inds,inp)
        !enddo
      !endif
        !------------------------------------------------------------------------------
        tend = MPI_WTIME()
        wtime_vof = wtime_vof + (tend - tstart)

        tstart = MPI_WTIME()
        call divg
        call phcalc 

        call update_both_ghosts(n1,n2+1,dph,kstart,kend)
        
        call updvp  ! SOLENOIDAL VEL FIELD
        call prcalc  ! PRESSURE FIELD

        tend = MPI_WTIME()
        pressure_wtime = pressure_wtime + (tend - tstart)


        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)
        call update_both_ghosts(n1,n2,pr,kstart,kend)
        call update_both_ghosts(n1,n2,temp,kstart,kend)

!     ======================================================
!     End pressure correction
!     ======================================================

        enddo

      end

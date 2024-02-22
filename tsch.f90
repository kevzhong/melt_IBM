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

        !do i=1,n1
        !   do j = 1,n2
        !     do k = kstart-1,kend+1
        !      VOFx(i,j,k) = 1.
        !      VOFy(i,j,k) = 1.
        !      VOFz(i,j,k) = 1.
        !      VOFp(i,j,k) = 1.
        !    enddo
        !  enddo
        !enddo
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

        do inp=1,Nparticle
          call calc_rot_matrix(quat(:,inp),AA)

          call get_bbox_inds(bbox_inds,inp)

          call convex_hull_q12(AA,bbox_inds,inp)
          call convex_hull_q22(AA,bbox_inds,inp)
          call convex_hull_q32(AA,bbox_inds,inp)
          call convex_hull_qc2(AA,bbox_inds,inp)

        enddo
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


        !do i=1,n1
        !   do j = 1,n2
        !     do k = kstart-1,kend+1
        !      VOFx(i,j,k) = 1.
        !      VOFy(i,j,k) = 1.
        !      VOFz(i,j,k) = 1.
        !      VOFp(i,j,k) = 1.
        !    enddo
        !  enddo
        !enddo
        VOFx(:,:,:) = 1.
        VOFy(:,:,:) = 1.
        VOFz(:,:,:) = 1.
        VOFp(:,:,:) = 1.

        tstart = MPI_WTIME()
        call particle
        tend = MPI_WTIME()
        mls_wtime = mls_wtime + (tend - tstart)

        tstart = MPI_WTIME()

        !------------ KZ: update VOF, remove later since called by Newton--Euler -----
        do inp=1,Nparticle
          call calc_rot_matrix(quat(:,inp),AA)
          call get_bbox_inds(bbox_inds,inp)
          call convex_hull_q12(AA,bbox_inds,inp)
          call convex_hull_q22(AA,bbox_inds,inp)
          call convex_hull_q32(AA,bbox_inds,inp)
          call convex_hull_qc2(AA,bbox_inds,inp)
        enddo
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

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
        
        VOFx(:,:,:) = 1.
        VOFy(:,:,:) = 1.
        VOFz(:,:,:) = 1.
        VOFp(:,:,:) = 1.

        solid_mask(:,:,:) = .false.

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

        call hdnl1
        call hdnl2
        call hdnl3
        call hdnlte

        call invtr1 
        call invtr2      
        call invtr3
        !write(*,*) "Starting invtrte"
        call invtrte

        VOFx(:,:,:) = 1.
        VOFy(:,:,:) = 1.
        VOFz(:,:,:) = 1.
        VOFp(:,:,:) = 1.

        call particle

        !------------ KZ: update VOF, remove later since called by Newton--Euler -----
        !if ( (imlsfor.eq.1) .and. (imlsstr .eq. 0 ) ) then

        VOFx(:,:,:) = 1.
        VOFy(:,:,:) = 1.
        VOFz(:,:,:) = 1.
        VOFp(:,:,:) = 1.
        
        do inp=1,Nparticle
         call get_bbox_inds(bbox_inds,inp)
         call convex_hull_q12(bbox_inds,inp)
         call convex_hull_q22(bbox_inds,inp)
         call convex_hull_q32(bbox_inds,inp)
         call convex_hull_qc2(bbox_inds,inp)
        enddo
      !endif
        !------------------------------------------------------------------------------
        call divg
        call phcalc 

        call update_both_ghosts(n1,n2+1,dph,kstart,kend)
        
        call updvp  ! SOLENOIDAL VEL FIELD
        call prcalc  ! PRESSURE FIELD
        
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

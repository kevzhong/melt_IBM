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
      real,dimension(3,Nparticle) :: vel_m1,pos_m1,pos_k,om_m1

      integer :: ns, inp, ntr, nsub
      integer :: i,j,k
      

      beta=dt/ren*0.5d0

      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        do i=1,n1
           do j = 1,n2
             do k = kstart-1,kend+1
              ax(i,j,k) = 1.
              ay(i,j,k) = 1.
              az(i,j,k) = 1.
        enddo
        enddo
        enddo

        do inp=1,Nparticle
        call calc_rot_matrix(quat(:,inp),AA)

        call convex_hull_q12(AA,inp)
        call convex_hull_q22(AA,inp)
        call convex_hull_q32(AA,inp)
        enddo

        !call hdnl1
        !call hdnl2
        !call hdnl3

        !call invtr1 
        !call invtr2      
        !call invtr3
        call invtrte

        do i=1,n1
           do j = 1,n2
             do k = kstart-1,kend+1
              ax(i,j,k) = 1.
              ay(i,j,k) = 1.
              az(i,j,k) = 1.
        enddo
        enddo
        enddo
        call particle

        !call divg
        !call phcalc 

        !call update_both_ghosts(n1,n2+1,dph,kstart,kend)
        
        !call updvp  ! SOLENOIDAL VEL FIELD
        !call prcalc  ! PRESSURE FIELD

        !call update_both_ghosts(n1,n2,vx,kstart,kend)
        !call update_both_ghosts(n1,n2,vy,kstart,kend)
        !call update_both_ghosts(n1,n2,vz,kstart,kend)
        !call update_both_ghosts(n1,n2,pr,kstart,kend)
        call update_both_ghosts(n1,n2,temp,kstart,kend)

!     ======================================================
!     End pressure correction
!     ======================================================

        enddo

      end

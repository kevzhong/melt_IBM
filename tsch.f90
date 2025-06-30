      subroutine tschem 
      use param
      use local_arrays
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_aux
      !use stat_arrays, only: vxvyvz_rms_vol
      implicit none
      !real,dimension(3,3)     :: AA, AAT
      integer,dimension(3,2)     :: bbox_inds
      real,dimension(3,Nparticle) :: vel_m1,pos_m1,pos_k,om_m1
      real :: tstart, tend

      integer :: ns, inp, ntr, nsub
      integer :: i,j,k

      beta=dt/ren*0.5d0
      ! Accumulate Structural loads across all sub-steps
      !Fp(1:3) = 0.0d0
      !Ftau(1:3) = 0.0d0
      !Torq_p(1:3) = 0.0d0
      !Torq_tau(1:3) = 0.0d0

      ! Code timing
      wtime_vof = 0.
      eul_solve_wtime = 0.
      mls_wtime = 0.
      pressure_wtime = 0.


      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)
        
        if (timeflag) call tic(tstart)

        if ( .not. is_stationarySolid ) then
          call tagCells
        endif

        if (timeflag) call toc(tstart,tend,wtime_vof)
 
        

        if (timeflag) call tic(tstart)

        !call hdnl1
        !call hdnl2
        !call hdnl3
        !call hdnlte

        !call invtr1 
        !call invtr2      
        !call invtr3
        call invtrte

        if (timeflag) call toc(tstart,tend,eul_solve_wtime)

        if (timeflag) call tic(tstart)

        call particle

        if (timeflag) call toc(tstart,tend,mls_wtime)


        !if (timeflag) call tic(tstart)
        !call divg
        !call phcalc 
        !call update_both_ghosts(n1,n2+1,dph,kstart,kend)
        !call updvp  ! SOLENOIDAL VEL FIELD
        !call prcalc  ! PRESSURE FIELD
        !if (timeflag) call toc(tstart,tend,pressure_wtime)
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

      subroutine tschem 
      use param
      use local_arrays
      use mpih
      use phasefield
      use mpi_param, only: kstart,kend
      use local_aux
      !use stat_arrays, only: vxvyvz_rms_vol
      implicit none
      !real,dimension(3,3)     :: AA, AAT
      !integer,dimension(3,2)     :: bbox_inds
      !real,dimension(3,Nparticle) :: vel_m1,pos_m1,pos_k,om_m1
      integer :: ns, inp, ntr, nsub
      integer :: i,j,k

      beta=dt/ren*0.5d0

      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)
        

        if (meltmode .eq. 1) then
          call hdnl_phase
          call invtr_phase
          call update_both_ghosts(n1,n2,phi,kstart,kend)

          call compute_local_vmelt
        endif

        ! ! Compute surface metrics (OPTIONAL)
        ! if (pfmode .eq. 1) then
        !   call compute_psi_from_phi
        !   call compute_normals
        !   call compute_curvature
        ! endif

        call hdnl1
        call hdnl2
        call hdnl3

        if (which_hit .eq. 3) call add_linearHITForce

        if (meltmode .eq. 1) then
          ! Temperature only updated if melting: otherwise, field kept frozen
          call hdnlte
          !call hdnlte_QUICK
          call invtrte
          call update_both_ghosts(n1,n2,temp,kstart,kend)
        endif

        call invtr1 
        call invtr2
        call invtr3



        ! For single-phase
        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)

        !call particle

        call divg
        call phcalc 

        call update_both_ghosts(n1,n2+1,dph,kstart,kend)
        
        call updvp  ! SOLENOIDAL VEL FIELD
        call prcalc  ! PRESSURE FIELD

        
        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)
        call update_both_ghosts(n1,n2,pr,kstart,kend)
        !call update_both_ghosts(n1,n2,temp,kstart,kend)

!     ======================================================
!     End pressure correction
!     ======================================================


        enddo

      end

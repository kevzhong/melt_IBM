!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none
!==========================================================
!       read from input file bou.in
!==========================================================
        integer   :: n2, n3,n1
        integer   :: nsst, nread,pread, ntst, ireset, temp_restart
        real      :: tframe,tpin,tmax,walltimemax
        real      :: xlen, ylen, zlen
        real      :: pra,dt,resid,cflmax,tsta
        integer   :: starea
        real      :: dtmax,cfllim,cflfix
        real      :: tl,epsstar,kf_on_kmin, C_HIT, a0, k0
        integer   :: nson,idtv,forcing, which_hit

        integer :: pfmode, meltmode, melt_icond
        real      :: Usolid, Tmelt, Tliq, Tsol, latHeat, cpliquid
        real      :: rad_sph, halfthick, grv_depth, grv_width
!=================================================
!       end of input file
!=================================================
        real :: time
!******* Grid parameters**************************
        real :: dx2,dx3,dx1
        real :: dx2q,dx3q,dx1q
         
        real, allocatable, dimension(:) :: xc,xm
        real, allocatable, dimension(:) :: yc,ym
        real, allocatable, dimension(:) :: zc,zm

        ! volume integral
        !real, allocatable, dimension(:,:,:) :: VOFx, VOFy, VOFz, VOFp

        !logical :: timeflag = .false.
        logical :: specflag = .false.


!==========================================================
!******* Grid indices**************************************
        integer, allocatable, dimension(:) :: jmv,jpv
        integer, allocatable, dimension(:) :: imv,ipv
        integer, allocatable, dimension(:) :: jmhv
        integer, allocatable, dimension(:) :: kmv,kpv
!============================================================
!******* Variables for FFTW and Poisson solver****************
        real, dimension(13) :: ifx1
        integer*8 :: fwd_plan,bck_plan
        integer*8 :: fwdplan_1d,bckplan_1d
        real, allocatable, dimension(:) :: ao,ap,af
        real, allocatable, dimension(:) :: ak1,ak2,ak3
        
!===========================================================
!******* Other variables ***********************************
        integer  :: n2m, n3m, n1m,n1mh,n2mh
        integer  :: iaxsy
        real :: cflm 
        real :: ren, prandtl, pec, betagz
        real :: pi
        real :: al,ga,ro
        real :: beta, betatemp
        real :: qqmax,qqtot
        real :: re
        integer :: ntime
        real, dimension(1:3) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        real :: tempmin, tempmax, phimin, phimax

        integer :: nmodes
        integer, allocatable :: waveN(:)
        complex(kind=kind(0d0)), allocatable :: exp_I_kl_xi(:,:), exp_I_km_yj(:,:), exp_I_kn_zk(:,:)
        complex(kind=kind(0d0)), allocatable :: exp_I_kl_xsi(:,:), exp_I_km_ysj(:,:), exp_I_kn_zsk(:,:)
        complex(kind=kind(0d0)), allocatable :: bhat(:,:,:,:)
        

        logical :: ismaster = .false.
        logical :: solvestructure = .false.
        logical :: mlsforcing = .false.

        
      end module param
      
!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******
      module local_arrays
      use param
        implicit none
        real,allocatable,dimension(:,:,:) :: vx,vy,vz
        real,allocatable,dimension(:,:,:) :: temp
        real,allocatable,dimension(:,:,:) :: qbuf,forcx,forcy,forcz
        real,allocatable,dimension(:,:,:) :: pr,rhs
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3,rut
        real,allocatable,dimension(:,:,:) :: qcap, htemp
        real,allocatable,dimension(:,:,:) :: dph,dq

        ! QUICK
        real,allocatable,dimension(:,:,:) :: temp2, xflux_imh, yflux_jmh, zflux_kmh
        
      end module local_arrays

!===============================================================
      module stat_arrays
       implicit none
       real,allocatable, dimension(:,:,:) :: vx_me,vx_rms,pr_me,pr_rms
       real,allocatable, dimension(:,:,:) :: vy_me,vz_me,vy_rms,vz_rms 
       integer :: timeint_cdsp
        real :: vxvyvz_rms_vol
      end module stat_arrays
!=====================================================       
      module mpih
        use mpi
        implicit none
        !include 'mpif.h'
        integer :: myid, numtasks, ierr
        integer :: p_row, p_col, my_p_row, my_p_col
        integer, parameter :: master=0
        integer, parameter :: lvlhalo=1
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer :: MCP = MPI_DOUBLE_COMPLEX
        integer :: STATUS(MPI_STATUS_SIZE,4)
        integer :: req(1:4)
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
      end module mpih
      
      module mpi_param
        implicit none
        integer :: istart,iend,jstart,jend, kstart,kend
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 

        real, allocatable, dimension(:,:) :: buf_n1n2
      end module mpi_param

      module local_aux
       use param
       implicit none
       real,allocatable,dimension(:,:,:) :: vorx, vory, vorz
       real,allocatable,dimension(:,:,:) :: diss, tke, chi
      end module local_aux

module phasefield
        implicit none
        real,allocatable,dimension(:,:,:) :: phi
        real,allocatable,dimension(:,:,:) :: nhat_x,nhat_y,nhat_z,psi,curv,vmelt
        real,allocatable,dimension(:,:,:) :: hphi, ruphi

        real :: pf_eps, pf_A, D_pf
end module phasefield
module modspec
        implicit none
        integer*8 :: specplan
        complex, allocatable :: uhat(:,:,:)
        !complex, allocatable :: uhat(:)
end module modspec


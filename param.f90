!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none
!==========================================================
!       read from input file bou.in
!==========================================================
        integer   :: n2, n3,n1
        integer   :: nsst, nwrit, nread,pread, ntst, ireset
        real      :: tframe,tpin,tmax,walltimemax
        real      :: xlen, ylen, zlen
        real      :: pra,dt,resid,cflmax,tsta
        integer   :: starea
        real      :: dtmax,cfllim
        real      :: tl,epsstar,kfmax
        integer   :: nson,idtv,forcing
        real      :: Tmelt
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
        real, allocatable, dimension(:,:,:) :: ax, ay, az
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
        real :: keta
        real :: ren, prandtl, pec
        real :: pi
        real :: al,ga,ro
        real :: beta, betatemp
        real :: qqmax,qqtot
        real :: re
        integer :: ntime
        real, dimension(1:3) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        complex,allocatable,dimension(:,:) :: term1a,term1b
        complex,allocatable,dimension(:,:) :: term2a,term2b
        complex,allocatable,dimension(:,:) :: term3a,term3b
        real, dimension(6,7,7,7) :: bcoefs
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
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3
        real,allocatable,dimension(:,:,:) :: qcap
        real,allocatable,dimension(:,:,:) :: dph,dq
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
        implicit none
        include 'mpif.h'
        integer :: myid, numtasks, ierr
        integer, parameter :: master=0
        integer, parameter :: lvlhalo=2
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
      end module mpi_param

      module local_aux
       use param
       implicit none
       !real, allocatable, dimension(:,:,:) :: vxc, vyc, vzc !KZ unused
       !real, allocatable, dimension(:,:,:) :: matderxc, matderx !KZ unused
       !real, allocatable, dimension(:,:,:) :: matderyc, matdery !KZ unused
       !real, allocatable, dimension(:,:,:) :: matderzc, matderz !KZ unused
       real,allocatable,dimension(:,:,:) :: vorx, vory, vorz
       !real,allocatable,dimension(:,:,:) :: vxo, vyo, vzo !KZ unused
      end module local_aux

      module mls_param
      implicit none
      
      integer :: Nparticle
      parameter( Nparticle=1)

      !=================================================
      !       read from input file part.in
      !==========================================================
      integer   :: imlsfor,imlsstr
      real      :: wcon, wscl, dens_ratio
      real      :: inert_fac, ref_pos_fac
      character(50)  gtsfx
      real      :: rad_p
      !=================================================
      !       end of input file
      !=================================================


!     --------variables for structural solver------------------------
      integer  :: max_n_edge_of_vert
      integer, parameter :: nel=27

      integer, dimension(:),   allocatable :: n_edge_of_vert
      integer, dimension(:,:), allocatable :: vert_of_edge
      integer, dimension(:,:), allocatable :: vert_of_face
      integer, dimension(:,:), allocatable :: edge_of_face
      integer, dimension(:,:), allocatable :: vert_of_vert
      integer, dimension(:,:), allocatable :: edge_of_vert
      integer, dimension(:,:), allocatable :: face_of_edge

      integer, dimension(:,:,:), allocatable :: pind
      integer, dimension(:,:), allocatable :: pind1
      real,dimension(:,:,:),     allocatable :: dismax

      real, dimension(:,:,:), allocatable :: tri_ver, vel_tri
      real, dimension(:,:,:), allocatable :: tri_bar, tri_nor
 
      
      real, dimension(:,:), allocatable :: dist
      real, dimension(:,:), allocatable :: sur, vol

      real :: celvol
      real :: chi, I11, I22, I33

      real, dimension(:,:),   allocatable :: xyz0
      real, dimension(:,:,:), allocatable :: xyzv
      real, dimension(:,:,:), allocatable :: dxyz_CM_b
      real, dimension(:,:,:), allocatable :: dxyz_CM_s



      integer :: maxnv,maxne,maxnf

      !-- Particle vars
      real, dimension(:,:), allocatable :: fpxyz,     ftxyz
      real, dimension(:,:), allocatable :: pos_CM,    vel_CM, a_CM
      real, dimension(:,:), allocatable :: omega_b,   omega_dot_b, alpha_b
      real, dimension(:,:), allocatable :: om_b_sqr,  om_b_sqr_m1
      real, dimension(:,:), allocatable :: u_tot,     u_tot_m1
      real, dimension(:,:), allocatable :: r_x_u_tot, r_x_u_tot_m1
      real, dimension(:,:), allocatable :: omega_s
      real, dimension(:,:), allocatable :: tail_head
      real, dimension(:,:), allocatable :: quat, quat_m1, quat_dot

      real, dimension(:),allocatable :: Surface
      real, dimension(:),allocatable :: Volume

      real   alph_dx, invdx1celvol, invdx1dt, invwcon, alph_dx_invwcon
      real, dimension(:), allocatable :: cfac
      real :: Hboxx
      real, dimension (3,3) :: i_inv, i_inv2

      !-- mlsWeight
      real, dimension(:,:,:), allocatable :: ptxAB_q1,ptxAB_q2,ptxAB_q3
      real, dimension(:,:,:), allocatable :: ptxAB_temp
      !real,dimension(:,:,:,:), allocatable :: ptxAB_pr

      end module mls_param

!====================================================
      module mls_local
        use param
        implicit none
        real, dimension(:,:,:), allocatable :: for_xc, for_yc, for_zc, for_temp
      end module mls_local


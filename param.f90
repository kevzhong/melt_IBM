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
        real      :: tl,epsstar,kf_on_kmin, C_HIT
        integer   :: nson,idtv,forcing, which_hit
        real      :: Tmelt, Tliq, Tsol, latHeat, cpliquid
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
        real, allocatable, dimension(:,:,:) :: VOFx, VOFy, VOFz, VOFp
        !real, allocatable, dimension(:,:,:) :: usolid_x, usolid_y, usolid_z
        !real, allocatable, dimension(:,:,:) :: d_UsolidT_dxj
        !logical, allocatable, dimension(:,:,:) :: solid_mask
        real, allocatable, dimension(:,:,:) :: sdf

        logical :: timeflag = .false.
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

        integer :: nmodes
        integer, allocatable :: waveN(:)
        complex(kind=kind(0d0)), allocatable :: exp_I_kl_xi(:,:), exp_I_km_yj(:,:), exp_I_kn_zk(:,:)
        complex(kind=kind(0d0)), allocatable :: exp_I_kl_xsi(:,:), exp_I_km_ysj(:,:), exp_I_kn_zsk(:,:)
        complex(kind=kind(0d0)), allocatable :: bhat(:,:,:,:)
        

        logical :: ismaster = .false.
        logical :: solvestructure = .false.
        logical :: mlsforcing = .false.


        ! For code timing / benchmarking
        real :: wtime_vof
        real :: eul_solve_wtime
        real :: mls_wtime
        real :: pressure_wtime
        real :: hit_wtime
        real :: wtime_total
        
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
      end module local_aux

      module mls_param
      implicit none
      
      integer :: Nparticle
      parameter( Nparticle=1)

      !=================================================
      !       read from input file part.in
      !==========================================================
      integer   :: imlsfor,imlsstr, imelt, iremesh
      real      :: wcon, wscl, dens_ratio
      real      :: inert_fac, ref_pos_fac
      character(50)  gtsfx
      real      :: rad_p
      integer :: tagType
      logical :: is_stationarySolid = .true.
      logical :: initial_tag = .true.
      !integer   :: VERTBUFFER !KZ Max no. of faces adjoining the vertices, should be precomputed in geom. pre-processing
      !=================================================
      !       end of input file
      !=================================================


!     --------variables for structural solver------------------------
      !integer  :: max_n_edge_of_vert
      integer, parameter :: nel=27

      !integer, dimension(:),   allocatable :: n_edge_of_vert
      integer, dimension(:,:,:), allocatable :: vert_of_edge
      integer, dimension(:,:,:), allocatable :: vert_of_face
      integer, dimension(:,:,:), allocatable :: edge_of_face
      integer, dimension(:,:,:), allocatable :: face_of_edge

      !Re-meshing data structures
      logical, dimension(:,:), allocatable :: isGhostFace
      logical, dimension(:,:), allocatable :: isGhostEdge
      logical, dimension(:,:), allocatable :: isGhostVert
      real, dimension(:,:), allocatable :: eLengths
      real, dimension(:,:), allocatable :: skewness
      !real :: PERC_Athresh, A_thresh, skew_thresh, V_ON_VE_PERC, V_thresh
      real :: PERC_Ethresh, E_thresh, V_ON_VE_PERC, V_thresh
      
      logical, dimension(:), allocatable :: rm_flag ! Remeshing flag
      logical, dimension(:,:), allocatable :: anchorVert
      logical, dimension(:,:), allocatable :: flagged_edge

      integer, dimension(:,:,:), allocatable :: pind
      !integer, dimension(:,:,:), allocatable :: pindv
      integer, dimension(:,:), allocatable :: pind1
      !real,dimension(:,:,:),     allocatable :: dismax

      real, dimension(:,:,:), allocatable :: tri_ver, vel_tri
      real, dimension(:,:,:), allocatable :: tri_bar, tri_nor
      real, dimension(:,:,:), allocatable :: vert_nor !KZ: normal vectors of vertices

      
      real, dimension(:,:), allocatable :: sur, vol

      real :: celvol
      !real :: chi, I11, I22, I33

      real, dimension(:,:),   allocatable :: xyz0
      real, dimension(:,:,:), allocatable :: xyzv

      !COM-relative vertex coordinates
      real, dimension(:,:,:), allocatable :: dxyzv_s
      real, dimension(:,:,:), allocatable :: dxyz_s



      integer :: maxnv,maxne,maxnf

      !-- Particle vars
      real, dimension(:,:), allocatable :: pos_CM,    vel_CM
      real, dimension(:,:), allocatable :: omega_c
      real, dimension(:,:), allocatable :: quat
      !quat_m1, quat_dot, quat_dot_m1

      real, dimension(3,3,Nparticle) :: InertTensor

      real, dimension(:),allocatable :: Surface
      real, dimension(:),allocatable :: Volume

      real    invdx1dt
      real, dimension(:,:), allocatable :: cfac
      real :: h_eulerian, A_eulerian

      !-- mlsWeight
      real, dimension(:,:,:), allocatable :: ptxAB_q1,ptxAB_q2,ptxAB_q3
      real, dimension(:,:,:), allocatable :: ptxAB_temp
      real, dimension(:,:), allocatable :: tau_n1, tau_n2, tau_n3 ! Structural loads
      real, dimension(:,:,:), allocatable ::  press_n_tri ! Structural loads
      real, dimension(:,:), allocatable :: r_x_tau_n1, r_x_tau_n2, r_x_tau_n3 ! Structural torques
      real, dimension(:,:,:), allocatable :: r_x_prn_tri ! Structural torques
      real, dimension(:,:), allocatable :: qw_o, qw_i ! Normal gradients at interface faces (outward and inward dirn)
      real, dimension(:,:), allocatable :: qw_oVert, qw_iVert ! Normal gradients at interface VERTICES (outward and inward dirn)
      real, dimension(:,:), allocatable :: Avert ! Area of each vertex
      real, dimension(:,:,:), allocatable :: vmelt, vmelt_m1 ! local melting velocity at vertices
      real, dimension(:,:),allocatable :: int_tau_dA, int_r_x_tau_dA,int_prn_dA, int_r_x_prn_dA
      real, dimension(:,:),allocatable :: int_tau_dA_m1, int_r_x_tau_dA_m1,int_prn_dA_m1, int_r_x_prn_dA_m1
      real, dimension(3) :: Fp, Ftau ! Integrated loads
      real, dimension(3) :: Torq_p, Torq_tau ! Integrated torques

      end module mls_param

!====================================================
      module mls_local
        use param
        implicit none
        real, dimension(:,:,:), allocatable :: for_xc, for_yc, for_zc, for_temp
      end module mls_local

module modspec
        implicit none
        integer*8 :: specplan
        complex, allocatable :: uhat(:,:,:)
        !complex, allocatable :: uhat(:)

end module modspec


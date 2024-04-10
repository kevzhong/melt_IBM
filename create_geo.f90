subroutine setup_particles
  use param
  use mls_param
  implicit none
  integer :: i,j,k,inp,count
  integer :: Npx,Npy,Npz
  real :: a,b,c
  real :: x1,x2,x3
  real,dimension(3,3):: AA,AAT
  real :: alpha
  real :: angle
  real :: I1, I2, I3
  real, dimension(4) :: Q_buffer

  ! Assign COMs for each particle and initialise rigid-body parameters (vel_COM, orientation, etc)

  Npx = 1; Npy = 1; Npz = 1

  ! KZ: for future with multiple particles, can transform the inertia tensor with just
  ! I' = A * I * A^T

  !--------Pre-rotate geometry------------------
  angle = 30.0
  Q_buffer(1) = cosd(angle / 2.0)
  Q_buffer(2) = 0.0
  Q_buffer(3) = sind(angle / 2.0)
  Q_buffer(4) = 0.0

  call calc_rot_matrix(Q_buffer,AA)
  AAT = transpose(AA)
  do i = 1,maxnv
    xyz0(:,i) = matmul(AAT, xyz0(:,i) )
  enddo
 !----------------------------------------------
  
  call calc_rigidBody_params(pos_CM(:,1),Volume(1),InertTensor(:,:,1),maxnv,maxnf,&
  xyz0(:,:),vert_of_face(:,:,1),isGhostFace(:,1) )

  ! KZ: note the hard-coded single particle for now
  pos_CM(1,1) = 0.5*xlen
  pos_CM(2,1) = 0.5*ylen
  pos_CM(3,1) = 0.5*zlen

  !write(*,*), "Volume is", Volume(1)
  
     do inp=1,Nparticle
     ! Initialise rigid body variables
      !call random_number(a)
      ! b=a ; c=a
      
!      a=pi/4  ; b=0.0 ; c=0.0
      !quat(1,inp) = cos(0.5*(a))*cos(0.5*(b+c))
      !quat(2,inp) = sin(0.5*(a))*cos(0.5*(b-c))
      !quat(3,inp) = sin(0.5*(a))*sin(0.5*(b-c))
      !quat(4,inp) = cos(0.5*(a))*sin(0.5*(b+c))

      ! Identity quaternion
      quat(1,inp) = 1.0
      quat(2,inp) = 0.0
      quat(3,inp) = 0.0
      quat(4,inp) = 0.0


      omega_c(:,inp) = 0.0
      alpha_b(:,inp) = 0.0
      omega_dot_b = 0.0

      quat_dot = 0.0
      quat_dot_m1 = 0.0
      
      om_b_sqr = 0.0

      u_tot_m1 = 0.
      r_x_u_tot_m1 = 0.

      vel_CM    = 0.
      a_CM      = 0.
      u_tot     = 0.
      r_x_u_tot = 0.

      !KZ : for initialising stationary
      vel_tri(:,:,inp) = 0.

     do j=1,Nparticle
     if((j.ne.inp).and.ismaster)then
     write(*,*)'dist part N',inp,'from part N',j, sqrt((pos_CM(1,inp)-pos_CM(1,j))**2+ &
                            (pos_CM(2,inp)-pos_CM(2,j))**2+ &
                            (pos_CM(3,inp)-pos_CM(3,j))**2)
     endif
     enddo
   end do

end subroutine

subroutine init_geomCoords
  use param
  use mls_param
  use geom
  use mpih
  implicit none
  integer :: inp


  ! Set the vertex coordinates of the geometry based on assigned COM

  do inp = 1,Nparticle
    ! KZ store vertex coordinates
    xyzv(1:3,:,inp) = xyz0(1:3,:) !+ pos_CM(:,inp)

    ! COM-relative assignment
    xyzv(1,:,inp) = xyzv(1,:,inp) + pos_CM(1,inp)
    xyzv(2,:,inp) = xyzv(2,:,inp) + pos_CM(2,inp)
    xyzv(3,:,inp) = xyzv(3,:,inp) + pos_CM(3,inp)

    call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) 
  enddo




end subroutine


 subroutine import_particles
  use param
  use mls_param
  use mpih
  implicit none
  character(200) :: dname 
  real :: var = 0

    dname = trim("pos_CM")
    call hdf_read_2d(pos_CM,shape(pos_CM),dname)
    dname = trim("vel_CM")
    call hdf_read_2d(vel_CM,shape(vel_CM),dname)
    dname = trim("quat")
    call hdf_read_2d(quat,shape(quat),dname)
    dname = trim("quat_dot")
    call hdf_read_2d(quat_dot,shape(quat_dot),dname)
    dname = trim("omega_c")
    call hdf_read_2d(omega_c,shape(omega_c),dname)
    dname = trim("om_b_sqr")
    call hdf_read_2d(om_b_sqr,shape(om_b_sqr),dname)


    dname = trim("u_tot")
    call hdf_read_2d(u_tot,shape(u_tot),dname)
    dname = trim("r_x_u_tot")
    call hdf_read_2d(r_x_u_tot,shape(r_x_u_tot),dname)

    dname = trim("volume")
    call hdf_read_1d(Volume,Nparticle,dname)
    dname = trim("total_area")
    call hdf_read_1d(Surface,Nparticle,dname)


    ! ------------------- GEOMETRY --------------------------------------
    ! ----------------------------   Edge information ----------------------------

    !----- Connectivity information ---------------
    dname = trim("vert_of_edge")
    call hdf_read_3dInt(vert_of_edge,shape(vert_of_edge),dname)

    dname = trim("face_of_edge")
    call hdf_read_3dInt(face_of_edge,shape(face_of_edge),dname)

    dname = trim("isGhostEdge")
    call hdf_read_2dInt(isGhostEdge,shape(isGhostEdge),dname)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------

    dname = trim("eLengths")
    call hdf_read_2d(eLengths,shape(eLengths),dname)

    !---------------------------------------------

    
    ! ----------------------------   Face information ----------------------------

    !----- Connectivity information ---------------
    dname = trim("vert_of_face")
    call hdf_read_3dInt(vert_of_face,shape(vert_of_face),dname)

    dname = trim("edge_of_face")
    call hdf_read_3dInt(edge_of_face,shape(edge_of_face),dname)

    dname = trim("isGhostFace")
    call hdf_read_2dInt(isGhostFace,shape(isGhostFace),dname)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------

    dname = trim("tri_bar")
    call hdf_read_3d(tri_bar,shape(tri_bar),dname)

    dname = trim("tri_nor")
    call hdf_read_3d(tri_nor,shape(tri_nor),dname)

    dname = trim("Atri")
    call hdf_read_2d(sur,shape(sur),dname)

    dname = trim("skewness")
    call hdf_read_2d(skewness,shape(skewness),dname)
    !---------------------------------------------

    ! ----------------------------   Vertex information --------------------------

    !----- Connectivity information ---------------

    dname = trim("isGhostVert")
    call hdf_read_2dInt(isGhostVert,shape(isGhostVert),dname)

    dname = trim("xyzv")
    call hdf_read_3d(xyzv,shape(xyzv),dname)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------
    dname = trim("vert_nor")
    call hdf_read_3d(vert_nor,shape(vert_nor),dname)

    dname = trim("Avert")
    call hdf_read_2d(Avert,shape(Avert),dname)

    !---------------------------------------------

    !--------------------------------------------------------------------


    call mpi_bcast(pos_CM,size(pos_cm),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(vel_CM,size(vel_cm),mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(quat,     size(quat),     mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(quat_dot, size(quat_dot), mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(omega_c,  size(omega_c),  mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(om_b_sqr, size(om_b_sqr), mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(u_tot,    size(u_tot),    mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(r_x_u_tot,size(r_x_u_tot),mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(Volume,size(Volume),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(Surface,size(Surface),mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(vert_of_edge,size(vert_of_edge),mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(face_of_edge,size(face_of_edge),mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(isGhostEdge,size(isGhostEdge),mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(eLengths,size(eLengths),mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(vert_of_face,size(vert_of_face),mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(edge_of_face,size(edge_of_face),mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(isGhostFace,size(isGhostFace),mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(tri_bar,size(tri_bar),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(tri_nor,size(tri_nor),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(sur,size(sur),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(skewness,size(skewness),mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(isGhostVert,size(isGhostVert),mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(xyzv,size(xyzv),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(vert_nor,size(vert_nor),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(Avert,size(Avert),mpi_double,0,mpi_comm_world,ierr)

 end subroutine

 subroutine continua_particle
  use param
  use mpih
  use mls_param
 
  implicit none

  character(70) filename
  character(30) dataset
  character(5) ipfi

  integer, dimension(2) :: dims

  dims(1)=n1m
  dims(2)=n2m

  filename = 'continuation/particles.h5'
   

  if (myid.eq.0) then 
    call hdf5_create_blank_file(filename)

    !---------------------- Newton--Euler, rigid-body quantities ------------
    dataset = trim("pos_CM")
    call hdf_write_2d(pos_CM,(/3,Nparticle/),dataset)
    dataset = trim("vel_CM")
    call hdf_write_2d(vel_CM,(/3,Nparticle/),dataset)
    dataset = trim("quat")
    call hdf_write_2d(quat,(/4,Nparticle/),dataset)
    dataset = trim("quat_dot")
    call hdf_write_2d(quat_dot,(/4,Nparticle/),dataset)
    dataset = trim("omega_c")
    call hdf_write_2d(omega_c,(/3,Nparticle/),dataset)
    dataset = trim("om_b_sqr")
    call hdf_write_2d(om_b_sqr,(/3,Nparticle/),dataset)

    dataset = trim("u_tot")
    call hdf_write_2d(u_tot,(/3,Nparticle/),dataset)
    dataset = trim("r_x_u_tot")
    call hdf_write_2d(r_x_u_tot,(/3,Nparticle/),dataset)


    ! Global quantities: Volume, surface area, 
    ! Inertia tensor, principal axes, etc (eventually)

    dataset = trim("volume")
    call hdf_write_1d(Volume,Nparticle,dataset)
    dataset = trim("total_area")
    call hdf_write_1d(Surface,Nparticle,dataset)

    !-------------- BEGIN GEOMETRY INFORMATION -----------------------------------

    ! ----------------------------   Edge information ----------------------------

    !----- Connectivity information ---------------
    dataset = trim("vert_of_edge")
    call hdf_write_3dInt(vert_of_edge,(/2,maxne,Nparticle/),dataset)

    dataset = trim("face_of_edge")
    call hdf_write_3dInt(face_of_edge,(/2,maxne,Nparticle/),dataset)

    dataset = trim("isGhostEdge")
    call hdf_write_2dInt(isGhostEdge,(/maxne,Nparticle/),dataset)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------

    dataset = trim("eLengths")
    call hdf_write_2d(eLengths,(/maxne,Nparticle/),dataset)

    !---------------------------------------------

    
    ! ----------------------------   Face information ----------------------------

    !----- Connectivity information ---------------
    dataset = trim("vert_of_face")
    call hdf_write_3dInt(vert_of_face,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("edge_of_face")
    call hdf_write_3dInt(edge_of_face,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("isGhostFace")
    call hdf_write_2dInt(isGhostFace,(/maxnf,Nparticle/),dataset)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------

    dataset = trim("tri_bar")
    call hdf_write_3d(tri_bar,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("tri_nor")
    call hdf_write_3d(tri_nor,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("Atri")
    call hdf_write_2d(sur,(/maxnf,Nparticle/),dataset)

    dataset = trim("skewness")
    call hdf_write_2d(skewness,(/maxnf,Nparticle/),dataset)
    !---------------------------------------------


    ! ----------------------------   Vertex information --------------------------

    !----- Connectivity information ---------------

    dataset = trim("isGhostVert")
    call hdf_write_2dInt(isGhostVert,(/maxnv,Nparticle/),dataset)

    dataset = trim("xyzv")
    call hdf_write_3d(xyzv,(/3,maxnv,Nparticle/),dataset)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------
    dataset = trim("vert_nor")
    call hdf_write_3d(vert_nor,(/3,maxnv,Nparticle/),dataset)

    dataset = trim("Avert")
    call hdf_write_2d(Avert,(/maxnv,Nparticle/),dataset)

    !---------------------------------------------



  end if
  end subroutine


subroutine print_particle_info
  use param
  use mls_param
  use mpih
  implicit none
  real :: ddx3

  ddx3 = n3m/zlen

  if(myid.eq.0) then
    write(*,*) 'Ave edge/dx: ',    (sum(eLengths(:,1))/float(maxne)) * ddx3
    write(*,*) 'Max edge/dx: ',    maxval(eLengths(:,1)) * ddx3
    write(*,*) 'Min edge/dx: ',    minval(eLengths(:,1)) * ddx3
  end if
end subroutine


subroutine write_h5_geo
  use param
  use mpih
  use mls_param
 
  implicit none

  character(70) filename
  character(30) dataset
  character(5) ipfi

  integer, dimension(2) :: dims
  real tprfi
  integer :: itime

  tprfi = 1/tframe
  itime=nint(time*tprfi)


  dims(1)=n1m
  dims(2)=n2m
  82 format(i7.7)
  write(ipfi,82) itime

  filename='flowmov/partHist_'//ipfi
  !filename = 'continuation/particles.h5'
   

  if (myid.eq.0) then 
    call hdf5_create_blank_file(filename)

    !---------------------- Newton--Euler, rigid-body quantities ------------
    dataset = trim("pos_CM")
    call hdf_write_2d(pos_CM,(/3,Nparticle/),dataset)
    dataset = trim("vel_CM")
    call hdf_write_2d(vel_CM,(/3,Nparticle/),dataset)
    dataset = trim("quat")
    call hdf_write_2d(quat,(/4,Nparticle/),dataset)
    dataset = trim("quat_dot")
    call hdf_write_2d(quat_dot,(/4,Nparticle/),dataset)
    dataset = trim("omega_c")
    call hdf_write_2d(omega_c,(/3,Nparticle/),dataset)
    dataset = trim("om_b_sqr")
    call hdf_write_2d(om_b_sqr,(/3,Nparticle/),dataset)

    dataset = trim("u_tot")
    call hdf_write_2d(u_tot,(/3,Nparticle/),dataset)
    dataset = trim("r_x_u_tot")
    call hdf_write_2d(r_x_u_tot,(/3,Nparticle/),dataset)


    ! Global quantities: Volume, surface area, 
    ! Inertia tensor, principal axes, etc (eventually)

    dataset = trim("volume")
    call hdf_write_1d(Volume,Nparticle,dataset)
    dataset = trim("total_area")
    call hdf_write_1d(Surface,Nparticle,dataset)

    !-------------- BEGIN GEOMETRY INFORMATION -----------------------------------

    ! ----------------------------   Edge information ----------------------------

    !----- Connectivity information ---------------
    dataset = trim("vert_of_edge")
    call hdf_write_3dInt(vert_of_edge,(/2,maxne,Nparticle/),dataset)

    dataset = trim("face_of_edge")
    call hdf_write_3dInt(face_of_edge,(/2,maxne,Nparticle/),dataset)

    dataset = trim("isGhostEdge")
    call hdf_write_2dInt(isGhostEdge,(/maxne,Nparticle/),dataset)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------

    dataset = trim("eLengths")
    call hdf_write_2d(eLengths,(/maxne,Nparticle/),dataset)

    !---------------------------------------------

    
    ! ----------------------------   Face information ----------------------------

    !----- Connectivity information ---------------
    dataset = trim("vert_of_face")
    call hdf_write_3dInt(vert_of_face,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("edge_of_face")
    call hdf_write_3dInt(edge_of_face,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("isGhostFace")
    call hdf_write_2dInt(isGhostFace,(/maxnf,Nparticle/),dataset)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------

    dataset = trim("tri_bar")
    call hdf_write_3d(tri_bar,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("tri_nor")
    call hdf_write_3d(tri_nor,(/3,maxnf,Nparticle/),dataset)

    dataset = trim("Atri")
    call hdf_write_2d(sur,(/maxnf,Nparticle/),dataset)

    dataset = trim("skewness")
    call hdf_write_2d(skewness,(/maxnf,Nparticle/),dataset)
    !---------------------------------------------


    ! ----------------------------   Vertex information --------------------------

    !----- Connectivity information ---------------

    dataset = trim("isGhostVert")
    call hdf_write_2dInt(isGhostVert,(/maxnv,Nparticle/),dataset)

    dataset = trim("xyzv")
    call hdf_write_3d(xyzv,(/3,maxnv,Nparticle/),dataset)
    !---------------------------------------------

    !----- Derived geoemtric quantities ----------
    dataset = trim("vert_nor")
    call hdf_write_3d(vert_nor,(/3,maxnv,Nparticle/),dataset)

    dataset = trim("Avert")
    call hdf_write_2d(Avert,(/maxnv,Nparticle/),dataset)

    !---------------------------------------------



  end if
  end subroutine

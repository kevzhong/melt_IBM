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

  Npx = 1; Npy = 1; Npz = 1

  angle = 0.0
  pos_CM(:,1)=0.5*xlen
!  pos_CM(3,1)=0.8*zlen
!  pos_CM(1,2)=0.7*xlen
!  pos_CM(2:3,:)=0.5*xlen
!  pos_CM(1,1)=0.3*xlen
!  if(Nparticle.gt.1)then
!  do inp=1,Nparticle
!     if(inp.eq.1)then
   
!     call random_number(x1)
!     call random_number(x2)
!     call random_number(x3)

!     pos_CM(1,inp) = x1*xlen
!     pos_CM(2,inp) = x2*ylen
!     pos_CM(3,inp) = x3*zlen
!     else
!     count=inp
     
!     call random_number(x1)
!     call random_number(x2)
!     call random_number(x3)

!     pos_CM(1,count) = x1*xlen
!     pos_CM(2,count) = x2*ylen
!     pos_CM(3,count) = x3*zlen
!     end if

!40 continue
!     do j=1,count
!     if((count.ne.j).and.(sqrt((pos_CM(1,count)-pos_CM(1,j))**2+ &
!                               (pos_CM(2,count)-pos_CM(2,j))**2+  &
!                               (pos_CM(3,count)-pos_CM(3,j))**2).le.2.5))then

!     call random_number(x1)
!     call random_number(x2)
!     call random_number(x3)

!     pos_CM(1,count) = x1*xlen
!     pos_CM(2,count) = x2*ylen
!     pos_CM(3,count) = x3*zlen
!     endif
!     enddo

!     do j=1,count
!     if((j.ne.count).and.(sqrt((pos_CM(1,count)-pos_CM(1,j))**2+ &
!                               (pos_CM(2,count)-pos_CM(2,j))**2+  &
!                               (pos_CM(3,count)-pos_CM(3,j))**2).le.2))then
!     go to 40
!     end if
!     enddo

!     end if !inp.ne.1
!     enddo  !!
!     end if !Nparticle

     do inp=1,Nparticle
     ! Initialise rigid body variables
      call random_number(a)
       b=a ; c=a
      
!      a=pi/4  ; b=0.0 ; c=0.0
      quat(1,inp) = cos(0.5*(a))*cos(0.5*(b+c))
      quat(2,inp) = sin(0.5*(a))*cos(0.5*(b-c))
      quat(3,inp) = sin(0.5*(a))*sin(0.5*(b-c))
      quat(4,inp) = cos(0.5*(a))*sin(0.5*(b+c))


      omega_b(:,inp) = 0.0
      alpha_b(:,inp) = 0.0
      omega_dot_b = 0.0

      quat_dot = 0.0
      om_b_sqr = 0.0

      u_tot_m1 = 0.
      r_x_u_tot_m1 = 0.

      vel_CM    = 0.
      a_CM      = 0.
      u_tot     = 0.
      r_x_u_tot = 0.

     do j=1,Nparticle
     if((j.ne.inp).and.ismaster)then
     write(*,*)'dist part N',inp,'from part N',j, sqrt((pos_CM(1,inp)-pos_CM(1,j))**2+ &
                            (pos_CM(2,inp)-pos_CM(2,j))**2+ &
                            (pos_CM(3,inp)-pos_CM(3,j))**2)
     endif
     enddo
   end do

end subroutine

subroutine set_particle_rad
  use param
  use mls_param
  use geom
  use mpih
  implicit none
  integer :: inp,i
  integer :: v1,v2,v3
  real    :: AAT_P(3,3)

  AAT_P = princ_axis_rotm()

  if(myid.eq.0) print *, AAT_P
  do inp = 1,Nparticle
    ! This is r
    do i=1,maxnf

      v1=vert_of_face(1,i)
      v2=vert_of_face(2,i)
      v3=vert_of_face(3,i)

      dxyz_CM_b(:,i,inp) = matmul(AAT_P, (xyz0(:,v1) + xyz0(:,v2) + xyz0(:,v3)) / 3.d0)
    end do

  enddo

end subroutine
 

subroutine set_xyz
  use param
  use mls_param
  implicit none
  real,dimension(3)  :: om_dCM, pos, vel
  integer :: i,inp
  real,dimension(3,3):: AA,AAT,AAR
  real, dimension(2,2) :: Rot
  real :: radius,angle, om,tp
  real ::zmin,zmax

   do inp = 1,Nparticle

      call calc_rot_matrix(quat(:,inp),AA)

      !  compute tranpose
      AAT = transpose(AA)
      ! x-rot pi/2
      tail_head(:,inp) = AAT(:,3)
!      AAR1(1,1) = 1
!      AAR1(1,2) = 0
!      AAR1(1,3) = 0
!      AAR1(2,1) = 0
!      AAR1(2,2) = 0
!      AAR1(2,3) = -1
!      AAR1(3,1) = 0
!      AAR1(3,2) = 1
!      AAR1(3,3) = 0
!      AAR1=transpose(AAR1)
      ! y-rot pi/2
!      AAR(1,1) = 1
!      AAR(1,2) = 0
!      AAR(1,3) = 0

!      AAR(2,1) = 0
!      AAR(2,2) = -1
!      AAR(2,3) = 0

!      AAR(3,1) = 0
!      AAR(3,2) = 0
!      AAR(3,3) = 1
     ! AAR=transpose(AAR)

      do i=1,maxnf
         !-- position
         dxyz_CM_s(:,i,inp) = matmul(AAT,dxyz_CM_b(:,i,inp))
!         dxyz_CM_s(:,i,inp) = matmul(AAT,dxyz_CM_s(:,i,inp))

         !tri_bar(:,i,inp) = pos_cm(:,inp) + dxyz_CM_s(:,i,inp) 
            ! KZ: above commented: don't store tri_centroids as COM-relative for now due to loss of sig figs.
         call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face,maxnf,maxnv) 

          !-- velocity
         omega_s(:,inp) = matmul(AAT,omega_b(:,inp))
         call cross(om_dCM(:), omega_s(:,inp), dxyz_CM_s(:,i,inp))
         vel_tri(:,i,inp) = vel_CM(:,inp) + om_dCM(:) 
      end do
   end do


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
    dname = trim("omega_b")
    call hdf_read_2d(omega_b,shape(omega_b),dname)
    dname = trim("om_b_sqr")
    call hdf_read_2d(om_b_sqr,shape(om_b_sqr),dname)


    dname = trim("u_tot")
    call hdf_read_2d(u_tot,shape(u_tot),dname)
    dname = trim("r_x_u_tot")
    call hdf_read_2d(r_x_u_tot,shape(r_x_u_tot),dname)


    call mpi_bcast(pos_CM,size(pos_cm),mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(vel_CM,size(vel_cm),mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(quat,     size(quat),     mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(quat_dot, size(quat_dot), mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(omega_b,  size(omega_b),  mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(om_b_sqr, size(om_b_sqr), mpi_double,0,mpi_comm_world,ierr)

    call mpi_bcast(u_tot,    size(u_tot),    mpi_double,0,mpi_comm_world,ierr)
    call mpi_bcast(r_x_u_tot,size(r_x_u_tot),mpi_double,0,mpi_comm_world,ierr)

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

    ! 2D
    dataset = trim("pos_CM")
    call hdf_write_2d(pos_CM,(/3,Nparticle/),dataset)
    dataset = trim("vel_CM")
    call hdf_write_2d(vel_CM,(/3,Nparticle/),dataset)
    dataset = trim("quat")
    call hdf_write_2d(quat,(/4,Nparticle/),dataset)
    dataset = trim("quat_dot")
    call hdf_write_2d(quat_dot,(/4,Nparticle/),dataset)
    dataset = trim("omega_b")
    call hdf_write_2d(omega_b,(/3,Nparticle/),dataset)
    dataset = trim("om_b_sqr")
    call hdf_write_2d(om_b_sqr,(/3,Nparticle/),dataset)

    dataset = trim("u_tot")
    call hdf_write_2d(u_tot,(/3,Nparticle/),dataset)
    dataset = trim("r_x_u_tot")
    call hdf_write_2d(r_x_u_tot,(/3,Nparticle/),dataset)

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
    write(*,*) 'Ave edge/dx: ',    (sum(dist(:,1))/float(maxne)) * ddx3
    write(*,*) 'Max edge/dx: ',    maxval(dist(:,1)) * ddx3
    write(*,*) 'Min edge/dx: ',    minval(dist(:,1)) * ddx3
  end if
end subroutine

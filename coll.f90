module coll_mod
use param
use mls_param
use mpih

real :: fp(3,Nparticle)
real :: tp(3,Nparticle)

real :: fp_m1(3,Nparticle)
real :: tp_m1(3,Nparticle)

real :: dt_p
real :: short_dist

logical :: is_coll(Nparticle)
logical :: coll_check

! -------------------------------------------------------------
! Updating position and orientation
! -------------------------------------------------------------
real, dimension(3,Nparticle) :: coll_vel,      coll_vel_m1
real, dimension(3,Nparticle) :: coll_om_b,     coll_om_b_m1
real, dimension(4,Nparticle) :: coll_quat_dot, coll_quat_dot_m1
 
 

contains



subroutine collision 
implicit none
integer :: i,j
real    :: p1(3),p2(3),dij(3),dist,dc,x1(3),x2(3)
real    :: f1(3), f2(3), t1(3), t2(3)
real    :: AA1(3,3)
real    :: AA2(3,3)


if (coll_check.eq. .false.) then
  fp = 0.d0
  tp = 0.d0
endif


dc = 1./16.

do i = 1,Nparticle
   do j = i+1,Nparticle

      call closest_distance(i,j,p1,p2,x1,x2,dist,AA1,AA2)

      if (dist .lt. dc) then

        is_coll(i) = .true.
        is_coll(j) = .true.

        dij = p2 - p1
        f1 = -1e4 * ( ((dist-dc)/dc)**2 ) * dij/norm2(dij)
        f2 = -f1

        call cross(t1, p1-x1, f1)
        call cross(t2, p2-x2, f2)

        t1 = matmul(AA1, t1)
        t2 = matmul(AA2, t2)

      if (coll_check.eq. .false.) then
        fp(1:3,i) = fp(1:3,i) + f1
        fp(1:3,j) = fp(1:3,j) + f2
        tp(1:3,i) = tp(1:3,i) + t1
        tp(1:3,j) = tp(1:3,j) + t2
      endif

      endif

      short_dist = dist

enddo
enddo

end subroutine


subroutine closest_distance(i,j,l1_short,l2_short,pos1,pos2,dist,AA1,AA2)
implicit none
integer :: i,j,np, inp1, inp2
real :: AA1(3,3), AAT1(3,3)
real :: AA2(3,3), AAT2(3,3)
real :: y1(3,3), y2(3,3)
real :: part_1_p1(3,3), part_1_p2(3,3)
real :: part_2_p1(3,3), part_2_p2(3,3)
real :: P1(3), P2(3), Q1(3),Q2(3),u(3),v(3),w0(3)
real :: d1(3),d2(3),pos1(3),pos2(3)
real :: l1_short(3), l2_short(3)
real :: a,b,c,d,e,den,s,t
real :: mydist, dist


call calc_rot_matrix(quat(:,i),AA1)
call calc_rot_matrix(quat(:,j),AA2)

AAT1 = transpose(AA1) 
AAT2 = transpose(AA2) 

call part_body_frame(y1,y2)

pos1 = pos_cm(1:3, i)
pos2 = closest_particle(pos_cm(1:3, i), pos_cm(1:3, j))

do np = 1,3
  part_1_p1(:,np) = pos_cm(1:3, i) + matmul(AAT1, y1(:,np))
  part_1_p2(:,np) = pos_cm(1:3, i) + matmul(AAT1, y2(:,np))

  part_2_p1(:,np) = pos2 + matmul(AAT2, y1(:,np))
  part_2_p2(:,np) = pos2 + matmul(AAT2, y2(:,np))
enddo

dist = 1.e6

! short dist calc
do inp1=1,3
  do inp2=1,3
    P1 = part_1_p1(:,inp1) 
    P2 = part_2_p1(:,inp2) 

    Q1 = part_1_p2(:,inp1)
    Q2 = part_2_p2(:,inp2)

    u = Q1-P1
    v = Q2-P2
    w0 = P1-P2

    a = dot_product(u, u)
    b = dot_product(u, v)
    c = dot_product(v, v)
    d = dot_product(w0, u)
    e = dot_product(w0, v)

    den = a*c - b**2 + 1e-10
    s = (b*e - c*d) / den
    t = (a*e - b*d) / den 

    s = min(s,1.);
    s = max(s,0.);
    t = min(t,1.);
    t = max(t,0.);

    d1 = P1 + s*u;
    d2 = P2 + t*v;

    mydist = norm2(d1 - d2)
    if (mydist.lt.dist) then
      l1_short = d1
      l2_short = d2
      dist = mydist
    endif

  enddo
enddo

dist = dist - 2* 0.2125855

!if(myid.eq.0) print *, dist 

!   ! -------------------------------------------------
!   ! brute force - just for verification
!   mydist = 1e6
!   do inp1=1,maxnf
!    do inp2=1,maxnf
!       mydist = min(mydist, norm2(tri_bar(1:3,inp1,1) -  tri_bar(1:3,inp2,2)))
!       enddo
!       enddo
!if(myid.eq.0) print *, "brute force", mydist

end subroutine


! y1(1,i) 1st x-value of particle i
! y2(1,i) 2nd x-value of particle i
subroutine part_body_frame(y1,y2) 
implicit none
real    :: y1(3,3), y2(3,3)
real    :: r

r = .544279548

y1 = 0.
y2 = 0.


! particle 1
y1(1,1) = -r
y2(1,1) =  r
! particle 2
y1(1,2) = -r
y2(1,2) = -r
y2(2,2) = 2.*r
! particle 3
y1(1,3) = r
y2(1,3) = r
y2(3,3) = -2.*r

! moviing centroid such that volumetric centre is at (0,0,0)
y1(2,:) = y1(2,:) - 0.2125855
y2(2,:) = y2(2,:) - 0.2125855

y1(3,:) = y1(3,:) + 0.2125855
y2(3,:) = y2(3,:) + 0.2125855
end subroutine


function closest_particle(p1,p2) result(p2_out)
integer :: i
real    :: p1(3), p2(3), p2_out(3)
real    :: xb(3), yb(3), zb(3)
real    :: y(3,27)
real    :: dist(27)

  xb = 0.; 
  yb = 0.;
  zb = 0.;

  xb(1) = xlen
  yb(2) = ylen
  zb(3) = zlen

  !mid plane
  y(1:3,1) = p2 
  y(1:3,2) = p2 - xb
  y(1:3,3) = p2 + xb
  y(1:3,4) = p2 - yb
  y(1:3,5) = p2 + yb
  y(1:3,6) = p2 - xb + yb
  y(1:3,7) = p2 - xb - yb
  y(1:3,8) = p2 + xb + yb
  y(1:3,9) = p2 + xb - yb

  !lower plane
  y(1:3,10) = p2 - zb 
  y(1:3,11) = p2 - xb - zb
  y(1:3,12) = p2 + xb - zb
  y(1:3,13) = p2 - yb - zb
  y(1:3,14) = p2 + yb - zb
  y(1:3,15) = p2 - xb + yb - zb
  y(1:3,16) = p2 - xb - yb - zb
  y(1:3,17) = p2 + xb + yb - zb
  y(1:3,18) = p2 + xb - yb - zb

  !upper plane
  y(1:3,19) = p2 + zb
  y(1:3,20) = p2 - xb + zb
  y(1:3,21) = p2 + xb + zb
  y(1:3,22) = p2 - yb + zb
  y(1:3,23) = p2 + yb + zb
  y(1:3,24) = p2 - xb + yb + zb
  y(1:3,25) = p2 - xb - yb + zb
  y(1:3,26) = p2 + xb + yb + zb
  y(1:3,27) = p2 + xb - yb + zb

  do i = 1,27
   dist(i) = norm2( p1 - y(1:3,i) )
  enddo

  ! min loc
  i = minloc(dist, dim=1)

  p2_out = y(1:3,i)

end function

! write collision continuation data
subroutine continua_collision
implicit none
character(30) dataset
integer, dimension(2) :: dims

dims(1)=n1m
dims(2)=n2m


if (myid.eq.0) then 
  ! 2D
  dataset = trim("fp")
  call hdf_write_2d(fp,(/3,Nparticle/),dataset)
  dataset = trim("tp")
  call hdf_write_2d(tp,(/3,Nparticle/),dataset)

  !-- RigidAuxRoutines
  dataset = trim("coll_vel")
  call hdf_write_2d(coll_vel,(/3,Nparticle/),dataset)
  dataset = trim("coll_om_b")
  call hdf_write_2d(coll_om_b,(/3,Nparticle/),dataset)
  dataset = trim("coll_quat_dot")
  call hdf_write_2d(coll_quat_dot,(/3,Nparticle/),dataset)
end if
end subroutine


! read collision continuation data

subroutine import_collision
implicit none
character(200) :: dname 

  dname = trim("fp")
  call hdf_read_2d(fp,shape(fp),dname)
  dname = trim("tp")
  call hdf_read_2d(tp,shape(tp),dname)
  !-- RigidAuxRoutines
  dname = trim("coll_vel")
  call hdf_read_2d(coll_vel,shape(coll_vel),dname)
  dname = trim("coll_om_b")
  call hdf_read_2d(coll_om_b,shape(coll_om_b),dname)
  dname = trim("coll_quat_dot")
  call hdf_read_2d(coll_quat_dot,shape(coll_quat_dot),dname)

  ! broadcast vars to all ranks
  call mpi_bcast(fp,   size(fp),   mpi_double,0,mpi_comm_world,ierr)
  call mpi_bcast(tp,   size(tp),   mpi_double,0,mpi_comm_world,ierr)
  !-- RigidAuxRoutines
  call mpi_bcast(coll_vel,       size(coll_vel),       mpi_double,0,mpi_comm_world,ierr)
  call mpi_bcast(coll_om_b,      size(coll_om_b),      mpi_double,0,mpi_comm_world,ierr)
  call mpi_bcast(coll_quat_dot,  size(coll_quat_dot),  mpi_double,0,mpi_comm_world,ierr)
end subroutine

end module

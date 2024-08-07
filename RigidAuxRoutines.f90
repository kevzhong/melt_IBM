subroutine update_part_pos()
use param
use mls_param
use mpih
!use coll_mod ! KZ: no collisions for now
implicit none
integer :: inp,i
real,dimension(3,3) :: AA, AAT

! set flag to false
!is_coll = .false.

! check for collision
!coll_check = .false.
!call collision
!coll_check = .false.


! KZ: ---- pre-compute COM-relative, world-aligned coordinates
do inp=1,Nparticle
  do i = 1,maxnv
    if (.not. isGhostVert(i,inp) ) then
      dxyzv_s(:,i,inp) =   xyzv(:,i,inp) - pos_CM(:,inp)
    endif
  enddo
enddo

! COM-relative of centroids
!call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) 
do inp=1,Nparticle
  do i = 1,maxnf
    if (.not. isGhostFace(i,inp) ) then
      dxyz_s(:,i,inp) = tri_bar(:,i,inp) - pos_CM(:,inp)
    endif
  enddo
enddo

call update_regular
!call update_substepping

do inp=1,nParticle
    ! KZ: need to re-orient inertia tensor / principal axes post-rotation
  !AAT = transpose(AA)
  !tail_head(:,inp) = AAT(:,3) 

 if (pos_cm(1,inp).lt.0.)   pos_cm(1,inp) = pos_cm(1,inp) + xlen
 if (pos_cm(1,inp).gt.xlen) pos_cm(1,inp) = pos_cm(1,inp) - xlen

 if (pos_cm(2,inp).lt.0.)   pos_cm(2,inp) = pos_cm(2,inp) + ylen
 if (pos_cm(2,inp).gt.ylen) pos_cm(2,inp) = pos_cm(2,inp) - ylen

 if (pos_cm(3,inp).lt.0.)   pos_cm(3,inp) = pos_cm(3,inp) + zlen
 if (pos_cm(3,inp).gt.zlen) pos_cm(3,inp) = pos_cm(3,inp) - zlen
enddo

! KZ replace this, and make use of updated xyzv_b above
call update_xyz

end subroutine

subroutine update_regular
use param
use mls_param
!use coll_mod
implicit none
integer :: i,inp
!real,dimension(3,3)     :: AA, AAT
real,dimension(3,2)     :: bbox_inds
real,dimension(Nparticle) :: error
real,dimension(3,Nparticle) :: vel_m1,pos_m1,pos_k,om_m1,acm_m1

do inp=1,Nparticle

  !   vel_m1(:,inp)       = vel_cm(:,inp)
  !   pos_m1(:,inp)       = pos_cm(:,inp)
  !   om_m1(:,inp)        = omega_c(:,inp)
  !   acm_m1(:,inp)       = a_cm(:,inp)
  !      ! cutvol
  !   u_tot_m1(:,inp)     = u_tot(:,inp)
  !   r_x_u_tot_m1(:,inp) = r_x_u_tot(:,inp)

  !   ! compute int(u)dV and int(r x u)dV
  !   u_tot(:,inp)     = 0.
  !   r_x_u_tot(:,inp) = 0.

  !   VOFx(:,:,:) = 1.
  !   VOFy(:,:,:) = 1.
  !   VOFz(:,:,:) = 1.

  !  call get_bbox_inds(bbox_inds,inp)
  !  call convex_hull_q1(bbox_inds,inp)
  !  call convex_hull_q2(bbox_inds,inp)
  !  call convex_hull_q3(bbox_inds,inp)

  !   call newton_euler(fpxyz(:,inp),       ftxyz(:,inp),                  &
  !   vel_CM(:,inp),      vel_m1(:,inp),                 &
  !   pos_CM(:,inp),      pos_m1(:,inp),                 &
  !   omega_c(:,inp),     om_m1(:,inp),                  &
  !   quat(:,inp),                                       &
  !   InertTensor(:,:,inp),                              &
  !   u_tot(:,inp),       u_tot_m1(:,inp),               &
  !   r_x_u_tot(:,inp),   r_x_u_tot_m1(:,inp),           &
  !   inp,                 a_CM(:,inp),                  &
  !   alpha_b(:,inp),      acm_m1(:,inp),                &
  !   Volume(inp) ) ! Volume pre-factor for IBM force term


  vel_m1(:,inp)       = vel_cm(:,inp)
  pos_m1(:,inp)       = pos_cm(:,inp)
  om_m1(:,inp)        = omega_c(:,inp)

  call NE_MLS(int_prn_dA(1:3,inp),        int_tau_dA(1:3,inp),        &
              int_prn_dA_m1(1:3,inp),     int_tau_dA_m1(1:3,inp),     &
              int_r_x_prn_dA(1:3,inp),    int_r_x_tau_dA(1:3,inp),    &
              int_r_x_prn_dA_m1(1:3,inp), int_r_x_tau_dA_m1(1:3,inp), &
              vel_CM(:,inp),              vel_m1(:,inp),              &
              pos_CM(:,inp),              pos_m1(:,inp),              &
              omega_c(:,inp),             om_m1(:,inp),               &
              InertTensor(:,:,inp),                                   &
              inp,                                                    &
              Volume(inp) )
enddo

end subroutine


! subroutine newton_euler(For_tot,  torq_surf,   &
!                         vel_CM,   vel_cmm1,    &
!                         pos_CM,   pos_cmm1,    &
!                         omega_c,  omega_c_m1,  &
!                         quat,     quat_m1,     &
!                         quat_dot, quat_dot_m1, &
!                         I_ij,                  &
!                         u_tot,    u_tot_m1,    &
!                         r_x_u,    r_x_u_m1,    &
!                         inp,      a_CM,        &
!                         alpha_b,  acm_m1,      &
!                         Volume)

subroutine NE_MLS(int_prn_dA,         int_tau_dA,        &
                  int_prn_dA_m1,      int_tau_dA_m1,     &
                  int_r_x_prn_dA,     int_r_x_tau_dA,    &
                  int_r_x_prn_dA_m1,  int_r_x_tau_dA_m1, &
                  vel_CM,   vel_cmm1,                    &
                  pos_CM,   pos_cmm1,                    &
                  omega_c,  omega_c_m1,                  &
                  I_ij,                                  &
                  inp,                                   &
                  Volume)

! Given total force and total torque acting on the body, it evolves
! the newton equation for center of mass and
! the Euler equation for rotation (quaternion notation Allen Tildesley pag 103)
! Breugem 2012 (JCP)
use param, only: pi,ntime,al,ga,ro,dt, ismaster
use mls_param,only: dens_ratio
use mpih
implicit none

real,dimension(3) :: int_prn_dA, int_tau_dA ! total forces acting on the rigid body
real,dimension(3) :: int_prn_dA_m1, int_tau_dA_m1 ! total forces acting on the rigid body
real,dimension(3) :: int_r_x_prn_dA, int_r_x_tau_dA ! total torques acting on the rigid body
real,dimension(3) :: int_r_x_prn_dA_m1, int_r_x_tau_dA_m1 ! total torques acting on the rigid body
real,dimension(3) :: vel_CM        ! velocity of the center of mass
real,dimension(3) :: pos_CM        ! position of the center of mass
real,dimension(4)   :: quat        ! quaternions
real,dimension(3)   :: omega_c     ! angular vel body frame
real,dimension(3,3) :: I_ij    ! rotation matrix "body  = AA  space "
real ,dimension(3)  :: e_z
! local variables
real,dimension(3) :: vel_CMm1            ! velocity of the center of mass
real,dimension(3) :: pos_CMm1            ! position of the center of mass
real,dimension(3) :: omega_c_m1          ! angular vel body frame
integer :: inp
real :: Volume, angle_buffer, magOM, pre_fac

! Iterative inertia tensor update (Ardekani et al. 2016)
real :: tolerance, Iresidual ! Iteration tolerance for the lhs inertia tensor
integer :: n_iter
real, dimension(3,3) :: I_prev, I_LHS, AA, buffer_3x3, I_inv
real, dimension(3)  :: om_buffer
! For LAPACK
integer :: LWORK = 3
real, dimension(3) :: WORK
integer :: IPIV(3), INFO

! Coefficienti
! RK3 coeff here, alm(ns), are two times alpha(k) in Rai & Moin (JCP) 1991
! In Rai and moin:     HERE               BREUGEM 2012 (JCP)
! alpha(1) = 4/15      alpha(1) = 8/15    a(1)+b(1) = 8/15
! alpha(2) = 1/15      alpha(2) = 2/15    a(2)+b(2) = 2/15
! alpha(3) = 1/6       alpha(3) = 2/6     a(3)+b(3) = 2/6

! -------------------------------------
!              Translation
! -------------------------------------

e_z = 0.; e_z(3) = -1.0

pre_fac = 1.0 / Volume / dens_ratio
! translation
!a_CM = - pre_fac * for_tot &
!      + e_z / dens_ratio &
!      + (u_tot - u_tot_m1) / (dens_ratio*dt)


!vel_CM = vel_CMm1 + al * dt * pre_fac * ( -int_prn_dA + int_tau_dA )         & ! Hydrodynamic loads
!                  + al * dt * (1.0 - ( 1.0 / dens_ratio ) ) * e_z        ! Gravity term

vel_CM = vel_CMm1 +  ga * dt * pre_fac * ( -int_prn_dA + int_tau_dA ) &
                  +  ro * dt * pre_fac * ( -int_prn_dA_m1 + int_tau_dA_m1 ) &
                   + al * dt * (1.0 - ( 1.0 / dens_ratio ) ) * e_z        ! Gravity term



!vel_CM=0.0d0
pos_CM = pos_CMm1 + 0.5 * al * dt * ( vel_CM + vel_CMm1 )


!if (ismaster) then
!  write(*,*) "POS_CM: ", pos_CM
!endif
! ------------------------------------- 
!               Rotation
! ------------------------------------- 


!dr_x_u_b = r_x_u - r_x_u_m1

!omega_b = omega_b_m1 + matmul(I_inv, -(dt/dens_ratio)*torq_b + & ! IBM force term
!                                 dr_x_u_b)  & ! Torque impulse term
!                              + dt*al*matmul(I_inv2, omega_b_m1_squared) ! Non-inertial reference frame correction
I_LHS(:,:) = I_ij
I_inv = I_LHS
call dgetrf(3, 3, I_inv, 3, IPIV, INFO) ! Compute the LU factorization of I_ij
if (info /= 0) then
print *, "Error in dgetrf"
stop
end if
call dgetri(3, I_inv, 3, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization

if (info /= 0) then
print *, "Error in dgetri"
stop
end if

!tolerance = 1e-13
!Iresidual = 1.0
!n_iter = 1
!do while(  (n_iter .le. 10000 ) .or.  (Iresidual .gt. tolerance)   )


do n_iter = 1,10 ! Iterative update of inertia tensor post-rotation

  !omega_c =  matmul(I_inv, matmul(I_ij,omega_c_m1) &
  !                  + al * dt * ( -int_r_x_prn_dA + int_r_x_tau_dA ) )

   omega_c =  matmul(I_inv, matmul(I_ij,omega_c_m1) &
                     +  ga * dt * ( -int_r_x_prn_dA + int_r_x_tau_dA )  &
                     +  ro * dt * ( -int_r_x_prn_dA_m1 + int_r_x_tau_dA_m1 ) ) 


  ! No rotation
  !omega_c = 0.0

  ! Ardekani et al. (2016) iterative update method
  om_buffer = 0.5* (omega_c + omega_c_m1 )
  magOM = norm2(om_buffer)

  angle_buffer = magOM * al* dt

  if (magOM .lt. EPSILON(1.0d0) ) then ! Axis vector from rotation is zero-vector, avoid div-by-zero
    om_buffer(1:3) = 0.0
  if (ismaster) then
    write(*,*) "Div. by zero avoided for rotation axis!"
  endif
  else
    om_buffer(1:3) = om_buffer(1:3) / magOM
  endif

  ! Set provisional quarternion from this new axis-angle representation
  quat(1) = cos( 0.5 * angle_buffer )
  quat(2) = sin( 0.5 * angle_buffer ) * om_buffer(1)
  quat(3) = sin( 0.5 * angle_buffer ) * om_buffer(2)
  quat(4) = sin( 0.5 * angle_buffer ) * om_buffer(3)
  quat = quat / norm2(quat)

  call calc_rot_matrix(quat,AA)

  !----- New iteratively-updated estimate of inertia tensor -------
  I_prev = I_LHS

  buffer_3x3 = matmul(I_ij,transpose(AA))
  I_LHS = matmul(AA, buffer_3x3)

  !buffer_3x3 = matmul(I_LHS,transpose(AA))
  !I_LHS = matmul(AA, buffer_3x3)

  !buffer_3x3 = matmul(I_LHS,AA)
  !I_LHS = matmul(transpose(AA), buffer_3x3)

  Iresidual = maxval (abs (I_LHS - I_prev) )

  I_inv = I_LHS
  call dgetrf(3, 3, I_inv, 3, IPIV, INFO) ! Compute the LU factorization of I_ij
  if (info /= 0) then
  print *, "Error in dgetrf"
  stop
  end if
  call dgetri(3, I_inv, 3, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization

  if (info /= 0) then
  print *, "Error in dgetri"
  stop
  end if

!n_iter = n_iter + 1
enddo

!vel_CM = vel_CMm1 - al * dt * ( -int_prn_dA + int_tau_dA )         & ! Hydrodynamic loads
!                  + al * dt * (1.0 - ( 1.0 / dens_ratio ) ) * e_z        ! Gravity term



! if (ismaster) then
! write(*,*) "Iresidual: ", Iresidual
! write(*,*) "int_prn_dA", int_prn_dA
! write(*,*) "int_tau_dA", int_tau_dA
! endif

! if(ismaster) then
!   open(112,file='flowmov/mlsLoads.txt',status='unknown', position='append')
!         write(112,'(40E17.5)') int_prn_dA, int_tau_dA
!   close(112)
! end if



end subroutine NE_MLS


! subroutine newton_euler(For_tot,  torq_surf,   &
!                         vel_CM,   vel_cmm1,    &
!                         pos_CM,   pos_cmm1,    &
!                         omega_c,  omega_c_m1,  &
!                         quat,                  &
!                         I_ij,                  &
!                         u_tot,    u_tot_m1,    &
!                         r_x_u,    r_x_u_m1,    &
!                         inp,      a_CM,        &
!                         alpha_b,  acm_m1,      &
!                         Volume)



! ! Given total force and total torque acting on the body, it evolves
! ! the newton equation for center of mass and
! ! the Euler equation for rotation (quaternion notation Allen Tildesley pag 103)
! ! Breugem 2012 (JCP)
!   use param, only: pi,ntime,al,ga,ro,dt, ismaster
!   use mls_param,only: dens_ratio
!   use mpih
!   implicit none

!   real,dimension(3) :: For_tot       ! total force actin on the rigid body
!   real,dimension(3) :: a_CM          ! acceleration of the center of mass
!   real,dimension(3) :: vel_CM        ! velocity of the center of mass
!   real,dimension(3) :: pos_CM        ! position of the center of mass
!   real              :: pre_fac       ! prefactor
!   real,dimension(3) :: torq_b        ! torque actin on the rigid body
!   real,dimension(3) :: dr_x_u_b      ! torque actin on the rigid body
!                                      ! (with respect to mass center and represented
!                                      !  in space frame)
!   real,dimension(4)   :: quat        ! quaternions
!   real,dimension(3)   :: alpha_b     ! angular acc body frame
!   real,dimension(3)   :: omega_c     ! angular vel body frame
!   real,dimension(3,3) :: I_ij    ! rotation matrix "body  = AA  space "
!   real,dimension(3)   :: u_tot       ! <u>_V over the particle
!   real,dimension(3)   :: torq_surf, r_x_u, r_x_u_m1
!   real ,dimension(3)  :: e_z
! ! local variables
!   real,dimension(3) :: acm_m1            ! velocity of the center of mass
!   real,dimension(3) :: vel_CMm1            ! velocity of the center of mass
!   real,dimension(3) :: pos_CMm1            ! position of the center of mass
!   real,dimension(3) :: omega_c_m1          ! angular vel body frame
!   real,dimension(3) :: u_tot_m1            ! <u>_V over the particle
!   integer :: inp
!   real :: Volume, angle_buffer, magOM

!   ! Iterative inertia tensor update (Ardekani et al. 2016)
!   real :: tolerance, Iresidual ! Iteration tolerance for the lhs inertia tensor
!   integer :: n_iter
!   real, dimension(3,3) :: I_prev, I_LHS, AA, buffer_3x3, I_inv
!   real, dimension(3)  :: om_buffer
!   ! For LAPACK
!   integer :: LWORK = 3
!   real, dimension(3) :: WORK
!   integer :: IPIV(3), INFO

! ! Coefficienti
! ! RK3 coeff here, alm(ns), are two times alpha(k) in Rai & Moin (JCP) 1991
! ! In Rai and moin:     HERE               BREUGEM 2012 (JCP)
! ! alpha(1) = 4/15      alpha(1) = 8/15    a(1)+b(1) = 8/15
! ! alpha(2) = 1/15      alpha(2) = 2/15    a(2)+b(2) = 2/15
! ! alpha(3) = 1/6       alpha(3) = 2/6     a(3)+b(3) = 2/6

! ! -------------------------------------
! !              Translation
! ! -------------------------------------

!   e_z = 0.; e_z(3) = -1.0
  
!   pre_fac = 1.0 / Volume / dens_ratio
! ! translation
!   !a_CM = - pre_fac * for_tot &
!    !      + e_z / dens_ratio &
!    !      + (u_tot - u_tot_m1) / (dens_ratio*dt)


!   vel_CM = vel_CMm1 - dt * pre_fac * For_tot           & ! IBM force term
!                    + dt * al * (1.0 - ( 1.0 / dens_ratio ) ) * e_z       & ! Gravity term
!                    + (u_tot - u_tot_m1) / dens_ratio    ! Impulse term, already normalised on VOF-computed particle volume


! !vel_CM=0.0d0
! pos_CM = pos_CMm1 + 0.5 * al * dt * ( vel_CM + vel_CMm1 )


!   !if (ismaster) then
!   !  write(*,*) "POS_CM: ", pos_CM
!   !endif
!   ! ------------------------------------- 
!   !               Rotation
!   ! ------------------------------------- 
  

!   dr_x_u_b = r_x_u - r_x_u_m1

!   !omega_b = omega_b_m1 + matmul(I_inv, -(dt/dens_ratio)*torq_b + & ! IBM force term
!   !                                 dr_x_u_b)  & ! Torque impulse term
!   !                              + dt*al*matmul(I_inv2, omega_b_m1_squared) ! Non-inertial reference frame correction
!   I_LHS(:,:) = I_ij
!   I_inv = I_LHS
!   call dgetrf(3, 3, I_inv, 3, IPIV, INFO) ! Compute the LU factorization of I_ij
!   if (info /= 0) then
!     print *, "Error in dgetrf"
!     stop
!   end if
!   call dgetri(3, I_inv, 3, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization

!   if (info /= 0) then
!     print *, "Error in dgetri"
!     stop
!   end if
  
!   !tolerance = 1e-13
!   !Iresidual = 1.0
!   !n_iter = 1
!   !do while(  (n_iter .le. 10000 ) .or.  (Iresidual .gt. tolerance)   )


!   do n_iter = 1,10 ! Iterative update of inertia tensor post-rotation

!     omega_c =  + matmul(I_inv, matmul(I_ij,omega_c_m1) &
!     - dt *torq_surf  &   ! IBM force term:  density ratio pre-factors are not needed since absorbed into inertia tensor
!     +  dr_x_u_b  )   ! Torque impulse term: density ratio pre-factors are not needed since absorbed into inertia tensor

!    ! omega_c =  + matmul(I_inv, matmul(I_ij,omega_c_m1) ) &
!    !            - matmul(I_inv, - dt *torq_surf  ) &   ! IBM force term:  density ratio pre-factors are not needed since absorbed into inertia tensor
!    !            + matmul(I_inv,   dr_x_u_b  )   ! Torque impulse term: density ratio pre-factors are not needed since absorbed into inertia tensor
  
!   ! No rotation
!   !omega_c = 0.0

!  ! Ardekani et al. (2016) iterative update method
!   om_buffer = 0.5* (omega_c + omega_c_m1 )
!   magOM = norm2(om_buffer)

!   angle_buffer = magOM * al* dt

!   if (magOM .lt. 1.0e-16 ) then ! Axis vector from rotation is zero-vector, avoid div-by-zero
!     om_buffer(1:3) = 0.0
!     if (ismaster) then
!       write(*,*) "Div. by zero avoided for rotation axis!"
!     endif
!   else
!     om_buffer(1:3) = om_buffer(1:3) / magOM
!   endif

!   ! Set provisional quarternion from this new axis-angle representation
!   quat(1) = cos( 0.5 * angle_buffer )
!   quat(2) = sin( 0.5 * angle_buffer ) * om_buffer(1)
!   quat(3) = sin( 0.5 * angle_buffer ) * om_buffer(2)
!   quat(4) = sin( 0.5 * angle_buffer ) * om_buffer(3)
!   quat = quat / norm2(quat)

!   call calc_rot_matrix(quat,AA)

!   !----- New iteratively-updated estimate of inertia tensor -------
!   I_prev = I_LHS

!   buffer_3x3 = matmul(I_ij,transpose(AA))
!   I_LHS = matmul(AA, buffer_3x3)

!   !buffer_3x3 = matmul(I_LHS,transpose(AA))
!   !I_LHS = matmul(AA, buffer_3x3)

!   !buffer_3x3 = matmul(I_LHS,AA)
!   !I_LHS = matmul(transpose(AA), buffer_3x3)

!   Iresidual = maxval (abs (I_LHS - I_prev) )

!   I_inv = I_LHS
!   call dgetrf(3, 3, I_inv, 3, IPIV, INFO) ! Compute the LU factorization of I_ij
!   if (info /= 0) then
!     print *, "Error in dgetrf"
!     stop
!   end if
!   call dgetri(3, I_inv, 3, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization

!   if (info /= 0) then
!     print *, "Error in dgetri"
!     stop
!   end if
  
!   !n_iter = n_iter + 1
!   enddo


!   !omega_c =  + matmul(I_inv, matmul(I_ij,omega_c_m1) &
!   !- (dt/dens_ratio)*torq_surf  &   ! IBM force term
!   !+ (1.0 / dens_ratio) * dr_x_u_b  )   ! Torque impulse term
!                              ! No rotation


!   if (ismaster) then
!     write(*,*) "Iresidual: ", Iresidual
!     !write(*,*) "-------------------------"
!     !write(*,*) -1.0*dt * pre_fac * For_tot, " IBM translation"
!     !write(*,*) (u_tot - u_tot_m1) / dens_ratio, " Impulse translation"
!     !write(*,*) "-------------------------"
!     !write(*,*) matmul(I_inv, matmul(I_ij,omega_c_m1) ), "I_ij x om_m1"
!     !rite(*,*) -1.0*matmul(I_inv, - dt *torq_surf  ), "IBM rotation"
!     !write(*,*) matmul(I_inv,   dr_x_u_b  ), "Rotation impulse"
!     !write(*,*) "-------------------------"
!   endif



! end subroutine newton_euler

!===========================================================================
! collision
!===========================================================================

! subroutine update_substepping
! use param
! use mls_param
! use coll_mod
! implicit none
! integer :: i,inp
! real,dimension(3,3)          :: AA, AAT
! real,dimension(Nparticle)    :: error
! integer,dimension(Nparticle) :: step
! real,dimension(3,Nparticle)  :: vel_m1,pos_m1,pos_k,om_m1
! real :: pos_new(3)
! integer :: iter,nIter



! nIter = 50
! dt_p = dt / dble(nIter)

! do iter=1,nIter

!   do inp=1,Nparticle
!       quat_m1(:,inp) = quat(:,inp)
!       vel_m1(:,inp)  = vel_cm(:,inp)
!       pos_m1(:,inp)  = pos_cm(:,inp)
!       om_m1(:,inp)   = omega_c(:,inp)
!       !-- force & torque
!       fp_m1(:,inp)   = fp(:,inp)
!       tp_m1(:,inp)   = tp(:,inp)
!       !-- collision vars
!       coll_vel_m1(:,inp)      = coll_vel(:,inp)
!       coll_om_b_m1(:,inp)     = coll_om_b(:,inp)
!       coll_quat_dot_m1(:,inp) = coll_quat_dot(:,inp)

!     if (is_coll(inp).eqv. .false.) then
!       !-- reset vars
!       coll_vel(:,inp)      = 0.d0
!       coll_om_b(:,inp)     = 0.d0
!       coll_quat_dot(:,inp) = 0.d0
!     endif

!   enddo

!   error = 1e6
!   step  = 0



!   do while ( maxval(error) .gt. zlen/n3m .and. maxval(step) .lt. 10)

!        do inp = 1,Nparticle
!            if (is_coll(inp).eqv. .true. ) then

!            ! update fp and tp
!            call collision

!            ! get latest rotation matrix
!            call calc_rot_matrix(quat(:,inp),AA)
!            AAT = transpose(AA)

!            ! store pos
!            pos_k(:,inp) = pos_cm(1:3,inp) + matmul(AAT, dxyz_CM_b(:,1,inp))

!            ! update vel, pos, omega and quat
!            call newton_coll(  vel_CM(:,inp), vel_m1(:,inp),                          &
!                               pos_CM(:,inp), pos_m1(:,inp),                          &
!                               omega_c(:,inp), om_m1(:,inp),                          &
!                               quat(:,inp), quat_m1(:,inp),                           & 
!                               AA(:,:),                                               &    
!                               coll_vel(:,inp), coll_vel_m1(:,inp),                   &
!                               coll_om_b(:,inp),                                      &
!                               coll_quat_dot(:,inp), coll_quat_dot_m1(:,inp),         &
!                               fp(:,inp), fp_m1(:,inp), tp(:,inp), tp_m1(:,inp),      & 
!                               dt_p) 


!            ! update rotation matrix
!            call calc_rot_matrix(quat(:,inp),AA)
!            AAT = transpose(AA)

!            pos_new = pos_cm(1:3,inp) + matmul(AAT, dxyz_CM_b(:,1,inp))

!            error(inp) = norm2(pos_k(:,inp) - pos_new)
!            step(inp) = step(inp) + 1

!            else

!              error(inp) = 0.d0

!            endif
!            enddo
!   enddo
! enddo
! end subroutine


! subroutine newton_coll ( vel_CM, vel_m1,                    &
!                          pos_CM, pos_m1,                    &
!                          omega_c, omega_c_m1,               &
!                          quat, quat_m1,                     &
!                          AA,                                &
!                          coll_vel, coll_vel_m1,             &
!                          coll_om_b,                         &
!                          coll_quat_dot, coll_quat_dot_m1,   &
!                          fp, fp_m1, tp, tp_m1,              &
!                          dt )

! ! Given total force and total torque acting on the body, it evolves
! ! the newton equation for center of mass and
! ! the Euler equation for rotation (quaternion notation Allen Tildesley pag 103)
! ! Breugem 2012 (JCP)
!   use param, only: al
!   use mls_param,only: dens_ratio, i_inv
!   implicit none

!   real,dimension(3)   :: vel_CM            ! velocity of the center of mass 
!   real,dimension(3)   :: pos_CM            ! position of the center of mass 
!   real,dimension(3)   :: omega_c           ! angular vel body frame
!   real,dimension(4)   :: quat              ! quaternions
!   real,dimension(3,3) :: AA                ! rotation matrix "body  = AA  space "
!   real,dimension(3)   :: fp,fp_m1
!   real,dimension(3)   :: tp,tp_m1
!   real                :: dt
!   real,dimension(3)   :: coll_vel
!   real,dimension(3)   :: coll_om_b
!   real,dimension(4)   :: coll_quat_dot
!   ! local variables 
!   real,dimension(3) :: vel_m1              ! velocity of the center of mass 
!   real,dimension(3) :: pos_m1              ! position of the center of mass 
!   real,dimension(3) :: omega_c_m1          ! angular vel body frame
!   real,dimension(4) :: quat_m1             ! quaternions previous time step
!   real,dimension(3) :: tp_b                ! collision torque in body frame
!   real,dimension(3) :: coll_vel_m1
!   real,dimension(4) :: coll_quat_dot_m1
 

!   ! ---------------------------------------------------- 
!   !              Translation
!   ! ---------------------------------------------------- 

!   coll_vel = 0.5*dt*al*(fp + fp_m1) 

!   pos_cm = pos_m1  +  0.5*dt*al*(coll_vel + coll_vel_m1)


!   ! ---------------------------------------------------- 
!   !               Rotation
!   ! ---------------------------------------------------- 

!   tp_b = tp + tp_m1

!   coll_om_b = 0.5*dt*al*matmul(I_inv, tp_b)


!   call quatMul(quat, coll_om_b, coll_quat_dot)

!   quat = quat_m1 + 0.25*dt*al*( coll_quat_dot + coll_quat_dot_m1 )

!   ! ---------------------------------------------------- 
!   !               Adding back to main vars
!   ! ---------------------------------------------------- 
!   vel_cm   = vel_m1      + coll_vel
!   omega_c  = omega_c_m1  + coll_om_b 

! end subroutine newton_coll






subroutine quatMul(quat,omega,result)
implicit none
real,dimension(4,4) :: QMAT            ! Allen Tildedley 3.37 ({qdot} = 0.5[Q]{wb})
real,dimension(4) :: result,quat,Wquat
real,dimension(3) :: omega 

quat = quat / norm2(quat)

QMAT(1,1) = quat(1)
QMAT(1,2) =-quat(2)
QMAT(1,3) =-quat(3)
QMAT(1,4) =-quat(4)

QMAT(2,1) = quat(2)
QMAT(2,2) = quat(1)
QMAT(2,3) =-quat(4)
QMAT(2,4) = quat(3)

QMAT(3,1) = quat(3)
QMAT(3,2) = quat(4)
QMAT(3,3) = quat(1)
QMAT(3,4) =-quat(2)

QMAT(4,1) = quat(4)
QMAT(4,2) =-quat(3)
QMAT(4,3) = quat(2)
QMAT(4,4) = quat(1)

Wquat = 0.d0
Wquat(2:4) = omega(1:3)
result = matmul(QMAT,Wquat)


end subroutine

subroutine quatMul2(quat,omega,result)
  implicit none
  real,dimension(4,4) :: QMAT            
  real,dimension(4) :: result,quat,Wquat
  real,dimension(3) :: omega 

  ! Allen Tildedley 3.37
  
  ! Space-fixed version!!!!
!      _
!    d q(t)     1   _     _
!  --------- = ---  W(t)  q(t)  
!     d t       2
!
  quat = quat / norm2(quat)

  QMAT(1,1) = quat(1)
  QMAT(1,2) =-quat(2)
  QMAT(1,3) =-quat(3)
  QMAT(1,4) =-quat(4)
  
  QMAT(2,1) = quat(2)
  QMAT(2,2) = quat(1)
  QMAT(2,3) = quat(4) ! Changed
  QMAT(2,4) = -quat(3) ! Changed
  
  QMAT(3,1) = quat(3)
  QMAT(3,2) = -quat(4) !Changed
  QMAT(3,3) = quat(1)
  QMAT(3,4) = quat(2) ! Changed
  
  QMAT(4,1) = quat(4) 
  QMAT(4,2) = quat(3) ! Changed
  QMAT(4,3) = -quat(2) ! Changed
  QMAT(4,4) = quat(1)
  
  Wquat = 0.d0
  Wquat(2:4) = omega(1:3)
  result = matmul(QMAT,Wquat)
  
  
  end subroutine


subroutine calc_rot_matrix(quat,AA)
! Given quaternions, calculate rotation mastrix
  implicit none

  real,dimension(4) :: quat        ! quaternions
  real,dimension(3,3) :: AA        ! rotation matrix "body  = AA  space "

! local variables
  real :: norm

! normalize quaternions
  norm = sqrt( quat(1)**2 + quat(2)**2 + quat(3)**2 + quat(4)**2 )
  quat = quat/norm

! rotation matrix (Allen Tildesley 3.36, but q_0 = quat(1)
  AA(1,1) = quat(1)**2 + quat(2)**2 - quat(3)**2 - quat(4)**2
  AA(1,2) = 2.0D0* ( quat(2)*quat(3) + quat(1)*quat(4) )
  AA(1,3) = 2.0D0* ( quat(2)*quat(4) - quat(1)*quat(3) )

  AA(2,1) = 2.0D0* ( quat(2)*quat(3) - quat(1)*quat(4) )
  AA(2,2) = quat(1)**2 - quat(2)**2 + quat(3)**2 - quat(4)**2
  AA(2,3) = 2.0D0* ( quat(3)*quat(4) + quat(1)*quat(2) )

  AA(3,1) = 2.0D0* ( quat(2)*quat(4) + quat(1)*quat(3) )
  AA(3,2) = 2.0D0* ( quat(3)*quat(4) - quat(1)*quat(2) )
  AA(3,3) = quat(1)**2 - quat(2)**2 - quat(3)**2 + quat(4)**2

end subroutine calc_rot_matrix

subroutine calc_quat_from_rotMatrix(quat,AA)
  ! Given rotation/orientation matrix, calculate quaternion
    implicit none
  
    real,dimension(4) :: quat        ! quaternions
    real,dimension(3,3) :: AA        ! rotation matrix "body  = AA  space "
  

    ! §6.5 of "Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors (Dibel, 2005)
    
    if ( ( AA(2,2) .gt. -AA(3,3) ) .and. ( AA(1,1) .gt. -AA(2,2) ) .and. ( AA(1,1) .gt. -AA(3,3) ) ) then
      !q0
      quat(1) = sqrt( 1.0 + AA(1,1) + AA(2,2) + AA(3,3) )
      quat(2) = ( AA(2,3) - AA(3,2) ) / quat(1)
      quat(3) = ( AA(3,1) - AA(1,3) ) / quat(1)
      quat(4) = ( AA(1,2) - AA(2,1) ) / quat(1)
      quat(1:4) = quat(1:4) * 0.5
    elseif ( ( AA(2,2) .lt. -AA(3,3) ) .and. ( AA(1,1) .gt. AA(2,2) ) .and. ( AA(1,1) .gt. AA(3,3) ) ) then
      !q1 
      quat(2) = sqrt( 1.0 + AA(1,1) - AA(2,2) - AA(3,3) ) 
      quat(1) = ( AA(2,3) - AA(3,2) ) / quat(2)
      quat(3) = ( AA(1,2) + AA(2,1) ) / quat(2)
      quat(4) = ( AA(3,1) + AA(1,3) ) / quat(2)
      quat(1:4) = quat(1:4) * 0.5
    elseif ( ( AA(2,2) .gt. AA(3,3) ) .and. ( AA(1,1) .lt. AA(2,2) ) .and. ( AA(1,1) .lt. -AA(3,3) ) ) then
      !q2
      quat(3) = sqrt( 1.0 - AA(1,1) + AA(2,2) - AA(3,3) )
      quat(1) = ( AA(3,1) - AA(1,3) ) / quat(3)
      quat(2) = ( AA(1,2) - AA(2,1) ) / quat(3)
      quat(4) = ( AA(2,3) - AA(3,2) ) / quat(3)
      quat(1:4) = quat(1:4) * 0.5
    elseif ( ( AA(2,2) .lt. AA(3,3) ) .and. ( AA(1,1) .lt. -AA(2,2) ) .and. ( AA(1,1) .lt. AA(3,3) ) ) then
      !q3
      quat(4) = sqrt( 1.0 - AA(1,1) - AA(2,2) + AA(3,3) )
      quat(1) = ( AA(1,2) - AA(2,1) ) / quat(4)
      quat(2) = ( AA(3,1) - AA(1,3) ) / quat(4)
      quat(3) = ( AA(2,3) - AA(3,2) ) / quat(4)
      quat(1:4) = quat(1:4) * 0.5
    else
      write(*,*) "No quat case, problem encountered?"
    endif

  ! normalize quaternions
    quat = quat / norm2 (quat)
  
  end subroutine calc_quat_from_rotMatrix


subroutine calc_rigidBody_params (COM,Vol,I_ij,nv,nf,xyz,vert_of_face,isGhostFace)

  ! Numerically compute "rigid"-body properties required for advancing Newton--Euler eqns.
  ! Moment of inertia tensor, centre of mass, volume, etc.
  ! The quantities are volume integrals, computed as surface integrals (sum over triangles)
  ! by way of the divergence theorem.
  ! Following:
  !   Computing the Moment of Inertia of a Solid Defined by a Triangle Mesh (Kallay 2006)
  !   https://github.com/erich666/jgt-code/blob/master/Volume_11/Number_2/Kallay2006/Moment_of_Inertia.cpp
  ! This is a Fortran port of the above C++ code
  ! For underlying theory, see also D. H. Eberly (2015, pp. 74) and Kallay (2006)

  use mls_param, only: dens_ratio

  implicit none

  integer :: i, nf, nv
  real,dimension(3,3) :: I_ij, Iij_buffer       ! Inertia tensor
  logical, dimension(nf) :: isGhostFace
  integer, dimension(3,nf) :: vert_of_face
  real, dimension(3,nv) :: xyz

  real :: Cx, Cy, Cz
  real :: Ixx, Iyy, Izz, Ixy, Ixz, Iyz
  real :: Vol
  real, dimension(3) :: I_princ, COM
  real :: x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4, v, r
  real :: xx, yy, zz

  ! For LAPACK eigenvalue routines
  integer :: LWORK = 3
  real, dimension(3) :: WORK
  integer :: IPIV(3), INFO


  ! COM coordinates in world coordinates
  Cx = 0.0d0
  Cy = 0.0d0
  Cz = 0.0d0
  
  !Inertia tensor components
  Ixx = 0.0d0
  Iyy = 0.0d0
  Izz = 0.0d0
  
  Ixy = 0.0d0
  Ixz = 0.0d0
  Iyz = 0.0d0

  Vol = 0.0d0
  
  do i = 1,nf
    if (.not. isGhostFace(i) ) then
       x1 = xyz(1, vert_of_face(1,i) ) 
       y1 = xyz(2, vert_of_face(1,i) ) 
       z1 = xyz(3, vert_of_face(1,i) ) 
  
       x2 = xyz(1, vert_of_face(2,i) ) 
       y2 = xyz(2, vert_of_face(2,i) ) 
       z2 = xyz(3, vert_of_face(2,i) )  
 
       x3 = xyz(1, vert_of_face(3,i) ) 
       y3 = xyz(2, vert_of_face(3,i) ) 
       z3 = xyz(3, vert_of_face(3,i) ) 
  
  
      ! Signed volume of this tetrahedron.
       v = x1*y2*z3 + y1*z2*x3 + x2*y3*z1 - (x3*y2*z1 + x2*y1*z3 + y3*z2*x1) 
  
       ! Accumulate volume
       Vol = Vol + v
  
      ! Contribution to the centroid
       x4 = x1 + x2 + x3 ;    Cx = Cx + (v * x4)
       y4 = y1 + y2 + y3 ;    Cy = Cy + (v * y4)
       z4 = z1 + z2 + z3 ;    Cz = Cz + (v * z4)
  
  
       ! Contribution to moment of inertia monomials
       ! Relative to world axis
       Ixx = Ixx + v * (x1*x1 + x2*x2 + x3*x3 + x4*x4)
       Iyy = Iyy + v * (y1*y1 + y2*y2 + y3*y3 + y4*y4)
       Izz = Izz + v * (z1*z1 + z2*z2 + z3*z3 + z4*z4)
       Ixy = Ixy + v * (y1*x1 + y2*x2 + y3*x3 + y4*x4)
       Ixz = Ixz + v * (z1*x1 + z2*x2 + z3*x3 + z4*x4)
       Iyz = Iyz + v * (z1*y1 + z2*y2 + z3*y3 + z4*y4)
    endif
  enddo
  
  
  !  "Post-process" step
  !  Centroid.
  !  The case _m = 0 needs to be addressed here.
   r = 1.0 / (4.0 * Vol)
   Cx = Cx * r
   Cy = Cy * r
   Cz = Cz * r
  
  !  Mass/volume: tetrahedral 1/6 correction
   Vol = Vol / 6.0
  
  !  Moment of inertia about the centroid
  !  Application of parallel axis theorem
  r = 1.0 / 120.0
  Ixy = Ixy * r - Vol * Cy*Cx
  Ixz = Ixz * r - Vol * Cz*Cx
  Iyz = Iyz * r - Vol * Cz*Cy
  
  xx = Ixx * r - Vol * Cx*Cx
  yy = Iyy * r - Vol * Cy*Cy
  zz = Izz * r - Vol * Cz*Cz
  
  Ixx = yy + zz
  Iyy = zz + xx
  Izz = xx + yy
  
  ! ------- Store -------------------------
  COM(1) = Cx ; COM(2) = Cy; COM(3) = Cz
  ! I = [ Ixx, -Ixy, -Ixz ;
  !      -Ixy,  Iyy, -Iyz ;
  !      -Ixz, -Iyz,  Izz ];
  I_ij(1,1) =  Ixx ;   I_ij(1,2) = -Ixy ;    I_ij(1,3) = -Ixz
  I_ij(2,1) = -Ixy ;   I_ij(2,2) =  Iyy ;    I_ij(2,3) = -Iyz
  I_ij(3,1) = -Ixz ;   I_ij(3,2) = -Iyz ;    I_ij(3,3) =  Izz

  ! Density prefactor: I_ij now stores the inertia matrix normalised on the (unit) fluid density
  I_ij = I_ij * dens_ratio

  ! Solve for principal values: 3x3 symmetric eigenproblem
  ! Eigenvalues I_princ(1:3), in ascending order by default
  ! Eigenvectors (principal axes) overwritten in I_ij(3,3)
  !call dsyev('V', 'U', 3, I_ij, 3, I_princ, WORK, LWORK, INFO)

  !Iij_buffer = I_ij

  ! Invert 3x3 symmetric I_ij matrix
  !call dgetrf(3, 3, Iij_buffer, 3, IPIV, INFO) ! Compute the LU factorization of I_ij

  !if (info /= 0) then
  !  print *, "Error in dgetrf"
  !  stop
  !end if

  !call dgetri(3, Iij_buffer, 3, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization

  !if (info /= 0) then
  !  print *, "Error in dgetri"
  !  stop
  !end if

  ! Assign to global variables
  !I_inv(:,:) = Iij_buffer(:,:)

  !write(*,*) "invI(1,:)", I_inv(1,:)
  !write(*,*) "invI(2,:)", I_inv(2,:)
  !write(*,*) "invI(3,:)", I_inv(3,:)

  !I_inv(1,1) = 1. / I_princ(1) 
  !I_inv(2,2) = 1. / I_princ(2) 
  !I_inv(3,3) = 1. / I_princ(3)

  !I_inv2 = 0.
  !I_inv2(1,1) = ( I_princ(2)-I_princ(3) ) / I_princ(1)
  !I_inv2(2,2) = ( I_princ(3)-I_princ(1) ) / I_princ(2)
  !I_inv2(3,3) = ( I_princ(1)-I_princ(2) ) / I_princ(3)

end subroutine calc_rigidBody_params

subroutine update_xyz
  use param
  use mls_param
  use mpih
  implicit none
  real,dimension(3)  :: om_dCM, pos, vel
  integer :: i,inp
  real,dimension(3,3):: AA,AAT,AAR
  real, dimension(2,2) :: Rot
  real :: radius,angle, om,tp
  real ::zmin,zmax


  do inp = 1,Nparticle
    call calc_rot_matrix(quat(:,inp),AA)
    AAT = transpose(AA) ! Space-aligned = A' * Body-fixed
    ! KZ: Should test whether this should be A or AAT

    do i = 1,maxnv
      if ( .not. isGhostVert(i,inp) ) then
        xyzv(:,i,inp) =  matmul( AAT, dxyzv_s(:,i,inp) ) + pos_CM(:,inp)
        !xyzv(:,i,inp) =  matmul( AA, dxyzv_s(:,i,inp) ) + pos_CM(:,inp)

      endif
    enddo

    ! COM-relative assignment
    !xyzv(2,:,inp) = xyzv(2,:,inp) + pos_CM(2,inp)
    !xyzv(3,:,inp) = xyzv(3,:,inp) + pos_CM(3,inp)

    ! New triangle centroid locations
    call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) 

    !-- velocity

    do i = 1,maxnf
      if (.not. isGhostFace(i,inp) ) then
        
        ! Add Urot = omega x r contribution to local surface velocity
        call cross(om_dCM(:), omega_c(:,inp), dxyz_s(:,i,inp)  )

        vel_tri(:,i,inp) = vel_CM(:,inp) + om_dCM(:)

      endif

    enddo

 end do

end subroutine update_xyz

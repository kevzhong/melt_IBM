subroutine convex_hull_q1(AA,inp)
  use mls_param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vx,vy,vz
  use param
  implicit none
  real, dimension(3,3) :: AA
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: volp
  real :: alpha
  real :: alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real tot_vol_1, r_x_u_1(3)

  ! get bounding box
  do i = 1,3
    lim(i,1) = minval(tri_bar(i,:,inp))
    lim(i,2) = maxval(tri_bar(i,:,inp))
  end do

  ind = floor(lim*dx1) + 1 ! compute indices cell centered
  ! expanding bounding box to be extra safe
  ind(:,1) = ind(:,1) - 7; ind(:,2) = ind(:,2) + 7;

  tot_vol_1    = 0.
  r_x_u_1(1:3) = 0.


  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)


           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices(k,x_gc)


           if (k.ge.kstart.and.k.le.kend) then 
             ! position cell centered location
             x_grid(1) = xc(i)
             x_grid(2) = ym(j)
             x_grid(3) = zm(k)
             r = x_grid - x_GC ! relative distance 

             call level_set1(i,j,k,AA,inp,x_gc,alpha) 

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1
             alpha_q = alpha*celvol*vx(ii,jj,k)
             if(VOFx(ii,jj,k).lt.1.0)then
             VOFx(ii,jj,k) = VOFx(ii,jj,k)
             else
             VOFx(ii,jj,k) = 1.0 - alpha
             end if
             ! compute int u over V
             u_tot(1,inp) = u_tot(1,inp) + alpha_q

             ! compute int r x u over V, r = xgrid - pos_cm
             r_x_u_1(2) = r_x_u_1(2) + alpha_q * r(3)
             r_x_u_1(3) = r_x_u_1(3) - alpha_q * r(2)

             tot_vol_1 = tot_vol_1 + celvol*alpha
             endif

      end do
    end do
  end do

  call mpi_globalsum_double_var(tot_vol_1)
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_globalsum_double_var(u_tot(1,inp))
  call mpi_globalsum_double_arr(r_x_u_1,3)

  u_tot(1,inp)       = u_tot(1,inp) / tot_vol_1
  r_x_u_tot(1:3,inp) = r_x_u_tot(1:3,inp) + r_x_u_1 / tot_vol_1

 !volp = (4./3.)*pi*0.5**3
 !if (myid.eq.0) print*, "volume fraction", tot_vol_1/ volp

! if(myid.eq.0)then
!   open(112,file='flowmov/vol.txt',status='unknown', position='append')
!         write(112,'(40E17.5)') time, tot_vol /volp
!   close(112)
! end if
end subroutine

subroutine level_set1(ic,jc,kc,AA,inp,x_gc,alpha)
  use param, only: xm,yc,zc
  use geom
  implicit none
  integer :: i,j,k    ! int running over corners
  integer :: ic,jc,kc ! center position 
  real :: phi,alpha,phi_tot,inva,invb,invc
  real, dimension(3)   :: x_GC
  real, dimension(3)   :: r, v, x_grid
  real, dimension(3,3) :: AA
  integer :: inp

  ! Compute alpha over cell, see Kempe 2012 (JCP)
  ! alpha = sum(-phi(m) * H(-phi(m))) / sum(|phi(m)|)
  ! alpha in [0,1]
  
  phi_tot = 0.
  alpha = 0.

  ! run over 8 corners of cell
  do i = ic-1,ic
    do j = jc,jc+1
      do k = kc,kc+1

        x_grid(1) = xm(i)
        x_grid(2) = yc(j)
        x_grid(3) = zc(k)

        phi = loopoverbeams(x_grid,x_gc,AA,inp) 

        if (phi.gt.0.) alpha = alpha + phi
        phi_tot = phi_tot + abs(phi) 

      end do
    end do
  end do

alpha = alpha / phi_tot
end subroutine


!=================================
!    q2
!=================================

subroutine convex_hull_q2(AA,inp)
  use mls_param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vy
  use geom
  implicit none
  real, dimension(3,3) :: AA
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: volp
  real :: alpha, alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real tot_vol_2, r_x_u_2(3)

  ! get bounding box
  do i = 1,3
    lim(i,1) = minval(tri_bar(i,:,inp))
    lim(i,2) = maxval(tri_bar(i,:,inp))
  end do

  ind = floor(lim*dx1) + 1 ! compute indices cell centered
  ! expanding bounding box to be extra safe
  ind(:,1) = ind(:,1) - 7; ind(:,2) = ind(:,2) + 7;

  tot_vol_2 = 0.
  r_x_u_2   = 0.


  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices(k,x_gc)


           if (k.ge.kstart.and.k.le.kend) then 
             ! position cell centered location
             x_grid(1) = xm(i)
             x_grid(2) = yc(j)
             x_grid(3) = zm(k)
             r = x_grid - x_GC ! relative distance 

             call level_set2(i,j,k,AA,x_GC,inp,alpha) 

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1
             alpha_q = alpha*celvol*vy(ii,jj,k)
             

             if(VOFy(ii,jj,k).lt.1.0)then
             VOFy(ii,jj,k) = VOFy(ii,jj,k)
             else
             VOFy(ii,jj,k) = 1.0 - alpha
             end if

             ! compute int u over V 
             u_tot(2,inp) = u_tot(2,inp) + alpha_q 

             ! compute int r x u over V, r = xgrid - pos_cm
             r_x_u_2(1) = r_x_u_2(1) - alpha_q * r(3)
             r_x_u_2(3) = r_x_u_2(3) + alpha_q * r(1)

             tot_vol_2 = tot_vol_2 + celvol *alpha
             endif

      end do
    end do
  end do

  call mpi_globalsum_double_var(tot_vol_2)
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_globalsum_double_var(u_tot(2,inp))
  call mpi_globalsum_double_arr(r_x_u_2,3)

  u_tot(2,inp)       =  u_tot(2,inp) / tot_vol_2
  r_x_u_tot(1:3,inp) = r_x_u_tot(1:3,inp) + r_x_u_2  / tot_vol_2


! volp = (4./3.)*pi*0.5**3
! if (myid.eq.0) print*, "volume fraction", tot_vol_2/ volp
! if(myid.eq.0)then
!   open(112,file='flowmov/vol.txt',status='unknown', position='append')
!         write(112,'(40E17.5)') time, tot_vol /volp
!   close(112)
! end if
end subroutine

subroutine level_set2(ic,jc,kc,AA,x_GC,inp,alpha)
  use param, only: xc,ym,zc
  use geom
  implicit none
  integer :: i,j,k    ! int running over corners
  integer :: ic,jc,kc ! center position 
  real :: phi,alpha,phi_tot,inva,invb,invc
  real, dimension(3)   :: x_GC
  real, dimension(3)   :: r, v, x_grid
  real, dimension(3,3) :: AA
  integer :: inp
  ! Compute alpha over cell, see Kempe 2012 (JCP)
  ! alpha = sum(-phi(m) * H(-phi(m))) / sum(|phi(m)|)
  ! alpha in [0,1]
  
  phi_tot = 0.
  alpha  =  0.
  ! run over 8 corners of cell
  do i = ic,ic+1
    do j = jc-1,jc
      do k = kc,kc+1

        ! set location of cell corner
        x_grid(1) = xc(i)
        x_grid(2) = ym(j)
        x_grid(3) = zc(k)


        phi = loopoverbeams(x_grid,x_gc,AA,inp) 

        if (phi.gt.0.) alpha = alpha + phi
        phi_tot = phi_tot + abs(phi) 

      end do
    end do
  end do

  alpha = alpha / phi_tot

end subroutine

!=================================
!    q3
!=================================

subroutine convex_hull_q3(AA,inp)
  use mls_param
  use param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vz
  implicit none
  real, dimension(3,3) :: AA
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: volp
  real :: alpha, alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real tot_vol_3, r_x_u_3(3)

  ! get bounding box
  do i = 1,3
    lim(i,1) = minval(tri_bar(i,:,inp))
    lim(i,2) = maxval(tri_bar(i,:,inp))
  end do

  ind = floor(lim*dx1) + 1 ! compute indices cell centered
  ! expanding bounding box to be extra safe
  ind(:,1) = ind(:,1) - 7; ind(:,2) = ind(:,2) + 7;

  tot_vol_3 = 0.
  r_x_u_3   = 0.


  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices(k,x_gc)

           if (k.ge.kstart.and.k.le.kend) then 
             ! position cell centered location
             x_grid(1) = xm(i)
             x_grid(2) = ym(j)
             x_grid(3) = zc(k)
             r = x_grid - x_GC ! relative distance 

             call level_set3(i,j,k,AA,x_GC,inp,alpha) 
             ! compute int u over V 

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1
             alpha_q = alpha*celvol*vz(ii,jj,k)

             if(VOFz(ii,jj,k).lt.1.0)then
             VOFz(ii,jj,k) = VOFz(ii,jj,k)
             else
             VOFz(ii,jj,k) = 1.0 - alpha
             end if

             ! compute int u over V 
             u_tot(3,inp) = u_tot(3,inp) + alpha_q 

             ! compute int r x u over V, r = xgrid - pos_cm
             r_x_u_3(1) = r_x_u_3(1) + alpha_q * r(2)
             r_x_u_3(2) = r_x_u_3(2) - alpha_q * r(1)

             tot_vol_3 = tot_vol_3 + celvol*alpha
             endif

      end do
    end do
  end do

  call mpi_globalsum_double_var(tot_vol_3)
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_globalsum_double_var(u_tot(3,inp))
  call mpi_globalsum_double_arr(r_x_u_3,3)

  u_tot(3,inp)       = u_tot(3,inp) / tot_vol_3
  r_x_u_tot(1:3,inp) = r_x_u_tot(1:3,inp) + r_x_u_3  / tot_vol_3

!volp = (4./3.)*pi*0.5**3
!if (myid.eq.0) print*, "volume fraction", tot_vol_3/ volp

! if(myid.eq.0)then
!   open(112,file='flowmov/vol.txt',status='unknown', position='append')
!         write(112,'(40E17.5)') time, tot_vol /volp
!   close(112)
! end if
end subroutine

subroutine level_set3(ic,jc,kc,AA,x_GC,inp,alpha)
  use param, only: xc,yc,zm
  use geom
  implicit none
  integer :: i,j,k    ! int running over corners
  integer :: ic,jc,kc ! center position 
  real :: phi,alpha,phi_tot,inva,invb,invc
  real, dimension(3)   :: x_GC
  real, dimension(3)   :: r, v, x_grid
  real, dimension(3,3) :: AA
  integer :: inp

  ! Compute alpha over cell, see Kempe 2012 (JCP)
  ! alpha = sum(-phi(m) * H(-phi(m))) / sum(|phi(m)|)
  ! alpha in [0,1]
  
  phi_tot = 0.
  alpha  =  0.
  ! run over 8 corners of cell
  do i = ic,ic+1
    do j = jc,jc+1
      do k = kc-1,kc

        ! set location of cell corner
        x_grid(1) = xc(i)
        x_grid(2) = yc(j)
        x_grid(3) = zm(k)

        phi = loopoverbeams(x_grid,x_gc,AA,inp) 

        if (phi.gt.0.) alpha = alpha + phi
        phi_tot = phi_tot + abs(phi) 



      end do
    end do
  end do

alpha = alpha / phi_tot
end subroutine

subroutine get_periodic_indices(k,x)
  use param
  implicit none
  integer :: k
  real    :: x(3)

  if (k .ge. n3) then
     k = k - n3m
    x(3) = x(3) - zlen
  end if

  if (k .lt. 1) then
     k = k + n3m
     x(3) = x(3) + zlen
  end if
end subroutine

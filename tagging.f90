subroutine convex_hull_q12(AA,inp)
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

  ! run over cell centered grid and find int(u)dV and int(r x u)dV
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)


           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices(k,x_gc)


           if (k.ge.kstart.and.k.le.kend) then 
             call level_set12(i,j,k,AA,inp,x_gc,alpha)

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1
             if(ax(ii,jj,k).lt.1.0)then
             ax(ii,jj,k) = ax(ii,jj,k)
             else
             ax(ii,jj,k) = 1.0 - alpha
             end if
             endif

      end do
    end do
  end do

end subroutine

subroutine level_set12(ic,jc,kc,AA,inp,x_gc,alpha)
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

subroutine convex_hull_q22(AA,inp)
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

  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices2(k,x_gc)


           if (k.ge.kstart.and.k.le.kend) then 
           call level_set22(i,j,k,AA,x_GC,inp,alpha)

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1
            
             if(ay(ii,jj,k).lt.1.0)then
             ay(ii,jj,k) = ay(ii,jj,k)
             else
             ay(ii,jj,k) = 1.0 - alpha
             end if

           endif

      end do
    end do
  end do

end subroutine

subroutine level_set22(ic,jc,kc,AA,x_GC,inp,alpha)
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

subroutine convex_hull_q32(AA,inp)
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

  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices2(k,x_gc)

           if (k.ge.kstart.and.k.le.kend) then 

             call level_set32(i,j,k,AA,x_GC,inp,alpha)
             ! compute int u over V 

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1

             if(az(ii,jj,k).lt.1.0)then
             az(ii,jj,k) = az(ii,jj,k)
             else
             az(ii,jj,k) = 1.0 - alpha
             end if
             endif

      end do
    end do
  end do

end subroutine

subroutine level_set32(ic,jc,kc,AA,x_GC,inp,alpha)
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

subroutine get_periodic_indices2(k,x)
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

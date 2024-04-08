  ! Calculate the volume-of-fluid / solid field given the bounding-box indices defined by ind(3,2) for particle inp
  ! Using level-set scheme of Kempe & Frohlich (2012) eqns 32-33
  ! The signed distance function, phi, is evaluated as the distance from the Eulerian-cell-corner point to the plane of the closest triangle centroid

  subroutine convex_hull_q12(ind,inp)
  use mls_param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vx,vy,vz
  use param
  implicit none
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: volp
  real :: alpha
  real :: alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real tot_vol_1, r_x_u_1(3)



  ! run over cell centered grid and find int(u)dV and int(r x u)dV
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices(k,x_gc)


           if (k.ge.kstart.and.k.le.kend) then 
             call level_set12(i,j,k,inp,x_gc,alpha)

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1
             if(VOFx(ii,jj,k).lt.1.0)then
             VOFx(ii,jj,k) = VOFx(ii,jj,k)
             else
            VOFx(ii,jj,k) = 1.0 - alpha
             end if
             endif

      end do
    end do
  end do

end subroutine

subroutine level_set12(ic,jc,kc,inp,x_gc,alpha)
  use param, only: xm,yc,zc
  use geom
  implicit none
  integer :: i,j,k    ! int running over corners
  integer :: ic,jc,kc ! center position 
  integer :: nf
  real :: phi,alpha,phi_tot,inva,invb,invc
  real, dimension(3)   :: x_GC
  real, dimension(3)   :: r, v, x_grid
  !real, dimension(3,nf) :: tri_bar, tri_nor
  integer :: inp, tri_ind

  ! Compute alpha over cell, see Kempe 2012 (JCP)
  ! alpha = sum(-phi(m) * H(-phi(m))) / sum(|phi(m)|)
  ! alpha in [0,1]
  
  phi_tot = 0.
  alpha = 0.

  !x_grid(1) = xm(ic)
  !x_grid(2) = ym(jc)
  !x_grid(3) = zm(kc)
  !call find_closestTri_ind(tri_ind,x_grid,tri_bar,nf)

  ! run over 8 corners of cell
  do i = ic-1,ic
    do j = jc,jc+1
      do k = kc,kc+1

        x_grid(1) = xm(i)
        x_grid(2) = yc(j)
        x_grid(3) = zc(k)

      ! KZ: If this proves to be expensive, can just fix tri_ind = const. for entire cell
        call find_closestTri_ind(tri_ind,x_grid,tri_bar(:,:,inp),isGhostFace(:,inp),maxnf)

        !phi = loopoverbeams(x_grid,x_gc,AA,inp) 
        phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))

       !if (phi.gt.0.) alpha = alpha + phi
        if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)

        phi_tot = phi_tot + abs(phi) 

      end do
    end do
  end do

alpha = alpha / phi_tot
end subroutine


!=================================
!    q2
!=================================

subroutine convex_hull_q22(ind,inp)
  use mls_param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vy
  use geom
  implicit none
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: volp
  real :: alpha, alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real tot_vol_2, r_x_u_2(3)

  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices2(k,x_gc)


           if (k.ge.kstart.and.k.le.kend) then 
           call level_set22(i,j,k,x_GC,inp,alpha)

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1
            
             if(VOFy(ii,jj,k).lt.1.0)then
             VOFy(ii,jj,k) = VOFy(ii,jj,k)
             else
             VOFy(ii,jj,k) = 1.0 - alpha
             end if

           endif

      end do
    end do
  end do

end subroutine

subroutine level_set22(ic,jc,kc,x_GC,inp,alpha)
  use param, only: xc,ym,zc
  use geom
  implicit none
  integer :: i,j,k    ! int running over corners
  integer :: ic,jc,kc ! center position 
  real :: phi,alpha,phi_tot,inva,invb,invc
  real, dimension(3)   :: x_GC
  real, dimension(3)   :: r, v, x_grid
  integer :: inp, tri_ind
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


        ! KZ: If this proves to be expensive, can just fix tri_ind = const. for entire cell
        call find_closestTri_ind(tri_ind,x_grid,tri_bar(:,:,inp),isGhostFace(:,inp),maxnf)

        phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))

        !phi = loopoverbeams(x_grid,x_gc,AA,inp) 

        !if (phi.gt.0.) alpha = alpha + phi
        if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)

        phi_tot = phi_tot + abs(phi) 

      end do
    end do
  end do

  alpha = alpha / phi_tot

end subroutine

!=================================
!    q3
!=================================

subroutine convex_hull_q32(ind,inp)
  use mls_param
  use param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vz
  implicit none
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: volp
  real :: alpha, alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real tot_vol_3, r_x_u_3(3)

  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices2(k,x_gc)

           if (k.ge.kstart.and.k.le.kend) then 

             call level_set32(i,j,k,x_GC,inp,alpha)
             ! compute int u over V 

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1

             if(VOFz(ii,jj,k).lt.1.0)then
             VOFz(ii,jj,k) = VOFz(ii,jj,k)
             else
             VOFz(ii,jj,k) = 1.0 - alpha
             end if
             endif

      end do
    end do
  end do

end subroutine


subroutine level_set32(ic,jc,kc,x_GC,inp,alpha)
  use param, only: xc,yc,zm
  use geom
  implicit none
  integer :: i,j,k    ! int running over corners
  integer :: ic,jc,kc ! center position 
  real :: phi,alpha,phi_tot,inva,invb,invc
  real, dimension(3)   :: x_GC
  real, dimension(3)   :: r, v, x_grid
  integer :: inp, tri_ind

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

        ! KZ: If this proves to be expensive, can just fix tri_ind = const. for entire cell
        call find_closestTri_ind(tri_ind,x_grid,tri_bar(:,:,inp),isGhostFace(:,inp),maxnf)

        phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))

        !phi = loopoverbeams(x_grid,x_gc,AA,inp) 

        !if (phi.gt.0.) alpha = alpha + phi
        if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)

        phi_tot = phi_tot + abs(phi) 


      end do
    end do
  end do

alpha = alpha / phi_tot
end subroutine

subroutine convex_hull_qc2(ind,inp)
  use mls_param
  use param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vz
  implicit none
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: volp
  real :: alpha, alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real tot_vol_3, r_x_u_3(3)

  ! run over cell centered grid and find int(u)dV and int(r x u)dV 
  do i = ind(1,1),ind(1,2)
    do j = ind(2,1),ind(2,2)
      do kk = ind(3,1),ind(3,2)

           ! periodic BC
           x_gc = pos_cm(1:3,inp) 
           k = kk
           call get_periodic_indices2(k,x_gc)

           if (k.ge.kstart.and.k.le.kend) then 

             call level_setc2(i,j,k,x_GC,inp,alpha)
             ! compute int u over V 

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1

             if(VOFp(ii,jj,k).lt.1.0)then
             VOFp(ii,jj,k) = VOFp(ii,jj,k)
             else
              VOFp(ii,jj,k) = 1.0 - alpha
             end if
             endif

      end do
    end do
  end do

end subroutine


subroutine level_setc2(ic,jc,kc,x_GC,inp,alpha)
  use param, only: xc,yc,zc
  use geom
  implicit none
  integer :: i,j,k    ! int running over corners
  integer :: ic,jc,kc ! center position 
  real :: phi,alpha,phi_tot,inva,invb,invc
  real, dimension(3)   :: x_GC
  real, dimension(3)   :: r, v, x_grid
  integer :: inp, tri_ind

  ! Compute alpha over cell, see Kempe 2012 (JCP)
  ! alpha = sum(-phi(m) * H(-phi(m))) / sum(|phi(m)|)
  ! alpha in [0,1]
  
  phi_tot = 0.
  alpha  =  0.
  ! run over 8 corners of cell
  do i = ic,ic+1
    do j = jc,jc+1
      do k = kc,kc+1

        ! set location of cell corner
        x_grid(1) = xc(i)
        x_grid(2) = yc(j)
        x_grid(3) = zc(k)

        ! KZ: If this proves to be expensive, can just fix tri_ind = const. for entire cell
        call find_closestTri_ind(tri_ind,x_grid,tri_bar(:,:,inp),isGhostFace(:,inp),maxnf)

        phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))

        !phi = loopoverbeams(x_grid,x_gc,AA,inp) 

        !if (phi.gt.0.) alpha = alpha + phi
        if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)

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

subroutine find_closestTri_ind(tri_ind,x_grid,tri_bar,isGhostFace,nf)

  ! Find the index of the triangle centroid which lies closest in space to some point in space, x_grid

  use param
  implicit none
  integer :: i, nf, tri_ind
  real :: mindist
  real    :: x_grid(3)
  real, dimension(3,nf) :: tri_bar
  logical, dimension(nf) :: isGhostFace
  real  :: dist

  mindist = 1.0e6
  tri_ind = 0

  do i = 1,nf
    if (isGhostFace(i) .eqv. .false.) then
      dist = norm2( x_grid - tri_bar(:,i) )

      if (dist .lt. mindist) then
        mindist = dist
        tri_ind = i
      endif
    endif
  enddo

  !if (tri_ind .eq. 0 ) then
  !  write(*,*) "Something went wrong finding the closest centroid"
  !  exit
  !endif

end subroutine

subroutine get_bbox_inds(bbox_inds,inp)
  ! Retrieve the xyz indices of the bounding box for a given particle
  use param
  use mls_param
  implicit none
  integer :: i, nf, tri_ind
  real, dimension(3,2) :: lim
  integer, dimension(3,2) :: bbox_inds
  integer  :: inp, padSize

  ! Padding size indices for safety
  padSize = 2

! get bounding box
  do i = 1,3
    lim(i,1) = minval( pack(xyzv(i,:,inp) , .not. isGhostVert(:,inp)  ) )
    lim(i,2) = maxval( pack(xyzv(i,:,inp) , .not. isGhostVert(:,inp)  ) )
  end do

  bbox_inds = floor(lim*dx1) + 1 ! compute indices cell centered

  ! expanding bounding box to be extra safe
  bbox_inds(:,1) = bbox_inds(:,1) - padSize
  bbox_inds(:,2) = bbox_inds(:,2) + padSize

end subroutine
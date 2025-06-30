  ! Calculate the volume-of-fluid / solid field given the bounding-box indices defined by ind(3,2) for particle inp
  ! Using level-set scheme of Kempe & Frohlich (2012) eqns 32-33
  ! The signed distance function, phi, is evaluated as the distance from the Eulerian-cell-corner point to the plane of the closest triangle centroid

  subroutine sphereTagging(ind,inp)
    use param, only: VOFx, VOFy, VOFz, VOFp, dx1, dx2,dx3
    use mls_param, only: celvol
    use mpih
    use mpi_param
    implicit none
    integer, dimension(3,2) :: ind
    integer :: inp
    real(8) :: vol_sphere


    VOFx(:,:,:) = 1.
    VOFy(:,:,:) = 1.
    VOFz(:,:,:) = 1.
    VOFp(:,:,:) = 1.

    call convex_hull_qc2(ind,inp)
    call convex_hull_q12(ind,inp)
    call convex_hull_q22(ind,inp)
    call convex_hull_q32(ind,inp)

    ! vol_sphere = sum( 1.0 - VOFp(:,:,kstart:kend) ) * celvol
    ! call MPI_ALLREDUCE(MPI_IN_PLACE,vol_sphere,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
    ! if (myid .eq. 0) then
    !  write(*,*) "Vsphere = ", vol_sphere
    ! endif

  end subroutine sphereTagging

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
              !x_grid(1) = xc(i)
              !x_grid(2) = ym(j)
              !x_grid(3) = zm(k)
              !r = x_grid - x_GC ! relative distance 

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
  
          phi = signDist_sphere(x_grid,x_gc,inp) 
  
          !phi = loopoverbeams(x_grid,x_gc,AA,inp) 
          !phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))
  
         !if (phi.gt.0.) alpha = alpha + phi
          !if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)

          alpha = alpha + (-phi) * heaviside(-phi)
  
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
             call get_periodic_indices(k,x_gc)
  
  
             if (k.ge.kstart.and.k.le.kend) then 
              !x_grid(1) = xm(i)
              !x_grid(2) = yc(j)
              !x_grid(3) = zm(k)
              !r = x_grid - x_GC ! relative distance 

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
          !call find_closestTri_ind(tri_ind,x_grid,tri_bar(:,:,inp),isGhostFace(:,inp),maxnf)
  
          !phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))
  
          !phi = loopoverbeams(x_grid,x_gc,AA,inp) 
          phi = signDist_sphere(x_grid,x_gc,inp) 

  
          !if (phi.gt.0.) alpha = alpha + phi
          !if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)

          alpha = alpha + (-phi) * heaviside(-phi)
  
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
             call get_periodic_indices(k,x_gc)
  
             if (k.ge.kstart.and.k.le.kend) then 
              !x_grid(1) = xm(i)
              !x_grid(2) = ym(j)
              !x_grid(3) = zc(k)
              !r = x_grid - x_GC ! relative distance 

  
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
          !call find_closestTri_ind(tri_ind,x_grid,tri_bar(:,:,inp),isGhostFace(:,inp),maxnf)
  
          !phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))
  
          !phi = loopoverbeams(x_grid,x_gc,AA,inp) 
  
          phi = signDist_sphere(x_grid,x_gc,inp) 


          !if (phi.gt.0.) alpha = alpha + phi
          !if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)
  
          alpha = alpha + (-phi) * heaviside(-phi)

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
    use local_arrays, only: temp
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: volp
    real :: alpha, alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real tot_vol_3, r_x_u_3(3)
    real :: u_imh, u_iph, v_jmh, v_jph, w_kmh, w_kph, h31, h32, h33
    real :: udx1, udx2, udx3
    integer :: jc,kc
    integer :: km,kp,jm,jp,ic,im,ip

  
    udx1=dx1*0.5d0
    udx2=dx2*0.5d0
    udx3=dx3*0.5d0

    ! run over cell centered grid and find int(u)dV and int(r x u)dV 
    do i = ind(1,1),ind(1,2)
      do j = ind(2,1),ind(2,2)
        do kk = ind(3,1),ind(3,2)
  
             ! periodic BC
             x_gc = pos_cm(1:3,inp) 
             k = kk
             call get_periodic_indices(k,x_gc)
  
             if (k.ge.kstart.and.k.le.kend) then 
  
               call level_setc2(i,j,k,x_GC,inp,alpha)
               ! compute int u over V 
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1

               kc = k
               km=kc-1
               kp=kc+1

               jc = jj
               jm=jmv(jj)
               jp=jpv(jj)

               ic = ii
               im=imv(ii)
               ip=ipv(ii)
  
               if(VOFp(ii,jj,k).lt.1.0)then
               VOFp(ii,jj,k) = VOFp(ii,jj,k)
               else
                VOFp(ii,jj,k) = 1.0 - alpha
               end if

              !  ! Solid cell d_Usolid_dxj
              !  if (VOFp(ii,jj,k).lt. 1.0e-6) then

              !   solid_mask(ii,jj,k) = .true.

              !   !                d  u T   |          1   [                              ]
              !   !             ----------- |  =     ----- |  uT |      -      uT |       |
              !   !                d   x    |i,j,k     dx  [     i+1/2            i-1/2   ]

              !   ! uT |_{i-1/2}
              !   x_grid(1) = xc(i)
              !   x_grid(2) = ym(j)
              !   x_grid(3) = zm(kk)
              !   r = x_grid - x_GC ! relative distance 
              !   u_imh = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

              !   ! uT |_{i+1/2}
              !   x_grid(1) = xc(i+1)
              !   r = x_grid - x_GC ! relative distance 
              !   u_iph = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

              !   h31=( u_iph*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
              !   -u_imh*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1

              !   !                d  v T   |          1   [                              ]
              !   !             ----------- |  =     ----- |  vT |      -      vT |       |
              !   !                d   y    |i,j,k     dy  [     j+1/2            j-1/2   ] 

              !   ! vT |_{j-1/2}
              !   x_grid(1) = xm(i)
              !   x_grid(2) = yc(j)
              !   x_grid(3) = zm(kk)
              !   r = x_grid - x_GC ! relative distance 
              !   v_jmh = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

              !   ! vT |_{j+1/2}
              !   x_grid(2) = yc(j+1)
              !   r = x_grid - x_GC ! relative distance 
              !   v_jph = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

              !   h32=( v_jph*(temp(ic,jp,kc)+temp(ic,jc,kc)) &
              !   -v_jmh*(temp(ic,jc,kc)+temp(ic,jm,kc)) )*udx2


              !   !                d  w T   |          1   [                              ]
              !   !             ----------- |  =     ----- |  wT |      -      wT |       |
              !   !                d   z    |i,j,k     dz  [     k+1/2            k-1/2   ]

              !   ! wT |_{k-1/2}
              !   x_grid(1) = xm(i)
              !   x_grid(2) = ym(j)
              !   x_grid(3) = zc(kk)
              !   r = x_grid - x_GC ! relative distance 
              !   w_kmh = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

              !   ! wT |_{k+1/2}
              !   x_grid(3) = zc(kk+1)
              !   r = x_grid - x_GC ! relative distance 
              !   w_kph = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

              !   h33=( w_kph*(temp(ic,jc,kp)+temp(ic,jc,kc)) &
              !   -w_kmh*(temp(ic,jc,kc)+temp(ic,jc,km)) )*udx3

              !   d_UsolidT_dxj(ii,jj,k) = (h31+h32+h33)

              !  endif




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
          !call find_closestTri_ind(tri_ind,x_grid,tri_bar(:,:,inp),isGhostFace(:,inp),maxnf)
  
          !phi = signed_distance(x_grid,tri_bar(1:3,tri_ind,inp),tri_nor(1:3,tri_ind,inp))
  
          !phi = loopoverbeams(x_grid,x_gc,AA,inp) 
  
          !if (phi.gt.0.) alpha = alpha + phi

          phi = signDist_sphere(x_grid,x_gc,inp) 


          !if (-phi .gt. 0.0d0 ) alpha = alpha + (-phi)
  
          alpha = alpha + (-phi) * heaviside(-phi)
          
          phi_tot = phi_tot + abs(phi) 
  
  
        end do
      end do
    end do
  
  alpha = alpha / phi_tot
  end subroutine
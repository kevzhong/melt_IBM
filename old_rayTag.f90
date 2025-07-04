 
 ! Tag fluid / solid cells in the Eulerian grid given the bounding-box indices defined by ind(3,2) for particle inp
 ! Bounding box indices are obtained from the get_bbox_inds(ind,inp) routine below
 ! Cells are identified as either fluid (VOF = 1), solid (VOF = 0.0), or interface (VOF between 0 and 1)
 ! The tagging operation (fluid, interface, or solid) is done through ray-tagging
 ! That is, by performing queries on ray-triangle intersections
 ! Ray-triangle intersections are efficiently computed using the Möller & Trumbore (2005) algorithm
 ! As noted in §7.5 of O'Rourke (1998), counting the no. of ray-triangle intersections will tell us if a cell is internal/external (solid/fluid)
 ! provided the geometry is 'watertight' (a closed, triangulated geometry)
 !      Cell is external (fluid): even number of intersections:  
 !      Cell is internal (solid): odd number of intersections:   
 !      Cell is interface: ray-triangle intersection point lies inside cell volume coordinates [x-dx/2, x+dx/2],[y-dy/2, y+dy/2],[z-dz/2, z+dz/2]
 ! Interface treatment for VOF somewhat ambiguous still
 !
 ! Besides Möller & Trumbore (2005) algorithm, some additional resources on ray-tracing are:
 !
 ! - IBM lecture notes (de Tullio et al., 2014) http://people.uniroma2.it/roberto.verzicco/IBnotes.pdf
 ! - "Computational Geometry in C" (O'Rourke, 1988, §7.5)
 ! - https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
 ! - "Essential Mathematics for Games and Interactive Applications" (Van Verth & Bishop, 2016, §12.3.5.2)
 ! - "Real-time Rendering" (Akenine-Möller et al., 2018, §22.8)

   subroutine tagCells(ind,inp)
    use param, only: VOFx, VOFy, VOFz, VOFp, solid_mask, dx1, dx2,dx3
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
    solid_mask(:,:,:) = .false.

    call convex_hull_qc2(ind,inp)
    call convex_hull_q12(ind,inp)
    call convex_hull_q22(ind,inp)
    call convex_hull_q32(ind,inp)

    !vol_sphere = sum( 1.0 - VOFp(:,:,kstart:kend) ) * celvol
    !call MPI_ALLREDUCE(MPI_IN_PLACE,vol_sphere,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
    !if (myid .eq. 0) then
    ! write(*,*) "Vsphere = ", vol_sphere
    !endif

  end subroutine tagCells

  subroutine convex_hull_q12(ind,inp)
    ! Tagging for x-staggered vx grid
  use mls_param
  use mpih
  use mpi_param, only: kstart,kend
  use local_arrays, only: vx,vy,vz
  use param
  implicit none
  real, dimension(3)   :: x_grid
  real, dimension(3)   :: r, x_GC
  real, dimension(3,2) :: lim
  real :: vof
  real :: alpha
  real :: alpha_q
  integer :: inp, i,j,k, ii,jj,kk
  integer, dimension(3,2) :: ind
  real ,dimension(3) :: Q, C

  ! Q stores the cell-centre coordinates and is the query point
  ! C is the control-point and stores the ray origin

  ! Random control-point well-outside computational domain
  C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]

  do i = ind(1,1),ind(1,2)
    Q(1) = xc(i)
    do j = ind(2,1),ind(2,2)
        Q(2) = ym(j)
      do kk = ind(3,1),ind(3,2)
           Q(3) = zm(kk)

           ! periodic BC
           k = modulo(kk-1,n3m) + 1

           if (k.ge.kstart.and.k.le.kend) then 

            !call rayTagQ(vof,C,Q,inp)
            call rayTagQ_faster(vof,Q,kk,inp)

             ii = modulo(i-1,n1m) + 1
             jj = modulo(j-1,n2m) + 1

             if(VOFx(ii,jj,k).lt.1.0)then
             VOFx(ii,jj,k) = VOFx(ii,jj,k)
             else
              VOFx(ii,jj,k) = vof
             end if

             endif

      end do
    end do
  end do

end subroutine

!=================================
!    q2
!=================================

subroutine convex_hull_q22(ind,inp)
      ! Tagging for y-staggered vy grid
    use mls_param
    use mpih
    use mpi_param, only: kstart,kend
    use local_arrays, only: vx,vy,vz
    use param
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: vof
    real :: alpha
    real :: alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real ,dimension(3) :: Q, C
  
    ! Q stores the cell-centre coordinates and is the query point
    ! C is the control-point and stores the ray origin

    ! Random control-point well-outside computational domain
    C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]
  
    do i = ind(1,1),ind(1,2)
      Q(1) = xm(i)
      do j = ind(2,1),ind(2,2)
          Q(2) = yc(j)
        do kk = ind(3,1),ind(3,2)
             Q(3) = zm(kk)
  
             ! periodic BC
             k = modulo(kk-1,n3m) + 1
  
             if (k.ge.kstart.and.k.le.kend) then 
  
              !call rayTagQ(vof,C,Q,inp)
              call rayTagQ_faster(vof,Q,kk,inp)
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1

               if(VOFy(ii,jj,k).lt.1.0)then
               VOFy(ii,jj,k) = VOFy(ii,jj,k)
               else
                VOFy(ii,jj,k) = vof
               end if

            endif
  
        end do
      end do
    end do
  
  end subroutine

!=================================
!    q3
!=================================

  subroutine convex_hull_q32(ind,inp)
    ! Tagging for z-staggered vz grid
    use mls_param
    use mpih
    use mpi_param, only: kstart,kend
    use local_arrays, only: vx,vy,vz
    use param
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: vof
    real :: alpha
    real :: alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real ,dimension(3) :: Q, C
  
    ! Q stores the cell-centre coordinates and is the query point
    ! C is the control-point and stores the ray origin
  
    ! Random control-point well-outside computational domain
    C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]

    do i = ind(1,1),ind(1,2)
      Q(1) = xm(i)
      do j = ind(2,1),ind(2,2)
          Q(2) = ym(j)
        do kk = ind(3,1),ind(3,2)
             Q(3) = zc(kk)
  
             ! periodic BC
             k = modulo(kk-1,n3m) + 1
  
             if (k.ge.kstart.and.k.le.kend) then 
  
              !call rayTagQ(vof,C,Q,inp)
              call rayTagQ_faster(vof,Q,kk,inp)
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1
              
               if(VOFz(ii,jj,k).lt.1.0)then
               VOFz(ii,jj,k) = VOFz(ii,jj,k)
               else
               VOFz(ii,jj,k) = vof
               end if

               endif
  
        end do
      end do
    end do
  
  end subroutine

  subroutine convex_hull_qc2(ind,inp)
    ! Tagging for cell-centred scalar grid
    use mls_param
    use mpih
    use mpi_param, only: kstart,kend
    use local_arrays, only: temp
    use param
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: vof
    real :: alpha
    real :: alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real ,dimension(3) :: Q, C
    real :: u_imh, u_iph, v_jmh, v_jph, w_kmh, w_kph, h31, h32, h33
    real :: udx1, udx2, udx3
    integer :: kc,km,kp,jm,jc,jp,ic,im,ip

  
    udx1=dx1*0.5d0
    udx2=dx2*0.5d0
    udx3=dx3*0.5d0
  
    ! Q stores the cell-centre coordinates and is the query point
    ! C is the control-point and stores the ray origin

    ! Random control-point well-outside computational domain
    !C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]
  
    do i = ind(1,1),ind(1,2)
      Q(1) = xm(i)
      do j = ind(2,1),ind(2,2)
          Q(2) = ym(j)
        do kk = ind(3,1),ind(3,2)
             Q(3) = zm(kk)
  
             ! periodic BC
             k = modulo(kk-1,n3m) + 1
              
  
             if (k.ge.kstart.and.k.le.kend) then 

              x_gc = pos_cm(1:3,inp)

              !call rayTagQ(vof,C,Q,inp)
              call rayTagQ_faster(vof,Q,kk,inp)
  
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
               VOFp(ii,jj,k) = vof
               end if
               
               if (vof .eq. 0.0) then
                !write(*,*) "Tagging solid cell",ii,jj,k
                solid_mask(ii,jj,k) = .true.

                !                d  u T   |          1   [                              ]
                !             ----------- |  =     ----- |  uT |      -      uT |       |
                !                d   x    |i,j,k     dx  [     i+1/2            i-1/2   ]

                ! uT |_{i-1/2}
                x_grid(1) = xc(i)
                x_grid(2) = ym(j)
                x_grid(3) = zm(kk)
                r = x_grid - x_GC ! relative distance 

                u_imh = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

                ! uT |_{i+1/2}
                x_grid(1) = xc(i+1)
                r = x_grid - x_GC ! relative distance 
                u_iph = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

                h31=( u_iph*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
                -u_imh*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1


                !                d  v T   |          1   [                              ]
                !             ----------- |  =     ----- |  vT |      -      vT |       |
                !                d   y    |i,j,k     dy  [     j+1/2            j-1/2   ] 

                ! vT |_{j-1/2}
                x_grid(1) = xm(i)
                x_grid(2) = yc(j)
                x_grid(3) = zm(kk)
                r = x_grid - x_GC ! relative distance 
                v_jmh = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

                ! vT |_{j+1/2}
                x_grid(2) = yc(j+1)
                r = x_grid - x_GC ! relative distance 
                v_jph = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

                h32=( v_jph*(temp(ic,jp,kc)+temp(ic,jc,kc)) &
                -v_jmh*(temp(ic,jc,kc)+temp(ic,jm,kc)) )*udx2


                !                d  w T   |          1   [                              ]
                !             ----------- |  =     ----- |  wT |      -      wT |       |
                !                d   z    |i,j,k     dz  [     k+1/2            k-1/2   ]

                ! wT |_{k-1/2}
                x_grid(1) = xm(i)
                x_grid(2) = ym(j)
                x_grid(3) = zc(kk)
                r = x_grid - x_GC ! relative distance 
                w_kmh = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

                ! wT |_{k+1/2}
                x_grid(3) = zc(kk+1)
                r = x_grid - x_GC ! relative distance 
                w_kph = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

                h33=( w_kph*(temp(ic,jc,kp)+temp(ic,jc,kc)) &
                -w_kmh*(temp(ic,jc,kc)+temp(ic,jc,km)) )*udx3

                d_UsolidT_dxj(ii,jj,k) = (h31+h32+h33)


               endif
               
               endif
  
        end do
      end do
    end do
  
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

  ! Hard code vertical for testing free-fall
  !bbox_inds(3,1) = 1
  !bbox_inds(3,2) = n3m


  !! Hard-code the full domain for testing
  !bbox_inds(:,1) = [1, 1, 1] 
  !bbox_inds(:,2) = [n1m, n2m, n3m]

end subroutine


subroutine rayTagQ(vof,C,Q,inp)
    ! Tag the cell Q, defined by its cell-centered coordinates as being either interface, fluid or solid
    use param
    use mls_param
    implicit none
    real, dimension(3) :: V0, V1, V2, Q,C, intPoint
    real, dimension(2) :: xx,yy,zz

    real :: vof
    integer :: i,inp, int_count
    logical :: intersect

    int_count = 0

    do i =1,maxnf
        if (.not. isGhostFace(i,inp) ) then
            V0 = xyzv(1:3,vert_of_face(1,i,inp),inp)
            V1 = xyzv(1:3,vert_of_face(2,i,inp),inp)
            V2 = xyzv(1:3,vert_of_face(3,i,inp),inp)

            call rayTriangle_intersect(intersect,intPoint,C,Q,V0,V1,V2)

            ! If intersected, check if cell is interface cell
            if (intersect) then
              if( ( ( Q(1) - 0.5/dx1 .lt. intPoint(1) ) .and.  ( intPoint(1)  .lt.   Q(1) + 0.5/dx1 )  ) .and. &
                ( ( Q(2) - 0.5/dx2 .lt. intPoint(2) ) .and.  ( intPoint(2)  .lt.   Q(2) + 0.5/dx2 )  ) .and. &
                 ( ( Q(3) - 0.5/dx3 .lt. intPoint(3) ) .and.  ( intPoint(3)  .lt.   Q(3) + 0.5/dx3 )  ) ) then

                !! Cell is interface, estimate VOF based on signed distance to the triangle plane
                !xx = [ Q(1) - 0.5/dx1 , Q(1) + 0.5/dx1 ]
                !yy = [ Q(2) - 0.5/dx2 , Q(2) + 0.5/dx2 ]
                !zz = [ Q(3) - 0.5/dx3 , Q(3) + 0.5/dx3 ]
                !call compute_interface_VOF(vof,xx,yy,zz,tri_nor(1:3,i,inp),V0)
                vof = 0.5

                ! Exit, no need to check other triangles
                return

              else ! Intersected: add to counter
                int_count = int_count + 1
              endif ! if interface
            endif !if intersect
        endif !ifGhost
    enddo

    ! Even intersections: external, odd intersections: internal
    if (modulo(int_count,2) .eq. 0) then !even = fluid
        vof = 1.0
    else ! odd = solid
        vof = 0.0
    endif

end subroutine

subroutine rayTagQ_faster(vof,Q,k,inp)
  ! Tag the cell Q, defined by its cell-centered coordinates as being either interface, fluid or solid
  use param
  use mls_param
  implicit none
  real, dimension(3) :: V0, V1, V2, Q,C, intPoint
  real, dimension(2) :: xx,yy,zz

  real :: vof, zcent
  integer :: i,inp, int_count, pad, k1,k
  logical :: intersect, insideBox

  ! Force C to be parallel to the computational grid lines: efficient

  C = [-xlen*10.0, -ylen*3.0, Q(3) ]

  int_count = 0
  
  pad = 2

  do i =1,maxnf
      if (.not. isGhostFace(i,inp) ) then
          V0 = xyzv(1:3,vert_of_face(1,i,inp),inp)
          V1 = xyzv(1:3,vert_of_face(2,i,inp),inp)
          V2 = xyzv(1:3,vert_of_face(3,i,inp),inp)

        ! If triangle in cell, cell is interface, break out
        call triangleBoxIntersect(insideBox,V0,V1,V2, Q, dx1, dx2, dx3)    
        if (insideBox .eqv. .true. ) then
          vof = 0.5
          return
        endif       

        zcent = (1.0/3.0)* ( V0(3) + V1(3) + V2(3) )
        k1 = floor(zcent*dx3) + 1

        if (abs(k1 - k) .le. pad ) then
          call rayTriangle_intersect(intersect,intPoint,C,Q,V0,V1,V2)

          ! If intersected, check if cell is interface cell
          if (intersect) then
            call pointBoxIntersect(insideBox,intPoint, Q, dx1, dx2, dx3)
            if (insideBox .eqv. .true. ) then

              !! Cell is interface, estimate VOF based on signed distance to the triangle plane
              !xx = [ Q(1) - 0.5/dx1 , Q(1) + 0.5/dx1 ]
              !yy = [ Q(2) - 0.5/dx2 , Q(2) + 0.5/dx2 ]
              !zz = [ Q(3) - 0.5/dx3 , Q(3) + 0.5/dx3 ]
              !call compute_interface_VOF(vof,xx,yy,zz,tri_nor(1:3,i,inp),V0)
              vof = 0.5

              ! Exit, no need to check other triangles
              return

            else ! Intersected: add to counter
              int_count = int_count + 1
            endif ! if interface
          endif !if intersect
        endif
      endif !ifGhost
  enddo

  ! Even intersections: external, odd intersections: internal
  if (modulo(int_count,2) .eq. 0) then !even = fluid
      vof = 1.0
  else ! odd = solid
      vof = 0.0
  endif

end subroutine


recursive subroutine rayTriangle_intersect(intersect,intPoint,C,Q,V0,V1,V2)
    ! Using the Möller & Trumbore (2005) algorithm to
    ! test for a ray-triangle intersection where the ray is cast from an
    ! origin C to the query point Q
    ! The triangle is defined by its three vertex coordinates V0,V1,V2

    use param
    !use mls_param
    implicit none
    real, dimension(3) :: V0, V1, V2, Q,C,Cprime, r, e1, e2, s, p, q1, intPoint
    real :: EPS, a, f, u , v, t
    logical :: intersect

    EPS = 1e-14

    intersect = .true.

    ! Ray direction vector
    r = Q - C
    
    ! Edge-vectors of triangle
    e1 = V1 - V0
    e2 = V2 - V0

    s = c - V0
    call cross(p, r, e2)

    a = dot_product(e1,p)

    !Parallel ray, change query point and recursively call 
    !if ( (a .gt. -EPS) .and. (a .le. EPS ) ) then
    if (abs(a) .lt. EPSILON(1.0d0) ) then
        Cprime = [C(2), C(3), C(1) ]
        write(*,*) "Warning: parallel ray detected!"
        call rayTriangle_intersect(intersect,intPoint,Cprime,Q,V0,V1,V2)
    endif

    f = 1.0d0 / a

    ! First Barycentric limit u[0,1]
    u = f * dot_product(s,p) 
    if ( (u < EPSILON(1.0d0) ) .or. (u .gt. 1.0) ) then
        intersect = .false.
    endif

    ! Second barycentric limit v[0,1] u + v <= 1
    call cross(q1, s, e1)
    v = f * dot_product(r, q1)

    if (  (v .lt. EPSILON(1.0d0)) .or. ( (u + v) .gt. 1.0  ) ) then
        intersect = .false.
    endif

    ! Line parameter
    ! Note in general: if above two u,v conditions were passed but below t is not in [0,1]
    ! This means the extruded ray will intersect the triangle, but not in the range specified by C-->Q

    t = f * dot_product(e2, q1)
    if ( (t .lt. 0) .or. (t .gt. 1.0) ) then
        intersect = .false.
    endif

    ! If intersected, compute the intersection point to be returned
    if (intersect) then
        intPoint = C + t * r
    endif

end subroutine

subroutine pointBoxIntersect(insideBox,intPoint, Q, dx1, dx2, dx3)
  ! Test if point intPoint is inside 3D box with centroid Q
  !use mls_param
  implicit none
  real, dimension(3) :: Q, intPoint
  real :: dx1, dx2, dx3
  logical :: insideBox

  insideBox = .false.


  if( ( ( Q(1) - 0.5/dx1 .lt. intPoint(1) ) .and.  ( intPoint(1)  .lt.   Q(1) + 0.5/dx1 )  ) .and. &
  ( ( Q(2) - 0.5/dx2 .lt. intPoint(2) ) .and.  ( intPoint(2)  .lt.   Q(2) + 0.5/dx2 )  ) .and. &
   ( ( Q(3) - 0.5/dx3 .lt. intPoint(3) ) .and.  ( intPoint(3)  .lt.   Q(3) + 0.5/dx3 )  ) ) then
    insideBox = .true.
   endif

end subroutine pointBoxIntersect

subroutine triangleBoxIntersect(insideBox,V0,V1,V2, Q, dx1, dx2, dx3)
  ! Test if point intPoint is inside 3D box with centroid Q
  !use mls_param
  implicit none
  real, dimension(3) :: Q, V0, V1, V2
  real :: dx1, dx2, dx3
  logical :: insideBox

  insideBox = .false.

  call pointBoxIntersect(insideBox,V0, Q, dx1, dx2, dx3)

  if (insideBox .eqv. .false.) then
    call pointBoxIntersect(insideBox,V1, Q, dx1, dx2, dx3)
  endif

  if (insideBox .eqv. .false.) then
    call pointBoxIntersect(insideBox,V2, Q, dx1, dx2, dx3)
  endif


end subroutine triangleBoxIntersect

subroutine interp_vof
  ! Interpolate VOFx, VOFv, VOFw from cell-centered VOFp information
  use param
  use mpi_param, only: kstart,kend
  implicit none
  integer :: ic,jc,kc
  integer :: km,kp,jm,jp,im,ip

  !call update_both_ghosts(n1,n2,VOFx,kstart,kend)
  !call update_both_ghosts(n1,n2,VOFy,kstart,kend)
  !call update_both_ghosts(n1,n2,VOFz,kstart,kend)
  call update_both_ghosts(n1,n2,VOFp,kstart,kend)

  do kc=kstart,kend
    km=kc-1
    kp=kc+1
      do jc=1,n2m
        jm=jmv(jc)
        jp=jpv(jc)
        do ic=1,n1m
          ip=ipv(ic)
          im=imv(ic)
          
          VOFx(ic,jc,kc) = 0.5 * (VOFp(im,jc,kc) + VOFp(ic,jc,kc) )
          VOFy(ic,jc,kc) = 0.5 * (VOFp(ic,jm,kc) + VOFp(ic,jc,kc) )
          VOFz(ic,jc,kc) = 0.5 * (VOFp(ic,jc,km) + VOFp(ic,jc,kc) )
        enddo
      enddo
  enddo
end subroutine interp_vof

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
 
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

            call rayTagQ(vof,C,Q,inp)

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
  
              call rayTagQ(vof,C,Q,inp)
  
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
  
              call rayTagQ(vof,C,Q,inp)
  
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
             Q(3) = zm(kk)
  
             ! periodic BC
             k = modulo(kk-1,n3m) + 1
  
             if (k.ge.kstart.and.k.le.kend) then 
  
              call rayTagQ(vof,C,Q,inp)
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1

               if(VOFp(ii,jj,k).lt.1.0)then
               VOFp(ii,jj,k) = VOFp(ii,jj,k)
               else
               VOFp(ii,jj,k) = vof
               end if
               
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

    EPS = 1e-10

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
    if ( (a .gt. -EPS) .and. (a .le. EPS ) ) then
        Cprime = [C(2), C(3), C(1) ]
        call rayTriangle_intersect(intersect,intPoint,Cprime,Q,V0,V1,V2)
    endif

    f = 1.0d0 / a

    ! First Barycentric limit u[0,1]
    u = f * dot_product(s,p) 
    if ( (u < EPS ) .or. (u .gt. 1.0) ) then
        intersect = .false.
    endif

    ! Second barycentric limit v[0,1] u + v <= 1
    call cross(q1, s, e1)
    v = f * dot_product(r, q1)

    if (  (v .lt. EPS) .or. ( (u + v) .gt. 1.0  ) ) then
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

subroutine compute_interface_VOF(VOF)
  !use param
  !use mls_param
  !use geom
  implicit none

  real :: VOF
  VOF = 0.5

end subroutine



! subroutine compute_interface_VOF(VOF,xx,yy,zz,nhat,V0)
!     !use param
!     !use mls_param
!     use geom
!     implicit none

!     real :: sumPhi, alpha, phi,VOF
!     real, dimension(2) :: xx, yy, zz
!     real, dimension(3) :: x0, V0, nhat
!     integer :: i,j,k

! ! Estimate volume-of-fluid of an interface cell based on the signed distance at each cell corner
! ! Signed distance computed as the signed distance to the triangle plane

!     sumPhi = 0.0
!     alpha = 0.0
!     ! Loop over cell corners
!     do i = 1,2
!         x0(1) = xx(i)
!         do j = 1,2
!             x0(2) = yy(j)
!             do k = 1,2
!                 x0(3) = zz(k)

!                 phi = -dot_product(nhat,x0 - V0)
!                 alpha = alpha - phi * heaviside(-phi)
!                 sumPhi = sumPhi + abs(phi)
!             enddo
!         enddo
!     enddo

!     alpha = alpha / sumPhi
!     VOF = 1.0 - alpha
! end subroutine
! Auxilary helper routines for ray-tagging
! Mostly contains primitive geometric operations like intersection-testing, etc
! with the help of chatGPT here on edge-intersection testing

module rayAux
  implicit none
  contains

  recursive subroutine rayTriangle_intersect(intersect,intPoint,C,Q,V0,V1,V2)
    ! Using the Möller & Trumbore (2005) algorithm to
    ! test for a ray-triangle intersection where the ray is cast from an
    ! origin C to the query point Q
    ! The triangle is defined by its three vertex coordinates V0,V1,V2
    use param, only: dx1, dx2, dx3
    !use mls_param
    implicit none
    real, dimension(3) :: V0, V1, V2, Q,C,Cprime, r, e1, e2, s, p, q1, intPoint
    real ::  a, f, u , v, t
    logical :: intersect

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
        Cprime = [C(1) , C(2) + 0.1/dx2, C(3) -0.1/dx3 ] ! Slightly perturb the C-point: still contained in pencil, but centre-aligned
        call rayTriangle_intersect(intersect,intPoint,Cprime,Q,V0,V1,V2)
        write(*,*) "Warning: parallel ray detected!"
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
    if ( (t .lt. 0.0 ) .or. (t .gt. 1.0) ) then
        intersect = .false.
    endif

    ! If intersected, compute the intersection point to be returned
    if (intersect) then
        intPoint = C + t * r
    endif

end subroutine

  
    subroutine intersects_yz_range(V0, V1, V2, ymin, ymax, zmin, zmax, intersects)
      ! Evaluate whether a triangle (V0,V1,V2) intersects with a [y,z]-range
      ! Used for x-aligned pencils which are defined by their [ymin,ymax], [zmin, zmax] square
      ! To test whether this intersection occurs, we test
      !   1) if the vertices of the triangle are inside the square
      !   2) if the bounding box of the triangle intersects with the square
      !   3) if the triangle edges intersect with the square edges
      real, dimension(3), intent(in) :: V0, V1, V2   ! Vertices of the triangle in 3D space (x, y, z)
      real, intent(in) :: ymin, ymax                  ! y bounds
      real, intent(in) :: zmin, zmax                  ! z bounds
      logical, intent(out) :: intersects                  ! Output: True if the triangle intersects the y-z range
  
      real, dimension(3) :: yvals, zvals              ! y and z coordinates of the triangle's vertices
  
      ! Extract y and z coordinates of the triangle vertices
      yvals = [ V0(2), V1(2), V2(2) ]
      zvals = [ V0(3), V1(3), V2(3) ]
  
      ! Check if any vertex of the triangle lies inside the y-z range
      if ( any((yvals >= ymin) .and. (yvals <= ymax) .and. (zvals >= zmin) .and. (zvals <= zmax)) ) then
        intersects = .true.
        return
      endif
  
      ! Check if the bounding box of the triangle in y-z intersects with the y-z range
      if (maxval(yvals) < ymin .or. minval(yvals) > ymax .or. maxval(zvals) < zmin .or. minval(zvals) > zmax) then
        intersects = .false.
        return
      endif
  
      ! Perform a more detailed intersection check for edges of the triangle with the bounding box edges
      call triangle_square_edges_intersection(yvals, zvals, ymin, ymax, zmin, zmax, intersects)
  
    end subroutine intersects_yz_range
  
    subroutine triangle_square_edges_intersection(yvals, zvals, ymin, ymax, zmin, zmax, intersects)
      ! This subroutine checks if any edge of the triangle intersects with the bounding square in the y-z plane
      real, dimension(3), intent(in) :: yvals, zvals
      real, intent(in) :: ymin, ymax, zmin, zmax
      logical, intent(out) :: intersects
      integer :: i, j
      real :: y1, y2, z1, z2
      logical :: e_int_box
  
      intersects = .false.
  
      ! Check all three edges of the triangle
      do i = 1, 3
        j = mod(i, 3) + 1
        y1 = yvals(i)
        y2 = yvals(j)
        z1 = zvals(i)
        z2 = zvals(j)
  
        ! Check if the edge intersects with the bounding square
        call edge_intersects_square(y1, z1, y2, z2, ymin, ymax, zmin, zmax, e_int_box)
  
        if (e_int_box) then
          intersects = .true.
          return
        endif
      enddo
    end subroutine triangle_square_edges_intersection
  
    subroutine edge_intersects_square(y1, z1, y2, z2, ymin, ymax, zmin, zmax, intersects)
      ! Test generated with the help of chatGPT
      ! Check if a line segment (y1, z1) to (y2, z2) intersects with the square defined by [ymin, ymax] and [zmin, zmax]
      real, intent(in) :: y1, z1, y2, z2, ymin, ymax, zmin, zmax
      logical, intent(out) :: intersects
  
      real :: t
      real :: y, z
  
      intersects = .false.
  
      ! Check against each side of the box
  
      ! Left (y = ymin) and Right (y = ymax)
      if (y1 .ne. y2) then
        t = (ymin - y1) / (y2 - y1)
        if (t .ge. 0.0 .and. t .le. 1.0 ) then
          z = z1 + t * (z2 - z1)
          if (z .ge. zmin .and. z .le. zmax) then
            intersects = .true.
            return
          endif
        endif
  
        t = (ymax - y1) / (y2 - y1)
        if (t .ge. 0.0 .and. t .le. 1.0 ) then
          z = z1 + t * (z2 - z1)
          if (z .ge. zmin .and. z .le. zmax) then
            intersects = .true.
            return
          endif
          endif
      endif
  
      ! Bottom (z = zmin) and Top (z = zmax)
      if (z1 .ne. z2) then
        t = (zmin - z1) / (z2 - z1)
        if (t .ge. 0.0 .and. t .le. 1.0 ) then
          y = y1 + t * (y2 - y1)
          if (y .ge. ymin .and. y .le. ymax) then
            intersects = .true.
            return
          endif
        endif
  
        t = (zmax - z1) / (z2 - z1)
        if (t .ge. 0.0 .and. t .le. 1.0 ) then
          y = y1 + t * (y2 - y1)
          if (y .ge. ymin .and. y .le. ymax) then
            intersects = .true.
            return
          endif
        endif
      endif
    end subroutine edge_intersects_square

    subroutine triangleXIntersect(xintersect,V0,V1,V2, xmin, xmax )
      ! Test if triangle (V0,V1,V2) is inside within (xmin,xmax)
      ! This supplements the pencil-intersection routine, which has already tested for [ymin,ymax] [zmin, zmax] intersection
      ! Pencil(j,k) is defined spatially by its cell-centred coordinate, XC_jk which stores the (ym,zm) coordinates
      implicit none
      real, dimension(3) :: V0, V1, V2
      real ::   xmin_tri, xmax_tri, xmin, xmax
      logical :: xintersect, point_in_range
    
      xintersect = .false.
  
      xmin_tri = min( V0(1), V1(1), V2(1) )
      xmax_tri = max( V0(1), V1(1), V2(1) )
  
      ! Check if the bounding box of the triangle overlaps
      if (xmax_tri < xmin .or. xmin_tri > xmax) return  ! No overlap in x-direction
  
  
      ! Check if any vertex of the triangle is inside the 2D range
      point_in_range = .false.
      point_in_range = point_in_range .or. (V0(1) >= xmin .and. V0(1) <= xmax )
      point_in_range = point_in_range .or. (V1(1) >= xmin .and. V1(1) <= xmax )
      point_in_range = point_in_range .or. (V2(1) >= xmin .and. V2(1) <= xmax )
  
      if (point_in_range) then
          xintersect = .true.
          return
      end if
  
  end subroutine triangleXIntersect


  subroutine pointBoxIntersect(insideBox,intPoint, xmin, xmax, ymin, ymax, zmin, zmax)
      ! Test point is inside a box
      implicit none
      real, dimension(3) :: intPoint
      real :: xmin, xmax, ymin, ymax, zmin, zmax
      logical :: insideBox
    
  
      if (intPoint(1) >= xmin .and. intPoint(1) <= xmax .and. &  ! x-direction check
          intPoint(2) >= ymin .and. intPoint(2) <= ymax .and. &  ! y-direction check
          intPoint(3) >= zmin .and. intPoint(3) <= zmax) then    ! z-direction check
  
          insideBox = .true.
      else
          insideBox = .false.
      end if
  
  end subroutine pointBoxIntersect

  
end module rayAux
  
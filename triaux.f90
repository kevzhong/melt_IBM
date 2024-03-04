!------------------------------------------------------
subroutine calculate_volume2 (Volume,nf,tri_nor,sur,tri_bar,isGhostFace)
! Evaluate the triangulated surface volume using the divergence theorem
! using the (pre-)computed triangle areas and (outward-facing) normal vectors
!
! ---     _           ---  _   _
! /// div F dV    =    //  F . n  dA
! ---                 ---
!  V                   S
!                _                            _
! Here we choose F = [x, 0, 0]' such that div F = 1, yielding the desired volume integral on the LHS and the RHS becomes
!
!       ---  _   _                ---                                    ---
!        //  F . n  dA    =       //  [x, 0, 0]' .  [nx, ny, nz]' dA =   //  x * nx  dA
!       ---                      ---                                    ---
!        S                        S                                      S
!
! RHS surface integral is then evaluated numerically as
!
! Nfaces
! ____
! \
!      |   nx_i  * x_i * A_i   |
! /
! ----  
! i=1
!
! Where x_i is taken to be the centroid of triangle i
!       A_i is the area of triangle i
!       nx_i is the x-component of the outward-facing normal vector for triangle FACE i
! Note the absolute value for the summation argument, this seems necessary
! Probably related to how the definition of outward-facing changes when we go from positive x-axis to negative x-axis

implicit none
integer :: nf,i
real,dimension (3,nf) :: tri_nor
real, dimension (3,nf) :: tri_bar
real, dimension (nf) :: sur
logical, dimension(nf) :: isGhostFace
real :: Volume
        
Volume=0.0
do i=1,nf
        if (isGhostFace(i) .eqv. .false. ) then
                Volume = Volume + tri_nor(1,i) * tri_bar(1,i) * sur(i) 
        endif
enddo        
return

end subroutine calculate_volume2
!------------------------------------------------------

subroutine calculate_volume (Volume,nv,nf,xyz,vert_of_face,vol)

implicit none
integer :: nv,nf,v1,v2,v3,i
integer, dimension (3,nf) :: vert_of_face
real, dimension (3,nv) ::xyz
real, dimension (nf) :: vol
real :: Volume

Volume=0.0
do i=1,nf
        v1=vert_of_face(1,i)
        v2=vert_of_face(2,i)
        v3=vert_of_face(3,i)

        vol(i) =           xyz(1,v1) * (xyz(2,v2)*xyz(3,v3) - xyz(3,v2)*xyz(2,v3)) &
                        +  xyz(1,v2) * (xyz(2,v3)*xyz(3,v1) - xyz(3,v3)*xyz(2,v1)) &
                        +  xyz(1,v3) * (xyz(2,v1)*xyz(3,v2) - xyz(3,v1)*xyz(2,v2))
        Volume=Volume+vol(i)
enddo
 
Volume=Volume/6.

return
end subroutine calculate_volume
!------------------------------------------------------
subroutine calculate_area (Surface,nv,nf,xyz,vert_of_face,sur,isGhostFace,rm_flag,A_thresh)

implicit none
integer :: nv,nf,v1,v2,v3,i
integer, dimension (3,nf) :: vert_of_face
logical, dimension(nf) :: isGhostFace
real, dimension (3,nv) ::xyz
real, dimension (nf) :: sur
real :: Surface,d12,d23,d31,sp
logical :: rm_flag
real :: A_thresh

Surface=0.0
do i=1,nf
        if (isGhostFace(i) .eqv. .false. ) then
                v1=vert_of_face(1,i)
                v2=vert_of_face(2,i)
                v3=vert_of_face(3,i)

                d12=sqrt( (xyz(1,v1)-xyz(1,v2))**2. &
                         +(xyz(2,v1)-xyz(2,v2))**2. &
                         +(xyz(3,v1)-xyz(3,v2))**2. )
                d23=sqrt( (xyz(1,v2)-xyz(1,v3))**2. &
                         +(xyz(2,v2)-xyz(2,v3))**2. &
                         +(xyz(3,v2)-xyz(3,v3))**2. )
                d31=sqrt( (xyz(1,v3)-xyz(1,v1))**2. &
                         +(xyz(2,v3)-xyz(2,v1))**2. &
                         +(xyz(3,v3)-xyz(3,v1))**2. )
                sp=(d12+d23+d31)/2.
                sur(i)=sqrt(sp*(sp-d12)*(sp-d23)*(sp-d31))
                Surface=Surface+sur(i)

                if (sur(i) .le. A_thresh) then
                        rm_flag = .true. ! set the remesh flag to true if any triangle less than the threshold area
                endif
        endif
enddo

return
end subroutine calculate_area

subroutine calculate_vert_area (Avert,nv,nf,vert_of_face,sur,isGhostFace)
        implicit none
integer :: nv,nf,v1,v2,v3,i
integer, dimension (3,nf) :: vert_of_face
logical, dimension(nf) :: isGhostFace
real, dimension (nv) ::Avert
real, dimension (nf) :: sur
        ! Vertex area defined as the total area of faces adjoining a vertex

Avert = 0.

do i = 1,nf
        if (isGhostFace(i) .eqv. .false. ) then
                v1=vert_of_face(1,i)
                v2=vert_of_face(2,i)
                v3=vert_of_face(3,i)

                Avert(v1) = Avert(v1) + sur(i)
                Avert(v2) = Avert(v2) + sur(i)
                Avert(v3) = Avert(v3) + sur(i)
        endif
enddo

return

end subroutine calculate_vert_area


!------------------------------------------------------
subroutine calculate_normal (tri_nor,nv,nf,xyz,vert_of_face)

implicit none
integer :: nv,nf,v1,v2,v3,i
integer, dimension (3,nf) :: vert_of_face
real,dimension (3,nf) :: tri_nor
real, dimension (3,nv) ::xyz
real, dimension(3) :: ve1,ve2

do i=1,nf
        v1=vert_of_face(1,i)
        v2=vert_of_face(2,i)
        v3=vert_of_face(3,i)

        ! Vectors ve1 = P2 - P1   ve2 = P3 - P1
        ve1(1:3) = xyz(1:3,v2) - xyz(1:3,v1)
        ve2(1:3) = xyz(1:3,v3) - xyz(1:3,v1)

        ! Compute cross-product cross(ve1, ve2)
        tri_nor(1,i) = ve1(2)*ve2(3) - ve1(3)*ve2(2)
        tri_nor(2,i) = ve1(3)*ve2(1) - ve1(1)*ve2(3)
        tri_nor(3,i) = ve1(1)*ve2(2) - ve1(2)*ve2(1)
        tri_nor(1:3,i) = tri_nor(1:3,i) / sqrt ( sum ( tri_nor(1:3,i)**2  )  )
enddo

return

end subroutine calculate_normal

subroutine update_tri_normal (tri_nor,nv,nf,xyz,vert_of_face,isGhostFace)
        ! "Safely" update the triangle normal vector using information from prev. iteration/timestep normal vector
        ! Take sign of dot product with newly-computed normal vector and previous normal vector
        ! Choose the one that is best-aligned with prev iteration 
        ! i.e. choose nhat_new such that sign( dot(nhat_old, nhat_new) ) = 1

        implicit none
        integer :: nv,nf,v1,v2,v3,i
        integer, dimension (3,nf) :: vert_of_face
        logical, dimension(nf) :: isGhostFace
        real,dimension (3,nf) :: tri_nor
        real, dimension (3,nv) ::xyz
        real, dimension(3) :: ve1,ve2
        real, dimension(3) :: nhat_old
        real :: sgn

        do i=1,nf
                if ( isGhostFace(i) .eqv. .false. ) then
                        nhat_old(1:3) = tri_nor(1:3,i)
                        
                        v1=vert_of_face(1,i)
                        v2=vert_of_face(2,i)
                        v3=vert_of_face(3,i)
                        
                        ! Vectors ve1 = P2 - P1   ve2 = P3 - P1
                        ve1(1:3) = xyz(1:3,v2) - xyz(1:3,v1)
                        ve2(1:3) = xyz(1:3,v3) - xyz(1:3,v1)
                        
                        ! Compute cross-product cross(ve1, ve2)
                        tri_nor(1,i) = ve1(2)*ve2(3) - ve1(3)*ve2(2)
                        tri_nor(2,i) = ve1(3)*ve2(1) - ve1(1)*ve2(3)
                        tri_nor(3,i) = ve1(1)*ve2(2) - ve1(2)*ve2(1)
                        tri_nor(1:3,i) = tri_nor(1:3,i) / sqrt ( sum ( tri_nor(1:3,i)**2  )  )
                        
                        ! flip if needed
                        sgn = dot_product( tri_nor(1:3,i) , nhat_old(1:3) ) 
                        sgn = sign(1.0, sgn)
                        
                        tri_nor(1:3,i) = tri_nor(1:3,i)*sgn
                endif
        enddo
        
        return
        
end subroutine update_tri_normal


!------------------------------------------------------

subroutine calculate_areaWeighted_vert_normal (tri_nor,vert_nor,nv,nf,Atri,vert_of_face,isGhostFace,isGhostVert)

! This routine computes the (normalised) normal vectors at each of the geometric vertices
! Computed as area-weighted average (from triangle areas Atri)

implicit none
integer :: i,nv,nf,v
real,dimension (3,nf) :: tri_nor
real,dimension (3,nv) :: vert_nor
logical, dimension(nf) :: isGhostFace
logical, dimension(nv) :: isGhostVert
real, dimension(nf) :: Atri
real :: magN
integer, dimension(3,nf) :: vert_of_face

vert_nor(:,:) = 0.

do i=1,nf
        if (isGhostFace(i) .eqv. .false. ) then
                do v=1,3 ! For each 3 vertices of each triangle
                        vert_nor( 1:3, vert_of_face(v,i) ) = vert_nor( 1:3, vert_of_face(v,i) ) + Atri(i)*tri_nor(1:3,i)
                enddo
        endif
enddo

! No need to track area of vertices, since we normalise in the end
! Normalise for all vertices
do i=1,nv
        if (isGhostVert(i) .eqv. .false. ) then
                magN = norm2 ( vert_nor(:,i) )
                vert_nor(1:3,i) = vert_nor(1:3,i) / magN
        endif
enddo


return

end subroutine calculate_areaWeighted_vert_normal
!------------------------------------------------------

subroutine calculate_eLengths(eLengths,nv,ne,xyz,vert_of_edge,isGhostEdge)

implicit none
integer :: nv,ne,v1,v2,i
integer, dimension (2,ne) :: vert_of_edge
logical, dimension(ne) :: isGhostEdge
real, dimension (3,nv) ::xyz
real, dimension (ne) :: eLengths
do i=1,ne
        if ( isGhostEdge(i) .eqv. .false. ) then
                v1=vert_of_edge(1,i)
                v2=vert_of_edge(2,i)
                eLengths(i)=sqrt( (xyz(1,v1)-xyz(1,v2))**2 + (xyz(2,v1)-xyz(2,v2))**2 + (xyz(3,v1)-xyz(3,v2))**2 ) 

        endif
enddo

return
end subroutine calculate_eLengths

!     ----------------------------------------------------------------
subroutine calc_centroids_from_vert(tri_cent,xyz,vert_of_face,nf,nv,isGhostFace)
        implicit none
        integer :: i, nf, nv
        real, dimension(3,nf) :: tri_cent
        real, dimension(3,nv) :: xyz
        logical, dimension(nf) :: isGhostFace
        integer, dimension(3, nf) :: vert_of_face 

        do i = 1,nf
                if (isGhostFace(i) .eqv. .false.) then
                        tri_cent(1,i) =  sum ( xyz( 1, vert_of_face(1:3,i)  )  )  / 3.0d0
                        tri_cent(2,i) =  sum ( xyz( 2, vert_of_face(1:3,i)  )  )  / 3.0d0
                        tri_cent(3,i) =  sum ( xyz( 3, vert_of_face(1:3,i)  )  )  / 3.0d0
                endif
        enddo
end subroutine

subroutine calculate_skewness (ne,nf,edge_of_face,sur,eLengths,skewness,isGhostFace,rm_flag,skew_thresh)

        implicit none
        integer :: ne,nf,e1,e2,e3,i
        integer, dimension (3,nf) :: edge_of_face
        logical, dimension(nf) :: isGhostFace
        real, dimension (nf) :: skewness, sur
        real, dimension(ne) :: eLengths
        real :: perim, sur_opt
        logical :: rm_flag
        real :: skew_thresh

        do i=1,nf
                if (isGhostFace(i) .eqv. .false. ) then
                        ! Compute triangle perimeter
                        e1=edge_of_face(1,i)
                        e2=edge_of_face(2,i)
                        e3=edge_of_face(3,i)

                        perim = eLengths(e1) + eLengths(e2) + eLengths(e3)

                        ! Equilateral triangle area with matched perimeter
                        sur_opt = sqrt(3.0) / 4.0 * (perim / 3.0)**2

                        skewness(i) = ( sur_opt - sur(i) ) / sur_opt

                        if (skewness(i) .gt. skew_thresh) then
                                rm_flag = .true. ! set the remesh flag to true if any edgeLength less than the threshold length
                        endif
                endif
        enddo
        
        return
end subroutine calculate_skewness

!------------------------------------------------------


subroutine dot(xy,x,y)

implicit none
real x(3),y(3),xy

xy = x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
return
end subroutine dot

!     ----------------------------------------------------------------

subroutine  sub(xmy,x,y)

implicit none
real x(3),y(3),xmy(3)

xmy = x-y
return
end subroutine sub

!     ----------------------------------------------------------------

subroutine  cross(xcy,x,y)
implicit none
real x(3),y(3),xcy(3)

xcy(1) = x(2)*y(3)-y(2)*x(3)
xcy(2) = x(3)*y(1)-y(3)*x(1)
xcy(3) = x(1)*y(2)-y(1)*x(2)
return
end subroutine cross

!     ----------------------------------------------------------------

subroutine inverseLU(a,c)
implicit none
real :: a(4,4), c(4,4)
real :: L(4,4), U(4,4)
real :: b(4), d(4), x(4)
real :: coeff
integer i,j,k

!L = 0.0 ; U = 0.0 ; b = 0.0

!forward elimination
do k=1,3
 do i=k+1,4
  coeff=a(i,k)/a(k,k)
  L(i,k) = coeff
  do j=k+1,4
   a(i,j)=a(i,j)-coeff*a(k,j)
  end do
 end do
end do

!prepare L U
do i=1,4
 L(i,i) = 1.0
end do
do j=1,4
 do i=1,j
  U(i,j) = a(i,j)
 end do
end do

!compute columns of c
do k=1,4
 b(k)=1.0
 d(1)=b(1)
 do i=2,4
  d(i)=b(i)
  do j=1,i-1
   d(i) = d(i)-L(i,j)*d(j)
  end do
 end do
 !solve ux=d with back subs.
 x(4)=d(4)/U(4,4)
 do i=3,1,-1
  x(i) = d(i)
   do j=4,i+1,-1
    x(i)=x(i)-U(i,j)*x(j)
   end do
   x(i) = x(i)/u(i,i)
  end do
 !fill solns of x(n) to k of C
  do i=1,4
   c(i,k)=x(i)
  end do
  b(k) = 0.0
end do

return
end


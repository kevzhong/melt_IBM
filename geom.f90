module geom
use param
use mls_param
implicit none

contains

function delta(r) result(phi)
implicit none
real :: r, phi

   if (r.le.1.5 .and. r.ge.0.5) then
     phi = 5.0 - 3.0*abs(r) - sqrt(-3.0*(1.0 - abs(r))**2 + 1.0)
     phi = phi / 6.0
   elseif (r.le. 0.5 .and. r.ge.0.0) then
     phi = 1.0 + sqrt(-3.0*r**2 + 1.0)
     phi = phi / 3.0
   else
     phi = 0.0
   endif

end



function loopoverbeams(x0,x_cm,AA,inp) result(phi)
implicit none
integer :: i,inp
real    :: mindist,t,d,r
real    :: y1(3,3), y2(3,3), AA(3,3), AAT(3,3), AAT_P(3,3)
real, dimension(3) :: x0,x1,x2,x_cm
real    :: phi

r = .544279548

y1 = 0.
y2 = 0.
!chiral
! particle 1
!y1(1,1) = -r
!y2(1,1) =  r
! particle 2
!y1(1,2) = -r
!y2(1,2) = -r
!y2(2,2) = r
!y2(3,2) = 0.942719830656624
! particle 3
!y1(1,3) = r
!y2(1,3) = r
!y2(2,3) = r
!y2(3,3) = -0.942719830656624

!chiral
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


!achiral
! particle 1
!y1(1,1) = -r
!y2(1,1) =  r
! particle 2
!y1(1,2) = -r
!y2(1,2) = -r
!y2(3,2) = -2.*r
! particle 3
!y1(1,3) = r
!y2(1,3) = r
!y2(2,3) = 2.*r

! moviing centroid such that volumetric centre is at (0,0,0)
y1(2,:) = y1(2,:) - 0.2125855
y2(2,:) = y2(2,:) - 0.2125855

y1(3,:) = y1(3,:) + 0.2125855
y2(3,:) = y2(3,:) + 0.2125855



mindist = 1.e6

AAT   = transpose(AA)
AAT_P = princ_axis_rotm()

do i = 1,3  
  x1 = matmul(AAT_P, y1(:,i))
  x2 = matmul(AAT_P, y2(:,i))
  x1 = matmul(AAT,x1) + x_cm 
  x2 = matmul(AAT,x2) + x_cm 

  t = step(x0,x1,x2)

  if (t.ge.0 .and. t.le.1) then
   d = shortdist(x0,x1,x2)
   mindist = min(d, mindist)
  endif
  ! spherical cap
    d = norm2(x1-x0)
    mindist = min(d, mindist)
    d = norm2(x2-x0)
    mindist = min(d, mindist)

enddo

  phi = 1 - mindist/(0.4*r) 


end

function shortdist(x0,x1,x2) result(d)
implicit none
real, dimension(3) :: x0, x1, x2, temp
real               :: d

call cross(temp, x2-x1, x1-x0)

d = norm2(temp) / norm2(x2-x1)


end

function step(x0,x1,x2) result(t)
implicit none
real, dimension(3) :: x0, x1, x2
real :: t

  t =  dot_product(x0-x1, x2-x1) / norm2(x2-x1)**2
end

function princ_axis_rotm() result(AAT)
 implicit none
 real :: a,b,c,q(4),AAT(3,3)
!chiral
!  AAT(1,1)=  -0.719319343
!  AAT(1,2)=  -0.0476315399
!  AAT(1,3)=  -0.693044673

!  AAT(2,1)=  -0.0324121773
!  AAT(2,2)=   0.998861272
!  AAT(2,3)=  -0.0350087021

!  AAT(3,1)=   0.693923002
!  AAT(3,2)=  -0.00271934974
!  AAT(3,3)=  -0.720044076

!chiral
AAT(1,1)= -0.773471492617178
AAT(1,2)=  0.448186260089139
AAT(1,3)=  0.448186263037887

AAT(2,1)=  0.000000012288949
AAT(2,2)= -0.707106772908658
AAT(2,3)=  0.707106789464437

AAT(3,1)= -0.633831089572415
AAT(3,2)= -0.546926949394537
AAT(3,3)= -0.546926925573637
       
!achiral
!AAT(1,1)= -0.773471492617178
!AAT(1,2)= -0.448186260089139
!AAT(1,3)= -0.448186263037887

!AAT(2,1)= -0.000000012288949
!AAT(2,2)= -0.707106772908658
!AAT(2,3)=  0.707106789464437

!AAT(3,1)=  0.633831089572415
!AAT(3,2)= -0.546926949394537
!AAT(3,3)= -0.546926925573637
end



end

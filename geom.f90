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

function mls_gaussian(r,rcoef) result(phi)
  ! Gaussian weight function of MLS
  ! E.g. eq. (3.150) in Liu & Gu (2005)
  implicit none
  real :: r, phi
  real, intent(in) :: rcoef
     if (r.gt.1.0d0) then
       phi = 0.0d0
     else 
       phi = exp( - (r / rcoef )**2  )
     endif
  
  end

  function mls_gauss_deriv(r,rcoef) result(dphidr)
    ! Gaussian weight function of MLS
    ! Analytical derivative
    ! Given the weight function W(r) , this computes dW/dr
    ! For Cartesian coordinates, still need to convert from dW/dr to, e.g. , dW/dx with chain rule
    ! for use in MLS
    implicit none
    real :: r, dphidr
    real, intent(in) :: rcoef
       if (r.gt.1.0d0) then
        dphidr = 0.0d0
       else 
        dphidr = - 2.0d0 * r / rcoef**2 * exp( - (r / rcoef )**2  )
       endif
    
    end


function signDist_sphere(x0,x_cm,inp) result(phi)
implicit none
integer :: i,inp
real    :: mindist,t,d,d1
real    :: y1(3), y2(3,3), AA(3,3), AAT(3,3), AAT_P(3,3)
real, dimension(3) :: x0,x1,x2,x_cm
real    :: phi

! mindist = 1.e6


!    x1 = matmul(AAT_P, y1)
!    x1 =  matmul(AAT,x1) + x_cm

!   ! spherical cap
!     d = norm2(x1-x0)
!     mindist = min(d, mindist)

!   phi = 1 - mindist/(rad_p)

phi = sqrt( ( x0(1) - x_cm(1)  )**2 +  ( x0(2) - x_cm(2)  )**2 + ( x0(3) - x_cm(3)  )**2  ) - rad_p


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
 
 AAT(1,1) = 1.0d0
 AAT(1,2) = 0.0d0
 AAT(1,3) = 0.0d0

 AAT(2,1) = 0.0d0
 AAT(2,2) = 1.0d0
 AAT(2,3) = 0.0d0

 AAT(3,1) = 0.0d0
 AAT(3,2) = 0.0d0
 AAT(3,3) = 1.0d0
end

function signed_distance(x0,x1,nhat) result(phi)
  implicit none
  real, dimension(3) :: x0,x1,nhat
  real :: phi
  ! Compute the SIGNED distance of a point x0 to a PLANE defined by its normal nhat and a point on the plane, x1

  phi = dot_product( nhat , x0 - x1   )
end 



function heaviside(x) result(answer)
  implicit none
  real :: x, answer
  if (x .ge. 0) then
    answer = 1.0d0
  else
    answer = 0.0
  endif

end 

end
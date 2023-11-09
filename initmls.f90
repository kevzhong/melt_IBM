subroutine init_mlsForce
use param
use mls_param
implicit none
integer :: ntr
integer :: k
real :: I1, I2, I3

!
!     Quantites for mls 
!
  Hboxx = (1.0d0/dx1+1.0d0/dx2+1.0d0/dx3)/3.0d0 !Average Eulerian grid spacing

  celvol = 1.0d0 / (dx1*dx2*dx3) ! (uniform) Eulerian cell volume

  alph_dx = dx1/wscl
  !invdx1celvol = 1.0d0/dx1/celvol
  invwcon      = 1.0d0/wcon

  alph_dx_invwcon = alph_dx*invwcon

  do ntr = 1,maxnf

     !cfac(ntr) = sur(ntr,1)*invdx1celvol ! For isotropic grid only
    cfac(ntr) = ( sur(ntr,1) * Hboxx ) / celvol

  enddo

  ! particle moment of intertia
  I1 = 0.0261201133
  I2 = 0.0261201133
  I3 = 0.0261201133

  I_inv = 0.
  I_inv(1,1) = 1. / I1 
  I_inv(2,2) = 1. / I2 
  I_inv(3,3) = 1. / I3 

  I_inv2 = 0.
  I_inv2(1,1) = (I2-I3) / I1
  I_inv2(2,2) = (I3-I1) / I2
  I_inv2(3,3) = (I1-I2) / I3
end

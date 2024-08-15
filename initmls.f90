subroutine init_mlsForce
use param
use mls_param
implicit none
integer :: ntr
integer :: k, inp

!
! Quantites for mls 
!

  wcon = 0.3 !MLS Gaussian weight coefficient
  wscl = 1.5 ! Support-cage half-width

  h_eulerian = (1.0d0/dx1+1.0d0/dx2+1.0d0/dx3)/3.0d0 !Average Eulerian grid spacing
  celvol = 1.0d0 / (dx1*dx2*dx3) ! (uniform) Eulerian cell volume
  
  !A_eulerian =  ( celvol**(1.0/3.0) )**2 
  !A_thresh = PERC_Athresh * A_eulerian ! Threshold triangle area for remesh flagging

  E_thresh = PERC_Ethresh * h_eulerian
  
  V_thresh = V_ON_VE_PERC * celvol ! Threshold volume for object

  cfac(:,:) = ( sur(:,:) * h_eulerian ) / celvol

end

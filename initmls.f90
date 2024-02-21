subroutine init_mlsForce
use param
use mls_param
implicit none
integer :: ntr
integer :: k, inp

!
!     Quantites for mls 
!
  h_eulerian = (1.0d0/dx1+1.0d0/dx2+1.0d0/dx3)/3.0d0 !Average Eulerian grid spacing
  celvol = 1.0d0 / (dx1*dx2*dx3) ! (uniform) Eulerian cell volume
  
  A_eulerian =  ( celvol**(1.0/3.0) )**2 
  A_thresh = PERC_Athresh * A_eulerian ! Threshold triangle area for remesh flagging
  
  !KZ: For future multi-body melting, should cfac should be size cfac(maxnf,Nparticle)
  !do inp = 1, Nparticle
  !  do ntr = 1,maxnf
  !    !cfac(ntr) = sur(ntr,1)*invdx1celvol ! For isotropic grid only
  !    cfac(ntr,inp) = ( sur(ntr,inp) * h_eulerian ) / celvol
  !  enddo
  !enddo

  cfac(:,:) = ( sur(:,:) * h_eulerian ) / celvol

end

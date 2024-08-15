subroutine mlsForce
USE param 
USE mls_param 
USE mpih
use mpi_param
USE mls_local
implicit none
real,dimension(4,4) :: ptx
real,dimension(3) :: tsur_xyz,fsur_xyz,pos_vec, pos
real              :: angle
integer :: inp,ntr

do inp = 1, Nparticle
 do ntr = 1, maxnf
      if (isGhostFace(ntr,inp) .eqv. .false. ) then

      if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then

         call forc1(ntr,inp,ptxAB_q1(1:nel,ntr,inp),vel_tri(1,ntr,inp) )
         call forc2(ntr,inp,ptxAB_q2(1:nel,ntr,inp),vel_tri(2,ntr,inp) )
         call forctemp(ntr,inp,ptxAB_temp(1:nel,ntr,inp))
      endif

      !if(pind(6,ntr,inp).ge.kstart-1 .and. pind(6,ntr,inp).le.kend+1) then
      if(pind(6,ntr,inp).ge.kstart .and. pind(6,ntr,inp).le.kend) then
         call forc3(ntr,inp,ptxAB_q3(1:nel,ntr,inp),vel_tri(3,ntr,inp) )

      endif

   endif

 enddo
enddo


end subroutine mlsForce

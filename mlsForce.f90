subroutine mlsForce
USE param 
USE mls_param 
USE mpih
use mpi_param
USE mls_local
implicit none
real,dimension(4) :: ptx
real,dimension(3) :: tsur_xyz,fsur_xyz,pos_vec, pos
real              :: angle
integer :: inp,ntr
real :: force(3,Nparticle), torque(3,Nparticle)

force = 0.
torque = 0.

do inp = 1, Nparticle
 do ntr = 1, maxnf


      if(pind(3,ntr,inp).ge.kstart-1 .and. pind(3,ntr,inp).le.kend+1) then
         pos_vec(1:3) = tri_bar(1:3,ntr,inp) - pos_cm(1:3,inp)

         call forc1(ntr,inp,ptxAB_q1(1:nel,ntr,inp),vel_tri(1,ntr,inp),pos_vec,torque(1:3,inp),force(1,inp))
         call forc2(ntr,inp,ptxAB_q2(1:nel,ntr,inp),vel_tri(2,ntr,inp),pos_vec,torque(1:3,inp),force(2,inp))
         call forc3(ntr,inp,ptxAB_q3(1:nel,ntr,inp),vel_tri(3,ntr,inp),pos_vec,torque(1:3,inp),force(3,inp))

      endif

 enddo
enddo

do inp=1,Nparticle
   call mpi_globalsum_double_arr(force(1:3,inp),3)
   call mpi_globalsum_double_arr(torque(1:3,inp),3)
   fpxyz(:,inp) = fpxyz(:,inp) + force(:,inp)
   ftxyz(:,inp) = ftxyz(:,inp) + torque(:,inp)
   enddo
!   if(ismaster)then
!   write(*,*)fpxyz(:,1)
!   endif
  
end subroutine mlsForce

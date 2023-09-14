subroutine forc3(ntr,inp,ptxAB,Vel,rad,torque,force)
USE param
USE mls_param
USE local_arrays, only: vz
USE mpi_param, only: kstart, kend
USE mls_local, only: for_zc
use mpih, only: myid
implicit none
real,dimension(3) :: rad,torque
real,dimension(nel) :: ui
real,dimension(nel) :: ptxAB(nel)
integer :: inp,ntr,inw,i,j,k
real Um, Vel
integer, dimension(3) :: pind_i, pind_o
real force,for
integer :: ii,jj
integer :: kstartp

if (kstart.eq.1) then
 kstartp = 2
else
 kstartp = kstart
endif

if (pind(6,ntr,inp).ge.kstartp .and. pind(6,ntr,inp).le.kend) then

  pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1
  pind_i(2)=pind(2,ntr,inp)-1;  pind_o(2)=pind(2,ntr,inp)+1
  pind_i(3)=pind(6,ntr,inp)-1;  pind_o(3)=pind(6,ntr,inp)+1

  inw = 1
  Um  = 0.

  ! spline weights for all three components
  do k=pind_i(3), pind_o(3)
    do j=pind_i(2), pind_o(2)

      jj = modulo(j-1,n2m) + 1

       do i=pind_i(1), pind_o(1)

          ii = modulo(i-1,n1m) + 1
 
          Um = Um + vz(ii,jj,k)*ptxAB(inw)

          inw = inw + 1
    enddo
   enddo
  enddo

  ! =========
  ! forcing

  inw = 1

  do k=pind_i(3), pind_o(3)
   do j=pind_i(2), pind_o(2)

      jj = modulo(j-1,n2m) + 1

    do i=pind_i(1), pind_o(1)

         ii = modulo(i-1,n1m) + 1

         for_zc(ii,jj,k) = for_zc(ii,jj,k) +  cfac(ntr) * (Vel-Um) * ptxAB(inw)
         inw = inw +1
         enddo
         enddo
         enddo

  ! =========
  ! particle force
  for = (Vel-Um) * sur(ntr,inp) * invdx1dt
  force = force + for

  !-- torque
  torque(1) = torque(1) + for*rad(2)
  torque(2) = torque(2) - for*rad(1)


endif
end

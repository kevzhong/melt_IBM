subroutine forc2(ntr,inp,ptxAB,Vel)
USE param
USE mls_param
USE local_arrays, only: vy
USE mpi_param, only: kstart, kend
USE mls_local, only: for_yc
use mpih, only: myid
implicit none
real,dimension(nel) :: ui
real,dimension(nel) :: ptxAB(nel)
integer :: inp,ntr,inw,i,j,k,k1
real Um, Vel
integer, dimension(3) :: pind_i, pind_o
real for
integer :: ii,jj,kk

if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then

  pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1
  pind_i(2)=pind(5,ntr,inp)-1;  pind_o(2)=pind(5,ntr,inp)+1
  pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1


  Um = 0.
  inw = 1

  ! spline weights for all three components
  do k=pind_i(3), pind_o(3)
!    kk = modulo(k-1,n3m) + 1 !
     do j=pind_i(2), pind_o(2)

       jj = modulo(j-1,n2m) + 1

       do i=pind_i(1), pind_o(1)

          ii = modulo(i-1,n1m) + 1

          Um = Um + vy(ii,jj,k)*ptxAB(inw)

          inw = inw + 1
    enddo
   enddo
  enddo

  ! =========
  ! forcing

  inw = 1

  do k=pind_i(3), pind_o(3)
!    kk = modulo(k-1,n3m) + 1 !
   do j=pind_i(2), pind_o(2)

     jj = modulo(j-1,n2m) + 1

    do i=pind_i(1), pind_o(1)

         ii = modulo(i-1,n1m) + 1

         for_yc(ii,jj,k) = for_yc(ii,jj,k) + cfac(ntr,inp) * (Vel-Um) * ptxAB(inw)

         inw = inw +1
         enddo
         enddo
         enddo


endif
end

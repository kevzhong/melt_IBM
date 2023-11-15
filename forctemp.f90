subroutine forctemp(ntr,inp,ptxAB,Tmelt)
USE param
USE mls_param
USE local_arrays, only: temp
USE mpi_param, only: kstart, kend
USE mls_local, only: for_temp
use mpih, only: myid
implicit none
real,dimension(nel) :: ui
real,dimension(nel) :: ptxAB(nel)
integer :: inp,ntr,inw,i,j,k
real Tm, Tmelt
integer, dimension(3) :: pind_i, pind_o
integer :: ii,jj,kk
integer kstartp

if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then

    ! Support cage indices
  pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1 !i1 i2 i3
  pind_i(2)=pind(2,ntr,inp)-1;  pind_o(2)=pind(2,ntr,inp)+1 !j1 j2 j3
  pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1 !k1 k2 k3

  inw = 1
  Tm  = 0. !MLS-interpolated value at triangle centroid

  ! spline weights for all three components
  do k=pind_i(3), pind_o(3)
!kk = modulo(k-1,n3m) + 1 !
   do j=pind_i(2), pind_o(2)

     jj = modulo(j-1,n2m) + 1

     do i=pind_i(1), pind_o(1)

        ii = modulo(i-1,n1m) + 1

        Tm = Tm + temp(ii,jj,k)*ptxAB(inw)

        inw = inw + 1
    enddo
   enddo
  enddo

  ! =========
  ! forcing

  inw = 1

  do k=pind_i(3), pind_o(3)
!kk = modulo(k-1,n3m) + 1 !
   do j=pind_i(2), pind_o(2)

     jj = modulo(j-1,n2m) + 1

    do i=pind_i(1), pind_o(1)

         ii = modulo(i-1,n1m) + 1

         for_temp(ii,jj,k) = for_temp(ii,jj,k) + cfac(ntr) * (Tmelt-Tm) * ptxAB(inw)

         inw = inw +1
         enddo
         enddo
         enddo


!  ! =========
!  ! particle force
!  for = (Vel-Um) * sur(ntr,inp) * invdx1dt
!  force = force + for
!
!  !-- torque
  !torque(2) = torque(2) + for*rad(3)
  !torque(3) = torque(3) - for*rad(2)

endif
end


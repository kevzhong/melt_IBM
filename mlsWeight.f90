subroutine mlsWeight
USE param
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real,dimension(4) :: ptx
real,dimension(3) :: pos
integer :: inp,ntr

invdx1dt = 1.0d0/dx1/dt


do inp=1,Nparticle
 do ntr = 1, maxnf

      if(pind(3,ntr,inp).ge.kstart-1 .and. pind(3,ntr,inp).le.kend+1) then
    
         pos(1:3) = tri_bar(1:3,ntr,inp)

         ! initialise pre-factor matrix
         ptx(1)   = 1.d0; 
         ptx(2:4) = pos(1:3)


         call wght1(ntr,inp,pos,ptx,ptxAB_q1(1:nel,ntr,inp))
         call wght2(ntr,inp,pos,ptx,ptxAB_q2(1:nel,ntr,inp))
         call wght3(ntr,inp,pos,ptx,ptxAB_q3(1:nel,ntr,inp))

      endif

 enddo
enddo

end subroutine mlsWeight


subroutine wght1(ntr,inp,pos,ptx,ptxAB)
USE param
USE mls_param
USE mpi_param, only: kstart, kend
use geom
implicit none
real,dimension(4) :: pxk,ptx,ptxA
real,dimension(3) :: pos,norp,Wt
real,dimension(4,4) :: pinvA,invA
real,dimension(4,nel) :: B
real,dimension(nel) :: ptxAB(nel)
real :: Wtx, Wt23
integer :: inp,ntr,inw,i,j,k,k1
integer, dimension(3) :: pind_i, pind_o

if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then

  !-------------FORCING FUNCTION------------------------
  ! volume of a face with a specific marker - thickness taken as average of grid spacing

  pind_i(1)=pind(4,ntr,inp)-1;  pind_o(1)=pind(4,ntr,inp)+1
  pind_i(2)=pind(2,ntr,inp)-1;  pind_o(2)=pind(2,ntr,inp)+1
  !pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1

  k1  = floor(pos(3)*dx3) + 1
  pind_i(3)=k1-1
  pind_o(3)=k1+1


  inw = 1

  ! spline weights for all three components
  do k=pind_i(3), pind_o(3)
     norp(3) = abs(zm(k)-pos(3))*dx3
     do j=pind_i(2), pind_o(2)
          norp(2) = abs(ym(j)-pos(2))*dx2
         do i=pind_i(1), pind_o(1)
            norp(1) = abs(xc(i)-pos(1))*dx1

           ptxab(inw) =   delta(norp(1))*delta(norp(2))*delta(norp(3))

           inw = inw + 1
    enddo
   enddo
  enddo

endif
end

subroutine wght2(ntr,inp,pos,ptx,ptxAB)
USE param
USE mls_param
USE mpi_param, only: kstart, kend
use geom
implicit none
real,dimension(4) :: pxk,ptx,ptxA
real,dimension(3) :: pos,norp,Wt
real,dimension(4,4) :: pinvA,invA
real,dimension(4,nel) :: B
real,dimension(nel) :: ptxAB(nel)
real :: Wtx, Wt23
integer :: inp,ntr,inw,i,j,k,k1
integer, dimension(3) :: pind_i, pind_o


if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then

  pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1
  pind_i(2)=pind(5,ntr,inp)-1;  pind_o(2)=pind(5,ntr,inp)+1
  !pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1

  k1  = floor(pos(3)*dx3) + 1
  pind_i(3)=k1-1
  pind_o(3)=k1+1

  inw = 1

  ! spline weights for all three components
  do k=pind_i(3), pind_o(3)
     norp(3) = abs(zm(k)-pos(3))*dx3
     do j=pind_i(2), pind_o(2)
          norp(2) = abs(yc(j)-pos(2))*dx2
         do i=pind_i(1), pind_o(1)
            norp(1) = abs(xm(i)-pos(1))*dx1

           ptxab(inw) =   delta(norp(1))*delta(norp(2))*delta(norp(3))

           inw = inw + 1
    enddo
   enddo
  enddo



endif
end

subroutine wght3(ntr,inp,pos,ptx,ptxAB)
USE param
USE mls_param
USE mpi_param, only: kstart, kend
use geom
implicit none
real,dimension(4) :: pxk,ptx,ptxA
real,dimension(3) :: pos,norp,Wt
real,dimension(4,4) :: pinvA,invA
real,dimension(4,nel) :: B
real,dimension(nel) :: ptxAB(nel)
real :: Wtx,Wt23
integer :: inp,ntr,inw,i,j,k,kst
integer, dimension(3) :: pind_i, pind_o

if (pind(6,ntr,inp).ge.kstart .and. pind(6,ntr,inp).le.kend) then

  pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1
  pind_i(2)=pind(2,ntr,inp)-1;  pind_o(2)=pind(2,ntr,inp)+1
  !pind_i(3)=pind(6,ntr,inp)-1;  pind_o(3)=pind(6,ntr,inp)+1

  kst = nint(pos(3)*dx3)  + 1
  pind_i(3) = kst-1
  pind_o(3) = kst+1


  inw = 1

  ! spline weights for all three components
  do k=pind_i(3), pind_o(3)
     norp(3) = abs(zc(k)-pos(3))*dx3
     do j=pind_i(2), pind_o(2)
          norp(2) = abs(ym(j)-pos(2))*dx2
         do i=pind_i(1), pind_o(1)
            norp(1) = abs(xm(i)-pos(1))*dx1

           ptxab(inw) =   delta(norp(1))*delta(norp(2))*delta(norp(3))

           inw = inw + 1
    enddo
   enddo
  enddo

endif

end

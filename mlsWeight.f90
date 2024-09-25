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
  if (isGhostFace(ntr,inp) .eqv. .false. ) then
      if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then
    
         pos(1:3) = tri_bar(1:3,ntr,inp)
         ! initialise pre-factor matrix
         ptx(1)   = 1.d0; 
         ptx(2:4) = pos(1:3)

         call wght1(ntr,inp,pos,ptx,ptxAB_q1(1:nel,ntr,inp))
         call wght2(ntr,inp,pos,ptx,ptxAB_q2(1:nel,ntr,inp))
         !call wght3(ntr,inp,pos,ptx,ptxAB_q3(1:nel,ntr,inp))
         call wghttemp(ntr,inp,pos,ptx,ptxAB_temp(1:nel,ntr,inp))

      endif

      if(pind(6,ntr,inp).ge.kstart .and. pind(6,ntr,inp).le.kend) then

        pos(1:3) = tri_bar(1:3,ntr,inp)
        ! initialise pre-factor matrix
        ptx(1)   = 1.d0; 
        ptx(2:4) = pos(1:3)

        call wght3(ntr,inp,pos,ptx,ptxAB_q3(1:nel,ntr,inp))
      endif

    endif

 enddo
enddo

end subroutine mlsWeight


subroutine wght1(ntr,inp,pos,ptx,ptxAB)
USE param
USE geom
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real,dimension(4) :: pxk,ptx,ptxA
real,dimension(3) :: pos,norp,Wt
real,dimension(4,4) :: pinvA,invA
real,dimension(4,nel) :: B
real,dimension(nel) :: ptxAB(nel)
real :: Wtx, Wt23
integer :: inp,ntr,inw,i,j,k,k1
integer, dimension(3) :: pind_i, pind_o

! For LAPACK
integer :: LWORK = 4
real, dimension(4) :: WORK
integer :: IPIV(4), INFO

if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then

  !-------------FORCING FUNCTION------------------------
  ! volume of a face with a specific marker - thickness taken as average of grid spacing

  pind_i(1)=pind(4,ntr,inp)-1;  pind_o(1)=pind(4,ntr,inp)+1
  pind_i(2)=pind(2,ntr,inp)-1;  pind_o(2)=pind(2,ntr,inp)+1
! pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1

  k1  = floor(pos(3)*dx3) + 1
  pind_i(3)=k1-1
  pind_o(3)=k1+1



  pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
  inw = 1


! Accumulate A(4,4)   , B(4,27) linear system
  do k=pind_i(3), pind_o(3)
  
    norp(3)=abs(zm(k)-pos(3)) / (wscl / dx3)
    Wt(3) = mls_gaussian( norp(3) , wcon )
  
    do j=pind_i(2), pind_o(2)
  
        norp(2)=abs(ym(j)-pos(2)) / (wscl / dx2)
        Wt(2) = mls_gaussian( norp(2) , wcon )
        Wt23 = Wt(2)*Wt(3)
  
        do i=pind_i(1), pind_o(1)
  
            norp(1)=abs(xc(i)-pos(1)) / (wscl / dx1)
            Wt(1) = mls_gaussian( norp(1) , wcon )
            Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)
  
            pxk(1)=1.0d0
            pxk(2)=xc(i)
            pxk(3)=ym(j)
            pxk(4)=zm(k)
  
            call DGEMM('N','T',4,4,1,Wtx,pxk,4,pxk,4, 1.0d0,pinvA,4)
            B(1:4,inw)=Wtx*pxk(1:4)
  
            inw = inw + 1
        enddo !end i
    enddo !end j
enddo !end k

  ! calling routine to compute inverse
  !call inverseLU(pinvA,invA)
  call dgetrf(4, 4, pinvA, 4, IPIV, INFO) ! Compute the LU factorization of I_ij
  if (info /= 0) then
  print *, "Error in dgetrf"
  stop
  end if
  call dgetri(4, pinvA, 4, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization
  invA = pinvA
        
  !------------------------------------------------------
  ! matrix multiplications for final interpolation
  ! DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
  ! C = alpha * A * B + beta * C
  !---------------Shape function calculation---------------
  call DGEMM('N','N',1,4  ,4,1.0d0,ptx ,1,invA,4,0.0d0,ptxA ,1) 
  call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B   ,4,0.0d0,ptxAB,1) 



endif
end

subroutine wght2(ntr,inp,pos,ptx,ptxAB)
USE param
USE geom
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real,dimension(4) :: pxk,ptx,ptxA
real,dimension(3) :: pos,norp,Wt
real,dimension(4,4) :: pinvA,invA
real,dimension(4,nel) :: B
real,dimension(nel) :: ptxAB(nel)
real :: Wtx, Wt23
integer :: inp,ntr,inw,i,j,k,k1
integer, dimension(3) :: pind_i, pind_o

! For LAPACK
integer :: LWORK = 4
real, dimension(4) :: WORK
integer :: IPIV(4), INFO


if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then

  pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1
  pind_i(2)=pind(5,ntr,inp)-1;  pind_o(2)=pind(5,ntr,inp)+1
! pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1

  k1  = floor(pos(3)*dx3) + 1
  pind_i(3)=k1-1
  pind_o(3)=k1+1


  pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
  inw = 1

! Accumulate A(4,4)   , B(4,27) linear system
  do k=pind_i(3), pind_o(3)
  
    norp(3)=abs(zm(k)-pos(3)) / (wscl / dx3)
    Wt(3) = mls_gaussian( norp(3) , wcon )
  
    do j=pind_i(2), pind_o(2)
  
        norp(2)=abs(yc(j)-pos(2)) / (wscl / dx2)
        Wt(2) = mls_gaussian( norp(2) , wcon )
        Wt23 = Wt(2)*Wt(3)
  
        do i=pind_i(1), pind_o(1)
  
            norp(1)=abs(xm(i)-pos(1)) / (wscl / dx1)
            Wt(1) = mls_gaussian( norp(1) , wcon )
            Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)
  
            pxk(1)=1.0d0
            pxk(2)=xm(i)
            pxk(3)=yc(j)
            pxk(4)=zm(k)
  
            call DGEMM('N','T',4,4,1,Wtx,pxk,4,pxk,4, 1.0d0,pinvA,4)
            B(1:4,inw)=Wtx*pxk(1:4)
  
            inw = inw + 1
        enddo !end i
    enddo !end j
enddo !end k

  ! calling routine to compute inverse
  !call inverseLU(pinvA,invA)
  call dgetrf(4, 4, pinvA, 4, IPIV, INFO) ! Compute the LU factorization of I_ij
  if (info /= 0) then
  print *, "Error in dgetrf"
  stop
  end if
  call dgetri(4, pinvA, 4, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization
  invA = pinvA       

  !------------------------------------------------------
  ! matrix multiplications for final interpolation
  ! DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
  ! C = alpha * A * B + beta * C
  !---------------Shape function calculation---------------
  call DGEMM('N','N',1,4  ,4,1.0d0,ptx ,1,invA,4,0.0d0,ptxA ,1) 
  call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B   ,4,0.0d0,ptxAB,1) 

endif
end

subroutine wght3(ntr,inp,pos,ptx,ptxAB)
USE param
USE geom
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real,dimension(4) :: pxk,ptx,ptxA
real,dimension(3) :: pos,norp,Wt
real,dimension(4,4) :: pinvA,invA
real,dimension(4,nel) :: B
real,dimension(nel) :: ptxAB(nel)
real :: Wtx,Wt23
integer :: inp,ntr,inw,i,j,k,kst
integer, dimension(3) :: pind_i, pind_o

! For LAPACK
integer :: LWORK = 4
real, dimension(4) :: WORK
integer :: IPIV(4), INFO

if (pind(6,ntr,inp).ge.kstart .and. pind(6,ntr,inp).le.kend) then

  pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1
  pind_i(2)=pind(2,ntr,inp)-1;  pind_o(2)=pind(2,ntr,inp)+1
 !pind_i(3)=pind(6,ntr,inp)-1;  pind_o(3)=pind(6,ntr,inp)+1

  kst = nint(pos(3)*dx3)  + 1
  pind_i(3) = kst-1
  pind_o(3) = kst+1


  pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
  inw = 1


! Accumulate A(4,4)   , B(4,27) linear system
  do k=pind_i(3), pind_o(3)
  
    norp(3)=abs(zc(k)-pos(3)) / (wscl / dx3)
    Wt(3) = mls_gaussian( norp(3) , wcon )
  
    do j=pind_i(2), pind_o(2)
  
        norp(2)=abs(ym(j)-pos(2)) / (wscl / dx2)
        Wt(2) = mls_gaussian( norp(2) , wcon )
        Wt23 = Wt(2)*Wt(3)
  
        do i=pind_i(1), pind_o(1)
  
            norp(1)=abs(xm(i)-pos(1)) / (wscl / dx1)
            Wt(1) = mls_gaussian( norp(1) , wcon )
            Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)
  
            pxk(1)=1.0d0
            pxk(2)=xm(i)
            pxk(3)=ym(j)
            pxk(4)=zc(k)
  
            call DGEMM('N','T',4,4,1,Wtx,pxk,4,pxk,4, 1.0d0,pinvA,4)
            B(1:4,inw)=Wtx*pxk(1:4)
  
            inw = inw + 1
        enddo !end i
    enddo !end j
enddo !end k

  ! calling routine to compute inverse
  !call inverseLU(pinvA,invA)

  call dgetrf(4, 4, pinvA, 4, IPIV, INFO) ! Compute the LU factorization of I_ij
  if (info /= 0) then
  print *, "Error in dgetrf"
  stop
  end if
  call dgetri(4, pinvA, 4, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization
  invA = pinvA
        
  !------------------------------------------------------
  ! matrix multiplications for final interpolation
  ! DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
  ! C = alpha * A * B + beta * C
  !---------------Shape function calculation---------------
  call DGEMM('N','N',1,4  ,4,1.0d0,ptx ,1,invA,4,0.0d0,ptxA ,1) 
  call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B   ,4,0.0d0,ptxAB,1) 
endif
end

subroutine wghttemp(ntr,inp,pos,ptx,ptxAB)

USE param
USE geom
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real,dimension(4) :: pxk,ptx,ptxA
real,dimension(3) :: pos,norp,Wt
real,dimension(4,4) :: pinvA,invA
real,dimension(4,nel) :: B
real,dimension(nel) :: ptxAB(nel)
real :: Wtx, Wt23
integer :: inp,ntr,inw,i,j,k,k1
integer, dimension(3) :: pind_i, pind_o

! For LAPACK
integer :: LWORK = 4
real, dimension(4) :: WORK
integer :: IPIV(4), INFO

!-------------Shape function for cell centres (temp. or pressure cells) -------------------------

  
if(pind(3,ntr,inp).ge.kstart .and. pind(3,ntr,inp).le.kend) then
  
!-------------FORCING FUNCTION------------------------
! volume of a face with a specific marker - thickness taken as average of grid spacing
  
!WGHT1
pind_i(1)=pind(1,ntr,inp)-1;  pind_o(1)=pind(1,ntr,inp)+1 ! xi indices
pind_i(2)=pind(2,ntr,inp)-1;  pind_o(2)=pind(2,ntr,inp)+1 !yj indices
! pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1
  
! zk indinces
k1  = floor(pos(3)*dx3) + 1
pind_i(3)=k1-1
pind_o(3)=k1+1
  
pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
inw = 1
  
  
! Accumulate A(4,4)   , B(4,27) linear system
do k=pind_i(3), pind_o(3)
  
    norp(3)=abs(zm(k)-pos(3)) / (wscl / dx3)
    Wt(3) = mls_gaussian( norp(3) , wcon )
  
    do j=pind_i(2), pind_o(2)
  
        norp(2)=abs(ym(j)-pos(2)) / (wscl / dx2)
        Wt(2) = mls_gaussian( norp(2) , wcon )
        Wt23 = Wt(2)*Wt(3)
  
        do i=pind_i(1), pind_o(1)
  
            norp(1)=abs(xm(i)-pos(1)) / (wscl / dx1)
            Wt(1) = mls_gaussian( norp(1) , wcon )
            Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)
  
            pxk(1)=1.0d0
            pxk(2)=xm(i)
            pxk(3)=ym(j)
            pxk(4)=zm(k)
  
            call DGEMM('N','T',4,4,1,Wtx,pxk,4,pxk,4, 1.0d0,pinvA,4)
            B(1:4,inw)=Wtx*pxk(1:4)
  
            inw = inw + 1
        enddo !end i
    enddo !end j
enddo !end k
  
    ! calling routine to compute inverse
    ! SPD matrix for uniform grids, we can use Cholesky decomp. instead: dpotrf
    !call inverseLU(pinvA,invA)

    call dgetrf(4, 4, pinvA, 4, IPIV, INFO) ! Compute the LU factorization of I_ij
    if (info /= 0) then
    print *, "Error in dgetrf"
    stop
    end if
    call dgetri(4, pinvA, 4, IPIV, WORK, LWORK, INFO) ! Invert the LU factorization
    invA = pinvA
          
    !------------------------------------------------------
    ! matrix multiplications for final interpolation
    ! DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
    ! C = alpha * A * B + beta * C
    !---------------Shape function calculation---------------
    call DGEMM('N','N',1,4  ,4,1.0d0,ptx ,1,invA,4,0.0d0,ptxA ,1) 
    call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B   ,4,0.0d0,ptxAB,1) 
  
  
  
  endif
  end
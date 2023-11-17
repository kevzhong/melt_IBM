       subroutine Weight_W3(dif,dsk,w,ndex,const,idir)

       implicit none
       integer :: ndex,i,idir
       real, dimension (3,ndex) :: dif,dsk,dswk
       real, dimension (ndex,4) :: w
       real, dimension (3) :: dsw
       real :: ep,difx,dify,difz,drdx,drdy,drdz,const
       real :: dwxdx,dwydy,dwzdz,rx,ry,rz,wx,wy,wz,alpha

       alpha=0.3
       dswk(1:3,1:ndex)=const*dsk(1:3,1:ndex)

!----- From 1 to 4 of the second dimension of w: w,dwdx,dwdy,dwdz

       ep=1.0e-20

       do i=1,ndex

          difx=dif(1,i)
          dify=dif(2,i)
          difz=dif(3,i)

          if (abs(difx).le.ep) then
             drdx=0.0
          else
             drdx=(difx/abs(difx))/dswk(1,i)
          endif

          if (abs(dify).le.ep) then
             drdy=0.0
          else
             drdy=(dify/abs(dify))/dswk(2,i)
          endif

          if (abs(difz).le.ep) then
             drdz=0.0
          else
             drdz=(difz/abs(difz))/dswk(3,i)
          endif

          rx=abs(dif(1,i))/dswk(1,i)
          ry=abs(dif(2,i))/dswk(2,i)
          rz=abs(dif(3,i))/dswk(3,i)

          if (rx.gt.1.0) then
             wx=0.0
             dwxdx=0.0
          else
             wx=exp(-(rx/alpha)**2)
             dwxdx=(-1./(alpha**2)*2.*rx*exp(-(rx/alpha)**2))*drdx
          endif
 
          if (ry.gt.1.0) then
             wy=0.0
             dwydy=0.0
          else
             wy=exp(-(ry/alpha)**2)
             dwydy=(-1./(alpha**2)*2.*ry*exp(-(ry/alpha)**2))*drdy
          endif

          if (rz.gt.1.0) then
             wz=0.0
             dwzdz=0.0
          else
             wz=exp(-(rz/alpha)**2)
             dwzdz=(-1./(alpha**2)*2.*rz*exp(-(rz/alpha)**2))*drdz
          endif

          w(i,1)=wx*wy*wz
          w(i,2)=wy*wz*dwxdx
          w(i,3)=wx*wz*dwydy
          w(i,4)=wx*wy*dwzdz

       enddo

       return
       end subroutine Weight_W3

       subroutine Compute_Basis(gpos,gp)
!-------------------------------------------------------------------
!      Compute basis functions and their derivatives
!
!      P: (1,x,y,z)^T
!
!      The matrix GP(4x4) contains:
!
!      | p(4)    | 
!      | dpdx(4) |
!      | dpdy(4) |
!      | dpdz(4) |
!-------------------------------------------------------------------
       implicit none
       real, dimension (3) :: gpos
       real, dimension (4,4) :: gp

       gp(1,1)=1.0
       gp(1,2)=gpos(1)
       gp(1,3)=gpos(2)
       gp(1,4)=gpos(3)

       gp(2,1)=0.0
       gp(2,2)=1.0
       gp(2,3)=0.0
       gp(2,4)=0.0

       gp(3,1)=0.0
       gp(3,2)=0.0
       gp(3,3)=1.0
       gp(3,4)=0.0

       gp(4,1)=0.0
       gp(4,2)=0.0
       gp(4,3)=0.0
       gp(4,4)=1.0

       return
       end subroutine Compute_Basis
!...................................................................
!===================================================================
       subroutine Compute_AB(const,gpos,posk,dsk,a,b,ndex,idir)

       implicit none
       integer :: ndex,i,j,k,iii,ik,jk,ikk,jkk,kkk,idir
       real, dimension (3) :: gpos
       real, dimension (3,ndex) :: posk,dsk
       real, dimension (4,4,4) :: A
       real, dimension (4,ndex,4) :: B
       real, dimension (3,ndex) :: dif
       real, dimension (4,ndex) :: pk
       real, dimension (ndex,4) :: w
       real, dimension (4,4) :: pp
       real :: const

!----- From 1 to 4 of the third dimension of A: a,dadx,dady,dadz
!----- From 1 to 4 of the third dimension of B: b,dbdx,dbdy,dbdz

       do i=1,ndex
          pk(1,i)=1.0
          pk(2,i)=posk(1,i)
          pk(3,i)=posk(2,i)
          pk(4,i)=posk(3,i)
          dif(1,i)=gpos(1)-posk(1,i)
          dif(2,i)=gpos(2)-posk(2,i)
          dif(3,i)=gpos(3)-posk(3,i)
       enddo
 
       call Weight_W3(dif,dsk,w,ndex,const,idir) ! Exponential


!----- Compute B and its derivatives

       do i=1,4
          do j=1,ndex
             do k=1,4
                B(i,j,k)=pk(i,j)*w(j,k)
             enddo
          enddo
       enddo

!----- Compute A and its derivatives 

       A(1:4,1:4,1:4)=0.

       do iii=1,ndex
          do ik=1,4
             do jk=1,4
                pp(ik,jk)=pk(ik,iii)*pk(jk,iii)
             enddo
          enddo

          do ikk=1,4
             do jkk=1,4
                do kkk=1,4
                   A(ikk,jkk,kkk)=A(ikk,jkk,kkk)+w(iii,kkk)*pp(ikk,jkk) 
                enddo
             enddo
          enddo
       enddo

       return
       end subroutine Compute_AB
!...................................................................
!===================================================================
       subroutine MLS_ShapeFunc_3D(const,gpos,posk,dsk,phi,ndex,kwji,idir)

       implicit none
       integer :: ndex,i,j,k,kwji,idir
       real, dimension (3) :: gpos
       real, dimension (3,ndex) :: posk,dsk
       real, dimension (4,4) :: gp
       real, dimension (4,4) :: gam
       real, dimension (4,4,4) :: A
       real, dimension (4,ndex,4) :: B
       real, dimension (4) :: c
       real, dimension (4,4) :: aa
       real, dimension (4,ndex) :: phi
       real :: ep,const

       call Compute_Basis(gpos,gp)

       call Compute_AB(const,gpos,posk,dsk,a,b,ndex,idir)

       ep=1.0e-20

       c(1:4)=gp(1,1:4)

       aa(1:4,1:4)=a(1:4,1:4,1)

       gam(1:4,1:4)=0.0

!----- Compute gam 

       call GaussEqSolver_Sym(aa,c,ep,kwji)

       gam(1:4,1)=c(1:4)

!----- Compute dgamdx 

       do i=1,4
          c(i)=0.0
          do j=1,4
             c(i)=c(i)+A(i,j,2)*gam(j,1)
          enddo
       enddo

       do k=1,4
          c(k)=gp(2,k)-c(k)
       enddo

       do i=1,4
          do j=1,4
             aa(i,j)=a(i,j,1)
          enddo
       enddo

       call GaussEqSolver_Sym(aa,c,ep,kwji)

       gam(1:4,2)=c(1:4)

!----- Compute dgamdy 

       do i=1,4
          c(i)=0.0
          do j=1,4
             c(i)=c(i)+a(i,j,3)*gam(j,1)
          enddo
       enddo

       do k=1,4
          c(k)=gp(3,k)-c(k)
       enddo

       do i=1,4
          do j=1,4
             aa(i,j)=a(i,j,1)
          enddo
       enddo

       call GaussEqSolver_Sym(aa,c,ep,kwji)

       gam(1:4,3)=c(1:4)

!----- Compute dgamdz 

       do i=1,4
          c(i)=0.0
          do j=1,4
             c(i)=c(i)+a(i,j,4)*gam(j,1)
          enddo
       enddo

       do k=1,4
          c(k)=gp(4,k)-c(k)
       enddo

       do i=1,4
          do j=1,4
             aa(i,j)=a(i,j,1)
          enddo
       enddo

       call GaussEqSolver_Sym(aa,c,ep,kwji)

       gam(1:4,4)=c(1:4)

!----- Compute Phi and its derivatives

       do i=1,ndex
          phi(1:4,i)=0.0
          do j=1,4
             phi(1,i)=phi(1,i)+gam(j,1)*b(j,i,1)
             phi(2,i)=phi(2,i)+gam(j,2)*b(j,i,1)+gam(j,1)*b(j,i,2)
             phi(3,i)=phi(3,i)+gam(j,3)*b(j,i,1)+gam(j,1)*b(j,i,3)
             phi(4,i)=phi(4,i)+gam(j,4)*b(j,i,1)+gam(j,1)*b(j,i,4)
          enddo
       enddo

       return
       end subroutine MLS_ShapeFunc_3D
!...................................................................
!===================================================================
       subroutine GaussEqSolver_Sym(a,b,ep,kwji)

       implicit none
       integer :: i,j,k,io,jo,kwji,in,i1
       integer, dimension (5) :: m
       real, dimension (4,4) :: a
       real, dimension (4) :: b
       real :: ep,p,t

       kwji=0

       do 10 i=1,4
10     m(i)=i
       do 120 k=1,4
          p=0.0
          do 20 i=k,4
             do 20 j=k,4
                if (abs(a(i,j)).gt.abs(p)) then
                   p=a(i,j)
                   io=i
                   jo=j
                endif
20        continue
          if (abs(p)-ep) 30,30,35
30        kwji=1
          return
35        continue
          if (jo.eq.k) go to 45
          do 40 i=1,4
             t=a(i,jo)
             a(i,jo)=a(i,k)
             a(i,k)=t
40        continue
          j=m(k)
          m(k)=m(jo)
          m(jo)=j
45        if(io.eq.k) go to 55
          do 50 j=k,4
             t=a(io,j)
             a(io,j)=a(k,j)
             a(k,j)=t
50        continue
          t=b(io)
          b(io)=b(k)
          b(k)=t
55        p=1./p
          in=4-1
          if(k.eq.4) go to 65
          do 60 j=k,in
60        a(k,j+1)=a(k,j+1)*p
65        b(k)=b(k)*p
          if(k.eq.4) go to 120
          do 80 i=k,in
             do 70 j=k,in
70              a(i+1,j+1)=a(i+1,j+1)-a(i+1,k)*a(k,j+1)
80              b(i+1)=b(i+1)-a(i+1,k)*b(k)
120       continue
          do 130 i1=2,4
             i=4+1-i1
             do 130 j=i,in
130       b(i)=b(i)-a(i,j+1)*b(j+1)
          do 140 k=1,4
             i=m(k)
140       a(1,i)=b(k)
          do 150 k=1,4
150       b(k)=a(1,k)
          kwji=0

          return
          end subroutine GaussEqSolver_Sym
!...................................................................
 
       subroutine get_indices_27(idir,ist,jst,kst,ndex,nx,ny,nz,i,j,k,idp,            &
                                posk,dsk,x,y,z,m1,m2,m3,dx1,dx2,dx3,g2,g3,constp,     &
                                alx1,alx2,alx3)

       implicit none
       integer :: ist,jst,kst,m1,m2,m3,m,nn
       integer :: idir,ndex,nx,ny,nz,i,j,k,kk,jj,ii,n,nb
       integer, dimension (3,ndex) :: idp
       integer, dimension (9) :: mim1,mip1,mi
       real, dimension (3,ndex) :: posk,dsk
       integer imin,imax,jmin,jmax,kmin,kmax
       real, dimension(m1) :: x
       real, dimension(m2) :: y,g2
       real, dimension(m3) :: z,g3
       real constp,dx1,dx2,dx3,alx1,alx2,alx3

       mim1(1)=1
       mim1(2)=4
       mim1(3)=7
       mim1(4)=10
       mim1(5)=13
       mim1(6)=16
       mim1(7)=19
       mim1(8)=22
       mim1(9)=25

       mi(1)=2
       mi(2)=5
       mi(3)=8
       mi(4)=11
       mi(5)=14
       mi(6)=17
       mi(7)=20
       mi(8)=23
       mi(9)=26

       mip1(1)=3
       mip1(2)=6
       mip1(3)=9
       mip1(4)=12
       mip1(5)=15
       mip1(6)=18
       mip1(7)=21
       mip1(8)=24
       mip1(9)=27



       if (i.eq.1) then
          if (ist.eq.0) then
             do nn=1,9
                m=mim1(nn)
                idp(1,m)=nx-1 
                posk(1,m)=x(nx-1)-alx1
                dsk(1,m)=1.0/dx1
             enddo
          elseif (ist.eq.1) then
             do nn=1,9
                m=mim1(nn)
                idp(1,m)=nx 
                posk(1,m)=x(nx)-alx1
                dsk(1,m)=1.0/dx1
             enddo
          endif  
          do nn=1,9
             m=mi(nn)
             idp(1,m)=i 
             posk(1,m)=x(i)
             dsk(1,m)=1.0/dx1
          enddo
          do nn=1,9
             m=mip1(nn)
             idp(1,m)=i+1 
             posk(1,m)=x(i+1)
             dsk(1,m)=1.0/dx1
          enddo
       endif 

       if (i.eq.nx) then
          if (ist.eq.0) then
             do nn=1,9
                m=mip1(nn)
                idp(1,m)=1+1 
                posk(1,m)=x(1+1)+alx1
                dsk(1,m)=1.0/dx1
             enddo
          elseif (ist.eq.1) then
             do nn=1,9
                m=mip1(nn)
                idp(1,m)=1 
                posk(1,m)=x(1)+alx1
                dsk(1,m)=1.0/dx1
             enddo
          endif   
          do nn=1,9
             m=mi(nn)
             idp(1,m)=i         
             posk(1,m)=x(i)
             dsk(1,m)=1.0/dx1
          enddo
          do nn=1,9
             m=mim1(nn)
             idp(1,m)=i-1
             posk(1,m)=x(i-1)
             dsk(1,m)=1.0/dx1
          enddo
       endif 

       if (i.ne.1.and.i.ne.nx) then
          m=0
          do kk=-1,1
             do jj=-1,1
                do ii=-1,1
                   m=m+1
                   idp(1,m)=i+ii 
                   posk(1,m)=x(i+ii)
                   dsk(1,m)=1.0/dx1
                enddo       
             enddo      
          enddo        
       endif 


       
       m=0
       do kk=-1,1
          do jj=-1,1
             do ii=-1,1
                m=m+1

                idp(2,m)=j+jj 
                posk(2,m)=y(j+jj)
                dsk(2,m)=g2(j+jj)/dx2
 
                idp(3,m)=k+kk 
                posk(3,m)=z(k+kk)
                dsk(3,m)=g3(k+kk)/dx3

             enddo       
          enddo      
       enddo        


       return
       end subroutine get_indices_27
!.......................................................................

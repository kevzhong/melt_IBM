!------------------------------------------------------------------
!     Routine to compute indices of all centroids on all particles
!     Only starting and ending indices of the support domain are 
!     reqiured. Then loop over it.
!------------------------------------------------------------------
subroutine findCentroidIndices
  USE mpih
  USE param
  USE mls_param, only : pind,Nparticle,maxnf,tri_bar,isGhostFace
  IMPLICIT NONE
  real pos(3)
  integer i1,j1,k1,ist,jst,kst,inp,ntr

  do inp=1,Nparticle
    do ntr=1,maxnf
      if (isGhostFace(ntr,inp) .eqv. .false. ) then
       pos(1:3) = tri_bar(1:3,ntr,inp)

       !++++++++Indices of the marker+++++++++++++++
       !     X - indices
       i1  = floor(pos(1)*dx1) + 1
       ist = nint(pos(1)*dx1)  + 1

       i1  = modulo(i1-1,n1m)  + 1
       ist = modulo(ist-1,n1m) + 1

       !     Y - indices
       j1  = floor(pos(2)*dx2) + 1
       jst = nint(pos(2)*dx2)  + 1

       j1  = modulo(j1-1,n2m)  + 1
       jst = modulo(jst-1,n2m) + 1

       !     Z - indices
       k1  = floor(pos(3)*dx3) + 1
       kst = nint(pos(3)*dx3)  + 1

       k1  = modulo(k1-1,n3m)  + 1
       kst = modulo(kst-1,n3m) + 1

       !-------------------------------------------------------------
       pind(1,ntr,inp)=i1 ; pind(2,ntr,inp)=j1 ; pind(3,ntr,inp)=k1
       pind(4,ntr,inp)=ist; pind(5,ntr,inp)=jst; pind(6,ntr,inp)=kst
       !-------------------------------------------------------------
      endif
    end do
  end do

end subroutine findCentroidIndices

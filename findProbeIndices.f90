!------------------------------------------------------------------
!     Routine to compute the indices of all positive and negative probes
!     extruded from the triangle centroid locations of each particle
!     extrusion is done in the +ve/-ve normal direction, with a length
!     equal to the average Euerlian grid spacing

!     Only required for cell-centered temperature field
!------------------------------------------------------------------
subroutine findProbeIndices
  USE mpih
  USE param
  USE mls_param, only : pind,Nparticle,maxnf,tri_bar, tri_nor, h_eulerian, isGhostFace
  IMPLICIT NONE
  real pos(3)
  integer i1,j1,k1,inp,nf
  real :: s

  ! Re-use pind() memory from centroid indices

  do inp=1,Nparticle
    do nf=1,maxnf
      if ( isGhostFace(nf,inp) .eqv. .false. ) then
        ! Get nhat vector length: should be unity but compute anyway for safety
        s = norm2( tri_nor(1:3,nf,inp) )

       !------------------ POSITIVE PROBE-------------------------
       pos(1:3) = tri_bar(1:3,nf,inp) + h_eulerian * tri_nor(1:3,nf,inp) / s !Positive probe location

       !++++++++Indices of the marker+++++++++++++++
       !     X - indices
       i1  = floor(pos(1)*dx1) + 1
       !ist = nint(pos(1)*dx1)  + 1
       !i1  = modulo(i1-1,n1m)  + 1

       !     Y - indices
       j1  = floor(pos(2)*dx2) + 1
       !jst = nint(pos(2)*dx2)  + 1
       !j1  = modulo(j1-1,n2m)  + 1

       !     Z - indices
       k1  = floor(pos(3)*dx3) + 1
       !kst = nint(pos(3)*dx3)  + 1

       k1  = modulo(k1-1,n3m)  + 1
       !kst = modulo(kst-1,n3m) + 1

       pind(1,nf,inp)=i1 ; pind(2,nf,inp)=j1 ; pind(3,nf,inp)=k1   ! indices for positive probe

       !------------------ NEGATIVE PROBE-------------------------
       pos(1:3) = tri_bar(1:3,nf,inp) - h_eulerian * tri_nor(1:3,nf,inp) / s ! negative probe location

       !++++++++Indices of the marker+++++++++++++++
       !     X - indices
       i1  = floor(pos(1)*dx1) + 1
       !ist = nint(pos(1)*dx1)  + 1
       !i1  = modulo(i1-1,n1m)  + 1

       !     Y - indices
       j1  = floor(pos(2)*dx2) + 1
       !jst = nint(pos(2)*dx2)  + 1
       !j1  = modulo(j1-1,n2m)  + 1

       !     Z - indices
       k1  = floor(pos(3)*dx3) + 1
       !kst = nint(pos(3)*dx3)  + 1

       k1  = modulo(k1-1,n3m)  + 1
       !kst = modulo(kst-1,n3m) + 1

       pind(4,nf,inp)=i1; pind(5,nf,inp)=j1; pind(6,nf,inp)=k1  ! indices for negative probe
       !------------------END NEGATIVE PROBE-------------------------
      endif
    end do
  end do

end subroutine findProbeIndices

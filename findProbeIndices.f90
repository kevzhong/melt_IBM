!------------------------------------------------------------------
!     Routine to compute the indices of all positive and negative probes
!     extruded from the vertex locations of each particle
!     extrusion is done in the +ve/-ve normal direction, with a length
!     equal to the average Euerlian grid spacing

!     Only required for cell-centered temperature field
!------------------------------------------------------------------
subroutine findProbeIndices
  USE mpih
  USE param
  USE mls_param, only : pindv,Nparticle,maxnv,xyzv, vert_nor, h_eulerian
  IMPLICIT NONE
  real pos(3)
  integer i1,j1,k1,,inp,nv
  real :: s

  do inp=1,Nparticle
    do nv=1,maxnv
        ! Get nhat vector length: should be unity but compute anyway for safety
        s = norm2( vert_nor(1:3,nv,inp) )

       !------------------ POSITIVE PROBE-------------------------
       pos(1:3) = xyzv(1:3,nv,inp) + h_eulerian * vert_nor(1:3,nv,inp) / s !Positive probe location

       !++++++++Indices of the marker+++++++++++++++
       !     X - indices
       i1  = floor(pos(1)*dx1) + 1
       !ist = nint(pos(1)*dx1)  + 1

       !     Y - indices
       j1  = floor(pos(2)*dx2) + 1
       !jst = nint(pos(2)*dx2)  + 1

       !     Z - indices
       k1  = floor(pos(3)*dx3) + 1
       !kst = nint(pos(3)*dx3)  + 1

       k1  = modulo(k1-1,n3m)  + 1
       !kst = modulo(kst-1,n3m) + 1

       pindv(1,nv,inp)=i1 ; pindv(2,nv,inp)=j1 ; pindv(3,nv,inp)=k1   ! indices for positive probe

       !------------------ NEGATIVE PROBE-------------------------
       pos(1:3) = xyzv(1:3,nv,inp) - h_eulerian * vert_nor(1:3,nv,inp) / s ! negative probe location

       !++++++++Indices of the marker+++++++++++++++
       !     X - indices
       i1  = floor(pos(1)*dx1) + 1
       !ist = nint(pos(1)*dx1)  + 1

       !     Y - indices
       j1  = floor(pos(2)*dx2) + 1
       !jst = nint(pos(2)*dx2)  + 1

       !     Z - indices
       k1  = floor(pos(3)*dx3) + 1
       !kst = nint(pos(3)*dx3)  + 1

       k1  = modulo(k1-1,n3m)  + 1
       !kst = modulo(kst-1,n3m) + 1

       pind(4,nv,inp)=i1; pind(5,nv,inp)=j1; pind(6,nv,inp)=k1  ! indices for negative probe
       !------------------END NEGATIVE PROBE-------------------------
    end do
  end do

end subroutine findProbeIndices

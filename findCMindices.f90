
subroutine findCMindices
  USE mpih
  USE param
  USE mls_param
  IMPLICIT NONE
  integer i2,j2,k2,inp

  do inp=1,Nparticle
  !++++++++Indices of the CM+++++++++++++++
  !     X - indices
  i2  = floor(pos_cm(1,inp)*dx1) + 1
  i2   = modulo(i2-1,n1m) + 1

  !     Y - indices
  j2  = floor(pos_cm(2,inp)*dx2) + 1
  j2 = modulo(j2-1,n2m) + 1

  !     Z - indices
  k2  = floor(pos_cm(3,inp)*dx3) + 1
  k2  = modulo(k2-1,n3m)  + 1

  !-------------------------------------------------------------
  pind1(1,inp)=i2 ; pind1(2,inp)=j2 ; pind1(3,inp)=k2
  !-------------------------------------------------------------

enddo

end subroutine findCMindices

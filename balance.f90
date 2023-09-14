subroutine balance
 USE param
 USE mls_param
 use mls_local
 use mpi_param
 use local_arrays

 implicit none
 integer ic,jc,kc,ip,jp,kp,inp
 real :: fpcoupling
 character(70) namfile
 fpcoupling=0.0d0
     do inp=1,Nparticle

     fpcoupling=fpxyz(1,inp)*vel_cm(1,inp)+ &
                fpxyz(2,inp)*vel_cm(2,inp)+ &
                fpxyz(3,inp)*vel_cm(3,inp)

     end do
 

      if(ismaster) then
      namfile='flowmov/fpcoupling.txt'
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time,fpcoupling
      close(92)
      end if
   

end subroutine balance

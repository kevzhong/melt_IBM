SUbroutine particle
  use param
  use mls_param
  use mls_local
  use local_arrays
  use mpi_param
  use mpih

  implicit none 

  integer :: mstep, inp
  integer :: my_up,my_down

  if(imlsfor.eq.1)then

    call findindices
    call mlsWeight

    fpxyz=0.0d0
    ftxyz=0.0d0

     my_down=myid-1
     my_up=myid+1

    do mstep=1,1
      call update_both_ghosts(n1,n2,vx,kstart,kend)
      call update_both_ghosts(n1,n2,vy,kstart,kend)
      call update_both_ghosts(n1,n2,vz,kstart,kend)
    

      for_xc = 0.
      for_yc = 0.
      for_zc = 0.

      call mlsForce
      call velforce
      end do


  endif
 
    call update_both_ghosts(n1,n2,vx,kstart,kend)
    call update_both_ghosts(n1,n2,vy,kstart,kend)
    call update_both_ghosts(n1,n2,vz,kstart,kend)



if (imlsstr.eq.1) then
      call update_part_pos
  endif
 end

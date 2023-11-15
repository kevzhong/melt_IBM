subroutine tempforce
 USE param
 USE mls_param
 use mls_local
 use mpi_param
 use local_arrays

 implicit none
 integer ic,jc,kc,ip,jp,kp
 real :: for_temp_mean
 character(70) namfile

 call update_add_upper_ghost(for_temp)
 call update_add_lower_ghost(for_temp)
 
 for_temp_mean=0.

 do kc=kstart, kend
  do jc=1,n2m
   do ic=1,n1m
      for_temp_mean = for_temp_mean + for_temp(ic,jc,kc)
   end do
  end do
 end do

call mpi_globalsum_double_var(for_temp_mean)

for_temp_mean = for_temp_mean / dble(n1m*n2m*n3m)

 do kc=kstart,kend
  do jc=1,n2m
   do ic=1,n1m
      temp(ic,jc,kc) = temp(ic,jc,kc) +  for_temp(ic,jc,kc) - for_temp_mean
   end do
  end do
 end do

end

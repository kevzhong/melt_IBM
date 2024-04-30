subroutine velforce
 USE param
 USE mls_param
 use mls_local
 use mpi_param
 use local_arrays

 implicit none
 integer ic,jc,kc,ip,jp,kp
 real :: for_x_mean, for_y_mean, for_z_mean
 character(70) namfile

 call update_add_upper_ghost(for_xc)
 call update_add_upper_ghost(for_yc)
 call update_add_upper_ghost(for_zc)

 call update_add_lower_ghost(for_xc)
 call update_add_lower_ghost(for_yc)
 call update_add_lower_ghost(for_zc)
 
 for_x_mean=0.
 for_y_mean=0.
 for_z_mean=0.

 do kc=kstart, kend
  do jc=1,n2m
   do ic=1,n1m
      for_x_mean = for_x_mean + for_xc(ic,jc,kc)
      for_y_mean = for_y_mean + for_yc(ic,jc,kc)
      for_z_mean = for_z_mean + for_zc(ic,jc,kc)
   end do
  end do
 end do

call mpi_globalsum_double_var(for_x_mean)
call mpi_globalsum_double_var(for_y_mean)
call mpi_globalsum_double_var(for_z_mean)

for_x_mean = for_x_mean / dble(n1m*n2m*n3m)
for_y_mean = for_y_mean / dble(n1m*n2m*n3m)
for_z_mean = for_z_mean / dble(n1m*n2m*n3m)

 do kc=kstart,kend
  do jc=1,n2m
   do ic=1,n1m
      vx(ic,jc,kc) = vx(ic,jc,kc) +  for_xc(ic,jc,kc) !- for_x_mean
      vy(ic,jc,kc) = vy(ic,jc,kc) +  for_yc(ic,jc,kc) !- for_y_mean
      vz(ic,jc,kc) = vz(ic,jc,kc) +  for_zc(ic,jc,kc) !- for_z_mean
   end do
  end do
 end do

 ! To obtain a statistically steady state, mean must be subtracted
 !if (forcing .eq. 1) then 
   vx(:,:,:) = vx(:,:,:) - for_x_mean
   vy(:,:,:) = vy(:,:,:) - for_y_mean
   vz(:,:,:) = vz(:,:,:) - for_z_mean
 !endif

end

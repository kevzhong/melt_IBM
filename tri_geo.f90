subroutine tri_geo
      
  use param
  use mls_param
  use coll_mod
  implicit none

  integer :: inp

  ! Set maxnv, maxne and maxnf
  call set_particle_array_sizes

  ! max n vert
  call get_maxn_n_edge_of_vert
  ! Allocations using vertices,edges and faces as parameters     
  call allocate_trigeo

  ! Read in the vertices, edges and faces from the file
  ! Also make connections between the three
  call set_connectivity

  if(pread.eq.0) then
    ! Reposition particles and set physical quantities
    call setup_particles

  else if(pread.eq.1)then
    ! read in particle data from file
    call import_particles
    call import_collision
  end if

  call set_particle_rad !KZ: hard-coded AAT matrix influences tri_bar (centroid) calculations
  call set_xyz

  do inp=1,Nparticle

        call calculate_distance(dist(:,inp),maxnv,maxne,xyz0(:,:), vert_of_edge(:,:))
        call calculate_area(Surface(inp),maxnv,maxnf,xyz0(:,:), vert_of_face(:,:),sur(:,inp))
        call calculate_normal(tri_nor(:,:,inp),maxnv,maxnf,xyz0(:,:), vert_of_face(:,:))
        call calculate_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,VERTBUFFER,maxnf,faces_of_vert)
        call calculate_volume(Volume(inp),maxnv,maxnf,xyz0(:,:),vert_of_face(:,:),vol(:,inp))
        !call calculate_volume2 (Volume(inp),maxnf,tri_nor(:,:,inp),sur(:,inp),tri_bar(:,:,inp))
  end do

       
        if(ismaster)then
        write(*,*)'Volume fraction: ',Volume(1)*Nparticle*1.0e2/(xlen*ylen*zlen)
        end if

  call print_particle_info

  call init_mlsForce

end subroutine tri_geo

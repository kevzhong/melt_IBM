subroutine tri_geo
      
  use param
  use mls_param
  use coll_mod
  implicit none

  integer :: inp

  ! Set maxnv, maxne and maxnf
  call set_particle_array_sizes

  ! max n vert
  ! call get_maxn_n_edge_of_vert
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
  call set_xyz ! KZ: Setting COM-relative xyz

  do inp=1,Nparticle

        call calculate_eLengths(eLengths(:,inp),maxnv,maxne,xyz0(:,:), vert_of_edge(:,:),isGhostEdge(:,inp))

        call calculate_area(Surface(inp),maxnv,maxnf,xyz0(:,:), vert_of_face(:,:),sur(:,inp),&
                            isGhostFace(:,inp),rm_flag(inp),A_thresh)

        call calculate_vert_area (Avert(:,inp),maxnv,maxnf,vert_of_face(:,:),sur(:,inp),isGhostFace(:,inp))

        call calculate_normal(tri_nor(:,:,inp),maxnv,maxnf,xyz0(:,:), vert_of_face(:,:))

        call calculate_areaWeighted_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,maxnf,sur(:,inp),&
        vert_of_face(:,:),isGhostFace(:,inp),isGhostVert(:,inp))

        !call calculate_volume(Volume(inp),maxnv,maxnf,xyz0(:,:),vert_of_face(:,:),vol(:,inp))

        call calculate_volume2 (Volume(inp),maxnf,tri_nor(:,:,inp),sur(:,inp),tri_bar(:,:,inp),isGhostFace(:,inp))

  end do

       
        if(ismaster)then
        write(*,*)'Volume fraction: ',Volume(1)*Nparticle*1.0e2/(xlen*ylen*zlen)
        write(*,*)'Volume of object: ',Volume(1)

        end if



  call print_particle_info

  call init_mlsForce

end subroutine tri_geo

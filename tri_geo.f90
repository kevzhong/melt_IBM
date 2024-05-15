subroutine tri_geo
      
  use param
  use mls_param
  !use coll_mod KZ: no collisions for now
  implicit none

  integer :: inp

  ! Set maxnv, maxne and maxnf
  call set_particle_array_sizes

  ! max n vert
  ! call get_maxn_n_edge_of_vert
  ! Allocations using vertices,edges and faces as parameters     
  call allocate_trigeo


  if(pread.eq.0) then ! Read in initial geometry
    call set_connectivity ! Read xyz vertex coordinates and triangle connectivity arrays

    ! Reposition particles and set physical quantities
    call setup_particles

    call init_geomCoords !KZ: Set vertex and centroid coordinates in xyz space
    !call set_xyz 
  
    do inp=1,Nparticle
      call calculate_eLengths(eLengths(:,inp),maxnv,maxne,xyz0(:,:), vert_of_edge(:,:,inp),isGhostEdge(:,inp))
      call calculate_area(Surface(inp),maxnv,maxnf,xyz0(:,:), vert_of_face(:,:,inp),sur(:,inp),&
                        isGhostFace(:,inp),rm_flag(inp),A_thresh)
      call calculate_vert_area (Avert(:,inp),maxnv,maxnf,vert_of_face(:,:,inp),sur(:,inp),isGhostFace(:,inp))
      call calculate_skewness (maxne,maxnf,edge_of_face(:,:,inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),isGhostFace(:,inp),&
                              rm_flag(inp),skew_thresh)
      call calculate_normal(tri_nor(:,:,inp),maxnv,maxnf,xyz0(:,:), vert_of_face(:,:,inp))
      call calculate_areaWeighted_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,maxnf,sur(:,inp),&
            vert_of_face(:,:,inp),isGhostFace(:,inp),isGhostVert(:,inp))

      ! Volume is pre-computed by rigidBody_calcs
      !call calculate_volume2 (Volume(inp),maxnf,tri_nor(:,:,inp),sur(:,inp),tri_bar(:,:,inp),isGhostFace(:,inp))

    enddo

  else if(pread.eq.1)then ! Otherwise, read from continuation files
    ! read in particle data from file
    call import_particles
    
    !if (Nparticle .gt. 1 )then
    !  call import_collision
    !endif
  end if

        if(ismaster)then
        write(*,*)'Volume fraction: ',sum( Volume )*1.0e2/(xlen*ylen*zlen)
        !write(*,*)'Volume of object: ',Volume(1)
        end if

  call print_particle_info

  call init_mlsForce

end subroutine tri_geo

Subroutine particle
use param
use mls_param
use mls_local
use local_arrays
use mpi_param
use mpih

implicit none 

integer :: mstep, inp, i,ntr
integer :: my_up,my_down

if(imlsfor.eq.1)then

    call findCentroidIndices
    call mlsWeight

    fpxyz=0.0d0
    ftxyz=0.0d0

    my_down=myid-1
    my_up=myid+1

    do mstep=1,1 !Multi-direct forcing iteration. cf. Breugem (2012) eqn. 10
        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)
        call update_both_ghosts(n1,n2,temp,kstart,kend)


        for_xc = 0.0d0
        for_yc = 0.0d0
        for_zc = 0.0d0
        for_temp = 0.0d0

        call mlsForce

        if (imelt .eq. 1) then 
            call findProbeIndices ! Indices of inward/outward probes extrapolated from triangle faces
            call mls_heatFlux ! Calculate heat flux at +/- faces, then interpolate to vertices
            
            call MPI_ALLREDUCE(MPI_IN_PLACE,qw_oVert,maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
            call MPI_ALLREDUCE(MPI_IN_PLACE,qw_iVert,maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    

        endif

        call velforce
        call tempforce
    end do
 
call update_both_ghosts(n1,n2,vx,kstart,kend)
call update_both_ghosts(n1,n2,vy,kstart,kend)
call update_both_ghosts(n1,n2,vz,kstart,kend)
call update_both_ghosts(n1,n2,temp,kstart,kend)


if (imelt .eq. 1) then    
    call apply_StefanCondition ! Update vertex locations xyzv
    
    ! Update triangulated geometry details
    ! KZ: Entirety of same computation done by each process since each process stores all the geo info, could be parallelised later if a bottleneck
    do inp = 1,Nparticle
        call calculate_area(Surface(inp),maxnv,maxnf,xyzv(1:3,:,inp),vert_of_face(:,:,inp),sur(:,inp),&
                            isGhostFace(:,inp),rm_flag(inp),A_thresh) ! Update sur
        call calculate_eLengths(eLengths(:,inp),maxnv,maxne,xyzv(1:3,:,inp), vert_of_edge(:,:,inp),isGhostEdge(:,inp))
        call update_tri_normal (tri_nor(:,:,inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
        call calculate_skewness (maxne,maxnf,edge_of_face(:,:,inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),isGhostFace(:,inp)&
                                ,rm_flag(inp),skew_thresh)

        ! Remesh if <= threshold_area detected from calculate_area()
        if ( (rm_flag(inp) .eqv. .true.) .and. (iremesh .eq. 1 ) ) then
            call main_remesh (Surface(inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),maxnf,maxne,maxnv,&
                        xyzv(:,:,inp),tri_nor(:,:,inp),&
                        A_thresh,skew_thresh,vert_of_face(:,:,inp),edge_of_face(:,:,inp),vert_of_edge(:,:,inp),&
                        face_of_edge(:,:,inp),isGhostFace(:,inp),isGhostEdge(:,inp),isGhostVert(:,inp),rm_flag(inp))
        endif

        call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) ! Update tri_bar
        call calculate_vert_area (Avert(:,inp),maxnv,maxnf,vert_of_face(:,:,inp),sur(:,inp),isGhostFace(:,inp)) ! Update vertex areas
        call calculate_volume2 (Volume(inp),maxnf,tri_nor(:,:,inp),sur(:,inp),tri_bar(:,:,inp),isGhostFace(:,inp))

        call calculate_areaWeighted_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,maxnf,sur(:,inp),vert_of_face(:,:,inp),&
                                        isGhostFace(:,inp), isGhostVert(:,inp) )
    enddo

    ! Update Eulerian < -- > Lagrangian forcing transfer coefficient
    cfac(:,:) = ( sur(:,:) * h_eulerian ) / celvol ! Note the hard-coded single-particle for cfac

endif

if (imlsstr.eq.1) then
    call update_part_pos
endif

endif

 end

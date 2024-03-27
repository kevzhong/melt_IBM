!------------------------------------------------------
subroutine main_remesh (Surface,sur,eLengths,skewness,nf,ne,nv,xyz,tri_nor,A_thresh,skew_thresh,&
                        vert_of_face,edge_of_face,vert_of_edge,face_of_edge,&
                        isGhostFace,isGhostEdge,isGhostVert,rm_flag,anchorVert,flagged_edge)

    use mpih
    use param, only: ismaster

    implicit none
    integer :: nf,nv,ne
    real :: A_thresh, Surface, skew_thresh
    integer, dimension (3,nf) :: vert_of_face, edge_of_face
    integer, dimension(2,ne) :: vert_of_edge, face_of_edge
    logical, dimension(nf) :: isGhostFace
    logical, dimension(ne) :: isGhostEdge, flagged_edge
    logical, dimension(nv) :: isGhostVert, anchorVert
    logical :: rm_flag
    real, dimension(3,nv) :: xyz
    real, dimension(3,nf) :: tri_nor
    real, dimension(nf) :: sur, skewness
    real, dimension(ne) :: eLengths
    real, dimension(4,4) :: Q1, Q2 ! Error quadrics of vertices v1, v2
    integer :: f, e, v1, v2
    integer :: F1, F2, e1, e2, re1, re2
    integer :: ecol_cnt, EN, flip_cnt, cnt, niter
    character*50 :: dsetname,filename

    !integer :: valence, i

    ! Iteratively remesh a triangulated geometry using quadric-error-metrics-guided edge collapses
    ! QEM algorithm from "Surface Simplification Using Quadric Error Metrics" (Garland & Heckert, 1997)
    ! Here, we perform edge collapses based on triangles being smaller than a threshold Eulerian mesh measure,
    ! rather than a greedy cost sorting as in the original Garland & Heckbert implementation

    ! After each edge-collapse operation we remove:
    ! -3 edges
    ! -1 vertex
    ! -2 faces

    ! For topology-preserving operations (which an edge-collapse is), Euler's polyhedral number, EN, should stay the same:
    ! EN = numVert - numEdge + numFace
    ! where we will have EN = 2 for closed-manifold surfaces
    ecol_cnt = 0
    flip_cnt = 0

    do while (rm_flag .eqv. .true.) 
    do f = 1,nf
        if ( (isGhostFace(f) .eqv. .false.) .and.  ( ( sur(f) .le. A_thresh) .or. (skewness(f) .ge. skew_thresh) ) ) then
            !write(*,*) "Collapsing face", f
            ! Select the smallest edge of the triangle to remove
            e = minloc(  eLengths( edge_of_face(1:3,f) ) , 1  )
            e = edge_of_face(e,f)

            !write(*,*) "Collapsing edge", e
            
            v1 = vert_of_edge(1,e)
            v2 = vert_of_edge(2,e)

            ! Compute new optimal position of (v1,v2) collapse, vbar, using QEM
            ! Reset quadrics
            Q1(1:4,1:4) = 0.0d0
            Q2(1:4,1:4) = 0.0d0
            call calc_errorQuadric_of_v(Q1,v1,e,ne,nf,nv,tri_nor,xyz,edge_of_face,face_of_edge,vert_of_edge)
            call calc_errorQuadric_of_v(Q2,v2,e,ne,nf,nv,tri_nor,xyz,edge_of_face,face_of_edge,vert_of_edge)

            call solve_QEM_system(Q1,Q2, xyz(1:3,v1) ) ! v1 -> vbar
            !write(*,*) "New position from QEM of ecol ", e, "is ", xyz(1:3,v1)

            !-------------- Apply ghost flags --------------------------------
            ! Faces to be removed: the faces of edge e
            F1 = face_of_edge(1,e)
            F2 = face_of_edge(2,e)
            !write(*,*) "Faces to ghost of ecol ", e, "is ", F1, F2

            ! Remove 1 edge from each of face F1 and F2, not counting the edge e, call this e1, e2
            ! Also retain 1 edge (not counting e) from each of F1, F2, call this re1, re2
            call get_ghost_and_retained_edges(F1,e,nf,edge_of_face,e1,re1)
            call get_ghost_and_retained_edges(F2,e,nf,edge_of_face,e2,re2)

            !write(*,*) "ecol ", e, "e2 is ", e2, "re2 is ", re2

            isGhostFace(F1) = .true.
            isGhostFace(F2) = .true.

            isGhostVert(v2) = .true.

            isGhostEdge(e) = .true.
            isGhostEdge(e1) = .true.
            isGhostEdge(e2) = .true.

            !-------------- Update connectivity --------------------------------
            call update_face_connectivity(v1,v2,e1,e2,re1,re2,nf,vert_of_face,edge_of_face)
            call update_edge_connectivity(v1,v2,e1,e2,re1,re2,F1,F2,ne,vert_of_edge,face_of_edge)

            ! Optimisation steps: tangential relaxation
            !call write_tecplot_geom ! for debugging
            !call debug_write(flip_cnt) ! for debugging

            !if ( (ismaster) .and. (cnt .ne. 0) ) then
            !    write(*,*) "Entering 1-ring edge-flip for vertex", v1
            !    call debug_write(flip_cnt) ! for debugging
            !endif

            !do niter = 1,4
                !call calculate_eLengths(eLengths,nv,ne,xyz,vert_of_edge,isGhostEdge)
                !call calculate_area(Surface,nv,nf,xyz,vert_of_face,sur,isGhostFace,rm_flag,A_thresh) ! Update sur
                !call calculate_skewness (ne,nf,edge_of_face,sur,eLengths,skewness,isGhostFace)
                !call optimiseSkewness_1ring(cnt,v1,re1,ne,nf,nv,vert_of_face,face_of_edge,edge_of_face,vert_of_edge,xyz)
                !flip_cnt = flip_cnt + cnt
            !enddo

            !flip_cnt = flip_cnt + cnt
            !if ( (ismaster) .and. (cnt .ne. 0) ) then
            !    write(*,*) "Exited 1-ring edge-flip with", cnt, "flips!"
            !    call debug_write(flip_cnt) ! for debugging
            !endif

            ! Flag un-anchored vertices, for coupling to smoothing
            call flag_neighbours_1ring(anchorVert,flagged_edge,v1,re1,nv,ne,nf,face_of_edge,vert_of_edge,edge_of_face)


            !-------------- Update relevant geometric information --------------

            !Reset flag after remesh,  will be modified to .true. again if small triangle area detected
            rm_flag = .false.

            call calculate_eLengths(eLengths,nv,ne,xyz,vert_of_edge,isGhostEdge)
            call calculate_area(Surface,nv,nf,xyz,vert_of_face,sur,isGhostFace,rm_flag,A_thresh) ! Update sur
            call update_tri_normal (tri_nor,nv,nf,xyz,vert_of_face,isGhostFace)
            call calculate_skewness (ne,nf,edge_of_face,sur,eLengths,skewness,isGhostFace,rm_flag,skew_thresh)

            !!---------------------- Equalize valences in 1-ring of v1 with edge-flips-----------------


            ! Euler's polyhedral number should be two for closed manifold
            EN = count(isGhostVert .eqv. .false.) -  count(isGhostEdge .eqv. .false.) +  count(isGhostFace .eqv. .false.)
            if (EN .ne. 2) then
                write(*,*) "Topology not preserved, error in remeshing, writing geom and exiting program now!"
                call write_tecplot_geom
                !stop
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
                call MPI_Finalize(ierr)
            endif

            if (count(isGhostFace .eqv. .false.) .le. 4 ) then
                write(*,*) "Geometry reduced to tetrahedron, writing geom and exiting program now!"
                call write_tecplot_geom
                !stop
                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
                call MPI_Finalize(ierr)
            endif

            ecol_cnt = ecol_cnt + 1

        endif
    enddo
enddo !while

! if (ismaster) then
!     write(*,*) ecol_cnt," edge collapses performned "
!     write(*,*) flip_cnt," edge flips performned "
! endif

! if (ismaster) then
!     filename = 'continuation/isGhostVert.h5'
!     dsetname = trim('isGhostVert')
!     call HdfWriteSerialInt2D(filename,dsetname,nv,1,isGhostVert)

!     filename = 'continuation/isGhostEdge.h5'
!     dsetname = trim('isGhostEdge')
!     call HdfWriteSerialInt2D(filename,dsetname,ne,1,isGhostEdge)

!     filename = 'continuation/isGhostFace.h5'
!     dsetname = trim('isGhostFace')
!     call HdfWriteSerialInt2D(filename,dsetname,nf,1,isGhostFace)

!     filename = 'continuation/anchorVert.h5'
!     dsetname = trim('anchorVert')
!     call HdfWriteSerialInt2D(filename,dsetname,nv,1,anchorVert)

!     filename = 'continuation/flagged_edge.h5'
!     dsetname = trim('flagged_edge')
!     call HdfWriteSerialInt2D(filename,dsetname,ne,1,flagged_edge)

!     filename = 'continuation/vert_of_edge.h5'
!     dsetname = trim('vert_of_edge')
!     call HdfWriteSerialInt2D(filename,dsetname,2,ne,vert_of_edge)

!     filename = 'continuation/vert_of_face.h5'
!     dsetname = trim('vert_of_face')
!     call HdfWriteSerialInt2D(filename,dsetname,3,nf,vert_of_face)

!     filename = 'continuation/edge_of_face.h5'
!     dsetname = trim('edge_of_face')
!     call HdfWriteSerialInt2D(filename,dsetname,3,nf,edge_of_face)

!     filename = 'continuation/face_of_edge.h5'
!     dsetname = trim('face_of_edge')
!     call HdfWriteSerialInt2D(filename,dsetname,2,ne,face_of_edge)

!     filename = 'continuation/xyz.h5'
!     dsetname = trim('xyz')
!     call HdfWriteSerialReal2D(filename,dsetname,3,nv,xyz)
! endif

! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
! call MPI_Finalize(ierr)

end subroutine main_remesh
!------------------------------------------------------
subroutine calc_errorQuadric_of_v(Q,v,e,ne,nf,nv,tri_nor,xyz,edge_of_face,face_of_edge,vert_of_edge)
    ! Accumualte the error quadric of a vertex, v, based on its adjacent faces
    ! Face-adjacency is queried based on an arbitrary incident edge e of vertex v
    ! I.e. sum the fundamental quadrics, eqn. (2) in Garland & Heckbert (1997)
    implicit none
    integer :: v, e, ne, nv, nf
    integer, dimension(3,nf) :: edge_of_face
    integer, dimension(2,ne) :: face_of_edge, vert_of_edge
    real, dimension(3,nf) :: tri_nor
    real, dimension(3,nv) :: xyz
    real, dimension(4:4) :: Q ! Error quadric of vertex v
    integer :: prevFace, currentFace, currentEdge, prevEdge, F1
    real, dimension(4) :: p


    ! Accmulate the face adjacency
    ! Arbitrary starting face
    F1 = face_of_edge(1,e)
    prevFace = F1 
    currentFace = 0

    currentEdge = 0
    prevEdge = e

    ! First, accumulate fundamental quadric contribution from face F1
    p(1:3) = tri_nor(1:3,F1)
    p(4) = -dot_product( tri_nor(1:3,F1) , xyz(1:3,v) )
    p(1:4) = p(1:4) / norm2(p)
    call DGEMM('N','T',4,4,1,1.0d0,p,4,p,4, 1.0d0,Q,4)

    do while (currentEdge .ne. e )

        if (face_of_edge(1,prevEdge) .ne. prevFace) then
            currentFace = face_of_edge(1,prevEdge)
        else
            currentFace = face_of_edge(2,prevEdge)
        endif

        if (currentFace .eq. F1) exit

        ! Accumulate fundamental quadric contribution from currentFace
        p(1:3) = tri_nor(1:3,currentFace)
        p(4) = -dot_product( tri_nor(1:3,currentFace) , xyz(1:3,v) )
        p(1:4) = p(1:4) / norm2(p)
        call DGEMM('N','T',4,4,1,1.0d0,p,4,p,4, 1.0d0,Q,4)

        ! Find next edge, store as currentEdge
        call get_next_edge_of_v(v,nf,ne,edge_of_face,currentFace,vert_of_edge,prevEdge,currentEdge)

        prevFace = currentFace
        prevEdge = currentEdge

    enddo 



end subroutine calc_errorQuadric_of_v

    subroutine solve_QEM_system(Q1,Q2,xyz_new)
        ! Solve the QEM linear system to obtain new optimal position, xyz_newnew
        implicit none
        real, dimension(4,4) :: Q1, Q2
        real, dimension(3) :: xyz_new
        real, dimension(4) :: b
        integer, dimension(4) :: IPIV
        integer :: INFO
        Q1 = Q1 + Q2 ! Qbar in Garland & Heckbert notation

        ! Setup eqn. 1 in Garland & Heckbert (1997)
        b = [0.0d0, 0.0d0, 0.0d0, 1.0d0] ! RHS vector
        Q1(4,1:4) = [0.0d0, 0.0d0, 0.0d0, 1.0d0]

        ! Solve 
        call  dgesv	(4,1,Q1,4,IPIV,b,4,INFO)

        ! ToDo: some exception-handling if INFO .ne. 0 (unsuccessful)

        xyz_new(1:3) = b(1:3)
    end subroutine solve_QEM_system


subroutine get_ghost_and_retained_edges(f,e,nf,edge_of_face,ghostEdge,retainedEdge)
    ! Given face f and collapsed edge e, select a further edge to ghost (that is not e) and edge to retain
    implicit none
    integer :: i,f,e,ghostEdge,retainedEdge, nf, cnt
    integer, dimension(3,nf) :: edge_of_face

    integer, dimension(3) :: ebuffer
    
    cnt = 1

    do i = 1,3
        if ( edge_of_face(i,f) .ne. e ) then
            ebuffer(cnt) = edge_of_face(i,f)
            cnt = cnt + 1
        endif
    enddo

    ! Assert for cnt

    ! Arbitrarily set as first to be ghosted, second to be retained
    ghostEdge = ebuffer(1)
    retainedEdge = ebuffer(2)
end

subroutine update_face_connectivity(v1,v2,e1,e2,re1,re2,nf,vert_of_face,edge_of_face)
    implicit none
    integer :: v1, v2, e1, e2, re1, re2, nf
    integer, dimension(3,nf) :: vert_of_face, edge_of_face


    ! All faces prev. connected to v2 now connected to v1
    where (vert_of_face .eq. v2 ) vert_of_face = v1

    ! All faces connected to the ghosted edges e1/e2 now connected to the retained edge re1/re2
    where (edge_of_face .eq. e1 ) edge_of_face = re1
    where (edge_of_face .eq. e2 ) edge_of_face = re2

end subroutine update_face_connectivity

subroutine update_edge_connectivity(v1,v2,e1,e2,re1,re2,F1,F2,ne,vert_of_edge,face_of_edge)
    implicit none
    integer :: v1, v2, e1, e2, re1, re2, ne, F1, F2
    integer :: i, cnt
    integer, dimension(2,ne) :: vert_of_edge, face_of_edge
    integer, dimension(4) :: fbuffer

    ! All edges prev. connected to v2 now connected to v1
    where (vert_of_edge .eq. v2 ) vert_of_edge = v1

    ! Face conncection: faces connected to the ghosted edges e1/e2 now connected to the retained edge re1/re2
    fbuffer(1:2) = face_of_edge(1:2,e1)
    fbuffer(3:4) = face_of_edge(1:2,re1)

    cnt = 1
    do i = 1,4
        if ( fbuffer(i) .ne. F1 ) then
            face_of_edge(cnt,re1) = fbuffer(i)
            cnt = cnt + 1
        endif
    enddo

    !check cnt should be 3 at end

    fbuffer(1:2) = face_of_edge(1:2,e2)
    fbuffer(3:4) = face_of_edge(1:2,re2)

    cnt = 1
    do i = 1,4
        if ( fbuffer(i) .ne. F2 ) then
            face_of_edge(cnt,re2) = fbuffer(i)
            cnt = cnt + 1
        endif
    enddo

    !check cnt should be 3 at end


end subroutine update_edge_connectivity

subroutine get_next_edge_of_v(v,nf,ne,edge_of_face,currentFace,vert_of_edge,prevEdge,nextEdge)
    ! Retrive the next edge that is not prevEdge joined to vertex v
    implicit none
    integer :: v, ne, nf,i
    integer, dimension(3,nf) :: edge_of_face
    integer, dimension(2,ne) ::  vert_of_edge
    integer :: currentFace, nextEdge, prevEdge
    integer, dimension(3) :: edges
    logical :: success
    
    success = .false.

    edges(1:3) = edge_of_face(1:3,currentFace)

    do i = 1,3
        if (edges(i) .ne. prevEdge) then
            if( ( vert_of_edge(1, edges(i) ) .eq. v ) .or. ( vert_of_edge(2, edges(i) ) .eq. v ) ) then
                nextEdge = edges(i)
                success = .true.
            endif
        endif
    enddo

    if (success .eqv. .false.) then
        write(*,*) "Error finding next edge of vertex", v
        stop
    endif
    end subroutine get_next_edge_of_v



subroutine flag_neighbours_1ring(anchorVert,flagged_edge,v,e,nv,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
        use param, only: ismaster
        ! Flag the neighbours (verts + edges) in the 1-ring neighbourhood
        ! of vertex v to be un-anchored for later smoothing
        implicit none
        integer :: v, e, nv,ne, nf
        integer ::  cnt
        logical, dimension(nv) :: anchorVert
        logical, dimension(ne) :: flagged_edge
        integer, dimension(3,nf) :: edge_of_face
        integer, dimension(2,ne) :: face_of_edge, vert_of_edge
        integer :: prevFace, currentFace, currentEdge, prevEdge, F1
    
        ! 1-ring center
        anchorVert(v) = .false.
        flagged_edge(e) = .true.

        ! Accmulate the face adjacency
        ! Arbitrary starting face
        F1 = face_of_edge(1,e)
        prevFace = F1 
        currentFace = 0
    
        currentEdge = 0
        prevEdge = e
    
        cnt = 1
    
        do while (currentEdge .ne. e )
    
            if (face_of_edge(1,prevEdge) .ne. prevFace) then
                currentFace = face_of_edge(1,prevEdge)
            else
                currentFace = face_of_edge(2,prevEdge)
            endif
    
            if (currentFace .eq. F1) exit
    
            cnt = cnt + 1
    
            ! Find next edge, store as currentEdge
            call get_next_edge_of_v(v,nf,ne,edge_of_face,currentFace,vert_of_edge,prevEdge,currentEdge)
            anchorVert( vert_of_edge(:,currentEdge)  ) = .false.
            flagged_edge(currentEdge) = .true.
            prevFace = currentFace
            prevEdge = currentEdge
        enddo 
end subroutine flag_neighbours_1ring


! subroutine optimiseSkewness_1ring(flip_cnt,v,e,ne,nf,nv,vert_of_face,face_of_edge,edge_of_face,vert_of_edge,xyz)
!     ! Do a sequence of edge flips around the 1-ring neighbourhood of vertex v
!     use param, only: ismaster
!     !use geom
!     implicit none
!     integer :: v, e, ne, nv, nf, i,j,k
!     integer, dimension(3,nf) :: vert_of_face, edge_of_face
!     integer, dimension(2,ne) :: face_of_edge, vert_of_edge
!     real, dimension(3,nv) :: xyz
!     real, dimension(3) :: vecBuffer
!     integer :: prevFace, currentFace, currentEdge, prevEdge, startFace
!     integer :: v1, v2, v3, v4, F1, F2
!     integer :: flip_cnt
!     real :: skew_pre, skew_post
!     real, dimension(2) :: perim, sur_opt, sur_buffer
!     logical :: flip_OK
!     integer :: v1_val, v2_val, v3_val, v4_val
!     integer :: e_of_v1, e_of_v2, e_of_v3, e_of_v4
!     integer :: debug_cnt
    
!     integer, allocatable, dimension(:) :: eNeighbours_of_v

!     flip_cnt = 0

!     skew_pre = 0.
!     skew_post = 0.
!     perim = 0.
!     sur_opt = 0.
!     sur_buffer = 0.

!     ! First start off by calculating the valence of the vertex v: no. of neighbours
!     call get_vertValence(v1_val,v,e,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!     allocate(eNeighbours_of_v(v1_val))

!     ! Populate array
!     call get_vertEdges(v1_val,eNeighbours_of_v,v,e,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!     ! Now loop over each edge
!     do j = 1,v1_val
!         currentEdge = eNeighbours_of_v(j)

!         F1 = face_of_edge(1,currentEdge)
!         F2 = face_of_edge(2,currentEdge)

!         call get_quadrilateral_vertices(v1,v2,v3,v4,v,F1,F2,currentEdge,ne,nf,vert_of_face,vert_of_edge)

!         ! Now do a bunch of tests to see if the flip is OK
!         call get_vertValence(v1_val,v1,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!         call get_vertValence(v2_val,v2,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!         do i = 1,3
!             if ( count ( vert_of_edge(:,edge_of_face(i,F1) ) .eq. v3) .gt. 0 ) then
!                 e_of_v3 = edge_of_face(i,F1)
!             endif

!             if ( count ( vert_of_edge(:,edge_of_face(i,F2) ) .eq. v4) .gt. 0 ) then
!                 e_of_v4 = edge_of_face(i,F2)
!             endif
!         enddo

!         call get_vertValence(v3_val,v3,e_of_v3,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!         call get_vertValence(v4_val,v4,e_of_v4,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!         ! Do tests to see if edge flip is legal
!         call test_edgeFlip(flip_OK,v1,v2,v3,v4,v1_val,v2_val,v3_val,v4_val,xyz,ne,nf,nv)

!         if (flip_OK .eqv. .true.) then
!             ! Get skewness of F1 and F2
!                 !F1 = (v1,v2,v3)
!                 perim(1) = norm2( xyz(1:3,v2) - xyz(1:3,v1)  ) + norm2( xyz(1:3,v3) - xyz(1:3,v1)  ) + &
!                            norm2( xyz(1:3,v3) - xyz(1:3,v2)  )

!                 perim(2) = norm2( xyz(1:3,v2) - xyz(1:3,v1)  ) + norm2( xyz(1:3,v4) - xyz(1:3,v1)  ) + &
!                 norm2( xyz(1:3,v4) - xyz(1:3,v2)  )

!                 call cross(vecBuffer, xyz(1:3,v2) - xyz(1:3,v1), xyz(1:3,v3) - xyz(1:3,v1) )
!                 sur_buffer(1) = 0.5 * norm2( vecBuffer)
!                 call cross(vecBuffer, xyz(1:3,v2) - xyz(1:3,v1), xyz(1:3,v4) - xyz(1:3,v1) )
!                 sur_buffer(2) = 0.5 * norm2( vecBuffer)

!                 !sur_buffer(1) = 0.5 * norm2(   cross(xyz(1:3,v2) - xyz(1:3,v1), xyz(1:3,v3) - xyz(1:3,v1))    )
!                 !sur_buffer(2) = 0.5 * norm2(   cross(xyz(1:3,v2) - xyz(1:3,v1), xyz(1:3,v4) - xyz(1:3,v1))    )

!                 sur_opt = sqrt(3.0) / 4.0 * (perim / 3.0)**2

!                 skew_pre = sum ( ( sur_opt - sur_buffer ) / sur_opt )

!             ! Perform edge flip
!             call edgeFlip(currentEdge,F1,F2,v1,v2,v3,v4,ne,nf,vert_of_edge,vert_of_face,edge_of_face,face_of_edge)

!             ! Compute new valences, original e is now connected to (v3,v4)
!             call get_vertValence(v3_val,v3,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!             call get_vertValence(v4_val,v4,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!             do i = 1,3
!                 if ( count ( vert_of_edge(:,edge_of_face(i,F1) ) .eq. v1) .gt. 0 ) then
!                     e_of_v1 = edge_of_face(i,F1)
!                 endif
            
!                 if ( count ( vert_of_edge(:,edge_of_face(i,F2) ) .eq. v2) .gt. 0 ) then
!                     e_of_v2 = edge_of_face(i,F2)
!                 endif
!             enddo

!             call get_vertValence(v1_val,v1,e_of_v1,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!             call get_vertValence(v2_val,v2,e_of_v2,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!             ! After flip: F1(v3,v4,v1)  F2(v3,v4,v2)

!             perim(1) = norm2( xyz(1:3,v3) - xyz(1:3,v4)  ) + norm2( xyz(1:3,v3) - xyz(1:3,v1)  ) + &
!             norm2( xyz(1:3,v4) - xyz(1:3,v1)  )

!             perim(2) = norm2( xyz(1:3,v3) - xyz(1:3,v4)  ) + norm2( xyz(1:3,v3) - xyz(1:3,v2)  ) + &
!             norm2( xyz(1:3,v4) - xyz(1:3,v2)  )

!             call cross(vecBuffer, xyz(1:3,v4) - xyz(1:3,v3), xyz(1:3,v4) - xyz(1:3,v1) )
!             sur_buffer(1) = 0.5 * norm2( vecBuffer)
!             call cross(vecBuffer, xyz(1:3,v4) - xyz(1:3,v3), xyz(1:3,v4) - xyz(1:3,v2) )
!             sur_buffer(2) = 0.5 * norm2( vecBuffer)
!             sur_opt = sqrt(3.0) / 4.0 * (perim / 3.0)**2

!             skew_post = sum ( ( sur_opt - sur_buffer ) / sur_opt )


!             if (skew_pre .le. skew_post) then ! Flip back if flip did not optimize valences
!                 call edgeFlip(currentEdge,F1,F2,v3,v4,v1,v2,ne,nf,vert_of_edge,vert_of_face,edge_of_face,face_of_edge)
!             else
!                 flip_cnt = flip_cnt + 1
!             endif
!         endif

!     enddo

!     deallocate(eNeighbours_of_v)


! end subroutine optimiseSkewness_1ring

! subroutine equalizeValences_1ring(flip_cnt,v,e,ne,nf,nv,vert_of_face,face_of_edge,edge_of_face,vert_of_edge,xyz)
!     ! Do a sequence of edge flips around the 1-ring neighbourhood of vertex v
!     use param, only: ismaster
!     implicit none
!     integer :: v, e, ne, nv, nf, i,j
!     integer, dimension(3,nf) :: vert_of_face, edge_of_face
!     integer, dimension(2,ne) :: face_of_edge, vert_of_edge
!     real, dimension(3,nv) :: xyz
!     integer :: prevFace, currentFace, currentEdge, prevEdge, startFace
!     integer :: v1, v2, v3, v4, F1, F2
!     integer :: flip_cnt
!     integer :: deviation_pre, deviation_post
!     logical :: flip_OK
!     integer :: v1_val, v2_val, v3_val, v4_val
!     integer :: e_of_v1, e_of_v2, e_of_v3, e_of_v4
!     integer :: debug_cnt
    
!     integer, allocatable, dimension(:) :: eNeighbours_of_v

!     flip_cnt = 0

!     ! First start off by calculating the valence of the vertex v: no. of neighbours
!     call get_vertValence(v1_val,v,e,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!     allocate(eNeighbours_of_v(v1_val))

!     ! Populate array
!     call get_vertEdges(v1_val,eNeighbours_of_v,v,e,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!     ! Now loop over each edge
!     do j = 1,v1_val
!         currentEdge = eNeighbours_of_v(j)

!         F1 = face_of_edge(1,currentEdge)
!         F2 = face_of_edge(2,currentEdge)

!         call get_quadrilateral_vertices(v1,v2,v3,v4,v,F1,F2,currentEdge,ne,nf,vert_of_face,vert_of_edge)

!         ! Now do a bunch of tests to see if the flip is OK
!         call get_vertValence(v1_val,v1,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!         call get_vertValence(v2_val,v2,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!         do i = 1,3
!             if ( count ( vert_of_edge(:,edge_of_face(i,F1) ) .eq. v3) .gt. 0 ) then
!                 e_of_v3 = edge_of_face(i,F1)
!             endif

!             if ( count ( vert_of_edge(:,edge_of_face(i,F2) ) .eq. v4) .gt. 0 ) then
!                 e_of_v4 = edge_of_face(i,F2)
!             endif
!         enddo

!         call get_vertValence(v3_val,v3,e_of_v3,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!         call get_vertValence(v4_val,v4,e_of_v4,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!         ! Do tests to see if edge flip is legal
!         call test_edgeFlip(flip_OK,v1,v2,v3,v4,v1_val,v2_val,v3_val,v4_val,xyz,ne,nf,nv)

!         if (flip_OK .eqv. .true.) then
!             deviation_pre = abs(v1_val - 6) + abs(v2_val - 6) + &
!                         abs(v3_val - 6) + abs(v4_val - 6) 

!             ! Perform edge flip
!             call edgeFlip(currentEdge,F1,F2,v1,v2,v3,v4,ne,nf,vert_of_edge,vert_of_face,edge_of_face,face_of_edge)

!             ! Compute new valences, original e is now connected to (v3,v4)
!             call get_vertValence(v3_val,v3,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!             call get_vertValence(v4_val,v4,currentEdge,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!             do i = 1,3
!                 if ( count ( vert_of_edge(:,edge_of_face(i,F1) ) .eq. v1) .gt. 0 ) then
!                     e_of_v1 = edge_of_face(i,F1)
!                 endif
            
!                 if ( count ( vert_of_edge(:,edge_of_face(i,F2) ) .eq. v2) .gt. 0 ) then
!                     e_of_v2 = edge_of_face(i,F2)
!                 endif
!             enddo

!             call get_vertValence(v1_val,v1,e_of_v1,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!             call get_vertValence(v2_val,v2,e_of_v2,ne,nf,face_of_edge,vert_of_edge,edge_of_face)

!             deviation_post = abs(v1_val - 6) + abs(v2_val - 6) + &
!                             abs(v3_val - 6) + abs(v4_val - 6) 

!             if (deviation_pre .le. deviation_post) then ! Flip back if flip did not optimize valences
!                 call edgeFlip(currentEdge,F1,F2,v3,v4,v1,v2,ne,nf,vert_of_edge,vert_of_face,edge_of_face,face_of_edge)
!             else
!                 flip_cnt = flip_cnt + 1
!             endif
!         endif

!     enddo

!     deallocate(eNeighbours_of_v)


! end subroutine equalizeValences_1ring


! subroutine get_quadrilateral_vertices(v1,v2,v3,v4,v,F1,F2,e,ne,nf,vert_of_face,vert_of_edge)
!     ! Given faces F1 and F2 which are defined to be adjacent by edge e, calculate the bounding quadrilateral vertices,
!     ! (v1,v2,v3,v4)
!     ! Edge e is defined as (v1,v2)
!     ! v3 is joined to F1, v4 is joined to F2

!         ! Index naming convetion used:
!     !          v3                                
!     !          O  
!     !        /   \
!     !       /     \
!     !      /   F1  \
!     !     /         \
!     ! v1 O-----e-----O v2
!     !     \ F2      / 
!     !      \       /  
!     !       \     /   
!     !        \   /   
!     !          O 
!     !          v4 

!     implicit none
!     integer :: v, e, ne, nf, i
!     integer, dimension(3,nf) :: vert_of_face
!     integer, dimension(2,ne) :: vert_of_edge
!     integer ::  v1, v2, v3, v4, F1, F2

!     ! Define v3 and v4 such that
!     ! F1 is (v1,v2,v3)
!     ! F2 is (v1,v2,v4)

!     ! Test edge flip
!     ! Fix v1 == v always
!     v1 = v 
!     if (vert_of_edge(1,e) .eq. v) then
!         v2 = vert_of_edge(2,e)
!     else
!         v2 = vert_of_edge(1,e)
!     endif

!     do i = 1,3
!         if ( (vert_of_face(i,F1) .ne. v) .and. (vert_of_face(i,F1) .ne. v2) ) then
!             v3 = vert_of_face(i,F1)
!         endif

!         if ( (vert_of_face(i,F2) .ne. v) .and. (vert_of_face(i,F2) .ne. v2) ) then
!             v4 = vert_of_face(i,F2)
!         endif
!     enddo

! end subroutine get_quadrilateral_vertices

! subroutine get_vertValence(valence,v,e,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!     use param, only: ismaster
!     ! Compute the valence of vertex v
!     implicit none
!     integer :: v, e, ne, nf
!     integer, dimension(3,nf) :: edge_of_face
!     integer, dimension(2,ne) :: face_of_edge, vert_of_edge
!     integer :: prevFace, currentFace, currentEdge, prevEdge, F1
!     integer :: valence

!     ! Accmulate the face adjacency
!     ! Arbitrary starting face
!     F1 = face_of_edge(1,e)
!     prevFace = F1 
!     currentFace = 0

!     currentEdge = 0
!     prevEdge = e

!     valence = 1

!     do while (currentEdge .ne. e )

!         if (face_of_edge(1,prevEdge) .ne. prevFace) then
!             currentFace = face_of_edge(1,prevEdge)
!         else
!             currentFace = face_of_edge(2,prevEdge)
!         endif

!         if (currentFace .eq. F1) exit

!         valence = valence + 1

!         ! Find next edge, store as currentEdge
!         call get_next_edge_of_v(v,nf,ne,edge_of_face,currentFace,vert_of_edge,prevEdge,currentEdge)

!         prevFace = currentFace
!         prevEdge = currentEdge
!     enddo 
! end subroutine get_vertValence

! subroutine get_vertEdges(valence,eNeighbours,v,e,ne,nf,face_of_edge,vert_of_edge,edge_of_face)
!     use param, only: ismaster
!     ! Compute the valence of vertex v
!     implicit none
!     integer :: v, e, ne, nf
!     integer :: valence, cnt
!     integer, dimension(3,nf) :: edge_of_face
!     integer, dimension(2,ne) :: face_of_edge, vert_of_edge
!     integer, dimension(valence) :: eNeighbours
!     integer :: prevFace, currentFace, currentEdge, prevEdge, F1

!     ! Accmulate the face adjacency
!     ! Arbitrary starting face
!     F1 = face_of_edge(1,e)
!     prevFace = F1 
!     currentFace = 0

!     currentEdge = 0
!     prevEdge = e

!     cnt = 1

!     eNeighbours(cnt) = e


!     do while (currentEdge .ne. e )

!         if (face_of_edge(1,prevEdge) .ne. prevFace) then
!             currentFace = face_of_edge(1,prevEdge)
!         else
!             currentFace = face_of_edge(2,prevEdge)
!         endif

!         if (currentFace .eq. F1) exit

!         cnt = cnt + 1

!         ! Find next edge, store as currentEdge
!         call get_next_edge_of_v(v,nf,ne,edge_of_face,currentFace,vert_of_edge,prevEdge,currentEdge)
!         eNeighbours(cnt) = currentEdge

!         prevFace = currentFace
!         prevEdge = currentEdge
!     enddo 
! end subroutine get_vertEdges

! subroutine test_edgeFlip(flip_OK,v1,v2,v3,v4,v1_val,v2_val,v3_val,v4_val,xyz,ne,nf,nv)
!     ! Compute the valence of vertex v
!     use param, only: ismaster
!     implicit none
!     integer :: ne, nf,nv
!     real, dimension(3,nv) :: xyz
!     integer :: v1,v2,v3,v4
!     integer :: v1_val, v2_val, v3_val, v4_val
!     logical :: flip_OK
!     real :: pi, angle1, angle2
!     real, dimension(3) :: d12, d13, d23, d14, d24


!     flip_OK = .true.

!     ! Do a series of tests to see if an edge-flip operation is permissable

!     ! First, flip is only permissiable if valence of v1,v2 vertices is greater than 3
!     if ( (v1_val .le. 3) .or. (v2_val .le. 3) ) then
!         !write(*,*) "Valence <=3 detected, no edge flip"
!         flip_OK = .false.
!     endif

!     ! Bounding quadrilaterial must be convex: so no internal angle >180 degrees allowed
!     if (flip_OK .eqv. .true.) then
!         pi = 4.0d0 * atan(1.0d0)
!         d12 =   xyz(:,v2) - xyz(:,v1) 

!         ! (angle between e(v1, v2) and e(v3,v1)
!         d13 =   xyz(:,v3) - xyz(:,v1) 
!         d23 =   xyz(:,v3) - xyz(:,v2) 

!         !Note the signs
!         angle1 = acos ( dot_product(d13,d12) / ( norm2(d13)*norm2(d12)  ) ) * (180.0d0 / pi)
!         angle2 = acos ( dot_product(d23,-d12) / ( norm2(d23)*norm2(d12)  ) ) * (180.0d0 / pi)

    
!         d14 = xyz(:,v4) - xyz(:,v1) 
!         d24 = xyz(:,v4) - xyz(:,v2) 

!         angle1 = angle1 + acos ( dot_product(d14,d12) / ( norm2(d14)*norm2(d12)  ) ) * (180.0d0/pi)
!         angle2 = angle2 + acos ( dot_product(d24,-d12) / ( norm2(d14)*norm2(d12)  ) ) * (180.0d0/pi)

!         if ( (angle1 .ge. 180.0d0) .or. (angle2 .ge. 180.0d0) ) then
!             !if (ismaster) then
!             !    write(*,*) "Concavity detected, no edge flip, angles are", angle1, angle2
!             !endif
!             flip_OK = .false.
!         endif
!     endif

!     ! Dihedral angle test: skip for now

! end subroutine test_edgeFlip


! subroutine edgeFlip(e,F1,F2,v1,v2,v3,v4,ne,nf,vert_of_edge,vert_of_face,edge_of_face,face_of_edge)
!     ! Perform an edge-flip for edge e
!     ! Index naming convetion used:

!     !          v3                                 v3
!     !          O                                  O
!     !        /   \                              / | \
!     !       /     \                            /  |  \
!     !      /   F1  \                          /   |   \
!     !     /         \          eflip         /    |    \
!     ! v1 O-----e-----O v2    -------->   v1 O  F1 e F2  O  v2
!     !     \ F2      /                        \    |    /
!     !      \       /                          \   |   /
!     !       \     /                            \  |  /
!     !        \   /                              \ | /
!     !          O                                  O
!     !          v4                                 v4

!     ! Rules:
!     ! - Edge to be flipped, e has original vertices (v1,v2)
!     ! - Original v3 and v4 belong to faces F1 and F2 respectively
!     ! Upon flipping:
!     !       - Flipped edge e, will be (v3,v4)
!     !       - New F1 will include vertex v1, new F2 will include vertex v2
!     use param, only: ismaster
!     implicit none
!     integer :: ne, nf, i,j
!     integer :: v1,v2,v3,v4,e, F1, F2
!     integer :: cnt_v1, cnt_v2
!     integer, dimension(2) :: edges_of_v1, edges_of_v2
!     integer, dimension(2,ne) :: vert_of_edge, face_of_edge
!     integer, dimension(3,nf) :: vert_of_face, edge_of_face
!     integer, dimension(6) :: edges
!     logical, dimension(2,2) :: mask

!     ! Full set of edges in the bounding quadrilateral
!     edges(1:3) = edge_of_face(1:3,F1)
!     edges(4:6) = edge_of_face(1:3,F2)

!     cnt_v1 = 1
!     cnt_v2 = 1


!     do i = 1,6
!         if ( ( edges(i) .ne. e ) .and. ( count( vert_of_edge(:,edges(i)) .eq. v1  ) .ne. 0 )    ) then
!             edges_of_v1(cnt_v1) = edges(i)
!             cnt_v1 = cnt_v1 + 1
!         endif

!         if ( ( edges(i) .ne. e ) .and. ( count( vert_of_edge(:,edges(i)) .eq. v2  ) .ne. 0 )    ) then
!             edges_of_v2(cnt_v2) = edges(i)
!             cnt_v2 = cnt_v2 + 1
!         endif
!     enddo

!     ! Connectivity updates after performing the flip
!     vert_of_edge(:,e) = [v3,v4]

!     vert_of_face(1:3,F1) = [v3,v4,v1]
!     vert_of_face(1:3,F2) = [v3,v4,v2]

!     edge_of_face(:,F1) = [e,edges_of_v1]
!     edge_of_face(:,F2) = [e,edges_of_v2]

!     !Peripherary edges in the bounding quadrilateral
!     do i = 1,2
!         do j = 1,2

!             if (face_of_edge(j, edges_of_v1(i)) .eq. F2) then
!                 face_of_edge(j, edges_of_v1(i)) = F1
!             endif

!             if (face_of_edge(j, edges_of_v2(i)) .eq. F1) then
!                 face_of_edge(j, edges_of_v2(i)) = F2
!             endif

!         enddo
!     enddo


! end subroutine edgeFlip
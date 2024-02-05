!------------------------------------------------------
subroutine main_remesh (Surface,sur,eLengths,nf,ne,nv,xyz,tri_nor,A_thresh,&
                        vert_of_face,edge_of_face,vert_of_edge,face_of_edge,&
                        isGhostFace,isGhostEdge,isGhostVert,rm_flag)
    !use param
    !use mls_param
    !use coll_mod
    implicit none
    integer :: nf,nv,ne,i
    real :: A_thresh, Surface
    integer, dimension (3,nf) :: vert_of_face, edge_of_face
    integer, dimension(2,ne) :: vert_of_edge, face_of_edge
    logical, dimension(nf) :: isGhostFace
    logical, dimension(ne) :: isGhostEdge
    logical, dimension(nv) :: isGhostVert
    logical :: rm_flag
    real, dimension(3,nv) :: xyz
    real, dimension(3,nf) :: tri_nor
    real, dimension(nf) :: sur
    real, dimension(ne) :: eLengths
    real, dimension(4,4) :: Q1, Q2 ! Error quadrics of vertices v1, v2
    integer :: f, e, v1, v2
    integer :: F1, F2, e1, e2, re1, re2

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

    do while (rm_flag .eqv. .true.) 
    do f = 1,nf
        if ( (isGhostFace(f) .neqv. .false.) .and.  ( sur(f) .le. A_thresh ) ) then
            ! Select the smallest edge of the triangle to remove
            e = minloc(  eLengths( edge_of_face(1:3,f) ) , 1  )
            e = edge_of_face(e,f)

            v1 = vert_of_edge(1,e)
            v2 = vert_of_edge(2,e)

            ! Compute new optimal position of (v1,v2) collapse, vbar, using QEM
            ! Reset quadrics
            Q1(1:4,1:4) = 0.0d0
            Q2(1:4,1:4) = 0.0d0
            call calc_errorQuadric_of_v(Q1,v1,e,ne,nf,nv,tri_nor,xyz,edge_of_face,face_of_edge,vert_of_edge)
            call calc_errorQuadric_of_v(Q2,v2,e,ne,nf,nv,tri_nor,xyz,edge_of_face,face_of_edge,vert_of_edge)

            call solve_QEM_system(Q1,Q2, xyz(1:3,v1) ) ! v1 -> vbar

            !-------------- Apply ghost flags --------------------------------
            ! Faces to be removed: the faces of edge e
            F1 = face_of_edge(1,e)
            F2 = face_of_edge(1,e)

            ! Remove 1 edge from each of face F1 and F2, not counting the edge e, call this e1, e2
            ! Also retain 1 edge (not counting e) from each of F1, F2, call this re1, re2
            call get_ghost_and_retained_edges(F1,e,nf,edge_of_face,e1,re1)
            call get_ghost_and_retained_edges(F2,e,nf,edge_of_face,e2,re2)

            isGhostFace(F1) = .true.
            isGhostFace(F2) = .true.

            isGhostVert(v2) = .true.
            isGhostEdge(e) = .true.
            isGhostEdge(e1) = .true.
            isGhostEdge(e2) = .true.

            !-------------- Update connectivity --------------------------------
            call update_face_connectivity(v1,v2,e1,e2,re1,re2,nf,vert_of_face,edge_of_face)
            call update_edge_connectivity(v1,v2,e1,e2,re1,re2,F1,F2,ne,vert_of_edge,face_of_edge)


            !-------------- Update relevant geometric information --------------

            !Reset flag after remesh,  will be modified to .true. again if small triangle area detected
            rm_flag = .false.

            ! Update edge lengths
            call calculate_eLengths(eLengths,nv,ne,xyz,vert_of_edge,isGhostEdge)
            call calculate_area(Surface,nv,nf,xyz,vert_of_face,sur,isGhostFace,rm_flag,A_thresh) ! Update sur
            call update_tri_normal (tri_nor,nv,nf,xyz,vert_of_face,isGhostFace)

            ! Optional: collapse-counter, assert topology-preservation

        endif
    enddo
enddo !while
!! reset remesh flag at end
!rm_flag = .false.

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

subroutine get_next_edge_of_v(v,nf,ne,edge_of_face,currentFace,vert_of_edge,prevEdge,nextEdge)
    ! Retrive the next edge that is not prevEdge joined to vertex v
    implicit none
    integer :: v, ne, nf,i
    integer, dimension(3,nf) :: edge_of_face
    integer, dimension(2,ne) ::  vert_of_edge
    integer :: currentFace, nextEdge, prevEdge
    integer, dimension(3) :: edges
    
    edges(1:3) = edge_of_face(1:3,currentFace)

    do i = 1,3
        if (edges(i) .ne. prevEdge) then
            if( ( vert_of_edge(1, edges(i) ) .eq. v ) .or. ( vert_of_edge(2, edges(i) ) .eq. v ) ) then
                nextEdge = edges(i)
            endif
        endif
    enddo
    end subroutine get_next_edge_of_v

    subroutine solve_QEM_system(Q1,Q2,xyz_new)
        ! Solve the QEM linear system to obtain new optimal position, xyz_newnew
        implicit none
        real, dimension(4,4) :: Q1, Q2
        real, dimension(3) :: xyz_new
        real, dimension(4) :: b = 0.0d0
        integer, dimension(4) :: IPIV
        integer :: INFO
        Q1 = Q1 + Q2 ! Qbar in Garland & Heckbert notation

        ! Setup eqn. 1 in Garland & Heckbert (1997)
        b(4) = 1.0d0 ! RHS vector
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


end
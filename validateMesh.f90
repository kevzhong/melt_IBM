subroutine repair_face_orientations(vert_of_edge, face_of_edge, vert_of_face, edge_of_face, isGhostEdge,isGhostFace,isGhostVert,&
                                     Nvert,Nedge, Nface)

    ! Repair/check face orientations with a breadth-first search

    implicit none
    integer, intent(in) :: Nvert,Nedge, Nface
    integer, dimension(2,Nedge), intent(in) :: vert_of_edge, face_of_edge
    integer, dimension(3,Nface), intent(inout) :: vert_of_face, edge_of_face
    logical, dimension(Nvert) :: isGhostVert
    logical, dimension(Nedge) :: isGhostEdge
    logical, dimension(Nface) :: isGhostFace
    logical, dimension(Nface) :: visited
    integer, dimension(3,Nface) :: face_neighbors
    integer :: i, j, e, f, f_curr, f_adj, stack_size
    integer, dimension(Nface) :: stack
    integer :: v1, v2, dir_curr, dir_adj
    integer, dimension(3) :: vf_curr, vf_adj
    integer :: seed, nflip

    visited = .false.
    face_neighbors = 0
    stack = 0

    ! Build face_neighbors: each face shares 3 edges with 3 neighbors
    do e = 1, Nedge
    if (.not. isGhostEdge(e)) then
        f = face_of_edge(1,e)
        if ( isGhostFace(f) .eqv. .true. ) error stop "Error: face f is ghosted!"

        if (f > 0) then
            do j = 1,3
                if (face_neighbors(j,f) == 0) then
                    face_neighbors(j,f) = face_of_edge(2,e)
                    exit
                endif
            enddo
        endif
        f = face_of_edge(2,e)
        if ( isGhostFace(f) .eqv. .true. ) error stop "Error: face f is ghosted!"
        if (f > 0) then
            do j = 1,3
                if (face_neighbors(j,f) == 0) then
                    face_neighbors(j,f) = face_of_edge(1,e)
                    exit
                endif
            enddo
        endif
    endif
    enddo


    nflip = 0

    ! Start DFS traversal from first non-ghost face
    do f = 1,Nface
        if (.not. isGhostFace(f ) ) seed = f
    enddo

    visited(seed) = .true.
    stack_size = 1
    stack(stack_size) = seed

    do while (stack_size > 0)
        f_curr = stack(stack_size)
        stack_size = stack_size - 1

        vf_curr = vert_of_face(:,f_curr)

        do i = 1,3
            e = edge_of_face(i, f_curr)
            !if ( isGhostEdge(e) .eqv. .true. ) error stop "Error: edge e is ghosted!"
            if ( isGhostEdge(e) .eqv. .true. ) then
                write(*,*) "edge ", e, "is ghosted for f_curr", f_curr
            endif

            f_adj = face_of_edge(1,e)
            if ( isGhostFace(f_adj) .eqv. .true. ) error stop "Error: face f_adj is ghosted!"

            if (f_adj == f_curr) f_adj = face_of_edge(2,e)

            if (f_adj == 0 .or. visited(f_adj)) cycle

            vf_adj = vert_of_face(:,f_adj)

            v1 = vert_of_edge(1,e)
            v2 = vert_of_edge(2,e)
            if ( isGhostVert(v1) .eqv. .true. ) error stop "Error: Vert v1 is ghosted!"
            if ( isGhostVert(v2) .eqv. .true. ) error stop "Error: Vert v2 is ghosted!"

            call get_edge_direction_in_face(dir_curr,vf_curr, v1, v2)
            call get_edge_direction_in_face(dir_adj,vf_adj, v1, v2)

            if (dir_curr == dir_adj) then
                ! Flip adjacent face
                call reverse_face(vert_of_face(:,f_adj))

                nflip = nflip + 1
            endif

            visited(f_adj) = .true.
            stack_size = stack_size + 1
            stack(stack_size) = f_adj
        enddo
    enddo


    write(*,*) "Finishing face orienations with", nflip, "face flips!"


end subroutine repair_face_orientations

! subroutine rebuild_edge_of_face(vert_of_face, vert_of_edge, edge_of_face, Nface, Nedge)
!     implicit none
!     integer, intent(in) :: Nface, Nedge
!     integer, dimension(3,Nface), intent(in) :: vert_of_face
!     integer, dimension(2,Nedge), intent(in) :: vert_of_edge
!     integer, dimension(3,Nface), intent(out) :: edge_of_face

!     integer :: f, i, e
!     integer :: v1, v2, a, b
!     logical :: found

!     do f = 1, Nface
!         do i = 1, 3
!             v1 = vert_of_face(i, f)
!             v2 = vert_of_face(mod(i,3)+1, f)

!             found = .false.
!             do e = 1, Nedge
!                 a = vert_of_edge(1,e)
!                 b = vert_of_edge(2,e)
!                 if ( (v1 == a .and. v2 == b) .or. (v1 == b .and. v2 == a) ) then
!                     edge_of_face(i, f) = e
!                     found = .true.
!                     exit
!                 endif
!             enddo

!             if (.not. found) then
!                 write(*,*) 'Error: edge not found for face ', f, ' vertices ', v1, v2
!                 stop
!             endif
!         enddo
!     enddo
! end subroutine rebuild_edge_of_face

! ! ----------------------------

! subroutine orient_edges_from_faces(vert_of_face, edge_of_face, face_of_edge, vert_of_edge, Nface, Nedge)
!     implicit none
!     integer, intent(in) :: Nface, Nedge
!     integer, dimension(3,Nface), intent(in) :: vert_of_face
!     integer, dimension(3,Nface), intent(in) :: edge_of_face
!     integer, dimension(2,Nedge), intent(inout) :: vert_of_edge
!     integer, dimension(2,Nedge), intent(in) :: face_of_edge

!     integer :: f, i, e, v1, v2
!     integer :: a, b
!     integer :: j

!     do f = 1, Nface
!         do i = 1, 3
!             e = edge_of_face(i, f)
!             v1 = vert_of_face(i, f)
!             v2 = vert_of_face(mod(i,3)+1, f)  ! CCW: i→i+1

!             ! Check if vert_of_edge matches face usage; if not, flip
!             a = vert_of_edge(1, e)
!             b = vert_of_edge(2, e)

!             if ((a == v2) .and. (b == v1)) then
!                 ! Edge used in reverse direction → flip it
!                 vert_of_edge(1,e) = v1
!                 vert_of_edge(2,e) = v2
!             endif

!             ! If edge does not match either (v1→v2 or v2→v1), something is wrong
!             if (.not. ((vert_of_edge(1,e) == v1 .and. vert_of_edge(2,e) == v2) .or. &
!                        (vert_of_edge(1,e) == v2 .and. vert_of_edge(2,e) == v1))) then
!                 write(*,*) "Error: edge", e, "does not match face", f
!                 stop
!             endif
!         enddo
!     enddo
! end subroutine orient_edges_from_faces



subroutine get_edge_direction_in_face(direction,vf, v1, v2)
    implicit none
    integer, intent(in) :: vf(3), v1, v2
    integer :: direction, i, j

    direction = 0
    do i = 1,3
        j = mod(i,3) + 1
        if (vf(i) == v1 .and. vf(j) == v2) then
            direction = +1
            return
        elseif (vf(i) == v2 .and. vf(j) == v1) then
            direction = -1
            return
        endif
    enddo
end subroutine get_edge_direction_in_face

! ----------------------------

subroutine reverse_face(vf)
    implicit none
    integer, dimension(3), intent(inout) :: vf
    integer :: tmp
    tmp = vf(2)
    vf(2) = vf(3)
    vf(3) = tmp
end subroutine reverse_face

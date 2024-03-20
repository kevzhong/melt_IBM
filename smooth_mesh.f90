!------------------------------------------------------
! subroutine main_smooth(nv,ne,nf,xyz,isGhostVert,isGhostEdge,isGhostFace,anchorVert,vert_of_edge,vert_of_face)

!     use mpih
!     use param, only: ismaster

!     implicit none
!     integer :: nv,ne,nf
!     integer :: nv_active, cnt, v1, v2, v1m, v2m, nua
!     integer :: i, icnt, jcnt, j
!     logical, dimension(ne) :: isGhostEdge
!     logical, dimension(nf) :: isGhostFace
!     logical, dimension(nv) :: isGhostVert, anchorVert
!     integer, dimension(2,ne) :: vert_of_edge
!     integer, dimension(3,nf) :: vert_of_face
!     !integer, dimension(nv) :: vert_mask
!     real, dimension(3,nv) :: xyz
!     real :: valence
!     integer :: INFO
!     real, allocatable, dimension(:,:) :: L, Lreduced ! Laplacian matrix
!     real, allocatable, dimension(:,:) :: b,buffer ! rhs vector
!     !real, allocatable, dimension(:) :: bx ! rhs vector
!     !real, allocatable, dimension(:) :: by ! rhs vector
!     !real, allocatable, dimension(:) :: bz ! rhs vector
!     integer, allocatable, dimension(:) :: IPIV ! rhs vector
!     integer, allocatable, dimension(:) :: vert_mask
!     character*50 :: dsetname,filename

!     ! Perform mesh-smoothing by minimising an energy functional
!     ! Resultings in a system of the form
!     ! L^n x = 0 where L is a Laplacian matrix
!     ! exponent n controls the fuctional to be minimised, n = 1 is the membrane/Dirichlet energy (minimise stretching)
!     ! n = 2 is the thin-plate energy (minimise bending)

!     ! subject to Dirichlet boundary conditions x = a_anchor for specified anchored vertices

!     ! See for example, Baerentzen et al. , §9.6.2

!     ! if (ismaster) then

!     !     filename = 'continuation/isGhostVert.h5'
!     !     dsetname = trim('isGhostVert')
!     !     call HdfWriteSerialInt2D(filename,dsetname,nv,1,isGhostVert)

!     !     filename = 'continuation/isGhostEdge.h5'
!     !     dsetname = trim('isGhostEdge')
!     !     call HdfWriteSerialInt2D(filename,dsetname,ne,1,isGhostEdge)

!     !     filename = 'continuation/isGhostFace.h5'
!     !     dsetname = trim('isGhostFace')
!     !     call HdfWriteSerialInt2D(filename,dsetname,nf,1,isGhostFace)

!     !     filename = 'continuation/anchorVert.h5'
!     !     dsetname = trim('anchorVert')
!     !     call HdfWriteSerialInt2D(filename,dsetname,nv,1,anchorVert)

!     !     filename = 'continuation/vert_of_edge.h5'
!     !     dsetname = trim('vert_of_edge')
!     !     call HdfWriteSerialInt2D(filename,dsetname,2,ne,vert_of_edge)

!     !     filename = 'continuation/vert_of_face.h5'
!     !     dsetname = trim('vert_of_face')
!     !     call HdfWriteSerialInt2D(filename,dsetname,3,nf,vert_of_face)

!     !     filename = 'continuation/xyz.h5'
!     !     dsetname = trim('xyz')
!     !     call HdfWriteSerialReal2D(filename,dsetname,3,nv,xyz)
!     !  endif


!     !No. of active vertices (non-ghosted)
!     nv_active = count(isGhostVert .eqv. .false.)

!     !Un-anchored should not count ghots by construction, but for safety
!     nua = count( ( .not. isGhostVert ) .and. ( .not. anchorVert ) ) 

!     !write(*,*) "No. of un-anchored vertices: ", count(anchorVert .eqv. .false.)


!     vert_mask = [(i, i = 1, nv)]
!     vert_mask = pack(vert_mask(:) , .not. isGhostVert  ) ! Array of all active vertices
!     ! vert_mask(i) stores the global vertex index of the i'th active (un-ghosted) vertex
    

!     ! Array of active mertices size should match nv_active
!     !if ( shape(vert_mask)  .ne. nv_active ) then
!     !    write(*,*) "Size of vert_mask does not match active vertices, exiting"
!     !    !stop
!     !    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
!     !    call MPI_Finalize(ierr)
!     !endif

!     allocate( L(nv_active,nv_active) )
!     allocate( buffer( nv_active,3) )

!     allocate( Lreduced(nua,nua) )
!     allocate( b( nua,3) )
!     allocate( IPIV(nua) )


!     ! Build mask

!     ! Build Laplacian matrix
!     L(:,:) = 0.0d0
!     buffer(:,:) = 0.0d0




!     do i = 1,ne
!         if ( (isGhostEdge(i) .eqv. .false.)  ) then
!             v1 = vert_of_edge(1,i)
!             v2 = vert_of_edge(2,i)

!             if ( (isGhostVert(v1) .eqv. .true.) .or. (isGhostVert(v2) .eqv. .true.) ) then
!                 write(*,*) "there is a problem"
!             endif

!             ! Mapping from global vert indices to active vert indices
!             v1m = findloc(vert_mask,v1 , 1 )
!             v2m = findloc(vert_mask,v2 , 1 )

!             if ( (v1m .gt. nv_active) .or. (v2m .gt. nv_active) ) then
!                 write(*,*) "there is a problem with v1m v2m"
!             endif

!             L(v1m,v2m) = 1.0d0
!             L(v2m,v1m) = 1.0d0


!             L(v1m,v1m) = 1.0d0
!             L(v2m,v2m) = 1.0d0

!         endif
!     enddo

!     ! Normalise rows by valence
!     ! vert_mask(i) stores the global vertex index of the i'th active (un-ghosted) vertex
!     do i =1,nv_active
!         valence = sum( L(i,:) ) - 1.0d0 ! minus 1, don't count self (from identity matrix)
!         L(i,:) = -L(i,:) / valence
!         L(i,i) = 1.0d0
!     enddo

!     ! Exponentiate for higher-order smoothing
!     L = matmul(L,L)

!     ! Apply Dirichlet boundary conditions to anchored vertices
!     do i =1,nv_active
!         if ( anchorVert( vert_mask(i) ) .eqv. .true. ) then
!             L(i,:) = 0.0d0
!             L(i,i) = 1.0d0
!         endif
!     enddo

!     Lreduced(:,:) = 0.0d0
!     b(:,:) = 0.0d0

!     icnt = 1
!     jcnt = 1

!     do i = 1,nv_active
!         if ( anchorVert( vert_mask(i) ) .eqv. .true.  ) then
!             buffer(:,1) = buffer(:,1) - xyz( 1,vert_mask(i) ) * L(:,i)
!             buffer(:,2) = buffer(:,2) - xyz( 2,vert_mask(i) ) * L(:,i)
!             buffer(:,3) = buffer(:,3) - xyz( 3,vert_mask(i) ) * L(:,i)
!         else ! Un-anchored, append to system
!             do j = 1,nv_active
!                 v1 = vert_mask(j)
!                 if ( anchorVert( vert_mask(v1) ) .eqv. .false.  ) then
!                     Lreduced(icnt,jcnt) = L(i,j)
!                     jcnt = jcnt + 1
!                 endif
!             enddo !end j
!             jcnt = 1
!             icnt = icnt + 1
!         endif
!     enddo ! end i

!     ! Now reduce the RHS vector from buffer
!     icnt = 1

!     do i = 1,nv_active
!         if ( anchorVert( vert_mask(i) ) .eqv. .false.  ) then
!             b(icnt,:) = b(icnt,:) + buffer(i,:)
!             icnt = icnt + 1
!         endif
!     enddo



!     if (ismaster) then
!          filename = 'continuation/Lmat.h5'
!          dsetname = trim('Lmat')
!          call HdfWriteSerialReal2D(filename,dsetname,nua,nua,Lreduced)
!       endif

!     ! Solve
!     call  dgesv	(nua,3,Lreduced,nua,IPIV,b,nua,INFO)


!     ! Now solve linear system
!     ! Test for now, can make much more efficient later


!     !call  dgesv	(nv_active,3,L,nv_active,IPIV,b,nv_active,INFO)

!     ! if (ismaster) then
!     !     write(*,*) "dgesv info flag:", INFO
!     ! endif

!     ! if (ismaster) then
!     !     filename = 'continuation/by.h5'
!     !     dsetname = trim('by')
!     !     call HdfWriteSerialReal2D(filename,dsetname,1,nv_active,by)

!     !     filename = 'continuation/bz.h5'
!     !     dsetname = trim('bz')
!     !     call HdfWriteSerialReal2D(filename,dsetname,1,nv_active,bz)
!     !  endif

!     ! Allocate to update
!     vert_mask = [(i, i = 1, nv)]
!     vert_mask = pack(vert_mask(:) , (.not. isGhostVert) .and. (.not. anchorVert)  ) ! Array of all active vertices
!     write(*,*) "shape(vert_mask)", shape(vert_mask)

!     do i = 1,nua
!         xyz( 1:3 , vert_mask(i) ) = b(i,1:3)
!     enddo


!     ! Reset anchor flags
!     anchorVert(:) = .true.

!     deallocate( L )
!     deallocate( buffer  )

!     deallocate( Lreduced )
!     deallocate( b  )

!     deallocate( IPIV )
!     deallocate( vert_mask )

!     write(*,*) "Reached end of smoothing"
!      call write_tecplot_geom

!     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
!     call MPI_Finalize(ierr)


! end subroutine main_smooth

!------------------------------------------------------
!!!!!!!!!!! LARGE SYSTEM SOLVE VERSION !!!!!!!!!!!!!!!!!!!!!
subroutine main_smooth(nv,ne,nf,xyz,isGhostVert,isGhostEdge,isGhostFace,anchorVert,vert_of_edge,vert_of_face)

    use mpih
    use param, only: ismaster

    implicit none
    integer :: nv,ne,nf
    integer :: nv_active, cnt, v1, v2, v1m, v2m
    integer :: i
    logical, dimension(ne) :: isGhostEdge
    logical, dimension(nf) :: isGhostFace
    logical, dimension(nv) :: isGhostVert, anchorVert
    integer, dimension(2,ne) :: vert_of_edge
    integer, dimension(3,nf) :: vert_of_face
    !integer, dimension(nv) :: vert_mask
    real, dimension(3,nv) :: xyz
    real :: valence
    integer :: INFO
    real, allocatable, dimension(:,:) :: L ! Laplacian matrix
    real, allocatable, dimension(:,:) :: L_squared !
    real, allocatable, dimension(:,:) :: b ! rhs vector
    !real, allocatable, dimension(:) :: bx ! rhs vector
    !real, allocatable, dimension(:) :: by ! rhs vector
    !real, allocatable, dimension(:) :: bz ! rhs vector
    integer, allocatable, dimension(:) :: IPIV ! rhs vector
    integer, allocatable, dimension(:) :: vert_mask
    character*50 :: dsetname,filename

    ! Perform mesh-smoothing by minimising an energy functional
    ! Resultings in a system of the form
    ! L^n x = 0 where L is a Laplacian matrix
    ! exponent n controls the fuctional to be minimised, n = 1 is the membrane/Dirichlet energy (minimise stretching)
    ! n = 2 is the thin-plate energy (minimise bending)

    ! subject to Dirichlet boundary conditions x = a_anchor for specified anchored vertices

    ! See for example, Baerentzen et al. , §9.6.2

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

    !     filename = 'continuation/vert_of_edge.h5'
    !     dsetname = trim('vert_of_edge')
    !     call HdfWriteSerialInt2D(filename,dsetname,2,ne,vert_of_edge)

    !     filename = 'continuation/vert_of_face.h5'
    !     dsetname = trim('vert_of_face')
    !     call HdfWriteSerialInt2D(filename,dsetname,3,nf,vert_of_face)

    !     filename = 'continuation/xyz.h5'
    !     dsetname = trim('xyz')
    !     call HdfWriteSerialReal2D(filename,dsetname,3,nv,xyz)
    !  endif


    !No. of active vertices (non-ghosted)
    nv_active = count(isGhostVert .eqv. .false.)

    !write(*,*) "No. of un-anchored vertices: ", count(anchorVert .eqv. .false.)


    vert_mask = [(i, i = 1, nv)]
    vert_mask = pack(vert_mask(:) , .not. isGhostVert  ) ! Array of all active vertices
    ! vert_mask(i) stores the global vertex index of the i'th active (un-ghosted) vertex
    

    ! Array of active mertices size should match nv_active
    !if ( shape(vert_mask)  .ne. nv_active ) then
    !    write(*,*) "Size of vert_mask does not match active vertices, exiting"
    !    !stop
    !    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    !    call MPI_Finalize(ierr)
    !endif

    allocate( L(nv_active,nv_active) )
    allocate( L_squared(nv_active,nv_active) )
    !allocate( bx(nv_active) )
    !allocate( by(nv_active) )
    !allocate( bz(nv_active) )
    allocate( b( nv_active,3) )
    allocate( IPIV(nv_active) )

    ! Build mask

    ! Build Laplacian matrix
    L(:,:) = 0.0d0
    !bx(:) = 0.0d0
    !by(:) = 0.0d0
    !bz(:) = 0.0d0
    b(:,:) = 0.0d0
    !if (ismaster) then
        !write(*,*) "nv_active is ", nv_active
        !write(*,*) "shape(vert_mask) is ", shape(vert_mask)

    !endif


    do i = 1,ne
        if ( (isGhostEdge(i) .eqv. .false.)  ) then
            v1 = vert_of_edge(1,i)
            v2 = vert_of_edge(2,i)

            if ( (isGhostVert(v1) .eqv. .true.) .or. (isGhostVert(v2) .eqv. .true.) ) then
                write(*,*) "there is a problem"
            endif

            ! Mapping from global vert indices to active vert indices
            v1m = findloc(vert_mask,v1 , 1 )
            v2m = findloc(vert_mask,v2 , 1 )

            if ( (v1m .gt. nv_active) .or. (v2m .gt. nv_active) ) then
                write(*,*) "there is a problem with v1m v2m"
            endif

            L(v1m,v2m) = 1.0d0
            L(v2m,v1m) = 1.0d0


            L(v1m,v1m) = 1.0d0
            L(v2m,v2m) = 1.0d0

        endif
    enddo

    ! Normalise rows by valence
    ! vert_mask(i) stores the global vertex index of the i'th active (un-ghosted) vertex
    do i =1,nv_active
        valence = sum( L(i,:) ) - 1.0d0 ! minus 1, don't count self (from identity matrix)
        L(i,:) = -L(i,:) / valence
        L(i,i) = 1.0d0
    enddo

    ! Exponentiate for higher-order smoothing
    !L = matmul(L,L)
    call DGEMM('N','N',nv_active,nv_active,nv_active,1.0d0,L,nv_active,L,nv_active, 0.0d0,L_squared,nv_active)

    ! Apply Dirichlet boundary conditions to anchored vertices
    do i =1,nv_active
        if ( anchorVert( vert_mask(i) ) .eqv. .true. ) then
            L_squared(i,:) = 0.0d0
            L_squared(i,i) = 1.0d0
            !bx(i) = xyz(1,vert_mask(i))
            !by(i) = xyz(2,vert_mask(i))
            !bz(i) = xyz(3,vert_mask(i))
            b(i,1:3) = xyz(1:3,vert_mask(i))
        endif
    enddo


    ! if (ismaster) then
    !      filename = 'continuation/Lmat.h5'
    !      dsetname = trim('Lmat')
    !      call HdfWriteSerialReal2D(filename,dsetname,nv_active,nv_active,L)
    !   endif


    ! Now solve linear system
    ! Test for now, can make much more efficient later
    call  dgesv	(nv_active,3,L_squared,nv_active,IPIV,b,nv_active,INFO)

    ! if (ismaster) then
    !     write(*,*) "dgesv info flag:", INFO
    ! endif

    ! if (ismaster) then
    !     filename = 'continuation/by.h5'
    !     dsetname = trim('by')
    !     call HdfWriteSerialReal2D(filename,dsetname,1,nv_active,by)

    !     filename = 'continuation/bz.h5'
    !     dsetname = trim('bz')
    !     call HdfWriteSerialReal2D(filename,dsetname,1,nv_active,bz)
    !  endif

    ! Allocate to update

    do i = 1,nv_active
        xyz( 1:3 , vert_mask(i) ) = b(i,1:3)
    enddo

    !if (ismaster) then
    !    write(*,*) "Reached end of smoothing"
    !    write(*,*) "max(b(:,1,2,3))", maxval(b(:,1),1), maxval(b(:,2),1), maxval(b(:,3),1) 
    !    write(*,*) "min(b(:,1,2,3))", minval(b(:,1),1), minval(b(:,2),1), minval(b(:,3),1) 
    !endif

    ! Reset anchor flags
    anchorVert(:) = .true.

    deallocate( L )
    deallocate( L_squared )
    deallocate( b  )
    deallocate( IPIV )
    deallocate( vert_mask )


    ! call write_tecplot_geom

    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    !call MPI_Finalize(ierr)


end subroutine main_smooth

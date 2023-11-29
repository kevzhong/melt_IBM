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

    do mstep=1,1 !KZ: iteration(?) over enforcing immersed-boundary condition
        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)
        call update_both_ghosts(n1,n2,temp,kstart,kend)


        for_xc = 0.0d0
        for_yc = 0.0d0
        for_zc = 0.0d0
        for_temp = 0.0d0

        call mlsForce

        ! if (imelt .eq. 1) then 
        !     call findProbeIndices
        !     ! Calculate dT/dn at immersed interface location (vertices), stored in dtdn_o, dtdn_i for outward/inward
        !     call mls_normDerivs
        !     ! TODO: optional step: dTdn at triangle centroids
        ! endif

        call velforce
        call tempforce
    end do
endif
 
call update_both_ghosts(n1,n2,vx,kstart,kend)
call update_both_ghosts(n1,n2,vy,kstart,kend)
call update_both_ghosts(n1,n2,vz,kstart,kend)
call update_both_ghosts(n1,n2,temp,kstart,kend)

if (imelt .eq. 1) then
    ! Calculate dT/dn at immersed interface location (vertices), stored in dtdn_o, dtdn_i for outward/inward
    call findProbeIndices
    call mls_normDerivs
    ! TODO: optional step: dTdn at triangle centroids
    
    call apply_StefanCondition ! Update vertex locations xyzv
    
    ! Update triangulated geometry details
    ! KZ: Entirety of same computation done by each process since each process stores all the geo info, could be parallelised later if a bottleneck
    do inp = 1,Nparticle
        call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face,maxnf,maxnv)
        call calculate_area(Surface,maxnv,maxnf,xyzv(1:3,:,inp),vert_of_face,sur(:,inp))

        !if (ismaster) then
        !    do i=1,maxnf
        !    write(*,*) "Updated area of triangle ", i, "is ", sur(i,1)
        !    enddo
        !endif
        ! Update Eulerian < -- > Lagrangian forcing transfer coefficient
          cfac = ( sur(:,inp) * h_eulerian ) / celvol ! Note the hard-coded single-particle for cfac
          
          call calculate_normal(tri_nor(:,:,inp),maxnv,maxnf,xyzv(:,:,inp), vert_of_face(:,:))
          call calculate_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,VERTBUFFER,maxnf,faces_of_vert)
    enddo
endif

if (imlsstr.eq.1) then
    call update_part_pos
endif

! Re-meshing check somewhere in this place

 end

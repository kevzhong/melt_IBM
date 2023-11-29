Subroutine particle
use param
use mls_param
use mls_local
use local_arrays
use mpi_param
use mpih

implicit none 

integer :: mstep, inp
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

        if (imelt .eq. 1) then 
            call findProbeIndices
            ! Calculate dT/dn at immersed interface location (vertices), stored in dtdn_o, dtdn_i for outward/inward
            call mls_normDerivs
            ! KZ consider applying the allreduce() to rhs of Stefan condition rather than dtdn?
            !call MPI_ALLREDUCE(MPI_IN_PLACE,dtdn_o,maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
            !call MPI_ALLREDUCE(MPI_IN_PLACE,dtdn_i,maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
            ! TODO: optional step: dTdn at triangle centroids
        endif

        call velforce
        call tempforce
    end do
endif
 
call update_both_ghosts(n1,n2,vx,kstart,kend)
call update_both_ghosts(n1,n2,vy,kstart,kend)
call update_both_ghosts(n1,n2,vz,kstart,kend)
call update_both_ghosts(n1,n2,temp,kstart,kend)

if (imelt .eq. 1) then
    call apply_StefanCondition
    if (ismaster) then
        do inp=1,maxnv
        write(*,*) "Updated vertex location ", inp, "is ", xyzv(:,inp,1)
        enddo
    endif
    ! Update centroids
    ! Update triAreas -> update cfac
    ! Update faceNormals
    ! Update vertexNormals

endif

if (imlsstr.eq.1) then
    call update_part_pos
endif

! Re-meshing check somewhere in this place

 end

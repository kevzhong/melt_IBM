Subroutine particle
use param
use mls_param
use mls_local
use local_arrays
use mpi_param
use mpih

! Loose / explicit coupling approach is adopted (Hu et al. 2001, Borazjani et al. 2008)
! That is, the melting, rigid-body motion, and re-meshing are all de-coupled in time
!
! RK3 substep l:  t_l --------------------------------------->  t_{l+1}
!                       (1)      (2)       (3)        (4)
!
! (1) Apply direct force, enforces T = Tmelt, u = u_object on immersed boundary, accumulate forces/torques
! (2) Melt object from heat-fluxes evaluated at (1)
! (3) Move/rotate body with Newton--Euler equations
! (4) Remesh body if needed

implicit none 

integer :: mstep, inp, i,ntr
integer :: my_up,my_down
logical :: did_remesh
real :: vol_pre, vol_melt, vol_coarse, vol_smooth , drift ! Store some volumes to track drift in volume
integer :: n_ecol,  n_erel

did_remesh = .false.
vol_pre = 0.0d0
vol_melt = 0.0d0
vol_coarse = 0.0d0
vol_smooth = 0.0d0


!------------------------------------ (1) BEGIN IBM FORCING ----------------------------------------------
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
endif

call update_both_ghosts(n1,n2,vx,kstart,kend)
call update_both_ghosts(n1,n2,vy,kstart,kend)
call update_both_ghosts(n1,n2,vz,kstart,kend)
call update_both_ghosts(n1,n2,temp,kstart,kend)

 !------------------------------------ (1) END IBM FORCING ----------------------------------------------


!------------------------------------ (2) BEGIN MELTING  ------------------------------------------------

if (imelt .eq. 1) then
    call apply_StefanCondition ! Update vertex locations xyzv

    do inp = 1,Nparticle
        call calculate_volume (Volume(inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
        ! Track volume change for volume-conserving remeshing
        vol_melt = Volume(1)
        !write(*,*) "Volume after melting:", Volume(1)
    enddo

endif
!------------------------------------ (2) END MELTING  --------------------------------------------------


!------------------------- (3) BEGIN NEWTON--EULER OBJECT MOTION  ---------------------------------------
if (imlsstr.eq.1) then

    ! Update: tri-centroid locations, object COM, Volume, Inertia tensor components (rotation matrix)
    do inp = 1,Nparticle
        call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) 
        call update_tri_normal (tri_nor(:,:,inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
        call calc_rigidBody_params(pos_CM(:,inp),Volume(inp),InertTensor(:,:,inp),maxnv,maxnf,&
        xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp) )
    enddo

    !if (ismaster) then
    !    write(*,*) "Volume fraction is ", Volume(1)*100.0
    !    write(*,*) "trace(I) is ", InertTensor(1,1,1) + InertTensor(2,2,1) + InertTensor(3,3,1)
    !endif    
    !KZ: check if anything else needs to be updated

    ! Move / rotate object by solving Newton--Euler
    call update_part_pos
endif
!------------------------- (3) END NEWTON--EULER OBJECT MOTION  -----------------------------------------

!-------------------------------- (4)  BEGIN REMESHING  -------------------------------------------------
if (iremesh .eq. 1 ) then
do inp = 1,Nparticle
    ! First, evaluate triangle areas and skewnesses -> scan for whether remeshing is needed
    call calculate_area(Surface(inp),maxnv,maxnf,xyzv(1:3,:,inp),vert_of_face(:,:,inp),sur(:,inp),&
                        isGhostFace(:,inp)) 
    call calculate_eLengths(eLengths(:,inp),maxnv,maxne,xyzv(1:3,:,inp), vert_of_edge(:,:,inp),isGhostEdge(:,inp),&
            rm_flag(inp),E_thresh)
    call calculate_skewness (maxne,maxnf,edge_of_face(:,:,inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),isGhostFace(:,inp))


     if ( (rm_flag(inp) .eqv. .true.) ) then

        did_remesh = .true.
        !----------Mesh coarsening----------
        call update_tri_normal (tri_nor(:,:,inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
        call remesh_coarsen (n_ecol,Surface(inp),sur(:,inp),eLengths(:,inp),maxnf,maxne,maxnv,&
        xyzv(:,:,inp),tri_nor(:,:,inp),&
        E_thresh,vert_of_face(:,:,inp),edge_of_face(:,:,inp),vert_of_edge(:,:,inp),&
        face_of_edge(:,:,inp),isGhostFace(:,inp),isGhostEdge(:,inp),isGhostVert(:,inp),rm_flag(inp),&
        anchorVert(:,inp),flagged_edge(:,inp) )

        ! call main_remesh (n_ecol,Surface(inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),maxnf,maxne,maxnv,&
        !                  xyzv(:,:,inp),tri_nor(:,:,inp),&
        !                  A_thresh,skew_thresh,vert_of_face(:,:,inp),edge_of_face(:,:,inp),vert_of_edge(:,:,inp),&
        !                  face_of_edge(:,:,inp),isGhostFace(:,inp),isGhostEdge(:,inp),isGhostVert(:,inp),rm_flag(inp),&
        !                  anchorVert(:,inp),flagged_edge(:,inp) )
        
        ! Update volume after mesh-coarsening to specify target volume change in smoothing
        call calculate_volume (Volume(inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
        vol_coarse = Volume(1)
    
        !----------Mesh smoothing----------
        call remesh_smooth( -(vol_coarse - vol_melt),n_erel,drift,maxnv,maxne,maxnf,xyzv(:,:,inp),isGhostVert(:,inp),&
        isGhostEdge(:,inp),isGhostFace(:,inp),flagged_edge(:,inp),vert_of_edge(:,:,inp), vert_of_face(:,:,inp),&
        face_of_edge(:,:,inp), edge_of_face(:,:,inp) ) 
        
        call calculate_volume (Volume(inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
        vol_smooth = Volume(1)
    endif
    
enddo
endif

!-------------------------------- (4)  END REMESHING  --------------------------------------------------

!--------------------------  PRIMITIVE GEOMETRY INFO  --------------------------------------------------
do inp = 1,Nparticle

    call calculate_area(Surface(inp),maxnv,maxnf,xyzv(1:3,:,inp),vert_of_face(:,:,inp),sur(:,inp),&
                        isGhostFace(:,inp)) ! Update sur
    call calculate_eLengths(eLengths(:,inp),maxnv,maxne,xyzv(1:3,:,inp), vert_of_edge(:,:,inp),isGhostEdge(:,inp),&
                rm_flag(inp),E_thresh)
    call update_tri_normal (tri_nor(:,:,inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
    call calculate_skewness (maxne,maxnf,edge_of_face(:,:,inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),isGhostFace(:,inp))
    call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) ! Update tri_bar
    call calculate_vert_area (Avert(:,inp),maxnv,maxnf,vert_of_face(:,:,inp),sur(:,inp),isGhostFace(:,inp)) ! Update vertex areas
    call calculate_volume (Volume(inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
    call calculate_areaWeighted_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,maxnf,sur(:,inp),vert_of_face(:,:,inp),&
                                    isGhostFace(:,inp), isGhostVert(:,inp) )    


    call calc_rigidBody_params(pos_CM(:,inp),Volume(inp),InertTensor(:,:,inp),maxnv,maxnf,&
    xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp) )
enddo

! Update Eulerian < -- > Lagrangian forcing transfer coefficient
if (imelt .eq. 1) then
    ! Tri area only changes if melting
    cfac(:,:) = ( sur(:,:) * h_eulerian ) / celvol ! Note the hard-coded single-particle for cfac
endif

!-------------------------- END PRIMITIVE GEOMETRY INFO  ---------------------------------------------


if (did_remesh .eqv. .true.) then
   call writeRemeshStats(n_ecol, n_erel, vol_melt - vol_smooth, drift)
    !if (ismaster) then
    !  write(*,*) "Re-meshing residual (Vmelt - Vsmooth)", vol_melt - vol_smooth , "Max vert_drift / dx = ", drift * dx1
    !endif
endif


!---------------- EXIT CONDITION FOR SMALL GEOMETRY ------------------------------------------
if (Volume(1) .lt. V_thresh ) then
    write(*,*) "Geometry V/VE is ", Volume(1)/celvol, "smaller than threshold, exiting now!"
    call write_tecplot_geom
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    call MPI_Finalize(ierr)
endif
!--------------------------------------------------------------------------------------------


! if (imelt .eq. 1) then
    
!     call apply_StefanCondition ! Update vertex locations xyzv
    
!     ! Update triangulated geometry details
!     ! KZ: Entirety of same computation done by each process since each process stores all the geo info, could be parallelised later if a bottleneck
!     do inp = 1,Nparticle
!         vol_pre = Volume(1) ! Track volume drift

!     call calculate_area(Surface(inp),maxnv,maxnf,xyzv(1:3,:,inp),vert_of_face(:,:,inp),sur(:,inp),&
!                         isGhostFace(:,inp),rm_flag(inp),A_thresh) ! Update sur
!     call calculate_eLengths(eLengths(:,inp),maxnv,maxne,xyzv(1:3,:,inp), vert_of_edge(:,:,inp),isGhostEdge(:,inp),&
!                             rm_flag(inp))
!     call update_tri_normal (tri_nor(:,:,inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
!     !call calculate_areaWeighted_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,maxnf,sur(:,inp),vert_of_face(:,:,inp),&
!     !isGhostFace(:,inp), isGhostVert(:,inp) )
!     call calculate_skewness (maxne,maxnf,edge_of_face(:,:,inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),isGhostFace(:,inp),&
!     rm_flag(inp), skew_thresh )

!     !call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) ! Update tri_bar
!     !call calculate_volume2 (Volume(inp),maxnf,tri_nor(:,:,inp),sur(:,inp),tri_bar(:,:,inp),isGhostFace(:,inp))
!     call calculate_volume (Volume(inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
!     ! Volume change due to melting
!     vol_melt = Volume(1)

!     call calculate_area(Surface(inp),maxnv,maxnf,xyzv(1:3,:,inp),vert_of_face(:,:,inp),sur(:,inp),&
!         isGhostFace(:,inp),rm_flag(inp),A_thresh) ! Update sur
!     call calculate_eLengths(eLengths(:,inp),maxnv,maxne,xyzv(1:3,:,inp), vert_of_edge(:,:,inp),isGhostEdge(:,inp),&
!                 rm_flag(inp))
!     call update_tri_normal (tri_nor(:,:,inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
!     call calculate_areaWeighted_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,maxnf,sur(:,inp),vert_of_face(:,:,inp),&
!         isGhostFace(:,inp), isGhostVert(:,inp) )
!     call calculate_skewness (maxne,maxnf,edge_of_face(:,:,inp),sur(:,inp),eLengths(:,inp),skewness(:,inp),isGhostFace(:,inp),&
!         rm_flag(inp), skew_thresh )
!     call calc_centroids_from_vert(tri_bar(1:3,:,inp),xyzv(1:3,:,inp),vert_of_face(:,:,inp),maxnf,maxnv,isGhostFace(:,inp)) ! Update tri_bar
!     call calculate_vert_area (Avert(:,inp),maxnv,maxnf,vert_of_face(:,:,inp),sur(:,inp),isGhostFace(:,inp)) ! Update vertex areas
!     call calculate_volume2 (Volume(inp),maxnf,tri_nor(:,:,inp),sur(:,inp),tri_bar(:,:,inp),isGhostFace(:,inp))
!     !call calculate_volume (Volume(inp),maxnv,maxnf,xyzv(:,:,inp),vert_of_face(:,:,inp),isGhostFace(:,inp))
!     call calculate_areaWeighted_vert_normal (tri_nor(:,:,inp),vert_nor(:,:,inp),maxnv,maxnf,sur(:,inp),vert_of_face(:,:,inp),&
!                                     isGhostFace(:,inp), isGhostVert(:,inp) )

!     enddo

!     ! Update Eulerian < -- > Lagrangian forcing transfer coefficient
!     cfac(:,:) = ( sur(:,:) * h_eulerian ) / celvol ! Note the hard-coded single-particle for cfac

! endif

! Update allocation of I_inv, etc, somewhere

 end

 ! Tag fluid / solid cells in the Eulerian grid given the bounding-box indices defined by ind(3,2) for particle inp
 ! Bounding box indices are obtained from the get_bbox_inds(ind,inp) routine below
 ! Cells are identified as either fluid (VOF = 1), solid (VOF = 0.0), or interface (VOF between 0 and 1)
 ! The tagging operation (fluid, interface, or solid) is done through ray-tagging
 ! That is, by performing queries on ray-triangle intersections
 ! Ray-triangle intersections are efficiently computed using the Möller & Trumbore (2005) algorithm
 ! As noted in §7.5 of O'Rourke (1998), counting the no. of ray-triangle intersections will tell us if a cell is internal/external (solid/fluid)
 ! provided the geometry is 'watertight' (a closed, triangulated geometry)
 !      Cell is external (fluid): even number of intersections:  
 !      Cell is internal (solid): odd number of intersections:   
 !      Cell is interface: ray-triangle intersection point lies inside cell volume coordinates [x-dx/2, x+dx/2],[y-dy/2, y+dy/2],[z-dz/2, z+dz/2]
 ! Interface treatment for VOF somewhat ambiguous still
 !
 ! Besides Möller & Trumbore (2005) algorithm, some additional resources on ray-tracing are:
 !
 ! - IBM lecture notes (de Tullio et al., 2014) http://people.uniroma2.it/roberto.verzicco/IBnotes.pdf
 ! - "Computational Geometry in C" (O'Rourke, 1988, §7.5)
 ! - https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
 ! - "Essential Mathematics for Games and Interactive Applications" (Van Verth & Bishop, 2016, §12.3.5.2)
 ! - "Real-time Rendering" (Akenine-Möller et al., 2018, §22.8)
 !
 ! ---------------------------------------------------------------------------------------------------------------
 ! ____  _____ _   _  ____ ___ _              _             _                   
 ! |  _ \| ____| \ | |/ ___|_ _| |         ___| |_ _ __ __ _| |_ ___  __ _ _   _ 
 ! | |_) |  _| |  \| | |    | || |   _____/ __| __| '__/ _` | __/ _ \/ _` | | | |
 ! |  __/| |___| |\  | |___ | || |__|_____\__ \ |_| | | (_| | ||  __/ (_| | |_| |
 ! |_|   |_____|_| \_|\____|___|_____|    |___/\__|_|  \__,_|\__\___|\__, |\__, |
 !                                                                   |___/ |___/ 
 !
 ! This version of the routine is a pencil-accelerated decomposition strategy for ray-tagging
 ! Implemeneted out of necessity since ray-tagging poses the main bottleneck for large triangulations unfortunately
 ! The strategy uses x-aligned pencils, decomposing the bounding box in the (y,z) directions with pencils p_row, p_col
 ! where p_row * p_col = numtasks
 ! Since the object can shrink, the bounding box size in y-z can be less than numtasks
 ! in this case, we switch back to the original naiive slab method

subroutine pencilTag(ind, inp)
    use mls_param
    use mpih
    use mpi_param
    use rayAux
    !use local_arrays, only: vx,vy,vz
    use param

    implicit none
    integer :: inp, f,j,k, Njk_TOT
    integer, dimension(3,2) :: ind
    integer :: jbegin, jfinish, kbegin, kfinish, jsize, ksize
    integer :: my_jbegin, my_jfinish, my_kbegin, my_kfinish
    integer :: BUFF_SIZE
    integer, allocatable :: pencil_triList(:,:,:), pencil_triCount(:,:)
    logical :: insidePencil, pencil_accelerated
    real, dimension(3) :: V0, V1, V2
    real :: ymin, ymax, zmin, zmax

    real :: vol_sphere

    VOFx(:,:,:) = 1.
    VOFy(:,:,:) = 1.
    VOFz(:,:,:) = 1.
    VOFp(:,:,:) = 1.
    

    ! Maximum number of pencils allowed based on the bounding-box size
    jbegin = ind(2,1) ; jfinish = ind(2,2)
    kbegin = ind(3,1) ; kfinish = ind(3,2)

    jsize = jfinish - jbegin + 1
    ksize = kfinish - kbegin + 1
    Njk_TOT = jsize * ksize

    
    !if (myid .eq. 0) then
    !    write(*,*) "jsize, ksize: ", jsize, ksize
    !    write(*,*) "jstart:jend, kstart:kend ", jbegin, jfinish, kbegin, kfinish
    !    write(*,*) "numtasks ", numtasks
    !endif

        ! Only do pencil-strategy if Njk_TOT >= numtasks, otherwise just do naiive slab method

        !---------------------------- BEGIN PENCIL-STRATEGY -----------------------------------
        
        !------ Allocate work distribution based on bounding-box size ------
        BUFF_SIZE = n1m ! Max. no. of triangles in a pencil

        ! Global bounding box arrays
        !allocate( pencil_triList(jbegin:jfinish, kbegin:kfinish, 1:BUFF_SIZE) )
        allocate( pencil_triList(1:BUFF_SIZE, jbegin:jfinish, kbegin:kfinish ) )
        allocate( pencil_triCount(jbegin:jfinish, kbegin:kfinish ) )
        pencil_triCount = 0
        pencil_triList = 0

        if ( (Njk_TOT .ge. numtasks ) .and. (jsize .ge. p_col) .and. (ksize .ge. p_row ) ) then ! DO PENCIL STRATEGY
          pencil_accelerated = .true.
        ! Retrieve local indices for each pencil
        call pencil_workAlloc(myid, numtasks, p_row, p_col, my_p_row, my_p_col,&
                                jbegin, jfinish, kbegin, kfinish,jsize,ksize, &
                                my_jbegin, my_jfinish, my_kbegin, my_kfinish)
        else ! Naive serial
          pencil_accelerated = .false.
          my_jbegin = jbegin ; my_jfinish = jfinish
          my_kbegin = kbegin ; my_kfinish = kfinish
        endif


        !-------- Begin pencil-accelerated counting of triangles ------------------
        do j = my_jbegin, my_jfinish
            ymin = ym(j) - 0.5/dx2
            ymax = ym(j) + 0.5/dx2
            do k = my_kbegin,my_kfinish
                zmin = zm(k) - 0.5/dx3
                zmax = zm(k) + 0.5/dx3
                do f =1,maxnf
                    if (.not. isGhostFace(f,inp) ) then
                        V0 = xyzv(1:3,vert_of_face(1,f,inp),inp)
                        V1 = xyzv(1:3,vert_of_face(2,f,inp),inp)
                        V2 = xyzv(1:3,vert_of_face(3,f,inp),inp)
                        !call trianglePencilIntersect(insidePencil,V0,V1,V2, ymin, ymax, zmin, zmax)
                        !insidePencil = intersects_yz_range(V0, V1, V2, ymin, ymax, zmin, zmax)
                        call intersects_yz_range(V0, V1, V2, ymin, ymax, zmin, zmax, insidePencil)
                        ! Store triangle in list if intersected
                        if (insidePencil) then
                            pencil_triCount(j,k) = pencil_triCount(j,k) + 1

                            ! Raise error if the triCount exceeds the buffer size
                            if (pencil_triCount(j,k) .gt. BUFF_SIZE) then
                                write(*,*) "triCount(",j,k,") = ", pencil_triCount(j,k), " has exceeded buffer size, aborting"
                                call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
                                call MPI_Finalize(ierr)
                            endif

                            !pencil_triList(j,k,pencil_triCount(j,k)) = f
                            pencil_triList(pencil_triCount(j,k),j,k) = f
                        endif ! end insidePencil

                    endif ! end ghostFace
                enddo ! end f
            enddo ! end k
        enddo !end j

        if (pencil_accelerated ) then
          ! Collect to global arrays if pencil-accelerated
          call MPI_ALLREDUCE(MPI_IN_PLACE,pencil_triCount,jsize*ksize,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)        
          call MPI_ALLREDUCE(MPI_IN_PLACE,pencil_triList,jsize*ksize*BUFF_SIZE,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)        
        endif

        ! Now perform the ray-tagging operations with the computed triangle lists
        call pencil_tagCells(ind,inp,&
        jbegin, jfinish, kbegin, kfinish, BUFF_SIZE, pencil_triList, pencil_triCount)

        deallocate(pencil_triList, pencil_triCount)
    !---------------------------- END PENCIL-STRATEGY -----------------------------------

  !vol_sphere = sum( 1.0 - VOFp(:,:,kstart:kend) ) * celvol
  !call MPI_ALLREDUCE(MPI_IN_PLACE,vol_sphere,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
  !if (myid .eq. 0) then
  !write(*,*) "Vsphere = ", vol_sphere
  !endif

  !call MPI_Barrier(MPI_COMM_WORLD, ierr)
  !call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  !call MPI_Finalize(ierr)

end subroutine pencilTag


subroutine pencil_workAlloc(myid, numtasks, p_row, p_col, my_p_row, my_p_col,&
                            jbegin, jfinish, kbegin, kfinish,jsize,ksize, &
                            my_jbegin, my_jfinish, my_kbegin, my_kfinish)

    ! Given the global work array of size ( jbegin:jfinish, kbegin:kfinish  )
    ! we divide this work evenly into p_row * p_col processes, each with their own
    ! my_jbegin, my_jfinish etc. work indices
    implicit none
    integer, intent(in) :: myid, numtasks, p_row, p_col, my_p_row, my_p_col
    integer, intent(in) :: jbegin, jfinish, kbegin, kfinish, jsize, ksize
    integer, intent(out) :: my_jbegin, my_jfinish, my_kbegin, my_kfinish

    ! Local
    integer :: j_per_proc, k_per_proc, j_rem, k_rem
    integer :: k_block_size, k_offset, j_block_size, j_offset


    ! First compute local work size
    k_per_proc = ksize / p_row
    j_per_proc = jsize / p_col

    k_rem = mod(ksize, p_row)
    j_rem = mod(jsize, p_col)

    ! Row-work division in k direction
    if (my_p_row < k_rem) then
        k_block_size = k_per_proc + 1
        k_offset = my_p_row * k_block_size
    else
        k_block_size = k_per_proc
        k_offset = my_p_row * k_block_size + k_rem
    end if

    ! Column-work division in j direction
    if (my_p_col < j_rem) then
        j_block_size = j_per_proc + 1
        j_offset = my_p_col * j_block_size
    else
        j_block_size = j_per_proc
        j_offset = my_p_col * j_block_size + j_rem
    end if

    ! Now compute local indices in world/global index coordinates

    ! Compute the local start and end indices for k and i for this processor
    my_kbegin = kbegin + k_offset
    my_kfinish   = my_kbegin + k_block_size - 1

    my_jbegin = jbegin + j_offset
    my_jfinish   = my_jbegin + j_block_size - 1

    !write(*,*) "myid: ", myid, " pencil location: ", '(', my_p_row, my_p_col, ')',&
    !     'indices: (', my_jbegin, my_jfinish, ')', '(', my_kbegin, my_kfinish, ')'

end subroutine pencil_workAlloc

subroutine pencil_tagCells(ind,inp,&
    jbegin, jfinish, kbegin, kfinish, BUFF_SIZE, pencil_triList, pencil_triCount)
    ! Tagging for cell-centred scalar grid
    use mls_param
    use mpih
    use mpi_param, only: kstart,kend
    use local_arrays, only: temp
    use param
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: vof
    real :: alpha
    real :: alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real ,dimension(3) :: Q, C
    !real :: u_imh, u_iph, v_jmh, v_jph, w_kmh, w_kph, h31, h32, h33
    real :: udx1, udx2, udx3
    integer :: kc,km,kp,jm,jc,jp,ic,im,ip

    integer :: jbegin, jfinish, kbegin, kfinish, BUFF_SIZE
    !integer ::  pencil_triList(jbegin:jfinish, kbegin:kfinish, 1:BUFF_SIZE) 
    integer ::  pencil_triList(1:BUFF_SIZE, jbegin:jfinish, kbegin:kfinish) 
    integer ::  pencil_triCount(jbegin:jfinish, kbegin:kfinish ) 

  
    udx1=dx1*0.5d0
    udx2=dx2*0.5d0
    udx3=dx3*0.5d0
  
    ! Q stores the cell-centre coordinates and is the query point
    ! C is the control-point and stores the ray origin

    ! Random control-point well-outside computational domain
    !C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]
  
    do i = ind(1,1),ind(1,2)
      Q(1) = xm(i)
      do j = jbegin, jfinish
          Q(2) = ym(j)
        do kk = kbegin, kfinish
             Q(3) = zm(kk)
  
             ! periodic BC
             k = modulo(kk-1,n3m) + 1
              
             if (k.ge.kstart.and.k.le.kend) then 

              x_gc = pos_cm(1:3,inp)

              !call rayTagQ(vof,C,Q,inp)
              !call rayTagQ_pencil(i,j,kk,vof,Q,inp,pencil_triList(j,kk,1:pencil_triCount(j,kk)),pencil_triCount(j,kk))
              call rayTagQ_pencil(i,j,kk,vof,Q,inp,pencil_triList(1:pencil_triCount(j,kk),j,kk),pencil_triCount(j,kk))

              !if ( (kk .eq. 64) .and. (j .eq. 64)  ) then
              !  write(*,*) "VOF: ", vof
              !endif
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1

               kc = k
               km=kc-1
               kp=kc+1

               jc = jj
               jm=jmv(jj)
               jp=jpv(jj)

               ic = ii
               im=imv(ii)
               ip=ipv(ii)

               if(VOFp(ii,jj,k).lt.1.0)then
               VOFp(ii,jj,k) = VOFp(ii,jj,k)
               else
               VOFp(ii,jj,k) = vof
               end if
               
              !  if (vof .eq. 0.0) then
              !   !write(*,*) "Tagging solid cell",ii,jj,k
              !   solid_mask(ii,jj,k) = .true.

              !   !                d  u T   |          1   [                              ]
              !   !             ----------- |  =     ----- |  uT |      -      uT |       |
              !   !                d   x    |i,j,k     dx  [     i+1/2            i-1/2   ]

              !   ! uT |_{i-1/2}
              !   x_grid(1) = xc(i)
              !   x_grid(2) = ym(j)
              !   x_grid(3) = zm(kk)
              !   r = x_grid - x_GC ! relative distance 

              !   u_imh = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

              !   ! uT |_{i+1/2}
              !   x_grid(1) = xc(i+1)
              !   r = x_grid - x_GC ! relative distance 
              !   u_iph = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

              !   h31=( u_iph*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
              !   -u_imh*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1


              !   !                d  v T   |          1   [                              ]
              !   !             ----------- |  =     ----- |  vT |      -      vT |       |
              !   !                d   y    |i,j,k     dy  [     j+1/2            j-1/2   ] 

              !   ! vT |_{j-1/2}
              !   x_grid(1) = xm(i)
              !   x_grid(2) = yc(j)
              !   x_grid(3) = zm(kk)
              !   r = x_grid - x_GC ! relative distance 
              !   v_jmh = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

              !   ! vT |_{j+1/2}
              !   x_grid(2) = yc(j+1)
              !   r = x_grid - x_GC ! relative distance 
              !   v_jph = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

              !   h32=( v_jph*(temp(ic,jp,kc)+temp(ic,jc,kc)) &
              !   -v_jmh*(temp(ic,jc,kc)+temp(ic,jm,kc)) )*udx2


              !   !                d  w T   |          1   [                              ]
              !   !             ----------- |  =     ----- |  wT |      -      wT |       |
              !   !                d   z    |i,j,k     dz  [     k+1/2            k-1/2   ]

              !   ! wT |_{k-1/2}
              !   x_grid(1) = xm(i)
              !   x_grid(2) = ym(j)
              !   x_grid(3) = zc(kk)
              !   r = x_grid - x_GC ! relative distance 
              !   w_kmh = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

              !   ! wT |_{k+1/2}
              !   x_grid(3) = zc(kk+1)
              !   r = x_grid - x_GC ! relative distance 
              !   w_kph = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

              !   h33=( w_kph*(temp(ic,jc,kp)+temp(ic,jc,kc)) &
              !   -w_kmh*(temp(ic,jc,kc)+temp(ic,jc,km)) )*udx3

              !   d_UsolidT_dxj(ii,jj,k) = (h31+h32+h33)


              !  endif
               
               endif
  
        end do
      end do
    end do
  
  end subroutine pencil_tagCells

  subroutine rayTagQ_pencil(ni,nj,nk,vof,Q,inp,&
    triList, numTri)
    ! Tag the cell Q, defined by its cell-centered coordinates as being either interface, fluid or solid
    use param
    use mls_param
    use rayAux
    implicit none
    real, dimension(3) :: V0, V1, V2, Q,C, intPoint
    real, dimension(2) :: xx,yy,zz
  
    real :: vof
    integer :: i,inp, int_count,nt
    logical :: intersect, insideBox

    integer ::  triList(1:numTri) 
    integer ::  numTri

    integer :: ni, nj, nk
  

    ! The ray-origin (control point) must be parallel to the x-aligned pencils
    C = [-xlen*10.0, Q(2) , Q(3) ]
  
    int_count = 0
      
    do nt =1,numTri
        i = triList(nt)
        if (.not. isGhostFace(i,inp) ) then ! Should always pass this condition anyway
            V0 = xyzv(1:3,vert_of_face(1,i,inp),inp)
            V1 = xyzv(1:3,vert_of_face(2,i,inp),inp)
            V2 = xyzv(1:3,vert_of_face(3,i,inp),inp)
  
        ! If triangle in cell, cell is interface, break out            
        ! Pencil already inside (j,k) by construction, so box-intersection testing only x
          call triangleXIntersect(insideBox,V0,V1,V2, Q(1) - 0.5/dx1, Q(1) + 0.5/dx1 )

          if (insideBox .eqv. .true. ) then
            vof = 0.5
            return
          endif       

  
        call rayTriangle_intersect(intersect,intPoint,C,Q,V0,V1,V2)
        

        ! if  ( (ni .eq. 64) .and. (nj .eq. 64) .and. (nk .eq. 64)  .and. (i .eq. 3677) ) then
        !   write(*,*) "V0: ", V0
        !   write(*,*) "V1: ", V1
        !   write(*,*) "V2: ", V2
        !   write(*,*) "intersect: ", intersect

        !  endif
  
            ! If intersected, check if cell is interface cell
            if (intersect) then
              call pointBoxIntersect(insideBox,intPoint, Q(1) - 0.5/dx1, Q(1) + 0.5/dx1 ,&
                                                         Q(2) - 0.5/dx2, Q(2) + 0.5/dx2 ,&
                                                         Q(3) - 0.5/dx3, Q(3) + 0.5/dx3 )


              if (insideBox .eqv. .true. ) then
                vof = 0.5
                ! Exit, no need to check other triangles
                return
              else ! Intersected: add to counter
                int_count = int_count + 1



              endif ! if interface
            endif !if intersect
        endif !ifGhost
    enddo
  
    ! Even intersections: external, odd intersections: internal
    if (modulo(int_count,2) .eq. 0) then !even = fluid
        vof = 1.0
    else ! odd = solid
        vof = 0.0
    endif
  end subroutine rayTagQ_pencil

  ! subroutine slab_tagCells(ind,inp)
  !   ! Tagging for cell-centred scalar grid
  !   use mls_param
  !   use mpih
  !   use mpi_param, only: kstart,kend
  !   use local_arrays, only: temp
  !   use param
  !   implicit none
  !   real, dimension(3)   :: x_grid
  !   real, dimension(3)   :: r, x_GC
  !   real, dimension(3,2) :: lim
  !   real :: vof
  !   real :: alpha
  !   real :: alpha_q
  !   integer :: inp, i,j,k, ii,jj,kk
  !   integer, dimension(3,2) :: ind
  !   real ,dimension(3) :: Q, C
  !   real :: u_imh, u_iph, v_jmh, v_jph, w_kmh, w_kph, h31, h32, h33
  !   real :: udx1, udx2, udx3
  !   integer :: kc,km,kp,jm,jc,jp,ic,im,ip

  
  !   udx1=dx1*0.5d0
  !   udx2=dx2*0.5d0
  !   udx3=dx3*0.5d0
  
  !   ! Q stores the cell-centre coordinates and is the query point
  !   ! C is the control-point and stores the ray origin

  !   ! Random control-point well-outside computational domain
  !   !C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]
  
  !   do i = ind(1,1),ind(1,2)
  !     Q(1) = xm(i)
  !     do j = ind(2,1),ind(2,2)
  !         Q(2) = ym(j)
  !       do kk = ind(3,1),ind(3,2)
  !            Q(3) = zm(kk)
  
  !            ! periodic BC
  !            k = modulo(kk-1,n3m) + 1
              
  
  !            if (k.ge.kstart.and.k.le.kend) then 

  !             x_gc = pos_cm(1:3,inp)

  !             !call rayTagQ(vof,C,Q,inp)
  !             call rayTagQ_slab(vof,Q,kk,inp)
  
  !              ii = modulo(i-1,n1m) + 1
  !              jj = modulo(j-1,n2m) + 1

  !              kc = k
  !              km=kc-1
  !              kp=kc+1

  !              jc = jj
  !              jm=jmv(jj)
  !              jp=jpv(jj)

  !              ic = ii
  !              im=imv(ii)
  !              ip=ipv(ii)

  !              if(VOFp(ii,jj,k).lt.1.0)then
  !              VOFp(ii,jj,k) = VOFp(ii,jj,k)
  !              else
  !              VOFp(ii,jj,k) = vof
  !              end if
               
  !              if (vof .eq. 0.0) then
  !               !write(*,*) "Tagging solid cell",ii,jj,k
  !               solid_mask(ii,jj,k) = .true.

  !               !                d  u T   |          1   [                              ]
  !               !             ----------- |  =     ----- |  uT |      -      uT |       |
  !               !                d   x    |i,j,k     dx  [     i+1/2            i-1/2   ]

  !               ! uT |_{i-1/2}
  !               x_grid(1) = xc(i)
  !               x_grid(2) = ym(j)
  !               x_grid(3) = zm(kk)
  !               r = x_grid - x_GC ! relative distance 

  !               u_imh = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

  !               ! uT |_{i+1/2}
  !               x_grid(1) = xc(i+1)
  !               r = x_grid - x_GC ! relative distance 
  !               u_iph = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

  !               h31=( u_iph*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
  !               -u_imh*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1


  !               !                d  v T   |          1   [                              ]
  !               !             ----------- |  =     ----- |  vT |      -      vT |       |
  !               !                d   y    |i,j,k     dy  [     j+1/2            j-1/2   ] 

  !               ! vT |_{j-1/2}
  !               x_grid(1) = xm(i)
  !               x_grid(2) = yc(j)
  !               x_grid(3) = zm(kk)
  !               r = x_grid - x_GC ! relative distance 
  !               v_jmh = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

  !               ! vT |_{j+1/2}
  !               x_grid(2) = yc(j+1)
  !               r = x_grid - x_GC ! relative distance 
  !               v_jph = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

  !               h32=( v_jph*(temp(ic,jp,kc)+temp(ic,jc,kc)) &
  !               -v_jmh*(temp(ic,jc,kc)+temp(ic,jm,kc)) )*udx2


  !               !                d  w T   |          1   [                              ]
  !               !             ----------- |  =     ----- |  wT |      -      wT |       |
  !               !                d   z    |i,j,k     dz  [     k+1/2            k-1/2   ]

  !               ! wT |_{k-1/2}
  !               x_grid(1) = xm(i)
  !               x_grid(2) = ym(j)
  !               x_grid(3) = zc(kk)
  !               r = x_grid - x_GC ! relative distance 
  !               w_kmh = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

  !               ! wT |_{k+1/2}
  !               x_grid(3) = zc(kk+1)
  !               r = x_grid - x_GC ! relative distance 
  !               w_kph = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

  !               h33=( w_kph*(temp(ic,jc,kp)+temp(ic,jc,kc)) &
  !               -w_kmh*(temp(ic,jc,kc)+temp(ic,jc,km)) )*udx3

  !               d_UsolidT_dxj(ii,jj,k) = (h31+h32+h33)


  !              endif
               
  !              endif
  
  !       end do
  !     end do
  !   end do
  
  !   call interp_vof
  ! end subroutine slab_tagCells

  ! subroutine rayTagQ_slab(vof,Q,k,inp)
  !   ! Tag the cell Q, defined by its cell-centered coordinates as being either interface, fluid or solid
  !   use param
  !   use mls_param
  !   implicit none
  !   real, dimension(3) :: V0, V1, V2, Q,C, intPoint
  !   real, dimension(2) :: xx,yy,zz
  
  !   real :: vof, zcent
  !   integer :: i,inp, int_count, pad, k1,k
  !   logical :: intersect, insideBox
  
  !   ! Force C to be parallel to the computational grid lines: efficient
  
  !   C = [-xlen*10.0, -ylen*3.0, Q(3) ]
  
  !   int_count = 0
    
  !   pad = 2
  
  !   do i =1,maxnf
  !       if (.not. isGhostFace(i,inp) ) then
  !           V0 = xyzv(1:3,vert_of_face(1,i,inp),inp)
  !           V1 = xyzv(1:3,vert_of_face(2,i,inp),inp)
  !           V2 = xyzv(1:3,vert_of_face(3,i,inp),inp)
  
  !           ! If triangle in cell, cell is interface, break out
  !           call triangleBoxIntersect(insideBox,V0,V1,V2, Q(1) - 0.5/dx1, Q(1) + 0.5/dx1,&
  !                                                         Q(2) - 0.5/dx2, Q(2) + 0.5/dx2,&
  !                                                         Q(3) - 0.5/dx3, Q(3) + 0.5/dx3)

  !         if (insideBox .eqv. .true. ) then
  !           vof = 0.5
  !           return
  !         endif       
  
  !         zcent = (1.0/3.0)* ( V0(3) + V1(3) + V2(3) )
  !         k1 = floor(zcent*dx3) + 1
  
  !         if (abs(k1 - k) .le. pad ) then
  !           call rayTriangle_intersect(intersect,intPoint,C,Q,V0,V1,V2)
  
  !           ! If intersected, check if cell is interface cell
  !           if (intersect) then
  !               call pointBoxIntersect(insideBox,intPoint, Q(1) - 0.5/dx1, Q(1) + 0.5/dx1 ,&
  !                                                          Q(2) - 0.5/dx2, Q(2) + 0.5/dx2 ,&
  !                                                          Q(3) - 0.5/dx3, Q(3) + 0.5/dx3 )

  !             if (insideBox .eqv. .true. ) then
  !               vof = 0.5
  
  !               ! Exit, no need to check other triangles
  !               return
  
  !             else ! Intersected: add to counter
  !               int_count = int_count + 1
  !             endif ! if interface
  !           endif !if intersect
  !         endif
  !       endif !ifGhost
  !   enddo
  
  !   ! Even intersections: external, odd intersections: internal
  !   if (modulo(int_count,2) .eq. 0) then !even = fluid
  !       vof = 1.0
  !   else ! odd = solid
  !       vof = 0.0
  !   endif
  
  ! end subroutine rayTagQ_slab


  ! subroutine interp_vof
  !   ! Interpolate VOFx, VOFv, VOFw from cell-centered VOFp information
  !   use param
  !   use mpi_param, only: kstart,kend
  !   implicit none
  !   integer :: ic,jc,kc
  !   integer :: km,kp,jm,jp,im,ip
  
  !   !call update_both_ghosts(n1,n2,VOFx,kstart,kend)
  !   !call update_both_ghosts(n1,n2,VOFy,kstart,kend)
  !   !call update_both_ghosts(n1,n2,VOFz,kstart,kend)
  !   call update_both_ghosts(n1,n2,VOFp,kstart,kend)
  
  !   do kc=kstart,kend
  !     km=kc-1
  !     kp=kc+1
  !       do jc=1,n2m
  !         jm=jmv(jc)
  !         jp=jpv(jc)
  !         do ic=1,n1m
  !           ip=ipv(ic)
  !           im=imv(ic)
            
  !           VOFx(ic,jc,kc) = 0.5 * (VOFp(im,jc,kc) + VOFp(ic,jc,kc) )
  !           VOFy(ic,jc,kc) = 0.5 * (VOFp(ic,jm,kc) + VOFp(ic,jc,kc) )
  !           VOFz(ic,jc,kc) = 0.5 * (VOFp(ic,jc,km) + VOFp(ic,jc,kc) )
  !         enddo
  !       enddo
  !   enddo
  ! end subroutine interp_vof


!   subroutine trianglePencilIntersect(insidePencil,V0,V1,V2, ymin, ymax, zmin, zmax)
!     ! Test if triangle (V0,V1,V2) is inside pencil (j,k)
!     ! Pencil(j,k) is defined spatially by its cell-centred coordinate, XC_jk which stores the (ym,zm) coordinates
!     implicit none
!     real, dimension(3) :: V0, V1, V2
!     real ::  ymin_tri, ymax_tri, zmin_tri, zmax_tri, ymin, ymax, zmin, zmax
!     logical :: insidePencil, point_in_range
  
!     insidePencil = .false.

!     ymin_tri = min( V0(2), V1(2), V2(2) )
!     ymax_tri = max( V0(2), V1(2), V2(2) )

!     zmin_tri = min( V0(3), V1(3), V2(3) )
!     zmax_tri = max( V0(3), V1(3), V2(3) )


!     ! Check if the bounding box of the triangle overlaps the 2D range
!     if (ymax_tri < ymin .or. ymin_tri > ymax) return  ! No overlap in y-direction
!     if (zmax_tri < zmin .or. zmin_tri > zmax) return  ! No overlap in z-direction


!     ! Check if any vertex of the triangle is inside the 2D range
!     point_in_range = .false.
!     point_in_range = point_in_range .or. (V0(2) >= ymin .and. V0(2) <= ymax .and. V0(3) >= zmin .and. V0(3) <= zmax)
!     point_in_range = point_in_range .or. (V1(2) >= ymin .and. V1(2) <= ymax .and. V1(3) >= zmin .and. V1(3) <= zmax)
!     point_in_range = point_in_range .or. (V2(2) >= ymin .and. V2(2) <= ymax .and. V2(3) >= zmin .and. V2(3) <= zmax)

!     if (point_in_range) then
!         insidePencil = .true.
!         return
!     end if


! end subroutine trianglePencilIntersect

! subroutine get_bbox_inds(bbox_inds,inp)
!     ! Retrieve the xyz indices of the bounding box for a given particle
!     use param
!     use mls_param
!     implicit none
!     integer :: i, nf, tri_ind
!     real, dimension(3,2) :: lim
!     integer, dimension(3,2) :: bbox_inds
!     integer  :: inp, padSize
  
!     ! Padding size indices for safety
!     padSize = 2
  
!   ! get bounding box
!     do i = 1,3
!       lim(i,1) = minval( pack(xyzv(i,:,inp) , .not. isGhostVert(:,inp)  ) )
!       lim(i,2) = maxval( pack(xyzv(i,:,inp) , .not. isGhostVert(:,inp)  ) )
!     end do
  
!     bbox_inds = floor(lim*dx1) + 1 ! compute indices cell centered
  
!     ! expanding bounding box to be extra safe
!     bbox_inds(:,1) = bbox_inds(:,1) - padSize
!     bbox_inds(:,2) = bbox_inds(:,2) + padSize
  
!     ! Hard code vertical for testing free-fall
!     !bbox_inds(3,1) = 1
!     !bbox_inds(3,2) = n3m
  
  
!     !! Hard-code the full domain for testing
!     !bbox_inds(:,1) = [1, 1, 1] 
!     !bbox_inds(:,2) = [n1m, n2m, n3m]
  
! end subroutine get_bbox_inds

! subroutine get_periodic_indices(k,x)
!   use param
!   implicit none
!   integer :: k
!   real    :: x(3)

!   if (k .ge. n3) then
!      k = k - n3m
!     x(3) = x(3) - zlen
!   end if

!   if (k .lt. 1) then
!      k = k + n3m
!      x(3) = x(3) + zlen
!   end if
! end subroutine


subroutine apply_StefanCondition
USE param
USE mls_param
use mpih
implicit none
real, dimension(3) :: nhat
integer :: inp,nv,nf
integer :: v1, v2, v3
real :: vmelt_fixed
real :: rand_noise


! Prescribe a fixed local melt-rate for testing (should be negative)
vmelt_fixed = -7.0

! RK3 memory swap
vmelt_m1 = vmelt

do inp=1,Nparticle
    do nv = 1,maxnv
        call random_number(rand_noise) ! different for each rank, be careful
        if (isGhostVert(nv,inp) .eqv. .false. ) then
            nhat(1:3) = vert_nor(1:3,nv,inp)
            !vmelt(1:3,nv,inp) = cpliquid / latHeat * ( qw_iVert(nv,inp) - qw_oVert(nv,inp) ) * nhat(1:3)


            !vmelt(1:3,nv,inp) = (vmelt_fixed - rand_noise * 0.5) * nhat(1:3)
            vmelt(1:3,nv,inp) = vmelt_fixed  * nhat(1:3)

            ! Update interface location
            !xyzv(1:3,nv,inp) = xyzv(1:3,nv,inp) +  al * dt * vmelt(1:3,nv,inp)
            xyzv(1:3,nv,inp) = xyzv(1:3,nv,inp) + dt * ( ga * vmelt(1:3,nv,inp) + ro * vmelt_m1(1:3,nv,inp) )

        endif
    enddo

    !KZ: interface motion does not contribute to fluid velocity
    
    ! ! Add vmelt contribution to vel_tri for each Lagrangian marker
    ! do nf = 1,maxnf
    !     if (.not. isGhostFace(nf,inp) ) then
    !         ! Vertex -> face interpolation, equal weigting
    !         v1 = vert_of_face(1,nf,inp)
    !         v2 = vert_of_face(2,nf,inp)
    !         v3 = vert_of_face(3,nf,inp)

    !         ! Equal load distribution (1/3) for each vertex on each triangle
    !         !vel_tri(1:3,nf,inp) = vel_tri(1:3,nf,inp) + (1.0/3.0) * vmelt(1:3,v1,inp)
    !         !vel_tri(1:3,nf,inp) = vel_tri(1:3,nf,inp) + (1.0/3.0) * vmelt(1:3,v2,inp)
    !         !vel_tri(1:3,nf,inp) = vel_tri(1:3,nf,inp) + (1.0/3.0) * vmelt(1:3,v3,inp)

    !         ! Equal load distribution (1/3) for each vertex on each triangle
    !         vel_tri(1:3,nf,inp) = (1.0/3.0) * vmelt(1:3,v1,inp)
    !         vel_tri(1:3,nf,inp) = (1.0/3.0) * vmelt(1:3,v2,inp)
    !         vel_tri(1:3,nf,inp) = (1.0/3.0) * vmelt(1:3,v3,inp)

    !         !Reset at end

    !     endif
    ! enddo

enddo


end subroutine apply_StefanCondition
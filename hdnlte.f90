    subroutine hdnlte
    use param
    use local_arrays, only: vy,vz,htemp,vx,temp,rhs
    use mpi_param, only: kstart,kend
    use param
    use phasefield
    !use mls_param,only: dens_ratio, omega_c, vel_CM, pos_CM, Nparticle
    implicit none
    integer :: jc,kc
    integer :: km,kp,jm,jp,ic,im,ip
    real    :: h32,h33,h31
    real    :: udx1,udx2,udx3
    real, dimension(3)   :: x_grid, r
    integer :: inp
    real    :: u_imh, v_jmh, w_kmh, u_iph
    integer :: im2


    ! Staggered grid arrangement, the stagger is by -1/2, index convention:

    !                     vy(i,j+1,k)
    !                      ^
    !                      |
    !                      |
    !    __________________|___________________
    !   |                                      |
    !   |                                      |
    !   |                                      |
    !   |                                      |
    !   |  vx(i,j,k)                           |
    ! ----->                               ------>  vx(i+1,j,k)
    !   |                  O                   |
    !   |              temp(i,j,k)             |
    !   |              or p(i,j,k)             |
    !   |                                      |
    !   |                                      |
    !   |                  ^                   |
    !   |__________________|___________________
    !                      |
    !                      |
    !                     vy(i,j,k)


    udx1=dx1*0.5d0
    udx2=dx2*0.5d0
    udx3=dx3*0.5d0

    do kc=kstart,kend
    km=kc-1
    kp=kc+1
    do jc=1,n2m
    jm=jmv(jc)
    jp=jpv(jc)
    do ic=1,n1m
    ip=ipv(ic)
    im=imv(ic)

    !---------------------  FLUID VELOCITY TERMS -------------------------------
!                d  u T   |          1   [                              ]
!             ----------- |  =     ----- |  uT |      -      uT |       |
!                d   x    |i,j,k     dx  [     i+1/2            i-1/2   ]
!
!

         h31=( vx(ip,jc,kc)*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
                -vx(ic,jc,kc)*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1

        !  h31=( (vx(ip,jc,kc)-Usolid)*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
        !         -(vx(ic,jc,kc)-Usolid)*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1
!
!                d  v T   |          1   [                              ]
!             ----------- |  =     ----- |  vT |      -      vT |       |
!                d   y    |i,j,k     dy  [     j+1/2            j-1/2   ]  
!

    h32=( vy(ic,jp,kc)*(temp(ic,jp,kc)+temp(ic,jc,kc)) &
               -vy(ic,jc,kc)*(temp(ic,jc,kc)+temp(ic,jm,kc)) )*udx2


!
!                d  w T   |          1   [                              ]
!             ----------- |  =     ----- |  wT |      -      wT |       |
!                d   z    |i,j,k     dz  [     k+1/2            k-1/2   ]    
!

    h33=( vz(ic,jc,kp)*(temp(ic,jc,kp)+temp(ic,jc,kc)) &
               -vz(ic,jc,kc)*(temp(ic,jc,kc)+temp(ic,jc,km)) )*udx3



    htemp(ic,jc,kc)=-(h31+h32+h33)
    enddo
    enddo
    enddo


    if (meltmode .eq. 1) then
    
    do kc=kstart,kend
    km=kc-1
    kp=kc+1
    do jc=1,n2m
    jm=jmv(jc)
    jp=jpv(jc)
    do ic=1,n1m
    ip=ipv(ic)
    im=imv(ic)
    im2 = 1 + modulo(ic - 1 - 2, n1m)

        ! rhs stores dphi at this point
        htemp(ic,jc,kc) = htemp(ic,jc,kc) + latHeat / cpliquid * rhs(ic,jc,kc) / (al * dt)
        
        ! CDS
        h31 = Usolid * udx1 * ( phi(ip,jc,kc)  -  phi(im,jc,kc) )  

        htemp(ic,jc,kc) = htemp(ic,jc,kc) + latHeat / cpliquid * h31

    enddo
    enddo
    enddo
        
    endif

    return
end subroutine hdnlte


! UPWINDING
subroutine hdnlte_QUICK
    use param
    use local_arrays, only: vy,vz,htemp,vx,temp,temp2,xflux_imh,yflux_jmh,zflux_kmh,rhs
    use mpi_param, only: kstart,kend
    use param
    use phasefield
    implicit none
    integer :: jc,kc
    integer :: km,kp,jm,jp,ic,im,ip
    real    :: h32,h33,h31
    real    :: udx1,udx2,udx3
    real, dimension(3)   :: x_grid, r
    integer :: inp
    real    :: u_imh, v_jmh, w_kmh, u_iph
    integer :: im2,jm2,km2


    ! Memory swap for larger ghost-cell array
    !temp2(1:n1,1:n2,kstart-2:kend+2)
    !temp(1:n1,1:n2,kstart-1:kend+1)
    temp2(:,:,kstart:kend) = temp(:,:,kstart:kend)
    call update_both_ghosts2(n1,n2,temp2,kstart,kend)

    udx1=dx1*0.5d0
    udx2=dx2*0.5d0
    udx3=dx3*0.5d0

    ! First accumulate fluxes at cell faces
    do kc=kstart,kend
    km=kc-1
    kp=kc+1
    km2 = kc-2
    do jc=1,n2m
    jm=jmv(jc)
    jp=jpv(jc)
    jm2 = 1 + modulo(jc - 1 - 2, n2m)

    do ic=1,n1m
    ip=ipv(ic)
    im=imv(ic)
    im2 = 1 + modulo(ic - 1 - 2, n1m)

    ! x-flux
    if (vx(ic,jc,kc) .ge. 0.0) then
        xflux_imh(ic,jc,kc) = vx(ic,jc,kc) * 0.125 * ( -1.0*temp2(im2,jc,kc) + 6.0*temp2(im,jc,kc) + 3.0*temp2(ic,jc,kc) )
    else
        xflux_imh(ic,jc,kc) = vx(ic,jc,kc) * 0.125 * (  3.0*temp2(im,jc,kc) + 6.0*temp2(ic,jc,kc) - 1.0*temp2(ip,jc,kc) )
    endif

    ! y-flux
    if (vy(ic,jc,kc) .ge. 0.0) then
        yflux_jmh(ic,jc,kc) = vy(ic,jc,kc) * 0.125 * ( -1.0*temp2(ic,jm2,kc) + 6.0*temp2(ic,jm,kc) + 3.0*temp2(ic,jc,kc) )
    else
        yflux_jmh(ic,jc,kc) = vy(ic,jc,kc) * 0.125 * (  3.0*temp2(ic,jm,kc) + 6.0*temp2(ic,jc,kc) - 1.0*temp2(ic,jp,kc) )
    endif

    ! z-flux
    if (vz(ic,jc,kc) .ge. 0.0) then
        zflux_kmh(ic,jc,kc) = vz(ic,jc,kc) * 0.125 * ( -1.0*temp2(ic,jc,km2) + 6.0*temp2(ic,jc,km) + 3.0*temp2(ic,jc,kc) )
    else
        zflux_kmh(ic,jc,kc) = vz(ic,jc,kc) * 0.125 * (  3.0*temp2(ic,jc,km) + 6.0*temp2(ic,jc,kc) - 1.0*temp2(ic,jc,kp) )
    endif

    enddo
    enddo
    enddo

    call update_both_ghosts(n1,n2,xflux_imh,kstart,kend)
    call update_both_ghosts(n1,n2,yflux_jmh,kstart,kend)
    call update_both_ghosts(n1,n2,zflux_kmh,kstart,kend)

    do kc=kstart,kend
    km=kc-1
    kp=kc+1
    do jc=1,n2m
    jm=jmv(jc)
    jp=jpv(jc)

    do ic=1,n1m
    ip=ipv(ic)
    im=imv(ic)


    h31 = ( xflux_imh(ip,jc,kc) - xflux_imh(ic,jc,kc) ) * dx1
    h32 = ( yflux_jmh(ic,jp,kc) - yflux_jmh(ic,jc,kc) ) * dx2
    h33 = ( zflux_kmh(ic,jc,kp) - zflux_kmh(ic,jc,kc) ) * dx3

    htemp(ic,jc,kc) = -(h31+h32+h33)
    enddo
    enddo
    enddo

    if (meltmode .eq. 1) then
    
    do kc=kstart,kend
    km=kc-1
    kp=kc+1
    do jc=1,n2m
    jm=jmv(jc)
    jp=jpv(jc)
    do ic=1,n1m
    ip=ipv(ic)
    im=imv(ic)
    im2 = 1 + modulo(ic - 1 - 2, n1m)

        ! rhs stores dphi at this point
        htemp(ic,jc,kc) = htemp(ic,jc,kc) + latHeat / cpliquid * rhs(ic,jc,kc) / (al * dt)
        
        ! CDS
        !h31 = Usolid * udx1 * ( phi(ip,jc,kc)  -  phi(im,jc,kc) )  

        h31 = Usolid * 0.125 * dx1 * ( phi(im2,jc,kc) - 7.0*phi(im,jc,kc) + 3.0*phi(ic,jc,kc) + 3.0*phi(ip,jc,kc) )

        htemp(ic,jc,kc) = htemp(ic,jc,kc) + latHeat / cpliquid * h31

    enddo
    enddo
    enddo
        
    endif

    return
end subroutine hdnlte_QUICK
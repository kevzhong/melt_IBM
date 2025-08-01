    subroutine hdnlte
    use param
    use local_arrays, only: vy,vz,htemp,vx,temp
    use mpi_param, only: kstart,kend
    use mls_param,only: dens_ratio, omega_c, vel_CM, pos_CM, Nparticle
    implicit none
    integer :: jc,kc
    integer :: km,kp,jm,jp,ic,im,ip
    real    :: h32,h33,h31
    real    :: udx1,udx2,udx3
    real, dimension(3)   :: x_grid, r
    integer :: inp
    real    :: u_imh, v_jmh, w_kmh


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

    do inp = 1,Nparticle
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



    htemp(ic,jc,kc)=-(h31+h32+h33)*VOFp(ic,jc,kc)

    !----------------------- SOLID VELOCITY TERMS --------------------------------------

    ! U_i+1/2 = U_i-1/2 since moment arm does not depend on x location

    ! uT |_{i-1/2}
    x_grid(1) = xc(ic)
    x_grid(2) = ym(jc)
    x_grid(3) = zm(kc)
    r = x_grid - pos_cm(1:3,inp) ! relative distance
    ! Periodicity: minimum image convection
    r(1) = r(1) - xlen * nint(r(1) / xlen)
    r(2) = r(2) - ylen * nint(r(2) / ylen)
    r(3) = r(3) - zlen * nint(r(3) / zlen)

    u_imh = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

    ! uT |_{i+1/2}
    !x_grid(1) = xc(i+1)
    !r(1) = x_grid(1) - pos_cm(1,inp) ! relative distance 
    ! Periodicity: minimum image convection
    !r(1) = r(1) - xlen * nint(r(1) / xlen)
    !u_iph = vel_CM(1,inp) + omega_c(2,inp)*r(3) - omega_c(3,inp)*r(2)

    ! h31=( u_iph*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
    ! -u_imh*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1
    h31 = u_imh * ( temp(ip,jc,kc) - temp(im,jc,kc)  )*udx1

    ! vT |_{j-1/2}
    x_grid(1) = xm(ic)
    x_grid(2) = yc(jc)
    !x_grid(3) = zm(kc)
    r = x_grid - pos_cm(1:3,inp) ! relative distance
    ! Periodicity: minimum image convection
    r(1) = r(1) - xlen * nint(r(1) / xlen)
    r(2) = r(2) - ylen * nint(r(2) / ylen)
    r(3) = r(3) - zlen * nint(r(3) / zlen)
    v_jmh = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

    ! vT |_{j+1/2}
    !x_grid(2) = yc(j+1)
    !r = x_grid - x_GC ! relative distance 
    !v_jph = vel_CM(2,inp) + omega_c(3,inp)*r(1) - omega_c(1,inp)*r(3)

    ! h32=( v_jph*(temp(ic,jp,kc)+temp(ic,jc,kc)) &
    ! -v_jmh*(temp(ic,jc,kc)+temp(ic,jm,kc)) )*udx2

    h32 = v_jmh * ( temp(ic,jp,kc) - temp(ic,jm,kc) )*udx2

    ! wT |_{k-1/2}
    !x_grid(1) = xm(ic)
    x_grid(2) = ym(jc)
    x_grid(3) = zc(kc)
    r = x_grid - pos_cm(1:3,inp) ! relative distance
    w_kmh = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

    ! wT |_{k+1/2}
    !x_grid(3) = zc(kk+1)
    !r = x_grid - x_GC ! relative distance 
    !w_kph = vel_CM(3,inp) + omega_c(1,inp)*r(2) - omega_c(2,inp)*r(1)

    h33 = w_kmh * ( temp(ic,jc,kp) - temp(ic,jc,km) )*udx3

    ! h33=( w_kph*(temp(ic,jc,kp)+temp(ic,jc,kc)) &
    ! -w_kmh*(temp(ic,jc,kc)+temp(ic,jc,km)) )*udx3
    htemp(ic,jc,kc)= htemp(ic,jc,kc) -(h31+h32+h33)*(1.0 - VOFp(ic,jc,kc) )

    enddo
    enddo
    enddo
enddo

    return
    end


    ! subroutine hdnlte
!     use param
!     use local_arrays, only: vy,vz,htemp,vx,temp
!     use mpi_param, only: kstart,kend
!     use mls_param,only: dens_ratio
!     implicit none
!     integer :: jc,kc
!     integer :: km,kp,jm,jp,ic,im,ip
!     real    :: h32,h33,h31
!     real    :: udx1,udx2,udx3

!     ! Staggered grid arrangement, the stagger is by -1/2, index convention:

!     !                     vy(i,j+1,k)
!     !                      ^
!     !                      |
!     !                      |
!     !    __________________|___________________
!     !   |                                      |
!     !   |                                      |
!     !   |                                      |
!     !   |                                      |
!     !   |  vx(i,j,k)                           |
!     ! ----->                               ------>  vx(i+1,j,k)
!     !   |                  O                   |
!     !   |              temp(i,j,k)             |
!     !   |              or p(i,j,k)             |
!     !   |                                      |
!     !   |                                      |
!     !   |                  ^                   |
!     !   |__________________|___________________
!     !                      |
!     !                      |
!     !                     vy(i,j,k)


!     udx1=dx1*0.5d0
!     udx2=dx2*0.5d0
!     udx3=dx3*0.5d0

!     do kc=kstart,kend
!     km=kc-1
!     kp=kc+1
!     do jc=1,n2m
!     jm=jmv(jc)
!     jp=jpv(jc)
!     do ic=1,n1m
!     ip=ipv(ic)
!     im=imv(ic)

!     if ( solid_mask(ic,jc,kc) .eqv. .false. ) then
! !                d  u T   |          1   [                              ]
! !             ----------- |  =     ----- |  uT |      -      uT |       |
! !                d   x    |i,j,k     dx  [     i+1/2            i-1/2   ]
! !
! !

!          h31=( vx(ip,jc,kc)*(temp(ip,jc,kc)+temp(ic,jc,kc)) & 
!                 -vx(ic,jc,kc)*(temp(ic,jc,kc)+temp(im,jc,kc)) )*udx1

! !
! !                d  v T   |          1   [                              ]
! !             ----------- |  =     ----- |  vT |      -      vT |       |
! !                d   y    |i,j,k     dy  [     j+1/2            j-1/2   ]  
! !

!     h32=( vy(ic,jp,kc)*(temp(ic,jp,kc)+temp(ic,jc,kc)) &
!                -vy(ic,jc,kc)*(temp(ic,jc,kc)+temp(ic,jm,kc)) )*udx2


! !
! !                d  w T   |          1   [                              ]
! !             ----------- |  =     ----- |  wT |      -      wT |       |
! !                d   z    |i,j,k     dz  [     k+1/2            k-1/2   ]    
! !

!     h33=( vz(ic,jc,kp)*(temp(ic,jc,kp)+temp(ic,jc,kc)) &
!                -vz(ic,jc,kc)*(temp(ic,jc,kc)+temp(ic,jc,km)) )*udx3



!     htemp(ic,jc,kc)=-(h31+h32+h33)

!     else ! solid_mask eqv. true 
!     htemp(ic,jc,kc) = -d_UsolidT_dxj(ic,jc,kc)

!     endif


!     enddo
!     enddo
!     enddo

!     return
!     end

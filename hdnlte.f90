subroutine hdnlte
    use param
    use local_arrays, only: vy,vz,htemp,vx,temp
    use mpi_param, only: kstart,kend
    use mls_param,only: dens_ratio
    implicit none
    integer :: jc,kc
    integer :: km,kp,jm,jp,ic,im,ip
    real    :: h32,h33,h31
    real    :: udx1,udx2,udx3

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

    if ( solid_mask(ic,jc,kc) .eqv. .false. ) then
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



    htemp(ic,jc,kc)=-(h31+h32+h33)

    else ! solid_mask eqv. true 
    htemp(ic,jc,kc) = -d_UsolidT_dxj(ic,jc,kc)

    endif


    enddo
    enddo
    enddo

    return
    end

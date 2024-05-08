      subroutine hdnl3 
      use param
      use local_arrays, only: vy,vz,qcap,vx,temp,forcz
      use mpi_param, only: kstart,kend
      use mls_param,only: dens_ratio
      implicit none
      integer :: jc,kc
      integer :: km,kp,jmm,jpp,ic,im,ip
      real    :: h32,h33,h31
      real    :: udx1,udx2,udx3, fbz, T_interp

      udx1=dx1*0.25
      udx2=dx2*0.25
      udx3=dx3*0.25

      do kc=kstart,kend
      km=kc-1
      kp=kc+1
      do jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
!
!
!    vz vx term
!
!
!                d  q_x q_z 
!             -----------
!                d   x      
!
!
      h31=(((vx(ip,jc,kc)+vx(ip,jc,km)) &
           *(vz(ip,jc,kc)+vz(ic,jc,kc))) &
          -((vx(ic,jc,kc)+vx(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(im,jc,kc))))*udx1
!
!    vz vy term
!
!
!                d  q_y q_z 
!             -----------
!                d   y      
!
      h32=(((vy(ic,jpp,kc)+vy(ic,jpp,km)) &
           *(vz(ic,jpp,kc)+vz(ic,jc,kc))) &
          -((vy(ic,jc,kc)+vy(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(ic,jmm,kc))))*udx2
!
!    vz vz term
!
!
!                 d  q_z q_z 
!                -----------
!                 d   z      
!
      h33=(((vz(ic,jc,kp)+vz(ic,jc,kc)) &
           *(vz(ic,jc,kp)+vz(ic,jc,kc))) &
          -((vz(ic,jc,kc)+vz(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(ic,jc,km))))*udx3


!   buoyancy term
      T_interp =  ( temp(ic,jc,kc) + temp(ic,jc,km) ) * 0.5d0 
      fbz = betagz * ( T_interp - Tliq ) ! Relative to ambient liquid
 

      !qcap(ic,jc,kc)=-(h31+h32+h33) + fbz + &
      !                  az(ic,jc,kc)*forcz(ic,jc,kc)/zlen+   &
      !                (1.0-az(ic,jc,kc))*forcz(ic,jc,kc)*dens_ratio/(zlen)

      ! qcap(ic,jc,kc)=( -(h31+h32+h33) + fbz ) * VOFz(ic,jc,kc) + &
      !                   VOFz(ic,jc,kc)*forcz(ic,jc,kc)/zlen+   &
      !                   (1.0-VOFz(ic,jc,kc))*forcz(ic,jc,kc)*dens_ratio/(zlen)

      qcap(ic,jc,kc)=( -(h31+h32+h33) + fbz )  + &
      VOFz(ic,jc,kc)*forcz(ic,jc,kc)/zlen+   &
      (1.0-VOFz(ic,jc,kc))*forcz(ic,jc,kc)*dens_ratio/(zlen)

      enddo
      enddo
      enddo

      return
      end

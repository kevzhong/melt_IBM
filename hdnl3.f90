      subroutine hdnl3 
      use param
      use local_arrays, only: vy,vz,qcap,vx,temp,forcz
      use mpi_param, only: kstart,kend
      use phasefield
      !use mls_param,only: dens_ratio
      implicit none
      integer :: jc,kc
      integer :: km,kp,jmm,jpp,ic,im,ip
      real    :: h32,h33,h31
      real    :: udx1,udx2,udx3, fbz, T_interp
      real :: phi_interp

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
 



      ! phi_interp = 0.5 * ( phi(ic,jc,km) + phi(ic,jc,kc) )

      ! qcap(ic,jc,kc)=( -(h31+h32+h33) + fbz )  + &
      !             (1.0-phi_interp)*forcz(ic,jc,kc)/zlen 



      qcap(ic,jc,kc)=( -(h31+h32+h33) + fbz ) 


      enddo
      enddo
      enddo

      ! HIT forcing for stochastic scheme or ABC
      if (forcing .eq. 1) then
      if ( (which_hit .eq. 1) .or. (which_hit .eq. 2) ) then
      do kc=kstart,kend
      km=kc-1
      !kp=kc+1
      do jc=1,n2m
      !jm=jmv(jc)
      !jp=jpv(jc)
      do ic=1,n1m
      !im=imv(ic)
      !ip=ipv(ic)
            phi_interp = 0.5 * ( phi(ic,jc,km) + phi(ic,jc,kc) )
            qcap(ic,jc,kc) = qcap(ic,jc,kc) + (1.0 - phi_interp)*forcz(ic,jc,kc)/zlen
      enddo
      enddo
      enddo

      endif
      endif


      ! Phase-field volume penalty
      if (pfmode .eq. 1) then

      do kc=kstart,kend
      km=kc-1
      kp=kc+1
      do jc=1,n2m
      do ic=1,n1m
            phi_interp = 0.5 * ( phi(ic,jc,km) + phi(ic,jc,kc) )

            qcap(ic,jc,kc) = qcap(ic,jc,kc)  - phi_interp**2 * vz(ic,jc,kc) / (al * dt)

      enddo
      enddo
      enddo

      endif



      return
      end

      subroutine hdnl2 
      use param
      use phasefield
      use local_arrays, only: vx,vy,vz,dph,forcy
      use mpi_param, only: kstart,kend
      !use mls_param,only: dens_ratio
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip,km
      real    :: h22,h23,udx1,udx2,h21,udx3
      real :: phi_interp,fmean

      udx1=dx1*0.25
      udx2=dx2*0.25
      udx3=dx3*0.25
!
!     h term for the vy momentum equation at i+1/2,j,k+1/2
!
      do kc=kstart,kend
      km=kc-1
      kp=kc+1
      do jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)

!     vx vy term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   t      
!
      h21=( (vy(ip,jc,kc)+vy(ic,jc,kc)) &
           *(vx(ip,jc,kc)+vx(ip,jm,kc)) &
           -(vy(ic,jc,kc)+vy(im,jc,kc)) &
           *(vx(ic,jc,kc)+vx(ic,jm,kc)) &
          )*udx1
      
!     vy vy term
!
!
!                 d  q_r q_r 
!                ------------
!                 d   r      
!
      h22=( (vy(ic,jp,kc)+vy(ic,jc,kc)) &
           *(vy(ic,jp,kc)+vy(ic,jc,kc)) &
           -(vy(ic,jm,kc)+vy(ic,jc,kc)) &
           *(vy(ic,jm,kc)+vy(ic,jc,kc)) &
          )*udx2
!
!     vy vz term
!
!
!                 d  q_x q_r 
!                -----------
!                 d   x      
!
      h23=((vz(ic,jc,kp)+vz(ic,jm,kp))*(vy(ic,jc,kp)+vy(ic,jc,kc)) &
          -(vz(ic,jc,kc)+vz(ic,jm,kc))*(vy(ic,jc,kc)+vy(ic,jc,km)) &
          )*udx3


      ! phi_interp = 0.5 * ( phi(ic,jm,kc) + phi(ic,jc,kc) )
      ! dph(ic,jc,kc)=-(h21+h22+h23) + (1.0-phi_interp)*forcy(ic,jc,kc)/ylen 

      dph(ic,jc,kc)=-(h21+h22+h23)

      enddo
      enddo
      enddo

      ! HIT forcing for stochastic scheme or ABC
      if (forcing .eq. 1) then
      if ( (which_hit .eq. 1) .or. (which_hit .eq. 2) ) then
      do kc=kstart,kend
      !km=kc-1
      !kp=kc+1
      do jc=1,n2m
      jm=jmv(jc)
      !jp=jpv(jc)
      do ic=1,n1m
      !im=imv(ic)
      !ip=ipv(ic)
            phi_interp = 0.5 * ( phi(ic,jm,kc) + phi(ic,jc,kc) )
            dph(ic,jc,kc) = dph(ic,jc,kc) + (1.0 - phi_interp)*forcy(ic,jc,kc)/ylen
      enddo
      enddo
      enddo

      endif
      endif


      ! Phase-field volume penalty
      if ( (pfmode .eq. 1) .and. (ibtype .eq. 1) ) then

      !fmean = 0.0

      do kc=kstart,kend
      !km=kc-1
      !kp=kc+1
      do jc=1,n2m
      jm=jmv(jc)
      !jp=jpv(jc)
      do ic=1,n1m
      !im=imv(ic)
      !ip=ipv(ic)

            phi_interp = 0.5 * ( phi(ic,jm,kc) + phi(ic,jc,kc) )
            !dph(ic,jc,kc) = dph(ic,jc,kc) -  phi_interp**2 * vy(ic,jc,kc) / (al * dt)
            !dph(ic,jc,kc) = dph(ic,jc,kc) -  phi_interp * vy(ic,jc,kc) / ( al * dt)

            !fmean = fmean - phi_interp * vy(ic,jc,kc) 

            fmean = 4.0 * phi_interp * (1.0 - phi_interp)
            dph(ic,jc,kc) = dph(ic,jc,kc) -  fmean * vy(ic,jc,kc) / ( al * dt)


      enddo
      enddo
      enddo

      !call MpiAllSumRealScalar(fmean)
      !fmean = fmean / (al * dt * n1m * n2m * n3m)

      ! ! subtract mean forcing
      ! do kc=kstart,kend
      ! do jc=1,n2m
      ! do ic=1,n1m
      !       dph(ic,jc,kc) = dph(ic,jc,kc) - fmean
      ! enddo
      ! enddo
      ! enddo

      endif

      
      return
      end


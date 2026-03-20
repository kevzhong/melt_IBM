      subroutine invtr2 
      use param
      use local_arrays, only: vy,pr,rhs,dph,ru2,forcy
      use mpi_param, only: kstart,kend
      use phasefield
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real    :: udx2
      real    :: dcvy,dpx22
      real    :: d22vy,d33vy,d11vy
      real    :: alre,udx1q,udx2q,udx3q
      real :: phi_interp

      alre=al/ren
      udx2=dx2*al
      udx1q=dx1q
      udx2q=dx2q
      udx3q=dx3q

!
!  compute the rhs of the factored equation
!  everything at i,j+1/2,k+1/2
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

!   viscid terms
!
!   x- second derivative of vy
!
            d11vy=(vy(ip,jc,kc)-2.0*vy(ic,jc,kc)+vy(im,jc,kc))*udx1q

!   y- second derivative of vy

            d22vy=(vy(ic,jp,kc)-2.0*vy(ic,jc,kc)+vy(ic,jm,kc))*udx2q

!   z- second derivative of vy

            d33vy=(vy(ic,jc,kp)-2.0*vy(ic,jc,kc)+vy(ic,jc,km))*udx3q

            dcvy=d11vy+d22vy+d33vy
 
!
!   component of grad(pr) along 2 direction
!
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2

            rhs(ic,jc,kc)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc) &
                          +alre*dcvy-dpx22)*dt

            ! HIT forcing
            !rhs(ic,jc,kc) = rhs(ic,jc,kc) + forcy(ic,jc,kc) * al * dt

            ru2(ic,jc,kc)=dph(ic,jc,kc)
         enddo
       enddo
      enddo

      ! HIT forcing, linear forcing gets aldt coefficient
      if (forcing .eq. 1) then
      if (which_hit .eq. 3) then
        do kc=kstart,kend
        do jc=1,n2m
        !jm=jmv(jc)
        do ic=1,n1m
              !phi_interp = 0.5 * ( phi(ic,jm,kc) + phi(ic,jc,kc) )
              !rhs(ic,jc,kc) = rhs(ic,jc,kc) + (1.0 - phi_interp) * forcy(ic,jc,kc) * al * dt
              rhs(ic,jc,kc) = rhs(ic,jc,kc) + forcy(ic,jc,kc) * al * dt
        enddo
        enddo
        enddo
      endif
      endif



      call solxi(beta*al*dx1q)
      call solxj(beta*al*dx2q)
      call solxk(vy(1:n1,1:n2,kstart:kend),beta*al*dx3q)


      
      !vy(:,n2,:) = vy(:,1,:)
     
      return
      end

      subroutine invtr3 
      use param
      use local_arrays, only: vz,qcap, pr,ru3,rhs,forcz
      use mpi_param, only: kstart,kend
      use phasefield
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: udx3
      real    :: dvz2,dvz3,dcvz,dpx33,dvz1
      real    :: alre,udx1q,udx2q,udx3q
      real :: phi_interp

      alre=al/ren
      udx1q=dx1q
      udx2q=dx2q
      udx3q=dx3q
      udx3 = al*dx3
!
!  compute the rhs of the factored equation
!  everything at i+1/2,j+1/2,k
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

!   viscous terms
!
!   x- second derivatives of vz

            dvz1=(vz(im,jc,kc)-2.0*vz(ic,jc,kc)+vz(ip,jc,kc))*udx1q

!   y- second derivatives of vz

            dvz2=(vz(ic,jm,kc)-2.0*vz(ic,jc,kc)+vz(ic,jp,kc))*udx2q

!   z- second derivatives of vz

            dvz3=(vz(ic,jc,kp)-2.0*vz(ic,jc,kc)+vz(ic,jc,km))*udx3q
 
            dcvz=dvz2+dvz3+dvz1
!
!  component of grad(pr) along x3 direction
!
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3

            rhs(ic,jc,kc)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc) &
                          +alre*dcvz-dpx33)*dt 

            ! HIT forcing
            !rhs(ic,jc,kc) = rhs(ic,jc,kc) + forcz(ic,jc,kc) * al * dt

!  updating of the explicit terms

            ru3(ic,jc,kc)=qcap(ic,jc,kc)
         enddo
       enddo
      enddo

      ! HIT forcing, linear forcing gets aldt coefficient
      if (forcing .eq. 1) then
      if (which_hit .eq. 3) then
        do kc=kstart,kend
        !km=kc-1
        !kp=kc+1
        do jc=1,n2m
        !jm=jmv(jc)
        !jp=jpv(jc)
        do ic=1,n1m
        !im=imv(ic)
        !ip=ipv(ic)
          !phi_interp = 0.5 * ( phi(ic,jc,km) + phi(ic,jc,kc) )
          !rhs(ic,jc,kc) = rhs(ic,jc,kc) + (1.0 - phi_interp) * forcz(ic,jc,kc) * al * dt
          rhs(ic,jc,kc) = rhs(ic,jc,kc) + forcz(ic,jc,kc) * al * dt

        enddo
        enddo
        enddo
      endif
      endif

      call solxi(beta*al*dx1q)
      call solxj(beta*al*dx2q)
      call solxk(vz(1:n1,1:n2,kstart:kend),beta*al*dx3q)
 
      return
      end

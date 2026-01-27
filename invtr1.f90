      subroutine invtr1 
      use param
      use local_arrays, only: pr,rhs,ru1,vx,dq,forcx
      use mpi_param, only: kstart,kend
      use phasefield

      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real    :: udx1
      real    :: dcvx,dpx11
      real    :: d22vx,d33vx,d11vx
      real    :: alre,udx1q,udx2q,udx3q
      real    :: phi_interp

      alre=al/ren

      udx1=dx1*al
      udx1q=dx1q
      udx2q=dx2q
      udx3q=dx3q
!
!  compute the rhs of the factored equation
!  everything at i,j+1/2,k+1/2
!
!    points inside the flowfield
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

!
!    viscid terms
!
!   x- second derivative of vx
!
            d11vx=(vx(ip,jc,kc)-2.0*vx(ic,jc,kc)+vx(im,jc,kc))*udx1q

!   y- second derivative of vx

            d22vx=(vx(ic,jp,kc)-2.0*vx(ic,jc,kc)+vx(ic,jm,kc))*udx2q
 
!   z- second derivative of vx
 
            d33vx=(vx(ic,jc,kp)-2.0*vx(ic,jc,kc)+vx(ic,jc,km))*udx3q

            dcvx=d11vx+d22vx+d33vx
 
!
!   component of grad(pr) along 2 direction
!
            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1

            

            rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc) &
                          +alre*dcvx-dpx11)*dt

                

            ! ! Uniform pressure gradient for Poiseiulle flow!
            ! rhs(ic,jc,kc) = rhs(ic,jc,kc) + 1.0 * al * dt

            ! HIT forcing
            rhs(ic,jc,kc) = rhs(ic,jc,kc) + forcx(ic,jc,kc) * al * dt


            ru1(ic,jc,kc)=dq(ic,jc,kc)
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
        im=imv(ic)
        !ip=ipv(ic)
              phi_interp = 0.5 * ( phi(im,jc,kc) + phi(ic,jc,kc) )
              rhs(ic,jc,kc) = rhs(ic,jc,kc) + (1.0 - phi_interp) * forcx(ic,jc,kc) * al * dt
        enddo
        enddo
        enddo
      endif
      endif

      call solxi(beta*al*dx1q )
      call solxj(beta*al*dx2q )
      call solxk(vx(1:n1,1:n2,kstart:kend),beta*al*dx3q )

      ! call solxi_FSI(beta*al*dx1q, usolid_x(1:n2,1:n2,kstart:kend ) )
      ! call solxj_FSI(beta*al*dx2q, usolid_x(1:n2,1:n2,kstart:kend ) )
      ! call solxk_FSI(vx(1:n1,1:n2,kstart:kend),beta*al*dx3q, usolid_x(1:n2,1:n2,kstart:kend ) )

      !vx(n1,:,:) = vx(1,:,:)
     
      return
      end

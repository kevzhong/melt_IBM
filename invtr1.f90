      subroutine invtr1 
      use param
      use local_arrays, only: pr,rhs,ru1,vx,dq
      use mpi_param, only: kstart,kend

      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real    :: udx1
      real    :: dcvx,dpx11
      real    :: d22vx,d33vx,d11vx
      real    :: alre,udx1q,udx2q,udx3q

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

 
            if ( abs( usolid_x(ic,jc,kc) ) .gt. 0.0 ) then
              rhs(ic,jc,kc) = ( usolid_x(ic,jc,kc) - vx(ic,jc,kc) ) / (al*dt)
            else
            rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc) &
                          +alre*dcvx-dpx11)*dt
            endif

            ru1(ic,jc,kc)=dq(ic,jc,kc)
         enddo
       enddo

      enddo

      call solxi(beta*al*dx1q, usolid_x(1:n2,1:n2,kstart:kend ) )
      call solxj(beta*al*dx2q, usolid_x(1:n2,1:n2,kstart:kend ) )
      call solxk(vx(1:n1,1:n2,kstart:kend),beta*al*dx3q, usolid_x(1:n2,1:n2,kstart:kend ) )
      
      !vx(n1,:,:) = vx(1,:,:)
     
      return
      end

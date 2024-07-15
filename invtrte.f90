      subroutine invtrte 
      use param
      use local_arrays, only: temp, pr,rhs,rut,htemp
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: udx3
      real    :: dte2,dte3,dcte,dpx33,dte1
      real    :: alre,udx1q,udx2q,udx3q

      alre=al/ren
      udx1q=dx1q
      udx2q=dx2q
      udx3q=dx3q
      udx3 = al*dx3
      
      do kc=kstart,kend
        km=kc-1
        kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do ic=1,n1m
              im=imv(ic)
              ip=ipv(ic)

!   diffusive terms
!
!   x- second derivatives of temperature

            dte1=(temp(im,jc,kc)-2.0*temp(ic,jc,kc)+temp(ip,jc,kc))*udx1q

!   y- second derivatives of temperature

            dte2=(temp(ic,jm,kc)-2.0*temp(ic,jc,kc)+temp(ic,jp,kc))*udx2q

!   z- second derivatives of temperature

            dte3=(temp(ic,jc,kp)-2.0*temp(ic,jc,kc)+temp(ic,jc,km))*udx3q
 
            dcte=dte2+dte3+dte1
!
!  component of grad(pr) along x3 direction
!
            !dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3

            !rhs(ic,jc,kc)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc) &
            !              +alre*dcvz-dpx33)*dt 

            rhs(ic,jc,kc) =  dt * (ga*htemp(ic,jc,kc) + ro*rut(ic,jc,kc) ) + &
             (al/pec) * dcte * dt

!  updating of the explicit terms

            rut(ic,jc,kc)=htemp(ic,jc,kc)

         enddo
       enddo
      enddo

      call solxi(dt/pec*0.5d0*al*dx1q )
      call solxj(dt/pec*0.5d0*al*dx2q )
      call solxk(temp(1:n1,1:n2,kstart:kend),dt/pec*0.5d0*al*dx3q )

 
      return
      end

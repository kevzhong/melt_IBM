      subroutine invtrte_dirichlet 
      use param
      use local_arrays, only: temp, pr,rhs,rut,htemp
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: udx3
      real    :: dte2,dte3,dcte,dpx33,dte1
      real    :: alre,udx1q,udx2q,udx3q
      real :: Tconst



      Tconst = 1.0

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
            if (ic .eq. 1) then
              dte1=(temp(ip,jc,kc)-3.0*temp(ic,jc,kc)+ 2.0*Tconst)*udx1q
            elseif (ic .eq. n1m) then
              dte1=(2.0*Tconst-3.0*temp(ic,jc,kc)+ temp(im,jc,kc) )*udx1q
            else
              dte1=(temp(im,jc,kc)-2.0*temp(ic,jc,kc)+temp(ip,jc,kc))*udx1q
            endif

!   y- second derivatives of temperature
            if (jc .eq. 1) then
              dte2=(temp(ic,jp,kc)-3.0*temp(ic,jc,kc) + 2.0*Tconst )*udx2q
            elseif (jc .eq. n2m) then
              dte2=(2.0*Tconst - 3.0*temp(ic,jc,kc) + temp(ic,jm,kc) )*udx2q
            else
              dte2=(temp(ic,jm,kc)-2.0*temp(ic,jc,kc)+temp(ic,jp,kc))*udx2q
            endif

!   z- second derivatives of temperature
            if (kc .eq. 1) then
              dte3=(temp(ic,jc,kp)-3.0*temp(ic,jc,kc)+ 2.0*Tconst ) * udx3q
    
            elseif (kc .eq. n3m) then
              dte3=(2.0*Tconst-3.0*temp(ic,jc,kc)+temp(ic,jc,km))*udx3q
            else
              dte3=(temp(ic,jc,kp)-2.0*temp(ic,jc,kc)+temp(ic,jc,km))*udx3q
            endif

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

      call solxi_dirichlet(dt/pec*0.5d0*al*dx1q)
      call solxj_dirichlet(dt/pec*0.5d0*al*dx2q)
      call solxk_dirichlet(temp(1:n1,1:n2,kstart:kend),dt/pec*0.5d0*al*dx3q)
 
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcDissipation.F90                            !
!    CONTAINS: subroutine CalcDissipation                 !
!                                                         ! 
!    PURPOSE: Calculates the instantaneous kinetic        !
!     energy dissipation and writes it in dissip.out      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcDissipation
      use mpih
      use param
      use local_arrays,only: vx,vy,vz
      use stat_arrays
      use mls_param
      use mpi_param, only: kstart,kend

      implicit none
      integer :: ic,jc,kc
      integer :: im,ip,jm,jp,kp,km
      integer :: nfluid
      real :: eta,re_lam
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: diss_volAvg, nu, kmax_eta, Re_L,lambda_t
      real :: udx1,udx2,dissipte,udx3, kenerg, urms, L_int
      character(70) namfile
      
      diss_volAvg = 0.0d0
      kenerg=0.0d0

      nfluid = 0

      udx1=dx1
      udx2=dx2
      udx3=dx3

      call update_both_ghosts(n1,n2,vx,kstart,kend)
      call update_both_ghosts(n1,n2,vy,kstart,kend)
      call update_both_ghosts(n1,n2,vz,kstart,kend)

!================================================================
!      Dissipation rates
!================================================================
!
!

       do kc=kstart,kend
        kp=kc+1
        km=kc-1

       do jc=1,n2m
        jm=jmv(jc)
        jp=jpv(jc)

       do ic=1,n1m
       im= imv(ic)
       ip= ipv(ic)

!       if(imlsfor.eq.1) then
      !  if((VOFx(ip,jc,kc).eq.1).and.(VOFx(ic,jc,kc).eq.1).and. &
      !     (VOFx(ic,jp,kc).eq.1).and.(VOFx(ip,jp,kc).eq.1).and. &
      !     (VOFx(ic,jm,kc).eq.1).and.(VOFx(ip,jm,kc).eq.1).and. &
      !     (VOFx(ic,jc,kp).eq.1).and.(VOFx(ip,jc,kp).eq.1).and. &
      !     (VOFx(ic,jc,km).eq.1).and.(VOFx(ip,jc,km).eq.1).and. &
      !     (VOFy(ip,jc,kc).eq.1).and.(VOFy(ip,jp,kc).eq.1).and. &
      !     (VOFy(im,jc,kc).eq.1).and.(VOFy(im,jp,kc).eq.1).and. &
      !     (VOFy(ic,jp,kc).eq.1).and.(VOFy(ic,jc,kc).eq.1).and. &
      !     (VOFy(ic,jc,kp).eq.1).and.(VOFy(ic,jp,kp).eq.1).and. &
      !     (VOFy(ic,jc,km).eq.1).and.(VOFy(ic,jp,km).eq.1).and. &
      !     (VOFz(ip,jc,kc).eq.1).and.(VOFz(ip,jc,kp).eq.1).and. &
      !     (VOFz(im,jc,kc).eq.1).and.(VOFz(im,jc,kp).eq.1).and. &
      !     (VOFz(ic,jp,kc).eq.1).and.(VOFz(ic,jp,kp).eq.1).and. &
      !     (VOFz(ic,jm,kc).eq.1).and.(VOFz(ic,jm,kp).eq.1).and. &
      !     (VOFz(ic,jc,kp).eq.1).and.(VOFz(ic,jc,kc).eq.1))then

      if ( solid_mask(ic,jc,kc) .eqv. .false. ) then

      nfluid = nfluid + 1

       !du/dx
       h11=( vx(ip,jc,kc)-vx(ic,jc,kc) )*udx1

       !du/dy
       h12=( (vx(ic,jp,kc)+vx(ip,jp,kc))- &
             (vx(ic,jm,kc)+vx(ip,jm,kc)) )*0.25d0*udx2

      !du/dz
       h13=( (vx(ic,jc,kp)+vx(ip,jc,kp))- &
               (vx(ic,jc,km)+vx(ip,jc,km)) )*0.25d0*udx3

       !---

      !dv/dx
       h21=( (vy(ip,jc,kc)+vy(ip,jp,kc))- &
             (vy(im,jc,kc)+vy(im,jp,kc)) )*0.25d0*udx1

      !dv/dy
       h22=( vy(ic,jp,kc)-vy(ic,jc,kc) )*udx2

       !dv/dz
       h23=( (vy(ic,jc,kp)+vy(ic,jp,kp))- &
             (vy(ic,jc,km)+vy(ic,jp,km)) )*0.25d0*udx3

       !---

      !dw/dx
       h31=( (vz(ip,jc,kc)+vz(ip,jc,kp))- &
             (vz(im,jc,kc)+vz(im,jc,kp)) )*0.25d0*udx1

      !dw/dy
       h32=( (vz(ic,jp,kc)+vz(ic,jp,kp))- &
             (vz(ic,jm,kc)+vz(ic,jm,kp)) )*0.25d0*udx2

      !dw/dz
       h33=( vz(ic,jc,kp)-vz(ic,jc,kc) )*udx3

      !  kenerg = kenerg +  ((vx(ic,jc,kc)+vx(ip,jc,kc))*0.5)**2+ &
      !  &                   ((vy(ic,jc,kc)+vy(ic,jp,kc))*0.5)**2 + &
      !  &                   ((vz(ic,jc,kc)+vz(ic,jc,kp)*0.5))**2

       kenerg = kenerg + (vx(ic,jc,kc)**2+vx(ip,jc,kc)**2+ &
       &                   vy(ic,jc,kc)**2+vy(ic,jp,kc)**2 + &
       &                   vz(ic,jc,kc)**2+vz(ic,jc,kp)**2)*0.5d0


       !dissipte = 2.0*(h11**2+h22**2+h33**2)+ &
       !        (h21+h12)**2+ (h31+h13)**2+ (h32+h23)**2

      !sij sij = s11^2 + s22^2 + s33^2 + 2*(  s12^2 + s13^2 + s23^2  )
      dissipte = h11**2 + h22**2 + h33**2 + & ! s11^2 + s22^2 + s33^2            
                  0.5d0 * (h21+h12)**2 + & ! 2 * s12^2
                  0.5d0 * (h31+h13)**2 + & ! 2 * s13^2
                  0.5d0 * (h32+h23)**2     ! 2 * s23^2

      diss_volAvg = diss_volAvg+dissipte

      end if

       end do
       end do
       end do
           
      call MpiAllSumRealScalar(diss_volAvg)
      call MpiAllSumRealScalar(kenerg)

      ! kenerg =0.5d0* kenerg / (dble(n1m*n2m*n3m))
      ! nu=1.0d0/ren
      ! diss_volAvg = 2.0d0 * nu*diss_volAvg / (dble(n1m*n2m*n3m))


      ! Fluid-domain average
      kenerg =0.5d0* kenerg / (dble(nfluid))
      nu=1.0d0/ren
      diss_volAvg = 2.0d0 * nu*diss_volAvg / (dble(nfluid))


      urms = sqrt( 2.0 / 3.0 * kenerg )

      eta = (nu**3/diss_volAvg)**(0.25)
      kmax_eta = pi * dx1 * eta ! kmax_eta where kmax = pi /  dx is the Nyquist limit
      Re_L = kenerg**2 / (nu * diss_volAvg)

      lambda_t = sqrt( 15.0 * nu * urms**2 / diss_volAvg )
      re_lam = urms * lambda_t / nu

      L_int = sqrt(kenerg**3) / diss_volAvg ! Large-eddy length scale
      !re_lam = sqrt(20.0 / 3.0 * kenerg**2 / (diss_volAvg * nu) ) ! cf. eqn (6.59--6.65), Pope (2000)
      if(ismaster) then
            write(6,*) "Re_lambda", re_lam 
      namfile='stringdata/diss.txt'
      
      open(unit=92,file=namfile, Access='append', Status='unknown')
      !write(92,'(100E15.7)') time,dt,diss_volAvg,kmax_eta,eta,re_lam!,nu,ren
      write(92,'(100E15.7)') time, diss_volAvg, kenerg, urms, eta, lambda_t, Re_L, re_lam, L_int

      close(92)
      end if
   
      return
      end

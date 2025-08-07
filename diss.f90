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
      use local_arrays,only: vx,vy,vz,temp
      use local_aux,only: diss, tke,chi
      use stat_arrays
      use mls_param
      use mpi_param, only: kstart,kend

      implicit none
      integer :: ic,jc,kc
      integer :: im,ip,jm,jp,kp,km
      real :: sumvof
      real :: eta,re_lam
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: diss_volAvg, nu, kmax_eta, Re_L,lambda_t, chi_volAvg
      real :: udx1,udx2,dissipte,udx3, kenerg, urms, L_int
      character(70) namfile
      
      diss_volAvg = 0.0d0
      chi_volAvg = 0.0d0
      kenerg=0.0d0

      sumvof = 0.0

      udx1=dx1
      udx2=dx2
      udx3=dx3

      call update_both_ghosts(n1,n2,vx,kstart,kend)
      call update_both_ghosts(n1,n2,vy,kstart,kend)
      call update_both_ghosts(n1,n2,vz,kstart,kend)
      call update_both_ghosts(n1,n2,temp,kstart,kend)

      !-------------------- Re-tag cells --------------------------
      if ( .not. is_stationarySolid ) then
        call tagCells
      endif
 !-------------------- End re-tag cells --------------------------

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

       ! KZ: Fluid-domain average
      tke(ic,jc,kc) =  (vx(ic,jc,kc)**2+vx(ip,jc,kc)**2+ &
       &                   vy(ic,jc,kc)**2+vy(ic,jp,kc)**2 + &
       &                   vz(ic,jc,kc)**2+vz(ic,jc,kp)**2)*0.5d0

       tke(ic,jc,kc) = tke(ic,jc,kc) * 0.5


       !dissipte = 2.0*(h11**2+h22**2+h33**2)+ &
       !        (h21+h12)**2+ (h31+h13)**2+ (h32+h23)**2

      !sij sij = s11^2 + s22^2 + s33^2 + 2*(  s12^2 + s13^2 + s23^2  )
      ! dissipte = h11**2 + h22**2 + h33**2 + & ! s11^2 + s22^2 + s33^2            
      !             0.5d0 * (h21+h12)**2 + & ! 2 * s12^2
      !             0.5d0 * (h31+h13)**2 + & ! 2 * s13^2
      !             0.5d0 * (h32+h23)**2     ! 2 * s23^2

      diss(ic,jc,kc) = h11**2 + h22**2 + h33**2 + & ! s11^2 + s22^2 + s33^2            
            0.5d0 * (h21+h12)**2 + & ! 2 * s12^2
            0.5d0 * (h31+h13)**2 + & ! 2 * s13^2
            0.5d0 * (h32+h23)**2     ! 2 * s23^2

      diss(ic,jc,kc) = diss(ic,jc,kc) * 2.0 / ren

      !diss_volAvg = diss_volAvg+dissipte



      ! Scalar dissipation rate
      chi(ic,jc,kc) = (0.5 * dx1 * ( temp(ip,jc,kc) - temp(im,jc,kc) ) )**2 + & ! (dT / dx)^2
                      (0.5 * dx2 * ( temp(ic,jp,kc) - temp(ic,jm,kc) ) )**2 + & ! (dT / dy)^2
                      (0.5 * dx3 * ( temp(ic,jc,kp) - temp(ic,jc,km) ) )**2     ! (dT / dz)^2

      chi(ic,jc,kc) = chi(ic,jc,kc) / pec

      ! Only accumulate fluid-domain averages in fully fluid regions
      if (VOFp(ic,jc,kc) .eq. 1.0) then

            sumvof = sumvof + VOFp(ic,jc,kc)
            
            kenerg = kenerg + tke(ic,jc,kc) !* VOFp(ic,jc,kc)
            diss_volAvg = diss_volAvg + diss(ic,jc,kc) !* VOFp(ic,jc,kc)
            chi_volAvg = chi_volAvg + chi(ic,jc,kc) !* VOFp(ic,jc,kc)
      endif
            

       end do
       end do
       end do
           
      call MpiAllSumRealScalar(diss_volAvg)
      call MpiAllSumRealScalar(chi_volAvg)
      call MpiAllSumRealScalar(kenerg)
      call MpiAllSumRealScalar(sumVOF)

      ! kenerg =0.5d0* kenerg / (dble(n1m*n2m*n3m))
      ! nu=1.0d0/ren
      ! diss_volAvg = 2.0d0 * nu*diss_volAvg / (dble(n1m*n2m*n3m))


      ! Fluid-domain average
      !kenerg =0.5d0* kenerg / sumVOF
      kenerg = kenerg / sumVOF
      nu=1.0d0/ren
      !diss_volAvg = 2.0d0 * nu*diss_volAvg / sumVOF
      diss_volAvg = diss_volAvg / sumVOF
      chi_volAvg = chi_volAvg / sumVOF

      !write(*,*) "kenerg", kenerg
      !write(*,*) "diss_volAvg", diss_volAvg



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
      write(92,'(100E15.7)') time, diss_volAvg, kenerg, urms, eta, lambda_t, Re_L, re_lam, L_int,chi_volAvg

      close(92)
      end if
   
      return
      end

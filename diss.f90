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
      real :: eta,re_lam
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: nute, nu
      real :: udx1,udx2,dissipte,udx3, kenerg
      character(70) namfile
      
      nute = 0.0d0
      kenerg=0.0d0

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
!                                   1  |         | 2
!                   dissipation:  ---- | nabla  u|
!                                  Re  |         |
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
       if((VOFx(ip,jc,kc).eq.1).and.(VOFx(ic,jc,kc).eq.1).and. &
          (VOFx(ic,jp,kc).eq.1).and.(VOFx(ip,jp,kc).eq.1).and. &
          (VOFx(ic,jm,kc).eq.1).and.(VOFx(ip,jm,kc).eq.1).and. &
          (VOFx(ic,jc,kp).eq.1).and.(VOFx(ip,jc,kp).eq.1).and. &
          (VOFx(ic,jc,km).eq.1).and.(VOFx(ip,jc,km).eq.1).and. &
          (VOFy(ip,jc,kc).eq.1).and.(VOFy(ip,jp,kc).eq.1).and. &
          (VOFy(im,jc,kc).eq.1).and.(VOFy(im,jp,kc).eq.1).and. &
          (VOFy(ic,jp,kc).eq.1).and.(VOFy(ic,jc,kc).eq.1).and. &
          (VOFy(ic,jc,kp).eq.1).and.(VOFy(ic,jp,kp).eq.1).and. &
          (VOFy(ic,jc,km).eq.1).and.(VOFy(ic,jp,km).eq.1).and. &
          (VOFz(ip,jc,kc).eq.1).and.(VOFz(ip,jc,kp).eq.1).and. &
          (VOFz(im,jc,kc).eq.1).and.(VOFz(im,jc,kp).eq.1).and. &
          (VOFz(ic,jp,kc).eq.1).and.(VOFz(ic,jp,kp).eq.1).and. &
          (VOFz(ic,jm,kc).eq.1).and.(VOFz(ic,jm,kp).eq.1).and. &
          (VOFz(ic,jc,kp).eq.1).and.(VOFz(ic,jc,kc).eq.1))then

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

       kenerg = kenerg +  ((vx(ic,jc,kc)+vx(ip,jc,kc))*0.5)**2+ &
       &                   ((vy(ic,jc,kc)+vy(ic,jp,kc))*0.5)**2 + &
       &                   ((vz(ic,jc,kc)+vz(ic,jc,kp)*0.5))**2

      end if

      ! dissipte = 2.0*(h11**2+h22**2+h33**2)+ &
      !         (h21+h12)**2+ (h31+h13)**2+ (h32+h23)**2

      !sij sij = s11^2 + s22^2 + s33^2 + 2*(  s12^2 + s13^2 + s23^2  )
      dissipte = h11**2 + h22**2 + h33**2 + & ! s11^2 + s22^2 + s33^2            
                  0.5d0 * (h21+h12)**2 + & ! 2 * s12^2
                  0.5d0 * (h31+h13)**2 + & ! 2 * s13^2
                  0.5d0 * (h32+h23)**2     ! 2 * s23^2

       nute = nute+dissipte


       end do
       end do
       end do
     
    !  dx1=1/dx1
      
      call MpiAllSumRealScalar(nute)
      call MpiAllSumRealScalar(kenerg)
      kenerg =0.5d0* kenerg*(xlen*ylen*zlen)/(dble(n1m*n2m*n3m))

      nu=1.0d0/ren
      nute = 2.0d0 * nu*nute/((dx1*dx2*dx3)*(xlen*ylen*zlen))

      eta = (nu**3/nute)**(0.25)
      keta = float(n1m/2)*eta*2*pi
      !re_lam = sqrt(15.d0*(vxvyvz_rms_vol/sqrt(3.d0))**4/(nute*nu))
      re_lam = sqrt(20.0 / 3.0 * kenerg**2 / (nute * nu) ) ! cf. eqn (6.59--6.65), Pope (2000)
      if(ismaster) then
      namfile='flowmov/diss.txt'
      
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time,dt,nute,keta,eta,re_lam!,nu,ren
      close(92)
      end if
   
      return
      end

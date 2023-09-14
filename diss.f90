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
      use mpi_param

      implicit none
      integer :: ic,jc,kc
      integer :: im,ip,jm,jp,kp,km
      real    :: ti(2), tin(3)
      real    :: dmax,tpc,texit
      integer :: errorcode

      real :: eta,re_lam
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: nute, nu
      real :: udx1,udx2,dissipte,udx3
      character(70) namfile
      
      nute = 0.0d0

      udx1=dx1
      udx2=dx2
      udx3=dx3
      
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

       if(imlsfor.eq.1) then
       if((ax(ip,jc,kc).eq.1).and.(ax(ic,jc,kc).eq.1).and. &
          (ax(ic,jp,kc).eq.1).and.(ax(ip,jp,kc).eq.1).and. &
          (ax(ic,jm,kc).eq.1).and.(ax(ip,jm,kc).eq.1).and. &
          (ax(ic,jc,kp).eq.1).and.(ax(ip,jc,kp).eq.1).and. &
          (ax(ic,jc,km).eq.1).and.(ax(ip,jc,km).eq.1).and. &
          (ay(ip,jc,kc).eq.1).and.(ay(ip,jp,kc).eq.1).and. &
          (ay(im,jc,kc).eq.1).and.(ay(im,jp,kc).eq.1).and. &
          (ay(ic,jp,kc).eq.1).and.(ay(ic,jc,kc).eq.1).and. &
          (ay(ic,jc,kp).eq.1).and.(ay(ic,jp,kp).eq.1).and. &
          (ay(ic,jc,km).eq.1).and.(ay(ic,jp,km).eq.1).and. &
          (az(ip,jc,kc).eq.1).and.(az(ip,jc,kp).eq.1).and. &
          (az(im,jc,kc).eq.1).and.(az(im,jc,kp).eq.1).and. &
          (az(ic,jp,kc).eq.1).and.(az(ic,jp,kp).eq.1).and. &
          (az(ic,jm,kc).eq.1).and.(az(ic,jm,kp).eq.1).and. &
          (az(ic,jc,kp).eq.1).and.(az(ic,jc,kc).eq.1))then

       h11=( vx(ip,jc,kc)-vx(ic,jc,kc) )*udx1

       h12=( (vx(ic,jp,kc)+vx(ip,jp,kc))- &
             (vx(ic,jm,kc)+vx(ip,jm,kc)) )*0.25d0*udx2

       h13=( (vx(ic,jc,kp)+vx(ip,jc,kp))- &
               (vx(ic,jc,km)+vx(ip,jc,km)) )*0.25d0*udx3

       !---
       h21=( (vy(ip,jc,kc)+vy(ip,jp,kc))- &
             (vy(im,jc,kc)+vy(im,jp,kc)) )*0.25d0*udx1
       h22=( vy(ic,jp,kc)-vy(ic,jc,kc) )*udx2
       h23=( (vy(ic,jc,kp)+vy(ic,jp,kp))- &
             (vy(ic,jc,km)+vy(ic,jp,km)) )*0.25d0*udx3

       !---
       h31=( (vz(ip,jc,kc)+vz(ip,jc,kp))- &
             (vz(im,jc,kc)+vz(im,jc,kp)) )*0.25d0*udx1
       h32=( (vz(ic,jp,kc)+vz(ic,jp,kp))- &
             (vz(ic,jm,kc)+vz(ic,jm,kp)) )*0.25d0*udx2
       h33=( vz(ic,jc,kp)-vz(ic,jc,kc) )*udx3
       end if
       else
       h11=( vx(ip,jc,kc)-vx(ic,jc,kc) )*udx1

       h12=( (vx(ic,jp,kc)+vx(ip,jp,kc))- &
             (vx(ic,jm,kc)+vx(ip,jm,kc)) )*0.25d0*udx2

       h13=( (vx(ic,jc,kp)+vx(ip,jc,kp))- &
               (vx(ic,jc,km)+vx(ip,jc,km)) )*0.25d0*udx3

       !---
       h21=( (vy(ip,jc,kc)+vy(ip,jp,kc))- &
             (vy(im,jc,kc)+vy(im,jp,kc)) )*0.25d0*udx1
       h22=( vy(ic,jp,kc)-vy(ic,jc,kc) )*udx2
       h23=( (vy(ic,jc,kp)+vy(ic,jp,kp))- &
             (vy(ic,jc,km)+vy(ic,jp,km)) )*0.25d0*udx3

       !---
       h31=( (vz(ip,jc,kc)+vz(ip,jc,kp))- &
             (vz(im,jc,kc)+vz(im,jc,kp)) )*0.25d0*udx1
       h32=( (vz(ic,jp,kc)+vz(ic,jp,kp))- &
             (vz(ic,jm,kc)+vz(ic,jm,kp)) )*0.25d0*udx2
       h33=( vz(ic,jc,kp)-vz(ic,jc,kc) )*udx3
       end if

       dissipte = 2.0*(h11**2+h22**2+h33**2)+ &
               (h21+h12)**2+ (h31+h13)**2+ (h32+h23)**2


       nute = nute+dissipte
       end do
       end do
       end do


      call MpiAllSumRealScalar(nute)
      nu=1/ren
      nute = nute/ren/float(n1m*n2m*n3m)
      eta = (nu**3/nute)**(0.25)
      keta = float(n1m/2)*eta*2*pi
      re_lam = sqrt(15.d0*(vxvyvz_rms_vol/sqrt(3.d0))**4/(nute*nu))
      if(ismaster) then
      namfile='flowmov/diss.txt'
      
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time,dt,nute,keta,eta,re_lam!,nu,ren
      close(92)
      end if

      !if(nute.ge.9)then
      !errorcode = 1
      !call QuitRoutine(tin,.true.,errorcode)
      !end if 
   
      return
      end

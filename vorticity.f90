      subroutine vorticity
      use mpih
      use param
      use local_arrays, only: vx,vy,vz
      use local_aux, only: vorx, vory, vorz
      use mpi_param, only: kstart,kend
      use mls_param
      use stat_arrays
      implicit none
      real :: dvxx1,dvxx2,dvxx3
      real :: dvyx1,dvyx2,dvyx3
      real :: dvzx1,dvzx2,dvzx3
      integer :: kc,kp,jp,jm,jc,ic,im,ip,km

      call update_both_ghosts(n1,n2,vx,kstart,kend)
      call update_both_ghosts(n1,n2,vy,kstart,kend)
      call update_both_ghosts(n1,n2,vz,kstart,kend)

      call update_add_upper_ghost(vorx)
      call update_add_upper_ghost(vory)
      call update_add_upper_ghost(vorz)

      call update_add_lower_ghost(vorx)
      call update_add_lower_ghost(vory)
      call update_add_lower_ghost(vorz)

      do kc=kstart,kend
      kp=kc+1
      km=kc-1
        do jc=1,n2m
        jp=jpv(jc)
        jm=jmv(jc)
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)

          dvxx1=vx(ip,jc,kc)-vx(ic,jc,kc)

          dvxx2=( (vx(ic,jp,kc)+vx(ip,jp,kc))- &
                  (vx(ic,jm,kc)+vx(ip,jm,kc)) )*0.25d0

          dvxx3=( (vx(ic,jc,kp)+vx(ip,jc,kp))- &
                  (vx(ic,jc,km)+vx(ip,jc,km)) )*0.25d0

          !---
          dvyx1=( (vy(ip,jc,kc)+vy(ip,jp,kc))- &
                  (vy(im,jc,kc)+vy(im,jp,kc)) )*0.25d0

          dvyx2=vy(ic,jp,kc)-vy(ic,jc,kc)

          dvyx3=( (vy(ic,jc,kp)+vy(ic,jp,kp))- &
                  (vy(ic,jc,km)+vy(ic,jp,km)) )*0.25d0

          !---
          dvzx1=( (vz(ip,jc,kc)+vz(ip,jc,kp))- &
                  (vz(im,jc,kc)+vz(im,jc,kp)) )*0.25d0

          dvzx2=( (vz(ic,jp,kc)+vz(ic,jp,kp))- &
                  (vz(ic,jm,kc)+vz(ic,jm,kp)) )*0.25d0

          dvzx3=vz(ic,jc,kp)-vz(ic,jc,kc)


          vorx(ic,jc,kc) = dvzx2*dx2 - dvyx3*dx3
          vory(ic,jc,kc) = dvxx3*dx3 - dvzx1*dx1
          vorz(ic,jc,kc) = dvyx1*dx1 - dvxx2*dx2
          end do
        end do     
     end do
      call update_both_ghosts(n1,n2,vorx,kstart,kend)
      call update_both_ghosts(n1,n2,vory,kstart,kend)
      call update_both_ghosts(n1,n2,vorz,kstart,kend)
      return
      end subroutine vorticity

      subroutine circulation_y
      use param
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr
      use local_aux, only: vory

      IMPLICIT none

      integer ic,jc,kc
      integer ip,jp,kp
      integer inp,boxdim
      integer ii,kk
      real, allocatable, dimension(:,:) :: vor2
      real :: ycirc,udx1,udx3
      character(70) namfile
      allocate(vor2(n1m,n3m))

      udx1=1/dx1
      udx3=1/dx3
      if(Nparticle.eq.1)then
         do inp=1,Nparticle
            jc = pind1(2,inp)
         end do
      else
      jc = n2m/2
      end if
      jp=jc+1
      do kc=kstart,kend
         kp=kc+1

       do ic=1,n1m
        ip=ipv(ic)

        vor2(ic,kc) = (vory(ic,jc,kc)+vory(ic,jp,kc))*0.5

        end do
      end do

      ycirc=0.
      boxdim=n1m/2
      if(Nparticle.eq.1) then
         do inp=1,Nparticle
           do ii=1,boxdim+1
           ic=pind1(1,inp)-boxdim/2+(ii-1)
           kc=pind1(3,inp)-boxdim/2+(ii-1)

           ic=modulo(ic-1,n1m) + 1
           kc=modulo(kc-1,n1m) + 1
        
           ycirc=ycirc+vor2(ic,kc)
           end do
           end do
      end if
      call mpi_globalsum_double_var(ycirc)
      ycirc=((boxdim-1)**2*udx1*udx3*ycirc)/(boxdim**2)

 
      if(ismaster) then
      namfile='flowmov/circ_y.txt'
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time,ycirc
      close(92)
      end if
      return
      end subroutine circulation_y

      subroutine circulation_z
      use param
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr
      use local_aux, only: vorx,vory,vorz

      IMPLICIT none

      integer ic,jc,kc
      integer ip,jp,kp
      integer inp,boxdim
      integer ii,kk
      real, allocatable, dimension(:,:) :: vor3
      real :: zcirc,udx1,udx2
      character(70) namfile
      allocate(vor3(n1m,n2m))

      call update_both_ghosts(n1,n2,vorx,kstart,kend)
      call update_both_ghosts(n1,n2,vory,kstart,kend)
      call update_both_ghosts(n1,n2,vorz,kstart,kend)

      udx1=1/dx1
      udx2=1/dx2

      if(Nparticle.eq.1)then
         do inp=1,Nparticle
            kc = pind1(3,inp)
         end do
      else
      kc = n3m/2
      end if
      kp=kc+1

      if (kc.ge.kstart .and. kc.le.kend) then

      do jc=1,n2m
         jp=jc+1

       do ic=1,n1m
        ip=ipv(ic)

        vor3(ic,jc) = (vorz(ic,jc,kc)+vorz(ic,jc,kp))*0.5

        end do
      end do
     end if
      zcirc=0.
      boxdim=n1m/2
      if(Nparticle.eq.1)then
         do inp=1,Nparticle
           do ii=1,boxdim+1
           ic=pind1(1,inp)-boxdim/2+(ii-1)
           jc=pind1(2,inp)-boxdim/2+(ii-1)

           ic=modulo(ic-1,n1m) + 1
           kc=modulo(jc-1,n1m) + 1
        
           zcirc=zcirc+vor3(ic,jc)
           end do
           end do
      end if

      call mpi_globalsum_double_var(zcirc)
      zcirc=((boxdim-1)**2*udx1*udx2*zcirc)/(boxdim**2)

      if(ismaster) then
      namfile='flowmov/circ_z.txt'
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time,zcirc
      close(92)
      end if

      return
      end subroutine circulation_z


      subroutine calc_helicity
      use param
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr
      use local_aux, only: vorx,vory,vorz

      IMPLICIT none

      integer ic,jc,kc
      integer ip,jp,kp
      integer im,jm,km
      integer ii,kk
      real :: helicity
      real :: udx1,udx2,udx3
      character(70) namfile

      call update_both_ghosts(n1,n2,vorx,kstart,kend)
      call update_both_ghosts(n1,n2,vory,kstart,kend)
      call update_both_ghosts(n1,n2,vorz,kstart,kend)

      call update_both_ghosts(n1,n2,vx,kstart,kend)
      call update_both_ghosts(n1,n2,vy,kstart,kend)
      call update_both_ghosts(n1,n2,vz,kstart,kend)
 

      udx1=1/dx1
      udx2=1/dx2
      udx3=1/dx3

     helicity=0.0
 
     do kc=kstart,kend
        if(kc.ge.(n3m/4).and.kc.le.(3*n3m/4))then
        kp=kc+1
        km=kc-1
        do jc=n2m/4,3*n2m/4
           jm=jmv(jc)
           jp=jpv(jc)
           do ic=n1m/4,3*n1m/4
              im= imv(ic)
              ip= ipv(ic)

       if((imlsfor.eq.1).and.(ax(ip,jc,kc).eq.1).and.(ax(ic,jc,kc).eq.1).and. &
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


    
       helicity=helicity+(vx(ic,jc,kc)*vorx(ic,jc,kc)+ &
                          vx(ip,jc,kc)*vorx(ip,jc,kc)+ &
                          vy(ic,jc,kc)*vory(ic,jc,kc)+ &
                          vy(ic,jp,kc)*vory(ic,jp,kc)+ &
                          vz(ic,jc,kc)*vorz(ic,jc,kc)+ &
                          vz(ic,jc,kp)*vorz(ic,jc,kp))*0.5d0

      else


    
       helicity=helicity+(vx(ic,jc,kc)*vorx(ic,jc,kc)+ &
                          vx(ip,jc,kc)*vorx(ip,jc,kc)+ &
                          vy(ic,jc,kc)*vory(ic,jc,kc)+ &
                          vy(ic,jp,kc)*vory(ic,jp,kc)+ &
                          vz(ic,jc,kc)*vorz(ic,jc,kc)+ &
                          vz(ic,jc,kp)*vorz(ic,jc,kp))*0.5d0
 
      endif
          end do
        end do
      end if
      end do

      call mpi_globalsum_double_var(helicity)

      helicity=helicity/dble(n1m/2*n2m/2*n3m/2)

      if(ismaster) then
      namfile='flowmov/helicity.txt'
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time,helicity
      close(92)
      end if

      return
      end subroutine calc_helicity

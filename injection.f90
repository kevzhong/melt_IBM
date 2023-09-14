      subroutine CalcInjection
      use mpih
      use param
      use local_arrays, only: vx,vy,vz,forcx,forcy,forcz
      use mpi_param, only: kstart,kend
      use mls_param
      implicit none
      integer :: kc,jc,ic,ip,jp,kp
      real :: eps_in,kenerg
      character(70) namfile
      
      eps_in=0.0d0
      kenerg=0.0d0
      do kc=kstart,kend
         kp=kc+1
       do jc=1,n2m
          jp=jpv(jc)
        do ic=1,n1m
           ip=ipv(ic)
       
       if(imlsfor.eq.1)then
       if((ax(ic,jc,kc).eq.1).and.(ax(ip,jc,kc).eq.1).and. &
          (ay(ic,jc,kc).eq.1).and.(ay(ic,jp,kc).eq.1).and. &
          (az(ic,jc,kc).eq.1).and.(az(ic,jc,kp).eq.1))then

       eps_in = eps_in + ((forcx(ic,jc,kc)*vx(ic,jc,kc)+forcx(ip,jc,kc)*vx(ip,jc,kc))/xlen+ &
      &                    (forcy(ic,jc,kc)*vy(ic,jp,kc)+forcy(ic,jp,kc)*vy(ic,jp,kc))/ylen + &
      &                    (forcz(ic,jc,kc)*vz(ic,jc,kc)+forcz(ic,jc,kp)*vz(ic,jc,kp))/zlen)*0.5d0
 
       kenerg = kenerg + (vx(ic,jc,kc)**2+vx(ip,jc,kc)**2+ &
      &                   vy(ic,jp,kc)**2+vy(ic,jp,kc)**2 + &
      &                   vz(ic,jc,kc)**2+vz(ic,jc,kp)**2)*0.5d0

      end if
      else
       eps_in = eps_in + ((forcx(ic,jc,kc)*vx(ic,jc,kc)+forcx(ip,jc,kc)*vx(ip,jc,kc))/xlen+ &
      &                    (forcy(ic,jc,kc)*vy(ic,jp,kc)+forcy(ic,jp,kc)*vy(ic,jp,kc))/ylen + &
      &                    (forcz(ic,jc,kc)*vz(ic,jc,kc)+forcz(ic,jc,kp)*vz(ic,jc,kp))/zlen)*0.5d0
 
       kenerg = kenerg + (vx(ic,jc,kc)**2+vx(ip,jc,kc)**2+ &
      &                   vy(ic,jp,kc)**2+vy(ic,jp,kc)**2 + &
      &                   vz(ic,jc,kc)**2+vz(ic,jc,kp)**2)*0.5d0
      end if
        enddo
       enddo
      enddo
 
      call MpiAllSumRealScalar(eps_in)
      call MpiAllSumRealScalar(kenerg)

      eps_in = eps_in/(dble(n1m*n2m*n3m))
      kenerg =0.5d0* kenerg/(dble(n1m*n2m*n3m))
      if(ismaster) then
       namfile='flowmov/eps_in.txt'

       open(unit=92,file=namfile, Access='append', Status='unknown')
       write(92,'(100E15.7)') time,eps_in,kenerg
       close(92)
      end if
      end subroutine CalcInjection

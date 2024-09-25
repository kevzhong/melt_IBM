      subroutine CalcInjection
      use mpih
      use param
      use local_arrays
      use mpi_param
      use mls_param
      use mls_local
      implicit none
      integer :: kc,jc,ic,ip,jp,kp,inp
      real :: eps_in,kenerg,en_ibm
      real,dimension(3,2)     :: bbox_inds
      character(70) namfile

call update_both_ghosts(n1,n2,vx,kstart,kend)
call update_both_ghosts(n1,n2,vy,kstart,kend)
call update_both_ghosts(n1,n2,vz,kstart,kend)

 call update_add_upper_ghost(for_xc)
 call update_add_upper_ghost(for_yc)
 call update_add_upper_ghost(for_zc)

 call update_add_lower_ghost(for_xc)
 call update_add_lower_ghost(for_yc)
 call update_add_lower_ghost(for_zc)

 !-------------------- Re-tag cells --------------------------
 if (  (imlsfor .eq. 1 ) ) then
 do inp=1,Nparticle
   call get_bbox_inds(bbox_inds,inp)
   call tagCells(bbox_inds, inp)
 enddo
 endif
 !-------------------- End re-tag cells --------------------------

      eps_in=0.0d0
      kenerg=0.0d0
      en_ibm=0.0d0
      do kc=kstart,kend
         kp=kc+1
       do jc=1,n2m
          jp=jpv(jc)
        do ic=1,n1m
           ip=ipv(ic)
       
 !      if(imlsfor.eq.1)then
      !  if((VOFx(ic,jc,kc).eq.1).and.(VOFx(ip,jc,kc).eq.1).and. &
      !     (VOFy(ic,jc,kc).eq.1).and.(VOFy(ic,jp,kc).eq.1).and. &
      !     (VOFz(ic,jc,kc).eq.1).and.(VOFz(ic,jc,kp).eq.1))then
      if (solid_mask(ic,jc,kc) .eqv. .false.) then

       eps_in = eps_in + ((forcx(ic,jc,kc)*vx(ic,jc,kc)+forcx(ip,jc,kc)*vx(ip,jc,kc))/xlen+ &
      &                    (forcy(ic,jc,kc)*vy(ic,jp,kc)+forcy(ic,jp,kc)*vy(ic,jp,kc))/ylen + &
      &                    (forcz(ic,jc,kc)*vz(ic,jc,kc)+forcz(ic,jc,kp)*vz(ic,jc,kp))/zlen)*0.5d0
 
 !      kenerg = kenerg +  ((vx(ic,jc,kc)+vx(ip,jc,kc))*0.5)**2+ &
 !     &                   ((vy(ic,jc,kc)+vy(ic,jp,kc))*0.5)**2 + &
 !     &                   ((vz(ic,jc,kc)+vz(ic,jc,kp)*0.5))**2

 !     en_p=en_p + -1*(vz(ic,jc,kc)+vz(ic,jc,kp))*0.5/(dens_ratio-1)

!      endif
!     eps_in = eps_in +   dens_ratio* ((forcx(ic,jc,kc)*vx(ic,jc,kc)+forcx(ip,jc,kc)*vx(ip,jc,kc))/xlen+ &
!      &                    (forcy(ic,jc,kc)*vy(ic,jp,kc)+forcy(ic,jp,kc)*vy(ic,jp,kc))/ylen + &
!      &                    (forcz(ic,jc,kc)*vz(ic,jc,kc)+forcz(ic,jc,kp)*vz(ic,jc,kp))/zlen)*0.5d0
 
       kenerg = kenerg +  ((vx(ic,jc,kc)+vx(ip,jc,kc))*0.5)**2+ &
      &                   ((vy(ic,jc,kc)+vy(ic,jp,kc))*0.5)**2 + &
      &                   ((vz(ic,jc,kc)+vz(ic,jc,kp)*0.5))**2
      endif
       ! for_xc dimensionally is a velocity so has to be divided by dt and divided by dx1*dx2*dx3
       en_ibm = en_ibm + &
      &        (for_xc(ic,jc,kc)*vx(ic,jc,kc))/dt + &
      &        (for_yc(ic,jc,kc)*vy(ic,jc,kc))/dt + &
      &        (for_zc(ic,jc,kc)*vz(ic,jc,kc))/dt


  !    end if
        enddo
       enddo
      enddo
 
      call MpiAllSumRealScalar(eps_in)
      call MpiAllSumRealScalar(kenerg)
      call MpiAllSumRealScalar(en_ibm)

      eps_in = eps_in*(xlen*ylen*zlen)/(dble(n1m*n2m*n3m))
      kenerg =0.5d0* kenerg*(xlen*ylen*zlen)/(dble(n1m*n2m*n3m))
      !en_ibm=2.9103*en_ibm/dble(n1m*n2m*n3m)
      en_ibm = en_ibm/dble(n1m*n2m*n3m)

      if(ismaster) then
       namfile='stringdata/eps_in.txt'

       open(unit=92,file=namfile, Access='append', Status='unknown')
       write(92,'(100E15.7)') time,eps_in,kenerg,en_ibm
       close(92)
      end if
      end subroutine CalcInjection

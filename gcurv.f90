subroutine gcurv
use mpih
use mpi_param
use param
use local_arrays, only: vy,vz,pr,vx
use local_aux, only: vorx,vory,vorz
use mls_param
use coll_mod
implicit none
real    :: ti(2), tin(3)
real    :: dmax,tpc,texit
integer :: i,j,coll_num,p_old1,p_old2
real    :: p1(3),p2(3),dij(3),dista,dc,x1(3),x2(3)
real    :: AA1(3,3)
real    :: AA2(3,3)
integer :: errorcode
character(70) namfile
  call mpi_workdistribution
  call InitArrays

  tin(1) = MPI_WTIME()

  call HdfStart
  call InitStats
  call cordin 
  call phini
  call tri_geo

  time=0.d0
  vmax=0.0d0

      if(myid.eq.0) then
      write(6,754)n1m,n2m,n3m
  754 format(/,5x,'grid resolution: ',' n1m= ',i5,' n2m= ',i5, &
       ' n3m= ',i5/)                       
      write(6,755) 1.d0/dx1,1.d0/dx2,1.d0/dx3,dt,ntst                  
  755 format(/,2x,' dx1=',e10.3,' dx2=',e10.3,' dx3=',e10.3,' dt=' &
       ,e10.3,' ntst=',i7,/)
      endif

!===================================                                                                       
!   Read or create initial fields                                       
!===================================  


!     create the initial conditions

      if(nread.eq.0) then

       if(ismaster) write(6,*)' nread=0 ---> create initial conditions'

       ntime=0                                                         
       time=0.d0
       cflm=0.d0
         
       call inqpr 

      else

       if(ismaster) write(6,*)' nread=1 ---> Read initial conditions'
         call inirea
      endif                                                             

       call update_both_ghosts(n1,n2,pr,kstart,kend)
       call update_both_ghosts(n1,n2,vx,kstart,kend)
       call update_both_ghosts(n1,n2,vy,kstart,kend)
       call update_both_ghosts(n1,n2,vz,kstart,kend)

       call cfl 

!=================================
! Check velocity divergence of initial conditions
!=================================

      call divgck(dmax)

      if(ismaster) write(6,*)' initial divg dmax  ',dmax

!m================================
      if(ismaster) then
       write(6,711) tframe,ntst,tpin
711    format(3x,'check in cond : tframe =',f10.1,  &
             '  ntst=',i8,2x,'tpin =',f10.1//)
      endif
!m================================ 
!       If variable timestep

    
      cflm=cflm*dt

      call InitRandomForce

      tin(2) = MPI_WTIME()

      if(ismaster) then
       write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif

      allocate(ax(n1,n2,kstart-1:kend+1))
      allocate(ay(n1,n2,kstart-1:kend+1))
      allocate(az(n1,n2,kstart-1:kend+1))
!
!  ********* starts the time dependent calculation ***

      do ntime=1,ntst

       ti(1) = MPI_WTIME()

       call cfl 

       if(idtv.eq.1) then
         if(ntime.ne.1) then
            dt=cflmax/cflm
            if(dt.gt.dtmax) dt=dtmax
         endif
       else
         cflm=cflm*dt
       endif

       if (forcing.eq.1) then
       call CalcHITRandomForce
       endif

       call tschem

        if(mod(time,tpin).lt.dt) then
          if(ismaster) then
          write(6,*) "---------------"
          write(6,'(A,I10,A,E10.3)')"nt  ", ntime," time  ",time
          write(6,'(A,E10.3,A,E10.3)')"dt  ", dt,   " cfl   ",cflm*dt
          endif
        endif

                             

 !       call spec
          if(mod(time,tframe).lt.dt) then
           call findCMindices
!           call vorticity
!           call mkmov_hdf_ycut
!           call mkmov_hdf_xcut
!           call circulation_y
!           call circulation_z
 !             if (imlsfor.eq.1) then
!              call writePind
 !             call writePPpartpos
!              call writePPquat
 !             call write_tail_head
!              call write_shortdist
!              endif
          end if

       if((time.gt.tsta).and.(mod(ntime,5).eq.0)) then   !make it run every 4/5 time steps
         call vorticity
         call UpdateStats
         if(forcing.eq.1)then
           call CalcInjection
           call CalcDissipation
           call balance
         endif
        endif
              call write_shortdist
              call write_partvel
              call write_partrot
              call calc_helicity
!              call collision_stat
              call writePPpartpos
              call writePPquat
              call write_partacm
              call write_partalp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~collision_stat~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(ntime.eq.1)then
      coll_num=0
      else
      coll_num=coll_num
      p_old1=p_old1
      p_old2=p_old2
      endif
      !dc = 1./16.
      dc=6*xlen/n1m

      do i = 1,Nparticle
        do j = i+1,Nparticle

         call closest_distance(i,j,p1,p2,x1,x2,dista,AA1,AA2)
         
         if((dista.lt.dc).and.(p_old1 .ne. i).and.(p_old1 .ne. j))then
            coll_num=coll_num+1
            if(ismaster)then
            namfile='flowmov/coll_stat.txt'
            open(unit=92,file=namfile, Access='append', Status='unknown')
            write(92,'(1E15.7, 3I10.1, 200E15.7)') time, coll_num, i, j, x1(1), x1(2), x1(3), x2(1), x2(2), x2(3), quat(:,i), quat(:,j)
            close(92)
            endif
            p_old1=i
            p_old2=j
         endif
        end do
      end do
      coll_num=coll_num
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~collision_stat~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      time=time+dt

!     if((ntime.eq.1).or.(mod(time,tpin).lt.dt)) then 

!      call vmaxv !Make sure velocities are not exceeding maximum

     !  call cfl !Recalculate CFL 
     !  if(idtv.ne.1) cflm=cflm*dt 

     !  call divgck(dmax) !Make sure velocity is solenoidal

!     end if
!      tpc=0.001
!      if(mod(time,tpc).lt.dt) then
!      if(imlsfor.eq.1) then
!      call mpi_write_field
!      call mkmov_hdf_ycut
!      call writePPpartpos
!      call writePPquat
!      else
!      call mpi_write_field_noParts
!      end if
!      end if
       tpc=10
       if((ntime.eq.ntst).or.(mod(time,tpc).lt.dt)) then          !to perform when needed not only at the end
!      call WriteStats
      call mpi_write_field
      call write_tecplot_geom
      call mpi_write_continua
      call WriteRandForcCoef
      call continua_particle
      call continua_collision
      end if
!     enddo

!      texit=108
      if(mod(time,texit).lt.dt)then
!       if(mod(time,tframe).lt.dt)then
      call mpi_write_continua
      call mpi_write_field
      call WriteRandForcCoef
      call continua_particle
      call continua_collision

!      errorcode = 1
!      call QuitRoutine(tin,.true.,errorcode)

      end if
enddo
end

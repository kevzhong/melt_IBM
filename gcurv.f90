subroutine gcurv
use mpih
use mpi_param
use param
use local_arrays
use local_aux
use mls_param
use mls_local

!use coll_mod ! KZ: no collision module for now
implicit none
real    :: ti(2), tin(3)
real    :: dmax,tpc
real :: tstart, tend
integer :: inp
real,dimension(3,2)     :: bbox_inds
character(70) namfile
  call mpi_workdistribution
  call get_prow_pcol ! KZ: pencils for bounding box in ray-tagging
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
         
       !call inqpr_rotated
       call inqpr
       !call inqpr_taylorGreen

      else

       if(ismaster) write(6,*)' nread=1 ---> Read initial conditions'
         call inirea
      endif   
      
      if (temp_restart .eq. 1) then
        call restart_temperature
      endif

       call update_both_ghosts(n1,n2,pr,kstart,kend)
       call update_both_ghosts(n1,n2,vx,kstart,kend)
       call update_both_ghosts(n1,n2,vy,kstart,kend)
       call update_both_ghosts(n1,n2,vz,kstart,kend)
       call update_both_ghosts(n1,n2,temp,kstart,kend)

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
      if (nread .eq. 1) then
        call hdf_read_bhat
      endif

      tin(2) = MPI_WTIME()

      if(ismaster) then
       write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif
! 
!  ********* starts the time dependent calculation ***

      do ntime=1,ntst

       ti(1) = MPI_WTIME()

       call cfl 

       if(idtv.eq.1) then !CFL mode
         if(ntime.ne.1) then
            dt=dtmax
         endif
        call get_dt

       else ! Constant DT mode
         cflm=cflm*dt
       endif


       if (forcing.eq.1) then
        if (which_hit .eq. 1) then
          call CalcHITRandomForce
        elseif (which_hit .eq. 2) then
          call CalcABC_HITForce
        endif
       endif

       call tschem

        !if(mod(time,tpin).lt.dt) then !KZ: commented to dump at every timestep
          if(ismaster) then
          write(6,*) "---------------"
          write(6,'(A,I10,A,E10.3)')"nt  ", ntime," time  ",time
          write(6,'(A,E10.3,A,E10.3)')"dt  ", dt,   " cfl   ",cflm*dt
          write(6,*) "Ntri", count(isGhostFace(:,1) .eqv. .false.)
          write(6,*) "V(t)/VE", Volume(1) / celvol
          write(6,'(A,F10.6)') "min elength/dx:", minval( pack(eLengths(:,:) , .not. isGhostEdge(:,:)  ) )*dx1 
          write(6,'(A,F10.6,F10.6,F10.6)') "vel_CM:", vel_CM(:,1)
          write(6,'(A,F10.6,F10.6,F10.6)') "omega_CM:", omega_c(:,1)
          endif
        !endif


 
          if(mod(time,tframe).lt.dt) then !KZ: comment to dump cuts at every timestep
           !call mkmov_hdf_xcut
           call mkmov_hdf_ycut
           call mkmov_hdf_zcut
           call write_tecplot_geom
           !call mpi_write_tempField
           !call mpi_write_vel
           !call mpi_write_field
         endif
         
        ! ASCII write
          call writePartVol
          call writeInertTens
          call write_partrot
          call write_partpos
          call write_partvel
          call writeStructLoads

          call writeTriMeshStats
          call writeClock

          call CalcInjection
          call CalcDissipation

      time=time+dt
      if((ntime.eq.ntst).or.(mod(ntime,1000).eq.0)) then          !to perform when needed not only at the end
      call mpi_write_continua
      call mpi_write_field
      call WriteRandForcCoef
      call continua_particle
     end if

     ti(2) = MPI_WTIME()
     wtime_total = ti(2) - ti(1)

      enddo
end

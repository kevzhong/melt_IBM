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
integer :: icut,jcut,kcut
real,dimension(3,2)     :: bbox_inds
character(70) namfile
  call mpi_workdistribution
  call get_prow_pcol ! KZ: pencils for bounding box in ray-tagging
  call InitArrays

  if (specflag) call initSpectra

  tin(1) = MPI_WTIME()

  call HdfStart
  call InitStats
  call cordin 
  call phini
  call tri_geo


  ! Initial cell-tagging operation at start of runtime if IBM is active
  ! the most expensive full tagging of ALL cells
  initial_tag = .true.
  if (imlsfor .eq. 1) call tagCells
  initial_tag = .false. ! Subsequent time-steps: only tag along a narrow band

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
       !call inqpr
       call inqpr_taylorGreen

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
          !write(6,'(A,F10.6,F10.6,F10.6)') "vel_CM:", vel_CM(:,1)
          write(6,'(A,F10.6,F10.6,F10.6)') "pos_CM:", pos_CM(:,1)
          write(6,'(A,F10.6,F10.6,F10.6)') "omega_CM:", omega_c(:,1)
          endif
        !endif


          !------ ASCII write -----------------
          call writePartVol
          call writeInertTens
          call write_partrot
          call write_partpos
          call write_partvel
          call writeStructLoads

          call writeTriMeshStats
          call writeClock

          !call CalcInjection
          call CalcDissipation

          ! KZ: relative Lagrangian motion tracking
          call calcFluidVelAvgs
          !call calcRelShellVel
          call calcLocalShellFlow
          call vorticity
          !------ END ASCII -----------------


          if(mod(time,tframe).lt.dt) then !KZ: comment to dump cuts at every timestep
           
            if (imlsstr .eq. 1) then
              ! For Lagrangian: cuts follow centroid of object
              icut = floor( pos_CM(1,1) * dx1 ) + 1
              jcut = floor( pos_CM(2,1) * dx2 ) + 1
              kcut = floor( pos_CM(3,1) * dx3 ) + 1

              icut = modulo(icut-1,n1m)  + 1
              jcut = modulo(jcut-1,n2m)  + 1
              kcut = modulo(kcut-1,n3m)  + 1
            else
              icut = n1m / 2
              jcut = n2m / 2
              kcut = n3m / 2
            endif

           call mkmov_hdf_xcut(icut)
           call mkmov_hdf_ycut(jcut)
           call mkmov_hdf_zcut(kcut)
           call write_tecplot_geom
           !call mpi_write_tempField
           !call mpi_write_vel
           !call mpi_write_field

           if (specflag) call compute_1d_spectra
         endif
         


      time=time+dt

      ! KZ: These are executed at the end by QuitRoutine() instead
      !-------------------------------------------------------------------------------------------------
      !if((ntime.eq.ntst).or.(mod(ntime,1000).eq.0)) then !to perform when needed not only at the end
    !   if( (ntime.eq.ntst) ) then !to perform when needed not only at the end
    !   call mpi_write_continua
    !   !call mpi_write_field
    !   call WriteRandForcCoef
    !   call continua_particle
    !  end if
    !-------------------------------------------------------------------------------------------------

     ti(2) = MPI_WTIME()
     wtime_total = ti(2) - ti(1)

      enddo
end

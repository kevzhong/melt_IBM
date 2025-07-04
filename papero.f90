      program papero
      use mpih
      use mpi_param
      use param
      use local_arrays, only: vy,vz,pr,vx
      use mls_param
      implicit none
      character(len=4) :: dummy
      integer :: n,ns,nt,errorcode

      integer :: ntstf
      real    :: dmax,minwtdt
      real    :: ti(2), tin(3)

      call InitializeMPI

      tin(1) = mpi_wtime()
      tin(2) = mpi_wtime()

!RO   Calculate number of OpenMP Threads

      nt = 0
      nt = nt + 1

      if (ismaster) then 
        write(6,*) 'MPI tasks=', numtasks
        write(6,*) 'OMP threads per task=', nt
        write(6,*) 'No. of processors=', nt*numtasks
      end if

      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1m,n2m,n3m
        read(15,301) dummy
        read(15,*) xlen,ylen,zlen
        read(15,301) dummy
        read(15,*) nread, pread 
        read(15,301) dummy
        read(15,*) ntst,nsst,tframe,tpin,ireset
        read(15,301) dummy
        read(15,*) idtv, dt, dtmax, cflfix
        read(15,301) dummy
        read(15,*) ren,prandtl,betagz
        read(15,301) dummy
        read(15,*) Tmelt , Tliq, Tsol, latHeat, cpliquid, temp_restart
        read(15,301) dummy
        read(15,*) forcing 
      close(15)

      open(unit=15,file='HITForcing.in',status='old')
      read(15,301) dummy
      read(15,*) which_hit
      read(15,301) dummy
      read(15,*) tl, epsstar, kf_on_kmin
      read(15,301) dummy
      read(15,*) C_HIT
    close(15)

      open(unit=15,file='part.in',status='old')
        read(15,301) dummy       
        read(15,*) imlsfor, imlsstr, imelt
        read(15,301) dummy
        read(15,*) dens_ratio
        read(15,301) dummy       
        read(15,*) gtsfx, rad_p, tagType
        read(15,301) dummy       
        read(15,*) iremesh, PERC_Ethresh, V_ON_VE_PERC
        !read(15,*) iremesh, PERC_Athresh, skew_thresh, V_ON_VE_PERC


      close(15)

301     format(a4)          

      ! KZ Verify correctness of bou.in,  part.in

      ! Can only melt if MLS forcing is enabled
      if (imelt.eq.1) then
            if(imlsfor.ne.1) then
            write(*,*) "Rank", myid, "Melting enabled but MLS forcing disabled, exiting"
            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
            endif
      endif

      ! Can only do FSI if MLS forcing is enabled
      if (imlsstr.eq.1) then
            if(imlsfor.ne.1) then
            write(*,*) "Rank", myid, "FSI enabled but MLS forcing disabled, exiting"
            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
            endif
      endif

      ! Logical flag for solid
      ! If object is stationary (true by default on on initialization), only need to tag cells once
      is_stationarySolid = .true.
      if (imlsfor .eq. 1) then
            if ( (imlsstr .eq. 1) .or. (imelt .eq. 1) ) then ! non-stationary if melting or FSI enabled
                  is_stationarySolid = .false.
            endif
      endif

      gtsfx = "gts/" // trim(gtsfx)

      starea = 0 ! KZ: fix to zero, not using stat.f90


      pec = ren * prandtl


      n1=n1m+1                                                          
      n2=n2m+1  
      n3=n3m+1
      n1mh = n1m/2 + 1
      n2mh = n2m/2 + 1

      nsst=3
      if(nsst.eq.3)then
      gam(1)=8.d0/15.d0
      gam(2)=5.d0/12.d0
      gam(3)=3.d0/4.d0
      rom(1)=0.d0
      rom(2)=-17.d0/60.d0
      rom(3)=-5.d0/12.d0
      endif
      if(nsst.eq.1)then
      gam(1)=3./2.
      gam(2)=0.d0
      gam(3)=0.d0
      rom(1)=-1./2.
      rom(2)=0.d0
      rom(3)=0.d0
      endif


      do ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
      end do

  call MpiBarrier


      pi=2.d0*dasin(1.d0)                          

      if(ismaster) then
!m====================================================                                                                             
      write(6,112)xlen/zlen,ylen/zlen
  112 format(//,20x,'H I T',//,10x, &
       '3D Cube with aspect-ratio:  L1/L3 = ',f5.2,'   L2/L3 = ',f5.2)
      write(6,120)nsst
  120 format(/,5x, &
       'nsst', i7,/)    
      write(6,202) ren, prandtl
  202 format(/,5x,'Parameters: ',' Re=',e10.3,'       Prandtl = ',f5.2)
      if(idtv.eq.1) then
         write(6,204) cflfix
  204 format(/,5x,'Variable dt and fixed cfl= ', &
       e11.4,/ )            
      else 
         write(6,205) dtmax,cfllim
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=', &
        e11.4,/ )            
        write(*,*) "IMLSFOR = ", imlsfor, "IMLSSTR = ", imlsstr, "IMELT = ",imelt
      endif
!m====================================================    
      endif
      
!m======================================================
      
#ifdef TIMED
      timeflag = .true.
#endif

#ifdef SPEC
      specflag = .true.
#endif

      call  gcurv

      errorcode = 1
      call QuitRoutine(tin,.true.,errorcode)
      end

subroutine CalcHITRandomForce
      use param
      use mpih
      use mpi_param, only: kstart, kend
      use local_arrays, only: forcx, forcy, forcz

      implicit none
      integer :: l, m, nn,ndir
      integer :: i, j, k
      real :: kl, km, kn, kmag, kf
      real :: u1, u2
      real :: kvec(3)
      complex :: rand_complex
      complex :: Fhat(3,nmodes,nmodes,nmodes)
      complex :: buff1
  
      real :: tstart, tend
          ! il, jm, knn
  
      kf = kf_on_kmin * 2.0 * pi / xlen

      Fhat = 0.0
  
      if (ismaster) then ! Only for 1 process to sync random numbers
          do l = 1,nmodes ! x wavenumber
              kl = float(waveN(l))*2.0*pi / xlen
              do m = 1,nmodes ! y wavenumber
                  km = float(waveN(m))*2.0*pi / ylen
                  do nn = 1,nmodes ! z wavenumber
                      kn = float(waveN(nn))*2.0*pi / zlen
                      kmag = sqrt(kl**2 + km**2 + kn**2)
  
                      if  ( (kmag .gt. 0.0 ) .and. (kmag .lt. kf) ) then
                          do ndir = 1,3
                              ! Box-Muller transform
                              call random_number(u1)
                              call random_number(u2)
                              rand_complex = cmplx( sqrt(-2.0*log(u1))*cos(2.0*pi*u2) , &
                                                     sqrt(-2.0*log(u1))*sin(2.0*pi*u2) )
  
                              !write(*,*) "factor:", 2.0*epsstar*dt/TL**2  
  
  
                              bhat(ndir,l,m,nn) = bhat(ndir,l,m,nn) * (1.0 - dt / TL) + &
                                                 rand_complex*sqrt(2.0*epsstar*dt/(TL**2))  
                                              
  
                          enddo ! end 1:3
  
                                            ! Subtract projection along wavenumber: enforce solenoidality
                      kvec = [kl,km,kn]
                      Fhat(1:3,l,m,nn) = bhat(1:3,l,m,nn) - kvec * dot_product(kvec,  bhat(1:3,l,m,nn) ) / kmag**2
                      !write(*,*) "Fhat(1:3,l,m,nn) is :", Fhat(1:3,l,m,nn)
  
                      endif ! end kmag
    
                  enddo ! end nn
              enddo !end m
          enddo !end l
      endif
  
      ! Hard-set the zero mode to have no contribution
      !Fhat(1:3,(nmodes-1)/2+1,(nmodes-1)/2+1,(nmodes-1)/2+1 ) = 0.0 !/
      call MPI_BCAST(Fhat, 3*nmodes**3, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  
  
      ! Begin accmulation of terms
  
      forcx = 0.0
      forcy = 0.0
      forcz = 0.0
  
  
  
      ! x forcing
      do nn = 1,nmodes
          do m = 1,nmodes
              do l = 1,nmodes
                  do k = kstart,kend
                      do j = 1,n2m
                          buff1 = exp_I_km_yj(j,m) * exp_I_kn_zk(k,nn)
                          do i = 1,n1m
                              forcx(i,j,k) = forcx(i,j,k) + &
                              real( Fhat(1,l,m,nn) * exp_I_kl_xsi(i,l) * buff1 )
                          enddo
                      enddo
                  enddo
              enddo
          enddo
      enddo
  
  
      ! y forcing
      do nn = 1,nmodes
          do m = 1,nmodes
              do l = 1,nmodes
                  do k = kstart,kend
                      do j = 1,n2m
                          buff1 = exp_I_km_ysj(j,m) * exp_I_kn_zk(k,nn)
                          do i = 1,n1m
                              forcy(i,j,k) = forcy(i,j,k) + &
                              real( Fhat(2,l,m,nn) * exp_I_kl_xi(i,l)  * buff1 )
                          enddo
                      enddo
                  enddo
              enddo
          enddo
      enddo
  
      ! z forcing
      do nn = 1,nmodes
          do m = 1,nmodes
              do l = 1,nmodes
                  do k = kstart,kend
                      do j = 1,n2m
                          buff1 = exp_I_km_yj(j,m) * exp_I_kn_zsk(k,nn)
                          do i = 1,n1m
                              forcz(i,j,k) = forcz(i,j,k) + &
                              real( Fhat(3,l,m,nn) * exp_I_kl_xi(i,l)  * buff1 )
                          enddo
                      enddo
                  enddo
              enddo
          enddo
      enddo
  
  end subroutine CalcHITRandomForce

      subroutine InitRandomForce
            use param
            use mpih
            use mpi_param, only: kstart, kend

            implicit none
            integer :: Nbuffer, i,j,k,nn, count
            real :: kx, ky, kz
            integer :: l,m
            complex :: im=(0.,1.)
            real :: kappa,kf, magK
        
        
            ! First count the number of modes 
        
            ! {0, k1, k2, ........,  k_nmodes}    where k_nmodes < kf
        
            !Nbuffer counts the no. of positive wavenumbers + the zero mode
            Nbuffer = floor(kf_on_kmin) + 1
            
            !nmodes counts all -ve and +ve wavenumbers in 1 direciton, including the zero mode
            nmodes = (Nbuffer-1)*2 + 1
        
            allocate( waveN(nmodes) )
        
            ! Coefficient in forcing for later
            allocate( bhat(3,nmodes,nmodes,nmodes) )
            bhat = 0.0
        
        
            ! Allocate zero wavneumbers
            do i = 1,Nbuffer-1
                waven(i) = -( Nbuffer - i)
            enddo
        
            ! Allocate zero and positive wavenumbers
            do i = 0,Nbuffer-1
                waveN(Nbuffer+i) = i
            enddo
        
            ! Exponential pre-factor terms in Fourier expansion
            allocate( exp_I_kl_xi(n1m, nmodes) ) ! Nx x k_l
            allocate( exp_I_km_yj(n2m, nmodes) ) ! Ny x k_m
            allocate( exp_I_kn_zk(kstart:kend, nmodes) ) ! Nz x k_n
        
            ! Staggered array copies
            allocate( exp_I_kl_xsi(n1m, nmodes) ) ! Nx x k_l
            allocate( exp_I_km_ysj(n2m, nmodes) ) ! Ny x k_m
            allocate( exp_I_kn_zsk(kstart:kend, nmodes) ) ! Nz x k_n
        
        
            ! Duplicates needed for staggered grid
        
        
            ! x-modes
            do i = 1,n1m
                do nn=1,nmodes
                    kappa = float(waveN(nn))*2.0*pi / xlen
                    exp_I_kl_xi(i,nn) = exp( im * kappa * xm(i) )
                    exp_I_kl_xsi(i,nn) = exp( im * kappa * xc(i) )
        
                enddo
            enddo
        
            ! y-modes
            do j = 1,n2m
                do nn=1,nmodes
                    kappa = float(waveN(nn))*2.0*pi / ylen
                    exp_I_km_yj(j,nn) = exp( im * kappa * ym(j) )
                    exp_I_km_ysj(j,nn) = exp( im * kappa * yc(j) )
        
                enddo
            enddo
        
        
            ! z-modes
            do k = kstart,kend
                do nn=1,nmodes
                    kappa = float(waveN(nn))*2.0*pi / zlen
                    exp_I_kn_zk(k,nn) = exp( im * kappa * zm(k) )
                    exp_I_kn_zsk(k,nn) = exp( im * kappa * zc(k) )
                enddo
            enddo
        
            !write(*,*) "waveN is: ", exp_I_kn_zsk(kstart,nmodes)
        
            ! Count the number of modes
            count = 0
        
            kf = kf_on_kmin * 2.0*pi/xlen
        
            do l = 1,nmodes
                kx = float( waveN(l) )*2.0*pi/xlen
                do m=1,nmodes
                    ky = float( waveN(m) )*2.0*pi/ylen
                    do nn = 1,nmodes
                        kz = float( waveN(nn) )*2.0*pi/zlen
                        magK = sqrt( kx**2 + ky**2 + kz**2)
                        if ( (magK .le. kf) .and. (magK .ne. 0)) then
                            count = count + 1
                        endif
                enddo
            enddo
        enddo
        
        if (myid .eq. 0) write(*,*) count," wavenumber modes will be forced!"
        
end subroutine InitRandomForce

subroutine CalcABC_HITForce
      use param
      use local_arrays, only: forcx,forcy,forcz
      use mpih
      use mpi_param
      use stat_arrays
      implicit none
      integer :: j,k,i
      real :: C

      C = C_HIT * (2.0*pi / xlen)**2

      do k=kstart,kend
            do j=1,n2m
                  do i=1,n1m
                        forcx(i,j,k) = C * sin(2.0 * pi * zm(k) / zlen ) + C * cos(2.0 * pi * ym(j) / ylen )
                        forcy(i,j,k) = C * sin(2.0 * pi * xm(i) / xlen ) + C * cos(2.0 * pi * zm(k) / zlen )
                        forcz(i,j,k) = C * sin(2.0 * pi * ym(j) / ylen ) + C * cos(2.0 * pi * xm(i) / xlen )
                  enddo
            enddo
      enddo


      return
end


  subroutine WriteRandForcCoef
    use mpih
    use param
    use hdf5
    IMPLICIT NONE
    integer hdf_error, ndims
    character(40) filnambc
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_coef, dset_coef2
    integer(HID_T) :: dspace_coef
    integer(HSIZE_T) :: dims_coef(4)
    
    
    
    ! Size of bhat is
    ! bhat(3,nmodes,nmodes,nmodes) ) as a complex double
    
    if(ismaster) then 
    
    ndims=4
    dims_coef(1)=3
    dims_coef(2)=nmodes
    dims_coef(3)=nmodes
    dims_coef(4)=nmodes
    
    filnambc = 'continuation/continua_bhat.h5'
    call h5fcreate_f(filnambc,H5F_ACC_TRUNC_F, file_id, hdf_error)

    call h5screate_simple_f(ndims, dims_coef, dspace_coef, hdf_error)

    call h5dcreate_f(file_id, 'bhat_re', H5T_NATIVE_DOUBLE, &
                    dspace_coef, dset_coef, hdf_error)

    call h5dcreate_f(file_id, 'bhat_im', H5T_NATIVE_DOUBLE, &
                    dspace_coef, dset_coef2, hdf_error)

    call h5dwrite_f(dset_coef, H5T_NATIVE_DOUBLE, real(bhat), &
           dims_coef,hdf_error)

    call h5dwrite_f(dset_coef2, H5T_NATIVE_DOUBLE, aimag(bhat), &
           dims_coef,hdf_error)

    call h5dclose_f(dset_coef, hdf_error)
    call h5sclose_f(dspace_coef, hdf_error)


    call h5fclose_f(file_id, hdf_error)

    !filnambc = 'continuation/continua_bhat_im.h5'
    !call h5fcreate_f(filnambc,H5F_ACC_TRUNC_F, file_id, hdf_error)
    !call h5screate_simple_f(ndims, dims_coef, dspace_coef, hdf_error)
    !call h5dcreate_f(file_id, 'bhat_im', H5T_NATIVE_DOUBLE, &
    !                dspace_coef, dset_coef, hdf_error)
    !call h5dwrite_f(dset_coef, H5T_NATIVE_DOUBLE, aimag(bhat), &
    !       dims_coef,hdf_error)
    !call h5dclose_f(dset_coef, hdf_error)
    !call h5sclose_f(dspace_coef, hdf_error)
    !call h5fclose_f(file_id, hdf_error)
    
    endif
    
    return                                    
end      

subroutine hdf_read_bhat
    use hdf5
    use mpih
    use param
    implicit none
    character(200) :: dsetname,filename
    integer(HID_T) :: file_id
    integer(HID_T) :: dset
    integer :: hdf_error
    integer(HSIZE_T) :: dims(4)
    real :: re_buffer(3,nmodes,nmodes,nmodes)
    real :: im_buffer(3,nmodes,nmodes,nmodes)

    
    
    if (myid.eq.0) then
        dims(1) = 3
        dims(2) = nmodes
        dims(3) = nmodes
        dims(4) = nmodes

    !write(*,*) "Attempting to read bhat_re!"

    filename = 'continuation/continua_bhat.h5'
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf_error)
    !write(*,*) "Finishing open_f!"

    call h5dopen_f(file_id, 'bhat_re', dset, hdf_error)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
    re_buffer, dims, hdf_error)
    call h5dclose_f(dset, hdf_error)


    call h5dopen_f(file_id, 'bhat_im', dset, hdf_error)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
    im_buffer, dims, hdf_error)
    call h5dclose_f(dset, hdf_error)


    call h5fclose_f(file_id, hdf_error)

    bhat = cmplx(1.0,0.0) * re_buffer + cmplx(0.0,1.0) * im_buffer
    end if
    end subroutine
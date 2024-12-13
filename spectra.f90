! This subroutine is for computing 1D spectra on-the-fly during run-time

! For now, we only compute Exx(kx) (average taken in x and z), requiring only an fft in the x-direction
subroutine compute_1d_spectra
    use local_arrays, only: vx
    use modspec
    use param
    use mpih
    use mpi_param, only: kstart, kend
    implicit none
    real, allocatable, dimension(:) :: Exx_kx
    !real, allocatable, dimension(:,:,:) :: Exx

    integer i,j,k
    real :: dkx

    allocate(Exx_kx(1:n1m/2+1))
    !allocate(Exx(1:n1m/2+1,1:n2m, kstart:kend ) )

    dkx = 2.0 * pi / xlen

    Exx_kx = 0.0
    ! Exx = 0.0

    do j = 1,n2m
        do k = kstart, kend
            call dfftw_execute_dft_r2c(specplan, vx(1:n1m,j,k), uhat(:,j,k))
            !call dfftw_execute_dft(specplan, vx(:,j,k), uhat(:,j,k))

            uhat(:,j,k) = uhat(:,j,k) / dble(n1m)
            !Exx_kx(:) = Exx_kx(:) + ( real(uhat)**2 + aimag(uhat)**2   )
            !Exx(:,j,k) = ( real(uhat(:,j,k))**2 + aimag(uhat(:,j,k))**2   )

            Exx_kx = Exx_kx + ( real(uhat(:,j,k))**2 + aimag(uhat(:,j,k))**2   )

        enddo
    enddo

    ! do j = 1,n2m
    !     do k = kstart, kend
    !         Exx_kx(:) = Exx_kx(:) + Exx(:,j,k)
    !     enddo
    ! enddo

        ! do j = 1,n2m
        !     do k = kstart,kend
        !         !call dfftw_execute_dft(specplan, vx(:,j,k), uhat(:,j,k))
        !         !Exx_kx(i) = Exx_kx(i) + abs( uhat(i,j,k) )**2

        !         call dfftw_execute_dft(specplan, vx(:,j,k), uhat)
        !         uhat = uhat / dble(n1m)
        !         Exx_kx = Exx_kx + ( real(uhat)**2 + aimag(uhat)**2   )

        !     enddo
        ! enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE,Exx_kx,n1m/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    Exx_kx = Exx_kx / dble(n2m * n3m) / dkx

    ! Write to file
    call writeSpec(Exx_kx , n1m/2+1 )

    deallocate( Exx_kx )
    !deallocate( Exx )

end subroutine compute_1d_spectra 


! Initialisation of working memory, FFTW plans etc used in computing spectra
subroutine initSpectra
    use local_arrays, only: vx
    use modspec
    use param
    use mpi_param, only: kstart, kend
    implicit none
    integer :: my_n3m
    integer FFTW_EXHAUSTIVE
    parameter(FFTW_EXHAUSTIVE=64)
    
    my_n3m = kend - kstart + 1
    allocate(uhat(1:(n1m/2+1), 1:n2m, kstart:kend))
    !allocate( uhat(1:(n1m/2+1) ) )


    !call dfftw_plan_many_dft_r2c(specplan,1, n1m, n2m*my_n3m, vx(:,:,kstart:kend), [1, n2m*my_n3m], n2m*my_n3m, 1, uhat, [1, n2m*my_n3m], n2m*my_n3m, 1, FFTW_EXHAUSTIVE)
    !call dfftw_plan_many_dft_r2c(specplan,1, n1m, n2m*my_n3m, vx(:,:,kstart:kend), n1m, 1, n1m, uhat, [1, n2m*my_n3m], n2m*my_n3m, 1, FFTW_EXHAUSTIVE)



    call dfftw_plan_dft_r2c_1d(specplan,n1m,vx(1:n1m,1,kstart),uhat(:,1,kstart),FFTW_EXHAUSTIVE)

    !call dfftw_plan_dft_r2c_1d(specplan,n1m,vx(:,1,kstart),uhat(:),FFTW_EXHAUSTIVE)


end subroutine initSpectra

subroutine dealloc_spec
    use modspec
    implicit none
    
    call dfftw_destroy_plan(specplan)

    if(allocated(uhat)) deallocate(uhat)

end subroutine dealloc_spec


subroutine writeSpec(Exx_kx, size_kx)
    use param
    use mpih
    use mls_param
    use hdf5
   
    implicit none
  
    character(70) filename
    character(30) dataset
    character(30) dsetname
    character(5) ipfi
    integer :: itime
    real :: tprfi
    integer     ::  rank = 1     ! Dataset rank
    integer(hid_t) :: file_id       ! File identifier
    integer(hid_t) :: dset_id       ! Dataset identifier
    integer(hid_t) :: dspace_id     ! Dataspace identifier
    integer(hsize_t) :: dims(1)
    integer     ::  error        ! Error flag
    logical :: fileexists


    integer :: size_kx
    real, dimension(size_kx) :: Exx_kx

  
    tprfi = 1/tframe
    itime=nint(time*tprfi)
    write(ipfi,82)itime
    82 format(i5.5)

    filename ='spectra/Exx_kx_'//ipfi//'.h5'
     
  
    if (myid.eq.0) then 
      call hdf5_create_blank_file(filename)
      dsetname = trim("Exx_kx")
      !call hdf_write_1d(Exx_kx,n1m/2+1,dataset)

     !hdf_write_1d(dataset,d,dsetname)

      dims(1) = n1m/2+1
     
      call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
      call h5screate_simple_f(rank, dims, dspace_id, error)
    
      call h5lexists_f(file_id,dsetname,fileexists,error)
      if(fileexists) call h5ldelete_f(file_id,dsetname,error)
    
      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Exx_kx, dims, error)
    
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)
      call h5fclose_f(file_id, error)

    end if

    end subroutine writeSpec
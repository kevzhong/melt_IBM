      subroutine spec
      use param
      use mpih
      use hdf5
      use mpi_param
      use local_arrays
      use local_aux
      use AuxiliaryRoutines

      implicit none
      real :: vx_mean,kref
      real,allocatable,dimension(:) :: spec_tot,specx
      real,allocatable,dimension(:,:,:) :: v,rv1
      real,allocatable,dimension(:,:,:) :: specx_tot
      complex,allocatable,dimension(:,:,:) :: cv1,re_uhat,im_uhat
      integer :: j,k,i
      complex,allocatable,dimension(:) :: xb,xf
      real,allocatable,dimension(:) :: xr
      complex,allocatable,dimension(:,:) :: xa
      complex,allocatable,dimension(:,:,:) :: buf
      complex,allocatable,dimension(:,:,:) :: buft

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count
      integer(HSSIZE_T), dimension(3) :: data_offset

      integer :: comm, info
      integer :: ndims

      character*70 :: filnam1
      character*10 :: dsetname
      character(70) namfile

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      vx_mean=0.0d0
      allocate(rv1(n1m,1,1),v(n1m,n2m,n3m),cv1(n1m,n2m,n3m))
      allocate(specx_tot(1:n1mh,1:n2m,1:n3m))
      allocate(spec_tot(1:n1mh),specx(1:n1mh))

      pi = 2.*asin(1.)
      kref=2.*pi/xlen
    
      call AllocateReal1DArray(xr,1,n1m)
      call AllocateCplx1DArray(xb,1,n3m)
      call AllocateCplx1DArray(xf,1,n1m)
      call AllocateCplx2DArray(xa,1,n2mh,1,n1m)
      call AllocateCplx3DArray(buf,1,n1,1,n2,1,n3)
      call AllocateCplx3DArray(buft,1,n3,1,n1,1,n2)
 
          dsetname = trim('Vx')
          filnam1 = trim('continuation/field_01683.h5')

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1m
      dims(2)=n2m
      dims(3)=n3m


      data_count(1) = n1m
      data_count(2) = n2m
      data_count(3) = n3m


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dopen_f(file_id, dsetname, dset_qua, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,  &
         v(1:n1m,1:n2m,1:n3m), dims,              &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

     

!        do k=kstart,kend
!         do j=1,n2m
!            do i=1,n1m
!               vx_mean = vx_mean + vx(i,j,k)
!            end do
!         end do
!      end do

!     call MpiAllSumRealScalar(vx_mean)
!     vx_mean=vx_mean/float(n1m*n2m*n3m)
!    compute velocity fluctuation
     rv1=v!-vx_mean


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the comments are related to FFT in all three directions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      do k=1,n3m
!       xr(1:n1m,1:n2m) = rv1(1:n1m,1:n2m,1)
!       call dfftw_execute_dft_r2c(fwd_plan,xr,xa)
!       buf(1:n1mh,1:n2m,1) = xa(1:n1mh,1:n2m)
!      end do

!      call PackZ_UnpackR_C(buf(:,:,1:n3m),buft(:,:,1:n2m))
 
!      do j=jstart,jend
!!!!! keep it 1D but over the entire domain  !!!!!!!!!!!!!!!!!!!!!!!!!!
          do j=1,n2m
          do k=1,n3m
          xb(1:n1m) = rv1(1:n1m,j,k)!buft(1,1:n1m,1)
          xr(1:n1m) =xb(1:n1m)
!     end do
!     end do
          call dfftw_execute_dft(fwdplan_1d,xb,xb)
 !       buft(1:n3m,1,1) = xb(1:n3m)
!      end do
!      end do
          xb=xb/float(n1m)!*n2m*n3m)
!      call PackR_UnpackZ_C(buft(:,:,1:n2m),buf(:,:,1:n3m))
          specx=(real(xb))**2+(aimag(xb))**2
          specx=specx(1:n1mh)
   
          specx(1)=specx(1)/2
          specx(n1mh)=specx(n1mh)/2
          specx=2.*specx/kref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            do i=1,(2*n2m)
            specx_tot(1:n1mh,j,k)=specx
!            enddo

        enddo
      enddo
!!!!! keep it 1D but over the entire domain  !!!!!!!!!!!!!!!!!!!!!!!!!!
      spec_tot=0.0d0
      do k=1,n3m
      do j=1,n2m
      spec_tot(1:n1mh)=spec_tot(1:n1mh)+specx_tot(1:n1mh,j,k)
      enddo
      enddo

      spec_tot=spec_tot/(n3m*n2m)
    
      do i=1,n1mh
      if(ismaster) then
      namfile='flowmov/spec.txt'

      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)')xr(i),rv1(i,1,1), real(xb(i)),aimag(xb(i)), specx(i), kref,specx_tot(i,20,20),spec_tot(i)
      close(92)
      end if
      enddo
    


      deallocate(rv1,cv1,v,specx_tot)
      deallocate(spec_tot,specx)
      call DestroyReal1DArray(xr)
      call DestroyCplx1DArray(xb)
      call DestroyCplx1DArray(xf)
      call DestroyCplx2DArray(xa)
      call DestroyCplx3DArray(buf)
      call DestroyCplx3DArray(buft)
      return

      end subroutine spec

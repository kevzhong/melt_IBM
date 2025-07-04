      subroutine UpdateStats
      use param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr
      use local_aux, only: vorx, vory, vorz
      use stat_arrays
      use mpih
      implicit none
      real :: vz_rms_vol,vy_rms_vol,vx_rms_vol,vx_m,vy_m,vz_m
      real :: vorx_me,vory_me,vorz_me,vorx_rms_vol,vory_rms_vol,vorz_rms_vol
      real :: vorx_mes,vory_mes,vorz_mes,vorx_rms_svol,vory_rms_svol,vorz_rms_svol

      integer :: i,j,k

      timeint_cdsp = timeint_cdsp + 1

768   format(10e20.10)

!      pi = 2.*asin(1.)

      vx_rms_vol = 0.0d0
      vy_rms_vol = 0.0d0
      vz_rms_vol = 0.0d0
      vxvyvz_rms_vol = 0.0d0
      vorx_me = 0.0d0
      vory_me = 0.0d0
      vorz_me = 0.0d0

      vx_m = 0.0d0
      vy_m = 0.0d0
      vz_m = 0.0d0
      vorx_rms_vol = 0.0d0
      vory_rms_vol = 0.0d0
      vorz_rms_vol = 0.0d0
      do k=kstart,kend
         do j=1,n2m
            do i=1,n1m
               vx_m = vx_m + vx(i,j,k)*VOFx(i,j,k)
               vy_m = vy_m + vy(i,j,k)*VOFy(i,j,k)
               vz_m = vz_m + vz(i,j,k)*VOFz(i,j,k)
               pr_me(i,j,k) = pr_me(i,j,k) + pr(i,j,k)
               vorx_me = vorx_me + vorx(i,j,k)*VOFx(i,j,k)
               vory_me = vory_me + vory(i,j,k)*VOFy(i,j,k)
               vorz_me = vorz_me + vorz(i,j,k)*VOFz(i,j,k)

               vorx_rms_vol = vorx_rms_vol + (vorx(i,j,k)*VOFx(i,j,k))**2
               vory_rms_vol = vory_rms_vol + (vory(i,j,k)*VOFy(i,j,k))**2
               vorz_rms_vol = vorz_rms_vol + (vorz(i,j,k)*VOFz(i,j,k))**2

               vx_rms(i,j,k) = vx_rms(i,j,k) + (vx(i,j,k)*VOFx(i,j,k))**2
               vy_rms(i,j,k) = vy_rms(i,j,k) + (vy(i,j,k)*VOFy(i,j,k))**2
               vz_rms(i,j,k) = vz_rms(i,j,k) + (vz(i,j,k)*VOFz(i,j,k))**2
               pr_rms(i,j,k) = pr_rms(i,j,k) + pr(i,j,k)**2

               vx_rms_vol = vx_rms_vol + (vx(i,j,k)*VOFx(i,j,k))**2
               vy_rms_vol = vy_rms_vol + (vy(i,j,k)*VOFy(i,j,k))**2
               vz_rms_vol = vz_rms_vol + (vz(i,j,k)*VOFz(i,j,k))**2
               vxvyvz_rms_vol = vxvyvz_rms_vol + ( &
                (vx(i+1,j,k)*VOFx(i+1,j,k)+vx(i,j,k)*VOFx(i,j,k))**2+ &
                (vy(i,j+1,k)*VOFy(i,j+1,k)+vy(i,j,k)*VOFy(i,j,k))**2+ &
                (vz(i,j,k+1)*VOFz(i,j,k+1)+vz(i,j,k)*VOFz(i,j,k))**2)*0.25
            end do
         end do
      end do

      call MpiAllSumRealScalar(vx_rms_vol)
      call MpiAllSumRealScalar(vy_rms_vol)
      call MpiAllSumRealScalar(vz_rms_vol)
      call MpiAllSumRealScalar(vxvyvz_rms_vol)
      call MpiAllSumRealScalar(vx_m)
      call MpiAllSumRealScalar(vy_m)
      call MpiAllSumRealScalar(vz_m)
      call MpiAllSumRealScalar(vorx_me)
      call MpiAllSumRealScalar(vory_me)
      call MpiAllSumRealScalar(vorz_me)
      call MpiAllSumRealScalar(vorx_rms_vol)
      call MpiAllSumRealScalar(vory_rms_vol)
      call MpiAllSumRealScalar(vorz_rms_vol)


      vx_m=vx_m/float(n1m*n2m*n3m)
      vy_m=vy_m/float(n1m*n2m*n3m)
      vz_m=vz_m/float(n1m*n2m*n3m)

      vorx_me=vorx_me/float(n1m*n2m*n3m)
      vory_me=vory_me/float(n1m*n2m*n3m)
      vorz_me=vorz_me/float(n1m*n2m*n3m)

      vorx_rms_vol=sqrt(vorx_rms_vol/float(n1m*n2m*n3m))
      vory_rms_vol=sqrt(vory_rms_vol/float(n1m*n2m*n3m))
      vorz_rms_vol=sqrt(vorz_rms_vol/float(n1m*n2m*n3m))

      vx_rms_vol=sqrt(vx_rms_vol/float(n1m*n2m*n3m))
      vy_rms_vol=sqrt(vy_rms_vol/float(n1m*n2m*n3m))
      vz_rms_vol=sqrt(vz_rms_vol/float(n1m*n2m*n3m))
      vxvyvz_rms_vol=sqrt(vxvyvz_rms_vol/float(n1m*n2m*n3m))
        
       if(ismaster) then
       open(unit=43,file='flowmov/rms.txt',Access = 'append', Status='unknown')
       write(43,768) time,vx_m,vy_m,vz_m,vx_rms_vol,vy_rms_vol,vz_rms_vol, vxvyvz_rms_vol
       close(43)
       open(unit=43,file='flowmov/vor_rms.txt',Access = 'append', Status='unknown')
       write(43,768) time,vorx_me,vory_me,vorz_me,vorx_rms_vol,vory_rms_vol,vorz_rms_vol
       close(43)
       endif

      vorx_mes = 0.0d0
      vory_mes = 0.0d0
      vorz_mes = 0.0d0
      vorx_rms_svol = 0.0d0
      vory_rms_svol = 0.0d0
      vorz_rms_svol = 0.0d0

   do k=kstart,kend
   if(k.ge.(n3m/4).and.k.le.(3*n3m/4))then

         do j=n2m/4,3*n2m/4

            do i=n1m/4,3*n1m/4
               vorx_mes = vorx_mes + vorx(i,j,k)*VOFx(i,j,k)
               vory_mes = vory_mes + vory(i,j,k)*VOFy(i,j,k)
               vorz_mes = vorz_mes + vorz(i,j,k)*VOFz(i,j,k)

               vorx_rms_svol = vorx_rms_svol + (vorx(i,j,k)*VOFx(i,j,k))**2
               vory_rms_svol = vory_rms_svol + (vory(i,j,k)*VOFy(i,j,k))**2
               vorz_rms_svol = vorz_rms_svol + (vorz(i,j,k)*VOFz(i,j,k))**2

            end do
         end do
      end if
      end do
      vorx_mes=vorx_mes/float(n1m/2*n2m/2*n3m/2)
      vory_mes=vory_mes/float(n1m/2*n2m/2*n3m/2)
      vorz_mes=vorz_mes/float(n1m/2*n2m/2*n3m/2)
      vorx_rms_svol=sqrt(vorx_rms_svol/float(n1m/2*n2m/2*n3m/2))
      vory_rms_svol=sqrt(vory_rms_svol/float(n1m/2*n2m/2*n3m/2))
      vorz_rms_svol=sqrt(vorz_rms_svol/float(n1m/2*n2m/2*n3m/2))
       if(ismaster) then
       open(unit=43,file='flowmov/vor_rms_svol.txt',Access = 'append', Status='unknown')
       write(43,768) time,vorx_mes,vory_mes,vorz_mes,vorx_rms_svol,vory_rms_svol,vorz_rms_svol
       close(43)
       endif



      end
!    
!***********************************************************************
      subroutine WriteStats
      use mpih
      use param
      use mpi_param, only: kstart,kend
      use stat_arrays
      use hdf5

      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_vxme
      integer(HID_T) :: dset_vyme
      integer(HID_T) :: dset_vzme
      integer(HID_T) :: dset_prme

      integer(HID_T) :: dset_vxrms
      integer(HID_T) :: dset_vyrms
      integer(HID_T) :: dset_vzrms
      integer(HID_T) :: dset_prrms

      integer(HSIZE_T) :: dims(3)

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character(40)filnam1
      character(40)filnamgrid

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam1 = 'continuation/stafield_data.h5'

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1m
      dims(2)=n2m
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims,  &
                              filespace, hdf_error)

      data_count(1) = n1m
      data_count(2) = n2m
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!RO   Open file and create dataspace

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, & 
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, &
        hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

!RO   Create the datasets with default properties.

      call h5dcreate_f(file_id, 'Vx_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vxme, hdf_error)

      call h5dcreate_f(file_id, 'Vy_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vyme, hdf_error)

      call h5dcreate_f(file_id, 'Vz_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vzme, hdf_error)

      call h5dcreate_f(file_id, 'pr_mean', H5T_NATIVE_DOUBLE, &
                      filespace, dset_prme, hdf_error)

      call h5dcreate_f(file_id, 'Vx_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vxrms, hdf_error)

      call h5dcreate_f(file_id, 'Vy_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vyrms, hdf_error)

      call h5dcreate_f(file_id, 'Vz_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vzrms, hdf_error)

      call h5dcreate_f(file_id, 'pr_rms', H5T_NATIVE_DOUBLE, &
                      filespace, dset_prrms, hdf_error)



      call h5sclose_f(filespace, hdf_error)

!RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 


!RO   Select hyperslab  and then write it

!RO   Q1me


      call h5dget_space_f(dset_vxme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vxme, H5T_NATIVE_DOUBLE, &
         vx_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2me

      call h5dget_space_f(dset_vyme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vyme, H5T_NATIVE_DOUBLE, &
         vy_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3me

      call h5dget_space_f(dset_vzme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vzme, H5T_NATIVE_DOUBLE, &
         vz_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)


!EP   prme

      call h5dget_space_f(dset_prme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_prme, H5T_NATIVE_DOUBLE, &
         pr_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q1rms

      call h5dget_space_f(dset_vxrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vxrms, H5T_NATIVE_DOUBLE, &
         vx_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2rms

      call h5dget_space_f(dset_vyrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vyrms, H5T_NATIVE_DOUBLE, &
         vy_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3rms

      call h5dget_space_f(dset_vzrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vzrms, H5T_NATIVE_DOUBLE, &
         vz_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!     prrms

      call h5dget_space_f(dset_prrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_prrms, H5T_NATIVE_DOUBLE, &
         pr_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)



!RO   Close properties and datasets

      call h5dclose_f(dset_vxme, hdf_error)
      call h5dclose_f(dset_vyme, hdf_error)
      call h5dclose_f(dset_vzme, hdf_error)
      call h5dclose_f(dset_prme, hdf_error)

      call h5dclose_f(dset_vxrms, hdf_error)
      call h5dclose_f(dset_vyrms, hdf_error)
      call h5dclose_f(dset_vzrms, hdf_error)
      call h5dclose_f(dset_prrms, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!RO   Write the grid & statistics information
!RO   only if master process

      if (myid.eq.0) then

      ndims=1

      filnamgrid = 'continuation/stafield_master.h5'
      call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_id, hdf_error)

!RO   Write amount of averages 

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'averaging_time', H5T_NATIVE_INTEGER, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Write Reynolds number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Re', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ren, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)
           

!EP   Write Prandtl number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, pra, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)
           


!RO   Write the grid information 

      dims_grid(1)=n1m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'X_Cordin', H5T_NATIVE_DOUBLE,&
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xm(1:n1m),&
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Y_cordin', H5T_NATIVE_DOUBLE,&
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ym(1:n2m),&
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_id, 'Z_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
              dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Close file

      call h5fclose_f(file_id, hdf_error)

      endif

      return  
      end

!***********************************************************************

      subroutine InitStats
      use param
      use mpi_param, only: kstart,kend
      use stat_arrays
      use mpih
      use hdf5
      implicit none
      integer :: i,j,k

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: memspace
      integer(HID_T) :: slabspace

      integer(HID_T) :: dset_vxme
      integer(HID_T) :: dset_vyme
      integer(HID_T) :: dset_vzme
      integer(HID_T) :: dset_prme

      integer(HID_T) :: dset_vxrms
      integer(HID_T) :: dset_vyrms
      integer(HID_T) :: dset_vzrms
      integer(HID_T) :: dset_prrms

      integer(HSIZE_T) :: dims(3)

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_id
      integer(HID_T) :: plist_full
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character(40) filnam1
      character(40) filnamgrid

!EP   Read or initialize stat arrays

      if(starea.eq.1) then

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam1 = 'continuation/stafield_data.h5'
      filnamgrid = 'continuation/stafield_master.h5'

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1m
      dims(2)=n2m
      dims(3)=n3m

      data_count(1)=n1m
      data_count(2)=n2m
      data_count(3)=kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1
         
!RO   Open file and create dataspace

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_full, hdf_error)
      call h5pset_fapl_mpio_f(plist_full, comm, info, hdf_error)

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, hdf_error, &
                       access_prp=plist_full)


!RO   Create the datasets with default properties.

      call h5dopen_f(file_id, 'Vx_mean', dset_vxme, hdf_error)
      call h5dopen_f(file_id, 'Vy_mean', dset_vyme, hdf_error)
      call h5dopen_f(file_id, 'Vz_mean', dset_vzme, hdf_error)
      call h5dopen_f(file_id, 'pr_mean', dset_prme, hdf_error)

      call h5dopen_f(file_id, 'Vx_rms', dset_vxrms, hdf_error)
      call h5dopen_f(file_id, 'Vy_rms', dset_vyrms, hdf_error)
      call h5dopen_f(file_id, 'Vz_rms', dset_vzrms, hdf_error)
      call h5dopen_f(file_id, 'pr_rms', dset_prrms, hdf_error)


!RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!RO   Select hyperslab  and then read it

!RO   Q1me

      call h5dget_space_f(dset_vxme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vxme, H5T_NATIVE_DOUBLE, &
         vx_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2me

      call h5dget_space_f(dset_vyme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vyme, H5T_NATIVE_DOUBLE, &
         vy_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3me

      call h5dget_space_f(dset_vzme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vzme, H5T_NATIVE_DOUBLE, &
         vz_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   prme

      call h5dget_space_f(dset_prme, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_prme, H5T_NATIVE_DOUBLE, &
         pr_me(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q1rms

      call h5dget_space_f(dset_vxrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vxrms, H5T_NATIVE_DOUBLE, &
         vz_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q2rms

      call h5dget_space_f(dset_vyrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vyrms, H5T_NATIVE_DOUBLE, &
         vy_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   Q3rms

      call h5dget_space_f(dset_vzrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_vzrms, H5T_NATIVE_DOUBLE, &
         vz_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

!EP   PRrms

      call h5dget_space_f(dset_prrms, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_prrms, H5T_NATIVE_DOUBLE, &
         pr_rms(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)


!RO   Close properties and datasets and file

      call h5dclose_f(dset_vxme, hdf_error)
      call h5dclose_f(dset_vyme, hdf_error)
      call h5dclose_f(dset_vzme, hdf_error)
      call h5dclose_f(dset_prme, hdf_error)

      call h5dclose_f(dset_vxrms, hdf_error)
      call h5dclose_f(dset_vyrms, hdf_error)
      call h5dclose_f(dset_vzrms, hdf_error)
      call h5dclose_f(dset_prrms, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


!RO   Read the grid & statistics information

      ndims=1
      dims_grid(1)=1

      call h5fopen_f(filnamgrid, H5F_ACC_RDONLY_F, file_id, hdf_error, &
                       access_prp=plist_full)

      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dopen_f(file_id, 'averaging_time', dset_grid, hdf_error)

      call h5dread_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5pclose_f(plist_full, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP starea.eq.0
      else
      
      timeint_cdsp = 0

!EP   Initialize to 0

      do k=kstart,kend
      do j=1,n2m
      do i=1,n1m
       vx_me(i,j,k)    =0.0d0
       vy_me(i,j,k)    =0.0d0
       vz_me(i,j,k)    =0.0d0
       pr_me(i,j,k)    =0.0d0
       vx_rms(i,j,k)   =0.0d0
       vy_rms(i,j,k)   =0.0d0
       vz_rms(i,j,k)   =0.0d0
       pr_rms(i,j,k)   =0.0d0
      enddo
      enddo
      enddo

      endif

      return
      end

      subroutine mkmov_hdf_xcut
      use param
      use hdf5
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr,forcx,forcy,forcz

      IMPLICIT none

      integer ic,jc,kc
      integer ip,jp,kp

      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      integer(HID_T) :: dset_prv

      integer(HSIZE_T) :: dims(2)


      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      real tprfi
      real prx(n2m,n3m),v1(n2m,n3m),v2(n2m,n3m),v3(n2m,n3m)
      integer itime

      character(70) namfile,xdmnam
      character(5) ipfi


      ndims=2
      ic = pind1(1,1) 
      ip=ic+1
      
      do kc=kstart,kend
       kp=kc+1 

       do jc=1,n2m                                                     
        jp=jpv(jc)

        v1(jc,kc) = (vx(ic,jc,kc)+vx(ip,jc,kc))*0.5
        v2(jc,kc) = (vy(ic,jc,kc)+vy(ic,jp,kc))*0.5
        v3(jc,kc) = (vz(ic,jc,kc)+vz(ic,jc,kp))*0.5
        prx(jc,kc)= pr(ic,jc,kc)
       end do
      end do

!RO   File writing part

!RO   Form the name of the file

      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)

      namfile='flowmov/frame_x_'//ipfi//'.h5'
      xdmnam='flowmov/frame_x_'//ipfi//'.xmf'

!RO   Sort out MPI definitions and open file

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

!RO   Create dataspace

      dims(1)=n2m
      dims(2)=n3m
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
         

!RO   Create the dataset with default properties.

      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, filespace, dset_q1v, hdf_error)

      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, filespace, dset_q2v, hdf_error)

      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, filespace, dset_q3v, hdf_error)

      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, filespace, dset_prv, hdf_error)


!RO   Set offsets and element counts

      data_count(1)=n2m
      data_count(2)=kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1

!RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!RO   Select hyperslab  and then write it

      call h5dget_space_f(dset_q1v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_q1v, H5T_NATIVE_DOUBLE, &
        v1(1:n2m,kstart:kend), dims, &
        hdf_error, file_space_id = filespace, mem_space_id = memspace, &
        xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dget_space_f(dset_q2v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_q2v, H5T_NATIVE_DOUBLE, &
         v2(1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dget_space_f(dset_q3v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

       call h5dwrite_f(dset_q3v, H5T_NATIVE_DOUBLE, &
         v3(1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

      call h5dget_space_f(dset_prv, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_prv, H5T_NATIVE_DOUBLE, & 
         prx(1:n2m,kstart:kend), dims, & 
         hdf_error, file_space_id = filespace, mem_space_id = memspace, & 
         xfer_prp = plist_id)
!RO   Close properties and file

      call h5dclose_f(dset_q1v, hdf_error)
      call h5dclose_f(dset_q2v, hdf_error)
      call h5dclose_f(dset_q3v, hdf_error)
      call h5dclose_f(dset_prv, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      if (myid.eq.0) then

      open(45,file=xdmnam,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""thetacut"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""2DCORECTMESH"" NumberOfElements=""",i4," ",i4,"""/>")') n3m,n2m
      write(45,'("<Geometry GeometryType=""ORIGIN_DXDY"">")')
      write(45,'("<DataItem Name=""Origin"" Dimensions=""2"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 0,0.0
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Name=""Spacing"" Dimensions=""2"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 1./dx3,1./dx2
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""X-velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_x_",i5.5,".h5:/Vx")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Y-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_x_",i5.5,".h5:/Vy")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Z-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_x_",i5.5,".h5:/Vz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Pressure"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_x_",i5.5,".h5:/Pr")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,"""/>")')time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      end if


      return
      end


      subroutine mkmov_hdf_ycut
      use param
      use hdf5
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr,temp
      use local_aux, only: vorx,vory,vorz

      IMPLICIT none

      integer ic,jc,kc
      integer ip,jp,kp
      integer inp
      real, allocatable, dimension(:,:) :: prx,v1,v2,v3,tempx
      real, allocatable, dimension(:,:) :: vor1,vor2,vor3,hel
      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      integer(HID_T) :: dset_prv
      integer(HID_T) :: dset_tempv
      integer(HID_T) :: dset_vor1
      integer(HID_T) :: dset_vor2
      integer(HID_T) :: dset_vor3
      integer(HID_T) :: dset_hel
 
      integer(HSIZE_T) :: dims(2)


      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 
      integer :: comm, info
      integer :: ndims
      
      real tprfi
      integer itime

      character(70) namfile,xdmnam
      character(5) ipfi

      allocate(prx(n1m,n3m),v1(n1m,n3m),v2(n1m,n3m),v3(n1m,n3m), tempx(n1m,n3m))
      allocate(vor1(n1m,n3m),vor2(n1m,n3m),vor3(n1m,n3m),hel(n1m,n3m))

      ndims=2
!      if((imlsfor.eq.1).and.(Nparticle.eq.1))then
!            jc = pind1(2,1)
!      else
      jc = n2m/2
!      end if
      jp=jc+1
      do kc=kstart,kend
         kp=kc+1

       do ic=1,n1m
        ip=ipv(ic)

        v1(ic,kc) = (vx(ic,jc,kc)+vx(ip,jc,kc))*0.5
        v2(ic,kc) = (vy(ic,jc,kc)+vy(ic,jp,kc))*0.5
        v3(ic,kc) = (vz(ic,jc,kc)+vz(ic,jc,kp))*0.5
        prx(ic,kc)= pr(ic,jc,kc)
        tempx(ic,kc) = temp(ic,jc,kc)
        vor1(ic,kc) = (vorx(ic,jc,kc)+vorx(ip,jc,kc))*0.5
        vor2(ic,kc) = (vory(ic,jc,kc)+vory(ic,jp,kc))*0.5
        vor3(ic,kc) = (vorz(ic,jc,kc)+vorz(ic,jc,kp))*0.5
        hel(ic,kc)  = (vx(ic,jc,kc)*vorx(ic,jc,kc)+ &
                       vx(ip,jc,kc)*vorx(ip,jc,kc)+ &
                       vy(ic,jc,kc)*vory(ic,jc,kc)+ &
                       vy(ic,jp,kc)*vory(ic,jp,kc)+ &
                       vz(ic,jc,kc)*vorz(ic,jc,kc)+ &
                       vz(ic,jc,kp)*vorz(ic,jc,kp))*0.5d0
        end do
      end do

!RO   File writing part

!RO   Form the name of the file

      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)

      namfile='flowmov/frame_y_'//ipfi//'.h5'
      xdmnam='flowmov/frame_y_'//ipfi//'.xmf'

!RO   Sort out MPI definitions and open file

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

!RO   Create dataspace

      dims(1)=n1m
      dims(2)=n3m
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
         

!RO   Create the dataset with default properties.

      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, filespace, dset_q1v, hdf_error)

      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, filespace, dset_q2v, hdf_error)

      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, filespace, dset_q3v, hdf_error)

      call h5dcreate_f(file_id, 'temperature', H5T_NATIVE_DOUBLE, filespace, dset_tempv, hdf_error)

      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, filespace, dset_prv, hdf_error)

      call h5dcreate_f(file_id, 'vorx', H5T_NATIVE_DOUBLE, filespace, dset_vor1, hdf_error)

      call h5dcreate_f(file_id, 'vory', H5T_NATIVE_DOUBLE, filespace, dset_vor2, hdf_error)

      call h5dcreate_f(file_id, 'vorz', H5T_NATIVE_DOUBLE, filespace, dset_vor3, hdf_error)

      call h5dcreate_f(file_id, 'hel', H5T_NATIVE_DOUBLE, filespace, dset_hel, hdf_error)


!RO   Set offsets and element counts

      data_count(1)=n1m
      data_count(2)=kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1

!RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!RO   Select hyperslab  and then write it

      call h5dget_space_f(dset_q1v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_q1v, H5T_NATIVE_DOUBLE, &
        v1(1:n1m,kstart:kend), dims, &
        hdf_error, file_space_id = filespace, mem_space_id = memspace, &
        xfer_prp = plist_id)

      call h5dget_space_f(dset_q2v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_q2v, H5T_NATIVE_DOUBLE, &
         v2(1:n1m,kstart:kend), dims, &
         hdf_error, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

      call h5dget_space_f(dset_q3v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

       call h5dwrite_f(dset_q3v, H5T_NATIVE_DOUBLE, &
         v3(1:n1m,kstart:kend), dims, &
         hdf_error, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

      call h5dget_space_f(dset_prv, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_prv, H5T_NATIVE_DOUBLE, & 
         prx(1:n1m,kstart:kend), dims, & 
         hdf_error, file_space_id = filespace, mem_space_id = memspace, & 
         xfer_prp = plist_id)

      ! TEMPERATURE
         call h5dget_space_f(dset_tempv, filespace, hdf_error)
         call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
   
         call h5dwrite_f(dset_tempv, H5T_NATIVE_DOUBLE, & 
            tempx(1:n1m,kstart:kend), dims, & 
            hdf_error, file_space_id = filespace, mem_space_id = memspace, & 
            xfer_prp = plist_id)

      ! END TEMPERATURE

      call h5dget_space_f(dset_vor1, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_vor1, H5T_NATIVE_DOUBLE, &
        vor1(1:n1m,kstart:kend), dims, &
        hdf_error, file_space_id = filespace, mem_space_id = memspace, &
        xfer_prp = plist_id)

      call h5dget_space_f(dset_vor2, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dwrite_f(dset_vor2, H5T_NATIVE_DOUBLE, &
         vor2(1:n1m,kstart:kend), dims, &
         hdf_error, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

      call h5dget_space_f(dset_vor3, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

       call h5dwrite_f(dset_vor3, H5T_NATIVE_DOUBLE, &
         vor3(1:n1m,kstart:kend), dims, &
         hdf_error, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

      call h5dget_space_f(dset_hel, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

       call h5dwrite_f(dset_hel, H5T_NATIVE_DOUBLE, &
         hel(1:n1m,kstart:kend), dims, &
         hdf_error, file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)
!RO   Close properties and file

      call h5dclose_f(dset_q1v, hdf_error)
      call h5dclose_f(dset_q2v, hdf_error)
      call h5dclose_f(dset_q3v, hdf_error)
      call h5dclose_f(dset_tempv, hdf_error)
      call h5dclose_f(dset_prv, hdf_error)
      call h5dclose_f(dset_vor1, hdf_error)
      call h5dclose_f(dset_vor2, hdf_error)
      call h5dclose_f(dset_vor3, hdf_error)
      call h5dclose_f(dset_hel, hdf_error)


      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      if (myid.eq.0) then

      open(45,file=xdmnam,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""thetacut"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""2DCORECTMESH"" NumberOfElements=""",i4," ",i4,"""/>")') n3m,n1m
      write(45,'("<Geometry GeometryType=""ORIGIN_DXDY"">")')
      write(45,'("<DataItem Name=""Origin"" Dimensions=""2"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 0.0,0.0
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Name=""Spacing"" Dimensions=""2"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 1./dx3,1./dx1

      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""X-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/Vx")') itime

      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Y-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/Vy")') itime

      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Z-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/Vz")') itime

      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Temperature"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/temperature")') itime

      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Pressure"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/Pr")') itime

      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""vorx"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/vorx")') itime

      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""vory"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/vory")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""vorz"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/vorz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""hel"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n1m
      write(45,'("frame_y_",i5.5,".h5:/hel")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,"""/>")')time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      end if
      deallocate(prx,v1,v2,v3,tempx)
      deallocate(vor1,vor2,vor3,hel)
      return
      end


      subroutine mkmov_hdf_zcut

      use param
      use hdf5
      use mpih
      use mls_param
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr
      use local_aux, only: vorx,vory,vorz

      IMPLICIT none

      integer ic,jc,kc
      integer ip,jp,kp
      integer inp
      real, allocatable, dimension(:,:) :: prx,v1,v2,v3
      real, allocatable, dimension(:,:) :: vor1,vor2,vor3
      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      integer(HID_T) :: dset_prv
      integer(HID_T) :: dset_vor1
      integer(HID_T) :: dset_vor2
      integer(HID_T) :: dset_vor3
 

      integer(HSIZE_T) :: dims(2)


      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      real tprfi
   !   real prx(n1m,n2m),v1(n1m,n2m),v2(n1m,n2m),v3(n1m,n2m)
      integer itime

      character(70) namfile,xdmnam
      character(5) ipfi
     
      allocate(prx(n1m,n2m),v1(n1m,n2m),v2(n1m,n2m),v3(n1m,n2m))
      allocate(vor1(n1m,n2m),vor2(n1m,n2m),vor3(n1m,n2m))
      
      ndims=2
      if(imlsfor.eq.1)then

            kc = pind1(3,1)

      else
      kc = n3m/2
      end if
      kp=kc+1
      if (kc.ge.kstart .and. kc.le.kend) then
      
      do jc= 1,n2m
        jp=jpv(jc)

       do ic=1,n1m                                                     
        ip=ipv(ic)

        v1(ic,jc) = (vx(ic,jc,kc)+vx(ip,jc,kc))*0.5
        v2(ic,jc) = (vy(ic,jc,kc)+vy(ic,jp,kc))*0.5
        v3(ic,jc) = (vz(ic,jc,kc)+vz(ic,jc,kp))*0.5
        prx(ic,jc)= pr(ic,jc,kc)
        vor1(ic,jc) = (vorx(ic,jc,kc)+vorx(ip,jc,kc))*0.5
        vor2(ic,jc) = (vory(ic,jc,kc)+vory(ic,jp,kc))*0.5
        vor3(ic,jc) = (vorz(ic,jc,kc)+vorz(ic,jc,kp))*0.5

       end do
      end do

!RO   File writing part

!RO   Form the name of the file

      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)

      namfile='flowmov/frame_z_'//ipfi//'.h5'
      xdmnam='flowmov/frame_z_'//ipfi//'.xmf'

!RO   Sort out MPI definitions and open file

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      dims(1)=n1m
      dims(2)=n2m

      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)


!RO   Create the dataset with default properties.

      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, filespace, &
                      dset_q1v, hdf_error)
      call h5dwrite_f(dset_q1v, H5T_NATIVE_DOUBLE, v1, dims, hdf_error)


      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, filespace, &
                      dset_q2v, hdf_error)
      call h5dwrite_f(dset_q2v, H5T_NATIVE_DOUBLE, v2, dims, hdf_error)


      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, filespace, &
                      dset_q3v, hdf_error)
      call h5dwrite_f(dset_q3v, H5T_NATIVE_DOUBLE, v3, dims, hdf_error)

      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, filespace, &
                      dset_prv, hdf_error)
      call h5dwrite_f(dset_prv, H5T_NATIVE_DOUBLE, prx, dims, hdf_error)

      call h5dcreate_f(file_id, 'vorx', H5T_NATIVE_DOUBLE, filespace, &
                      dset_vor1, hdf_error)
      call h5dwrite_f(dset_vor1, H5T_NATIVE_DOUBLE, vor1, dims, hdf_error)


      call h5dcreate_f(file_id, 'vory', H5T_NATIVE_DOUBLE, filespace, &
                      dset_vor2, hdf_error)
      call h5dwrite_f(dset_vor2, H5T_NATIVE_DOUBLE, vor2, dims, hdf_error)


      call h5dcreate_f(file_id, 'vorz', H5T_NATIVE_DOUBLE, filespace, &
                      dset_vor3, hdf_error)
      call h5dwrite_f(dset_vor3, H5T_NATIVE_DOUBLE, vor3, dims, hdf_error)


!RO   Close properties and file

      call h5dclose_f(dset_q1v, hdf_error)
      call h5dclose_f(dset_q2v, hdf_error)
      call h5dclose_f(dset_q3v, hdf_error)
      call h5dclose_f(dset_prv, hdf_error)
      call h5dclose_f(dset_vor1, hdf_error)
      call h5dclose_f(dset_vor2, hdf_error)
      call h5dclose_f(dset_vor3, hdf_error)
      
      call h5sclose_f(filespace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      open(45,file=xdmnam,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""thetacut"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""2DCORECTMESH"" NumberOfElements=""",i4," ",i4,"""/>")') n1m,n2m
      write(45,'("<Geometry GeometryType=""ORIGIN_DXDY"">")')
      write(45,'("<DataItem Name=""Origin"" Dimensions=""2"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 0,0.0
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Name=""Spacing"" Dimensions=""2"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 1./dx1,1./dx2
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""X-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m,n2m
      write(45,'("frame_z_",i5.5,".h5:/Vx")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Y-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m,n2m
      write(45,'("frame_z_",i5.5,".h5:/Vy")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Z-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m,n2m
      write(45,'("frame_z_",i5.5,".h5:/Vz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Pressure"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m,n2m
      write(45,'("frame_z_",i5.5,".h5:/Pr")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""vorx"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m,n2m
      write(45,'("frame_z_",i5.5,".h5:/vorx")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""vory"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m,n2m
      write(45,'("frame_z_",i5.5,".h5:/vory")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""vorz"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m,n2m
      write(45,'("frame_z_",i5.5,".h5:/vorz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,"""/>")')time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      end if


      return
      end



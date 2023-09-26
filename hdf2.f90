subroutine hdf5_create_blank_file(filename)
use hdf5
implicit none
character(70),intent(in) :: filename
integer(HID_T) :: file_id
integer :: error


call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)
call h5fclose_f(file_id, error)

end subroutine

subroutine hdf_write_1d(dataset,d,dsetname)
use hdf5

implicit none

character(70) filename
character(30) dsetname

integer(hid_t) :: file_id       ! File identifier

integer(hid_t) :: dset_id       ! Dataset identifier
integer(hid_t) :: dspace_id     ! Dataspace identifier

integer(hsize_t) :: dims(1)
integer          :: d
real        :: dataset(d)
integer     ::  rank = 1     ! Dataset rank
integer     ::  error        ! Error flag
logical :: fileexists


filename = 'continuation/particles.h5'
dsetname = trim(dsetname)

  dims(1) = d
 

  call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
  call h5screate_simple_f(rank, dims, dspace_id, error)

  call h5lexists_f(file_id,dsetname,fileexists,error)
  if(fileexists) call h5ldelete_f(file_id,dsetname,error)

  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dims, error)

  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  call h5fclose_f(file_id, error)

end subroutine

subroutine hdf_write_2d(dataset,d,dsetname)
use hdf5

implicit none

character(70) filename
character(30) dsetname

integer(hid_t) :: file_id       ! File identifier

integer(hid_t) :: dset_id       ! Dataset identifier
integer(hid_t) :: dspace_id     ! Dataspace identifier

integer(hsize_t) :: dims(2)
integer          :: d(2)
real,dimension(d(1),d(2))  :: dataset
integer     ::  rank = 2     ! Dataset rank
integer     ::  error        ! Error flag
logical :: fileexists



filename = 'continuation/particles.h5'
dsetname = trim(dsetname)
!allocate( dataset( dims(1),dims(2),dims(3) ) )

  dims(1:2) = d(1:2)
 
  call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
  call h5screate_simple_f(rank, dims, dspace_id, error)

  call h5lexists_f(file_id,dsetname,fileexists,error)
  if(fileexists) call h5ldelete_f(file_id,dsetname,error)

  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dims, error)

  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  call h5fclose_f(file_id, error)

end subroutine

subroutine hdf_write_3d(dataset,d,dsetname)
use hdf5

implicit none

character(70) filename
character(30) dsetname

integer(hid_t) :: file_id       ! File identifier

integer(hid_t) :: dset_id       ! Dataset identifier
integer(hid_t) :: dspace_id     ! Dataspace identifier

integer(hsize_t) :: dims(3)
integer          :: d(3)
real,dimension(d(1),d(2),d(3))  :: dataset
integer     ::  rank = 3     ! Dataset rank
integer     ::  error        ! Error flag
logical :: fileexists



filename = 'continuation/particles.h5'
dsetname = trim(dsetname)

  ! NEED rank, dims, dset, dsetname 
  dims(1:3) = d(1:3)
 
  call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
  call h5screate_simple_f(rank, dims, dspace_id, error)

  call h5lexists_f(file_id,dsetname,fileexists,error)
  if(fileexists) call h5ldelete_f(file_id,dsetname,error)

  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dims, error)

  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  call h5fclose_f(file_id, error)

!deallocate(dataset)

end subroutine

subroutine hdf_read_1d(var,d,dsetname)
use hdf5
use mpih
implicit none
character(200) :: dsetname,filename
integer, intent(in) :: d
real, dimension(d), intent(out) :: var
integer(HID_T) :: file_id
integer(HID_T) :: dset
integer :: hdf_error
integer(HSIZE_T) :: dims(1)

filename = 'continuation/particles.h5'


if (myid.eq.0) then
call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf_error)

dims(1)=d

call h5dopen_f(file_id, dsetname, dset, hdf_error)

call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
   var, dims, hdf_error)

call h5dclose_f(dset, hdf_error)

call h5fclose_f(file_id, hdf_error)

end if
end subroutine

subroutine hdf_read_2d(var,d,dsetname)
use hdf5
use mpih
implicit none
character(200) :: dsetname,filename
integer, intent(in) :: d(2)
real, dimension(d), intent(out) :: var(d(1),d(2))
integer(HID_T) :: file_id
integer(HID_T) :: dset
integer :: hdf_error
integer(HSIZE_T) :: dims(2)

filename = 'continuation/particles.h5'


if (myid.eq.0) then
call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf_error)

dims(1)=d(1)
dims(2)=d(2)

call h5dopen_f(file_id, dsetname, dset, hdf_error)

call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
   var, dims, hdf_error)

call h5dclose_f(dset, hdf_error)

call h5fclose_f(file_id, hdf_error)

end if
end subroutine

subroutine hdf_read_3d(var,d,dsetname)
use hdf5
use mpih
implicit none
character(200) :: dsetname,filename
integer, intent(in) :: d(3)
real, dimension(d), intent(out) :: var(d(1),d(2),d(2))
integer(HID_T) :: file_id
integer(HID_T) :: dset
integer :: hdf_error
integer(HSIZE_T) :: dims(3)

filename = 'continuation/particles.h5'


if (myid.eq.0) then
call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf_error)

dims(1)=d(1)
dims(2)=d(2)
dims(3)=d(3)

call h5dopen_f(file_id, dsetname, dset, hdf_error)

call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
   var, dims, hdf_error)

call h5dclose_f(dset, hdf_error)

call h5fclose_f(file_id, hdf_error)

end if
end subroutine


subroutine init_quat_flowmov
  use param
  use mls_param
  use hdf5
  use mpih

  implicit none

  character(len=16), parameter :: filename   = "vtkfiles/vtk.h5"
  integer :: hdferr

  integer(hid_t) :: file_id

  logical :: file_exists

  if (myid.eq.0) then

  inquire(file=filename, exist=file_exists)

  if (file_exists .eq. .false. .or. nread.eq.0) then

  call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, hdferr)
  call h5fclose_f(file_id, hdferr)

  endif

  endif

end subroutine

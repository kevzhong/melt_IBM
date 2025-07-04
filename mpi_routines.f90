!===============================================
      subroutine block(n, p, irank, istart, iend, blcsz)
      implicit none
      integer,intent(in) :: n,p,irank
      integer,intent(out) :: istart,iend
      integer :: i
      integer,dimension(0:p-1),intent(out) :: blcsz
      
      do i=0,p-1
      blcsz(i) = floor(real((n+p-i-1)/p))
      enddo
      istart = sum(blcsz(0:irank))-blcsz(irank)+1
      iend = istart+blcsz(irank)-1

      end subroutine block

!=================================================           
      subroutine mpi_workdistribution
      use param
      use mpih 
      use mpi_param
      implicit none
      integer :: i
      
      if(.not. allocated(countj)) allocate(countj(0:numtasks-1))
      if(.not. allocated(countjp)) allocate(countjp(0:numtasks-1))
      if(.not. allocated(countk)) allocate(countk(0:numtasks-1))

!EP   For PERIODIC pressure solver
      call block(n2+1, numtasks, myid, jstartp, jendp, countjp)
      djp=jendp-jstartp+1

      call block(n2m, numtasks, myid, jstart, jend, countj)
      dj=jend-jstart+1
      
      call block(n3m, numtasks, myid, kstart, kend, countk)
      dk=kend-kstart+1

      !KZ FOR DEBUG
      !write(*,*) "myid:  ",myid, "kstart is: ", kstart, "kend is: ", kend

#ifdef DEBUG
      write(*,*) "myid:  ",myid, "(jstart,jend): ", jstart, jend, "(kstart,kend): ", kstart, kend
      !write(*,*) "jstart: ",jstart
      !write(*,*) "jend: ",jend
      !write(*,*) "jstartp: ",jstart
      !write(*,*) "jendp: ",jend
      !write(*,*) "kstart: ",kstart
      !write(*,*) "kend: ",kend
#endif

      if( dj .lt. 1 ) then            
       write(6,*)'process ',myid,' has work load <1 cell in j direction'
       write(6,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif

      if( dk .lt. 1 ) then            
       write(6,*)'process ',myid,' has work load <1 cell in k direction'
       write(6,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
  
      if(.not. allocated(offsetjp)) allocate(offsetjp(0:numtasks-1))
      if(.not. allocated(offsetj)) allocate(offsetj(0:numtasks-1))
      if(.not. allocated(offsetk)) allocate(offsetk(0:numtasks-1))
      
      offsetjp(:)=0
      offsetj(:)=0
      offsetk(:)=0
      do i=1,numtasks-1
        offsetjp(i)= offsetjp(i-1) + countjp(i-1)
        offsetj(i)= offsetj(i-1) + countj(i-1)
        offsetk(i)= offsetk(i-1) + countk(i-1)
      end do
      
      !-------For MPI-IO--------------------------------
      mydata= n2*dk*n1
      mydatam = n2m*dk*n1m

      if(myid .eq. numtasks-1) mydata = n2*(dk+1)*n1
      
      if(.not. allocated(countf)) allocate(countf(0:numtasks-1))
      if(.not. allocated(offsetf)) allocate(offsetf(0:numtasks-1))
       
      call MPI_ALLGATHER(mydata, 1, MPI_INTEGER, countf, 1, MPI_INTEGER, &
       MPI_COMM_WORLD,ierr)
    
      offsetf(:)=0
      do i=1,numtasks-1
        offsetf(i)= offsetf(i-1) + countf(i-1)
      end do
      
      !------------------------------------------------
      
      end subroutine mpi_workdistribution 
!===============================================
      subroutine update_both_ghosts(n1,n2,q,ks,ke)
      use mpih
      implicit none
      integer, intent(in) :: ks,ke
      real,intent(inout) :: q(n1,n2,ks-lvlhalo:ke+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      
      mydata= n1*n2*lvlhalo
      
      my_down=myid-1
      
      my_up=myid+1

      if(myid .eq. 0) my_down=numtasks-1
      if(myid .eq. numtasks-1) my_up=0

      tag=1
      call MPI_ISEND(q(1,1,ke-lvlhalo+1), mydata, MDP, &
       my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      
      call MPI_ISEND(q(1,1,ks), mydata,  MDP, &
       my_down,tag,MPI_COMM_WORLD,req(2), ierr)
     
      call MPI_IRECV(q(1,1,ks-lvlhalo), mydata,  MDP,  &
       my_down,tag,MPI_COMM_WORLD,req(3),ierr)
     
      call MPI_IRECV(q(1,1,ke+1), mydata,  MDP, &
       my_up, tag,MPI_COMM_WORLD,req(4),ierr)
     
      call MPI_Waitall(4,req,status,ierr)

      end subroutine update_both_ghosts
!=========================================
      subroutine update_upper_ghost(n1,n2,q)
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout) :: q(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down=numtasks-1
      if(myid .eq. numtasks-1) my_up= 0
     
      tag=1
      
      call MPI_ISEND(q(1,1,kstart), mydata, MDP, &
       my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(q(1,1,kend+1), mydata, MDP, &
       my_up,tag, MPI_COMM_WORLD, req(2), ierr)
       
      call MPI_Waitall(2,req,status,ierr)

     
      end subroutine update_upper_ghost
!=========================================
      subroutine update_lower_ghost(n1,n2,q)
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout) :: q(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down= numtasks-1
      if(myid .eq. numtasks-1) my_up= 0
      
      tag=1
      
      call MPI_ISEND(q(1,1,kend-lvlhalo+1), mydata,  MDP, &
       my_up, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(q(1,1,kstart-lvlhalo), mydata,  MDP, &
       my_down,tag, MPI_COMM_WORLD, req(2), ierr)
       
      call MPI_Waitall(2,req,status,ierr)

      end subroutine update_lower_ghost

!=========================================
subroutine update_add_lower_ghost(q1)
  use param, only: n1, n2
  use mpi_param, only: kstart,kend, buf_n1n2
  use mpih
  implicit none
  real,intent(inout) :: q1(n1,n2,kstart-1:kend+1)
  !real :: buf(n1,n2)
  !real, allocatable :: buf(:,:)

  integer :: mydata
  integer :: my_down, my_up,tag
  integer :: kc,ii

  mydata= n1*n2

  !allocate(buf(n1,n2))    ! <<--- now on the heap

  
  my_down= myid-1
  
  my_up= myid+1

  buf_n1n2=0.0d0

  if(myid .eq. 0) my_down= numtasks-1
  if(myid .eq. numtasks-1) my_up= 0
 
      
  tag=1
  
  call MPI_ISEND(q1(1,1,kend+1),mydata,MDP, my_up, tag, MPI_COMM_WORLD, req(1), ierr)

  call MPI_IRECV(buf_n1n2(1,1), mydata, MDP, my_down,tag, MPI_COMM_WORLD, req(2), ierr)

  call MPI_Waitall(2,req,status,ierr)

  kc = kstart
  do ii = 1,n2
    q1(:,ii,kc) = q1(:,ii,kc) + buf_n1n2(:,ii)
  end do

     !deallocate(buf)

end subroutine update_add_lower_ghost

subroutine update_add_upper_ghost(q1)
  use param, only: n1, n2
  use mpi_param, only: kstart,kend,buf_n1n2
  use mpih
  implicit none
  real,intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo)
  !real :: buf(n1,n2)
  !real, allocatable :: buf(:,:)
  integer :: mydata
  integer :: my_down, my_up,tag
  integer :: kc,ii


  !write(*,*) "proc", myid, "entering upperGhost routine"
  mydata= n1*n2
  
  !allocate(buf(n1,n2))    ! <<--- now on the heap

  my_down= myid-1
  
  my_up= myid+1

  buf_n1n2=0.0d0

  if(myid .eq. 0) my_down= numtasks-1
  if(myid .eq. numtasks-1) my_up= 0

   tag=1

   call MPI_ISEND(q1(1,1,kstart-1),mydata,MDP, my_down, tag, MPI_COMM_WORLD, req(1), ierr)
   
   call MPI_IRECV(buf_n1n2(1,1), mydata, MDP, my_up,tag, MPI_COMM_WORLD, req(2), ierr)

   call MPI_Waitall(2,req,status,ierr)

   !write(*,*) "myid starting loop", myid
   kc=kend
   do ii=1,n2
     q1(:,ii,kc) = q1(:,ii,kc) + buf_n1n2(:,ii)
   enddo

   !deallocate(buf)

end subroutine update_add_upper_ghost
!
subroutine get_prow_pcol
      ! KZ: auxilary pencil routine for pencil-accelerated ray-tagging
      use mpih
      implicit none
      integer :: i, factor1, factor2
      real :: aspect_ratio, best_aspect_ratio
  
      ! First, find the best choice of p_row * p_col given numtasks
  
      ! Initialize best aspect ratio to a large value
      best_aspect_ratio = 1.0e30
      p_row = 1
      p_col = numtasks
  
      ! Loop over possible factors of numtasks
      do i = 1, int(sqrt(real(numtasks, 8)))
          if (mod(numtasks, i) == 0) then
              ! i is a factor, so numtasks / i is the corresponding pair
              factor1 = i
              factor2 = numtasks / i
  
              ! Compute the aspect ratio and find the best one
              aspect_ratio = abs(real(factor1, 8) / real(factor2, 8) - 1.0)
              if (aspect_ratio < best_aspect_ratio) then
                  best_aspect_ratio = aspect_ratio
                  p_row = factor1
                  p_col = factor2
              end if
          end if
      end do
  
      my_p_row  = myid / p_col
      my_p_col = mod(myid, p_col)
  
      !write(*,*) "Total prow, pcol: ", p_row, p_col
  
  end subroutine get_prow_pcol

!==============================================




      !   subroutine mpi_globalsum_double_arr(var,nvar)
      !   use mpih
      !     implicit none
      !     real,intent(inout),dimension(nvar) :: var
      !     real,dimension(nvar) :: var2
      !     integer,intent(in) :: nvar

      !     call MPI_ALLREDUCE(var,var2,nvar,MPI_DOUBLE_PRECISION, &
      !         MPI_SUM,MPI_COMM_WORLD,ierr)          
      !     var = var2

      !   end subroutine mpi_globalsum_double_arr

subroutine mpi_globalsum_double_var(var)
use mpih
  implicit none
  real,intent(inout) :: var
  real :: var2

  call MPI_ALLREDUCE([var],[var2],1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,ierr)
  var = var2

end subroutine mpi_globalsum_double_var
!==============================================
        subroutine mpi_globalsum_double_forc(var)
        use mpih
          implicit none
          real,intent(inout),dimension(6,7,7,7) :: var
          real,dimension(6,7,7,7) :: var2
!          integer :: nvar=1715
          integer :: nvar=2058!6*7*7*7

          call MPI_ALLREDUCE(var,var2,nvar,MPI_DOUBLE_PRECISION, &
              MPI_SUM,MPI_COMM_WORLD,ierr)          
          var = var2

        end subroutine mpi_globalsum_double_forc

!==================================================      
      
      subroutine mem_dealloc
      use local_arrays
      use local_aux
      use mpi_param
      use stat_arrays
      use mls_param
      use mls_local
      implicit none
      
      if(allocated(vx)) deallocate(vx)
      if(allocated(vy)) deallocate(vy)
      if(allocated(vz)) deallocate(vz)
      if(allocated(temp)) deallocate(temp)

      if(allocated(qcap)) deallocate(qcap)
      
      if(allocated(pr)) deallocate(pr)
      
      if(allocated(rhs)) deallocate(rhs)
      
      if(allocated(dph)) deallocate(dph)
      
      if(allocated(ru1)) deallocate(ru1)
      if(allocated(ru2)) deallocate(ru2)
      if(allocated(ru3)) deallocate(ru3)
      
      if(allocated(vorx)) deallocate(vorx)
      if(allocated(vory)) deallocate(vory)
      if(allocated(vorz)) deallocate(vorz)
      
      if(allocated(diss)) deallocate(diss)
      if(allocated(tke)) deallocate(tke)
      if(allocated(chi)) deallocate(chi)

      !---------------------------------------
      if(allocated(countj)) deallocate(countj)
      if(allocated(countk)) deallocate(countk)

      if(allocated(offsetj)) deallocate(offsetj)
      if(allocated(offsetk)) deallocate(offsetk)
      
      if(allocated(countf)) deallocate(countf)
      
      if(allocated(offsetf)) deallocate(offsetf)
      
      if(allocated(vx_me)) deallocate(vx_me)
      if(allocated(vy_me)) deallocate(vy_me)
      if(allocated(vz_me)) deallocate(vz_me)
      if(allocated(vx_rms)) deallocate(vx_rms)
      if(allocated(vy_rms)) deallocate(vy_rms)
      if(allocated(vz_rms)) deallocate(vz_rms)
      if(allocated(pr_rms)) deallocate(pr_rms)
      if(allocated(pr_me)) deallocate(pr_me)


      !------------------------------------------
        ! KZ HIT stuff
      deallocate( waveN )
      deallocate( bhat )
      deallocate( exp_I_kl_xi, exp_I_km_yj, exp_I_kn_zk)
      deallocate( exp_I_kl_xsi, exp_I_km_ysj, exp_I_kn_zsk)


      ! KZ: IBM stuff
      ! Cell-tagging masks
      if(allocated(VOFx)) deallocate(VOFx)
      if(allocated(VOFy)) deallocate(VOFy)
      if(allocated(VOFz)) deallocate(VOFz)
      if(allocated(VOFp)) deallocate(VOFp)
      !if(allocated(d_UsolidT_dxj)) deallocate(d_UsolidT_dxj)
      !if(allocated(solid_mask)) deallocate(solid_mask)
      if(allocated(sdf)) deallocate(sdf)

      ! HIT forcing field
      if(allocated(forcx)) deallocate(forcx)
      if(allocated(forcy)) deallocate(forcy)
      if(allocated(forcz)) deallocate(forcz)

      ! IBM force
      if(allocated(for_xc)) deallocate(for_xc)
      if(allocated(for_yc)) deallocate(for_yc)
      if(allocated(for_zc)) deallocate(for_zc)
      if(allocated(for_temp)) deallocate(for_temp)


      ! MPI buffer
      if(allocated(buf_n1n2)) deallocate(buf_n1n2)

      end subroutine mem_dealloc
!================================================
!       subroutine mpi_write_continua
!       use param
!       use mpih
!       use mpi_param, only: kstart,kend
!       use local_arrays, only: vy,vz,vx,pr,temp
!       use hdf5
!       implicit none

!       integer hdf_error

!       integer(HID_T) :: file_id
!       integer(HID_T) :: filespace
!       integer(HID_T) :: slabspace
!       integer(HID_T) :: memspace

!       integer(HID_T) :: dset_vx
!       integer(HID_T) :: dset_vy
!       integer(HID_T) :: dset_vz
!       integer(HID_T) :: dset_pr
!       integer(HID_T) :: dset_temp
!       integer(HID_T) :: dset_enst

!       integer(HSIZE_T) :: dims(3)

!       integer(HID_T) :: plist_id
!       integer(HSIZE_T), dimension(3) :: data_count  
!       integer(HSSIZE_T), dimension(3) :: data_offset 

!       integer(HSIZE_T) :: dims_grid(1)
!       integer(HID_T) :: dset_grid
!       integer(HID_T) :: dspace_grid

!       integer :: comm, info
!       integer :: ndims

!       character(40) filnam2,filnam3,filnam4
!       character(40) filnamgrid,filnam5,filnam6

! !RO   Sort out MPI definitions

!       comm = MPI_COMM_WORLD
!       info = MPI_INFO_NULL

! !RO   Form the name of the file

!       filnam2 = 'continuation/continua_vx.h5'
!       filnam3 = 'continuation/continua_vy.h5'
!       filnam4 = 'continuation/continua_vz.h5'
!       filnam5 = 'continuation/continua_pr.h5'
!       filnam6 = 'continuation/continua_temp.h5'

! !RO   Set offsets and element counts
   
!       ndims = 3

!       dims(1)=n1
!       dims(2)=n2
!       dims(3)=n3m

!       call h5screate_simple_f(ndims, dims, filespace, hdf_error)

!       data_count(1) = n1
!       data_count(2) = n2
!       data_count(3) = kend-kstart+1

!       data_offset(1) = 0
!       data_offset(2) = 0
!       data_offset(3) = kstart-1

! !RO   pressure

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
!           hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, &
!         hdf_error)

!       call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'pr', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_pr, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_pr, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
!          pr(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_pr, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! !EP   vx

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vx, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_vx, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
!          vx(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_vx, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! !EP   vy

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vy, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_vy, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
!          vy(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_vy, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! !EP   vz

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vz, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_vz, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
!          vz(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_vz, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! ! Temperature

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam6, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Temp', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_temp, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_temp, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_temp, H5T_NATIVE_DOUBLE, &
!          temp(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_temp, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)
 
!       if (myid .eq. 0) then
!        open(13,file='continuation/continua_grid.dat',status='unknown')
!        rewind(13)                                                      
!        write(13,*) n1,n2,n3,time
!        close(13)
!       endif
      
! !RO   Write the grid & statistics information
! !RO   only if master process

!       if (myid.eq.0) then

!       ndims=1

!       filnamgrid = 'continuation/continua_master.h5'
!       call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_id, hdf_error)

! !RO   Write Reynolds number

!       dims_grid(1)=1
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'Re', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ren, &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)
           

! !EP   Write Prandtl number

!       dims_grid(1)=1
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, pra, &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)
           

! !RO   Write the grid information 

!       dims_grid(1)=n1
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'X_cordin', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xc(1:n1), &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)

!       dims_grid(1)=n2
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'Y_cordin', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, yc(1:n2), &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)

!       dims_grid(1)=n3m
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
!       call h5dcreate_f(file_id, 'Z_cordin', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
!               dims_grid, hdf_error)


!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)

! !RO   Close file

!       call h5fclose_f(file_id, hdf_error)

!       endif
!       end subroutine mpi_write_continua

! ! THIS VERSION CORRECTS FOR LARGE-MEMORY BUFFERS
! subroutine mpi_write_continua
!       use param
!       use mpih
!       use mpi_param, only: kstart,kend
!       use local_arrays, only: vy,vz,vx,pr,temp
!       use hdf5
!       implicit none

!       integer hdf_error

!       integer(HID_T) :: file_id
!       integer(HID_T) :: filespace
!       integer(HID_T) :: slabspace
!       integer(HID_T) :: memspace

!       integer(HID_T) :: dset_vx
!       integer(HID_T) :: dset_vy
!       integer(HID_T) :: dset_vz
!       integer(HID_T) :: dset_pr
!       integer(HID_T) :: dset_temp
!       integer(HID_T) :: dset_enst

!       integer(HSIZE_T) :: dims(3)

!       integer(HID_T) :: plist_id
!       integer(HSIZE_T), dimension(3) :: data_count  
!       integer(HSSIZE_T), dimension(3) :: data_offset 

!       integer(HSIZE_T) :: dims_grid(1)
!       integer(HID_T) :: dset_grid
!       integer(HID_T) :: dspace_grid

!       integer :: comm, info
!       integer :: ndims

!       character(40) filnam2,filnam3,filnam4
!       character(40) filnamgrid,filnam5,filnam6

! !RO   Sort out MPI definitions

!       comm = MPI_COMM_WORLD
!       info = MPI_INFO_NULL

! !RO   Form the name of the file

!       filnam2 = 'continuation/continua_vx.h5'
!       filnam3 = 'continuation/continua_vy.h5'
!       filnam4 = 'continuation/continua_vz.h5'
!       filnam5 = 'continuation/continua_pr.h5'
!       filnam6 = 'continuation/continua_temp.h5'

! !RO   Set offsets and element counts
   
!       ndims = 3

!       !dims(1)=n1
!       !dims(2)=n2
!       !dims(3)=n3m

!       data_count(1) = n1
!       data_count(2) = n2
!       data_count(3) = kend-kstart+1

!       data_offset(1) = 0
!       data_offset(2) = 0
!       data_offset(3) = kstart-1

!       dims = data_count ! KZ local slice size

!       call h5screate_simple_f(ndims, dims, filespace, hdf_error)



! !RO   pressure

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
!           hdf_error)


!       call h5pset_fapl_mpio_f(plist_id, comm, info, &
!         hdf_error)

!       call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)


!       call h5pclose_f(plist_id, hdf_error)


!       call h5dcreate_f(file_id, 'pr', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_pr, hdf_error)


!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_pr, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!       !  call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
!       !    pr(1:n1,1:n2,kstart:kend), dims,  &
!       !    hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!       !    xfer_prp = plist_id)
!       call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
!          pr(1,1,kstart), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_pr, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! !EP   vx

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vx, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_vx, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!       !  call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
!       !    vx(1:n1,1:n2,kstart:kend), dims,  &
!       !    hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!       !    xfer_prp = plist_id)
!       call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
!          vx(1,1,kstart), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_vx, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! !EP   vy

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vy, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_vy, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!       !  call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
!       !    vy(1:n1,1:n2,kstart:kend), dims,  &
!       !    hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!       !    xfer_prp = plist_id)
!       call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
!          vy(1,1,kstart), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_vy, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! !EP   vz

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vz, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_vz, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!       !  call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
!       !    vz(1:n1,1:n2,kstart:kend), dims,  &
!       !    hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!       !    xfer_prp = plist_id)
!       call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
!          vz(1,1,kstart), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_vz, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)

! ! Temperature

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam6, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Temp', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_temp, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!       call h5dget_space_f(dset_temp, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!       !  call h5dwrite_f(dset_temp, H5T_NATIVE_DOUBLE, &
!       !    temp(1:n1,1:n2,kstart:kend), dims,  &
!       !    hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!       !    xfer_prp = plist_id)
!       call h5dwrite_f(dset_temp, H5T_NATIVE_DOUBLE, &
!          temp(1,1,kstart), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_temp, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)
 
!       if (myid .eq. 0) then
!        open(13,file='continuation/continua_grid.dat',status='unknown')
!        rewind(13)                                                      
!        write(13,*) n1,n2,n3,time
!        close(13)
!       endif
      
! !RO   Write the grid & statistics information
! !RO   only if master process

!       if (myid.eq.0) then

!       ndims=1

!       filnamgrid = 'continuation/continua_master.h5'
!       call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_id, hdf_error)

! !RO   Write Reynolds number

!       dims_grid(1)=1
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'Re', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ren, &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)
           

! !EP   Write Prandtl number

!       dims_grid(1)=1
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, pra, &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)
           

! !RO   Write the grid information 

!       dims_grid(1)=n1
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'X_cordin', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xc(1:n1), &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)

!       dims_grid(1)=n2
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

!       call h5dcreate_f(file_id, 'Y_cordin', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, yc(1:n2), &
!              dims_grid,hdf_error)

!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)

!       dims_grid(1)=n3m
!       call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
!       call h5dcreate_f(file_id, 'Z_cordin', H5T_NATIVE_DOUBLE, &
!                       dspace_grid, dset_grid, hdf_error)

!       call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
!               dims_grid, hdf_error)


!       call h5dclose_f(dset_grid, hdf_error)
!       call h5sclose_f(dspace_grid, hdf_error)

! !RO   Close file

!       call h5fclose_f(file_id, hdf_error)

!       endif
!       end subroutine mpi_write_continua


subroutine mpi_write_continua
  use param
  use mpih
  use mpi_param, only: kstart, kend
  use local_arrays, only: vx, vy, vz, pr, temp
  use hdf5
  implicit none

  integer            :: hdf_error, comm, info, ndims
  integer(HID_T)     :: file_id, filespace, slabspace, memspace
  integer(HID_T)     :: dset_vx, dset_vy, dset_vz, dset_pr, dset_temp
  integer(HID_T)     :: plist_id
  integer(HSIZE_T)   :: file_dims(3), mem_dims(3)
  integer(HSIZE_T)   :: data_count(3)
  integer(HSSIZE_T)  :: data_offset(3)

  character(40)      :: filnam_vx, filnam_vy, filnam_vz
  character(40)      :: filnam_pr, filnam_temp
  character(40)      :: filnamgrid

  integer(HSIZE_T) :: dims_grid(1)
  integer(HID_T) :: dset_grid
  integer(HID_T) :: dspace_grid

  ! MPI‐IO setup
  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  ! file names
  filnam_vx   = 'continuation/continua_vx.h5'
  filnam_vy   = 'continuation/continua_vy.h5'
  filnam_vz   = 'continuation/continua_vz.h5'
  filnam_pr   = 'continuation/continua_pr.h5'
  filnam_temp = 'continuation/continua_temp.h5'

  ndims = 3

  ! full on-disk dimensions 
  file_dims = (/ n1, n2, n3m /)

  ! local slab dimensions in k
  data_count  = (/ n1, n2, kend - kstart + 1 /)
  data_offset = (/ 0,  0,  kstart - 1       /)

  ! memory dims = local slab size
  mem_dims = data_count

!   ! Optional debug print
!   if (myid .eq. 0) then
!     write(*,*) ' file_dims(3)=', file_dims(3),        &
!                ' offset=',   data_offset(3),         &
!                ' count=',    data_count(3),          &
!                ' sum=',      data_offset(3)+data_count(3)
!   end if

  ! ————————————————
  ! Pressure dataset
  ! ————————————————
  ! create file dataspace
  call h5screate_simple_f(ndims, file_dims, filespace, hdf_error)

  ! setup MPI‐IO FAPL and open file
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
  call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
  call h5fcreate_f(filnam_pr, H5F_ACC_TRUNC_F, file_id, hdf_error, &
                   access_prp=plist_id)
  call h5pclose_f(plist_id, hdf_error)

  ! create the dataset
  call h5dcreate_f(file_id, 'pr', H5T_NATIVE_DOUBLE, filespace, dset_pr, hdf_error)

  ! create memory dataspace
  call h5screate_simple_f(ndims, mem_dims, memspace, hdf_error)

  ! select slab on file
  call h5dget_space_f(dset_pr, slabspace, hdf_error)
  call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)

  ! collective‐IO setup
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

  ! **write only the local slab** into the file hyperslab
  call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE,        &
                  pr(1,1,kstart),                    & ! buffer start
                  mem_dims,                          & ! local slab shape
                  hdf_error,                         &
                  file_space_id=slabspace,           &
                  mem_space_id=memspace,             &
                  xfer_prp=plist_id)
  call h5pclose_f(plist_id, hdf_error)

  call h5dclose_f(dset_pr, hdf_error)
  call h5sclose_f(memspace, hdf_error)
  call h5fclose_f(file_id, hdf_error)
  call h5sclose_f(filespace, hdf_error)

  ! — Repeat the above block for vx, vy, vz, temp —  
  ! simply replace filnam_pr→filnam_vx, dset_pr→dset_vx, 
  ! pr→vx, etc., and you get parallel writes of each dataset.

  ! … (VX) …
  call h5screate_simple_f(ndims, file_dims, filespace, hdf_error)
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
  call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
  call h5fcreate_f(filnam_vx, H5F_ACC_TRUNC_F, file_id, hdf_error, &
                   access_prp=plist_id)
  call h5pclose_f(plist_id, hdf_error)
  call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, filespace, dset_vx, hdf_error)
  call h5screate_simple_f(ndims, mem_dims,   memspace,   hdf_error)
  call h5dget_space_f(dset_vx, slabspace, hdf_error)
  call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
  call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, vx(1,1,kstart), mem_dims, hdf_error, &
                  file_space_id=slabspace, mem_space_id=memspace, xfer_prp=plist_id)
  call h5pclose_f(plist_id, hdf_error)
  call h5dclose_f(dset_vx, hdf_error)
  call h5sclose_f(memspace, hdf_error)
  call h5fclose_f(file_id, hdf_error)
  call h5sclose_f(filespace, hdf_error)

!———————————————————————————————————————————————————————————————
! Vy
!———————————————————————————————————————————————————————————————
call h5screate_simple_f(ndims, file_dims,  filespace, hdf_error)
call h5pcreate_f    (H5P_FILE_ACCESS_F, plist_id, hdf_error)
call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
call h5fcreate_f    (filnam_vy, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
call h5pclose_f     (plist_id, hdf_error)

call h5dcreate_f    (file_id, 'Vy', H5T_NATIVE_DOUBLE, filespace, dset_vy, hdf_error)
call h5screate_simple_f(ndims, mem_dims, memspace, hdf_error)

call h5dget_space_f (dset_vy, slabspace, hdf_error)
call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)

call h5pcreate_f    (H5P_DATASET_XFER_F, plist_id, hdf_error)
call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

call h5dwrite_f     (dset_vy, H5T_NATIVE_DOUBLE, &
                     vy(1,1,kstart), mem_dims, hdf_error, &
                     file_space_id=slabspace, mem_space_id=memspace, xfer_prp=plist_id)

call h5pclose_f     (plist_id, hdf_error)
call h5dclose_f     (dset_vy, hdf_error)
call h5sclose_f     (memspace, hdf_error)
call h5fclose_f     (file_id, hdf_error)
call h5sclose_f     (filespace, hdf_error)


!———————————————————————————————————————————————————————————————
! Vz
!———————————————————————————————————————————————————————————————
call h5screate_simple_f(ndims, file_dims,  filespace, hdf_error)
call h5pcreate_f    (H5P_FILE_ACCESS_F, plist_id, hdf_error)
call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
call h5fcreate_f    (filnam_vz, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
call h5pclose_f     (plist_id, hdf_error)

call h5dcreate_f    (file_id, 'Vz', H5T_NATIVE_DOUBLE, filespace, dset_vz, hdf_error)
call h5screate_simple_f(ndims, mem_dims, memspace, hdf_error)

call h5dget_space_f (dset_vz, slabspace, hdf_error)
call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)

call h5pcreate_f    (H5P_DATASET_XFER_F, plist_id, hdf_error)
call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

call h5dwrite_f     (dset_vz, H5T_NATIVE_DOUBLE, &
                     vz(1,1,kstart), mem_dims, hdf_error, &
                     file_space_id=slabspace, mem_space_id=memspace, xfer_prp=plist_id)

call h5pclose_f     (plist_id, hdf_error)
call h5dclose_f     (dset_vz, hdf_error)
call h5sclose_f     (memspace, hdf_error)
call h5fclose_f     (file_id, hdf_error)
call h5sclose_f     (filespace, hdf_error)


!———————————————————————————————————————————————————————————————
! Temp
!———————————————————————————————————————————————————————————————
call h5screate_simple_f(ndims, file_dims,  filespace, hdf_error)
call h5pcreate_f    (H5P_FILE_ACCESS_F, plist_id, hdf_error)
call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
call h5fcreate_f    (filnam_temp, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
call h5pclose_f     (plist_id, hdf_error)

call h5dcreate_f    (file_id, 'Temp', H5T_NATIVE_DOUBLE, filespace, dset_temp, hdf_error)
call h5screate_simple_f(ndims, mem_dims, memspace, hdf_error)

call h5dget_space_f (dset_temp, slabspace, hdf_error)
call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)

call h5pcreate_f    (H5P_DATASET_XFER_F, plist_id, hdf_error)
call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

call h5dwrite_f     (dset_temp, H5T_NATIVE_DOUBLE, &
                     temp(1,1,kstart), mem_dims, hdf_error, &
                     file_space_id=slabspace, mem_space_id=memspace, xfer_prp=plist_id)

call h5pclose_f     (plist_id, hdf_error)
call h5dclose_f     (dset_temp, hdf_error)
call h5sclose_f     (memspace, hdf_error)
call h5fclose_f     (file_id, hdf_error)
call h5sclose_f     (filespace, hdf_error)

      if (myid .eq. 0) then
       open(13,file='continuation/continua_grid.dat',status='unknown')
       rewind(13)                                                      
       write(13,*) n1,n2,n3,time
       close(13)
      endif
      
!RO   Write the grid & statistics information
!RO   only if master process

      if (myid.eq.0) then

      ndims=1

      filnamgrid = 'continuation/continua_master.h5'
      call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_id, hdf_error)

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

      dims_grid(1)=n1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'X_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xc(1:n1), &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Y_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, yc(1:n2), &
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

end subroutine mpi_write_continua

      
!================================================


    subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
      use mpih
      use param
      use hdf5
      implicit none
      integer, intent(in) :: ks,ke,n2o,n1o,n3o
      real, dimension(1:n1o,1:n2o,ks-1:ke+1)::qua

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

      integer, intent(in) :: intvar
      character*70 :: filnam1
      character*10 :: dsetname

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!EP   Select file and dataset based on intvar

      select case (intvar)
        case (1)
          dsetname = trim('Vx')
          filnam1 = trim('continuation/continua_vx.h5')
        case (2)
          dsetname = trim('Vy')
          filnam1 = trim('continuation/continua_vy.h5')
        case (3)
          dsetname = trim('Vz')
          filnam1 = trim('continuation/continua_vz.h5')
        case (4)
          dsetname = trim('pr')
          filnam1 = trim('continuation/continua_pr.h5')
        case (5)
          dsetname = trim('Temp')
          filnam1 = trim('continuation/continua_temp.h5')
      end select

!RO   Set offsets and element counts
   
      ndims = 3

      data_count(1) = n1o
      data_count(2) = n2o
      data_count(3) = ke-ks+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = ks-1

      dims = data_count ! KZ local slice size



!     call h5open_f(hdf_error)

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
      !  call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,  &
      !    qua(1:n1o,1:n2o,ks:ke), dims,              &
      !    hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
      !    xfer_prp = plist_id)

             call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,  &
         qua(1,1,ks), dims,              &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
         xfer_prp = plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
!     call h5close_f(hdf_error)

      if(myid.eq.0)write(*,'(5x,a)')'reading complete: '//filnam1

      end subroutine mpi_read_continua
!================================================


!     subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
!       use mpih
!       use param
!       use hdf5
!       implicit none
!       integer, intent(in) :: ks,ke,n2o,n1o,n3o
!       real, dimension(1:n1o,1:n2o,ks-1:ke+1)::qua

!       integer hdf_error

!       integer(HID_T) :: file_id
!       integer(HID_T) :: slabspace
!       integer(HID_T) :: memspace

!       integer(HID_T) :: dset_qua

!       integer(HSIZE_T) :: dims(3)

!       integer(HID_T) :: plist_id
!       integer(HSIZE_T), dimension(3) :: data_count
!       integer(HSSIZE_T), dimension(3) :: data_offset

!       integer :: comm, info
!       integer :: ndims

!       integer, intent(in) :: intvar
!       character*70 :: filnam1
!       character*10 :: dsetname

!       comm = MPI_COMM_WORLD
!       info = MPI_INFO_NULL

! !EP   Select file and dataset based on intvar

!       select case (intvar)
!         case (1)
!           dsetname = trim('Vx')
!           filnam1 = trim('continuation/continua_vx.h5')
!         case (2)
!           dsetname = trim('Vy')
!           filnam1 = trim('continuation/continua_vy.h5')
!         case (3)
!           dsetname = trim('Vz')
!           filnam1 = trim('continuation/continua_vz.h5')
!         case (4)
!           dsetname = trim('pr')
!           filnam1 = trim('continuation/continua_pr.h5')
!         case (5)
!           dsetname = trim('Temp')
!           filnam1 = trim('continuation/continua_temp.h5')
!       end select

! !RO   Set offsets and element counts
   
!       ndims = 3

!       dims(1)=n1o
!       dims(2)=n2o
!       dims(3)=n3o-1


!       data_count(1) = n1o
!       data_count(2) = n2o
!       data_count(3) = ke-ks+1

!       data_offset(1) = 0
!       data_offset(2) = 0
!       data_offset(3) = ks-1



! !     call h5open_f(hdf_error)

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, hdf_error, access_prp=plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dopen_f(file_id, dsetname, dset_qua, hdf_error)
!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

!       call h5dget_space_f(dset_qua, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
!        call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,  &
!          qua(1:n1o,1:n2o,ks:ke), dims,              &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
!          xfer_prp = plist_id)
!       call h5pclose_f(plist_id, hdf_error)

!       call h5dclose_f(dset_qua, hdf_error)

!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)
! !     call h5close_f(hdf_error)

!       if(myid.eq.0)write(*,'(5x,a)')'reading complete: '//filnam1

!       end subroutine mpi_read_continua
! !================================================

      subroutine mpi_write_field
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: vy,vz,vx,pr,temp
      use hdf5
      implicit none
      integer ic,jc,kc
      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_vx
      integer(HID_T) :: dset_vy
      integer(HID_T) :: dset_vz
      integer(HID_T) :: dset_pr
      integer(HID_T) :: dset_temp
      integer(HID_T) :: dset_VOFx
      integer(HID_T) :: dset_VOFy
      integer(HID_T) :: dset_VOFz
      integer(HID_T) :: dset_VOFp
      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count
      integer(HSSIZE_T), dimension(3) :: data_offset

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer :: comm, info
      integer :: ndims
      real tprfi
      integer itime
      character(70) filnam2,filnam3,filnam4,xdmnam
      character(70) filnamgrid,filnam5,filnam6
      character(7) ipfi
      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i7.7)

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

!       filnam2 = 'continuation/field.h5'
       filnam2 = 'continuation/field_'//ipfi//'.h5'
       xdmnam  = 'continuation/field_'//ipfi//'.xmf'

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m
 

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1


!EP   vx

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vx, hdf_error)
      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vy, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vz, hdf_error)
      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
                      filespace, dset_pr, hdf_error)
      call h5dcreate_f(file_id, 'Temp', H5T_NATIVE_DOUBLE, &
      filespace, dset_temp, hdf_error)
      call h5dcreate_f(file_id, 'VOFx', H5T_NATIVE_DOUBLE, &
                      filespace, dset_VOFx, hdf_error)
      call h5dcreate_f(file_id, 'VOFy', H5T_NATIVE_DOUBLE, &
                      filespace, dset_VOFy, hdf_error)
      call h5dcreate_f(file_id, 'VOFz', H5T_NATIVE_DOUBLE, &
                      filespace, dset_VOFz, hdf_error)
      call h5dcreate_f(file_id, 'VOFp', H5T_NATIVE_DOUBLE, &
      filespace, dset_VOFp, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
! vx
      call h5dget_space_f(dset_vx, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
         vx(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
! vy
      call h5dget_space_f(dset_vy, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
         vy(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
!vz
      call h5dget_space_f(dset_vz, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
         vz(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
!pr
      call h5dget_space_f(dset_pr, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
         pr(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
!temp
         call h5dget_space_f(dset_temp, slabspace, hdf_error)
         call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                               data_offset, data_count, hdf_error)
         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                                 hdf_error)
          call h5dwrite_f(dset_temp, H5T_NATIVE_DOUBLE, &
            temp(1:n1,1:n2,kstart:kend), dims,  &
            hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
            xfer_prp = plist_id)
! ax
      call h5dget_space_f(dset_VOFx, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_VOFx, H5T_NATIVE_DOUBLE, &
       VOFx(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
! ay
      call h5dget_space_f(dset_VOFy, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_VOFy, H5T_NATIVE_DOUBLE, &
       VOFy(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
!az
      call h5dget_space_f(dset_VOFz, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_VOFz, H5T_NATIVE_DOUBLE, &
       VOFz(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
!vofp
         call h5dget_space_f(dset_VOFp, slabspace, hdf_error)
         call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                               data_offset, data_count, hdf_error)
         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                                 hdf_error)
          call h5dwrite_f(dset_VOFp, H5T_NATIVE_DOUBLE, &
          VOFp(1:n1,1:n2,kstart:kend), dims,  &
            hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
            xfer_prp = plist_id)

      call h5dclose_f(dset_vx, hdf_error)
      call h5dclose_f(dset_vy, hdf_error)
      call h5dclose_f(dset_vz, hdf_error)
      call h5dclose_f(dset_pr, hdf_error)
      call h5dclose_f(dset_temp, hdf_error)
      call h5dclose_f(dset_VOFx, hdf_error)
      call h5dclose_f(dset_VOFy, hdf_error)
      call h5dclose_f(dset_VOFz, hdf_error)
      call h5dclose_f(dset_VOFp, hdf_error)


      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      if (myid.eq.0) then

      open(45,file=xdmnam,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""thetacut"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""ORIGIN_DXDYDZ"">")')
      write(45,'("<DataItem Name=""Origin"" Dimensions=""3"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 0.0,0.0,0.0
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Name=""Spacing"" Dimensions=""3"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
      write(45,'(2E15.7)') 1./dx3,1./dx2,1./dx1
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""X-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/Vx")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Y-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/Vy")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Z-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/Vz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Pressure"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/Pr")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Temperature"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/Temp")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""VOFx"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/VOFx")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""VOFy"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/VOFy")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""VOFz"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/VOFz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""VOFp"" AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
      write(45,'("field_",i7.7,".h5:/VOFp")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      !write(45,'("<Time Value=""",e12.5,"""/>")')time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)
      endif
      end subroutine mpi_write_field

!================================================

!================================================

      subroutine mpi_write_vel
            use param
            use mpih
            use mpi_param, only: kstart,kend
            use local_arrays, only: vy,vz,vx,pr,temp
            use hdf5
            implicit none
            integer ic,jc,kc
            integer hdf_error
      
            integer(HID_T) :: file_id
            integer(HID_T) :: filespace
            integer(HID_T) :: slabspace
            integer(HID_T) :: memspace
      
            integer(HID_T) :: dset_vx
            integer(HID_T) :: dset_vy
            integer(HID_T) :: dset_vz
            integer(HSIZE_T) :: dims(3)
      
            integer(HID_T) :: plist_id
            integer(HSIZE_T), dimension(3) :: data_count
            integer(HSSIZE_T), dimension(3) :: data_offset
      
            integer(HSIZE_T) :: dims_grid(1)
            integer(HID_T) :: dset_grid
            integer(HID_T) :: dspace_grid
      
            integer :: comm, info
            integer :: ndims
            real tprfi
            integer itime
            character(70) filnam2,filnam3,filnam4,xdmnam
            character(70) filnamgrid,filnam5,filnam6
            character(7) ipfi
            tprfi = 1/tframe
            itime=nint(time*tprfi)
            write(ipfi,82)itime
         82 format(i7.7)
      
      !RO   Sort out MPI definitions
      
            comm = MPI_COMM_WORLD
            info = MPI_INFO_NULL
      
      !RO   Form the name of the file
      
      !       filnam2 = 'continuation/field.h5'
             filnam2 = 'continuation/velfield_'//ipfi//'.h5'
             xdmnam  = 'continuation/velfield_'//ipfi//'.xmf'
      
      !RO   Set offsets and element counts
         
            ndims = 3
      
            dims(1)=n1
            dims(2)=n2
            dims(3)=n3m
       
      
            call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      
            data_count(1) = n1
            data_count(2) = n2
            data_count(3) = kend-kstart+1
      
            data_offset(1) = 0
            data_offset(2) = 0
            data_offset(3) = kstart-1
      
      
      !EP   vx
      
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      
            call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      
            call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
             hdf_error, access_prp=plist_id)
      
            call h5pclose_f(plist_id, hdf_error)
      
            call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
                            filespace, dset_vx, hdf_error)
            call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
                            filespace, dset_vy, hdf_error)
            call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
                            filespace, dset_vz, hdf_error)
            call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      ! vx
            call h5dget_space_f(dset_vx, slabspace, hdf_error)
            call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                                  data_offset, data_count, hdf_error)
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                                    hdf_error)
             call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
               vx(1:n1,1:n2,kstart:kend), dims,  &
               hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
               xfer_prp = plist_id)
      ! vy
            call h5dget_space_f(dset_vy, slabspace, hdf_error)
            call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                                  data_offset, data_count, hdf_error)
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                                    hdf_error)
             call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
               vy(1:n1,1:n2,kstart:kend), dims,  &
               hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
               xfer_prp = plist_id)
      !vz
            call h5dget_space_f(dset_vz, slabspace, hdf_error)
            call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                                  data_offset, data_count, hdf_error)
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                                    hdf_error)
             call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
               vz(1:n1,1:n2,kstart:kend), dims,  &
               hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
               xfer_prp = plist_id)
      
            call h5dclose_f(dset_vx, hdf_error)
            call h5dclose_f(dset_vy, hdf_error)
            call h5dclose_f(dset_vz, hdf_error)
      
            call h5sclose_f(filespace, hdf_error)
            call h5pclose_f(plist_id, hdf_error)
            call h5sclose_f(memspace, hdf_error)
            call h5fclose_f(file_id, hdf_error)
      
            if (myid.eq.0) then
      
            open(45,file=xdmnam,status='unknown')
            rewind(45)
            write(45,'("<?xml version=""1.0"" ?>")')
            write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
            write(45,'("<Xdmf Version=""2.0"">")')
            write(45,'("<Domain>")')
            write(45,'("<Grid Name=""thetacut"" GridType=""Uniform"">")')
            write(45,'("<Topology TopologyType=""3DRectMesh"" NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
            write(45,'("<Geometry GeometryType=""ORIGIN_DXDYDZ"">")')
            write(45,'("<DataItem Name=""Origin"" Dimensions=""3"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
            write(45,'(2E15.7)') 0.0,0.0,0.0
            write(45,'("</DataItem>")')
            write(45,'("<DataItem Name=""Spacing"" Dimensions=""3"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
            write(45,'(2E15.7)') 1./dx3,1./dx2,1./dx1
            write(45,'("</DataItem>")')
            write(45,'("</Geometry>")')
            write(45,'("<Attribute Name=""X-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
            write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
            write(45,'("velfield_",i7.7,".h5:/Vx")') itime
            write(45,'("</DataItem>")')
            write(45,'("</Attribute>")')
            write(45,'("<Attribute Name=""Y-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
            write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
            write(45,'("velfield_",i7.7,".h5:/Vy")') itime
            write(45,'("</DataItem>")')
            write(45,'("</Attribute>")')
            write(45,'("<Attribute Name=""Z-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
            write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
            write(45,'("velfield_",i7.7,".h5:/Vz")') itime
            write(45,'("</DataItem>")')
            write(45,'("</Attribute>")')
            !write(45,'("<Time Value=""",e12.5,"""/>")')time
            write(45,'("</Grid>")')
            write(45,'("</Domain>")')
            write(45,'("</Xdmf>")')
            close(45)
            endif
      end subroutine mpi_write_vel
      
      !================================================


      subroutine mpi_write_tempField
            use param
            use mpih
            use mpi_param, only: kstart,kend
            use local_arrays, only: temp
            use hdf5
            implicit none
            integer ic,jc,kc
            integer hdf_error
      
            integer(HID_T) :: file_id
            integer(HID_T) :: filespace
            integer(HID_T) :: slabspace
            integer(HID_T) :: memspace
      
            integer(HID_T) :: dset_temp
            integer(HSIZE_T) :: dims(3)
      
            integer(HID_T) :: plist_id
            integer(HSIZE_T), dimension(3) :: data_count
            integer(HSSIZE_T), dimension(3) :: data_offset
      
            integer(HSIZE_T) :: dims_grid(1)
            integer(HID_T) :: dset_grid
            integer(HID_T) :: dspace_grid
      
            integer :: comm, info
            integer :: ndims
            real tprfi
            integer itime
            character(70) filnam2,filnam3,filnam4,xdmnam
            character(70) filnamgrid,filnam5,filnam6
            character(7) ipfi
            tprfi = 1/tframe
            itime=nint(time*tprfi)
            write(ipfi,82)itime
         82 format(i7.7)
      
      !RO   Sort out MPI definitions
      
            comm = MPI_COMM_WORLD
            info = MPI_INFO_NULL
      
      !RO   Form the name of the file
      
      !       filnam2 = 'continuation/field.h5'
             filnam2 = 'continuation/tempfield_'//ipfi//'.h5'
             xdmnam  = 'continuation/tempfield_'//ipfi//'.xmf'
      
      !RO   Set offsets and element counts
         
            ndims = 3
      
            dims(1)=n1
            dims(2)=n2
            dims(3)=n3m
       
      
            call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      
            data_count(1) = n1
            data_count(2) = n2
            data_count(3) = kend-kstart+1
      
            data_offset(1) = 0
            data_offset(2) = 0
            data_offset(3) = kstart-1
            
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      
            call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      
            call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
             hdf_error, access_prp=plist_id)
      
            call h5pclose_f(plist_id, hdf_error)
      
            call h5dcreate_f(file_id, 'Temp', H5T_NATIVE_DOUBLE, &
            filespace, dset_temp, hdf_error)
 
            call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
 
      !temp
               call h5dget_space_f(dset_temp, slabspace, hdf_error)
               call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                                     data_offset, data_count, hdf_error)
               call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
               call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                                       hdf_error)
                call h5dwrite_f(dset_temp, H5T_NATIVE_DOUBLE, &
                  temp(1:n1,1:n2,kstart:kend), dims,  &
                  hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
                  xfer_prp = plist_id)
      
            call h5dclose_f(dset_temp, hdf_error)

            call h5sclose_f(filespace, hdf_error)
            call h5pclose_f(plist_id, hdf_error)
            call h5sclose_f(memspace, hdf_error)
            call h5fclose_f(file_id, hdf_error)
      
            if (myid.eq.0) then
      
            open(45,file=xdmnam,status='unknown')
            rewind(45)
            write(45,'("<?xml version=""1.0"" ?>")')
            write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
            write(45,'("<Xdmf Version=""2.0"">")')
            write(45,'("<Domain>")')
            write(45,'("<Grid Name=""thetacut"" GridType=""Uniform"">")')
            write(45,'("<Topology TopologyType=""3DRectMesh"" NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
            write(45,'("<Geometry GeometryType=""ORIGIN_DXDYDZ"">")')
            write(45,'("<DataItem Name=""Origin"" Dimensions=""3"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
            write(45,'(2E15.7)') 0.0,0.0,0.0
            write(45,'("</DataItem>")')
            write(45,'("<DataItem Name=""Spacing"" Dimensions=""3"" NumberType=""Float"" Precision=""4"" Format=""XML"">")')
            write(45,'(2E15.7)') 1./dx3,1./dx2,1./dx1
            write(45,'("</DataItem>")')
            write(45,'("</Geometry>")')
            write(45,'("<Attribute Name=""Temperature"" AttributeType=""Scalar"" Center=""Node"">")')
            write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2,n1
            write(45,'("tempfield_",i7.7,".h5:/Temp")') itime
            write(45,'("</DataItem>")')
            write(45,'("</Attribute>")')
            !write(45,'("<Time Value=""",e12.5,"""/>")')time
            write(45,'("</Grid>")')
            write(45,'("</Domain>")')
            write(45,'("</Xdmf>")')
            close(45)
            endif
            end subroutine mpi_write_tempField


! KZ auxilary routines for timing code
subroutine tic(tstart)
      use mpih
      implicit none
      real :: tstart
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      tstart = MPI_WTIME()
end subroutine tic

subroutine toc(tstart,tend,elapsed)
      use mpih
      implicit none
      real :: tstart, tend, elapsed
          call MPI_Barrier(MPI_COMM_WORLD, ierr)
          tend = MPI_WTIME()
          elapsed = elapsed + (tend - tstart)
end subroutine toc

!======================================================

!       subroutine mpi_write_field_noParts
!       use param
!       use mpih
!       use mpi_param, only: kstart,kend
!       use local_arrays, only: vy,vz,vx,pr
!       use hdf5
!       implicit none
!       integer ic,jc,kc
!       integer hdf_error

!       integer(HID_T) :: file_id
!       integer(HID_T) :: filespace
!       integer(HID_T) :: slabspace
!       integer(HID_T) :: memspace

!       integer(HID_T) :: dset_vx
!       integer(HID_T) :: dset_vy
!       integer(HID_T) :: dset_vz
!       integer(HSIZE_T) :: dims(3)

!       integer(HID_T) :: plist_id
!       integer(HSIZE_T), dimension(3) :: data_count
!       integer(HSSIZE_T), dimension(3) :: data_offset

!       integer(HSIZE_T) :: dims_grid(1)
!       integer(HID_T) :: dset_grid
!       integer(HID_T) :: dspace_grid

!       integer :: comm, info
!       integer :: ndims
!       real tprfi
!       integer itime
!       character(70) filnam2,filnam3,filnam4
!       character(70) filnamgrid,filnam5,filnam6
!       character(7) ipfi
!       tprfi = 1/tframe
!       itime=nint(time*tprfi)
!       write(ipfi,82)itime
!    82 format(i7.7)

! !RO   Sort out MPI definitions

!       comm = MPI_COMM_WORLD
!       info = MPI_INFO_NULL

! !RO   Form the name of the file

!       filnam2 = 'continuation/velfield.h5'

!       !filnam2 = 'continuation/velfield_'//ipfi//'.h5'
!       !xdmnam  = 'continuation/velfield_'//ipfi//'.xmf'

! !RO   Set offsets and element counts
   
!       ndims = 3

!       dims(1)=n1
!       dims(2)=n2
!       dims(3)=n3m
 

!       call h5screate_simple_f(ndims, dims, filespace, hdf_error)

!       data_count(1) = n1
!       data_count(2) = n2
!       data_count(3) = kend-kstart+1

!       data_offset(1) = 0
!       data_offset(2) = 0
!       data_offset(3) = kstart-1


! !EP   vx

!       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

!       call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

!       call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
!        hdf_error, access_prp=plist_id)

!       call h5pclose_f(plist_id, hdf_error)

!       call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vx, hdf_error)
!       call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vy, hdf_error)
!       call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
!                       filespace, dset_vz, hdf_error)

!       call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
! ! vx
!       call h5dget_space_f(dset_vx, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
!          vx(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
! ! vy
!       call h5dget_space_f(dset_vy, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
!          vy(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)
! !vz
!       call h5dget_space_f(dset_vz, slabspace, hdf_error)
!       call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
!                             data_offset, data_count, hdf_error)
!       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
!                               hdf_error)
!        call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
!          vz(1:n1,1:n2,kstart:kend), dims,  &
!          hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
!          xfer_prp = plist_id)

!       call h5dclose_f(dset_vx, hdf_error)
!       call h5dclose_f(dset_vy, hdf_error)
!       call h5dclose_f(dset_vz, hdf_error)

!       call h5sclose_f(filespace, hdf_error)
!       call h5pclose_f(plist_id, hdf_error)
!       call h5sclose_f(memspace, hdf_error)
!       call h5fclose_f(file_id, hdf_error)
!       end subroutine mpi_write_field_noParts



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write particles for paraview

 !     subroutine write_parts


  !    do inp=1,Nparticle

!      call write_geom(maxnv,maxnf,xyzv(:,:,inp),vel_tri(:,:,inp))
!      enddo

!      end subroutine write_parts



! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! subroutine write_tecplot_geom
! use mpih
! use mpi_param
! use param
! use local_arrays, only: vy,vz,pr,vx
! use local_aux, only: vorx,vory,vorz
! use mls_param
! use coll_mod

! character(70) namfile,namfi2
! character(7) ipfi
! character*2 ibod
! real tprfi
! integer itime
! tprfi = 1/tframe
! itime=nint(time*tprfi)
! write(ipfi,82)itime
! 82 format(i7.7)
! 98 format(i2.2)

! if(ismaster)then
!     do inp=1,Nparticle

!    write(ibod,98) inp
!    write(ipfi,82) itime

! namfi2='continuation/G_'//ibod//'_'//ipfi

! call write_geom (maxnv,maxnf,xyzv(:,:,inp),tri_bar(:,:,inp),namfi2)
! end do
! end if

! end subroutine write_tecplot_geom

! !---------------------------------------------------------------------------------------------
!       subroutine write_geom (nv,nf,xyz,tri_vel,filename)
!       use param
!       use mpih
!       use mls_param
!       use mpi_param, only: kstart,kend
!       use local_arrays, only: vy,vz,vx,pr
!       implicit none
!       character(70) filename,geotecfile
!       integer ic,jc,kc,i,nv,nf
!       integer :: v1,v2,v3
!       real, dimension (3,nv) :: xyz
!       real, dimension (3,nf) :: tri_vel,node
!       real tprfi
!       integer itime
!       character(7) ipfi
!       tprfi = 1/tframe
!       itime=nint(time*tprfi)
!       write(ipfi,82)itime
!    82 format(i7.7)

     

!         geotecfile=trim(filename)//'.dat'
! !        write(*,*)' Write file ',trim(geotecfile)

!         open(11,file=geotecfile)

!         write(11,*)'TITLE = "Geo"'
!         write(11,*)'VARIABLES = X Y Z Vx Vy Vz'
!     !    write(11,*)'ZONE T="DOMAIN 0", N=',nv,' E=',nf,' F=FEBLOCK, ET=TRIANGLE'
!         write(11,*)'ZONE T="FETri" N=',nv,' E=',nf,' ZONETYPE=FETriangle'
!         ! write(11,*)'ZONE T=FETri N=',nvc,' E=',ntri,' ZONETYPE=FETriangle'
!         write(11,*)'DATAPACKING=BLOCK                                       '
!         write(11,*)'VARLOCATION=([4-6]=CELLCENTERED)'

!         do i=1,nv
!         write(11,*)xyz(1,i)
!         end do

!         do i=1,nv
!         write(11,*)xyz(2,i)
!         end do

!         do i=1,nv
!         write(11,*)xyz(3,i)
!         end do

!         do i=1,nf
!         write(11,*)tri_vel(1,i)
!         end do

!         do i=1,nf
!         write(11,*)tri_vel(2,i)
!         end do

!         do i=1,nf
!         write(11,*)tri_vel(3,i)
!         end do

!         do i=1,nf
!         write(11,*)vert_of_face(1:3,i)
!         end do


!         close(11)
!         return
!         end subroutine write_geom
! !---------------------------------------------------------------------------------------------

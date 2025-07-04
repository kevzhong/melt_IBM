!***********************************************************************
      subroutine InitArrays
      use param
      use local_arrays
      use local_aux, only: vorx, vory, vorz, diss, tke,chi
      use mpi_param, only: kstart,kend, buf_n1n2
      use mpih, only: lvlhalo
      use stat_arrays
      use mls_local
      use AuxiliaryRoutines
      implicit none
      integer :: j,k,kc,i

      !Solution arrays, u,v,w, temperature, pressure
      call AllocateReal3DArray(vx,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vy,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vz,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(temp,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(pr,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)

      ! Phase-indicator, internal-solid treatment
      allocate(VOFx(n1,n2,kstart-1:kend+1))
      allocate(VOFy(n1,n2,kstart-1:kend+1))
      allocate(VOFz(n1,n2,kstart-1:kend+1))
      allocate(VOFp(n1,n2,kstart-1:kend+1))
      !call AllocateReal3DArray(d_UsolidT_dxj,1,n1,1,n2,kstart,kend)
      !call AllocateLogical3DArray(solid_mask,1,n1,1,n2,kstart,kend)

      ! Default values for single-phase
      VOFx(:,:,:) = 1.
      VOFy(:,:,:) = 1.
      VOFz(:,:,:) = 1.
      VOFp(:,:,:) = 1.
      !solid_mask(:,:,:) = .false.

     allocate(sdf(n1,n2,kstart:kend))

      ! Auxilary fractional-step pseudo-pressure
      call AllocateReal3DArray(dph,1,n1,1,n2+1, &
          kstart-lvlhalo,kend+lvlhalo)

     ! Vorticity
      call AllocateReal3DArray(vorx,1,n1,1,n2,  &
           kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vory,1,n1,1,n2, &
           kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vorz,1,n1,1,n2,  &
           kstart-lvlhalo,kend+lvlhalo)

     ! dissipation, tke
     call AllocateReal3DArray(diss,1,n1,1,n2,kstart,kend)
     call AllocateReal3DArray(tke,1,n1,1,n2,kstart,kend)
     call AllocateReal3DArray(chi,1,n1,1,n2,kstart,kend)

     !HIT forcing
      call AllocateReal3DArray(forcx,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(forcy,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(forcz,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(rhs,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(dq,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(qcap,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(htemp,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru1,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru2,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru3,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(rut,1,n1,1,n2,kstart,kend)


      call AllocateReal3DArray(vx_me,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vy_me,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vz_me,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(pr_me,1,n1,1,n2,kstart,kend)

      call AllocateReal3DArray(vx_rms,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vy_rms,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vz_rms,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(pr_rms,1,n1,1,n2,kstart,kend)

      call AllocateReal3DArray(for_xc,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(for_yc,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(for_zc,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(for_temp,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)

      ! KZ: for add ghost routines: store memory in heap rather than stack for performance and large problems
      call AllocateReal2DArray(buf_n1n2,1,n1,1,n2)

      end   

subroutine apply_StefanCondition
USE param
USE mls_param
USE mpi_param, only: kstart, kend
use mpih
implicit none
real,dimension(4) :: ptx
real,dimension(3) :: pos
integer, dimension(3) :: probe_inds
real, dimension(3) :: gradT
integer :: inp,nv
real :: s
rhs_stefan(:,:,:) = 0.0d0

call MPI_ALLREDUCE(MPI_IN_PLACE,dtdn_o,maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
call MPI_ALLREDUCE(MPI_IN_PLACE,dtdn_i,maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        

! KZ: for now, therm. diffusivity kappa = 1/prandtl
do inp=1,Nparticle
    do nv = 1,maxnv

        ! ! Accumulate outward-normal term, outward = liquid
        ! if(pindv(3,nv,inp).ge.kstart-1 .and. pindv(3,nv,inp).le.kend+1) then
        !     rhs_stefan(:,nv,inp) = rhs_stefan(:,nv,inp) + al*dt * cpsolid / latHeat * ( -(1.0d0/prandtl)*dtdn_o(nv,inp) ) * vert_nor(:,nv,inp)
        ! endif

        ! ! Accumulate inward-normal term, inward = solid
        ! if(pindv(6,nv,inp).ge.kstart-1 .and. pindv(6,nv,inp).le.kend+1) then
        !     rhs_stefan(:,nv,inp) = rhs_stefan(:,nv,inp) + al*dt * cpsolid / latHeat * ( (1.0d0/prandtl)*dtdn_i(nv,inp) ) * vert_nor(:,nv,inp)
        ! endif

        ! Accumulate outward-normal term, outward = liquid
        vmelt_n(:,nv,inp) = cpsolid / latHeat * ( (1.0d0/prandtl)*dtdn_i(nv,inp) -(1.0d0/prandtl)*dtdn_o(nv,inp) ) * vert_nor(:,nv,inp)

    enddo
enddo
! Accumulate rhs across all processes
!call MPI_ALLREDUCE(MPI_IN_PLACE,rhs_stefan,3*maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
! Update vertice locations
! if (ismaster) then
!     do inp = 1,maxnv
!         write(*,*) "RHS stefan_x" , inp, "is ", rhs_stefan(1,inp,1)
!         !write(*,*) "dtdn_o " , inp, "is ", dtdn_o(inp,1)

!     enddo
! endif
xyzv(:,:,:) = xyzv(:,:,:) + al*dt*vmelt_n
!  if (ismaster) then
!      do inp = 1,maxnv
!          write(*,*) "RHS xv" , inp, "is ", xyzv(1,inp,1)
!          !write(*,*) "dtdn_o " , inp, "is ", dtdn_o(inp,1)

!      enddo
! endif
end subroutine apply_StefanCondition
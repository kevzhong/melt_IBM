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

! KZ: for now, therm. diffusivity kappa = 1/prandtl
do inp=1,Nparticle
    do nv = 1,maxnv

        ! Accumulate outward-normal term, outward = liquid
        if(pindv(3,nv,inp).ge.kstart-1 .and. pindv(3,nv,inp).le.kend+1) then
            rhs_stefan(:,nv,inp) = rhs_stefan(:,nv,inp) + dt * cpsolid / latHeat * ( -(1.0d0/prandtl)*dtdn_o(nv,inp) ) * vert_nor(:,nv,inp)
        endif

        ! Accumulate inward-normal term, inward = solid
        if(pindv(6,nv,inp).ge.kstart-1 .and. pindv(6,nv,inp).le.kend+1) then
            rhs_stefan(:,nv,inp) = rhs_stefan(:,nv,inp) + dt * cpsolid / latHeat * ( (1.0d0/prandtl)*dtdn_i(nv,inp) ) * vert_nor(:,nv,inp)
        endif

    enddo
enddo
! Accumulate rhs across all processes
call MPI_ALLREDUCE(MPI_IN_PLACE,rhs_stefan,3*maxnv*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        

! Update vertice locations
xyzv(:,:,:) = xyzv(:,:,:) + rhs_stefan
end subroutine apply_StefanCondition
subroutine apply_StefanCondition
USE param
USE mls_param
use mpih
implicit none
real, dimension(3) :: rhs_stefan
integer :: inp,nv

! KZ: for now, therm. diffusivity kappa = 1/prandtl
do inp=1,Nparticle
    do nv = 1,maxnv

        ! Accumulate outward-normal term, outward = liquid
        rhs_stefan(1:3) =  cpsolid / latHeat * ( (1.0d0/prandtl)*dtdn_iVert(nv,inp) -(1.0d0/prandtl)*dtdn_oVert(nv,inp) ) * vert_nor(1:3,nv,inp)
        vmelt(nv,inp) = dot_product( rhs_stefan,vert_nor(1:3,nv,inp) ) ! Store the scalar local melting velocity

        ! Update interface location
        xyzv(1:3,nv,inp) = xyzv(1:3,nv,inp) +  al * dt * rhs_stefan(1:3)
    enddo
enddo

end subroutine apply_StefanCondition
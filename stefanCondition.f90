subroutine apply_StefanCondition
USE param
USE mls_param
use mpih
implicit none
real, dimension(3) :: rhs_stefan, nhat
integer :: inp,nv

! KZ: for now, therm. diffusivity kappa = 1/prandtl
do inp=1,Nparticle
    do nv = 1,maxnv
        if (isGhostVert(nv,inp) .eqv. .false. ) then
            nhat(1:3) = vert_nor(1:3,nv,inp)
            ! Accumulate outward-normal term, outward = liquid
            vmelt(nv,inp) = cpliquid / latHeat * ( (1.0d0/prandtl)*dtdn_iVert(nv,inp) -(1.0d0/prandtl)*dtdn_oVert(nv,inp) )
            rhs_stefan(1:3) = vmelt(nv,inp)*nhat(1:3)
            !rhs_stefan(1:3) =  cpliquid / latHeat * ( (1.0d0/prandtl)*dtdn_iVert(nv,inp) -(1.0d0/prandtl)*dtdn_oVert(nv,inp) ) * nhat
            !vmelt(nv,inp) = dot_product( rhs_stefan, nhat ) ! Store the scalar local melting velocity
            
            ! Update interface location
            xyzv(1:3,nv,inp) = xyzv(1:3,nv,inp) +  al * dt * rhs_stefan(1:3)
        endif
    enddo
enddo

end subroutine apply_StefanCondition
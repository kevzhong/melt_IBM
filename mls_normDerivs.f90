subroutine mls_normDerivs
USE param
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real,dimension(4) :: ptx
real,dimension(3) :: pos
integer, dimension(3) :: probe_inds
real, dimension(3) :: gradT
integer :: inp,nf
real :: s

! Evaluate dTdn at probe locations extrapolated from triangle centroids

dtdn_o(:,:) = 0.0d0
dtdn_i(:,:) = 0.0d0

dtdn_oVert(:,:) = 0.0d0
dtdn_iVert(:,:) = 0.0d0


gradT = 0.0d0

do inp=1,Nparticle
 do nf = 1, maxnf

      if(pind(3,nf,inp).ge.kstart .and. pind(3,nf,inp).le.kend) then
        s = norm2( tri_nor(1:3,nf,inp) )
        !------------------ POSITIVE PROBE-------------------------
        pos(1:3) = tri_bar(1:3,nf,inp) + h_eulerian * tri_nor(1:3,nf,inp) / s !Positive probe location

         ! initialise pre-factor matrix
         ptx(1)   = 1.0d0 
         ptx(2:4) = pos(1:3)

         probe_inds(1:3) = pind(1:3,nf,inp) ! Positive probe cage-indices

         call wght_gradT(pos,ptx,gradT,probe_inds)
         dtdn_o(nf,inp) = tri_nor(1,nf,inp) * gradT(1) & ! nx dTdx
                        + tri_nor(2,nf,inp) * gradT(2) & ! ny dTdy
                        + tri_nor(3,nf,inp) * gradT(3)   ! nz dTdz

         ! Transfer the face A*dTdn to its vertices, to dtdn_oVert and area, Avert 
         call faceToVert_interp(vert_of_face(1:3,nf),sur(nf,inp),Avert(:,inp),dtdn_o(nf,inp),dtdn_oVert(:,inp),maxnv)

         !write(*,*) "centroid loc", nf, " is ", tri_bar(1:3,nf,1)
      endif !end if pind(3...)

        !------------------ NEGATIVE PROBE-------------------------
      if(pind(6,nf,inp).ge.kstart .and. pind(6,nf,inp).le.kend) then
        s = norm2( tri_nor(1:3,nf,inp) )
        pos(1:3) = tri_bar(1:3,nf,inp) - h_eulerian * tri_nor(1:3,nf,inp) / s !Negative probe location

         ! initialise pre-factor matrix
         ptx(1)   = 1.0d0 
         ptx(2:4) = pos(1:3)

         probe_inds(1:3) = pind(4:6,nf,inp) ! Negative probe cage-indices

         call wght_gradT(pos,ptx,gradT,probe_inds)
         dtdn_i(nf,inp) = tri_nor(1,nf,inp) * gradT(1) & ! nx dTdx
                        + tri_nor(2,nf,inp) * gradT(2) & ! ny dTdy
                        + tri_nor(3,nf,inp) * gradT(3)   ! nz dTdz

        ! Note the same normal-sign convetion
          
        call faceToVert_interp(vert_of_face(1:3,nf),sur(nf,inp),Avert(:,inp),dtdn_i(nf,inp),dtdn_iVert(:,inp),maxnv)

      endif !end if pindv(6...)

 enddo
enddo
end subroutine mls_normDerivs


subroutine wght_gradT(pos,ptx,gradT,probe_inds)
USE param
USE geom
use local_arrays, only: temp
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real, dimension(3), intent(inout) :: gradT !dTdx, dTdy, dTdz at the probe
real,dimension(nel) :: ddx_PtxAB(nel), ddy_PtxAB(nel), ddz_PtxAB(nel) !Shape function derivatives
real,dimension(nel) :: Tnel ! Eulerian support domain values of temperature
real,dimension(4) :: pxk ! [1,x,y,z] Basis function evaluated at Eulerian cage coordinates
real,dimension(4), intent(in) :: ptx ! Basis function [1, x, y, z] where (x,y,z) evaluated at the probe location
real,dimension(3) :: pos,norp,Wt ! Probe location, normalised cage distances and weights
real,dimension(4,4) :: pinvA,invA, dAdx, dAdy, dAdz
real,dimension(4,nel) :: B, dBdx, dBdy, dBdz
real, dimension(4) :: gamvec(4), dgamdx(4), dgamdy(4), dgamdz(4) ! Auxilary vectors for computing shape function derivatives
real :: Wtx, Wt23
real :: dWx_dx, dWy_dy, dWz_dz
real :: dWdx, dWdy, dWdz
integer :: inp,inw,i,j,k,k1
integer, dimension(3) :: pind_i, pind_o, probe_inds

!-------------Shape function for cell centres (temp. or pressure cells) -------------------------

if(probe_inds(3).ge.kstart .and. probe_inds(3).le.kend) then
  
!-------------FORCING FUNCTION------------------------
! volume of a face with a specific marker - thickness taken as average of grid spacing
  
!WGHT1
!pind_i(1)=pindv(1,nv,inp)-1;  pind_o(1)=pindv(1,nv,inp)+1
!pind_i(2)=pindv(2,nv,inp)-1;  pind_o(2)=pindv(2,nv,inp)+1
!pind_i(3)=pind(3,ntr,inp)-1;  pind_o(3)=pind(3,ntr,inp)+1

pind_i(1)=probe_inds(1)-1;  pind_o(1)=probe_inds(1)+1
pind_i(2)=probe_inds(2)-1;  pind_o(2)=probe_inds(2)+1
pind_i(3)=probe_inds(3)-1;  pind_o(3)=probe_inds(3)+1


  
k1  = floor(pos(3)*dx3) + 1
pind_i(3)=k1-1
pind_o(3)=k1+1
  
pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
dAdx(1:4,1:4)=0.0d0
dAdy(1:4,1:4)=0.0d0
dAdz(1:4,1:4)=0.0d0

inw = 1
  
  
! Accumulate A(4,4)   , B(4,27) linear system
! Likewise for derivatives of A, B
do k=pind_i(3), pind_o(3)
  
    norp(3)=abs( pos(3) - zm(k) ) / (wscl / dx3)
    Wt(3) = mls_gaussian( norp(3) , wcon )
    dWz_dz = mls_gauss_deriv(norp(3),wcon) ! dWz / dr_z
    dWz_dz = dWz_dz * ( ( pos(3) - zm(k) ) / abs( pos(3) - zm(k)  ) /  (wscl/dx3) ) ! dWz/dz = ( dWz / dr_z ) * ( dr_z / dz )
    do j=pind_i(2), pind_o(2)
  
        norp(2)=abs( pos(2) - ym(j) ) / (wscl / dx2)
        Wt(2) = mls_gaussian( norp(2) , wcon )
        Wt23 = Wt(2)*Wt(3)
  
        dWy_dy = mls_gauss_deriv(norp(2),wcon) ! dWy / dr_y
        dWy_dy = dWy_dy * ( ( pos(2) - ym(j) ) / abs( pos(2) - ym(j)  ) /  (wscl/dx2) ) ! dWy/dy = ( dWy / dr_y ) * ( dr_y / dy )
        do i=pind_i(1), pind_o(1)
  
            norp(1)=abs(pos(1) - xm(i) ) / (wscl / dx1)
            Wt(1) = mls_gaussian( norp(1) , wcon )
            Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)

            dWx_dx = mls_gauss_deriv(norp(1),wcon) ! dWx / dr_x
            dWx_dx = dWx_dx * ( ( pos(1) - xm(i) ) / abs( pos(1) - xm(i)  ) /  (wscl/dx1) ) ! dWx/dx = ( dWx / dr_x ) * ( dr_x / dx )

            pxk(1)=1.0d0
            pxk(2)=xm(i)
            pxk(3)=ym(j)
            pxk(4)=zm(k)
            ! C = alpha * A * B + beta * C
            ! Note how the result of the computation is stored in C
            !DGEMM(transA,transB,numRowsA,numColsB,numColsA,alpha,  A   ,numRowsA,  B     , numRowsB, beta , C/result, numRowsC)
            call DGEMM('N','T',4,4,1,Wtx,pxk,4,pxk,4, 1.0d0,pinvA,4) ! Accumulate A matrix ( eventually inverted to pinvA )
            B(1:4,inw)=Wtx*pxk(1:4)

            ! Eq. (3.165) and Program 3.5, 3.6 in Liu & Gu 2005
            dWdx = dWx_dx * Wt(2)  * Wt(3)
            dWdy = Wt(1)  * dWy_dy * Wt(3)
            dWdz = Wt(1)  * Wt(2)  * dWz_dz
  
            call DGEMM('N','T',4,4,1,dWdx,pxk,4,pxk,4, 1.0d0,dAdx,4) 
            call DGEMM('N','T',4,4,1,dWdy,pxk,4,pxk,4, 1.0d0,dAdy,4) 
            call DGEMM('N','T',4,4,1,dWdz,pxk,4,pxk,4, 1.0d0,dAdz,4) 

            dBdx(1:4,inw)=dWdx * pxk(1:4)
            dBdy(1:4,inw)=dWdy * pxk(1:4)
            dBdz(1:4,inw)=dWdz * pxk(1:4)

            Tnel(inw) = temp(i,j,k) ! Store Eulerian support domain values for summation later

            inw = inw + 1
        enddo !end i
    enddo !end j
enddo !end k
  
    ! calling routine to compute inverse
    ! SPD matrix for uniform grids, we can use Cholesky decomp. instead: dpotrf
    call inverseLU(pinvA,invA)

    !---------------------LIN-ALG ROUTINES----------------------------
    ! matrix multiplications for final interpolation
    !DGEMM(transA,transB,numRowsA,numColsB,numColsA,alpha,  A   ,numRowsA,  B     , numRowsB, beta , C/result, numRowsC)
    ! C = alpha * A * B + beta * C
    ! Note how the result of the computation is stored in C
    ! 
    !
    ! Matrix-vector operations: y = alpha * Ax + beta * y   ; result stored in y   incx incy are strides (should be 1)
    !DGEMV(TPOSE?, nRowA, nColA, alpha, A, nRowA, x, incx, beta, y, incy) 
    !
    ! NB: for the intel MKL implementations, must have B =/= C ! But this is not the case for LAPACK

    !---------------SHAPE FUNCTION DERIVATIVE CALCULATIONS---------------
    ! For evaluation of shape function derivatives, ddx_ptxAB, ddy_ptxAB, ddz_ptxAB
    ! See equations (3.138--3.141 , 3.144) in Liu & Gu (2005)
    ! Our MLS shape function ptxAB := Phi in their notation
    ! Notably, only inversion of matrix A is required in this framework (stored in invA)
    ! Also no need for evaluating the shape function (i.e. ptxAB), as we only need its derivatives here

    ! Eq. (3.139 / 3.140) Compute gamma = invA * p  , result gamvec = size(4) vector
    call DGEMV('N',    4, 4, 1.0d0, invA, 4, ptx, 1, 0.0d0, gamvec, 1 ) 

    ! Now compute derivatives of gamma for use in (3.141)

    !---------------------- x derivative ddx_ptxAB-------------------------------------------
    pxk = [0.0d0,1.0d0,0.0d0,0.0d0] ! Store d/dx(ptx) initially, re-use pxk memory
    call DGEMV('N',    4, 4, -1.0d0, dAdx, 4, gamvec, 1, 1.0d0, pxk, 1 ) ! Get -dAdx*gamma + dpdx first, store in pxk
    !Now get dgamdx = Ainv * (-dAdx*gamma + dpdx) where (-dAdx*gamma + dpdx) was computed above and is stored in pxk
    call DGEMV('N',    4, 4, 1.0d0, invA, 4, pxk, 1, 0.0d0, dgamdx, 1 ) 

    ! Compute MLS shape function derivative ddx_ptxAB = dgamdx^T*B + gamma^T*dBdx
    call DGEMM('N', 'N', 1, nel, 4, 1.0d0, dgamdx, 1, B, 4, 0.0d0, ddx_ptxAB, 1) ! Term dgamdx^T*B
    call DGEMM('N', 'N', 1, nel, 4, 1.0d0, gamvec, 1, dBdx, 4 , 1.0d0 ,ddx_ptxAB, 1) ! += gamma^T*dBdx



    !---------------------- y derivative ddy_ptxAB-------------------------------------------
    pxk = [0.0d0,0.0d0,1.0d0,0.0d0] ! Store d/dy(ptx) initially
    call DGEMV('N',    4, 4, -1.0d0, dAdy, 4, gamvec, 1, 1.0d0, pxk, 1 )
    call DGEMV('N',    4, 4, 1.0d0, invA, 4, pxk, 1, 0.0d0, dgamdy, 1 ) 
    call DGEMM('N', 'N', 1, nel, 4, 1.0d0, dgamdy, 1, B, 4, 0.0d0, ddy_ptxAB, 1)
    call DGEMM('N', 'N', 1, nel, 4, 1.0d0, gamvec, 1, dBdy, 4 , 1.0d0 ,ddy_ptxAB, 1)

    !---------------------- z derivative ddz_ptxAB-------------------------------------------
    pxk = [0.0d0,0.0d0,0.0d0,1.0d0] ! Store d/dz(ptx) initially
    call DGEMV('N',    4, 4, -1.0d0, dAdz, 4, gamvec, 1, 1.0d0, pxk, 1 )
    call DGEMV('N',    4, 4, 1.0d0, invA, 4, pxk, 1, 0.0d0, dgamdz, 1 ) 
    call DGEMM('N', 'N', 1, nel, 4, 1.0d0, dgamdz, 1, B, 4, 0.0d0, ddz_ptxAB, 1)
    call DGEMM('N', 'N', 1, nel, 4, 1.0d0, gamvec, 1, dBdz, 4 , 1.0d0 ,ddz_ptxAB, 1)
    
    !-------------------- Now compute derivatives at the marker location---------------------
    ! e.g. Vanella & Balaras (2009) equation (27)

    gradT(1) = sum( ddx_PtxAB * Tnel) !dT dx
    gradT(2) = sum( ddy_PtxAB * Tnel) !dT dy
    gradT(3) = sum( ddz_PtxAB * Tnel) !dT dz

  endif
  end subroutine wght_gradT

  subroutine faceToVert_interp(vert_of_face,Atri,Avert,dtdn_F,dtdn_v,nv)

    ! For a SINGLE triangle/face, spreads the area-weighted heat flux and area to its 3 vertices
    implicit none

    integer :: i, nv
    integer, dimension(3) :: vert_of_face ! 3 vertices of the specified triangle
    real, intent(in) :: Atri, dtdn_F
    real, dimension(nv), intent(out)  :: dtdn_v ! Vertex heat flux
    real, dimension(nv), intent(in)  :: Avert ! Area associated with vertices

    do i = 1,3 !vertices v1, v2, v3
      dtdn_v( vert_of_face(i) ) = dtdn_v( vert_of_face(i) ) + dtdn_F * ( Atri / Avert( vert_of_face(i) ) )
    enddo

  end subroutine faceToVert_interp

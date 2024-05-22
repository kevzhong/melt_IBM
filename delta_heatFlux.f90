subroutine mls_heatFlux
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

! Evaluate heat flux at probe locations extrapolated from triangle centroids

qw_o(:,:) = 0.0d0
qw_i(:,:) = 0.0d0

qw_oVert(:,:) = 0.0d0
qw_iVert(:,:) = 0.0d0


gradT = 0.0d0

do inp=1,Nparticle
 do nf = 1, maxnf
  if (isGhostFace(nf,inp) .eqv. .false. ) then
      if(pind(3,nf,inp).ge.kstart .and. pind(3,nf,inp).le.kend) then
        s = norm2( tri_nor(1:3,nf,inp) )
        !------------------ POSITIVE PROBE-------------------------
        pos(1:3) = tri_bar(1:3,nf,inp) + h_eulerian * tri_nor(1:3,nf,inp) / s !Positive probe location

         ! initialise pre-factor matrix
         ptx(1)   = 1.0d0 
         ptx(2:4) = pos(1:3)

         probe_inds(1:3) = pind(1:3,nf,inp) ! Positive probe cage-indices

         call wght_gradT(pos,gradT,probe_inds)
         qw_o(nf,inp) = tri_nor(1,nf,inp) * gradT(1) & ! nx dTdx
                        + tri_nor(2,nf,inp) * gradT(2) & ! ny dTdy
                        + tri_nor(3,nf,inp) * gradT(3)   ! nz dTdz

         qw_o(nf,inp) = qw_o(nf,inp) / pec

         ! Transfer the face A*qw to its vertices, to qw_oVert and area, Avert 
         call faceToVert_interp(vert_of_face(1:3,nf,inp),sur(nf,inp),Avert(:,inp),qw_o(nf,inp),qw_oVert(:,inp),maxnv)

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

         call wght_gradT(pos,gradT,probe_inds)
         qw_i(nf,inp) = tri_nor(1,nf,inp) * gradT(1) & ! nx dTdx
                        + tri_nor(2,nf,inp) * gradT(2) & ! ny dTdy
                        + tri_nor(3,nf,inp) * gradT(3)   ! nz dTdz

        qw_i(nf,inp) = qw_i(nf,inp) / pec

        ! Note the same normal-sign convetion
          
        call faceToVert_interp(vert_of_face(1:3,nf,inp),sur(nf,inp),Avert(:,inp),qw_i(nf,inp),qw_iVert(:,inp),maxnv)

      endif !end if pindv(6...)
    endif
 enddo
enddo
end subroutine mls_heatFlux


subroutine wght_gradT(pos,gradT,probe_inds)
USE param
USE geom
use local_arrays, only: temp
USE mls_param
USE mpi_param, only: kstart, kend
implicit none
real, dimension(3), intent(inout) :: gradT !dTdx, dTdy, dTdz at the probe
real,dimension(3) :: pos,norp,Wt ! Probe location, normalised cage distances and weights
real :: Wtx, Wt23
real :: dWx_dx, dWy_dy, dWz_dz
real :: dWdx, dWdy, dWdz
integer :: inp,inw,i,j,k,k1
integer :: ii,jj,kk,nk

integer, dimension(3) :: pind_i, pind_o, probe_inds, keul_inds, k_inds

!-------------Shape function for cell centres (temp. or pressure cells) -------------------------

if(probe_inds(3).ge.kstart .and. probe_inds(3).le.kend) then
  
!-------------FORCING FUNCTION------------------------
! volume of a face with a specific marker - thickness taken as average of grid spacing

pind_i(1)=probe_inds(1)-1;  pind_o(1)=probe_inds(1)+1
pind_i(2)=probe_inds(2)-1;  pind_o(2)=probe_inds(2)+1
pind_i(3)=probe_inds(3)-1;  pind_o(3)=probe_inds(3)+1


  
k1  = floor(pos(3)*dx3) + 1
pind_i(3)=k1-1
pind_o(3)=k1+1

! For the spatial grid
k_inds = [k1-1,k1,k1]

! For the Eulerian field
k1  = modulo(k1-1,n3m)  + 1
keul_inds = [k1-1,k1,k1+1]

inw = 1
  
  
gradT(:) = 0.0

!do k=pind_i(3), pind_o(3)
 do nk = 1,3 ! Different to apply periodicty in z
    k = k_inds(nk)
    norp(3) = abs(zm(k)-pos(3))*dx3
    !norp(3) = (zm(k)-pos(3))*dx3
    !norp(3) = (pos(3) - zm(k) )*dx3

    Wt(3) = delta(norp(3))!*dx3
    dWz_dz = delta_deriv(norp(3)) ! dWz / dr_z
    !dWz_dz = dWz_dz * dx3 ! dWz/dz = ( dWz / dr_z ) * ( dr_z / dz )
    dWz_dz = dWz_dz * dx3 * ( ( zm(k) - pos(3) ) / abs(  zm(k) - pos(3) )  )
    do j=pind_i(2), pind_o(2)
  
        norp(2) = abs(ym(j)-pos(2))*dx2  
        !norp(2) = (ym(j)-pos(2))*dx2  
        !norp(2) = (pos(2) - ym(j) )*dx2

        Wt(2) = delta(norp(2))!*dx2
        dWy_dy = delta_deriv(norp(2)) ! dWy / dr_y
        !dWy_dy = dWy_dy * dx2 ! dWy/dy = ( dWy / dr_y ) * ( dr_y / dy )
        dWy_dy = dWy_dy * dx2 * ( (  ym(j) - pos(2) ) / abs(  ym(j) - pos(2) )  )

        do i=pind_i(1), pind_o(1)
  
            norp(1) = abs(xm(i)-pos(1))*dx1  
            !norp(1) = (xm(i)-pos(1))*dx1  
            !norp(1) = (pos(1) - xm(i) )*dx1
            Wt(1) = delta(norp(1))!*dx1
            dWx_dx = delta_deriv(norp(1)) ! dWx / dr_x
            !dWx_dx = dWx_dx * dx1 ! dWx/dx = ( dWx / dr_x ) * ( dr_x / dx )
            dWx_dx = dWx_dx * dx1 * ( (  xm(i) - pos(1) ) / abs(  xm(i) - pos(1) )  )

            ! Eq. (21) Roma et al. (1999) differentiated
            dWdx = dWx_dx * Wt(2)  * Wt(3)  !/ (dx1 * dx2 * dx3)
            dWdy = Wt(1)  * dWy_dy * Wt(3)  !/ (dx1 * dx2 * dx3)
            dWdz = Wt(1)  * Wt(2)  * dWz_dz !/ (dx1 * dx2 * dx3)

            !Tnel(inw) = temp(i,j,k) ! Store Eulerian support domain values for summation later

            ii = modulo(i-1,n1m) + 1
            jj = modulo(j-1,n2m) + 1
            !kk = modulo(k-1,n3m) + 1

            kk = keul_inds(nk)

            gradT(1) = gradT(1) + dWdx * temp(ii,jj,kk)
            gradT(2) = gradT(2) + dWdy * temp(ii,jj,kk)
            gradT(3) = gradT(3) + dWdz * temp(ii,jj,kk)

            inw = inw + 1
        enddo !end i
    enddo !end j
enddo !end k

  endif
  end subroutine wght_gradT

  subroutine faceToVert_interp(vert_of_face,Atri,Avert,qw_F,qw_v,nv)

    ! For a SINGLE triangle/face, spreads the area-weighted heat flux and area to its 3 vertices
    implicit none

    integer :: i, nv
    integer, dimension(3) :: vert_of_face ! 3 vertices of the specified triangle
    real, intent(in) :: Atri, qw_F
    real, dimension(nv), intent(out)  :: qw_v ! Vertex heat flux
    real, dimension(nv), intent(in)  :: Avert ! Area associated with vertices

    do i = 1,3 !vertices v1, v2, v3
      qw_v( vert_of_face(i) ) = qw_v( vert_of_face(i) ) + qw_F * ( Atri / Avert( vert_of_face(i) ) )
    enddo

  end subroutine faceToVert_interp

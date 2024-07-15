subroutine mls_structLoads
    USE param
    USE mls_param
    USE mpi_param, only: kstart, kend
    use mpih
    implicit none
    real,dimension(4) :: ptx
    real,dimension(3) :: pos, rad
    integer, dimension(3) :: probe_inds
    real, dimension(3) :: gradP, gradU, gradV, gradW
    integer :: inp,nf
    integer :: icent, istag, jcent, jstag, kcent, kstag
    real :: s, press_probe
    real :: S11, S12, S13
    real :: S21, S22, S23
    real :: S31, S32, S33
    real :: sum_buffer

    
    ! Evaluate pressures and viscous stresses at Lagrangian nodes. Surface values extrapolated from probe locations
    ! Pre-multiplications by triangle areas also done to obtain surface-integrated loads, so result is effectively a force on triangles
    press_tri(:,:) = 0.0d0
    tau_n1(:,:) = 0.0d0
    tau_n2(:,:) = 0.0d0
    tau_n3(:,:) = 0.0d0

    
    gradP = 0.0d0
    gradU = 0.0d0
    gradV = 0.0d0
    gradW = 0.0d0

    do inp=1,Nparticle
     do nf = 1, maxnf
      if (isGhostFace(nf,inp) .eqv. .false. ) then
        ! Reset
        S11 = 0.0 ; S12 = 0.0; S13 = 0.0
        S21 = 0.0 ; S22 = 0.0; S23 = 0.0
        S31 = 0.0 ; S32 = 0.0; S33 = 0.0
          !if(pind(3,nf,inp).ge.kstart .and. pind(3,nf,inp).le.kend) then
            s = norm2( tri_nor(1:3,nf,inp) )
            !------------------ POSITIVE PROBE-------------------------
            pos(1:3) = tri_bar(1:3,nf,inp) + h_eulerian * tri_nor(1:3,nf,inp) / s !Positive probe location

            rad(1:3) = tri_bar(1:3,nf,inp) - pos_cm(1:3,inp) ! Moment arm for torques
    
             kcent  = floor(pos(3)*dx3) + 1
             kstag = nint(pos(3)*dx3)  + 1

             kcent  = modulo(kcent-1,n3m)  + 1    
             kstag = modulo(kstag-1,n3m)  + 1 

            ! initialise pre-factor matrix
             ptx(1)   = 1.0d0 
             ptx(2:4) = pos(1:3)

             jcent  = floor(pos(2)*dx2) + 1
             jstag = nint(pos(2)*dx2)  + 1
             
             icent  = floor(pos(1)*dx1) + 1
             istag = nint(pos(1)*dx1)  + 1
    
            
            if(kcent.ge.kstart .and. kcent.le.kend) then
             !----------------- Pressure: cell-centered grid -------------------------

            ! Extrapolation from probe value x_p to L-node location X:
            !   _     _       _
            !   x   = X   + h n 
            !    p

            ! Applying Taylor series expansion for P(X):
            !
            ! P(X) = P(x_p) - h (dp/dn)_{x_p}
            
             probe_inds(1:3) = [icent, jcent, kcent]

             call wght_press(pos,ptx,press_probe,probe_inds) !Calc. pressure at probe

             call wght_gradP(pos,ptx,gradP,probe_inds) ! Calc dP/dxi at probe

             press_tri(nf,inp) = press_probe - h_eulerian * dot_product( gradP, tri_nor(1:3,nf,inp) )

             ! Convert pressure to force
             press_tri(nf,inp) = press_tri(nf,inp) * sur(nf,inp) 

            !----------------- Viscous stress calculations -------------------------

            ! Denote S_ij := dui / dxj as the velocity gradient tensor
             !Viscous stress tensor is given by
             ! tau_ij := nu * ( dui/dxj + duj / dxi ) := nu * (S_ij + S_ji)

             !----------------- Project velocity gradients to normal direction ----------------
             !                [   2 S_11    ]   [  nhat_x  ]
             ! tau_n1 := nu * |  S12 + S21  | . [  nhat_y  ]
             !                [  S13 + S31  ]   [  nhat_z  ]


             !                [  S21 + S12  ]   [  nhat_x  ]
             ! tau_n2 := nu * |   2 S_22    | . [  nhat_y  ]
             !                [  S23 + S32  ]   [  nhat_z  ]

             !                [  S31 + S13  ]   [  nhat_x  ]
             ! tau_n3 := nu * |  S32 + S23  | . [  nhat_y  ]
             !                [   2 S_33    ]   [  nhat_z  ]
             !
             ! Note the separate if condition for W-velocity gradients due to staggered grid
             ! Also note the premultiplication by triangle area, sur, to obtain a differential force

            !----------------- U velocity gradients ---------------------------------
             probe_inds(1:3) = [istag, jcent, kcent]
             call wght_gradU(pos,ptx,gradU,probe_inds) ! Calc dU/dxi at probe
             S11 = gradU(1) ; S12 = gradU(2) ; S13 = gradU(3)

            !----------------- V velocity gradients ---------------------------------
             probe_inds(1:3) = [icent, jstag, kcent]
             call wght_gradV(pos,ptx,gradV,probe_inds) ! Calc dV/dxi at probe
             S21 = gradV(1) ; S22 = gradV(2) ; S23 = gradV(3)

             !--------------- Accumulate (U,V) contributions to viscous stresses
             tau_n1(nf,inp) = tau_n1(nf,inp) + (sur(nf,inp)/ren) * ( 2.0*S11    * tri_nor(1,nf,inp) + & 
                                                                    (S12 + S21) * tri_nor(2,nf,inp) + & 
                                                                    S13         * tri_nor(3,nf,inp) )

             tau_n2(nf,inp) = tau_n2(nf,inp) + (sur(nf,inp)/ren) * ((S21 + S12)    * tri_nor(1,nf,inp) + & 
                                                                      2.0*S22      * tri_nor(2,nf,inp) + & 
                                                                        S23        * tri_nor(3,nf,inp) )    
                                                                
            tau_n3(nf,inp) = tau_n3(nf,inp) + (sur(nf,inp)/ren) * (    S13      * tri_nor(1,nf,inp) + & 
                                                                       S23      * tri_nor(2,nf,inp) )
            endif

            if(kstag.ge.kstart .and. kstag.le.kend) then
                !----------------- W velocity gradients ---------------------------------
                 probe_inds(1:3) = [icent, jcent, kstag]
                 call wght_gradW(pos,ptx,gradW,probe_inds) ! Calc dW/dxi at probe
                 S31 = gradW(1) ; S32 = gradW(2) ; S33 = gradW(3)

                !--------------- Accumulate W contributions to viscous stresses
                 tau_n1(nf,inp) = tau_n1(nf,inp) + (sur(nf,inp)/ren) * (S31 * tri_nor(3,nf,inp) )

                 tau_n2(nf,inp) = tau_n2(nf,inp) + (sur(nf,inp)/ren) * (S32 * tri_nor(3,nf,inp) )

                 tau_n3(nf,inp) = tau_n3(nf,inp) + (sur(nf,inp)/ren) * ( S31    * tri_nor(1,nf,inp) + & 
                                                                         S32    * tri_nor(2,nf,inp) + &
                                                                        2.0*S33 * tri_nor(3,nf,inp) )

            endif


        
        endif ! end if GhostFace
     enddo
    enddo

    ! Sum across all processors (triangles)
    do inp=1,Nparticle

        int_pr_dA(inp) = sum( press_tri(:,inp) )
        call MPI_ALLREDUCE(MPI_IN_PLACE,int_pr_dA(inp),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)     
        
        int_tau_dA(1,inp) = sum( tau_n1(:,inp) )
        call MPI_ALLREDUCE(MPI_IN_PLACE,int_tau_dA(1,inp),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)     

        int_tau_dA(2,inp) = sum( tau_n2(:,inp) )
        call MPI_ALLREDUCE(MPI_IN_PLACE,int_tau_dA(2,inp),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  

        int_tau_dA(3,inp) = sum( tau_n3(:,inp) )
        call MPI_ALLREDUCE(MPI_IN_PLACE,int_tau_dA(3,inp),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
    enddo

    !call MPI_ALLREDUCE(MPI_IN_PLACE,press_tri,maxnf*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
    !call MPI_ALLREDUCE(MPI_IN_PLACE,tau_n1,maxnf*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
    !call MPI_ALLREDUCE(MPI_IN_PLACE,tau_n2,maxnf*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        
    !call MPI_ALLREDUCE(MPI_IN_PLACE,tau_n3,maxnf*Nparticle,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)        

    end subroutine mls_structLoads



subroutine wght_press(pos,ptx,press_probe,probe_inds)

        USE param
        USE geom
        USE mls_param
        use local_arrays, only: pr
        USE mpi_param, only: kstart, kend
        implicit none
        real,dimension(4) :: pxk,ptx,ptxA
        real,dimension(3) :: pos,norp,Wt
        real,dimension(4,4) :: pinvA,invA
        real,dimension(4,nel) :: B
        real,dimension(nel) :: ptxAB(nel)
        real :: press_probe
        real :: Wtx, Wt23
        integer :: inp,ntr,inw,i,j,k,k1, ii, jj, kk
        integer, dimension(3) :: pind_i, pind_o, probe_inds, keul_inds
        
        !-------------Shape function for cell centres (temp. or pressure cells) -------------------------
        
          
        if(probe_inds(3).ge.kstart .and. probe_inds(3).le.kend) then
          
        !-------------FORCING FUNCTION------------------------
        ! volume of a face with a specific marker - thickness taken as average of grid spacing
          
        !WGHT1
        pind_i(1)=probe_inds(1)-1;  pind_o(1)=probe_inds(1)+1 ! xi indices
        pind_i(2)=probe_inds(2)-1;  pind_o(2)=probe_inds(2)+1 ! yj indices
          
        k1  = floor(pos(3)*dx3) + 1
        pind_i(3)=k1-1
        pind_o(3)=k1+1
        
        ! For the spatial grid
        !k_inds = [k1-1,k1,k1]
        
        ! For the Eulerian field
        k1  = modulo(k1-1,n3m)  + 1
        keul_inds = [k1-1,k1,k1+1]


        pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
        inw = 1
          
          
        ! Accumulate A(4,4)   , B(4,27) linear system
        do k=pind_i(3), pind_o(3)
          
            norp(3)=abs(zm(k)-pos(3)) / (wscl / dx3)
            Wt(3) = mls_gaussian( norp(3) , wcon )
          
            do j=pind_i(2), pind_o(2)
          
                norp(2)=abs(ym(j)-pos(2)) / (wscl / dx2)
                Wt(2) = mls_gaussian( norp(2) , wcon )
                Wt23 = Wt(2)*Wt(3)
          
                do i=pind_i(1), pind_o(1)
          
                    norp(1)=abs(xm(i)-pos(1)) / (wscl / dx1)
                    Wt(1) = mls_gaussian( norp(1) , wcon )
                    Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)
          
                    pxk(1)=1.0d0
                    pxk(2)=xm(i)
                    pxk(3)=ym(j)
                    pxk(4)=zm(k)
          
                    call DGEMM('N','T',4,4,1,Wtx,pxk,4,pxk,4, 1.0d0,pinvA,4)
                    B(1:4,inw)=Wtx*pxk(1:4)
          
                    inw = inw + 1
                enddo !end i
            enddo !end j
        enddo !end k
          
            ! calling routine to compute inverse
            ! SPD matrix for uniform grids, we can use Cholesky decomp. instead: dpotrf
            call inverseLU(pinvA,invA)
                  
            !------------------------------------------------------
            ! matrix multiplications for final interpolation
            ! DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
            ! C = alpha * A * B + beta * C
            !---------------Shape function calculation---------------
            call DGEMM('N','N',1,4  ,4,1.0d0,ptx ,1,invA,4,0.0d0,ptxA ,1) 
            call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B   ,4,0.0d0,ptxAB,1) 
          
          
            inw = 1
            press_probe  = 0.0d0 !MLS-interpolated value at probe
          
        do kk=1,3
            k = keul_inds(kk)
          !kk = modulo(k-1,n3m) + 1 !
             do j=pind_i(2), pind_o(2)
               jj = modulo(j-1,n2m) + 1
               do i=pind_i(1), pind_o(1)
                  ii = modulo(i-1,n1m) + 1
                  press_probe = press_probe + pr(ii,jj,k)*ptxAB(inw)
                  inw = inw + 1
              enddo
             enddo
            enddo


          
          endif
end subroutine wght_press
    
          
    
    
subroutine wght_gradP(pos,ptx,gradP,probe_inds)
    USE param
    USE geom
    use local_arrays, only: pr
    USE mls_param
    USE mpi_param, only: kstart, kend
    implicit none
    real, dimension(3), intent(inout) :: gradP !dPdx, dPdy, dPdz at the probe
    real,dimension(nel) :: ddx_PtxAB(nel), ddy_PtxAB(nel), ddz_PtxAB(nel) !Shape function derivatives
    real,dimension(nel) :: Pnel ! Eulerian support domain values of temperature
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
    integer :: ii,jj,kk,nk
    
    integer, dimension(3) :: pind_i, pind_o, probe_inds, keul_inds, k_inds
    
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
    
    ! For the spatial grid
    k_inds = [k1-1,k1,k1]
    
    ! For the Eulerian field
    k1  = modulo(k1-1,n3m)  + 1
    keul_inds = [k1-1,k1,k1+1]
    
      
    pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
    dAdx(1:4,1:4)=0.0d0
    dAdy(1:4,1:4)=0.0d0
    dAdz(1:4,1:4)=0.0d0
    
    inw = 1
      
      
    ! Accumulate A(4,4)   , B(4,27) linear system
    ! Likewise for derivatives of A, B
    !do k=pind_i(3), pind_o(3)
     do nk = 1,3 ! Different to apply periodicty in z
        k = k_inds(nk)
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
    
                !Tnel(inw) = temp(i,j,k) ! Store Eulerian support domain values for summation later
    
                ii = modulo(i-1,n1m) + 1
                jj = modulo(j-1,n2m) + 1
                !kk = modulo(k-1,n3m) + 1
    
                kk = keul_inds(nk)
    
                !Tnel(inw) = temp(ii,jj,kk) ! Store Eulerian support domain values for summation later
                Pnel(inw) = pr(ii,jj,kk) ! Store Eulerian support domain values for summation later
    
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
    
        gradP(1) = sum( ddx_PtxAB * Pnel) !dT dx
        gradP(2) = sum( ddy_PtxAB * Pnel) !dT dy
        gradP(3) = sum( ddz_PtxAB * Pnel) !dT dz
    
      endif
      end subroutine wght_gradP

      subroutine wght_gradU(pos,ptx,gradU,probe_inds)
        USE param
        USE geom
        use local_arrays, only: vx
        USE mls_param
        USE mpi_param, only: kstart, kend
        implicit none
        real, dimension(3), intent(inout) :: gradU !dPdx, dPdy, dPdz at the probe
        real,dimension(nel) :: ddx_PtxAB(nel), ddy_PtxAB(nel), ddz_PtxAB(nel) !Shape function derivatives
        real,dimension(nel) :: Unel ! Eulerian support domain values of temperature
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
        integer :: ii,jj,kk,nk
        
        integer, dimension(3) :: pind_i, pind_o, probe_inds, keul_inds, k_inds
        
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
        
        ! For the spatial grid
        k_inds = [k1-1,k1,k1]
        
        ! For the Eulerian field
        k1  = modulo(k1-1,n3m)  + 1
        keul_inds = [k1-1,k1,k1+1]
        
          
        pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
        dAdx(1:4,1:4)=0.0d0
        dAdy(1:4,1:4)=0.0d0
        dAdz(1:4,1:4)=0.0d0
        
        inw = 1
          
          
        ! Accumulate A(4,4)   , B(4,27) linear system
        ! Likewise for derivatives of A, B
        !do k=pind_i(3), pind_o(3)
         do nk = 1,3 ! Different to apply periodicty in z
            k = k_inds(nk)
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
          
                    norp(1)=abs(pos(1) - xc(i) ) / (wscl / dx1)
                    Wt(1) = mls_gaussian( norp(1) , wcon )
                    Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)
        
                    dWx_dx = mls_gauss_deriv(norp(1),wcon) ! dWx / dr_x
                    dWx_dx = dWx_dx * ( ( pos(1) - xc(i) ) / abs( pos(1) - xc(i)  ) /  (wscl/dx1) ) ! dWx/dx = ( dWx / dr_x ) * ( dr_x / dx )
        
                    pxk(1)=1.0d0
                    pxk(2)=xc(i)
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
        
                    !Tnel(inw) = temp(i,j,k) ! Store Eulerian support domain values for summation later
        
                    ii = modulo(i-1,n1m) + 1
                    jj = modulo(j-1,n2m) + 1
                    !kk = modulo(k-1,n3m) + 1
        
                    kk = keul_inds(nk)
        
                    !Tnel(inw) = temp(ii,jj,kk) ! Store Eulerian support domain values for summation later
                    Unel(inw) = vx(ii,jj,kk) ! Store Eulerian support domain values for summation later
        
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
        
            gradU(1) = sum( ddx_PtxAB * Unel) !dU dx
            gradU(2) = sum( ddy_PtxAB * Unel) !dU dy
            gradU(3) = sum( ddz_PtxAB * Unel) !dU dz
        
          endif
          end subroutine wght_gradU


          subroutine wght_gradV(pos,ptx,gradV,probe_inds)
            USE param
            USE geom
            use local_arrays, only: vy
            USE mls_param
            USE mpi_param, only: kstart, kend
            implicit none
            real, dimension(3), intent(inout) :: gradV !dPdx, dPdy, dPdz at the probe
            real,dimension(nel) :: ddx_PtxAB(nel), ddy_PtxAB(nel), ddz_PtxAB(nel) !Shape function derivatives
            real,dimension(nel) :: Vnel ! Eulerian support domain values of temperature
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
            integer :: ii,jj,kk,nk
            
            integer, dimension(3) :: pind_i, pind_o, probe_inds, keul_inds, k_inds
            
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
            
            ! For the spatial grid
            k_inds = [k1-1,k1,k1]
            
            ! For the Eulerian field
            k1  = modulo(k1-1,n3m)  + 1
            keul_inds = [k1-1,k1,k1+1]
            
              
            pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
            dAdx(1:4,1:4)=0.0d0
            dAdy(1:4,1:4)=0.0d0
            dAdz(1:4,1:4)=0.0d0
            
            inw = 1
              
              
            ! Accumulate A(4,4)   , B(4,27) linear system
            ! Likewise for derivatives of A, B
            !do k=pind_i(3), pind_o(3)
             do nk = 1,3 ! Different to apply periodicty in z
                k = k_inds(nk)
                norp(3)=abs( pos(3) - zm(k) ) / (wscl / dx3)
                Wt(3) = mls_gaussian( norp(3) , wcon )
                dWz_dz = mls_gauss_deriv(norp(3),wcon) ! dWz / dr_z
                dWz_dz = dWz_dz * ( ( pos(3) - zm(k) ) / abs( pos(3) - zm(k)  ) /  (wscl/dx3) ) ! dWz/dz = ( dWz / dr_z ) * ( dr_z / dz )
                do j=pind_i(2), pind_o(2)
              
                    norp(2)=abs( pos(2) - yc(j) ) / (wscl / dx2)
                    Wt(2) = mls_gaussian( norp(2) , wcon )
                    Wt23 = Wt(2)*Wt(3)
              
                    dWy_dy = mls_gauss_deriv(norp(2),wcon) ! dWy / dr_y
                    dWy_dy = dWy_dy * ( ( pos(2) - yc(j) ) / abs( pos(2) - yc(j)  ) /  (wscl/dx2) ) ! dWy/dy = ( dWy / dr_y ) * ( dr_y / dy )
                    do i=pind_i(1), pind_o(1)
              
                        norp(1)=abs(pos(1) - xm(i) ) / (wscl / dx1)
                        Wt(1) = mls_gaussian( norp(1) , wcon )
                        Wtx = Wt(1)*Wt23 !Eq. (3.165) Liu & Gu (2005)
            
                        dWx_dx = mls_gauss_deriv(norp(1),wcon) ! dWx / dr_x
                        dWx_dx = dWx_dx * ( ( pos(1) - xm(i) ) / abs( pos(1) - xm(i)  ) /  (wscl/dx1) ) ! dWx/dx = ( dWx / dr_x ) * ( dr_x / dx )
            
                        pxk(1)=1.0d0
                        pxk(2)=xm(i)
                        pxk(3)=yc(j)
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
            
                        !Tnel(inw) = temp(i,j,k) ! Store Eulerian support domain values for summation later
            
                        ii = modulo(i-1,n1m) + 1
                        jj = modulo(j-1,n2m) + 1
                        !kk = modulo(k-1,n3m) + 1
            
                        kk = keul_inds(nk)
            
                        !Tnel(inw) = temp(ii,jj,kk) ! Store Eulerian support domain values for summation later
                        Vnel(inw) = vy(ii,jj,kk) ! Store Eulerian support domain values for summation later
            
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
            
                gradV(1) = sum( ddx_PtxAB * Vnel) !dV dx
                gradV(2) = sum( ddy_PtxAB * Vnel) !dV dy
                gradV(3) = sum( ddz_PtxAB * Vnel) !dV dz
            
              endif
end subroutine wght_gradV
    

subroutine wght_gradW(pos,ptx,gradW,probe_inds)
                USE param
                USE geom
                use local_arrays, only: vz
                USE mls_param
                USE mpi_param, only: kstart, kend
                implicit none
                real, dimension(3), intent(inout) :: gradW !dPdx, dPdy, dPdz at the probe
                real,dimension(nel) :: ddx_PtxAB(nel), ddy_PtxAB(nel), ddz_PtxAB(nel) !Shape function derivatives
                real,dimension(nel) :: Wnel ! Eulerian support domain values of temperature
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
                integer :: ii,jj,kk,nk
                
                integer, dimension(3) :: pind_i, pind_o, probe_inds, keul_inds, k_inds
                
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
                
                ! For the spatial grid
                k_inds = [k1-1,k1,k1]
                
                ! For the Eulerian field
                k1  = modulo(k1-1,n3m)  + 1
                keul_inds = [k1-1,k1,k1+1]
                
                  
                pinvA(1:4,1:4)=0.0d0 ! Is summed in the loop below
                dAdx(1:4,1:4)=0.0d0
                dAdy(1:4,1:4)=0.0d0
                dAdz(1:4,1:4)=0.0d0
                
                inw = 1
                  
                  
                ! Accumulate A(4,4)   , B(4,27) linear system
                ! Likewise for derivatives of A, B
                !do k=pind_i(3), pind_o(3)
                 do nk = 1,3 ! Different to apply periodicty in z
                    k = k_inds(nk)
                    norp(3)=abs( pos(3) - zc(k) ) / (wscl / dx3)
                    Wt(3) = mls_gaussian( norp(3) , wcon )
                    dWz_dz = mls_gauss_deriv(norp(3),wcon) ! dWz / dr_z
                    dWz_dz = dWz_dz * ( ( pos(3) - zc(k) ) / abs( pos(3) - zc(k)  ) /  (wscl/dx3) ) ! dWz/dz = ( dWz / dr_z ) * ( dr_z / dz )
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
                            pxk(3)=yc(j)
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
                
                            !Tnel(inw) = temp(i,j,k) ! Store Eulerian support domain values for summation later
                
                            ii = modulo(i-1,n1m) + 1
                            jj = modulo(j-1,n2m) + 1
                            !kk = modulo(k-1,n3m) + 1
                
                            kk = keul_inds(nk)
                
                            !Tnel(inw) = temp(ii,jj,kk) ! Store Eulerian support domain values for summation later
                            Wnel(inw) = vz(ii,jj,kk) ! Store Eulerian support domain values for summation later
                
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
                
                    gradW(1) = sum( ddx_PtxAB * Wnel) !dW dx
                    gradW(2) = sum( ddy_PtxAB * Wnel) !dW dy
                    gradW(3) = sum( ddz_PtxAB * Wnel) !dW dz
                
                  endif
end subroutine wght_gradW
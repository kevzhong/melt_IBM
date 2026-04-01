      subroutine solxi_ib(betadx,bcval,field)
!EP   Solves tridiagonal system in i direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      use phasefield
      use ieee_arithmetic !debug
      implicit none
      integer :: jc,kc,ic
      real,intent(in) :: betadx
      real, allocatable, dimension(:) :: amil,apil,acil,fil
      real :: ackl_b

      ! IB variables
      real :: bcval,alpha
      real, intent(in), dimension(1:n1,1:n2,kstart:kend) :: field

      allocate(amil(1:n1))
      allocate(apil(1:n1))
      allocate(acil(1:n1))
      allocate(fil(1:n1))


      ! TRIDIAGONAL SYSTEM TO INVERT
      ! For { u_i-1, u_i, u_i+1 }
      !
      !
      ! [       G     ]           [     ]         [       G     ]           [       1     ]
      ! | -  -------- | u      +  |  1  | u     + | -  -------- | u      =  |   --------- | RHS
      ! [    1 + 2*G  ]  i-1      [     ]  i      [    1 + 2*G  ]  i+1      [    1 + 2*G  ]    i



      ! WORK IN PROGRESS

      
      ! Here, G := betadx

      do kc=kstart,kend
          do jc=1,n2m
             do ic=1,n1m
                ! ackl_b = 1.0/(1.0+2.0*betadx)
                ! apil(ic)=-betadx*ackl_b ! Super-diagonal elements
                ! acil(ic)=1.0d0          ! Diagonal elements
                ! amil(ic)=-betadx*ackl_b ! Sub-diagonal elements
                ! fil(ic)=rhs(ic,jc,kc)*ackl_b !RHS vector

                ! ackl_b = 1.0 / ( 1.0 + (1.0 - phi(ic,jc,kc) ) * 2.0 * betadx )
                ! apil(ic) = - (1.0 - phi(ic,jc,kc) ) * betadx * ackl_b ! Super-diagonal elements
                ! acil(ic) = 1.0d0          ! Diagonal elements
                ! amil(ic) = - (1.0 - phi(ic,jc,kc) ) * betadx * ackl_b ! Sub-diagonal elements

                ! fil(ic) = rhs(ic,jc,kc) * ( 1.0 - phi(ic,jc,kc) )
                ! fil(ic) = fil(ic) + phi(ic,jc,kc) * ( bcval - field(ic,jc,kc) )
                ! fil(ic) = fil(ic) * ackl_b

                alpha =  1.0 - phi(ic,jc,kc) 

                apil(ic) = -betadx*alpha
                acil(ic) = 1.0 + 2.0 * betadx * alpha
                amil(ic) = -betadx * alpha

                fil(ic) = rhs(ic,jc,kc) * alpha
                fil(ic) = fil(ic) + phi(ic,jc,kc) * ( bcval - field(ic,jc,kc) )

                ! !-------- IB volume-penalty interpolation ----------------
                ! fil(ic)  = fil(ic) * ( 1.0 - phi(ic,jc,kc) ) + phi(ic,jc,kc) * ( bcval - field(ic,jc,kc) ) * ackl_b
                ! apil(ic) = apil(ic) * ( 1.0 - phi(ic,jc,kc) )
                ! acil(ic) = acil(ic) * ( 1.0 - phi(ic,jc,kc) )
                ! amil(ic) = amil(ic) * ( 1.0 - phi(ic,jc,kc) )

             enddo
                call tridiag_periodic(amil,acil,apil,fil,1,n1m,n1)
             do ic=1,n1m
                rhs(ic,jc,kc) = fil(ic)  
             enddo
          end do
      end do 


write(*,*) "min(rhs) after solxi: ", minval(rhs)
write(*,*) "max(rhs) after solxi: ", maxval(rhs)
write(*,*) "any NaN  in rhs:      ", any(ieee_is_nan(rhs))
write(*,*) "any Inf  in rhs:      ", any(.not. ieee_is_finite(rhs) .and. .not. ieee_is_nan(rhs))

 
      if(allocated(amil)) deallocate(amil)
      if(allocated(acil)) deallocate(apil)
      if(allocated(apil)) deallocate(acil)
      if(allocated(fil)) deallocate(fil)

      return
      end

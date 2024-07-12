      subroutine solxi(betadx)
!EP   Solves tridiagonal system in i direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      real,intent(in) :: betadx
      real, allocatable, dimension(:) :: amil,apil,acil,fil
      real :: ackl_b

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


      ! Here, G := betadx


      do kc=kstart,kend
          do jc=1,n2m
             do ic=1,n1m
                ackl_b = 1.0/(1.0+2.0*betadx)
                apil(ic)=-betadx*ackl_b ! Super-diagonal elements
                acil(ic)=1.0d0          ! Diagonal elements
                amil(ic)=-betadx*ackl_b ! Sub-diagonal elements
                fil(ic)=rhs(ic,jc,kc)*ackl_b !RHS vector
             enddo
                call tripvmyline(amil,acil,apil,fil,1,n1m,n1)
             do ic=1,n1m
                rhs(ic,jc,kc) = fil(ic)  
             enddo
          end do
      end do 
 
      if(allocated(amil)) deallocate(amil)
      if(allocated(acil)) deallocate(apil)
      if(allocated(apil)) deallocate(acil)
      if(allocated(fil)) deallocate(fil)

      return
      end

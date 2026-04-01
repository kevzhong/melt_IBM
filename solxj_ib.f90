      subroutine solxj_ib(betadx,bcval,field)
!EP   Solves tridiagonal system in j direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      use phasefield
      use ieee_arithmetic
      implicit none
      integer :: jc,kc,ic
      real,intent(in) :: betadx
      real, allocatable, dimension(:) :: amjl,apjl,acjl,fjl
      real :: ackl_b

      ! IB variables
      real :: bcval,alpha
      real, intent(in), dimension(1:n1,1:n2,kstart:kend) :: field


      allocate(amjl(1:n2))
      allocate(apjl(1:n2))
      allocate(acjl(1:n2))
      allocate(fjl(1:n2))

      do kc=kstart,kend
          do ic=1,n1m
             do jc=1,n2m
                ! ackl_b = 1.0/(1.0+2.0*betadx)
                ! apjl(jc)=-betadx*ackl_b
                ! acjl(jc)=1.0d0
                ! amjl(jc)=-betadx*ackl_b
                ! fjl(jc)=rhs(ic,jc,kc)*ackl_b

                !ackl_b = 1.0 / ( 1.0 + (1.0 - phi(ic,jc,kc) ) * 2.0 * betadx )
                !apjl(jc) = - (1.0 - phi(ic,jc,kc) ) * betadx * ackl_b ! Super-diagonal elements
                !acjl(jc) = 1.0d0          ! Diagonal elements
                !amjl(jc) = - (1.0 - phi(ic,jc,kc) ) * betadx * ackl_b ! Sub-diagonal elements

                alpha = 1.0 - phi(ic,jc,kc)
                apjl(jc)=-betadx*alpha
                acjl(jc)= 1.0 + 2.0 * betadx * alpha
                amjl(jc)=-betadx * alpha

                fjl(jc) = rhs(ic,jc,kc)

             enddo
                call tridiag_periodic(amjl,acjl,apjl,fjl,1,n2m,n2)
             do jc=1,n2m
                rhs(ic,jc,kc) = fjl(jc)  
             enddo
          end do
      end do 

      write(*,*) "min(rhs) after solxj: ", minval(rhs)
      write(*,*) "max(rhs) after solxj: ", maxval(rhs)
      write(*,*) "any NaN  in rhs:      ", any(ieee_is_nan(rhs))
      write(*,*) "any Inf  in rhs:      ", any(.not. ieee_is_finite(rhs) .and. .not. ieee_is_nan(rhs))

      if(allocated(amjl)) deallocate(amjl)
      if(allocated(acjl)) deallocate(apjl)
      if(allocated(apjl)) deallocate(acjl)
      if(allocated(fjl)) deallocate(fjl)

      return
      end

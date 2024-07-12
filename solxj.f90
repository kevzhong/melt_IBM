      subroutine solxj(betadx)
!EP   Solves tridiagonal system in j direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      real,intent(in) :: betadx
      real, allocatable, dimension(:) :: amjl,apjl,acjl,fjl
      real :: ackl_b

      allocate(amjl(1:n2))
      allocate(apjl(1:n2))
      allocate(acjl(1:n2))
      allocate(fjl(1:n2))

      do kc=kstart,kend
          do ic=1,n1m
             do jc=1,n2m
                ackl_b = 1.0/(1.0+2.0*betadx)
                apjl(jc)=-betadx*ackl_b
                acjl(jc)=1.0d0
                amjl(jc)=-betadx*ackl_b
                fjl(jc)=rhs(ic,jc,kc)*ackl_b
             enddo
                call tripvmyline(amjl,acjl,apjl,fjl,1,n2m,n2)
             do jc=1,n2m
                rhs(ic,jc,kc) = fjl(jc)  
             enddo
          end do
      end do 


      if(allocated(amjl)) deallocate(amjl)
      if(allocated(acjl)) deallocate(apjl)
      if(allocated(apjl)) deallocate(acjl)
      if(allocated(fjl)) deallocate(fjl)

      return
      end

      subroutine solxi(betadx,solid_mask)
!EP   Solves tridiagonal system in i direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic, dirn
      real,intent(in) :: betadx
      real, allocatable, dimension(:) :: amil,apil,acil,fil
      real :: ackl_b
      real, dimension(3) :: x_grid
      real, dimension(n1m, n2m, kstart:kend), optional, intent(in)  :: solid_mask

      allocate(amil(1:n1))
      allocate(apil(1:n1))
      allocate(acil(1:n1))
      allocate(fil(1:n1))

      do kc=kstart,kend
          do jc=1,n2m
             do ic=1,n1m
               if ( present(solid_mask) ) then
                  if ( abs( solid_mask(ic,jc,kc) ) .gt. 0.0 ) then
                     ! Solid cell, force cell to translational + rotational component
                     apil(ic) = 0.0d0
                     acil(ic)=1.0d0
                     amil(ic) = 0.0d0
                     fil(ic)=rhs(ic,jc,kc)
                  else
                     ackl_b = 1.0/(1.0+2.0*betadx)
                     apil(ic)=-betadx*ackl_b
                     acil(ic)=1.0d0
                     amil(ic)=-betadx*ackl_b
                     fil(ic)=rhs(ic,jc,kc)*ackl_b
                  endif
               else
                  ackl_b = 1.0/(1.0+2.0*betadx)
                  apil(ic)=-betadx*ackl_b
                  acil(ic)=1.0d0
                  amil(ic)=-betadx*ackl_b
                  fil(ic)=rhs(ic,jc,kc)*ackl_b
               endif
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

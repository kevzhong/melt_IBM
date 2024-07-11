      subroutine solxk(q,betadx,solid_mask)
      use param
      use local_arrays, only : rhs
      use mpi_param
      use mpih
      implicit none
      real, intent(inout) :: q(1:n1,1:n2,kstart:kend)
      real, allocatable, dimension(:) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,ic
      real :: betadx,ackl_b
      real,allocatable :: rhst(:,:,:)
      real, dimension(n1m, n2m, kstart:kend), optional, intent(in)  :: solid_mask

      allocate(rhst(1:n3,1:n1,jstart:jend))

      allocate(amkl(1:n3))
      allocate(apkl(1:n3))
      allocate(ackl(1:n3))
      allocate(fkl(1:n3))

      call PackZ_UnpackR(rhs(:,:,kstart:kend),rhst(:,:,jstart:jend))

      do jc=jstart,jend
         do ic=1,n1m
          do kc=1,n3m
          if ( present(solid_mask) ) then
            if ( abs( solid_mask(ic,jc,kc) ) .gt. 0.0 ) then
              amkl(kc)=0.0d0
              ackl(kc)=1.0d0
              apkl(kc)=0.0d0
              fkl(kc)=rhst(kc,ic,jc)
            else
              ackl_b=1.0d0/(1.+2.0*betadx)
              amkl(kc)=-betadx*ackl_b
              ackl(kc)=1.0d0
              apkl(kc)=-betadx*ackl_b
              fkl(kc)=rhst(kc,ic,jc)*ackl_b
            endif
          else
            ackl_b=1.0d0/(1.+2.0*betadx)
            amkl(kc)=-betadx*ackl_b
            ackl(kc)=1.0d0
            apkl(kc)=-betadx*ackl_b
            fkl(kc)=rhst(kc,ic,jc)*ackl_b
          endif
          end do

          call tripvmyline(amkl,ackl,apkl,fkl,1,n3m,n3)

          do kc=1,n3m
            rhst(kc,ic,jc)=fkl(kc)
          end do
         enddo
      end do

      call PackR_UnpackZ(rhst(:,:,jstart:jend),rhs(:,:,kstart:kend))

!     Update velocities

      do kc=kstart,kend
      do jc=1,n2m
      do ic=1,n1m
      q(ic,jc,kc) = q(ic,jc,kc) + rhs(ic,jc,kc)
      enddo
      enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)

      if(allocated(amkl)) deallocate(amkl)
      if(allocated(ackl)) deallocate(apkl)
      if(allocated(apkl)) deallocate(ackl)
      if(allocated(fkl)) deallocate(fkl)

      return
      end

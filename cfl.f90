      subroutine cfl
      use param
      use local_arrays, only: vx,vy,vz
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      integer :: j,k,jp,kp,i,ip
      real :: qcf
      
      cflm=0.00000001d0
                                                                       
      do k=kstart,kend
        kp=k+1
        do j=1,n2m
          jp=j+1
          do i=1,n1m
            ip=i+1
            qcf=( abs((vx(i,j,k)+vx(ip,j,k))*0.5d0*dx1)  &
                 +abs((vy(i,j,k)+vy(i,jp,k))*0.5d0*dx2)  &
                 +abs((vz(i,j,k)+vz(i,j,kp))*0.5d0*dx3))

            cflm = max(cflm,qcf)
      enddo
      enddo
      enddo
            
      call MpiAllMaxRealScalar(cflm)


      return  
      end                                                               

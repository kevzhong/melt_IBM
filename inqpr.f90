!
!     Creates Initial condition (all zero for now)
!

      subroutine inqpr
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, n1m, n2m
      use mpi_param
      use mls_param, only: rad_p
      implicit none
      integer :: ic,jc,kc
      real :: rr

      vx=0.0d0
      vy=0.0d0
      vz=0.0d0

      ! Hard-coded for unit length, r = 0.1 for now
      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        rr = sqrt ( (xm(ic) - 0.5d0)**2 + (ym(jc) - 0.5d0)**2 + (zm(kc) - 0.5d0)**2 )
                        if (rr .ge. rad_p) then !Liquid exterior
                              temp(ic,jc,kc) = 1.0d0
                        else !Solid interior
                              temp(ic,jc,kc) = 0.0d0 
                        endif
                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
      end                                                               

!
!     Creates Initial condition (all zero for now)
!

      subroutine inqpr
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm
      integer :: ic,jc,kc
      real :: rr
      implicit none

      vx=0.0d0
      vy=0.0d0
      vz=0.0d0

      ! Hard-coded for unit length, r = 0.1 for now
      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        rr = sqrt ( (xm(ic) - 0.5d0)**2 + (ym(jc) - 0.5d0)**2 + (zm(jc) - 0.5d0)**2 )
                        if (rr .ge. 0.1d0) then !Liquid exterior
                              temp[ic,jc,kc] = 1.0d0
                        else !Solid interior
                              temp[ic,jc,kc] = 0.0d0 
                        endif !end if
                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
      end                                                               

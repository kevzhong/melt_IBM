!
!     Creates Initial condition (all zero for now)
!

      subroutine inqpr
      use local_arrays, only: vy,vz,vx
      implicit none

      vx=0.0d0
      vy=0.0d0
      vz=0.0d0

      return                                                            
      end                                                               

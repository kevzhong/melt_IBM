!     Creates Initial condition
!
subroutine inqpr
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1
      use mpi_param
      use mls_param, only: rad_p
      implicit none
      integer :: ic,jc,kc
      real :: rr

      vx=0.0d0
      vy=0.0d0
      vz=0.0d0
      temp= Tsol
      
      !For temperature: temp = Tsol in solid interior, otherwise Tliq in liquid exterior

      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        rr = sqrt ( (xm(ic) - 0.5d0*xlen )**2 + (ym(jc) - 0.5d0*ylen )**2 + (zm(kc) - 0.5d0*zlen )**2 )
                        ! Sigmoid fit
                        ! temp(ic,jc,kc) = Tliq - (Tliq - Tsol) / ( 1 + exp(2.0d0 / dx1 * (rr - rad_p)  )

                        !Tanh
                        temp(ic,jc,kc) = Tliq - (Tliq - Tsol)*0.5* (1.0d0 - tanh( (rr - rad_p) * dx1 / 2.0 )  )
                        !if (rr .ge. rad_p) then !Liquid exterior
                        !      temp(ic,jc,kc) = Tliq
                        !endif
                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end                                                               

! Quick and dirty routine to reset the temperature field for uniform/tanh Tsolid, Tliquid
subroutine restart_temperature
      use local_arrays, only: temp
      use param, only: xm, ym, zm, n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1
      use mpi_param
      use mls_param, only: rad_p
      implicit none
      integer :: ic,jc,kc
      real :: rr

      temp= Tsol
            
      !For temperature: temp = Tsol in solid interior, otherwise Tliq in liquid exterior
      
      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        rr = sqrt ( (xm(ic) - 0.5d0*xlen )**2 + (ym(jc) - 0.5d0*ylen )**2 + (zm(kc) - 0.5d0*zlen )**2 )
                        ! Sigmoid fit
                        ! temp(ic,jc,kc) = Tliq - (Tliq - Tsol) / ( 1 + exp(2.0d0 / dx1 * (rr - rad_p)  )
      
                        !Tanh
                        temp(ic,jc,kc) = Tliq - (Tliq - Tsol)*0.5* (1.0d0 - tanh( (rr - rad_p) * dx1 / 2.0 )  )
                        !if (rr .ge. rad_p) then !Liquid exterior
                        !      temp(ic,jc,kc) = Tliq
                        !endif
                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end                      
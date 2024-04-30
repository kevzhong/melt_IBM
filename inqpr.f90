!     Creates Initial condition
!
subroutine inqpr
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1
      use mpi_param
      use mls_param, only: rad_p, pos_CM
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
                        !rr = sqrt ( (xm(ic) - 0.5d0*xlen )**2 + (ym(jc) - 0.5d0*ylen )**2 + (zm(kc) - 0.5d0*zlen )**2 )
                        rr = norm2 (  [ xm(ic),ym(jc), zm(kc) ]  - pos_CM(:,1)  )
                        ! Sigmoid fit
                        ! temp(ic,jc,kc) = Tliq - (Tliq - Tsol) / ( 1 + exp(2.0d0 / dx1 * (rr - rad_p)  )

                        !Tanh
                        temp(ic,jc,kc) = Tliq - (Tliq - Tsol)*0.5* (1.0d0 - tanh( (rr - rad_p) * dx1 / 2.0 )  )
                        !if (rr .ge. rad_p) then !Liquid exterior
                        !      temp(ic,jc,kc) = Tliq
                        !endif
                        
                        ! Ellipsoid
                        !rr =  ( (xm(ic) - 0.5d0*xlen) / 0.20 )**2 + ( (ym(jc) - 0.5d0*ylen)/0.10 )**2 + &
                        !( ( zm(kc) - 0.5d0*zlen) / 0.1 )**2
                        !if (rr .gt. 1.0) then
                        !      temp(ic,jc,kc) = Tliq
                        !endif

                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end   

subroutine inqpr_rotated
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1
      use mpi_param
      use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr, angle
      real, dimension(4) :: quaternion
      real, dimension(3,3) :: rotmat
      real, dimension(3) :: r, rprime

      ! 30 degree CCW rotation about the y axis
      angle = 30.0 
      quaternion(1) = cosd(angle / 2.0)
      quaternion(2) = 0.0
      quaternion(3) = sind(angle / 2.0)
      quaternion(4) = 0.0

      call calc_rot_matrix(quaternion,rotmat)

      ! (Hard-coded) rotated geometry initial condition

      vx=0.0d0
      vy=0.0d0
      vz=0.0d0
      temp= Tsol
      
      !For temperature: temp = Tsol in solid interior, otherwise Tliq in liquid exterior

      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        !rr = sqrt ( (xm(ic) - 0.5d0*xlen )**2 + (ym(jc) - 0.5d0*ylen )**2 + (zm(kc) - 0.5d0*zlen )**2 )

                        r = [ xm(ic) , ym(jc), zm(kc) ] - pos_CM(:,1)
                        rprime = matmul(rotmat, r)

                        rr = (rprime(1) / 0.2)**2  + (rprime(2)/0.1)**2  + ( rprime(3) / 0.1)**2

                        
                        ! Ellipsoid
                        !rr =  ( (xm(ic) - 0.5d0*xlen) / 0.20 )**2 + ( (ym(jc) - 0.5d0*ylen)/0.10 )**2 + &
                        !( ( zm(kc) - 0.5d0*zlen) / 0.1 )**2
                        if (rr .gt. 1.0) then
                              temp(ic,jc,kc) = Tliq
                        endif

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
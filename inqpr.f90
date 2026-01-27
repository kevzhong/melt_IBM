!     Creates Initial condition
!
subroutine ICOND_zeroVelocity
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1
      use mpi_param
      !use mls_param, only: rad_p, pos_CM
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
                        rr = norm2 (  [ xm(ic),ym(jc), zm(kc) ]  - [0.5,0.5,0.5]  )
                        ! Sigmoid fit
                        ! temp(ic,jc,kc) = Tliq - (Tliq - Tsol) / ( 1 + exp(2.0d0 / dx1 * (rr - rad_p)  )

                        !Tanh
                        temp(ic,jc,kc) = Tliq - (Tliq - Tsol)*0.5* (1.0d0 - tanh( (rr - 0.1) * dx1 / 2.0 )  )
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


! Quick and dirty routine to reset the temperature field for uniform/tanh Tsolid, Tliquid
subroutine restart_temperature
      use local_arrays, only: temp
      use param, only: xm, ym, zm, n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1
      use mpi_param
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr, rx_periodic,ry_periodic, rz_periodic 

      temp= Tsol
            
      !For temperature: temp = Tsol in solid interior, otherwise Tliq in liquid exterior
      
      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m

                        rx_periodic =  xm(ic)  - 0.5*xlen
                        ry_periodic =  ym(jc)  - 0.5*ylen
                        rz_periodic =  zm(kc)  - 0.5*zlen

                        rr = sqrt( rx_periodic**2 + ry_periodic**2 + rz_periodic**2 )

                        ! Sigmoid fit
                        ! temp(ic,jc,kc) = Tliq - (Tliq - Tsol) / ( 1 + exp(2.0d0 / dx1 * (rr - rad_p)  )
      
                        !Tanh
                        temp(ic,jc,kc) = Tliq - (Tliq - Tsol)*0.5* (1.0d0 - tanh( (rr - 0.1) * dx1 / 2.0 )  )
                        !if (rr .ge. rad_p) then !Liquid exterior
                        !      temp(ic,jc,kc) = Tliq
                        !endif
                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end                      


subroutine ICOND_TaylorGreen
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, xc,yc,zc,n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1, pi
      use mpi_param
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr

      temp= Tsol
      
      !For temperature: temp = Tsol in solid interior, otherwise Tliq in liquid exterior

      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        !rr = sqrt ( (xm(ic) - 0.5d0*xlen )**2 + (ym(jc) - 0.5d0*ylen )**2 + (zm(kc) - 0.5d0*zlen )**2 )
                        rr = norm2 (  [ xm(ic),ym(jc), zm(kc) ]  - [0.5, 0.5, 0.5]  )
                        ! Sigmoid fit
                        ! temp(ic,jc,kc) = Tliq - (Tliq - Tsol) / ( 1 + exp(2.0d0 / dx1 * (rr - rad_p)  )

                        !Tanh
                        temp(ic,jc,kc) = Tliq - (Tliq - Tsol)*0.5* (1.0d0 - tanh( (rr - 0.1) * dx1 / 2.0 )  )

                        ! Taylor-green vortices IC

                        !if ( (rr .gt. 0.1) ) then

                        vx(ic,jc,kc) = sin(2.0 * pi * xc(ic) / xlen ) * cos(2.0 * pi * ym(jc) / ylen ) &
                                     * cos(2.0 * pi * zm(kc) / zlen )
                        vy(ic,jc,kc) = -cos(2.0 * pi * xc(ic) / xlen ) * sin(2.0 * pi * yc(jc) / ylen ) &
                                     * cos(2.0 * pi * zm(kc) / zlen )
                        vz(ic,jc,kc) = 0.0
                        !endif

                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end


subroutine ICOND_Phasefield
      use param, only: xm, ym, zm, xc,yc,zc,n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1, pi
      use mpi_param
      use local_arrays, only: temp
      use phasefield
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr,rrp
      real :: phip,phic,phic2
      real :: halfthick,rad,x01,x02,y0


      halfthick = 0.2
      rad = 0.1

      x01 = 0.5*xlen - 0.5*halfthick - rad
      x02 = 0.5*xlen + 0.5*halfthick + rad

      y0 = 0.5*ylen

      !For temperature: temp = Tsol in solid interior, otherwise Tliq in liquid exterior

      ! do kc = kstart, kend
      !       do ic = 1, n1m
      !             do jc = 1, n2m

      !                   !! Slabs at top/bottom
      !                   !rr = zm(kc) - 0.5*zlen
      !                   !phi(ic,jc,kc) =  1.0 - 0.5 * ( tanh(0.5 * (rr + 0.4) / pf_eps) - tanh(0.5 * (rr - 0.4) / pf_eps) )


      !                   !! Sphere
      !                   !rr = norm2 (  [ xm(ic),ym(jc), zm(kc) ]  - [0.5, 0.5, 0.5]  )
      !                   !phi(ic,jc,kc) = 0.5* (1.0d0 - tanh( (rr - 0.1) * dx1 / 2.0 )  )


      !                   ! Groove
      !                   rrp = xm(ic) -  0.3
      !                   phip = 0.5*(1.0 - tanh( 0.5*rrp / pf_eps ) )
      !                   rr = norm2 (  [ xm(ic), zm(kc) ]  - [0.5, 0.5]  )
      !                   phic = 0.5*(1.0 + tanh( 0.5*(rr-0.1) / pf_eps ) )

      !                   if ( (rr .le. 0.1) .and. (rrp .lt. 4.0 * pf_eps) ) then
      !                         phi(ic,jc,kc) = phic
      !                   else
      !                         phi(ic,jc,kc) = phip
      !                   endif

      !             enddo !end j
      !       enddo !end i
      ! enddo !end k

            do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m

                        !! Slabs at top/bottom
                        !rr = zm(kc) - 0.5*zlen
                        !phi(ic,jc,kc) =  1.0 - 0.5 * ( tanh(0.5 * (rr + 0.4) / pf_eps) - tanh(0.5 * (rr - 0.4) / pf_eps) )


                        !! Sphere
                        !rr = norm2 (  [ xm(ic),ym(jc), zm(kc) ]  - [0.5, 0.5, 0.5]  )
                        !phi(ic,jc,kc) = 0.5* (1.0d0 - tanh( (rr - 0.1) * dx1 / 2.0 )  )


                        ! Plane slab
                        rrp = abs( xm(ic) - 0.5*xlen )
                        phip = 0.5 * ( 1.0 - tanh(0.5*(rrp - halfthick) / pf_eps) )

                        rr = norm2 (  [ xm(ic) , zm(kc) ]  - [x01, y0 ]  )
                        phic = 0.5 * ( 1.0 + tanh(0.5*(rr - rad) / pf_eps  ))

                        rr = norm2 (  [ xm(ic) , zm(kc) ]  - [x02, y0 ]  )
                        phic2 = 0.5 * ( 1.0 + tanh(0.5*(rr - rad) / pf_eps  ))

                        phi(ic,jc,kc) = phip * phic * phic2

                        ! Map to temperature
                        temp(ic,jc,kc) = (1.0 - phi(ic,jc,kc) ) * Tliq + phi(ic,jc,kc) * Tsol
                  
                  enddo !end j
            enddo !end i
      enddo !end k


      !write(*,*) " sumphi,pf_eps", sumphi, pf_eps
      return                                                            
end

subroutine ICOND_random
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, xc,yc,zc,n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1, pi
      use mpi_param
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rn1, rn2, rn3,maxU

      temp= Tsol
      
      !For temperature: temp = Tsol in solid interior, otherwise Tliq in liquid exterior


      maxU = 5.0

      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m

                        call random_number(rn1)
                        call random_number(rn2)
                        call random_number(rn3)

                        vx(ic,jc,kc) = (rn1 - 0.5) * 2.0 * maxU
                        vy(ic,jc,kc) = (rn2 - 0.5) * 2.0 * maxU
                        vz(ic,jc,kc) = (rn3 - 0.5) * 2.0 * maxU

                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end
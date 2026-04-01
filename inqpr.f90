!     Creates Initial condition
!
subroutine ICOND_zeroVelocity
      use local_arrays, only: vy,vz,vx
      use param
      use mpi_param
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr

      vx=0.0d0
      vy=0.0d0
      vz=0.0d0

      ! ! Solid body velocity
      ! do kc = kstart, kend
      ! do jc = 1, n2m
      ! do ic = 1, n1m
      !       rr = norm2 (  [ xc(ic),ym(jc), zm(kc) ]  - [0.5*xlen, 0.5*ylen, 0.5*zlen]  )
      !       vx(ic,jc,kc) = 0.0 - (0.0 - Usolid)*0.5* (1.0d0 - tanh( 0.5 * (rr - rad_sph) * dx1  )  )
      ! enddo
      ! enddo
      ! enddo            

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
      
      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        vx(ic,jc,kc) = sin(2.0 * pi * xc(ic) / xlen ) * cos(2.0 * pi * ym(jc) / ylen ) &
                                     * cos(2.0 * pi * zm(kc) / zlen )
                        vy(ic,jc,kc) = -cos(2.0 * pi * xm(ic) / xlen ) * sin(2.0 * pi * yc(jc) / ylen ) &
                                     * cos(2.0 * pi * zm(kc) / zlen )
                        vz(ic,jc,kc) = 0.0
                        !endif
                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end

subroutine ICOND_TaylorGreen2D
      use local_arrays, only: vy,vz,vx,temp
      use param, only: xm, ym, zm, xc,yc,zc,n1m, n2m, xlen, ylen, zlen, Tsol, Tliq,dx1, pi
      use mpi_param
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr
      
      do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m

                        ! xy plane
                        vx(ic,jc,kc) = sin(2.0 * pi * xc(ic) / xlen ) * cos(2.0 * pi * ym(jc) / ylen ) 
                        vy(ic,jc,kc) = -cos(2.0 * pi * xm(ic) / xlen ) * sin(2.0 * pi * yc(jc) / ylen )
                        vz(ic,jc,kc) = 0.0

                        ! ! xz plane
                        ! vx(ic,jc,kc) = sin(2.0 * pi * xc(ic) / xlen ) * cos(2.0 * pi * zm(kc) / zlen ) 
                        ! vy(ic,jc,kc) = 0.0
                        ! vz(ic,jc,kc) = -cos(2.0 * pi * xm(ic) / xlen ) * sin(2.0 * pi * zc(kc) / zlen )

                        ! yz plane
                        !vx(ic,jc,kc) = 0.0
                        !vy(ic,jc,kc) = sin(2.0 * pi * yc(jc) / ylen ) * cos(2.0 * pi * zm(kc) / zlen ) 
                        !vz(ic,jc,kc) = -cos(2.0 * pi * ym(jc) / ylen ) * sin(2.0 * pi * zc(kc) / zlen )

                        !endif
                  enddo !end j
            enddo !end i
      enddo !end k
      return                                                            
end

subroutine ICOND_SPHERE
      use param
      use mpi_param
      use local_arrays, only: temp
      use phasefield
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr,rrp
      real :: phip,phic,phic2
      real :: rad,x01,x02,y0




            do kc = kstart, kend
            do jc = 1, n2m
                  do ic = 1, n1m
                        ! Sphere
                        rr = norm2 (  [ xm(ic),ym(jc), zm(kc) ]  - [0.5*xlen, 0.5*ylen, 0.5*zlen]  )

                        phi(ic,jc,kc) = 0.5* (1.0d0 - tanh( 0.5 * (rr - rad_sph) / pf_eps  )  )

                        !if (phi(ic,jc,kc) .ge. 0.99 ) then !clamp
                        !      phi(ic,jc,kc) = 1.0
                        !endif

                        ! Map to temperature
                        temp(ic,jc,kc) = (1.0 - phi(ic,jc,kc) ) * Tliq + phi(ic,jc,kc) * Tsol
                  
                  enddo !end j
            enddo !end i
      enddo !end k


      !write(*,*) " sumphi,pf_eps", sumphi, pf_eps
      return                                                            
end

subroutine ICOND_SPHERE_SHARP
      ! Step-function-like initialization
      use param
      use mpi_param
      use local_arrays, only: temp
      use phasefield
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr,rrp
      real :: phip,phic,phic2
      real :: rad,x01,x02,y0



            do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m
                        ! Sphere
                        rr = norm2 (  [ xm(ic),ym(jc), zm(kc) ]  - [0.5*xlen, 0.5*ylen, 0.5*zlen]  )

                        if (rr .ge. rad_sph) then
                              phi(ic,jc,kc) = 0.0
                        else
                              phi(ic,jc,kc) = 1.0
                        endif
                        !phi(ic,jc,kc) = 0.5* (1.0d0 - tanh( 0.5 * (rr - rad_sph) * dx1  )  )

                        !if (phi(ic,jc,kc) .ge. 0.99 ) then !clamp
                        !      phi(ic,jc,kc) = 1.0
                        !endif

                        ! Map to temperature
                        temp(ic,jc,kc) = (1.0 - phi(ic,jc,kc) ) * Tliq + phi(ic,jc,kc) * Tsol
                  
                  enddo !end j
            enddo !end i
      enddo !end k


      !write(*,*) " sumphi,pf_eps", sumphi, pf_eps
      return                                                            
end



subroutine ICOND_GROOVE
      use param
      use mpi_param
      use local_arrays, only: temp
      use phasefield
      !use mls_param, only: rad_p, pos_CM
      implicit none
      integer :: ic,jc,kc
      real :: rr,rrp
      real :: phip,phic,phic2
      real :: rad,z01,z02,x0
      real :: Fx, Fz, F, gradF, sdf,a,b


      z01 = 0.5*zlen - halfthick
      z02 = 0.5*zlen + halfthick
      x0 = 0.5*xlen

      ! Ellipse semi-axes for convenience
      a = grv_depth
      b = 0.5 * grv_width

            do kc = kstart, kend
            do ic = 1, n1m
                  do jc = 1, n2m

                        ! Plane slab
                        rrp = abs( zm(kc) - 0.5*zlen )
                        phip = 0.5 * ( 1.0 - tanh(0.5*(rrp - halfthick) / pf_eps) )

                        ! LHS ellipse, computed based on tanh(sdf)
                        Fz = zm(kc) - z01
                        Fx = xm(ic) - 0.5*xlen
                        F = (Fx/a)**2 + (Fz/b)**2 - 1.0
                        gradF = 2.0*sqrt( (Fx**2)/(a**4) + (Fz**2)/(b**4) )
                        sdf = F / (gradF + 1e-15)
                        phic = 0.5 * (1.0 + tanh(0.5 * sdf / pf_eps) )

                        ! RHS ellipse, computed based on tanh(sdf)
                        Fz = zm(kc) - z02
                        Fx = xm(ic) - 0.5*xlen
                        F = (Fx/a)**2 + (Fz/b)**2 - 1.0
                        gradF = 2.0*sqrt( (Fx**2)/(a**4) + (Fz**2)/(b**4) )
                        sdf = F / (gradF + 1e-15)
                        phic2 = 0.5 * (1.0 + tanh(0.5 * sdf / pf_eps) )

                        phi(ic,jc,kc) = phip * phic * phic2

                        ! Map to temperature
                        temp(ic,jc,kc) = (1.0 - phi(ic,jc,kc) ) * Tliq + phi(ic,jc,kc) * Tsol
                  
                  enddo !end j
            enddo !end i
      enddo !end k


      !write(*,*) " sumphi,pf_eps", sumphi, pf_eps
      return                                                            
end

! subroutine ICOND_GROOVE
!       use param
!       use mpi_param
!       use local_arrays, only: temp
!       use phasefield
!       !use mls_param, only: rad_p, pos_CM
!       implicit none
!       integer :: ic,jc,kc
!       real :: rr,rrp
!       real :: phip,phic,phic2
!       real :: rad,x01,x02,z0
!       real :: Fx, Fz, F, gradF, sdf,a,b


!       x01 = 0.5*xlen - halfthick
!       x02 = 0.5*xlen + halfthick
!       z0 = 0.5*zlen

!       ! Ellipse semi-axes for convenience
!       a = grv_depth
!       b = 0.5 * grv_width

!             do kc = kstart, kend
!             do ic = 1, n1m
!                   do jc = 1, n2m

!                         ! Plane slab
!                         rrp = abs( xm(ic) - 0.5*xlen )
!                         phip = 0.5 * ( 1.0 - tanh(0.5*(rrp - halfthick) / pf_eps) )

!                         ! LHS ellipse, computed based on tanh(sdf)
!                         Fx = xm(ic) - x01
!                         Fz = zm(kc) - 0.5*zlen
!                         F = (Fx/a)**2 + (Fz/b)**2 - 1.0
!                         gradF = 2.0*sqrt( (Fx**2)/(a**4) + (Fz**2)/(b**4) )
!                         sdf = F / (gradF + 1e-15)
!                         phic = 0.5 * (1.0 + tanh(0.5 * sdf / pf_eps) )

!                         ! RHS ellipse, computed based on tanh(sdf)
!                         Fx = xm(ic) - x02
!                         Fz = zm(kc) - 0.5*zlen
!                         F = (Fx/a)**2 + (Fz/b)**2 - 1.0
!                         gradF = 2.0*sqrt( (Fx**2)/(a**4) + (Fz**2)/(b**4) )
!                         sdf = F / (gradF + 1e-15)
!                         phic2 = 0.5 * (1.0 + tanh(0.5 * sdf / pf_eps) )

!                         phi(ic,jc,kc) = phip * phic * phic2

!                         ! Map to temperature
!                         temp(ic,jc,kc) = (1.0 - phi(ic,jc,kc) ) * Tliq + phi(ic,jc,kc) * Tsol
                  
!                   enddo !end j
!             enddo !end i
!       enddo !end k


!       !write(*,*) " sumphi,pf_eps", sumphi, pf_eps
!       return                                                            
! end

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
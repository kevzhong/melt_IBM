    
    subroutine invtr_phase 
      use param
      use local_arrays, only: rhs
      use mpi_param, only: kstart,kend
      use phasefield
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: udx3
      real    :: dte2,dte3,dcte,dpx33,dte1
      real    :: alre,udx1q,udx2q,udx3q

      alre=al/ren
      udx1q=dx1q
      udx2q=dx2q
      udx3q=dx3q
      udx3 = al*dx3
      
      do kc=kstart,kend
        km=kc-1
        kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do ic=1,n1m
              im=imv(ic)
              ip=ipv(ic)

!   diffusive terms
!
!   x- second derivatives of phi

            dte1=(phi(im,jc,kc)-2.0*phi(ic,jc,kc)+phi(ip,jc,kc))*udx1q

!   y- second derivatives of phi

            dte2=(phi(ic,jm,kc)-2.0*phi(ic,jc,kc)+phi(ic,jp,kc))*udx2q

!   z- second derivatives of phi

            dte3=(phi(ic,jc,kp)-2.0*phi(ic,jc,kc)+phi(ic,jc,km))*udx3q
 
            dcte=dte2+dte3+dte1

            rhs(ic,jc,kc) =  dt * (ga*hphi(ic,jc,kc) + ro*ruphi(ic,jc,kc) ) + &
             (al*D_pf) * dcte * dt

!  updating of the explicit terms

            ruphi(ic,jc,kc)=hphi(ic,jc,kc)
         enddo
       enddo
      enddo

      call solxi(dt*D_pf*0.5d0*al*dx1q )
      call solxj(dt*D_pf*0.5d0*al*dx2q )
      call solxk(phi(1:n1,1:n2,kstart:kend),dt*D_pf*0.5d0*al*dx3q )

 
      return
end subroutine invtr_phase


      subroutine hdnl_phase
      use param
      use local_arrays, only: temp
      use mpi_param, only: kstart,kend
      use phasefield
      !use mls_param,only: dens_ratio
      implicit none
      integer :: ic,jc,kc,ip,im
      real :: nl_phi
      real :: u_iph, u_imh, udx1, h31

      !udx1 = 0.25 * dx1

      do kc=kstart,kend
      do jc=1,n2m
      do ic=1,n1m
        !im=imv(ic)
        !ip=ipv(ic)

        ! Beckermann--Hester formulation
        nl_phi = -D_pf / pf_eps**2 * &
                             phi(ic,jc,kc)  * &
                             ( 1.0 - phi(ic,jc,kc) ) * &
                             (1.0 - 2.0 * phi(ic,jc,kc) + pf_A*(temp(ic,jc,kc) - Tmelt )  )

        hphi(ic,jc,kc) = nl_phi


        ! ! Inclusion of solid advection term, TEST!
        ! u_imh =  ( phi(im,jc,kc) + phi(ic,jc,kc) ) * Usolid
        ! u_iph =  ( phi(ic,jc,kc) + phi(ip,jc,kc) ) * Usolid

        !  h31 = ( u_iph * ( phi(ip,jc,kc) + phi(ic,jc,kc) ) & 
        !        - u_imh * ( phi(ic,jc,kc) + phi(im,jc,kc) ) )*udx1

        ! hphi(ic,jc,kc) = nl_phi - h31

      enddo
      enddo
      enddo
      return
    end subroutine hdnl_phase



    ! Write instantaneous turbulence metrics at current time instant into a text file

      subroutine calcPhaseStats
      use mpih
      use param
      use phasefield,only: phi
      use local_arrays, only: temp
      use stat_arrays
      use mpi_param, only: kstart,kend

      implicit none
      integer :: ic,jc,kc
      real :: vol, tfluid, sumvof
      character(70) namfile

      vol = 0.0
      tfluid = 0.0
      sumvof = 0.0
       do kc=kstart,kend
       do jc=1,n2m
       do ic=1,n1m
            vol = vol + phi(ic,jc,kc)
            sumvof = sumvof + ( 1.0 - phi(ic,jc,kc) )
            tfluid = tfluid + ( 1.0 - phi(ic,jc,kc) ) * temp(ic,jc,kc)
       end do
       end do
       end do
           
      call MpiAllSumRealScalar(vol)
      call MpiAllSumRealScalar(sumvof)
      call MpiAllSumRealScalar(tfluid)

      tfluid = tfluid / sumvof
      vol = vol / (dx1 * dx2 * dx3)

      if(ismaster) then
        write(6,*) "vol:", vol 
        namfile='stringdata/phasefield.txt'
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time, vol, tfluid
      close(92)
      end if

      return
    end subroutine calcPhaseStats

    subroutine compute_psi_from_phi
      use param
      use mpi_param, only: kstart,kend
      use phasefield
      implicit none
      integer :: ic,jc,kc
      real :: small = 1.e-15
      real :: tmp

      ! Calculate (approximate) signed distance function psi from phi

      ! First

      do kc=kstart,kend
      do jc=1,n2m
      do ic=1,n1m

        ! Clip phi for robustness
        tmp = max( phi(ic,jc,kc) , 0.0 )
        tmp = min(tmp, 1.0)

        psi(ic,jc,kc) = pf_eps * log( (tmp + small) / (1.0 - tmp + small ) )

      enddo
      enddo
      enddo
      
      call update_both_ghosts(n1,n2,psi,kstart,kend)

      return
    end subroutine compute_psi_from_phi


    subroutine compute_normals
      use param
      use mpi_param, only: kstart,kend
      use phasefield
      implicit none
      integer :: ic,jc,kc,ip,im, jp,jm, kp,km
      real :: udx1, udx2, udx3
      real :: dpsidx, dpsidy, dpsidz, mag
      real :: small = 1.e-15

      udx1 = 0.5 * dx1
      udx2 = 0.5 * dx2
      udx3 = 0.5 * dx3

      do kc=kstart,kend
        km = kc - 1
        kp = kc + 1
      do jc=1,n2m
        jm = jmv(jc)
        jp = jpv(jc)
      do ic=1,n1m
        im=imv(ic)
        ip=ipv(ic)

        dpsidx = ( psi(ip,jc,kc) - psi(im,jc,kc) ) * udx1
        dpsidy = ( psi(ic,jp,kc) - psi(ic,jm,kc) ) * udx2
        dpsidz = ( psi(ic,jc,kp) - psi(ic,jc,km) ) * udx3

        mag = sqrt(dpsidx**2 + dpsidy**2 + dpsidz**2)

        nhat_x(ic,jc,kc) = dpsidx / (mag + small)
        nhat_y(ic,jc,kc) = dpsidy / (mag + small)
        nhat_z(ic,jc,kc) = dpsidz / (mag + small)

      enddo
      enddo
      enddo

      call update_both_ghosts(n1,n2,nhat_x,kstart,kend)
      call update_both_ghosts(n1,n2,nhat_y,kstart,kend)
      call update_both_ghosts(n1,n2,nhat_z,kstart,kend)

      return

    end subroutine compute_normals


    subroutine compute_curvature
      use param
      use mpi_param, only: kstart,kend
      use phasefield
      implicit none
      integer :: ic,jc,kc,ip,im, jp,jm, kp,km
      real :: udx1, udx2, udx3
      real :: dndx, dndy, dndz
      real :: small = 1.e-15

      udx1 = 0.5 * dx1
      udx2 = 0.5 * dx2
      udx3 = 0.5 * dx3

      do kc=kstart,kend
        km = kc - 1
        kp = kc + 1
      do jc=1,n2m
        jm = jmv(jc)
        jp = jpv(jc)
      do ic=1,n1m
        im=imv(ic)
        ip=ipv(ic)

        dndx = ( nhat_x(ip,jc,kc) - nhat_x(im,jc,kc) ) * udx1
        dndy = ( nhat_y(ic,jp,kc) - nhat_y(ic,jm,kc) ) * udx2
        dndz = ( nhat_z(ic,jc,kp) - nhat_z(ic,jc,km) ) * udx3

        curv(ic,jc,kc) = dndx + dndy + dndz

      enddo
      enddo
      enddo

      call update_both_ghosts(n1,n2,curv,kstart,kend)

      return

    end subroutine compute_curvature

subroutine compute_local_vmelt
      use param
      use mpi_param, only: kstart,kend
      use local_arrays, only: rhs
      use phasefield
      implicit none
      integer :: ic,jc,kc,ip,im, jp,jm, kp,km
      real :: udx1, udx2, udx3
      real :: dphidx, dphidy, dphidz, mag
      real :: small = 1.e-15

      udx1 = 0.5 * dx1
      udx2 = 0.5 * dx2
      udx3 = 0.5 * dx3

      do kc=kstart,kend
        km = kc - 1
        kp = kc + 1
      do jc=1,n2m
        jm = jmv(jc)
        jp = jpv(jc)
      do ic=1,n1m
        im=imv(ic)
        ip=ipv(ic)

        ! rhs(ic,jc,kc) / (al * dt) stores dphi/dt

        dphidx = ( phi(ip,jc,kc) - phi(im,jc,kc) ) * udx1
        dphidy = ( phi(ic,jp,kc) - phi(ic,jm,kc) ) * udx2
        dphidz = ( phi(ic,jc,kp) - phi(ic,jc,km) ) * udx3

        mag = sqrt(dphidx**2 + dphidy**2 + dphidz**2)

        vmelt(ic,jc,kc) = rhs(ic,jc,kc) / (al * dt * mag + small) 
      enddo
      enddo
      enddo

      return

  end subroutine compute_local_vmelt

subroutine minmax_scalars
use param
use local_arrays, only: temp
use mpi_param, only: kstart,kend
use mpih
use phasefield
implicit none
integer :: jc,kc,ic
character(70) namfile

tempmin =  huge(tempmin)
tempmax = -huge(tempmax)

phimin =  huge(phimin)
phimax = -huge(phimax)

 do kc=kstart,kend
  do jc=1,n2m
    do ic=1,n1m
      
     tempmin = min( tempmin, temp(ic,jc,kc ) )
     tempmax = max( tempmax, temp(ic,jc,kc ) )

     phimin = min( phimin, phi(ic,jc,kc ) )
     phimax = max( phimax, phi(ic,jc,kc ) )

   enddo
  enddo
 enddo

call MpiAllMinRealScalar(tempmin)
call MpiAllMinRealScalar(phimin)
call MpiAllMaxRealScalar(tempmax)
call MpiAllMaxRealScalar(phimax)

if (myid .eq. 0) then
  write(*,*) "tempmin, tempmax: ", tempmin, tempmax
  write(*,*) "phimin, phimax: ", phimin, phimax
endif

if(ismaster) then
  namfile='stringdata/minmax.txt'
open(unit=92,file=namfile, Access='append', Status='unknown')
write(92,'(100E15.7)') time, tempmin, tempmax, phimin, phimax
close(92)
end if

end subroutine minmax_scalars


    
    subroutine init_phaseParams
        use param
      !use local_arrays
      use phasefield
      use mpi_param, only: kstart,kend
      use mpih, only: lvlhalo
      !use stat_arrays
      !use mls_local
      use AuxiliaryRoutines
      implicit none
      !integer :: j,k,kc,i
      real :: ste


      pf_eps = 1.0 / dx1

      ste = latHeat / ( cpliquid * (Tliq - Tmelt) )

      D_pf = 6.0 / (5.0 * pf_A * pec * ste ) 
    end 
    
    subroutine init_phaseMemory
      use param
      !use local_arrays
      use phasefield
      use mpi_param, only: kstart,kend
      use mpih, only: lvlhalo
      !use stat_arrays
      !use mls_local
      use AuxiliaryRoutines
      implicit none
      !integer :: j,k,kc,i

      call AllocateReal3DArray(phi,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(nhat_x,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(nhat_y,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(nhat_z,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(psi,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(curv,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vmelt,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)


      call AllocateReal3DArray(hphi,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ruphi,1,n1,1,n2,kstart,kend)
    end

    subroutine dealloc_phaseMem
        use phasefield
      implicit none
      
      if(allocated(phi)) deallocate(phi)
      if(allocated(nhat_x)) deallocate(nhat_x)
      if(allocated(nhat_y)) deallocate(nhat_y)
      if(allocated(nhat_z)) deallocate(nhat_z)
      if(allocated(psi)) deallocate(psi)
      if(allocated(curv)) deallocate(curv)
      if(allocated(vmelt)) deallocate(vmelt)

      if(allocated(hphi)) deallocate(hphi)
      if(allocated(ruphi)) deallocate(ruphi)


      end subroutine dealloc_phaseMem

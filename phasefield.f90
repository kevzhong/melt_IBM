    
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
      integer :: ic,jc,kc
      real :: nl_phi

      do kc=kstart,kend
      do jc=1,n2m
      do ic=1,n1m

        ! Beckermann--Hester formulation
        nl_phi = -D_pf / pf_eps**2 * &
                             phi(ic,jc,kc)  * &
                             ( 1.0 - phi(ic,jc,kc) ) * &
                             (1.0 - 2.0 * phi(ic,jc,kc) + pf_A*(temp(ic,jc,kc) - Tmelt )  )

        hphi(ic,jc,kc) = nl_phi

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
      use stat_arrays
      use mpi_param, only: kstart,kend

      implicit none
      integer :: ic,jc,kc
      real :: vol
      character(70) namfile

      vol = 0.0
       do kc=kstart,kend
       do jc=1,n2m
       do ic=1,n1m
            vol = vol + phi(ic,jc,kc)
       end do
       end do
       end do
           
      call MpiAllSumRealScalar(vol)

      vol = vol / (dx1 * dx2 * dx3)

      if(ismaster) then
        write(6,*) "vol:", vol 
        namfile='stringdata/phasefield.txt'
      open(unit=92,file=namfile, Access='append', Status='unknown')
      write(92,'(100E15.7)') time, vol
      close(92)
      end if

      return
    end subroutine calcPhaseStats
    
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
      
      if(allocated(hphi)) deallocate(hphi)
      if(allocated(ruphi)) deallocate(ruphi)


      end subroutine dealloc_phaseMem

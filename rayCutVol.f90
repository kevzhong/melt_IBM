! Numerical integration scheme for equivalent viscous stresses based on volume-integrated impulses and IBM force
! Underlying equation in Breugem (2012)
! The volume-integration VOF, level-set scheme is as in Kempe & Frohlich (2012)
!
! Variable definitions follow as:
!
!
!           1     /// _
!u_tot =  -----  ///  u dV 
!           V   ///
!
!
!             /// _     _
!r_x_u_1 =   ///  r  x  u  dV 
!           ///
!
! Cell-tagging operations are done as in found in the rayTagging.f90 routines

subroutine convex_hull_q1(ind,inp)
    use mls_param
    use mpih
    use mpi_param, only: kstart,kend
    use local_arrays, only: vx,vy,vz
    use param
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: vof
    real :: alpha
    real :: alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real ,dimension(3) :: Q, C
    real tot_vol_1, r_x_u_1(3)

    tot_vol_1    = 0.
    r_x_u_1(1:3) = 0.

    ! Q stores the cell-centre coordinates and is the query point
    ! C is the control-point and stores the ray origin
  
    ! Random control-point well-outside computational domain
    C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]
  
    do i = ind(1,1),ind(1,2)
      Q(1) = xc(i)
      do j = ind(2,1),ind(2,2)
          Q(2) = ym(j)
        do kk = ind(3,1),ind(3,2)
             Q(3) = zm(kk)
  
             ! Torque/moment arm
             r = pos_CM(1:3,inp) - Q(1:3)

             ! periodic BC
             k = modulo(kk-1,n3m) + 1
  
             if (k.ge.kstart.and.k.le.kend) then 
  
              call rayTagQ(vof,C,Q,inp)
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1

               !if(VOFx(ii,jj,k).lt.1.0)then
               !VOFx(ii,jj,k) = VOFx(ii,jj,k)
               !else
                VOFx(ii,jj,k) = vof
               !end if
                alpha = 1.0 - vof ! alpha := volume of solid = 1 - vof
                alpha_q = alpha*celvol*vx(ii,jj,k)

                ! compute int u over V
                u_tot(1,inp) = u_tot(1,inp) + alpha_q

                ! compute int r x u over V
                r_x_u_1(2) = r_x_u_1(2) + alpha_q * r(3)
                r_x_u_1(3) = r_x_u_1(3) - alpha_q * r(2)
   
                tot_vol_1 = tot_vol_1 + celvol*alpha
               endif
  
        end do
      end do
    end do

    call mpi_globalsum_double_var(tot_vol_1)
    call mpi_barrier(mpi_comm_world,ierr)
    call mpi_globalsum_double_var(u_tot(1,inp))
    call mpi_globalsum_double_arr(r_x_u_1,3)


  u_tot(1,inp)       = u_tot(1,inp) !/ tot_vol_1
  r_x_u_tot(1:3,inp) = r_x_u_tot(1:3,inp) + r_x_u_1 ! / tot_vol_1  ! KZ: volume normalisation not needed
  end subroutine

!=================================
!    q2
!=================================

  subroutine convex_hull_q2(ind,inp)
    use mls_param
    use mpih
    use mpi_param, only: kstart,kend
    use local_arrays, only: vx,vy,vz
    use param
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: vof
    real :: alpha
    real :: alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real ,dimension(3) :: Q, C
    real tot_vol_2, r_x_u_2(3)

    tot_vol_2 = 0.
    r_x_u_2   = 0.

    ! Q stores the cell-centre coordinates and is the query point
    ! C is the control-point and stores the ray origin

    ! Random control-point well-outside computational domain
    C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]
  
    do i = ind(1,1),ind(1,2)
      Q(1) = xm(i)
      do j = ind(2,1),ind(2,2)
          Q(2) = yc(j)
        do kk = ind(3,1),ind(3,2)
             Q(3) = zm(kk)
  
             ! Torque/moment arm
             r = pos_CM(1:3,inp) - Q(1:3)

             ! periodic BC
             k = modulo(kk-1,n3m) + 1
  
             if (k.ge.kstart.and.k.le.kend) then 
  
              call rayTagQ(vof,C,Q,inp)
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1
               !if(VOFx(ii,jj,k).lt.1.0)then
               !VOFy(ii,jj,k) = VOFy(ii,jj,k)
               !else
                VOFy(ii,jj,k) = vof

                alpha = 1.0 - vof ! alpha := volume of solid = 1 - vof
                alpha_q = alpha*celvol*vy(ii,jj,k)

                ! compute int u over V 
                u_tot(2,inp) = u_tot(2,inp) + alpha_q 

                ! compute int r x u over V, r = xgrid - pos_cm
                r_x_u_2(1) = r_x_u_2(1) - alpha_q * r(3)
                r_x_u_2(3) = r_x_u_2(3) + alpha_q * r(1)
   
                tot_vol_2 = tot_vol_2 + celvol *alpha

               !end if
               endif
  
        end do
      end do
    end do

    call mpi_globalsum_double_var(tot_vol_2)
    call mpi_barrier(mpi_comm_world,ierr)
    call mpi_globalsum_double_var(u_tot(2,inp))
    call mpi_globalsum_double_arr(r_x_u_2,3)
  
    u_tot(2,inp)       =  u_tot(2,inp) !/ tot_vol_2
    r_x_u_tot(1:3,inp) = r_x_u_tot(1:3,inp) + r_x_u_2  !/ tot_vol_2 ! KZ: volume normalisation not needed
  
  end subroutine

!=================================
!    q3
!=================================

  subroutine convex_hull_q3(ind,inp)
    use mls_param
    use mpih
    use mpi_param, only: kstart,kend
    use local_arrays, only: vx,vy,vz
    use param
    implicit none
    real, dimension(3)   :: x_grid
    real, dimension(3)   :: r, x_GC
    real, dimension(3,2) :: lim
    real :: vof
    real :: alpha
    real :: alpha_q
    integer :: inp, i,j,k, ii,jj,kk
    integer, dimension(3,2) :: ind
    real ,dimension(3) :: Q, C
    real tot_vol_3, r_x_u_3(3)

    tot_vol_3 = 0.
    r_x_u_3   = 0.
    
    ! Q stores the cell-centre coordinates and is the query point
    ! C is the control-point and stores the ray origin
  
    ! Random control-point well-outside computational domain
    C = [0.0, 0.0, 0.0 ] - [xlen*5.0, ylen*4.0, zlen * 3.0]

    do i = ind(1,1),ind(1,2)
      Q(1) = xm(i)
      do j = ind(2,1),ind(2,2)
          Q(2) = ym(j)
        do kk = ind(3,1),ind(3,2)
             Q(3) = zc(kk)

             ! Torque/moment arm
             r = pos_CM(1:3,inp) - Q(1:3)
  
             ! periodic BC
             k = modulo(kk-1,n3m) + 1
  
             if (k.ge.kstart.and.k.le.kend) then 
  
              call rayTagQ(vof,C,Q,inp)
  
               ii = modulo(i-1,n1m) + 1
               jj = modulo(j-1,n2m) + 1
               !if(VOFz(ii,jj,k).lt.1.0)then
               !VOFz(ii,jj,k) = VOFz(ii,jj,k)
               !else
               VOFz(ii,jj,k) = vof
               !end if
               alpha = 1.0 - vof ! alpha := volume of solid = 1 - vof
               alpha_q = alpha*celvol*vz(ii,jj,k)

                ! compute int u over V 
                u_tot(3,inp) = u_tot(3,inp) + alpha_q 

                ! compute int r x u over V, r = xgrid - pos_cm
                r_x_u_3(1) = r_x_u_3(1) + alpha_q * r(2)
                r_x_u_3(2) = r_x_u_3(2) - alpha_q * r(1)

                tot_vol_3 = tot_vol_3 + celvol*alpha
               endif
  
        end do
      end do
    end do
  
    call mpi_globalsum_double_var(tot_vol_3)
    call mpi_barrier(mpi_comm_world,ierr)
    call mpi_globalsum_double_var(u_tot(3,inp))
    call mpi_globalsum_double_arr(r_x_u_3,3)
  
    u_tot(3,inp)       = u_tot(3,inp) !/ tot_vol_3
    r_x_u_tot(1:3,inp) = r_x_u_tot(1:3,inp) + r_x_u_3  !/ tot_vol_3 ! KZ: volume normalisation not needed
  end subroutine

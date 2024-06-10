      subroutine CalcHITRandomForce
      use param
      use local_arrays, only: forcx,forcy,forcz
      use mpih
      use mpi_param
      use stat_arrays
      implicit none
      integer :: j,k,i,l,n,m
      integer :: kc,jc,ic
      real :: kl,kn,km
      real :: wlow,whigh
      real :: waven
      real :: scalpr, scalpc
      real :: u1, u2
      real, dimension(6,7,7,7) :: fcoefs
      complex :: fcoefsx,fcoefsy,fcoefsz
      complex :: dummy1,dummy3,dummy1b,dummy2b,dummy3b

      wlow=0.01 !kmin=0.0159 -> kmin*2*pi/xlen=0.01
      !wlow = 2.0d0 * pi / xlen ! Lowest wavenumber <--> largest length scale to force set by box size
      !wlow = 1.0e-10 ! Avoid zero wavenumber
      whigh=kf_on_kmin*2*pi/xlen ! Prescribed highest threshold wavenumber <---> smallest length scale to force
      fcoefs=0.0d0

      if(myid.eq.0) then
      do n=1,7
       kn = float(n-4)*2*pi/xlen
       do m=1,7
        km = float(m-4)*2*pi/xlen
        do l=1,4
         kl = float(l-4)*2*pi/xlen
         waven=sqrt(kl**2+km**2+kn**2)
         if((waven.gt.wlow).and.(waven.lt.(whigh))) then
         do i=1,3
            ! Box-Muller transform
          call random_number(u1)
          call random_number(u2)
          bcoefs(2*i-1,l,m,n) = bcoefs(2*i-1,l,m,n)*(1-dt/tl) +  & 
     &     sqrt(2*epsstar*dt/(tl*tl))*sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
          bcoefs(2*i,l,m,n) = bcoefs(2*i,l,m,n)*(1-dt/tl) +  & 
     &     sqrt(2*epsstar*dt/(tl*tl))*sqrt(-2.0*log(u1))*sin(2.0*pi*u2)
         end do
         scalpr=bcoefs(1,l,m,n)*kl + bcoefs(2,l,m,n)*km + &
     &          bcoefs(3,l,m,n)*kn
         scalpc=bcoefs(4,l,m,n)*kl + bcoefs(5,l,m,n)*km + &
     &          bcoefs(6,l,m,n)*kn
         fcoefs(1,l,m,n) =bcoefs(1,l,m,n)-scalpr*kl/(waven**2)
         fcoefs(2,l,m,n) =bcoefs(2,l,m,n)-scalpr*km/(waven**2)
         fcoefs(3,l,m,n) =bcoefs(3,l,m,n)-scalpr*kn/(waven**2)
         fcoefs(4,l,m,n) =bcoefs(4,l,m,n)-scalpc*kl/(waven**2)
         fcoefs(5,l,m,n) =bcoefs(5,l,m,n)-scalpc*km/(waven**2)
         fcoefs(6,l,m,n) =bcoefs(6,l,m,n)-scalpc*kn/(waven**2)
         end if
        end do
       end do
      end do

!RO     Complex conjugate
        do k=1,7
        kc = 8-k
        do j=1,7
        jc = 8-j
        do i=5,7
         ic = 8-i
         fcoefs(1,i,j,k) =fcoefs(1,ic,jc,kc)
         fcoefs(2,i,j,k) =fcoefs(2,ic,jc,kc)
         fcoefs(3,i,j,k) =fcoefs(3,ic,jc,kc)
         fcoefs(4,i,j,k)=-fcoefs(4,ic,jc,kc)
         fcoefs(5,i,j,k)=-fcoefs(5,ic,jc,kc)
         fcoefs(6,i,j,k)=-fcoefs(6,ic,jc,kc)
        end do
        end do
        end do

!RO     Special case i=4 (k=0)
        i=4 
        do k=1,7
         do j=1,7
         fcoefs(4,i,j,k)=0.
         fcoefs(5,i,j,k)=0.
         fcoefs(6,i,j,k)=0.
         end do
        end do


        do k=1,7
         kc=8-k
         do j=5,7
         jc=8-j
         fcoefs(1,i,j,k)=fcoefs(1,i,jc,kc)
         fcoefs(2,i,j,k)=fcoefs(2,i,jc,kc)
         fcoefs(3,i,j,k)=fcoefs(3,i,jc,kc)
         end do
        end do

        j=4
        do k=5,7
         kc=8-k
         fcoefs(1,i,j,k)=fcoefs(1,i,j,kc)
         fcoefs(2,i,j,k)=fcoefs(2,i,j,kc)
         fcoefs(3,i,j,k)=fcoefs(3,i,j,kc)
        end do

!RO    0 constant force

       fcoefs(:,4,4,4) = 0.0d0

       end if

!RS    Check the data range that actually needs to be send in this routine
       call mpi_globalsum_double_forc(fcoefs)

  
       forcx=0.0d0 ! Initial force as we are using these matrix elements to sum
       do n=2,6 ! The ring values (n=1,n=7) are zero, so exclude from summation
       do m=2,6 ! The ring values (m=1,m=7) are zero, so exclude from summation
       do l=2,6 ! The ring values (l=1,l=7) are zero, so exclude from summation
        fcoefsx=cmplx(fcoefs(1,l,m,n),fcoefs(4,l,m,n))
         do k=kstart,kend
         dummy1=term3a(k,n) ! use pre-allocated axial terms
         do j=1,n2m
          dummy1b=dummy1*term2a(j,m)!*term3a(k,n) ! use pre-allocated radial terms
          do i=1,n1m
           forcx(i,j,k) = forcx(i,j,k)+real(fcoefsx*dummy1b*term1b(i,l)) ! sum force
          end do 
         end do
         end do
        end do
        end do
        end do

!      Seperate loops per component, seem to benefit data locality, could be system dependent
       forcy=0.0d0 ! Initial force as we are using these matrix elements to sum
       do n=2,6 ! The ring values (n=1,n=7) are zero, so exclude from summation
       do m=2,6 ! The ring values (m=1,m=7) are zero, so exclude from summation
       do l=2,6 ! The ring values (l=1,l=7) are zero, so exclude from summation
        fcoefsy=cmplx(fcoefs(2,l,m,n),fcoefs(5,l,m,n))
         do k=kstart,kend
         dummy1=term3a(k,n) ! use pre-allocated axial terms
         do j=1,n2m 
         dummy2b=dummy1*term2b(j,m) ! use pre-allocated radial terms
          do i=1,n1m
           forcy(i,j,k)=forcy(i,j,k)+real(fcoefsy*dummy2b*term1a(i,l)) ! sum force
          end do
         end do
         end do
       end do
       end do
       end do

!      Seperate loops per component, seem to benefit data locality, could be system dependent
       forcz=0.0d0 ! Initial force as we are using these matrix elements to sum
       do n=2,6 ! The ring values (n=1,n=7) are zero, so exclude from summation
       do m=2,6 ! The ring values (m=1,m=7) are zero, so exclude from summation
       do l=2,6 ! The ring values (l=1,l=7) are zero, so exclude from summation
        fcoefsz=cmplx(fcoefs(3,l,m,n),fcoefs(6,l,m,n))
         do k=kstart,kend
         dummy3=term3b(k,n) ! use pre-allocated axial terms
         do j=1,n2m
          dummy3b=dummy3*term2a(j,m) ! use pre-allocated radial terms
          do i=1,n1m
           forcz(i,j,k)=forcz(i,j,k)+real(fcoefsz*dummy3b*term1a(i,l)) ! sum force
          end do
         end do
         end do
      end do
      end do
      end do

      return
      end


      subroutine InitRandomForce
      use mpih
      use mpi_param, only: kstart,kend
      use param
      use hdf5
      IMPLICIT NONE
      integer i,j,k,l,m,n
      integer hdf_error, ndims
      integer(HID_T) :: file_id
      integer(HID_T) :: dset_coef
      integer(HID_T) :: dspace_coef
      integer(HSIZE_T) :: dims(4)
      character(40) filnambc
      real kl, km, kn
      complex :: im=(0.,1.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Pre-allocation of Geometrical terms used in the forcing (calcforc.F90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! Azimuthal terms
      do i=1,n1m
       do l=1,7
       kl = float(l-4)*2*pi
       term1a(i,l)=exp(im*kl*xm(i)/xlen)
       term1b(i,l)=exp(im*kl*xc(i)/xlen)
       enddo
      enddo

!! Radial terms
      do j=1,n2m
       do m=1,7
       km = float(m-4)*2*pi
       term2a(j,m)=exp(im*km*ym(j)/ylen)
       term2b(j,m)=exp(im*km*yc(j)/ylen)
       enddo
      enddo

!! Axial terms
      do k=kstart,kend
       do n=1,7
       kn = float(n-4)*2*pi
       term3a(k,n)=exp(im*kn*zm(k)/zlen)
       term3b(k,n)=exp(im*kn*zc(k)/zlen)
       end do
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End pre-allocation of Geometrical terms used in the forcing (calcforc.F90)
end

subroutine WriteRandForcCoef
use mpih
use param
use hdf5
IMPLICIT NONE
integer hdf_error, ndims
character(40) filnambc
integer(HID_T) :: file_id
integer(HID_T) :: dset_coef
integer(HID_T) :: dspace_coef
integer(HSIZE_T) :: dims_coef(4)


filnambc = 'continuation/continua_bcoefs.h5'

if(ismaster) then 

ndims=4
dims_coef(1)=6
dims_coef(2)=7
dims_coef(3)=7
dims_coef(4)=7

call h5fcreate_f(filnambc,H5F_ACC_TRUNC_F, file_id, hdf_error)

call h5screate_simple_f(ndims, dims_coef, dspace_coef, hdf_error)

call h5dcreate_f(file_id, 'bcoefs', H5T_NATIVE_DOUBLE, &
                dspace_coef, dset_coef, hdf_error)

call h5dwrite_f(dset_coef, H5T_NATIVE_DOUBLE, bcoefs, &
       dims_coef,hdf_error)

call h5dclose_f(dset_coef, hdf_error)
call h5sclose_f(dspace_coef, hdf_error)

call h5fclose_f(file_id, hdf_error)

endif


return                                    
end      

subroutine CalcABC_HITForce
      use param
      use local_arrays, only: forcx,forcy,forcz
      use mpih
      use mpi_param
      use stat_arrays
      implicit none
      integer :: j,k,i
      real :: C

      C = 1.0

      do k=kstart,kend
            do j=1,n2m
                  do i=1,n1m
                        forcx(i,j,k) = C * sin(2.0 * pi * zm(k) / zlen ) + C * cos(2.0 * pi * ym(j) / ylen )
                        forcy(i,j,k) = C * sin(2.0 * pi * xm(i) / xlen ) + C * cos(2.0 * pi * zm(k) / zlen )
                        forcz(i,j,k) = C * sin(2.0 * pi * ym(j) / ylen ) + C * cos(2.0 * pi * xm(i) / xlen )

                        !forcx(i,j,k) = C * sin(2.0 * pi * zm(k) / zlen ) + C * cos(2.0 * pi * ym(j) / ylen )
                        !forcy(i,j,k) = C * sin(2.0 * pi * xm(i) / xlen ) + C * cos(2.0 * pi * zm(k) / zlen )
                        !forcz(i,j,k) = C * sin(2.0 * pi * ym(j) / ylen ) + C * cos(2.0 * pi * xm(i) / xlen )
                  enddo
            enddo
      enddo


      return
end

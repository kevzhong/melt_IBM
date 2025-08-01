!-----------------------------------------------------
!     Auxilliary particle statistics
!-----------------------------------------------------

subroutine write_partpos
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(Nparticle*3) :: pos
  integer                      :: i,idx,inp

  character(70) namfile


  if (myid.eq.0) then

  namfile='stringdata/part_pos.txt'

  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  write(43,'(100E15.7)')time,dt, pos_cm 

  close(43)
  end if

end subroutine write_partpos

subroutine write_partvel 
  use param
  use mls_param
  use mpih

  IMPLICIT none

  integer i,idx,inp
  real vel(3,Nparticle)

  character(70) namfile


  ! first of all, fill up velocity arrays

  if (myid.eq.0) then

    namfile='stringdata/part_vel.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
      write(43,'(100E15.7)') time,dt, vel_cm 
    close(43)
  end if

end subroutine write_partvel 

! subroutine write_shortdist 
!   use param
!   use mls_param
!   use mpih
!   use coll_mod

!   IMPLICIT none

!   integer i,idx,inp

!   character(70) namfile


!   ! first of all, fill up velocity arrays

!   if (myid.eq.0) then

!     namfile='flowmov/shortdist.txt'

!     open(unit=43,file=namfile,Access = 'append', Status='unknown')
!       write(43,'(100E15.7)') time,dt, short_dist 
!     close(43)
!   end if

! end subroutine 




subroutine write_quat 
  use param
  use mls_param
  use mpih

  IMPLICIT none

  integer i,idx,inp
  real vel(3,Nparticle)

  character(70) namfile


  ! first of all, fill up velocity arrays

  if (myid.eq.0) then


    namfile='flowmov/quat.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
      write(43,'(200E15.7)') time,dt, quat 
    close(43)
  end if

end subroutine write_quat 

subroutine write_partrot 
  use param
  use mls_param
  use mpih

  implicit none

  integer i,idx,inp

  character(70) namfile


  ! first of all, fill up velocity arrays

  if (myid.eq.0) then

    namfile='stringdata/part_rot.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
    write(43,'(100E15.7)')time,dt, omega_c 
   
    close(43)
  end if

end subroutine write_partrot 


subroutine set_particle_array_sizes
  use param
  use mpih
  use mls_param
  implicit none

  open(109,file=gtsfx)

    read(109,*) maxnv, maxne, maxnf

  close(109)

  if(myid.eq.0) then
    write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(6,*)'reading mesh', gtsfx
    write(6,*)'N vertices, N edges, N faces', maxnv, maxne, maxnf
  end if


end subroutine

subroutine set_connectivity
  use param
  use mpih
  use mls_param
  implicit none
  integer :: i,j,v1,v2,v3,e1,e2,e3
  integer :: nv,ne,nf, inp
  integer :: dummyInd

    inp = 1

    open(109,file=gtsfx)

      read(109,*) nv, ne, nf

      do i=1,maxnv
        read(109,*)xyz0(1,i),xyz0(2,i),xyz0(3,i)
      end do

      !n_edge_of_vert(1:maxnv)=0

      do i=1,maxne
        read(109,*)v1,v2
        vert_of_edge(1,i,inp)=v1
        vert_of_edge(2,i,inp)=v2
      enddo

      do i=1,maxnf
          read(109,*)edge_of_face(1,i,inp),edge_of_face(2,i,inp), edge_of_face(3,i,inp)
      end do

      do i=1,maxnf
        e1=edge_of_face(1,i,inp)
        e2=edge_of_face(2,i,inp)

        if (vert_of_edge(2,e1,inp).eq.vert_of_edge(1,e2,inp)) then
           v1=vert_of_edge(1,e1,inp)
           v2=vert_of_edge(2,e1,inp)
           v3=vert_of_edge(2,e2,inp)
        elseif(vert_of_edge(2,e1,inp).eq.vert_of_edge(2,e2,inp)) then
           v1=vert_of_edge(1,e1,inp)
           v2=vert_of_edge(2,e1,inp)
           v3=vert_of_edge(1,e2,inp)
        elseif(vert_of_edge(1,e1,inp).eq.vert_of_edge(1,e2,inp)) then
           v1=vert_of_edge(2,e1,inp)
           v2=vert_of_edge(1,e1,inp)
           v3=vert_of_edge(2,e2,inp)
        else 
           v1=vert_of_edge(2,e1,inp)
           v2=vert_of_edge(1,e1,inp)
           v3=vert_of_edge(1,e2,inp)
        endif 

           vert_of_face(1,i,inp)=v1
           vert_of_face(2,i,inp)=v2
           vert_of_face(3,i,inp)=v3
      end do

    close(109)  !VS   Completed reading the gts file

    ! Check - for faces and edges
    face_of_edge(1:2,1:maxne,inp)=0
    do i=1,maxnf
       e1=edge_of_face(1,i,inp)
       e2=edge_of_face(2,i,inp)
       e3=edge_of_face(3,i,inp)
       if (face_of_edge(1,e1,inp).eq.0) then
          face_of_edge(1,e1,inp)=i
       elseif (face_of_edge(2,e1,inp).eq.0) then
          face_of_edge(2,e1,inp)=i
       else
          write(*,*)'Edge error', i,e1,e2,e3
       endif
       if (face_of_edge(1,e2,inp).eq.0) then
          face_of_edge(1,e2,inp)=i
       elseif (face_of_edge(2,e2,inp).eq.0) then
          face_of_edge(2,e2,inp)=i
       else
          write(*,*)'Edge error'
       endif
       if (face_of_edge(1,e3,inp).eq.0) then
          face_of_edge(1,e3,inp)=i
       elseif (face_of_edge(2,e3,inp).eq.0) then
          face_of_edge(2,e3,inp)=i
       else
          write(*,*)'Edge error'
       endif
    enddo 

    !Check
    do i=1,maxne
       if (face_of_edge(1,i,inp).eq.face_of_edge(2,i,inp)) then
          write(*,*)'Error on edges '
       endif
    enddo

    ! Copy initial connectivity arrays for every particle
    if (Nparticle .gt. 1) then
      do inp = 2,Nparticle
        face_of_edge(:,:,inp) = face_of_edge(:,:,1)
        vert_of_edge(:,:,inp) = vert_of_edge(:,:,1)

        edge_of_face(:,:,inp) = edge_of_face(:,:,1)
        vert_of_face(:,:,inp) = vert_of_face(:,:,1)
      enddo
    endif

end subroutine


subroutine writePind
  use param
  use mls_param
  use mpih

  IMPLICIT none

  character(70) namfile

  if (myid.eq.0) then

    namfile='flowmov/partPind.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')

      write(43,'(1E15.7,200I4.1)') time,pind1
    close(43)
  end if

end subroutine writePind

subroutine writePPpartpos
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(Nparticle*3) :: pos
  integer                      :: i,idx,inp

  character(70) namfile


  if (myid.eq.0) then

  namfile='flowmov/partPPpos.txt'

  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  write(43,'(100E15.7)')time,dt, pos_cm

  close(43)
  end if

end subroutine writePPpartpos

subroutine writePPquat
  use param
  use mls_param
  use mpih

  IMPLICIT none

  integer i,idx,inp
  real vel(3,Nparticle)

  character(70) namfile


  ! first of all, fill up velocity arrays

  if (myid.eq.0) then


    namfile='flowmov/partPPquat.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
      write(43,'(200E15.7)') time,dt, quat
    close(43)
  end if

end subroutine writePPquat

subroutine writeInertTens
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(Nparticle*3) :: pos
  integer                      :: i,idx,inp

  character(70) namfile


  if (myid.eq.0) then

  namfile='stringdata/inertTensor.txt'
 !KZ: note hard-coded single particle for now
  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  write(43,'(100E15.7)') InertTensor(1,1,1), InertTensor(1,2,1), InertTensor(1,3,1) , &
                        InertTensor(2,2,1), InertTensor(2,3,1), InertTensor(3,3,1) 

  close(43)
  end if
end subroutine writeInertTens

subroutine writeStructLoads
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(Nparticle*3) :: pos
  integer                      :: i,idx,inp

  character(70) namfile


if (IMLSSTR .eq. 1) then
  if (myid.eq.0) then
  namfile='stringdata/mlsLoads.txt'

  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  write(43,'(100E15.7)') Fp, Ftau, torq_p, torq_tau

  close(43)
  end if
endif
end subroutine writeStructLoads

subroutine writePartVol
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(Nparticle*3) :: pos
  integer                      :: i,idx,inp

  character(70) namfile


  if (myid.eq.0) then

  namfile='stringdata/part_vol.txt'
 !KZ: note hard-coded single particle for now
  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  !write(43,'(100E15.7)')time, Volume(1), Surface(1), maxval( pack(skewness(:,:) , .not. isGhostFace(:,:)  ) ) , vel_CM(3,1) 
  write(43,'(100E15.7)')time, Volume(1), Surface(1)

  close(43)
  end if
end subroutine writePartVol

subroutine writeTriMeshStats
  use param
  use mls_param
  use mpih

  IMPLICIT none

  integer                      :: i,idx,inp

  character(70) namfile


  if (myid.eq.0) then

  namfile='stringdata/triStats.txt'
 !KZ: note hard-coded single particle for now

  ! Time, Nv, Ne, Nf, min(skew), max(skew), min(elength), max(elength), min(atri), max(atri)
  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  !write(43,'(100E15.7)')time,&
  write(43, '(F10.6, 3(I10), 6(F10.6))')time,&
  count(isGhostVert(:,1) .eqv. .false.),count(isGhostEdge(:,1) .eqv. .false.),count(isGhostFace(:,1) .eqv. .false.),&
  minval( pack(skewness(:,:) , .not. isGhostFace(:,:)  ) ),   maxval( pack(skewness(:,:) , .not. isGhostFace(:,:)  ) ) , &
  minval( pack(eLengths(:,:) , .not. isGhostEdge(:,:)  ) ),   maxval( pack(eLengths(:,:) , .not. isGhostEdge(:,:)  ) ) , &
  minval( pack(sur(:,:) , .not. isGhostFace(:,:)  ) ),   maxval( pack(sur(:,:) , .not. isGhostFace(:,:)  ) ) 

  close(43)
  end if
end subroutine writeTriMeshStats

subroutine writeClock
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(Nparticle*3) :: pos
  integer                      :: i,idx,inp

  character(70) namfile

  !wtime_vof = 0.
  !eul_solve_wtime = 0.
  !mls_wtime = 0.
  !pressure_wtime = 0.
  !hit_wtime = 0.


  if (myid.eq.0) then

  namfile='stringdata/clock.txt'

  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  write(43,'(100E15.7)')time, wtime_vof, eul_solve_wtime, mls_wtime, pressure_wtime, hit_wtime,wtime_total

  close(43)
  end if
end subroutine writeClock



subroutine write_tecplot_geom
  use mpih
  use mpi_param
  use param
  use mls_param
  
  character(70) namfile,namfi2
  character(7) ipfi
  character*2 ibod
  real tprfi
  integer itime
  tprfi = 1/tframe
  itime=nint(time*tprfi)
  write(ipfi,82)itime
  82 format(i7.7)
  98 format(i2.2)
  
  if(ismaster)then
      do inp=1,Nparticle
  
        write(ibod,98) inp
        write(ipfi,82) itime
  
        namfi2='flowmov/geom_'//ibod//'_'//ipfi
  
        !call write_geom (maxnv,maxnf,xyzv(:,:,inp),sur(:,inp),namfi2)
        call write_VertGeom (inp,maxnv,maxnf,xyzv(:,:,inp),vmelt(:,:,inp),qw_oVert(:,inp),qw_iVert(:,inp),&
        vert_nor(:,:,inp),namfi2,isGhostFace(:,inp),isGhostVert(:,inp))

      end do
  end if
  
  end subroutine write_tecplot_geom
   !---------------------------------------------------------------------------------------------
  subroutine write_VertGeom (inp,nv,nf,xyz,vmelt,qw_o,qw_i,vert_nor,filename,isGhostFace,isGhostVert)
    use param
    use mpih
    use mls_param, only: vert_of_face!, isGhostFace, isGhostVert
    implicit none
    character(70) filename,geotecfile
    integer i,nv,nf,inp
    real, dimension (3,nv) :: xyz , vert_nor 
    real, dimension (nv) :: qw_o, qw_i
    real, dimension(3,nv) :: vmelt
    logical, dimension(nv) :: isGhostVert
    logical, dimension(nf) :: isGhostFace
    integer, dimension(nv) :: vert_mask
    integer :: vcnt
    real tprfi
    integer itime
    character(7) ipfi


    tprfi = 1/tframe
    itime=nint(time*tprfi)
    write(ipfi,82)itime
 82 format(i7.7)

   

      geotecfile=trim(filename)//'.dat'
!        write(*,*)' Write file ',trim(geotecfile)

      open(11,file=geotecfile)

      write(11,*)'TITLE = "Geo"'
      !write(11,*)'VARIABLES = X Y Z Vx Vy Vz'
      write(11,*)'VARIABLES = X Y Z vmelt qw_liquid qw_solid nhat_x nhat_y nhat_z'
  !    write(11,*)'ZONE T="DOMAIN 0", N=',nv,' E=',nf,' F=FEBLOCK, ET=TRIANGLE'
      write(11,*)'ZONE T="FETri" N=',count(isGhostVert .eqv. .false.),' E=',count(isGhostFace .eqv. .false.),' ZONETYPE=FETriangle'
      ! write(11,*)'ZONE T=FETri N=',nvc,' E=',ntri,' ZONETYPE=FETriangle'
      write(11,*)'DATAPACKING=POINT                                       '
      !write(11,*)'VARLOCATION=([4-6]=CELLCENTERED)' ! KZ: lists 4-6 are centroid data, can specify as nodal for vertices
      !write(11,*)'VARLOCATION=([4-9]=NODAL)' 

      vcnt = 1

      do i=1,nv
        if (isGhostVert(i) .eqv. .false.) then
        write(11,*)xyz(1,i), xyz(2,i), xyz(3,i), norm2(vmelt(1:3,i)), qw_o(i), qw_i(i), vert_nor(1,i), vert_nor(2,i), vert_nor(3,i)

        ! For tracking of cumulative non-ghost vertices
        vert_mask(i) = vcnt
        vcnt = vcnt + 1
        endif
      end do

      do i=1,nf
        if (isGhostFace(i) .eqv. .false.) then
          write(11,*) vert_mask( vert_of_face(1,i,inp) ), vert_mask( vert_of_face(2,i,inp) ), vert_mask( vert_of_face(3,i,inp) )
        endif
      end do


      close(11)
      return
      end subroutine write_VertGeom

subroutine writeRemeshStats(n_ecol, n_erel, DV_residual, maxdrift,minA)
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real :: DV_residual, maxdrift,minA
  integer :: n_ecol, n_erel
  character(70) namfile


  if (myid.eq.0) then

  namfile='stringdata/remeshStats.txt'
 !KZ: note hard-coded single particle for now
  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  !write(43,'(100E15.7)')vol_pre, vol_melt, vol_coarse, vol_smooth 
  write(43,'(2I6, 3E15.7)')n_ecol, n_erel ,DV_residual, maxdrift*dx1,minA
  close(43)
  end if
end subroutine writeRemeshStats

subroutine write_triDebug
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real :: DV_residual, maxdrift
  integer :: n_ecol, n_erel,i
  character(70) namfile


  if (myid.eq.0) then

  namfile='flowmov/tribar_debug.txt'
 !KZ: note hard-coded single particle for now
  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  !write(43,'(100E15.7)')vol_pre, vol_melt, vol_coarse, vol_smooth 
  do i = 1,maxnf
    write(43,'(3E15.7)') tri_bar(1,i,1), tri_bar(2,i,1), tri_bar(3,i,1)
  enddo
  close(43)
  end if
end subroutine write_triDebug


subroutine writePrincAxes(xCOM,Ib_ij)
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(3) :: xCOM
  real, dimension(3,3) :: Ib_ij
  character(70) namfile,namfi2
  character(7) ipfi
  character*2 ibod
  real tprfi
  integer itime, inp


  if(ismaster)then

  inp = 1

  tprfi = 1/tframe
  itime=nint(time*tprfi)
  write(ipfi,82)itime
  82 format(i7.7)
  98 format(i2.2)

  write(ibod,98) inp
  write(ipfi,82) itime
  
  namfi2='flowmov/partAxes_'//ibod//'_'//ipfi
  namfi2=trim(namfi2)//'.txt'

 !KZ: note hard-coded single particle for now
  open(unit=43,file=namfi2,Access = 'append', Status='unknown')
  !write(43,'(100E15.7)')vol_pre, vol_melt, vol_coarse, vol_smooth 
  write(43,*)'X  Y  Z  e11  e21  e31  e12  e22  e32  e13  e23  e33'
  write(43,'(12E15.7)')xCOM(1), xCOM(2), xCOM(3), Ib_ij(1,1) , Ib_ij(2,1), Ib_ij(3,1) , &
                                                 Ib_ij(1,2) , Ib_ij(2,2), Ib_ij(3,2) , &
                                                 Ib_ij(1,3) , Ib_ij(2,3), Ib_ij(3,3) 
  close(43)
  end if
end subroutine writePrincAxes
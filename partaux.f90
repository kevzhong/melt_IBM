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

  namfile='flowmov/part_pos.txt'

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

    namfile='flowmov/part_vel.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
      write(43,'(100E15.7)') time,dt, vel_cm 
    close(43)
  end if

end subroutine write_partvel 

subroutine write_shortdist 
  use param
  use mls_param
  use mpih
  use coll_mod

  IMPLICIT none

  integer i,idx,inp

  character(70) namfile


  ! first of all, fill up velocity arrays

  if (myid.eq.0) then

    namfile='flowmov/shortdist.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
      write(43,'(100E15.7)') time,dt, short_dist 
    close(43)
  end if

end subroutine 




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

    namfile='flowmov/part_rot.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
    write(43,'(100E15.7)')time,dt, omega_b 
   
    close(43)
  end if

end subroutine write_partrot 

subroutine write_tail_head 
  use param
  use mls_param
  use mpih

  implicit none

  integer i,idx,inp

  character(70) namfile

  ! first of all, fill up velocity arrays

  if (myid.eq.0) then

    namfile='flowmov/tail_head.txt'

    open(unit=43,file=namfile,Access = 'append', Status='unknown')
    write(43,'(100E15.7)')time,dt, tail_head 
   
    close(43)
  end if

end subroutine write_tail_head



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
  integer :: nv,ne,nf
  integer :: dummyInd

    open(109,file=gtsfx)

      read(109,*) nv, ne, nf

      do i=1,maxnv
        read(109,*)xyz0(1,i),xyz0(2,i),xyz0(3,i)
      end do

      ! n_edge_of_vert(1:maxnv)=0

      ! do i=1,maxne
      !   read(109,*)v1,v2
      !   vert_of_edge(1,i)=v1
      !   vert_of_edge(2,i)=v2
      !   n_edge_of_vert(v1)=n_edge_of_vert(v1)+1
      !   n_edge_of_vert(v2)=n_edge_of_vert(v2)+1
      !   vert_of_vert(n_edge_of_vert(v1),v1)=v2           
      !   vert_of_vert(n_edge_of_vert(v2),v2)=v1           
      !   edge_of_vert(n_edge_of_vert(v1),v1)=i
      !   edge_of_vert(n_edge_of_vert(v2),v2)=i
      ! enddo

      do i=1,maxnf
          read(109,*)edge_of_face(1,i),edge_of_face(2,i), edge_of_face(3,i)
      end do

      do i=1,maxnf
        e1=edge_of_face(1,i)
        e2=edge_of_face(2,i)

        if (vert_of_edge(2,e1).eq.vert_of_edge(1,e2)) then
           v1=vert_of_edge(1,e1)
           v2=vert_of_edge(2,e1)
           v3=vert_of_edge(2,e2)
        elseif(vert_of_edge(2,e1).eq.vert_of_edge(2,e2)) then
           v1=vert_of_edge(1,e1)
           v2=vert_of_edge(2,e1)
           v3=vert_of_edge(1,e2)
        elseif(vert_of_edge(1,e1).eq.vert_of_edge(1,e2)) then
           v1=vert_of_edge(2,e1)
           v2=vert_of_edge(1,e1)
           v3=vert_of_edge(2,e2)
        else 
           v1=vert_of_edge(2,e1)
           v2=vert_of_edge(1,e1)
           v3=vert_of_edge(1,e2)
        endif 

           vert_of_face(1,i)=v1
           vert_of_face(2,i)=v2
           vert_of_face(3,i)=v3
      end do

    close(109)  !VS   Completed reading the gts file

    ! Check - for faces and edges
    face_of_edge(1:2,1:maxne)=0
    do i=1,maxnf
       e1=edge_of_face(1,i)
       e2=edge_of_face(2,i)
       e3=edge_of_face(3,i)
       if (face_of_edge(1,e1).eq.0) then
          face_of_edge(1,e1)=i
       elseif (face_of_edge(2,e1).eq.0) then
          face_of_edge(2,e1)=i
       else
          write(*,*)'Edge error', i,e1,e2,e3
       endif
       if (face_of_edge(1,e2).eq.0) then
          face_of_edge(1,e2)=i
       elseif (face_of_edge(2,e2).eq.0) then
          face_of_edge(2,e2)=i
       else
          write(*,*)'Edge error'
       endif
       if (face_of_edge(1,e3).eq.0) then
          face_of_edge(1,e3)=i
       elseif (face_of_edge(2,e3).eq.0) then
          face_of_edge(2,e3)=i
       else
          write(*,*)'Edge error'
       endif
    enddo 

    ! ! Check - vertex cannot be connected to itself
    ! do i=1,maxnv
    !   do j=1,n_edge_of_vert(i)
    !      if (vert_of_vert(j,i).eq.i)  write(*,*)'Error ',vert_of_vert(j,i),i
    !   enddo
    ! enddo

    !Check
    do i=1,maxne
       if (face_of_edge(1,i).eq.face_of_edge(2,i)) then
          write(*,*)'Error on edges '
       endif
    enddo


    ! ! KZ: Compute the face-of-vertices connectivity
    ! faces_of_vert = 0

    ! do i=1,maxnf

    !   ! Scan for vertex to assign for vertex 1
    !   do j = 1,VERTBUFFER
    !       if ( faces_of_vert( j, vert_of_face(1,i)  ) .eq. 0 ) then
    !           dummyInd = j !Next-in-line index of face to store
    !           exit
    !       endif
    !   enddo
    !   faces_of_vert(dummyInd, vert_of_face(1,i)  ) = i

    !   ! Scan for vertex to assign for vertex 2
    !   do j = 1,VERTBUFFER
    !       if ( faces_of_vert( j, vert_of_face(2,i)  ) .eq. 0 ) then
    !           dummyInd = j !Next-in-line index of face to store
    !           exit
    !       endif
    !   enddo
    !   faces_of_vert(dummyInd, vert_of_face(2,i)  ) = i


    !   ! Scan for vertex to assign for vertex 3
    !   do j = 1,VERTBUFFER
    !       if ( faces_of_vert( j, vert_of_face(3,i)  ) .eq. 0 ) then
    !           dummyInd = j !Next-in-line index of face to store
    !           exit
    !       endif
    !   enddo
    !   faces_of_vert(dummyInd, vert_of_face(3,i)  ) = i
      
    ! enddo

    ! KZ store vertex coordinates
    ! Eventually should be superceded by a COM-relative formulation
       xyzv(1:3,:,1)= xyz0(1:3,:)
       xyzv(1,:,1) = xyzv(1,:,1) + 0.5d0*xlen
       xyzv(2,:,1) = xyzv(2,:,1) + 0.5d0*ylen
       xyzv(3,:,1) = xyzv(3,:,1) + 0.5d0*zlen

end subroutine



! subroutine get_maxn_n_edge_of_vert
!   use param
!   use mpih
!   use mls_param
!   implicit none
!   real    :: dummy
!   integer :: i,j,v1,v2,v3,e1,e2,e3
!   integer :: nv,ne,nf
!   integer, dimension(:), allocatable :: num_edge_of_vert

!     open(109,file=gtsfx)

!       read(109,*) nv, ne, nf

!       do i=1,nv
!         read(109,*) dummy, dummy, dummy
!       end do

!       allocate(num_edge_of_vert(nv))
!       num_edge_of_vert(1:nv)=0

!       do i=1,maxne
!         read(109,*) v1,v2
!         num_edge_of_vert(v1) = num_edge_of_vert(v1) + 1
!         num_edge_of_vert(v2) = num_edge_of_vert(v2) + 1
!       enddo

!      max_n_edge_of_vert = -1e6
!      do i=1,nv
!        max_n_edge_of_vert = max(num_edge_of_vert(i), max_n_edge_of_vert )
!      enddo

!     close(109)
! end subroutine

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

subroutine writePPpartVol
  use param
  use mls_param
  use mpih

  IMPLICIT none

  real, dimension(Nparticle*3) :: pos
  integer                      :: i,idx,inp

  character(70) namfile


  if (myid.eq.0) then

  namfile='flowmov/partPPVol.txt'

  open(unit=43,file=namfile,Access = 'append', Status='unknown')
  write(43,'(100E15.7)')time, Volume(1) !KZ: note hard-coded single particle for now

  close(43)
  end if
end subroutine writePPpartVol



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
        call write_VertGeom (maxnv,maxnf,xyzv(:,:,inp),vmelt(:,inp),dtdn_oVert(:,inp),dtdn_iVert(:,inp),vert_nor(:,:,inp),namfi2)

      end do
  end if
  
  end subroutine write_tecplot_geom
  
  !---------------------------------------------------------------------------------------------
        subroutine write_VertGeom (nv,nf,xyz,vmelt,dtdn_o,dtdn_i,vert_nor,filename)
        use param
        use mpih
        use mls_param, only: vert_of_face
        implicit none
        character(70) filename,geotecfile
        integer i,nv,nf
        real, dimension (3,nv) :: xyz , vert_nor 
        real, dimension (nv) :: vmelt, dtdn_o, dtdn_i
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
          write(11,*)'VARIABLES = X Y Z vmelt dtdn_o dtdn_i nhat_x nhat_y nhat_z'
      !    write(11,*)'ZONE T="DOMAIN 0", N=',nv,' E=',nf,' F=FEBLOCK, ET=TRIANGLE'
          write(11,*)'ZONE T="FETri" N=',nv,' E=',nf,' ZONETYPE=FETriangle'
          ! write(11,*)'ZONE T=FETri N=',nvc,' E=',ntri,' ZONETYPE=FETriangle'
          write(11,*)'DATAPACKING=POINT                                       '
          !write(11,*)'VARLOCATION=([4-6]=CELLCENTERED)' ! KZ: lists 4-6 are centroid data, can specify as nodal for vertices
          !write(11,*)'VARLOCATION=([4-9]=NODAL)' 

          do i=1,nv
          write(11,*)xyz(1,i), xyz(2,i), xyz(3,i), vmelt(i), dtdn_o(i), dtdn_i(i), vert_nor(1,i), vert_nor(2,i), vert_nor(3,i)
          end do

          do i=1,nf
          write(11,*)vert_of_face(1:3,i)
          end do
  
  
          close(11)
          return
          end subroutine write_VertGeom
  !---------------------------------------------------------------------------------------------
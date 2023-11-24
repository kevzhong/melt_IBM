Subroutine particle
use param
use mls_param
use mls_local
use local_arrays
use mpi_param
use mpih

implicit none 

integer :: mstep, inp
integer :: my_up,my_down

if(imlsfor.eq.1)then

	call findCentroidIndices
	call mlsWeight

	fpxyz=0.0d0
	ftxyz=0.0d0

	my_down=myid-1
	my_up=myid+1

	do mstep=1,1 !KZ: iteration(?) over enforcing immersed-boundary condition
		call update_both_ghosts(n1,n2,vx,kstart,kend)
		call update_both_ghosts(n1,n2,vy,kstart,kend)
		call update_both_ghosts(n1,n2,vz,kstart,kend)
		call update_both_ghosts(n1,n2,temp,kstart,kend)


		for_xc = 0.0d0
		for_yc = 0.0d0
		for_zc = 0.0d0
		for_temp = 0.0d0

		call mlsForce

		if (imelt .eq. 1) then 
			call findProbeIndices
		! mlsMelt
		! call mlsVertWeight
		! call mlsdtdn
		endif

		call velforce
		call tempforce
	end do
endif
 
call update_both_ghosts(n1,n2,vx,kstart,kend)
call update_both_ghosts(n1,n2,vy,kstart,kend)
call update_both_ghosts(n1,n2,vz,kstart,kend)
call update_both_ghosts(n1,n2,temp,kstart,kend)



if (imlsstr.eq.1) then
	call update_part_pos
endif

! if IMELT 
	! meltVertices (Stefan condition update)
	! Update triangle centroids
	! Update triAreas
		! Update cfac -> triAreas

	! ---- after re-meshing (check) sometime above ----

	! Update faceNormals
	! Update vertexNormals

 end

! A wrapper function for tagging cells (computing the indicator functions), VOF, switching based on the tagging type specified in the input
! tagType = 0   sphereTagging, for a rigid sphere
! tagType = 1   rayTagging with signed-distance treatment at interface, general purpose for arbitrary geometry
subroutine tagCells
    use param
    use mls_param
    implicit none
    integer :: inp
    integer,dimension(3,2)     :: bbox_inds

    if (tagType .eq. 0) then
           do inp=1,Nparticle
             call get_bbox_inds(bbox_inds,inp)
             call sphereTagging(bbox_inds,inp)
           enddo
        elseif (tagType .eq. 1) then
           do inp=1,Nparticle
            if (initial_tag) then ! Ray-tagging only needed for initialization
              call get_bbox_inds(bbox_inds,inp)
              call pencilTag(bbox_inds,inp)
            endif
            call narrowBandTagging(inp)
           enddo
        else
            write(*,*) "Invalid tag type specified?"
    endif
end subroutine tagCells


! SHARED ROUTINES AMONGST ALL TAGGING TYPES

  subroutine get_bbox_inds(bbox_inds,inp)
    ! Retrieve the xyz indices of the bounding box for a given particle
    use param
    use mls_param
    implicit none
    integer :: i, nf, tri_ind
    real, dimension(3,2) :: lim
    integer, dimension(3,2) :: bbox_inds
    integer  :: inp, padSize
  
    ! Padding size indices for safety
    padSize = 2
  
  ! get bounding box
    do i = 1,3
      lim(i,1) = minval( pack(xyzv(i,:,inp) , .not. isGhostVert(:,inp)  ) )
      lim(i,2) = maxval( pack(xyzv(i,:,inp) , .not. isGhostVert(:,inp)  ) )
    end do
  
    bbox_inds = floor(lim*dx1) + 1 ! compute indices cell centered
  
    ! expanding bounding box to be extra safe
    bbox_inds(:,1) = bbox_inds(:,1) - padSize
    bbox_inds(:,2) = bbox_inds(:,2) + padSize
  
    ! Hard code vertical for testing free-fall
    !bbox_inds(3,1) = 1
    !bbox_inds(3,2) = n3m
  
  
    !! Hard-code the full domain for testing
    !bbox_inds(:,1) = [1, 1, 1] 
    !bbox_inds(:,2) = [n1m, n2m, n3m]
  
end subroutine

subroutine get_periodic_indices(k,x)
  use param
  implicit none
  integer :: k
  real    :: x(3)

  if (k .ge. n3) then
     k = k - n3m
    x(3) = x(3) - zlen
  end if

  if (k .lt. 1) then
     k = k + n3m
     x(3) = x(3) + zlen
  end if
end subroutine

subroutine interp_vof
    ! Interpolate VOFx, VOFv, VOFw from cell-centered VOFp information
    use param
    use mpi_param, only: kstart,kend
    implicit none
    integer :: ic,jc,kc
    integer :: km,kp,jm,jp,im,ip
  
    call update_both_ghosts(n1,n2,VOFp,kstart,kend)
  
    do kc=kstart,kend
      km=kc-1
      kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
            
            VOFx(ic,jc,kc) = 0.5 * (VOFp(im,jc,kc) + VOFp(ic,jc,kc) )
            VOFy(ic,jc,kc) = 0.5 * (VOFp(ic,jm,kc) + VOFp(ic,jc,kc) )
            VOFz(ic,jc,kc) = 0.5 * (VOFp(ic,jc,km) + VOFp(ic,jc,kc) )
          enddo
        enddo
    enddo
  end subroutine interp_vof
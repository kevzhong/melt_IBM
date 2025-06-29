! A routine for backing out a local color function from the genSDF algorithm of Patil et al. (2025)
! The strategy computes the signed distance only along a narrow band surrouding the triangles of the surface geometry
! A mapping from the sdf to color-function / VOF is done using the standard sin() mapping as in level-set methods
! See e.g. Olsson & Kreiss 2005 JCP
! Since this tagging strategy is only done for the narrow band, this routine is meant to be supplemented with a full
! cell-tagging at initial runtime (e.g. with ray-tagging)
! For subsequent timesteps after initialization: only neighbouring cells to the interface cells need to be updated
! That is, we should only need to re-tag cells along the narrow band at subsequent timesteps

subroutine narrowBandTagging(inp)
    use param
    use mls_param
    use mpi_param
    implicit none
    integer, intent(in) :: inp
    integer :: f, i, j, k,  kk_mod, jj_mod, ii_mod
    integer :: ii, jj, kk
    real, dimension(3) :: pos, V0, V1, V2, P, Q
    integer :: sb

    real :: sgnd, dotval, sgn_dot_nhat, vof
    real :: BIGNUMBER = 999


    sb = 2 ! The narrow band size surrouding each triangle

    ! Reset SDF array
    sdf = BIGNUMBER

    !------------------- COMPUTE SIGNED DISTANCE ALONG NARROW BAND ----------------------
    do f = 1, maxnf
        if (.not. isGhostFace(f,inp) ) then

            pos(1:3) = tri_bar(1:3,f,inp)
            i  = floor(pos(1)*dx1) + 1
            j  = floor(pos(2)*dx2) + 1
            k  = floor(pos(3)*dx3) + 1

            V0 = xyzv(1:3,vert_of_face(1,f,inp),inp)
            V1 = xyzv(1:3,vert_of_face(2,f,inp),inp)
            V2 = xyzv(1:3,vert_of_face(3,f,inp),inp)

            do kk = k-sb, k+sb
                kk_mod = modulo(kk-1,n3m)  + 1
                if ( (kk_mod .ge. kstart) .and. (kk_mod .le. kend) ) then
                    do jj = j-sb, j+sb
                        do ii = i-sb, i+sb

                            P = [ xm(ii), ym(jj), zm(kk) ]

                            call distance_point_triangle_3d(P, V0, V1, V2, sgnd)

                            ii_mod = modulo(ii-1,n1m)  + 1
                            jj_mod = modulo(jj-1,n1m)  + 1

                            !Only update if |distance| smaller than value to beoverried
                            if ( abs(sdf(ii_mod,jj_mod,kk_mod)) .gt.  sgnd ) then
                                ! Convert from unsigned to signed distance
                                Q = pos(1:3) - P
                                dotval = dot_product( tri_nor(1:3,f,inp), Q )
                                sgn_dot_nhat = sign( 1.0, dotval )
                                sgnd = sgnd * (-sgn_dot_nhat)



                                sdf(ii_mod,jj_mod,kk_mod) = sgnd

                                ! Map the SDF to a VOF
                                call sdf_to_VOF(sgnd, vof)
                                VOFp(ii_mod,jj_mod,kk_mod) = vof

                            endif
                        enddo
                    enddo
                endif
            enddo
        endif
    enddo
!--------------------- END NARROW-BAND SIGNED DISTANCE CALCULATION -----------------------------
call interp_vof

end subroutine narrowBandTagging


  subroutine distance_point_triangle_3d(P, V0, V1, V2, magDist)
    ! Algorithm here https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    ! For computing the sign-unaware (Euclidean distance) of the cell-centre P to triangle (V0, V1, V2), magDist
    implicit none
    real(8), intent(in)  :: P(3), V0(3), V1(3), V2(3)
    real(8), intent(out) :: magDist
    real(8) :: ClosestPoint(3)

    real(8) :: E0(3), E1(3), DE(3)
    real(8) :: a, b, c, d,e, f, det, s, t
    real(8) :: invDet, numer, denom
    real(8) :: tmp0, tmp1

    integer :: i

    E0 = V1 - V0
    E1 = V2 - V0
    DE  = V0 - P

    a = dot_product(E0, E0)
    b = dot_product(E0, E1)
    c = dot_product(E1, E1)
    d = dot_product(E0, DE)
    e = dot_product(E1, DE)
    f = dot_product(DE, DE)

    det = a*c - b*b
    s = b*e - c*d
    t = b*d - a*e

    if (s + t <= det) then
      if (s < 0.0d0) then
        if (t < 0.0d0) then
          ! Region 4
          if (d < 0.0d0) then
            t = 0.0d0
            if (-d >= a) then
              s = 1.0d0
            else
              s = -d / a
            end if
          else
            s = 0.0d0
            if (e >= 0.0d0) then
              t = 0.0d0
            elseif (-e >= c) then
              t = 1.0d0
            else
              t = -e / c
            end if
          end if
        else
          ! Region 3
          s = 0.0d0
          if (e >= 0.0d0) then
            t = 0.0d0
          elseif (-e >= c) then
            t = 1.0d0
          else
            t = -e / c
          end if
        end if
      elseif (t < 0.0d0) then
        ! Region 5
        t = 0.0d0
        if (d >= 0.0d0) then
          s = 0.0d0
        elseif (-d >= a) then
          s = 1.0d0
        else
          s = -d / a
        end if
      else
        ! Region 0
        invDet = 1.0d0 / det
        s = s * invDet
        t = t * invDet
      end if
    else
      if (s < 0.0d0) then
        ! Region 2
        tmp0 = b + d
        tmp1 = c + e
        if (tmp1 > tmp0) then
          numer = tmp1 - tmp0
          denom = a - 2.0d0*b + c
          if (numer >= denom) then
            s = 1.0d0
            t = 0.0d0
          else
            s = numer / denom
            t = 1.0d0 - s
          end if
        else
          s = 0.0d0
          if (tmp1 <= 0.0d0) then
            t = 1.0d0
          elseif (e >= 0.0d0) then
            t = 0.0d0
          else
            t = -e / c
          end if
        end if
      elseif (t < 0.0d0) then
        ! Region 6
        tmp0 = b + e
        tmp1 = a + d
        if (tmp1 > tmp0) then
          numer = tmp1 - tmp0
          denom = a - 2.0d0*b + c
          if (numer >= denom) then
            t = 1.0d0
            s = 0.0d0
          else
            t = numer / denom
            s = 1.0d0 - t
          end if
        else
          t = 0.0d0
          if (tmp1 <= 0.0d0) then
            s = 1.0d0
          elseif (d >= 0.0d0) then
            s = 0.0d0
          else
            s = -d / a
          end if
        end if
      else
        ! Region 1
        numer = c + e - b - d
        denom = a - 2.0d0*b + c
        if (numer <= 0.0d0) then
          s = 0.0d0
        elseif (numer >= denom) then
          s = 1.0d0
        else
          s = numer / denom
        end if
        t = 1.0d0 - s
      end if
    end if

    ! Compute closest point
    ClosestPoint = V0 + s*E0 + t*E1

    ! Compute squared distance
    DE = ClosestPoint - P
    magDist = dot_product(DE, DE)
    magDist = sqrt(magDist)
  end subroutine distance_point_triangle_3d


  subroutine sdf_to_VOF(sgnd, vof)
    use param, only: dx1,pi
    implicit none
    real, intent(inout) :: vof
    real, intent(in) :: sgnd
    real :: eps
    eps = 1.0 / dx1

    ! Map the signed distance, sgnd to a color function/VOF
    ! Common mapping employed in level-set formulations (See Olsson & Kreiss for example)
    if (sgnd .le. -eps) then
        vof = 0.0
    elseif (sgnd .ge. eps) then
        vof = 1.0
    else
        vof = 0.5 * ( 1.0 + sgnd * dx1 + (1.0/pi)*sin(pi*sgnd * dx1 )  )
    endif

  end subroutine sdf_to_VOF
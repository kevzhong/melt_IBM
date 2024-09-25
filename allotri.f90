subroutine allocate_trigeo
use param
use mls_param
use mpih
implicit none 

  allocate(vert_of_edge(2,maxne,Nparticle))
  allocate(face_of_edge(2,maxne,Nparticle))
  allocate(vert_of_face(3,maxnf,Nparticle))
  allocate(edge_of_face(3,maxnf,Nparticle))

  ! KZ: remeshing data structures
  allocate(isGhostFace(maxnf,Nparticle))
  allocate(isGhostEdge(maxne,Nparticle))
  allocate(isGhostVert(maxnv,Nparticle))
  allocate(anchorVert(maxnv,Nparticle))
  allocate(flagged_edge(maxne,Nparticle))

  allocate(eLengths(maxne,Nparticle))
  allocate(skewness(maxnf,Nparticle))
  allocate(rm_flag(Nparticle))

  isGhostFace(:,:) = .false.
  isGhostEdge(:,:) = .false.
  isGhostVert(:,:) = .false.
  rm_flag(:) = .false.
  anchorVert(:,:) = .true.
  flagged_edge(:,:) = .false.


  allocate(pind(6,maxnf,Nparticle))
  allocate(pind1(3,Nparticle))

  allocate(tri_ver(9,maxnf,Nparticle))
  allocate(vel_tri(3,maxnf,Nparticle))
  allocate(tri_bar(3,maxnf,Nparticle))
  allocate(tri_nor(3,maxnf,Nparticle))
  allocate(vert_nor(3,maxnv,Nparticle)) !KZ: vertex normal vectors

  allocate(sur(maxnf,Nparticle))
  allocate(vol(maxnf,Nparticle))
  allocate(Avert(maxnv,Nparticle))


  ! Initial geometry
  allocate(xyz0(3,maxnv))
  allocate(xyzv(3,maxnv,Nparticle))

  ! Vertex coordinates relative to COM
  allocate(dxyzv_s(3,maxnv,Nparticle))
  allocate(dxyz_s(3,maxnf,Nparticle))

  allocate(ptxAB_q1(nel,maxnf,Nparticle),ptxAB_q2(nel,maxnf,Nparticle),ptxAB_q3(nel,maxnf,Nparticle))
  allocate(  ptxAB_temp(nel,maxnf,Nparticle)  )

  allocate(  qw_o(maxnf,Nparticle) , qw_i(maxnf,Nparticle)  ) !Normal gradients at vertices
  allocate(  qw_oVert(maxnv,Nparticle) , qw_iVert(maxnv,Nparticle)  ) !Normal gradients at vertices
  allocate(vmelt(3,maxnv,Nparticle))
  allocate(vmelt_m1(3,maxnv,Nparticle))

  !-- particle
  allocate( pos_CM(3, Nparticle),    vel_CM(3, Nparticle) )
  allocate( omega_c(3, Nparticle) )
  allocate( quat(4, Nparticle)  )


    ! Structural loads at Lagrangian nodes
  allocate(  tau_n1(maxnf,Nparticle) ,tau_n2(maxnf,Nparticle), tau_n3(maxnf,Nparticle)  ) 
  allocate(  press_n_tri(3,maxnf,Nparticle)  ) 

  allocate(  r_x_tau_n1(maxnf,Nparticle) ,r_x_tau_n2(maxnf,Nparticle), r_x_tau_n3(maxnf,Nparticle)  ) 
  allocate(  r_x_prn_tri(3,maxnf,Nparticle)  ) 

  allocate(int_prn_dA(3,Nparticle))
  allocate(int_r_x_prn_dA(3,Nparticle))
  allocate(int_tau_dA(3,Nparticle))
  allocate(int_r_x_tau_dA(3,Nparticle))

  allocate(int_prn_dA_m1(3,Nparticle))
  allocate(int_r_x_prn_dA_m1(3,Nparticle))
  allocate(int_tau_dA_m1(3,Nparticle))
  allocate(int_r_x_tau_dA_m1(3,Nparticle))

  allocate(Surface(Nparticle))
  allocate(Volume(Nparticle))

  allocate(cfac(maxnf,Nparticle))
end subroutine allocate_trigeo

subroutine dealloc_trigeo
  use param
  use mls_param
  use mpih
  implicit none 

  ! LAGRANGIAN memory
  if(allocated(vert_of_edge)) deallocate(vert_of_edge)
  if(allocated(face_of_edge)) deallocate(face_of_edge)
  if(allocated(vert_of_face)) deallocate(vert_of_face)
  if(allocated(edge_of_face)) deallocate(edge_of_face)

  if(allocated(isGhostFace)) deallocate(isGhostFace)
  if(allocated(isGhostEdge)) deallocate(isGhostEdge)
  if(allocated(isGhostVert)) deallocate(isGhostVert)
  if(allocated(anchorVert)) deallocate(anchorVert)
  if(allocated(flagged_edge)) deallocate(flagged_edge)

  if(allocated(eLengths)) deallocate(eLengths)
  if(allocated(skewness)) deallocate(skewness)
  if(allocated(rm_flag)) deallocate(rm_flag)

  if(allocated(pind)) deallocate(pind)
  if(allocated(pind1)) deallocate(pind1)

  if(allocated(tri_ver)) deallocate(tri_ver)
  if(allocated(vel_tri)) deallocate(vel_tri)
  if(allocated(tri_bar)) deallocate(tri_bar)
  if(allocated(tri_nor)) deallocate(tri_nor)
  if(allocated(vert_nor)) deallocate(vert_nor)

  if(allocated(sur)) deallocate(sur)
  if(allocated(vol)) deallocate(vol)
  if(allocated(Avert)) deallocate(Avert)

  if(allocated(xyz0)) deallocate(xyz0)
  if(allocated(xyzv)) deallocate(xyzv)

  if(allocated(dxyzv_s)) deallocate(dxyzv_s)
  if(allocated(dxyz_s)) deallocate(dxyz_s)


  if(allocated(ptxAB_q1)) deallocate(ptxAB_q1)
  if(allocated(ptxAB_q2)) deallocate(ptxAB_q2)
  if(allocated(ptxAB_q3)) deallocate(ptxAB_q3)
  if(allocated(ptxAB_temp)) deallocate(ptxAB_temp)

  if(allocated(qw_o)) deallocate(qw_o)
  if(allocated(qw_i)) deallocate(qw_i)
  if(allocated(qw_oVert)) deallocate(qw_oVert)
  if(allocated(qw_iVert)) deallocate(qw_iVert)
  if(allocated(vmelt)) deallocate(vmelt)
  if(allocated(vmelt_m1)) deallocate(vmelt_m1)

  if(allocated(pos_CM)) deallocate(pos_CM)
  if(allocated(vel_CM)) deallocate(vel_CM)
  if(allocated(omega_c)) deallocate(omega_c)
  if(allocated(quat)) deallocate(quat)

  if(allocated(tau_n1)) deallocate(tau_n1)
  if(allocated(tau_n2)) deallocate(tau_n2)
  if(allocated(tau_n3)) deallocate(tau_n3)
  if(allocated(press_n_tri)) deallocate(press_n_tri)

  if(allocated(r_x_tau_n1)) deallocate(r_x_tau_n1)
  if(allocated(r_x_tau_n2)) deallocate(r_x_tau_n2)
  if(allocated(r_x_tau_n3)) deallocate(r_x_tau_n3)
  if(allocated(r_x_prn_tri)) deallocate(r_x_prn_tri)

  if(allocated(int_prn_dA)) deallocate(int_prn_dA)
  if(allocated(int_r_x_prn_dA)) deallocate(int_r_x_prn_dA)
  if(allocated(int_tau_dA)) deallocate(int_tau_dA)
  if(allocated(int_r_x_tau_dA)) deallocate(int_r_x_tau_dA)

  if(allocated(int_prn_dA_m1)) deallocate(int_prn_dA_m1)
  if(allocated(int_r_x_prn_dA_m1)) deallocate(int_r_x_prn_dA_m1)
  if(allocated(int_tau_dA_m1)) deallocate(int_tau_dA_m1)
  if(allocated(int_r_x_tau_dA_m1)) deallocate(int_r_x_tau_dA_m1)

  if(allocated(Surface)) deallocate(Surface)
  if(allocated(Volume)) deallocate(Volume)
  if(allocated(cfac)) deallocate(cfac)

end subroutine dealloc_trigeo
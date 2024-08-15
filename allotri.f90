subroutine allocate_trigeo
use param
use mls_param
use mpih
implicit none 

  !allocate(n_edge_of_vert(maxnv))
  allocate(vert_of_edge(2,maxne,Nparticle))
  allocate(face_of_edge(2,maxne,Nparticle))
  allocate(vert_of_face(3,maxnf,Nparticle))
  allocate(edge_of_face(3,maxnf,Nparticle))
  !allocate(vert_of_vert(max_n_edge_of_vert,maxnv))
  !allocate(edge_of_vert(max_n_edge_of_vert,maxnv))

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

  !allocate(faces_of_vert(VERTBUFFER,maxnv)) !KZ face-vertex connectvity for vertex normals

  allocate(pind(6,maxnf,Nparticle))
  !allocate(pindv(6,maxnv,Nparticle))
  allocate(pind1(3,Nparticle))
  !allocate(dismax(3,maxnf,Nparticle))

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

  ! Probably not needed anymore
  allocate(dxyz_CM_b(3,maxnf,Nparticle))
  allocate(dxyz_CM_s(3,maxnf,Nparticle))

  ! Vertex coordinates relative to COM
  allocate(dxyzv_s(3,maxnv,Nparticle))
  allocate(dxyz_s(3,maxnf,Nparticle))

  allocate(ptxAB_q1(nel,maxnf,Nparticle),ptxAB_q2(nel,maxnf,Nparticle),ptxAB_q3(nel,maxnf,Nparticle))
  allocate(  ptxAB_temp(nel,maxnf,Nparticle)  )

  !allocate(  ddx_ptxAB(nel,maxnv,Nparticle) , ddy_ptxAB(nel,maxnv,Nparticle), ddz_ptxAB(nel,maxnv,Nparticle) )

  allocate(  qw_o(maxnf,Nparticle) , qw_i(maxnf,Nparticle)  ) !Normal gradients at vertices
  allocate(  qw_oVert(maxnv,Nparticle) , qw_iVert(maxnv,Nparticle)  ) !Normal gradients at vertices
  allocate(vmelt(3,maxnv,Nparticle))
  allocate(vmelt_m1(3,maxnv,Nparticle))

  !-- particle
  allocate( fpxyz(3, Nparticle),     ftxyz(3, Nparticle) )
  allocate( pos_CM(3, Nparticle),    vel_CM(3, Nparticle),  a_CM(3, Nparticle) )
  allocate( omega_c(3, Nparticle) )
  allocate( u_tot(3, Nparticle),     u_tot_m1(3, Nparticle) )
  allocate( r_x_u_tot(3, Nparticle), r_x_u_tot_m1(3, Nparticle) )
  allocate( omega_s(3, Nparticle) )
  allocate( tail_head(3, Nparticle) )
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


  !allocate( quat(4, Nparticle), quat_m1(4, Nparticle), quat_dot(4, Nparticle)  )

  !allocate(  quat_dot_m1(4, Nparticle)  ) !KZ: needed?

  allocate(Surface(Nparticle))
  allocate(Volume(Nparticle))

  allocate(cfac(maxnf,Nparticle))
end 

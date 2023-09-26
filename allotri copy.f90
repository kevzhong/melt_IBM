subroutine allocate_trigeo
use param
use mls_param
use mpih
implicit none 

  allocate(n_edge_of_vert(maxnv))
  allocate(vert_of_edge(2,maxne))
  allocate(face_of_edge(2,maxne))
  allocate(vert_of_face(3,maxnf))
  allocate(edge_of_face(3,maxnf))
  allocate(vert_of_vert(max_n_edge_of_vert,maxnv))
  allocate(edge_of_vert(max_n_edge_of_vert,maxnv))

  allocate(pind(6,maxnf,Nparticle))
  allocate(pind1(3,Nparticle))
  allocate(dismax(3,maxnf,Nparticle))

  allocate(tri_ver(9,maxnf,Nparticle))
  allocate(vel_tri(3,maxnf,Nparticle))
  allocate(tri_bar(3,maxnf,Nparticle))
  allocate(tri_nor(3,maxnf,Nparticle))

  allocate(dist(maxne,Nparticle))
  allocate(sur(maxnf,Nparticle))
  allocate(vol(maxnf,Nparticle))


  allocate(xyz0(3,maxnv))
  allocate(xyzv(3,maxnv,Nparticle))

  allocate(dxyz_CM_b(3,maxnf,Nparticle))
  allocate(dxyz_CM_s(3,maxnf,Nparticle))

  allocate(ptxAB_q1(nel,maxnf,Nparticle),ptxAB_q2(nel,maxnf,Nparticle),ptxAB_q3(nel,maxnf,Nparticle))

  !-- particle
  allocate( fpxyz(3, Nparticle),     ftxyz(3, Nparticle) )
  allocate( pos_CM(3, Nparticle),    vel_CM(3, Nparticle) )
  allocate( omega_b(3, Nparticle),   omega_dot_b(3, Nparticle) )
  allocate( u_tot(3, Nparticle),     u_tot_m1(3, Nparticle) )
  allocate( r_x_u_tot(3, Nparticle), r_x_u_tot_m1(3, Nparticle) )
  allocate( omega_s(3, Nparticle) )
  allocate( om_b_sqr(3, Nparticle), om_b_sqr_m1(3, Nparticle) )
  allocate( tail_head(3, Nparticle) )
  allocate( quat(4, Nparticle), quat_m1(4, Nparticle), quat_dot(4, Nparticle) )

  allocate(Surface(Nparticle))
  allocate(Volume(Nparticle))

  allocate(cfac(maxnf))
end 

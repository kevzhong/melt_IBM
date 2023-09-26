!=====================================================================================================
        subroutine Loads (tri_ver,tri_bar,tri_sup,tri_nor,ntri,bbox,            &
                          m1,m2,m3,n1,n2,n3,x1c,x1m,x2c,x2m,x3c,x3m,            &
                          g2c,g2m,g3c,g3m,dt,dx1,dx2,dx3,const,ndex,            &
                          tri_cell,tri_cell_probe,h_probe,tri_bar_probe,        &
                          q1,q2,q3,pr,fxk,fyk,fzk,IB_force,dto,dtoo,            &
                          acc_tri,vel_tri,Pc,Uc,Vc,Wc,dUdxyz,dVdxyz,dWdxyz,     &
                          ntime,ifsi,time,ren,alx1,alx2,alx3,                   &
                          For_pres,For_visc,Mom_pres,Mom_visc,                  &
                          cen,cen_g,grav,rhom,Vol,For_grav,Mom_grav,            &
                          n_max_data,nsamp,tad,dp,rib,visct,iles,               &
                          nv,ne,nf,vert_of_face,fxyze,iopen,loadext)


        implicit none

        integer :: ndex,ntri,l,kwji,n,nst,n1m,n2m,n3m,iles
        integer :: i,j,k,m1,m2,m3,n1,n2,n3,ntime,ifsi,n_max_data,nsamp
        integer,dimension(ntri,3,4) :: tri_cell
        integer,dimension(ntri,3,4,2) :: tri_cell_probe
        integer, dimension(3,ntri,ndex) :: idp
        real, dimension (3) :: gpos,cen
        real, dimension (3,ndex) :: posk,dsk
        real, dimension (ndex) :: uk,vk,wk,pk,muk
        real, dimension (ntri,4,ndex) :: phis
        real,dimension(ntri,9) :: tri_ver
        real,dimension(ntri,3) :: tri_bar,tri_nor
        real,dimension(ntri,3,2) :: tri_bar_probe
        real,dimension(ntri) :: tri_sup
        real,dimension(3,2) :: bbox
        real,dimension(3) :: IB_force,F_Lag,F_Eul
        real, dimension(m1) :: x1c,x1m
        real, dimension(m2) :: x2c,x2m,g2c,g2m
        real, dimension(m3) :: x3c,x3m,g3c,g3m
        real,dimension(ntri) :: Uc,Vc,Wc,hll,cl,DVEl
        real,dimension(ntri,3) :: Pc
        real,dimension(ntri) :: Fxc,Fyc,Fzc,noloads,muc
        real,dimension(ntri) :: const_tri
        real,dimension(ntri,3) :: press
        real,dimension(ntri,3) :: acc_tri,vel_tri
        real, dimension(m1,m2,m3) :: q1,q2,q3,pr,fxk,fyk,fzk,rib,visct
        real :: ck01,ck02,dt,dx1,dx2,dx3,const
        real :: summa,fconst,fzlag,fzeul,dto,dtoo,time,h_probe
        real, dimension (ndex,ntri) :: phile
        real,dimension(ntri,3,3) :: dUdxyz,dVdxyz,dWdxyz
        real :: rhom,Vol,e11,e12,e13,e21,e22,e23,e31,e32,e33,visc,ren,dpdn
        real, dimension(3) :: For_pres,For_visc,Mom_pres,Mom_visc
        real, dimension(3) :: tau,rho,cen_g,grav,For_grav,Mom_grav
        real :: alx1,alx2,alx3,control,overpressure
        real, dimension (n_max_data) :: tad,dp
        integer :: nv,ne,nf,iopen,v1,v2,v3
        real, dimension (3,nv) :: fxyze
        integer, dimension (3,nf) :: vert_of_face
        real, dimension (2,nf) :: loadext
    

        n1m=n1-1
        n2m=n2-1
        n3m=n3-1

        const_tri(1:ntri)=const
!       rib(1:n1m,1:n2m,1:n3) = 0.0

!====== U positive normal direction =========================================

        do l=1,ntri
           gpos(1:3)=tri_bar_probe(l,1:3,1)

           i=tri_cell_probe(l,1,1,1)
           j=tri_cell_probe(l,2,1,1)
           k=tri_cell_probe(l,3,1,1)

           dUdxyz(l,1:3,1)=0.0

           if (i+j+k.ne.0) then

               call get_indices_27(1,0,1,1,ndex,n1,n2m,n3m,i,j,k,idp(1:3,l,1:ndex),       &
                                  posk(1:3,1:ndex),dsk(1:3,1:ndex),                       &
                                  x1c,x2m,x3m,m1,m2,m3,dx1,dx2,dx3,g2m,g3m,const_tri(l),  &
                                  alx1,alx2,alx3)

               do n=1,ndex
                  uk(n)=q1(idp(1,l,n),idp(2,l,n),idp(3,l,n))
               enddo

               !---- Find Phi and its derivatives ----------------------------------

               call MLS_ShapeFunc_3D(const_tri(l),gpos,posk,dsk(1:3,1:ndex),phis(l,1:4,1:ndex),ndex,kwji,1)

               phile(1:ndex,l)=phis(l,1,1:ndex)

               do n=1,ndex
                  dUdxyz(l,1,1)=dUdxyz(l,1,1)+phis(l,2,n)*uk(n)
                  dUdxyz(l,2,1)=dUdxyz(l,2,1)+phis(l,3,n)*uk(n)
                  dUdxyz(l,3,1)=dUdxyz(l,3,1)+phis(l,4,n)*uk(n)
               enddo

           endif

        enddo


!====== V positive normal direction =========================================

        do l=1,ntri

           gpos(1:3)=tri_bar_probe(l,1:3,1)

           i=tri_cell_probe(l,1,2,1)
           j=tri_cell_probe(l,2,2,1)
           k=tri_cell_probe(l,3,2,1)

           dVdxyz(l,1:3,1)=0.0

           if (i+j+k.ne.0) then

               call get_indices_27(2,1,0,1,ndex,n1m,n2,n3m,i,j,k,idp(1:3,l,1:ndex),        &
                                  posk(1:3,1:ndex),dsk(1:3,1:ndex),                       &
                                  x1m,x2c,x3m,m1,m2,m3,dx1,dx2,dx3,g2c,g3m,const_tri(l),  &
                                  alx1,alx2,alx3)

               do n=1,ndex
                  vk(n)=q2(idp(1,l,n),idp(2,l,n),idp(3,l,n))
               enddo

               !---- Find Phi and its derivatives ----------------------------------

               call MLS_ShapeFunc_3D(const_tri(l),gpos,posk,dsk(1:3,1:ndex),phis(l,1:4,1:ndex),ndex,kwji,2)

               phile(1:ndex,l)=phis(l,1,1:ndex)

               do n=1,ndex
                  dVdxyz(l,1,1)=dVdxyz(l,1,1)+phis(l,2,n)*vk(n)
                  dVdxyz(l,2,1)=dVdxyz(l,2,1)+phis(l,3,n)*vk(n)
                  dVdxyz(l,3,1)=dVdxyz(l,3,1)+phis(l,4,n)*vk(n)
               enddo

           endif

        enddo


!====== W positive normal direction =========================================

        do l=1,ntri

           gpos(1:3)=tri_bar_probe(l,1:3,1)

           i=tri_cell_probe(l,1,3,1)
           j=tri_cell_probe(l,2,3,1)
           k=tri_cell_probe(l,3,3,1)

           dWdxyz(l,1:3,1)=0.0

           if (i+j+k.ne.0) then

               call get_indices_27(3,1,1,0,ndex,n1m,n2m,n3,i,j,k,idp(1:3,l,1:ndex),        &
                                  posk(1:3,1:ndex),dsk(1:3,1:ndex),                       &
                                  x1m,x2m,x3c,m1,m2,m3,dx1,dx2,dx3,g2m,g3c,const_tri(l),  &
                                  alx1,alx2,alx3)

               do n=1,ndex
                  wk(n)=q3(idp(1,l,n),idp(2,l,n),idp(3,l,n))
               enddo

               !---- Find Phi and its derivatives ----------------------------------

               call MLS_ShapeFunc_3D(const_tri(l),gpos,posk,dsk(1:3,1:ndex),phis(l,1:4,1:ndex),ndex,kwji,3)

               phile(1:ndex,l)=phis(l,1,1:ndex)

               do n=1,ndex
                  dWdxyz(l,1,1)=dWdxyz(l,1,1)+phis(l,2,n)*wk(n)
                  dWdxyz(l,2,1)=dWdxyz(l,2,1)+phis(l,3,n)*wk(n)
                  dWdxyz(l,3,1)=dWdxyz(l,3,1)+phis(l,4,n)*wk(n)
               enddo

           endif

        enddo

!====== Viscosity at cell center ============================================

        do l=1,ntri

           gpos(1:3)=tri_bar(l,1:3)

           i=tri_cell(l,1,4)
           j=tri_cell(l,2,4)
           k=tri_cell(l,3,4)

!           Pc(l,2)=0.0

           if (iles.ne.0) muc(l)=0.0

           if (i+j+k.ne.0) then

               call get_indices_27(4,1,1,1,ndex,n1m,n2m,n3m,i,j,k,idp(1:3,l,1:ndex),       &
                                  posk(1:3,1:ndex),dsk(1:3,1:ndex),                       &
                                  x1m,x2m,x3m,m1,m2,m3,dx1,dx2,dx3,g2m,g3m,const_tri(l),  &
                                  alx1,alx2,alx3)

!               do n=1,ndex
!                  pk(n)=pr(idp(1,l,n),idp(2,l,n),idp(3,l,n))
!               enddo

               if (iles.ne.0) then
                  do n=1,ndex
                     muk(n)=visct(idp(1,l,n),idp(2,l,n),idp(3,l,n))
                  enddo
               endif

               !---- Find Phi and its derivatives ----------------------------------

               call MLS_ShapeFunc_3D(const_tri(l),gpos,posk,dsk(1:3,1:ndex),phis(l,1:4,1:ndex),ndex,kwji,4)

!               do n=1,ndex
!                  Pc(l,2)=Pc(l,2)+phis(l,1,n)*pk(n)
!               enddo

               if (iles.eq.3) then
                  do n=1,ndex
                     muc(l)=muc(l)+phis(l,1,n)*muk(n)
                  enddo
               endif

           endif

        enddo

!====== P positive normal direction =========================================

        do l=1,ntri

           gpos(1:3)=tri_bar_probe(l,1:3,1)

           i=tri_cell_probe(l,1,4,1)
           j=tri_cell_probe(l,2,4,1)
           k=tri_cell_probe(l,3,4,1)

           Pc(l,1)=0.0

           if (i+j+k.ne.0) then

               call get_indices_27(4,1,1,1,ndex,n1m,n2m,n3m,i,j,k,idp(1:3,l,1:ndex),       &
                                  posk(1:3,1:ndex),dsk(1:3,1:ndex),                       &
                                  x1m,x2m,x3m,m1,m2,m3,dx1,dx2,dx3,g2m,g3m,const_tri(l),  &
                                  alx1,alx2,alx3)

               do n=1,ndex
!                 rib(idp(1,l,n),idp(2,l,n),idp(3,l,n)) = 1.0
                  pk(n)=pr(idp(1,l,n),idp(2,l,n),idp(3,l,n))
               enddo

               !---- Find Phi and its derivatives ----------------------------------

               call MLS_ShapeFunc_3D(const_tri(l),gpos,posk,dsk(1:3,1:ndex),phis(l,1:4,1:ndex),ndex,kwji,4)


               do n=1,ndex
                  Pc(l,1)=Pc(l,1)+phis(l,1,n)*pk(n)
               enddo

           endif

        enddo


301     format(50(1x,1pe18.10))

!------ Calculate Forces --------------------------------------------------------

        do l=1,ntri
           ! positive probe
           dpdn = -1.0*(acc_tri(l,1)*tri_nor(l,1)+acc_tri(l,2)*tri_nor(l,2)+acc_tri(l,3)*tri_nor(l,3))!0.0
           press(l,1)=Pc(l,1)-dpdn*h_probe

           loadext(1,l) = Pc(l,1)
        enddo

        For_pres(1:3)=0.0 
        For_visc(1:3)=0.0 
        Mom_pres(1:3)=0.0 
        Mom_visc(1:3)=0.0 
        fxyze(1:3,1:nv)=0.0

        visc=1.0/ren

        do l=1,ntri
           if (iles.ne.0) visc=visc+muc(l)

           rho(1:3)=tri_bar(l,1:3)-cen(1:3)

           press(l,3)=press(l,1)
           For_pres(1:3)=For_pres(1:3)-press(l,3)*tri_nor(l,1:3)*tri_sup(l)

           v1=vert_of_face(1,l)
           v2=vert_of_face(2,l)
           v3=vert_of_face(3,l)

           fxyze(1:3,v1)=fxyze(1:3,v1)-press(l,3)*tri_nor(l,1:3)*tri_sup(l)/3.
           fxyze(1:3,v2)=fxyze(1:3,v2)-press(l,3)*tri_nor(l,1:3)*tri_sup(l)/3.
           fxyze(1:3,v3)=fxyze(1:3,v3)-press(l,3)*tri_nor(l,1:3)*tri_sup(l)/3.

           e11=dUdxyz(l,1,1)
           e12=0.5*(dUdxyz(l,2,1)+dVdxyz(l,1,1))
           e13=0.5*(dUdxyz(l,3,1)+dWdxyz(l,1,1))
           e21=e12
           e22=dVdxyz(l,2,1)
           e23=0.5*(dVdxyz(l,3,1)+dWdxyz(l,2,1))
           e31=e13
           e32=e23
           e33=dWdxyz(l,3,1)

           tau(1)=2.0*visc*(e11*tri_nor(l,1)+e12*tri_nor(l,2)+e13*tri_nor(l,3))
           tau(2)=2.0*visc*(e21*tri_nor(l,1)+e22*tri_nor(l,2)+e23*tri_nor(l,3))
           tau(3)=2.0*visc*(e31*tri_nor(l,1)+e32*tri_nor(l,2)+e33*tri_nor(l,3))

           fxyze(1:3,v1)=fxyze(1:3,v1)+tau(1:3)*tri_sup(l)/3.
           fxyze(1:3,v2)=fxyze(1:3,v2)+tau(1:3)*tri_sup(l)/3.
           fxyze(1:3,v3)=fxyze(1:3,v3)+tau(1:3)*tri_sup(l)/3.

           For_visc(1)=For_visc(1)+tri_sup(l)*tau(1)
           For_visc(2)=For_visc(2)+tri_sup(l)*tau(2)
           For_visc(3)=For_visc(3)+tri_sup(l)*tau(3)

!          loadext(1,l) = sqrt( (tau(1)-press(l,3)*tri_nor(l,1))**2+ &
!                             (tau(2)-press(l,3)*tri_nor(l,2))**2+ &
!                             (tau(3)-press(l,3)*tri_nor(l,3))**2 )
!          loadext(2,l) = press(l,3)
!          loadext(1,l) = press(l,1)
!          loadext(2,l) = press(l,2)

           Mom_pres(1)=Mom_pres(1)-press(l,1)*tri_sup(l)*(tri_nor(l,3)*rho(2)-tri_nor(l,2)*rho(3))
           Mom_pres(2)=Mom_pres(2)-press(l,1)*tri_sup(l)*(tri_nor(l,1)*rho(3)-tri_nor(l,3)*rho(1))
           Mom_pres(3)=Mom_pres(3)-press(l,1)*tri_sup(l)*(tri_nor(l,2)*rho(1)-tri_nor(l,1)*rho(2))
           Mom_visc(1)=Mom_visc(1)+tri_sup(l)*(tau(3)*rho(2)-tau(2)*rho(3))
           Mom_visc(2)=Mom_visc(2)+tri_sup(l)*(tau(1)*rho(3)-tau(3)*rho(1))
           Mom_visc(3)=Mom_visc(3)+tri_sup(l)*(tau(2)*rho(1)-tau(1)*rho(2))

        enddo

        !gravity
        Mom_grav(1:3)=0.0


        return
        end subroutine Loads
!.....................................................................................................
!=====================================================================================================

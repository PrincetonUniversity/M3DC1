module physical_mesh
  use element
  use read_vmec

  implicit none

  real :: xcenter, zcenter ! coords of the center of the logical mesh
  integer :: igeometry   ! 0 = identity; 1 = prescribed; 2 = solve Laplace equation  
  integer :: iread_vmec 
  integer :: nperiods
#ifdef USEST
  ! 0 = physical basis;, 1 = logical basis & DoF; 2= logical basis, physical DoF. 
  integer :: ilog        
  real :: mesh_period 
  ! arrays to store dofs of rst & zst fields 
  real, allocatable :: rstnode(:,:)
  real, allocatable :: zstnode(:,:)
  real :: m_max, n_max, dphi, mf
  real :: rm0, rm1, rm2, zm0, zm1

contains

  subroutine physical_mesh_setup(period)
    use math
    implicit none

    real, intent(in) :: period 

    mesh_period = period/nperiods 
    ! metric factor for phi when itor=0
    mf = twopi/(mesh_period*nperiods)
    ! phase shift
    dphi = 0.*twopi/(2*nperiods)
    !dphi = 1.*mesh_period/2
    ! maxium m & n numbers
    m_max = 10.5
    n_max = 18.5
    ! 2D analytical boundary parameters
    rm0 = 6.
    rm1 = 2.5
    rm2 = 0.
    zm0 = 0.
    zm1 = 2.
  
  end subroutine physical_mesh_setup


  ! Return rst and zst given x, phi, z 
  elemental subroutine physical_geometry(rout, zout, x, phi, z)
    use math
    implicit none

    real, intent(in) :: x, phi, z
    real, intent(out) :: rout, zout 
    real :: r, theta, ds, r2n
    integer :: i, js
    real, dimension(mn_mode) :: rstc, zsts, co, sn 
    real :: phis

    phis = phi*mf+dphi

    r = sqrt((x - xcenter)**2 + (z - zcenter)**2 + 0e-6)
    theta = atan2(z - zcenter, x - xcenter)
    rout = 0
    zout = 0
    if (iread_vmec==1) then ! use VMEC surfaces to calculate geometry
      co = cos(xmv*theta+xnv*phis)
      sn = sin(xmv*theta+xnv*phis)
      r2n = r**2*(ns-1)
      js = ceiling(r2n)
      if (js>(ns-1)) js = ns-1 
!      if (js>1) then ! interpolate between surfaces
      if (js>ns) then ! do not interpolate 
        ds = js - r2n 
        rstc = rmnc(:,js+1)*(1-ds) + rmnc(:,js)*ds
        zsts = zmns(:,js+1)*(1-ds) + zmns(:,js)*ds
        do i = 1, mn_mode 
          if (xmv(i)<m_max .and. abs(xnv(i))<n_max) then
            rout = rout + rstc(i)*co(i)
            zout = zout + zsts(i)*sn(i)
          end if
        end do
      else ! near axis, use routine 
!        ds = sqrt(r2n) 
!        rstc = rmnc(:,js+1)
!        zsts = zmns(:,js+1)
        ds = r 
        rstc = rbc
        zsts = zbs
        do i = 1, mn_mode 
          if (xmv(i)<m_max .and. abs(xnv(i))<n_max) then
            rout = rout + rstc(i)*co(i)*ds**xmv(i)
            zout = zout + zsts(i)*sn(i)*ds**xmv(i)
          end if
        end do
        do i = 1, n_tor ! shift axis to magnetic axis 
            if (xmv(i)<m_max .and. abs(xnv(i))<n_max) then
              rout = rout + (rmnc(i,1)-rbc(i))*co(i)*(1-ds**2)
              zout = zout + (zmns(i,1)-zbs(i))*sn(i)*(1-ds**2)
            end if
        end do
      end if
    else ! use routine to generate from boundary geometry
      co = cos(mb*theta-nb*phis*nperiods)
      sn = sin(mb*theta-nb*phis*nperiods)
      do i = 1, mn_mode 
        if ((mb(i)<m_max)) then
          rout = rout + rbc(i)*co(i)*r**mb(i)
          zout = zout + zbs(i)*sn(i)*r**mb(i)
        end if
      end do
    end if
!    rout = rm0 + rm1*r*cos(theta+rm2*sin(theta)) 
!    zout = zm0 + zm1*r*sin(theta) 
  end subroutine physical_geometry

  ! Calculate curvature and normal vector on physical boundary. 
  ! Works only for circular logical meshes.
  subroutine get_boundary_curv(normal, curv, inode, curv3)
    use math
    implicit none

    integer, intent(in) :: inode 
    real, intent(out) :: normal(2), curv
    real, intent(out), optional :: curv3(2)
    real :: dr, ddr, dz, ddz, theta, x, phi, z, phis
    real, dimension(mn_mode) :: co, sn 
#ifdef USE3D
    real :: dr1, ddr1, dz1, ddz1 
#endif
    real :: coords(3)
    integer :: i
    real, dimension(dofs_per_node) :: rn, zn 

    ! get logical coordinates    
    call m3dc1_node_getcoord(inode-1,coords)
    x = coords(1)
    z = coords(2)
    phi = coords(3)
    phis = phi*mf + dphi

    theta = atan2(z - zcenter, x - xcenter)
    dr = 0
    dz = 0
    ddr = 0
    ddz = 0
#ifdef USE3D
    dr1 = 0
    dz1 = 0
    ddr1 = 0
    ddz1 = 0
#endif

    if (igeometry.eq.1 .and. ilog.ne.1) then ! use rst, zst node data 
       rn = rstnode(:,inode) 
       zn = zstnode(:,inode) 
       dr = rn(DOF_DR)*(-z) + rn(DOF_DZ)*x
       dz = zn(DOF_DR)*(-z) + zn(DOF_DZ)*x
       ddr = rn(DOF_DRR)*(-z)**2 + rn(DOF_DZZ)*x**2 &
            +2*rn(DOF_DRZ)*(-x*z) &
            +rn(DOF_DR)*(-x) + rn(DOF_DZ)*(-z)
       ddz = zn(DOF_DRR)*(-z)**2 + zn(DOF_DZZ)*x**2 &
            +2*zn(DOF_DRZ)*(-x*z) &
            +zn(DOF_DR)*(-x) + zn(DOF_DZ)*(-z)
#ifdef USE3D
       dr1 = rn(DOF_DRP)*(-z) + rn(DOF_DZP)*x
       dz1 = zn(DOF_DRP)*(-z) + zn(DOF_DZP)*x
       ddr1 = rn(DOF_DRRP)*(-z)**2 + rn(DOF_DZZP)*x**2 &
             +2*rn(DOF_DRZP)*(-x*z) &
             +rn(DOF_DRP)*(-x) + rn(DOF_DZP)*(-z)
       ddz1 = zn(DOF_DRRP)*(-z)**2 + zn(DOF_DZZP)*x**2 &
             +2*zn(DOF_DRZP)*(-x*z) &
             +zn(DOF_DRP)*(-x) + zn(DOF_DZP)*(-z)
#endif
       !if (inode.eq.1) print *, 'this is a test' 
    else ! use Fourier coefficients
       if (iread_vmec==1) then
         co = cos(xmv*theta+xnv*phis)
         sn = sin(xmv*theta+xnv*phis)
         do i = 1, mn_mode 
           if (xmv(i)<m_max .and. abs(xnv(i))<n_max) then
             dr = dr - rbc(i)*sn(i)*xmv(i)
             dz = dz + zbs(i)*co(i)*xmv(i)
             ddr = ddr - rbc(i)*co(i)*xmv(i)**2
             ddz = ddz - zbs(i)*sn(i)*xmv(i)**2
#ifdef USE3D
             dr1 = dr1 - rbc(i)*co(i)*xmv(i)*xnv(i)*mf
             dz1 = dz1 - zbs(i)*sn(i)*xmv(i)*xnv(i)*mf
             ddr1 = ddr1 + rbc(i)*sn(i)*xmv(i)**2*xnv(i)*mf
             ddz1 = ddz1 - zbs(i)*co(i)*xmv(i)**2*xnv(i)*mf
#endif
           end if
         end do
       else
         co = cos(mb*theta-nb*phis*nperiods)
         sn = sin(mb*theta-nb*phis*nperiods)
         do i = 1, mn_mode 
           if ((mb(i)<m_max)) then
             dr = dr - rbc(i)*sn(i)*mb(i)
             dz = dz + zbs(i)*co(i)*mb(i)
             ddr = ddr - rbc(i)*co(i)*mb(i)**2
             ddz = ddz - zbs(i)*sn(i)*mb(i)**2
#ifdef USE3D
             dr1 = dr1 - rbc(i)*co(i)*mb(i)*(-nb(i))*nperiods*mf
             dz1 = dz1 - zbs(i)*sn(i)*mb(i)*(-nb(i))*nperiods*mf
             ddr1 = ddr1 + rbc(i)*sn(i)*mb(i)**2*(-nb(i))*nperiods*mf
             ddz1 = ddz1 - zbs(i)*co(i)*mb(i)**2*(-nb(i))*nperiods*mf
#endif
           end if 
         end do
       end if 
    end if 
!    dr = -rm1*sin(theta+rm2*sin(theta))*(1+rm2*cos(theta)) 
!    ddr = -rm1*cos(theta+rm2*sin(theta))*(1+rm2*cos(theta))**2 &
!          +rm1*sin(theta+rm2*sin(theta))*rm2*sin(theta)
!    dz = zm1*cos(theta) 
!    ddz = -zm1*sin(theta) 
#ifdef USE3D
!    dr1 = 0.
!    dz1 = 0.
!    ddr1 = 0.
!    ddz1 = 0.
#endif
    curv = (dr*ddz - dz*ddr)/((dr**2 + dz**2)*sqrt(dr**2 + dz**2))
    normal(1) = dz/sqrt(dr**2 + dz**2)
    normal(2) = -dr/sqrt(dr**2 + dz**2)
#ifdef USE3D
    if(present(curv3)) then
       curv3(1) = (dr*dz1 - dz*dr1)/(dr**2 + dz**2)
       curv3(2) = ((dr1*ddz + dr*ddz1 - dz1*ddr - dz*ddr1)/(sqrt(dr**2 + dz**2))&
                  - 3.*curv*(dr*dr1 + dz*dz1))/(dr**2 + dz**2) 
    end if
#endif
  end subroutine get_boundary_curv

  ! calculate matrix that transforms logical dofs to physical
  pure subroutine l2p_matrix(l2p,inode)
    implicit none

    integer, intent(in) :: inode 
    real, intent(out) :: l2p(dofs_per_node, dofs_per_node) 
    real :: di, di2, di3, dx, dy, f, g
    real, dimension(dofs_per_node) :: rn, zn 

    rn = rstnode(:,inode) 
    zn = zstnode(:,inode) 

    ! calculate expressions needed
    ! inverse of D = Rx*Zy - Ry*Zx
    di = 1./(rn(DOF_DR)*zn(DOF_DZ) - zn(DOF_DR)*rn(DOF_DZ))
    di2 = di*di
    di3 = di*di2
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, di(1) 
    ! Dx = Rx*Zxy + Rxx*Zy - Ry*Zxx - Rxy*Zx  
    dx = rn(DOF_DR)*zn(DOF_DRZ) + rn(DOF_DRR)*zn(DOF_DZ)&
       - rn(DOF_DZ)*zn(DOF_DRR) - rn(DOF_DRZ)*zn(DOF_DR)
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, dx(1) 
    ! Dy = Rx*Zyy + Rxy*Zy - Ry*Zxy - Ryy*Zx 
    dy = rn(DOF_DR)*zn(DOF_DZZ) + rn(DOF_DRZ)*zn(DOF_DZ)&
       - rn(DOF_DZ)*zn(DOF_DRZ) - rn(DOF_DZZ)*zn(DOF_DR)
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, dy(1) 
    ! F = Rx*Dy - Ry*Dx
    f = rn(DOF_DR)*dy - rn(DOF_DZ)*dx
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, f(1) 
    ! G = Zx*Dy - Zy*Dx
    g = zn(DOF_DR)*dy - zn(DOF_DZ)*dx
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, g(1) 

    l2p = 0
    l2p(DOF_1,DOF_1) = 1.
    ! fR = (Zy/D)*fx - (Zx/D)*fy
    l2p(DOF_DR,DOF_DR) = di*zn(DOF_DZ)
    l2p(DOF_DR,DOF_DZ) = -di*zn(DOF_DR)
    ! fZ = (Rx/D)*fy - (Ry/D)*fx
    l2p(DOF_DZ,DOF_DZ) = di*rn(DOF_DR)
    l2p(DOF_DZ,DOF_DR) = -di*rn(DOF_DZ)
    ! fRR = (Zy/D)^2*fxx + (Zx/D)^2*fyy - 2(Zx*Zy/D^2)*fxy
    !     + [(Zy*Zxy - Zx*Zyy)/D^2 + G*Zy/D^3]*fx      
    !     + [(Zx*Zxy - Zy*Zxx)/D^2 - G*Zx/D^3]*fy      
    l2p(DOF_DRR,DOF_DRR) = di2*zn(DOF_DZ)**2
    l2p(DOF_DRR,DOF_DZZ) = di2*zn(DOF_DR)**2
    l2p(DOF_DRR,DOF_DRZ) = -2*di2*zn(DOF_DR)*zn(DOF_DZ)
    l2p(DOF_DRR,DOF_DR) = (zn(DOF_DZ)*zn(DOF_DRZ) - zn(DOF_DR)*zn(DOF_DZZ))*di2&
                          + g*zn(DOF_DZ)*di3
    l2p(DOF_DRR,DOF_DZ) = (zn(DOF_DR)*zn(DOF_DRZ) - zn(DOF_DZ)*zn(DOF_DRR))*di2&
                          - g*zn(DOF_DR)*di3
    ! fZZ = (Ry/D)^2*fxx + (Rx/D)^2*fyy - 2(Rx*Ry/D^2)*fxy
    !     + [(Ry*Rxy - Rx*Ryy)/D^2 + F*Ry/D^3]*fx      
    !     + [(Rx*Rxy - Ry*Rxx)/D^2 - F*Rx/D^3]*fy      
    l2p(DOF_DZZ,DOF_DRR) = di2*rn(DOF_DZ)**2
    l2p(DOF_DZZ,DOF_DZZ) = di2*rn(DOF_DR)**2
    l2p(DOF_DZZ,DOF_DRZ) = -2*di2*rn(DOF_DR)*rn(DOF_DZ)
    l2p(DOF_DZZ,DOF_DR) = (rn(DOF_DZ)*rn(DOF_DRZ) - rn(DOF_DR)*rn(DOF_DZZ))*di2&
                          + f*rn(DOF_DZ)*di3
    l2p(DOF_DZZ,DOF_DZ) = (rn(DOF_DR)*rn(DOF_DRZ) - rn(DOF_DZ)*rn(DOF_DRR))*di2&
                          - f*rn(DOF_DR)*di3
    ! fRZ = [(Rx*Zy + Ry*Zx)/D^2]*fxy - (Ry*Zy/D^2)*fxx - (Rx*Zx/D^2)*fyy 
    !     - [(Zy*Rxy - Zx*Ryy)/D^2 + G*Ry/D^3]*fx      
    !     - [(Zx*Rxy - Zy*Rxx)/D^2 - G*Rx/D^3]*fy      
    l2p(DOF_DRZ,DOF_DRZ) = di2*(rn(DOF_DR)*zn(DOF_DZ) + zn(DOF_DR)*rn(DOF_DZ)) 
    l2p(DOF_DRZ,DOF_DRR) = -di2*rn(DOF_DZ)*zn(DOF_DZ)
    l2p(DOF_DRZ,DOF_DZZ) = -di2*rn(DOF_DR)*zn(DOF_DR) 
    l2p(DOF_DRZ,DOF_DR) = -((zn(DOF_DZ)*rn(DOF_DRZ) - zn(DOF_DR)*rn(DOF_DZZ))*di2&
                          + g*rn(DOF_DZ)*di3)
    l2p(DOF_DRZ,DOF_DZ) = -((zn(DOF_DR)*rn(DOF_DRZ) - zn(DOF_DZ)*rn(DOF_DRR))*di2&
                          - g*rn(DOF_DR)*di3)

  end subroutine l2p_matrix

  ! calculate matrix that transforms physical dofs to logical
  pure subroutine p2l_matrix(p2l,inode)
    implicit none

    integer, intent(in) :: inode 
    real, intent(out) :: p2l(dofs_per_node, dofs_per_node) 
    real, dimension(dofs_per_node) :: rn, zn 

    rn = rstnode(:,inode) 
    zn = zstnode(:,inode) 

!    rn(DOF_DP:DOF_DZZP) = 0.
!    zn(DOF_DP:DOF_DZZP) = 0.

    p2l = 0
    p2l(DOF_1,DOF_1) = 1.
    ! fx = Rx*fR + Zx*fZ
    p2l(DOF_DR,DOF_DR) = rn(DOF_DR)
    p2l(DOF_DR,DOF_DZ) = zn(DOF_DR)
    ! fy = Ry*fR + Zy*fZ
    p2l(DOF_DZ,DOF_DR) = rn(DOF_DZ)
    p2l(DOF_DZ,DOF_DZ) = zn(DOF_DZ)
    ! fxx = Rx^2*fRR + Zx^2*fZZ + 2(Rx*Zx)*fRZ
    !     + Rxx*fR + Zxx*fZ 
    p2l(DOF_DRR,DOF_DRR) = rn(DOF_DR)**2
    p2l(DOF_DRR,DOF_DZZ) = zn(DOF_DR)**2
    p2l(DOF_DRR,DOF_DRZ) = 2*zn(DOF_DR)*rn(DOF_DR)
    p2l(DOF_DRR,DOF_DR) = rn(DOF_DRR)
    p2l(DOF_DRR,DOF_DZ) = zn(DOF_DRR)
    ! fyy = Ry^2*fRR + Zy^2*fZZ + 2(Ry*Zy)*fRZ
    !     + Ryy*fR + Zyy*fZ
    p2l(DOF_DZZ,DOF_DRR) = rn(DOF_DZ)**2
    p2l(DOF_DZZ,DOF_DZZ) = zn(DOF_DZ)**2
    p2l(DOF_DZZ,DOF_DRZ) = 2*zn(DOF_DZ)*rn(DOF_DZ)
    p2l(DOF_DZZ,DOF_DR) = rn(DOF_DZZ)
    p2l(DOF_DZZ,DOF_DZ) = zn(DOF_DZZ)
    ! fxy = Ry*Rx*fRR + Zy*Zx*fZZ + (Rx*Zy+Ry*Zx)*fRZ
    !     + Rxy*fR + Zxy*fZ 
    p2l(DOF_DRZ,DOF_DRR) = rn(DOF_DR)*rn(DOF_DZ)
    p2l(DOF_DRZ,DOF_DZZ) = zn(DOF_DR)*zn(DOF_DZ) 
    p2l(DOF_DRZ,DOF_DRZ) = rn(DOF_DR)*zn(DOF_DZ) + zn(DOF_DR)*rn(DOF_DZ) 
    p2l(DOF_DRZ,DOF_DR) = rn(DOF_DRZ) 
    p2l(DOF_DRZ,DOF_DZ) = zn(DOF_DRZ)
#ifdef USE3D
    p2l(DOF_DP,DOF_DP) = 1.
    ! fxz = Rx*fRz + Zx*fZz + Rxz* fR + Zxz * fZ
    p2l(DOF_DRP,DOF_DRP) = rn(DOF_DR)
    p2l(DOF_DRP,DOF_DZP) = zn(DOF_DR)
    p2l(DOF_DRP,DOF_DR) = rn(DOF_DRP)
    p2l(DOF_DRP,DOF_DZ) = zn(DOF_DRP)
    ! fyz = Ry*fRz + Zy*fZz + Ryz* fR + Zyz * fZ
    p2l(DOF_DZP,DOF_DRP) = rn(DOF_DZ)
    p2l(DOF_DZP,DOF_DZP) = zn(DOF_DZ)
    p2l(DOF_DZP,DOF_DR) = rn(DOF_DZP)
    p2l(DOF_DZP,DOF_DZ) = zn(DOF_DZP)
    ! fxxz = Rx^2*fRRz + Zx^2*fZZz + 2(Rx*Zx)*fRZz
    !     + Rxx*fRz + Zxx*fZz + Rxxz*fR + Zxxz*fZ 
    !     + 2*Rx*Rxz*fRR + 2*Zx*Zxz*fZZ + 2(Rxz*Zx+Rx*Zxz)*fRZ
    p2l(DOF_DRRP,DOF_DRRP) = rn(DOF_DR)**2
    p2l(DOF_DRRP,DOF_DZZP) = zn(DOF_DR)**2
    p2l(DOF_DRRP,DOF_DRZP) = 2*zn(DOF_DR)*rn(DOF_DR)
    p2l(DOF_DRRP,DOF_DRP) = rn(DOF_DRR)
    p2l(DOF_DRRP,DOF_DZP) = zn(DOF_DRR)
    p2l(DOF_DRRP,DOF_DRR) = 2*rn(DOF_DR)*rn(DOF_DRP)
    p2l(DOF_DRRP,DOF_DZZ) = 2*zn(DOF_DR)*zn(DOF_DRP)
    p2l(DOF_DRRP,DOF_DRZ) = 2*(zn(DOF_DRP)*rn(DOF_DR)+zn(DOF_DR)*rn(DOF_DRP))
    p2l(DOF_DRRP,DOF_DR) = rn(DOF_DRRP)
    p2l(DOF_DRRP,DOF_DZ) = zn(DOF_DRRP)
    ! fyyz = Ry^2*fRRz + Zy^2*fZZz + 2(Ry*Zy)*fRZz
    !     + Ryy*fRz + Zyy*fZz + Ryyz*fR + Zyyz*fZ 
    !     + 2*Ry*Ryz*fRR + 2*Zy*Zyz*fZZ + 2(Ry*Zyz+Ryz*Zy)*fRZ
    p2l(DOF_DZZP,DOF_DRRP) = rn(DOF_DZ)**2
    p2l(DOF_DZZP,DOF_DZZP) = zn(DOF_DZ)**2
    p2l(DOF_DZZP,DOF_DRZP) = 2*zn(DOF_DZ)*rn(DOF_DZ)
    p2l(DOF_DZZP,DOF_DRP) = rn(DOF_DZZ)
    p2l(DOF_DZZP,DOF_DZP) = zn(DOF_DZZ)
    p2l(DOF_DZZP,DOF_DRR) = 2*rn(DOF_DZ)*rn(DOF_DZP)
    p2l(DOF_DZZP,DOF_DZZ) = 2*zn(DOF_DZ)*zn(DOF_DZP)
    p2l(DOF_DZZP,DOF_DRZ) = 2*(zn(DOF_DZP)*rn(DOF_DZ)+zn(DOF_DZ)*rn(DOF_DZP))
    p2l(DOF_DZZP,DOF_DR) = rn(DOF_DZZP)
    p2l(DOF_DZZP,DOF_DZ) = zn(DOF_DZZP)
    ! fxyz = Ry*Rx*fRRz + Zy*Zx*fZZz + (Rx*Zy+Ry*Zx)*fRZz
    !     + Rxy*fRz + Zxy*fZz + Rxyz*fR + Zxyz*fZ
    !     + (Ry*Rxz+Ryz*Rx)*fRR + (Zy*Zxz+Zyz*Zx)*fZZ 
    !     + (Rx*Zyz+Rxz*Zy+Ry*Zxz+Ryz*Zx)*fRZ
    p2l(DOF_DRZP,DOF_DRRP) = rn(DOF_DR)*rn(DOF_DZ)
    p2l(DOF_DRZP,DOF_DZZP) = zn(DOF_DR)*zn(DOF_DZ) 
    p2l(DOF_DRZP,DOF_DRZP) = rn(DOF_DR)*zn(DOF_DZ) + zn(DOF_DR)*rn(DOF_DZ) 
    p2l(DOF_DRZP,DOF_DRP) = rn(DOF_DRZ) 
    p2l(DOF_DRZP,DOF_DZP) = zn(DOF_DRZ)
    p2l(DOF_DRZP,DOF_DRR) = rn(DOF_DRP)*rn(DOF_DZ)+rn(DOF_DR)*rn(DOF_DZP)
    p2l(DOF_DRZP,DOF_DZZ) = zn(DOF_DRP)*zn(DOF_DZ)+zn(DOF_DR)*zn(DOF_DZP)  
    p2l(DOF_DRZP,DOF_DRZ) = rn(DOF_DRP)*zn(DOF_DZ) + zn(DOF_DRP)*rn(DOF_DZ) & 
                          + rn(DOF_DR)*zn(DOF_DZP) + zn(DOF_DR)*rn(DOF_DZP) 
    p2l(DOF_DRZP,DOF_DR) = rn(DOF_DRZP) 
    p2l(DOF_DRZP,DOF_DZ) = zn(DOF_DRZP)
#endif
  end subroutine p2l_matrix

#endif 

end module physical_mesh

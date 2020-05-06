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
  real :: rm1, rm2, zm1
  ! arrays to store dofs of rst & zst fields 
  real, allocatable :: rstnode(:,:)
  real, allocatable :: zstnode(:,:)

contains

  ! Return rst and zst given x, phi, z 
  elemental subroutine physical_geometry(rout, zout, x, phi, z)
    use math
    implicit none

    real, intent(in):: x, phi, z
    vectype, intent(out):: rout, zout 
    real :: r, theta, ds, r2n
    integer :: i, js
    real, dimension(mn_mode) :: rstc, zsts
    real :: m_max, dphi

    dphi = 1.*mesh_period/2
    m_max = 15.5
    r = sqrt((x - xcenter)**2 + (z - zcenter)**2 + 0e-6)
    theta = atan2(z - zcenter, x - xcenter)
!    rout = 5.4 + rm1*r*cos(theta+rm2*sin(theta)) 
!    zout = 0 + zm1*r*sin(theta) 
    rout = 0
    zout = 0
    if (iread_vmec==1) then ! use VMEC surfaces to calculate geometry
      r2n = r**2*(ns-1)
      js = ceiling(r2n)
      if (js>(ns-1)) js = ns-1 
      if (js>0) then ! interpolate between surfaces
        ds = js - r2n 
        rstc = rmnc(:,js+1)*(1-ds) + rmnc(:,js)*ds
        zsts = zmns(:,js+1)*(1-ds) + zmns(:,js)*ds
        do i = 1, mn_mode 
          if (xmv(i)<m_max) then
            rout = rout + rstc(i)*cos(xmv(i)*theta+xnv(i)*(phi+dphi))
            zout = zout + zsts(i)*sin(xmv(i)*theta+xnv(i)*(phi+dphi))
          end if
        end do
      else ! near axis, use routine 
        ds = sqrt(r2n) 
        rstc = rmnc(:,js+1)
        zsts = zmns(:,js+1)
        do i = 1, mn_mode 
          if (xmv(i)<m_max) then
            rout = rout + rstc(i)*cos(xmv(i)*theta+xnv(i)*(phi+dphi))*ds**xmv(i)
            zout = zout + zsts(i)*sin(xmv(i)*theta+xnv(i)*(phi+dphi))*ds**xmv(i)
          end if
        end do
      end if
    else ! use routine to generate from boundary geometry
      do i = 1, mn_mode 
        if ((mb(i)<m_max).and.(abs(nb(i))<5)) then
          rout = rout + rbc(i)*cos(mb(i)*theta-nb(i)*(phi+dphi)*nperiods)*r**mb(i)
          zout = zout + zbs(i)*sin(mb(i)*theta-nb(i)*(phi+dphi)*nperiods)*r**mb(i)
        end if
      end do
    end if
  end subroutine physical_geometry

  ! Calculate curvature and normal vector on physical boundary
  subroutine get_boundary_curv(normal, curv, inode)
    implicit none

    integer, intent(in) :: inode 
    real, intent(out) :: normal(2), curv
    real :: dr, ddr, dz, ddz, theta, x, phi, z
    real :: coords(3)
    integer :: i
    real :: m_max, dphi

    ! get logical coordinates    
    call m3dc1_node_getcoord(inode-1,coords)
    x = coords(1)
    z = coords(2)
    phi = coords(3)
   
!    rm1 = 1.
!    rm2 = .0
!    zm1 = 1.
    dphi = 1.*mesh_period/2
    m_max = 15.5
    theta = atan2(z - zcenter, x - xcenter)
!    dr = -rm1*sin(theta+rm2*sin(theta))*(1+rm2*cos(theta)) 
!    ddr = -rm1*cos(theta+rm2*sin(theta))*(1+rm2*cos(theta))**2 &
!          +rm1*sin(theta+rm2*sin(theta))*rm2*sin(theta)
!    dz = zm1*cos(theta) 
!    ddz = -zm1*sin(theta) 
    dr = 0
    dz = 0
    ddr = 0
    ddz = 0
    if (iread_vmec==1) then
      do i = 1, mn_mode 
        if (xmv(i)<m_max) then
          dr = dr - rbc(i)*sin(xmv(i)*theta+xnv(i)*(phi+dphi))*xmv(i)
          dz = dz + zbs(i)*cos(xmv(i)*theta+xnv(i)*(phi+dphi))*xmv(i)
          ddr = ddr - rbc(i)*cos(xmv(i)*theta+xnv(i)*(phi+dphi))*xmv(i)**2
          ddz = ddz - zbs(i)*sin(xmv(i)*theta+xnv(i)*(phi+dphi))*xmv(i)**2
        end if
      end do
    else
      do i = 1, mn_mode 
        if ((mb(i)<m_max).and.(abs(nb(i))<5)) then
          dr = dr - rbc(i)*sin(mb(i)*theta-nb(i)*(phi+dphi)*nperiods)*mb(i)
          dz = dz + zbs(i)*cos(mb(i)*theta-nb(i)*(phi+dphi)*nperiods)*mb(i)
          ddr = ddr - rbc(i)*cos(mb(i)*theta-nb(i)*(phi+dphi)*nperiods)*mb(i)**2
          ddz = ddz - zbs(i)*sin(mb(i)*theta-nb(i)*(phi+dphi)*nperiods)*mb(i)**2
        end if 
      end do
    end if 
    curv = (dr*ddz - dz*ddr)/((dr**2 + dz**2)*sqrt(dr**2 + dz**2))
    normal(1) = dz/sqrt(dr**2 + dz**2)
    normal(2) = -dr/sqrt(dr**2 + dz**2)

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
    p2l(DOF_DRZ,DOF_DRR) = rn(DOF_DZ)*rn(DOF_DR)
    p2l(DOF_DRZ,DOF_DZZ) = zn(DOF_DR)*zn(DOF_DZ) 
    p2l(DOF_DRZ,DOF_DRZ) = rn(DOF_DR)*zn(DOF_DZ) + zn(DOF_DR)*rn(DOF_DZ) 
    p2l(DOF_DRZ,DOF_DR) = rn(DOF_DRZ) 
    p2l(DOF_DRZ,DOF_DZ) = zn(DOF_DRZ)

  end subroutine p2l_matrix

#endif 

end module physical_mesh

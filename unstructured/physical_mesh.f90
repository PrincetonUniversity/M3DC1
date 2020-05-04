module physical_mesh
  use element
  use read_vmec

  implicit none

  real :: xcenter, zcenter ! coords of the center of the logical mesh
  integer :: igeometry   ! 0 = identity; 1 = prescribed; 2 = solve Laplace equation  
  integer :: iread_vmec 
  integer :: nperiods
#ifdef USEST
  real :: mesh_period 
  real :: rm1, rm2, zm1
  ! arrays to store dofs of rst & zst fields 
  real, allocatable :: rnode(:,:)
  real, allocatable :: znode(:,:)

contains

  ! Return rst and zst given x, phi, z 
  subroutine physical_geometry(rout, zout, x, phi, z)
    use math
    implicit none

    real, intent(in):: x, phi, z
    vectype, intent(out):: rout, zout 
    real :: r, theta, ds, r2n, m_max, dphi
    integer :: i, js

    dphi = 1.*mesh_period/2
    m_max = 12.5
!    rm1 = .8
!    rm2 = .3
!    zm1 = 1.2
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
      if (js>1) then ! interpolate between surfaces
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
    real :: dr, ddr, dz, ddz, theta, x, phi, z, m_max, dphi 
    real :: coords(3)
    integer :: i

    ! get logical coordinates    
    call m3dc1_node_getcoord(inode-1,coords)
    x = coords(1)
    z = coords(2)
    phi = coords(3)
   
    dphi = 1.*mesh_period/2
    m_max = 12.5
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
#endif 

end module physical_mesh

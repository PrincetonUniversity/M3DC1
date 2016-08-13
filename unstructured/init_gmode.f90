!============================================================================
! Gravitational Instability Equilibrium (itor = 0, itaylor = 5)
!============================================================================
module grav

contains

subroutine grav_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_local_vals(l)

     call get_node_pos(l, x, phi, z)

     call grav_equ(x, z)
     call grav_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine grav_init

subroutine grav_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: fac1, n0, pn0, ppn0

  n0 = exp(z/ln)
  pn0 = n0/ln
  ppn0 = pn0/ln

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  psi0_l = 0.

  fac1 = gravz*ln+gam*p0
  bz0_l(1) = sqrt(bzero**2 - 2.*fac1*(n0-1.))
  bz0_l(2) = 0.
  bz0_l(3) = -pn0*fac1/bz0_l(1)
  bz0_l(4) = 0.
  bz0_l(5) = 0.
  bz0_l(6) = (fac1/bz0_l(1))*(pn0*bz0_l(3)/bz0_l(1) - ppn0)

  
  fac1 = p0 - pi0
  pe0_l(1) = fac1+gam*fac1*(n0-1.)
  pe0_l(2) = 0.
  pe0_l(3) = gam*fac1*pn0
  pe0_l(4) = 0.
  pe0_l(5) = 0.
  pe0_l(6) = gam*fac1*ppn0

  den0_l(1) = n0
  den0_l(2) = 0.
  den0_l(3) = pn0
  den0_l(4) = 0.
  den0_l(5) = 0.
  den0_l(6) = ppn0

  fac1 = p0
  p0_l(1) = fac1+gam*fac1*(n0-1.)
  p0_l(2) = 0.
  p0_l(3) = gam*fac1*pn0
  p0_l(4) = 0.
  p0_l(5) = 0.
  p0_l(6) = gam*fac1*ppn0

end subroutine grav_equ


subroutine grav_per(x, z)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in) :: x, z

  real :: kx, kz, alx, alz

  call get_bounding_box_size(alx, alz)
  kx = twopi/alx
  kz = pi/alz

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.

  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.
  p1_l = 0.

  den1_l(1) =  eps*sin(kx*(x-xzero))*sin(kz*(z-zzero))
  den1_l(2) =  eps*cos(kx*(x-xzero))*sin(kz*(z-zzero))*kx
  den1_l(3) =  eps*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kz
  den1_l(4) = -eps*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kx**2
  den1_l(5) =  eps*cos(kx*(x-xzero))*cos(kz*(z-zzero))*kx*kz
  den1_l(6) = -eps*sin(kx*(x-xzero))*sin(kz*(z-zzero))*kz**2

end subroutine grav_per

end module grav

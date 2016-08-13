!==============================================================================
! Force-Free Taylor State (itaylor = 2)
!==============================================================================
module force_free_state

contains

subroutine force_free_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz

  call get_bounding_box_size(alx, alz)

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)

     call get_node_pos(l, x, phi, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call get_local_vals(l)

     call force_free_equ(x, z)
     call force_free_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine force_free_init

subroutine force_free_equ(x, z)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in) :: x, z
  
  real :: kx, kz, alx, alz, alam
  
  call get_bounding_box_size(alx, alz)

  kx = pi/alx
  kz = pi/alz

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  alam = sqrt(kx**2+kz**2)

  psi0_l(1) = sin(kx*x)*sin(kz*z)
  psi0_l(2) = kx*cos(kx*x)*sin(kz*z)
  psi0_l(3) = kz*sin(kx*x)*cos(kz*z)
  psi0_l(4) = -kx**2*sin(kx*x)*sin(kz*z)
  psi0_l(5) =  kx*kz*cos(kx*x)*cos(kz*z)
  psi0_l(6) = -kz**2*sin(kz*x)*sin(kz*z)

  bz0_l = alam*psi0_l

  call constant_field( pe0_l, p0-pi0)
  call constant_field(  p0_l, p0)
  call constant_field(den0_l, 1.)
  
end subroutine force_free_equ

subroutine force_free_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.

  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.

  p1_l = 0.
  den1_l = 0.

end subroutine force_free_per

end module force_free_state

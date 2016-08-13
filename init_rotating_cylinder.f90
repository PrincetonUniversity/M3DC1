!==============================================================================
! Rotating cylinder
! ~~~~~~~~~~~~~~~~~
!
! This is a rotating equilibrium with a radially increasing density to offset
! the centrifugal force.  In the isothermal case, with do density diffusion or
! thermal diffusion, this is a static equilibrium 
! (with or without gyroviscosity).
!==============================================================================
module rotate

  real, private :: kx, kz

contains

subroutine rotate_init()
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

     z = z - alz*.5

     call get_local_vals(l)

     call rotate_equ(x, z)
     call rotate_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine rotate_init

subroutine rotate_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z


  u0_l = 0.
  vz0_l(1) =    vzero*x**2
  vz0_l(2) = 2.*vzero*x
  vz0_l(3) = 0.
  vz0_l(4) = 2.*vzero
  vz0_l(5) = 0.
  vz0_l(6) = 0.
  chi0_l = 0.

  psi0_l = 0.
  bz0_l(1) =    bzero*x**2
  bz0_l(2) = 2.*bzero*x
  bz0_l(3) = 0.
  bz0_l(4) = 2.*bzero
  bz0_l(5) = 0.
  bz0_l(6) = 0.
  pe0_l(1) = p0 - pi0
  pe0_l(2:6) = 0.
  p0_l(1) = p0
  p0_l(2:6) = 0.
  den0_l(1) =  2.*(bzero/vzero)**2
  den0_l(2:6) = 0.
  

end subroutine rotate_equ


subroutine rotate_per(x, z)
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
  den1_l = 0.
  p1_l = 0.

end subroutine rotate_per

end module rotate

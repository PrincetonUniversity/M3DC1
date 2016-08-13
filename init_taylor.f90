!==============================================================================
! Taylor Reconnection (itaylor = 1)
!==============================================================================
module taylor_reconnection

contains

subroutine taylor_reconnection_init()
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

     call taylor_reconnection_equ(x, z)
     call taylor_reconnection_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine taylor_reconnection_init



subroutine taylor_reconnection_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  psi0_l(1) = -0.5*z**2
  psi0_l(2) =  0.
  psi0_l(3) = -z
  psi0_l(4) =  0.
  psi0_l(5) =  0.
  psi0_l(6) = -1.

  call constant_field( bz0_l, bzero)
  call constant_field( pe0_l, p0-pi0)
  call constant_field(den0_l, 1.)
  
end subroutine taylor_reconnection_equ

subroutine taylor_reconnection_per(x, z)
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

end subroutine taylor_reconnection_per


end module taylor_reconnection

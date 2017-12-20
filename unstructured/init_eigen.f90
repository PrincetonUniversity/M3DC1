!==============================================================================
! eigen
! ~~~~~~~~~~~~
!==============================================================================
module eigen

!    real, private :: kx, kz

contains

!========================================================
! init
!========================================================
subroutine eigen_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z


  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     call eigen_equ(x-rzero, z)
     call eigen_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)
end subroutine eigen_init


!========================================================
! equ
!========================================================
subroutine eigen_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
 


  call constant_field(den0_l, 1.)
  call constant_field(p0_l, 0.1)
  call constant_field(pe0_l, 0.05)
  call constant_field(psi0_l, 0.)
  call constant_field(bz0_l, rzero)
  


end subroutine eigen_equ


!========================================================
! per
!========================================================
subroutine eigen_per(x, phi, z)

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask
!!$
!!$
!!$
!!$
     vmask = 1.
!     vmask(1:6) = /0.05
!     vmask(1) = vmask(1) !- pedge/p0
     
     ! initial parallel rotation
     u1_l = phizero*vmask

     ! allow for initial toroidal rotation
     vz1_l = 0.
!     if(vzero.ne.0) call add_angular_velocity(vz1_l, x+xzero, vzero*vmask)

     ! add random perturbations
     if(nonrect.eq.0) then
        vmask(1) = 1.
        vmask(2:6) = 0.
     endif
     call random_per(x,phi,z,vmask)
end subroutine eigen_per

end module eigen

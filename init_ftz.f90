!==============================================================================
! ftz
! ~~~~~~~~~~~~
!==============================================================================
module ftz

!    real, private :: kx, kz

contains

!========================================================
! init
!========================================================
subroutine ftz_init()
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

     call ftz_equ(x-rzero, z)
     call ftz_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)
end subroutine ftz_init


!========================================================
! equ
!========================================================
subroutine ftz_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r, j0

  call constant_field(den0_l, 1.)
  
  
  j0=4./rzero/2.8

  r=sqrt(x**2 + z**2)

  
  p0_l(1)=p0*(1.- r**2)
  p0_l(2)=p0*(- 2.*x)
  p0_l(3)=p0*(- 2.*z)
  p0_l(4)=p0*(- 2.)
  p0_l(5)=0.
  p0_l(6)=p0*(- 2.)

  pe0_l(1)=0.5*p0_l(1)
  pe0_l(2)=0.5*p0_l(2)
  pe0_l(3)=0.5*p0_l(3)
  pe0_l(4)=0.5*p0_l(4)
  pe0_l(5)=0.5*p0_l(5)
  pe0_l(6)=0.5*p0_l(6)

  psi0_l(1) = -rzero*j0*(r**2/4.-r**4/16.-3./16.)
  psi0_l(2) = -rzero*j0*(1./2.-r**2/4.)*x
  psi0_l(3) = -rzero*j0*(1./2.-r**2/4.)*z
  psi0_l(4) = -rzero*j0*(-1./2.*x**2+1./2.-r**2/4.)
  psi0_l(5) = 1./2.*rzero*j0*x*z
  psi0_l(6) =  -rzero*j0*(-1./2.*z**2+1./2.-r**2/4.)
  




  bz0_l(1)=rzero
  bz0_l(2)=0.
  bz0_l(3)=0.
  bz0_l(4)=0.
  bz0_l(5)=0.
  bz0_l(6)=0.

  

end subroutine ftz_equ


!========================================================
! per
!========================================================
subroutine ftz_per(x, phi, z)

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
end subroutine ftz_per

end module ftz

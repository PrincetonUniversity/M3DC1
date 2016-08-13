!==============================================================================
! 3D diffusion Test
! ~~~~~~~~~~~~
! This is a test case which initializes a 2-dimensional equilibrium
! With a 3-dimensional initial perturbation
!==============================================================================
module threed_diffusion_test

  real, private :: kx, kz

contains

subroutine threed_diffusion_test_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, x1, x2, z1, z2
  real :: phi0, x0, z0

  call get_bounding_box(x1, z1, x2, z2)
  kx = pi/(x2-x1)
  kz = pi/(z2-z1)

#ifdef USE3D
  phi0 = pi
#else
  phi0 = 0
#endif
  z0 = (z1 + z2)/2.
  x0 = (x1 + 2.*x2)/3.

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     call threed_diffusion_test_equ(x-x1, phi, z-z1)
     call threed_diffusion_test_per(x-x0, phi-phi0, z-z0)

     call set_local_vals(l)
  enddo
end subroutine threed_diffusion_test_init

subroutine threed_diffusion_test_equ(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z

  psi0_l(1) = bx0*sin(kx*x)*sin(kz*z)
  psi0_l(2) = bx0*cos(kx*x)*sin(kz*z)*kx
  psi0_l(3) = bx0*sin(kx*x)*cos(kz*z)*kz
  psi0_l(4) =-bx0*sin(kx*x)*sin(kz*z)*kx**2
  psi0_l(5) = bx0*cos(kx*x)*cos(kz*z)*kx*kz
  psi0_l(6) =-bx0*sin(kx*x)*sin(kz*z)*kz**2

  call constant_field(bz0_l, bzero)
  call constant_field(den0_l, 1.)
  call constant_field(p0_l, p0)
  call constant_field(pe0_l, p0-pi0)

end subroutine threed_diffusion_test_equ


subroutine threed_diffusion_test_per(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z
  real :: ex

  ex = eps*exp(-(rzero*phi/ln)**2-(x/ln)**2-(z/ln)**2)

  p1_l(1) = ex
  p1_l(2) = -2.*ex*x/ln**2
  p1_l(3) = -2.*ex*z/ln**2
  p1_l(4) = 4.*ex*x**2/ln**4 - 2.*ex*ln**2
  p1_l(5) = 4.*ex*x*z/ln**4
  p1_l(6) = 4.*ex*z**2/ln**4 - 2.*ex*ln**2
#ifdef USE3D
  p1_l(7:12) = -2.*phi*p1_l(1:6)*(rzero/ln)**2
#endif

  if(idens.eq.1) den1_l = p1_l
  pe1_l = p1_l*pefac

end subroutine threed_diffusion_test_per

end module threed_diffusion_test

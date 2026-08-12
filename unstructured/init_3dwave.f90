!==============================================================================
! 3D wave Test
! ~~~~~~~~~~~~
! This is a test case which initializes a 2-dimensional equilibrium
! With a 3-dimensional initial perturbation
!==============================================================================
module threed_wave_test

  real, private :: kx, kz, kphi, omega, k2, kp2, kB, b2

contains

subroutine threed_wave_test_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, x1, x2, z1, z2

  call get_bounding_box(x1, z1, x2, z2)
  kx = pi/(x2-x1)
  kz = pi/(z2-z1)
  kphi = ntor/rzero
  kp2 = kx**2 + kz**2
  k2 = kp2 + kphi**2
  kB = kx*Bx0 + kphi*bzero
  b2 = bzero**2 + Bx0**2

  select case(numvar)
  case(1)
     select case(iwave)
     case(0)
        omega = kB

     case default
        if(myrank.eq.0) then
           print *, 'Error: iwave must be one of the following:'
           print *, ' iwave = 0: Shear Alfven wave'
        end if
        call safestop(4)
     end select
     
  case(2)
     select case(iwave)
     case(0)
        if(myrank.eq.0) print *, 'Shear Alfven wave'
        
        omega = bx0*kx*kB + 0.5*b2*kphi**2 &
             + 0.5*sqrt((2.*bx0*kx*kB + B2*kphi**2)**2 &
             -4.*(bx0*kx*kB)**2*k2/kp2)
        omega = sqrt(omega)

     case(1)
        if(myrank.eq.0) print *, 'Fast magnetosonic wave'

        omega = bx0*kx*kB + 0.5*b2*kphi**2 &
             - 0.5*sqrt((2.*bx0*kx*kB + B2*kphi**2)**2 &
             -4.*(bx0*kx*kB)**2*k2/kp2)
        omega = sqrt(omega)

     case default
        if(myrank.eq.0) then
           print *, 'Error: iwave must be one of the following:'
           print *, ' iwave = 0: Shear Alfven wave'
           print *, ' iwave = 1: Fast wave'
        end if
        call safestop(4)
     end select

  case(3)
     if(bx0 .ne. 0.) then
        if(myrank.eq.0) print *, 'Only Bx0=0 supported for numvar=3'
        call safestop(4)
     endif

     select case(iwave)
     case(0)
        if(myrank.eq.0) print *, 'Shear Alfven wave'
        omega = bzero*kphi

     case(1)
        if(myrank.eq.0) print *, 'Fast magnetosonic wave'
        omega = 0.5*k2*(bzero**2 + gam*p0) &
             + 0.5*sqrt((k2*(bzero**2 + gam*p0))**2 &
             -4.*(bzero*kphi)**2*k2*gam*p0)
        omega = sqrt(omega)

     case(2)
        if(myrank.eq.0) print *, 'Slow magnetosonic wave'
        omega = 0.5*k2*(bzero**2 + gam*p0) &
             - 0.5*sqrt((k2*(bzero**2 + gam*p0))**2 &
             -4.*(bzero*kphi)**2*k2*gam*p0)
        omega = sqrt(omega)

     case default
        if(myrank.eq.0) then
           print *, 'Error: iwave must be one of the following:'
           print *, ' iwave = 0: Shear Alfven wave'
           print *, ' iwave = 1: Fast wave'
           print *, ' iwave = 2: Slow wave'
        end if
        call safestop(4)
     end select
  end select

  if(myrank.eq.0) then
     print *, ' kphi = ', myrank, kphi
     print *, ' k2 = ', myrank, k2
     print *, ' kp2 = ', myrank, kp2
     print *, ' omega = ', myrank, omega
     print *, ' rzero = ', myrank, rzero
     print *, ' bzero = ', myrank, bzero
     print *, ' p0 = ', myrank, p0
  end if

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     x = x - x1
     z = z - z1

     call get_local_vals(l)

     call threed_wave_test_equ(x, phi, z)
     call threed_wave_test_per(x, phi, z)

     call set_local_vals(l)
  enddo
end subroutine threed_wave_test_init

subroutine threed_wave_test_equ(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z

  psi0_l = 0.

  call constant_field(bz0_l, bzero)
  call constant_field(den0_l, 1.)
  call constant_field(p0_l, p0)
  call constant_field(pe0_l, p0-pi0)

  psi0_l(:) = 0.
  psi0_l(1) = -Bx0*z
  psi0_l(3) = -Bx0

end subroutine threed_wave_test_equ


subroutine threed_wave_test_per(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z
  vectype, dimension(dofs_per_node) :: val

  val(1) = eps*sin(kx*x)*sin(kz*z)*cos(kphi*phi)
  val(2) = eps*cos(kx*x)*sin(kz*z)*cos(kphi*phi)*kx
  val(3) = eps*sin(kx*x)*cos(kz*z)*cos(kphi*phi)*kz
  val(4) =-eps*sin(kx*x)*sin(kz*z)*cos(kphi*phi)*kx**2
  val(5) = eps*cos(kx*x)*cos(kz*z)*cos(kphi*phi)*kx*kz
  val(6) =-eps*sin(kx*x)*sin(kz*z)*cos(kphi*phi)*kz**2
#ifdef USE3D
  val(7)  =-eps*sin(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi
  val(8)  =-eps*cos(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi*kx
  val(9)  =-eps*sin(kx*x)*cos(kz*z)*sin(kphi*phi)*kphi*kz
  val(10) = eps*sin(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi*kx**2
  val(11) =-eps*cos(kx*x)*cos(kz*z)*sin(kphi*phi)*kphi*kx*kz
  val(12) = eps*sin(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi*kz**2
#endif

  select case(numvar)
  case(1)
     psi1_l = val
     u1_l = -val*kB/omega

  case(2)
     psi1_l = val
     u1_l = -val*kB/omega
     vz1_l = -val*(0,1)*Bx0*kz*kphi*kp2 &
          / (k2*(Bx0*kx)**2 - kp2*omega**2)
     bz1_l = (-kp2)* (val*(0,1)*Bx0**2*kx*kz*kphi) &
          / (k2*(Bx0*kx)**2 - kp2*omega**2)
     
  case(3)
     if(iwave.eq.0) then
        psi1_l = val
        u1_l = -val*bzero*kphi/omega
     else
        bz1_l =-val*kp2
        vz1_l = val*bzero*kphi*kp2*k2*p0*gam / &
             (omega*(k2*p0*gam - omega**2))
#if defined(USE3D)
        ! unlike other fields, chi has -sin(phi) dependence
        chi1_l(1:6) = (val(7:12)/kphi)* &
             (0,1)*bzero*k2*(kphi**2*p0*gam - omega**2) / &
             (omega*(k2*p0*gam - omega**2))
        ! toroidal derivative has -cos(phi) dependence
        chi1_l(7:12) = (-val(1:6)*kphi)* &
             (0,1)*bzero*k2*(kphi**2*p0*gam - omega**2) / &
             (omega*(k2*p0*gam - omega**2))
#elif defined(USECOMPLEX)
        chi1_l = val*(0,1)*bzero*k2*(kphi**2*p0*gam - omega**2) / &
             (omega*(k2*p0*gam - omega**2))
#endif
        p1_l = val*bzero*kp2*k2*p0*gam / &
             (k2*p0*gam - omega**2)
     endif
  end select

end subroutine threed_wave_test_per

end module threed_wave_test

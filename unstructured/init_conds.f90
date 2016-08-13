!==============================
subroutine constant_field(outarr, val)
  use element

  implicit none

  vectype, dimension(dofs_per_node), intent(out) :: outarr
  real, intent(in) :: val

  outarr(1) = val
  outarr(2:6) = 0.
#ifdef USE3D
  outarr(7:12) = 0.
#endif

end subroutine constant_field
!==============================
subroutine plane_wave(outarr, x, z, kx, kz, amp, phase)
  implicit none

  real, intent(in) :: x, z
  vectype, dimension(6), intent(out) :: outarr
  real, intent(in)  :: kx, kz, amp, phase

  real :: arg
  arg = kx*x + kz*z + phase

  outarr(1) =  amp*cos(arg)
  outarr(2) = -amp*sin(arg)*kx
  outarr(3) = -amp*sin(arg)*kz
  outarr(4) = -amp*cos(arg)*kx*kx
  outarr(5) = -amp*cos(arg)*kx*kz
  outarr(6) = -amp*cos(arg)*kz*kz
end subroutine plane_wave
!==============================
subroutine plane_wave2(outarr,x,phi,z,kx,kphi,kz,amp,phasex,phasephi,phasez)
  use element

  implicit none

  vectype, dimension(dofs_per_node), intent(out) :: outarr
  real, intent(in)  :: x, phi, z, kx, kphi, kz, amp, phasex, phasephi, phasez

  real :: argx,argp,argz,cox,cop,coz,six,sip,siz
  argx = kx*x     + phasex
  argp = kphi*phi + phasephi
  argz = kz*z     + phasez

  six = sin(argx)
  cox = cos(argx)
  sip = sin(argp)
  cop = cos(argp)
  siz = sin(argz)
  coz = cos(argz)

  outarr(1) =  amp*cop*six*siz
  outarr(2) =  amp*cop*cox*siz*kx
  outarr(3) =  amp*cop*six*coz*kz
  outarr(4) = -amp*cop*six*siz*kx*kx
  outarr(5) =  amp*cop*cox*coz*kx*kz
  outarr(6) = -amp*cop*six*siz*kz*kz
#ifdef USE3D
  outarr(7 ) = -amp*sip*six*siz*kphi
  outarr(8 ) = -amp*sip*cox*siz*kx*kphi
  outarr(9 ) = -amp*sip*six*coz*kz*kphi
  outarr(10) =  amp*sip*six*siz*kx*kx*kphi
  outarr(11) = -amp*sip*cox*coz*kx*kz*kphi
  outarr(12) =  amp*sip*six*siz*kz*kz*kphi
#endif
end subroutine plane_wave2
!==============================
subroutine add_angular_velocity(outarr, x,omega)
  use arrays

  implicit none

  vectype, dimension(6), intent(inout) :: outarr
  vectype, dimension(6), intent(in) :: omega
  real, intent(in) :: x

  vectype, dimension(6) :: temp

  temp(1) = x**2 * omega(1)
  temp(2) = x**2 * omega(2) + 2.*x*omega(1)
  temp(3) = x**2 * omega(3)
  temp(4) = x**2 * omega(4) + 4.*x*omega(2) + 2.*omega(1)
  temp(5) = x**2 * omega(5) + 2.*x*omega(3)
  temp(6) = x**2 * omega(6)

  outarr = outarr + temp

end subroutine add_angular_velocity
!=============================
subroutine add_product(x,a,b)
  vectype, dimension(6), intent(in) :: a,b
  vectype, dimension(6), intent(inout) :: x

  x(1) = x(1) + a(1)*b(1)
  x(2) = x(2) + a(1)*b(2) + a(2)*b(1)
  x(3) = x(3) + a(1)*b(3) + a(3)*b(1)
  x(4) = x(4) + a(1)*b(4) + a(4)*b(1) + 2.*a(2)*b(2)
  x(5) = x(5) + a(1)*b(5) + a(5)*b(1) + a(2)*b(3) + a(3)*b(2)
  x(6) = x(6) + a(1)*b(6) + a(6)*b(1) + 2.*a(3)*b(3)
end subroutine add_product
!=============================
subroutine random_per(x,phi,z,fac)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in) :: x, phi, z
  vectype, intent(in), dimension(dofs_per_node) :: fac
  integer, allocatable :: seed(:)
  integer :: i, j, n
  real :: alx, alz, kx, kp, kz, xx, zz, random, rsq, r,ri,roundoff,ri3,rexp,co,sn
  vectype, dimension(dofs_per_node) :: temp

  call get_bounding_box_size(alx, alz)
!
! changed to be consistent with fortran95 7/12/2011 ... scj
  call random_seed(SIZE = n)
  allocate(seed(n))
  seed = 23
  call random_seed(PUT = seed)
  deallocate(seed)

  temp = 0.
  roundoff = 1.e-12

  xx = x - xzero
  zz = z - zzero
  rsq = xx**2 + zz**2 + roundoff
  r = sqrt(rsq)
  ri = 1./sqrt(rsq + roundoff)
  ri3 = ri/rsq
  rexp = exp(-rsq/ln)
  co = cos(phi)
  sn = sin(phi)

  do i=1,maxn
     kx = pi*i/alx
     select case (icsym)

     case (0)   !  original option...no symmetry imposed
     do j=1, maxn
        kz = j*pi/alz
        kp = j
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(psi1_l,fac,temp)
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(u1_l,fac,temp)
     end do

     case (1)  !   make U odd symmetry about midplane:  perturb only U
     do j=1, maxn/2
        kz = 2.*j*pi/alz
        kp = j
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(u1_l,fac,temp)
     end do

     case (2)  !   make U even  symmetry about midplane:  perturb only U
     do j=1, maxn/2
        kz = (2.*j-1)*pi/alz
        kp = j
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(u1_l,fac,temp)
     end do

     case (3)  !   NOT RANDOM....start in (1,1) eigenfunction
     temp(1) = eps* r * rexp*(zz*co - xx*sn)
     temp(2) = eps* ri * rexp*(zz*xx*co - xx*xx*sn)   &
             - eps*(2./ln)* r * rexp*(zz*xx*co - xx*xx*sn)    &
             - eps* r * rexp*sn 
     temp(3) = eps* ri * rexp*(zz*zz*co - xx*zz*sn)   &
             - eps*(2./ln)* r * rexp*(zz*zz*co - xx*zz*sn)    &
             + eps* r * rexp*co
     temp(4) = -eps* ri3 * rexp*(zz*xx*xx*co - xx*xx*xx*sn)   &
             + eps*(2./ln)* r * rexp*((2.*zz*xx*xx/ln - zz)*co - (2*xx*xx*xx/ln - 3.*xx)*sn)    &
             + eps* ri * rexp*((zz - 4.*zz*xx*xx/ln)*co - (3*xx-4*xx**3/ln)*sn)               
          
     temp(5) = -eps* ri3 * rexp*(zz*zz*xx*co - zz*xx*xx*sn)   &
             + eps*(2./ln)* r * rexp*((2.*zz*zz*xx/ln - xx)*co - (2*zz*xx*xx/ln - zz)*sn)    &
             + eps* ri * rexp*((xx - 4.*zz*zz*xx/ln)*co - (zz-4*zz*xx**2/ln)*sn)               
         
     temp(6) = -eps* ri3 * rexp*(zz*zz*zz*co - xx*zz*zz*sn)   &
             + eps*(2./ln)* r * rexp*((2.*zz*zz*zz/ln - 3*zz)*co - (2*xx*zz*zz/ln - xx)*sn)    &
             + eps* ri * rexp*((3*zz - 4.*zz*zz*zz/ln)*co - (xx-4*xx*zz*zz/ln)*sn)   
#ifdef USE3D
     temp(7) = eps* r * rexp*(-zz*sn - xx*co)
     temp(8) = eps* ri * rexp*(-zz*xx*sn - xx*xx*co)   &
             - eps*(2./ln)* r * rexp*(-zz*xx*sn - xx*xx*co)    &
             - eps* r * rexp*co 
     temp(9) = eps* ri * rexp*(-zz*zz*sn - xx*zz*co)   &
             - eps*(2./ln)* r * rexp*(-zz*zz*sn - xx*zz*co)    &
             - eps* r * rexp*sn
     temp(10) = -eps* ri3 * rexp*(-zz*xx*xx*sn - xx*xx*xx*co)   &
             + eps*(2./ln)* r * rexp*(-(2.*zz*xx*xx/ln - zz)*sn - (2*xx*xx*xx/ln - 3.*xx)*co)    &
             + eps* ri * rexp*(-(zz - 4.*zz*xx*xx/ln)*sn - (3*xx-4*xx**3/ln)*co)               
          
     temp(11) = -eps* ri3 * rexp*(-zz*zz*xx*sn - zz*xx*xx*co)   &
             + eps*(2./ln)* r * rexp*(-(2.*zz*zz*xx/ln - xx)*sn - (2*zz*xx*xx/ln - zz)*co)    &
             + eps* ri * rexp*(-(xx - 4.*zz*zz*xx/ln)*sn - (zz-4*zz*xx**2/ln)*co)               
         
     temp(12) = -eps* ri3 * rexp*(-zz*zz*zz*sn - xx*zz*zz*co)   &
             + eps*(2./ln)* r * rexp*(-(2.*zz*zz*zz/ln - 3*zz)*sn - (2*xx*zz*zz/ln - xx)*co)    &
             + eps* ri * rexp*(-(3*zz - 4.*zz*zz*zz/ln)*sn - (xx-4*xx*zz*zz/ln)*co)   
#endif               
 
     call add_product(u1_l,fac,temp)

     end select
  end do

end subroutine random_per
!===========================
subroutine cartesian_to_cylindrical(x,vec)
  implicit none

  vectype, dimension(6), intent(inout) :: vec
  real, intent(in) :: x

  vec(6) = vec(6) * x
  vec(5) = vec(5) * x +    vec(3)
  vec(4) = vec(4) * x + 2.*vec(2)
  vec(3) = vec(3) * x
  vec(2) = vec(2) * x +    vec(1)
  vec(1) = vec(1) * x 
end subroutine cartesian_to_cylindrical
!===========================
subroutine cartesian_to_cylindrical_all()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: inode, numnodes, icounter_tt
  real :: x, phi, z

  numnodes = owned_nodes()

   do icounter_tt=1,numnodes
     inode = nodes_owned(icounter_tt)

     call get_node_pos(inode, x, phi, z)

     call get_local_vals(inode)

     call cartesian_to_cylindrical(x,psi0_l)
     call cartesian_to_cylindrical(x,psi1_l)
     call cartesian_to_cylindrical(x,  u0_l)
     call cartesian_to_cylindrical(x,  u1_l)
     
     if(numvar.ge.2) then
        call cartesian_to_cylindrical(x,bz0_l)
        call cartesian_to_cylindrical(x,bz1_l)
        call cartesian_to_cylindrical(x,vz0_l)
        call cartesian_to_cylindrical(x,vz1_l)
     endif
     
     call set_local_vals(inode)
  end do
end subroutine cartesian_to_cylindrical_all


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



!==============================================================================
! frs
! ~~~~~~~~~~~~
!==============================================================================
module frs

!    real, private :: kx, kz

contains

!========================================================
! init
!========================================================
subroutine frs_init()
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

     call frs_equ(x-rzero, z)
     call frs_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)

end subroutine frs_init


!========================================================
! equ
!========================================================
subroutine frs_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r,Bp0,integral,rs, Bz_edge,Bz,dBzdx,dBzdz,r0
  integer :: m,n


  r0=xlim ! Current channel with
  rs=zlim! Position of singular surface

  m=mpol 
  n=ntor
  Bz_edge=1. 


  Bp0=sqrt((Bz_edge**2 -2.*p0*(1.-rs**2)) &
       /((real(m)/real(n)*rzero*(rs/r0)/(1+(rs/r0)**2)/rs )**2 &
       +2.*(-1./2.)*((1.+(rs/r0)**2)**(-2)-(1.+(1./r0)**2)**(-2))))

  call constant_field(den0_l, 1.)

  r=sqrt(x**2 + z**2)

  
  p0_l(1)=p0*(1.- r**2)
  p0_l(2)=p0*(- 2.*x)
  p0_l(3)=p0*(- 2.*z)
  p0_l(4)=p0*(- 2.)
  p0_l(5)=0.
  p0_l(6)=p0*(- 2.)



  pe0_l(1)=0.
  pe0_l(2)=0.
  pe0_l(3)=0.
  pe0_l(4)=0.
  pe0_l(5)=0.
  pe0_l(6)=0.


  psi0_l(1) = (-rzero*Bp0*r0/2.*log(1.+(r/r0)**2) &
       +rzero*Bp0*r0/2.*log(1.+(1./r0)**2))
  psi0_l(2) = -rzero*Bp0*(x/r0)/(1.+(r/r0)**2)
  psi0_l(3) = -rzero*Bp0*(z/r0)/(1.+(r/r0)**2)
  psi0_l(4) = -rzero*Bp0*r0/2.*(-(2.*x/r0**2)**2/(1+(r/r0)**2)**2 &
       +2./r0**2/(1+(r/r0)**2))
  psi0_l(5) = -rzero*Bp0*r0/2.*(-4.*z*x)/r0**4/(1+(r/r0)**2)**2
  psi0_l(6) = -rzero*Bp0*r0/2.*(-(2.*z/r0**2)**2/(1+(r/r0)**2)**2 &
       +2./r0**2/(1+(r/r0)**2))
  


  integral=(-1./2.)*((1.+(r/r0)**2)**(-2)-(1.+(1./r0)**2)**(-2))

  Bz=sqrt(Bz_edge**2-2.*p0_l(1)-2.*Bp0**2*integral)

  dBzdx=(1./Bz)*(-Bp0**2*2.*x/r0**2/(1.+(r/r0)**2)**3+2.*p0*x)
  dBzdz=(1./Bz)*(-Bp0**2*2.*z/r0**2/(1.+(r/r0)**2)**3+2.*p0*z)
  bz0_l(1)=rzero*Bz
  bz0_l(2)=rzero*dBzdx
  bz0_l(3)=rzero*dBzdz
  bz0_l(4)=rzero* &
       (-1./Bz**2*dBzdx*(-Bp0**2*2.*x/r0**2/(1.+(r/r0)**2)**3+2.*p0*x) &
       +1./Bz*(-Bp0**2*(2./r0**2/(1.+(r/r0)**2)**3 &
                        - 3.*(2.*x/r0**2)**2/(1.+(r/r0)**2)**4  )+2.*p0))
  bz0_l(5)=rzero* &
       (-1./Bz**2*dBzdz*(-Bp0**2*2.*x/r0**2/(1.+(r/r0)**2)**3+2.*p0*x) &
       +1./Bz*(-Bp0**2*(-3.)*4.*x*z/r0**4/(1.+(r/r0)**2)**4))
  bz0_l(6)=rzero* &
       (-1./Bz**2*dBzdz*(-Bp0**2*2.*z/r0**2/(1.+(r/r0)**2)**3+2.*p0*z) &
       +1./Bz*(-Bp0**2*(2./r0**2/(1.+(r/r0)**2)**3 &
                        - 3.*(2.*z/r0**2)**2/(1.+(r/r0)**2)**4  )+2.*p0))

end subroutine frs_equ

subroutine frs1_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes
  real :: x, phi, z


  numnodes = owned_nodes()
  do l=1, numnodes
     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     call frs1_equ(x-rzero, z)
     call frs_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo

end subroutine frs1_init

!========================================================
! equ
!========================================================
subroutine frs1_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r, Bz, Bz_edge

  Bz_edge=rzero

  call constant_field(den0_l, 1.)

  r=sqrt(x**2 + z**2)

  
  call constant_field(p0_l,p0)

  psi0_l(1) = -.5/q0*(1.-.5*r**2)**2
  psi0_l(2) =   x/q0*(1.-.5*r**2)
  psi0_l(3) =   z/q0*(1.-.5*r**2)
  psi0_l(4) =  1./q0*(1.-.5*r**2-x**2)
  psi0_l(5) = -1./q0*x*z
  psi0_l(6) =  1./q0*(1.-.5*r**2-z**2)
  


  Bz=sqrt(Bz_edge**2+4./q0**2*(5./24. - r**2/2. + 3.*r**4/8. - r**6/12))

  bz0_l(1)= Bz
  bz0_l(2)= 2*x/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)
  bz0_l(3)= 2*z/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)
  bz0_l(4)= -4*x**2/(bz**3*q0**4)*(-1+3*r**2/2 -r**4/2.)**2   &
          + 2/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)   &
          + 2*x**2/(bz*q0**2)*(3-2*r)
  bz0_l(5)=  -4*x*z/(bz**3*q0**4)*(-1+3*r**2/2 -r**4/2.)**2   &
          + 2*x*z/(bz*q0**2)*(3-2*r)
  bz0_l(6)= -4*z**2/(bz**3*q0**4)*(-1+3*r**2/2 -r**4/2.)**2   &
          + 2/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)   &
          + 2*z**2/(bz*q0**2)*(3-2*r)

end subroutine frs1_equ

!========================================================
! per
!========================================================
subroutine frs_per(x, phi, z)

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  integer :: i
  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask

     vmask = 1.
     vmask(1:6) = p0_l(1:6)/0.05

     
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
end subroutine frs_per

end module frs

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

!==============================================================================
! int_kink
! ~~~~~~~~~~~~
!==============================================================================
module int_kink


contains

!========================================================
! init
!========================================================
subroutine int_kink_init()
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

     call int_kink_equ(x, z)
     call int_kink_per(x, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)

end subroutine int_kink_init


!========================================================
! equ
!========================================================
subroutine int_kink_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r, preszero
  real, parameter :: alpha=0.6579, rmin=1., pi=3.1415926536
  
  preszero=p0
  r=sqrt(x**2 + z**2)

  if (NUMVAR.eq.1) then
    call constant_field(den0_l, 1.)
    call constant_field(p0_l,p0)
    call constant_field(pe0_l,p0/2.)
    call constant_field(bz0_l,bzero)
  else ! only for rmin=1.

    bz0_l(1) = 0.225079* &
         Sqrt(1.80383 - 0.957984*preszero - 0.108208*(x**2 + z**2) &
         + 3.22684*preszero*(x**2 + z**2) + 0.162312*(x**2 + z**2)**2 &
         - 4.40951*preszero*(x**2 + z**2)**2 - 0.120231*(x**2 + z**2)**3 &
         + 3.28718*preszero*(x**2 + z**2)**3 + 0.0450867*(x**2 + z**2)**4 &
         - 1.5192*preszero*(x**2 + z**2)**4 - 0.00721387*(x**2 + z**2)**5 &
         + 0.45427*preszero*(x**2 + z**2)**5 &
         - 0.0952705*preszero*(x**2 + z**2)**6 &
         + 0.0159547*preszero*(x**2 + z**2)**7 &
         - 0.00245458*preszero*(x**2 + z**2)**8 &
         + 0.00018182*preszero*(x**2 + z**2)**9)

    bz0_l(2) = (0.11253953951963827* &
         (0. - 0.21641620500000003*x + 6.45367005*preszero*x &
         + 0.6492486150000001*x*(x**2 + z**2) &
         - 17.638059524400003*preszero*x*(x**2 + z**2) &
         - 0.7213873500000001*x*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*x*(x**2 + z**2)**2 &
         + 0.36069367500000005*x*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*x*(x**2 + z**2)**3 &
         - 0.07213873500000001*x*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*x*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*x*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*x*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*x*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*x*(x**2 + z**2)**8)) &
         / Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(3) = (0.11253953951963827* &
         (0. - 0.21641620500000003*z + 6.45367005*preszero*z &
         + 0.6492486150000001*z*(x**2 + z**2) &
         - 17.638059524400003*preszero*z*(x**2 + z**2) &
         - 0.7213873500000001*z*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*z*(x**2 + z**2)**2 &
         + 0.36069367500000005*z*(x**2 + z**2)**3 & 
         - 12.153599550456903*preszero*z*(x**2 + z**2)**3 &
         - 0.07213873500000001*z*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*z*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*z*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*z*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*z*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*z*(x**2 + z**2)**8))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(4) = (-0.056269769759819135* &
         (0. - 0.21641620500000003*x + 6.45367005*preszero*x &
         + 0.6492486150000001*x*(x**2 + z**2) &
         - 17.638059524400003*preszero*x*(x**2 + z**2) &
         - 0.7213873500000001*x*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*x*(x**2 + z**2)**2 &
         + 0.36069367500000005*x*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*x*(x**2 + z**2)**3 &
         - 0.07213873500000001*x*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*x*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*x*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*x*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*x*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*x*(x**2 + z**2)**8)**2)/ &
         (1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 & 
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)**1.5 &
         + (0.11253953951963827*&
         (-0.21641620500000003 + 6.45367005*preszero &
         + 1.2984972300000002*x**2 - 35.276119048800005*preszero*x**2 &
         + 0.6492486150000001*(x**2 + z**2) &
         - 17.638059524400003*preszero*(x**2 + z**2) &
         - 2.8855494000000004*x**2*(x**2 + z**2) &
         + 78.89236034011381*preszero*x**2*(x**2 + z**2) &
         - 0.7213873500000001*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*(x**2 + z**2)**2 &
         + 2.1641620500000003*x**2*(x**2 + z**2)**2 &
         - 72.92159730274142*preszero*x**2*(x**2 + z**2)**2 &
         + 0.36069367500000005*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*(x**2 + z**2)**3 &
         - 0.5771098800000001*x**2*(x**2 + z**2)**3 &
         + 36.341613059299846*preszero*x**2*(x**2 + z**2)**3 &
         - 0.07213873500000001*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*(x**2 + z**2)**4 &
         - 11.432456270094134*preszero*x**2*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*(x**2 + z**2)**5 &
         + 2.680396653036631*preszero*x**2*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*(x**2 + z**2)**6 &
         - 0.5498249544690526*preszero*x**2*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*(x**2 + z**2)**7 &
         + 0.052364281378004994*preszero*x**2*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*(x**2 + z**2)**8))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(5) = (-0.056269769759819135* &
         (0. - 0.21641620500000003*x + 6.45367005*preszero*x &
         + 0.6492486150000001*x*(x**2 + z**2) &
         - 17.638059524400003*preszero*x*(x**2 + z**2) &
         - 0.7213873500000001*x*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*x*(x**2 + z**2)**2 &
         + 0.36069367500000005*x*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*x*(x**2 + z**2)**3 &
         - 0.07213873500000001*x*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*x*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*x*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*x*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*x*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*x*(x**2 + z**2)**8)* &
         (0. - 0.21641620500000003*z + 6.45367005*preszero*z &
         + 0.6492486150000001*z*(x**2 + z**2) &
         - 17.638059524400003*preszero*z*(x**2 + z**2) &
         - 0.7213873500000001*z*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*z*(x**2 + z**2)**2 &
         + 0.36069367500000005*z*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*z*(x**2 + z**2)**3 &
         - 0.07213873500000001*z*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*z*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*z*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*z*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*z*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*z*(x**2 + z**2)**8))/ &
         (1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)**1.5 &
         +(0.11253953951963827* &
         (1.2984972300000002*x*z - 35.276119048800005*preszero*x*z &
         - 2.8855494000000004*x*z*(x**2 + z**2) &
         + 78.89236034011381*preszero*x*z*(x**2 + z**2) &
         + 2.1641620500000003*x*z*(x**2 + z**2)**2 &
         - 72.92159730274142*preszero*x*z*(x**2 + z**2)**2 &
         - 0.5771098800000001*x*z*(x**2 + z**2)**3 &
         + 36.341613059299846*preszero*x*z*(x**2 + z**2)**3 &
         - 11.432456270094134*preszero*x*z*(x**2 + z**2)**4 &
         + 2.680396653036631*preszero*x*z*(x**2 + z**2)**5 &
         - 0.5498249544690526*preszero*x*z*(x**2 + z**2)**6 &
         + 0.052364281378004994*preszero*x*z*(x**2 + z**2)**7))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(6) = (-0.056269769759819135* &
         (0. - 0.21641620500000003*z + 6.45367005*preszero*z &
         + 0.6492486150000001*z*(x**2 + z**2) &
         - 17.638059524400003*preszero*z*(x**2 + z**2) &
         - 0.7213873500000001*z*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*z*(x**2 + z**2)**2 &
         + 0.36069367500000005*z*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*z*(x**2 + z**2)**3 &
         - 0.07213873500000001*z*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*z*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*z*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*z*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*z*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*z*(x**2 + z**2)**8)**2)/ &
         (1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)**1.5 &
         + (0.11253953951963827* &
         (-0.21641620500000003 + 6.45367005*preszero &
         + 1.2984972300000002*z**2 - 35.276119048800005*preszero*z**2 &
         + 0.6492486150000001*(x**2 + z**2) &
         - 17.638059524400003*preszero*(x**2 + z**2) &
         - 2.8855494000000004*z**2*(x**2 + z**2) &
         + 78.89236034011381*preszero*z**2*(x**2 + z**2) &
         - 0.7213873500000001*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*(x**2 + z**2)**2 &
         + 2.1641620500000003*z**2*(x**2 + z**2)**2 &
         - 72.92159730274142*preszero*z**2*(x**2 + z**2)**2 &
         + 0.36069367500000005*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*(x**2 + z**2)**3 &
         - 0.5771098800000001*z**2*(x**2 + z**2)**3 &
         + 36.341613059299846*preszero*z**2*(x**2 + z**2)**3 &
         - 0.07213873500000001*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*(x**2 + z**2)**4 &
         - 11.432456270094134*preszero*z**2*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*(x**2 + z**2)**5 &
         + 2.680396653036631*preszero*z**2*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*(x**2 + z**2)**6 &
         - 0.5498249544690526*preszero*z**2*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*(x**2 + z**2)**7 &
         + 0.052364281378004994*preszero*z**2*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*(x**2 + z**2)**8))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 & 
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    p0_l(1) =preszero*(1. - 3.106*(x**2 + z**2) + 3.54*(x**2 + z**2)**2 - 1.408*(x**2 + z**2)**3)                         
    p0_l(2) =preszero*(-6.212*x + 14.16*x*(x**2 + z**2) - 8.448*x*(x**2 + z**2)**2)                                       
    p0_l(3) =preszero*(-6.212*z + 14.16*z*(x**2 + z**2) - 8.448*z*(x**2 + z**2)**2)                                       
    p0_l(4) =preszero*(-6.212 + 28.32*(x**2) + 14.16*(x**2 + z**2) - 33.792*(x**2)*(x**2 + z**2) - 8.448*(x**2 + z**2)**2) 
    p0_l(5) =preszero*(28.32*x*z - 33.792*x*z*(x**2 + z**2))
    p0_l(6) =preszero*(-6.212 + 28.32*(z**2) + 14.16*(x**2 + z**2) - 33.792*(z**2)*(x**2 + z**2) - 8.448*(x**2 + z**2)**2) 

    pe0_l(1) = p0_l(1)/2.
    pe0_l(2) = p0_l(2)/2.
    pe0_l(3) = p0_l(3)/2.
    pe0_l(4) = p0_l(4)/2.
    pe0_l(5) = p0_l(5)/2.
    pe0_l(6) = p0_l(6)/2.

    den0_l(1) =1. - 1.6*(x**2 + z**2) + 0.8*(x**2 + z**2)**2 - 3.677e-8*(x**2 + z**2)**3
    den0_l(2) =-3.2*x + 3.2*x*(x**2 + z**2) - 2.2062e-7*x*(x**2 + z**2)**2
    den0_l(3) =-3.2*z + 3.2*z*(x**2 + z**2) - 2.2062e-7*z*(x**2 + z**2)**2
    den0_l(4) =-3.2 + 6.4*(x**2) + 3.2*(x**2 + z**2) - 8.8248e-7*(x**2)*(x**2 + z**2) - 2.2062e-7*(x**2 + z**2)**2
    den0_l(5) =6.4*x*z - 8.8248e-7*x*z*(x**2 + z**2)
    den0_l(6) =-3.2 + 6.4*(z**2) + 3.2*(x**2 + z**2) - 8.8248e-7*(z**2)*(x**2 + z**2) - 2.2062e-7*(x**2 + z**2)**2

  end if

  call constant_field(u0_l, 0.)
  call constant_field(vz0_l, 0.)
  call constant_field(chi0_l, 0.)

  psi0_l(1) = alpha * ( r**2/4. - r**4/(8. * rmin**2)               + r**6/(36. * rmin**4)    )
  psi0_l(2) = alpha * ( x/2.    - r**2 * x/(2. * rmin**2)           + r**4 * x/(6. * rmin**4) )
  psi0_l(3) = alpha * ( z/2.    - r**2 * z/(2. * rmin**2)           + r**4 * z/(6. * rmin**4) )
  psi0_l(4) = alpha * ( 1./2.   - (3. * x**2 + z**2)/(2. * rmin**2) + (5. * x**4 + 6. * x**2 * z**2 + z**4)/(6. * rmin**4) )
  psi0_l(5) = alpha * (         - x * z / rmin**2                   + (x**3 * z + x * z**3) * 2./(3. * rmin**4) )
  psi0_l(6) = alpha * ( 1./2.   - (x**2 + 3. * z**2)/(2. * rmin**2) + (x**4 + 6. * x**2 * z**2 + 5. * z**4)/(6. * rmin**4) )
 

end subroutine int_kink_equ


!========================================================
! per
!========================================================
subroutine int_kink_per(x, phi, z)

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask

     ! add random perturbations
        vmask(1) = 1.
        vmask(2:dofs_per_node) = 0.
     call random_per(x,phi,z,vmask)





end subroutine int_kink_per

end module int_kink


!=========================================
! set_neo_vel
!=========================================
subroutine set_neo_vel
  use basic
  use arrays
  use m3dc1_nint
  use newvar_mod
  use neo
  use math
  use sparse
  use model

  implicit none

  integer :: i, j, itri, nelms, ier

  type(field_type) :: vz_vec, u_f, chi_f, diamag
  type(vector_type) :: vp_vec
  type(matrix_type) :: vpol_mat
  vectype, dimension(dofs_per_element, dofs_per_element, 2, 2) :: temp
  vectype, dimension(dofs_per_element) :: temp2, temp4
  vectype, dimension(dofs_per_element, 2) :: temp3

  real, dimension(MAX_PTS) :: theta, vtor, vpol, psival
  vectype, dimension(MAX_PTS) :: vz, vp, dia

  integer, dimension(dofs_per_element) :: imask_vor, imask_chi
  integer, dimension(MAX_PTS) :: iout
  integer :: imag, magnetic_region

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Setting velocity from NEO data"

     if(db.ne.0. .and. ineo_subtract_diamag.eq.1) call create_field(diamag)
  call create_field(vz_vec)
  call create_vector(vp_vec, 2)
  call associate_field(u_f, vp_vec, 1)
  call associate_field(chi_f, vp_vec, 2)

  call set_matrix_index(vpol_mat, vpol_mat_index)
  call create_mat(vpol_mat, 2, 2, icomplex, 1)

  nelms = local_elements()
  do itri=1,nelms
        
     call define_element_quadrature(itri,int_pts_main,5)
     call define_fields(itri,0,1,0)
     call eval_ops(itri, psi_field(0), ps079)
     if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then 
        call eval_ops(itri, den_field(0), n079)
        call eval_ops(itri, p_field(0), p079)
        call eval_ops(itri, pe_field(0), pe079)
        pi079 = p079 - pe079

        ! zeff does not appear here because db includes zeff
        dia = db*(pi079(:,OP_DR)*ps079(:,OP_DR)+pi079(:,OP_DZ)*ps079(:,OP_DZ))&
             / n079(:,OP_1)
     end if

     theta = atan2(z_79-zmag,x_79-xmag)
     psival = -(ps079(:,OP_1) - psimin)
     call neo_eval_vel(int_pts_main, psival, theta, vpol, vtor, iout)
     vz = vtor / (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)
     vp = vpol / (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)

     ! NEO coordinates are (r, theta, phi_neo) = (r, phi, theta)
     ! --> theta = grad(r)xgrad(phi)
     ! if grad(psi).grad(r) < 0 then vpol is opposite direction to Bpol
     if(psimin.gt.psibound) vp = -vp
     ! phi_neo = -phi
     vz = -vz

     temp79e = sqrt((ps079(:,OP_DR)**2 + ps079(:,OP_DZ)**2)*ri2_79)

     do i=1, int_pts_main
        imag = magnetic_region(ps079(i,:), x_79(i), z_79(i))
        if(imag.ne.0) then
           vz(i) = 0.
           vp(i) = 0.
           iout(i) = 1
        endif
     end do

     where(iout.eq.1 .or. abs(temp79e).lt.1e-2)
        temp79f = 0.
        dia = 0.
     elsewhere
        temp79f = 1./temp79e
        dia = dia / &
             (ps079(:,OP_DR)*ps079(:,OP_DR) + ps079(:,OP_DZ)*ps079(:,OP_DZ))
     end where

     call get_vor_mask(itri, imask_vor)
     call get_chi_mask(itri, imask_chi)

     do i=1,dofs_per_element          
        ! assemble matrix
        do j=1, dofs_per_element

           ! vorticity equation
           if(imask_vor(i).eq.0) then
              temp(i,j,1,:) = 0.
           else
              temp(i,j,1,1) = &
                   int3(r2_79,mu79(:,OP_DR,i),nu79(:,OP_DR,j)) + &
                   int3(r2_79,mu79(:,OP_DZ,i),nu79(:,OP_DZ,j))
              temp(i,j,1,2) = &
                   int3(ri_79,mu79(:,OP_DR,i),nu79(:,OP_DZ,j)) - &
                   int3(ri_79,mu79(:,OP_DZ,i),nu79(:,OP_DR,j))
           end if
           
           ! compression equation
           if(imask_chi(i).eq.0) then
              temp(i,j,2,:) = 0.
           else
              temp(i,j,2,1) = &
                   int3(ri_79,mu79(:,OP_DZ,i),nu79(:,OP_DR,j)) - &
                   int3(ri_79,mu79(:,OP_DR,i),nu79(:,OP_DZ,j))
              temp(i,j,2,2) = &
                   int3(ri4_79,mu79(:,OP_DR,i),nu79(:,OP_DR,j)) + &
                   int3(ri4_79,mu79(:,OP_DZ,i),nu79(:,OP_DZ,j))
           end if
        end do

        ! assemble RHS
        ! toroidal rotation
        select case(ivform)
        case(0)
           temp2(i) = int3(r_79,mu79(:,OP_1,i),vz)
        case(1)
           temp2(i) = int3(ri_79,mu79(:,OP_1,i),vz)
        end select

        ! vorticity
        if(imask_vor(i).eq.0) then
           temp3(i,1) = 0.
        else
           temp3(i,1) = &
                int4(temp79f,vp,mu79(:,OP_DR,i),ps079(:,OP_DR)) + &
                int4(temp79f,vp,mu79(:,OP_DZ,i),ps079(:,OP_DZ))
        endif

        ! compression
        if(imask_vor(i).eq.0) then
           temp3(i,2) = 0.
        else
           temp3(i,2) = &
                int5(ri3_79,temp79f,vp,mu79(:,OP_DZ,i),ps079(:,OP_DR)) - &
                int5(ri3_79,temp79f,vp,mu79(:,OP_DR,i),ps079(:,OP_DZ))
        endif

        ! diamagnetic term
        if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then
           temp4(i) = int2(mu79(:,OP_1,i), dia)
        end if
     end do

     call insert_block(vpol_mat, itri, 1, 1, temp(:,:,1,1), MAT_ADD)
     call insert_block(vpol_mat, itri, 1, 2, temp(:,:,1,2), MAT_ADD)
     call insert_block(vpol_mat, itri, 2, 1, temp(:,:,2,1), MAT_ADD)
     call insert_block(vpol_mat, itri, 2, 2, temp(:,:,2,2), MAT_ADD)
     call vector_insert_block(vp_vec, itri, 1, temp3(:,1), VEC_ADD)
     call vector_insert_block(vp_vec, itri, 2, temp3(:,2), VEC_ADD)
     call vector_insert_block(vz_vec%vec, itri, 1, temp2(:), VEC_ADD)

     if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then
        call vector_insert_block(diamag%vec, itri, 1, temp4(:), VEC_ADD)
     end if
  end do
  call sum_shared(vz_vec%vec)
  call sum_shared(vp_vec)
  if(db.ne.0. .and. ineo_subtract_diamag.eq.1) call sum_shared(diamag%vec)

  ! solve vz
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving vz..."
  call newsolve(mass_mat_lhs%mat,vz_vec%vec,ier)

  ! add neoclassical rotation to base rotation
  call add_field_to_field(vz_field(0), vz_vec)

  ! subtract diamagnetic term from rotation (because neo is adding this in)
  if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then
     ! solve diamagnetic term
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving diamag..."
     call newsolve(mass_mat_lhs%mat,diamag%vec,ier)
     call add_field_to_field(vz_field(0), diamag, -1.)
  end if

  ! solve vpol
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving vpol..."
  call boundary_vpol(vp_vec, u_f, chi_f, vpol_mat)
  call finalize(vpol_mat)
  call newsolve(vpol_mat, vp_vec, ier)
  u_field(0) = u_f
  chi_field(0) = chi_f

  call destroy_field(vz_vec)
  call destroy_vector(vp_vec)
  call destroy_mat(vpol_mat)
  if(db.ne.0. .and. ineo_subtract_diamag.eq.1) call destroy_field(diamag)

  if(myrank.eq.0 .and. iprint.ge.1) print *, " Done calculating velocity"
end subroutine set_neo_vel


!=====================================
subroutine initial_conditions()
  use basic
  use vector_mod
  use arrays

  use tilting_cylinder
  use taylor_reconnection
  use force_free_state
  use gem_reconnection
  use wave_propagation
  use gradshafranov
  use mri
  use grav
  use strauss
  use circular_field
  use rotate
  use eqdsk_eq
  use dskbal_eq
  use jsolver_eq
  use biharmonic
  use circ_shell_only
  use resistive_wall_test
  use threed_wave_test
  use threed_diffusion_test
  use frs
  use ftz
  use eigen
  use neo
  use int_kink
  use rwm
  use solovev
  use auxiliary_fields
  use basicq
  use basicj
  use rmp
  use init_common

  implicit none

  integer :: ierr

  if(iread_neo.eq.1) then
     call read_neo(ierr)
     if(ierr.ne.0) return
  end if

  if(iread_eqdsk.ge.1) then
     call eqdsk_init()
  else if(iread_dskbal.ge.1) then
     call dskbal_init()
  else if(iread_jsolver.ge.1) then
     call jsolver_init()
  else
     if(itor.eq.0) then
        ! slab equilibria
        select case(itaylor)
        case(0)
           call tilting_cylinder_init()
        case(1)
           call taylor_reconnection_init()
        case(2)
           call force_free_init()
        case(3)
           call gem_reconnection_init()
        case(4)
           call wave_init()
        case(5)
           call grav_init()
        case(6)
           call strauss_init()
        case(7)
           call circular_field_init()
        case(8)
           call biharmonic_init(1)
        case(9)
           call biharmonic_init(0)
        case(10,11)
           call circ_shell_only_init()
        case(12,13)
           call resistive_wall_test_init()
        case(14)
           call threed_wave_test_init()
        case(15)
           call threed_diffusion_test_init()
        case(16)
           call frs_init()
        case(17)
           call ftz_init()
        case(18)
           call eigen_init()
        case(19)
           call int_kink_init()
        case(20)
           call kstar_profiles()
        case(21,22,25,26,27,28)
           call fixed_q_profiles()
        case(23)
           call frs1_init()
        case(24)
           call rwm_init()
        case(29)
           call basicj_init()
        end select
     else
        ! toroidal equilibria
        select case(itaylor)
        case(0)
           call tilting_cylinder_init()
           call cartesian_to_cylindrical_all()
        case(1)
           call gradshafranov_init()
        case(2)
           call mri_init()
        case(3)
           call rotate_init()
        case(7)
           call circular_field_init()
           call cartesian_to_cylindrical_all()
        case(10,11)
           call circ_shell_only_init()
           call cartesian_to_cylindrical_all()
        case(12,13)
           call resistive_wall_test_init()
           call cartesian_to_cylindrical_all()
        case(14)
           call threed_wave_test_init()
        case(15)
           call threed_diffusion_test_init()
        case(16)
           call frs_init()
        case(17)
           call ftz_init()
        case(18)
           call eigen_init()
        case(19)
           call solovev_init()
        case(24)
           call rwm_init()
        end select
     endif
  end if

  if(iread_neo.eq.1) then
     call set_neo_vel
     call unload_neo
  end if

  call den_eq()
  call den_per()

  if(irmp.ge.1 .or. iread_ext_field.ge.1 .or. &
       tf_tilt.ne.0. .or. tf_shift.ne.0. .or. &
       any(pf_tilt.ne.0.) .or. any(pf_shift.ne.0.)) call rmp_per()

  ! calculate equilibrium and perturbed temperature profiles
  call calculate_temperatures(0, te_field(0),ti_field(0), 1)
  call calculate_temperatures(1, te_field(1),ti_field(1), 1)

  if(iflip_b.eq.1) call mult(bz_field(0), -1.)
  if(iflip_j.gt.0) call mult(psi_field(0), -1.)
  if(iflip_v.eq.1) call mult(vz_field(0), -1.)
  if(iflip_v.eq.-1) call mult(vz_field(0), 0.)
end subroutine initial_conditions
                                                                     
subroutine kstar_profiles()

  use basic
  use math
  use mesh_mod
  use sparse
  use arrays
  use m3dc1_nint
  use newvar_mod
  use boundary_conditions
  use model
  use gradshafranov
  use int_kink

  vectype, dimension (dofs_per_element,dofs_per_element) :: temp
  vectype, dimension (dofs_per_element) :: temp2
  vectype, dimension (MAX_PTS) :: co, sn, r, theta, rdpsidr
  real :: x, phi, z
  real, parameter :: e=2.7183
  real :: a, r1, r2, u, fa, ra
  real :: b0,r0
  integer :: m,n
  integer :: numnodes, nelms, l, itri, i, j, ier, icounter_tt
  integer :: imask(dofs_per_element)
  type (field_type) :: psi_vec
  type(matrix_type) :: psi_mat
  
  !call create_field(dpsi_dr)
  call create_field(psi_vec)
  
  call set_matrix_index(psi_mat, psi_mat_index)
  call create_mat(psi_mat,1,1,icomplex, 1)

  psi_vec = 0.

  !input variables: B0,fA,r1,r2,q0,R0,m,n,u,rA
  b0 = bzero
  r0 = rzero
  m = mpol
  n = igs
  a = alpha0
  r1 = alpha1
  r2 = alpha2
  u = alpha3
  fa = p1
  ra = p2
  
  if(myrank.eq.0 .and. iprint.ge.1) write(*,*)   &
                                    'b0,r0,m,n,e,a,r1,r2,u,fa,ra'   &
                                    ,b0,r0,m,n,e,a,r1,r2,u,fa,ra

  numnodes = owned_nodes()

  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l,x,phi,z)
     
     call constant_field(den0_l,1.)
     call constant_field(bz0_l,bzero)
     call constant_field(p0_l,p0)
     call constant_field(pe0_l,p0/2.)
     call int_kink_per(x,phi,z)
     
     call set_node_data(den_field(0),l,den0_l)
     call set_node_data(bz_field(0),l,bz0_l)
     call set_node_data(p_field(0),l,p0_l)
     call set_node_data(pe_field(0),l,pe0_l)
     call set_node_data(u_field(1),l,u1_l)
  enddo

  nelms = local_elements()
  do itri=1,nelms
     
     call define_element_quadrature(itri,int_pts_diag, int_pts_tor)
     call define_fields(itri,0,1,0) ! defines x_79,z_79,mu,nu
     !   call eval_ops(itri,dpsi_dr,dpt79)
     
     r = sqrt((x_79-xmag)**2 + (z_79-zmag)**2)
     theta = atan2(z_79-zmag,x_79-xmag)
     co = cos(theta)
     sn = sin(theta)

     rdpsidr = (B0*r**2)/ &
          ((1+fA/E**((-r1+r/a)**2/r2**2))*q0*R0* &
          (1+((r*Abs(-1+(m/(n*q0))**u)**(1/(2.*u)))/(a*rA))**(2*u))**(1/u))
      
     call get_boundary_mask(itri, BOUNDARY_DIRICHLET, imask, domain_boundary)

     !  assemble matrix    
     do i=1,dofs_per_element
        if(imask(i).eq.0) then
           temp(i,:) = 0.
        else
           do j=1,dofs_per_element
              temp(i,j) = int4(mu79(:,OP_1,i),nu79(:,OP_DR,j),co(:),r(:)) &
                   +      int4(mu79(:,OP_1,i),nu79(:,OP_DZ,j),sn(:),r(:))
           enddo
        end if
        !  assemble rhs
        temp2(i) = int2(mu79(:,OP_1,i),rdpsidr)
     enddo
     
     call insert_block(psi_mat, itri, 1,1, temp(:,:), MAT_ADD)
     
     call vector_insert_block(psi_vec%vec, itri, 1, temp2(:), MAT_ADD)
  enddo
  
  call sum_shared(psi_vec%vec)
  call flush(psi_mat)
  call boundary_gs(psi_vec%vec, 0., psi_mat)
  call finalize(psi_mat)

! solve for psi
  if(myrank.eq.0 .and. iprint.ge.1) print *, "solving psi"
 
  call newsolve(psi_mat,psi_vec%vec,ier)
  if(eqsubtract.eq.1) then
     psi_field(0) = psi_vec
  else
     psi_field(1) = psi_vec
  endif
  
  call destroy_mat(psi_mat)
  call destroy_field(psi_vec)
  !call destroy_field(dpsi_dr)
  
end subroutine kstar_profiles

module LZeqbm
  IMPLICIT NONE

  !Parameters
  REAL, PARAMETER :: aspect_ratio = 18.0  ! R_major / r_wall_b
  REAL, PARAMETER :: B0 = 1.0             ! Toroidal field at r=0
  REAL, PARAMETER :: dfrac = 0.1          ! Bndry layer width / r_plas_a
  REAL, PARAMETER :: Trat = 100.0         ! Ratio of plas to vac temperature
  REAL, PARAMETER :: nrat = 1.0           ! Ratio of plas to vac density
  REAL, PARAMETER :: q_interior = 1.0     ! Safety factor across bulk plasma
  REAL, PARAMETER :: q_a = 0.75           ! Safety factor on plasma surface
  REAL, PARAMETER :: r_plas_a = 0.6       ! Plasma minor radius
  REAL, PARAMETER :: r_wall_b = 1.0       ! Minor radius of ideal wall
  REAL, PARAMETER :: r_tile_c = 0.7       ! Minor radius of tile surface

  !Derived quantities
  REAL R_major                            ! Major radius of magnetic axis
  REAL delta                              ! Width of boundary layer
  REAL r_interior                         ! Radius to interior of bndry layer
  REAL prat                               ! Ratio pf plas to vac pressure
  REAL p_interior                         ! Pressure in bulk plasma

  REAL p_vacuum                           ! Pressure in vacuum region
  REAL q0R0                               ! Central safety factor * major radius
  REAL Jdenom                             ! Intermediate derived quantity
  REAL Bbz                                ! Toroidal field in bndry layer
  REAL gamma                              ! Pressure gradient in bndry layer
end module LZeqbm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Profiles for a pseudo-VMEC equilibrium structure resembling L. Zakharov's    !
! circular cross-section, straight-cylinder, flat q, ideal equilibrium with a  !
! surface current.                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setupLZeqbm
    use lzeqbm
    IMPLICIT NONE

    delta = dfrac * r_plas_a
    r_interior = r_plas_a - delta
    R_major = aspect_ratio * r_wall_b
    q0R0 = q_interior * R_major
    Jdenom = q0R0**2 + r_interior**2
    Bbz = q0R0**2 * B0 / Jdenom
    gamma = (1.5 * (q0R0**2 * B0**2 * r_interior**4 / Jdenom**2 - &
         r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / &
         (r_interior**3 - r_plas_a**3)) ! /aspect_ratio
    prat = nrat * Trat
    p_vacuum = gamma*delta/(prat - 1.0)
    p_interior = p_vacuum + gamma*delta
  end subroutine setupLZeqbm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function p_LZ(psi)
    use lzeqbm
    implicit none
    REAL p_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       p_LZ = p_interior
    ELSE IF (rmin < r_plas_a) THEN
       p_LZ = p_interior - gamma*(rmin - r_interior)
    ELSE
       p_LZ = p_vacuum
    endIF
  end function p_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pprime_LZ(psi)
    use lzeqbm
    implicit none
    REAL pprime_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, dpdr

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       pprime_LZ = 0.0
       RETURN
    endIF

    IF (rmin < r_plas_a) THEN
       dpdr = -gamma
    ELSE
       dpdr = 0.0
    endIF

    pprime_LZ = 0.5*dpdr/rmin
  end function pprime_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function q_LZ(psi)
    use lzeqbm
    implicit none
    intrinsic sqrt
    REAL q_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, Btheta


    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       q_LZ = q_interior
    ELSE IF (rmin < r_plas_a) THEN
       Btheta = SQRT(2.0*gamma*(rmin**3 - r_plas_a**3)/3.0 + &
            r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / rmin
       q_LZ = rmin * Bbz / (R_major * Btheta)
    ELSE
       q_LZ = q_a * (rmin/r_plas_a)**2
    endIF
  end function q_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function qprime_LZ(psi)
    use lzeqbm
    implicit none
    REAL qprime_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, Btheta, Bthetaprime, dqdr

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       qprime_LZ = 0.0
       RETURN
    endIF

    IF (rmin < r_plas_a) THEN
       Btheta = SQRT(2.0*gamma*(rmin**3 - r_plas_a**3)/3.0 + &
            r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / rmin
       dqdr = Bbz * (2.0 - gamma*rmin/Btheta**2) / (R_major * Btheta)
    ELSE
       dqdr = 2.0*rmin*q_a / r_plas_a**2
    endIF

    qprime_LZ = 0.5*dqdr/rmin
  end function qprime_LZ

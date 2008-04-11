!==============================
subroutine constant_field(outarr, val)
  implicit none

  vectype, dimension(6), intent(out) :: outarr
  vectype, intent(in) :: val

  outarr(1) = val
  outarr(2:6) = 0.

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
subroutine plane_wave2(outarr,x,z,kx,kz,amp,phasex,phasez)
  implicit none

  vectype, dimension(6), intent(out) :: outarr
  real, intent(in)  :: x, z, kx, kz, amp, phasex, phasez

  real :: argx,argz,cox,coz,six,siz
  argx = kx*x + phasex
  argz = kz*z + phasez

  six = sin(argx)
  cox = cos(argx)
  siz = sin(argz)
  coz = cos(argz)

  outarr(1) =  amp*six*siz
  outarr(2) =  amp*cox*siz*kx
  outarr(3) =  amp*six*coz*kz
  outarr(4) = -amp*six*siz*kx*kx
  outarr(5) =  amp*cox*coz*kx*kz
  outarr(6) = -amp*six*siz*kz*kz
end subroutine plane_wave2
!==============================
subroutine random_per(x,z,seed)
  use basic
  use arrays

  implicit none

  real :: rand

  real, intent(in) :: x, z
  integer, intent(in) :: seed
  integer :: i, j
  integer, parameter :: maxn = 10
  real :: alx, alz, a, kx, kz
  vectype, dimension(6) :: temp

  call getboundingboxsize(alx, alz)

  call srand(seed)

  psi1_l = 0
  u1_l = 0
  temp = 0

#ifdef _AIX
#define RAND_ARG
#else
#define RAND_ARG 0
#endif

  do i=1,maxn
     kx = pi*i/alx
     do j=1, maxn
        kz = pi*j/alz
        call plane_wave2(temp,x,z,kx,kz,2.*eps*(rand(RAND_ARG)-.5),0.,0.)
        psi1_l = psi1_l + temp
        call plane_wave2(temp,x,z,kx,kz,2.*eps*(rand(RAND_ARG)-.5),0.,0.)
        u1_l = u1_l + temp
     end do
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

  implicit none

  integer :: inode, numnodes
  real :: x
  double precision :: coords(3)

  call numnod(numnodes)

  do inode=1, numnodes
     call xyznod(inode, coords)
     x = coords(1) + xzero

     call assign_local_pointers(inode)

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
     
  end do
end subroutine cartesian_to_cylindrical_all


!==============================================================================
! Tilting Cylinder (itaylor = 0)
!==============================================================================
module tilting_cylinder

  implicit none

  real, parameter :: k = 3.8317059702

contains

subroutine tilting_cylinder_init()
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     call assign_local_pointers(l)

     call cylinder_equ(x, z)
     call cylinder_per(x, z)
  enddo

end subroutine tilting_cylinder_init

subroutine cylinder_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: rr, ri, arg, befo, ff, fp, fpp, j0, j1, kb
  real :: s17aef, s17aff
  integer :: ifail1, ifail2, ifail3
  
  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  !.....Use the equilibrium given in R. Richard, et al,
  !     Phys. Fluids B, 2 (3) 1990 p 489
  !.....(also Strauss and Longcope)
  !
  rr = sqrt(x**2 + z**2)
  ri = 0.
  if(rr.ne.0) ri = 1./rr
  arg = k*rr
  if(rr.le.1) then
     befo = 2./s17aef(k,ifail1)
     j0 = s17aef(arg,ifail2)
     j1 = s17aff(arg,ifail3)
     ff = .5*befo
     fp = 0
     fpp = -.125*k**2*befo
     if(arg.ne.0)  then
        ff   = befo*j1/arg
        fp  = befo*(j0 - 2.*j1/arg)/rr
        fpp = befo*(-3*j0 + 6.*j1/arg - arg*j1)/rr**2
     endif
  else
     ff   = (1.-1./rr**2)
     fp  = 2./rr**3
     fpp = -6/rr**4
  endif
   
  psi0_l(1) = ff* z
  psi0_l(2) = fp*z*x *ri
  psi0_l(3) = fp*z**2*ri + ff
  psi0_l(4) = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)
  psi0_l(5) = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
  psi0_l(6) = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)

  if(rr.gt.1) then
     call constant_field(bz0_l, bzero   )
     call constant_field(pe0_l, p0-pi0*ipres)
     call constant_field( p0_l, p0)
  else 
     kb = k**2*(1.-beta)
     bz0_l(1) = sqrt(kb*psi0_l(1)**2+bzero**2)
     bz0_l(2) = kb/bz0_l(1)*psi0_l(1)*psi0_l(2)
     bz0_l(3) = kb/bz0_l(1)*psi0_l(1)*psi0_l(3)
     bz0_l(4) = kb/bz0_l(1)                                &
          *(-kb/bz0_l(1)**2*(psi0_l(1)*psi0_l(2))**2     &
          + psi0_l(2)**2+psi0_l(1)*psi0_l(4))
     bz0_l(5) = kb/bz0_l(1)*(-kb/bz0_l(1)**2*       &
          psi0_l(1)**2*psi0_l(2)*psi0_l(3)                &
          + psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
     bz0_l(6) = kb/bz0_l(1)                               &
          *(-kb/bz0_l(1)**2*(psi0_l(1)*psi0_l(3))**2     &
          + psi0_l(3)**2+psi0_l(1)*psi0_l(6))
        
     kb = k**2*beta*(p0 - pi0*ipres)/p0
     pe0_l(1) = 0.5*kb*psi0_l(1)**2 + p0 - pi0*ipres
     pe0_l(2) = kb*psi0_l(1)*psi0_l(2)
     pe0_l(3) = kb*psi0_l(1)*psi0_l(3)
     pe0_l(4) = kb*(psi0_l(2)**2+psi0_l(1)*psi0_l(4))
     pe0_l(5) = kb*(psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
     pe0_l(6) = kb*(psi0_l(3)**2+psi0_l(1)*psi0_l(6))

     kb = k**2*beta
     p0_l(1) = 0.5*kb*psi0_l(1)**2 + p0
     p0_l(2) = kb*psi0_l(1)*psi0_l(2)
     p0_l(3) = kb*psi0_l(1)*psi0_l(3)
     p0_l(4) = kb*(psi0_l(2)**2+psi0_l(1)*psi0_l(4))
     p0_l(5) = kb*(psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
     p0_l(6) = kb*(psi0_l(3)**2+psi0_l(1)*psi0_l(6))
  endif

  call constant_field(den0_l, 1.)

end subroutine cylinder_equ

subroutine cylinder_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: befo, uzero, czero
  
  befo = eps*exp(-(x**2+z**2))
  uzero = befo*cos(z)
  czero = befo*sin(z)
  u1_l(1) = uzero
  u1_l(2) = - 2.*x*uzero
  u1_l(3) = - 2.*z*uzero - czero
  u1_l(4) = (-2.+4.*x**2)*uzero
  u1_l(5) = 4.*x*z*uzero + 2.*x*czero
  u1_l(6) = (-3.+4.*z**2)*uzero + 4.*z*czero
  
  vz1_l =  0.

  chi1_l(1) = czero
  chi1_l(2) = -2*x*czero
  chi1_l(3) = -2*z*czero + uzero
  chi1_l(4) = (-2.+4*x**2)*czero
  chi1_l(5) = 4.*x*z*czero - 2.*x*uzero
  chi1_l(6) = (-3.+4.*z**2)*czero - 4*z*uzero

  psi1_l =  0.
  bz1_l = 0.
  pe1_l = 0.
  p1_l = 0.
  den1_l = 0.

end subroutine cylinder_per

end module tilting_cylinder



!==============================================================================
! Taylor Reconnection (itaylor = 1)
!==============================================================================

module taylor_reconnection

contains


subroutine taylor_reconnection_init()
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     call assign_local_pointers(l)

     call taylor_reconnection_equ(x, z)
     call taylor_reconnection_per(x, z)
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
  call constant_field( pe0_l, p0-pi0*ipres)
  call constant_field(den0_l, 1)
  
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


!==============================================================================
! Force-Free Taylor State (itaylor = 2)
!==============================================================================
module force_free_state

contains

subroutine force_free_init()
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     call assign_local_pointers(l)

     call force_free_equ(x, z)
     call force_free_per(x, z)
  enddo

end subroutine force_free_init

subroutine force_free_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  
  real :: kx, kz, alx, alz, alam
  
  call getboundingboxsize(alx, alz)

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

  call constant_field( pe0_l, p0-pi0*ipres)
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

!===========================================================================
! GEM Reconnection (itaylor = 3)
!===========================================================================
module gem_reconnection

  real, private :: akx, akz

contains

subroutine gem_reconnection_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  akx = 2.*pi/alx
  akz = pi/alz

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     call assign_local_pointers(l)

     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     call gem_reconnection_equ(x, z)
     call gem_reconnection_per(x, z)
  enddo

end subroutine gem_reconnection_init

subroutine gem_reconnection_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: sech, pezero

  call constant_field(  u0_l, 0.)
  if(numvar.ge.2) call constant_field( vz0_l, 0.)
  if(numvar.ge.3) call constant_field(chi0_l, 0.)

  psi0_l(1) = 0.5*alog(cosh(2.*z))
  psi0_l(2) = 0.0
  psi0_l(3) = tanh(2.*z)  
  psi0_l(4) = 0.0
  psi0_l(5) = 0.0
  psi0_l(6) = 2.*sech(2.*z)**2

  ! if numvar = 2, then use Bz to satisfy force balance
  if(numvar.eq.2) then
     bz0_l(1) = sqrt(bzero**2 + sech(2.*z)**2)
     bz0_l(2) = 0.
     bz0_l(3) = -2.*tanh(2.*z)*sech(2.*z)**2/bz0_l(1)
     bz0_l(4) = 0.
     bz0_l(5) = 0.
     bz0_l(6) = (2.*sech(2.*z)**2/bz0_l(1))                      &
          *(bz0_l(3)*tanh(2.*z)/bz0_l(1)                          &
          - 2.*sech(2.*z)**2 + 4.*tanh(2.*z)**2)
  endif

  ! if numvar >= 3, then use pressure to satisfy force balance
  if(numvar.ge.3) then

     if(bzero.eq.0) then
        bz0_l(1) = sqrt(1.-2.*p0)*sech(2.*z)
        bz0_l(2) = 0.
        bz0_l(3) = -2.*bz0_l(1)*tanh(2.*z)
        bz0_l(4) = 0.
        bz0_l(5) = 0.
        bz0_l(6) = -2.*bz0_l(3)*tanh(2.*z) &
             -4.*bz0_l(1)*(1.-tanh(2.*z)**2)
     else
        bz0_l(1) = sqrt(bzero**2 + (1.-2.*p0)*sech(2.*z)**2)
        bz0_l(2) = 0.
        bz0_l(3) = -(1.-2.*p0)*2.*tanh(2.*z)*sech(2.*z)**2/bz0_l(1)
        bz0_l(4) = 0.
        bz0_l(5) = 0.
        bz0_l(6) = (1.-2.*p0)*(2.*sech(2.*z)**2/bz0_l(1))            &
             *(bz0_l(3)*tanh(2.*z)/bz0_l(1)                          &
             - 2.*sech(2.*z)**2 + 4.*tanh(2.*z)**2)
     endif

     pezero = p0 - pi0*ipres
    
     pe0_l(1) = pezero*(sech(2.*z)**2 + 0.2)
     pe0_l(2) = 0.
     pe0_l(3) = pezero*(-4.*sech(2.*z)**2*tanh(2.*z))
     pe0_l(4) = 0.
     pe0_l(5) = 0.
     pe0_l(6) = pezero*(-8.*sech(2.*z)**2*(sech(2.*z)**2-2.*tanh(2.*z)**2))
  endif

  if(ipres.eq.1) then
     p0_l(1) = p0*(sech(2.*z)**2 + 0.2)
     p0_l(2) = 0.
     p0_l(3) = p0*(-4.*sech(2.*z)**2*tanh(2.*z))
     p0_l(4) = 0.
     p0_l(5) = 0.
     p0_l(6) = p0*(-8.*sech(2.*z)**2*(sech(2.*z)**2-2.*tanh(2.*z)**2)) 
  else 
     call constant_field(p0_l, 1.)
  endif

  if(idens.eq.1) then
     den0_l(1) = sech(2*z)**2 + 0.2
     den0_l(2) = 0.
     den0_l(3) = -4.*sech(2*z)**2*tanh(2*z)
     den0_l(4) = 0.
     den0_l(5) = 0.
     den0_l(6) = -8.*sech(2*z)**2*(sech(2*z)**2-2.*tanh(2*z)**2)
  else
     call constant_field(den0_l, 1.)
  endif

end subroutine gem_reconnection_equ

subroutine gem_reconnection_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.
  
  psi1_l(1) =  eps*cos(akx*x)*cos(akz*z)
  psi1_l(2) = -eps*sin(akx*x)*cos(akz*z)*akx
  psi1_l(3) = -eps*cos(akx*x)*sin(akz*z)*akz
  psi1_l(4) = -eps*cos(akx*x)*cos(akz*z)*akx**2
  psi1_l(5) =  eps*sin(akx*x)*sin(akz*z)*akx*akz
  psi1_l(6) = -eps*cos(akx*x)*cos(akz*z)*akz**2

  bz1_l = 0.
  pe1_l = 0. 
  p1_l = 0.
  den1_l = 0.

end subroutine gem_reconnection_per

end module gem_reconnection


!==============================================================================
! Wave Propagation (itaylor = 4)
!==============================================================================
module wave_propagation

  implicit none
  real, private :: alx, alz, akx, akx2, omega

contains

subroutine wave_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  akx = 2.*pi/alx
  akx2 = akx**2
  omega = sqrt(akx2 + .5*db**2*akx2**2                                 &
       + db*akx2*sqrt(akx2 + .25*db**2*akx2**2))

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     call assign_local_pointers(l)

     call wave_equ(x, z)
     call wave_per(x, z)
  enddo

end subroutine wave_init

subroutine wave_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.
  
  psi0_l(1) = z
  psi0_l(2) = 0.0
  psi0_l(3) = 1.0
  psi0_l(4) = 0.0
  psi0_l(5) = 0.0
  psi0_l(6) = 0.0

  call constant_field( bz0_l, bzero)
  call constant_field( pe0_l, p0-pi0*ipres)
  call constant_field(  p0_l, p0)
  call constant_field(den0_l, 1.)

end subroutine wave_equ


subroutine wave_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: omega

  call plane_wave(u1_l, x, z,  akx, 0, eps*akx/omega, pi/2.)
  call plane_wave(vz1_l, x, z, akx, 0, eps*(akx2/omega**2-1), pi/2.)
  chi1_l = 0.
  
  call plane_wave(psi1_l, x, z, akx, 0, eps, pi/2.)
  call plane_wave( bz1_l, x, z, akx, 0, eps*(akx/omega-omega/akx), pi/2.)
  pe1_l = 0.
  p1_l = 0.

  den1_l = 0.
end subroutine wave_per

end module wave_propagation



!============================================================================
! Gravitational Instability Equilibrium (itor = 0, itaylor = 5)
!============================================================================
module grav

contains

subroutine grav_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call assign_local_pointers(l)

     call xyznod(l, coords)

     x = coords(1) + xzero - xmin
     z = coords(2) + zzero - zmin

     call grav_equ(x, z)
     call grav_per(x, z)
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

  
  fac1 = p0 - pi0*ipres
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
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: kx, kz, alx, alz

  call getboundingboxsize(alx, alz)
  kx = 2.*pi/alx
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


!==============================================================================
! Strauss Equilibrium (itor = 0, itaylor = 6)
!
! H.R. Strauss, Phys. Fluids 19(1):134 (1976)
!==============================================================================
module strauss

  implicit none

  real :: alx,alz

contains

subroutine strauss_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) - xmin - alx/2.
     z = coords(2) - zmin - alz/2.

     call assign_local_pointers(l)

     call strauss_equ(x, z)
     call strauss_per(x, z)
  enddo

end subroutine strauss_init

subroutine strauss_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: kx, kz, a0

    u0_l = 0.
  if(numvar.ge.2) vz0_l = 0.
  if(numvar.ge.3) chi0_l = 0.

  kx = pi/alx
  kz = pi/alz

  a0 = 2.*alz*alx*bzero/(pi*ln*q0)

  psi0_l(1) = -a0*cos(kx*x)*cos(kz*z)
  psi0_l(2) =  a0*sin(kx*x)*cos(kz*z)*kx
  psi0_l(3) =  a0*cos(kx*x)*sin(kz*z)*kz
  psi0_l(4) =  a0*cos(kx*x)*cos(kz*z)*kx**2
  psi0_l(5) = -a0*sin(kx*x)*sin(kz*z)*kx*kz
  psi0_l(6) =  a0*cos(kx*x)*cos(kz*z)*kz**2

  call constant_field(bz0_l, bzero)
  call constant_field(pe0_l, p0 - ipres*pi0)
  call constant_field(den0_l, 1.)
  call constant_field(p0_l, p0)

end subroutine strauss_equ


subroutine strauss_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: akx, akz

  akx = pi/alx
  akz = 2.*pi/alz

  psi1_l(1) = -eps*cos(akx*x)*sin(akz*z)
  psi1_l(2) =  eps*sin(akx*x)*sin(akz*z)*akx
  psi1_l(3) = -eps*cos(akx*x)*cos(akz*z)*akz
  psi1_l(4) =  eps*cos(akx*x)*sin(akz*z)*akx**2
  psi1_l(5) =  eps*sin(akx*x)*cos(akz*z)*akx*akz
  psi1_l(6) =  eps*cos(akx*x)*sin(akz*z)*akz**2
  u1_l = 0.
  
  bz1_l = 0.
  vz1_l = 0.
  pe1_l = 0.
  chi1_l = 0.
  den1_l = 0.
  p1_l = 0.

end subroutine strauss_per

end module strauss


!==============================================================================
! Circular magnetic field equilibrium (for parallel conduction tests)
!==============================================================================
module circular_field

  real :: alx, alz, xmin, zmin

contains

subroutine circular_field_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) - xmin - alx/2.
     z = coords(2) - zmin - alz/2.

     call assign_local_pointers(l)

     call circular_field_equ(x, z)
     call circular_field_per(x, z)
  enddo

end subroutine circular_field_init

subroutine circular_field_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: ss

  ss = min(alx,alz)/10.

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  psi0_l(1) = x**2 + z**2
  psi0_l(2) = 2.*x
  psi0_l(3) = 2.*z
  psi0_l(4) = 2.
  psi0_l(5) = 0.
  psi0_l(6) = 2.
!!$  psi0_l(1) = exp(-(x**2+z**2)/(2.*ss**2))
!!$  psi0_l(2) = -x*psi0_l(1)/ss**2
!!$  psi0_l(3) = -z*psi0_l(1)/ss**2
!!$  psi0_l(4) = ((x/ss)**2 - 1.)*psi0_l(1)/ss**2
!!$  psi0_l(5) =  x*z*psi0_l(1)/ss**4
!!$  psi0_l(6) = ((z/ss)**2 - 1.)*psi0_l(1)/ss**2


  call constant_field(bz0_l , bzero)
  call constant_field(p0_l  , p0)
  call constant_field(pe0_l , p0-pi0*idens)
  call constant_field(den0_l, 1.)

end subroutine circular_field_equ


subroutine circular_field_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: x0,z0

  x0 = 1.
  z0 = 0.

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.

  psi1_l = 0.
  bz1_l = 0.

  p1_l(1) = eps*exp(-((x-x0)**2+z**2)/(2.*ln**2))
  p1_l(2) = -(x-x0)*p1_l(1)/ln**2
  p1_l(3) = -(z-z0)*p1_l(1)/ln**2
  p1_l(4) = (((x-x0)/ln)**2 - 1.)*p1_l(1)/ln**2
  p1_l(5) =  (x-x0)*(z-z0)*p1_l(1)/ln**4
  p1_l(6) = (((z-z0)/ln)**2 - 1.)*p1_l(1)/ln**2

  ! for viscosity test..
  vz1_l = p1_l
  p1_l = 0.

  if(ipres.eq.1) then
     pe1_l = pefac*p1_l
  else
     pe1_l = p1_l
  endif

  den1_l = 0.

end subroutine circular_field_per

end module circular_field


!==============================================================================
! Magnetorotational Equilibrium (itor = 1, itaylor = 2)
!==============================================================================
module mri

  real, private :: kx, kz

contains

subroutine mri_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  kx = pi/alx
  kz = 2.*pi/alz

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) + xzero - xmin
     z = coords(2) + zzero - zmin - alz*.5

     call assign_local_pointers(l)

     call mri_equ(x, z)
     call mri_per(x, z)
  enddo

end subroutine mri_init

subroutine mri_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: nu, fac1, fac2

  nu = db*0.5*gyro*pi0/bzero
  fac1 = sqrt(gravr)
  fac2 = (9./128.)*nu**2/fac1

  call constant_field(  u0_l,0.)

  psi0_l(1) = bzero * x**2 / 2.
  psi0_l(2) = bzero * x
  psi0_l(3) = 0.
  psi0_l(4) = bzero
  psi0_l(5) = 0.
  psi0_l(6) = 0.

  if(numvar.ge.2) then
     bz0_l = 0.

     vz0_l(1)  = fac1*sqrt(x) + (3./8.)*nu - fac2/sqrt(x)
     vz0_l(2)  = 0.5*fac1/sqrt(x) + 0.5*fac2/(sqrt(x)**3)
     vz0_l(3)  = 0.
     vz0_l(4)  = -0.25*fac1/(sqrt(x)**3) - 0.75*fac2/(sqrt(x)**5)
     vz0_l(5) = 0.
     vz0_l(6) = 0.
  end if

  if(numvar.ge.3) then
     call constant_field( pe0_l,p0 - pi0*ipres)
     call constant_field(chi0_l,0.)
  end if

  call constant_field(den0_l, 1.)
  call constant_field(p0_l,p0)

end subroutine mri_equ


subroutine mri_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: fac1

  u1_l = 0.
  if(numvar.ge.2)  then
     fac1 = eps*sqrt(gravr*x)
     vz1_l(1) =  fac1*sin(kx*(x-xzero))*sin(kz*(z-zzero))
     vz1_l(2) =  fac1*cos(kx*(x-xzero))*sin(kz*(z-zzero))*kx
     vz1_l(3) =  fac1*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kz
     vz1_l(4) = -fac1*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kx**2
     vz1_l(5) =  fac1*cos(kx*(x-xzero))*cos(kz*(z-zzero))*kx*kz
     vz1_l(6) = -fac1*sin(kx*(x-xzero))*sin(kz*(z-zzero))*kz**2
  endif
  chi1_l = 0.
  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.

  den1_l = 0.
  p1_l = 0.

end subroutine mri_per

end module mri


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

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     x = coords(1) + xzero - xmin
     z = coords(2) + zzero - zmin - alz*.5

     call assign_local_pointers(l)

     call rotate_equ(x, z)
     call rotate_per(x, z)
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
  pe0_l(1) = p0 - ipres*pi0
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

  real :: fac1

  u1_l = 0.
  vz1_l = 0
  chi1_l = 0.
  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.
  den1_l = 0.
  p1_l = 0.

end subroutine rotate_per

end module rotate




!=====================================
subroutine initial_conditions()
  use basic

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

  implicit none

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
     end select
  endif

end subroutine initial_conditions

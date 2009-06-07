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
subroutine random_per(x,z,seed,fac)
  use basic
  use arrays

  implicit none

  real :: rand

  real, intent(in) :: x, z
  vectype, intent(in), dimension(6) :: fac
  integer, intent(in) :: seed
  integer :: i, j
  integer, parameter :: maxn = 10
  real :: alx, alz, kx, kz, xx, zz
  vectype, dimension(6) :: temp

  call getboundingboxsize(alx, alz)

  call srand(seed)

  temp = 0.

#ifdef _AIX
#define RAND_ARG
#else
#define RAND_ARG 0
#endif

  xx = x - xzero
  zz = z - zzero

  do i=1,maxn
     kx = pi*i/alx
     select case (icsym)
!
     case (0)   !  original option...no symmetry imposed
     do j=1, maxn
        kz = j*pi/alz
        call plane_wave2(temp,xx,zz,kx,kz,2.*eps*(rand(RAND_ARG)-.5),0.,0.)
        call add_product(psi1_l,fac,temp)
        call plane_wave2(temp,xx,zz,kx,kz,2.*eps*(rand(RAND_ARG)-.5),0.,0.)
        call add_product(u1_l,fac,temp)
     end do
!
     case (1)  !   make U odd symmetry about midplane:  perturb only U
     do j=1, maxn/2
        kz = 2.*j*pi/alz
        call plane_wave2(temp,xx,zz,kx,kz,2.*eps*(rand(RAND_ARG)-.5),0.,0.)
        call add_product(u1_l,fac,temp)
     end do
!
     case (2)  !   make U even  symmetry about midplane:  perturb only U
     do j=1, maxn/2
        kz = (2.*j-1)*pi/alz
        call plane_wave2(temp,xx,zz,kx,kz,2.*eps*(rand(RAND_ARG)-.5),0.,0.)
        call add_product(u1_l,fac,temp)
     end do
!
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

  implicit none

  integer :: inode, numnodes
  real :: x, z

  call numnod(numnodes)

  do inode=1, numnodes
     call nodcoord(inode, x, z)

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
subroutine den_eq
  use basic
  use arrays
  use diagnostics

  integer :: numnodes, inode
  real :: temp(6), k, kx
  real :: sech
  
  if(idenfunc.eq.0) return

  call numnod(numnodes)

  do inode=1, numnodes 

     call assign_local_pointers(inode)
     
     select case(idenfunc)
     case(1)      ! added 08/05/08 for stability benchmarking
        den0_l(1) = den0*.5* &
             (1. + &
             tanh((real(psi0_l(1))-(psibound+denoff*(psibound-psimin)))&
             /(dendelt*(psibound-psimin))))
        
     case(2)
        temp(1) = real((psi0_l(1)-psimin)/(psibound-psimin))
        temp(2:6) = real(psi0_l(2:6))/(psibound-psimin)

        k = 1./dendelt
        kx = k*(temp(1) - denoff)

        den0_l(1) = 1. + tanh(kx)
        den0_l(2) = k*sech(kx)**2 * temp(2)
        den0_l(3) = k*sech(kx)**2 * temp(3)
        den0_l(4) = k*sech(kx)**2 * temp(4) &
             -2.*k**2*sech(kx)**2*tanh(kx) * temp(2)**2
        den0_l(5) = k*sech(kx)**2 * temp(5) &
             -2.*k**2*sech(kx)**2*tanh(kx) * temp(2)*temp(3)
        den0_l(6) = k*sech(kx)**2 * temp(6) &
             -2.*k**2*sech(kx)**2*tanh(kx) * temp(3)**2

        den0_l = den0_l * 0.5*(den_edge-den0)
        den0_l(1) = den0_l(1) + den0
     end select

  end do

end subroutine den_eq
!==============================================================================
subroutine rmp_per
  use basic
  use arrays
  use coils

  implicit none

  logical :: is_boundary
  integer :: izone, izonedim, numnodes, l
  real :: normal(2), curv, x, z

  call load_field_from_coils('rmp_coil.dat', 'rmp_current.dat', &
       field,num_fields,psi_g)

  ! leave perterbation only on the boundary
  if(irmp.eq.2) then
     call numnod(numnodes)
     do l=1, numnodes
        call boundary_node(l,is_boundary,izone,izonedim,normal,curv,x,z)
        if(.not.is_boundary) then
           call assign_local_pointers(l)
           psi1_l = 0.
        endif
     end do
  endif
end subroutine rmp_per

!==============================================================================
! Tilting Cylinder (itaylor = 0)
!==============================================================================
module tilting_cylinder

  implicit none

  real, parameter :: k = 3.8317059702

contains

subroutine tilting_cylinder_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero 

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
     fp = 0.
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
     if(numvar.ge.2) then
        kb = k**2*(1.-beta)
        bz0_l(1) = sqrt(kb*psi0_l(1)**2+bzero**2)
        bz0_l(2) = kb/bz0_l(1)*psi0_l(1)*psi0_l(2)
        bz0_l(3) = kb/bz0_l(1)*psi0_l(1)*psi0_l(3)
        bz0_l(4) = kb/bz0_l(1)                              &
             *(-kb/bz0_l(1)**2*(psi0_l(1)*psi0_l(2))**2     &
             + psi0_l(2)**2+psi0_l(1)*psi0_l(4))
        bz0_l(5) = kb/bz0_l(1)*(-kb/bz0_l(1)**2*            &
             psi0_l(1)**2*psi0_l(2)*psi0_l(3)               &
             + psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
        bz0_l(6) = kb/bz0_l(1)                              &
             *(-kb/bz0_l(1)**2*(psi0_l(1)*psi0_l(3))**2     &
             + psi0_l(3)**2+psi0_l(1)*psi0_l(6))
     else 
       call constant_field(bz0_l, bzero)
     end if

     if(numvar.ge.3) then
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
     else
        call constant_field(pe0_l, p0-pi0*ipres)
        call constant_field( p0_l, p0)
     end if
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
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

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
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z, alx, alz

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x,z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

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
  real :: x, z, alx, alz

  call getboundingboxsize(alx, alz)

  akx = 2.*pi/alx
  akz = pi/alz

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     call assign_local_pointers(l)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

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


!=============================================================================
! Wave Propagation (itaylor = 4)
!=============================================================================
module wave_propagation

  implicit none
  real, private :: alx, alz, akx, akx2, omega
  real, private :: psiper, phiper, bzper, vzper, peper, chiper, nper, pper

contains

subroutine wave_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z
  real :: b2,a2
  real :: kp,km,t1,t2,t3
  real :: coef(4)
  real :: root(3)
  real :: error(3)
  real :: bi

  call getboundingboxsize(alx, alz)

  ! for itaylorw=3, set up a phi perturbation only
  if(iwave.eq.3) then
     phiper = eps
     vzper = 0.
     chiper = 0.
     psiper = 0.
     bzper = 0.
     nper = 0.
     pper = 0.
     peper = 0.
     goto 1
  endif

  akx = 2.*pi/alx      
  akx2 = akx**2
  b2 = bzero*bzero + bx0*bx0
  a2 = bx0**2/b2
  bi = 2.*pi0*gyro
  
  ! numvar=2 =================
  if(numvar.eq.2) then
     kp = akx * (b2 + bi*(3.*a2-1.)/4.)
     km = akx * (b2 - bi*(3.*a2-1.)/4.)

     if(iwave.eq.0) then ! fast wave
        omega = (akx/2.)*sqrt(a2/b2)*(sqrt(4.*b2**2+km**2)+kp)
        psiper = eps
        phiper = -eps*2.*b2 / (sqrt(4.*b2**2+km**2)+km)
        bzper = akx*psiper
        vzper = akx*phiper
        
     else ! slow wave
        omega = (akx/2.)*sqrt(a2/b2)*(sqrt(4.*b2**2+km**2)-kp)
        psiper = eps
        phiper = -eps*2.*b2 / (sqrt(4.*b2**2+km**2)-km)
        bzper = -akx*psiper
        vzper = -akx*phiper
     endif
     chiper = 0.
     nper = 0.
     peper = 0.
     pper = 0.
     
  ! numvar=3 =================
  elseif(numvar.ge.3) then
     
     coef(4) = 96.*b2**4
     coef(3) = -6.*akx2*b2**3*(16.*b2**2*(1.+a2+akx2*a2)            &
          + akx2*bi**2*(1.+6.*a2-3.*a2**2)                          &
          + 16.*gam*p0*b2)
     coef(2) = 6.*a2*akx2**2*b2**3                                  &
          *(b2*((4.*b2-bi*akx2*(1.+a2))**2                          &
          +4.*bi**2*akx2*(1.-a2)*(1.+akx2*a2))                      &
          + gam*p0*(32.*b2**2+16.*akx2*b2**2                        &
          +bi**2*akx2*(1.-3.*a2)**2))
     coef(1) = -6.*gam*p0*a2**2*akx2**3*b2**4                       &
          *(4.*b2+akx2*bi*(1.-3.*a2))**2

     if(myrank.eq.0 .and. iprint.eq.1) then
        write(*,*) "Coefs: ", coef(1), coef(2), coef(3), coef(4)
     endif

     call cubic_roots(coef, root, error)

     if(myrank.eq.0 .and. iprint.eq.1) then
        write(*,*) "Coefs: ", coef(1), coef(2), coef(3), coef(4)
        write(*,*) "Roots: ", root(1), root(2), root(3)
        write(*,*) "Error: ", error(1), error(2), error(3)
     endif

     select case(iwave)
     case(1)
        omega=sqrt(root(1))
     case(2)
        omega=sqrt(root(2))
     case default 
        omega=sqrt(root(3))
     end select

     t1 = akx2*bi/(4.*b2)
     t2 = (akx2/omega**2)*                                          &
          (bzero**2/(1.-gam*p0*akx2/(omega**2)))
     t3 = bx0 * akx / omega

     psiper = eps
     phiper = psiper*(akx2*t3*(1.+(t3**2/(1.-t2-t3**2)))-(1.-t2)/t3)&
          / (1. - t2 + t1*t2*(1.+3.*a2)                             &
          + t1*t3**2*(2.*t2-(1.-3.*a2))/(1-t2-t3**2))
     vzper = phiper*t1*t3*(2.*t2-(1.-3.*a2))/(1-t2-t3**2)           &
          - psiper*akx2*t3**2/(1.-t2-t3**2)
     bzper = (-phiper*t1*t2*(1.+3.*a2) - vzper*t3 + psiper*akx2*t3) &
          / (1.-t2)
     chiper = bzero / (1.-gam*p0*akx2/(omega**2)) / omega           &
          * ((1.+3.*a2)*t1*phiper - bzper)
     peper = -chiper*gam*(p0-pi0)*akx2 / omega
     nper = -chiper*akx2/omega
     pper = -chiper*gam*p0*akx2 / omega
  endif

  if(myrank.eq.0) then
     print *, "Wave angular frequency: ", omega
     print *, "Wave phase velocity: ", omega/akx
     print *, "Wave transit time: ", 2.*pi/omega
  end if

1 continue

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

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
  
  psi0_l(1) = z*bx0
  psi0_l(2) = 0.0
  psi0_l(3) = bx0
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

  call plane_wave(u1_l, x, z,  akx, 0, phiper, pi/2.)
  call plane_wave(vz1_l, x, z, akx, 0, vzper, pi/2.)
  call plane_wave(chi1_l, x, z, akx, 0, chiper, pi)
  
  call plane_wave(psi1_l, x, z, akx, 0, psiper, pi/2.)
  call plane_wave( bz1_l, x, z, akx, 0, bzper, pi/2.)

  if(ipres.eq.1) then 
     call plane_wave(pe1_l, x, z, akx, 0, peper, pi/2.)
     call plane_wave( p1_l, x, z, akx, 0,  pper, pi/2.)
  else
     call plane_wave(pe1_l, x, z, akx, 0, pper, pi/2.)
     call plane_wave( p1_l, x, z, akx, 0, pper, pi/2.)    
  endif

  if(idens.eq.1) then
     call plane_wave(den1_l, x, z, akx, 0, nper, pi/2.)    
  endif

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
  real :: x, z

  call numnod(numnodes)
  do l=1, numnodes
     call assign_local_pointers(l)

     call nodcoord(l, x, z)

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
  real :: x, z

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     x = x - alx/2.
     z = z - alz/2.

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

  real :: alx, alz

contains

subroutine circular_field_init()
  use basic
  use arrays

  implicit none

  integer :: l, numnodes
  real :: x, z

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     x = x - alx/2.
     z = z - alz/2.

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
  real :: ss

  ss = min(alx,alz)/4.

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.

  psi1_l = 0.
  bz1_l = 0.
  p1_l = 0.

  p1_l(1) = eps*exp(-((x-x0)**2+z**2)/(2.*ln**2))
  p1_l(2) = -(x-x0)*p1_l(1)/ln**2
  p1_l(3) = -(z-z0)*p1_l(1)/ln**2
  p1_l(4) = (((x-x0)/ln)**2 - 1.)*p1_l(1)/ln**2
  p1_l(5) =  (x-x0)*(z-z0)*p1_l(1)/ln**4
  p1_l(6) = (((z-z0)/ln)**2 - 1.)*p1_l(1)/ln**2

  ! for viscosity test..
!  vz1_l = p1_l
!  p1_l = 0.

  ! for parallel viscosity test...

  u1_l(1) = eps*exp(-(x**2+z**2)/(2.*ss**2))
  u1_l(2) = -x*u1_l(1)/ss**2
  u1_l(3) = -z*u1_l(1)/ss**2
  u1_l(4) = ((x/ss)**2 - 1.)*u1_l(1)/ss**2
  u1_l(5) =  x*z*u1_l(1)/ss**4
  u1_l(6) = ((z/ss)**2 - 1.)*u1_l(1)/ss**2
  vz1_l(1) = vzero*exp(-(x**2+z**2)/(2.*ss**2))
  vz1_l(2) = -x*vz1_l(1)/ss**2
  vz1_l(3) = -z*vz1_l(1)/ss**2
  vz1_l(4) = ((x/ss)**2 - 1.)*vz1_l(1)/ss**2
  vz1_l(5) =  x*z*psi0_l(1)/ss**4
  vz1_l(6) = ((z/ss)**2 - 1.)*vz1_l(1)/ss**2
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
  real :: x, z, alx, alz

  call getboundingboxsize(alx, alz)

  kx = pi/alx
  kz = 2.*pi/alz

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     z = z - alz*.5

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
  real :: x, z, alx, alz

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     z = z - alz*.5

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


!==============================================================================
! Eqdsk_eq
! ~~~~~~~~
!
! Loads eqdsk equilibrium
!==============================================================================
module eqdsk_eq

contains

subroutine eqdsk_init()
  use basic
  use arrays
  use eqdsk
  use gradshafranov
  use newvar_mod
  use sparse
  use diagnostics

  implicit none

  integer :: l, numnodes
  real :: x, z
  real, parameter :: amu0 = pi*4.e-7
  real, allocatable :: flux(:)

  print *, "eqdsk_init called"

  call load_eqdsk
  press = press*amu0
  pprime = pprime*amu0
  current = current*amu0

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     call assign_local_pointers(l)

     call eqdsk_equ(x, z)
     call eqdsk_per(x, z)
  enddo

  tcuro = current
  xmag = rmaxis
  zmag = zmaxis
  rzero = rmaxis
  bzero = fpol(nw)/rzero

  if(igs.gt.0) then
     fieldi = field0
     separatrix_top = abs(znull)
     separatrix_bottom = -separatrix_top

     if(iread_eqdsk.eq.2) then
        call default_profiles
     else
        allocate(flux(nw))
        do l=1,nw
           flux(l) = l*(sibry-simag)/(nw-1.)
        end do
        call create_profile(nw,press,pprime,fpol,ffprim,flux)
        deallocate(flux)
     end if

     ! initial guess
!!$     call newvar(mass_matrix_lhs_dc,jphi,field0,psi_g,num_fields, &
!!$          gs_matrix_rhs_dc,NV_DCBOUND)
     call deltafun(xmag,zmag,tcuro,jphi,1,1)

     call gradshafranov_solve
     call gradshafranov_per
  else
     psibound = sibry
     psimin = simag
  endif

  call unload_eqdsk

end subroutine eqdsk_init

subroutine eqdsk_equ(x, z)
  use basic
  use arrays

  use eqdsk

  implicit none

  real, intent(in) :: x, z
  real :: p, q
  integer :: i, j, n, m

  real :: dx, dz, temp, dpsi
  real, dimension(4,4) :: a
  real, dimension(4) :: b, c
  
  dx = rdim/(nw - 1.)
  dz = zdim/(nh - 1.)
  p = (x-rleft)/dx + 1.
  q = (z-zmid)/dz + nh/2. + 1.
  i = p
  j = q

  call bicubic_interpolation_coeffs(psirz,nw,nh,i,j,a)

  psi0_l = 0.
  do n=1, 4
     do m=1, 4
        temp = a(n,m)
        if(n.gt.1) temp = temp*(p-i)**(n-1)
        if(m.gt.1) temp = temp*(q-j)**(m-1)
        psi0_l(1) = psi0_l(1) + temp

        temp = a(n,m)*(n-1)
        if(n.gt.2) temp = temp*(p-i)**(n-2)
        if(m.gt.1) temp = temp*(q-j)**(m-1)
        psi0_l(2) = psi0_l(2) + temp/dx

        temp = a(n,m)*(m-1)
        if(n.gt.1) temp = temp*(p-i)**(n-1)
        if(m.gt.2) temp = temp*(q-j)**(m-2)
        psi0_l(3) = psi0_l(3) + temp/dz

        temp = a(n,m)*(n-1)*(n-2)
        if(n.gt.3) temp = temp*(p-i)**(n-3)
        if(m.gt.1) temp = temp*(q-j)**(m-1)
        psi0_l(4) = psi0_l(4) + temp/dx**2

        temp = a(n,m)*(n-1)*(m-1)
        if(n.gt.2) temp = temp*(p-i)**(n-2)
        if(m.gt.2) temp = temp*(q-j)**(m-2)
        psi0_l(5) = psi0_l(5) + temp/(dx*dz)

        temp = a(n,m)*(m-1)*(m-2)
        if(n.gt.1) temp = temp*(p-i)**(n-1)
        if(m.gt.3) temp = temp*(q-j)**(m-3)
        psi0_l(6) = psi0_l(6) + temp/dz**2
     end do
  end do

  ! calculation of p
  ! ~~~~~~~~~~~~~~~~
  dpsi = (sibry - simag)/(nw - 1.)
  p = (psi0_l(1) - simag)/dpsi + 1.
  i = p

  if(i.gt.nw) then
     call constant_field(p0_l, press(nw))
     call constant_field(bz0_l, fpol(nw))
  else

     ! use press and fpol to calculate values of p and I
     call cubic_interpolation_coeffs(press,nw,i,b)
     call cubic_interpolation_coeffs(fpol,nw,i,c)

     do n=1,4
        temp = b(n)
        if(n.gt.1) temp = temp*(p-i)**(n-1)
        p0_l(1) = p0_l(1) + temp

        temp = c(n)
        if(n.gt.1) temp = temp*(p-i)**(n-1)
        bz0_l(1) = bz0_l(1) + temp
     end do


     ! use pprime and ffprime to calculate derivatives
     call cubic_interpolation_coeffs(pprime,nw,i,b)
     call cubic_interpolation_coeffs(ffprim,nw,i,c)

     do n=1,4
        temp = b(n)
        if(n.gt.1) temp = temp*(p-i)**(n-1)
        p0_l(2) = p0_l(2) + temp*psi0_l(2)
        p0_l(3) = p0_l(3) + temp*psi0_l(3)
        
        temp = b(n)*(n-1)
        if(n.gt.2) temp = temp*(p-i)**(n-2)
        p0_l(4) = p0_l(4) + temp*psi0_l(4)
        p0_l(5) = p0_l(5) + temp*psi0_l(5)
        p0_l(6) = p0_l(6) + temp*psi0_l(6)

        temp = c(n)/bz0_l(1)
        if(n.gt.1) temp = temp*(p-i)**(n-1)
        bz0_l(2) = bz0_l(2) + temp*psi0_l(2)
        bz0_l(3) = bz0_l(3) + temp*psi0_l(3)
        
        temp = c(n)*(n-1)/bz0_l(1)
        if(n.gt.2) temp = temp*(p-i)**(n-2)
        bz0_l(4) = bz0_l(4) + temp*psi0_l(4)
        bz0_l(5) = bz0_l(5) + temp*psi0_l(5)
        bz0_l(6) = bz0_l(6) + temp*psi0_l(6)
     end do
  endif

  p0_l = p0_l + pedge

  ! Set electron pressure and density
  pe0_l = (1. - ipres*pi0/p0)*p0_l

  if(expn.eq.0) then
     call constant_field(den0_l,1.)
  else
     den0_l = (p0_l/p0)**expn
  end if

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

end subroutine eqdsk_equ


subroutine eqdsk_per(x, z)
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

end subroutine eqdsk_per
  
end module eqdsk_eq


!==============================================================================
! Dskbal_eq
! ~~~~~~~~~
!
! Loads dskbal equilibrium
!==============================================================================
module dskbal_eq

  implicit none

contains

subroutine dskbal_init()
  use basic
  use arrays
  use dskbal
  use gradshafranov
  use newvar_mod
  use sparse
  use diagnostics

  implicit none

  integer :: i, inode, numnodes
  real, parameter :: amu0 = pi*4.e-7
  real :: dp, a(4), minden, minte, minti

  print *, "dskbal_init called"

  call load_dskbal
  p_bal = p_bal*amu0
  pprime_bal = pprime_bal*amu0

  ! find "vacuum" values
  do i=1,npsi_bal
     if(ne_bal(i) .le. 0.) exit
     minden = ne_bal(i)
     minte = te_bal(i)
     minti = ti_bal(i)
  end do

  ! replace zeroed-out values with "vacuum" values
  do i=1,npsi_bal
     if(ne_bal(i) .le. 0.) ne_bal(i) = minden
     if(ti_bal(i) .le. 0.) ti_bal(i) = minti
     if(te_bal(i) .le. 0.) te_bal(i) = minte
     if(p_bal(i) .le. 0.) p_bal(i) = minden*(minte + minti)*1.6022e-12*amu0/10.
  end do

  xmag = x_bal(1,1)
  zmag = z_bal(1,1)
  rzero = xmag
  bzero = f_bal(npsi_bal)/rzero

  ifixedb = 1

  fieldi = field0
  separatrix_top = abs(znull)
  separatrix_bottom = -separatrix_top

  if(iread_dskbal.eq.2) then
     call default_profiles
  else
     call create_profile(npsi_bal,p_bal,pprime_bal,f_bal,ffprime_bal,psi_bal)
  end if
  
  ! initial plasma current filament
  call deltafun(xmag,zmag,tcuro,jphi,1,1)

  call gradshafranov_solve
  call gradshafranov_per  

  ! set density profile
  call numnod(numnodes)
  do inode=1, numnodes
     call assign_local_pointers(inode)
     do i=1, npsi_bal-1
        if((psi_bal(i+1)-psi_bal(1))/(psi_bal(npsi_bal)-psi_bal(1)) &
             .gt. (real(psi0_l(1))-psimin)/(psilim-psimin)) exit
     end do

     call cubic_interpolation_coeffs(ne_bal,npsi_bal,i,a)

     dp = ((real(psi0_l(1))-psimin) &
          *(psi_bal(npsi_bal) - psi_bal(1))/(psilim - psimin) &
          -(psi_bal(i)-psi_bal(1))) / (psi_bal(i+1) - psi_bal(i))

     den0_l(1) = a(1) + a(2)*dp + a(3)*dp**2 + a(4)*dp**3
     den0_l(2) = (a(2) + 2.*a(3)*dp + 3.*a(4)*dp**2)*psi0_l(2)
     den0_l(3) = (a(2) + 2.*a(3)*dp + 3.*a(4)*dp**2)*psi0_l(3)
     den0_l(4) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(4)
     den0_l(5) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(5)
     den0_l(6) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(6)
     den0_l = den0_l / n0_norm
  enddo

  call unload_dskbal

end subroutine dskbal_init
  
end module dskbal_eq


!==============================================================================
! Jsolver_eq
! ~~~~~~~~~~
!
! Loads eqdsk equilibrium
!==============================================================================
module jsolver_eq

contains

subroutine jsolver_init()
  use basic
  use arrays
  use eqdsk
  use gradshafranov
  use newvar_mod
  use sparse
  use coils
  use diagnostics
  use jsolver

  implicit none

  integer :: l, numnodes
  real :: x, z
  real, parameter :: amu0 = pi*4.e-7
  real, allocatable :: ffprime(:)

  print *, "jsolver_init called"

  call load_jsolver

  call numnod(numnodes)
  do l=1, numnodes
     call nodcoord(l, x, z)

     call assign_local_pointers(l)

     call jsolver_equ(x, z)
     call jsolver_per(x, z)
  enddo

  xmag = x_jsv(1,1)
  zmag = z_jsv(1,1)
  rzero = xzero_jsv
  bzero = gxx_jsv(1)

  ! remove factor of 1/xzero_jsv in definition of g
  gxx_jsv = gxx_jsv*xzero_jsv
  gpx_jsv = gpx_jsv*xzero_jsv

  ! shift gxx and p to be defined at psi_j locations
  ! instead of at psi_{j+1/2} locations
  do l=npsi_jsv, 2, -1
     gxx_jsv(l) = 0.5*(gxx_jsv(l) + gxx_jsv(l-1))
     p_jsv(l) = 0.5*(p_jsv(l) + p_jsv(l-1))
  end do

  if(igs.gt.0) then
     fieldi = field0
     separatrix_top = abs(znull)
     separatrix_bottom = -separatrix_top

     if(iread_jsolver.eq.2) then
        call default_profiles
     else
        allocate(ffprime(npsi_jsv))
        ffprime = gpx_jsv*gxx_jsv
        call create_profile(npsi_jsv,p_jsv,ppxx_jsv,gxx_jsv,ffprime,psival_jsv)
        deallocate(ffprime)
     end if

     call deltafun(xmag,zmag,tcuro,jphi,1,1)

     call gradshafranov_solve
     call gradshafranov_per
  else
     psibound = psival_jsv(npsi_jsv)
     psimin = psival_jsv(1)
  endif

  call unload_jsolver

end subroutine jsolver_init

subroutine jsolver_equ(x, z)
  use basic
  use arrays

  use jsolver

  implicit none

  real, intent(in) :: x, z

  integer :: i, j, n
  real :: i0, j0, d, dmin, temp
  real, dimension(4) :: a, b, c

  ! determine jsolver index (i0,j0) of this node
  dmin = (x - x_jsv(2,1))**2 + (z - z_jsv(2,1))**2
  do i=2, nthe
     do j=1, npsi_jsv
        d = (x - x_jsv(i,j))**2 + (z - z_jsv(i,j))**2
        if(d.lt.dmin) then
           dmin = d
           i0 = i
           j0 = j
        end if
     end do
  end do
  
  j = j0
  call cubic_interpolation_coeffs(psival_jsv,npsi_jsv,j,a)
  call cubic_interpolation_coeffs(gxx_jsv,npsi_jsv,j,b)
  call cubic_interpolation_coeffs(p_jsv,npsi_jsv,j,c)
  
  do n=1,4
     temp = a(n)
     if(n.gt.1) temp = temp*(j-j0)**(n-1)
     psi0_l(1) = psi0_l(1) + temp

     temp = b(n)
     if(n.gt.1) temp = temp*(j-j0)**(n-1)
     bz0_l(1) = bz0_l(1) + temp

     temp = c(n)
     if(n.gt.1) temp = temp*(j-j0)**(n-1)
     p0_l(1) = p0_l(1) + temp
  end do
  
  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

end subroutine jsolver_equ


subroutine jsolver_per(x, z)
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

end subroutine jsolver_per
  
end module jsolver_eq



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
  use eqdsk_eq
  use dskbal_eq
  use jsolver_eq

  implicit none

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
        end select
     endif
  end if
     
  call den_eq()
  if(irmp.ge.1) call rmp_per()

end subroutine initial_conditions

! ===============================================================
subroutine cubic_roots(coef, root, error)

!     calculate the roots of a 3rd degree polynomial
!     coef(4)*x**3 + coef(3)*x**2 + coef(2)*x + coef(1) = 0

  use basic
  implicit none
  real, intent(out), dimension(3) :: root, error
  real, intent(inout), dimension(4) :: coef
  real qqs, rrs, ths, fun, dfun
  integer r, ll

  !     normalize the coefficients to the 3rd degree coefficient
  coef(1) = coef(1)/coef(4)
  coef(2) = coef(2)/coef(4)
  coef(3) = coef(3)/coef(4)
  coef(4) = 1.0
       
  !     solve using method in Numerical Recipes
  qqs = ((coef(3)**2 - 3.*coef(2))/9.)
  rrs = (2.*coef(3)**3 - 9.*coef(3)*coef(2) + 27.*coef(1))/54.
  ths = acos(rrs/sqrt(qqs**3))

  root(1) = -2.*sqrt(qqs)*cos(ths/3.) - coef(3)/3.
  root(2) = -2.*sqrt(qqs)*cos((ths+6.283185307)/3.) - coef(3)/3.
  root(3) = -2.*sqrt(qqs)*cos((ths-6.283185307)/3.) - coef(3)/3.
  
  !........refine with Newton's method
  do r=1,3
     do ll=1,3
        fun  = root(r)**3 + coef(3)*root(r)**2 + coef(2)*root(r)    &
             + coef(1)
        dfun = 3.*root(r)**2 + 2.*coef(3)*root(r) + coef(2)
        root(r) = root(r) - fun/dfun
     enddo
     error(r) =  root(r)**3 + coef(3)*root(r)**2 + coef(2)*root(r)  &
          + coef(1)
  enddo
!      
end subroutine cubic_roots

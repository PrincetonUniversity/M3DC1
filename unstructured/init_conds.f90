!=====================================
subroutine initial_conditions()
  use basic

  implicit none

  integer :: ibegin, iendplusone, l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  do l=1, numnodes
     call xyznod(l, coords)

     if(itor.eq.0) then

        x = coords(1) - xmin - alx*.5
        z = coords(2) - zmin - alz*.5

        select case(itaylor)
        case(0)
           call cylinder_equ(x, z, l)
           call cylinder_per(x, z, l)
        case(1)
           call taylor_reconnection_equ(x, z, l)
           call taylor_reconnection_per(x, z, l)
        case(2)
           call force_free_equ(x, z, l)
           call force_free_per(x, z, l)
        case(3)
           call gem_reconnection_equ(x, z, l)
           call gem_reconnection_per(x, z, l)
        case(4)
           call wave_equ(x, z, l)
           call wave_per(x, z, l)
        end select
     else
        x = coords(1) + xzero - xmin
        z = coords(2) + zzero - zmin - alz*.5

        select case(itaylor)
        case(0)
           call solovev_equ(x, z, l)
           call solovev_per(x, z, l)
        end select
     endif
  enddo

end subroutine initial_conditions

!==============================
subroutine constant_field(outarr, val)
  implicit none

  real, dimension(6), intent(out) :: outarr
  real, intent(in) :: val

  outarr(1) = val
  outarr(2:6) = 0.

end subroutine constant_field
!==============================
subroutine static_equ(ibegin)
  use basic
  use arrays

  implicit none

  integer, intent(in) :: ibegin

  call constant_field(vel0(ibegin   :ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(vel0(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(vel0(ibegin+12:ibegin+17), 0.)
end subroutine static_equ
!==============================
subroutine plane_wave(x, z, outarr, kx, kz, amp, phase)
  implicit none

  real, intent(in) :: x, z
  real, dimension(6), intent(out) :: outarr
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



!==============================================================================
! Tilting Cylinder (itaylor = 0)
!==============================================================================
subroutine cylinder_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real, parameter :: k = 3.8317059702
  real :: rr, ri, arg, befo, ff, fp, fpp, j0, j1, kb
  real :: d1, d2, d3, d4, d5, d6
  real :: s17aef, s17aff
  integer :: ifail1, ifail2, ifail3

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call static_equ(ibegin)

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
   
  phi0(ibegin) = ff* z
  phi0(ibegin+1) = fp*z*x *ri
  phi0(ibegin+2) = fp*z**2*ri + ff
  phi0(ibegin+3) = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)
  phi0(ibegin+4) = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
  phi0(ibegin+5) = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)

  if(rr.gt.1) then
     if(numvar.ge.2) call constant_field(phi0(ibegin+6 :ibegin+11), bzero   )
     if(numvar.ge.3) call constant_field(phi0(ibegin+12:ibegin+17), p0-pi0*ipres)
  else 
     if(numvar.ge.2) then
        kb = k**2*(1.-beta)
        phi0(ibegin+6) = sqrt(kb*phi0(ibegin)**2+bzero**2)
        phi0(ibegin+7) = kb/phi0(ibegin+6)*phi0(ibegin)*phi0(ibegin+1)
        phi0(ibegin+8) = kb/phi0(ibegin+6)*phi0(ibegin)*phi0(ibegin+2)
        phi0(ibegin+9) = kb/phi0(ibegin+6)                                &
             *(-kb/phi0(ibegin+6)**2*(phi0(ibegin)*phi0(ibegin+1))**2     &
             + phi0(ibegin+1)**2+phi0(ibegin)*phi0(ibegin+3))
        phi0(ibegin+10) = kb/phi0(ibegin+6)*(-kb/phi0(ibegin+6)**2*       &
             phi0(ibegin)**2*phi0(ibegin+1)*phi0(ibegin+2)                &
             + phi0(ibegin+1)*phi0(ibegin+2)+phi0(ibegin)*phi0(ibegin+4))
        phi0(ibegin+11) = kb/phi0(ibegin+6)                               &
             *(-kb/phi0(ibegin+6)**2*(phi0(ibegin)*phi0(ibegin+2))**2     &
             + phi0(ibegin+2)**2+phi0(ibegin)*phi0(ibegin+5))
        
        if(numvar.ge.3) then
           kb = k**2*beta
           phi0(ibegin+12) = 0.5*kb*phi0(ibegin)**2 + p0 - pi0*ipres
           phi0(ibegin+13) = kb*phi0(ibegin)*phi0(ibegin+1)
           phi0(ibegin+14) = kb*phi0(ibegin)*phi0(ibegin+2)
           phi0(ibegin+15) = kb*(phi0(ibegin+1)**2+                      &
                phi0(ibegin)*phi0(ibegin+3))
           phi0(ibegin+16) = kb*(phi0(ibegin+1)*phi0(ibegin+2)+          &
                phi0(ibegin)*phi0(ibegin+4))
           phi0(ibegin+17) = kb*(phi0(ibegin+2)**2+                      &
                phi0(ibegin)*phi0(ibegin+5))
        endif
     endif
  endif

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)

     d1 = ff* z
     d3 = fp*z**2*ri + ff
     d2 = fp*z*x *ri
     d6 = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)
     d5 = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
     d4 = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)

     kb = k**2*beta
     den0(ibegin) = 0.5*kb*d1**2 + 1.
     den0(ibegin+1) = kb*d1*d2
     den0(ibegin+2) = kb*d1*d3
     den0(ibegin+3) = kb*(d2**2+d1*d4)
     den0(ibegin+4) = kb*(d2*d3+d1*d5)
     den0(ibegin+5) = kb*(d3**2+d1*d6)
  endif

end subroutine cylinder_equ

subroutine cylinder_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  
  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call constant_field(vel(ibegin   :ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(vel(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(vel(ibegin+12:ibegin+17), 0.)
  
  call constant_field(phi(ibegin   :ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(phi(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(phi(ibegin+12:ibegin+17), 0.)

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den(ibegin:ibegin+5 ), 0.)
  endif

end subroutine cylinder_per



!==============================================================================
! Taylor Reconnection (itaylor = 1)
!==============================================================================

subroutine taylor_reconnection_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call static_equ(ibegin)

  phi0(ibegin  ) = -0.5*z**2
  phi0(ibegin+1) =  0.
  phi0(ibegin+2) = -z
  phi0(ibegin+3) =  0.
  phi0(ibegin+4) =  0.
  phi0(ibegin+5) = -1.

  if(numvar.ge.2) call constant_field(phi0(ibegin+6 :ibegin+11), bzero   )
  if(numvar.ge.3) call constant_field(phi0(ibegin+12:ibegin+17), p0-pi0*ipres)
  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den0 (ibegin:ibegin+5), 1)
  endif
     
end subroutine taylor_reconnection_equ

subroutine taylor_reconnection_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call constant_field(vel(ibegin   :ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(vel(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(vel(ibegin+12:ibegin+17), 0.)
  
  call constant_field(phi(ibegin   :ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(phi(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(phi(ibegin+12:ibegin+17), 0.)  

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den(ibegin:ibegin+5 ), 0.)
  end if

end subroutine taylor_reconnection_per



!==============================================================================
! Force-Free Taylor State (itaylor = 2)
!==============================================================================
subroutine force_free_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode
  
  integer :: ibegin, iendplusone
  real :: kx, kz, alx, alz, alam
  
  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call getboundingboxsize(alx, alz)

  kx = pi/alx
  kz = pi/alz

  call static_equ(ibegin)

  alam = sqrt(kx**2+kz**2)

  phi0(ibegin) = sin(kx*x)*sin(kz*z)
  phi0(ibegin+1) = kx*cos(kx*x)*sin(kz*z)
  phi0(ibegin+2) = kz*sin(kx*x)*cos(kz*z)
  phi0(ibegin+3) = -kx**2*sin(kx*x)*sin(kz*z)
  phi0(ibegin+4) =  kx*kz*cos(kx*x)*cos(kz*z)
  phi0(ibegin+5) = -kz**2*sin(kz*x)*sin(kz*z)

  if(numvar.ge.2) then
     phi0(ibegin+6) = alam*phi0(ibegin)
     phi0(ibegin+7) = alam*phi0(ibegin+1)
     phi0(ibegin+8) = alam*phi0(ibegin+2)
     phi0(ibegin+9) = alam*phi0(ibegin+3)
     phi0(ibegin+10) = alam*phi0(ibegin+4)
     phi0(ibegin+11) = alam*phi0(ibegin+5)
  endif

  if(numvar.ge.3) call constant_field(phi0(ibegin+12:ibegin+17), p0-pi0*ipres)
  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den0 (ibegin:ibegin+5), 1.)
  endif

end subroutine force_free_equ

subroutine force_free_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call constant_field(vel(ibegin:ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(vel(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(vel(ibegin+12:ibegin+17), 0.)
  
  call constant_field(phi(ibegin:ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(phi(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(phi(ibegin+12:ibegin+17), 0.)  

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den(ibegin:ibegin+5 ), 0.)
  end if

end subroutine force_free_per

!==============================================================================
! GEM Reconnection (itaylor = 3)
!==============================================================================
subroutine gem_reconnection_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: sech, pezero

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call static_equ(ibegin)

  phi0(ibegin)   = 0.5*alog(cosh(2.*z))
  phi0(ibegin+1) = 0.0
  phi0(ibegin+2) = tanh(2.*z)  
  phi0(ibegin+3) = 0.0
  phi0(ibegin+4) = 0.0
  phi0(ibegin+5) = 2.*sech(2.*z)**2

  ! if numvar = 2, then use Bz to satisfy force balance
  if(numvar.eq.2) then
     phi0(ibegin+6) =  sqrt(bzero**2 + sech(2.*z)**2)
     phi0(ibegin+7) =  0.
     phi0(ibegin+8) =  -2.*tanh(2.*z)*sech(2.*z)**2/phi0(ibegin+6)
     phi0(ibegin+9) = 0.
     phi0(ibegin+10) = 0.
     phi0(ibegin+11) = (2.*sech(2.*z)**2/phi0(ibegin+6))                      &
          *(phi0(ibegin+8)*tanh(2.*z)/phi0(ibegin+6)                          &
          - 2.*sech(2.*z)**2 + 4.*tanh(2.*z)**2)
  endif

  ! if numvar >= 3, then use pressure to satisfy force balance
  if(numvar.ge.3) then
     call constant_field(phi0(ibegin+6 :ibegin+11), bzero)

     pezero = p0-ipres*pi0
    
     phi0(ibegin+12) = pezero*(sech(2.*z)**2 + 0.2)
     phi0(ibegin+13) = 0.
     phi0(ibegin+14) = pezero*(-4.*sech(2.*z)**2*tanh(2.*z))
     phi0(ibegin+15) = 0.
     phi0(ibegin+16) = 0.
     phi0(ibegin+17) = pezero*(-8.*sech(2.*z)**2*(sech(2.*z)**2-2.*tanh(2.*z)**2))
  endif

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     den0(ibegin) = sech(2*z)**2 + 0.2
     den0(ibegin+1) = 0.
     den0(ibegin+2) = -4.*sech(2*z)**2*tanh(2*z)
     den0(ibegin+3) = 0.
     den0(ibegin+4) = 0.
     den0(ibegin+5) = -8.*sech(2*z)**2*(sech(2*z)**2-2.*tanh(2*z)**2)
  endif
end subroutine gem_reconnection_equ

subroutine gem_reconnection_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: akx, akz, alx, alz

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call getboundingboxsize(alx, alz)

  akx = 2.*pi/alx
  akz = pi/alz

  call constant_field(vel(ibegin:ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(vel(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(vel(ibegin+12:ibegin+17), 0.)
  
  phi(ibegin)   =  eps*cos(akx*x)*cos(akz*z)
  phi(ibegin+1) = -eps*sin(akx*x)*cos(akz*z)*akx
  phi(ibegin+2) = -eps*cos(akx*x)*sin(akz*z)*akz
  phi(ibegin+3) = -eps*cos(akx*x)*cos(akz*z)*akx**2
  phi(ibegin+4) =  eps*sin(akx*x)*sin(akz*z)*akx*akz
  phi(ibegin+5) = -eps*cos(akx*x)*cos(akz*z)*akz**2

  if(numvar.ge.2)  call constant_field(phi(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(phi(ibegin+12:ibegin+17), 0.)  

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den(ibegin:ibegin+5 ), 0.)
  end if

end subroutine gem_reconnection_per


!==============================================================================
! Wave Propagation (itaylor = 4)
!==============================================================================
subroutine wave_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call static_equ(ibegin)

  phi0(ibegin) = z
  phi0(ibegin+1) = 0.0
  phi0(ibegin+2) = 1.0
  phi0(ibegin+3) = 0.0
  phi0(ibegin+4) = 0.0
  phi0(ibegin+5) = 0.0

  if(numvar.ge.2) call constant_field(phi0(ibegin+6 :ibegin+11), bzero)
  if(numvar.ge.3) call constant_field(phi0(ibegin+12:ibegin+17), p0-pi0*ipres)
  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den0(ibegin:ibegin+5),     1.)
  end if

end subroutine wave_equ


subroutine wave_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: alx, alz, akx, akx2, omega

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call getboundingboxsize(alx, alz)

  akx = 2.*pi/alx
  akx2 = akx**2
  omega = sqrt(akx2 + .5*db**2*akx2**2                                 &
       + db*akx2*sqrt(akx2 + .25*db**2*akx2**2))
  
  call plane_wave(x, z, vel(ibegin:ibegin+5), akx, 0, &
       eps*akx/omega, pi/2.)
  if(numvar.ge.2)  call plane_wave(x, z, vel(ibegin+6:ibegin+11),      &
       akx, 0, eps*(akx2/omega**2-1), pi/2.)
  if(numvar.ge.3)  call constant_field(vel(ibegin+12:ibegin+17), 0.)
  
  call plane_wave(x, z, phi(ibegin:ibegin+5), akx, 0, eps, pi/2.)
  if(numvar.ge.2) call plane_wave(x, z, phi(ibegin+6:ibegin+11), &
       akx, 0, eps*(akx/omega-omega/akx), pi/2.)
  if(numvar.ge.3)  call constant_field(phi(ibegin+12:ibegin+17), 0.)  

  if(idens.eq.1) then
     call entdofs(numvar, inode, 0, ibegin, iendplusone)
     call constant_field(den(ibegin:ibegin+5 ), 0.)
  end if

end subroutine wave_per


!==============================================================================
! Solov'ev Equilibrium (itor = 1, itaylor = 0)
!==============================================================================
subroutine solovev_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: psi_psib, psi_psib_x, psi_psib_xx
  real :: psi_psib_z, psi_psib_xz, psi_psib_zz
  real :: yb, e, d, y, r0, psib, f12, pp0, aspect_ratio

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call static_equ(ibegin)

  aspect_ratio = xmag / (xmag - xlim)
  yb = 2.*aspect_ratio / (1. + aspect_ratio**2)
  d = 0.1
  e = 3.
  psib = 1.
  r0 = xmag
  y = (x/r0)**2 - 1.

  psi_psib    = ((4./e**2)*(1.+(1.-d)*y)*(z/r0)**2+y**2)/yb**2
  psi_psib_x  = ((4./e**2)*(1.-d)*(2.*x/r0**2)*(z/r0)**2+4.*y*x/r0**2)/yb**2
  psi_psib_z  = ((8./e**2)*(1.+(1.-d)*y)*z/r0**2)/yb**2
  psi_psib_xx = ((4./e**2)*(1.-d)*(2./r0**2)*(z/r0)**2+4.*y/r0**2+8.*(x/r0**2)**2)/yb**2
  psi_psib_xz = ((8./e**2)*(1.-d)*(2.*x/r0**2)*z/r0**2)/yb**2
  psi_psib_zz = ((8./e**2)*(1.+(1.-d)*y)/r0**2)/yb**2

  phi0(ibegin)   = psib * psi_psib
  phi0(ibegin+1) = psib * psi_psib_x
  phi0(ibegin+2) = psib * psi_psib_z
  phi0(ibegin+3) = psib * psi_psib_xx
  phi0(ibegin+4) = psib * psi_psib_xz
  phi0(ibegin+5) = psib * psi_psib_zz

  if(numvar.ge.2) then
     f12 = d*(4.*psib/(yb*e*r0))**2

     phi0(ibegin+6) = sqrt(bzero**2 + f12*(1.- psi_psib))
     phi0(ibegin+7) = -0.5*f12*psi_psib_x/phi0(ibegin+6)
     phi0(ibegin+8) = -0.5*f12*psi_psib_z/phi0(ibegin+6)
     phi0(ibegin+9) = -0.5*f12*psi_psib_xx/phi0(ibegin+6) &
          + (0.5*f12*psi_psib_x)**2/phi0(ibegin+6)**3
     phi0(ibegin+10) = -0.5*f12*psi_psib_xz/phi0(ibegin+6) &
          + (0.5*f12)**2*psi_psib_x*psi_psib_z/phi0(ibegin+6)**3
     phi0(ibegin+11) = -0.5*f12*psi_psib_zz/phi0(ibegin+6) &
          + (0.5*f12*psi_psib_z)**2/phi0(ibegin+6)**3
  end if

  if(numvar.ge.3) then
     pp0 = 8*(psib/(yb*e*r0**2))**2 * (e**2 + 1 - d)

     phi0(ibegin+12) = pp0 * (1. - psi_psib)
     phi0(ibegin+13) = - pp0 * psi_psib_x
     phi0(ibegin+14) = - pp0 * psi_psib_z
     phi0(ibegin+15) = - pp0 * psi_psib_xx
     phi0(ibegin+16) = - pp0 * psi_psib_xz
     phi0(ibegin+17) = - pp0 * psi_psib_zz
  end if


  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den0(ibegin:ibegin+5),1.)
  endif

end subroutine solovev_equ


subroutine solovev_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call constant_field(vel(ibegin:ibegin+5), 0.)
  if(numvar.ge.2)  call constant_field(vel(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(vel(ibegin+12:ibegin+17), 0.)
  
  call constant_field(phi(ibegin:ibegin+5), 0.)
  if(numvar.ge.2)  call constant_field(phi(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(phi(ibegin+12:ibegin+17), 0.)

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den(ibegin+12:ibegin+17), 0.)
  endif

end subroutine solovev_per







!===========================================================================
! OLD STUFF
! ==========================================================================

!=====================================
subroutine velequ(dum, numberingid)
  
  use basic
  
  implicit none
  integer :: numnodes, l, ibegin, iendplusone, i, numberingid
  real :: dum(*)

  real :: alx, alz, akx, akz, x, z, xmin, zmin
  double precision :: coords(3)

  call getboundingboxsize(alx, alz)

  ! initialize equilibrium velocity variables to zero 
  call numnod(numnodes)
  do l=1, numnodes
     call entdofs(numberingid, l, 0, ibegin, iendplusone)

     call xyznod(l, coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     dum(ibegin:ibegin+5) = 0
     if(numvar.ge.2) dum(ibegin+6 :ibegin+11) = 0
     if(numvar.ge.3) dum(ibegin+12:ibegin+17) = 0
  enddo

  return
end subroutine velequ

!============================================================
subroutine phiequ(dum, numberingid)
  
  use basic
  implicit none
  integer :: l, i, ifail1, ifail2, ifail3, numnodes, ibegin, iendplusone, numberingid
  real :: dum(*), x, z, r, ri, arg, befo, s17aef, s17aff
  real :: ff, fp, fpp, alam, sech, alx, alz, xmin, zmin
  real :: k, j0, j1, kb, akx, akz
  double precision :: coords(3)

  k = 3.8317059702

  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0

  call numnod(numnodes)
  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)
  
  do l=1,numnodes
     call entdofs(numberingid, l, 0, ibegin, iendplusone)
     call xyznod(l, coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5
     
     select case(itaylor)
     
     case(0)
!
!.....for itaylor = 0, use the equilibrium given in R. Richard, et al,
!     Phys. Fluids B, 2 (3) 1990 p 489
!.....(also Strauss and Longcope)
!
        r = sqrt(x**2 + z**2)
        ri = 0.
        if(r.ne.0) ri = 1./r
        arg = k*r
        if(r.le.1) then
           befo = 2./s17aef(k,ifail1)
           j0 = s17aef(arg,ifail2)
           j1 = s17aff(arg,ifail3)
           ff = .5*befo
           fp = 0
           fpp = -.125*k**2*befo
           if(arg.ne.0)  then
              ff   = befo*j1/arg
              fp  = befo*(j0 - 2.*j1/arg)/r
              fpp = befo*(-3*j0 + 6.*j1/arg - arg*j1)/r**2
           endif
        else
           ff   = (1.-1./r**2)
           fp  = 2./r**3
           fpp = -6/r**4
        endif
   
        dum(ibegin) = ff* z
        dum(ibegin+1) = fp*z*x *ri
        dum(ibegin+2) = fp*z**2*ri + ff
        dum(ibegin+3) = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)
        dum(ibegin+4) = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
        dum(ibegin+5) = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)


        if(numvar.ge.2) then

           if(r.gt.1) then

              if(numvar.ge.2) then
                 dum(ibegin+6) = bzero
                 do i=ibegin+7,ibegin+11
                    dum(i) = 0.
                 enddo
              endif
              if(numvar.ge.3) then
                 if(ipres.eq.0) then
                    dum(ibegin+12) = p0
                 else
                    dum(ibegin+12) = p0 - pi0
                 endif
                 do i=ibegin+13,ibegin+17
                    dum(i) = 0.
                 enddo
              endif

           else 

              kb = k**2*(1.-beta)
              dum(ibegin+6) = sqrt(kb*dum(ibegin)**2+bzero**2)
              dum(ibegin+7) = kb/dum(ibegin+6)*dum(ibegin)*dum(ibegin+1)
              dum(ibegin+8) = kb/dum(ibegin+6)*dum(ibegin)*dum(ibegin+2)
              dum(ibegin+9) = kb/dum(ibegin+6)                               &
                   *(-kb/dum(ibegin+6)**2*(dum(ibegin)*dum(ibegin+1))**2     &
                   + dum(ibegin+1)**2+dum(ibegin)*dum(ibegin+3))
              dum(ibegin+10) = kb/dum(ibegin+6)*(-kb/dum(ibegin+6)**2*       &
                   dum(ibegin)**2*dum(ibegin+1)*dum(ibegin+2)                &
                   + dum(ibegin+1)*dum(ibegin+2)+dum(ibegin)*dum(ibegin+4))
              dum(ibegin+11) = kb/dum(ibegin+6)                              &
                   *(-kb/dum(ibegin+6)**2*(dum(ibegin)*dum(ibegin+2))**2     &
                   + dum(ibegin+2)**2+dum(ibegin)*dum(ibegin+5))
              
              if(numvar.ge.3) then
                 kb = k**2*beta
                 if(ipres.eq.0) then
                    dum(ibegin+12) = 0.5*kb*dum(ibegin)**2 + p0
                 else
                    dum(ibegin+12) = 0.5*kb*dum(ibegin)**2 + p0 - pi0
                 endif
                 dum(ibegin+13) = kb*dum(ibegin)*dum(ibegin+1)
                 dum(ibegin+14) = kb*dum(ibegin)*dum(ibegin+2)
                 dum(ibegin+15) = kb*(dum(ibegin+1)**2+                      &
                      dum(ibegin)*dum(ibegin+3))
                 dum(ibegin+16) = kb*(dum(ibegin+1)*dum(ibegin+2)+           &
                      dum(ibegin)*dum(ibegin+4))
                 dum(ibegin+17) = kb*(dum(ibegin+2)**2+                      &
                      dum(ibegin)*dum(ibegin+5))
              endif
           endif
        endif

      case(1) !  Taylor reconnection (itaylor=1)

         dum(ibegin) = -0.5*z**2
         dum(ibegin+1) = 0.
         dum(ibegin+2) = -z
         dum(ibegin+4) = 0
         dum(ibegin+4) = 0
         dum(ibegin+5) = -1

      case(2)  !  force-free taylor state (itaylor=2)
         call xyznod(l, coords)
         x = coords(1)
         z = coords(2)
         alam = sqrt((pi/alx)**2+(pi/alz)**2)

         dum(ibegin) = sin(pi*x/alx)*sin(pi*z/alz)
         dum(ibegin+1) = (pi/alx)*cos(pi*x/alx)*sin(pi*z/alz)
         dum(ibegin+2) = (pi/alz)*sin(pi*x/alx)*cos(pi*z/alz)
         dum(ibegin+3) = -(pi/alx)**2*sin(pi*x/alx)*sin(pi*z/alz)
         dum(ibegin+4) = (pi/alx)*(pi/alz)*cos(pi*x/alx)*cos(pi*z/alz)
         dum(ibegin+5) = -(pi/alz)**2*sin(pi*x/alx)*sin(pi*z/alz)

         dum(ibegin+6) = alam*dum(ibegin)
         dum(ibegin+7) = alam*dum(ibegin+1)
         dum(ibegin+8) = alam*dum(ibegin+2)
         dum(ibegin+9) = alam*dum(ibegin+3)
         dum(ibegin+10) = alam*dum(ibegin+4)
         dum(ibegin+11) = alam*dum(ibegin+5)

      case(3)  ! GEM reconnection (itaylor=3)
         
         dum(ibegin) = 0.5*alog(cosh(2.*z))
         dum(ibegin+1) = 0.0
         dum(ibegin+2) = tanh(2.*z)  
         dum(ibegin+3) = 0.0
         dum(ibegin+4) = 0.0
         dum(ibegin+5) = 2.*(1.-tanh(2.*z)**2)

         if(numvar.ge.2) then
            dum(ibegin+6) = bzero
            do i=ibegin+7,ibegin+11
               dum(i) = 0.
            enddo
         endif

         if(numvar.ge.3) then
            dum(ibegin+12) = 0.5*sech(2*z)**2 + 0.1
            dum(ibegin+13) = 0.
            dum(ibegin+14) = -2.*sech(2*z)**2*tanh(2*z)
            dum(ibegin+15) = 0.
            dum(ibegin+16) = 0.
            dum(ibegin+17) = -4.*sech(2*z)**2*             &
                 &        (sech(2*z)**2-2.*tanh(2*z)**2)
         endif ! if on numvar.ge.3
         
      case(4)  ! wave propagation (itaylor=4)
         
         dum(ibegin) = z
         dum(ibegin+1) = 0.0
         dum(ibegin+2) = 1.0
         dum(ibegin+3) = 0.0
         dum(ibegin+4) = 0.0
         dum(ibegin+5) = 0.0
         
         if(numvar.ge.2) then
            dum(ibegin+6) = bzero
            do i=ibegin+7,ibegin+11
               dum(i) = 0.
            enddo
         endif

      case(7)  ! heat conduction test

         dum(ibegin) = z
         dum(ibegin+1) = 0.0
         dum(ibegin+2) = 1.0
         dum(ibegin+3) = 0.0
         dum(ibegin+4) = 0.0
         dum(ibegin+5) = 0.0
         
         if(numvar.ge.2) then
            dum(ibegin+6) = bzero
            do i=ibegin+7,ibegin+11
               dum(i) = 0.
            enddo
         endif

      case(100) ! test case
         
         akx = 2.*pi/alx
         akz = 2.*2.*pi/alz

         dum(ibegin  ) = cos(akx*x)*cos(akz*z)
         dum(ibegin+1) =-sin(akx*x)*cos(akz*z)*akx
         dum(ibegin+2) =-cos(akx*x)*sin(akz*z)*akz
         dum(ibegin+3) =-cos(akx*x)*cos(akz*z)*akx**2
         dum(ibegin+4) = sin(akx*x)*sin(akz*z)*akx*akz
         dum(ibegin+5) =-cos(akx*x)*cos(akz*z)*akz**2
      end select
   end do

end subroutine phiequ
!============================================================
subroutine denequ(dum, numberingid)
  use basic
  implicit none

  real, intent(inout) :: dum(*)
  integer, intent(in) :: numberingid

  integer :: l, numnodes, ii, ifail1, ifail2
  integer :: ifail3, ibegin, iendplusone

  real :: x, z
  real :: j0, j1, k, kb, r, ri, arg, befo, ff, s17aef, s17aff
  real :: fp, fpp, d1, d3, d2, d6, d5, d4, sech
  double precision :: coords(3)
  real :: alx, alz, xmin, zmin

  k = 3.8317059702

  call getboundingboxsize(alx, alz)
  call getmincoord(xmin, zmin)

  ! initialize equilibrium to be the Strauss and Longcope solution
  
  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0
  call numnod(numnodes)

  print *, 'reached dd', numnodes

  do l=1,numnodes
     call entdofs(numberingid, l, 0, ibegin, iendplusone)
     call xyznod(l,coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     if(idens.eq.0) then
        dum(ibegin) = 1.
        dum(ibegin+1) = 0.
        dum(ibegin+2) = 0.
        dum(ibegin+3) = 0.
        dum(ibegin+4) = 0.
        dum(ibegin+5) = 0.
        cycle
     endif

     select case(itaylor)

     case(0)         

        ! for itaylor = 0, use the equilibrium given in R. Richard, et al,
        !     Phys. Fluids B, 2 (3) 1990 p 489

        r = sqrt(x**2 + z**2)
        ri = 0.
        if(r.ne.0) ri = 1./r
        arg = k*r
        if(r.le.1) then
           befo = 2./s17aef(k,ifail1)
           j0 = s17aef(arg,ifail2)
           j1 = s17aff(arg,ifail3)
           ff = .5*befo
           fp = 0
           fpp = -.125*k**2*befo
           if(arg.ne.0)  then
              ff   = befo*j1/arg
              fp  = befo*(j0 - 2.*j1/arg)/r
              fpp = befo*(-3*j0 + 6.*j1/arg - arg*j1)/r**2
           endif
        else
           ff   = (1.-1./r**2)
           fp  = 2./r**3
           fpp = -6/r**4
        endif

        d1 = ff* z
        d3 = fp*z**2*ri + ff
        d2 = fp*z*x *ri
        d6 = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)
        d5 = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
        d4 = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)

        kb = k**2*beta
        dum(ibegin) = 0.5*kb*d1**2 + 1.

        dum(ibegin+1) = kb*d1*d2
        dum(ibegin+2) = kb*d1*d3
        dum(ibegin+3) = kb*(d2**2+d1*d4)
        dum(ibegin+4) = kb*(d2*d3+d1*d5)
        dum(ibegin+5) = kb*(d3**2+d1*d6)

     case(3) ! GEM reconnection (itaylor=3)
        dum(ibegin) = sech(2*z)**2 + 0.2
        dum(ibegin+1) = 0.
        dum(ibegin+2) = -4.*sech(2*z)**2*tanh(2*z)
        dum(ibegin+3) = 0.
        dum(ibegin+4) = 0.
        dum(ibegin+5) = -8.*sech(2*z)**2*(sech(2*z)**2-2.*tanh(2*z)**2)

     case default
        dum(ibegin) = 1.
        do ii=1,5
           dum(ibegin+ii) = 0.
        enddo
        
     end select
     
     
  end do

end subroutine denequ
!========================
subroutine velinit(dum)

  use basic
  implicit none
  real dum(*), x, z, akx, akz, akx2, omega, befo, alx, alz, xmin, zmin
  integer numnodes, l, ibegin, iendplusone
  double precision coords(3) 

  call numnod(numnodes)
  call getboundingboxsize(alx, alz)
  call getmincoord(xmin, zmin)

  ! initialize perturbed velocity variables

  do l=1,numnodes
     call entdofs(numvar, l, 0, ibegin, iendplusone)
     call xyznod(l,coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5
   
     dum(ibegin:iendplusone-1) = 0.

     select case(itaylor)

     case(4) ! wave propagation (itaylor=4)
        akx = 2.*pi/alx
        akx2 = akx**2
        omega = sqrt(akx2 + .5*db**2*akx2**2                                 &
             + db*akx2*sqrt(akx2 + .25*db**2*akx2**2))
        befo = eps*akx/omega
        
        dum(ibegin) = befo*sin(akx*x)
        dum(ibegin+1) = befo*akx*cos(akx*x)
        dum(ibegin+2) =  0.
        dum(ibegin+3) =  -befo*akx**2*sin(akx*x)
        dum(ibegin+4) = 0.
        dum(ibegin+5) =  0.
        
        if(numvar.ge.2) then
           befo = eps*(akx2/omega**2-1)
           
           dum(ibegin+6) = befo*sin(akx*x)
           dum(ibegin+7) = befo*akx*cos(akx*x)
           dum(ibegin+8) =  0.
           dum(ibegin+9) =  -befo*akx**2*sin(akx*x)
           dum(ibegin+10) = 0.
           dum(ibegin+11) =  0.
        endif

     case(100) ! test case
         
        akx = 3.*2.*pi/alx
        akz = 4.*2.*pi/alz
         
        dum(ibegin  ) = cos(akx*x)*cos(akz*z)
        dum(ibegin+1) =-sin(akx*x)*cos(akz*z)*akx
        dum(ibegin+2) =-cos(akx*x)*sin(akz*z)*akz
        dum(ibegin+3) =-cos(akx*x)*cos(akz*z)*akx**2
        dum(ibegin+4) = sin(akx*x)*sin(akz*z)*akx*akz
        dum(ibegin+5) =-cos(akx*x)*cos(akz*z)*akz**2
        
     end select

  end do

end subroutine velinit
!===========================================================================
subroutine phiinit(dum)
  use basic
  implicit none
  
  integer :: l, numnodes, ibegin, iendplusone
  real x, z, akx, akz, akx2, omega, befo, alx, alz, xmin, zmin
  real dum(*)
  real :: peper, r, r0, loc, locp
  double precision :: coords(3)

! start a loop over the number of vertices to define initial conditions

! note coordinates are centered at 0,0
  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)
  call numnod(numnodes)

 
  do l=1,numnodes
     call xyznod(l,coords)
     call entdofs(numvar, l, 0, ibegin, iendplusone)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     dum(ibegin:iendplusone-1) = 0.

     select case(itaylor)

     case(0,1,2)
     
     case(3) ! GEM reconnection (itaylor=3)
        akx = 2.*pi/alx
        akz = pi/alz

        dum(ibegin) = eps*cos(akx*x)*cos(akz*z)
        dum(ibegin+1) = -eps*akx*sin(akx*x)*cos(akz*z)
        dum(ibegin+2) =  - eps*akz*cos(akx*x)*sin(akz*z)
        dum(ibegin+3) = -eps*akx**2*cos(akx*x)*cos(akz*z)
        dum(ibegin+4) = eps*akx*akz*sin(akx*x)*sin(akz*z)
        dum(ibegin+5) =  - eps*akz**2*cos(akx*x)*cos(akz*z)

     case(4) ! wave propagation (itaylor=4)
        akx = 2.*pi/alx
        akx2 = akx**2
        omega = sqrt(akx2 + .5*db**2*akx2**2                                 &
             + db*akx2*sqrt(akx2 + .25*db**2*akx2**2))

        dum(ibegin) =  eps*sin(akx*x)
        dum(ibegin+1) =  eps*akx*cos(akx*x)
        dum(ibegin+2) =  0.
        dum(ibegin+3) = -eps*akx**2*sin(akx*x)
        dum(ibegin+4) = 0.
        dum(ibegin+5) =  0.
        
        if(numvar.ge.2) then
           befo = eps*(akx/omega-omega/akx)
           dum(ibegin+6) =  befo*sin(akx*x)
           dum(ibegin+7) =  befo*akx*cos(akx*x)
           dum(ibegin+8) =  0.
           dum(ibegin+9) = -befo*akx**2*sin(akx*x)
           dum(ibegin+10) = 0.
           dum(ibegin+11) =  0.
        endif ! on numvar.ge.2

     case(7) ! heat conduction test
        if(ipres.eq.1) then
           peper = eps*(p0-pi0)        
        else
           peper = eps*p0
        endif

        r = sqrt(x**2+z**2+eps)
        r0 = 0.0
        
        loc = alx/ln
        locp = alz/ln

        dum(ibegin+12) = peper*exp(-((x-r0)/loc)**2-(z/locp)**2)
        dum(ibegin+13) = -2.*(x-r0)*dum(ibegin+12)/loc**2
        dum(ibegin+14) = -2.*z*dum(ibegin+12)/locp**2
        dum(ibegin+15) = 2.*(2.*((x-r0)/loc)**2 - 1.)*dum(ibegin+12)/loc**2
        dum(ibegin+16) = 4.*(x-r0)*z*dum(ibegin+12)/(loc*locp)**2
        dum(ibegin+17) = 2.*(2.*(z/locp)**2 - 1.)*dum(ibegin+12)/locp**2
   
     end select

  end do

end subroutine phiinit
!============================================================
subroutine deninit(dum)
  use basic

  implicit none

  real, intent(out) :: dum(*)
  
  integer :: l, numnodes, ibegin, iendplusone
  real :: x, z, alx, alz, xmin, zmin

  double precision :: coords(3)


  ! note coordinates are centered at 0,0
  call numnod(numnodes)
  call getboundingboxsize(alx, alz)
  call getmincoord(xmin, zmin)

  do l=1,numnodes
     call xyznod(l,coords)
     call entdofs(1, l, 0, ibegin, iendplusone)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     dum(ibegin) = 0.
     dum(ibegin+1) = 0.
     dum(ibegin+2) = 0.
     dum(ibegin+3) = 0.
     dum(ibegin+4) = 0.
     dum(ibegin+5) = 0.
  end do

end subroutine deninit

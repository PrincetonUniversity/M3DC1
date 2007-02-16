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
module tilting_cylinder

  implicit none

  real, parameter :: k = 3.8317059702

contains

subroutine tilting_cylinder_init()
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

     call cylinder_equ(x, z, l)
     call cylinder_per(x, z, l)
  enddo

end subroutine tilting_cylinder_init

subroutine cylinder_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: rr, ri, arg, befo, uzero, czero, ff, fp, fpp, j0, j1, kb
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
     call constant_field(den0(ibegin:ibegin+5), 1.)
  endif

end subroutine cylinder_equ

subroutine cylinder_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: befo, uzero, czero
  
  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  befo = eps*exp(-(x**2+z**2))
  uzero = befo*cos(z)
  czero = befo*sin(z)
  vel(ibegin  ) = uzero
  vel(ibegin+1) = - 2.*x*uzero
  vel(ibegin+2) = - 2.*z*uzero - czero
  vel(ibegin+3) = (-2.+4.*x**2)*uzero
  vel(ibegin+4) = 4.*x*z*uzero + 2.*x*czero
  vel(ibegin+5) = (-3.+4.*z**2)*uzero + 4.*z*czero
  
  if(numvar.ge.2) call constant_field(vel(ibegin+6:ibegin+11), 0.)

  if(numvar.ge.3) then
     vel(ibegin+12) = czero
     vel(ibegin+13) = -2*x*czero
     vel(ibegin+14) = -2*z*czero + uzero
     vel(ibegin+15) = (-2.+4*x**2)*czero
     vel(ibegin+16) = 4.*x*z*czero - 2.*x*uzero
     vel(ibegin+17) = (-3.+4.*z**2)*czero - 4*z*uzero
  endif

  call constant_field(phi(ibegin   :ibegin+5 ), 0.)
  if(numvar.ge.2)  call constant_field(phi(ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(phi(ibegin+12:ibegin+17), 0.)

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(den(ibegin:ibegin+5 ), 0.)
  endif

end subroutine cylinder_per

end module tilting_cylinder



!==============================================================================
! Taylor Reconnection (itaylor = 1)
!==============================================================================

module taylor_reconnection

contains


subroutine taylor_reconnection_init()
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

     call taylor_reconnection_equ(x, z, l)
     call taylor_reconnection_per(x, z, l)
  enddo

end subroutine taylor_reconnection_init



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

end subroutine taylor_reconnection_per


end module taylor_reconnection


!==============================================================================
! Force-Free Taylor State (itaylor = 2)
!==============================================================================
module force_free_state

contains

subroutine force_free_init()
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

     call force_free_equ(x, z, l)
     call force_free_per(x, z, l)
  enddo

end subroutine force_free_init

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

end module force_free_state

!==============================================================================
! GEM Reconnection (itaylor = 3)
!==============================================================================
module gem_reconnection

contains

subroutine gem_reconnection_init()
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

     call gem_reconnection_equ(x, z, l)
     call gem_reconnection_per(x, z, l)
  enddo

end subroutine gem_reconnection_init

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

     pezero = p0 - pi0*ipres
    
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

  call constant_field(vel(ibegin:ibegin+5), 0.)
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
     call constant_field(den(ibegin:ibegin+5), 0.)
  end if

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

     call wave_equ(x, z, l)
     call wave_per(x, z, l)
  enddo

end subroutine wave_init

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
  real :: alx, alz, omega

  call entdofs(numvar, inode, 0, ibegin, iendplusone)
  
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

end module wave_propagation


!==============================================================================
! Solov'ev Equilibrium (itor = 1, itaylor = 0)
!==============================================================================
module solovev

contains

subroutine solovev_init()
  use basic

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

     call solovev_equ(x, z, l)
     call solovev_per(x, z, l)
  enddo

end subroutine solovev_init

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
  e = 1.5
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
          - (0.5*f12*psi_psib_x)**2/phi0(ibegin+6)**3
     phi0(ibegin+10) = -0.5*f12*psi_psib_xz/phi0(ibegin+6) &
          - (0.5*f12)**2*psi_psib_x*psi_psib_z/phi0(ibegin+6)**3
     phi0(ibegin+11) = -0.5*f12*psi_psib_zz/phi0(ibegin+6) &
          - (0.5*f12*psi_psib_z)**2/phi0(ibegin+6)**3
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
     call constant_field(den(ibegin:ibegin+5), 0.)
  endif

end subroutine solovev_per

end module solovev




!=====================================
subroutine initial_conditions()
  use basic

  use tilting_cylinder
  use taylor_reconnection
  use force_free_state
  use gem_reconnection
  use wave_propagation
  use solovev
  use gradshafranov

  implicit none

  integer :: ibegin, iendplusone, l, numnodes
  real :: x, z, alx, alz, xmin, zmin
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

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
     end select
  else
     ! toroidal equilibria
     select case(itaylor)
     case(0)
        call solovev_init()
     case(1)
        call gradshafranov_init()
     end select
  endif

end subroutine initial_conditions

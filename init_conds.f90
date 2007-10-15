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

  call constant_field(phi0_v(ibegin+phi_off:ibegin+phi_off+5 ), 0.)
  if(numvar.ge.2) &
       call constant_field(vz0_v(ibegin+vz_off:ibegin+vz_off+5), vzero)
  if(numvar.ge.3) &
       call constant_field(chi0_v(ibegin+chi_off:ibegin+chi_off+5), 0.)
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

subroutine cartesian_to_cylindrical(x,vec)
  implicit none

  real, dimension(6), intent(inout) :: vec
  real, intent(in) :: x

  vec(6) = vec(6) * x
  vec(5) = vec(5) * x +    vec(3)
  vec(4) = vec(4) * x + 2.*vec(2)
  vec(3) = vec(3) * x
  vec(2) = vec(2) * x +    vec(1)
  vec(1) = vec(1) * x 
end subroutine cartesian_to_cylindrical

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

     call assign_vectors(inode)

     call cartesian_to_cylindrical(x,psi0_l(1:6))
     call cartesian_to_cylindrical(x,psi1_l(1:6))
     call cartesian_to_cylindrical(x,phi0_l(1:6))
     call cartesian_to_cylindrical(x,phi1_l(1:6))
     
     if(numvar.ge.2) then
        call cartesian_to_cylindrical(x,bz0_l(1:6))
        call cartesian_to_cylindrical(x,bz1_l(1:6))
        call cartesian_to_cylindrical(x,vz0_l(1:6))
        call cartesian_to_cylindrical(x,vz1_l(1:6))
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

     call assign_vectors(l)

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

  call constant_field(phi1_l(1:6), 0.)
  if(numvar.ge.2)  call constant_field( vz1_l(1:6), 0.)
  if(numvar.ge.3)  call constant_field(chi1_l(1:6), 0.)

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
     if(numvar.ge.2) call constant_field(bz0_l(1:6), bzero   )
     if(numvar.ge.3) call constant_field(pe0_l(1:6), p0-pi0*ipres)
     if(ipres.eq.1)  call constant_field( p0_l(1:6), p0)
  else 
     if(numvar.ge.2) then
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
        
        if(numvar.ge.3) then
           kb = k**2*beta*(p0 - pi0*ipres)/p0
           pe0_l(1) = 0.5*kb*psi0_l(1)**2 + p0 - pi0*ipres
           pe0_l(2) = kb*psi0_l(1)*psi0_l(2)
           pe0_l(3) = kb*psi0_l(1)*psi0_l(3)
           pe0_l(4) = kb*(psi0_l(2)**2+psi0_l(1)*psi0_l(4))
           pe0_l(5) = kb*(psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
           pe0_l(6) = kb*(psi0_l(3)**2+psi0_l(1)*psi0_l(6))
        endif
     endif

     if(ipres.eq.1) then
        kb = k**2*beta
        p0_l(1) = 0.5*kb*psi0_l(1)**2 + p0
        p0_l(2) = kb*psi0_l(1)*psi0_l(2)
        p0_l(3) = kb*psi0_l(1)*psi0_l(3)
        p0_l(4) = kb*(psi0_l(2)**2+psi0_l(1)*psi0_l(4))
        p0_l(5) = kb*(psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
        p0_l(6) = kb*(psi0_l(3)**2+psi0_l(1)*psi0_l(6))
     endif
  endif

  if(idens.eq.1) then
     call constant_field(den0_l(1:6), 1.)
  endif

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
  phi1_l(1) = uzero
  phi1_l(2) = - 2.*x*uzero
  phi1_l(3) = - 2.*z*uzero - czero
  phi1_l(4) = (-2.+4.*x**2)*uzero
  phi1_l(5) = 4.*x*z*uzero + 2.*x*czero
  phi1_l(6) = (-3.+4.*z**2)*uzero + 4.*z*czero
  
  if(numvar.ge.2) call constant_field(vz1_l(1:6), 0.)

  if(numvar.ge.3) then
     chi1_l(1) = czero
     chi1_l(2) = -2*x*czero
     chi1_l(3) = -2*z*czero + uzero
     chi1_l(4) = (-2.+4*x**2)*czero
     chi1_l(5) = 4.*x*z*czero - 2.*x*uzero
     chi1_l(6) = (-3.+4.*z**2)*czero - 4*z*uzero
  endif

  call constant_field(psi1_l(1:6), 0.)
  if(numvar.ge.2)  call constant_field( bz1_l(1:6), 0.)
  if(numvar.ge.3)  call constant_field( pe1_l(1:6), 0.)
  if(ipres.eq.1)   call constant_field(  p1_l(1:6), 0.)
  if(idens.eq.1)   call constant_field(den1_l(1:6), 0.)

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

  integer :: ibegin, iendplusone, ibegin1, iendplusone1
  real :: sech, pezero

  call entdofs(numvar, inode, 0, ibegin, iendplusone)
  if(idens.eq.1 .or. ipres.eq.1) then
     call entdofs(1, inode, 0, ibegin1, iendplusone1)
  endif

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

  if(ipres.eq.1) then
     pres0(ibegin1  ) = p0*(sech(2.*z)**2 + 0.2)
     pres0(ibegin1+1) = 0.
     pres0(ibegin1+2) = p0*(-4.*sech(2.*z)**2*tanh(2.*z))
     pres0(ibegin1+3) = 0.
     pres0(ibegin1+4) = 0.
     pres0(ibegin1+5) = p0*(-8.*sech(2.*z)**2*(sech(2.*z)**2-2.*tanh(2.*z)**2)) 
  endif

  if(idens.eq.1) then
     den0(ibegin1  ) = sech(2*z)**2 + 0.2
     den0(ibegin1+1) = 0.
     den0(ibegin1+2) = -4.*sech(2*z)**2*tanh(2*z)
     den0(ibegin1+3) = 0.
     den0(ibegin1+4) = 0.
     den0(ibegin1+5) = -8.*sech(2*z)**2*(sech(2*z)**2-2.*tanh(2*z)**2)
  endif

end subroutine gem_reconnection_equ

subroutine gem_reconnection_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone, ibegin1, iendplusone1
  real :: akx, akz, alx, alz

  call entdofs(numvar, inode, 0, ibegin, iendplusone)
  if(idens.eq.1 .or. ipres.eq.1) then
     call entdofs(1, inode, 0, ibegin1, iendplusone1)
  endif

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

  if(numvar.ge.2)  call constant_field(phi (ibegin+6 :ibegin+11), 0.)
  if(numvar.ge.3)  call constant_field(phi (ibegin+12:ibegin+17), 0.)  
  if(ipres.eq.1)   call constant_field(pres(ibegin1  :ibegin1+5), 0.)
  if(idens.eq.1)   call constant_field(den (ibegin1  :ibegin1+5), 0.)

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
  real :: omega

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
! Gravitational Instability Equilibrium (itor = 0, itaylor = 5)
!==============================================================================
module grav

contains

subroutine grav_init()
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
     z = coords(2) + zzero - zmin

     call grav_equ(x, z, l)
     call grav_per(x, z, l)
  enddo

end subroutine grav_init

subroutine grav_equ(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: fac1, n0, pn0, ppn0

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  n0 = exp(z/ln)
  pn0 = n0/ln
  ppn0 = pn0/ln

  call static_equ(ibegin)

  call constant_field(phi0(ibegin:ibegin+5), 0.)

  if(numvar.ge.2) then
     fac1 = gravz*ln+gam*p0
     phi0(ibegin+6 ) = sqrt(bzero**2 - 2.*fac1*(n0-1.))
     phi0(ibegin+7 ) = 0.
     phi0(ibegin+8 ) = -pn0*fac1/phi0(ibegin+6)
     phi0(ibegin+9 ) = 0.
     phi0(ibegin+10) = 0.
     phi0(ibegin+11) = (fac1/phi0(ibegin+6))*(pn0*phi0(ibegin+8)/phi0(ibegin+6) - ppn0)
  end if

  if(numvar.ge.3) then
     fac1 = p0 - pi0*ipres
     phi0(ibegin+12) = fac1+gam*fac1*(n0-1.)
     phi0(ibegin+13) = 0.
     phi0(ibegin+14) = gam*fac1*pn0
     phi0(ibegin+15) = 0.
     phi0(ibegin+16) = 0.
     phi0(ibegin+17) = gam*fac1*ppn0
  end if

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     den0(ibegin  ) = n0
     den0(ibegin+1) = 0.
     den0(ibegin+2) = pn0
     den0(ibegin+3) = 0.
     den0(ibegin+4) = 0.
     den0(ibegin+5) = ppn0
  endif

  if(ipres.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     fac1 = p0
     pres0(ibegin  ) = fac1+gam*fac1*(n0-1.)
     pres0(ibegin+1) = 0.
     pres0(ibegin+2) = gam*fac1*pn0
     pres0(ibegin+3) = 0.
     pres0(ibegin+4) = 0.
     pres0(ibegin+5) = gam*fac1*ppn0
  endif

end subroutine grav_equ


subroutine grav_per(x, z, inode)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  integer, intent(in) :: inode

  integer :: ibegin, iendplusone
  real :: kx, kz, alx, alz

  call getboundingboxsize(alx, alz)
  kx = 2.*pi/alx
  kz = pi/alz

  call entdofs(numvar, inode, 0, ibegin, iendplusone)

  call constant_field(phi(ibegin:ibegin+5), 0.)
  call constant_field(vel(ibegin:ibegin+5), 0.)
  
  if(numvar.ge.2)  then
     call constant_field(phi(ibegin+6 :ibegin+11), 0.)
     call constant_field(vel(ibegin+6 :ibegin+11), 0.)
  endif
  if(numvar.ge.3)  then
     call constant_field(phi(ibegin+12:ibegin+17), 0.)
     call constant_field(vel(ibegin+12:ibegin+17), 0.)
  endif

  if(idens.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     den(ibegin  ) =  eps*sin(kx*(x-xzero))*sin(kz*(z-zzero))
     den(ibegin+1) =  eps*cos(kx*(x-xzero))*sin(kz*(z-zzero))*kx
     den(ibegin+2) =  eps*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kz
     den(ibegin+3) = -eps*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kx**2
     den(ibegin+4) =  eps*cos(kx*(x-xzero))*cos(kz*(z-zzero))*kx*kz
     den(ibegin+5) = -eps*sin(kx*(x-xzero))*sin(kz*(z-zzero))*kz**2
  endif

  if(ipres.eq.1) then
     call entdofs(1, inode, 0, ibegin, iendplusone)
     call constant_field(pres(ibegin:ibegin+5),0.)
  endif

end subroutine grav_per

end module grav




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

  aspect_ratio = xmag / abs(xmag - xlim)
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

     call assign_vectors(l)

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

  call constant_field(phi0_l(1:6),0.)

  psi0_l(1) = bzero * x**2 / 2.
  psi0_l(2) = bzero * x
  psi0_l(3) = 0.
  psi0_l(4) = bzero
  psi0_l(5) = 0.
  psi0_l(6) = 0.

  if(numvar.ge.2) then
     call constant_field(bz0_l(1:6),0.)
     vz0_l(1)  = fac1*sqrt(x) + (3./8.)*nu - fac2/sqrt(x)
     vz0_l(2)  = 0.5*fac1/sqrt(x) + 0.5*fac2/(sqrt(x)**3)
     vz0_l(3)  = 0.
     vz0_l(4)  = -0.25*fac1/(sqrt(x)**3) - 0.75*fac2/(sqrt(x)**5)
     vz0_l(5) = 0.
     vz0_l(6) = 0.
  end if

  if(numvar.ge.3) then
     call constant_field( pe0_l(1:6),p0 - pi0*ipres)
     call constant_field(chi0_l(1:6),0.)
  end if

  if(idens.eq.1) call constant_field(den0_l(1:6),1.)

  if(ipres.eq.1) call constant_field(p0_l(1:6),p0)

end subroutine mri_equ


subroutine mri_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: fac1

  call constant_field(phi1_l(1:6), 0.)
  if(numvar.ge.2)  then
     fac1 = eps*sqrt(gravr*x)
     vz1_l(1) =  fac1*sin(kx*(x-xzero))*sin(kz*(z-zzero))
     vz1_l(2) =  fac1*cos(kx*(x-xzero))*sin(kz*(z-zzero))*kx
     vz1_l(3) =  fac1*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kz
     vz1_l(4) = -fac1*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kx**2
     vz1_l(5) =  fac1*cos(kx*(x-xzero))*cos(kz*(z-zzero))*kx*kz
     vz1_l(6) = -fac1*sin(kx*(x-xzero))*sin(kz*(z-zzero))*kz**2
  endif
  if(numvar.ge.3)  call constant_field(chi1_l(1:6), 0.)
  
  call constant_field(psi1_l(1:6), 0.)
  if(numvar.ge.2) call constant_field(bz1_l(1:6), 0.)
  if(numvar.ge.3) call constant_field(pe1_l(1:6), 0.)

  if(idens.eq.1) call constant_field(den1_l(1:6), 0.)
  if(ipres.eq.1) call constant_field(p1_l(1:6),0.)

end subroutine mri_per

end module mri




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
  use mri
  use grav

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
     end select
  else
     ! toroidal equilibria
     select case(itaylor)
     case(0)
!        call solovev_init()
        call tilting_cylinder_init()
        call cartesian_to_cylindrical_all()
     case(1)
        call gradshafranov_init()
     case(2)
        call mri_init()
     end select
  endif

end subroutine initial_conditions

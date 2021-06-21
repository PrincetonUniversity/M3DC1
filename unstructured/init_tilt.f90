!==============================================================================
! Tilting Cylinder (itaylor = 0)
!==============================================================================
module tilting_cylinder

  implicit none

  real, parameter, private :: k = 3.8317059702

contains

subroutine tilting_cylinder_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz
#ifdef USEST
  real :: xout, zout
#endif

  call get_bounding_box_size(alx, alz)
!  alx = 5.
!  alz = 4.

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)
#ifdef USEST
     if(igeometry.eq.1) then
       call physical_geometry(xout, zout, x, phi, z) 
       x = xout
       z = zout
     end if
#endif
     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call get_local_vals(l)

     call cylinder_equ(x, z)
     call cylinder_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine tilting_cylinder_init

subroutine cylinder_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: rr, ri, arg, befo, ff, fp, fpp, j0, j1, kb
  real :: dbesj0, dbesj1
  
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
     befo = 2./dbesj0(dble(k))
     j0 = dbesj0(dble(arg))
     j1 = dbesj1(dble(arg))
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
     call constant_field(pe0_l, p0-pi0)
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
        kb = k**2*beta*pefac
        pe0_l(1) = 0.5*kb*psi0_l(1)**2 + p0 - pi0
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
        call constant_field(pe0_l, p0-pi0)
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

!=============================================================================
! Wave Propagation (itaylor = 4)
!=============================================================================
module wave_propagation

  implicit none
  real, private :: alx, alz, akx, akx2, omega
  real, private :: psiper, phiper, bzper, vzper, peper, chiper, nper, pper

contains

subroutine wave_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z
  real :: b2,a2
  real :: kp,km,t1,t2,t3
  real :: coef(4)
  real :: root(3)
  real :: error(3)
  real :: bi

  call get_bounding_box_size(alx, alz)

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

  akx = twopi/alx      
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

     if(myrank.eq.0 .and. iprint.ge.1) then
        write(*,*) "Coefs: ", coef(1), coef(2), coef(3), coef(4)
     endif

     call cubic_roots(coef, root, error)

     if(myrank.eq.0 .and. iprint.ge.1) then
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

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call get_local_vals(l)

     call wave_equ(x, z)
     call wave_per(x, z)

     call set_local_vals(l)
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
  call constant_field( pe0_l, p0-pi0)
  call constant_field(  p0_l, p0)
  call constant_field(den0_l, 1.)

end subroutine wave_equ


subroutine wave_per(x, z)
  use math
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  call plane_wave(u1_l, x, z,  akx, 0., phiper, pi/2.)
  call plane_wave(vz1_l, x, z, akx, 0., vzper, pi/2.)
  call plane_wave(chi1_l, x, z, akx, 0., chiper, pi)
  
  call plane_wave(psi1_l, x, z, akx, 0., psiper, pi/2.)
  call plane_wave( bz1_l, x, z, akx, 0., bzper, pi/2.)

  call plane_wave(pe1_l, x, z, akx, 0., peper, pi/2.)
  call plane_wave( p1_l, x, z, akx, 0.,  pper, pi/2.)

  if(idens.eq.1) then
     call plane_wave(den1_l, x, z, akx, 0., nper, pi/2.)    
  endif

end subroutine wave_per

end module wave_propagation

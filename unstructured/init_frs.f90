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

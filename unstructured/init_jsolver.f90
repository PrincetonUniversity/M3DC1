!==============================================================================
! Jsolver_eq
! ~~~~~~~~~~
!
! Loads jsolver equilibrium
!==============================================================================
module jsolver_eq

contains

subroutine jsolver_init()
  use basic
  use arrays
  use gradshafranov
  use newvar_mod
  use sparse
  use coils
  use diagnostics
  use jsolver
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, pzero_jsv, gzero_jsv
  real, allocatable :: ffprime(:),ppxx_jsv2(:),gpx_jsv2(:)

  if(myrank.eq.0 .and. iprint .ge. 1) print *, "jsolver_init called"

  call load_jsolver

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     call jsolver_equ(x, z)
     call jsolver_per(x, z)

     call set_local_vals(l)
  enddo

  xmag = x_jsv(1,1)
  zmag = z_jsv(1,1)

if(myrank.eq.0 .and. iprint .ge. 1) print *, "xmag,zmag=",xmag,zmag

  ! remove factor of 1/xzero_jsv in definition of g
  gxx_jsv = gxx_jsv*xzero_jsv
  gpx_jsv = gpx_jsv*xzero_jsv
  pzero_jsv = 1.5*p_jsv(1)   - .5*p_jsv(2)
  gzero_jsv = 1.5*gxx_jsv(1) - .5*gxx_jsv(2)
!
!.calculate ppxx_jsv and gpx_jsv another way for comparison
  allocate(ppxx_jsv2(npsi_jsv),gpx_jsv2(npsi_jsv))
  do l=2,npsi_jsv-1
     ppxx_jsv2(l) =  (p_jsv(l) - p_jsv(l-1)) /(psival_jsv(l) - psival_jsv(l-1))
     gpx_jsv2(l) = (gxx_jsv(l) - gxx_jsv(l-1))/(psival_jsv(l) - psival_jsv(l-1))
  enddo
! shift gxx and p to be defined at psi_j locations
! instead of at psi_{j+1/2} locations
  do l=npsi_jsv,2,-1
     p_jsv(l) = (p_jsv(l-1) + p_jsv(l))/2.
     gxx_jsv(l) = (gxx_jsv(l-1) + gxx_jsv(l))/2.
  end do
  p_jsv(1)   = pzero_jsv
  gxx_jsv(1) = gzero_jsv

!...Adopt JSOLVER convention that R*B_T = xzero_jsv
  rzero = xzero_jsv
  bzero = 1.0
  if(myrank.eq.0 .and. iprint .ge. 1) print *, "rzero,bzero=",rzero,bzero

  if(igs.gt.0) then
     if(iread_jsolver.eq.2) then
        call default_profiles
     else
        allocate(ffprime(npsi_jsv))
        ffprime = gpx_jsv*gxx_jsv
        call create_profile(npsi_jsv,p_jsv(1:npsi_jsv),ppxx_jsv(1:npsi_jsv), &
             gxx_jsv(1:npsi_jsv),ffprime,psival_jsv(1:npsi_jsv))
!
      if(myrank.eq.0 .and. iprint.ge.1) then
        open(unit=77,file="debug-out",status="unknown")
        write(77,2010) xmag,zmag,tcuro
  2010 format("xmag,zmag,tcuro =",1p3e12.4,/,  &
     "i   press       pprime      fpol        ffprim      flux")
        do l=1,npsi_jsv
        write(77,2011) l,p_jsv(l),ppxx_jsv(l),gxx_jsv(l),ffprime(l),psival_jsv(l)
        enddo
  2011 format(i3,1p5e12.4)
       close(77)
      endif
!
        deallocate(ffprime)
     end if
     if(sigma0.eq.0) then
        call deltafun(xmag,zmag,tcuro,jphi_field)
     else
        call gaussianfun(xmag,zmag,tcuro,sigma0,jphi_field)
     endif
     psibound = psival_jsv(npsi_jsv)
     psimin = psival_jsv(1)

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

  integer :: i, j
  real :: i0, j0, d, dmin, dj
  real, dimension(4) :: a

  ! determine jsolver index (i0,j0) of this node
  i0 = 3
  j0 = 1
  dmin = (x - x_jsv(3,1))**2 + (z - z_jsv(3,1))**2
  do i=3, nthe+2
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
  dj = j - j0
  call cubic_interpolation_coeffs(psival_jsv,npsi_jsv,j,a)
  psi0_l(1) = a(1) + a(2)*dj + a(3)*dj**2 + a(4)*dj**3

  call cubic_interpolation_coeffs(gxx_jsv,npsi_jsv,j,a)
  bz0_l(1) = a(1) + a(2)*dj + a(3)*dj**2 + a(4)*dj**3

  call cubic_interpolation_coeffs(p_jsv,npsi_jsv,j,a)
  p0_l(1) = a(1) + a(2)*dj + a(3)*dj**2 + a(4)*dj**3
  
  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  if(pedge.ge.0.) p0_l = p0_l + pedge

  ! Set electron pressure and density
  pe0_l = pefac*p0_l


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

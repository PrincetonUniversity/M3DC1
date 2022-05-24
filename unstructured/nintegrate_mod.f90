!======================================================================
! nintegrate
! ~~~~~~~~~~
! This module defines general procedures for defining gaussian
! integration quadratures over interior and surface domains
!======================================================================
module nintegrate

  use element

implicit none

! NOTE: All element-specific variables should be declared OMP THREADPRIVATE

  integer :: npoints        ! number of points in Gaussian quadrature
  integer :: npoints_pol
  integer :: npoints_tor
  logical :: surface_int
!$OMP THREADPRIVATE(npoints, npoints_pol, npoints_tor, surface_int)

  real, dimension(MAX_PTS) :: x_79, phi_79, z_79
!$OMP THREADPRIVATE(x_79,phi_79,z_79)
  real, dimension(MAX_PTS) :: xi_79, zi_79, eta_79, weight_79
!$OMP THREADPRIVATE(xi_79,zi_79,eta_79,weight_79)
  vectype, dimension(MAX_PTS,2) :: norm79
!$OMP THREADPRIVATE(norm79)

  real, private, dimension(12) :: alpha_12, beta_12, gamma_12, area_weight_12
  real, private, dimension(25) :: alpha_25, beta_25, gamma_25, area_weight_25
  real, private, dimension(79) :: alpha_79, beta_79, gamma_79, area_weight_79

  real, private, dimension(2) :: delta_2, line_weight_2
  real, private, dimension(3) :: delta_3, line_weight_3
  real, private, dimension(5) :: delta_5, line_weight_5

data delta_2        / -0.577350269, 0.577350269 /
data line_weight_2  /  1.,          1.          /

data delta_3        / -0.774596669, 0.,            0.774596669 /
data line_weight_3  /  0.555555556, 0.888888889  , 0.555555556 /

data delta_5 &
     / -0.906179846, -0.53846931, 0.,         0.53846931, 0.906179846 /
data line_weight_5 &
     /  0.236926885,  0.47862867, 0.56888889, 0.47862867, 0.236926885 /

data alpha_12 &
/ 0.501426509658179, 0.249286745170910, 0.249286745170910, 0.873821971016996, &
  0.063089014491502, 0.063089014491502, 0.053145049844817, 0.310352351033784, &
  0.636502499121399, 0.053145049844817, 0.310352351033784, 0.636502499121399 /

data beta_12 &
/ 0.249286745170910, 0.501426509658179, 0.249286745170910, 0.063089014491502, &
  0.873821971016996, 0.063089014491502, 0.310352351033784, 0.636502499121399, &
  0.053145049844817, 0.636502499121399, 0.053145049844817, 0.310352351033784 /

data gamma_12 &
/ 0.249286745170910, 0.249286745170910, 0.501426509658179, 0.063089014491502, &
  0.063089014491502, 0.873821971016996, 0.636502499121399, 0.053145049844817, &
  0.310352351033784, 0.310352351033784, 0.636502499121399, 0.053145049844817 /

data area_weight_12 &
/ 0.116786275726379, 0.116786275726379, 0.116786275726379, 0.050844906370207, &
  0.050844906370207, 0.050844906370207, 0.082851075618374, 0.082851075618374, &
  0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374/


data alpha_25 &
/ 0.333333333333333, 0.028844733232685, 0.485577633383657, 0.485577633383657, &
  0.781036849029926, 0.109481575485037, 0.109481575485037, 0.141707219414880, &
  0.307939838764121, 0.550352941820999, 0.307939838764121, 0.141707219414880, &
  0.550352941820999, 0.025003534762686, 0.246672560639903, 0.728323904597411, &
  0.246672560639903, 0.025003534762686, 0.728323904597411, 0.009540815400299, &
  0.066803251012200, 0.923655933587500, 0.066803251012200, 0.009540815400299, &
  0.923655933587500 /

data beta_25 &
/ 0.333333333333333, 0.485577633383657, 0.485577633383657, 0.028844733232685, &
  0.109481575485037, 0.109481575485037, 0.781036849029926, 0.307939838764121, &
  0.550352941820999, 0.141707219414880, 0.141707219414880, 0.550352941820999, &
  0.307939838764121, 0.246672560639903, 0.728323904597411, 0.025003534762686, &
  0.025003534762686, 0.728323904597411, 0.246672560639903, 0.066803251012200, &
  0.923655933587500, 0.009540815400299, 0.009540815400299, 0.923655933587500, &
  0.066803251012200 /

data gamma_25 &
/ 0.333333333333333, 0.485577633383657, 0.028844733232685, 0.485577633383657, &
  0.109481575485037, 0.781036849029926, 0.109481575485037, 0.550352941820999, &
  0.141707219414880, 0.307939838764121, 0.550352941820999, 0.307939838764121, &
  0.141707219414880, 0.728323904597411, 0.025003534762686, 0.246672560639903, &
  0.728323904597411, 0.246672560639903, 0.025003534762686, 0.923655933587500, &
  0.009540815400299, 0.066803251012200, 0.923655933587500, 0.066803251012200, &
  0.009540815400299 /

data area_weight_25 &
/ 0.090817990382754, 0.036725957756467, 0.036725957756467, 0.036725957756467, &
  0.045321059435528, 0.045321059435528, 0.045321059435528, 0.072757916845420, &
  0.072757916845420, 0.072757916845420, 0.072757916845420, 0.072757916845420, &
  0.072757916845420, 0.028327242531057, 0.028327242531057, 0.028327242531057, &
  0.028327242531057, 0.028327242531057, 0.028327242531057, 0.009421666963733, &
  0.009421666963733, 0.009421666963733, 0.009421666963733, 0.009421666963733, &
  0.009421666963733 /


data alpha_79 &
/ 0.333333333333333,-0.001900928704400, 0.500950464352200, 0.500950464352200, &
  0.023574084130543, 0.488212957934729, 0.488212957934729, 0.089726636099435, &
  0.455136681950283, 0.455136681950283, 0.196007481363421, 0.401996259318289, &
  0.401996259318289, 0.488214180481157, 0.255892909759421, 0.255892909759421, &
  0.647023488009788, 0.176488255995106, 0.176488255995106, 0.791658289326483, &
  0.104170855336758, 0.104170855336758, 0.893862072318140, 0.053068963840930, &
  0.053068963840930, 0.916762569607942, 0.041618715196029, 0.041618715196029, &
  0.976836157186356, 0.011581921406822, 0.011581921406822, 0.048741583664839, &
  0.606402646106160, 0.344855770229001, 0.606402646106160, 0.048741583664839, &
  0.344855770229001, 0.006314115948605, 0.615842614456541, 0.377843269594854, &
  0.615842614456541, 0.006314115948605, 0.377843269594854, 0.134316520547348, &
  0.559048000390295, 0.306635479062357, 0.559048000390295, 0.134316520547348, &
  0.306635479062357, 0.013973893962392, 0.736606743262866, 0.249419362774742, &
  0.736606743262866, 0.013973893962392, 0.249419362774742, 0.075549132909764, &
  0.711675142287434, 0.212775724802802, 0.711675142287434, 0.075549132909764, &
  0.212775724802802,-0.008368153208227, 0.861402717154987, 0.146965436053239, &
  0.861402717154987,-0.008368153208227, 0.146965436053239, 0.026686063258714, &
  0.835586957912363, 0.137726978828923, 0.835586957912363, 0.026686063258714, &
  0.137726978828923, 0.010547719294141, 0.929756171556853, 0.059696109149007, &
  0.929756171556853, 0.010547719294141, 0.059696109149007 /

data beta_79 &
/ 0.333333333333333, 0.500950464352200,-0.001900928704400, 0.500950464352200, &
  0.488212957934729, 0.023574084130543, 0.488212957934729, 0.455136681950283, &
  0.089726636099435, 0.455136681950283, 0.401996259318289, 0.196007481363421, &
  0.401996259318289, 0.255892909759421, 0.488214180481157, 0.255892909759421, &
  0.176488255995106, 0.647023488009788, 0.176488255995106, 0.104170855336758, &
  0.791658289326483, 0.104170855336758, 0.053068963840930, 0.893862072318140, &
  0.053068963840930, 0.041618715196029, 0.916762569607942, 0.041618715196029, &
  0.011581921406822, 0.976836157186356, 0.011581921406822, 0.344855770229001, &
  0.048741583664839, 0.606402646106160, 0.344855770229001, 0.606402646106160, &
  0.048741583664839, 0.377843269594854, 0.006314115948605, 0.615842614456541, &
  0.377843269594854, 0.615842614456541, 0.006314115948605, 0.306635479062357, &
  0.134316520547348, 0.559048000390295, 0.306635479062357, 0.559048000390295, &
  0.134316520547348, 0.249419362774742, 0.013973893962392, 0.736606743262866, &
  0.249419362774742, 0.736606743262866, 0.013973893962392, 0.212775724802802, &
  0.075549132909764, 0.711675142287434, 0.212775724802802, 0.711675142287434, &
  0.075549132909764, 0.146965436053239,-0.008368153208227, 0.861402717154987, &
  0.146965436053239, 0.861402717154987,-0.008368153208227, 0.137726978828923, &
  0.026686063258714, 0.835586957912363, 0.137726978828923, 0.835586957912363, &
  0.026686063258714, 0.059696109149007, 0.010547719294141, 0.929756171556853, &
  0.059696109149007, 0.929756171556853, 0.010547719294141 /

data gamma_79 &
/ 0.333333333333333, 0.500950464352200, 0.500950464352200,-0.001900928704400, &
  0.488212957934729, 0.488212957934729, 0.023574084130543, 0.455136681950283, &
  0.455136681950283, 0.089726636099435, 0.401996259318289, 0.401996259318289, &
  0.196007481363421, 0.255892909759421, 0.255892909759421, 0.488214180481157, &
  0.176488255995106, 0.176488255995106, 0.647023488009788, 0.104170855336758, &
  0.104170855336758, 0.791658289326483, 0.053068963840930, 0.053068963840930, &
  0.893862072318140, 0.041618715196029, 0.041618715196029, 0.916762569607942, &
  0.011581921406822, 0.011581921406822, 0.976836157186356, 0.606402646106160, &
  0.344855770229001, 0.048741583664839, 0.048741583664839, 0.344855770229001, &
  0.606402646106160, 0.615842614456541, 0.377843269594854, 0.006314115948605, &
  0.006314115948605, 0.377843269594854, 0.615842614456541, 0.559048000390295, &
  0.306635479062357, 0.134316520547348, 0.134316520547348, 0.306635479062357, &
  0.559048000390295, 0.736606743262866, 0.249419362774742, 0.013973893962392, &
  0.013973893962392, 0.249419362774742, 0.736606743262866, 0.711675142287434, &
  0.212775724802802, 0.075549132909764, 0.075549132909764, 0.212775724802802, &
  0.711675142287434, 0.861402717154987, 0.146965436053239,-0.008368153208227, &
 -0.008368153208227, 0.146965436053239, 0.861402717154987, 0.835586957912363, &
  0.137726978828923, 0.026686063258714, 0.026686063258714, 0.137726978828923, &
  0.835586957912363, 0.929756171556853, 0.059696109149007, 0.010547719294141, &
  0.010547719294141, 0.059696109149007, 0.929756171556853 /

data area_weight_79 &
/ 0.033057055541624, 0.000867019185663, 0.000867019185663, 0.000867019185663, &
  0.011660052716448, 0.011660052716448, 0.011660052716448, 0.022876936356421, &
  0.022876936356421, 0.022876936356421, 0.030448982673938, 0.030448982673938, &
  0.030448982673938, 0.030624891725355, 0.030624891725355, 0.030624891725355, &
  0.024368057676800, 0.024368057676800, 0.024368057676800, 0.015997432032024, &
  0.015997432032024, 0.015997432032024, 0.007698301815602, 0.007698301815602, &
  0.007698301815602,-0.000632060497488,-0.000632060497488,-0.000632060497488, &
  0.001751134301193, 0.001751134301193, 0.001751134301193, 0.016465839189576, &
  0.016465839189576, 0.016465839189576, 0.016465839189576, 0.016465839189576, &
  0.016465839189576, 0.004839033540485, 0.004839033540485, 0.004839033540485, &
  0.004839033540485, 0.004839033540485, 0.004839033540485, 0.025804906534650, &
  0.025804906534650, 0.025804906534650, 0.025804906534650, 0.025804906534650, &
  0.025804906534650, 0.008471091054441, 0.008471091054441, 0.008471091054441, &
  0.008471091054441, 0.008471091054441, 0.008471091054441, 0.018354914106280, &
  0.018354914106280, 0.018354914106280, 0.018354914106280, 0.018354914106280, &
  0.018354914106280, 0.000704404677908, 0.000704404677908, 0.000704404677908, &
  0.000704404677908, 0.000704404677908, 0.000704404677908, 0.010112684927462, &
  0.010112684927462, 0.010112684927462, 0.010112684927462, 0.010112684927462, &
  0.010112684927462, 0.003573909385950, 0.003573909385950, 0.003573909385950, &
  0.003573909385950, 0.003573909385950, 0.003573909385950 /

contains

!==================================================
! check_npoints
! ~~~~~~~~~~~~~
!==================================================
subroutine check_npoints()
  implicit none
  if(npoints .gt. MAX_PTS) then
     print *, 'ERROR: npoints > MAX_PTS.  npoints, MAX_PTS = ', &
          npoints, MAX_PTS
     call safestop(13)
  endif
end subroutine check_npoints

!==================================================
! quadrature_implemented
! ~~~~~~~~~~~~~~~~~~~~~~
! returns true if an i-point numerical integration
! quadrature is implemented
!==================================================
logical function quadrature_implemented(i)
  implicit none
  integer, intent(in) :: i
  
  quadrature_implemented = (i.eq.12).or.(i.eq.25).or.(i.eq.79)
  return
end function quadrature_implemented


!==============================================
! edge_to_local
! -------------
!
! Calculates linear transformation from [-1,1] 
! to local coordinates (si, eta).
!==============================================
subroutine edge_to_local(ngauss, delta, line_weight, &
     si1, eta1, si2, eta2, si, eta, local_weight, &
     n1, n2, co, sn)

  implicit none

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: delta, line_weight
  real, intent(in) :: si1, si2, eta1, eta2, co, sn
  real, dimension(ngauss), intent(out) :: si, eta, local_weight
  real, dimension(2), intent(in) :: n1, n2

  real :: l

!!$  real :: m1(2), m2(2)    ! m are the normal vectors in the local coord sys
  real, parameter :: epsilon = 1.-1.e-6

  l = sqrt((si2-si1)**2 + (eta2-eta1)**2)

  si =  0.5*(( si2- si1)*delta +  si2 +  si1)
  eta = 0.5*((eta2-eta1)*delta + eta2 + eta1)
  local_weight = 0.5*line_weight*l

  norm79 = 0.
  norm79(1:ngauss,1) = 0.5*(n2(1)*(1.+delta) + n1(1)*(1.-delta))
  norm79(1:ngauss,2) = 0.5*(n2(2)*(1.+delta) + n1(2)*(1.-delta))

!!$  ! calculate expected normal vector (in global coordinates)
!!$  m1(1) = (eta2 - eta1)*co + (si2 - si1)*sn
!!$  m1(2) = (eta2 - eta1)*sn - (si2 - si1)*co
!!$  ! calculate actual normal vector
!!$  m2 = 0.5*(n1 + n2)
!!$
!!$  ! if actual normal vector is opposite sign of expected vector,
!!$  ! flip integration limits (i.e. flip sign of Jacobian)
!!$  if(m1(1)*m2(1) + m1(2)*m2(2) .lt. 0.) then
!!$     write(*,'(A,6f10.4)') 'Flipping edge.', m1,n1,n2
!!$     local_weight = -local_weight
!!$  endif
 
end subroutine edge_to_local


!==============================================
! area_to_local
! -------------
!
! Calculates linear transformation from area coordinates (alpha, beta, gamma) 
! to local coordinates (si, eta).
!==============================================
subroutine area_to_local(ngauss, alpha, beta, gamma, area_weight, &
     a, b, c, si, eta, local_weight)

  implicit none

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: alpha, beta, gamma, area_weight
  real, intent(in) :: a, b, c
  real, dimension(ngauss), intent(out) :: si, eta, local_weight

  si = (a+b)*(beta - gamma)/2. + (a-b)*(1.-alpha)/2.
  eta = c*alpha
  local_weight = area_weight*(a+b)*c/2.

end subroutine area_to_local


!======================================================================
! extrude_quadrature
! ~~~~~~~~~~~~~~~~~~
! Extends a 2D npol-point quadrature to 3D 
! by taking external product with ntor-point 1D quadrature
!======================================================================
subroutine extrude_quadrature(d, npol, ntor)
  use math
  use mesh_mod
  implicit none

  real, intent(in) :: d        ! toroidal extent of element
  integer, intent(in) :: npol  ! number of poloidal quadrature points
  integer, intent(in) :: ntor  ! number of toroidal quadrature points

  real, dimension(ntor) :: phi, wt
  integer :: i, j

  select case(ntor)

  ! if ntor==0, axisymmetry is assumed
  case(0)
     npoints_pol = npol
     npoints_tor = 1
     npoints = npol
     zi_79(1:npol) = 0.
     weight_79(1:npol) = weight_79(1:npol)*d
     return

  case(2)
     phi = d*(delta_2 + 1.)/2.
     wt = d*line_weight_2/2.

  case(3)
     phi = d*(delta_3 + 1.)/2.
     wt = d*line_weight_3/2.

  case(5)
     phi = d*(delta_5 + 1.)/2.
     wt = d*line_weight_5/2.
     
  case default
     print *, "Error! ", ntor, "-point line quadrature not defined."
     call safestop(44)
  end select

  ! replicate xi, eta, and weight arrays for subsequent toroidal planes
  j = npol + 1
  do i=2, ntor
     xi_79(j:j+npol-1) = xi_79(1:npol)
     eta_79(j:j+npol-1) = eta_79(1:npol)
     weight_79(j:j+npol-1) = weight_79(1:npol)
     if(surface_int) then
        norm79(j:j+npol-1,:) = norm79(1:npol,:)
     endif
     j = j + npol
  end do

  ! set zi and modify weight functions
  j = 1
  do i=1, ntor
     zi_79(j:j+npol-1) = phi(i)
     weight_79(j:j+npol-1) = weight_79(j:j+npol-1)*wt(i)
     j = j + npol
  end do

  npoints_pol = npol
  npoints_tor = ntor
  npoints = npol*ntor
end subroutine extrude_quadrature


!======================================================================
! define_boundary_quadrature
! ~~~~~~~~~~~~~~~~~~~~~~~~~
! Defines the quadrature over a part of the boundary of element ielm.
! In 2D, this boundary is the edge between vertices iedge and iedge+1 
! In 3D, this boundary is the same edge as in 2D, extruded in the 
! toroidal direction to form a rectangular surface.
!======================================================================
subroutine define_boundary_quadrature(ielm, iedge, npol, ntor, normal, idim)
  use mesh_mod
  
  implicit none

  integer, intent(in) :: ielm     ! which element
  integer, intent(in) :: iedge    ! which edge
  integer, intent(in) :: npol     ! number of sampling points in poloidal dir
  integer, intent(in) :: ntor     ! number of sampling points in toroidal dir
  integer, intent(in) :: idim(3)  ! dim of each node.
                                  ! 2=interior, 1=smooth bdry, 0=cusp
  real, intent(in), dimension(2,3) :: normal  ! outward normal vector at
                                              ! each vertex in ielm

  type(element_data) :: d
  real :: si1, si2, eta1, eta2, n1(2), n2(2)

  call get_element_data(ielm,d)
  select case(iedge)
  case(1)
     si1 = -d%b; eta1 = 0.
     si2 =  d%a; eta2 = 0.
  case(2)
     si1 =  d%a; eta1 = 0.
     si2 =  0.;  eta2 = d%c
  case(3)
     si1 =  0.;  eta1 = d%c
     si2 = -d%b; eta2 = 0.    
  case default
     print *, "Error: invalid iedge. ", iedge
     call safestop(45)
  end select
  n1 = normal(:,iedge)
  n2 = normal(:,mod(iedge,3)+1)

  ! Set normal vector of corner nodes
  ! equal to the the normal vector of its adjacent node on this edge.
  if(idim(iedge).eq.0 .and. idim(mod(iedge,3)+1).eq.0) then
     print *, 'Dropping double cusp'
     n1 = 0.
     n2 = 0.
  else
     if(idim(iedge).eq.0) n1 = n2
     if(idim(mod(iedge,3)+1).eq.0) n2 = n1
  end if

  select case(npol)
  case(2)
     call edge_to_local(npol, delta_2, line_weight_2, &
     si1, eta1, si2, eta2, xi_79, eta_79, weight_79, n1, n2, d%co, d%sn)

  case(3)
     call edge_to_local(npol, delta_3, line_weight_3, &
     si1, eta1, si2, eta2, xi_79, eta_79, weight_79, n1, n2, d%co, d%sn)

  case(5)
     call edge_to_local(npol, delta_5, line_weight_5, &
     si1, eta1, si2, eta2, xi_79, eta_79, weight_79, n1, n2, d%co, d%sn)
  case default 
     print *, "Error: ", npol, &
          "-point quadrature not defined for line integration"
  end select

  surface_int = .true.

#ifdef USE3D
  call extrude_quadrature(d%d,npol,ntor)
#else
  npoints_pol = npol
  npoints_tor = 1
  npoints = npol
  zi_79 = 0.
#endif
end subroutine define_boundary_quadrature


!=====================================================
! define_element_quadrature
!=====================================================
subroutine define_element_quadrature(ielm, pol_gauss, tor_gauss)
  use mesh_mod

  integer, intent(in) :: ielm, pol_gauss, tor_gauss

  type(element_data) :: d

  call get_element_data(ielm,d)
  select case(pol_gauss)
  case(12)
     call area_to_local(12,                          &
          alpha_12,beta_12,gamma_12,area_weight_12,  &
          d%a, d%b, d%c, xi_79, eta_79, weight_79)
  case(25)
     call area_to_local(25,                          &
          alpha_25,beta_25,gamma_25,area_weight_25,  &
          d%a, d%b, d%c, xi_79, eta_79, weight_79)
  case(79)
     call area_to_local(79,                          &
          alpha_79,beta_79,gamma_79,area_weight_79,  &
          d%a, d%b, d%c, xi_79, eta_79, weight_79)
  case default
     print *, "Error! ", pol_gauss, "-point triangle quadrature not defined."
     call safestop(44)
  end select

  surface_int = .false.
#ifdef USE3D
  call extrude_quadrature(d%d,pol_gauss,tor_gauss)
#else
  npoints_pol = pol_gauss
  npoints_tor = 1
  npoints = pol_gauss
  zi_79 = 0.0
#endif
end subroutine define_element_quadrature


!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: int0
pure vectype function int0()

  implicit none

  integer :: k

  int0 = 0.
  do k=1, npoints
     int0 = int0 + weight_79(k)
  enddo
end function int0
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: int1
pure vectype function int1(vari)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari

  integer :: k

  int1 = 0.
  do k=1, npoints
     int1 = int1 + vari(k)*weight_79(k)
  enddo
end function int1
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: int2
pure vectype function int2(vari,varj)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj

  integer :: k

  int2 = 0.
  do k=1, npoints
     int2 = int2 + vari(k)*varj(k)*weight_79(k)
  enddo
end function int2
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: int3
pure vectype function int3(vari,varj,vark)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj, vark

  integer :: k

  int3 = 0.
  do k=1, npoints
     int3 = int3 + vari(k)*varj(k)*vark(k)*weight_79(k)
  enddo
end function int3
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: int4
pure vectype function int4(vari,varj,vark,varl)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj, vark, varl

  integer :: k

  int4 = 0.
  do k=1, npoints
     int4 = int4 + vari(k)*varj(k)*vark(k)*varl(k)*weight_79(k)
  enddo
end function int4
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: int5
pure vectype function int5(vari,varj,vark,varl,varm)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj, vark, varl, varm

  integer :: k

  int5 = 0.
  do k=1, npoints
     int5 = int5 + vari(k)*varj(k)*vark(k)*varl(k)*varm(k)*weight_79(k)
  enddo
end function int5
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx0
pure function intx0()

  implicit none

  vectype, dimension(dofs_per_element) :: intx0

  integer :: k

  intx0 = 0.
  do k=1, npoints
     intx0 = intx0 + weight_79(k)
  enddo
end function intx0
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx1
pure function intx1(vari)

  implicit none

  vectype, dimension(dofs_per_element) :: intx1
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari

  integer :: k

  intx1 = 0.
  do k=1, npoints
     intx1 = intx1 + vari(:,k)*weight_79(k)
  enddo
end function intx1
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx2
pure function intx2(vari,varj)

  implicit none

  vectype, dimension(dofs_per_element) :: intx2
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari
  vectype, dimension(npoints), intent(in) :: varj

  integer :: k

  intx2 = 0.
  do k=1, npoints
     intx2 = intx2 + vari(:,k)*varj(k)*weight_79(k)
  enddo
end function intx2
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx3
pure function intx3(vari,varj,vark)

  implicit none

  vectype, dimension(dofs_per_element) :: intx3
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari
  vectype, dimension(npoints), intent(in) :: varj, vark

  integer :: k

  intx3 = 0.
  do k=1, npoints
     intx3 = intx3 + vari(:,k)*varj(k)*vark(k)*weight_79(k)
  enddo
end function intx3
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx4
pure function intx4(vari,varj,vark,varl)

  implicit none

  vectype, dimension(dofs_per_element) :: intx4
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari
  vectype, dimension(npoints), intent(in) :: varj, vark, varl

  integer :: k

  intx4 = 0.
  do k=1, npoints
     intx4 = intx4 + vari(:,k)*varj(k)*vark(k)*varl(k)*weight_79(k)
  enddo
end function intx4
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx5
pure function intx5(vari,varj,vark,varl,varm)

  implicit none

  vectype, dimension(dofs_per_element) :: intx5
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari
  vectype, dimension(npoints), intent(in) :: varj, vark, varl, varm

  integer :: k

  intx5 = 0.
  do k=1, npoints
     intx5 = intx5 + vari(:,k)*varj(k)*vark(k)*varl(k)*varm(k)*weight_79(k)
  enddo
end function intx5
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx6
pure function intx6(vari,varj,vark,varl,varm,varn)

  implicit none

  vectype, dimension(dofs_per_element) :: intx6
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari
  vectype, dimension(npoints), intent(in) :: varj, vark, varl, varm, varn

  integer :: k

  intx6 = 0.
  do k=1, npoints
     intx6 = intx6 + vari(:,k)*varj(k)*vark(k)*varl(k)*varm(k)*varn(k)*weight_79(k)
  enddo
end function intx6
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intx7
pure function intx7(vari,varj,vark,varl,varm,varn,varo)

  implicit none

  vectype, dimension(dofs_per_element) :: intx7
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari
  vectype, dimension(npoints), intent(in) :: varj, vark, varl, varm, varn, varo

  integer :: k

  intx7 = 0.
  do k=1, npoints
     intx7 = intx7 + vari(:,k)*varj(k)*vark(k)*varl(k)*varm(k)*varn(k)*varo(k)*weight_79(k)
  enddo
end function intx7
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intxx0
pure function intxx0()

  implicit none

  vectype, dimension(dofs_per_element, dofs_per_element) :: intxx0

  integer :: k

  intxx0 = 0.
  do k=1, npoints
     intxx0 = intxx0 + weight_79(k)
  enddo
end function intxx0
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intxx1
pure function intxx1(vari)

  implicit none

  vectype, dimension(dofs_per_element, dofs_per_element) :: intxx1
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari

  integer :: i,k

  intxx1 = 0.
  do k=1, npoints
     do i=1, dofs_per_element
        intxx1(i,:) = intxx1(i,:) + vari(i,k)*weight_79(k)
     end do
  enddo
end function intxx1
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intxx2
pure function intxx2(vari,varj)

  implicit none

  vectype, dimension(dofs_per_element, dofs_per_element) :: intxx2
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari, varj

  integer :: i,j,k

  intxx2 = 0.
  do k=1, npoints
     do i=1, dofs_per_element
        do j=1, dofs_per_element
           intxx2(i,j) = intxx2(i,j) + vari(i,k)*varj(j,k)*weight_79(k)
        end do
     end do
  end do
end function intxx2
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intxx3
pure function intxx3(vari,varj,vark)

  implicit none

  vectype, dimension(dofs_per_element,dofs_per_element) :: intxx3
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari, varj
  vectype, dimension(npoints), intent(in) :: vark

  integer :: i,j,k

  intxx3 = 0.
  do k=1, npoints
     do i=1, dofs_per_element
        do j=1, dofs_per_element
           intxx3(i,j) = intxx3(i,j) + vari(i,k)*varj(j,k)*vark(k)*weight_79(k)
        end do
     end do
  enddo
end function intxx3
!DEC$ ATTRIBUTES FORCEINLINE :: intxx4
!==============================================
pure function intxx4(vari,varj,vark,varl)

  implicit none

  vectype, dimension(dofs_per_element, dofs_per_element) :: intxx4
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari, varj
  vectype, dimension(npoints), intent(in) :: vark, varl

  integer :: i,j,k

  intxx4 = 0.
  do k=1, npoints
     do i=1, dofs_per_element
        do j=1, dofs_per_element
           intxx4(i,j) = intxx4(i,j) + &
                vari(i,k)*varj(j,k)*vark(k)*varl(k)*weight_79(k)
        end do
     end do
  end do
end function intxx4
!==============================================
!DEC$ ATTRIBUTES FORCEINLINE :: intxx5
pure function intxx5(vari,varj,vark,varl,varm)

  implicit none

  vectype, dimension(dofs_per_element, dofs_per_element) :: intxx5
  vectype, dimension(dofs_per_element, npoints), intent(in) :: vari, varj
  vectype, dimension(npoints), intent(in) :: vark, varl, varm

  integer :: i,j,k

  intxx5 = 0.
  do k=1, npoints
     do i=1, dofs_per_element
        do j=1, dofs_per_element
           intxx5(i,j) = intxx5(i,j) + &
                vari(i,k)*varj(j,k)*vark(k)*varl(k)*varm(k)*weight_79(k)
        end do
     end do
  enddo
end function intxx5

end module nintegrate

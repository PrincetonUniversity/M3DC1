module nintegrate_mod

implicit none

integer :: npoints        ! number of points in Gaussian quadrature
logical :: surface_int

integer, parameter :: OP_1    = 1
integer, parameter :: OP_DR   = 2
integer, parameter :: OP_DZ   = 3
integer, parameter :: OP_DRR  = 4
integer, parameter :: OP_DRZ  = 5
integer, parameter :: OP_DZZ  = 6
integer, parameter :: OP_LP   = 7
integer, parameter :: OP_GS   = 8
#ifndef USECOMPLEX
integer, parameter :: OP_NUM  = 8
#else
integer, parameter :: OP_DP    = 9
integer, parameter :: OP_DRP   = 10
integer, parameter :: OP_DZP   = 11
integer, parameter :: OP_DRRP  = 12
integer, parameter :: OP_DRZP  = 13
integer, parameter :: OP_DZZP  = 14
integer, parameter :: OP_LPP   = 15
integer, parameter :: OP_GSP   = 16
integer, parameter :: OP_DPP   = 17
integer, parameter :: OP_DRPP  = 18
integer, parameter :: OP_DZPP  = 19
integer, parameter :: OP_DRRPP = 20
integer, parameter :: OP_DRZPP = 21
integer, parameter :: OP_DZZPP = 22
integer, parameter :: OP_LPPP  = 23
integer, parameter :: OP_GSPP  = 24
integer, parameter :: OP_NUM   = 24
#endif

integer, parameter :: FIELD_PHI =     1
integer, parameter :: FIELD_PSI =     2
integer, parameter :: FIELD_V   =     4
integer, parameter :: FIELD_I   =     8
integer, parameter :: FIELD_CHI =    16
integer, parameter :: FIELD_PE  =    32
integer, parameter :: FIELD_P   =    64
integer, parameter :: FIELD_N   =   128
integer, parameter :: FIELD_J   =   256
integer, parameter :: FIELD_VOR =   512
integer, parameter :: FIELD_COM =  1024
integer, parameter :: FIELD_NI  =  2048
integer, parameter :: FIELD_B2I =  4096
integer, parameter :: FIELD_ETA =  8192
integer, parameter :: FIELD_KAP = 16384
integer, parameter :: FIELD_SIG = 32768
integer, parameter :: FIELD_SRC = 65536
integer, parameter :: FIELD_MU  =131072


vectype, dimension(MAX_PTS, OP_NUM, 18) :: g79
real, dimension(MAX_PTS) :: x_79, z_79
vectype, dimension(MAX_PTS) :: r_79, r2_79, r3_79, &
     ri_79, ri2_79, ri3_79, ri4_79, ri5_79, ri6_79, ri7_79, ri8_79
vectype, dimension(MAX_PTS) :: temp79a, temp79b, temp79c, &
                               temp79d, temp79e, temp79f
vectype, dimension(MAX_PTS, OP_NUM) :: sz79
vectype, dimension(MAX_PTS, OP_NUM) :: tm79, ni79, b2i79, sb179, sb279, sp179
vectype, dimension(MAX_PTS, OP_NUM) :: ps179, bz179, pe179, n179, & 
                                       ph179, vz179, ch179, p179
vectype, dimension(MAX_PTS, OP_NUM) :: pst79, bzt79, pet79, nt79, &
                                       pht79, vzt79, cht79, pt79
vectype, dimension(MAX_PTS, OP_NUM) :: vis79, vic79, vip79
vectype, dimension(MAX_PTS, OP_NUM) :: jt79, cot79, vot79, pit79, &
                                       eta79, sig79, bf79
vectype, dimension(MAX_PTS, OP_NUM) :: kap79, kar79, kax79
vectype, dimension(MAX_PTS, OP_NUM) :: ps079, bz079, pe079, n079, &
                                       ph079, vz079, ch079, p079
vectype, dimension(MAX_PTS, OP_NUM) :: pss79, bzs79, phs79, vzs79, chs79

real, dimension(MAX_PTS) :: si_79, eta_79, weight_79

real, private, dimension(12) :: alpha_12, beta_12, gamma_12, area_weight_12
real, private, dimension(25) :: alpha_25, beta_25, gamma_25, area_weight_25
real, private, dimension(79) :: alpha_79, beta_79, gamma_79, area_weight_79

real, private, dimension(5) :: delta_5, line_weight_5
vectype, dimension(MAX_PTS,2) :: norm79

data delta_5        / -0.906180, -0.538469, 0.,       0.538469, 0.906180 /
data line_weight_5  /  0.236927,  0.478629, 0.568889, 0.478629, 0.236927 /

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
     n1, n2, theta)

  use basic

  implicit none

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: delta, line_weight
  real, intent(in) :: si1, si2, eta1, eta2, theta
  real, dimension(ngauss), intent(out) :: si, eta, local_weight
  real, dimension(2), intent(in) :: n1, n2

  real :: l

  real :: m1(2), m2(2)    ! m are the normal vectors in the local coord sys
  real :: co, sn
  real, parameter :: epsilon = 1.-1.e-6

  l = sqrt((si2-si1)**2 + (eta2-eta1)**2)

  si =  0.5*(( si2- si1)*delta +  si2 +  si1)
  eta = 0.5*((eta2-eta1)*delta + eta2 + eta1)
  local_weight = 0.5*line_weight*l

  if(icurv.eq.0) then
     norm79(1:ngauss,1) = 0.5*(n2(1) + n1(1))
     norm79(1:ngauss,2) = 0.5*(n2(2) + n1(2))
  else
     norm79(1:ngauss,1) = 0.5*(n2(1)*(1.+delta) + n1(1)*(1.-delta))
     norm79(1:ngauss,2) = 0.5*(n2(2)*(1.+delta) + n1(2)*(1.-delta))
  end if

  co = cos(theta)
  sn = sin(theta)

  ! calculate expected normal vector (in global coordinates)
  m1(1) = (eta2 - eta1)*co - (si1 - si2)*sn
  m1(2) = (eta2 - eta1)*sn + (si1 - si2)*co
  ! calculate actual normal vector
  m2 = 0.5*(n1 + n2)

  ! if actual normal vector is opposite sign of expected vector,
  ! flip integration limits (i.e. flip sign of Jacobian)
  if(m1(1)*m2(1) + m1(2)*m2(2) .lt. 0.) then
     write(*,'(A,6f10.4)') 'Flipping edge.', m1,m2
     local_weight = -local_weight
  endif
 
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



!=====================================================
! calcpos
! -------
!
! Calculates global coordinates at each sampling point
!=====================================================
subroutine calcpos(itri,si,eta,ngauss,x,z)

  use basic
  use t_data

  implicit none

  integer, intent(in) :: itri, ngauss
  real, dimension(ngauss), intent(in) :: si, eta
  real, dimension(ngauss), intent(out) :: x, z

  integer, dimension(4) :: nodeids
  integer :: i
  real :: xoff, zoff, b, co, sn

  call nodfac(itri,nodeids)
  call nodcoord(nodeids(1), xoff, zoff)

  b = btri(itri)
  co = cos(ttri(itri))
  sn = sin(ttri(itri))

  do i=1, ngauss
     x(i) = (si(i) + b)*co - eta(i)*sn + xoff
     z(i) = (si(i) + b)*sn + eta(i)*co + zoff
  end do

end subroutine calcpos



!===============================================
! eval_ops
! --------
!
! evaluates linear unitary operators
!===============================================
subroutine eval_ops(avector,si,eta,theta,rinv,ngauss,outarr)

  use basic

  implicit none
      
  integer, intent(in) :: ngauss
  vectype, dimension(20), intent(in) :: avector
  real, dimension(ngauss), intent(in) :: si, eta
  vectype, dimension(ngauss), intent(in) :: rinv
  vectype, dimension(MAX_PTS, OP_NUM), intent(out) :: outarr
  real, intent(in) :: theta

  integer :: k,p,op
  real, dimension(OP_NUM) :: sum
  real :: co, sn, co2, sn2, cosn, temp

  co = cos(theta)
  sn = sin(theta)

  co2 = co*co
  sn2 = sn*sn
  cosn = co*sn

  outarr = 0.
 
  ! calculate the answer at each sampling point
  do k=1,ngauss

     do p=1,20
         
        sum = 0.

        sum(OP_1) = si(k)**mi(p)*eta(k)**ni(p)

        if(mi(p).ge.1) then
           ! d_si terms
           temp = mi(p)*si(k)**(mi(p)-1) * eta(k)**ni(p)
           sum(OP_DR) = sum(OP_DR) + co*temp
           sum(OP_DZ) = sum(OP_DZ) + sn*temp           

           if(mi(p).ge.2) then
              ! d_si^2 terms
              temp = si(k)**(mi(p)-2)*(mi(p)-1)*mi(p) * eta(k)**ni(p)
              sum(OP_DRR) = sum(OP_DRR) + co2*temp
              sum(OP_DZZ) = sum(OP_DZZ) + sn2*temp
              sum(OP_DRZ) = sum(OP_DRZ) + cosn*temp
              sum(OP_LP) = sum(OP_LP) + temp
           endif
        endif
        if(ni(p).ge.1) then
           ! d_eta terms
           temp = si(k)**mi(p) * eta(k)**(ni(p)-1)*ni(p)
           sum(OP_DR) = sum(OP_DR) - sn*temp
           sum(OP_DZ) = sum(OP_DZ) + co*temp
           
           if(ni(p).ge.2) then
              ! d_eta^2 terms
              temp = si(k)**mi(p) * eta(k)**(ni(p)-2)*(ni(p)-1)*ni(p)
              sum(OP_DRR) = sum(OP_DRR) + sn2*temp
              sum(OP_DZZ) = sum(OP_DZZ) + co2*temp
              sum(OP_DRZ) = sum(OP_DRZ) - cosn*temp
              sum(OP_LP) = sum(OP_LP) + temp
           endif

           if(mi(p).ge.1) then
              ! d_eta_si terms
              temp = si(k)**(mi(p)-1)*mi(p) * eta(k)**(ni(p)-1)*ni(p)

              sum(OP_DRR) = sum(OP_DRR) - 2.*cosn*temp
              sum(OP_DZZ) = sum(OP_DZZ) + 2.*cosn*temp
              sum(OP_DRZ) = sum(OP_DRZ) + (co2-sn2)*temp
           endif
        endif
               
        ! Grad-Shafranov operator, and
        ! cylindrical correction to Laplacian
        sum(OP_GS) = sum(OP_LP)
        if(itor.eq.1) then
           sum(OP_GS) = sum(OP_GS) - sum(OP_DR)*rinv(k)
           sum(OP_LP) = sum(OP_LP) + sum(OP_DR)*rinv(k)
        endif

        do op=1,OP_NUM
           outarr(k, op) = outarr(k, op) + avector(p)*sum(op)
        enddo

     enddo
  enddo

end subroutine eval_ops


!=====================================================
! define_triangle_quadrature
!=====================================================
subroutine define_triangle_quadrature(itri, ngauss)
  use t_data

  integer, intent(in) :: itri, ngauss

  select case(ngauss)
  case(12)
     call area_to_local(12,                          &
          alpha_12,beta_12,gamma_12,area_weight_12,  &
          atri(itri), btri(itri), ctri(itri),        &
          si_79, eta_79, weight_79)
  case(25)
     call area_to_local(25,                          &
          alpha_25,beta_25,gamma_25,area_weight_25,  &
          atri(itri), btri(itri), ctri(itri),        &
          si_79, eta_79, weight_79)
  case(79)
     call area_to_local(79,                          &
          alpha_79,beta_79,gamma_79,area_weight_79,  &
          atri(itri), btri(itri), ctri(itri),        &
          si_79, eta_79, weight_79)
  case default
     print *, "Error! ", ngauss, "-point triangle quadrature not defined."
     call safestop(44)
  end select
  npoints = ngauss

  surface_int = .false.
end subroutine define_triangle_quadrature


!=====================================================
! define_edge_quadrature
! ~~~~~~~~~~~~~~~~~~~~~~
! 
!=====================================================
subroutine define_edge_quadrature(itri, ivertex, ngauss, normal, idim)
  use t_data 
  integer, intent(in) :: itri, ivertex, ngauss, idim(3)
  real, intent(in), dimension(2,3) :: normal
  real :: si1, si2, eta1, eta2, n1(2), n2(2)

  select case(ivertex)
  case(1)
     si1 = -btri(itri); eta1 = 0.
     si2 =  atri(itri); eta2 = 0.
  case(2)
     si1 =  atri(itri); eta1 = 0.
     si2 =  0.;         eta2 = ctri(itri)
  case(3)
     si1 =  0.;         eta1 = ctri(itri)
     si2 = -btri(itri); eta2 = 0.
  case default
     print *, "Error: invalid ivertex. ", ivertex
     call safestop(45)
  end select
  n1 = normal(:,ivertex)
  n2 = normal(:,mod(ivertex,3)+1)

  ! Set normal vector of corner nodes
  ! equal to the the normal vector of its adjacent node on this edge.
  if(idim(ivertex).eq.0) n1 = n2
  if(idim(mod(ivertex,3)+1).eq.0) n2 = n1

  select case(ngauss)
  case(5)
     call edge_to_local(ngauss, delta_5, line_weight_5, &
     si1, eta1, si2, eta2, si_79, eta_79, weight_79, n1, n2, ttri(itri))
  case default 
     print *, "Error: ", ngauss, &
          "-point quadrature not defined for line integration"
  end select

  npoints = ngauss
  surface_int = .true.
end subroutine define_edge_quadrature


!=====================================================
! define_fields
!=====================================================
subroutine define_fields(itri, fields, gdef, ilin)

  use basic
  use t_data
  use arrays

  implicit none
  
  integer, intent(in) :: itri, fields, gdef, ilin

  integer :: i
  vectype, dimension(20) :: avec

  ! calculate the hyperviscosity coefficients and
  ! the size field for this element.
  hypf = hyper *deex**2
  hypi = hyperi*deex**2
  hypv = hyperv*deex**2
  hypc = hyperc*deex**2
  hypp = hyperp*deex**2
  call interpolate_size_field(itri)

  ! calculate the major radius, and useful powers
  call calcpos(itri, si_79, eta_79, npoints, x_79, z_79)
  if(itor.eq.1) then 
     r_79 = x_79 
  else 
     r_79 = 1.
  endif
  ri_79 = 1./r_79
  ri2_79 = ri_79*ri_79
  ri3_79 = ri2_79*ri_79
  ri4_79 = ri2_79*ri2_79
  ri5_79 = ri3_79*ri2_79
  ri6_79 = ri3_79*ri3_79
  ri7_79 = ri4_79*ri3_79
  ri8_79 = ri4_79*ri4_79
  r2_79 = r_79*r_79
  r3_79 = r2_79*r_79

  if(ijacobian.eq.1) weight_79 = weight_79 * r_79

 
  ! PHI
  ! ~~~
  if(iand(fields, FIELD_PHI).eq.FIELD_PHI) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   U..."
     
     if(ilin.eq.0) then
        call calcavector(itri, field, u_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ph179)
#ifdef USECOMPLEX
        ph179(:,OP_DP :OP_GSP ) = (0,1)*ntor*ph179(:,OP_1:OP_GS)
        ph179(:,OP_DPP:OP_GSPP) =   -ntor**2*ph179(:,OP_1:OP_GS)
#endif
     else
        ph179 = 0.
     endif

     if(eqsubtract.eq.1) then
        call calcavector(itri, field0, u_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ph079)
        phs79 = ph079 + ph179/2.
        pht79 = ph079 + ph179
     else
        pht79 = ph179
        phs79 = ph179/2.
     endif

  endif

  ! PSI
  ! ~~~
  if(iand(fields, FIELD_PSI).eq.FIELD_PSI) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   psi..."

     if(ilin.eq.0) then
        call calcavector(itri, field, psi_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ps179)
#ifdef USECOMPLEX
        ps179(:,OP_DP :OP_GSP ) = (0,1)*ntor*ps179(:,OP_1:OP_GS)
        ps179(:,OP_DPP:OP_GSPP) =   -ntor**2*ps179(:,OP_1:OP_GS)
#endif
     else
        ps179 = 0.
     end if

     if(eqsubtract.eq.1) then
        call calcavector(itri, field0, psi_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ps079)
        pst79 = ps079 + ps179
        pss79 = ps079 + ps179/2.
     else
        pst79 = ps179
        pss79 = ps179/2.
     endif
  endif

  ! V
  ! ~
  if(iand(fields, FIELD_V).eq.FIELD_V) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   V..."

     if(ilin.eq.0) then
        call calcavector(itri, field, vz_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vz179)
#ifdef USECOMPLEX
        vz179(:,OP_DP :OP_GSP ) = (0,1)*ntor*vz179(:,OP_1:OP_GS)
        vz179(:,OP_DPP:OP_GSPP) = -ntor**2*vz179(:,OP_1:OP_GS)
#endif
     else
        vz179 = 0.
     end if

     if(eqsubtract.eq.1) then
        call calcavector(itri, field0, vz_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vz079)
        vzs79 = vz079 + vz179/2.
        vzt79 = vz079 + vz179
     else
        vzt79 = vz179
        vzs79 = vz179/2.
     endif
  endif

  ! I
  ! ~
  if(iand(fields, FIELD_I).eq.FIELD_I) then
     
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   I..."

     if(ilin.eq.0) then
        call calcavector(itri, field, bz_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, bz179)
#ifdef USECOMPLEX
        bz179(:,OP_DP :OP_GSP ) = (0,1)*ntor   *bz179(:,OP_1:OP_GS)
        bz179(:,OP_DPP:OP_GSPP) =      -ntor**2*bz179(:,OP_1:OP_GS)

        call calcavector(itri, bf, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, bf79)
        bf79(:,OP_DP :OP_GSP ) = (0,1)*ntor   *bf79(:,OP_1:OP_GS)
        bf79(:,OP_DPP:OP_GSPP) =      -ntor**2*bf79(:,OP_1:OP_GS)
#endif
     else
        bz179 = 0.
        bf79 = 0.
     endif
       
     if(eqsubtract.eq.1) then
        call calcavector(itri, field0, bz_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, bz079)
        bzt79 = bz079 + bz179
        bzs79 = bz079 + bz179/2.
     else
        bzt79 = bz179
        bzs79 = bz179/2.
     endif
  endif

  ! CHI
  ! ~~~
  if(iand(fields, FIELD_CHI).eq.FIELD_CHI) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   chi..."

     if(ilin.eq.0) then
        call calcavector(itri, field, chi_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ch179)
#ifdef USECOMPLEX
        ch179(:,OP_DP :OP_GSP ) = (0,1)*ntor*ch179(:,OP_1:OP_GS)
        ch179(:,OP_DPP:OP_GSPP) = -ntor**2*ch179(:,OP_1:OP_GS)
#endif
     else
        ch179 = 0.
     end if

     if(eqsubtract.eq.1) then
        call calcavector(itri, field0, chi_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ch079)
        chs79 = ch079 + ch179/2.
        cht79 = ch079 + ch179
     else
        cht79 = ch179
        chs79 = ch179/2.
     endif
  endif

  ! P & PE
  ! ~~~~~~
  if((iand(fields, FIELD_PE).eq.FIELD_PE) .or. &
       (iand(fields, FIELD_P).eq.FIELD_P)) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   P..."

     if(ilin.eq.0) then
        if(ipres.eq.1) then
           call calcavector(itri, field, p_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, p179)
           call calcavector(itri, field, pe_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, pe179)
        else
           call calcavector(itri, field, pe_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, p179)
           pe179 = p179
        endif
#ifdef USECOMPLEX
        p179(:,OP_DP :OP_GSP ) = (0,1)*ntor*p179(:,OP_1:OP_GS)
        p179(:,OP_DPP:OP_GSPP) = -ntor**2*p179(:,OP_1:OP_GS)
        pe179(:,OP_DP :OP_GSP ) = (0,1)*ntor*pe179(:,OP_1:OP_GS)
        pe179(:,OP_DPP:OP_GSPP) = -ntor**2*pe179(:,OP_1:OP_GS)
#endif
     else
        p179 = 0.
        pe179 = 0.
     end if
        
     if(eqsubtract.eq.1) then
        if(ipres.eq.1) then
           call calcavector(itri, field0, p_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, p079)
           call calcavector(itri, field0, pe_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, pe079)
        else
           call calcavector(itri, field0, pe_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, p079)
           pe079 = p079
        endif
        pet79 = pe079 + pe179
        pt79  =  p079 +  p179
     else
        pet79 = pe179
        pt79  =  p179
     endif

     pit79 = pt79 - pefac*pet79
  endif

  
  ! N
  ! ~
  if(iand(fields, FIELD_N).eq.FIELD_N) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   n..."

     if(ilin.eq.0) then
        call calcavector(itri, field, den_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, n179)
#ifdef USECOMPLEX
        n179(:,OP_DP :OP_GSP ) = (0,1)*ntor*n179(:,OP_1:OP_GS)
        n179(:,OP_DPP:OP_GSPP) = -ntor**2*n179(:,OP_1:OP_GS)
#endif    
     else
        n179 = 0.
     end if

     if(eqsubtract.eq.1) then
        if(idenfunc.eq.3) then
           temp79a = (pst79(:,OP_1) - psimin)/(psibound - psimin)
           temp79b = (pst79(:,OP_DR)*(x_79 - xmag) &
                +     pst79(:,OP_DZ)*(z_79 - zmag))*(psibound-psimin)
           n079 = 0.
           where(real(temp79a).lt.denoff .and. real(temp79b).gt.0.)
              n079(:,OP_1) = den0
           elsewhere
              n079(:,OP_1) = den_edge
           end where
        else
           call calcavector(itri, field0, den_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, n079)
        end if
        nt79 = n079 + n179
     else
        nt79 = n179
     endif

  endif

  ! NI
  ! ~~
  if(iand(fields, FIELD_NI).eq.FIELD_NI) then
     ni79(:,OP_1  ) = 1./nt79(:,OP_1)
     ni79(:,OP_DR ) = -ni79(:,OP_1)**2 * nt79(:,OP_DR)
     ni79(:,OP_DZ ) = -ni79(:,OP_1)**2 * nt79(:,OP_DZ)
     ni79(:,OP_DRR) = 2.*ni79(:,OP_1)**3 * nt79(:,OP_DR )**2             &
                      -  ni79(:,OP_1)**2 * nt79(:,OP_DRR)
     ni79(:,OP_DRZ) = 2.*ni79(:,OP_1)**3 * nt79(:,OP_DR )*nt79(:,OP_DZ ) &
                      -  ni79(:,OP_1)**2 * nt79(:,OP_DRZ)
     ni79(:,OP_DZZ) = 2.*ni79(:,OP_1)**3 * nt79(:,OP_DZ )**2             &
                      -  ni79(:,OP_1)**2 * nt79(:,OP_DZZ)
     ni79(:,OP_LP)  = ni79(:,OP_DRR) + ni79(:,OP_DZZ)
     ni79(:,OP_GS)  = ni79(:,OP_LP)
     if(itor.eq.1) then
        ni79(:,OP_LP) = ni79(:,OP_LP) + ri_79*ni79(:,OP_DR)
        ni79(:,OP_GS) = ni79(:,OP_GS) - ri_79*ni79(:,OP_DR)
     endif
#ifdef USECOMPLEX
     ni79(:,OP_DP :OP_GSP ) = 0.
     ni79(:,OP_DPP:OP_GSPP) = 0.
#endif    
  endif
  
  ! J
  ! ~
  if(iand(fields, FIELD_J).eq.FIELD_J) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   j..."

     if(ilin.eq.0) then
        call calcavector(itri, jphi, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, jt79)
     else
        jt79 = 0.
     end if
  endif

  ! VOR
  ! ~~~
  if(iand(fields, FIELD_VOR).eq.FIELD_VOR) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   vor..."

     if(ilin.eq.0) then
        call calcavector(itri, vor, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vot79)
     else
        vot79 = 0.
     end if
  endif

  ! COM
  ! ~~~
  if(iand(fields, FIELD_COM).eq.FIELD_COM) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   com..."

     if(ilin.eq.0) then
        call calcavector(itri, com, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, cot79)
     else
        cot79 = 0.
     end if
  endif


  ! SRC
  ! ~~~
  if(iand(fields, FIELD_SRC).eq.FIELD_SRC) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   sources..."

     call calcavector(itri, sb1, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, sb179)
     if(numvar.ge.2) then
        call calcavector(itri, sb2, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, sb279)
     endif
     if(numvar.ge.3) then
        call calcavector(itri, sp1, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, sp179)
     endif
  endif

  ! B2I
  ! ~~~
  if(iand(fields, FIELD_B2I).eq.FIELD_B2I) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   B^-2..."

     temp79a = ri2_79* &
          (pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2 + bzt79(:,OP_1)**2)

     b2i79(:,OP_1 ) = 1./temp79a
     b2i79(:,OP_DR) = -2.*b2i79(:,OP_1)**2 * ri2_79 * &
          (pst79(:,OP_DR)*pst79(:,OP_DRR) + pst79(:,OP_DZ)*pst79(:,OP_DRZ) &
          +bzt79(:,OP_1 )*bzt79(:,OP_DR ))
     b2i79(:,OP_DZ) = -2.*b2i79(:,OP_1)**2 * ri2_79 * &
          (pst79(:,OP_DR)*pst79(:,OP_DRZ) + pst79(:,OP_DZ)*pst79(:,OP_DZZ) &
          +bzt79(:,OP_1 )*bzt79(:,OP_DZ ))

     if(itor.eq.1) then 
        b2i79(:,OP_DR) = b2i79(:,OP_DR) + 2.*b2i79(:,OP_1)*ri_79
     endif

  endif


  ! ETA
  ! ~~~
  if(iand(fields, FIELD_ETA).eq.FIELD_ETA) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   eta..."

     if(iresfunc.eq.3) then
        temp79a = (pst79(:,OP_1) - psimin)/(psibound - psimin)
        temp79b = (pst79(:,OP_DR)*(x_79 - xmag) &
             +     pst79(:,OP_DZ)*(z_79 - zmag))*(psibound-psimin)
        eta79 = 0.
        where(real(temp79a).lt.etaoff .and. real(temp79b).gt.0.)
           eta79(:,OP_1) = etar
        elsewhere
           eta79(:,OP_1) = eta0
        end where
     else
        call calcavector(itri, resistivity, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, eta79)
     end if
  end if

  ! KAP
  ! ~~~
  if(iand(fields, FIELD_KAP).eq.FIELD_KAP) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   kappa..."

     call calcavector(itri, kappa, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, kap79)

     if(ikapscale.eq.1) then
        kar79 = kappar*kap79
     else
        kar79 = 0.
        kar79(:,OP_1) = kappar
     endif
     kax79 = 0.
     kax79(:,OP_1) = kappax
  end if

  ! SIG
  ! ~~~
  if((iand(fields, FIELD_SIG).eq.FIELD_SIG) .and. idens.eq.1) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   sigma..."

     call calcavector(itri, sigma, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, sig79)
  else
     sig79 = 0.
  end if

  ! MU
  ! ~~
  if(iand(fields, FIELD_MU).eq.FIELD_MU) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   vis..."

     if(ivisfunc.eq.3) then
        temp79a = (pst79(:,OP_1) - psimin)/(psibound - psimin)
        temp79b = (pst79(:,OP_DR)*(x_79 - xmag) &
             +     pst79(:,OP_DZ)*(z_79 - zmag))*(psibound-psimin)
        vis79 = 0.
        vic79 = 0.
        where(real(temp79a).lt.amuoff .and. real(temp79b).gt.0.)
           vis79(:,OP_1) = amu
           vic79(:,OP_1) = amuc
        elsewhere
           vis79(:,OP_1) = amu_edge
           vic79(:,OP_1) = amu_edge
        end where
     else
        call calcavector(itri, visc, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vis79)

        if(numvar.ge.3) then
           call calcavector(itri, visc_c, 1, 1, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vic79)
        endif
     endif

     if(amupar.ne.0.) vip79 = amupar*pit79/2.
  end if

  if(gdef.eq.1) then
     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, &
             ttri(itri), ri_79, npoints, g79(:,:,i))
#ifdef USECOMPLEX
        g79(:,OP_DP :OP_GSP, i) = (0,1)*ntor*g79(:,OP_1:OP_GS,i)
        g79(:,OP_DPP:OP_GSPP,i) =   -ntor**2*g79(:,OP_1:OP_GS,i)
#endif
     end do
  endif
end subroutine define_fields


!==============================================
vectype function int0()

  implicit none

  integer :: k
  vectype :: ksum

  ksum = 0.
  do k=1, npoints
     ksum = ksum + weight_79(k)
  enddo

  int0 = ksum

end function int0
!==============================================
vectype function int1(vari)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari

  integer :: k
  vectype :: ksum

  ksum = 0.
  do k=1, npoints
     ksum = ksum + vari(k)*weight_79(k)
  enddo

  int1 = ksum

end function int1
!==============================================
vectype function int2(vari,varj)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj

  integer :: k
  vectype :: ksum

  ksum = 0.
  do k=1, npoints
     ksum = ksum + vari(k)*varj(k)*weight_79(k)
  enddo

  int2 = ksum

end function int2
!==============================================
vectype function int3(vari,varj,vark)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj, vark

  integer :: k
  vectype :: ksum

  ksum = 0.
  do k=1, npoints
     ksum = ksum + vari(k)*varj(k)*vark(k)*weight_79(k)
  enddo

  int3 = ksum

end function int3
!==============================================
vectype function int4(vari,varj,vark,varl)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj, vark, varl

  integer :: k
  vectype :: ksum

  ksum = 0.
  do k=1, npoints
     ksum = ksum + vari(k)*varj(k)*vark(k)*varl(k)*weight_79(k)
  enddo

  int4 = ksum

end function int4
!==============================================
vectype function int5(vari,varj,vark,varl,varm)

  implicit none

  vectype, dimension(npoints), intent(in) :: vari, varj, vark, varl, varm

  integer :: k
  vectype :: ksum

  ksum = 0.
  do k=1, npoints
     ksum = ksum + vari(k)*varj(k)*vark(k)*varl(k)*varm(k)*weight_79(k)
  enddo

  int5 = ksum

end function int5

subroutine interpolate_size_field(itri)

  use basic
  use t_data

  implicit none

  integer, intent(in) :: itri

  double precision, dimension(3) :: node_sz
  real :: a,b,c,theta,k,l,m,d

  if(ihypdx.eq.0) then
     sz79(:,OP_1  ) = 1.
     sz79(:,OP_DR ) = 0.
     sz79(:,OP_DZ ) = 0.
     sz79(:,OP_DRR) = 0.
     sz79(:,OP_DRZ) = 0.
     sz79(:,OP_DZZ) = 0.
  end if

  a = atri(itri)
  b = btri(itri)
  c = ctri(itri)
  theta = ttri(itri)
  call getelmsizes(itri, node_sz)

  ! use size**2 field
  node_sz = node_sz**ihypdx

  d = b / (a + b)

  m = (node_sz(3) - node_sz(1) - d*(node_sz(2) - node_sz(1))) / c
  l = (node_sz(2) - node_sz(1)) / (a + b)
  k = node_sz(1) + l*b

  sz79(:,OP_1  ) = k + l*si_79 + m*eta_79
  sz79(:,OP_DR ) = k + l*cos(theta) + m*sin(theta)
  sz79(:,OP_DZ ) = k - l*sin(theta) + m*cos(theta)
  sz79(:,OP_DRR) = 0.
  sz79(:,OP_DRZ) = 0.
  sz79(:,OP_DZZ) = 0. 

end subroutine interpolate_size_field

end module nintegrate_mod

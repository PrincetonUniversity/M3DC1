module nintegrate_mod

implicit none

integer :: npoints                 ! number of points in Gaussian quadrature
!integer, parameter :: MAX_PTS = 79 ! maximum number of quad. points allowed

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
integer, parameter :: OP_DPP   = 15
integer, parameter :: OP_DRPP  = 16
integer, parameter :: OP_DZPP  = 17
integer, parameter :: OP_DRRPP = 18
integer, parameter :: OP_DRZPP = 19
integer, parameter :: OP_DZZPP = 20
integer, parameter :: OP_NUM   = 20
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
vectype, dimension(MAX_PTS) :: temp79a, temp79b, temp79c, temp79d, temp79e, temp79f
vectype, dimension(MAX_PTS, OP_NUM) :: sz79
vectype, dimension(MAX_PTS, OP_NUM) :: tm79, ni79, b2i79, sb179, sb279, sp179
vectype, dimension(MAX_PTS, OP_NUM) :: ps179, bz179, pe179, n179, p179, ph179, vz179, ch179
vectype, dimension(MAX_PTS, OP_NUM) :: pst79, bzt79, pet79, nt79, pt79, pht79, vzt79, cht79
vectype, dimension(MAX_PTS, OP_NUM) :: vis79, vic79, vip79
vectype, dimension(MAX_PTS, OP_NUM) :: jt79, cot79, vot79, pit79, eta79, sig79, bf79
vectype, dimension(MAX_PTS, OP_NUM) :: kap79, kar79, kax79
vectype, dimension(MAX_PTS, OP_NUM) :: ps079, bz079, pe079, n079, p079, ph079, vz079, ch079
vectype, dimension(MAX_PTS, OP_NUM) :: pss79, bzs79, phs79, vzs79, chs79

real, dimension(MAX_PTS) :: si_79, eta_79, weight_79

real, dimension(12) :: alpha_12, beta_12, gamma_12, area_weight_12
real, dimension(25) :: alpha_25, beta_25, gamma_25, area_weight_25
real, dimension(79) :: alpha_79, beta_79, gamma_79, area_weight_79

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

  integer :: i

  do i=1, ngauss
     si(i) = (a+b)*(beta(i) - gamma(i))/2. + (a-b)*(1.-alpha(i))/2.
     eta(i) = c*alpha(i)
     local_weight(i) = area_weight(i)*(a+b)*c/2.
  end do

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
! define_fields
!=====================================================
subroutine define_fields(itri, fields, ngauss, gdef)

  use basic
  use t_data
  use arrays

  implicit none
  
  integer, intent(in) :: itri, fields, gdef, ngauss
  
  integer :: i
  vectype, dimension(20) :: avec

  ! calculate the local sampling points and weights for numerical integration
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
     print *, "Error! ", ngauss, "-point quadrature not defined."
     call safestop(44)
  end select
  npoints = ngauss

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
     
     call calcavector(itri, field, u_g, num_fields, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ph179)
#ifdef USECOMPLEX
     ph179(:,OP_DP :OP_DZZP ) = (0,1)*ntor*ph179(:,OP_1:OP_DZZ)
     ph179(:,OP_DPP:OP_DZZPP) =   -ntor**2*ph179(:,OP_1:OP_DZZ)
#endif

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

     call calcavector(itri, field, psi_g, num_fields, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ps179)
#ifdef USECOMPLEX
     ps179(:,OP_DP :OP_DZZP ) = (0,1)*ntor*ps179(:,OP_1:OP_DZZ)
     ps179(:,OP_DPP:OP_DZZPP) =   -ntor**2*ps179(:,OP_1:OP_DZZ)
#endif

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

     call calcavector(itri, field, vz_g, num_fields, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vz179)
#ifdef USECOMPLEX
     vz179(:,OP_DP :OP_DZZP ) = (0,1)*ntor*vz179(:,OP_1:OP_DZZ)
     vz179(:,OP_DPP:OP_DZZPP) = -ntor**2*vz179(:,OP_1:OP_DZZ)
#endif
    
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

     if(numvar.ge.2) then
        
        call calcavector(itri, field, bz_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, bz179)
#ifdef USECOMPLEX
        bz179(:,OP_DP :OP_DZZP ) = (0,1)*ntor   *bz179(:,OP_1:OP_DZZ)
        bz179(:,OP_DPP:OP_DZZPP) =      -ntor**2*bz179(:,OP_1:OP_DZZ)

        call calcavector(itri, bf, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, bf79)
        bf79(:,OP_DP :OP_DZZP ) = (0,1)*ntor   *bf79(:,OP_1:OP_DZZ)
        bf79(:,OP_DPP:OP_DZZPP) =      -ntor**2*bf79(:,OP_1:OP_DZZ)
#endif
       
        if(eqsubtract.eq.1) then
           call calcavector(itri, field0, bz_g, num_fields, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, bz079)
           bzt79 = bz079 + bz179
           bzs79 = bz079 + bz179/2.
        else
           bzt79 = bz179
           bzs79 = bz179/2.
        endif
     else

        bz079 = 0.
        bz179 = 0.
        bzt79 = 0.
        bzs79 = 0.

        ! negative signs adjust for left-handed coordinates.
        bzt79(:,OP_1) = -bzero
        bzs79(:,OP_1) = -bzero/2.
        
        if(eqsubtract.eq.1) then
           bz079(:,OP_1) = -bzero
        else
           bz179(:,OP_1) = -bzero
        endif

        if(itor.eq.1) then
           bzt79(:,OP_1) = bzt79(:,OP_1)*rzero !cjdebug xzero
           bzs79(:,OP_1) = bzs79(:,OP_1)*rzero !cjdebug xzero
           bz079(:,OP_1) = bz079(:,OP_1)*rzero !cjdebug xzero
           bz179(:,OP_1) = bz179(:,OP_1)*rzero !cjdebug xzero
        endif
     endif
  endif

  ! CHI
  ! ~~~
  if(iand(fields, FIELD_CHI).eq.FIELD_CHI) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   chi..."

     call calcavector(itri, field, chi_g, num_fields, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ch179)
#ifdef USECOMPLEX
     ch179(:,OP_DP :OP_DZZP ) = (0,1)*ntor*ch179(:,OP_1:OP_DZZ)
     ch179(:,OP_DPP:OP_DZZPP) = -ntor**2*ch179(:,OP_1:OP_DZZ)
#endif

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
     p179(:,OP_DP :OP_DZZP ) = (0,1)*ntor*p179(:,OP_1:OP_DZZ)
     p179(:,OP_DPP:OP_DZZPP) = -ntor**2*p179(:,OP_1:OP_DZZ)
     pe179(:,OP_DP :OP_DZZP ) = (0,1)*ntor*pe179(:,OP_1:OP_DZZ)
     pe179(:,OP_DPP:OP_DZZPP) = -ntor**2*pe179(:,OP_1:OP_DZZ)
#endif
        
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

     call calcavector(itri, field, den_g, num_fields, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, n179)
#ifdef USECOMPLEX
     n179(:,OP_DP :OP_DZZP ) = (0,1)*ntor*n179(:,OP_1:OP_DZZ)
     n179(:,OP_DPP:OP_DZZPP) = -ntor**2*n179(:,OP_1:OP_DZZ)
#endif    
     if(eqsubtract.eq.1) then
        call calcavector(itri, field0, den_g, num_fields, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, n079)
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
  endif
  
  ! J
  ! ~
  if(iand(fields, FIELD_J).eq.FIELD_J) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   j..."

     call calcavector(itri, jphi, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, jt79)
  endif

  ! VOR
  ! ~~~
  if(iand(fields, FIELD_VOR).eq.FIELD_VOR) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   vor..."

     call calcavector(itri, vor, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vot79)
  endif

  ! COM
  ! ~~~
  if(iand(fields, FIELD_COM).eq.FIELD_COM) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.1) print *, "   com..."

     call calcavector(itri, com, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, cot79)
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

     call calcavector(itri, resistivity, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, eta79)
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

     call calcavector(itri, visc, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vis79)

     if(numvar.ge.3) then
        call calcavector(itri, visc_c, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, vic79)
     endif

     if(amupar.ne.0.) vip79 = amupar*pit79/2.
  end if

  if(gdef.eq.1) then
     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, &
             ttri(itri), ri_79, npoints, g79(:,:,i))
#ifdef USECOMPLEX
        g79(:,OP_DP :OP_DZZP, i) = (0,1)*ntor*g79(:,OP_1:OP_DZZ,i)
        g79(:,OP_DPP:OP_DZZPP,i) =   -ntor**2*g79(:,OP_1:OP_DZZ,i)
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

end module nintegrate_mod

!==========================================================
subroutine calcavector(itri, inarr, itype, numvare, avector)

  use t_data

  implicit none

  integer, intent(in) :: itri, itype, numvare
  vectype, dimension(*), intent(in) :: inarr
  vectype, dimension(20), intent(out) :: avector

  integer :: ibegin, iendplusone
    
  integer :: i, ii, iii, k
  vectype, dimension(18) :: wlocal

  ! construct the 18 vector corresponding to this triangle
  ! calculate the index and local coordinates for this triangle
  do iii=1,3  
     call entdofs(numvare, ist(itri,iii)+1, 0, ibegin, iendplusone)
     do ii=1,6
        i = (iii-1)*6 + ii
        wlocal(i) = inarr(ibegin+ii-1+(itype-1)*6)
     enddo
  enddo

  ! calculate the function value corresponding to this point
  do i=1,20
     avector(i) = 0.
     do k=1,18
        avector(i) = avector(i) + gtri(i,k,itri)*wlocal(k)
     enddo
  enddo

end subroutine calcavector
!==========================================================
subroutine calcrvector(itri, x1, z1, rvector)

  use t_data
  use basic

  implicit none

  integer, intent(in) :: itri
  real, intent(in) :: x1, z1
  real, dimension(20), intent(out) :: rvector
  
  integer :: i, iii, k
  real, dimension(18) :: rlocal

  ! construct the 6-vector corresponding to 1/r at the nodes of the triangle
  do iii=1,3  
     i = (iii-1)*6 + 1
     rlocal(i)   =  1./(x1+xzero)
     rlocal(i+1) = -1./(x1+xzero)**2
     rlocal(i+2) = 0.
     rlocal(i+3) =  2./(x1+xzero)**3
     rlocal(i+4) = 0.
     rlocal(i+5) = 0.
     rlocal(i+6) = 0.
  enddo

  ! calculate the value of 1/r corresponding to this point
  do i=1,20
     rvector(i) = 0.
     do k=1,18
        rvector(i) = rvector(i) + gtri(i,k,itri)*rlocal(k)
     enddo
  enddo

end subroutine calcrvector

!============================================================
subroutine evaluate(x,z,ans,ans2,dum,itype,numvare,itri)
  
  use p_data
  use t_data
  use basic

  use nintegrate_mod

  implicit none

  include 'mpif.h'

  integer, intent(in) :: itype, numvare
  integer, intent(inout) :: itri
  real, intent(in) :: x, z
  vectype, intent(in) :: dum(*)
  real, intent(out) :: ans, ans2

  integer :: p, nodeids(4), ier
  real :: x1, z1
  vectype, dimension(20) :: avector
  real :: ri, si, eta, co, sn
  real :: term1, term2
  real, dimension(2) :: temp1, temp2
  integer :: hasval, tothasval

  ! evaluate the solution to get the value [ans] at one point (x,z)

  ! first find out what triangle x,z is in.  whattri
  ! returns itri, x1, and z1 with x1 and z1 being
  ! the coordinates of the first node/vertex

  if(itri.eq.0) then
     call whattri(x,z,itri,x1,z1)
  else
     call nodfac(itri,nodeids)
     call nodcoord(nodeids(1), x1, z1)
  endif

  ans = 0.
  ans2 = 0.


  ! If this process contains the point, evaluate the field at that point.
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(itri.gt.0) then

     ! calculate local coordinates
     co = cos(ttri(itri))
     sn = sin(ttri(itri))
  
     si  = (x-x1)*co + (z-z1)*sn - btri(itri)
     eta =-(x-x1)*sn + (z-z1)*co

     ! calculate the inverse radius
     if(itor.eq.1) then
        ri = 1./x
     else
        ri = 1.
     endif

     ! calculate the value of the function
     call calcavector(itri, dum, itype, numvare, avector)
     
     do p=1,20
     
        term1 = si**mi(p)*eta**ni(p)
        term2 = 0.
        
        if(mi(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 - 2.*co*(mi(p)*si**(mi(p)-1) * eta**ni(p))*ri
           endif
           
           if(mi(p).ge.2) then
              term2 = term2 + si**(mi(p)-2)*(mi(p)-1)*mi(p) * eta**ni(p)
           endif
        endif
     
        if(ni(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 + 2.*sn*(si**mi(p) * eta**(ni(p)-1)*ni(p))*ri
           endif
           
           if(ni(p).ge.2) then
              term2 = term2 + si**mi(p) * eta**(ni(p)-2)*(ni(p)-1)*ni(p)
           endif
        endif
     
        ans = ans + avector(p)*term1
        ans2 = ans2 + avector(p)*term2
        hasval = 1
     enddo
  else
     hasval = 0
  endif
     

  ! Distribute the result if necessary
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(maxrank.gt.1) then
     ! Determine number of processes whose domains contain this point
     call mpi_allreduce(hasval, tothasval, 1, MPI_INTEGER, MPI_SUM, &
          MPI_COMM_WORLD, ier)

     ! Find the average value at this point over all processes containing
     ! the point.  (Each value should be identical.)
     temp1(1) = ans
     temp1(2) = ans2
     call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_WORLD, ier)
     ans = temp2(1)/tothasval
     ans2 = temp2(2)/tothasval
  endif

end subroutine evaluate
!============================================================

subroutine interpolate_size_field(itri)

  use basic
  use nintegrate_mod

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

  call getelmparams(itri, a, b, c, theta)
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

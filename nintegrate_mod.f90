module nintegrate_mod

implicit none

integer, parameter :: OP_1    = 1
integer, parameter :: OP_DR   = 2
integer, parameter :: OP_DZ   = 3
integer, parameter :: OP_DRR  = 4
integer, parameter :: OP_DRZ  = 5
integer, parameter :: OP_DZZ  = 6
integer, parameter :: OP_LP   = 7
integer, parameter :: OP_GS   = 8

integer, parameter :: OP_NUM  = 8

real, dimension(25) :: r_25, r2_25, ri_25, ri2_25, ri3_25, ri4_25
real, dimension(25, OP_NUM, 18) :: g25
real, dimension(25, OP_NUM) :: ps025, bz025, pe025, n025, p025, ph025, vz025, ch025
real, dimension(25, OP_NUM) :: ps125, bz125, pe125, n125, p125, ph125, vz125, ch125
real, dimension(25, OP_NUM) :: pst25, bzt25, pet25, nt25, pt25, pht25, vzt25, cht25

real, dimension(79) :: r_79, r2_79, ri_79, ri2_79, ri3_79, ri4_79, ri5_79, ri6_79, ri7_79
real, dimension(79, OP_NUM, 18) :: g79
real, dimension(79, OP_NUM) :: tm79, ni79, b2i79
real, dimension(79, OP_NUM) :: ps079, bz079, pe079, n079, p079, ph079, vz079, ch079
real, dimension(79, OP_NUM) :: ps179, bz179, pe179, n179, p179, ph179, vz179, ch179
real, dimension(79, OP_NUM) :: pst79, bzt79, pet79, nt79, pt79, pht79, vzt79, cht79
real, dimension(79, OP_NUM) :: pss79, bzs79, phs79, vzs79, chs79
real, dimension(79, OP_NUM) :: sb179, sb279
real, dimension(79) :: temp79a, temp79b, temp79c, temp79d, temp79e, temp79f

real, dimension(25) :: si_25, eta_25, weight_25
real, dimension(79) :: si_79, eta_79, weight_79

real, dimension(25) ::  alpha_25, beta_25, gamma_25, area_weight_25
real, dimension(79) ::  alpha_79, beta_79, gamma_79, area_weight_79

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
     eta(i) = 2.*c*alpha(i)/2.
     local_weight(i) = area_weight(i)*(a+b)*c/2.
  end do

end subroutine area_to_local


!============================================
! calcr
! -----
!
! Calculates r at each sampling point
!============================================
subroutine calcr(itri,si,eta,ngauss,r)

  use basic
  use t_data

  implicit none

  integer, intent(in) :: itri, ngauss
  real, dimension(ngauss), intent(in) :: si, eta
  real, dimension(ngauss), intent(out) :: r

  integer, dimension(4) :: nodeids(4)
  integer :: i
  real :: xoff, b, co, sn, xmin, zmin
  double precision :: coords(3)

  if(itor.eq.1) then

     call getmincoord(xmin,zmin)
     call nodfac(itri,nodeids)
     call xyznod(nodeids(1), coords)

     xoff = coords(1) - xmin + xzero

     b = btri(itri)
     co = cos(ttri(itri))
     sn = sin(ttri(itri))

     do i=1, ngauss
        r(i) = (si(i) + b)*co - eta(i)*sn + xoff
     end do
  else
     r = 1.
  endif

end subroutine calcr


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
  real, dimension(20), intent(in) :: avector
  real, dimension(ngauss), intent(in) :: si, eta, rinv
  real, dimension(ngauss, OP_NUM), intent(out) :: outarr
  real, intent(in) :: theta

  integer :: k,p,op
  real, dimension(OP_NUM) :: sum
  real :: co, sn, co2, sn2, cosn, temp

!!$  if(itor.eq.1) then
     co = cos(theta)
     sn = sin(theta)
!!$  else
!!$     co = 1.
!!$     sn = 0.
!!$  endif
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
           sum(OP_GS) = sum(OP_GS) - 2.*sum(OP_DR)*rinv(k)
           sum(OP_LP) = sum(OP_LP) +    sum(OP_DR)*rinv(k)
        endif

        do op=1,OP_NUM
           outarr(k, op) = outarr(k, op) + avector(p)*sum(op)
        enddo

     enddo
  enddo

end subroutine eval_ops

!=====================================================
! define_fields_25
!=====================================================
subroutine define_fields_25(itri)

  use basic
  use t_data
  use arrays

  implicit none
 
  integer, intent(in) :: itri
  integer :: i
  real, dimension(20) :: avec

  ! calculate the local sampling points and weights for numerical integration
  call area_to_local(25,                                            &
       alpha_25,beta_25,gamma_25,area_weight_25,                    &
       atri(itri), btri(itri), ctri(itri),                          &
       si_25, eta_25, weight_25)

  call calcr(itri, si_25, eta_25, 25, r_25)
  ri_25 = 1./r_25
  ri2_25 = ri_25*ri_25
  ri3_25 = ri2_25*ri_25
  ri4_25 = ri2_25*ri2_25
  r2_25 = r_25*r_25

  call calcavector(itri, vel, 1, numvar, avec)
  call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, ph125)
  call calcavector(itri, phi, 1, numvar, avec)
  call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, ps125)

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     call calcavector(itri, vel0, 1, numvar, avec)
     call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, ph025)
     call calcavector(itri, phi0, 1, numvar, avec)
     call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, ps025)
     pht25 = ph025 + ph125
     pst25 = ps025 + ps125
  else
     pht25 = ph125
     pst25 = ps125
  endif

  if(numvar.ge.2) then
     call calcavector(itri, vel, 2, numvar, avec)
     call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, vz125)
     call calcavector(itri, phi, 2, numvar, avec)
     call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, bz125)
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, vel0, 2, numvar, avec)
        call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, vz025)
        call calcavector(itri, phi0, 2, numvar, avec)
        call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, bz025)
        vzt25 = vz025 + vz125
        bzt25 = bz025 + bz125
     else
        vzt25 = vz125
        bzt25 = bz125
     endif
    
     if(numvar.ge.3) then
        call calcavector(itri, vel, 3, numvar, avec)
        call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, ch125)

        if(ipres.eq.1) then
           call calcavector(itri, pres, 1, 1, avec)
           call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, p125)
           call calcavector(itri, phi, 3, numvar, avec)
           call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, pe125)
        else
           call calcavector(itri, phi, 3, numvar, avec)
           call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, p125)
!           pe125 = ((p0-pi0)/p0)*p125
           pe125 = p125
        endif
           
        if(linear.eq.1 .or. eqsubtract.eq.1) then
           call calcavector(itri, vel0, 3, numvar, avec)
           call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, ch025)

           if(ipres.eq.1) then
              call calcavector(itri, pres0, 1, 1, avec)
              call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, p025)
              call calcavector(itri, phi0, 3, numvar, avec)
              call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, pe025)
           else
              call calcavector(itri, phi0, 3, numvar, avec)
              call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, p025)
!              pe025 = ((p0-pi0)/p0)*p025
              pe025 = p025
           endif

           cht25 = ch025 + ch125
           pet25 = pe025 + pe125
           pt25  =  p025 +  p125
        else
           cht25 = ch125
           pet25 = pe125
           pt25  =  p125
        endif
     endif

  endif
  
  if(idens.eq.1) then
     call calcavector(itri, den, 1, 1, avec)
     call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, n125)
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, den0, 1, 1, avec)
        call eval_ops(avec, si_25, eta_25, ttri(itri), ri_25,25, n025)
        nt25 = n025 + n125
     else
        nt25 = n125
     endif
  endif

  do i=1,18
     call eval_ops(gtri(:,i,itri), si_25, eta_25, ttri(itri), ri_25, 25, g25(:,:,i))
  end do

end subroutine define_fields_25


!=====================================================
! define_fields_79
!=====================================================
subroutine define_fields_79(itri)

  use basic
  use t_data
  use arrays

  implicit none
  
  integer, intent(in) :: itri
  
  integer :: i
  real, dimension(20) :: avec

  ! calculate the local sampling points and weights for numerical integration
  call area_to_local(79,                                            &
       alpha_79,beta_79,gamma_79,area_weight_79,                    &
       atri(itri), btri(itri), ctri(itri),                          &
       si_79, eta_79, weight_79)

  call calcr(itri, si_79, eta_79, 79, r_79)
  ri_79 = 1./r_79
  ri2_79 = ri_79*ri_79
  ri3_79 = ri2_79*ri_79
  ri4_79 = ri2_79*ri2_79
  ri5_79 = ri3_79*ri2_79
  ri6_79 = ri3_79*ri3_79
  ri7_79 = ri4_79*ri3_79
  r2_79 = r_79*r_79

  call calcavector(itri, sb1, 1, 1, avec)
  call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, sb179)
  
  call calcavector(itri, vel, 1, numvar, avec)
  call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ph179)
  call calcavector(itri, phi, 1, numvar, avec)
  call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps179)

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     call calcavector(itri, vel0, 1, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ph079)
     call calcavector(itri, phi0, 1, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps079)
     pht79 = ph079 + ph179
     pst79 = ps079 + ps179
     phs79 = ph079 + ph179/2.
     pss79 = ps079 + ps179/2.
  else
     pht79 = ph179
     pst79 = ps179
     phs79 = ph179/2.
     pss79 = ps179/2.
  endif

  if(numvar.ge.2) then
     call calcavector(itri, sb2, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, sb279)

     call calcavector(itri, vel, 2, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, vz179)
     call calcavector(itri, phi, 2, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz179)
    
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, vel0, 2, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, vz079)
        call calcavector(itri, phi0, 2, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz079)
        vzt79 = vz079 + vz179
        bzt79 = bz079 + bz179
        bzs79 = bz079 + bz179/2.
     else
        vzt79 = vz179
        bzt79 = bz179
        bzs79 = bz179/2.
     endif
    
     if(numvar.ge.3) then
        call calcavector(itri, vel, 3, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ch179)

        if(ipres.eq.1) then
           call calcavector(itri, pres, 1, 1, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, p179)
           call calcavector(itri, phi, 3, numvar, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, pe179)
        else
           call calcavector(itri, phi, 3, numvar, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, p179)
!           pe179 = ((p0-pi0)/p0)*p179
           pe179 = p179
        endif
           
        if(linear.eq.1 .or. eqsubtract.eq.1) then
           call calcavector(itri, vel0, 3, numvar, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ch079)

           if(ipres.eq.1) then
              call calcavector(itri, pres0, 1, 1, avec)
              call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, p079)
              call calcavector(itri, phi0, 3, numvar, avec)
              call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, pe079)
           else
              call calcavector(itri, phi0, 3, numvar, avec)
              call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, p079)
!              pe079 = ((p0-pi0)/p0)*p079
              pe079 = p079
           endif

           cht79 = ch079 + ch179
           pet79 = pe079 + pe179
           pt79  =  p079 +  p179
           chs79 = ch079 + ch179/2.
        else
           cht79 = ch179
           pet79 = pe179
           pt79  =  p179
           chs79 = ch179/2.
        endif
     endif

  endif
  
  if(idens.eq.1) then
     call calcavector(itri, den, 1, 1, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, n179)
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, den0, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, n079)
        nt79 = n079 + n179
     else
        nt79 = n179
     endif

     ! calculate 1/n
     ni79(:,OP_1  ) = 1./nt79(:,OP_1)
     temp79a = -ni79(:,OP_1)**2
     ni79(:,OP_DR ) = temp79a*nt79(:,OP_DR)
     ni79(:,OP_DZ ) = temp79a*nt79(:,OP_DZ)
     temp79b = -2.*temp79a*ni79(:,OP_1)
     ni79(:,OP_DRR) = temp79b*nt79(:,OP_DR)**2           +temp79a*nt79(:,OP_DRR)
     ni79(:,OP_DRZ) = temp79b*nt79(:,OP_DR)*nt79(:,OP_DZ)+temp79a*nt79(:,OP_DRZ)
     ni79(:,OP_DZZ) = temp79b*nt79(:,OP_DZ)**2           +temp79a*nt79(:,OP_DZZ)
  endif

  if(numvar.eq.1) then
     b2i79(:,OP_1) = 1./(pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)
  else
     b2i79(:,OP_1) = 1./(pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2 + bzt79(:,OP_1)**2)
  endif

  do i=1,18
     call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
  end do

end subroutine define_fields_79



!==============================================
real function int1(vari,weight,ngauss)

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: vari, weight

  integer :: k
  real :: ksum

  ksum = 0.
  do k=1, ngauss
     ksum = ksum + vari(k)*weight(k)
  enddo

  int1 = ksum

end function int1
!==============================================
real function int2(vari,varj,weight,ngauss)

  implicit none

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: vari, varj, weight

  integer :: k
  real :: ksum

  ksum = 0.
  do k=1, ngauss
     ksum = ksum + vari(k)*varj(k)*weight(k)
  enddo

  int2 = ksum

end function int2
!==============================================
real function int3(vari,varj,vark, weight,ngauss)

  implicit none

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: vari, varj, vark, weight

  integer :: k
  real :: ksum

  ksum = 0.
  do k=1, ngauss
     ksum = ksum + vari(k)*varj(k)*vark(k)*weight(k)
  enddo

  int3 = ksum

end function int3
!==============================================
real function int4(vari,varj,vark,varl,weight,ngauss)

  use t_data

  implicit none

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: vari, varj, vark, varl,weight

  integer :: k
  real :: ksum

  ksum = 0.
  do k=1, ngauss
     ksum = ksum + vari(k)*varj(k)*vark(k)*varl(k)*weight(k)
  enddo

  int4 = ksum

end function int4
!==============================================
real function int5(vari,varj,vark,varl,varm,weight,ngauss)

  use t_data

  implicit none

  integer, intent(in) :: ngauss
  real, dimension(ngauss), intent(in) :: vari, varj, vark, varl, varm, weight

  integer :: k
  real :: ksum

  ksum = 0.
  do k=1, ngauss
     ksum = ksum + vari(k)*varj(k)*vark(k)*varl(k)*varm(k)*weight(k)
  enddo

  int5 = ksum

end function int5

end module nintegrate_mod

!==========================================================
subroutine calcavector(itri, inarr, itype, numvare, avector)

  use t_data

  implicit none

  integer, intent(in) :: itri, itype, numvare
  real, dimension(*), intent(in) :: inarr
  real, dimension(20), intent(out) :: avector

  integer :: ibegin, iendplusone
    
  integer :: i, ii, iii, k
  real, dimension(18) :: wlocal

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
  integer, intent(in) :: itype, numvare, itri
  real, intent(in) :: x, z, dum(*)
  real, intent(out) :: ans, ans2

  integer :: p, nodeids(4)
  real :: x1, z1, xmin, zmin
  real, dimension(20) :: avector
  real :: ri, si, eta, co, sn
  real :: term1, term2

  double precision :: coords(3)

  ! itype=1 for first velocity (flux) variable

  ! evaluate the solution to get the value [ans] at one point (x,z)

  ! first find out what triangle x,z is in.  whattri
  ! returns itri, x1, and z1 with x1 and z1 being
  ! the coordinates of the first node/vertex

  if(itri.eq.0) then
     call whattri(x,z,itri,x1,z1)
  else
     call nodfac(itri,nodeids)
     call xyznod(nodeids(1), coords)
     x1 = coords(1)
     z1 = coords(2)
  endif

  call getmincoord(xmin, zmin)

  ! calculate local coordinates

  co = cos(ttri(itri))
  sn = sin(ttri(itri))

  si  = (x-x1)*co + (z-z1)*sn - btri(itri)
  eta =-(x-x1)*sn + (z-z1)*co

  ! calculate the inverse radius
  if(itor.eq.1) then
     ri = 1./(x - xmin + xzero)
  else
     ri =1.
  endif

  ! calculate the value of the function
  call calcavector(itri, dum, itype, numvare, avector)

  ans = 0.
  ans2 = 0.

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
  enddo

end subroutine evaluate
!============================================================

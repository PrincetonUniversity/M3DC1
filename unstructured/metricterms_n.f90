module metricterms_n

implicit none

contains

!===============================================================================
! V1 TERMS
!===============================================================================


! V1umu
! =====
real function v1umu(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = -int2(e(:,OP_GS),f(:,OP_LP),weight_79,79)

  v1umu = temp
  return
end function v1umu


! V1un
! ====
real function v1un(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  if(idens.eq.0) then
     temp = int2(e(:,OP_DR),f(:,OP_DR),weight_79,79) &
          + int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79)
  else
     temp = int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  endif

  v1un = temp
  return
end function v1un


! V1chin
! ======
real function v1chin(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  if(idens.eq.0) then
     temp = 0.
  else
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79)
  endif

  v1chin = temp
  return
end function v1chin


! V1uun
! =====
real function v1uun(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp
           
  if(idens.eq.0) then

     temp = int4(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_LP),weight_79,79) &
          - int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_LP),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + int3(e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_LP),weight_79,79) &
             + int3(e(:,OP_1 ),f(:,OP_DZZ),g(:,OP_DZ),weight_79,79) &
             + int3(e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DR),weight_79,79) &
             + int3(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZ),weight_79,79) &
             + int3(e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DR),weight_79,79)
     endif
  else

     temp = -int5(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DRR),h(:,OP_DR),weight_79,79) &
          +  int5(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DRZ),h(:,OP_DR),weight_79,79) &
          -  int5(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DRZ),h(:,OP_DZ),weight_79,79) &
          +  int5(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DZZ),h(:,OP_DZ),weight_79,79) &
          +  int5(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_LP ),h(:,OP_DR),weight_79,79) &
          -  int5(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_LP ),h(:,OP_DZ),weight_79,79) &
          +  int5(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_LP ),h(:,OP_1 ),weight_79,79) &
          -  int5(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_LP ),h(:,OP_1 ),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          +  int4(e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_LP),h(:,OP_1),weight_79,79) &
          +  int4(e(:,OP_1 ),f(:,OP_DZZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          +  int4(e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          +  int4(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          +  int4(e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     endif
  endif

  v1uun = temp
  return
end function v1uun


! v1uchin
! ~~~~~~~
real function v1uchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  real :: temp

  if(idens.eq.0) then
     temp =-int3(e(:,OP_DZ),f(:,OP_LP),g(:,OP_DZ),weight_79,79) &
          - int3(e(:,OP_DR),f(:,OP_LP),g(:,OP_DR),weight_79,79)
  else 
     temp =-int4(e(:,OP_DZ),f(:,OP_LP ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
          - int4(e(:,OP_DR),f(:,OP_LP ),g(:,OP_DR ),h(:,OP_1 ),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_LP ),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_LP ),g(:,OP_DR ),h(:,OP_DR),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_DR),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_DZ),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRR),h(:,OP_DZ),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_DR),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_DZ),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DZZ),h(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
             + int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_DR),weight_79,79)
     endif
  endif

  v1uchin = temp
end function v1uchin


! v1chichin
! =========
real function v1chichin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  real :: temp

  if(idens.eq.0) then
     temp = 0
  else
     temp = -0.5* &
          (int5(ri_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DR),weight_79,79) &
          +int5(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DR),h(:,OP_DR),weight_79,79) &
          -int5(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          -int5(ri_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),h(:,OP_DZ),weight_79,79))
  endif
  
  v1chichin = temp

  return
end function v1chichin



! V1psisb1
! ========
real function v1psisb1(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),weight_79,79) &
       - int4(ri3_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
       + int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS),weight_79,79) &
       - int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS),weight_79,79)

  v1psisb1 = temp
  return
end function v1psisb1


! V1psipsi
! ========
real function v1psipsi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp
          
  temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),weight_79,79) &
       - int4(ri3_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),weight_79,79)

  v1psipsi = temp
  return
end function v1psipsi


! V1upsipsi
! =========
real function v1upsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  ! |psi(1), u|
  temp79a = f(:,OP_DR )*g(:,OP_DZ ) - f(:,OP_DZ )*g(:,OP_DR )
  ! |psi(1), u|,r
  temp79b = f(:,OP_DRR)*g(:,OP_DZ ) - f(:,OP_DRZ)*g(:,OP_DR ) &
       +    f(:,OP_DR )*g(:,OP_DRZ) - f(:,OP_DZ )*g(:,OP_DRR)
  ! |psi(1), u|,z
  temp79c = f(:,OP_DRZ)*g(:,OP_DZ ) - f(:,OP_DZZ)*g(:,OP_DR ) &
       +    f(:,OP_DR )*g(:,OP_DZZ) - f(:,OP_DZ )*g(:,OP_DRZ)

  ! |psi(2), nu|
  temp79d = e(:,OP_DR )*h(:,OP_DZ ) - e(:,OP_DZ )*h(:,OP_DR)
  ! |psi(2), nu|,r
  temp79e = e(:,OP_DRR)*h(:,OP_DZ ) - e(:,OP_DRZ)*h(:,OP_DR ) &
       +    e(:,OP_DR )*h(:,OP_DRZ) - e(:,OP_DZ )*h(:,OP_DRR)
  ! |psi(2), nu|,z
  temp79f = e(:,OP_DRZ)*h(:,OP_DZ ) - e(:,OP_DZZ)*h(:,OP_DR ) &
       +    e(:,OP_DR )*h(:,OP_DZZ) - e(:,OP_DZ )*h(:,OP_DRZ)
  
  temp = int4(ri2_79,e(:,OP_DR),temp79c,h(:,OP_GS),weight_79,79) &
       - int4(ri2_79,e(:,OP_DZ),temp79b,h(:,OP_GS),weight_79,79) &
       - int3(ri2_79,temp79b,temp79e,weight_79,79) &
       - int3(ri2_79,temp79c,temp79f,weight_79,79)
  if(itor.eq.1) then
     temp = temp &
       + int3(ri3_79,temp79b,temp79d,weight_79,79) &
       - int3(ri3_79,temp79a,temp79e,weight_79,79) &
       - int4(ri3_79,e(:,OP_DZ),temp79a,h(:,OP_GS),weight_79,79) &
       + int3(ri4_79,temp79a,temp79d,weight_79,79)
  endif

  v1upsipsi = temp
  return
end function v1upsipsi


! V1vpsib
! =======
real function v1vpsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     temp = int5(ri6_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          - int5(ri6_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          + int5(ri7_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1),weight_79,79)
  endif

  v1vpsib = temp
  return
end function v1vpsib


! V1bb 
! ====
real function v1bb(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     temp = -2.*int4(ri4_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79) 
  endif

  v1bb = temp
  return
end function v1bb


! V1ubb 
! =====
real function v1ubb(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     temp = 2.* &
          (int5(ri5_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),g(:,OP_1),weight_79,79) &
          -int5(ri5_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),g(:,OP_1),weight_79,79))
  endif
  
  v1ubb = temp
  return
end function v1ubb


! V1chipsipsi
! ===========
real function v1chipsipsi(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
  temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
  temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)

  temp = int4(ri2_79,h(:,OP_GS ),e(:,OP_DR ),temp79c,weight_79,79) &
       - int4(ri2_79,h(:,OP_GS ),e(:,OP_DZ ),temp79b,weight_79,79) &
       + int4(ri3_79,e(:,OP_DZZ),h(:,OP_DR ),temp79c,weight_79,79) &
       - int4(ri3_79,e(:,OP_DRZ),h(:,OP_DZ ),temp79c,weight_79,79) &
       + int4(ri3_79,e(:,OP_DZ ),h(:,OP_DRZ),temp79c,weight_79,79) &
       - int4(ri3_79,e(:,OP_DR ),h(:,OP_DZZ),temp79c,weight_79,79) &
       + int4(ri3_79,e(:,OP_DRZ),h(:,OP_DR ),temp79b,weight_79,79) &
       - int4(ri3_79,e(:,OP_DRR),h(:,OP_DZ ),temp79b,weight_79,79) &
       + int4(ri3_79,e(:,OP_DZ ),h(:,OP_DRR),temp79b,weight_79,79) &
       - int4(ri3_79,e(:,OP_DR ),h(:,OP_DRZ),temp79b,weight_79,79)

  if(itor.eq.1) then
     temp = temp &
          + int4(ri4_79,e(:,OP_DR),h(:,OP_DZ),temp79b,weight_79,79) &
          - int4(ri4_79,e(:,OP_DZ),h(:,OP_DR),temp79b,weight_79,79) &
          + int4(ri4_79,e(:,OP_DR),h(:,OP_DZ),temp79b,weight_79,79) &
          - int4(ri4_79,e(:,OP_DZ),h(:,OP_DR),temp79b,weight_79,79)
  endif

  v1chipsipsi = temp
  return
end function v1chipsipsi



! V1chibb
! =======
real function v1chibb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(itor.eq.0) then
     v1chibb = 0
     return
  endif

  temp = -2.* &
       (int5(ri6_79,e(:,OP_DZ),g(:,OP_1 ),f(:,OP_GS),h(:,OP_1),weight_79,79) &
       +int5(ri6_79,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_DZ),h(:,OP_1),weight_79,79) &
       +int5(ri6_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DR),h(:,OP_1),weight_79,79))
  
  v1chibb = temp
  return
end function v1chibb



! V1bsb2
! ======
real function v1bsb2(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  if(itor.eq.0) then 
     temp = 0.
  else
     temp = 2.*int4(ri4_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_1),weight_79,79)
  endif

  v1bsb2 = temp
  return
end function v1bsb2


! V1ngrav
! =======
real function v1ngrav(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int3(ri_79,e(:,OP_1),f(:,OP_DR),weight_79,79)

  v1ngrav = temp
  return
end function v1ngrav


! V1ungrav
! ========
real function v1ungrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int3(e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
       temp = temp - 2.*int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  endif

  v1ungrav = temp
  return
end function v1ungrav


! V1chingrav
! ==========
real function v1chingrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       + int4(ri_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
       + int4(ri_79,e(:,OP_DR),f(:,OP_LP),g(:,OP_1 ),weight_79,79)

  v1chingrav = temp
  return
end function v1chingrav


! V1ndenmgrav
! ===========
real function v1ndenmgrav(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = -int3(ri_79,e(:,OP_DR),f(:,OP_LP),weight_79,79)

  v1ndenmgrav = temp
  return
end function v1ndenmgrav


!===============================================================================
! V2 TERMS
!===============================================================================


! V2vn
! ====
real function v2vn(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  if(idens.eq.0) then
     temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)
  else
     temp = int3(e(:,OP_1),f(:,OP_1),g(:,OP_1),weight_79,79)
  endif

  v2vn = temp
  return
end function v2vn


! V2vmu
! =====
real function v2vmu(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = -int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       -  int2(e(:,OP_DR),f(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp - int3(ri3_79,e(:,OP_1),f(:,OP_1),weight_79,79)
  endif

  v2vmu = temp
  return
end function v2vmu


! V2vun
! =====
real function v2vun(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(idens.eq.0) then
     temp = int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  else
     temp = int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
  endif

  v2vun = temp
  return
end function v2vun


! V2psib
! ======
real function v2psib(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  v2psib = temp
  return
end function v2psib


! V2vpsipsi
! =========
real function v2vpsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp


  ! [nu,psi(2)] + nu*psi(2)_z/r
  temp79a = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
  if(itor.eq.1) temp79a = temp79a + ri_79*e(:,OP_1)*h(:,OP_DZ)

  ! [v,psi(1)] + v*psi(1)_z/r
  temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
  if(itor.eq.1) temp79b = temp79b + ri_79*f(:,OP_1)*g(:,OP_DZ)

  temp = -int3(ri4_79, temp79a, temp79b, weight_79, 79)

  v2vpsipsi = temp
  return
end function v2vpsipsi


! V2upsib
! =========
real function v2upsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp


  ! [nu,bz] + nu*bz_z/r
  temp79a = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
  if(itor.eq.1) temp79a = temp79a + ri_79*e(:,OP_1)*h(:,OP_DZ)

  ! [nu,psi] + nu*psi_z/r
  temp79b = e(:,OP_DZ)*g(:,OP_DR) - e(:,OP_DR)*g(:,OP_DZ)
  if(itor.eq.1) temp79b = temp79b + ri_79*e(:,OP_1)*g(:,OP_DZ)

  temp = int4(ri_79, f(:,OP_DZ), g(:,OP_DR), temp79a, weight_79, 79)  &
       - int4(ri_79, f(:,OP_DR), g(:,OP_DZ), temp79a, weight_79, 79)  &
       + int4(ri3_79, f(:,OP_DR), h(:,OP_DZ), temp79b, weight_79, 79) &
       - int4(ri3_79, f(:,OP_DZ), h(:,OP_DR), temp79b, weight_79, 79)

  v2upsib = temp
  return
end function v2upsib


! v2upsisb2
! =========
real function v2upsisb2(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  temp = int5(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
       - int5(ri_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
       - int5(ri_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DR),weight_79,79) &
       + int5(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp &
          + int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          - int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DZ),weight_79,79)
  endif

  v2upsisb2 = temp
  return
end function v2upsisb2


! v2chipsib
! =========
real function v2chipsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp


  temp79a = h(:,OP_1 )*f(:,OP_GS) &
       +    h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
  temp79b = e(:,OP_DZ)*g(:,OP_DR) - e(:,OP_DR)*g(:,OP_DZ)

  temp = int5(ri2_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
       + int5(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
       - int5(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DR),weight_79,79) &
       - int5(ri2_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DR),h(:,OP_DR),weight_79,79) &
       + int3(ri4_79,temp79a,temp79b,weight_79,79)
  if(itor.eq.1) then
     temp = temp &
          - int5(ri3_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          - int5(ri3_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          + int4(ri5_79,temp79a,e(:,OP_1),g(:,OP_DZ),weight_79,79)
  endif

  v2chipsib = temp
  return
end function v2chipsib


! v2vchin
! =======
real function v2vchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(idens.eq.0) then
     temp =-int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)
  else
     temp =-int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)
  endif

  v2vchin = temp
  return
end function v2vchin



!===============================================================================
! V3 TERMS
!===============================================================================

! V3chin
! ======
real function v3chin(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  
  real :: temp

  if(idens.eq.0) then
     temp = - int2(e(:,OP_DR),f(:,OP_DR),weight_79,79) &
          -   int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79)
  else
     temp = - int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
          -   int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  endif

  v3chin = temp
  return
end function v3chin



! V3chimu
! =======
real function v3chimu(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_LP),f(:,OP_LP),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*int3(ri_79,e(:,OP_DR),f(:,OP_LP),weight_79,79)
  endif

  v3chimu = temp
  return
end function v3chimu


! V3umu
! =====
real function v3umu(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     temp = -4.*int3(ri_79,e(:,OP_DZ),f(:,OP_DR),weight_79,79)
  endif

  v3umu = temp
  return
end function v3umu


! V3un
! ====
real function v3un(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  if(idens.eq.0) then
     temp = 0.
  else
     temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) 
  endif

  v3un = temp
  return
end function v3un


! V3p
! ===
real function v3p(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       + int2(e(:,OP_DR),f(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.*int3(ri_79,e(:,OP_1),f(:,OP_DR),weight_79,79)
  endif

  v3p = temp
  return
end function v3p



! V3up
! ====
real function v3up(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(r_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(r_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          (int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          -int3(e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          -2.*gam* &
          (int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          +int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79)))
  endif

  v3up = temp
  
  return
end function v3up


! V3chip
! ======
real function v3chip(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       + int3(e(:,OP_LP),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
       + gam*int3(e(:,OP_LP),f(:,OP_LP),g(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.* &
          (int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          +int4(ri_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          +gam*int4(ri_79,e(:,OP_DR),f(:,OP_LP),g(:,OP_1),weight_79,79))
  endif

  v3chip = temp

  return
end function v3chip



! V3psipsi
! ========
real function v3psipsi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.*int4(ri3_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DR),weight_79,79)
  endif

  v3psipsi = temp

  return
end function v3psipsi


! V3bb
! ====
real function v3bb(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          int4(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),weight_79,79)
  endif

  v3bb = temp
  return
end function v3bb


! V3psisb1
! ========
real function v3psisb1(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          int4(ri3_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DR),weight_79,79)
  endif

  v3psisb1 = temp

end function v3psisb1


! V3bsb2
! ======
real function v3bsb2(e,f,g)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = -int4(ri2_79,e(:,OP_LP),f(:,OP_1),g(:,OP_1),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 4.*int4(ri4_79,e(:,OP_1),f(:,OP_1),g(:,OP_1),weight_79,79)
  endif

  v3bsb2 = temp
  return
end function v3bsb2


! V3upsipsi
! =========
real function v3upsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  ! Linearize in u
  ! ~~~~~~~~~~~~~~

  ! [b,c]
  temp79a = f(:,OP_DZ )*g(:,OP_DR ) - f(:,OP_DR )*g(:,OP_DZ)

  ! [b,c],r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)

  ! [b,c],z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)

  temp =-int4(ri_79,temp79b,e(:,OP_DRZ),h(:,OP_DRZ),weight_79,79) & 
       - int4(ri_79,temp79b,e(:,OP_DRR),h(:,OP_DRR),weight_79,79) &
       - int4(ri_79,temp79c,e(:,OP_DZZ),h(:,OP_DZZ),weight_79,79) & 
       - int4(ri_79,temp79c,e(:,OP_DRZ),h(:,OP_DRZ),weight_79,79) &
       + int4(ri_79,e(:,OP_DR),temp79b,h(:,OP_GS),weight_79,79) &
       + int4(ri_79,e(:,OP_DZ),temp79c,h(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     temp = temp &
          + int4(ri2_79,e(:,OP_DR ),temp79a,h(:,OP_GS ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DRZ),temp79a,h(:,OP_DZ ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DRR),temp79a,h(:,OP_DR ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DZ ),temp79a,h(:,OP_DRZ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DR ),temp79a,h(:,OP_DRR),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_1  ),temp79b,h(:,OP_DRR),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_1  ),temp79c,h(:,OP_DRZ),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_DR ),temp79b,h(:,OP_DR ),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_DZ ),temp79c,h(:,OP_DR ),weight_79,79) &
          + 2.*int4(ri2_79,e(:,OP_1  ),temp79b,h(:,OP_GS ),weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_1  ),temp79b,h(:,OP_DR ),weight_79,79) &
          - 2.*int4(ri3_79,e(:,OP_DR ),temp79a,h(:,OP_DR ),weight_79,79) &
          - 2.*int4(ri3_79,e(:,OP_1  ),temp79a,h(:,OP_DRR),weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_1  ),temp79a,h(:,OP_GS ),weight_79,79) &
          + 2.*int4(ri4_79,e(:,OP_1  ),temp79a,h(:,OP_DR ),weight_79,79)
  endif

  v3upsipsi = temp
  return
end function v3upsipsi


! V3ubb
! =====
real function v3ubb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  real :: temp

  temp = int5(ri3_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       - int5(ri3_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)

  if(itor.eq.1) then
     temp = temp - 4.* &
          (int5(ri5_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          -int5(ri5_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79))
  endif

  v3ubb = temp
  return
end function v3ubb


! v3vpsib
! =======
real function v3vpsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  temp = int5(ri4_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       - int5(ri4_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp &
          - int5(ri5_79,e(:,OP_LP),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          -4.* &
          (int5(ri6_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          -int5(ri6_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          -int5(ri7_79,e(:,OP_1),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1),weight_79,79))
  endif

  v3vpsib = temp
  return
end function v3vpsib


! V3chipsipsi
! ===========
real function v3chipsipsi(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  ! <f,g>
  temp79a = f(:,OP_DR )*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DZ )

  ! <f,g>,r
  temp79b = f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DRZ) &
       +    f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DRZ)

  ! <f,g>,z
  temp79c = f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DZZ) &
       +    f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DZZ)

  temp = int4(ri2_79,temp79b,e(:,OP_DRZ),h(:,OP_DZ ),weight_79,79) &
       + int4(ri2_79,temp79b,e(:,OP_DRR),h(:,OP_DR ),weight_79,79) &
       + int4(ri2_79,temp79b,e(:,OP_DZ ),h(:,OP_DRZ),weight_79,79) &
       + int4(ri2_79,temp79b,e(:,OP_DR ),h(:,OP_DRR),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DZZ),h(:,OP_DZ ),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DRZ),h(:,OP_DR ),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DZ ),h(:,OP_DZZ),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DR ),h(:,OP_DRZ),weight_79,79) &
       - int4(ri2_79,temp79c,e(:,OP_DZ ),h(:,OP_GS ),weight_79,79) &
       - int4(ri2_79,temp79b,e(:,OP_DR ),h(:,OP_GS ),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.* &
          (int4(ri3_79,temp79c,e(:,OP_1  ),h(:,OP_DRZ),weight_79,79) &
          +int4(ri3_79,temp79b,e(:,OP_1  ),h(:,OP_DRR),weight_79,79) &
          +int4(ri3_79,temp79c,e(:,OP_DRZ),h(:,OP_1  ),weight_79,79) &
          +int4(ri3_79,temp79b,e(:,OP_DRR),h(:,OP_1  ),weight_79,79) &
          -int4(ri3_79,temp79b,e(:,OP_1  ),h(:,OP_GS ),weight_79,79) &
          -int4(ri4_79,temp79b,e(:,OP_1  ),h(:,OP_DR ),weight_79,79))
  endif

  v3chipsipsi = temp

  return
end function v3chipsipsi


! V3chibb
! =======
real function v3chibb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  temp = int5(ri4_79,e(:,OP_LP),f(:,OP_GS),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
       + int5(ri4_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       + int5(ri4_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp - 4.*&
          (int5(ri6_79,e(:,OP_1),f(:,OP_GS),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
          +int5(ri6_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          +int5(ri6_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79))
  endif

  v3chibb = temp

  return
end function v3chibb


! V3uun
! =====
real function v3uun(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(idens.eq.0) then
     temp = int4(r2_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DRZ),weight_79,79) &
          - int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DRR),weight_79,79) &
          + int4(r2_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DRZ),weight_79,79) &
          - int4(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZZ),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ ),weight_79,79) &
             - int4(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ ),weight_79,79) &
             - int4(r_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DZ ),weight_79,79) &
             + int4(r_79,e(:,OP_DR),g(:,OP_DZ),f(:,OP_DZ ),weight_79,79) &
             + int4(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DRZ),weight_79,79) &
             - int4(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DZZ),weight_79,79)
     endif
  else
     temp = int5(r2_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          - int5(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DRR),h(:,OP_1),weight_79,79) &
          + int5(r2_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          - int5(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZZ),h(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + int5(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             + int5(r_79,e(:,OP_DR),g(:,OP_DZ),f(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             + int5(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DRZ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DZZ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
             + int5(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79)
     endif
  endif

  v3uun = temp
  return
end function v3uun


! V3uchin
! =======
real function v3uchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(idens.eq.0) then
     temp = 2.* &
          (int3(e(:,OP_1),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
          -int3(e(:,OP_1),f(:,OP_DZZ),g(:,OP_DZ ),weight_79,79) &
          -int3(e(:,OP_1),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
          -int3(e(:,OP_1),f(:,OP_DZ ),g(:,OP_DRR),weight_79,79))
     if(itor.eq.1) then
        temp = temp &
             -2.*int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     endif
  else
     temp = 2.* &
          (int4(e(:,OP_1),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          -int4(e(:,OP_1),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          -int4(e(:,OP_1),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          -int4(e(:,OP_1),f(:,OP_DZ ),g(:,OP_DRR),h(:,OP_1),weight_79,79))
     if(itor.eq.1) then
        temp = temp &
             -2.*int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
     endif
  endif

  v3uchin = temp
  return
end function v3uchin



! V3chichin
! =========
real function v3chichin(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(idens.eq.0) then
     temp = int3(e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          - int4(ri_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DZ ),weight_79,79) &
          - int4(ri_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DR ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DRR),g(:,OP_DR ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DRR),weight_79,79)
     endif
  else
     temp = int4(e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),h(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DR ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DRR),h(:,OP_1 ),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DR ),h(:,OP_DR),weight_79,79)
     endif
  endif

  temp = temp / 2.
  v3chichin = temp

  return
end function v3chichin


! V3ngrav
! =======
real function v3ngrav(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = -int2(e(:,OP_1),f(:,OP_DZ),weight_79,79)

  v3ngrav = temp
  return
end function v3ngrav


! V3ungrav
! ========
real function v3ungrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(r_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
       - int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.*int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  endif

  v3ungrav = temp
  return
end function v3ungrav


! V3chingrav
! ==========
real function v3chingrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp =-int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       - int3(e(:,OP_DZ),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
       - int3(e(:,OP_DZ),f(:,OP_LP),g(:,OP_1 ),weight_79,79)

  v3chingrav = temp
  return
end function v3chingrav


! V3ndenmgrav
! ===========
real function v3ndenmgrav(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_DZ),f(:,OP_LP),weight_79,79)

  v3ndenmgrav = temp
  return
end function v3ndenmgrav


!===============================================================================
! B1 TERMS
!===============================================================================


! B1psi
! =====
real function b1psi(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)

  b1psi = temp
  return
end function b1psi


! B1psiu
! ======
real function b1psiu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  
  b1psiu = temp
  return
end function b1psiu


! B1psichi
! =======
real function b1psichi(e,f,g)

  use basic
  use nintegrate_mod
  
  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = -int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       -  int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)

  b1psichi = temp
  return
end function


! B1psieta
! ========
real function b1psieta(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp =-int2(e(:,OP_DR),f(:,OP_DR),weight_79,79) &
       - int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79)
  if(itor.eq.1) then
     temp = temp - 2.*int3(ri_79,e(:,OP_1),f(:,OP_DR),weight_79,79)
  endif

  b1psieta = temp
  return
end function b1psieta


! B1psibd
! =======
real function b1psibd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  ! Linearize in psi
  ! ~~~~~~~~~~~~~~~~
  if(idens.eq.0) then
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79)
  else
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79)
  endif

  b1psibd = temp
  return
end function b1psibd


!===============================================================================
! B2 TERMS
!===============================================================================

! B2b
! ===
real function b2b(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)

  b2b = temp
  return
end function b2b


! B2beta
! ======
real function b2beta(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = -int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       -  int3(ri2_79,e(:,OP_DR),f(:,OP_DR),weight_79,79)

  b2beta = temp
  return
end function b2beta


! B2bu
! ====
real function b2bu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  b2bu = temp
  return
end function b2bu


! B2bchi
! ======
real function b2bchi(e,f,g)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)

  b2bchi = temp
  return
end function b2bchi


! B2psiv
! ======
real function b2psiv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp + int4(ri3_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  endif

  b2psiv = temp
  return
end function b2psiv


! B2psipsid
! =========
real function b2psipsid(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(idens.eq.0) then
     temp = int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS),weight_79,79) &
          - int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS),weight_79,79)
  else
     temp = int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS),ni79(:,OP_1),weight_79,79) &
          - int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS),ni79(:,OP_1),weight_79,79)
  endif
  
  b2psipsid = temp
  return 
end function b2psipsid


! B2bbd
! =====
real function b2bbd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(idens.eq.0) then
     temp = 0.
     if(itor.eq.1) then
        temp = temp &
             + 2.*int4(ri4_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
     endif
  else
     temp = int5(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          - int5(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + 2.*int5(ri4_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     endif
  endif
  
  b2bbd = temp
  return 
end function b2bbd



! B2ped
! =====
real function b2ped(e,f,g)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp
  
  if(idens.eq.0) then
     temp = 0.
  else
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  endif

  b2ped = temp
  return
end function b2ped

!===============================================================================
! B3 TERMS
!===============================================================================

! B3pe
! ====
real function b3pe(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)

  b3pe = temp
  return
end function b3pe


! B3psipsieta
! ~~~~~~~~~~~
real function b3psipsieta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  
  real :: temp

  if(itor.eq.0) then
     b3psipsieta = 0.
     return
  endif

  temp = (gam-1.)*int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),weight_79,79)

  b3psipsieta = temp
  
  return
end function b3psipsieta


! B3bbeta
! ~~~~~~~
real function b3bbeta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  
  real :: temp

  if(itor.eq.0) then
     b3bbeta = 0.
     return
  endif

  temp = (gam-1.)* &
       (int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       +int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79))

  b3bbeta = temp
  
  return
end function b3bbeta


! B3pebd
! ~~~~~~
real function b3pebd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  
  real :: temp

  if(idens.eq.0) then
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79)
  else 
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          + gam* &
          (int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          -int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR),weight_79,79))
  endif

  b3pebd = temp
  
  return
end function b3pebd



!===============================================================================
! N1 TERMS
!===============================================================================

! N1n
! ===
real function n1n(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)

  n1n = temp
  return
end function n1n


! N1ndenm
! =======
real function n1ndenm(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp =-int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       - int2(e(:,OP_DR),f(:,OP_DR),weight_79,79)

  n1ndenm = temp
  return
end function n1ndenm



! N1nu
! ====
real function n1nu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(r_79, e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(r_79, e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*int3(e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
  endif

  n1nu = temp
  return
end function n1nu


! N1nchi
! ======
real function n1nchi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
       + int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)

  n1nchi = temp
  return
end function n1nchi


!===============================================================================
! P1 TERMS
!===============================================================================

! P1pu
! ====
real function p1pu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*gam*int4(r_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
  endif

  p1pu = temp

  return
end function p1pu


! P1pchi
! ======
real function p1pchi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = gam* &
       (int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79)  &
       +int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)) &
       +(gam-1.)* &
       (int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79)  &
       +int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79))

  p1pchi = temp

  return
end function p1pchi

end module metricterms_n

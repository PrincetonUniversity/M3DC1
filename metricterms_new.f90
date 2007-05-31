module metricterms_new

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

  temp = -int2(e(:,OP_GS),f(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     temp = temp &
          - 4.*int3(ri_79,e(:,OP_DR),f(:,OP_GS),weight_79,79)
  endif

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
     if(itor.eq.1) then
        temp = temp + 2.*int3(ri_79,e(:,OP_1),f(:,OP_DR),weight_79,79)
     endif
  else
     temp = int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),weight_79,79)
     endif
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
     temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  endif

  v1chin = temp
  return
end function v1chin


! V1psipsi
! ========
real function v1psipsi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp
          
  temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS),weight_79,79) &
       - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     temp = temp - 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_GS),weight_79,79)
  endif

  v1psipsi = temp
  return
end function v1psipsi


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
     temp = -2.*int4(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79) 
  endif

  v1bb = temp
  return
end function v1bb


! V1uun
! =====
real function v1uun(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp
           

  if(idens.eq.0) then
     temp =  int4(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
          -  int4(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.*int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZ),weight_79,79) 
     endif
  else
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

     temp = int5(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          + 0.5*(int4(ri_79,temp79a,e(:,OP_DR),h(:,OP_DZ),weight_79,79) &
                -int4(ri_79,temp79a,e(:,OP_DZ),h(:,OP_DR),weight_79,79))

     if(itor.eq.1) then
        temp = temp + &
             2.*int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
             +  int4(ri2_79,temp79a,e(:,OP_1),h(:,OP_DZ),weight_79,79)
     endif

  endif

  v1uun = temp
  return
end function v1uun


! v1vvn
! =====
real function v1vvn(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp  

  if(itor.eq.0) then
     temp = 0.
  else
     if(idens.eq.0) then
        temp = -int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_1),weight_79,79)
     else
        temp = -int5(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_1),h(:,OP_1),weight_79,79)
     endif
  endif

  v1vvn = temp
  return
end function v1vvn


! v1uchin
! =======
real function v1uchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  real :: temp

  if(idens.eq.0) then
     temp = -int3(e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
          -  int3(e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DR),weight_79,79) 
     endif
  else
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

     temp =-int4(e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int4(e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          + int3(temp79a,e(:,OP_DZ),h(:,OP_DR),weight_79,79) &
          - int3(temp79a,e(:,OP_DR),h(:,OP_DZ),weight_79,79)

     if(itor.eq.1) then
        temp = temp - 2.*&
             (int5(ri_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79) &
             +int4(ri_79,temp79a,e(:,OP_1),h(:,OP_DZ),weight_79,79))
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
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

     temp = 0.5* &
          (int4(r_79,e(:,OP_DR),temp79a,h(:,OP_DZ),weight_79,79) &
          -int4(r_79,e(:,OP_DZ),temp79a,h(:,OP_DR),weight_79,79))

     if(itor.eq.1) then
        temp = temp + int3(e(:,OP_1),temp79a,h(:,OP_DZ),weight_79,79)
     endif
         
  endif
  
  v1chichin = temp

  return
end function v1chichin



! V1upsipsi
! =========
real function v1upsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  ! |u, psi(1)|,r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
  ! |u, psi(1)|,z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)

  ! |nu, psi(2)|,r
  temp79e = e(:,OP_DRZ)*h(:,OP_DR ) - e(:,OP_DRR)*h(:,OP_DZ ) &
       +    e(:,OP_DZ )*h(:,OP_DRR) - e(:,OP_DR )*h(:,OP_DRZ)
  ! |nu, psi(2)|,z
  temp79f = e(:,OP_DZZ)*h(:,OP_DR ) - e(:,OP_DRZ)*h(:,OP_DZ ) &
       +    e(:,OP_DZ )*h(:,OP_DRZ) - e(:,OP_DR )*h(:,OP_DZZ)
  
  temp =-int3(ri2_79,temp79b,temp79e,weight_79,79) &
       - int3(ri2_79,temp79c,temp79f,weight_79,79) &
       + int4(ri2_79,e(:,OP_DZ),temp79b,h(:,OP_GS),weight_79,79) &
       - int4(ri2_79,e(:,OP_DR),temp79c,h(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     ! |u, psi(1)|
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

     ! |nu,psi(2)|
     temp79d = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)

     temp = temp &
          -    int4(ri3_79,e(:,OP_DZ),temp79a,h(:,OP_GS),weight_79,79) &
          - 2.*int4(ri3_79,e(:,OP_1 ),temp79c,h(:,OP_GS),weight_79,79) &
          +    int3(ri3_79,temp79e,temp79a,weight_79,79) &
          -    int3(ri3_79,temp79d,temp79b,weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_1 ),temp79b,h(:,OP_DRZ),weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_1 ),temp79c,h(:,OP_DZZ),weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_DR),temp79b,h(:,OP_DZ ),weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_DZ),temp79c,h(:,OP_DZ ),weight_79,79) &
          +    int3(ri4_79,temp79a,temp79d,weight_79,79) &
          - 2.*int4(ri4_79,e(:,OP_DR),temp79a,h(:,OP_DZ ),weight_79,79) &
          - 2.*int4(ri4_79,e(:,OP_1 ),temp79a,h(:,OP_DRZ),weight_79,79)
  endif

  v1upsipsi = temp
  return
end function v1upsipsi

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
          (   int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          -   int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          -2.*int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1),weight_79,79))
  endif
  
  v1ubb = temp
  return
end function v1ubb


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
     temp = 2.*(   int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
               -   int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
               +2.*int5(ri4_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1),weight_79,79))
  endif

  v1vpsib = temp
  return
end function v1vpsib



! V1chipsipsi
! ===========
real function v1chipsipsi(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
  temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)

  temp79d = temp79b * &
        (e(:,OP_DRZ)*h(:,OP_DR ) - e(:,OP_DRR)*h(:,OP_DZ )  &
        +e(:,OP_DZ )*h(:,OP_DRR) - e(:,OP_DR )*h(:,OP_DRZ)) & 
           +temp79c * &
        (e(:,OP_DZZ)*h(:,OP_DR ) - e(:,OP_DRZ)*h(:,OP_DZ )  &
        +e(:,OP_DZ )*h(:,OP_DRZ) - e(:,OP_DR )*h(:,OP_DZZ))


  temp = int2(ri_79,temp79d,weight_79,79) &
       + int4(ri_79,e(:,OP_DR),temp79c,h(:,OP_GS),weight_79,79) &
       - int4(ri_79,e(:,OP_DZ),temp79b,h(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

     temp = temp + 2.* &
          (int4(ri2_79,e(:,OP_1 ),temp79c,h(:,OP_GS ),weight_79,79) &
          -int4(ri2_79,e(:,OP_1 ),temp79b,h(:,OP_DRZ),weight_79,79) &
          -int4(ri2_79,e(:,OP_1 ),temp79c,h(:,OP_DZZ),weight_79,79) &
          -int4(ri2_79,e(:,OP_DR),temp79b,h(:,OP_DZ ),weight_79,79) &
          -int4(ri2_79,e(:,OP_DZ),temp79c,h(:,OP_DZ ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DZ),temp79b,h(:,OP_DR ),weight_79,79) &
          -int4(ri2_79,e(:,OP_DR),temp79b,h(:,OP_DZ ),weight_79,79))
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
     v1chibb = 0.
     return
  endif

  temp = -2.* &
       (int5(ri2_79,e(:,OP_DZ),g(:,OP_1 ),f(:,OP_GS),h(:,OP_1),weight_79,79) &
       +int5(ri2_79,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_DZ),h(:,OP_1),weight_79,79) &
       +int5(ri2_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DR),h(:,OP_1),weight_79,79))
  
  v1chibb = temp
  return
end function v1chibb


! V1psisb1
! ========
real function v1psisb1(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),weight_79,79) &
       - int4(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
       + int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS),weight_79,79) &
       - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     temp = temp - 2.* &
          (int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
          +int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_GS),weight_79,79))
  endif

  v1psisb1 = temp
  return
end function v1psisb1



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
     temp = 2.*int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_1),weight_79,79)
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

  temp = int3(r_79,e(:,OP_1),f(:,OP_DR),weight_79,79)

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

  temp79a = g(:,OP_DZ)*f(:,OP_DR) - f(:,OP_DZ)*g(:,OP_DR)

  temp = int2(temp79a,e(:,OP_DR),weight_79,79)
  if(itor.eq.1) then 
     temp = temp + 2.*int3(ri_79,temp79a,e(:,OP_1),weight_79,79)
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

  temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) + f(:,OP_LP)*g(:,OP_1)

  temp = int3(r_79,temp79a,e(:,OP_DR),weight_79,79)
  if(itor.eq.1) then 
     temp = temp + 2.*int2(temp79a,e(:,OP_1),weight_79,79)
  endif

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

  temp = -int3(r_79,e(:,OP_DR),f(:,OP_LP),weight_79,79)

  if(itor.eq.1) then
     temp = temp - 2.*int2(e(:,OP_1),f(:,OP_LP),weight_79,79)
  endif

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
real function v2vmu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real, intent(in) :: g
  real :: temp

  temp = -g* &
       (int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       +int2(e(:,OP_DR),f(:,OP_DR),weight_79,79))

  if(itor.eq.1) then
     temp = temp - 2.*g*int3(ri_79,e(:,OP_1 ),f(:,OP_DR),weight_79,79)
  endif

  v2vmu = temp
  return
end function v2vmu


! V2vhypv
! =======
real function v2vhypv(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real, intent(in) :: g,h
  real :: temp

  temp = -g*h*int2(e(:,OP_GS),f(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     temp = temp - 4.*g*h* &
          int3(ri_79,e(:,OP_DR),f(:,OP_GS),weight_79,79)
  endif

  v2vhypv = temp
  return
end function v2vhypv


! V2vun
! =====
real function v2vun(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(idens.eq.0) then
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  else
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
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

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

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


  ! [nu,psi(2)]
  temp79a = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)

  temp = int4(ri2_79,f(:,OP_DR),g(:,OP_DZ),temp79a,weight_79, 79) &
       - int4(ri2_79,f(:,OP_DZ),g(:,OP_DR),temp79a,weight_79, 79)
  if(itor.eq.1) then
     temp = temp - 2.*int4(ri3_79,f(:,OP_1),g(:,OP_DZ),temp79a,weight_79,79)
  endif

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

  temp79a = e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ)

  temp = int4(ri2_79,temp79a,g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
       - int4(ri2_79,temp79a,g(:,OP_DZ),h(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.* &
          (int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          -int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79))
  endif

  v2upsib = temp
  return
end function v2upsib


! v2upsisb2
! ========
real function v2upsisb2(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  v2upsisb2 = 0.
  return
end function v2upsisb2


! v2ubsb1
! =======
real function v2ubsb1(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  v2ubsb1 = 0.
  return
end function v2ubsb1


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

  temp79b = e(:,OP_DR)*h(:,OP_DZ) - e(:,OP_DZ)*h(:,OP_DR)
  
  temp = int4(ri_79,temp79a,e(:,OP_DZ),g(:,OP_DR),weight_79,79) &
       - int4(ri_79,temp79a,e(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       + int4(ri_79,temp79b,f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       + int4(ri_79,temp79b,f(:,OP_DR),g(:,OP_DR),weight_79,79) 

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
     temp =-int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          - int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)
  else
     temp =-int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)
  endif

  v2vchin = temp
  return
end function v2vchin


! v2chibsb1
! =========
real function v2chibsb1(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  v2chibsb1 = 0.
  return
end function v2chibsb1

! v2psisb2
! ========
real function v2psisb2(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  v2psisb2 = temp
  return
end function v2psisb2

! v2bsb1
! ======
real function v2bsb1(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79)

  v2bsb1 = temp
  return
end function v2bsb1




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

  v3chimu = temp
  return
end function v3chimu


! V3umu
! =====
real function v3umu(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real, intent(in) :: g,h
  real :: temp

  v3umu = 0.
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
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) 
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

  temp = int4(ri_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

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

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS),weight_79,79)

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

  temp = int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       + int3(ri2_79,e(:,OP_DR),f(:,OP_DR),weight_79,79)

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

  ! [f,g],r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)

  ! [f,g],z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)

  temp = int4(ri3_79,e(:,OP_DZ ),temp79c,h(:,OP_GS ),weight_79,79) &
       + int4(ri3_79,e(:,OP_DR ),temp79b,h(:,OP_GS ),weight_79,79) &
       - int4(ri3_79,e(:,OP_DZZ),temp79c,h(:,OP_DZ ),weight_79,79) &
       - int4(ri3_79,e(:,OP_DR ),temp79c,h(:,OP_DRZ),weight_79,79) &
       - int4(ri3_79,e(:,OP_DZ ),temp79c,h(:,OP_DZZ),weight_79,79) &
       - int4(ri3_79,e(:,OP_DRZ),temp79c,h(:,OP_DR ),weight_79,79) &
       - int4(ri3_79,e(:,OP_DRZ),temp79b,h(:,OP_DZ ),weight_79,79) &
       - int4(ri3_79,e(:,OP_DR ),temp79b,h(:,OP_DRR),weight_79,79) &
       - int4(ri3_79,e(:,OP_DZ ),temp79b,h(:,OP_DRZ),weight_79,79) &
       - int4(ri3_79,e(:,OP_DRR),temp79b,h(:,OP_DR ),weight_79,79)

  if(itor.eq.1) then
     ! [f,g]
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

     temp = temp &
          - int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_GS ),weight_79,79) &
          + int4(ri4_79,e(:,OP_DRZ),temp79a,h(:,OP_DZ ),weight_79,79) &
          + int4(ri4_79,e(:,OP_DZ ),temp79a,h(:,OP_DRZ),weight_79,79) &
          + int4(ri4_79,e(:,OP_DRR),temp79a,h(:,OP_DR ),weight_79,79) &
          + int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_DRR),weight_79,79)
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

  temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       - int5(ri3_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          int5(ri4_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_1),h(:,OP_1),weight_79,79)
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

  temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       - int5(ri3_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)

  if(itor.eq.1) then
     temp = temp - 2.* &
          int5(ri4_79,e(:,OP_GS),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
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

  ! <f,g>,r
  temp79b = f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DRZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRR) + f(:,OP_DZ )*g(:,OP_DRZ)

  ! <f,g>,z
  temp79c = f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRZ) + f(:,OP_DZ )*g(:,OP_DZZ)

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

  temp = int5(ri2_79,e(:,OP_GS),f(:,OP_GS),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
       + int5(ri2_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       + int5(ri2_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)

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
     temp = &
          - int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),weight_79,79) &
          + 0.5* &
          (int4(ri2_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),weight_79,79) &
          +int4(ri2_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),weight_79,79))

     if(itor.eq.1) then
        temp = temp - &
             (int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
             +int4(ri3_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),weight_79,79))
     endif
  else
     temp = &
          - int5(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          + 0.5* &
          (int5(ri2_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),h(:,OP_1),weight_79,79))

     if(itor.eq.1) then
        temp = temp - &
             (int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
             +int5(ri3_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79))
     endif
  endif

  v3uun = temp
  return
end function v3uun


! V3vvn
! =====
real function v3vvn(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h
  real :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     if(idens.eq.0) then
        temp = -int4(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_1),weight_79,79)
     else 
        temp = -int5(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_1),h(:,OP_1),weight_79,79)
     endif
  endif

  v3vvn = temp
  return
end function v3vvn


! V3uchin
! =======
real function v3uchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

  if(idens.eq.0) then
     temp = int3(ri_79,e(:,OP_LP),temp79a,weight_79,79) &
          + int4(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),weight_79,79) & 
          - int4(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),weight_79,79)
  else
     temp = int4(ri_79,e(:,OP_LP),temp79a,h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79) & 
          - int5(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) & 
          + int4(ri_79,e(:,OP_DZ),temp79a,h(:,OP_DZ),weight_79,79) & 
          + int4(ri_79,e(:,OP_DR),temp79a,h(:,OP_DR),weight_79,79)
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
     temp = 0.5* &
          (int3(e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),weight_79,79) &
          +int3(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),weight_79,79) &
          +int3(e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
          +int3(e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
          +int3(e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
          +int3(e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
          +int3(e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),weight_79,79) &
          +int3(e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),weight_79,79))
  else
     temp = 0.5* &
          (int4(e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          +int4(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),h(:,OP_1),weight_79,79) &
          +int4(e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          +int4(e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          +int4(e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          +int4(e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          +int4(e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          +int4(e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),h(:,OP_1),weight_79,79))
  endif

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

  temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
       - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),weight_79,79)

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

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  
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
real function b1psieta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real, intent(in) :: h
  real :: temp

  temp = -(int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1 ),weight_79,79) &
          +int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),weight_79,79) &
          +int3(e(:,OP_1 ),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          +int3(e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DZ),weight_79,79))
  if(itor.eq.1) then
     temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),weight_79,79)
  endif

  if(h.ne.0) then
     temp79a = e(:,OP_1)*g(:,OP_LP) + e(:,OP_LP)*g(:,OP_1) &
          + 2.*(e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))

!!$     if(itor.eq.1) temp79a = temp79a + 2.*ri_79* &
!!$          (e(:,OP_DR)*g(:,OP_1) + e(:,OP_1)*g(:,OP_DR))

     temp = temp - h*int2(temp79a,f(:,OP_GS),weight_79,79)
     if(itor.eq.1) temp = temp + 2.*h*int3(ri_79,temp79a,f(:,OP_DR),weight_79,79)
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
real function b2beta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real, intent(in) :: h
  real :: temp

  temp = -(int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          +int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79))

  if(itor.eq.1) then
     temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),weight_79,79)
  endif

  if(h.ne.0.) then
     temp79a = (e(:,OP_LP)*g(:,OP_1) + &
          e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))
     if(itor.eq.1) temp79a = temp79a + 2.*ri_79* &
          (e(:,OP_DR)*g(:,OP_1) + e(:,OP_1)*g(:,OP_DR))

     temp = temp - h*int2(temp79a,f(:,OP_GS),weight_79,79)
  endif

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

  if(itor.eq.1) then
     temp = temp - 2.*int4(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
  endif

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

  temp = int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
       + int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.*int4(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),weight_79,79)
  endif

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

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_1),weight_79,79)
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
     temp = int4(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.* &
             int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZ),weight_79,79)
     endif
  else
     temp = int5(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.* &
             int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     endif
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
             + 2.*int4(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
     endif
  else
     temp = int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + 2.*int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
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
     temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
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
! ===========
real function b3psipsieta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  
  real :: temp

  if(itor.eq.0) then
     b3psipsieta = 0.
     return
  endif

  temp = (gam-1.)*int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),h(:,OP_1),weight_79,79)

  b3psipsieta = temp
  
  return
end function b3psipsieta


! B3bbeta
! =======
real function b3bbeta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  
  real :: temp

  if(itor.eq.0) then
     b3bbeta = 0.
     return
  endif

  temp = (gam-1.)* &
       (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79))

  b3bbeta = temp
  
  return
end function b3bbeta


! B3pebd
! ======
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


! B3pedkappa
! ==========
real function b3pedkappa(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real, intent(in) :: h, i
  real :: temp

  if(idens.eq.0) then
     temp = -h* &
          (int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
          +int2(e(:,OP_DR),f(:,OP_DR),weight_79,79))
  else 
     temp = -h* &
          (int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),weight_79,79) &
          +int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1 ),weight_79,79) &
          +int3(e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),weight_79,79) &
          +int3(e(:,OP_DR),f(:,OP_1 ),g(:,OP_DR),weight_79,79))
  endif

  if(i.ne.0) then
     ! Laplacian[f g]
     temp79a = f(:,OP_LP)*g(:,OP_1) + f(:,OP_1)*g(:,OP_LP) &
          + 2.*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR))
     temp = temp - h*i*int2(temp79a,e(:,OP_LP),weight_79,79)
  endif

  b3pedkappa = temp  
  return
end function b3pedkappa



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
real function n1ndenm(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real, intent(in) :: g, h
  real :: temp

  temp = -g* &
       (int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       +int2(e(:,OP_DR),f(:,OP_DR),weight_79,79))

  if(h .ne. 0) then
     temp = temp - h*g*int2(e(:,OP_LP),f(:,OP_LP),weight_79,79)
  endif

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

  temp = int4(ri_79, e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79, e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

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

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

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


! P1kappar
! ========
real function p1kappar(e,f,g,h,i,j)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i,j
  real :: temp

  temp79a = ri2_79*(e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)

  if(idens.eq.0) then
     temp = int3(temp79a,g(:,OP_DZ),h(:,OP_DR),weight_79,79) &
          - int3(temp79a,g(:,OP_DR),h(:,OP_DZ),weight_79,79)
  else
     temp = int4(temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1 ),weight_79,79) &
          - int4(temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1 ),weight_79,79) &
          + int4(temp79a,g(:,OP_DZ),h(:,OP_1 ),i(:,OP_DR),weight_79,79) &
          - int4(temp79a,g(:,OP_DR),h(:,OP_1 ),i(:,OP_DZ),weight_79,79)
  endif

  p1kappar = temp
  return
end function p1kappar




!======================================================================
! Gyroviscous terms
!======================================================================

! g1u
! ===
real function g1u(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  temp79b = e(:,OP_DZZ) - e(:,OP_DRR)
  if(itor.eq.1) temp79b = temp79b - 3.0*ri_79*e(:,OP_DR)
  temp79c = f(:,OP_DRZ)
  if(itor.eq.1) temp79c = temp79c - 0.5*ri_79*f(:,OP_DZ)
  temp79d = f(:,OP_DZZ) - f(:,OP_DRR)
  if(itor.eq.1) temp79d = temp79d +     ri_79*f(:,OP_DR)
  temp79e = e(:,OP_DRZ)
  if(itor.eq.1) temp79e = temp79e + 1.5*ri_79*e(:,OP_DZ)

  temp79a = &
       ri_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
       (1.+1.5*ri2_79*b2i79(:,OP_1)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2)) * &
       (temp79b*temp79c - temp79d*temp79e)

  if(itor.eq.1) then 
     temp79a = temp79a &
          +2.25*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
           ((pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) * &
            (e(:,OP_DZ)*temp79d - f(:,OP_DZ)*temp79b) &
           +4.*pst79(:,OP_DR)*pst79(:,OP_DZ) * &
            (e(:,OP_DZ)*temp79c - f(:,OP_DZ)*temp79e))
  endif

  g1u = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g1u

! g1v
! ===
real function g1v(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  temp79b = f(:,OP_DR)
  if(itor.eq.1) temp79b = temp79b - 2.0*ri_79*f(:,OP_1)
  temp79c = e(:,OP_DRR) - e(:,OP_DZZ)
  if(itor.eq.1) temp79c = temp79c + 3.0*ri_79*e(:,OP_DR)
  temp79d = e(:,OP_DRZ)
  if(itor.eq.1) temp79d = temp79d + 3.0*ri_79*e(:,OP_DZ)

  temp79a = &
       0.25*ri_79*b2i79(:,OP_1)*(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) * &
        ((pst79(:,OP_DR)*f(:,OP_DZ) + pst79(:,OP_DZ)*temp79b)*temp79c &
         +2.*f(:,OP_DZ )*pst79(:,OP_DZ)*temp79d &
         -2.*e(:,OP_DRZ)*pst79(:,OP_DR)*temp79b)

  temp79b = pst79(:,OP_DZ)*f(:,OP_DR) - pst79(:,OP_DR)*f(:,OP_DZ)
  if(itor.eq.1) temp79b = temp79b - 2.*ri_79*f(:,OP_1)*pst79(:,OP_DZ)
  temp79c = e(:,OP_DZZ)
  if(itor.eq.1) temp79c = temp79c - 2.0*ri_79*e(:,OP_DR)
  temp79d = e(:,OP_DRZ)
  if(itor.eq.1) temp79d = temp79d + 1.5*ri_79*e(:,OP_DZ)
  temp79e = e(:,OP_DRR)
  if(itor.eq.1) temp79e = temp79e +     ri_79*e(:,OP_DR)

  temp79a = temp79a &
       -0.75*ri3_79*b2i79(:,OP_1)**2*temp79b * &
         (e(:,OP_GS)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2) - &
         2.*(   temp79c*pst79(:,OP_DR)**2 &
            -2.*temp79d*pst79(:,OP_DR)*pst79(:,OP_DZ) &
            +   temp79e*pst79(:,OP_DZ)**2))

  if(itor.eq.1) then
     temp79a = temp79a &
          + 4.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2*e(:,OP_DZ) * &
          (pst79(:,OP_DZ)*f(:,OP_DZ) + pst79(:,OP_DR)*f(:,OP_DR) &
          -2.*ri_79*f(:,OP_1)*pst79(:,OP_DR))
  endif
       
  g1v = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g1v

! g1chi
! =====
real function g1chi(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  temp79b = e(:,OP_DRZ)
  if(itor.eq.1) temp79b = temp79b + 3.0*ri_79*e(:,OP_DZ)
  temp79c = e(:,OP_DRR) - e(:,OP_DZZ)
  if(itor.eq.1) temp79c = temp79c + 3.0*ri_79*e(:,OP_DR)
  temp79d = f(:,OP_DZZ)
  if(itor.eq.1) temp79d = temp79d -     ri_79*f(:,OP_DR)
  temp79e = f(:,OP_DRR)
  if(itor.eq.1) temp79e = temp79e -     ri_79*f(:,OP_DR)

  temp79a = &
       0.5*b2i79(:,OP_1)*bzt79(:,OP_1) * &
        (temp79c*(f(:,OP_DRR) - f(:,OP_DZZ)) &
        +2.*(e(:,OP_DRZ) + temp79b)*f(:,OP_DRZ)) &
       +1.5*ri2_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
        (temp79c*(f(:,OP_GS)*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
        -f(:,OP_DZZ)*pst79(:,OP_DZ)**2 &
        +f(:,OP_DRR)*pst79(:,OP_DR)**2) &
        +2.*(f(:,OP_DRZ)*(pst79(:,OP_DZ)**2*e(:,OP_DRZ) &
                         +pst79(:,OP_DR)**2*temp79b) &
            -pst79(:,OP_DR)*pst79(:,OP_DZ)* &
             (e(:,OP_DRZ)*temp79d + temp79b*temp79e)))
       
  g1chi = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g1chi


! g2u
! ===
real function g2u(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
  if(itor.eq.1) temp79b = temp79b +     ri_79*f(:,OP_DR)
  temp79c = f(:,OP_DRZ)
  if(itor.eq.1) temp79c = temp79c +     ri_79*f(:,OP_DZ)
  temp79d = f(:,OP_DRZ)
  if(itor.eq.1) temp79d = temp79d - 2.0*ri_79*f(:,OP_DZ)
  temp79e = f(:,OP_DRZ)
  if(itor.eq.1) temp79e = temp79e - 0.5*ri_79*f(:,OP_DZ)

  temp79a = 0.25*ri_79*b2i79(:,OP_1)* &
       (1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)* &
        ((e(:,OP_DR)*pst79(:,OP_DZ) + e(:,OP_DZ)*pst79(:,OP_DR))*temp79b &
        -2.*e(:,OP_DZ)*pst79(:,OP_DZ)*temp79c &
        +2.*e(:,OP_DR)*pst79(:,OP_DR)*temp79d) &
       + 0.75*ri3_79*b2i79(:,OP_1)**2 * &
        (pst79(:,OP_DZ)*e(:,OP_DR) - pst79(:,OP_DR)*e(:,OP_DZ)) * &
        (   temp79b*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
        +4.*temp79e* pst79(:,OP_DR)    * pst79(:,OP_DZ)   )

  if(itor.eq.1) then
     temp79a = temp79a &
          -4.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2 * &
           f(:,OP_DZ)*(e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR))
  endif

  g2u = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g2u


! g2v
! ===
real function g2v(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  temp79b = e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ)
  if(itor.eq.1) temp79b = temp79b - 2.*ri_79*e(:,OP_DZ)*f(:,OP_1)

  temp79a = 0.25*ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)*temp79b* &
       (1.-3.*ri2_79*b2i79(:,OP_1)* &
        (pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 - bzt79(:,OP_1)**2))

  g2v = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g2v


! g2chi
! =====
real function g2chi(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  temp79b = f(:,OP_DZZ)
  if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)
  temp79c = f(:,OP_DRR)
  if(itor.eq.1) temp79c = temp79c - ri_79*f(:,OP_DR)

  temp79a = &
       -0.5*b2i79(:,OP_1)* &
        (e(:,OP_DZ)*temp79b*pst79(:,OP_DZ) &
        +e(:,OP_DR)*temp79c*pst79(:,OP_DR) &
        +f(:,OP_DRZ)*(e(:,OP_DR)*pst79(:,OP_DZ) + e(:,OP_DZ)*pst79(:,OP_DR))) &
       +1.5*ri2_79*b2i79(:,OP_1)**2 * &
        (pst79(:,OP_DZ)*e(:,OP_DR) - pst79(:,OP_DR)*e(:,OP_DZ)) * &
        (pst79(:,OP_DZ)*pst79(:,OP_DR)*(f(:,OP_DZZ) - f(:,OP_DRR)) &
        -(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*f(:,OP_DRZ)) &
       -1.5*ri2_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2 * &
        (e(:,OP_DZ)*temp79c*pst79(:,OP_DZ) &
        +e(:,OP_DR)*temp79b*pst79(:,OP_DR) &
        -f(:,OP_DRZ)*(e(:,OP_DR)*pst79(:,OP_DZ) + e(:,OP_DZ)*pst79(:,OP_DR)))

  g2chi = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g2chi


! g3u
! ===
real function g3u(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  if(itor.eq.1) then
     temp79b = ri_79*e(:,OP_DR)
     temp79c = ri_79*e(:,OP_DZ)
     temp79d = ri_79*f(:,OP_DR)
     temp79e = ri_79*f(:,OP_DZ)
  else
     temp79b = 0.
     temp79c = 0.
     temp79d = 0.
     temp79e = 0.
  endif

  temp79a = &
       0.5*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
        ((e(:,OP_DZZ) - e(:,OP_DRR))* &
         (f(:,OP_DZZ) - f(:,OP_DRR) + temp79d) &
        +4.*e(:,OP_DRZ)*(f(:,OP_DRZ) - 0.5*temp79e)) &
       +1.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
        ((f(:,OP_DZZ) - f(:,OP_DRR) + temp79d) * &
         (pst79(:,OP_DR)**2*(e(:,OP_DZZ) - temp79b) &
         -pst79(:,OP_DZ)**2*(e(:,OP_DRR) - temp79b)) &
        +2.*(e(:,OP_DRZ)*pst79(:,OP_DR)**2*(f(:,OP_DRZ) +    temp79e) &
            +e(:,OP_DRZ)*pst79(:,OP_DZ)**2*(f(:,OP_DRZ) - 2.*temp79e) &
            -pst79(:,OP_DR)*pst79(:,OP_DZ)* &
             ( f(:,OP_DRZ)           *(e(:,OP_DRR) - temp79b    ) &
             +(f(:,OP_DRZ) - temp79e)*(e(:,OP_DZZ) - temp79b    ) &
             + temp79e               *(e(:,OP_DRR) - e(:,OP_DZZ)))))

  g3u = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g3u


! g3v
! ===
real function g3v(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  if(itor.eq.1) then
     temp79b = ri_79*e(:,OP_DR)
     temp79c = ri_79*f(:,OP_1)
  else
     temp79b = 0.
     temp79c = 0.
  endif

  temp79a = &
       -0.5*ri2_79*b2i79(:,OP_1) * &
        (pst79(:,OP_DZ)*(e(:,OP_DZZ) - temp79b)* f(:,OP_DZ) &
        +pst79(:,OP_DR)*(e(:,OP_DRR) - temp79b)*(f(:,OP_DR) - 2.*temp79c) &
        +e(:,OP_DRZ)*(pst79(:,OP_DZ)*(f(:,OP_DR) - 2.*temp79c) &
                     +pst79(:,OP_DR)* f(:,OP_DZ))) &
       -1.5*ri4_79*b2i79(:,OP_1)**2 * &
        (f(:,OP_DZ)*pst79(:,OP_DR) - f(:,OP_DR)*pst79(:,OP_DZ) &
        +2.*temp79c*pst79(:,OP_DZ)) * &
        (pst79(:,OP_DZ)*pst79(:,OP_DR)*(e(:,OP_DZZ) - e(:,OP_DRR)) &
        -e(:,OP_DRZ)*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)) &
       +1.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2 * &
        (e(:,OP_DRZ)*(pst79(:,OP_DZ)*(f(:,OP_DR)-2.*temp79c) &
                     +pst79(:,OP_DR)* f(:,OP_DZ)) &
        -(e(:,OP_DZZ) - temp79b)*(f(:,OP_DR)-2.*temp79c)*pst79(:,OP_DR) &
        -(e(:,OP_DRR) - temp79b)* f(:,OP_DZ)            *pst79(:,OP_DZ))
       
  g3v = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g3v


! g3chi
! =====
real function g3chi(e,f)

  use basic
  use nintegrate_mod

  real, intent(in), dimension(79,OP_NUM) :: e,f

  if(itor.eq.1) then
     temp79b = ri_79*e(:,OP_DR)
     temp79c = ri_79*f(:,OP_DR)
  else
     temp79b = 0.
     temp79c = 0.
  endif

  temp79a = &
       -ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)* &
        ((e(:,OP_DZZ) - e(:,OP_DRR))*f(:,OP_DRZ) &
        -(f(:,OP_DZZ) - f(:,OP_DRR))*e(:,OP_DRZ)) &
       -3.*ri3_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
        (pst79(:,OP_DR)**2 * &
         ((e(:,OP_DZZ) - temp79b)*f(:,OP_DRZ) &
         -(f(:,OP_DZZ) - temp79c)*e(:,OP_DRZ)) &
        +pst79(:,OP_DZ)**2 * &
         ((f(:,OP_DRR) - temp79c)*e(:,OP_DRZ) &
         -(e(:,OP_DRR) - temp79b)*f(:,OP_DRZ)) &
        +pst79(:,OP_DZ)*pst79(:,OP_DR) * &
         ((e(:,OP_DRR) - temp79b)*(f(:,OP_DZZ) - temp79c) &
         -(e(:,OP_DZZ) - temp79b)*(f(:,OP_DRR) - temp79c)))
       
  g3chi = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g3chi





! ==============================================================
! Ohmic heating terms
! ==============================================================

! qpsipsieta
! ==========
real function qpsipsieta(e,hypf)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e
  real, intent(in) :: hypf
  real :: temp

  if(hypf.eq.0.) then
     temp = 0.
  else
     temp79a = ri2_79* &
          (jt79(:,OP_DZ)**2 + jt79(:,OP_DR)**2)
     if(itor.eq.1) then
        temp79a = temp79a - 2.*jt79(:,OP_1)* &
             (ri3_79*jt79(:,OP_DR) - ri4_79*jt79(:,OP_1))
     endif

     temp = hypf*int3(e(:,OP_1),eta79(:,OP_1),temp79a,weight_79,79)

     if(idens.eq.1) then
        temp79a = jt79(:,OP_DZ)*ni79(:,OP_DZ) + jt79(:,OP_DR)*ni79(:,OP_DR)
        if(itor.eq.1) then
           temp79a = temp79a - ri_79*jt79(:,OP_1)*ni79(:,OP_DR)
        endif
        temp79a = temp79a * ri2_79*nt79(:,OP_1)*jt79(:,OP_1)

        temp = temp + hypf*int3(e(:,OP_1),eta79(:,OP_1),temp79a,weight_79,79)
     endif
  endif

  qpsipsieta = temp
  return
end function qpsipsieta

! qbbeta
! ======
real function qbbeta(e,hypi)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e
  real, intent(in) :: hypi
  real :: temp

  if(hypi.eq.0.) then
     temp = 0.
  else
     temp79a = ri2_79*(bzt79(:,OP_GS)*bzt79(:,OP_GS) &
          + 2.*(bzt79(:,OP_DRZ)**2 - bzt79(:,OP_DRR)*bzt79(:,OP_DZZ)))

     if(itor.eq.1) then 
        temp79a = temp79a + 2.* &
             (ri3_79*bzt79(:,OP_DZZ)*bzt79(:,OP_DR) &
             -ri3_79*bzt79(:,OP_DRZ)*bzt79(:,OP_DZ) &
             +ri4_79*bzt79(:,OP_DZ )*bzt79(:,OP_DZ))
     endif

     temp = hypi*int3(e(:,OP_1),eta79(:,OP_1),temp79a,weight_79,79)

     if(idens.eq.1) then
        temp79a = &
               bzt79(:,OP_DZ)*bzt79(:,OP_DZZ)*ni79(:,OP_DZ) &
             + bzt79(:,OP_DR)*bzt79(:,OP_DRZ)*ni79(:,OP_DZ) &
             + bzt79(:,OP_DZ)*bzt79(:,OP_DRZ)*ni79(:,OP_DR) &
             + bzt79(:,OP_DR)*bzt79(:,OP_DRR)*ni79(:,OP_DR)
        if(itor.eq.1) then
           temp79a = temp79a - ri_79*ni79(:,OP_DR)* &
                (bzt79(:,OP_DZ)**2 + bzt79(:,OP_DR)**2)
        endif
        temp79a = temp79a * ri2_79*nt79(:,OP_1)

        temp = temp + hypi*int3(e(:,OP_1),eta79(:,OP_1),temp79a,weight_79,79)
     endif

  endif

  qbbeta = temp
  return
end function qbbeta


! ==============================================================
! Viscous heating terms
! ==============================================================

! quumu
! =====
real function quumu(e,f,g,h,i,j)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real, intent(in) :: h, i, j
  real :: temp

  temp = -h*int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),weight_79,79)

  if(j.ne.0) then
     temp = temp - h*j* &
          (int4(ri2_79,e(:,OP_1),vot79(:,OP_DZ),vot79(:,OP_DZ),weight_79,79) &
          +int4(ri2_79,e(:,OP_1),vot79(:,OP_DR),vot79(:,OP_DR),weight_79,79))
  endif

  quumu = temp
  return
end function quumu


! quchimu
! =======
real function quchimu(e,f,g,h,i,j)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real, intent(in) :: h, i, j
  real :: temp

  quchimu = 0.
  return
end function quchimu


! qvvmu
! =====
real function qvvmu(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real, intent(in) :: h, i
  real :: temp

  temp = -h* &
       (int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       +int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79))
  
  if(i.ne.1) then
     temp = temp - h*i*int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),weight_79,79)
  endif

  qvvmu = temp
  return
end function qvvmu


! qchichimu
! =========
real function qchichimu(e,f,g,h,i,j)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real, intent(in) :: h, i, j
  real :: temp

  temp = -2.*i*int4(ri2_79,e(:,OP_1),cot79(:,OP_1),cot79(:,OP_1),weight_79,79)
  
  if(j.ne.1) then
     temp = temp - 2.*i*j* &
          (int4(ri2_79,e(:,OP_1),cot79(:,OP_DZ),cot79(:,OP_DZ),weight_79,79) &
          +int4(ri2_79,e(:,OP_1),cot79(:,OP_DR),cot79(:,OP_DR),weight_79,79))
  endif

  qchichimu = temp
  return
end function qchichimu


!======================================================================
! FLUXES
!======================================================================

! Diffusive flux
! --------------
real function flux_diffusive

  use basic
  use nintegrate_mod

  if(idens.eq.0 .or. denm.eq.0.) then
     flux_diffusive = 0.
     return
  endif

  temp79a = ri2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2)
  if(numvar.ge.2) then
     temp79a = temp79a + ri2_79*vzt79(:,OP_1)**2
  endif
  if(numvar.ge.3) then
      temp79a = temp79a &
           + cht79(:,OP_DZ)**2 + cht79(:,OP_DR)**2 &
           + 2.*ri_79* &
           ( cht79(:,OP_DZ)*pht79(:,OP_DR) - cht79(:,OP_DR)*pht79(:,OP_DZ))
  endif

  flux_diffusive = 0.5*denm*int2(nt79(:,OP_LP),temp79a,weight_79,79)
  return
end function flux_diffusive


! Pressure convection
! -------------------
real function flux_pressure(dbf)

  use basic
  use nintegrate_mod

  real, intent(in) :: dbf

  if(numvar.lt.3 .or. gam.eq.1.) then
     flux_pressure = 0.
     return
  endif

  temp79a = pt79(:,OP_1)*cht79(:,OP_LP) &
       + pt79(:,OP_DZ)*cht79(:,OP_DZ) + pt79(:,OP_DR)*cht79(:,OP_DR) &
       + ri_79* &
       ( pt79(:,OP_DZ)*pht79(:,OP_DR) - pt79(:,OP_DR)*pht79(:,OP_DZ) &
       - pefac*dbf* &
         ( ni79(:,OP_1)*(pet79(:,OP_DZ)*bzt79(:,OP_DR) - pet79(:,OP_DR)*bzt79(:,OP_DZ)) &
         +pet79(:,OP_1)*( ni79(:,OP_DZ)*bzt79(:,OP_DR) -  ni79(:,OP_DR)*bzt79(:,OP_DZ))))

  flux_pressure = -gam/(gam-1.)*int1(temp79a,weight_79,79)
  return
end function flux_pressure


! Kinetic Energy Convection
! -------------------------
real function flux_ke

  use basic
  use nintegrate_mod

  real :: temp


  ! numvar = 1
  temp79a = ri2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2)
  if(idens.eq.1) then
     temp79b = ri_79*(nt79(:,OP_DZ)*pht79(:,OP_DR) - nt79(:,OP_DR)*pht79(:,OP_DZ))
  else
     temp79b = 0.
  endif
  temp79c = 0.5*ri3_79* &
       (pht79(:,OP_DR) * &
        (pht79(:,OP_DZZ)*pht79(:,OP_DZ ) + pht79(:,OP_DRZ)*pht79(:,OP_DR ) &
        +pht79(:,OP_DZ )*pht79(:,OP_DZZ) + pht79(:,OP_DR )*pht79(:,OP_DRZ)) &
       -pht79(:,OP_DZ) * &
        (pht79(:,OP_DRZ)*pht79(:,OP_DZ ) + pht79(:,OP_DRR)*pht79(:,OP_DR ) &
        +pht79(:,OP_DZ )*pht79(:,OP_DRZ) + pht79(:,OP_DR )*pht79(:,OP_DRR)))
  if(itor.eq.1) then
     temp79c = temp79c + ri4_79*pht79(:,OP_DZ)* &
       (pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2)
  endif

  ! numvar = 2
  if(numvar.ge.2) then
     temp79a = temp79a + ri2_79*vzt79(:,OP_1)**2
     temp79c = temp79c + &
          ri3_79*vzt79(:,OP_1) * &
          (vzt79(:,OP_DZ)*pht79(:,OP_DR) - vzt79(:,OP_DR)*pht79(:,OP_DZ))
     if(itor.eq.1) then
        temp79c = temp79c + ri4_79*vzt79(:,OP_1)**2*pht79(:,OP_DZ)
     endif
  endif

  ! numvar = 3
  if(numvar.ge.3) then
      temp79a = temp79a &
           + cht79(:,OP_DZ)*cht79(:,OP_DZ) + cht79(:,OP_DR)*cht79(:,OP_DR) &
           + 2.*ri_79* &
           ( cht79(:,OP_DZ)*pht79(:,OP_DR) - cht79(:,OP_DR)*pht79(:,OP_DZ))
      temp79b = temp79b + nt79(:,OP_1)*cht79(:,OP_LP) &
            + nt79(:,OP_DZ)*cht79(:,OP_DZ) + nt79(:,OP_DR)*cht79(:,OP_DR)
      temp79c = temp79c + 0.5*ri_79* &
           (pht79(:,OP_DR) * &
            (cht79(:,OP_DZZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRZ)*cht79(:,OP_DR ) &
            +cht79(:,OP_DZ )*cht79(:,OP_DZZ) + cht79(:,OP_DR )*cht79(:,OP_DRZ)) &
           -pht79(:,OP_DZ) * &
            (cht79(:,OP_DRZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRR)*cht79(:,OP_DR ) &
            +cht79(:,OP_DZ )*cht79(:,OP_DRZ) + cht79(:,OP_DR )*cht79(:,OP_DRR))) &
           + ri_79* &
           (cht79(:,OP_DZ) * &
            (cht79(:,OP_DZZ)*pht79(:,OP_DR ) - cht79(:,OP_DRZ)*pht79(:,OP_DZ ) &
            +cht79(:,OP_DZ )*pht79(:,OP_DRZ) - cht79(:,OP_DR )*pht79(:,OP_DZZ)) &
           +cht79(:,OP_DR) * &
            (cht79(:,OP_DRZ)*pht79(:,OP_DR ) - cht79(:,OP_DRR)*pht79(:,OP_DZ ) &
            +cht79(:,OP_DZ )*pht79(:,OP_DRR) - cht79(:,OP_DR )*pht79(:,OP_DRZ))) &
           + ri2_79* &
           (pht79(:,OP_DR) * &
            (cht79(:,OP_DZZ)*pht79(:,OP_DR ) - cht79(:,OP_DRZ)*pht79(:,OP_DZ ) &
            +cht79(:,OP_DZ )*pht79(:,OP_DRZ) - cht79(:,OP_DR )*pht79(:,OP_DZZ)) &
           -pht79(:,OP_DZ) * &
            (cht79(:,OP_DRZ)*pht79(:,OP_DR ) - cht79(:,OP_DRR)*pht79(:,OP_DZ ) &
            +cht79(:,OP_DZ )*pht79(:,OP_DRR) - cht79(:,OP_DR )*pht79(:,OP_DRZ))) &
           + 0.5*ri2_79* &
           (cht79(:,OP_DZ) * &
            (pht79(:,OP_DZZ)*pht79(:,OP_DZ ) + pht79(:,OP_DRZ)*pht79(:,OP_DR ) &
            +pht79(:,OP_DZ )*pht79(:,OP_DZZ) + pht79(:,OP_DR )*pht79(:,OP_DRZ)) &
           +cht79(:,OP_DR) * &
            (pht79(:,OP_DRZ)*pht79(:,OP_DZ ) + pht79(:,OP_DRR)*pht79(:,OP_DR ) &
            +pht79(:,OP_DZ )*pht79(:,OP_DRZ) + pht79(:,OP_DR )*pht79(:,OP_DRR))) &
           + ri2_79*vzt79(:,OP_1)* &
            (vzt79(:,OP_DZ )*cht79(:,OP_DZ ) + vzt79(:,OP_DR )*cht79(:,OP_DR )) &
           + 0.5* &
           (cht79(:,OP_DZ) * &
            (cht79(:,OP_DZZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRZ)*cht79(:,OP_DR ) &
            +cht79(:,OP_DZ )*cht79(:,OP_DZZ) + cht79(:,OP_DR )*cht79(:,OP_DRZ)) &
           +cht79(:,OP_DR) * &
            (cht79(:,OP_DRZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRR)*cht79(:,OP_DR ) &
            +cht79(:,OP_DZ )*cht79(:,OP_DRZ) + cht79(:,OP_DR )*cht79(:,OP_DRR)))
      
      if(itor.eq.1) then
         temp79c = temp79c &
              + ri2_79*cht79(:,OP_DR)* &
              (pht79(:,OP_DZ)*cht79(:,OP_DR) - pht79(:,OP_DR)*cht79(:,OP_DZ)) &
              - ri3_79*pht79(:,OP_DZ)* &
              (pht79(:,OP_DZ)*cht79(:,OP_DR) - pht79(:,OP_DR)*cht79(:,OP_DZ)) &
              - ri3_79*cht79(:,OP_DR)* &
              (pht79(:,OP_DZ)*pht79(:,OP_DZ) + pht79(:,OP_DR)*pht79(:,OP_DR)) &
              - ri3_79*vzt79(:,OP_1)**2*cht79(:,OP_DR)
      endif
  endif

  temp = 0.5*int2(temp79a,temp79b,weight_79,79)

  if(idens.eq.1) then
     temp = temp + int2(nt79(:,OP_1),temp79c,weight_79,79)
  else
     temp = temp + int1(temp79c,weight_79,79)
  endif

  flux_ke = -temp
  return
end function flux_ke


! Poynting flux
! -------------
real function flux_poynting(dbf)

  use basic
  use nintegrate_mod

  real, intent(in) :: dbf

  if(idens.eq.0) then
     ni79 = 0.
     ni79(:,OP_1) = 1.
  endif

  temp79a = ri3_79*(pst79(:,OP_GS)* &
        (pst79(:,OP_DZ )*pht79(:,OP_DR ) - pst79(:,OP_DR )*pht79(:,OP_DZ )) &
       +pst79(:,OP_DZ)* &
        (pst79(:,OP_DZZ)*pht79(:,OP_DR ) - pst79(:,OP_DRZ)*pht79(:,OP_DZ ) &
        +pst79(:,OP_DZ )*pht79(:,OP_DRZ) - pst79(:,OP_DR )*pht79(:,OP_DZZ)) &
       +pst79(:,OP_DR)* &
        (pst79(:,OP_DRZ)*pht79(:,OP_DR ) - pst79(:,OP_DRR)*pht79(:,OP_DZ) &
        +pst79(:,OP_DZ )*pht79(:,OP_DRR) - pst79(:,OP_DR )*pht79(:,OP_DRZ))) &
       - ri2_79* &
       +(eta79(:,OP_1 )*  jt79(:,OP_1 )* jt79(:,OP_1 ) &
        +eta79(:,OP_1 )*(pst79(:,OP_DZ)* jt79(:,OP_DZ) + pst79(:,OP_DR)* jt79(:,OP_DR)) &
        + jt79(:,OP_1 )*(pst79(:,OP_DZ)*eta79(:,OP_DZ) + pst79(:,OP_DR)*eta79(:,OP_DR)))
  
  if(itor.eq.1) then
     temp79a = temp79a - ri4_79*pst79(:,OP_DR)* &
          (pst79(:,OP_DZ)*pht79(:,OP_DR) - pst79(:,OP_DR)*pht79(:,OP_DZ))
  endif


  if(numvar.ge.2) then
     temp79a = temp79a &
          - ri2_79* &
          (eta79(:,OP_1)* bzt79(:,OP_1)*bzt79(:,OP_GS) &
          +eta79(:,OP_1)*(bzt79(:,OP_DZ)*bzt79(:,OP_DZ) + bzt79(:,OP_DR)*bzt79(:,OP_DR)) &
          +bzt79(:,OP_1)*(bzt79(:,OP_DZ)*eta79(:,OP_DZ) + bzt79(:,OP_DR)*eta79(:,OP_DR))) &
          - ri3_79* &
          (vzt79(:,OP_1)*(bzt79(:,OP_DZ)*pst79(:,OP_DR) - bzt79(:,OP_DR)*pst79(:,OP_DZ)) &
          +bzt79(:,OP_1)*(vzt79(:,OP_DZ)*pst79(:,OP_DR) - vzt79(:,OP_DR)*pst79(:,OP_DZ))) &
          - ri3_79*dbf* &
           (bzt79(:,OP_1)* &
           ( jt79(:,OP_1 )*(ni79(:,OP_DZ)*pst79(:,OP_DR) - ni79(:,OP_DR)*pst79(:,OP_DZ)) &
           + ni79(:,OP_1 )*(jt79(:,OP_DZ)*pst79(:,OP_DR) - jt79(:,OP_DR)*pst79(:,OP_DZ))) &
           +( ni79(:,OP_DZ)*pst79(:,OP_DZ) +  ni79(:,OP_DR)*pst79(:,OP_DR))* &
            (pst79(:,OP_DZ)*bzt79(:,OP_DR) - pst79(:,OP_DR)*bzt79(:,OP_DZ)) &
           +ni79(:,OP_1)* &
           (pst79(:,OP_DZ)* &
            (pst79(:,OP_DZZ)*bzt79(:,OP_DR ) - pst79(:,OP_DRZ)*bzt79(:,OP_DZ ) &
            +pst79(:,OP_DZ )*bzt79(:,OP_DRZ) - pst79(:,OP_DR )*bzt79(:,OP_DZZ)) &
           +pst79(:,OP_DR)* &
            (pst79(:,OP_DRZ)*bzt79(:,OP_DR ) - pst79(:,OP_DRR)*bzt79(:,OP_DZ ) &
            +pst79(:,OP_DZ )*bzt79(:,OP_DRR) - pst79(:,OP_DR )*bzt79(:,OP_DRZ))) &
           +bzt79(:,OP_1)**2*(ni79(:,OP_DZ)*bzt79(:,OP_DR) - ni79(:,OP_DR)*bzt79(:,OP_DZ))) &
          +2.*ri3_79*bzt79(:,OP_1)* &
           (bzt79(:,OP_DZ)*pht79(:,OP_DR) - bzt79(:,OP_DR)*pht79(:,OP_DZ))

     if(itor.eq.1) then
        temp79a = temp79a &
             - 2.*ri4_79*pst79(:,OP_DZ)*vzt79(:,OP_1)*bzt79(:,OP_1) &
             +dbf*ri4_79*ni79(:,OP_1)* &
              (pst79(:,OP_DR)*(pst79(:,OP_DZ)*bzt79(:,OP_DR) - pst79(:,OP_DR)*bzt79(:,OP_DZ)) &
              -2.*pst79(:,OP_DZ)*pst79(:,OP_GS)*bzt79(:,OP_1) &
              -2.*bzt79(:,OP_1)**2*bzt79(:,OP_DZ)) &
             +2.*ri4_79*bzt79(:,OP_1)**2*pht79(:,OP_DZ)
     endif
  endif

  if(numvar.ge.3) then
     temp79a = temp79a &
          - ri_79*pefac*dbf* &
            (bzt79(:,OP_1)*( ni79(:,OP_DZ)*pet79(:,OP_DR) -  ni79(:,OP_DR)*pet79(:,OP_DZ)) &
            + ni79(:,OP_1)*(bzt79(:,OP_DZ)*pet79(:,OP_DR) - bzt79(:,OP_DR)*pet79(:,OP_DZ))) &
          + ri2_79*bzt79(:,OP_1)* &
            (bzt79(:,OP_1)*cht79(:,OP_GS) + &
            2.*(bzt79(:,OP_DZ)*cht79(:,OP_DZ) + bzt79(:,OP_DR)*cht79(:,OP_DR))) &
          + ri2_79* &
            (pst79(:,OP_GS)*(pst79(:,OP_DZ )*cht79(:,OP_DZ ) + pst79(:,OP_DR )*cht79(:,OP_DR )) &
            +pst79(:,OP_DZ)*(pst79(:,OP_DZZ)*cht79(:,OP_DZ ) + pst79(:,OP_DRZ)*cht79(:,OP_DR ) &
                            +pst79(:,OP_DZ )*cht79(:,OP_DZZ) + pst79(:,OP_DR )*cht79(:,OP_DRZ)) &
            +pst79(:,OP_DR)*(pst79(:,OP_DRZ)*cht79(:,OP_DZ ) + pst79(:,OP_DRR)*cht79(:,OP_DR ) &
                            +pst79(:,OP_DZ )*cht79(:,OP_DRZ) + pst79(:,OP_DR )*cht79(:,OP_DRR)))
          
  endif

!!$ MISSING HYPER-RESISTIVE TERMS

  flux_poynting = -int1(temp79a,weight_79,79)
  return
end function flux_poynting


! Heat flux
! ---------
real function flux_heat

  use basic
  use nintegrate_mod

  real :: temp

  if(numvar.lt.3) then
     flux_heat = 0.
     return
  endif

  ! Isotropic heat flux
  if(idens.eq.0) then
     temp79a = kappat*pt79(:,OP_LP)
  else
     temp79a = kappat*(ni79(:,OP_1)*pt79(:,OP_LP) + ni79(:,OP_LP)*pt79(:,OP_1) &
          + 2.*(ni79(:,OP_DZ)*pt79(:,OP_DZ) + ni79(:,OP_DR)*pt79(:,OP_DR)))
  endif

  temp = int1(temp79a,weight_79,79)


  ! Parallel heat flux
  if(kappar.ne.0.) then

     temp79c = kappar*(b2i79(:,OP_DZ)*pst79(:,OP_DR) - b2i79(:,OP_DR)*pst79(:,OP_DZ))
     if(itor.eq.1) temp79c = temp79c + ri_79*kappar*b2i79(:,OP_1)*pst79(:,OP_DZ)

     if(idens.eq.0) then
        temp79b = pt79(:,OP_DZ)*pst79(:,OP_DR) - pt79(:,OP_DR)*pst79(:,OP_DZ)

        temp79d =-pst79(:,OP_DZ)* &
              (pt79(:,OP_DRZ)*pst79(:,OP_DR ) - pt79(:,OP_DRR)*pst79(:,OP_DZ ) &
              +pt79(:,OP_DZ )*pst79(:,OP_DRR) - pt79(:,OP_DR )*pst79(:,OP_DRZ)) &
                 +pst79(:,OP_DR)* &
              (pt79(:,OP_DZZ)*pst79(:,OP_DR ) - pt79(:,OP_DRZ)*pst79(:,OP_DZ ) &
              +pt79(:,OP_DZ )*pst79(:,OP_DRZ) - pt79(:,OP_DR )*pst79(:,OP_DZZ))
     else
        temp79b = ni79(:,OP_1)*(pt79(:,OP_DZ)*pst79(:,OP_DR) - pt79(:,OP_DR)*pst79(:,OP_DZ)) &
             +    pt79(:,OP_1)*(ni79(:,OP_DZ)*pst79(:,OP_DR) - ni79(:,OP_DR)*pst79(:,OP_DZ))
     
        temp79d = -ni79(:,OP_1)* &
             (pst79(:,OP_DZ)* &
              (pt79(:,OP_DRZ)*pst79(:,OP_DR ) - pt79(:,OP_DRR)*pst79(:,OP_DZ ) &
              +pt79(:,OP_DZ )*pst79(:,OP_DRR) - pt79(:,OP_DR )*pst79(:,OP_DRZ)) &
             -pst79(:,OP_DR)* &
              (pt79(:,OP_DZZ)*pst79(:,OP_DR ) - pt79(:,OP_DRZ)*pst79(:,OP_DZ ) &
              +pt79(:,OP_DZ )*pst79(:,OP_DRZ) - pt79(:,OP_DR )*pst79(:,OP_DZZ))) &
             - pt79(:,OP_1)* &
             (pst79(:,OP_DZ)* &
              (ni79(:,OP_DRZ)*pst79(:,OP_DR ) - ni79(:,OP_DRR)*pst79(:,OP_DZ ) &
              +ni79(:,OP_DZ )*pst79(:,OP_DRR) - ni79(:,OP_DR )*pst79(:,OP_DRZ)) &
             -pst79(:,OP_DR)* &
              (ni79(:,OP_DZZ)*pst79(:,OP_DR ) - ni79(:,OP_DRZ)*pst79(:,OP_DZ ) &
              +ni79(:,OP_DZ )*pst79(:,OP_DRZ) - ni79(:,OP_DR )*pst79(:,OP_DZZ))) &
             + 2.*(pt79(:,OP_DZ)*pst79(:,OP_DR) - pt79(:,OP_DR)*pst79(:,OP_DZ))* &
                  (ni79(:,OP_DZ)*pst79(:,OP_DR) - ni79(:,OP_DR)*pst79(:,OP_DZ))
     endif

     temp = temp &
          + int3(ri2_79,temp79b,temp79c,weight_79,79) &
          + kappar*int3(ri2_79,b2i79(:,OP_1),temp79d,weight_79,79)
  endif

!!$ MISSING HYPER-DIFFUSIVE TERMS

  flux_heat = temp
  return
end function flux_heat

end module metricterms_new






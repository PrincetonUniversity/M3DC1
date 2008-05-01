module metricterms_new

implicit none

contains      

!=============================================================================
! V1 TERMS
!=============================================================================
  

! V1umu 
! =====
vectype function v1umu(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)

     temp79a = e(:,OP_GS)*g(:,OP_1) &
          + e(:,OP_DR)*g(:,OP_DR) + e(:,OP_DZ)*g(:,OP_DZ)
     if(itor.eq.1) then 
        temp79a = temp79a  &
             + 4.*ri_79*e(:,OP_DR)*g(:,OP_1 ) &
             + 2.*ri_79*e(:,OP_1) *g(:,OP_DR)
     endif

     temp79b = e(:,OP_DR)*f(:,OP_DR) + e(:,OP_DZ)*f(:,OP_DZ)
     temp79c = g(:,OP_DR)*f(:,OP_DR) + g(:,OP_DZ)*f(:,OP_DZ)

     temp = &
          - int2(temp79a,f(:,OP_GS),weight_79,79) &
          - int2(temp79b,g(:,OP_GS),weight_79,79) &
          - int2(temp79c,e(:,OP_GS),weight_79,79)

     if(itor.eq.1) then
        temp = temp &
             - 2.*int3(ri_79,temp79b,g(:,OP_DR),weight_79,79) &
             - 4.*int3(ri_79,temp79c,e(:,OP_DR),weight_79,79) &
             - 2.*int4(ri_79,e(:,OP_1),g(:,OP_LP),f(:,OP_DR),weight_79,79)
     endif

#ifdef USECOMPLEX
     temp = temp + &
          (int4(ri2_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),weight_79,79) &
          +int4(ri2_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),weight_79,79))
     if(itor.eq.1) then
        temp = temp + &
             2.*int4(ri3_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_1),weight_79,79)
     endif
#endif

     
  case(1)

     temp79a = 4.*e(:,OP_DRZ)*f(:,OP_DRZ) &
          +(e(:,OP_DZZ) - e(:,OP_DRR))*(f(:,OP_DZZ) - f(:,OP_DRR))

     if(itor.eq.1) then
        temp79a = temp79a - ri_79* &
             (e(:,OP_DR)*(f(:,OP_DZZ) - f(:,OP_DRR)) &
             +f(:,OP_DR)*(e(:,OP_DZZ) - e(:,OP_DRR)) &
             -2.*(e(:,OP_DZ)*f(:,OP_DRZ) + &
                  f(:,OP_DZ)*e(:,OP_DRZ))) &
             + ri2_79*(e(:,OP_DR)*f(:,OP_DR) - 4.*e(:,OP_DZ)*f(:,OP_DZ))
     endif

     temp = -int3(r2_79,temp79a,g(:,OP_1),weight_79,79)
     
     if(itor.eq.1) then
        temp = temp - &
             8.*int3(e(:,OP_DZ),f(:,OP_DZ),h(:,OP_1),weight_79,79)
     endif

  end select

  v1umu = temp
  return
end function v1umu



! V1chimu
! =======
vectype function v1chimu(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     temp79a = g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ)
     temp79b = r_79*(e(:,OP_DZ)*g(:,OP_DR) - e(:,OP_DR)*g(:,OP_DZ))
     if(itor.eq.1) temp79b = temp79b - 2.*e(:,OP_1)*g(:,OP_DZ)

     temp = int3(r_79,e(:,OP_LP),temp79a,weight_79,79) &
          + int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS),weight_79,79) &
          - int4(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS),weight_79,79) &
          + int2(temp79b,f(:,OP_GS),weight_79,79)

     if(itor.eq.1) then
        temp = temp &
             + 4.*int2(      e(:,OP_DR),temp79a,weight_79,79) &
             + 4.*int3(ri_79,e(:,OP_1 ),temp79a,weight_79,79) &
             + 4.*int3(ri_79,f(:,OP_DR),temp79b,weight_79,79) &
             - 2.*int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_GS),weight_79,79)
     endif

  case(1)

     temp79a = e(:,OP_DRZ)*(f(:,OP_DZZ) - f(:,OP_DRR)) &
          -    f(:,OP_DRZ)*(e(:,OP_DZZ) - e(:,OP_DRR))

     if(itor.eq.1) then
        temp79a = temp79a - ri_79* &
             (e(:,OP_DZ)*f(:,OP_DRR) - e(:,OP_DR)*f(:,OP_DRZ) &
             -2.*e(:,OP_DRZ)*f(:,OP_DR) &
             -f(:,OP_DZ)*(e(:,OP_DZZ) - e(:,OP_DRR))) &
             - ri2_79* &
             (e(:,OP_DR)*f(:,OP_DZ) - e(:,OP_DZ)*f(:,OP_DR))
     endif

     temp = temp - 2.*int3(ri_79,temp79a,g(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp = temp &
             +4.*int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),h(:,OP_1),weight_79,79) &
             -4.*int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_1),weight_79,79)
     endif
     
  end select

  v1chimu = temp
  return
end function v1chimu




! V1un
! ====
vectype function v1un(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),weight_79,79)
     endif

  case(1)
     temp = (int4(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
          +  int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79))/rzero**2
  end select

  v1un = temp
  return
end function v1un


! V1chin
! ======
vectype function v1chin(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     
  case(1)
     temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),weight_79,79)
  end select

  v1chin = temp
  return
end function v1chin


! V1psipsi
! ========
vectype function v1psipsi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
          
  temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS),weight_79,79) &
       - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS),weight_79,79)

  select case(ivform)
  case(0)
     if(itor.eq.1) then
        temp = temp - 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_GS),weight_79,79)
     endif
  case(1)
     
  end select
  
  v1psipsi = temp
  return
end function v1psipsi


! V1psib
! ======
vectype function v1psib(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
#ifdef USECOMPLEX
  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1),weight_79,79)

  select case(ivform)
  case(0)
  if(itor.eq.1) then
     temp = temp &
          + 2.*int4(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_1),weight_79,79)
  endif
  case(1)
!
  end select

  v1psib = temp
#else
  v1psib = 0.
#endif
end function v1psib


! V1bb
! ====
vectype function v1bb(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(itor.eq.0) then
        temp = 0.
     else
        temp = -2.*int4(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79) 
     endif

  case(1)

     temp = 0.

!!$     CHANGE: NMF 4/30/08
!!$     temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),weight_79,79) &
!!$          - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  end select

  v1bb = temp
  return
end function v1bb


! V1uun 
! =====
vectype function v1uun(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp
           

  select case(ivform)
  case(0)
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

  case(1)
     temp79a = e(:,OP_DZ)*(f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ)) &
             + e(:,OP_DR)*(f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ))
     temp = -int3(r3_79,temp79a,h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp = temp &
             + int5(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
             + int5(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     end if

  end select

  v1uun = temp
  return
end function v1uun


! v1vvn
! =====
vectype function v1vvn(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     select case(ivform)
     case(0)
        temp = -int5(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_1),h(:,OP_1),weight_79,79)

     case(1)
        temp = -int5(r2_79, e(:,OP_DZ),f(:,OP_1),g(:,OP_1),h(:,OP_1),weight_79,79)
     end select
  endif

  v1vvn = temp
  return
end function v1vvn


! v1uchin
! =======
vectype function v1uchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  vectype :: temp

  select case(ivform)
  case(0)
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

  case(1)
     temp79a = e(:,OP_DR)*(f(:,OP_DR)*g(:,OP_DZZ) + f(:,OP_DRZ)*g(:,OP_DZ) &
                          -f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DRR)*g(:,OP_DR)) &
             + e(:,OP_DR)*(f(:,OP_DZ)*g(:,OP_DRR) + f(:,OP_DRZ)*g(:,OP_DR) &
                          -f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZZ)*g(:,OP_DZ))
     temp = -int2(temp79a,h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp79a = e(:,OP_DR)*(f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ)) &
                + f(:,OP_DZ)*(e(:,OP_DR)*g(:,OP_DZ) - e(:,OP_DZ)*g(:,OP_DR))
        temp = temp - int3(ri_79,temp79a,h(:,OP_1),weight_79,79)
     end if
     
  end select

  v1uchin = temp
end function v1uchin


! v1chichin
! =========
vectype function v1chichin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  vectype :: temp

  select case(ivform)
  case(0)
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

     temp = 0.5* &
          (int4(r_79,e(:,OP_DR),temp79a,h(:,OP_DZ),weight_79,79) &
          -int4(r_79,e(:,OP_DZ),temp79a,h(:,OP_DR),weight_79,79))
     
     if(itor.eq.1) then
        temp = temp + int3(e(:,OP_1),temp79a,h(:,OP_DZ),weight_79,79)
     endif

  case(1)

     temp79a = g(:,OP_DZ)*(e(:,OP_DR)*f(:,OP_DZZ) - e(:,OP_DZ)*f(:,OP_DRZ)) &
             + g(:,OP_DR)*(e(:,OP_DR)*f(:,OP_DRZ) - e(:,OP_DZ)*f(:,OP_DRR))

     temp = -int3(ri3_79,temp79a,h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp79a = g(:,OP_DR)*(e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))

        temp = temp - 2.*int3(ri4_79,temp79a,h(:,OP_1),weight_79,79)
     end if

  end select
  
  v1chichin = temp

  return
end function v1chichin



! V1upsipsi
! =========
vectype function v1upsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

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

  select case(ivform)
  case(0)  

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

  case(1)
     temp = -int2(temp79b,temp79e,weight_79,79)/rzero**2 &
          - int2(temp79c,temp79f,weight_79,79)/rzero**2 &
          + int3(e(:,OP_DZ),temp79b,h(:,OP_GS),weight_79,79)/rzero**2 &
          - int3(e(:,OP_DR),temp79c,h(:,OP_GS),weight_79,79)/rzero**2

     if(itor.eq.1) then
        ! |u, psi(1)|
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
        
        ! |nu,psi(2)|
        temp79d = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
        temp = temp             &
             - int3(ri_79,temp79d,temp79b,weight_79,79)/rzero**2  &
             - int3(ri_79,temp79a,temp79e,weight_79,79)/rzero**2  &
             - int3(ri2_79,temp79a,temp79d,weight_79,79)/rzero**2  &
             + int4(ri_79,h(:,OP_GS),temp79a,e(:,OP_DZ),weight_79,79)/rzero**2 
      endif

  end select


  v1upsipsi = temp
  return
end function v1upsipsi

! V1upsib
! =======
vectype function v1upsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp
#ifdef USECOMPLEX

  temp79a = g(:,OP_DR)*(e(:,OP_DZZ)*f(:,OP_DZP)-e(:,OP_DZ)*f(:,OP_DZZP) &
                       +e(:,OP_DRZ)*f(:,OP_DRP)-e(:,OP_DR)*f(:,OP_DRZP)) &
           -g(:,OP_DZ)*(e(:,OP_DRR)*f(:,OP_DRP)-e(:,OP_DR)*f(:,OP_DRRP) &
                       +e(:,OP_DRZ)*f(:,OP_DZP)-e(:,OP_DZ)*f(:,OP_DRZP))
  select case(ivform)
  case(0)

    temp = -int3(ri3_79,temp79a,h(:,OP_1),weight_79,79)

    if(itor.eq.1) then
       temp79a = g(:,OP_DZ)*(e(:,OP_1)*f(:,OP_DRRP) - e(:,OP_DZ)*f(:,OP_DZP) &
                                               - 2.*e(:,OP_DR)*f(:,OP_DRP)) &
              +g(:,OP_DR)*(e(:,OP_DZ)*f(:,OP_DRP) - e(:,OP_1 )*f(:,OP_DRZP))
       temp = temp &
          - 2.*int3(ri4_79,temp79a,h(:,OP_1),weight_79,79) &
          + 2.*int5(ri5_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1),weight_79,79)
    endif

  case(1)
    temp = -int3(ri_79,temp79a,h(:,OP_1),weight_79,79)/rzero**2
    if(itor.eq.1) then
       temp79a = f(:,OP_DZP)*g(:,OP_DR)-f(:,OP_DRP)*g(:,OP_DZ)
       temp = temp + (int5(ri2_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_DR),h(:,OP_1),weight_79,79) &
                     -int5(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1),weight_79,79))/rzero**2

    endif

  end select

  v1upsib = temp
#else
  v1upsib = 0.
#endif

end function v1upsib


! V1ubb 
! =====
vectype function v1ubb(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
  select case(ivform)
  case(0)
     if(itor.eq.1) then
        temp = 2.* &
             (   int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
             -   int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
             -2.*int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1),weight_79,79))
     endif

#ifdef USECOMPLEX
     temp = temp + &
          (int5(ri4_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),h(:,OP_1),weight_79,79) &
          +int5(ri4_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1),weight_79,79))
     if(itor.eq.1) then
        temp = temp + &
             2.*int5(ri5_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1),weight_79,79)
     endif

#endif
  case(1)

!!$ CHANGED NMF 4/30/08
!!$     if(itor.eq.1) then
!!$        temp = 2.* &
!!$             (   int5(ri_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
!!$             -   int5(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79))/rzero**2
!!$     endif

#ifdef USECOMPLEX
     temp = temp + &
          (int5(ri2_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1),weight_79,79))/rzero**2

#endif

  end select
  v1ubb = temp
  return
end function v1ubb


! V1up
! ====
vectype function v1up(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = 0.
     
  case(1)

     temp = 0.
     if(itor.eq.1) then 
        temp = 2.* &
             (int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
             -int4(r_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),weight_79,79)) &
             -4.*gam*int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
     end if

  end select

  v1up = temp
  return
end function v1up



! V1vpsipsi
! =========
vectype function v1vpsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX

  select case(ivform)
  case(0)
     temp79a = g(:,OP_DZ)*( h(:,OP_DZ)*e(:,OP_DRZ) + h(:,OP_DR)*e(:,OP_DRR)  &
                          - e(:,OP_DZ)*h(:,OP_DRZ) - e(:,OP_DR)*h(:,OP_DRR)) &
             + g(:,OP_DR)*( e(:,OP_DR)*h(:,OP_DRZ) + e(:,OP_DZ)*h(:,OP_DZZ)  &
                          - h(:,OP_DR)*e(:,OP_DRZ) - h(:,OP_DZ)*e(:,OP_DZZ))
     temp = -int3(ri3_79,f(:,OP_DP),temp79a,weight_79,79)

     if(itor.eq.1) then
        temp79a = e(:,OP_DZ)*(g(:,OP_DZ)*h(:,OP_DZ ) + g(:,OP_DR)*h(:,OP_DR)) &
             + 2.*e(:,OP_DR)*(g(:,OP_DZ)*h(:,OP_DR ))                         &
             +    e(:,OP_1) *(g(:,OP_DR)*h(:,OP_DZZ) - g(:,OP_DZ)*h(:,OP_DRR))
        temp = temp &
             -2.*int3(ri4_79,f(:,OP_DP),temp79a,weight_79,79) &
             -2.*int5(ri5_79,e(:,OP_1),f(:,OP_DP),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
     endif

  case(1)
     temp79a = g(:,OP_DZ)*(-h(:,OP_DZ)*e(:,OP_DRZ) - h(:,OP_DR)*e(:,OP_DRR)  &
                          + e(:,OP_DZ)*h(:,OP_DRZ) + e(:,OP_DR)*h(:,OP_DRR)) &
             + g(:,OP_DR)*(-e(:,OP_DR)*h(:,OP_DRZ) - e(:,OP_DZ)*h(:,OP_DZZ)  &
                          + h(:,OP_DR)*e(:,OP_DRZ) + h(:,OP_DZ)*e(:,OP_DZZ))
     temp = int3(ri_79,f(:,OP_DP),temp79a,weight_79,79)

    
  end select

  v1vpsipsi = temp
#else
  v1vpsipsi = 0.
#endif
  return
end function v1vpsipsi



! V1vpsib
! =======
vectype function v1vpsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     select case(ivform)
     case(0)
        temp = 2.*(   int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
                  -   int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
                  +2.*int5(ri4_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1),weight_79,79))
#ifdef USECOMPLEX
        temp = temp - &
             (int5(ri4_79,e(:,OP_DZ),f(:,OP_DPP),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
             +int5(ri4_79,e(:,OP_DR),f(:,OP_DPP),g(:,OP_DR),h(:,OP_1),weight_79,79))
        if(itor.eq.1) then
           temp = temp - &
                2.*int5(ri5_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DR),h(:,OP_1),weight_79,79)
        endif
#endif

     case(1)
        if(itor.eq.1) then
        temp = 2.*(int5(ri_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
                  -int5(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79))
        endif
#ifdef USECOMPLEX
        temp = temp - &
             (int5(ri2_79,e(:,OP_DZ),f(:,OP_DPP),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
             +int5(ri2_79,e(:,OP_DR),f(:,OP_DPP),g(:,OP_DR),h(:,OP_1),weight_79,79))
#endif
     end select
  endif

  v1vpsib = temp
  return
end function v1vpsib



! V1chipsipsi
! ===========
vectype function v1chipsipsi(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  vectype :: temp

  temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
  temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)

  temp79e = e(:,OP_DRZ)*h(:,OP_DR ) - e(:,OP_DRR)*h(:,OP_DZ ) &
       +    e(:,OP_DZ )*h(:,OP_DRR) - e(:,OP_DR )*h(:,OP_DRZ)
  temp79f = e(:,OP_DZZ)*h(:,OP_DR ) - e(:,OP_DRZ)*h(:,OP_DZ ) &
       +    e(:,OP_DZ )*h(:,OP_DRZ) - e(:,OP_DR )*h(:,OP_DZZ)


  select case(ivform)
  case(0)

     temp = int3(ri_79,temp79b,temp79e,weight_79,79) &
          + int3(ri_79,temp79c,temp79f,weight_79,79) &
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

  case(1)
     temp = int3(ri3_79,temp79b,temp79e,weight_79,79) &
          + int3(ri3_79,temp79c,temp79f,weight_79,79) &
          + int4(ri3_79,e(:,OP_DR),temp79c,h(:,OP_GS),weight_79,79) &
          - int4(ri3_79,e(:,OP_DZ),temp79b,h(:,OP_GS),weight_79,79)
     
     if(itor.eq.1) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
        temp79d = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
        
        temp = temp + &
             (2.*int4(ri4_79,e(:,OP_DZ),temp79a,h(:,OP_GS),weight_79,79) &
             -   int3(ri4_79,temp79e,temp79a,weight_79,79) &
             +   int3(ri4_79,temp79d,temp79b,weight_79,79) &
             -2.*int3(ri5_79,temp79d,temp79a,weight_79,79))
     endif
  end select

  v1chipsipsi = temp
  return
end function v1chipsipsi


! V1chipsib
! =========
vectype function v1chipsib(e,f,g,h)
  
  use basic
  use arrays
  use nintegrate_mod

  implicit none
  
  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h
  vectype :: temp

   temp = 0.

#ifdef USECOMPLEX
   temp79a = f(:,OP_DRP )*e(:,OP_DZZ)-f(:,OP_DZP )*e(:,OP_DRZ) &
            -f(:,OP_DRRP)*e(:,OP_DZ )-f(:,OP_DRZP)*e(:,OP_DZ )
   temp79b = f(:,OP_DZP )*e(:,OP_DRR)-f(:,OP_DRP )*e(:,OP_DRZ) &
            -f(:,OP_DRZP)*e(:,OP_DR )-f(:,OP_DZZP)*e(:,OP_DZ )

   temp = temp +   &
        int4(ri2_79,g(:,OP_DR),temp79a,h(:,OP_1),weight_79,79) &
       +int4(ri2_79,g(:,OP_DZ),temp79b,h(:,OP_1),weight_79,79)

   if(itor.eq.1) then
      temp = temp - &
               int5(ri3_79,f(:,OP_DRP), e(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)  &
           -2.*int5(ri3_79,f(:,OP_DRRP),e(:,OP_1), g(:,OP_DR),h(:,OP_1),weight_79,79)  &
           -2.*int5(ri4_79,f(:,OP_DRP), e(:,OP_1), g(:,OP_DR),h(:,OP_1),weight_79,79)  &
           +   int5(ri3_79,f(:,OP_DZP), e(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79)  &
           -2.*int5(ri3_79,f(:,OP_DRP), e(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79)  &
           -2.*int5(ri3_79,f(:,OP_DRZP),e(:,OP_1), g(:,OP_DZ),h(:,OP_1),weight_79,79)
   endif
#endif

   v1chipsib = temp
   return
 end function v1chipsib



! V1chibb
! =======
vectype function v1chibb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  vectype :: temp

  select case(ivform)
  case(0)
     temp = 0.
     if(itor.eq.1) then
        temp = temp -2.* &
             (int5(ri2_79,e(:,OP_DZ),g(:,OP_1 ),f(:,OP_GS),h(:,OP_1),weight_79,79) &
             +int5(ri2_79,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_DZ),h(:,OP_1),weight_79,79) &
             +int5(ri2_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DR),h(:,OP_1),weight_79,79))
     endif
!
!-----------------start of code added 04/09/08 (scj)
#ifdef USECOMPLEX
     temp = temp +   &
          int5(ri3_79,f(:,OP_DZPP),e(:,OP_DR),g(:,OP_1),h(:,OP_1),weight_79,79) &
          -int5(ri3_79,f(:,OP_DRPP),e(:,OP_DZ),g(:,OP_1),h(:,OP_1),weight_79,79)
     
     if(itor.eq.1) then
        temp = temp + 2.*int5(ri4_79,f(:,OP_DZPP),e(:,OP_1),g(:,OP_1),h(:,OP_1),weight_79,79)
     endif
#endif
!------------------end of code added 04/09/08

  case(1)
     temp = 0.

  end select
  
  v1chibb = temp
  return
end function v1chibb


! V1chip
! ======
vectype function v1chip(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = 0.
     
  case(1)
     temp = 0.

     if(itor.eq.1) then
        temp = 2.* &
             (int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
             +int4(ri2_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
             +int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_1 ),weight_79,79))
     end if
  end select

  v1chip = temp
  return
end function v1chip


! V1psisb1
! ========
vectype function v1psisb1(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

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
vectype function v1bsb2(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

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
vectype function v1ngrav(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  temp = gravz*int3( r_79,e(:,OP_1),f(:,OP_DR),weight_79,79) &
       - gravr*int3(ri_79,e(:,OP_1),f(:,OP_DZ),weight_79,79)

  v1ngrav = temp
  return
end function v1ngrav


! V1ungrav
! ========
vectype function v1ungrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp79a = f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)

  temp = gravz*int2(       e(:,OP_DR),temp79a,weight_79,79) &
       - gravr*int3(ri2_79,e(:,OP_DZ),temp79a,weight_79,79)

  if(itor.eq.1) &
       temp = temp + 2.*gravz*int3(ri_79,e(:,OP_1),temp79a,weight_79,79)

  v1ungrav = temp
  return
end function v1ungrav


! V1chingrav
! ==========
vectype function v1chingrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp79a = r_79*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
       + g(:,OP_1)*f(:,OP_LP))
     
  temp = gravz*int2(       e(:,OP_DR),temp79a,weight_79,79) &
       - gravr*int3(ri2_79,e(:,OP_DZ),temp79a,weight_79,79)

  if(itor.eq.1) &
       temp = temp + 2.*gravz*int3(ri_79,e(:,OP_1),temp79a,weight_79,79)

  v1chingrav = temp
  return
end function v1chingrav


! V1ndenmgrav
! ===========
vectype function v1ndenmgrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  real, intent(in) :: g
  vectype :: temp

  temp79a = -g*r_79*f(:,OP_LP)

  temp = gravz*int2(       e(:,OP_DR),temp79a,weight_79,79) &
       - gravr*int3(ri2_79,e(:,OP_DZ),temp79a,weight_79,79)

  if(itor.eq.1) &
       temp = temp + 2.*gravz*int3(ri_79,e(:,OP_1),temp79a,weight_79,79)

  v1ndenmgrav = temp
  return
end function v1ndenmgrav


! V1us
! ====
vectype function v1us(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = -int3(e(:,OP_DZ),f(:,OP_DZ),temp79a,weight_79,79) &
            -int3(e(:,OP_DR),f(:,OP_DR),temp79a,weight_79,79)

     if(itor.eq.1) then 
        temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),weight_79,79)
     endif

  case(1)
     temp = -int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),temp79a,weight_79,79) &
            -int4(r2_79,e(:,OP_DR),f(:,OP_DR),temp79a,weight_79,79)

  end select

  v1us = temp
  return
end function v1us


! V1chis
! ======
vectype function v1chis(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)


  select case(ivform)
  case(0)
     temp = int4(r_79,e(:,OP_DZ),f(:,OP_DR),temp79a,weight_79,79) &
          - int4(r_79,e(:,OP_DR),f(:,OP_DZ),temp79a,weight_79,79)

     if(itor.eq.1) then 
        temp = temp - 2.*int3(e(:,OP_1),f(:,OP_DZ),temp79a,weight_79,79)
     endif

  case(1)
     temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),temp79a,weight_79,79) &
          - int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),temp79a,weight_79,79)
     
  end select

  v1chis = temp
  return
end function v1chis


! V1psif
! ======
vectype function v1psif(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
#ifdef USECOMPLEX
  temp = &
       - int3(e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZP),weight_79,79) &
       - int3(e(:,OP_DR),f(:,OP_GS),g(:,OP_DRP),weight_79,79)
  if(itor.eq.1) then
     temp = temp - 2.* & 
          int4(ri_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DRP),weight_79,79)
  endif
  
  v1psif = temp
#else
  v1psif = 0.
#endif

  return
end function v1psif


! V1bf
! ====
vectype function v1bf(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
#ifdef USECOMPLEX
  temp = &
       + int4(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DRPP),weight_79,79) &
       - int4(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZPP),weight_79,79)
  if(itor.eq.1) then
     temp = temp - 2.* & 
          int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZPP),weight_79,79)
  endif

  v1bf = temp
#else
  v1bf = 0.
#endif

  return
end function v1bf


! V1p
! ===
vectype function v1p(e,f)
  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f

  vectype :: temp
  temp = 0.

  select case(ivform)
  case(0)
    temp = 0.

  case(1)
    temp = &
       + int3(r_79,e(:,OP_DZ),f(:,OP_DR),weight_79,79) &
       - int3(r_79,e(:,OP_DR),f(:,OP_DZ),weight_79,79)

  end select

  v1p = temp
  return
end function v1p



!============================================================================
! V2 TERMS
!============================================================================


! V2vn
! ====
vectype function v2vn(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int3(e(:,OP_1),f(:,OP_1),g(:,OP_1),weight_79,79)
  case(1)
     temp = int4(r2_79,e(:,OP_1),f(:,OP_1),g(:,OP_1),weight_79,79)
  end select

  v2vn = temp
  return
end function v2vn


! V2vmu
! =====
vectype function v2vmu(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int3(e(:,OP_1),f(:,OP_GS),g(:,OP_1),weight_79,79) &
          + int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          + int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)
     
     if(itor.eq.1) then
        temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),weight_79,79)
     endif

#ifdef USECOMPLEX
     temp = temp + 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DPP),h(:,OP_1),weight_79,79)
#endif

     ! hyperviscous
     if(hypv.ne.0.) then
        temp79a = e(:,OP_GS)*g(:,OP_1) + &
             e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR)
        if(itor.eq.1) temp79a = temp79a + 4.*ri_79*e(:,OP_DR)*g(:,OP_1)
        temp = temp - int3(temp79a,f(:,OP_GS),i(:,OP_1),weight_79,79)  
     endif

  case(1)

     temp = -int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          -  int4(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79)
     
     ! hyperviscous
     if(hypv.ne.0.) then
        temp79a = e(:,OP_GS)*g(:,OP_1) + &
             e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR)
        if(itor.eq.1) temp79a = temp79a + 4.*ri_79*e(:,OP_DR)*g(:,OP_1)
        
        temp = temp - int4(r2_79,temp79a,f(:,OP_GS),i(:,OP_1),weight_79,79)
        if(itor.eq.1) then
           temp = temp - 4.*int4(r_79,temp79a,f(:,OP_DR),i(:,OP_1),weight_79,79)
        endif
     end if

#ifdef USECOMPLEX
     temp = temp + 2.*int3(e(:,OP_1),f(:,OP_DPP),h(:,OP_1),weight_79,79)
#endif

  end select

  v2vmu = temp
  return
end function v2vmu


! V2vun
! =====
vectype function v2vun(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
  case(1)
     temp = int5(r3_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(r3_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp = temp + &
             2.*int5(r2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     end if

  end select

  v2vun = temp
  return
end function v2vun


! V2up
! ====
vectype function v2up(e,f,g)
  use basic
  use arrays
  use nintegrate_mod
  
  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.

#ifdef USECOMPLEX
  select case(ivform)
    case(0)
      temp = -int4(ri_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR),weight_79,79) &
           +  int4(ri_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),weight_79,79)
    case(1)
      temp = ( -int4(r_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR),weight_79,79) &
           +    int4(r_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),weight_79,79))/rzero**2
    end select
#endif

  v2up = temp
  return
end function v2up

! V2vp
! =======
vectype function v2vp(e,f,g)
  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.

#ifdef USECOMPLEX
  temp = gam*int3(e(:,OP_1),f(:,OP_DPP),g(:,OP_1),weight_79,79)
#endif

  v2vp = temp
  return
end function v2vp


! V2chip
! =======
vectype function v2chip(e,f,g)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.

#ifdef USECOMPLEX
  temp =     int3(e(:,OP_1),f(:,OP_DRP),g(:,OP_DR),weight_79,79)    &
           + int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ),weight_79,79)    &
       + gam*int3(e(:,OP_1),f(:,OP_DRRP),g(:,OP_1),weight_79,79)    &
       + gam*int3(e(:,OP_1),f(:,OP_DZZP),g(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp +   &
          + gam*int4(ri_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_1),weight_79,79)
  end if
#endif

  v2chip = temp
  return
end function v2chip

! V2p
! ===
vectype function v2p(e,f)
  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f

  vectype :: temp
  temp = 0.

#ifdef USECOMPLEX
  temp = -int2(e(:,OP_1),f(:,OP_DP),weight_79,79)
#endif

  v2p = temp
  return
end function v2p


! V2psipsi
! ========
vectype function v2psipsi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  temp = 0.5* &
       (int4(ri2_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ),weight_79,79) &
       +int4(ri2_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DR),weight_79,79))
  v2psipsi = -temp
#else
  v2psipsi = 0.
#endif
  return
end function v2psipsi


! V2psib
! ======
vectype function v2psib(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  v2psib = temp
  return
end function v2psib


! V2vpsipsi
! =========
vectype function v2vpsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  ! [nu,psi(2)]
  temp79a = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)

  select case(ivform)
  case(0)
     temp = int4(ri2_79,f(:,OP_DR),g(:,OP_DZ),temp79a,weight_79, 79) &
          - int4(ri2_79,f(:,OP_DZ),g(:,OP_DR),temp79a,weight_79, 79)
     if(itor.eq.1) then
        temp = temp - 2.*int4(ri3_79,f(:,OP_1),g(:,OP_DZ),temp79a,weight_79,79)
     endif

#ifdef USECOMPLEX
     temp = temp + &
          (int5(ri4_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          +int5(ri4_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DR),h(:,OP_DR),weight_79,79))
#endif

  case(1)
     temp = int3(f(:,OP_DR),g(:,OP_DZ),temp79a,weight_79, 79) &
          - int3(f(:,OP_DZ),g(:,OP_DR),temp79a,weight_79, 79)

#ifdef USECOMPLEX
     temp = temp + &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DR),h(:,OP_DR),weight_79,79))
#endif

  end select

  v2vpsipsi = temp
  return
end function v2vpsipsi


! V2vpsib
! =======
vectype function v2vpsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  select case(ivform)
  case(0)
     temp = int5(ri3_79,e(:,OP_1),f(:,OP_DP),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          - int5(ri3_79,e(:,OP_1),f(:,OP_DP),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
  case(1)
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DP),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DP),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
  end select
  v2vpsib = -temp
#else
  v2vpsib = 0.
#endif
  return
end function v2vpsib



! V2upsipsi
! =========
vectype function v2upsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX  
  select case(ivform)
  case(0)
     temp79a = h(:,OP_DZ)*(f(:,OP_DZZP)*g(:,OP_DR )-f(:,OP_DRZP)*g(:,OP_DZ )  &
                       +f(:,OP_DZP )*g(:,OP_DRZ)-f(:,OP_DRP )*g(:,OP_DZZ)) &
            + h(:,OP_DR)*(f(:,OP_DRZP)*g(:,OP_DR )-f(:,OP_DRRP)*g(:,OP_DZ )  &
                       +f(:,OP_DZP )*g(:,OP_DRR)-f(:,OP_DRP )*g(:,OP_DRZ))
    temp = -int3(ri3_79,e(:,OP_1),temp79a,weight_79,79)

    if(itor.eq.1) then
       temp79a = (f(:,OP_DZP)*g(:,OP_DR)-f(:,OP_DRP)*g(:,OP_DZ))*h(:,OP_DR)
       temp = temp + int3(ri4_79,e(:,OP_1),temp79a,weight_79,79)
    endif
 
  case(1)
   temp79a = h(:,OP_DZ)*(f(:,OP_DZZP)*g(:,OP_DR )-f(:,OP_DRZP)*g(:,OP_DZ )  &
                       +f(:,OP_DZP )*g(:,OP_DRZ)-f(:,OP_DRP )*g(:,OP_DZZ)) &
            + h(:,OP_DR)*(f(:,OP_DRZP)*g(:,OP_DR )-f(:,OP_DRRP)*g(:,OP_DZ )  &
                       +f(:,OP_DZP )*g(:,OP_DRR)-f(:,OP_DRP )*g(:,OP_DRZ))
    temp = -int3(ri_79,e(:,OP_1),temp79a,weight_79,79)/rzero**2

    if(itor.eq.1) then
       temp79a = (f(:,OP_DZP)*g(:,OP_DR)-f(:,OP_DRP)*g(:,OP_DZ))*h(:,OP_DR)
       temp = temp - int3(ri2_79,e(:,OP_1),temp79a,weight_79,79)/rzero**2
    endif

  end select
#else
  temp = 0.
#endif

  v2upsipsi = temp
  return
end function v2upsipsi



! V2upsib
! =======
vectype function v2upsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp79a = e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ)
  select case(ivform)
  case(0)

    temp = int4(ri2_79,temp79a,g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
         - int4(ri2_79,temp79a,g(:,OP_DZ),h(:,OP_DR),weight_79,79)
    if(itor.eq.1) then
       temp = temp + 2.* &
            (int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
            -int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79))
    endif

#ifdef USECOMPLEX
    temp = temp - &
         (int5(ri4_79,e(:,OP_1),f(:,OP_DZPP),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
         +int5(ri4_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_DR),h(:,OP_1),weight_79,79))
#endif

  case(1)
     temp = (int3(temp79a,g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          -  int3(temp79a,g(:,OP_DZ),h(:,OP_DR),weight_79,79))/rzero**2

#ifdef USECOMPLEX
     temp = temp - &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZPP),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_DR),h(:,OP_1),weight_79,79))/rzero**2
#endif

  end select

  v2upsib = temp
  return
end function v2upsib


! V2ubb
! =====
vectype function v2ubb(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  select case(ivform)
  case(0)
    temp = &
         (int5(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
         -int5(ri3_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR),h(:,OP_1),weight_79,79))
  case(1)
   temp = &
         (int5(ri_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
         -int5(ri_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR),h(:,OP_1),weight_79,79))/rzero**2
  end select

#else
  temp = 0.
#endif
  v2ubb = temp
  return
end function v2ubb



! v2upsisb2
! ========
vectype function v2upsisb2(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  v2upsisb2 = 0.
  return
end function v2upsisb2


! v2ubsb1
! =======
vectype function v2ubsb1(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  v2ubsb1 = 0.
  return
end function v2ubsb1


! v2chipsipsi
! =========
vectype function v2chipsipsi(e,f,g,h)
  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#ifdef USECOMPLEX
  temp = temp +   &
       +int5(ri2_79,f(:,OP_DRRP),e(:,OP_1),g(:,OP_DR ),h(:,OP_DR),weight_79,79) &
       +int5(ri2_79,f(:,OP_DRP ),e(:,OP_1),g(:,OP_DRR),h(:,OP_DR),weight_79,79) &
       +int5(ri2_79,f(:,OP_DRZP),e(:,OP_1),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
       +int5(ri2_79,f(:,OP_DZP ),e(:,OP_1),g(:,OP_DRZ),h(:,OP_DR),weight_79,79) &
       +int5(ri2_79,f(:,OP_DRZP),e(:,OP_1),g(:,OP_DR ),h(:,OP_DZ),weight_79,79) &
       +int5(ri2_79,f(:,OP_DRP ),e(:,OP_1),g(:,OP_DRZ),h(:,OP_DZ),weight_79,79) &
       +int5(ri2_79,f(:,OP_DZZP),e(:,OP_1),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79) &
       +int5(ri2_79,f(:,OP_DZP ),e(:,OP_1),g(:,OP_DZZ),h(:,OP_DZ),weight_79,79)
#endif

  v2chipsipsi = temp
  return
end function v2chipsipsi


! v2chipsib
! =========
vectype function v2chipsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp


  temp79a = h(:,OP_1 )*f(:,OP_GS) &
       +    h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)

  temp79b = e(:,OP_DR)*h(:,OP_DZ) - e(:,OP_DZ)*h(:,OP_DR)

  temp79c = e(:,OP_DZ)*g(:,OP_DR) - e(:,OP_DR)*g(:,OP_DZ)
  
  select case(ivform)
  case(0)
     temp = int3(ri_79,temp79a,temp79c,weight_79,79) &
          + int4(ri_79,temp79b,f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          + int4(ri_79,temp79b,f(:,OP_DR),g(:,OP_DR),weight_79,79) 
#ifdef USECOMPLEX
     temp = temp +   &
          -int5(ri3_79,f(:,OP_DZPP),e(:,OP_1),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          +int5(ri3_79,f(:,OP_DRPP),e(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
#endif
     
  case(1)

     temp = int3(ri3_79,temp79a,temp79c,weight_79,79) &
          + int4(ri3_79,temp79b,f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          + int4(ri3_79,temp79b,f(:,OP_DR),g(:,OP_DR),weight_79,79)

     if(itor.eq.1) then
        temp = temp - &
             2.*int4(ri4_79,temp79c,f(:,OP_DR),h(:,OP_1),weight_79,79)
     endif
  end select

  v2chipsib = temp
  return
end function v2chipsib


! v2chibb
! =========
vectype function v2chibb(e,f,g,h)
  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#ifdef USECOMPLEX
  temp = temp +   &
      +int5(ri2_79,f(:,OP_DRP),e(:,OP_1),g(:,OP_DR),h(:,OP_1),weight_79,79) &
      +int5(ri2_79,f(:,OP_DZP),e(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
#endif

  v2chibb = temp
  return
end function v2chibb



! v2vchin
! =======
vectype function v2vchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     temp =-int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)
  case(1)
     temp =-int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp = temp - 2.*int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_1),weight_79,79)
     endif
  end select

  v2vchin = temp
  return
end function v2vchin


! v2chibsb1
! =========
vectype function v2chibsb1(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  v2chibsb1 = 0.
  return
end function v2chibsb1

! v2psisb2
! ========
vectype function v2psisb2(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  v2psisb2 = temp
  return
end function v2psisb2

! v2bsb1
! ======
vectype function v2bsb1(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79)

  v2bsb1 = temp
  return
end function v2bsb1


! V2vs
! ====
vectype function v2vs(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v2vs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = -int3(e(:,OP_1),f(:,OP_1),temp79a,weight_79,79)
  case(1)
     temp = -int4(r2_79,e(:,OP_1),f(:,OP_1),temp79a,weight_79,79)
  end select

  v2vs = temp
  return
end function v2vs


! V2psif
! ======
vectype function v2psif(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
#ifdef USECOMPLEX
  temp = &
       + int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZPP),weight_79,79) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DRPP),weight_79,79)
  
  v2psif = temp
#else
  v2psif = 0.
#endif

  return
end function v2psif


! V2bf
! ====
vectype function v2bf(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
#ifdef USECOMPLEX
  temp = &
       - int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZP),weight_79,79) &
       - int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DRP),weight_79,79)

  v2bf = temp
#else
  v2bf = 0.
#endif

  return
end function v2bf



!===============================================================================
! V3 TERMS
!===============================================================================

! V3chin
! ======
vectype function v3chin(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  
  vectype :: temp

  select case(ivform)
  case(0)
     temp = - int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
            - int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)

  case(1)
     temp = - int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
            - int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  end select

  v3chin = temp
  return
end function v3chin



! V3chimu
! =======
vectype function v3chimu(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     temp = 2.*int3(e(:,OP_LP),f(:,OP_LP),h(:,OP_1 ),weight_79,79) &
          +    int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          +    int3(e(:,OP_LP),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          +    int3(e(:,OP_DZ),f(:,OP_LP),g(:,OP_DZ),weight_79,79) &
          +    int3(e(:,OP_DR),f(:,OP_LP),g(:,OP_DR),weight_79,79) &
          +    int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_LP),weight_79,79) &
          +    int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_LP),weight_79,79)
  case(1)
     temp79a = e(:,OP_DZZ)*f(:,OP_DZZ) + e(:,OP_DRR)*f(:,OP_DRR) &
          + 2.*e(:,OP_DRZ)*f(:,OP_DRZ)

     if(itor.eq.1) then
        temp79a = temp79a - 2.*ri_79* &
             (e(:,OP_DR)*f(:,OP_DRR) + e(:,OP_DZ)*f(:,OP_DRZ) &
             +f(:,OP_DR)*e(:,OP_DRR) + f(:,OP_DZ)*e(:,OP_DRZ)) &
             + 5.*ri2_79*e(:,OP_DR)*f(:,OP_DR) &
             + 2.*ri2_79*e(:,OP_DZ)*f(:,OP_DZ)
     end if

     temp = 2.*int3(ri4_79,temp79a,g(:,OP_1),weight_79,79) &
          + 2.*int4(ri4_79,e(:,OP_GS),f(:,OP_GS),h(:,OP_1),weight_79,79) &
          - 2.*int4(ri4_79,e(:,OP_GS),f(:,OP_GS),g(:,OP_1),weight_79,79)
  end select

  v3chimu = temp
  return
end function v3chimu


! V3umu
! =====
vectype function v3umu(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(ri_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          + int4(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),weight_79,79) &
          - int4(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
          + int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_LP),weight_79,79) &
          - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_LP),weight_79,79)

  case(1)
     temp79a = f(:,OP_DRZ)*(e(:,OP_DZZ)-e(:,OP_DRR)) &
          -    e(:,OP_DRZ)*(f(:,OP_DZZ)-f(:,OP_DRR))
     temp = 2.*int3(ri_79,temp79a,g(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp79a = e(:,OP_DZ)*(f(:,OP_DZZ)-f(:,OP_DRR)) &
             + 2.*e(:,OP_DR)*f(:,OP_DRZ) &
             - e(:,OP_DRR)*f(:,OP_DZ) + e(:,OP_DRZ)*f(:,OP_DR) &
             + ri_79*(e(:,OP_DR)*f(:,OP_DZ) - e(:,OP_DZ)*f(:,OP_DR))
        
        temp79b = h(:,OP_1) - g(:,OP_1)
        
        temp = temp &
             + 2.*int3(ri2_79,temp79a,g(:,OP_1),weight_79,79) &
             - 4.*int4(ri2_79,e(:,OP_GS),f(:,OP_DZ),temp79b,weight_79,79)
     endif
  end select

  v3umu = temp
  return
end function v3umu


! V3un
! ====
vectype function v3un(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) 
  case(1)
     temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),weight_79,79) 
  end select

  v3un = temp
  return
end function v3un


! V3p
! ===
vectype function v3p(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  select case(ivform)
  case(0)
     temp = -int2(e(:,OP_1),f(:,OP_LP),weight_79,79)
     
  case(1)
     temp = int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
          + int3(ri2_79,e(:,OP_DR),f(:,OP_DR),weight_79,79)
  end select

  v3p = temp
  return
end function v3p



! V3up
! ====
vectype function v3up(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(ri_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  case(1)
     temp = int4(ri_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp - 2.*gam* &
             int4(ri2_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_1),weight_79,79)
     endif
  end select

  v3up = temp
  
  return
end function v3up


! V3vp
! =======
vectype function v3vp(e,f,g)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.

#ifdef USECOMPLEX
  select case(ivform)
  case(0)
  temp = gam*int4(ri2_79,e(:,OP_LP),f(:,OP_DP),g(:,OP_1),weight_79,79)

  case(1)
  temp = gam*int3(e(:,OP_LP),f(:,OP_DP),g(:,OP_1),weight_79,79)

  end select
#endif

  v3vp = temp
  return
end function v3vp



! V3chip
! ======
vectype function v3chip(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          + int3(e(:,OP_LP),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          + gam*int3(e(:,OP_LP),f(:,OP_LP),g(:,OP_1),weight_79,79)

  case(1)
     temp = int3(e(:,OP_GS),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          + int3(e(:,OP_GS),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          + gam*int3(e(:,OP_GS),f(:,OP_GS),g(:,OP_1),weight_79,79)
  end select

  v3chip = temp

  return
end function v3chip



! V3psipsi
! ========
vectype function v3psipsi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS),weight_79,79) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS),weight_79,79)

  case(1)
     temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS),weight_79,79) &
          + int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS),weight_79,79)

  end select

  v3psipsi = temp

  return
end function v3psipsi


! V3bb
! ====
vectype function v3bb(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     temp = - &
          (int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          +int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          +int4(ri2_79,e(:,OP_1),f(:,OP_1 ),g(:,OP_GS),weight_79,79))
  case(1)
     temp = int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
          + int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  end select

  v3bb = temp
  return
end function v3bb


! V3psisb1
! ========
vectype function v3psisb1(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),weight_79,79) &
       + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS),weight_79,79)

  v3psisb1 = temp

end function v3psisb1


! V3bsb2
! ======
vectype function v3bsb2(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1 ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_1 ),g(:,OP_DR),weight_79,79)

  v3bsb2 = temp
  return
end function v3bsb2


! V3upsipsi
! =========
vectype function v3upsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  vectype :: temp

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
  
     select case(ivform)
     case(0)
        temp = temp &
             - int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_GS ),weight_79,79) &
             + int4(ri4_79,e(:,OP_DRZ),temp79a,h(:,OP_DZ ),weight_79,79) &
             + int4(ri4_79,e(:,OP_DZ ),temp79a,h(:,OP_DRZ),weight_79,79) &
             + int4(ri4_79,e(:,OP_DRR),temp79a,h(:,OP_DR ),weight_79,79) &
             + int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_DRR),weight_79,79)

     case(1)

        temp79d = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)
        
        temp = temp &
             -2.*int3(ri4_79,temp79b,temp79d,weight_79,79) &
             -   int4(ri4_79,e(:,OP_DRZ),temp79a,h(:,OP_DZ ),weight_79,79) &
             -   int4(ri4_79,e(:,OP_DZ ),temp79a,h(:,OP_DRZ),weight_79,79) &
             -   int4(ri4_79,e(:,OP_DRR),temp79a,h(:,OP_DR ),weight_79,79) &
             -   int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_DRR),weight_79,79) &
             +   int4(ri4_79,e(:,OP_DR),temp79a,h(:,OP_GS),weight_79,79) &
             +2.*int3(ri5_79,temp79d,temp79a,weight_79,79)
     end select
  endif

  v3upsipsi = temp
  return
end function v3upsipsi


! V3upsib
! =======
vectype function v3upsib(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none
  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  vectype :: temp
  temp = 0.

#ifdef USECOMPLEX
  temp = temp    &
       +int5(ri4_79,f(:,OP_DZZP),e(:,OP_DR ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
       -int5(ri4_79,f(:,OP_DRZP),e(:,OP_DZ ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
       -int5(ri4_79,f(:,OP_DRP ),e(:,OP_DRR),g(:,OP_DR),h(:,OP_1),weight_79,79) &
       -int5(ri4_79,f(:,OP_DZP ),e(:,OP_DRZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
       +int5(ri4_79,f(:,OP_DRRP),e(:,OP_DZ ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       -int5(ri4_79,f(:,OP_DRZP),e(:,OP_DR ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       -int5(ri4_79,f(:,OP_DZP ),e(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       -int5(ri4_79,f(:,OP_DRP ),e(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1),weight_79,79)
 
  if(itor.eq.1) then
     temp = temp    &
          -int5(ri5_79,f(:,OP_DRP),e(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          -int5(ri5_79,f(:,OP_DRP),e(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79)
  endif
#endif
  v3upsib = temp
  return
end function v3upsib



! V3ubb
! =====
vectype function v3ubb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri3_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp = temp + 2.* &
             int5(ri4_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_1),h(:,OP_1),weight_79,79)
     endif
#ifdef USECOMPLEX
     temp = temp -   &
          int5(ri5_79,f(:,OP_DRPP),e(:,OP_DZ),g(:,1),h(:,OP_1),weight_79,79) &
          +int5(ri5_79,f(:,OP_DZPP),e(:,OP_DR),g(:,1),h(:,OP_1),weight_79,79)
#endif

  case(1)
     temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

     temp = int3(ri3_79,e(:,OP_GS),temp79a,weight_79,79)

     if(itor.eq.1) then
        temp = temp - &
             2.*int3(ri4_79,e(:,OP_DR),temp79a,weight_79,79)
     endif


  end select


  v3ubb = temp
  return
end function v3ubb


! v3vpsipsi
! =========
vectype function v3vpsipsi(e,f,g,h)
!
!  e trial
!  f lin
!  g psi
!  h psi
  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#ifdef USECOMPLEX

  select case(ivform)
  case(0)
  temp = temp +   &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_DRR),h(:,OP_DR),weight_79,79) &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DRR),g(:,OP_DR ),h(:,OP_DR),weight_79,79) &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_DR),weight_79,79) &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_DRZ),h(:,OP_DZ),weight_79,79) &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DRZ),g(:,OP_DR ),h(:,OP_DZ),weight_79,79) &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_DZZ),h(:,OP_DZ),weight_79,79) &
      +int5(ri4_79,f(:,OP_DP),e(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79) &
      -int5(ri4_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_GS ),h(:,OP_DR),weight_79,79)  &
      -int5(ri4_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DZ),weight_79,79)
  case(1)
  temp = temp +   &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_DRR),h(:,OP_DR),weight_79,79) &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DRR),g(:,OP_DR ),h(:,OP_DR),weight_79,79) &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_DR),weight_79,79) &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_DRZ),h(:,OP_DZ),weight_79,79) &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DRZ),g(:,OP_DR ),h(:,OP_DZ),weight_79,79) &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_DZZ),h(:,OP_DZ),weight_79,79) &
      +int5(ri2_79,f(:,OP_DP),e(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79) &
      -int5(ri2_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_GS ),h(:,OP_DR),weight_79,79)  &
      -int5(ri2_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DZ),weight_79,79)
  end select
#endif

  v3vpsipsi = temp
  return
end function v3vpsipsi



! v3vpsib
! =======
vectype function v3vpsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri3_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp = temp - 2.* &
             int5(ri4_79,e(:,OP_GS),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     endif

#ifdef USECOMPLEX
     temp = temp    &
          +int5(ri5_79,f(:,OP_DPP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          -int5(ri5_79,f(:,OP_DPP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) 
#endif
     
  case(1)
     temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

     temp = int3(ri3_79,e(:,OP_GS),temp79a,weight_79,79)

     if(itor.eq.1) then
        temp = temp - &
             2.*int3(ri4_79,e(:,OP_DR),temp79a,weight_79,79)
     endif

  end select

  v3vpsib = temp
  return
end function v3vpsib


! V3chipsipsi
! ===========
vectype function v3chipsipsi(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h
  vectype :: temp

  ! <f,g>,r
  temp79b = f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DRZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRR) + f(:,OP_DZ )*g(:,OP_DRZ)
  
  ! <f,g>,z
  temp79c = f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRZ) + f(:,OP_DZ )*g(:,OP_DZZ)

  ! <e,h>,r
  temp79e = e(:,OP_DRR)*h(:,OP_DR ) + e(:,OP_DRZ)*h(:,OP_DZ ) &
       +    e(:,OP_DR )*h(:,OP_DRR) + e(:,OP_DZ )*h(:,OP_DRZ)
  
  ! <e,h>,z
  temp79f = e(:,OP_DRZ)*h(:,OP_DR ) + e(:,OP_DZZ)*h(:,OP_DZ ) &
       +    e(:,OP_DR )*h(:,OP_DRZ) + e(:,OP_DZ )*h(:,OP_DZZ)
  
  select case(ivform)
  case(0)

     temp = int3(ri2_79,temp79b,temp79e,weight_79,79) &
          + int3(ri2_79,temp79c,temp79f,weight_79,79) &
          - int4(ri2_79,e(:,OP_DZ),temp79c,h(:,OP_GS),weight_79,79) &
          - int4(ri2_79,e(:,OP_DR),temp79b,h(:,OP_GS),weight_79,79)

  case(1)
     temp = int3(ri6_79,temp79b,temp79e,weight_79,79) &
          + int3(ri6_79,temp79c,temp79f,weight_79,79) &
          - int4(ri6_79,e(:,OP_DZ),temp79c,h(:,OP_GS),weight_79,79) &
          - int4(ri6_79,e(:,OP_DR),temp79b,h(:,OP_GS),weight_79,79)

     if(itor.eq.1) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
        temp79d = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)

        temp = temp &
             -2.*int3(ri7_79,temp79b,temp79d,weight_79,79) &
             -2.*int3(ri7_79,temp79a,temp79e,weight_79,79) &
             +2.*int4(ri7_79,e(:,OP_DR),temp79a,h(:,OP_GS),weight_79,79) &
             +4.*int3(ri8_79,temp79a,temp79d,weight_79,79)
     endif
     
  end select

  v3chipsipsi = temp
  return
end function v3chipsipsi


! V3chipsib
! =====
vectype function v3chipsib(e,f,g,h)
  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  vectype :: temp
  temp = 0.

#ifdef USECOMPLEX
  temp = temp    &
      +int5(ri3_79,f(:,OP_DRRP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
      +int5(ri3_79,f(:,OP_DRP),e(:,OP_DRZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
      -int5(ri3_79,f(:,OP_DRZP),e(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79) &
      -int5(ri3_79,f(:,OP_DZP),e(:,OP_DRR),g(:,OP_DR),h(:,OP_1),weight_79,79) &
      +int5(ri3_79,f(:,OP_DRZP),e(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
      +int5(ri3_79,f(:,OP_DRP),e(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
      -int5(ri3_79,f(:,OP_DZZP),e(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
      -int5(ri3_79,f(:,OP_DZP),e(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1),weight_79,79)

  if(itor.eq.1) then
  temp = temp    &
       -int5(ri4_79,f(:,OP_DZP),e(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79) &
       +int5(ri4_79,f(:,OP_DRP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
  endif
#endif
  v3chipsib = temp
  return
end function v3chipsib


! V3chibb
! =======
vectype function v3chibb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int5(ri2_79,e(:,OP_GS),f(:,OP_GS),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
          + int5(ri2_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          + int5(ri2_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)

#ifdef USECOMPLEX
     temp = temp    &
          - int5(ri4_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1),weight_79,79) &
          - int5(ri4_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),h(:,OP_1),weight_79,79)
#endif
     
  case(1)
     temp79a = h(:,OP_1)*(g(:,OP_1)*f(:,OP_GS) &
          + g(:,OP_DZ)*f(:,OP_DZ) + g(:,OP_DR)*f(:,OP_DR))

     if(itor.eq.1) then
        temp79a = temp79a - &
             2.*ri_79*f(:,OP_DR)*g(:,OP_1)*h(:,OP_1)
     end if

     temp = int3(ri6_79,e(:,OP_GS),temp79a,weight_79,79)

     if(itor.eq.1) then
        temp = temp - &
             2.*int3(ri7_79,e(:,OP_DR),temp79a,weight_79,79)
     endif

  end select

  v3chibb = temp

  return
end function v3chibb


! V3uun
! =====
vectype function v3uun(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  vectype :: temp

  select case(ivform)
  case(0)
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

  case(1)
     temp79a = e(:,OP_DZ)*(f(:,OP_DR)*g(:,OP_DRZ) - f(:,OP_DZ)*g(:,OP_DRR)) &
          +    e(:,OP_DR)*(f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ))

     temp = int2(temp79a,h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp79a = g(:,OP_DZ)*(e(:,OP_DR)*f(:,OP_DZ) - e(:,OP_DZ)*f(:,OP_DR))
        temp = temp + int3(ri_79,temp79a,h(:,OP_1),weight_79,79)
     endif
  end select
     

  v3uun = temp
  return
end function v3uun


! V3vvn
! =====
vectype function v3vvn(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h
  vectype :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     select case(ivform)
     case(0)
        temp = -int5(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_1),h(:,OP_1),weight_79,79)
     case(1)
        temp = -int5(ri_79,e(:,OP_DR),f(:,OP_1),g(:,OP_1),h(:,OP_1),weight_79,79)
     end select
  endif

  v3vvn = temp
  return
end function v3vvn


! V3uchin
! =======
vectype function v3uchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  vectype :: temp

  temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

  select case(ivform)
  case(0)
     temp = int4(ri_79,e(:,OP_LP),temp79a,h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79) & 
          - int5(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) & 
          + int4(ri_79,e(:,OP_DZ),temp79a,h(:,OP_DZ),weight_79,79) & 
          + int4(ri_79,e(:,OP_DR),temp79a,h(:,OP_DR),weight_79,79)
  case(1)
     temp79a = e(:,OP_DZ)*(f(:,OP_DR)*g(:,OP_DZZ) + f(:,OP_DRZ)*g(:,OP_DZ)  &
                          -f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DRR)*g(:,OP_DR)) &
              +e(:,OP_DR)*(f(:,OP_DZ)*g(:,OP_DRR) + f(:,OP_DRZ)*g(:,OP_DR)  &
                          -f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZZ)*g(:,OP_DZ))
     temp = int3(ri3_79,temp79a,h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp79a = e(:,OP_DZ)*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)) &
             +    f(:,OP_DZ)*(e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))

        temp = temp + int3(ri4_79,temp79a,h(:,OP_1),weight_79,79)
     end if
  end select
  
  v3uchin = temp
  return
end function v3uchin



! V3chichin
! =========
vectype function v3chichin(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_1),weight_79,79)

  case(1)
     temp79a = e(:,OP_DZ)*(f(:,OP_DZ)*g(:,OP_DZZ) + f(:,OP_DR)*g(:,OP_DRZ)) &
          +    e(:,OP_DR)*(f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DR)*g(:,OP_DRR))

     temp = int3(ri6_79,temp79a,h(:,OP_1),weight_79,79)

     if(itor.eq.1) then
        temp79a = g(:,OP_DR)*(e(:,OP_DZ)*f(:,OP_DZ) + e(:,OP_DR)*f(:,OP_DR))

        temp = temp - 2.*int3(ri7_79,temp79a,h(:,OP_1),weight_79,79)
     endif

  end select

  v3chichin = temp
  return
end function v3chichin


! V3ngrav
! =======
vectype function v3ngrav(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  temp = gravz*int2(       e(:,OP_DZ),f(:,OP_1),weight_79,79) & 
       + gravr*int3(ri2_79,e(:,OP_DR),f(:,OP_1),weight_79,79) 

  v3ngrav = temp
  return
end function v3ngrav


! V3ungrav
! ========
vectype function v3ungrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
  
  temp = gravz*int3( ri_79,e(:,OP_DZ),temp79a,weight_79,79) &
       + gravr*int3(ri3_79,e(:,OP_DR),temp79a,weight_79,79)

  v3ungrav = temp
  return
end function v3ungrav


! V3chingrav
! ==========
vectype function v3chingrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  temp79a = -(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
       + f(:,OP_LP)*g(:,OP_1))

  temp = gravz*int2(       e(:,OP_DZ),temp79a,weight_79,79) &
       + gravr*int3(ri2_79,e(:,OP_DR),temp79a,weight_79,79)

  v3chingrav = temp
  return
end function v3chingrav


! V3ndenmgrav
! ===========
vectype function v3ndenmgrav(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  real, intent(in) :: g
  vectype :: temp

  temp = gravz*int2(       e(:,OP_DZ),f(:,OP_LP),weight_79,79) &
       + gravr*int3(ri2_79,e(:,OP_DR),f(:,OP_LP),weight_79,79)

  v3ndenmgrav = g*temp
  return
end function v3ndenmgrav


! V3us
! ====
vectype function v3us(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DR),temp79a,weight_79,79) &
        -int4(ri_79,e(:,OP_DR),f(:,OP_DZ),temp79a,weight_79,79)

  v3us = temp
  return
end function v3us


! V3chis
! ======
vectype function v3chis(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0) 
     temp = int3(e(:,OP_DZ),f(:,OP_DZ),temp79a,weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DR),temp79a,weight_79,79)

  case(1)
     temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),temp79a,weight_79,79) &
          + int4(ri4_79,e(:,OP_DR),f(:,OP_DR),temp79a,weight_79,79)
  end select

  v3chis = temp
  return
end function v3chis


!===============================================================================
! B1 TERMS
!===============================================================================


! B1psi
! =====
vectype function b1psi(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  if(jadv.eq.0) then
     temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)
  else
     temp79a = e(:,OP_DZ)*f(:,OP_DZ) + e(:,OP_DR)*f(:,OP_DR)
     temp = -int2(ri2_79,temp79a,weight_79,79)
  endif

  b1psi = temp
  return
end function b1psi


! B1psiu
! ======
vectype function b1psiu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
  select case(ivform)
  case(0)
     if(jadv.eq.0) then
        temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
             - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     else
        temp = int4(ri3_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
             + int4(ri3_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRR),weight_79,79) &
             - int4(ri3_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DZ ),weight_79,79) &
             - int4(ri3_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
             + int4(ri3_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DR ),weight_79,79) &
             + int4(ri3_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
             - int4(ri3_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
             - int4(ri3_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DZZ),weight_79,79) 
        if(itor.eq.1) then
           temp = temp &
                - int4(ri4_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
                + int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),weight_79,79) 
        endif
     endif
  case(1)
     if(jadv.eq.0) then
!!$     CHANGED NMF 4/30/08
        temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
             - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     else
        temp = (int4(ri_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
             + int4(ri_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRR),weight_79,79) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DZ ),weight_79,79) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
             + int4(ri_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DR ),weight_79,79) &
             + int4(ri_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
             - int4(ri_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
             - int4(ri_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DZZ),weight_79,79))/rzero**2 
        if(itor.eq.1) then
           temp = temp &
                + (int4(ri2_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
                - int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),weight_79,79))/rzero**2 
        endif
     endif
  end select
  
  b1psiu = temp
  return
end function b1psiu


! B1psiv
! ======
vectype function b1psiv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  if(jadv.eq.0) then
    temp = 0.
  else	
    select case(ivform)
    case(0)
       temp = int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DP),weight_79,79) &
            + int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DP),weight_79,79)
    case(1)
       temp = int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DP),weight_79,79) &
            + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DP),weight_79,79)
    end select
  endif
#else
  temp = 0
#endif

  b1psiv = temp
  return
end function b1psiv


! B1bu
! ====
vectype function b1bu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  if(jadv.eq.0) then
     temp = 0.
  else
     select case(ivform)
     case(0)    
     temp = -(int4(ri4_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZP),weight_79,79) &
            + int4(ri4_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DRP),weight_79,79)) 
     case(1)
     temp = -(int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZP),weight_79,79) &
            + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DRP),weight_79,79))/rzero**2
     end select
  endif

#else
  temp  = 0.
#endif
  b1bu = temp
end function b1bu


! B1bv
! ====
vectype function b1bv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  if(jadv.eq.0) then
     temp = 0.
  else
     select case(ivform)
     case(0)
        temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DP),weight_79,79) &
             + int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DP),weight_79,79)  
     case(1)
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DP),weight_79,79) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DP),weight_79,79)  
     end select
  endif

  b1bv = temp
#else
  b1bv = 0.
#endif
end function b1bv


! B1psichi
! ========
vectype function b1psichi(e,f,g)

  use basic
  use nintegrate_mod
  
  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(jadv.eq.0) then
        temp = -int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
             -  int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)
     else
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),weight_79,79)
     endif

  case(1)
     if(jadv.eq.0) then
        temp = -int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
               -int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)
     else

     endif

  end select

  b1psichi = temp
  return
end function b1psichi


! B1bchi
! =======
vectype function b1bchi(e,f,g)

  use basic
  use nintegrate_mod
  
  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

#ifdef USECOMPLEX
  if(jadv.eq.0) then
     temp = 0
  else
     temp =  &
       int4(ri3_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1),weight_79,79)  &
     + int4(ri3_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1),weight_79,79)
  endif

#else
  temp = 0.
#endif
  b1bchi = temp
  return
end function b1bchi


! B1psieta
! ========
vectype function b1psieta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(jadv.eq.0) then
!!$     temp = -(int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1 ),weight_79,79) &
!!$             +int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),weight_79,79) &
!!$             +int3(e(:,OP_1 ),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
!!$             +int3(e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DZ),weight_79,79))
!!$     if(itor.eq.1) then
!!$        temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),weight_79,79)
!!$     endif
     temp = int3(e(:,OP_1),f(:,OP_GS),g(:,OP_1),weight_79,79)

     if(hypf.ne.0.) then
        if(ihypeta.eq.1) then
           temp79a = e(:,OP_1)*g(:,OP_LP) + e(:,OP_LP)*g(:,OP_1) &
                + 2.*(e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))       
        else
           temp79a = e(:,OP_LP)
        end if

        temp = temp - int3(temp79a,f(:,OP_GS),h(:,OP_1),weight_79,79)
        if(itor.eq.1) temp = temp + &
             2.*int4(ri_79,temp79a,f(:,OP_DR),h(:,OP_1),weight_79,79)
     endif

  else
     temp = int4(ri2_79,g(:,OP_1),e(:,OP_GS),f(:,OP_GS),weight_79,79)
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),weight_79,79) &
          +int4(ri4_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),weight_79,79))
#endif
  endif

  b1psieta = temp
  return
end function b1psieta


! B1beta
! ======
vectype function b1beta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1),weight_79,79) &
       - int4(ri3_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1),weight_79,79)
  b1beta = -temp
#else
  b1beta = 0.
#endif

  return
end function b1beta



! B1psipsid
! =========
vectype function b1psipsid(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  if(jadv.eq.0) then
    temp = 0.
  else
     temp79a = e(:,OP_DZ)*(f(:,OP_DZZP)*g(:,OP_DZ ) + f(:,OP_DRZP)*g(:,OP_DR )  &
          +f(:,OP_DZP )*g(:,OP_DZZ) + f(:,OP_DRP )*g(:,OP_DRZ)) &
          +    e(:,OP_DR)*(f(:,OP_DRZP)*g(:,OP_DZ ) + f(:,OP_DRRP)*g(:,OP_DR )  &
          +f(:,OP_DZP )*g(:,OP_DRZ) + f(:,OP_DRP )*g(:,OP_DRR))
     if(itor.eq.1) then
        temp79a = temp79a -2.*ri_79*e(:,OP_DR)* &
             (f(:,OP_DZP)*g(:,OP_DZ) + f(:,OP_DRP)*g(:,OP_DR))
     endif

     temp = -0.5*int3(ri4_79,temp79a,h(:,OP_1),weight_79,79) &
          + int5(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_GS),h(:,OP_1),weight_79,79) &
          + int5(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_GS),h(:,OP_1),weight_79,79)
  endif
#else
  temp = 0
#endif
  b1psipsid = temp
  return
  
end function b1psipsid


! B1psibd
! =======
vectype function b1psibd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(jadv.eq.0) then
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79)
  else
     temp79a = f(:,OP_DZ )*g(:,OP_DR ) - f(:,OP_DR )*g(:,OP_DZ )
     temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ ) &
             + f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
     temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ ) &
             + f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)
     temp = -int4(ri3_79,e(:,OP_DZ),temp79c,h(:,OP_1 ),weight_79,79) &
          -  int4(ri3_79,e(:,OP_DR),temp79b,h(:,OP_1 ),weight_79,79) &
          -  int4(ri3_79,temp79a,e(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          -  int4(ri3_79,temp79a,e(:,OP_DR),h(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp + int4(ri4_79,e(:,OP_DR),temp79a,h(:,OP_1 ),weight_79,79)
     endif
  endif
#ifdef USECOMPLEX
  if(jadv.eq.0) then
!
  else
     temp = temp - &
          (int5(ri5_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1),weight_79,79) &
          -int5(ri5_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1),h(:,OP_1),weight_79,79))
  endif
#endif

  b1psibd = temp
  return
end function b1psibd


! B1bbd
! =====
vectype function b1bbd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  temp = 0.5* &
       (int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DP),h(:,OP_1),weight_79,79) &
       +int5(ri4_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DP),h(:,OP_1),weight_79,79) &
       +int5(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
       +int5(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1 ),h(:,OP_1),weight_79,79))

  b1bbd = temp
#else
  b1bbd = 0.
#endif
  return
end function b1bbd
! B1ped
! ====
vectype function b1ped(e,f,g)

  use basic
  use nintegrate_mod

  implicit none
  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  if(jadv.eq.0) then
    temp =  int3(e(:,OP_1),f(:,OP_DP),g(:,OP_1),weight_79,79)
  else
    temp = int4(ri2_79,e(:,OP_DR),f(:,OP_DP),g(:,OP_DR),weight_79,79) &
         + int4(ri2_79,e(:,OP_DZ),f(:,OP_DP),g(:,OP_DZ),weight_79,79)
  endif
#else
  temp = 0.
#endif

  b1ped = temp
  return
end function b1ped


! B1feta
! ======
vectype function b1feta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  if(jadv.eq.0) then
  else
  endif

  b1feta = temp
#else
  b1feta = 0.
#endif
end function b1feta


!===============================================================================
! B2 TERMS
!===============================================================================

! B2b
! ===
vectype function b2b(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  select case(ibform)
  case(0)
     temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)
     
  case(1)
     temp = int3(ri2_79,e(:,OP_1),f(:,OP_1),weight_79,79)
  end select

  b2b = temp
  return
end function b2b


! B2psieta
! ========
vectype function b2psieta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1),weight_79,79) &
       - int4(ri_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp - 2.* &
          int4(ri2_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_1),weight_79,79)
  endif

  ! hyper-resistive
  if(hypi.ne.0.) then
     ! del*(f')
     temp79a = f(:,OP_DZZP) + f(:,OP_DRRP)
     if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DRP)
     
     temp = temp  &
          + int5(ri_79,e(:,OP_DZ),temp79a,g(:,OP_1 ),h(:,OP_DR),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),temp79a,g(:,OP_1 ),h(:,OP_DZ),weight_79,79) &
          + int5(ri_79,e(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_1 ),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_1 ),weight_79,79)

     if(itor.eq.1) then
        temp = temp - 2.* &
             (int5(ri2_79,e(:,OP_1),temp79a,g(:,OP_1 ),h(:,OP_DZ),weight_79,79) &
             +int5(ri2_79,e(:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_1 ),weight_79,79))
     endif
  endif

  b2psieta = temp
#else
  b2psieta = 0.
#endif
  return
end function b2psieta



! B2beta
! ======
vectype function b2beta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ibform)
  case(0)
     temp = int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          + int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          + int3(e(:,OP_1),f(:,OP_GS),g(:,OP_1 ),weight_79,79)
     
     if(hypi.ne.0.) then
        if(ihypeta.eq.1) then
           temp79a = (e(:,OP_LP)*g(:,OP_1) + &
                e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))
           if(itor.eq.1) temp79a = temp79a + 2.*ri_79* &
                (e(:,OP_DR)*g(:,OP_1) + e(:,OP_1)*g(:,OP_DR))
        else
           temp79a = e(:,OP_LP)
           if(itor.eq.1) temp79a = temp79a + 2.*ri_79*e(:,OP_DR)
        end if
  
        temp = temp - int3(temp79a,f(:,OP_GS),h(:,OP_1),weight_79,79)

#ifdef USECOMPLEX
        temp = temp + &
             (int4(e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),h(:,OP_1),weight_79,79) &
             +int4(e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1),weight_79,79))
        if(itor.eq.1) then
           temp = temp + 2.* &
                (int5(ri_79,e(:,OP_1 ),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1),weight_79,79) &
                -int5(ri_79,e(:,OP_DR),f(:,OP_DPP ),g(:,OP_1),h(:,OP_1),weight_79,79) &
                -2.*int5(ri2_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_1),h(:,OP_1),weight_79,79))
           
        endif
#endif
     endif

  case(1)
     temp = -int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          -  int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) 

  end select

  b2beta = temp
  return
end function b2beta


! B2bu
! ====
vectype function b2bu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ibform)
  case(0)

     select case(ivform)
     case(0)
        
        temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
             - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
        
        if(itor.eq.1) then
           temp = temp - 2.*int4(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
        endif
     case(1)
        temp = (int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
             - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79))/rzero**2
     end select

  case(1)
     select case(ivform)
     case(0)
        
        temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR),weight_79,79) &
             - int4(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ),weight_79,79)
        
     case(1)
        temp = (int4(ri_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR),weight_79,79) &
             -  int4(ri_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ),weight_79,79))/rzero**2
     end select

  end select
 
  b2bu = temp
  return
end function b2bu


! B2bchi
! ======
vectype function b2bchi(e,f,g)

  use basic
  use nintegrate_mod

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ibform)
  case(0)
     temp = int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)

     if(itor.eq.1) then
        temp = temp + 2.*int4(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),weight_79,79)
     endif

  case(1)
     select case(ivform)
     case(0)
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)
     case(1)
        temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
             + int4(ri4_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)        
     end select

  end select

  b2bchi = temp
  return
end function b2bchi


! B2psiv
! ======
vectype function b2psiv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ibform)
  case(0)
     select case(ivform)
     case(0)
        temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
             - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
        if(itor.eq.1) then
           temp = temp + 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_1),weight_79,79)
        endif
     case(1)
        temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
             - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     end select

  case(1)
     select case(ivform)
     case(0)
        temp = int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
             - int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),weight_79,79)

     case(1)
        temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
             - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),weight_79,79)
     end select

  end select

  b2psiv = temp
  return
end function b2psiv


! B2psipsid
! =========
vectype function b2psipsid(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ibform)
  case(0)
     temp = int5(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.* &
             int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     endif

  case(1)
     temp = int5(ri3_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri3_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1),weight_79,79)
  end select

  b2psipsid = temp
  return 
end function b2psipsid


! B2psibd
! =======
vectype function b2psibd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  temp = &
       -(int5(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1),h(:,OP_1),weight_79,79) &
       +int5(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1),h(:,OP_1),weight_79,79))
  if(itor.eq.1) temp = temp &
       - 2.*int5(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_1),h(:,OP_1),weight_79,79)

  b2psibd = -temp
#else
  b2psibd = 0.
#endif
  return 
end function b2psibd


! B2bbd
! =====
vectype function b2bbd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ibform)
  case(0)
     temp = int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + 2.*int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79)
     endif

  case(1)
     temp = int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),h(:,OP_1),weight_79,79) &
          - int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),h(:,OP_1),weight_79,79)
     
  end select
  
  b2bbd = temp
  return 
end function b2bbd



! B2ped
! =====
vectype function b2ped(e,f,g)

  use basic
  use nintegrate_mod

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ibform)
  case(0)
     temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  case(1)
     temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),weight_79,79)
  end select

  b2ped = temp
  return
end function b2ped


! B2psif
! ======
vectype function b2psif(e,f,g)

  use basic
  use nintegrate_mod

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
  
#ifdef USECOMPLEX
  temp = int3(e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZP),weight_79,79) &
       + int3(e(:,OP_DR),f(:,OP_GS),g(:,OP_DRP),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          int4(ri_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DRP),weight_79,79)
  endif

  b2psif = temp
#else
  b2psif = 0.
#endif
  return
end function b2psif


! B2bf
! ====
vectype function b2bf(e,f,g)

  use basic
  use nintegrate_mod

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp
  
#ifdef USECOMPLEX
  temp = - &
       (int4(ri3_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DRPP),weight_79,79) &
       -int4(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZPP),weight_79,79))

  if(itor.eq.1) then
     temp = temp + 2.* &
          int4(ri4_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZPP),weight_79,79)
  endif

  b2bf = temp
#else
  b2bf = 0.
#endif
  return
end function b2bf


!=============================================================================
! B3 TERMS
!=============================================================================

! B3pe
! ====
vectype function b3pe(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)

  b3pe = temp
  return
end function b3pe


! B3psipsieta
! ===========
vectype function b3psipsieta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.eq.1.) then
     temp = 0.
  else
     temp = (gam-1.)* &
          int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),h(:,OP_1),weight_79,79)
  end if


  b3psipsieta = temp
  
  return
end function b3psipsieta


! B3bbeta
! =======
vectype function b3bbeta(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.eq.1) then
     temp = 0.
  else 
     temp = (gam-1.)* &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79))
  end if

  b3bbeta = temp
  
  return
end function b3bbeta


! B3pebd
! ======
vectype function b3pebd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
        -int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       + gam* &
       (int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
       -int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR),weight_79,79))

  b3pebd = temp
  
  return
end function b3pebd


! B3pedkappa
! ==========
vectype function b3pedkappa(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  if(gam.eq.1.) then
     b3pedkappa = 0.
     return
  end if

  temp = &
       - int4(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
       - int4(e(:,OP_DR),f(:,OP_DR),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
       - int4(e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       - int4(e(:,OP_DR),f(:,OP_1 ),g(:,OP_DR),h(:,OP_1),weight_79,79)
  
  if(hypp.ne.0.) then
     ! Laplacian[f g]
     temp79a = f(:,OP_LP)*g(:,OP_1) + f(:,OP_1)*g(:,OP_LP) &
          + 2.*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR))
     
     temp = temp - &
          (int4(temp79a,e(:,OP_LP),h(:,OP_1 ),i(:,OP_1),weight_79,79) &
          +int4(temp79a,e(:,OP_DZ),h(:,OP_DZ),i(:,OP_1),weight_79,79) &
          +int4(temp79a,e(:,OP_DR),h(:,OP_DR),i(:,OP_1),weight_79,79))
  endif

  if(iupwind.eq.1) then
     temp79a = 0.5*ri2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2) &
          +    0.5*       (cht79(:,OP_DZ)**2 + cht79(:,OP_DR)**2) &
          + ri_79*(cht79(:,OP_DZ)*pht79(:,OP_DR) &
          -cht79(:,OP_DR)*pht79(:,OP_DZ))
     temp79a = sqrt(sz79(:,OP_1)*temp79a)
     temp = temp + 0.5*int3(temp79a,f(:,OP_LP),g(:,OP_1 ),weight_79,79) &
          +        0.5*int3(temp79a,f(:,OP_1 ),g(:,OP_LP),weight_79,79) &
          +            int3(temp79a,f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          +            int3(temp79a,f(:,OP_DR),g(:,OP_DR),weight_79,79)

  endif 

  b3pedkappa = (gam-1.)*temp  
  return
end function b3pedkappa



!===============================================================================
! N1 TERMS
!===============================================================================

! N1n
! ===
vectype function n1n(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)

  n1n = temp
  return
end function n1n


! N1ndenm
! =======
vectype function n1ndenm(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,h
  real, intent(in) :: g
  vectype :: temp

!!$  temp = -g* &
!!$       (int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
!!$       +int2(e(:,OP_DR),f(:,OP_DR),weight_79,79))
  temp = g*int2(e(:,OP_1),f(:,OP_LP),weight_79,79)

#ifdef USECOMPLEX
  temp = temp + g*int3(ri2_79,e(:,OP_1),f(:,OP_DPP),weight_79,79)
#endif

  if(hypp.ne.0.) then
     temp = temp - g*int3(e(:,OP_LP),f(:,OP_LP),h(:,OP_1),weight_79,79)
  endif

  if(iupwind.eq.1) then
     temp79a = 0.5*ri2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2) &
          +    0.5*       (cht79(:,OP_DZ)**2 + cht79(:,OP_DR)**2) &
          + ri_79*(cht79(:,OP_DZ)*pht79(:,OP_DR) &
          -cht79(:,OP_DR)*pht79(:,OP_DZ))
     temp79a = sqrt(sz79(:,OP_1)*temp79a)
     temp = temp + 0.5*int2(temp79a,f(:,OP_LP),weight_79,79)
  endif 

  n1ndenm = temp
  return
end function n1ndenm



! N1nu
! ====
vectype function n1nu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(ri_79, e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79, e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     
  case(1)
     temp = int4(r_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR),weight_79,79) &
          - int4(r_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ),weight_79,79)
  end select

  n1nu = temp
  return
end function n1nu


! N1nv
! ====
vectype function n1nv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  select case(ivform)
  case(0)
     temp = -int4(ri2_79,e(:,OP_1),f(:,OP_1 ),g(:,OP_DP),weight_79,79) &
          -  int4(ri2_79,e(:,OP_1),f(:,OP_DP),g(:,OP_1 ),weight_79,79)
  case(1)
     temp = -int3(e(:,OP_1),f(:,OP_1 ),g(:,OP_DP),weight_79,79) &
          -  int3(e(:,OP_1),f(:,OP_DP),g(:,OP_1 ),weight_79,79)
  end select

  n1nv = temp
#else
  n1nv = 0.
#endif
  return
end function n1nv

! N1nchi
! ======
vectype function n1nchi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     temp = int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)

  case(1)
     temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)
  end select

  n1nchi = temp
  return
end function n1nchi


! N1s
! ===
vectype function n1s(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f
  vectype :: temp

  temp = int2(e(:,OP_1),f(:,OP_1),weight_79,79)

  n1s = temp
  return
end function n1s



!===============================================================================
! P1 TERMS
!===============================================================================

! P1pu
! ====
vectype function p1pu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  case(1)
     temp = int4(r2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR),weight_79,79) &
          - int4(r2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ),weight_79,79)

     if(itor.eq.1) then
        temp = temp + &
             2.*(gam-1.)*int3(e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
     endif
  end select

  p1pu = temp

  return
end function p1pu


! P1pv
! ====
vectype function p1pv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

#ifdef USECOMPLEX
  temp = -gam*int3(e(:,OP_1),f(:,OP_1),g(:,OP_DP),weight_79,79)
#else
  temp = 0.
#endif

  p1pv = temp

  return
end function p1pv


! P1pchi
! ======
vectype function p1pchi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     temp = gam* &
          (int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79)  &
          +int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)) &
          +(gam-1.)* &
          (int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79)  &
          +int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79))
  case(1)
     temp79a = e(:,OP_DR)*g(:,OP_DR) + e(:,OP_DZ)*g(:,OP_DZ) &
          - (gam-1.)*e(:,OP_1)*g(:,OP_GS)
     
     temp = int3(ri2_79,temp79a,f(:,OP_1),weight_79,79)
  end select

  p1pchi = temp

  return
end function p1pchi


! P1kappar
! ========
vectype function p1kappar(e,f,g,h,i,j,k)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i,j,k
  vectype :: temp

  if(gam.eq.1.) then
     p1kappar = 0.
     return
  end if

  temp79a = k(:,OP_1)*ri2_79* &
       (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)

  temp = int4(temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1 ),weight_79,79) &
       - int4(temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1 ),weight_79,79) &
       + int4(temp79a,g(:,OP_DZ),h(:,OP_1 ),i(:,OP_DR),weight_79,79) &
       - int4(temp79a,g(:,OP_DR),h(:,OP_1 ),i(:,OP_DZ),weight_79,79)

  p1kappar = (gam - 1.) * temp
  return
end function p1kappar


! P1kappax
! ========
vectype function p1kappax(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  if(gam.eq.1.) then
     p1kappax = 0.
     return
  end if

  temp79a = ri_79*i(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))

  temp79b = ri_79*i(:,OP_1)*(e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ))

  temp = (gam-1.)* &
       (int3(g(:,OP_1),h(:,OP_1),temp79a,weight_79,79)  &
       +int3(f(:,OP_1),g(:,OP_1),temp79b,weight_79,79))

  p1kappax = (gam - 1.) * temp
  return
end function p1kappax



! P1uus
! =====
vectype function p1uus(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.eq.1.) then
     p1uus = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = 0.5*(gam-1.)* &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a,weight_79,79) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a,weight_79,79))
  case(1)
     temp = 0.5*(gam-1.)* &
          (int5(r2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a,weight_79,79) &
          +int5(r2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a,weight_79,79))
  end select

  p1uus = temp
  return
end function p1uus


! P1vvs
! =====
vectype function p1vvs(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.eq.1.) then
     p1vvs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = 0.5*(gam-1.)* &
          int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_1),temp79a,weight_79,79)
  case(1)
     temp = 0.5*(gam-1.)* &
          int5(r2_79,e(:,OP_1),f(:,OP_1),g(:,OP_1),temp79a,weight_79,79)
  end select

  p1vvs = temp
  return
end function p1vvs


! P1chichis
! =========
vectype function p1chichis(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.eq.1.) then
     p1chichis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = 0.5*(gam-1.)* &
          (int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a,weight_79,79) &
          +int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a,weight_79,79))

  case(1)
     temp = 0.5*(gam-1.)* &
          (int5(ri4_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a,weight_79,79) &
          +int5(ri4_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a,weight_79,79))
  end select
     
  p1chichis = temp
  return
end function p1chichis


! P1uchis
! =======
vectype function p1uchis(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.eq.1.) then
     p1uchis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  temp = -(gam-1.)* & 
       (int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),temp79a,weight_79,79) &
       -int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),temp79a,weight_79,79))

  p1uchis = temp
  return
end function p1uchis


!======================================================================
! Parallel Viscous Terms
!======================================================================

! PVS1
! ====
subroutine PVS1(i,o)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: i
  vectype, intent(out), dimension(79) :: o
  
  o = ri_79* &
       (i(:,OP_DRZ) * (pst79(:,OP_DR)**2 - pst79(:,OP_DZ)**2) &
       +pst79(:,OP_DR)*pst79(:,OP_DZ)*(i(:,OP_DZZ) - i(:,OP_DRR)))

  if(itor.eq.1) then
     o = o + ri2_79*(pst79(:,OP_DZ) * &
          (i(:,OP_DZ)*pst79(:,OP_DZ) + i(:,OP_DR)*pst79(:,OP_DR)) &
          - i(:,OP_DZ)*bzt79(:,OP_1)**2)
  end if

  o = o * ri2_79*b2i79(:,OP_1)

end subroutine PVS1

! PVS2
! ====
subroutine PVS2(i,o)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: i
  vectype, intent(out), dimension(79) :: o

  select case(ivform)
  case(0)
     o = ri_79*bzt79(:,OP_1)* &
          (i(:,OP_DZ)*pst79(:,OP_DR) - i(:,OP_DR)*pst79(:,OP_DZ))

     if(itor.eq.1) then
        o = o + 2.*ri2_79*bzt79(:,OP_1)*pst79(:,OP_DZ)*i(:,OP_1)
     endif

  case(1)
     o = r_79*bzt79(:,OP_1)* &
          (i(:,OP_DZ)*pst79(:,OP_DR) - i(:,OP_DR)*pst79(:,OP_DZ))
  end select

  o = o * ri2_79*b2i79(:,OP_1)

end subroutine PVS2

! PVS3
! ====
subroutine PVS3(i,o)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: i
  vectype, intent(out), dimension(79) :: o

  o = i(:,OP_DZZ)*pst79(:,OP_DR)**2 + i(:,OP_DRR)*pst79(:,OP_DZ)**2 &
    - 2.*i(:,OP_DRZ)*pst79(:,OP_DR)*pst79(:,OP_DZ)

  if(itor.eq.1) then
     o = o + ri_79*i(:,OP_DR)*bzt79(:,OP_1)**2
  endif

  o = o * ri2_79*b2i79(:,OP_1) - (1./3.)*i(:,OP_LP)
 
end subroutine PVS3


! PVV1
! ====
subroutine PVV1(e,o)
  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e
  vectype, intent(out), dimension(79) :: o
  
  o =   (e(:,OP_DZZ) - e(:,OP_DRR))*pst79(:,OP_DR)*pst79(:,OP_DZ) &
       + e(:,OP_DRZ)*(pst79(:,OP_DR)**2 - pst79(:,OP_DZ)**2)

  if(itor.eq.1) then
     o = o + ri_79* &
          (pst79(:,OP_DZ)* &
          (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR)) &
          -2.*(e(:,OP_DZ)*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
               +2.*e(:,OP_DR)*pst79(:,OP_DR)*pst79(:,OP_DZ)) &
          +e(:,OP_DZ)*bzt79(:,OP_1)**2)
  endif

  o = 3.*ri_79*b2i79(:,OP_1)*o
end subroutine  PVV1

! PVV2
! ====
subroutine PVV2(e,o)
  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e
  vectype, intent(out), dimension(79) :: o

  o = 3.*ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)* &
       (e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ))
     
end subroutine  PVV2


! PVV3
! ====
subroutine PVV3(e,o)
  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e
  vectype, intent(out), dimension(79) :: o

  o = (1. - 3.*ri2_79*b2i79(:,OP_1)* &
       (pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)) * e(:,OP_LP) &
      - 3.*ri2_79*b2i79(:,OP_1)* &
      (e(:,OP_DZ)* &
       (pst79(:,OP_DZ)*pst79(:,OP_DZZ) + pst79(:,OP_DR)*pst79(:,OP_DRZ)) &
      +e(:,OP_DR)* &
       (pst79(:,OP_DZ)*pst79(:,OP_DRZ) + pst79(:,OP_DR)*pst79(:,OP_DRR)) &
      -pst79(:,OP_DZ)* &
       (pst79(:,OP_DZZ)*e(:,OP_DZ ) + pst79(:,OP_DRZ)*e(:,OP_DR ) &
       +pst79(:,OP_DZ )*e(:,OP_DZZ) + pst79(:,OP_DR )*e(:,OP_DRZ)) &
      -pst79(:,OP_DR)* &
       (pst79(:,OP_DRZ)*e(:,OP_DZ ) + pst79(:,OP_DRR)*e(:,OP_DR ) &
       +pst79(:,OP_DZ )*e(:,OP_DRZ) + pst79(:,OP_DR )*e(:,OP_DRR)))

  if(itor.eq.1) then
     o = o + 3.*ri3_79*b2i79(:,OP_1)*e(:,OP_DR)* &
          (pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2 - bzt79(:,OP_1)**2)
  endif

end subroutine  PVV3


! P1vip
! =====
vectype function P1vip(e)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e

  if(amupar.eq.0) then
     P1vip = 0.
     return
  endif

  call PVS1(pht79,temp79a)

  if(numvar.ge.2) then
     call PVS2(vzt79,temp79b)
     temp79a = temp79a + temp79b
  endif

  if(numvar.ge.3) then
     call PVS3(cht79,temp79c)
     temp79a = temp79a + temp79c
  endif

  temp79d = 3.*temp79a - cht79(:,OP_LP)

  P1vip = int4(e(:,OP_1),vip79(:,OP_1),temp79a,temp79d,weight_79,79)
end function P1vip


!======================================================================
! Gyroviscous terms
!======================================================================

! g1u
! ===
vectype function g1u(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

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
vectype function g1v(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

  select case(ivform)
  case(0)
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

  case(1)

     temp79c = e(:,OP_DRR) - e(:,OP_DZZ)
     if(itor.eq.1) temp79c = temp79c + 3.0*ri_79*e(:,OP_DR)
     temp79d = e(:,OP_DRZ)
     if(itor.eq.1) temp79d = temp79d + 3.0*ri_79*e(:,OP_DZ)

     temp79a = &
          0.25*r_79*b2i79(:,OP_1)*(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) * &
          ((pst79(:,OP_DR)*f(:,OP_DZ) + pst79(:,OP_DZ)*f(:,OP_DR))*temp79c &
          +2.*f(:,OP_DZ)*pst79(:,OP_DZ)*temp79d &
          -2.*f(:,OP_DR)*pst79(:,OP_DR)*e(:,OP_DRZ))

     temp79b = pst79(:,OP_DZ)*f(:,OP_DR) - pst79(:,OP_DR)*f(:,OP_DZ)
     temp79c = e(:,OP_DZZ)
     if(itor.eq.1) temp79c = temp79c - 2.0*ri_79*e(:,OP_DR)
     temp79d = e(:,OP_DRZ)
     if(itor.eq.1) temp79d = temp79d + 1.5*ri_79*e(:,OP_DZ)
     temp79e = e(:,OP_DRR)
     if(itor.eq.1) temp79e = temp79e +     ri_79*e(:,OP_DR)
     

     temp79a = temp79a &
          -0.75*ri_79*b2i79(:,OP_1)**2*temp79b * &
          (e(:,OP_GS)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2) - &
           2.*(   temp79c*pst79(:,OP_DR)**2 &
              -2.*temp79d*pst79(:,OP_DR)*pst79(:,OP_DZ) &
              +   temp79e*pst79(:,OP_DZ)**2))

     if(itor.eq.1) then
        temp79a = temp79a &
             + 4.5*ri2_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2*e(:,OP_DZ) * &
             (pst79(:,OP_DZ)*f(:,OP_DZ) + pst79(:,OP_DR)*f(:,OP_DR))
     endif
  end select
       
  g1v = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g1v

! g1chi
! =====
vectype function g1chi(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

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
vectype function g2u(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

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
vectype function g2v(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

  select case(ivform)
  case(0)
     temp79b = e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ)
     if(itor.eq.1) temp79b = temp79b - 2.*ri_79*e(:,OP_DZ)*f(:,OP_1)
     
     temp79a = 0.25*ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)*temp79b* &
          (1.-3.*ri2_79*b2i79(:,OP_1)* &
          (pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 - bzt79(:,OP_1)**2))

  case(1)
     temp79b = e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ)
     
     temp79a = 0.25*r_79*b2i79(:,OP_1)*bzt79(:,OP_1)*temp79b* &
          (1.-3.*ri2_79*b2i79(:,OP_1)* &
          (pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 - bzt79(:,OP_1)**2))

  end select

  g2v = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g2v


! g2chi
! =====
vectype function g2chi(e,f)

  use basic
  use nintegrate_mod

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

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
vectype function g3u(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

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
vectype function g3v(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

  select case(ivform)
  case(0)
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
     
  case(1)
     if(itor.eq.1) then
        temp79b = ri_79*e(:,OP_DR)
     else
        temp79b = 0.
     endif

   
     temp79a = &
          -0.5*b2i79(:,OP_1) * &
          (pst79(:,OP_DZ)*(e(:,OP_DZZ) - temp79b)*f(:,OP_DZ) &
          +pst79(:,OP_DR)*(e(:,OP_DRR) - temp79b)*f(:,OP_DR) &
          +e(:,OP_DRZ)*(pst79(:,OP_DZ)*f(:,OP_DR) &
                       +pst79(:,OP_DR)*f(:,OP_DZ))) &
          -1.5*ri2_79*b2i79(:,OP_1)**2 * &
          (f(:,OP_DZ)*pst79(:,OP_DR) - f(:,OP_DR)*pst79(:,OP_DZ)) * &
          (pst79(:,OP_DZ)*pst79(:,OP_DR)*(e(:,OP_DZZ) - e(:,OP_DRR)) &
          -e(:,OP_DRZ)*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)) &
          +1.5*ri2_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2 * &
          (e(:,OP_DRZ)*(pst79(:,OP_DZ)*f(:,OP_DR)  &
                       +pst79(:,OP_DR)*f(:,OP_DZ)) &
          -(e(:,OP_DZZ) - temp79b)*f(:,OP_DR)*pst79(:,OP_DR) &
          -(e(:,OP_DRR) - temp79b)*f(:,OP_DZ)*pst79(:,OP_DZ))
  end select
  g3v = int2(pit79(:,OP_1),temp79a,weight_79,79)
  return
end function g3v


! g3chi
! =====
vectype function g3chi(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f

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
vectype function qpsipsieta(e)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e
  vectype :: temp

  if(hypf.eq.0. .or. jadv.eq.1) then
     qpsipsieta = 0.
     return
  end if

  temp79a = ri2_79* &
       (jt79(:,OP_DZ)**2 + jt79(:,OP_DR)**2)
  if(itor.eq.1) then
     temp79a = temp79a - 2.*jt79(:,OP_1)* &
          (ri3_79*jt79(:,OP_DR) - ri4_79*jt79(:,OP_1))
  endif

  if(ihypeta.eq.1) then
     temp = hypf*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  else
     temp = hypf*int3(e(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  endif


  temp79a = jt79(:,OP_DZ)*ni79(:,OP_DZ) + jt79(:,OP_DR)*ni79(:,OP_DR)
  if(itor.eq.1) then
     temp79a = temp79a - ri_79*jt79(:,OP_1)*ni79(:,OP_DR)
  endif
  temp79a = temp79a * ri2_79*nt79(:,OP_1)*jt79(:,OP_1)
     
  if(ihypeta.eq.1) then
     temp = temp + hypf*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  else
     temp = temp + hypf*int3(e(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  endif

  qpsipsieta = temp
  return
end function qpsipsieta

! qbbeta
! ======
vectype function qbbeta(e)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e
  vectype :: temp

  if(hypi.eq.0.) then
     qbbeta = 0.
     return
  end if

  temp79a = ri2_79*(bzt79(:,OP_GS)*bzt79(:,OP_GS) &
       + 2.*(bzt79(:,OP_DRZ)**2 - bzt79(:,OP_DRR)*bzt79(:,OP_DZZ)))
  
  if(itor.eq.1) then 
     temp79a = temp79a + 2.* &
          (ri3_79*bzt79(:,OP_DZZ)*bzt79(:,OP_DR) &
          -ri3_79*bzt79(:,OP_DRZ)*bzt79(:,OP_DZ) &
          +ri4_79*bzt79(:,OP_DZ )*bzt79(:,OP_DZ))
  endif

  if(ihypeta.eq.1) then
     temp = hypi*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  else
     temp = hypi*int3(e(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  endif

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
     
  if(ihypeta.eq.1) then
     temp = temp + hypi*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  else
     temp = temp + hypi*int3(e(:,OP_1),temp79a,sz79(:,OP_1),weight_79,79)
  endif

  qbbeta = temp
  return
end function qbbeta


! ==============================================================
! Viscous heating terms
! ==============================================================

! quumu
! =====
vectype function quumu(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  temp = -int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),h(:,OP_1),weight_79,79) &
       + 4.* int5(ri2_79,e(:,OP_1),f(:,OP_DRR),g(:,OP_DZZ),h(:,OP_1),weight_79,79) &
       - 4.* int5(ri2_79,e(:,OP_1),f(:,OP_DRZ),g(:,OP_DRZ),h(:,OP_1),weight_79,79)

  if(itor.eq.1) then
     temp = temp &
          + 4.*int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - 4.*int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          + 4.*int5(ri4_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79)
  endif

  if(hypc.ne.0.) then
     temp79a = h(:,OP_1)*i(:,OP_1)
     temp = temp - &
          (int5(ri2_79,e(:,OP_1),vot79(:,OP_DZ),vot79(:,OP_DZ),temp79a,weight_79,79) &
          +int5(ri2_79,e(:,OP_1),vot79(:,OP_DR),vot79(:,OP_DR),temp79a,weight_79,79))
  endif

  quumu = temp
  return
end function quumu


! quchimu
! =======
vectype function quchimu(e,f,g,h,i,j)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i,j
  vectype :: temp

  quchimu = 0.
  return
end function quchimu


! qvvmu
! =====
vectype function qvvmu(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  select case(ivform)
  case(0)
     temp = - &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79))

     if(itor.eq.1) then
        temp = temp &
             + 4.*int5(ri3_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),h(:,OP_1),weight_79,79) &
             - 4.*int5(ri3_79,e(:,OP_1),f(:,OP_1 ),g(:,OP_1),h(:,OP_1),weight_79,79)
     endif
  
     if(hypv.ne.0.) then
        temp79a = h(:,OP_1)*i(:,OP_1)
        temp = temp - int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),temp79a,weight_79,79)
     endif
  case(1)
     temp = - &
          (int5(r2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          +int5(r2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79))
 
     if(hypv.ne.0.) then
        temp79a = h(:,OP_1)*i(:,OP_1)
        temp = temp - int5(r2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),temp79a,weight_79,79)
     endif

  end select

  qvvmu = temp
  return
end function qvvmu


! qchichimu
! =========
vectype function qchichimu(e,f,g,h,i)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(79,OP_NUM) :: e,f,g,h,i
  vectype :: temp

!!$  temp = -2.*int4(e(:,OP_1),cot79(:,OP_1),cot79(:,OP_1),h(:,OP_1),weight_79,79)

  temp = 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DZZ),g(:,OP_DRZ),weight_79,79) &
       - 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DRZ),g(:,OP_DZZ),weight_79,79) &
       + 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DRZ),g(:,OP_DRR),weight_79,79) &
       - 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DRR),g(:,OP_DRZ),weight_79,79)

  if(itor.eq.1) then
     temp = temp &
          + 4.*int5(ri2_79,e(:,OP_DR),h(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          - 4.*int5(ri2_79,e(:,OP_DZ),h(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)
  endif
  
  if(hypc.ne.0.) then
     temp79a = h(:,OP_1)*i(:,OP_1)
     temp = temp - 2.* &
          (int4(e(:,OP_1),cot79(:,OP_DZ),cot79(:,OP_DZ),temp79a,weight_79,79) &
          +int4(e(:,OP_1),cot79(:,OP_DR),cot79(:,OP_DR),temp79a,weight_79,79))
  endif

  qchichimu = temp
  return
end function qchichimu


!======================================================================
! ENERGY
!======================================================================

#ifdef USECOMPLEX
#define CONJUGATE(x) conjg(x)
#else
#define CONJUGATE(x) x
#endif

! Poloidal magnetic
! -----------------
real function energy_mp()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = .5* &
          (int3(ri2_79,pss79(:,OP_DZ),CONJUGATE(ps179(:,OP_DZ)),weight_79,79) &
          +int3(ri2_79,pss79(:,OP_DR),CONJUGATE(ps179(:,OP_DR)),weight_79,79) &
          +int3(ri2_79,ps179(:,OP_DZ),CONJUGATE(pss79(:,OP_DZ)),weight_79,79) &
          +int3(ri2_79,ps179(:,OP_DR),CONJUGATE(pss79(:,OP_DR)),weight_79,79))
  else
     temp = .5* &
          (int3(ri2_79,pst79(:,OP_DZ),CONJUGATE(pst79(:,OP_DZ)),weight_79,79) &
          +int3(ri2_79,pst79(:,OP_DR),CONJUGATE(pst79(:,OP_DR)),weight_79,79))
  endif

  energy_mp = temp
  return
end function energy_mp


! Toroidal magnetic
! -----------------
real function energy_mt()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = .5* &
          (int3(ri2_79,bzs79(:,OP_1),CONJUGATE(bz179(:,OP_1)),weight_79,79) &
          +int3(ri2_79,bz179(:,OP_1),CONJUGATE(bzs79(:,OP_1)),weight_79,79))
  else
     temp = .5*int3(ri2_79,bzt79(:,OP_1),CONJUGATE(bzt79(:,OP_1)),weight_79,79)
  endif

  energy_mt = temp
  return
end function energy_mt


! Pressure
! --------
real function energy_p()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(gam.eq.1.) then 
     temp = 0.
  else
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = int1(p179,weight_79,79) / (gam - 1.)
     else
        temp = int1(pt79,weight_79,79) / (gam - 1.)
     endif
  endif

  energy_p = temp
  return
end function energy_p



! Poloidal kinetic
! ----------------
real function energy_kp()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1) then
        temp = .5* &
             (int4(ri2_79,phs79(:,OP_DZ),CONJUGATE(ph179(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri2_79,phs79(:,OP_DR),CONJUGATE(ph179(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int4(ri2_79,ph179(:,OP_DZ),CONJUGATE(phs79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri2_79,ph179(:,OP_DR),CONJUGATE(phs79(:,OP_DR)),nt79(:,OP_1),weight_79,79))
     else
        temp = .5* &
             (int4(ri2_79,pht79(:,OP_DZ),CONJUGATE(pht79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri2_79,pht79(:,OP_DR),CONJUGATE(pht79(:,OP_DR)),nt79(:,OP_1),weight_79,79))
     endif

  case(1)
     if(linear.eq.1) then
        temp = .5* &
             (int4(r2_79,phs79(:,OP_DZ),CONJUGATE(ph179(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(r2_79,phs79(:,OP_DR),CONJUGATE(ph179(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int4(r2_79,ph179(:,OP_DZ),CONJUGATE(phs79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(r2_79,ph179(:,OP_DR),CONJUGATE(phs79(:,OP_DR)),nt79(:,OP_1),weight_79,79))
     else
        temp = .5* &
             (int4(r2_79,pht79(:,OP_DZ),CONJUGATE(pht79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(r2_79,pht79(:,OP_DR),CONJUGATE(pht79(:,OP_DR)),nt79(:,OP_1),weight_79,79))
     endif
  end select

  energy_kp = temp
  return
end function energy_kp


! Toroidal kinetic
! ----------------
real function energy_kt()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1) then
        temp = .5* &
             (int4(ri2_79,vzs79(:,OP_1),CONJUGATE(vz179(:,OP_1)),nt79(:,OP_1),weight_79,79) &
             +int4(ri2_79,vz179(:,OP_1),CONJUGATE(vzs79(:,OP_1)),nt79(:,OP_1),weight_79,79))
     else
        temp = .5*int4(ri2_79,vzt79(:,OP_1),CONJUGATE(vzt79(:,OP_1)),nt79(:,OP_1),weight_79,79)
     endif

  case(1)
     if(linear.eq.1) then
        temp = .5* &
             (int4(r2_79,vzs79(:,OP_1),CONJUGATE(vz179(:,OP_1)),nt79(:,OP_1),weight_79,79) &
             +int4(r2_79,vz179(:,OP_1),CONJUGATE(vzs79(:,OP_1)),nt79(:,OP_1),weight_79,79))
     else
        temp = .5*int4(r2_79,vzt79(:,OP_1),CONJUGATE(vzt79(:,OP_1)),nt79(:,OP_1),weight_79,79)
     endif
  end select

  energy_kt = temp
  return
end function energy_kt


! Compressional kinetic
! ---------------------
real function energy_k3()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1) then
        temp = .5* &
             (int3(chs79(:,OP_DZ),CONJUGATE(ch179(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int3(chs79(:,OP_DR),CONJUGATE(ch179(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int3(ch179(:,OP_DZ),CONJUGATE(chs79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int3(ch179(:,OP_DR),CONJUGATE(chs79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,chs79(:,OP_DZ),CONJUGATE(ph179(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,chs79(:,OP_DR),CONJUGATE(ph179(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,ch179(:,OP_DZ),CONJUGATE(phs79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,ch179(:,OP_DR),CONJUGATE(phs79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,CONJUGATE(chs79(:,OP_DZ)),ph179(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,CONJUGATE(chs79(:,OP_DR)),ph179(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,CONJUGATE(ch179(:,OP_DZ)),phs79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,CONJUGATE(ch179(:,OP_DR)),phs79(:,OP_DZ),nt79(:,OP_1),weight_79,79))
     else
        temp = .5* &
             (int3(cht79(:,OP_DZ),CONJUGATE(cht79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int3(cht79(:,OP_DR),CONJUGATE(cht79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,cht79(:,OP_DZ),CONJUGATE(pht79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,cht79(:,OP_DR),CONJUGATE(pht79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,CONJUGATE(cht79(:,OP_DZ)),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,CONJUGATE(cht79(:,OP_DR)),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79))
     endif
     
  case(1)
     if(linear.eq.1) then
        temp = .5* &
             (int4(ri4_79,chs79(:,OP_DZ),CONJUGATE(ch179(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri4_79,chs79(:,OP_DR),CONJUGATE(ch179(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int4(ri4_79,ch179(:,OP_DZ),CONJUGATE(chs79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri4_79,ch179(:,OP_DR),CONJUGATE(chs79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,chs79(:,OP_DZ),CONJUGATE(ph179(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,chs79(:,OP_DR),CONJUGATE(ph179(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,ch179(:,OP_DZ),CONJUGATE(phs79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,ch179(:,OP_DR),CONJUGATE(phs79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,CONJUGATE(chs79(:,OP_DZ)),ph179(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,CONJUGATE(chs79(:,OP_DR)),ph179(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,CONJUGATE(ch179(:,OP_DZ)),phs79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,CONJUGATE(ch179(:,OP_DR)),phs79(:,OP_DZ),nt79(:,OP_1),weight_79,79))
     else
        temp = .5* &
             (int4(ri4_79,cht79(:,OP_DZ),CONJUGATE(cht79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri4_79,cht79(:,OP_DR),CONJUGATE(cht79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,cht79(:,OP_DZ),CONJUGATE(pht79(:,OP_DR)),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,cht79(:,OP_DR),CONJUGATE(pht79(:,OP_DZ)),nt79(:,OP_1),weight_79,79) &
             +int4(ri_79,CONJUGATE(cht79(:,OP_DZ)),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             -int4(ri_79,CONJUGATE(cht79(:,OP_DR)),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79))
     endif
  end select

  energy_k3 = temp
  return
end function energy_k3


! Poloidal resistive
! ------------------
real function energy_mpd()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = &
          -int4(ri2_79,pss79(:,OP_GS),CONJUGATE(ps179(:,OP_GS)),eta79(:,OP_1),weight_79,79) &
          -int4(ri2_79,ps179(:,OP_GS),CONJUGATE(pss79(:,OP_GS)),eta79(:,OP_1),weight_79,79)
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,pss79(:,OP_DZP),CONJUGATE(ps179(:,OP_DZP)),eta79(:,OP_1),weight_79,79) &
          +int4(ri4_79,pss79(:,OP_DRP),CONJUGATE(ps179(:,OP_DRP)),eta79(:,OP_1),weight_79,79) &
          +int4(ri4_79,ps179(:,OP_DZP),CONJUGATE(pss79(:,OP_DZP)),eta79(:,OP_1),weight_79,79) &
          +int4(ri4_79,ps179(:,OP_DRP),CONJUGATE(pss79(:,OP_DRP)),eta79(:,OP_1),weight_79,79))
#endif
  else
     temp = -int4(ri2_79,pst79(:,OP_GS),CONJUGATE(pst79(:,OP_GS)),eta79(:,OP_1),weight_79,79)
  endif

  energy_mpd = temp
  return
end function energy_mpd


! Toroidal resistive
! ------------------
real function energy_mtd()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = - &
          (int4(ri2_79,bzs79(:,OP_DZ),CONJUGATE(bz179(:,OP_DZ)),eta79(:,OP_1),weight_79,79) &
          +int4(ri2_79,bzs79(:,OP_DR),CONJUGATE(bz179(:,OP_DR)),eta79(:,OP_1),weight_79,79) &
          +int4(ri2_79,bz179(:,OP_DZ),CONJUGATE(bzs79(:,OP_DZ)),eta79(:,OP_1),weight_79,79) &
          +int4(ri2_79,bz179(:,OP_DR),CONJUGATE(bzs79(:,OP_DR)),eta79(:,OP_1),weight_79,79))
  else
     temp = - &
          (int4(ri2_79,bzt79(:,OP_DZ),CONJUGATE(bzt79(:,OP_DZ)),eta79(:,OP_1),weight_79,79) &
          +int4(ri2_79,bzt79(:,OP_DR),CONJUGATE(bzt79(:,OP_DR)),eta79(:,OP_1),weight_79,79))
  end if

  energy_mtd = temp
  return
end function energy_mtd


! Poloidal viscous
! ----------------
real function energy_kpd()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = - &
          (int4(ri2_79,phs79(:,OP_GS),CONJUGATE(ph179(:,OP_GS)),vis79(:,OP_1),weight_79,79) &
          +int4(ri2_79,ph179(:,OP_GS),CONJUGATE(phs79(:,OP_GS)),vis79(:,OP_1),weight_79,79))
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,phs79(:,OP_DZP),CONJUGATE(ph179(:,OP_DZP)),vis79(:,OP_1),weight_79,79) &
          +int4(ri4_79,phs79(:,OP_DRP),CONJUGATE(ph179(:,OP_DRP)),vis79(:,OP_1),weight_79,79) &
          +int4(ri4_79,ph179(:,OP_DZP),CONJUGATE(phs79(:,OP_DZP)),vis79(:,OP_1),weight_79,79) &
          +int4(ri4_79,ph179(:,OP_DRP),CONJUGATE(phs79(:,OP_DRP)),vis79(:,OP_1),weight_79,79))
#endif
  else
     temp = -int4(ri2_79,pht79(:,OP_GS),CONJUGATE(pht79(:,OP_GS)),vis79(:,OP_1),weight_79,79)
  endif

  energy_kpd = temp
  return
end function energy_kpd


! Toroidal viscous
! ----------------
real function energy_ktd()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = - &
             (int4(ri2_79,vzs79(:,OP_DZ),CONJUGATE(vz179(:,OP_DZ)),vis79(:,OP_1),weight_79,79) &
             +int4(ri2_79,vzs79(:,OP_DR),CONJUGATE(vz179(:,OP_DR)),vis79(:,OP_1),weight_79,79) &
             +int4(ri2_79,vz179(:,OP_DZ),CONJUGATE(vzs79(:,OP_DZ)),vis79(:,OP_1),weight_79,79) &
             +int4(ri2_79,vz179(:,OP_DR),CONJUGATE(vzs79(:,OP_DR)),vis79(:,OP_1),weight_79,79))
     else
        temp = - &
             (int4(ri2_79,vzt79(:,OP_DZ),CONJUGATE(vzt79(:,OP_DZ)),vis79(:,OP_1),weight_79,79) &
             +int4(ri2_79,vzt79(:,OP_DR),CONJUGATE(vzt79(:,OP_DR)),vis79(:,OP_1),weight_79,79))
     endif
  case(1)
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = - &
             (int4(r2_79,vzs79(:,OP_DZ),CONJUGATE(vz179(:,OP_DZ)),vis79(:,OP_1),weight_79,79) &
             +int4(r2_79,vzs79(:,OP_DR),CONJUGATE(vz179(:,OP_DR)),vis79(:,OP_1),weight_79,79) &
             +int4(r2_79,vz179(:,OP_DZ),CONJUGATE(vzs79(:,OP_DZ)),vis79(:,OP_1),weight_79,79) &
             +int4(r2_79,vz179(:,OP_DR),CONJUGATE(vzs79(:,OP_DR)),vis79(:,OP_1),weight_79,79))
     else
        temp = - &
             (int4(r2_79,vzt79(:,OP_DZ),CONJUGATE(vzt79(:,OP_DZ)),vis79(:,OP_1),weight_79,79) &
             +int4(r2_79,vzt79(:,OP_DR),CONJUGATE(vzt79(:,OP_DR)),vis79(:,OP_1),weight_79,79))
     endif
  end select

  energy_ktd = temp
  return
end function energy_ktd

! Compressional viscous
! ---------------------
real function energy_k3d()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = - 2.* &
          (int3(chs79(:,OP_LP),CONJUGATE(ch179(:,OP_LP)),vic79(:,OP_1),weight_79,79) &
          +int3(ch179(:,OP_LP),CONJUGATE(chs79(:,OP_LP)),vic79(:,OP_1),weight_79,79))
  else
     temp = - 2.*int3(cht79(:,OP_LP),CONJUGATE(cht79(:,OP_LP)),vic79(:,OP_1),weight_79,79)
  endif

  energy_k3d = temp
  return
end function energy_k3d


! Poloidal hyper-viscous
! ----------------------
real function energy_kph()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  temp = - hypc* &
       (int5(ri2_79,vot79(:,OP_DZ),CONJUGATE(vot79(:,OP_DZ)),vis79(:,OP_1),sz79(:,OP_1),weight_79,79) &
       +int5(ri2_79,vot79(:,OP_DR),CONJUGATE(vot79(:,OP_DR)),vis79(:,OP_1),sz79(:,OP_1),weight_79,79))

  energy_kph = temp
  return
end function energy_kph


! Toroidal hyper-viscous
! ----------------------
real function energy_kth()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     temp = -hypv*int5(ri2_79,vzt79(:,OP_GS),CONJUGATE(vzt79(:,OP_GS)),vis79(:,OP_1),sz79(:,OP_1),weight_79,79)
  case(1)
     temp = -hypv*int5(r2_79,vzt79(:,OP_GS),CONJUGATE(vzt79(:,OP_GS)),vis79(:,OP_1),sz79(:,OP_1),weight_79,79)
  end select

  energy_kth = temp
  return
end function energy_kth

! Compressional hyper-viscous
! ---------------------------
real function energy_k3h()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  temp = -2.*hypc* &
       (int4(cot79(:,OP_DZ),CONJUGATE(cot79(:,OP_DZ)),vic79(:,OP_1),sz79(:,OP_1),weight_79,79) &
       +int4(cot79(:,OP_DR),CONJUGATE(cot79(:,OP_DR)),vic79(:,OP_1),sz79(:,OP_1),weight_79,79))

  energy_k3h = temp
  return
end function energy_k3h



!======================================================================
! FLUXES
!======================================================================


! Pressure convection
! -------------------
real function flux_pressure()

  use basic
  use nintegrate_mod

  implicit none

  if(numvar.lt.3 .or. gam.eq.1.) then
     flux_pressure = 0.
     return
  endif

  temp79a = pt79(:,OP_1)*cht79(:,OP_LP) &
       + pt79(:,OP_DZ)*cht79(:,OP_DZ) + pt79(:,OP_DR)*cht79(:,OP_DR) &
       + ri_79* &
       ( pt79(:,OP_DZ)*pht79(:,OP_DR) - pt79(:,OP_DR)*pht79(:,OP_DZ))

  temp79a = temp79a - ri_79*pefac*dbf* &
       ( ni79(:,OP_1)*(pet79(:,OP_DZ)*bzt79(:,OP_DR) - pet79(:,OP_DR)*bzt79(:,OP_DZ)) &
       +pet79(:,OP_1)*( ni79(:,OP_DZ)*bzt79(:,OP_DR) -  ni79(:,OP_DR)*bzt79(:,OP_DZ)))

  flux_pressure = real(-gam/(gam-1.)*int1(temp79a,weight_79,79))
  return
end function flux_pressure


! Kinetic Energy Convection
! -------------------------
real function flux_ke()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp


  ! numvar = 1
  temp79a = ri2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2)

  temp79b = ri_79*(nt79(:,OP_DZ)*pht79(:,OP_DR) - nt79(:,OP_DR)*pht79(:,OP_DZ))

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
     select case(ivform)
     case(0)
        temp79a = temp79a + ri2_79*vzt79(:,OP_1)**2
        temp79c = temp79c + &
             ri3_79*vzt79(:,OP_1) * &
             (vzt79(:,OP_DZ)*pht79(:,OP_DR) - vzt79(:,OP_DR)*pht79(:,OP_DZ))
        if(itor.eq.1) then
           temp79c = temp79c + ri4_79*vzt79(:,OP_1)**2*pht79(:,OP_DZ)
        endif
     case(1)
        temp79a = temp79a + r2_79*vzt79(:,OP_1)**2
        temp79c = temp79c + &
             r_79*vzt79(:,OP_1) * &
             (vzt79(:,OP_DZ)*pht79(:,OP_DR) - vzt79(:,OP_DR)*pht79(:,OP_DZ))
        if(itor.eq.1) then
           temp79c = temp79c - vzt79(:,OP_1)**2*pht79(:,OP_DZ)
        endif
     end select
  endif

  ! numvar = 3
  if(numvar.ge.3) then
      temp79a = temp79a &
           + cht79(:,OP_DZ)*cht79(:,OP_DZ) + cht79(:,OP_DR)*cht79(:,OP_DR) &
           + 2.*ri_79* &
           ( cht79(:,OP_DZ)*pht79(:,OP_DR) - cht79(:,OP_DR)*pht79(:,OP_DZ))
      temp79b = temp79b + nt79(:,OP_1)*cht79(:,OP_LP) &
            + nt79(:,OP_DZ)*cht79(:,OP_DZ) + nt79(:,OP_DR)*cht79(:,OP_DR)
      temp79c = temp79c + ri_79* &
           (pht79(:,OP_DR) * &
            (cht79(:,OP_DZZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRZ)*cht79(:,OP_DR )) &
           -pht79(:,OP_DZ) * &
            (cht79(:,OP_DRZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRR)*cht79(:,OP_DR ))) &
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
           + ri2_79* &
           (cht79(:,OP_DZ) * &
            (pht79(:,OP_DZZ)*pht79(:,OP_DZ ) + pht79(:,OP_DRZ)*pht79(:,OP_DR )) &
           +cht79(:,OP_DR) * &
            (pht79(:,OP_DRZ)*pht79(:,OP_DZ ) + pht79(:,OP_DRR)*pht79(:,OP_DR ))) &
           +cht79(:,OP_DZ) * &
            (cht79(:,OP_DZZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRZ)*cht79(:,OP_DR )) &
           +cht79(:,OP_DR) * &
            (cht79(:,OP_DRZ)*cht79(:,OP_DZ ) + cht79(:,OP_DRR)*cht79(:,OP_DR ))
      
      if(itor.eq.1) then
         temp79c = temp79c &
              + ri2_79*cht79(:,OP_DR)* &
              (pht79(:,OP_DZ)*cht79(:,OP_DR) - pht79(:,OP_DR)*cht79(:,OP_DZ)) &
              - ri3_79*pht79(:,OP_DZ)* &
              (pht79(:,OP_DZ)*cht79(:,OP_DR) - pht79(:,OP_DR)*cht79(:,OP_DZ)) &
              - ri3_79*cht79(:,OP_DR)* &
              (pht79(:,OP_DZ)*pht79(:,OP_DZ) + pht79(:,OP_DR)*pht79(:,OP_DR))
      endif

      select case(ivform)
      case(0)
         temp79c = temp79c + ri2_79*vzt79(:,OP_1)* &
              (vzt79(:,OP_DZ)*cht79(:,OP_DZ) + vzt79(:,OP_DR)*cht79(:,OP_DR))
         if(itor.eq.1) then
            temp79c = temp79c - ri3_79*vzt79(:,OP_1)**2*cht79(:,OP_DR)
         endif

      case(1)
         temp79c = temp79c + r2_79*vzt79(:,OP_1)* &
              (vzt79(:,OP_DZ)*cht79(:,OP_DZ) + vzt79(:,OP_DR)*cht79(:,OP_DR))
         if(itor.eq.1) then
            temp79c = temp79c + r_79*vzt79(:,OP_1)**2*cht79(:,OP_DR)
         endif

      end select
  endif

  temp = 0.5*int2(temp79a,temp79b,weight_79,79)

  temp = temp + int2(nt79(:,OP_1),temp79c,weight_79,79)

  flux_ke = real(-temp)
  return
end function flux_ke


! Poynting flux
! -------------
real function flux_poynting()

  use basic
  use nintegrate_mod

  implicit none

  temp79a = -ri2_79*pst79(:,OP_GS)*vloop/(2.*pi)

!!$  if(idens.eq.0) then
!!$     ni79 = 0.
!!$     ni79(:,OP_1) = 1.
!!$  endif

!!$  temp79a = ri3_79*(pst79(:,OP_GS)* &
!!$        (pst79(:,OP_DZ )*pht79(:,OP_DR ) - pst79(:,OP_DR )*pht79(:,OP_DZ )) &
!!$       +pst79(:,OP_DZ)* &
!!$        (pst79(:,OP_DZZ)*pht79(:,OP_DR ) - pst79(:,OP_DRZ)*pht79(:,OP_DZ ) &
!!$        +pst79(:,OP_DZ )*pht79(:,OP_DRZ) - pst79(:,OP_DR )*pht79(:,OP_DZZ)) &
!!$       +pst79(:,OP_DR)* &
!!$        (pst79(:,OP_DRZ)*pht79(:,OP_DR ) - pst79(:,OP_DRR)*pht79(:,OP_DZ) &
!!$        +pst79(:,OP_DZ )*pht79(:,OP_DRR) - pst79(:,OP_DR )*pht79(:,OP_DRZ))) &
!!$       - ri2_79* &
!!$       +(eta79(:,OP_1 )* pst79(:,OP_GS)* jt79(:,OP_1 ) &
!!$        +eta79(:,OP_1 )*(pst79(:,OP_DZ)* jt79(:,OP_DZ) + pst79(:,OP_DR)* jt79(:,OP_DR)) &
!!$        + jt79(:,OP_1 )*(pst79(:,OP_DZ)*eta79(:,OP_DZ) + pst79(:,OP_DR)*eta79(:,OP_DR)))
!!$  
!!$  if(itor.eq.1) then
!!$     temp79a = temp79a - ri4_79*pst79(:,OP_DR)* &
!!$          (pst79(:,OP_DZ)*pht79(:,OP_DR) - pst79(:,OP_DR)*pht79(:,OP_DZ))
!!$  endif
 
!!$  if(numvar.ge.2) then
!!$     temp79a = temp79a &
!!$          - ri2_79* &
!!$          (eta79(:,OP_1)* bzt79(:,OP_1) *bzt79(:,OP_GS) &
!!$          +eta79(:,OP_1)*(bzt79(:,OP_DZ)*bzt79(:,OP_DZ) + bzt79(:,OP_DR)*bzt79(:,OP_DR)) &
!!$          +bzt79(:,OP_1)*(bzt79(:,OP_DZ)*eta79(:,OP_DZ) + bzt79(:,OP_DR)*eta79(:,OP_DR))) &
!!$          - ri3_79* &
!!$          (vzt79(:,OP_1)*(bzt79(:,OP_DZ)*pst79(:,OP_DR) - bzt79(:,OP_DR)*pst79(:,OP_DZ)) &
!!$          +bzt79(:,OP_1)*(vzt79(:,OP_DZ)*pst79(:,OP_DR) - vzt79(:,OP_DR)*pst79(:,OP_DZ))) &
!!$          - ri3_79*dbf* &
!!$           (bzt79(:,OP_1)* &
!!$           ( jt79(:,OP_1 )*(ni79(:,OP_DZ)*pst79(:,OP_DR) - ni79(:,OP_DR)*pst79(:,OP_DZ)) &
!!$           + ni79(:,OP_1 )*(jt79(:,OP_DZ)*pst79(:,OP_DR) - jt79(:,OP_DR)*pst79(:,OP_DZ))) &
!!$           +( ni79(:,OP_DZ)*pst79(:,OP_DZ) +  ni79(:,OP_DR)*pst79(:,OP_DR))* &
!!$            (pst79(:,OP_DZ)*bzt79(:,OP_DR) - pst79(:,OP_DR)*bzt79(:,OP_DZ)) &
!!$           +ni79(:,OP_1)* &
!!$           (pst79(:,OP_DZ)* &
!!$            (pst79(:,OP_DZZ)*bzt79(:,OP_DR ) - pst79(:,OP_DRZ)*bzt79(:,OP_DZ ) &
!!$            +pst79(:,OP_DZ )*bzt79(:,OP_DRZ) - pst79(:,OP_DR )*bzt79(:,OP_DZZ)) &
!!$           +pst79(:,OP_DR)* &
!!$            (pst79(:,OP_DRZ)*bzt79(:,OP_DR ) - pst79(:,OP_DRR)*bzt79(:,OP_DZ ) &
!!$            +pst79(:,OP_DZ )*bzt79(:,OP_DRR) - pst79(:,OP_DR )*bzt79(:,OP_DRZ))) &
!!$           +bzt79(:,OP_1)**2*(ni79(:,OP_DZ)*bzt79(:,OP_DR) - ni79(:,OP_DR)*bzt79(:,OP_DZ))) &
!!$          +2.*ri3_79*bzt79(:,OP_1)* &
!!$           (bzt79(:,OP_DZ)*pht79(:,OP_DR) - bzt79(:,OP_DR)*pht79(:,OP_DZ))
!!$
!!$     if(itor.eq.1) then
!!$        temp79a = temp79a &
!!$             - 2.*ri4_79*pst79(:,OP_DZ)*vzt79(:,OP_1)*bzt79(:,OP_1) &
!!$             +dbf*ri4_79*ni79(:,OP_1)* &
!!$              (pst79(:,OP_DR)*(pst79(:,OP_DZ)*bzt79(:,OP_DR) - pst79(:,OP_DR)*bzt79(:,OP_DZ)) &
!!$              -2.*pst79(:,OP_DZ)*pst79(:,OP_GS)*bzt79(:,OP_1) &
!!$              -2.*bzt79(:,OP_1)**2*bzt79(:,OP_DZ))  &
!!$             +2.*ri4_79*bzt79(:,OP_1)**2*pht79(:,OP_DZ)
!!$     endif
!!$  endif

!!$  if(numvar.ge.3) then
!!$     temp79a = temp79a &
!!$          - ri_79*pefac*dbf* &
!!$            (bzt79(:,OP_1)*( ni79(:,OP_DZ)*pet79(:,OP_DR) -  ni79(:,OP_DR)*pet79(:,OP_DZ)) &
!!$            + ni79(:,OP_1)*(bzt79(:,OP_DZ)*pet79(:,OP_DR) - bzt79(:,OP_DR)*pet79(:,OP_DZ))) &
!!$          + ri2_79*bzt79(:,OP_1)* &
!!$            (bzt79(:,OP_1)*cht79(:,OP_GS) + &
!!$            2.*(bzt79(:,OP_DZ)*cht79(:,OP_DZ) + bzt79(:,OP_DR)*cht79(:,OP_DR))) &
!!$          + ri2_79* &
!!$            (pst79(:,OP_GS)*(pst79(:,OP_DZ )*cht79(:,OP_DZ ) + pst79(:,OP_DR )*cht79(:,OP_DR )) &
!!$            +pst79(:,OP_DZ)*(pst79(:,OP_DZZ)*cht79(:,OP_DZ ) + pst79(:,OP_DRZ)*cht79(:,OP_DR ) &
!!$                            +pst79(:,OP_DZ )*cht79(:,OP_DZZ) + pst79(:,OP_DR )*cht79(:,OP_DRZ)) &
!!$            +pst79(:,OP_DR)*(pst79(:,OP_DRZ)*cht79(:,OP_DZ ) + pst79(:,OP_DRR)*cht79(:,OP_DR ) &
!!$                            +pst79(:,OP_DZ )*cht79(:,OP_DRZ) + pst79(:,OP_DR )*cht79(:,OP_DRR)))
!!$          
!!$  endif

! MISSING HYPER-RESISTIVE TERMS

  flux_poynting = real(-int1(temp79a,weight_79,79))
  return
end function flux_poynting


! Heat flux
! ---------
real function flux_heat()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(numvar.lt.3) then
     flux_heat = 0.
     return
  endif

  tm79(:,OP_1  ) = pt79(:,OP_1)/nt79(:,OP_1)
  tm79(:,OP_DR ) = (pt79(:,OP_DR ) - nt79(:,OP_DR )*tm79(:,OP_1))/nt79(:,OP_1)
  tm79(:,OP_DZ ) = (pt79(:,OP_DZ ) - nt79(:,OP_DZ )*tm79(:,OP_1))/nt79(:,OP_1)
  tm79(:,OP_DRR) = (pt79(:,OP_DRR) - nt79(:,OP_DRR)*tm79(:,OP_1) &
       -2.*nt79(:,OP_DR)*tm79(:,OP_DR))/nt79(:,OP_1)
  tm79(:,OP_DRZ) = (pt79(:,OP_DRZ) - nt79(:,OP_DRZ)*tm79(:,OP_1) &
       -nt79(:,OP_DZ)*tm79(:,OP_DR)-nt79(:,OP_DR)*tm79(:,OP_DZ)) &
       /nt79(:,OP_1)
  tm79(:,OP_DZZ) = (pt79(:,OP_DZZ) - nt79(:,OP_DZZ)*tm79(:,OP_1) &
       -2.*nt79(:,OP_DZ)*tm79(:,OP_DZ))/nt79(:,OP_1)
  tm79(:,OP_LP ) = tm79(:,OP_DRR) + ri_79*tm79(:,OP_DR) + tm79(:,OP_DZZ)

  ! Isotropic heat flux
  temp79a = kap79(:,OP_1)*tm79(:,OP_LP) &
       + kap79(:,OP_DZ)*tm79(:,OP_DZ) + kap79(:,OP_DR)*tm79(:,OP_DR)

  temp = int1(temp79a,weight_79,79)


  ! Parallel heat flux
  if(kappar.ne.0.) then

     temp79b = tm79(:,OP_DZ)*pst79(:,OP_DR) - tm79(:,OP_DR)*pst79(:,OP_DZ)

     temp79c = &
          kar79(:,OP_1)*(b2i79(:,OP_DZ)*pst79(:,OP_DR)-b2i79(:,OP_DR)*pst79(:,OP_DZ)) &
         +b2i79(:,OP_1)*(kar79(:,OP_DZ)*pst79(:,OP_DR)-kar79(:,OP_DR)*pst79(:,OP_DZ))
     if(itor.eq.1) temp79c = temp79c + ri_79*kar79(:,OP_1)*b2i79(:,OP_1)*pst79(:,OP_DZ)

     temp79d =-pst79(:,OP_DZ)* &
          (tm79(:,OP_DRZ)*pst79(:,OP_DR ) - tm79(:,OP_DRR)*pst79(:,OP_DZ ) &
          +tm79(:,OP_DZ )*pst79(:,OP_DRR) - tm79(:,OP_DR )*pst79(:,OP_DRZ)) &
          +pst79(:,OP_DR)* &
          (tm79(:,OP_DZZ)*pst79(:,OP_DR ) - tm79(:,OP_DRZ)*pst79(:,OP_DZ ) &
          +tm79(:,OP_DZ )*pst79(:,OP_DRZ) - tm79(:,OP_DR )*pst79(:,OP_DZZ))

     temp = temp &
          + int3(ri2_79,temp79b,temp79c,weight_79,79) &
          + int4(kar79(:,OP_1),ri2_79,b2i79(:,OP_1),temp79d,weight_79,79)
  endif

  ! Cross heat flux
  if(kappax.ne.0.) then
     temp79a = bzt79(:,OP_1) * &
          (kax79(:,OP_DZ)*tm79(:,OP_DR) - kax79(:,OP_DR)*tm79(:,OP_DZ)) &
          +    kax79(:,OP_1) * &
          (bzt79(:,OP_DZ)*tm79(:,OP_DR) - bzt79(:,OP_DR)*tm79(:,OP_DZ))

     temp = temp + int2(ri_79,temp79a,weight_79,79)
  end if

!!$ MISSING HYPER-DIFFUSIVE TERMS

  flux_heat = real(temp)
  return
end function flux_heat


! Grav_pot
! --------
real function grav_pot()

  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     grav_pot = 0.
     return
  endif

  temp = gravr*int3(ri3_79,pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
       - gravz*int3( ri_79,pht79(:,OP_DR),nt79(:,OP_1),weight_79,79)
     
  if(numvar.ge.3) then
     temp = -gravr*int3(ri2_79,cht79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
          -gravz*int2(       cht79(:,OP_DZ),nt79(:,OP_1),weight_79,79)
  endif

  grav_pot = real(temp)
  return
end function grav_pot


!======================================================================
! Toroidal (angular) momentum
!======================================================================

! torque_em
! ~~~~~~~~~
vectype function torque_em()
  use nintegrate_mod

  implicit none

  vectype :: temp

  temp = int3(ri_79,bzt79(:,OP_DZ),pst79(:,OP_DR),weight_79,79) &
       - int3(ri_79,bzt79(:,OP_DR),pst79(:,OP_DZ),weight_79,79)
  
  torque_em = temp
end function torque_em

! torque_sol
! ~~~~~~~~~~
vectype function torque_sol()
  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int4(ri_79,pht79(:,OP_DZ),vzt79(:,OP_DR),nt79(:,OP_1 ),weight_79,79) &
          - int4(ri_79,pht79(:,OP_DR),vzt79(:,OP_DZ),nt79(:,OP_1 ),weight_79,79) &
          + int4(ri_79,pht79(:,OP_DZ),vzt79(:,OP_1 ),nt79(:,OP_DR),weight_79,79) &
          - int4(ri_79,pht79(:,OP_DR),vzt79(:,OP_1 ),nt79(:,OP_DZ),weight_79,79)
  case(1)
     temp = int4(r_79,pht79(:,OP_DZ),vzt79(:,OP_DR),nt79(:,OP_1 ),weight_79,79) &
          - int4(r_79,pht79(:,OP_DR),vzt79(:,OP_DZ),nt79(:,OP_1 ),weight_79,79) &
          + int4(r_79,pht79(:,OP_DZ),vzt79(:,OP_1 ),nt79(:,OP_DR),weight_79,79) &
          - int4(r_79,pht79(:,OP_DR),vzt79(:,OP_1 ),nt79(:,OP_DZ),weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.* &
             int3(pht79(:,OP_DZ),vzt79(:,OP_1),nt79(:,OP_1),weight_79,79)
     end if
  end select
  
  torque_sol = temp
end function torque_sol

! torque_com
! ~~~~~~~~~~
vectype function torque_com()
  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  if(numvar.lt.3) then
     torque_com = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp = int3(cht79(:,OP_DZ),vzt79(:,OP_DZ),nt79(:,OP_1 ),weight_79,79) &
          + int3(cht79(:,OP_DR),vzt79(:,OP_DR),nt79(:,OP_1 ),weight_79,79) &
          + int3(cht79(:,OP_DZ),vzt79(:,OP_1 ),nt79(:,OP_DZ),weight_79,79) &
          + int3(cht79(:,OP_DR),vzt79(:,OP_1 ),nt79(:,OP_DR),weight_79,79) &
          + int3(cht79(:,OP_LP),vzt79(:,OP_1 ),nt79(:,OP_1 ),weight_79,79)
  case(1)
     temp = int4(r2_79,cht79(:,OP_DZ),vzt79(:,OP_DZ),nt79(:,OP_1 ),weight_79,79) &
          + int4(r2_79,cht79(:,OP_DR),vzt79(:,OP_DR),nt79(:,OP_1 ),weight_79,79) &
          + int4(r2_79,cht79(:,OP_DZ),vzt79(:,OP_1 ),nt79(:,OP_DZ),weight_79,79) &
          + int4(r2_79,cht79(:,OP_DR),vzt79(:,OP_1 ),nt79(:,OP_DR),weight_79,79) &
          + int4(r2_79,cht79(:,OP_LP),vzt79(:,OP_1 ),nt79(:,OP_1 ),weight_79,79)

     if(itor.eq.1) then
        temp = temp - 2.* &
             int4(r_79,cht79(:,OP_DR),vzt79(:,OP_1),nt79(:,OP_1 ),weight_79,79) 
     endif
  end select
  
  torque_com = -temp
end function torque_com


! torque_visc
! ~~~~~~~~~~~
vectype function torque_visc()
  use basic
  use nintegrate_mod

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int2(vzt79(:,OP_DZ),vis79(:,OP_DZ),weight_79,79) &
          + int2(vzt79(:,OP_DR),vis79(:,OP_DR),weight_79,79) &
          + int2(vzt79(:,OP_GS),vis79(:,OP_1 ),weight_79,79)

     if(itor.eq.1) then
        temp = temp - 2.*int3(ri_79,vzt79(:,OP_1),vis79(:,OP_DR),weight_79,79)
     endif
  case(1)
     temp = int3(r2_79,vzt79(:,OP_DZ),vis79(:,OP_DZ),weight_79,79) &
          + int3(r2_79,vzt79(:,OP_DR),vis79(:,OP_DR),weight_79,79) &
          + int3(r2_79,vzt79(:,OP_GS),vis79(:,OP_1 ),weight_79,79)

     if(itor.eq.1) then
        temp = temp + 4.*int3(r_79,vzt79(:,OP_DR),vis79(:,OP_1),weight_79,79)
     endif
  end select

  torque_visc = temp
end function torque_visc


! torque_gyro
! ~~~~~~~~~~~
vectype function torque_gyro(itri)
  use basic
  use arrays
  use t_data
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri
  vectype, dimension(20) :: avec
  vectype :: temp

  if(gyro.eq.0 .or. numvar.lt.2) then
     torque_gyro = 0.
     return
  endif

  call calcavector(itri, gyro_tau, 1, 1, avec)
  call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, temp79a)

  torque_gyro = int1(temp79a,weight_79,79)

end function torque_gyro


! torque_denm
! ~~~~~~~~~~~
vectype function torque_denm()
  use basic
  use nintegrate_mod

  implicit none

  if(denm.eq.0 .or. idens.eq.0) then
     torque_denm = 0.
     return
  endif

  select case(ivform)
  case(0)
     torque_denm = denm*int2(nt79(:,OP_LP),vzt79(:,OP_1),weight_79,79)
  case(1)
     torque_denm = denm*int3(r2_79,nt79(:,OP_LP),vzt79(:,OP_1),weight_79,79)
  end select

end function torque_denm


end module metricterms_new



module two_fluid

implicit none

contains
vectype function v1hupsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = int3(g(:,OP_GS),e(:,OP_DZ),f(:,OP_DZP)) &
             + int3(g(:,OP_GS),e(:,OP_DR),f(:,OP_DRP)) &
             - int3(f(:,OP_LP),e(:,OP_DZ),g(:,OP_DZP)) &
             - int3(f(:,OP_LP),e(:,OP_DR),g(:,OP_DRP)) 
     end if
#endif

  v1hupsi = temp
  return
end function v1hupsi

vectype function v1hub(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = int4(g(:,OP_1),ri_79,f(:,OP_DZPP),e(:,OP_DR))   &
             - int4(g(:,OP_1),ri_79,f(:,OP_DRPP),e(:,OP_DZ))
     end if
#endif

  v1hub = temp
  return
end function v1hub

vectype function v1huf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = -int4(r_79,f(:,OP_LP),e(:,OP_DZ),g(:,OP_DRPP))  &
               +int4(r_79,f(:,OP_LP),e(:,OP_DR),g(:,OP_DZPP))
     end if
#endif

  v1huf = temp
  return
end function v1huf

vectype function v1hvpsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

     if(surface_int) then 
        temp = 0.
     else
        temp79a = r_79*f(:,OP_LP) 
        if(itor.eq.1) temp79a = temp79a + 2.*f(:,OP_DR)
        temp79b = r_79*(e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))
        if(itor.eq.1) temp79b = temp79b + 2.*f(:,OP_1)*e(:,OP_DZ)
        temp = int3(temp79a,e(:,OP_DZ),g(:,OP_DR)) &
             - int3(temp79a,e(:,OP_DR),g(:,OP_DZ)) &
             + int2(g(:,OP_GS),temp79b)
     end if

  v1hvpsi = temp
  return
end function v1hvpsi

vectype function v1hvb(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = int3(g(:,OP_1),f(:,OP_DRP),e(:,OP_DR))  &
              +int3(g(:,OP_1),f(:,OP_DZP),e(:,OP_DZ))  &
          + 2.*int4(g(:,OP_1),f(:,OP_DP),ri_79,e(:,OP_DR))
     end if
#endif

  v1hvb = temp
  return
end function v1hvb

vectype function v1hvf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp79a = r2_79*f(:,OP_LP) + 2.*r_79*f(:,OP_DR)
        temp = -int3(temp79a,e(:,OP_DR),g(:,OP_DRP))   &
               -int3(temp79a,e(:,OP_DZ),g(:,OP_DZP))
     end if
#endif

  v1hvf = temp
  return
end function v1hvf

vectype function v1hchipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = -int4(ri3_79,f(:,OP_GSP),e(:,OP_DZ),g(:,OP_DR )) &
               +int4(ri3_79,f(:,OP_GSP),e(:,OP_DR),g(:,OP_DZ )) &
               -int4(ri3_79,g(:,OP_GS) ,e(:,OP_DZ),f(:,OP_DRP)) &
               +int4(ri3_79,g(:,OP_GS) ,e(:,OP_DR),f(:,OP_DZP)) &
            + 2.*int4(ri4_79,f(:,OP_DZ) ,e(:,OP_DR),g(:,OP_DRP))  &
            + 2.*int4(ri4_79,f(:,OP_DZ) ,e(:,OP_DZ),g(:,OP_DZP))
               
     end if
#endif

  v1hchipsi = temp
  return
end function v1hchipsi

vectype function v1hchib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp =  -int4(ri4_79,g(:,OP_1),f(:,OP_DRPP),e(:,OP_DR)) &
                -int4(ri4_79,g(:,OP_1),f(:,OP_DZPP),e(:,OP_DZ))

     end if
#endif
  v1hchib = temp
  return
end function v1hchib

vectype function v1hchif(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = 2.*int4(ri3_79,f(:,OP_DZ),e(:,OP_DZ),g(:,OP_DRPP))  &
              -2.*int4(ri3_79,f(:,OP_DZ),e(:,OP_DR),g(:,OP_DZPP))  &
                 +int4(ri2_79,f(:,OP_GSP),e(:,OP_DR),g(:,OP_DRP))  &
                 +int4(ri2_79,f(:,OP_GSP),e(:,OP_DZ),g(:,OP_DZP))
     end if
#endif

  v1hchif = temp
  return
end function v1hchif

vectype function v2hupsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.

     if(surface_int) then 
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp = int4(e(:,OP_1 ),ri_79,g(:,OP_DZP),f(:,OP_DRP)) &
             - int4(e(:,OP_1 ),ri_79,g(:,OP_DRP),f(:,OP_DZP)) &
             + int4(e(:,OP_1 ),ri_79,g(:,OP_DZ),f(:,OP_DRPP)) &
             - int4(e(:,OP_1 ),ri_79,g(:,OP_DR),f(:,OP_DZPP)) 
#endif
        temp = temp                                           &
             + int4(f(:,OP_LP),r_79,e(:,OP_DZ),g(:,OP_DR)) &
             - int4(f(:,OP_LP),r_79,e(:,OP_DR),g(:,OP_DZ)) 
     end if

  v2hupsi = temp
  return
end function v2hupsi

vectype function v2hub(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = int3(e(:,OP_1),f(:,OP_DRP),g(:,OP_DR))   &
             + int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ))
     end if
#endif

  v2hub = temp
  return
end function v2hub

vectype function v2huf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then 
        temp = 0.
     else
        temp = int3(e(:,OP_1),f(:,OP_DRPP),g(:,OP_DRP))   &
             + int3(e(:,OP_1),f(:,OP_DZPP),g(:,OP_DZP))   &
             + int3(e(:,OP_1),f(:,OP_DRP),g(:,OP_DRPP))   &
             + int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_DZPP))   &
             - int4(r2_79,f(:,OP_LP),e(:,OP_DR),g(:,OP_DRP)) &
             - int4(r2_79,f(:,OP_LP),e(:,OP_DZ),g(:,OP_DZP)) &
             - int4(r2_79,f(:,OP_LP),e(:,OP_1), g(:,OP_LPP)) 

     end if
#endif

  v2huf = temp
  return
end function v2huf

vectype function v2hvpsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = - int3(e(:,OP_1),f(:,OP_DRP),g(:,OP_DR )) &
               - int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ )) &
               - int3(e(:,OP_1),f(:,OP_DR ),g(:,OP_DRP)) &
               - int3(e(:,OP_1),f(:,OP_DZ ),g(:,OP_DZP)) &
               - 2*int4(e(:,OP_1),ri_79,f(:,OP_DP),g(:,OP_DR )) &
               - 2*int4(e(:,OP_1),ri_79,f(:,OP_1) ,g(:,OP_DRP)) 

     end if
#endif
  v2hvpsi = temp
  return
end function v2hvpsi

vectype function v2hvb(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

     if(surface_int) then 
        temp = 0.
     else
        temp = int4(e(:,OP_1),r_79,g(:,OP_DZ),f(:,OP_DR))  &
              -int4(e(:,OP_1),r_79,g(:,OP_DR),f(:,OP_DZ))
        if(itor.eq.1) temp = temp - 2.*int3(e(:,OP_1),f(:,OP_1),g(:,OP_DZ))
     end if

  v2hvb = temp
  return
end function v2hvb

vectype function v2hvf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = int4(e(:,OP_1),r_79,g(:,OP_DZPP),f(:,OP_DR))  &
              -int4(e(:,OP_1),r_79,g(:,OP_DRPP),f(:,OP_DZ))  &
           -2.*int3(e(:,OP_1),f(:,OP_1),g(:,OP_DZPP))        &
             + int4(e(:,OP_1),r_79,g(:,OP_DZP),f(:,OP_DRP))  &
              -int4(e(:,OP_1),r_79,g(:,OP_DRP),f(:,OP_DZP))  &
           -2.*int3(e(:,OP_1),f(:,OP_DP),g(:,OP_DZP))
     end if
#endif
  v2hvf = temp
  return
end function v2hvf

vectype function v2hchipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
      temp = 0.

     if(surface_int) then 
        temp = 0.
     else
        if(itor.eq.1) then
          temp = -2*int4(ri3_79,f(:,OP_DZ),e(:,OP_DZ),g(:,OP_DR )) &
                 +2*int4(ri3_79,f(:,OP_DZ),e(:,OP_DR),g(:,OP_DZ )) 
        endif
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp                                             &
               +int4(ri4_79,e(:,OP_1) ,g(:,OP_DR),f(:,OP_DRPP)) &
               +int4(ri4_79,e(:,OP_1) ,g(:,OP_DZ),f(:,OP_DZPP)) &
               +int4(ri4_79,e(:,OP_1) ,g(:,OP_DRP),f(:,OP_DRP)) &
               +int4(ri4_79,e(:,OP_1) ,g(:,OP_DZP),f(:,OP_DZP)) 
#endif
               
     end if

  v2hchipsi = temp
  return
end function v2hchipsi

vectype function v2hchib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp =  -int4(ri3_79,e(:,OP_1),g(:,OP_DZ),f(:,OP_DRP)) &
                +int4(ri3_79,e(:,OP_1),g(:,OP_DR),f(:,OP_DZP))

     end if
#endif

  v2hchib = temp
  return
end function v2hchib

vectype function v2hchif(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = 2.*int4(ri2_79,f(:,OP_DZ) ,e(:,OP_DR),g(:,OP_DRP))  &
              +2.*int4(ri2_79,f(:,OP_DZ) ,e(:,OP_DZ),g(:,OP_DZP))  &
              +2.*int4(ri2_79,f(:,OP_DZ) ,e(:,OP_1), g(:,OP_LPP))  &
                 -int4(ri3_79,e(:,OP_1),g(:,OP_DZPP),f(:,OP_DRP))  &
                 +int4(ri3_79,e(:,OP_1),g(:,OP_DRPP),f(:,OP_DZP))  &
                 -int4(ri3_79,e(:,OP_1),g(:,OP_DZP),f(:,OP_DRPP))  &
                 +int4(ri3_79,e(:,OP_1),g(:,OP_DRP),f(:,OP_DZPP))
     end if
#endif

  v2hchif = temp
  return
end function v2hchif
vectype function v3hupsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = int4(ri3_79,f(:,OP_LP),e(:,OP_DZ),g(:,OP_DRP)) &
             - int4(ri3_79,f(:,OP_LP),e(:,OP_DR),g(:,OP_DZP)) &
             + int4(ri3_79,f(:,OP_GS),f(:,OP_DZP),e(:,OP_DR)) &
             - int4(ri3_79,f(:,OP_GS),f(:,OP_DRP),e(:,OP_DZ)) 
     end if
#endif

  v3hupsi = temp
  return
end function v3hupsi

vectype function v3hub(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

     if(surface_int) then 
        temp = 0.
     else
        temp =  int4(g(:,OP_1),ri2_79,f(:,OP_LP),e(:,OP_LP))   
        if(itor.eq.1) temp = temp - 4.*int4(g(:,OP_1),ri3_79,f(:,OP_LP),e(:,OP_DR))   
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp                                            &
             +  int4(g(:,OP_DP),ri4_79,f(:,OP_DRP),e(:,OP_DR)) &
             +  int4(g(:,OP_DP),ri4_79,f(:,OP_DZP),e(:,OP_DZ))
#endif
     end if

  v3hub = temp
  return
end function v3hub

vectype function v3huf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = -int4(ri2_79,f(:,OP_LP),e(:,OP_DR),g(:,OP_DRPP))  &
               -int4(ri2_79,f(:,OP_LP),e(:,OP_DZ),g(:,OP_DZPP))
     end if
#endif
  v3huf = temp
  return
end function v3huf

vectype function v3hvpsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

     if(surface_int) then 
        temp = 0.
     else
        temp79a = ri2_79*f(:,OP_LP) 
        if(itor.eq.1) temp79a = temp79a + 2.*ri3_79*f(:,OP_DR)
        temp79b = ri2_79*(e(:,OP_DR)*f(:,OP_DR) + e(:,OP_DZ)*f(:,OP_DZ))  
        if(itor.eq.1) temp79b = temp79b + 2.*ri3_79*f(:,OP_1)*e(:,OP_DR)
        temp = int3(temp79a,e(:,OP_DR),g(:,OP_DR)) &
             + int3(temp79a,e(:,OP_DZ),g(:,OP_DZ)) &
             + int2(g(:,OP_GS),temp79b)
     end if

  v3hvpsi = temp
  return
end function v3hvpsi

vectype function v3hvb(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = int4(ri3_79,g(:,OP_DP),e(:,OP_DZ),f(:,OP_DR))  &
              -int4(ri3_79,g(:,OP_DP),e(:,OP_DR),f(:,OP_DZ))  &
          + 2.*int4(ri4_79,g(:,OP_DP),e(:,OP_DZ),f(:,OP_1))
     end if
#endif

  v3hvb = temp
  return
end function v3hvb

vectype function v3hvf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then 
        temp = 0.
     else
        temp79a = f(:,OP_DR) + 2.*ri_79*f(:,OP_1)
        temp = -int4(temp79a,ri_79 ,e(:,OP_DRZ),g(:,OP_DRP))   &
               +int4(temp79a,ri_79 ,e(:,OP_DRR),g(:,OP_DZP))   &
               -int4(temp79a,ri_79 ,e(:,OP_DZ),g(:,OP_DRRP))   &
               +int4(temp79a,ri_79 ,e(:,OP_DR),g(:,OP_DRZP))   &

               +int4(temp79a,ri2_79,e(:,OP_DZ),g(:,OP_DRP))   &
               -int4(temp79a,ri2_79,e(:,OP_DR),g(:,OP_DZP))   &

               -int3(f(:,OP_DZ),e(:,OP_DZZ),g(:,OP_DRP))  &
               +int3(f(:,OP_DZ),e(:,OP_DRZ),g(:,OP_DZP))  &
               -int3(f(:,OP_DZ),e(:,OP_DZ),g(:,OP_DRZP))  &
               +int3(f(:,OP_DZ),e(:,OP_DR),g(:,OP_DZZP))  
     end if
#endif
  v3hvf = temp
  return
end function v3hvf

vectype function v3hchipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then 
        temp = 0.
     else
        temp = -int4(ri6_79,f(:,OP_GSP),e(:,OP_DR),g(:,OP_DR )) &
               -int4(ri6_79,f(:,OP_GSP),e(:,OP_DZ),g(:,OP_DZ )) &
               -int4(ri6_79,g(:,OP_GS) ,e(:,OP_DZ),f(:,OP_DZP)) &
               -int4(ri6_79,g(:,OP_GS) ,e(:,OP_DR),f(:,OP_DRP)) &
            - 2*int4(ri7_79,f(:,OP_DZ),e(:,OP_DZ),g(:,OP_DRP))  &
            + 2*int4(ri7_79,f(:,OP_DZ),e(:,OP_DR),g(:,OP_DZP))
               
     end if
#endif

  v3hchipsi = temp
  return
end function v3hchipsi

vectype function v3hchib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.

     if(surface_int) then 
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp =  -int4(ri7_79,g(:,OP_DP),e(:,OP_DZ),f(:,OP_DRP)) &
                +int4(ri7_79,g(:,OP_DP),e(:,OP_DR),f(:,OP_DZP)) 
#endif
        if(itor.eq.1) then
          temp = temp                                              &
                  -2.*int4(ri6_79,g(:,OP_1),f(:,OP_DZ),e(:,OP_LP)) &
                  +8.*int4(ri7_79,g(:,OP_1),f(:,OP_DZ),e(:,OP_DR))
        endif

     end if

  v3hchib = temp
  return
end function v3hchib

vectype function v3hchif(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

     if(surface_int) then 
        temp = 0.
     else
        temp = 2.*int4(ri6_79,f(:,OP_DZ),e(:,OP_DR),g(:,OP_DRPP))  &
              +2.*int4(ri6_79,f(:,OP_DZ),e(:,OP_DZ),g(:,OP_DZPP))  &

                 +int4(ri5_79,f(:,OP_DRP),e(:,OP_DRZ),g(:,OP_DRP))  &
                 -int4(ri5_79,f(:,OP_DRP),e(:,OP_DRR),g(:,OP_DZP))  &
                 +int4(ri5_79,f(:,OP_DRP),e(:,OP_DZ),g(:,OP_DRRP))  &
                 -int4(ri5_79,f(:,OP_DRP),e(:,OP_DR),g(:,OP_DRZP))  &

                 -int4(ri6_79,f(:,OP_DRP),e(:,OP_DZ),g(:,OP_DRP))  &
                 +int4(ri6_79,f(:,OP_DRP),e(:,OP_DR),g(:,OP_DZP))  &

                 +int4(ri5_79,f(:,OP_DZP),e(:,OP_DZZ),g(:,OP_DRP))  &
                 -int4(ri5_79,f(:,OP_DZP),e(:,OP_DRZ),g(:,OP_DZP))  &
                 +int4(ri5_79,f(:,OP_DZP),e(:,OP_DZ),g(:,OP_DRZP))  &
                 -int4(ri5_79,f(:,OP_DZP),e(:,OP_DR),g(:,OP_DZZP))  
     end if
#endif
  v3hchif = temp
  return
end function v3hchif
end module two_fluid

module two_fluid

implicit none

contains

function v1hupsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hupsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx3(e(:,:,OP_DZ),g(:,OP_GS),f(:,OP_DZP)) &
       + intx3(e(:,:,OP_DR),g(:,OP_GS),f(:,OP_DRP)) &
       - intx3(e(:,:,OP_DZ),f(:,OP_LP),g(:,OP_DZP)) &
       - intx3(e(:,:,OP_DR),f(:,OP_LP),g(:,OP_DRP)) 

  v1hupsi = temp
#else
  v1hupsi = 0.
#endif

end function v1hupsi

function v1hub(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hub
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),g(:,OP_1),ri_79,f(:,OP_DZPP)) &
       - intx4(e(:,:,OP_DZ),g(:,OP_1),ri_79,f(:,OP_DRPP))

  v1hub = temp
#else
  v1hub = 0.
#endif
end function v1hub


function v1huf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1huf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  
  temp = -intx4(e(:,:,OP_DZ),r_79,f(:,OP_LP),g(:,OP_DRP))  &
       +  intx4(e(:,:,OP_DR),r_79,f(:,OP_LP),g(:,OP_DZP))

  v1huf = temp
#else
  v1huf = 0.
#endif
end function v1huf

function v1hvpsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hvpsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp79a = r_79*f(:,OP_LP) 
  if(itor.eq.1) temp79a = temp79a + 2.*f(:,OP_DR)

  temp = intx3(e(:,:,OP_DZ),temp79a,g(:,OP_DR)) &
       - intx3(e(:,:,OP_DR),temp79a,g(:,OP_DZ)) &
       + intx4(e(:,:,OP_DZ),f(:,OP_DR),r_79,g(:,OP_GS)) &
       - intx4(e(:,:,OP_DR),f(:,OP_DZ),r_79,g(:,OP_GS))
  if(itor.eq.1) then
     temp = temp + 2.*intx3(e(:,:,OP_DZ),f(:,OP_1),g(:,OP_GS))
  end if
  
  v1hvpsi = temp
end function v1hvpsi

function v1hvb(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hvb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx3(e(:,:,OP_DR),g(:,OP_1),f(:,OP_DRP))  &
       + intx3(e(:,:,OP_DZ),g(:,OP_1),f(:,OP_DZP))  &
       + 2.*intx4(e(:,:,OP_DR),g(:,OP_1),f(:,OP_DP),ri_79)

  v1hvb = temp
#else
  v1hvb = 0.
#endif
end function v1hvb

function v1hvf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hvf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp79a = r2_79*f(:,OP_LP) + 2.*r_79*f(:,OP_DR)
  temp = -intx3(e(:,:,OP_DR),temp79a,g(:,OP_DR))   &
       -  intx3(e(:,:,OP_DZ),temp79a,g(:,OP_DZ))

  v1hvf = temp
#else
  v1hvf = 0.
#endif
end function v1hvf

function v1hchipsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hchipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = -  intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_GSP),g(:,OP_DR )) &
       +    intx4(e(:,:,OP_DR),ri3_79,f(:,OP_GSP),g(:,OP_DZ )) &
       + 2.*intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DR),f(:,OP_DRP))  &
       - 2.*intx4(e(:,:,OP_DR),ri4_79,g(:,OP_DZ),f(:,OP_DRP))  &
       -    intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_GS),f(:,OP_DRP))  &
       +    intx4(e(:,:,OP_DR),ri3_79,g(:,OP_GS),f(:,OP_DZP))  &
       + 2.*intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DZ),g(:,OP_DRP))  &
       + 2.*intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_DZP))
     
  v1hchipsi = temp
#else
  v1hchipsi = 0.
#endif
end function v1hchipsi

function v1hchib(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hchib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = -intx4(e(:,:,OP_DR),ri4_79,g(:,OP_1),f(:,OP_DRPP)) &
       -  intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_1),f(:,OP_DZPP))

  v1hchib = temp
#else
  v1hchib = 0.
#endif

end function v1hchib

function v1hchif(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1hchif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = 2.*intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DZ),g(:,OP_DRP))  &
       - 2.*intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZ),g(:,OP_DZP))  &
       +    intx4(e(:,:,OP_DR),ri2_79,f(:,OP_GSP),g(:,OP_DR))  &
       +    intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_GSP),g(:,OP_DZ))

  v1hchif = temp
#else
  v1hchif = 0.
#endif

end function v1hchif

function v2hupsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hupsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx4(e(:,:,OP_1 ),ri_79,g(:,OP_DZP),f(:,OP_DRP)) &
       - intx4(e(:,:,OP_1 ),ri_79,g(:,OP_DRP),f(:,OP_DZP)) &
       + intx4(e(:,:,OP_1 ),ri_79,g(:,OP_DZ),f(:,OP_DRPP)) &
       - intx4(e(:,:,OP_1 ),ri_79,g(:,OP_DR),f(:,OP_DZPP)) 
#endif
  temp = temp                                           &
       + intx4(e(:,:,OP_DZ),f(:,OP_LP),r_79,g(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),f(:,OP_LP),r_79,g(:,OP_DZ)) 

  v2hupsi = temp
end function v2hupsi

function v2hub(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hub
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx3(e(:,:,OP_1),f(:,OP_DRP),g(:,OP_DR))   &
       + intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_DZ))

  v2hub = temp
#else
  v2hub = 0.
#endif
end function v2hub

function v2huf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2huf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx3(e(:,:,OP_1),f(:,OP_DRPP),g(:,OP_DR))   &
       + intx3(e(:,:,OP_1),f(:,OP_DZPP),g(:,OP_DZ))   &
       + intx3(e(:,:,OP_1),f(:,OP_DRP),g(:,OP_DRP))   &
       + intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_DZP))   &
       - intx4(e(:,:,OP_DR),r2_79,f(:,OP_LP),g(:,OP_DR)) &
       - intx4(e(:,:,OP_DZ),r2_79,f(:,OP_LP),g(:,OP_DZ)) &
       - intx4(e(:,:,OP_1),r2_79,f(:,OP_LP), g(:,OP_LP)) 

  v2huf = temp
#else
  v2huf = 0.
#endif
end function v2huf

function v2hvpsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hvpsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = - intx3(e(:,:,OP_1),f(:,OP_DRP),g(:,OP_DR )) &
       - intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_DZ )) &
       - intx3(e(:,:,OP_1),f(:,OP_DR ),g(:,OP_DRP)) &
       - intx3(e(:,:,OP_1),f(:,OP_DZ ),g(:,OP_DZP)) &
       - 2*intx4(e(:,:,OP_1),ri_79,f(:,OP_DP),g(:,OP_DR )) &
       - 2*intx4(e(:,:,OP_1),ri_79,f(:,OP_1) ,g(:,OP_DRP)) 

  v2hvpsi = temp
#else
  v2hvpsi = 0.
#endif
end function v2hvpsi

function v2hvb(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hvb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_1),r_79,g(:,OP_DZ),f(:,OP_DR))  &
       -intx4(e(:,:,OP_1),r_79,g(:,OP_DR),f(:,OP_DZ))
  if(itor.eq.1) temp = temp + 2.*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DZ))

  v2hvb = temp
end function v2hvb

function v2hvf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hvf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp =   intx4(e(:,:,OP_1),r_79,g(:,OP_DZP),f(:,OP_DR))  &
       -   intx4(e(:,:,OP_1),r_79,g(:,OP_DRP),f(:,OP_DZ))  &
       +2.*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DZP))        &
       +   intx4(e(:,:,OP_1),r_79,g(:,OP_DZ),f(:,OP_DRP))  &
       -   intx4(e(:,:,OP_1),r_79,g(:,OP_DR),f(:,OP_DZP))  &
       +2.*intx3(e(:,:,OP_1),f(:,OP_DP),g(:,OP_DZ))

  v2hvf = temp
#else
  v2hvf = 0.
#endif
end function v2hvf

function v2hchipsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hchipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  if(itor.eq.1) then
     temp = -2*intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DZ),g(:,OP_DR )) &
          +  2*intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZ),g(:,OP_DZ )) 
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp                                             &
       +intx4(e(:,:,OP_1),ri4_79,g(:,OP_DR),f(:,OP_DRPP)) &
       +intx4(e(:,:,OP_1),ri4_79,g(:,OP_DZ),f(:,OP_DZPP)) &
       +intx4(e(:,:,OP_1),ri4_79,g(:,OP_DRP),f(:,OP_DRP)) &
       +intx4(e(:,:,OP_1),ri4_79,g(:,OP_DZP),f(:,OP_DZP)) 
#endif

  v2hchipsi = temp
end function v2hchipsi

function v2hchib(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hchib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g


#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp =-intx4(e(:,:,OP_1),ri3_79,g(:,OP_DZ),f(:,OP_DRP)) &
       + intx4(e(:,:,OP_1),ri3_79,g(:,OP_DR),f(:,OP_DZP))

  v2hchib = temp
#else
  v2hchib = 0.
#endif
end function v2hchib

function v2hchif(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2hchif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp =2.*intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DZ), g(:,OP_DR))  &
       +2.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ) ,g(:,OP_DZ))  &
       +2.*intx4(e(:,:,OP_1), ri2_79,f(:,OP_DZ) ,g(:,OP_LP))  &
       -   intx4(e(:,:,OP_1),ri3_79,g(:,OP_DZP),f(:,OP_DRP))  &
       +   intx4(e(:,:,OP_1),ri3_79,g(:,OP_DRP),f(:,OP_DZP))  &
       -   intx4(e(:,:,OP_1),ri3_79,g(:,OP_DZ),f(:,OP_DRPP))  &
       +   intx4(e(:,:,OP_1),ri3_79,g(:,OP_DR),f(:,OP_DZPP))

  v2hchif = temp
#else
  v2hchif = 0.
#endif
end function v2hchif


function v3hupsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hupsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_LP),g(:,OP_DRP)) &
       - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_LP),g(:,OP_DZP)) &
       + intx4(e(:,:,OP_DR),ri3_79,f(:,OP_GS),f(:,OP_DZP)) &
       - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_GS),f(:,OP_DRP)) 

  v3hupsi = temp
#else
  v3hupsi = 0.
#endif
end function v3hupsi

function v3hub(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hub
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp =  intx4(e(:,:,OP_LP),g(:,OP_1),ri2_79,f(:,OP_LP))   
  if(itor.eq.1) temp = temp &
       - 4.*intx4(e(:,:,OP_DR),g(:,OP_1),ri3_79,f(:,OP_LP))   
#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp                                            &
       +  intx4(e(:,:,OP_DR),g(:,OP_DP),ri4_79,f(:,OP_DRP)) &
       +  intx4(e(:,:,OP_DZ),g(:,OP_DP),ri4_79,f(:,OP_DZP))
#endif

  v3hub = temp
end function v3hub

function v3huf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3huf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = -intx4(e(:,:,OP_DR),ri2_79,f(:,OP_LP),g(:,OP_DRP)) &
       -  intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_LP),g(:,OP_DZP))

  v3huf = temp
#else
  v3huf = 0.
#endif
end function v3huf

function v3hvpsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hvpsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp79a = ri2_79*f(:,OP_LP) 
  if(itor.eq.1) temp79a = temp79a + 2.*ri3_79*f(:,OP_DR)
  temp = intx3(e(:,:,OP_DR),temp79a,g(:,OP_DR)) &
       + intx3(e(:,:,OP_DZ),temp79a,g(:,OP_DZ)) &
       + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_GS)) &
       + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_GS))
  if(itor.eq.1) then
     temp = temp &
          + 2.*intx4(e(:,:,OP_DR),ri3_79,f(:,OP_1),g(:,OP_GS))
  end if

  v3hvpsi = temp
end function v3hvpsi

function v3hvb(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hvb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DP),f(:,OP_DR))  &
       - intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DP),f(:,OP_DZ))  &
       + 2.*intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DP),f(:,OP_1))

  v3hvb = temp
#else
  v3hvb = 0.
#endif
end function v3hvb

function v3hvf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hvf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp79a = f(:,OP_DR) + 2.*ri_79*f(:,OP_1)
  temp = -intx4(e(:,:,OP_DRZ),temp79a,ri_79 ,g(:,OP_DR))   &
       +intx4(e(:,:,OP_DRR),temp79a,ri_79 ,g(:,OP_DZ))   &
       -intx4(e(:,:,OP_DZ),temp79a,ri_79 ,g(:,OP_DRR))   &
       +intx4(e(:,:,OP_DR),temp79a,ri_79 ,g(:,OP_DRZ))   &
       +intx4(e(:,:,OP_DZ),temp79a,ri2_79,g(:,OP_DR))   &
       -intx4(e(:,:,OP_DR),temp79a,ri2_79,g(:,OP_DZ))   &      
       -intx3(e(:,:,OP_DZZ),f(:,OP_DZ),g(:,OP_DR))  &
       +intx3(e(:,:,OP_DRZ),f(:,OP_DZ),g(:,OP_DZ))  &
       -intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_DRZ))  &
       +intx3(e(:,:,OP_DR),f(:,OP_DZ),g(:,OP_DZZ))  

  v3hvf = temp
#else
  v3hvf = 0.
#endif
end function v3hvf

function v3hchipsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hchipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = - intx4(e(:,:,OP_DR),ri6_79,f(:,OP_GSP),g(:,OP_DR )) &
       -   intx4(e(:,:,OP_DZ),ri6_79,f(:,OP_GSP),g(:,OP_DZ )) &
       -   intx4(e(:,:,OP_DZ),ri6_79,g(:,OP_GS) ,f(:,OP_DZP)) &
       -   intx4(e(:,:,OP_DR),ri6_79,g(:,OP_GS) ,f(:,OP_DRP)) &
       - 2*intx4(e(:,:,OP_DZ),ri7_79,f(:,OP_DZ),g(:,OP_DRP))  &
       + 2*intx4(e(:,:,OP_DR),ri7_79,f(:,OP_DZ),g(:,OP_DZP))
     
  v3hchipsi = temp
#else
  v3hchipsi = 0.
#endif
end function v3hchipsi

function v3hchib(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hchib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp = -intx4(e(:,:,OP_DZ),ri7_79,g(:,OP_DP),f(:,OP_DRP)) &
       +  intx4(e(:,:,OP_DR),ri7_79,g(:,OP_DP),f(:,OP_DZP)) 
#endif
  if(itor.eq.1) then
     temp = temp                                             &
          -2.*intx4(e(:,:,OP_LP),ri6_79,g(:,OP_1),f(:,OP_DZ)) &
          +8.*intx4(e(:,:,OP_DR),ri7_79,g(:,OP_1),f(:,OP_DZ))
  endif

  v3hchib = temp
end function v3hchib

function v3hchif(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3hchif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = 2.*intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DZ),g(:,OP_DRP))  &
       +2.*intx4(e(:,:,OP_DZ),ri6_79,f(:,OP_DZ),g(:,OP_DZP))  &
       +intx4(e(:,:,OP_DRZ),ri5_79,f(:,OP_DRP),g(:,OP_DR))  &
       -intx4(e(:,:,OP_DRR),ri5_79,f(:,OP_DRP),g(:,OP_DZ))  &
       +intx4(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP),g(:,OP_DRR))  &
       -intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DRP),g(:,OP_DRZ))  &
       -intx4(e(:,:,OP_DZ),ri6_79,f(:,OP_DRP),g(:,OP_DR))  &
       +intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DRP),g(:,OP_DZ))  &
       +intx4(e(:,:,OP_DZZ),ri5_79,f(:,OP_DZP),g(:,OP_DR))  &
       -intx4(e(:,:,OP_DRZ),ri5_79,f(:,OP_DZP),g(:,OP_DZ))  &
       +intx4(e(:,:,OP_DZ),ri5_79,f(:,OP_DZP),g(:,OP_DRZ))  &
       -intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DZP),g(:,OP_DZZ))  

  v3hchif = temp
#else
  v3hchif = 0.
#endif

end function v3hchif

function b1vphi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1vphi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),r_79,f(:,OP_DZ),g(:,OP_DRR))   &
       + intx4(e(:,:,OP_DZ),r_79,f(:,OP_DZ),g(:,OP_DRZ))   &
       + intx4(e(:,:,OP_DR),r_79,g(:,OP_DR),f(:,OP_DRZ))   &
       + intx4(e(:,:,OP_DZ),r_79,g(:,OP_DR),f(:,OP_DZZ))   &
       - intx4(e(:,:,OP_DR),r_79,f(:,OP_DR),g(:,OP_DRZ))   &
       - intx4(e(:,:,OP_DZ),r_79,f(:,OP_DR),g(:,OP_DZZ))   &
       - intx4(e(:,:,OP_DR),r_79,g(:,OP_DZ),f(:,OP_DRR))   &
       - intx4(e(:,:,OP_DZ),r_79,g(:,OP_DZ),f(:,OP_DRZ))   
  if(itor.eq.1) then
     temp = temp + 3.*intx3(e(:,:,OP_DR),f(:,OP_DZ),g(:,OP_DR)) &
          - 3.*intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_DZ)) &
          - 2.*intx3(e(:,:,OP_DR),f(:,OP_1),g(:,OP_DRZ)) &
          - 2.*intx3(e(:,:,OP_DZ),f(:,OP_1),g(:,OP_DZZ)) &
          - 2.*intx3(e(:,:,OP_DR),g(:,OP_DZ),f(:,OP_DR)) &
          - 2.*intx3(e(:,:,OP_DZ),g(:,OP_DZ),f(:,OP_DZ)) &
          - 4.*intx4(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZ))
  endif

  b1vphi = temp
end function b1vphi

function b1vchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1vchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRR),g(:,OP_DR))  &
       +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_DRR))  &
       +intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DRZ),g(:,OP_DR))  &
       +intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DR),g(:,OP_DRZ))  &
       +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRZ),g(:,OP_DZ))  &
       +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DZ),g(:,OP_DRZ))  &
       +intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZZ),g(:,OP_DZ))  &
       +intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_DZZ))  

  if(itor.eq.1) then
     temp = temp + 2.*intx4(e(:,:,OP_DR),ri3_79,f(:,OP_1),g(:,OP_DRR))  &
          + 2.*intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_1),g(:,OP_DRZ))  &
          + 2.*intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DR),f(:,OP_DR))  &
          + 2.*intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DR),f(:,OP_DZ))  &
          - 2.*intx4(e(:,:,OP_DR),ri4_79,f(:,OP_1),g(:,OP_DR))

  endif

  b1vchi = temp
end function b1vchi

function b1vmun(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1vmun
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_LP),g(:,OP_1),h(:,OP_1),f(:,OP_LP))

  if(itor.eq.1) then

  endif

  b1vmun = temp
end function b1vmun

vectype function b2uu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  temp = -int4(r_79,f(:,OP_LP),e(:,OP_DZ),g(:,OP_DR))  &
       +  int4(r_79,f(:,OP_LP),e(:,OP_DR),g(:,OP_DZ))

  b2uu = temp
end function b2uu

vectype function b2uchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  temp = -int4(ri2_79,f(:,OP_LP),e(:,OP_DR),g(:,OP_DR))  &
       -  int4(ri2_79,f(:,OP_LP),e(:,OP_DZ),g(:,OP_DZ))

  b2uchi = temp
end function b2uchi

vectype function b2phimun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = -int4(g(:,OP_1),f(:,OP_LP),e(:,OP_DR),h(:,OP_DR)) &
       -  int4(g(:,OP_1),f(:,OP_LP),e(:,OP_DZ),h(:,OP_DZ)) &
       -  int4(g(:,OP_1),f(:,OP_LP),e(:,OP_LP),h(:,OP_1)) 

  if(itor.eq.1) then
     
  endif

  b2phimun = temp
end function b2phimun

vectype function b2chimun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = -2.*int5(ri3_79,g(:,OP_1),f(:,OP_LP),h(:,OP_DZ),e(:,OP_DR)) &
       +  2.*int5(ri3_79,g(:,OP_1),f(:,OP_LP),h(:,OP_DR),e(:,OP_DZ)) 

  if(itor.eq.1) then
  endif

  b2chimun = temp
end function b2chimun

end module two_fluid

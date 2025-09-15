module metricterms_new


implicit none

contains

!============================================================================
! V1 TERMS
!============================================================================
  

! V1umu 
! =====
function v1umu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1umu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

     temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

     temp79d = 2.*f(:,OP_DRZ)
     if(itor.eq.1) then
        temp79d = temp79d + ri_79*f(:,OP_DZ)
     endif

     temp = &
          - intx4(e(:,:,OP_DZZ),r2_79,temp79b,g(:,OP_1)) &
          + intx4(e(:,:,OP_DRR),r2_79,temp79b,g(:,OP_1)) &
          - 2.*intx4(e(:,:,OP_DRZ),r2_79,temp79d,g(:,OP_1))
        
     if(itor.eq.1) then
        temp = temp &
             + intx4(e(:,:,OP_DR),r_79,temp79b,g(:,OP_1)) &
             - intx4(e(:,:,OP_DZ),r_79,temp79d,g(:,OP_1)) &
             + 5.*intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
             - 8.*intx3(e(:,:,OP_DZ),f(:,OP_DZ),h(:,OP_1))
     endif
     
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp &
#ifdef USEST
          - intx3(e(:,:,OP_DZP),f(:,OP_DZP),g(:,OP_1)) &
          - intx3(e(:,:,OP_DRP),f(:,OP_DRP),g(:,OP_1)) &
          - intx3(e(:,:,OP_DZ),f(:,OP_DZP),g(:,OP_DP)) &
          - intx3(e(:,:,OP_DR),f(:,OP_DRP),g(:,OP_DP))
#else
          + intx3(e(:,:,OP_DZ),f(:,OP_DZPP),g(:,OP_1)) &
          + intx3(e(:,:,OP_DR),f(:,OP_DRPP),g(:,OP_1))
#endif
#endif

  v1umu = temp
end function v1umu


! V1vmu
! =====
function v1vmu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vmu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx4(e(:,:,OP_DR),r_79,f(:,OP_DZP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),r_79,f(:,OP_DRP),g(:,OP_1))
        
  if(itor.eq.1) then
     temp = temp &
          + 4.*intx3(e(:,:,OP_DZ),f(:,OP_DP),h(:,OP_1)) &
          - 2.*intx3(e(:,:,OP_DZ),f(:,OP_DP),g(:,OP_1))
  endif
#endif
  v1vmu = temp
end function v1vmu


! V1chimu
! =======
function v1chimu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chimu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
  temp79c = f(:,OP_DRZ)
  if(itor.eq.1) then
     temp79b = temp79b + 2.*ri_79*f(:,OP_DR)
     temp79c = temp79c -    ri_79*f(:,OP_DZ)
  endif
     
  temp = -2.* &
       (intx4(e(:,:,OP_DRZ),ri_79,temp79b,g(:,OP_1)) &
       -intx4(e(:,:,OP_DZZ),ri_79,temp79c,g(:,OP_1)) &
       +intx4(e(:,:,OP_DRR),ri_79,temp79c,g(:,OP_1)))
        
  if(itor.eq.1) then
     temp = temp &
          -    intx4(e(:,:,OP_DZ),ri2_79,temp79b,g(:,OP_1)) &
          - 2.*intx4(e(:,:,OP_DR),ri2_79,temp79c,g(:,OP_1)) &
          + 4.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),h(:,OP_1)) &
          - 3.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),g(:,OP_1))
  endif
     
#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp &
#ifdef USEST
       - intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZP),g(:,OP_1)) &
       + intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DRP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_DP)) &
       + intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_DP))
#else
       + intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZPP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRPP),g(:,OP_1))
#endif
#endif

  v1chimu = temp
end function v1chimu




! V1un
! ====
function v1un(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1un
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),r2_79,f(:,OP_DR),g(:,OP_1)) &
       + intx4(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),g(:,OP_1))

  v1un = temp
end function v1un


! V1chin
! ======
function v1chin(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1))

  v1chin = temp
end function v1chin


! V1psipsi
! ========
function v1psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1psipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_GS)) &
       - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_GS))

  v1psipsi = temp
end function v1psipsi


! V1psib
! ======
function v1psib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1psib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
       + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1))

  v1psib = temp
#else
  v1psib = 0.
#endif
end function v1psib


! V1bb
! ====
function v1bb(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1bb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  v1bb = 0.
end function v1bb


! V1uun 
! =====
function v1uun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1uun
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1uun = 0.
     return
  end if

  temp = -intx5(e(:,:,OP_DZ),r3_79,f(:,OP_DR),g(:,OP_LP),h(:,OP_1)) &
       +  intx5(e(:,:,OP_DR),r3_79,f(:,OP_DZ),g(:,OP_LP),h(:,OP_1)) &
       -  intx5(e(:,:,OP_DR),r3_79,f(:,OP_DRZ),g(:,OP_DR),h(:,OP_1)) &
       -  intx5(e(:,:,OP_DR),r3_79,f(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1)) &
       +  intx5(e(:,:,OP_DZ),r3_79,f(:,OP_DRR),g(:,OP_DR),h(:,OP_1)) &
       +  intx5(e(:,:,OP_DZ),r3_79,f(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1))

  if(itor.eq.1) then
     temp = temp &
          + intx5(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          + intx5(e(:,:,OP_DZ),r2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  end if

  v1uun = temp
end function v1uun


! V1uvn
! =====
function v1uvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1uvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1uvn = 0.
     return
  end if

  temp = &
       - intx5(e(:,:,OP_DZ),r2_79,f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
       - intx5(e(:,:,OP_DR),r2_79,f(:,OP_DRP),g(:,OP_1),h(:,OP_1))
#else
  temp = 0.
#endif
  v1uvn = temp
end function v1uvn


! v1uchin
! =======
function v1uchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1uchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1uchin = 0.
     return
  end if


  temp79a = (f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &
       +     f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ))
  temp79b = (f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &
       +     f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ))
  temp = -intx4(e(:,:,OP_DR),f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) &
       -  intx4(e(:,:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
       -  intx3(e(:,:,OP_DZ),temp79a,h(:,OP_1)) &
       +  intx3(e(:,:,OP_DR),temp79b,h(:,OP_1))

  if(itor.eq.1) then
     temp = temp &
          - 2.*intx5(e(:,:,OP_DR),ri_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))&
          - 2.*intx5(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1))&
          + intx5(e(:,:,OP_DZ),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
          - intx5(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
  end if

  v1uchin = temp
end function v1uchin


! V1vvn
! =====
function v1vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1vvn = 0.
     return
  end if

  temp = 0.
  if(itor.eq.1) then
     temp = -intx5(e(:,:,OP_DZ),r2_79,f(:,OP_1),g(:,OP_1),h(:,OP_1))
  endif

  v1vvn = temp
end function v1vvn


! V1vchin
! =======
function v1vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1vchin = 0.
     return
  end if

  temp = intx5(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1)) &
       - intx5(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1))
#else
  temp = 0.
#endif
  v1vchin = temp
end function v1vchin


! v1chichin
! =========
function v1chichin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chichin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1chichin = 0.
     return
  end if

  temp = -intx5(e(:,:,OP_DR),ri3_79,f(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1)) &
       -  intx5(e(:,:,OP_DR),ri3_79,f(:,OP_DRZ),g(:,OP_DR),h(:,OP_1)) &
       +  intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1)) &
       -  intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_DRR),g(:,OP_DR),h(:,OP_1))

  if(itor.eq.1) then
     temp = temp + 2.*&
          (intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) & 
          +intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)))
  end if
  
  v1chichin = temp
end function v1chichin



! V1upsipsi
! =========
function v1upsipsi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1upsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp
  
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe, tempf
  integer :: j
  

  ! |u, psi(1)|,r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
  ! |u, psi(1)|,z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)

  do j=1, dofs_per_element
     ! |nu, psi(2)|,r
     tempe(j,:) = e(j,:,OP_DRZ)*h(:,OP_DR ) - e(j,:,OP_DRR)*h(:,OP_DZ ) &
          +       e(j,:,OP_DZ )*h(:,OP_DRR) - e(j,:,OP_DR )*h(:,OP_DRZ)
     ! |nu, psi(2)|,z
     tempf(j,:) = e(j,:,OP_DZZ)*h(:,OP_DR ) - e(j,:,OP_DRZ)*h(:,OP_DZ ) &
          +       e(j,:,OP_DZ )*h(:,OP_DRZ) - e(j,:,OP_DR )*h(:,OP_DZZ)
  end do

  temp = -intx2(tempe,temp79b)  &
       - intx2(tempf,temp79c)  &
       + intx3(e(:,:,OP_DZ),temp79b,h(:,OP_GS))  &
       - intx3(e(:,:,OP_DR),temp79c,h(:,OP_GS))

  if(itor.eq.1) then
     ! |u, psi(1)|
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           
     do j=1, dofs_per_element
        ! |nu,psi(2)|
        tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
     end do
     temp = temp             &
          - intx3(tempd,ri_79,temp79b)   &
          - intx3(tempe,ri_79,temp79a)   &
          - intx3(tempd,ri2_79,temp79a)   &
          + intx4(e(:,:,OP_DZ),ri_79,h(:,OP_GS),temp79a)
  endif

  v1upsipsi = temp
end function v1upsipsi


! V1upsib
! =======
function v1upsib(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1upsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempc, tempd
  integer :: j

  do j=1, dofs_per_element
     tempa(j,:) = h(:,OP_1)*e(j,:,OP_GS)  &
          + h(:,OP_DZ)*e(j,:,OP_DZ) + h(:,OP_DR)*e(j,:,OP_DR)
     tempc(j,:) = e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP)
     tempd(j,:) = h(:,OP_DP)* &
          (e(j,:,OP_DZ)*f(:,OP_DR )-e(j,:,OP_DR)*f(:,OP_DZ )) &
          +    h(:,OP_1 )* &
          (e(j,:,OP_DZ)*f(:,OP_DRP)-e(j,:,OP_DR)*f(:,OP_DZP))
  end do
  temp79b = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)
#ifdef USEST
  temp79e = h(:,OP_1)*f(:,OP_GS) &
       + h(:,OP_DZ)*f(:,OP_DZ ) + h(:,OP_DR)*f(:,OP_DR ) 
#else
  temp79e = h(:,OP_DP)*f(:,OP_GS) + h(:,OP_1)*f(:,OP_GSP) &
       +    h(:,OP_DZP)*f(:,OP_DZ ) + h(:,OP_DRP)*f(:,OP_DR ) &
       +    h(:,OP_DZ )*f(:,OP_DZP) + h(:,OP_DR )*f(:,OP_DRP)
#endif
     
  temp = intx3(tempa,ri_79,temp79b) &
       + intx4(tempc,ri_79,f(:,OP_DR),h(:,OP_DZ)) &
       - intx4(tempc,ri_79,f(:,OP_DZ),h(:,OP_DR)) &
       - intx3(tempd,ri_79,g(:,OP_GS)) &
#ifdef USEST
       + intx4(e(:,:,OP_DZP),ri_79,g(:,OP_DR),temp79e) &
       - intx4(e(:,:,OP_DRP),ri_79,g(:,OP_DZ),temp79e) &
       + intx4(e(:,:,OP_DZ),ri_79,g(:,OP_DRP),temp79e) &
       - intx4(e(:,:,OP_DR),ri_79,g(:,OP_DZP),temp79e)
#else
       - intx4(e(:,:,OP_DZ),ri_79,g(:,OP_DR),temp79e) &
       + intx4(e(:,:,OP_DR),ri_79,g(:,OP_DZ),temp79e)
#endif
  temp = -temp

  v1upsib = temp
#else
  v1upsib = 0.
#endif

end function v1upsib


! V1ubb 
! =====
function v1ubb(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ubb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
  integer :: j

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST 
  do j=1, dofs_per_element
     tempa(j,:) = &
          -(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
          *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) &
          -(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
          *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP))
  end do
  temp = intx2(tempa,ri2_79)
        !temp = intx5(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1),ri2_79,h(:,OP_DPP))&
        !     + intx5(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1),ri2_79,h(:,OP_DPP))
#else
  do j=1, dofs_per_element
     tempa(j,:) = &
          (e(j,:,OP_DZ)*f(:,OP_DZPP) + e(j,:,OP_DR)*f(:,OP_DRPP)) &
          *g(:,OP_1) &
          + 2.*(e(j,:,OP_DZ)*f(:,OP_DZP) + e(j,:,OP_DR)*f(:,OP_DRP)) &
          *g(:,OP_DP) &
          +    (e(j,:,OP_DZ)*f(:,OP_DZ) + e(j,:,OP_DR)*f(:,OP_DR)) &
          *g(:,OP_DPP)
  end do
  temp = intx3(tempa,ri2_79,h(:,OP_1))
#endif
#endif

  v1ubb = temp
end function v1ubb

#ifdef USE3D
! V1upsif
! =====
function v1upsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1upsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe
  integer :: j

  ! [u, psi]_R*R
  temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
       + f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)
  if(itor.eq.1) then
     temp79a = temp79a - ri_79* &
          (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
  end if
  ! [u, psi]_Z*R
  temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
          + f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ) 
  ! (u, f')_R
  temp79c = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
          + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
  ! (u, f')_Z
  temp79d = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
          + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
  if(itor.eq.1) then
     ! (u, f')
     temp79e = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ)           
  endif

  do j=1, dofs_per_element
     ! (nu, f')_R
     tempa(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) &
          + e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) 
     ! (nu, f')_Z
     tempb(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) &
          + e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) 
     ! [nu, psi]_R*R
     tempc(j,:) = &
            e(j,:,OP_DZ)*g(:,OP_DRR) - e(j,:,OP_DR)*g(:,OP_DRZ) &
          + e(j,:,OP_DRZ)*g(:,OP_DR) - e(j,:,OP_DRR)*g(:,OP_DZ)
     if(itor.eq.1) then
        tempc(j,:) = tempc(j,:) - ri_79* &
             (e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)) 
     end if
     ! [nu, psi]_Z*R
     tempd(j,:) = &
            e(j,:,OP_DZ)*g(:,OP_DRZ) - e(j,:,OP_DR)*g(:,OP_DZZ) &
          + e(j,:,OP_DZZ)*g(:,OP_DR) - e(j,:,OP_DRZ)*g(:,OP_DZ) 
     ! [nu, R^2*(U, f')]/R
     tempe(j,:) = e(j,:,OP_DZ)*temp79c - e(j,:,OP_DR)*temp79d
     if(itor.eq.1) then
        tempe(j,:) = tempe(j,:) + 2*e(j,:,OP_DZ)*temp79e*ri_79 
     endif
  end do

  temp = intx3(tempa,temp79a,r_79) &
       + intx3(tempb,temp79b,r_79) &
       + intx3(tempc,temp79c,r_79) &
       + intx3(tempd,temp79d,r_79) &
       - intx3(tempe,g(:,OP_GS),r_79) 

  v1upsif = temp
end function v1upsif

! V1ubf
! =====
function v1ubf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ubf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j


  do j=1, dofs_per_element
     ! (nu, f')'
     tempa(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) &
          + e(j,:,OP_DZP)*h(:,OP_DZ) + e(j,:,OP_DRP)*h(:,OP_DR) 
     ! [nu,f'']*R
     tempb(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
  end do
  ! F*u_LP + R^2*(u, F/R^2)
  temp79a = g(:,OP_1)*f(:,OP_LP)  & 
       + g(:,OP_DR)*f(:,OP_DR) + g(:,OP_DZ)*f(:,OP_DZ)
  if(itor.eq.1) then
     temp79a = temp79a - 2*ri_79*g(:,OP_1)*f(:,OP_DR)
  end if
  ! (R^2(u, f'))_R/R^2
  temp79b = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
       + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79b = temp79b + 2*ri_79* &
          (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
  end if
  ! (R^2(u, f'))_Z/R^2
  temp79c = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
       +    f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
  ![u,F]*R
  temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ) 

  temp = intx2(tempa,temp79a) &
       + intx3(e(:,:,OP_DRP),temp79b,g(:,OP_1)) &
       + intx3(e(:,:,OP_DZP),temp79c,g(:,OP_1)) &
       + intx3(e(:,:,OP_DR),temp79b,g(:,OP_DP)) &
       + intx3(e(:,:,OP_DZ),temp79c,g(:,OP_DP)) &
       + intx2(tempb,temp79d)

  v1ubf = temp
end function v1ubf

! V1uff
! =====
function v1uff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1uff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     !(nu, f')_R
     tempa(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) &
          + e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) 
     !(nu, f')_Z
     tempb(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) &
          + e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) 
  end do
  !(u, f')_R
  temp79a = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &           
       +    f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ)            
  !(u, f')_Z
  temp79b = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &           
       +    f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ)            

  temp = -intx3(tempa,temp79a,r2_79) &
       -  intx3(tempb,temp79b,r2_79) 

  v1uff = temp
end function v1uff
#endif

! V1up
! ====
function v1up(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1up
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     temp = 2.* &
          (intx4(e(:,:,OP_DZ),r_79,f(:,OP_DR),g(:,OP_DZ)) &
          -intx4(e(:,:,OP_DZ),r_79,f(:,OP_DZ),g(:,OP_DR))) &
          -4.*gam*intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
  end if

  v1up = temp
end function v1up



! V1vpsipsi
! =========
function v1vpsipsi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vpsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempc
  integer :: j

  do j=1, dofs_per_element
     tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP)
     tempc(j,:) = &
          f(:,OP_DP)* &
          (e(j,:,OP_DZ)*h(:,OP_DR ) - e(j,:,OP_DR)*h(:,OP_DZ )) &
          +  f(:,OP_1 )* &
          (e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP))
  end do
#ifdef USEST
  temp79b = f(:,OP_1)*g(:,OP_GS) &
       + f(:,OP_DZ )*g(:,OP_DZ) + f(:,OP_DR )*g(:,OP_DR)
#else
  temp79b = f(:,OP_DP)*g(:,OP_GS) + f(:,OP_1)*g(:,OP_GSP) &
       + f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
       + f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)
#endif
  temp = intx4(tempa,ri_79,g(:,OP_DZ),f(:,OP_DR)) &
       - intx4(tempa,ri_79,g(:,OP_DR),f(:,OP_DZ)) &
#ifdef USEST
       - intx4(e(:,:,OP_DZP),ri_79,h(:,OP_DR),temp79b) &
       + intx4(e(:,:,OP_DRP),ri_79,h(:,OP_DZ),temp79b) &
       - intx4(e(:,:,OP_DZ),ri_79,h(:,OP_DRP),temp79b) &
       + intx4(e(:,:,OP_DR),ri_79,h(:,OP_DZP),temp79b) &
#else
       + intx4(e(:,:,OP_DZ),ri_79,h(:,OP_DR),temp79b) &
       - intx4(e(:,:,OP_DR),ri_79,h(:,OP_DZ),temp79b) &
#endif
       + intx3(tempc,ri_79,g(:,OP_GS))

  temp= -temp

  v1vpsipsi = temp
#else
  v1vpsipsi = 0.
#endif
end function v1vpsipsi



! V1vpsib
! =======
function v1vpsib(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vpsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa
  integer :: j

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST 
        !temp = - intx5(e(:,:,OP_DR),f(:,OP_1),g(:,OP_DR),ri2_79,h(:,OP_DPP))&
        !       - intx5(e(:,:,OP_DZ),f(:,OP_1),g(:,OP_DZ),ri2_79,h(:,OP_DPP))
  do j=1, dofs_per_element
     tempa(j,:) = &
          +(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
          *(g(:,OP_DZP)*f(:,OP_1) + g(:,OP_DZ)*f(:,OP_DP)) &
          +(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
          *(g(:,OP_DRP)*f(:,OP_1) + g(:,OP_DR)*f(:,OP_DP))
  end do
  temp = intx2(tempa,ri2_79)
#else
  do j=1, dofs_per_element
     tempa(j,:) = &
          f(:,OP_DPP)* &
          (e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR)) &
          +2.*f(:,OP_DP)* &
          (e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP)) &
          +   f(:,OP_1)* &
          (e(j,:,OP_DZ)*g(:,OP_DZPP) + e(j,:,OP_DR)*g(:,OP_DRPP))
  end do
  temp = -intx3(tempa,ri2_79,h(:,OP_1))
#endif
#endif

  v1vpsib = temp
end function v1vpsib

#ifdef USE3D
! V1vpsif
! =====
function v1vpsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vpsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe, tempf
  integer :: j


  do j=1, dofs_per_element
     ! [nu, psi]*R
     tempa(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
     ! [nu, f'']*R
     tempb(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
     ! (nu, psi')
     tempc(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP) 
     ! (nu, f')
     tempd(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
     ! (nu, f'')
     tempe(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
     ! (nu', f')
     tempf(j,:) = &
          e(j,:,OP_DZP)*h(:,OP_DZ) + e(j,:,OP_DRP)*h(:,OP_DR) 
  end do

  ! [v, f']'*R
  temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
  ! [v, psi]*R
  temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)            
  ! (v, f') + v*f'_LP
  temp79c = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
       + f(:,OP_1)*h(:,OP_LP)
  ! (v, psi)
  temp79d = f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ) 
  
  temp = - intx2(tempa,temp79a) &
       + intx2(tempb,temp79b) &
       - intx2(tempc,temp79c) &
       - intx2(tempe,temp79d) & 
       - intx2(tempf,temp79d) & 
       + intx3(tempd,g(:,OP_GS),f(:,OP_DP)) &
       - intx3(tempf,g(:,OP_GS),f(:,OP_1)) 

  v1vpsif = temp
end function v1vpsif

! V1vbf
! =====
function v1vbf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vbf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
  integer :: j


  do j=1, dofs_per_element
     ! [nu, f']*R
     tempa(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ) 
     ! [nu, f'']*R
     tempb(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
     ! [nu', f'']*R
     tempc(j,:) = &
          e(j,:,OP_DZP)*h(:,OP_DRP) - e(j,:,OP_DRP)*h(:,OP_DZP) 
  end do

  temp = - intx4(tempa,g(:,OP_1),f(:,OP_DPP),ri_79) &
       - intx4(tempb,g(:,OP_1),f(:,OP_DP),ri_79) &
       + intx4(tempc,g(:,OP_1),f(:,OP_1),ri_79) 

  v1vbf = temp
end function v1vbf

! V1vff
! =====
function v1vff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     ! (nu, f')
     tempa(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR) 
     ! [nu, f'']*R
     tempb(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DRP) - e(j,:,OP_DR)*g(:,OP_DZP) 
  end do
  ! [v, f']'*R 
  temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
       + f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
  ! (v, f')
  temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) 

  temp = intx3(tempa,temp79a,r_79) &
       - intx3(tempb,temp79b,r_79) 

  v1vff = temp
end function v1vff
#endif


! V1vp
! ====
function v1vp(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vp
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  if(itor.eq.1) then
     temp79a =  f(:,OP_1 )*g(:,OP_DP) &
          + gam*f(:,OP_DP)*g(:,OP_1 )
     temp = 2.*intx2(e(:,:,OP_DZ),temp79a)
  end if
#endif

  v1vp = temp
end function v1vp



! V1chipsipsi
! ===========
function v1chipsipsi(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempd, tempe, tempf
  integer :: j

  temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
  temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)
        
  do j=1, dofs_per_element
     tempe(j,:) = e(j,:,OP_DRZ)*h(:,OP_DR ) - e(j,:,OP_DRR)*h(:,OP_DZ ) &
          +       e(j,:,OP_DZ )*h(:,OP_DRR) - e(j,:,OP_DR )*h(:,OP_DRZ)
     tempf(j,:) = e(j,:,OP_DZZ)*h(:,OP_DR ) - e(j,:,OP_DRZ)*h(:,OP_DZ ) &
          +       e(j,:,OP_DZ )*h(:,OP_DRZ) - e(j,:,OP_DR )*h(:,OP_DZZ)
  end do

  temp = intx3(tempe,ri3_79,temp79b) &
       + intx3(tempf,ri3_79,temp79c) &
       + intx4(e(:,:,OP_DR),ri3_79,temp79c,h(:,OP_GS)) &
       - intx4(e(:,:,OP_DZ),ri3_79,temp79b,h(:,OP_GS))
     
  if(itor.eq.1) then
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
     do j=1, dofs_per_element
        tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
     end do
           
     temp = temp + &
          (2.*intx4(e(:,:,OP_DZ),ri4_79,temp79a,h(:,OP_GS)) &
          -2.*intx3(tempe,ri4_79,temp79a) &
          +   intx3(tempd,ri4_79,temp79b) &
          -2.*intx3(tempd,ri5_79,temp79a))
  endif

  v1chipsipsi = temp
end function v1chipsipsi


! V1chipsib
! =========
function v1chipsib(e,f,g,h) 
  use basic
  use arrays
  use m3dc1_nint

  implicit none
  
  vectype, dimension(dofs_per_element) :: v1chipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  temp79a = h(:,OP_DZP)*f(:,OP_DR ) - h(:,OP_DRP)*f(:,OP_DZ ) &
       +    h(:,OP_DZ )*f(:,OP_DRP) - h(:,OP_DR )*f(:,OP_DZP)
  do j=1, dofs_per_element
     tempb(j,:) = (e(j,:,OP_DZ)*f(:,OP_DZ )+e(j,:,OP_DR)*f(:,OP_DR )) &
          *h(:,OP_DP) &
          +    (e(j,:,OP_DZ)*f(:,OP_DZP)+e(j,:,OP_DR)*f(:,OP_DRP)) &
          *h(:,OP_1 )
  end do
  temp79c = f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)
  temp79d = h(:,OP_1)*f(:,OP_GS) &
       + h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
  if(itor.eq.1) then
     temp79a = temp79a + 4.*ri_79* &
          (h(:,OP_DP)*f(:,OP_DZ) + h(:,OP_1)*f(:,OP_DZP))
     temp79d = temp79d - 2.*ri_79*h(:,OP_1)*f(:,OP_DR)
  endif
     
  temp = intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DR),temp79a)  &
       - intx4(e(:,:,OP_DR),ri4_79,g(:,OP_DZ),temp79a)  &
       - intx3(tempb,ri4_79,g(:,OP_GS))               &
       - intx4(e(:,:,OP_GS),ri4_79,h(:,OP_1 ),temp79c)  &
       - intx4(e(:,:,OP_DZ),ri4_79,h(:,OP_DZ),temp79c)  &
       - intx4(e(:,:,OP_DR),ri4_79,h(:,OP_DR),temp79c)  &
       + intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),temp79d) &
       + intx4(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),temp79d)
  temp = -temp

  v1chipsib = temp
#else
  v1chipsib = 0.
#endif

end function v1chipsib


! V1chibb
! =======
function v1chibb(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chibb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa
  integer :: j

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST 
        !temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),g(:,OP_1),ri5_79,h(:,OP_DPP))&
        !     - intx5(e(:,:,OP_DR),f(:,OP_DZ),g(:,OP_1),ri5_79,h(:,OP_DPP))
  do j=1, dofs_per_element
     tempa(j,:) = &
          (e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
          *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP)) &
          -(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
          *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP))
  end do
  temp = intx2(tempa,ri5_79)
#else
  do j=1, dofs_per_element
     tempa(j,:) = &
          (e(j,:,OP_DZ)*f(:,OP_DR) - e(j,:,OP_DR)*f(:,OP_DZ))*g(:,OP_DPP) &
          + 2.*(e(j,:,OP_DZ)*f(:,OP_DRP) - e(j,:,OP_DR)*f(:,OP_DZP))*g(:,OP_DP) &
          +    (e(j,:,OP_DZ)*f(:,OP_DRPP) - e(j,:,OP_DR)*f(:,OP_DZPP))*g(:,OP_1)
  end do
  temp = -intx3(tempa,ri5_79,h(:,OP_1))
#endif
#endif

  v1chibb = temp
end function v1chibb

#ifdef USE3D
! V1chipsif
! =====
function v1chipsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chipsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempc, tempd
  integer :: j

  do j=1, dofs_per_element
     ! ([nu, psi]*R^2)_R / R
     tempa(j,:) = &
            e(j,:,OP_DRZ)*g(:,OP_DR) - e(j,:,OP_DRR)*g(:,OP_DZ) &
          + e(j,:,OP_DZ)*g(:,OP_DRR) - e(j,:,OP_DR)*g(:,OP_DRZ)
     if(itor.eq.1) then 
        tempa(j,:) = tempa(j,:) + ri_79* &
             (e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)) 
     end if
     ! ([nu, psi]*R^2)_Z / R
     tempb(j,:) = &
            e(j,:,OP_DZZ)*g(:,OP_DR) - e(j,:,OP_DRZ)*g(:,OP_DZ) &
          + e(j,:,OP_DZ)*g(:,OP_DRZ) - e(j,:,OP_DR)*g(:,OP_DZZ)
     ! (R^2*(nu, f'))_R / R^2
     tempc(j,:) = &
            e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) & 
          + e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) 
     if(itor.eq.1) then 
        tempc(j,:) = tempc(j,:) + 2*ri_79* &
             (e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)) 
     end if
     ! (R^2*(nu, f'))_Z / R^2
     tempd(j,:) = &
            e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) & 
          + e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) 
  end do

  ! [chi, f']_R*R 
  temp79a = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
  if(itor.eq.1) then 
     temp79a = temp79a - ri_79* &        
          (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
  end if
  ! [chi, f']_Z*R 
  temp79b = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
  ! ((chi, psi)/R^2)_R * R^2
  temp79c = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &
       +    f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ) 
  if(itor.eq.1) then 
     temp79c = temp79c - 2*ri_79*  & 
          (f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ)) 
  end if
        ! ((chi, psi)/R^2)_Z * R^2
  temp79d = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &
       +    f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ) 

  temp = intx3(tempa,temp79a,ri2_79) &
       + intx3(tempb,temp79b,ri2_79) &
       - intx3(tempc,temp79c,ri2_79) &
       - intx3(tempd,temp79d,ri2_79) &
       - intx4(e(:,:,OP_DZ),temp79a,ri2_79,g(:,OP_GS)) &
       + intx4(e(:,:,OP_DR),temp79b,ri2_79,g(:,OP_GS)) 

  v1chipsif = temp
end function v1chipsif

! V1chibf
! =====
function v1chibf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chibf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     ! (nu, f')
     tempa(j,:) = &
          + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
     ! [nu, f'']*R
     tempb(j,:) = &
          + e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
  end do

  ! [chi, F/R^4]'*R^5
  temp79a = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP) 
  if(itor.eq.1) then 
     temp79a = temp79a - 4*ri_79* &        
          (f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) 
  end if
  ! [chi, f']_R*R 
  temp79b = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
  if(itor.eq.1) then 
     temp79b = temp79b - ri_79* &        
          (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
  end if
  ! [chi, f']_Z*R 
  temp79c = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
  ! (chi, F/R^2)*R^2 + F*chi_GS
  temp79d = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
       + f(:,OP_GS)*g(:,OP_1) 
  if(itor.eq.1) then 
     temp79d = temp79d - 2*ri_79*f(:,OP_DR)*g(:,OP_1) 
  end if

  temp = - intx3(tempa,temp79a,ri3_79) &
       - intx3(tempb,temp79d,ri3_79) &
       + intx4(e(:,:,OP_DZ),temp79c,ri3_79,g(:,OP_DP)) &
       + intx4(e(:,:,OP_DR),temp79b,ri3_79,g(:,OP_DP)) &
       + intx4(e(:,:,OP_DZP),temp79c,ri3_79,g(:,OP_1)) &
       + intx4(e(:,:,OP_DRP),temp79b,ri3_79,g(:,OP_1)) 

  v1chibf = temp
end function v1chibf

! V1chiff
! =====
function v1chiff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chiff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     ! (R^2*(nu, f'))_R / R^2
     tempa(j,:) = &
            e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) & 
          + e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) 
     if(itor.eq.1) then 
        tempa(j,:) = tempa(j,:) + 2*ri_79* &
             (e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)) 
     end if
     ! (R^2*(nu, f'))_Z / R^2
     tempb(j,:) = &
            e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) & 
          + e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) 
  end do

  ! [chi, f']_R*R
  temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
  if(itor.eq.1) then 
     temp79a = temp79a - ri_79* &        
          (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
  end if
  ! [chi, f']_Z*R
  temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            

  temp = - intx3(tempa,temp79a,ri_79) &
       -   intx3(tempb,temp79b,ri_79) 

  v1chiff = temp
end function v1chiff
#endif


! V1chip
! ======
function v1chip(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chip
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     temp79a = gam*f(:,OP_GS)*g(:,OP_1) + &
          f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
     
     temp = 2.*intx3(e(:,:,OP_DZ),ri2_79,temp79a)
  end if

  v1chip = temp
end function v1chip


! V1ngrav
! =======
function v1ngrav(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ngrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ngrav = 0.
     return
  endif

  temp = gravz*intx3(e(:,:,OP_1), r_79,f(:,OP_DR)) &
       - gravr*intx3(e(:,:,OP_1),ri_79,f(:,OP_DZ))

  v1ngrav = temp
end function v1ngrav


! V1ungrav
! ========
function v1ungrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ungrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ungrav = 0.
     return
  endif

  temp79a = f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)

  temp = gravz*intx2(e(:,:,OP_DR),temp79a) &
       - gravr*intx3(e(:,:,OP_DZ),ri2_79,temp79a)
     
  if(itor.eq.1) &
       temp = temp + 2.*gravz*intx3(e(:,:,OP_1),ri_79,temp79a)

  v1ungrav = temp
end function v1ungrav


! V1chingrav
! ==========
function v1chingrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chingrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1chingrav = 0.
     return
  endif

  temp79a = r_79*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
       + g(:,OP_1)*f(:,OP_LP))

  temp = gravz*intx2(e(:,:,OP_DR),temp79a) &
       - gravr*intx3(e(:,:,OP_DZ),ri2_79,temp79a)

  if(itor.eq.1) &
       temp = temp + 2.*gravz*intx3(e(:,:,OP_1),ri_79,temp79a)

  v1chingrav = temp
end function v1chingrav


! V1ndenmgrav
! ===========
function v1ndenmgrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ndenmgrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  real, intent(in) :: g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ndenmgrav = 0.
     return
  endif

  temp79a = -g*r_79*f(:,OP_LP)

  temp = gravz*intx2(e(:,:,OP_DR),temp79a) &
       - gravr*intx3(e(:,:,OP_DZ),ri2_79,temp79a)

  if(itor.eq.1) &
       temp = temp + 2.*gravz*intx3(e(:,:,OP_1),ri_79,temp79a)

  v1ndenmgrav = temp
end function v1ndenmgrav


! V1us
! ====
function v1us(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1us
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = -intx4(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),temp79a) &
       -  intx4(e(:,:,OP_DR),r2_79,f(:,OP_DR),temp79a)

  v1us = temp
end function v1us


! V1chis
! ======
function v1chis(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DR),temp79a) &
       - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZ),temp79a)

  v1chis = temp
end function v1chis


! V1psif
! ======
function v1psif(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1psif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = &
       - intx3(e(:,:,OP_DZ),f(:,OP_GS),g(:,OP_DZ)) &
       - intx3(e(:,:,OP_DR),f(:,OP_GS),g(:,OP_DR)) 
#else
  temp = 0.
#endif
  v1psif = temp
end function v1psif


! V1bf
! ====
function v1bf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1bf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = &
       + intx4(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DRP)) &
       - intx4(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZP))
#else
  temp = 0.
#endif
  v1bf = temp
end function v1bf


! V1p
! ===
function v1p(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1p
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  vectype, dimension(dofs_per_element) :: temp

  if(itor.eq.0) then
     v1p = 0.
     return
  end if

  temp = intx3(e(:,:,OP_DZ),r_79,f(:,OP_DR)) &
       - intx3(e(:,:,OP_DR),r_79,f(:,OP_DZ))

  v1p = temp
end function v1p

! V1psiforce
! ===
vectype function v1psiforce(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g
  vectype :: temp

  temp = int3(g(:,OP_1),e(:,OP_DR),f(:,OP_DR)) &
       + int3(g(:,OP_1),e(:,OP_DZ),f(:,OP_DZ))

  v1psiforce = temp
end function v1psiforce

! V1be
! ===

function v1par(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1par
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = - intx2(e(:,:,OP_DZ),f(:,OP_1))

  v1par = temp
end function v1par

function v1parb2ipsipsi(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1parb2ipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i

  vectype, dimension(dofs_per_element) :: temp

  temp79a = f(:,OP_1)*g(:,OP_1)*ri_79
  temp =  intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DR),i(:,OP_DRR))   &
       +  intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DZ),i(:,OP_DRZ))   &
       -  intx4(e(:,:,OP_DR),temp79a,h(:,OP_DR),i(:,OP_DRZ))   &
       -  intx4(e(:,:,OP_DR),temp79a,h(:,OP_DZ),i(:,OP_DZZ))   

  temp79b = -f(:,OP_1)*ri_79
  temp = temp                                                         &
       +  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DR))  &
       +  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ))  &
       -  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DR))  &
       -  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

  temp79c = -g(:,OP_1)*ri_79
  temp = temp                                                         &
       +  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DR)) &
       +  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ)) &
       -  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DR)) &
       -  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

  temp79d = -f(:,OP_1)*g(:,OP_1)*h(:,OP_GS)*ri_79
  temp = temp                                   &
       +  intx3(e(:,:,OP_DZ),temp79d,i(:,OP_DR))   &
       -  intx3(e(:,:,OP_DR),temp79d,i(:,OP_DZ))   

  v1parb2ipsipsi = temp
end function v1parb2ipsipsi

function v1parb2ipsib(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1parb2ipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = - f(:,OP_DP)*g(:,OP_1)*i(:,OP_1)*ri2_79
  temp = +intx3(e(:,:,OP_DR),temp79a,h(:,OP_DR)) &
       +  intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DZ))
#endif

  v1parb2ipsib = temp
end function v1parb2ipsib

#ifdef USEPARTICLES
function v1p_2(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1p_2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        if(inoslip_pol.eq.1 .or. iconst_p.eq.1) then
           temp = 0.
        else
           temp = &
                + intx5(e(:,:,OP_1),r_79,norm79(:,1),f(:,OP_DZ),g) &
                - intx5(e(:,:,OP_1),r_79,norm79(:,2),f(:,OP_DR),g)
        endif
     else
        temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_DR),g) &
             - intx4(e(:,:,OP_DR),r_79,f(:,OP_DZ),g)
     end if

  v1p_2 = temp
end function v1p_2

function v1pbb(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1pbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79b = -g*ri_79
        temp =  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DR)) &
              + intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DZ)) &
              - intx5(e(:,:,OP_DR),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DR)) &
              - intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = - f(:,OP_DP)*bzt79(:,OP_1)*g*ri2_79
        temp = temp +intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DR)) &
              +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ))
        temp79a = g
        temp = temp +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
        temp79a = -g
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DZ)) 
        temp79a = -f(:,OP_DP)*g*bzt79(:,OP_1)*ri_79
        temp = temp +intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR)) &
              -intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ)) 
        temp79a = g*r_79
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              -intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
#endif
    end if

  v1pbb = temp
end function v1pbb

function v1jxb(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1jxb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79a = -f*pst79(:,OP_GS)*ri_79
        temp = temp                                   &
             +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR))   &
             -  intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ))   
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = - f*bzt79(:,OP_1)*ri_79
     temp = temp + intx3(e(:,:,OP_DZ),temp79a,(bzt79(:,OP_DR)+bfpt79(:,OP_DRP))) &
          - intx3(e(:,:,OP_DR),temp79a,bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))
     temp79a = -f*bzt79(:,OP_1)*ri2_79
     temp = temp+intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DRP)) &
          +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZP))
     temp79a = f*pst79(:,OP_GS)
     temp = temp+intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR)) &
          +  intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ))
#endif
    end if

  v1jxb = temp
end function v1jxb

function v1pbb1psi(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1pbb1psi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79b = -g*ri_79
        temp =  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DZ),h(:,OP_DR),pst79(:,OP_DR)) &
              + intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DZ),h(:,OP_DR),pst79(:,OP_DZ)) &
              - intx5(e(:,:,OP_DR),temp79b,f(:,OP_DR),h(:,OP_DZ),pst79(:,OP_DR)) &
              - intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DR),h(:,OP_DZ),pst79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = -g
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,h(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DZ),temp79a,h(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,h(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              +intx5(e(:,:,OP_DR),temp79a,h(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DZ)) 
#endif
    end if

  v1pbb1psi = temp
end function v1pbb1psi

function v1pbb1f(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1pbb1f
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
    if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = g
        temp = temp +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DR),h(:,OP_DR)) &
              +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DZ),h(:,OP_DZ)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DR),h(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DZ),h(:,OP_DZ)) 
        temp79a = g*r_79
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR),f(:,OP_DR),h(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR),f(:,OP_DZ),h(:,OP_DZ)) &
              -intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ),f(:,OP_DR),h(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ),f(:,OP_DZ),h(:,OP_DZ)) 
#endif
    end if

  v1pbb1f = temp
end function v1pbb1f

#endif
!============================================================================
! V2 TERMS
!============================================================================


! V2vn
! ====
function v2vn(e,f,g)

  use basic
  use m3dc1_nint

  implicit none
  
  vectype, dimension(dofs_per_element) :: v2vn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  v2vn = intx4(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_1))
end function v2vn


! V2umu
! =====
function v2umu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2umu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx4(e(:,:,OP_DR),r_79,f(:,OP_DZP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),r_79,f(:,OP_DRP),g(:,OP_1))
  if(itor.eq.1) then
     temp = temp &
          + 2.*intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_1)) &
          - 4.*intx3(e(:,:,OP_1),f(:,OP_DZP),h(:,OP_1))
  endif
#endif
  v2umu = temp
end function v2umu


! V2vmu
! =====
function v2vmu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vmu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = -intx4(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),g(:,OP_1)) &
       -  intx4(e(:,:,OP_DR),r2_79,f(:,OP_DR),g(:,OP_1))

#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp + 2.*intx3(e(:,:,OP_1),f(:,OP_DPP),h(:,OP_1))
#endif

  v2vmu = temp
end function v2vmu


! V2chimu
! =======
function v2chimu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chimu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = h(:,OP_1) - g(:,OP_1)
  temp = &
       - intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1)) &
#ifdef USEST
       - 2.*intx4(e(:,:,OP_DP),ri2_79,f(:,OP_GS),temp79a) &
       - 2.*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GS),h(:,OP_DP)) &
       + 2.*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GS),g(:,OP_DP))
#else
       + 2.*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GSP),temp79a)
#endif
  if(itor.eq.1) then
     temp = temp &
          +2.*intx4(e(:,:,OP_1),ri3_79,f(:,OP_DRP),g(:,OP_1))
  endif

#endif
  v2chimu = temp
end function v2chimu



! V2vun
! =====
function v2vun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vun
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v2vun = 0.
     return
  end if

  temp = intx5(e(:,:,OP_1),r3_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
       - intx5(e(:,:,OP_1),r3_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1))

  if(itor.eq.1) then
     temp = temp + &
          2.*intx5(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_1))
  end if

  v2vun = temp
end function v2vun


! V2vvn
! =====
function v2vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v2vvn = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp = -intx5(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_DP),h(:,OP_1))
#else
  temp = 0.
#endif
  v2vvn = temp
end function v2vvn


! V2up
! ====
function v2up(e,f,g)
  use basic
  use arrays
  use m3dc1_nint
  
  implicit none

  vectype, dimension(dofs_per_element) :: v2up
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx4(e(:,:,OP_1),r_79,f(:,OP_DRP),g(:,OP_DZ)) &
       - intx4(e(:,:,OP_1),r_79,f(:,OP_DZP),g(:,OP_DR)) &
       + intx4(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_DZP)) &
       - intx4(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_DRP))
  if(itor.eq.1) then
     temp = temp - 2.*gam* &
          (intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_1)) &
          +intx3(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DP)))
  endif
#endif

  v2up = temp
end function v2up

! V2vp
! ====
function v2vp(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vp
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp =          intx3(e(:,:,OP_1),g(:,OP_DPP),f(:,OP_1)) &
       + (1.+gam)*intx3(e(:,:,OP_1),g(:,OP_DP),f(:,OP_DP)) &
       + gam*     intx3(e(:,:,OP_1),g(:,OP_1),f(:,OP_DPP))
#endif

  v2vp = temp
end function v2vp


! V2chip
! ======
function v2chip(e,f,g)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chip
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp =     intx4(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR))    &
       + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ))    &
       + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP))    &
       + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP))    &
#ifdef USEST
       - gam*intx4(e(:,:,OP_DP),ri2_79,f(:,OP_GS),g(:,OP_1 ))    
#else
       + gam*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GSP),g(:,OP_1 ))    &
       + gam*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GS),g(:,OP_DP))
#endif
#endif

  v2chip = temp
end function v2chip


! V2p
! ===
function v2p(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2p
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  v2p = -intx2(e(:,:,OP_1),f(:,OP_DP))
#else
  v2p = 0.
#endif

end function v2p


! V2psipsi
! ========
function v2psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  v2psipsi = - &
       (intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ)) &
       +intx4(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR)))
#else
  v2psipsi = 0.
#endif
end function v2psipsi


! V2psib
! ======
function v2psib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  v2psib = intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
       -   intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR))
end function v2psib


! V2vpsipsi
! =========
function v2vpsipsi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vpsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
  integer :: j

  do j=1, dofs_per_element
     ! [nu,psi(2)]
     tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
  end do

  temp = intx3(tempa,f(:,OP_DR),g(:,OP_DZ)) &
       - intx3(tempa,f(:,OP_DZ),g(:,OP_DR))

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
  temp79b = &
         f(:,OP_DP)*(g(:,OP_DZ )*h(:,OP_DZ ) + g(:,OP_DR )*h(:,OP_DR )) &
       + f(:,OP_1 )*(g(:,OP_DZP)*h(:,OP_DZ ) + g(:,OP_DRP)*h(:,OP_DR )) 
  temp = temp - intx3(e(:,:,OP_DP),ri2_79,temp79b)
#else
  temp79b = &
       f(:,OP_DPP)*(g(:,OP_DZ )*h(:,OP_DZ ) + g(:,OP_DR )*h(:,OP_DR )) &
       + 2.*f(:,OP_DP)*(g(:,OP_DZP)*h(:,OP_DZ ) + g(:,OP_DRP)*h(:,OP_DR )) &
       +   f(:,OP_DP )*(g(:,OP_DZ )*h(:,OP_DZP) + g(:,OP_DR )*h(:,OP_DRP)) &
       + f(:,OP_1 )*(g(:,OP_DZPP)*h(:,OP_DZ ) + g(:,OP_DRPP)*h(:,OP_DR )) &
       + f(:,OP_1 )*(g(:,OP_DZP )*h(:,OP_DZP) + g(:,OP_DRP )*h(:,OP_DRP))
  temp = temp + intx3(e(:,:,OP_1),ri2_79,temp79b)
#endif
#endif

  v2vpsipsi = temp
end function v2vpsipsi


! V2vpsib
! =======
function v2vpsib(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vpsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp79a = f(:,OP_DP)*(g(:,OP_DZ )*h(:,OP_DR)-g(:,OP_DR )*h(:,OP_DZ)) &
       +    f(:,OP_1 )*(g(:,OP_DZP)*h(:,OP_DR)-g(:,OP_DRP)*h(:,OP_DZ))
  temp = intx3(e(:,:,OP_1),ri_79,temp79a)

  v2vpsib = temp
#else
  v2vpsib = 0.
#endif
end function v2vpsib

#ifdef USE3D
! V2vpsif
! =====
function v2vpsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vpsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j


  do j=1, dofs_per_element
     ! (nu, f') + nu*f'_LP
     tempa(j,:) = e(j,:,OP_1)*h(:,OP_LP) &
          + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
     ! [nu, psi]*R
     tempb(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
  end do

  ! [v, psi]*R
  temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)            
  ! (v, f') + v*f'_LP
  temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
       + f(:,OP_1)*h(:,OP_LP)
  ! [psi,f']*R
  temp79c = g(:,OP_DZ)*h(:,OP_DR) - g(:,OP_DR)*h(:,OP_DZ)            
  ! [psi,f']'*R
  temp79d = g(:,OP_DZP)*h(:,OP_DR) - g(:,OP_DRP)*h(:,OP_DZ) & 
       +    g(:,OP_DZ)*h(:,OP_DRP) - g(:,OP_DR)*h(:,OP_DZP)        

  temp = + intx3(tempa,temp79a,r_79) &
       + intx3(tempb,temp79b,r_79) &
       - 2*intx4(e(:,:,OP_DP),f(:,OP_DP),temp79c,ri_79) &
       - intx4(e(:,:,OP_DP),f(:,OP_1),temp79d,ri_79) 

  v2vpsif = temp
end function v2vpsif

! V2vbf
! =====
function v2vbf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vbf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp

  ! (F,f')
  temp79a = g(:,OP_DZ)*h(:,OP_DZ) + g(:,OP_DR)*h(:,OP_DR)
  ! (F,f'')
  temp79b = g(:,OP_DZ)*h(:,OP_DZP) + g(:,OP_DR)*h(:,OP_DRP)        

  temp = intx3(e(:,:,OP_1),f(:,OP_DP),temp79a) &
       + intx3(e(:,:,OP_1),f(:,OP_1),temp79b) 

  v2vbf = temp
end function v2vbf

! V2vff
! =====
function v2vff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
  integer :: j


  do j=1, dofs_per_element
     ! (nu, f') + nu*f'_LP
     tempa(j,:) = e(j,:,OP_1)*g(:,OP_LP) &
          + e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR) 
  end do
  ! (v, f') + v*f'_LP
  temp79a = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
       + f(:,OP_1)*h(:,OP_LP)
  ! (f',f')
  temp79b = g(:,OP_DZ)*h(:,OP_DZ) + g(:,OP_DR)*h(:,OP_DR)            
  ! (f',f'') 
  temp79c = g(:,OP_DZ)*h(:,OP_DZP) + g(:,OP_DR)*h(:,OP_DRP) 

  temp = - intx3(tempa,temp79a,r2_79) &
       - intx3(e(:,:,OP_DP),f(:,OP_DP),temp79b) &
       - intx3(e(:,:,OP_DP),f(:,OP_1),temp79c) 

  v2vff = temp
end function v2vff
#endif

! V2upsipsi
! =========
function v2upsipsi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2upsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)  
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempc, tempd
  integer :: j

  temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
  temp79b = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)

  do j=1, dofs_per_element
#ifdef USEST
     tempc(j,:) =  &
#else
     tempc(j,:) = e(j,:,OP_1)*h(:,OP_GS) &
#endif
             +    e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)

#ifdef USEST
     tempd(j,:) = - e(j,:,OP_DP)*h(:,OP_GS) &
#else
     tempd(j,:) = e(j,:,OP_1)*h(:,OP_GSP) &
#endif
             +    e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP)
  end do
  temp = intx3(tempd,ri_79,temp79a) &
       + intx3(tempc,ri_79,temp79b)

  v2upsipsi = temp
#else
  v2upsipsi = 0.
#endif
end function v2upsipsi



! V2upsib
! =======
function v2upsib(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2upsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa
  integer :: j

  do j=1, dofs_per_element
     tempa(j,:) = e(j,:,OP_DZ)*f(:,OP_DR) - e(j,:,OP_DR)*f(:,OP_DZ)
  end do

  temp = (intx3(tempa,g(:,OP_DR),h(:,OP_DZ)) &
       -  intx3(tempa,g(:,OP_DZ),h(:,OP_DR)))

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
  temp79b = &
       h(:,OP_DP)*(f(:,OP_DZ )*g(:,OP_DZ ) + f(:,OP_DR )*g(:,OP_DR )) &
     + h(:,OP_1 )*(f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR )) 
  temp = temp + intx3(e(:,:,OP_DP),ri2_79,temp79b)
#else
  temp79b = &
       2.*(f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ))*h(:,OP_DP ) &
       +  (f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP))*h(:,OP_DP ) &
       +  (f(:,OP_DZ )*g(:,OP_DZ ) + f(:,OP_DR )*g(:,OP_DR ))*h(:,OP_DPP) &
       +  (f(:,OP_DZPP)*g(:,OP_DZ ) + f(:,OP_DRPP)*g(:,OP_DR ))*h(:,OP_1) &
       +  (f(:,OP_DZP )*g(:,OP_DZP) + f(:,OP_DRP )*g(:,OP_DRP))*h(:,OP_1)
  temp = temp - intx3(e(:,:,OP_1),ri2_79,temp79b)
#endif
#endif

  v2upsib = temp
end function v2upsib


! V2ubb
! =====
function v2ubb(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2ubb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp79a = h(:,OP_DP)*(f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)) &
       +    h(:,OP_1)*(f(:,OP_DRP)*g(:,OP_DZ) - f(:,OP_DZP)*g(:,OP_DR))
  temp = intx3(e(:,:,OP_1),ri_79,temp79a)

  v2ubb = temp
#else
  v2ubb = 0.
#endif
end function v2ubb

#ifdef USE3D
! V2upsif
! =====
function v2upsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2upsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp


  ! ([u, psi]*R^2)_R/R
  temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
          + f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79a = temp79a + ri_79* &
          (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
  end if
  ! ([u, psi]*R^2)_Z/R
  temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
          + f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            
  ! ((u, f')*R^2)_R/R^2
  temp79c = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
          + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79c = temp79c + 2*ri_79* &
          (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ)) 
  end if
  ! ((u, f')*R^2)_Z/R^2
  temp79d = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
          + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            

  temp = - intx3(e(:,:,OP_DP),h(:,OP_DZ),temp79a) &
       + intx3(e(:,:,OP_DP),h(:,OP_DR),temp79b) &
       - intx3(e(:,:,OP_DP),g(:,OP_DR),temp79c) &
       - intx3(e(:,:,OP_DP),g(:,OP_DZ),temp79d) 

  v2upsif = temp
end function v2upsif

! V2ubf
! =====
function v2ubf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2ubf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
  integer :: j


  do j=1, dofs_per_element
     ! (nu, f') + nu*f'_LP
     tempa(j,:) = e(j,:,OP_1)*h(:,OP_LP) &
          + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
  end do
  ! (F'[u,f'] + F[u',f'])*R 
  temp79a = g(:,OP_1)* &
       (h(:,OP_DR)*f(:,OP_DZP) - h(:,OP_DZ)*f(:,OP_DRP)) &
       + g(:,OP_DP)* &
       (h(:,OP_DR)*f(:,OP_DZ) - h(:,OP_DZ)*f(:,OP_DR))
  ! (R^2(u, f'))_R/R^2
  temp79b = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
       + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79b = temp79b + 2*ri_79* &
          (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
  end if
  ! (R^2(u, f'))_Z/R^2
  temp79c = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
       + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
  ![u,F]*R
  temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ) 

  temp = intx3(e(:,:,OP_DP),temp79a,ri_79) &
       + intx4(e(:,:,OP_1),temp79c,g(:,OP_DR),r_79) &
       - intx4(e(:,:,OP_1),temp79b,g(:,OP_DZ),r_79) &
       + intx3(tempa,temp79d,r_79)

  v2ubf = temp
end function v2ubf

! V2uff
! =====
function v2uff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2uff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp

  !((u, f')*R^2)_R/R^2
  temp79a = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &           
          + f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79a = temp79a + 2*ri_79* &
          (f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ))
  end if
  !((u, f')*R^2)_Z/R^2
  temp79b = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &           
          + f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ)            

  temp = intx4(e(:,:,OP_DP),h(:,OP_DZ),temp79a,r_79) &
       - intx4(e(:,:,OP_DP),h(:,OP_DR),temp79b,r_79) 

  v2uff = temp
end function v2uff
#endif

! v2upsisb2
! ========
function v2upsisb2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2upsisb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  v2upsisb2 = 0.
end function v2upsisb2


! v2ubsb1
! =======
function v2ubsb1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2ubsb1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  v2ubsb1 = 0.
end function v2ubsb1


! v2chipsipsi
! ===========
function v2chipsipsi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempc, tempd
  integer :: j

  temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
  temp79b = f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)

  do j=1, dofs_per_element
#ifdef USEST
     tempc(j,:) =  &
#else
     tempc(j,:) = e(j,:,OP_1)*h(:,OP_GS) &
#endif
          + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
#ifdef USEST
     tempd(j,:) = -e(j,:,OP_DP)*h(:,OP_GS) &
#else
     tempd(j,:) = e(j,:,OP_1)*h(:,OP_GSP) &
#endif
          + e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP)
  end do
        
  temp = -intx3(tempd,ri4_79,temp79a) &
       -  intx3(tempc,ri4_79,temp79b)

  v2chipsipsi = temp
#else
  v2chipsipsi = 0.
#endif
end function v2chipsipsi


! v2chipsib
! =========
function v2chipsib(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempb, tempc
  integer :: j


  temp79a = h(:,OP_1 )*f(:,OP_GS) &
       +    h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)

  do j=1, dofs_per_element
     tempb(j,:) = e(j,:,OP_DR)*h(:,OP_DZ) - e(j,:,OP_DZ)*h(:,OP_DR)

     tempc(j,:) = e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)
  end do
  

  temp = intx3(tempc,ri3_79,temp79a) &
       + intx4(tempb,ri3_79,f(:,OP_DZ),g(:,OP_DZ)) &
       + intx4(tempb,ri3_79,f(:,OP_DR),g(:,OP_DR))

  if(itor.eq.1) then
     temp = temp - &
          2.*intx4(tempc,ri4_79,f(:,OP_DR),h(:,OP_1))
  endif
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
  temp79d = &
         (f(:,OP_DZ)*g(:,OP_DR ) - f(:,OP_DR)*g(:,OP_DZ ))*h(:,OP_DP ) &
       + (f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ))*h(:,OP_1 ) 
  temp = temp + intx3(e(:,:,OP_DP),ri5_79,temp79d)
#else
  temp79d = &
       2.*(f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ))*h(:,OP_DP ) &
       +  (f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP))*h(:,OP_DP ) &
       +  (f(:,OP_DZ )*g(:,OP_DR ) - f(:,OP_DR )*g(:,OP_DZ ))*h(:,OP_DPP) &
       +  (f(:,OP_DZPP)*g(:,OP_DR ) - f(:,OP_DRPP)*g(:,OP_DZ ))*h(:,OP_1) &
       +  (f(:,OP_DZP )*g(:,OP_DRP) - f(:,OP_DRP )*g(:,OP_DZP))*h(:,OP_1)
  temp = temp - intx3(e(:,:,OP_1),ri5_79,temp79d)
#endif
#endif

  v2chipsib = temp
end function v2chipsib


! v2chibb
! =======
function v2chibb(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chibb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx5(e(:,:,OP_1),ri4_79,g(:,OP_1),h(:,OP_DR),f(:,OP_DRP)) &
       + intx5(e(:,:,OP_1),ri4_79,g(:,OP_1),h(:,OP_DZ),f(:,OP_DZP)) &
       + intx5(e(:,:,OP_1),ri4_79,g(:,OP_DP),h(:,OP_DR),f(:,OP_DR)) &
       + intx5(e(:,:,OP_1),ri4_79,g(:,OP_DP),h(:,OP_DZ),f(:,OP_DZ))

  v2chibb = temp
#else
  v2chibb = 0.
#endif
end function v2chibb

#ifdef USE3D
! V2chipsif
! =====
function v2chipsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chipsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp

  ! [chi, f']_R*R 
  temp79a = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
          + f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
  if(itor.eq.1) then 
     temp79a = temp79a - ri_79* &        
          (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
  end if
  ! [chi, f']_Z*R 
  temp79b = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
  ! ((chi, psi)/R^2)_R * R^2
  temp79c = f(:,OP_DRZ)*g(:,OP_DZ) + f(:,OP_DRR)*g(:,OP_DR) &
       +    f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DR)*g(:,OP_DRR) 
  if(itor.eq.1) then 
     temp79c = temp79c - 2*ri_79*  & 
          (f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)) 
  end if
  ! ((chi, psi)/R^2)_Z * R^2
  temp79d = f(:,OP_DZZ)*g(:,OP_DZ) + f(:,OP_DRZ)*g(:,OP_DR) &
       +    f(:,OP_DZ)*g(:,OP_DZZ) + f(:,OP_DR)*g(:,OP_DRZ) 

  temp = &
       - intx4(e(:,:,OP_DP),temp79a,ri3_79,g(:,OP_DR)) &
       - intx4(e(:,:,OP_DP),temp79b,ri3_79,g(:,OP_DZ)) &
       - intx4(e(:,:,OP_DP),temp79d,ri3_79,h(:,OP_DR)) &
       + intx4(e(:,:,OP_DP),temp79c,ri3_79,h(:,OP_DZ)) 

  v2chipsif = temp
end function v2chipsif

! V2chibf
! =====
function v2chibf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chibf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j


  do j=1, dofs_per_element
     ! (nu, f') + nu*f'_LP
     tempa(j,:) = e(j,:,OP_1)*h(:,OP_LP) &
          + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
     ! [nu, F]*R
     tempb(j,:) = &
          + e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
  end do

  ! [chi, f']*R
  temp79a = f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)  
  ! (chi, f')
  temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) 
  ! (chi', f')
  temp79c = f(:,OP_DRP)*h(:,OP_DR) + f(:,OP_DZP)*h(:,OP_DZ) 
  ! (chi, F/R^4)*R^4 + F*chi_LP
  temp79d = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
       + f(:,OP_LP)*g(:,OP_1) 
  if(itor.eq.1) then 
     temp79d = temp79d - 4*ri_79*f(:,OP_DR)*g(:,OP_1) 
  end if

  temp = - intx3(tempb,temp79a,ri2_79) &
       - intx3(tempa,temp79d,ri2_79) &
       - intx4(e(:,:,OP_DP),temp79b,ri4_79,g(:,OP_DP)) &
       - intx4(e(:,:,OP_DP),temp79c,ri4_79,g(:,OP_1)) 

  v2chibf = temp
end function v2chibf

! V2chiff
! =====
function v2chiff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chiff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     ! [nu, f']*R
     tempa(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)  
     ! [nu, f'']*R
     tempb(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP)  
  end do

  ! [chi, f']*R
  temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)        
  ! [chi, f']'*R
  temp79b = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP)            

  temp = - intx3(tempa,temp79b,ri2_79) &
       - intx3(tempb,temp79a,ri2_79) 

  v2chiff = temp
end function v2chiff
#endif

! v2vchin
! =======
function v2vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v2vchin = 0.
     return
  end if

  temp =-intx4(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
       - intx4(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  if(itor.eq.1) then
     temp = temp &
          - 2.*intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DR),h(:,OP_1))
  endif

  v2vchin = temp
end function v2vchin


! v2chibsb1
! =========
function v2chibsb1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chibsb1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  v2chibsb1 = 0.
end function v2chibsb1

! v2psisb2
! ========
function v2psisb2(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psisb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
       - intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR))

  v2psisb2 = temp
end function v2psisb2

! v2bsb1
! ======
function v2bsb1(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2bsb1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR)) &
       - intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ))

  v2bsb1 = temp
end function v2bsb1


! V2vs
! ====
function v2vs(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vs
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(idens.eq.0 .or. nosig.eq.1) then
     v2vs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  v2vs = -intx4(e(:,:,OP_1),r2_79,f(:,OP_1),temp79a)

end function v2vs


! V2psif1
! =======
function v2psif1(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psif1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif1 = intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZP)) &
       -    intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DRP))
#else
  v2psif1 = 0.
#endif
end function v2psif1

! V2psif2
! =======
function v2psif2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psif2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif2 = intx4(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ)) &
       -    intx4(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR))
#else
  v2psif2 = 0.
#endif
end function v2psif2



! V2bf
! ====
function v2bf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2bf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  
#if defined(USE3D) || defined(USECOMPLEX)
  v2bf = &
       - intx3(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
       - intx3(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR))
#else
  v2bf = 0.
#endif
end function v2bf

! V2ff
! ====
function v2ff(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2ff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  
#if defined(USE3D) || defined(USECOMPLEX)
  v2ff = &
       - intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_DZ)) &
       - intx3(e(:,:,OP_1),f(:,OP_DRP),g(:,OP_DR))
#else
  v2ff = 0.
#endif
end function v2ff


! V2be
! ===

function v2parpb2ipsipsi(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2parpb2ipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = -f(:,OP_DP)*g(:,OP_1)*ri2_79
  temp = intx4(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR))    &
       + intx4(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ)) 
#endif

  v2parpb2ipsipsi = temp
end function v2parpb2ipsipsi

function v2parpb2ibb(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2parpb2ibb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = -f(:,OP_DP)*g(:,OP_1)*h(:,OP_1)*i(:,OP_1)*ri2_79
  temp = temp + intx2(e(:,:,OP_1),temp79a)  
#endif

  v2parpb2ibb = temp
end function v2parpb2ibb


function v2parpb2ipsib(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2parpb2ipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i

  temp79a =  f(:,OP_1)*g(:,OP_1)*i(:,OP_1)*ri_79

  v2parpb2ipsib = intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DR))    &
       -          intx3(e(:,:,OP_DR),temp79a,h(:,OP_DZ)) 
end function v2parpb2ipsib

#ifdef USEPARTICLES
function v2p_2(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2p_2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  if(surface_int) then
     v2p_2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  v2p_2 = -intx3(e(:,:,OP_1),f(:,OP_DP),g)
#else
  v2p_2 = 0.
#endif

end function v2p_2

function v2pbb(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2pbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = -f(:,OP_DP)*g*bzt79(:,OP_1)*bzt79(:,OP_1)*ri2_79
  temp = temp + intx2(e(:,:,OP_1),temp79a)  
  temp79a = g*bzt79(:,OP_1)
  temp = temp + intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),bfpt79(:,OP_DR))    &
       +  intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),bfpt79(:,OP_DZ)) 
#endif
  temp79a = -g*bzt79(:,OP_1)*ri_79
  temp = temp + intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),pst79(:,OP_DR))    &
       -          intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),pst79(:,OP_DZ)) 
    end if

  v2pbb = temp
end function v2pbb

function v2jxb(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2jxb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = f*ri_79
        temp = temp                                   &
             +  intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DZ),bzt79(:,OP_DR)+bfpt79(:,OP_DRP))   &
             -  intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DR),bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))   
  temp79a = f*ri2_79
  temp = temp + intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DR),pst79(:,OP_DRP))    &
       + intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DZ),pst79(:,OP_DZP)) 
  temp79a = f
  temp = temp + intx4(e(:,:,OP_1),temp79a,bfpt79(:,OP_DR),(bzt79(:,OP_DR)+bfpt79(:,OP_DRP)))    &
       + intx4(e(:,:,OP_1),temp79a,bfpt79(:,OP_DZ),(bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))) 
  temp79a = f*ri_79
  temp = temp + intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DZP),bfpt79(:,OP_DR))    &
       - intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DRP),bfpt79(:,OP_DZ)) 
#endif

    end if

  v2jxb = temp
end function v2jxb

function v2pbb1psi(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2pbb1psi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
    if(surface_int) then
        temp = 0.
     else
  temp79a = -g*bzt79(:,OP_1)*ri_79
  temp = temp + intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),h(:,OP_DR))    &
       -          intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),h(:,OP_DZ)) 
    end if

  v2pbb1psi = temp
end function v2pbb1psi


function v2pbb1f(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2pbb1f
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = g*bzt79(:,OP_1)
  temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),h(:,OP_DR))    &
       +  intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),h(:,OP_DZ)) 
#endif
     end if

  v2pbb1f = temp

end function v2pbb1f

#endif
!==============================================================================
! V3 TERMS
!==============================================================================

! V3chin
! ======
function v3chin(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = - intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_1)) &
       -   intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_1))

  v3chin = temp
end function v3chin



! V3chimu
! =======
function v3chimu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chimu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp79b = f(:,OP_DRR)
  temp79d = f(:,OP_DRZ)
  if(itor.eq.1) then
     temp79b = temp79b - 2.*ri_79*f(:,OP_DR)
     temp79d = temp79d -    ri_79*f(:,OP_DZ)
  endif
  temp = 2.* &
       (intx4(e(:,:,OP_DZZ),ri4_79,f(:,OP_DZZ),g(:,OP_1)) &
       +intx4(e(:,:,OP_DRR),ri4_79,temp79b,g(:,OP_1)) &
       +2.*intx4(e(:,:,OP_DRZ),ri4_79,temp79d,g(:,OP_1)) & 
       +intx4(e(:,:,OP_GS),ri4_79,f(:,OP_GS),h(:,OP_1)) &
       -intx4(e(:,:,OP_GS),ri4_79,f(:,OP_GS),g(:,OP_1)))
  if(itor.eq.1) then
     temp = temp &
          + 2.*intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DR),g(:,OP_1)) &
          - 4.*intx4(e(:,:,OP_DR),ri5_79,temp79b,g(:,OP_1)) &
          - 4.*intx4(e(:,:,OP_DZ),ri5_79,temp79d,g(:,OP_1))
  endif
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
  temp = temp + &
       (intx4(e(:,:,OP_DZP),ri6_79,f(:,OP_DZP),g(:,OP_1)) &
       +intx4(e(:,:,OP_DRP),ri6_79,f(:,OP_DRP),g(:,OP_1))) &
       +(intx4(e(:,:,OP_DZ),ri6_79,f(:,OP_DZP),g(:,OP_DP)) &
       +intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DRP),g(:,OP_DP)))
#else
  temp = temp - &
       (intx4(e(:,:,OP_DZ),ri6_79,f(:,OP_DZPP),g(:,OP_1)) &
       +intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DRPP),g(:,OP_1)))
#endif
#endif

  v3chimu = temp
end function v3chimu


! V3umu
! =====
function v3umu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3umu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

     temp79c = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79c = temp79c - ri_79*f(:,OP_DR)

     temp = 2.* &
          (intx4(e(:,:,OP_DZZ),ri_79,f(:,OP_DRZ),g(:,OP_1)) &
          -intx4(e(:,:,OP_DRR),ri_79,f(:,OP_DRZ),g(:,OP_1)) &
          -intx4(e(:,:,OP_DRZ),ri_79,temp79c,g(:,OP_1)))
     if(itor.eq.1) then
        temp = temp - 2.* &
             (intx4(e(:,:,OP_DRR),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
             -intx4(e(:,:,OP_DR ),ri3_79,f(:,OP_DZ),g(:,OP_1))) &
             + 4.*intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRZ),g(:,OP_1)) &
             + 2.*intx4(e(:,:,OP_DZ),ri2_79,temp79c,g(:,OP_1))
        
        temp79d = g(:,OP_1)-h(:,OP_1)
        temp = temp + 4.*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DZ),temp79d)
     endif
        
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
     temp = temp &
          - intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZP),g(:,OP_1)) &
          + intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DRP),g(:,OP_1)) &
          - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_DP)) &
          + intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_DP))
#else
     temp = temp &
          + intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZPP),g(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRPP),g(:,OP_1))
#endif
#endif

  v3umu = temp
end function v3umu


! V3vmu
! =====
function v3vmu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vmu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp = &
       - intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1)) &
       + 2.*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),h(:,OP_1)) &
       - 2.*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),g(:,OP_1))
  if(itor.eq.1) then
     temp = temp &
          + 2.*intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DP),g(:,OP_1))
  endif

#endif
  v3vmu = temp
end function v3vmu

! V3un
! ====
function v3un(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3un
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1)) 

  v3un = temp
end function v3un


! V3p
! ===
function v3p(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3p
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = intx3(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ)) &
       + intx3(e(:,:,OP_DR),ri2_79,f(:,OP_DR))

  v3p = temp
end function v3p



! V3up
! ====
function v3up(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3up
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_GS),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
       - intx4(e(:,:,OP_GS),ri_79,f(:,OP_DZ),g(:,OP_DR))
  if(itor.eq.1) then
     temp = temp - 2.*gam* &
          intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DZ),g(:,OP_1))
  endif

  v3up = temp
end function v3up


! V3vp
! ====
function v3vp(e,f,g)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vp
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g

  vectype, dimension(dofs_per_element) :: temp
  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp = gam*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),g(:,OP_1)) &
       + intx4(e(:,:,OP_GS),ri2_79,f(:,OP_1),g(:,OP_DP))
#endif

  v3vp = temp
end function v3vp



! V3chip
! ======
function v3chip(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chip
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DZ),g(:,OP_DZ)) &
       + intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DR),g(:,OP_DR)) &
       + gam*intx4(e(:,:,OP_GS),ri4_79,f(:,OP_GS),g(:,OP_1))

  v3chip = temp
end function v3chip



! V3psipsi
! ========
function v3psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3psipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_GS)) &
       + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_GS))

  v3psipsi = temp
end function v3psipsi


! V3psib
! ======
function v3psib(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3psib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = - &
       (intx4(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP),g(:,OP_1)) &
       -intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DZP),g(:,OP_1)))
#else
  temp = 0.
#endif
  v3psib = temp
end function v3psib


! V3psif
! ======
function v3psif(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3psif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DR),f(:,OP_GS)) &
       - intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DZ),f(:,OP_GS))
#else
  temp = 0.
#endif

  v3psif = temp
end function v3psif


! V3bb
! ====
function v3bb(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3bb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_1)) &
       + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_1))

  v3bb = temp
end function v3bb


! V3bf
! ====
function v3bf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3bf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx4(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_1)) &
       + intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_1))
#else
  temp = 0.
#endif

  v3bf = temp
end function v3bf


! V3psisb1
! ========
vectype function v3psisb1(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ)) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR)) &
       + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS)) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS))

  v3psisb1 = temp

end function v3psisb1


! V3bsb2
! ======
vectype function v3bsb2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 )) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1 )) &
       + int4(ri2_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ)) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_1 ),g(:,OP_DR))

  v3bsb2 = temp
  return
end function v3bsb2


! V3upsipsi
! =========
function v3upsipsi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3upsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempd
  integer :: j
 
  ! [f,g],r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
  
  ! [f,g],z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)
  
  temp = intx4(e(:,:,OP_DZ ),ri3_79,temp79c,h(:,OP_GS )) &
       + intx4(e(:,:,OP_DR ),ri3_79,temp79b,h(:,OP_GS )) &
       - intx4(e(:,:,OP_DZZ),ri3_79,temp79c,h(:,OP_DZ )) &
       - intx4(e(:,:,OP_DR ),ri3_79,temp79c,h(:,OP_DRZ)) &
       - intx4(e(:,:,OP_DZ ),ri3_79,temp79c,h(:,OP_DZZ)) &
       - intx4(e(:,:,OP_DRZ),ri3_79,temp79c,h(:,OP_DR )) &
       - intx4(e(:,:,OP_DRZ),ri3_79,temp79b,h(:,OP_DZ )) &
       - intx4(e(:,:,OP_DR ),ri3_79,temp79b,h(:,OP_DRR)) &
       - intx4(e(:,:,OP_DZ ),ri3_79,temp79b,h(:,OP_DRZ)) &
       - intx4(e(:,:,OP_DRR),ri3_79,temp79b,h(:,OP_DR ))
    
  if(itor.eq.1) then
     ! [f,g]
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
  

        do j=1, dofs_per_element
           tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
        end do
        
        temp = temp &
             +2.*intx3(tempd,ri4_79,temp79b) &
             -   intx4(e(:,:,OP_DRZ),ri4_79,temp79a,h(:,OP_DZ )) &
             -   intx4(e(:,:,OP_DZ ),ri4_79,temp79a,h(:,OP_DRZ)) &
             -   intx4(e(:,:,OP_DRR),ri4_79,temp79a,h(:,OP_DR )) &
             -   intx4(e(:,:,OP_DR ),ri4_79,temp79a,h(:,OP_DRR)) &
             +   intx4(e(:,:,OP_DR),ri4_79,temp79a,h(:,OP_GS)) &
             +2.*intx3(tempd,ri5_79,temp79a)
  endif

  v3upsipsi = temp
end function v3upsipsi


! V3upsib
! =======
function v3upsib(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3upsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX) 
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempd, tempf
  integer :: j

  do j=1, dofs_per_element
     tempa(j,:) = e(j,:,OP_DZ)*g(:,OP_DRP) - e(j,:,OP_DR)*g(:,OP_DZP)
     tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DR)  - e(j,:,OP_DR)*h(:,OP_DZ)
     tempf(j,:) = h(:,OP_DP)* &
          (e(j,:,OP_DZ)*f(:,OP_DZ )+e(j,:,OP_DR)*f(:,OP_DR )) &
          +    h(:,OP_1 )* &
          (e(j,:,OP_DZ)*f(:,OP_DZP)+e(j,:,OP_DR)*f(:,OP_DRP))
     if(itor.eq.1) then 
        tempd(j,:) = tempd(j,:) - 4.*ri_79*e(j,:,OP_DZ)*h(:,OP_1)
     endif
  end do
  
  temp79b = h(:,OP_DZ)*f(:,OP_DR) - h(:,OP_DR)*f(:,OP_DZ)
  temp79c = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &
       +    f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP)
#ifdef USEST
  temp79e =  h(:,OP_1)*f(:,OP_GS) &
       +     f(:,OP_DZ )*h(:,OP_DZ) + f(:,OP_DR )*h(:,OP_DR)
#else
  temp79e = h(:,OP_DP)*f(:,OP_GS) + h(:,OP_1)*f(:,OP_GSP) &
       +    f(:,OP_DZP)*h(:,OP_DZ ) + f(:,OP_DRP)*h(:,OP_DR ) &
       +    f(:,OP_DZ )*h(:,OP_DZP) + f(:,OP_DR )*h(:,OP_DRP)
#endif
  temp = intx3(tempa,ri4_79,temp79b) &
       + intx3(tempd,ri4_79,temp79c) &
#ifdef USEST
       - intx4(e(:,:,OP_DZP),ri4_79,temp79e,g(:,OP_DZ)) &
       - intx4(e(:,:,OP_DRP),ri4_79,temp79e,g(:,OP_DR)) &
       - intx4(e(:,:,OP_DZ),ri4_79,temp79e,g(:,OP_DZP)) &
       - intx4(e(:,:,OP_DR),ri4_79,temp79e,g(:,OP_DRP)) &
#else
       + intx4(e(:,:,OP_DZ),ri4_79,temp79e,g(:,OP_DZ)) &
       + intx4(e(:,:,OP_DR),ri4_79,temp79e,g(:,OP_DR)) &
#endif
       + intx3(tempf,ri4_79,g(:,OP_GS))

  v3upsib = temp
#else
  v3upsib = 0.
#endif
end function v3upsib



! V3ubb
! =====
function v3ubb(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ubb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb
  integer :: j

  temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

  temp = intx3(e(:,:,OP_GS),ri3_79,temp79a)

  if(itor.eq.1) then
     temp = temp - &
          2.*intx3(e(:,:,OP_DR),ri4_79,temp79a)
  endif

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
  do j=1, dofs_per_element 
     tempb(j,:) = &
          (e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
          *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP)) &
          -(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
          *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP))
  end do
  temp = temp + intx2(tempb,ri5_79)
#else
  do j=1, dofs_per_element
     tempb(j,:) = &
          (e(j,:,OP_DZ)*f(:,OP_DR)-e(j,:,OP_DR)*f(:,OP_DZ))*g(:,OP_DPP) &
          + 2.*(e(j,:,OP_DZ)*f(:,OP_DRP)-e(j,:,OP_DR)*f(:,OP_DZP))*g(:,OP_DP) &
          +    (e(j,:,OP_DZ)*f(:,OP_DRPP)-e(j,:,OP_DR)*f(:,OP_DZPP))*g(:,OP_1)
  end do
  temp = temp - intx3(tempb,ri5_79,h(:,OP_1))
#endif
#endif

  v3ubb = temp
end function v3ubb

#ifdef USE3D
! V3upsif
! =====
function v3upsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3upsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempd
  integer :: j


  do j=1, dofs_per_element
     ! ((nu, psi)/R^2)_R*R^2
     tempa(j,:) = &
            e(j,:,OP_DZ)*g(:,OP_DRZ) + e(j,:,OP_DR)*g(:,OP_DRR) &
          + e(j,:,OP_DRZ)*g(:,OP_DZ) + e(j,:,OP_DRR)*g(:,OP_DR) 
     if(itor.eq.1) then
        tempa(j,:) = tempa(j,:) - 2*ri_79*&
             (e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DZ)) 
     end if
     ! ((nu, psi)/R^2)_Z*R^2
     tempb(j,:) = &
            e(j,:,OP_DZ)*g(:,OP_DZZ) + e(j,:,OP_DR)*g(:,OP_DRZ) &
          + e(j,:,OP_DZZ)*g(:,OP_DZ) + e(j,:,OP_DRZ)*g(:,OP_DR) 
     ! [nu, f']_R*R
     tempc(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DRR) - e(j,:,OP_DR)*h(:,OP_DRZ) &
          + e(j,:,OP_DRZ)*h(:,OP_DR) - e(j,:,OP_DRR)*h(:,OP_DZ)
     if(itor.eq.1) then
        tempc(j,:) = tempc(j,:) - ri_79*&
             (e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)) 
     end if
     ! [nu, f']_Z*R
     tempd(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DRZ) - e(j,:,OP_DR)*h(:,OP_DZZ) &
          + e(j,:,OP_DZZ)*h(:,OP_DR) - e(j,:,OP_DRZ)*h(:,OP_DZ)
  end do

  ! ([u, psi]*R^2)_R/R
  temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79a = temp79a + ri_79* &
          (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
  end if
  ! ([u, psi]*R^2)_Z/R
  temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            
  ! (R^2*(u, f'))_R/R^2
  temp79c = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79c = temp79c + 2*ri_79* &           
          (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
  endif
  ! (R^2*(u, f'))_Z/R^2
  temp79d = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
       +    f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            

  temp = intx3(tempa,temp79c,ri2_79) &
       + intx3(tempb,temp79d,ri2_79) &
       - intx3(tempc,temp79a,ri2_79) &
       - intx3(tempd,temp79b,ri2_79) &
       - intx4(e(:,:,OP_DZ),temp79d,g(:,OP_GS),ri2_79) & 
       - intx4(e(:,:,OP_DR),temp79c,g(:,OP_GS),ri2_79) 

  v3upsif = temp
end function v3upsif

! V3ubf
! =====
function v3ubf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ubf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j


  do j=1, dofs_per_element
     ! [nu, f']'*R
     tempa(j,:) = &
            e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) &
          + e(j,:,OP_DZP)*h(:,OP_DR) - e(j,:,OP_DRP)*h(:,OP_DZ) 
     ! (nu,f'')
     tempb(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
  end do
  ! F*u_GS + (u, F)
  temp79a = g(:,OP_1)*f(:,OP_GS) & 
       +    g(:,OP_DR)*f(:,OP_DR) + g(:,OP_DZ)*f(:,OP_DZ)
  ! (R^2(u, f'))_R/R^2
  temp79b = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79b = temp79b + 2*ri_79* &
          (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
  end if
  ! (R^2(u, f'))_Z/R^2
  temp79c = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
       +    f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
  ![u,F]*R
  temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ) 

  temp = - intx3(tempa,temp79a,ri3_79) &
       + intx3(tempb,temp79d,ri3_79) &
       - intx4(e(:,:,OP_DZP),temp79b,g(:,OP_1),ri3_79) &
       + intx4(e(:,:,OP_DRP),temp79c,g(:,OP_1),ri3_79) &
       - intx4(e(:,:,OP_DZ),temp79b,g(:,OP_DP),ri3_79) &
       + intx4(e(:,:,OP_DR),temp79c,g(:,OP_DP),ri3_79) 

  v3ubf = temp
end function v3ubf

! V3uff
! =====
function v3uff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3uff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     ![nu, f']_R*R
     tempa(j,:) = &
            e(j,:,OP_DZ)*g(:,OP_DRR) - e(j,:,OP_DR)*g(:,OP_DRZ) &
          + e(j,:,OP_DRZ)*g(:,OP_DR) - e(j,:,OP_DRR)*g(:,OP_DZ)
     if(itor.eq.1) then
        tempa(j,:) = tempa(j,:) - ri_79*&
             (e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)) 
     end if
     ![nu, f']_Z*R
     tempb(j,:) = &
            e(j,:,OP_DZ)*g(:,OP_DRZ) - e(j,:,OP_DR)*g(:,OP_DZZ) &
          + e(j,:,OP_DZZ)*g(:,OP_DR) - e(j,:,OP_DRZ)*g(:,OP_DZ)
  end do
  ! (R^2*(u, f'))_R/R^2
  temp79a = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
  if(itor.eq.1) then
     temp79a = temp79a + 2*ri_79* &
          (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
  endif
  ! (R^2*(u, f'))_Z/R^2
  temp79b = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
       +    f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)     

  temp = intx3(tempa,temp79a,ri_79) &
       + intx3(tempb,temp79b,ri_79) 

  v3uff = temp
end function v3uff
#endif

! v3vpsipsi
! =========
function v3vpsipsi(e,f,g,h)
!
!  e trial
!  f lin
!  g psi
!  h psi
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vpsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempf
  integer :: j

  do j=1, dofs_per_element
     tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP)
     tempf(j,:) = f(:,OP_DP)* &
          (e(j,:,OP_DZ)*g(:,OP_DZ )+e(j,:,OP_DR)*g(:,OP_DR )) &
          +       f(:,OP_1 )* &
          (e(j,:,OP_DZ)*g(:,OP_DZP)+e(j,:,OP_DR)*g(:,OP_DRP))
  end do

  temp79b = g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ)
#ifdef USEST
  temp79e =  f(:,OP_1)*g(:,OP_GS) &
       + g(:,OP_DZ )*f(:,OP_DZ) + g(:,OP_DR )*f(:,OP_DR)
#else
  temp79e = f(:,OP_DP)*g(:,OP_GS) + f(:,OP_1)*g(:,OP_GSP) &
       + g(:,OP_DZP)*f(:,OP_DZ ) + g(:,OP_DRP)*f(:,OP_DR ) &
       + g(:,OP_DZ )*f(:,OP_DZP) + g(:,OP_DR )*f(:,OP_DRP)
#endif
  temp = intx3(tempa,ri4_79,temp79b) &
#ifdef USEST
       + intx4(e(:,:,OP_DZ),ri4_79,temp79e,h(:,OP_DZP)) &
       + intx4(e(:,:,OP_DR),ri4_79,temp79e,h(:,OP_DRP)) &
       + intx4(e(:,:,OP_DZP),ri4_79,temp79e,h(:,OP_DZ)) &
       + intx4(e(:,:,OP_DRP),ri4_79,temp79e,h(:,OP_DR)) &
#else
       - intx4(e(:,:,OP_DZ),ri4_79,temp79e,h(:,OP_DZ)) &
       - intx4(e(:,:,OP_DR),ri4_79,temp79e,h(:,OP_DR)) &
#endif
       - intx3(tempf,ri4_79,h(:,OP_GS))

  v3vpsipsi = temp
#else
  v3vpsipsi = 0.
#endif
end function v3vpsipsi



! v3vpsib
! =======
function v3vpsib(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vpsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb
  integer :: j

  temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

  temp = intx3(e(:,:,OP_GS),ri3_79,temp79a)

  if(itor.eq.1) then
     temp = temp - &
          2.*intx3(e(:,:,OP_DR),ri4_79,temp79a)
  endif

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
  do j=1, dofs_per_element 
     tempb(j,:) = &
          -(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
          *(g(:,OP_DRP)*f(:,OP_1) + g(:,OP_DR)*f(:,OP_DP)) &
          +(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
          *(g(:,OP_DZP)*f(:,OP_1) + g(:,OP_DZ)*f(:,OP_DP))
  end do
  temp = temp + intx2(tempb,ri5_79)
#else
  do j=1, dofs_per_element
     tempb(j,:) = f(:,OP_DPP)* &
          (e(j,:,OP_DZ)*g(:,OP_DR)-e(j,:,OP_DR)*g(:,OP_DZ)) &
          + 2.*f(:,OP_DP)* &
          (e(j,:,OP_DZ)*g(:,OP_DRP)-e(j,:,OP_DR)*g(:,OP_DZP)) &
          +    f(:,OP_1)* &
          (e(j,:,OP_DZ)*g(:,OP_DRPP)-e(j,:,OP_DR)*g(:,OP_DZPP))
  end do
  temp = temp + intx3(tempb,ri5_79,h(:,OP_1))
#endif
#endif

  v3vpsib = temp
end function v3vpsib


! V3vbb
! =====
function v3vbb(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  v3vbb= 0.
end function v3vbb

#ifdef USE3D
! V3vpsif
! =====
function v3vpsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vpsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe, tempf
  integer :: j

  do j=1, dofs_per_element
     ! [nu, psi']*R
     tempa(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DRP) - e(j,:,OP_DR)*g(:,OP_DZP) 
     ! (nu, f'')
     tempb(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
     ! (nu, psi)
     tempc(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR) 
     ! [nu, f']*R
     tempd(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ) 
     ! [nu, f'']*R
     tempe(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
     ! [nu', f']*R
     tempf(j,:) = &
          e(j,:,OP_DZP)*h(:,OP_DR) - e(j,:,OP_DRP)*h(:,OP_DZ) 
  end do

  ! [v, f']'*R 
  temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
  ! [v, psi]*R
  temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)            
  ! (v, f') + v*f'_LP
  temp79c = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
       +    f(:,OP_1)*h(:,OP_LP)
  ! (v, psi)
  temp79d = f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ) 

  temp = - intx3(tempc,temp79a,ri3_79) &
       + intx3(tempb,temp79b,ri3_79) &
       + intx3(tempa,temp79c,ri3_79) &
       + intx3(tempe,temp79d,ri3_79) &
       + intx3(tempf,temp79d,ri3_79) &
       - intx4(tempd,g(:,OP_GS),f(:,OP_DP),ri3_79) &
       + intx4(tempf,g(:,OP_GS),f(:,OP_1),ri3_79) 

  v3vpsif = temp
end function v3vpsif

! V3vbf
! =====
function v3vbf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vbf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
  integer :: j


  do j=1, dofs_per_element
     ! (nu, f')
     tempa(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
     ! (nu, f'')
     tempb(j,:) = &
          e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
     ! (nu', f'')
     tempc(j,:) = &
          e(j,:,OP_DZP)*h(:,OP_DZP) + e(j,:,OP_DRP)*h(:,OP_DRP) 
  end do
  ! (v, f') + v*f'_LP
  temp79a = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
       +    f(:,OP_1)*h(:,OP_LP)

  temp = - intx4(tempa,g(:,OP_1),f(:,OP_DPP),ri4_79) &
       - intx4(tempb,g(:,OP_1),f(:,OP_DP),ri4_79) &
       + intx4(tempc,g(:,OP_1),f(:,OP_1),ri4_79) &
       + intx4(e(:,:,OP_GS),g(:,OP_1),temp79a,ri2_79)
  if(itor.eq.1) then
     temp = temp &
          - 2*intx4(e(:,:,OP_DR),g(:,OP_1),temp79a,ri3_79)
  end if

  v3vbf = temp
end function v3vbf

! V3vff
! =====
function v3vff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     ! (nu, f'')
     tempa(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP) 
     ! [nu, f']*R
     tempb(j,:) = &
          e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
  end do
  ! [v, f']'*R
  temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
  ! (v, f') 
  temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ)  
  
  temp = - intx3(tempa,temp79b,ri2_79) &
       -   intx3(tempb,temp79a,ri2_79) 

  v3vff = temp
end function v3vff
#endif


! V3chipsipsi
! ===========
function v3chipsipsi(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempd, tempe, tempf
  integer :: j

  ! <f,g>,r
  temp79b = f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DRZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRR) + f(:,OP_DZ )*g(:,OP_DRZ)
  
  ! <f,g>,z
  temp79c = f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRZ) + f(:,OP_DZ )*g(:,OP_DZZ)

  do j=1, dofs_per_element
     ! <e,h>,r
     tempe(j,:) = e(j,:,OP_DRR)*h(:,OP_DR ) + e(j,:,OP_DRZ)*h(:,OP_DZ ) &
          +       e(j,:,OP_DR )*h(:,OP_DRR) + e(j,:,OP_DZ )*h(:,OP_DRZ)
  
     ! <e,h>,z
     tempf(j,:) = e(j,:,OP_DRZ)*h(:,OP_DR ) + e(j,:,OP_DZZ)*h(:,OP_DZ ) &
          +       e(j,:,OP_DR )*h(:,OP_DRZ) + e(j,:,OP_DZ )*h(:,OP_DZZ)
  end do
  
  temp = intx3(tempe,ri6_79,temp79b) &
       + intx3(tempf,ri6_79,temp79c) &
       - intx4(e(:,:,OP_DZ),ri6_79,temp79c,h(:,OP_GS)) &
       - intx4(e(:,:,OP_DR),ri6_79,temp79b,h(:,OP_GS))

  if(itor.eq.1) then
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

     do j=1, dofs_per_element
        tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
     end do

     temp = temp &
          -2.*intx3(tempd,ri7_79,temp79b) &
          -2.*intx3(tempe,ri7_79,temp79a) &
          +2.*intx4(e(:,:,OP_DR),ri7_79,temp79a,h(:,OP_GS)) &
          +4.*intx3(tempd,ri8_79,temp79a)
  endif

  v3chipsipsi = temp
end function v3chipsipsi


! V3chipsib
! =========
function v3chipsib(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb, tempd
  integer :: j

  temp79a = h(:,OP_1)*f(:,OP_GS) &
       + h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
  temp79c = f(:,OP_DZP)*h(:,OP_DR ) - f(:,OP_DRP)*h(:,OP_DZ ) &
       +    f(:,OP_DZ )*h(:,OP_DRP) - f(:,OP_DR )*h(:,OP_DZP)
  if(itor.eq.1) then
     temp79a = temp79a - 2.*ri_79*f(:,OP_DR)*h(:,OP_1)
     temp79c = temp79c - 4.*ri_79* &
          (f(:,OP_DZP)*h(:,OP_1) + f(:,OP_DZ)*h(:,OP_DP))
  endif


  do j=1, dofs_per_element
     tempb(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
     tempd(j,:) = (f(:,OP_DZ )* &
          e(j,:,OP_DR)-f(:,OP_DR )*e(j,:,OP_DZ))*h(:,OP_DP) &
          +    (f(:,OP_DZP)* &
          e(j,:,OP_DR)-f(:,OP_DRP)*e(j,:,OP_DZ))*h(:,OP_1 )
     if(itor.eq.1) then
        tempb(j,:) = tempb(j,:) - 4.*ri_79*e(j,:,OP_DZ)*h(:,OP_1)
     end if
  end do

  temp = intx4(e(:,:,OP_DZ),ri7_79,g(:,OP_DRP),temp79a) &
       - intx4(e(:,:,OP_DR),ri7_79,g(:,OP_DZP),temp79a) &
       - intx4(tempb,ri7_79,f(:,OP_DZP),g(:,OP_DZ)) &
       - intx4(tempb,ri7_79,f(:,OP_DRP),g(:,OP_DR)) &
       - intx4(tempb,ri7_79,f(:,OP_DZ),g(:,OP_DZP)) &
       - intx4(tempb,ri7_79,f(:,OP_DR),g(:,OP_DRP)) &
       + intx4(e(:,:,OP_DZ),ri7_79,temp79c,g(:,OP_DZ)) &
       + intx4(e(:,:,OP_DR),ri7_79,temp79c,g(:,OP_DR)) &
       + intx3(tempd,ri7_79,g(:,OP_GS))

  v3chipsib = temp
#else
  v3chipsib = 0.
#endif
end function v3chipsib


! V3chibb
! =======
function v3chibb(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chibb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb
  integer :: j

  temp79a = g(:,OP_1)*f(:,OP_GS) &
       + g(:,OP_DZ)*f(:,OP_DZ) + g(:,OP_DR)*f(:,OP_DR)
        
  if(itor.eq.1) then
     temp79a = temp79a - &
          2.*ri_79*f(:,OP_DR)*g(:,OP_1)
  end if
        
  temp = intx4(e(:,:,OP_GS),ri6_79,temp79a,h(:,OP_1))

  if(itor.eq.1) then
     temp = temp - &
          2.*intx4(e(:,:,OP_DR),ri7_79,temp79a,h(:,OP_1))
  endif
        
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
  do j=1, dofs_per_element 
     tempb(j,:) = &
          +(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
          *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) &
          +(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
          *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP))
  end do
  temp = temp + intx2(tempb,ri8_79)
#else
  do j=1, dofs_per_element
     tempb(j,:) = &
          (e(j,:,OP_DZ)*f(:,OP_DZPP) + e(j,:,OP_DR)*f(:,OP_DRPP))*g(:,OP_1  ) &
          + 2.*(e(j,:,OP_DZ)*f(:,OP_DZP ) + e(j,:,OP_DR)*f(:,OP_DRP ))*g(:,OP_DP ) &
          +    (e(j,:,OP_DZ)*f(:,OP_DZ  ) + e(j,:,OP_DR)*f(:,OP_DR  ))*g(:,OP_DPP)
  end do
  temp = temp - intx3(tempb,ri8_79,h(:,OP_1))
#endif
#endif

  v3chibb = temp
end function v3chibb

#ifdef USE3D
! V3chipsif
! =====
function v3chipsif(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chipsif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempc, tempd
  integer :: j


  do j=1, dofs_per_element
     ! [nu, f']_R * R
     tempa(j,:) = &
            e(j,:,OP_DRZ)*h(:,OP_DR) - e(j,:,OP_DRR)*h(:,OP_DZ) &
          + e(j,:,OP_DZ)*h(:,OP_DRR) - e(j,:,OP_DR)*h(:,OP_DRZ)
     if(itor.eq.1) then 
        tempa(j,:) = tempa(j,:) - ri_79* &
             (e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)) 
     end if
     ! [nu, f']_Z * R
     tempb(j,:) = &
            e(j,:,OP_DZZ)*h(:,OP_DR) - e(j,:,OP_DRZ)*h(:,OP_DZ) &
          + e(j,:,OP_DZ)*h(:,OP_DRZ) - e(j,:,OP_DR)*h(:,OP_DZZ)
     ! ((nu, psi)/R^2)_R * R^2
     tempc(j,:) = &
            e(j,:,OP_DRZ)*g(:,OP_DZ) + e(j,:,OP_DRR)*g(:,OP_DR) & 
          + e(j,:,OP_DZ)*g(:,OP_DRZ) + e(j,:,OP_DR)*g(:,OP_DRR) 
     if(itor.eq.1) then 
        tempc(j,:) = tempc(j,:) - 2*ri_79* &
             (e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR)) 
     end if
     ! ((nu, psi)/R^2)_Z * R^2
     tempd(j,:) = &
            e(j,:,OP_DZZ)*g(:,OP_DZ) + e(j,:,OP_DRZ)*g(:,OP_DR) & 
          + e(j,:,OP_DZ)*g(:,OP_DZZ) + e(j,:,OP_DR)*g(:,OP_DRZ) 
  end do

  ! [chi, f']_R*R 
  temp79a = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
  if(itor.eq.1) then 
     temp79a = temp79a - ri_79* &        
          (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
  end if
  ! [chi, f']_Z*R 
  temp79b = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
  ! ((chi, psi)/R^2)_R * R^2
  temp79c = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &
       +    f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ) 
  if(itor.eq.1) then 
     temp79c = temp79c - 2*ri_79* & 
          (f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ)) 
  end if
  ! ((chi, psi)/R^2)_Z * R^2
  temp79d = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &
       +    f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ) 

  temp = intx3(tempa,temp79c,ri5_79) &
       + intx3(tempb,temp79d,ri5_79) &
       + intx3(tempc,temp79a,ri5_79) &
       + intx3(tempd,temp79b,ri5_79) &
       - intx4(e(:,:,OP_DZ),temp79b,ri5_79,g(:,OP_GS)) &
       - intx4(e(:,:,OP_DR),temp79a,ri5_79,g(:,OP_GS)) 

  v3chipsif = temp
end function v3chipsif

! V3chibf
! =====
function v3chibf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chibf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j


  do j=1, dofs_per_element
     ! [nu, f']*R
     tempa(j,:) = &
          + e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ) 
     ! (nu, f'')
     tempb(j,:) = &
          + e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
  end do

  ! [chi, F/R^4]'*R^5
  temp79a = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &           
       + f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP) 
  if(itor.eq.1) then 
     temp79a = temp79a - 4*ri_79* &        
          (f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) 
  end if
  ! [chi, f']_R*R 
  temp79b = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
  if(itor.eq.1) then 
     temp79b = temp79b - ri_79* &        
          (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
  end if
  ! [chi, f']_Z*R 
  temp79c = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
       +    f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
  ! (chi, F/R^2)*R^2 + F*chi_GS
  temp79d = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
       +    f(:,OP_GS)*g(:,OP_1) 
  if(itor.eq.1) then 
     temp79d = temp79d - 2*ri_79*f(:,OP_DR)*g(:,OP_1) 
  end if

  temp = intx3(tempa,temp79a,ri6_79) &
       - intx3(tempb,temp79d,ri6_79) &
       - intx4(e(:,:,OP_DZ),temp79b,ri6_79,g(:,OP_DP)) &
       + intx4(e(:,:,OP_DR),temp79c,ri6_79,g(:,OP_DP)) &
       - intx4(e(:,:,OP_DZP),temp79b,ri6_79,g(:,OP_1)) &
       + intx4(e(:,:,OP_DRP),temp79c,ri6_79,g(:,OP_1)) 

  v3chibf = temp
end function v3chibf

! V3chiff
! =====
function v3chiff(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chiff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  integer :: j

  do j=1, dofs_per_element
     ! [nu, f']_R*R
     tempa(j,:) = &
            e(j,:,OP_DRZ)*h(:,OP_DR) - e(j,:,OP_DRR)*h(:,OP_DZ) & 
          + e(j,:,OP_DZ)*h(:,OP_DRR) - e(j,:,OP_DR)*h(:,OP_DRZ) 
     if(itor.eq.1) then 
        tempa(j,:) = tempa(j,:) - ri_79* &
             (e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)) 
     end if
     ! [nu, f']_Z*R
     tempb(j,:) = &
            e(j,:,OP_DZZ)*h(:,OP_DR) - e(j,:,OP_DRZ)*h(:,OP_DZ) & 
          + e(j,:,OP_DZ)*h(:,OP_DRZ) - e(j,:,OP_DR)*h(:,OP_DZZ) 
  end do

  ! [chi, f']_R*R
  temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
  if(itor.eq.1) then 
     temp79a = temp79a - ri_79* &        
          (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
  end if
  ! [chi, f']_Z*R 
  temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
       +    f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            

  temp = intx3(tempa,temp79a,ri4_79) &
       + intx3(tempb,temp79b,ri4_79) 

  v3chiff = temp
end function v3chiff
#endif

! V3uun
! =====
function v3uun(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3uun
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v3uun = 0.
     return
  end if

  temp = intx4(e(:,:,OP_DZ),f(:,OP_DR),g(:,OP_DRZ),h(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_DRR),h(:,OP_1)) &
       + intx4(e(:,:,OP_DR),f(:,OP_DZ),g(:,OP_DRZ),h(:,OP_1)) &
       - intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_DZZ),h(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_LP),h(:,OP_1)) &
       - intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_LP),h(:,OP_1))
  if(itor.eq.1) then
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
     temp = temp + intx4(e(:,:,OP_DZ),ri_79,temp79a,h(:,OP_1))
  endif

  v3uun = temp
end function v3uun


! V3uvn
! =====
function v3uvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3uvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v3uvn = 0.
     return
  end if

  temp = intx5(e(:,:,OP_DZ),ri_79,f(:,OP_DRP),g(:,OP_1),h(:,OP_1)) &
       - intx5(e(:,:,OP_DR),ri_79,f(:,OP_DZP),g(:,OP_1),h(:,OP_1))

#else
  temp = 0.
#endif

  v3uvn = temp
end function v3uvn



! V3vvn
! =====
function v3vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v3vvn = 0.
     return
  end if

  if(itor.eq.0) then
     temp = 0.
  else
     temp = -intx5(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_1),h(:,OP_1))
  endif

  v3vvn = temp
end function v3vvn


! V3uchin
! =======
function v3uchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3uchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v3uchin = 0.
     return
  end if


     temp79a = g(:,OP_DZZ)*f(:,OP_DR) - g(:,OP_DRZ)*f(:,OP_DZ) &
          +    g(:,OP_DZ)*f(:,OP_DRZ) - g(:,OP_DR)*f(:,OP_DZZ)
     temp79b = g(:,OP_DRZ)*f(:,OP_DR) - g(:,OP_DRR)*f(:,OP_DZ) &
          +    g(:,OP_DZ)*f(:,OP_DRR) - g(:,OP_DR)*f(:,OP_DRZ)

     temp = intx4(e(:,:,OP_DZ),ri3_79,temp79a,h(:,OP_1)) &
          + intx4(e(:,:,OP_DR),ri3_79,temp79b,h(:,OP_1)) &
          + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) &
          - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1))

     if(itor.eq.1) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
        temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

        temp = temp + &
             2.*intx4(e(:,:,OP_DZ),ri4_79,temp79a,h(:,OP_1)) &
             +  intx4(e(:,:,OP_DR),ri4_79,temp79b,h(:,OP_1))
     end if
  
  v3uchin = temp
end function v3uchin


! V3vchin
! =======
function v3vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v3vchin = 0.
     return
  end if

  temp = intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1)) &
       + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1))
#else
  temp = 0.
#endif
  v3vchin = temp
end function v3vchin


! V3chichin
! =========
function v3chichin(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chichin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v3chichin = 0.
     return
  end if

  temp79a = f(:,OP_DZ)*g(:,OP_DZZ) + f(:,OP_DR)*g(:,OP_DRZ)
  temp79b = f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DR)*g(:,OP_DRR)

  temp = intx4(e(:,:,OP_DZ),ri6_79,temp79a,h(:,OP_1)) &
       + intx4(e(:,:,OP_DR),ri6_79,temp79b,h(:,OP_1))

  if(itor.eq.1) then
     temp = temp &
          - 2.*intx5(e(:,:,OP_DZ),ri7_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
          - 2.*intx5(e(:,:,OP_DR),ri7_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  endif

  v3chichin = temp
end function v3chichin


! V3ngrav
! =======
function v3ngrav(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ngrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ngrav = 0.
     return
  endif

  temp = gravz*intx2(e(:,:,OP_DZ),f(:,OP_1)) & 
       + gravr*intx3(e(:,:,OP_DR),ri2_79,f(:,OP_1)) 

  v3ngrav = temp
end function v3ngrav


! V3ungrav
! ========
function v3ungrav(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ungrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ungrav = 0.
     return
  endif

  temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
     
  temp = gravz*intx3(e(:,:,OP_DZ), ri_79,temp79a) &
       + gravr*intx3(e(:,:,OP_DR),ri3_79,temp79a)

  v3ungrav = temp
end function v3ungrav


! V3chingrav
! ==========
function v3chingrav(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chingrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3chingrav = 0.
     return
  endif

  temp79a = -(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
       + f(:,OP_LP)*g(:,OP_1))

  temp = gravz*intx2(e(:,:,OP_DZ),temp79a) &
       + gravr*intx3(e(:,:,OP_DR),ri2_79,temp79a)

  v3chingrav = temp
end function v3chingrav


! V3ndenmgrav
! ===========
function v3ndenmgrav(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ndenmgrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  real, intent(in) :: g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ndenmgrav = 0.
     return
  endif

  temp = gravz*intx2(e(:,:,OP_DZ),f(:,OP_LP)) &
       + gravr*intx3(e(:,:,OP_DR),ri2_79,f(:,OP_LP))

  v3ndenmgrav = g*temp
end function v3ndenmgrav


! V3us
! ====
function v3us(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3us
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),temp79a) &
       - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),temp79a)

  v3us = temp
end function v3us


! V3chis
! ======
function v3chis(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),temp79a) &
       + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),temp79a)

  v3chis = temp
end function v3chis


! V3psiforce
! ===
vectype function v3psiforce(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g

  vectype :: temp

  temp = -int4(ri3_79,g(:,OP_1),e(:,OP_DZ),f(:,OP_DR)) &
       +  int4(ri3_79,g(:,OP_1),e(:,OP_DR),f(:,OP_DZ))

  v3psiforce = temp
end function v3psiforce



! V3par
! =====
function v3par(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3par
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  if(itor.eq.1) then
     temp = - intx3(e(:,:,OP_DR),ri3_79,f(:,OP_1))
  endif

  v3par = temp
end function v3par

function v3parb2ipsipsi(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3parb2ipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  temp79a = f(:,OP_1)*g(:,OP_1)*ri4_79
  temp =  temp                                                 &
       +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DR),i(:,OP_DRR))   &
       +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DZ),i(:,OP_DRZ))   &
       +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DRR),i(:,OP_DR))   &
       +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DRZ),i(:,OP_DZ))   &
       +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DR),i(:,OP_DRZ))   &
       +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DZ),i(:,OP_DZZ))   &
       +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DRZ),i(:,OP_DR))   &
       +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DZZ),i(:,OP_DZ))   

  temp79b = f(:,OP_1)*ri4_79
  temp = temp                                                          &
       +  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DR)) &
       -  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ)) &
       -  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DR)) &
       +  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

  temp79c = g(:,OP_1)*ri4_79
  temp = temp                                                         &
       +  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DR)) &
       -  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ)) &
       -  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DR)) &
       +  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

  temp79d = -f(:,OP_1)*g(:,OP_1)*h(:,OP_GS)*ri4_79
  temp = temp                                   &
       +  intx3(e(:,:,OP_DR),temp79d,i(:,OP_DR)) &
       +  intx3(e(:,:,OP_DZ),temp79d,i(:,OP_DZ))   

  v3parb2ipsipsi = temp
end function v3parb2ipsipsi

function v3parb2ipsib(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3parb2ipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i

  vectype,dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = f(:,OP_DP)*g(:,OP_1)*i(:,OP_1)*ri5_79
  temp = intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DR)) &
       - intx3(e(:,:,OP_DR),temp79a,h(:,OP_DZ))
#endif

  v3parb2ipsib = temp
end function v3parb2ipsib

#ifdef USEPARTICLES
function v3p_2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3p_2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        temp = 0.
!!$        temp = &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR)) &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ))
     else
        temp = intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g) &
             + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g)
     end if

  v3p_2 = temp
end function v3p_2

function v3pbb(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3pbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79b = g*ri4_79
        temp = temp                                                         &
             +  intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DR)) &
             -  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DZ)) &
             -  intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DR)) &
             +  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DZ))

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = f(:,OP_DP)*g*bzt79(:,OP_1)*ri5_79
     temp = temp + intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR)) &
          - intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ))
     temp79a = -g*ri3_79
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
        temp79a = -g*ri3_79
        temp = temp +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              -intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DZ)) 
        temp79a = -f(:,OP_DP)*g*bzt79(:,OP_1)*ri4_79
        temp = temp +intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR)) &
              +intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ)) 
        temp79a = g*ri2_79
        temp = temp +intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
#endif
    end if

  v3pbb = temp
end function v3pbb

function v3jxb(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3jxb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79a = -f*pst79(:,OP_GS)*ri4_79
        temp = temp                                   &
             +  intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DR)) &
             +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ))   
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = - f*bzt79(:,OP_1)*ri4_79
     temp = temp + intx3(e(:,:,OP_DR),temp79a,bzt79(:,OP_DR)+bfpt79(:,OP_DRP)) &
          - intx3(e(:,:,OP_DZ),temp79a,bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))
     temp79a = -f*bzt79(:,OP_1)*ri5_79
     temp = temp+intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DZP)) &
          -  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DRP))
     temp79a = -f*pst79(:,OP_GS)*ri3_79
     temp = temp+intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR)) &
          -  intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ))
#endif
    end if

  v3jxb = temp
end function v3jxb


function v3pbb1psi(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3pbb1psi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79b = g*ri4_79
        temp = temp                                                         &
             +  intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DZ),h(:,OP_DR),pst79(:,OP_DR)) &
             -  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DZ),h(:,OP_DR),pst79(:,OP_DZ)) &
             -  intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DR),h(:,OP_DZ),pst79(:,OP_DR)) &
             +  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DR),h(:,OP_DZ),pst79(:,OP_DZ))

#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = -g*ri3_79
        temp = temp +intx5(e(:,:,OP_DR),temp79a,h(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,h(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,h(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              -intx5(e(:,:,OP_DZ),temp79a,h(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DZ)) 
#endif
    end if

  v3pbb1psi = temp
end function v3pbb1psi

function v3pbb1f(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3pbb1f
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = -g*ri3_79
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DR),h(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DZ),h(:,OP_DZ)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DR),h(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DZ),h(:,OP_DZ)) 
        temp79a = g*ri2_79
        temp = temp +intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR),f(:,OP_DR),h(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ),f(:,OP_DR),h(:,OP_DR)) &
              +intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR),f(:,OP_DZ),h(:,OP_DZ)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ),f(:,OP_DZ),h(:,OP_DZ)) 
#endif
    end if

  v3pbb1f = temp
end function v3pbb1f

#endif
!==============================================================================
! B1 TERMS
!==============================================================================


! B1psi
! =====
function b1psi(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

!!$  if(jadv.eq.0) then
!!$     temp = int2(e(:,OP_1),f(:,OP_1))
!!$  else
!!$!     temp79a = e(:,OP_DZ)*f(:,OP_DZ) + e(:,OP_DR)*f(:,OP_DR)
!!$!     temp = -int2(ri2_79,temp79a) 
!!$!...changed 7/22/08    scj
!!$     temp = -int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ)) &
!!$            -int3(ri2_79,e(:,OP_DR),f(:,OP_DR))
!!$  endif

  if(jadv.eq.0) then
     temp = intx2(e(:,:,OP_1),f(:,OP_1))
  else
     temp = intx3(e(:,:,OP_GS),ri2_79,f(:,OP_1))
  endif

  b1psi = temp
end function b1psi


! B1psiu
! ======
function b1psiu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psiu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     ! surface terms
  if(jadv.eq.0) then
     temp = intx4(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_DZ)) &
          - intx4(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_DR))
  else
     temp = intx4(e(:,:,OP_GS),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
          - intx4(e(:,:,OP_GS),ri_79,f(:,OP_DZ),g(:,OP_DR))
  endif
  
  b1psiu = temp
end function b1psiu


! B1psiv
! ======
function b1psiv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psiv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
    temp = 0.
  else
     temp = intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1 )) &
          + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1 )) &
          + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR ),g(:,OP_DP)) &
          + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ ),g(:,OP_DP))
  endif
#else
  temp = 0
#endif

  b1psiv = temp
end function b1psiv


! B1psid
! ======
function b1psid(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b1psid
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(mass_ratio.eq.0. .or. db.eq.0.) then
     b1psid = 0.
     return
  endif

  temp = intx3(e(:,:,OP_1),f(:,OP_GS),g(:,OP_1))

  b1psid = temp*me_mp*mass_ratio*db**2
end function b1psid


! B1bu
! ====
function b1bu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     temp = -(intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZP)) &
          + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DRP)) &
          + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DP),g(:,OP_DZ)) &
          + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DP),g(:,OP_DR)))
  endif

#else
  temp  = 0.
#endif
  b1bu = temp
end function b1bu

! B1jrepsieta
! ====
function b1jrepsieta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jrepsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     temp = -(intx6(e(:,:,OP_DZ),ri3_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
          -intx6(e(:,:,OP_DR),ri3_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)) &
          +intx6(e(:,:,OP_DZ),ri3_79,f(:,OP_DP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
          -intx6(e(:,:,OP_DR),ri3_79,f(:,OP_DP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
  endif

#else
  temp  = 0.
#endif
  b1jrepsieta = temp
end function b1jrepsieta

! B1jrebeta
! ====
function b1jrebeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jrebeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if (jadv .eq. 0) then
     temp = -intx5(e(:,:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_1),i(:,OP_1))
  else
     temp = -intx6(e(:,:,OP_GS),ri2_79,f(:,OP_1),g(:,OP_1),h(:,OP_1),i(:,OP_1))
  endif
  temp = temp*1.000
  b1jrebeta = temp
end function b1jrebeta

! B1jrefeta
! ====
function b1jrefeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jrefeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     temp = (intx6(e(:,:,OP_DZ),ri2_79,f(:,OP_DP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1))&
          +intx6(e(:,:,OP_DR),ri2_79,f(:,OP_DP),g(:,OP_DR),h(:,OP_1),i(:,OP_1))&
          +intx6(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1),i(:,OP_1))&
          +intx6(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)))
  end if

#else
  temp  = 0.
#endif
  b1jrefeta = temp
end function b1jrefeta

! B2jrepsieta
! ====
function b2jrepsieta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2jrepsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = 0.
  else
     temp = -(intx6(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)) &
          +intx6(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DR),h(:,OP_1),i(:,OP_1)))
  endif

  b2jrepsieta = temp
end function b2jrepsieta

! B2jrefeta
! ====
function b2jrefeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2jrefeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     temp = -(intx6(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
          -   intx6(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
  endif

#else
  temp  = 0.
#endif
  b2jrefeta = temp
end function b2jrefeta

! B1bv
! ====
function b1bv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  b1bv = 0.
end function b1bv


! B1psichi
! ========
function b1psichi(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b1psichi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = -intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ)) &
          -  intx4(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR))
  else
     temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZZ),g(:,OP_DZ )) &
          + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_DZZ)) &
          + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DRZ),g(:,OP_DR )) &
          + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DR ),g(:,OP_DRZ)) &
          + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRZ),g(:,OP_DZ )) &
          + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DZ ),g(:,OP_DRZ)) &
          + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRR),g(:,OP_DR )) &
          + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_DRR))
     if(itor.eq.1) then
        temp = temp - 2.* &
             (intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DR),g(:,OP_DR)) &
             +intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DZ),g(:,OP_DZ)))
     endif
  end if

  b1psichi = temp
end function b1psichi


! B1bchi
! ======
function b1bchi(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b1bchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     temp = intx4(e(:,:,OP_DZ),ri5_79,g(:,OP_DRP),f(:,OP_1))  &
          - intx4(e(:,:,OP_DR),ri5_79,g(:,OP_DZP),f(:,OP_1))  &
          + intx4(e(:,:,OP_DZ),ri5_79,g(:,OP_DR),f(:,OP_DP))  &
          - intx4(e(:,:,OP_DR),ri5_79,g(:,OP_DZ),f(:,OP_DP))
  endif

#else
  temp = 0.
#endif
  b1bchi = temp
end function b1bchi


! B1psieta
! ========
function b1psieta1(e,f,g,h,imod)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psieta1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp
  logical, intent(in) :: imod

  if(jadv.eq.0) then

     temp = intx3(e(:,:,OP_1),f(:,OP_GS),g(:,OP_1))
#if defined(USE3D) || defined(USECOMPLEX)
     if(iupstream.eq.1) then 
        temp79a = abs(h(:,OP_1))*magus
        temp = temp + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a)
     elseif(iupstream.eq.2) then
        temp79a = abs(h(:,OP_1))*magus
        temp = temp - intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),temp79a)
     endif
#endif

     if(hypf.ne.0) then
        if(ihypeta.eq.1) then
           temp = temp - hypf* &
                (intx3(e(:,:,OP_1),g(:,OP_LP),f(:,OP_GS)) &
                +intx3(e(:,:,OP_LP),g(:,OP_1),f(:,OP_GS)) &
                +2.*intx3(e(:,:,OP_DZ),g(:,OP_DZ),f(:,OP_GS)) &
                +2.*intx3(e(:,:,OP_DR),g(:,OP_DR),f(:,OP_GS)))
           if(itor.eq.1) temp = temp  &
                - 2.*hypf*intx4(e(:,:,OP_1 ),ri_79 ,f(:,OP_GS),g(:,OP_DR)) &
                - 2.*hypf*intx4(e(:,:,OP_DR),ri_79 ,f(:,OP_GS),g(:,OP_1))
        else
           temp = temp - hypf*intx2(e(:,:,OP_LP),f(:,OP_GS))
           
           if(itor.eq.1) then
              temp = temp - 2.*hypf*intx3(e(:,:,OP_DR),ri_79,f(:,OP_GS))
           endif
        end if
     end if
  else
     temp = intx4(e(:,:,OP_GS),ri2_79,g(:,OP_1),f(:,OP_GS))

#if defined(USE3D) || defined(USECOMPLEX)
     if(.not.imod) then
        temp = temp 
        if(iupstream.eq.1) then   
           temp79a = abs(h(:,OP_1))*magus
           temp = temp - &
                (intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZPP),temp79a) &
                +intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRPP),temp79a))
           if(itor.eq.1) then
              temp = temp + intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DPP),temp79a)
           endif
        elseif(iupstream.eq.2) then
           temp79a = abs(h(:,OP_1))*magus
           temp = temp + &
                (intx4(e(:,:,OP_DZPP),ri6_79,f(:,OP_DZPP),temp79a) &
                +intx4(e(:,:,OP_DRPP),ri6_79,f(:,OP_DRPP),temp79a))
           if(itor.eq.1) then
              temp = temp - intx4(e(:,:,OP_DRPP),ri7_79,f(:,OP_DPP),temp79a)
           endif
        endif
     end if
#endif
  endif

  b1psieta1 = temp
end function b1psieta1


function b1psieta2(e,f,g,h,imod)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psieta2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp
  logical, intent(in) :: imod

  temp = 0.
  if(jadv.ne.0) then
#if defined(USECOMPLEX)
     if(.not.imod) then
        temp = - &
             (intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZPP),g(:,OP_1)) &
             +intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRPP),g(:,OP_1)) &
             +intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_DP)) &
             +intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_DP)))
     end if
#endif
#if defined(USE3D)
     if(.not.imod) then
        temp =  &
             intx4(e(:,:,OP_DZP),ri4_79,f(:,OP_DZP),g(:,OP_1)) &
             +intx4(e(:,:,OP_DRP),ri4_79,f(:,OP_DRP),g(:,OP_1)) 
     end if
#endif
  endif
  b1psieta2 = temp
end function b1psieta2


! B1jeta
! ========
function b1jeta(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  ! note: f = -delstar(psi)

  if(ihypeta.eq.1) then

     temp = hypf*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_GS),g(:,OP_1))

#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + hypf* &
          (intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DP) ,g(:,OP_DP)) &
          +intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DPP),g(:,OP_1)))
#endif
   
  else
     temp = hypf*intx3(e(:,:,OP_GS),ri2_79,f(:,OP_GS))
     
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + hypf* &
          intx3(e(:,:,OP_GS),ri4_79,f(:,OP_DPP))
#endif
  endif

  b1jeta = temp
end function b1jeta


! B1beta
! ======
function b1beta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1beta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else 
#if defined(USE3D)
     temp = - (intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZ),g(:,OP_1 )) &
          - intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DR),g(:,OP_1 )))
#endif
#if defined(USECOMPLEX)
     temp = intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_1 )) &
          - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_1 )) &
          + intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZ ),g(:,OP_DP)) &
          - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DR ),g(:,OP_DP))
#endif
#ifndef USEST
     if(hypf.gt.0. .and. imp_hyper.le.1) then
        if(ihypeta.eq.0) then
           temp = temp - hypf*intx3(e(:,:,OP_DZP),ri5_79,f(:,OP_DRPP)) &
                + hypf*intx3(e(:,:,OP_DRP),ri5_79,f(:,OP_DZPP))
           if(itor.eq.1) then
              temp = temp + 4.*hypf*intx3(e(:,:,OP_DZP),ri6_79,f(:,OP_DPP))
           endif
        else
           temp = temp + hypf* &
                (-intx4(e(:,:,OP_DZP),ri5_79,g(:,OP_1),f(:,OP_DRPP))   &
                + intx4(e(:,:,OP_DRP),ri5_79,g(:,OP_1),f(:,OP_DZPP))   &
                - intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DRP),f(:,OP_GS))   &
                + intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DZP),f(:,OP_GS))   &
                - intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DR),f(:,OP_GSP))   &
                + intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DZ),f(:,OP_GSP)))
           if(itor.eq.1) then
              temp = temp + &
                   4.*hypf*intx4(e(:,:,OP_DZP),ri6_79,g(:,OP_1),f(:,OP_DPP))
           endif

        endif
     endif
#endif !USEST
  endif

#else
  temp = 0.
#endif

  b1beta = temp
end function b1beta

! B1psij
! ======
function b1psij(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psij
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = hypf*(g(:,OP_DRR)+ri2_79*g(:,OP_DPP) +g(:,OP_DZZ))
  if(itor.eq.1) temp79a = temp79a + hypf*ri_79*g(:,OP_DR)
  if     (ihypeta.eq.1) then
     temp79a =     eta79(:,OP_1)*temp79a          &
          + hypf*( eta79(:,OP_DR)*g(:,OP_DR)      &
          + ri2_79*eta79(:,OP_DP)*g(:,OP_DP)      &
          +        eta79(:,OP_DZ)*g(:,OP_DZ))
  else if(ihypeta.eq.2) then
     temp79a = pt79(:,OP_1)*temp79a          &
          + hypf*( pt79(:,OP_DR)*g(:,OP_DR)      &
          + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
          +        pt79(:,OP_DZ)*g(:,OP_DZ))
  else if(ihypeta.gt.2) then
     temp79a = (bharhypeta)**beta*pt79(:,OP_1)*temp79a          &
          + hypf*(bharhypeta)**beta*(pt79(:,OP_DR)*g(:,OP_DR)      &
          + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
          +        pt79(:,OP_DZ)*g(:,OP_DZ))
  endif
  temp = -intx5(e(:,:,OP_DZP),ri3_79,b2i79(:,OP_1),f(:,OP_DR),temp79a)    &
       +  intx5(e(:,:,OP_DRP),ri3_79,b2i79(:,OP_1),f(:,OP_DZ),temp79a)
  
#else
  temp = 0.
#endif

  b1psij = temp
end function b1psij

! B1bj
! ====
function b1bj(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bj
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp
  real :: hypfm


     temp79a = hypf*(g(:,OP_DRR)+g(:,OP_DZZ))
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a + hypf*ri2_79*g(:,OP_DPP)
#endif
     if(itor.eq.1) temp79a = temp79a + hypf*ri_79*g(:,OP_DR)
     if     (ihypeta.eq.1) then             
         temp79a = eta79(:,OP_1)*temp79a          &
          + hypf*(eta79(:,OP_DR)*g(:,OP_DR)      &
                 + eta79(:,OP_DZ)*g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp79a = temp79a + hypf*ri2_79*eta79(:,OP_DP)*g(:,OP_DP)
#endif
     else if(ihypeta.eq.2) then             
         temp79a = pt79(:,OP_1)*temp79a          &
          + hypf*(pt79(:,OP_DR)*g(:,OP_DR)      &
                 + pt79(:,OP_DZ)*g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp79a = temp79a + hypf*ri2_79*pt79(:,OP_DP)*g(:,OP_DP)
#endif
     else if(ihypeta.gt.2) then
         hypfm = hypf*(bharhypeta)**beta
         temp79a = pt79(:,OP_1)*temp79a*(bharhypeta)**beta &
          + hypfm*(pt79(:,OP_DR)*g(:,OP_DR)      &
                 + pt79(:,OP_DZ)*g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp79a = temp79a + hypfm*ri2_79*pt79(:,OP_DP)*g(:,OP_DP)
#endif
     endif

     temp = intx5(e(:,:,OP_GS),ri2_79,b2i79(:,OP_1),f(:,OP_1),temp79a)

  b1bj = temp
end function b1bj

! B1fj
! ======
function b1fj(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fj
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = hypf*(g(:,OP_DRR)+ri2_79*g(:,OP_DPP) +g(:,OP_DZZ))
     if(itor.eq.1) temp79a = temp79a + hypf*ri_79*g(:,OP_DR)
     if     (ihypeta.eq.1) then
        temp79a = eta79(:,OP_1)*temp79a          &
                + hypf*(eta79(:,OP_DR)*g(:,OP_DR)      &
               + ri2_79*eta79(:,OP_DP)*g(:,OP_DP)      &
                      + eta79(:,OP_DZ)*g(:,OP_DZ))
     else if(ihypeta.eq.2) then
        temp79a = pt79(:,OP_1)*temp79a          &
                + hypf*(pt79(:,OP_DR)*g(:,OP_DR)      &
               + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
                      + pt79(:,OP_DZ)*g(:,OP_DZ))
     else if(ihypeta.gt.2) then
        temp79a = pt79(:,OP_1)*temp79a*(bharhypeta)**beta  &
                + hypf*(bharhypeta)**beta*(pt79(:,OP_DR)*g(:,OP_DR)      &
                                   + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
                                          + pt79(:,OP_DZ)*g(:,OP_DZ))
     endif
     temp =  intx5(e(:,:,OP_DRP),ri2_79,b2i79(:,OP_1),f(:,OP_DR),temp79a)    &
            +intx5(e(:,:,OP_DZP),ri2_79,b2i79(:,OP_1),f(:,OP_DZ),temp79a)
#else
  temp = 0.
#endif

  b1fj = temp
end function b1fj



! B1psipsid
! =========
function b1psipsid(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psipsid
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psipsid = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
          + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1))
  else
     temp = intx5(e(:,:,OP_GS),ri4_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
          + intx5(e(:,:,OP_GS),ri4_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
#ifdef USEST
          - intx5(e(:,:,OP_DZP),ri4_79,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_1)) &
          - intx5(e(:,:,OP_DRP),ri4_79,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_1))
#else
          + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_GS ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_GS ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_GSP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_GSP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DP)) &
          + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_DP))
#endif
  endif
  b1psipsid = temp
#else
  b1psipsid = 0.
#endif
end function b1psipsid


! B1psibd1
! ========
function b1psibd1(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibd1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibd1 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
          - intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
  else
     temp = intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
          - intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
  endif

  b1psibd1 = temp
end function b1psibd1

! B1psibd2
! ========
function b1psibd2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibd2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibd2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     temp = &
#ifdef USEST
          - intx5(e(:,:,OP_DRP),ri5_79,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_1)) &
          + intx5(e(:,:,OP_DZP),ri5_79,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_1))
#else
          + intx5(e(:,:,OP_DR),ri5_79,f(:,OP_DZPP),g(:,OP_1 ),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DZ),ri5_79,f(:,OP_DRPP),g(:,OP_1 ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),ri5_79,f(:,OP_DZP ),g(:,OP_DP),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP ),g(:,OP_DP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),ri5_79,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_DP)) &
          - intx5(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_DP))
#endif
  endif
#else
  temp = 0.
#endif

  b1psibd2 = temp
end function b1psibd2



! B1psifd1
! ========
function b1psifd1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifd1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifd1 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR ),h(:,OP_1)) &
          - intx5(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ ),h(:,OP_1))
  else
     temp = intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DZP),g(:,OP_DR ),h(:,OP_1)) &
          - intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DRP),g(:,OP_DZ ),h(:,OP_1)) &
#ifdef USEST
          - intx5(e(:,:,OP_DZP),ri3_79,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_1)) &
          + intx5(e(:,:,OP_DRP),ri3_79,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_1))
#else
          + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GSP),g(:,OP_DR ),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GSP),g(:,OP_DZ ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS ),g(:,OP_DRP),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS ),g(:,OP_DZP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_DP)) &
          - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_DP))
#endif
  endif
  b1psifd1 = temp
#else
  b1psifd1 = 0.
#endif
end function b1psifd1

! B1psifd2
! ========
function b1psifd2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifd2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifd2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ ),g(:,OP_DRP),h(:,OP_1)) &
          - intx5(e(:,:,OP_1),ri_79,f(:,OP_DR ),g(:,OP_DZP),h(:,OP_1))
  else
     temp = intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DZ ),g(:,OP_DRP),h(:,OP_1)) &
          - intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DR ),g(:,OP_DZP),h(:,OP_1))
  endif
  b1psifd2 = temp
#else
  b1psifd2 = 0.
#endif
end function b1psifd2




! B1bbd
! =====
function b1bbd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bbd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bbd = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
     temp = intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_1 ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_1 ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_DP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_DP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_1 ),h(:,OP_DP)) &
          + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_1 ),h(:,OP_DP))
  endif
  b1bbd = temp
#else
  b1bbd = 0.
#endif
end function b1bbd


! B1bfd1
! ======
function b1bfd1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfd1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfd1 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = intx4(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          + intx4(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  else
     temp = intx5(e(:,:,OP_GS),ri2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          + intx5(e(:,:,OP_GS),ri2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  endif
  b1bfd1 = temp
#else
  b1bfd1 = 0.
#endif
end function b1bfd1

! B1bfd2
! ======
function b1bfd2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfd2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfd2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
#ifdef USECOMPLEX
     temp = intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_DP),h(:,OP_1 ))  &
          + intx5(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_DP),h(:,OP_1 ))  &
          + intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_1 ),h(:,OP_DP))  &
          + intx5(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_1 ),h(:,OP_DP))

     ! f''' term hack
     temp = temp + rfac* &
          (intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
          +intx5(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#elif defined(USE3D)
        ! here, we can integrate by parts
     temp = - &
          (intx5(e(:,:,OP_DZP),ri4_79,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
          +intx5(e(:,:,OP_DRP),ri4_79,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#endif
  endif
  b1bfd2 = temp
#else
  b1bfd2 = 0.
#endif
end function b1bfd2



! B1ped
! =====
function b1ped(e,f,g)
  use basic
  use m3dc1_nint

  implicit none
  vectype, dimension(dofs_per_element) :: b1ped
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1ped = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = intx3(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1))
  else
     temp = intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),g(:,OP_1)) &
          + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1 )) &
          + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1 )) &
          + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ ),g(:,OP_DP)) &
          + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR ),g(:,OP_DP))
  endif
  b1ped = temp
#else
  b1ped = 0.
#endif
end function b1ped

! B1psipsin
! =========
function b1psipsin(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psipsin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psipsin = 0.
     return
  end if

  if(jadv.eq.0) then
     temp79a = ri2_79*h(:,OP_1)*ni79(:,OP_1)**2
     temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),g(:,OP_DZP)) &
          + intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),g(:,OP_DRP))
  else
     temp79a = ri4_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
          + intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
#ifdef USEST
          - intx5(e(:,:,OP_DZP),temp79a,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_1)) &
          - intx5(e(:,:,OP_DRP),temp79a,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_1)) &
          - 2*intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DP)) &
          - 2*intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_DP))
#else
          + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_GS ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_GS ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_GSP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_GSP),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DP)) &
          - intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_DP))
#endif
  endif
  b1psipsin = temp
#else
  b1psipsin = 0.
#endif
end function b1psipsin


! B1psibn1
! ========
function b1psibn1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibn1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibn1 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp79a = ri_79*h(:,OP_1)*ni79(:,OP_1)**2
     temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),g(:,OP_DR)) &
          - intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),g(:,OP_DZ))
  else
     temp79a = ri3_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
          - intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
  endif

  b1psibn1 = temp
end function b1psibn1

! B1psibd2
! ========
function b1psibn2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibn2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibn2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
     temp79a = ri5_79*ni79(:,OP_1)**2
     temp = &
#ifdef USEST
          - intx5(e(:,:,OP_DRP),temp79a,f(:,OP_DZP ),g(:,OP_1),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZP),temp79a,f(:,OP_DRP ),g(:,OP_1),h(:,OP_1 )) &
          - 2*intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_DP)) &
          + 2*intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_DP))
#else
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZPP),g(:,OP_1 ),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRPP),g(:,OP_1 ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZP ),g(:,OP_DP),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRP ),g(:,OP_DP),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_DP)) &
          + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_DP))
#endif
  endif
  b1psibn2 = temp
#else
  b1psibn2 = 0.
#endif
end function b1psibn2



! B1psifn1
! ========
function b1psifn1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifn1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifn1 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp79a =ri_79*h(:,OP_1)*ni79(:,OP_1)**2
     temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZP),g(:,OP_DR )) &
          - intx4(e(:,:,OP_1),temp79a,f(:,OP_DRP),g(:,OP_DZ ))
  else
     temp79a = ri3_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZP),g(:,OP_DR ),h(:,OP_1)) &
          - intx5(e(:,:,OP_GS),temp79a,f(:,OP_DRP),g(:,OP_DZ ),h(:,OP_1)) &
#ifdef USEST
          - intx5(e(:,:,OP_DZP),temp79a,f(:,OP_GS ),g(:,OP_DR),h(:,OP_1 ))&
          + intx5(e(:,:,OP_DRP),temp79a,f(:,OP_GS ),g(:,OP_DZ),h(:,OP_1 ))&
          - 2.*intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_DP))&
          + 2.*intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_DP))
#else
          + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GSP),g(:,OP_DR ),h(:,OP_1 ))&
          - intx5(e(:,:,OP_DR),temp79a,f(:,OP_GSP),g(:,OP_DZ ),h(:,OP_1 ))&
          + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS ),g(:,OP_DRP),h(:,OP_1 ))&
          - intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS ),g(:,OP_DZP),h(:,OP_1 ))&
          - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_DP))&
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_DP))
#endif
  endif
  b1psifn1 = temp
#else
  b1psifn1 = 0.
#endif
end function b1psifn1

! B1psifn2
! ========
function b1psifn2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifn2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifn2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp79a = ri_79*h(:,OP_1)*ni79(:,OP_1)**2
     temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ ),g(:,OP_DRP)) &
          - intx4(e(:,:,OP_1),temp79a,f(:,OP_DR ),g(:,OP_DZP))
  else
     temp79a = ri3_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ ),g(:,OP_DRP),h(:,OP_1)) &
          - intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR ),g(:,OP_DZP),h(:,OP_1))
  endif
  b1psifn2 = temp
#else
  b1psifn2 = 0.
#endif
end function b1psifn2




! B1bbn
! =====
function b1bbn(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bbn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bbn = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
     temp79a = ri4_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_1 ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_1 ),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_DP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_DP),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_1 ),h(:,OP_DP)) &
          - intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_1 ),h(:,OP_DP))
  endif
  b1bbn = temp
#else
  b1bbn = 0.
#endif
end function b1bbn


! B1bfn1
! ======
function b1bfn1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfn1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfn1 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp79a = h(:,OP_1)*ni79(:,OP_1)**2
     temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),g(:,OP_DZ)) &
          + intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),g(:,OP_DR))
  else
     temp79a = ri2_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          + intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  endif
  b1bfn1 = temp
#else
  b1bfn1 = 0.
#endif
end function b1bfn1

! B1bfn2
! ======
function b1bfn2(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfn2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfn2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
#ifdef USECOMPLEX
     temp79a = ri4_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_DZ),temp79a,g(:,OP_DZP),f(:,OP_DP),h(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),temp79a,g(:,OP_DRP),f(:,OP_DP),h(:,OP_1 )) &
          - intx5(e(:,:,OP_DZ),temp79a,g(:,OP_DZP),f(:,OP_1 ),h(:,OP_DP)) &
          - intx5(e(:,:,OP_DR),temp79a,g(:,OP_DRP),f(:,OP_1 ),h(:,OP_DP))

        ! f''' term hack
     temp = temp + rfac* &
          (intx5(e(:,:,OP_DZ),temp79a,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
          +intx5(e(:,:,OP_DR),temp79a,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#elif defined(USE3D)
     temp79a = ri4_79*ni79(:,OP_1)**2
        ! here, we can integrate by parts
     temp = - &
          (intx5(e(:,:,OP_DZP),temp79a,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
          +intx5(e(:,:,OP_DRP),temp79a,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#endif
  endif
  b1bfn2 = temp
#else
  b1bfn2 = 0.
#endif
end function b1bfn2



! B1pen
! =====
function b1pen(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1pen
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1pen = 0.
     return
  end if

  if(jadv.eq.0) then
     temp79a = ni79(:,OP_1)**2
     temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DP),g(:,OP_1))
  else
     temp79a = ri2_79*ni79(:,OP_1)**2
     temp = intx4(e(:,:,OP_GS),temp79a,f(:,OP_DP),g(:,OP_1)) &
          + intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_1 )) &
          + intx4(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_1 )) &
          - intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_DP)) &
          - intx4(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_DP))
  endif
  b1pen = temp
#else
  b1pen = 0.
#endif
end function b1pen



! B1feta
! ======
function b1feta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1feta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = 0.
  else
#ifdef USECOMPLEX
     temp = intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZPP),g(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRPP),g(:,OP_1))
#else
     temp = -intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZP),g(:,OP_1)) &
          +  intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DRP),g(:,OP_1))
#endif
  end if

  b1feta = temp
#else
  b1feta = 0.
#endif
end function b1feta


! b1fu
! ====
function b1fu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = &
          - intx4(e(:,:,OP_1),r2_79,f(:,OP_DZ),g(:,OP_DZ)) &
          - intx4(e(:,:,OP_1),r2_79,f(:,OP_DR),g(:,OP_DR))
  else
     temp = &
          - intx3(e(:,:,OP_GS),f(:,OP_DZ),g(:,OP_DZ)) &
          - intx3(e(:,:,OP_GS),f(:,OP_DR),g(:,OP_DR))
  endif

  b1fu = temp
#else
  b1fu = 0.
#endif
end function b1fu


! b1fv
! ====
function b1fv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = 0.
  else
     temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DRP),g(:,OP_1)) &
          - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZP),g(:,OP_1)) &
          + intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR ),g(:,OP_DP)) &
          - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ ),g(:,OP_DP))
  endif

  b1fv = temp
#else
  b1fv = 0.
#endif
end function b1fv


! b1fchi
! ======
function b1fchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR)) &
          - intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ))
  else
     temp = intx4(e(:,:,OP_GS),ri3_79,f(:,OP_DZ),g(:,OP_DR)) &
          - intx4(e(:,:,OP_GS),ri3_79,f(:,OP_DR),g(:,OP_DZ))
  endif

  b1fchi = temp
#else
  b1fchi = 0.
#endif
end function b1fchi


! B1e
! ===
function b1e(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1e
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.1) then
     temp = 0.
  else
     temp = -intx2(e(:,:,OP_1),f(:,OP_DP))
  endif
  b1e = temp
#else
  b1e = 0.
#endif
end function b1e

function b1vzdot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1vzdot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = intx3(e(:,:,OP_1),r2_79,f(:,OP_1))
  else
     temp = -intx2(e(:,:,OP_DR),f(:,OP_DR))    &
          -intx2(e(:,:,OP_DZ),f(:,OP_DZ))
     if(itor.eq.1) temp = temp - 2.*intx3(e(:,:,OP_DR),ri_79,f(:,OP_1))
  endif

  b1vzdot = temp
end function b1vzdot

function b1chidot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1chidot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = 0.
  else
     temp = -intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRP))    &
          -intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP))   
  endif

  b1chidot = temp
#else
  b1chidot = 0.
#endif
end function b1chidot


function b1psi2bfpe(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psi2bfpe
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j
  vectype, dimension(dofs_per_element) :: temp
  temp = 0.

  temp79a = ri_79*(j(:,OP_DZ)*f(:,OP_DR) - j(:,OP_DR)*f(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = temp79a - (j(:,OP_DR)*i(:,OP_DR) + j(:,OP_DZ)*i(:,OP_DZ))   &
                    + ri2_79*h(:,OP_1)*j(:,OP_DP)
#endif
  temp79a = temp79a*b2i79(:,OP_1)*ri2_79*ni79(:,OP_1)


  if(jadv.eq.0) then
     temp = 0.
  else
     temp = intx3(e(:,:,OP_GS),temp79a,h(:,OP_1))
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp &
          - intx4(e(:,:,OP_DZP),ri_79,temp79a,g(:,OP_DR))   &
          + intx4(e(:,:,OP_DRP),ri_79,temp79a,g(:,OP_DZ))   &
          + intx3(e(:,:,OP_DRP),temp79a,i(:,OP_DR))        &
          + intx3(e(:,:,OP_DZP),temp79a,i(:,OP_DZ))
#endif
  end if

  b1psi2bfpe = temp
end function b1psi2bfpe


!==============================================================================
! B2 TERMS
!==============================================================================

! B2b
! ===
function b2b(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2b
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = intx3(e(:,:,OP_1),ri2_79,f(:,OP_1))

  b2b = temp
end function b2b


! B2psieta
! ========
function b2psieta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_1))

#ifndef USEST
  if(hypi.ne.0 .and. imp_hyper.le.1) then
     if(ihypeta.eq.0) then          
        temp = temp + 2.*hypi* &
             (intx3(e(:,:,OP_DZZ),ri3_79,f(:,OP_DRZP)) &
             -intx3(e(:,:,OP_DRR),ri3_79,f(:,OP_DRZP)) &
             -intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZZP)) &
             +intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DRRP)))
           
        if(itor.eq.1) then
           temp = temp - 2.*hypi* &
                (   intx3(e(:,:,OP_DZZ),ri4_79,f(:,OP_DZP)) &
                -   intx3(e(:,:,OP_DRR),ri4_79,f(:,OP_DZP)) &
                +   intx3(e(:,:,OP_DR),ri5_79,f(:,OP_DZP)) &
                +2.*intx3(e(:,:,OP_DRZ),ri4_79,f(:,OP_DRP)) &
                -4.*intx3(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP)) &
                -   intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRZP)) &
                +   intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZZP)) &
                +2.*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DRRP)))
        endif
     endif
  end if
#endif !USEST

  b2psieta = temp
#else
  b2psieta = 0.
#endif
end function b2psieta


! B2psimue
! ========
vectype function b2psimue(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
  
  temp = &
       - int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
       - int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR)) &
       - int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_1 ))

  b2psimue = temp
end function b2psimue



! B2beta
! ======
function b2beta(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2beta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = &
       - intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
       - intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_1)) 
#if defined(USE3D) || defined(USECOMPLEX)
  if(iupstream.eq.1) then    
     temp79a = abs(h(:,OP_1))*magus
     temp = temp + intx4(e(:,:,OP_1),ri4_79,f(:,OP_DPP),temp79a)
  elseif(iupstream.eq.2) then
     temp79a = abs(h(:,OP_1))*magus
     temp = temp - intx4(e(:,:,OP_DPP),ri6_79,f(:,OP_DPP),temp79a)
  endif
#endif     

  if(imp_hyper.le.1) then
!    the following coding should be checked.  It does not agree with my derivation  scj 4/30/14
     if(hypi.ne.0.) then
        if(ihypeta.eq.1) then
           temp = temp - hypi* &
                (intx4(e(:,:,OP_GS),ri2_79,f(:,OP_GS),g(:,OP_1)) &
                +intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),g(:,OP_DZ)) &
                +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_GS),g(:,OP_DR)))
        else
           temp = temp + hypi* &
                (-intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DZZ)) &
                + intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DZZ)) &
                + intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DRR)) &
                - intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DRR)) &
                - 4.*intx3(e(:,:,OP_DRZ),ri2_79,f(:,OP_DRZ)))

           if(itor.eq.1) then
              temp = temp + hypi*&
                   (-intx3(e(:,:,OP_DZZ),ri3_79,f(:,OP_DR)) &
                   + intx3(e(:,:,OP_DRR),ri3_79,f(:,OP_DR)) &
                   - intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DR)) &
                   - 2.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZ)) &
                   + 4.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZ)) &
                   - 4.*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ)) &
                   - intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DZZ)) &
                   + intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DRR)) &
                   + 2.*intx3(e(:,:,OP_DZ),ri3_79,f(:,OP_DRZ)))
           endif
           
#if defined(USE3D) || defined(USECOMPLEX)
           temp = temp &
#ifdef USEST
                - hypi*intx3(e(:,:,OP_DZP),ri4_79,f(:,OP_DZP)) &
                - hypi*intx3(e(:,:,OP_DRP),ri4_79,f(:,OP_DRP))
#else
                + hypi*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZPP)) &
                + hypi*intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRPP))
#endif
#endif
        endif
     endif
  end if
       
  b2beta = temp
end function b2beta


! B2feta
! ======
function b2feta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2feta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  temp = &
       - intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
       - intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1))

#ifndef USEST
  if(imp_hyper.le.1) then

!   the following coding should be checked.  does not agree with my derivation scj 4/30/2014
     if(hypi.ne.0.) then
        if(ihypeta.eq.0) then
           temp = temp + hypi*&
                (-intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DZZP)) &
                + intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DZZP)) &
                + intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DRRP)) &
                - intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DRRP)) &
                - 4.*intx3(e(:,:,OP_DRZ),ri2_79,f(:,OP_DRZP)))

           if(itor.eq.1) then
              temp = temp + hypi* &
                   (-intx3(e(:,:,OP_DZZ),ri3_79,f(:,OP_DRP)) &
                   + intx3(e(:,:,OP_DRR),ri3_79,f(:,OP_DRP)) &
                   - intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRP)) &
                   - 2.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZP)) &
                   + 4.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZP)) &
                   - 4.*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP)) &
                   - intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DZZP)) &
                   + intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DRRP)) &
                   + 2.*intx3(e(:,:,OP_DZ),ri3_79,f(:,OP_DRZP)))
           endif
        endif
     endif
  end if
#endif !USEST

  b2feta = temp
#else
  b2feta = 0.
#endif
end function b2feta

! B2FJ
! ====
function b2fj(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2fj
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  real :: hypfm

  temp = 0.
  if     (ihypeta.eq.1) then
     temp79a = hypf*eta79(:,OP_1)*g(:,OP_DR)
     temp79b = hypf*eta79(:,OP_1)*g(:,OP_DZ)
     temp79c = hypf*(eta79(:,OP_DP)*g(:,OP_DP) + eta79(:,OP_1)*g(:,OP_DPP))
  else if(ihypeta.eq.2) then
     temp79a = hypf*pt79(:,OP_1)*g(:,OP_DR)
     temp79b = hypf*pt79(:,OP_1)*g(:,OP_DZ)
     temp79c = hypf*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
  else if(ihypeta.gt.2) then
     hypfm = hypf*(bharhypeta)**beta
     temp79a = hypfm*pt79(:,OP_1)*g(:,OP_DR)
     temp79b = hypfm*pt79(:,OP_1)*g(:,OP_DZ)
     temp79c = hypfm*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
  else
     temp79a = hypf*g(:,OP_DR)
     temp79b = hypf*g(:,OP_DZ)
     temp79c = hypf*g(:,OP_DPP)
  endif
     
  temp = -intx5(e(:,:,OP_DZ),b2i79(:,OP_DR),ri_79,f(:,OP_DR ),temp79a) &
       +intx5(e(:,:,OP_DR ),b2i79(:,OP_DR),ri_79,f(:,OP_DZ ),temp79a) &
       -intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri_79,f(:,OP_DR ),temp79a) &
       +intx5(e(:,:,OP_DRR),b2i79(:,OP_1 ),ri_79,f(:,OP_DZ ),temp79a) &
       -intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri_79,f(:,OP_DRR),temp79a) &
       +intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri_79,f(:,OP_DRZ),temp79a) &
       -intx5(e(:,:,OP_DZ ),b2i79(:,OP_DZ),ri_79,f(:,OP_DR ),temp79b) &
       +intx5(e(:,:,OP_DR ),b2i79(:,OP_DZ),ri_79,f(:,OP_DZ ),temp79b) &
       -intx5(e(:,:,OP_DZZ),b2i79(:,OP_1 ),ri_79,f(:,OP_DR ),temp79b) &
       +intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri_79,f(:,OP_DZ ),temp79b) &
       -intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri_79,f(:,OP_DRZ),temp79b) &
       +intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri_79,f(:,OP_DZZ),temp79b) 
  if(itor.eq.1) temp = temp &
       +intx5(e(:,:,OP_DZ ),b2i79(:,OP_1),ri2_79,f(:,OP_DR ),temp79a) &
       -intx5(e(:,:,OP_DR ),b2i79(:,OP_1),ri2_79,f(:,OP_DZ ),temp79a) 
  temp = temp                                                         &
       +intx5(e(:,:,OP_DZ ),b2i79(:,OP_1),ri3_79,f(:,OP_DR ),temp79c) &
       -intx5(e(:,:,OP_DR ),b2i79(:,OP_1),ri3_79,f(:,OP_DZ ),temp79c) 

  b2fj = temp
#else
  b2fj = 0.
#endif
end function b2fj

! B2PSIJ
! ======
function b2psij(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psij
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp
  real :: hypfm

  if     (ihypeta.eq.1) then
     temp79a = hypf*eta79(:,OP_1)*g(:,OP_DR)
     temp79b = hypf*eta79(:,OP_1)*g(:,OP_DZ)
  else if(ihypeta.eq.2) then
     temp79a = hypf*pt79(:,OP_1)*g(:,OP_DR)
     temp79b = hypf*pt79(:,OP_1)*g(:,OP_DZ)
  else if(ihypeta.gt.2) then
     hypfm = hypf*(bharhypeta)**beta
     temp79a = hypfm*pt79(:,OP_1)*g(:,OP_DR)
     temp79b = hypfm*pt79(:,OP_1)*g(:,OP_DZ)
  else
     temp79a = hypf*g(:,OP_DR)
     temp79b = hypf*g(:,OP_DZ)
  endif

  temp = -intx5(e(:,:,OP_DR ),b2i79(:,OP_DR),ri2_79,f(:,OP_DR ),temp79a) &
       -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_DR),ri2_79,f(:,OP_DZ ),temp79a) &
       -  intx5(e(:,:,OP_DRR),b2i79(:,OP_1 ),ri2_79,f(:,OP_DR ),temp79a) &
       -  intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DZ ),temp79a) &
       -  intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DRR),temp79a) &
       -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DRZ),temp79a) &
       -  intx5(e(:,:,OP_DR ),b2i79(:,OP_DZ),ri2_79,f(:,OP_DR ),temp79b) &
       -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_DZ),ri2_79,f(:,OP_DZ ),temp79b) &
       -  intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DR ),temp79b) &
       -  intx5(e(:,:,OP_DZZ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DZ ),temp79b) &
       -  intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DRZ),temp79b) &
       -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DZZ),temp79b) 
  if(itor.eq.1) temp = temp &
       +2.*intx5(e(:,:,OP_DR ),b2i79(:,OP_1),ri3_79,f(:,OP_DR ),temp79a) &
       +2.*intx5(e(:,:,OP_DZ ),b2i79(:,OP_1),ri3_79,f(:,OP_DZ ),temp79a) 
#if defined(USE3D) || defined(USECOMPLEX)
  if     (ihypeta.eq.1) then
     temp79c = hypf*(eta79(:,OP_DP)*g(:,OP_DP) + eta79(:,OP_1)*g(:,OP_DPP))
  else if(ihypeta.eq.2) then
     temp79c = hypf*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
  else if(ihypeta.gt.2) then
     temp79c = hypfm*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
  else
     temp79c = hypf*g(:,OP_DPP)
  endif
  temp = temp                                                        &
       +intx5(e(:,:,OP_DR),b2i79(:,OP_1),ri4_79,f(:,OP_DR ),temp79c) &
       +intx5(e(:,:,OP_DZ),b2i79(:,OP_1),ri4_79,f(:,OP_DZ ),temp79c)
#endif

  b2psij = temp
end function b2psij

! B2bu
! ====
function b2bu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2bu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZ))

  b2bu = temp
end function b2bu


! B2bchi
! ======
function b2bchi(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2bchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_1),g(:,OP_DZ)) &
       + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_1),g(:,OP_DR))        

  b2bchi = temp
end function b2bchi


! B2bd
! ====
function b2bd(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2bd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(mass_ratio.eq.0. .or. db.eq.0.) then
     b2bd = 0.
     return
  endif

  temp = - &
       (intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
       +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_1)))

  b2bd = temp*me_mp*mass_ratio*db**2
end function b2bd



! B2psiv
! ======
function b2psiv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psiv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1))

  b2psiv = temp
end function b2psiv


! B2fv
! ====
function b2fv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2fv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
       + intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1))
#else
  temp = 0.
#endif

  b2fv = temp
end function b2fv


! B2psipsid
! =========
function b2psipsid(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psipsid
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psipsid = 0.
     return
  end if

  temp = intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
       - intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))

  b2psipsid = temp
end function b2psipsid


! B2psibd
! =======
function b2psibd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psibd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psibd = 0.
     return
  end if

  temp = &
       -(intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
       +intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_1),h(:,OP_1)))

  b2psibd = temp
#else
  b2psibd = 0.
#endif
end function b2psibd


! B2bbd
! =====
function b2bbd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2bbd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bbd = 0.
     return
  end if

  temp = intx5(e(:,:,OP_DR),ri3_79,f(:,OP_DZ),g(:,OP_1),h(:,OP_1)) &
       - intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_DR),g(:,OP_1),h(:,OP_1))
  
  b2bbd = temp
end function b2bbd



! B2ped
! =====
function b2ped(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2ped
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2ped = 0.
     return
  end if

  temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1))

  b2ped = temp
end function b2ped


! B2psifd
! =======
function b2psifd(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2psifd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psifd = 0.
     return
  end if
  
  temp = intx5(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
       + intx5(e(:,:,OP_DR),ri2_79,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))

  b2psifd = temp
#else
  b2psifd = 0.
#endif
end function b2psifd


! B2bfd
! =====
function b2bfd(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2bfd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bfd = 0.
     return
  end if
  
  temp = - &
       (intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1)) &
       -intx5(e(:,:,OP_DR),ri3_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1)))
  
  b2bfd = temp
#else
  b2bfd = 0.
#endif
end function b2bfd

! B2psipsin
! =========
function b2psipsin(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psipsin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psipsin = 0.
     return
  end if

  temp79a = ri3_79*ni79(:,OP_1)**2
  temp = intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
       - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))

  b2psipsin = temp
end function b2psipsin


! B2psibn
! =======
function b2psibn(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psibn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psibn = 0.
     return
  end if

  temp79a = ri4_79*ni79(:,OP_1)**2
  temp = &
       -(intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
       +intx5(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_1),h(:,OP_1)))

  b2psibn = temp
#else
  b2psibn = 0.
#endif
end function b2psibn


! B2bbn
! =====
function b2bbn(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2bbn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bbn = 0.
     return
  end if

  temp79a = ri3_79*ni79(:,OP_1)**2
  temp = intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZ),g(:,OP_1),h(:,OP_1)) &
       - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DR),g(:,OP_1),h(:,OP_1))
  
  b2bbn = temp
end function b2bbn



! B2pen
! =====
function b2pen(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2pen
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2pen = 0.
     return
  end if

  temp79a = ri_79*ni79(:,OP_1)**2
  temp = intx4(e(:,:,OP_DR),temp79a,f(:,OP_DZ),g(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DR),g(:,OP_1))

  b2pen = temp
end function b2pen


! B2psifn
! =======
function b2psifn(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2psifn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psifn = 0.
     return
  end if
  
  temp79a = ri2_79*ni79(:,OP_1)**2
  temp = intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
       + intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))

  b2psifn = temp
#else
  b2psifn = 0.
#endif
end function b2psifn


! B2bfn
! =====
function b2bfn(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2bfn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bfn = 0.
     return
  end if
  
  temp79a = ri3_79*ni79(:,OP_1)**2
  temp = - &
       (intx5(e(:,:,OP_DZ),temp79a,f(:,OP_1),g(:,OP_DRP),h(:,OP_1)) &
       -intx5(e(:,:,OP_DR),temp79a,f(:,OP_1),g(:,OP_DZP),h(:,OP_1)))
  
  b2bfn = temp
#else
  b2bfn = 0.
#endif
end function b2bfn

function b2phidot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2phidot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = intx2(e(:,:,OP_DR),f(:,OP_DR))   &
       + intx2(e(:,:,OP_DZ),f(:,OP_DZ))

  b2phidot = temp
end function b2phidot

function b2chidot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2chidot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  if(itor.eq.1) temp = 2.*intx3(e(:,:,OP_1),ri4_79,f(:,OP_DZ))

  b2chidot = temp
end function b2chidot

function b2psi2bfpe(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psi2bfpe
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j
  vectype, dimension(dofs_per_element) :: temp
  temp = 0.

  temp79a = ri_79*(j(:,OP_DZ)*f(:,OP_DR) - j(:,OP_DR)*f(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = temp79a - (j(:,OP_DR)*i(:,OP_DR) + j(:,OP_DZ)*i(:,OP_DZ))   &
                    + ri2_79*h(:,OP_1)*j(:,OP_DP)
#endif
  temp79a = temp79a*b2i79(:,OP_1)*ni79(:,OP_1)


  temp = intx4(e(:,:,OP_DR),ri2_79,temp79a,g(:,OP_DR))    &
       + intx4(e(:,:,OP_DZ),ri2_79,temp79a,g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp + intx4(e(:,:,OP_DZ),ri_79,temp79a,i(:,OP_DR))   &
       - intx4(e(:,:,OP_DR),ri_79,temp79a,i(:,OP_DZ))
#endif

  b2psi2bfpe = temp
end function b2psi2bfpe
!=============================================================================
! B3 TERMS
!=============================================================================

! B3pe
! ====
function b3pe(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pe
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = intx2(e(:,:,OP_1),f(:,OP_1))

  b3pe = temp
end function b3pe


! B3pe27
! ======
function b3pe27(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pe27
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp
  real, dimension(MAX_PTS) :: r

  r = sqrt((x_79-xmag)**2 + (z_79-zmag)**2)
  temp79b = 1. + tanh((r-libetap)/p1)

  temp = intx3(e(:,:,OP_1),f(:,OP_1),temp79b)

  b3pe27 = temp
end function b3pe27


! B3q
! ===
function b3q(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3q
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = intx2(e(:,:,OP_1),f(:,OP_1))

  b3q = temp
end function b3q

! B3psipsieta
! ===========
function b3psipsieta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3psipsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     temp = 0.
  else
     temp = (gam-1.)* &
           intx5(e(:,:,OP_1),ri2_79,f(:,OP_GS), g(:,OP_GS), h(:,OP_1))   
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + (gam-1)*   &
           (intx5(e(:,:,OP_1),ri4_79,f(:,OP_DRP),g(:,OP_DRP),h(:,OP_1))   &
         +  intx5(e(:,:,OP_1),ri4_79,f(:,OP_DZP),g(:,OP_DZP),h(:,OP_1)))
     if(irunaway .gt. 2) then
        temp = temp + 1.*(gam-1.) * &
              (-intx6(e(:,:,OP_1),ri3_79,f(:,OP_DZ),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
              + intx6(e(:,:,OP_1),ri3_79,f(:,OP_DR),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)))
     endif
#endif
  end if

  b3psipsieta = temp
end function b3psipsieta


! B3psibeta
! ===========
function b3psibeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3psibeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     temp = 0.
  else
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (intx5(e(:,:,OP_1),ri3_79,f(:,OP_DZP),g(:,OP_DR),h(:,OP_1))  &
          -intx5(e(:,:,OP_1),ri3_79,f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1)))
#endif
  end if

  if(irunaway .gt. 2) then
     temp = temp + 1.*(gam-1.) * &
                   (-intx6(e(:,:,OP_1),ri2_79,f(:,OP_GS),g(:,OP_1),h(:,OP_1),i(:,OP_1)) &
                  + intx6(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
                  + intx6(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
  end if

  b3psibeta = temp
end function b3psibeta


! B3psifeta
! =========
function b3psifeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3psifeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     temp = 0.
  else
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (intx5(e(:,:,OP_1),ri3_79,f(:,OP_DZP),g(:,OP_DRP),h(:,OP_1))  &
          -intx5(e(:,:,OP_1),ri3_79,f(:,OP_DRP),g(:,OP_DZP),h(:,OP_1)))
      if(irunaway .gt. 2) then
         temp = temp + 1.*(gam-1.) * &
                       (intx6(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
                      + intx6(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)) &
                      - intx6(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
                      - intx6(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
      end if 
#endif
  end if

  b3psifeta = temp
end function b3psifeta


! B3bbeta
! =======
function b3bbeta(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3bbeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     temp = 0.
  else 
     temp = (gam-1.)* &
          (intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          +intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))
  end if

  b3bbeta = temp
end function b3bbeta


! B3bfeta
! =======
function b3bfeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3bfeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     temp = 0.
  else 
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
          +intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)))
#endif
     if(irunaway .gt. 2) then
        temp = temp + 1.*(gam-1.)* &
                      (intx6(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
                     - intx6(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     endif
  endif

  b3bfeta = temp
end function b3bfeta


! B3ffeta
! =======
function b3ffeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3ffeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     temp = 0.
  else 
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZP),h(:,OP_1)) &
          +intx5(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DRP),h(:,OP_1)))
     if(irunaway .gt. 0) then
        temp = temp + 1.*(gam-1.)* &
               (intx6(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
              - intx6(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     endif
#endif
  end if

  b3ffeta = temp
end function b3ffeta



! B3pepsid
! ========
function b3pepsid(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pepsid
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b3pepsid = 0.
     return
  end if

  temp = intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
       + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
       - intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_GS),h(:,OP_1)) &
       + gam* &
       (intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_DZ)) &
       +intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_DR)) &
       -intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_GS),h(:,OP_DP)))

  b3pepsid = temp
#else
  b3pepsid = 0.
#endif
end function b3pepsid


! B3pebd
! ======
function b3pebd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pebd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b3pebd = 0.
     return
  end if

  temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
        -intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
       + gam* &
       (intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DR),h(:,OP_DZ)) &
       -intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_DR)))

  b3pebd = temp
end function b3pebd


! B3pefd
! ======
function b3pefd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pefd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b3pefd = 0.
     return
  end if

  temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DRP),h(:,OP_1)) &
        -intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZP),h(:,OP_1)) &
       + gam* &
       (intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_DZ)) &
       -intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_DR)))

  b3pefd = temp
#else
  b3pefd = 0.
#endif
end function b3pefd



! B3pedkappa
! ==========
function b3pedkappa(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pedkappa
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3pedkappa = 0.
     return
  end if

  temp = &
       - intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1)) &
       - intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1 ),h(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1)) &
       - intx4(e(:,:,OP_DR),f(:,OP_1 ),g(:,OP_DR),h(:,OP_1))
  
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = h(:,OP_1)
  if(iupstream.eq.1) then    
     temp79a = temp79a + abs(i(:,OP_1))*magus
  endif
  temp = temp +                       &
       intx5(e(:,:,OP_1),ri2_79,f(:,OP_DPP),g(:,OP_1),temp79a)
  if(iupstream.eq.2) then
     temp79a = abs(i(:,OP_1))*magus
     temp = temp -                       &
          intx5(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),g(:,OP_1),temp79a)
  endif
#endif
  if(hypp.ne.0.) then
     ! Laplacian[f g]
     temp79a = f(:,OP_LP)*g(:,OP_1) + f(:,OP_1)*g(:,OP_LP) &
          + 2.*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR))

#ifdef USE3D
     temp79a = temp79a + ri2_79* &
          (f(:,OP_DPP)*g(:,OP_1) + f(:,OP_1)*g(:,OP_DPP) &
          + 2.*f(:,OP_DP)*g(:,OP_DP))
#endif

     if(ihypkappa.eq.1) then        
        temp = temp - hypp* &
             (intx3(e(:,:,OP_LP),temp79a,h(:,OP_1 )) &
             +intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DZ)) &
             +intx3(e(:,:,OP_DR),temp79a,h(:,OP_DR)))
#ifdef USE3D
        temp = temp - hypp* &
             (intx4(e(:,:,OP_DPP),ri2_79,temp79a,h(:,OP_1 )) &
             +intx4(e(:,:,OP_DP ),ri2_79,temp79a,h(:,OP_DP)))
#endif
     else
        temp = temp - hypp*intx2(e(:,:,OP_LP),temp79a)
#ifdef USE3D
        temp = temp - hypp*intx3(e(:,:,OP_DPP),ri2_79,temp79a)
#endif
     endif
  endif

  b3pedkappa = (gam-1.)*temp  
end function b3pedkappa


! B3tekappa
! ==========
function b3tekappa(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3tekappa
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3tekappa = 0.
     return
  end if

  temp = &
       - intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
       - intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1))
  
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = g(:,OP_1)
  if(iupstream.eq.1) then    
     temp79a = temp79a + abs(h(:,OP_1))*magus
  endif
  temp = temp +                       &
       intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a)
  if(iupstream.eq.2) then    
     temp79a = abs(h(:,OP_1))*magus
     temp = temp -                    &
          intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),temp79a)
  endif
#endif
  if(hypp.ne.0.) then

     if(ihypkappa.eq.1) then        
        temp = temp - hypp* &
             (intx3(e(:,:,OP_LP),f(:,OP_LP),g(:,OP_1 )) &
             +intx3(e(:,:,OP_DZ),f(:,OP_LP),g(:,OP_DZ)) &
             +intx3(e(:,:,OP_DR),f(:,OP_LP),g(:,OP_DR)))
#ifdef USE3D
        temp = temp - hypp* &
             (intx4(e(:,:,OP_LP ),ri2_79,f(:,OP_DPP),g(:,OP_1 )) &
             +intx4(e(:,:,OP_DPP),ri2_79,f(:,OP_LP ),g(:,OP_1 )) &
             +intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),g(:,OP_1 )) &
             +intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DPP),g(:,OP_DZ)) &
             +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DPP),g(:,OP_DR)) &
             +intx4(e(:,:,OP_DP),ri2_79,f(:,OP_LP ),g(:,OP_DP)) &
             +intx4(e(:,:,OP_DP),ri4_79,f(:,OP_DPP),g(:,OP_DP)))

#endif
     else
        temp = temp - hypp*intx2(e(:,:,OP_LP),f(:,OP_LP))
#ifdef USE3D
        temp = temp - hypp* &
             (intx3(e(:,:,OP_LP ),ri2_79,f(:,OP_DPP)) &
             +intx3(e(:,:,OP_DPP),ri2_79,f(:,OP_LP )) &
             +intx3(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP)))
#endif
     endif
  endif

  b3tekappa = (gam-1.)*temp
end function b3tekappa

! B3pedkappag
! ===========
function b3pppkappag(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pppkappag
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, intent(in), dimension(MAX_PTS) :: i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3pppkappag = 0.
     return
  end if

#ifdef USECOMPLEX
  temp79a = h(:,OP_DR)*conjg(g(:,OP_DR)) + h(:,OP_DZ)*conjg(g(:,OP_DZ)) &
       + h(:,OP_DP)*conjg(g(:,OP_DP))*ri2_79
#else
  temp79a = h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)
#endif

#ifdef USE3D
  temp79a = temp79a + h(:,OP_DP)*g(:,OP_DP)*ri2_79
#endif

  temp = &
       - intx4(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,i) &
       - intx4(e(:,:,OP_DR),f(:,OP_DR),temp79a,i)
  
#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a,i)
#endif

  b3pppkappag = (gam-1.)*kappag*temp
end function b3pppkappag

! B3pkappag
! =========
function b3pkappag(e,f,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pkappag
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3pkappag = 0.
     return
  end if

  temp = &
       - intx3(e(:,:,OP_DZ),f(:,OP_DZ),i) &
       - intx3(e(:,:,OP_DR),f(:,OP_DR),i)
  
#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),i)
#endif

  b3pkappag = -gradp_crit**2*(gam-1.)*kappag*temp
end function b3pkappag



!============================================================================
! N1 TERMS
!============================================================================

! N1n
! ===
function n1n(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element,dofs_per_element) :: n1n
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e,f

  n1n = intxx2(e(:,:,OP_1),f(:,:,OP_1))
end function n1n


! N1ndenm
! =======
function n1ndenm(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1ndenm
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = -(intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
       +   intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1)))

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = g(:,OP_1)
  if(iupstream .eq. 1) then   
     temp79a = temp79a+abs(h(:,OP_1))*magus
  endif
  temp = temp + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a) &
       + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_DP))
  if(iupstream .eq. 2) then   
     temp79a = abs(h(:,OP_1))*magus
     temp = temp - intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),temp79a)
  endif
#endif

  if(hypp.ne.0.) then
     if(ihypkappa.eq.1) then
        temp = temp - hypp*intx3(e(:,:,OP_LP),f(:,OP_LP),g(:,OP_1))
#ifdef USE3D
        temp = temp - hypp* &
             (intx4(e(:,:,OP_LP ),ri2_79,f(:,OP_DPP),g(:,OP_1)) &
             +intx4(e(:,:,OP_DPP),ri2_79,f(:,OP_LP ),g(:,OP_1)) &
             +intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),g(:,OP_1)))
#endif
     else
        temp = temp - hypp*intx2(e(:,:,OP_LP),f(:,OP_LP))
#ifdef USE3D
        temp = temp - hypp* &
             (intx3(e(:,:,OP_LP ),ri2_79,f(:,OP_DPP)) &
             +intx3(e(:,:,OP_DPP),ri2_79,f(:,OP_LP )) &
             +intx3(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP)))
#endif
     endif
  endif

  n1ndenm = temp
end function n1ndenm


! N1nu
! ====
function n1nu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1nu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ))

  n1nu = temp
end function n1nu


! N1nv
! ====
function n1nv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1nv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = -intx3(e(:,:,OP_1),f(:,OP_1 ),g(:,OP_DP)) &
       -  intx3(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1 ))

  n1nv = temp
#else
  n1nv = 0.
#endif
end function n1nv

! N1nchi
! ======
function n1nchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1nchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZ)) &
       + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DR))

  n1nchi = temp
end function n1nchi


! N1s
! ===
function n1s(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1s
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  n1s = intx2(e(:,:,OP_1),f(:,OP_1))
end function n1s

function t3tndenm(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tndenm
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp =  &
       (intx4(e(:,:,OP_DZ),f(:,OP_1),g(:,OP_DZ),h(:,OP_1)) &
       +intx4(e(:,:,OP_DR),f(:,OP_1),g(:,OP_DR),h(:,OP_1)) &
       +intx4(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
       +intx4(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))

#if defined(USE3D) || defined(USECOMPLEX)
  temp = temp - &
       (   intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DPP),h(:,OP_1 )) &
       +   intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DP ),h(:,OP_DP)))
#endif

  if(hypp.ne.0.) then
     if(ihypkappa.eq.1) then
        temp = temp + hypp* &
             ( intx4(e(:,:,OP_LP),f(:,OP_1),g(:,OP_LP),h(:,OP_1))   &
             + intx4(e(:,:,OP_1),f(:,OP_LP),g(:,OP_LP),h(:,OP_1))   &
             + 2.*intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_LP),h(:,OP_1))&
             + 2.*intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_LP),h(:,OP_1)))
#ifdef USE3D
        temp = temp + hypp* &
             ( intx5(e(:,:,OP_DPP),ri2_79,f(:,OP_1),g(:,OP_LP ),h(:,OP_1))    &
             + intx5(e(:,:,OP_LP ),ri2_79,f(:,OP_1),g(:,OP_DPP),h(:,OP_1))    &
             + intx5(e(:,:,OP_DPP),ri4_79,f(:,OP_1),g(:,OP_DPP),h(:,OP_1))    &
             + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DPP),g(:,OP_LP ),h(:,OP_1))    &
             + intx5(e(:,:,OP_1),ri2_79,f(:,OP_LP ),g(:,OP_DPP),h(:,OP_1))    &
             + intx5(e(:,:,OP_1),ri4_79,f(:,OP_DPP),g(:,OP_DPP),h(:,OP_1))    &
             + 2.*intx5(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_DPP),h(:,OP_1)) &
             + 2.*intx5(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_DPP),h(:,OP_1)) &
             + 2.*intx5(e(:,:,OP_DP),ri2_79,f(:,OP_DP),g(:,OP_LP ),h(:,OP_1)) &
             + 2.*intx5(e(:,:,OP_DP),ri4_79,f(:,OP_DP),g(:,OP_DPP),h(:,OP_1)))
#endif

     else
        temp = temp + hypp* &
             ( intx3(e(:,:,OP_LP),f(:,OP_1),g(:,OP_LP))   &
             + intx3(e(:,:,OP_1),f(:,OP_LP),g(:,OP_LP))   &
             + 2.*intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_LP)) &
             + 2.*intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_LP)))
#ifdef USE3D
        temp = temp + hypp* &
             ( intx4(e(:,:,OP_DPP),ri2_79,f(:,OP_1),g(:,OP_LP ))    &
             + intx4(e(:,:,OP_LP ),ri2_79,f(:,OP_1),g(:,OP_DPP))    &
             + intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_1),g(:,OP_DPP))    &
             + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),g(:,OP_LP))     &
             + intx4(e(:,:,OP_1),ri2_79,f(:,OP_LP ),g(:,OP_DPP))    &
             + intx4(e(:,:,OP_1),ri4_79,f(:,OP_DPP),g(:,OP_DPP))    &
             + 2.*intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_DPP)) &
             + 2.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_DPP)) &
             + 2.*intx4(e(:,:,OP_DP),ri2_79,f(:,OP_DP),g(:,OP_LP )) &
             + 2.*intx4(e(:,:,OP_DP),ri4_79,f(:,OP_DP),g(:,OP_DPP)))
#endif
     endif
  endif

  t3tndenm = temp
end function t3tndenm


function t3ts(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3ts
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = -intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_1))

  t3ts = temp
end function t3ts

!============================================================================
! NRE1 TERMS
!============================================================================

! NRE1nrediff
! =======
function nre1nrediff(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nrediff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D)
  ! B . Grad(p)
  ! [ p, psi] / R + 1/R^2 F p' - <f', p>
  temp79a = ri_79*(f(:,OP_DZ)*b2i79(:,OP_1)*pstx79(:,OP_DR) - &
            f(:,OP_DR)*b2i79(:,OP_1)*pstx79(:,OP_DZ))
  temp79a = temp79a+0.5*ri_79*(f(:,OP_1)*b2i79(:,OP_DZ)*pstx79(:,OP_DR) - &
            f(:,OP_1)*b2i79(:,OP_DR)*pstx79(:,OP_DZ))
  temp79a = temp79a &
          + ri2_79*bztx79(:,OP_1)*b2i79(:,OP_1)*f(:,OP_DP) &
          - b2i79(:,OP_1)*(bfptx79(:,OP_DZ)*f(:,OP_DZ) & 
          + bfptx79(:,OP_DR)*f(:,OP_DR))

  temp = intx4(e(:,:,OP_DZ),ri_79,-temp79a,pstx79(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),ri_79,-temp79a,pstx79(:,OP_DZ))
  temp = temp &
       + intx4(e(:,:,OP_DP),ri2_79,-temp79a,bztx79(:,OP_1)) &
       - intx3(e(:,:,OP_DZ),-temp79a,bfptx79(:,OP_DZ)) &
       - intx3(e(:,:,OP_DR),-temp79a,bfptx79(:,OP_DR))
#else
  temp79a = ri_79*(f(:,OP_DZ)*b2i79(:,OP_1)*pstx79(:,OP_DR) - &
            f(:,OP_DR)*b2i79(:,OP_1)*pstx79(:,OP_DZ))
  temp79a = temp79a+0.5*ri_79*(f(:,OP_1)*b2i79(:,OP_DZ)*pstx79(:,OP_DR) - &
            f(:,OP_1)*b2i79(:,OP_DR)*pstx79(:,OP_DZ))

  temp = intx4(e(:,:,OP_DZ),ri_79,-temp79a,pstx79(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),ri_79,-temp79a,pstx79(:,OP_DZ))
#endif

  nre1nrediff = temp
end function nre1nrediff

! NRE1nrepsi
! ====
function nre1nrepsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nrepsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = intx5(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_1),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_1),h(:,OP_DZ))
  temp = temp * 1.000 * cre 

  nre1nrepsi = temp
end function nre1nrepsi

! NRE1nreb
! ====
function nre1nreb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nreb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = 0
#if defined(USE3D)
  temp = -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_1),h(:,OP_1)) - &
          intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DP),h(:,OP_1)) - &
          intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_1),h(:,OP_DP))
#elif defined(USECOMPLEX)
  temp = -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_1),h(:,OP_1))
#endif
  temp = temp * 1.000 * cre

  nre1nreb = temp
end function nre1nreb

! NRE1nreu
! ====
function nre1nreu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nreu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ))

  nre1nreu = temp
end function nre1nreu

!============================================================================
! P1 TERMS
!============================================================================

! P1pu
! ====
function p1pu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1pu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ))

  if(itor.eq.1) then
     temp = temp + &
          2.*(gam-1.)*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DZ))
  endif

  p1pu = temp
end function p1pu


! P1pv
! ====
function p1pv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1pv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = - intx3(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1)) &
       - gam*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DP))
#else
  temp = 0.
#endif

  p1pv = temp
end function p1pv


! P1pchi
! ======
function p1pchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1pchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = intx4(e(:,:,OP_DR),ri2_79,g(:,OP_DR),f(:,OP_1)) &
       + intx4(e(:,:,OP_DZ),ri2_79,g(:,OP_DZ),f(:,OP_1)) &
       - (gam-1.)*intx4(e(:,:,OP_1),ri2_79,g(:,OP_GS),f(:,OP_1))

  p1pchi = temp
end function p1pchi


! P1psipsikappar
! ==============
function p1psipsikappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psipsikappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psipsikappar = 0.
     return
  end if

  temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)*i(:,OP_1)
  temp79b = ri2_79*k(:,OP_1)*j(:,OP_1)*h(:,OP_1)

  temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DZ),i(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DZ),i(:,OP_DR)) &
       - intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DR),i(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DR),i(:,OP_DZ))

  p1psipsikappar = (gam - 1.) * temp
end function p1psipsikappar

! P1psipsipnkappar
! ================
function p1psipsipnkappar(e,f,g,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psipsipnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. fac1.eq.0) then
     p1psipsipnkappar = 0.
     return
  end if

  temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! [T,psi]*n = [n p/n^2,psi]*n
  temp79b = ni79(:,OP_1)* &
       (i(:,OP_1)*(h(:,OP_DZ)*g(:,OP_DR)-h(:,OP_DR)*g(:,OP_DZ)) &
       +h(:,OP_1)*(i(:,OP_DZ)*g(:,OP_DR)-i(:,OP_DR)*g(:,OP_DZ))) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*g(:,OP_DR) - ni79(:,OP_DR)*g(:,OP_DZ))

  temp = intx4(e(:,:,OP_DZ),temp79a,temp79b,f(:,OP_DR)) &
       - intx4(e(:,:,OP_DR),temp79a,temp79b,f(:,OP_DZ))

  p1psipsipnkappar = (gam - 1.) * temp
end function p1psipsipnkappar


! P1psibkappar
! ============
function p1psibkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psibkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psibkappar = 0.
     return
  end if

  temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)*j(:,OP_1)  

  temp79b = f(:,OP_DR)*(h(:,OP_DZ)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DZ)) &
       -    f(:,OP_DZ)*(h(:,OP_DR)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DR))

  temp79c = f(:,OP_DRP)*(h(:,OP_DZ )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DZ )) &
       -    f(:,OP_DZP)*(h(:,OP_DR )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DR )) &
       +    f(:,OP_DR )*(h(:,OP_DZP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DZ )) &
       -    f(:,OP_DZ )*(h(:,OP_DRP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DR )) &
       +    f(:,OP_DR )*(h(:,OP_DZ )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DZP)) &
       -    f(:,OP_DZ )*(h(:,OP_DR )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DRP))

  temp79d = temp79c*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ) &
       +    temp79b*g(:,OP_DP)*j(:,OP_1 )*k(:,OP_1 ) &
       +    temp79b*g(:,OP_1 )*j(:,OP_DP)*k(:,OP_1 ) &
       +    temp79b*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_DP)

  temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,h(:,OP_DP),i(:,OP_1 )) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,h(:,OP_DP),i(:,OP_1 )) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
       + intx3(e(:,:,OP_1),ri3_79,temp79d)

  p1psibkappar = (gam - 1.) * temp
#else
  p1psibkappar = 0.
#endif
end function p1psibkappar

! P1psibpnkappar
! ==============
function p1psibpnkappar(e,f,g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psibpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psibpnkappar = 0.
     return
  end if

  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1)

  ! n*dT/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  ! d(temp79a)/dphi
  temp79b = kar79(:,OP_DP)*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_1 ) &
       +    kar79(:,OP_1 )*b2i79(:,OP_DP)*ni79(:,OP_1 )*g(:,OP_1 ) &
       +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_DP)*g(:,OP_1 ) &
       +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_DP)

  ! [T,psi]*n = [n p/n^2,psi]*n
  temp79c = ni79(:,OP_1)* &
       (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
       +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))

  ! d(temp79c)/dphi
  temp79d = ni79(:,OP_DP)* &
       (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
       +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZP)*f(:,OP_DR) - ni79(:,OP_DRP)*f(:,OP_DZ)) &
       +    ni79(:,OP_1)* &
       (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DRP)-h(:,OP_DR)*f(:,OP_DZP)) &
       +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DRP)-i(:,OP_DR)*f(:,OP_DZP))) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*f(:,OP_DRP) - ni79(:,OP_DR)*f(:,OP_DZP)) &
       +    ni79(:,OP_1)* &
       (i(:,OP_1 )*(h(:,OP_DZP)*f(:,OP_DR)-h(:,OP_DRP)*f(:,OP_DZ)) &
       +h(:,OP_DP)*(i(:,OP_DZ )*f(:,OP_DR)-i(:,OP_DR )*f(:,OP_DZ))) &
       + 2.*h(:,OP_DP)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ)) &
       +    ni79(:,OP_1)* &
       (i(:,OP_DP)*(h(:,OP_DZ )*f(:,OP_DR)-h(:,OP_DR )*f(:,OP_DZ)) &
       +h(:,OP_1 )*(i(:,OP_DZP)*f(:,OP_DR)-i(:,OP_DRP)*f(:,OP_DZ))) &
       + 2.*h(:,OP_1)*i(:,OP_DP)* &
       (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))
  
  temp = fac2*intx4(e(:,:,OP_1),ri3_79,temp79a,temp79d) &
       + fac2*intx4(e(:,:,OP_1),ri3_79,temp79b,temp79c) &
       + fac1*intx5(e(:,:,OP_DR),ri3_79,temp79a,f(:,OP_DZ),temp79e) &
       - fac1*intx5(e(:,:,OP_DZ),ri3_79,temp79a,f(:,OP_DR),temp79e)

  p1psibpnkappar = (gam - 1.) * temp
#else
  p1psibpnkappar = 0.
#endif
end function p1psibpnkappar


! P1bbkappar
! ==========
function p1bbkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bbkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1bbkappar = 0.
     return
  end if

  temp79a = h(:,OP_DP)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DP)
  temp79b = h(:,OP_DPP)*i(:,OP_1  ) &
       + 2.*h(:,OP_DP )*i(:,OP_DP ) &
       +    h(:,OP_1  )*i(:,OP_DPP)

  temp79c = f(:,OP_DP)*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 ) &
       +    f(:,OP_1 )*g(:,OP_DP)*temp79a*j(:,OP_1 )*k(:,OP_1 ) &
       +    f(:,OP_1 )*g(:,OP_1 )*temp79b*j(:,OP_1 )*k(:,OP_1 ) &
       +    f(:,OP_1 )*g(:,OP_1 )*temp79a*j(:,OP_DP)*k(:,OP_1 ) &
       +    f(:,OP_1 )*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_DP)

  temp = intx3(e(:,:,OP_1),ri4_79,temp79c)

  p1bbkappar = (gam - 1.) * temp
#else
  p1bbkappar = 0.
#endif
end function p1bbkappar


! P1bbpnkappar
! ===========
function p1bbpnkappar(e,f,g,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bbpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  
  if(gam.le.1. .or. fac1.eq.0) then
     p1bbpnkappar = 0.
     return
  end if

  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1)

  ! d(temp79a)/dphi
  temp79b = kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
       +    kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
       +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
       +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
       +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

  ! n*dT/dphi
  temp79c = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)
  ! d(temp79c)/dphi
  temp79d = ni79(:,OP_DP)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DPP) &
       +    ni79(:,OP_1)*(h(:,OP_DP)*i(:,OP_DP) + h(:,OP_DPP)*i(:,OP_1)) &
       + 2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DP) &
       +    ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DPP) + h(:,OP_DP)*i(:,OP_DP)) &
       + 2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DP)

  temp = intx4(e(:,:,OP_1),ri4_79,temp79a,temp79d) &
       + intx4(e(:,:,OP_1),ri4_79,temp79b,temp79c)

  p1bbpnkappar = (gam - 1.) * temp
#else
  p1bbpnkappar = 0.
#endif
end function p1bbpnkappar



! P1psifkappar
! ============
function p1psifkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psifkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psifkappar = 0.
     return
  end if

  temp79a = k(:,OP_1)*ri_79*j(:,OP_1)*i(:,OP_1)
  temp79b = k(:,OP_1)*ri_79*j(:,OP_1)*h(:,OP_1)

  temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DR),i(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DR),i(:,OP_DR)) &
       + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
       - intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DZ ),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DZ ),h(:,OP_DR)) &
       + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79b,f(:,OP_DR ),i(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79b,f(:,OP_DR ),i(:,OP_DZ)) &
       - intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79b,f(:,OP_DZ ),i(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),g(:,OP_DR),temp79b,f(:,OP_DZ ),i(:,OP_DR))

  p1psifkappar = (gam - 1.) * temp
#else
  p1psifkappar = 0.
#endif
end function p1psifkappar


! P1psifpnkappar
! ==============
function p1psifpnkappar(e,f,g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psifpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psifpnkappar = 0.
     return
  end if

  temp79a = ri_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! n*<T,f'>
  temp79b = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))
  ! n*[T,psi]
  temp79c = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DZ)*f(:,OP_DR) - i(:,OP_DR)*f(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DZ)*f(:,OP_DR) - h(:,OP_DR)*f(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))
  
  temp = fac1*intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DR),temp79b) &
       - fac1*intx4(e(:,:,OP_DR),temp79a,f(:,OP_DZ),temp79b) &
       + fac2*intx4(e(:,:,OP_DR),temp79a,g(:,OP_DR),temp79c) &
       + fac2*intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DZ),temp79c)

  p1psifpnkappar = (gam - 1.) * temp
#else
  p1psifpnkappar = 0.
#endif
end function p1psifpnkappar

! P1qpsikappar
! ==============
function p1qpsikappar(e,f,g,i,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1qpsikappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,k
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1qpsikappar = 0.
     return
  end if

  temp79a = ri_79*k(:,OP_1)*i(:,OP_1)
  
  temp = intx4(e(:,:,OP_DR),temp79a,g(:,OP_DZ),f(:,OP_1)) &
       - intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DR),f(:,OP_1)) 

  p1qpsikappar = (gam - 1.) * temp
  return
end function p1qpsikappar

! P1bfkappar
! ==========
function p1bfkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bfkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1bfkappar = 0.
     return
  end if

  temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)*j(:,OP_1)

  temp79b = g(:,OP_DZ)*(h(:,OP_DZ)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DZ)) &
       +    g(:,OP_DR)*(h(:,OP_DR)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DR))

  temp79c = g(:,OP_DZP)*(h(:,OP_DZ )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DZ )) &
       +    g(:,OP_DRP)*(h(:,OP_DR )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DR )) &
       +    g(:,OP_DZ )*(h(:,OP_DZP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DZ )) &
       +    g(:,OP_DR )*(h(:,OP_DRP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DR )) &
       +    g(:,OP_DZ )*(h(:,OP_DZ )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DZP)) &
       +    g(:,OP_DR )*(h(:,OP_DR )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DRP))

  temp79d = temp79c*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ) &
       +    temp79b*f(:,OP_DP)*j(:,OP_1 )*k(:,OP_1 ) &
       +    temp79b*f(:,OP_1 )*j(:,OP_DP)*k(:,OP_1 ) &
       +    temp79b*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_DP)

  temp = intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,h(:,OP_DP),i(:,OP_1 )) &
       + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,h(:,OP_DP),i(:,OP_1 )) &
       + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
       + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
       - intx3(e(:,:,OP_1),ri2_79,temp79d)

  p1bfkappar = (gam - 1.) * temp
#else
  p1bfkappar = 0.
#endif
end function p1bfkappar


! P1qbkappar
! ==========
function p1qbkappar(e,f,g,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1qbkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1qbkappar = 0.
     return
  end if

  temp79a =  ri2_79*i(:,OP_1)*j(:,OP_1)*g(:,OP_1)

  temp = -intx3(e(:,:,OP_DP),temp79a,f(:,OP_1 )) 

  p1qbkappar = (gam - 1.) * temp
#else
  p1qbkappar = 0.
#endif

end function p1qbkappar

! ==========
function p1bfpnkappar(e,f,g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bfpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1bfpnkappar = 0.
     return
  end if

  temp79a = ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)

  ! n*dT/dphi = n*d(n p/n^2)/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  ! d(temp79a)/dphi
  temp79b = ri2_79 * &
       (kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1) &
       +kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1) &
       +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1) &
       +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP))

  ! n*<T, f'> = n*<n p/n^2, f'>
  temp79c = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))

  ! d(temp79c)/dphi
  temp79d = ni79(:,OP_DP)*h(:,OP_1)* &
       (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
       + ni79(:,OP_DP)*i(:,OP_1)* &
       (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DRP)*g(:,OP_DR) + ni79(:,OP_DZP)*g(:,OP_DZ)) &
       +    ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DR)*g(:,OP_DRP) + i(:,OP_DZ)*g(:,OP_DZP)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DR)*g(:,OP_DRP) + ni79(:,OP_DZ)*g(:,OP_DZP)) &
       +    ni79(:,OP_1)*h(:,OP_DP)* &
       (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DRP)*g(:,OP_DR) + h(:,OP_DZP)*g(:,OP_DZ)) &
       + 2.*h(:,OP_DP)*i(:,OP_1)* &
       (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ)) &
       +    ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DRP)*g(:,OP_DR) + i(:,OP_DZP)*g(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_DP)* &
       (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_DP)* &
       (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))

  temp = fac2*intx4(e(:,:,OP_DR),temp79a,g(:,OP_DR),temp79e) &
       + fac2*intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DZ),temp79e) &
       - fac1*intx3(e(:,:,OP_1),temp79a,temp79d) &
       - fac1*intx3(e(:,:,OP_1),temp79b,temp79c)

  p1bfpnkappar = (gam - 1.) * temp
#else
  p1bfpnkappar = 0.
#endif
end function p1bfpnkappar


! P1ffkappar
! ==========
function p1ffkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1ffkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1ffkappar = 0.
     return
  end if

  temp79a = -k(:,OP_1)*j(:,OP_1)*i(:,OP_1)
  temp79b = -k(:,OP_1)*j(:,OP_1)*h(:,OP_1)

  temp = intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79b,g(:,OP_DR),i(:,OP_DR)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79b,g(:,OP_DR),i(:,OP_DR)) 

  p1ffkappar = (gam - 1.) * temp
#else
  p1ffkappar = 0.
#endif
end function p1ffkappar


! P1ffpnkappar
! ============
function p1ffpnkappar(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1ffpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1ffpnkappar = 0.
     return
  end if

  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! <T, f'>*n = <n p/n^2, f'>*n
  temp79c = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))

  temp = -intx4(e(:,:,OP_DR),temp79a,f(:,OP_DR),temp79c) &
       -  intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DZ),temp79c)

  p1ffpnkappar = (gam - 1.) * temp
#else
  p1ffpnkappar = 0.
#endif
end function p1ffpnkappar


! P1qfkappar
! ============
vectype function p1qfkappar(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  if(gam.le.1.) then
     p1qfkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = i(:,OP_1)*h(:,OP_1)

  temp = int4(temp79a,e(:,OP_DR),g(:,OP_DR),f(:,OP_1)) &
       + int4(temp79a,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_1))
#else
  temp = 0.
#endif

  p1qfkappar = (gam - 1.) * temp
end function p1qfkappar

! P1kappax
! ========
function p1kappax(e,f,g,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1kappax
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,j
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1kappax = 0.
     return
  end if

  temp79a = ri_79*i(:,OP_1)*g(:,OP_1)*j(:,OP_1)

  temp = intx3(e(:,:,OP_DZ),f(:,OP_DR),temp79a) &
       - intx3(e(:,:,OP_DR),f(:,OP_DZ),temp79a)

  p1kappax = (gam - 1.) * temp
end function p1kappax



! P1uus
! =====
function p1uus(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1uus
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1.) then
     p1uus = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = 0.5*(gam-1.)* &
       (intx5(e(:,:,OP_1),r2_79,f(:,OP_DZ),g(:,OP_DZ),temp79a) &
       +intx5(e(:,:,OP_1),r2_79,f(:,OP_DR),g(:,OP_DR),temp79a))

  p1uus = temp
end function p1uus


! P1vvs
! =====
function p1vvs(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1vvs
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1.) then
     p1vvs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = 0.5*(gam-1.)* &
       intx5(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_1),temp79a)

  p1vvs = temp
end function p1vvs


! P1chichis
! =========
function p1chichis(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1chichis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1.) then
     p1chichis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = 0.5*(gam-1.)* &
       (intx5(e(:,:,OP_1),ri4_79,f(:,OP_DZ),g(:,OP_DZ),temp79a) &
       +intx5(e(:,:,OP_1),ri4_79,f(:,OP_DR),g(:,OP_DR),temp79a))
     
  p1chichis = temp
end function p1chichis


! P1uchis
! =======
function p1uchis(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1uchis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1.) then
     p1uchis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = -(gam-1.)* & 
       (intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),temp79a) &
       -intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),temp79a))

  p1uchis = temp
end function p1uchis

!
! Extra diffusion to model upstream differencing
vectype function p1uspu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp 

  temp79a = abs(g(:,OP_DZ))
  temp79b = abs(g(:,OP_DR))

  temp = - int4(r_79,e(:,OP_DR),f(:,OP_DR),temp79a)  &
       - int4(r_79,e(:,OP_DZ),f(:,OP_DZ),temp79b)

  p1uspu = temp
end function p1uspu

!
! Extra diffusion to model upstream differencing
vectype function p1uspchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp 

  temp79a = abs(g(:,OP_DR))
  temp79b = abs(g(:,OP_DZ))

  temp = - int4(ri2_79,e(:,OP_DR),f(:,OP_DR),temp79a)  &
       -   int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),temp79b)

  p1uspchi = temp
end function p1uspchi

! Extra diffusion to model upstream differencing
vectype function p1uspv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype :: temp 

  temp79a = abs(g(:,OP_1))

  temp =  int3(e(:,OP_1),f(:,OP_DPP),temp79a)  

  p1uspv = temp
#else

  p1uspv = 0.
#endif
end function p1uspv


!======================================================================
! Parallel Viscous Terms
!======================================================================

! PVS1
! ====
subroutine PVS1(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

     o = 0.
     if(itor.eq.1) then
        o = o - (1./3.)*i(:,OP_DZ)
     endif

end subroutine PVS1


! PVS1psipsi
! ==========
subroutine PVS1psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o


     o = r_79* &
          (j(:,OP_DR)*k(:,OP_DZ)*(i(:,OP_DZZ) - i(:,OP_DRR)) &
          -(j(:,OP_DZ)*k(:,OP_DZ) - j(:,OP_DR)*k(:,OP_DR))*i(:,OP_DRZ))
     if(itor.eq.1) then
        o = o + j(:,OP_DR)* &
             (i(:,OP_DZ)*k(:,OP_DR) - i(:,OP_DR)*k(:,OP_DZ))
     endif

     o = o * ri2_79*b2i79(:,OP_1)

end subroutine PVS1psipsi

! PVS1psib
! ========
subroutine PVS1psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

     o = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + k(:,OP_1)* &
          (j(:,OP_DZ)*i(:,OP_DZP) + j(:,OP_DR)*i(:,OP_DRP))
     o = o * ri2_79*b2i79(:,OP_1)
#endif     

end subroutine PVS1psib


! PVS1bb
! ======
subroutine PVS1bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k 
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = 0.

end subroutine PVS1bb


! PVS2
! ====
subroutine PVS2(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

#if defined(USE3D) || defined(USECOMPLEX)
  o = -i(:,OP_DP)/3.
#else
  o = 0.
#endif

end subroutine PVS2


! PVS2psipsi
! ==========
subroutine PVS2psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = 0.
end subroutine PVS2psipsi


! PVS2psib
! ========
subroutine PVS2psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = r_79*k(:,OP_1)* &
       (i(:,OP_DZ)*j(:,OP_DR) - i(:,OP_DR)*j(:,OP_DZ))

  o = o * ri2_79*b2i79(:,OP_1)

end subroutine PVS2psib

! PVS2bb
! ======
subroutine PVS2bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k
  vectype, intent(out), dimension(MAX_PTS) :: o

#if defined(USE3D) || defined(USECOMPLEX)
  o = j(:,OP_1)*k(:,OP_1)*i(:,OP_DP)
  o = o * ri2_79*b2i79(:,OP_1)
#else
  o = 0.
#endif

end subroutine PVS2bb


! PVS3
! ====
subroutine PVS3(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = -(1./3.)*ri2_79*i(:,OP_GS)

  if(itor.eq.1) then
     o = o + (1./3.)*ri3_79*i(:,OP_DR)
  end if
 
end subroutine PVS3


! PVS3psipsi
! ==========
subroutine PVS3psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = -ri2_79* &
       (j(:,OP_DZ)*k(:,OP_DZ)*i(:,OP_DZZ) &
       +j(:,OP_DR)*k(:,OP_DR)*i(:,OP_DRR) &
       +2.*j(:,OP_DZ)*k(:,OP_DR)*i(:,OP_DRZ) &
       -i(:,OP_GS)*(j(:,OP_DZ)*k(:,OP_DZ) + j(:,OP_DR)*k(:,OP_DR))) 
  
  if(itor.eq.1) then
     o = o + 2.*ri3_79*j(:,OP_DZ) * &
          (i(:,OP_DZ)*k(:,OP_DR) - i(:,OP_DR)*k(:,OP_DZ))
     
  endif

  o = o * ri2_79*b2i79(:,OP_1)
 
end subroutine PVS3psipsi


! PVS3psib
! ========
subroutine PVS3psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

#if defined(USE3D) || defined(USECOMPLEX)
  o = ri3_79*k(:,OP_1) * &
       (i(:,OP_DZP)*j(:,OP_DR) - i(:,OP_DRP)*j(:,OP_DZ))
  o = o * ri2_79*b2i79(:,OP_1)
#else
  o = 0.
#endif
 
end subroutine PVS3psib


! PVS3bb
! ======
subroutine PVS3bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = 0.
 
end subroutine PVS3bb



! PVV1
! ====
subroutine PVV1(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

  o =   (e(:,OP_DZZ) - e(:,OP_DRR))*pst79(:,OP_DR)*pst79(:,OP_DZ) &
       + e(:,OP_DRZ)*(pst79(:,OP_DR)**2 - pst79(:,OP_DZ)**2)
        
  if(itor.eq.1) then
     o = o - ri_79* &
          (pst79(:,OP_DZ)* &
          (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR)) &
          +e(:,OP_DZ)*bzt79(:,OP_1)**2)
  endif
        
  o = 3.*ri_79*b2i79(:,OP_1)*o
        
  if(itor.eq.1) then
     o = o + 2.*e(:,OP_DZ)
  endif
        
#if defined(USECOMPLEX)
  o = o - rfac*3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
       (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR))
#endif

  o = -o
end subroutine  PVV1

! PVV2
! ====
subroutine PVV2(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = ri_79*(e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ))

  o = -3.*b2i79(:,OP_1)*bzt79(:,OP_1)*o

#if defined(USECOMPLEX)
  ! This term is a total derivative in phi
  ! and therefore integrates out of the 3D case

  o = o - rfac*e(:,OP_1) * &
       (1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)
#endif
end subroutine  PVV2


! PVV3
! ====
subroutine PVV3(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

  o =      e(:,OP_DZZ)*pst79(:,OP_DR)**2 &
       +   e(:,OP_DRR)*pst79(:,OP_DZ)**2 &
       -2.*e(:,OP_DRZ)*pst79(:,OP_DZ)*pst79(:,OP_DR)
        
  if(itor.eq.1) then
     o = o + 2.*ri_79*pst79(:,OP_DZ)* &
          (e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ)) &
          + ri_79*e(:,OP_DR)*bzt79(:,OP_1)**2
  endif
        
  o = 3.*ri4_79*b2i79(:,OP_1)*o - ri2_79*e(:,OP_GS)
        
#if defined(USECOMPLEX)
  o = o - rfac*3.*ri5_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
       (e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ))
#endif
     
end subroutine  PVV3


! P1vip
! =====
vectype function P1vip(e)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e

  if(amupar.eq.0) then
     P1vip = 0.
     return
  endif

  call PVS1      (pht79,temp79b)
  call PVS1psipsi(pht79,pst79,pst79,temp79c)
  call PVS1psib  (pht79,pst79,bzt79,temp79d)
  call PVS1bb    (pht79,bzt79,bzt79,temp79e)
  temp79a = temp79b + temp79c + temp79d + temp79e

  if(numvar.ge.2) then
     call PVS2      (vzt79,temp79b)
     call PVS2psipsi(vzt79,pst79,pst79,temp79c)
     call PVS2psib  (vzt79,pst79,bzt79,temp79d)
     call PVS2bb    (vzt79,bzt79,bzt79,temp79e)
     temp79a = temp79a + temp79b + temp79c + temp79d + temp79e
  endif

  if(numvar.ge.3) then
     call PVS3      (cht79,temp79b)
     call PVS3psipsi(cht79,pst79,pst79,temp79c)
     call PVS3psib  (cht79,pst79,bzt79,temp79d)
     call PVS3bb    (cht79,bzt79,bzt79,temp79e)
     temp79a = temp79a + temp79b + temp79c + temp79d + temp79e
  endif

  P1vip = 3.*int4(e(:,OP_1),vip79(:,OP_1),temp79a,temp79a)
end function P1vip


!======================================================================
! Gyroviscous terms
!======================================================================

! g1u
! ===
vectype function g1u(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g1u = gyro_vor_u(e,f)

end function g1u

! g1v
! ===
vectype function g1v(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g1v = gyro_vor_v(e,f)

end function g1v

! g1chi
! =====
vectype function g1chi(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g1chi = gyro_vor_x(e,f)

end function g1chi


! g2u
! ===
vectype function g2u(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g2u = gyro_tor_u(e,f)

end function g2u


! g2v
! ===
vectype function g2v(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g2v = gyro_tor_v(e,f)

end function g2v


! g2chi
! =====
vectype function g2chi(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g2chi = gyro_tor_x(e,f)

end function g2chi


! g3u
! ===
vectype function g3u(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g3u = gyro_com_u(e,f)

end function g3u


! g3v
! ===
vectype function g3v(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g3v = gyro_com_v(e,f)

end function g3v


! g3chi
! =====
vectype function g3chi(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  g3chi = gyro_com_x(e,f)

end function g3chi





! ==============================================================
! Ohmic heating terms
! ==============================================================

! qpsipsieta
! ==========
function qpsipsieta(e)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qpsipsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, dimension(dofs_per_element) :: temp

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
     temp = hypf*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = hypf*intx2(e(:,:,OP_1),temp79a)
  endif


  temp79a = jt79(:,OP_DZ)*ni79(:,OP_DZ) + jt79(:,OP_DR)*ni79(:,OP_DR)
  if(itor.eq.1) then
     temp79a = temp79a - ri_79*jt79(:,OP_1)*ni79(:,OP_DR)
  endif
  temp79a = temp79a * ri2_79*nt79(:,OP_1)*jt79(:,OP_1)
     
  if(ihypeta.eq.1) then
     temp = temp + hypf*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = temp + hypf*intx2(e(:,:,OP_1),temp79a)
  endif

  qpsipsieta = temp
end function qpsipsieta

! qbbeta
! ======
function qbbeta(e)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qbbeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, dimension(dofs_per_element) :: temp

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
     temp = hypi*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = hypi*intx2(e(:,:,OP_1),temp79a)
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
     temp = temp + hypi*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = temp + hypi*intx2(e(:,:,OP_1),temp79a)
  endif

  qbbeta = temp
end function qbbeta


! ==============================================================
! Viscous heating terms
! ==============================================================

! quumu
! =====
function quumu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: quumu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
        
  ! Not yet implemented

  quumu = 0.
end function quumu


! quchimu
! =======
function quchimu(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: quchimu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i

  quchimu = 0.
end function quchimu


! qvvmu
! =====
function qvvmu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qvvmu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = - &
       (intx5(e(:,:,OP_1),r2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
       +intx5(e(:,:,OP_1),r2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))
        
  if(hypv.ne.0.) then
     temp = temp - &
          hypv*intx5(e(:,:,OP_1),r2_79,f(:,OP_GS),g(:,OP_GS),h(:,OP_1))
  endif

  qvvmu = temp
end function qvvmu


! qchichimu
! =========
function qchichimu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qchichimu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  ! Not yet implemented

  qchichimu = 0.
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
real function energy_mp(mask)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), optional, dimension(MAX_PTS) :: mask

  vectype :: temp

  if(present(mask)) then
     temp79a = mask
  else
     temp79a = 1.
  end if


  if(linear.eq.1) then
     temp = .5* &
          (int4(ri2_79,ps179(:,OP_DZ),CONJUGATE(ps179(:,OP_DZ)),temp79a) &
          +int4(ri2_79,ps179(:,OP_DR),CONJUGATE(ps179(:,OP_DR)),temp79a))
#ifdef USECOMPLEX
     if(numvar.gt.1) then
        temp = temp + .5* &
             (int3(bfp179(:,OP_DZ),CONJUGATE(bfp179(:,OP_DZ)),temp79a) &
             +int3(bfp179(:,OP_DR),CONJUGATE(bfp179(:,OP_DR)),temp79a) &
             +int4(ri_79,ps179(:,OP_DZ),CONJUGATE(bfp179(:,OP_DR)),temp79a) &
             -int4(ri_79,ps179(:,OP_DR),CONJUGATE(bfp179(:,OP_DZ)),temp79a) &
             +int4(ri_79,CONJUGATE(ps179(:,OP_DZ)),bfp179(:,OP_DR),temp79a) &
             -int4(ri_79,CONJUGATE(ps179(:,OP_DR)),bfp179(:,OP_DZ),temp79a))
     endif
#endif
  else
!    nonlinear:  do not subtract off equilibrium piece
     temp = .5* &
          (int4(ri2_79,pst79(:,OP_DZ),pst79(:,OP_DZ),temp79a) &
          +int4(ri2_79,pst79(:,OP_DR),pst79(:,OP_DR),temp79a)) ! &
!          - .5* &
!          (int4(ri2_79,ps079(:,OP_DZ),ps079(:,OP_DZ),temp79a) &
!          +int4(ri2_79,ps079(:,OP_DR),ps079(:,OP_DR),temp79a))
#if defined(USE3D)
     if(numvar.gt.1) then
        temp = temp   &
             + .5* &
             (int3(bfpt79(:,OP_DZ),bfpt79(:,OP_DZ),temp79a) &
             +int3(bfpt79(:,OP_DR),bfpt79(:,OP_DR),temp79a) &
             +2.*int4(ri_79,pst79(:,OP_DZ),bfpt79(:,OP_DR),temp79a) &
             -2.*int4(ri_79,pst79(:,OP_DR),bfpt79(:,OP_DZ),temp79a) )
     endif
#endif

#ifdef USECOMPLEX
     if(numvar.gt.1) then
        temp = temp + .5* &
             (int3(bfpt79(:,OP_DZ),CONJUGATE(bfpt79(:,OP_DZ)),temp79a) &
             +int3(bfpt79(:,OP_DR),CONJUGATE(bfpt79(:,OP_DR)),temp79a) &
             +int4(ri_79,pst79(:,OP_DZ),CONJUGATE(bfpt79(:,OP_DR)),temp79a) &
             -int4(ri_79,pst79(:,OP_DR),CONJUGATE(bfpt79(:,OP_DR)),temp79a) &
             +int4(ri_79,CONJUGATE(pst79(:,OP_DZ)),bfpt79(:,OP_DR),temp79a) &
             -int4(ri_79,CONJUGATE(pst79(:,OP_DR)),bfpt79(:,OP_DZ),temp79a))
     endif
#endif
  endif

  energy_mp = temp
  return
end function energy_mp


! Toroidal magnetic
! -----------------
real function energy_mt()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = .5*int3(ri2_79,bz179(:,OP_1),CONJUGATE(bz179(:,OP_1)))
  else
!....nonlinear:  do not subtract off equilibrium piece
     temp = .5*int3(ri2_79,bzt79(:,OP_1),bzt79(:,OP_1))!   &
!          - .5*int3(ri2_79,bz079(:,OP_1),bz079(:,OP_1))
  endif

  energy_mt = temp
  return
end function energy_mt


! Pressure
! --------
real function energy_p(mask)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), optional, dimension(MAX_PTS) :: mask

  vectype :: temp

  if(present(mask)) then
     temp79a = mask
  else
     temp79a = 1.
  end if

  if(gam.le.1.) then 
     temp = 0.
  else
     if(linear.eq.1) then
        temp = int2(p179,temp79a) / (gam - 1.)
     else
!.......nonlinear: subtract off equilibrium piece
        !temp = (int2(pt79,temp79a) - int2(p079,temp79a))/ (gam - 1.)
        temp = int2(pt79,temp79a) / (gam - 1.)
     endif
  endif

  energy_p = temp
  return
end function energy_p


! Electron Pressure
! -----------------
real function energy_pe()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gam.le.1.) then 
     temp = 0.
  else
     if(linear.eq.1) then
        temp = int1(pe179) / (gam - 1.)
     else
!.......nonlinear: subtract off equilibrium piece
        temp = int1(pet79) / (gam - 1.)
     endif
  endif

  energy_pe = temp
  return
end function energy_pe



! Poloidal kinetic
! ----------------
real function energy_kp()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = .5* &
             (int4(r2_79,ph179(:,OP_DZ),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
             +int4(r2_79,ph179(:,OP_DR),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)))
     else
        temp = .5* &
             (int4(r2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),rho79(:,OP_1)) &
             +int4(r2_79,pht79(:,OP_DR),pht79(:,OP_DR),rho79(:,OP_1)))
     endif

  energy_kp = temp
  return
end function energy_kp


! Toroidal kinetic
! ----------------
real function energy_kt()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = .5*int4(r2_79,vz179(:,OP_1),CONJUGATE(vz179(:,OP_1)),n079(:,OP_1))
     else
        temp = .5*int4(r2_79,vzt79(:,OP_1),vzt79(:,OP_1),rho79(:,OP_1))
     endif

  energy_kt = temp
  return
end function energy_kt


! Compressional kinetic
! ---------------------
real function energy_k3()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = .5* &
          (int4(ri4_79,ch179(:,OP_DZ),CONJUGATE(ch179(:,OP_DZ)),n079(:,OP_1)) &
          +int4(ri4_79,ch179(:,OP_DR),CONJUGATE(ch179(:,OP_DR)),n079(:,OP_1)) &
          +int4(ri_79,ch179(:,OP_DZ),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)) &
          -int4(ri_79,ch179(:,OP_DR),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
          +int4(ri_79,CONJUGATE(ch179(:,OP_DZ)),ph179(:,OP_DR),n079(:,OP_1)) &
          -int4(ri_79,CONJUGATE(ch179(:,OP_DR)),ph179(:,OP_DZ),n079(:,OP_1)))
     else
        temp = .5* &
          (int4(ri4_79,cht79(:,OP_DZ),cht79(:,OP_DZ),rho79(:,OP_1)) &
          +int4(ri4_79,cht79(:,OP_DR),cht79(:,OP_DR),rho79(:,OP_1)) &
          +int4(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),rho79(:,OP_1)) &
          -int4(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),rho79(:,OP_1)) &
          +int4(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),rho79(:,OP_1)) &
          -int4(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),rho79(:,OP_1)))
     endif

  energy_k3 = temp
  return
end function energy_k3


! Poloidal resistive
! ------------------
real function energy_mpd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = -int4(ri2_79,ps179(:,OP_GS),CONJUGATE(ps179(:,OP_GS)),eta79(:,OP_1))
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,ps179(:,OP_DZP),CONJUGATE(ps179(:,OP_DZP)),eta79(:,OP_1)) &
          +int4(ri4_79,ps179(:,OP_DRP),CONJUGATE(ps179(:,OP_DRP)),eta79(:,OP_1)))
#endif
  else
     temp = -int4(ri2_79,pst79(:,OP_GS),CONJUGATE(pst79(:,OP_GS)),eta79(:,OP_1))
  endif

  energy_mpd = temp
  return
end function energy_mpd


! Toroidal resistive
! ------------------
real function energy_mtd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = - &
          (int4(ri2_79,bz179(:,OP_DZ),CONJUGATE(bz179(:,OP_DZ)),eta79(:,OP_1))&
          +int4(ri2_79,bz179(:,OP_DR),CONJUGATE(bz179(:,OP_DR)),eta79(:,OP_1)))
  else
     temp = - &
          (int4(ri2_79,bzt79(:,OP_DZ),CONJUGATE(bzt79(:,OP_DZ)),eta79(:,OP_1))&
          +int4(ri2_79,bzt79(:,OP_DR),CONJUGATE(bzt79(:,OP_DR)),eta79(:,OP_1)))
  end if

  energy_mtd = temp
  return
end function energy_mtd


! Poloidal viscous
! ----------------
real function energy_kpd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = -int4(ri2_79,ph179(:,OP_GS),CONJUGATE(ph179(:,OP_GS)),vis79(:,OP_1))
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,ph179(:,OP_DZP),CONJUGATE(ph179(:,OP_DZP)),vis79(:,OP_1)) &
          +int4(ri4_79,ph179(:,OP_DRP),CONJUGATE(ph179(:,OP_DRP)),vis79(:,OP_1)))
#endif
  else
     temp = -int4(ri2_79,pht79(:,OP_GS),CONJUGATE(pht79(:,OP_GS)),vis79(:,OP_1))
  endif

  energy_kpd = temp
  return
end function energy_kpd


! Toroidal viscous
! ----------------
real function energy_ktd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = - &
             (int4(r2_79,vz179(:,OP_DZ),CONJUGATE(vz179(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(r2_79,vz179(:,OP_DR),CONJUGATE(vz179(:,OP_DR)),vis79(:,OP_1)))
     else
        temp = - &
             (int4(r2_79,vzt79(:,OP_DZ),CONJUGATE(vzt79(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(r2_79,vzt79(:,OP_DR),CONJUGATE(vzt79(:,OP_DR)),vis79(:,OP_1)))
     endif

  energy_ktd = temp
  return
end function energy_ktd

! Compressional viscous
! ---------------------
real function energy_k3d()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = - 2.*int3(ch179(:,OP_LP),CONJUGATE(ch179(:,OP_LP)),vic79(:,OP_1))
  else
     temp = - 2.*int3(cht79(:,OP_LP),CONJUGATE(cht79(:,OP_LP)),vic79(:,OP_1))
  endif

  energy_k3d = temp
  return
end function energy_k3d


! Poloidal hyper-viscous
! ----------------------
real function energy_kph()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

!!$  if(hypc.ne.0.) then
!!$     temp = - hypc* &
!!$          (int4(ri2_79,vot79(:,OP_DZ),CONJUGATE(vot79(:,OP_DZ)),vis79(:,OP_1)) &
!!$          +int4(ri2_79,vot79(:,OP_DR),CONJUGATE(vot79(:,OP_DR)),vis79(:,OP_1)))
!!$  else
     temp = 0.
!!$  end if
  energy_kph = temp
  return
end function energy_kph


! Toroidal hyper-viscous
! ----------------------
real function energy_kth()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = -hypv*int4(r2_79,vzt79(:,OP_GS),CONJUGATE(vzt79(:,OP_GS)),vis79(:,OP_1))

  energy_kth = temp
  return
end function energy_kth

! Compressional hyper-viscous
! ---------------------------
real function energy_k3h()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

!!$  if(hypc.ne.0.) then
!!$     temp = -2.*hypc* &
!!$          (int3(cot79(:,OP_DZ),CONJUGATE(cot79(:,OP_DZ)),vic79(:,OP_1)) &
!!$          +int3(cot79(:,OP_DR),CONJUGATE(cot79(:,OP_DR)),vic79(:,OP_1)))
!!$  else
     temp = 0.
!!$  end if
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
  use m3dc1_nint

  implicit none

  vectype :: temp

  if((numvar.lt.3 .and. ipres.eq.0) .or. gam.eq.1. .or. gam.eq.0.) then
     flux_pressure = 0.
     return
  endif

  temp = 0.5* &
       (int4(r_79,pt79(:,OP_1),norm79(:,2),pht79(:,OP_DR)) &
       -int4(r_79,pt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
       +int4(ri2_79,pt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
       +int4(ri2_79,pt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ)))

  if(db .ne. 0.) then
     temp = temp + db* &
          (int5(ri_79,pet79(:,OP_1),ni79(:,OP_1),norm79(:,1),bzt79(:,OP_DZ)) &
          -int5(ri_79,pet79(:,OP_1),ni79(:,OP_1),norm79(:,2),bzt79(:,OP_DR)))
  endif

  flux_pressure = gam*real(temp)/(gam-1.)
end function flux_pressure


! Kinetic Energy Convection
! -------------------------
real function flux_ke()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp79a = r2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2) &
       +    r2_79*vzt79(:,OP_1)**2 &
       +   ri4_79*(cht79(:,OP_DZ)**2 + cht79(:,OP_DR)**2) &
       + 2.*ri_79* &
       (cht79(:,OP_DZ)*pht79(:,OP_DR)-cht79(:,OP_DR)*pht79(:,OP_DZ))

  temp = 0.5* &
       (int5(r_79,nt79(:,OP_1),norm79(:,2),pht79(:,OP_DR),temp79a) &
       -int5(r_79,nt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ),temp79a) &
       +int5(ri2_79,nt79(:,OP_1),norm79(:,1),cht79(:,OP_DR),temp79a) &
       +int5(ri2_79,nt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ),temp79a))

  flux_ke = real(temp)
end function flux_ke


! Poynting flux
! -------------
real function flux_poynting()

  use math
  use basic
  use m3dc1_nint

  implicit none
  vectype :: temp

  temp = -vloop/toroidal_period * &
       (int3(ri2_79,norm79(:,1),pst79(:,OP_DR)) &
       +int3(ri2_79,norm79(:,2),pst79(:,OP_DZ)))

  flux_poynting = real(temp)
end function flux_poynting


! Heat flux
! ---------
real function flux_heat()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp
  
  if(numvar.lt.3 .and. ipres.eq.0) then
     flux_heat = 0.
     return
  endif

!!$  temp = int4(kap79(:,OP_1),norm79(:,1),pt79(:,OP_DR),ni79(:,OP_1)) &
!!$       + int4(kap79(:,OP_1),norm79(:,2),pt79(:,OP_DZ),ni79(:,OP_1)) &
!!$       + int4(kap79(:,OP_1),norm79(:,1),pt79(:,OP_1),ni79(:,OP_DR)) &
!!$       + int4(kap79(:,OP_1),norm79(:,2),pt79(:,OP_1),ni79(:,OP_DZ))
!!$
!!$  if(kappar.ne.0.) then
!!$     temp79a = ni79(:,OP_1)* &
!!$          (pt79(:,OP_DZ)*pst79(:,OP_DR) - pt79(:,OP_DR)*pst79(:,OP_DZ)) &
!!$          +    pt79(:,OP_1)* &
!!$          (ni79(:,OP_DZ)*pst79(:,OP_DR) - ni79(:,OP_DR)*pst79(:,OP_DZ))
!!$     temp79b = norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR)
!!$     temp = temp &
!!$          + int5(ri2_79,kar79(:,OP_1),b2i79(:,OP_1),temp79a,temp79b)
!!$  endif

  temp = -int3(kap79(:,OP_1),norm79(:,1),tet79(:,OP_DR)) &
       - int3(kap79(:,OP_1),norm79(:,2),tet79(:,OP_DZ))
  temp = temp &
       - kappai_fac*(  int3(kap79(:,OP_1),norm79(:,1),tit79(:,OP_DR)) &
                     + int3(kap79(:,OP_1),norm79(:,2),tit79(:,OP_DZ)))

  if(kappar.ne.0.) then
     temp79a = (tet79(:,OP_DZ)*pst79(:,OP_DR)-tet79(:,OP_DR)*pst79(:,OP_DZ))
     temp79c = (tit79(:,OP_DZ)*pst79(:,OP_DR)-tit79(:,OP_DR)*pst79(:,OP_DZ))
     temp79b = norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR)
     temp = temp &
          + int5(ri2_79,kar79(:,OP_1),b2i79(:,OP_1),temp79a,temp79b) 
     temp = temp &
          + int5(ri2_79,kar79(:,OP_1),b2i79(:,OP_1),temp79c,temp79b)
  endif

  flux_heat = real(temp)
  return
end function flux_heat

! Grav_pot
! --------
real function grav_pot()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     grav_pot = 0.
     return
  endif

  temp = gravr*int3(ri3_79,pht79(:,OP_DZ),nt79(:,OP_1)) &
       - gravz*int3( ri_79,pht79(:,OP_DR),nt79(:,OP_1))
     
  if(numvar.ge.3) then
     temp = -gravr*int3(ri2_79,cht79(:,OP_DR),nt79(:,OP_1)) &
          -gravz*int2(       cht79(:,OP_DZ),nt79(:,OP_1))
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
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = int4(ri_79,bzt79(:,OP_1),norm79(:,2),pst79(:,OP_DR)) &
       - int4(ri_79,bzt79(:,OP_1),norm79(:,1),pst79(:,OP_DZ))
  
  torque_em = temp
end function torque_em

! torque_sol
! ~~~~~~~~~~
vectype function torque_sol()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = int5(r3_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
       - int5(r3_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),pht79(:,OP_DR))
  
  torque_sol = temp
end function torque_sol

! torque_com
! ~~~~~~~~~~
vectype function torque_com()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(numvar.lt.3) then
     torque_com = 0.
     return
  end if

  temp = &
       - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
       - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ))     

  torque_com = temp
end function torque_com


! torque_visc
! ~~~~~~~~~~~
vectype function torque_visc()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = int4(r2_79,vis79(:,OP_1),norm79(:,1),vzt79(:,OP_DR)) &
       + int4(r2_79,vis79(:,OP_1),norm79(:,2),vzt79(:,OP_DZ))

  torque_visc = temp
end function torque_visc


! torque_parvisc
! ~~~~~~~~~~~~~~
vectype function torque_parvisc()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp
  
  if(amupar.eq.0.) then
     torque_parvisc = 0.
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

  temp79d = 3.*vip79(:,OP_1)*bzt79(:,OP_1)*b2i79(:,OP_1)*temp79a
  temp = int4(ri_79,temp79d,norm79(:,2),pst79(:,OP_DR)) &
       - int4(ri_79,temp79d,norm79(:,1),pst79(:,OP_DZ))

  torque_parvisc = temp
end function torque_parvisc



! torque_gyro
! ~~~~~~~~~~~
vectype function torque_gyro()
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gyro.eq.0 .or. numvar.lt.2) then
     torque_gyro = 0.
     return
  endif

  temp79a = 0.25*db*b2i79(:,OP_1)
  temp79b = temp79a*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)

  temp79c = norm79(:,1)*pst79(:,OP_DZ) + norm79(:,2)*pst79(:,OP_DR)
  temp79d = norm79(:,1)*pst79(:,OP_DR) - norm79(:,2)*pst79(:,OP_DZ)
  temp79e = 3.*temp79a*b2i79(:,OP_1)* &
       (norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR))
  temp79f = pht79(:,OP_DZZ) - pht79(:,OP_DRR)
  if(itor.eq.1) then
     temp79f = temp79f - ri_79*pht79(:,OP_DR)
  endif

  ! U contribution
  temp = int4(r_79,temp79b,temp79c,temp79f) &
       +2.*int4(r_79,temp79b,temp79d,pht79(:,OP_DRZ)) &
       +   int5(ri_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),temp79f) &
       -   int5(ri_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),temp79f) &
       +4.*int5(ri_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),pht79(:,OP_DRZ))

  if(itor.eq.1) then
     temp = temp &
          + 2.*int4(temp79b,norm79(:,1),pst79(:,OP_DR),pht79(:,OP_DZ)) &
          - 2.*int4(temp79a,pht79(:,OP_DZ),norm79(:,1),pst79(:,OP_DR)) &
          - 2.*int4(temp79a,pht79(:,OP_DZ),norm79(:,2),pst79(:,OP_DZ)) &
          + 2.*int5(ri2_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),pht79(:,OP_DZ))
  endif
     
  ! omega contribution
  temp = temp + 2.* &
       (int5(r_79,temp79b,bzt79(:,OP_1),norm79(:,1),vzt79(:,OP_DZ)) &
       -int5(r_79,temp79b,bzt79(:,OP_1),norm79(:,2),vzt79(:,OP_DR)))

  temp79f = cht79(:,OP_DZZ) - cht79(:,OP_DRR)
  if(itor.eq.1) then
     temp79f = temp79f + 2.*ri_79*cht79(:,OP_DR)
  endif

  ! chi constribution
  temp = temp &
       +    int4(ri2_79,temp79b,temp79d,temp79f) &
       - 2.*int4(ri2_79,temp79b,temp79c,cht79(:,OP_DRZ)) &
       +    int5(ri2_79,temp79b,cht79(:,OP_GS),norm79(:,1),pst79(:,OP_DR)) &
       +    int5(ri2_79,temp79b,cht79(:,OP_GS),norm79(:,2),pst79(:,OP_DZ)) &
       - 2.*int4(ri2_79,cht79(:,OP_GS),norm79(:,1),pst79(:,OP_DR)) &
       - 2.*int4(ri2_79,cht79(:,OP_GS),norm79(:,2),pst79(:,OP_DZ)) &
       + 2.*int5(ri4_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),temp79f) &
       - 2.*int5(ri4_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),cht79(:,OP_DRZ)) &
       + 2.*int5(ri4_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),cht79(:,OP_DRZ))

  if(itor.eq.1) then
     temp = temp &
          + 2.*int4(ri3_79,temp79b,temp79c,cht79(:,OP_DZ)) &
          + 3.*int5(ri3_79,temp79b,cht79(:,OP_DR),norm79(:,1),pst79(:,OP_DR)) &
          + 3.*int5(ri3_79,temp79b,cht79(:,OP_DR),norm79(:,2),pst79(:,OP_DZ)) &
          - 6.*int4(ri3_79,cht79(:,OP_DR),norm79(:,1),pst79(:,OP_DR)) &
          - 6.*int4(ri3_79,cht79(:,OP_DR),norm79(:,2),pst79(:,OP_DZ)) &
          + 2.*int5(ri5_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),cht79(:,OP_DZ)) &
          - 2.*int5(ri5_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),cht79(:,OP_DZ))
  endif

  torque_gyro = 0.
end function torque_gyro


! torque_denm
! ~~~~~~~~~~~
vectype function torque_denm()
  use basic
  use m3dc1_nint

  implicit none

  torque_denm = 0.

end function torque_denm

! Volume of parallel diffusion (Paul, Hudson & Helander, JPP 2022) 
! ---------
real function volume_pd(mask)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), optional, dimension(MAX_PTS) :: mask

  vectype :: temp

  temp = 0.
#ifdef USE3D
  if(kappar.ne.0.) then
    if(present(mask)) then
      temp79a = mask
    else
      temp79a = 1.
    end if
    temp79b = kar79(:, OP_1) * b2i79 (:, OP_1) * &
              (ri_79 * (tet79(:, OP_DZ) * pst79(:, OP_DR) &
                       -tet79(:, OP_DR) * pst79(:, OP_DZ))&
              + ri2_79 * bzt79(:, OP_1) * tet79(:, OP_DP) & 
               - (tet79(:, OP_DZ) * bfpt79(:, OP_DZ)      &
                 +tet79(:, OP_DR) * bfpt79(:, OP_DR)))**2 &
    ! Assuming that B0.grad(T0) = 0
    !          (ri_79 * (te179(:, OP_DZ) * pst79(:, OP_DR) &
    !                   -te179(:, OP_DR) * pst79(:, OP_DZ) &
    !                   +te079(:, OP_DZ) * ps179(:, OP_DR) &
    !                   -te079(:, OP_DR) * ps179(:, OP_DZ))&
    !          + ri2_79 *(bzt79(:, OP_1) * te179(:, OP_DP) & 
    !                   + bz179(:, OP_1) * te079(:, OP_DP))& 
    !           - (te179(:, OP_DZ) * bfpt79(:, OP_DZ)      &
    !             +te179(:, OP_DR) * bfpt79(:, OP_DR)      &
    !             +te079(:, OP_DZ) * bfp179(:, OP_DZ)      &
    !             +te079(:, OP_DR) * bfp179(:, OP_DR)))**2 &
             -kap79(:, OP_1) * & 
               (tet79(:, OP_DR)**2 + tet79(:, OP_DZ)**2   &
               +tet79(:, OP_DP)**2 * ri2_79) 
    ! Heaviside function
    temp79c = sign(1., temp79b) * 0.5 + 0.5
    temp = int1(temp79c*temp79a)
  end if  
#endif
  volume_pd = real(temp)
  return
end function volume_pd 

function tepsipsikappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsipsikappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsipsikappar = 0.
     return
  end if

  temp79a = k(:,OP_1)*ri2_79*j(:,OP_1)

  temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DZ))

  tepsipsikappar = (gam - 1.) * temp
end function tepsipsikappar


function tepsibkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsibkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsibkappar = 0.
     return
  end if

  temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)*j(:,OP_1)  

  temp79b = f(:,OP_DR)*(h(:,OP_DZ) ) &
       -    f(:,OP_DZ)*(h(:,OP_DR) )

  temp79d = temp79b*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1)

  temp = intx4(e(:,:,OP_DZ),f(:,OP_DR),temp79a,h(:,OP_DP)) &
       - intx4(e(:,:,OP_DR),f(:,OP_DZ),temp79a,h(:,OP_DP)) &
       - intx3(e(:,:,OP_DP),ri3_79,temp79d)

  tepsibkappar = (gam - 1.) * temp
#else
  tepsibkappar = 0.
#endif
end function tepsibkappar

function tepsibkapparl(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsibkapparl
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsibkapparl = 0.
     return
  end if

  temp79a = i(:,OP_1)*j(:,OP_1)*g(:,OP_1)

  ! d(temp79a)/dphi
  temp79b = i(:,OP_DP)*j(:,OP_1 )*g(:,OP_1 ) &
       +    i(:,OP_1 )*j(:,OP_DP)*g(:,OP_1 ) &
       +    i(:,OP_1 )*j(:,OP_1 )*g(:,OP_DP)

  ! [T,psi]
  temp79c = (h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ))

  ! d(temp79c)/dphi
  temp79d =      &
       (h(:,OP_DZ)*f(:,OP_DRP)-h(:,OP_DR)*f(:,OP_DZP)) &
       + (h(:,OP_DZP)*f(:,OP_DR)-h(:,OP_DRP)*f(:,OP_DZ))

  temp = intx4(e(:,:,OP_1),ri3_79,temp79a,temp79d) &
       + intx4(e(:,:,OP_1),ri3_79,temp79b,temp79c) &
       + intx5(e(:,:,OP_DR),ri3_79,temp79a,f(:,OP_DZ),h(:,OP_DP)) &
       - intx5(e(:,:,OP_DZ),ri3_79,temp79a,f(:,OP_DR),h(:,OP_DP))
  tepsibkapparl = (gam - 1.) * temp
#else
  tepsibkapparl = 0.
#endif
end function tepsibkapparl


function tebbkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebbkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tebbkappar = 0.
     return
  end if

  temp79a = h(:,OP_DP)

  temp79c = f(:,OP_1)*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 )

  temp = -intx3(e(:,:,OP_DP),ri4_79,temp79c)

  tebbkappar = (gam - 1.) * temp
#else
  tebbkappar = 0.
#endif
end function tebbkappar
!
!...the following function must replace tebbkappar for linear runs.
function tebbkapparl(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebbkapparl
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  
  if(gam.le.1.) then
     tebbkapparl = 0.
     return
  end if

  temp79a = i(:,OP_1)*j(:,OP_1)*f(:,OP_1)*g(:,OP_1)

  ! d(temp79a)/dphi
  temp79b = i(:,OP_DP)*j(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
       +    i(:,OP_1)*j(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
       +    i(:,OP_1)*j(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
       +    i(:,OP_1)*j(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

  temp = intx4(e(:,:,OP_1),ri4_79,temp79a,h(:,OP_DPP)) &
       + intx4(e(:,:,OP_1),ri4_79,temp79b,h(:,OP_DP))

  tebbkapparl = (gam - 1.) * temp
#else
  tebbkapparl = 0.
#endif
end function tebbkapparl

function tepsifkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsifkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsifkappar = 0.
     return
  end if

  temp79a = k(:,OP_1)*ri_79*j(:,OP_1)

  temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
       + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
       - intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DZ ),h(:,OP_DR)) &
       - intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DZ ),h(:,OP_DR))

  tepsifkappar = (gam - 1.) * temp
#else
  tepsifkappar = 0.
#endif
end function tepsifkappar

function tebfkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebfkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tebfkappar = 0.
     return
  end if

  temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)*j(:,OP_1)

  temp79b = g(:,OP_DZ)*(h(:,OP_DZ) )&
       +    g(:,OP_DR)*(h(:,OP_DR) )

  temp79d = temp79b*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 )

  temp = intx4(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,h(:,OP_DP)) &
       + intx4(e(:,:,OP_DR),g(:,OP_DR),temp79a,h(:,OP_DP)) &
       + intx3(e(:,:,OP_DP),ri2_79,temp79d)

  tebfkappar = (gam - 1.) * temp
#else
  tebfkappar = 0.
#endif
end function tebfkappar

function tebfkapparl(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebfkapparl
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tebfkapparl = 0.
     return
  end if

  temp79a = ri2_79*j(:,OP_1)*i(:,OP_1)* f(:,OP_1)


  ! d(temp79a)/dphi
  temp79b = ri2_79 * &
       (j(:,OP_DP)*i(:,OP_1)* f(:,OP_1) &
       +j(:,OP_1)*i(:,OP_DP)* f(:,OP_1) &
       +j(:,OP_1)*i(:,OP_1)* f(:,OP_DP))

  !  <T, f'> 
  temp79c =  (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) 
    

  ! d(temp79c)/dphi
  temp79d =  &
       +(h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
       +(h(:,OP_DRP)*g(:,OP_DR) + h(:,OP_DZP)*g(:,OP_DZ)) 

  temp = intx4(e(:,:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DP)) &
       + intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DP)) &
       - intx3(e(:,:,OP_1),temp79a,temp79d) &
       - intx3(e(:,:,OP_1),temp79b,temp79c)

  tebfkapparl = (gam - 1.) * temp
#else
  tebfkapparl = 0.
#endif
  return
end function tebfkapparl


function teffkappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: teffkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     teffkappar = 0.
     return
  end if

  temp79a = -k(:,OP_1)*j(:,OP_1)

  temp = intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
       + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
       + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR))

  teffkappar = (gam - 1.) * temp
#else
  teffkappar = 0.
#endif
  return
end function teffkappar


function q_delta(e)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: q_delta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e

  q_delta = intx4(e(:,:,OP_1),net79(:,OP_1),tit79(:,OP_1),qd79) &
       -    intx4(e(:,:,OP_1),net79(:,OP_1),tet79(:,OP_1),qd79)
end function q_delta

function q_delta1(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: q_delta1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  q_delta1 = intx4(e(:,:,OP_1),f(:,OP_1),net79(:,OP_1),qd79) 
end function q_delta1

vectype function q1ppsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
       - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))

  q1ppsi =   temp
end function q1ppsi

vectype function q1pb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


#if defined(USE3D) || defined(USECOMPLEX)
  temp = int5(ri2_79,e(:,OP_1),f(:,OP_DP),g(:,OP_1),h(:,OP_1))
#else
  temp = 0.
#endif

  q1pb =   temp
end function q1pb

vectype function q1pf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype :: temp

  temp = - int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
       - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))

  q1pf = temp
#else
  q1pf = 0.
#endif
end function q1pf


function t3tn(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  temp = intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_1))

  t3tn = temp
end function t3tn

function t3t(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3t
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = intx2(e(:,:,OP_1),f(:,OP_1))

  t3t = temp
end function t3t

function t3tnu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tnu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = intx5(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_1),h(:,OP_DZ)) &
       - intx5(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_1),h(:,OP_DR))
  if(itor.eq.1) then
     temp = temp + &
          2.*(gam-1.)*intx4(e(:,:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_DZ))
  endif

  t3tnu = temp
end function t3tnu


function t3tnv(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tnv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = - intx4(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1),h(:,OP_1)) &
       - (gam-1.)*intx4(e(:,:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_DP))
#else
  temp = 0.
#endif

  t3tnv = temp
end function t3tnv

function t3tnchi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tnchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_1),h(:,OP_DR))  &
       -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_1),h(:,OP_DZ)) &
       -(gam-1.)* &
       intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_1),h(:,OP_GS))

  t3tnchi = temp
end function t3tnchi

vectype function j1b2ipsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

  temp = - int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_GS),h(:,OP_1))

  j1b2ipsib = temp
  return
end function j1b2ipsib
vectype function j1b2ibpsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

  temp = int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DR))    &
       + int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DZ))

  j1b2ibpsi = temp
  return
end function j1b2ibpsi

vectype function j1b2ipsif(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp
#if defined(USE3D) || defined(USECOMPLEX)
  temp = int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DRP))    &
       + int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DZP))    &
       - int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DRP),h(:,OP_DR))    &
       - int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZP),h(:,OP_DZ))
#else
  temp = 0
#endif

  j1b2ipsif = temp
  return
end function j1b2ipsif

vectype function j1b2ifb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = - int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR))    &
         + int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ))
#else
  temp = 0
#endif


  j1b2ifb = temp
  return
end function j1b2ifb

vectype function j1b2iff(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = - int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DRP))    &
         + int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZP))
#else
  temp = 0
#endif

  j1b2iff = temp
  return
end function j1b2iff

vectype function j1b2ipsipsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = + int5(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZP),h(:,OP_DR))    &
         - int5(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DRP),h(:,OP_DZ))
#else
  temp = 0
#endif

  j1b2ipsipsi = temp
  return
end function j1b2ipsipsi


vectype function pparpu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  temp = int4(r_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR)) &
       - int4(r_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ))

  pparpu = temp
end function pparpu

vectype function pparpupsipsib2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

  temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
  temp = - 2.*int5(ri_79,temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DR))  &
       + 2.*int5(ri_79,temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DZ))
  temp = temp                                                    &
       + int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DR)) &
       + int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DR)) &
       + int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
       + int5(ri_79,temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DR)) &
       - int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DZ)) &
       - int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DZ)) &
       - int5(ri_79,temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DZ)) &
       - int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
  temp = temp                                                    &
       -2*int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DRZ)) &
       -2*int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
       -2*int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DR),g(:,OP_DZZ)) &
       -2*int5(ri_79,temp79a,h(:,OP_DZZ),i(:,OP_DR),g(:,OP_DZ)) &
       +2*int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DRR)) &
       +2*int5(ri_79,temp79a,h(:,OP_DRR),i(:,OP_DZ),g(:,OP_DR)) &
       +2*int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DRZ)) &
       +2*int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
  if(itor.eq.1) then
     temp = temp                                                    &
          + 2.*int5(ri2_79,temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DZ)) &
          + 2.*int5(ri2_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DZ)) 
  endif

  pparpupsipsib2 = temp
end function pparpupsipsib2

vectype function pparpupsibb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

  temp = 0
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
  temp =  -2.*int5(ri2_79,temp79a,i(:,OP_1),g(:,OP_DRP),h(:,OP_DR))  &
       -2.*int5(ri2_79,temp79a,i(:,OP_1),g(:,OP_DZP),h(:,OP_DZ))
#else
  temp = 0
#endif

  pparpupsibb2 = temp
end function pparpupsibb2

vectype function pparpubbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

  temp = 0.
  if(itor.eq.1) then
     temp79a = 2.*e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
     temp = temp                                                    &
          + int5(ri2_79,temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DZ)) 
  endif

  pparpubbb2 = temp
end function pparpubbb2

! ===========
vectype function pparpvpsibb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

  temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
  temp = - 2.*int5(ri_79,temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1))  &
       + 2.*int5(ri_79,temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1))

  pparpvpsibb2 = temp
end function pparpvpsibb2

vectype function pparpvbbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp
  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp = -int3(e(:,OP_1),f(:,OP_1),g(:,OP_DP)) 
  temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
  temp = temp                                                    &
       - 2.*int5(ri2_79,temp79a,g(:,OP_DP),h(:,OP_1),i(:,OP_1))  
#endif

  pparpvbbb2 = temp
end function pparpvbbb2

vectype function pparpchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  temp = int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR)) &
       + int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ))

  pparpchi = temp
end function pparpchi

vectype function pparpchipsipsib2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

  temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
  temp = - 2.*int5(ri4_79,temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DR))  &
       - 2.*int5(ri4_79,temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DZ))
  temp = temp                                                     &
       + int5(ri4_79,temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DR)) &
       + int5(ri4_79,temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DR)) &
       + int5(ri4_79,temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DR)) &
       + int5(ri4_79,temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DR)) &
       + int5(ri4_79,temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DZ)) &
       + int5(ri4_79,temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DZ)) &
       + int5(ri4_79,temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DZ)) &
       + int5(ri4_79,temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DZ))
  temp = temp                                                       &
       -2.*int5(ri4_79,temp79a,g(:,OP_DZZ),h(:,OP_DR),i(:,OP_DR)) &
       -2.*int5(ri4_79,temp79a,g(:,OP_DZ),h(:,OP_DRZ),i(:,OP_DR)) &
       +2.*int5(ri4_79,temp79a,g(:,OP_DRZ),h(:,OP_DZ),i(:,OP_DR)) &
       +2.*int5(ri4_79,temp79a,g(:,OP_DR),h(:,OP_DZZ),i(:,OP_DR)) &
       +2.*int5(ri4_79,temp79a,g(:,OP_DRZ),h(:,OP_DR),i(:,OP_DZ)) &
       +2.*int5(ri4_79,temp79a,g(:,OP_DZ),h(:,OP_DRR),i(:,OP_DZ)) &
       -2.*int5(ri4_79,temp79a,g(:,OP_DRR),h(:,OP_DZ),i(:,OP_DZ)) &
       -2.*int5(ri4_79,temp79a,g(:,OP_DR),h(:,OP_DRZ),i(:,OP_DZ))
  if(itor.eq.1) then
     temp = temp                                                       &
          - 2.*int5(ri5_79,temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DR)) &
          - 2.*int5(ri5_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR))
     temp = temp                                                       &
          -6.*int5(ri5_79,temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DZ))  &
          +6.*int5(ri5_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR)) 
  endif

  pparpchipsipsib2 = temp
end function pparpchipsipsib2

vectype function pparpchipsibb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp
  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
  temp = temp                                                      &
       - 2.*int5(ri2_79,temp79a,g(:,OP_DZP),h(:,OP_DR),i(:,OP_1))  &
       + 2.*int5(ri2_79,temp79a,g(:,OP_DRP),h(:,OP_DZ),i(:,OP_1))  
#endif

  pparpchipsibb2 = temp
end function pparpchipsibb2

vectype function pparpchibbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

  temp = 0
  if(itor.eq.1) then
     temp79a = 2.*e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
     temp = temp                                                    &
          - int5(ri5_79,temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DR)) 
  endif

  pparpchibbb2 = temp
end function pparpchibbb2


function pperpu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = 2.*intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
       - 2.*intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ)) &
       + intx4(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_DR)) &
       - intx4(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_DZ))

  pperpu = temp
end function pperpu

function pperpupsipsib2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpupsipsib2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp79a = ri_79*f(:,OP_1)*j(:,OP_1)
  temp =   intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DR))  &
       -   intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DZ))
  temp = temp                                                    &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DR)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DR)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DR)) &
       + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DZ)) &
       + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DZ)) &
       + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DZ)) &
       + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
  temp = temp                                                    &
       + intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DRZ)) &
       + intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
       + intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DR),g(:,OP_DZZ)) &
       + intx5(e(:,:,OP_1),temp79a,h(:,OP_DZZ),i(:,OP_DR),g(:,OP_DZ)) &
       - intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DRR)) &
       - intx5(e(:,:,OP_1),temp79a,h(:,OP_DRR),i(:,OP_DZ),g(:,OP_DR)) &
       - intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DRZ)) &
       - intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
  if(itor.eq.1) then
     temp79a = ri2_79*f(:,OP_1)*j(:,OP_1)

     temp = temp                                                    &
          - intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DZ)) &
          - intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DZ)) 
  endif
  
  pperpupsipsib2 = temp
end function pperpupsipsib2

function pperpupsibb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpupsibb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = ri2_79*f(:,OP_1)*j(:,OP_1)
  temp = temp                                                      &
       + intx5(e(:,:,OP_1),temp79a,g(:,OP_DRP),h(:,OP_DR),i(:,OP_1))  &
       + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZP),h(:,OP_DZ),i(:,OP_1))  
#endif

  pperpupsibb2 = temp
end function pperpupsibb2


function pperpubbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpubbb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp = 0
  if(itor.eq.1) then
     temp79a =   ri2_79*f(:,OP_1)*j(:,OP_1)
     temp = temp                                                    &
          - intx5(e(:,:,OP_1),temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DZ)) 
  endif

  pperpubbb2 = temp
end function pperpubbb2

! ===========
function pperpvpsibb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpvpsibb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp79a =   ri_79*f(:,OP_1)*j(:,OP_1)
  temp =   intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1))  &
       - intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1))

  pperpvpsibb2 = temp
end function pperpvpsibb2

function pperpvbbb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpvbbb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp = -2.*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DP)) 
  temp79a =   ri2_79*f(:,OP_1)*j(:,OP_1)
  temp = temp                                       &
       + intx5(e(:,:,OP_1),temp79a, g(:,OP_DP),h(:,OP_1),i(:,OP_1))
#endif

  pperpvbbb2 = temp
end function pperpvbbb2

function pperpchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchi

  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  temp = 2.*intx4(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DR)) &
       + 2.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZ)) &
       + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR)) &
       + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ))

  pperpchi = temp
end function pperpchi


function pperpchipsipsib2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchipsipsib2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp79a =   ri4_79*f(:,OP_1)*j(:,OP_1)
  temp =   intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DR))  &
       + intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DZ))
  temp = temp                                                    &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DR)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DR)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DR)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DR)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DZ)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DZ)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DZ)) &
       - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DZ))

  temp = temp                                                       &
       + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZZ),h(:,OP_DR),i(:,OP_DR)) &
       + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DRZ),i(:,OP_DR)) &
       - intx5(e(:,:,OP_1),temp79a,g(:,OP_DRZ),h(:,OP_DZ),i(:,OP_DR)) &
       - intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DZZ),i(:,OP_DR)) &
       - intx5(e(:,:,OP_1),temp79a,g(:,OP_DRZ),h(:,OP_DR),i(:,OP_DZ)) &
       - intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DRR),i(:,OP_DZ)) &
       + intx5(e(:,:,OP_1),temp79a,g(:,OP_DRR),h(:,OP_DZ),i(:,OP_DZ)) &
       + intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DRZ),i(:,OP_DZ))
  if(itor.eq.1) then
     temp79a =   ri5_79*f(:,OP_1)*j(:,OP_1)
     temp = temp                                                    &
          + intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DR)) &
          + intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR))
     temp = temp                                                       &
          +3.*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DZ))  &
          -3.*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR))  
  endif

  pperpchipsipsib2 = temp
end function pperpchipsipsib2


function pperpchipsibb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchipsibb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp
  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = ri2_79*f(:,OP_1)*j(:,OP_1)
  temp = temp                                                      &
       + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZP),h(:,OP_DR),i(:,OP_1))  &
       - intx5(e(:,:,OP_1),temp79a,g(:,OP_DRP),h(:,OP_DZ),i(:,OP_1))  
#endif

  pperpchipsibb2 = temp
end function pperpchipsibb2

function pperpchibbb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchibbb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp = 0
  if(itor.eq.1) then
     temp79a =  ri5_79*f(:,OP_1)*j(:,OP_1)
     temp = temp                                                    &
          + intx5(e(:,:,OP_1),temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DR)) 
  endif

  pperpchibbb2 = temp
end function pperpchibbb2


function incvb(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: incvb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  incvb = intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_1))
end function incvb

function incupsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: incupsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  incupsi = intx3(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR)) &
       +    intx3(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ))
end function incupsi

function incchipsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: incchipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  incchipsi = intx4(e(:,:,OP_1),ri3_79,f(:,OP_DZ),g(:,OP_DR)) &
       -      intx4(e(:,:,OP_1),ri3_79,f(:,OP_DR),g(:,OP_DZ))
end function incchipsi

subroutine JxB_r(o, opol)
  use m3dc1_nint

  implicit none

  vectype, intent(out), dimension(MAX_PTS) :: o, opol
  vectype, dimension(MAX_PTS) :: otor

  otor = - ri2_79*pst79(:,OP_GS)*pstx79(:,OP_DR)
  opol = - ri2_79*bztx79(:,OP_1)*bzt79(:,OP_DR)
#if defined(USE3D) || defined(USECOMPLEX)
  otor = otor + ri_79*pst79(:,OP_GS)*bfptx79(:,OP_DZ)
  opol = opol &
       - ri2_79*bztx79(:,OP_1)*bfpt79(:,OP_DRP) &
       - ri3_79*bztx79(:,OP_1)*pst79(:,OP_DZP)
#endif
  o = opol + otor
end subroutine JxB_r

subroutine JxB_phi(o)
  use m3dc1_nint

  implicit none
  
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = &
       + ri2_79*bzt79(:,OP_DZ)*pstx79(:,OP_DR) &
       - ri2_79*bzt79(:,OP_DR)*pstx79(:,OP_DZ)

#if defined(USE3D) || defined(USECOMPLEX)
  o = o &
       + ri2_79*bfpt79(:,OP_DZP)*pstx79(:,OP_DR) &
       - ri2_79*bfpt79(:,OP_DRP)*pstx79(:,OP_DZ) &
       - ri_79*bzt79(:,OP_DZ)*bfptx79(:,OP_DZ) &
       - ri_79*bzt79(:,OP_DR)*bfptx79(:,OP_DR) &
       - ri_79*bfpt79(:,OP_DZP)*bfptx79(:,OP_DZ) &
       - ri_79*bfpt79(:,OP_DRP)*bfptx79(:,OP_DR) &
       - ri3_79*pst79(:,OP_DZP)*pstx79(:,OP_DZ) &
       - ri3_79*pst79(:,OP_DRP)*pstx79(:,OP_DR) &
       - ri2_79*pst79(:,OP_DZP)*bfptx79(:,OP_DR) &
       + ri2_79*pst79(:,OP_DRP)*bfptx79(:,OP_DZ)
#endif       

end subroutine JxB_phi
       
subroutine JxB_z(o, opol)
  use m3dc1_nint

  implicit none

  vectype, intent(out), dimension(MAX_PTS) :: o, opol
  vectype, dimension(MAX_PTS) :: otor

  otor = - ri2_79*pst79(:,OP_GS)*pstx79(:,OP_DZ)
  opol = - ri2_79*bztx79(:,OP_1)*bzt79(:,OP_DZ)

#if defined(USE3D) || defined(USECOMPLEX)
  otor = otor - ri_79*pst79(:,OP_GS)*bfptx79(:,OP_DR)
  opol = opol  &
       - ri2_79*bztx79(:,OP_1)*bfpt79(:,OP_DZP) &
       + ri3_79*bztx79(:,OP_1)*pst79(:,OP_DRP)
#endif
  o = opol + otor
  
end subroutine JxB_z

end module metricterms_new

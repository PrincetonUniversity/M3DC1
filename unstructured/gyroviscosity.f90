module gyroviscosity
  implicit none

  vectype, private, dimension(MAX_PTS) :: psir, psiz, gbb, gzz, grr, grz
  vectype, private, dimension(MAX_PTS) :: mucross

  private :: gyro_RR_u, gyro_RR_v, gyro_RR_x
  private :: gyro_RP_u, gyro_RP_v, gyro_RP_x
  private :: gyro_RZ_u, gyro_RZ_v, gyro_RZ_x
  private :: gyro_PP_u, gyro_PP_v, gyro_PP_x
  private :: gyro_ZP_u, gyro_ZP_v, gyro_ZP_x
  private :: gyro_ZZ_u, gyro_ZZ_v, gyro_ZZ_x

contains

  subroutine gyro_common()
    use m3dc1_nint

    implicit none

    psir = pst79(:,OP_DR)
    psiz = pst79(:,OP_DZ)
#ifdef USE3D
    ! f' = 0 in complex case, so this is only needed for USE3D
    psir = psir - r_79*bft79(:,OP_DZP)
    psiz = psiz + r_79*bft79(:,OP_DRP)
#endif

    gbb = 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2
    grr = 3.*ri2_79*b2i79(:,OP_1)*psir**2
    gzz = 3.*ri2_79*b2i79(:,OP_1)*psiz**2
    grz = 3.*ri2_79*b2i79(:,OP_1)*psir*psiz

    mucross = 0.25*pit79(:,OP_1)*b2i79(:,OP_1)
  end subroutine gyro_common

  !======================================================================
  ! Vorticity contributions
  !======================================================================
  vectype function gyro_vor_u(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RZ_u(f,temp79f)
    gyro_vor_u = int3(r_79,e(:,OP_DZZ),temp79f) &
         -       int3(r_79,e(:,OP_DRR),temp79f)
    if(itor.eq.1) gyro_vor_u = gyro_vor_u - int2(e(:,OP_DR),temp79f)
    
    call gyro_RR_u(f,temp79f)
    gyro_vor_u = gyro_vor_u + int3(r_79,e(:,OP_DRZ),temp79f)
    if(itor.eq.1) gyro_vor_u = gyro_vor_u + int2(e(:,OP_DZ),temp79f)

    call gyro_ZZ_u(f,temp79f)
    gyro_vor_u = gyro_vor_u - int3(r_79,e(:,OP_DRZ),temp79f)

    if(itor.eq.1) then 
       call gyro_PP_u(f,temp79f)
       gyro_vor_u = gyro_vor_u + int2(e(:,OP_DZ),temp79f)
    end if

#ifdef USE3D
    call gyro_RP_u(f,temp79f)
    gyro_vor_u = gyro_vor_u + int2(e(:,OP_DZP),temp79f)
    call gyro_ZP_u(f,temp79f)
    gyro_vor_u = gyro_vor_u - int2(e(:,OP_DRP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_RP_u(f,temp79f)
    gyro_vor_u = gyro_vor_u - rfac*int2(e(:,OP_DZ),temp79f)
    call gyro_ZP_u(f,temp79f)
    gyro_vor_u = gyro_vor_u + rfac*int2(e(:,OP_DR),temp79f)
#endif

    gyro_vor_u = -gyro_vor_u
  end function gyro_vor_u

  vectype function gyro_vor_v(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RZ_v(f,temp79f)
    gyro_vor_v = int3(r_79,e(:,OP_DZZ),temp79f) &
         -       int3(r_79,e(:,OP_DRR),temp79f)
    if(itor.eq.1) gyro_vor_v = gyro_vor_v - int2(e(:,OP_DR),temp79f)
    
    call gyro_RR_v(f,temp79f)
    gyro_vor_v = gyro_vor_v + int3(r_79,e(:,OP_DRZ),temp79f)
    if(itor.eq.1) gyro_vor_v = gyro_vor_v + int2(e(:,OP_DZ),temp79f)
     
    call gyro_ZZ_v(f,temp79f)
    gyro_vor_v = gyro_vor_v - int3(r_79,e(:,OP_DRZ),temp79f)

    if(itor.eq.1) then 
       call gyro_PP_v(f,temp79f)
       gyro_vor_v = gyro_vor_v + int2(e(:,OP_DZ),temp79f)
    end if

#ifdef USE3D
    call gyro_RP_v(f,temp79f)
    gyro_vor_v = gyro_vor_v + int2(e(:,OP_DZP),temp79f)
    call gyro_ZP_v(f,temp79f)
    gyro_vor_v = gyro_vor_v - int2(e(:,OP_DRP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_RP_v(f,temp79f)
    gyro_vor_v = gyro_vor_v - rfac*int2(e(:,OP_DZ),temp79f)
    call gyro_ZP_v(f,temp79f)
    gyro_vor_v = gyro_vor_v + rfac*int2(e(:,OP_DR),temp79f)
#endif

    gyro_vor_v = -gyro_vor_v
  end function gyro_vor_v

  vectype function gyro_vor_x(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RZ_x(f,temp79f)
    gyro_vor_x = int3(r_79,e(:,OP_DZZ),temp79f) &
         -       int3(r_79,e(:,OP_DRR),temp79f)
    if(itor.eq.1) gyro_vor_x = gyro_vor_x - int2(e(:,OP_DR),temp79f)
    
    call gyro_RR_x(f,temp79f)
    gyro_vor_x = gyro_vor_x + int3(r_79,e(:,OP_DRZ),temp79f)
    if(itor.eq.1) gyro_vor_x = gyro_vor_x + int2(e(:,OP_DZ),temp79f)
     
    call gyro_ZZ_x(f,temp79f)
    gyro_vor_x = gyro_vor_x - int3(r_79,e(:,OP_DRZ),temp79f)

    if(itor.eq.1) then 
       call gyro_PP_x(f,temp79f)
       gyro_vor_x = gyro_vor_x + int2(e(:,OP_DZ),temp79f)
    end if

#ifdef USE3D
    call gyro_RP_x(f,temp79f)
    gyro_vor_x = gyro_vor_x + int2(e(:,OP_DZP),temp79f)
    call gyro_ZP_x(f,temp79f)
    gyro_vor_x = gyro_vor_x - int2(e(:,OP_DRP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_RP_x(f,temp79f)
    gyro_vor_x = gyro_vor_x - rfac*int2(e(:,OP_DZ),temp79f)
    call gyro_ZP_x(f,temp79f)
    gyro_vor_x = gyro_vor_x + rfac*int2(e(:,OP_DR),temp79f)
#endif

    gyro_vor_x = -gyro_vor_x
  end function gyro_vor_x

  !======================================================================
  ! Toroidal contributions
  !======================================================================
  vectype function gyro_tor_u(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RP_u(f,temp79f)
    gyro_tor_u = int3(r_79,e(:,OP_DR),temp79f)
    call gyro_ZP_u(f,temp79f)
    gyro_tor_u = gyro_tor_u + int3(r_79,e(:,OP_DZ),temp79f)

#ifdef USE3D
    call gyro_PP_u(f,temp79f)
    gyro_tor_u = gyro_tor_u + int2(e(:,OP_DP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_PP_u(f,temp79f)
    gyro_tor_u = gyro_tor_u - rfac*int2(e(:,OP_1),temp79f)
#endif
  end function gyro_tor_u

  vectype function gyro_tor_v(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RP_v(f,temp79f)
    gyro_tor_v = int3(r_79,e(:,OP_DR),temp79f)
    call gyro_ZP_v(f,temp79f)
    gyro_tor_v = gyro_tor_v + int3(r_79,e(:,OP_DZ),temp79f)

#ifdef USE3D
    call gyro_PP_v(f,temp79f)
    gyro_tor_v = gyro_tor_v + int2(e(:,OP_DP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_PP_v(f,temp79f)
    gyro_tor_v = gyro_tor_v - rfac*int2(e(:,OP_1),temp79f)
#endif
  end function gyro_tor_v

  vectype function gyro_tor_x(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RP_x(f,temp79f)
    gyro_tor_x = int3(r_79,e(:,OP_DR),temp79f)
    call gyro_ZP_x(f,temp79f)
    gyro_tor_x = gyro_tor_x + int3(r_79,e(:,OP_DZ),temp79f)

#ifdef USE3D
    call gyro_PP_x(f,temp79f)
    gyro_tor_x = gyro_tor_x + int2(e(:,OP_DP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_PP_x(f,temp79f)
    gyro_tor_x = gyro_tor_x - rfac*int2(e(:,OP_1),temp79f)
#endif
  end function gyro_tor_x

  !======================================================================
  ! Compressional contributions
  !======================================================================
  vectype function gyro_com_u(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RR_u(f,temp79f)
    gyro_com_u = -int3(ri2_79,e(:,OP_DRR),temp79f)
    if(itor.eq.1) gyro_com_u = gyro_com_u + 2.*int3(ri3_79,e(:,OP_DR),temp79f)

    call gyro_ZZ_u(f,temp79f)
    gyro_com_u = gyro_com_u - int3(ri2_79,e(:,OP_DZZ),temp79f)

    call gyro_RZ_u(f,temp79f)
    gyro_com_u = gyro_com_u - 2.*int3(ri2_79,e(:,OP_DRZ),temp79f)
    if(itor.eq.1) gyro_com_u = gyro_com_u + 2.*int3(ri3_79,e(:,OP_DZ),temp79f)

    call gyro_PP_u(f,temp79f)
    gyro_com_u = gyro_com_u - int3(ri3_79,e(:,OP_DR),temp79f)

#ifdef USE3D
    call gyro_RP_u(f,temp79f)
    gyro_com_u = gyro_com_u - int3(ri3_79,e(:,OP_DRP),temp79f)
    call gyro_ZP_u(f,temp79f)
    gyro_com_u = gyro_com_u - int3(ri3_79,e(:,OP_DZP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_RP_u(f,temp79f)
    gyro_com_u = gyro_com_u + rfac*int3(ri3_79,e(:,OP_DR),temp79f)
    call gyro_ZP_u(f,temp79f)
    gyro_com_u = gyro_com_u + rfac*int3(ri3_79,e(:,OP_DZ),temp79f)
#endif
  end function gyro_com_u

  vectype function gyro_com_v(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RR_v(f,temp79f)
    gyro_com_v = -int3(ri2_79,e(:,OP_DRR),temp79f)
    if(itor.eq.1) gyro_com_v = gyro_com_v + 2.*int3(ri3_79,e(:,OP_DR),temp79f)

    call gyro_ZZ_v(f,temp79f)
    gyro_com_v = gyro_com_v - int3(ri2_79,e(:,OP_DZZ),temp79f)

    call gyro_RZ_v(f,temp79f)
    gyro_com_v = gyro_com_v - 2.*int3(ri2_79,e(:,OP_DRZ),temp79f)
    if(itor.eq.1) gyro_com_v = gyro_com_v + 2.*int3(ri3_79,e(:,OP_DZ),temp79f)

    call gyro_PP_v(f,temp79f)
    gyro_com_v = gyro_com_v - int3(ri3_79,e(:,OP_DR),temp79f)

#ifdef USE3D
    call gyro_RP_v(f,temp79f)
    gyro_com_v = gyro_com_v - int3(ri3_79,e(:,OP_DRP),temp79f)
    call gyro_ZP_v(f,temp79f)
    gyro_com_v = gyro_com_v - int3(ri3_79,e(:,OP_DZP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_RP_v(f,temp79f)
    gyro_com_v = gyro_com_v + rfac*int3(ri3_79,e(:,OP_DR),temp79f)
    call gyro_ZP_v(f,temp79f)
    gyro_com_v = gyro_com_v + rfac*int3(ri3_79,e(:,OP_DZ),temp79f)
#endif
  end function gyro_com_v

  vectype function gyro_com_x(e,f)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

    call gyro_RR_x(f,temp79f)
    gyro_com_x = -int3(ri2_79,e(:,OP_DRR),temp79f)
    if(itor.eq.1) gyro_com_x = gyro_com_x + 2.*int3(ri3_79,e(:,OP_DR),temp79f)

    call gyro_ZZ_x(f,temp79f)
    gyro_com_x = gyro_com_x - int3(ri2_79,e(:,OP_DZZ),temp79f)

    call gyro_RZ_x(f,temp79f)
    gyro_com_x = gyro_com_x - 2.*int3(ri2_79,e(:,OP_DRZ),temp79f)
    if(itor.eq.1) gyro_com_x = gyro_com_x + 2.*int3(ri3_79,e(:,OP_DZ),temp79f)

    call gyro_PP_x(f,temp79f)
    gyro_com_x = gyro_com_x - int3(ri3_79,e(:,OP_DR),temp79f)

#ifdef USE3D
    call gyro_RP_x(f,temp79f)
    gyro_com_x = gyro_com_x - int3(ri3_79,e(:,OP_DRP),temp79f)
    call gyro_ZP_x(f,temp79f)
    gyro_com_x = gyro_com_x - int3(ri3_79,e(:,OP_DZP),temp79f)
#elif defined(USECOMPLEX)
    call gyro_RP_x(f,temp79f)
    gyro_com_x = gyro_com_x + rfac*int3(ri3_79,e(:,OP_DR),temp79f)
    call gyro_ZP_x(f,temp79f)
    gyro_com_x = gyro_com_x + rfac*int3(ri3_79,e(:,OP_DZ),temp79f)
#endif
  end function gyro_com_x



  !======================================================================
  ! Pi_RR
  ! ~~~~~
  !======================================================================
  subroutine gyro_RR_u(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! u_ZZ - u_RR - u_R / R
    temp79a = f(:,OP_DZZ) - f(:,OP_DRR)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DR)

    ! u_RZ + u_Z / R
    temp79b = f(:,OP_DRZ)
    if(itor.eq.1) temp79b = temp79b + ri_79*f(:,OP_DZ)

    o = -(1.+gzz)*bzt79(:,OP_1)*temp79a - grz*2.*bzt79(:,OP_1)*temp79b

#if defined(USE3D) || defined(USECOMPLEX)
    o = o + ri_79*psir*f(:,OP_DZP) &
         + grz*ri_79*(psir*f(:,OP_DRP) + psiz*f(:,OP_DZP)) &
         - gbb*ri_79*psiz*f(:,OP_DRP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_RR_u

  subroutine gyro_RR_v(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    o = -psir*f(:,OP_DR) - gbb*psiz*f(:,OP_DZ) + &
         grz*(psir*f(:,OP_DZ) - psiz*f(:,OP_DR))

#if defined(USE3D) || defined(USECOMPLEX)
    o = o + grz*2.*ri_79*bzt79(:,OP_1)*f(:,OP_DP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_RR_v

  subroutine gyro_RR_x(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! chi_RZ - chi_Z / R
    temp79a = f(:,OP_DRZ)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DZ)
    
    ! chi_ZZ - chi_R / R
    temp79b = f(:,OP_DZZ)
    if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

    o = 2.*ri3_79*bzt79(:,OP_1)*((1.+gzz)*temp79a - grz*temp79b)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - ri4_79*psir*f(:,OP_DRP) &
         + grz*ri4_79*(psir*f(:,OP_DZP) - psiz*f(:,OP_DRP)) &
         - gbb*ri4_79*psiz*f(:,OP_DZP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_RR_x


  !======================================================================
  ! Pi_RP
  ! ~~~~~
  !======================================================================
  subroutine gyro_RP_u(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! u_ZZ - u_RR - u_R / R
    temp79a = f(:,OP_DZZ) - f(:,OP_DRR)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DR)

    ! u_RZ + u_Z / R
    temp79b = f(:,OP_DRZ)
    if(itor.eq.1) temp79b = temp79b + ri_79*f(:,OP_DZ)

    o = (1.-gzz)*psiz*temp79a &
         - (1.+gzz)*psir*f(:,OP_DRZ) &
         + (gbb-gzz)*psir*temp79b

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - ri_79*(1.-gbb)*bzt79(:,OP_1)*f(:,OP_DRP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_RP_u

  subroutine gyro_RP_v(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    o = -(1.-gbb)*bzt79(:,OP_1)*f(:,OP_DZ)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - (1.+gbb)*ri_79*psir*f(:,OP_DP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_RP_v

  subroutine gyro_RP_x(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! chi_RZ - chi_Z / R
    temp79a = f(:,OP_DRZ)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DZ)
    
    ! chi_ZZ - chi_R / R
    temp79b = f(:,OP_DZZ)
    if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

    ! chi_RR - 3*chi_R / R
    temp79c = f(:,OP_DRR)
    if(itor.eq.1) temp79c = temp79c - 3.*ri_79*f(:,OP_DR)

    o = ri3_79*(1.+gzz)*psir*temp79c &
         - 2.*ri3_79*(1.-gzz)*psiz*temp79a &
         + ri3_79*(gbb-gzz)*psir*temp79b

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - ri4_79*(1.-gbb)*bzt79(:,OP_1)*f(:,OP_DZP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_RP_x


  !======================================================================
  ! Pi_RZ
  ! ~~~~~
  !======================================================================
  subroutine gyro_RZ_u(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! u_RZ + u_Z / R
    temp79b = f(:,OP_DRZ)
    if(itor.eq.1) temp79b = temp79b + ri_79*f(:,OP_DZ)

    o = 2.*bzt79(:,OP_1)*((1.+gzz)*f(:,OP_DRZ) + (1.+grr)*temp79b)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o + (1.-gbb)*ri_79*(psiz*f(:,OP_DZP) - psir*f(:,OP_DRP)) &
         + (gzz-grr)*ri_79*(psir*f(:,OP_DRP) + psiz*f(:,OP_DZP))
#endif
    
    o = o*mucross
  end subroutine gyro_RZ_u

  subroutine gyro_RZ_v(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    o = -(1.-gbb)*(psiz*f(:,OP_DR) + psir*f(:,OP_DZ)) &
         + (gzz-grr)*(psir*f(:,OP_DZ) - psiz*f(:,OP_DR))

#if defined(USE3D) || defined(USECOMPLEX)
    o = o + (gzz-grr)*2.*ri_79*bzt79(:,OP_1)*f(:,OP_DP)
#endif
    o = o*mucross
  end subroutine gyro_RZ_v

  subroutine gyro_RZ_x(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! chi_ZZ - chi_R / R
    temp79b = f(:,OP_DZZ)
    if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

    ! chi_RR - 3*chi_R / R
    temp79c = f(:,OP_DRR)
    if(itor.eq.1) temp79c = temp79c - 3.*ri_79*f(:,OP_DR)

    o = 2.*ri3_79*bzt79(:,OP_1)* &
         ((1.+grr)*temp79b - (1.+gzz)*temp79c)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - ri4_79*(1.-gbb)*(psiz*f(:,OP_DRP) + psir*f(:,OP_DZP)) &
         + ri4_79*(gzz-grr)*(psir*f(:,OP_DZP) - psiz*f(:,OP_DRP))
#endif
    o = o*mucross
  end subroutine gyro_RZ_x

  !======================================================================
  ! Pi_PP
  ! ~~~~~
  !======================================================================
  subroutine gyro_PP_u(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! u_ZZ - u_RR - u_R / R
    temp79a = f(:,OP_DZZ) - f(:,OP_DRR)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DR)

    ! 2*u_RZ + u_Z / R
    temp79b = 2.*f(:,OP_DRZ)
    if(itor.eq.1) temp79b = temp79b + ri_79*f(:,OP_DZ)

    o = (gzz-grr)*bzt79(:,OP_1)*temp79a + grz*2.*bzt79(:,OP_1)*temp79b

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - (1.+gbb)*ri_79*(psir*f(:,OP_DZP) - psiz*f(:,OP_DRP))
#endif

    o = 2.*mucross*o
  end subroutine gyro_PP_u

  subroutine gyro_PP_v(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    o = (1.+gbb)*(psir*f(:,OP_DR) + psiz*f(:,OP_DZ))

    o = 2.*mucross*o
  end subroutine gyro_PP_v

  subroutine gyro_PP_x(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! chi_RZ - chi_Z / R
    temp79a = f(:,OP_DRZ)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DZ)
    
    ! chi_ZZ - chi_R / R
    temp79b = f(:,OP_DZZ)
    if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

    ! chi_RR - 3*chi_R / R
    temp79c = f(:,OP_DRR)
    if(itor.eq.1) temp79c = temp79c - 3.*ri_79*f(:,OP_DR)

    o = 2.*ri3_79*bzt79(:,OP_1)*(-(gzz-grr)*temp79a + grz*(temp79b-temp79c))

#if defined(USE3D) || defined(USECOMPLEX)
    o = o + (1.+gbb)*ri4_79*(psir*f(:,OP_DRP) + psiz*f(:,OP_DZP))
#endif

    o = 2.*mucross*o
  end subroutine gyro_PP_x


  !======================================================================
  ! Pi_ZP
  ! ~~~~~
  !======================================================================
  subroutine gyro_ZP_u(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! u_ZZ - u_RR - u_R / R
    temp79a = f(:,OP_DZZ) - f(:,OP_DRR)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DR)

    ! u_RZ + u_Z / R
    temp79b = f(:,OP_DRZ)
    if(itor.eq.1) temp79b = temp79b + ri_79*f(:,OP_DZ)

    o = (1.-grr)*psir*temp79a + (1.+grr)*psiz*temp79b &
         -(gbb-grr)*psiz*f(:,OP_DRZ)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - (1.-gbb)*ri_79*bzt79(:,OP_1)*f(:,OP_DZP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_ZP_u

  subroutine gyro_ZP_v(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    o = (1.-gbb)*bzt79(:,OP_1)*f(:,OP_DR)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - (1.+gbb)*ri_79*psiz*f(:,OP_DP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_ZP_v

  subroutine gyro_ZP_x(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! chi_RZ - chi_Z / R
    temp79a = f(:,OP_DRZ)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DZ)
    
    ! chi_ZZ - chi_R / R
    temp79b = f(:,OP_DZZ)
    if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

    ! chi_RR - 3*chi_R / R
    temp79c = f(:,OP_DRR)
    if(itor.eq.1) temp79c = temp79c - 3.*ri_79*f(:,OP_DR)

    o = ri3_79*((1.+grr)*psiz*temp79b + (gbb-grr)*psiz*temp79c &
         - 2.*(1.-grr)*psir*temp79a)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o + (1.-gbb)*ri4_79*bzt79(:,OP_1)*f(:,OP_DRP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_ZP_x


  !======================================================================
  ! Pi_ZZ
  ! ~~~~~
  !======================================================================
  subroutine gyro_ZZ_u(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! u_ZZ - u_RR - u_R / R
    temp79a = f(:,OP_DZZ) - f(:,OP_DRR)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DR)

    o = (1.+grr)*bzt79(:,OP_1)*temp79a &
         - grz*2.*bzt79(:,OP_1)*f(:,OP_DRZ)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - ri_79*psiz*f(:,OP_DRP) &
         - ri_79*grz*(psir*f(:,OP_DRP) + psiz*f(:,OP_DZP)) &
         + gbb*ri_79*psir*f(:,OP_DZP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_ZZ_u

  subroutine gyro_ZZ_v(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    o = -psiz*f(:,OP_DZ) - grz*(psir*f(:,OP_DZ) - psiz*f(:,OP_DR)) &
         - gbb*psir*f(:,OP_DR)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - 2.*ri_79*grz*bzt79(:,OP_1)*f(:,OP_DP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_ZZ_v

  subroutine gyro_ZZ_x(f,o)
    use basic
    use m3dc1_nint
    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, intent(out), dimension(MAX_PTS) :: o

    if(surface_int) then
       o = 0.
       return
    end if

    ! chi_RZ - chi_Z / R
    temp79a = f(:,OP_DRZ)
    if(itor.eq.1) temp79a = temp79a - ri_79*f(:,OP_DZ)
    
    ! chi_RR - 3*chi_R / R
    temp79c = f(:,OP_DRR)
    if(itor.eq.1) temp79c = temp79c - 3.*ri_79*f(:,OP_DR)

    o = 2.*ri3_79*bzt79(:,OP_1)*(-(1.+grr)*temp79a + grz*temp79c)

#if defined(USE3D) || defined(USECOMPLEX)
    o = o - ri4_79*psiz*f(:,OP_DZP) &
         - ri4_79*grz*(psir*f(:,OP_DZP) - psiz*f(:,OP_DRP)) &
         - ri4_79*gbb*psir*f(:,OP_DRP)
#endif

    o = 2.*mucross*o
  end subroutine gyro_ZZ_x

end module gyroviscosity

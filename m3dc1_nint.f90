!======================================================================
! m3dc1_nint
! ~~~~~~~~~~
! This module defines arrays of the values of fields at numerical
! integration quadrature sampling points, and routines for populating
! these arrays.
!======================================================================
module m3dc1_nint

  use nintegrate

  implicit none

  ! The following give the meaning of the value array returned by local_value
  integer, parameter :: OP_1    = 1
  integer, parameter :: OP_DR   = 2
  integer, parameter :: OP_DZ   = 3
  integer, parameter :: OP_DRR  = 4
  integer, parameter :: OP_DRZ  = 5
  integer, parameter :: OP_DZZ  = 6
  integer, parameter :: OP_LP   = 7
  integer, parameter :: OP_GS   = 8
  integer, parameter :: OP_NUM_POL = 8
#if defined(USECOMPLEX) || defined(USE3D) 
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
  integer, parameter :: OP_LPR   = 25
  integer, parameter :: OP_LPZ   = 26
  integer, parameter :: OP_NUM   = 26
#else
  integer, parameter :: OP_LPR  = 9
  integer, parameter :: OP_LPZ  = 10
  integer, parameter :: OP_NUM  = 10
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
  integer, parameter :: FIELD_MU  = 65536

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element) :: mu79, nu79
  vectype, dimension(MAX_PTS) :: r_79, r2_79, r3_79, &
     ri_79, ri2_79, ri3_79, ri4_79, ri5_79, ri6_79, ri7_79, ri8_79
  vectype, dimension(MAX_PTS) :: temp79a, temp79b, temp79c, &
       temp79d, temp79e, temp79f
  vectype, dimension(MAX_PTS, OP_NUM) :: sz79
  vectype, dimension(MAX_PTS, OP_NUM) :: tm79, ni79, b2i79
  vectype, dimension(MAX_PTS, OP_NUM) :: ps179, bz179, pe179, n179, & 
       ph179, vz179, ch179, p179
  vectype, dimension(MAX_PTS, OP_NUM) :: pst79, bzt79, pet79, nt79, &
       pht79, vzt79, cht79, pt79
  vectype, dimension(MAX_PTS, OP_NUM) :: bzx79, bzi79, psx79, psi79
  vectype, dimension(MAX_PTS, OP_NUM) :: vis79, vic79, vip79
  vectype, dimension(MAX_PTS, OP_NUM) :: jt79, cot79, vot79, pit79, &
       eta79, sig79
  vectype, dimension(MAX_PTS, OP_NUM) :: bf079, bf179, bft79, bfx79, bfi79
  vectype, dimension(MAX_PTS, OP_NUM) :: kap79, kar79, kax79
  vectype, dimension(MAX_PTS, OP_NUM) :: ps079, bz079, pe079, n079, &
       ph079, vz079, ch079, p079
  vectype, dimension(MAX_PTS, OP_NUM) :: pss79, bzs79

  ! precalculated terms
   real, private :: fterm(MAX_PTS, coeffs_per_element, OP_NUM)
  
contains

  !==================================================
  ! precalculate_terms
  ! ~~~~~~~~~~~~~~~~~~
  ! precalculates the values of each term in the
  ! finite element expansion at each sampling point
  !==================================================
  subroutine precalculate_terms(xi,zi,eta,co,sn,ri,npoints)
    use basic

    implicit none
      
    integer, intent(in) :: npoints
    real, dimension(MAX_PTS), intent(in) :: xi, zi, eta
    real, intent(in) :: co, sn
    vectype, dimension(MAX_PTS), intent(in) :: ri

    integer :: p
    real, dimension(MAX_PTS) :: temp
    real :: co2, sn2, cosn
    real :: xpow(MAX_PTS,-3:5), ypow(MAX_PTS,-3:5)
#ifdef USE3D
    integer :: i, j, op
    real :: zpow(MAX_PTS,-2:3)
#endif

    co2 = co*co
    sn2 = sn*sn
    cosn = co*sn

    xpow(:,-3:-1) = 0.
    ypow(:,-3:-1) = 0.
    xpow(:,0) = 1.
    ypow(:,0) = 1.

#ifdef USE3D
    zpow(:,-2:-1) = 0.
    zpow(:,0) = 1.
#endif
 
    do p=1, 5
       xpow(:,p) = xpow(:,p-1)*xi(:)
       ypow(:,p) = ypow(:,p-1)*eta(:)
    end do
#ifdef USE3D
    do p=1, 3
       zpow(:,p) = zpow(:,p-1)*zi(:)
    end do
#endif
    
    fterm = 0.
    do p=1, coeffs_per_tri
       fterm(:,p,OP_1) = xpow(:,mi(p))*ypow(:,ni(p))
          
       if(mi(p).ge.1) then
          ! d_si terms
          temp = mi(p)*xpow(:,mi(p)-1) * ypow(:,ni(p))
          fterm(:,p,OP_DR) = fterm(:,p,OP_DR) + co*temp
          fterm(:,p,OP_DZ) = fterm(:,p,OP_DZ) + sn*temp           
             
          if(mi(p).ge.2) then
             ! d_si^2 terms
             temp = xpow(:,mi(p)-2)*(mi(p)-1)*mi(p) * ypow(:,ni(p))
             fterm(:,p,OP_DRR) = fterm(:,p,OP_DRR) + co2*temp
             fterm(:,p,OP_DZZ) = fterm(:,p,OP_DZZ) + sn2*temp
             fterm(:,p,OP_DRZ) = fterm(:,p,OP_DRZ) + cosn*temp
             fterm(:,p,OP_LP) = fterm(:,p,OP_LP) + temp
          endif
       endif
       if(ni(p).ge.1) then
          ! d_eta terms
          temp = xpow(:,mi(p)) * ypow(:,ni(p)-1)*ni(p)
          fterm(:,p,OP_DR) = fterm(:,p,OP_DR) - sn*temp
          fterm(:,p,OP_DZ) = fterm(:,p,OP_DZ) + co*temp
          
          if(ni(p).ge.2) then
             ! d_eta^2 terms
             temp = xpow(:,mi(p)) * ypow(:,ni(p)-2)*(ni(p)-1)*ni(p)
             fterm(:,p,OP_DRR) = fterm(:,p,OP_DRR) + sn2*temp
             fterm(:,p,OP_DZZ) = fterm(:,p,OP_DZZ) + co2*temp
             fterm(:,p,OP_DRZ) = fterm(:,p,OP_DRZ) - cosn*temp
             fterm(:,p,OP_LP) = fterm(:,p,OP_LP) + temp
          endif
          
          if(mi(p).ge.1) then
             ! d_eta_si terms
             temp = xpow(:,mi(p)-1)*mi(p) * ypow(:,ni(p)-1)*ni(p)
             
             fterm(:,p,OP_DRR) = fterm(:,p,OP_DRR) - 2.*cosn*temp
             fterm(:,p,OP_DZZ) = fterm(:,p,OP_DZZ) + 2.*cosn*temp
             fterm(:,p,OP_DRZ) = fterm(:,p,OP_DRZ) + (co2-sn2)*temp
          endif
       endif
       
       ! for surface terms, higher derivatives may be taken
       if(surface_int) then
          if(mi(p).ge.2) then
             if(ni(p).ge.1) then
                ! d_si^2 d_eta terms
                temp = xpow(:,mi(p)-2)*ypow(:,ni(p)-1)*(mi(p)-1)*mi(p)*ni(p)
                fterm(:,p,OP_LPR) = fterm(:,p,OP_LPR) - sn*temp
                fterm(:,p,OP_LPZ) = fterm(:,p,OP_LPZ) + co*temp
             endif
          endif
          if(ni(p).ge.2) then
             if(mi(p).ge.1) then
                ! d_eta^2 d_si terms
                temp = xpow(:,mi(p)-1)*ypow(:,ni(p)-2)*mi(p)*(ni(p)-1)*ni(p)
                fterm(:,p,OP_LPR) = fterm(:,p,OP_LPR) + co*temp
                fterm(:,p,OP_LPZ) = fterm(:,p,OP_LPZ) + sn*temp
             endif
          endif
          
          if(mi(p).ge.3) then
             ! d_si^3 terms
             temp = xpow(:,mi(p)-3)*ypow(:,ni(p))*(mi(p)-2)*(mi(p)-1)*mi(p)
             fterm(:,p,OP_LPR) = fterm(:,p,OP_LPR) + co*temp
             fterm(:,p,OP_LPZ) = fterm(:,p,OP_LPZ) + sn*temp
          endif
          if(ni(p).ge.3) then
             ! d_eta^3 terms
             temp = xpow(:,mi(p))*ypow(:,ni(p)-3)*(ni(p)-2)*(ni(p)-1)*ni(p)
             fterm(:,p,OP_LPR) = fterm(:,p,OP_LPR) - sn*temp
             fterm(:,p,OP_LPZ) = fterm(:,p,OP_LPZ) + co*temp
          endif
       endif
       
       ! Grad-Shafranov operator, and
       ! cylindrical correction to Laplacian
       fterm(:,p,OP_GS) = fterm(:,p,OP_LP)
       if(itor.eq.1) then
          fterm(:,p,OP_GS) = fterm(:,p,OP_GS) - fterm(:,p,OP_DR)*ri(:)
          fterm(:,p,OP_LP) = fterm(:,p,OP_LP) + fterm(:,p,OP_DR)*ri(:)
          
          if(surface_int) then
             fterm(:,p,OP_LPR) = fterm(:,p,OP_LPR) + fterm(:,p,OP_DRR)*ri(:) &
                  - fterm(:,p,OP_DR)*ri(:)*ri(:)
             fterm(:,p,OP_LPZ) = fterm(:,p,OP_LPZ) + fterm(:,p,OP_DRZ)*ri(:)
          endif
       endif      

#ifdef USE3D
       do op=1, OP_NUM_POL
          temp(:) = fterm(:,p,op)

          do i=1, coeffs_per_dphi
             j = p + (i-1)*coeffs_per_tri

             fterm(:,j,op) = temp(:)*zpow(:,li(i))

             ! first toroidal derivative
             if(li(i).ge.1) then
                fterm(:,j,op+OP_NUM_POL) = temp(:) &
                     *zpow(:,li(i)-1)*li(i)
             endif
             ! second toroidal derivative
             if(li(i).ge.2) then
                fterm(:,j,op+2*OP_NUM_POL) = temp(:) &
                     *zpow(:,li(i)-2)*(li(i)-1)*li(i)
             endif
          end do
       end do
#endif
    end do
  end subroutine precalculate_terms

  !===============================================
  ! eval_ops
  ! --------
  !
  ! evaluates linear unitary operators
  !===============================================
  subroutine eval_ops(avector, npoints, outarr)
    implicit none
    
    integer, intent(in) :: npoints
    vectype, dimension(coeffs_per_element), intent(in) :: avector
    vectype, dimension(MAX_PTS, OP_NUM), intent(out) :: outarr

    integer :: p, i
    
    outarr = 0.
    do i=1, OP_NUM
       do p=1, coeffs_per_element
          outarr(:,i) = outarr(:,i) + avector(p)*fterm(:,p,i)
       end do
    end do
  end subroutine eval_ops

  !=====================================================
  ! define_fields
  !=====================================================
  subroutine define_fields(itri, fields, gdef, ilin)
    use basic
    use mesh_mod
    use arrays

    implicit none
  
    integer, intent(in) :: itri, fields, gdef, ilin

    real :: fac
    integer :: i
    type(element_data) :: d
    vectype, dimension(coeffs_per_element) :: avec
    vectype, dimension(dofs_per_element,coeffs_per_element) :: cl

    ! calculate the hyperviscosity coefficients and
    ! the size field for this element.
    if(ihypdx.eq.0) then
       fac = 1.
    else
       fac = deex**ihypdx
    endif
    hypf = hyper *fac
    hypi = hyperi*fac
    hypv = hyperv*fac
    hypc = hyperc*fac
    hypp = hyperp*fac

    call interpolate_size_field(itri)

    ! calculate the major radius, and useful powers
    call get_element_data(itri, d)
    call local_to_global(d, xi_79, zi_79, eta_79, x_79, phi_79, z_79)
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

    call precalculate_terms(xi_79,zi_79,eta_79,d%co,d%sn,ri_79,npoints)

    ! PHI
    ! ~~~
    if(iand(fields, FIELD_PHI).eq.FIELD_PHI) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   U..."
       
       if(ilin.eq.0) then
          call calcavector(itri, u_field(1), avec)
          call eval_ops(avec, npoints, ph179)
#ifdef USECOMPLEX
          ph179(:,OP_DP :OP_GSP ) = ph179(:,OP_1:OP_GS)*rfac
          ph179(:,OP_DPP:OP_GSPP) = ph179(:,OP_1:OP_GS)*rfac**2
#endif
       else
          ph179 = 0.
       endif
       
       if(eqsubtract.eq.1) then
          call calcavector(itri, u_field(0), avec)
          call eval_ops(avec, npoints, ph079)
          pht79 = ph079 + ph179
       else
          ph079 = 0.
          pht79 = ph179
       endif
    endif

    ! PSI
    ! ~~~
    if(iand(fields, FIELD_PSI).eq.FIELD_PSI) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   psi..."

       if(use_external_fields) then 
          call calcavector(itri, psi_ext, avec)
          call eval_ops(avec, npoints, psx79)
#ifdef USECOMPLEX
          psx79(:,OP_DP :OP_GSP ) = psx79(:,OP_1:OP_GS)*rfac
          psx79(:,OP_DPP:OP_GSPP) = psx79(:,OP_1:OP_GS)*rfac**2
#endif
       else
          psx79 = 0.
       end if
       
       if(ilin.eq.0) then
          call calcavector(itri, psi_field(1), avec)
          call eval_ops(avec, npoints, psi79)
#ifdef USECOMPLEX
          psi79(:,OP_DP :OP_GSP ) = psi79(:,OP_1:OP_GS)*rfac
          psi79(:,OP_DPP:OP_GSPP) = psi79(:,OP_1:OP_GS)*rfac**2
#endif
       else
          psi79 = 0.
       end if
       ps179 = psi79 + psx79
       
       if(eqsubtract.eq.1) then
          call calcavector(itri, psi_field(0), avec)
          call eval_ops(avec, npoints, ps079)
          pst79 = ps079 + ps179
          pss79 = ps079 + ps179/2.
       else
          ps079 = 0.
          pst79 = ps179
          pss79 = ps179/2.
       endif
    endif

    ! V
    ! ~
    if(iand(fields, FIELD_V).eq.FIELD_V) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   V..."
       
       if(ilin.eq.0) then
          call calcavector(itri, vz_field(1), avec)
          call eval_ops(avec, npoints, vz179)
#ifdef USECOMPLEX
          vz179(:,OP_DP :OP_GSP ) = vz179(:,OP_1:OP_GS)*rfac
          vz179(:,OP_DPP:OP_GSPP) = vz179(:,OP_1:OP_GS)*rfac**2
#endif
       else
          vz179 = 0.
       end if
       
       if(eqsubtract.eq.1) then
          call calcavector(itri, vz_field(0), avec)
          call eval_ops(avec, npoints, vz079)
          vzt79 = vz079 + vz179
       else
          vz079 = 0.
          vzt79 = vz179
       endif
    endif

    ! I
    ! ~
    if(iand(fields, FIELD_I).eq.FIELD_I) then
       
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   I..."
      
       if(use_external_fields) then
          call calcavector(itri, bz_ext, avec)
          call eval_ops(avec, npoints, bzx79)
#ifdef USECOMPLEX
          bzx79(:,OP_DP :OP_GSP ) = bzx79(:,OP_1:OP_GS)*rfac
          bzx79(:,OP_DPP:OP_GSPP) = bzx79(:,OP_1:OP_GS)*rfac**2
#endif
      
#if defined(USECOMPLEX) || defined(USE3D)    
          call calcavector(itri, bf_ext, avec)
          call eval_ops(avec, npoints, bfx79)
#endif
#ifdef USECOMPLEX
          bfx79(:,OP_DP :OP_GSP ) = bfx79(:,OP_1:OP_GS)*rfac
          bfx79(:,OP_DPP:OP_GSPP) = bfx79(:,OP_1:OP_GS)*rfac**2
#endif
       else
          bzx79 = 0.
          bfx79 = 0.
       endif
 
       if(ilin.eq.0) then
          call calcavector(itri, bz_field(1), avec)
          call eval_ops(avec, npoints, bzi79)
#ifdef USECOMPLEX
          bzi79(:,OP_DP :OP_GSP ) = bzi79(:,OP_1:OP_GS)*rfac
          bzi79(:,OP_DPP:OP_GSPP) = bzi79(:,OP_1:OP_GS)*rfac**2
#endif
      
#if defined(USECOMPLEX) || defined(USE3D)    
          call calcavector(itri, bf_field(1), avec)
          call eval_ops(avec, npoints, bfi79)
#endif
#ifdef USECOMPLEX
          bfi79(:,OP_DP :OP_GSP ) = bfi79(:,OP_1:OP_GS)*rfac
          bfi79(:,OP_DPP:OP_GSPP) = bfi79(:,OP_1:OP_GS)*rfac**2
#endif
       else
          bzi79 = 0.
          bfi79 = 0.
       endif
       bz179 = bzi79 + bzx79
       bf179 = bfi79 + bfx79
       
       if(eqsubtract.eq.1) then
          call calcavector(itri, bz_field(0), avec)
          call eval_ops(avec, npoints, bz079)
          bzt79 = bz079 + bz179
          bzs79 = bz079 + bz179/2.
          
#if defined(USECOMPLEX) || defined(USE3D)
          call calcavector(itri, bf_field(0), avec)
          call eval_ops(avec, npoints, bf079)
          bft79 = bf079 + bf179
#endif
       else
          bz079 = 0.
          bzt79 = bz179
          bzs79 = bz179/2.
          
#ifdef USECOMPLEX
          bf079 = 0.
          bft79 = bf179
#endif
       endif

       if(numvar.eq.1) bzs79 = bzt79
    endif

    ! CHI
    ! ~~~
    if(iand(fields, FIELD_CHI).eq.FIELD_CHI) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   chi..."
       
       if(ilin.eq.0) then
          call calcavector(itri, chi_field(1), avec)
          call eval_ops(avec, npoints, ch179)
#ifdef USECOMPLEX
          ch179(:,OP_DP :OP_GSP ) = ch179(:,OP_1:OP_GS)*rfac
          ch179(:,OP_DPP:OP_GSPP) = ch179(:,OP_1:OP_GS)*rfac**2
#endif
       else
          ch179 = 0.
       end if
       
       if(eqsubtract.eq.1) then
          call calcavector(itri, chi_field(0), avec)
          call eval_ops(avec, npoints, ch079)
          cht79 = ch079 + ch179
       else
          ch079 = 0.
          cht79 = ch179
       endif
    endif
    
    ! P & PE
    ! ~~~~~~
    if((iand(fields, FIELD_PE).eq.FIELD_PE) .or. &
         (iand(fields, FIELD_P).eq.FIELD_P)) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   P..."
       
       if(ilin.eq.0) then
          call calcavector(itri, p_field(1), avec)
          call eval_ops(avec, npoints, p179)
          call calcavector(itri, pe_field(1), avec)
          call eval_ops(avec, npoints, pe179)
#ifdef USECOMPLEX
          p179(:,OP_DP :OP_GSP ) = p179(:,OP_1:OP_GS)*rfac
          p179(:,OP_DPP:OP_GSPP) = p179(:,OP_1:OP_GS)*rfac**2
          pe179(:,OP_DP :OP_GSP ) = pe179(:,OP_1:OP_GS)*rfac
          pe179(:,OP_DPP:OP_GSPP) = pe179(:,OP_1:OP_GS)*rfac**2
#endif
       else
          p179 = 0.
          pe179 = 0.
       end if
       
       if(eqsubtract.eq.1) then
          call calcavector(itri, p_field(0), avec)
          call eval_ops(avec, npoints, p079)
          call calcavector(itri, pe_field(0), avec)
          call eval_ops(avec, npoints, pe079)

          pet79 = pe079 + pe179
          pt79  =  p079 +  p179
       else
          pe079 = 0.
          p079 = 0.
          pet79 = pe179
          pt79  =  p179
       endif
       pit79 = pt79 - pet79
    endif
    
    
    ! N
    ! ~
    if(iand(fields, FIELD_N).eq.FIELD_N) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   n..."
       
       if(ilin.eq.0) then
          call calcavector(itri, den_field(1), avec)
          call eval_ops(avec, npoints, n179)
#ifdef USECOMPLEX
          n179(:,OP_DP :OP_GSP ) = n179(:,OP_1:OP_GS)*rfac
          n179(:,OP_DPP:OP_GSPP) = n179(:,OP_1:OP_GS)*rfac**2
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
             call calcavector(itri, den_field(0), avec)
             call eval_ops(avec, npoints, n079)
          end if
          nt79 = n079 + n179
       else
          n079 = 0.
          nt79 = n179
       endif
    endif

  ! NI
  ! ~~
  if(iand(fields, FIELD_NI).eq.FIELD_NI) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   n^-1..."
     ni79(1:npoints,OP_1  ) = 1./nt79(1:npoints,OP_1)
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
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   j..."

     if(ilin.eq.0) then
        call calcavector(itri, jphi_field, avec)
        call eval_ops(avec, npoints, jt79)
     else
        jt79 = 0.
     end if
  endif

  ! VOR
  ! ~~~
  if(iand(fields, FIELD_VOR).eq.FIELD_VOR) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   vor..."

     if(ilin.eq.0) then
        call calcavector(itri, vor_field, avec)
        call eval_ops(avec, npoints, vot79)
     else
        vot79 = 0.
     end if
  endif

  ! COM
  ! ~~~
  if(iand(fields, FIELD_COM).eq.FIELD_COM) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   com..."

     if(ilin.eq.0) then
        call calcavector(itri, com_field, avec)
        call eval_ops(avec, npoints, cot79)
     else
        cot79 = 0.
     end if
  endif


  ! B2I
  ! ~~~
  if(iand(fields, FIELD_B2I).eq.FIELD_B2I) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   B^-2..."

     b2i79 = 0.

     temp79a = ri2_79* &
          (pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2 + bzt79(:,OP_1)**2)

#if defined(USECOMPLEX) || defined(USE3D)
     temp79b = &
          (bft79(:,OP_DRP)**2 + bft79(:,OP_DZP)**2) &
          + 2.*ri_79* &
          (pst79(:,OP_DZ)*bft79(:,OP_DRP) - pst79(:,OP_DR)*bft79(:,OP_DZP))

     b2i79(1:npoints,OP_1 ) = 1./(temp79a(1:npoints) + temp79b(1:npoints))
#else
     b2i79(1:npoints,OP_1 ) = 1./temp79a(1:npoints)
#endif
     b2i79(:,OP_DR) = ri2_79 * &
          (pst79(:,OP_DR)*pst79(:,OP_DRR) + pst79(:,OP_DZ)*pst79(:,OP_DRZ) &
          +bzt79(:,OP_1 )*bzt79(:,OP_DR ))
     b2i79(:,OP_DZ) = ri2_79 * &
          (pst79(:,OP_DR)*pst79(:,OP_DRZ) + pst79(:,OP_DZ)*pst79(:,OP_DZZ) &
          +bzt79(:,OP_1 )*bzt79(:,OP_DZ ))

     if(itor.eq.1) then 
        b2i79(:,OP_DR) = b2i79(:,OP_DR) - ri_79*temp79a
     endif

#if defined(USECOMPLEX) || defined(USE3D)
     b2i79(:,OP_DR) = b2i79(:,OP_DR) + ri_79* &
          (pst79(:,OP_DZ )*bft79(:,OP_DRRP)-pst79(:,OP_DR )*bft79(:,OP_DRZP) &
          +pst79(:,OP_DRZ)*bft79(:,OP_DRP )-pst79(:,OP_DRR)*bft79(:,OP_DZP ))&
          +bft79(:,OP_DRRP)*bft79(:,OP_DRP)+bft79(:,OP_DRZP)*bft79(:,OP_DZP)
     b2i79(:,OP_DZ) = b2i79(:,OP_DZ) + ri_79* &
          (pst79(:,OP_DZ )*bft79(:,OP_DRZP)-pst79(:,OP_DR )*bft79(:,OP_DZZP) &
          +pst79(:,OP_DZZ)*bft79(:,OP_DRP )-pst79(:,OP_DRZ)*bft79(:,OP_DZP ))&
          +bft79(:,OP_DRZP)*bft79(:,OP_DRP)+bft79(:,OP_DZZP)*bft79(:,OP_DZP)

     if(itor.eq.1) then
        b2i79(:,OP_DR) = b2i79(:,OP_DR) - ri2_79* &
          (pst79(:,OP_DZ)*bft79(:,OP_DRP) - pst79(:,OP_DR)*bft79(:,OP_DZP))
     endif
#endif

     b2i79(:,OP_DR) = -2.*b2i79(:,OP_DR)*b2i79(:,OP_1)**2
     b2i79(:,OP_DZ) = -2.*b2i79(:,OP_DZ)*b2i79(:,OP_1)**2
  endif

  ! ETA
  ! ~~~
  if(iand(fields, FIELD_ETA).eq.FIELD_ETA) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   eta..."

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
     else if(iresfunc.eq.4) then
        eta79 = 0.
        temp79b = pet79(:,OP_1)/nt79(:,OP_1)
        where(real(temp79b).gt.0.) 
           temp79a = sqrt(temp79b)
           eta79(:,OP_1 ) = 1. / temp79a**3
           eta79(:,OP_DR) = (-3./2.) / temp79a**5 * &
                (pet79(:,OP_DR)/nt79(:,OP_1) &
                -pet79(:,OP_1)*nt79(:,OP_DR)/nt79(:,OP_1)**2)
           eta79(:,OP_DZ) = (-3./2.) / temp79a**5 * &
                (pet79(:,OP_DZ)/nt79(:,OP_1) &
                -pet79(:,OP_1)*nt79(:,OP_DZ)/nt79(:,OP_1)**2)
           eta79(:,OP_DRR) = (15./4.) / temp79a**7 * &
                (pet79(:,OP_DR)/nt79(:,OP_1) &
                -pet79(:,OP_1)*nt79(:,OP_DR)/nt79(:,OP_1)**2)**2 &
                + (-3./2.) / temp79a**5 * &
                (pet79(:,OP_DRR)/nt79(:,OP_1) &
                -2.*pet79(:,OP_DR)*nt79(:,OP_DR)/nt79(:,OP_1)**2 &
                -pet79(:,OP_1)*nt79(:,OP_DRR)/nt79(:,OP_1)**2 &
                +2.*pet79(:,OP_1)*nt79(:,OP_DR)**2/nt79(:,OP_1)**3)
           eta79(:,OP_DRZ) = (15./4.) / temp79a**7 * &
                (pet79(:,OP_DR)/nt79(:,OP_1) &
                -pet79(:,OP_1)*nt79(:,OP_DR)/nt79(:,OP_1)**2) &
                *(pet79(:,OP_DZ)/nt79(:,OP_1) &
                -pet79(:,OP_1)*nt79(:,OP_DZ)/nt79(:,OP_1)**2) &
                + (-3./2.) / temp79a**5 * &
                (pet79(:,OP_DRZ)/nt79(:,OP_1) &
                -pet79(:,OP_DR)*nt79(:,OP_DZ)/nt79(:,OP_1)**2 &
                -pet79(:,OP_DZ)*nt79(:,OP_DR)/nt79(:,OP_1)**2 &
                -pet79(:,OP_1)*nt79(:,OP_DRZ)/nt79(:,OP_1)**2 &
                +2.*pet79(:,OP_1)*nt79(:,OP_DR)*nt79(:,OP_DZ)/nt79(:,OP_1)**3)
           eta79(:,OP_DZZ) = (15./4.) / temp79a**7 * &
                (pet79(:,OP_DZ)/nt79(:,OP_1) &
                -pet79(:,OP_1)*nt79(:,OP_DZ)/nt79(:,OP_1)**2)**2 &
                + (-3./2.) / temp79a**5 * &
                (pet79(:,OP_DZZ)/nt79(:,OP_1) &
                -2.*pet79(:,OP_DZ)*nt79(:,OP_DZ)/nt79(:,OP_1)**2 &
                -pet79(:,OP_1)*nt79(:,OP_DZZ)/nt79(:,OP_1)**2 &
                +2.*pet79(:,OP_1)*nt79(:,OP_DZ)**2/nt79(:,OP_1)**3)
        end where

        eta79 = eta79 * 3.4e-22*n0_norm**2/(b0_norm**4*l0_norm)*17.
     else
        call calcavector(itri, resistivity_field, avec)
        call eval_ops(avec, npoints, eta79)
     end if
  end if

  ! KAP
  ! ~~~
  if(iand(fields, FIELD_KAP).eq.FIELD_KAP) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   kappa..."

     call calcavector(itri, kappa_field, avec)
     call eval_ops(avec, npoints, kap79)

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
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   sigma..."

     call calcavector(itri, sigma_field, avec)
     call eval_ops(avec, npoints, sig79)
  else
     sig79 = 0.
  end if

  ! MU
  ! ~~
  if(iand(fields, FIELD_MU).eq.FIELD_MU) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   vis..."

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
        call calcavector(itri, visc_field, avec)
        call eval_ops(avec, npoints, vis79)

        if(numvar.ge.3) then
           call calcavector(itri, visc_c_field, avec)
           call eval_ops(avec, npoints, vic79)
        endif
     endif

!     if(amupar.ne.0.) vip79 = amupar*pit79/2.
     if(amupar.ne.0.) vip79 = amupar
  end if

  if(gdef.eq.1) then
     call local_coeff_vector(itri, cl, .false.)
     do i=1, dofs_per_element
        avec = cl(i,:)
        call eval_ops(avec, npoints, mu79(:,:,i))
#ifdef USECOMPLEX
        mu79(:,OP_DP :OP_GSP, i) = mu79(:,OP_1:OP_GS,i)*rfac
        mu79(:,OP_DPP:OP_GSPP,i) = mu79(:,OP_1:OP_GS,i)*rfac**2
#endif
     end do
     nu79 = mu79

     if(equilibrate.eq.2) then
        if(itri.eq.1 .and. myrank.eq.0) then 
           print *, 1./int2(mu79(:,OP_1,i),mu79(:,OP_1,i)), &
                equil_fac(i,itri)
        end if
!        mu79(:,:,i) = mu79(:,:,i) / int2(mu79(:,OP_1,i),mu79(:,OP_1,i))
     endif

     if(equilibrate.ne.0) then 
        do i=1, dofs_per_element
           mu79(:,:,i) = mu79(:,:,i)*equil_fac(i,itri)
        end do
     end if
  endif

end subroutine define_fields

subroutine interpolate_size_field(itri)

  use basic
  use mesh_mod

  implicit none

  integer, intent(in) :: itri

  type(element_data) :: d
  double precision, dimension(3) :: node_sz
  real :: k,l,m,h

  if(ihypdx.eq.0) then
     sz79(:,OP_1  ) = 1.
     sz79(:,OP_DR ) = 0.
     sz79(:,OP_DZ ) = 0.
     sz79(:,OP_DRR) = 0.
     sz79(:,OP_DRZ) = 0.
     sz79(:,OP_DZZ) = 0.
     return
  end if

  call get_element_data(itri, d)
  h = d%b / (d%a + d%b)

#ifdef USESCOREC
  call getelmsizes(itri, node_sz)
#else
  node_sz = h
#endif
  ! use size**2 field
  node_sz = node_sz**ihypdx

  m = (node_sz(3) - node_sz(1) - h*(node_sz(2) - node_sz(1))) / d%c
  l = (node_sz(2) - node_sz(1)) / (d%a + d%b)
  k = node_sz(1) + l*d%b

  sz79(:,OP_1  ) = k + l*xi_79 + m*eta_79
  sz79(:,OP_DR ) = k + l*d%co + m*d%sn
  sz79(:,OP_DZ ) = k - l*d%sn + m*d%co
  sz79(:,OP_DRR) = 0.
  sz79(:,OP_DRZ) = 0.
  sz79(:,OP_DZZ) = 0. 

end subroutine interpolate_size_field

end module m3dc1_nint

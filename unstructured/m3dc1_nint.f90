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
  integer, parameter :: FIELD_XX1 =   512  ! UNUSED
  integer, parameter :: FIELD_XX2 =  1024  ! UNUSED
  integer, parameter :: FIELD_NI  =  2048
  integer, parameter :: FIELD_B2I =  4096
  integer, parameter :: FIELD_ETA =  8192
  integer, parameter :: FIELD_KAP = 16384
  integer, parameter :: FIELD_SIG = 32768
  integer, parameter :: FIELD_MU  = 65536
  integer, parameter :: FIELD_TE  =131072
  integer, parameter :: FIELD_TI  =262144
  integer, parameter :: FIELD_Q   =524288
  integer, parameter :: FIELD_F   =1048576
  integer, parameter :: FIELD_PF  =2097152
  integer, parameter :: FIELD_ES  =4194304
  integer, parameter :: FIELD_CD  =8388608
  integer, parameter :: FIELD_RAD  =16777216
  integer, parameter :: FIELD_KIN  =33554432
  integer, parameter :: FIELD_RE   =67108864 
  integer, parameter :: FIELD_WALL =134217728  ! 2^27
  integer, parameter :: FIELD_DENM =268435456  ! 2^28
  integer, parameter :: FIELD_JBS = 536870912  ! 2^29

! NOTE: All element-specific variables should be declared OMP THREADPRIVATE

  vectype, dimension(dofs_per_element, MAX_PTS, OP_NUM) :: mu79, nu79
!$OMP THREADPRIVATE(mu79,nu79)
!$acc declare create(mu79,nu79)
#ifdef USEST 
! logical basis functions must and nust
  vectype, dimension(dofs_per_element, MAX_PTS, OP_NUM) :: must79, nust79
!$OMP THREADPRIVATE(must79,nust79)
  vectype, dimension(MAX_PTS,OP_NUM) :: rst79, zst79 
!$OMP THREADPRIVATE(rst79,zst79)
! logical coords of quadrature points
  vectype, dimension(MAX_PTS) :: xl_79, zl_79 
!$OMP THREADPRIVATE(xl_79,zl_79)
#endif
  vectype, dimension(MAX_PTS) :: r_79, r2_79, r3_79, &
     ri_79, ri2_79, ri3_79, ri4_79, ri5_79, ri6_79, ri7_79, ri8_79
!$OMP THREADPRIVATE(r_79,r2_79,r3_79)
!$OMP THREADPRIVATE(ri_79,ri2_79,ri3_79,ri4_79,ri5_79,ri6_79,ri7_79,ri8_79)
  vectype, dimension(MAX_PTS) :: temp79a, temp79b, temp79c, &
       temp79d, temp79e, temp79f, temp79g
!$OMP THREADPRIVATE(temp79a,temp79b,temp79c,temp79d,temp79e,temp79f,temp79g)
  vectype, dimension(MAX_PTS, OP_NUM) :: tm79, ni79, b2i79, bi79
!$OMP THREADPRIVATE(tm79,ni79,b2i79,bi79)
  vectype, dimension(MAX_PTS, OP_NUM) :: ps179, bz179, pe179, n179, & 
       ph179, vz179, ch179, p179, ne179, pi179
!$OMP THREADPRIVATE(ps179,bz179,pe179,n179,ph179,vz179,ch179,p179,ne179,pi179)
  vectype, dimension(MAX_PTS, OP_NUM) :: pst79, bzt79, pet79, nt79, &
       pht79, vzt79, cht79, pt79, net79
!$OMP THREADPRIVATE(pst79,bzt79,pet79,nt79,pht79,vzt79,cht79,pt79,net79)
  vectype, dimension(MAX_PTS, OP_NUM) :: rho79, nw79
!$OMP THREADPRIVATE(rho79, nw79)
  vectype, dimension(MAX_PTS, OP_NUM) :: vis79, vic79, vip79, for79, es179
!$OMP THREADPRIVATE(vis79,vic79,vip79,for79,es179)
  vectype, dimension(MAX_PTS, OP_NUM) :: jt79, pit79, &
       eta79, etaRZ79,sig79, fy79, q79, cd79, totrad79, linerad79, bremrad79, ionrad79, reckrad79, recprad79, sie79, sii79, sir79
!$OMP THREADPRIVATE(jt79,pit79,eta79,etaRZ79,sig79,fy79,cd79,totrad79,linerad79,bremrad79,ionrad79,reckrad79,recprad79,sie79,sii79,sir79)
  vectype, dimension(MAX_PTS, OP_NUM) :: bfp079, bfp179, bfpt79
!$OMP THREADPRIVATE(bfp079,bfp179,bfpt79)
  vectype, dimension(MAX_PTS, OP_NUM) :: bf079, bf179, bft79
!$OMP THREADPRIVATE(bf079,bf179,bft79)
  vectype, dimension(MAX_PTS, OP_NUM) :: kap79, kar79, kax79
!$OMP THREADPRIVATE(kap79,kar79,kax79)
  vectype, dimension(MAX_PTS, OP_NUM) :: denm79
!$OMP THREADPRIVATE(denm79)
  vectype, dimension(MAX_PTS, OP_NUM) :: ps079, bz079, pe079, n079, &
       ph079, vz079, ch079, p079, ne079, pi079
!$OMP THREADPRIVATE(ps079,bz079,pe079,n079,ph079,vz079,ch079,p079,ne079,pi079)
  vectype, dimension(MAX_PTS, OP_NUM) :: pss79, bzs79
!$OMP THREADPRIVATE(pss79,bzs79)
  vectype, dimension(MAX_PTS, OP_NUM) :: bzx79, psx79, bfpx79, bfx79, psc79
!$OMP THREADPRIVATE(bzx79,psx79,bfpx79,bfx79,psc79)
  vectype, dimension(MAX_PTS, OP_NUM) :: pstx79, bztx79, bfptx79, bftx79
!$OMP THREADPRIVATE(pstx79,bztx79,bfptx79,bftx79)
  vectype, dimension(MAX_PTS, OP_NUM) :: te179, te079, tet79
!$OMP THREADPRIVATE(te179,te079,tet79)
  vectype, dimension(MAX_PTS, OP_NUM) :: ti179, ti079, tit79
!$OMP THREADPRIVATE(ti179,ti079,tit79)
  vectype, dimension(MAX_PTS, OP_NUM) :: q179, q079, qt79, qe179, qe079, qet79
!$OMP THREADPRIVATE(q179,q079,qt79,qe079,qet79)
#ifdef USEPARTICLES
  vectype, dimension(MAX_PTS, OP_NUM) :: pfpar79, pfper79, pf079
!$OMP THREADPRIVATE(pfpar79,pfper79,pf079)
  vectype, dimension(MAX_PTS, OP_NUM) :: pipar79, piper79, pfi079
!$OMP THREADPRIVATE(pipar79,piper79,pfi079)
  vectype, dimension(MAX_PTS, OP_NUM) :: vfpar79, vfpar079, vipar79
!$OMP THREADPRIVATE(vfpar79,vfpar079,vipar79)
  vectype, dimension(MAX_PTS, OP_NUM) :: nf79, nf079
!$OMP THREADPRIVATE(nf79, nf079)
  vectype, dimension(MAX_PTS, OP_NUM) :: nfi79, nfi079
!$OMP THREADPRIVATE(nfi79, nfi079)
  vectype, dimension(MAX_PTS, OP_NUM) :: rhof79
!$OMP THREADPRIVATE(rhof79)
#endif
  vectype, dimension(MAX_PTS, OP_NUM) :: nre079, nre179
!$OMP THREADPRIVATE(nre079,nre179)
  vectype, dimension(MAX_PTS, OP_NUM) :: wall79
!$OMP THREADPRIVATE(wall79)
  vectype, dimension(MAX_PTS) :: qd79
!$OMP THREADPRIVATE(qd79)
  vectype, dimension(MAX_PTS, OP_NUM) :: jbsl3179,jbsl3279,jbsl3479,jbsalpha79,jbsfluxavg_iBsq_B79,jbsfluxavg_G79,jbs_dtedpsit79
!$OMP THREADPRIVATE(jbsl3179,jbsl3279,jbsl3479,jbsalpha79,jbsfluxavg_iBsq_B79,jbsfluxavg_G79,jbs_dtedpsit79)
  vectype, dimension(MAX_PTS, OP_NUM) :: jbs_ftrap79,jbs_qR79,jbs_invAspectRatio79
!$OMP THREADPRIVATE(jbs_ftrap79,jbs_qR79,jbs_invAspectRatio79)

  ! precalculated terms
   real, private :: fterm(MAX_PTS, OP_NUM, coeffs_per_element)
!$OMP THREADPRIVATE(fterm)
  
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
       fterm(:,OP_1,p) = xpow(:,mi(p))*ypow(:,ni(p))
          
       if(mi(p).ge.1) then
          ! d_si terms
          temp = mi(p)*xpow(:,mi(p)-1) * ypow(:,ni(p))
          fterm(:,OP_DR,p) = fterm(:,OP_DR,p) + co*temp
          fterm(:,OP_DZ,p) = fterm(:,OP_DZ,p) + sn*temp           
             
          if(mi(p).ge.2) then
             ! d_si^2 terms
             temp = xpow(:,mi(p)-2)*(mi(p)-1)*mi(p) * ypow(:,ni(p))
             fterm(:,OP_DRR,p) = fterm(:,OP_DRR,p) + co2*temp
             fterm(:,OP_DZZ,p) = fterm(:,OP_DZZ,p) + sn2*temp
             fterm(:,OP_DRZ,p) = fterm(:,OP_DRZ,p) + cosn*temp
             fterm(:,OP_LP,p) = fterm(:,OP_LP,p) + temp
          endif
       endif
       if(ni(p).ge.1) then
          ! d_eta terms
          temp = xpow(:,mi(p)) * ypow(:,ni(p)-1)*ni(p)
          fterm(:,OP_DR,p) = fterm(:,OP_DR,p) - sn*temp
          fterm(:,OP_DZ,p) = fterm(:,OP_DZ,p) + co*temp
          
          if(ni(p).ge.2) then
             ! d_eta^2 terms
             temp = xpow(:,mi(p)) * ypow(:,ni(p)-2)*(ni(p)-1)*ni(p)
             fterm(:,OP_DRR,p) = fterm(:,OP_DRR,p) + sn2*temp
             fterm(:,OP_DZZ,p) = fterm(:,OP_DZZ,p) + co2*temp
             fterm(:,OP_DRZ,p) = fterm(:,OP_DRZ,p) - cosn*temp
             fterm(:,OP_LP,p) = fterm(:,OP_LP,p) + temp
          endif
          
          if(mi(p).ge.1) then
             ! d_eta_si terms
             temp = xpow(:,mi(p)-1)*mi(p) * ypow(:,ni(p)-1)*ni(p)
             
             fterm(:,OP_DRR,p) = fterm(:,OP_DRR,p) - 2.*cosn*temp
             fterm(:,OP_DZZ,p) = fterm(:,OP_DZZ,p) + 2.*cosn*temp
             fterm(:,OP_DRZ,p) = fterm(:,OP_DRZ,p) + (co2-sn2)*temp
          endif
       endif
       
       ! for surface terms, higher derivatives may be taken
       if(surface_int) then
          if(mi(p).ge.2) then
             if(ni(p).ge.1) then
                ! d_si^2 d_eta terms
                temp = xpow(:,mi(p)-2)*ypow(:,ni(p)-1)*(mi(p)-1)*mi(p)*ni(p)
                fterm(:,OP_LPR,p) = fterm(:,OP_LPR,p) - sn*temp
                fterm(:,OP_LPZ,p) = fterm(:,OP_LPZ,p) + co*temp
             endif
          endif
          if(ni(p).ge.2) then
             if(mi(p).ge.1) then
                ! d_eta^2 d_si terms
                temp = xpow(:,mi(p)-1)*ypow(:,ni(p)-2)*mi(p)*(ni(p)-1)*ni(p)
                fterm(:,OP_LPR,p) = fterm(:,OP_LPR,p) + co*temp
                fterm(:,OP_LPZ,p) = fterm(:,OP_LPZ,p) + sn*temp
             endif
          endif
          
          if(mi(p).ge.3) then
             ! d_si^3 terms
             temp = xpow(:,mi(p)-3)*ypow(:,ni(p))*(mi(p)-2)*(mi(p)-1)*mi(p)
             fterm(:,OP_LPR,p) = fterm(:,OP_LPR,p) + co*temp
             fterm(:,OP_LPZ,p) = fterm(:,OP_LPZ,p) + sn*temp
          endif
          if(ni(p).ge.3) then
             ! d_eta^3 terms
             temp = xpow(:,mi(p))*ypow(:,ni(p)-3)*(ni(p)-2)*(ni(p)-1)*ni(p)
             fterm(:,OP_LPR,p) = fterm(:,OP_LPR,p) - sn*temp
             fterm(:,OP_LPZ,p) = fterm(:,OP_LPZ,p) + co*temp
          endif
       endif
       
       ! Grad-Shafranov operator, and
       ! cylindrical correction to Laplacian
       fterm(:,OP_GS,p) = fterm(:,OP_LP,p)
       if(itor.eq.1) then
          fterm(:,OP_GS,p) = fterm(:,OP_GS,p) - fterm(:,OP_DR,p)*ri(:)
          fterm(:,OP_LP,p) = fterm(:,OP_LP,p) + fterm(:,OP_DR,p)*ri(:)
          
          if(surface_int) then
             fterm(:,OP_LPR,p) = fterm(:,OP_LPR,p) + fterm(:,OP_DRR,p)*ri(:) &
                  - fterm(:,OP_DR,p)*ri(:)*ri(:)
             fterm(:,OP_LPZ,p) = fterm(:,OP_LPZ,p) + fterm(:,OP_DRZ,p)*ri(:)
          endif
       endif      

#ifdef USE3D
       do op=1, OP_NUM_POL
          temp(:) = fterm(:,op,p)

          do i=1, coeffs_per_dphi
             j = p + (i-1)*coeffs_per_tri

             fterm(:,op,j) = temp(:)*zpow(:,li(i))

             ! first toroidal derivative
             if(li(i).ge.1) then
                fterm(:,op+OP_NUM_POL,j) = temp(:) &
                     *zpow(:,li(i)-1)*li(i)
             endif
             ! second toroidal derivative
             if(li(i).ge.2) then
                fterm(:,op+2*OP_NUM_POL,j) = temp(:) &
                     *zpow(:,li(i)-2)*(li(i)-1)*li(i)
             endif
          end do
       end do
#endif
    end do

  end subroutine precalculate_terms

  subroutine define_basis(itri)
    use basic
    implicit none

    integer, intent(in) :: itri
    real, dimension(dofs_per_element,coeffs_per_element) :: cl
    integer :: i
#ifndef USEBLAS
    integer :: p, op
#elif defined(USECOMPLEX)
    real, dimension(dofs_per_element, MAX_PTS, OP_NUM) :: real_mu
#endif

    if(iprecompute_metric.eq.1) then
       cl = ctri(:,:,itri)
    else
       call local_coeff_vector(itri, cl)
    endif 

#ifdef USEBLAS

#ifdef USECOMPLEX
    call dgemm('N','T',dofs_per_element,MAX_PTS*OP_NUM,coeffs_per_element, &
               1.,cl,dofs_per_element,fterm,MAX_PTS*OP_NUM, &
               0.,real_mu,dofs_per_element)
    mu79 = real_mu
#else
    call dgemm('N','T',dofs_per_element,MAX_PTS*OP_NUM,coeffs_per_element, &
               1.,cl,dofs_per_element,fterm,MAX_PTS*OP_NUM, &
               0.,mu79,dofs_per_element)
#endif

#else ! USEBLAS not defined
    mu79 = 0.
    do op=1, OP_NUM
       do i=1, dofs_per_element
          do p=1, coeffs_per_element
             mu79(i,:,op) = mu79(i,:,op) + fterm(:,op,p)*cl(i,p)
          end do
       end do
    end do
#endif

#ifdef USECOMPLEX
    mu79(:,:,OP_DP :OP_GSP ) = mu79(:,:,OP_1:OP_GS)*rfac
    mu79(:,:,OP_DPP:OP_GSPP) = mu79(:,:,OP_1:OP_GS)*rfac**2
#endif
    nu79 = mu79

    if(equilibrate.ne.0) then
       do i=1, dofs_per_element
          mu79(i,:,:) = mu79(i,:,:)*equil_fac(i, itri)
       end do
    end if

  end subroutine define_basis

#ifdef USEST 
  ! calculate physical basis functions from logical basis functions
  subroutine define_physical_basis(itri)
    use basic
    use arrays
    implicit none

    integer, intent(in) :: itri
    integer :: i
    vectype, dimension(MAX_PTS) :: di_79, di2_79
    real :: negone
    real :: value_nan ! NaN to be assigned to forbidden operations
    
    ! calculate logical derivatives of geometry
    call eval_ops(itri, rst, rst79)
    call eval_ops(itri, zst, zst79)

    ! save logical coordinates
    xl_79 = x_79
    zl_79 = z_79
    ! update physical element data
    x_79 = rst79(:,OP_1)
    z_79 = zst79(:,OP_1)
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

    ! calculate expressions needed
    ! inverse of D = Rx*Zy - Ry*Zx
    di_79 = 1./(rst79(:,OP_DR)*zst79(:,OP_DZ) - zst79(:,OP_DR)*rst79(:,OP_DZ))
    di2_79 = di_79*di_79
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, di_79(1) 
    ! Dx = Rx*Zxy + Rxx*Zy - Ry*Zxx - Rxy*Zx  
    temp79a = rst79(:,OP_DR)*zst79(:,OP_DRZ) + rst79(:,OP_DRR)*zst79(:,OP_DZ)&
            - rst79(:,OP_DZ)*zst79(:,OP_DRR) - rst79(:,OP_DRZ)*zst79(:,OP_DR)
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, temp79a(1) 
    ! Dy = Rx*Zyy + Rxy*Zy - Ry*Zxy - Ryy*Zx 
    temp79b = rst79(:,OP_DR)*zst79(:,OP_DZZ) + rst79(:,OP_DRZ)*zst79(:,OP_DZ)&
            - rst79(:,OP_DZ)*zst79(:,OP_DRZ) - rst79(:,OP_DZZ)*zst79(:,OP_DR)
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, temp79b(1) 
    ! F = Rx*Dy - Ry*Dx
    temp79c = rst79(:,OP_DR)*temp79b - rst79(:,OP_DZ)*temp79a
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, temp79c(1) 
    ! G = Zy*Dx - Zx*Dy
    temp79d = zst79(:,OP_DZ)*temp79a - zst79(:,OP_DR)*temp79b
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, temp79d(1) 
#ifdef USE3D
    ! Dz = Rx*Zyz + Rxz*Zy - Ry*Zxz - Ryz*Zx  
    temp79e = rst79(:,OP_DR)*zst79(:,OP_DZP) + rst79(:,OP_DRP)*zst79(:,OP_DZ)&
            - rst79(:,OP_DZ)*zst79(:,OP_DRP) - rst79(:,OP_DZP)*zst79(:,OP_DR)
#endif    
    do i=1, dofs_per_element
      ! fR = (Zy/D)*fx - (Zx/D)*fy
      mu79(i,:,OP_DR) = di_79*zst79(:,OP_DZ)*must79(i,:,OP_DR)&
                      - di_79*zst79(:,OP_DR)*must79(i,:,OP_DZ) 
      ! fZ = (Rx/D)*fy - (Ry/D)*fx
      mu79(i,:,OP_DZ) = di_79*rst79(:,OP_DR)*must79(i,:,OP_DZ)&
                      - di_79*rst79(:,OP_DZ)*must79(i,:,OP_DR) 
      ! fRR = (Zy/D)^2*fxx + (Zx/D)^2*fyy - 2(Zx*Zy/D^2)*fxy
      !     + [(Zy*Zxy - Zx*Zyy)/D^2 ]*fx      
      !     + [(Zx*Zxy - Zy*Zxx)/D^2 ]*fy - (G/D^2)*fR      
      mu79(i,:,OP_DRR) = (zst79(:,OP_DZ)**2*must79(i,:,OP_DRR)& 
                       + zst79(:,OP_DR)**2*must79(i,:,OP_DZZ)& 
                       - 2*zst79(:,OP_DR)*zst79(:,OP_DZ)*must79(i,:,OP_DRZ)& 
                       + (zst79(:,OP_DZ)*zst79(:,OP_DRZ) - zst79(:,OP_DR)*zst79(:,OP_DZZ))&
                       * must79(i,:,OP_DR) + must79(i,:,OP_DZ) &  
                       * (zst79(:,OP_DR)*zst79(:,OP_DRZ) - zst79(:,OP_DZ)*zst79(:,OP_DRR))&
                       - temp79d*mu79(i,:,OP_DR))*di2_79
      ! fZZ = (Ry/D)^2*fxx + (Rx/D)^2*fyy - 2(Rx*Ry/D^2)*fxy
      !     + [(Ry*Rxy - Rx*Ryy)/D^2 ]*fx      
      !     + [(Rx*Rxy - Ry*Rxx)/D^2 ]*fy -(F/D^2)*fZ 
      mu79(i,:,OP_DZZ) = (rst79(:,OP_DZ)**2*must79(i,:,OP_DRR)& 
                       + rst79(:,OP_DR)**2*must79(i,:,OP_DZZ)& 
                       - 2*rst79(:,OP_DR)*rst79(:,OP_DZ)*must79(i,:,OP_DRZ)& 
                       + (rst79(:,OP_DZ)*rst79(:,OP_DRZ) - rst79(:,OP_DR)*rst79(:,OP_DZZ))&
                       * must79(i,:,OP_DR) + must79(i,:,OP_DZ) &
                       * (rst79(:,OP_DR)*rst79(:,OP_DRZ) - rst79(:,OP_DZ)*rst79(:,OP_DRR))&
                       - temp79c*mu79(i,:,OP_DZ))*di2_79
      ! fRZ = [(Rx*Zy + Ry*Zx)/D^2]*fxy - (Ry*Zy/D^2)*fxx - (Rx*Zx/D^2)*fyy 
      !     - [(Zy*Rxy - Zx*Ryy)/D^2 ]*fx      
      !     - [(Zx*Rxy - Zy*Rxx)/D^2 ]*fy - (G/D^2)*fZ     
      mu79(i,:,OP_DRZ) = ((rst79(:,OP_DR)*zst79(:,OP_DZ)&
                       + zst79(:,OP_DR)*rst79(:,OP_DZ))*must79(i,:,OP_DRZ)& 
                       - rst79(:,OP_DR)*zst79(:,OP_DR)*must79(i,:,OP_DZZ)& 
                       - rst79(:,OP_DZ)*zst79(:,OP_DZ)*must79(i,:,OP_DRR)& 
                       - (zst79(:,OP_DZ)*rst79(:,OP_DRZ) - zst79(:,OP_DR)*rst79(:,OP_DZZ))&
                       * must79(i,:,OP_DR) - must79(i,:,OP_DZ) &
                       * (zst79(:,OP_DR)*rst79(:,OP_DRZ) - zst79(:,OP_DZ)*rst79(:,OP_DRR))&
                       - temp79d*mu79(i,:,OP_DZ))*di2_79
      mu79(i,:,OP_LP) = mu79(i,:,OP_DRR) + mu79(i,:,OP_DZZ) 
      mu79(i,:,OP_GS) = mu79(i,:,OP_LP) 
      if(itor.eq.1) then ! toroidal corrections
        mu79(i,:,OP_LP) = mu79(i,:,OP_LP) + mu79(i,:,OP_DR)*ri_79 
        mu79(i,:,OP_GS) = mu79(i,:,OP_GS) - mu79(i,:,OP_DR)*ri_79 
      end if
#ifdef USE3D
      ! fP = - Rz*fR - Zz*fZ + fz
      mu79(i,:,OP_DP) = - rst79(:,OP_DP)*mu79(i,:,OP_DR)&
                        - zst79(:,OP_DP)*mu79(i,:,OP_DZ)& 
                        + must79(i,:,OP_DP)
      ! fRz = (Zy/D)*fxz - (Zx/D)*fyz - (Dz/D)*fR + (Zyz/D)*fx - (Zxz/D)*fy
      mu79(i,:,OP_DRP) = (zst79(:,OP_DZ)*must79(i,:,OP_DRP)&
                       - zst79(:,OP_DR)*must79(i,:,OP_DZP)& 
                       + zst79(:,OP_DZP)*must79(i,:,OP_DR)&
                       - zst79(:,OP_DRP)*must79(i,:,OP_DZ)& 
                       - temp79e*mu79(i,:,OP_DR))*di_79 
      ! fZz = (Rx/D)*fyz - (Ry/D)*fxz - (Dz/D)*fZ + (Rxz/D)*fy - (Ryz/D)*fx
      mu79(i,:,OP_DZP) = (rst79(:,OP_DR)*must79(i,:,OP_DZP)&
                       - rst79(:,OP_DZ)*must79(i,:,OP_DRP)& 
                       + rst79(:,OP_DRP)*must79(i,:,OP_DZ)& 
                       - rst79(:,OP_DZP)*must79(i,:,OP_DR)&
                       - temp79e*mu79(i,:,OP_DZ))*di_79 
      ! fPz = - Rz*fRz - Zz*fZz + fzz - Rzz*fR -Zzz*fZ
      mu79(i,:,OP_DPP) = - rst79(:,OP_DP)*mu79(i,:,OP_DRP)&
                         - zst79(:,OP_DP)*mu79(i,:,OP_DZP)& 
                         - rst79(:,OP_DPP)*mu79(i,:,OP_DR)&
                         - zst79(:,OP_DPP)*mu79(i,:,OP_DZ)& 
                         + must79(i,:,OP_DPP)
      ! fRP = - Rz*fRR - Zz*fZR + fRz
      mu79(i,:,OP_DRP) = - rst79(:,OP_DP)*mu79(i,:,OP_DRR)&
                         - zst79(:,OP_DP)*mu79(i,:,OP_DRZ)& 
                         + mu79(i,:,OP_DRP)
      ! fZP = - Rz*fRZ - Zz*fZZ + fZz
      mu79(i,:,OP_DZP) = - rst79(:,OP_DP)*mu79(i,:,OP_DRZ)&
                         - zst79(:,OP_DP)*mu79(i,:,OP_DZZ)& 
                         + mu79(i,:,OP_DZP)
      ! fPP = - Rz*fRP - Zz*fZP + fPz
      mu79(i,:,OP_DPP) = - rst79(:,OP_DP)*mu79(i,:,OP_DRP)&
                         - zst79(:,OP_DP)*mu79(i,:,OP_DZP)& 
                         + mu79(i,:,OP_DPP)
#endif
    end do
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, mu79(1,1,:) 
    !if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, must79(1,1,:) 
#ifdef USE3D
    ! set forbidden operations to NaN
    negone = -1.
    value_nan = sqrt(negone)
    mu79(:,:,OP_DRRP:OP_GSP) = value_nan
    mu79(:,:,OP_DRPP:OP_GSPP) = value_nan 
    !mu79(:,:,OP_GSP) = value_nan 
#endif
    nu79 = mu79
    ! modify Jabobian
    if(ijacobian.eq.1) weight_79 = weight_79/di_79 
  end subroutine define_physical_basis
#endif

  !===============================================
  ! eval_ops
  ! --------
  !
  ! evaluates linear unitary operators
  !===============================================
  subroutine eval_ops(itri,fin,outarr,rfac)
    use field
    implicit none

    integer, intent(in) :: itri
    type(field_type), intent(in) :: fin
    vectype, dimension(MAX_PTS, OP_NUM), intent(out) :: outarr
    complex, optional :: rfac

    vectype, dimension(dofs_per_element) :: dofs

#ifndef USEBLAS
    integer :: i, op
#endif

    call get_element_dofs(fin, itri, dofs)

    ! Calculate outarr = dot(nu79, dofs)
#ifdef USEBLAS

#ifdef USECOMPLEX
    call zgemv('T',dofs_per_element,MAX_PTS*OP_NUM,&
         (1.,0.),nu79,dofs_per_element,dofs,1,(0.,0.),outarr,1)
    outarr(:,OP_NUM_POL+1:) = 0.
#else
    call dgemv('T',dofs_per_element,MAX_PTS*OP_NUM,&
         1.,nu79,dofs_per_element,dofs,1,0.,outarr,1)
#endif

#else ! USEBLAS not defined
    outarr = 0.
#ifdef USECOMPLEX
    do op=1, OP_NUM_POL
       do i=1, dofs_per_element
          outarr(:,op) = outarr(:,op) + dofs(i)*nu79(i,:,op)
       end do
    end do
#else
    do op=1, OP_NUM
       do i=1, dofs_per_element
          outarr(:,op) = outarr(:,op) + dofs(i)*nu79(i,:,op)
       end do
    end do
#endif

#endif ! on USEBLAS


#ifdef USECOMPLEX
    if(present(rfac)) then
       outarr(:,OP_DP :OP_GSP ) = outarr(:,OP_1 :OP_GS)*rfac
       outarr(:,OP_DPP:OP_GSPP) = outarr(:,OP_1 :OP_GS)*rfac**2
    end if
#endif

  end subroutine eval_ops

  !=====================================================
  ! define_fields
  !=====================================================
  subroutine define_fields(itri, fieldi, gdef, ilin, ieqs)
    use basic
    use mesh_mod
    use arrays
    use math
    use resistive_wall

    implicit none
  
    integer, intent(in) :: itri, fieldi, gdef, ilin
    integer, intent(in), optional :: ieqs

    real :: fac
    real :: p_floor
    integer :: izone, ieqsub, fields, i, iz
    integer, dimension(MAX_PTS) :: izarr
    type(element_data) :: d
    real :: kr_tmin, kr_tmax
    integer::ipoint

    fields = fieldi

    if(present(ieqs)) then
       ieqsub = ieqs
    else
       ieqsub = eqsubtract
    end if

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

    call get_zone(itri, izone)

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
    
    call precalculate_terms(xi_79,zi_79,eta_79,d%co,d%sn,ri_79,npoints)
    call define_basis(itri)

#ifdef USEST 
    if (igeometry.eq.1) then
       ! copy logical basis functions
       must79 = mu79
       nust79 = mu79
       if(ilog.eq.0) then  ! use physical basis when ilog==0 
          call define_physical_basis(itri)
       else
          if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, &
          "Use logical basis..."
       end if
    end if
#endif

    if(ijacobian.eq.1) weight_79 = weight_79 * r_79

    ! some field calculations require other field calculation first
    if(iand(fields, FIELD_N).eq.FIELD_N) then
       if(idenfunc.eq.3) fields = ior(fields,FIELD_PSI)
    end if
    if(iand(fields, FIELD_ETA).eq.FIELD_ETA) then
       if(iresfunc.eq.2 .or. iresfunc.eq.3) fields = ior(fields,FIELD_PSI)
       if(iresfunc.eq.4) then
          fields = ior(ior(fields,FIELD_N),FIELD_P)
          if(itemp.eq.1) fields = ior(fields,FIELD_TE)
       end if
    end if
    if(iand(fields, FIELD_KAP).eq.FIELD_KAP) then
       if(ikappafunc.eq.5) then
          fields = ior(ior(ior(fields,FIELD_N),FIELD_P),FIELD_PSI)
          if(itemp.eq.1) fields = ior(fields,FIELD_TE)
       end if
       if(ikapparfunc.eq.1) fields = ior(fields,FIELD_TE)
       if(ikapparfunc.eq.2) then
          fields =ior(fields,FIELD_TE)
          !fields = ior(ior(fields,FIELD_N),FIELD_P)
          !if(itemp.eq.1) then
          !   fields =ior(fields,FIELD_TE)
          !   if(ipres.eq.1) fields = ior(fields,FIELD_TI)
          !end if
       end if
    end if
    if(iand(fields, FIELD_DENM).eq.FIELD_DENM) then
       if(idenmfunc.eq.1) then
          fields = ior(ior(ior(fields,FIELD_N),FIELD_P),FIELD_PSI)
          if(itemp.eq.1) fields = ior(fields,FIELD_TE)
       end if
    end if
    if(iand(fields, FIELD_MU).eq.FIELD_MU) then
       if(ivisfunc.eq.3) fields = ior(fields,FIELD_PSI)
       if(ivisfunc.eq.4) fields = ior(fields,FIELD_WALL)
    end if

    ! PHI
    ! ~~~
    if(iand(fields, FIELD_PHI).eq.FIELD_PHI) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   U..."
       
       if(ilin.eq.0) then
          call eval_ops(itri, u_field(1), ph179, rfac)
       else
          ph179 = 0.
       endif
       
       if(ieqsub.eq.1) then
          call eval_ops(itri, u_field(0), ph079)
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
          call eval_ops(itri, psi_ext, psx79, rfac)
       else
          psx79 = 0.
       end if
       
       if(ilin.eq.0) then
          call eval_ops(itri, psi_field(1), ps179, rfac)
       else
          ps179 = 0.
       end if
       
       if(ieqsub.eq.1) then
          call eval_ops(itri, psi_field(0), ps079)
       else
          ps079 = 0.
       endif

       if(icsubtract.eq.1) then 
          call eval_ops(itri, psi_coil_field, psc79)
          ps079 = ps079 + psc79
       else
          psc79 = 0.
       end if

       pst79 = ps079 + ps179
       pss79 = ps079 + ps179/2.
       pstx79 = pst79 + psx79
    endif

    ! V
    ! ~
    if(iand(fields, FIELD_V).eq.FIELD_V) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   V..."
       
       if(ilin.eq.0) then
          call eval_ops(itri, vz_field(1), vz179, rfac)
       else
          vz179 = 0.
       end if
       
       if(ieqsub.eq.1) then
          call eval_ops(itri, vz_field(0), vz079)
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
          call eval_ops(itri, bz_ext, bzx79, rfac)
#if defined(USECOMPLEX) || defined(USE3D)    
          call eval_ops(itri, bf_ext, bfx79, rfac)
          call eval_ops(itri, bfp_ext, bfpx79, rfac)
#endif
       else
          bzx79 = 0.
          bfx79 = 0.
          bfpx79 = 0.
       endif
 
       if(ilin.eq.0) then
          call eval_ops(itri, bz_field(1), bz179, rfac)
#if defined(USECOMPLEX) || defined(USE3D)    
          call eval_ops(itri, bf_field(1), bf179, rfac)
          call eval_ops(itri, bfp_field(1), bfp179, rfac)
#endif
       else
          bz179 = 0.
          bf179 = 0.
          bfp179 = 0.
       endif
       
       if(ieqsub.eq.1) then
          call eval_ops(itri, bz_field(0), bz079)
          bzt79 = bz079 + bz179
          bzs79 = bz079 + bz179/2.
          
#if defined(USECOMPLEX) || defined(USE3D)
          call eval_ops(itri, bf_field(0), bf079)
          bft79 = bf079 + bf179
          call eval_ops(itri, bfp_field(0), bfp079)
          bfpt79 = bfp079 + bfp179
#endif
       else
          bz079 = 0.
          bzt79 = bz179
          bzs79 = bz179/2.
          
#if defined(USECOMPLEX) || defined(USE3D)
          bf079 = 0.
          bft79 = bf179
          bfp079 = 0.
          bfpt79 = bfp179
#endif
       endif

       bztx79 = bzt79 + bzx79
#if defined(USECOMPLEX) || defined(USE3D)
       bftx79 = bft79 + bfx79
       bfptx79 = bfpt79 + bfpx79
#endif

       if(numvar.eq.1) bzs79 = bzt79
    endif

    ! CHI
    ! ~~~
    if(iand(fields, FIELD_CHI).eq.FIELD_CHI) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   chi..."
       
       if(ilin.eq.0) then
          call eval_ops(itri, chi_field(1), ch179, rfac)
       else
          ch179 = 0.
       end if
       
       if(ieqsub.eq.1) then
          call eval_ops(itri, chi_field(0), ch079)
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
          call eval_ops(itri, p_field(1), p179, rfac)
          call eval_ops(itri, pe_field(1), pe179, rfac)
       else
          p179 = 0.
          pe179 = 0.
       end if
       
       if(ieqsub.eq.1) then
          call eval_ops(itri, p_field(0), p079)
          call eval_ops(itri, pe_field(0), pe079)

          pet79 = pe079 + pe179
          pt79  =  p079 +  p179
       else
          pe079 = 0.
          p079 = 0.
          pet79 = pe179
          pt79  =  p179
       endif
       pi179 = p179 - pe179
       pi079 = p079 - pe079
       pit79 = pt79 - pet79

       p_floor = iset_pe_floor*pe_floor + iset_pi_floor*pi_floor

       if(iset_pe_floor.eq.1) then
          if(ilin.eq.0) then
             where(real(pet79(:,OP_1)).lt.pe_floor)
                pe179(:,OP_1) = pe_floor - pe079(:,OP_1)
             end where
          end if
          where(real(pet79(:,OP_1)).lt.pe_floor)
             pet79(:,OP_1) = pe_floor
          end where
       end if

       if(iset_pi_floor.eq.1) then
          if(ilin.eq.0) then
             where(real(pit79(:,OP_1)).lt.pi_floor)
                pi179(:,OP_1) = pi_floor - pi079(:,OP_1)
             end where
          end if
          where(real(pit79(:,OP_1)).lt.pi_floor)
             pit79(:,OP_1) = pi_floor
          end where
       end if

       if(iset_pe_floor.eq.1 .or. iset_pi_floor.eq.1) then
          if(ilin.eq.0) then
             where(real(pt79(:,OP_1)).lt.p_floor)
                p179(:,OP_1) = p_floor - p079(:,OP_1)
             end where
          end if
          where(real(pt79(:,OP_1)).lt.p_floor)
             pt79(:,OP_1) = p_floor
          end where
       end if
       
    endif
   
    
    ! N
    ! ~
    if(iand(fields, FIELD_N).eq.FIELD_N) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   n..."
       
       if(ilin.eq.0) then
          call eval_ops(itri, den_field(1), n179, rfac)
          call eval_ops(itri, ne_field(1), ne179, rfac)
       else
          n179 = 0.
          ne179 = 0.
       end if

       if(ieqsub.eq.1) then
          if(idenfunc.eq.3) then
             n079 = 0.
             if(izone.gt.1) then
                n079(:,OP_1) = den_edge
             else 
                temp79a = (pst79(:,OP_1) - psimin)/(psibound - psimin)
                temp79b = (pst79(:,OP_DR)*(x_79 - xmag) &
                     +     pst79(:,OP_DZ)*(z_79 - zmag))*(psibound-psimin)
                where(real(temp79a).lt.denoff .and. real(temp79b).gt.0.)
                   n079(:,OP_1) = den0
                elsewhere
                   n079(:,OP_1) = den_edge
                end where
             end if
             ne079 = n079*z_ion
          else
             call eval_ops(itri, den_field(0), n079)
             call eval_ops(itri, ne_field(0), ne079)
          end if
          nt79 = n079 + n179
          net79 = ne079 + ne179
       else
          n079 = 0.
          nt79 = n179
          ne079 = 0.
          net79 = ne179
       endif

     if(iset_ne_floor.eq.1) then
        if(ilin.eq.0) then
           where(real(net79(:,OP_1)).lt.ne_floor)
              ne179(:,OP_1) = ne_floor - ne079(:,OP_1)
           end where
        end if
        where(real(net79(:,OP_1)).lt.ne_floor)
           net79(:,OP_1) = ne_floor
        end where
     end if
     
     if(iset_ni_floor.eq.1) then
        if(ilin.eq.0) then
           where(real(nt79(:,OP_1)).lt.ni_floor)
              n179(:,OP_1) = ni_floor - n079(:,OP_1)
           end where
        end if
        where(real(nt79(:,OP_1)).lt.ni_floor)
           nt79(:,OP_1) = ni_floor
        end where
     end if

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
#if defined(USECOMPLEX) || defined(USE3D)
     ni79(:,OP_DP) = -nt79(:,OP_DP)*ni79(:,OP_1)**2
     ni79(:,OP_DRP) = -nt79(:,OP_DRP)*ni79(:,OP_1)**2
     ni79(:,OP_DZP) = -nt79(:,OP_DZP)*ni79(:,OP_1)**2
     ni79(:,OP_DRRP) = -nt79(:,OP_DRRP)*ni79(:,OP_1)**2
     ni79(:,OP_DRZP) = -nt79(:,OP_DRZP)*ni79(:,OP_1)**2
     ni79(:,OP_DZZP) = -nt79(:,OP_DZZP)*ni79(:,OP_1)**2
     ni79(:,OP_LPP) = -nt79(:,OP_LPP)*ni79(:,OP_1)**2
     ni79(:,OP_GSP) = -nt79(:,OP_GSP)*ni79(:,OP_1)**2

     ni79(:,OP_DPP) = -nt79(:,OP_DPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_DP)*ni79(:,OP_1)*ni79(:,OP_DP)
     ni79(:,OP_DRPP) = -nt79(:,OP_DRPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_DRP)*ni79(:,OP_1)*ni79(:,OP_DP)
     ni79(:,OP_DZPP) = -nt79(:,OP_DZPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_DZP)*ni79(:,OP_1)*ni79(:,OP_DP)
     ni79(:,OP_DRRPP) = -nt79(:,OP_DRRPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_DRRP)*ni79(:,OP_1)*ni79(:,OP_DP)
     ni79(:,OP_DRZPP) = -nt79(:,OP_DRZPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_DRZP)*ni79(:,OP_1)*ni79(:,OP_DP)
     ni79(:,OP_DZZPP) = -nt79(:,OP_DZZPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_DZZP)*ni79(:,OP_1)*ni79(:,OP_DP)
     ni79(:,OP_LPPP) = -nt79(:,OP_LPPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_LPP)*ni79(:,OP_1)*ni79(:,OP_DP)
     ni79(:,OP_GSPP) = -nt79(:,OP_GSPP)*ni79(:,OP_1)**2 &
          - 2.*nt79(:,OP_GSP)*ni79(:,OP_1)*ni79(:,OP_DP)
#endif


     if(linear.eq.0) then
        where(ni79.ne.ni79) ni79 = 0.
        where(real(ni79(:,OP_1)).lt.0.) ni79(:,OP_1) = 0.
     end if
  endif

  ! TE
  ! ~~~
  if(iand(fields, FIELD_TE).eq.FIELD_TE) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   TE..."
     
     if(ilin.eq.0) then
        call eval_ops(itri, te_field(1), te179, rfac)
     else
        te179 = 0.
     endif
     
     if(ieqsub.eq.1) then
        call eval_ops(itri, te_field(0), te079)
        tet79 = te079 + te179
     else
        te079 = 0.
        tet79 = te179
     endif

     if(iset_te_floor.eq.1) then
        if(ilin.eq.0) then
           where(real(tet79(:,OP_1)).lt.te_floor)
              te179(:,OP_1) = te_floor - te079(:,OP_1)
           end where
        end if
        where(real(tet79(:,OP_1)).lt.te_floor)
           tet79(:,OP_1) = te_floor
        end where
     end if

  endif

  ! TI
  ! ~~~
  if(iand(fields, FIELD_TI).eq.FIELD_TI) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   TI..."
     
     if(ilin.eq.0) then
        call eval_ops(itri, ti_field(1), ti179, rfac)
     else
        ti179 = 0.
     endif
     
     if(ieqsub.eq.1) then
        call eval_ops(itri, ti_field(0), ti079)
        tit79 = ti079 + ti179
     else
        ti079 = 0.
        tit79 = ti179
     endif

     if(iset_ti_floor.eq.1) then
        if(ilin.eq.0) then
           where(real(tit79(:,OP_1)).lt.ti_floor)
              ti179(:,OP_1) = ti_floor - ti079(:,OP_1)
           end where
        end if
        where(real(tit79(:,OP_1)).lt.ti_floor)
           tit79(:,OP_1) = ti_floor
        end where
     end if

  endif
  
  ! J
  ! ~
  if(iand(fields, FIELD_J).eq.FIELD_J) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   j..."

     if(ilin.eq.0) then
        call eval_ops(itri, jphi_field, jt79)
     else
        jt79 = 0.
     end if
  endif


  ! B2I
  ! ~~~
  if(iand(fields, FIELD_B2I).eq.FIELD_B2I) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   B^-2..."

     b2i79 = 0.

     temp79a = ri2_79* &
          (pstx79(:,OP_DR)**2 + pstx79(:,OP_DZ)**2 + k_fac*bztx79(:,OP_1)**2)

#if defined(USECOMPLEX) || defined(USE3D)
     temp79b = &
          (bfptx79(:,OP_DR)**2 + bfptx79(:,OP_DZ)**2) &
          + 2.*ri_79* &
          (pstx79(:,OP_DZ)*bfptx79(:,OP_DR) - pstx79(:,OP_DR)*bfptx79(:,OP_DZ))

     b2i79(1:npoints,OP_1 ) = 1./(temp79a(1:npoints) + temp79b(1:npoints))
#else
     b2i79(1:npoints,OP_1 ) = 1./temp79a(1:npoints)
#endif
     bi79(1:npoints,OP_1)  = sqrt(b2i79(1:npoints,OP_1))
     b2i79(:,OP_DR) = ri2_79 * &
          (pstx79(:,OP_DR)*pstx79(:,OP_DRR)+pstx79(:,OP_DZ)*pstx79(:,OP_DRZ) &
          +k_fac*bztx79(:,OP_1 )*bztx79(:,OP_DR ))
     b2i79(:,OP_DZ) = ri2_79 * &
          (pstx79(:,OP_DR)*pstx79(:,OP_DRZ)+pstx79(:,OP_DZ)*pstx79(:,OP_DZZ) &
          +k_fac*bztx79(:,OP_1 )*bztx79(:,OP_DZ ))

     if(itor.eq.1) then 
        b2i79(:,OP_DR) = b2i79(:,OP_DR) - ri_79*temp79a
     endif

#if defined(USECOMPLEX) || defined(USE3D)
     b2i79(:,OP_DR) = b2i79(:,OP_DR) + ri_79* &
          (pstx79(:,OP_DZ )*bfptx79(:,OP_DRR) &
          -pstx79(:,OP_DR )*bfptx79(:,OP_DRZ) &
          +pstx79(:,OP_DRZ)*bfptx79(:,OP_DR ) &
          -pstx79(:,OP_DRR)*bfptx79(:,OP_DZ ))&
          +bfptx79(:,OP_DRR)*bfptx79(:,OP_DR) &
          +bfptx79(:,OP_DRZ)*bfptx79(:,OP_DZ)
     b2i79(:,OP_DZ) = b2i79(:,OP_DZ) + ri_79* &
          (pstx79(:,OP_DZ )*bfptx79(:,OP_DRZ) &
          -pstx79(:,OP_DR )*bfptx79(:,OP_DZZ) &
          +pstx79(:,OP_DZZ)*bfptx79(:,OP_DR ) &
          -pstx79(:,OP_DRZ)*bfptx79(:,OP_DZ ))&
          +bfptx79(:,OP_DRZ)*bfptx79(:,OP_DR) &
          +bfptx79(:,OP_DZZ)*bfptx79(:,OP_DZ)

     if(itor.eq.1) then
        b2i79(:,OP_DR) = b2i79(:,OP_DR) - ri2_79* &
          (pstx79(:,OP_DZ)*bfptx79(:,OP_DR) - pstx79(:,OP_DR)*bfptx79(:,OP_DZ))
     endif
#endif

     b2i79(:,OP_DR) = -2.*b2i79(:,OP_DR)*b2i79(:,OP_1)**2
     b2i79(:,OP_DZ) = -2.*b2i79(:,OP_DZ)*b2i79(:,OP_1)**2
  endif

  ! ETA
  ! ~~~
  if(iand(fields, FIELD_ETA).eq.FIELD_ETA) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   eta..."

     if(izone.eq.ZONE_VACUUM) then
        eta79 = 0.
        eta79(:,OP_1) = eta_vac
     else if(izone.eq.ZONE_CONDUCTOR) then
        call get_zone_index(itri,iz)
        izarr = iz
        eta79 = 0.
        eta79(:,OP_1) = wall_resistivity(x_79,phi_79,z_79,izarr)
        etaRZ79 = 0
        etaRZ79(:,OP_1) = wall_resistivityRZ(x_79,phi_79,z_79,izarr)
     else 
        if(iresfunc.eq.2) then
           if(linear.eq.1) then
              tm79 = ((ps079-psimin)/(psibound-psimin) - etaoff) / etadelt
           else
              tm79 = ((pst79-psimin)/(psibound-psimin) - etaoff) / etadelt
           end if
           eta79 = 0.
           eta79(:,OP_1 )  = 1. + tanh(real(tm79(:,OP_1)))
           eta79(:,OP_DR)  = sech(real(tm79(:,OP_1)))**2 * tm79(:,OP_DR)
           eta79(:,OP_DZ)  = sech(real(tm79(:,OP_1)))**2 * tm79(:,OP_DZ)
           eta79(:,OP_DRR) = sech(real(tm79(:,OP_1)))**2 * tm79(:,OP_DRR) &
                - 2.*tanh(real(tm79(:,OP_1)))*sech(real(tm79(:,OP_1)))**2 &
                * tm79(:,OP_DR)**2
           eta79(:,OP_DRZ) = sech(real(tm79(:,OP_1)))**2 * tm79(:,OP_DRZ) &
                - 2.*tanh(real(tm79(:,OP_1)))*sech(real(tm79(:,OP_1)))**2 &
                * tm79(:,OP_DR)*tm79(:,OP_DZ)
           eta79(:,OP_DZZ) = sech(real(tm79(:,OP_1)))**2 * tm79(:,OP_DRR) &
                - 2.*tanh(real(tm79(:,OP_1)))*sech(real(tm79(:,OP_1)))**2 &
                * tm79(:,OP_DZ)**2

           eta79 = (0.5*eta79*eta0 + etar)*eta_fac
        else if(iresfunc.eq.3) then
           temp79a = (pst79(:,OP_1) - psimin)/(psibound - psimin)
           temp79b = (pst79(:,OP_DR)*(x_79 - xmag) &
                +     pst79(:,OP_DZ)*(z_79 - zmag))*(psibound-psimin)
           eta79 = 0.
           where(real(temp79a).lt.etaoff .and. real(temp79b).gt.0.)
              eta79(:,OP_1) = etar
           elsewhere
              eta79(:,OP_1) = eta0
           end where
           eta79 = eta79*eta_fac
           !     else if(iresfunc.eq.0 .or. iresfunc.eq.4) then
        else if(iresfunc.eq.4) then
           ! Here eta79 = 1/T^1.5 .  Factor of efac is included later
           eta79 = 0.
           eta79(:,OP_1) = eta_max / efac

           ! Te
           if(itemp.eq.1) then
              temp79b = tet79(:,OP_1) - eta_te_offset
           else
              temp79b = pet79(:,OP_1)/net79(:,OP_1) - eta_te_offset
           end if

#ifdef USE3D
           ! dTe/dphi
           temp79c = pet79(:,OP_DP)/net79(:,OP_1) - &
                pet79(:,OP_1)*net79(:,OP_DP)/net79(:,OP_1)**2
#endif
           where(real(temp79b).gt.(eta_max/efac)**(-2./3.))
              temp79a = sqrt(temp79b)
              eta79(:,OP_1 ) = 1. / temp79a**3
              eta79(:,OP_DR) = (-3./2.) / temp79a**5 * &
                   (pet79(:,OP_DR)/net79(:,OP_1) &
                   -pet79(:,OP_1)*net79(:,OP_DR)/net79(:,OP_1)**2)
              eta79(:,OP_DZ) = (-3./2.) / temp79a**5 * &
                   (pet79(:,OP_DZ)/net79(:,OP_1) &
                   -pet79(:,OP_1)*net79(:,OP_DZ)/net79(:,OP_1)**2)
              eta79(:,OP_DRR) = (15./4.) / temp79a**7 * &
                   (pet79(:,OP_DR)/net79(:,OP_1) &
                   -pet79(:,OP_1)*net79(:,OP_DR)/net79(:,OP_1)**2)**2 &
                   + (-3./2.) / temp79a**5 * &
                   (pet79(:,OP_DRR)/net79(:,OP_1) &
                   -2.*pet79(:,OP_DR)*net79(:,OP_DR)/net79(:,OP_1)**2 &
                   -pet79(:,OP_1)*net79(:,OP_DRR)/net79(:,OP_1)**2 &
                   +2.*pet79(:,OP_1)*net79(:,OP_DR)**2/net79(:,OP_1)**3)
              eta79(:,OP_DRZ) = (15./4.) / temp79a**7 * &
                   (pet79(:,OP_DR)/net79(:,OP_1) &
                   -pet79(:,OP_1)*net79(:,OP_DR)/net79(:,OP_1)**2) &
                   *(pet79(:,OP_DZ)/net79(:,OP_1) &
                   -pet79(:,OP_1)*net79(:,OP_DZ)/net79(:,OP_1)**2) &
                   + (-3./2.) / temp79a**5 * &
                   (pet79(:,OP_DRZ)/net79(:,OP_1) &
                   -pet79(:,OP_DR)*net79(:,OP_DZ)/net79(:,OP_1)**2 &
                   -pet79(:,OP_DZ)*net79(:,OP_DR)/net79(:,OP_1)**2 &
                   -pet79(:,OP_1)*net79(:,OP_DRZ)/net79(:,OP_1)**2 &
                   +2.*pet79(:,OP_1)*net79(:,OP_DR)*net79(:,OP_DZ) &
                   /net79(:,OP_1)**3)
              eta79(:,OP_DZZ) = (15./4.) / temp79a**7 * &
                   (pet79(:,OP_DZ)/net79(:,OP_1) &
                   -pet79(:,OP_1)*net79(:,OP_DZ)/net79(:,OP_1)**2)**2 &
                   + (-3./2.) / temp79a**5 * &
                   (pet79(:,OP_DZZ)/net79(:,OP_1) &
                   -2.*pet79(:,OP_DZ)*net79(:,OP_DZ)/net79(:,OP_1)**2 &
                   -pet79(:,OP_1)*net79(:,OP_DZZ)/net79(:,OP_1)**2 &
                   +2.*pet79(:,OP_1)*net79(:,OP_DZ)**2/net79(:,OP_1)**3)
#ifdef USE3D
              eta79(:,OP_DP) = -(3./2.)*temp79c / temp79a**5
#endif
           end where

           eta79 = eta79 * efac

           call calculate_zeff(itri,temp79a)
           do i=1, OP_NUM
              eta79(:,i) = eta79(:,i) * temp79a
           end do
        else
           call eval_ops(itri, resistivity_field, eta79)
        end if

        where(eta79.ne.eta79) eta79 = 0.
        where(real(eta79(:,OP_1)).lt.eta_min) eta79(:,OP_1) = eta_min
        where(real(eta79(:,OP_1)).gt.eta_max) eta79(:,OP_1) = eta_max
     end if

  end if

  ! denm
  ! ~~~~
  if(iand(fields, FIELD_DENM).eq.FIELD_DENM) then

     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   denm..."
     
     if(idenmfunc.eq.1) then
        denm79 = 0.
        if(izone.eq.1. .and. denmt.gt.0) then
           ! Te
           if(itemp.eq.1) then
              temp79b = tet79(:,OP_1)
           else
              temp79b = pet79(:,OP_1)/net79(:,OP_1)
           end if

           where(real(temp79b).gt.(denmt/(denmmax-denm)))
              denm79(:,OP_1) = net79(:,OP_1)/pet79(:,OP_1)
              denm79(:,OP_DR) = net79(:,OP_DR)/pet79(:,OP_1) &
                   - net79(:,OP_1)*pet79(:,OP_DR)/pet79(:,OP_1)**2
              denm79(:,OP_DZ) = net79(:,OP_DZ)/pet79(:,OP_1) &
                   - net79(:,OP_1)*pet79(:,OP_DZ)/pet79(:,OP_1)**2
              denm79(:,OP_DRR) = net79(:,OP_DRR)/pet79(:,OP_1) &
                   - 2.*net79(:,OP_DR)*pet79(:,OP_DR)/pet79(:,OP_1)**2 &
                   + 2.*net79(:,OP_1)*pet79(:,OP_DR)**2/pet79(:,OP_1)**3 &
                   - net79(:,OP_1)*pet79(:,OP_DRR)/pet79(:,OP_1)**2
              denm79(:,OP_DRZ) = net79(:,OP_DRZ)/pet79(:,OP_1) &
                   - (net79(:,OP_DR)*pet79(:,OP_DZ) + &
                      net79(:,OP_DZ)*pet79(:,OP_DR))/pet79(:,OP_1)**2 &
                   + 2.*net79(:,OP_1)*pet79(:,OP_DR)*pet79(:,OP_DZ)/pet79(:,OP_1)**3 &
                   - net79(:,OP_1)*pet79(:,OP_DRZ)/pet79(:,OP_1)**2
              denm79(:,OP_DZZ) = net79(:,OP_DZZ)/pet79(:,OP_1) &
                   - 2.*net79(:,OP_DZ)*pet79(:,OP_DZ)/pet79(:,OP_1)**2 &
                   + 2.*net79(:,OP_1)*pet79(:,OP_DZ)**2/pet79(:,OP_1)**3 &
                   - net79(:,OP_1)*pet79(:,OP_DZZ)/pet79(:,OP_1)**2
#ifdef USE3D
              denm79(:,OP_DP) = net79(:,OP_DP)/pet79(:,OP_1) &
                   - net79(:,OP_1)*pet79(:,OP_DP)/pet79(:,OP_1)**2
#endif
           end where

           denm79 = denm79*denmt
           denm79(:,OP_1) = denm79(:,OP_1) + denm

           where((denmmax-denm)*real(temp79b).le.denmt)
              denm79(:,OP_1) = denmmax
           end where
           where((denmmin-denm)*real(temp79b).ge.denmt)
              denm79(:,OP_1) = denmmin
           end where

           where(denm79.ne.denm79) denm79 = 0.
           temp79a = denm79(:,OP_1)
           do i=1, OP_NUM
              where(real(temp79a).lt.0.) denm79(:,i) = 0.
              where(real(temp79a).gt.denmmax) denm79(:,i) = 0.
              where(real(temp79a).lt.denmmin) denm79(:,i) = 0.
           end do
           where(real(temp79a).gt.denmmax) denm79(:,OP_1) = denmmax
           where(real(temp79a).lt.denmmin) denm79(:,OP_1) = denmmin

        end if
     else
        call eval_ops(itri, denm_field, denm79)
     end if

  end if


  ! KAP
  ! ~~~
  if(iand(fields, FIELD_KAP).eq.FIELD_KAP) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   kappa..."

     if(ikappafunc.eq.5) then
        kap79 = 0.
        if(izone.eq.1. .and. kappa0.gt.0) then
           ! Te
           if(itemp.eq.1) then
              temp79b = tet79(:,OP_1)
           else
              temp79b = pet79(:,OP_1)/net79(:,OP_1)
           end if

           kap79 = 0.
           where(real(temp79b).gt.(kappa0/(kappa_max-kappat)))
              kap79(:,OP_1) = net79(:,OP_1)/pet79(:,OP_1)
              kap79(:,OP_DR) = net79(:,OP_DR)/pet79(:,OP_1) &
                   - net79(:,OP_1)*pet79(:,OP_DR)/pet79(:,OP_1)**2
              kap79(:,OP_DZ) = net79(:,OP_DZ)/pet79(:,OP_1) &
                   - net79(:,OP_1)*pet79(:,OP_DZ)/pet79(:,OP_1)**2
              kap79(:,OP_DRR) = net79(:,OP_DRR)/pet79(:,OP_1) &
                   - 2.*net79(:,OP_DR)*pet79(:,OP_DR)/pet79(:,OP_1)**2 &
                   + 2.*net79(:,OP_1)*pet79(:,OP_DR)**2/pet79(:,OP_1)**3 &
                   - net79(:,OP_1)*pet79(:,OP_DRR)/pet79(:,OP_1)**2
              kap79(:,OP_DRZ) = net79(:,OP_DRZ)/pet79(:,OP_1) &
                   - (net79(:,OP_DR)*pet79(:,OP_DZ) + &
                      net79(:,OP_DZ)*pet79(:,OP_DR))/pet79(:,OP_1)**2 &
                   + 2.*net79(:,OP_1)*pet79(:,OP_DR)*pet79(:,OP_DZ)/pet79(:,OP_1)**3 &
                   - net79(:,OP_1)*pet79(:,OP_DRZ)/pet79(:,OP_1)**2
              kap79(:,OP_DZZ) = net79(:,OP_DZZ)/pet79(:,OP_1) &
                   - 2.*net79(:,OP_DZ)*pet79(:,OP_DZ)/pet79(:,OP_1)**2 &
                   + 2.*net79(:,OP_1)*pet79(:,OP_DZ)**2/pet79(:,OP_1)**3 &
                   - net79(:,OP_1)*pet79(:,OP_DZZ)/pet79(:,OP_1)**2
#ifdef USE3D
              kap79(:,OP_DP) = net79(:,OP_DP)/pet79(:,OP_1) &
                   - net79(:,OP_1)*pet79(:,OP_DP)/pet79(:,OP_1)**2
#endif
           end where

           kap79 = kap79*kappa0
           kap79(:,OP_1) = kap79(:,OP_1) + kappat

           where(real(temp79b).le.(kappa0/(kappa_max-kappat)))
              kap79(:,OP_1) = kappa_max
           end where

           if(kappaf.ge.0. .and. gradp_crit.ne.0) then
              temp79a = pt79(:,OP_DR)**2 + pt79(:,OP_DZ)**2
#ifdef USE3D
              temp79a = temp79a + ri2_79*pt79(:,OP_DP)**2
#endif
              do i=1, OP_NUM
                 where(real(temp79a).lt.gradp_crit**2) kap79(:,i) = kap79(:,i) * kappaf
              end do
           end if

           if(kappah.ne.0.) then
              tm79 = (pst79 - psimin)/(psibound - psimin)
              temp79a = tanh((real(tm79(:,OP_1))-1.)/.2)
              temp79b = cosh((real(tm79(:,OP_1))-1.)/.2)**(-1)
              kap79(:,OP_1) = kap79(:,OP_1) + kappah*temp79a**2
              kap79(:,OP_DR) = kap79(:,OP_DR) &
                   + kappah*temp79a*temp79b**2*tm79(:,OP_DR)
              kap79(:,OP_DZ) = kap79(:,OP_DZ) &
                   + kappah*temp79a*temp79b**2*tm79(:,OP_DZ)
              kap79(:,OP_DRR) = kap79(:,OP_DRR) &
                   + kappah*(temp79a*temp79b**2*tm79(:,OP_DRR) &
                             + 0.5*(temp79b**4)*(tm79(:,OP_DR)**2) &
                             - (temp79a**2)*(temp79b**2)*(tm79(:,OP_DR)**2))
              kap79(:,OP_DRZ) = kap79(:,OP_DRZ) &
                   + kappah*(temp79a*temp79b**2*tm79(:,OP_DRZ) &
                             + 0.5*(temp79b**4)*(tm79(:,OP_DR)**2) &
                             - (temp79a**2)*(temp79b**2)*(tm79(:,OP_DR)**2))
              kap79(:,OP_DZZ) = kap79(:,OP_DZZ) &
                   + kappah*(temp79a*temp79b**2*tm79(:,OP_DZZ) &
                             + 0.5*(temp79b**4)*tm79(:,OP_DR)*tm79(:,OP_DZ) &
                             - (temp79a**2)*(temp79b**2)*tm79(:,OP_DR)*tm79(:,OP_DZ))
#ifdef USE3D
              kap79(:,OP_DP) = kap79(:,OP_DP) &
                   + kappah*temp79a*temp79b**2*tm79(:,OP_DP)
#endif
           end if

           where(kap79.ne.kap79) kap79 = 0.
           temp79a = kap79(:,OP_1)
           do i=1, OP_NUM
              where(real(temp79a).lt.0.) kap79(:,i) = 0.
              where(real(temp79a).gt.kappa_max) kap79(:,i) = 0.
           end do
           where(real(temp79a).gt.kappa_max) kap79(:,OP_1) = kappa_max

        end if
     else
        call eval_ops(itri, kappa_field, kap79)
     end if

     if(ikapscale.eq.1) then
        kar79 = kappar*kap79
     else if(ikapparfunc.eq.2) then

        ! Here eta79 = 1/T^1.5 .  Factor of efac is included later
        kar79 = 0.
        if(izone.eq.1) then

           kar79(:,OP_1) = kappar_max / krfac

           ! Te
           temp79b = tet79(:,OP_1)
           
           kr_tmin = (kappar_min/krfac)**0.4
           kr_tmax = (kappar_max/krfac)**0.4
           if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "kr_tlims", kr_tmin, kr_tmax
           where(real(temp79b).lt.(kappar_max/krfac)**0.4 .and. &
                 real(temp79b).gt.(kappar_min/krfac)**0.4)
              temp79a = sqrt(temp79b)
              kar79(:,OP_1 ) = temp79a**5
              kar79(:,OP_DR) = (5./2.) * temp79a**3 * tet79(:,OP_DR)
              kar79(:,OP_DZ) = (5./2.) * temp79a**3 * tet79(:,OP_DZ)
              kar79(:,OP_DRR) = (15./4.) * temp79a * tet79(:,OP_DR)**2 &
                   + (5./2.) * temp79a**3 * tet79(:,OP_DRR)
              kar79(:,OP_DRZ) = (15./4.) * temp79a * tet79(:,OP_DR) * tet79(:,OP_DZ)  &
                   + (5./2.) * temp79a**3 * tet79(:,OP_DRZ)
              kar79(:,OP_DZZ) = (15./4.) * temp79a *  tet79(:,OP_DZ)**2 &
                   + (5./2.) * temp79a**3 * tet79(:,OP_DZZ)
              
#ifdef USE3D
              kar79(:,OP_DP) = (5./2.) * tet79(:,OP_DP) * temp79a**3
#endif
           end where

           kar79 = kar79 * krfac
           where(kar79.ne.kar79) kar79 = 0.
           where(real(kar79(:,OP_1)).lt.kappar_min) kar79(:,OP_1) = kappar_min
           where(real(kar79(:,OP_1)).gt.kappar_max) kar79(:,OP_1) = kappar_max
        end if
     else
        call eval_ops(itri, kappar_field, kar79)
     end if

     kax79 = 0.
     kax79(:,OP_1) = kappax
  end if

  ! SIG
  ! ~~~
  if((iand(fields, FIELD_SIG).eq.FIELD_SIG) .and. idens.eq.1 &
       .and. density_source) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   sigma..."

     call eval_ops(itri, sigma_field, sig79)
  else
     sig79 = 0.
  end if

  ! F
  ! ~
  if((iand(fields, FIELD_F).eq.FIELD_F) &
       .and. momentum_source) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   F..."

     call eval_ops(itri, Fphi_field, fy79)
  else
     fy79 = 0.
  end if

  ! Q
  ! ~
  if((iand(fields, FIELD_Q).eq.FIELD_Q) &
       .and. heat_source) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   Q..."

     call eval_ops(itri, Q_field, q79)
  else
     q79 = 0.
  end if

  ! Rad
  ! ~
  if((iand(fields, FIELD_RAD).eq.FIELD_RAD) &
       .and. rad_source) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   Rad..."

     call eval_ops(itri, Totrad_field,  totrad79)
     call eval_ops(itri, Linerad_field, linerad79)
     call eval_ops(itri, Bremrad_field, bremrad79)
     call eval_ops(itri, Ionrad_field,  ionrad79)
     call eval_ops(itri, Reckrad_field, reckrad79)
     call eval_ops(itri, Recprad_field, recprad79)
     
  else
     totrad79 = 0.
     linerad79 = 0.
     bremrad79 = 0.
     ionrad79 = 0.
     reckrad79 = 0.
     recprad79 = 0.
  end if

  ! cd
  ! ~
  if((iand(fields, FIELD_CD).eq.FIELD_CD) &
       .and. icd_source.gt.0) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   cd..."

     call eval_ops(itri, cd_field, cd79)
  else
     cd79 = 0.
  end if


  ! Jbs Coefs
  ! ~
  if((iand(fields, FIELD_JBS).eq.FIELD_JBS) &
       .and. ibootstrap.gt.0) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   Jbs_coefs..."

     call eval_ops(itri, Jbs_fluxavg_iBsq_field, jbsfluxavg_iBsq_B79)
     call eval_ops(itri, Jbs_fluxavg_G_field, jbsfluxavg_G79)
     if(ibootstrap.eq.2 .or. ibootstrap .eq.1) then
      call eval_ops(itri, Jbs_L31_field, jbsl3179)
      call eval_ops(itri, Jbs_L32_field, jbsl3279)
      call eval_ops(itri, Jbs_L34_field, jbsl3479)
      call eval_ops(itri, Jbs_alpha_field, jbsalpha79)
      call eval_ops(itri, Jbs_dtedpsit_field, jbs_dtedpsit79)
     endif
     if(ibootstrap.eq.3) then
      jbsl3179 = 0.
      jbsl3279 = 0.
      jbsl3479 = 0.
      jbsalpha79 = 0.
      call eval_ops(itri, Jbs_dtedpsit_field, jbs_dtedpsit79)
      call eval_ops(itri, Jbs_ftrap_field, jbs_ftrap79)
      call eval_ops(itri, Jbs_qR_field, jbs_qR79)
      call eval_ops(itri, Jbs_invAspectRatio_field, jbs_invAspectRatio79)
     endif
     
  else
     jbsl3179 = 0.
     jbsl3279 = 0.
     jbsl3479 = 0.
     jbsalpha79 = 0.
     jbsfluxavg_iBsq_B79 = 0.
     jbsfluxavg_G79 = 0.
     if(ibootstrap.eq.2) then
      jbs_dtedpsit79 =0.
     endif
     if(ibootstrap.eq.3) then
      jbs_dtedpsit79 =0.
      jbs_ftrap79 =0.
      jbs_qR79 =0.
      jbs_invAspectRatio79 =0.
     endif
  end if

  ! Wall dist field
  ! ~~~~~~~~~~~~~~~
  if((iand(fields, FIELD_WALL).eq.FIELD_WALL)) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   Wall dist..."

     call eval_ops(itri, wall_dist, wall79)
  endif

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
     else if(ivisfunc.eq.4) then
        where(real(wall79(:,OP_1)).gt.amuoff)
           vis79(:,OP_1) = amu
           vic79(:,OP_1) = amuc
        elsewhere
           vis79(:,OP_1) = amu_edge
           vic79(:,OP_1) = amu_edge
        end where
     else
        call eval_ops(itri, visc_field, vis79)

        if(numvar.ge.3) then
           call eval_ops(itri, visc_c_field, vic79)
        endif
     endif

!     if(amupar.ne.0.) vip79 = amupar*pit79/2.
     if(amupar.ne.0.) vip79 = amupar
  end if


  ! Poloidal Momentum Force
  ! ~~~
  if((iand(fields, FIELD_PF).eq.FIELD_PF)   &
      .and. ipforce .gt. 0) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   PFORCE..."

      call eval_ops(itri, pforce_field, for79)
  endif


  ! Electrostatic Potential
  ! ~~~
  if((iand(fields, FIELD_ES).eq.FIELD_ES)   &
      .and. jadv.eq.0) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   potential..."

      call eval_ops(itri, e_field(1), es179)
  endif

#ifdef USEPARTICLES
    ! Kinetic Pressure Terms
    ! ~~~
    if((iand(fields, FIELD_KIN).eq.FIELD_KIN)   &
        .and. kinetic .eq. 1) then
       if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   kinetic..."
        call eval_ops(itri, p_f_par, pfpar79, rfac)
        call eval_ops(itri, p_f_perp, pfper79, rfac)
        call eval_ops(itri, pf_field, pf079, rfac)

        if (kinetic_fast_ion.eq.1) then
           call eval_ops(itri, p_f_par, pfpar79, rfac)
           call eval_ops(itri, p_f_perp, pfper79, rfac)
           call eval_ops(itri, pf_field, pf079, rfac)

           call eval_ops(itri, den_f_1, nf79, rfac)
           call eval_ops(itri, nf_field, nf079, rfac)

           call eval_ops(itri, v_f_par, vfpar79, rfac)
           !call eval_ops(itri, v_f_par_0, vfpar079, rfac)
           !pf079 = 0.
        else
           pfpar79 = 0.
           pfper79 = 0.
           pf079 = 0.
           nf79 = 0.
           nf079 = 0.
           vfpar79 = 0.
        endif

        if (kinetic_thermal_ion.eq.1) then
           call eval_ops(itri, p_i_par, pipar79, rfac)
           call eval_ops(itri, p_i_perp, piper79, rfac)
           call eval_ops(itri, pfi_field, pfi079, rfac)

           call eval_ops(itri, den_i_1, nfi79, rfac)
           call eval_ops(itri, nfi_field, nfi079, rfac)

           call eval_ops(itri, v_i_par, vipar79, rfac)
        else
           pipar79 = 0.
           piper79 = 0.
           pfi079 = 0.
           nfi79 = 0.
           nfi079 = 0.
           vipar79 = 0.
        endif
        call eval_ops(itri, rho_field, rhof79, rfac)
           ! pipar79 = 0.
           ! piper79 = 0.
           ! pfi079 = 0.
         
        !do ipoint=1,MAX_PTS
        !      if (real(rhof79(ipoint,OP_1))>0.85) then
        !   pipar79(ipoint,:)=0.
        !   piper79(ipoint,:)=0.
        !   !nfi79(ipoint,:)=0.
        !   endif
        !enddo

    endif
#endif

  ! Runaway Electron Density
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  if((iand(fields, FIELD_RE).eq.FIELD_RE) .and. irunaway.gt.0) then
     if(itri.eq.1 .and. myrank.eq.0 .and. iprint.ge.2) print *, "   RE density..."

     if(ilin.eq.0) then
       call eval_ops(itri, nre_field(1), nre179, rfac)
     else
       nre179 = 0.
     endif
     if(ieqsub.eq.1) then
        if(idenfunc.eq.3) then
           if(izone.gt.1) then
              nre079(:,OP_1) = 0.
           else
              temp79a = (pst79(:,OP_1) - psimin)/(psibound - psimin)
              temp79b = (pst79(:,OP_DR)*(x_79 - xmag) &
                   +     pst79(:,OP_DZ)*(z_79 - zmag))*(psibound-psimin)
              where(real(temp79a).ge.denoff .and. real(temp79b).le.0.)
                 nre079(:,OP_1) = 0.
              end where
           end if
        else
           call eval_ops(itri, nre_field(0), nre079)
        end if
     else
        nre079 = 0.
     endif
  endif

end subroutine define_fields

end module m3dc1_nint

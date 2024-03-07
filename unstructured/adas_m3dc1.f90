!========================
! M3D-C1 ADAS connection
!========================
module adas_m3dc1

  character(len=256) :: adas_adf11 ! root folder location to ADAS adf11 data

  integer, parameter :: izdimd = 20, itdimd = 50, iddimd = 50
  integer :: iz0, iz1min, iz1max

  ! for radiation from plt data
  integer :: iblmx_plt, itmax_plt, idmax_plt
  integer :: isstgr_plt(izdimd)
  real :: ddens_plt(iddimd), dtev_plt(itdimd), drcof_plt(izdimd, itdimd, iddimd)

  ! for ionization from scd data
  integer :: iblmx_scd, itmax_scd, idmax_scd
  integer :: isstgr_scd(izdimd)
  real :: ddens_scd(iddimd), dtev_scd(itdimd), drcof_scd(izdimd, itdimd, iddimd)

contains

  subroutine load_adf11(Z)
    use basic
    implicit none
    integer, intent(in) :: Z

    integer :: iclass
    character(len=19) :: pltname, scdname
    character(len=255) :: dsname

    integer, parameter :: iunit = 67
    integer, parameter :: idptn = 128, idptnc = 256, idcnct = 100, idptnl = 4
    logical :: lexist, lres, lstan, lptn
    integer :: nptnl, ncnct, ismax
    integer :: iptnca(idptnl, idptn, idptnc), ispbr(izdimd), isppr(izdimd)
    integer :: icnctv(idcnct), iptnla(idptnl), iptna(idptnl, idptn)
    integer :: nptn(idptnl), nptnc(idptnl, idptn)
    character(len=12) :: dnr_ele
    real :: dnr_ams

    ! using high quality data unless otherwise noted by necessity
    select case (Z)
       
    case (2) ! HELIUM
       pltname = '/plt96/plt96_he.dat'
       scdname = '/scd96/scd96_he.dat'
    case (4) ! BERYLLIUM
       pltname = '/plt96/plt96_be.dat'
       scdname = '/scd96/scd96_be.dat'
    case (5) ! BORON
       pltname = '/plt89/plt89_b.dat' ! low quality (only available)
       scdname = '/scd89/scd89_b.dat' ! low quality (only available)
    case (6) ! CARBON
       pltname = '/plt96/plt96_c.dat'
       scdname = '/scd96/scd96_c.dat'
    case (10) ! NEON
       pltname = '/plt96/plt96_ne.dat'
       scdname = '/scd96/scd96_ne.dat'
    case (18) ! ARGON
       pltname = '/plt89/plt89_ar.dat' ! low quality (only available)
       scdname = '/scd85/scd85_ar.dat' ! medium quality (could use Summers [scd74] or this one)
    case DEFAULT
       if(myrank.eq.0) write(*,*) 'THIS ELEMENT NOT YET IMPLEMENTED WITH ADAS'
       call safestop(1001)
    end select

    call xx0000

    ! Load radiation data
    dsname = trim(adas_adf11) // trim(pltname)
    iclass = 8
    
    inquire(FILE=dsname, EXIST=lexist)
    if (lexist) then
       if(myrank.eq.0) print *, trim(dsname), " Found!"
       open(unit=iunit, file=dsname, status='old')
    else
       if(myrank.eq.0) print *, trim(dsname), " does not exist"
       call safestop(1001)
    endif

    call xxdata_11(iunit, iclass, izdimd, iddimd, itdimd, idptnl, idptn, &
         idptnc, idcnct, iz0, iz1min, iz1max, nptnl, nptn, nptnc, &
         iptnla, iptna, iptnca, ncnct, icnctv, iblmx_plt, ismax, dnr_ele, &
         dnr_ams, isppr, ispbr, isstgr_plt, idmax_plt, itmax_plt, ddens_plt, &
         dtev_plt, drcof_plt, lres, lstan, lptn)

    close(iunit)

    ! Load ionization data
    dsname = trim(adas_adf11) // trim(scdname)
    iclass = 2
    
    inquire(FILE=dsname, EXIST=lexist)
    if (lexist) then
       if(myrank.eq.0) print *, trim(dsname), " Found!"
       open(unit=iunit, file=dsname, status='old')
    else
       if(myrank.eq.0) print *, trim(dsname), " does not exist"
       call safestop(1001)
    endif

    call xxdata_11(iunit, iclass, izdimd, iddimd, itdimd, idptnl, idptn, &
         idptnc, idcnct, iz0, iz1min, iz1max, nptnl, nptn, nptnc, &
         iptnla, iptna, iptnca, ncnct, icnctv, iblmx_scd, ismax, dnr_ele, &
         dnr_ams, isppr, ispbr, isstgr_scd, idmax_scd, itmax_scd, ddens_scd, &
         dtev_scd, drcof_scd, lres, lstan, lptn)

    close(iunit)

    return
  end subroutine load_adf11


  subroutine interp_adf11(iclass, N, iz1, te, dens, coeff)
    use basic
    implicit none
    integer, parameter :: lck = 51
    integer, intent(in) :: iclass, N, iz1
    real, intent(in) :: te(N), dens(N)
    real, intent(out) :: coeff(N)

    integer :: it, id, isel, indsel, isdat
    integer :: nin
    logical :: lsetx
    real :: ytot
    logical :: ltrng(lck), ldrng(lck)
    real :: dtev(lck), ddens(lck)
    real :: yin(lck), yout(lck), ypass(lck,lck), dy(lck), xin(lck)
    real :: r8fun1

    integer :: iblmx, itmax, idmax
    integer :: isstgr(izdimd)
    real :: ddens_arr(iddimd), dtev_arr(itdimd), drcof_arr(izdimd, itdimd, iddimd)
    external :: r8fun1

    if (N.gt.lck) then
       if(myrank.eq.0) print *, 'ADAS ERROR: N > LCK - INCREASE LCK'
       call safestop(1002)
    endif
    if (iz1 > iz0 + 1) then
       if(myrank.eq.0) print *, 'ADAS ERROR: charge requested > nuclear charge'
       call safestop(1002)
    endif

    if (iz1 < iz1min) then
       if(myrank.eq.0) print *, 'ADAS ERROR: charge outside range (too low)'
       call safestop(1002)
    endif

    if (iz1 > iz1max) then
       if(myrank.eq.0) print *, 'ADAS ERROR: charge outside range (too high)'
       call safestop(1002)
    endif

    ! ionization from scd data
    if (iclass == 2) then
       iblmx  = iblmx_scd
       itmax  = itmax_scd
       idmax  = idmax_scd
       isstgr = isstgr_scd
       ddens_arr  = ddens_scd
       dtev_arr   = dtev_scd
       drcof_arr  = drcof_scd

    ! line radiation from plt data
    elseif (iclass == 8) then
       iblmx  = iblmx_plt
       itmax  = itmax_plt
       idmax  = idmax_plt
       isstgr = isstgr_plt
       ddens_arr  = ddens_plt
       dtev_arr   = dtev_plt
       drcof_arr  = drcof_plt
    else
       if(myrank.eq.0) print *, 'ADAS ERROR: invalid ADAS class type. Must be 2 (scd) or 8 (plt)'
       call safestop(1002)
    end if

    indsel = 0
    do isel = 1, iblmx
       isdat = isstgr(isel)
       if (isdat == iz1) indsel = isel
    end do

    if (indsel == 0) then
       if(myrank.eq.0) print *, 'ADAS ERROR: no valid data block'
       call safestop(1002)
    endif

    ! Transform to log space
    do it = 1, N
       dtev(it) = dlog10(te(it))
       ddens(it) = dlog10(dens(it))
    end do

    lsetx = .true.

    ! At every density location in the file, do a 1D temperature interpolation
    do id = 1, idmax
       ytot = 0.0D0
       do it = 1, itmax
          xin(it) = dtev_arr(it)
          yin(it) = drcof_arr(indsel, it, id)
          ytot = ytot + yin(it)
       end do
       nin = itmax

       nin = 0
       do it =1, itmax
          if (drcof_arr(indsel,it,id).GT.-40) then
             xin(nin+1) = dtev_arr(it)
             yin(nin+1) = drcof_arr(indsel,it,id)
             nin = nin + 1
          endif
       end do

       if (nin > 1) then
          call xxsple(lsetx, 0, r8fun1, nin, xin, yin, N, dtev, yout, dy, ltrng)
       else
          yout = -74.0d0
       endif

       do it = 1, N
          ypass(it,id) = yout(it)
       end do

    end do

    lsetx = .true.

    ! Now do the 1D density interpolation
    do it = 1, N
       do id = 1, idmax
          yin(id) = ypass(it,id)
       end do
       call xxsple(lsetx, 0, r8fun1, idmax, ddens_arr, yin, 1, ddens(it), yout(it), dy, ldrng(it))
    end do
    coeff = yout(1:N)
    return

  end subroutine
  

end module adas_m3dc1

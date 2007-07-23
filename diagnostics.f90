module diagnostics

  implicit none

  real :: tflux0, totcur0

  ! scalars integrated over entire computational domain
  real :: tflux, area, totcur, totden, tmom, tvor, psilim

  ! scalars integrated within lcfs
  real :: pflux, parea, pcur, pden, pmom, pvor

  real :: chierror

  ! energy diagnostics
  real :: ekin, emag, ekind, emagd, ekino, emago, ekindo, emagdo,      &
       ekint,emagt,ekintd,emagtd,ekinto,emagto,ekintdo,emagtdo,        &
       ekinp,emagp,ekinpd,emagpd,ekinpo,emagpo,ekinpdo,emagpdo,        &
       ekinph,ekinth,emagph,emagth,ekinpho,ekintho,emagpho,emagtho,    &
       ekin3,ekin3d,ekin3h,emag3,ekin3o,ekin3do,ekin3ho,emag3o,        &
       emag3h,emag3d,emag3ho,emag3do
  real :: efluxd,efluxp,efluxk,efluxs,efluxt,epotg,etot,ptot,eerr,ptoto

  ! timing diagnostics
  real :: t_ludefall, t_sources, t_smoother, t_aux, t_onestep
  real :: t_solve_v, t_solve_n, t_solve_p, t_solve_b, t_mvm
  real :: t_output_cgm, t_output_hdf5, t_output_reset
  real :: t_gs, t_gs_magaxis, t_gs_fundef, t_gs_solve, t_gs_init

contains

  ! ======================================================================
  ! reset_timings
  ! 
  ! resents timing vaiables to 0
  ! ======================================================================
  subroutine reset_timings
    implicit none
    t_ludefall = 0.
    t_sources = 0.
    t_smoother = 0.
    t_aux = 0.
    t_solve_v = 0.
    t_solve_n = 0.
    t_solve_p = 0.
    t_solve_b = 0.
    t_output_cgm = 0.
    t_output_hdf5 = 0.
    t_output_reset = 0.
    t_mvm = 0.
    t_onestep = 0.
  end subroutine reset_timings


  ! ======================================================================
  ! distribute_timings
  ! 
  ! distributes timing vaiables to all processors
  ! ======================================================================
  subroutine distribute_timings

    implicit none

    include 'mpif.h'
    integer :: ier
    real, dimension(13) :: vin, vout

    vin(1) =  t_ludefall
    vin(2) =  t_sources
    vin(3) =  t_smoother
    vin(4) =  t_aux
    vin(5) =  t_solve_v
    vin(6) =  t_solve_n
    vin(7) =  t_solve_p
    vin(8) =  t_solve_b
    vin(9) =  t_output_cgm
    vin(10) = t_output_hdf5
    vin(11) = t_output_reset
    vin(12) = t_mvm
    vin(13) = t_onestep
    call MPI_ALLREDUCE(vin, vout, 13, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, ier)
    t_ludefall      = vout(1)
    t_sources       = vout(2)
    t_smoother      = vout(3)
    t_aux           = vout(4)
    t_solve_v       = vout(5)
    t_solve_n       = vout(6)
    t_solve_p       = vout(7)
    t_solve_b       = vout(8)
    t_output_cgm    = vout(9)
    t_output_hdf5   = vout(10)
    t_output_reset  = vout(11)
    t_mvm           = vout(12)
    t_onestep       = vout(13)
    
  end subroutine distribute_timings


  ! ======================================================================
  ! reset_scalars
  ! -------------
  ! 
  ! resets diagnostic energy and scalar quantities to zero
  ! ======================================================================
  subroutine reset_scalars()
    implicit none

    ekin = 0.
    emag = 0.
    ekind = 0.
    emagd = 0.
    ekinp = 0.
    emagp = 0.
    ekinpd = 0.
    emagpd = 0.      
    ekint = 0.
    emagt = 0.
    ekintd = 0.
    emagtd = 0.
    ekinph = 0.
    ekinth = 0.
    emagph = 0.
    emagth = 0.
    ekin3 = 0.
    ekin3d = 0.
    ekin3h = 0.
    emag3 = 0.
    emag3d = 0.
    emag3h = 0.
    efluxd = 0.
    efluxp = 0.
    efluxk = 0.
    efluxs = 0.
    efluxt = 0.
    epotg = 0.

    area = 0.
    totcur = 0.
    tflux = 0.
    totden = 0.
    tmom = 0.
    tvor = 0.
    parea = 0.
    pcur = 0.
    pflux = 0.
    pden = 0.
    pmom = 0.
    pvor = 0.

  end subroutine reset_scalars


  ! ======================================================================
  ! distribute_scalars
  ! 
  ! distributes diagnostic energy and scalar quantities
  ! ======================================================================
  subroutine distribute_scalars()    
    use basic

    implicit none

    include 'mpif.h'

    integer, parameter :: num_scalars = 36
    integer :: ier
    double precision, dimension(num_scalars) :: temp, temp2

    ! Allreduce energy terms
    if(maxrank .gt. 1) then
       temp(1) = ekinp
       temp(2) = emagp
       temp(3) = ekinpd
       temp(4) = emagpd      
       temp(5) = ekint
       temp(6) = emagt
       temp(7) = ekintd
       temp(8) = emagtd
       temp(9) = ekinph
       temp(10) = ekinth
       temp(11) = emagph
       temp(12) = emagth
       temp(13) = ekin3
       temp(14) = ekin3d
       temp(15) = ekin3h
       temp(16) = emag3
       temp(17) = emag3d
       temp(18) = emag3h
       temp(19) = efluxd
       temp(20) = efluxp
       temp(21) = efluxk
       temp(22) = efluxs
       temp(23) = efluxt
       temp(24) = epotg
       temp(25) = area
       temp(26) = totcur
       temp(27) = totden
       temp(28) = tflux
       temp(29) = tmom
       temp(30) = tvor
       temp(31) = parea
       temp(32) = pcur
       temp(33) = pflux
       temp(34) = pden
       temp(35) = pmom
       temp(36) = pvor
         
       !checked that this should be MPI_DOUBLE_PRECISION
       call mpi_allreduce(temp, temp2, num_scalars, MPI_DOUBLE_PRECISION,  &
            MPI_SUM, MPI_COMM_WORLD, ier) 
         
       ekinp =  temp2( 1)
       emagp =  temp2( 2)
       ekinpd = temp2( 3)
       emagpd = temp2( 4)      
       ekint =  temp2( 5)
       emagt =  temp2( 6)
       ekintd = temp2( 7)
       emagtd = temp2( 8)
       ekinph = temp2( 9)
       ekinth = temp2(10)
       emagph = temp2(11)
       emagth = temp2(12)
       ekin3 =  temp2(13)
       ekin3d = temp2(14)
       ekin3h = temp2(15)
       emag3 =  temp2(16)
       emag3d = temp2(17)
       emag3h = temp2(18)
       efluxd = temp2(19)
       efluxp = temp2(20)
       efluxk = temp2(21)
       efluxs = temp2(22)
       efluxt = temp2(23)
       epotg =  temp2(24)
       area =   temp2(25)
       totcur = temp2(26)
       totden = temp2(27)
       tflux =  temp2(28)
       tmom =   temp2(29)
       tvor =   temp2(30)
       parea =  temp2(31)
       pcur  =  temp2(32)
       pflux =  temp2(33)
       pden =   temp2(34)
       pmom =   temp2(35)
       pvor =   temp2(36)

    endif !if maxrank .gt. 1

  end subroutine distribute_scalars


  ! ======================================================================
  ! reconnected_flux
  !
  ! calculates the reconnected flux assuming an x-point at the center of 
  ! the computational domain, and an o-point at the center of the right
  ! boundary.
  ! ======================================================================
  real function reconnected_flux()
    
    use basic
    use arrays
  
    implicit none

    real :: alx, alz
    real, dimension(2) :: temp, temp2
    integer :: itri

    call getboundingboxsize(alx,alz)

    itri = 0
    call evaluate(alx   ,alz/2.,temp(1),temp2(1),phi,1,numvar,itri)
    itri = 0
    call evaluate(alx/2.,alz/2.,temp(2),temp2(2),phi,1,numvar,itri)

    reconnected_flux = 0.5*(temp(2)-temp(1))

    return
  end function reconnected_flux

! ======================================================================
! second
!
! returns the current time in seconds
! ======================================================================
  subroutine second(tcpu)
    implicit none
!    real*4 etime
!    external etime
!    dimension tarray(2)
!    tcpu = etime(tarray) 
    real :: tcpu
!    integer*8 nticks, tickspersec
    integer :: nticks, tickspersec
    intrinsic system_clock
    call system_clock(count=nticks, count_rate=tickspersec)
    tcpu = 0.
    if(tickspersec .gt. 0) tcpu = nticks/(1.0*tickspersec)
    return
  end subroutine second


! ======================================================================
! calc_chi_error
! --------------
!
! calculates the error in the equation Laplacian(chi) = com
! ======================================================================
subroutine calc_chi_error(chierror)

  use arrays
  use t_data

  implicit none

  include "mpif.h"

  real, intent(out) :: chierror

  integer :: ier, itri, j, jone, numelms
  real :: chierror_local
  real :: fintl(-6:maxi,-6:maxi), d2term(18)

  call numfac(numelms)

  chierror_local = 0

  do itri=1,numelms
     call calcfint(fintl, maxi, atri(itri),btri(itri), ctri(itri))
     call calcd2term(itri, d2term, fintl)
     do j=1,18
        jone = isval1(itri,j)
        chierror_local = chierror_local + d2term(j)*com(jone)
     enddo
  enddo                  ! loop over itri
  call mpi_allreduce(chierror_local, chierror, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, ier)

end subroutine calc_chi_error


! ======================================================================
! calculate scalars
! -----------------
!
! calculates energy, energy fluxes, and other scalars
! ======================================================================
subroutine calculate_scalars()

  use basic
  use t_data
  use arrays
  use nintegrate_mod

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none
 
  integer :: itri, numelms, i, ione, ndof, def_fields
  real :: x, z, xmin, zmin, hypf, hypi, hypv, hypc, hypp, dbf, factor

  double precision, dimension(3)  :: cogcoords

  ptoto = ptot

  ekino = ekin
  emago = emag
  ekindo = ekind
  emagdo = emagd
!
  ekinpo = ekinp
  emagpo = emagp
  ekinpdo = ekinpd
  emagpdo = emagpd
!
  ekinto = ekint
  emagto = emagt
  ekintdo = ekintd
  emagtdo = emagtd
!
  ekinpho = ekinph
  ekintho = ekinth
  emagpho = emagph
  emagtho = emagth
!
  ekin3o = ekin3 
  ekin3do = ekin3d
  ekin3ho = ekin3h 
  emag3o = emag3
  emag3do = emag3d
  emag3ho = emag3h
!  
  call reset_scalars()

  hypf = hyper *deex**2
  hypi = hyperi*deex**2
  hypv = hyperv*deex**2
  hypc = hyperc*deex**2
  hypp = hyperp*deex**2

  ! Specify which fields need to be calculated
  def_fields = FIELD_PSI + FIELD_PHI + FIELD_J + FIELD_ETA
  if(numvar.ge.2) def_fields = def_fields + FIELD_I + FIELD_V
  if(numvar.ge.3) then
     def_fields = def_fields + FIELD_CHI + &
          FIELD_PE + FIELD_P + FIELD_KAP
     if(kappar.ne.0) def_fields = def_fields + FIELD_B2I
  endif
  if(idens.eq.1) def_fields = def_fields + FIELD_N + FIELD_NI

  if(hypc.ne.0.) then 
     def_fields = def_fields + FIELD_VOR
     if(numvar.ge.3) def_fields = def_fields + FIELD_COM
  end if

  tm79 = 0.
  tm79(:,OP_1) = 1.

  call getmincoord(xmin,zmin)
  
  ! Define RHS of equation
  call numfac(numelms)
  do itri=1,numelms

     if(imask.eq.1) then
        call cogfac(itri,cogcoords)
        x = cogcoords(1)-xmin
        z = cogcoords(2)-zmin
        call mask(x,z,factor)
     else
        factor = 1.
     endif
     dbf = db*factor

     call define_fields_79(itri, def_fields)

     do i=1,79 
        if(pst79(i,OP_1).ge.psilim) then
           temp79a(i) = 1.
        else
           temp79a(i) = 0.
        endif
     end do

     ! Calculate Scalars
     ! ~~~~~~~~~~~~~~~~~
     !  (extra factor of 1/r comes from delta function in toroidal coordinate)
     area   = area   + int1( ri_79,                       weight_79,79)
     parea  = parea  + int2( ri_79,               temp79a,weight_79,79)
     totcur = totcur - int2(ri2_79,pst79(:,OP_GS),        weight_79,79)
     pcur   = pcur   - int3(ri2_79,pst79(:,OP_GS),temp79a,weight_79,79)
     tvor   = tvor   - int2(ri2_79,pht79(:,OP_GS),        weight_79,79)
     pvor   = pvor   - int3(ri2_79,pht79(:,OP_GS),temp79a,weight_79,79)
     if(numvar.ge.2) then
        tflux = tflux+ int2(ri2_79,bzt79(:,OP_1 ),        weight_79,79)
        pflux = pflux+ int3(ri2_79,bzt79(:,OP_1 ),temp79a,weight_79,79)
     endif
     if(idens.eq.1) then
        totden = totden + int2( ri_79, nt79(:,OP_1),        weight_79,79)
        pden   = pden   + int3( ri_79, nt79(:,OP_1),temp79a,weight_79,79)
     endif
     if(numvar.ge.2) then
        if(idens.eq.0) then
           tmom = tmom &
                + int2(ri_79,vzt79(:,OP_1),                     weight_79,79)
           pmom = pmom &
                + int3(ri_79,vzt79(:,OP_1),temp79a,             weight_79,79)
        else
           tmom = tmom &
                + int3(ri_79,vzt79(:,OP_1),nt79(:,OP_1),        weight_79,79)
           pmom = pmom &
                + int4(ri_79,vzt79(:,OP_1),nt79(:,OP_1),temp79a,weight_79,79)
        endif
     endif


     ! Calculate energy
     ! ~~~~~~~~~~~~~~~~
     ekinp  = ekinp  + energy_kp ()
     ekinpd = ekinpd + energy_kpd()
     ekinph = ekinph + energy_kph(hypc)

     emagp  = emagp  + energy_mp ()
     emagpd = emagpd + energy_mpd()
     emagpd = emagpd - qpsipsieta(tm79,hypf)

     if(numvar.ge.2) then       
        ekint  = ekint  + energy_kt ()
        ekintd = ekintd + energy_ktd()
        ekinth = ekinth + energy_kth(hypv)

        emagt  = emagt  + energy_mt ()
        emagtd = emagtd + energy_mtd()
        emagth = emagth - qbbeta(tm79,hypi)
     endif

     if(numvar.ge.3) then
        ekin3  = ekin3  + energy_k3 ()
        ekin3d = ekin3d + energy_k3d()
        ekin3h = ekin3h + energy_k3h(hypc)
                   
        emag3 = emag3 + energy_p()
     endif

     ! Calculate fluxes
     ! ~~~~~~~~~~~~~~~~
     efluxd = efluxd + flux_diffusive()
     efluxk = efluxk + flux_ke()
     efluxp = efluxp + flux_pressure(dbf)
     efluxs = efluxs + flux_poynting(dbf)
     efluxt = efluxt + flux_heat()

     epotg = epotg + grav_pot()
  end do

  if(idens.eq.0) totden = area

  call distribute_scalars

  ekin = ekinp + ekint + ekin3
  emag = emagp + emagt + emag3
  ekind = ekinpd + ekintd + ekin3d
  emagd = emagpd + emagtd + emag3d

  ! sum all fluxes to get total energy lost through boundary
  ! (remove poynting flux due to conducting boundary)
  ptot = ptot + (efluxd + efluxk + efluxp + 0.*efluxs + efluxt + epotg &
       - vloop*totcur/(2.*pi))*dt
  if(numvar.lt.3) ptot = ptot + (ekind + emagd)*dt

  ! total energy, including energy lost through boundary flux and
  ! internal dissipation
  etot = ekin + emag - ptoto

  if(myrank.eq.0 .and. iprint.ge.1) then 
     print *, "Total energy = ", etot
     print *, "Total energy lost = ", ptot
 
     print *, "Scalars:"
     print *, "  Area = ", area
     print *, "  Toroidal current = ", totcur
     print *, "  Total electrons = ", totden
  endif

end subroutine calculate_scalars


!=====================================================
! magaxis
!
! locates the magnetic axis and the value of psi there
!=====================================================
subroutine magaxis(xguess,zguess,phin,numvari,psimin)
  use basic
  use t_data
  use nintegrate_mod

  implicit none

  include 'mpif.h'

  real, intent(inout) :: xguess, zguess
  real, intent(in), dimension(*) :: phin
  integer, intent(in) :: numvari
  real, intent(out) :: psimin

  integer, parameter :: iterations = 5

  integer :: itri, inews
  integer :: i, ier
  real :: x1, z1, x, z, theta, b, co, sn, si, eta
  real :: sum, sum1, sum2, sum3, sum4, sum5
  real :: term1, term2, term3, term4, term5
  real :: pt, pt1, pt2, p11, p22, p12
  real :: xnew, znew, denom, sinew, etanew
  real :: alx, alz
  real, dimension(20) :: avector
  real, dimension(3) :: temp1, temp2


  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " magaxis: guess=", xguess, zguess

  call getboundingboxsize(alx, alz)

  x = xguess
  z = zguess
  
  do inews=1, iterations

     call whattri(x,z,itri,x1,z1)

     ! calculate position of minimum
     if(itri.gt.0) then
        call calcavector(itri, phin, 1, numvari, avector)
         
        ! calculate local coordinates
        theta = ttri(itri)
        b = btri(itri)
        co = cos(theta)
        sn = sin(theta)
        si  = (x-x1)*co + (z-z1)*sn - b
        eta =-(x-x1)*sn + (z-z1)*co
  
        ! evaluate the polynomial and second derivative
        sum = 0.
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        sum5 = 0.
        do i=1,20
           sum = sum + avector(i)*si**mi(i)*eta**ni(i)
           term1 = 0.
           if(mi(i).ge.1) term1 = mi(i)*si**(mi(i)-1)*eta**ni(i)
           term2 = 0.
           if(ni(i).ge.1) term2 = ni(i)*si**mi(i)*eta**(ni(i)-1)
           term3 = 0.
           if(mi(i).ge.2) term3 = mi(i)*(mi(i)-1)*si**(mi(i)-2)*eta**ni(i)
           term4 = 0.
           if(ni(i).ge.2) term4 = ni(i)*(ni(i)-1)*si**mi(i)*eta**(ni(i)-2)
           term5 = 0.
           if(ni(i)*mi(i) .ge. 1)                                          &
                term5 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
           
           sum1 = sum1 + avector(i)*term1
           sum2 = sum2 + avector(i)*term2
           sum3 = sum3 + avector(i)*term3
           sum4 = sum4 + avector(i)*term4
           sum5 = sum5 + avector(i)*term5
        enddo
        pt  = sum
        pt1 = sum1
        pt2 = sum2
        p11 = sum3
        p22 = sum4
        p12 = sum5

        denom = p22*p11 - p12**2
        sinew = si -  ( p22*pt1 - p12*pt2)/denom
        etanew= eta - (-p12*pt1 + p11*pt2)/denom

        xnew = x1 + co*(b+sinew) - sn*etanew
        znew = z1 + sn*(b+sinew) + co*etanew
     else
        xnew = 0.
        znew = 0.
        pt   = 0.
     endif  ! on itri.gt.0
     
     ! communicate new minimum to all processors
     if(maxrank.gt.1) then
        temp1(1) = xnew
        temp1(2) = znew
        temp1(3) = pt
        call mpi_allreduce(temp1, temp2, 3, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
        xnew = temp2(1)
        znew = temp2(2)
        pt   = temp2(3)
     endif

     ! check to see whether the new minimum is outside the simulation domain
     if(xnew .lt. 0 .or. xnew.gt.alx .or. &
          znew .lt. 0 .or. znew.gt.alz .or. &
          xnew.ne.xnew .or. znew.ne.znew) then
        ! if not within the domain, safestop.

        write(*,3333) inews,x,z,xnew,znew
3333    format("magaxis: new minimum outside domain. ",i3,1p4e12.4)
        call safestop(27)       
     else
        x = xnew
        z = znew
     endif
  end do

  xguess = x
  zguess = z
  psimin = pt

  if(myrank.eq.0 .and. iprint.gt.0) print *, " magaxis: minimum at ", xguess, zguess
  
end subroutine magaxis



!=====================================================
! lcfs
!
! locates the magnetic axis and the value of psi there
!=====================================================
subroutine lcfs(phin, numvari)
  use basic
  use t_data
  use nintegrate_mod

  implicit none

  include 'mpif.h'

  real, intent(in), dimension(*) :: phin
  integer, intent(in) :: numvari

  integer :: itri
  real :: ajlim

  itri = 0.
  call evaluate(xlim-xzero,zlim-zzero,psilim,ajlim,phin,1,numvari,itri)

!!$  integer :: i, ifirst, ierr, numnodes, ibegin, iendplusone
!!$  integer :: izone, izonedim, ibottom, itop, ileft, iright
!!$  real :: numelms, normalderiv
!!$  real, dimension(2) :: vin, vout
!!$
!!$  call numnod(numnodes)
!!$
!!$  call getmodeltags(ibottom, iright, itop, ileft)
!!$
!!$  ifirst = 1
!!$
!!$  do i=1, numnodes
!!$     call zonenod(i,izone,izonedim)
!!$
!!$     if(izonedim.ne.1) cycle
!!$
!!$     call entdofs(numvari, i, 0, ibegin, iendplusone)
!!$
!!$     if(izone.eq.iright .and. iper.eq.0) then
!!$        normalderiv = phin(ibegin+1)
!!$     else if(izone.eq.ileft .and. iper.eq.0) then
!!$        normalderiv = -phin(ibegin+1)
!!$     else if(izone.eq.itop .and. jper.eq.0) then
!!$        normalderiv = phin(ibegin+2)
!!$     else if(izone.eq.ibottom .and. jper.eq.0) then
!!$        normalderiv = -phin(ibegin+2)
!!$     else
!!$        cycle
!!$     end if
!!$
!!$     if(normalderiv .lt. 0 .and. &
!!$          (phin(ibegin).gt.psilim .or. ifirst.eq.1)) then
!!$        psilim = phin(ibegin)
!!$        ifirst = 0
!!$     endif
!!$  end do
!!$
!!$  if(maxrank.gt.1) then
!!$     vin(1) = psilim
!!$     call MPI_allreduce(vin, vout, 1, MPI_2DOUBLE_PRECISION, &
!!$          MPI_MAXLOC, MPI_COMM_WORLD, ierr)
!!$     psilim = vout(1)
!!$  endif
  
end subroutine lcfs


end module diagnostics

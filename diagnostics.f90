module diagnostics

  implicit none

  ! scalar diagnostics
  real :: tflux0, totcur0
  real :: tflux, area, totcur, totden, tmom, tvor

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

    integer, parameter :: num_scalars = 30
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


end module diagnostics


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
  use diagnostics

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

     ! Calculate Scalars
     ! ~~~~~~~~~~~~~~~~~
     !  (extra factor of 1/r comes from delta function in toroidal coordinate)
     area   = area   + int1( ri_79,               weight_79,79)
     totcur = totcur - int2(ri2_79,pst79(:,OP_GS),weight_79,79)
     tvor   = tvor   - int2(ri2_79,pht79(:,OP_GS),weight_79,79)
     if(numvar.ge.2) tflux = tflux  + int2(ri2_79,bzt79(:,OP_1),weight_79,79)
     if(idens.eq.1) totden = totden + int2( ri_79, nt79(:,OP_1),weight_79,79)
     if(numvar.ge.2) then
        if(idens.eq.0) then
           tmom = tmom + int2(ri_79,vzt79(:,OP_1),             weight_79,79)
        else
           tmom = tmom + int3(ri_79,vzt79(:,OP_1),nt79(:,OP_1),weight_79,79)
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

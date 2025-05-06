#ifdef USECOMPLEX
#define CONJUGATE(x) conjg(x)
#else
#define CONJUGATE(x) x
#endif


module diagnostics

  implicit none

  real :: tflux0

  integer :: itri_magaxis = 0
  integer :: itri_te_max  = 0
  integer :: itri_te_max2 = 0

  ! scalars integrated over entire computational domain
  real :: tflux, area, volume, totcur, wallcur, totden, tmom, tvor, bwb2, totne, pinj, totkprad, totkprad0
  real :: totrad, linerad, bremrad, ionrad, reckrad, recprad
  real :: w_pe   ! electron thermal energy
  real :: w_m    ! totoidal poloidal magnetic energy inside plasma
  real :: w_p    ! totoidal poloidal magnetic energy inside plasma
  real :: totre  ! total number of runaway electrons
  real :: helicity ! total helicity
  
  ! wall forces in R, phi, and Z directions
  real :: wall_force_n0_x, wall_force_n0_y, wall_force_n0_z
  real :: wall_force_n0_x_halo, wall_force_n0_z_halo

  ! wall forces in x, y, and z directions
  ! note: x = R cos(phi), y = R sin(phi)
  real :: wall_force_n1_x, wall_force_n1_y, wall_force_n1_z

  ! scalars integrated within lcfs
  real :: pflux, parea, pcur, pcur_co, pcur_sn, pden, pmom, pvol, m_iz, m_iz_co, m_iz_sn, volpd
  real :: jbs
  real :: chierror, psi0

  ! density diagnostics
  real :: nfluxd, nfluxv, nsource

  ! energy diagnostics
  real :: ekin, emag, ekind, emagd, ekino, emago, ekindo, emagdo,      &
       ekint,emagt,ekintd,emagtd,ekinto,emagto,ekintdo,emagtdo,        &
       ekinp,emagp,ekinpd,emagpd,ekinpo,emagpo,ekinpdo,emagpdo,        &
       ekinph,ekinth,emagph,emagth,ekinpho,ekintho,emagpho,emagtho,    &
       ekin3,ekin3d,ekin3h,emag3,ekin3o,ekin3do,ekin3ho,emag3o,        &
       emag3h,emag3d,emag3ho,emag3do, avep,emagpv,emagtv,emagpc,emagtc
  real :: efluxp,efluxk,efluxs,efluxt,epotg,etot,ptot,eerr,ptoto

  ! in subroutine calculate_ke()
  real, allocatable:: keharmonic(:)
  integer :: NMAX

  ! in subroutine calculate_keB()
  real, allocatable:: bharmonic(:)
  integer :: BNMAX

  ! momentum diagnostics
  real :: tau_em, tau_sol, tau_com, tau_visc, tau_gyro, tau_parvisc

  logical :: ifirsttime = .true.

  ! xray diagnostic
  real :: xray_signal

  ! timing diagnostics
  real :: t_ludefall, t_sources, t_smoother, t_aux, t_onestep
  real :: t_solve_v, t_solve_n, t_solve_p, t_solve_b, t_mvm
  real :: t_output_cgm, t_output_hdf5, t_output_reset
  real :: t_gs, t_kprad
#ifdef USEPARTICLES
  real :: t_particle
#endif

  integer, parameter :: imag_probes_max = 100
  integer :: imag_probes
  real, dimension(imag_probes_max) :: mag_probe_x, mag_probe_phi, mag_probe_z
  real, dimension(imag_probes_max) :: mag_probe_nx, mag_probe_nphi, mag_probe_nz
  real, dimension(imag_probes_max) :: mag_probe_val
  integer, dimension(imag_probes_max) :: mag_probe_itri

  integer, parameter :: iflux_loops_max = 100
  integer :: iflux_loops
  real, dimension(iflux_loops_max) :: flux_loop_x, flux_loop_z
  real, dimension(iflux_loops_max) :: flux_loop_val
  integer, dimension(iflux_loops_max) :: flux_loop_itri

  vectype, dimension(MAX_PTS) :: rhop79
  vectype, dimension(MAX_PTS) :: Lorentz_pel

  integer, parameter :: REGION_PLASMA = 0
  integer, parameter :: REGION_SOL = 1
  integer, parameter :: REGION_PF  = 2

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
    t_kprad = 0.
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
    integer, parameter :: num_scalars = 14
    real, dimension(num_scalars) :: vin, vout

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
    vin(14) = t_kprad
    call MPI_ALLREDUCE(vin, vout, num_scalars, MPI_DOUBLE_PRECISION, &
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
    t_kprad         = vout(14)
    
  end subroutine distribute_timings


  ! ======================================================================
  ! reset_scalars
  ! -------------
  ! 
  ! resets diagnostic energy and scalar quantities to zero
  ! ======================================================================
  subroutine reset_scalars()
    use pellet

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
    efluxp = 0.
    efluxk = 0.
    efluxs = 0.
    efluxt = 0.
    epotg = 0.
    emagpc = 0.
    emagtc = 0.
    emagpv = 0.
    emagtv = 0.

    ptot = 0.

    area = 0.
    volume = 0.
    totcur = 0.
    wallcur = 0.
    tflux = 0.
    totden = 0.
    totne = 0.
    totrad = 0.
    linerad = 0.
    bremrad = 0.
    ionrad = 0.
    reckrad = 0.
    recprad = 0.
    tmom = 0.
    tvor = 0.
    parea = 0.
    pcur = 0.
    pcur_co = 0.
    pcur_sn = 0.
    m_iz = 0.
    m_iz_co = 0.
    m_iz_sn = 0.
    pflux = 0.
    pden = 0.
    pmom = 0.
    pvol = 0.
    pinj = 0.
    totkprad = 0.
    totkprad0 = 0.
    jbs = 0.
    volpd = 0.

    tau_em = 0.
    tau_sol = 0.
    tau_com = 0.
    tau_visc = 0.
    tau_gyro = 0.
    tau_parvisc = 0.

    nfluxd = 0.
    nfluxv = 0.
    nsource = 0.
    if(ipellet_abl.gt.0) then
       nsource_pel = 0. ! this is an array
       temp_pel = 0.    ! this is an array
    end if

    bwb2 = 0.

    xray_signal = 0.

    psi0 = 0.

    totre = 0.
    w_pe = 0.
    w_m = 0.
    w_p = 0.

    wall_force_n0_x = 0.
    wall_force_n0_y = 0.
    wall_force_n0_z = 0.
    wall_force_n1_x = 0.
    wall_force_n1_y = 0.
    wall_force_n1_z = 0.
    wall_force_n0_x_halo = 0.
    wall_force_n0_z_halo = 0.

    helicity = 0.

  end subroutine reset_scalars


  ! ======================================================================
  ! distribute_scalars
  ! 
  ! distributes diagnostic energy and scalar quantities
  ! ======================================================================
  subroutine distribute_scalars()    
    use basic
    use pellet

    implicit none

    include 'mpif.h'

    integer, parameter :: num_scalars = 82
    integer :: ier
    double precision, dimension(num_scalars) :: temp, temp2
    double precision, allocatable  :: ptemp(:)

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
       temp(19) = efluxp
       temp(20) = efluxk
       temp(21) = efluxs
       temp(22) = efluxt
       temp(23) = epotg
       temp(24) = area
       temp(25) = totcur
       temp(26) = totden
       temp(27) = tflux
       temp(28) = tmom
       temp(29) = tvor
       temp(30) = parea
       temp(31) = pcur
       temp(32) = pflux
       temp(33) = pden
       temp(34) = pmom
       temp(35) = pvol
       temp(36) = nfluxd
       temp(37) = nfluxv
       temp(38) = nsource
       temp(39) = tau_em
       temp(40) = tau_sol
       temp(41) = tau_com
       temp(42) = tau_visc
       temp(43) = tau_gyro
       temp(44) = tau_parvisc
       temp(45) = bwb2
       temp(46) = volume
       temp(47) = xray_signal
       temp(48) = wallcur
       temp(49) = totrad         
       temp(50) = linerad        
       temp(51) = bremrad        
       temp(52) = ionrad         
       temp(53) = reckrad        
       temp(54) = recprad        
       temp(55) = totre          
       temp(56) = m_iz           
       temp(57) = wall_force_n0_x
       temp(58) = wall_force_n0_y
       temp(59) = wall_force_n0_z
       temp(60) = wall_force_n1_x
       temp(61) = wall_force_n1_y
       temp(62) = wall_force_n1_z
       temp(63) = totne          
       temp(64) = w_pe           
       temp(65) = pcur_co        
       temp(66) = pcur_sn        
       temp(67) = m_iz_co        
       temp(68) = m_iz_sn        
       temp(69) = w_m
       temp(70) = w_p
       temp(71) = wall_force_n0_x_halo
       temp(72) = wall_force_n0_z_halo
       temp(73) = helicity
       temp(74) = pinj
       temp(75) = totkprad
       temp(76) = totkprad0
       temp(77) = emagpc
       temp(78) = emagtc
       temp(79) = emagpv
       temp(80) = emagtv
       temp(81) = jbs
       temp(82) = volpd

       !checked that this should be MPI_DOUBLE_PRECISION
       call mpi_allreduce(temp, temp2, num_scalars, MPI_DOUBLE_PRECISION,  &
            MPI_SUM, MPI_COMM_WORLD, ier) 
         
       ekinp           = temp2( 1)
       emagp           = temp2( 2)
       ekinpd          = temp2( 3)
       emagpd          = temp2( 4)      
       ekint           = temp2( 5)
       emagt           = temp2( 6)
       ekintd          = temp2( 7)
       emagtd          = temp2( 8)
       ekinph          = temp2( 9)
       ekinth          = temp2(10)
       emagph          = temp2(11)
       emagth          = temp2(12)
       ekin3           = temp2(13)
       ekin3d          = temp2(14)
       ekin3h          = temp2(15)
       emag3           = temp2(16)
       emag3d          = temp2(17)
       emag3h          = temp2(18)
       efluxp          = temp2(19)
       efluxk          = temp2(20)
       efluxs          = temp2(21)
       efluxt          = temp2(22)
       epotg           = temp2(23)
       area            = temp2(24)
       totcur          = temp2(25)
       totden          = temp2(26)
       tflux           = temp2(27)
       tmom            = temp2(28)
       tvor            = temp2(29)
       parea           = temp2(30)
       pcur            = temp2(31)
       pflux           = temp2(32)
       pden            = temp2(33)
       pmom            = temp2(34)
       pvol            = temp2(35)
       nfluxd          = temp2(36)
       nfluxv          = temp2(37)
       nsource         = temp2(38)
       tau_em          = temp2(39)
       tau_sol         = temp2(40)
       tau_com         = temp2(41)
       tau_visc        = temp2(42)
       tau_gyro        = temp2(43)
       tau_parvisc     = temp2(44)
       bwb2            = temp2(45)
       volume          = temp2(46)
       xray_signal     = temp2(47)
       wallcur         = temp2(48)
       totrad          = temp2(49)
       linerad         = temp2(50)
       bremrad         = temp2(51)
       ionrad          = temp2(52)
       reckrad         = temp2(53)
       recprad         = temp2(54)
       totre           = temp2(55)
       m_iz            = temp2(56)
       wall_force_n0_x = temp2(57)
       wall_force_n0_y = temp2(58)
       wall_force_n0_z = temp2(59)
       wall_force_n1_x = temp2(60)
       wall_force_n1_y = temp2(61)
       wall_force_n1_z = temp2(62)
       totne           = temp2(63)
       w_pe            = temp2(64)
       pcur_co         = temp2(65)
       pcur_sn         = temp2(66)
       m_iz_co         = temp2(67)
       m_iz_sn         = temp2(68)
       w_m             = temp2(69)
       w_p             = temp2(70)
       wall_force_n0_x_halo = temp2(71)
       wall_force_n0_z_halo = temp2(72)
       helicity        = temp2(73)
       pinj            = temp2(74)
       totkprad        = temp2(75)
       totkprad0       = temp2(76)
       emagpc          = temp2(77)
       emagtc          = temp2(78)
       emagpv          = temp2(79)
       emagtv          = temp2(80)
       jbs             = temp2(81)
       volpd           = temp2(82)

       if(ipellet_abl.gt.0) then
          allocate(ptemp(npellets))

          ptemp = nsource_pel
          call mpi_allreduce(ptemp, nsource_pel, npellets, MPI_DOUBLE_PRECISION,  &
                             MPI_SUM, MPI_COMM_WORLD, ier)
          ptemp = temp_pel
          call mpi_allreduce(ptemp, temp_pel, npellets, MPI_DOUBLE_PRECISION,  &
                            MPI_SUM, MPI_COMM_WORLD, ier)

          deallocate(ptemp)
       end if

    endif

  end subroutine distribute_scalars


!============================================================
! evaluate
! ~~~~~~~~
! calculates the value ans of field dum at global coordinates
! (x,z).  itri is the element containing (x,z).  (If this
! element does not reside on this process, itri=-1).
!============================================================
subroutine evaluate(x,phi,z,ans,fin,itri,ierr)
  
  use mesh_mod
  use basic
  use m3dc1_nint
  use field

  implicit none

  include 'mpif.h'

  integer, intent(inout) :: itri
  real, intent(in) :: x, phi, z
  type(field_type), intent(in) :: fin
  integer, intent(out) :: ierr ! = 0 on success

  real, intent(out), dimension(OP_NUM) :: ans

  type(element_data) :: d
  integer :: nodeids(nodes_per_element), ier
  real :: x1, phi1, z1
  vectype, dimension(OP_NUM) :: temp1, temp2
  integer :: hasval, tothasval

  ! evaluate the solution to get the value [ans] at one point (x,z)

  ! first find out what triangle x,z is in.  whattri
  ! returns itri, x1, and z1 with x1 and z1 being
  ! the coordinates of the first node/vertex

  if(itri.eq.0) then
     call whattri(x,phi,z,itri,x1,z1)
  else if(itri.gt.0) then
     call get_element_nodes(itri,nodeids)
     call get_node_pos(nodeids(1), x1, phi1, z1)
  endif

  temp1 = 0.
  temp2 = 0.
  ans = 0.

  ! if this process contains the point, evaluate the field at that point.
  if(itri.gt.0) then

     call get_element_data(itri, d)

     ! calculate local coordinates
     call global_to_local(d, x, phi, z, xi_79(1), zi_79(1), eta_79(1))

     ! calculate the inverse radius
     if(itor.eq.1) then
        ri_79(1) = 1./x
     else
        ri_79(1) = 1.
     endif

     call precalculate_terms(xi_79, zi_79, eta_79, d%co, d%sn, ri_79, 1)
     call define_basis(itri)

     ! calculate the value of the function
     call eval_ops(itri, fin, tm79, rfac)
     
     temp1 = tm79(1,:)

     hasval = 1
  else
     hasval = 0
  endif
     

  ! Distribute the result if necessary
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(maxrank.gt.1) then
     ! Determine number of processes whose domains contain this point
     call mpi_allreduce(hasval, tothasval, 1, MPI_INTEGER, MPI_SUM, &
          MPI_COMM_WORLD, ier)
  else
     tothasval = hasval
  end if

  if(tothasval.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          write(*,'(A,3f12.4)') 'Point not found in domain: ', x, phi, z
     ierr = 1
     return
  end if

  if(maxrank.gt.1) then
     ! Find the average value at this point over all processes containing
     ! the point.  (Each value should be identical.)
#ifdef USECOMPLEX
     call mpi_allreduce(temp1, temp2, OP_NUM, MPI_DOUBLE_COMPLEX, MPI_SUM, &
          MPI_COMM_WORLD, ier)
#else
     call mpi_allreduce(temp1, temp2, OP_NUM, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_WORLD, ier)
#endif
  else
     temp2 = temp1
  endif

#ifdef USECOMPLEX
  ans = real(temp2*exp(rfac*phi))/tothasval
#else
  ans = temp2/tothasval
#endif

  ierr = 0
end subroutine evaluate


  ! ======================================================================
  ! reconnected_flux
  !
  ! calculates the reconnected flux assuming an x-point at the center of
  ! the computational domain, and an o-point at the center of the
  ! left boundary
  ! ======================================================================
  real function reconnected_flux()
    
    use basic
    use arrays
    use mesh_mod
    use m3dc1_nint
  
    implicit none

    real :: alx, alz, xrel, zrel
    real, dimension(OP_NUM) :: temp, temp2
    integer :: itri, ierr

    call get_bounding_box_size(alx,alz)


    itri = 0
    xrel = xzero
    zrel = alz/2. + zzero
    call evaluate(xrel,0.,zrel,temp,psi_field(1),itri,ierr)

    itri = 0
    xrel = alx/2. + xzero
    zrel = alz/2. + zzero
    call evaluate(xrel,0.,zrel,temp2,psi_field(1),itri,ierr)

    reconnected_flux = 0.5*(temp2(OP_1)-temp(OP_1))

    return
  end function reconnected_flux

! ======================================================================
! second
!
! returns the current time in seconds
! ======================================================================

  subroutine second(tcpu)
    implicit none
    real :: tcpu
  intrinsic cpu_time
    call cpu_time(tcpu)
    return
  end subroutine second

  subroutine second_new(time)
    implicit none
    real, intent(out) :: time
    integer :: count, rate
    call system_clock(count,count_rate=rate)
    time=real(count)/rate
    return
  end subroutine second_new

!   Added 1/1/2016 to get consistency between 2D,3D,Cyl,Tor
subroutine tpi_factors(tpifac,tpirzero)
  use basic
  use math
  implicit none
  real, intent(out) :: tpifac, tpirzero
  if(nplanes.eq.1) then
     if(itor.eq.1) then
        tpifac = 1.
        tpirzero = 1.
     else
        tpifac = 1./rzero
        tpirzero = 1.
     endif
  else
     if(itor.eq.1) then
        tpifac = twopi
        tpirzero = twopi
     else
        tpifac = twopi
        tpirzero = twopi*rzero
     endif
     if(ifull_torus.eq.0) then
        tpifac = twopi/nperiods
        tpirzero = tpirzero/nperiods
     endif
  endif
end subroutine tpi_factors

! ======================================================================
! calculate scalars
! -----------------
!
! calculates energy, energy fluxes, and other scalars
! ======================================================================
subroutine calculate_scalars()

  use basic
  use mesh_mod
  use arrays
  use m3dc1_nint
  use newvar_mod
  use sparse
  use metricterms_new
  use boundary_conditions
  use math
  use gyroviscosity
  use pellet
  use kprad_m3dc1
  use bootstrap

  implicit none
 
  include 'mpif.h'

  integer :: itri, numelms, def_fields, ier
  integer :: is_edge(3)  ! is inode on boundary
  real :: n(2,3),tpifac,tpirzero, t0
  integer :: iedge, idim(3), izone, izonedim, izone_ind, i, j
  real, dimension(OP_NUM) :: dum1
  vectype, dimension(MAX_PTS) :: mr
  vectype, dimension(MAX_PTS) :: co, sn

  integer :: ip

  call tpi_factors(tpifac,tpirzero)

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

  ! Specify which fields need to be calculated
  if(ike_only.eq.1) then
     def_fields = FIELD_PHI + FIELD_N
     if(numvar.ge.2) def_fields = def_fields + FIELD_V
     if(numvar.ge.3) def_fields = def_fields + FIELD_CHI
  else
     def_fields = FIELD_PSI + FIELD_PHI + FIELD_ETA + FIELD_MU &
          + FIELD_N + FIELD_NI + FIELD_SIG
     if(numvar.ge.2) def_fields = def_fields + FIELD_I + FIELD_V
     if(numvar.ge.3) then
        def_fields = def_fields + FIELD_CHI 
     endif

     if(gyro.eq.1 .or. amupar.ne.0) then
        def_fields = def_fields + FIELD_B2I
     endif

     if(numvar.ge.3 .or. ipres.eq.1) then
        def_fields = def_fields + FIELD_P + FIELD_KAP + FIELD_TE + FIELD_TI + FIELD_Q
        if(hyper.eq.0.) def_fields = def_fields + FIELD_J
        if(rad_source) def_fields = def_fields + FIELD_RAD
     end if

     if(irunaway.gt.0) then 
        def_fields = def_fields + FIELD_RE
     end if

     if(ibootstrap.gt.0) then
        def_fields = def_fields + FIELD_JBS
     endif
  endif

!  tm79 = 0.
!  tm79(:,OP_1) = 1.
  mr = 0.

  call finalize(field0_vec)
  call finalize(field_vec)

  numelms = local_elements()
  
  if(ipellet.ne.0) call calculate_Lor_vol

  ! BCL Warning: nsource_pel and temp_pel are now vectors
  !              this compiles, but may break at runtime for OpenMP (OMP=1)
!!$OMP PARALLEL DO PRIVATE(mr,dum1,ier,is_edge,n,iedge,idim,izone,izonedim,izone_ind,i) &
!!$OMP& REDUCTION(+:ekinp,ekinpd,ekinph,ekint,ekintd,ekinth,ekin3,ekin3d,ekin3h,wallcur,emagp,emagpd,emagph,emagt,emagtd,emagth,emag3,area,parea,totcur,pcur,m_iz,tflux,pflux,tvor,volume,pvol,volpd,totden,pden,totrad,linerad,bremrad,ionrad,reckrad,recprad,totre,nsource,epotg,tmom,pmom,bwb2,efluxp,efluxt,efluxs,efluxk,tau_em,tau_sol,tau_com,tau_visc,tau_gyro,tau_parvisc,nfluxd,nfluxv,xray_signal,Lor_vol,nsource_pel,temp_pel,wall_force_n0_x,wall_force_n0_y,wall_force_n0_z,wall_force_n1_x,wall_force_n1_y,wall_force_n1_z,totne,w_pe,pcur_co,pcur_sn,m_iz_co,m_iz_sn,w_m,w_p,wall_force_n0_x_halo,wall_force_n0_z_halo,helicity,pinj,totkprad,totkprad0,emagpc,emagtc,emagpv,emagtv,jbs)
  do itri=1,numelms

     !call zonfac(itri, izone, izonedim)
     call m3dc1_ent_getgeomclass(2, itri-1,izonedim,izone_ind)
     izone = zone_type(izone_ind)

     call define_element_quadrature(itri, int_pts_diag, int_pts_tor)
     call define_fields(itri, def_fields, 0, 0)
     call calculate_rho(itri)
     if(gyro.eq.1) call gyro_common

#ifdef USE3D
     if(ifull_torus.eq.1) then
        co = cos(phi_79*twopi/toroidal_period)
        sn = sin(phi_79*twopi/toroidal_period)
     else
        co = cos(phi_79*twopi/(toroidal_period*nperiods))
        sn = sin(phi_79*twopi/(toroidal_period*nperiods))
     end if
#endif

     if(imulti_region.eq.1 .and. izone.eq.ZONE_CONDUCTOR) then
        wallcur = wallcur - int2(ri2_79,pst79(:,OP_GS))/tpirzero

        call jxb_r(temp79a, temp79d)
        call jxb_phi(temp79b)
        call jxb_z(temp79c, temp79e)

        wall_force_n0_x = wall_force_n0_x + int1(temp79a)*twopi/tpifac
        wall_force_n0_y = wall_force_n0_y + int1(temp79b)*twopi/tpifac
        wall_force_n0_z = wall_force_n0_z + int1(temp79c)*twopi/tpifac

        wall_force_n0_x_halo = wall_force_n0_x_halo + &
             int1(temp79d)*twopi/tpifac
        wall_force_n0_z_halo = wall_force_n0_z_halo + &
             int1(temp79e)*twopi/tpifac

#ifdef USE3D
        wall_force_n1_x = wall_force_n1_x &
             + int2(temp79a,co)*twopi/tpifac &
             - int2(temp79b,sn)*twopi/tpifac
        wall_force_n1_y = wall_force_n1_y &
             + int2(temp79a,sn)*twopi/tpifac &
             + int2(temp79b,co)*twopi/tpifac
#endif
     end if

     if(izone.eq.ZONE_CONDUCTOR) then
        emagpc  = emagpc + twopi*energy_mp()/tpifac
        emagtc  = emagtc + twopi*energy_mt()/tpifac
     end if

     if(izone.eq.ZONE_VACUUM) then
        emagpv  = emagpv + twopi*energy_mp()/tpifac
        emagtv  = emagtv + twopi*energy_mt()/tpifac
     end if

     if(izone.ne.ZONE_PLASMA) cycle

     do i=1, npoints
        if(linear.eq.1) then
           call magnetic_region(ps079(i,OP_1),ps079(i,OP_DR),ps079(i,OP_DZ), &
                x_79(i),z_79(i),j)
           if(j.eq.REGION_PLASMA) then
              mr(i) = 1.
           else
              mr(i) = 0.
           end if
        else
           call magnetic_region(pst79(i,OP_1),pst79(i,OP_DR),pst79(i,OP_DZ), &
                x_79(i),z_79(i),j)
           if(j.eq.REGION_PLASMA) then
              mr(i) = 1.
           else
              mr(i) = 0.
           end if
        endif
     end do

     ! Calculate energy
     ! ~~~~~~~~~~~~~~~~
     ekinp  = ekinp  + twopi*energy_kp ()/tpifac
     ekinpd = ekinpd + twopi*energy_kpd()/tpifac
     ekinph = ekinph + twopi*energy_kph()/tpifac

     ekint  = ekint  + twopi*energy_kt ()/tpifac
     ekintd = ekintd + twopi*energy_ktd()/tpifac
     ekinth = ekinth + twopi*energy_kth()/tpifac

     ekin3  = ekin3  + twopi*energy_k3 ()/tpifac
     ekin3d = ekin3d + twopi*energy_k3d()/tpifac
     ekin3h = ekin3h + twopi*energy_k3h()/tpifac

     if(ike_only.eq.1) cycle

     emagp  = emagp  + twopi*energy_mp ()/tpifac
     emagpd = emagpd + twopi*energy_mpd()/tpifac
!     emagph = emagph - twopi*qpsipsieta(tm79)/tpifac

     emagt  = emagt  + twopi*energy_mt ()/tpifac
     emagtd = emagtd + twopi*energy_mtd()/tpifac
!     emagth = emagth - twopi*qbbeta(tm79)/tpifac

     emag3 = emag3 + twopi*energy_p()/tpifac
     w_pe = w_pe + twopi*energy_pe()/tpifac
     w_p  = w_p +  twopi*energy_p(mr)/tpifac
     w_m  = w_m + twopi*energy_mp(mr)/tpifac

     ! Calculate Scalars
     ! ~~~~~~~~~~~~~~~~~
     ! extra factor of 1/r comes from delta function in toroidal coordinate)
     area   = area   + int1(ri_79)/tpirzero
     parea  = parea  + int2(ri_79,mr)/tpirzero

     ! toroidal current
     totcur = totcur - int2(ri2_79,pst79(:,OP_GS))/tpirzero
     pcur   = pcur   - int3(ri2_79,pst79(:,OP_GS),mr)/tpirzero
#ifdef USE3D
     pcur_co = pcur_co - int4(ri2_79,pst79(:,OP_GS),mr,co)/tpirzero * 2.
     pcur_sn = pcur_sn - int4(ri2_79,pst79(:,OP_GS),mr,sn)/tpirzero * 2.
#endif
    ! bootstrap current
    if (ibootstrap.eq.1)then
      call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
     ! Jp_BS_Phi = intx3(ri_79,bzt79(:,OP_1),temp79a) 
      
      jbs = jbs + int4(ri2_79,bzt79(:,OP_1),temp79a,mr) 
    elseif (ibootstrap.eq.2)then
      call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
      jbs = jbs + int4(ri2_79,bzt79(:,OP_1),temp79a,mr) 
    endif
     ! M_iz = int(dV Z*J)
     ! This is used for calculating the vertical "center" of the plasma current
     temp79a = z_79
     m_iz   = m_iz   - int4(ri2_79,temp79a,pst79(:,OP_GS),mr)/tpirzero
#ifdef USE3D
     m_iz_co = m_iz_co - int5(ri2_79,temp79a,pst79(:,OP_GS),mr,co)/tpirzero * 2.
     m_iz_sn = m_iz_sn - int5(ri2_79,temp79a,pst79(:,OP_GS),mr,sn)/tpirzero * 2.
#endif

     ! toroidal flux
     tflux = tflux + int2(ri2_79,bzt79(:,OP_1))/tpirzero
     pflux = pflux + int3(ri2_79,bzt79(:,OP_1),mr)/tpirzero
     
     ! enstrophy
     tvor   = tvor &
          -    (int1(pht79(:,OP_LP)) &
          + 2.*int2(ri4_79,cht79(:,OP_DZ)))/tpirzero

     ! volume
     volume = volume + twopi*int0()/tpifac
     pvol = pvol + twopi*int1(mr)/tpifac
     volpd = volpd + twopi*volume_pd(mr)/tpifac

     ! particle number
     totden = totden + twopi*int1(nt79(:,OP_1))/tpifac
     totne = totne + twopi*int1(net79(:,OP_1))/tpifac
     pden = pden + twopi*int2(nt79(:,OP_1),mr)/tpifac
      ! radiation
     totrad = totrad + twopi*int1(totrad79(:,OP_1))/tpifac
     linerad = linerad + twopi*int1(linerad79(:,OP_1))/tpifac
     bremrad = bremrad + twopi*int1(bremrad79(:,OP_1))/tpifac
     ionrad = ionrad + twopi*int1(ionrad79(:,OP_1))/tpifac
     reckrad = reckrad + twopi*int1(reckrad79(:,OP_1))/tpifac
     recprad = recprad + twopi*int1(recprad79(:,OP_1))/tpifac

     if(ikprad.ne.0) then
        call eval_ops(itri, kprad_n(0), tm79, rfac)
        t0 = twopi*int1(tm79(:,OP_1))/tpifac
        totkprad0 = totkprad0 + t0
        totkprad = totkprad + t0
        do i=1, kprad_z
           call eval_ops(itri, kprad_n(i), tm79, rfac)
           totkprad = totkprad + twopi*int1(tm79(:,OP_1))/tpifac
        end do
     end if
     
     if(irunaway.gt.0) then
        totre = totre - int2(ri_79,nre179(:,OP_1))/tpirzero
        totre = totre - int2(ri_79,nre079(:,OP_1))/tpirzero
     end if

     helicity = helicity &
          + twopi*int3(ri2_79,pstx79(:,OP_1),bztx79(:,OP_1))/tpifac
#if defined(USE3D) || defined(USECOMPLEX)
     helicity = helicity &
          + twopi*int2(bftx79(:,OP_1),pst79(:,OP_GS))/tpifac
#endif

     ! particle source
     if(idens.eq.1) then        
        nsource = nsource - twopi*int1(sig79)/tpifac
     
        ! Pellet radius and density/temperature at the pellet surface
        if(ipellet_abl.gt.0) then
           do ip=1,npellets
              if((pellet_state(ip).eq.1).and.(r_p(ip).ge.1e-8)) then
                 ! weight density/temp by pellet distribution (normalized)
                 temp79a = pellet_distribution(ip, x_79, phi_79, z_79, real(pt79(:,OP_1)), 1, izone)
                 nsource_pel(ip) = nsource_pel(ip) + twopi*int2(net79(:,OP_1),temp79a)/tpifac
                 temp79b = pet79(:,OP_1)/net79(:,OP_1)
                 if(ikprad_te_offset .gt. 0) temp79b = temp79b - eta_te_offset
                 temp_pel(ip) = temp_pel(ip) + twopi*int2(temp79b,temp79a)*p0_norm/(1.6022e-12*n0_norm*tpifac)
              else
                 nsource_pel(ip) = 0.
                 temp_pel(ip) = 0.
              end if
           end do
       endif
     endif

     ! gravitational potential energy
     epotg = epotg + grav_pot()

     ! toroidal (angular) momentum
     if(numvar.ge.2) then
        tmom = tmom &
             + twopi*int3(r2_79,vzt79(:,OP_1),rho79(:,OP_1))/tpifac
        pmom = pmom &
             + twopi*int4(r2_79,vzt79(:,OP_1),rho79(:,OP_1),mr)/tpifac
     endif

     ! injected power
     pinj = pinj + twopi*int1(q79)/tpifac

     if(amupar.ne.0.) then
        call PVS1(pht79,temp79a)

        if(numvar.ge.2) then
           call PVS2(vzt79,temp79b)
           temp79a = temp79a + temp79b
        endif

        if(numvar.ge.3) then
           call PVS3(cht79,temp79c)
           temp79a = temp79a + temp79c
        endif

        bwb2 = bwb2 + &
             3.*int3(vip79(:,OP_1),temp79a,CONJUGATE(temp79a))
     end if

     ! add surface terms
     call boundary_edge(itri, is_edge, n, idim, BOUND_FIRSTWALL)

     do iedge=1,3
        if(is_edge(iedge).eq.0) cycle

        call define_boundary_quadrature(itri, iedge, 5, 5, n, idim)
        call define_fields(itri, def_fields, 1, 0)
        call calculate_rho(itri)
        if(gyro.eq.1) call gyro_common

        ! Energy fluxes
        ! ~~~~~~~~~~~~~
        efluxp = efluxp + twopi*flux_pressure()/tpifac
        efluxt = efluxt + twopi*flux_heat()/tpifac
        efluxs = efluxs + twopi*flux_poynting()/tpifac
        efluxk = efluxk + twopi*flux_ke()/tpifac

        ! Toroidal momentum fluxes
        ! ~~~~~~~~~~~~~~~~~~~~~~~~
        if(numvar.ge.2) then
           tau_em   = tau_em   + torque_em()
           tau_sol  = tau_sol  + torque_sol()
           tau_com  = tau_com  + torque_com()
           tau_visc = tau_visc + torque_visc()
           tau_gyro = tau_gyro + torque_gyro()
           tau_parvisc = tau_parvisc + torque_parvisc()
        endif

        ! Particle fluxes
        ! ~~~~~~~~~~~~~~~
        if(idens.eq.1) then
           nfluxd = nfluxd - twopi*denm* &
                (int2(norm79(:,1),nt79(:,OP_DR)) &
                +int2(norm79(:,2),nt79(:,OP_DZ)))/tpifac

           nfluxv = nfluxv &
                + int4(r_79,nt79(:,OP_1),norm79(:,2),pht79(:,OP_DR)) &
                - int4(r_79,nt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ))
           
           if(numvar.ge.3) then
              nfluxv = nfluxv &
                   + int4(ri2_79,nt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
                   + int4(ri2_79,nt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ))
           endif

           nfluxv = nfluxv*twopi/tpifac
        end if

        ! xray signal
        temp79b = bremsstrahlung(nt79(:,OP_1), pet79(:,OP_1))
        call get_chord_mask(xray_r0, xray_phi0*pi/180., xray_z0, &
             x_79, phi_79, z_79, npoints, &
             xray_theta*pi/180., xray_sigma*pi/180., temp79a)
        xray_signal = xray_signal + int2(temp79a, temp79b)
     end do
  end do
!!$OMP END PARALLEL DO

  call distribute_scalars

  if(ipellet_abl.gt.0) call calculate_ablation

  ekin = ekinp + ekint + ekin3
  emag = emagp + emagt + emag3
  ekind = ekinpd + ekintd + ekin3d
  emagd = emagpd + emagtd + emag3d

  ! sum all fluxes to get total energy lost through boundary
  ptot = ptot + (efluxk + efluxp + efluxs + efluxt + epotg)*dt
  if(numvar.lt.3) ptot = ptot + (ekind + emagd)*dt

  ! total energy, including energy lost through boundary flux and
  ! internal dissipation
  etot = ekin + emag - ptoto
!
!   volume averaged pressure for beta calculation
    avep = (gam - 1.)*(w_p / pvol)

    ! psi on axis
    itri = 0
    call evaluate(xmag,0.,zmag,dum1,psi_field(1),itri,ier)
    psi0 = dum1(OP_1)

#ifdef USE3D
  if(ike_harmonics .gt. 0) call calculate_ke()
  if(ibh_harmonics .gt. 0) call calculate_bh()
#endif

  call evaluate_mag_probes
  call evaluate_flux_loops

  if(myrank.eq.0 .and. iprint.ge.1) then 
     print *, "Total energy = ", etot, emagp, emagt,emag3
     print *, "Total energy lost = ", ptot
 
     print *, "Scalars:"
     print *, "  Area = ", area
     print *, "  Toroidal current = ", totcur
     print *, "  Toroidal flux = ", tflux
     print *, "  Volume = ", volume
     print *, "  Total particles = ", totden
     print *, "  Total radiation = ", totrad
     print *, "  Line radiation = ", linerad
     print *, "  Bremsstrahlung radiation = ", bremrad
     print *, "  Ionization loss = ", ionrad
     print *, "  Recombination radiation (kinetic) = ", reckrad
     print *, "  Recombination radiation (potential) = ", recprad
     if(ipellet_abl.gt.0) then
        if(iprint.ge.3 .or. npellets.eq.1) then
           do ip=1,npellets
              print *, "  Pellet #", ip
              print *, "    state = ", pellet_state(ip)
              print *, "    particles injected = ", pellet_rate(ip)*dt*(n0_norm*l0_norm**3)
              print *, "    radius (in cm) = ", r_p(ip)*l0_norm
              print *, "    local electron temperature (in eV) = ", temp_pel(ip)
              print *, "    local electron density (in ne14) = ", nsource_pel(ip)
              print *, "    rpdot (in cm/s) = ", rpdot(ip)*l0_norm/t0_norm
              print *, "    Lor_vol = ", Lor_vol(ip)
              print *, "    R position: ", pellet_r(ip)*l0_norm
              print *, "    phi position: ", pellet_phi(ip)
              print *, "    Z position: ", pellet_z(ip)*l0_norm
           end do
        elseif(iprint.ge.2) then
           print *, "  Pellet diagnostic output suppressed for multi-pellet"
           print *, "    Set iprint >= 3 to output multi-pellet diagnostics"
        endif
     endif
  endif

end subroutine calculate_scalars


subroutine calculate_Lor_vol()

  use basic
  use mesh_mod
  use m3dc1_nint
  use math
  use pellet

  implicit none
 
  include 'mpif.h'

  integer :: itri, numelms, ier
  integer :: is_edge(3)  ! is inode on boundary
  real :: tpifac,tpirzero
  integer :: izone, izonedim, izone_ind
  real, allocatable :: temp(:)
  integer :: ip

  call tpi_factors(tpifac,tpirzero)

  numelms = local_elements()

  allocate(temp(npellets))
  temp = 0.

  do itri=1,numelms

     call m3dc1_ent_getgeomclass(2, itri-1,izonedim,izone_ind)
     izone = zone_type(izone_ind)
         
     if(izone.ne.1) cycle
     call define_element_quadrature(itri, int_pts_diag, int_pts_tor)
     call define_fields(itri, FIELD_P, 0, 0)

     ! perform volume integral of pellet cloud (without normalization)
     do ip=1,npellets
        temp79a  = pellet_distribution(ip, x_79, phi_79, z_79, real(pt79(:,OP_1)), 0, izone)
        temp(ip) = temp(ip) + twopi*int1(temp79a)/tpifac
     end do

  end do

  call mpi_allreduce(temp, Lor_vol, npellets, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier )

  deallocate(temp)

end subroutine calculate_Lor_vol


!======================================================================
! magnetic_region
! ~~~~~~~~~~~~~~~
! determines what magnetic region the point x, z is in
! 0: inside plasma
! 1: scrape-off layer
! 2: private flux
!======================================================================
elemental subroutine magnetic_region(psi, psix, psiz, x, z, mr, psinb)
  use basic

  implicit none

  vectype, intent(in) :: psi, psix, psiz
  real, intent(in) :: x, z 
  integer, intent(out) :: mr
  real, intent(out), optional :: psinb    ! psi_norm at the "boundary" of this region
  real :: psii, psii_null, psii_null2, dpsii, pl, rl, al

  mr = REGION_PLASMA
  if(present(psinb)) psinb = 1.

#ifdef USEST
  if (igeometry.ge.1) return 
#endif

  dpsii = psibound - psimin
  if(dpsii.eq.0.) return

  psii = (real(psi) - psimin)/dpsii
  if(psii .lt. 0.) then
     mr = REGION_PLASMA
     return
  end if

  if(xnull.ne.0) psii_null = (psinull - psimin)/dpsii
  if(xnull2.ne.0) psii_null2 = (psinull2 - psimin)/dpsii

  if(psii.gt.1. .and. &
       (xnull.eq.0. .or. psii.gt.psii_null) .and. &
       (xnull2.eq.0. .or. psii.gt.psii_null2)) then
     mr = REGION_SOL
     return
  end if

  pl = sqrt(real(psix)**2 + real(psiz)**2)
  rl = sqrt((x-xmag)**2 + (z-zmag)**2)
  if(pl.eq.0. .or. rl.eq.0.) then
     mr = REGION_PLASMA
     return
  end if

  al = (real(psix)*(x-xmag) + real(psiz)*(z-zmag))/(pl*rl)
  if(xnull.gt.0) then
     if(psii.lt.psii_null) then
        ! Make sure we're in the direction of this x-point
        if((x-xmag)*(xnull-xmag) + (z-zmag)*(znull-zmag) .gt. 0.) then
           ! Make sure either flux gradient has flipped, or we're 
           ! farther from axis than x-point
           if(al*dpsii/abs(dpsii) .lt. 0.3 .or. &
                rl.gt.sqrt((xmag-xnull)**2 + (zmag-znull)**2)) then
              if(present(psinb)) psinb = psii_null
              mr = REGION_PF
              return
           end if
        end if
     end if
  end if
  if(xnull2.gt.0) then
     if(psii.lt.psii_null2) then
        ! Make sure we're in the direction of this x-point
        if((x-xmag)*(xnull2-xmag) + (z-zmag)*(znull2-zmag) .gt. 0.) then
           ! Make sure either flux gradient has flipped, or we're 
           ! farther from axis than x-point
           if(al*dpsii/abs(dpsii) .lt. 0.3 .or. &
                rl.gt.sqrt((xmag-xnull2)**2 + (zmag-znull2)**2)) then
              if(present(psinb)) psinb = psii_null2
              mr = REGION_PF
              return
           end if
        end if
     end if
  end if

  if(psii .gt. 1.) then
     ! if Psi > 1, and wer're not in a private flux region, 
     ! we are in scrape-off layer
     mr = REGION_SOL
     return
  end if

  mr = REGION_PLASMA
end subroutine magnetic_region

subroutine reset_itris()

  implicit none

  itri_magaxis = 0
  itri_te_max  = 0
  itri_te_max2 = 0
  return

end subroutine reset_itris

!=====================================================
! magaxis
! ~~~~~~~
! locates the magnetic axis and the value of psi there
! imethod = 0 finds the local minimum of psi
! imethod = 1 finds the local zero of <psi,psi>
!=====================================================
subroutine magaxis(xguess,zguess,psi,psim,imethod,ier)
  use basic
  use mesh_mod
  use m3dc1_nint
  use field

  implicit none

  include 'mpif.h'

  real, intent(inout) :: xguess, zguess
  type(field_type), intent(in) :: psi
  integer, intent(in) :: imethod
  real, intent(out) :: psim

  integer, parameter :: iterations = 50  !  max number of Newton iterations
  real, parameter :: bfac = 0.1  !max zone fraction for movement each iteration
  real, parameter :: tol = 1e-3   ! convergence tolorance (fraction of h)

  type(element_data) :: d
  integer :: inews
  integer :: i, ier, in_domain, converged, izone
  real :: x1, z1, x, z, si, zi, eta, h
  real :: sum, sum1, sum2, sum3, sum4, sum5
  real :: term1, term2, term3, term4, term5
  real :: pt, pt1, pt2, p11, p22, p12
  real :: xnew, znew, denom, sinew, etanew
  real :: xtry, ztry, rdiff
  vectype, dimension(coeffs_per_element) :: avector
  real, dimension(6) :: temp1, temp2
  integer :: itri

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,2E12.4)') '  magaxis: guess = ', xguess, zguess

  converged = 0
  izone = 0

  x = xguess
  z = zguess
  itri = itri_magaxis
  
  newton :  do inews=1, iterations

     call whattri(x,0.,z,itri,x1,z1)

     ! calculate position of minimum
     if(itri.gt.0) then
        call get_zone(itri,izone)
        call calcavector(itri, psi, avector)
        call get_element_data(itri, d)

        ! calculate local coordinates
        call global_to_local(d, x, 0., z, si, zi, eta)

        ! calculate mesh size
        h = sqrt((d%a+d%b)*d%c)

        ! evaluate the polynomial and second derivative
        sum = 0.
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        sum5 = 0.
        do i=1, coeffs_per_tri
           sum = sum + avector(i)*si**mi(i)*eta**ni(i)
           term1 = 0.
           if(mi(i).ge.1) term1 = mi(i)*si**(mi(i)-1)*eta**ni(i)
           term2 = 0.
           if(ni(i).ge.1) term2 = ni(i)*si**mi(i)*eta**(ni(i)-1)
           term3 = 0.
           if(mi(i).ge.2) term3 = mi(i)*(mi(i)-1)*si**(mi(i)-2)*eta**ni(i)
           term4 = 0.
           if(ni(i)*mi(i) .ge. 1)                                          &
                term4 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
           term5 = 0.
           if(ni(i).ge.2) term5 = ni(i)*(ni(i)-1)*si**mi(i)*eta**(ni(i)-2)
           
           sum1 = sum1 + avector(i)*term1
           sum2 = sum2 + avector(i)*term2
           sum3 = sum3 + avector(i)*term3
           sum4 = sum4 + avector(i)*term4
           sum5 = sum5 + avector(i)*term5
        enddo

        select case(imethod)
        case(0)  ! find local minimum of psi
           pt  = sum
           pt1 = sum1
           pt2 = sum2
           p11 = sum3
           p12 = sum4
           p22 = sum5

           denom = p22*p11 - p12**2
           if(denom.ne.0.) then
              sinew = si -  ( p22*pt1 - p12*pt2)/denom
              etanew= eta - (-p12*pt1 + p11*pt2)/denom
           else
              sinew = si
              etanew= eta
           endif

        case(1)  ! find local zero of <psi,psi>
           pt = sum1**2 + sum2**2
           pt1 = 2.*(sum1*sum3 + sum2*sum4)
           pt2 = 2.*(sum1*sum4 + sum2*sum5)

           denom = pt1**2 + pt2**2
           if(denom.ne.0.) then
              sinew = si - pt*pt1/denom
              etanew = eta - pt*pt2/denom
           else
              sinew = si
              etanew = eta
           end if
        case default
           print *, 'Error: unknown null-finding method: ', imethod
           sinew = si
           etanew = eta
        end select

        xtry = x1 + d%co*(d%b+sinew) - d%sn*etanew
        ztry = z1 + d%sn*(d%b+sinew) + d%co*etanew

!....limit movement to bfac times zone spacing per iteration
        rdiff = sqrt((x-xtry)**2 + (z-ztry)**2)
        if(rdiff .lt. bfac*h) then
          xnew = xtry
          znew = ztry
        else
          xnew = x + bfac*h*(xtry-x)/rdiff
          znew = z + bfac*h*(ztry-z)/rdiff
        endif

        in_domain = 1
        if(rdiff/h .lt. tol) converged = 1
     else
        xnew = 0.
        znew = 0.
        sum   = 0.
        rdiff = 0.
        in_domain = 0
        converged = 0
        izone = 0
     endif  ! on itri.gt.0
     
     ! communicate new minimum to all processors
     if(maxrank.gt.1) then
        temp1(1) = xnew
        temp1(2) = znew
        temp1(3) = sum
        temp1(4) = in_domain
        temp1(5) = converged
        temp1(6) = izone
        call mpi_allreduce(temp1, temp2, 6, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
        xnew  = temp2(1)
        znew  = temp2(2)
        sum   = temp2(3)
        in_domain = temp2(4)
        converged = temp2(5)
        izone = temp2(6)

        if(in_domain .gt. 1) then
           if(myrank.eq.0 .and. iprint.ge.1) &
                print *, "In multiple domains.", in_domain, inews
           if(myrank.le.2 .and. iprint.ge.2) &
                print *,"myrank,x,z", myrank,xnew,znew

           xnew = xnew / in_domain
           znew = znew / in_domain
           sum = sum / in_domain
        end if
     endif
     ! check to see whether the new minimum is outside the simulation domain
     if(in_domain.eq.0) then
        ! if not within the domain, safestop.
        if(myrank.eq.0 .and. iprint.ge.1)   &
             write(*,'(A,2E12.4)') '  magaxis: guess outside domain ', &
             xnew, znew
        ier = 1
        return
     else
        x = xnew
        z = znew
     endif

     if(converged.ge.1) exit newton
  end do newton

  xguess = x
  zguess = z
  psim = sum
  ier = 0
  if(izone.ne.1 .and. iwall_is_limiter.eq.1) ier = 2

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,I12,2E12.4)') '  magaxis: iterations, x, z = ', inews, x, z

  itri_magaxis = itri

end subroutine magaxis

! te_max
! ~~~~~~~
! locates the extreemum in te and the value of te there
! imethod = 0 finds the local maximum
! imethod = 1 finds the local zero of <te,te>
!=====================================================
subroutine te_max(xguess,zguess,te,tem,imethod,ier)
  use basic
  use mesh_mod
  use m3dc1_nint
  use field

  implicit none

  include 'mpif.h'

  real, intent(inout) :: xguess, zguess
  type(field_type), intent(in) :: te
  integer, intent(in) :: imethod
  real, intent(out) :: tem

  integer, parameter :: iterations = 20  !  max number of Newton iterations
  real, parameter :: bfac = 0.1  !max zone fraction for movement each iteration
  real, parameter :: tol = 1e-3   ! convergence tolorance (fraction of h)

  type(element_data) :: d
  integer :: inews
  integer :: i, ier, in_domain, converged
  real :: x1, z1, x, z, si, zi, eta, h
  real :: sum, sum1, sum2, sum3, sum4, sum5
  real :: term1, term2, term3, term4, term5
  real :: pt, pt1, pt2, p11, p22, p12
  real :: xnew, znew, denom, sinew, etanew
  real :: xtry, ztry, rdiff
  vectype, dimension(coeffs_per_element) :: avector
  real, dimension(5) :: temp1, temp2
  integer :: itri

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,2E12.4)') '  te_max: guess = ', xguess, zguess

  converged = 0

  x = xguess
  z = zguess
  itri = itri_te_max

  newton :  do inews=1, iterations

     call whattri(x,0.,z,itri,x1,z1)

     ! calculate position of maximum
     if(itri.gt.0) then
        call calcavector(itri, te, avector)
        call get_element_data(itri, d)

        ! calculate local coordinates
        call global_to_local(d, x, 0., z, si, zi, eta)

        ! calculate mesh size
        h = sqrt((d%a+d%b)*d%c)

        ! evaluate the polynomial and second derivative
        sum = 0.
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        sum5 = 0.
        do i=1, coeffs_per_tri
           sum = sum + avector(i)*si**mi(i)*eta**ni(i)
           term1 = 0.
           if(mi(i).ge.1) term1 = mi(i)*si**(mi(i)-1)*eta**ni(i)
           term2 = 0.
           if(ni(i).ge.1) term2 = ni(i)*si**mi(i)*eta**(ni(i)-1)
           term3 = 0.
           if(mi(i).ge.2) term3 = mi(i)*(mi(i)-1)*si**(mi(i)-2)*eta**ni(i)
           term4 = 0.
           if(ni(i)*mi(i) .ge. 1)                                          &
                term4 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
           term5 = 0.
           if(ni(i).ge.2) term5 = ni(i)*(ni(i)-1)*si**mi(i)*eta**(ni(i)-2)
           
           sum1 = sum1 + avector(i)*term1
           sum2 = sum2 + avector(i)*term2
           sum3 = sum3 + avector(i)*term3
           sum4 = sum4 + avector(i)*term4
           sum5 = sum5 + avector(i)*term5
        enddo

        select case(imethod)
        case(0)  ! find local maximum of te
           pt  = sum
           pt1 = sum1
           pt2 = sum2
           p11 = sum3
           p12 = sum4
           p22 = sum5

           denom = p22*p11 - p12**2
           if(denom.ne.0.) then
              sinew = si -  ( p22*pt1 - p12*pt2)/denom
              etanew= eta - (-p12*pt1 + p11*pt2)/denom
           else
              sinew = si
              etanew= eta
           endif

        case(1)  ! find local zero of <te,te>
           pt = sum1**2 + sum2**2
           pt1 = 2.*(sum1*sum3 + sum2*sum4)
           pt2 = 2.*(sum1*sum4 + sum2*sum5)

           denom = pt1**2 + pt2**2
           if(denom.ne.0.) then
              sinew = si - pt*pt1/denom
              etanew = eta - pt*pt2/denom
           else
              sinew = si
              etanew = eta
           end if
        case default
           print *, 'Error: unknown null-finding method: ', imethod
           sinew = si
           etanew = eta
        end select

        xtry = x1 + d%co*(d%b+sinew) - d%sn*etanew
        ztry = z1 + d%sn*(d%b+sinew) + d%co*etanew

!....limit movement to bfac times zone spacing per iteration
        rdiff = sqrt((x-xtry)**2 + (z-ztry)**2)
        if(rdiff .lt. bfac*h) then
          xnew = xtry
          znew = ztry
        else
          xnew = x + bfac*h*(xtry-x)/rdiff
          znew = z + bfac*h*(ztry-z)/rdiff
        endif

        in_domain = 1
        if(iprint.ge.2) then
           write(*,'(A,4E12.4)') &
                '  te_max: rdiff/h, tol, xnew,znew', rdiff/h, tol, xnew, znew
        end if
        if(rdiff/h .lt. tol) converged = 1
     else
        xnew = 0.
        znew = 0.
        sum   = 0.
        rdiff = 0.
        in_domain = 0
        converged = 0
     endif  ! on itri.gt.0
     
     ! communicate new maximum to all processors
     if(maxrank.gt.1) then
        temp1(1) = xnew
        temp1(2) = znew
        temp1(3) = sum
        temp1(4) = in_domain
        temp1(5) = converged
        call mpi_allreduce(temp1, temp2, 5, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
        xnew  = temp2(1)
        znew  = temp2(2)
        sum   = temp2(3)
        in_domain = temp2(4)
        converged = temp2(5)

        if(in_domain .gt. 1) then
           if(myrank.eq.0 .and. iprint.ge.1) &
                print *, "In multiple domains.", in_domain

           xnew = xnew / in_domain
           znew = znew / in_domain
           sum = sum / in_domain
        end if
     endif
     ! check to see whether the new maximum is outside the simulation domain
     if(in_domain.eq.0) then
        ! if not within the domain, safestop.
        if(myrank.eq.0 .and. iprint.ge.1)   &
             write(*,'(A,2E12.4)') '  te_max: guess outside domain ', &
             xnew, znew
        ier = 1
        return
     else
        x = xnew
        z = znew
     endif

     if(converged.ge.1) exit newton
  end do newton

 ! xguess = x
 ! zguess = z
  tem = sum
  ier = 0

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,I12,2E12.4)') '  te_max: iterations, x, z = ', inews, x, z

  itri_te_max = itri

end subroutine te_max

subroutine te_max2(xguess,zguess,te,tem,imethod,ier)
  use basic
  use mesh_mod
  use m3dc1_nint
  use field

  implicit none

  include 'mpif.h'

  real, intent(inout) :: xguess, zguess
  type(field_type), intent(in) :: te
  integer, intent(in) :: imethod
  real, intent(out) :: tem
  type(element_data) :: d
  integer :: i, ier
  real :: x1, z1, x, z, si, zi, eta
  real :: sum
  vectype, dimension(coeffs_per_element) :: avector
  real, dimension(5) :: temp1, temp2
  integer :: itri

  x = xguess
  z = zguess
  itri = itri_te_max2
  sum = 0.

  call whattri(x,0.,z,itri,x1,z1)

  ! calculate te at x,0,z
  if(itri.gt.0) then
     call calcavector(itri, te, avector)
     call get_element_data(itri, d)
     ! calculate local coordinates
     call global_to_local(d, x, 0., z, si, zi, eta)
     ! evaluate the polynomial
     sum = 0.
     do i=1, coeffs_per_tri
        sum = sum + avector(i)*si**mi(i)*eta**ni(i)
     enddo
  endif  ! on itri.gt.0
  ! communicate new maximum to all processors
  if(maxrank.gt.1) then
     temp1(1) = sum
     call mpi_allreduce(temp1, temp2, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
     sum   = temp2(1)
  endif
  tem = sum
  ier = 0
  itri_te_max2 = itri
end subroutine te_max2



!=====================================================
! lcfs
!
! locates the magnetic axis and the value of psi there
!=====================================================
subroutine lcfs(psi, test_wall, findx)
  use arrays
  use basic
  use mesh_mod
  use m3dc1_nint
  use field
  use boundary_conditions

  implicit none

  include 'mpif.h'

  type(field_type), intent(in) :: psi
  logical, intent(in), optional :: test_wall
  logical, intent(in), optional :: findx

  type(field_type) :: temp_field
  real :: psib, psim
  real :: x, z, temp1, temp2, temp_min, temp_max, phi
  integer :: ier, ier2, numnodes, inode, izone, izonedim, itri, icounter_t
  logical :: is_boundary, first_point
  real, dimension(2) :: normal
  real, dimension(OP_NUM) :: dum1
  real :: curv(3)
  logical :: tw, fx
  vectype, dimension(dofs_per_node) :: data

  if(present(test_wall)) then 
     tw = test_wall
  else
     tw = .true.
  endif
  if(present(findx)) then 
     fx = findx
  else
     fx = .true.
  end if

  call create_field(temp_field)
  temp_field = psi
  if(icsubtract.eq.1) call add(temp_field, psi_coil_field)

  ! Find magnetic axis
  ! ~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Finding magnetic axis.  FX = ', fx
  call magaxis(xmag,zmag,temp_field,psim,0,ier)
  if(ier.eq.0) then
     psimin = psim
     
     if(myrank.eq.0 .and. iprint.ge.1) then 
        write(*,'(A,2E12.4)') '  magnetic axis: ', xmag, zmag
        write(*,'(A, E12.4)') '  psi at magnetic axis: ', psimin
     end if
  else
     if(myrank.eq.0 .and. iprint.ge.1) then 
        write(*,'(A,2E12.4)') '  no magnetic axis found near ', xmag, zmag
     end if
  endif

  if(ifixedb.eq.1 .and. ntime.le.0) then
     psilim = 0.
     psilim2 = 0.
     psibound = 0.
     call destroy_field(temp_field)     
     return
  endif

  ! Find the maximum value of psi at the boundary 
  ! such that psi is not increasing outward 
  ! (as in a private flux region)
  first_point = .true.
  numnodes = owned_nodes()
  if(tw .and. (iwall_is_limiter.eq.1)) then
     do icounter_t=1,numnodes
        inode = nodes_owned(icounter_t) 
        call boundary_node(inode,is_boundary,izone,izonedim,normal,curv, &
             x,phi,z, BOUND_FIRSTWALL)
        if(.not.is_boundary) cycle
        
        call get_node_data(temp_field,inode,data)
        
        if(((x-xmag)*real(data(2)) + (z-zmag)*real(data(3))) &
             *(real(data(1))-psimin).lt.0.) cycle
        if(z*zmag .lt. 0) cycle    !  added 11/21/2018 scj        
        if(first_point) then
           psib = real(data(1))
           first_point =.false.
        else
           if(abs(real(data(1)) - psimin).lt.abs(psib - psimin)) &
                psib = real(data(1))
        endif
     end do
  end if

  if(first_point) then
     temp1 = -1e10
     temp2 =  1e10
  else 
     temp1 = psib
     temp2 = psib
  end if
  call mpi_allreduce(temp1, temp_max, 1, &
       MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call mpi_allreduce(temp2, temp_min, 1, &
       MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  if(abs(temp_max - psimin).lt.abs(temp_min - psimin)) then
     psib = temp_max
  else
     psib = temp_min
  endif

  psibound = psib
  is_diverted=.false.


  ! Calculate psi at the x-point(s)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  psinull = psib
  psinull2 = psib
  if(fx) then 
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Finding X-point'
     if(xnull .gt. 0.) call magaxis(xnull,znull,temp_field,psinull,1,ier)
     if(xnull2 .gt. 0.) call magaxis(xnull2,znull2,temp_field,psinull2,1,ier2)
  else
     itri = 0.
     if(xnull.gt.0.) then
        call evaluate(xnull,0.,znull,dum1,temp_field,itri,ier)
        psinull = dum1(OP_1)
     end if
     if(xnull2.gt.0.) then
        call evaluate(xnull2,0.,znull2,dum1,temp_field,itri,ier2)
        psinull2 = dum1(OP_1)
     end if
  end if

  if(xnull.gt.0) then
     if(ier.eq.0) then
        if(myrank.eq.0 .and. iprint.ge.1) then
           write(*,'(A,2E12.4)') '  X-point 1 found at ', xnull, znull
           write(*,'(A,2E12.4)') '   Psi_x = ', psinull
        end if
        if(abs(psinull - psimin).lt.abs(psibound - psimin)) then
           is_diverted = .true.
           psibound = psinull
        end if
     else 
        if(myrank.eq.0 .and. iprint.ge.1) then 
           write(*,'(A,2E12.4)') '  X-point 1 NOT found near ', xnull, znull
        end if
     endif
  end if
  if(xnull2 .gt. 0) then
     if(ier2.eq.0) then
        if(myrank.eq.0 .and. iprint.ge.1) then
           write(*,'(A,2E12.4)') '  X-point 2 found at ', xnull2, znull2
           write(*,'(A,2E12.4)') '   Psi_x = ', psinull2
        end if
        if(abs(psinull2 - psimin).lt.abs(psibound - psimin)) then
           is_diverted = .true.
           psibound = psinull2
        end if
     else 
        if(myrank.eq.0 .and. iprint.ge.1) then 
           write(*,'(A,2E12.4)') '  X-point 2 NOT found near ', xnull2, znull2
        end if
     endif
  end if


  ! Calculate psi at the limiter
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Finding LCFS'
  if(xlim.eq.0.) then
     ! when xlim = 0, use the lcfs as the limiting flux
     psilim = psibound
     psilim2 = psilim
  else
     itri = 0
     call evaluate(xlim,0.,zlim,dum1,temp_field,itri,ier)
     if(ier.eq.0) then
        psilim = dum1(OP_1)
     else
        psilim = psibound
        if(myrank.eq.0 .and. iprint.ge.1) print *, 'Limiter #1 not found.'
     end if
     
     ! calculate psi at a second limiter point as a diagnostic
     if(xlim2.gt.0) then
        itri = 0
        call evaluate(xlim2,0.,zlim2,dum1,temp_field,itri,ier)
        if(ier.eq.0) then
           psilim2 = dum1(OP_1)
        else
           psilim2 = psibound
           if(myrank.eq.0 .and. iprint.ge.1) print *, 'Limiter #2 not found.'
        end if
     else
        psilim2 = psilim
     endif

     if(abs(psilim - psimin) .lt. abs(psibound - psimin)) then
        is_diverted = .false.
        psibound = psilim
     endif
     if(abs(psilim2 - psimin) .lt. abs(psibound - psimin)) then
        is_diverted = .false.
        psibound = psilim2
     endif
  endif

  if(myrank.eq.0 .and. iprint.ge.1) then
     write(*,'(1A10,6A11)') 'psi at:', &
          'axis', 'wall', 'divertor', 'lim1', 'lim2', 'lcfs'
     write(*,'(1I10,1p6e11.3)') myrank,  &
          psimin, psib, psinull, psilim, psilim2, psibound
  endif

  ! daignostic output
  if(myrank.eq.0 .and. iprint.ge.1) then
     if(psibound.eq.psib) then
        print *, ' Plasma is limited by wall.'
     else if(psibound.eq.psinull) then
        print *, ' Plasma is diverted by x-point #1.'
     else if(psibound.eq.psinull2) then
        print *, ' Plasma is diverted by x-point #2.'
     else if(psibound.eq.psilim) then
        print *, ' Plasma is limited by internal limiter #1.'
     else if(psibound.eq.psilim2) then
        print *, ' Plasma is limited by internal limiter #2.'
     else 
        print *, ' Plasma limiter is unknown!'
     end if
  end if

  call destroy_field(temp_field)
     
end subroutine lcfs

!======================================================================
! get_chord_mask
! ~~~~~~~~~~~~~~
! calculates the transfer function
! exp(-t^2/(2*sigma)^2) / |r - r0|^2
! where cos(t) = d.(r - r0) and 
! d is the unit vector in the director of the chord
!======================================================================
subroutine get_chord_mask(r0, phi0, z0, r, phi, z, npts, theta, sigma, mask)
  implicit none

  integer, intent(in) :: npts      ! number of source points
  real, intent(in) :: r0, phi0, z0 ! location of detector
  real, intent(in), dimension(npts) :: r, phi, z    ! location of source
  real, intent(in) :: theta        ! angle of chord w.r.t. horizontal (radians)
  real, intent(in) :: sigma        ! variance of chord (radians)
  vectype, intent(out), dimension(npts) :: mask

  real, dimension(npts) :: l, t

  ! distance of source to detector
  l = sqrt(r**2 + r0**2 - 2.*r*r0*cos(phi-phi0) + (z-z0)**2)

  ! angle of source to detector relative to angle of chord
  t = acos(((r*cos(phi-phi0) - r0)*cos(theta) + (z-z0)*sin(theta))/l)

  ! assume that signal falls off as l**2
  ! and shape function is guassian in angle from chord
  mask = exp(-t**2/(2.*sigma**2))/l**2
end subroutine get_chord_mask


subroutine calculate_rho(itri)
  use basic
  use kprad
  use kprad_m3dc1
  use m3dc1_nint

  implicit none

  integer, intent(in) :: itri
  integer :: i

  rho79 = nt79

  if(ikprad.ne.0) then 
     do i=1, kprad_z
        call eval_ops(itri, kprad_n(i), tm79, rfac)
        rho79 = rho79 + tm79*kprad_mz/ion_mass
     end do
  end if
  
end subroutine calculate_rho


!======================================================================
! bremsstrahlung
! ~~~~~~~~~~~~~~
! calculates the power per volume of bremsstrahlung given a local
! density n and electron pressure p
!======================================================================
elemental vectype function bremsstrahlung(n, p)
  implicit none
  vectype, intent(in) :: n, p

#ifdef USECOMPLEX
  bremsstrahlung = sqrt(p/n)
#else
  if(p.le.0. .or. n.le.0.) then
     bremsstrahlung = 0.
  else
     bremsstrahlung = sqrt(n*p)
  end if
#endif
end function bremsstrahlung

!======================================================================
! calculate_ke
! ~~~~~~~~~~~~~~
! calculates each Fourer harmonics for kinetic energy
!======================================================================
subroutine calculate_ke()
#ifdef USE3D
  use basic
  use mesh_mod
  use arrays
  use m3dc1_nint
  use newvar_mod
  use sparse
  use metricterms_new
  use boundary_conditions
  use math

  implicit none
  include 'mpif.h'
  integer :: itri, numelms, def_fields
  real :: ke_N, ketotal, fac
  integer :: ier, k, l, numnodes, N, icounter_t
  integer :: izone, izonedim, izone_ind
  vectype, dimension(dofs_per_node) :: vec_l

  real, allocatable :: i1ck(:,:), i1sk(:,:)
  real, allocatable :: i2ck(:,:), i2sk(:,:)

!  type(vector_type) :: transform_field
  type(field_type) :: u_transformc, vz_transformc, chi_transformc
  type(field_type) :: u_transforms, vz_transforms, chi_transforms
#ifdef USEST
  type(field_type) :: vx_field, vy_field, vp_field
  vectype, dimension(dofs_per_element) :: dofs
#endif

  NMAX = ike_harmonics
  numnodes = owned_nodes()
  if(.not.allocated(keharmonic)) allocate(keharmonic(0:NMAX))
  allocate(i1ck(0:nplanes-1,0:NMAX))
  allocate(i1sk(0:nplanes-1,0:NMAX))
  allocate(i2ck(0:nplanes-1,0:NMAX))
  allocate(i2sk(0:nplanes-1,0:NMAX))

  ! create the sin and cos arrays
  do N = 0, NMAX
     do k = 0, nplanes-1
        call ke_I1(NMAX, k, N, i1ck(k,N), i1sk(k,N))
        call ke_I2(NMAX, k, N, i2ck(k,N), i2sk(k,N))
     enddo
  enddo

  if(myrank.eq.0 .and. iprint.eq.1) then
     write(*,900) ntime, numnodes, NMAX
900  format("calculate_ke called-1,   ntime  numnodes  NMAX=",3i6)
  endif

  
  call create_field(  u_transformc)
  call create_field(  u_transforms)
  call create_field( vz_transformc)
  call create_field( vz_transforms)
  call create_field(chi_transformc)
  call create_field(chi_transforms)

#ifdef USEST
  call create_field(vx_field)
  call create_field(vy_field)
  call create_field(vp_field)
  def_fields = FIELD_PHI 
  if(numvar.ge.2) def_fields = def_fields + FIELD_V
  if(numvar.ge.3) def_fields = def_fields + FIELD_CHI
  ! calculate R, phi, Z components of velocity
  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_diag, int_pts_tor)
     call define_fields(itri,def_fields,1,0)

     ! R component: -phi_Z*R + chi_R/R^2
     dofs = - intx3(mu79(:,:,OP_1),pht79(:,OP_DZ),r_79) &
            + intx3(mu79(:,:,OP_1),cht79(:,OP_DR),ri2_79)
     call vector_insert_block(vx_field%vec, itri, 1, dofs, VEC_ADD)

     ! Z component: phi_R*R + chi_Z/R^2
     dofs =   intx3(mu79(:,:,OP_1),pht79(:,OP_DR),r_79) &
            + intx3(mu79(:,:,OP_1),cht79(:,OP_DZ),ri2_79)
     call vector_insert_block(vy_field%vec, itri, 1, dofs, VEC_ADD)

     ! Phi component: omega*R 
     dofs =   intx3(mu79(:,:,OP_1),vzt79(:,OP_1),r_79) 
     call vector_insert_block(vp_field%vec, itri, 1, dofs, VEC_ADD)
  end do

  call newvar_solve(vx_field%vec,mass_mat_lhs)
  call newvar_solve(vy_field%vec,mass_mat_lhs)
  call newvar_solve(vp_field%vec,mass_mat_lhs)
#endif


  if(myrank.eq.0 .and. iprint.eq.1) then
     write(*,901) ntime
901  format("calculate_ke called-2,   ntime=",i6)
  endif

  ! for each Fourier mode
  do N=0,NMAX

     k = local_plane()

     !eq 12: U cos
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(vx_field, l , u1_l )! u1_l is R component of velocity 
#else
        call get_node_data(u_field(1), l , u1_l )! u1_l is “U” (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(u_field(0), l , u0_l )
           u1_l = u1_l + u0_l
        end if
#endif

        vec_l(1)= u1_l(1) * i1ck(k,N) + u1_l( 7)*i2ck(k,N)
        vec_l(2)= u1_l(2) * i1ck(k,N) + u1_l( 8)*i2ck(k,N)
        vec_l(3)= u1_l(3) * i1ck(k,N) + u1_l( 9)*i2ck(k,N)
        vec_l(4)= u1_l(4) * i1ck(k,N) + u1_l(10)*i2ck(k,N)
        vec_l(5)= u1_l(5) * i1ck(k,N) + u1_l(11)*i2ck(k,N)
        vec_l(6)= u1_l(6) * i1ck(k,N) + u1_l(12)*i2ck(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        
        call set_node_data(u_transformc,l,vec_l)
     enddo
     call finalize(u_transformc%vec)
     call m3dc1_field_sum_plane(u_transformc%vec%id) ! sum vec%datator at each (R,Z) node over k

     
     !eq 12: U sin
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(vx_field, l , u1_l )! u1_l is R component of velocity 
#else
        call get_node_data(u_field(1), l , u1_l )! u1_l is “U” (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(u_field(0), l , u0_l )
           u1_l = u1_l + u0_l
        end if
#endif

        vec_l(1)= u1_l(1) * i1sk(k,N) + u1_l( 7)*i2sk(k,N)
        vec_l(2)= u1_l(2) * i1sk(k,N) + u1_l( 8)*i2sk(k,N)
        vec_l(3)= u1_l(3) * i1sk(k,N) + u1_l( 9)*i2sk(k,N)
        vec_l(4)= u1_l(4) * i1sk(k,N) + u1_l(10)*i2sk(k,N)
        vec_l(5)= u1_l(5) * i1sk(k,N) + u1_l(11)*i2sk(k,N)
        vec_l(6)= u1_l(6) * i1sk(k,N) + u1_l(12)*i2sk(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(u_transforms,l,vec_l)
     enddo
     call finalize(u_transforms%vec)
     call m3dc1_field_sum_plane(u_transforms%vec%id) ! sum vec%datator at each (R,Z) node over k


     !eq 12: omega cos
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(vp_field, l , vz1_l )! vz1_l is phi component of velocity 
#else
        call get_node_data(vz_field(1), l , vz1_l) ! vz1_l is “ω” ( dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(vz_field(0), l , vz0_l)
           vz1_l = vz1_l + vz0_l
        end if
#endif

        vec_l(1)= vz1_l(1) * i1ck(k,N) + vz1_l( 7)*i2ck(k,N)
        vec_l(2)= vz1_l(2) * i1ck(k,N) + vz1_l( 8)*i2ck(k,N)
        vec_l(3)= vz1_l(3) * i1ck(k,N) + vz1_l( 9)*i2ck(k,N)
        vec_l(4)= vz1_l(4) * i1ck(k,N) + vz1_l(10)*i2ck(k,N)
        vec_l(5)= vz1_l(5) * i1ck(k,N) + vz1_l(11)*i2ck(k,N)
        vec_l(6)= vz1_l(6) * i1ck(k,N) + vz1_l(12)*i2ck(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(vz_transformc,l,vec_l)
     enddo
     call finalize(vz_transformc%vec)
     call m3dc1_field_sum_plane(vz_transformc%vec%id) ! sum vec%datator at each (R,Z) node over k


     !eq 12: omega sin
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(vp_field, l , vz1_l )! vz1_l is phi component of velocity 
#else
        call get_node_data(vz_field(1), l , vz1_l) ! vz1_l is “ω” ( dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(vz_field(0), l , vz0_l)
           vz1_l = vz1_l + vz0_l
        end if
#endif

        vec_l(1)= vz1_l(1) * i1sk(k,N) + vz1_l( 7)*i2sk(k,N)
        vec_l(2)= vz1_l(2) * i1sk(k,N) + vz1_l( 8)*i2sk(k,N)
        vec_l(3)= vz1_l(3) * i1sk(k,N) + vz1_l( 9)*i2sk(k,N)
        vec_l(4)= vz1_l(4) * i1sk(k,N) + vz1_l(10)*i2sk(k,N)
        vec_l(5)= vz1_l(5) * i1sk(k,N) + vz1_l(11)*i2sk(k,N)
        vec_l(6)= vz1_l(6) * i1sk(k,N) + vz1_l(12)*i2sk(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(vz_transforms,l,vec_l)
     enddo
     call finalize(vz_transforms%vec)
     call m3dc1_field_sum_plane(vz_transforms%vec%id) ! sum vec%datator at each (R,Z) node over k


     !eq 12: chi cos
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(vy_field, l , chi1_l )! chi1_l is Z component of velocity 
#else
        call get_node_data(chi_field(1), l , chi1_l ) ! chi1_l is “χ” (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(chi_field(0), l , chi0_l )
           chi1_l = chi1_l + chi0_l
        end if
#endif

        vec_l(1)= chi1_l(1) * i1ck(k,N) + chi1_l( 7)*i2ck(k,N)
        vec_l(2)= chi1_l(2) * i1ck(k,N) + chi1_l( 8)*i2ck(k,N)
        vec_l(3)= chi1_l(3) * i1ck(k,N) + chi1_l( 9)*i2ck(k,N)
        vec_l(4)= chi1_l(4) * i1ck(k,N) + chi1_l(10)*i2ck(k,N)
        vec_l(5)= chi1_l(5) * i1ck(k,N) + chi1_l(11)*i2ck(k,N)
        vec_l(6)= chi1_l(6) * i1ck(k,N) + chi1_l(12)*i2ck(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(chi_transformc,l,vec_l)
     enddo
     call finalize(chi_transformc%vec)
     call m3dc1_field_sum_plane(chi_transformc%vec%id) ! sum vec%datator of size 6 at each (R,Z) node over k
     

     ! eq 12: chi sin
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(vy_field, l , chi1_l )! chi1_l is Z component of velocity 
#else
        call get_node_data(chi_field(1), l , chi1_l ) ! chi1_l is “χ” (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(chi_field(0), l , chi0_l )
           chi1_l = chi1_l + chi0_l
        end if
#endif

        vec_l(1)= chi1_l(1) * i1sk(k,N) + chi1_l( 7)*i2sk(k,N)
        vec_l(2)= chi1_l(2) * i1sk(k,N) + chi1_l( 8)*i2sk(k,N)
        vec_l(3)= chi1_l(3) * i1sk(k,N) + chi1_l( 9)*i2sk(k,N)
        vec_l(4)= chi1_l(4) * i1sk(k,N) + chi1_l(10)*i2sk(k,N)
        vec_l(5)= chi1_l(5) * i1sk(k,N) + chi1_l(11)*i2sk(k,N)
        vec_l(6)= chi1_l(6) * i1sk(k,N) + chi1_l(12)*i2sk(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(chi_transforms,l,vec_l)
     enddo
     call finalize(chi_transforms%vec)
     call m3dc1_field_sum_plane(chi_transforms%vec%id) ! sum vec%datator at each (R,Z) node over k
     

     !eq 4b: Calculate energy for each Fourier Harminics N
   
     ke_N = 0.
     def_fields = FIELD_N
     numelms = local_elements()
     
!!$OMP PARALLEL DO REDUCTION(+:ke_N)
     do itri=1,numelms
        call m3dc1_ent_getgeomclass(2, itri-1,izonedim,izone_ind)
        izone = zone_type(izone_ind)

        if(izone.ne.1) cycle

        call define_element_quadrature(itri, int_pts_diag, int_pts_tor)
        call define_fields(itri, def_fields, 1, 0)
        call calculate_rho(itri)
!
!       cosine harmonics
        call eval_ops(itri,  u_transformc,pht79)
        call eval_ops(itri, vz_transformc,vzt79)
        call eval_ops(itri,chi_transformc,cht79)

#ifdef USEST
        ke_N = ke_N + int3(rho79(:,OP_1), pht79(:,OP_1), pht79(:,OP_1))
        ke_N = ke_N + int3(rho79(:,OP_1), cht79(:,OP_1), cht79(:,OP_1))
        ke_N = ke_N + int3(rho79(:,OP_1), vzt79(:,OP_1), vzt79(:,OP_1))
#else
        ke_N = ke_N + int4(r2_79, rho79(:,OP_1), pht79(:,OP_DR), pht79(:,OP_DR))
        ke_N = ke_N + int4(r2_79, rho79(:,OP_1), pht79(:,OP_DZ), pht79(:,OP_DZ))

        ke_N = ke_N + int4(r2_79, rho79(:,OP_1), vzt79(:,OP_1), vzt79(:,OP_1))

        ke_N = ke_N + int4(ri4_79, rho79(:,OP_1), cht79(:,OP_DR), cht79(:,OP_DR))
        ke_N = ke_N + int4(ri4_79, rho79(:,OP_1), cht79(:,OP_DZ), cht79(:,OP_DZ))

        ke_N = ke_N + 2.*int4(ri_79, rho79(:,OP_1), pht79(:,OP_DR), cht79(:,OP_DZ))
        ke_N = ke_N - 2.*int4(ri_79, rho79(:,OP_1), pht79(:,OP_DZ), cht79(:,OP_DR))
#endif

!
!       sine harmonics
        call eval_ops(itri,  u_transforms,pht79)
        call eval_ops(itri, vz_transforms,vzt79)
        call eval_ops(itri,chi_transforms,cht79)

#ifdef USEST
        ke_N = ke_N + int3(rho79(:,OP_1), pht79(:,OP_1), pht79(:,OP_1))
        ke_N = ke_N + int3(rho79(:,OP_1), cht79(:,OP_1), cht79(:,OP_1))
        ke_N = ke_N + int3(rho79(:,OP_1), vzt79(:,OP_1), vzt79(:,OP_1))
#else
        ke_N = ke_N + int4(r2_79, rho79(:,OP_1), pht79(:,OP_DR), pht79(:,OP_DR))
        ke_N = ke_N + int4(r2_79, rho79(:,OP_1), pht79(:,OP_DZ), pht79(:,OP_DZ))

        ke_N = ke_N + int4(r2_79, rho79(:,OP_1), vzt79(:,OP_1), vzt79(:,OP_1))

        ke_N = ke_N + int4(ri4_79, rho79(:,OP_1), cht79(:,OP_DR), cht79(:,OP_DR))
        ke_N = ke_N + int4(ri4_79, rho79(:,OP_1), cht79(:,OP_DZ), cht79(:,OP_DZ))

        ke_N = ke_N + 2.*int4(ri_79, rho79(:,OP_1), pht79(:,OP_DR), cht79(:,OP_DZ))
        ke_N = ke_N - 2.*int4(ri_79, rho79(:,OP_1), pht79(:,OP_DZ), cht79(:,OP_DR))
#endif

     end do
!!$OMP END PARALLEL DO

     call mpi_allreduce(ke_N, ketotal, 1, MPI_DOUBLE_PRECISION, &
          MPI_SUM, mpi_comm_world, ier)

! BCL 11/6/19: All transform vectors are constant in phi
!              so integral picks up a 2*pi
!              this is correct for n=0, but twice size for n>0
     if(N.gt.0) ketotal = 0.5*ketotal

     keharmonic(N) = 0.5*ketotal ! 0.5 for 1/2 rho v^2

  end do

!!!!!....we need to save keharmonic for output <===
  if(myrank.eq.0 .and. iprint.ge.1) then
     write(*,1001) ntime
     write(*,1002) (keharmonic(N),N=0,NMAX)
1001 format(" keharmonics at cycle",i6)
1002 format(1p5e12.4)
  endif

  deallocate(i1ck, i1sk, i2ck, i2sk)
  call destroy_field(  u_transformc)
  call destroy_field(  u_transforms)
  call destroy_field( vz_transformc)
  call destroy_field( vz_transforms)
  call destroy_field(chi_transformc)
  call destroy_field(chi_transforms)
#ifdef USEST
  call destroy_field(vx_field)
  call destroy_field(vy_field)
  call destroy_field(vp_field)
#endif
#endif
end subroutine calculate_ke



!======================================================================
! calculate_bh
! ~~~~~~~~~~~~~~
! calculates each Fourer harmonics for magnetic energy
!======================================================================
subroutine calculate_bh()
#ifdef USE3D
  use basic
  use mesh_mod
  use arrays
  use m3dc1_nint
  use newvar_mod
  use sparse
  use metricterms_new
  use boundary_conditions
  use math
  implicit none
  include 'mpif.h'
  integer :: itri, numelms, def_fields
  real:: bh_N, bhtotal, fac
  integer :: ier, k, l, numnodes, N, icounter_t
  integer :: izone, izonedim, izone_ind
  vectype, dimension(dofs_per_node) :: vec_l

  real, allocatable :: i1ck(:,:), i1sk(:,:)
  real, allocatable :: i2ck(:,:), i2sk(:,:)

  type(field_type) :: psi_transformc, F_transformc, fp_transformc
  type(field_type) :: psi_transforms, F_transforms, fp_transforms
#ifdef USEST
  type(field_type) :: bx_field, by_field, bp_field
  vectype, dimension(dofs_per_element) :: dofs
#endif

  BNMAX = ibh_harmonics
  numnodes = owned_nodes()
  if(.not.allocated(bharmonic)) allocate(bharmonic(0:BNMAX))
  allocate(i1ck(0:nplanes-1,0:BNMAX))
  allocate(i1sk(0:nplanes-1,0:BNMAX))
  allocate(i2ck(0:nplanes-1,0:BNMAX))
  allocate(i2sk(0:nplanes-1,0:BNMAX))
  
  ! create the sin and cos arrays
  do N = 0, BNMAX
     do k = 0, nplanes-1
        call ke_I1(BNMAX, k, N, i1ck(k,N), i1sk(k,N))
        call ke_I2(BNMAX, k, N, i2ck(k,N), i2sk(k,N))
     enddo
  enddo
  
  if(myrank.eq.0 .and. iprint.eq.1) then
     write(*,900) ntime, numnodes, BNMAX
900  format("calculate_bh called-1,   ntime  numnodes  BNMAX=",3i6)
  endif
  
  call create_field(psi_transformc)
  call create_field(psi_transforms)
  call create_field(F_transformc)
  call create_field(F_transforms)
  call create_field(fp_transformc)
  call create_field(fp_transforms)

#ifdef USEST
  call create_field(bx_field)
  call create_field(by_field)
  call create_field(bp_field)
  def_fields = FIELD_PSI 
  if(numvar.ge.2) def_fields = def_fields + FIELD_I 
  ! calculate R, phi, Z components of magnetic field 
  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_diag,int_pts_tor)
     call define_fields(itri,def_fields,1,0)

     ! R component: -psi_Z/R - f'_R
     dofs = - intx3(mu79(:,:,OP_1),pst79(:,OP_DZ),ri_79) &
            - intx2(mu79(:,:,OP_1),bfpt79(:,OP_DR))
     call vector_insert_block(bx_field%vec, itri, 1, dofs, VEC_ADD)

     ! Z component: psi_R/R - f'_Z
     dofs =   intx3(mu79(:,:,OP_1),pst79(:,OP_DR),ri_79) &
            - intx2(mu79(:,:,OP_1),bfpt79(:,OP_DZ))
     call vector_insert_block(by_field%vec, itri, 1, dofs, VEC_ADD)

     ! Phi component: F/R 
     dofs =   intx3(mu79(:,:,OP_1),bzt79(:,OP_1),ri_79) 
     call vector_insert_block(bp_field%vec, itri, 1, dofs, VEC_ADD)
  end do

  call newvar_solve(bx_field%vec,mass_mat_lhs)
  call newvar_solve(by_field%vec,mass_mat_lhs)
  call newvar_solve(bp_field%vec,mass_mat_lhs)
#endif
 
  if(myrank.eq.0 .and. iprint.eq.1) then
     write(*,901) ntime
901  format("calculate_bh called-2,   ntime=",i6)
  endif
  

  ! for each Fourier mode
  do N=0,BNMAX

     k = local_plane()

     !eq 12: psi cos
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(bx_field, l , psi1_l )! psi1_l is R component of field 
#else
        call get_node_data(psi_field(1), l, psi1_l) ! psi1_1 is ψ (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(psi_field(0), l, psi0_l)
           psi1_l = psi1_l + psi0_l
        end if
        if(icsubtract.eq.1) then
           call get_node_data(psi_coil_field, l, psi0_l)
           psi1_l = psi1_l + psi0_l
        end if
#endif

        vec_l(1)= psi1_l(1) * i1ck(k,N) + psi1_l( 7)*i2ck(k,N)
        vec_l(2)= psi1_l(2) * i1ck(k,N) + psi1_l( 8)*i2ck(k,N)
        vec_l(3)= psi1_l(3) * i1ck(k,N) + psi1_l( 9)*i2ck(k,N)
        vec_l(4)= psi1_l(4) * i1ck(k,N) + psi1_l(10)*i2ck(k,N)
        vec_l(5)= psi1_l(5) * i1ck(k,N) + psi1_l(11)*i2ck(k,N)
        vec_l(6)= psi1_l(6) * i1ck(k,N) + psi1_l(12)*i2ck(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        
        call set_node_data(psi_transformc,l,vec_l)
     enddo
     call finalize(psi_transformc%vec)
     call m3dc1_field_sum_plane(psi_transformc%vec%id) ! sum vec%datator at each (R,Z) node over k
     
     !eq 12: psi sin
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(bx_field, l , psi1_l )! psi1_l is R component of field 
#else
        call get_node_data(psi_field(1), l, psi1_l) ! psi1_1 is ψ (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(psi_field(0), l, psi0_l)
           psi1_l = psi1_l + psi0_l
        end if
        if(icsubtract.eq.1) then
           call get_node_data(psi_coil_field, l, psi0_l)
           psi1_l = psi1_l + psi0_l
        end if
#endif

        vec_l(1)= psi1_l(1) * i1sk(k,N) + psi1_l( 7)*i2sk(k,N)
        vec_l(2)= psi1_l(2) * i1sk(k,N) + psi1_l( 8)*i2sk(k,N)
        vec_l(3)= psi1_l(3) * i1sk(k,N) + psi1_l( 9)*i2sk(k,N)
        vec_l(4)= psi1_l(4) * i1sk(k,N) + psi1_l(10)*i2sk(k,N)
        vec_l(5)= psi1_l(5) * i1sk(k,N) + psi1_l(11)*i2sk(k,N)
        vec_l(6)= psi1_l(6) * i1sk(k,N) + psi1_l(12)*i2sk(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(psi_transforms,l,vec_l)
     enddo
     call finalize(psi_transforms%vec)
     call m3dc1_field_sum_plane(psi_transforms%vec%id) ! sum vec%datator at each (R,Z) node over k

     !eq 12: F cos
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(bp_field, l , bz1_l )! bz1_l is phi component of field 
#else
        call get_node_data(bz_field(1), l, bz1_l) ! bz1_l is F (dimension 12)
        if(eqsubtract.eq.1) then
          call get_node_data(bz_field(0), l, bz0_l)
          bz1_l = bz1_l + bz0_l
       end if
#endif

        vec_l(1)= bz1_l(1) * i1ck(k,N) + bz1_l( 7)*i2ck(k,N)
        vec_l(2)= bz1_l(2) * i1ck(k,N) + bz1_l( 8)*i2ck(k,N)
        vec_l(3)= bz1_l(3) * i1ck(k,N) + bz1_l( 9)*i2ck(k,N)
        vec_l(4)= bz1_l(4) * i1ck(k,N) + bz1_l(10)*i2ck(k,N)
        vec_l(5)= bz1_l(5) * i1ck(k,N) + bz1_l(11)*i2ck(k,N)
        vec_l(6)= bz1_l(6) * i1ck(k,N) + bz1_l(12)*i2ck(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(F_transformc,l,vec_l)
     enddo
     call finalize(F_transformc%vec)
     call m3dc1_field_sum_plane(F_transformc%vec%id) ! sum vec%datator at each (R,Z) node over k


     !eq 12: F sin
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(bp_field, l , bz1_l )! bz1_l is phi component of field 
#else
        call get_node_data(bz_field(1), l, bz1_l) ! bz1_l is F (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(bz_field(0), l, bz0_l)
           bz1_l = bz1_l + bz0_l
        end if
#endif

        vec_l(1)= bz1_l(1) * i1sk(k,N) + bz1_l( 7)*i2sk(k,N)
        vec_l(2)= bz1_l(2) * i1sk(k,N) + bz1_l( 8)*i2sk(k,N)
        vec_l(3)= bz1_l(3) * i1sk(k,N) + bz1_l( 9)*i2sk(k,N)
        vec_l(4)= bz1_l(4) * i1sk(k,N) + bz1_l(10)*i2sk(k,N)
        vec_l(5)= bz1_l(5) * i1sk(k,N) + bz1_l(11)*i2sk(k,N)
        vec_l(6)= bz1_l(6) * i1sk(k,N) + bz1_l(12)*i2sk(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(F_transforms,l,vec_l)
     enddo
     call finalize(F_transforms%vec)
     call m3dc1_field_sum_plane(F_transforms%vec%id) ! sum vec%datator at each (R,Z) node over k
     
     !eq 12: f' cos
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(by_field, l , bfp1_l )! bfp1_l is Z component of field 
#else
        call get_node_data(bfp_field(1), l, bfp1_l) ! bfp1_l is f (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(bfp_field(0), l, bfp0_l)
           bfp1_l = bfp1_l + bfp0_l
        end if
#endif

        vec_l(1)= bfp1_l(1) * i1ck(k,N) + bfp1_l( 7)*i2ck(k,N)
        vec_l(2)= bfp1_l(2) * i1ck(k,N) + bfp1_l( 8)*i2ck(k,N)
        vec_l(3)= bfp1_l(3) * i1ck(k,N) + bfp1_l( 9)*i2ck(k,N)
        vec_l(4)= bfp1_l(4) * i1ck(k,N) + bfp1_l(10)*i2ck(k,N)
        vec_l(5)= bfp1_l(5) * i1ck(k,N) + bfp1_l(11)*i2ck(k,N)
        vec_l(6)= bfp1_l(6) * i1ck(k,N) + bfp1_l(12)*i2ck(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(fp_transformc,l,vec_l)
     enddo
     call finalize(fp_transformc%vec)
     call m3dc1_field_sum_plane(fp_transformc%vec%id) ! sum vec%datator of size 6 at each (R,Z) node over k
     
     ! eq 12: f' sin
     do icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
#ifdef USEST
        call get_node_data(by_field, l , bfp1_l )! bfp1_l is Z component of field 
#else
        call get_node_data(bfp_field(1), l, bfp1_l) ! bfp1_l is f (dimension 12)
        if(eqsubtract.eq.1) then
           call get_node_data(bfp_field(0), l, bfp0_l)
           bfp1_l = bfp1_l + bfp0_l
        end if
#endif

        vec_l(1)= bfp1_l(1) * i1sk(k,N) + bfp1_l( 7)*i2sk(k,N)
        vec_l(2)= bfp1_l(2) * i1sk(k,N) + bfp1_l( 8)*i2sk(k,N)
        vec_l(3)= bfp1_l(3) * i1sk(k,N) + bfp1_l( 9)*i2sk(k,N)
        vec_l(4)= bfp1_l(4) * i1sk(k,N) + bfp1_l(10)*i2sk(k,N)
        vec_l(5)= bfp1_l(5) * i1sk(k,N) + bfp1_l(11)*i2sk(k,N)
        vec_l(6)= bfp1_l(6) * i1sk(k,N) + bfp1_l(12)*i2sk(k,N)
        vec_l(7:12) = 0. ! pad with zeros
        call set_node_data(fp_transforms,l,vec_l)
     enddo
     call finalize(fp_transforms%vec)
     call m3dc1_field_sum_plane(fp_transforms%vec%id) ! sum vec%datator at each (R,Z) node over k
     
     
     !eq 4b: Calculate energy for each Fourier Harminics N
     
     bh_N = 0.
     def_fields = 0
     numelms = local_elements()
     
!!$OMP PARALLEL DO REDUCTION(+:bh_N)
     do itri=1,numelms
        call m3dc1_ent_getgeomclass(2, itri-1,izonedim,izone_ind)
        izone = zone_type(izone_ind)

        if(izone.ne.1) cycle

        call define_element_quadrature(itri, int_pts_diag, int_pts_tor)
        call define_fields(itri, def_fields, 1, 0)
        
!       cosine harmonics
        call eval_ops(itri,psi_transformc,pst79)
        call eval_ops(itri,F_transformc,bzt79)
        call eval_ops(itri,fp_transformc,bfpt79)

#ifdef USEST
        bh_N = bh_N + int2(pst79(:,OP_1), pst79(:,OP_1))   
        bh_N = bh_N + int2(bfpt79(:,OP_1), bfpt79(:,OP_1))
        bh_N = bh_N + int2(bzt79(:,OP_1), bzt79(:,OP_1))
#else
        bh_N = bh_N + int3(ri2_79, pst79(:,OP_DR), pst79(:,OP_DR))   &
                    + int3(ri2_79, pst79(:,OP_DZ), pst79(:,OP_DZ))

        bh_N = bh_N + int3(ri2_79, bzt79(:,OP_1), bzt79(:,OP_1))

        bh_N = bh_N + int2(bfpt79(:,OP_DR), bfpt79(:,OP_DR))   &
                    + int2(bfpt79(:,OP_DZ), bfpt79(:,OP_DZ))
        bh_N = bh_N - 2.*int3(ri_79, pst79(:,OP_DR), bfpt79(:,OP_DZ)) &
                    + 2.*int3(ri_79, pst79(:,OP_DZ), bfpt79(:,OP_DR))
#endif

!       sine harmonics
        call eval_ops(itri,psi_transforms,pst79)
        call eval_ops(itri,F_transforms,bzt79)
        call eval_ops(itri,fp_transforms,bfpt79)

#ifdef USEST
        bh_N = bh_N + int2(pst79(:,OP_1), pst79(:,OP_1))   
        bh_N = bh_N + int2(bfpt79(:,OP_1), bfpt79(:,OP_1))
        bh_N = bh_N + int2(bzt79(:,OP_1), bzt79(:,OP_1))
#else
        bh_N = bh_N + int3(ri2_79,  pst79(:,OP_DR), pst79(:,OP_DR))   &
                    + int3(ri2_79,  pst79(:,OP_DZ), pst79(:,OP_DZ))

        bh_N = bh_N + int3(ri2_79,  bzt79(:,OP_1), bzt79(:,OP_1))

        bh_N = bh_N + int2(bfpt79(:,OP_DR), bfpt79(:,OP_DR))   &
                    + int2(bfpt79(:,OP_DZ), bfpt79(:,OP_DZ))
        bh_N = bh_N - 2.*int3(ri_79, pst79(:,OP_DR), bfpt79(:,OP_DZ)) &
                    + 2.*int3(ri_79, pst79(:,OP_DZ), bfpt79(:,OP_DR))
#endif

     end do
!!$OMP END PARALLEL DO

     call mpi_allreduce(bh_N, bhtotal, 1, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, mpi_comm_world, ier)

! BCL 11/6/19: All transform vectors are constant in phi
!              so integral picks up a 2*pi
!              this is correct for n=0, but twice size for n>0
     if(N.gt.0) bhtotal = 0.5*bhtotal

     ! 0.5 for proper normalization
     bharmonic(N) = 0.5*bhtotal
  end do
!    NOTE:  bharmonic must be divided by (2 pi)**2 mu_0 to get actual SI magnetic energy
!           This is done in the idl routine plot_bhmn.pro
!    BCL 11/7/19: Not sure this ^ is correct

!  save one harmonic to scale hyper for ihypeta .gt. 2
  if(ihypeta.gt.2) bharhypeta = bharmonic(ihypeta)


  if(myrank.eq.0 .and. iprint.ge.1) then
     write(*,1001) ntime
     write(*,1002) (bharmonic(N),N=0,BNMAX)
1001 format(" calculate_bh called-2,   bharmonics at cycle",i6)
1002 format(1p5e12.4)
  endif
  
  deallocate(i1ck, i1sk, i2ck, i2sk)

  call destroy_field(psi_transformc)
  call destroy_field(psi_transforms)
  call destroy_field(F_transformc)
  call destroy_field(F_transforms)
  call destroy_field(fp_transformc)
  call destroy_field(fp_transforms)
#ifdef USEST
  call destroy_field(bx_field)
  call destroy_field(by_field)
  call destroy_field(bp_field)
#endif
#endif
end subroutine calculate_bh


!...........................................................
! input: k, N
! output: i1ck, i1sk
! eq 10
subroutine ke_I1(NMAX, k, N, i1ck, i1sk)
  use basic, ONLY: myrank, itor
  use math
  use mesh_mod
  implicit none
  integer:: k, N, NMAX
  real:: i1ck, i1sk
  
  real:: hm, hp
  real:: xm, x0, xp
  
  if(k==0) then
     call m3dc1_plane_getphi(nplanes-1, xm)
     if(itor.eq.0) xm = xm*2.*pi/toroidal_period
     xm = xm - 2.*pi
  else
     call m3dc1_plane_getphi(k-1, xm)
     if(itor.eq.0) xm = xm*2.*pi/toroidal_period
  endif
  call m3dc1_plane_getphi(k, x0)
  if(itor.eq.0) x0 = x0*2.*pi/toroidal_period
  if(k==nplanes-1) then
     call m3dc1_plane_getphi(0, xp)
     if(itor.eq.0) xp = xp*2.*pi/toroidal_period
     xp = xp + 2.*pi
  else
     call m3dc1_plane_getphi(k+1, xp)
     if(itor.eq.0) xp = xp*2.*pi/toroidal_period
  endif
  hm = x0 - xm
  hp = xp - x0
  
  if(N .eq. 0) then
     i1ck = (hm + hp)/(4.*pi)
     i1sk = 0.
     return
  endif

  i1ck = ( hm*N*((6 + hm**2*N**2)*Sin(N*x0) + 6*Sin(N*(-hm + x0))) & 
          + 12*(Cos(N*x0) - Cos(N*(-hm + x0))))/(hm**3*N**4) &
       + (-hp*N*((6 + hp**2*N**2)*Sin(N*x0) + 6*Sin(N*(hp + x0))) &
          + 12*(Cos(N*x0) - Cos(N*( hp + x0))))/(hp**3*N**4)
  i1ck = i1ck/pi

  i1sk = (-hm*N*(6 + hm**2*N**2)*Cos(N*x0) - 6*hm*N*Cos(N*(-hm + x0)) &
          + 12*(Sin(N*(hm - x0)) + Sin(N*x0)))/(hm**3*N**4) &
       + ( hp*N*(6 + hp**2*N**2)*Cos(N*x0) + 6*hp*N*Cos(N*( hp + x0)) &
          + 12*(Sin(N*x0) - Sin(N*(hp + x0))))/(hp**3*N**4)
  i1sk = i1sk/pi
  
end subroutine ke_I1


!...........................................................
! input: N, k, delta_phi
! output: i2ck, i2sk
! eq 10
subroutine ke_I2(NMAX, k, N, i2ck, i2sk)
  use basic, ONLY: myrank, itor
  use math
  use mesh_mod
  implicit none
  integer:: k, N, NMAX
  real:: i2ck, i2sk
  
  real:: hm, hp
  real:: xm, x0, xp
  
  if(k==0) then
     call m3dc1_plane_getphi(nplanes-1, xm)
     if(itor.eq.0) xm = xm*2.*pi/toroidal_period
     xm = xm - 2.*pi
  else
     call m3dc1_plane_getphi(k-1, xm)
     if(itor.eq.0) xm = xm*2.*pi/toroidal_period
  endif
  call m3dc1_plane_getphi(k, x0)
  if(itor.eq.0) x0 = x0*2.*pi/toroidal_period
  if(k==nplanes-1) then
     call m3dc1_plane_getphi(0, xp)
     if(itor.eq.0) xp = xp*2.*pi/toroidal_period
     xp = xp + 2.*pi
  else
     call m3dc1_plane_getphi(k+1, xp)
     if(itor.eq.0) xp = xp*2.*pi/toroidal_period
  endif
  hm = x0 - xm
  hp = xp - x0
  
  if(N .le. 0) then
     i2ck = (hp**2 - hm**2)/(24.*pi)
     i2sk = 0.
     return
  endif

  i2ck = ((-6 + hm**2*N**2)*Cos(N*x0) + 6*Cos(N*(-hm + x0)) & 
          -2*hm*N*(2*Sin(N*x0) + Sin(N*(-hm + x0))))/(hm**2*N**4) &
       + (( 6 - hp**2*N**2)*Cos(N*x0) - 6*Cos(N*( hp + x0)) &
          -2*hp*N*(2*Sin(N*x0) + Sin(N*( hp + x0))))/(hp**2*N**4)
  i2ck = i2ck/pi

  i2sk =  (hm*N*(4*Cos(N*x0) + 2*Cos(N*(-hm + x0)) + hm*N*Sin(N*x0)) &
           - 6*(Sin(N*(hm - x0)) + Sin(N*x0)))/(hm**2*N**4) &
         +(hp*N*(4*Cos(N*x0) + 2*Cos(N*( hp + x0)) - hp*N*Sin(N*x0)) &
           + 6*(Sin(N*x0) - Sin(N*(hp + x0))))/(hp**2*N**4)
  i2sk = i2sk/pi
  
end subroutine ke_I2


! te_max
! ~~~~~~~
! searches each cell for the extreemum in te and the value of te there
!  finds the local global maximum of te
!=====================================================
subroutine te_max_dev(xguess,zguess,te,tem,imethod,ier)
  use basic
  use mesh_mod
  use m3dc1_nint
  use field

  implicit none

  include 'mpif.h'

  real, intent(inout) :: xguess, zguess
  type(field_type), intent(in) :: te
  real, intent(out) :: tem

  integer, parameter :: iterations = 20  !  max number of Newton iterations
  real, parameter :: bfac = 0.1  !max zone fraction for movement each iteration
  real, parameter :: tol = 1e-3   ! convergence tolorance (fraction of h)

  type(element_data) :: d
  integer :: inews
  integer :: i, ier, in_domain, converged
  real :: x1, z1, x, z, si, zi, eta, h, phi
  real :: sum, sum1, sum2, sum3, sum4, sum5
  real :: term1, term2, term3, term4, term5
  real :: pt, pt1, pt2, p11, p22, p12
  real :: xnew, znew, denom, sinew, etanew
  real :: xtry, ztry, rdiff
  vectype, dimension(coeffs_per_element) :: avector
  real :: temp1, temp2
  real, dimension(nodes_per_element) :: xnode , phinode, znode
  real :: xsum, zsum, psum, summax
  integer :: itri, itri1, numelms, inode, imethod
  integer :: nodeids(nodes_per_element)

! search over all triangles local to this processor
!
  numelms = local_elements()
  summax = 0.
  sum = 0.
   if(myrank.eq.0 .and. iprint.ge.2) &
     write(*,'(A,5i5)') 'parameters',iterations,coeffs_per_element,nodes_per_element, &
                         numelms      !,coeffs_per_tri
  triangles : do itri = 1,numelms

    call get_element_nodes(itri,nodeids)

!   calculate the triangle center
      xsum = 0
      zsum = 0
      psum = 0
      do inode=1,3
         call get_node_pos(nodeids(inode),xnode(inode),phinode(inode),znode(inode))
         xsum = xsum + xnode(inode)
         zsum = zsum + znode(inode)
         psum = psum + phinode(inode)
      enddo
      x = xsum/3.
      z = zsum/3.
      phi = psum/3.

                 if(myrank.eq.0 .and. iprint.ge.2) &
                   write(*,'(A,i5,1p3e12.4)') '  itri ,x,z,phi = ', itri, x, z,phi

      converged = 0
  
    newton :  do inews=1, iterations



       call whattri(x,phi,z,itri1,x1,z1)
       if(myrank.eq.0 .and. iprint.ge.2) &
         write(*,'(A,2i5,1p2e12.4)') 'itri,itri1,x1,z1',itri,itri1,x1,z1

     ! calculate position of maximum
       if(itri1.eq.itri) then
        call calcavector(itri, te, avector)
        call get_element_data(itri, d)

        ! calculate local coordinates
        call global_to_local(d, x, phi , z, si, zi, eta)

        ! calculate mesh size
        h = sqrt((d%a+d%b)*d%c)

       if(myrank.eq.0) write(*,'(A,1p6e12.4)') 'triangle data',d%a,d%b,d%c,si,zi,eta

        ! evaluate the polynomial and second derivative
        sum = 0.
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        sum5 = 0.
    exit newton
        do i=1, coeffs_per_tri
           sum = sum + avector(i)*si**mi(i)*eta**ni(i)
           term1 = 0.
           if(mi(i).ge.1) term1 = mi(i)*si**(mi(i)-1)*eta**ni(i)
           term2 = 0.
           if(ni(i).ge.1) term2 = ni(i)*si**mi(i)*eta**(ni(i)-1)
           term3 = 0.
           if(mi(i).ge.2) term3 = mi(i)*(mi(i)-1)*si**(mi(i)-2)*eta**ni(i)
           term4 = 0.
           if(ni(i)*mi(i) .ge. 1)                                          &
                term4 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
           term5 = 0.
           if(ni(i).ge.2) term5 = ni(i)*(ni(i)-1)*si**mi(i)*eta**(ni(i)-2)
           
           sum1 = sum1 + avector(i)*term1
           sum2 = sum2 + avector(i)*term2
           sum3 = sum3 + avector(i)*term3
           sum4 = sum4 + avector(i)*term4
           sum5 = sum5 + avector(i)*term5
        enddo


          ! find local maximum of te
           pt  = sum
           pt1 = sum1
           pt2 = sum2
           p11 = sum3
           p12 = sum4
           p22 = sum5

           denom = p22*p11 - p12**2
           if(denom.ne.0.) then
              sinew = si -  ( p22*pt1 - p12*pt2)/denom
              etanew= eta - (-p12*pt1 + p11*pt2)/denom
           else
              sinew = si
              etanew= eta
           endif


        xtry = x1 + d%co*(d%b+sinew) - d%sn*etanew
        ztry = z1 + d%sn*(d%b+sinew) + d%co*etanew
        rdiff = sqrt((x-xtry)**2 + (z-ztry)**2)

        if(rdiff/h .lt. tol) converged = 1
        x = xtry
        z = ztry
     else ! on itri.eq.itri1
        xnew = 0.
        znew = 0.
        sum   = 0.
        rdiff = 0.
        in_domain = 0
        converged = 0
        exit newton
     endif  ! on itri.eq.itri1
  

     if(converged.ge.1) exit newton
    end do newton
!
    summax = max(sum,summax)
  end do triangles

  ! select maximum over all processors
  if(maxrank.gt.1) then
     temp1 = summax
     call mpi_allreduce(temp1, temp2, 1, &
          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
     summax   = temp2
  endif

  tem = summax
  ier = 0

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,E12.4)') '  te_max_dev: summax', summax
  
end subroutine te_max_dev

  subroutine evaluate_mag_probes()
    use basic
    use arrays
    use m3dc1_nint
    
    implicit none

    integer :: ierr, i
    real, dimension(OP_NUM) :: val
    real :: r

    do i=1, imag_probes
       mag_probe_val(i) = 0.

       if(itor.eq.1) then
          r = mag_probe_x(i)
       else
          r = 1.
       end if

       ! Read poloidal field
       if(mag_probe_nx(i).ne.0. .or. mag_probe_nz(i).ne.0.) then
          ! psi
          call evaluate(mag_probe_x(i),mag_probe_phi(i),mag_probe_z(i), &
               val,psi_field(1),mag_probe_itri(i),ierr)
          if(ierr.ne.0) then
             if(myrank.eq.0) print *, 'Error evaluating flux loop ', i
             cycle
          end if

          mag_probe_val(i) = mag_probe_val(i) &
               - mag_probe_nx(i)*val(OP_DZ) / r &
               + mag_probe_nz(i)*val(OP_DR) / r

          ! f
#if defined(USE3D) || defined(USECOMPLEX)
          call evaluate(mag_probe_x(i),mag_probe_phi(i),mag_probe_z(i), &
               val,bf_field(1),mag_probe_itri(i),ierr)      
          if(ierr.ne.0) then
             if(myrank.eq.0) print *, 'Error evaluating flux loop ', i
             cycle
          end if

          mag_probe_val(i) = mag_probe_val(i) &
               + mag_probe_nx(i)*val(OP_DRP) &
               + mag_probe_nz(i)*val(OP_DZP)
#endif
       end if
       
       ! Read toroidal field
       if(mag_probe_nphi(i).ne.0.) then
#ifdef USEPARTICLES
          ! phi
          call evaluate(mag_probe_x(i),mag_probe_phi(i),mag_probe_z(i), &
               val,u_field(1),mag_probe_itri(i),ierr)
#else
          ! bz
          call evaluate(mag_probe_x(i),mag_probe_phi(i),mag_probe_z(i), &
               val,bz_field(1),mag_probe_itri(i),ierr)
#endif
          if(ierr.ne.0) then
             if(myrank.eq.0) print *, 'Error evaluating flux loop ', i
             cycle
          end if

#ifdef USEPARTICLES
          mag_probe_val(i) = mag_probe_val(i) + mag_probe_nphi(i)*val(OP_1)
#else
          mag_probe_val(i) = mag_probe_val(i) + mag_probe_nphi(i)*val(OP_1)*r
#endif
       end if
    end do
  end subroutine evaluate_mag_probes


  subroutine evaluate_flux_loops()
    use basic
    use arrays
    
    implicit none

    integer :: ierr, i

    do i=1, iflux_loops
       call phi_int(flux_loop_x(i),flux_loop_z(i),flux_loop_val(i),psi_field(1), &
            flux_loop_itri(i),ierr)
       if(ierr .ne. 0) then
          flux_loop_val(i) = 0
       end if
    end do

  end subroutine evaluate_flux_loops


!============================================================
! phi_int
! ~~~~~~~
! calculates the integral of a field at (R,Z) over phi
!============================================================
subroutine phi_int(x,z,ans,fin,itri,ierr)
  
  use mesh_mod
  use basic
  use m3dc1_nint
  use field

  implicit none

  include 'mpif.h'

  integer, intent(inout) :: itri
  real, intent(in) :: x, z
  type(field_type), intent(in) :: fin
  integer, intent(out) :: ierr ! = 0 on success

  real, intent(out) :: ans

#if defined(USECOMPLEX)
  ans = 0.
  ierr = 0

#elif defined(USE3D)

  type(element_data) :: d
  integer :: nodeids(nodes_per_element), ier, iplane
  real :: phi, x1, phi1, z1, my_ans
  integer :: hasval, tothasval

  iplane = local_plane()
  call m3dc1_plane_getphi(iplane, phi)
  if(myrank.eq.0) print *, 'diagnostics Plane ', iplane, 'at angle ', phi

  if(itri.eq.0) then
     call whattri(x,phi,z,itri,x1,z1)
  else if(itri.gt.0) then
     call get_element_nodes(itri,nodeids)
     call get_node_pos(nodeids(1), x1, phi1, z1)
  endif

  ans = 0.
  my_ans = 0.

  ! if this process contains the point, evaluate the field at that point.
  if(itri.gt.0) then

     call get_element_data(itri, d)

     ! calculate local coordinates
     call global_to_local(d, x, phi, z, xi_79(1), zi_79(1), eta_79(1))

     weight_79 = 1
     call extrude_quadrature(d%d, 1, 5)

     ! calculate the inverse radius
     if(itor.eq.1) then
        ri_79(1) = 1./x
     else
        ri_79(1) = 1.
     endif

     call precalculate_terms(xi_79, zi_79, eta_79, d%co, d%sn, ri_79, npoints)
     call define_basis(itri)

     ! calculate the value of the function
     call eval_ops(itri, fin, tm79, rfac)
     
     ! integrate
     my_ans = int1(tm79(:,OP_1))

     hasval = 1
  else
     hasval = 0
  endif


  ! Distribute the result if necessary
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! NOT YET DONE -- REDUCE VALUE WITHIN PLANE!
  
  if(maxrank.gt.1) then
     ! Determine number of processes whose domains contain this point
     call mpi_allreduce(hasval, tothasval, 1, MPI_INTEGER, MPI_SUM, &
          MPI_COMM_WORLD, ier)
  else
     tothasval = hasval
  end if

  if(tothasval.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          write(*,'(A,3f12.4)') 'Point not found in domain: ', x, phi, z
     ierr = 1
     return
  end if

  if(maxrank.gt.1) then
     ! Sum the integrals from each plane
     call mpi_allreduce(my_ans, ans, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_WORLD, ier)
  endif

  ierr = 0

#else
  
  real, dimension(OP_NUM) :: temp
  call evaluate(x,0.,z,temp,fin,itri,ierr)

  if(ifull_torus.eq.1) then
     ans = temp(OP_1)*toroidal_period
  else
     ans = temp(OP_1)*toroidal_period*nperiods
  endif
#endif

end subroutine phi_int

  
end module diagnostics

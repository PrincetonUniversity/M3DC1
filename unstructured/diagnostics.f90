#ifdef USECOMPLEX
#define CONJUGATE(x) conjg(x)
#else
#define CONJUGATE(x) x
#endif


module diagnostics

  implicit none

  real :: tflux0, totcur0

  ! scalars integrated over entire computational domain
  real :: tflux, area, volume, totcur, totden, tmom, tvor, bwb2

  ! scalars integrated within lcfs
  real :: pflux, parea, pcur, pden, pmom, pvor

  real :: chierror

  ! density diagnostics
  real :: nfluxd, nfluxv, nsource

  ! energy diagnostics
  real :: ekin, emag, ekind, emagd, ekino, emago, ekindo, emagdo,      &
       ekint,emagt,ekintd,emagtd,ekinto,emagto,ekintdo,emagtdo,        &
       ekinp,emagp,ekinpd,emagpd,ekinpo,emagpo,ekinpdo,emagpdo,        &
       ekinph,ekinth,emagph,emagth,ekinpho,ekintho,emagpho,emagtho,    &
       ekin3,ekin3d,ekin3h,emag3,ekin3o,ekin3do,ekin3ho,emag3o,        &
       emag3h,emag3d,emag3ho,emag3do, avep
  real :: efluxp,efluxk,efluxs,efluxt,epotg,etot,ptot,eerr,ptoto

  ! in subroutine calculate_ke()
  real, allocatable:: keharmonic(:)
  integer :: NMAX

  ! momentum diagnostics
  real :: tau_em, tau_sol, tau_com, tau_visc, tau_gyro, tau_parvisc

  logical :: ifirsttime = .true.

  ! xray diagnostic
  real :: xray_signal

  ! timing diagnostics
  real :: t_ludefall, t_sources, t_smoother, t_aux, t_onestep
  real :: t_solve_v, t_solve_n, t_solve_p, t_solve_b, t_mvm
  real :: t_output_cgm, t_output_hdf5, t_output_reset
  real :: t_gs

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
    efluxp = 0.
    efluxk = 0.
    efluxs = 0.
    efluxt = 0.
    epotg = 0.

    ptot = 0.

    area = 0.
    volume = 0.
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

    tau_em = 0.
    tau_sol = 0.
    tau_com = 0.
    tau_visc = 0.
    tau_gyro = 0.
    tau_parvisc = 0.

    nfluxd = 0.
    nfluxv = 0.
    nsource = 0.

    bwb2 = 0.

    xray_signal = 0.
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

    integer, parameter :: num_scalars = 47
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
       temp(35) = pvor
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
       efluxp = temp2(19)
       efluxk = temp2(20)
       efluxs = temp2(21)
       efluxt = temp2(22)
       epotg =  temp2(23)
       area =   temp2(24)
       totcur = temp2(25)
       totden = temp2(26)
       tflux =  temp2(27)
       tmom =   temp2(28)
       tvor =   temp2(29)
       parea =  temp2(30)
       pcur  =  temp2(31)
       pflux =  temp2(32)
       pden =   temp2(33)
       pmom =   temp2(34)
       pvor =   temp2(35)
       nfluxd = temp2(36)
       nfluxv = temp2(37)
       nsource= temp2(38)
       tau_em  =temp2(39)
       tau_sol =temp2(40)
       tau_com =temp2(41)
       tau_visc=temp2(42)
       tau_gyro=temp2(43)
       tau_parvisc=temp2(44)
       bwb2    =temp2(45)
       volume  =temp2(46)
       xray_signal=temp2(47)

    endif !if maxrank .gt. 1

  end subroutine distribute_scalars


!============================================================
! evaluate
! ~~~~~~~~
! calculates the value ans of field dum at global coordinates
! (x,z).  itri is the element containing (x,z).  (If this
! element does not reside on this process, itri=-1).
!============================================================
subroutine evaluate(x,phi,z,ans,ans2,fin,itri,ierr)
  
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

  real, intent(out) :: ans, ans2

  type(element_data) :: d
  integer :: p, nodeids(nodes_per_element), ier
  real :: x1, phi1, z1
  vectype, dimension(coeffs_per_element) :: avector
  real :: ri, si, zi, eta
  real :: term1, term2
  real, dimension(2) :: temp1, temp2
  integer :: hasval, tothasval

  ! evaluate the solution to get the value [ans] at one point (x,z)

  ! first find out what triangle x,z is in.  whattri
  ! returns itri, x1, and z1 with x1 and z1 being
  ! the coordinates of the first node/vertex

  if(itri.eq.0) then
     call whattri(x,phi,z,itri,x1,z1)
  else
     call get_element_nodes(itri,nodeids)
     call get_node_pos(nodeids(1), x1, phi1, z1)
  endif

  ans = 0.
  ans2 = 0.


  ! if this process contains the point, evaluate the field at that point.
  if(itri.gt.0) then

     call get_element_data(itri, d)

     ! calculate local coordinates
     call global_to_local(d, x, phi, z, si, zi, eta)

     ! calculate the inverse radius
     if(itor.eq.1) then
        ri = 1./x
     else
        ri = 1.
     endif

     ! calculate the value of the function
     call calcavector(itri, fin, avector)
     
     do p=1,20
     
        term1 = si**mi(p)*eta**ni(p)
        term2 = 0.
        
        if(mi(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 - 2.*d%co*(mi(p)*si**(mi(p)-1) * eta**ni(p))*ri
           endif
           
           if(mi(p).ge.2) then
              term2 = term2 + si**(mi(p)-2)*(mi(p)-1)*mi(p) * eta**ni(p)
           endif
        endif
     
        if(ni(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 + 2.*d%sn*(si**mi(p) * eta**(ni(p)-1)*ni(p))*ri
           endif
           
           if(ni(p).ge.2) then
              term2 = term2 + si**mi(p) * eta**(ni(p)-2)*(ni(p)-1)*ni(p)
           endif
        endif
     
        ans = ans + avector(p)*term1
        ans2 = ans2 + avector(p)*term2
        hasval = 1
     enddo
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
     temp1(1) = ans
     temp1(2) = ans2
     call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_WORLD, ier)
     ans = temp2(1)/tothasval
     ans2 = temp2(2)/tothasval
  endif

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
  
    implicit none

    real :: alx, alz, xrel, zrel
    real, dimension(2) :: temp, temp2
    integer :: itri, ierr

    call get_bounding_box_size(alx,alz)


    itri = 0
    xrel = xzero
    zrel = alz/2. + zzero
    call evaluate(xrel,0.,zrel,temp(1),temp2(1),psi_field(1),itri,ierr)

    itri = 0
    xrel = alx/2. + xzero
    zrel = alz/2. + zzero
    call evaluate(xrel,0.,zrel,temp(2),temp2(2),psi_field(1),itri,ierr)

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
    real :: tcpu
  intrinsic cpu_time
    call cpu_time(tcpu)
    return
  end subroutine second


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

  implicit none
 
  include 'mpif.h'

  integer :: itri, numelms, def_fields
  integer :: is_edge(3)  ! is inode on boundary
  real :: n(2,3),tpifac
  integer :: iedge, idim(3)
  real :: delta_t

  tpifac = 1.
  if(nplanes.gt.1) tpifac = twopi

  print *, 'Setting ptoto!', ptot
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
     def_fields = FIELD_PHI
     if(numvar.ge.2) def_fields = def_fields + FIELD_V
     if(numvar.ge.3) def_fields = def_fields + FIELD_CHI
  else
     def_fields = FIELD_PSI + FIELD_PHI + FIELD_ETA + FIELD_MU &
          + FIELD_N + FIELD_NI + FIELD_SIG
     if(numvar.ge.2) def_fields = def_fields + FIELD_I + FIELD_V
     if(numvar.ge.3) then
        def_fields = def_fields + FIELD_CHI + FIELD_P + FIELD_KAP
     endif

     if(gyro.eq.1 .or. amupar.ne.0) then
        def_fields = def_fields + FIELD_B2I
     endif

     if(numvar.ge.3 .or. ipres.eq.1) then
        if(hyper.eq.0.) def_fields = def_fields + FIELD_J
        if(hyperc.ne.0.) def_fields = def_fields + FIELD_VOR + FIELD_COM
     end if
  endif

  tm79 = 0.
  tm79(:,OP_1) = 1.

  call finalize(field0_vec)
  call finalize(field_vec)

  numelms = local_elements()
  do itri=1,numelms

     dbf = db

     call define_element_quadrature(itri, int_pts_diag, int_pts_tor)
     call define_fields(itri, def_fields, 0, 0)
     if(gyro.eq.1) call gyro_common

     ! Calculate energy
     ! ~~~~~~~~~~~~~~~~
     ekinp  = ekinp  + energy_kp ()
     ekinpd = ekinpd + energy_kpd()
     ekinph = ekinph + energy_kph()

     ekint  = ekint  + energy_kt ()
     ekintd = ekintd + energy_ktd()
     ekinth = ekinth + energy_kth()

     ekin3  = ekin3  + energy_k3 ()
     ekin3d = ekin3d + energy_k3d()
     ekin3h = ekin3h + energy_k3h()

     if(ike_only.eq.1) cycle

     emagp  = emagp  + energy_mp ()
     emagpd = emagpd + energy_mpd()
     emagph = emagph - qpsipsieta(tm79)

     emagt  = emagt  + energy_mt ()
     emagtd = emagtd + energy_mtd()
     emagth = emagth - qbbeta(tm79)

     emag3 = emag3 + energy_p()


     ! Calculate Scalars
     ! ~~~~~~~~~~~~~~~~~
     ! extra factor of 1/r comes from delta function in toroidal coordinate)
     area   = area   + int1(ri_79)/tpifac

     ! toroidal current
     totcur = totcur - int2(ri2_79,pst79(:,OP_GS))/tpifac

     ! toroidal flux
     tflux = tflux + int2(ri2_79,bzt79(:,OP_1 ))/tpifac
     
     ! enstrophy
     select case(ivform)
     case(0)
        tvor   = tvor   - int2(ri2_79,pht79(:,OP_GS))/tpifac
     case(1)
        tvor   = tvor &
             -    (int1(pht79(:,OP_LP)) &
             + 2.*int2(ri4_79,cht79(:,OP_DZ)))/tpifac
     end select

     ! volume
     volume = volume + int0()/tpifac

     ! particle number
     totden = totden + int1(nt79(:,OP_1))/tpifac

     ! particle source
     if(idens.eq.1) then        
        nsource = nsource - int1(sig79)/tpifac
     endif

     ! gravitational potential energy
     epotg = epotg + grav_pot()

     ! toroidal (angular) momentum
     if(numvar.ge.2) then
        select case(ivform)
        case(0)
           tmom = tmom &
                + int2(vzt79(:,OP_1),nt79(:,OP_1))
!           pmom = pmom &
!                + int3(vzt79(:,OP_1),nt79(:,OP_1),temp79a)
        case(1)
           tmom = tmom &
                + int3(r2_79,vzt79(:,OP_1),nt79(:,OP_1))
!           pmom = pmom &
!                + int4(r2_79,vzt79(:,OP_1),nt79(:,OP_1),temp79a)
        end select
     endif


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
     call boundary_edge(itri, is_edge, n, idim)

     do iedge=1,3
        if(is_edge(iedge).eq.0) cycle

        call define_boundary_quadrature(itri, iedge, 5, 5, n, idim)
        call define_fields(itri, def_fields, 1, 0)
        if(gyro.eq.1) call gyro_common

        ! Energy fluxes
        ! ~~~~~~~~~~~~~
        efluxp = efluxp + flux_pressure()
        efluxt = efluxt + flux_heat()
        efluxs = efluxs + flux_poynting()
        efluxk = efluxk + flux_ke()

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
           nfluxd = nfluxd - denm* &
                (int2(norm79(:,1),nt79(:,OP_DR)) &
                +int2(norm79(:,2),nt79(:,OP_DZ)))

           select case(ivform)
           case(0)
              nfluxv = nfluxv &
                   + int4(ri_79,nt79(:,OP_1),norm79(:,2),pht79(:,OP_DR)) &
                   - int4(ri_79,nt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ))
              
              if(numvar.ge.3) then
                 nfluxv = nfluxv &
                      + int3(nt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
                      + int3(nt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ))
              endif
           case(1)
              nfluxv = nfluxv &
                   + int4(r_79,nt79(:,OP_1),norm79(:,2),pht79(:,OP_DR)) &
                   - int4(r_79,nt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ))
              
              if(numvar.ge.3) then
                 nfluxv = nfluxv &
                      + int4(ri2_79,nt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
                      + int4(ri2_79,nt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ))
              endif
           end select
        end if

        ! xray signal
        temp79b = bremsstrahlung(nt79(:,OP_1), pet79(:,OP_1))
        call get_chord_mask(xray_r0, xray_phi0*pi/180., xray_z0, &
             x_79, phi_79, z_79, npoints, &
             xray_theta*pi/180., xray_sigma*pi/180., temp79a)
        xray_signal = xray_signal + int2(temp79a, temp79b)
     end do
  end do

  call distribute_scalars

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
    avep = (gam - 1.)*(emag3 / (volume*tpifac))
!
#ifdef USE3D
  if(ike_harmonics .gt. 0) call calculate_ke()
#endif

  if(myrank.eq.0 .and. iprint.ge.1) then 
     print *, "Total energy = ", etot
     print *, "Total energy lost = ", ptot
 
     print *, "Scalars:"
     print *, "  Area = ", area
     print *, "  Toroidal current = ", totcur
     print *, "  Toroidal flux = ", tflux
     print *, "  Volume = ", volume
     print *, "  Total particles = ", totden
  endif

end subroutine calculate_scalars


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
  integer, save :: itri = 0

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,2E12.4)') '  magaxis: guess = ', xguess, zguess

  converged = 0

  x = xguess
  z = zguess
  
  newton :  do inews=1, iterations

     call whattri(x,0.,z,itri,x1,z1)

     ! calculate position of minimum
     if(itri.gt.0) then
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
        if(iprint.ge.2) then
           write(*,'(A,4E12.4)') &
                '  magaxis: rdiff/h, tol, xnew,znew', rdiff/h, tol, xnew, znew
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
     
     ! communicate new minimum to all processors
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

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,I12,2E12.4)') '  magaxis: iterations, x, z = ', inews, x, z
  
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
  integer, save :: itri = 0

  if(myrank.eq.0 .and. iprint.ge.2) &
       write(*,'(A,2E12.4)') '  te_max: guess = ', xguess, zguess

  converged = 0

  x = xguess
  z = zguess
  
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
  
end subroutine te_max



!=====================================================
! lcfs
!
! locates the magnetic axis and the value of psi there
!=====================================================
subroutine lcfs(psi)
  use arrays
  use basic
  use mesh_mod
  use m3dc1_nint
  use field
  use boundary_conditions

  implicit none

  include 'mpif.h'

  type(field_type), intent(in) :: psi

  type(field_type) :: temp_field
  real :: psix, psib, psim
  real :: x, z, temp1, temp2, temp_min, temp_max, ajlim
  integer :: ier, numnodes, inode, izone, izonedim, itri
  logical :: is_boundary, first_point
  real, dimension(2) :: normal
  real :: curv
  vectype, dimension(dofs_per_node) :: data

  call create_field(temp_field)
  temp_field = psi
  if(icsubtract.eq.1) call add(temp_field, psi_coil_field)

  ! Find magnetic axis
  ! ~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Finding magnetic axis'
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
  do inode=1, numnodes
     call boundary_node(inode,is_boundary,izone,izonedim,normal,curv,x,z, &
          inner_wall)
     if(.not.is_boundary) cycle

     call get_node_data(temp_field,inode,data)

     if(((x-xmag)*real(data(2)) + (z-zmag)*real(data(3))) &
          *(real(data(1))-psimin).lt.0.) cycle

     if(first_point) then
        psib = real(data(1))
        first_point =.false.
     else
        if(abs(real(data(1)) - psimin).lt.abs(psib - psimin)) &
             psib = real(data(1))
     endif
  end do

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


  ! Calculate psi at the x-point
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Finding X-point'
  call magaxis(xnull,znull,temp_field,psix,1,ier)
  if(ier.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) then
        write(*,'(A,2E12.4)') '  X-point found at ', xnull, znull
     end if
  else 
     psix = psib

     if(myrank.eq.0 .and. iprint.ge.1) then 
        write(*,'(A,2E12.4)') '  no X-point found near ', xnull, znull
     end if
  endif


  if(abs(psix - psimin).lt.abs(psib - psimin)) then
     is_diverted = .true.
     psibound = psix
  else
     is_diverted = .true.
     psibound = psib
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
     call evaluate(xlim,0.,zlim,psilim,ajlim,temp_field,itri,ier)
     if(ier.ne.0) psilim = psibound
     
     ! calculate psi at a second limiter point as a diagnostic
     if(xlim2.gt.0) then
        itri = 0
        call evaluate(xlim2,0.,zlim2,psilim2,ajlim,temp_field,itri,ier)
        if(ier.ne.0) psilim = psibound
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

  ! daignostic output
  if(myrank.eq.0 .and. iprint.ge.1) then
     if(iprint.ge.1) then
        write(*,'(1A10,6A11)') 'psi at:', &
             'axis', 'wall', 'divertor', 'lim1', 'lim2', 'lcfs'
        write(*,'(1A10,1p6e11.3)') '',  &
             psimin, psib, psix, psilim, psilim2, psibound
     endif

     if(psibound.eq.psib) then
        print *, ' Plasma is limited by wall'
     else if(psibound.eq.psix) then
        print *, ' Plasma is diverted'
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
  real:: ke_N, ketotal, fac, delta_phi
  integer :: ier, k, l, numnodes, N
  vectype, dimension(dofs_per_node) :: vec_l

  real, allocatable :: i1ck(:,:), i1sk(:,:)
  real, allocatable :: i2ck(:,:), i2sk(:,:)

!  type(vector_type) :: transform_field
  type(field_type) :: u_transformc, vz_transformc, chi_transformc
  type(field_type) :: u_transforms, vz_transforms, chi_transforms


     NMAX = ike_harmonics
     numnodes = owned_nodes()
     if(.not.allocated(keharmonic)) allocate(keharmonic(0:NMAX))
     allocate(i1ck(nplanes,0:NMAX))
     allocate(i1sk(nplanes,0:NMAX))
     allocate(i2ck(nplanes,0:NMAX))
     allocate(i2sk(nplanes,0:NMAX))

!    create the sin and cos arrays
     do N = 0, NMAX
        do k = 1, nplanes
       call ke_I1(nplanes, NMAX, k, N, i1ck(k,N), i1sk(k,N))
       call ke_I2(nplanes, NMAX, k, N, i2ck(k,N), i2sk(k,N))
        enddo
     enddo

    if(myrank.eq.0 .and. iprint.eq.1) then
       write(*,900) ntime, numnodes, NMAX
       900 format("calculate_ke called-1,   ntime  numnodes  NMAX=",3i6)
     endif

!    call create_vector(transform_field,6)
!   
     call create_field(  u_transformc)
     call create_field(  u_transforms)
     call create_field( vz_transformc)
     call create_field( vz_transforms)
     call create_field(chi_transformc)
     call create_field(chi_transforms)

!     call associate_field(  u_transformc,transform_field,1)
!     call associate_field(  u_transforms,transform_field,2)
!     call associate_field( vz_transformc,transform_field,3)
!     call associate_field( vz_transforms,transform_field,4)
!     call associate_field(chi_transformc,transform_field,5)
!     call associate_field(chi_transforms,transform_field,6)

     if(myrank.eq.0 .and. iprint.eq.1) then
       write(*,901) ntime
       901 format("calculate_ke called-2,   ntime=",i6)
     endif


! for each Fourier mode
  do N=0,NMAX

  fac = 2.
  If (N.eq.0) fac = 1.

        k = local_plane() + 1

!test   delta_phi = 2. * pi / nplanes
!test1  u_field(1)=10.
!test2  u_field(1)= cos(2. * k * delta_phi)
!test3  u_field(1)= sin(2. * k * delta_phi)
        !eq 12: U cos
        do l =1,numnodes
           call get_node_data(u_field(1), l , u1_l ) ! u1_l is “U” (dimension 12)

           vec_l(1)= fac*(u1_l(1) * i1ck(k,N) + u1_l( 7)*i2ck(k,N))
           vec_l(2)= fac*(u1_l(2) * i1ck(k,N) + u1_l( 8)*i2ck(k,N))
           vec_l(3)= fac*(u1_l(3) * i1ck(k,N) + u1_l( 9)*i2ck(k,N))
           vec_l(4)= fac*(u1_l(4) * i1ck(k,N) + u1_l(10)*i2ck(k,N))
           vec_l(5)= fac*(u1_l(5) * i1ck(k,N) + u1_l(11)*i2ck(k,N))
           vec_l(6)= fac*(u1_l(6) * i1ck(k,N) + u1_l(12)*i2ck(k,N))
           vec_l(7:12) = 0. ! pad with zeros

           call set_node_data(u_transformc,l,vec_l)
        enddo
        call finalize(u_transformc%vec)
        call sum_vec_planes(u_transformc%vec%data) ! sum vec%datator at each (R,Z) node over k

!
!test1  u_field(1)=10.
!test2  u_field(1)= cos(2. * k * delta_phi)
!test3  u_field(1)= sin(2. * k * delta_phi)
        !eq 12: U sin
        do l =1,numnodes
           call get_node_data(u_field(1), l , u1_l ) ! u1_l is “U” (dimension 12)

           fac = 2.
           If (N.eq.0) fac = 1.
           vec_l(1)= fac*(u1_l(1) * i1sk(k,N) + u1_l( 7)*i2sk(k,N))
           vec_l(2)= fac*(u1_l(2) * i1sk(k,N) + u1_l( 8)*i2sk(k,N))
           vec_l(3)= fac*(u1_l(3) * i1sk(k,N) + u1_l( 9)*i2sk(k,N))
           vec_l(4)= fac*(u1_l(4) * i1sk(k,N) + u1_l(10)*i2sk(k,N))
           vec_l(5)= fac*(u1_l(5) * i1sk(k,N) + u1_l(11)*i2sk(k,N))
           vec_l(6)= fac*(u1_l(6) * i1sk(k,N) + u1_l(12)*i2sk(k,N))
           vec_l(7:12) = 0. ! pad with zeros
           call set_node_data(u_transforms,l,vec_l)
        enddo
        call finalize(u_transforms%vec)
        call sum_vec_planes(u_transforms%vec%data) ! sum vec%datator at each (R,Z) node over k

!
!test1  vz_field(1)=10.
!test2  vz_field(1)= cos(2. * k * delta_phi)
!test3  vz_field(1)= sin(2. * k * delta_phi)
        !eq 12: omega cos
        do l =1,numnodes
           call get_node_data(vz_field(1), l , vz1_l) ! vz1_l is “ω” ( dimension 12)

           vec_l(1)= fac*(vz1_l(1) * i1ck(k,N) + vz1_l( 7)*i2ck(k,N))
           vec_l(2)= fac*(vz1_l(2) * i1ck(k,N) + vz1_l( 8)*i2ck(k,N))
           vec_l(3)= fac*(vz1_l(3) * i1ck(k,N) + vz1_l( 9)*i2ck(k,N))
           vec_l(4)= fac*(vz1_l(4) * i1ck(k,N) + vz1_l(10)*i2ck(k,N))
           vec_l(5)= fac*(vz1_l(5) * i1ck(k,N) + vz1_l(11)*i2ck(k,N))
           vec_l(6)= fac*(vz1_l(6) * i1ck(k,N) + vz1_l(12)*i2ck(k,N))
           vec_l(7:12) = 0. ! pad with zeros
          call set_node_data(vz_transformc,l,vec_l)
        enddo
        call finalize(vz_transformc%vec)
        call sum_vec_planes(vz_transformc%vec%data) ! sum vec%datator at each (R,Z) node over k

!
!test1  vz_field(1)=10.
!test2  vz_field(1)= cos(2. * k * delta_phi)
!test3  vz_field(1)= sin(2. * k * delta_phi)
        !eq 12: omega sin
        do l =1,numnodes
           call get_node_data(vz_field(1), l , vz1_l) ! vz1_l is “ω” ( dimension 12)

           vec_l(1)= fac*(vz1_l(1) * i1sk(k,N) + vz1_l( 7)*i2sk(k,N))
           vec_l(2)= fac*(vz1_l(2) * i1sk(k,N) + vz1_l( 8)*i2sk(k,N))
           vec_l(3)= fac*(vz1_l(3) * i1sk(k,N) + vz1_l( 9)*i2sk(k,N))
           vec_l(4)= fac*(vz1_l(4) * i1sk(k,N) + vz1_l(10)*i2sk(k,N))
           vec_l(5)= fac*(vz1_l(5) * i1sk(k,N) + vz1_l(11)*i2sk(k,N))
           vec_l(6)= fac*(vz1_l(6) * i1sk(k,N) + vz1_l(12)*i2sk(k,N))
           vec_l(7:12) = 0. ! pad with zeros
           call set_node_data(vz_transforms,l,vec_l)
        enddo
        call finalize(vz_transforms%vec)
        call sum_vec_planes(vz_transforms%vec%data) ! sum vec%datator at each (R,Z) node over k

!
!test1  chi_field(1)=10.
!test2  chi_field(1)= cos(2. * k * delta_phi)
!test3  chi_field(1)= sin(2. * k * delta_phi)
        !eq 12: chi cos
        do l =1,numnodes
           call get_node_data(chi_field(1), l , chi1_l ) ! chi1_l is “χ” (dimension 12)

           fac = 2.
           If (N.eq.0) fac = 1.
           vec_l(1)= fac*(chi1_l(1) * i1ck(k,N) + chi1_l( 7)*i2ck(k,N))
           vec_l(2)= fac*(chi1_l(2) * i1ck(k,N) + chi1_l( 8)*i2ck(k,N))
           vec_l(3)= fac*(chi1_l(3) * i1ck(k,N) + chi1_l( 9)*i2ck(k,N))
           vec_l(4)= fac*(chi1_l(4) * i1ck(k,N) + chi1_l(10)*i2ck(k,N))
           vec_l(5)= fac*(chi1_l(5) * i1ck(k,N) + chi1_l(11)*i2ck(k,N))
           vec_l(6)= fac*(chi1_l(6) * i1ck(k,N) + chi1_l(12)*i2ck(k,N))
           vec_l(7:12) = 0. ! pad with zeros
           call set_node_data(chi_transformc,l,vec_l)
        enddo
        call finalize(chi_transformc%vec)
        call sum_vec_planes(chi_transformc%vec%data) ! sum vec%datator of size 6 at each (R,Z) node over k

!
!test1  chi_field(1)=10.
!test2  chi_field(1)= cos(2. * k * delta_phi)
!test3  chi_field(1)= sin(2. * k * delta_phi)
        ! eq 12: chi sin
        do l =1,numnodes
           call get_node_data(chi_field(1), l , chi1_l ) ! chi1_l is “χ” (dimension 12)

           vec_l(1)= fac*(chi1_l(1) * i1sk(k,N) + chi1_l( 7)*i2sk(k,N))
           vec_l(2)= fac*(chi1_l(2) * i1sk(k,N) + chi1_l( 8)*i2sk(k,N))
           vec_l(3)= fac*(chi1_l(3) * i1sk(k,N) + chi1_l( 9)*i2sk(k,N))
           vec_l(4)= fac*(chi1_l(4) * i1sk(k,N) + chi1_l(10)*i2sk(k,N))
           vec_l(5)= fac*(chi1_l(5) * i1sk(k,N) + chi1_l(11)*i2sk(k,N))
           vec_l(6)= fac*(chi1_l(6) * i1sk(k,N) + chi1_l(12)*i2sk(k,N))
           vec_l(7:12) = 0. ! pad with zeros
           call set_node_data(chi_transforms,l,vec_l)
        enddo
        call finalize(chi_transforms%vec)
        call sum_vec_planes(chi_transforms%vec%data) ! sum vec%datator at each (R,Z) node over k

!    enddo
!    call finalize(transform_field)


!eq 4b: Calculate energy for each Fourier Harminics N
   
!
  ke_N = 0.
  def_fields = 0
     numelms = local_elements()
     do itri=1,numelms

        call define_element_quadrature(itri, int_pts_diag, int_pts_tor)
        call define_fields(itri, def_fields, 1, 0)
!
!       cosine harmonics
        call eval_ops(itri,  u_transformc,pht79)
        call eval_ops(itri, vz_transformc,vzt79)
        call eval_ops(itri,chi_transformc,cht79)
!
        ke_N = ke_N + int3(r2_79,  pht79(:,OP_DR), pht79(:,OP_DR))   &
                    + int3(r2_79,  pht79(:,OP_DZ), pht79(:,OP_DZ))

        ke_N = ke_N + int3(r2_79,  vzt79(:,OP_1), vzt79(:,OP_1))

        ke_N = ke_N + int3(ri4_79,  cht79(:,OP_DR), cht79(:,OP_DR))   &
                    + int3(ri4_79,  cht79(:,OP_DZ), cht79(:,OP_DZ))
!
!       sine harmonics
        call eval_ops(itri,  u_transforms,pht79)
        call eval_ops(itri, vz_transforms,vzt79)
        call eval_ops(itri,chi_transforms,cht79)
!
        ke_N = ke_N + int3(r2_79,  pht79(:,OP_DR), pht79(:,OP_DR))   &
                    + int3(r2_79,  pht79(:,OP_DZ), pht79(:,OP_DZ))

        ke_N = ke_N + int3(r2_79,  vzt79(:,OP_1), vzt79(:,OP_1))

        ke_N = ke_N + int3(ri4_79,  cht79(:,OP_DR), cht79(:,OP_DR))   &
                    + int3(ri4_79,  cht79(:,OP_DZ), cht79(:,OP_DZ))
     end do

     call mpi_allreduce(ke_N, ketotal, 1, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, mpi_comm_world, ier)

     keharmonic(N) = ketotal / 4.
  end do

!!!!!....we need to save keharmonic for output <===
  if(myrank.eq.0 .and. iprint.ge.1) then
      write(*,1001) ntime
      write(*,1002) (keharmonic(N),N=0,NMAX)
 1001 format(" keharmonics at cycle",i6)
 1002 format(1p5e12.4)
  endif

! deallocate(keharmonic)
  deallocate(i1ck, i1sk, i2ck, i2sk)
!   
     call destroy_field(  u_transformc)
     call destroy_field(  u_transforms)
     call destroy_field( vz_transformc)
     call destroy_field( vz_transforms)
     call destroy_field(chi_transformc)
     call destroy_field(chi_transforms)
!    call destroy_vector(transform_field)

end subroutine calculate_ke



!...........................................................
! input: k, N
! output: i1ck, i1sk
! eq 10
  subroutine ke_I1(nplanes, NMAX, k, N, i1ck, i1sk)
    use math
    implicit none
    integer:: k, N, nplanes, NMAX
    real:: i1ck, i1sk

    real:: delta_phi

#ifdef OLD
    real:: I1_cNk_minus, I1_cNk_plus, I1_sNk_minus, I1_sNk_plus
    real:: Phi_1_plus, Phi_1_zero, Phi_1_minus, &
           Phi_1p_plus, Phi_1p_zero, Phi_1p_minus, &
           Phi_1pp_plus, Phi_1pp_zero, Phi_1pp_minus, &
           Phi_1ppp_plus, Phi_1ppp_zero_up, Phi_1ppp_zero_dn, Phi_1ppp_minus
    real:: sin_plus, sin_zero, sin_minus, &
           cos_plus, cos_zero, cos_minus 

#else
!  use 3/8 Simpson's rule
    real:: hh, x1, x2, x3, x4, f1, f2, f3, f4
    real:: phi1, phi2, phi3, phi4
#endif


   delta_phi = 2. * pi / nplanes

! $ \Phi_1(x) = ( |x|-1)^2( 2|x|+1), |x| \leq 1 $
   if(N .le. 0) then
      i1ck = delta_phi
      i1sk = 0.
      return
   endif

#ifdef OLD
   Phi_1_plus  =0. ! x=1 \\
   Phi_1_zero  =1. ! x=0 \\
   Phi_1_minus =0. ! x=-1

   Phi_1p_plus  =0. ! x=1 \\
   Phi_1p_zero  =0. ! x=0 \\
   Phi_1p_minus =0. ! x=-1

   Phi_1pp_plus  = 6.  ! x=1 \\
   Phi_1pp_zero  =-6.  ! x=0 \\
   Phi_1pp_minus = 6.  ! x=-1

   Phi_1ppp_plus    = 12.  ! x=1 \\
   Phi_1ppp_zero_up = 12.  ! x=0 \mbox{  in interval } [0,1]\\
   Phi_1ppp_zero_dn =-12.  ! x=0 \mbox{  in interval } [-1,0]\\
   Phi_1ppp_minus   =-12.  ! x=-1

   sin_plus  = N*( 1 + k-1) * delta_phi
   sin_zero  = N*( 0 + k-1) * delta_phi
   sin_minus = N*(-1 + k-1) * delta_phi
   cos_plus  = N*( 1 + k-1) * delta_phi
   cos_zero  = N*( 0 + k-1) * delta_phi
   cos_minus = N*(-1 + k-1) * delta_phi

   I1_cNk_minus =                                                                    &
     (                                                                               &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_1_zero * sin_zero                       &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_zero * cos_zero                      &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_zero * sin_zero                     &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_zero_dn * cos_zero                 &
     )                                                                               &
     -                                                                               &
     (                                                                               &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_1_minus * sin_minus                     &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_minus * cos_minus                    &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_minus * sin_minus                   &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_minus * cos_minus                  &
     )

   I1_cNk_plus =                                                                     &
     (                                                                               &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_1_plus * sin_plus                       &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_plus * cos_plus                      &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_plus * sin_plus                     &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_plus * cos_plus                    &
     )                                                                               &
     -                                                                               &
     (                                                                               &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_1_zero * sin_zero                       &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_zero * cos_zero                      &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_zero * sin_zero                     &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_zero_up * cos_zero                 &
     )

   i1ck =  I1_cNk_minus + I1_cNk_plus

#else
!  use 3/8 Simpson's rule
   hh = 2. / 3.
   x1 = -1.
   x2 = -1. + hh
   x3 =  1. - hh
   x4 =  1.
   phi1 = ( abs(x1) - 1. ) * ( abs(x1) - 1. ) * ( 2. * abs(x1) + 1. )
   phi2 = ( abs(x2) - 1. ) * ( abs(x2) - 1. ) * ( 2. * abs(x2) + 1. )
   phi3 = ( abs(x3) - 1. ) * ( abs(x3) - 1. ) * ( 2. * abs(x3) + 1. )
   phi4 = ( abs(x4) - 1. ) * ( abs(x4) - 1. ) * ( 2. * abs(x4) + 1. )
   f1 = cos(N*(x1+k-1)*delta_phi) * phi1
   f2 = cos(N*(x2+k-1)*delta_phi) * phi2
   f3 = cos(N*(x3+k-1)*delta_phi) * phi3
   f4 = cos(N*(x4+k-1)*delta_phi) * phi4
   i1ck =  hh*3./8. * ( f1 + 3.*f2 + 3.*f3 + f4)
   i1ck =  i1ck * delta_phi
#endif


#ifdef OLD
   I1_sNk_minus  =                                                                   &
     (                                                                               &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_1_zero * cos_zero                       &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_zero * sin_zero                      &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_zero * cos_zero                     &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_zero_dn * sin_zero                 &
     )                                                                               &
    -                                                                                &
     (                                                                               &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_1_minus * cos_minus                     &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_minus * sin_minus                    &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_minus * cos_minus                   &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_minus * sin_minus                  &
     )

   I1_sNk_plus =                                                                     &
     (                                                                               &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_1_plus * cos_plus                       &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_plus * sin_plus                      &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_plus * cos_plus                     &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_plus * sin_plus                    &
     )                                                                               &
     -                                                                               &
     (                                                                               &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_1_zero * cos_zero                       &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_1p_zero * sin_zero                      &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_1pp_zero * cos_zero                     &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_1ppp_zero_up * sin_zero                 &
     )

   i1sk =  I1_sNk_minus + I1_sNk_plus

#else
!  use 3/8 Simpson's rule
   f1 = sin(N*(x1+k-1)*delta_phi) * phi1
   f2 = sin(N*(x2+k-1)*delta_phi) * phi2
   f3 = sin(N*(x3+k-1)*delta_phi) * phi3
   f4 = sin(N*(x4+k-1)*delta_phi) * phi4
   i1sk =  hh*3./8. * ( f1 + 3.*f2 + 3.*f3 + f4)
   i1sk =  i1sk * delta_phi
#endif


end subroutine ke_I1


!...........................................................
! input: N, k, delta_phi
! output: i2ck, i2sk
! eq 10
  subroutine ke_I2(nplanes, NMAX, k, N, i2ck, i2sk)
    use math
    implicit none
    integer:: k, N, nplanes, NMAX
    real:: i2ck, i2sk

    real:: delta_phi

#ifdef OLD
    real:: I2_cNk_minus, I2_cNk_plus, I2_sNk_minus, I2_sNk_plus
    real:: Phi_2_plus, Phi_2_zero, Phi_2_minus, &
           Phi_2p_plus, Phi_2p_zero, Phi_2p_minus, &
           Phi_2pp_plus, Phi_2pp_zero_up, Phi_2pp_zero_dn, Phi_2pp_minus, &
           Phi_2ppp_plus, Phi_2ppp_zero, Phi_2ppp_minus
    real:: sin_plus, sin_zero, sin_minus, &
           cos_plus, cos_zero, cos_minus

#else
!  use 3/8 Simpson's rule
    real:: hh, x1, x2, x3, x4, f1, f2, f3, f4
    real:: phi1, phi2, phi3, phi4
#endif

   delta_phi = 2. * pi / nplanes

!  $ \Phi_2(x) = x ( |x|-1)^2, |x| \leq 1 $
   if(N .le. 0) then
      i2ck = 0.
      i2sk = 0.
      return
   endif

#ifdef OLD
   Phi_2_plus  =0. ! x=1 \\
   Phi_2_zero  =0. ! x=0 \\
   Phi_2_minus =0. ! x=-1 

   Phi_2p_plus  =0. ! x=1 \\
   Phi_2p_zero  =1. ! x=0 \\
   Phi_2p_minus =0. ! x=-1 

   Phi_2pp_plus    =  2. ! x=1 \\
   Phi_2pp_zero_up = -4. ! x=0 \mbox{  in interval  } [0, 1]\\
   Phi_2pp_zero_dn =  4. ! x=-1 \mbox{  in interval  } [-1, 0] \\
   Phi_2pp_minus   = -2. ! x=-1

   Phi_2ppp_plus  = 6.  ! x=1 \\
   Phi_2ppp_zero  = 6.  ! x=0 \\
   Phi_2ppp_minus = 6.  ! x=-1 \\

   sin_plus  = N*( 1 + k-1) * delta_phi
   sin_zero  = N*( 0 + k-1) * delta_phi
   sin_minus = N*(-1 + k-1) * delta_phi
   cos_plus  = N*( 1 + k-1) * delta_phi
   cos_zero  = N*( 0 + k-1) * delta_phi
   cos_minus = N*(-1 + k-1) * delta_phi

   I2_cNk_minus =                                                            &
     (                                                                       &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_2_zero * sin_zero               &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_zero * cos_zero              &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_zero_dn * sin_zero          &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_zero * cos_zero            &
     )                                                                       &
     -                                                                       &
     (                                                                       &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_2_minus * sin_minus             &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_minus * cos_minus            &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_minus * sin_minus           &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_minus * cos_minus          &
     )          

   I2_cNk_plus =                                                             & 
     (                                                                       &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_2_plus * sin_plus               &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_plus * cos_plus              &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_plus * sin_plus             &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_plus * cos_plus            &
     )                                                                       &
     -                                                                       &
     (                                                                       &
     + delta_phi * 1./(N*delta_phi)**1 * Phi_2_zero * sin_zero               &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_zero * cos_zero              &
     - delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_zero_up * sin_zero          &  
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_zero * cos_zero            &
     )

   i2ck =  I2_cNk_minus + I2_cNk_plus

#else
!  use 3/8 Simpson's rule
   hh = 2. / 3.
   x1 = -1.
   x2 = -1. + hh
   x3 =  1. - hh
   x4 =  1.
   phi1 = x1 * ( abs(x1) - 1. ) * ( abs(x1) - 1. )
   phi2 = x2 * ( abs(x2) - 1. ) * ( abs(x2) - 1. )
   phi3 = x3 * ( abs(x3) - 1. ) * ( abs(x3) - 1. )
   phi4 = x4 * ( abs(x4) - 1. ) * ( abs(x4) - 1. )
   f1 = cos(N*(x1+k-1)*delta_phi) * phi1
   f2 = cos(N*(x2+k-1)*delta_phi) * phi2
   f3 = cos(N*(x3+k-1)*delta_phi) * phi3
   f4 = cos(N*(x4+k-1)*delta_phi) * phi4
   i2ck =  hh*3./8. * ( f1 + 3.*f2 + 3.*f3 + f4)
   i2ck =  i2ck * delta_phi
#endif


#ifdef OLD
   I2_sNk_minus  =                                                           &
     (                                                                       &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_2_zero * cos_zero               &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_zero * sin_zero              &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_zero_dn * cos_zero          &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_zero * sin_zero            &
     )                                                                       &
    -                                                                        &
     (                                                                       &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_2_minus * cos_minus             &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_minus * sin_minus            &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_minus * cos_minus           &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_minus * sin_minus          &
     )

   I2_sNk_plus =                                                             &
     (                                                                       &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_2_plus * cos_plus               &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_plus * sin_plus              &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_plus * cos_plus             &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_plus * sin_plus            &
     )                                                                       &
     -                                                                       &
     (                                                                       &
     - delta_phi * 1./(N*delta_phi)**1 * Phi_2_zero * cos_zero               &
     + delta_phi * 1./(N*delta_phi)**2 * Phi_2p_zero * sin_zero              &
     + delta_phi * 1./(N*delta_phi)**3 * Phi_2pp_zero_up * cos_zero          &
     - delta_phi * 1./(N*delta_phi)**4 * Phi_2ppp_zero * sin_zero            &
     )

   i2sk =  I2_sNk_minus + I2_sNk_plus

#else
!  use 3/8 Simpson's rule
   f1 = sin(N*(x1+k-1)*delta_phi) * phi1
   f2 = sin(N*(x2+k-1)*delta_phi) * phi2
   f3 = sin(N*(x3+k-1)*delta_phi) * phi3
   f4 = sin(N*(x4+k-1)*delta_phi) * phi4
   i2sk =  hh*3./8. * ( f1 + 3.*f2 + 3.*f3 + f4)
   i2sk =  i2sk * delta_phi
#endif


end subroutine ke_I2


end module diagnostics

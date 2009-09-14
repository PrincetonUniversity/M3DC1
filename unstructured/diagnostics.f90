module diagnostics

  implicit none

  real :: tflux0, totcur0

  ! scalars integrated over entire computational domain
  real :: tflux, area, totcur, totden, tmom, tvor, bwb2

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
       emag3h,emag3d,emag3ho,emag3do
  real :: efluxp,efluxk,efluxs,efluxt,epotg,etot,ptot,eerr,ptoto


  ! momentum diagnostics
  real :: tau_em, tau_sol, tau_com, tau_visc, tau_gyro, tau_parvisc

  logical :: ifirsttime = .true.

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

    integer, parameter :: num_scalars = 45
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

    endif !if maxrank .gt. 1

  end subroutine distribute_scalars


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
  
    implicit none

    real :: alx, alz, xrel, zrel
    real, dimension(2) :: temp, temp2
    integer :: itri

    call getboundingboxsize(alx,alz)


    itri = 0
    xrel = xzero
    zrel = alz/2. + zzero
    call evaluate(xrel,zrel,temp(1),temp2(1),field,psi_g,num_fields,itri)

    itri = 0
    xrel = alx/2. + xzero
    zrel = alz/2. + zzero
    call evaluate(xrel,zrel,temp(2),temp2(2),field,psi_g,num_fields,itri)

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
  use newvar_mod
  use sparse
  use metricterms_new

  implicit none
#include "finclude/petsc.h" 
 
  integer :: itri, numelms, i, ione, def_fields
  real :: x, z, xmin, zmin, factor
  logical :: is_edge(3)  ! is inode on boundary
  real :: n(2,3)
  integer :: iedge, idim(3)


  double precision, dimension(3)  :: cogcoords

  integer :: ier
  PetscTruth :: flg_petsc, flg_solve2, flg_solve1
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier)

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
     def_fields = FIELD_PSI + FIELD_PHI + FIELD_J + FIELD_ETA + FIELD_MU &
          + FIELD_N + FIELD_NI + FIELD_SIG
     if(numvar.ge.2) def_fields = def_fields + FIELD_I + FIELD_V
     if(numvar.ge.3) then
        def_fields = def_fields + FIELD_CHI + &
             FIELD_PE + FIELD_P + FIELD_KAP
     endif
     if((numvar.ge.3. .and. kappar.ne.0.) &
          .or. amupar.ne.0.) def_fields = def_fields + FIELD_B2I
   
     if(hypc.ne.0.) then 
        def_fields = def_fields + FIELD_VOR
        if(numvar.ge.3) def_fields = def_fields + FIELD_COM
     end if
  endif

  tm79 = 0.
  tm79(:,OP_1) = 1.

  if(isources.eq.1) then
     sb1 = 0.
     sb2 = 0.
     sp1 = 0.
  end if

  call getmincoord(xmin,zmin)

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

     call define_triangle_quadrature(itri, int_pts_diag)
     call define_fields(itri, def_fields, isources, 0)


     ! Define Source terms
     ! ~~~~~~~~~~~~~~~~~~~
     if(isources.eq.1) then

        do i=1,18
           ione = isval1(itri,i)
          
           ! Definition of Source Terms
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
           sb1(ione) = sb1(ione) + b1psieta(g79(:,:,i),pst79,eta79,hypf*sz79)
           
           if(numvar.ge.2) then
              sb1(ione) = sb1(ione) + b1psibd(g79(:,:,i),pst79,bzt79,ni79)*dbf
              
              sb2(ione) = sb2(ione)  &
                   + b2psipsid(g79(:,:,i),pst79,pst79,ni79)*dbf &
                   + b2bbd    (g79(:,:,i),bzt79,bzt79,ni79)*dbf &
                   + b2beta   (g79(:,:,i),bzt79,eta79,hypi*sz79)
           endif
        
           if(numvar.ge.3) then
              sb2(ione) = sb2(ione) + b2ped(g79(:,:,i),pet79,ni79)*dbf*pefac
              
              sp1(ione) = sp1(ione) &
                   + b3psipsieta(g79(:,:,i),pst79,pst79,eta79)   &
                   + b3bbeta    (g79(:,:,i),bzt79,bzt79,eta79)   &
!                   + b3pedkappa (g79(:,:,i),pt79,ni79,kappat,hypp*sz79) &
!                   + p1kappar   (g79(:,:,i),pst79,pst79,pet79,ni79,b2i79,kar79) &
!                   + p1kappax   (g79(:,:,i),pet79,bzt79,ni79,b2i79,kar79) &
                   + b3pebd(g79(:,:,i),pet79,bzt79,ni79)*dbf*pefac
              
              ! ohmic heating
              sp1(ione) = sp1(ione) + (gam-1.)* &
                   (qpsipsieta(g79(:,:,i)) &
                   +qbbeta    (g79(:,:,i)))
           
              ! viscous heating
              sp1(ione) = sp1(ione) - (gam-1.)* &
                   (quumu    (g79(:,:,i),pht79,pht79,vis79,      hypc*sz79) &
                   +qvvmu    (g79(:,:,i),vzt79,vzt79,vis79,      hypv*sz79) &
                   +quchimu  (g79(:,:,i),pht79,cht79,vis79,vic79,hypc*sz79) &
                   +qchichimu(g79(:,:,i),cht79,cht79,      vic79,hypc*sz79))
           endif ! on numvar.ge.3
        end do

     endif ! on isources


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
     !  (extra factor of 1/r comes from delta function in toroidal coordinate)
     area   = area   + int1(ri_79)
     totcur = totcur - int2(ri2_79,pst79(:,OP_GS))
     
     ! enstrophy
     select case(ivform)
     case(0)
        tvor   = tvor   - int2(ri2_79,pht79(:,OP_GS))
     case(1)
        tvor   = tvor &
             -    int1(pht79(:,OP_LP)) &
             - 2.*int2(ri4_79,cht79(:,OP_DZ))
     end select

     ! toroidal flux
     tflux = tflux + int2(ri2_79,bzt79(:,OP_1 ))

     ! particle number
     totden = totden + int1(nt79(:,OP_1))

     ! particle source
     if(idens.eq.1) then        
        nsource = nsource - int1(sig79)
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
             3.*int3(vip79(:,OP_1),temp79a,temp79a)
     end if


     ! add surface terms
     call boundary_edge(itri, is_edge, n, idim)

     do iedge=1,3
        if(.not.is_edge(iedge)) cycle

        call define_edge_quadrature(itri, iedge, 5, n, idim)
        call define_fields(itri, def_fields, 1, 0)

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
     end do
  end do


  if(isources.eq.1) then
     call solve_newvar(sb1, NV_DCBOUND, mass_matrix_lhs_dc)
     if(numvar.ge.2) call solve_newvar(sb2, NV_DCBOUND, mass_matrix_lhs_dc)
     if(numvar.ge.3) call solve_newvar(sp1, NV_DCBOUND, mass_matrix_lhs_dc)
  endif

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

  if(myrank.eq.0 .and. iprint.ge.1) then 
     print *, "Total energy = ", etot
     print *, "Total energy lost = ", ptot
 
     print *, "Scalars:"
     print *, "  Area = ", area
     print *, "  Toroidal current = ", totcur
     print *, "  Total particles = ", totden
  endif

end subroutine calculate_scalars


!=====================================================
! magaxis
!
! locates the magnetic axis and the value of psi there
! imethod = 0 finds the local minimum of psi
! imethod = 1 finds the local zero of <psi,psi>
!=====================================================
subroutine magaxis(xguess,zguess,phin,iplace,numvari,psim,imethod,ier)
  use basic
  use t_data
  use nintegrate_mod

  implicit none

  include 'mpif.h'

  real, intent(inout) :: xguess, zguess
  vectype, intent(in), dimension(*) :: phin
  integer, intent(in) :: iplace, numvari, imethod
  real, intent(out) :: psim

  integer, parameter :: iterations = 20  !  max number of Newton iterations
  real, parameter :: bfac = 0.5  !max zone fraction for movement each iteration
  real, parameter :: tol = 1e-3   ! convergence tolorance (fraction of h)

  integer :: itri, inews
  integer :: i, ier, in_domain, converged
  real :: x1, z1, x, z, theta, a, b, c, co, sn, si, eta, h
  real :: sum, sum1, sum2, sum3, sum4, sum5
  real :: term1, term2, term3, term4, term5
  real :: pt, pt1, pt2, p11, p22, p12
  real :: xnew, znew, denom, sinew, etanew
  real :: xtry, ztry, rdiff
  vectype, dimension(20) :: avector
  real, dimension(5) :: temp1, temp2


  if(myrank.eq.0 .and. iprint.gt.0) &
       write(*,'(A,2E12.4)') '  magaxis: guess = ', xguess, zguess

  converged = 0

  x = xguess
  z = zguess
  
  newton :  do inews=1, iterations

     call whattri(x,z,itri,x1,z1)

     ! calculate position of minimum
     if(itri.gt.0) then
        call calcavector(itri, phin, iplace, numvari, avector)

        ! calculate local coordinates
        theta = ttri(itri)
        a = atri(itri)
        b = btri(itri)
        c = ctri(itri)
        h = sqrt((a+b)*c)
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
           sinew = si -  ( p22*pt1 - p12*pt2)/denom
           etanew= eta - (-p12*pt1 + p11*pt2)/denom

        case(1)  ! find local zero of <psi,psi>
           pt = sum1**2 + sum2**2
           pt1 = 2.*(sum1*sum3 + sum2*sum4)
           pt2 = 2.*(sum1*sum4 + sum2*sum5)

           denom = pt1**2 + pt2**2
           sinew = si - pt*pt1/denom
           etanew = eta - pt*pt2/denom
        end select

        xtry = x1 + co*(b+sinew) - sn*etanew
        ztry = z1 + sn*(b+sinew) + co*etanew

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
        if(iprint.ge.2) &
             write(*,'(A,4E12.4)') '  magaxis: rdiff/h, tol, xnew,znew', rdiff/h, tol, xnew,znew
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
             write(*,'(A,2E12.4)') '  magaxis: guess outside domain ', x, z
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

  if(myrank.eq.0 .and. iprint.ge.1) &
       write(*,'(A,I12,2E12.4)') '  magaxis: iterations, x, z = ', inews, x, z
  
end subroutine magaxis



!=====================================================
! lcfs
!
! locates the magnetic axis and the value of psi there
!=====================================================
subroutine lcfs(phin, iplace, numvari)
  use basic
  use t_data
  use nintegrate_mod

  implicit none

  include 'mpif.h'

  vectype, intent(in), dimension(*) :: phin
  integer, intent(in) :: iplace,numvari

  real :: psix, psib, psim
  real :: x, z, temp1, temp2, temp_min, temp_max, ajlim
  integer :: ier, numnodes, inode, izone, izonedim, itri
  integer :: ibegin, iendplusone, index
  logical :: is_boundary, first_point
  real, dimension(2) :: normal
  real :: curv

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Finding LCFS:'

  ! Find magnetic axis
  ! ~~~~~~~~~~~~~~~~~~
  call magaxis(xmag,zmag,phin,iplace,numvari,psim,0,ier)
  if(ier.eq.0) then
     psimin = psim
     
     if(myrank.eq.0 .and. iprint.ge.1) then 
        write(*,'(A,2E12.4)') '  magnetic axis found at ', xmag, zmag
        write(*,'(A, E12.4)') '  value of psi at magnetic axis ', psimin
     end if
  else
     if(myrank.eq.0 .and. iprint.ge.1) then 
        write(*,'(A,2E12.4)') '  no magnetic axis found near ', xmag, zmag
     end if
  endif

  if(ifixedb.eq.1) then
     psilim = 0.
     psilim2 = 0.
     psibound = 0.
     return
  endif

  ! Find the maximum value of psi at the boundary 
  ! such that psi is not increasing outward 
  ! (as in a private flux region)
  first_point = .true.
  call numnod(numnodes)
  do inode=1, numnodes
     call boundary_node(inode,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call entdofs(numvari,inode,0,ibegin,iendplusone)
     index = ibegin+(iplace-1)*6
     if(((x-xmag)*real(phin(index+1)) + &
          (z-zmag)*real(phin(index+2)))*(real(phin(index))-psimin).lt.0.) cycle

     if(first_point) then
        psib = real(phin(index))
        first_point =.false.
     else
        if(abs(real(phin(index)) - psimin).lt.abs(psib - psimin)) &
             psib = real(phin(index))
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
  call magaxis(xnull,znull,phin,iplace,numvari,psix,1,ier)
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
  if(xlim.eq.0.) then
     ! when xlim = 0, use the lcfs as the limiting flux
     psilim = psibound
     psilim2 = psilim
  else
     itri = 0
     call evaluate(xlim,zlim,psilim,ajlim,phin,iplace,numvari,itri)
     
     ! calculate psi at a second limiter point as a diagnostic
     if(xlim2.gt.0) then
        itri = 0
        call evaluate(xlim2,zlim2,psilim2,ajlim,phin,iplace,numvari,itri)
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
     write(*,'(1A10,6A11)') 'psi at:', &
          'axis', 'wall', 'divertor', 'lim1', 'lim2', 'lcfs'
     write(*,'(1A10,1p6e11.3)') '',  &
          psimin, psib, psix, psilim, psilim2, psibound

     if(psibound.eq.psib) then
        print *, 'Plasma is limited by wall'
     else if(psibound.eq.psix) then
        print *, 'Plasma is diverted'
     else if(psibound.eq.psilim) then
        print *, 'Plasma is limited by internal limiter #1.'
     else if(psibound.eq.psilim2) then
        print *, 'Plasma is limited by internal limiter #2.'
     else 
        print *, 'Plasma limiter is unknown!'
     end if
  end if

     
end subroutine lcfs


end module diagnostics

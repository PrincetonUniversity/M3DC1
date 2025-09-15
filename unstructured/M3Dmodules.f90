module basic
  use mesh_mod
  use pid_controller
  use spline

  implicit none

  integer, parameter :: ijacobian = 1

  integer, parameter :: version = 45

#if defined(USE3D) || defined(USECOMPLEX)
  integer, parameter :: i3d = 1
#else
  integer, parameter :: i3d = 0
#endif
#ifdef USECOMPLEX
  integer, parameter :: icomplex = 1
#else
  integer, parameter :: icomplex = 0
#endif

  real, parameter :: c_light = 2.9979e10
  real, parameter :: e_c = 4.8032e-10
  real, parameter :: m_p = 1.6726219e-24
  real, parameter :: m_e = 9.1094e-28
  real, parameter :: me_mp = m_e / m_p
  real, parameter :: mp_me = m_p / m_e
  real, parameter :: N_Avo = 6.022140857e23

  logical :: print_help

  logical :: density_source
  logical :: momentum_source
  logical :: heat_source
  logical :: rad_source

  ! normalizations
  real :: b0_norm     ! magnetic field normalization (in Gauss)
  real :: l0_norm     ! length normalization (in centimeters)
  real :: n0_norm     ! density normalization (in cm^-3)
  real :: v0_norm
  real :: t0_norm
  real :: p0_norm
  real :: e0_norm
  real :: j0_norm
  real :: m0_norm

  ! transport coefficients
  real :: amu         ! incompressible viscosity
  real :: amuc        ! compressible viscosity
  real :: amue        ! bootstrap viscosity coefficient
  real :: amupar      ! parallel viscosity coefficient
  real :: amu_wall
  real :: amu_wall_off
  real :: amu_wall_delt
  integer :: iresfunc   ! if 1, use new resistivity function
  integer :: ivisfunc   ! if 1, use new resistivity function
  integer :: ikappafunc ! select electron thermal conductivity function
  integer :: ikapparfunc ! select electron parallel thermal conductivity function
  integer :: idenmfunc
  integer :: ikappar_ni
  real :: etar, eta0  ! iresfunc=0:  resistivity = etar + eta0/T^(3/2)
  real :: eta_fac
  real :: eta_max, eta_min
  integer :: eta_mod
  real :: eta_te_offset  ! offset in Te when calculating eta
  integer :: ikprad_te_offset  ! if 1, eta_te_offset also applied to kprad and ablation
  real :: etaoff, etadelt !iresfunc=1: = etar + .5 eta0 (1+tanh(psi-psilim(1+etaoff*DP)/etadelt*DP))
  !                                                      DP = psilim - psimin
  real :: amuoff, amudelt, amuoff2, amudelt2
  real :: kappaoff, kappadelt
  real :: lambdae     ! multiplier of electron mass term in psi equation
  real :: kappat      ! isotropic temperature conductivity
  real :: kappa0      ! kappa = kappat + kappa0*n/T^(1/2)
  real :: kappah      ! phenomenological model for H-mode
  real :: kappar      ! coefficient of field-aligned temperature diffusion
  real :: kappari_fac  ! ion parallel thermal conductivity is kappari_fac x electron value
  real :: tcrit ! for ikapparfunc=1, parallel TC is kappar/( (tcrit/te)**2.5 + 1)
  real :: k_fac   ! factor by which the toroidal field is multiplied in 1/B^2 that appears in kappa_parallel
  real :: kappax      ! coefficient of B x Grad[T] temperature diffusion
  real :: kappag
  real :: kappaf
  real :: kappai_fac
  real :: kappa_max
  real :: kappar_max, kappar_min
  real :: denm        ! artificial density diffusion used in idenmfunc = 0,1
  real :: denmt       ! temperature dependent density diffusion used in idenmfunc = 1
  real :: denmmin     ! Minimum value of density diffusion
  real :: denmmax     ! Maximum value of density diffusion
  real :: deex        ! scale length of hyperviscosity term
  real :: hyper,hyperi,hyperv,hyperc,hyperp
  real :: gradp_crit
  real :: temin_qd    ! minimum temperature used in equipartition term for ipres=1
  real :: efac        ! eta = efac / T^(3/2)
  real :: nufac       ! nu = nufac * n / T^(3/2)
  real :: krfac       ! kappar = krfac * T^5/2


  ! physical parameters
  integer :: itor     ! 1 = cylindrical coordinates; 0 = cartesian coordinates
  real :: db          ! ion skin depth
  real :: db_fac      ! factor to scale physical value of db
  real :: gam         ! ratio of specific heats
  real :: gravr,gravz ! gravitational acceleration
  real :: vloop0      ! initial loop voltage
  real :: vloop       ! loop voltage
  real :: vloopRZ     ! rate at which boundary TF changes
  real :: vloop_freq  ! frequency of loop voltage
  real :: mass_ratio  ! me/mi (in units of me/mp)
  real :: z_ion       ! Z of main ion species
  real :: ion_mass    ! Effective mass of ions (in proton mass/particle)
  real :: lambda_coulomb ! coulomb logarithm
  real :: thermal_force_coeff

  ! domain parameters
  real :: rzero    ! nominal major radius of the device
  real :: libetap   !  li/2 + beta_poloidal estimate for equili/2brium calculation
  real :: xlim2   !  x-position of second limiter point as a diagnostic
  real :: zlim2   !  z-position of second limiter point as a diagnostic
  integer :: iwall_is_limiter ! 1 = wall acts as limiter
  integer :: nonrect     ! 1 = non-rectangular boundary; 0 = rectangular boundary
  integer :: ifixedb   !  1 = plasma boundary is mesh boundary (for nonrect=1);   0 = free boundary

  ! boundary conditions
  integer :: inonormalflow   ! 1 = no normal flow
  integer :: inoslip_pol     ! 1 = no slip (poloidal flow)
  integer :: inoslip_tor     ! 1 = no slip (toroidal flow)
  integer :: inostress_tor   ! 1 = no stress (toroidal flow)
  integer :: inocurrent_pol  ! 1 = no tangential current
  integer :: inocurrent_tor  ! 1 = no toroidal current
  integer :: inocurrent_norm ! 1 = no toroidal current
  integer :: iconst_p        ! 1 = pressure held constant, 2 = pressure=pedge
  integer :: iconst_t        ! 1 = temperature held constant
  integer :: iconst_n        ! 1 = density held constant
  integer :: iconst_bz       ! 1 = toroidal field held constant
  integer :: iconst_bn       ! 1 = normal magnetic field held constant
  integer :: inograd_p       ! 1 = no normal pressure gradient
  integer :: inograd_t       ! 1 = no normal temperature gradient
  integer :: inograd_n       ! 1 = no normal density gradient
  integer :: com_bc          ! 1 = forces del^2(chi) = 0 on boundary
  integer :: vor_bc          ! 1 = forces del*(phi) = 0 on boundary
  integer :: ifbound         ! bc on f.  0=none, 1=dirichlet, 2=neumann
  real :: amu_edge    ! factor by which to increase viscosity at boundaries

  real :: eta_vac     ! resistivity of vacuum region

  ! density sources
  integer :: ionization     ! 1 = include edge reionization
  real :: ionization_rate   ! rate of ionization
  real :: coolrate    !  cooling rate for iheat_sink=1 and itaylor=27
  real :: ionization_temp   ! temperature above which ionization occurs
  real :: ionization_depth  ! temperature scale of neutral burnout
  integer :: nosig          ! 1 = drop sigma terms from momentum eqn

  integer :: isink               ! number of density sinks
  real :: sink1_x, sink2_x       ! x coordinates of sinks
  real :: sink1_z, sink2_z       ! z coordinates of sinks
  real :: sink1_rate, sink2_rate ! rate of sinks
  real :: sink1_var, sink2_var   ! spatial dispersion of sinks

  integer :: iarc_source         ! density source due to halo current
  real :: arc_source_alpha
  real :: arc_source_eta

  integer :: idenfloor    ! defines sources/sinks to keep density constant in vacuum region
  real :: alphadenfloor  ! proportionality constant for setting density constant

  integer :: igaussian_heat_source
  real :: ghs_x, ghs_z, ghs_phi, ghs_rate, ghs_var, ghs_var_tor

  integer :: iprad    ! Use Qian Teng's prad module
  integer :: prad_z   ! Z of impurity species in prad module
  real    :: prad_fz  ! density of impurity species as fraction of electron density 

  ! general equilibrium parameters
  integer :: irestart ! 1 = reads restart file as initial condition
                      ! 2 = reads restart file to initialize GS solve
                      ! 3 = reads 2D RL=1 restart file o initialize 2D COM=1 run
  integer :: irestart_slice   ! field output timeslice from which to restart
  integer :: version_in  ! Version of restart file
  integer :: itaylor  ! equilibrium
  integer :: idevice  ! for itor=1, itaylor=1, selects tokamak configuration
                      !  0 = generic
                      !  1 = CDX-U
                      !  2 = NSTX
  integer :: iupstream  !  if 1, adds 2nd order toroidal derivative to pressure, density, and fields
                        !  if 2, adds 4th order    "         "       "    "         "      "    "
  real :: magus       ! magnitude of upstream term
  real :: bzero       ! guide field
  real :: bx0         ! initial field in x-direction
  real :: vzero       ! initial toroidal velocity
  real :: phizero     ! initial poloidal velocity
  real :: verzero     ! initial vertical velocity
  real :: v0_cyl    !  initial central axial velocity for cylinder equilibrium
  real :: v1_cyl    !  vz = v0_cyl + v1_cyl*psi**beta          0  < psi < 1
  real :: p0, pi0     ! total, ion pressures
  real :: pscale      ! factor by which to scale equilibrium pressure
  real :: bscale      ! factor by which to scale equilibrium toroidal field
  real :: vscale      ! factor by which to scale equilibrium toroidal rotation
  real :: bpscale     ! factor by which to scale F'
  real :: batemanscale! Bateman scaling of toroidal field, keeping current density fixed
  integer :: iread_bscale
  integer :: iread_pscale
  real :: ln          ! length of equilibrium gradient
  real :: eps         ! size of initial perturbation
  integer :: iwave    ! which wave to initialize in wave prop. equilibrium
  integer :: irmp     ! 1 = read rmp coil/currents from rmp_coil.dat, rmp_current.dat
  integer :: iread_ext_field
  integer :: isample_ext_field
  integer :: isample_ext_field_pol

  real :: scale_ext_field
  integer :: type_ext_field ! 0 = text schaffer field; 1,2 = fieldlines or mgrid file.
  character(len=256) :: file_ext_field ! External field (to be subtracted for ST)
  character(len=256) :: file_total_field ! Stellarator field (plasma+coils) to be read for itaylor=41
  real, dimension(8) :: shift_ext_field
  integer :: maxn     ! maximum frequency in random initial conditions

  ! grad-shafranov options
  integer :: divertors! number of divertors
  integer :: igs      ! number of grad-shafranov iterations
  real    :: eta_gs   ! factor for smoothing nonaxisymmetric psi in GS solve
  integer :: igs_pp_ffp_rescale ! rescale p' and FF' based on p and F
  integer :: igs_extend_p ! extend p past psi=1 using te and ne profiles
  integer :: igs_extend_diamag ! extend diamagnetic rotation past psi=1
  integer :: nv1equ   ! if set to 1, use numvar equilibrium for numvar > 1
  real :: xmag, zmag  ! position of magnetic axis
#ifdef USEST
  real :: xmagp, zmagp  ! physical position of magnetic axis
#endif
  real :: xlim, zlim  ! position of limiter
  real :: xdiv, zdiv  ! position of divertor
  real :: tcuro       ! initial toroidal current
  real :: divcur      ! current in divertor (as fraction of tcuro)
  real :: djdpsi
  real :: p1, p2
  real :: pedge       ! pressure in SOL
  real :: tedge       ! electron temperature in SOL
  real :: tebound      ! boundary condition for electron temperature
  real :: tibound      ! boundary condition for ion temperature
  real :: expn        ! density = pressure**expn
  real :: q0          ! safety factor at magnetic axis
  real :: th_gs       ! relaxation factor
  real :: tol_gs      ! error tolorance for GS solver
  real :: psiscale    ! profile scale-factor (psiscale < 1 throws out edge pts)
  real :: sigma0      ! minor radius of initial current distribution to initialize GS
  real :: elongation  ! elongation of solovev equilibrium

  integer :: idenfunc ! specifies a specific form for equilibrium density
  real :: den_edge
  real :: den0
  real :: dendelt
  real :: denoff

  integer :: irot  !  if irot=1, include toroidal rotation in equilibrium
  integer :: iscale_rot_by_p      ! if 0, don't scale rotation by pressure
  real :: alpha0, alpha1, alpha2, alpha3 !  rotation profile is
!                                   (alpha0 + alpha1*psin + alpha2*psin**2)*pressure

  ! model options
  integer :: linear      ! 1 = linear simulation; 0 = nonlinear simulation
  integer :: eqsubtract  ! 1 = subtract equilibrium in noninear simulations
  integer :: extsubtract ! 1 = subtract external fields from solutions
  integer :: numvar      ! 1 = 2-field; 2 = reduced MHD; 3 = compressible MHD
  integer :: idens       ! evolve density
  integer :: ipres       ! evolve total and electron pressures separately
  integer :: ipressplit  ! separate the pressure (or temperature) solves from the field solves
  integer :: itemp       ! advance pressures for itemp=0, Temperatures for itemp=1
  integer :: iadiabat    ! 1:  takes into account density variation when itemp=1
  integer :: imode       ! specifies which of the 4 modes of treating pressure (and temperature) for ipressplit=1
  integer :: imp_bf      ! include bf implicitly
  integer :: imp_temp    ! include implicit equation for temperature
  integer :: gyro        ! include gyroviscosity
  integer :: jadv        ! 1 = use current density equation, not flux equation
  integer :: isources    ! 1 = include "source" terms in velocity advance
  integer :: istatic     ! 1 = do not advance velocity
                         ! 3 = zero out chi only
  integer :: iestatic    ! 1 = do not advance fields
  integer :: igauge
  integer :: ihypeta     ! 1 = scale hyper-resistivity with eta
                         ! 2 = scale hyper-resistivity with pressure for imp_hyper=2
                         ! >2 hyper-resistivity also scaled by keharmonic(ihypeta)
  real :: bharhypeta    ! bharmonic(ihypeta)
  integer :: ihypkappa   ! 1 = scale hyper-diffusivity with kappa
  integer :: ihypdx      ! scale hyper-resistivity with dx**ihypdx
  integer :: imp_hyper   ! 1 = include hyper-resistivity implicitly in psi equation
                         ! 2 = implicit hyper-resistivity is of the Boozer form
  integer :: ikapscale   ! 1 = scale kappar with kappa
  integer :: inertia     ! 1 = include ion inertial terms (v.grad(v))
  integer :: itwofluid   ! 1 = include two-fluid terms in ohm's law (electron form)
                         ! 2 = ion form of 2F equations
                         ! 3 = parallel electron pressure gradient only
  integer :: no_vdg_T    ! 1 = do not include the V dot Grad(T) terms in temperature equation (for debug) 
  integer :: ibootstrap  ! bootstrap current model
  integer :: irunaway    ! runaway electron model
  integer :: cre         ! runaway speed
  integer :: ra_cyc      ! runaway subcycle
  real :: radiff         ! runaway diffusion
  real :: rjra           ! jra/j0
  integer :: runaway_characteristics           ! use method of characteristics
  integer :: iflip       ! 1 = flip handedness
  integer :: iflip_b     ! 1 = flip equilibrium toroidal field
  integer :: iflip_j     ! 1 = flip equilibrium toroidal current density
  integer :: iflip_v     ! 1 = flip equilibrium toroidal velocity
  integer :: iflip_z     ! 1 = flip equilibrium across z=0 plane
  integer :: ieq_bdotgradt ! 1 = include equilibrium parallel T gradient term
  integer :: icsubtract  ! 1 = subtract fields from poloidal field coils
  integer :: kinetic     ! 1 = use kinetic PIC hot ion pressure tensor
                         ! 2 = CGL form for the pressure tensor (incompressible)
                         ! 3 = CGL form for pressure tensor (full)
#ifdef USEPARTICLES
  integer :: kinetic_fast_ion, kinetic_thermal_ion
  integer :: igyroaverage
  integer :: particle_linear
  integer :: particle_substeps
  integer :: particle_subcycles
  integer :: particle_couple
  integer :: particle_nodelete
  integer :: iconst_f0
  real :: fast_ion_mass, fast_ion_z
  integer :: fast_ion_dist
  real :: fast_ion_max_energy
  integer :: num_par_max
  real, dimension(2) :: num_par_scale
  real, dimension(2) :: kinetic_nrmfac_scale
  integer :: ikinetic_vpar
  real :: kinetic_rhomax
  real :: vpar_reduce
  real :: smooth_par, smooth_dens_parallel
#endif

  integer :: iohmic_heating  ! 1 = include ohmic heating
  integer :: irad_heating  ! 1 = include radiation heat source

  ! numerical parameters
  integer :: ntimemax    ! number of timesteps
  integer :: nskip       ! number of timesteps per matrix recalculation
  integer :: pskip       ! number of timesteps per perconditioner recalculation
  integer :: iskippc     ! number of times preconditioner is reused
  integer :: iconstflux  ! 1 = conserve toroidal flux
  integer :: integrator  ! 0 = Crank-Nicholson, 1 = BDF2
  integer :: isplitstep  ! 1 = do timestep splitting
  integer :: imp_mod
  integer :: iteratephi  ! 1 = iterate field solve
  integer :: icsym  ! symmetry of initial conditions (0) no; (1) even in U: (2) odd in U
  integer :: inumgs  ! if 1, use numerical definition of p and g from profile-p and profile-g files
  integer :: irecalc_eta ! 1 = recalculate transport coeffs after den solve
  integer :: iconst_eta  ! 1 = don't evolve resistivity
  integer :: int_pts_main! pol points in int. quad. for time advance matrices
  integer :: int_pts_aux ! pol points in int. quad. for aux. var, definitions
  integer :: int_pts_diag! pol points in int. quad. for diagnostic calculations
  integer :: int_pts_tor ! tor points in int. quad.
  integer :: isurface    ! include surface terms
  integer :: equilibrate ! 1 = scale trial functions so L2 norms = 1
  integer :: itime_independent ! 1 = exclude d/dt terms
  integer :: iset_pe_floor, iset_pi_floor
  integer :: iset_te_floor, iset_ti_floor
  integer :: iset_ne_floor, iset_ni_floor
  integer :: idiff       ! 1 = solve for difference in solution in B,p from n to n+1
  integer :: idifv       ! 1 = solve for difference in solution in v from n to n+1
  integer :: ksp_max     ! if.gt.0  max number of petsc iterations before time step is repeated
  integer :: max_repeat  ! max number of times time-step is repeated
  integer :: ksp_warn    ! time step is  reduced  if max Petsc iterations > ksp_warn
  integer :: ksp_min     ! time step is increased if max Petsc iterations < ksp_min
  real :: dt, dtold      ! timestep (present and previous)
  real :: dtmin,dtmax,dtkecrit,dtfrac  ! quantities used in variable_timestep option
  real :: ddt            ! change in timestep per timestep
  real :: thimp          ! implicitness parameter (for Crank-Nicholson)
  real :: thimpsm        ! implicitness parameter for smoothers
  real :: regular        ! regularization constant in chi equation
  real :: max_ke         ! max KE before fields are re-scaled when linear==1
  real :: chiiner        ! factor to multiply chi inertial terms
  real :: harned_mikic   ! coefficient of harned-mikic 2f stabilization term
  real :: gamma_gr       ! growth rate based on kinetic energy -- used in variable_timestep
  real :: pe_floor, pi_floor
  real :: te_floor, ti_floor
  real :: ne_floor, ni_floor
  real :: frequency      ! frequency in time-independent calculation
  real :: caramana_fac   ! coefficient for the explicit term in Caramana method (imp_mod=1)

  ! poloidal force parameters
  integer :: ipforce     ! 1 = include poloidal momentum source
  integer :: nforce         !  exponent
  real :: dforce         !  half width
  real :: xforce         !  location of peak (0 to 1(edge))
  real :: aforce         !  magnitude

  integer :: ifixed_temax   !  if nonzero, evaluate temax at xmag0,zmag0

  ! curent drive source
  integer :: icd_source  ! 1 = include current drive in flux equation
  real :: j_0cd          ! amplitude of Gaussian
  real :: r_0cd          ! r coordinate of peak
  real :: z_0cd          ! z coordinate of peak
  real :: w_cd           ! width of Gaussian
  real :: delta_cd       ! shift of Gaussian

  ! xray diagnostic parameters
  integer :: xray_detector_enabled ! 1 = enable xray diagnostic
  real :: xray_r0        ! R coordinate of xray detector
  real :: xray_phi0      ! Phi coordinate of xray detector
  real :: xray_z0        ! Z coordinate of xray detector
  real :: xray_theta     ! angle of xray detector chord (degrees)
  real :: xray_sigma     ! spread of xray detector chord (degrees)

  ! current controller parameters
  real :: tcur           ! target toroidal current (constant in time)
  real :: tcuri          ! initial current for variable target current
  real :: tcurf          ! final current for variable target current
  real :: tcur_t0        ! transition time for variable target current
  real :: tcur_tw        ! transition time width for variable target current
  real :: control_p      ! proportionality constant
  real :: control_i      ! integral control inverse time-scale
  real :: control_d      ! derivative control time-scale
  integer :: control_type  ! type of current control  0:orig  1:standard
  ! density controller parameters
  real :: n_target       ! target particle number
  real :: n_control_p    ! proportionality constant
  real :: n_control_i    ! integral control inverse time-scale
  real :: n_control_d    ! derivative control time-scale
  integer :: n_control_type  ! type of density control  0:orig 1:standard


  ! input/output parameters
  integer :: iprint        ! print extra debugging info
  integer :: itimer        ! print timing info
  integer :: ntimepr       ! number of timesteps per output  
  integer :: ntimers       ! number of timesteps per restart output
  integer :: icalc_scalars ! 1 = calculate scalars
  integer :: ike_only      ! 1 = only calculate kinetic energy
  integer :: ike_harmonics  ! number of Fourier harmonics of ke to be calculated and output
  integer :: ibh_harmonics  ! number of Fourier harmonics of magnetic field perturbation to be calculated and output
  integer :: ifout         ! 1 = output f field
  integer :: irestart_fp   ! -1 = default; 1 = fp present at restart; 0 = absent; 
  integer :: itemp_plot    ! 1 =output vdotgradt, deldotq_perp, deldotq_par,eta_jsq
  integer :: ibdgp         ! option to make partial plots of b dot grad potential
  integer :: iveldif         ! option to make partial plots of V x B - grad(potential)
  integer :: iread_eqdsk   ! 1 = read geqdsk input
                           ! 2 = read geqdsk for psi, but use default profiles
  integer :: iread_dskbal  ! 1 = read dskbal input
  integer :: iread_jsolver
  integer :: iread_omega
  integer :: iread_omega_e
  integer :: iread_omega_ExB
  integer :: iread_ne      
  integer :: iread_te
  integer :: iread_p
  integer :: iread_f
  integer :: iread_j
  integer :: iread_heatsource ! 1 = read heat source profile (in terms of Psi normalized), source is scaled with ghs_rate
  integer :: iread_particlesource ! 1 = read particle source profile (in terms of Psi normalized), source is scaled with pellet_rate
  integer :: iheat_sink   !  add a sink term in p equation (initially for itaylor=27)
  integer :: iread_neo      ! 1 = read velocity profiles from NEO output
  integer :: ineo_subtract_diamag ! 1 = subtract v* from input v profile

  ! adaptation options
  integer :: iadapt     ! 1,2 = adapts mesh after initialization
  real :: adapt_factor
  real :: adapt_hmin
  real :: adapt_hmax  
  real :: adapt_smooth  ! value which controls smoothing of size field
  real :: adapt_psin_vacuum ! value of psin in vacuum region for adaptation
  real :: adapt_psin_wall   ! value of psin in wall region for adaptation
  integer :: iadapt_pack_rationals
  real :: adapt_pack_factor

  integer :: ispradapt
  integer :: isprrefinelevel
  integer :: isprcoarsenlevel
  integer :: isprntime
  real :: isprmaxsize
  real :: isprweight

  real :: beta
  real :: pefac

  ! complex options
  integer :: ntor     ! toroidal mode number
  integer :: mpol     ! poloidal mode number for certain test problems
  complex :: rfac

  !     derived quantities
  real :: bdf,hypv,hypc,hypf,hypi,hypp,   &
       time,                                     &
       gbound

  logical :: use_external_fields = .false.

  ! magnetic diagnostics
  real :: psimin            ! flux value at magnetic axis
  real :: psilim,psilim2    ! flux at the limiter
  real :: psibound          ! flux at the lcfs
  logical :: is_diverted    ! whether plasma is diverted or not
  real :: xnull, znull      ! coordinates of the limiting x-point
  real :: xnull2, znull2    ! coordinates of the limiting x-point
  real :: psinull, psinull2
  integer :: mod_null_rs, mod_null_rs2  ! if 1, modify xnull,znull or xnull2,znull2 at restart
  real :: temax, temax_readin      ! maximum temperature, reading in for ibootstrap=3

  integer :: isolve_with_guess=0 ! (=0; use zero initial guess); (=1; use previous step value as non-zero initial guess)

  ! PID controllers
  type(pid_control), save :: i_control, n_control

  integer :: ntime, ntime0

  ! Deprecated
  real :: zeff_xxx       ! Effective Z of ion fluid

  ! MPI variable(s)
  integer myrank, maxrank
#if defined(_OPENACC) || defined(_OPENMP)
  integer :: num_devices, igpu
#endif

  type(spline1d) :: q_spline

end module basic

module arrays
  use field
  use element

  integer, parameter :: r8 = selected_real_kind(12,100)

  ! arrays defined at all vertices
  ! any change to this list of variables needs to be taken into
  ! account in the arrayresizevec subroutine

  ! Arrays containing physical fields
  type(vector_type), target :: field_vec, field0_vec
  type(vector_type), target :: field_vec_pre

  ! Arrays containing external fields
  type(field_type) :: psi_ext, bz_ext, bf_ext, bfp_ext

  ! Arrays containing auxiliary variables
  type(field_type) :: jphi_field
  type(field_type) :: resistivity_field, kappa_field, kappar_field, denm_field
  type(field_type) :: sigma_field, Fphi_field, Q_field, cd_field
  type(field_type) :: Totrad_field, Linerad_field, Bremrad_field, Ionrad_field, Reckrad_field, Recprad_field
  type(field_type) :: visc_field, visc_c_field, visc_e_field, pforce_field, pmach_field
  type(field_type) :: Jbs_L31_field, Jbs_L32_field, Jbs_L34_field, Jbs_alpha_field, Jbs_fluxavg_iBsq_field &
  , Jbs_fluxavg_G_field, Jbs_dtedpsit_field
  type(field_type) :: Jbs_ftrap_field,Jbs_qR_field,Jbs_invAspectRatio_field 
  type(field_type) :: temporary_field
  
  type(field_type) :: psi_coil_field

  ! the indicies of the named fields within the field vector
  integer, parameter :: u_g = 1
  integer, parameter :: psi_g = 2
  integer, parameter :: vz_g = 3
  integer, parameter :: bz_g = 4
  integer, parameter :: chi_g = 5
  integer, parameter :: pe_g = 6
  integer, parameter :: den_g = 7
  integer, parameter :: p_g = 8
  integer, parameter :: te_g = 9
  integer, parameter :: ti_g = 10
  integer, parameter :: e_g = 11
  integer, parameter :: ne_g = 12
  integer, parameter :: nre_g = 13
  integer, parameter :: num_fields = 13


  type(field_type) :: u_field(0:1), vz_field(0:1), chi_field(0:1)
  type(field_type) :: psi_field(0:1), bz_field(0:1), pe_field(0:1)
  type(field_type) :: den_field(0:1), p_field(0:1), ne_field(0:1)
  type(field_type) :: bf_field(0:1), bfp_field(0:1), e_field(0:1)
  type(field_type) :: te_field(0:1), ti_field(0:1)
  type(field_type) :: u_field_pre, psi_field_pre
  type(field_type) :: nre_field(0:1)  ! runaway electron density
  type(field_type) :: wall_dist
#ifdef USEST
  type(field_type) :: rst, zst ! Stellarator geometry field
#endif
#ifdef USEPARTICLES
  type(field_type) :: rho_field, nf_field, tf_field, pf_field, vfpar0_field
  type(field_type) :: nfi_field, tfi_field, pfi_field, densmooth_field, vparsmooth_field
  type(field_type) :: epar_field, den2_field

  type(field_type) :: p_f_par, p_f_perp  !Kinetic pressure tensor components
  type(field_type) :: p_i_par, p_i_perp  !Kinetic pressure tensor components
  type(field_type) :: den_i_0, den_i_1, den_f_0, den_f_1
  type(field_type) :: v_i_par
  type(field_type) :: v_f_par
#endif

  ! the following pointers point to the locations of the named field within
  ! the respective vector.  set by assign_local_pointers()
  vectype, dimension(dofs_per_node) ::   u1_l,   u0_l
  vectype, dimension(dofs_per_node) ::  vz1_l,  vz0_l
  vectype, dimension(dofs_per_node) :: chi1_l, chi0_l
  vectype, dimension(dofs_per_node) :: psi1_l, psi0_l
  vectype, dimension(dofs_per_node) ::  bz1_l,  bz0_l
  vectype, dimension(dofs_per_node) ::  bfp1_l,  bfp0_l
  vectype, dimension(dofs_per_node) ::  e1_l,  e0_l
  vectype, dimension(dofs_per_node) ::  pe1_l,  pe0_l
  vectype, dimension(dofs_per_node) :: den1_l, den0_l
  vectype, dimension(dofs_per_node) ::   p1_l,   p0_l
  vectype, dimension(dofs_per_node) ::  te1_l,  te0_l
  vectype, dimension(dofs_per_node) ::  ti1_l,  ti0_l
  vectype, dimension(dofs_per_node) ::  qe1_l,  qe0_l
  vectype, dimension(dofs_per_node) ::  qi1_l,  qi0_l
  vectype, dimension(dofs_per_node) :: nre1_l, nre0_l


contains

  !======================================================
  ! get_local_vals
  ! ~~~~~~~~~~~~~~
  !======================================================
  subroutine get_local_vals(inode)

    use basic

    implicit none

    integer, intent(in) :: inode

    call get_node_data(  u_field(0), inode,   u0_l)
    call get_node_data( vz_field(0), inode,  vz0_l)
    call get_node_data(chi_field(0), inode, chi0_l)
    call get_node_data(psi_field(0), inode, psi0_l)
    call get_node_data( bz_field(0), inode,  bz0_l)
    call get_node_data( pe_field(0), inode,  pe0_l)
    call get_node_data(  p_field(0), inode,   p0_l)
    call get_node_data(den_field(0), inode, den0_l)
    call get_node_data( te_field(0), inode,  te0_l)
    call get_node_data( ti_field(0), inode,  ti0_l)
    call get_node_data(nre_field(0), inode, nre0_l)
    call get_node_data(  u_field(1), inode,   u1_l)
    call get_node_data( vz_field(1), inode,  vz1_l)
    call get_node_data(chi_field(1), inode, chi1_l)
    call get_node_data(psi_field(1), inode, psi1_l)
    call get_node_data( bz_field(1), inode,  bz1_l)
    call get_node_data( pe_field(1), inode,  pe1_l)
    call get_node_data(  p_field(1), inode,   p1_l)
    call get_node_data(den_field(1), inode, den1_l)
    call get_node_data( te_field(1), inode,  te1_l)
    call get_node_data( ti_field(1), inode,  ti1_l)
    call get_node_data(nre_field(1), inode, nre1_l)

  end subroutine get_local_vals

  subroutine set_local_vals(inode)
    use basic
    implicit none

    integer, intent(in) :: inode

    call set_node_data(  u_field(0), inode,   u0_l)
    call set_node_data( vz_field(0), inode,  vz0_l)
    call set_node_data(chi_field(0), inode, chi0_l)
    call set_node_data(psi_field(0), inode, psi0_l)
    call set_node_data( bz_field(0), inode,  bz0_l)
    call set_node_data( pe_field(0), inode,  pe0_l)
    call set_node_data(  p_field(0), inode,   p0_l)
    call set_node_data(den_field(0), inode, den0_l)
    call set_node_data( te_field(0), inode,  te0_l)
    call set_node_data( ti_field(0), inode,  ti0_l)
    call set_node_data(nre_field(0), inode, nre0_l)
    call set_node_data(  u_field(1), inode,   u1_l)
    call set_node_data( vz_field(1), inode,  vz1_l)
    call set_node_data(chi_field(1), inode, chi1_l)
    call set_node_data(psi_field(1), inode, psi1_l)
    call set_node_data( bz_field(1), inode,  bz1_l)
    call set_node_data( pe_field(1), inode,  pe1_l)
    call set_node_data(  p_field(1), inode,   p1_l)
    call set_node_data(den_field(1), inode, den1_l)
    call set_node_data( te_field(1), inode,  te1_l)
    call set_node_data( ti_field(1), inode,  ti1_l)
    call set_node_data(nre_field(1), inode, nre1_l)

  end subroutine set_local_vals
end module arrays
  
module sparse
  use matrix_mod

  integer, parameter :: numvar1_numbering = 1
  integer, parameter :: numvar2_numbering = 2
  integer, parameter :: numvar3_numbering = 3
  integer, parameter :: numvar4_numbering = 4
  integer, parameter :: numvar5_numbering = 5
  integer, parameter :: numvar6_numbering = 6
  integer, parameter :: numvar7_numbering = 7
  integer, parameter :: numvar8_numbering = 8

  integer, parameter :: s8_mat_index = 1
  integer, parameter :: s7_mat_index = 2
  integer, parameter :: s4matrix_sm = 3
  integer, parameter :: s5_mat_index = 4
  integer, parameter :: s1_mat_index = 5
  integer, parameter :: s2_mat_index = 6
  integer, parameter :: d1_mat_index = 7
  integer, parameter :: d2_mat_index = 8
  integer, parameter :: d4matrix_sm = 9
  integer, parameter :: d8_mat_index = 10
  integer, parameter :: q1_mat_index = 11
  integer, parameter :: r2_mat_index = 12
  integer, parameter :: r8_mat_index = 13
  integer, parameter :: q2_mat_index = 14
  integer, parameter :: q8_mat_index = 15
  integer, parameter :: gsmatrix_sm = 16
  integer, parameter :: s9_mat_index = 17
  integer, parameter :: d9_mat_index = 18
  integer, parameter :: r9_mat_index = 19
  integer, parameter :: q9_mat_index = 20
  integer, parameter :: r14_mat_index = 21
  integer, parameter :: mass_mat_lhs_index = 22
  integer, parameter :: mass_mat_lhs_dc_index = 23
  integer, parameter :: mass_mat_rhs_dc_index = 24
  integer, parameter :: o1_mat_index = 25
  integer, parameter :: o2_mat_index = 26
  integer, parameter :: gs_mat_rhs_dc_index = 27
  integer, parameter :: lp_mat_rhs_index = 28
  integer, parameter :: lp_mat_rhs_dc_index = 29
  integer, parameter :: bf_mat_lhs_dc_index = 30
  integer, parameter :: q42_mat_index = 31
  integer, parameter :: d5_mat_index = 32
  integer, parameter :: s10_mat_index = 33
  integer, parameter :: d10_mat_index = 34
  integer, parameter :: d7_mat_index = 36
  integer, parameter :: ppmatrix_lhs = 37
  integer, parameter :: br_mat_index = 38
  integer, parameter :: bfp_mat_rhs_index = 39
  integer, parameter :: dp_mat_lhs_index = 40
  integer, parameter :: mass_mat_rhs_index = 41
  integer, parameter :: rwpsi_mat_index = 42
  integer, parameter :: rwbf_mat_index = 43
  integer, parameter :: ecpsi_mat_index = 44
  integer, parameter :: ecbf_mat_index = 45
  integer, parameter :: rw_rhs_mat_index = 46
  integer, parameter :: rw_lhs_mat_index = 47
  integer, parameter :: o9_mat_index = 48
  integer, parameter :: p1_mat_index = 49
  integer, parameter :: gs_mat_rhs_index = 50
  integer, parameter :: r42_mat_index = 51
  integer, parameter :: bf_mat_index = 52
  integer, parameter :: s11_mat_index = 53
  integer, parameter :: s12_mat_index = 54
  integer, parameter :: d11_mat_index = 55
  integer, parameter :: d12_mat_index = 56
  integer, parameter :: vpol_mat_index = 57
  integer, parameter :: o3_mat_index = 58
  integer, parameter :: psi_mat_index = 59
  integer, parameter :: qp42_mat_index = 60
  integer, parameter :: rp42_mat_index = 61
  integer, parameter :: dr_mat_index = 62
  integer, parameter :: lp_mat_index = 63
  integer, parameter :: wall_mat_index = 64
  integer, parameter :: kprad_lhs_index = 65
  integer, parameter :: kprad_rhs_index = 66
  integer, parameter :: s15_mat_index = 67
  integer, parameter :: d15_mat_index = 68
  integer, parameter :: r15_mat_index = 69
  integer, parameter :: q15_mat_index = 70
  integer, parameter :: k15_mat_index = 71
  integer, parameter :: q43_mat_index = 72
  integer, parameter :: r43_mat_index = 73
  integer, parameter :: pot2_mat_lhs_index = 74
  integer, parameter :: st_mat_index = 75
  integer, parameter :: hypv_lhs_index = 76
  integer, parameter :: hypv_rhs_index = 77
  integer, parameter :: num_matrices = 77

  type(matrix_type) :: rwpsi_mat, rwbf_mat, ecpsi_mat, ecbf_mat
  type(matrix_type), save :: rw_rhs_mat, rw_lhs_mat

  logical :: sparse_initialized = .false.

contains
  subroutine delete_matrices
    implicit none

#ifdef USESCOREC
    integer :: i

    if(.not. sparse_initialized) return

    do i=1, num_matrices
#ifdef M3DC1_TRILINOS
       call m3dc1_epetra_delete(i)
#else
       call m3dc1_matrix_delete(i)
#endif
    end do
#endif
  end subroutine delete_matrices
  
end module sparse

!cjomp
module m3dc1_omp
#ifdef _OPENMP
  implicit none

  include 'omp_lib.h'
  integer :: ithread, nthreads
#endif
end module m3dc1_omp

!cj velocity equation profiling
module m3dc1_vel_prof
#if PETSC_VERSION >= 39
#include <petsc/finclude/petscksp.h>
      use petscksp
#else
#include <petsc/finclude/petscsysdef.h>
#endif
      implicit none
      PetscLogStage  stageA,stageS
end module m3dc1_vel_prof

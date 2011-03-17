module p_data
  implicit none
 
  integer, parameter :: ijacobian = 1
end module p_data

module basic
  use p_data
  use mesh_mod
  use pid_controller

  integer, parameter :: version = 4

  integer :: i3d
#ifdef USECOMPLEX
  integer, parameter :: icomplex = 1
#else
  integer, parameter :: icomplex = 0
#endif

  real, parameter :: me_mi = 1./1836.2

  ! normalizations
  real :: b0_norm     ! magnetic field normalization (in Gauss)
  real :: l0_norm     ! length normalization (in centimeters)
  real :: n0_norm     ! density normalization (in cm^-3)

  ! transport coefficients
  real :: amu         ! incompressible viscosity
  real :: amuc        ! compressible viscosity
  real :: amue        ! bootstrap viscosity coefficient
  real :: amupar      ! parallel viscosity coefficient
  integer :: iresfunc   ! if 1, use new resistivity function
  integer :: ivisfunc   ! if 1, use new resistivity function
  integer :: ikappafunc ! if 1, use new resistivity function
  real :: etar, eta0  ! iresfunc=0:  resistivity = etar + eta0/T^(3/2)
  real :: etaoff, etadelt !iresfunc=1: = etar + .5 eta0 (1+tanh(psi-psilim(1+etaoff*DP)/etadelt*DP))
  !                                                      DP = psilim - psimin
  real :: amuoff, amudelt, amuoff2, amudelt2
  real :: kappaoff, kappadelt
  real :: lambdae     ! multiplier of electron mass term in psi equation
  real :: kappat      ! isotropic temperature conductivity
  real :: kappa0      ! kappa = kappat + kappa0*n/T^(1/2)
  real :: kappah      ! phenomenological model for H-mode
  real :: kappar      ! coefficient of field-aligned temperature diffusion
  real :: kappax      ! coefficient of B x Grad[T] temperature diffusion
  real :: denm        ! artificial density diffusion
  real :: deex        ! scale length of hyperviscosity term
  real :: hyper,hyperi,hyperv,hyperc,hyperp


  ! physical parameters
  integer :: itor     ! 1 = cylindrical coordinates; 0 = cartesian coordinates
  real :: db          ! ion skin depth
  real :: gam         ! ratio of specific heats
  real :: gravr,gravz ! gravitational acceleration
  real :: vloop       ! loop voltage
  real :: mass_ratio  ! me/mi (in units of me/mp)

  ! domain parameters
  real :: rzero    ! nominal major radius of the device
  real :: libetap   !  li/2 + beta_poloidal estimate for equili/2brium calculation
  real :: xlim2   !  x-position of second limiter point as a diagnostic
  real :: zlim2   !  z-position of second limiter point as a diagnostic
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
  integer :: iconst_p        ! 1 = pressure held constant
  integer :: iconst_t        ! 1 = pressure held constant
  integer :: iconst_n        ! 1 = density held constant
  integer :: iconst_bz       ! 1 = toroidal field held constant
  integer :: inograd_p       ! 1 = no normal pressure gradient
  integer :: inograd_n       ! 1 = no normal density gradient
  integer :: com_bc          ! 1 = forces del^2(chi) = 0 on boundary
  integer :: vor_bc          ! 1 = forces del*(phi) = 0 on boundary
  integer :: ifbound         ! bc on f.  0=none, 1=dirichlet, 2=neumann
  real :: amu_edge    ! factor by which to increase viscosity at boundaries

  real :: eta_wall    ! resistivity of boundary
  real :: delta_wall  ! thickness of boundary

  ! density sources
  integer :: ipellet  ! 1 = include pellet injection density source
  real :: pellet_x    ! x coordinate of pellet injection
  real :: pellet_z    ! z coordinate of pellet injection
  real :: pellet_rate ! amplitude of pellet density source
  real :: pellet_var  ! spatial dispersion of density source 
  integer :: ionization     ! 1 = include edge reionization
  real :: ionization_rate   ! rate of ionization
  real :: ionization_temp   ! temperature above which ionization occurs
  real :: ionization_depth  ! temperature scale of neutral burnout
  integer :: nosig          ! 1 = drop sigma terms from momentum eqn

  integer :: isink               ! number of density sinks
  real :: sink1_x, sink2_x       ! x coordinates of sinks
  real :: sink1_z, sink2_z       ! z coordinates of sinks
  real :: sink1_rate, sink2_rate ! rate of sinks
  real :: sink1_var, sink2_var   ! spatial dispersion of sinks


  ! general equilibrium parameters
  integer :: irestart ! 1 = reads restart file as initial condition
                      ! 2 = reads restart file to initialize GS solve
  integer :: itaylor  ! equilibrium
  integer :: idevice  ! for itor=1, itaylor=1, selects tokamak configuration
                      !  0 = generic
                      !  1 = CDX-U
                      !  2 = NSTX
  integer :: iupstream  !  if 1, adds diffusion term to pressure like upstream differencing
  real :: bzero       ! guide field
  real :: bx0         ! initial field in x-direction
  real :: vzero       ! initial toroidal velocity
  real :: phizero     ! initial poloidal velocity
  real :: p0, pi0     ! total, ion pressures
  real :: pscale      ! factor by which to scale equilibrium pressure
  real :: bscale      ! factor by which to scale equilibrium toroidal field
  real :: ln          ! length of equilibrium gradient
  real :: eps         ! size of initial perturbation
  integer :: iwave    ! which wave to initialize in wave prop. equilibrium
  integer :: irmp     ! 1 = read rmp coil/currents from rmp_coil.dat, rmp_current.dat
  integer :: maxn     ! maximum frequency in random initial conditions


  ! grad-shafranov options
  integer :: divertors! number of divertors
  integer :: igs      ! number of grad-shafranov iterations
  integer :: nv1equ   ! if set to 1, use numvar equilibrium for numvar > 1
  integer :: igs_method  ! 1 = use node-based method (fastest, least accurate)
                         ! 2 = use element-based method, and calculate p from
                         !     input p profile (closest fit to input equil.)
                         ! 3 = use element-based method, and calculate p from
                         !     input p' profile (best gs solution)
  real :: xmag, zmag  ! position of magnetic axis
  real :: xlim, zlim  ! position of limiter
  real :: xdiv, zdiv  ! position of divertor
  real :: tcuro       ! initial toroidal current
  real :: divcur      ! current in divertor (as fraction of tcuro)
  real :: djdpsi
  real :: p1, p2, pedge
  real :: expn        ! density = pressure**expn
  real :: q0          ! safety factor at magnetic axis
  real :: th_gs       ! relaxation factor
  real :: tol_gs      ! error tolorance for GS solver
  real :: psiscale    ! profile scale-factor (psiscale < 1 throws out edge pts)

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
  integer :: numvar      ! 1 = 2-field; 2 = reduced MHD; 3 = compressible MHD
  integer :: idens       ! evolve density
  integer :: ipres       ! evolve total and electron pressures separately
  integer :: imp_bf      ! include bf implicitly
  integer :: gyro        ! include gyroviscosity
  integer :: jadv        ! 1 = use current density equation, not flux equation
  integer :: isources    ! 1 = include "source" terms in velocity advance
  integer :: istatic     ! 1 = do not advance velocity
  integer :: iestatic    ! 1 = do not advance fields
  integer :: igauge
  integer :: ivform      ! 0: V = v grad(phi).  1: V = R^2 v grad(phi)
  integer :: ibform      ! 0: multiply bz equation by r^2
  integer :: ihypeta     ! 1 = scale hyper-resistivity with eta
  integer :: ihypamu     ! 1 = scale hyper-viscosity with visc
  integer :: ihypkappa   ! 1 = scale hyper-diffusivity with kappa
  integer :: ihypdx      ! scale hyper-resistivity with dx**ihypdx
  integer :: ikapscale   ! 1 = scale kappar with kappa
  integer :: inertia     ! 1 = include ion inertial terms (v.grad(v))
  integer :: itwofluid   ! 1 = include two-fluid terms in ohm's law
  integer :: ibootstrap  ! bootstrap current model
  integer :: iflip       ! 1 = flip handedness
  integer :: iflip_b     ! 1 = flip equilibrium toroidal field
  integer :: iflip_j     ! 1 = flip equilibrium toroidal current density
  integer :: iflip_v     ! 1 = flip equilibrium toroidal velocity
  integer :: iflip_z     ! 1 = flip equilibrium across z=0 plane

  ! numerical parameters
  integer :: ntimemax    ! number of timesteps
  integer :: nskip       ! number of timesteps per matrix recalculation
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
  real :: dt             ! timestep
  real :: ddt            ! change in timestep per timestep
  real :: thimp          ! implicitness parameter (for Crank-Nicholson)
  real :: thimpsm        ! implicitness parameter for smoothers
  real :: regular        ! regularization constant in chi equation
  real :: max_ke         ! max KE before fields are re-scaled when linear==1
  real :: chiiner        ! factor to multiply chi inertial terms
  real :: harned_mikic   ! coefficient of harned-mikic 2f stabilization term


  ! current controller parameters
  real :: tcur           ! target toroidal current
  real :: control_p      ! proportionality constant
  real :: control_i      ! integral control inverse time-scale
  real :: control_d      ! derivative control time-scale
  ! density controller parameters
  real :: n_target       ! target toroidal current
  real :: n_control_p    ! proportionality constant
  real :: n_control_i    ! integral control inverse time-scale
  real :: n_control_d    ! derivative control time-scale


  ! input/output parameters
  integer :: iprint        ! print extra debugging info
  integer :: itimer        ! print timing info
  integer :: ntimepr       ! number of timesteps per output  
  integer :: iglobalout    ! 1 = write global restart files
  integer :: iglobalin     ! 1 = read global restart files
  integer :: icalc_scalars ! 1 = calculate scalars
  integer :: ike_only      ! 1 = only calculate kinetic energy
  integer :: ifout         ! 1 = output f field
  integer :: iread_eqdsk   ! 1 = read geqdsk input
                           ! 2 = read geqdsk for psi, but use default profiles
  integer :: iread_dskbal  ! 1 = read dskbal input
  integer :: iread_jsolver
  integer :: iread_omega   ! 2 = read omega from cerfit file
  integer :: iread_ne      
  integer :: iread_te
  integer :: iwrite_restart ! 0 = don't write restart files

  ! adaptation options
  integer :: iadapt     ! 1,2 = adapts mesh after initialization
  real :: adapt_factor
  real :: adapt_hmin
  real :: adapt_hmax  

  real :: beta
  real :: pefac

  ! complex options
  integer :: ntor     ! toroidal mode number
  integer :: mpol     ! poloidal mode number for certain test problems
  complex :: rfac

!
!.....input quantities---defined in subroutine input or in namelist
!
  namelist / inputnl/                                          &
       adapt_factor, adapt_hmax, adapt_hmin, &
       alpha0, alpha1, alpha2, alpha3,  &
       amu, amuc, amudelt, amudelt2, amue, amu_edge, amuoff, amuoff2, amupar, &
       b0_norm, beta, bscale, bx0, bzero, &
       chiiner, com_bc, &
       control_d, control_i, control_p, &
       db, ddt, &
       deex, delta_wall,                                   &
       den_edge, &
       den0, dendelt, denm, denoff, &
       divcur, divertors, &
       djdpsi, &
       dt, &
       eps, eqsubtract, equilibrate, &
       eta_wall, eta0, etadelt, etaoff, etar, &
       expn, &
       gam, gravr, gravz, gyro, &
       harned_mikic,                                   &
       hyper, hyperc, hyperi, hyperp, hyperv, &
       iadapt, ibform, ibootstrap, icalc_scalars, &
       iconst_bz, iconst_eta, iconst_n, iconst_p, iconst_t, iconstflux, &
       icsym, icurv, &
       idenfunc, idens, &
       idevice, iestatic, ifbound, ifixedb, &
       iflip_b, iflip_j, iflip_v, iflip_z, iflip, &
       ifout, igauge, iupstream,                                         &
       iglobalin, iglobalout, &
       igs_method, igs, &
       ihypamu, ihypdx, ihypeta, ihypkappa, &
       ikappafunc, ikapscale, &
       ike_only, imp_bf, imp_mod, inertia, &
       inocurrent_norm, inocurrent_pol, inocurrent_tor, &
       inograd_n, inograd_p, &
       inonormalflow, inoslip_pol, inoslip_tor, inostress_tor, &
       int_pts_aux, int_pts_diag, int_pts_main, int_pts_tor,   &
       integrator, inumgs, &
       ionization_depth, ionization_rate, ionization_temp, ionization, &
       ipellet, iper, ipres, iprint, &
       iread_dskbal, iread_eqdsk, iread_jsolver, &
       iread_ne, iread_omega, iread_te,                                     &
       irecalc_eta, iresfunc, irestart, irmp, irot, iscale_rot_by_p, &
       isink, isources, isplitstep, istatic, isurface, itaylor, &
       iteratephi, itime_independent, itimer, itor, itwofluid,     &
       ivform, ivisfunc, iwave, iwrite_restart, &
       jadv, jper, &
       kappa0, kappadelt, kappah, kappaoff, kappar, kappat, kappax, &
       l0_norm, lambdae, libetap, linear, ln, &
       mass_ratio, max_ke, maxn, mpol, &
       n_control_d, n_control_i, n_control_p, n_target, &
       n0_norm, nonrect, nosig, nplanes, &
       nskip, ntimemax, ntimepr, &
       ntor, numvar, nv1equ, &
       p0, p1, p2, pedge, &
       pellet_rate, pellet_var, pellet_x, pellet_z, &
       phizero, pi0, pscale, psiscale, &
       q0, &
       regular, rzero, &
       sink1_rate, sink1_var, sink1_x, sink1_z, &
       sink2_rate, sink2_var, sink2_x, sink2_z, &
       tcur, tcuro, &
       th_gs, thimp, thimpsm, &
       tiltangled, tol_gs, &
       vloop, vor_bc, vzero, &
       xdiv, xlim, xlim2, xmag, xnull, xzero, &
       zdiv, zlim, zlim2, zmag, znull, zzero

  !     derived quantities
  real :: dbf,bdf,hypv,hypc,hypf,hypi,hypp,   &
       time,                                     &
       gbound

  ! magnetic diagnostics
  real :: psimin            ! flux value at magnetic axis
  real :: psilim,psilim2    ! flux at the limiter
  real :: psibound          ! flux at the lcfs
  logical :: is_diverted    ! whether plasma is diverted or not
  real :: xnull, znull      ! coordinates of the limiting x-point

  ! PID controllers
  type(pid_control), save :: i_control, n_control

  integer :: ntime, ntime0

  ! MPI variable(s)
  integer myrank, maxrank

end module basic

module arrays
  use p_data
  use field
  use element

  integer, parameter :: r8 = selected_real_kind(12,100)

  ! arrays defined at all vertices
  ! any change to this list of variables needs to be taken into
  ! account in the arrayresizevec subroutine

  ! Arrays containing physical fields
  type(vector_type), target :: field_vec, field0_vec

  ! Arrays containing auxiliary variables
  type(field_type) :: jphi_field, vor_field, com_field
  type(field_type) :: resistivity_field, kappa_field, sigma_field
  type(field_type) :: visc_field, visc_c_field, visc_e_field
  type(field_type) :: tempvar_field

  type(field_type) :: temporary_field
  
  type(vector_type) :: external_field
  type(field_type) :: external_psi_field, external_bf_field, external_bz_field

  ! the indicies of the named fields within the field vector
  integer, parameter :: u_g = 1
  integer, parameter :: psi_g = 2
  integer, parameter :: vz_g = 3
  integer, parameter :: bz_g = 4
  integer, parameter :: chi_g = 5
  integer, parameter :: pe_g = 6
  integer, parameter :: den_g = 7
  integer, parameter :: p_g = 8
  integer, parameter :: num_fields = 8

  type(field_type) :: u_field(0:1), vz_field(0:1), chi_field(0:1)
  type(field_type) :: psi_field(0:1), bz_field(0:1), pe_field(0:1)
  type(field_type) :: den_field(0:1), p_field(0:1)
  type(field_type) :: bf_field(0:1)

  ! the following pointers point to the locations of the named field within
  ! the respective vector.  set by assign_local_pointers()
  vectype, dimension(dofs_per_node) ::   u1_l,   u0_l
  vectype, dimension(dofs_per_node) ::  vz1_l,  vz0_l
  vectype, dimension(dofs_per_node) :: chi1_l, chi0_l
  vectype, dimension(dofs_per_node) :: psi1_l, psi0_l
  vectype, dimension(dofs_per_node) ::  bz1_l,  bz0_l
  vectype, dimension(dofs_per_node) ::  pe1_l,  pe0_l
  vectype, dimension(dofs_per_node) :: den1_l, den0_l
  vectype, dimension(dofs_per_node) ::   p1_l,   p0_l


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
    call get_node_data(  u_field(1), inode,   u1_l)
    call get_node_data( vz_field(1), inode,  vz1_l)
    call get_node_data(chi_field(1), inode, chi1_l)
    call get_node_data(psi_field(1), inode, psi1_l)
    call get_node_data( bz_field(1), inode,  bz1_l)
    call get_node_data( pe_field(1), inode,  pe1_l)
    call get_node_data(  p_field(1), inode,   p1_l)
    call get_node_data(den_field(1), inode, den1_l)

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
    call set_node_data(  u_field(1), inode,   u1_l)
    call set_node_data( vz_field(1), inode,  vz1_l)
    call set_node_data(chi_field(1), inode, chi1_l)
    call set_node_data(psi_field(1), inode, psi1_l)
    call set_node_data( bz_field(1), inode,  bz1_l)
    call set_node_data( pe_field(1), inode,  pe1_l)
    call set_node_data(  p_field(1), inode,   p1_l)
    call set_node_data(den_field(1), inode, den1_l)

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
  integer, parameter :: bf_mat_rhs_index = 39
  integer, parameter :: dp_mat_lhs_index = 40
  integer, parameter :: mass_mat_rhs_index = 41
  integer, parameter :: rwpsi_mat_index = 42
  integer, parameter :: rwbf_mat_index = 43
  integer, parameter :: ecpsi_mat_index = 44
  integer, parameter :: ecbf_mat_index = 45
  integer, parameter :: rw_rhs_mat_index = 46
  integer, parameter :: rw_lhs_mat_index = 47
  integer, parameter :: o9_mat_index = 48
  integer, parameter :: num_matrices = 48

  type(matrix_type), target :: s1_mat, d1_mat, q1_mat, r14_mat, o1_mat
  type(matrix_type), target :: q42_mat
  type(matrix_type), target :: s2_mat, d2_mat, r2_mat, q2_mat, o2_mat
  type(matrix_type), target :: s8_mat, d8_mat, r8_mat, q8_mat
  type(matrix_type), target :: s9_mat, d9_mat, r9_mat, q9_mat, o9_mat
  type(matrix_type) :: rwpsi_mat, rwbf_mat, ecpsi_mat, ecbf_mat
  type(matrix_type), save :: rw_rhs_mat, rw_lhs_mat

contains
  subroutine delete_matrices
    implicit none

#ifdef USESCOREC
    integer :: i

    do i=1, num_matrices
       call deletematrix(i)
    end do
#endif
  end subroutine delete_matrices
  
end module sparse

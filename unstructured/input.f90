subroutine add_var_double(name, var, default, desc, grp)
  implicit none

  character(len=*), intent(in) :: name, desc
  real :: var
  real, intent(in) :: default
  integer, intent(in) :: grp

  call add_variable_double(name//char(0), var, default, desc//char(0), grp)
end subroutine add_var_double

subroutine add_var_double_array(name, var, size, default, desc, grp)
  implicit none

  character(len=*), intent(in) :: name, desc
  integer, intent(in) :: size
  real, dimension(*) :: var
  real, intent(in) :: default
  integer, intent(in) :: grp

  call add_variable_double_array(name//char(0), var, size, default, &
       desc//char(0), grp)
end subroutine add_var_double_array

subroutine add_var_int(name, var, default, desc, grp)
  implicit none

  character(len=*), intent(in) :: name, desc
  integer :: var
  integer, intent(in) :: default
  integer, intent(in) :: grp

  call add_variable_int(name//char(0), var, default, desc//char(0), grp)
end subroutine add_var_int

subroutine add_var_string(name, var, len, default, desc, grp)
  implicit none

  character(len=*), intent(in) :: name, desc
  character(len=*) :: var
  character(len=*), intent(in) :: default
  integer, intent(in) :: len
  integer, intent(in) :: grp

  call add_variable_string(name//char(0), var, len, default//char(0), &
       desc//char(0), grp, ' ')
end subroutine add_var_string

subroutine add_group(name, handle)
  implicit none
  character(len=*), intent(in) :: name
  integer, intent(out) :: handle

  call create_group(name//char(0), handle)
end subroutine add_group


!=========================
! input
! ~~~~~
! reads input namelist
!=========================
subroutine input
  use basic

  implicit none

  integer :: ierr

#include "mpif.h"

  call set_defaults

  ! Read input file
  ! ~~~~~~~~~~~~~~~
  if(.not.print_help) then
     call read_namelist("C1input"//char(0), ierr)
     if(ierr.ne.0) call safestop(3)
  end if

  if(myrank.eq.0 .and. iprint.ge.1) print *, " validating input"
  call validate_input

  if(print_help) then
     if(myrank.eq.0) call print_variables(3)
     call safestop(0)
  else
     if(myrank.eq.0) call print_variables(0)
  end if
end subroutine input


!==========================
! set_defaults
! ~~~~~~~~~~~~
! sets defa
!==========================
subroutine set_defaults
  use basic
  use m3dc1_output
  use neutral_beam
  use element
  use pellet
  use mesh_mod
  use gradshafranov

  implicit none

  integer :: model_grp
  integer :: eq_grp 
  integer :: transp_grp
  integer :: bc_grp
  integer :: norm_grp
  integer :: gs_grp
  integer :: num_grp
  integer :: hyper_grp
  integer :: time_grp
  integer :: mesh_grp
  integer :: adapt_grp
  integer :: input_grp
  integer :: output_grp
  integer :: diagnostic_grp
  integer :: source_grp
  integer :: misc_grp
  integer :: deprec_grp


  call add_group("Model Options", model_grp)
  call add_group("Equilibrium", eq_grp)
  call add_group("Grad-Shafranov Solver", gs_grp)
  call add_group("Transport Coefficients", transp_grp)
  call add_group("Hyper Diffusivity", hyper_grp)
  call add_group("Normalizations", norm_grp)
  call add_group("Boundary Conditions", bc_grp)
  call add_group("Time Step", time_grp)
  call add_group("Mesh", mesh_grp)
  call add_group("Mesh Adaptation", adapt_grp)
  call add_group("Numerical Options", num_grp)
  call add_group("Input", input_grp)
  call add_group("Output", output_grp)
  call add_group("Diagnostics", diagnostic_grp)
  call add_group("Sources/Sinks", source_grp)
  call add_group("Miscellaneous", misc_grp)
  call add_group("Deprecated", deprec_grp)


  ! Normalizations
  call add_var_double("b0_norm", b0_norm, 1.e4, &
       "Normalization magnetic field (in G)", norm_grp)
  call add_var_double("n0_norm", n0_norm, 1.e14, &
       "Normalization density (in e-/cm3)", norm_grp)
  call add_var_double("l0_norm", l0_norm, 100.,  &
       "Normalization length (in cm)", norm_grp)
 
  ! Input
  call add_var_int("iread_eqdsk", iread_eqdsk, 0, "", input_grp)
  call add_var_int("iread_dskbal", iread_dskbal, 0, "", input_grp)
  call add_var_int("iread_jsolver", iread_jsolver, 0, "", input_grp)
  call add_var_int("iread_omega", iread_omega, 0, "", input_grp)
  call add_var_int("iread_omega_e", iread_omega_e, 0, &
       "Read electron rotation (same options as iread_omega)", input_grp)
  call add_var_int("iread_omega_ExB", iread_omega_ExB, 0, &
       "Read ExB rotation (same options as iread_omega)", input_grp)
  call add_var_int("iread_ne", iread_ne, 0, "", input_grp)
  call add_var_int("iread_te", iread_te, 0, "", input_grp)
  call add_var_int("iread_p", iread_p, 0, "", input_grp)
  call add_var_int("iread_neo", iread_neo, 0, &
       "Read velocity data from NEO output", input_grp)
  call add_var_int("ineo_subtract_diamag", ineo_subtract_diamag, 0, &
       "Subtract diamag. term from input vel. when reading NEO vel.", &
       input_grp)

  ! Transport parameters
  call add_var_int("ivisfunc", ivisfunc, 0, "", transp_grp)
  call add_var_double("amuoff", amuoff, 0., "", transp_grp)
  call add_var_double("amudelt", amudelt, 0., "", transp_grp)
  call add_var_double("amuoff2", amuoff2, 0., "", transp_grp)
  call add_var_double("amudelt2", amudelt2, 0., "", transp_grp)
  call add_var_double("amu", amu, 0., &
       "Isotropic viscosity", transp_grp)
  call add_var_double("amuc", amuc, 0., &
       "Compressional viscosity", transp_grp)
  call add_var_double("amue", amue, 0., "", transp_grp)
  call add_var_double("amupar", amupar, 0., &
       "Parallel viscosity", transp_grp)
  call add_var_double("amu_edge", amu_edge, 0., "", transp_grp)

  call add_var_int("iresfunc", iresfunc, 0, "", transp_grp)
  call add_var_double("etaoff", etaoff, 0., "", transp_grp)
  call add_var_double("etadelt", etadelt, 0., "", transp_grp)
  call add_var_double("etar", etar, 0., &
       "Isotropic resistivity", transp_grp)
  call add_var_double("eta0", eta0, 0., "", transp_grp)
  call add_var_double("eta_fac", eta_fac, 1., &
       "Uniform resistivity multiplier", transp_grp)

  call add_var_int("ikappafunc", ikappafunc, 0, "", transp_grp)
  call add_var_int("ikapscale", ikapscale, 0, "", transp_grp)
  call add_var_int("ikappar_ni", ikappar_ni, 0, &
       "Include 1/n terms in parallel heat flux", transp_grp)
  call add_var_double("kappaoff", kappaoff, 0., "", transp_grp)
  call add_var_double("kappadelt", kappadelt, 0., "", transp_grp)
  call add_var_double("kappat", kappat, 0., &
       "Isotropic thermal conductivity", transp_grp)
  call add_var_double("kappa0", kappa0, 0., "", transp_grp)
  call add_var_double("kappar", kappar, 0., &
       "Parallel thermal conductivity", transp_grp)
  call add_var_double("kappax", kappax, 0., "", transp_grp)
  call add_var_double("kappah", kappah, 0., "", transp_grp)

  call add_var_double("denm", denm, 0., &
       "Density hyperdiffusion coefficient", transp_grp)
  
  call add_var_double("gam", gam, 5./3., &
       "Ratio of specific heats", misc_grp)
  call add_var_double("db", db, 0., &
       "Collisionless ion skin depth", misc_grp)
  call add_var_double("mass_ratio", mass_ratio, 0., "", misc_grp)
  call add_var_double("lambdae", lambdae, 0., "", misc_grp)
  call add_var_double("zeff", zeff, 1., "Z effective", misc_grp)
  call add_var_double("ion_mass", ion_mass, 1., &
       "Ion mass (in units of m_p)", misc_grp)
  call add_var_double("lambda_coulomb", lambda_coulomb, 17., &
       "Coulomb logarithm", misc_grp)


  ! Model options
  call add_var_int("numvar", numvar, 3, &
       "1: 2-Field;  2: 4-Field;  3: 6-Field", model_grp)
  call add_var_int("linear", linear, 0, &
       "1: Use linearized equations", model_grp)
  call add_var_int("eqsubtract", eqsubtract, 0, &
       "1: Subtract equilibrium fields", model_grp)
  call add_var_int("extsubtract", extsubtract, 0, &
       "1: Subtract fields from non-axisymmetric coils", model_grp)
  call add_var_int("icsubtract", icsubtract, 0, &
       "1: Subtract fields from poloidal field coils", model_grp)
  call add_var_int("idens", idens, 0, &
       "1: Include density equation", model_grp)
  call add_var_int("ipres", ipres, 0, &
       "1: Include total pressure equation", model_grp)
  call add_var_int("ipressplit", ipressplit, 0, &
       "1: Separate pressure solves from field solves", model_grp)
  call add_var_int("itemp", itemp, 0, &
       "1: Advance Temperatures rather than Pressures", model_grp)
  call add_var_int("gyro", gyro, 0, &
       "1: Include Braginskii gyroviscosity", model_grp)
  call add_var_int("igauge", igauge, 0, "", model_grp)
  call add_var_int("inertia", inertia, 1, &
       "1: Include V.Grad(V) terms", model_grp)
  call add_var_int("itwofluid", itwofluid, 1, &
       "1: -electron 2F,  2: ion 2F", model_grp)
  call add_var_int("ibootstrap", ibootstrap, 0, "", model_grp)
  call add_var_int("imp_bf", imp_bf, 0, &
       "1: Include implicit equation for f", model_grp)
  call add_var_int("nosig", nosig, 0, "", model_grp)
  call add_var_int("itor", itor, 0, &
       "1: Use toroidal geometry", model_grp)

  call add_var_double("gravr", gravr, 0., "", model_grp)
  call add_var_double("gravz", gravz, 0., "", model_grp)
  call add_var_int("istatic", istatic, 0, &
       "1: Do not advance velocity fields", model_grp)
  call add_var_int("iestatic", iestatic, 0, &
       "1: Do not advance magnetic fields", model_grp)
  call add_var_double("chiiner", chiiner, 1., "", model_grp)
  call add_var_int("ieq_bdotgradt", ieq_bdotgradt, 1, "", model_grp)
  call add_var_int("iwall_is_limiter", iwall_is_limiter, 1, &
       "1 = Wall acts as limiter", model_grp)
  call add_var_int("no_vdg_T", no_vdg_T,0, &
       "1: do not include V dot grad T in Temp equation (debug)",model_grp)
    
  ! Time step options
  call add_var_int("ntimemax", ntimemax, 20, &
       "Total number of timesteps", time_grp)
  call add_var_int("integrator", integrator, 0, "", time_grp)
  call add_var_int("isplitstep", isplitstep, 1, &
       "0: Unsplit time step;  1: Split time step", time_grp)
  call add_var_int("iteratephi", iteratephi, 0, "", time_grp)
  call add_var_int("imp_mod", imp_mod, 1, &
       "Type of split step.  0: Standard;  1: Caramana", time_grp)
  call add_var_int("idiff", idiff, 0, "only solve for difference in B,p", time_grp)
  call add_var_int("idifv", idifv, 0, "only solve for difference in v", time_grp)
  call add_var_int("irecalc_eta", irecalc_eta, 0, "", time_grp)
  call add_var_int("iconst_eta", iconst_eta, 0, "", time_grp)
  call add_var_int("itime_independent", itime_independent, 0, "", time_grp)
  call add_var_double("thimp", thimp, 0.5, &
       "Implicitness of timestep (.5<thimp<1)", time_grp)
  call add_var_double("thimpsm", thimpsm, 1., "", time_grp)
  call add_var_double("harned_mikic", harned_mikic, 0., "", time_grp)
  call add_var_int("isources", isources, 0, "", time_grp)
  call add_var_int("nskip", nskip, 1, "", time_grp)
  call add_var_int("iskippc", iskippc, 1, "", time_grp)
  call add_var_double("dt", dt, 0.1, &
       "Size of time step", time_grp)
  call add_var_double("ddt", ddt, 0., "", time_grp)

  ! variable_timestep parameters

  call add_var_double("dtmin",dtmin,4.0,"minimum time step",time_grp)
  call add_var_double("dtmax",dtmax,40.,"maximum time step",time_grp)
  call add_var_double("dtkecrit",dtkecrit,0.0,"ekin limit on timestep",time_grp)
  call add_var_double("dtfrac",dtfrac,0.1,"fractional change of time step",time_grp)
  call add_var_int("max_repeat", max_repeat, 3, &
       "maximum number of times a time step can be attempted", time_grp)
  call add_var_int("ksp_max", ksp_max, 10000, &
       "maximum number of ksp iterations without repeating time step", time_grp)
  call add_var_int("ksp_min", ksp_min, 500, &
       "time step is increased if max ksp iterations is less than this", time_grp)
  call add_var_int("ksp_warn", ksp_warn, 1000, &
       "time step is reduced if max ksp iterations exceeds this", time_grp)

  ! Numerical methods
  call add_var_int("jadv", jadv, 1, &
       "Use Del*(psi) eqn. instead of psi eqn.", num_grp)
  call add_var_int("ivform", ivform, 1, &
       "V = R^J Grad(U)XGrad(phi) + R^K V Grad(phi) + R^L Grad(chi) |&
       &0: J=0, K=0, L=0;  1: J=2, K=2, L=-2", num_grp)

  call add_var_int("int_pts_main", int_pts_main, 25, "", num_grp)
  call add_var_int("int_pts_aux", int_pts_aux, 25, "", num_grp)
  call add_var_int("int_pts_diag", int_pts_diag, 25, "", num_grp)
  call add_var_int("int_pts_tor", int_pts_tor, 5, "", num_grp)
  call add_var_double("max_ke", max_ke, 1., &
       "Value of ke at which linear sims are rescaled|(ignore if 0)", num_grp)
  call add_var_int("equilibrate", equilibrate, 0, "", num_grp)
  call add_var_double("regular", regular, 0., "", num_grp)
  call add_var_int("iset_pe_floor", iset_pe_floor, 0, &
       "1: Do not let pe drop below pe_floor", num_grp)
  call add_var_double("pe_floor", pe_floor, 0., &
       "Minimum allowed value for pe when iset_pe_floor=1", num_grp)
  call add_var_int("iprecompute_metric", iprecompute_metric, 0, &
       "1: precompute full metric tensor", num_grp)

  ! Equilibrium 
  call add_var_int("itaylor", itaylor, 0, "", eq_grp)
  call add_var_int("iupstream", iupstream, 0, "", eq_grp)
  call add_var_int("iflip", iflip, 0, "", eq_grp)
  call add_var_int("iflip_b", iflip_b, 0, &
       "Reverse equilibrium toroidal field", eq_grp)
  call add_var_int("iflip_j", iflip_j, 0, &
       "Reverse equilibrium toroidal current", eq_grp)
  call add_var_int("iflip_v", iflip_v, 0, &
       "Reverse equilibrium toroidal velocity", eq_grp)
  call add_var_int("iflip_z", iflip_z, 0, "", eq_grp)
  call add_var_int("icsym", icsym, 0, "", eq_grp)
  call add_var_double("bzero", bzero, 1., "", eq_grp)
  call add_var_double("bx0", bx0, 0., "", eq_grp)
  call add_var_double("vzero", vzero, 0., "", eq_grp)
  call add_var_double("phizero", phizero, 0., "", eq_grp)
  call add_var_double("verzero", verzero, 0., "", eq_grp)
  call add_var_int("idevice", idevice, 0, "", eq_grp)
  call add_var_int("iwave", iwave, 0, "", eq_grp)
  call add_var_double("eps", eps, 0.01, &
       "Magnitude of initial perturbations*", eq_grp)
  call add_var_int("maxn", maxn, 200, "", eq_grp)
  call add_var_int("irmp", irmp, 0, &
       "1: Apply nonaxisym. fields throughout plasma|&
       &2: Apply nonaxisym. fields only at boundaries", eq_grp)

  call add_var_int("iread_ext_field", iread_ext_field, 0, &
       "1: Read external field", eq_grp)
  call add_var_int("isample_ext_field", isample_ext_field, 1, &
       "Factor to down-sample external field data toroidally", eq_grp)
  call add_var_int("isample_ext_field_pol", isample_ext_field_pol, 1, &
       "Factor to down-sample external field data poloidally", eq_grp)
  call add_var_double("scale_ext_field", scale_ext_field, 1., &
       "Factor to scale external field", eq_grp)
  call add_var_double_array("shift_ext_field", shift_ext_field, 8, 0., &
       "Toroidal shift (in deg) of external fields", eq_grp)
  call add_var_double("beta", beta, 0., "", eq_grp)
  call add_var_double("ln", ln, 0., "", eq_grp)
  call add_var_double("elongation", elongation, 1., "", eq_grp)
  
  ! Grad-Shafranov
  call add_var_int("inumgs", inumgs, 0, "", gs_grp)
  call add_var_int("igs", igs, 80, "", gs_grp)
  call add_var_int("igs_pp_ffp_rescale", igs_pp_ffp_rescale, 0, &
       "Rescale p' and FF' to match p and F", gs_grp)
  call add_var_int("igs_extend_p", igs_extend_p, 0, &
       "Extend p past Psi=1 using ne and Te profiles", gs_grp)
  call add_var_int("igs_start_xpoint_search", igs_start_xpoint_search, 0, &
       "Number of GS its. before searching for xpoint", gs_grp)
  call add_var_int("igs_forcefree_lcfs", igs_forcefree_lcfs, -1, &
       "Ensure that GS solution is force-free at LCFS", gs_grp)
  call add_var_int("nv1equ", nv1equ, 0, "", gs_grp)
  call add_var_double("eta_gs", eta_gs, 1e3, &
       "Factor for smoothing nonaxisymmetries in psi in GS solve", gs_grp)
  call add_var_double("tcuro", tcuro, 1., &
       "Total current in initial current filament", gs_grp)
  call add_var_double("xmag", xmag, 1., &
       "R-coordinate of initial current filament", gs_grp)
  call add_var_double("zmag", zmag, 0., &
       "Z-coordinate of initial current filament", gs_grp)
  call add_var_double("xlim", xlim, 0., &
       "R-coordinate of limiter #1", gs_grp)
  call add_var_double("zlim", zlim, 0., &
       "Z-coordinate of limiter #1", gs_grp)
  call add_var_double("xlim2", xlim2, 0., &
       "R-coordinate of limiter #2", gs_grp)
  call add_var_double("zlim2", zlim2, 0., &
       "Z-coordinate of limiter #2", gs_grp)
  call add_var_double("rzero", rzero, -1., "", gs_grp)
  call add_var_double("libetap", libetap, 1.2, "", gs_grp)
  call add_var_double("p0", p0, 0.01, "", gs_grp)
  call add_var_double("pi0", pi0, 0.005, "", gs_grp)
  call add_var_double("p1", p1, 0., "", gs_grp)
  call add_var_double("p2", p2, 0., "", gs_grp)
  call add_var_double("pedge", pedge, -1., &
       "Pressure outside separatrix (ignore if < 0)", gs_grp)
  call add_var_double("tedge", tedge, -1., &
       "Temperature outside separatrix (ignore if < 0)", gs_grp)
  call add_var_double("expn", expn, 0., &
       "Density profile = p^expn", gs_grp)
  call add_var_double("q0", q0, 1., "", gs_grp)
  call add_var_double("sigma0", sigma0, 0., "", gs_grp)
  call add_var_double("djdpsi", djdpsi, 0., "", gs_grp)
  call add_var_double("th_gs", th_gs, 0.8, &
       "Implicitness of GS Picard iterations", gs_grp)
  call add_var_double("tol_gs", tol_gs, 1.e-8, "", gs_grp)
  call add_var_double("psiscale", psiscale, 1., "", gs_grp)
  call add_var_double("pscale", pscale, 1., &
       "Factor multiplying pressure profile", gs_grp)
  call add_var_double("bscale", bscale, 1., &
       "Factor multiplying toroidal field profile", gs_grp)
  call add_var_double("batemanscale", batemanscale, 1., &
       "Bateman scaling factor for TF (keeping current density fixed)", gs_grp)
  call add_var_double("bpscale", bpscale, 1., &
       "Factor multiplying F' (keeping F0 constant)", gs_grp)
  call add_var_int("iread_bscale", iread_bscale, 0, &
       "1: read profile_bscale for factor to scale F", gs_grp)
  call add_var_int("iread_pscale", iread_pscale, 0, &
       "1: read profile_pscale for factor to scale p and p'", gs_grp)
  call add_var_double("vscale", vscale, 1., &
       "Factor multiplying toroidal rotation profile", gs_grp)
  call add_var_double_array("gs_vertical_feedback", gs_vertical_feedback, &
       maxcoils, 0., &
       "Proportional feedback of each coil to vertical displacements", gs_grp)
  call add_var_double_array("gs_radial_feedback", gs_radial_feedback, &
       maxcoils, 0., &
       "Proportional feedback of each coil to radial displacements", gs_grp)

  call add_var_int("irot", irot, 0, &
       "Include toroidal rotation", gs_grp)
  call add_var_int("iscale_rot_by_p", iscale_rot_by_p, 1, &
       "0: omega^2 = 2.*p0*(alphai * Psi^i)/n0|&
       &1: omega^2 = 2.*(alphai * Psi^i)/n0,&
       &2: omega^2 = 2.*(alphai * Psi^i), alphai = a0 + a1*exp(-((psii-a2)/a3)**2) ", gs_grp)
  call add_var_double("alpha0", alpha0, 0., &
       "Constant term in analytic rotation profile", gs_grp)
  call add_var_double("alpha1", alpha1, 0., &
       "Linear term in analytic rotation profile", gs_grp)
  call add_var_double("alpha2", alpha2, 0., &
       "Quadratic term in analytic rotation profile", gs_grp)
  call add_var_double("alpha3", alpha3, 0., &
       "Cubic term in analytic rotation profile", gs_grp)
  
  call add_var_int("idenfunc", idenfunc, 0, "", gs_grp)
  call add_var_double("den_edge", den_edge, 0., "", gs_grp)
  call add_var_double("den0", den0, 1., "", gs_grp)
  call add_var_double("dendelt", dendelt, 0.1, "", gs_grp)
  call add_var_double("denoff", denoff, 1., "", gs_grp)

  call add_var_int("divertors", divertors, 0, "", gs_grp)
  call add_var_double("xdiv", xdiv, 0., "", gs_grp)
  call add_var_double("zdiv", zdiv, 0., "", gs_grp)
  call add_var_double("divcur", divcur, 0.1, "", gs_grp)

  call add_var_double("xnull", xnull, 0., &
       "Guess for R-coordinate of active x-point", gs_grp)
  call add_var_double("znull", znull, 0., &
       "Guess for Z-coordinate of axtive x-point", gs_grp)
  call add_var_double("xnull2", xnull2, 0., &
       "Guess for R-coordinate of inactive x-point", gs_grp)
  call add_var_double("znull2", znull2, 0., &
       "Guess for Z-coordinate of inaxtive x-point", gs_grp)


  ! Hyper diffusion
  call add_var_double("deex", deex, 1., "", hyper_grp)
  call add_var_double("hyper", hyper, 0., "", hyper_grp)
  call add_var_double("hyperc", hyperc, 0., "", hyper_grp)
  call add_var_double("hyperi", hyperi, 0., "", hyper_grp)
  call add_var_double("hyperp", hyperp, 0., "", hyper_grp)
  call add_var_double("hyperv", hyperv, 0., "", hyper_grp)
  call add_var_int("ihypdx", ihypdx, 0, "", hyper_grp)
  call add_var_int("ihypeta", ihypeta, 1, "", hyper_grp)
  call add_var_int("ihypamu", ihypamu, 1, "", hyper_grp)
  call add_var_int("ihypkappa", ihypkappa, 1, "", hyper_grp)
  call add_var_int("imp_hyper", imp_hyper, 0,     &
        "1: implicit hyper-resistivity in psi equation", hyper_grp)


  ! Boundary conditions
  call add_var_int("isurface", isurface, 1, "", bc_grp)
  call add_var_int("icurv", icurv, 2, "", bc_grp)
  call add_var_int("nonrect", nonrect, 0, "", bc_grp)
  call add_var_int("ifixedb", ifixedb, 0, &
       "1: Force psi=0 on boundary", bc_grp)
  call add_var_int("com_bc", com_bc, 0, "", bc_grp)
  call add_var_int("vor_bc", vor_bc, 0, "", bc_grp)
  call add_var_int("iconst_p", iconst_p, 1, &
       "1: Hold pressure constant on boundary", bc_grp)
  call add_var_int("iconst_n", iconst_n, 1, &
       "1: Hold density constant on boundary", bc_grp)
  call add_var_int("iconst_t", iconst_t, 1, &
       "1: Hold temperature constant on boundary", bc_grp)
  call add_var_int("iconst_bn", iconst_bn, 1, &
       "1: Hold normal field constant on boundary", bc_grp)
  call add_var_int("iconst_bz", iconst_bz, 1, &
       "1: Hold toroidal field constant on boundary", bc_grp)
  call add_var_int("inograd_p", inograd_p, 0, "", bc_grp)
  call add_var_int("inograd_t", inograd_t, 0, "", bc_grp)
  call add_var_int("inograd_n", inograd_n, 0, "", bc_grp)
  call add_var_int("inonormalflow", inonormalflow, 1, &
       "1: No-normal-flow boundary condition", bc_grp)
  call add_var_int("inoslip_pol", inoslip_pol, 1, &
       "1: No-slip boundary condition on pol. velocity", bc_grp)
  call add_var_int("inoslip_tor", inoslip_tor, 1, &
       "1: No-slip boundary condition on tor. velocity", bc_grp)
  call add_var_int("inostress_tor", inostress_tor, 0, "", bc_grp)
  call add_var_int("inocurrent_pol", inocurrent_pol, 0, "", bc_grp)
  call add_var_int("inocurrent_tor", inocurrent_tor, 0, "", bc_grp)
  call add_var_int("inocurrent_norm", inocurrent_norm, 0, "", bc_grp)
  call add_var_int("ifbound", ifbound, 1, "", bc_grp)
  call add_var_int("iconstflux", iconstflux, 0, "", bc_grp)
  call add_var_int("iper", iper, 0, &
       "1: Periodic boundary condition in R direction", bc_grp)
  call add_var_int("jper", jper, 0, &
       "1: Preiodic boundary condition in Z direction", bc_grp)

  
  ! resistive wall
  call add_var_double("eta_wall", eta_wall, 1e-3, &
       "Resistivity of conducting wall region", misc_grp)
  call add_var_double("delta_wall", delta_wall, 1., "", misc_grp)


  ! loop voltage
  call add_var_double("vloop", vloop, 0., "", source_grp)
  call add_var_double("tcur", tcur, 0., "", source_grp)
  call add_var_double("control_p", control_p, 0., "", source_grp)
  call add_var_double("control_i", control_i, 0., "", source_grp)
  call add_var_double("control_d", control_d, 0., "", source_grp)
  call add_var_int("control_type", control_type, 0, "", source_grp)

  
  ! density source
  call add_var_int("ipellet", ipellet, 0, &
       "1 = include a gaussian pellet source", source_grp)
  call add_var_double("pellet_x", pellet_x, 0., &
       "Initial radial position of the pellet", source_grp)
  call add_var_double("pellet_phi", pellet_phi, 0., &
       "Initial toroidal position of the pellet", source_grp)
  call add_var_double("pellet_z", pellet_z, 0., &
       "Initial vertical position of the pellet", source_grp)
  call add_var_double("pellet_rate", pellet_rate, 0., "", source_grp)
  call add_var_double("pellet_var", pellet_var, 1., "", source_grp)
  call add_var_double("pellet_velx", pellet_velx, 0., &
       "Radial velocity of the pellet", source_grp)
  call add_var_double("pellet_velphi", pellet_velphi, 0., &
       "Toroidal velocity of the pellet", source_grp)
  call add_var_double("pellet_velz", pellet_velz, 0., &
       "Vertical velocity of the pellet", source_grp)

  ! beam source
  call add_var_int("ibeam", ibeam, 0, &
       "1: Include neutral beam source", source_grp)
  call add_var_double("beam_x", beam_x, 0., &
       "R-coordinate of beam center", source_grp)
  call add_var_double("beam_z", beam_z, 0., &
       "Z-coordinate of beam center", source_grp)
  call add_var_double("beam_v", beam_v, 1.e4, &
       "Beam voltage (in volts)", source_grp)
  call add_var_double("beam_rate", beam_rate, 0., &
       "Ions/second deposited by beam", source_grp)
  call add_var_double("beam_dr", beam_dr, 0.1, &
       "Dispersion of beam deposition", source_grp)
  call add_var_double("beam_dv", beam_dv, 100., &
       "Dispersion of beam voltage (in volts)", source_grp)
  call add_var_double("beam_fracpar", beam_fracpar, 1.0, &
       "Cosine of beam angle relative to parallel", source_grp)

  ! current drive source
  call add_var_int("icd_source",icd_source,0, &
       "1: Include current drive source",source_grp)
  call add_var_double("J_0cd", j_0cd, 0., &
       "amplitude of current drive", source_grp)
  call add_var_double("R_0cd", r_0cd, 0., &
       "R-coordinate of cd maximum", source_grp)
  call add_var_double("Z_0cd", z_0cd, 0., &
       "Z-coordinate of cd maximum", source_grp)
  call add_var_double("W_cd", w_cd, 0., &
       "width of cd gaussian", source_grp)
  call add_var_double("delta_cd", delta_cd, 0., &
       "shift of cd gaussian", source_grp)

  ! poloidal momentum source
  call add_var_int("ipforce", ipforce, 0, &
       "1: Include Poloidal momentum source", source_grp)
  call add_var_double("dforce", dforce, 0., &
       "half-width of poloidal momentum source", source_grp)
  call add_var_double("xforce", xforce, 0., &
       "location [0,1] of poloidal momentum source", source_grp) 
  call add_var_int("nforce", nforce, 0, &
       "exponent of (1-x) multiplying poloidal mom. source", source_grp) 
  call add_var_double("aforce", aforce, 0., &
       "magnitude of poloidal momentum source", source_grp)


  ! gaussian heat source
  call add_var_int("igaussian_heat_source", igaussian_heat_source, 0, &
       "Include gaussian heat source", source_grp)
  call add_var_double("ghs_x", ghs_x, 0., &
       "R coordinate of gaussian heat source", source_grp)
  call add_var_double("ghs_z", ghs_z, 0., &
       "Z coordinate of gaussian heat source", source_grp)
  call add_var_double("ghs_rate", ghs_rate, 0., &
       "Amplitude of gaussian heat source", source_grp)
  call add_var_double("ghs_var", ghs_var, 1., &
       "Variance of gaussian heat source", source_grp)

  call add_var_int("ionization", ionization, 0, "", source_grp)
  call add_var_double("ionization_rate", ionization_rate, 0., "", source_grp)
  call add_var_double("ionization_temp", ionization_temp, 0.01, "", source_grp)
  call add_var_double("ionization_depth", ionization_depth, 0.01, "", source_grp)
  
  call add_var_int("isink", isink, 0, "", source_grp)
  call add_var_double("sink1_x", sink1_x, 0., "", source_grp)
  call add_var_double("sink1_z", sink1_z, 0., "", source_grp)
  call add_var_double("sink1_rate", sink1_rate, 0., "", source_grp)
  call add_var_double("sink1_var", sink1_var, 1., "", source_grp)
  call add_var_double("sink2_x", sink2_x, 0., "", source_grp)
  call add_var_double("sink2_z", sink2_z, 0., "", source_grp)
  call add_var_double("sink2_rate", sink2_rate, 0., "", source_grp)
  call add_var_double("sink2_var", sink2_var, 1., "", source_grp)

  call add_var_int("idenfloor", idenfloor, 0, "", source_grp)
  call add_var_double("alphadenfloor", alphadenfloor, 0., "", source_grp)

  call add_var_double("n_target", n_target, 1., "", source_grp)
  call add_var_double("n_control_p", n_control_p, 0., "", source_grp)
  call add_var_double("n_control_i", n_control_i, 0., "", source_grp)
  call add_var_double("n_control_d", n_control_d, 0., "", source_grp)
  call add_var_int("n_control_type", n_control_type, 0, "", source_grp)
  

  ! Output
  call add_var_int("iprint", iprint, 0, "", output_grp)
  call add_var_int("ntimepr", ntimepr, 5, &
       "Number of time steps per field/restart output", output_grp)
  call add_var_int("iglobalout", iglobalout, 0, "", output_grp)
  call add_var_int("iglobalin", iglobalin, 0, "", output_grp)
  call add_var_int("iwrite_restart", iwrite_restart, 1, &
       "1: Write restart files", output_grp)
  call add_var_int("ifout",  ifout, -1, "", output_grp)
  call add_var_int("icalc_scalars", icalc_scalars, 1, &
       "1: Calculate scalar diagnostics", output_grp)
  call add_var_int("ike_only", ike_only, 0, &
       "1: Only calculate ke scalar diagnostic", output_grp)
  call add_var_int("ike_harmonics", ike_harmonics, 0, &
       "Number of Fourier harmonics of ke to be calculated and output", output_grp)
  call add_var_int("ibh_harmonics", ibh_harmonics, 0, &
       "Number of Fourier harmonics of magnetic perturbation to be calculated and output", output_grp)
  call add_var_int("irestart", irestart, 0, "", output_grp)
  call add_var_int("itimer", itimer, 0, &
       "1: Output internal timer data", output_grp)
  call add_var_int("iwrite_transport_coeffs", iwrite_transport_coeffs, 1, &
       "1: Output transport coefficient fields", output_grp)
  call add_var_int("iwrite_aux_vars", iwrite_aux_vars, 1, &
       "1: Output auxiliary variable fields", output_grp)
  call add_var_int("itemp_plot", itemp_plot, 0, &
       "1: Output additional temperature plots", output_grp)
  call add_var_int("ibdgp", ibdgp, 0, &
       "ne.0: bdgp plot contains only partial results ", output_grp)

  ! diagnostics
  call add_var_int("xray_detector_enabled", xray_detector_enabled, 0, &
       "1: enable xray detector", diagnostic_grp)
  call add_var_double("xray_r0", xray_r0, 0., &
       "R coordinate of xray detector", diagnostic_grp)
  call add_var_double("xray_phi0", xray_phi0, 0., &
       "Phi coordinate of xray detector", diagnostic_grp)
  call add_var_double("xray_z0", xray_z0, 0., &
       "Z coordinate of xray detector", diagnostic_grp)
  call add_var_double("xray_theta", xray_theta, 0., &
       "Angle of xray detector chord (degrees)", diagnostic_grp)
  call add_var_double("xray_sigma", xray_sigma, 1., &
       "Spread of xray detector chord (degrees)", diagnostic_grp)

  ! 3-D options
  call add_var_int("ntor", ntor, 0, &
       "Toroidal mode number", misc_grp)
  call add_var_int("mpol", mpol, 0, "", misc_grp)


  ! Mesh adaptation
  call add_var_int("iadapt", iadapt, 0, "", adapt_grp)
  call add_var_double("adapt_factor", adapt_factor, 1., "", adapt_grp)
  call add_var_double("adapt_hmin", adapt_hmin, 0.001, "", adapt_grp)
  call add_var_double("adapt_hmax", adapt_hmax, 0.1, "", adapt_grp)
  call add_var_double("adapt_smooth", adapt_smooth, 2./3., "", adapt_grp)
  call add_var_double("adapt_psin_vacuum", adapt_psin_vacuum, 0., &
       "", adapt_grp)

  ! Mesh
  call add_var_int("nplanes", nplanes, 1, &
       "Number of toroidal planes", mesh_grp)
  call add_var_double("xzero", xzero, 0., "", mesh_grp)
  call add_var_double("zzero", zzero, 0., "", mesh_grp)
  call add_var_double("tiltangled", tiltangled, 0., "", mesh_grp)
  call add_var_string("mesh_filename", mesh_filename, 256, "struct-dmg.sms", &
       "", mesh_grp)
  call add_var_string("mesh_model", mesh_model, 256, "struct.dmg", &
       "", mesh_grp)
  call add_var_int("ipartitioned",ipartitioned,0,&
       "1 = the input mesh is partitioned", mesh_grp)
  call add_var_int("imatassemble", imatassemble, 0, &
       "0: use scorec matrix parallel assembly; 1 use petsc", mesh_grp)

  call add_var_int("imulti_region", imulti_region, 0, &
       "1 = Mesh has multiple physical regions", mesh_grp)
  
  ! Deprecated
  call add_var_int("ibform", ibform, -1, "", deprec_grp)
  call add_var_int("igs_method", igs_method, -1, "", gs_grp)
end subroutine set_defaults


subroutine validate_input
  use basic
  use mesh_mod
  use m3dc1_nint
  use transport_coefficients
  use neutral_beam
  use pellet
  use math
  use gradshafranov

  implicit none

#include "finclude/petsc.h"
#ifdef PetscDEV
  PetscBool :: flg_petsc, flg_solve2, flg_pdslin
#else
  PetscTruth :: flg_petsc, flg_solve2, flg_pdslin
#endif
  integer :: ier

  if(myrank.eq.0) then
     print *, "============================================="
     print *, " VALIDATING INPUT"
     print *, " ----------------"
  end if
  
!...check if correct code version is being used
#if defined(USE3D)
    if(linear.eq.1 .or. nplanes.le.1) then
      if(myrank.eq.0) print *, "must have linear=0 and nplanes>1 for 3D version"
      call safestop(1)
    endif
#endif

#if defined(USECOMPLEX)
    if(linear.eq.0 .or. nplanes.gt.1) then
      if(myrank.eq.0) print *, "must have linear=1 and nplanes=1 for complex version"
      call safestop(1)
    endif
#endif
    if(linear.eq.0 .and. nplanes.eq.1) then
      ier = 0
#if defined(USECOMPLEX) || (USE3D)
      ier = 1
#endif
      if(ier.ne.0) then
        if(myrank.eq.0) print *,"must use RL version for linear.eq.0 .and. nplanes.eq.1"
        call safestop(1)
      endif
    endif

  if(amuc.eq.0.) amuc = amu

  if(linear.eq.1) then
     eqsubtract = 1
     if(iteratephi.eq.1) then
        if(myrank.eq.0) print *, "iteratephi=1 is not allowed with linear=1."
        call safestop(1)
     endif
  endif

  if(amupar.ne.0 .and. ivform.eq.0) then
     if(myrank.eq.0) print *, "Parallel viscosity not implemented for ivform=0"
     call safestop(1)
  end if

  if(ipressplit.eq.0 .and. itemp.eq.1) then
     if(myrank.eq.0) print *, "itemp=1 not allowed with ipressplit=0"
     call safestop(1)
  endif

  if(isplitstep.eq.0 .and. ipressplit.eq.1) then
     if(myrank.eq.0) print *, "ipressplit=1 not allowed with isplitstep=0"
     call safestop(1)
  endif

  if(isplitstep.eq.0 .and. idiff.eq.1) then
     if(myrank.eq.0) print *, "idiff=1 not allowed with isplitstep=0"
     call safestop(1)
  endif

  if(igs_method.ne.-1) then 
     if(myrank.eq.0) print *, "WARNING: igs_method is now deprecated"
  end if
  
  ! calculate pfac (pe*pfac = electron pressure)
  if(p0.gt.0.) then
     pefac = (p0-pi0)/p0
  else
     pefac = 0.
  endif
  if(myrank.eq.0 .and. iprint.ge.1) print *, "pefac = ", pefac

#if defined(USE3D) || defined(USECOMPLEX)
  i3d = 1
#else
  i3d = 0
#endif

  if(iadapt.gt.0) then
#if defined(USECOMPLEX)
      if(myrank.eq.0) print *, "ERROR:  must use real version of code for iadapt.gt.0"
      call safestop(1)
#endif
   endif

  if(ifout.eq.-1) ifout = i3d
  if(ibform.ne.-1) then
     if(myrank.eq.0) print *, 'WARNING: ibform input parameter deprecated'
  endif
  if(i3d.eq.1 .and. jadv.eq.0) then
     if(myrank.eq.0) &
          print *, 'WARNING: nonaxisymmetric cases should use jadv=1'
  endif

  if(isplitstep.eq.0) then
     imp_mod = 0

     if(hyperc.ne.0.) &
          print *, 'WARNING: poloidal velocity smoothing not available with isplitstep=0'
     if(jadv.eq.1 .and. hyper.ne.0) &
          print *, 'WARNING: poloidal flux smoothing not available with isplitstep=0 and jadv=1'
  end if

  if(i3d.eq.0 .and. imp_bf.ne.0) then
     imp_bf = 0
  endif

  if(rzero.eq.-1) then
     if(itor.eq.1) then 
        rzero = xzero
     else
        rzero = 1.
     endif
  endif

  if(rzero.le.0) then
     if(myrank.eq.0) print *, 'WARNING: rzero <= 0'
  endif

  if(pefac.eq.0. .and. eta0.ne.0 .and. iresfunc.eq.0) then
     if(myrank.eq.0) print *, 'ERROR: Te = 0, but eta0 != 0.'
     call safestop(1)
  endif

  if(amuc.lt.(2./3.)*amu) then
     if(myrank.eq.0) &
          print *, 'ERROR: Constraint amuc >= (2/3)*amu violated.'
     call safestop(1)
  endif
  if(icalc_scalars.eq.0) then
     if(isources.eq.1) then
        if(myrank.eq.0) print *, 'ERROR: isources=1 requires icalcscalars=1'
        call safestop(1)
     endif
  endif
#ifndef USE3D
  int_pts_tor = 1
#endif
  if(int_pts_main*int_pts_tor .gt. MAX_PTS) then
     if(myrank.eq.0) print*, 'ERROR: int_pts_max > MAX_PTS = ', MAX_PTS
     call safestop(1)
  endif
  if(int_pts_aux*int_pts_tor .gt. MAX_PTS) then
     if(myrank.eq.0) print*, 'ERROR: int_pts_aux > MAX_PTS = ', MAX_PTS
     call safestop(1)
  endif
  if(int_pts_diag*int_pts_tor .gt. MAX_PTS) then
     if(myrank.eq.0) print*, 'ERROR: int_pts_diag > MAX_PTS = ', MAX_PTS
     call safestop(1)
  endif
  if((.not.quadrature_implemented(int_pts_main)) .or. &
       (.not.quadrature_implemented(int_pts_aux)) .or. &
       (.not.quadrature_implemented(int_pts_diag))) then 
     if(myrank.eq.0) print*, 'ERROR: integration quadrature not implemented.'
     call safestop(1)
  endif
  if(psiscale.gt.1.) then
     if(myrank.eq.0) print*, 'Warning: psiscale > 1 not supported'
     psiscale = 1.
  endif

  if(integrator.eq.1) then
     thimp = 1.
  endif

  if(itime_independent.eq.1) then
     thimp = 1.
  end if
      
  if(iflip.eq.1) then
     vloop = -vloop
     tcur = -tcur
  endif

  if(igs_forcefree_lcfs.eq.-1) then
     if(iread_eqdsk.ne.0 .and. igs_extend_p.eq.0 .and. irot.le.0) then
        igs_forcefree_lcfs = 2
     else
        igs_forcefree_lcfs = 0
     end if
  end if

!!$  if(eta_wall.ne.0 .and. iconst_bn.eq.1) then
!!$     if(myrank.eq.0) &
!!$          print *, 'Error: eta_wall!=0 is incompatible with iconst_bn==1'
!!$     call safestop(1)
!!$  endif
!!$
!  if(eta_wall.ne.0.) then
!     if(maxrank.gt.1) then
!        print *, 'Error: resistive wall only supported for single-process runs'
!        call safestop(1)
!     end if
!  endif

  if(numvar.eq.1 .and. imp_bf.eq.1) then
     imp_bf = 0
     if(myrank.eq.0) print *, 'WARNING: numvar==1; setting imp_bf = 0'
  end if

  if(ntimepr.lt.1) ntimepr = 1

  if(nplanes.lt.1) nplanes = 1

#ifndef USE3D
  if(nplanes.ne.1) then
     print *, "Compile option '3D=1' must be set to use nplanes>1"
     call safestop(1)
  end if
#endif

#if defined(USE3D) && defined(USEPETSC)
  if(maxrank.ne.nplanes) then 
     print *, 'Must run with procs = nplanes'
     call safestop(1)
  end if
#endif

  if(iread_omega_e .ne. 0) then
     if(iread_omega .ne. 0) then
        if(myrank.eq.0) print *, "Error, can't read multiple rotation profiles"
        call safestop(1)
     end if
     iread_omega = iread_omega_e
  end if
  if(iread_omega_ExB .ne. 0) then
     if(iread_omega .ne. 0) then
        if(myrank.eq.0) print *, "Error, can't read multiple rotation profiles"
        call safestop(1)
     end if
     iread_omega = iread_omega_ExB
  end if

  ! Read PETSc options
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-pdslin', flg_pdslin,ier)
  
  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, "petsc arguments: ipetsc, solve2, solve1", flg_petsc, flg_solve2, flg_pdslin
     print *, "petsc true/false", PETSC_TRUE, PETSC_FALSE
     if(flg_petsc.eq.PETSC_TRUE) print*, 'Using SCOREC PETSc.'
     if(flg_solve2.eq.PETSC_TRUE) print*, 'Using PPPL solve2.'
     if(flg_pdslin.eq.PETSC_TRUE) print*, 'Using PDSLin.'
  endif

  is_rectilinear = (nonrect.eq.0)

  density_source = idens.eq.1 .and. &
       (ipellet.ge.1 .or. ionization.ge.1 .or. isink.gt.0 &
                           .or. idenfloor.gt.0 .or. ibeam.eq.1 .or. ibeam.eq.2)
  momentum_source = (ibeam.eq.1 .or. ibeam.eq.4)
  heat_source = (numvar.ge.3 .or. ipres.eq.1) .and. &
       (igaussian_heat_source.eq.1 .or. ibeam.ge.1)

  if(myrank.eq.0 .and. iprint.ge.1) then 
     print *, 'Density source: ', density_source
     print *, 'Momentum source: ', momentum_source
     print *, 'Heat source: ', heat_source
  end if

  if(den_edge .eq.0) den_edge = den0*(pedge/p0)**expn

  if(irmp.eq.0 .and. iread_ext_field.eq.0) then
     if(extsubtract.ne.0) then
        print *, 'Error: with no external fields, set extsubtract=0'
        call safestop(1)
     end if
  end if

  v0_norm = b0_norm / sqrt(4.*pi*ion_mass*m_p*n0_norm)
  t0_norm = l0_norm / v0_norm
  p0_norm = b0_norm**2/(4.*pi)
  e0_norm = v0_norm*b0_norm / c_light
  
  if(ibeam.ge.1) call neutral_beam_init
  if(ipellet.ge.1) call pellet_init

  if(myrank.eq.0) then
     print *, "============================================="
  end if

end subroutine validate_input


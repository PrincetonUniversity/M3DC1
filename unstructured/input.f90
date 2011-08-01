subroutine add_var_double(name, var, default)
  implicit none

  character*(*), intent(in) :: name
  real :: var
  real, intent(in) :: default

  call add_variable_double(name//char(0), var, default)
end subroutine add_var_double

subroutine add_var_int(name, var, default)
  implicit none

  character*(*), intent(in) :: name
  integer :: var
  integer, intent(in) :: default

  call add_variable_int(name//char(0), var, default)
end subroutine add_var_int


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

  if(myrank.eq.0) print *, " setting defaults"
  call set_defaults

  ! Read input file
  ! ~~~~~~~~~~~~~~~
  if(myrank.eq.0) print *, " reading input"
  call read_namelist("C1input", ierr)

  if(ierr.ne.0) call safestop(3)

  if(myrank.eq.0) call print_namelist

  if(myrank.eq.0 .and. iprint.ge.1) print *, " validating input"
  call validate_input
end subroutine input


!==========================
! set_defaults
! ~~~~~~~~~~~~
! sets defa
!==========================
subroutine set_defaults
  use basic

  implicit none

  ! normalizations
  call add_var_double("b0_norm", b0_norm, 1.e-4)
  call add_var_double("n0_norm", n0_norm, 1.e14)
  call add_var_double("l0_norm", l0_norm, 100.)
  
  ! equilibria input options
  call add_var_int("iread_eqdsk", iread_eqdsk, 0)
  call add_var_int("iread_dskbal", iread_dskbal, 0)
  call add_var_int("iread_jsolver", iread_jsolver, 0)
  call add_var_int("iread_omega", iread_omega, 0)
  call add_var_int("iread_ne", iread_ne, 0)
  call add_var_int("iread_te", iread_te, 0)

  ! transport coefficients
  call add_var_int("ivisfunc", ivisfunc, 0)
  call add_var_double("amuoff", amuoff, 0.)
  call add_var_double("amudelt", amudelt, 0.)
  call add_var_double("amuoff2", amuoff2, 0.)
  call add_var_double("amudelt2", amudelt2, 0.)
  call add_var_double("amu", amu, 0.)
  call add_var_double("amuc", amuc, 0.)
  call add_var_double("amue", amue, 0.)
  call add_var_double("amupar", amupar, 0.)
  call add_var_double("amu_edge", amu_edge, 0.)

  call add_var_int("iresfunc", iresfunc, 0)
  call add_var_double("etaoff", etaoff, 0.)
  call add_var_double("etadelt", etadelt, 0.)
  call add_var_double("etar", etar, 0.)
  call add_var_double("eta0", eta0, 0.)

  call add_var_int("ikappafunc", ikappafunc, 0)
  call add_var_int("ikapscale", ikapscale, 0)
  call add_var_double("kappaoff", kappaoff, 0.)
  call add_var_double("kappadelt", kappadelt, 0.)
  call add_var_double("kappat", kappat, 0.)
  call add_var_double("kappa0", kappa0, 0.)
  call add_var_double("kappar", kappar, 0.)
  call add_var_double("kappax", kappax, 0.)
  call add_var_double("kappah", kappah, 0.)

  call add_var_double("denm", denm, 0.)
  
  call add_var_double("lambdae", lambdae, 0.)

  call add_var_double("gam", gam, 5./3.)
  call add_var_double("db", db, 0.)
  call add_var_double("mass_ratio", mass_ratio, 0.)

  ! model options
  call add_var_int("numvar", numvar, 3)
  call add_var_int("linear", linear, 0)
  call add_var_int("eqsubtract", eqsubtract, 0)
  call add_var_int("idens", idens, 1)
  call add_var_int("ipres", ipres, 0)
  call add_var_int("gyro", gyro, 0)
  call add_var_int("igauge", igauge, 0)
  call add_var_int("inertia", inertia, 1)
  call add_var_int("itwofluid", itwofluid, 1)
  call add_var_int("ibootstrap", ibootstrap, 0)
  call add_var_int("imp_bf", imp_bf, 0)
  call add_var_int("nosig", nosig, 0)
  call add_var_int("itor", itor, 0)

  call add_var_double("gravr", gravr, 0.)
  call add_var_double("gravz", gravz, 0.)
  call add_var_int("istatic", istatic, 0)
  call add_var_int("iestatic", iestatic, 0)
  call add_var_double("chiiner", chiiner, 1.)
  call add_var_int("ieq_bdotgradt", ieq_bdotgradt, 1)
    
  ! time-step options
  call add_var_int("integrator", integrator, 0)
  call add_var_int("isplitstep", isplitstep, 1)
  call add_var_int("iteratephi", iteratephi, 0)
  call add_var_int("imp_mod", imp_mod, 1)
  call add_var_int("irecalc_eta", irecalc_eta, 0)
  call add_var_int("iconst_eta", iconst_eta, 0)
  call add_var_int("itime_independent", itime_independent, 0)
  call add_var_double("thimp", thimp, 0.5)
  call add_var_double("thimpsm", thimpsm, 1.)
  call add_var_double("harned_mikic", harned_mikic, 0.)
  call add_var_int("isources", isources, 0)
  call add_var_int("nskip", nskip, 1)
  call add_var_double("dt", dt, 0.1)
  call add_var_double("ddt", ddt, 0.)
  
  ! numerical options
  call add_var_int("jadv", jadv, 0)
  call add_var_int("ivform", ivform, 0)
  call add_var_int("ibform", ibform, -1)
  call add_var_int("int_pts_main", int_pts_main, 25)
  call add_var_int("int_pts_aux", int_pts_aux, 25)
  call add_var_int("int_pts_diag", int_pts_diag, 25)
  call add_var_int("int_pts_tor", int_pts_tor, 5)
  call add_var_double("max_ke", max_ke, 1.)
  call add_var_int("equilibrate", equilibrate, 0)
  call add_var_double("regular", regular, 0.)

  ! equilibrium options
  call add_var_int("itaylor", itaylor, 0)
  call add_var_int("iupstream", iupstream, 0)
  call add_var_int("iflip", iflip, 0)
  call add_var_int("iflip_b", iflip_b, 0)
  call add_var_int("iflip_j", iflip_j, 0)
  call add_var_int("iflip_v", iflip_v, 0)
  call add_var_int("iflip_z", iflip_z, 0)
  call add_var_int("icsym", icsym, 0)
  call add_var_double("bzero", bzero, 1.)
  call add_var_double("bx0", bx0, 0.)
  call add_var_double("vzero", vzero, 0.)
  call add_var_double("phizero", phizero, 0.)
  call add_var_int("idevice", idevice, 0)
  call add_var_int("iwave", iwave, 0)
  call add_var_double("eps", eps, 0.01)
  call add_var_int("maxn", maxn, 200)
  call add_var_int("irmp", irmp, 0)
  call add_var_double("beta", beta, 0.)
  call add_var_double("ln", ln, 0.)

  ! grad-shafranov options
  call add_var_int("inumgs", inumgs, 0)
  call add_var_int("igs", igs, 80)
  call add_var_int("igs_method", igs_method, 2)
  call add_var_int("nv1equ", nv1equ, 0)
  call add_var_double("tcuro", tcuro, 1.)
  call add_var_double("xmag", xmag, 1.)
  call add_var_double("zmag", zmag, 0.)
  call add_var_double("xlim", xlim, 0.)
  call add_var_double("zlim", zlim, 0.)
  call add_var_double("xlim2", xlim2, 0.)
  call add_var_double("zlim2", zlim2, 0.)
  call add_var_double("rzero", rzero, -1.)
  call add_var_double("libetap", libetap, 1.2)
  call add_var_double("p0", p0, 0.01)
  call add_var_double("pi0", pi0, 0.005)
  call add_var_double("p1", p1, 0.)
  call add_var_double("p2", p2, 0.)
  call add_var_double("pedge", pedge, -1.) ! If pedge < 0, don't use pedge
  call add_var_double("tedge", tedge, -1.) ! ditto for tedge
  call add_var_double("expn", expn, 0.)
  call add_var_double("q0", q0, 1.)
  call add_var_double("djdpsi", djdpsi, 0.)
  call add_var_double("th_gs", th_gs, 0.8)
  call add_var_double("tol_gs", tol_gs, 1.e-8)
  call add_var_double("psiscale", psiscale, 1.)
  call add_var_double("pscale", pscale, 1.)
  call add_var_double("bscale", bscale, 1.)

  call add_var_int("irot", irot, 0)
  call add_var_int("iscale_rot_by_p", iscale_rot_by_p, 1)
  call add_var_double("alpha0", alpha0, 0.)
  call add_var_double("alpha1", alpha1, 0.)
  call add_var_double("alpha2", alpha2, 0.)
  call add_var_double("alpha3", alpha3, 0.)
  
  call add_var_int("idenfunc", idenfunc, 0)
  call add_var_double("den_edge", den_edge, 0.)
  call add_var_double("den0", den0, 1.)
  call add_var_double("dendelt", dendelt, 0.1)
  call add_var_double("denoff", denoff, 1.)

  call add_var_int("divertors", divertors, 0)
  call add_var_double("xdiv", xdiv, 0.)
  call add_var_double("zdiv", zdiv, 0.)
  call add_var_double("divcur", divcur, 0.1)

  call add_var_double("xnull", xnull, 0.)
  call add_var_double("znull", znull, 0.)

  ! hyper-diffusivity
  call add_var_double("deex", deex, 1.)
  call add_var_double("hyper", hyper, 0.)
  call add_var_double("hyperi", hyperi, 0.)
  call add_var_double("hyperv", hyperv, 0.)
  call add_var_double("hyperp", hyperp, 0.)
  call add_var_double("hyperv", hyperv, 0.)
  call add_var_double("hyperc", hyperc, 0.)
  call add_var_int("ihypdx", ihypdx, 2)
  call add_var_int("ihypeta", ihypeta, 1)
  call add_var_int("ihypamu", ihypamu, 1)
  call add_var_int("ihypkappa", ihypkappa, 1)

  ! boundary conditions
  call add_var_int("isurface", isurface, 1)
  call add_var_int("icurv", icurv, 2)
  call add_var_int("nonrect", nonrect, 0)
  call add_var_int("ifixedb", ifixedb, 0)
  call add_var_int("com_bc", com_bc, 0)
  call add_var_int("vor_bc", vor_bc, 0)
  call add_var_int("iconst_p", iconst_p, 1)
  call add_var_int("iconst_n", iconst_n, 1)
  call add_var_int("iconst_t", iconst_t, 0)
  call add_var_int("iconst_bz", iconst_bz, 1)
  call add_var_int("inograd_p", inograd_p, 0)
  call add_var_int("inograd_n", inograd_n, 0)
  call add_var_int("inonormalflow", inonormalflow, 1)
  call add_var_int("inoslip_pol", inoslip_pol, 1)
  call add_var_int("inoslip_tor", inoslip_tor, 1)
  call add_var_int("inostress_tor", inostress_tor, 0)
  call add_var_int("inocurrent_pol", inocurrent_pol, 0)
  call add_var_int("inocurrent_tor", inocurrent_tor, 0)
  call add_var_int("inocurrent_norm", inocurrent_norm, 0)
  call add_var_int("ifbound", ifbound, 1)
  call add_var_int("iconstflux", iconstflux, 0)
  
  ! resistive wall
  call add_var_double("eta_wall", eta_wall, 0.)
  call add_var_double("delta_wall", delta_wall, 1.)

  ! loop voltage
  call add_var_double("vloop", vloop, 0.)
  call add_var_double("tcur", tcur, 0.)
  call add_var_double("control_p", control_p, 0.)
  call add_var_double("control_i", control_i, 0.)
  call add_var_double("control_d", control_d, 0.)
  
  ! density source
  call add_var_int("ipellet", ipellet, 0)
  call add_var_double("pellet_x", pellet_x, 0.)
  call add_var_double("pellet_z", pellet_z, 0.)
  call add_var_double("pellet_rate", pellet_rate, 0.)
  call add_var_double("pellet_var", pellet_var, 1.)

  call add_var_int("ionization", ionization, 0)
  call add_var_double("ionization_rate", ionization_rate, 0.)
  call add_var_double("ionization_temp", ionization_temp, 0.01)
  call add_var_double("ionization_depth", ionization_depth, 0.01)
  
  call add_var_int("isink", isink, 0)
  call add_var_double("sink1_x", sink1_x, 0.)
  call add_var_double("sink1_z", sink1_z, 0.)
  call add_var_double("sink1_rate", sink1_rate, 0.)
  call add_var_double("sink1_var", sink1_var, 1.)
  call add_var_double("sink2_x", sink2_x, 0.)
  call add_var_double("sink2_z", sink2_z, 0.)
  call add_var_double("sink2_rate", sink2_rate, 0.)
  call add_var_double("sink2_var", sink2_var, 1.)

  call add_var_double("n_target", n_target, 1.)
  call add_var_double("n_control_p", n_control_p, 0.)
  call add_var_double("n_control_i", n_control_i, 0.)
  call add_var_double("n_control_d", n_control_d, 0.)
  
  ! I/O options
  call add_var_int("iprint", iprint, 0)
  call add_var_int("ntimemax", ntimemax, 20)
  call add_var_int("ntimepr", ntimepr, 5)
  call add_var_int("iglobalout", iglobalout, 0)
  call add_var_int("iglobalin", iglobalin, 0)
  call add_var_int("iwrite_restart", iwrite_restart, 1)
  call add_var_int("ifout",  ifout, -1)
  call add_var_int("icalc_scalrs", icalc_scalars, 1)
  call add_var_int("ike_only", ike_only, 0)
  call add_var_int("irestart", irestart, 0)
  call add_var_int("itimer", itimer, 0)

  ! 3-D options
  call add_var_int("nplanes", nplanes, 1)
  call add_var_int("ntor", ntor, 0)
  call add_var_int("mpol", mpol, 0)

  ! adaptation options
  call add_var_int("iadapt", iadapt, 0)
  call add_var_double("adapt_factor", adapt_factor, 1.)
  call add_var_double("adapt_hmin", adapt_hmin, 0.001)
  call add_var_double("adapt_hmax", adapt_hmax, 0.1)

  ! mesh options
  call add_var_int("iper", iper, 0)
  call add_var_int("jper", jper, 0)
  call add_var_double("xzero", xzero, 0.)
  call add_var_double("zzero", zzero, 0.)
  call add_var_double("tiltangled", tiltangled, 0.)
end subroutine set_defaults


subroutine validate_input
  use basic
  use mesh_mod
  use m3dc1_nint

  implicit none

#include "finclude/petsc.h"
#ifdef PetscDEV
  PetscBool :: flg_petsc, flg_solve2, flg_pdslin
#else
  PetscTruth :: flg_petsc, flg_solve2, flg_pdslin
#endif
  integer :: ier

  if(amuc.eq.0.) amuc = amu

  if(linear.eq.1) then
     eqsubtract = 1
     if(iteratephi.eq.1) then
        if(myrank.eq.0) print *, "iteratephi=1 is not allowed with linear=1."
        call safestop(1)
     endif
  endif
  
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

#if defined(USECOMPLEX)
    if(iadapt.gt.0) then
      if(myrank.eq.0) print *, "ERROR:  must use real version of code for iadapt.gt.0"
      call safestop(1)
    endif
#endif
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

  if(rzero.eq.-1) then
     if(itor.eq.1) then 
        rzero = xzero
     else
        rzero = 1.
     endif
  endif

  if(rzero.le.0) then
     print *, 'WARNING: rzero <= 0'
  endif

  if(pefac.eq.0. .and. eta0.ne.0) then
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
  if(int_pts_main .gt. MAX_PTS) then
     if(myrank.eq.0) print*, 'ERROR: int_pts_max > MAX_PTS = ', MAX_PTS
     call safestop(1)
  endif
  if(int_pts_aux .gt. MAX_PTS) then
     if(myrank.eq.0) print*, 'ERROR: int_pts_aux > MAX_PTS = ', MAX_PTS
     call safestop(1)
  endif
  if(int_pts_diag .gt. MAX_PTS) then
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
      
  if(iflip.eq.1) then
     vloop = -vloop
     tcur = -tcur
  endif

#ifdef USE3D
  if(gyro.eq.1) then
     print *, 'Error: gyroviscosity not yet implemented in 3D'
     call safestop(1)
  endif
#endif

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

  ! Read PETSc options
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-pdslin', flg_pdslin,ier)
  
  if(myrank.eq.0) then
     print *, "petsc arguments: ipetsc, solve2, solve1", flg_petsc, flg_solve2, flg_pdslin
     print *, "petsc true/false", PETSC_TRUE, PETSC_FALSE
     if(flg_petsc.eq.PETSC_TRUE) print*, 'Using SCOREC PETSc.'
     if(flg_solve2.eq.PETSC_TRUE) print*, 'Using PPPL solve2.'
     if(flg_pdslin.eq.PETSC_TRUE) print*, 'Using PDSLin.'
  endif

  is_rectilinear = (nonrect.eq.0)

  if(myrank.eq.0) call print_namelist

end subroutine validate_input


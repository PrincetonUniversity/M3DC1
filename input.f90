!=========================
! input
! ~~~~~
! reads input namelist
!=========================
subroutine input
  use basic

  implicit none

  if(myrank.eq.0 .and. iprint.eq.1) print *, " setting defaults"
  call set_defaults

  ! Read input file
  ! ~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.eq.1) print *, " reading input"
  open(5,file='C1input',form='formatted',status='old')
  read(5,nml=inputnl)
  close(5)
 
  if(myrank.eq.0 .and. iprint.eq.1) print *, " validating input"
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
  b0_norm = 1.e4
  n0_norm = 1.e14
  l0_norm = 100.
  
  ! equilibria input options
  iread_eqdsk = 0
  iread_dskbal = 0
  iread_jsolver = 0

  ! transport coefficients
  ivisfunc = 0
  amuoff = 0.
  amudelt = 0.
  amuoff2 = 0.
  amudelt2 = 0.
  amu = 0.
  amuc = 0.
  amue = 0.
  amupar = 0.
  amu_edge = 0.

  iresfunc = 0
  etaoff = 0.
  etadelt = 0.
  etar = 0.
  eta0 = 0.

  ikappafunc = 0
  kappaoff = 0.
  kappadelt = 0.
  kappat = 0.
  kappa0 = 0.

  kappar = 0.
  kappax = 0.
  kappah = 0.
  ikapscale = 0

  lambdae = 0.

  gam = 5./3.
  db = 0.
  mass_ratio = 0.

  ! model options
  numvar = 3
  linear = 0
  idens = 1
  ipres = 0
  igauge = 0
  inertia = 1
  itwofluid = 1
  ibootstrap = 0
  nosig = 0
  itor = 0
  gravr = 0.
  gravz = 0.
  istatic = 0
  iestatic = 0
  chiiner = 1.
  
  ! time-step options
  integrator = 0
  isplitstep = 1
  iteratephi = 0
  imp_mod = 1
  irecalc_eta = 0
  iconst_eta = 0
  thimp = 0.5
  thimp_ohm = -1.
  thimpsm = 1.
  harned_mikic = 0.
  isources = 0
  nskip = 1
  dt = 0.1
  ddt = 0.
  
  ! numerical options
  jadv = 0
  ivform = 0
  ibform = -1
  int_pts_main = 79
  int_pts_aux = 79
  int_pts_diag = 79
  max_ke = 1.

  ! equilibrium options
  itaylor = 0
  iflip = 0
  iflip_b = 0
  iflip_v = 0
  icsym = 0
  bzero = 1.
  bx0 = 0.
  vzero = 0.
  phizero = 0.
  idevice = 0
  iwave = 0
  eps = .01
  maxn = 200
  irmp = 0
  bscale = 1.

  ! grad-shafranov options
  inumgs = 0
  igs = 80
  igs_method = 1
  nv1equ = 0
  tcuro = 1.
  xmag = 1.
  zmag = 0.
  xlim = 0.
  zlim = 0.
  rzero = -1.
  xlim2 = 0.
  zlim2 = 0.
  libetap = 1.2
  p0 = 0.01      
  p1 = -1.
  p2 = -2.
  pedge = 0.
  expn = 0.
  q0 = 1.
  djdpsi = 0.0
  th_gs = 0.8
  tol_gs = 1.e-8
  psiscale = 1.

  irot = 0.
  alpha0 = 0.
  alpha1 = 0.
  alpha2 = 0.

  idenfunc = 0
  den_edge = 0.
  den0 = 1.
  dendelt = 0.1
  denoff = 1.

  divertors = 0
  xdiv = 0.
  zdiv = 0.
  divcur = 0.1

  xnull = 0.
  znull = 0.


  ! hyper-diffusivity
  deex = 1.
  hyper = 0.0
  hyperi= 0.0
  hyperv= 0.0
  hyperp= 0.0
  ihypdx = 2
  ihypeta = 1
  ihypamu = 1
  ihypkappa = 1


  ! boundary conditions
  isurface = 1
  icurv = 2
  nonrect = 0
  ifixedb = 0
  com_bc = 0
  vor_bc = 0
  iconst_p = 1
  iconst_n = 1
  iconst_t = 0
  iconst_bz = 1
  inograd_p = 0
  inograd_n = 0
  inonormalflow = 1
  inoslip_pol = 1
  inoslip_tor = 1
  inostress_tor = 0
  inocurrent_pol = 1
  inocurrent_tor = 0
  ifbound = 1
  iconstflux = 0
  
  ! resistive wall
  eta_wall = 0.
  delta_wall = 1.

  ! loop voltage
  vloop = 0.
  tcur = 0.
  control_p = 0.
  control_i = 0.
  control_d = 0.
  
  ! density source
  ipellet = 0
  pellet_x = 0.
  pellet_z = 0.
  pellet_rate = 0. 
  pellet_var = 1.

  ionization = 0
  ionization_rate = 0.
  ionization_temp = 0.01
  ionization_depth = 0.01

  isink = 0
  sink1_x = 0.
  sink1_z = 0.
  sink2_rate = 0. 
  sink2_var = 1.
  sink2_x = 0.
  sink2_z = 0.
  sink2_rate = 0. 
  sink2_var = 1.

  n_target = 1.
  n_control_p = 0.
  n_control_i = 0.
  n_control_d = 0.
  
  ! I/O options
  ntimemax = 20
  ntimepr   = 5
  iglobalout = 0
  iglobalin = 0
  iwrite_restart = 1
  ifout = -1
  icalc_scalars = 1
  ike_only = 0
  irestart = 0

  ! 3-D options
  ntor = 0

  ! adaptation options
  iadapt = 0
  adapt_factor = 1.
  adapt_hmin = 0.001
  adapt_hmax = 0.1

  ! mesh options
  xzero = 0.
  zzero = 0.
  tiltangled = 0.

  imask = 0

end subroutine set_defaults


subroutine validate_input
  use basic
  use nintegrate_mod

  implicit none

#include "finclude/petsc.h"
  PetscTruth :: flg_petsc, flg_solve2
  integer :: ier

  if(amuc.eq.0.) amuc = amu
  if(thimp_ohm.lt.0) thimp_ohm = thimp

  if(linear.eq.1) then
     eqsubtract = 1
     if(iteratephi.eq.1) then
        if(myrank.eq.0) print *, "iteratephi=1 is not allowed with linear=1."
        call safestop(1)
     endif
  endif
  
  ! calculate pfac (pe*pfac = electron pressure)
  if(ipres.eq.1) then
     pefac = 1.
  else
     if(p0.gt.0.) then
        pefac = (p0-pi0)/p0
     else
        pefac = 0.
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "pefac = ", pefac
  endif

  if(icomplex.eq.1) then
     i3d = 1
  else 
     i3d = 0
  endif
  if(ifout.eq.-1) ifout = i3d
  if(ibform.ne.-1) then
     if(myrank.eq.0) print *, 'WARNING: ibform input parameter deprecated'
  endif
  if(i3d.eq.1 .and. jadv.eq.0) then
     if(myrank.eq.0) &
          print *, 'WARNING: nonaxisymmetric cases should use jadv=1'
  endif

  if(isplitstep.eq.0) imp_mod = 0
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
     thimp_ohm = 1.
  endif
      
  if(iflip.eq.1) then
     vloop = -vloop
     tcur = -tcur
  endif

  if(ntimepr.lt.1) ntimepr = 1

  ! Read PETSc options
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  
  if(myrank.eq.0) then
     print *, "petsc arguments: ipetsc, solve2", flg_petsc, flg_solve2
     print *, "petsc true/false", PETSC_TRUE, PETSC_FALSE
     if(flg_petsc.eq.PETSC_TRUE) print*, 'Using PETSc.'
     if(flg_solve2.eq.PETSC_TRUE) print*, 'Using solve2.'
  endif

  if(flg_petsc.eq.PETSC_TRUE .and. flg_solve2.eq.PETSC_TRUE) then 
     drop_zeroes = 0 
  else 
     select case(nonrect)
     case(0)
        drop_zeroes = 1
     case(1)
        drop_zeroes = 0 
     end select
  endif

  if(myrank.eq.0) write(*,nml=inputnl)

end subroutine validate_input


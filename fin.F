      subroutine rinvdef(rinv)

      use basic

      implicit none
      integer i, numnodes, ibegin, iendplusone
      real, intent(out) :: rinv(*)

      real :: x, z

      call numnod(numnodes)
      do i=1,numnodes
         call nodcoord(i, x, z)

         call entdofs(1, i, 0, ibegin, iendplusone)
         rinv(ibegin) =  1./x
         rinv(ibegin+1) = -1./x**2
         rinv(ibegin+3) =  2./x**3
      enddo

      return
      end

!============================================================
      subroutine input
      use basic
      use arrays
      use nintegrate_mod

      implicit none
#include "finclude/petsc.h"

      PetscTruth :: flg_petsc, flg_solve2
      integer :: ier

      if(myrank.eq.0 .and. iprint.eq.1) print *, " entering input.."

      icurv = 2

      b0_norm = 0.
      n0_norm = 0.
      l0_norm = 0.

      iread_eqdsk = 0
      iread_dskbal = 0
      iread_jsolver = 0
!
!.....linear parameter  0-nonlinear,  1-linear
      linear = 0
!
!.....density advance parameter: 0 no advance,  1-advance density
      idens = 1
!
!.....electron pressure advance parameter: 0 no advance, 1-advance pe
      ipres = 0
!
!.....use current density equation instead of poloidal flux equation
      jadv = 0
!
!.....restart parameter 0-no restart, 1 restart
      irestart = 0
!
!.....associated restart parameter: 0-normal, 1-start time from zero
      istart = 0
!
!.....maximum KE before solution is re-scaled in linear case
      max_ke = 1.
!
!.....boundary parameter  0-nonperiodic, 1-periodic
!      iper = 0
!      jper = 0
!
!.....switch for Taylor Problem
      itaylor = 0
!
!.....time integrator: 0 = Crank-Nicholson, 1 = BDF2
      integrator = 0
!
!.....do timestep splitting
      isplitstep = 1
!
!.....iterate field solve
      iteratephi = 0
      icsym = 0
      inumgs = 0
!
!.....iterate field solve
      irecalc_eta = 0
!
!.....modifications to implicitization
      imp_mod = 1
      igauge = 0
!
!.....form of toroidal velocity
      ivform = 0
!
!.....flip handedness
      iflip = 0
!
!.....form of toroidal field equation
      ibform = -1
!
!.....number of grad-shafranov iterations
      igs = 80
      igs_method = 1
!
!.....switch to use numvar equilibrium when numvar > 1
      nv1equ = 0
!
!.....CN implicitness parameter (thimp > 0.5 for stability)
      thimp = 0.5
      thimp_ohm = -1.
      thimpsm = 1.

      harned_mikic = 0.
!
!.....hyper-resistivity coefficients
      deex = 1.
      hyper = 0.0
      hyperi= 0.0
      hyperv= 0.0
      hyperp= 0.0
!
!.....scale hyper-resistivity with dx**ihypdx
      ihypdx = 2
!
!.....scale hyper-resistivity with eta
      ihypeta = 1
!
!.....scale hyper-viscosity with visc
      ihypamu = 1
!
!.....don't evolve resistivity
      iconst_eta = 0
!
!.....scale kappar with kappa
      ikapscale = 0
!
      inertia = 1
      itwofluid = 1
!
      isurface = 0
!
!.....bootstrap current model (0 = none)
      ibootstrap = 0
!
!.....viscosity, resistivity, and heat conduction
      amu = 0.
      amuc = 0.
      amue = 0.
      amupar = 0.
      iresfunc = 0     ! use old definition of resistivity
      etaoff = 0.
      etadelt = 0.
      etar = 0.
      eta0   = 0.
      lambdae = 0.
      kappat = 0.
      kappar = 0.
      kappa0 = 0.
      kappax = 0.
      kappah = 0.
      mass_ratio = 0.
!
      ivisfunc = 0
      amuoff = 0.
      amudelt = 0.
      amuoff2 = 0.
      amudelt2 = 0.
      ikappafunc = 0
      kappaoff = 0.
      kappadelt = 0.
!
!
!.....ratio of specific heats
      gam = 5./3.
!
!.....loop voltage
      vloop = 0.
!
!......toroidal magnetic field
      bzero = 1.
      bx0 = 0.
!
!.....initial toroidal velocity
      vzero = 0.
      phizero = 0.
!
!.....rectangular mesh or curved mesh
      nonrect = 0
!
!.....fixed boundary option
      ifixedb = 0
!
!.....boundary conditions
!     com_bc: 1 = sets laplacian(chi) = 0 on boundary
      com_bc = 1
      iconst_p = 1
      iconst_n = 1
      iconst_t = 0
      iconst_bz = 1
      inograd_p = 0
      inograd_n = 0
      inonormalflow = 1
      inoslip_pol = 0
      inoslip_tor = 1
      inostress_tor = 0
      inocurrent_pol = 1
      inocurrent_tor = 1
!
!.....density source
      ipellet = 0
      pellet_x = 0.
      pellet_z = 0.
      pellet_rate = 0. 
      pellet_var = 1.
!
      ionization = 0
      ionization_rate = denm
      ionization_temp = 0.01
      ionization_depth = 0.01
!
      isink = 0
      sink1_x = 0.
      sink1_z = 0.
      sink2_rate = 0. 
      sink2_var = 1.
      sink2_x = 0.
      sink2_z = 0.
      sink2_rate = 0. 
      sink2_var = 1.

!
      nosig = 0
!
!.....ion skin depth
      db = 0.0
!
!     masking switch (0 no mask,  1 mask hyper terms at boundary
      imask = 0
!
      amu_edge = 0.
!
      isources = 0
!
!.....how many cycles to skip before inverting matrix
      nskip = 1
!
!.....timestep
      dt = 0.1
      ddt = 0.
!.....ntimemax is max time cycles, nimepr is cycles between print cycle
!
      ntimemax = 20
      ntimepr   = 5
!
!.....itor=0 for cylinder, 1 for torus
      itor=0
      idevice = 0
      iwave = 0
!
!.....parameters needed for toroidal plasma
      tcuro = 1.
      xmag = 1.
      zmag = 0.
      xlim = 0.
      zlim = 0.
      rzero = -1.
      xlim2 = 0.
      zlim2 = 0.
      libetap = 1.2
!
!..... toroidal current controller parameters
      tcur = 0.
      control_p = 0.
      control_i = 0.
      control_d = 0.
!..... density controller parameters
      n_target = 1.
      n_control_p = 0.
      n_control_i = 0.
      n_control_d = 0.

!
!.....pressure = p0*(1 + p1*psi + p2*psi**2 + ...)
!                where psi is the normalized poloidal flux
      p0 = 0.01      
      p1 = -1.
      p2 = -2.
      pedge = 0.
      expn = 0.
      q0 = 1.
      th_gs = 0.5
      tol_gs = 1.e-8

      idenfunc = 0
      den_edge = 0.
      den0 = 1.
      dendelt = 0.1
      denoff = 1.

!....equilibrium rotation parameters
      irot = 0.
      alpha0 = 0.
      alpha1 = 0.
      alpha2 = 0.
!
!.....divertor information
      divertors = 0
      xdiv = xmag
      zdiv = zzero-1.
      divcur = 0.1
!
!.....approximate position of x-point
      xnull = 0.
      znull = 0.
!
!     derivative of current wrt normalized poloidal flux
      djdpsi = 0.0
!
!....number of velocity and field variables
      numvar = 3
!
!.....parameters for taylor problem (itaylor=1)
      eps = .01

      maxn = 200

      irmp = 0
!
!.....gravitational acceleration
      gravr = 0.
      gravz = 0.
!
!.....1 = don't advance velocity
      istatic = 0
      iestatic = 0
!
!.....1 = conserve toroidal flux
      iconstflux = 0
!
!.....Toroidal mode number
      ntor = 0
!
!.....1 = output/input global restart files
      iglobalout = 0
      iglobalin = 0
      iwrite_restart = 1
!     
!.....1 = output f field
      ifout = -1
!
      icalc_scalars = 1
      ike_only = 0
!
      int_pts_main = 79
      int_pts_aux = 79
      int_pts_diag = 79
!
!.....1 = adapt mesh upon initialization
      iadapt = 0
      adapt_factor = 1.
      adapt_hmin = 0.001
      adapt_hmax = 0.1
!
!......factor multiplying chi inertial term
      chiiner = 1.
!
!......tilt angle of rectangular mesh (in degrees)
      tiltangled = 0.
!
!.....Bateman scaling factor
      bscale = 1.

!
      open(5,file='C1input',form='formatted',status='old')
      read(5,nml=inputnl)
      close(5)
      if(amuc.eq.0.) amuc = amu
      if(thimp_ohm.lt.0) thimp_ohm = thimp
!
      if(linear.eq.1) then
         eqsubtract = 1
         if(iteratephi.eq.1) then
            if(myrank.eq.0) 
     &           print *, "iteratephi=1 is not allowed with linear=1."
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




      ! Validate input
      ! ~~~~~~~~~~~~~~

      if(icomplex.eq.1) then
         i3d = 1
      else 
         i3d = 0
      endif
      if(ifout.eq.-1) ifout = i3d
      if(ibform.ne.-1) then
         if(myrank.eq.0) 
     &        print *, 'WARNING: ibform input parameter deprecated'
      endif
      if(i3d.eq.1 .and. jadv.eq.0) then
         if(myrank.eq.0)
     &        print *, 'WARNING: nonaxisymmetric cases should &
     &        use jadv=1'
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
         if(myrank.eq.0) 
     &        print *, 'ERROR: Constraint amuc >= (2/3)*amu violated.'
         call safestop(1)
      endif
      if(icalc_scalars.eq.0) then
         if(isources.eq.1) then
            if(myrank.eq.0) 
     &           print *, 'ERROR: isources=1 requires icalcscalars=1'
            call safestop(1)
         endif
      endif
      if(int_pts_main .gt. MAX_PTS) then
         if(myrank.eq.0) 
     &        print*, 'ERROR: int_pts_max > MAX_PTS = ', MAX_PTS
         call safestop(1)
      endif
      if(int_pts_aux .gt. MAX_PTS) then
         if(myrank.eq.0) 
     &        print*, 'ERROR: int_pts_aux > MAX_PTS = ', MAX_PTS
         call safestop(1)
      endif
      if(int_pts_diag .gt. MAX_PTS) then
         if(myrank.eq.0) 
     &        print*, 'ERROR: int_pts_diag > MAX_PTS = ', MAX_PTS
         call safestop(1)
      endif
      if((.not.quadrature_implemented(int_pts_main)) .or.
     &     (.not.quadrature_implemented(int_pts_aux)) .or.
     &     (.not.quadrature_implemented(int_pts_diag))) then 
         if(myrank.eq.0) 
     &        print*, 'ERROR: integration quadrature not implemented.'
         call safestop(1)
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
      call PetscOptionsHasName(PETSC_NULL_CHARACTER,
     &     '-ipetsc', flg_petsc,ier)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER,
     &     '-solve2', flg_solve2,ier)

      if(myrank.eq.0) then
         print *, "petsc arguments: ipetsc, solve2",
     &        flg_petsc, flg_solve2
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

      if(myrank.eq.0 .and. iprint.eq.1) print *, " exiting input.."

      return
      end
!============================================================
      real function fact(n)

      implicit none

      integer, intent(in) :: n
      integer :: i
      real :: ans

! acbauer - my implict none changes here caused problems, probably from real double precision integer things
!
!.....fact calculates the factorial of n...an integer
      ans = 1.
      if (n.le.1) go to 100
      do i=1,n
        ans = ans*i
      enddo
 100  fact = ans
      return
      end
!============================================================
      subroutine scaleback
!
      use basic
      use arrays
      use diagnostics

      implicit none

      real, parameter :: scalefac = 1.e-10

!.....scale solution back if a linear calculation

      if(ekin.lt.max_ke .or. max_ke.eq.0) return
      if(myrank.eq.0) write(*,*) " =>solution scaled back at time", time

      field = scalefac*field
      fieldi = scalefac*fieldi
      phiold = scalefac*phiold
      if(isplitstep.eq.1) then
         velold = scalefac*velold
         veloldn = scalefac*veloldn
         if(idens.eq.1) denold = scalefac*denold
         if(ipres.eq.1) presold = scalefac*presold
      endif
      if(i3d.eq.1) bf = scalefac*bf

      end
!==================================
      real function sech(x)
      implicit none
      real, intent(in) :: x
      
      sech = 1./cosh(x)
      return
      end
!============================================================








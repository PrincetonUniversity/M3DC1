module p_data

  implicit none
  
  integer :: maxi     ! highest degree polynomial kept in nonlinear calculations
  integer :: ires     ! linear resolution of the plot files
  integer :: ntri
  integer :: maxplots ! maximum dimension of the graph array
  integer :: ijacobian

  parameter(maxplots=50, maxi=20, ires=201, ijacobian=1)
  
  integer, dimension(ires, ires) :: whichtri

end module p_data

module basic
  use p_data
  use pid_controller

#ifdef USECOMPLEX
  integer, parameter :: icomplex = 1
  integer, parameter :: i3d = 1
#else
  integer, parameter :: icomplex = 0
  integer, parameter :: i3d = 0
#endif

  ! transport coefficients
  real :: amu         ! incompressible viscosity
  real :: amuc        ! compressible viscosity
  real :: amupar      ! parallel viscosity coefficient
  real :: etar, eta0  ! resistivity = etar + eta0/T^(3/2)
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

  ! domain parameters
  real :: xzero, zzero  ! cooridinates of lower left corner of domain
  real :: rzero    ! nominal major radius of the device
  real :: libetap   !  li/2 + beta_poloidal estimate for equili/2brium calculation
  real :: xlim2   !  x-position of second limiter point as a diagnostic
  real :: zlim2   !  z-position of second limiter point as a diagnostic

  ! boundary conditions
  integer :: iper, jper ! periodic boundary conditions
  integer :: imask      ! 1 = ignore 2-fluid terms near boundaries
  integer :: inonormalflow ! 1 = no normal flow
  integer :: inoslip_pol   ! 1 = no slip (poloidal flow)
  integer :: inoslip_tor   ! 1 = no slip (toroidal flow)
  integer :: inostress_tor ! 1 = no stress (toroidal flow)
  integer :: inocurrent_pol! 1 = no tangential current
  integer :: inocurrent_tor! 1 = no toroidal current
  integer :: iconst_t      ! 1 = temperature held constant
  integer :: iconst_bz     ! 1 = toroidal field held constant
  integer :: inograd_t     ! 1 = no normal temperature gradient
  integer :: com_bc   ! 1 = forces div(V) = 0 on boundary
  real :: amu_edge    ! factor by which to increase viscosity at boundaries

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
  integer :: itaylor  ! equilibrium
  integer :: idevice  ! for itor=1, itaylor=1, selects tokamak configuration
                      !  0 = generic
                      !  1 = CDX-U
                      !  2 = NSTX
  real :: bzero       ! guide field
  real :: vzero       ! initial toroidal velocity
  real :: p0, pi0     ! total, ion pressures
  real :: ln          ! length of equilibrium gradient
  real :: eps         ! size of initial perturbation

  ! grad-shafranov options
  integer :: divertors! number of divertors
  integer :: igs      ! number of grad-shafranov iterations
  integer :: nv1equ   ! if set to 1, use numvar equilibrium for numvar > 1
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
  real :: tol_gs  ! error tolorance for GS solver

  ! model options
  integer :: linear      ! 1 = linear simulation; 0 = nonlinear simulation
  integer :: eqsubtract  ! 1 = subtract equilibrium in noninear simulations
  integer :: numvar      ! 1 = 2-field; 2 = reduced MHD; 3 = compressible MHD
  integer :: idens       ! evolve density
  integer :: ipres       ! evolve total and electron pressures separately
  integer :: gyro        ! include gyroviscosity
  integer :: jadv        ! 1 = use current density equation, not flux equation
  integer :: isources    ! 1 = include "source" terms in velocity advance
  integer :: istatic     ! 1 = do not advance velocity
  integer :: iestatic    ! 1 = do not advance fields
  integer :: igauge
  integer :: ivform      ! 0: V = v grad(phi).  1: V = R^2 v grad(phi)
  integer :: ibform      ! 0: multiply bz equation by r^2
  integer :: ihypeta     ! 1 = scale hyper-resistivity with eta
  integer :: ihypdx      ! scale hyper-resistivity with dx**ihypdx
  integer :: ikapscale   ! 1 = scale kappar with kappa
  integer :: iupwind     ! 1 = include upwinding term in pressure eqn

  ! numerical parameters
  integer :: ntimemax    ! number of timesteps
  integer :: nskip       ! number of timesteps per matrix recalculation
  integer :: iconstflux  ! 1 = conserve toroidal flux
  integer :: integrator  ! 0 = Crank-Nicholson, 1 = BDF2
  integer :: isplitstep  ! 1 = do timestep splitting
  integer :: imp_mod
  integer :: iteratephi  ! 1 = iterate field solve
  integer :: irecalc_eta ! 1 = recalculate transport coeffs after den solve
  integer :: iconst_eta  ! 1 = don't evolve resistivity
  real :: dt             ! timestep
  real :: thimp          ! implicitness parameter (for Crank-Nicholson)
  real :: thimp_ohm      ! implicitness parameter for ohmic heating
  real :: regular        ! regularization constant in chi equation
  real :: max_ke         ! max KE before fields are re-scaled when linear==1
  real :: eta_djdt       ! eta = (1+eta_djdt*dj/dt)*eta

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


  ! output parameters
  integer :: iprint     ! print extra debugging info
  integer :: itimer     ! print timing info
  integer :: ntimepr    ! number of timesteps per output  
  integer :: iglobalout ! 1 = write global restart files
  integer :: iglobalin  ! 1 = read global restart files

  ! general behavior
  integer :: iadapt     ! 1 = adapts mesh after initialization

  real :: dndt_fac      ! factor to multiply dn/dt term
  real :: dvdt_fac      ! factor to multiply dv/dt terms
  real :: dbdt_fac      ! factor to multiply db/dt terms

  
  integer :: istart
  real :: beta
  real :: pefac

  ! complex options
  integer :: ntor     ! toroidal mode number

!
!.....input quantities---defined in subroutine input or in namelist
!
  namelist / inputnl/                                          &
       itaylor,                                                &
       xzero,zzero,beta,rzero, libetap, xlim2, zlim2,     &
       numvar,idens,ipres,gyro,isources,nosig,itor,jadv,       &
       gam,db,gravr,gravz,                                     &
       p0,pi0,bzero,vzero,                                     &
       etar,eta0,amu,amuc,amupar,denm,                         &
       kappat,kappa0,kappar,kappax,kappah,                     &
       hyper,hyperi,hyperv,hyperc,hyperp,deex,                 &
       iper,jper,imask,amu_edge,com_bc,pedge,                  &
       eps,ln,                                                 &
       vloop,control_p,control_i,control_d,tcur,               &
       ipellet, pellet_x, pellet_z, pellet_rate, pellet_var,   &
       ionization, ionization_rate, ionization_temp, ionization_depth, &
       isink, sink1_x, sink2_x, sink1_z, sink2_z,              &
       sink1_rate, sink2_rate, sink1_var, sink2_var,           &
       ntimemax,dt,integrator,thimp,thimp_ohm,imp_mod,igauge,  &
       isplitstep,                                             &
       linear,nskip,eqsubtract,                                &
       itimer,iprint,ntimepr,iglobalout,iglobalin,             &
       irestart,istart,                                        &
       tcuro,djdpsi,xmag,zmag,xlim,zlim,                       &
       expn,q0,divertors,xdiv,zdiv,divcur,th_gs,p1,p2,         &
       idevice,igs,nv1equ,tol_gs,                              &
       iconstflux,regular,max_ke,                              &
       ntor,iadapt,istatic,iestatic,ivform, ibform,            &
       ihypeta,ikapscale,                                      &
       iteratephi,                                             &
       inonormalflow, inoslip_pol, inoslip_tor, inostress_tor, &
       iconst_t, inograd_t, inocurrent_pol, inocurrent_tor,    &
       iconst_bz,                                              &
       irecalc_eta,ihypdx, iconst_eta,                         &
       iupwind, dndt_fac, dvdt_fac, dbdt_fac,                  &
       eta_djdt,                                               &
       n_target, n_control_p, n_control_i, n_control_d


  !     derived quantities
  real :: pi,dbf,bdf,hypv,hypc,hypf,hypi,hypp,   &
       time,                                     &
       gbound,fbound

  ! PID controllers
  type(pid_control) :: i_control, n_control

  real :: eta_fac

  integer :: ntime, ntime0
  character*10 :: datec, timec

  integer :: ni(20),mi(20)  
  data mi /0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0/
  data ni /0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5/

! MPI variable(s)
  integer myrank, maxrank
end module basic

module t_data
  use p_data

  ! variables used to define the triangles
  real, allocatable :: atri(:),btri(:),ctri(:),ttri(:)
  vectype, allocatable :: gtri(:,:,:)
  real, allocatable :: rinv(:)
  integer, allocatable :: ist(:,:)

end module t_data

module arrays
  use p_data
  integer, parameter :: r8 = selected_real_kind(12,100)
  ! indices
  integer :: p,q,r,s
  integer :: maxdofs1, maxdofs2, maxdofsn
  integer, allocatable :: isvaln(:,:),isval1(:,:),isval2(:,:)
  real :: fint(-6:maxi,-6:maxi), xi(3),zi(3),df(0:4,0:4)

  ! arrays defined at all vertices
  ! any change to this list of variables needs to be taken into
  ! account in the arrayresizevec subroutine

  ! Arrays containing physical fields
  vectype, allocatable, target :: &
       field(:), field0(:), fieldi(:) 

  ! Arrays containing auxiliary variables
  vectype, allocatable :: &
       jphi(:), vor(:), com(:),                     &
       resistivity(:), tempvar(:),                  &
       kappa(:),sigma(:), sb1(:), sb2(:), sp1(:),   &
       visc(:), visc_c(:), bf(:), gyro_tau(:)
       


  ! Arrays for advance
  vectype, allocatable, target :: &
       phi(:), phiold(:),  &
       vel(:), velold(:),  &
       den(:), denold(:),  &
       pres(:), presold(:),  &
       q4(:), r4(:), qn4(:), qp4(:)  

  vectype, allocatable :: &
       veln(:), veloldn(:), phip(:),                          &
       b1_vel(:), b2_vel(:), b1_phi(:), b2_phi(:),            &
       tempcompare(:)


  ! the following pointers point to the vector containing the named field.
  ! set by assign_variables()
  vectype, pointer ::   u_v(:),   uo_v(:)
  vectype, pointer ::  vz_v(:),  vzo_v(:)
  vectype, pointer :: chi_v(:), chio_v(:)
  vectype, pointer :: psi_v(:), psio_v(:)
  vectype, pointer ::  bz_v(:),  bzo_v(:)
  vectype, pointer ::  pe_v(:),  peo_v(:)
  vectype, pointer :: den_v(:), deno_v(:)
  vectype, pointer ::   p_v(:),   po_v(:)


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

  ! the indicies of the named fields within their respective vectors
  integer :: u_i, vz_i, chi_i
  integer :: psi_i, bz_i, pe_i
  integer :: den_i, p_i
  integer :: bf_i

  ! the offset (relative to the node offset) of the named field within
  ! their respective vectors
  integer :: u_off, vz_off, chi_off
  integer :: psi_off, bz_off, pe_off
  integer :: den_off, p_off
  integer :: bf_off
  integer :: vecsize_vel, vecsize_phi, vecsize_n, vecsize_p
  
  ! the following pointers point to the locations of the named field within
  ! the respective vector.  set by assign_local_pointers()
  vectype, pointer ::   u1_l(:),   u0_l(:),   us_l(:)
  vectype, pointer ::  vz1_l(:),  vz0_l(:),  vzs_l(:)
  vectype, pointer :: chi1_l(:), chi0_l(:), chis_l(:) 
  vectype, pointer :: psi1_l(:), psi0_l(:), psis_l(:)
  vectype, pointer ::  bz1_l(:),  bz0_l(:),  bzs_l(:) 
  vectype, pointer ::  pe1_l(:),  pe0_l(:),  pes_l(:)
  vectype, pointer :: den1_l(:), den0_l(:), dens_l(:) 
  vectype, pointer ::   p1_l(:),   p0_l(:),   ps_l(:)

  contains


!==========================================================
! assign_variables
! ~~~~~~~~~~~~~~~
! Assigns variables to appropriate vectors for time advance
!==========================================================
    subroutine assign_variables()

      use basic

      implicit none

      if(isplitstep.eq.1) then

         u_v => vel
         uo_v => velold    
         psi_v => phi
         psio_v => phiold

         if(numvar.ge.2) then
            vz_v => vel
            vzo_v => velold
            bz_v => phi
            bzo_v => phiold
         endif

         if(numvar.ge.3) then
            chi_v => vel
            chio_v => velold
            pe_v => phi
            peo_v => phiold
            if(ipres.eq.1) then
               p_v => pres
               po_v => presold
            end if
         endif

         if(idens.eq.1) then
            den_v => den
            deno_v => denold
         end if

         u_i = 1
         psi_i = 1
         vz_i = 2
         bz_i = 2
         chi_i = 3
         pe_i = 3
         den_i = 1
         bf_i = 1

         if(ipres.eq.1) then
            p_i = 1
         else
            p_i = 3
         end if

      else
         u_v => phi
         uo_v => phiold
         psi_v => phi
         psio_v => phiold

         if(numvar.ge.2) then
            vz_v => phi
            vzo_v => phiold
            bz_v => phi
            bzo_v => phiold
         endif

         if(numvar.ge.3) then
            chi_v => phi
            chio_v => phiold
            pe_v => phi
            peo_v => phiold
            if(ipres.eq.1) then
               p_v => phi
               po_v => phiold
            endif
         endif
     
         if(idens.eq.1) then
            den_v => phi
            deno_v => phiold
         end if
         
         u_i = 1
         psi_i = 2
         vz_i = 3
         bz_i = 4
         chi_i = 5
         pe_i = 6    
         den_i = 2*numvar+1
         p_i = 2*numvar+idens+1
         bf_i = 1
         
      endif
      
      u_off = (u_i-1)*6
      psi_off = (psi_i-1)*6
      vz_off = (vz_i-1)*6
      bz_off = (bz_i-1)*6
      chi_off = (chi_i-1)*6
      pe_off = (pe_i-1)*6
      den_off = (den_i-1)*6
      p_off = (p_i-1)*6
      bf_off = (bf_i-1)*6

    end subroutine assign_variables

!======================================================
! assign_local_pointers
! ~~~~~~~~~~~~~~~~~~~~~
! Assigns local field pointers to appropriate locations
! in global field vectors.
!
!======================================================
    subroutine assign_local_pointers(inode)

      use basic

      implicit none

      integer, intent(in) :: inode
      integer :: ibegin, iendplusone, iend

      call entdofs(num_fields, inode, 0, ibegin, iendplusone)
      iend = ibegin+5

      psi1_l => field (ibegin+(psi_g-1)*6:iend+(psi_g-1)*6)
      psi0_l => field0(ibegin+(psi_g-1)*6:iend+(psi_g-1)*6)
      psis_l => fieldi(ibegin+(psi_g-1)*6:iend+(psi_g-1)*6)
        u1_l => field (ibegin+(  u_g-1)*6:iend+(  u_g-1)*6)
        u0_l => field0(ibegin+(  u_g-1)*6:iend+(  u_g-1)*6)
        us_l => fieldi(ibegin+(  u_g-1)*6:iend+(  u_g-1)*6)
       vz1_l => field (ibegin+( vz_g-1)*6:iend+( vz_g-1)*6)
       vz0_l => field0(ibegin+( vz_g-1)*6:iend+( vz_g-1)*6)
       vzs_l => fieldi(ibegin+( vz_g-1)*6:iend+( vz_g-1)*6)
       bz1_l => field (ibegin+( bz_g-1)*6:iend+( bz_g-1)*6)
       bz0_l => field0(ibegin+( bz_g-1)*6:iend+( bz_g-1)*6)
       bzs_l => fieldi(ibegin+( bz_g-1)*6:iend+( bz_g-1)*6)
      chi1_l => field (ibegin+(chi_g-1)*6:iend+(chi_g-1)*6)
      chi0_l => field0(ibegin+(chi_g-1)*6:iend+(chi_g-1)*6)
      chis_l => fieldi(ibegin+(chi_g-1)*6:iend+(chi_g-1)*6)
       pe1_l => field (ibegin+( pe_g-1)*6:iend+( pe_g-1)*6)
       pe0_l => field0(ibegin+( pe_g-1)*6:iend+( pe_g-1)*6)
       pes_l => fieldi(ibegin+( pe_g-1)*6:iend+( pe_g-1)*6)
        p1_l => field (ibegin+(  p_g-1)*6:iend+(  p_g-1)*6)
        p0_l => field0(ibegin+(  p_g-1)*6:iend+(  p_g-1)*6)
        ps_l => fieldi(ibegin+(  p_g-1)*6:iend+(  p_g-1)*6)
      den1_l => field (ibegin+(den_g-1)*6:iend+(den_g-1)*6)
      den0_l => field0(ibegin+(den_g-1)*6:iend+(den_g-1)*6)
      dens_l => fieldi(ibegin+(den_g-1)*6:iend+(den_g-1)*6)
    
    end subroutine assign_local_pointers
!================================
    subroutine createvec(vec, numberingid)
      implicit none
      integer :: numberingid, ndof
      vectype, allocatable :: vec(:)

      if(allocated(vec)) call deletevec(vec)

      call numdofs(numberingid, ndof)
      allocate(vec(ndof))

#ifdef USECOMPLEX
      call createppplvec(vec, numberingid, 1)
#else
      call createppplvec(vec, numberingid, 0)
#endif
      vec = 0.
    end subroutine createvec
!===============================
    subroutine createrealvec(vec, numberingid)
      implicit none
      integer :: numberingid, ndof
      real, allocatable ::  vec(:)
      
      if(allocated(vec)) call deleterealvec(vec)

      call numdofs(numberingid, ndof)
      allocate(vec(ndof))

      call createppplvec(vec, numberingid, 0)

      vec = 0.
    end subroutine createrealvec
!================================
    subroutine deletevec(vec)
      implicit none
      integer :: i
      vectype, allocatable :: vec(:) 
      
      call checkppplveccreated(vec, i)
      if(i .ne. 0) then
         call deleteppplvec(vec)
         deallocate(vec)
      else
         write(*,*) 'vector does not exist'
      endif
    end subroutine deletevec
!================================
    subroutine deleterealvec(vec)
      implicit none
      integer :: i
      real, allocatable :: vec(:) 
      
      call checkppplveccreated(vec, i)
      if(i .ne. 0) then
         call deleteppplvec(vec)
         deallocate(vec)
      else
         write(*,*) 'vector does not exist'
      endif
    end subroutine deleterealvec
!================================
end module arrays

  
module sparse
  integer, parameter :: numvar1_numbering = 1
  integer, parameter :: numvar2_numbering = 2
  integer, parameter :: numvar3_numbering = 3
  integer, parameter :: numvar4_numbering = 4
  integer, parameter :: numvar5_numbering = 5
  integer, parameter :: numvar6_numbering = 6
  integer, parameter :: numvar7_numbering = 7
  integer, parameter :: numvar8_numbering = 8

  integer, parameter :: s8matrix_sm = 1
  integer, parameter :: s7matrix_sm = 2
  integer, parameter :: s4matrix_sm = 3
  integer, parameter :: s5matrix_sm = 4
  integer, parameter :: s1matrix_sm = 5
  integer, parameter :: s2matrix_sm = 6
  integer, parameter :: d1matrix_sm = 7
  integer, parameter :: d2matrix_sm = 8
  integer, parameter :: d4matrix_sm = 9
  integer, parameter :: d8matrix_sm = 10
  integer, parameter :: q1matrix_sm = 11
  integer, parameter :: r2matrix_sm = 12
  integer, parameter :: r8matrix_sm = 13
  integer, parameter :: q2matrix_sm = 14
  integer, parameter :: q8matrix_sm = 15
  integer, parameter :: gsmatrix_sm = 16
  integer, parameter :: s9matrix_sm = 17
  integer, parameter :: d9matrix_sm = 18
  integer, parameter :: r9matrix_sm = 19
  integer, parameter :: q9matrix_sm = 20
  integer, parameter :: r14matrix_sm = 21
  integer, parameter :: mass_matrix_lhs = 22
  integer, parameter :: mass_matrix_lhs_dc = 23
  integer, parameter :: poisson_matrix_lhs = 24
  integer, parameter :: o1matrix_sm = 25
  integer, parameter :: o2matrix_sm = 26
  integer, parameter :: gs_matrix_rhs_dc = 27
  integer, parameter :: lp_matrix_rhs = 28
  integer, parameter :: lp_matrix_rhs_dc = 29
  integer, parameter :: bf_matrix_rhs_dc = 30
  integer, parameter :: gyro_torque_sm = 31
  integer, parameter :: q42matrix_sm = 32


  
end module sparse


subroutine insval(imatrix, val, icomplex, i, j, iop)

  implicit none
  
  integer, intent(in) :: imatrix, i, j, iop, icomplex
  vectype, intent(in) :: val

  if(val.eq.0.) return

  call insertval(imatrix, val, icomplex, i, j, iop)

end subroutine insval



subroutine resizevec(vec, ivecsize)
  use arrays
  implicit none
  integer ivecsize
  double precision :: vec

  call arrayresizevec(vec, ivecsize)

  return
end subroutine resizevec

subroutine arrayresizevec(vec, ivecsize)
  use arrays

  implicit none
  integer :: ivecsize, i
  double precision :: vec

  print *, "In arrayresizevec!", ivecsize

  call checkppplveccreated(vec, i)
  if(i .eq. 0) then
     call printfpointer(vec)
     write(*,*) 'trying to resize a vector that has not been created'
     call safestop(8844)
  endif


  call checksameppplvec(field, vec, i)
  if(i .eq. 1) then
     print *, "field"
     if(allocated(field)) deallocate(field, STAT=i)
     allocate(field(ivecsize))
     field = 0.
     call updateids(vec, field)
     print *, 'done'
     return
  endif

  call checksameppplvec(field0, vec, i)
  if(i .eq. 1) then
     print *, "field0"
     if(allocated(field0)) deallocate(field0, STAT=i)
     allocate(field0(ivecsize))
     field0 = 0.
     call updateids(vec, field0)
     return
  endif

  call checksameppplvec(fieldi, vec, i)
  if(i .eq. 1) then
     print *, "fieldi"
     if(allocated(fieldi)) deallocate(fieldi, STAT=i)
     allocate(fieldi(ivecsize))
     fieldi = 0.
     call updateids(vec, fieldi)
     return
  endif

  call checksameppplvec(jphi, vec, i)
  if(i .eq. 1) then
     print *, "jphi"
     if(allocated(jphi)) deallocate(jphi, STAT=i)
     allocate(jphi(ivecsize))
     jphi = 0.
     call updateids(vec, jphi)
     return
  endif
  
  call checksameppplvec(vor, vec, i)
  if(i .eq. 1) then
     print *, "vor"
     if(allocated(vor)) deallocate(vor, STAT=i)
     allocate(vor(ivecsize))
     vor = 0.
     call updateids(vec, vor)
     return
  endif
  
  call checksameppplvec(com, vec, i)
  if(i .eq. 1) then
     print *, "com"
     if(allocated(com)) deallocate(com, STAT=i)
     allocate(com(ivecsize))
     com = 0.
     call updateids(vec, com)
     return
  endif
  
  call checksameppplvec(resistivity, vec, i)
  if(i .eq. 1) then
     print *, "resistivity"
     if(allocated(resistivity)) deallocate(resistivity, STAT=i)
     allocate(resistivity(ivecsize))
     resistivity = 0.
     call updateids(vec, resistivity)
     return
  endif
  
  call checksameppplvec(tempvar, vec, i)
  if(i .eq. 1) then
     print *, "tempvar"
     if(allocated(tempvar)) deallocate(tempvar, STAT=i)
     allocate(tempvar(ivecsize))
     tempvar = 0.
     call updateids(vec, tempvar)
     return
  endif
  
  call checksameppplvec(kappa, vec, i)
  if(i .eq. 1) then
     print *, "kappa"
     if(allocated(kappa)) deallocate(kappa, STAT=i)
     allocate(kappa(ivecsize))
     kappa = 0.
     call updateids(vec, kappa)
     return
  endif

  call checksameppplvec(sigma, vec, i)
  if(i .eq. 1) then
     print *, "sigma"
     if(allocated(sigma)) deallocate(sigma, STAT=i)
     allocate(sigma(ivecsize))
     sigma = 0.
     call updateids(vec, sigma)
     return
  endif
      
  call checksameppplvec(sb1, vec, i)
  if(i .eq. 1) then
     print *, "sb1"
     if(allocated(sb1)) deallocate(sb1, STAT=i)
     allocate(sb1(ivecsize))
     sb1 = 0.
     call updateids(vec, sb1)
     return
  endif
  
  call checksameppplvec(sb2, vec, i)
  if(i .eq. 1) then
     print *, "sb2"
     if(allocated(sb2)) deallocate(sb2, STAT=i)
     allocate(sb2(ivecsize))
     sb2 = 0.
     call updateids(vec, sb2)
     return
  endif
  
  call checksameppplvec(sp1, vec, i)
  if(i .eq. 1) then
     print *, "sp1"
     if(allocated(sp1)) deallocate(sp1, STAT=i)
     allocate(sp1(ivecsize))
     sp1 = 0.
     call updateids(vec, sp1)
     return
  endif
  
  call checksameppplvec(visc, vec, i)
  if(i .eq. 1) then
     print *, "visc"
     if(allocated(visc)) deallocate(visc, STAT=i)
     allocate(visc(ivecsize))
     visc = 0.
     call updateids(vec, visc)
     return
  endif
  
  call checksameppplvec(visc_c, vec, i)
  if(i .eq. 1) then
     print *, "visc_c"
     if(allocated(visc_c)) deallocate(visc_c, STAT=i)
     allocate(visc_c(ivecsize))
     visc_c = 0.
     call updateids(vec, visc_c)
     return
  endif

  call checksameppplvec(bf, vec, i)
  if(i .eq. 1) then
     print *, "bf"
     if(allocated(bf)) deallocate(bf, STAT=i)
     allocate(bf(ivecsize))
     bf = 0.
     call updateids(vec, bf)
     return
  endif

  call checksameppplvec(gyro_tau, vec, i)
  if(i .eq. 1) then
     print *, "gyro_tau"
     if(allocated(gyro_tau)) deallocate(gyro_tau, STAT=i)
     allocate(gyro_tau(ivecsize))
     gyro_tau = 0.
     call updateids(vec, gyro_tau)
     return
  endif

  call checksameppplvec(phi, vec, i)
  if(i .eq. 1) then
     print *, "phi"
     if(allocated(phi)) deallocate(phi, STAT=i)
     allocate(phi(ivecsize))
     phi = 0.
     call updateids(vec, phi)
     return
  endif
     
  call checksameppplvec(phiold, vec, i)
  if(i .eq. 1) then
     print *, "phiold"
     if(allocated(phiold)) deallocate(phiold, STAT=i)
     allocate(phiold(ivecsize))
     phiold = 0.
     call updateids(vec, phiold)
     return
  endif

  call checksameppplvec(vel, vec, i)
  if(i .eq. 1) then
     print *, "vel"
     if(allocated(vel)) deallocate(vel, STAT=i)
     allocate(vel(ivecsize))
     vel = 0.
     call updateids(vec, vel)
     return
  endif
      
  call checksameppplvec(velold, vec, i)
  if(i .eq. 1) then
     print *, "velold"
     if(allocated(velold)) deallocate(velold, STAT=i)
     allocate(velold(ivecsize))
     velold = 0.
     call updateids(vec, velold)
     return
  endif
   
     
  call checksameppplvec(den, vec, i)
  if(i .eq. 1) then
     print *, "den"
     if(allocated(den)) deallocate(den, STAT=i)
     allocate(den(ivecsize))
     den = 0.
     call updateids(vec, den)
     return
  endif
  
  call checksameppplvec(denold, vec, i)
  if(i .eq. 1) then
     print *, "denold"
     if(allocated(denold)) deallocate(denold, STAT=i)
     allocate(denold(ivecsize))
     denold = 0.
     call updateids(vec, denold)
     return
  endif
  
  call checksameppplvec(pres, vec, i)
  if(i .eq. 1) then
     print *, "pres"
     if(allocated(pres)) deallocate(pres, STAT=i)
     allocate(pres(ivecsize))
     pres = 0.
     call updateids(vec, pres)
     return
  endif
  
  call checksameppplvec(presold, vec, i)
  if(i .eq. 1) then
     print *, "presold"
     if(allocated(presold)) deallocate(presold, STAT=i)
     allocate(presold(ivecsize))
     presold = 0.
     call updateids(vec, presold)
     return
  endif
     
  call checksameppplvec(q4, vec, i)
  if(i .eq. 1) then
     print *, "q4"
     if(allocated(q4)) deallocate(q4, STAT=i)
     allocate(q4(ivecsize))
     q4 = 0.
     call updateids(vec, q4)
     return
  endif

  call checksameppplvec(r4, vec, i)
  if(i .eq. 1) then
     print *, "r4"
     if(allocated(r4)) deallocate(r4, STAT=i)
     allocate(r4(ivecsize))
     r4 = 0.
     call updateids(vec, r4)
     return
  endif

  call checksameppplvec(qn4, vec, i)
  if(i .eq. 1) then
     print *, "qn4"
     if(allocated(qn4)) deallocate(qn4, STAT=i)
     allocate(qn4(ivecsize))
     qn4 = 0.
     call updateids(vec, qn4)
     return
  endif
  
  call checksameppplvec(qp4, vec, i)
  if(i .eq. 1) then
     print *, "qp4"
     if(allocated(qp4)) deallocate(qp4, STAT=i)
     allocate(qp4(ivecsize))
     qp4 = 0.
     call updateids(vec, qp4)
     return
  endif

  call checksameppplvec(veln, vec, i)
  if(i .eq. 1) then
     print *, "veln"
     if(allocated(veln)) deallocate(veln, STAT=i)
     allocate(veln(ivecsize))
     veln = 0.
     call updateids(vec, veln)
     return
  endif
      
  call checksameppplvec(veloldn, vec, i)
  if(i .eq. 1) then
     print *, "veloldn"
     if(allocated(veloldn)) deallocate(veloldn, STAT=i)
     allocate(veloldn(ivecsize))
     veloldn = 0.
     call updateids(vec, veloldn)
     return
  endif

  call checksameppplvec(phip, vec, i)
  if(i .eq. 1) then
     print *, "phip"
     if(allocated(phip)) deallocate(phip, STAT=i)
     allocate(phip(ivecsize))
     phip = 0.
     call updateids(vec, phip)
     return
  endif
  
  call checksameppplvec(b1_phi, vec, i)
  if(i .eq. 1) then
     print *, "b1_phi"
     if(allocated(b1_phi)) deallocate(b1_phi, STAT=i)
     allocate(b1_phi(ivecsize))
     b1_phi = 0.
     call updateids(vec, b1_phi)
     return
  endif
  
  call checksameppplvec(b2_phi, vec, i)
  if(i .eq. 1) then
     print *, "b2_phi"
     if(allocated(b2_phi)) deallocate(b2_phi, STAT=i)
     allocate(b2_phi(ivecsize))
     b2_phi = 0.
     call updateids(vec, b2_phi)
     return
  endif
  
  call checksameppplvec(b1_vel, vec, i)
  if(i .eq. 1) then
     print *, "b1_vel"
     if(allocated(b1_vel)) deallocate(b1_vel, STAT=i)
     allocate(b1_vel(ivecsize))
     b1_vel = 0.
     call updateids(vec, b1_vel)
     return
  endif
  
  call checksameppplvec(b2_vel, vec, i)
  if(i .eq. 1) then
     print *, "b2_vel"
     if(allocated(b2_vel)) deallocate(b2_vel, STAT=i)
     allocate(b2_vel(ivecsize))
     b2_vel = 0.
     call updateids(vec, b2_vel)
     return
  endif
    
  call checksameppplvec(tempcompare, vec, i)
  if(i .eq. 1) then
     print *, "tempcompare"
     if(allocated(tempcompare)) deallocate(tempcompare, STAT=i)
     allocate(tempcompare(ivecsize))
     tempcompare = 0.
     call updateids(vec, tempcompare)
     return
  endif  

  print *, "Error: unknown vector"

end subroutine arrayresizevec

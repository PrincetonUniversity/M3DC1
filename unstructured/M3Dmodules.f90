module p_data

  implicit none
  
  integer :: maxi     ! highest degree polynomial kept in nonlinear calculations
  integer :: ntimep   ! maximum number of timesteps
  integer :: ires     ! linear resolution of the plot files
  integer :: ntri     ! maximum number of HDF5 files
  integer :: maxplots ! maximum dimension of the graph array
  integer :: ijacobian

  parameter(maxplots=50, maxi=20, ntimep=10000, ires=201, ijacobian=1)
  
  integer, dimension(ires, ires) :: whichtri

end module p_data

module basic
  use p_data

#ifdef USECOMPLEX
  integer, parameter :: icomplex = 1
#else
  integer, parameter :: icomplex = 0
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

  ! boundary conditions
  integer :: iper, jper ! periodic boundary conditions
  integer :: imask      ! 1 = ignore 2-fluid terms near boundaries
  integer :: v_bc     ! bc on angular momentum.  
                      ! 0 = no-slip, 1 = no normal stress
  integer :: p_bc     ! bc on pressure.
                      !   0 = constant pressure, 1 = insulating
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

  

  ! numerical parameters
  integer :: linear      ! 1 = linear simulation; 0 = nonlinear simulation
  integer :: eqsubtract  ! 1 = subtract equilibrium in noninear simulations
  integer :: numvar      ! 1 = 2-field; 2 = reduced MHD; 3 = compressible MHD
  integer :: idens       ! evolve density
  integer :: ipres       ! evolve total and electron pressures separately
  integer :: gyro        ! include gyroviscosity
  integer :: jadv        ! 1 = use current density equation, not flux equation
  integer :: isources    ! 1 = include "source" terms in velocity advance
  integer :: ntimemax    ! number of timesteps
  integer :: nskip       ! number of timesteps per matrix recalculation
  integer :: iconstflux  ! 1 = conserve toroidal flux
  integer :: integrator  ! 0 = Crank-Nicholson, 1 = BDF2
  integer :: isplitstep  ! 1 = do timestep splitting
  integer :: imp_mod
  integer :: igauge
  real :: dt             ! timestep
  real :: thimp          ! implicitness parameter (for Crank-Nicholson)
  real :: thimp_ohm      ! implicitness parameter for ohmic heating
  real :: facw, facd
  real :: regular        ! regularization constant in chi equation

  ! current controller parameters
  real :: tcur        ! target toroidal current
  real :: control_p      ! proportionality constant
  real :: control_i      ! integral control inverse time-scale
  real :: control_d      ! derivative control time-scale

  ! output parameters
  integer :: iprint   ! print extra debugging info
  integer :: itimer   ! print timing info
  integer :: ntimepr  ! number of timesteps per output  

  ! domain parameters
  real :: xzero, zzero  ! cooridinates of lower left corner of domain
  
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
       xzero,zzero,beta,                                       &
       numvar,idens,ipres,gyro,isources,nosig,itor,jadv,       &
       gam,db,gravr,gravz,                                     &
       p0,pi0,bzero,vzero,                                     &
       etar,eta0,amu,amuc,amupar,denm,                         &
       kappat,kappa0,kappar,kappax,kappah,                     &
       hyper,hyperi,hyperv,hyperc,hyperp,deex,                 &
       iper,jper,imask,amu_edge,v_bc,p_db,com_bc,pedge,        &
       eps,ln,                                                 &
       vloop,control_p,control_i,control_d,tcur,               &
       ipellet, pellet_x, pellet_z, pellet_rate, pellet_var,   &
       ionization, ionization_rate, ionization_temp, ionization_depth, &
       ntimemax,dt,integrator,thimp,thimp_ohm,imp_mod,igauge,  &
       isplitstep,                                             &
       linear,nskip,eqsubtract,                                &
       itimer,iprint,ntimepr,                                  &
       irestart,istart,                                        &
       tcuro,djdpsi,xmag,zmag,xlim,zlim,facw,facd,             &
       expn,q0,divertors,xdiv,zdiv,divcur,th_gs,p1,p2,p_edge,  &
       idevice,igs,th_gs,                                      &
       iconstflux,regular,                                     &
       ntor

  !     derived quantities
  real :: pi,dbf,bdf,hypv,hypc,hypf,hypi,hypp,                         &
       time,timer,ajmax,errori,enormi,ratioi,                          &
       gbound,fbound
  integer ::  ni(20),mi(20),                                           &
       ntime,ntimer,nrank,ntimemin,ntensor,idebug, islutype

  character*8 :: filename(50)
  character*10 :: datec, timec
  
! initialization quantities
  integer :: ifirsts4_lu, ifirsts5_lu, ifirsts7_lu

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
       vtemp(:), resistivity(:), tempvar(:),        &
       kappa(:),sigma(:), sb1(:), sb2(:), sp1(:),   &
       visc(:), visc_c(:)


  ! Arrays for advance
  vectype, allocatable, target :: &
       phi(:), phiold(:),  &
       vel(:), velold(:),  &
       den(:), denold(:),  &
       pres(:), presold(:),  &
       q4(:), r4(:), qn4(:), qp4(:)  

  vectype, allocatable :: &
       veln(:), veloldn(:), phip(:),                              &
       b1vector(:), b2vector(:), b3vector(:), b4vector(:),        &
       b5vector(:), tempcompare(:)

!!$  vectype, allocatable, target :: &
!!$       vel(:), velold(:), vel0(:),  vels(:),    &
!!$       phi(:), phiold(:), phi0(:),  phis(:),    &
!!$       den(:), denold(:), den0(:),  dens(:),    &
!!$       pres(:), presold(:), pres0(:), press(:), &
!!$       q4(:), r4(:), qn4(:), qp4(:)    



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

  ! the offset (relative to the node offset) of the named field within
  ! their respective vectors
  integer :: u_off, vz_off, chi_off
  integer :: psi_off, bz_off, pe_off
  integer :: den_off, p_off
  integer :: vecsize, vecsize1
  
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
         p_i = 1

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
         p_i = 2*numvar+2
      endif
      
      u_off = (u_i-1)*6
      psi_off = (psi_i-1)*6
      vz_off = (vz_i-1)*6
      bz_off = (bz_i-1)*6
      chi_off = (chi_i-1)*6
      pe_off = (pe_i-1)*6
      den_off = (den_i-1)*6
      p_off = (p_i-1)*6

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
      integer :: numberingid, i, ndof
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
      integer :: numberingid, i, ndof
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
  integer, parameter :: s6matrix_sm = 1
  integer, parameter :: s8matrix_sm = 2
  integer, parameter :: s7matrix_sm = 3
  integer, parameter :: s4matrix_sm = 4
  integer, parameter :: s3matrix_sm = 5
  integer, parameter :: s5matrix_sm = 6
  integer, parameter :: s1matrix_sm = 7
  integer, parameter :: s2matrix_sm = 8
  integer, parameter :: d1matrix_sm = 9
  integer, parameter :: d2matrix_sm = 10
  integer, parameter :: d4matrix_sm = 11
  integer, parameter :: d8matrix_sm = 12
  integer, parameter :: q1matrix_sm = 13
  integer, parameter :: r2matrix_sm = 14
  integer, parameter :: r8matrix_sm = 15
  integer, parameter :: q2matrix_sm = 16
  integer, parameter :: q8matrix_sm = 17
  integer, parameter :: gsmatrix_sm = 18
  integer, parameter :: s9matrix_sm = 19
  integer, parameter :: d9matrix_sm = 20
  integer, parameter :: r9matrix_sm = 21
  integer, parameter :: q9matrix_sm = 22
  integer, parameter :: r14matrix_sm = 23
  
end module sparse


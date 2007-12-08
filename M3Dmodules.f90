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
  integer :: iresolve    ! 1 = do second velocity solve after field solve
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
       isplitstep,iresolve,                                    &
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
  real, allocatable :: atri(:),btri(:),ctri(:),ttri(:),gtri(:,:,:)
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
  real :: xsep(5), zsep(5), graphit(0:ntimep,maxplots)

  ! arrays defined at all vertices
  ! any change to this list of variables needs to be taken into
  ! account in the arrayresizevec subroutine
  vectype, allocatable, target :: &
       vel(:), velold(:),  &
       phi(:), phiold(:),  &
       den(:), denold(:),  &
       pres(:), presold(:), &
       q4(:), r4(:), qn4(:), qp4(:)
  real, allocatable, target :: &
       vel0(:),  vels(:),  &
       phi0(:),  phis(:),  &
       den0(:),  dens(:),  &
       pres0(:), press(:)

  vectype, allocatable :: &
       veln(:), veloldn(:),                                       &
       phip(:), phioldn(:),                                       &
       b1vector(:), b2vector(:), b3vector(:), b4vector(:),        &
       jphi(:),vor(:),com(:),                                     &
       presoldn(:),                                               &
       b5vector(:), vtemp(:), resistivity(:), tempvar(:),         &
       kappa(:),sigma(:), sb1(:), sb2(:), sp1(:),                 &
       visc(:), visc_c(:), tempcompare(:)


  ! the following pointers point to the vector containing the named field.
  ! set by assign_variables()
  vectype, pointer :: phi1_v(:), phio_v(:)
  vectype, pointer ::  vz1_v(:),  vzo_v(:)
  vectype, pointer :: chi1_v(:), chio_v(:)
  vectype, pointer :: psi1_v(:), psio_v(:)
  vectype, pointer ::  bz1_v(:),  bzo_v(:)
  vectype, pointer ::  pe1_v(:),  peo_v(:)
  vectype, pointer :: den1_v(:), deno_v(:)
  vectype, pointer ::   p1_v(:),   po_v(:)
  real, pointer :: phi0_v(:), phis_v(:)
  real, pointer ::  vz0_v(:),  vzs_v(:)
  real, pointer :: chi0_v(:), chis_v(:)
  real, pointer :: psi0_v(:), psis_v(:)
  real, pointer ::  bz0_v(:),  bzs_v(:)
  real, pointer ::  pe0_v(:),  pes_v(:)
  real, pointer :: den0_v(:), dens_v(:)
  real, pointer ::   p0_v(:),   ps_v(:)


  ! the indicies of the named fields within their respective vectors
  integer :: phi_i, vz_i, chi_i
  integer :: psi_i, bz_i, pe_i
  integer :: den_i, p_i

  ! the offset (relative to the node offset) of the named field within
  ! their respective vectors
  integer :: phi_off, vz_off, chi_off
  integer :: psi_off, bz_off, pe_off
  integer :: den_off, p_off
  integer :: vecsize, vecsize1
  
  ! the following pointers point to the locations of the named field within
  ! the respective vector.  set by assign_vectors()
  vectype, pointer :: phi1_l(:)
  vectype, pointer ::  vz1_l(:)
  vectype, pointer :: chi1_l(:)
  vectype, pointer :: psi1_l(:)
  vectype, pointer ::  bz1_l(:)
  vectype, pointer ::  pe1_l(:)
  vectype, pointer :: den1_l(:)
  vectype, pointer ::   p1_l(:)
  real, pointer :: phi0_l(:)
  real, pointer ::  vz0_l(:)
  real, pointer :: chi0_l(:)
  real, pointer :: psi0_l(:)
  real, pointer ::  bz0_l(:)
  real, pointer ::  pe0_l(:)
  real, pointer :: den0_l(:)
  real, pointer ::   p0_l(:)

  contains
!================================
    subroutine assign_variables()

      use basic

      implicit none

      if(isplitstep.eq.1) then

         phi0_v => vel0
         phi1_v => vel
         phis_v => vels
         phio_v => velold
         
         psi0_v => phi0
         psi1_v => phi
         psis_v => phis
         psio_v => phiold

         if(numvar.ge.2) then
            vz0_v => vel0
            vz1_v => vel
            vzs_v => vels
            vzo_v => velold

            bz0_v => phi0
            bz1_v => phi
            bzs_v => phis
            bzo_v => phiold
         endif

         if(numvar.ge.3) then
            chi0_v => vel0
            chi1_v => vel
            chis_v => vels
            chio_v => velold

            pe0_v => phi0
            pe1_v => phi
            pes_v => phis
            peo_v => phiold
         endif

         if(idens.eq.1) then
            den0_v => den0
            den1_v => den
            dens_v => dens
            deno_v => denold
         endif

         if(ipres.eq.1) then
            p0_v => pres0
            p1_v => pres
            ps_v => press
            po_v => presold
         endif

         phi_i = 1
         psi_i = 1
         vz_i = 2
         bz_i = 2
         chi_i = 3
         pe_i = 3
         den_i = 1
         p_i = 1

      else
         phi0_v => phi0
         phi1_v => phi
         phis_v => phis
         phio_v => phiold
 
         psi0_v => phi0
         psi1_v => phi
         psis_v => phis
         psio_v => phiold

         if(numvar.ge.2) then
            vz0_v => phi0
            vz1_v => phi
            vzs_v => phis
            vzo_v => phiold

            bz0_v => phi0
            bz1_v => phi
            bzs_v => phis
            bzo_v => phiold

         endif
         if(numvar.ge.3) then
            chi0_v => phi0
            chi1_v => phi
            chis_v => phis
            chio_v => phiold

            pe0_v => phi0
            pe1_v => phi
            pes_v => phis
            peo_v => phiold
         endif
     
         if(idens.eq.1) then
            den0_v => phi0
            den1_v => phi
            dens_v => phis
            deno_v => phiold
         endif

         if(ipres.eq.1) then
            p0_v => phi0
            p1_v => phi
            ps_v => phis
            po_v => phiold
         endif
         
         phi_i = 1
         psi_i = 2
         vz_i = 3
         bz_i = 4
         chi_i = 5
         pe_i = 6    
         den_i = 2*numvar+1
         p_i = 2*numvar+2
      endif
      
      phi_off = (phi_i-1)*6
      psi_off = (psi_i-1)*6
      vz_off = (vz_i-1)*6
      bz_off = (bz_i-1)*6
      chi_off = (chi_i-1)*6
      pe_off = (pe_i-1)*6
      den_off = (den_i-1)*6
      p_off = (p_i-1)*6

    end subroutine assign_variables

    subroutine assign_vectors(inode)

      use basic

      implicit none

      integer, intent(in) :: inode
      integer :: ibegin, iendplusone

      call entdofs(vecsize, inode, 0, ibegin, iendplusone)
      
      phi0_l => phi0_v(ibegin+phi_off:ibegin+phi_off+5)
      phi1_l => phi1_v(ibegin+phi_off:ibegin+phi_off+5)
      psi0_l => psi0_v(ibegin+psi_off:ibegin+psi_off+5)
      psi1_l => psi1_v(ibegin+psi_off:ibegin+psi_off+5)

      if(numvar.ge.2) then
         vz0_l => vz0_v(ibegin+vz_off:ibegin+vz_off+5)
         vz1_l => vz1_v(ibegin+vz_off:ibegin+vz_off+5)
         bz0_l => bz0_v(ibegin+bz_off:ibegin+bz_off+5)
         bz1_l => bz1_v(ibegin+bz_off:ibegin+bz_off+5)
      endif
      
      if(numvar.ge.3) then
         chi0_l => chi0_v(ibegin+chi_off:ibegin+chi_off+5)
         chi1_l => chi1_v(ibegin+chi_off:ibegin+chi_off+5)
         pe0_l => pe0_v(ibegin+pe_off:ibegin+pe_off+5)
         pe1_l => pe1_v(ibegin+pe_off:ibegin+pe_off+5)
      endif

      if(isplitstep.eq.1) call entdofs(1, inode, 0, ibegin, iendplusone)

      if(idens.eq.1) then
         den0_l => den0_v(ibegin+den_off:ibegin+den_off+5)
         den1_l => den1_v(ibegin+den_off:ibegin+den_off+5)
      endif
      
      if(ipres.eq.1) then
         p0_l => p0_v(ibegin+p_off:ibegin+p_off+5)
         p1_l => p1_v(ibegin+p_off:ibegin+p_off+5)
      endif
      
    end subroutine assign_vectors
!================================
    subroutine createvec(vec, numberingid)
      implicit none
      integer numberingid, i, ndof
      vectype, allocatable ::  vec(:)
      
!      write(*,*) 'numbering id is ',numberingid
      ! if vec has not been created vec == 0
      call checkveccreated(vec, i)
      if(i .eq. 0) then
         call numdofs(numberingid, ndof)
         allocate(vec(ndof))
         call createppplvec(vec, numberingid)
         vec = 0.
      else
         write(*,*) 'vector is already created'
      endif
    end subroutine createvec
!===============================
    subroutine createrealvec(vec, numberingid)
      implicit none
      integer numberingid, i, ndof
      real, allocatable ::  vec(:)
      
!      write(*,*) 'numbering id is ',numberingid
      ! if vec has not been created vec == 0
      call checkveccreated(vec, i)
      if(i .eq. 0) then
         call numdofs(numberingid, ndof)
         allocate(vec(ndof))
         call createppplvec(vec, numberingid)
         vec = 0.
      else
         write(*,*) 'vector is already created'
      endif
    end subroutine createrealvec
!================================
    subroutine deletevec(vec)
      implicit none
      integer :: i
      vectype, allocatable :: vec(:) 
      
      call checkveccreated(vec, i)
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
      
      call checkveccreated(vec, i)
      if(i .ne. 0) then
         call deleteppplvec(vec)
         deallocate(vec)
      else
         write(*,*) 'vector does not exist'
      endif
    end subroutine deleterealvec
!================================
    subroutine arrayresizevec(vec, ivecsize)
      implicit none
      integer numberingid, ivecsize, i
      double precision :: vec

      call checkveccreated(vec, i)
      if(i .eq. 0) then
         call printfpointer(vec)
         write(*,*) 'trying to resize a vector that has not been created'
         call safestop(8844)
      endif

      call checksamevec(vel, vec, i)
      if(i .eq. 1) then
         if(allocated(vel)) deallocate(vel, STAT=i)
         allocate(vel(ivecsize))
         vel = 0.
         call updateids(vec, vel)
         return
      endif

      call checksamevec(vels, vec, i)
      if(i .eq. 1) then
         if(allocated(vels)) deallocate(vels, STAT=i)
         allocate(vels(ivecsize))
         vels = 0.
         call updateids(vec, vels)
         return
      endif
      
      call checksamevec(veln, vec, i)
      if(i .eq. 1) then
         if(allocated(veln)) deallocate(veln, STAT=i)
         allocate(veln(ivecsize))
         veln = 0.
         call updateids(vec, veln)
         return
      endif
      
      call checksamevec(veloldn, vec, i)
      if(i .eq. 1) then
         if(allocated(veloldn)) deallocate(veloldn, STAT=i)
         allocate(veloldn(ivecsize))
         veloldn = 0.
         call updateids(vec, veloldn)
         return
      endif
      
      call checksamevec(velold, vec, i)
      if(i .eq. 1) then
         if(allocated(velold)) deallocate(velold, STAT=i)
         allocate(velold(ivecsize))
         velold = 0.
         call updateids(vec, velold)
         return
      endif
      
      call checksamevec(vel0, vec, i)
      if(i .eq. 1) then
         if(allocated(vel0)) deallocate(vel0, STAT=i)
         allocate(vel0(ivecsize))
         vel0 = 0.
         call updateids(vec, vel0)
         return
      endif
            
      call checksamevec(phi, vec, i)
      if(i .eq. 1) then
         if(allocated(phi)) deallocate(phi, STAT=i)
         allocate(phi(ivecsize))
         phi = 0.
         call updateids(vec, phi)
         return
      endif

      call checksamevec(phis, vec, i)
      if(i .eq. 1) then
         if(allocated(phis)) deallocate(phis, STAT=i)
         allocate(phis(ivecsize))
         phis = 0.
         call updateids(vec, phis)
         return
      endif

      call checksamevec(phip, vec, i)
      if(i .eq. 1) then
         if(allocated(phip)) deallocate(phip, STAT=i)
         allocate(phip(ivecsize))
         phip = 0.
         call updateids(vec, phip)
         return
      endif

      call checksamevec(phiold, vec, i)
      if(i .eq. 1) then
         if(allocated(phiold)) deallocate(phiold, STAT=i)
         allocate(phiold(ivecsize))
         phiold = 0.
         call updateids(vec, phiold)
         return
      endif

      call checksamevec(phi0, vec, i)
      if(i .eq. 1) then
         if(allocated(phi0)) deallocate(phi0, STAT=i)
         allocate(phi0(ivecsize))
         phi0 = 0.
         call updateids(vec, phi0)
         return
      endif

      call checksamevec(jphi, vec, i)
      if(i .eq. 1) then
         if(allocated(jphi)) deallocate(jphi, STAT=i)
         allocate(jphi(ivecsize))
         jphi = 0.
         call updateids(vec, jphi)
         return
      endif

      call checksamevec(vor, vec, i)
      if(i .eq. 1) then
         if(allocated(vor)) deallocate(vor, STAT=i)
         allocate(vor(ivecsize))
         vor = 0.
         call updateids(vec, vor)
         return
      endif

      call checksamevec(com, vec, i)
      if(i .eq. 1) then
         if(allocated(com)) deallocate(com, STAT=i)
         allocate(com(ivecsize))
         com = 0.
         call updateids(vec, com)
         return
      endif

      call checksamevec(den, vec, i)
      if(i .eq. 1) then
         if(allocated(den)) deallocate(den, STAT=i)
         allocate(den(ivecsize))
         den = 0.
         call updateids(vec, den)
         return
      endif

      call checksamevec(den0, vec, i)
      if(i .eq. 1) then
         if(allocated(den0)) deallocate(den0, STAT=i)
         allocate(den0(ivecsize))
         den0 = 0.
         call updateids(vec, den0)
         return
      endif

      call checksamevec(dens, vec, i)
      if(i .eq. 1) then
         if(allocated(dens)) deallocate(dens, STAT=i)
         allocate(dens(ivecsize))
         dens = 0.
         call updateids(vec, dens)
         return
      endif

      call checksamevec(denold, vec, i)
      if(i .eq. 1) then
         if(allocated(denold)) deallocate(denold, STAT=i)
         allocate(denold(ivecsize))
         denold = 0.
         call updateids(vec, denold)
         return
      endif

      call checksamevec(pres, vec, i)
      if(i .eq. 1) then
         if(allocated(pres)) deallocate(pres, STAT=i)
         allocate(pres(ivecsize))
         pres = 0.
         call updateids(vec, pres)
         return
      endif

      call checksamevec(pres0, vec, i)
      if(i .eq. 1) then
         if(allocated(pres0)) deallocate(pres0, STAT=i)
         allocate(pres0(ivecsize))
         pres0 = 0.
         call updateids(vec, pres0)
         return
      endif

      call checksamevec(press, vec, i)
      if(i .eq. 1) then
         if(allocated(press)) deallocate(press, STAT=i)
         allocate(press(ivecsize))
         press = 0.
         call updateids(vec, press)
         return
      endif

      call checksamevec(presold, vec, i)
      if(i .eq. 1) then
         if(allocated(presold)) deallocate(presold, STAT=i)
         allocate(presold(ivecsize))
         presold = 0.
         call updateids(vec, presold)
         return
      endif

      call checksamevec(r4, vec, i)
      if(i .eq. 1) then
         if(allocated(r4)) deallocate(r4, STAT=i)
         allocate(r4(ivecsize))
         r4 = 0.
         call updateids(vec, r4)
         return
      endif

      call checksamevec(q4, vec, i)
      if(i .eq. 1) then
         if(allocated(q4)) deallocate(q4, STAT=i)
         allocate(q4(ivecsize))
         q4 = 0.
         call updateids(vec, q4)
         return
      endif

      call checksamevec(qn4, vec, i)
      if(i .eq. 1) then
         if(allocated(qn4)) deallocate(qn4, STAT=i)
         allocate(qn4(ivecsize))
         qn4 = 0.
         call updateids(vec, qn4)
         return
      endif

      call checksamevec(qp4, vec, i)
      if(i .eq. 1) then
         if(allocated(qp4)) deallocate(qp4, STAT=i)
         allocate(qp4(ivecsize))
         qp4 = 0.
         call updateids(vec, qp4)
         return
      endif

      call checksamevec(b1vector, vec, i)
      if(i .eq. 1) then
         if(allocated(b1vector)) deallocate(b1vector, STAT=i)
         allocate(b1vector(ivecsize))
         b1vector = 0.
         call updateids(vec, b1vector)
         return
      endif

      call checksamevec(b2vector, vec, i)
      if(i .eq. 1) then
         if(allocated(b2vector)) deallocate(b2vector, STAT=i)
         allocate(b2vector(ivecsize))
         b2vector = 0.
         call updateids(vec, b2vector)
         return
      endif

      call checksamevec(b3vector, vec, i)
      if(i .eq. 1) then
         if(allocated(b3vector)) deallocate(b3vector, STAT=i)
         allocate(b3vector(ivecsize))
         b3vector = 0.
         call updateids(vec, b3vector)
         return
      endif

      call checksamevec(b4vector, vec, i)
      if(i .eq. 1) then
         if(allocated(b4vector)) deallocate(b4vector, STAT=i)
         allocate(b4vector(ivecsize))
         b4vector = 0.
         call updateids(vec, b4vector)
         return
      endif

      call checksamevec(b5vector, vec, i)
      if(i .eq. 1) then
         if(allocated(b5vector)) deallocate(b5vector, STAT=i)
         allocate(b5vector(ivecsize))
         b5vector = 0.
         call updateids(vec, b5vector)
         return
      endif

      call checksamevec(vtemp, vec, i)
      if(i .eq. 1) then
         if(allocated(vtemp)) deallocate(vtemp, STAT=i)
         allocate(vtemp(ivecsize))
         vtemp = 0.
         call updateids(vec, vtemp)
         return
      endif

      call checksamevec(resistivity, vec, i)
      if(i .eq. 1) then
         if(allocated(resistivity)) deallocate(resistivity, STAT=i)
         allocate(resistivity(ivecsize))
         resistivity = 0.
         call updateids(vec, resistivity)
         return
      endif

      call checksamevec(tempvar, vec, i)
      if(i .eq. 1) then
         if(allocated(tempvar)) deallocate(tempvar, STAT=i)
         allocate(tempvar(ivecsize))
         tempvar = 0.
         call updateids(vec, tempvar)
         return
      endif

      call checksamevec(kappa, vec, i)
      if(i .eq. 1) then
         if(allocated(kappa)) deallocate(kappa, STAT=i)
         allocate(kappa(ivecsize))
         tempvar = 0.
         call updateids(vec, kappa)
         return
      endif

      call checksamevec(tempcompare, vec, i)
      if(i .eq. 1) then
         if(allocated(tempcompare)) deallocate(tempcompare, STAT=i)
         allocate(tempcompare(ivecsize))
         tempvar = 0.
         call updateids(vec, tempcompare)
         return
      endif

      call checksamevec(sigma, vec, i)
      if(i .eq. 1) then
         if(allocated(sigma)) deallocate(sigma, STAT=i)
         allocate(sigma(ivecsize))
         tempvar = 0.
         call updateids(vec, sigma)
         return
      endif
    end subroutine arrayresizevec
end module arrays

!below is so we can call arrayresizevec from c++
subroutine resizevec(vec, ivecsize)
  use arrays
  implicit none
  integer ivecsize
  double precision :: vec

  call arrayresizevec(vec, ivecsize)

  return
end subroutine resizevec
  

  
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
  integer, parameter :: s10matrix_sm = 24
  integer, parameter :: d10matrix_sm = 25
  integer, parameter :: q10matrix_sm = 26
  integer, parameter :: r10matrix_sm = 27
  
end module sparse


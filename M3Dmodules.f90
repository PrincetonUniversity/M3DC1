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
  real :: etar, eta0  ! resistivity = etar + eta0/T^(3/2)
  real :: kappa       ! pressure diffusion
  real :: kappat      ! isotropic temperature conductivity
  real :: kappar      ! anisotropic (field-aligned) temperature conductivity
  real :: denm        ! artificial density diffusion
  real :: deex        ! scale length of hyperviscosity term
  real :: hyper,hyperi,hyperv,hyperc,hyperp

  ! physical parameters
  integer :: itor     ! 1 = cylindrical coordinates; 0 = cartesian coordinates
  real :: db          ! ion skin depth
  real :: gam         ! ratio of specific heats
  real :: grav        ! gravitational acceleration
  real :: vloop       ! loop voltage

  ! general equilibrium parameters
  integer :: irestart ! 1 = reads restart file as initial condition
  integer :: itaylor  ! equilibrium
  real :: bzero       ! guide field
  real :: p0, pi0     ! total, ion pressures
  real :: ln          ! length of equilibrium gradient
  real :: eps         ! size of initial perturbation

  ! toroidal equilibrium parameters
  integer :: divertors! number of divertors
  real :: xmag, zmag  ! position of magnetic axis
  real :: xlim, zlim  ! position of limiter
  real :: xdiv, zdiv  ! position of divertor
  real :: tcuro       ! toroidal current
  real :: divcur      ! current in divertor (as fraction of tcuro)
  real :: djdpsi
  real :: p1, p2, pedge
  real :: expn        ! density = pressure**expn
  

  ! numerical parameters
  integer :: linear      ! 1 = linear simulation; 0 = nonlinear simulation
  integer :: eqsubtract  ! 1 = subtract equilibrium in noninear simulations
  integer :: numvar      ! 1 = 2-field; 2 = reduced MHD; 3 = compressible MHD
  integer :: idens       ! evolve density
  integer :: ipres       ! evolve total and electron pressures separately
  integer :: gyro        ! include gyroviscosity
  integer :: imask       ! 1 = ignore 2-fluid terms near boundaries
  integer :: ntimemax    ! number of timesteps
  integer :: nskip       ! number of timesteps per matrix recalculation
  integer :: iconstflux  ! 1 = conserve toroidal flux
  integer :: isources    ! 1 = include source terms in velocity advance
  integer :: integrator  ! 0 = Crank-Nicholson, 1 = BDF2
  real :: dt             ! timestep
  real :: thimp          ! implicitness parameter (for Crank-Nicholson)
  real :: facw, facd
  real :: regular        ! regularization constant in chi equation

  ! output parameters
  integer :: iprint   ! print extra debugging info
  integer :: itimer   ! print timing info
  integer :: ntimepr  ! number of timesteps per output  

  ! domain parameters
  integer :: iper, jper ! periodic boundary conditions
  real :: xzero, zzero  ! cooridinates of lower left corner of domain
  
  integer :: maxs,itest,isecondorder,istart
  real :: beta
  real :: pefac
!
!.....input quantities---defined in subroutine input or in namelist
!
  namelist / inputnl/                                          &
       linear,maxs,ntimemax,ntimepr,itor,                      &
       irestart,itaylor,itest,isecondorder,imask,nskip,        &
       numvar,istart,idens,ipres,thimp,amu,etar,dt,p1,p2,p0,   &
       tcuro,djdpsi,xmag,zmag,xlim,zlim,facw,facd,db,cb,       &
       bzero,hyper,hyperi,hyperv,hyperc,hyperp,gam,eps,        &
       kappa,iper,jper,iprint,itimer,xzero,zzero,beta,pi0,     &
       eqsubtract,denm,grav,kappat,kappar,ln,amuc,iconstflux,  &
       regular,deex,gyro,vloop,eta0,isources,pedge,integrator, &
       expn,divertors,xdiv,zdiv,divcur

  !     derived quantities
  real :: tt,pi,                                                       &
       time,timer,ajmax,errori,enormi,ratioi,                          &
       gbound,fbound
  integer ::  ni(20),mi(20),                                           &
       ntime,ntimer,nrank,ntimemin,ntensor,idebug, islutype
  real :: ekin, emag, ekind, emagd, ekino, emago, ekindo, emagdo,      &
       ekint,emagt,ekintd,emagtd,ekinto,emagto,ekintdo,emagtdo,        &
       ekinp,emagp,ekinpd,emagpd,ekinpo,emagpo,ekinpdo,emagpdo,        &
       ekinph,ekinth,emagph,emagth,ekinpho,ekintho,emagpho,emagtho,    &
       ekin3,ekin3d,ekin3h,emag3,ekin3o,ekin3do,ekin3ho,emag3o,        &
       emag3h,emag3d,emag3ho,emag3do,chierror,tflux0,totcur0,          &
       ekappar,ekappat,eeta
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
  integer :: maxdofs1, maxdofs2, maxdofs3
  integer, allocatable :: isvaln(:,:),isval1(:,:),isval2(:,:)
  real :: fint(-6:maxi,-6:maxi), xi(3),zi(3),df(0:4,0:4)
  real :: xsep(5), zsep(5), graphit(0:ntimep,maxplots)

  ! arrays defined at all vertices
  ! any change to this list of variables needs to be taken into
  ! account in the arrayresizevec subroutine
  real, allocatable::                                             &
       vel(:), vels(:), veln(:), veloldn(:),                      &
       velold(:), vel0(:), vel1(:),                               &
       phi(:), phis(:), phip(:),                                  &
       phiold(:), phi0(:), phi1(:),                               &
       jphi(:),sb1(:),sb2(:),sp1(:),                              &
       vor(:),com(:),                                             &
       den(:),den0(:),denold(:),deni(:),dens(:),                  &
       pres(:),pres0(:),presold(:),press(:),                      &
       r4(:),q4(:),qn4(:),qp4(:),                                 &
       b1vector(:), b2vector(:), b3vector(:), b4vector(:),        &
       b5vector(:), vtemp(:), resistivity(:), tempvar(:)

  contains
!================================
    subroutine createvec(vec, numberingid)
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
    end subroutine createvec
!================================
    subroutine deletevec(vec)
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
    end subroutine deletevec
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
      
      call checksamevec(vel1, vec, i)
      if(i .eq. 1) then
         if(allocated(vel1)) deallocate(vel1, STAT=i)
         allocate(vel1(ivecsize))
         vel1 = 0.
         call updateids(vec, vel1)
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

      call checksamevec(phi1, vec, i)
      if(i .eq. 1) then
         if(allocated(phi1)) deallocate(phi1, STAT=i)
         allocate(phi1(ivecsize))
         phi1 = 0.
         call updateids(vec, phi1)
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

      call checksamevec(sb1, vec, i)
      if(i .eq. 1) then
         if(allocated(sb1)) deallocate(sb1, STAT=i)
         allocate(sb1(ivecsize))
         sb1 = 0.
         call updateids(vec, sb1)
         return
      endif

      call checksamevec(sb2, vec, i)
      if(i .eq. 1) then
         if(allocated(sb2)) deallocate(sb2, STAT=i)
         allocate(sb2(ivecsize))
         sb2 = 0.
         call updateids(vec, sb2)
         return
      endif

      call checksamevec(sp1, vec, i)
      if(i .eq. 1) then
         if(allocated(sp1)) deallocate(sp1, STAT=i)
         allocate(sp1(ivecsize))
         sp1 = 0.
         call updateids(vec, sp1)
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

      call checksamevec(deni, vec, i)
      if(i .eq. 1) then
         if(allocated(deni)) deallocate(deni, STAT=i)
         allocate(deni(ivecsize))
         deni = 0.
         call updateids(vec, deni)
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
  integer, parameter :: r1matrix_sm = 13
  integer, parameter :: r2matrix_sm = 14
  integer, parameter :: r8matrix_sm = 15
  integer, parameter :: q2matrix_sm = 16
  integer, parameter :: q8matrix_sm = 17
  integer, parameter :: gsmatrix_sm = 18
  integer, parameter :: s9matrix_sm = 19
  integer, parameter :: d9matrix_sm = 20
  integer, parameter :: r9matrix_sm = 21
  integer, parameter :: q9matrix_sm = 22
  
end module sparse


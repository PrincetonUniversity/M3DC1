module p_data

  implicit none
  
  integer :: maxi     ! highest degree polynomial kept in nonlinear calculations
  integer :: ntimep   ! maximum number of timesteps
  integer :: ires     ! linear resolution of the plot files
  integer :: ntri     ! maximum number of HDF5 files
  integer :: ntridim  ! 
  integer :: maxhdf5
  integer :: maxplots ! maximum dimension of the graph array

  parameter(maxplots=50, maxi=20, ntimep=1000, ires=201, maxhdf5=30)
  
  integer, dimension(ires, ires) :: whichtri

end module p_data

module basic
  use p_data

  ! transport coefficients
  real :: amu         ! incompressible viscosity
  real :: amuc        ! compressible viscosity
  real :: etar        ! resistivity
  real :: kappa       ! pressure diffusion
  real :: kappat      ! isotropic temperature conductivity
  real :: kappar      ! anisotropic (field-aligned) temperature conductivity
  real :: denm        ! artificial density diffusion
  real :: hyper,hyperi,hyperv,hyperc,hyperp

  ! physical parameters
  integer :: itor     ! 1 = cylindrical coordinates; 0 = cartesian coordinates
  real :: db          ! ion skin depth
  real :: cb
  real :: gam         ! ratio of specific heats
  real :: grav        ! gravitational acceleration

  ! general equilibrium parameters
  integer :: irestart ! 1 = reads restart file as initial condition
  integer :: itaylor  ! equilibrium
  real :: bzero       ! guide field
  real :: p0, pi0     ! total, ion pressures
  real :: ln          ! length of equilibrium gradient
  real :: eps         ! size of initial perturbation

  ! toroidal equilibrium parameters
  real :: xmag, zmag  ! position of magnetic axis
  real :: xlim, zlim  ! position of limiter
  real :: tcuro       ! toroidal current
  real :: djdpsi
  real :: p1, p2

  ! numerical parameters
  integer :: linear      ! 1 = linear simulation; 0 = nonlinear simulation
  integer :: eqsubtract  ! 1 = subtract equilibrium in noninear simulations
  integer :: numvar      ! 1 = 2-field; 2 = reduced MHD; 3 = compressible MHD
  integer :: idens       ! evolve density
  integer :: ipres       ! evolve total and electron pressures separately
  integer :: imask       ! 1 = ignore 2-fluid terms near boundaries
  integer :: ntimemax    ! number of timesteps
  integer :: nskip       ! number of timesteps per matrix recalculation
  integer :: iconstflux  ! 1 = conserve toroidal flux
  real :: dt             ! timestep
  real :: thimp          ! implicitness parameter
  real :: facw, facd

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
  namelist / inputnl/                                        &
       linear,maxs,ntimemax,ntimepr,itor,                    &
       irestart,itaylor,itest,isecondorder,imask,nskip,      &
       numvar,istart,idens,ipres,thimp,amu,etar,dt,p1,p2,p0, &
       tcuro,djdpsi,xmag,zmag,xlim,zlim,facw,facd,db,cb,     &
       bzero,hyper,hyperi,hyperv,hyperc,hyperp,gam,eps,      &
       kappa,iper,jper,iprint,itimer,xzero,zzero,beta,pi0,   &
       eqsubtract,denm,grav,kappat,kappar,ln,amuc,iconstflux

  !     derived quantities
  real :: tt,gamma4,gamma2,gamma3,dpsii,psimin,psilim,pi,              &
       time,timer,ajmax,errori,enormi,ratioi,                          &
       tmesh, tsetup, tfirst,tsolve,tsecond,tzero,tthird,gbound
  integer ::  ni(20),mi(20) ,nbcgs,nbcp,nbcv,nbcn,iboundmax,           &
       ntime,ntimer,nrank,ntimemin,ntensor,idebug, islutype
  real :: ekin, emag, ekind, emagd, ekino, emago, ekindo, emagdo,      &
       ekint,emagt,ekintd,emagtd,ekinto,emagto,ekintdo,emagtdo,        &
       ekinp,emagp,ekinpd,emagpd,ekinpo,emagpo,ekinpdo,emagpdo,        &
       ekinph,ekinth,emagph,emagth,ekinpho,ekintho,emagpho,emagtho,    &
       ekin3,ekin3d,ekin3h,emag3,ekin3o,ekin3do,ekin3ho,emag3o,        &
       emag3h,emag3d,emag3ho,emag3do,tflux,chierror,totcur, area
  real :: xary(0:ires-1),yary(0:ires-1),mif(0:maxhdf5-1),maf(0:maxhdf5-1)
  character*8 :: filename(50)
  character*10 :: datec, timec
  
! initialization quantities
  integer ifirstd1_lu, ifirsts1_lu, ifirstd2_lu, ifirsts2_lu,    &
       ifirstr1_lu, ifirstr2_lu, ifirstq2_lu, ifirsts3_lu,       &
       ifirsts4_lu, ifirsts5_lu, ifirsts6_lu, ifirsts7_lu,       &
       ifirsts8_lu, ifirstd8_lu, ifirstq8_lu, ifirstr8_lu

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
  integer, allocatable :: iboundgs(:), iboundv(:), iboundp(:), iboundn(:)
  integer, allocatable :: isvaln(:,:),isval1(:,:),isval2(:,:)
  real :: fint(-6:maxi,-6:maxi), xi(3),zi(3),df(0:4,0:4)
  real :: xsep(5), zsep(5), graphit(0:ntimep,maxplots)
  real, allocatable :: psibounds(:), velbounds(:), combounds(:)

  ! arrays defined at all vertices
  real(r8), allocatable::                                         &
       b1vecini(:), vel(:), vels(:), veln(:),                     &
       velold(:), vel0(:), vel1(:),                               &
       b2vecini(:), phi(:), phis(:),                              &
       phiold(:), phi0(:), phi1(:),                               &
       jphi(:),jphi0(:),sb1(:),sb2(:),sb3(:),sp1(:),              &
       vor(:),vor0(:),com(:),com0(:),                             &
       den(:),den0(:),denold(:),deni(:),                          &
       pres(:),pres0(:),r4(:),q4(:),qn4(:),                       &
       b1vector(:), b2vector(:), b3vector(:), b4vector(:),        &
       b5vector(:), vtemp(:),                                     &
       fun1(:),fun4(:),fun2(:),fun3(:)

  contains
!================================
    subroutine createvec(vec, numberingid)
      implicit none
      integer numberingid, i, ndof
      real(r8), allocatable ::  vec(:)
      
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
      integer numberingid, i
      real(r8), allocatable :: vec(:) 
      
      call checkveccreated(vec, i)
      if(i .ne. 0) then
         call deleteppplvec(vec)
         deallocate(vec)
      else
         write(*,*) 'vector does not exist'
      endif
    end subroutine deletevec
end module arrays

  
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
  
end module sparse
!!$

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
  real :: amu         ! viscosity
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
  integer :: linear   ! 1 = linear simulation; 0 = nonlinear simulation
  integer :: eqsubtract  ! 1 = subtract equilibrium in noninear simulations
  integer :: numvar
  integer :: idens    ! evolve density
  integer :: ipres    ! evolve total and electron pressures separately
  integer :: imask    ! 1 = ignore 2-fluid terms near boundaries
  integer :: ntimemax ! number of timesteps
  integer :: nskip    ! number of timesteps per matrix recalculation
  real :: dt          ! timestep
  real :: thimp       ! implicitness parameter
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
       eqsubtract,denm,grav,kappat,kappar,ln

  !     derived quantities
  real :: tt,gamma4,gamma2,gamma3,dpsii,psimin,psilim,pi,              &
       time,timer,ajmax,errori,enormi,ratioi,                          &
       tmesh, tsetup, tfirst,tsolve,tsecond,tzero,tthird,gbound
  integer ::  ni(20),mi(20) ,nbcgs,nbcp,nbcv,nbcn,iboundmax,           &
       ntime,ntimer,nrank,ntimemin,ntensor, iframe,                    &
       ihdf5, idebug, islutype
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
  data iframe /0/

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
       jphi(:),jphi0(:),sb1(:),sb2(:),sb3(:),vor(:),vor0(:),      &
       com(:),com0(:),den(:),den0(:),denold(:),                   &
       pres(:),pres0(:),r4(:),q4(:),qn4(:),                       &
       b1vector(:), b2vector(:), b3vector(:), b4vector(:),        &
       b5vector(:), vtemp(:),                                     &
       fun1(:),fun4(:),fun2(:),fun3(:)

end module arrays

module sparse
  use sparse_matrix
  type(sparse_params_obj), pointer :: spo_numvar1, spo_numvar2, spo_numvar3
  type(sparse_matrix_obj) :: s6matrix_sm, s8matrix_sm, s7matrix_sm, &
       s4matrix_sm, s3matrix_sm, s5matrix_sm, s1matrix_sm,    &
       s2matrix_sm, d1matrix_sm, d2matrix_sm, d4matrix_sm,    &
       d8matrix_sm, r1matrix_sm, r2matrix_sm, r8matrix_sm,    &   
       q2matrix_sm, q8matrix_sm, gsmatrix_sm
  
end module sparse

#ifdef IS_LIBRARY
subroutine reducedquintic(isfirst, inmyrank, inmaxrank)
#else
Program Reducedquintic
#endif

!   Ref:  [1] Strang and Fix, An Analysis of the Finite Element Method, page 83
!         [2] G.R. Cowper, et al, AIAA Journal, Vol 7, no. 10 page 19
!         [3] JCP Vol 147 p 318-336 (1998)

  use p_data
  use t_data
  use basic
  use arrays
  use newvar_mod
  use sparse
  use sparse_params
#ifdef mpi
  use supralu_dist_mod
  implicit none
  include 'mpif.h'
#else
  implicit none
#endif

#ifdef IS_LIBRARY
  integer, intent(in) :: isfirst, inmyrank, inmaxrank
#endif
  integer :: j, i, ier, numvars, ifail, maxts, numelms, numnodes
  integer :: index1, index2
  real :: dens, dtmin, ratemin, ratemax, temp(5)

  real :: tstart, tend

  double precision :: coords(3)
  integer :: nodeids(4)

  integer :: ibegin, iendplusone

#ifdef mpi
#ifndef IS_LIBRARY
  ! Start up message passing, SuperLU process grid
  call MPI_Init(ier)
  if (ier /= 0) then
     print *,'Error initializing MPI:',ier
     call safestop(1)
  endif
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_rank:',ier
     call safestop(1)
  endif
  ! initialize the SUPERLU process grid
  call SLUD_init
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_size:',ier
     call safestop(1)
  endif
#endif
#ifdef IS_LIBRARY
  write(*,*) 'input is myrank, maxrank, isfirst', inmyrank, &
       inmaxrank, isfirst
      
  if(isfirst .eq. 1) then
     call SLUD_init
  endif
  myrank = inmyrank
  maxrank = inmaxrank
#endif
#else
  myrank = 0
#endif
  temp = 0.
#ifndef IS_LIBRARY
  print *, 'starting library'
  if(myrank.eq.0 .and. itimer.ge.1) call second(tstart)
  call loadmesh("struct.dmg", "struct-dmg.sms")
  if(myrank.eq.0 .and. itimer.ge.1) then
     call second(tend)
     print *, 'Time spent in loadmesh: ', tend-tstart
  endif
#else
  print *, 'starting program'
#endif
  call numnod(numnodes)
  call numfac(numelms)
  write(*,*) 'numnodes and numfaces',numnodes,numelms

  ! arrays for triangle parameters
  allocate(atri(numelms),btri(numelms),ctri(numelms),ttri(numelms),      &
       gtri(20,18,numelms)) 
  
  call precalc_whattri()

  if(myrank.eq.0) then
     call date_and_time( datec, timec)
     write(*,1001) datec(1:4),datec(5:6),datec(7:8),                        &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("M3D-C1 VERSION 1.0    DATE: "a4,1x,a2,1x,a2,3x,               &
          "TIME: "a2,":",a2,":",a4,/)
  endif

  pi = acos(-1.)
  
  ! special switch installed 12/03/04 to write slu matrices
  idebug = 0

  ! Choose which test problem to run (0 indicates reading input file C1input)
  itest=0

  ! initialize needed variables and define geometry and triangles
  if(myrank.eq.0 .and. itimer.ge.1) call second(tstart)
  call init
  if(myrank.eq.0 .and. itimer.ge.1) then
     call second(tend)
     print *, 'Time spent in init: ', tend-tstart
  endif

  if(maxrank .gt. 1) then
     write(*,*) 'Currently analysis can only be done with 1 proc'
     call safestop(42)
  endif

  ! sparse_params_objects
  allocate(spo_numvar1, spo_numvar2, spo_numvar3)
  call init_spo(spo_numvar1,1)
  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'finished with init_spo'
  if(numvar .ge. 2) call init_spo(spo_numvar2,2)
  if(numvar .ge. 3) call init_spo(spo_numvar3,3)

  ! calculate the RHS (forcing function)
  call rhsdef
  
  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'finished rhsdef'

  ! note that the matrix indices now refer to vertices, not triangles.
  ! ...............................................................

  ! initialize the solution to an equilibrium and save
  if(irestart.eq.1) then
     call rdrestart
     if(istart.ne.0) then
        ntimer = 0
        timer = 0.
     endif
  else
     ntimer = 0
     timer = 0.
     if(idens.eq.1) then
        call denequ(den0, 1)
        if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'finished denequ'
        if(myrank.eq.0) call oneplot(den0,1,1,"den0",1)
        if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'finished oneplot'
     endif

     if(itor.eq.1 .and. itaylor.eq.1) then
        numvars = numvar

        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        call gradshafranov
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           write(*,*) "Time spent in gradshafranov: ", tend - tstart
        endif

        numvar = numvars
!        phiold = phi
        do i=1,numnodes
           do j=1,6
              index2 = (i-1)*6*numvar + j
              index1 = (i-1)*6 + j
              phiold(index2) = phi(index1)
           enddo
        enddo
     else
        velold = 0.
        phiold = 0.
        call velequ(velold, numvar)
        call phiequ(phiold, numvar)
     endif

     vel0 = velold
     phi0 = phiold
     denold = den0

     ! calculate initial perturbed fields
     call velinit(vel)
     call phiinit(phi)
     if(idens.eq.1) call deninit(den)

     phis = 0
     vels = 0
     dens = 0
  endif                     !  end of the branch on restart/no restart

  ! correct for left-handed coordinates
  do i=1,numnodes
     call entdofs(numvar, i, 0, ibegin, iendplusone)
     if(numvar.eq.1) then
        j = 5
     else
        j = 11
     endif
     phi(ibegin:ibegin+j) = -phi(ibegin:ibegin+j)
     vel(ibegin:ibegin+j) = -vel(ibegin:ibegin+j)
     phi0(ibegin:ibegin+j) = -phi0(ibegin:ibegin+j)
     vel0(ibegin:ibegin+j) = -vel0(ibegin:ibegin+j)
     phiold(ibegin:ibegin+j) = -phiold(ibegin:ibegin+j)
     velold(ibegin:ibegin+j) = -velold(ibegin:ibegin+j)
  enddo

  if(myrank.eq.0) call plotit(vel,phi,0)

  ! calculate the equilibrium current density
  if(itaylor.ne.4) then
     call newvar(phi0,jphi0,numnodes,numvar,1,VAR_J,1)
  else
     jphi0 = 0
  endif

  ! calculate the equilibrium vorticity and velocity divergence
  vor0 = 0.
  com0 = 0.


  ! combine the equilibrium and perturbed fields of linear=0
  ! unless eqsubtract = 1
  if(linear.eq.0 .and. eqsubtract.eq.0) then
     phi = phi + phi0
     phi0 = 0.

     vel = vel + vel0
     vel0 = 0.

     if(idens.eq.1) then
        den = den + den0
        den0 = 0.
     endif
  endif

  call plotit(vel+vel0,phi+phi0,1)

!  call axis(phi,xsep,zsep,0)

! start the time dependent loop


  ifail=0
  ier = 0

  time = timer
  errori = 0.
  enormi = 0.
  ratioi = 0.
  dtmin = 0.001*dt
  ntime = ntimer
  call energy
!      if (myrank.eq.0) call output
  if(ntimemax.le.ntimer) go to 101

  do ntime=ntimer+1,ntimemax
     time = time + dt

     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call onestep
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        write(*,*) "Time spent in onestep: ", tend - tstart
     endif
!     call exportfield2(1,numvar,phi, ntime)
     
     if(myrank.eq.0) then
        if(itimer.eq.1) call second(tstart)
        call output
        if(itimer.eq.1) then
           call second(tend)
           write(*,*) "Time spent in output: ", tend - tstart
        endif
     endif

!     if(linear.eq.1) call scaleback

     if(ekin .gt. 100.) then
        write(*,*) 'ekin is greater than 100'
        go to 100
     endif
  enddo ! ntime

 100  continue
  call second(tsolve)
  
  maxts = ntime-1

101 continue

  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'about to export field'
!      call exportfield2(1,numvar,phi, 0) ! 0 for now for adaptive-loop.sh script
!      call exportfield2(1,numvar,jphi, 1)
  ratemin = 0.
  ratemax = 0.
  do ntime=ntimemin,maxts
     ratemin = min(ratemin,graphit(ntime,26))
     ratemax = max(ratemax,graphit(ntime,26))
  enddo

!  call errorcalc(numvar, phi, 1)
!
  if (myrank.eq.0) then
     write(*,5002) tsolve-tfirst,                                      &
          numvar,amu,etar,dt,thimp,cb,db,hyper,hyperi,hyperv,hyperc,      &
          ratemin,ratemax,ajmax,graphit(ntimemax,25)
     write(*,5003) linear, itaylor, 0, imask, irestart,                &
          facd,bzero,eps

     write(9,5002) tsolve-tfirst,                                      &
          numvar,amu,etar,dt,thimp,cb,db,hyper,hyperi,hyperv,hyperc,      &
          ratemin,ratemax,ajmax,graphit(ntimemax,25)
     write(9,5003) linear, itaylor, 0, imask, irestart,                &
          facd,bzero,eps
  endif

  if(ntime.gt.1 .and. myrank.eq.0) then
     call plotenergy(graphit,ntimep,maxts,ntimemin,numvar)
     open(99,file='C1graphit',form='formatted',status='unknown')
     write(99,8001) maxts
     do i=1,maxts
        write(99,8002) (graphit(i,j),j=1,30)
     enddo
8001 format(i5)
8002 format(1p10e12.4)
  endif
999 continue
  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'writing the restart file'
  call wrrestart
  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'done writing the restart file'
      
!     free memory from sparse matrices
  call free_smo(gsmatrix_sm)
  call free_smo(s6matrix_sm)
  call free_smo(s8matrix_sm)
  call free_smo(s7matrix_sm)
  call free_smo(s4matrix_sm)
  call free_smo(s3matrix_sm)
  call free_smo(s5matrix_sm)
  call free_smo(s1matrix_sm)
  call free_smo(s2matrix_sm)
  call free_smo(d1matrix_sm)
  call free_smo(d2matrix_sm)
  call free_smo(d4matrix_sm)
  call free_smo(d8matrix_sm)
  call free_smo(r1matrix_sm)
  call free_smo(r2matrix_sm)
  call free_smo(r8matrix_sm)
  call free_smo(q2matrix_sm)
  call free_smo(q8matrix_sm)
  call free_spo(spo_numvar1)
  if(numvar .ge. 2) call free_spo(spo_numvar2)
  if(numvar .ge. 3) call free_spo(spo_numvar3)
  deallocate(spo_numvar1, spo_numvar2, spo_numvar3)
  call deletesearchstructure()
  if (myrank.eq.0) call plote
  
5002 format(" tsolve =", 1pe11.4,   "  numvar =", 0p1i4,               &
          "   amu,etar =", 1p2e10.2,  /,"  dt,thimp =",1p2e10.2,          &
          "  cb,db=",1p2e12.2,/,                                          &
          "  hyper,hyperi,hyperv,hyperc = ",  1p4e12.4,  /,               &
          "  ratemin,ratemax,ajmax,reconflux = ",1p4e12.4) 
5003 format(" linear, itaylor, isetup, imask, irestart: ",             &
          5i4, / ," facd, bzero, eps ", 1p3e12.4)
2323 format(1p5e12.4)
2121 format(" start of time dependent loop")
5001 format("ifail .ne. 0 after call to f04aaf")
7011 format(1p6e20.12)
7012 format(6e20.12)

#ifdef IS_LIBRARY
  return
end subroutine reducedquintic

#else
  write(*,*) myrank
  call safestop(2)

end Program Reducedquintic
#endif



!============================================================
! ONESTEP
! ~~~~~~~
! advance solution to time ntime+1
!============================================================
subroutine onestep

  use p_data
  use t_data
  use basic
  use arrays
  use sparse
  use newvar_mod

  implicit none

  integer :: l, jer, i, jone, j, itri
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  real :: eerror, ediff, etotd, sum, fintl(-6:maxi,-6:maxi), d2term(18)
  integer, allocatable:: itemp(:)
  integer :: ndofs, numelms, numnodes

  real :: tstart, tend
  
  ! acbauer -- these variables are not initialized before use
  ediff = 1.
  etotd = 1.

  call numnod(numnodes)
  call numfac(numelms)

  ! define the inverse density array deni
!!$  if(idens.eq.1) then
!!$     if(linear.eq.1 .or. eqsubtract.eq.1) then
!!$        call inverse(den+den0,deni)
!!$     else
!!$        call inverse(den,deni)
!!$     endif
!!$  endif

  ! define the current vector jphi and the RHS vectors sb1 and sb2
  !      call newvar(phi,jphi,numnodes,numvar,1,VAR_J,1)
  call newvar(phi,sb1,numnodes,numvar,1,VAR_SB1,1)
!!$  if(numvar.ge.2) call newvar(phi,sb2,numnodes,numvar,1,VAR_SB2,1)
  sb2 = 0.
  if(numvar.ge.3) call newvar(phi,sb3,numnodes,numvar,1,VAR_SB3,1)

!  call conserve_tflux()

  if(ntime.le.ntimer+1.or. (linear.eq.0 .and. mod(ntime,nskip).eq.0)) then
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call ludefall
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        write(*,*) " onestep: Time spent in ludefall:", tend - tstart
     endif
  endif

  veln = vel
  
  ! Advance Velocity
  ! ================

  ! Calculate LU decomposition of velocity matrix if needed

  ! b1vector = r1matrix_sm * phi(n)
  if(iprint.ge.1) write(*,*) "before sparseR8_A_dot_X"
  call matrix_vector_mult(r1matrix_sm, phi, b1vector, r1matrix_sm%spo_ptr%numcols)
  
  ! vtemp = d1matrix_sm * vel(n)
  vtemp = 0.
  if(iprint.eq.1)write(*,*) "before second sparseR8_A_dot_X"
  call matrix_vector_mult(d1matrix_sm,vel,vtemp, d1matrix_sm%spo_ptr%numcols)

  vtemp = vtemp + dt*facd*b1vecini + b1vector + r4

  ! apply boundary conditions
  do l=1,nbcv
     vtemp(iboundv(l)) = velbounds(l)
  enddo

  ! solve linear system with rhs in vtemp (note LU-decomp done first time)
  if(myrank.eq.0 .and. iprint.eq.1) write(*,*) "before dsupralu_solve_s1handle"
  call solve(s1matrix_sm, vtemp, jer)
  if(myrank.eq.0 .and. iprint.eq.1) write(*,*) "after dsupralu_solve_s1handle"
  if(jer.ne.0) then
     write(*,*) 'after sparseR8d_solve', jer
     call safestop(42)
  endif
! ok to here     call printarray(vtemp, 150, 0, 'vtemp aa')

!
!.....coding to calculate the error in the delsquared chi equation
  chierror = 0
  sum = 0
  if(numvar.ge.3) then
     call newvar(vtemp,com,numnodes,numvar,3,VAR_COM,0)
     
     do itri=1,numelms
        call calcfint(fintl, maxi, atri(itri),btri(itri), ctri(itri))
        call calcd2term(itri, d2term, fintl)
        do j=1,18
           jone = isval1(itri,j)
           chierror = chierror + d2term(j)*com(jone)
        enddo
     enddo                  ! loop over itri
     if(myrank.eq.0 .and. iprint.ge.1) then
        print *, "Error in com = ", chierror 
     endif
     
     if(hyperc.gt.0) then
        call smoother3(com,vtemp,numnodes,numvar,3)
        call newvar(vtemp,com,numnodes,numvar,3,VAR_COM,0)
     endif
     call oneplot(com,1,1,"com",0)
     
  endif

  call newvar(vtemp,vor,numnodes,numvar,1,VAR_VOR,1)
  
  if(hyperc.gt.0) then
     !
     !.......calculate vorticity, apply smoothing operator, and redefine vor array
     call smoother1(vor,vtemp,numnodes,numvar,1)
     call newvar(vtemp,vor,numnodes,numvar,1,VAR_VOR,1)
     call oneplot(com,1,1,"vor",0)
  endif

!.....new velocity solution at time n+1 (or n* for second order advance)
!!$      call numdofs(numvar, ndofs)
!!$      do i=1,ndofs
!!$         veln(i) = vel(i)
!!$         vel(i) = vtemp(i)
!!$      enddo
  vel = vtemp

!
! Advance Density
! ===============
  if(idens.eq.1) then
     if(iprint.ge.1) write(*,*) "s8handle"

     ! b2vector = r8matrix_lu * vel(n+1)
     call matrix_vector_mult(r8matrix_sm,vel,b2vector, r8matrix_sm%spo_ptr%numcols)

     ! b3vector = q8matrix_sm * vel(n)
     call matrix_vector_mult(q8matrix_sm,veln,b3vector, q8matrix_sm%spo_ptr%numcols)

     ! vtemp = d8matrix_sm * phi(n)
     vtemp = 0.
     call matrix_vector_mult(d8matrix_sm,den,vtemp, d8matrix_sm%spo_ptr%numcols)
     allocate(itemp(numnodes*6)) ! this is used to make sure that we don't double count the sum for periodic dofs
     do l=1,numnodes*6
        itemp(l) = 1
     enddo
     do l=1,numnodes
        call entdofs(1, l, 0, ibegin, iendplusone)
        call entdofs(numvar, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           vtemp(ibegin+i) = vtemp(ibegin+i) + itemp(ibegin+i) * &
                (b2vector(ibeginnv+i) + b3vector(ibeginnv+i) + qn4(ibegin+i))
           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)

     do l=1,nbcn
        vtemp(iboundn(l)) = 0.
        if(linear.eq.0 .and. eqsubtract.eq.0) then
           vtemp(iboundn(l)) = vtemp(iboundn(l)) + denold(iboundn(l))
        endif
     enddo

     ! solve linear system...LU decomposition done first time
! -- okay to here      call printarray(vtemp, 150, 0, 'vtemp on')
     call solve(s8matrix_sm, vtemp, jer)
     if(jer.ne.0) then
        write(*,*) 'after 2nd sparseR8d_solve', jer
        call safestop(29)
     endif

     ! new field solution at time n+1 (or n* for second order advance)
     den = vtemp
     !      call printarray(den, 150, 0, 'den vtem')
  endif

!
! Advance Pressure
! ================
!
! Advance Fields
! ==============
  
  ! Calculate LU decomposition of field matrix if needed

#ifdef mpi

  ! b2vector = r2matrix_lu * vel(n+1)
  call matrix_vector_mult(r2matrix_sm,vel,b2vector, r2matrix_sm%spo_ptr%numcols)

  ! b3vector = q2matrix_lu * vel(n)
  call matrix_vector_mult(q2matrix_sm,veln,b3vector, q2matrix_sm%spo_ptr%numcols)

  ! vtemp = d2matrix_sm * phi(n)
  vtemp = 0.
  call matrix_vector_mult(d2matrix_sm,phi,vtemp, d2matrix_sm%spo_ptr%numcols)

  vtemp = vtemp + dt*facd*b2vecini + b2vector + b3vector + q4
  
  ! Insert boundary conditions
  do l=1,nbcp
     vtemp(iboundp(l)) = psibounds(l)
     if(linear.eq.0 .and. eqsubtract.eq.0) then
        vtemp(iboundp(l)) = vtemp(iboundp(l)) + phiold(iboundp(l))
     endif
  enddo

  ! solve linear system...LU decomposition done first time

  call solve(s2matrix_sm, vtemp, jer)
  if(jer.ne.0) then
     write(*,*) 'after 2nd sparseR8d_solve', jer
     call safestop(29)
  endif
#endif

! new field solution at time n+1 (or n* for second order advance)
!      phiold = phi
  phi = vtemp

  call energy

  eerror = 0.
  if(ntime.gt.5) then
     eerror = 2.*abs(ediff-etotd)/(abs(ediff)+abs(etotd))
  endif
!      if(eerror .gt. 0.10 .and. linear.eq.0) dt = 0.9*dt
!      if(dt.lt.dtmin) then
!        write(*,*) 'timestep too small, dt,dtmin =  ' ,dt, dtmin
!        go to 100
!      endif
!      if(eerror .gt. 1.00  .and. linear.eq.0) then
!        write(*,*) 'eerror too large' , eerror
!        go to 100
!      endif
!
!     NOTE:  if facw=1., facd is zeroed after the first cycle
  if(facd .ne. 0 .and. facw.eq.1 .and. time.ge.0.1) then
     facd = 0.
     ntimemin = max(ntimemin,ntime+1)
  endif

end subroutine onestep



subroutine conserve_tflux()

  use basic
  use t_data
  use arrays

  implicit none

  integer :: numelms, numnodes, index, itri, j, jone, j2, ivertex
  real :: correction
  real :: d2term(18), fintl(-6:maxi,-6:maxi)

  call numnod(numnodes)
  call numfac(numelms)

  totcur = 0
  area = 0.
  tflux = 0.
  ! calculate the total perturbed current and area, and toroidal flux
  do itri=1,numelms
     call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
     call calcd2term(itri, d2term, fintl)
     do j=1,18
        jone = isval1(itri,j)
        totcur = totcur + d2term(j)*jphi(jone)

        if(numvar.ge.2) then
           j2 = isvaln(itri,j) + 6
           tflux = tflux + d2term(j)*phi(j2)
        endif

     enddo
     do j=1,13,6
        area = area + d2term(j)
     enddo
  enddo                     ! loop over itri

  ! adjust toroidal field and boundary condition to conserve toroidal flux
  if(numvar.ge.2) then
     correction = tflux/area
     gbound = gbound - correction
     do ivertex=1,numnodes
        index = 6*numvar*(ivertex-1) + 7
        phi(index) = phi(index) - correction
     enddo
  endif
  
  if(myrank.eq.0 .and. iprint.eq.1) write(*,*) "tflux, gbound",tflux,gbound
  if(myrank.eq.0) write(9,*) "tflux, gbound",tflux,gbound

end subroutine conserve_tflux

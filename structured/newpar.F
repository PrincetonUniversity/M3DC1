!123456789-123456789-123456789-123456789-123456789-123456789-123456789-23
!
!     VERSION 12.1   10/01/06
!     runs with three velocity variables (numvar=3)
!     added terms for electron pressure in ludefphi
!       write HDF5 files
!       runs and restarts under mpi with disk read/writes
!       compiles and runs on sgi,NERSC,petrel
!       accepts namelist input
!       smoother for vorticity and compressibility
!       parallel ludefvel and ludefphi
!       density advance equation added, inverse density defined
!       runs with idens=1 with constant density..same results for idens=0
!       1/ne multiplying hyperresistivity terms.
!       jper periodicity (in addition to iper)
!       itaylor=4 wave propagation tests
!       imask coordinates centered in triangle
!       sign error corrected in ludefden
!       Grad-Shafranov operator corrected
!       Ohmic and Viscous heating terms added for numvar=3
!       Order of loops in ludef changed to help disk reads
!       Viscous terms corrected for symmetric stress tensor form
!       eqsubtract term added
!       kappat term added (including g13erm)
!       section timings added to check parallel scaling
!       performed energy check with Samtaney code
!       added new iboundn array for density
!       new calculation of C10-C13 for improved numerical stability
!       added density diffusion (proportional to dt) namelist variable denm
!       added most of pressure coding from Nate
!       added semi-implicit terms bzerosi,isemii
!       added new subroutine newvarbv as needed for symmetry boundary conditions
!       jsym option tested for reconnection problem
!       hyper resistivity does not contain density dependence
!       hyper resistivity term added to ohmic heating
!       g14,g15,g16 added for viscous heating
!       boundary condition for CHI to fix one value removed in boundaryv
!       smoother3 changed to set com=0 on boundary (and fix one value of CHI)
!       added itaylor=5 for blob propagation
!       added amuc for compressible viscosity
!       fixed bugs with symmetry option and computer time graphic
!       uncommented lines setting Fact=2 for density and pressure advance
!       went back to old form for viscous heating term
!       o6 option for GEM test problem
!       sphie defined as full pressure source term for ipres=0
!       moved allocation of superlu arrays to after opdef
!       removed b1vecini and b1vecini arrays
!       added regularization term "regular" for chi equation
!
!     NMF contributions (4/24/06):
!     ----------------------------
!     added full pressure evolution
!     added boundary conditions for pressure 
!     parallelized and reordered do loops in newvarb, newvar
!     gam defaults to 5/3, but input file definition no longer overridden.
!     added gyroviscosity
!     moved calculation of itaylor=4 eigenvectors to wave_perturbation
!     includes gravity
!     split r4 into linear (r4) and nonlinear (r40) parts
!     added gravitational instability test (itaylor=5)
!     removing implicit varible declarations...
!     moved superlu allocations/setup out of onestep
!     direct writing to superlu matrix in ludef* and newvar*
!     changed toroidal flux correction
!
#define REAL64 real
!
      module p_data
!
!.....contains all parameters
!     maxnumvar is maximum number of velocity variables
!     maxi is highest degree polynomial kept in nonlinear calculations
!     ntimep is maximum number of timesteps
!     irs is the linear resolution of the plot files
!     ioddm is the maximum number of triangles in cylinder calculation
!     maxfiles is maximum number of scratch disk files (now 50 are used)
!     maxhdf5 is the maximum number of HDF5 files
!     maxplots is the maximum dimension of the graph array
!
      implicit none
      integer  maxi,ntimep,irs,ntri,ntridim,ioddm,maxfiles,             &
     &         maxhdf5,maxplots
      parameter(maxplots=60,maxi=20,ntimep=10000,irs=201,ntridim=1,     &
     &          ioddm=2, maxfiles=100, maxhdf5=30)
      character*11, parameter :: version="M3D-C1 12.2"
      end module p_data
!
      module basic
      use p_data

      real alx,alxp,alz,alzp,xzero,zzero,thimp,amu,etar,dt,p1,p2,p0,pi0, &
     &     tcuro,djdpsi,xmag,zmag,xlim,zlim,facw,facd,db,cb,            &
     &     bzero,hyper,hyperi,hyperv,hyperc,gam,beta,eps,tau,           &
     &     kappa,kappat,hyperp,del,denm,bzerosi,grav,ln,amuc,chipar,    &
     &     regular
      integer n,m,linear,maxs,ntimemax,ntimepr,itor,ioddmx,isetup,      &
     &   irestart,iper,jper,itaylor,itest,isecondorder,imask,nskip,     &
     &        numvar,istart,idens,ipres,itaylorw,isym,jsym,eqsubtract,  &
     &        isemii,elvis, gyro
      character*40 netcdffilename
!
!.....input quantities---defined in subroutine input or in namelist
!
      namelist / inputnl/                                               &
     &        n,m,linear,maxs,ntimemax,ntimepr,itor,ioddmx,isetup,      &
     &   irestart,iper,jper,itaylor,itest,isecondorder,imask,nskip,     &
     &        numvar,istart,idens,ipres,itaylorw,isym,jsym,eqsubtract,  &
     &        isemii, iprint,elvis,gyro,netcdffilename,                 &
     &        alx,alz,xzero,zzero,thimp,amu,amuc,etar,dt,p1,p2,p0,pi0,  &
     &        tcuro,djdpsi,xmag,zmag,xlim,zlim,facw,facd,db,cb,         &
     &        bzero,hyper,hyperi,hyperv,hyperc,gam,beta,eps,tau,        &
     &        kappa,kappat,hyperp,del,denm,bzerosi,grav,ln,chipar,      &
     &        regular
!
!     derived quantities
      real    deex,deez,tt,gamma4,gamma2,gamma3,dpsii,psimin,psilim,pi, &
     &  time,timer,ajmax,thdts,ththmdts,thsqudts,errori,enormi,ratioi,  &
     &  tread, telements, tfirst,ttotal,tsolve,tonestep,tmpi,gbound,    &
     &  tadotb,tinit,tnewvarb,pefac
      integer ni(20),mi(20),nreg,nbcgs,nbcp,nbcv,nbcn,                  &
     &  iboundmax ,ntime,ntimer,nrank,msize,ntimemin,                   &
     &  nlast, iread(maxfiles), ntensor, iframe, ihdf5, idebug,         &
     &  islutype, iprint,itimer, m_loc, nnz_loc, ifs_row
      real    ekin, emag, ekind, emagd, ekino, emago, ekindo, emagdo,   &
     &     ekint,emagt,ekintd,emagtd,ekinto,emagto,ekintdo,emagtdo,     &
     &     ekinp,emagp,ekinpd,emagpd,ekinpo,emagpo,ekinpdo,emagpdo,     &
     &     ekinph,ekinth,emagph,emagth,ekinpho,ekintho,emagpho,emagtho, &
     &     ekin3,ekin3d,ekin3h,emag3,ekin3o,ekin3do,ekin3ho,emag3o,     &
     &     emag3h,emag3d,emag3ho,emag3do,tflux,chierror,totcur, area
      real   xary(0:irs-1),yary(0:irs-1),mif(0:maxhdf5-1),              &
     & maf(0:maxhdf5-1)
!.....define some arrays for transmitting to ElVis files
      real curarray(irs),psiarray(irs),ephiarray(irs),jcbarray(irs),    &
     &     vcbarray(irs),etajarray(irs),hyprarray(irs),densarray(irs),  &
     &     pearray(irs)
      character*8 filename(79)
      character*10 datec, timec
      data filename /8hscratch1, 8hscratch2, 8hscratch3, 8hscratch4,    &
     &               8hscratch5, 8hscratch6, 8hscratch7, 8hscratch8,    & 
     &               8hscratch9, 8hscratc10, 8hscratc11, 8hscratc12,    &
     &               8hscratc13, 8hscratc14, 8hscratc15, 8hscratc16,    & 
     &               8hscratc17, 8hscratc18, 8hscratc19, 8hscratc20,    &
     &               8hscratc21, 8hscratc22, 8hscratc23, 8hscratc24,    &
     &               8hscratc25, 8hscratc26, 8hscratc27, 8hscratc28,    &
     &               8hscratc29, 8hscratc30, 8hscratc31, 8hscratc32,    &
     &               8hscratc33, 8hscratc34, 8hscratc35, 8hscratc36,    &
     &               8hscratc37, 8hscratc38, 8hscratc39, 8hscratc40,    &
     &               8hscratc41, 8hscratc42, 8hscratc43, 8hscratc44,    &
     &               8hscratc45, 8hscratc46, 8hscratc47, 8hscratc48,    &
     &               8hscratc49, 8hscratc50, 8hscratc51, 8hscratc52,    &
     &               8hscratc53, 8hscratc54, 8hscratc55, 8hscratc56,    &
     &               8hscratc57, 8hscratc58, 8hscratc59, 8hscratc60,    &
     &               8hscratc61, 8hscratc62, 8hscratc63, 8hscratc64,    &
     &               8hscratc65, 8hscratc66, 8hscratc67, 8hscratc68,    &
     &               8hscratc69, 8hscratc70, 8hscratc71, 8hscratc72,    &
     &               8hscratc73, 8hscratc74, 8hscratc75, 8hscratc76,    &
     &               8hscratc77, 8hscratc78, 8hscratc79/
!
!.....initialization quantities
      integer ifirstd1_lu, ifirsts1_lu, ifirstd2_lu, ifirsts2_lu,       &
     &        ifirstr1_lu, ifirstr2_lu, ifirstq2_lu, ifirsts3_lu,       &
     &        ifirsts4_lu, ifirsts5_lu, ifirsts6_lu, ifirsts7_lu,       &
     &        ifirsts9_lu, ifirstd9_lu, ifirstq9_lu, ifirstr9_lu,       &
     &        ifirsts8_lu, ifirstd8_lu, ifirstq8_lu, ifirstr8_lu,       &
     &        ifirstsa_lu
!
      data mi /0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0/
      data ni /0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5/
      data iframe /0/
!
!.....MPI variable(s)
      integer myrank, maxrank
      end module basic
!
      module t_data
      use p_data
!.....variables used to define the triangles
      real   atri(ioddm),btri(ioddm),ctri(ioddm),ttri(ioddm),           &
     &  gtri(20,18,ioddm)
      real     , allocatable:: xcord(:),zcord(:),rinv(:)
      integer  , allocatable:: ist(:,:),ist0(:,:),incb(:,:),itypemax(:)
      integer  , allocatable:: jbiglist(:,:)
      integer  , allocatable:: origcolind(:), origcolind1(:)
      end module t_data
!
      module arrays
      use p_data

!     indices
      integer p,q,r,s
      integer, allocatable:: iboundgs(:), iboundv(:), iboundp(:),       &
     & iboundn(:), isvaln(:,:),isval0(:,:),isval1(:,:),isval2(:,:)  
      real    fint(-6:maxi,-6:maxi), xi(3),zi(3),df(0:4,0:4),           &
     &        xsep(5), zsep(5), graphit(0:ntimep,maxplots)
      real, allocatable:: psibounds(:), velbounds(:), combounds(:) 
!
!.....arrays defined at all vertices
#ifndef BIT64
      integer, parameter :: r8 = selected_real_kind(12,100)
      real(r8) ,allocatable::                                           &
#else
      REAL64, allocatable::
#endif
     &        vel(:), vels(:), veln(:),                                 &
     &        velold(:), vel0(:), vel1(:),                              &
     &        phi(:), phis(:),                                           &
     &        phiold(:), phi0(:), phi1(:), phip(:),                     &
     &       jphi(:),jphi0(:),sphi1(:),sphi2(:),vor(:),vor0(:),         &
     &       com(:),com0(:),den(:),den0(:),dent(:),deni(:),denold(:),   &
     &       pres(:),pres0(:),r4(:),q4(:),r4temp(:),bsi(:),r40(:),      &
     &       b1vector(:), b2vector(:), b3vector(:),q34(:),q34temp(:),   &
     &       vtemp(:),phistar(:),ohmic(:),viscous(:),                   &
     &       phitimen(:), veltimen(:),                                  &
     &       fun1(:),fun4(:),fun2(:),fun3(:),sphip(:),sphie(:),         &
     &       galph(:),ggam(:),gmu(:),bsqr(:),bsqri(:),ephi(:),          &
     &       eph2(:),eph3(:),eph4(:),eph5(:),eph6(:),eph7(:),           &
     &       gw1(:),gw2(:),gw3(:),tw1(:),tw2(:),tw3(:),tw4(:),tw5(:),   &
     &       tw40(:),tw50(:),
!
!.....metric terms with 1-2 indices
     &        aterm(:,:,:), bterm(:,:,:), termbb(:,:),                  &
     &        dterm(:,:,:), sterm(:,:,:), d2term(:,:),                  &
     &        xterm(:,:,:), yterm(:,:,:),                               &
!
!.....metric terms involving 3 indices
     &        k0term (:,:,:,:), k1term (:,:,:,:), k2term (:,:,:,:),     &
     &        g0term (:,:,:,:), g2term (:,:,:,:), g3term (:,:,:,:),     &
     &        g4term (:,:,:,:), g5term (:,:,:,:),                       &
     &        g6term (:,:,:,:), g7term (:,:,:,:),                       &
     &        g8term (:,:,:,:), g9term (:,:,:,:),                       &
     &        g10erm (:,:,:,:), g11erm (:,:,:,:),                       &
     &        g12erm (:,:,:,:), g13erm (:,:,:,:), g14erm (:,:,:,:),     &
     &        g15erm (:,:,:,:), g16erm (:,:,:,:),                       &
     &        h3term (:,:,:,:), h5term (:,:,:,:),                       &
     &        x0term (:,:,:,:), y0term (:,:,:,:),                       &
     &        x1term (:,:,:,:), y1term (:,:,:,:),                       &
     &        x2term (:,:,:,:), y2term (:,:,:,:),                       &
!
!     local metric terms involving 2 indices
     &        u1terml(:,:), u2terml(:,:), u4terml(:,:),                 &
     &        u3terml(:,:), u6terml(:,:), u5terml(:,:),                 &
     &        u7terml(:,:), u8terml(:,:), c10erml(:,:),                 &
     &        c11erml(:,:), q0terml(:,:), c12erml(:,:),                 &
     &        c13erml(:,:), c14erml(:,:), c15erml(:,:),                 &
     &        c16erml(:,:), c1term2(:,:), c3term2(:,:),                 &
     &        c4term2(:,:), c5term2(:,:), c6term2(:,:),                 &
     &        c7term2(:,:), c8term2(:,:), c9term2(:,:),                 &
     &        c10erm2(:,:), c11erm2(:,:), c12erm2(:,:),                 &
     &        c13erm2(:,:), c14erm2(:,:), c15erm2(:,:),                 &
     &        c16erm2(:,:),                                             &
     &        v1terml(:,:), v2terml(:,:), v3terml(:,:), v4terml(:,:),   &
     &        v5terml(:,:), v6terml(:,:), v7terml(:,:), v8terml(:,:),   &
     &        v9terml(:,:), v10erml(:,:), v11erml(:,:), v12erml(:,:),   &
     &        v13erml(:,:), v14erml(:,:), v15erml(:,:), v16erml(:,:),   &
!..........gyroviscous
     &     gv11arml(:,:), gv12arml(:,:), gv12brml(:,:), gv12crml(:,:),  &
     &     gv13arml(:,:), gv13brml(:,:), gv13crml(:,:), gv22arml(:,:),  &
     &     gv23arml(:,:), gv23brml(:,:), gv33brml(:,:), gv33crml(:,:),  &
     &     gv12axml(:,:), gv12bxml(:,:), gv12cxml(:,:), gv13bxml(:,:),  &
     &     gv13cxml(:,:), gv23axml(:,:), gv23bxml(:,:),                 &
!..........anisotropic heat flux
     &     t1terml(:,:), t2terml(:,:), t3terml(:,:),                    &
     &     t4terml(:,:), t5terml(:,:),                                  &
!
!.....metric terms involved in toroidal geometry
     &        hterm(:,:,:), termh2(:,:,:,:),                            &
!
!     intermediate metric terms with 1-3 indices                        
     &       d2matrix(:),     bmatrix(:,:),  amatrix(:,:),              &
     &        dmatrix(:,:),  qmatrix(:,:), ttermij(:,:),                &
     &       xmatrix(:,:), ymatrix(:,:),                                &
     &       termg0(:,:,:), termg2(:,:,:), termg3(:,:,:),               &
     &       termg4(:,:,:), termg5(:,:,:), termg6(:,:,:),               &
     &       termg7(:,:,:), termg8(:,:,:), termg9(:,:,:),               &
     &       terg10(:,:,:), terg11(:,:,:), terg12(:,:,:),               &
     &       terg13(:,:,:), terg14(:,:,:), terg15(:,:,:),               &
     &       terg16(:,:,:), termh (:,:,:), termh3(:,:,:),               &
     &       termh5(:,:,:), termk0(:,:,:), termk1(:,:,:),               &
     &       termk2(:,:,:), tensor(:,:,:,:), tterm(:,:,:,:),            &
     &       termx0(:,:,:), termx1(:,:,:), termx2(:,:,:),               &
     &       termy0(:,:,:), termy1(:,:,:), termy2(:,:,:),
!     these contain the data that is given to superlu
!
     &       ss(:,:,:,:), sstemp(:,:,:,:) !,dd(:,:,:,:), rr(:,:,:,:),   &
!     &       qq(:,:,:,:)
      end module arrays
!
      Program Reducedquintic
!
!
!...Program to test the "Reduced Quintic" 2D Finite Elements for the
!
!   Ref:  [1] Strang and Fix, An Analysis of the Finite Element Method, page 83
!         [2] G.R. Cowper, et al, AIAA Journal, Vol 7, no. 10 page 19
!         [3] JCP Vol 147 p 318-336 (1998)
!         [4] JCP Vol 200 p 133-152 (2004)
!         [5] Phys Plasma May 2005
!
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
      use opdef_mod
      use newvar_mod
      use init_conds_mod
      use inout_mod
#ifdef mpi
      use supralu_dist_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif

      integer ier, lx, lz, l, ione, i, j, index1, index2, ifail
      integer maxts, numvars
      real dtmin, ratemin, ratemax, actualhyper
      integer  nran1, msize1, m_lo1, nnz_lo1, ifs_ro1

#ifdef mpi
!....Start up message passing, SuperLU process grid
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
!.....initialize the SUPERLU process grid
      call SLUD_init
      call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
      if (ier /= 0) then
         print *,'Error in MPI_Comm_size:',ier
         call safestop(1)
      endif
#else
      myrank = 0
#endif
!
      if(myrank.eq.0) then
         call date_and_time( datec, timec)
         write(*,1001) version, datec(1:4),datec(5:6),datec(7:8),                &
     &        timec(1:2),timec(3:4),timec(5:8)
 1001    format("VERSION ",a11,  "  DATE: ",a4,1x,a2,1x,a2,3x,           &
     &        "TIME: ",a2,":",a2,":",a4,/)
      endif
!
      pi = acos(-1.)
!.....special switch installed 12/03/04 to write slu matrices
      idebug = 0
      islutype = 0
!
!.....special switch to provide timing information if set
      itimer = 1
!
!.....debug printout
      iprint = 0
!
!.....Choose which test problem to run (0 indicates reading input file C1input)
      itest=0
!
!.....initialize needed variables and define geometry and triangles
      call init
!
!.....define the operators in terms of local coordinates
      call ireaddef
!
      if(ipres.eq.1) then
         pefac = 1.
      else
         pefac = (p0-pi0)/p0
      endif
!
      if(isetup.ne.0) then
         call opdef1
         call opdef2
         if (myrank.eq.0) then
            call writeit
         endif
         if(ntimemax.eq.0) go to 999
      else
         call readit
      endif
!
!     initialize superlu matrices
      nrank = n*m*6*numvar
!.....total number of entries in the sparse matrix (approximate)
      msize = 9*n*m*(6*numvar)**2
      call jbdecomp1(n*m, 6*numvar, m_loc, nnz_loc, ifs_row)

      if(myrank.eq.0) write(*,*) "allocating..."

#ifdef mpi
      if(idens.eq.1 .or. ipres.eq.1) then
         nran1 = n*m*6 
         msize1 = 9*n*m*6**2
         call jbdecomp1(n*m, 6, m_lo1, nnz_lo1, ifs_ro1)
         allocate(origcolind1(nnz_lo1))
      endif
      allocate(jbiglist(m*n,9))
      allocate(origcolind(nnz_loc))

      call sparseR8d_init(s1matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
      call sparseR8d_init(d1matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
      call sparseR8d_init(r1matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
      call sparseR8d_init(s2matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
      call sparseR8d_init(d2matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
      call sparseR8d_init(r2matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
      call sparseR8d_init(q2matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
      call sparseR8d_setup(s1matrix_lu,numvar,jbiglist)
      call sparseR8d_setup(d1matrix_lu,numvar,jbiglist)
      call sparseR8d_setup(r1matrix_lu,numvar,jbiglist)
      call sparseR8d_setup(s2matrix_lu,numvar,jbiglist)
      call sparseR8d_setup(d2matrix_lu,numvar,jbiglist)
      call sparseR8d_setup(r2matrix_lu,numvar,jbiglist)
      call sparseR8d_setup(q2matrix_lu,numvar,jbiglist)
      call sparseR8d_new(s1matrix_lu, ier)
      call sparseR8d_new(s2matrix_lu, ier)
      if(idens.eq.1) then
        call sparseR8d_init(s8matrix_lu,nran1,nnz_lo1,m_lo1,ifs_ro1,ier)
        call sparseR8d_init(d8matrix_lu,nran1,nnz_lo1,m_lo1,ifs_ro1,ier)
        call sparseR8d_init(r8matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
        call sparseR8d_init(q8matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
        call sparseR8d_setup(s8matrix_lu,1,jbiglist)
        call sparseR8d_setup(d8matrix_lu,1,jbiglist)
        call sparseR8d_setup(r8matrix_lu,numvar,jbiglist)
        call sparseR8d_setup(q8matrix_lu,numvar,jbiglist)
      call sparseR8d_new(s8matrix_lu, ier)
      endif
      if(ipres.eq.1) then
        call sparseR8d_init(s9matrix_lu,nran1,nnz_lo1,m_lo1,ifs_ro1,ier)
        call sparseR8d_init(d9matrix_lu,nran1,nnz_lo1,m_lo1,ifs_ro1,ier)
        call sparseR8d_init(r9matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
        call sparseR8d_init(q9matrix_lu,nrank,nnz_loc,m_loc,ifs_row,ier)
        call sparseR8d_setup(s9matrix_lu,1,jbiglist)
        call sparseR8d_setup(d9matrix_lu,1,jbiglist)
        call sparseR8d_setup(r9matrix_lu,numvar,jbiglist)
        call sparseR8d_setup(q9matrix_lu,numvar,jbiglist)
        call sparseR8d_new(s9matrix_lu, ier)
      endif

      origcolind = s1matrix_lu%colind
      if(idens.eq.1) then 
         origcolind1 = s8matrix_lu%colind
      else if(ipres.eq.1) then
         origcolind1 = s9matrix_lu%colind
      endif
#endif

      if(myrank.eq.0) write(*,*) "done allocating."

      if (ier.ne.0) call safestop(5)

!
!.....calculate the RHS (forcing function)
      call rhsdef
!  
!     end of loop over all triangles defining matrices   
!.....note that the matrix indices now refer to vertices, not triangles.
!................................................................
!
!     initialize the solution to an equilibrium and save
!
!     for the wave propagation test, output calculated frequency and 
!     eigenvectors
      if((itaylor.eq.4).and.(myrank.eq.0)) then
         write(*,*) "calculating eigenmodes"
         call wave_perturbation(0,0,0,5)
      endif
!
      if(irestart.eq.1) then
         call rdrestart
         if(istart.ne.0) then
            ntimer = 0
            timer = 0.
         else
            tfirst = tfirst -graphit(ntimer,33)
            tread = graphit(ntimer,35)
            telements = graphit(ntimer,36) - tread
            tsolve = graphit(ntimer,37) - telements - tread
            tonestep = graphit(ntimer,38)
            tmpi = graphit(ntimer,39) - tsolve - telements - tread           &
     &           - tadotb - tinit
            tadotb = graphit(ntimer,40) - tsolve - telements - tread
            tinit  = graphit(ntimer,41) - tsolve - telements - tread         &
     &           - tadotb
!>>>DIAG
            write(*,6001) tfirst,tread,telements,tsolve,tonestep
 6001       format("tfirst,tread,telements,tsolve,tonestep",/,1p5e12.4)
         endif
      else                      !  start of branch to initialize (no restart)
         ntimer = 0
         timer = 0.
         if(idens.eq.1) then
            call denequ(den0)
            if(myrank.eq.0) call oneplot(0,0,den0,1,1,"den0 ")
         endif                  ! on idens .eq.1
         if(ipres.eq.1) then
            call presequ(pres0)
            if(myrank.eq.0) call oneplot(0,0,pres0,1,1,"p0   ")
         endif                  ! on ipres.eq.1
         if(itor.eq.0) then
            velold = 0
            phiold = 0
            call velequ(velold)
            call phiequ(phiold)
         else
            numvars = numvar
            call gradshafranov
            numvar = numvars
            do i=1,n*m
               do j=1,6
                  index2 = (i-1)*6*numvar + j
                  index1 = (i-1)*6 + j
                  phiold(index2) = vel(index1)
               enddo
            enddo
         endif
         vel0 = velold
         phi0 = phiold 
         denold = den0

!......calculate initial perturbed fields
         call velinit(vel)
         call phiinit(phi)
         if(idens.eq.1) then
            call deninit(den)
         endif
         if(ipres.eq.1) then
            call presinit(pres)
         endif
!
!......for nonlinear problem with eqsubtract=0, put equilibrium in perturbation
         if(linear.eq.0 .and. eqsubtract.eq.0) then
            vel = vel + vel0
            vel0 = 0
            phi = phi + phi0
            phi0 = 0
            if(idens.eq.1) then
               den = den + den0
               den0 = 0
            endif
         endif     
!     
!        calculate the equilibrium current density
         if(itaylor.ne.4 .and. itaylor.ne.6) then
            call newvarb(phi0,jphi0,m*n,numvar,1,1)
            if(myrank.eq.0) call oneplot(0,0,jphi0,1,1,"jphi0")
         else
            jphi0 = 0
         endif

!     
!        calculate the equilibrium vorticity and velocity divergence
         vor0 = 0.
         com0 = 0.

!......save initial conditions for analysis
         phis = phi
         vels = vel
      endif                     !  end of the branch on restart/no restart
!     
!     define ficticious field array for semi-implicit method
         if(isemii.gt.0) then
            do lx=1,n
               do lz=1,m
                  l = lx + (lz-1)*n
                  ione = (l-1)*6
                  bsi(ione+1) = bzerosi
                  do i=2,6
                     bsi(ione+i) = 0
                  enddo
               enddo
            enddo
         endif                  ! on isemii
!
      dent = den + den0

      if(chipar.ne.0) then
         call newvar(phi,bsqr,m*n,numvar,1,8)
         call inverse(bsqr,bsqri)
         call inverse(dent,deni)

         if(idens.eq.1) then
            call newvar(phi,tw1,m*n,numvar,1,15)
         else
            tw1 = 0.
         endif
         call newvar(phi,tw2,m*n,numvar,1,16)
         call newvar(phi,tw3,m*n,numvar,1,17)
         call newvar(phi+phi0,tw4,m*n,numvar,3,18)
         call newvar(phi+phi0,tw5,m*n,numvar,3,19)
         call newvar(phi0,tw40,m*n,numvar,3,18)
         call newvar(phi0,tw50,m*n,numvar,3,19) 

         call q34def
      endif
!
!
!      call axis(phi,xsep,zsep,0)   
!
!.....start the time dependent loop
!
      ier = 0
      nlast = 0
      ifail = 0
!
      time = timer
      errori = 0.
      enormi = 0.
      ratioi = 0.
      dtmin = 0.001*dt
      ntime = ntimer
      if(irestart.ne.1) then
        call energy
        if (myrank.eq.0) call output
      endif
      if(ntimemax.le.ntimer) go to 101
      if(myrank.eq.0) then
         call plotit(vel,phi,0)
         call plotit(vel0,phi0,0)
      endif

      do ntime=ntimer+1,ntimemax
         time = time + dt
         if(myrank.eq.0 .and. iprint.eq.1)
     &        write(*,*) "calling onestep..."
         call onestep
         if(myrank.eq.0) call second(ttotal)
         if(myrank.eq.0) call output
         if(linear.eq.1) call scaleback
c$$$         if(ekin .gt. 100.) go to 100
         if(myrank.eq.0 .and. elvis.gt.0) call netcdfout
      enddo ! ntime
!
!
      go to 102
 100  continue
      ntimemax = ntime
      if(myrank.eq.0 .and. elvis.gt.0) call netcdfout
  102 continue

      if(myrank.eq.0 .and. iprint.eq.1) then
         print *, 'Deallocating matrices...'
      endif
#ifdef mpi
      deallocate(jbiglist)
      deallocate(origcolind)

      if(idens.eq.1 .or. ipres.eq.1) deallocate(origcolind1)
      
      call sparseR8d_del(s1matrix_lu,ier)
      call sparseR8d_free(s1matrix_lu,ier)
      call sparseR8d_del(s2matrix_lu,ier)
      call sparseR8d_free(s2matrix_lu,ier)
      if(ifirsts3_lu.ne.0) then
         call sparseR8d_del(s3matrix_lu,ier)
         call sparseR8d_free(s3matrix_lu,ier)
      endif
      if(ifirsts5_lu.ne.0) then
         call sparseR8d_del(s5matrix_lu,ier)
         call sparseR8d_free(s5matrix_lu,ier)
      endif
      if(ifirsts6_lu.ne.0) then
         call sparseR8d_del(s6matrix_lu,ier)
         call sparseR8d_free(s6matrix_lu,ier)
      endif
      if(ifirsts7_lu.ne.0) then
         call sparseR8d_del(s7matrix_lu,ier)
         call sparseR8d_free(s7matrix_lu,ier)
      endif
      if(ifirstsa_lu.ne.0) then
         call sparseR8d_del(samatrix_lu,ier)
         call sparseR8d_free(samatrix_lu,ier)
      endif
      call sparseR8d_free(d1matrix_lu,ier)
      call sparseR8d_free(d2matrix_lu,ier)
      call sparseR8d_free(d3matrix_lu,ier)
      call sparseR8d_free(r1matrix_lu,ier)
      call sparseR8d_free(r2matrix_lu,ier)
      call sparseR8d_free(q2matrix_lu,ier)
      if(idens.eq.1) then
         call sparseR8d_del(s8matrix_lu,ier)
         call sparseR8d_free(s8matrix_lu,ier)
         call sparseR8d_free(d8matrix_lu,ier)
         call sparseR8d_free(r8matrix_lu,ier)
         call sparseR8d_free(q8matrix_lu,ier)
      endif
      if(ipres.eq.1) then
         call sparseR8d_del(s9matrix_lu,ier)
         call sparseR8d_free(s9matrix_lu,ier)
         call sparseR8d_free(d9matrix_lu,ier)
         call sparseR8d_free(r9matrix_lu,ier)
         call sparseR8d_free(q9matrix_lu,ier)
      endif
#else 
      call dsupralu_del(s1handle, jer)
      call sparseR8_free(s1matrix_lu,ier)
      call dsupralu_del(s2handle, jer)
      call sparseR8_free(s2matrix_lu,ier)
      call dsupralu_del(s3handle, jer)
      call sparseR8_free(s3matrix_lu,ier)
      call sparseR8_free(d1matrix_lu,ier)
      call sparseR8_free(d2matrix_lu,ier)
      call sparseR8_free(d3matrix_lu,ier)
      call sparseR8_free(r1matrix_lu,ier)
      call sparseR8_free(r2matrix_lu,ier)
      call sparseR8_free(q2matrix_lu,ier)
#endif
      if(myrank.eq.0 .and. iprint.eq.1) then
         print *, 'Done deallocating matrices'
      endif
!
!     
      maxts = ntime-1
!     if(maxts .gt. nlast) then
!........plot linear solution
!        if (myrank.eq.0) call plotit(vel,phi,0)
!        vel1 = vel + vel0
!        phi1 = phi + phi0
!
!........plot full solution
!      if (myrank.eq.0 .and. linear.eq.0 .and. eqsubtract.ne.0) then
!         call plotit(vel1,phi1,1)
!      endif
!      endif
      if(itaylor.eq.4 .and. myrank.eq.0) then
         call errorcheck
      endif
!
 101  continue
!
      ratemin = 0.
      ratemax = 0.
      do ntime=ntimemin,maxts
         ratemin = min(ratemin,graphit(ntime,26))
         ratemax = max(ratemax,graphit(ntime,26))
      enddo
      actualhyper = hyper*deex*deex*etar
!
      if (myrank.eq.0) then
      write(*,5002) ttotal-tfirst,                                      &
     &  n,m,numvar,amu,etar,dt,thimp,cb,db,hyper,hyperi,hyperv,hyperc,  &
     &  actualhyper,ratemin,ratemax,ajmax,graphit(ntimemax,25)
      write(*,5003) linear, itaylor, itaylorw, iper, jper,              &
     &     isetup, imask, irestart, alx,alz,facd,bzero,eps,amuc
!
      write(9,5002) ttotal-tfirst,                                      &
     &  n,m,numvar,amu,etar,dt,thimp,cb,db,hyper,hyperi,hyperv,hyperc,  &
     &  actualhyper,ratemin,ratemax,ajmax,graphit(ntimemax,25)
      write(9,5003) linear, itaylor, itaylorw, iper, jper,              &
     &     isetup, imask,  irestart,alx,alz,facd,bzero,eps,amuc
      endif
!
      if(myrank.eq.0) write(*,*) "telements, tsolve, tnewvarb",         &
     &     telements, tsolve, tnewvarb
!
      if(ntime.gt.1 .and. myrank.eq.0) then
         call plotenergy(graphit,ntimep,maxts,ntimemin,numvar)
         open(99,file='C1graphit',form='formatted',status='unknown')
         write(99,8001) maxts
         do i=1,maxts
           write(99,8002) (graphit(i,j),j=1,maxplots)
         enddo
 8001    format(i5)
 8002    format(1p10e12.4)
      endif
 999  continue
!
!
      if (myrank.eq.0) call plote()
      call safestop(2)
      write(*,*) myrank
!
 5002 format(" ttotal =", 1pe11.4,   "  n,m,numvar =", 0p3i4,           &
     &  "   amu,etar =", 1p2e10.2,  /,"  dt,thimp =",1p2e10.2,          &
     &  "  cb,db=",1p2e12.2,/,                                          &
     &  "  hyper,hyperi,hyperv,hyperc,actual hyper = ",  1p5e12.4,  /,  &
     &  "  ratemin,ratemax,ajmax,reconflux = ",1p4e12.4) 
 5003 format(" linear, itaylor(2),iper,jper,isetup, imask, irestart: ", &
     & 8i4, / ," alx, alz, facd, bzero, eps ", 1p5e12.4)
 2323 format(1p5e12.4)
 2121 format(" start of time dependent loop")
 5001 format("ifail .ne. 0 after call to f04aaf")
 7011 format(1p6e20.12)
 7012 format(6e20.12)
!
!
      end program Reducedquintic
!============================================================
      subroutine onestep
!
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
      use newvar_mod

      implicit none
      real tstartonestep, tendonestep, tbeforeadotb, tafteradotb
      real taftersolve
      integer itri, iodd, irect, jrect, i, j, jone, ivertex
      integer i0, ione, l, ier, jer, lx, lz
      real eerror, ediff, etotd, sum
      real correction, index, j2

      if(myrank.eq.0 .and. itimer.gt.0) call second(tstartonestep)
!
!.....advance solution to time ntime+1
!
      ththmdts = thimp*(thimp-1)*dt**2
      thsqudts = thimp**2*dt**2
      thdts = thimp*dt**2
!
      if(iprint.eq.1 .and. myrank.eq.0) write(*,*) "starting onestep"
!
!.....define the inverse density array deni
      if(idens.eq.1) then
         dent = den + den0
         call inverse(dent,deni)
      endif
!
!.....define the current vector jphi and the RHS vectors sphi1 and sphi2
      call newvarb(phi,jphi,m*n,numvar,1,1)
      call newvarb(phi,sphi1,m*n,numvar,1,4)
      if(numvar.ge.2) call newvarb(phi,sphi2,m*n,numvar,1,5)
!
!.....Define magnetic field strength variables
      if(gyro.eq.1 .or. chipar.ne.0) then
         call newvar(phi,bsqr,m*n,numvar,1,8)
         call inverse(bsqr,bsqri)
      endif

!.....Define gyroviscous quantities
!.....-----------------------------
      if(gyro.eq.1) then
!.....   define the other relevant gyroviscous quantities
         call newvar(phi,galph,m*n,numvar,1,9)
         call newvar(phi,ggam,m*n,numvar,1,10)
         call newvar(phi,gw1,m*n,numvar,1,12)
         if(numvar.ge.2) then
            call newvar(phi,gmu,m*n,numvar,1,11)
            call newvar(phi,gw2,m*n,numvar,1,13)
            call newvar(phi,gw3,m*n,numvar,1,14)
         endif
      endif ! gyro
!.....----------------------------

!
      totcur = 0
      area = 0.
      tflux = 0.
!.....calculate the total perturbed current and area, and toroidal flux
      itri = 0
      do iodd=1,ioddmx
         do jrect=1,m-1+jper
            do irect=1,n-1+iper
               itri = itri + 1
               do j=1,18
                  jone = isval1(itri,j)
                  totcur = totcur + d2term(iodd,j)*jphi(jone)*(1+jsym)
!     
                  if(numvar.ge.2) then
                     j2 = isvaln(itri,j) + 6
                     tflux = tflux + d2term(iodd,j)*(1+jsym)            &
     &                    *(phi(j2) + phi0(j2) - phiold(j2))
                  endif
!     
               enddo
               do j=1,13,6
                  area = area + d2term(iodd,j)
               enddo
            enddo               ! loop over irect
         enddo                  ! loop over jrect
      enddo                     ! loop over iodd
     
!.....adjust toroidal field and boundary condition to conserve toroidal flux
!     if(numvar.ge.2) then
!        correction = tflux/area
!        gbound = gbound - correction
!        ivertex = 0
!        do jrect=1,m
!           do irect=1,n
!              ivertex = ivertex + 1
!              index = 6*numvar*(ivertex-1) + 7
!              phi(index) = phi(index) - correction
!           enddo
!        enddo
!     endif
!
! Advance VEL ===========================================
! 
      if(iprint.eq.1 .and. myrank.eq.0) write(*,*) "advancing vel"
!.....Calculate LU decomposition of velocity matrix if needed
      if(ntime.le.ntimer+1.or.                                          &
     &     (linear.eq.0 .and. mod(ntime,nskip).eq.0)) then
         call ludefvel
!
!......special option to output matrix and quit
       if(idebug.eq.1 .and. islutype.eq.1) then
       call sparseR8d_output(s1matrix_lu)
       call safestop(0)
       endif
!
      endif
      if(grav.ne.0. .and. idens.gt.0) call r4def
      if(myrank.eq.0) call second(tbeforeadotb)
!
!     create new phi vector with total pressure (instead of pe)
      phip = phi
      if((numvar.ge.3).and.(ipres.eq.1)) then
         ivertex = 0
         do jrect=1,m
            do irect=1,n
               ivertex = ivertex + 1
               i0 = 6*numvar*(ivertex-1)
               ione = 6*(ivertex-1) + 1
               phip(i0+13) = pres(ione)
               phip(i0+14) = pres(ione+1)
               phip(i0+15) = pres(ione+2)
               phip(i0+16) = pres(ione+3)
               phip(i0+17) = pres(ione+4)
               phip(i0+18) = pres(ione+5)
            enddo
         enddo
      endif
!     
!.....multiply r1matrix_lu x time n fields, put result in b1vector
      if(iprint.ge.1) write(*,*) "before sparseR8_A_dot_X"
#ifdef mpi
      call sparseR8d_A_dot_X(r1matrix_lu,phip,b1vector,ier)
#else
      call sparseR8_A_dot_X(r1matrix_lu,phip,b1vector,ier)
#endif
      if(ier.ne.0) then
         write(*,*) "ier after A_dot,   ",ier
         call safestop(34)
      endif
!     
!.....multiply d1matrix_lu x time n velocity, add RHS vector to obtain vtemp
      vtemp = 0.
      if(iprint.eq.1)write(*,*) "before second sparseR8_A_dot_X"
#ifdef mpi
      call sparseR8d_A_dot_X(d1matrix_lu,vel,vtemp,ier)
#else
      call sparseR8_A_dot_X(d1matrix_lu,vel,vtemp,ier)
#endif
      if(ier.ne.0) then
         write(*,*) "ier after A_dot,   ",ier
         call safestop(35)
      endif

      vtemp = vtemp +   b1vector + r4 + r40

      do l=1,nbcv
         vtemp(iboundv(l)) = velbounds(l)
      enddo
      if(myrank.eq.0) call second(tafteradotb)
      tadotb = tadotb + tafteradotb-tbeforeadotb
!     
!.....solve linear system with rhs in vtemp (note LU-decomp done first time)
      if(iprint.eq.1) write(*,*) "before dsupralu_solve_s1handle" 
#ifdef mpi
      call sparseR8d_solve(s1matrix_lu, vtemp, jer)

      if(jer.ne.0) then
         write(*,*) 'after sparseR8d_solve', jer
         call safestop(42)
      endif
#else
      call dsupralu_solve(s1handle, vtemp, jer)
      if(jer.ne.0) then
         write(*,*) 'after dsupralu_solve', jer
         call safestop(43)
      endif
#endif
      if(myrank.eq.0) call second(taftersolve)
      tsolve = tsolve + taftersolve - tafteradotb
!
!
!.....coding to calculate the error in the delsquared chi equation
      chierror = 0
      sum = 0
      if(numvar.ge.3) then
         call newvar(vtemp,com,m*n,numvar,3,3)
         itri = 0
         do iodd=1,ioddmx
            do jrect=1,m-1+jper
               do irect=1,n-1+iper
                  itri = itri + 1
                  do j=1,18
                     jone = isval1(itri,j)
                     chierror = chierror + d2term(iodd,j)*com(jone)
                  enddo
               enddo            ! loop over irect
            enddo               ! loop over jrect
         enddo                  ! loop over iodd
!
!
         if(hyperc.gt.0) then
            call smoother3(com,vtemp,m*n,numvar,3)
            call newvar(vtemp,com,m*n,numvar,3,3)
         endif
!
      endif
!
      call newvarbv(vtemp,vor,m*n,numvar,1,2)
!     
      if(hyperc.gt.0) then
!     
!.......calculate vorticity, apply smoothing operator, and redefine vor array
         call smoother1(vor,vtemp,m*n,numvar,1)
         call newvarbv(vtemp,vor,m*n,numvar,1,2)
!
      endif
!     
!.....new velocity solution at time n+1 (or n* for second order advance)
      velold = vel
      if(itaylor.ne.7 .or. itaylorw.ne.1) then
        vel = vtemp
      endif
!     
! Advance DENSITY ========================================================
      if(iprint.eq.1 .and. myrank.eq.0) write(*,*) "advancing den"
      if(idens.eq.1) then
         if(ntime.le.ntimer+1.or.                                           &
     &        (linear.eq.0 .and. mod(ntime,nskip).eq.0)) then
            call ludefden
         endif
         if(iprint.ge.1) write(*,*) "s8handle"
         if(myrank.eq.0) call second(tbeforeadotb)
#ifdef mpi
!     
!.....multiply r8matrix_lu x time n+1 velocity, put result in b2vector
         call sparseR8d_A_dot_X(r8matrix_lu,vel,b2vector,ier)
!
!.....multiply q8matrix_lu x time n velocity, put result in b3vector
         call sparseR8d_A_dot_X(q8matrix_lu,velold,b3vector,ier)

!.....multiply d8matrix_lu x time n fields, add RHS vector to obtain vtemp
         vtemp = 0.
!
         call sparseR8d_A_dot_X(d8matrix_lu,den,vtemp,ier)
#else
!     
!.....multiply r8matrix_lu x time n+1 velocity, put result in b2vector
         call sparseR8_A_dot_X(r8matrix_lu,vel,b2vector,ier)
!     
!.....multiply q8matrix_lu x time n velocity, put result in b3vector
         call sparseR8_A_dot_X(q8matrix_lu,velold,b3vector,ier)
!
!.....multiply d8matrix_lu x time n fields, add RHS vector to obtain vtemp
         vtemp = 0.
!
         call sparseR8_A_dot_X(d8matrix_lu,den,vtemp,ier)
!
#endif
         do lx=1,n
            do lz=1,m
               l = lx + (lz-1)*n
               do i=1,6
                  vtemp(6*(l-1)+i) =  vtemp(6*(l-1)+i)                  &
     &                 + b2vector(6*((l-1)*numvar)+i)                   &
     &                 + b3vector(6*((l-1)*numvar)+i)
               enddo
            enddo
         enddo
!     set density boundary conditions
         do l=1,nbcn
            vtemp(iboundn(l)) = 0.
            if(linear.eq.0 .and. eqsubtract.eq.0) then
               vtemp(iboundn(l)) = vtemp(iboundn(l))+denold(iboundn(l))
            endif
         enddo
         if(myrank.eq.0) call second(tafteradotb)
         tadotb = tadotb + tafteradotb-tbeforeadotb
!
#ifdef mpi
!     
!.....solve linear system...LU decomposition done first time
         call sparseR8d_solve(s8matrix_lu, vtemp, jer)

         if(jer.ne.0) then
            write(*,*) 'after 2nd sparseR8d_solve', jer
            call safestop(29)
         endif
#else
!
!.....solve linear system using LU decomp encoded in s8handle, rhs in vtemp
      if(iprint.eq.1) write(*,*) "before dsupralu_solve ... s8handle",  &
     &        s8handle
         call dsupralu_solve(s8handle, vtemp, jer)
         if(jer.ne.0) then
            write(*,*) 'after 2nd dsupralu_solve', jer
            call safestop(29)
         endif
#endif
         if(myrank.eq.0) call second(taftersolve)
         tsolve = tsolve + taftersolve-tafteradotb
!
!.....new field solution at time n+1 (or n* for second order advance)
         den = vtemp
      endif                     ! on density advance
!
! Advance PRESSURE =====================================================
      if(ipres.eq.1) then
         if(iprint.eq.1 .and. myrank.eq.0) write(*,*) "advancing pres"
         call newvarb(phi,sphip,m*n,numvar,1,7)
         if(ntime.le.ntimer+1.or.                                       &
     &        (linear.eq.0 .and. mod(ntime,nskip).eq.0)) then
            call ludefpres   
         endif
         if(iprint.ge.1) write(*,*) "s9handle"
#ifdef mpi
!
!........multiply r9matrix_lu x time n+1 velocity, put result in b2vector
         call sparseR8d_A_dot_X(r9matrix_lu,vel,b2vector,ier)
!
!........multiply q9matrix_lu x time n velocity, put result in b3vector
         call sparseR8d_A_dot_X(q9matrix_lu,velold,b3vector,ier)
!
!........multiply d9matrix_lu x time n fields, add RHS vector to obtain vtemp
         vtemp = 0.

         call sparseR8d_A_dot_X(d9matrix_lu,pres,vtemp,ier)
#else
!
!........multiply r9matrix_lu x time n+1 velocity, put result in b2vector
         call sparseR8_A_dot_X(r9matrix_lu,vel,b2vector,ier)
!
!........multiply q9matrix_lu x time n velocity, put result in b3vector
         call sparseR8_A_dot_X(q9matrix_lu,velold,b3vector,ier)
!
!........multiply d9matrix_lu x time n fields, add RHS vector to obtain vtemp
         vtemp = 0.
!
         call sparseR8_A_dot_X(d9matrix_lu,pres,vtemp,ier)
#endif
         do lx=1,n
            do lz=1,m
               l = lx + (lz-1)*n
               do i=1,6
                  vtemp(6*(l-1)+i) =  vtemp(6*(l-1)+i)                  &
     &                 + b2vector(6*((l-1)*numvar)+i)                   & 
     &                 + b3vector(6*((l-1)*numvar)+i)                   &
     &                 + q4(6*((l-1)*numvar)+i)
               enddo
            enddo
         enddo
!
!........set boundary conditions
         do l=1, nbcn
            vtemp(iboundn(l)) = 0.
         enddo
!
#ifdef mpi         
!
!........solve linear system...LU decomposition done first time
         call sparseR8d_solve(s9matrix_lu, vtemp, jer)
         if(jer.ne.0) then
            write(*,*) 'after 2nd sparseR8d_solve', jer
            call safestop(29)
         endif
#else
!
!........solve linear system using LU decomp encoded in s9handle, rhs in vtemp
         if(iprint.eq.1) write(*,*)                                     &
     &        "before dsupralu_solve ... s9handle", s9handle
         call dsupralu_solve(s9handle, vtemp, jer)
         if(jer.ne.0) then
            write(*,*) 'after 2nd dsupralu_solve', jer
            call safestop(29)
         endif
#endif
!
!........new field solution at time n+1 (or n* for second order advance)
         pres = vtemp
      endif
!

! Advance PHI =========================================================
      if(iprint.eq.1 .and. myrank.eq.0) write(*,*) "advancing phi"

      if(numvar.ge.3) then
!        anisotropic heat flux terms
         if(chipar.ne.0) then
            if(idens.eq.1) then
               call newvar(phi,tw1,m*n,numvar,1,15)
            else
               tw1 = 0.
            endif
            call newvar(phi,tw2,m*n,numvar,1,16)
            call newvar(phi,tw3,m*n,numvar,1,17)
            call newvar(phi+phi0,tw4,m*n,numvar,3,18)
            call newvar(phi+phi0,tw5,m*n,numvar,3,19)
            call newvar(phi0,tw40,m*n,numvar,3,18)
            call newvar(phi0,tw50,m*n,numvar,3,19)
         endif

!        electron pressure source term
         if(ipres.eq.0) then
            ! sphie = full pressure source term
            call newvar(phi,sphie,m*n,numvar,1,7)
         else
            ! sphie = electron pressure source term
            call newvar(phi,sphie,m*n,numvar,1,6)
         endif

!        assemble q34 matrix term
         call q34def
      endif
!
!.....Calculate LU decomposition of field matrix if needed
      if(ntime.le.ntimer+1.or.                                             &
     &     (linear.eq.0 .and. mod(ntime,nskip).eq.0)) then
         call ludefphi         
      endif
      if(iprint.ge.1) write(*,*) "s2handle" 
      if(myrank.eq.0) call second(tbeforeadotb)
#ifdef mpi
!
!.....multiply r2matrix_lu x time n+1 velocity, put result in b2vector
      call sparseR8d_A_dot_X(r2matrix_lu,vel,b2vector,ier)
!
!.....multiply q2matrix_lu x time n velocity, put result in b3vector
      call sparseR8d_A_dot_X(q2matrix_lu,velold,b3vector,ier)
!
!.....multiply d2matrix_lu x time n fields, add RHS vector to obtain vtemp
      vtemp = 0.
!
      call sparseR8d_A_dot_X(d2matrix_lu,phi,vtemp,ier)
#else
!
!.....multiply r2matrix_lu x time n+1 velocity, put result in b2vector
      call sparseR8_A_dot_X(r2matrix_lu,vel,b2vector,ier)
!
!.....multiply q2matrix_lu x time n velocity, put result in b3vector
      call sparseR8_A_dot_X(q2matrix_lu,velold,b3vector,ier)
!
!.....multiply d2matrix_lu x time n fields, add RHS vector to obtain vtemp
      vtemp = 0.
!
      call sparseR8_A_dot_X(d2matrix_lu,phi,vtemp,ier)
!

#endif
      if(myrank.eq.0) call second(tafteradotb)
      tadotb = tadotb + tafteradotb-tbeforeadotb
      vtemp = vtemp +  b2vector + b3vector                              &
     &     + dt*(ohmic + viscous)
     &     + dt*q34
      do l=1,nbcp
!......added 01/23/06..for nonlinear problem, add equilibrium to boundary
         vtemp(iboundp(l)) = psibounds(l) 
         if(linear.eq.0 .and. eqsubtract.eq.0) then
            vtemp(iboundp(l)) = vtemp(iboundp(l)) + phiold(iboundp(l))
         endif
      enddo

#ifdef mpi
!
!.....solve linear system...LU decomposition done first time
      call sparseR8d_solve(s2matrix_lu, vtemp, jer)
      if(jer.ne.0) then
         write(*,*) 'after 2nd sparseR8d_solve', jer
         call safestop(29)
      endif
#else
!.....solve linear system using LU decomp encoded in s2handle, rhs in vtemp
      if(iprint.eq.1) write(*,*) "before dsupralu_solve ... s2handle",     &
     &                s2handle
      call dsupralu_solve(s2handle, vtemp, jer)
      if(jer.ne.0) then
         write(*,*) 'after 2nd dsupralu_solve', jer
         call safestop(29)
      endif
#endif

      if(myrank.eq.0) call second(taftersolve)
      tsolve = tsolve + taftersolve-tafteradotb
!
!.....new field solution at time n+1 (or n* for second order advance)
      phi = vtemp
!
!
      call energy
!
!.....calculate some arrays for diagnostic plots
      if(mod(ntime,ntimepr).eq.0 .or. elvis.gt.0) then
      if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "b calls to newvarb"
        call newvarb(phi,ephi,m*n,numvar,1,8)
        call newvarb(phi,eph2,m*n,numvar,1,9)
        call newvarb(phi,eph3,m*n,numvar,1,10)
        call newvarb(phi,eph4,m*n,numvar,1,11)
        call newvarb(phi,eph5,m*n,numvar,1,12)
        call newvarb(phi,eph6,m*n,numvar,1,13)
      if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "a calls to newvarb"
        eph7 = eph2 + eph3
      endif  ! on mod(ntime,ntimepr)
!
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
!
      if(myrank.eq.0 .and. itimer.gt.0) then
         call second(tendonestep)
         tonestep = tonestep +tendonestep - tstartonestep
         write(*,3001) tonestep
 3001    format(" time spent in onestep: ",1pe12.4)
      endif

      return
      end
!============================================================
      real function fact(n)
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
      real function f(m,n,a,b,c)
!
!.....f calculates the double integral over a triangle
!     with base a + b and height c  [see ref 2, figure 1]
!
!     f = Int{si**m*eta**n}d(si)d(eta)
!

        anum = c**(n+1)*(a**(m+1)-(-b)**(m+1))*fact(m)*fact(n)
        denom = fact(m+n+2)
        f = anum/denom
        return
      end function f
!============================================================
      subroutine mask(x,z,factor)
      use basic
      factor = 1.
      if(iper.eq.1 .and. jper.eq.1) return
!
!     radius where the factor becomes 0
      r1 = 2.0
!
!     radius where the factor becomes 1
      r2 = 1.5
      alpha = 6./(r1 - r2)
      r0 = (r1+r2)/2.
!
      if(iper.eq.0) then
        rnorm = sqrt(alx*alz)/4.
        r = sqrt( (x-alx/2.)**2 + (z-alz/2.)**2)/rnorm
      else
        rnorm = alz/4.
        r = abs(z-alz/2.)/rnorm
      endif
!
      factor = (1 + tanh(-alpha*(r-r0)))/2.
!
 1001 format(1p4e12.4)
      return
      end
!============================================================
      subroutine ludefvel
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
      use mesh_mod
#ifdef mpi
      use supralu_dist_mod
      use superlu_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer icount, irect, jrect, iodd, itri, ll, itype
      integer i, i1, i2, i3
      integer j, j01, j1, j2, j3
      integer l, lone, l01, l1, l2, l3
      integer k, kone, k01, k1, k2, k3
      integer ir1, ir2, ir3, ier, jer
      real tbeforeread, tafterread, tafterelements, tbeforeinit
      real tafterinit, tstartludefvel
      real hypvf, hypv, thimpv, factor, dbf
      real ssterm(3,3),ddterm(3,3),rrterm(3,3)
      real terma, termb, termd, gmix, qmix, qmix2, c6ave, pmix, tempvar
      real x, z, pressure

!
      if(myrank.eq.0 .and. itimer.gt.0) call second(tstartludefvel)
!>>>  zero, count, and icon appear never to be used.
!      zero = 0.
!      count = zero
      hypv = hyperv*deex**2
      thimpv = 1.
      icount = 0
!      icon = 1
!
c$$$      if(ifirsts1_lu.eq.1) then
c$$$         call f_Destroy_LU(s1matrix_lu%n, SLUprocgrid,                     &
c$$$     &        s1matrix_lu%LUstruct)
c$$$      endif
      call sparseR8d_reset(s1matrix_lu, ier)
      call sparseR8d_reset(d1matrix_lu, ier)
      call sparseR8d_reset(r1matrix_lu, ier)
      s1matrix_lu%colind = origcolind
      d1matrix_lu%colind = origcolind
      r1matrix_lu%colind = origcolind
!     
!.....form velocity matrices
      r40 = 0.
      veln = vel
      itri = 0

!.......position files with metric terms
      do itype=1,ntensor
         if(iread(itype).eq.1) rewind(itype+100)
      enddo
!     
      do iodd=1,ioddmx
      do i=1,18

      do j=1,18

!
         if(myrank.eq.0 .and. itimer.gt.0) call second(tbeforeread)
! ..........read the 4th rank tensor terms from the disk
         if(iread(05).eq.1) then                                     
            read(105,err=1001,end=1001)                                  &
     &           ir1,ir2,ir3,((v7terml(k,l),l=1,18),k=1,18)
            if(ir1.ne.iodd .or. ir2.ne.i .or. ir3.ne.j) go to 1001
         endif
!     
        if(iread(06).eq.1)                                              &
     &       read(106) ir1,ir2,ir3,((v8terml(k,l),l=1,18),k=1,18)
        if(iread(07).eq.1)                                              &
     &       read(107) ir1,ir2,ir3,((v9terml(k,l),l=1,18),k=1,18)
        if(iread(08).eq.1)                                              &
     &       read(108) ir1,ir2,ir3,((v10erml(k,l),l=1,18),k=1,18)
        if(iread(11).eq.1)                                              &
     &       read(111) ir1,ir2,ir3,((v1terml(k,l),l=1,18),k=1,18)
        if(iread(12).eq.1)                                              &
     &       read(112) ir1,ir2,ir3,((v2terml(k,l),l=1,18),k=1,18)
        if(iread(13).eq.1)                                              &
     &       read(113) ir1,ir2,ir3,((c12erml(k,l),l=1,18),k=1,18)
        if(iread(15).eq.1)                                              &
     &       read(115) ir1,ir2,ir3,((c13erml(k,l),l=1,18),k=1,18)
        if(iread(17).eq.1)                                              &
     &       read(117) ir1,ir2,ir3,((c14erml(k,l),l=1,18),k=1,18)
        if(iread(19).eq.1)                                              &
     &       read(119) ir1,ir2,ir3,((c15erml(k,l),l=1,18),k=1,18)
        if(iread(21).eq.1)                                              &
     &       read(121) ir1,ir2,ir3,((c16erml(k,l),l=1,18),k=1,18)
        if(iread(23).eq.1)                                              &
     &       read(123) ir1,ir2,ir3,((c1term2(k,l),l=1,18),k=1,18)
        if(iread(25).eq.1)                                              &      
     &       read(125) ir1,ir2,ir3,((c6term2(k,l),l=1,18),k=1,18)
        if(iread(27).eq.1)                                              &
     &       read(127) ir1,ir2,ir3,((c7term2(k,l),l=1,18),k=1,18)
        if(iread(29).eq.1)                                              &
     &       read(129) ir1,ir2,ir3,((c8term2(k,l),l=1,18),k=1,18)
        if(iread(31).eq.1)                                              &
     &       read(131) ir1,ir2,ir3,((c9term2(k,l),l=1,18),k=1,18)
        if(iread(33).eq.1)                                              &
     &       read(133) ir1,ir2,ir3,((c10erm2(k,l),l=1,18),k=1,18)
        if(iread(35).eq.1)                                              &
     &       read(135) ir1,ir2,ir3,((c11erm2(k,l),l=1,18),k=1,18)
        if(iread(37).eq.1)                                              &
     &       read(137) ir1,ir2,ir3,((v3terml(k,l),l=1,18),k=1,18)
        if(iread(38).eq.1)                                              &
     &       read(138) ir1,ir2,ir3,((v4terml(k,l),l=1,18),k=1,18)
        if(iread(39).eq.1)                                              &
     &       read(139) ir1,ir2,ir3,((v11erml(k,l),l=1,18),k=1,18)
        if(iread(40).eq.1)                                              &
     &       read(140) ir1,ir2,ir3,((v12erml(k,l),l=1,18),k=1,18)
        if(iread(41).eq.1)                                              &
     &       read(141) ir1,ir2,ir3,((v13erml(k,l),l=1,18),k=1,18)
        if(iread(42).eq.1)                                              &
     &       read(142) ir1,ir2,ir3,((v14erml(k,l),l=1,18),k=1,18)
        if(iread(43).eq.1)                                              &
     &       read(143) ir1,ir2,ir3,((v15erml(k,l),l=1,18),k=1,18)
        if(iread(44).eq.1)                                              &
     &       read(144) ir1,ir2,ir3,((v16erml(k,l),l=1,18),k=1,18)
        if(iread(45).eq.1)                                              &
     &       read(145) ir1,ir2,ir3,((v5terml(k,l),l=1,18),k=1,18)
        if(iread(46).eq.1)                                              &
     &       read(146) ir1,ir2,ir3,((v6terml(k,l),l=1,18),k=1,18)
        if(iread(54).eq.1)                                              &
     &       read(154) ir1,ir2,ir3,((gv11arml(k,l),l=1,18),k=1,18)
        if(iread(55).eq.1)                                              &
     &       read(155) ir1,ir2,ir3,((gv12arml(k,l),l=1,18),k=1,18)
        if(iread(56).eq.1)                                              &
     &       read(156) ir1,ir2,ir3,((gv12brml(k,l),l=1,18),k=1,18)
        if(iread(57).eq.1)                                              &
     &       read(157) ir1,ir2,ir3,((gv12crml(k,l),l=1,18),k=1,18)
        if(iread(58).eq.1)                                              &
     &       read(158) ir1,ir2,ir3,((gv13arml(k,l),l=1,18),k=1,18)
        if(iread(59).eq.1)                                              &
     &       read(159) ir1,ir2,ir3,((gv13brml(k,l),l=1,18),k=1,18)
        if(iread(60).eq.1)                                              &
     &       read(160) ir1,ir2,ir3,((gv13crml(k,l),l=1,18),k=1,18)
        if(iread(61).eq.1)                                              &
     &       read(161) ir1,ir2,ir3,((gv22arml(k,l),l=1,18),k=1,18)
        if(iread(62).eq.1)                                              &
     &       read(162) ir1,ir2,ir3,((gv23arml(k,l),l=1,18),k=1,18)
        if(iread(63).eq.1)                                              &
     &       read(163) ir1,ir2,ir3,((gv23brml(k,l),l=1,18),k=1,18)
        if(iread(64).eq.1)                                              &
     &       read(164) ir1,ir2,ir3,((gv33brml(k,l),l=1,18),k=1,18)
        if(iread(65).eq.1)                                              &
     &       read(165) ir1,ir2,ir3,((gv33crml(k,l),l=1,18),k=1,18)
        if(iread(66).eq.1)                                              &
     &       read(166) ir1,ir2,ir3,((gv12axml(k,l),l=1,18),k=1,18)
        if(iread(67).eq.1)                                              &
     &       read(167) ir1,ir2,ir3,((gv12bxml(k,l),l=1,18),k=1,18)
        if(iread(68).eq.1)                                              &
     &       read(168) ir1,ir2,ir3,((gv12cxml(k,l),l=1,18),k=1,18)
        if(iread(69).eq.1)                                              &
     &       read(169) ir1,ir2,ir3,((gv13bxml(k,l),l=1,18),k=1,18)
        if(iread(70).eq.1)                                              &
     &       read(170) ir1,ir2,ir3,((gv13cxml(k,l),l=1,18),k=1,18)
        if(iread(71).eq.1)                                              &
     &       read(171) ir1,ir2,ir3,((gv23axml(k,l),l=1,18),k=1,18)
        if(iread(72).eq.1)                                              &
     &       read(172) ir1,ir2,ir3,((gv23bxml(k,l),l=1,18),k=1,18)
!     
        go to 1002
 1001   continue
        write(*,*) ir1,ir2,ir3,iodd,i,j
        call safestop(1001)
 1002   continue
!
!
        if(myrank.eq.0 .and. itimer.gt.0) call second(tafterread)
        tread = tread + (tafterread-tbeforeread)
!
        do jrect=1,m-1+jper
        do irect=1,n-1+iper
!
           ll = irect + (jrect-1)*(n-1+iper)
           itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll

           i1 = isvaln(itri,i)
           
           ! check whether this node is assigned to this processor
           if(sparseR8d_is_local_row(s1matrix_lu, i1).eq.0) goto 800

           j01= isval0(itri,j)         
           j1 = isvaln(itri,j)
           j2 = j1 + 6
           i2 = i1 + 6
           j3 = j2 + 6
           i3 = i2 + 6

           x = (irect-iodd/3.)*deex
           z = (jrect-1+iodd/3.)*deez
           factor = 1.
           if(imask.eq.1) call mask(x,z,factor)
           hypvf = hypv
!     
           dbf = db*factor
!
           ssterm = 0
           ddterm = 0
           rrterm = 0
!
           terma = aterm(iodd,i,j)
           termb = bterm(iodd,i,j)
           termd = dterm(iodd,i,j)
!
!...NOTE: s1matrix_lu (formed from ss) multiplies the solution at time NTIME+1
!
! LINEAR TERMS

           ssterm(1,1) = terma*(1-idens)      - thimp*dt*amu*termb
           ddterm(1,1) = terma*(1-idens) + (1.-thimp)*dt*amu*termb
!     
           if(numvar.ge.2) then
              ssterm(2,2) = termd*(1-idens) -     thimp *dt*amu*terma   &
     &                 + thimpv *dt*amu*hypvf*termb
              ddterm(2,2) = termd*(1-idens) + (1.-thimp)*dt*amu*terma   &
     &             - (1.-thimpv)*dt*amu*hypvf*termb
           endif
           if(numvar.ge.3) then
              ssterm(3,3) = terma*(1-idens) -   thimp *dt*2.*amuc*termb
     &               - regular*termd
              ddterm(3,3) = terma*(1-idens)+(1.-thimp)*dt*2.*amuc*termb
              rrterm(3,3) = - dt*terma
           endif
!
!
! SIMPLE NONLINEAR TERMS
           do k=1,18
              k1= isvaln(itri,k)
              k01 = isval0(itri,k)
              k2= k1 + 6
              k3= k2 + 6
              kone = isval1(itri,k)
              gmix = g0term(iodd,i,j,k)+g0term(iodd,i,k,j)
              if(ipres.eq.0) then 
                 if(numvar.eq.3) then
                    pressure = phi(k3) + phi0(k3)
                 else
                    pressure = p0 / 18.
                 endif
              else 
                 pressure = pres(kone) + pres0(kone)
              endif
!
              ssterm(1,1) = ssterm(1,1) + dt*gmix*(1-idens)*            &
     &             thimp*(vel(k1)+vel0(k1))
!     
              ddterm(1,1) = ddterm(1,1) + dt*gmix*(1-idens)*            &
     &             (thimp*(vel(k1)+vel0(k1))-(.5*vel(k1) + vel0(k1)))
!             
              if(idens.gt.0) then
                 ssterm(1,1) = ssterm(1,1)+dent(kone)*g2term(iodd,i,j,k)    
                 ddterm(1,1) = ddterm(1,1)+dent(kone)*g2term(iodd,i,j,k)
!     
!                gravity
!                =======
                 ssterm(1,1) = ssterm(1,1) + grav*thsqudts                 &
     &                *x1term(iodd,i,j,k)*dent(kone)
                 ddterm(1,1) = ddterm(1,1) + grav*ththmdts                 &
     &                *x1term(iodd,i,j,k)*dent(kone)
!                =======
              endif
!
              rrterm(1,1)=rrterm(1,1) + dt*gmix*(.5*phi(k1) + phi0(k01)) &
     &             + thdts*gmix*sphi1(kone)
!     
              r40(i1) = r40(i1) + thdts*gmix*sphi1(kone)*phi0(j01)
!
              if(numvar.ge.2) then
!
              ssterm(2,1) = ssterm(2,1)-dt*k0term(iodd,i,j,k)*(1-idens) &
     &          *thimp*(vel(k2)+vel0(k2))
!
              ddterm(2,1) = ddterm(2,1)+dt*k0term(iodd,i,j,k)*(1-idens) &
     &         *(-thimp*(vel(k2)+vel0(k2))+(0.5*vel(k2)+vel0(k2)))
!
              rrterm(2,1) = rrterm(2,1) - dt*k0term(iodd,i,j,k)         &
     &          *(.5*phi(k2)+phi0(k2))                                  &
     &               - thdts*k0term(iodd,i,j,k)*sphi2(kone)  
!
              ssterm(2,2) = ssterm(2,2)-(1-idens)*dt*k0term(iodd,i,k,j) &
     &          *thimp*(vel(k1)+vel0(k1))
!
              ddterm(2,2) = ddterm(2,2)+(1-idens)*dt*k0term(iodd,i,k,j) &
     &          *(-thimp*(vel(k1)+vel0(k1))+(0.5*vel(k1)+vel0(k1)))

              if(idens.gt.0) then
                ssterm(2,2) = ssterm(2,2)+dent(kone)*k1term(iodd,i,j,k)
                ddterm(2,2) = ddterm(2,2)+dent(kone)*k1term(iodd,i,j,k)
              endif
!
              rrterm(2,2) = rrterm(2,2) + dt*k0term(iodd,i,j,k)         &
     &             *(.5*phi(k1)+phi0(k01))                              &
     &             - thdts*k0term(iodd,i,k,j)*sphi1(kone)              
!
              r40(i2) = r40(i2)                                         &
     &             - thdts*(k0term(iodd,i,j,k)*sphi2(kone)*phi0(j01)    &
     &             + k0term(iodd,i,k,j)*sphi1(kone)*phi0(j2))
              endif
!
!
              if(numvar.ge.3) then
!
              ssterm(1,1) = ssterm(1,1) + dt*(1-idens)*                 &
     &              g8term(iodd,i,j,k)*thimp*(vel(k3)+vel0(k3))
              ddterm(1,1) = ddterm(1,1) + dt*(1-idens)*                 &
     &              g8term(iodd,i,j,k)*(thimp*(vel(k3)+vel0(k3))        &
     &                                - (0.5*vel(k3)+vel0(k3)))
!
              ssterm(1,3) = ssterm(1,3)+dt*g8term(iodd,i,k,j)*(1-idens) &
     &                    *thimp*(vel(k1)+vel0(k1))
!
              ddterm(1,3) = ddterm(1,3)+dt*g8term(iodd,i,k,j)*(1-idens) &
     &            *(thimp*(vel(k1)+vel0(k1))-(0.5*vel(k1)+vel0(k1)))

              if(idens.gt.0) then
!....changed 11/19/06
                ssterm(1,3) = ssterm(1,3)-dent(kone)*k0term(iodd,k,j,i)
                ddterm(1,3) = ddterm(1,3)-dent(kone)*k0term(iodd,k,j,i)
!               ssterm(1,3) = ssterm(1,3)+dent(kone)*k0term(iodd,i,j,k)
!               ddterm(1,3) = ddterm(1,3)+dent(kone)*k0term(iodd,i,j,k)
              endif                  
!
              ssterm(2,2) = ssterm(2,2)-dt*g2term(iodd,j,k,i)*(1-idens) &
     &              *thimp*(vel(k3)+vel0(k3))
!
              ddterm(2,2) = ddterm(2,2)+dt*g2term(iodd,j,k,i)*(1-idens) &
     &             *(-thimp*(vel(k3)+vel0(k3))+(0.5*vel(k3)+vel0(k3)))
!
              ssterm(2,3) = ssterm(2,3)-dt*g2term(iodd,k,j,i)*(1-idens) &
     &              *thimp*(vel(k2)+vel0(k2))
!
              ddterm(2,3) = ddterm(2,3)-dt*g2term(iodd,k,j,i)*(1-idens) &
     &            *(thimp*(vel(k2)+vel0(k2))-(0.5*vel(k2)+vel0(k2)))
!
              ssterm(3,1) = ssterm(3,1) - thsqudts*g7term(iodd,i,j,k)   &
     &                                         *pressure                &
     &         +dt*(1-idens)*g9term(iodd,i,k,j)*thimp*(vel(k3)+vel0(k3)) &
     &         +dt*(1-idens)*(g11erm(iodd,i,k,j)+g11erm(iodd,i,j,k))*    &
     &                thimp*(vel(k1)+vel0(k1))
!
              ddterm(3,1) = ddterm(3,1) - ththmdts*g7term(iodd,i,j,k)   &
     &                                          *pressure               &
     &          +dt*(1-idens)*g9term(iodd,i,k,j)*                       &
     &              (thimp*(vel(k3)+vel0(k3))-(0.5*vel(k3)+vel0(k3)))   &
     &          +dt*(1-idens)*(g11erm(iodd,i,k,j)+g11erm(iodd,i,j,k))*  &
     &              (thimp*(vel(k1)+vel0(k1))-(0.5*vel(k1)+vel0(k1)))
!
              if(idens.gt.0) then
!....changed 11/19/06   scj
                ssterm(3,1) = ssterm(3,1)+dent(kone)*k0term(iodd,k,j,i)
                ddterm(3,1) = ddterm(3,1)+dent(kone)*k0term(iodd,k,j,i)
!               ssterm(3,1) = ssterm(3,1)-dent(kone)*k0term(iodd,i,j,k)
!               ddterm(3,1) = ddterm(3,1)-dent(kone)*k0term(iodd,i,j,k)
              endif 
!
!
              rrterm(3,1) = rrterm(3,1) + dt*(g4term(iodd,i,j,k)        &
     &           +g4term(iodd,i,k,j))*(0.5*phi(k1)+phi0(k01))           &
     &           -thdts*h3term(iodd,i,j,k)*sphi1(kone)
!
              rrterm(3,2) = rrterm(3,2) + dt*(g5term(iodd,i,j,k)        &
     &           +g5term(iodd,i,k,j))*(0.5*phi(k2)+phi0(k2))            &
     &           -thdts*h5term(iodd,i,j,k)*sphi2(kone)
!
              r40(i3) = r40(i3) - thdts*                                &
     &                 (h3term(iodd,i,j,k)*sphi1(kone)*phi0(j01)        &
     &                 +h5term(iodd,i,j,k)*sphi2(kone)*phi0(j2))
!
              ssterm(3,3) = ssterm(3,3) - thsqudts*g6term(iodd,i,j,k)   &
     &                                          *pressure               &
     &         +(1-idens)*dt*g9term(iodd,i,j,k)*thimp*(vel(k1)+vel0(k1)) &
     &         +(1-idens)*dt*(g10erm(iodd,i,k,j)+g10erm(iodd,i,j,k))*   &
     &                thimp*(vel(k3)+vel0(k3))
!
              ddterm(3,3) = ddterm(3,3) - ththmdts*g6term(iodd,i,j,k)   &
     &                                          *pressure               &
     &         +(1-idens)*dt*g9term(iodd,i,j,k)*                        &
     &              (thimp*(vel(k1)+vel0(k1))-(0.5*vel(k1)+vel0(k1)))   &
     &        +(1-idens)*dt*(g10erm(iodd,i,k,j)+g10erm(iodd,i,j,k))*    &
     &             (thimp*(vel(k3)+vel0(k3))-(0.5*vel(k3)+vel0(k3)))
!
              if(idens.gt.0) then
                ssterm(3,3) = ssterm(3,3)+dent(kone)*g2term(iodd,i,j,k)    
                ddterm(3,3) = ddterm(3,3)+dent(kone)*g2term(iodd,i,j,k)
!
!               gravity
!               =======
                ssterm(1,3) = ssterm(1,3) + grav*thsqudts               &
     &               *x2term(iodd,i,j,k)*dent(kone)
                ddterm(1,3) = ddterm(1,3) + grav*ththmdts               &
     &               *x2term(iodd,i,j,k)*dent(kone)
                ssterm(3,1) = ssterm(3,1) - grav*thsqudts               &
     &               *y1term(iodd,i,j,k)*dent(kone)
                ddterm(3,1) = ddterm(3,1) - grav*ththmdts               &
     &               *y1term(iodd,i,j,k)*dent(kone)
                ssterm(3,3) = ssterm(3,3) - grav*thsqudts               &
     &               *y2term(iodd,i,j,k)*dent(kone)
                ddterm(3,3) = ddterm(3,3) - grav*ththmdts               &
     &               *y2term(iodd,i,j,k)*dent(kone)
!               =======
             endif
             endif

! QUADRATIC NONLINEAR TERMS
          do l=1,18
             l1 = isvaln(itri,l)
             l01= isval0(itri,l)
             l2 = l1 + 6
             l3 = l2 + 6
             lone = isval1(itri,l)
!
             pmix = c1term2(k,l)*                                       &
     &            (phi0(l01)+phi(l1))*(phi0(k01)+phi(k1))
!
             ssterm(1,1) = ssterm(1,1) + thsqudts*pmix   
             ddterm(1,1) = ddterm(1,1) + ththmdts*pmix
!
             if(idens.gt.0) then
                ssterm(1,1) = ssterm(1,1) + dt*thimp*dent(lone)         &
     &               *(v1terml(k,l)+v2terml(k,l))*(vel(k1)+vel0(k1))
!
                ddterm(1,1) = ddterm(1,1) + dt*dent(lone)               &
     &            *(v1terml(k,l)+v2terml(k,l))*(thimp*(vel(k1)+vel0(k1))&
     &                       -(0.5*vel(k1)+vel0(k1)))
             endif
!   
!            gyroviscous contribution:
!            -------------------------
             if(gyro.eq.1) then
!               (1, 1)
                tempvar = gv11arml(k,l)*(galph(lone)+gw1(lone))
                if(numvar.ge.2) then
                   tempvar = tempvar*(phi(k2)+phi0(k2))
                else 
                   tempvar = tempvar*bzero/18.
                endif
                tempvar = tempvar
                ssterm(1,1) = ssterm(1,1) + dt* thimp    *tempvar
                ddterm(1,1) = ddterm(1,1) + dt*(thimp-1.)*tempvar
             endif
!            -------------------------
!
             if(numvar.ge.2) then
!
                qmix = c7term2(k,l)*                                    &
     &               (phi(k1)+phi0(k01))*(phi(l2)+phi0(l2))
! 
                qmix2 =   c9term2(k,l)                                  &
     &               *(phi(k1)+phi0(k01))*(phi(l1)+phi0(l01))
!
                ssterm(2,1) = ssterm(2,1) - thsqudts*qmix 
                ddterm(2,1) = ddterm(2,1) - ththmdts*qmix
!
                ssterm(2,2) = ssterm(2,2) - thsqudts*qmix2
                ddterm(2,2) = ddterm(2,2) - ththmdts*qmix2

                if(idens.gt.0) then
                   ssterm(2,1) = ssterm(2,1) + dt*dent(lone)            &
     &                  *v8terml(k,l)  *thimp*(vel(k2)+vel0(k2))
                   ddterm(2,1) = ddterm(2,1) + dt*dent(lone)            &
     &                  *v8terml(k,l)*(thimp*(vel(k2)+vel0(k2))         &
     &                  - (0.5*vel(k2)+vel0(k2)))
!
                   ssterm(2,2) = ssterm(2,2) + dt*dent(lone)            &
     &                  *v7terml(k,l)  *thimp*(vel(k1)+vel0(k1))
                   ddterm(2,2) = ddterm(2,2) + dt*dent(lone)            &
     &                  *v7terml(k,l)*(thimp*(vel(k1)+vel0(k1))         &
     &                  - (0.5*vel(k1)+vel0(k1)))
                endif           ! on idens
!
!             gyroviscous contribution:
!             -------------------------
              if(gyro.eq.1) then
!               (1,2)
                 tempvar = (phi(k1)+phi0(k01))                          &
     &                *(gv12arml(k,l)*(galph(lone)-gmu(lone))           &
     &                 +gv12brml(k,l)*gw2(lone)                         &
     &                 +gv12crml(k,l)*gw3(lone))
                 tempvar = tempvar
                 ssterm(1,2) = ssterm(1,2) + dt* thimp    *tempvar
                 ddterm(1,2) = ddterm(1,2) + dt*(thimp-1.)*tempvar
!     
!               (2,1)
                 tempvar = (phi(k1)+phi0(k01))                          &
     &                *(gv12axml(k,l)*(galph(lone)-gmu(lone))           &
     &                 +gv12bxml(k,l)*gw2(lone)                         &
     &                 +gv12cxml(k,l)*gw3(lone))
                 tempvar = tempvar
                 ssterm(2,1) = ssterm(2,1) + dt* thimp    *tempvar
                 ddterm(2,1) = ddterm(2,1) + dt*(thimp-1.)*tempvar
!     
!               (2,2)
                 tempvar = gv22arml(k,l)*(phi(k2)+phi0(k2))             &
     &                *(galph(lone)+gmu(lone)                           &
     &                 - 2.*gw1(lone))    
                 tempvar = tempvar
                 ssterm(2,2) = ssterm(2,2) + dt* thimp    *tempvar
                 ddterm(2,2) = ddterm(2,2) + dt*(thimp-1.)*tempvar
              endif
!             -------------------------
!
            endif ! on numvar

            if(numvar.ge.3) then
               if(idens.gt.0) then
                  ssterm(1,1) = ssterm(1,1) + dt*thimp*dent(lone)       &
     &                 *v3terml(k,l)*(vel(k3)+vel0(k3))
!     
                  ddterm(1,1) = ddterm(1,1) + dt*dent(lone)             &
     &                 *v3terml(k,l)*(thimp*(vel(k3)+vel0(k3))          &
     &                 -(0.5*vel(k3)+vel0(k3)))
               endif            ! on idens .gt. 0
!
               c6ave =  c6term2(k,l)
               ssterm(1,3) = ssterm(1,3) + thsqudts*                    &
     &              c6ave*(phi0(l01)+phi(l1))*(phi0(k01)+phi(k1))
!     
               ddterm(1,3) = ddterm(1,3) + ththmdts*                    &
     &              c6ave*(phi0(l01)+phi(l1))*(phi0(k01)+phi(k1))
!     
               if(idens.gt.0) then
                  ssterm(1,3) = ssterm(1,3) + dt*thimp*dent(lone)       &
     &                 *( (v5terml(k,l)+v6terml(k,l))*(vel(k3)+vel0(k3)) &
     &                 +v4terml(k,l)*(vel(k1)+vel0(k1)))
                  
                  ddterm(1,3) = ddterm(1,3) + dt*dent(lone)             &
     &          *( (v5terml(k,l)+v6terml(k,l))*(thimp*(vel(k3)+vel0(k3))&
     &                       -(0.5*vel(k3)+vel0(k3)))                   &
     &             +v4terml(k,l)*(thimp*(vel(k1)+vel0(k1))              &
     &                       -(0.5*vel(k1)+vel0(k1))))
!
                  ssterm(2,2) = ssterm(2,2) + dt*dent(lone)             &
     &                 *v9terml(k,l)  *thimp*(vel(k3)+vel0(k3))
                  ddterm(2,2) = ddterm(2,2) + dt*dent(lone)             &
     &                 *v9terml(k,l)*(thimp*(vel(k3)+vel0(k3))          &
     &                 - (0.5*vel(k3)+vel0(k3)))
               endif            ! on idens .gt. 0
!
!
               ssterm(2,3) = ssterm(2,3) - thsqudts*                    &
     &              c8term2(k,l)*(phi0(k01)+phi(k1))*(phi0(l2)+phi(l2))
!
               ddterm(2,3) = ddterm(2,3) - ththmdts*                    &
     &              c8term2(k,l)*(phi0(k01)+phi(k1))*(phi0(l2)+phi(l2))
!
               if(idens.gt.0) then
                  ssterm(2,3) = ssterm(2,3) + dt*thimp*dent(lone)       &
     &                 *v10erml(k,l)*(vel(k2)+vel0(k2))
                  
                  ddterm(2,3) = ddterm(2,3) + dt*dent(lone)             &
     &                 *v10erml(k,l)*(thimp*(vel(k2)+vel0(k2))          &
     &                 -(0.5*vel(k2)+vel0(k2)))                         &
               endif            ! on idens .gt. 0
!
               ssterm(3,1) = ssterm(3,1) + thsqudts*(                   &
     &        c10erm2(k,l)*(phi0(l01)+phi(l1))*(phi0(k01)+phi(k1))      &
c$$$     &  +     c11erm2(k,l)*(phi0(l01)+phi(l1))*(jphi0(kone)+jphi(kone)) &
     &       +c14erml(k,l)*(phi0(l2)+phi(l2))*(phi0(k2)+phi(k2)))
!
               ddterm(3,1) = ddterm(3,1) + ththmdts*(                   &
     &        c10erm2(k,l)*(phi0(l01)+phi(l1))*(phi0(k01)+phi(k1))      &
c$$$     &  +     c11erm2(k,l)*(phi0(l01)+phi(l1))*(jphi0(kone)+jphi(kone)) &
     &       +c14erml(k,l)*(phi0(l2)+phi(l2))*(phi0(k2)+phi(k2)))
!
               if(idens.gt.0) then
                  ssterm(3,1) = ssterm(3,1) + dt*thimp*dent(lone)       &
     &                *( (v15erml(k,l)+v16erml(k,l))*(vel(k1)+vel0(k1)) &
     &                +v12erml(k,l)*(vel(k3)+vel0(k3)))
                 
                  ddterm(3,1) = ddterm(3,1) + dt*dent(lone)             &
     &          *( (v15erml(k,l)+v16erml(k,l))*(thimp*(vel(k1)+vel0(k1))&
     &                -(0.5*vel(k1)+vel0(k1)))                          &
     &                +v12erml(k,l)*(thimp*(vel(k3)+vel0(k3))           &
     &                -(0.5*vel(k3)+vel0(k3))))
               endif            ! on idens.gt.0
!
               ssterm(3,2) = ssterm(3,2) + thsqudts*c16erml(k,l)        &
     &              *(phi(k1)+phi0(k01))*(phi(l2)+phi0(l2))

               ddterm(3,2) = ddterm(3,2) + ththmdts*c16erml(k,l)        &
     &              *(phi(k1)+phi0(k01))*(phi(l2)+phi0(l2))
!
               ssterm(3,3) = ssterm(3,3) + thsqudts*                    &
     &        (c13erml(k,l)*(phi0(l01)+phi(l1))*(phi0(k01)+phi(k1))     &
c$$$     &        +c12erml(k,l)*(phi0(k01)+phi(k1))*(jphi0(lone)+jphi(lone))&
     &        +c15erml(k,l)*(phi0(l2)+phi(l2))*(phi0(k2)+phi(k2)))
!
               ddterm(3,3) = ddterm(3,3) + ththmdts*                    &
     &        (c13erml(k,l)*(phi0(l01)+phi(l1))*(phi0(k01)+phi(k1))     &
c$$$     &        +c12erml(k,l)*(phi0(k01)+phi(k1))*(jphi0(lone)+jphi(lone))&
     &        +c15erml(k,l)*(phi0(l2)+phi(l2))*(phi0(k2)+phi(k2)))
!
!......semi-implicit terms added 03/12/06
               if(isemii.gt.0) then
                  ssterm(3,3) = ssterm(3,3) + thsqudts*                 &
     &                 c15erml(k,l)*bsi(lone)*bsi(kone)
                  ddterm(3,3) = ddterm(3,3) + thsqudts*                 &
     &                 c15erml(k,l)*bsi(lone)*bsi(kone)
               endif            ! on isemii .gt. 0
!
!              gyroviscous contribution:
!              -------------------------
               if(gyro.eq.1) then
!                (1,3)
                  tempvar = (phi(k2)+phi0(k2))                          &
     &                 *(gv13arml(k,l)*(galph(lone)+gw1(lone))          &
     &                  +gv13brml(k,l)*gw2(lone)                        &
     &                  +gv13crml(k,l)*gw3(lone))
                  tempvar = tempvar
                  ssterm(1,3) = ssterm(1,3) + dt* thimp    *tempvar
                  ddterm(1,3) = ddterm(1,3) + dt*(thimp-1.)*tempvar
!
!                (2,3)
                  tempvar = (phi(k1)+phi0(k01))                         &
     &                 *(gv23arml(k,l)*galph(lone)                      &
     &                  +gv23brml(k,l)*gmu(lone)                        &
     &                  -gv12cxml(k,l)*gw2(lone)                        &
     &                  +gv12bxml(k,l)*gw3(lone))
                  tempvar = tempvar
                  ssterm(2,3) = ssterm(2,3) + dt* thimp    *tempvar
                  ddterm(2,3) = ddterm(2,3) + dt*(thimp-1.)*tempvar
!
!                (3,1)
                  tempvar = -(phi(k2)+phi0(k2))                         &
     &                 *(gv13arml(k,l)*(galph(lone)+gw1(lone))          &
     &                  +gv13bxml(k,l)*gw2(lone)                        &
     &                  +gv13cxml(k,l)*gw3(lone))
                  tempvar = tempvar
                  ssterm(3,1) = ssterm(3,1) + dt* thimp    *tempvar
                  ddterm(3,1) = ddterm(3,1) + dt*(thimp-1.)*tempvar
!
!                (3,2)
                  tempvar = (phi(k1)+phi0(k01))                         &
     &                 *(gv23axml(k,l)*galph(lone)                      &
     &                  +gv23bxml(k,l)*gmu(lone)                        &
     &                  -gv12crml(k,l)*gw2(lone)                        &
     &                  +gv12brml(k,l)*gw3(lone))
                  tempvar = tempvar
                  ssterm(3,2) = ssterm(3,2) + dt* thimp    *tempvar
                  ddterm(3,2) = ddterm(3,2) + dt*(thimp-1.)*tempvar
!
!                (3,3)
                  tempvar = (phi(k2)+phi0(k2))                          &
     &                 *(gv11arml(k,l)*(galph(lone)+gw1(lone))          &
     &                  +gv33brml(k,l)*gw2(lone)                        &
     &                  +gv33crml(k,l)*gw3(lone))
                  tempvar = tempvar
                  ssterm(3,3) = ssterm(3,3) + dt* thimp    *tempvar
                  ddterm(3,3) = ddterm(3,3) + dt*(thimp-1.)*tempvar
               endif
!              -------------------------
!
               if(idens.gt.0) then
                  ssterm(3,3) = ssterm(3,3) + dt*thimp*dent(lone)       &
     &                 *( (v13erml(k,l)+v14erml(k,l))*(vel(k3)+vel0(k3)) &
     &                 +v11erml(k,l)*(vel(k1)+vel0(k1)))
                  
                  ddterm(3,3) = ddterm(3,3) + dt*dent(lone)             &
     &          *( (v13erml(k,l)+v14erml(k,l))*(thimp*(vel(k3)+vel0(k3))&
     &                 -(0.5*vel(k3)+vel0(k3)))                         &
     &                 +v11erml(k,l)*(thimp*(vel(k1)+vel0(k1))          &
     &                 -(0.5*vel(k1)+vel0(k1))))
             endif              ! on idens.gt.0
!
          endif                 ! on numvar.ge.3
       enddo                    ! on l
      enddo                     ! on k

      call sparseR8d_increment_array(s1matrix_lu,i1,j1,ssterm,numvar)
      call sparseR8d_increment_array(d1matrix_lu,i1,j1,ddterm,numvar)
      call sparseR8d_increment_array(r1matrix_lu,i1,j1,rrterm,numvar)
 800  continue
      enddo                     ! on irect
      enddo                     ! on jrect
      if(myrank.eq.0 .and. itimer.gt.0) call second(tafterelements)
      telements = telements + tafterelements-tafterread
      enddo                     ! on j
      enddo                     ! on i
      enddo                     ! on iodd

#ifdef mpi 
      icount = 6*numvar*m*n
      vtemp = r40
      call MPI_ALLREDUCE(vtemp,r40,icount,MPI_DOUBLE_PRECISION,         &
     &     MPI_SUM, MPI_COMM_WORLD,ier)
      if(ier /= 0) then
         print *,'Error with MPI_allreduce-r40',ier
         call safestop(1)
      endif
#endif
      if(myrank.eq.0) call second(tbeforeinit)
      tmpi = tmpi + tbeforeinit-tafterelements
!
!.....modify the s-matrix, inserting the boundary conditions
!
!.....calculate boundary conditions
      call boundaryv(iboundv,nbcv)
      if(nbcv .gt. iboundmax) then
         write(*,4888) nbcv, iboundmax
 4888    format(" ERROR: nbcv > iboundmax", 2i5)
         call safestop(3)
      endif

      do l=1,nbcv

         i = iboundv(l)

         if(sparseR8d_is_local_row(s1matrix_lu, i).eq.1) then
            if(i.le.0 .or. i.ge.m*n*6*numvar+1) then
               write(*,3341) i,l,iboundv(l)
 3341          format(" error in iboundv index",3i5)
               call safestop(3341)
            endif
            call sparseR8d_zero_row(s1matrix_lu,i)
            call sparseR8d_increment(s1matrix_lu,i,i,1.,numvar,1)
         endif
       enddo
    
      jer = 0
      if(iprint.ge.1) write(*,*) "before save s1matrix" 
      if(myrank.eq.0) call second(tafterinit)
      tinit = tinit + tafterinit-tbeforeinit

      if(ifirsts1_lu.eq.0) then
         ifirsts1_lu = 1
         call set_superlu_options(s1matrix_lu%options, Fact=0)
      else
         call set_superlu_options(s1matrix_lu%options, Fact=2)
      endif

      return
      end
!============================================================
      subroutine ludefphi
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
#ifdef mpi
      use supralu_dist_mod
      use superlu_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer i, i1, i2, i3
      integer j, jone, j01, j1, j2, j3
      integer k, kone, k01, k1, k2, k3
      integer l, lone, l1, l2, l3
      integer itype, iodd, itri, ll, irect, jrect, icount
      integer ir1, ir2, ir3, ier
      real ssterm(3,3),ddterm(3,3),rrterm(3,3),qqterm(3,3)
      real gmix, terma, termb, termd, u12terml, u78terml, termbf
      real t24terml, t35terml, temp
      real x, z, dbf, factor, hypf, hypi, hypp, denf
      real tstartludefphi, tendludefphi, tbeforeread, tafterread
      real tafterelements, tbeforeinit, tafterinit
      real o2, chif

      ! this factor multiplies pe terms
      o2 = pefac
!
!     implicit parameter for hyperviscosity terms
      if(myrank.eq.0 .and. itimer.gt.0) call second(tstartludefphi)
!
c$$$      if(ifirsts2_lu.eq.1) then
c$$$         call f_Destroy_LU(s2matrix_lu%n, SLUprocgrid,                     &
c$$$     &        s2matrix_lu%LUstruct)
c$$$      endif
      call sparseR8d_reset(s2matrix_lu, ier)
      call sparseR8d_reset(d2matrix_lu, ier)
      call sparseR8d_reset(r2matrix_lu, ier)
      call sparseR8d_reset(q2matrix_lu, ier)
      s2matrix_lu%colind = origcolind
      d2matrix_lu%colind = origcolind
      r2matrix_lu%colind = origcolind
      q2matrix_lu%colind = origcolind

!
!.....form field matrices
      ohmic = 0
      viscous = 0.
      q34 = 0.
      itri = 0

!.......position files with metric terms
      do itype=1,ntensor
         if(iread(itype).eq.1) rewind(itype+100)
      enddo
!
      do iodd=1,ioddmx
      do i=1,18
      do j=1,18
!
        if(myrank.eq.0 .and. itimer.gt.0) call second(tbeforeread)
!........read the 4th rank tensor terms from the disk
        if(iread(01).eq.1) then
           read(101,err=1001,end=1001)                                  &
     &          ir1,ir2,ir3,((u3terml(k,l),l=1,18),k=1,18)
           if(ir1.ne.iodd .or. ir2.ne.i .or. ir3.ne.j) go to 1001
        endif
!     
        if(iread(02).eq.1)                                              &
     &       read(102) ir1,ir2,ir3,((u4terml(k,l),l=1,18),k=1,18)
        if(iread(03).eq.1)                                              &
     &       read(103) ir1,ir2,ir3,((u1terml(k,l),l=1,18),k=1,18)
        if(iread(04).eq.1)                                              &
     &       read(104) ir1,ir2,ir3,((u2terml(k,l),l=1,18),k=1,18)
        if(iread(09).eq.1)                                              &
     &       read(109) ir1,ir2,ir3,((u5terml(k,l),l=1,18),k=1,18)
        if(iread(10).eq.1)                                              &
     &       read(110) ir1,ir2,ir3,((u6terml(k,l),l=1,18),k=1,18)
        if(iread(47).eq.1)                                              &
     &       read(147) ir1,ir2,ir3,((u7terml(k,l),l=1,18),k=1,18)
        if(iread(48).eq.1)                                              &
     &       read(148) ir1,ir2,ir3,((u8terml(k,l),l=1,18),k=1,18)
        if(iread(73).eq.1)                                              &
     &       read(173) ir1,ir2,ir3,((t1terml(k,l),l=1,18),k=1,18)
        if(iread(74).eq.1)                                              &
     &       read(174) ir1,ir2,ir3,((t2terml(k,l),l=1,18),k=1,18)
        if(iread(75).eq.1)                                              &
     &       read(175) ir1,ir2,ir3,((t3terml(k,l),l=1,18),k=1,18)
        if(iread(76).eq.1)                                              &
     &       read(176) ir1,ir2,ir3,((t4terml(k,l),l=1,18),k=1,18)
        if(iread(77).eq.1)                                              &
     &       read(177) ir1,ir2,ir3,((t5terml(k,l),l=1,18),k=1,18)
!     
        go to 1002
 1001   continue
        write(*,*) ir1,ir2,ir3,iodd,i,j
        call safestop(1001)
 1002   continue
        if(myrank.eq.0 .and. itimer.gt.0) call second(tafterread)
        tread = tread + tafterread-tbeforeread

        do jrect=1,m-1+jper
        do irect=1,n-1+iper
!
           ll = irect + (jrect-1)*(n-1+iper)
           itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll

           i1 = isvaln(itri,i)

           ! check whether this node is assigned to this processor
           if(sparseR8d_is_local_row(s2matrix_lu, i1).eq.0) goto 800

           j1 = isvaln(itri,j)
           j01 = isval0(itri,j)
           j2 = j1 + 6
           i2 = i1 + 6
           j3 = j2 + 6
           i3 = i2 + 6
           jone = isval1(itri,j)

           x = (irect-iodd/3.)*deex
           z = (jrect-1+iodd/3.)*deez
           factor = 1
           if(imask.eq.1) call mask(x,z,factor)
!
           ssterm = 0
           ddterm = 0
           rrterm = 0
           qqterm = 0
!
           terma = aterm(iodd,i,j)
           termb = bterm(iodd,i,j)
           termd = dterm(iodd,i,j)
           hypf = hyper *deex**2
           hypi = hyperi *deex**2
           hypp = hyperp *deex**2
           termbf = hypf*termb
!     
           dbf = db*factor
           chif = (gam-1.)*chipar*factor
!
!
!...NOTE: s2matrix_lu (formed from ss) multiplies the solution at time NTIME+1
!
! LINEAR TERMS
            ssterm(1,1) =  termd -     thimp  *dt*etar*(terma-termbf)  
            ddterm(1,1) =  termd + (1.-thimp) *dt*etar*(terma-termbf)   
!
            if(linear.eq.0) then
               ohmic(i1) = ohmic(i1) + etar*(terma-termbf)*phi0(j01)
            endif               ! on linear.eq.0
!
            if(numvar.ge.2) then
!
              ssterm(2,2) =  termd   - thimp *dt*etar*(terma-hypi*termb) 
              ddterm(2,2) =  termd+(1.-thimp)*dt*etar*(terma-hypi*termb) 
              if(linear.eq.0) then
                 ohmic(i2) = ohmic(i2)+etar*(terma-hypi*termb)*phi0(j2)
              endif            ! on linear.eq.0
!
              if(numvar.ge.3) then
!     
                 ssterm(3,3) =  termd -     thimp *dt*kappa             &
     &                                            *(terma-hypp*termb)
                 ddterm(3,3) =  termd + (1.-thimp)*dt*kappa             &
     &                                            *(terma-hypp*termb)
!
                 if(linear.eq.0) then
                    ohmic(i3) = ohmic(i3) + kappa*(terma-hypp*termb)    &
     &                                            *phi0(j3)
                 endif          ! on linear.eq.0
!
              endif ! numvar.ge.3
            endif ! numvar.ge.2
!
!
!
!
! SIMPLE NONLINEAR TERMS
            do k=1,18
              k1= isvaln(itri,k)
              k01 = isval0(itri,k)
              k2= k1 + 6
              k3= k2 + 6
              kone = isval1(itri,k)

              gmix = g0term(iodd,i,j,k)+g0term(iodd,i,k,j)
!
              ssterm(1,1)=ssterm(1,1) - dt*thimp                        &
     &           *(k0term(iodd,i,k,j)*(veln(k1)+vel0(k1)))
!
              ddterm(1,1)=ddterm(1,1) + dt*k0term(iodd,i,k,j)*          &
     &          (-thimp*(veln(k1)+vel0(k1))+0.5*veln(k1)+vel0(k1))
!
              rrterm(1,1)=rrterm(1,1) + thimp*dt*k0term(iodd,i,j,k)*    &
     &                      (phi(k1)+phi0(k01))
!
              qqterm(1,1)=qqterm(1,1) +       dt*k0term(iodd,i,j,k)*    &
     &          (-thimp*(phi(k1)+phi0(k01))+0.5*phi(k1)+phi0(k01))
!
              if(numvar.ge.2) then
                 if(idens.eq.0) then
!
                 ssterm(1,1)=ssterm(1,1) - dbf*dt*k0term(iodd,i,j,k)    &
     &                *thimp*(phi(k2)+phi0(k2))
!     
                 ddterm(1,1)=ddterm(1,1) - dbf*dt*k0term(iodd,i,j,k)    &
     &                *(thimp*(phi(k2)+phi0(k2))-(0.5*phi(k2)+phi0(k2)))
!
                 ssterm(1,2) = ssterm(1,2) - dt*dbf*k0term(iodd,i,k,j)  &
     &                *thimp*(phi(k1)+phi0(k01))
                 ddterm(1,2) = ddterm(1,2) -   dt*dbf*k0term(iodd,i,k,j)&
     &              *(thimp*(phi(k1)+phi0(k01))-(0.5*phi(k1)+phi0(k01)))
!
                 ssterm(2,1) = ssterm(2,1) -dt*thimp*                   &
     &                dbf*gmix*(phi(k1)+phi0(k01))
                 ddterm(2,1) = ddterm(2,1) - dt*                        &
     &                               dbf*gmix*(thimp*(phi(k1)+phi0(k01))&  
     &                                          -(.5*phi(k1)+phi0(k01)))
                 endif ! on idens.eq.0
!
!
                 ssterm(2,1) = ssterm(2,1)                              &
     &                -dt*k0term(iodd,i,k,j)*thimp*(veln(k2)+vel0(k2)) 
!
                 ddterm(2,1) = ddterm(2,1)                              &
     &                -dt*k0term(iodd,i,k,j)*(thimp*(veln(k2)+vel0(k2)) &
     &                                    -(0.5*veln(k2)+vel0(k2)) )
!
                 rrterm(2,1)=rrterm(2,1) + thimp*dt*k0term(iodd,i,j,k)* &
     &                (phi(k2)+phi0(k2))
!
                 qqterm(2,1)=qqterm(2,1) + dt*k0term(iodd,i,j,k)*       &
     &                (-thimp*(phi(k2)+phi0(k2))+0.5*phi(k2)+phi0(k2))
!
!
                 ssterm(2,2)=ssterm(2,2) - dt*k0term(iodd,i,k,j)*       &
     &                thimp*(veln(k1)+vel0(k1))
!
                 ddterm(2,2)=ddterm(2,2) + dt*k0term(iodd,i,k,j)*       &
     &                (-thimp*(veln(k1)+vel0(k1))+0.5*veln(k1)+vel0(k1))
!
                 rrterm(2,2) = rrterm(2,2) + dt*thimp*                  &
     &                k0term(iodd,i,j,k)*(phi(k1)+phi0(k01))
!
                 qqterm(2,2) = qqterm(2,2) + dt*k0term(iodd,i,j,k)*     &
     &                (-thimp*(phi(k1)+phi0(k01))+0.5*phi(k1)+phi0(k01))
!
                 if(numvar.ge.3) then
!
                    ssterm(1,1)=ssterm(1,1) - dt*thimp*                 &
     &                   g2term(iodd,j,k,i)*(veln(k3)+vel0(k3))
!
                    ddterm(1,1)=ddterm(1,1)                             &
     &                 -dt*g2term(iodd,j,k,i)*(thimp*(veln(k3)+vel0(k3)) &
     &                               -(0.5*veln(k3)+vel0(k3)))
!
!
                 rrterm(1,3)=rrterm(1,3) + thimp*dt*g2term(iodd,k,j,i)* &
     &                      (phi(k1)+phi0(k01))
!     
                    qqterm(1,3)=qqterm(1,3) - dt*g2term(iodd,k,j,i)*    &
     &               (thimp*(phi(k1)+phi0(k01))-(0.5*phi(k1)+phi0(k01)))
!
!
                    ssterm(2,2)=ssterm(2,2) + dt*thimp*                 &
     &                   g2term(iodd,i,k,j)*(veln(k3)+vel0(k3))
!
                    ddterm(2,2)=ddterm(2,2) +                           &
     &                  dt*g2term(iodd,i,k,j)*(thimp*(veln(k3)+vel0(k3)) &
     &                  -(0.5*veln(k3)+vel0(k3)))
!
!
                    rrterm(2,3)=rrterm(2,3)- thimp*dt*g2term(iodd,i,j,k) &
     &                   *(phi(k2)+phi0(k2))
!
                    qqterm(2,3)=qqterm(2,3) + dt*g2term(iodd,i,j,k)*    &
     &            (thimp*(phi(k2)+phi0(k2))-(0.5*phi(k2)+phi0(k2)))
!
                   rrterm(3,1)=rrterm(3,1) - dt*thimp*k0term(iodd,i,k,j) &
     &           *(phi(k3)+phi0(k3))
!
                    qqterm(3,1)=qqterm(3,1) + dt*k0term(iodd,i,k,j)*    &
     &                 (thimp*(phi(k3)+phi0(k3))-(0.5*phi(k3)+phi0(k3)))
!
                    ssterm(3,3) = ssterm(3,3) + thimp*dt*               &
     &                   (k0term(iodd,i,j,k)*(veln(k1)+vel0(k1))         &
     &                  +(g2term(iodd,i,k,j)-(gam-1)*k2term(iodd,i,k,j)) &
     &                   *(veln(k3)+vel0(k3))                            &
     &          - (gam-1.)*kappat*(g2term(iodd,i,k,j)+g2term(iodd,i,j,k)  &
     &                           - hypp*g13erm(iodd,i,k,j))*deni(kone))
!
                    ddterm(3,3) = ddterm(3,3) +       dt*               &
     &                  (k0term(iodd,i,j,k)*(thimp*(  veln(k1)+vel0(k1)) &
     &                   -(0.5*veln(k1)+vel0(k1)))                       &
     &                  +(g2term(iodd,i,k,j)-(gam-1)*k2term(iodd,i,k,j)) &
     &              *(thimp*(veln(k3)+vel0(k3))-(0.5*veln(k3)+vel0(k3))) &
     &          + (gam-1.)*kappat*(g2term(iodd,i,k,j)+g2term(iodd,i,j,k) &
     &                - hypp*g13erm(iodd,i,k,j))*(1.-thimp)*deni(kone))
!
              if(linear.eq.0) then
                 ohmic(i3) = ohmic(i3) + (gam-1.)*kappat                  &
     &          *(g2term(iodd,i,k,j)+g2term(iodd,i,j,k)                  &
     &                  - hypp*g13erm(iodd,i,k,j))*deni(kone)*phi0(j3)
              endif ! on linear.eq.0
!
!      Ohmic heating term  (added 1/12/06) Viscous added 1/19/06
              if(linear.eq.1) then
              ohmic(i3) = ohmic(i3) + (gam-1.)*etar*                      &
     &           (k1term(iodd,i,j,k)*((jphi(kone)+jphi0(kone))           &
     &                              *(jphi(jone)+jphi0(jone))            &
     &                               -jphi0(kone)*jphi0(jone))           &
     &        -g2term(iodd,k,j,i)*((phi(k2)+phi0(k2))*(phi(j2)+phi0(j2)) &
     &                                -phi0(k2)*phi0(j2))                &
     &        -hypf*g2term(iodd,k,j,i)*((jphi(kone)+jphi0(kone))         &
     &                                 *(jphi(jone)+jphi0(jone))         &
     &                                      -jphi0(kone)*jphi0(jone))    &
     &        +hypi*g12erm(iodd,i,j,k)*((phi(k2)+phi0(k2))               &
     &                                 *(phi(j2)+phi0(j2))               &
     &                                 -phi0(k2)*phi0(j2)))
              else ! on linear
              ohmic(i3) = ohmic(i3) + (gam-1.)*etar*(                     &
     &            k1term(iodd,i,j,k)*(jphi(kone)+jphi0(kone))            &
     &                              *(jphi(jone)+jphi0(jone))            &
     &        -g2term(iodd,k,j,i)*(phi(k2)+phi0(k2))*(phi(j2)+phi0(j2))  &
     &        -hypf*g2term(iodd,k,j,i)*(jphi(kone)+jphi0(kone))          &
     &                                 *(jphi(jone)+jphi0(jone))         &
     &        +hypi*g12erm(iodd,i,j,k)*(phi(k2)+phi0(k2))                &
     &                                *(phi(j2)+phi0(j2)))
              endif ! on linear
!
!.....changed def of viscous term on 08/08/06
!             viscous(i3) = viscous(i3) + (gam-1.)*amu*                  &
!    &             (k1term(iodd,i,j,k)* vor(kone)*vor(jone)              &
!    &             - g2term(iodd,k,j,i)*veln(k2)*veln(j2))                 &
!    &         + (gam-1.)*amuc*k1term(iodd,i,j,k)*com(kone)*com(jone)
!
              viscous(i3) = viscous(i3) + (gam-1.)*amu*                  &
     &         (  g14erm(iodd,i,j,k)*veln(j1)*veln(k1)                     &
     &           -g2term(iodd,k,j,i)*veln(k2)*veln(j2))
              viscous(i3) = viscous(i3) + (gam-1.)*.5*(amu+amuc)*        &
     &            g16erm(iodd,i,j,k)*veln(j1)*veln(k3)
              viscous(i3) =viscous(i3) + (gam-1.)*amuc*                  &
     &         (  g15erm(iodd,i,j,k)*veln(j3)*veln(k3))

              if(idens.eq.0) then
!
              ssterm(3,2)=ssterm(3,2) - dt*dbf*k0term(iodd,i,k,j)        &
     &           *thimp*(phi(k3)+phi0(k3))
!
              ddterm(3,2)=ddterm(3,2) - dt*dbf*k0term(iodd,i,k,j)*       &
     &           (thimp*(phi(k3)+phi0(k3))-(0.5*phi(k3)+phi0(k3)))
!
              ssterm(3,3) = ssterm(3,3)                                  &
     &            -dt*dbf*k0term(iodd,i,j,k)*thimp*(phi(k2)+phi0(k2))    &
!
              ddterm(3,3) = ddterm(3,3)                                  &
     &            -dt*dbf*k0term(iodd,i,j,k)*(thimp*(phi(k2)+phi0(k2))  &
     &                                 -(0.5*phi(k2)+phi0(k2)))          &
!
              else ! on idens.eq.0
                 ssterm(2,3)=ssterm(2,3) - dt*dbf*k0term(iodd,i,k,j)    &
     &                *thimp *deni(kone)
!
                 ddterm(2,3)=ddterm(2,3) + dt*dbf*k0term(iodd,i,k,j)    &
     &                *(1-thimp)*deni(kone)
!
!NATE: Why is the following commented out?
!                ssterm(3,3) = ssterm(3,3) -    thimp*dt*kappa*deni(kone) &
!     &                *(g2term(iodd,i,j,k)+g2term(iodd,i,k,j))
!                ddterm(3,3) = ddterm(3,3) +(1-thimp)*dt*kappa*deni(kone) &
!     &                *(g2term(iodd,i,j,k)+g2term(iodd,i,k,j))
!
              endif ! on idens
!
              rrterm(3,3) = rrterm(3,3) -dt*thimp*                      &
     &                  (g2term(iodd,i,j,k)-(gam-1)*k2term(iodd,i,j,k)) &
     &                  *(phi(k3)+phi0(k3))
!
              qqterm(3,3) = qqterm(3,3) +dt*                            &
     &                  (g2term(iodd,i,j,k)-(gam-1)*k2term(iodd,i,j,k)) &
     &          *(thimp*(phi(k3)+phi0(k3))-(0.5*phi(k3)+phi0(k3)))
!
           endif                ! on numvar.ge.3
!
        endif                   ! on numvar.ge.2
!

! QUADRATIC NONLINEAR TERMS
        do l=1,18
              l1 = isvaln(itri,l)
              l2 = l1 + 6
              l3 = l2 + 6
              lone = isval1(itri,l)

              if(idens.gt.0) then
              if(numvar.ge.2) then
                denf = deni(lone)*dbf
                ssterm(1,1) = ssterm(1,1) - dt*u3terml(k,l)*denf        &
     &                *thimp*(phi(k2)+phi0(k2))
                ddterm(1,1) = ddterm(1,1) - dt*u3terml(k,l)*denf        &
     &               *(thimp*(phi(k2)+phi0(k2))-(0.5*phi(k2)+phi0(k2)))
!
                ssterm(1,2) = ssterm(1,2) - dt*u4terml(k,l)*denf        &
     &                *thimp*(phi(k1)+phi0(k01))
                ddterm(1,2) = ddterm(1,2) - dt*u4terml(k,l)*denf        &
     &              *(thimp*(phi(k1)+phi0(k01))-(0.5*phi(k1)+phi0(k01)))
!
                u12terml = u1terml(k,l)+u2terml(k,l)
                u78terml = u7terml(k,l)+u8terml(k,l)
                ssterm(2,1) = ssterm(2,1) - dt*u12terml*denf            &
     &                *thimp*(phi(k1)+phi0(k01))
                ddterm(2,1) = ddterm(2,1) - dt*u12terml*denf            &
     &              *(thimp*(phi(k1)+phi0(k01))-(0.5*phi(k1)+phi0(k01)))
!
!>>> Changed +dt's to -dt's 10/19/06.  NMF.
                ssterm(2,2) = ssterm(2,2) - dt*u78terml*denf            &
     &                *thimp*(phi(k2)+phi0(k2))
                ddterm(2,2) = ddterm(2,2) - dt*u78terml*denf            &
     &               *(thimp*(phi(k2)+phi0(k2))-(0.5*phi(k2)+phi0(k2)))
!
              if(numvar.ge.3) then
!
              ssterm(3,2)=ssterm(3,2) - dt*denf*o2*u6terml(k,l)         &
     &           *thimp*(phi(k3)+phi0(k3))
!
              ddterm(3,2)=ddterm(3,2) - dt*denf*o2*u6terml(k,l)         &
     &           *(thimp*(phi(k3)+phi0(k3))-(0.5*phi(k3)+phi0(k3)))
!
              ssterm(3,3) = ssterm(3,3)                                 &
     &            -dt*denf*o2*u5terml(k,l)*thimp*(phi(k2)+phi0(k2))         
!
              ddterm(3,3) = ddterm(3,3)                                 &
     &            -dt*denf*o2*u5terml(k,l)*(thimp*(phi(k2)+phi0(k2))    &
     &                                 -(0.5*phi(k2)+phi0(k2)))         &
!
!

              endif ! on numvar.gt.3
              endif ! on numvar.gt.2
              endif ! on idens.gt.0

              ! anisotropic heat flux
 
              if (numvar.ge.3. .and. chipar.ne.0.) then
                 t24terml = t2terml(k,l) + t4terml(k,l)
                 t35terml = t3terml(k,l) + t5terml(k,l)

                 temp = chif*                                           &
     &                (t35terml*tw5(lone)-t24terml*tw4(lone))*          &
     &                (phi(k1)+phi0(k01))
                 ssterm(3,1) = ssterm(3,1) + thimp*dt*temp              &
                 ddterm(3,1) = ddterm(3,1) + thimp*dt*temp -            &
     &                chif*dt*                                          &
     &                (t35terml*tw50(lone)-t24terml*tw40(lone))*        &
     &                (phi(k1)/2.+phi0(k01))

                 temp = chif*                                           &
     &                (-t1terml(k,l)*tw1(lone)                          &
     &                +t2terml(k,l)*tw2(lone)                           &
     &                -t3terml(k,l)*tw3(lone))*                         &
     &                (phi(k1)+phi0(k01))
                 ssterm(3,3) = ssterm(3,3) +  thimp    *dt*temp
                 ddterm(3,3) = ddterm(3,3) + (thimp-1.)*dt*temp         &
              endif

              enddo ! on l
            enddo ! on k

            if(i1.eq.0 .or. j1.eq.0) then
              write(*,2001) i,j,i1,j1,itri,isvaln(itri,i),isvaln(itri,j)
 2001          format(" before call increment,4",7i5 )
               call safestop(409)
            endif
            call sparseR8d_increment_array(s2matrix_lu,i1,j1,           &
     &           ssterm,numvar)
            call sparseR8d_increment_array(d2matrix_lu,i1,j1,           &
     &           ddterm,numvar)
            call sparseR8d_increment_array(r2matrix_lu,i1,j1,           &
     &           rrterm,numvar)
            call sparseR8d_increment_array(q2matrix_lu,i1,j1,           &
     &           qqterm,numvar)
 800        continue
         enddo  
      enddo  
      if(myrank.eq.0 .and. itimer.gt.0) call second(tafterelements)
      telements = telements + tafterelements-tafterread
      enddo                     ! on j
      enddo                     ! on i
      enddo                     ! on iodd
#ifdef mpi
      icount = 6*numvar*m*n
      vtemp = ohmic
      call MPI_ALLREDUCE(vtemp,ohmic,icount,MPI_DOUBLE_PRECISION,       &
     &     MPI_SUM, MPI_COMM_WORLD,ier)
      if(ier /= 0) then
         print *,'Error with MPI_allreduce-ohmic',ier
         call safestop(1)
      endif
      vtemp = viscous
      call MPI_ALLREDUCE(vtemp,viscous,icount,MPI_DOUBLE_PRECISION,     &
     &     MPI_SUM, MPI_COMM_WORLD,ier)
      if(ier /= 0) then
         print *,'Error with MPI_allreduce-viscous',ier
         call safestop(1)
      endif
#endif
!     
      if(myrank.eq.0 .and. itimer.gt.0) call second(tbeforeinit)
      tmpi = tmpi + tbeforeinit-tafterelements
!.....modify the s-matrix, inserting the boundary conditions
!     
!.....calculate boundary conditions

      call boundaryp(iboundp,nbcp)
      if(nbcp .gt. iboundmax) then
         write(*,4888) nbcp, iboundmax
 4888    format(" ERROR: nbcp > iboundmax", 2i5)
         call safestop(9) 
      endif

      do l=1,nbcp
         i = iboundp(l)
         if(sparseR8d_is_local_row(s2matrix_lu, i).eq.1) then
            if(i.le.0 .or. i.ge.m*n*6*numvar+1) then
               write(*,3342) i,l,iboundp(l)
 3342          format(" error in iboundp index",3i5)
               call safestop(3342)
            endif
            call sparseR8d_zero_row(s2matrix_lu,i)
            call sparseR8d_increment(s2matrix_lu,i,i,1.,numvar,1)
         endif
      enddo
      
      if(myrank.eq.0) call second(tafterinit)
      tinit = tinit + tafterinit-tbeforeinit

      if(ifirsts2_lu.eq.0) then
         ifirsts2_lu = 1
         call set_superlu_options(s2matrix_lu%options, Fact=0)
      else
         call set_superlu_options(s2matrix_lu%options, Fact=2)
      endif

      if(myrank.eq.0 .and. itimer.gt.0) then
         if(myrank.eq.0) call second(tendludefphi)
         tsolve = tsolve + tendludefphi-tafterinit
      endif

      return
      end
!======================================================
      subroutine ludefden
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
#ifdef mpi
      use supralu_dist_mod
      use superlu_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      real  ssterm(1,1),ddterm(1,1),rrterm(3,3),qqterm(3,3)
      integer itri, irect, jrect, ll, iodd
      integer i, ione, i1
      integer j, jone, j1
      integer k, kone, k1, k3
      integer l
      integer ier
      real x, z, terma, termb, termd
      real hypp, factor, coef
      real tafterelements, taftermpi
      real tludefstart !, tludefend
!      integer icount
!
c$$$      if(ifirsts8_lu.eq.1) then
c$$$         call f_Destroy_LU(s8matrix_lu%n, SLUprocgrid,                     &
c$$$     &        s8matrix_lu%LUstruct)
c$$$      endif
      call sparseR8d_reset(s8matrix_lu, ier)
      call sparseR8d_reset(d8matrix_lu, ier)
      call sparseR8d_reset(r8matrix_lu, ier)
      call sparseR8d_reset(q8matrix_lu, ier)
      s8matrix_lu%colind = origcolind1
      d8matrix_lu%colind = origcolind1
      r8matrix_lu%colind = origcolind
      q8matrix_lu%colind = origcolind

      itri = 0
      if(myrank.eq.0 .and. itimer.gt.0) call second(tludefstart)

      do iodd=1,ioddmx
!     
      do i=1,18
      do j=1,18
!
      do jrect=1,m-1+jper
      do irect=1,n-1+iper
!
         ll = irect + (jrect-1)*(n-1+iper)
         itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
         i1 = isvaln(itri,i)         
         ione = isval1(itri,i)

         ! check whether this node is assigned to this processor
         if(sparseR8d_is_local_row(r8matrix_lu, i1).eq.0) goto 800

         j1 = isvaln(itri,j)
         jone = isval1(itri,j)

         x = (irect-iodd/3.)*deex
         z = (jrect-1+iodd/3.)*deez
         factor = 1
         if(imask.eq.1) call mask(x,z,factor)
!     
         ssterm = 0.
         ddterm = 0.
         rrterm = 0.
         qqterm = 0.
!     
         terma = aterm(iodd,i,j)
         termb = bterm(iodd,i,j)
         termd = dterm(iodd,i,j)
         hypp = hyperp*deex**2*termb
         
!     
!...NOTE: s8matrix_lu (formed from ss) multiplies the solution at time NTIME+1
!
! LINEAR TERMS
!
!....note "artificial diffusion" proportional to dt and denm in density
            coef =dt*denm
            ssterm(1,1) =  termd - thimp*dt*coef*terma                  &
     &                           + thimp*dt*coef*hypp
            ddterm(1,1) = termd + (1.-thimp)*dt*coef*terma              &
     &                           - (1.-thimp)*dt*coef*hypp
!
         do k=1,18
            k1= isvaln(itri,k)
            kone = isval1(itri,k)
            k3 = k1+12
!
! NONLINEAR TERMS
!
            ssterm(1,1)=ssterm(1,1) - dt*thimp                          &
     &           *(k0term(iodd,i,k,j)*(veln(k1)+vel0(k1)))
!     
            ddterm(1,1)=ddterm(1,1) + dt*k0term(iodd,i,k,j)*            &
     &           (-thimp*(veln(k1)+vel0(k1))+0.5*veln(k1)+vel0(k1))
!     
            rrterm(1,1)=rrterm(1,1) + thimp*dt*k0term(iodd,i,j,k)*      &
     &           (den(kone)+den0(kone))
!     
            qqterm(1,1)=qqterm(1,1) + dt*k0term(iodd,i,j,k)*            &
     &          (-thimp*(den(kone)+den0(kone))+0.5*den(kone)+den0(kone))
!     
!     
            if(numvar.ge.3) then
!     
               ssterm(1,1)=ssterm(1,1) + dt*thimp*                      &
     &              g2term(iodd,i,k,j)*(veln(k3)+vel0(k3))
!     
               ddterm(1,1)=ddterm(1,1)                                  &
     &              +dt*g2term(iodd,i,k,j)*(thimp*(veln(k3)+vel0(k3))    &
     &              -(0.5*veln(k3)+vel0(k3)))
!
!     
               rrterm(1,3)=rrterm(1,3) - thimp*dt*g2term(iodd,i,j,k)*   &
     &              (den(kone)+den0(kone))
!---> NOTE:  Plus sign in preceeding equation changed to minus on 1/09/06   SCJ
!     
               qqterm(1,3)=qqterm(1,3) + dt*g2term(iodd,i,j,k)*         &
     &         (thimp*(den(kone)+den0(kone))-(0.5*den(kone)+den0(kone)))
!
!
            endif
!
         enddo                  ! on k
         if(i1.eq.0 .or. j1.eq.0) then
            write(*,2001) i,j,i1,j1,itri,isvaln(itri,i),isvaln(itri,j)
 2001       format(" before call increment,4",7i5 )
            call safestop(409)
         endif
         call sparseR8d_increment(s8matrix_lu,ione,jone,ssterm(1,1),1,0)
         call sparseR8d_increment(d8matrix_lu,ione,jone,ddterm(1,1),1,0)
         call sparseR8d_increment_array(r8matrix_lu,i1,j1,rrterm,numvar)
         call sparseR8d_increment_array(q8matrix_lu,i1,j1,qqterm,numvar)
 800     continue
      enddo                     ! on irect
      enddo                     ! on jrect
      enddo                     ! on j
      enddo                     ! on i
      enddo                     ! on iodd

      if(myrank.eq.0 .and. itimer.gt.0) call second(taftermpi)
      tmpi = tmpi + taftermpi - tafterelements
!     
!.....modify the s-matrix, inserting the boundary conditions
!     
!.....calculate boundary conditions
      call boundaryds(iboundn,nbcn,1)
      if(nbcn .gt. iboundmax) then
         write(*,4888) nbcn, iboundmax
 4888    format(" ERROR: nbcn > iboundmax", 2i5)
         call safestop(9) 
      endif
      do l=1,nbcn
         i = iboundn(l)
         if(sparseR8d_is_local_row(s8matrix_lu, i).eq.1) then
            if(i.le.0 .or. i.ge.m*n*6+1) then
               write(*,3342) i,l,iboundp(l)
 3342          format(" error in iboundp index",3i5)
               call safestop(3342)
            endif
            call sparseR8d_zero_row(s8matrix_lu, i)
            call sparseR8d_increment(s8matrix_lu,i,i,1.,1,1)
         endif
      enddo

      if(ifirsts8_lu.eq.0) then
         ifirsts8_lu = 1
         call set_superlu_options(s8matrix_lu%options, Fact=0)
      else
         call set_superlu_options(s8matrix_lu%options, Fact=2)
      endif

      return
      end
!======================================================
      subroutine ludefpres
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
#ifdef mpi
      use supralu_dist_mod
      use superlu_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer itri, iodd, irect, jrect, ll, ier
      integer i, ione, i1
      integer j, jone, j01, j1
      integer k, kone, k01, k1, k3
      integer l
      real ssterm(1,1),ddterm(1,1),rrterm(3,3),qqterm(3,3)
      real terma, termb, termd, factor, x, z, hypp
      integer icount, ir1, ir2, ir3

c$$$      if(ifirsts9_lu.eq.1) then
c$$$         call f_Destroy_LU(s9matrix_lu%n, SLUprocgrid,                     &
c$$$     &        s9matrix_lu%LUstruct)
c$$$      endif
      call sparseR8d_reset(s9matrix_lu, ier)
      call sparseR8d_reset(d9matrix_lu, ier)
      call sparseR8d_reset(r9matrix_lu, ier)
      call sparseR8d_reset(q9matrix_lu, ier)
      s9matrix_lu%colind = origcolind1
      d9matrix_lu%colind = origcolind1
      r9matrix_lu%colind = origcolind
      q9matrix_lu%colind = origcolind

      q4 = 0.
      itri = 0

      ! prepare files for anisotropic heat flux
      if (chipar.ne.0 .and. numvar.ge.3) then
         rewind(173)
         rewind(174)
         rewind(175)
         rewind(176)
         rewind(177)
      endif

      do iodd=1,ioddmx

      do i=1,18
      do j=1,18

      if(iread(73).eq.1)                                                &
     &     read(173) ir1,ir2,ir3,((t1terml(k,l),l=1,18),k=1,18)
      if(iread(74).eq.1)                                                &
     &     read(174) ir1,ir2,ir3,((t2terml(k,l),l=1,18),k=1,18)
      if(iread(75).eq.1)                                                &
     &     read(175) ir1,ir2,ir3,((t3terml(k,l),l=1,18),k=1,18)
      if(iread(76).eq.1)                                                &
     &     read(176) ir1,ir2,ir3,((t4terml(k,l),l=1,18),k=1,18)
      if(iread(77).eq.1)                                                &
     &     read(177) ir1,ir2,ir3,((t5terml(k,l),l=1,18),k=1,18)

      do jrect=1,m-1+jper
      do irect=1,n-1+iper

         ll = irect + (jrect-1)*(n-1+iper)
         itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
         i1 = isvaln(itri,i)
         ione = isval1(itri,i)

         ! check whether this node is assigned to this processor
         if(sparseR8d_is_local_row(r9matrix_lu, i1).eq.0) goto 800
         if(sparseR8d_is_local_row(s9matrix_lu, ione).eq.0) then
            write(*,*) "weirdness, i1, ione", i1, ione
            goto 800
         endif

         x = (irect-iodd/3.)*deex
         z = (jrect-1+iodd/3.)*deez
         factor = 1
         if(imask.eq.1) call mask(x,z,factor)

         ssterm = 0.
         ddterm = 0.
         rrterm = 0.
         qqterm = 0.
     
         terma = aterm(iodd,i,j)
         termb = bterm(iodd,i,j)
         termd = dterm(iodd,i,j)
         hypp = hyperp*deex**2*termb
         j1 = isvaln(itri,j)
         j01 = isval0(itri, j)
         jone = isval1(itri,j)
!     
!...NOTE: s9matrix_lu (formed from ss) multiplies the solution at time NTIME+1
!     
! LINEAR TERMS

         ssterm(1,1) =  termd -     thimp *dt*kappa*(terma-hypp*termb)
         ddterm(1,1) =  termd + (1.-thimp)*dt*kappa*(terma-hypp*termb)
!     
         q4(i1) = q4(i1) + dt*termd*sphip(jone)
!     
         do k=1,18
            k01 = isval0(itri,k)
            k1= isvaln(itri,k)
            kone = isval1(itri,k)
            k3 = k1+12
!     
! NONLINEAR TERMS
!     
            ssterm(1,1)=ssterm(1,1) - dt*thimp                          &
     &           *(k0term(iodd,i,k,j)*(veln(k1)+vel0(k1))               &
     &           +(gam-1.)*kappat*(g2term(iodd,i,k,j)+g2term(iodd,i,j,k)&
     &                           - hypp*g13erm(iodd,i,k,j))*deni(kone))
!     
            ddterm(1,1)=ddterm(1,1) + dt*(k0term(iodd,i,k,j)*           &
     &           (-thimp*(veln(k1)+vel0(k1))+0.5*veln(k1)+vel0(k1))     &
     &          +(gam-1.)*kappat*(g2term(iodd,i,k,j)+g2term(iodd,i,j,k) &
     &                - hypp*g13erm(iodd,i,k,j))*(1.-thimp)*deni(kone))
!     
            rrterm(1,1)=rrterm(1,1) + thimp*dt*k0term(iodd,i,j,k)*      &
     &           (pres(kone)+pres0(kone))
!     
            qqterm(1,1)=qqterm(1,1) + dt*k0term(iodd,i,j,k)*            &
     &           (-thimp*(pres(kone)+pres0(kone))                       &
     &           +0.5*pres(kone)+pres0(kone))
!     
            if(numvar.ge.3) then
!     
               ssterm(1,1)=ssterm(1,1) + dt*thimp*                      &
     &              (g2term(iodd,i,k,j)-(gam-1.)*k2term(iodd,i,k,j))    &
     &              *(veln(k3)+vel0(k3))
!     
               ddterm(1,1)=ddterm(1,1)                                  &
     &              +dt*(g2term(iodd,i,k,j)-(gam-1.)*k2term(iodd,i,k,j))&
     &              *(thimp*(veln(k3)+vel0(k3))-(0.5*veln(k3)+vel0(k3)))
!     
               rrterm(1,3)=rrterm(1,3) - thimp*dt
     &              *(g2term(iodd,i,j,k)-(gam-1.)*k2term(iodd,i,j,k))*  &
     &              (pres(kone)+pres0(kone))
!     
               qqterm(1,3)=qqterm(1,3)                                  &
     &              +dt*(g2term(iodd,i,j,k)-(gam-1.)*k2term(iodd,i,j,k))&
     &              *(thimp*(pres(kone)+pres0(kone))                    &
     &              -(0.5*pres(kone)+pres0(kone)))
            endif

         enddo                  ! on k
         if(i1.eq.0 .or. j1.eq.0) then
            write(*,2001) i,j,i1,j1,itri,isvaln(itri,i),isvaln(itri,j)
 2001       format(" before call increment,4",7i5 )
            call safestop(409)
         endif
         call sparseR8d_increment(s9matrix_lu,ione,jone,ssterm(1,1),1,0)
         call sparseR8d_increment(d9matrix_lu,ione,jone,ddterm(1,1),1,0)
         call sparseR8d_increment_array(r9matrix_lu,i1,j1,rrterm,numvar)
         call sparseR8d_increment_array(q9matrix_lu,i1,j1,qqterm,numvar)
 800     continue
      enddo                     ! on irect
      enddo                     ! on jrect
      enddo                     ! on j
      enddo                     ! on i
      enddo                     ! on iodd

#ifdef mpi
      icount = m*n*numvar*6
      vtemp = q4
      call MPI_ALLREDUCE(vtemp,q4,icount,MPI_DOUBLE_PRECISION,          &
     &     MPI_SUM, MPI_COMM_WORLD,ier)
      if(ier /= 0) then
         print *,'Error with MPI_allreduce-q4',ier
         call safestop(1)
      endif
#endif
      if(myrank.eq.0) write(*,*) "inserting bcs"
!     
!.....modify the s-matrix, inserting the boundary conditions
!     
!.....calculate boundary conditions
      call boundaryds(iboundn,nbcn,1)
      if(nbcn .gt. iboundmax) then
         write(*,4888) nbcn, iboundmax
 4888    format(" ERROR: nbcn > iboundmax", 2i5)
         call safestop(9) 
      endif
      do l=1,nbcn
         i = iboundn(l)
         if(sparseR8d_is_local_row(s9matrix_lu, i).eq.1) then
            if(i.le.0 .or. i.ge.m*n*6+1) then
               write(*,3342) i,l,iboundp(l)
 3342          format(" error in iboundp index",3i5)
               call safestop(3342)
            endif
            call sparseR8d_zero_row(s9matrix_lu, i)
            call sparseR8d_increment(s9matrix_lu,i,i,1.,1,1)
         endif
      enddo

      if(ifirsts9_lu.eq.0) then
         ifirsts9_lu = 1
         call set_superlu_options(s9matrix_lu%options, Fact=0)
      else
         call set_superlu_options(s9matrix_lu%options, Fact=2)
      endif

      return
      end
!==============================================================
      subroutine boundaryv(ibound,nbc)
      use basic
      use arrays
      implicit none
      real x, ak, epsdot
      integer index, indext, lx, lz, nbc
!
!.....defines the velocity boundary conditions for a rectangular region
      integer ibound(*)
!
      velbounds = 0
!
!.....FIRST EQUATION:  Velocity stream function PHI
      index = 0
      indext=0
      if(jper.eq.1) go to 506
!.....start loop on top and bottom boundaries (y-direction)
      do lx=1,n
!
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 1
!.....tangential derivative (x)
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 2
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 2
!.....tangential second derivatives (xx)
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 4
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 4
!.....normal second derivatives (yy) (on bottom, for non-symmetric problem only)
!         if(jsym.eq.0) then
         index = index + 1
         ibound(index  ) = 6*numvar*(          lx-1) + 6
!         endif
         index = index + 1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 6
!
      enddo                     ! loop along top and bottom boundaries
!
!.....special for setting up boundary conditions for Taylor Problem
      if(itaylor.eq.1) then
         do lx=1,n
            x = (lx-1)*deex - alx/2.*(1-isym)
            if(tau.eq.0) tau = 1.
            epsdot = .01*(time/tau**2)*exp(-time/tau)
            ak = 2.*pi/alx
            indext = indext+2
            velbounds(indext-1) = (1./ak)*epsdot*sin(ak*x)
            velbounds(indext  ) =-(1./ak)*epsdot*sin(ak*x)
!.....tangential derivative (x)
            indext = indext+2
            velbounds(indext-1) =-epsdot*cos(ak*x)
            velbounds(indext  ) = epsdot*cos(ak*x)
!.....tangential second derivatives (xx)
            indext = indext+2
            velbounds(indext-1) = (ak)*epsdot*sin(ak*x)
            velbounds(indext  ) =-(ak)*epsdot*sin(ak*x)
         enddo ! loop over top and bottom boundaries for Taylor problem
      endif ! for Taylor Problem boundary conditions
 506  continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 511
      if(m.le.2) go to 501
      do lz=1,m-2
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 1
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 1
!.....tangential derivative (y)
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 3
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 3
!.....tangential second derivatives (yy)
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 6
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 6
!.....normal second derivative (xx)
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 4
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 4
      enddo
 501  continue
!.....corner points
      index = index+4
      ibound(index-3) =                 3
      ibound(index-2) = 6*numvar*(n-1     ) + 3
      ibound(index-1) = 6*numvar*((m-1)*n ) + 3
      ibound(index  ) = 6*numvar*( m*n - 1) + 3
 511  continue ! end of branch on iper
      if(numvar.le.1) then
        nbc = index
        return
      endif
!.....SECOND EQUATION:  Toroidal Velocity Vz
!.....start loop on top and bottom boundaries (x-direction)
      if(jper.eq.1) go to 516
      do lx=1,n
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 7
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 7
!.....tangential derivative (x)
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 8
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 8
!.....tangential second derivative (xx)
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 10
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 10
        if(hyperv.ne.0 .and. imask.ne.1) then
!......normal derivative (y) and (xy)
           if(jsym.eq.0) then
              index = index+1
              ibound(index  ) = 6*numvar*(       lx-1) + 9
           endif
           index = index+1
           ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 9
           if(jsym.eq.0) then
              index = index+1
              ibound(index  ) = 6*numvar*(       lx-1) + 11
           endif
           index = index+1
           ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 11
        endif
      enddo
  516 continue

!.....start loop on left and right boundaries (z-direction)
      if(iper.eq.1) go to 512
      if(m.le.2) go to 502
      do lz=1,m-2
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 7
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 7
!.....tangential derivative (y)
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 9
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 9
!.....tangential second derivative (yy)
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 12
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 12
         if(hyperv.ne.0 .and. imask.ne.1) then
!......normal second derivative (xx)
            index = index+2
            ibound(index-1) = 6*numvar*( lz*n      ) + 8
            ibound(index  ) = 6*numvar*( lz*n + n-1) + 8
            index = index+2
            ibound(index-1) = 6*numvar*( lz*n      ) + 11
            ibound(index  ) = 6*numvar*( lz*n + n-1) + 11
         endif
      enddo ! end of loop and left and right boundaries
 502  continue
!.....corner points
      index = index+4
      ibound(index-3) =                       12
      ibound(index-2) = 6*numvar*(n-1     ) + 12
      ibound(index-1) = 6*numvar*((m-1)*n ) + 12
      ibound(index  ) = 6*numvar*( m*n - 1) + 12
      index = index+4
      ibound(index-3) =                       9
      ibound(index-2) = 6*numvar*(n-1     ) + 9
      ibound(index-1) = 6*numvar*((m-1)*n ) + 9
      ibound(index  ) = 6*numvar*( m*n - 1) + 9
 512  continue
      if(numvar.le.2) then
        nbc = index
        return
      endif
!.....THIRD EQUATION:  Velocity potential Chi
!.....start loop on top and bottom boundaries (z-direction)
      if(jper.eq.1) go to 526
      do lx=1,n
!.....normal derivative
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 15
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 15
!.....cross derivative
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 17
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 17
      enddo
  526 continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 513
      if(m.le.2) go to 503
      do lz=1,m-2
!.....normal derivative
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 14
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 14
!.....cross derivative
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 17
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 17
      enddo
 503  continue
!.....X derivatives
      index = index+4
      ibound(index-3) =                 14
      ibound(index-2) = 6*numvar*(n-1     ) + 14
      ibound(index-1) = 6*numvar*((m-1)*n ) + 14
      ibound(index  ) = 6*numvar*( m*n - 1) + 14
 513  continue
!
!     Fix one point to zero:  Note, this will lead to different offset
!     for symmetric and non-symmetric problem, but other answers should
!     be ok....
!      index = index + 1
!      ibound(index) = 6*numvar*(m*n-1) + 13
      nbc = index
      return
      end
! BOUNDARYP ==================================================================
      subroutine boundaryp(ibound,nbc)
      use basic
      use arrays
      implicit none
      real x, z, ak
      integer index, lx, lz, nbc
!
!.....defines the boundary conditions for a rectangular region
      integer ibound(*)
      psibounds = 0
      if(tau.eq.0) tau = 1.
      eps = .01*(1.-(1.+time/tau)*exp(-time/tau))
      ak = 2.*pi/alx
!
!.....FIRST EQUATION:  poloidal flux psi
      index = 0
      if(jper.eq.1) go to 506
!.....start loop on top and bottom boundaries (y-direction)
      do lx=1,n
         x = (lx-1)*deex - alx/2.*(1-isym)
!.....TOP BOUNDARY
!
         index = index+1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 1
         if(itaylor.eq.1) then
            z = alz/2
            psibounds(index  )  =-eps*z*cos(ak*x)
         endif
!     
!.....tangential derivative (x)
         index = index+1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 2
         if(itaylor.eq.1) then
            z = alz/2
            psibounds(index  )  =  ak*eps*z*sin(ak*x)
         endif
!.....tangential second derivatives (xx)
         index = index+1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 4
         if(itaylor.eq.1) then
            z = alz/2
            psibounds(index  )  = +ak**2*eps*z*cos(ak*x)
         endif
         
         if(hyper.ne.0 .and. imask.ne.1) then
!.......normal derivative (y)
            index = index+1
            ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 3
            if(itaylor.eq.1) then
               psibounds(index  )  =-eps*cos(ak*x)
            endif
!.......cross derivative xy
            index = index+1
            ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 5
            if(itaylor.eq.1) then
               psibounds(index  )  =  ak*eps*sin(ak*x)
            endif
         endif
!     
         if(jsym.eq.0) then
!.......BOTTOM BOUNDARY
            index = index+1
            ibound(index  ) = 6*numvar*(          lx-1) + 1
            if(itaylor.eq.1) then
               z = -alz/2.*(1-jsym)
               psibounds(index  )  = eps*z*cos(ak*x)
            endif
!.......tangential derivative (x)
            index = index+1
            ibound(index  ) = 6*numvar*(          lx-1) + 2
            if(itaylor.eq.1) then
               z = -alz/2.*(1-jsym)
               psibounds(index  )  = -ak*eps*z*sin(ak*x)
            endif
!.......tangential second derivatives (xx)
            index = index+1
            ibound(index  ) = 6*numvar*(          lx-1) + 4
            if(itaylor.eq.1) then
               z = -alz/2.*(1-jsym)
               psibounds(index  )  = -ak**2*eps*z*cos(ak*x)
            endif
            
            if(hyper.ne.0 .and. imask.ne.1) then
!.......normal derivative (y)
               index = index+1
               ibound(index  ) = 6*numvar*(          lx-1) + 3
               if(itaylor.eq.1) then
                  psibounds(index  )  = eps*cos(ak*x)
               endif
!.......cross derivative xy
               
               index = index+1
               ibound(index  ) = 6*numvar*(          lx-1) + 5
               if(itaylor.eq.1) then
                  psibounds(index  )  = -ak*eps*sin(ak*x)
               endif
            endif               ! on hyper .ne.0 .and. imask.ne.1
!     
!.....FOR SYMMETRY OPTION
         else                   ! do the following for jsym.eq.1
            index = index+1
            ibound(index) = 6*numvar*(          lx-1) + 3
            index = index+1
            ibound(index) = 6*numvar*(          lx-1) + 5
         endif                  ! branch on symmetry options
!     
!     
!     
      enddo                     ! on bottom boundary
 506  continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 511
      if(m.le.2) go to 501
      do lz=1,m-2
!.......RIGHT BOUNDARY
         index = index+1
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 1
!.......tangential derivative (y)
         index = index+1
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 3
!.......tangential second derivatives (yy)
         index = index+1
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 6
         if(hyper.ne.0 .and. imask.ne.1) then
!.......normal derivative (x)
            index = index+1
            ibound(index  ) = 6*numvar*( lz*n + n-1) + 2
!.......cross derivative (xy)
            index = index+1
            ibound(index  ) = 6*numvar*( lz*n + n-1) + 5
         endif
         if(isym.eq.0) then
!.........LEFT BOUNDARY
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 1
!.........tangential derivative (y)
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 3
!.........tangential second derivatives (yy)
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 6
            if(hyper.ne.0 .and. imask.ne.1) then
!.........normal derivative (x)
               index = index+1
               ibound(index) = 6*numvar*( lz*n      ) + 2
!...........cross derivative (xy)
               index = index+1
               ibound(index) = 6*numvar*( lz*n      ) + 5
            endif
         else
            index = index+2
            ibound(index-1) = 6*numvar*( lz*n      ) + 2
            ibound(index  ) = 6*numvar*( lz*n      ) + 5
         endif
      enddo
 501  continue
      
!.....corner points
!     UPPER RIGHT
      index = index+1
      ibound(index) = 6*numvar*( m*n - 1) + 6
!     UPPER LEFT
      if(isym.eq.0) then
         index = index+1
         ibound(index) = 6*numvar*((m-1)*n ) + 6
      endif
!     BOTTOM RIGHT
      if(jsym.eq.0) then
         index = index+1
         ibound(index) = 6*numvar*(n-1     ) + 6
      endif
!     BOTTOM LEFT
      if(isym.eq.0 .and. jsym.eq.0) then
         index = index+1
         ibound(index) =                 6
      endif
!     
 511  continue
      if(numvar.le.1) then
         nbc = index
         return
      endif
!.....SECOND EQUATION: toroidal field function I
!.....start loop on top and bottom boundaries (z-direction)
      if(jper.eq.1) go to 516
      do lx=1,n
!.....function and tangential derivatives
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 7
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 7
         psibounds(index-1) = gbound
         psibounds(index  ) = gbound
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 8
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 8
         index = index+2
         ibound(index-1) = 6*numvar*(          lx-1) + 10
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 10
      enddo
 516  continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 512
      if(m.le.2) go to 502
      do lz=1,m-2
!.....function
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 7
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 7
         psibounds(index-1) = gbound
         psibounds(index  ) = gbound
!.....tangential derivative (y)
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 9
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 9
!.....tangential second derivative (yy)
         index = index+2
         ibound(index-1) = 6*numvar*( lz*n      ) + 12
         ibound(index  ) = 6*numvar*( lz*n + n-1) + 12
      enddo
 502  continue
!.....corner points
      index = index+4
      ibound(index-3) =                       12
      ibound(index-2) = 6*numvar*(n-1     ) + 12
      ibound(index-1) = 6*numvar*((m-1)*n ) + 12
      ibound(index  ) = 6*numvar*( m*n - 1) + 12
      index = index+4
      ibound(index-3) =                       9
      ibound(index-2) = 6*numvar*(n-1     ) + 9
      ibound(index-1) = 6*numvar*((m-1)*n ) + 9
      ibound(index  ) = 6*numvar*( m*n - 1) + 9
 512  continue
      if(numvar.le.2) then
         nbc = index
         return
      endif
!.....THIRD EQUATION:  Electron Pressure Pe
!.....start loop on top and bottom boundaries (z-direction)
      if(jper.eq.1) go to 526
      do lx=1,n
!.......TOP BOUNDARY
         index = index+1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 13
!.......tangential derivative
         index = index+1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 14
!.......tangential second derivative
         index = index+1
         ibound(index  ) = 6*numvar*((m-1)*n + lx-1) + 16
!.......BOTTOM BOUNDARY
         if(jsym.eq.0)  then
            index = index+1
            ibound(index) = 6*numvar*(          lx-1) + 13
!.........tangential derivative
            index = index+1
            ibound(index) = 6*numvar*(          lx-1) + 14
!.........tangential second derivative
            index = index+1
            ibound(index) = 6*numvar*(          lx-1) + 16
         else
!.........normal derivative
            index = index+1
            ibound(index) = 6*numvar*(          lx-1) + 15
!.........cross derivative
            index = index+1
            ibound(index) = 6*numvar*(          lx-1) + 17
         endif
         
      enddo
 526  continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 513
      if(m.le.2) go to 503
      do lz=1,m-2
!.....RIGHT BOUNDARY
         index = index+1
         ibound(index) = 6*numvar*( lz*n + n-1) + 13
!.......tangential derivative
         index = index+1
         ibound(index) = 6*numvar*( lz*n + n-1) + 15
!.......tangential second derivative
         index = index+1
         ibound(index) = 6*numvar*( lz*n + n-1) + 18
!.....LEFT BOUNDARY
         if(isym.eq.0) then
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 13
!.......tangential derivative
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 15
!.......tangential second derivative
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 18
         else
!     
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 14
            index = index+1
            ibound(index) = 6*numvar*( lz*n      ) + 17
         endif
      enddo
 503  continue
!.....corner points
!     
!.....UPPER RIGHT
      index = index+2
      ibound(index  ) = 6*numvar*( m*n - 1) + 15
      ibound(index-1) = 6*numvar*( m*n - 1) + 18
!.....UPPER LEFT
      if(isym.eq.0) then
         index = index+2
         ibound(index  ) = 6*numvar*((m-1)*n ) + 15
         ibound(index-1) = 6*numvar*((m-1)*n ) + 18
      endif
!.....LOWER RIGHT
      if(jsym.eq.0) then
         index = index+2
         ibound(index  ) = 6*numvar*(n-1     ) + 15
         ibound(index-1) = 6*numvar*(n-1     ) + 18
      endif
!.....LOWER LEFT
      if(isym.eq.0 .and. jsym.eq.0) then
         index = index+2
         ibound(index  ) =                 15
         ibound(index-1) =                 18
      endif
 513  continue
      nbc = index
      return
      end
!============================================================
      subroutine init
      use p_data
      use t_data
      use basic
      use arrays
      use inout_mod
      use mesh_mod
      implicit none
      integer nhalf, mhalf

      if(myrank.eq.0 .and. iprint.gt.0) then
         print *, "Entering init."
      endif

!
      if(myrank.eq.0) call second(tfirst) !Get current system time in seconds
      tread = 0.
      telements = 0.
      tsolve = 0.
      tonestep = 0
      tinit = 0
      tmpi = 0.
      tnewvarb = 0.
!
      ifirstd1_lu = 0
      ifirsts1_lu = 0
      ifirstd2_lu = 0
      ifirsts2_lu = 0
      ifirstr1_lu = 0
      ifirstr2_lu = 0
      ifirstq2_lu = 0
      ifirsts3_lu = 0
      ifirsts4_lu = 0
      ifirsts5_lu = 0
      ifirsts6_lu = 0
      ifirstsa_lu = 0
      ifirsts7_lu = 0
      ifirstd8_lu = 0
      ifirsts8_lu = 0
      ifirstr8_lu = 0
      ifirstq8_lu = 0
      ifirstd9_lu = 0
      ifirsts9_lu = 0
      ifirstr9_lu = 0
      ifirstq9_lu = 0
!
      gbound = 0.
!
!.....define all input parameters (choose by test problem)
      call input

      ntimemin = min(2,ntimemax-1)
      iboundmax = 12*numvar*(n+m-2)
!
      tt = 2.*thimp
!
!.....parameter that specifies grid type
      ioddmx = 2
!
!.....redefine n and m to be odd for ioddmx=8
      if(ioddmx.eq.8) then
        nhalf = (n-1)/2
        n = nhalf*2 + 1
        mhalf = (m-1)/2
        m = mhalf*2 + 1
      endif
!
!      allocate big arrays
      call space
!
!
!.....define the arrays needed to define the sparse matrices:
!     jbig = ibig + incb(ibig,l),l=1,itypemax(ibig)
      call incbdef

!.................................................................
!     PART 1:
!.....initialize needed variables and define geometry and triangles
      alxp = alx*(1./(1.+isym))
      alzp = alz*(1./(1.+jsym))
      deex = alx/(n+iper-1.)*(1./(1.+isym))
      deez = alz/(m+jper-1.)*(1./(1.+jsym))
      if(ntimemax .gt. ntimep) then
            write(*,6665) ntimemax,ntimep
 6665 format(" ntimemax .gt. ntimep", 2i5)
            call safestop(18) 
      endif
!
!
!.....number of rectangular regions (ioddmx triangles per region)
      nreg = (n+iper-1)*(m+jper-1)*2/ioddmx
!
!.....zero out big arrays
!     b1vecini = 0.
!     b2vecini = 0.
      d2term = 0
      aterm = 0
      bterm = 0
      dterm = 0
!      eterm = 0          ! eterm does not appear to be used
      g0term = 0
      g4term = 0
      g5term = 0
      g6term = 0
      g7term = 0.
      g8term = 0.
      g9term = 0.
      g10erm = 0.
      g11erm = 0.
      g12erm = 0.
      g13erm = 0.
      g14erm = 0.
      g15erm = 0.
      g16erm = 0.
      h3term = 0
      h5term = 0
      k0term = 0.
      hterm = 0.
      sterm = 0.
      rinv = 0.
      xterm = 0.
      yterm = 0.
      x0term = 0.
      y0term = 0.
      x1term = 0.
      x2term = 0.
      y1term = 0.
      y2term = 0.
!
!.....define properties of triangles
      call tridef
!
!.....define new terms needed for boundary integrals
!.....define 1/R for toroidal problem
      if(itor.eq.1) call rinvdef(rinv)
!
!     open all needed output files
      call openf
!
!.....end of section defining geometrical parameters for triangles

      if(myrank.eq.0 .and. iprint.gt.0) then
         print *, "Exiting init."
      endif

      return
      end
!============================================================
      subroutine rhsdef
!
      use p_data
      use t_data
      use basic
      use arrays
      use mesh_mod
!
      implicit none
      integer :: iodd, itri, i, ii, iii, j,  k, l, ll, i1, i2, i3
      real :: fac1, fac2, fac3, sum, term
      real :: f, fact
!
      itri = 0
      do 505 iodd=1,ioddmx
!     calculate matrix elements of the integration matrix fint
         do i=-4,maxi
            do j=-4,maxi
               fint(i,j) = 0
               if(i.ge.0 .and.j.ge.0)                                   &
     &              fint(i,j) = f(i,j,atri(iodd),btri(iodd),ctri(iodd))
            enddo
         enddo
      do 499 ll=1,nreg
      itri = itri + 1

!.....x,z dependent part follows
!
!     expand x and z in terms of si and eta for this element
      call local(xi,zi,ll,btri(iodd),ttri(iodd),iodd)
!
!
!     calculate the RHS function and it's first 2 derivatives
      call fun(xi(1),zi(1),df,alx,alz)  
!
!.....start loops to compute the RHS vector for this triangular element,
!
       do iii=1,3
         do ii=1,6
           i = (iii-1)*6 + ii
           isval1(itri,i) =         6*ist(itri,iii)+ii
           isval2(itri,i) =        12*ist(itri,iii)+ii
           isvaln(itri,i) =  numvar*6*ist(itri,iii)+ii
           isval0(itri,i) =  numvar*6*ist0(itri,iii)+ii
           i1 = isvaln(itri,i)
           i2 = i1 + 6
           i3 = i2 + 6
!
            term = 0
            do k=1,20
!
             sum = 0
             do j=0,4
             do l=0,j
               fac1 = df(j,l)
               do p=0,l
               fac2 = fac1*xi(2)**(l-p)  *xi(3)**p/(fact(p)*fact(l-p))
               do q=0,j-l
               fac3 = fac2*zi(2)**(j-l-q)*zi(3)**q/(fact(q)*fact(j-l-q))
!
               sum = sum + fac3*fint(mi(k)+j-p-q,ni(k)+p+q)
!
               enddo
               enddo
             enddo
             enddo
!
             term =  term + gtri(k,i,iodd)*sum
            enddo
!
!............velocity perturbations
!            b1vecini(i1) = b1vecini(i1) - term                         &
!    &                    *((pi/alx)**2+(pi/alz)**2)
!            if(numvar.ge.2) b1vecini(i2) = b1vecini(i2)
!    &              - term*((pi/alx)**2+(pi/alz)**2)
!            if(numvar.ge.3) b1vecini(i3) = 0.
!
!............fluxes
!            b2vecini(i1) = 0.
!            if(numvar.ge.2) b2vecini(i2) = 0.
!            if(numvar.ge.3) b2vecini(i3) = 0.
!
           enddo
           enddo
!.....end of do loops on i
!
 499  continue
 505  continue
!
      return
      end
!============================================================
      subroutine energy 
      use p_data
      use t_data
      use basic
      use arrays
!
      implicit none
      integer :: iodd, itri, irect, jrect, ll
      integer :: i, ione, i01, i1, i2, i3 
      integer :: j, jone, j01, j1, j2, j3
      integer :: k, kone
      real :: x, z, terma, termb, termd, hypf, hypi, hypv, hypc
      real :: gamfac, factor

      ekino = ekin
      emago = emag
      ekindo = ekind
      emagdo = emagd
!
      ekinpo = ekinp
      emagpo = emagp
      ekinpdo = ekinpd
      emagpdo = emagpd
!
      ekinto = ekint
      emagto = emagt
      ekintdo = ekintd
      emagtdo = emagtd
!
      ekinpho = ekinph
      ekintho = ekinth
      emagpho = emagph
      emagtho = emagth
!
      ekin3o = ekin3 
      ekin3do = ekin3d
      ekin3ho = ekin3h 
      emag3o = emag3
      emag3do = emag3d
      emag3ho = emag3h
!
      ekin = 0.
      emag = 0.
      ekind = 0.
      emagd = 0.
      ekinp = 0.
      emagp = 0.
      ekinpd = 0.
      emagpd = 0.      
      ekint = 0.
      emagt = 0.
      ekintd = 0.
      emagtd = 0.
      ekinph = 0
      ekinth = 0.
      emagph = 0.
      emagth = 0.
      ekin3 = 0.
      ekin3d = 0.
      ekin3h = 0.
      emag3 = 0.
      emag3d = 0.
      emag3h = 0.
!
!     volume terms
      itri = 0
      do iodd=1,ioddmx
!      do ll=1,nreg
       do jrect=1,m-1+jper
       do irect=1,n-1+iper
       ll = irect + (jrect-1)*(n-1+iper)
       x = (irect-iodd/3.)*deex
       z = (jrect-1+iodd/3.)*deez
        itri = itri + 1
        factor = 1
        if(imask.eq.1) call mask(x,z,factor)
        do j=1,18
        do i=1,18
          terma = aterm(iodd,i,j)
          termb = bterm(iodd,i,j)
          termd = dterm(iodd,i,j)
          hypf = hyper*deex**2
          hypi = hyperi*deex**2
          hypv = hyperv*deex**2
          hypc = hyperc*deex**2
!
          j1 = isvaln(itri,j)
          i1 = isvaln(itri,i)
          ione = isval1(itri,i)
          jone = isval1(itri,j)
          j01 = isval0(itri,j)
          i01 = isval0(itri,i)
!
          if(idens.eq.0) then
            ekinp = ekinp - .5*terma*(1+jsym)*                          &
     &      ((vel(i1)+vel0(i1))*(vel(j1)+vel0(j1))-vel0(i1)*vel0(j1))
          else
            do k=1,18
            kone = isval1(itri,k)
            ekinp = ekinp - .5*g2term(iodd,i,j,k)*dent(kone)*(1+jsym)   &
     &      *((vel(i1)+vel0(i1))*(vel(j1)+vel0(j1))-vel0(i1)*vel0(j1))
            enddo
          endif
!
          emagp = emagp - .5*terma*(1+jsym)*                            &
     &    ((phi(i1)+phi0(i01))*(phi(j1)+phi0(j01))-phi0(i01)*phi0(j01))
!
          ekinpd= ekinpd - termb*amu*(1+jsym)*                          &
     &    ((vel(i1)+vel0(i1))*(vel(j1)+vel0(j1))-vel0(i1)*vel0(j1))
!
          if(linear.eq.1) then
          emagpd= emagpd - termb*etar*(1+jsym)*                         &
     &    ((phi(i1)+phi0(i01))*(phi(j1)+phi0(j01))-phi0(i01)*phi0(j01))
          else
          emagpd= emagpd - termb*etar*(1+jsym)*                         &
     &    (phi(i1)+phi0(i01))*(phi(j1)+phi0(j01))
          endif
!
          ekinph= ekinph + terma*hypc*amu*(1+jsym)*                     &
     &    ((vor(ione)+vor0(ione))*(vor(jone)+vor0(jone))-               &
     &                vor0(ione)*vor0(jone))
!
          if(linear.eq.1) then
          emagph= emagph + terma*hypf*etar*(1+jsym)*                    &
     &    ((jphi(ione)+jphi0(ione))*(jphi(jone)+jphi0(jone))-           &
     &                jphi0(ione)*jphi0(jone))
          else
          emagph= emagph + terma*hypf*etar*(1+jsym)*                    &
     &    (jphi(ione)+jphi0(ione))*(jphi(jone)+jphi0(jone))       
          endif
!
!
          if(numvar.ge.2) then
          j2 = j1 + 6
          i2 = i1 + 6
          if(idens.eq.0) then
            ekint = ekint + .5*termd*(1+jsym)*                          &
     &      ((vel(i2)+vel0(i2))*(vel(j2)+vel0(j2))-vel0(i2)*vel0(j2))
          else
            do k=1,18
            kone = isval1(itri,k)
            ekint = ekint + 0.5*k1term(iodd,i,j,k)*dent(kone)*(1+jsym)  &
     &      *((vel(i2)+vel0(i2))*(vel(j2)+vel0(j2))-vel0(i2)*vel0(j2))
            enddo
          endif
!
          emagt = emagt + .5*termd*(1+jsym)*                            &
     &    ((phi(i2)+phi0(i2))*(phi(j2)+phi0(j2))-phi0(i2)*phi0(j2))
!
          ekintd= ekintd + amu*terma*(1+jsym)*                          &
     &    ((vel(i2)+vel0(i2))*(vel(j2)+vel0(j2))-vel0(i2)*vel0(j2))
          ekinth= ekinth - amu*termb*hypv*(1+jsym)*                     &
     &    ((vel(i2)+vel0(i2))*(vel(j2)+vel0(j2))-vel0(i2)*vel0(j2))
!
          if(linear.eq.1) then
          emagtd= emagtd + etar*terma*(1+jsym)*                         &
     &    ((phi(i2)+phi0(i2))*(phi(j2)+phi0(j2))-phi0(i2)*phi0(j2))
          emagth= emagth - etar*termb*hypi*(1+jsym)*                    &
     &    ((phi(i2)+phi0(i2))*(phi(j2)+phi0(j2))-phi0(i2)*phi0(j2))
          else
          emagtd= emagtd + etar*terma*(1+jsym)*                         &
     &    (phi(i2)+phi0(i2))*(phi(j2)+phi0(j2))
          emagth= emagth - etar*termb*hypi*(1+jsym)*                    &
     &    (phi(i2)+phi0(i2))*(phi(j2)+phi0(j2))
          endif
!
          endif
          if(numvar.ge.3) then
          j3 = j1 + 12
          i3 = i1 + 12
          if(idens.eq.0) then
            ekin3 = ekin3 - 0.5*terma*(1+jsym)*                         &
     &      ((vel(i3)+vel0(i3))*(vel(j3)+vel0(j3))-vel0(i3)*vel0(j3))
          else
            do k=1,18
            kone = isval1(itri,k)
            ekin3 = ekin3 - .5*g2term(iodd,i,j,k)*dent(kone)*(1+jsym)   &
     &      *((vel(i3)+vel0(i3))*(vel(j3)+vel0(j3))-vel0(i3)*vel0(j3))  &
     &                    + k0term(iodd,k,j,i)*dent(kone)*(1+jsym)      &
     &      *((vel(i1)+vel0(i1))*(vel(j3)+vel0(j3))-vel0(i1)*vel0(j3))
            enddo
          endif
!
!
          ekin3d= ekin3d - amuc*termb*(1+jsym)*                         &
     &    ((vel(i3)+vel0(i3))*(vel(j3)+vel0(j3))-vel0(i3)*vel0(j3))
!
          ekin3h= ekin3h + terma*hypc*amuc*(1+jsym)*                    &
     &    ((com(ione)+com0(ione))*(com(jone)+com0(jone))-               &
     &                com0(ione)*com0(jone))
!
          endif

        enddo   ! end of loop on i
!
            if(numvar.ge.3) then
              gamfac = 1.
              if(gam-1.0 .ne. 0.) gamfac = 1./(gam-1.)
              emag3 = emag3 + gamfac*d2term(iodd,j)*phi(j3)*(1+jsym)
            endif
!
        enddo   ! end of loop on j
      enddo     ! end of loop on irect
      enddo     ! end of loop on jrect
      enddo     ! end of loop on iodd
!    
 8900 format(/," ntime=",i5)
 8901 format(1p5e12.4)
      ekin = ekinp + ekint + ekin3
      emag = emagp + emagt + emag3
      ekind = ekinpd + ekintd + ekin3d
      emagd = emagpd + emagtd + emag3d
!
!
      return
      end
!============================================================
      subroutine rotation(rot,ndim,theta)

      implicit none

      real, intent(out), dimension(ndim, *) ::  rot
      integer, intent(in) :: ndim
      real, intent(in):: theta

      real, dimension(6, 6) :: r1
      real :: co, sn
      integer :: i, j

!
!     calculate the rotation matrix R from Table 2 of Ref 2
!
      co = cos(theta)
      sn = sin(theta)
      do i=1,6
        do j=1,6
          r1(i,j) = 0
        enddo
      enddo
!
      r1(1,1) = 1.
!
      r1(2,2) = co
      r1(2,3) = sn
!
      r1(3,2) = -sn
      r1(3,3) = co
!
      r1(4,4) = co**2
      r1(4,5) = 2.*sn*co
      r1(4,6) = sn**2
!
      r1(5,4) = -sn*co
      r1(5,5) = co**2-sn**2
      r1(5,6) = sn*co
!
      r1(6,4) = sn**2
      r1(6,5) = -2*sn*co
      r1(6,6) = co**2
! 
      do i=1,18
        do j=1,18
          rot(i,j) = 0
        enddo
      enddo
!
      do i=1,6
        do j=1,6
          rot(i,j)       = r1(i,j)
          rot(i+6,j+6)   = r1(i,j)
          rot(i+12,j+12) = r1(i,j)
        enddo
      enddo
!
      return
      end
!============================================================
      subroutine fun(x,z,df,alx,alz)

      implicit none
      
      real, intent(in) :: x, z, alx, alz
      real, intent(out), dimension(0:4, 0:4) :: df

      integer :: i, j
      real :: akx, akz, pi, anorm, snx, cox, snz, coz
      
!
!.....source term and derivatives for initial perturbation
!     first dimension of df is number of derivatives, second is
!     number of x derivatives
!
      pi = acos(-1.)
      anorm = ((pi/alx)**2+(pi/alz)**2)
!
      akx = pi/alx
      akz = pi/alz
!
      snx = sin(akx*x)
      cox = cos(akx*x)
      snz = sin(akz*z)
      coz = cos(akz*z)
!
      do i=0,4
        do j=0,4
          df(i,j) = 0.
        enddo
      enddo
!.....first index is total derivatives, second is x derivative
      df(0,0) = anorm       *snx*snz
!
      df(1,0) = anorm       *akz   *snx*coz
      df(1,1) = anorm*akx          *cox*snz
!
      df(2,0) =-anorm       *akz**2*snx*snz
      df(2,1) = anorm*akx*   akz   *cox*coz
      df(2,2) =-anorm*akx**2       *snx*snz      
!
      df(3,0) =-anorm*       akz**3*snx*coz   
      df(3,1) =-anorm*akx*   akz**2*cox*snz
      df(3,2) =-anorm*akx**2*akz   *snx*coz
      df(3,3) =-anorm*akx**3       *cox*snz
!
      df(4,0) = anorm*       akz**4*snx*snz
      df(4,1) =-anorm*akx*   akz**3*cox*coz
      df(4,2) = anorm*akx**2*akz**2*snx*snz
      df(4,3) =-anorm*akx**3*akz   *cox*coz
      df(4,4) = anorm*akx**4       *snx*snz
!
      return
      end

!============================================================
      subroutine assign(aa,mmnn,i1,j1,val,numvara)
      use p_data
      use t_data
      use basic

      implicit none

      real, intent(out), dimension(9,mmnn,6*numvara,6*numvara) :: aa
      integer, intent(in) :: mmnn, i1, j1, numvara
      real, intent(in) :: val

      integer :: nlocal, ibig, jbig, l, iterm, isl, jsl, itype
!
      nlocal = 6*numvara
      ibig = (i1-1)/nlocal + 1
      jbig = (j1-1)/nlocal + 1
      itype = 0
      do l=1,itypemax(ibig)
         iterm = ibig+incb(ibig,l)
         if(iterm.lt.1 .or. iterm.gt.n*m) then
            write(*,*) l,iterm,ibig
            call safestop(20)
         endif
         if(jbig.eq.ibig+incb(ibig,l)) itype = l
      enddo
      if(itype.eq.0) return
!
      isl = i1 - (ibig-1)*nlocal
      jsl = j1 - (jbig-1)*nlocal
!
      aa(itype,ibig,isl,jsl) = val
!
      return
      end
!============================================================
      subroutine define(slumatrix,aa,mmnn,numvard)
      use p_data
      use t_data
      use basic
#ifdef mpi      
      use supralu_dist_mod
#else
      use supralu_mod
#endif
      implicit none

      real, intent(in), dimension(9,mmnn,6*numvard,6*numvard) :: aa
      integer, intent(in) :: mmnn, numvard

#ifndef mpi
!      type(sparseR8_obj) :: slumatrix

!      call olddefine(slumatrix,aa,mmnn,numvard)
      return
#else
      type(sparseR8d_obj) :: slumatrix

#ifndef BIT64
      integer, parameter :: r8d = selected_real_kind(12,100)
      real(r8d) :: val
#else
      REAL64 :: val
#endif
      integer ::  globsize, nlocal, remainder, firstrow, blocksize
      logical ::  newrow_flag
      integer :: ibig, jbig, il, jl, itype, icol, jcol, irow, ier

!.....Determine local subarray
      globsize = n*m
      nlocal = globsize/maxrank
      remainder = globsize - nlocal*maxrank
      if (myrank .lt. remainder) then
         nlocal = nlocal + 1
         firstrow = myrank*nlocal
      else
         firstrow = myrank*nlocal + remainder
      endif

!.....Outer loop over local rows
      blocksize = 6*numvard
      do ibig=firstrow+1,firstrow+nlocal
         do il=1,blocksize
            newrow_flag = .true.

!...........Inner loop over local columns
            do jbig=1,globsize
               do itype=1,itypemax(ibig)
                  if(incb(ibig,itype) + ibig .eq. jbig) go to 102
               enddo
               cycle
 102           icol = (jbig-1)*blocksize
               do jl=1,blocksize
                  val = aa(itype,ibig,il,jl)
                  call sparseR8d_set_next(slumatrix, icol, newrow_flag,  &
     &                 val, ier)
                  if(idebug.ge.1 .and. islutype .eq. 1) then
                     irow = (ibig-1)*blocksize + il - 1
                     if(myrank.eq.0) write(31,1001) irow,jcol,val
 1001                format(2i10,1pe20.12)
                  endif
                  newrow_flag = .false.
                  icol = icol + 1
               enddo            !jl
            enddo               !jbig
         enddo                  !il
      enddo                     !ibig
!
      if(idebug.ge.1 .and. islutype.eq. 1 .and. myrank.eq.0) then
        close(31)
      endif
!
      return
#endif
      end
!======================================================
      subroutine increment(aa,mmnn,i1,j1,term,numvari,id)
      use p_data
      use t_data
      use basic

      implicit none

      real, dimension(9,mmnn,6*numvari,6*numvari) :: aa
      integer, intent(in) :: mmnn, i1, j1, numvari, id
      real, intent(in), dimension(3,3) :: term
      
      integer :: ibig, isl, itl, iul, jbig, jsl, jtl, jul
      integer :: nlocal, itype, iterm, l

!
      nlocal = 6*numvari
      ibig = (i1-1)/nlocal + 1
      isl = i1 - (ibig-1)*nlocal
      itl = isl + 6 
      iul = itl + 6
      jbig = (j1-1)/nlocal + 1
      jsl = j1 - (jbig-1)*nlocal
      jtl = jsl + 6
      jul = jtl + 6
!
      itype = 0
      do l=1,itypemax(ibig)
         iterm = ibig+incb(ibig,l)
         if(iterm.lt.1 .or. iterm.gt.n*m) then
            write(*,*) l,iterm,ibig
            call safestop(21)
         endif
         if(jbig.eq.ibig+incb(ibig,l)) itype = l
      enddo
!
!     if(isl.eq.1 .and. jsl.eq.1) write(*,9998)ibig,jbig,itype,term(1,1)
!c9998 format(3i5,1pe12.4)
      if(itype.eq.0) then
         write(*,9999) ibig,jbig,i1,j1,isl,jsl
 9999    format("ERROR in increment: ibig,jbig,i1,j1,isl,jsl=",6i5)
         write(*,9998) nlocal,itypemax(ibig),id
 9998    format(" nlocal, itypemax(ibig),id =",3i5)
!     
         call safestop(22)
      endif
!     
      aa(itype,ibig,isl,jsl) = aa(itype,ibig,isl,jsl) + term(1,1)
      if(numvari.eq.1) return
      aa(itype,ibig,isl,jtl) = aa(itype,ibig,isl,jtl) + term(1,2)
      aa(itype,ibig,itl,jsl) = aa(itype,ibig,itl,jsl) + term(2,1)
      aa(itype,ibig,itl,jtl) = aa(itype,ibig,itl,jtl) + term(2,2)
      if(numvari.eq.2) return
      aa(itype,ibig,isl,jul) = aa(itype,ibig,isl,jul) + term(1,3)
      aa(itype,ibig,itl,jul) = aa(itype,ibig,itl,jul) + term(2,3)
      aa(itype,ibig,iul,jsl) = aa(itype,ibig,iul,jsl) + term(3,1)
      aa(itype,ibig,iul,jtl) = aa(itype,ibig,iul,jtl) + term(3,2)
      aa(itype,ibig,iul,jul) = aa(itype,ibig,iul,jul) + term(3,3)
      return
      end
!============================================================
      subroutine inverse(inarray,outarray)
!
!     calculates the inverse of the numvar=1 arrary inarray, and puts
!     it into the numvar=1 array outarray
!
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
!
      real inarray(*),outarray(*)
!
      do lx=1,n
      do lz=1,m
         l = lx + (lz-1)*n
         i0 = 6*(l-1)
         outarray(i0+1) = 1./inarray(i0+1)
         outarray(i0+2) = -inarray(i0+2)/inarray(i0+1)**2
         outarray(i0+3) = -inarray(i0+3)/inarray(i0+1)**2
         outarray(i0+4) = (2.*inarray(i0+2)**2-inarray(i0+1)            &
     &                    *inarray(i0+4))/inarray(i0+1)**3
         outarray(i0+5) = (2.*inarray(i0+2)*inarray(i0+3)               &
     &         -inarray(i0+1)*inarray(i0+5))/inarray(i0+1)**3
         outarray(i0+6) = (2.*inarray(i0+3)**2-inarray(i0+1)            &
     &                    *inarray(i0+6))/inarray(i0+1)**3
      enddo ! on lz
      enddo ! on lx
!
      return
      end
!============================================================
      subroutine boundaryds(ibound,nbc,jsymtype)
      use basic
!
!.....defines the boundary conditions for a rectangular region
      dimension ibound(*)
!
!
!
!.....Define Dirichlet Boundary arrays for the del-squared equation
!
!      jsymtype = 1 for even up-down symmetry
!               = 2 for odd   "        "
      index = 0
!.....start loop on top and bottom boundaries (x-direction)
      if(jper.eq.1) go to 506
      if(jsym.eq.0 .or. jsymtype.eq.2) then
        do lx=1,n
          index = index+2
          ibound(index-1) = 6*(          lx-1) + 1
          ibound(index  ) = 6*((m-1)*n + lx-1) + 1
!.........tangential derivative
          index = index+2
          ibound(index-1) = 6*(          lx-1) + 2
          ibound(index  ) = 6*((m-1)*n + lx-1) + 2
!.........tangential second derivatives
          index = index+2
          ibound(index-1) = 6*(          lx-1) + 4
          ibound(index  ) = 6*((m-1)*n + lx-1) + 4
!.........normal second derivatives
!         index = index+2
!         ibound(index-1) = 6*(          lx-1) + 6
!         ibound(index  ) = 6*((m-1)*n + lx-1) + 6
        enddo
      else
        do lx=1,n
          index = index+2
          ibound(index-1) = 6*(          lx-1) + 3
          ibound(index  ) = 6*((m-1)*n + lx-1) + 1
!.........tangential derivative
          index = index+2
          ibound(index-1) = 6*(          lx-1) + 5
          ibound(index  ) = 6*((m-1)*n + lx-1) + 2
!.........tangential second derivatives
          index = index+1
          ibound(index  ) = 6*((m-1)*n + lx-1) + 4
        enddo
      endif

 506  continue
!.....start loop on left and right boundaries (z-direction)
      if(iper.eq.1) go to 511
      if(m.le.2) go to 501
      if(isym.eq.0) then
        do lz=1,m-2
          index = index+2
          ibound(index-1) = 6*( lz*n      ) + 1
          ibound(index  ) = 6*( lz*n + n-1) + 1
!.........tangential derivative
          index = index+2
          ibound(index-1) = 6*( lz*n      ) + 3
          ibound(index  ) = 6*( lz*n + n-1) + 3
!.........tangential second derivatives
          index = index+2
          ibound(index-1) = 6*( lz*n      ) + 6
          ibound(index  ) = 6*( lz*n + n-1) + 6
        enddo
      else
        do lz=1,m-2
          index = index+2
          ibound(index-1) = 6*( lz*n      ) + 2
          ibound(index  ) = 6*( lz*n + n-1) + 1
!.........tangential derivative
          index = index+2
          ibound(index-1) = 6*( lz*n      ) + 5
          ibound(index  ) = 6*( lz*n + n-1) + 3
!.........tangential second derivatives
          index = index+1
          ibound(index  ) = 6*( lz*n + n-1) + 6
        enddo
      endif
 501  continue
!
!.....corner points
      index = index+3
      ibound(index-2) = 6*( m*n - 1) + 3
      ibound(index-1) = 6*( m*n - 1) + 6
      ibound(index  ) = 6*( m*n - 1) + 5
      if(isym.eq.0) then
        index = index+3
        ibound(index-2) = 6*((m-1)*n ) + 3
        ibound(index-1) = 6*((m-1)*n ) + 6
        ibound(index  ) = 6*((m-1)*n ) + 5
      endif
      if(jsym.eq.0) then
        index = index+3
        ibound(index-2) = 6*(n-1     ) + 3
        ibound(index-1) = 6*(n-1     ) + 6
        ibound(index  ) = 6*(n-1     ) + 5
      endif
      if(isym.eq.0 .and. jsym.eq.0) then
        index = index+3
        ibound(index-2) =                3
        ibound(index-1) =                6
        ibound(index  ) =                5
      endif
!
 511  continue
!
      nbc = index
!
      return
      end
!============================================================
      subroutine scaleback
!
      use p_data
      use t_data
      use basic
      use arrays
!.....scale solution back if a linear calculation
      alpha = 1.
      if(ekin .gt. .001) alpha = sqrt(.001/ekin)
!
      vel = alpha*vel
      phi = alpha*phi
      den = alpha*den
      pres = alpha*pres
!
      alphas = alpha*alpha
      ekin = ekin*alphas
      emag = emag*alphas
      ekind = ekind*alphas
      emagd = emagd*alphas
!
      ekinp = ekinp*alphas
      emagp = emagp*alphas
      ekinpd = ekinpd*alphas
      emagpd = emagpd*alphas
!
!
      ekint = ekint*alphas
      emagt = emagt*alphas
      ekintd = ekintd*alphas
      emagtd = emagtd*alphas
!
      ekinph = ekinph*alphas
      ekinth = ekinth*alphas
      emagph = emagph*alphas
      emagth = emagth*alphas
!
      return
      end
!============================================================
      subroutine evaluate(x,z,ans,ans2,dum,itype,numvare)
      use p_data
      use t_data
      use basic
      use mesh_mod

      implicit none
      real, intent(in) :: x, z
      integer, intent(in) :: itype, numvare
      real, dimension(*), intent(in) :: dum
      real, intent(out) :: ans, ans2
      

!      integer ii, iii, index, k
      integer i, itri, iodd
      real x1, z1, a, b, c, theta, si, eta
      real sum, sum2, term1, term2, small
      real avector(20)!, avector2(20)
!      real wlocal(18)
!      real sumr, rlocal(18),rvector(20)
!     
!.....itype=1 for first velocity (flux) variable
!
!.....evaluate the solution to get the value [ans] at one point (x,z)
!
!
!     first find out what triangle x,z is in
      call whattri(x,z,itri,iodd,x1,z1)
!
!.....calculate local coordinates
      theta = ttri(iodd)
      b = btri(iodd) 
      a = atri(iodd)
      c = ctri(iodd)
      si  = (x-x1)*cos(theta) + (z-z1)*sin(theta) - b
      eta =-(x-x1)*sin(theta) + (z-z1)*cos(theta)
      small = 1.e-6
      if( si.gt.a+small .or. si.lt.-(b+small) .or.                      &
     &    eta .lt.-small .or.eta.gt.c+small) then
      write(*,2222) x,x1,z,z1,si,eta,a,b,c,iodd,itri
 2222 format("error in evaluate",/,1p9e12.4,2i5)
      call safestop(24)
      endif

      call calcavector(itri,iodd,dum,itype,numvare,avector)
!
!.....evaluate the polynomial and second derivative
      sum = 0.
      sum2= 0.
      do i=1,20
         sum = sum + avector(i)*si**mi(i)*eta**ni(i)
!         if(itor.eq.0) then
!.... calculate del-squared
         term1 = 0.
         if(mi(i).ge.2) term1 = mi(i)*(mi(i)-1)*si**(mi(i)-2)*eta**ni(i)
         term2 = 0.
         if(ni(i).ge.2) term2 = ni(i)*(ni(i)-1)*si**mi(i)*eta**(ni(i)-2)
         sum2= sum2 + avector(i)*(term1+term2)
!         else
!.... calculate (1/R) * del-star
!            do k=1,20
!               term1 = 0.
!               if(mi(i)+mi(k).ge.2) term1 = mi(i)*(mi(i)+mi(k)-1)       &
!     &              *si**(mi(i)+mi(k)-2)*eta**(ni(i)+ni(k))
!               term2 = 0.
!               if(ni(i)+ni(k).ge.2) term2 = ni(i)*(ni(i)+ni(k)-1)       &
!     &              *si**(mi(i)+mi(k))*eta**(ni(i)+ni(k)-2)
!               sum2= sum2 + avector(i)*rvector(k)*(term1+term2)
!            enddo
!         endif
      enddo
      ans = sum
      ans2= sum2
!     return
      end
!============================================================
      subroutine second(tcpu)
!     real*4 etime
!     external etime
!     dimension tarray(2)
!     tcpu = etime(tarray)       
#ifndef BIT64
      integer*8 nticks, tickspersec
#else
      integer nticks, tickspersec
#endif
      intrinsic system_clock
      call system_clock(count=nticks, count_rate=tickspersec)

      tcpu = 0.
      if(tickspersec .gt. 0) tcpu = nticks/(1.0*tickspersec)
      return
      end
!============================================================
      subroutine space
!
!     allocates space for big arrays
      use p_data
      use t_data
      use basic
      use arrays
      use superlu

      implicit none

      integer :: mmnn, mmnn18, mmnn18s, mmnn6, nbound, n6

      if(myrank.eq.0 .and. iprint.gt.0) then
         print *, "Entering space."
      endif

      mmnn = m*n
      mmnn18 = m*n*numvar*6
      mmnn18s = (m+1)*n*numvar*6
      mmnn6  = m*n*6
      ntri   = 2*m*(n+1)
      nbound = 12*numvar*(n+m-2)
!
!.....arrays defined at all vertices
      allocate( vel(mmnn18), vels(mmnn18),veln(mmnn18),                 &
     &        velold(mmnn18), vel0(mmnn18), vel1(mmnn18),               &
     &        phi(mmnn18), phis(mmnn18),                                &
     &      phiold(mmnn18s), phi0(mmnn18s), phi1(mmnn18), phip(mmnn18), &
     & b1vector(mmnn18), b2vector(mmnn18),ohmic(mmnn18),viscous(mmnn18), &
     &       b3vector(mmnn18), vtemp(mmnn18),                           &
     &       r4(mmnn18),q4(mmnn18),                                     &
     &       r40(mmnn18),q34(mmnn18),q34temp(mmnn18))
      allocate(r4temp(mmnn18))
!
      allocate(jphi(mmnn6),jphi0(mmnn6), sphi1(mmnn6), sphi2(mmnn6),    &
     &     vor(mmnn6), vor0(mmnn6),  com(mmnn6), com0(mmnn6),           &
     &    den(mmnn6),den0(mmnn6),dent(mmnn6),deni(mmnn6),denold(mmnn6), &
     &    sphip(mmnn6),sphie(mmnn6),pres(mmnn6),pres0(mmnn6),bsi(mmnn6),&
     &    ephi(mmnn6),eph2(mmnn6),eph3(mmnn6),eph4(mmnn6),eph5(mmnn6),   &
     &    eph6(mmnn6),eph7(mmnn6))
!
      if(itor.eq.1) then
         allocate(fun1(mmnn6),fun2(mmnn6),fun3(mmnn6),fun4(mmnn6))
      endif
!
      allocate(bsqr(mmnn6), bsqri(mmnn6))

      if(gyro.eq.1) then
         allocate(galph(mmnn6), ggam(mmnn6), gmu(mmnn6),                &
     &        gw1(mmnn6), gw2(mmnn6), gw3(mmnn6))
      endif
      if(chipar.ne.0) then
         allocate(tw1(mmnn6),tw2(mmnn6),tw3(mmnn6),tw4(mmnn6),          &
     &        tw5(mmnn6), tw40(mmnn6), tw50(mmnn6))
      endif
!
      allocate(incb(mmnn,9),itypemax(mmnn),xcord(mmnn),zcord(mmnn))
!
!.....boundary arrays
      allocate(iboundgs(nbound),iboundv(nbound),iboundp(nbound),        &
     &     psibounds(nbound),velbounds(nbound),combounds(nbound),       &
     &     iboundn(12*(n+m-2)))
!
!.....arrays associated with the triangles
      allocate(ist(ntri,3),isvaln(ntri,18),isval1(ntri,18),             &
     &         isval2(ntri,18),ist0(ntri,3),isval0(ntri,18))
!
!.....metric terms with 1-2 indices
      allocate(aterm(ioddm,18,18), bterm(ioddm,18,18), termbb(18,18),   &
     &     dterm(ioddm,18,18), sterm(ioddm,18,18), d2term(ioddm,18),    &
     &     xterm(ioddm,18,18), yterm(ioddm,18,18))
!
!.....metric terms involving 3 indices
      allocate(k0term (ioddm,18,18,18), k1term (ioddm,18,18,18),        &
     &         k2term (ioddm,18,18,18), g0term (ioddm,18,18,18),        &
     &         g2term (ioddm,18,18,18), g3term (ioddm,18,18,18),        &
     &         g4term (ioddm,18,18,18), g5term (ioddm,18,18,18),        &
     &         g6term (ioddm,18,18,18), g7term (ioddm,18,18,18),        &
     &         g8term (ioddm,18,18,18), g9term (ioddm,18,18,18),        &
     &         g10erm (ioddm,18,18,18), g11erm (ioddm,18,18,18),        &
     &         g12erm (ioddm,18,18,18), g13erm (ioddm,18,18,18),        &
     &         g14erm (ioddm,18,18,18), g15erm (ioddm,18,18,18),        &
     &         g16erm (ioddm,18,18,18),                                 &
     &         h3term (ioddm,18,18,18), h5term (ioddm,18,18,18),        &
     &         x0term (ioddm,18,18,18), y0term (ioddm,18,18,18),        &
     &         x1term (ioddm,18,18,18), y1term (ioddm,18,18,18),        &
     &         x2term (ioddm,18,18,18), y2term (ioddm,18,18,18))
!
!     local metric terms involving 2 indices
      allocate(u1terml(18,18), u2terml(18,18), u4terml(18,18),          &
     &        u3terml(18,18), u6terml(18,18), u5terml(18,18),           &
     &        u7terml(18,18), u8terml(18,18), c10erml(18,18),           &
     &        c11erml(18,18), c12erm2(18,18), c13erm2(18,18),           &
     &        c14erm2(18,18), c15erm2(18,18), c16erm2(18,18),           &
     &        v1terml(18,18), v2terml(18,18), v3terml(18,18),           &
     &        v4terml(18,18), v5terml(18,18), v6terml(18,18),           &
     &        v7terml(18,18), v8terml(18,18), v9terml(18,18),           &
     &        v10erml(18,18), v11erml(18,18), v12erml(18,18),           &
     &        v13erml(18,18), v14erml(18,18), v15erml(18,18),           &
     &        v16erml(18,18),                                           &
     &        q0terml(18,18), c12erml(18,18),                           &
     &        c13erml(18,18), c14erml(18,18), c15erml(18,18),           &
     &        c16erml(18,18), c1term2(18,18), c3term2(18,18),           &
     &        c4term2(18,18), c5term2(18,18), c6term2(18,18),           &
     &        c7term2(18,18), c8term2(18,18), c9term2(18,18),           &
     &        c10erm2(18,18), c11erm2(18,18))
      if(gyro.eq.1) then
         allocate(                                                      &
     &        gv11arml(18,18), gv12arml(18,18), gv12brml(18,18),        &
     &        gv12crml(18,18), gv13arml(18,18), gv13brml(18,18),        &
     &        gv13crml(18,18), gv22arml(18,18), gv23arml(18,18),        &
     &        gv23brml(18,18), gv33brml(18,18), gv33crml(18,18),        &
     &        gv12axml(18,18), gv12bxml(18,18), gv12cxml(18,18),        &
     &        gv13bxml(18,18), gv13cxml(18,18), gv23axml(18,18),        &
     &        gv23bxml(18,18))
      endif
      if(chipar.ne.0) then
         allocate(                                                      &
     &        t1terml(18,18), t2terml(18,18), t3terml(18,18),           &
     &        t4terml(18,18), t5terml(18,18))
      endif
!
!.....metric terms involved in toroidal geometry
      if(itor.eq.1) then
      allocate(hterm(ntridim,18,18),termh2(ioddm,18,18,20),rinv(mmnn6))
      endif
!
!     intermediate metric terms with 1-3 indices                        
      allocate(d2matrix(20),     bmatrix(20,20),  amatrix(20,20),       &
     &        dmatrix(20,20),  qmatrix(20,20),                          &
     &       xmatrix(20,20),   ymatrix(20,20),                          &
     &       termg0(20,20,20), termg2(20,20,20), termg3(20,20,20),      &
     &       termg4(20,20,20), termg5(20,20,20), termg6(20,20,20),      &
     &       termg7(20,20,20), termg8(20,20,20), termg9(20,20,20),      &
     &       terg10(20,20,20), terg11(20,20,20), terg12(20,20,20),      &
     &       terg13(20,20,20), terg14(20,20,20), terg15(20,20,20),      &
     &       terg16(20,20,20), termh(20,20,20),  termh3(20,20,20),      &
     &       termh5(20,20,20), termk0(20,20,20), termk1(20,20,20),      &
     &       termk2(20,20,20), termx0(20,20,20), termy0(20,20,20),      &
     &       termx1(20,20,20), termy1(20,20,20), termx2(20,20,20),      &
     &       termy2(20,20,20))
!
!
!     these contain the data that is given to superlu
!
      n6=max(numvar,2)*6
!      allocate(ss(9,mmnn,n6,n6),dd(9,mmnn,n6,n6), rr(9,mmnn,n6,n6),     &
!     &       qq(9,mmnn,n6,n6),sstemp(9,mmnn,n6,n6))
      allocate(sstemp(9,mmnn,n6,n6))
      allocate(ss(9,mmnn,n6,n6))

      if(myrank.eq.0 .and. iprint.gt.0) then
         print *, "Exiting space."
      endif

      return
      end
!==============================================================
      subroutine safestop(iarg)
#ifdef mpi
      use supralu_dist_mod
      include 'mpif.h'
      integer ier
      call SLUD_exit
      call MPI_Finalize(ier)
      if (ier /= 0) print *,'Error terminating MPI:',ier
#endif
      close(UNIT=9, STATUS='KEEP')
      write(*,*) "stopped at", iarg
      stop
      end
!==============================================================
      subroutine smoother1(inarray,outarray,mmnn,numvard,iplace)
!
!     This routine takes as input the Laplacian of the scalar field in
!     the array inarray, smooths it, and returns the smoothed scalar array
!     in location iplace in the stand-alone array outarray(mmnn6).  The LU
!     decomposition takes place the first-time only.
!
!
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
      dimension ssterm(3,3)
      real  outarray(mmnn*6*numvard),inarray(mmnn*6)
#ifndef BIT64
      integer, parameter :: r9a = selected_real_kind(12,100)
      real(r9a), allocatable:: temp(:)
#else
      REAL64, allocatable :: temp(:)
#endif
      integer,  save, allocatable:: iboundsm(:)
      integer, save :: nbcsm
!
      mmnn6 = 6*mmnn
      allocate(temp(2*mmnn6))
      temp = 0
      hyp = hyperc*deex**2*amu
      thimpsm = 1.
! 
      if(ifirsts5_lu.ne.0) go to 500
      ifirsts5_lu = 1
      allocate(iboundsm(24*(n+m-2)))
!
      numvarsm = 2
      nrasm = n*m*6*numvarsm
      msizesm = 9*n*m*(6*numvarsm)**2
      call jbdecomp1(n*m, 6*numvarsm, m_lsm, nnz_lsm, ifs_rsm)
      base = 0
!
!.....compute LU decomposition only once
!
!.....form matrix
      ss = 0
      itri = 0
      do iodd=1,ioddmx
      do ll=1,nreg
        itri = itri + 1
        do i=1,18
        do j=1,18
            jtwo = isval2(itri,j)
            itwo = isval2(itri,i)
            ssterm(1,1) = dterm(iodd,i,j)
            ssterm(1,2) =-aterm(iodd,i,j)
            ssterm(2,1) =  +dt*hyp*thimpsm*bterm(iodd,i,j)
            ssterm(2,2) = aterm(iodd,i,j)
            call increment(ss,m*n,itwo,jtwo,ssterm,numvarsm,10)
!
        enddo
        enddo
      enddo  !  on nreg
      enddo  !  on iodd
!
!......define indices for boundary arrays
      call boundarysm(iboundsm,nbcsm,iplace)
!
!.....modify the s-matrix, inserting the boundary conditions
      do l=1,nbcsm
        i = iboundsm(l)
        do j= 1,nrasm
            if(i.le.0 .or. i.ge.m*n*6*numvarsm+1) then
            write(*,3343) i,l,iboundsm(l)
 3343       format(" error in iboundsm index",3i5)
            call safestop(3343)
          endif
         call assign(ss,n*m,i,j,0.,numvarsm)
        enddo
      enddo
      do l=1,nbcsm
         call assign(ss,n*m,iboundsm(l),iboundsm(l),1.,numvarsm)
      enddo
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "before init"
!  
#ifdef mpi
      call sparseR8d_init(s5matrix_lu,nrasm,nnz_lsm,m_lsm,ifs_rsm,ier)
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "before 2nd init"
#else
      call sparseR8_init(s5matrix_lu,nrasm,msizesm,base,ier)
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "between inits"
      s5matrix_lu %permc_spec = 0
#endif
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "after init"
      call define(s5matrix_lu,ss,m*n,numvarsm)
!
!.....perform LU decomposition of sparse matrix "s5matrix"
!     store result in opaque object "s5handle"
      jer = 0
#ifdef mpi
       call sparseR8d_new(s5matrix_lu, jer)
       if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "after new"
       if(jer.ne.0) then
         write(*,*) 'after 3rd sparseR8d_new', jer
         call safestop(46)
       endif
#else
       call dsupralu_new(s5handle, s5matrix_lu%values(1),               &
     &      s5matrix_lu%irow(1), s5matrix_lu%jcol_ptr(1),               &
     &      s5matrix_lu%nnz, nrasm , jer)
       if(jer.ne.0) then
         write(*,*) 'after dsupralu_new', jer
         call safestop(47)
       endif
       if(myrank.eq.0 .and.iprint.ge.1) write(*,*) "before _colperm-s3"
       call dsupralu_colperm(s5handle, s5matrix_lu% permc_spec, jer)
       if(myrank.eq.0 .and.iprint.ge.1) write(*,*) "after  _colperm-s3"
       if(jer.ne.0) then
         write(*,*) 'after dsupralu_new', jer
         call safestop(48)
       endif
       call dsupralu_lu(s5handle, jer)
       if(jer.ne.0) then
         write(*,*) 'after dsupralu_lu', jer
         call safestop(49)
       endif
#endif
!
 500   continue
!
!.....define RHS vector
      itri = 0
      do iodd=1,ioddmx
      do ll=1,nreg
        itri = itri + 1
        do i=1,18
        do j=1,18
            jone = isval1(itri,j)
            itwo = isval2(itri,i)
!
            temp(itwo+6) = temp(itwo+6)+(dterm(iodd,i,j)                &
     &             - dt*hyp*(1.-thimpsm)*bterm(iodd,i,j))*inarray(jone)
!
        enddo
        enddo
      enddo
      enddo
      do l=1,nbcsm
      if(iboundsm(l).le.0 .or. iboundsm(l).gt.mmnn6*numvarsm) then
        write(*,*) "error in smoother 1"
        call safestop(98)
        endif
      temp(iboundsm(l)) = 0.
      enddo
!
!.....perform LU backsubstitution to get outarray solution
      if(myrank.eq.0) call second(tbeforesolve)

#ifdef mpi
      call sparseR8d_solve(s5matrix_lu,temp,ier)
#else
      call dsupralu_solve(s5handle,temp,ier)
#endif
      if(myrank.eq.0) call second(taftersolve)
      tsolve = tsolve + taftersolve-tbeforesolve
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "solve 3"
!
!.....store solution
      do lx=1,n
      do lz=1,m
         l = lx + (lz-1)*n
         do i=1,6
         outarray(6*((l-1)*numvard +(iplace-1))+i)                      &
     &     = temp(6*((l-1)*numvarsm+        1 )+i)
         enddo
      enddo
      enddo
!
      deallocate(temp)
      return
!
      end
!============================================================
      subroutine smoother3(inarray,outarray,mmnn,numvard,iplace)
!
!     This routine takes as input the Laplacian of the scalar field in
!     the array inarray, smooths it, and returns the smoothed scalar array
!     in location iplace in the stand-alone array outarray(mmnn6).  The LU
!     decomposition takes place the first-time only.
!
!
      use p_data
      use t_data
      use basic
      use arrays
      use superlu

      implicit none

      real, dimension(mmnn*6*numvard), intent(out) :: outarray
      real, dimension(mmnn*6), intent(in) :: inarray
      integer, intent(in) :: numvard, iplace, mmnn

      real, dimension(3,3) :: ssterm
#ifndef BIT64
      integer, parameter :: r9a = selected_real_kind(12,100)
      real(r9a), allocatable:: temp(:)
#else
      REAL64, allocatable :: temp(:)
#endif
      integer,  save, allocatable:: iboundsm(:)
      integer, save :: nbcsm
      integer :: itri, iodd, l, lx, lz, ier, jer, jone, i, j, ll
      integer :: itwo, jtwo, nrasm, mmnn6
      integer :: m_lsm, nnz_lsm, ifs_rsm, msizesm, numvarsm
      real :: tafterinit, tbeforeinit, tbeforesolve, taftersolve
      real :: hyp, thimpsm
!
      mmnn6 = 6*mmnn
      allocate(temp(2*mmnn6))
      temp = 0
      hyp = hyperc*deex**2*amuc
      thimpsm = 1.
      numvarsm = 2
! 
      if(ifirsts7_lu.ne.0) go to 500
      ifirsts7_lu = 1
      allocate(iboundsm(24*(n+m-2)))
!
      nrasm = n*m*6*numvarsm
      msizesm = 9*n*m*(6*numvarsm)**2
      call jbdecomp1(n*m, 6*numvarsm, m_lsm, nnz_lsm, ifs_rsm)
      base = 0
!
!.....compute LU decomposition only once
!
!.....form matrix
      ss = 0
      itri = 0
      do iodd=1,ioddmx
      do ll=1,nreg
        itri = itri + 1
        do i=1,18
        do j=1,18
            jtwo = isval2(itri,j)
            itwo = isval2(itri,i)
            ssterm(1,1) = dterm(iodd,i,j)
            ssterm(1,2) =-aterm(iodd,i,j)
            ssterm(2,1) =  +dt*hyp*thimpsm*bterm(iodd,i,j)
            ssterm(2,2) = aterm(iodd,i,j) - regular*dterm(iodd,i,j)
            call increment(ss,m*n,itwo,jtwo,ssterm,numvarsm,10)
!
        enddo
        enddo
      enddo
      enddo
!
!......define indices for boundary arrays
      call boundarysm(iboundsm,nbcsm,iplace)
!
!.....modify the s-matrix, inserting the boundary conditions
      do l=1,nbcsm
        i = iboundsm(l)
        do j= 1,nrasm
            if(i.le.0 .or. i.ge.m*n*6*numvarsm+1) then
            write(*,3343) i,l,iboundsm(l)
 3343       format(" error in iboundsm index",3i5)
            call safestop(3343)
          endif
         call assign(ss,n*m,i,j,0.,numvarsm)
        enddo
      enddo
      do l=1,nbcsm
         call assign(ss,n*m,iboundsm(l),iboundsm(l),1.,numvarsm)
      enddo
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "before init"
      if(myrank.eq.0) call second(tbeforeinit)
!  
#ifdef mpi
      call sparseR8d_init(s7matrix_lu,nrasm,nnz_lsm,m_lsm,ifs_rsm,ier)
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "before 2nd init"
#else
      call sparseR8_init(s7matrix_lu,nrasm,msizesm,base,ier)
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "between inits"
      s7matrix_lu %permc_spec = 0
#endif
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "after init"
      call define(s7matrix_lu,ss,m*n,numvarsm)
      if(myrank.eq.0) call second(tafterinit)
!      tinit = tinit + tafterinit-tbeforeinit
!
!.....perform LU decomposition of sparse matrix "s7matrix"
!     store result in opaque object "s7handle"
      jer = 0
#ifdef mpi
       call sparseR8d_new(s7matrix_lu, jer)
       if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "after new"
       if(jer.ne.0) then
         write(*,*) 'after 3rd sparseR8d_new', jer
         call safestop(46)
       endif
#else
       call dsupralu_new(s7handle, s7matrix_lu%values(1),               &
     &      s7matrix_lu%irow(1), s7matrix_lu%jcol_ptr(1),               &
     &      s7matrix_lu%nnz, nrasm , jer)
       if(jer.ne.0) then
         write(*,*) 'after dsupralu_new', jer
         call safestop(47)
       endif
       if(myrank.eq.0 .and.iprint.ge.1) write(*,*) "before _colperm-s3"
       call dsupralu_colperm(s7handle, s7matrix_lu% permc_spec, jer)
       if(myrank.eq.0 .and.iprint.ge.1) write(*,*) "after  _colperm-s3"
       if(jer.ne.0) then
         write(*,*) 'after dsupralu_new', jer
         call safestop(48)
       endif
       call dsupralu_lu(s7handle, jer)
       if(jer.ne.0) then
         write(*,*) 'after dsupralu_lu', jer
         call safestop(49)
       endif
#endif
!
 500   continue
!
!.....define RHS vector
      itri = 0
      do iodd=1,ioddmx
      do ll=1,nreg
        itri = itri + 1
        do i=1,18
        do j=1,18
            jone = isval1(itri,j)
            itwo = isval2(itri,i)
!
            temp(itwo+6) = temp(itwo+6)+(dterm(iodd,i,j)                &
     &             - dt*hyp*(1.-thimpsm)*bterm(iodd,i,j))*inarray(jone)
!
        enddo
        enddo
      enddo
      enddo
!!.....added 5/26/06 SCJ....(should remove combounds altogether)
      combounds = 0
      do l=1,nbcsm
      if(iboundsm(l).le.0 .or. iboundsm(l).gt.mmnn6*numvarsm) then
        write(*,*) "error in smoother 3"
        call safestop(98)
      endif
      temp(iboundsm(l)) = combounds(l)
      enddo
!
!.....perform LU backsubstitution to get outarray solution

      if(myrank.eq.0) call second(tbeforesolve)
#ifdef mpi
      call sparseR8d_solve(s7matrix_lu,temp,ier)
#else
      call dsupralu_solve(s7handle,temp,ier)
#endif
      if(myrank.eq.0) call second(taftersolve)
      tsolve = tsolve + taftersolve-tbeforesolve
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "after solve 4"
!
!.....store solution
      do lx=1,n
      do lz=1,m
         l = lx + (lz-1)*n
         do i=1,6
         outarray(6*((l-1)*numvard +(iplace-1))+i)                      &
     &     = temp(6*((l-1)*numvarsm+        1 )+i)
         enddo
      enddo
      enddo
!
      deallocate(temp)
      return
!
      end
!============================================================
      subroutine boundarysm(ibound,nbc,iplace)
      use basic
      use arrays
!
!.....defines the velocity boundary conditions for a rectangular region
      dimension ibound(*)
      velbounds = 0
      combounds = 0
      numvarsm = 2
      index = 0
      if(iplace.ge.3) go to 600
!
!.....FIRST EQUATION: ..... for vorticity
      if(jper.eq.1) go to 506
!.....start loop on top and bottom boundaries (y-direction)
      do lx=1,n
!
      index = index+1
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 1
!.....tangential derivative (x)
      index = index+1
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 2
!.....tangential second derivatives (xx)
      index = index+1
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 4
!     if(jsym.eq.0) then
        index = index+1
        ibound(index  ) = 6*numvarsm*(          lx-1) + 1
!.......tangential derivative (x)
        index = index+1
        ibound(index  ) = 6*numvarsm*(          lx-1) + 2
!.......tangential second derivatives (xx)
        index = index+1
        ibound(index  ) = 6*numvarsm*(          lx-1) + 4
!     endif
!
      enddo
  506 continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 511
      if(m.le.2) go to 501
      do lz=1,m-2
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 1
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 1
!.....tangential derivative (y)
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 3
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 3
!.....tangential second derivatives (yy)
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 6
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 6
      enddo
 501  continue
!.....corner points
      index = index+4
      ibound(index-3) =                 6
      ibound(index-2) = 6*numvarsm*(n-1     ) + 6
      ibound(index-1) = 6*numvarsm*((m-1)*n ) + 6
      ibound(index  ) = 6*numvarsm*( m*n - 1) + 6
      index = index+4
      ibound(index-3) =                 3
      ibound(index-2) = 6*numvarsm*(n-1     ) + 3
      ibound(index-1) = 6*numvarsm*((m-1)*n ) + 3
      ibound(index  ) = 6*numvarsm*( m*n - 1) + 3
 511  continue
      if(numvarsm.le.1) then
        nbc = index
        return
      endif
!.....SECOND EQUATION: ... potential
!.....start loop on top and bottom boundaries (x-direction)
      if(jper.ge.1) go to 516
      do lx=1,n
      index = index+2
      ibound(index-1) = 6*numvarsm*(          lx-1) + 7
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 7
!.....tangential derivative (x)
      index = index+2
      ibound(index-1) = 6*numvarsm*(          lx-1) + 8
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 8
!.....tangential second derivative (xx)
      index = index+2
      ibound(index-1) = 6*numvarsm*(          lx-1) + 10
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 10
!.....normal second derivative (yy)
!     if(jsym.eq.0) then
        index = index+1
        ibound(index  ) = 6*numvarsm*(          lx-1) + 12
!     endif
      index = index+1
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 12
      enddo
  516 continue

!.....start loop on left and right boundaries (z-direction)
      if(iper.eq.1) go to 512
      if(m.le.2) go to 502
      do lz=1,m-2
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 7
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 7
!.....tangential derivative (y)
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 9
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 9
!.....tangential second derivative (yy)
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 12
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 12
!.....normal second derivative (xx)
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 10
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 10
      enddo
 502  continue
!.....corner points
      index = index+4
      ibound(index-3) =                       9
      ibound(index-2) = 6*numvarsm*(n-1     ) + 9
      ibound(index-1) = 6*numvarsm*((m-1)*n ) + 9
      ibound(index  ) = 6*numvarsm*( m*n - 1) + 9
 512  continue
      nbc = index
      return
!
!.....start of boundary conditions for smoother3
 600  continue
      combounds = 0
      if(jper.ge.1) go to 526
!.....start loop on top and bottom boundaries (y-direction)
      do lx=1,n
!
!.....FIRST EQUATION...del squared chi
!.....top boundary
      index = index+1
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 1
      combounds(index  ) = com(6*((m-1)*n + lx-1) + 1)
!.....tangential derivative (x)
      index = index+1
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 2
      combounds(index  ) = com(6*((m-1)*n + lx-1) + 2)
!.....tangential second derivatives (xx)
      index = index+1
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 4
      combounds(index  ) = com(6*((m-1)*n + lx-1) + 4)
!
!.....bottom boundary
      if(jsym.eq.0) then
        index = index+1
        ibound(index) = 6*numvarsm*(          lx-1) + 1
        combounds(index) = com(6*(          lx-1) + 1)
!.......tangential derivative (x)
        index = index+1
        ibound(index) = 6*numvarsm*(          lx-1) + 2
        combounds(index) = com(6*(          lx-1) + 2)
!.......tangential second derivatives (xx)
        index = index+1
        ibound(index) = 6*numvarsm*(          lx-1) + 4
        combounds(index) = com(6*(          lx-1) + 4)
      else
!
!.......normal derivative (y)
!        index = index+1
!        ibound(index) = 6*numvarsm*(          lx-1) + 3
!        combounds(index) = 0
!.......cross second derivatives (xy)
!        index = index+1
!        ibound(index) = 6*numvarsm*(          lx-1) + 5
!        combounds(index) = 0
      endif
!

      enddo
  526 continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 611
      if(m.le.2) go to 601
      do lz=1,m-2
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 1
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 1
      combounds(index-1) = com(6*( lz*n      ) + 1)
      combounds(index  ) = com(6*( lz*n + n-1) + 1)
!.....tangential derivative (y)
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 3
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 3
      combounds(index-1) = com(6*( lz*n      ) + 3)
      combounds(index  ) = com(6*( lz*n + n-1) + 3)
!.....tangential second derivatives (yy)
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 6
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 6
      combounds(index-1) = com(6*( lz*n      ) + 6)
      combounds(index  ) = com(6*( lz*n + n-1) + 6)
      enddo
 601  continue
!.....corner points
      index = index+4
      ibound(index-3) =                 6
      ibound(index-2) = 6*numvarsm*(n-1     ) + 6
      ibound(index-1) = 6*numvarsm*((m-1)*n ) + 6
      ibound(index  ) = 6*numvarsm*( m*n - 1) + 6
      combounds(index-3) = com(6)
      combounds(index-2) = com(6*(n-1     ) + 6)
      combounds(index-1) = com(6*((m-1)*n ) + 6)
      combounds(index  ) = com(6*( m*n - 1) + 6)
      index = index+4
      ibound(index-3) =                 3
      ibound(index-2) = 6*numvarsm*(n-1     ) + 3
      ibound(index-1) = 6*numvarsm*((m-1)*n ) + 3
      ibound(index  ) = 6*numvarsm*( m*n - 1) + 3
      combounds(index-3) = com(3)
      combounds(index-2) = com(6*(n-1     ) + 3)
      combounds(index-1) = com(6*((m-1)*n ) + 3)
      combounds(index  ) = com(6*( m*n - 1) + 3)
 611  continue
      if(numvarsm.le.1) then
        nbc = index
        return
      endif
!.....SECOND EQUATION .... Chi
      if(jper.ge.1) go to 536
      do lx=1,n
!.....normal derivative
      index = index+2
      ibound(index-1) = 6*numvarsm*(          lx-1) + 9
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 9
!.....cross derivative
      index = index+2
      ibound(index-1) = 6*numvarsm*(          lx-1) + 11
      ibound(index  ) = 6*numvarsm*((m-1)*n + lx-1) + 11
      enddo
  536 continue
!.....start loop on left and right boundaries (x-direction)
      if(iper.eq.1) go to 513
      if(m.le.2) go to 503
      do lz=1,m-2
!.....normal derivative
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 8
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 8
!.....cross derivative
      index = index+2
      ibound(index-1) = 6*numvarsm*( lz*n      ) + 11
      ibound(index  ) = 6*numvarsm*( lz*n + n-1) + 11
      enddo
 503  continue
!.....corner points
      index = index+4
      ibound(index-3) =                 8
      ibound(index-2) = 6*numvarsm*(n-1     ) + 8
      ibound(index-1) = 6*numvarsm*((m-1)*n ) + 8
      ibound(index  ) = 6*numvarsm*( m*n - 1) + 8
 513  continue
!     index = index + 1
!     ibound(index) = 6*numvarsm*(m*n-2) + 7
      nbc = index
      return
      end
      subroutine MPE_DECOMP1D(nproc,mrows,myrank,jfirst,jlast)
      ilast = 0
      ifirst = 1
      mleft = mrows
      nleft = nproc + 1
      do i=0,myrank
        mleft = mleft - (ilast-ifirst+1)
        nleft = nleft - 1
        ifirst = ilast + 1
        ifrac = mleft/nleft
        if(nleft*ifrac.eq.mleft) then
          ilast = ifirst + ifrac - 1
        else
          ilast = ifirst + ifrac 
        endif
      enddo
      jfirst = ifirst
      jlast  = ilast
      return
      end
!===========================================================
      subroutine errorcheck
      use p_data
      use t_data
      use basic
      use arrays
!
      dimension array(irs,3,4),errorp(3),errorv(3),                     &
     &          ave(3,4),errorpn(3),errorvn(3),xval(irs)
!
      ave = 0.
      errorp = 0.
      errorv = 0.
      errorpn = 0.
      errorvn = 0.
      do ii=1,numvar
!
        do ix=1,irs
          x = (ix-1)*alx*(1./(1.+isym))/(irs-1)
          xval(ix) = x
!
          call evaluate(x,alz/2.,array(ix,ii,1),dum,phi,ii,numvar)
          call evaluate(x,alz/2.,array(ix,ii,2),dum,phis,ii,numvar)
          call evaluate(x,alz/2.,array(ix,ii,3),dum,vel,ii,numvar)
          call evaluate(x,alz/2.,array(ix,ii,4),dum,vels,ii,numvar)
!
          do i=1,4
            ave(ii,i) = ave(ii,i) + array(ix,ii,i)
          enddo
!
        enddo ! over irs
          do i=1,4
            ave(ii,i) = ave(ii,i)/irs
          enddo
        do ix=1,irs
          errorp(ii) = errorp(ii) + (array(ix,ii,1)-ave(ii,1)            &
     &                             -(array(ix,ii,2)-ave(ii,2)))**2
          errorv(ii) = errorv(ii) + (array(ix,ii,3)-ave(ii,3)            &
     &                             -(array(ix,ii,4)-ave(ii,4)))**2
        enddo
!
!
!.......calculate root mean square error
        errorp(ii)= sqrt(errorp(ii)/irs)
        errorv(ii)= sqrt(errorv(ii)/irs)
!
!.......calculate normalization
        do ix=1,irs
          errorpn(ii) = errorpn(ii)+ abs(array(ix,ii,2)-ave(ii,2)) 
          errorvn(ii) = errorvn(ii)+ abs(array(ix,ii,4)-ave(ii,4)) 
        enddo ! over irs
        errorpn(ii)= (errorpn(ii)/irs)
        errorvn(ii)= (errorvn(ii)/irs)
!
!       normalize errors
        errorp(ii) = errorp(ii)/errorpn(ii)
        errorv(ii) = errorv(ii)/errorvn(ii)
!

      enddo ! over ii
!
      write(*,1000) 
      do ii=1,3
         write(*,1001) ii,errorp(ii),errorv(ii),errorpn(ii),errorvn(ii)
         write(9,1001) ii,errorp(ii),errorv(ii),errorpn(ii),errorvn(ii)
      enddo
      do ii=1,3
         do i=1,4
            write(50,2001) ii,i
 2001       format("ii,i = ",2i3)
            write(50,2002) (array(ix,ii,i),ix=1,irs)
 2002       format(1p9e10.2)
         enddo                  ! on i
      enddo                     ! on ii
!
 1000 format(" ii  errorp      errorv      errorpn      ",
     &       "errorvn")
 1001 format(i3,1p4e12.4)
      return
      end
!============================================================
      subroutine jbdecomp1(meshpoints, blocksize, localrows, localnz,   &
     &     firstrow)
      use basic
      implicit none
!
      integer, intent(in)  :: meshpoints, blocksize
      integer, intent(out) :: localrows, localnz, firstrow
!
      integer localpoints, remainder, firstpoint
!
      localpoints = meshpoints / maxrank
      remainder = meshpoints - localpoints*maxrank
      if (myrank .lt. remainder) then
         localpoints = localpoints + 1
         firstpoint = myrank*localpoints
      else
         firstpoint = myrank*localpoints + remainder
      endif
!
      localrows = localpoints*blocksize
      localnz = 9*localpoints*(blocksize**2)
      firstrow = firstpoint*blocksize
      if(iprint.eq.1) then
         print *,'jbdecomp:',myrank,meshpoints,blocksize,localrows,     &
     &        localnz, firstrow
      endif
      end
!============================================================================
      subroutine wave_perturbation(per1,per2,per3, itype)
      use basic
      implicit none
!
!     Calculates the magnitudes of the initial perturbations for
!     the wave propagation test.  Puts the perturbation in dum
!
!     itype | field
!     --------------
!     1     | psi, bz, pe
!     2     | phi, vz, chi
!     3     | n
!     4     | p
!     5     | diagnostic output
!
      real, intent(out) :: per1, per2, per3
      integer, intent(in) :: itype
      real akx,akx2,b2,a2
      real kp,km,t1,t2,t3
      real omega
      real coef(4)
      real root(3)
      real error(3)
      real psiper, phiper, bzper, vzper, peper, chiper, nper, pper
      real bi
!
      akx = 2.*pi/alx
!
!     for itaylorw=3, set up a phi perturbation only
      if(itaylorw.eq.3) then
         phiper = eps
         vzper = 0.
         chiper = 0.
         psiper = 0.
         bzper = 0.
         nper = 0.
         pper = 0.
         peper = 0.
         goto 1
      endif
!
      akx2 = akx**2
      b2 = bzero*bzero + del*del
      a2 = del**2/b2
!
      if(gyro.eq.1) then
         bi = 2.*pi0
      else 
         bi = 0.0
      endif
!
!.....numvar=2 =================
      if(numvar.eq.2) then
         kp = akx * (b2 + bi*(3.*a2-1.)/4.)
         km = akx * (b2 - bi*(3.*a2-1.)/4.)
!........fast wave
         if(itaylorw.eq.0) then         
            omega = (akx/2.)*sqrt(a2/b2)*(sqrt(4.*b2**2+km**2)+kp)
            psiper = eps
            phiper = -eps*2.*b2 / (sqrt(4.*b2**2+km**2)+km)
            bzper = akx*psiper
            vzper = akx*phiper
!........slow wave
         else
            omega = (akx/2.)*sqrt(a2/b2)*(sqrt(4.*b2**2+km**2)-kp)
            psiper = eps
            phiper = -eps*2.*b2 / (sqrt(4.*b2**2+km**2)-km)
            bzper = -akx*psiper
            vzper = -akx*phiper
         endif
!
!.....numvar=3 =================
      elseif(numvar.ge.3) then
!
         coef(4) = 96.*b2**4
         coef(3) = -6.*akx2*b2**3*(16.*b2**2*(1.+a2+akx2*a2)            &
     &        + akx2*bi**2*(1.+6.*a2-3.*a2**2)                          &
     &        + 16.*gam*p0*b2)
         coef(2) = 6.*a2*akx2**2*b2**3                                  &
     &        *(b2*((4.*b2-bi*akx2*(1.+a2))**2                          &
     &        +4.*bi**2*akx2*(1.-a2)*(1.+akx2*a2))                      &
     &        + gam*p0*(32.*b2**2+16.*akx2*b2**2                        &
     &        +bi**2*akx2*(1.-3.*a2)**2))
         coef(1) = -6.*gam*p0*a2**2*akx2**3*b2**4                       &
     &        *(4.*b2+akx2*bi*(1.-3.*a2))**2
!
         if((itype.eq.5).and.(myrank.eq.0)) then
            write(*,*) "Coefs: ", coef(1), coef(2), coef(3), coef(4)
         endif
!
         call cubic_roots(coef, root, error)
!
         if((itype.eq.5).and.(myrank.eq.0)) then
            write(*,*) "Coefs: ", coef(1), coef(2), coef(3), coef(4)
            write(*,*) "Roots: ", root(1), root(2), root(3)
            write(*,*) "Error: ", error(1), error(2), error(3)
         endif
!
         omega = sqrt(root(3))
         if(itaylorw.eq.1) omega=sqrt(root(1))
         if(itaylorw.eq.2) omega=sqrt(root(2))
!
         t1 = akx2*bi/(4.*b2)
         t2 = (akx2/omega**2)*                                          &
     &        (bzero**2/(1.-gam*p0*akx2/(omega**2)))
         t3 = del * akx / omega
!
         psiper = eps
         phiper = psiper*(akx2*t3*(1.+(t3**2/(1.-t2-t3**2)))-(1.-t2)/t3)&
     &        / (1. - t2 + t1*t2*(1.+3.*a2)                             &
     &        + t1*t3**2*(2.*t2-(1.-3.*a2))/(1-t2-t3**2))
         vzper = phiper*t1*t3*(2.*t2-(1.-3.*a2))/(1-t2-t3**2)           &
     &        - psiper*akx2*t3**2/(1.-t2-t3**2)
         bzper = (-phiper*t1*t2*(1.+3.*a2) - vzper*t3 + psiper*akx2*t3) &
     &        / (1.-t2)
         chiper = bzero / (1.-gam*p0*akx2/(omega**2)) / omega           &
     &        * ((1.+3.*a2)*t1*phiper - bzper)
         peper = -chiper*gam*(p0-pi0)*akx2 / omega
         nper = -chiper*akx2/omega
         pper = -chiper*gam*p0*akx2 / omega
      endif
!
 1    goto (10,20,30,40,2), itype
      goto 50
!
!.....diagnosic output
 2    continue
      write(*,1020) omega,akx,omega/akx
      write(*,*) "psi1 = ", psiper
      write(*,*) "phi1 = ", phiper
      write(*,*) "vz1 = ", vzper
      write(*,*) "bz1 = ", bzper
      write(*,*) "pe1 = ", peper
      write(*,*) "chi1 = ", chiper
      write(*,*) "n1 = ", nper
      write(*,*) "p1 = ", pper
 1020 format(" omega, akx, omega/akx = ",1p4e12.4)
      go to 50
!
!.....psi, Bz, pe (itype=1)
 10   continue
      per1 = psiper
      per2 = bzper
      per3 = peper
      go to 50
!
!.....phi, vz, chi (itype=2)
 20   continue
      per1 = phiper
      per2 = vzper
      per3 = chiper
      go to 50
!
!.....n (itype=3)
 30   continue
!     
      per1 = nper
      per2 = nper
      per3 = nper
      go to 50
!
!.....p (itype=4)
 40   continue
      per1 = pper
      per2 = pper
      per3 = pper
      go to 50 
!     
 50   continue
!   
      end
!
! ===============================================================
      subroutine cubic_roots(coef, root, error)
!
!     calculate the roots of a 3rd degree polynomial
!     coef(4)*x**3 + coef(3)*x**2 + coef(2)*x + coef(1) = 0
!
      use basic
      implicit none
      real, intent(out), dimension(3) :: root, error
      real, intent(inout), dimension(4) :: coef
      real qqs, rrs, ths, fun, dfun
      integer r, ll
!
!     normalize the coefficients to the 3rd degree coefficient
      coef(1) = coef(1)/coef(4)
      coef(2) = coef(2)/coef(4)
      coef(3) = coef(3)/coef(4)
      coef(4) = 1.0
!     
!     solve using method in Numerical Recipes
      qqs = ((coef(3)**2 - 3.*coef(2))/9.)
      rrs = (2.*coef(3)**3 - 9.*coef(3)*coef(2) + 27.*coef(1))/54.
      ths = acos(rrs/sqrt(qqs**3))

      root(1) = -2.*sqrt(qqs)*cos(ths/3.) - coef(3)/3.
      root(2) = -2.*sqrt(qqs)*cos((ths+6.283185307)/3.) - coef(3)/3.
      root(3) = -2.*sqrt(qqs)*cos((ths-6.283185307)/3.) - coef(3)/3.
!
!........refine with Newton's method
      do r=1,3
         do ll=1,3
            fun  = root(r)**3 + coef(3)*root(r)**2 + coef(2)*root(r)    &
     &           + coef(1)
            dfun = 3.*root(r)**2 + 2.*coef(3)*root(r) + coef(2)
            root(r) = root(r) - fun/dfun
         enddo
         error(r) =  root(r)**3 + coef(3)*root(r)**2 + coef(2)*root(r)  &
     &        + coef(1)
      enddo
!      
      end
!======================================================
      subroutine newterm(termbbdum)
      use basic
      use t_data
      integer p,q
      dimension termbbdum(18,18)
!
      do i=1,18
      do j=1,18
      sum = 0.
      do p=1,20
      do q=1,20
      term = 0
      if(ni(p)+ni(q)-1 .ne. 0) go to 100
      if(mi(p)+mi(q)-1 .le. 0) go to 100
      nexp = mi(p)+mi(q)-1
      term = ni(p)*mi(q)*(mi(q)-1)*btri(1)**nexp/nexp
  100 continue
      if(ni(p)+ni(q)-3 .ne. 0) go to 200
      nexp = mi(p)+mi(q)+1
      term = term + ni(p)*ni(q)*(ni(q)-1)*btri(1)**nexp/nexp
  200 continue
      sum = sum + gtri(p,i,1)*gtri(q,j,1)*term
!
      enddo ! on q
      enddo ! on p
      termbbdum(i,j) = sum
      enddo ! on j
      enddo ! on i
!
      return
      end
!============================================================
      subroutine r4def
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
#ifdef mpi
      use supralu_dist_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer icount, irect, jrect, iodd, itri, ll, jfirst, jlast
      integer i, i1, i3
      integer j, jone
      integer ier

      r4 = 0.

      itri = 0
#ifdef mpi
      call MPE_DECOMP1D(maxrank,m-1+jper,myrank,jfirst,jlast)
      if(jfirst.le.0 .or. jlast.gt.m-1+jper) then
         write(*,*) myrank,maxrank,jfirst,jlast
         call safestop(1)
      endif
#else
      jfirst = 1
      jlast = m-1+jper
#endif
!
!
      do iodd=1,ioddmx
      do i=1,18
      do j=1,18
!
!
      do jrect=jfirst,jlast
      do irect=1,n-1+iper
         ll = irect + (jrect-1)*(n-1+iper)
         itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
!
         i1 = isvaln(itri,i)
         i3 = i1 + 12
!
         jone = isval1(itri,j)
!
!
!        gravity
!        =======
         if(idens.eq.1) then
            r4(i1) = r4(i1) + grav*dt*xterm(iodd,i,j)*den(jone)
            if(numvar.ge.3) then 
               r4(i3) = r4(i3) - grav*dt*yterm(iodd,i,j)*den(jone)
            endif
         endif
!
      enddo                     ! on irect
      enddo                     ! on jrect
      enddo                     ! on j
      enddo                     ! on i
      enddo                     ! on iodd
!

      if(myrank.eq.0 .and. iprint.eq.1) then
         write(*,*) "before r4 allreduce"
      endif

#ifdef mpi
      icount = m*n*numvar*6
      r4temp = r4
      call MPI_ALLREDUCE(r4temp,r4,icount,MPI_DOUBLE_PRECISION,         &
     &     MPI_SUM, MPI_COMM_WORLD,ier)
      if(ier /= 0) then
         print *,'Error with MPI_allreduce-r4',ier
         call safestop(1)
      endif
#endif
      return
      end

!============================================================
      subroutine q34def
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
#ifdef mpi
      use supralu_dist_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer icount, irect, jrect, iodd, itri, ll, jfirst, jlast
      integer i, i1, i3
      integer j, jone, j01
      integer k, k01
      integer l, lone
      integer ier, ir1, ir2, ir3
      real :: factor, chif, x, z

      q34 = 0.

      itri = 0
#ifdef mpi
      call MPE_DECOMP1D(maxrank,m-1+jper,myrank,jfirst,jlast)
      if(jfirst.le.0 .or. jlast.gt.m-1+jper) then
         write(*,*) myrank,maxrank,jfirst,jlast
         call safestop(1)
      endif
#else
      jfirst = 1
      jlast = m-1+jper
#endif
!
      ! prepare files for anisotropic heat flux
      if (chipar.ne.0) then
         rewind(174)
         rewind(175)
         rewind(176)
         rewind(177)
      endif
!
      do iodd=1,ioddmx
      do i=1,18
      do j=1,18
!
      if(iread(74).eq.1)                                                &
     &     read(174) ir1,ir2,ir3,((t2terml(k,l),l=1,18),k=1,18)
      if(iread(75).eq.1)                                                &
     &     read(175) ir1,ir2,ir3,((t3terml(k,l),l=1,18),k=1,18)
      if(iread(76).eq.1)                                                &
     &     read(176) ir1,ir2,ir3,((t4terml(k,l),l=1,18),k=1,18)
      if(iread(77).eq.1)                                                &
     &     read(177) ir1,ir2,ir3,((t5terml(k,l),l=1,18),k=1,18)
!
      do jrect=jfirst,jlast
      do irect=1,n-1+iper
         ll = irect + (jrect-1)*(n-1+iper)
         itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
!
         i1 = isvaln(itri,i)
         i3 = i1 + 12
!
         jone = isval1(itri,j)
         j01 = isval0(itri,j)

         x = (irect-iodd/3.)*deex
         z = (jrect-1+iodd/3.)*deez
         factor = 1.
         if(imask.eq.1) call mask(x,z,factor)

         chif = (gam-1.)*chipar*factor


         ! electron pressure source term
         q34(i3) = q34(i3) + dterm(iodd,i,j)*sphie(jone)
!
         ! anisotropic heat flux
         if(chipar.ne.0.) then
            do k=1,18
               k01 = isval0(itri,k)
               
               do l=1,18
                  lone = isval1(itri,l)
                  
                  q34(i3) = q34(i3) - 0.5*chif*                          &
     &                 ((t3terml(k,l) + t5terml(k,l))*tw50(lone)         &
     &                 -(t2terml(k,l) + t4terml(k,l))*tw40(lone))*       &
     &                 phi0(j01)*phi0(k01)
                  
               enddo
            enddo
         endif
 !
      enddo ! on irect
      enddo ! on jrect
      enddo ! on j
      enddo ! on i
      enddo ! on iodd
!

      if(myrank.eq.0 .and. iprint.eq.1) then
         write(*,*) "before q34 allreduce"
      endif

#ifdef mpi
      icount = m*n*numvar*6
      q34temp = q34
      call MPI_ALLREDUCE(q34temp,q34,icount,MPI_DOUBLE_PRECISION,         &
     &     MPI_SUM, MPI_COMM_WORLD,ier)
      if(ier /= 0) then
         print *,'Error with MPI_allreduce-q34',ier
         call safestop(1)
      endif
#endif
      return
      end

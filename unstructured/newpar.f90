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
  implicit none
#ifdef mpi
  include 'mpif.h'
#endif

#ifdef IS_LIBRARY
  integer, intent(in) :: isfirst, inmyrank, inmaxrank
#endif
  integer :: j, i, ier, numvars, ifail, maxts, numelms, numnodes
  integer :: index1, index2, ndofs
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
  call sludinit
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_size:',ier
     call safestop(1)
  endif
  ! initialize autopack
  call AP_INIT()
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
  if(myrank.eq.0) then
     call date_and_time( datec, timec)
     write(*,1001) datec(1:4),datec(5:6),datec(7:8),                        &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("M3D-C1 VERSION 1.0    DATE: "a4,1x,a2,1x,a2,3x,               &
          "TIME: "a2,":",a2,":",a4,/)
  endif

  pi = acos(-1.)

  call numnod(numnodes)
  call numfac(numelms)
  write(*,*) 'numnodes and numfaces',numnodes,numelms

  ! arrays for triangle parameters
  allocate(atri(numelms),btri(numelms),ctri(numelms),ttri(numelms),      &
       gtri(20,18,numelms)) 
  
  if(maxrank .eq. 1) call precalc_whattri()
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

  if(ipres.eq.1) then
     pefac = 1.
  else
     if(p0.gt.0) then
        pefac = (p0-pi0)/p0
     else
        pefac = 0.
     endif
     print *, "pefac = ", pefac
  endif


  ! calculate the RHS (forcing function)
  call rhsdef
  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'finished rhsdef'
  call numprocdofs(1, j)
  call numdofs(1, ndofs)
  write(*,*) 'proc, owned dofs, needed dofs',myrank, j, ndofs
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
        den0 = 1. ! acbauer - temp fix for periodic bc's on 11/14/06
        if(myrank.eq.0 .and. iprint.gt.0) write(*,*) 'finished denequ'
        if(maxrank.eq.1) call oneplot(den0,1,1,"den0",1)
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
        velold = 1. ! acbauer - temp fix for periodic bc's on 11/14/06
        call phiequ(phiold, numvar)
        phiold = 0. ! acbauer - temp fix for periodic bc's on 11/14/06
     endif

     vel0 = velold
     phi0 = phiold
     denold = den0

     ! calculate initial perturbed fields
     call velinit(vel)
     vel = 1. ! acbauer - temp fix for periodic bc's on 11/14/06
     call phiinit(phi)
     phi = 1. ! acbauer - temp fix for periodic bc's on 11/14/06
     if(idens.eq.1) then
        call deninit(den)
        den = 1. ! acbauer - temp fix for periodic bc's on 11/14/06
     endif

     phis = 0
     vels = 0
     dens = 0
  endif                     !  end of the branch on restart/no restart

  ! correct for left-handed coordinates
!!$  do i=1,numnodes
!!$     call entdofs(numvar, i, 0, ibegin, iendplusone)
!!$     if(numvar.eq.1) then
!!$        j = 5
!!$     else
!!$        j = 11
!!$     endif
!!$     phi(ibegin:ibegin+j) = -phi(ibegin:ibegin+j)
!!$     vel(ibegin:ibegin+j) = -vel(ibegin:ibegin+j)
!!$     phi0(ibegin:ibegin+j) = -phi0(ibegin:ibegin+j)
!!$     vel0(ibegin:ibegin+j) = -vel0(ibegin:ibegin+j)
!!$     phiold(ibegin:ibegin+j) = -phiold(ibegin:ibegin+j)
!!$     velold(ibegin:ibegin+j) = -velold(ibegin:ibegin+j)
!!$  enddo
! the above code looks like it just makes all of the dofs/unknowns negative
! which is what the code below does while also taking into account
! properly that multiple nodes may share the same dof/unknown
  call numdofs(numvar, ndofs)
  do i=1,ndofs
     phi(i) = -phi(i)
     vel(i) = -vel(i)
     phi0(i) = -phi0(i)
     vel0(i) = -vel0(i)
     phiold(i) = -phiold(i)
     velold(i) = -velold(i)
  enddo

  if(maxrank.eq.1) call plotit(vel,phi,0)

  ! calculate the equilibrium current density
  if(itaylor.ne.4) then
     call newvar(phi0,jphi0,numvar,1,VAR_J,1)
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

  if(maxrank .eq. 1) call plotit(vel+vel0,phi+phi0,1)

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
        if(maxrank .eq. 1) call output
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
  call freesmo(gsmatrix_sm)
  call freesmo(s6matrix_sm)
  call freesmo(s8matrix_sm)
  call freesmo(s7matrix_sm)
  call freesmo(s4matrix_sm)
  call freesmo(s3matrix_sm)
  call freesmo(s5matrix_sm)
  call freesmo(s1matrix_sm)
  call freesmo(s2matrix_sm)
  call freesmo(d1matrix_sm)
  call freesmo(d2matrix_sm)
  call freesmo(d4matrix_sm)
  call freesmo(d8matrix_sm)
  call freesmo(r1matrix_sm)
  call freesmo(r2matrix_sm)
  call freesmo(r8matrix_sm)
  call freesmo(q2matrix_sm)
  call freesmo(q8matrix_sm)  
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
  real(r8), allocatable:: temp(:)
  
  ! acbauer -- these variables are not initialized before use
  ediff = 1.
  etotd = 1.

  call numnod(numnodes)
  call numfac(numelms)
  call numdofs(1, ndofs)

  ! define the inverse density array deni
  if(idens.eq.1) then
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call inverse(den+den0,deni)
     else
        call inverse(den,deni)
     endif
  endif

  ! define the current vector jphi and the RHS vectors sb1 and sb2
  !      call newvar(phi,jphi,numvar,1,VAR_J,1)
  call newvar(phi,sb1,numvar,1,VAR_SB1,1)
  if(numvar.ge.2) call newvar(phi,sb2,numvar,1,VAR_SB2,1)
  if(numvar.ge.3) then
     call newvar(phi,sp1,numvar,1,VAR_SP1,1)
     if(ipres.eq.1) then
        call newvar(phi,sb3,numvar,1,VAR_SB3,1)
     else
        sb3 = sp1
     endif
  end if

  if(numvar.ge.2 .and. iconstflux.eq.1) call conserve_tflux()

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
  call matrixvectormult(r1matrix_sm, phi, b1vector)
  
  ! vtemp = d1matrix_sm * vel(n)
  vtemp = 0.
  if(iprint.eq.1)write(*,*) "before second sparseR8_A_dot_X"
  call matrixvectormult(d1matrix_sm,vel,vtemp)

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
     call newvar(vtemp,com,numvar,3,VAR_COM,0)
     
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
        call newvar(vtemp,com,numvar,3,VAR_COM,0)
     endif
     if(maxrank .eq. 1) call oneplot(com,1,1,"com",0)
     
  endif

  call newvar(vtemp,vor,numvar,1,VAR_VOR,1)
  
  if(hyperc.gt.0) then
     !
     !.......calculate vorticity, apply smoothing operator, and redefine vor array
     call smoother1(vor,vtemp,numnodes,numvar,1)
     call newvar(vtemp,vor,numvar,1,VAR_VOR,1)
     if(maxrank .eq. 1) call oneplot(com,1,1,"vor",0)
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
     call matrixvectormult(r8matrix_sm,vel,b2vector)

     ! b3vector = q8matrix_sm * vel(n)
     call matrixvectormult(q8matrix_sm,veln,b3vector)

     ! temp = d8matrix_sm * phi(n)
     call createvec(temp, 1)
     temp = 0.
     call matrixvectormult(d8matrix_sm,den,temp)
     call numdofs(numvar,ndofs)
     allocate(itemp(ndofs)) ! this is used to make sure that we don't double count the sum for periodic dofs
     itemp = 1

     do l=1,numnodes
        call entdofs(1, l, 0, ibegin, iendplusone)
        call entdofs(numvar, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           temp(ibegin+i) = temp(ibegin+i) + itemp(ibegin+i) * &
                (b2vector(ibeginnv+i) + b3vector(ibeginnv+i) + qn4(ibegin+i))
           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)

     do l=1,nbcn
        temp(iboundn(l)) = 0.
        if(linear.eq.0 .and. eqsubtract.eq.0) then
           temp(iboundn(l)) = temp(iboundn(l)) + denold(iboundn(l))
        endif
     enddo

     ! solve linear system...LU decomposition done first time
! -- okay to here      call printarray(temp, 150, 0, 'vtemp on')
     call solve(s8matrix_sm, temp, jer)
     if(jer.ne.0) then
        write(*,*) 'after 2nd sparseR8d_solve', jer
        call safestop(29)
     endif

     ! new field solution at time n+1 (or n* for second order advance)
     den = temp
     call deletevec(temp)
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
  call matrixvectormult(r2matrix_sm,vel,b2vector)

  ! b3vector = q2matrix_lu * vel(n)
  call matrixvectormult(q2matrix_sm,veln,b3vector)

  ! vtemp = d2matrix_sm * phi(n)
  vtemp = 0.
  call matrixvectormult(d2matrix_sm,phi,vtemp)

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
           tflux = tflux + d2term(j)*(phi(j2) + phi0(j2) - phiold(j2))
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

  return
end subroutine conserve_tflux

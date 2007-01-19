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
  use hdf5_output
  implicit none
#ifdef mpi
  include 'mpif.h'
#endif

#ifdef IS_LIBRARY
  integer, intent(in) :: isfirst, inmyrank, inmaxrank
#endif
  integer :: j, i, ier, numvars, ifail, maxts, numelms, numnodes
  integer :: index1, index2, ndofs, ibegin, iendplusone
  integer :: ibegin1, iendplusone1

  real :: dens, dtmin, ratemin, ratemax
  real :: tstart, tend

  double precision :: coords(3)
  integer :: nodeids(4)

  integer, allocatable ::  itemp(:)

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
1001 format("M3D-C1 VERSION 1.0    DATE: ", a4,1x,a2,1x,a2,3x,               &
          "TIME: ",a2,":",a2,":",a4,/)
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

  ! initialize hdf5
  call hdf5_initialize(ier)
  if(ier.lt.0) then 
     print *, "Error initializing HDF5"
  end if

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
  
  ! output simulation parameters
  call hdf5_write_parameters(ier)

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
     call initial_conditions()

     if(idens.eq.1 .and. maxrank.eq.1) call oneplot(den0,1,1,"den0",1)

     velold = vel0
     phiold = phi0
     denold = den0

     ! correct for left-handed coordinates
     call numdofs(numvar, ndofs)
     allocate(itemp(ndofs))
     itemp = 1
     if(numvar.eq.1) then
        j = 5
     else
        j = 11
     endif
     do i=1,numnodes
        call entdofs(numvar, i, 0, ibegin, iendplusone)
        if(itemp(ibegin) .eq. 1) then
           phi(ibegin:ibegin+j) = -phi(ibegin:ibegin+j)
           vel(ibegin:ibegin+j) = -vel(ibegin:ibegin+j)
           phi0(ibegin:ibegin+j) = -phi0(ibegin:ibegin+j)
           vel0(ibegin:ibegin+j) = -vel0(ibegin:ibegin+j)
           phiold(ibegin:ibegin+j) = -phiold(ibegin:ibegin+j)
           velold(ibegin:ibegin+j) = -velold(ibegin:ibegin+j)
           itemp(ibegin) = 0
        endif
     enddo
     deallocate(itemp)
     
     if(maxrank.eq.1) call plotit(vel,phi,0)

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
  endif                     !  end of the branch on restart/no restart

  
  ! create the newvar matrices
  call create_newvar_matrix(s6matrix_sm, 1)
  call create_newvar_matrix(s3matrix_sm, 0)

  ! calculate the equilibrium vorticity and velocity divergence
  call newvar_gs(vel, vor, numvar,1,1)
  call newvar_gs(phi+phi0, jphi, numvar,1,1)
  if(numvar.ge.3) then 
     call newvar_gs(vel+vel0, com, numvar,3,0)
  else
     com = 0.
  endif

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
  if (maxrank .eq. 1) call output

  call hdf5_write_scalars(ier)
  call hdf5_write_time_slice(ier)
  call hdf5_flush(ier)

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
     
     ! Write ictrans output
     if(maxrank .eq. 1) then
        if(itimer.eq.1) call second(tstart)
        call output
        if(itimer.eq.1) then
           call second(tend)
           write(*,*) "Time spent in output: ", tend - tstart
        endif
     endif

     ! Write HDF5 output
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call hdf5_write_scalars(ier)
     if(mod(ntime,ntimepr).eq.0) call hdf5_write_time_slice(ier)
     call hdf5_flush(ier)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        write(*,*) "Time spent in hdf5 output: ", tend - tstart
     end if

!     if(linear.eq.1) call scaleback

!!$     if(ekin .gt. 100.) then
!!$        write(*,*) 'ekin is greater than 100'
!!$        go to 100
!!$     endif
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
  ntime = max(ntime, ntimer)-1
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

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Calling plote..."
  if (myrank.eq.0) call plote
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Calling plote..."
  
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
  print *, myrank, "stopping"
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

  integer :: l, jer, i
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  integer, allocatable:: itemp(:)
  integer :: ndofs, numnodes

  real :: tstart, tend
  real, allocatable :: temp(:)
  
  call numnod(numnodes)

  ! define auxiliary variables
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  !   inverse density
  if(idens.eq.1) then
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call inverse(den+den0,deni)
     else
        call inverse(den,deni)
     endif
  endif
  !   toroidal current
  call newvar_gs(phi+phi0, jphi, numvar,1,1)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent defining auxiliary varibles:", tend - tstart
  endif


  ! define source terms
  ! ~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call define_sources
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent defining sources:", tend - tstart
  endif


  ! conserve toroidal flux
  ! ~~~~~~~~~~~~~~~~~~~~~~
  if(numvar.ge.2 .and. iconstflux.eq.1) call conserve_tflux()


  ! calculate matrices for time advance
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

  vtemp = vtemp + b1vector + r4

  ! apply boundary conditions
  do l=1,nbcv
!!$     if(iboundv2(l).gt.0) then
!!$        vtemp(iboundv2(l)) = vtemp(iboundv2(l)) - vtemp(iboundv(l))
!!$     endif
     vtemp(iboundv(l)) = velbounds(l)
  enddo

  ! solve linear system with rhs in vtemp (note LU-decomp done first time)
  if(myrank.eq.0 .and. iprint.eq.1) write(*,*) "before dsupralu_solve_s1handle"
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call solve(s1matrix_sm, vtemp, jer)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent in vel solve:", tend - tstart
  endif
  if(myrank.eq.0 .and. iprint.eq.1) write(*,*) "after dsupralu_solve_s1handle"
  if(jer.ne.0) then
     write(*,*) 'after sparseR8d_solve', jer
     call safestop(42)
  endif
! ok to here     call printarray(vtemp, 150, 0, 'vtemp aa')

!
!.....coding to calculate the error in the delsquared chi equation
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(numvar.ge.3) then
     call newvar_gs(vtemp,com,numvar,3,0)

!!$     call calc_chi_error(chierror)
!!$     if(myrank.eq.0 .and. iprint.ge.1) then
!!$        print *, "Error in com = ", chierror 
!!$     endif
     
     if(hyperc.gt.0) then
        call smoother3(com,vtemp,numnodes,numvar,3)
        call newvar_gs(vtemp,com,numvar,3,0)
     endif
  endif

  call newvar_gs(vtemp,vor,numvar,1,1)
  
  if(hyperc.gt.0) then

     ! calculate vorticity, apply smoothing operator, and redefine vor array
     call smoother1(vor,vtemp,numnodes,numvar,1)
!!$     if(maxrank .eq. 1) call oneplot(vtemp,1,numvar,"vte3",0)
     call newvar_gs(vtemp,vor,numvar,1,1)
!!$     if(maxrank .eq. 1) call oneplot(vor,1,1,"vor",0)
  endif
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     print *, " onestep: Time spent in smoothers", tend - tstart
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
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(s8matrix_sm, temp, jer)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        write(*,*) " onestep: Time spent in den solve:", tend - tstart
     endif
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

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  ! b2vector = r2matrix_lu * vel(n+1)
  call matrixvectormult(r2matrix_sm,vel,b2vector)

  ! b3vector = q2matrix_lu * vel(n)
  call matrixvectormult(q2matrix_sm,veln,b3vector)

  ! vtemp = d2matrix_sm * phi(n)
  vtemp = 0.
  call matrixvectormult(d2matrix_sm,phi,vtemp)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent in field matrixvectormult:", tend - tstart
  endif

  vtemp = vtemp + b2vector + b3vector + q4
  
  ! Insert boundary conditions
  do l=1,nbcp
     vtemp(iboundp(l)) = psibounds(l)
     if(linear.eq.0 .and. eqsubtract.eq.0) then
        vtemp(iboundp(l)) = vtemp(iboundp(l)) + phiold(iboundp(l))
     endif
  enddo

  ! solve linear system...LU decomposition done first time
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call solve(s2matrix_sm, vtemp, jer)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent in field solve:", tend - tstart
  endif
  if(jer.ne.0) then
     write(*,*) 'after 2nd sparseR8d_solve', jer
     call safestop(29)
  endif
#endif

! new field solution at time n+1 (or n* for second order advance)
!      phiold = phi
  phi = vtemp


  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call energy
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent in energy:", tend - tstart
  endif

!!$  eerror = 0.
!!$  if(ntime.gt.5) then
!!$     eerror = 2.*abs(ediff-etotd)/(abs(ediff)+abs(etotd))
!!$  endif
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
!!$  if(facd .ne. 0 .and. facw.eq.1 .and. time.ge.0.1) then
!!$     facd = 0.
!!$     ntimemin = max(ntimemin,ntime+1)
!!$  endif

end subroutine onestep


! ======================================================================
! conserve_tflux
! --------------
!
! adjusts the boundary conditions to conserve toroidal flux
!
! ======================================================================
subroutine conserve_tflux

  use basic
  use t_data
  use arrays

  implicit none
  include "mpif.h"
  
  integer :: numelms, index, itri, j, jone, j2, ivertex, ier, ndofs
  real :: correction
  real :: d2term(18), fintl(-6:maxi,-6:maxi)
  double precision :: valsin(3), valsout(3)

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

  if(maxrank .gt. 1) then
     valsin(1) = totcur
     valsin(2) = area
     valsin(3) = tflux
     call MPI_ALLREDUCE(valsin, valsout, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
     totcur = valsout(1)
     area = valsout(2)
     tflux = valsout(3)
  endif


  ! adjust toroidal field and boundary condition to conserve toroidal flux
  if(numvar.ge.2) then
     correction = tflux/area
     gbound = gbound - correction
     call numdofs(numvar, ndofs)
     do j=7,ndofs,6*numvar
        phi(j) = phi(j) - correction
     enddo
  endif

  return
end subroutine conserve_tflux


! ======================================================================
! calc_chi_error
! --------------
!
! calculates the error in the equation Laplacian(chi) = com
! ======================================================================
subroutine calc_chi_error(chierror)

  use arrays
  use t_data

  implicit none

  include "mpif.h"

  real, intent(out) :: chierror

  integer :: ier, itri, j, jone, numelms
  real :: chierror_local
  real :: fintl(-6:maxi,-6:maxi), d2term(18)

  call numfac(numelms)

  chierror_local = 0

  do itri=1,numelms
     call calcfint(fintl, maxi, atri(itri),btri(itri), ctri(itri))
     call calcd2term(itri, d2term, fintl)
     do j=1,18
        jone = isval1(itri,j)
        chierror_local = chierror_local + d2term(j)*com(jone)
     enddo
  enddo                  ! loop over itri
  call mpi_allreduce(chierror_local, chierror, 1, MPI_DOUBLE, &
       MPI_SUM, MPI_COMM_WORLD, ier)

end subroutine calc_chi_error


! ======================================================================
! energy
! ------
!
! calculates the various components of the energy
! ====================================================================== 

subroutine energy

  use p_data
  use t_data
  use basic
  use arrays
  use nintegrate_mod

  implicit none
#ifdef mpi
  include "mpif.h"
#endif
  integer :: itri, i, def_fields 
  real :: terma, termb, termd, hypf, hypi, hypv, hypc
  real :: gamfac, deex
  double precision temp(18), temp2(18)
  integer :: numelms
  real, dimension(20) :: avec

  call numfac(numelms)

  ekino = ekin
  emago = emag
  ekindo = ekind
  emagdo = emagd

  ekinpo = ekinp
  emagpo = emagp
  ekinpdo = ekinpd
  emagpdo = emagpd

  ekinto = ekint
  emagto = emagt
  ekintdo = ekintd
  emagtdo = emagtd

  ekinpho = ekinph
  ekintho = ekinth
  emagpho = emagph
  emagtho = emagth

  ekin3o = ekin3 
  ekin3do = ekin3d
  ekin3ho = ekin3h 
  emag3o = emag3
  emag3do = emag3d
  emag3ho = emag3h
  
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


  ! Determine which fields need to be calculated
  def_fields = FIELD_PSI + FIELD_PHI + FIELD_J + FIELD_VOR
  if(numvar.ge.2) def_fields = def_fields + FIELD_V + FIELD_I
  if(numvar.ge.3) def_fields = def_fields + FIELD_CHI + FIELD_PE + FIELD_COM
  if(idens.eq.1) def_fields = def_fields + FIELD_N
  if(ipres.eq.1) def_fields = def_fields + FIELD_P

  ! volume terms
  do itri=1,numelms
     ! calculate the field values and derivatives at the sampling points
     call define_fields_79(itri, def_fields)

     call getdeex(itri,deex)
     hypf = hyper *deex**2
     hypi = hyperi*deex**2
     hypv = hyperv*deex**2
     hypc = hyperc*deex**2

     if(idens.eq.0) then
        ekinp = ekinp + .5* &
             (int3(r2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),weight_79,79) &
             +int3(r2_79,pht79(:,OP_DR),pht79(:,OP_DR),weight_79,79))
     else
        ekinp = ekinp + .5* &
             (int4(r2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             +int4(r2_79,pht79(:,OP_DR),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79))
     endif

     ekinpd = ekinpd - amu*int2(pht79(:,OP_LP),pht79(:,OP_LP),weight_79,79)

     ekinph = ekinph - hypc*amu* &
          (int2(vot79(:,OP_DZ),vot79(:,OP_DZ),weight_79,79) &
          +int2(vot79(:,OP_DR),vot79(:,OP_DR),weight_79,79))

     emagp = emagp + .5* &
             (int3(ri2_79,pst79(:,OP_DZ),pst79(:,OP_DZ),weight_79,79) &
             +int3(ri2_79, pst79(:,OP_DR),pst79(:,OP_DR),weight_79,79))

     emagpd = emagpd - etar*int3(ri2_79,pst79(:,OP_GS),pst79(:,OP_GS),weight_79,79)

     emagph = emagph - hypf*etar* &
          (int2(jt79(:,OP_DZ),jt79(:,OP_DZ),weight_79,79) &
          +int2(jt79(:,OP_DR),jt79(:,OP_DR),weight_79,79))

     if(numvar.ge.2) then
        if(idens.eq.0) then
           ekint = ekint + .5*int2(vzt79(:,OP_1),vzt79(:,OP_1),weight_79,79)
        else
           ekint = ekint + .5*int3(vzt79(:,OP_1),vzt79(:,OP_1),nt79(:,OP_1),weight_79,79)
        endif

        ekintd = ekintd - amu* &
             (int2(vzt79(:,OP_DZ),vzt79(:,OP_DZ),weight_79,79) &
             +int2(vzt79(:,OP_DR),vzt79(:,OP_DR),weight_79,79))

        ekinth = ekinth - amu*hypv*int2(vzt79(:,OP_LP),vzt79(:,OP_LP),weight_79,79)

        emagt = emagt + .5*int3(ri2_79,bzt79(:,OP_1),bzt79(:,OP_1),weight_79,79)

        emagtd = emagtd - etar* &
             (int2(bzt79(:,OP_DZ),bzt79(:,OP_DZ),weight_79,79) &
             +int2(bzt79(:,OP_DR),bzt79(:,OP_DR),weight_79,79))

        emagth = emagth - etar*hypi*int2(bzt79(:,OP_LP),bzt79(:,OP_LP),weight_79,79)
     endif

     if(numvar.ge.3) then
        if(idens.eq.0) then
           ekin3 = ekin3 + .5* &
                (int2(cht79(:,OP_DZ),cht79(:,OP_DZ),weight_79,79) &
                +int2(cht79(:,OP_DR),cht79(:,OP_DR),weight_79,79))
           if(itor.eq.1) then
              ekin3 = ekin3 + &
                   (int3(r_79,cht79(:,OP_DZ),pht79(:,OP_DR),weight_79,79) &
                   -int3(r_79,cht79(:,OP_DR),pht79(:,OP_DZ),weight_79,79))
           endif
        else
           ekin3 = ekin3 + .5* &
                (int3(cht79(:,OP_DZ),cht79(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
                +int3(cht79(:,OP_DR),cht79(:,OP_DR),nt79(:,OP_1),weight_79,79))
           if(itor.eq.1) then
              ekin3 = ekin3 + &
                   (int4(r_79,cht79(:,OP_DZ),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
                   -int4(r_79,cht79(:,OP_DR),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79))
           endif
        endif

        ekin3d = ekin3d - amuc*int2(pht79(:,OP_LP),pht79(:,OP_LP),weight_79,79)
        
        ekin3h = ekin3h - hypc*amuc* &
             (int2(cot79(:,OP_DZ),cot79(:,OP_DZ),weight_79,79) &
             +int2(cot79(:,OP_DR),cot79(:,OP_DR),weight_79,79))
        
        emag3 = emag3 + int1(pt79,weight_79,79) / (gam - 1.)
        
     endif

  enddo     ! end of loop on itri

! parallel updating of values summed over elements
#ifdef mpi
  if(maxrank .gt. 1) then
     temp(1) = ekinp
     temp(2) = emagp
     temp(3) = ekinpd
     temp(4) = emagpd      
     temp(5) = ekint
     temp(6) = emagt
     temp(7) = ekintd
     temp(8) = emagtd
     temp(9) = ekinph
     temp(10) = ekinth
     temp(11) = emagph
     temp(12) = emagth
     temp(13) = ekin3
     temp(14) = ekin3d
     temp(15) = ekin3h
     temp(16) = emag3
     temp(17) = emag3d
     temp(18) = emag3h
         
     call mpi_allreduce(temp, temp2, 18, MPI_DOUBLE_PRECISION,  &
          MPI_SUM, MPI_COMM_WORLD, i) !checked that this should be MPI_DOUBLE_PRECISION
         
     ekinp = temp2(1)
     emagp = temp2(2)
     ekinpd = temp2(3)
     emagpd = temp2(4)      
     ekint = temp2(5)
     emagt = temp2(6)
     ekintd = temp2(7)
     emagtd = temp2(8)
     ekinph = temp2(9)
     ekinth = temp2(10)
     emagph = temp2(11)
     emagth = temp2(12)
     ekin3 = temp2(13)
     ekin3d = temp2(14)
     ekin3h = temp2(15)
     emag3 = temp2(16)
     emag3d = temp2(17)
     emag3h = temp2(18)
  endif !if maxrank .gt. 1
#endif

  ekin = ekinp + ekint + ekin3
  emag = emagp + emagt + emag3
  ekind = ekinpd + ekintd + ekin3d
  emagd = emagpd + emagtd + emag3d

  if(myrank.eq.0) then
     print *, "Energy at ntime = ", ntime
     print *, "ekinp, ekint, ekin3 = ", ekinp, ekint, ekin3
     print *, "ekinpd, ekintd, ekin3d = ", ekinpd, ekintd, ekin3d
     print *, "ekinph, ekinth, ekin3h = ", ekinph, ekinth, ekin3h
     print *, "emagp, emagt, emag3 = ", emagp, emagt, emag3
     print *, "emagpd, emagd, emag3d = ", emagpd, emagtd, emag3d
     print *, "emagph, emagth, emag3h = ", emagph, emagth, emag3h
  endif

  return
      
end subroutine energy


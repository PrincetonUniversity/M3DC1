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
  use diagnostics

  implicit none
#ifdef mpi
  include 'mpif.h'
#endif

#ifdef IS_LIBRARY
  integer, intent(in) :: isfirst, inmyrank, inmaxrank
#endif
  integer :: j, i, ier, ifail, maxts, numelms, numnodes
  integer :: ndofs, ibegin, iendplusone

  real :: dtmin, ratemin, ratemax
  real :: tstart, tend

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
  call init

  if(ipres.eq.1) then
     pefac = 1.
  else
     if(p0.gt.0) then
        pefac = (p0-pi0)/p0
     else
        pefac = 0.
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "pefac = ", pefac
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

  ! initialize hdf5
  call hdf5_initialize(ntimer, ier)
  if(ier.lt.0) then 
     print *, "Error initializing HDF5"
     call safestop(5)
  end if
  
  ! output simulation parameters
  if(irestart.eq.0) call hdf5_write_parameters(ier)
  
  ! create the newvar matrices
  call create_newvar_matrix(s6matrix_sm, 1)
  call create_newvar_matrix(s3matrix_sm, 0)

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

  if(itimer.eq.1) call reset_timings

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
  !   resistivity
  call newvar_eta
  !   toroidal current
  call newvar_gs(phi+phi0,jphi,1,1)
  !   vorticity
  call newvar_gs(vel+vel0, vor,1,1)
  !   compression
  if(numvar.ge.3) then 
     call newvar_gs(vel+vel0,com,3,0)
  else
     com = 0.
  endif
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif


  ! calculate other quantities of interest
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call total_flux
  tflux0 = tflux
  totcur0 = totcur
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent other quantities:", tend - tstart
  endif


  ! define source terms and calculate energy
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call define_sources
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_sources = t_sources + tend - tstart
  endif


  ! output initial conditions
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(irestart.eq.0) then
     if (maxrank .eq. 1) call output
     call hdf5_write_scalars(ier)
     call hdf5_write_time_slice(ier)
     if(itimer.eq.1) then 
        call hdf5_write_timings(ier)
        call reset_timings
     endif
     call hdf5_flush(ier)

  endif

  ! if there are no timesteps to calculate, then skip time loop
  if(ntimemax.le.ntimer) go to 101

  ! main time loop
  ! ~~~~~~~~~~~~~~
  do ntime=ntimer+1,ntimemax

     ! check for error
     if(ekin.ne.ekin .or. emag.ne.emag) then
        print *, "Error: energy is NaN"
        call safestop(3)
     endif

     !advance time
     time = time + dt

     if(myrank.eq.1 .and. iprint.ge.1) print *, "Before onestep"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call onestep
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_onestep = t_onestep + tend - tstart
     endif
     if(myrank.eq.1 .and. iprint.ge.1) print *, "After onestep"
!     call exportfield2(1,numvar,phi, ntime)


!     if(linear.eq.1) call scaleback

!!$     if(ekin .gt. 100.) then
!!$        write(*,*) 'ekin is greater than 100'
!!$        go to 100
!!$     endif

     
     ! Write ictrans output
     if(maxrank .eq. 1) then
        if(myrank.eq.1 .and. iprint.ge.1) print *, "Before output"
        if(itimer.eq.1) call second(tstart)
        call output
        if(itimer.eq.1) then
           call second(tend)
           t_output_cgm = t_output_cgm + tend - tstart
        endif
        if(myrank.eq.1 .and. iprint.ge.1) print *, "After onestep"
     endif

     ! Write HDF5 output
     if(myrank.eq.1 .and. iprint.ge.1) print *, "Before hdf5 output"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call hdf5_write_scalars(ier)
     if(mod(ntime,ntimepr).eq.0) then
        call hdf5_write_time_slice(ier)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_output_hdf5 = t_output_hdf5 + tend - tstart
     end if

     if(itimer.eq.1) then
        if(myrank.eq.0) call second(tstart)
        call hdf5_write_timings(ier)
        call reset_timings
     end if

     ! flush hdf5 data to disk
     call hdf5_flush(ier)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_output_hdf5 = t_output_hdf5 + tend - tstart
     end if
     if(myrank.eq.1 .and. iprint.ge.1) print *, "After hdf5 output"

  enddo ! ntime

 100  continue
  
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
!!$  if(maxrank.gt.1 .and. myrank.eq.0) then
!!$     write(*,5002) tsolve-tfirst,                                      &
!!$          numvar,amu,etar,dt,thimp,db,hyper,hyperi,hyperv,hyperc,      &
!!$          ratemin,ratemax,ajmax,graphit(ntimemax,25)
!!$     write(*,5003) linear, itaylor, 0, imask, irestart,                &
!!$          facd,bzero,eps
!!$
!!$     write(FILE__C1NEW,5002) tsolve-tfirst,                            &
!!$          numvar,amu,etar,dt,thimp,db,hyper,hyperi,hyperv,hyperc,      &
!!$          ratemin,ratemax,ajmax,graphit(ntimemax,25)
!!$     write(FILE__C1NEW,5003) linear, itaylor, 0, imask, irestart,      &
!!$          facd,bzero,eps
!!$  endif

  if(ntime.gt.1 .and. myrank.eq.0 .and. maxrank.eq.1) then
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

  if (myrank.eq.0 .and. maxrank.eq.1) call plote
  
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
  use diagnostics

  implicit none

  integer :: l, jer, i
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  integer, allocatable:: itemp(:)
  integer :: ndofs, numnodes

  real :: tstart, tend
  real, allocatable :: temp(:)
  
  call numnod(numnodes)

  ! apply loop voltage
  ! ~~~~~~~~~~~~~~~~~~
  fbound = fbound + dt*vloop/(2.*pi)


  ! calculate matrices for time advance
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ntime.le.ntimer+1.or. (linear.eq.0 .and. mod(ntime,nskip).eq.0)) then
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call ludefall
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif
  endif

  veln = vel
  

  ! Advance Velocity
  ! ================

  ! Calculate LU decomposition of velocity matrix if needed

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! b1vector = r1matrix_sm * phi(n)
  if(iprint.ge.1) write(*,*) "before sparseR8_A_dot_X"
  call matrixvectormult(r1matrix_sm, phi, b1vector)
  
  ! vtemp = d1matrix_sm * vel(n)
  vtemp = 0.
  if(iprint.eq.1)write(*,*) "before second sparseR8_A_dot_X"
  call matrixvectormult(d1matrix_sm,vel,vtemp)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_mvm = t_mvm + tend - tstart
  endif

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
     t_solve_v = t_solve_v + tend - tstart
  endif
  if(myrank.eq.0 .and. iprint.eq.1) write(*,*) "after dsupralu_solve_s1handle"
  if(jer.ne.0) then
     write(*,*) 'after sparseR8d_solve', jer
     call safestop(42)
  endif


  ! apply smoothing operators
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(hyperc.gt.0) then
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! smooth vorticity
     call newvar_gs(vtemp,vor,1,1)
     call smoother1(vor,vtemp,numnodes,numvar,1)

     ! smooth compression
     if(numvar.ge.3) then
        call newvar_gs(vtemp,com,3,0)

!!$        !
!!$        !.....coding to calculate the error in the delsquared chi equation
!!$        call calc_chi_error(chierror)
!!$        if(myrank.eq.0) then
!!$           print *, "Error in com = ", chierror 
!!$        endif

        call smoother3(com,vtemp,numnodes,numvar,3)     
     endif

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_smoother = t_smoother + tend - tstart
     endif
  endif

!!$  ! define vorticity and compression
!!$  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
!!$  call newvar_gs(vtemp,vor,1,1)
!!$  if(numvar.ge.3) call newvar_gs(vtemp,com,3,0)
!!$  if(myrank.eq.0 .and. itimer.eq.1) then
!!$     call second(tend)
!!$     t_aux = t_aux + tend - tstart
!!$  endif

!.....new velocity solution at time n+1 (or n* for second order advance)
  vel = vtemp

!
! Advance Density
! ===============
  if(idens.eq.1) then
     if(iprint.ge.1) write(*,*) "s8handle"

     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! b2vector = r8matrix_lu * vel(n+1)
     call matrixvectormult(r8matrix_sm,vel,b2vector)

     ! b3vector = q8matrix_sm * vel(n)
     call matrixvectormult(q8matrix_sm,veln,b3vector)

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif


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
        t_solve_n = t_solve_n + tend - tstart
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
     t_mvm = t_mvm + tend - tstart
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
     t_solve_b = t_solve_b + tend - tstart
  endif
  if(jer.ne.0) then
     write(*,*) 'after 2nd sparseR8d_solve', jer
     call safestop(29)
  endif
#endif

! new field solution at time n+1 (or n* for second order advance)
  phi = vtemp


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
  !   resistivity
  call newvar_eta
  !   toroidal current
  call newvar_gs(phi+phi0, jphi,1,1)
  !   vorticity
  call newvar_gs(vtemp,vor,1,1)
  !   compression
  if(numvar.ge.3) call newvar_gs(vtemp,com,3,0)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif


  ! calculate other quantities of interest
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call total_flux
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     call total_flux
  endif
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) " onestep: Time spent other quantities:", tend - tstart
  endif


  ! define source terms
  ! ~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call define_sources
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_sources = t_sources + tend - tstart
  endif 

end subroutine onestep


! ======================================================================
! conserve_flux
! --------------
!
! adjusts the boundary conditions to conserve flux
!
! ======================================================================
subroutine conserve_flux

  use basic
  use arrays
  use diagnostics

  implicit none
  
  integer :: ndofs, j

  ! adjust toroidal field and boundary condition to conserve toroidal flux
  if(numvar.ge.2 .and. iconstflux.eq.1) then
     gbound = (tflux0-tflux)/area
     call numdofs(numvar, ndofs)
     do j=7,ndofs,6*numvar
        phi(j) = phi(j) + gbound
     enddo
     if(myrank.eq.0) then
        print *, "Correction to toroidal flux: ", gbound*area
     end if
  endif

  return
end subroutine conserve_flux


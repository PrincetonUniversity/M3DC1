Program Reducedquintic

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
#ifdef _AIX
  include 'mpif.h'
#endif
#include "finclude/petsc.h"

  integer :: j, i, ier, ifail, maxts, numelms, numnodes
  integer :: ndofs, ibegin, iendplusone

  real :: dtmin, ratemin, ratemax
  real :: tstart, tend
  real :: tolerance
  real :: factor, hmin, hmax  

  integer, allocatable ::  itemp(:)
  PetscTruth :: flg
  PetscInt :: mpetscint,npetscint

  ! Start up message passing, SuperLU process grid
  call MPI_Init(ier)
  if (ier /= 0) then
     print *,'Error initializing MPI:',ier
     call safestop(1)
  endif
  call PetscInitialize(PETSC_NULL_CHARACTER, ier)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-m',mpetscint,flg,ier)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',npetscint,flg,ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_rank:',ier
     call safestop(1)
  endif
  ! initialize the SUPERLU process grid
  call initsolvers
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_size:',ier
     call safestop(1)
  endif
  ! initialize autopack
  call AP_INIT()

  print *, 'starting M3D-C1'
  if(myrank.eq.0 .and. itimer.ge.1) call second(tstart)
  call loadmesh("struct.dmg", "struct-dmg.sms")
  if(myrank.eq.0 .and. itimer.ge.1) then
     call second(tend)
     print *, 'Time spent in loadmesh: ', tend-tstart
  endif
  if(myrank.eq.0) then
     call date_and_time( datec, timec)
     write(*,1001) datec(1:4),datec(5:6),datec(7:8),                        &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("M3D-C1 VERSION 1.0    DATE: ", a4,1x,a2,1x,a2,3x,               &
          "TIME: ",a2,":",a2,":",a4,/)
  endif

  if(myrank.eq.0) print*, &
#ifdef NEW_VELOCITY
       "V = grad(U)xgrad(phi) + V grad(phi) + grad(chi)"
#else
       "V = r^2 grad(U)xgrad(phi) + r V grad(phi) + grad(chi)"
#endif

  pi = acos(-1.)

  call numnod(numnodes)
  call numfac(numelms)
  write(*,*) 'numnodes and numfaces',numnodes,numelms

  ! arrays for triangle parameters
  allocate(atri(numelms),btri(numelms),ctri(numelms),ttri(numelms),      &
       gtri(20,18,numelms)) 
  
!!$  if(maxrank .eq. 1) call precalc_whattri()
  ! special switch installed 12/03/04 to write slu matrices
  idebug = 0

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

  ! check time-integration options
  select case(integrator)
  case(1)
     ! For BDF2 integration, first timestep is Crank-Nicholson with thimp=1,
     ! and subsequent timesteps are BDF2.
     if(myrank.eq.0) print *, "Time integration: BDF2."
     thimp = 1.
     thimp_ohm = 1.
  case default
     if(myrank.eq.0) print *, "Time integration: Crank-Nicholson."
  end select

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
     time = timer
     ntime = ntimer
     if(istart.ne.0) then
        ntimer = 0
        timer = 0.
     endif
  else
     ptot = 0.
     ntimer = 0
     timer = 0.
     call initial_conditions()

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        phis = phi
        if(isplitstep.eq.1) then
           vels = vel
           if(idens.eq.1) dens = den
           if(ipres.eq.1) press = pres
        endif
     else
        phis = phi0
        if(isplitstep.eq.1) then
           vels = vel0
           if(idens.eq.1) dens = den0
           if(ipres.eq.1) press = pres0
        endif
     endif

     ! correct for left-handed coordinates
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "adjusting fields for left-handed coordinates"
     call numdofs(vecsize, ndofs)
     allocate(itemp(ndofs))
     itemp = 1
     do i=1,numnodes
        call entdofs(vecsize, i, 0, ibegin, iendplusone)
        if(itemp(ibegin) .eq. 1) then
           psi0_v(ibegin+psi_off:ibegin+psi_off+5) = &
                -psi0_v(ibegin+psi_off:ibegin+psi_off+5)
           psi1_v(ibegin+psi_off:ibegin+psi_off+5) = &
                -psi1_v(ibegin+psi_off:ibegin+psi_off+5)
           psis_v(ibegin+psi_off:ibegin+psi_off+5) = &
                -psis_v(ibegin+psi_off:ibegin+psi_off+5)

           phi0_v(ibegin+phi_off:ibegin+phi_off+5) = &
                -phi0_v(ibegin+phi_off:ibegin+phi_off+5)
           phi1_v(ibegin+phi_off:ibegin+phi_off+5) = &
                -phi1_v(ibegin+phi_off:ibegin+phi_off+5)
           phis_v(ibegin+phi_off:ibegin+phi_off+5) = &
                -phis_v(ibegin+phi_off:ibegin+phi_off+5)

           if(numvar.ge.2) then
              bz0_v(ibegin+bz_off:ibegin+bz_off+5) = &
                   -bz0_v(ibegin+bz_off:ibegin+bz_off+5)
              bz1_v(ibegin+bz_off:ibegin+bz_off+5) = &
                   -bz1_v(ibegin+bz_off:ibegin+bz_off+5)
              bzs_v(ibegin+bz_off:ibegin+bz_off+5) = &
                   -bzs_v(ibegin+bz_off:ibegin+bz_off+5)

              vz0_v(ibegin+vz_off:ibegin+vz_off+5) = &
                   -vz0_v(ibegin+vz_off:ibegin+vz_off+5)
              vz1_v(ibegin+vz_off:ibegin+vz_off+5) = &
                   -vz1_v(ibegin+vz_off:ibegin+vz_off+5)
              vzs_v(ibegin+vz_off:ibegin+vz_off+5) = &
                   -vzs_v(ibegin+vz_off:ibegin+vz_off+5)
           endif

           itemp(ibegin) = 0
        endif
     enddo
     deallocate(itemp)
     
     ! combine the equilibrium and perturbed fields of linear=0
     ! unless eqsubtract = 1
     if(linear.eq.0 .and. eqsubtract.eq.0) then
        phi = phi + phi0
        phi0 = 0.
        
        if(isplitstep.eq.1) then
           vel = vel + vel0
           vel0 = 0.
        
           if(idens.eq.1) then
              den = den + den0
              den0 = 0.
           endif

           if(ipres.eq.1) then
              pres = pres + pres0
              pres0 = 0.
           endif
        endif
     endif

     ! initialize t(n-1) values
     phiold = phi

     if(isplitstep.eq.1) then
        velold = vel
        if(idens.eq.1) denold = den
        if(ipres.eq.1) presold = pres
     endif
     
  endif                     !  end of the branch on restart/no restart

  ! initialize hdf5
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Initializing HDF5."
  call hdf5_initialize(ier)
  if(ier.lt.0) then 
     print *, "Error initializing HDF5"
     call safestop(5)
  end if
  
  ! output simulation parameters
  if(irestart.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "Writing simulation parameters."
     call hdf5_write_parameters(ier)

     ! Output the equilibrium
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call hdf5_write_time_slice(1,ier)
     endif
  end if
  
  ! create the newvar matrices
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Generating newvar matrices..."
  call create_newvar_matrix(s6matrix_sm, NV_DCBOUND)
  call create_newvar_matrix(s3matrix_sm, NV_NOBOUND)
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done generating newvar matrices."


  ifail=0
  ier = 0

  errori = 0.
  enormi = 0.
  ratioi = 0.
  dtmin = 0.001*dt

  if(itimer.eq.1) call reset_timings

  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call derived_quantities

  if(irestart.eq.0) then
     tflux0 = tflux
     totcur0 = totcur
  endif

  ! output initial conditions
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(irestart.eq.0) then
!!$     if (maxrank .eq. 1) call output
     call hdf5_write_scalars(ier)
     call hdf5_write_time_slice(0,ier)
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

     if(myrank.eq.0 .and. iprint.ge.1) print *, "ntime = ", ntime

     if(myrank.eq.0 .and. iprint.ge.1) print *, "Before onestep"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call onestep
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_onestep = t_onestep + tend - tstart
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "After onestep"
!     call exportfield2(1,numvar,phi, ntime)


!     if(linear.eq.1) call scaleback

!!$     if(ekin .gt. 100.) then
!!$        write(*,*) 'ekin is greater than 100'
!!$        go to 100
!!$     endif

     
!!$     ! Write ictrans output
!!$     if(maxrank .eq. 1) then
!!$        if(myrank.eq.0 .and. iprint.ge.1) print *, "Before output"
!!$        if(itimer.eq.1) call second(tstart)
!!$        call output
!!$        if(itimer.eq.1) then
!!$           call second(tend)
!!$           t_output_cgm = t_output_cgm + tend - tstart
!!$        endif
!!$        if(myrank.eq.0 .and. iprint.ge.1) print *, "After onestep"
!!$     endif

     ! Write output
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Writing output"
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Scalars (HDF5)."
     call hdf5_write_scalars(ier)
     if(mod(ntime,ntimepr).eq.0) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, "-Time slice (HDF5)."
        call hdf5_write_time_slice(0,ier)
        if(myrank.eq.0 .and. iprint.ge.1) print *, "-Restart file(s)."
        call wrrestart(time, max(maxts, ntimer))
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_output_hdf5 = t_output_hdf5 + tend - tstart
     end if

     if(itimer.eq.1) then
        if(myrank.eq.0) call second(tstart)
        if(myrank.eq.0 .and. iprint.ge.1) print *, "-Timings (HDF5)."
        call hdf5_write_timings(ier)
        call reset_timings
     end if

     ! flush hdf5 data to disk
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Flushing data to HDF5 file."
     call hdf5_flush(ier)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_output_hdf5 = t_output_hdf5 + tend - tstart
     end if
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Done output."

     ! feedback control on toroidal current
     if(itor.eq.1 .and. itaylor.eq.1) call control_pid
  enddo ! ntime

 100  continue
  
  maxts = ntime-1

101 continue

! below is for mesh adaptation
  tolerance = .00005 
  call outputfield(phi, numvar, 0, ntime, 123) 
  if(maxrank .eq. 1) then
!     call outputfield(phi, numvar, 0, ntime, 123) 
!     call writefieldatnodes(resistivity, 1, 1) 
     factor = .3
     hmin = .1
     hmax = .4
!     call hessianadapt(resistivity, 1, factor, hmin, hmax, ntime) 
  endif
  ratemin = 0.
  ratemax = 0.

!  call errorcalc(numvar, phi, 1)

  call wrrestart(time, max(maxts, ntimer))

!     free memory from sparse matrices
  call deletematrix(gsmatrix_sm)
  call deletematrix(s6matrix_sm)
  call deletematrix(s7matrix_sm)
  call deletematrix(s4matrix_sm)
  call deletematrix(s3matrix_sm)
  call deletematrix(s5matrix_sm)
  call deletematrix(s1matrix_sm)
  call deletematrix(s2matrix_sm)
  call deletematrix(d1matrix_sm)
  call deletematrix(d2matrix_sm)
  call deletematrix(d4matrix_sm)
  call deletematrix(q1matrix_sm)
  call deletematrix(r2matrix_sm)
  call deletematrix(q2matrix_sm)
  call deletematrix(r14matrix_sm)
  if(idens.eq.1) then
     call deletematrix(s8matrix_sm)
     call deletematrix(d8matrix_sm)
     call deletematrix(q8matrix_sm) 
     call deletematrix(r8matrix_sm)
  end if
  if(ipres.eq.1) then 
     call deletematrix(s9matrix_sm)
     call deletematrix(d9matrix_sm)
     call deletematrix(r9matrix_sm)
     call deletematrix(q9matrix_sm)
  endif
  call deletesearchstructure()
!  free memory for numberings
  call deletedofnumbering(1)
  call deletedofnumbering(2)
  if(vecsize.gt.2)  call deletedofnumbering(vecsize)

!!$  if (myrank.eq.0 .and. maxrank.eq.1) call plote

5003 format(" linear, itaylor, isetup, imask, irestart: ",             &
          5i4, / ," facd, bzero, eps ", 1p3e12.4)
2323 format(1p5e12.4)
2121 format(" start of time dependent loop")
5001 format("ifail .ne. 0 after call to f04aaf")
7011 format(1p6e20.12)
7012 format(6e20.12)

  call safestop(2)

end Program Reducedquintic



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
  use diagnostics
  use gradshafranov

  implicit none

  integer :: l, jer, i
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  integer, allocatable:: itemp(:)
  integer :: ndofs, numnodes

  integer :: calc_matrices

  real :: tstart, tend
  vectype, allocatable :: temp(:), temp2(:)
  
  call numnod(numnodes)

  ! apply loop voltage
  ! ~~~~~~~~~~~~~~~~~~
  fbound = fbound + dt*vloop/(2.*pi)


  ! Determine whether matrices should be re-calculated
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ntime.le.ntimer+1 &
       .or. (linear.eq.0 .and. mod(ntime,nskip).eq.0) &
       .or. (integrator.eq.1 .and. ntime.eq.1)) then
     calc_matrices = 1
  else
     calc_matrices = 0
  endif

  ! calculate matrices for time advance
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(calc_matrices.eq.1) then 
     if(myrank.eq.0 .and. iprint.eq.1) print *, "Defining matrices"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call ludefall
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif
     if(myrank.eq.0 .and. iprint.eq.1) print *, "Done defining matrices."
  endif

  if(isplitstep.eq.1) then
     ! Store current-time velocity matrices for use in field advance
     veln = vel
     veloldn = velold
 
     ! Advance Velocity
     ! ================
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Velocity"
     
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! b1vector = q1matrix_sm * phi(n)
     if(ipres.eq.1 .and. numvar.ge.3) then
        ! replace electron pressure with total pressure
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(numvar, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv   :ibeginnv+11) = phi(ibeginnv:ibeginnv+11)
           phip(ibeginnv+12:ibeginnv+17) = pres(ibegin:ibegin+5)
        enddo
        call matrixvectormult(q1matrix_sm, phip, b1vector)
     else
        call matrixvectormult(q1matrix_sm, phi , b1vector)
     endif
  
     ! vtemp = d1matrix_sm * vel(n)
     vtemp = 0.
     call matrixvectormult(d1matrix_sm,vel,vtemp)
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif

     vtemp = vtemp + b1vector + r4
     
     ! Include linear density terms
     if(idens.eq.1 .or. (gravr.ne.0 .or. gravz.ne.0)) then
        ! b2vector = r14 * den(n)
        
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(numvar, l, 0, ibeginnv, iendplusonenv)
           
           if(idens.eq.1) then
              phip(ibeginnv  :ibeginnv+5) = den(ibegin:ibegin+5)
           else
              phip(ibeginnv) = 1.
           end if
        enddo
        call matrixvectormult(r14matrix_sm,phip,b2vector)
        vtemp = vtemp + b2vector
     endif

     ! apply boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_vel(s1matrix_sm, vtemp)
        call finalizematrix(s1matrix_sm)
     else
        call boundary_vel(0, vtemp)
     endif
     
     ! solve linear system with rhs in vtemp (note LU-decomp done first time)
     if(myrank.eq.0) print *, "solving velocity advance..."
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(s1matrix_sm, vtemp, jer)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_v = t_solve_v + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in velocity solve', jer
        call safestop(42)
     endif
     
     if(integrator.eq.1 .and. ntime.gt.1) then
        vtemp = (2.*vtemp - velold)/3.
     endif

     !.....new velocity solution at time n+1 (or n* for second order advance)
     velold = vel
     vel = vtemp


     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth



     !
     ! Advance Density
     ! ===============
     if(idens.eq.1) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Density"
        
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

        ! b2vector = r8matrix_sm * vel(n+1)
        call matrixvectormult(r8matrix_sm,vel,b2vector)
        
        ! b3vector = q8matrix_sm * vel(n)
        call matrixvectormult(q8matrix_sm,veln,b3vector)
        
        ! b1vector = r8matrix_sm * vel(n-1)
        if(integrator.eq.1 .and. ntime.gt.1) then
           b2vector = 1.5*b2vector
           call matrixvectormult(r8matrix_sm,veloldn,b1vector)
           b1vector = 0.5*b1vector
        else
           b1vector = 0.
        endif

        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_mvm = t_mvm + tend - tstart
        endif

        ! temp = d8matrix_sm * phi(n)
        call createvec(temp, numvar1_numbering)
        temp = 0.

        call matrixvectormult(d8matrix_sm,den,temp)
        
        call numdofs(numvar,ndofs)
        allocate(itemp(ndofs)) ! this is used to make sure that we 
                               ! don't double count the sum for periodic dofs
        itemp = 1

        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(numvar, l, 0, ibeginnv, iendplusonenv)
           do i=0,iendplusone-ibegin-1
              temp(ibegin+i) = temp(ibegin+i) + itemp(ibegin+i) * &
                   (b2vector(ibeginnv+i) + b3vector(ibeginnv+i) + qn4(ibegin+i)) &
                   +b1vector(ibeginnv+i)
              
              itemp(ibegin+i) = 0
           enddo
        enddo
        deallocate(itemp)

        ! Insert boundary conditions
        if(calc_matrices.eq.1) then
           call boundary_den(s8matrix_sm, temp)
           call finalizematrix(s8matrix_sm)
        else
           call boundary_den(0, temp)
        endif
        
        if(myrank.eq.0 .and. iprint.ge.1) print *, " solving..."


        ! solve linear system...LU decomposition done first time
        ! -- okay to here      call printarray(temp, 150, 0, 'vtemp on')
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        call solve(s8matrix_sm, temp, jer)
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_solve_n = t_solve_n + tend - tstart
        endif
        if(jer.ne.0) then
           write(*,*) 'Error in density solve', jer
           call safestop(29)
        endif
        
        ! new field solution at time n+1 (or n* for second order advance)
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = (2.*temp - denold)/3.
        endif
        denold = den
        den = temp
        call deletevec(temp)
     endif
     
     !
     ! Advance Pressure
     ! ================
     if(ipres.eq.1) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure"
        
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        
        ! b2vector = r9matrix_sm * vel(n+1)
        call matrixvectormult(r9matrix_sm,vel,b2vector)
        
        ! b3vector = q9matrix_sm * vel(n)
        call matrixvectormult(q9matrix_sm,veln,b3vector)
        
        ! b1vector = r9matrix_sm * vel(n-1)
        if(integrator.eq.1 .and. ntime.gt.1) then
           b2vector = 1.5*b2vector
           call matrixvectormult(r9matrix_sm,veloldn,b1vector)
           b1vector = 0.5*b1vector
        else
           b1vector = 0.
        endif
        
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_mvm = t_mvm + tend - tstart
        endif
        
        ! temp = d8matrix_sm * pres(n)
        call createvec(temp, numvar1_numbering)
        temp = 0.
        call matrixvectormult(d9matrix_sm,pres,temp)
        
        call numdofs(numvar,ndofs)
        allocate(itemp(ndofs)) ! this is used to make sure that we 
                               ! don't double count the sum for periodic dofs
        itemp = 1
        
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(numvar, l, 0, ibeginnv, iendplusonenv)
           do i=0,iendplusone-ibegin-1
              temp(ibegin+i) = temp(ibegin+i) + itemp(ibegin+i) * &
                   (b2vector(ibeginnv+i) + b3vector(ibeginnv+i) + qp4(ibegin+i)) &
                   +b1vector(ibeginnv+i)
              itemp(ibegin+i) = 0
           enddo
        enddo
        deallocate(itemp)
        
        ! Insert boundary conditions
        if(calc_matrices.eq.1) then
           call boundary_pres(s9matrix_sm, temp)
           call finalizematrix(s9matrix_sm)
        else
           call boundary_pres(0, temp)
        endif
        
        ! solve linear system...LU decomposition done first time
        ! -- okay to here      call printarray(temp, 150, 0, 'vtemp on')
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        call solve(s9matrix_sm, temp, jer)
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_solve_p = t_solve_p + tend - tstart
        endif
        if(jer.ne.0) then
           write(*,*) 'Error in pressure solve', jer
           call safestop(29)
        endif
        
        ! new field solution at time n+1 (or n* for second order advance)
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = (2.*temp - presold)/3.
        endif
        presold = pres
        pres = temp
        call deletevec(temp)
     endif

     
     !
     ! Advance Fields
     ! ==============
     
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Fields"
     
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     
     ! b2vector = r2matrix_sm * vel(n+1)
     call matrixvectormult(r2matrix_sm,vel,b2vector)
     
     ! b3vector = q2matrix_sm * vel(n)
     call matrixvectormult(q2matrix_sm,veln,b3vector)
     
     ! b1vector = r2matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        b2vector = 1.5*b2vector
        call matrixvectormult(r2matrix_sm,veloldn,b1vector)
        b1vector = 0.5*b1vector
     else
        b1vector = 0.
     endif
     
     ! vtemp = d2matrix_sm * phi(n)
     vtemp = 0.
     call matrixvectormult(d2matrix_sm,phi,vtemp)
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
    
     vtemp = vtemp + b2vector + b3vector + q4  + b1vector   

     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_mag(s2matrix_sm, vtemp)
        call finalizematrix(s2matrix_sm)
     else 
        call boundary_mag(0, vtemp)
     endif
     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     call solve(s2matrix_sm, vtemp, jer)

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_b = t_solve_b + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in field solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1 (or n* for second order advance)
     if(integrator.eq.1 .and. ntime.gt.1) then
        vtemp = (2.*vtemp - phiold)/3.
     endif
     phiold = phi
     phi = vtemp
     
  else    
     ! ====================
     ! UNSPLIT TIME ADVANCE
     ! ====================

     if(myrank.eq.0 .and. iprint.ge.1) print *, "Solving matrix equation."
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
         
     ! vtemp = d1matrix_sm * phi(n)
     call matrixvectormult(d1matrix_sm,phi,vtemp)

     vtemp = vtemp + q4

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
     
     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_mag(s1matrix_sm, vtemp)
        call boundary_vel(s1matrix_sm, vtemp)
        if(idens.eq.1) call boundary_den(s1matrix_sm, vtemp)
        if(ipres.eq.1) call boundary_pres(s1matrix_sm, vtemp)
        call finalizematrix(s1matrix_sm)
     else 
        call boundary_mag(0, vtemp)
        call boundary_vel(0, vtemp)
        if(idens.eq.1) call boundary_den(0, vtemp)
        if(ipres.eq.1) call boundary_pres(0, vtemp)
     endif
     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(s1matrix_sm, vtemp, jer)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_b = t_solve_b + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in field solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1
     if(integrator.eq.1 .and. ntime.gt.1) then
        vtemp = (2.*vtemp - phiold)/3.
     endif
     phiold = phi
     phi = vtemp
     
  end if


  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call derived_quantities


  ! Conserve toroidal flux
  ! ~~~~~~~~~~~~~~~~~~~~~~
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     tflux = tflux + gbound*area
  endif

end subroutine onestep


! ======================================================================
! smooth
! ------
!
! applies smoothing operators
!
! ======================================================================
subroutine smooth
  use basic
  use arrays
  use newvar_mod
  use diagnostics

  implicit none

  real :: tstart, tend
  integer :: numnodes

  if(hyperc.eq.0) return

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  call numnod(numnodes)
     
  ! smooth vorticity
  call newvar_d2(vel,vor,1,NV_DCBOUND,NV_GS)
  call smoother1(vor,vel,numnodes,numvar,1)
     
  ! smooth compression
  if(numvar.ge.3) then
     if(com_bc.eq.1) then
        call newvar_d2(vel,com,3,NV_DCBOUND,NV_LP)
     else
        call newvar_d2(vel,com,3,NV_NOBOUND,NV_LP)
     endif
     
!!$        !
!!$        !.....coding to calculate the error in the delsquared chi equation
!!$        call calc_chi_error(chierror)
!!$        if(myrank.eq.0) then
!!$           print *, "Error in com = ", chierror 
!!$        endif
        
     call smoother3(com,vel,numnodes,numvar,3)     
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_smoother = t_smoother + tend - tstart
  endif

end subroutine smooth


! ======================================================================
! derived_quantities
! ------------------
!
! calculates all derived quantities, including auxiliary fields
! and scalars
!
! ======================================================================
subroutine derived_quantities
  use basic
  use arrays
  use newvar_mod
  use diagnostics

  implicit none

  real :: tstart, tend

  ! Define auxiliary fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining auxiliary fields."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  if(myrank.eq.0 .and. iprint.ge.1) print *, "-Transport coefficients"
  call define_transport_coefficients

  !   toroidal current
  if(myrank.eq.0 .and. iprint.ge.1) print *, "-Toroidal current"
  call newvar_d2(psi1_v,jphi,psi_i,NV_DCBOUND,NV_GS)

  if(hyperc.ne.0) then

     !   vorticity
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Vorticity"
     call newvar_d2(phi1_v,vor,phi_i,NV_DCBOUND,NV_GS)

     !   compression
     if(numvar.ge.3) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, "-Compression"
        if(com_bc.eq.1) then
           call newvar_d2(chi1_v,com,chi_i,NV_DCBOUND,NV_LP)
        else
           call newvar_d2(chi1_v,com,chi_i,NV_NOBOUND,NV_LP)
        endif
     else
        com = 0.
     endif
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif
  

  ! find lcfs
  ! ~~~~~~~~~
  call lcfs(psi1_v+psi0_v,vecsize) 


  ! calculate scalars
  ! ~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Calculating scalars"
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call calculate_scalars
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_sources = t_sources + tend - tstart
  endif

end subroutine derived_quantities

! ======================================================================
! conserve_flux
! --------------
!
! adjusts the boundary conditions to conserve toroidal flux
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
     call numdofs(vecsize, ndofs)
     do j=1,ndofs,6*vecsize
        bz1_v(j+bz_off) = bz1_v(j+bz_off) + gbound
     enddo
     if(myrank.eq.0) then
        print *, "Correction to toroidal flux: ", gbound*area
     end if
  endif

  return
end subroutine conserve_flux



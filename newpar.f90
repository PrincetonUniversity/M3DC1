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
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'defining initial conditions...'
     call initial_conditions
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'done initial conditions'

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        fieldi = field
     else
        fieldi = field0
     endif

     ! correct for left-handed coordinates
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "adjusting fields for left-handed coordinates"
     call numdofs(num_fields, ndofs)
     allocate(itemp(ndofs))
     itemp = 1
     do i=1,numnodes
        call entdofs(num_fields, i, 0, ibegin, iendplusone)
        if(itemp(ibegin) .eq. 1) then
           call assign_local_pointers(i)

           psi0_l = -psi0_l
           psi1_l = -psi1_l
           psis_l = -psis_l
             u0_l = -  u0_l
             u1_l = -  u1_l
             us_l = -  us_l
            bz0_l = - bz0_l
            bz1_l = - bz1_l
            bzs_l = - bzs_l
            vz0_l = - vz0_l
            vz1_l = - vz1_l
            vzs_l = - vzs_l

           itemp(ibegin) = 0
        endif
     enddo
     deallocate(itemp)

     ! combine the equilibrium and perturbed fields of linear=0
     ! unless eqsubtract = 1
     if(linear.eq.0 .and. eqsubtract.eq.0) then
        field = field + field0
        field0 = 0.
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
  call newvar_d2(vel,vor,1,vecsize,NV_DCBOUND,NV_GS)
  call smoother1(vor,vel,numnodes,numvar,1)
     
  ! smooth compression
  if(numvar.ge.3) then
     if(com_bc.eq.1) then
        call newvar_d2(vel,com,3,vecsize,NV_DCBOUND,NV_LP)
     else
        call newvar_d2(vel,com,3,vecsize,NV_NOBOUND,NV_LP)
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


subroutine copyvec(inarr, inpos, insize, outarr, outpos, outsize)

  real, intent(in) :: inarr(*)
  integer, intent(in) :: insize, inpos
  real, intent(out) :: outarr(*)
  integer, intent(in) :: outsize, outpos

  integer :: ibegini, iendplusonei
  integer :: ibegino, iendplusoneo
  integer :: l, numnodes, in_i, out_i

  in_i = (inpos-1)*6
  out_i = (output-1)*6

  call numnod(numnodes)

  do l=1,numnodes
     call entdofs(insize, l, 0, ibegini, iendplusonei)
     call entdofs(outsize, l, 0, ibegino, iendplusoneo)
    
     outarr(ibegino+out_i:ibegino+out_i+5) = &
          inarr(ibegini+in_i:ibegini+in_i+5)
  enddo
end subroutine copyvec

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
  call newvar_d2(field,jphi,psi_g,num_fields,NV_DCBOUND,NV_GS)

  if(hyperc.ne.0) then

     !   vorticity
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Vorticity"
     call newvar_d2(field,vor,u_g,num_fields,NV_DCBOUND,NV_GS)

     !   compression
     if(numvar.ge.3) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, "-Compression"
        if(com_bc.eq.1) then
           call newvar_d2(field,com,chi_g,num_fields,NV_DCBOUND,NV_LP)
        else
           call newvar_d2(field,com,chi_g,num_fields,NV_NOBOUND,NV_LP)
        endif
     else
        com = 0.
     endif
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif
  

!!$  ! find lcfs
!!$  ! ~~~~~~~~~
!!$  call lcfs(psi+psi0,1)

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
     call numdofs(1, ndofs)
     do j=1+(bz_g-1)*6,ndofs,6*num_fields
        field(j) = field(j) + gbound
     enddo
     if(myrank.eq.0) then
        print *, "Correction to toroidal flux: ", gbound*area
     end if
  endif

  return
end subroutine conserve_flux



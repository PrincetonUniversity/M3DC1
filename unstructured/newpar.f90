Program Reducedquintic

!   Ref:  [1] Strang and Fix, An Analysis of the Finite Element Method, page 83
!         [2] G.R. Cowper, et al, AIAA Journal, Vol 7, no. 10 page 19
!         [3] JCP Vol 147 p 318-336 (1998)

  use basic
  use arrays
  use newvar_mod
  use sparse
  use hdf5_output
  use diagnostics
  use boundary_conditions
  use time_step
  use m3dc1_output
  use auxiliary_fields
  use pellet
  use scorec_mesh_mod
  use adapt
  use particles
  use math
  use m3dc1_omp
  use restart_hdf5
  use wall
  use geometry 
  use neutral_beam
  use kprad_m3dc1
  use transport_coefficients
  use m3dc1_vel_prof
  use hypervisc
  use runaway_advection
#ifdef _OPENACC
  use openacc
#endif

#if PETSC_VERSION >= 38
  use petsc
  implicit none
#elif PETSC_VERSION >= 36
  implicit none
#include "petsc/finclude/petsc.h"
#else
  implicit none
#include "finclude/petsc.h"
#endif

  integer :: ier, i, adapt_flag
  real :: tstart, tend, dtsave, t_solve, t_compute
  character*10 :: datec, timec
  character*256 :: arg, solveroption_filename
  integer :: ip
  character(len=32) :: mesh_file_name
  logical :: update_mesh

  ! Initialize MPI
#ifdef _OPENMP
  integer :: omp_provided, omp_requested

  omp_requested = MPI_THREAD_SERIALIZED
  call MPI_Init_Thread(omp_requested, omp_provided, ier)
  if(omp_provided .lt. omp_requested) then
     print *, 'Error: MPI implementation does not support required level of OpenMP'
     stop
  end if
  
#else
  call MPI_Init(ier)

#endif

  if (ier /= 0) then
     print *, 'Error in MPI_Init', ier
     call safestop(1)
  endif
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_rank:',ier
     call safestop(1)
  endif
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_size:',ier
     call safestop(1)
  endif

#ifdef _OPENACC
  num_devices = acc_get_num_devices(acc_device_nvidia)
  if (num_devices == 0) num_devices = 1
  igpu=mod(myrank,num_devices)
  !!$acc set device_num(igpu)
  call acc_set_device_num(igpu,acc_device_default)
#endif

  print_help = .false.
  do i=1, command_argument_count()
     call get_command_argument(i, arg)
     if(trim(arg) == '--help') print_help = .true.
#ifdef USE3D
     if(trim(arg) == '-options_file') then
        call get_command_argument(i+1, solveroption_filename)
        !if(myrank==0) print '(2a)', 'solver option file : ', trim(solveroption_filename)
     endif
#endif
  end do

  ! Write version information
  if(myrank.eq.0) then
     print *, '=============================================================='
     print *, 'RELEASE VERSION: ', RELEASE_VERSION
     print *, 'BUILD DATE: ', DATE_BUILT
!     print *, BUILD_INFO
     call date_and_time( datec, timec)
     write(*,1001) datec(1:4),datec(5:6),datec(7:8), &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("START DATE: ", a4,1x,a2,1x,a2,3x,"TIME: ",a2,":",a2,":",a4,/)
#ifdef USECOMPLEX
     print *, 'COMPLEX VERSION'
#else
     print *, 'REAL VERSION'
#endif
#ifdef USE3D
     print *, '3D VERSION'
#else
     print *, '2D VERSION'
#endif
  endif

  ! OPENMP information
#ifdef _OPENMP
!$OMP PARALLEL
     !total number of threads in the group
     nthreads=OMP_GET_NUM_THREADS()
     !id of this thread in the group
     ithread=OMP_GET_THREAD_NUM()
     if(myrank.eq.0) &
     write(*,1004) myrank,ithread,nthreads
1004 format("rank",1x,i4,1x,'has OPENMP: the',1x,i4,1x,'th thread of',1x,i4,1x,'threads in this group.'/)
!$OMP END PARALLEL
#endif

#ifdef USESCOREC
  call m3dc1_domain_init()
#endif

#ifndef M3DC1_TRILINOS
  ! Initialize PETSc
  if(myrank.eq.0) print *, 'Initializing PETSc...'
  call PetscInitialize(PETSC_NULL_CHARACTER, ier)
  if (ier /= 0) then
     print *,'Error in PetscInitialize:',ier
     call safestop(1)
  endif
#endif

!profiling velocity solve
  call PetscLogStageRegister("Mat A",stageA,ier)
  call PetscLogStageRegister("KSP S",stageS,ier)

  ! read input file
  if(myrank.eq.0) print *, ' Reading input'
  call input

!if using SCOREC set adapt verbosity output if iprint.ge.1
  if (iprint.ge.1) then
    call m3dc1_domain_verbosity(1) ! 0 for non-verbose outputs
  end if

  ! load mesh
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Loading mesh nplane='
  !if(myrank==0 .and. nplanes.gt.1) call parse_solver_options(nplanes, trim(solveroption_filename)//PETSC_NULL_CHARACTER)

#ifndef M3DC1_TRILINOS
  call m3dc1_matrix_setassembleoption(imatassemble)
#endif

  call load_mesh

  call print_normal_curv()
!  call print_node_data
!  call safestop(1)

  ! allocate arrays
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Allocating arrays'
  call space(1)

  sparse_initialized = .true.

  ! initialize variables
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Initializing variables'
  call init

  ! output info about simulation to be run
  call print_info

#ifdef USEST 
  if (igeometry.eq.1) then
     ! calculate rst & zst fields
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating geometry'
     call calc_geometry 
     ! recalculate ctri such that DoFs are in physical derivatives 
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Redoing tridef'
     call tridef
     ! recalculate rst & zst fields using new ctri 
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Recalculating geometry'
     call calc_geometry 
  end if
#endif

  ! create the newvar matrices
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Generating newvar matrices'
  call create_newvar_matrices

  call calc_wall_dist

  ! Set initial conditions either from restart file
  ! or from initialization routine
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! initialize output for HDF5 and C1ke
  call initialize_output ()
 
  select case (irestart)
  case(0)
     ! Initialize from routine

     version_in = -1
     ptot = 0.
     ntime = 0
     time = 0.

     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Defining initial conditions'
     call initial_conditions

     ! correct for left-handed coordinates
     if(iflip.eq.1) then
        if(myrank.eq.0 .and. iprint.ge.1) &
             print *, ' Flipping coordinate system handedness'
        call flip_handedness
     endif

     ! initialize feedback systems
     i_control%err_i = 0.
     i_control%err_p_old = 0.
     n_control%err_i = 0.
     n_control%err_p_old = 0.

  case(1)
!
!....save timestep from input file (needed if not a variable timestep run)
     dtsave = dt
     ! Read restart file(s)

     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Reading restart file(s)'
     call rdrestart_hdf5
!
!....use timestep from input file if not a variable timestep run
    if(dtkecrit.eq.0) dt = dtsave

  case(3)
!....read 2D real RL=1 restart file to start 2D complex COM=1 run
!
!....save timestep from input file
     dtsave = dt
     call rdrestart_hdf5
     dt = dtsave

  end select                     !  end of the branch on restart/no restart

  if(irestart.ne.0) then
     call hdf5_reconcile_version(version_in, ier)
  end if

  ntime0 = ntime
  vloop0 = vloop

  ! zero-out scalar data
  call reset_scalars

  mag_probe_itri = 0
  flux_loop_itri = 0
  if(itimer.eq.1) call reset_timings

  ! output simulation parameters
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(irestart.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, " Writing simulation parameters"
     
     call hdf5_write_parameters(ier)
  end if

  ! output equilibrium time slice
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ntime.eq.0 .or. (ntime.eq.ntime0 .and. eqsubtract.eq.1)) then

     if(eqsubtract.eq.1) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
        call find_lcfs()
        call define_transport_coefficients
        call derived_quantities(0)
        if(iwrite_aux_vars.eq.1) call calculate_auxiliary_fields(0)
     end if

     call hdf5_write_time_slice(1,ier)
  end if


  ! combine the equilibrium and perturbed fields of linear=0
  ! unless eqsubtract = 1
  if(eqsubtract.eq.0) then
     call add(field_vec, field0_vec)
     field0_vec = 0.
  endif


  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
  call find_lcfs()
  call define_transport_coefficients
  call derived_quantities(1)

  ! Adapt the mesh
  ! ~~~~~~~~~~~~~~
#ifdef USESCOREC
  ! iadapt_writevtk=1, write the initial mesh before adaptation
  if (iadapt.gt.0 .and. iadapt_writevtk .eq. 1) then
    write(mesh_file_name,"(A7,A)") 'initial', 0
    call m3dc1_mesh_write (mesh_file_name,0,0)
  endif

  !if (iadapt .eq. 1) then
  if (mod (iadapt, 2) .eq. 1) then
    if (iprint.ge.1 .and. myrank.eq.0) write(*,*) "before adapt_by_psi:  psibound, psimin", psibound, psimin
    call adapt_by_psi
  endif
  if (iadapt .eq. 4) call adapt_by_error
#endif

  if(irestart.eq.0  .or. iadapt.gt.0) then
     tflux0 = tflux
  endif

  ! mark the fields necessary for solution transfer
  if (ispradapt .eq. 1) call marker
  
#ifdef USEPARTICLES
  if (kinetic.eq.1) then
     call particle_test
     !call safestop(0)
  endif
#endif

  ! output initial conditions
  call output

  ! if there are no timesteps to calculate, then skip time loop
  if(ntimemax.le.ntime) call safestop(0)

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Initializing timestep'
  call initialize_timestep

  if ((irunaway.ge.1).and.(runaway_characteristics.eq.1)) call runaway_advection_initialize

  ! main time loop
  ! ~~~~~~~~~~~~~~
  do ntime=ntime+1,ntimemax

     if(myrank.eq.0) print *, 'TIME STEP: ', ntime

     call init_hyperv_mat

     ! check for error
     if(ekin.ne.ekin .or. emag.ne.emag) then
        print *, "Error: energy is NaN"
        exit
     endif

     ! re-scale solution if energy is too large
     if(linear.eq.1) call scaleback

     ! take time step
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Calling onestep"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call onestep
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_onestep = t_onestep + tend - tstart
      t_solve = t_solve_v + t_solve_n + t_solve_p + t_solve_b
      t_compute = t_onestep - t_solve
      write(*,1002) ntime, t_onestep, t_compute, t_solve
 1002 format(" LOOP TIME",i5, "   Tot",1pe12.4, "   compute",1pe12.4,"   solve",1pe12.4)
     endif

     if(linear.eq.0 .and. eqsubtract.eq.0) then
        if(i_control%icontrol_type .ge. 0) then
     ! feedback control on toroidal current
          if(tcurf .ne. tcuri) then
          ! time varying target current
            call variable_tcur(tcuri,tcurf,tcur_t0,tcur_tw,time,tcur)
            i_control%target_val = tcur
          endif

          if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " Applying current feedback", &
             vloop, totcur, i_control%p, &
             i_control%target_val, i_control%err_p_old, i_control%err_i

          call control(totcur, vloop,       i_control, dt)

          if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " After current feedback", &
             vloop, totcur, i_control%p, &
             i_control%target_val, i_control%err_p_old, i_control%err_i
       else
          vloop = vloop0*cos(2.*pi*vloop_freq*time)
       endif
     endif

     if(linear.eq.0 .and. eqsubtract.eq.0 .and. n_control%icontrol_type .ge. 0) then
     ! feedback control on density source
          if(myrank.eq.0 .and. iprint.ge.1) print *, " Applying density feedback"
             do ip=1,npellets
                if(myrank.eq.0 .and. (iprint.ge.3 .or. (iprint.ge.1 .and. npellets.eq.1))) then
                   print *, "   ", pellet_rate(ip), totden, n_control%p, &
                        n_control%target_val, n_control%err_p_old, n_control%err_i
                end if
                call control(totden, pellet_rate(ip), n_control, dt) ! ???
             end do

          if(myrank.eq.0 .and. iprint.ge.1) then
             print *, " After density feedback"
             if(iprint.ge.3 .or. npellets.eq.1) then
                do ip=1,npellets
                   print *, "   ", pellet_rate(ip), totden, n_control%p, &
                        n_control%target_val, n_control%err_p_old, n_control%err_i,pellet_state
                end do
             end if
          end if
     endif

     ! Write output
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Writing output."
     call output

    ! for now call spr adapt every 10 time steps
    if (ispradapt .eq. 1) then
      if (mod(ntime, isprntime) .eq. 0) then
        write(mesh_file_name,"(A11,A)") 'beforeadapt', 0
        call m3dc1_mesh_write (mesh_file_name,0,ntime)
        ! if update_mesh is true
        ! the (2D) part.smb will be updated (overwritten) inside adapt_by_spr
        ! the frequency of update is the same as frequency of output time slices
        update_mesh = .false.
        if(mod(ntime-ntime0,ntimepr).eq.0) then
          update_mesh = .true.
        end if
        call adapt_by_spr(field_vec%id, psi_g, ntime, &
             isprweight, isprmaxsize, isprrefinelevel, isprcoarsenlevel, update_mesh)
      endif
    endif

      if (iadapt .gt. 1 .and. ispradapt .eq.0) then
      ! adapt_flag=1 if
      !(1) iadapt_ntime(N)>0 -- run adapt_by_error at the end of every N time steps
      !(2) non-linear & iadapt_ntime=0 -- run adapt_by_error at the end of every time step
      !(3) linear, adapt_ke>0 & ekin>adapt_ke -- run adapt_by_error in this time step  
        call diagnose_adapt(adapt_flag)
       if(adapt_flag .eq. 1) call adapt_by_error
     endif
  enddo ! ntime

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done time loop."

  call safestop(0)

end Program Reducedquintic


!============================================================
! init
! ~~~~
! initialize variables
!============================================================
subroutine init
  use basic
  use mesh_mod
  use basicq
  use runaway_mod
  use kprad_m3dc1
  use resistive_wall
  use pellet
  
  implicit none

  integer :: ierr
 
  rfac = (0,1)*ntor
  if(itor.eq.0) rfac = rfac/rzero

  ! define properties of triangles
  call tridef

  gbound = 0.

  ! Set up PID controllers
  i_control%p = control_p
  i_control%i = control_i
  i_control%d = control_d
  i_control%target_val = tcur
  i_control%icontrol_type = control_type

  n_control%p = n_control_p
  n_control%i = n_control_i
  n_control%d = n_control_d
  n_control%target_val = n_target
  n_control%icontrol_type = n_control_type

  dtold = dt

  call init_qp

  if(irunaway .eq. 2) call runaway_init

  call kprad_init(ierr)
  if(ierr.ne.0) call safestop(601)

  call init_resistive_wall(ierr)
  if(ierr.ne.0) call safestop(602)

  call pellet_domain

end subroutine init


!============================================
! print_info
! ~~~~~~~~~~
! print basic info about simulation options
!============================================
subroutine print_info
  use basic

  implicit none

  ! velocity form
  if(myrank.eq.0) then
     print*, "V = R^2 grad(U)xgrad(phi) + R^2 V grad(phi) + grad(chi)/R^2"
  endif

  ! check time-integration options
  select case(integrator)
  case(1)
     ! For BDF2 integration, first timestep is Crank-Nicholson with thimp=1,
     ! and subsequent timesteps are BDF2.
     if(myrank.eq.0) print *, "Time integration: BDF2."
  case default
     if(myrank.eq.0) print *, "Time integration: Crank-Nicholson."
  end select
end subroutine print_info

!=========================================
! safestop
! ~~~~~~~~
! stop program
!=========================================
subroutine safestop(iarg)

  use basic
  use sparse
  use m3dc1_output
  use time_step
  use auxiliary_fields
  use runaway_mod
  use wall
  use geometry 
  use kprad_m3dc1
  use particles
  use resistive_wall

#if PETSC_VERSION >= 38
  use petsc
  implicit none
#elif PETSC_VERSION >= 36
  implicit none
#include "petsc/finclude/petsc.h"
#else
  implicit none
#include "finclude/petsc.h"
#endif
      
  integer, intent(in) :: iarg
  integer :: ier
  character*10 :: datec, timec

#ifdef USEST
  if (igeometry.eq.1 .and. ilog.ne.-1) then
     if(myrank.eq.0 .and. iprint.ge.2) print *,"  destroying geometry..."
     call destroy_geometry
  end if
#endif

#ifdef USEPARTICLES
  if (kinetic.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.2) print *,"  finalizing particles..."
     call finalize_particles
  endif
#endif

  call destroy_auxiliary_fields
  if(irunaway .eq. 2) call runaway_deallocate
  call kprad_destroy
  call destroy_resistive_wall

  call destroy_wall_dist

  call finalize_timestep

  ! close hdf5 file
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing output..."
  call finalize_output

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  deleting matrices..."
  call delete_matrices

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  unloading mesh..."
  call unload_mesh

#ifndef M3DC1_TRILINOS
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing PETSC..."
  call PetscFinalize(ier)
  if(ier.ne.0) print *, 'Error finalizing PETSC:', ier
#endif

 ! Write time information
  if(myrank.eq.0) then
     print *, '=============================================================='
     call date_and_time( datec, timec)
     write(*,1002) datec(1:4),datec(5:6),datec(7:8), &
          timec(1:2),timec(3:4),timec(5:8)
1002 format(" END DATE: ", a4,1x,a2,1x,a2,3x,"TIME: ",a2,":",a2,":",a4,/)
  endif

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing MPI..."
  call MPI_Finalize(ier)
  if (ier.ne.0) print *, 'Error finalizing MPI:', ier
  
  if(myrank.eq.0) print *, "Stopped at", iarg
  stop
end subroutine safestop


! ======================================================================
! smooth_fields
! ~~~~~~~~~~~~~
!
! applies smoothing operators to fields
! ======================================================================
subroutine smooth_fields(psiin)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  type(field_type), intent(inout) :: psiin

  type(vector_type) :: temp_vec
  type(field_type) :: j_new, psi_new
  real :: tstart, tend

  if(jadv.eq.0 .or. hyper.eq.0. .or. (jadv.eq.1 .and. imp_hyper.ge.1)) return

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' smoothing fields...'
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! smooth current density
  call create_vector(temp_vec,2)
  call associate_field(j_new, temp_vec, 1)
  call associate_field(psi_new, temp_vec, 2)

  j_new = 0.
  psi_new = psiin
  call solve_newvar_axby(s10_mat,temp_vec,d10_mat,temp_vec)
  psiin = psi_new

  call destroy_vector(temp_vec)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_smoother = t_smoother + tend - tstart
  endif
     
end subroutine smooth_fields


!============================================
! find_lcfs
! ~~~~~~~~~
!============================================
subroutine find_lcfs()
  use basic
  use arrays
  use diagnostics
  implicit none

  type(field_type) :: psi_temp, te_temp

  ! Find lcfs
  ! ~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finding lcfs"
  if(eqsubtract.eq.1) then
     if(linear.eq.1) then 
        if(ntime.eq.ntime0) call lcfs(psi_field(0))
     else
if (ispradapt .eq. 1) then
        call create_field(psi_temp, "psi_temp")
else
        call create_field(psi_temp)
endif
        psi_temp = psi_field(0)
        call add_field_to_field(psi_temp, psi_field(1))
        call lcfs(psi_temp)
        call destroy_field(psi_temp)
     endif
  else
     call lcfs(psi_field(1))
  endif
end subroutine find_lcfs

! ======================================================================
! derived_quantities
! ~~~~~~~~~~~~~~~~~~
! calculates all derived quantities, including auxiliary fields
! and scalars
! ======================================================================
subroutine derived_quantities(ilin)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse
  use transport_coefficients
  use auxiliary_fields
  use gradshafranov

  implicit none

  type(field_type) :: psi_temp, te_temp
  integer, intent(in) :: ilin    ! 0 for equilibrium fields, 1 for perturbed
  integer :: ier
  real :: tstart, tend

  vectype :: temp

  ! Find maximum temperature:  te_max
  ! ~~~~~~~~~
  ier = 0
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finding temax"
  if(eqsubtract.eq.1) then
     if(linear.eq.1) then 
        if(ntime.eq.ntime0) then
          if(ifixed_temax .eq. 0) then
            call te_max(xmag,zmag,te_field(0),temax,0,ier)
          else
            call te_max2(xmag0,zmag0,te_field(0),temax,0,ier)
          endif
        endif
     else
if (ispradapt .eq. 1) then
        call create_field(te_temp, "te_temp")
else
        call create_field(te_temp)
endif
        te_temp = te_field(0)
        call add_field_to_field(te_temp, te_field(1))
        if(ifixed_temax .eq. 0) then
           call te_max(xmag,zmag,te_temp,temax,0,ier)
        else
           call te_max2(xmag0,zmag0,te_temp,temax,0,ier)
        endif
        call destroy_field(te_temp)
     endif
  else
     if(ifixed_temax .eq. 0) then
       call te_max(xmag,zmag,te_field(1),temax,0,ier)
     else
       call te_max2(xmag0,zmag0,te_field(1),temax,0,ier)
     endif
  endif
  if(ier.eq.0) then
    if(myrank.eq.0 .and. iprint.ge.1) then
      write(*,'(A, E12.4)') 'max te', temax
    endif
  else
    if(myrank.eq.0 .and. iprint.ge.1) then
      write(*,'(A,2e12.4)') ' no temperatue maximum found near ',xmag,zmag
    endif
  endif

  ! Electron density
  call calculate_ne(ilin, den_field(ilin), ne_field(ilin), eqsubtract)


  ! Define auxiliary fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  if(itemp.eq.0 .and. (numvar.eq.3 .or. ipres.gt.0) .and. imp_temp.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  temperatures"
     call calculate_temperatures(ilin, te_field(ilin), ti_field(ilin), &
          pe_field(ilin), p_field(ilin), ne_field(ilin), den_field(ilin), &
          eqsubtract)
  end if

  !   toroidal current
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  toroidal current"
  if(inocurrent_tor.eq.1) then
     call solve_newvar1(mass_mat_lhs_dc,jphi_field,gs_mat_rhs_dc, &
          psi_field(ilin))
  else
     call solve_newvar1(mass_mat_lhs,jphi_field,gs_mat_rhs, &
          psi_field(ilin))
  endif

  
  ! vector potential stream function
  !if(imp_bf.eq.0 .or. ilin.eq.0 .or. ntime.eq.0) then
     if((i3d.eq.1 .or. ifout.eq.1) .and. numvar.ge.2) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, "  f", ilin
        if((ilin.eq.0 .and. eqsubtract.eq.1) &
            .or. eqsubtract.eq.0) then
           if(itor.eq.0) then
              temp = bzero
           else
              temp = bzero*rzero
           end if
           call add(bz_field(ilin),-temp)
        endif
        call solve_newvar1(bf_mat_lhs,bf_field(ilin),mass_mat_rhs_bf, &
             bz_field(ilin), bf_field(ilin))
        if((ilin.eq.0 .and. eqsubtract.eq.1) &
            .or. eqsubtract.eq.0) call add(bz_field(ilin), temp)
     endif
  !end if

  ! toroidal derivative of vector potential stream function
  if(imp_bf.eq.0 .or. ntime.eq.ntime0) then
     if(i3d.eq.1 .and. ilin.eq.1 .and. numvar.ge.2) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, "  fp", ilin
        ! solve fp = df/dphi when restarting absent fp 
        if(irestart_fp.eq.0 .and. ntime.eq.ntime0) then 
           call solve_newvar1(mass_mat_lhs,bfp_field(ilin),dp_mat_rhs_bfp, &
               bf_field(ilin))
           if(extsubtract.eq.1) then ! also, external bfp
              call solve_newvar1(mass_mat_lhs,bfp_ext,dp_mat_rhs_bfp, &
                  bf_ext)
           endif
        else 
           call solve_newvar1(bf_mat_lhs,bfp_field(ilin),dp_mat_rhs_bfp, &
                bz_field(ilin), bfp_field(ilin))
        endif
     endif
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif

  ! calculate scalars
  ! ~~~~~~~~~~~~~~~~~
  if(icalc_scalars.eq.1 .and. ilin.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  scalars"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call calculate_scalars
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_sources = t_sources + tend - tstart
     endif
  end if

end subroutine derived_quantities

!======================================================================
! conserve_flux
! ~~~~~~~~~~~~~
! adjusts the boundary conditions to conserve toroidal flux
!======================================================================
subroutine conserve_flux

  use basic
  use arrays
  use diagnostics

  implicit none
  
  ! adjust toroidal field and boundary condition to conserve toroidal flux
  if(numvar.ge.2 .and. iconstflux.eq.1) then
     gbound = (tflux0-tflux)/area
     call add(bz_field(1), gbound)
     if(myrank.eq.0) then
        print *, "Correction to toroidal flux: ", gbound*area
     end if
  endif

  return
end subroutine conserve_flux


!======================================================================
! flip_handedness
! ~~~~~~~~~~~~~~~
! Flips coordinate system handedness by flipping sign of
! psi, u, vz, and bz
!======================================================================
subroutine flip_handedness

  use basic
  use arrays
  use diagnostics

  implicit none

  integer :: ilin

  vectype, parameter :: temp = -1

  do ilin=0,1
     call mult(psi_field(ilin), temp)
     call mult(bz_field(ilin),  temp)
     call mult(u_field(ilin),   temp)
     call mult(vz_field(ilin),  temp)
     if(icsubtract.eq.1) call mult(psi_coil_field, temp)
  end do

  psimin = -psimin
  psilim = -psilim
  psibound = -psibound

end subroutine flip_handedness


!============================================================
! rotation
! ~~~~~~~~
! calculates the rotation matrix rot given angle theta
!============================================================
subroutine rotation(rot,ndim,theta)
  implicit none

  integer, intent(in) :: ndim
  real, intent(in) :: theta
  real, intent(out) :: rot(ndim,*)

  integer :: i, j
  real :: r1(6,6), co, sn

  co = cos(theta)
  sn = sin(theta)
  do i=1,6
     do j=1,6
        r1(i,j) = 0.
     enddo
  enddo

  r1(1,1) = 1.

  r1(2,2) = co
  r1(2,3) = sn

  r1(3,2) = -sn
  r1(3,3) = co

  r1(4,4) = co**2
  r1(4,5) = 2.*sn*co
  r1(4,6) = sn**2

  r1(5,4) = -sn*co
  r1(5,5) = co**2-sn**2
  r1(5,6) = sn*co

  r1(6,4) = sn**2
  r1(6,5) = -2.*sn*co
  r1(6,6) = co**2

  do i=1,18
     do j=1,18
        rot(i,j) = 0.
     enddo
  enddo

  do i=1,6
     do j=1,6
        rot(i,j)       = r1(i,j)
        rot(i+6,j+6)   = r1(i,j)
        rot(i+12,j+12) = r1(i,j)
     enddo
  enddo

  return
end subroutine rotation


  !============================================================
  ! tridef
  ! ~~~~~~
  ! populates the *tri arrays
  !============================================================
  subroutine tridef
    use basic
    use math
    use mesh_mod

    implicit none

    include 'mpif.h'
  
    type(element_data) :: d
    integer :: itri, i, j, k, ii, jj, numelms, numnodes, ndofs, ierr
    real, dimension(coeffs_per_tri,coeffs_per_tri) :: ti, t
    real, dimension(dofs_per_tri, dofs_per_tri) :: rot
#ifdef USEST
    real, dimension(dofs_per_element, dofs_per_element) :: newrot
    integer :: l, m, n, idof, ip
!    real, dimension(dofs_per_node, dofs_per_node) :: rot1, rot2, rot3
#else
    real, dimension(dofs_per_tri, dofs_per_tri) :: newrot
#endif
    real :: sum, theta, mean_area, tot_area, mean_len
    real :: norm(2), curv(3), x, phi, z
    integer :: inode(nodes_per_element)
    logical :: is_boundary
    integer :: izone, izonedim
    integer :: tot_elms
    integer :: info1, info2
    real :: wkspce(9400)
    integer :: ipiv(coeffs_per_tri)

!    real, dimension(dofs_per_tri) :: temp_vec

    real, dimension(nodes_per_element) :: node_sz
#ifdef USEST
    real, dimension(dofs_per_element, dofs_per_element) :: p2l_mat 
#endif

    numelms = local_elements()
    numnodes = local_nodes()
    ndofs = numnodes*dofs_per_node

    ! start the loop over triangles within a rectangular region
    do itri=1,numelms

       ! define a,b,c and theta
       call get_element_data(itri,d)
       
       ! define the Inverse Transformation Matrix that enforces the 
       ! condition that the normal slope between triangles has only 
       ! cubic variation
       call tmatrix(t,d%a,d%b,d%c)
       
       ! calculate the inverse of t using lapack routines
       info1 = 0
       info2 = 0
       ti = t
       call dgetrf(coeffs_per_tri,coeffs_per_tri,ti,coeffs_per_tri,ipiv,info1)
       call dgetri(coeffs_per_tri,ti,coeffs_per_tri,ipiv,wkspce,400,info2)
       if(info1.ne.0.or.info2.ne.0) write(*,'(3I5)') info1,info2
    
       ! calculate the rotation matrix rot
       theta = atan2(d%sn,d%co)
       call rotation(rot,dofs_per_tri,theta)
       
       newrot = 0.
       call get_element_nodes(itri, inode)

#ifdef USEST
       p2l_mat=0.
       do i=1, nodes_per_element 
          k = (i-1)*dofs_per_node + 1

          if(igeometry.eq.1.and.ilog.eq.2) then
             call p2l_matrix(&
                  p2l_mat(k:(k+dofs_per_node-1), k:(k+dofs_per_node-1)),&
                  inode(i))
          end if
#else
       do i=1, 3
          k = (i-1)*6 + 1
#endif           
          call boundary_node(inode(i), &
               is_boundary, izone, izonedim, norm, curv, x, phi, z, &
               BOUND_ANY)
          if(is_boundary) then
#ifdef USEST
             call newrot_matrix(&
                  newrot(k:(k+dofs_per_node-1),k:(k+dofs_per_node-1)),&
                  norm,curv,2)
          else
             do j=1, dofs_per_node
#else
             newrot(k  ,k  ) = 1.
             newrot(k+1,k+1) =  norm(1)
             newrot(k+1,k+2) =  norm(2)
             newrot(k+1,k+3) =  curv(1)*norm(2)**2
             newrot(k+1,k+4) = -curv(1)*norm(1)*norm(2)
             newrot(k+1,k+5) =  curv(1)*norm(1)**2
             newrot(k+2,k+1) = -norm(2)
             newrot(k+2,k+2) =  norm(1)
             newrot(k+2,k+3) =  2.*curv(1)*norm(1)*norm(2)
             newrot(k+2,k+4) = -curv(1)*(norm(1)**2 - norm(2)**2) 
             newrot(k+2,k+5) = -2.*curv(1)*norm(1)*norm(2)
             newrot(k+3,k+3) =  norm(1)**2 
             newrot(k+3,k+4) =  norm(1)*norm(2)
             newrot(k+3,k+5) =  norm(2)**2
             newrot(k+4,k+3) = -2.*norm(1)*norm(2)
             newrot(k+4,k+4) =  norm(1)**2 - norm(2)**2
             newrot(k+4,k+5) =  2.*norm(1)*norm(2)
             newrot(k+5,k+3) =  norm(2)**2
             newrot(k+5,k+4) = -norm(1)*norm(2)
             newrot(k+5,k+5) =  norm(1)**2
          else
             do j=1, 6
#endif
                newrot(k+j-1,k+j-1) = 1.
             end do
          end if
       end do     
#ifdef USEST
       ! newrot now includes transformation between physical & logical DoFs
       if(igeometry.eq.1.and.ilog.eq.2) then
          newrot = matmul(newrot,transpose(p2l_mat))
       end if
       
       gtri(:,:,itri) = matmul(ti(:,1:dofs_per_tri),rot)

#else
       ! form the matrix g by multiplying ti and rot
       do k=1, coeffs_per_tri
          do j=1, dofs_per_tri
             sum = 0.
             do ii = 1, dofs_per_tri
                do jj=1, dofs_per_tri
                   sum = sum + newrot(j,jj)*ti(k,ii)*rot(ii,jj)
                enddo
             enddo
             gtri(k,j,itri) = sum
          enddo
       enddo
#endif

       htri(1,1,itri) = 1.
#ifdef USE3D
       htri(2,1,itri) = 0.
       htri(3,1,itri) =-3./d%d**2
       htri(4,1,itri) = 2./d%d**3

       htri(1,2,itri) = 0.
       htri(2,2,itri) = 1.
       htri(3,2,itri) =-2./d%d
       htri(4,2,itri) = 1./d%d**2

       htri(1,3,itri) = 0.
       htri(2,3,itri) = 0.
       htri(3,3,itri) = 3./d%d**2
       htri(4,3,itri) =-2./d%d**3

       htri(1,4,itri) = 0.
       htri(2,4,itri) = 0.
       htri(3,4,itri) =-1./d%d
       htri(4,4,itri) = 1./d%d**2
#endif
       if(iprecompute_metric.eq.1) then
          call local_coeff_vector(itri,ctri(:,:,itri))
#ifdef USEST
          if((igeometry.eq.1.and.ilog.eq.2).or.igeometry.eq.0) then
             ctri(:,:,itri) = matmul((newrot),ctri(:,:,itri))  
          end if
#endif
       end if
    end do
#ifdef USEST
    if((igeometry.eq.1.and.ilog.eq.2).or.igeometry.eq.0) then
       if(iprecompute_metric.eq.1) then
          deallocate(gtri,htri)
       end if
    end if
#endif

    select case(equilibrate)
    case(1)
       tot_area = 0.
       do itri=1, numelms
          call get_element_data(itri,d)
          tot_area = tot_area + (d%c)*(d%a + d%b)/2.
       end do

       call mpi_allreduce(numelms, tot_elms, 1, MPI_INTEGER, &
            MPI_SUM, MPI_COMM_WORLD, ierr)      
       call mpi_allreduce(tot_area, mean_area, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
       if(nplanes.le.1) then 
          mean_len = 1.
       else
          mean_len = 2.*pi/(nplanes-1)
       endif
       
       if(myrank.eq.0 .and. iprint.ge.1) then 
          print *, ' Total mesh area: ', mean_area
          print *, ' Total elements: ', tot_elms
       endif
       
       mean_area = mean_area/tot_elms
    
       if(myrank.eq.0 .and. iprint.ge.1) then 
          print *, ' Average area: ', mean_area
          print *, ' 1/mean_area = ', 1./mean_area
          print *, ' 1/mean_area**2 = ', 1./mean_area**2
       end if
       
       do i=1, pol_nodes_per_element
          equil_fac(1+dofs_per_node*(i-1),:) = 1./mean_area
          equil_fac(2+dofs_per_node*(i-1),:) = 1./sqrt(mean_area**3)
          equil_fac(3+dofs_per_node*(i-1),:) = 1./sqrt(mean_area**3)
          equil_fac(4+dofs_per_node*(i-1),:) = 1./mean_area**2
          equil_fac(5+dofs_per_node*(i-1),:) = 1./mean_area**2
          equil_fac(6+dofs_per_node*(i-1),:) = 1./mean_area**2
#ifdef USE3D
          do j=1, 6
             equil_fac(j+6+dofs_per_node*(i-1),:) = &
                  equil_fac(j+dofs_per_node*(i-1),:)/mean_len**2
          end do
#endif
       end do
       
    case(2)

!       temp_vec = 0.
       do itri=1, numelms
#ifdef USESCOREC
          !call getelmsizes(itri, node_sz)
          node_sz = 1.
#else
          node_sz = 1.
#endif
          if(myrank.eq.0 .and. itri.eq.1) print *, 'node_sz = ', node_sz

          do i=1, nodes_per_element
             equil_fac(1+dofs_per_node*(i-1),:) = 1./node_sz(i)**2
             equil_fac(2+dofs_per_node*(i-1),:) = 1./node_sz(i)**3
             equil_fac(3+dofs_per_node*(i-1),:) = 1./node_sz(i)**3
             equil_fac(4+dofs_per_node*(i-1),:) = 1./node_sz(i)**4
             equil_fac(5+dofs_per_node*(i-1),:) = 1./node_sz(i)**4
             equil_fac(6+dofs_per_node*(i-1),:) = 1./node_sz(i)**4
          end do
          if(myrank.eq.0 .and. itri.eq.1) print *, equil_fac(:,itri)
       end do
    end select
  end subroutine tridef

!============================================================
! space
! ~~~~~
! allocates space for big arrays
!
! ifirstcall = 1 if this is the first call and 0 otherwise
!============================================================
subroutine space(ifirstcall)

  use element
  use mesh_mod
  use basic
  use arrays
  use sparse
  use time_step
  use auxiliary_fields
  use transport_coefficients

  implicit none

  integer, intent(in) :: ifirstcall

  integer :: numelms

#ifdef USESCOREC
  integer :: i
  character(len=32) :: field_name
#endif

  if (myrank.eq.0 .and. iprint.ge.1) print *, " Entering space..."

!.....create numberings
#ifdef USESCOREC
  if(ifirstcall .eq. 1) then
     do i=1, num_fields
if (ispradapt .eq. 1) then
       write(field_name,"(A3,I0,A)")  "mat", i, 0
else
       write(field_name,"(I2,A)")  i,0
endif

#ifdef USECOMPLEX
       call m3dc1_field_create (i, trim(field_name), i, 1, dofs_per_node)
#else
       call m3dc1_field_create (i, trim(field_name), i, 0, dofs_per_node)
#endif
     end do
  endif ! on firstcall
#endif
  
  numelms = local_elements()

! arrays defined at all vertices
! createvec will delete the arrays if they have already been allocated
  if(ifirstcall.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'Allocating...'

if (ispradapt .eq. 1) then
     ! Physical Variables
     call create_vector(field_vec , num_fields, "field_vec")
     call create_vector(field0_vec, num_fields, "field_vec0")
     !if(iadapt .ne. 0) then
        call create_vector(field_vec_pre, 2, "field_vec_pre")
     !end if

     call mark_vector_for_solutiontransfer(field_vec)
     call mark_vector_for_solutiontransfer(field0_vec)
     call mark_vector_for_solutiontransfer(field_vec_pre)

  ! Auxiliary Variables
     call create_field(jphi_field, "jphi")
     call create_field(resistivity_field, "resistivity")
     call create_field(kappa_field, "kappa")
     call create_field(kappar_field, "kappar")
     call create_field(denm_field, "denm")
     call create_field(visc_field, "visc")
     call create_field(visc_c_field, "visc_c")
     if(ipforce.gt.0) call create_field(pforce_field, "pforce")
     if(ipforce.gt.0) call create_field(pmach_field, "pmach")
     if(density_source) call create_field(sigma_field, "sigma")
     if(momentum_source) call create_field(Fphi_field, "Fphi")
     if(heat_source) call create_field(Q_field, "Q")
     if(icd_source.gt.0) call create_field(cd_field, "cd")
     if(rad_source) then
        call create_field(Totrad_field, "Torad")
        call create_field(Linerad_field, "Linerad")
        call create_field(Bremrad_field, "Bremrad")
        call create_field(Ionrad_field, "Ionrad")
        call create_field(Reckrad_field, "Reckrad")
        call create_field(Recprad_field, "Recprad")
     end if
     call create_field(bf_field(0), "bf0")
     call create_field(bf_field(1), "bf1")
     call create_field(bfp_field(0), "bfp0")
     call create_field(bfp_field(1), "bfp1")
     if(ibootstrap.gt.0) call create_field(visc_e_field, "visc_e")
     !Bootstrap Coeff Fields
     if(ibootstrap.gt.0) call create_field(Jbs_L31_field, "Jbs_L31")
     if(ibootstrap.gt.0) call create_field(Jbs_L32_field, "Jbs_L32")
     if(ibootstrap.gt.0) call create_field(Jbs_L34_field, "Jbs_L34")
     if(ibootstrap.gt.0) call create_field(Jbs_alpha_field, "Jbs_alpha")
     if(ibootstrap.gt.0) call create_field(Jbs_fluxavg_iBsq_field, "Jbs_fluxavg_iBsq")
     if(ibootstrap.gt.0) call create_field(Jbs_fluxavg_G_field, "Jbs_fluxavg_G")
     if(ibootstrap.eq.2) call create_field(Jbs_dtedpsit_field, "Jbs_dtedpsit")

     call create_field(psi_coil_field, "psi_coil")

     ! create external fields
     if(extsubtract.eq.1) then
        call create_field(psi_ext, "pis_ext")
        call create_field(bz_ext, "bz_ext")
        call create_field(bf_ext, "bf_ext")
        call create_field(bfp_ext, "bfp_ext")
        use_external_fields = .true.
     end if
else
     ! Physical Variables
     call create_vector(field_vec , num_fields)
     call create_vector(field0_vec, num_fields)
     !if(iadapt .ne. 0) then
        call create_vector(field_vec_pre, 2)
     !end if

     ! Auxiliary Variables
     call create_field(jphi_field)
     call create_field(resistivity_field)
     call create_field(kappa_field)
     call create_field(kappar_field)
     call create_field(denm_field)
     call create_field(visc_field)
     call create_field(visc_c_field)
     if(ipforce.gt.0) call create_field(pforce_field)
     if(ipforce.gt.0) call create_field(pmach_field)
     if(density_source) call create_field(sigma_field)
     if(momentum_source) call create_field(Fphi_field)
     if(heat_source) call create_field(Q_field)
     if(icd_source.gt.0) call create_field(cd_field)
     if(rad_source) then
        call create_field(Totrad_field)
        call create_field(Linerad_field)
        call create_field(Bremrad_field)
        call create_field(Ionrad_field)
        call create_field(Reckrad_field)
        call create_field(Recprad_field)
     end if
     call create_field(bf_field(0))
     call create_field(bf_field(1))
     call create_field(bfp_field(0))
     call create_field(bfp_field(1))
     if(ibootstrap.gt.0) call create_field(visc_e_field)
     if(ibootstrap.gt.0) call create_field(Jbs_L31_field)
     if(ibootstrap.gt.0) call create_field(Jbs_L32_field)
     if(ibootstrap.gt.0) call create_field(Jbs_L34_field)
     if(ibootstrap.gt.0) call create_field(Jbs_alpha_field)
     if(ibootstrap.gt.0) call create_field(Jbs_fluxavg_iBsq_field)
     if(ibootstrap.gt.0) call create_field(Jbs_fluxavg_G_field)
     if(ibootstrap.eq.2) call create_field(Jbs_dtedpsit_field)

     call create_field(psi_coil_field)

     ! create external fields
     if(extsubtract.eq.1) then
        call create_field(psi_ext)
        call create_field(bz_ext)
        call create_field(bf_ext)
        call create_field(bfp_ext)
        use_external_fields = .true.
     end if
endif

#ifdef USEPARTICLES
     call create_field(p_f_par)
     call create_field(p_f_perp)
     call create_field(den_f_0)
     call create_field(den_f_1)
     call create_field(v_f_par)
     call create_field(p_i_par)
     call create_field(p_i_perp)
     call create_field(den_i_0)
     call create_field(den_i_1)
     call create_field(v_i_par)
     call create_field(rho_field)
     call create_field(nf_field)
     call create_field(tf_field)
     call create_field(pf_field)
     call create_field(nfi_field)
     call create_field(tfi_field)
     call create_field(pfi_field)
     call create_field(densmooth_field)
     call create_field(vparsmooth_field)
#endif

     call create_auxiliary_fields
  endif

  ! arrays associated with the triangles
  if(ifirstcall.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' deallocating...'
     deallocate(gtri,htri)
     if(iprecompute_metric.eq.1) deallocate(ctri)
     if(equilibrate.ne.0) deallocate(equil_fac)
  endif
  
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Allocating tri...'
  allocate(gtri(coeffs_per_tri,dofs_per_tri,numelms))
  allocate(htri(coeffs_per_dphi,dofs_per_dphi,numelms))
  if(iprecompute_metric.eq.1) &
       allocate(ctri(dofs_per_element,coeffs_per_element,numelms))
  if(equilibrate.ne.0) allocate(equil_fac(dofs_per_element,numelms))

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' associating...'
  call associate_field(u_field(1),   field_vec, u_g)
  call associate_field(vz_field(1),  field_vec, vz_g)
  call associate_field(chi_field(1), field_vec, chi_g)
  call associate_field(psi_field(1), field_vec, psi_g)
  call associate_field(bz_field(1),  field_vec, bz_g)
  call associate_field(pe_field(1),  field_vec, pe_g)
  call associate_field(den_field(1), field_vec, den_g)
  call associate_field(p_field(1),   field_vec, p_g)
  call associate_field(te_field(1),  field_vec, te_g)
  call associate_field(ti_field(1),  field_vec, ti_g)
  call associate_field(e_field(1),   field_vec, e_g)
  call associate_field(ne_field(1),  field_vec, ne_g)
  call associate_field(nre_field(1), field_vec, nre_g)

  call associate_field(u_field(0),   field0_vec, u_g)
  call associate_field(vz_field(0),  field0_vec, vz_g)
  call associate_field(chi_field(0), field0_vec, chi_g)
  call associate_field(psi_field(0), field0_vec, psi_g)
  call associate_field(bz_field(0),  field0_vec, bz_g)
  call associate_field(pe_field(0),  field0_vec, pe_g)
  call associate_field(den_field(0), field0_vec, den_g)
  call associate_field(p_field(0),   field0_vec, p_g)
  call associate_field(te_field(0),  field0_vec, te_g)
  call associate_field(ti_field(0),  field0_vec, ti_g)
  call associate_field(e_field(0),   field0_vec, e_g )
  call associate_field(ne_field(0),  field0_vec, ne_g)
  call associate_field(nre_field(0), field0_vec, nre_g)

  call allocate_kspits

  !if (iadapt .ne. 0)  then
     call associate_field(u_field_pre,   field_vec_pre, u_g)
     call associate_field(psi_field_pre, field_vec_pre, psi_g)
  !end if


  if(myrank.eq.0 .and. iprint.ge.1) print *, " Exiting space."

  return
end subroutine space

subroutine calculate_zeff(itri, z)
  use basic
  use kprad
  use kprad_m3dc1
  use m3dc1_nint

  implicit none

  integer, intent(in) :: itri

  integer :: i
  vectype, dimension(MAX_PTS), intent(out) :: z

  z = z_ion**2*nt79(:,OP_1)

  if(ikprad.ne.0) then 
     do i=1, kprad_z
        call eval_ops(itri, kprad_n(i), tm79, rfac)
        z = z + i**2*tm79(:,OP_1)
     end do
  end if

  z = z / net79(:,OP_1)
end subroutine calculate_zeff

! Calculate the factor that multiplies ne * (Ti - Te) in Q_Delta term
subroutine calculate_qdfac(itri, z)
  use basic
  use kprad
  use kprad_m3dc1
  use m3dc1_nint

  implicit none

  integer, intent(in) :: itri
  vectype, dimension(MAX_PTS), intent(out) :: z

  integer :: i

  z = z_ion**2 * nt79(:,OP_1) / ion_mass

  if(ikprad.ne.0) then 
     do i=1, kprad_z
        call eval_ops(itri, kprad_n(i), tm79, rfac)
        z = z + i**2 * tm79(:,OP_1) / kprad_mz
     end do
  end if
  temp79a = max(temin_qd,real(tet79(:,OP_1)))
  where(real(temp79a).gt.0.)
     z = z * 3. * me_mp * nufac / temp79a**(3./2.)
  elsewhere
     z = 0.
  end where
end subroutine calculate_qdfac

subroutine print_normal_curv()
!  use mpi
  use basic
  use mesh_mod
  implicit none

    include 'mpif.h'

  integer :: ierr, i, icounter_t, numnodes
  integer :: izone, izonedim
  real :: norm(2), curv(3), x, z, phi
  logical :: is_boundary
  character(len=255) :: buf
  character(len=:), allocatable :: local_buffer, remote_buffer
  integer :: local_length
  integer :: tag
  integer :: ifile
  integer :: source_rank
  integer :: remote_length
  integer, dimension(MPI_STATUS_SIZE) :: status

  ! Each process builds its own diagnostic output in a local buffer.
  local_buffer = ""
  numnodes = owned_nodes()
  do icounter_t = 1, numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i, is_boundary, izone, izonedim, norm, curv, x, phi, z, BOUND_ANY)
     if (.not. is_boundary) cycle
     buf = ""
     write(buf, '(7F10.5,A)') x, z, norm(1), norm(2), curv(1), curv(2), curv(3), new_line('a')
     local_buffer = local_buffer // trim(buf)
  end do
  local_length = len_trim(local_buffer)

  tag = 0
  if (myrank == 0) then
     ! Rank 0 performs all file I/O.
     call MPI_File_open(MPI_COMM_SELF, "normcurv", ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, ifile, ierr)
     if (local_length > 0) then
        call MPI_File_write(ifile, local_buffer, local_length, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
     end if

     ! Receive output from other ranks one at a time (in rank order).
     do source_rank = 1, maxrank - 1
        call MPI_Recv(remote_length, 1, MPI_INTEGER, source_rank, tag, MPI_COMM_WORLD, status, ierr)
        if (remote_length > 0) then
           allocate(character(len=remote_length) :: remote_buffer)
           call MPI_Recv(remote_buffer, remote_length, MPI_CHARACTER, source_rank, tag, MPI_COMM_WORLD, status, ierr)
           call MPI_File_write(ifile, remote_buffer, remote_length, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
           deallocate(remote_buffer)
        end if
     end do
     call MPI_File_close(ifile, ierr)
  else
     ! All nonzero ranks send their buffer length and then the buffer.
     call MPI_Send(local_length, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, ierr)
     if (local_length > 0) then
        call MPI_Send(local_buffer, local_length, MPI_CHARACTER, 0, tag, MPI_COMM_WORLD, ierr)
     end if
  end if

end subroutine print_normal_curv

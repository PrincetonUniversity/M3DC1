module time_step
  use time_step_split
  use time_step_unsplit
  implicit none

  integer :: meshAdapted
  data meshAdapted /0/

contains

  subroutine initialize_timestep
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call initialize_timestep_unsplit
       call assign_variables_unsplit
    case(1:2)
       call initialize_timestep_split
       call assign_variables_split
    end select

    if(myrank.eq.0 .and. iprint.ge.1) then
       print *, 'Index of U: ', u_i
       print *, 'Index of V: ', vz_i
       print *, 'Index of Chi: ', chi_i
       print *, 'Index of Psi: ', psi_i
       print *, 'Index of Bz: ', bz_i
       print *, 'Index of P: ', p_i
       print *, 'Index of n: ', den_i
       print *, 'Index of Pe: ', pe_i
       print *, 'Index of f: ', bf_i
       print *, 'Index of E: ', e_i
    end if
  end subroutine initialize_timestep

  subroutine finalize_timestep
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call finalize_timestep_unsplit
    case(1:2)
       call finalize_timestep_split
    end select
  end subroutine finalize_timestep

  subroutine clear_matrices
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call clear_matrices_unsplit
    case(1:2)
       call clear_matrices_split
    end select
  end subroutine clear_matrices
  
  subroutine finalize_matrices
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call finalize_matrices_unsplit
    case(1:2)
       call finalize_matrices_split
    end select
  end subroutine finalize_matrices


!============================================================
! ONESTEP
! ~~~~~~~
! advance solution to time ntime+1
!============================================================
subroutine onestep

  use basic
  use diagnostics
  use arrays
  use pellet
  use runaway_mod
  use kprad_m3dc1
  use transport_coefficients
  use particles
  use boundary_conditions

  implicit none

  integer :: calc_matrices, ivel_def, irepeat
  logical, save :: first_time = .true.

  real :: tstart, tend

#ifdef USE3D
  integer :: icount, maxiter
#endif

  ! Determine whether matrices should be re-calculated
  if(first_time &
       .or. (linear.eq.0 .and. mod(ntime,nskip).eq.0) &
       .or. (integrator.eq.1 .and. ntime.eq.1) .or. meshAdapted .eq. 1) then
     calc_matrices = 1
  else
     calc_matrices = 0
  endif

  ! Advance impurity charge states
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call kprad_ionize(dt)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_kprad = t_kprad + tend - tstart
  endif

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
  call define_transport_coefficients

  ! start of loop to repeat timestep if max iterations exceeded in 3D
  do irepeat = 1, max_repeat

    ! calculate matrices for time advance
    if(calc_matrices.eq.1) then
       if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining matrices"
       if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

       ! in linear case, eliminate second-order terms from matrix
       ivel_def = 1
       if(istatic.eq.1 .and. isplitstep.eq.1) ivel_def = 0
       call ludefall(ivel_def, idens, ipres, ipressplit, 1-iestatic)
       if(myrank.eq.0 .and. itimer.eq.1) then
          call second(tend)
          t_ludefall = t_ludefall + tend - tstart
       endif
       if(myrank.eq.0 .and. iprint.ge.1) print *, "Done defining matrices."
    endif

#ifdef USEPARTICLES
    if ((kinetic .eq. 1) .and. (linear .eq. 1)) then
       if (myrank .eq. 0 .and. itimer .eq. 1) call second(tstart)
       call ludefvel_nolin
       if (myrank .eq. 0 .and. itimer .eq. 1) then
          call second(tend)
          t_ludefall = t_ludefall + tend - tstart
       end if
    end if
#endif

    ! copy field data to time-advance vectors
    if(myrank.eq.0 .and. iprint.ge.1) print *,"Importing time advance vectors.."
    call import_time_advance_vectors

    ! advance time
    if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing time..."
    if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
    if(isplitstep.ge.1) then
       call step_split(calc_matrices)
    else
       call step_unsplit(calc_matrices)
    end if
    if(myrank.eq.0 .and. itimer.eq.1) then 
       call second(tend)
       if(iprint.ge.1) print *, "Time spent in *_step: ", tend-tstart
    end if

    if(nplanes.le.1) exit
    ! check if time step should be repeated

#ifdef USE3D
    maxiter = 0
    do icount=1,maxnumofsolves
       maxiter = max(maxiter,int(kspits(icount)))
    enddo
    if(myrank.eq.0) write(*,'(A,2i5)') 'ksp_max, maxiter =', ksp_max, maxiter
    if((maxiter .lt. ksp_max) .or. dtkecrit .le. 0) exit

    dt = dt/2.
    if(myrank.eq.0) write(*,'(A,e12.4)') 'Time step reduced by 2:  dt =', dt
#else
    exit
#endif
  enddo


  time = time + dt
  dtold = dt
  if(ntime.gt.1 .and. linear.eq.0) call variable_timestep

  call pellet_advance

!  call runaway_advance

  ! copy time advance vectors to field data
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Exporting time advance vectors.."

! if(eqsubtract.eq.0) call subtract_axi    !DEBUG
  call export_time_advance_vectors
  !if (irunaway .eq. 2) then
  !   call smooth_runaway
  !endif

!#ifdef USEPARTICLES
!  if(kinetic.eq.0) then
!     call set_parallel_velocity
!  endif
!#endif

#ifdef USEPARTICLES
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(kinetic.eq.1) call particle_step(dt*t0_norm)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_particle = t_particle + tend - tstart
  endif
#endif

  ! Calculate all quantities derived from basic fields
  call find_lcfs()
  call derived_quantities(1)

  if(ipellet_abl.gt.0) call pellet_shrink

  ! Advect impurity charge states
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call kprad_advect(dt)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_kprad = t_kprad + tend - tstart
  endif

  ! Conserve toroidal flux
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     tflux = tflux + gbound*area
  endif


  first_time = .false.
  meshAdapted = 0
end subroutine onestep

!======================================================================
! import_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! construct vectors for time-advance
!======================================================================
subroutine import_time_advance_vectors
  use basic
  implicit none

  select case(isplitstep)
  case(0)
     call import_time_advance_vectors_unsplit
  case(1:2)
     call import_time_advance_vectors_split
  end select

end subroutine import_time_advance_vectors


!======================================================================
! export_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! extract values from time-advance vectors
!======================================================================
subroutine export_time_advance_vectors
  use basic
  implicit none

  select case(isplitstep)
  case(0)
     call export_time_advance_vectors_unsplit
  case(1:2)
     call export_time_advance_vectors_split
  end select
end subroutine export_time_advance_vectors



!=============================
! scaleback
! ~~~~~~~~~
! rescale eigenfunction
!=============================
subroutine scaleback

  use basic
  use arrays
  use diagnostics
  use particles

  implicit none

  vectype, parameter :: scalefac = 1.e-10

  if(ekin.lt.max_ke .or. max_ke.eq.0) return
  if(myrank.eq.0) write(*,*) " =>solution scaled back at time", time

  call mult(field_vec, scalefac)
  if(i3d.eq.1) call mult(bf_field(1), scalefac)
  if(i3d.eq.1) call mult(bfp_field(1), scalefac)

#ifdef USEPARTICLES
  if(kinetic.eq.1) then
     call particle_scaleback(scalefac)
     ! epar=epar*scalefac**2
     ! if (hostrank==0) then
     !    do ipart=ipart_begin,ipart_end
     !       pdata(ipart)%wt=pdata(ipart)%wt*scalefac
     !       pdata(ipart)%wt2=pdata(ipart)%wt2*scalefac**2
     !    end do
     !    do ielm=ielm_min,ielm_max
     !       elfieldcoefs(ielm)%psiv1=elfieldcoefs(ielm)%psiv1*scalefac
     !       elfieldcoefs(ielm)%Bzv1=elfieldcoefs(ielm)%Bzv1*scalefac
     !       elfieldcoefs(ielm)%Bfv=elfieldcoefs(ielm)%Bfv*scalefac
     !    end do
     ! endif
     ! call update_particle_pressure
  end if
#endif
end subroutine scaleback


subroutine variable_timestep

  use basic
  use arrays
  use diagnostics

  implicit none
  include 'mpif.h'
  integer :: ierr

#ifdef USE3D
  integer :: maxiter, icount
#endif

if(dtkecrit.gt.0) then
  if(myrank.eq.0) then
!
! decrease timestep based on kinetic energy or ksp iterations
! but limit change to fraction dtfrac and bound by dtmin and dtmax
!
!   
       if(ekin.gt.dtkecrit) then
         dt = dtold/(1. + dtfrac)
       else

#ifdef USE3D
         maxiter = 0
         do icount=1,maxnumofsolves
           maxiter = max(maxiter,int(kspits(icount)))
         enddo

         if(maxiter .gt. ksp_warn) dt = dtold/(1. + dtfrac)
         if(maxiter .lt. ksp_min)  dt = dtold*(1. + dtfrac)
         if(iprint.ge.1) write(*,'("maxiter, ksp_warn, ksp_min",3i5)') maxiter, ksp_warn, ksp_min
#endif

         dt = max(dt,dtmin)
         dt = min(dt,dtmax)

         if(iprint.ge.1) write(*,'("dtold,dt,dtkecrit,ekin",1p4e12.4)') dtold,dt,dtkecrit,ekin
       endif



  endif ! on myrank.eq.0
  call MPI_bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

endif   ! on dtkecrit

end subroutine variable_timestep


end module time_step

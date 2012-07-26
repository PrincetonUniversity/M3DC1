module time_step
  use time_step_split
  use time_step_unsplit
  implicit none

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

  implicit none

  integer :: calc_matrices, ivel_def
  logical, save :: first_time = .true.

  real :: tstart, tend

  ! Determine whether matrices should be re-calculated
  if(first_time &
       .or. (linear.eq.0 .and. mod(ntime,nskip).eq.0) &
       .or. (integrator.eq.1 .and. ntime.eq.1)) then
     calc_matrices = 1
  else
     calc_matrices = 0
  endif

  ! calculate matrices for time advance
  if(calc_matrices.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining matrices"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! in linear case, eliminate second-order terms from matrix
     ivel_def = 1
     if(istatic.eq.1) ivel_def = 0
     call ludefall(ivel_def, idens, ipres, ipressplit, 1-iestatic)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Done defining matrices."
  endif


  ! copy field data to time-advance vectors
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Importing time advance vectors.."
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

  time = time + dt
  if(ntime.gt.1 .and. linear.eq.0) call variable_timestep

  ! copy time advance vectors to field data
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Exporting time advance vectors.."
  call export_time_advance_vectors


  ! Calculate all quantities derived from basic fields
  call derived_quantities(1)


  ! Conserve toroidal flux
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     tflux = tflux + gbound*area
  endif


  first_time = .false.

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

  implicit none

  vectype, parameter :: scalefac = 1.e-10

  if(ekin.lt.max_ke .or. max_ke.eq.0) return
  if(myrank.eq.0) write(*,*) " =>solution scaled back at time", time

  call mult(field_vec, scalefac)
  if(i3d.eq.1) call mult(bf_field(1), scalefac)
  
end subroutine scaleback


subroutine variable_timestep

  use basic
  use arrays
  use diagnostics

  implicit none
  include 'mpif.h'
  integer :: ierr
!
! increase or decrease timestep based on kinetic energy and gamma_gr,
! but limit change to fraction dtfrac and bound by dtmin and dtmax
!
  dtold = dt
  if(dtkecrit.eq.0 .or. dtgamma.eq.0) return
  if(myrank.eq.0) then
!
!
    if(ekin.lt.dtkecrit) then
       if(gamma_gr.gt.0) then
          if(dt .gt.dtgamma/abs(gamma_gr)) then
            dt = dtold/(1. + dtfrac)
          else
            dt = dtold*(1. + dtfrac)
          endif
       else
            dt = dtold*(1. + dtfrac)
       endif
    else
            dt = dtold/(1. + dtfrac)
    endif
    dt = max(dt,dtmin)
    dt = min(dt,dtmax)
    if(iprint.ge.1) write(*,1001) dtold,dt,gamma_gr,ekin
  endif
  call MPI_bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
 1001 format("dtold,dt,gamma_gr,ekin",1p4e12.4)
end subroutine variable_timestep


end module time_step

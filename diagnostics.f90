module diagnostics

  implicit none

  ! scalar diagnostics
  real :: tflux, area, totcur, totden

  ! timing diagnostics
  real :: t_ludefall, t_sources, t_smoother, t_aux, t_onestep
  real :: t_solve_v, t_solve_n, t_solve_p, t_solve_b, t_mvm
  real :: t_output_cgm, t_output_hdf5, t_output_reset
  real :: t_gs, t_gs_magaxis, t_gs_fundef, t_gs_solve, t_gs_init



contains

  ! ======================================================================
  ! reset_timings
  ! 
  ! resents timing vaiables to 0
  ! ======================================================================
  subroutine reset_timings
    implicit none
    t_ludefall = 0.
    t_sources = 0.
    t_smoother = 0.
    t_aux = 0.
    t_solve_v = 0.
    t_solve_n = 0.
    t_solve_p = 0.
    t_solve_b = 0.
    t_output_cgm = 0.
    t_output_hdf5 = 0.
    t_output_reset = 0.
    t_mvm = 0.
    t_onestep = 0.
  end subroutine reset_timings


  ! ======================================================================
  ! distribute_timings
  ! 
  ! distributes timing vaiables to all processors
  ! ======================================================================
  subroutine distribute_timings

    implicit none

    include 'mpif.h'
    integer :: ier
    real, dimension(13) :: vin, vout

    vin(1) =  t_ludefall
    vin(2) =  t_sources
    vin(3) =  t_smoother
    vin(4) =  t_aux
    vin(5) =  t_solve_v
    vin(6) =  t_solve_n
    vin(7) =  t_solve_p
    vin(8) =  t_solve_b
    vin(9) =  t_output_cgm
    vin(10) = t_output_hdf5
    vin(11) = t_output_reset
    vin(12) = t_mvm
    vin(13) = t_onestep
    call MPI_ALLREDUCE(vin, vout, 13, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, ier)
    t_ludefall      = vout(1)
    t_sources       = vout(2)
    t_smoother      = vout(3)
    t_aux           = vout(4)
    t_solve_v       = vout(5)
    t_solve_n       = vout(6)
    t_solve_p       = vout(7)
    t_solve_b       = vout(8)
    t_output_cgm    = vout(9)
    t_output_hdf5   = vout(10)
    t_output_reset  = vout(11)
    t_mvm           = vout(12)
    t_onestep       = vout(13)
    
  end subroutine distribute_timings


  ! ======================================================================
  ! total_flux
  ! 
  ! calculates the total area, toroidal flux, toroidal current, and 
  ! electron number (indegrated density) in the computational domain
  ! ======================================================================
  subroutine total_flux

    use basic
    use nintegrate_mod

    implicit none
    
    include 'mpif.h'

    integer :: numelms, itri, ier, def_fields
    integer, parameter :: num_scalars = 4
    real, dimension(num_scalars) :: valsin, valsout

    area = 0.
    tflux = 0.
    totcur = 0.
    totden = 0.

    call numfac(numelms)

    def_fields = FIELD_PSI
    if(numvar.ge.2) def_fields = def_fields + FIELD_I
    if(idens.eq.1) def_fields = def_fields + FIELD_N

    ! calculate the toroidal current, area, toroidal flux
    ! and electron number
    do itri=1,numelms

       call define_fields_79(itri, def_fields)

       if(ijacobian.eq.1) weight_79 = weight_79*ri_79

       area = area + int0(weight_79,79)
       totcur = totcur - int2(ri_79,pst79(:,OP_GS),weight_79,79)
       if(numvar.ge.2) tflux = tflux + int2(ri_79,bzt79(:,OP_1),weight_79,79)
       if(idens.eq.1) totden = totden + int1(nt79(:,OP_1),weight_79,79)

    enddo                     ! loop over itri

    if(idens.eq.0) totden = area

    if(maxrank .gt. 1) then
       valsin(1) = area
       valsin(2) = tflux
       valsin(3) = totcur
       valsin(4) = totden
       call MPI_ALLREDUCE(valsin, valsout, num_scalars, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ier)
       area = valsout(1)
       tflux = valsout(2)
       totcur = valsout(3)
       totden = valsout(4)
    endif

    if(myrank.eq.0 .and. iprint.ge.1) then 
       print *, "Scalars:"
       print *, "  Area = ", area
       print *, "  Toroidal current = ", totcur
       print *, "  Total electrons = ", totden
    endif

  end subroutine total_flux


  ! ======================================================================
  ! reconnected_flux
  !
  ! calculates the reconnected flux assuming an x-point at the center of 
  ! the computational domain, and an o-point at the center of the right
  ! boundary.
  ! ======================================================================
  real function reconnected_flux()
    
    use basic
    use arrays
  
    implicit none

    real :: alx, alz
    real, dimension(2) :: temp, temp2
    integer :: itri

    call getboundingboxsize(alx,alz)

    itri = 0
    call evaluate(alx   ,alz/2.,temp(1),temp2(1),phi,1,numvar,itri)
    itri = 0
    call evaluate(alx/2.,alz/2.,temp(2),temp2(2),phi,1,numvar,itri)

    reconnected_flux = 0.5*(temp(2)-temp(1))

    return
  end function reconnected_flux

! ======================================================================
! second
!
! returns the current time in seconds
! ======================================================================
  subroutine second(tcpu)
    implicit none
!    real*4 etime
!    external etime
!    dimension tarray(2)
!    tcpu = etime(tarray) 
    real :: tcpu
!    integer*8 nticks, tickspersec
    integer :: nticks, tickspersec
    intrinsic system_clock
    call system_clock(count=nticks, count_rate=tickspersec)
    tcpu = 0.
    if(tickspersec .gt. 0) tcpu = nticks/(1.0*tickspersec)
    return
  end subroutine second


end module diagnostics


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
  call mpi_allreduce(chierror_local, chierror, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, ier)

end subroutine calc_chi_error

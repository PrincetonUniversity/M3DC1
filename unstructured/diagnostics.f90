module diagnostics

  implicit none

  real :: tflux, area, totcur

contains

  
  subroutine total_flux

    use basic
    use nintegrate_mod

    implicit none
    
    include 'mpif.h'

    integer :: numelms, itri, ier, def_fields
    real, dimension(3) :: valsin, valsout

    area = 0.
    tflux = 0.
    totcur = 0.

    call numfac(numelms)

    def_fields = FIELD_J
    if(numvar.ge.2) def_fields = def_fields + FIELD_I

    ! calculate the total perturbed current and area, and toroidal flux
    do itri=1,numelms

       call define_fields_79(itri, def_fields)

       area = area + int1(ri_79,weight_79,79)
       totcur = totcur - int2(ri2_79,jt79,weight_79,79)
       if(numvar.ge.2) tflux = tflux + int2(ri2_79,bzt79,weight_79,79)

    enddo                     ! loop over itri

    if(maxrank .gt. 1) then
       valsin(1) = area
       valsin(2) = tflux
       valsin(3) = totcur
       call MPI_ALLREDUCE(valsin, valsout, 3, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ier)
       area = valsout(1)
       tflux = valsout(2)
       totcur = valsout(3)
    endif

    if(myrank.eq.0 .and. iprint.ge.1) print *, "Area = ", area

  end subroutine total_flux

  
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

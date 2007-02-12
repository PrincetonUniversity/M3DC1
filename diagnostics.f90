module diagnostics

  implicit none

contains

  real function reconnected_flux()
    
    use basic
    use arrays
  
    implicit none

    include "mpif.h"

    real :: alx, alz
    real, dimension(2) :: temp, temp2
    integer :: ier, itri

    call getboundingboxsize(alx,alz)

    itri = 0
    call evaluate(alx   ,alz/2.,temp(1),temp2(1),phi,1,numvar,itri)
    itri = 0
    call evaluate(alx/2.,alz/2.,temp(2),temp2(2),phi,1,numvar,itri)

    print *, "Collecting data...", temp(1), temp(2)

    if(maxrank.eq.1) then
       temp2 = temp
    else
       call MPI_ALLREDUCE(temp, temp2, 2, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ier)
    endif

    reconnected_flux = 0.5*(temp2(2)-temp2(1))

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

module diagnostics

  implicit none

  real :: pflux, tflux, area, totcur

contains

  
  subroutine total_flux

    use basic
    use p_data
    use t_data
    use arrays

    implicit none
    
    include 'mpif.h'

    integer :: numelms, itri, j, jone, j1, j2, ier
    real :: d2term(18), fintl(-6:maxi,-6:maxi)
    real, dimension(4) :: valsin, valsout

    area = 0.
    pflux = 0.
    tflux = 0.
    totcur = 0.

    call numfac(numelms)

    ! calculate the total perturbed current and area, and toroidal flux
    do itri=1,numelms
       call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
       call calcd2term(itri, d2term, fintl)
       do j=1,18
          jone = isval1(itri,j)
          
          j1 = isvaln(itri,j)
          pflux = pflux + d2term(j)*(phi(j1) + phi0(j1))
          totcur = totcur + d2term(j)*jphi(jone)

          if(numvar.ge.2) then
             j2 = j1 + 6
             tflux = tflux + d2term(j)*(phi(j1) + phi0(j1))
          endif
       enddo
       do j=1,13,6
          area = area + d2term(j)
       enddo
    enddo                     ! loop over itri

    if(maxrank .gt. 1) then
       valsin(1) = area
       valsin(2) = pflux
       valsin(3) = tflux
       valsin(4) = totcur
       call MPI_ALLREDUCE(valsin, valsout, 4, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ier)
       area = valsout(1)
       pflux = valsout(2)
       tflux = valsout(3)
       totcur = valsout(4)
    endif

  end subroutine total_flux

  
  real function reconnected_flux()
    
    use basic
    use arrays
  
    implicit none

    real :: alx, alz
    real, dimension(2) :: temp, temp2
    integer :: ier, itri

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

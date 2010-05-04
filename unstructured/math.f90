!============================================================================
! math.f90
! ~~~~~~~~
!
! contains common mathematical functions
!============================================================================

module math

contains

  !==============
  ! factorial
  ! ~~~~~~~~~
  ! returns n!
  !==============
  integer function factorial(n)
    
    implicit none
    
    integer, intent(in) :: n
    integer :: i
    integer :: ans
    
    if (n.le.1) then
       factorial = 1
    else
       ans = 1
       do i=1,n
          ans = ans*i
       enddo
       factorial = ans
    end if
    return
  end function factorial
  
  !==============
  ! sech
  ! ~~~~
  ! returns sech(x)
  !================
  real function sech(x)
    implicit none
    real, intent(in) :: x
    
    sech = 1./cosh(x)
    return
  end function sech
  
end module math

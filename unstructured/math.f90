!============================================================================
! math.f90
! ~~~~~~~~
!
! contains common mathematical functions
!============================================================================

module math

  real, parameter :: pi    = 3.14159265358979323846
  real, parameter :: twopi = 2.*pi
  real, parameter :: amu0  = pi*4.e-7

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

  !====================================
  ! binomial
  ! ~~~~~~~~
  ! returns n choose m = n!/(m!(n-m)!)
  !====================================
  real function binomial(n, m)
    implicit none

    integer, intent(in) :: n, m

    binomial = factorial(n)/(factorial(m)*factorial(n-m))
  end function binomial


!============================================================
! rq_tri_int 
! ~~~~~~~~~~
! calculates the double integral over a triangle
! with base a + b and height c  [see ref 2, figure 1]
!
! rq_tri_int = Int{si**m*eta**n}d(si)d(eta)
!============================================================
real function rq_tri_int(m,n,a,b,c)
  implicit none

  integer, intent(in) :: m,n
  real, intent(in) :: a, b, c
  real :: anum, denom
  
  anum = c**(n+1)*(a**(m+1)-(-b)**(m+1))*factorial(m)*factorial(n)
  denom = factorial(m+n+2)
  rq_tri_int = anum/denom
end function rq_tri_int

! ===============================================================
! cubic_roots
! ~~~~~~~~~~~
! calculate the roots of a 3rd degree polynomial
! coef(4)*x**3 + coef(3)*x**2 + coef(2)*x + coef(1) = 0
!================================================================
subroutine cubic_roots(coef, root, error)

  implicit none
  real, intent(out), dimension(3) :: root, error
  real, intent(inout), dimension(4) :: coef
  real qqs, rrs, ths, fun, dfun
  integer r, ll

  ! normalize the coefficients to the 3rd degree coefficient
  coef(1) = coef(1)/coef(4)
  coef(2) = coef(2)/coef(4)
  coef(3) = coef(3)/coef(4)
  coef(4) = 1.0
       
  ! solve using method in Numerical Recipes
  qqs = (coef(3)**2 - 3.*coef(2))/9.
  rrs = (2.*coef(3)**3 - 9.*coef(3)*coef(2) + 27.*coef(1))/54.
  ths = acos(rrs/sqrt(qqs**3))

  root(1) = -2.*sqrt(qqs)*cos(ths/3.) - coef(3)/3.
  root(2) = -2.*sqrt(qqs)*cos((ths+twopi)/3.) - coef(3)/3.
  root(3) = -2.*sqrt(qqs)*cos((ths-twopi)/3.) - coef(3)/3.
  
  ! refine with Newton's method
  do r=1,3
     do ll=1,3
        fun  = root(r)**3 + coef(3)*root(r)**2 + coef(2)*root(r)    &
             + coef(1)
        dfun = 3.*root(r)**2 + 2.*coef(3)*root(r) + coef(2)
        root(r) = root(r) - fun/dfun
     enddo
     error(r) =  root(r)**3 + coef(3)*root(r)**2 + coef(2)*root(r)  &
          + coef(1)
  enddo
      
end subroutine cubic_roots

!======================================================================
! bessel functions
! ~~~~~~~~~~~~~~~~
! Courtesy of Gnu Scientific Library
!======================================================================
real function bessel_I(n, x)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x

  call gsl_bessel_i(n,x,bessel_I)
end function bessel_I

  
end module math

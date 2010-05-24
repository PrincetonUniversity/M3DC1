!=====================================================
! bicubic_interpolation_coeffs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! calculates bicubic polynomial coefficients a of
! x, an array of dimension m x n,
! about index (i,j)
!=====================================================
subroutine bicubic_interpolation_coeffs(x,m,n,i,j,a)
  implicit none

  integer, intent(in) :: m, n, i, j
  real, intent(in), dimension(m,n) :: x
  real, intent(out), dimension(4,4) :: a

  a = 0.

  if(i.lt.1 .or. i.gt.m) return
  if(j.lt.1 .or. j.gt.n) return

  a(1,1) = x(i,j)
  a(1,2) = (-2.*x(i,j-1)-3.*x(i,j)+6.*x(i,j+1)-x(i,j+2))/6.
  a(1,3) = (    x(i,j-1)-2.*x(i,j)+   x(i,j+1)         )/2.
  a(1,4) = (   -x(i,j-1)+3.*x(i,j)-3.*x(i,j+1)+x(i,j+2))/6.
  
  a(2,1) = (-2.*x(i-1,j)- 3.*x(i,j)+6.*x(i+1,j)-x(i+2,j))/6.
  a(2,2) = ( 4.*x(i-1,j-1) +6.*x(i-1,j)-12.*x(i-1,j+1)+2.*x(i-1,j+2) &
            +6.*x(i  ,j-1) +9.*x(i  ,j)-18.*x(i  ,j+1)+3.*x(i  ,j+2) &
           -12.*x(i+1,j-1)-18.*x(i+1,j)+36.*x(i+1,j+1)-6.*x(i+1,j+2) &
            +2.*x(i+2,j-1) +3.*x(i+2,j) -6.*x(i+2,j+1)   +x(i+2,j+2))/36.
  a(2,3) = (-2.*x(i-1,j-1) +4.*x(i-1,j) -2.*x(i-1,j+1) &
            -3.*x(i  ,j-1) +6.*x(i  ,j) -3.*x(i  ,j+1) &
            +6.*x(i+1,j-1)-12.*x(i+1,j) +6.*x(i+1,j+1) &
               -x(i+2,j-1)+ 2.*x(i+2,j)    -x(i+2,j+1))/12.
  a(2,4) = (2.*x(i-1,j-1) -6.*x(i-1,j) +6.*x(i-1,j+1)-2.*x(i-1,j+2) &
           +3.*x(i  ,j-1) -9.*x(i  ,j) +9.*x(i  ,j+1)-3.*x(i  ,j+2) &
           -6.*x(i+1,j-1)+18.*x(i+1,j)-18.*x(i+1,j+1)+6.*x(i+1,j+2) &
              +x(i+2,j-1) -3.*x(i+2,j) +3.*x(i+2,j+1)   -x(i+2,j+2))/36.
  
  a(3,1) = (x(i-1,j)-2.*x(i,j)+x(i+1,j))/2.
  a(3,2) = (-2.*x(i-1,j-1)-3.*x(i-1,j) +6.*x(i-1,j+1)   -x(i-1,j+2) &
            +4.*x(i  ,j-1)+6.*x(i  ,j)-12.*x(i  ,j+1)+2.*x(i  ,j+2) &
            -2.*x(i+1,j-1)-3.*x(i+1,j) +6.*x(i+1,j+1)   -x(i+1,j+2))/12.
  a(3,3) = (   x(i-1,j-1)-2.*x(i-1,j)   +x(i-1,j+1) &
           -2.*x(i  ,j-1)+4.*x(i  ,j)-2.*x(i  ,j+1) &
              +x(i+1,j-1)-2.*x(i+1,j)   +x(i+1,j+1))/4.
  a(3,4) = (  -x(i-1,j-1)+3.*x(i-1,j)-3.*x(i-1,j+1)+   x(i-1,j+2) &
           +2.*x(i  ,j-1)-6.*x(i  ,j)+6.*x(i  ,j+1)-2.*x(i  ,j+2) &
              -x(i+1,j-1)+3.*x(i+1,j)-3.*x(i+1,j+1)+   x(i+1,j+2))/12.

  a(4,1) = (-x(i-1,j)+3.*x(i,j)-3.*x(i+1,j)+x(i+2,j))/6.
  a(4,2) = (2.*x(i-1,j-1)+3.*x(i-1,j) -6.*x(i-1,j+1)   +x(i-1,j+2) &
           -6.*x(i  ,j-1)-9.*x(i  ,j)+18.*x(i  ,j+1)-3.*x(i  ,j+2) &
           +6.*x(i+1,j-1)+9.*x(i+1,j)-18.*x(i+1,j+1)+3.*x(i+1,j+2) &
           -2.*x(i+2,j-1)-3.*x(i+2,j) +6.*x(i+2,j+1)   -x(i+2,j+2))/36.
  a(4,3) = (  -x(i-1,j-1)+2.*x(i-1,j)   -x(i-1,j+1) &
           +3.*x(i  ,j-1)-6.*x(i  ,j)+3.*x(i  ,j+1) &
           -3.*x(i+1,j-1)+6.*x(i+1,j)-3.*x(i+1,j+1) &
              +x(i+2,j-1)-2.*x(i+2,j)   +x(i+2,j+1))/12.
  a(4,4) = (   x(i-1,j-1)-3.*x(i-1,j)+3.*x(i-1,j+1)   -x(i-1,j+2) &
           -3.*x(i  ,j-1)+9.*x(i  ,j)-9.*x(i  ,j+1)+3.*x(i  ,j+2) &
           +3.*x(i+1,j-1)-9.*x(i+1,j)+9.*x(i+1,j+1)-3.*x(i+1,j+2) &
              -x(i+2,j-1)+3.*x(i+2,j)-3.*x(i+2,j+1)   +x(i+2,j+2))/36.  
end subroutine bicubic_interpolation_coeffs

!=====================================================
! cubic_interpolation_coeffs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! calculates cubic polynomial coefficients a of
! x, an array of dimension m,
! about index i
!=====================================================
subroutine cubic_interpolation_coeffs(x,m,i,a)

  implicit none

  integer, intent(in) :: m, i
  real, intent(in), dimension(m) :: x
  real, intent(out), dimension(4) :: a
  integer :: ihermite
  real  :: xpi,xpip

  ihermite = 1
  select case(ihermite)
  case(0)
     a(1) = x(i)
     if(i.eq.1) then
        a(2) = (-3.*x(i) + 4.*x(i+1) - x(i+2))/2.
        a(3) = (    x(i) - 2.*x(i+1) + x(i+2))/2.
        a(4) = 0.
     else if(i.eq.m-1) then
        a(2) = (-2.*x(i-1) - 3.*x(i) + 5.*x(i+1))/6.
        a(3) = (    x(i-1) - 2.*x(i)    + x(i+1))/2.
        a(4) = (   -x(i-1) + 3.*x(i) - 2.*x(i+1))/6.
     else if(i.eq.m) then
        a(2) = (-x(i-1) + x(i))/3.
        a(3) = ( x(i-1) - x(i))/2.
        a(4) = (-x(i-1) + x(i))/6.
!!$        a(2) = (-x(i-1) + x(i))/4.
!!$        a(3) = ( x(i-1) - x(i))/2.
!!$        a(4) = (-x(i-1) + x(i))/4.
     else
        a(2) = (-2.*x(i-1) - 3.*x(i) + 6.*x(i+1) - x(i+2))/6.
        a(3) = (    x(i-1) - 2.*x(i)    + x(i+1)         )/2.
        a(4) = (   -x(i-1) + 3.*x(i) - 3.*x(i+1) + x(i+2))/6.
     end if
  case(1)
    a(1) = x(i)
    if(i.eq.1) then
      xpip= .5*(x(i+2)-x(i))
      a(2) = -2*x(i) +2*x(i+1) - xpip     
      a(3) =    x(i) -  x(i+1) + xpip
      a(4) = 0.
    else if(i.eq.m-1) then
      xpi = .5*(x(i+1)-x(i-1))
      xpip=    (x(i+1)-x(i))
      a(2) = xpi
      a(3) = -.5*x(i+1) + 2.0*x(i) - 2.5*x(i-1) +     x(i-2)
      a(4) =  .5*x(i+1) - 1.5*x(i) + 1.5*x(i-1) - 0.5*x(i-2)
    else if(i.eq.m) then
      xpi = .5*(x(i+1)-x(i-1))
      xpip=     x(i) - x(i-1)
      a(2) = xpip
      a(3) = -.5*x(i) + 2.0*x(i-1) - 2.5*x(i-2) +     x(i-3)
      a(4) =  .5*x(i) - 1.5*x(i-1) + 1.5*x(i-2) - 0.5*x(i-3)
    else
      xpi = .5*(x(i+1)-x(i-1))
      xpip= .5*(x(i+2)-x(i))
      a(2) = xpi
      a(3) = -3*x(i) +3*x(i+1) - 2*xpi - xpip
      a(4) =  2*x(i) -2*x(i+1) +   xpi + xpip
    endif
  endselect
end subroutine cubic_interpolation_coeffs

!=====================================================
! cubic_interpolation
! ~~~~~~~~~~~~~~~~~~~
!
! interpolates function f(p) at point p0
! where f and p are arrays of length m
!=====================================================
subroutine cubic_interpolation(m, p, p0, f, f0)
  implicit none

  integer, intent(in) :: m
  real, intent(in), dimension(m) :: f, p
  real, intent(in) :: p0
  real, intent(out) :: f0

  real, dimension(4) :: a
  integer :: i, iguess
  real :: dp

  if(p0.lt.p(1)) then 
     f0 = f(1)
     return
  end if
  if(p0.gt.p(m)) then
     f0 = f(m)
     return
  end if

  iguess = (m-1)*(p0-p(1))/(p(m)-p(1)) + 1

  if(p(iguess).lt.p0) goto 90
  if(p(iguess).gt.p0) goto 80
  i = iguess
  goto 100

  ! search backward
80 continue
  do i=iguess, 1, -1
     if(p(i).le.p0) goto 100
  end do
  goto 100

90 continue 
  ! search forward
  do i=iguess, m-1
     if(p(i+1).gt.p0) goto 100
  end do

100 continue
  if(i .lt. 1) i = 1
  if(i .gt. m) i = m

  call cubic_interpolation_coeffs(f,m,i,a)

  if(i.ge.m) then
     f0 = a(1)
  else
     dp = (p0-p(i))/(p(i+1)-p(i))
     f0 = a(1) + a(2)*dp + a(3)*dp**2 + a(4)*dp**3
  endif
  
end subroutine cubic_interpolation

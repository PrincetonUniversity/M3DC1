!=====================================================
! bicubic_interpolation_coeffs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! calculates bicubic polynomial coefficients a of
! x, an array of dimension m x n,
! about index (i,j)
!=====================================================
subroutine bicubic_interpolation_coeffs(x,m,n,i0,j0,a,ierr)
  implicit none

  integer, intent(in) :: m, n
  real, intent(in) :: i0, j0
  real, intent(in), dimension(m,n) :: x
  real, intent(out), dimension(4,4) :: a
  integer, intent(out) :: ierr

  integer :: i, j

  i = i0
  j = j0

  ierr = 0
  if(i.lt.2) ierr = 1
  if(j.lt.2) ierr = 1
  if(i.gt.m-2) ierr = 1
  if(j.gt.n-2) ierr = 1
  if(ierr.ne.0) then
!     print *, 'Out of bounds ', i, j, m, n
     return
  end if

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
! bicubic_interpolation
! ~~~~~~~~~~~~~~~~~~~~~
!
! interpolates function f^(d)(x,y) at point (x0, y0)
!=====================================================
subroutine bicubic_interpolation(m, n, x, y, x0, y0, f, f0, dfdx, dfdy, ierr)
  implicit none

  integer, intent(in) :: m, n
  real, intent(in) :: f(m,n), x(m), y(n), x0, y0
  real, intent(out) :: f0, dfdx, dfdy
  integer, intent(out) :: ierr

  real, dimension(4,4) :: a
  real :: i0, j0
  integer :: i, j
  real :: di, dj, dxdi, dydj

  i0 = (x0 - x(1))/(x(m) - x(1))
  j0 = (y0 - y(1))/(y(n) - y(1))

  call bicubic_interpolation_coeffs(f,m,n,i0,j0,a,ierr)
  if(ierr.ne.0) return

  i = i0
  j = j0
  di = i0 - i
  dj = j0 - j

  dxdi = x(i+1)-x(i)
  dydj = y(i+1)-y(i)

  f0 =       a(1,1) + di*(a(2,1) + di*(a(3,1) + di*a(4,1)))    &
       + dj*(a(1,2) + di*(a(2,2) + di*(a(3,2) + di*a(4,2)))    &
       + dj*(a(1,3) + di*(a(2,3) + di*(a(3,3) + di*a(4,3)))    &
       + dj*(a(1,4) + di*(a(2,4) + di*(a(3,4) + di*a(4,4))))))
  dfdx =    (a(2,1) + di*(2.*a(3,1) + di*3.*a(4,1))            &
       + dj*(a(2,2) + di*(2.*a(3,2) + di*3.*a(4,2))            &
       + dj*(a(2,3) + di*(2.*a(3,3) + di*3.*a(4,3))            &
       + dj*(a(2,4) + di*(2.*a(3,4) + di*3.*a(4,4))))))/dxdi
  dfdy =    (a(1,2) + di*(a(2,2) + di*(a(3,2) + di*a(4,2)))    &
       + dj*(a(1,3) + di*(a(2,3) + di*(a(3,3) + di*a(4,3)))    &
       + dj*(a(1,4) + di*(a(2,4) + di*(a(3,4) + di*a(4,4))))))/dydj

  ierr = 0
end subroutine bicubic_interpolation


!=====================================================
! cubic_interpolation_coeffs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! calculates cubic polynomial coefficients a of
! x, an array of dimension m,
! about index i
!=====================================================
subroutine cubic_interpolation_coeffs(x,m,i0,a)

  implicit none

  integer, intent(in) :: m, i0
  real, intent(in), dimension(m) :: x
  real, intent(out), dimension(4) :: a
  integer :: i

  i = i0
  if(i.lt.1) i=1
  if(i.gt.m) i=m

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

!  if(p(m).eq.p(1)) print *, 'ERROR 1!!!', p(1), p(m) 
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
!     if(p(i+1).eq.p(i)) print *, 'ERROR 2!!!', p(i+1), p(i) 
     dp = (p0-p(i))/(p(i+1)-p(i))
     f0 = a(1) + a(2)*dp + a(3)*dp**2 + a(4)*dp**3
  endif

end subroutine cubic_interpolation

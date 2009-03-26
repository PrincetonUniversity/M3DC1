!=====================================================
! bicubic_interpolation
! ~~~~~~~~~~~~~~~~~~~~~
!
! calculates bicubic polynomial coefficients a of
! x, an array of dimension m x n,
! about index (i,j)
!=====================================================
subroutine bicubic_interpolation(x,m,n,i,j,a)
  implicit none

  integer, intent(in) :: m, n, i, j
  real, intent(in), dimension(m,n) :: x
  real, intent(out), dimension(4,4) :: a

  if(i.lt.1 .or. i.gt.m) return
  if(j.lt.1 .or. j.gt.n) return

  a(1,1) = x(i,j)
  a(1,2) = (-2.*x(i,j-1)-3.*x(i,j)+6.*x(i,j+1)-x(i,j+2))/6.
  a(1,3) = (    x(i,j-1)-2.*x(i,j)+   x(i,j+1))/2.
  a(1,4) = (   -x(i,j-1)+3.*x(i,j)-3.*x(i,j+1)+x(i,j+2))/6.
  
  a(2,1) = (-2.*x(i-1,j)-3.*x(i,j)+6.*x(i+1,j)-x(i+2,j))/6.
  a(2,2) = (4.*x(i-1,j-1)+6.*x(i-1,j)-12.*x(i-1,j+1) &
           +2.*x(i-1,j+2)+6.*x(i,j-1)+9.*x(i,j) &
           -18.*x(i,j+1)+3.*x(i,j+2)-12.*x(i+1,j-1) &
           -18.*x(i+1,j)+36.*x(i+1,j+1)-6.*x(i+1,j+2) &
           +2.*x(i+2,j-1)+3.*x(i+2,j)-6.*x(i+2,j+1) &
           +x(i+2,j+2))/36.
  a(2,3) = (-2.*x(i-1,j-1)+4.*x(i-1,j)-2.*x(i-1,j+1) &
            -3.*x(i,j-1)+6.*x(i,j)-3.*x(i,j+1) &
            +6*x(i+1,j-1)-12.*x(i+1,j)+6.*x(i+1,j+1) &
            -x(i+2,j-1)+2.*x(i+2,j)-x(i+2,j+1))/12.
  a(2,4) = (2.*x(i-1,j-1)-6.*x(i-1,j)+6.*x(i-1,j+1) &
           -2.*x(i-1,j+2)+3.*x(i,j-1)-9.*x(i,j) &
           +9.*x(i,j+1)-3.*x(i,j+2)-6.*x(i+1,j-1) &
           +18.*x(i+1,j)-18.*x(i+1,j+1)+6.*x(i+1,j+2) &
           +x(i+2,j-1)-3.*x(i+2,j)+3.*x(i+2,j+1) &
           -x(i+2,j+2))/36.
  
  a(3,1) = (x(i-1,j)-2.*x(i,j)+x(i+1,j))/2.
  a(3,2) = (-2.*x(i-1,j-1)-3.*x(i-1,j)+6.*x(i-1,j+1) &
            -x(i-1,j+2)+4.*x(i,j-1)+6.*x(i,j) &
            -12.*x(i,j+1)+2.*x(i,j+2)-2.*x(i+1,j-1) &
            -3.*x(i+1,j)+6.*x(i+1,j+1)-x(i+1,j+2))/12.
  a(3,3) = (x(i-1,j-1)-2.*x(i-1,j)+x(i-1,j+1) &
           -2.*x(i,j-1)+4.*x(i,j)-2.*x(i,j+1) &
           +x(i+1,j-1)-2.*x(i+1,j)+x(i+1,j+1))/4.
  a(3,4) = (-x(i-1,j-1)+3.*x(i-1,j)-3.*x(i-1,j+1) &
            +x(i-1,j+2)+2.*x(i,j-1)-6.*x(i,j) &
            +6.*x(i,j+1)-2.*x(i,j+2)-x(i+1,j-1) &
            +3.*x(i+1,j)-3.*x(i+1,j+1)+x(i+1,j+2))/12.

  a(4,1) = (-x(i-1,j)+3.*x(i,j)-3.*x(i+1,j)+x(i+2,j))/6.
  a(4,2) = (2.*x(i-1,j-1)+3.*x(i-1,j)-6.*x(i-1,j+1) &
           +x(i-1,j+2)-6.*x(i,j-1)-9.*x(i,j) &
           +18.*x(i,j+1)-3.*x(i,j+2)+6.*x(i+1,j-1) &
           +9.*x(i+1,j)-18.*x(i+1,j+1)+3.*x(i+1,j+2) &
           -2.*x(i+2,j-1)-3.*x(i+2,j)+6.*x(i+2,j+1) &
           -x(i+2,j+2))/36.
  a(4,3) = (-x(i-1,j-1)+2.*x(i-1,j)-x(i-1,j+1) &
            +3.*x(i,j-1)-6.*x(i,j)+3.*x(i,j+1) &
            -3.*x(i+1,j-1)+6.*x(i+1,j)-3.*x(i+1,j+1) &
            +x(i+2,j-1)-2.*x(i+2,j)+x(i+2,j+1))/12.
  a(4,4) = (x(i-1,j-1)-3.*x(i-1,j)+3.*x(i-1,j+1) &
           -x(i-1,j+2)-3.*x(i,j-1)+9.*x(i,j) &
           -9.*x(i,j+1)+3.*x(i,j+2)+3.*x(i+1,j-1) &
           -9.*x(i+1,j)+9.*x(i+1,j+1)-3.*x(i+1,j+2) &
           -x(i+2,j-1)+3.*x(i+2,j)-3.*x(i+2,j+1) &
           +x(i+2,j+2))/36.  
end subroutine bicubic_interpolation

!=====================================================
! cubic_interpolation
! ~~~~~~~~~~~~~~~~~~~
!
! calculates cubic polynomial coefficients a of
! x, an array of dimension m x n,
! about index i
!=====================================================
subroutine cubic_interpolation(x,m,i,a)

  implicit none

  integer, intent(in) :: m, i
  real, intent(in), dimension(m) :: x
  real, dimension(4) :: a

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
     a(2) = (-2.*x(i-1) + 2.*x(i))/6.
     a(3) = (    x(i-1) -    x(i))/2.
     a(4) = (   -x(i-1) +    x(i))/6.
  else
     a(2) = (-2.*x(i-1) - 3.*x(i) + 6.*x(i+1) - x(i+2))/6.
     a(3) = (    x(i-1) - 2.*x(i)    + x(i+1)         )/2.
     a(4) = (   -x(i-1) + 3.*x(i) - 3.*x(i+1) + x(i+2))/6.
  end if
end subroutine cubic_interpolation

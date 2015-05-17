module spline

  implicit none

  type spline1d
     real, allocatable :: x(:)
     real, allocatable :: y(:)
     integer :: n
  end type spline1d

contains

  subroutine copy_spline(s1, s2)
    implicit none

    type(spline1d) :: s1
    type(spline1d), intent(in) :: s2

    if(s1%n .ne. s2%n) then
       if(allocated(s1%x)) deallocate(s1%x)
       if(allocated(s1%y)) deallocate(s1%y)
       s1%n = s2%n
       allocate(s1%x(s1%n), s1%y(s1%n))
    end if

    s1%x = s2%x
    s1%y = s2%y
  end subroutine copy_spline

  subroutine create_spline(s, n, x, y)
    implicit none
    
    type(spline1d) :: s
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: x, y

    integer :: i

    if(allocated(s%x)) deallocate(s%x)
    if(allocated(s%y)) deallocate(s%y)

    allocate(s%x(n), s%y(n))
    s%n = n

    ! put x in increasing order
    if(x(1) .lt. x(n)) then 
       s%x = x
       s%y = y
    else
       print *, 'REVERSING', x(1), x(n)
       s%x = x(n:1:-1)
       s%y = y(n:1:-1)
    end if

    ! check for monotonicity
    do i=1, n-1
       if(s%x(i) .ge. s%x(i+1)) then
          print *, 'Warning: X not monotonic!', i, s%x(i), s%x(i+1)
       end if
    end do
  end subroutine create_spline

  subroutine destroy_spline(s)
    implicit none

    type(spline1d) :: s
    if(s%n.eq.0) return
    if(allocated(s%x)) deallocate(s%x)
    if(allocated(s%y)) deallocate(s%y)
    s%n = 0
  end subroutine destroy_spline

  subroutine which_domain(s, x, i, dx)
    implicit none
    type(spline1d), intent(in) :: s
    real, intent(in) :: x
    integer, intent(out) :: i
    real, intent(out) :: dx

    ! check if point is outside domain
    if(x.lt.s%x(1)) then
       i = 0
       dx = x - s%x(1)
       return
    end if
    if(x.ge.s%x(s%n)) then
       i = s%n
       dx = x - s%x(s%n-1)
       return
    end if

    ! find appropriate domain
    do i=1, s%n-1
       if(x.lt.s%x(i+1)) exit
    end do
!    if(i.lt.1 .or. i.ge.s%n) print *, 'ERROR!!!', i, s%n, x, s%x(1), s%x(s%n)
    dx = x - s%x(i)
  end subroutine which_domain

  subroutine get_hermite_coeffs(s, ii, a, extrapolate)
    implicit none

    type(spline1d), intent(in) :: s
    integer, intent(in) :: ii
    real, dimension(4), intent(out) :: a
    integer, intent(in), optional :: extrapolate

    integer :: ex, i
    real :: dx, dydx, dydx0, dydx1
    
    if(present(extrapolate)) then
       ex = extrapolate
    else
       ex = 0
    end if

    ! make sure that i is inside spline domain
    if(ii.lt.1) then
       i = 1
    else if(ii.ge.s%n) then
       i = s%n-1
    else
       i = ii
    end if

    ! calculate polynomial coefficients using hermite interpolation
    dx = s%x(i+1) - s%x(i)
    dydx = (s%y(i+1) - s%y(i))/dx
    if(i.eq.1) then
       dydx0 = dydx
       dydx1 = (s%y(i+2) - s%y(i  ))/(s%x(i+2) - s%x(i  ))
    else if(i.eq.s%n-1) then
       dydx0 = (s%y(i+1) - s%y(i-1))/(s%x(i+1) - s%x(i-1))
       dydx1 = dydx
    else
       dydx0 = (s%y(i+1) - s%y(i-1))/(s%x(i+1) - s%x(i-1))
       dydx1 = (s%y(i+2) - s%y(i  ))/(s%x(i+2) - s%x(i  ))
    end if
    a(1) = s%y(i)
    a(2) = dydx0
    a(3) = (3.*dydx - (2.*dydx0+dydx1))/dx
    a(4) = (-2.*dydx + (dydx0+dydx1))/dx**2

    ! if outside of domain, restrict coeffs based on extrapolation param
    if(ii.lt.1) then
       if(ex.lt.0) a(1) = 0.
       if(ex.lt.1) a(2) = 0.
       if(ex.lt.2) a(3) = 0.
       if(ex.lt.3) a(4) = 0.
    else if(ii.ge.s%n) then
       a(1) = s%y(s%n)
       if(ex.lt.0) a(1) = 0.
       if(ex.lt.1) a(2) = 0.
       if(ex.lt.2) a(3) = 0.
       if(ex.lt.3) a(4) = 0.
    end if

  end subroutine get_hermite_coeffs

  subroutine evaluate_spline(s, x, y, yp, ypp, yppp, iout, extrapolate)
    implicit none

    type(spline1d), intent(in) :: s
    real, intent(in) :: x
    real, intent(out) :: y
    real, intent(out), optional :: yp, ypp, yppp
    integer, intent(in), optional :: extrapolate
    integer, intent(out), optional :: iout

    real, dimension(4) :: a
    real :: dx
    integer :: i

    if(present(iout)) then 
       if(x .gt. s%x(s%n)) then
          iout = 1
       elseif(x .lt. s%x(1)) then
          iout = -1
       else
          iout = 0
       end if
    end if

    call which_domain(s, x, i, dx)
    call get_hermite_coeffs(s, i, a, extrapolate)
    
    ! calculate values and derivatives
    y =             a(1) + a(2)*dx +    a(3)*dx**2 +    a(4)*dx**3
    if(present(yp)) yp   = a(2)    + 2.*a(3)*dx    + 3.*a(4)*dx**2
    if(present(ypp)) ypp =           2.*a(3)       + 6.*a(4)*dx   
    if(present(yppp)) yppp =                         6.*a(4)      
  end subroutine evaluate_spline

  subroutine integrate_spline(s, x0, x1, y, extrapolate)
    implicit none

    type(spline1d), intent(in) ::s
    real, intent(in) :: x0, x1
    real, intent(out) :: y
    integer, intent(in), optional :: extrapolate

    integer :: i, i0, i1
    real :: fac, dx0, dx1, dx
    real, dimension(4) :: a
  
    call which_domain(s, x0, i0, dx0)
    call which_domain(s, x1, i1, dx1)

    ! if limits are backwards, reverse them
    if(i1.lt.i0) then 
       i = i0
       i0 = i1
       i1 = i
       fac = dx0
       dx0 = dx1
       dx1 = fac
       fac = -1.
    else
       fac = 1.
    end if

    y = 0.
    do i=i0, i1
       call get_hermite_coeffs(s, i, a, extrapolate)
       ! lower limit
       if(i.eq.i0) then 
          dx = dx
       else
          dx = 0.
       endif
       y = y - a(1)*dx - a(2)*dx**2/2. - a(3)*dx**3/3. - a(4)*dx**4/4.

       ! upper limit
       if(i.eq.i1) then
          dx = dx1
       else
          if(i.lt.1 .or. i.ge.s%n) then
             print *, 'ERROR in spline integration!', i, s%n
             return
          end if
          dx = s%x(i+1) - s%x(i)
       endif
       y = y + a(1)*dx + a(2)*dx**2/2. + a(3)*dx**3/3. + a(4)*dx**4/4.
    end do

    y = y*fac
  end subroutine integrate_spline


end module spline

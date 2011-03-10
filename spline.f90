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

    allocate(s%x(n), s%y(n))
    s%x = x
    s%y = y
    s%n = n
  end subroutine create_spline

  subroutine destroy_spline(s)
    implicit none

    type(spline1d) :: s
    if(s%n.eq.0) return
    deallocate(s%x, s%y)
    s%n = 0
  end subroutine destroy_spline

  subroutine evaluate_spline(s, x, y, yp, ypp, yppp)
    implicit none

    type(spline1d), intent(in) :: s
    real, intent(in) :: x
    real, intent(out) :: y
    real, intent(out), optional :: yp, ypp, yppp

    real, dimension(4) :: a
    real :: dx, di, didx
    integer :: i

    ! If x is ouside domain, extrapolate using edge value only
    i = 0
    if(s%x(1) .lt. s%x(s%n)) then
       if(x.le.s%x(1)) then
          i = 1
       else if(x.ge.s%x(s%n)) then
          i = s%n
       end if
    else
       if(x.ge.s%x(1)) then
          i = 1
       else if(x.le.s%x(s%n)) then
          i = s%n
       end if
    end if
    if(i.ne.0) then
       y = s%y(i)
       if(present(yp)) yp = 0.
       if(present(ypp)) ypp = 0.
       if(present(yppp)) yppp = 0.
       return
    end if

    ! otherwise, find appropriate domain and interpolate
    if(s%x(1) .lt. s%x(s%n)) then
       do i=1, s%n-1
          if(x.lt.s%x(i+1)) exit
       end do
    else
       do i=1, s%n-1
          if(x.gt.s%x(i+1)) exit
       end do
    end if
    if(i.lt.1 .or. i.gt.s%n) print *, 'ERROR!!!'

    a(1) = s%y(i)
    if(i.eq.1) then
       a(2) = (-3.*s%y(i) + 4.*s%y(i+1) - s%y(i+2))/2.
       a(3) = (    s%y(i) - 2.*s%y(i+1) + s%y(i+2))/2.
       a(4) = 0.
    else if(i.eq.s%n-1) then
       a(2) = (-2.*s%y(i-1) - 3.*s%y(i) + 5.*s%y(i+1))/6.
       a(3) = (    s%y(i-1) - 2.*s%y(i)    + s%y(i+1))/2.
       a(4) = (   -s%y(i-1) + 3.*s%y(i) - 2.*s%y(i+1))/6.
    else if(i.eq.s%n) then
       a(2) = (-s%y(i-1) + s%y(i))/3.
       a(3) = ( s%y(i-1) - s%y(i))/2.
       a(4) = (-s%y(i-1) + s%y(i))/6.
    else
       a(2) = (-2.*s%y(i-1) - 3.*s%y(i) + 6.*s%y(i+1) - s%y(i+2))/6.
       a(3) = (    s%y(i-1) - 2.*s%y(i)    + s%y(i+1)           )/2.
       a(4) = (   -s%y(i-1) + 3.*s%y(i) - 3.*s%y(i+1) + s%y(i+2))/6.
    end if

    dx = x - s%x(i)
    didx = 1./(s%x(2) - s%x(1))
    di = dx*didx

    y =              a(1) + a(2)*di +    a(3)*di**2 +    a(4)*di**3
    if(present(yp)) yp   = (a(2)    + 2.*a(3)*di    + 3.*a(4)*di**2)*didx
    if(present(ypp)) ypp = (          2.*a(3)       + 6.*a(4)*di   )*didx**2
    if(present(yppp)) yppp =                          6.*a(4)       *didx**3

  end subroutine evaluate_spline

end module spline

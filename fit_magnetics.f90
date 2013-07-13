module fit_magnetic
  implicit none

  ! Maximum radial and vertical modes to keep
  integer, private :: n_max, m_max

  ! Total number of degrees of freedom
  integer, private :: ntot

  ! Number of measurements
  integer, private :: num_points

  ! (R, phi, Z) location of each measurement
  real, private, allocatable :: point(:,:)

  ! Direction of each measurement
  real, private, allocatable :: norm(:,:)

  ! Value of each measurement
  real, private, allocatable :: value(:)

  ! Largest vertical scale-length allowed
  real, private :: Lz

  complex, private, allocatable :: mat(:,:)

contains

  subroutine fit_read_data(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer, parameter :: ifile = 15
    integer :: n, i

    open(unit=ifile, file=filename, action='read', status='old')    

    ! count lines
    n = 0
    do
       read(ifile, *, end=100, err=100)
       n = n + 1
    end do
    
100 continue
    print *, 'Read ', n, ' lines'

    num_points = n

    ! allocate data
    allocate(point(3,num_points), norm(3,num_points), value(num_points))

    rewind(ifile)

    ! read data
    do i=1, n
       
    end do

    close(ifile)
  end subroutine fit_read_data

  subroutine fit_finalize
    implicit none

    deallocate(point, norm, value, mat)
  end subroutine fit_finalize

  subroutine mat_element(n, m, pt, norm, val)
    use math
    implicit none
    
    integer, intent(in) :: n, m
    real, intent(in), dimension(3) :: pt, norm
    complex, intent(out) :: val

    real :: k

    k = 2.*pi*m/Lz
    val = exp(cmplx(0.,1.)*(n*pt(2) + k*pt(3))) * &
         (0.5*norm(1)*k*(Bessel_I(n-1, k*pt(1)) + Bessel_I(n+1, k*pt(1))) &
         + (norm(2)*n/pt(1) + norm(3)*k)*cmplx(0.,1.)*Bessel_I(n, k*pt(1)))
  end subroutine mat_element

  subroutine mat_row(pt, norm, val)
    implicit none

    real, intent(in), dimension(3) :: pt, norm
    complex, intent(out), dimension(ntot) :: val

    integer :: n, m, i

    i = 1
    do n=-n_max, n_max
       do m=-m_max, m_max
          call mat_element(n, m, pt, norm, val(i))
          i = i + 1 
       end do
    end do
  end subroutine mat_row

  subroutine mat_create(pt, norm)
    implicit none
    
    real, intent(in), dimension(3, num_points) :: pt, norm
    integer :: i

    ntot = (2.*n_max + 1)*(2.*m_max + 1)
    allocate(mat(ntot,num_points)) 

    do i=1, num_points
       call mat_row(pt(:,i), norm(:,i), mat(:,i))
    end do
    
  end subroutine mat_create


end module fit_magnetic

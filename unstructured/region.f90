module region
  implicit none

  type plane_type
     integer :: n
     real, dimension(:), allocatable :: x   ! coordinates of points projected onto plane with unit normal 
     real, dimension(:), allocatable :: y   ! coordinates of points projected onto plane with unit normal 
     real, dimension(3) :: norm             ! unit normal direction (R,phi,Z)
  end type plane_type

  ! Axisymmetric regions have nplanes = 1
  type region_type
     integer :: nplanes
     type(plane_type), dimension(:), allocatable :: plane
  end type region_type
  
  contains

    subroutine create_plane(p, nr)
      implicit none
      
      type(plane_type), intent(inout) :: p
      integer, intent(in) :: nr

      p%n = nr
      allocate(p%x(nr))
      allocate(p%y(nr))

      if(nr.eq.1) then 
         p%norm(1) = 0.
         p%norm(2) = 1.
         p%norm(3) = 0.
      end if
    end subroutine create_plane

    subroutine create_region(r, np)
      implicit none

      type(region_type), intent(inout) :: r
      integer, intent(in) :: np

      r%nplanes = np
      allocate(r%plane(np))
    end subroutine create_region

    subroutine destroy_region(r)
      implicit none
      
      type(region_type), intent(inout) :: r
      integer :: i

      do i=1, r%nplanes
         call destroy_plane(r%plane(i))
      end do
      deallocate(r%plane)
      r%nplanes = 0
    end subroutine destroy_region

    subroutine destroy_plane(p)
      implicit none
      type(plane_type), intent(inout) :: p

      if(allocated(p%x)) deallocate(p%x)
      if(allocated(p%y)) deallocate(p%y)
      p%n = 0
    end subroutine destroy_plane

    subroutine create_region_from_file(r, filename, ierr)
      implicit none

      include 'mpif.h'

      type(region_type), intent(inout) :: r
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ierr

      integer, parameter :: ifile = 45
      integer :: i, n, rank
      real :: x, y

      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

      ierr = 0

      if(rank.eq.0) then
         open(unit=ifile, file=filename, action='READ', status='OLD')

         read(ifile, *) n

      end if

      call mpi_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call create_region(r, 1)
      call create_plane(r%plane(1), n+1)

      if(rank.eq.0) then
         do i=1, n
            read(ifile, *)  x, y
            r%plane(1)%x(i) = x
            r%plane(1)%y(i) = y
         end do

         r%plane(1)%x(n+1) = r%plane(1)%x(1)
         r%plane(1)%y(n+1) = r%plane(1)%y(1)

         close(ifile)
      end if

      call mpi_bcast(r%plane(1)%x, n+1, MPI_DOUBLE_PRECISION, &
           0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(r%plane(1)%y, n+1, MPI_DOUBLE_PRECISION, &
           0, MPI_COMM_WORLD, ierr)

    end subroutine create_region_from_file

    pure logical function point_in_region(reg, r, phi, z)
      implicit none

      type(region_type), intent(in) :: reg
      real, intent(in) :: r, phi, z

      if(reg%nplanes.eq.1) then
         point_in_region = point_in_plane(reg%plane(1), r, z)
      else
         point_in_region = .false.
      endif
    end function point_in_region


    pure logical function point_in_plane(p, x, y)
      implicit none

      type(plane_type), intent(in) :: p
      real, intent(in) :: x, y

      integer :: i
      real :: m, b
      
      point_in_plane = .false.

      ! For each line segment, determine if line running from
      ! (x,y) to (infinity, y) passes through the segment.
      ! If an odd number of segments are crossed, the point is inside the plane
      do i=1, p%n-1
         if(y.gt.p%y(i) .and. y.gt.p%y(i+1)) cycle

         if(y.lt.p%y(i) .and. y.lt.p%y(i+1)) cycle

         if(p%y(i) .eq. p%y(i+1)) then
            point_in_plane = .not.point_in_plane
            cycle
         end if

         m = (p%x(i+1) - p%x(i)) / (p%y(i+1) - p%y(i))
         b = p%x(i) - m*p%y(i)
         if(x .le. m*y + b) point_in_plane = .not.point_in_plane
      end do

    end function point_in_plane
end module region

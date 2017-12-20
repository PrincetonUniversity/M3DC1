module read_ascii

contains

  !======================================================================
  ! read_ascii_column
  ! ~~~~~~~~~~~~~~~~~
  ! reads a column of an ascii file into x
  ! n = number of rows to read on input, rows read on output
  !     if n<=0 on input, number of rows will be determined automatically
  ! skip = number of header rows to skip (default = 0)
  ! xrow = column to read (default = 1)
  !======================================================================
  subroutine read_ascii_column(filename, x, n, skip, icol, read_until)
    implicit none

    include 'mpif.h'

    character(len=*), intent(in) :: filename
    real, allocatable :: x(:)
    integer, intent(inout) :: n
    integer, optional :: skip, icol
    character(len=*), intent(in), optional :: read_until

    integer :: i, ix, myrank, ierr
    integer, parameter :: ifile = 112
    real, allocatable :: val(:)
    character(len=256) :: buff

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    ! Determine which columns to read
    if(present(icol)) then
       ix = icol
    else
       ix = 1
    end if

    if(myrank.eq.0) then
       print *, 'Reading column ', ix, ' of ', filename, '...'
    end if

    ! if we don't know how many rows there are, count them
    if(n.le.0) then
       if(myrank.eq.0) then
          open(unit=ifile, file=filename, status='old', action='read', err=101)

          ! read until
          if(present(read_until)) then
             do 
                read(ifile, '(A)', err=100) buff
                if(buff(1:len(read_until)) .eq. read_until) exit
             end do
          end if

          ! skip header rows
          if(present(skip)) then
             do i=1, skip
                read(ifile, *)
             end do
          end if

          ! count lines until EOF
          n = 0
          do
             read(ifile, *, err=100, end=100)
             n = n+1
          end do

100       close(ifile)
          goto 102
101       n = 0
102       print *, ' Read ', n, ' lines'
       end if
       
       call MPI_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end if

    ! If no valid lines found, quit
    if(n.le.0) return

    ! Allocate array to read row data
    allocate(val(ix))

    ! Allocate profile arrays
    allocate(x(n))

    ! Read data
    if(myrank.eq.0) then
       open(unit=ifile, file=filename, status='old', action='read')

       ! read until
       if(present(read_until)) then
          do 
             read(ifile, '(A)') buff
             if(buff(1:len(read_until)) .eq. read_until) exit
          end do
       end if
       
       ! skip header rows
       if(present(skip)) then
          do i=1, skip
             read(ifile, *)
          end do
       end if

       ! read rows
       do i=1, n
          read(ifile, *) val
          x(i) = val(ix)
       end do
       
       close(ifile)
    end if
    
    deallocate(val)

    ! Share data with other processes
    call MPI_bcast(x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_ascii_column
end module read_ascii

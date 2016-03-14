module read_ascii
  implicit none

  interface read_ascii_column
     module procedure read_ascii_column_real
     module procedure read_ascii_column_int
  end interface

contains

  subroutine get_token(line, m, tok, ierr)
    implicit none
    
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: tok
    integer, intent(in) :: m
    integer, intent(out) :: ierr

    integer :: i, n, i0, nc

    i = 1
    n = len_trim(line)

    nc = 0
    ierr = 1
    
    do while(i.le.n)
       do while(line(i:i).eq.' ') 
          i = i + 1
          if(i.gt.n) return
       enddo
       i0 = i
       nc = nc + 1
       do while(line(i:i).ne.' ')
          i = i + 1
          if(i.gt.n) exit
       enddo
       if(m.eq.nc) then
          tok = trim(line(i0:i))
          ierr = 0
          return
       end if
    enddo
  end subroutine get_token
  
  !======================================================================
  ! read_ascii_column
  ! ~~~~~~~~~~~~~~~~~
  ! reads a column of an ascii file into x
  ! n = number of rows to read on input, rows read on output
  !     if n<=0 on input, number of rows will be determined automatically
  ! skip = number of header rows to skip (default = 0)
  ! xrow = column to read (default = 1)
  !======================================================================
  subroutine read_ascii_column_real(filename, x, n, skip, icol, read_until)
    implicit none

    include 'mpif.h'

    character(len=*), intent(in) :: filename
    real, allocatable :: x(:)
    integer, intent(inout) :: n
    integer, optional :: skip, icol
    character(len=*), intent(in), optional :: read_until

    integer :: i, ix, myrank, ierr
    integer, parameter :: ifile = 112
    character(len=256) :: buff, tok

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

    if(n.le.0) then
       ! count the number of lines and columns
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
    end if

    call MPI_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! If no valid lines found, quit
    if(n.le.0) return

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
          read(ifile, '(A)') buff
          call get_token(buff, ix, tok, ierr)
          if(ierr.eq.0) then
             read(tok, *) x(i)
          else
             x(i) = 0.
          end if
       end do
       
       close(ifile)
    end if

    ! Share data with other processes
    call MPI_bcast(x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_ascii_column_real

  subroutine read_ascii_column_int(filename, x, n, skip, icol, read_until)
    implicit none

    character(len=*), intent(in) :: filename
    integer, allocatable :: x(:)
    integer, intent(inout) :: n
    integer, optional :: skip, icol
    character(len=*), intent(in), optional :: read_until

    real, allocatable :: y(:)

    call read_ascii_column_real(filename, y, n, skip, icol, read_until)

    allocate(x(n))
    x = nint(y)
    deallocate(y)
  end subroutine read_ascii_column_int
end module read_ascii

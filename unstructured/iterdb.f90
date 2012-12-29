module iterdb

  integer :: idb_ishot
  integer :: idb_nj
  integer :: idb_nion
  integer :: idb_nprim
  integer :: idb_nimp
  integer :: idb_nneu
  integer :: idb_ibion
  integer :: idb_nplasbdry

  character(len=1), allocatable :: idb_namep(:)
  character(len=1), allocatable :: idb_namei(:)
  character(len=1), allocatable :: idb_namen(:)

  integer, parameter :: idb_nreal = 20
  real, dimension(idb_nreal) :: idb_real_data

  integer, parameter :: idb_nprof = 73
  real, allocatable :: idb_profile_data(:,:)

  real, allocatable :: idb_psi(:), idb_Te(:), idb_ne(:), idb_omega(:)

  real, allocatable :: idb_plasbdry(:,:)

  logical, private :: initialized = .false.

contains

  subroutine load_iterdb(filename, ierr)

    implicit none

    include 'mpif.h'

    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr
    
    integer :: ifile = 12
    integer :: size, rank, i

    call MPI_comm_size(MPI_COMM_WORLD,size,ierr)
    call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)

    if(rank.eq.0) then
       open(unit=ifile,file=filename,status='old',action='read')

       call read_next_integer(ifile, idb_ishot, ierr)
       call read_next_integer(ifile, idb_nj, ierr)
       call read_next_integer(ifile, idb_nion, ierr)
       call read_next_integer(ifile, idb_nprim, ierr)
       call read_next_integer(ifile, idb_nimp, ierr)
       call read_next_integer(ifile, idb_nneu, ierr)
       call read_next_integer(ifile, idb_ibion, ierr)
    endif

    call mpi_bcast(idb_ishot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_nj,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_nion,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_nprim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_nimp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_nneu,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_ibion, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    allocate(idb_namep(idb_nprim), idb_namei(idb_nimp), idb_namen(idb_nneu))
    allocate(idb_psi(idb_nj),idb_ne(idb_nj),idb_Te(idb_nj),idb_omega(idb_nj))
!!$    allocate(idb_profile_data(idb_nj,idb_nprof))

    if(rank.eq.0) then
       call read_next_char(ifile, idb_nprim, idb_namep, ierr)
       call read_next_char(ifile, idb_nimp,  idb_namei, ierr)
       call read_next_char(ifile, idb_nneu,  idb_namen, ierr)
    end if

    call mpi_bcast(idb_namep, idb_nprim, &
         MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_namei, idb_nimp,  &
         MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_namen, idb_nneu,  &
         MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    if(rank.eq.0) then
       call read_next_real(ifile, idb_nj, idb_psi, ierr, &
            "*  psi on rho grid, volt*second/radian")
       call read_next_real(ifile, idb_nj, idb_Te, ierr, &
            "*  electron temperature, keV")
       call read_next_real(ifile, idb_nj, idb_ne, ierr, &
            "*  electron density, #/meter**3")
       call read_next_real(ifile, idb_nj, idb_omega, ierr, &
            "*  angular rotation speed profile, rad/sec")
    end if

    call mpi_bcast(idb_psi, idb_nj, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_Te, idb_nj, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_ne, idb_nj, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(idb_omega, idb_nj, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!!$    if(rank.eq.0) then
!!$       do i=1, 10
!!$          call read_next_real(ifile, 1, idb_real_data(i), ierr)
!!$       end do
!!$       call read_next_real(ifile, 5, idb_real_data(11:15), ierr)
!!$       do i=16, 20
!!$          call read_next_real(ifile, 1, idb_real_data(i), ierr)
!!$       end do
!!$    end if
!!$    call mpi_bcast(idb_real_data, idb_nreal, &
!!$         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!!$    
!!$    if(rank.eq.0) then
!!$       do i=1, 70
!!$          call read_next_real(ifile, idb_nj, idb_profile_data(:,i), ierr)
!!$       end do
!!$       call read_next_integer(ifile, idb_nplasbdry, ierr)
!!$    end if
!!$
!!$    call mpi_bcast(idb_nplasbdry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!!$
!!$    allocate(idb_plasbdry(idb_nplasbdry, 2))
!!$
!!$    if(rank.eq.0) then
!!$       call read_next_real(ifile, idb_nplasbdry, idb_plasbdry(:,1), ierr) 
!!$       call read_next_real(ifile, idb_nplasbdry, idb_plasbdry(:,2), ierr)
!!$       do i=71, idb_nprof
!!$          call read_next_real(ifile, idb_nj, idb_profile_data(:,i), ierr)
!!$       end do
!!$    end if
!!$    call mpi_bcast(idb_profile_data, idb_nprof*idb_nj, &
!!$         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if(rank.eq.0) then 
       close(ifile)
    end if   

    initialized = .true.
  end subroutine load_iterdb
  
  subroutine unload_iterdb
    implicit none

    if(.not.initialized) return

    deallocate(idb_namep, idb_namei, idb_namen)
    deallocate(idb_psi, idb_ne, idb_Te, idb_omega)
!!$    deallocate(idb_profile_data)
!!$    deallocate(idb_plasbdry)

    initialized = .false.
    
  end subroutine unload_iterdb

  subroutine read_next_integer(ifile, ival, ierr)
    implicit none
    
    integer, intent(in) :: ifile
    integer, intent(out) :: ival, ierr

    character(len=256) :: linebuf

    ierr = 0
    do
       read(ifile, '(A)', end=100) linebuf
       if(linebuf(1:1) .ne. '*') exit
    end do

    read(linebuf, *) ival
    return

100 ierr=1
  end subroutine read_next_integer

  subroutine read_next_real(ifile, n, val, ierr, name)
    implicit none
    
    integer, intent(in) :: ifile, n
    real, intent(out), dimension(n) :: val
    integer, intent(out) :: ierr
    character(len=*), optional :: name

    character(len=256) :: linebuf

    ierr = 0

    if(present(name)) then
       do 
          read(ifile, '(A)', end=100) linebuf
          if(trim(linebuf) .eq. trim(name)) exit
       end do
    end if
    do
       read(ifile, '(A)', end=100) linebuf
       if(linebuf(1:1) .ne. '*') exit
    end do

90  format (5e16.4)
    if(n.gt.5) then
       read(linebuf, 90) val(1:5)
       read(ifile, 90) val(6:n)
    else
       read(linebuf, 90) val
    end if
    return

100 ierr=1
  end subroutine read_next_real

  subroutine read_next_char(ifile, n, ch, ierr)
    implicit none
    
    integer, intent(in) :: ifile, n
    character(len=1), intent(out), dimension(n) :: ch
    integer, intent(out) :: ierr

    character(len=256) :: linebuf

    ierr = 0
    do
       read(ifile, '(A)', end=100) linebuf
       if(linebuf(1:1) .ne. '*') exit
    end do

    read(linebuf, *) ch
    return

100 ierr=1
  end subroutine read_next_char



end module iterdb

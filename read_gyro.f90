module gyro

  integer :: gyro_nexp
  real :: gyro_bt_exp
  real :: gyro_arho_exp

  integer, parameter :: gyro_nprofiles = 40
  real, allocatable :: gyro_profile(:,:)

  integer, parameter :: gyro_rho      =  1
  integer, parameter :: gyro_rmin     =  2
  integer, parameter :: gyro_rmaj     =  3
  integer, parameter :: gyro_q        =  4
  integer, parameter :: gyro_kappa    =  5
  integer, parameter :: gyro_delta    =  6
  integer, parameter :: gyro_te       =  7
  integer, parameter :: gyro_ne       =  8
  integer, parameter :: gyro_zeff     =  9
  integer, parameter :: gyro_omega0   = 10
  integer, parameter :: gyro_flow_mom = 11
  integer, parameter :: gyro_pow_e    = 12
  integer, parameter :: gyro_pow_i    = 13
  integer, parameter :: gyro_pow_ei   = 14
  integer, parameter :: gyro_zeta     = 15
  integer, parameter :: gyro_flow_beam= 16
  integer, parameter :: gyro_flow_wall= 17
  integer, parameter :: gyro_zmag     = 18
  integer, parameter :: gyro_ptot     = 19
  integer, parameter :: gyro_polflux  = 20
  integer, parameter :: gyro_ni(5)    = (/ 21, 22, 23, 24, 25 /)
  integer, parameter :: gyro_ti(5)    = (/ 26, 27, 28, 29, 30 /)
  integer, parameter :: gyro_vtor(5)  = (/ 31, 32, 33, 34, 35 /)
  integer, parameter :: gyro_vpol(5)  = (/ 36, 37, 38, 39, 40 /)

  logical, private :: initialized = .false.

contains

  subroutine load_gyro(ierr)

    implicit none

    include 'mpif.h'

    integer, intent(out) :: ierr

    integer :: size, rank, ier, eqpos, i, j
    integer, parameter :: ifile=24
    character(len=256) :: linebuf

    ierr = 0
    call MPI_comm_size(MPI_COMM_WORLD,size,ierr)
    call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)

    if(rank.eq.0) then
       print *, 'Loading gyro data'

       open(ifile, file='input.profiles', action='read', status='old', err=100)

       ! skip header
       do
          read(ifile, '(A)', end=100) linebuf
          if(linebuf(1:1).ne.'#') exit
       end do

       eqpos = scan(linebuf, '=')
       read(linebuf(eqpos+1:), *) gyro_nexp
       read(ifile, '(A)', end=100) linebuf
       eqpos = scan(linebuf, '=')
       read(linebuf(eqpos+1:), *) gyro_bt_exp
       read(ifile, '(A)', end=100) linebuf
       eqpos = scan(linebuf, '=')
       read(linebuf(eqpos+1:), *) gyro_arho_exp

       ! skip labels

       write(*,'(A,I6,2f14.5)') 'nexp, bt, arho = ', &
            gyro_nexp, gyro_bt_exp, gyro_arho_exp

       goto 101
100    ierr = 1
       close(ifile)
101    continue
    end if
    call mpi_bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
    if(ierr.ne.0) then
       print *, 'error ', rank
       if(rank.eq.0) print *, 'Error reading input.profiles'
       return
    endif

    call mpi_bcast(gyro_nexp,    1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gyro_bt_exp,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gyro_arho_exp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    allocate(gyro_profile(gyro_nexp, gyro_nprofiles))

    if(rank.eq.0) then
       do i=1, gyro_nprofiles/5
          read(ifile, *)
          read(ifile, *)
          do j=1, gyro_nexp
             read(ifile, '(e14.7, 4e16.7)') & 
                  gyro_profile(j, 5*(i-1)+1:5*(i-1)+5)
          end do
       end do

       close(ifile)
    end if

    call mpi_bcast(gyro_profile, gyro_nexp*gyro_nprofiles, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    initialized = .true.
  end subroutine load_gyro

  subroutine unload_gyro
    implicit none

    if(.not.initialized) return

    deallocate(gyro_profile)
  end subroutine unload_gyro

end module gyro

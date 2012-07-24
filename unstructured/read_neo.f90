module neo
  use gyro
  use spline

  integer, private :: ntheta, nr, nspecies, mtheta

  real, allocatable, private :: vpol_c(:,:,:), vpol_s(:,:,:)
  real, allocatable, private :: vtor_c(:,:,:), vtor_s(:,:,:)

  real, allocatable, private :: theta(:)
  real, allocatable, private :: r(:)

  real, private :: anorm
  real, private, allocatable :: vnorm(:)

  type(spline1d), private :: index_vs_psi

contains

  subroutine read_neo(ierr)

    implicit none

    integer, intent(out) :: ierr
    integer :: i

    ierr = 0

    call read_neo_grid(ierr)
    if(ierr.ne.0) return

    call read_neo_expnorm(ierr)
    if(ierr.ne.0) return

    mtheta = (ntheta-1)/2

    allocate(vpol_c(mtheta, nr, nspecies), vpol_s(mtheta, nr, nspecies))
    allocate(vtor_c(mtheta, nr, nspecies), vtor_s(mtheta, nr, nspecies))
    call read_neo_vel_fourier(ierr)
    if(ierr.ne.0) return

    ! convert r, vpol, and vtor into mks units
    r = r*anorm 
    do i=1, nr
       vpol_c(:,i,:) = vpol_c(:,i,:)*vnorm(i)
       vpol_s(:,i,:) = vpol_s(:,i,:)*vnorm(i)
       vtor_c(:,i,:) = vtor_c(:,i,:)*vnorm(i)
       vtor_s(:,i,:) = vtor_s(:,i,:)*vnorm(i)
    end do

    ! read gyro profiles to get psi vs. r
    call load_gyro(ierr)

    call neo_setup_spline()
  end subroutine read_neo

  subroutine unload_neo()
    implicit none

    if(allocated(vnorm)) deallocate(vnorm)

    if(allocated(theta)) deallocate(theta)
    if(allocated(r)) deallocate(r)

    if(allocated(vpol_c)) deallocate(vpol_c)
    if(allocated(vpol_s)) deallocate(vpol_s)
    if(allocated(vtor_c)) deallocate(vtor_c)
    if(allocated(vtor_s)) deallocate(vtor_s)
  end subroutine unload_neo

  subroutine read_neo_expnorm(ierr)
    implicit none

    include 'mpif.h'

    integer, intent(out) :: ierr
    integer, parameter :: ifile=22
    integer :: size, rank, ier, i
    real :: dum

    ierr = 0
    call MPI_comm_size(MPI_COMM_WORLD,size,ierr)
    call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)

    allocate(vnorm(nr))

    anorm = 0
    vnorm = 0
  
    if(rank.eq.0) then
       print *, 'Reading neo normalizations'

       open(ifile, file='out.neo.expnorm',action='read',status='old',err=100)
       
       do i=1, nr
          read(ifile, '(7E16.8)', err=100, end=100) &
               dum, anorm, dum, dum, dum, vnorm(i), dum
       end do

       goto 101
100    ierr = 1
101    close(ifile)
    end if

    call mpi_bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
    if(ierr.ne.0) then
       if(rank.eq.0) print *, 'Error reading out.neo.expnorm'
       return
    endif

    call mpi_bcast(anorm,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vnorm, nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_neo_expnorm

  subroutine read_neo_grid(ierr)
    implicit none

    include 'mpif.h'

    integer, intent(out) :: ierr
    integer, parameter :: ifile=22
    integer :: size, rank, ier, i

    ierr = 0
    call MPI_comm_size(MPI_COMM_WORLD,size,ierr)
    call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)

    nspecies = 0
    nr = 0
    ntheta = 0

    if(rank.eq.0) then
       print *, 'Reading neo grid'

       open(ifile, file='out.neo.grid', action='read', status='old', err=100)

       read(ifile, *, err=100, end=100) nspecies
       read(ifile, *, err=100, end=100) ! n_energy
       read(ifile, *, err=100, end=100) ! n_xi
       read(ifile, *, err=100, end=100) ntheta
       
       allocate(theta(ntheta))

       do i=1, ntheta
          read(ifile, *, err=100, end=100) theta(i)
       end do
         
       read(ifile, *, err=100, end=100) nr
       
       allocate(r(nr))

       do i=1, nr
          read(ifile, *, err=100, end=100) r(i)
       end do

       goto 101
100    ierr = 1
101    close(ifile)
    end if

    call mpi_bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
    if(ierr.ne.0) then
       if(rank.eq.0) print *, 'Error reading out.neo.grid'
       return
    endif

    call mpi_bcast(ntheta,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nr,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nspecies, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if(rank.ne.0) then
       allocate(theta(ntheta))
       allocate(r(nr))
    end if
    
    call mpi_bcast(theta, ntheta, MPI_DOUBLE_PRECISION, &
         0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(r, nr, MPI_DOUBLE_PRECISION, &
         0, MPI_COMM_WORLD, ierr)

    if(rank.eq.0) then
       write(*, '(A,3I5)') '  ntheta, nr, nspecies = ', ntheta, nr, nspecies
    end if
  end subroutine read_neo_grid


  subroutine read_neo_vel_fourier(ierr)
    implicit none

    include 'mpif.h'

    integer, intent(out) :: ierr
    integer, parameter :: ifile=23
    integer :: size, rank, ier
    integer :: n, i, j, pos
    real, allocatable :: buffer(:)
    
    n = nspecies*nr*mtheta

    ierr = 0
    call MPI_comm_size(MPI_COMM_WORLD,size,ierr)
    call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)

    if(rank.eq.0) then
       open(ifile, file='out.neo.vel_fourier', action='read', status='old', &
            err=100)

       allocate(buffer(6*mtheta*nspecies))

       do i=1, nr
          read(ifile, *, err=100, end=100) buffer
          do j=1, nspecies
             pos = 6*mtheta*(j-1)+1
             vpol_c(:,i,j) = buffer(pos+mtheta*2:pos+mtheta*3-1)
             vpol_s(:,i,j) = buffer(pos+mtheta*3:pos+mtheta*4-1)
             vtor_c(:,i,j) = buffer(pos+mtheta*4:pos+mtheta*5-1)
             vtor_s(:,i,j) = buffer(pos+mtheta*5:pos+mtheta*6-1)
          end do
       end do

       deallocate(buffer)

       close(ifile)
       goto 101
100    ierr = 1
101    continue
    end if
    call mpi_bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
    if(ierr.ne.0) then
       if(rank.eq.0) print *, 'Error reading out.neo.vel_fourier'
       return
    endif

    call mpi_bcast(vpol_c, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vpol_s, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vtor_c, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(vtor_s, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_neo_vel_fourier

  subroutine neo_setup_spline()
    implicit none

    type(spline1d) :: index_vs_r
    type(spline1d) :: r_vs_psi

    integer :: i
    real, allocatable :: ind(:), rtmp(:)
    real :: rx

    ! the +1 is for an extra point at the separatrix
    allocate(ind(nr+1), rtmp(nr+1))
    do i=1, nr+1
       ind(i) = i
    end do
    rtmp(1:nr) = r
    rtmp(nr+1) = (gyro_profile(gyro_nexp, gyro_rmin) + r(nr))/2.
    call create_spline(index_vs_r, nr+1, rtmp, ind)
    deallocate(ind, rtmp)

    call create_spline(r_vs_psi, gyro_nexp, &
         gyro_profile(:,gyro_polflux), gyro_profile(:,gyro_rmin))

    ! for each psi in gyro_polflux, calculate corresponding index
    allocate(ind(gyro_nexp))
    do i=1, gyro_nexp
       call evaluate_spline(r_vs_psi, gyro_profile(i, gyro_polflux), rx)
       call evaluate_spline(index_vs_r, rx, ind(i))
    end do

    call create_spline(index_vs_psi, gyro_nexp, &
         gyro_profile(:,gyro_polflux), ind)
    deallocate(ind)

    ! shift bottom index to give vpol = 0 at axis
!    index_vs_psi%x(1) = index_vs_psi%x(2)/2.

    call destroy_spline(index_vs_r)
    call destroy_spline(r_vs_psi)
  end subroutine neo_setup_spline


  subroutine neo_eval_vel(n, psi, theta, vpol, vtor, iout)
    implicit none

    integer, intent(in) :: n
    real, dimension(n), intent(in) :: psi, theta
    real, dimension(n), intent(out) :: vpol, vtor
    integer, dimension(n), intent(out) :: iout

    integer, parameter :: species = 1
    integer :: i, j, i1, i2
    real :: xind, di, co, sn
    real :: vp1, vp2, vt1, vt2

    iout = 0
    do i=1, n
       ! determine the appropriate index and offset for given psi value
       call evaluate_spline(index_vs_psi, psi(i), xind)
       i1 = xind
       if(i1.ge.nr+1 .or. i1.lt.0) then
          iout(i) = 1
          vtor(i) = 0.
          vpol(i) = 0.
          cycle
       end if
       if(i1.lt.1) i1 = 1
       di = xind - i1
       if(i1.gt.nr) i1 = nr
       if(di.gt.1.) di = 1.
       i2 = i1 + 1

       vp1 = 0.
       vp2 = 0.
       vt1 = 0.
       vt2 = 0.
       do j=1, mtheta
          co = cos(theta(i)*(j-1))
          sn = sin(theta(i)*(j-1))

          if(i1.gt.1) then
             vp1 = vp1 &
                  + vpol_c(j,i1,species)*co &
                  + vpol_s(j,i1,species)*sn
             vt1 = vt1 &
                  + vtor_c(j,i1,species)*co &
                  + vtor_s(j,i1,species)*sn
          else
             vt1 = vt1 &
                  + vtor_c(j,1,species)*co &
                  + vtor_s(j,1,species)*sn
          end if
             
          if(i2.le.nr) then
             vp2 = vp2 &
                  + vpol_c(j,i2,species)*co &
                  + vpol_s(j,i2,species)*sn
             vt2 = vt2 &
                  + vtor_c(j,i2,species)*co &
                  + vtor_s(j,i2,species)*sn
          end if
       end do

       ! do linear interpolation between surfaces
       vpol(i) = vp1*(1.-di) + vp2*di
       vtor(i) = vt1*(1.-di) + vt2*di
    end do

  end subroutine neo_eval_vel

end module neo

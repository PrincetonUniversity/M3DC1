module read_schaffer_field
  implicit none

  integer, private :: nr, nz, nphi

  real, private, allocatable :: aphi(:,:,:)
  real, private, allocatable :: bphi(:,:,:)
  real, private, allocatable :: r(:), phi(:), z(:)

contains

  subroutine load_schaffer_field(filename, ierr)
    implicit none

    include 'mpif.h'

    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr

    integer, parameter :: header_lines = 5
    integer, parameter :: ifile = 37
    integer :: i, j, k, ier
    integer :: rank
    real :: phi1, r1, z1, a, b
    real :: phi0, r0, z0
    character(len=10) :: dummy
    

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ier)

    ! First, let rank zero parse the file to determine the size of the data
    ierr = 0
    if(rank.eq.0) then 
       open(ifile, name=filename, status='old', action='read', err=200)

       do i=1, header_lines
          read(ifile,*) dummy
          print *, dummy
       end do

       ! determine number of R's, Z's, and Phi's
       nz = 1
       nr = 1
       nphi = 1
       read(ifile,'(3F20.11)') phi0, r0, z0 
       do
          read(ifile,'(3F20.11)',end=100) phi1, r1, z1 
          if(phi1 .ne. phi0) then
             phi0 = phi1
             nphi = nphi + 1
             nz = 1
             nr = 1
          end if
          if(z1 .ne. z0) then
             z0 = z1
             nz = nz + 1
             nr = 1
          end if
          nr = nr + 1
       end do

100    continue
       nr = nr-1
       nz = nz-1
    end if

    goto 300
200 continue
    ierr = 1
    print *, 'Error reading ', filename
300 continue

    ! Tell all processors whether rank 0 had an error
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    if(ierr.ne.0) return
       
    ! Send size data to all processors
    call MPI_Bcast(nr,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(nz,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(nphi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    write(*,'(A,3I5)') 'nr, nz, nphi = ', nr, nz, nphi

    ! Allocate space for data
    allocate(aphi(nphi,nr,nz),bphi(nphi,nr,nz))
    allocate(r(nr),phi(nphi),z(nz))

    ! Now read in the data
    if(rank.eq.0) then
       ! position read at start of file again
       rewind(ifile,err=2000)

       do i=1, header_lines 
          read(ifile,*) dummy
          print *, dummy
       end do

       do k=1, nphi
          do j=1, nz
             do i=1, nr
                read(ifile,'(5F20.11)',end=1000) phi1, r1, z1, a, b
                aphi(k,i,j) = a
                bphi(k,i,j) = b
                r(i) = r1
             end do
             z(j) = z1
          end do
          phi(k) = phi1
       end do

       goto 1100
1000   continue
       print *, "Shouldn't be here!!!", i, j, k
1100   close(ifile)
    end if

    goto 3000
2000 continue
    print *, 'Error rewinding ', filename
    ierr = 1
3000 continue

    ! Tell all processes whether rank 0 had an error
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    if(ierr.ne.0) return

    ! Send data to all processes
    call MPI_Bcast(r,   nr,   MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(z,   nz,   MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(phi, nphi, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(aphi, nr*nz*nphi, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bphi, nr*nz*nphi, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

    if(rank.eq.0) call calculate_external_field_ft(1)
  end subroutine load_schaffer_field



  subroutine unload_schaffer_field
    implicit none

    deallocate(aphi,bphi,r,phi,z)
  end subroutine unload_schaffer_field


  subroutine calculate_external_field_ft(ntor)
    implicit none

    integer, intent(in) :: ntor

    complex, dimension(nphi) :: out
    real, dimension(nphi) :: in

    in = aphi(:,20,20)

    call fftw(in, out, nphi)
    
    print *, out

  end subroutine calculate_external_field_ft

end module read_schaffer_field

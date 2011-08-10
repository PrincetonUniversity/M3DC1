module read_schaffer_field
  implicit none

  integer, private :: nr, nz, nphi

  real, private, allocatable :: br(:,:,:), bphi(:,:,:), bz(:,:,:)
  real, private, allocatable :: r(:), phi(:), z(:)
  complex, private, allocatable :: br_ft(:,:), bphi_ft(:,:), bz_ft(:,:)

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
    real :: phi1, r1, z1, br1, bz1, bphi1
    real :: phi0, r0, z0
    character(len=10) :: dummy
    

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ier)

    ! First, let rank zero parse the file to determine the size of the data
    ierr = 0
    if(rank.eq.0) then 
       open(ifile, name=filename, status='old', action='read', err=200)

       do i=1, header_lines
          read(ifile,*) dummy
       end do

       ! determine number of R's, Z's, and Phi's
       nz = 1
       nr = 1
       nphi = 1
       read(ifile,'(3F20.11)') phi0, r0, z0 
       do
          ! read line
          read(ifile,'(3F20.11)',end=100) phi1, r1, z1 
          
          ! skip empty lines
          if(phi1.eq.0. .and. r1.eq.0. .and. z1.eq.0.) cycle

          ! if phi has changed, we have another plane
          if(phi1 .ne. phi0) then
             phi0 = phi1
             nphi = nphi + 1
             nz = 1
             nr = 1
          end if
          
          ! if z has changed, we have another row
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

       write(*, '(A,3I5)') &
            'Reading external fields with nr, nz, nphi =', nr, nz, nphi
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

    ! Allocate space for data
    allocate(br(nphi,nr,nz),bphi(nphi,nr,nz),bz(nphi,nr,nz))
    allocate(r(nr),phi(nphi),z(nz))

    ! Now read in the data
    if(rank.eq.0) then
       ! position read at start of file again
       rewind(ifile,err=2000)

       do i=1, header_lines 
          read(ifile,*) dummy
       end do

       do k=1, nphi
          do j=1, nz
             do i=1, nr
                ! read line
999             read(ifile,'(6F20.11)',end=1000) phi1, r1, z1, bphi1, br1, bz1

                ! skip empty lines
                if(phi1.eq.0. .and. r1.eq.0. .and. z1.eq.0.) goto 999

                ! put data in arrays
                br(k,i,j) = br1
                bphi(k,i,j) = bphi1
                bz(k,i,j) = bz1
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
    call MPI_Bcast(br,   nr*nz*nphi, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bz,   nr*nz*nphi, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bphi, nr*nz*nphi, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  end subroutine load_schaffer_field



  subroutine unload_schaffer_field
    implicit none

    deallocate(br,bz,bphi,r,phi,z)
    if(allocated(br_ft)) deallocate(br_ft)
    if(allocated(bz_ft)) deallocate(bz_ft)
    if(allocated(bphi_ft)) deallocate(bphi_ft)
  end subroutine unload_schaffer_field


  subroutine calculate_external_field_ft(ntor)
    implicit none

    integer, intent(in) :: ntor

    integer :: i, j
    real, dimension(nphi) :: in
    complex, dimension(nphi/2+1) :: out

    if(ntor .ge. nphi/2+1) then
       print *, 'Not enough toroidal data for ntor = ', ntor
       return
    end if

    allocate(br_ft(nr,nz))
    allocate(bz_ft(nr,nz))
    allocate(bphi_ft(nr,nz))

    do i=1, nr
       do j=1, nz
          in = br(:,i,j)
          call fftw(in, out, nphi)
          br_ft(i,j) = out(abs(ntor)+1)/nphi
          in = bz(:,i,j)
          call fftw(in, out, nphi)
          bz_ft(i,j) = out(abs(ntor)+1)/nphi
          in = bphi(:,i,j)
          call fftw(in, out, nphi)
          bphi_ft(i,j) = out(abs(ntor)+1)/nphi
       end do
    end do

    if(ntor .lt. 0) then
       br_ft = conjg(br_ft)
       bz_ft = conjg(bz_ft)
       bphi_ft = conjg(bphi_ft)
    end if
  end subroutine calculate_external_field_ft

  subroutine get_external_field_ft(r1, z1, br_out, bphi_out, bz_out, npts)
    implicit none

    integer, intent(in) :: npts
    real, intent(in), dimension(npts) :: r1, z1
    complex, intent(out), dimension(npts) :: br_out, bphi_out, bz_out

    integer :: i0, j0, p
    real :: ai, aj, di, dj, didr, djdz
    real, dimension(4,4) :: are, aim
    complex, dimension(4,4) :: a

    didr = (nr-1)/(r(nr) - r(1))
    djdz = (nz-1)/(z(nz) - z(1))

    do p=1,npts
       ai = (nr-1)*(r1(p) - r(1))/(r(nr) - r(1)) + 1
       aj = (nz-1)*(z1(p) - z(1))/(z(nz) - z(1)) + 1
       i0 = ai
       j0 = aj

       di = ai - i0
       dj = aj - j0

       call bicubic_interpolation_coeffs(real(br_ft),nr,nz,i0,j0,are)
       call bicubic_interpolation_coeffs(aimag(br_ft),nr,nz,i0,j0,aim)
       a = cmplx(are, aim)
       
       br_out(p) = (a(1,1) + a(1,2)*dj + a(1,3)*dj**2 + a(1,4)*dj**3)       &
            +      (a(2,1) + a(2,2)*dj + a(2,3)*dj**2 + a(2,4)*dj**3)*di    &
            +      (a(3,1) + a(3,2)*dj + a(3,3)*dj**2 + a(3,4)*dj**3)*di**2 &
            +      (a(4,1) + a(4,2)*dj + a(4,3)*dj**2 + a(4,4)*dj**3)*di**3

       call bicubic_interpolation_coeffs(real(bz_ft),nr,nz,i0,j0,are)
       call bicubic_interpolation_coeffs(aimag(bz_ft),nr,nz,i0,j0,aim)
       a = cmplx(are, aim)
       
       bz_out(p) = (a(1,1) + a(1,2)*dj + a(1,3)*dj**2 + a(1,4)*dj**3)       &
            +      (a(2,1) + a(2,2)*dj + a(2,3)*dj**2 + a(2,4)*dj**3)*di    &
            +      (a(3,1) + a(3,2)*dj + a(3,3)*dj**2 + a(3,4)*dj**3)*di**2 &
            +      (a(4,1) + a(4,2)*dj + a(4,3)*dj**2 + a(4,4)*dj**3)*di**3

       call bicubic_interpolation_coeffs(real(bphi_ft),nr,nz,i0,j0,are)
       call bicubic_interpolation_coeffs(aimag(bphi_ft),nr,nz,i0,j0,aim)
       a = cmplx(are, aim)
       
       bphi_out(p) = (a(1,1) + a(1,2)*dj + a(1,3)*dj**2 + a(1,4)*dj**3)       &
            +        (a(2,1) + a(2,2)*dj + a(2,3)*dj**2 + a(2,4)*dj**3)*di    &
            +        (a(3,1) + a(3,2)*dj + a(3,3)*dj**2 + a(3,4)*dj**3)*di**2 &
            +        (a(4,1) + a(4,2)*dj + a(4,3)*dj**2 + a(4,4)*dj**3)*di**3
    end do
  end subroutine get_external_field_ft

end module read_schaffer_field

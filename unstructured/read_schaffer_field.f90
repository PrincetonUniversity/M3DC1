module read_schaffer_field
  implicit none

  integer, private :: nr, nz, nphi

  real, private, allocatable :: br(:,:,:), bphi(:,:,:), bz(:,:,:)
  real, private, allocatable :: r(:), phi(:), z(:)
  complex, private, allocatable :: br_ft(:,:), bphi_ft(:,:), bz_ft(:,:)

contains

  subroutine load_schaffer_field(filename, ierr)
    use math
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
          phi(k) = phi1*pi/180.
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

  subroutine get_external_field(r1, phi1, z1, br_out, bphi_out, bz_out, npts)
    implicit none

    include 'mpif.h'

    integer, intent(in) :: npts
    real, intent(in), dimension(npts) :: r1, phi1, z1
    real, intent(out), dimension(npts) :: br_out, bphi_out, bz_out

    integer :: j, k, l, m, n, ierr, rank
    integer :: i(3)
    real :: x(3), dx(3)
    real :: f(3,4), g(3,4), h(3,4)

    logical :: out_of_bounds = .false.

    do k=1, npts
       x(1) = (nr   - 1)*(r1(k) - r(1))/(r(nr) - r(1)) + 1
       x(2) = (nphi - 1)*(phi1(k) - phi(1))/(phi(nphi) - phi(1)) + 1
       x(3) = (nz   - 1)*(z1(k) - z(1))/(z(nz) - z(1)) + 1

       ! choose the indicies 
       i(:) = x(:) + 1e-5

       ! if the index is out of bounds, set flag
       if(i(1).lt.1 .or. i(1).gt.nr .or. i(3).lt.1 .or. i(3).gt.nz) then
          out_of_bounds = .true.
       end if

       ! if the index is out of bounds, extrapolate
       if(i(1).lt.2) i(1) = 2
       if(i(1).gt.nr-2) i(1) = nr - 2
       if(i(3).lt.2) i(3) = 2
       if(i(3).gt.nz-2) i(3) = nz - 2

       dx(:) = x(:) - i(:)

       ! do tri-cubic interpolation
       do l=1, 4
          do m=1, 4
             j = i(2)+m-2
             if(j.lt.1) j = j+nphi
             if(j.gt.nphi) j = j-nphi

             do n=1, 4
                if(j.lt.1 .or. j.gt.nphi) print *, 'Error! j=', j
                if(i(1)+n-2.lt.1 .or. i(1)+n-2.gt.nr) &
                     print *, 'Error! i(1)=', i(1)
                if(i(3)+l-2.lt.1 .or. i(3)+l-2.gt.nz) &
                     print *, 'Error! i(3)=', i(3)
                f(1,n) = br  (j,i(1)+n-2,i(3)+l-2)
                f(2,n) = bphi(j,i(1)+n-2,i(3)+l-2)
                f(3,n) = bz  (j,i(1)+n-2,i(3)+l-2)
             end do
!             g(:,m) = (1.-dx(1))*f(:,2) + dx(1)*f(:,3)
!             g(:,m) = f(:,2) &
!                  + dx(1)*((f(:,3)-f(:,1))/2. &
!                  + dx(1)*((f(:,3)+f(:,1))/2. - f(:,2))) 
             g(:,m) = f(:,1) &
                  + dx(1)*((-2.*f(:,1) - 3.*f(:,2) + 6.*f(:,3) - f(:,4))/6. &
                  + dx(1)*((f(:,1) - 2.*f(:,2) + f(:,3))/2. &
                  + dx(1)*(-f(:,1) + 3.*f(:,2) - 3.*f(:,3) + f(:,4))/6.))
          end do
          h(:,l) = (1.-dx(2))*g(:,2) + dx(2)*g(:,3)
!          h(:,l) = g(:,2) &
!               + dx(2)*((g(:,3)-g(:,1))/2. &
!               + dx(2)*((g(:,3)+g(:,1))/2. - g(:,2))) 
!          h(:,l) = g(:,1) &
!               + dx(2)*((-2.*g(:,1) - 3.*g(:,2) + 6.*g(:,3) - g(:,4))/6. &
!               + dx(2)*((g(:,1) - 2.*g(:,2) + g(:,3))/2. &
!               + dx(2)*(-g(:,1) + 3.*g(:,2) - 3.*g(:,3) + g(:,4))/6.))
       end do
!       f(:,1) = (1.-dx(3))*h(:,2) + dx(3)*h(:,3)
!       f(:,1)= h(:,2) &
!            + dx(3)*((h(:,3)-h(:,1))/2. &
!            + dx(3)*((h(:,3)+h(:,1))/2. - h(:,2))) 
       f(:,1) = h(:,1) &
            + dx(3)*((-2.*h(:,1) - 3.*h(:,2) + 6.*h(:,3) - h(:,4))/6. &
            + dx(3)*((h(:,1) - 2.*h(:,2) + h(:,3))/2. &
            + dx(3)*(-h(:,1) + 3.*h(:,2) - 3.*h(:,3) + h(:,4))/6.))
       br_out(k)   = f(1,1)
       bphi_out(k) = f(2,1)
       bz_out(k)   = f(3,1)
    end do

    if(out_of_bounds) then
       call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
       if(rank.eq.0) &
            print *, 'Warning: some external field values extrapolated.'
    endif
  end subroutine get_external_field

end module read_schaffer_field

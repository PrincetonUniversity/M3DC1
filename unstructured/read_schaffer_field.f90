module read_schaffer_field
  implicit none

  type schaffer_field
     integer :: nr, nz, nphi

     real, allocatable :: br(:,:,:), bphi(:,:,:), bz(:,:,:), p(:,:,:)
     real, allocatable :: r(:), phi(:), z(:)
     complex, allocatable :: br_ft(:,:), bphi_ft(:,:), bz_ft(:,:)

     logical :: initialized = .false.
     logical :: vmec
  end type schaffer_field

contains

  subroutine load_schaffer_field(sf, filename, isamp, isamp_pol, ierr)
    use math
    implicit none

    include 'mpif.h'

    type(schaffer_field), intent(inout) :: sf
    character(len=*), intent(in) :: filename
    integer, intent(in) :: isamp, isamp_pol
    integer, intent(out) :: ierr

    integer, parameter :: header_lines = 0
    integer, parameter :: ifile = 37
    integer :: i, j, k, l, m, n, ier
    integer :: rank
    real :: phi1, r1, z1, br1, bz1, bphi1, p1
    real :: phi0, r0, z0
    integer, parameter :: catch = 100
    character(len=20) :: dummy
    integer :: file_nr, file_nz, file_nphi
    integer :: lines_read

    sf%vmec = .false.
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ier)

    ! First, let rank zero parse the file to determine the size of the data
    ierr = 0
    if(rank.eq.0) then 
       print *, 'Reading ASCII probe_g file ', filename
       open(ifile, file=filename, status='old', action='read', err=200)
       do i=1, catch
          read(ifile, '(A)', err=200, end=200) dummy
          if(dummy.eq.'             phi_tor') goto 90
          if(dummy.eq.'            %phi_tor') goto 90
          if(dummy.eq.'       %phi_tor(deg)') goto 90
          if(dummy.eq.'# R[m]              ') then
             sf%vmec = .true.
             goto 90
          end if
       end do
       goto 200
90     continue
       if(sf%vmec) then 
          print *, 'External field format is VMEC style'
       else
          print *, 'External field format is PROBE_g style'
       end if
       do i=1, header_lines-1
          read(ifile,*) dummy
       end do
       print *, 'last line = ', dummy

       ! determine number of R's, Z's, and Phi's
       file_nz = 1
       file_nr = 1
       file_nphi = 1
       if(sf%vmec) then
          read(ifile,*) r0, phi0, z0 
       else
          read(ifile,'(3F20.11)') phi0, r0, z0 
       end if
       print *, 'phi0, r0, z0: ', phi0, r0, z0
       do
          ! read line
          if(sf%vmec) then
             read(ifile,*,err=100,end=100) r1, phi1, z1 
             phi1 = phi1*180./pi
          else
             read(ifile,'(3F20.11)',err=100,end=100) phi1, r1, z1 
          end if
          
          ! skip empty lines
          if(phi1.eq.0. .and. r1.eq.0. .and. z1.eq.0.) cycle

          ! if phi has changed, we have another plane
          if(phi1 .ne. phi0) then
             phi0 = phi1
             file_nphi = file_nphi + 1
             file_nz = 1
             file_nr = 1
          end if
          
          ! if z has changed, we have another row
          if(z1 .ne. z0) then
             z0 = z1
             file_nz = file_nz + 1
             file_nr = 1
          end if
          file_nr = file_nr + 1
       end do

100    continue
       file_nr = file_nr-1
       file_nz = file_nz-1
       if(sf%vmec) file_nphi = file_nphi-1

       write(*, '(A,3I5)') &
            'Reading external fields with nr, nz, nphi =', &
            file_nr, file_nz, file_nphi
    end if

    goto 300
200 continue
    ierr = 1
    print *, 'Error reading ', filename
300 continue

    ! Tell all processors whether rank 0 had an error
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    if(ierr.ne.0) return
    
    if(rank.eq.0) then 
       sf%nphi = (file_nphi-1) / isamp + 1
       sf%nr = (file_nr-1) / isamp_pol + 1
       sf%nz = (file_nz-1) / isamp_pol + 1
       write(*, '(A,3I5)') &
            'Downsampling to ', sf%nr, sf%nz, sf%nphi
    end if
       
    ! Send size data to all processors
    call MPI_Bcast(sf%nr,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%nz,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%nphi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%vmec, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)

    ! Allocate space for data
    if(.not. sf%initialized) then
       allocate(sf%r(sf%nr),sf%phi(sf%nphi),sf%z(sf%nz))
       allocate(sf%br(sf%nphi,sf%nr,sf%nz), &
            sf%bphi(sf%nphi,sf%nr,sf%nz),sf%bz(sf%nphi,sf%nr,sf%nz))
       sf%br = 0.
       sf%bphi = 0.
       sf%bz = 0.

       if(sf%vmec) then
          allocate(sf%p(sf%nphi,sf%nr,sf%nz))
          sf%p = 0.
       end if
    endif

    ! Now read in the data
    if(rank.eq.0) then
       ! position read at start of file again
       lines_read = 0
       rewind(ifile,err=2000)

       do
          read(ifile, '(A)', err=2000, end=2000) dummy
          lines_read = lines_read + 1
          if(dummy.eq.'             phi_tor') goto 900
          if(dummy.eq.'            %phi_tor') goto 900
          if(dummy.eq.'       %phi_tor(deg)') goto 900
          if(dummy.eq.'# R[m]              ') goto 900
       end do
       goto 2000
900    continue
       do i=1, header_lines-1
          read(ifile,*) dummy
          lines_read = lines_read + 1
       end do
       print *, 'last line = ', dummy

       do k=1, file_nphi
          do j=1, file_nz
             do i=1, file_nr
                ! read line
999             if(sf%vmec) then
                   read(ifile,*,err=1000,end=1000) &
                        r1, phi1, z1, br1, bphi1, bz1, p1
                   phi1 = phi1*180./pi
                else
                   read(ifile,'(6F20.11)',err=1000,end=1000) &
                        phi1, r1, z1, bphi1, br1, bz1
                end if
                lines_read = lines_read + 1

                ! skip empty lines
                if(phi1.eq.0. .and. r1.eq.0. .and. z1.eq.0.) goto 999

                if(mod(k-1,isamp).eq.0 .and. &
                     mod(j-1,isamp_pol).eq.0 .and. &
                     mod(i-1,isamp_pol).eq.0) then 
                   l = (k-1)/isamp + 1
                   m = (i-1)/isamp_pol + 1
                   n = (j-1)/isamp_pol + 1

                   ! put data in arrays
                   sf%br(l,m,n) = sf%br(l,m,n) + br1
                   sf%bphi(l,m,n) = sf%bphi(l,m,n) + bphi1
                   sf%bz(l,m,n) = sf%bz(l,m,n) + bz1
                   if(sf%vmec) sf%p(l,m,n) = sf%p(l,m,n) + p1
                   sf%r(m) = r1
                end if
             end do
             if(mod(j-1,isamp_pol).eq.0) then 
                n = (j-1)/isamp_pol + 1
                sf%z(n) = z1
             end if
          end do
          if(mod(k-1,isamp).eq.0) then 
             l = (k-1)/isamp + 1
             sf%phi(l) = phi1*pi/180.
          end if
       end do

       goto 1100
1000   continue
       print *, "Shouldn't be here!!!", i, j, k
1100   close(ifile)
       if(sf%vmec) sf%phi = sf%phi - sf%phi(1)
       print *, 'Number of lines read: ', lines_read
       write(*, '(A,3I5)') 'NR, NZ, NPHI', sf%nr, sf%nz, sf%nphi
       write(*, '(A,F12.4," -- ",F12.4)') 'R:   ', sf%r(1), sf%r(sf%nr)
       write(*, '(A,F12.4," -- ",F12.4)') 'Z:   ', sf%z(1), sf%z(sf%nz)
       write(*, '(A,F12.4," -- ",F12.4)') 'Phi: ', sf%phi(1), sf%phi(sf%nphi)
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
    call MPI_Bcast(sf%r,   sf%nr,   MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%z,   sf%nz,   MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%phi, sf%nphi, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%br,   sf%nr*sf%nz*sf%nphi, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%bz,   sf%nr*sf%nz*sf%nphi, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_Bcast(sf%bphi, sf%nr*sf%nz*sf%nphi, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    if(sf%vmec) then
       call MPI_Bcast(sf%p,   sf%nr*sf%nz*sf%nphi, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    end if

    sf%initialized = .true.
    if(rank.eq.0) print *, 'Done reading fields.'
  end subroutine load_schaffer_field



  subroutine unload_schaffer_field(sf)
    implicit none

    type(schaffer_field) :: sf

    deallocate(sf%br,sf%bz,sf%bphi,sf%r,sf%phi,sf%z)
    if(allocated(sf%br_ft)) deallocate(sf%br_ft)
    if(allocated(sf%bz_ft)) deallocate(sf%bz_ft)
    if(allocated(sf%bphi_ft)) deallocate(sf%bphi_ft)
    if(allocated(sf%p)) deallocate(sf%p)

    sf%initialized = .false.
  end subroutine unload_schaffer_field


  subroutine calculate_external_field_ft(sf, ntor)
    implicit none

    type(schaffer_field) :: sf
    integer, intent(in) :: ntor

    integer :: i, j
    real, dimension(sf%nphi) :: in
    complex, dimension(sf%nphi/2+1) :: out

    if(ntor .ge. sf%nphi/2+1) then
       print *, 'Not enough toroidal data for ntor = ', ntor
       return
    end if

    allocate(sf%br_ft(sf%nr,sf%nz))
    allocate(sf%bz_ft(sf%nr,sf%nz))
    allocate(sf%bphi_ft(sf%nr,sf%nz))

    call fftw_init(in, out, sf%nphi)
    do i=1, sf%nr
       do j=1, sf%nz
          in = sf%br(:,i,j)
          call fftw_run
          sf%br_ft(i,j) = out(abs(ntor)+1)/sf%nphi
          in = sf%bz(:,i,j)
          call fftw_run
          sf%bz_ft(i,j) = out(abs(ntor)+1)/sf%nphi
          in = sf%bphi(:,i,j)
          call fftw_run
          sf%bphi_ft(i,j) = out(abs(ntor)+1)/sf%nphi
       end do
    end do
!    call fftw_destroy

    if(ntor .lt. 0) then
       sf%br_ft = conjg(sf%br_ft)
       sf%bz_ft = conjg(sf%bz_ft)
       sf%bphi_ft = conjg(sf%bphi_ft)
    end if
  end subroutine calculate_external_field_ft

  subroutine get_external_field_ft(sf, r1, z1, br_out, bphi_out, bz_out, npts)
    implicit none

    type(schaffer_field), intent(in) :: sf
    integer, intent(in) :: npts
    real, intent(in), dimension(npts) :: r1, z1
    complex, intent(out), dimension(npts) :: br_out, bphi_out, bz_out

    integer :: i0, j0, p
    real :: ai, aj, di, dj, didr, djdz
    complex, dimension(4,4) :: a

    didr = (sf%nr-1)/(sf%r(sf%nr) - sf%r(1))
    djdz = (sf%nz-1)/(sf%z(sf%nz) - sf%z(1))

    do p=1,npts
       ai = (sf%nr-1)*(r1(p) - sf%r(1))/(sf%r(sf%nr) - sf%r(1)) + 1
       aj = (sf%nz-1)*(z1(p) - sf%z(1))/(sf%z(sf%nz) - sf%z(1)) + 1
       i0 = ai
       j0 = aj

       di = ai - i0
       dj = aj - j0

       call bicubic_interpolation_coeffs_complex(sf%br_ft,sf%nr,sf%nz,i0,j0,a)

       
       br_out(p) = (a(1,1) + a(1,2)*dj + a(1,3)*dj**2 + a(1,4)*dj**3)       &
            +      (a(2,1) + a(2,2)*dj + a(2,3)*dj**2 + a(2,4)*dj**3)*di    &
            +      (a(3,1) + a(3,2)*dj + a(3,3)*dj**2 + a(3,4)*dj**3)*di**2 &
            +      (a(4,1) + a(4,2)*dj + a(4,3)*dj**2 + a(4,4)*dj**3)*di**3

       call bicubic_interpolation_coeffs_complex(sf%bz_ft,sf%nr,sf%nz,i0,j0,a)
       
       bz_out(p) = (a(1,1) + a(1,2)*dj + a(1,3)*dj**2 + a(1,4)*dj**3)       &
            +      (a(2,1) + a(2,2)*dj + a(2,3)*dj**2 + a(2,4)*dj**3)*di    &
            +      (a(3,1) + a(3,2)*dj + a(3,3)*dj**2 + a(3,4)*dj**3)*di**2 &
            +      (a(4,1) + a(4,2)*dj + a(4,3)*dj**2 + a(4,4)*dj**3)*di**3

       call bicubic_interpolation_coeffs_complex(sf%bphi_ft,sf%nr,sf%nz,i0,j0,a)
       
       bphi_out(p) = (a(1,1) + a(1,2)*dj + a(1,3)*dj**2 + a(1,4)*dj**3)       &
            +        (a(2,1) + a(2,2)*dj + a(2,3)*dj**2 + a(2,4)*dj**3)*di    &
            +        (a(3,1) + a(3,2)*dj + a(3,3)*dj**2 + a(3,4)*dj**3)*di**2 &
            +        (a(4,1) + a(4,2)*dj + a(4,3)*dj**2 + a(4,4)*dj**3)*di**3
    end do
  end subroutine get_external_field_ft

  subroutine get_external_field(sf,r1,phi1,z1,br_out,bphi_out,bz_out,p_out, &
       npts)
    implicit none

    include 'mpif.h'

    type(schaffer_field), intent(in) :: sf
    integer, intent(in) :: npts
    real, intent(in), dimension(npts) :: r1, phi1, z1
    real, dimension(npts) :: phim
    real, intent(out), dimension(npts) :: br_out, bphi_out, bz_out, p_out

    integer :: j, k, l, m, n, ierr, rank
    integer :: i(3)
    real :: x(3), dx(3)
    real :: f(4,4), g(4,4), h(4,4)

    logical :: out_of_bounds = .false.
    real :: obr, obz

    ! Modulo on phi in case data is only one field period
    phim = mod(phi1 - sf%phi(1),(sf%phi(sf%nphi) - sf%phi(1))) + sf%phi(1)

    do k=1, npts
       x(1) = (sf%nr   - 1)*(r1(k) - sf%r(1))/(sf%r(sf%nr) - sf%r(1)) + 1
       x(2) = (sf%nphi - 1)*(phim(k) - sf%phi(1)) &
            /(sf%phi(sf%nphi) - sf%phi(1)) + 1
       x(3) = (sf%nz   - 1)*(z1(k) - sf%z(1))/(sf%z(sf%nz) - sf%z(1)) + 1

       ! choose the indicies 
       i(:) = x(:)

       ! if the index is out of bounds, set flag
       if(i(1).lt.1 .or. i(1).gt.sf%nr .or. i(3).lt.1 .or. i(3).gt.sf%nz) then
          out_of_bounds = .true.
          obr = r1(k)
          obz = z1(k)
       end if

       ! if the index is out of bounds, extrapolate
       if(i(1).lt.2) i(1) = 2
       if(i(1).gt.sf%nr-2) i(1) = sf%nr - 2
       if(i(3).lt.2) i(3) = 2
       if(i(3).gt.sf%nz-2) i(3) = sf%nz - 2

       dx(:) = x(:) - i(:)

       ! do tri-cubic interpolation
       do l=1, 4
          do m=1, 4
             j = i(2)+m-2
             if(j.lt.1) j = j+sf%nphi
             if(j.gt.sf%nphi) j = j-sf%nphi

             do n=1, 4
                if(j.lt.1 .or. j.gt.sf%nphi) &
                     print *, 'Error! j=', j, sf%nphi
                if(i(1)+n-2.lt.1 .or. i(1)+n-2.gt.sf%nr) &
                     print *, 'Error! i(1)=', i(1)
                if(i(3)+l-2.lt.1 .or. i(3)+l-2.gt.sf%nz) &
                     print *, 'Error! i(3)=', i(3)
                f(1,n) = sf%br  (j,i(1)+n-2,i(3)+l-2)
                f(2,n) = sf%bphi(j,i(1)+n-2,i(3)+l-2)
                f(3,n) = sf%bz  (j,i(1)+n-2,i(3)+l-2)
                if(sf%vmec) f(4,n) = sf%p(j,i(1)+n-2,i(3)+l-2)
             end do
!             g(:,m) = (1.-dx(1))*f(:,2) + dx(1)*f(:,3)
             g(:,m) = f(:,2) &
                  + dx(1)*((f(:,3)-f(:,1))/2. &
                  + dx(1)*((f(:,3)+f(:,1))/2. - f(:,2))) 
!             g(:,m) = f(:,1) &
!                  + dx(1)*((-2.*f(:,1) - 3.*f(:,2) + 6.*f(:,3) - f(:,4))/6. &
!                  + dx(1)*((f(:,1) - 2.*f(:,2) + f(:,3))/2. &
!                  + dx(1)*(-f(:,1) + 3.*f(:,2) - 3.*f(:,3) + f(:,4))/6.))
          end do
!          h(:,l) = (1.-dx(2))*g(:,2) + dx(2)*g(:,3)
          h(:,l) = g(:,2) &
               + dx(2)*((g(:,3)-g(:,1))/2. &
               + dx(2)*((g(:,3)+g(:,1))/2. - g(:,2))) 
!          h(:,l) = g(:,1) &
!               + dx(2)*((-2.*g(:,1) - 3.*g(:,2) + 6.*g(:,3) - g(:,4))/6. &
!               + dx(2)*((g(:,1) - 2.*g(:,2) + g(:,3))/2. &
!               + dx(2)*(-g(:,1) + 3.*g(:,2) - 3.*g(:,3) + g(:,4))/6.))
       end do
!       f(:,1) = (1.-dx(3))*h(:,2) + dx(3)*h(:,3)
       f(:,1)= h(:,2) &
            + dx(3)*((h(:,3)-h(:,1))/2. &
            + dx(3)*((h(:,3)+h(:,1))/2. - h(:,2))) 
!       f(:,1) = h(:,1) &
!            + dx(3)*((-2.*h(:,1) - 3.*h(:,2) + 6.*h(:,3) - h(:,4))/6. &
!            + dx(3)*((h(:,1) - 2.*h(:,2) + h(:,3))/2. &
!            + dx(3)*(-h(:,1) + 3.*h(:,2) - 3.*h(:,3) + h(:,4))/6.))
       br_out(k)   = f(1,1)
       bphi_out(k) = f(2,1)
       bz_out(k)   = f(3,1)
       if(sf%vmec) then
          p_out(k) = max(f(4,1),0.)
       else
          p_out(k) = 0
       end if
    end do

    if(out_of_bounds) then
       call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
       if(rank.eq.0) &
            print *, 'Warning: some external field values extrapolated. ', &
            obr, obz
    endif
  end subroutine get_external_field

  subroutine load_fieldlines_field(sf, filename, error)
    use hdf5
    implicit none

    type(schaffer_field), intent(inout) :: sf
    character(len=*), intent(in) :: filename
    integer, intent(out) :: error
    integer :: k, lpres 
    integer(HID_T) :: file_id, dset_id, attr_id
    integer(HSIZE_T), dimension(1) :: dim0 = 1
    integer(HSIZE_T), dimension(1) :: dim1 
    integer(HSIZE_T), dimension(3) :: dim3
    real, allocatable :: br1(:,:,:), bphi1(:,:,:), bz1(:,:,:), p1(:,:,:)

    error = 0

    call h5open_f(error)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
    if(error.lt.0) then
       return
    end if

    call h5dopen_f(file_id, 'lpres', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, lpres, dim0, error)

    if(lpres.eq.0) then
       sf%vmec = .false.
    else 
       sf%vmec = .true.
    end if

    ! read array sizes nr, nphi, and nz 
    call h5dopen_f(file_id, 'nr', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, sf%nr, dim0, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'nphi', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, sf%nphi, dim0, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'nz', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, sf%nz, dim0, error)
    call h5dclose_f(dset_id, error)

    if(.not. sf%initialized) then
       allocate(sf%r(sf%nr))
       allocate(sf%z(sf%nz))
       allocate(sf%phi(sf%nphi))
       allocate(sf%br(sf%nphi,sf%nr,sf%nz))
       allocate(sf%bphi(sf%nphi,sf%nr,sf%nz))
       allocate(sf%bz(sf%nphi,sf%nr,sf%nz))
       allocate(br1(sf%nr,sf%nphi,sf%nz))
       allocate(bphi1(sf%nr,sf%nphi,sf%nz))
       allocate(bz1(sf%nr,sf%nphi,sf%nz))
       if(lpres.eq.1) then
          allocate(sf%p(sf%nphi,sf%nr,sf%nz))
          allocate(p1(sf%nr,sf%nphi,sf%nz))
       end if
    end if
    ! read 1d arrays r, phi, and z
    dim1(1) = sf%nr 
    call h5dopen_f(file_id, 'raxis', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sf%r, dim1, error)
    call h5dclose_f(dset_id, error)
    dim1(1) = sf%nz 
    call h5dopen_f(file_id, 'zaxis', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sf%z, dim1, error)
    call h5dclose_f(dset_id, error)
    dim1(1) = sf%nphi 
    call h5dopen_f(file_id, 'phiaxis', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sf%phi, dim1, error)
    call h5dclose_f(dset_id, error)
    ! read 3d arrays br, bphi, and bz
    dim3(1) = sf%nr 
    dim3(2) = sf%nphi 
    dim3(3) = sf%nz
    call h5dopen_f(file_id, 'B_R', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, br1, dim3, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'B_PHI', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, bphi1, dim3, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'B_Z', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, bz1, dim3, error)
    call h5dclose_f(dset_id, error)
    if(lpres.eq.1) then
       call h5dopen_f(file_id, 'PRES', dset_id, error)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, p1, dim3, error)
       call h5dclose_f(dset_id, error)

       if(error .lt. 0) then
          print *, 'ERROR: lpres = 1 but PRES not found'
          sf%vmec = .false.
          error = 0
       end if
    end if
    do k = 1, sf%nz
       sf%br(:,:,k) = transpose(br1(:,:,k))
       sf%bphi(:,:,k) = transpose(bphi1(:,:,k))
       sf%bz(:,:,k) = transpose(bz1(:,:,k))
       if(lpres.eq.1) sf%p(:,:,k) = transpose(p1(:,:,k))
    end do 
    deallocate(br1,bz1,bphi1)
    if(lpres.eq.1) deallocate(p1)

    sf%initialized = .true.

  end subroutine load_fieldlines_field

  subroutine load_mips_field(sf, filename, error)
    use math
    implicit none

    type(schaffer_field), intent(inout) :: sf
    character(len=*), intent(in) :: filename
    integer, intent(out) :: error
    real, allocatable :: br1(:,:,:), bphi1(:,:,:), bz1(:,:,:), p1(:,:,:)
    integer::lrnet, lznet, lphinet
    integer::lphinet4
    integer::i, j, k
    real(8)::rmin, rmax, zmin, zmax, phimin, phimax
    real(8)::dr, dz, dphi
    real(8), allocatable::br(:,:,:), bz(:,:,:), bphi(:,:,:), prs(:,:,:)

    open(60,file=filename, form='unformatted',convert='BIG_ENDIAN')

    read(60) lrnet, lznet, lphinet4 
    close(60)
    lphinet = lphinet4 - 4 

    allocate(br(lrnet,lznet,lphinet4))
    allocate(bz(lrnet,lznet,lphinet4))
    allocate(bphi(lrnet,lznet,lphinet4))
    allocate(prs(lrnet,lznet,lphinet4))

    open(60,file=filename, form='unformatted',convert='BIG_ENDIAN')
    read(60) lrnet, lznet, lphinet4, &
         rmin, rmax, zmin, zmax, phimin, phimax, &
         dr, dz, dphi, &
         br, bz, bphi, prs
    close(60)
    prs = prs/(pi*4e-7)

    sf%nr = lrnet
    sf%nz = lrnet
    sf%nphi = lphinet+1
    if(.not. sf%initialized) then
       allocate(sf%r(sf%nr))
       allocate(sf%z(sf%nz))
       allocate(sf%phi(sf%nphi))
       allocate(sf%br(sf%nphi,sf%nr,sf%nz))
       allocate(sf%bphi(sf%nphi,sf%nr,sf%nz))
       allocate(sf%bz(sf%nphi,sf%nr,sf%nz))
       allocate(br1(sf%nr,sf%nphi,sf%nz))
       allocate(bphi1(sf%nr,sf%nphi,sf%nz))
       allocate(bz1(sf%nr,sf%nphi,sf%nz))
       allocate(sf%p(sf%nphi,sf%nr,sf%nz))
       allocate(p1(sf%nr,sf%nphi,sf%nz))
    end if

    do i = 1, sf%nr
       sf%r(i) = rmin + dr*(i-1)
    end do

    do j = 1, sf%nphi
       sf%phi(j) = phimin + dphi*(j-1) - .5*(phimax-phimin)
       br1(:,j,:) = br(:,:,j+2)
       bphi1(:,j,:) = bphi(:,:,j+2)
       bz1(:,j,:) = bz(:,:,j+2)
       p1(:,j,:) = prs(:,:,j+2)
    end do

    do k = 1, sf%nz
       sf%z(k) = zmin + dz*(k-1)
       sf%br(:,:,k) = transpose(br1(:,:,k))
       sf%bphi(:,:,k) = transpose(bphi1(:,:,k))
       sf%bz(:,:,k) = transpose(bz1(:,:,k))
       sf%p(:,:,k) = transpose(p1(:,:,k))
    end do 
    deallocate(br1,bz1,bphi1,p1)
    deallocate(br,bphi,bz,prs)

    sf%vmec = .true.

    sf%initialized = .true.

  end subroutine load_mips_field

#ifdef USEST
  subroutine check(istatus)
  implicit none

  integer, intent (in) :: istatus
  if (istatus.ne.0) then
     call safestop(53)
  end if
  end subroutine check

  subroutine load_mgrid_field(sf, mgrid_filename, vmec_filename, error)
    use iso_c_binding
    use math
    implicit none

    type(schaffer_field), intent(inout) :: sf
    character(len=*), intent(in) :: mgrid_filename, vmec_filename
    integer, intent(out) :: error

    integer :: ll, ii, kk
    integer(c_int) :: ncid, ncid_vmec, istatus
    real :: curfac, dr, dz, dphi
    real :: per

! Dimension IDs
    integer(c_int) :: radDimID, phiDimID, zeeDimID
! Variable IDs
    integer(c_int) :: rminID, rmaxID, zminID, zmaxID, nfpID, nextcurID, modeID, extcurID, rawcurID
! Variable values
    real :: rmin, rmax, zmin, zmax
    real, allocatable :: extcur(:)
    integer :: nfp, nextcur
    character :: mgrid_mode
! Magnetic field
    integer(c_int) :: varID
    character(len=30) :: var
    real, allocatable :: brtemp(:,:,:), bphitemp(:,:,:), bztemp(:,:,:), temp(:,:,:)

    integer(c_int), parameter :: NC_NOWRITE = 0
    integer(c_size_t) :: dimlen
    character(kind=c_char, len=257) :: c_filename

    interface
      function nc_open(path, omode, ncidp) bind(C, name="nc_open") result(r)
        use iso_c_binding
        integer(c_int) :: r
        character(c_char), intent(in) :: path(*)
        integer(c_int), value :: omode
        integer(c_int), intent(out) :: ncidp
      end function

      function nc_inq_dimid(ncid, name, idp) bind(C, name="nc_inq_dimid") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid
        character(c_char), intent(in) :: name(*)
        integer(c_int), intent(out) :: idp
      end function

      function nc_inq_dimlen(ncid, dimid, lenp) bind(C, name="nc_inq_dimlen") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, dimid
        integer(c_size_t), intent(out) :: lenp
      end function

      function nc_inq_varid(ncid, name, varidp) bind(C, name="nc_inq_varid") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid
        character(c_char), intent(in) :: name(*)
        integer(c_int), intent(out) :: varidp
      end function

      function nc_get_var_int(ncid, varid, ip) bind(C, name="nc_get_var_int") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, varid
        integer(c_int), intent(out) :: ip
      end function

      function nc_get_var_double(ncid, varid, dp) bind(C, name="nc_get_var_double") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, varid
        real(c_double), intent(out) :: dp
      end function

      function nc_get_var_text(ncid, varid, tp) bind(C, name="nc_get_var_text") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, varid
        character(c_char), intent(out) :: tp(*)
      end function
    end interface

    error = 0

! Check extension
    ll = len_trim(mgrid_filename)
    if (mgrid_filename(ll-2:ll).eq.'.nc') then
       c_filename = trim(mgrid_filename) // c_null_char
       call check(int(nc_open(c_filename, NC_NOWRITE, ncid)))

       ! Get dimensions
       call check(int(nc_inq_dimid(ncid, "rad"//c_null_char, radDimID)))
       call check(int(nc_inq_dimid(ncid, "phi"//c_null_char, phiDimID)))
       call check(int(nc_inq_dimid(ncid, "zee"//c_null_char, zeeDimID)))
       call check(int(nc_inq_dimlen(ncid, radDimID, dimlen))); sf%nr   = int(dimlen)
       call check(int(nc_inq_dimlen(ncid, phiDimID, dimlen))); sf%nphi = int(dimlen)
       call check(int(nc_inq_dimlen(ncid, zeeDimID, dimlen))); sf%nz   = int(dimlen)

       ! Get variable values
       call check(int(nc_inq_varid(ncid, "rmin"//c_null_char, rminID)))
       call check(int(nc_inq_varid(ncid, "rmax"//c_null_char, rmaxID)))
       call check(int(nc_inq_varid(ncid, "zmin"//c_null_char, zminID)))
       call check(int(nc_inq_varid(ncid, "zmax"//c_null_char, zmaxID)))
       call check(int(nc_inq_varid(ncid, "nfp"//c_null_char, nfpID)))
       call check(int(nc_inq_varid(ncid, "nextcur"//c_null_char, nextcurID)))
       call check(int(nc_inq_varid(ncid, "mgrid_mode"//c_null_char, modeID)))

       call check(int(nc_get_var_double(ncid, rminID, rmin)))
       call check(int(nc_get_var_double(ncid, rmaxID, rmax)))
       call check(int(nc_get_var_double(ncid, zminID, zmin)))
       call check(int(nc_get_var_double(ncid, zmaxID, zmax)))
       call check(int(nc_get_var_int(ncid, nfpID, nfp)))
       call check(int(nc_get_var_int(ncid, nextcurID, nextcur)))
       call check(int(nc_get_var_text(ncid, modeID, mgrid_mode)))

       ! Check mgrid_mode
       allocate(extcur(nextcur))
       if (mgrid_mode.eq.'S') then
          print *, 'Coil currents are SCALED...'
          print *, 'Opening vmec file'
          c_filename = trim(vmec_filename) // c_null_char
          call check(int(nc_open(c_filename, NC_NOWRITE, ncid_vmec)))
          istatus = nc_inq_varid(ncid_vmec, "extcur"//c_null_char, extcurID)
          if(istatus.eq.0) then
             print *, 'Reading external currents from VMEC file.'
             call check(int(nc_get_var_double(ncid_vmec, extcurID, extcur(1))))
          else
             print *, 'External currents not found in VMEC file.'
             print *, 'Looking for external currents in MGRID file.'
             call check(int(nc_inq_varid(ncid, "raw_coil_cur"//c_null_char, rawcurID)))
             print *, 'Reading external currents from MGRID file.'
             call check(int(nc_get_var_double(ncid, rawcurID, extcur(1))))
             print *, 'extcur = ', extcur
          end if
       else if (mgrid_mode.eq.'R') then
          print *, 'Actual currents supplied'
          extcur = 1.0
       else
          print *, 'ERROR: Invalid coil currents'
          call safestop(54)
       end if

! Allocate temporary arrays
       allocate(brtemp(sf%nr,sf%nz,sf%nphi)) ! According to libstell (r,z,phi)
       allocate(bphitemp(sf%nr,sf%nz,sf%nphi))
       allocate(bztemp(sf%nr,sf%nz,sf%nphi))
       brtemp = 0.0
       bphitemp = 0.0
       bztemp = 0.0

! Get fields
       allocate(temp(sf%nr,sf%nz,sf%nphi))
8001   format(a,'_',i3.3)
       do ii = 1, nextcur, 1
          write(var,8001) 'br', ii
          call check(int(nc_inq_varid(ncid, trim(var)//c_null_char, varID)))
          call check(int(nc_get_var_double(ncid, varID, temp(1,1,1))))
          brtemp = brtemp + extcur(ii)*temp

          write(var,8001) 'bp', ii
          call check(int(nc_inq_varid(ncid, trim(var)//c_null_char, varID)))
          call check(int(nc_get_var_double(ncid, varID, temp(1,1,1))))
          bphitemp = bphitemp + extcur(ii)*temp

          write(var,8001) 'bz', ii
          call check(int(nc_inq_varid(ncid, trim(var)//c_null_char, varID)))
          call check(int(nc_get_var_double(ncid, varID, temp(1,1,1))))
          bztemp = bztemp + extcur(ii)*temp
       end do
       deallocate(temp)

! Schaffer field needs 1 more toroidal grid point than MGRID!
       sf%nphi = sf%nphi + 1
       if(.not. sf%initialized) then
         allocate(sf%br(sf%nphi,sf%nr,sf%nz))
         allocate(sf%bphi(sf%nphi,sf%nr,sf%nz))
         allocate(sf%bz(sf%nphi,sf%nr,sf%nz))
         allocate(sf%r(sf%nr))
         allocate(sf%z(sf%nz))
         allocate(sf%phi(sf%nphi))
       end if

! Transpose (r,z,phi) -> (phi,r,z)
       do kk = 1, sf%nphi-1
          sf%br(kk,:,:) = brtemp(:,:,kk)
          sf%bphi(kk,:,:) = bphitemp(:,:,kk)
          sf%bz(kk,:,:) = bztemp(:,:,kk)
       end do
! Data on extra toroidal grid point
       sf%br(sf%nphi,:,:) = brtemp(:,:,1)
       sf%bphi(sf%nphi,:,:) = bphitemp(:,:,1)
       sf%bz(sf%nphi,:,:) = bztemp(:,:,1)
       deallocate(brtemp,bztemp,bphitemp)

! Make grid
       dr = (rmax-rmin)/(sf%nr-1)
       dz = (zmax-zmin)/(sf%nz-1)
       per = twopi/nfp
       dphi = per/(sf%nphi-1)

       do kk = 1, sf%nr, 1
          sf%r(kk) = rmin + (kk-1)*dr
       end do
       do kk = 1, sf%nphi, 1
          sf%phi(kk) = 0.0 + (kk-1)*dphi
       end do
       do kk = 1, sf%nz, 1
          sf%z(kk) = zmin + (kk-1)*dz
       end do
       sf%initialized = .true.
    else
       print *, 'ERROR: Invalid extension. MGRID must be .nc'
       call safestop(52)
    end if
  end subroutine load_mgrid_field

  subroutine load_hint_field(sf, hint_filename, error)
    use iso_c_binding
    use math
    implicit none

    type(schaffer_field), intent(inout) :: sf
    character(len=*), intent(in) :: hint_filename
    integer, intent(out) :: error

    integer :: ll, ii, kk
    integer(c_int) :: ncid
    real :: dr, dz, dphi
    real :: per

! Dimension IDs
    integer(c_int) :: radDimID, phiDimID, zeeDimID, timeDimID
! Variable IDs
    integer(c_int) :: rminID, rmaxID, zminID, zmaxID, nfpID
! Variable values
    real :: rmin, rmax, zmin, zmax
    integer :: nfp, kstep
! Magnetic field
    integer(c_int) :: br_varid, bp_varid, bz_varid, p_varid
    real, allocatable :: brtemp(:,:,:), bphitemp(:,:,:), bztemp(:,:,:), ptemp(:,:,:)

    integer(c_int), parameter :: NC_NOWRITE = 0
    integer(c_size_t) :: dimlen
    integer(c_size_t) :: start_c(4), count_c(4)
    character(kind=c_char, len=257) :: c_filename

    interface
      function nc_open(path, omode, ncidp) bind(C, name="nc_open") result(r)
        use iso_c_binding
        integer(c_int) :: r
        character(c_char), intent(in) :: path(*)
        integer(c_int), value :: omode
        integer(c_int), intent(out) :: ncidp
      end function

      function nc_close(ncid) bind(C, name="nc_close") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid
      end function

      function nc_inq_dimid(ncid, name, idp) bind(C, name="nc_inq_dimid") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid
        character(c_char), intent(in) :: name(*)
        integer(c_int), intent(out) :: idp
      end function

      function nc_inq_dimlen(ncid, dimid, lenp) bind(C, name="nc_inq_dimlen") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, dimid
        integer(c_size_t), intent(out) :: lenp
      end function

      function nc_inq_varid(ncid, name, varidp) bind(C, name="nc_inq_varid") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid
        character(c_char), intent(in) :: name(*)
        integer(c_int), intent(out) :: varidp
      end function

      function nc_get_var_int(ncid, varid, ip) bind(C, name="nc_get_var_int") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, varid
        integer(c_int), intent(out) :: ip
      end function

      function nc_get_var_double(ncid, varid, dp) bind(C, name="nc_get_var_double") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, varid
        real(c_double), intent(out) :: dp
      end function

      function nc_get_vara_double(ncid, varid, startp, countp, ip) &
          bind(C, name="nc_get_vara_double") result(r)
        use iso_c_binding
        integer(c_int) :: r
        integer(c_int), value :: ncid, varid
        integer(c_size_t), intent(in) :: startp(*), countp(*)
        real(c_double), intent(out) :: ip(*)
      end function
    end interface

    error = 0

! Check extension
    ll = len_trim(hint_filename)
    if (hint_filename(ll-2:ll).eq.'.nc') then
       c_filename = trim(hint_filename) // c_null_char
       call check(int(nc_open(c_filename, NC_NOWRITE, ncid)))

! Get dimensions
       call check(int(nc_inq_dimid(ncid, "R"//c_null_char, radDimID)))
       call check(int(nc_inq_dimid(ncid, "phi"//c_null_char, phiDimID)))
       call check(int(nc_inq_dimid(ncid, "Z"//c_null_char, zeeDimID)))
       call check(int(nc_inq_dimid(ncid, "time"//c_null_char, timeDimID)))
       call check(int(nc_inq_dimlen(ncid, radDimID,  dimlen))); sf%nr   = int(dimlen)
       call check(int(nc_inq_dimlen(ncid, phiDimID,  dimlen))); sf%nphi = int(dimlen)
       call check(int(nc_inq_dimlen(ncid, zeeDimID,  dimlen))); sf%nz   = int(dimlen)
       call check(int(nc_inq_dimlen(ncid, timeDimID, dimlen))); kstep   = int(dimlen)

! Get variable values
       call check(int(nc_inq_varid(ncid, "rminb"//c_null_char, rminID)))
       call check(int(nc_inq_varid(ncid, "rmaxb"//c_null_char, rmaxID)))
       call check(int(nc_inq_varid(ncid, "zminb"//c_null_char, zminID)))
       call check(int(nc_inq_varid(ncid, "zmaxb"//c_null_char, zmaxID)))
       call check(int(nc_inq_varid(ncid, "mtor"//c_null_char, nfpID)))

       call check(int(nc_get_var_double(ncid, rminID, rmin)))
       call check(int(nc_get_var_double(ncid, rmaxID, rmax)))
       call check(int(nc_get_var_double(ncid, zminID, zmin)))
       call check(int(nc_get_var_double(ncid, zmaxID, zmax)))
       call check(int(nc_get_var_int(ncid, nfpID, nfp)))

! Allocate temporary arrays
       allocate(brtemp(sf%nr,sf%nz,sf%nphi))
       allocate(bphitemp(sf%nr,sf%nz,sf%nphi))
       allocate(bztemp(sf%nr,sf%nz,sf%nphi))
       allocate(ptemp(sf%nr,sf%nz,sf%nphi))

! Get fields (adapted from HINT source code)
       call check(int(nc_inq_varid(ncid, "B_R"//c_null_char,   br_varid)))
       call check(int(nc_inq_varid(ncid, "B_phi"//c_null_char, bp_varid)))
       call check(int(nc_inq_varid(ncid, "B_Z"//c_null_char,   bz_varid)))
       call check(int(nc_inq_varid(ncid, "P"//c_null_char,     p_varid)))

       ! C API uses C-order (reversed) dimensions and 0-based indices.
       ! Fortran dims: (nr, nz, nphi, time) -> C dims: (time, nphi, nz, nr)
       count_c = [1_c_size_t, int(sf%nphi,c_size_t), int(sf%nz,c_size_t), int(sf%nr,c_size_t)]
       start_c = [int(kstep-1,c_size_t), 0_c_size_t, 0_c_size_t, 0_c_size_t]

       call check(int(nc_get_vara_double(ncid, br_varid, start_c, count_c, brtemp)))
       call check(int(nc_get_vara_double(ncid, bp_varid, start_c, count_c, bphitemp)))
       call check(int(nc_get_vara_double(ncid, bz_varid, start_c, count_c, bztemp)))
       call check(int(nc_get_vara_double(ncid, p_varid,  start_c, count_c, ptemp)))
       call check(int(nc_close(ncid)))
! Pressure normalization
       ptemp = ptemp/(pi*4e-7)
! Schaffer field needs 1 more toroidal grid point than HINT!
       sf%nphi = sf%nphi + 1

       if(.not. sf%initialized) then
          allocate(sf%r(sf%nr))
          allocate(sf%z(sf%nz))
          allocate(sf%phi(sf%nphi))
          allocate(sf%br(sf%nphi,sf%nr,sf%nz))
          allocate(sf%bphi(sf%nphi,sf%nr,sf%nz))
          allocate(sf%bz(sf%nphi,sf%nr,sf%nz))
          allocate(sf%p(sf%nphi,sf%nr,sf%nz))
       end if

! Transpose (r,z,phi) -> (phi,r,z)
       do kk = 1, sf%nphi-1
          sf%br(kk,:,:) = brtemp(:,:,kk)
          sf%bphi(kk,:,:) = bphitemp(:,:,kk)
          sf%bz(kk,:,:) = bztemp(:,:,kk)
          sf%p(kk,:,:) = ptemp(:,:,kk)
       end do
! Data on extra toroidal grid point
       sf%br(sf%nphi,:,:) = brtemp(:,:,1)
       sf%bphi(sf%nphi,:,:) = bphitemp(:,:,1)
       sf%bz(sf%nphi,:,:) = bztemp(:,:,1)
       sf%p(sf%nphi,:,:) = ptemp(:,:,1)
       deallocate(brtemp,bztemp,bphitemp,ptemp)

! Make grid
       dr = (rmax-rmin)/(sf%nr-1)
       dz = (zmax-zmin)/(sf%nz-1)
       per = twopi/nfp
       dphi = per/(sf%nphi-1)

       do kk = 1, sf%nr, 1
          sf%r(kk) = rmin + (kk-1)*dr
       end do
       do kk = 1, sf%nphi, 1
          sf%phi(kk) = 0.0 + (kk-1)*dphi
       end do
       do kk = 1, sf%nz, 1
          sf%z(kk) = zmin + (kk-1)*dz
       end do
! Finalize
       sf%initialized = .true.
       sf%vmec = .true.
    else
       print *, 'ERROR: Invalid extension. Only support .nc format'
       call safestop(52)
    end if
  end subroutine load_hint_field

#endif
end module read_schaffer_field

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
       write(*, '(A,3I0)') 'NR, NZ, NPHI', sf%nr, sf%nz, sf%nphi
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

    do i=1, sf%nr
       do j=1, sf%nz
          in = sf%br(:,i,j)
          call fftw(in, out, sf%nphi)
          sf%br_ft(i,j) = out(abs(ntor)+1)/sf%nphi
          in = sf%bz(:,i,j)
          call fftw(in, out, sf%nphi)
          sf%bz_ft(i,j) = out(abs(ntor)+1)/sf%nphi
          in = sf%bphi(:,i,j)
          call fftw(in, out, sf%nphi)
          sf%bphi_ft(i,j) = out(abs(ntor)+1)/sf%nphi
       end do
    end do

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
    real, intent(out), dimension(npts) :: br_out, bphi_out, bz_out, p_out

    integer :: j, k, l, m, n, ierr, rank
    integer :: i(3)
    real :: x(3), dx(3)
    real :: f(4,4), g(4,4), h(4,4)

    logical :: out_of_bounds = .false.
    real :: obr, obz

    do k=1, npts
       x(1) = (sf%nr   - 1)*(r1(k) - sf%r(1))/(sf%r(sf%nr) - sf%r(1)) + 1
       x(2) = (sf%nphi - 1)*(phi1(k) - sf%phi(1)) &
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
          p_out(k) = f(4,1)
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

end module read_schaffer_field

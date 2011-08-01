program m3dc1_fortran_test
  implicit none

  integer :: i, j, k, ierr
  real :: r, phi, z, br, bphi, bz, n, p
  integer :: n_handle, p_handle
  integer :: timeslice = 1
  integer :: nr, nz, nphi
  integer, parameter :: MAX_SIZE = 100000000

  integer, parameter :: bout = 51
  integer, parameter :: pout = 52
  integer, parameter :: nout = 53
  
  real :: rmin, rmax, zmin, zmax, phimin, phimax
  real :: dr, dphi, dz
  real :: factor = 1

  character(len=32) :: argv

  common nr, nz, nphi, rmin, rmax, zmin, zmax, phimin, phimax

  integer :: argc
  argc = command_argument_count()

  if(argc .lt. 3) then
     call print_usage
     stop
  end if

  call get_command_argument(1,argv)
  read(argv,'(i10)') nr
  call get_command_argument(2,argv)
  read(argv,'(i10)') nz
  call get_command_argument(3,argv)
  read(argv,'(i10)') nphi

  if((nr.lt.2) .or. (nz.lt.2) .or. (nphi.lt.2)) then
     print *, 'Error: invalid number of points'
     stop
  end if
  if(nr*nz*nphi.gt.MAX_SIZE) then
     print *, 'Error: number of points exceeds MAX_SIZE'
     print *, '       To increase limit, raise MAX_SIZE and recompile.'
     stop
  end if


  ! Open the file
  ! (note that we're calling a c function 
  !  so the string must be null-terminated)
  call m3dc1_open_file("C1.h5"//char(0), ierr)
  if(ierr .ne. 0) then
     print *, 'Error opening file.'
     stop
  end if

  call m3dc1_load_magnetic_field(timeslice, ierr)
  if(ierr .ne. 0) then
     print *, 'Error loading magnetic field.'
     goto 100
  end if

  call m3dc1_load_field("den"//char(0), timeslice, n_handle, ierr)
  if(ierr.ne.0) then
     print *, 'Error loading density field.'
     goto 100
  end if

  call m3dc1_load_field("P"//char(0), timeslice, p_handle, ierr)
  if(ierr.ne.0) then
     print *, 'Error loading pressure field.'
     goto 100
  end if

  write(*,'(A,3I6)') "NR, NZ, NPHI = ", NR, NZ, NPHI
  write(*,'(A,1E13.5)') "Factor = ", factor

  call m3dc1_extent(timeslice, rmin, rmax, phimin, phimax, zmin, zmax, ierr)

  phimin = 0.
  phimax = 2.*3.14159265*(nphi-1)/nphi
  dr = (rmax - rmin)/(nr-1)
  dphi = (phimax - phimin)/(nphi-1)
  dz = (zmax - zmin)/(nz-1)

  open(unit=bout, file='B.out', action='WRITE', status='REPLACE')
  open(unit=pout, file='p.out', action='WRITE', status='REPLACE')
  open(unit=nout, file='n.out', action='WRITE', status='REPLACE')

  call write_header(bout)
  call write_header(pout)
  call write_header(nout)

  do i=1, nphi
     phi = dphi*(i-1) + phimin
     write(*, '("Plane ",I5,";  Phi = ",E13.5)') i, phi
     do j=1, nz
        z = dz*(j-1) + zmin
        do k=1, nr
           r = dr*(k-1) + rmin

           n = 0
           p = 0
           br = 0
           bphi = 0
           bz = 0

           call m3dc1_eval_magnetic_field(r, phi, z, br, bphi, bz, ierr)
           write(bout, '(3E13.5)') br, bz, bphi

           call m3dc1_eval_field(n_handle, r, phi, z, n, ierr)
           write(nout, '(1E13.5)') n

           call m3dc1_eval_field(p_handle, r, phi, z, p, ierr)
           write(nout, '(1E13.5)') p
        end do
     end do
  end do

  close(bout)
  close(pout)
  close(nout)
  
100 continue
  call m3dc1_close_file()

end program m3dc1_fortran_test

subroutine write_header(ifile)
  implicit none
  integer, intent(in) :: ifile

  integer :: nr, nz, nphi
  real :: rmin, rmax, zmin, zmax, phimin, phimax
  common nr, nz, nphi, rmin, rmax, zmin, zmax, phimin, phimax

  write(ifile,'(3I5)')   nr, nz, nphi
  write(ifile,'(6E13.5)') rmin, rmax, zmin, zmax, phimin, phimax
end subroutine write_header

subroutine print_usage
  implicit none

  print *, "Usage:  m3dc1_convert NR NZ NPHI"
  print *, " NR, NZ, NPHI: Number of points in R, Z, and PHI directions"
  print *, ""
  print *, "Output (B.out) is in following format: "
  print *, " NR NZ NPHI"
  print *, " RMIN RMAX ZMIN ZMAX PHIMIN PHIMAX"
  print *, " BR1 BZ1 BPHI1"
  print *, " BR2 BZ2 BPHI2"
  print *, " BR3 BZ3 BPHI3"
  print *, " ..."
  print *, ""
  print *, "where"
  print *, " RMIN, RMAX: minimum and maximum R coordinate"
  print *, " ZMIN, ZMAX: minimum and maximum Z coordinate"
  print *, " PHIMIN, PHIMAX: minimum and maximum PHI coordinate"
  print *, " BR1 BZ1 BPHI1: components of field at point 1"
  print *, ""
  print *, "points are ordered in the following way:"
  print *, " i = ir + iz*NR + iphi*NR*NZ"
  print *, "coordinates of point i are"
  print *, " (R, Phi, Z) = "
  print *, " (ir  *(  RMAX-  RMIN)/(NR  -1) +   RMIN, "
  print *, "  iphi*(PHIMAX-PHIMIN)/(NPHI-1) + PHIMIN, "
  print *, "  iz  *(  ZMAX-  ZMIN)/(NZ  -1) +   ZMIN)"
end subroutine print_usage

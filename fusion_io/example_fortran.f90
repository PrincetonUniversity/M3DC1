!==========================================================================
! fio_example
! ~~~~~~~~~~~
! In this example, several files from linear calculations are loaded.
! The fields are evaluated at several points, and the results are summed.
!==========================================================================
program fio_example
  use fusion_io

  ! factor by which to multiply perturbed (linear) part of solution
  real, parameter :: factor = 1.

  ! time slice to read
  integer, parameter :: timeslice = 1

  integer, parameter :: nfiles = 4
  character(len=64) :: filename(nfiles)
  integer, dimension(nfiles) :: isrc, ipres, idens, imag

  real :: p(1), n(1), b(3), x(3)
  real :: t1(1), t3(3)

  integer, parameter ::  npts = 10
  real :: R0, R1, Z0, Z1, phi0, phi1
  integer :: i, j, ierr

  filename(1) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=1/C1.h5'
  filename(2) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=2/C1.h5'
  filename(3) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=3/C1.h5'
  filename(4) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=4/C1.h5'

  ! read files and fields
  do i=1, nfiles
     print *, 'Reading ', filename(i)
     call fio_open_source_f(FIO_M3DC1_SOURCE, trim(filename(i)), isrc(i), ierr)
     if(ierr.ne.0) stop

     ! Set options appropriate to this source
     call fio_get_options_f(isrc(i), ierr)
     call fio_set_int_option_f(FIO_TIMESLICE, timeslice, ierr)
     call fio_set_real_option_f(FIO_LINEAR_SCALE, factor, ierr)

     ! For first file, read full fields (equilibrium + perturbed)
     ! For subsequent file, read only perturbed parts
     if(i.gt.1) call fio_set_int_option_f(FIO_PERTURBED_ONLY, 1, ierr)

     ! read fields
     call fio_get_field_f(isrc(i), FIO_TOTAL_PRESSURE, 0, ipres(i), ierr);
     call fio_get_field_f(isrc(i), FIO_DENSITY, FIO_ELECTRON, idens(i),ierr);
     call fio_get_field_f(isrc(i), FIO_MAGNETIC_FIELD, 0, imag(i), ierr);
  end do
 
  R0 = 1.6;
  R1 = 2.1;
  Z0 = 0.0;
  Z1 = 0.0;
  phi0 = 0.;
  phi1 = 0.;

  do i=1, npts
     x(1) = R0 + (R1-R0)*(i-1)/(npts-1);
     x(2) = phi0 + (phi1-phi0)*(i-1)/(npts-1);
     x(3) = Z0 + (Z1-Z0)*(i-1)/(npts-1);

     write(*, '("(",3F12.4,"):")') x

     p = 0.
     n = 0.
     b = 0.

     do j=1, nfiles
        call fio_eval_field_f(ipres(j), x, t1, ierr)
        p = p + t1
        call fio_eval_field_f(idens(j), x, t1, ierr)
        n = n + t1
        call fio_eval_field_f(imag(j), x, t3, ierr)
        b = b + t3
     end do

     write(*, '("        pressure = ",1pE12.4)') p
     write(*, '("        density = ",1pE12.4)') n
     write(*, '("        final field = ",1p3E12.4)') b
  end do

  do i=1, nfiles
     call fio_close_field_f(idens(i), ierr)
     call fio_close_field_f(ipres(i), ierr)
     call fio_close_field_f(imag(i), ierr)
     call fio_close_source_f(isrc(i), ierr)
  end do
end program fio_example

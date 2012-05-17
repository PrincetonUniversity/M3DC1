program fio_example
  use fusion_io

#include "fusion_io.h"

  integer :: isrc
  integer :: pressure, density, magnetic_field
  integer :: i, npts
  real :: R0, R1, Z0, Z1, phi0, phi1
  integer :: ierr

  real :: p(1), n(1), b(3), x(3)

  call fio_open_source(FIO_M3DC1_SOURCE, "C1.h5", isrc, ierr)
  if(ierr.ne.0) return

  call fio_get_options(isrc, ierr)

  call fio_set_int_option(FIO_TIMESLICE, 1, ierr)
  call fio_set_int_option(FIO_PERTURBED_ONLY, 1, ierr)
  call fio_set_real_option(FIO_LINEAR_SCALE, 10., ierr)

  call fio_get_field(isrc, FIO_PRESSURE, pressure, ierr);
  call fio_get_field(isrc, FIO_DENSITY, density, ierr);
  call fio_get_field(isrc, FIO_MAGNETIC_FIELD, magnetic_field, ierr);

  npts = 10;
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

     call fio_eval_field(pressure, x, p, ierr)
     write(*, '("        pressure = ",1pE12.4)') p

     call fio_eval_field(density, x, n, ierr)
     write(*, '("        pressure = ",1pE12.4)') n

     call fio_eval_field(magnetic_field, x, b, ierr)
     write(*, '("        pressure = ",1p3E12.4)') b
  end do

  call fio_close_field(density, ierr)
  call fio_close_field(pressure, ierr)
  call fio_close_field(magnetic_field, ierr)
  call fio_close_source(isrc, ierr)

end program fio_example

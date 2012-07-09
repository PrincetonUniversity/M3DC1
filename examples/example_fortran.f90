!==========================================================================
! fio_example
! ~~~~~~~~~~~
! In this example, several files from linear calculations are loaded.
! The fields are evaluated at several points, and the results are summed.
!==========================================================================
program fio_example
  use fusion_io

  implicit none

  ! factor by which to multiply perturbed (linear) part of solution
  real, parameter :: factor = 1.

  ! time slice to read
  integer, parameter :: timeslice = 1

  integer, parameter :: nfiles = 4
  character(len=64) :: filename(nfiles), filename_efit, filename_gato
  integer, dimension(nfiles) :: isrc, ipres, ine, ini, imag, ij
  integer :: isrc_efit, imag_efit, ipres_efit, ij_efit
  integer :: isrc_gato, imag_gato, ipres_gato

  real :: p(1), ne(1), ni(1), b(3), x(3), curr(3), curr2(3)
  real :: t1(1), t3(3)

  integer, parameter ::  npts = 10
  real :: R0, R1, Z0, Z1, phi0, phi1
  integer :: i, j, ierr

  filename(1) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=1/C1.h5'
  filename(2) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=2/C1.h5'
  filename(3) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=3/C1.h5'
  filename(4) = '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=4/C1.h5'

  filename_efit = &
       '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=3/geqdsk'

  ! read files and fields
  do i=1, nfiles
     print *, 'Reading ', filename(i)
     call fio_open_source_f(FIO_M3DC1_SOURCE, trim(filename(i)), isrc(i), ierr)
     if(ierr.ne.0) goto 100

     ! Set options appropriate to this source
     call fio_get_options_f(isrc(i), ierr)
     call fio_set_int_option_f(FIO_TIMESLICE, timeslice, ierr)
     call fio_set_real_option_f(FIO_LINEAR_SCALE, factor, ierr)

     ! For first file, read full fields (equilibrium + perturbed)
     ! For subsequent file, read only perturbed parts
     if(i.gt.1) call fio_set_int_option_f(FIO_PERTURBED_ONLY, 1, ierr)

     ! read fields
     ! magnetic field and total pressure are species-independent
     call fio_get_field_f(isrc(i), FIO_TOTAL_PRESSURE, ipres(i), ierr);
     call fio_get_field_f(isrc(i), FIO_MAGNETIC_FIELD, imag(i), ierr);
     call fio_get_field_f(isrc(i), FIO_CURRENT_DENSITY, ij(i), ierr);

     ! density is species dependent; specify electrons
     call fio_set_int_option_f(FIO_SPECIES, FIO_ELECTRON, ierr)
     call fio_get_field_f(isrc(i), FIO_DENSITY, ine(i),ierr);

     call fio_set_int_option_f(FIO_SPECIES, FIO_MAIN_ION, ierr)
     call fio_get_field_f(isrc(i), FIO_DENSITY, ini(i),ierr);
  end do

  ! open efit file
  call fio_open_source_f(FIO_GEQDSK_SOURCE,trim(filename_efit),isrc_efit,ierr)
  if(ierr.ne.0) goto 100

  call fio_get_options_f(isrc_efit, ierr)
  call fio_get_field_f(isrc_efit, FIO_MAGNETIC_FIELD, imag_efit, ierr);
  call fio_get_field_f(isrc_efit, FIO_TOTAL_PRESSURE, ipres_efit, ierr);
  call fio_get_field_f(isrc_efit, FIO_CURRENT_DENSITY, ij_efit, ierr);


  ! open gato file
  filename_gato = '/Users/ferraro/data/GATO/diagnostics.dat'
  call fio_open_source_f(FIO_GATO_SOURCE,trim(filename_gato),isrc_gato,ierr)
  if(ierr.ne.0) goto 100
  call fio_get_options_f(isrc_gato, ierr)
  call fio_get_field_f(isrc_gato, FIO_MAGNETIC_FIELD, imag_gato, ierr);
  call fio_get_field_f(isrc_gato, FIO_TOTAL_PRESSURE, ipres_gato, ierr);
 
  R0 = 1.6;
  R1 = 2.1;
  Z0 = 0.4;
  Z1 = 0.4;
  phi0 = 0.;
  phi1 = 0.;

  do i=1, npts
     x(1) = R0 + (R1-R0)*(i-1)/(npts-1);
     x(2) = phi0 + (phi1-phi0)*(i-1)/(npts-1);
     x(3) = Z0 + (Z1-Z0)*(i-1)/(npts-1);

     write(*, '("(",3F12.4,"):")') x

     p = 0.
     ne = 0.
     ni = 0.
     b = 0.
     curr = 0.

     do j=1, nfiles
        call fio_eval_field_f(ipres(j), x, t1, ierr)
        p = p + t1
        call fio_eval_field_f(ine(j), x, t1, ierr)
        ne = ne + t1
        call fio_eval_field_f(ini(j), x, t1, ierr)
        ni = ni + t1
        call fio_eval_field_f(imag(j), x, t3, ierr)
        b = b + t3
        call fio_eval_field_f(ij(j), x, t3, ierr)
        curr = curr + t3
     end do

     write(*, '("        pressure = ",1pE12.4)') p
     write(*, '("        electron density = ",1pE12.4)') ne
     write(*, '("        ion density = ",1pE12.4)') ni
     write(*, '("        total b = ",1p3E12.4)') b
     write(*, '("        total j = ",1p3E12.4)') curr

     b = 0.
     p = 0.
     curr2 = 0.
     call fio_eval_field_f(ipres_efit, x, p, ierr)
     call fio_eval_field_f(imag_efit, x, b, ierr)
     call fio_eval_field_f(ij_efit, x, curr2, ierr)
     write(*, '("        efit press = ", 1pE12.4)') p
     write(*, '("        efit b = ", 1p3E12.4)') b
     write(*, '("        efit j = ", 1p3E12.4)') curr2

     call fio_eval_field_f(ipres_gato, x, p, ierr)
     call fio_eval_field_f(imag_gato, x, b, ierr)
     write(*, '("        gato press = ", 1pE12.4)') p
     write(*, '("        gato b = ", 1p3E12.4)') b
  end do

100 continue 
  do i=1, nfiles
     call fio_close_field_f(ine(i), ierr)
     call fio_close_field_f(ini(i), ierr)
     call fio_close_field_f(ipres(i), ierr)
     call fio_close_field_f(imag(i), ierr)
     call fio_close_field_f(ij(i), ierr)
     call fio_close_source_f(isrc(i), ierr)
  end do

  call fio_close_field_f(imag_efit, ierr)
  call fio_close_field_f(ipres_efit, ierr)
  call fio_close_field_f(ij_efit, ierr)
  call fio_close_source_f(isrc_efit, ierr)

  call fio_close_field_f(imag_gato, ierr)
  call fio_close_field_f(ipres_gato, ierr)
  call fio_close_source_f(isrc_gato, ierr)
end program fio_example

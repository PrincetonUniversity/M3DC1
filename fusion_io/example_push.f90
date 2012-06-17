program example_push
  use fio_push

  implicit none

  ! type of source file to read
  integer, parameter :: itype = FIO_M3DC1_SOURCE

  ! name of source file
  character(len=*), parameter :: filename = 'C1.h5'

  ! for linear calculations, amount to scale perturbed part
  real, parameter :: linfac = 5.

  ! time (0 = vacuum fields, >0 = plasma response)
  real, parameter :: t = 1.
  
  integer, parameter ::  npts = 10
  real :: R0, R1, Z0, Z1, phi0, phi1
  real, dimension(3) :: q
  type(field_type) :: field
  integer :: i

  ! read source file
  call fio_push_initialize(itype, filename, linfac)

  R0 = 1.6;
  R1 = 2.1;
  Z0 = 0.0;
  Z1 = 0.0;
  phi0 = 0.;
  phi1 = 0.;

  do i=1, npts
     q(1) = R0 + (R1-R0)*(i-1)/(npts-1)           ! R
     q(2) = phi0 + (phi1-phi0)*(i-1)/(npts-1)     ! Phi
     q(3) = Z0 + (Z1-Z0)*(i-1)/(npts-1)           ! Z

     write(*, '("(",3F12.4,"):")') q

     call fio_push_field_eval(t, q, field)

     write(*, '(" Phi       = ",1pE12.4)')  field%phi
     write(*, '(" Grad(Phi) = ",1p3E12.4)') field%gradphi
     write(*, '(" A         = ",1p3E12.4)') field%a
     write(*, '(" Grad(A)   = ",1p3E12.4)') field%grada
  end do

  ! close source file
  call fio_push_finalize()

end program example_push

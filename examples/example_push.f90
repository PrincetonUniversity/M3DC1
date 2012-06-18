program example_push
  use fio_push

  implicit none

  ! type of source file to read
  integer, parameter :: itype = FIO_M3DC1_SOURCE

  ! name of source file
  character(len=*), parameter :: filename = &
       '/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=1/C1.h5'

  ! for linear calculations, amount to scale perturbed part
  real, parameter :: linfac = 5.

  ! time (0 = vacuum fields, >0 = plasma response)
  real, parameter :: t = 1.
  
  real, parameter :: epsilon = 1e-5
  real, parameter, dimension(3) :: q0 = (/ 1.6, 0., 0. /)
  real, dimension(3) :: q
  real, dimension(3) :: a0, a1
  real, dimension(3, 3) :: grada
  type(field_type) :: field
  integer :: i

  ! read source file
  call fio_push_initialize(itype, filename, linfac)

  ! evaluate fields at a point
  write(*, '("At (R, Phi, Z) = (",3F12.4,"):")') q0

  call fio_push_field_eval(t, q0, field)
  
  write(*, '("  Phi       = ",1pE12.4)')  field%phi
  write(*, '("  Grad(Phi) = ",1p3E12.4)') field%gradphi
  write(*, '("  A         = ",1p3E12.4)') field%a
  write(*, '("  Grad(A)   = ",1p3E12.4)') field%grada

  ! numerically test derivatives
  print *, 'Numerical derivatives :'
  do i=1, 3
     q = q0
     q(i) = q0(i) + epsilon
     call fio_push_field_eval(t, q, field)
     a1 = field%a
     q(i) = q0(i) - epsilon
     call fio_push_field_eval(t, q, field)
     a0 = field%a
     grada(i,:) = (a1 - a0)/(2.*epsilon)
  end do

  write(*, '("  Grad(A)   = ",1p3E12.4)') grada

  ! close source file
  call fio_push_finalize()

end program example_push

Program test_region
  use mpi
  use region

  implicit none

  type(region_type) :: reg
  integer :: ierr
  character(len=256) :: boundary_file
  
  integer, parameter :: n = 5
  real, dimension(n) :: r, phi, z
  logical :: in
  integer :: i

  r(1) = 3.;   z(1) = 2.
  r(2) = 5.;   z(2) = 2.
  r(3) = 4.;   z(3) = 5.
  r(4) = 8.;   z(4) = 1.
  r(5) = -2.;  z(5) = -12.
  phi = 0.

  call MPI_Init(ierr)

  boundary_file = "unstructured/templates/ITER/rw2_adapt/outer_boundary.pts"
  if(command_argument_count().ge.1) call get_command_argument(1, boundary_file)

  call create_region_from_file(reg, trim(boundary_file), ierr)
  
  do i=1, n
     in = point_in_region(reg, r(i), phi(i), z(i))
     write(*,'("Is (", G0, ", ", G0, " in region? ", L1)') r(i), z(i), in
  end do

  call destroy_region(reg)

  call MPI_Finalize(ierr)

end Program test_region

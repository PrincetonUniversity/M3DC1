Program test_region
  use region

  implicit none

  type(region_type) :: reg
  integer :: ierr
  
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

  call create_region_from_file(reg, "outer_boundary.pts", ierr)
  
  do i=1, n
     in = point_in_region(reg, r(i), phi(i), z(i))
     write(*,'("Is (", G, ", ", G," in region? ", L)') r(i), z(i), in
  end do

  call destroy_region(reg)

end Program test_region

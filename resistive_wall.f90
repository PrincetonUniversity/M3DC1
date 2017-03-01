module resistive_wall

  use basic

  implicit none

  integer, parameter :: imax_wall_breaks = 20

  real :: eta_wall
  real, dimension(imax_wall_breaks) :: eta_break
  real, dimension(imax_wall_breaks) :: wall_break_xmin, wall_break_xmax
  real, dimension(imax_wall_breaks) :: wall_break_zmin, wall_break_zmax
  real, dimension(imax_wall_breaks) :: wall_break_phimin, wall_break_phimax
  integer :: iwall_breaks

contains

  elemental real function wall_resistivity(x, phi, z)
    implicit none

    real, intent(in) :: x, phi, z

    integer :: i

    wall_resistivity = eta_wall

    do i=1, iwall_breaks
#ifdef USE3D
       if(phi.ge.wall_break_phimin(i) .and. &
            phi.le.wall_break_phimax(i)) then
#endif
          if(x.ge.wall_break_xmin(i) .and. &
               x.le.wall_break_xmax(i) .and. &
               z.ge.wall_break_zmin(i) .and. &
               z.le.wall_break_zmax(i)) then
             wall_resistivity = eta_break(i)
             continue
          end if
#ifdef USE3D
       endif
#endif
    end do

  end function wall_resistivity

end module resistive_wall

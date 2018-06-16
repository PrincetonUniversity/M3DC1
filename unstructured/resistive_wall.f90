module resistive_wall

  use basic
  use region

  implicit none

  integer, parameter :: imax_wall_breaks = 20

  real :: eta_wall
  real, dimension(imax_wall_breaks) :: eta_break
  real, dimension(imax_wall_breaks) :: wall_break_xmin, wall_break_xmax
  real, dimension(imax_wall_breaks) :: wall_break_zmin, wall_break_zmax
  real, dimension(imax_wall_breaks) :: wall_break_phimin, wall_break_phimax
  integer :: iwall_breaks

  integer, parameter :: imax_wall_regions = 20
  integer :: iwall_regions
  character(len=256), dimension(imax_wall_regions) :: wall_region_filename
  real, dimension(imax_wall_regions) :: wall_region_eta
  type(region_type), dimension(imax_wall_regions), private :: wall_region

contains

  subroutine init_resistive_wall(ierr)
    implicit none

    integer, intent(out) :: ierr

    integer :: i

    ierr = 0
    do i=1, iwall_regions
       call create_region_from_file(wall_region(i), wall_region_filename(i), &
            ierr)
       if(ierr.ne.0) return
    end do
  end subroutine init_resistive_wall

  subroutine destroy_resistive_wall
    implicit none

    integer :: i

    do i=1, iwall_regions
       call destroy_region(wall_region(i))
    end do
    iwall_regions = 0
  end subroutine destroy_resistive_wall

  elemental real function wall_resistivity(x, phi, z)
    implicit none

    real, intent(in) :: x, phi, z

    integer :: i

    wall_resistivity = eta_wall

    do i=iwall_regions, 1, -1
       if(point_in_region(wall_region(i), x, phi, z)) then
          wall_resistivity = wall_region_eta(i)
          exit
       end if
    end do

    do i=iwall_breaks, 1, -1
#ifdef USE3D
       if(phi.ge.wall_break_phimin(i) .and. &
            phi.le.wall_break_phimax(i)) then
#endif
          if(x.ge.wall_break_xmin(i) .and. &
               x.le.wall_break_xmax(i) .and. &
               z.ge.wall_break_zmin(i) .and. &
               z.le.wall_break_zmax(i)) then
             wall_resistivity = eta_break(i)
             exit
          end if
#ifdef USE3D
       endif
#endif
    end do

  end function wall_resistivity

end module resistive_wall

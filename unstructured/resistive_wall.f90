module resistive_wall

  use basic
  use region
  use mesh_mod

  implicit none

  integer, parameter :: imax_wall_breaks = 20

  real :: eta_wall, eta_wallRZ
  real, dimension(imax_wall_breaks) :: eta_break
  real, dimension(imax_wall_breaks) :: wall_break_xmin, wall_break_xmax
  real, dimension(imax_wall_breaks) :: wall_break_zmin, wall_break_zmax
  real, dimension(imax_wall_breaks) :: wall_break_phimin, wall_break_phimax
  integer :: iwall_breaks

  integer, parameter :: imax_wall_regions = 20
  integer :: iwall_regions
  character(len=256), dimension(imax_wall_regions) :: wall_region_filename
  real, dimension(imax_wall_regions) :: wall_region_eta, wall_region_etaRZ
  type(region_type), dimension(imax_wall_regions) :: wall_region

  real, dimension(max_zones) :: eta_zone, etaRZ_zone

  real :: eta_rekc
  integer :: ntor_rekc, mpol_rekc, isym_rekc
  real :: phi_rekc, theta_rekc, sigma_rekc
  real :: rzero_rekc, zzero_rekc

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

  elemental real function wall_resistivity(x, phi, z, izone)

    use math

    implicit none

    real, intent(in) :: x, phi, z
    integer, intent(in) :: izone
    real :: theta, f, f1, f2
    integer :: i

    if(eta_zone(izone).gt.0.) then
       wall_resistivity = eta_zone(izone)
    else
       wall_resistivity = eta_wall
    end if

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

    if(eta_rekc.gt.0) then
       theta = atan2(z-zzero_rekc, x-rzero_rekc)
      f = 0
#ifdef USE3D
       f = ntor_rekc*(phi-phi_rekc)*twopi/toroidal_period
#endif
       if(isym_rekc .eq. 0) then
          f = cos(f - mpol_rekc*(theta-theta_rekc))
          f = exp((f-1.)/sigma_rekc**2)
          wall_resistivity = 10.**(log10(wall_resistivity)*(1.-f) + log10(eta_rekc)*f)
       else
          f1 = cos(f - mpol_rekc*(theta-theta_rekc))
          f1 = exp((f1-1.)/sigma_rekc**2)
          f2 = cos(f + mpol_rekc*(theta-theta_rekc))
          f2 = exp((f2-1.)/sigma_rekc**2)
          wall_resistivity = 10.**(log10(wall_resistivity)*(1.-max(f1,f2)) + log10(eta_rekc)*max(f1,f2))
      endif
    end if

  end function wall_resistivity

  elemental real function wall_resistivityRZ(x, phi, z, izone)

    use math

    implicit none

    real, intent(in) :: x, phi, z
    integer, intent(in) :: izone
    real :: theta, f, f1, f2
    integer :: i

    if(etaRZ_zone(izone).gt.0.) then
       wall_resistivityRZ = etaRZ_zone(izone)
    else
       wall_resistivityRZ = eta_wallRZ
    end if

    do i=iwall_regions, 1, -1
       if(point_in_region(wall_region(i), x, phi, z)) then
          wall_resistivityRZ = wall_region_etaRZ(i)
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
             wall_resistivityRZ = eta_break(i)
             exit
          end if
#ifdef USE3D
       endif
#endif
    end do

    if(eta_rekc.gt.0) then
       theta = atan2(z-zzero_rekc, x-rzero_rekc)
       f = 0
#ifdef USE3D
       f = ntor_rekc*(phi-phi_rekc)*twopi/toroidal_period
#endif
       if(isym_rekc .eq. 0) then
          f = cos(f - mpol_rekc*(theta-theta_rekc))
          f = exp((f-1.)/sigma_rekc**2)
          wall_resistivityRZ = 10.**(log10(wall_resistivityRZ)*(1.-f) + log10(eta_rekc)*f)
       else
          f1 = cos(f - mpol_rekc*(theta-theta_rekc))
          f1 = exp((f1-1.)/sigma_rekc**2)
          f2 = cos(f + mpol_rekc*(theta-theta_rekc))
          f2 = exp((f2-1.)/sigma_rekc**2)
          wall_resistivityRZ = 10.**(log10(wall_resistivityRZ)*(1.-max(f1,f2)) + log10(eta_rekc)*max(f1,f2))
      endif
    end if
  end function wall_resistivityRZ


end module resistive_wall

module neutral_beam

  implicit none

  integer :: ibeam  ! 1 = include neutral beam source
  real :: beam_x    ! x coordinate of pellet injection
  real :: beam_z    ! z coordinate of pellet injection
  real :: beam_v    ! beam voltage (in volts)
  real :: beam_rate ! amplitude of pellet density source
  real :: beam_dr   ! spatial dispersion of beam source
  real :: beam_dv   ! voltage dispersion of beam source (in volts)
  real, parameter :: beam_zeff = 1

  real :: nb_r, nb_z, nb_v, nb_dr, nb_dz, nb_dv, nb_n

contains

  subroutine neutral_beam_init()
    use basic
    implicit none

    nb_n = beam_rate*t0_norm/(n0_norm*l0_norm**3)
    nb_r = beam_x/(l0_norm/100.)
    nb_z = beam_z/(l0_norm/100.)
    nb_dr = beam_dr/(l0_norm/100.)
    nb_v = sqrt(2.*beam_zeff*e_c*(abs(beam_v)*1.e8/c_light)/(ion_mass*m_p)) &
         /v0_norm
    if(beam_v.lt.0.) nb_v = -nb_v
    nb_dv =sqrt(2.*beam_zeff*e_c*(beam_dv*1.e8/c_light)/(ion_mass*m_p))/v0_norm
  end subroutine neutral_beam_init
  

  vectype elemental function neutral_beam_deposition(r, z)
    use math
    use basic
    implicit none
    real, intent(in) :: r, z
    
    neutral_beam_deposition = nb_n/(2.*r*(pi**2)*nb_dr**2) &
          *exp(-((r-nb_r)**2 + (z-nb_z)**2)/(2.*nb_dr**2))
  end function neutral_beam_deposition

end module neutral_beam

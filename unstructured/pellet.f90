module pellet
  implicit none

  integer :: ipellet  ! 1 = include pellet injection density source
                      ! 2 = distributed density source
  real :: pellet_x    ! x coordinate of pellet
  real :: pellet_phi  ! phi coordinate of pellet
  real :: pellet_z    ! z coordinate of pellet
  real :: pellet_rate ! amplitude of pellet density source
  real :: pellet_var  ! spatial dispersion of density source 

  real :: pellet_velx 
  real :: pellet_velphi
  real :: pellet_velz
  
contains

  subroutine pellet_init()
    use basic
    implicit none

    ! convert to m/s
    pellet_velx = pellet_velx/v0_norm * l0_norm
    pellet_velphi = pellet_velphi/v0_norm * l0_norm
    pellet_velz = pellet_velz/v0_norm * l0_norm
  end subroutine pellet_init
  
  vectype elemental function pellet_deposition(r, phi, z, pres, den)
    use math
    use basic
    implicit none
    real, intent(in) :: r, phi, z, pres, den

    select case(ipellet)

    ! gaussian pellet source
    case(1)
#ifdef USE3D
       pellet_deposition = pellet_rate/(sqrt(2.*pi)*pellet_var)**3 & 
            *exp(-(r**2 + pellet_x**2 -2.*x*pellet_x*cos(phi-pellet_phi) &
            + (z - pellet_z)**2) / (2.*pellet_var**2))
#else
       pellet_deposition = pellet_rate/(r*2.*pi*pellet_var**2) & 
            *exp(-((r - pellet_x)**2 + (z - pellet_z)**2) &
            /(2.*pellet_var**2))
#endif

    !......distributed source added 11/23/2011   (scj)
    case(2)
       pellet_deposition = pellet_rate*den0*(max(pedge,real(pres))/p0)**expn

    ! gaussian pellet source
    case(3)
#ifdef USE3D
       pellet_deposition = pres*pellet_rate/(sqrt(2.*pi)*pellet_var)**3 & 
            *exp(-(r**2 + pellet_x**2 -2.*x*pellet_x*cos(phi-pellet_phi) &
            + (z - pellet_z)**2) / (2.*pellet_var**2))
#else
       pellet_deposition = pres*pellet_rate/(r*2.*pi*pellet_var**2) & 
            *exp(-((r - pellet_x)**2 + (z - pellet_z)**2) &
            /(2.*pellet_var**2))
#endif

    case default
       pellet_deposition = 0
    end select
  end function pellet_deposition

  subroutine pellet_advance
    use basic
    implicit none
    
    pellet_x = pellet_x + pellet_velx*dt
    pellet_z = pellet_z + pellet_velz*dt
    pellet_phi = pellet_phi + pellet_velphi*dt
  end subroutine pellet_advance

end module pellet

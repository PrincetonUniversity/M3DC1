module pellet
  implicit none

 ! integer :: ipellet  ! 1 = include pellet injection density source
                       ! 2 = distributed density source
                       ! 3 = Gaussian source
                       ! 4 = Gaussian source including realistic ablation model  
 ! real :: pellet_x    ! x coordinate of pellet
 ! real :: pellet_phi  ! phi coordinate of pellet
 ! real :: pellet_z    ! z coordinate of pellet
 ! real :: pellet_rate ! amplitude of pellet density source
 ! real :: pellet_var  ! spatial dispersion of density source 
 ! real :: pellet_var_tor  ! spatial dispersion of density source 

  real :: pellet_velx 
  real :: pellet_velphi
  real :: pellet_velz
  
contains

  subroutine pellet_init()
     use basic
     use diagnostics
     implicit none

    if(pellet_var_tor.le.0) pellet_var_tor = pellet_var

    ! convert to m/s
    pellet_velx = pellet_velx/v0_norm * l0_norm
    pellet_velphi = pellet_velphi/v0_norm * l0_norm
    pellet_velz = pellet_velz/v0_norm * l0_norm

  end subroutine pellet_init
  
  vectype elemental function pellet_deposition(r, phi, z, pres, den, rate)
    use math
    use basic
    use diagnostics
    implicit none
    real, intent(in) :: r, phi, z, pres, den, rate
    real             :: rate_norm
    select case(abs(ipellet))

    ! gaussian pellet source
    case(1)

#ifdef USE3D
!!$       pellet_deposition = pellet_rate/(sqrt(2.*pi)*pellet_var)**3 & 
!!$            *exp(-(r**2 + pellet_x**2 -2.*r*pellet_x*cos(phi-pellet_phi) &
!!$            + (z - pellet_z)**2) / (2.*pellet_var**2))
       pellet_deposition = pellet_rate/ &
            (sqrt(2.*pi)**3*pellet_var**2*pellet_var_tor) & 
            *exp(-((r-pellet_x)**2 + (z-pellet_z)**2) &
                  /(2.*pellet_var**2) &
                 -2.*r*pellet_x*(1.-cos(phi-pellet_phi)) &
                  /(2.*pellet_var_tor**2))

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
            *exp(-(r**2 + pellet_x**2 -2.*r*pellet_x*cos(phi-pellet_phi) &
            + (z - pellet_z)**2) / (2.*pellet_var**2))
#else
       pellet_deposition = pres*pellet_rate/(r*2.*pi*pellet_var**2) & 
            *exp(-((r - pellet_x)**2 + (z - pellet_z)**2) &
            /(2.*pellet_var**2))
#endif

    case(4)

       select case(ipellet_abl)
       case(1)
          rate_norm = rate
       case(2)
#ifdef USE3D
          rate_norm = 6.022e23*rate*dt*t0_norm/(1.e6*pellet_volume*n0_norm)
#else
          rate_norm = 6.022e23*rate*dt*t0_norm/(1.e6*pellet_volume_2D*n0_norm)
#endif

       end select

#ifdef USE3D
       pellet_deposition = rate_norm/ &
            (sqrt(2.*pi)**3*(pellet_var)**2*pellet_var_tor) &
            *exp(-((r-pellet_x)**2 + (z-pellet_z)**2) &
            /(2.*(pellet_var)**2) &
            -2.*r*pellet_x*(1.-cos(phi-pellet_phi)) &
            /(2.*pellet_var_tor**2))
#else
       pellet_deposition = rate_norm/(r*2.*pi*(pellet_var)**2) &
            *exp(-((r - pellet_x)**2 + (z - pellet_z)**2) &
            /(2.*(pellet_var)**2))
#endif

    case default
       pellet_deposition = 0
    end select
  end function pellet_deposition

  subroutine pellet_advance
    use basic
    use diagnostics
    implicit none
    
    pellet_x = pellet_x + pellet_velx*dt
    pellet_z = pellet_z + pellet_velz*dt
    pellet_phi = pellet_phi + pellet_velphi*dt

    ! Pellet cloud radius which contains the same number of particles as the realistic pellet                                                                                
    if(ipellet.eq.4) then

       select case(ipellet_abl)
       case(1)
          pellet_var = cloud_pel*r_p
       case(2)
          pellet_var = cloud_pel*r_p2
       end select
    
    endif

  end subroutine pellet_advance

end module pellet

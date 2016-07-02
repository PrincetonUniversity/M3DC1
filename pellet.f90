module pellet
  implicit none

  integer :: ipellet  ! 1 = include pellet injection density source
                       ! 2 = distributed density source
                       ! 3 = Gaussian source
                       ! 4 = Gaussian source including realistic ablation model  
  real :: pellet_x    ! x coordinate of pellet
  real :: pellet_phi  ! phi coordinate of pellet
  real :: pellet_z    ! z coordinate of pellet
  real :: pellet_rate ! amplitude of pellet density source
  real :: pellet_var  ! spatial dispersion of density source 
  real :: pellet_var_tor  ! spatial dispersion of density source 

  real :: pellet_velx 
  real :: pellet_velphi
  real :: pellet_velz

  integer :: ipellet_abl, ipellet_spe

  real :: nsource_pel, temp_pel, Lor_vol, q_s, shield_p, f_b
  real :: n_g, subl, T_S, Mach, f_l, r_p, pellet_ablrate, pellet_rate1
  real :: pellet_rate2, C_abl, Xp_abl, Xn_abl, r_p2, a_Te, b_Te, c_Te, d_Te
  real :: B_Li, te_norm, pellet_volume, pellet_volume_2D, cloud_pel
  
contains

  subroutine pellet_init()
     use basic
!     use diagnostics
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
!    use diagnostics
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
          rate_norm = 6.022e23*rate*dt*t0_norm/(1.e6*n0_norm)
       end select

#ifdef USE3D
       pellet_deposition = rate_norm/ &
            (sqrt(2.*pi)**3*(pellet_var)**2*pellet_var_tor) &
            *exp(-((r-pellet_x)**2 + (z-pellet_z)**2) &
            /(2.*(pellet_var)**2) &
            -2.*r*pellet_x*(1.-cos(phi-pellet_phi)) &
            /(2.*pellet_var_tor**2))
#else
       pellet_deposition = rate_norm/sqrt(2.*pi*(pellet_var)**2) &
            *exp(-((r - pellet_x)**2 + (z - pellet_z)**2) &
            /(2.*(pellet_var)**2))
#endif

    case default
       pellet_deposition = 0
    end select
  end function pellet_deposition

  subroutine pellet_advance
    use basic
!    use diagnostics
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

  subroutine calculate_parks_model
    use basic
    use math

    implicit none
         ! First model: Parks NF 94 + Lunsford

     shield_p = 0.3
     f_b = 0.5  ! From Parks NF 94
     f_l = 0.16
     T_S = 0.14 !in eV
     subl = 1.6 !in eV/atom
     n_g = 0.534 !in g.cm-3
     Mach = 1.

     q_s = 0.5*nsource_pel*temp_pel*te_norm*sqrt(8.*temp_pel*te_norm/(PI*1.e3*m_p*me_mi))
     pellet_rate1 = 4.*PI*(l0_norm*pellet_var)**2*q_s*shield_p*f_b*0.906!/(1.e-3*n_g*(subl+10./3.*T_s))                                                                                            
     pellet_rate1 = dt*t0_norm*pellet_rate1/n0_norm

     r_p = r_p-dt*t0_norm*shield_p*f_b*f_l*1.e-6*1.e-5*sqrt(1.e-5)*q_s*0.906!/&     !(n_g*(subl+T_S*(2.5+0.833*Mach**2)))                                   

     ! Second model: Parks 2015 with multienergetic electrons
     B_Li = (1./3.)*sqrt(1./(2.*log(7.69e1*1.97836e-3*sqrt(2.*temp_pel*te_norm)*3.**(-1./3.)/(6e-1))*&
          log((2.*temp_pel*te_norm)/(3.33e1)*sqrt(exp(1.)/2.))))

     f_l = 0.2*(1.-0.0946*log((4.**1.3878+1.9155)/(1.9155)))

     Xn_abl = 8.1468e-9*(5./3.-1.)**(1./3.)*f_l**(1./3.)*6.941**(-1./3.)*(r_p2*1.e2)**(4./3.)&
          *(n0_norm*nsource_pel)**(1./3.)*(te_norm*temp_pel)**(11./6.)*B_Li**(2./3.)

     Xp_abl = (8.1468e-9/(4.*4.*atan(1.)*0.534))*(5./3.-1.)**(1./3.)*f_l**(1./3.)*6.941**(2./3.)*(r_p2*1.e2)**(-2./3.)&
          *(n0_norm*nsource_pel)**(1./3.)*(te_norm*temp_pel)**(11./6.)*B_Li**(2./3.)

     a_Te = 9.0624343*(te_norm*temp_pel/800.)**(log(9.0091993/9.0624343)/log(2.5))
     b_Te = 1.5915421*(te_norm*temp_pel/800.)**(log(1.1109535/1.5915421)/log(2.5))
     c_Te = 1.9177788*(te_norm*temp_pel/800.)**(log(1.7309993/1.9177788)/log(2.5))
     d_Te = 6.6657409*(te_norm*temp_pel/800.)**(log(4.10722254/6.6657409)/log(2.5))

     C_abl = a_Te*log(1.+b_Te*(r_p2*1.e2)**(2./3.)*(nsource_pel/0.45)**(2./3.))/&
          log(c_Te+d_Te*(r_p2*1.e2)**(2./3.)*(nsource_pel/0.45)**(2./3.))

     pellet_rate2 = C_abl*Xn_abl
     pellet_ablrate = pellet_rate2

     r_p2 = r_p2-dt*t0_norm*C_abl*Xp_abl*1.e-2

     pellet_volume = 16.*atan(1.)*(pellet_var)**2*pellet_var_tor

     pellet_volume_2D = 16.*atan(1.)*(pellet_var)**2*8.*atan(1.)*pellet_x

  end subroutine calculate_parks_model

end module pellet

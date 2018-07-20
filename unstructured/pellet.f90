module pellet
  implicit none

  integer :: ipellet   ! 1 = include pellet injection density source
                       ! 2 = distributed density source
                       ! 3 = Gaussian source
                       ! 4 = Gaussian source including realistic ablation model
  integer :: ipellet_z ! Atomic number of pellet source (0 = main ion)

  real :: pellet_r    ! x coordinate of pellet
  real :: pellet_phi  ! phi coordinate of pellet
  real :: pellet_z    ! z coordinate of pellet
  real :: pellet_rate ! amplitude of pellet density source
  real :: pellet_var  ! spatial dispersion of density source
  real :: pellet_var_tor  ! spatial dispersion of density source

  real :: pellet_velr
  real :: pellet_velphi
  real :: pellet_velz
  real :: pellet_velx, pellet_vely

  integer :: ipellet_abl
  real :: r_p
  real :: cloud_pel
  real :: pellet_mix      ! (moles D2)/(moles D2 + moles impurity)
  real :: temin_abl
  real :: pellet_rate_D2  ! rate of deuterium deposition from mixed pellets

  real :: nsource_pel, temp_pel, Lor_vol
  real :: rpdot

contains

  subroutine pellet_init()
    use basic
!    use diagnostics
    implicit none
    
    ! if we're ablating, pellet_var set by pellet & cloud size
    if(ipellet_abl.gt.0) pellet_var = cloud_pel*r_p
    
    if(pellet_var_tor.le.0) pellet_var_tor = pellet_var

    ! initialize Cartesian velocities
    pellet_velx = pellet_velr*cos(pellet_phi) - pellet_velphi*sin(pellet_phi)
    pellet_vely = pellet_velr*sin(pellet_phi) + pellet_velphi*cos(pellet_phi)
    
  end subroutine pellet_init

  vectype elemental function pellet_distribution(r, phi, z, pres, inorm)
    use math
    use basic
!    use diagnostics
    implicit none
    real, intent(in) :: r, phi, z, pres
    integer, intent(in) :: inorm

    real :: x, y, px, py
    
    select case(abs(ipellet))

#ifdef USE3D
    ! gaussian pellet source
    case(1, 4, 11)
       pellet_distribution = 1./ &
            (sqrt(2.*pi)**3*pellet_var**2*pellet_var_tor) &
            *exp(-((r-pellet_r)**2 + (z-pellet_z)**2) &
                  /(2.*pellet_var**2) &
                 -2.*r*pellet_r*(1.-cos(phi-pellet_phi)) &
                  /(2.*pellet_var_tor**2))

    !......distributed source added 11/23/2011   (scj)
    case(2)
       pellet_distribution = den0*(max(pedge,real(pres))/p0)**expn

    ! gaussian pellet source
    case(3)
       pellet_distribution = pres/(sqrt(2.*pi)*pellet_var)**3 &
            *exp(-(r**2 + pellet_r**2 -2.*r*pellet_r*cos(phi-pellet_phi) &
            + (z - pellet_z)**2) / (2.*pellet_var**2))


    ! spherical, cartesian gaussian
    case(12)
       x  = r*cos(phi)
       y  = r*sin(phi)
       px = pellet_r*cos(pellet_phi)
       py = pellet_r*sin(pellet_phi)

       pellet_distribution = 1./ &
            (sqrt(2.*pi*pellet_var)**3) &
            *exp(-((x-px)**2 + (y-py)**2 + (z-pellet_z)**2) &
                  /(2.*pellet_var**2))

    ! toroidal, axisymmetric gaussian
    case(13)
       pellet_distribution = 1./(2.*pi*pellet_var**2) &
            *exp(-((r - pellet_r)**2 + (z - pellet_z)**2) &
            /(2.*pellet_var**2))
       if(itor.eq.1) pellet_distribution = pellet_distribution / r

    
#else

    ! axisymmetric gaussian pellet source
    case(1, 11, 13)
       pellet_distribution = 1./(2.*pi*pellet_var**2) &
            *exp(-((r - pellet_r)**2 + (z - pellet_z)**2) &
            /(2.*pellet_var**2))
       if(itor.eq.1) pellet_distribution = pellet_distribution / r

    !......distributed source added 11/23/2011   (scj)
    case(2)
       pellet_distribution = den0*(max(pedge,real(pres))/p0)**expn

    ! pressure-weighted gaussian pellet source
    case(3)
       pellet_distribution = pres/(2.*pi*pellet_var**2) &
            *exp(-((r - pellet_r)**2 + (z - pellet_z)**2) &
            /(2.*pellet_var**2))
       if(itor.eq.1) pellet_distribution = pellet_distribution / r

    ! different normalization of axisymmetric gaussian
    case(4)
       pellet_distribution = 1./sqrt(2.*pi*(pellet_var)**2) &
            *exp(-((r - pellet_r)**2 + (z - pellet_z)**2) &
            /(2.*(pellet_var)**2))

    ! circular, cartesian gaussian
    case(12)
       pellet_distribution = 1./ &
            (2.*pi*pellet_var**2) &
            *exp(-((r-pellet_r)**2 + (z-pellet_z)**2) &
                  /(2.*pellet_var**2))

#endif

    case default
       pellet_distribution = 0
    end select

    ! normalize distribution to Lor_vol (for distributions with ipellet double digits)
    if(inorm.ne.0 .and. ipellet.ge.10) pellet_distribution = pellet_distribution/Lor_vol

  end function pellet_distribution

  subroutine pellet_advance
    use basic
!    use diagnostics
    implicit none
    real :: x, y

    x = pellet_r*cos(pellet_phi)
    y = pellet_r*sin(pellet_phi)

    x        = x        + pellet_velx*dt
    y        = y        + pellet_vely*dt
    pellet_z = pellet_z + pellet_velz*dt

    pellet_r   = sqrt(x**2 + y**2)
    pellet_phi = atan2(y,x)

    ! Pellet cloud radius which contains the same number of particles as the realistic pellet
    if(ipellet_abl.gt.0) then
       pellet_var = cloud_pel*r_p
    endif

  end subroutine pellet_advance

  subroutine calculate_ablation
    use basic
    use math

    implicit none

    real :: dr_p
    real :: q_s, shield_p, f_b
    integer :: z_abl
    real :: rho_z, M_z  ! density (g/cm^3) and molar weight (g/mol) for atom Z
    real :: subl, T_S, Mach, f_l
    real :: C_abl, Xp_abl, Xn_abl, a_Te, b_Te, c_Te, d_Te, B_Li
    real :: G, lambda, rho0
    real :: temin_eV
    real, parameter :: n_D2 = 0.2    ! density of solid D2
    real, parameter :: M_D2 = 4.0282 ! molar weight of D2
    real, parameter :: N_A  = 6.022140857e23  ! Avogadro's number
    real, parameter :: inv3 = 1./3.


    pellet_rate_D2 = 0. ! no mixture by default

    temin_eV = temin_abl*p0_norm/(1.6022e-12*n0_norm)
    if((r_p*l0_norm).lt.1e-8 .or. temp_pel.lt.temin_eV) then
       if((r_p*l0_norm).lt.1e-8) then
          if(myrank.eq.0 .and. iprint.ge.1) print *, "No pellet left to ablate"
          r_p = 0.
       else
          if(myrank.eq.0 .and. iprint.ge.1) print *, "Temperature too low for pellet ablation"
       end if
       pellet_rate = 0.
       pellet_rate_D2 = 0.
       rpdot = 0.
       return
    end if


    ! Define density and molar mass for various Z

    z_abl = ipellet_z
    if(z_abl.eq.0) then
       select case(ipellet_abl)
       case(1,2)
          ! Lithium by default for backward compatibility
          z_abl = 3
       case(3)
          ! Only valid for neon for now
          z_abl = 10
       case default
          ! Diatomic deuterium
          z_abl = 1
       end select
    end if

    select case(z_abl)
    case(1)
       ! Assume diatomic deuterium
       rho_z = n_D2
       M_z = M_D2
    case(3)
       ! Lithium
       rho_z = 0.534
       M_z = 6.941
    case(4)
       ! Beryllium
       rho_z = 1.85
       M_z = 9.012182
    case(6)
       ! Carbon (graphite)
       rho_z = 2.267
       M_z = 12.0107
    case(10)
       ! Neon
       rho_z = 1.444
       M_z = 20.1797
    case(18)
       ! Argon
       rho_z = 1.623
       M_z = 39.948
    case default
       if(myrank.eq.0) print *, "Cannot ablate for this ipellet_z"
       ipellet_abl = 0
       return
    end select

    select case(ipellet_abl)
    case(1)

       if(z_abl.ne.3 .and. myrank.eq.0) print *, "Warning: ipellet_abl=1 only valid for lithium"

       if(pellet_mix.gt.0) then
          if(myrank.eq.0) print*, "Warning: setting pellet_mix=0. for ipellet_abl=1"
          pellet_mix = 0.
       end if

       rho0 = rho_z

       ! First model: Parks NF 94 + Lunsford
       shield_p = 0.3
       f_b = 0.5  ! From Parks NF 94
       f_l = 0.16
       T_S = 0.14 !in eV
       subl = 1.6 !in eV/atom
       Mach = 1.

       q_s = 0.5*nsource_pel*temp_pel*sqrt(8.*temp_pel/(pi*1.e3*m_p*me_mp))
       pellet_rate = 4.*pi*(l0_norm*pellet_var)**2*q_s*shield_p*f_b*0.906!/(1.e-3*rho_z*(subl+10./3.*T_s))
       pellet_rate = t0_norm*pellet_rate/n0_norm

       rpdot = shield_p*f_b*f_l*1.e-6*1.e-5*sqrt(1.e-5)*q_s*0.906!/&     !(rho_z*(subl+T_S*(2.5+0.833*Mach**2)))

    case(2)

       if(z_abl.ne.3 .and. myrank.eq.0) print *, "Warning: ipellet_abl=2 only valid for lithium"
       if(pellet_mix.gt.0) then
          if(myrank.eq.0) print*, "Warning: setting pellet_mix=0. for ipellet_abl=2"
          pellet_mix=0.
       end if

       rho0 = rho_z

       ! Second model: Parks 2015 with multienergetic electrons
       B_Li = inv3*sqrt(1./(2.*log(7.69e1*1.97836e-3*sqrt(2.*temp_pel)*3.**(-inv3)/(6e-1))*&
            log((2.*temp_pel)/(3.33e1)*sqrt(exp(1.)/2.))))

       f_l = 0.2*(1.-0.0946*log((4.**1.3878+1.9155)/(1.9155)))

       Xn_abl = 8.1468e-9*(5.*inv3-1.)**(inv3)*f_l**inv3*M_z**(-inv3)*(r_p*l0_norm)**(4.*inv3)&
            *(n0_norm*nsource_pel)**inv3*temp_pel**(5.5*inv3)*B_Li**(2.*inv3)

       Xp_abl = Xn_abl*M_z/(4.*pi*rho_z*(r_p*l0_norm)**2)

       a_Te = 9.0624343*(1.25e-3*temp_pel)**(log(9.0091993/9.0624343)/log(2.5))
       b_Te = 1.5915421*(1.25e-3*temp_pel)**(log(1.1109535/1.5915421)/log(2.5))
       c_Te = 1.9177788*(1.25e-3*temp_pel)**(log(1.7309993/1.9177788)/log(2.5))
       d_Te = 6.6657409*(1.25e-3*temp_pel)**(log(4.10722254/6.6657409)/log(2.5))

       C_abl = a_Te*log(1.+b_Te*(r_p*l0_norm)**(2.*inv3)*(nsource_pel/0.45)**(2.*inv3))/&
            log(c_Te+d_Te*(r_p*l0_norm)**(2.*inv3)*(nsource_pel/0.45)**(2.*inv3))

       pellet_rate = N_A*C_abl*Xn_abl*t0_norm/(n0_norm*l0_norm**3)

       rpdot = C_abl*Xp_abl*1.e-2

    case(3)
       ! Parks composite neon-deuterium model from 6/20/2017 (implemented 6/11/18 BCL)
       if(z_abl.ne.10 .and. myrank.eq.0) print *, "Warning: ipellet_abl=3 only valid for neon"

       lambda = 27.0837 + tan(1.48709*pellet_mix)
       G = lambda*(temp_pel*5e-4)**(5.*inv3)*(5.*r_p*l0_norm)**(4.*inv3)*nsource_pel**(inv3)  ! g/s

       ! impurity number
       Xn_abl = (1.-pellet_mix)*G/(M_z*(1.-pellet_mix) + pellet_mix*M_D2) ! mole/s
       pellet_rate = N_A*Xn_abl*t0_norm/(n0_norm*l0_norm**3) ! particles injected

       ! D2 number
       Xn_abl = pellet_mix*G/(M_z*(1.-pellet_mix) + pellet_mix*M_D2) ! mole/s
       pellet_rate_D2 = N_A*Xn_abl*t0_norm/(n0_norm*l0_norm**3) ! particles injected

       ! pellet surface recession speed
       rho0 = ((1.-pellet_mix)*M_z + pellet_mix*M_D2)/((1.-pellet_mix)*(M_z/rho_z) + pellet_mix*(M_D2/n_D2)) ! g/cm^3
       rpdot = (G/(4.*pi*rho0*(r_p*l0_norm)**2))*(t0_norm/l0_norm)

    end select

    dr_p = dt*rpdot  ! change in pellet radius

    if(dr_p.gt.r_p) then
       ! we've ablated the whole pellet
       if(myrank.eq.0 .and. iprint.ge.1) print *, "Pellet fully ablated at radius ", r_p
       pellet_rate    = (N_A/(n0_norm*l0_norm**3*dt))*(4.*inv3*pi*(r_p*l0_norm)**3)*rho0*(1.-pellet_mix)/(M_z*(1.-pellet_mix)+M_D2*pellet_mix)
       pellet_rate_D2 = (N_A/(n0_norm*l0_norm**3*dt))*(4.*inv3*pi*(r_p*l0_norm)**3)*rho0*pellet_mix/(M_z*(1.-pellet_mix)+M_D2*pellet_mix)
       r_p = 0.0
    else
       r_p = r_p-dr_p
    end if

  end subroutine calculate_ablation

end module pellet

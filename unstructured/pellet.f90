module pellet
  implicit none

  integer :: ipellet   ! positive: pellet source at each time step
                       ! negative: pellet source as initial condition
                       ! See User Guide for description of each option
  integer :: ipellet_z ! Atomic number of pellet source (0 = main ion)

  integer :: iread_pellet  ! 1 = read pellet info from pellet.dat

  integer, parameter :: maxpellets = 128
  integer :: npellets

  real, dimension(maxpellets) :: pellet_r    ! x coordinate of pellet
  real, dimension(maxpellets) :: pellet_phi  ! phi coordinate of pellet
  real, dimension(maxpellets) :: pellet_z    ! z coordinate of pellet
  real, dimension(maxpellets) :: pellet_rate ! amplitude of pellet density source
  real, dimension(maxpellets) :: pellet_var  ! spatial dispersion of density source
  real, dimension(maxpellets) :: pellet_var_tor  ! spatial dispersion of density source

  real, dimension(maxpellets) :: pellet_velr
  real, dimension(maxpellets) :: pellet_velphi
  real, dimension(maxpellets) :: pellet_velz
  real, dimension(maxpellets) :: pellet_vx, pellet_vy

  integer :: ipellet_abl
  real, dimension(maxpellets) :: r_p
  real, dimension(maxpellets) :: cloud_pel
  real, dimension(maxpellets) :: pellet_mix      ! (moles D2)/(moles D2 + moles impurity)
  real :: temin_abl
  real, dimension(maxpellets) :: pellet_rate_D2  ! rate of deuterium deposition from mixed pellets

  real, dimension(maxpellets) :: nsource_pel, temp_pel, Lor_vol
  real, dimension(maxpellets) :: rpdot

  real :: pellet_r_scl, pellet_phi_scl, pellet_z_scl
  real :: pellet_rate_scl, pellet_var_scl, pellet_var_tor_scl
  real :: pellet_velr_scl, pellet_velphi_scl, pellet_velz_scl
  real :: r_p_scl, cloud_pel_scl, pellet_mix_scl

contains

  subroutine pellet_init()
    use basic
!    use diagnostics
    implicit none
    character(LEN=10), parameter :: pellet_filename = 'pellet.dat'

    if(iread_pellet.eq.0) then
       ! use the scalar values
       npellets = 1
       pellet_r(1)       = pellet_r_scl
       pellet_phi(1)     = pellet_phi_scl
       pellet_z(1)       = pellet_z_scl
       pellet_rate(1)    = pellet_rate_scl
       pellet_var(1)     = pellet_var_scl
       pellet_var_tor(1) = pellet_var_tor_scl
       pellet_velr(1)    = pellet_velr_scl
       pellet_velphi(1)  = pellet_velphi_scl
       pellet_velz(1)    = pellet_velz_scl
       r_p(1)            = r_p_scl
       cloud_pel(1)      = cloud_pel_scl
       pellet_mix(1)     = pellet_mix_scl

    else
       npellets = 0
       call read_ascii_column(pellet_filename, pellet_r,       npellets, icol=1)
       call read_ascii_column(pellet_filename, pellet_phi,     npellets, icol=2)
       call read_ascii_column(pellet_filename, pellet_z,       npellets, icol=3)
       call read_ascii_column(pellet_filename, pellet_rate,    npellets, icol=4)
       call read_ascii_column(pellet_filename, pellet_var,     npellets, icol=5)
       call read_ascii_column(pellet_filename, pellet_var_tor, npellets, icol=6)
       call read_ascii_column(pellet_filename, pellet_velr,    npellets, icol=7)
       call read_ascii_column(pellet_filename, pellet_velphi,  npellets, icol=8)
       call read_ascii_column(pellet_filename, pellet_velz,    npellets, icol=9)
       call read_ascii_column(pellet_filename, r_p,            npellets, icol=10)
       call read_ascii_column(pellet_filename, cloud_pel,      npellets, icol=11)
       call read_ascii_column(pellet_filename, pellet_mix,     npellets, icol=12)
    end if

    ! if we're ablating, pellet_var set by pellet & cloud size
    if(ipellet_abl.gt.0) pellet_var = cloud_pel*r_p

    where(pellet_var_tor.le.0) pellet_var_tor = pellet_var

    ! initialize Cartesian velocities
    pellet_vx = pellet_velr*cos(pellet_phi) - pellet_velphi*sin(pellet_phi)
    pellet_vy = pellet_velr*sin(pellet_phi) + pellet_velphi*cos(pellet_phi)

  end subroutine pellet_init

  vectype elemental function pellet_distribution(ip, r, phi, z, pres, inorm)
    use math
    use basic
!    use diagnostics
    implicit none
    integer, intent(in) :: ip
    real, intent(in) :: r, phi, z, pres
    integer, intent(in) :: inorm

    real :: x, y, px, py

    select case(abs(ipellet))

#ifdef USE3D
    ! gaussian pellet source
    case(1, 4, 11)
       pellet_distribution = 1./ &
            (sqrt(2.*pi)**3*pellet_var(ip)**2*pellet_var_tor(ip)) &
            *exp(-((r-pellet_r(ip))**2 + (z-pellet_z(ip))**2) &
                  /(2.*pellet_var(ip)**2) &
                 -2.*r*pellet_r(ip)*(1.-cos(phi-pellet_phi(ip))) &
                  /(2.*pellet_var_tor(ip)**2))

    !......distributed source added 11/23/2011   (scj)
    case(2)
       pellet_distribution = den0*(max(pedge,real(pres))/p0)**expn

    ! gaussian pellet source
    case(3)
       pellet_distribution = pres/(sqrt(2.*pi)*pellet_var(ip))**3 &
            *exp(-(r**2 + pellet_r(ip)**2 &
            - 2.*r*pellet_r(ip)*cos(phi-pellet_phi(ip)) &
            + (z - pellet_z(ip))**2) / (2.*pellet_var(ip)**2))


    ! spherical, cartesian gaussian
    case(12)
       x  = r*cos(phi)
       y  = r*sin(phi)
       px = pellet_r(ip)*cos(pellet_phi(ip))
       py = pellet_r(ip)*sin(pellet_phi(ip))

       pellet_distribution = 1./ &
            (sqrt(2.*pi*pellet_var(ip))**3) &
            *exp(-((x-px)**2 + (y-py)**2 + (z-pellet_z(ip))**2) &
                  /(2.*pellet_var(ip)**2))

    ! toroidal, axisymmetric gaussian
    case(13)
       pellet_distribution = 1./(2.*pi*pellet_var(ip)**2) &
            *exp(-((r - pellet_r(ip))**2 + (z - pellet_z(ip))**2) &
            /(2.*pellet_var(ip)**2))
       if(itor.eq.1) pellet_distribution = pellet_distribution / r


#else

    ! axisymmetric gaussian pellet source
    case(1, 11, 13)
       pellet_distribution = 1./(2.*pi*pellet_var(ip)**2) &
            *exp(-((r - pellet_r(ip))**2 + (z - pellet_z(ip))**2) &
            /(2.*pellet_var(ip)**2))
       if(itor.eq.1) pellet_distribution = pellet_distribution / r

    !......distributed source added 11/23/2011   (scj)
    case(2)
       pellet_distribution = den0*(max(pedge,real(pres))/p0)**expn

    ! pressure-weighted gaussian pellet source
    case(3)
       pellet_distribution = pres/(2.*pi*pellet_var(ip)**2) &
            *exp(-((r - pellet_r(ip))**2 + (z - pellet_z(ip))**2) &
            /(2.*pellet_var(ip)**2))
       if(itor.eq.1) pellet_distribution = pellet_distribution / r

    ! different normalization of axisymmetric gaussian
    case(4)
       pellet_distribution = 1./sqrt(2.*pi*(pellet_var(ip))**2) &
            *exp(-((r - pellet_r(ip))**2 + (z - pellet_z(ip))**2) &
            /(2.*(pellet_var(ip))**2))

    ! circular, cartesian gaussian
    case(12)
       pellet_distribution = 1./ &
            (2.*pi*pellet_var(ip)**2) &
            *exp(-((r-pellet_r(ip))**2 + (z-pellet_z(ip))**2) &
                  /(2.*pellet_var(ip)**2))

#endif

    case default
       pellet_distribution = 0
    end select

    ! normalize distribution to Lor_vol (for distributions with ipellet double digits)
    if(inorm.ne.0 .and. ipellet.ge.10) pellet_distribution = pellet_distribution/Lor_vol(ip)

  end function pellet_distribution

  subroutine pellet_advance
    use basic
!    use diagnostics
    implicit none
    real, dimension(maxpellets) :: x, y

    where((pellet_vx**2  + pellet_vy**2 + pellet_vz**2).gt.0.)
       x = pellet_r*cos(pellet_phi)
       y = pellet_r*sin(pellet_phi)

       x        = x        + pellet_vx*dt
       y        = y        + pellet_vy*dt
       pellet_z = pellet_z + pellet_velz*dt

       pellet_r   = sqrt(x**2 + y**2)
       pellet_phi = atan2(y,x)
    end where

    ! Pellet cloud radius which contains the same number of particles as the realistic pellet
    if(ipellet_abl.gt.0) then
       pellet_var = cloud_pel*r_p
       where(pellet_var.lt.1e-8) pellet_var = 1e-8
    endif

 end do

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

    do ip=1, npellets
       pellet_rate_D2(ip) = 0. ! no mixture by default

       temin_eV = temin_abl*p0_norm/(1.6022e-12*n0_norm)
       if((r_p(ip)*l0_norm).lt.1e-8 .or. temp_pel(ip).lt.temin_eV) then
          if((r_p(ip)*l0_norm).lt.1e-8) then
             if(myrank.eq.0 .and. iprint.ge.1) print *, "No pellet left to ablate"
             r_p(ip) = 0.
          else
             if(myrank.eq.0 .and. iprint.ge.1) print *, "Temperature too low for pellet ablation"
          end if
          pellet_rate(ip) = 0.
          pellet_rate_D2(ip) = 0.
          rpdot(ip) = 0.
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

          if(pellet_mix(ip).gt.0) then
             if(myrank.eq.0) print*, "Warning: setting pellet_mix=0. for ipellet_abl=1"
             pellet_mix(ip) = 0.
          end if

          rho0 = rho_z

          ! First model: Parks NF 94 + Lunsford
          shield_p = 0.3
          f_b = 0.5  ! From Parks NF 94
          f_l = 0.16
          T_S = 0.14 !in eV
          subl = 1.6 !in eV/atom
          Mach = 1.

          q_s = 0.5*nsource_pel(ip)*temp_pel(ip)*sqrt(8.*temp_pel(ip)/(pi*1.e3*m_p*me_mp))
          pellet_rate(ip) = 4.*pi*(l0_norm*pellet_var(ip))**2*q_s*shield_p*f_b*0.906!/(1.e-3*rho_z*(subl+10./3.*T_s))
          pellet_rate(ip) = t0_norm*pellet_rate(ip)/n0_norm

          rpdot(ip) = shield_p*f_b*f_l*1.e-6*1.e-5*sqrt(1.e-5)*q_s*0.906!/&     !(rho_z*(subl+T_S*(2.5+0.833*Mach**2)))

       case(2)

          if(z_abl.ne.3 .and. myrank.eq.0) print *, "Warning: ipellet_abl=2 only valid for lithium"
          if(pellet_mix(ip).gt.0) then
             if(myrank.eq.0) print*, "Warning: setting pellet_mix=0. for ipellet_abl=2"
             pellet_mix(ip) = 0.
          end if

          rho0 = rho_z

          ! Second model: Parks 2015 with multienergetic electrons
          B_Li = inv3*sqrt(1./(2.*log(7.69e1*1.97836e-3*sqrt(2.*temp_pel(ip))*3.**(-inv3)/(6e-1))*&
               log((2.*temp_pel(ip))/(3.33e1)*sqrt(exp(1.)/2.))))

          f_l = 0.2*(1.-0.0946*log((4.**1.3878+1.9155)/(1.9155)))

          Xn_abl = 8.1468e-9*(5.*inv3-1.)**(inv3)*f_l**inv3*M_z**(-inv3)*(r_p(ip)*l0_norm)**(4.*inv3)&
               *(n0_norm*nsource_pel(ip))**inv3*temp_pel(ip)**(5.5*inv3)*B_Li**(2.*inv3)

          Xp_abl = Xn_abl*M_z/(4.*pi*rho_z*(r_p(ip)*l0_norm)**2)

          a_Te = 9.0624343*(1.25e-3*temp_pel(ip))**(log(9.0091993/9.0624343)/log(2.5))
          b_Te = 1.5915421*(1.25e-3*temp_pel(ip))**(log(1.1109535/1.5915421)/log(2.5))
          c_Te = 1.9177788*(1.25e-3*temp_pel(ip))**(log(1.7309993/1.9177788)/log(2.5))
          d_Te = 6.6657409*(1.25e-3*temp_pel(ip))**(log(4.10722254/6.6657409)/log(2.5))

          C_abl = a_Te*log(1.+b_Te*(r_p(ip)*l0_norm)**(2.*inv3)*(nsource_pel(ip)/0.45)**(2.*inv3))/&
               log(c_Te+d_Te*(r_p(ip)*l0_norm)**(2.*inv3)*(nsource_pel(ip)/0.45)**(2.*inv3))

          pellet_rate(ip) = N_A*C_abl*Xn_abl*t0_norm/(n0_norm*l0_norm**3)

          rpdot(ip) = C_abl*Xp_abl*1.e-2

       case(3)
          ! Parks composite neon-deuterium model from 6/20/2017 (implemented 6/11/18 BCL)
          if(z_abl.ne.10 .and. myrank.eq.0) print *, "Warning: ipellet_abl=3 only valid for neon"

          lambda = 27.0837 + tan(1.48709*pellet_mix(ip))
          G = lambda*(temp_pel(ip)*5e-4)**(5.*inv3)*(5.*r_p(ip)*l0_norm)**(4.*inv3)*nsource_pel(ip)**(inv3)  ! g/s

          ! impurity number
          Xn_abl = (1.-pellet_mix(ip))*G/(M_z*(1.-pellet_mix(ip)) + pellet_mix(ip)*M_D2) ! mole/s
          pellet_rate(ip) = N_A*Xn_abl*t0_norm/(n0_norm*l0_norm**3) ! particles injected

          ! D2 number
          Xn_abl = pellet_mix(ip)*G/(M_z*(1.-pellet_mix(ip)) + pellet_mix(ip)*M_D2) ! mole/s
          pellet_rate_D2(ip) = N_A*Xn_abl*t0_norm/(n0_norm*l0_norm**3) ! particles injected

          ! pellet surface recession speed
          rho0 = ((1.-pellet_mix(ip))*M_z + pellet_mix(ip)*M_D2)/((1.-pellet_mix(ip))*(M_z/rho_z) + pellet_mix(ip)*(M_D2/n_D2)) ! g/cm^3
          rpdot(ip) = (G/(4.*pi*rho0*(r_p(ip)*l0_norm)**2))*(t0_norm/l0_norm)

       end select

       dr_p = dt*rpdot(ip)  ! change in pellet radius

       if(dr_p.gt.r_p(ip)) then
          ! we've ablated the whole pellet
          if(myrank.eq.0 .and. iprint.ge.1) print *, "Pellet fully ablated at radius ", r_p(ip)
          pellet_rate(ip)    = (N_A/(n0_norm*l0_norm**3*dt))*(4.*inv3*pi*(r_p(ip)*l0_norm)**3)*rho0*(1.-pellet_mix(ip))/&
                               (M_z*(1.-pellet_mix(ip))+M_D2*pellet_mix(ip))
          pellet_rate_D2(ip) = (N_A/(n0_norm*l0_norm**3*dt))*(4.*inv3*pi*(r_p(ip)*l0_norm)**3)*rho0*pellet_mix(ip)/&
                               (M_z*(1.-pellet_mix(ip))+M_D2*pellet_mix(ip))
          r_p(ip) = 0.0
       else
          r_p(ip) = r_p(ip) - dr_p
       end if
    end do

  end subroutine calculate_ablation

end module pellet

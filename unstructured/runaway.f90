module runaway_mod

  use basic
  use arrays
  use m3dc1_nint
  use field

  implicit none

  type(field_type), private :: dnre_field1,dndt_field,f_field

  type(field_type), private :: dnre_field2,jre_field,ere_field

  type(field_type), private :: depar_field,ecrit_field
    
  real, private, parameter :: eps0 = 8.854187817D-12 ![F/m]
  real, private, parameter :: c = 299792458.D0 ![m/s]
  real, private, parameter :: ec = 1.60217662D-19 ![C]
  real, private, parameter :: hbar = 1.05457180D-34
  real, private, parameter :: me = 9.10938356D-31 ![kg]
  real, private, parameter :: va = 1542350 ![m/s]

  real, dimension(MAX_PTS) :: re_j79

contains

  subroutine runaway_deallocate()
    implicit none
    if(irunaway.eq.0) return

    call destroy_field(dnre_field1)
    call destroy_field(dnre_field2)
    call destroy_field(depar_field)
    call destroy_field(ere_field)
    call destroy_field(ecrit_field)
    call destroy_field(dndt_field)
    call destroy_field(jre_field)
    call destroy_field(f_field)
  end subroutine runaway_deallocate
    
  subroutine runaway_init()
    use basic
    use math
    
    implicit none
    if(irunaway.eq.0) return

    print *, 'Estimated Ecrit for runaways = ', &
         ec**3*(n0_norm*1e6)*17./(4.*pi*eps0**2*me*c**2), ' V/m'

    call create_field(dnre_field1)
    call create_field(dnre_field2)
    call create_field(depar_field)
    call create_field(ere_field)
    call create_field(ecrit_field)
    call create_field(dndt_field)
    call create_field(jre_field)
    call create_field(f_field)

    dnre_field1 = 0.
    dnre_field2 = 0.
    depar_field = 0.
    ere_field = 0.
    ecrit_field = 0.
    dndt_field = 0.
    jre_field = 0.
    f_field = 0.
  end subroutine runaway_init
  
  elemental subroutine runaway_current(nre,epar,Temp,Dens,&
                                       Zeff,eta,Ecrit,re_79,&
                                       re_j79,re_epar,dndt,&
                                       mr,bz,bi,ri,z,&
                                       Dens_ion,Dens_imp)
    use math
    use basic
    use kprad
	use kprad_m3dc1

    implicit none

    real, intent(in) :: nre  ! [1/m^3]
    real, intent(in) :: epar,eta
    real, intent(in) :: Temp ! [eV]
    real, intent(in) :: Dens ! Total Elec. Density [1/m^3]
    real, intent(in) :: Dens_ion ! Ion Density [1/m^3]
    real, intent(in) :: Dens_imp ! Impurity Density [1/m^3]
    real, intent(in) :: Zeff ! [1]
    real, intent(in) :: bz
    real, intent(in) :: ri,z
    real, intent(in) :: bi
    real, intent(out) :: re_79,dndt,re_j79,re_epar,Ecrit
    real :: Clog,x,nu,vth,esign,teval,jpar,nrel,a,r,Ed, &
            f,dt_si,tmp,nretmp,Dens1,sd,sa,nra,gamma,tau
    real :: sbeta, scomp ! RiD: Terms for Tritium and Compton sources
    integer, intent(in) :: mr
    integer :: l, nl
    
    ! RiD: Additional variables for partially scrrned Dreicer term
    real :: ZImp, ZmaxImp
    
    
    dndt = 0.
    sd = 0. ! Dreicer
    sa = 0. ! Avalanche
    sbeta = 0. ! Tritium Beta
	scomp = 0. ! Gamma Compton source
	
    
    
    dt_si = dt * t0_norm
    nrel = nre
    nl = 10
    tmp = 0.
    nretmp = 0.
    Dens1 = Dens 
    
    
    if(mr.ne.0) then 
      re_79 = 0.D0
      re_j79 = 0.D0 
      return
    endif
             
!   based on formula in Stahl, et al, PRL 114 115002 (2015)
    teval = max(1.,Temp)   ! note: this sets a minimum temperature of 1 eV for runaway production
!   Note that nre is the density of RE flowing in the direction of B field,
!   thus their growth is linked to -epar (due to the negative charge)
    esign = sign(1., -epar)
       jpar = nre!*ec*cre*va
       nra = nre/ec/cre/va 
       re_epar = epar! - 1.0*abs(eta*jpar/b1/bz/ri)
       if (Temp.ge.teval .and. mr.eq.0) then
          Clog = 14.78D0-0.5*log(Dens1/1.d20)+log(Temp/1.d3)
          Ecrit = ec**3*Dens1*Clog/(4*pi*eps0**2*me*c**2) ! RiD: Connor Hastie Field 
          vth = sqrt(2*ec*Temp/me) ! Thermal Velocity
          nu = Dens1*ec**4*Clog/(4*pi*eps0**2*me**2*vth**3)
          tau = me*c/ec/Ecrit
          a = sqrt((1/ri-xmag)**2+(z-zmag)**2)
          r = sqrt(1/ri**2+z**2)
          gamma = 1/(1+1.46*sqrt(a/r)+1.72*a/r) 
          x = (abs(re_epar)*ec*Temp)/(Ecrit*me*c**2) ! Epar / ED
          Ed = abs(re_epar)/Ecrit ! RiD: Normalized electric field / CH field = E_star
          if(Ed < 1) Ed = 1
          if(abs(re_epar).gt.Ecrit) then
			! Dreicer Source 1 = Classical, 2 = Parial Screening, 0 = off
			  if (iDreicer.eq.1) then ! Classical Dreicer
				sd = Dens1*nu*x**(-3.D0*(1.D0+Zeff)/1.6D1) &
                 *exp(-1.D0/(4*x)-sqrt((1.D0+Zeff)/x)) ! Dreicer Source [per cubic m per s]
              endif
              if (iDreicer.eq.2) then ! Partially Screened Dreicer
				! For Partially Screened Dreicer
				ZImp = 1.0 ! Impurity Charge [-]
				ZmaxImp = 1.0 ! Atomic no. of Impurity [-]
				! Update nI, neI, ZI and ZmaxI if using KPRAD
				IF(ikprad.ne.0) then
							ZmaxImp = 1.0 *  kprad_z ! Impurity atomic no.
							Zimp = (Dens1 - 1.0*Dens_ion)/Dens_imp ! Impurity Ionization
				END IF
                                IF(x>0.13) THEN ! Upper bound of NN 
                                        x = 0.13
                                END IF
				sd = getDreicerHesslow(Dens_ion,Dens_imp,ZmaxImp,ZImp,x,Temp) ! [per cubic m per s]
              endif
              ! Tritium Beta Source
              if (iTritBeta.eq.1) then ! Tritium Beta Calculation
				sbeta = 0.5 * Dens_ion * beta_source(Ed) ! Tritium beta source [per cubic m per s]
              endif
              ! Compton Source
              if (iCompton.eq.1) then ! COmpton Source for SPARC
				scomp = Dens1 * compton_source(Ed) ! Compton Source [per cubic m per s]
				if (Temp < 2e3) then ! Non-Prompt Compton
					scomp = 1e-3 * scomp ! RiD: Reduce compton contribution when temp. < 2 keV 
				endif
			  endif
              ! avalanche growth only when epar and nre have opposite sign
              if (re_epar*nre<0) then
				if (iAvalanche.eq.1) then ! Avalanche Term Calculation
					sa = nra/tau/Clog*sqrt(pi*gamma/3/(Zeff+5))*(Ed-1)* &
                    1/sqrt(1-1/Ed+4*pi*(Zeff+1)**2/3/gamma/(Zeff+5)&
                    /(Ed**2+4/gamma**2-1)) 
                endif
              else
                 sa = 0.
              endif
              dndt = (1.0*sd*esign + &
					1.0*sbeta*esign + &
					1.0*scomp*esign + sa)*cre*ec*va ! Combining all sources
                 nrel = nre + dndt*dt_si
          else
              nrel = nre 
              dndt = 0.
          end if

       else
           nrel = nre
       endif
    
       re_79 = nrel
       re_j79 = nre
  
       if (abs(re_j79) .ge. abs(ri*bz)) then
            dndt = 0.
       endif

  end subroutine runaway_current

  subroutine eval_runaway(itri,izone)
    use basic
    use m3dc1_nint
    use electric_field
    use diagnostics
    use math
    use kprad
	use kprad_m3dc1
    
    implicit none

    integer, intent(in) :: itri,izone
    vectype, dimension(MAX_PTS) :: epar, te, ne, nre, eta, &
									n_ion, kp_den, kp_z
    vectype, dimension(MAX_PTS) :: dndt, re_79, re_j79, &
                                   re_epar, ecrit, bz, bi  
    vectype, dimension(dofs_per_element) :: dofs
    integer, dimension(MAX_PTS) :: mr, tmp
    
    ! For Partially Screened Dreicer
    ! vectype, dimension(MAX_PTS) :: nImp, neImp, ZImp


    if(irunaway.eq.0 .or. izone.ne.1) return
#ifdef USECOMPLEX
    re_79 = 0.
    re_j79 = 0.
    return
#else
    call electric_field_par(0, epar, izone)

    dofs = intx2(mu79(:,:,OP_1),epar)
    call vector_insert_block(ere_field%vec,itri,1,dofs,VEC_ADD)
    
    ! convert to SI units
    epar = epar*e0_norm*c*1e-4 ! Parallel E-field [V/m]
    te = tet79(:,OP_1)*(p0_norm/n0_norm)/1.6022e-12 ! Elec. Temp [eV]
    ne = net79(:,OP_1)*n0_norm*1e6 ! tot. electron density [per cubic m]
    nre = nre179(:,OP_1)*j0_norm/(c*1e-3)!n0_norm*1e6
    eta = eta79(:,OP_1)*e0_norm/j0_norm*c**2*1e-7 ! Resitivity [Ohm m]
    call magnetic_region(pst79(:,OP_1),&
         pst79(:,OP_DR),pst79(:,OP_DZ),x_79,Z_79,mr)
    bz = pst79(:,OP_GS)/((c*1e-3)/j0_norm)
    bi = bi79(:,OP_1)
    
    ! RiD: Adding quantities required for new sources
    n_ion = nt79(:,OP_1)*n0_norm*1e6 ! Ion density [per cubic m]
    kp_den = 0.0 * n_ion ! Impurity density [per cubic m]
    IF(ikprad.ne.0) THEN 
		kp_den = kprad_fz*nt79(:,OP_1)*n0_norm*1e6
	END IF
    
    ! Call RE subroutine
    call runaway_current(nre,epar,te,ne,z_ion,eta,ecrit,&
                         re_79,re_j79,re_epar,dndt,mr,bz,bi,ri_79,z_79,&
                         n_ion,kp_den) ! RiD: added ion density and impurity density

    ! convert back to normalized units
    dndt = dndt*t0_norm/(j0_norm/c/1e-3)

    re_79 = re_79*(c*1e-3)/j0_norm!/(n0_norm*1e6)

    re_j79 = re_j79*(c*1e-3)/j0_norm
   
    re_j79 = re_j79+dndt*dt

    re_epar = re_epar/(e0_norm*c*1e-4)

    ecrit = ecrit/(e0_norm*c*1e-4)
    
    dofs = intx2(mu79(:,:,OP_1),dndt) 
    call vector_insert_block(dndt_field%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),re_79)
    call vector_insert_block(dnre_field1%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),re_epar)
    call vector_insert_block(depar_field%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),ecrit)
    call vector_insert_block(ecrit_field%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),re_j79)
    call vector_insert_block(jre_field%vec,itri,1,dofs,VEC_ADD)



#endif

  end subroutine eval_runaway

  subroutine runaway_advance
    use basic
    use newvar_mod
    use m3dc1_nint
    implicit none

    integer :: ier

    if(irunaway.eq.0) return

    call newvar_solve(jre_field%vec, mass_mat_lhs)

    nre_field(1) = jre_field
    dnre_field2 = 0.
    dnre_field1 = 0.
    jre_field = 0.

  end subroutine runaway_advance

  subroutine smooth_runaway
    use basic
    use m3dc1_nint
    use newvar_mod
    use diagnostics
    use math
    implicit none

    vectype, dimension(MAX_PTS) :: jre
    vectype, dimension(MAX_PTS, OP_NUM) :: jre79
    vectype, dimension(dofs_per_element) :: dofs
    integer :: numelms, itri
    numelms = local_elements()    
    jre_field = 0.     
    do itri=1, numelms
       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)
       call eval_ops(itri, nre_field(1), jre79,rfac)

       where(real(jre79(:,OP_1)) .ge. 0.0)
               jre79(:,OP_1) = 0.0
       endwhere

       jre = jre79(:,OP_1)
       dofs = intx2(mu79(:,:,OP_1),jre)
       call vector_insert_block(jre_field%vec,itri,1,dofs,VEC_ADD)
    enddo

    call newvar_solve(jre_field%vec, mass_mat_lhs)
    nre_field(1) = jre_field
    jre_field= 0.




  end subroutine
  
  ! RiD: Adding Functions Required for Compton and Beta Tritium Sources
  pure function get_Wc(E_star) result(outval)
          ! Returns the critical RE energy in [eV]
          ! E_star = Par. E-field normalized with CH electric field
          real, intent(in) :: E_star
          real :: outval
          real alpha
          IF (E_star <= 1.0) THEN
                  alpha = 1.0 + 1.0e-6 ! Set to min. value for calculation
          ELSE
                  alpha = 1.0 * E_star
          ENDIF
          outval = 9.1093837015e-31 * (3.0e8)**2 * ((1-alpha**(-1))**(-0.5)-1) / 1.602e-19 ! Critical energy in [eV]
	end

	pure function beta_source(E_star) result(outval)
			! Returns the beta source rate in per second
			! Multiply with tritium denisty to get generation rate in per
			! cubic m
			! per second
			real, intent(in) :: E_star  !E-field normalized by Connor-hastie field
			real :: outval
			real A, mu, sig ! Fit parameters
			real Wc ! Critical Runaway energy in keV
			Wc = get_Wc(E_star) * 1.e-3 ! keV
			A =2.09023917e-09
			mu = -3.80424810
			sig = 6.81315848
			
			IF (Wc > 18.6) THEN
					outval = 0.0
			ELSEIF (Wc < 0.1) THEN
					outval = 1.8e-9 ! Maximum Rate
			ELSE
					outval = A * exp(-(Wc-mu)**2/(2*sig**2)) 
			ENDIF

	end

	pure function compton_source(E_star) result(outval)
			! Returns compton rate in per second
			! Multiply with total electron denisty to get rate in per cubic
			! m per
			! second
			real, intent(in) :: E_star ! E-field normalized by Connor-hastie field
			real :: outval
			real Wc, x
			real A, x0, L, c, B ! Fitting Parameters
			
			Wc = get_Wc(E_star) * 1.0e-3 ! Critical Runaway energy in [keV]
			x = log(Wc)
			IF (Wc < 0.1) THEN
					outval = 4.3e-11 ! MAX rate for SPARC
			ELSEIF (Wc > 10**4.5) THEN
					outval = 0.0
			ELSE
					! Fit parameters for SPARC
					A = 6.00578227e-01
					L = 1.25488546e+00
					x0 = 4.44639341e+00
					c = 2.12986727e-11
					B = -2.17893889e-11
					outval = B * tanh(A *(x-x0)/L) + c
			ENDIF
	end
	
	! RiD: Adding functions required for caluclation of partially screened
	! Dreicer term
	
	pure function getDreicerHesslow(nD,nI,ZmaxI,ZI,EoED,T) result(outval)
        implicit none
        include  'NN_params.inc'
        real, intent(in) :: nD, nI, ZI, EoED, T, ZmaxI
        real :: outval
        ! Returns the partially screened Dreicer source term [m^-3/s]
        ! based on original MATLAB implementation in
        !https://github.com/unnerfelt/dreicer-nn
    
        ! nD = Diuterium denisty [per cubic m]
        ! nI = Impurity Denisty  [per cubic m]
        ! ZmaxI = Atomic no. of Impurity [-]
        ! ZI = Impurity Charge [-]
        ! EoED = Normalized Electric Field [-]
        ! T = Temp. [eV]
        real :: nfree, ntot, Zeff, Zeff0, Z0oZ, ZZ0, lnnfree, nfreeontot,tau

        ! Derived Parameters
        nfree = nD * 1 + nI * ZI
        ntot = nD * 1 + nI * ZmaxI ! total electron density
        Zeff = (nD * 1**2 + nI * ZI**2) / nfree ! effective ioniztaion
        Zeff0 = (nD * 0 + nI * (ZmaxI**2-ZI**2)) / ntot
        Z0oZ = (nD * 1 + nI * ZI / ZmaxI) / ntot
        ZZ0 = (nD * 1 + nI * ZI * ZmaxI) / ntot
        lnnfree = log(nfree)
        nfreeontot = nfree / ntot

       ! INPUT
        X = (/ Zeff, &
               Zeff0, &
               Z0oZ, &
               ZZ0, &
               lnnfree,&
               nfreeontot,&
               EoED,&
               log(T/510998.946)/)


        ! Normalize your Inputs
        X = (X - input_mean) / input_std

        ! Neural Net Operations

        Y = dot_product(W5, &
             tanh(matmul(W4, &
                 tanh(matmul(W3, &
                     tanh(matmul(W2, &
                        tanh(matmul(W1, X) + b1)) &
                     + b2)) &
                + b3)) &
             + b4)) &
        + b5

        ! De-Normalize your output
        Y = exp(Y * output_std + output_mean)
        Y = Y * 4 / (3 * sqrt(3.14159265))
        tau = get_tau_ee(nfree,T) ! El-El collison time [s]
        outval = nfree * Y / tau ! [per cubic m per s]
        
	end function getDreicerHesslow  

	pure function get_tau_ee(density, temp_e) result(outval)
	  implicit none
	  real, intent(in) :: density, temp_e
	  double precision  :: e, eps0, me, lnLambda
	  double precision :: num, den
	  real :: outval
	  ! Constants in double precision
	  e = 1.602176634D-19      ! Elementary charge [C]
	  me = 9.1093837015D-31    ! Electron mass [kg]
	  eps0 = 8.8541878128D-12  ! Vacuum permittivity [F/m]

	  ! Coulomb logarithm
	  lnLambda = 14.9d0 - 0.5d0 * log(density / 1.0d20) + log(temp_e / 1.0d3)
	  ! Calculate tau_ee
	  num = 8D0*sqrt(2.D0)*3.14159265D0*(e*temp_e)**1.5*eps0**2*sqrt(me)
	  den = lnLambda* e**4 * density   
	  outval = num/den
	end function get_tau_ee



     
end module runaway_mod

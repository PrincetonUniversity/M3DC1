module runaway_mod

  use basic
  use arrays
  use m3dc1_nint
  use field
  use kprad
  use kprad_m3dc1
  use auxiliary_fields

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
  real, private, parameter :: pi = 3.14159265358979323846d0 ![-]
  integer, parameter :: dp = kind(1.0d0)

  real, dimension(MAX_PTS) :: re_j79
  
  
  ! CONSTANTS FOR AVALNCHING CALUCLATION (RiD) !
  
  ! --- NEON DATA (Size 10) ---
  real(dp), parameter, dimension(10) :: Ne_ln_aBar_const = (/ &
       4.7064_dp, 4.6097_dp, 4.5028_dp, 4.3858_dp, 4.2638_dp, &
       4.1266_dp, 3.9568_dp, 3.6810_dp, 3.1777_dp, 3.1278_dp /)
       
  real(dp), parameter, dimension(10) :: Ne_ln_I_const = (/ &
      -8.2227_dp, -8.0370_dp, -7.8614_dp, -7.6837_dp, -7.4994_dp, &
      -7.2788_dp, -6.9808_dp, -6.5976_dp, -5.8933_dp, -5.8320_dp /)

  ! --- ARGON DATA (Size 18) ---
  real(dp), parameter, dimension(18) :: Ar_ln_aBar_const = (/ &
       4.5677_dp, 4.5007_dp, 4.4289_dp, 4.3513_dp, 4.2698_dp, &
       4.1815_dp, 4.0829_dp, 3.9748_dp, 3.8511_dp, 3.7882_dp, &
       3.7203_dp, 3.6436_dp, 3.5579_dp, 3.4525_dp, 3.3093_dp, &
       3.0503_dp, 2.5578_dp, 2.5332_dp /)

  real(dp), parameter, dimension(18) :: Ar_ln_I_const = (/ &
      -7.9050_dp, -7.7532_dp, -7.6076_dp, -7.4626_dp, -7.3178_dp, &
      -7.1665_dp, -7.0055_dp, -6.8020_dp, -6.5538_dp, -6.4647_dp, &
      -6.3644_dp, -6.2465_dp, -6.1070_dp, -5.9219_dp, -5.6535_dp, &
      -5.3213_dp, -4.6937_dp, -4.6598_dp /)

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
  
  impure elemental subroutine runaway_current(nre,epar,Temp,Dens,&
                                       Zeff,eta,Ecrit,re_79,&
                                       re_j79,re_epar,dndt,&
                                       mr,bz,bi,ri,z,&
                                       Dens_ion,Dens_imp,btoroidal)
    use math
    use basic
    use kprad
    use kprad_m3dc1
    use auxiliary_fields

    implicit none

    real, intent(in) :: nre  ! [1/m^3]
    real, intent(in) :: epar,eta
    real, intent(in) :: Temp ! [eV]
    real, intent(in) :: Dens ! Total Elec. Density [1/m^3]
    real, intent(in) :: Dens_ion ! Ion Density [1/m^3]
    real, intent(in) :: Dens_imp ! Impurity Density [1/m^3]
    real, intent(in) :: btoroidal ! Toroidal Field * R [Tesla - m]
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
    
    ! RiD: Additional variables for partially screened Dreicer term
    real :: ZImp, ZmaxImp
    ! RiD: Additional variables for partially screened avalnching model
    real :: Zval(2), Z0high(2), Z0low(2), Z0val(2), nj(2)
    real :: ne_free, Ech, Estar, lnAc, Ylow, Yhigh, weight , ntot, Zeff_val
    
    
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
          a = sqrt((1/ri-xmag)**2+(z-zmag)**2) ! RiD: ri = 1/R
          r = sqrt(1/ri**2+z**2)
          gamma = 1/(1+1.46*sqrt(a/r)+1.72*a/r) 
          x = (abs(re_epar)*ec*Temp)/(Ecrit*me*c**2) ! Epar / ED
          Ed = abs(re_epar)/Ecrit ! RiD: Normalized electric field / CH field = E_star
          if(Ed < 1) Ed = 1
          if(abs(re_epar).gt.Ecrit) then ! E-field must be greater than Critical Field for RE generation
			! ***** DREICER GENERATION ****** !
			! Dreicer Source 1 = Classical, 2 = Parial Screening, 0 = off
			  if (iDreicer.eq.1) then ! Classical Dreicer
				sd = Dens1*nu*x**(-3.D0*(1.D0+Zeff)/1.6D1) &
                 *exp(-1.D0/(4*x)-sqrt((1.D0+Zeff)/x)) ! Dreicer Source [per cubic m per s]
              endif ! End Classical Dreicer 
              if (iDreicer.eq.2) then ! Partially Screened Dreicer
                    if ((x < 0.01) .or. (Dens_ion<1.D16)) then ! Lower Limit of Neural Network is x = 0.01
                              sd = 0.0
                     else
						! For Partially Screened Dreicer
						ZImp = 1.0 ! Impurity Charge [-]
						ZmaxImp = 1.0 ! Atomic no. of Impurity [-]
						! Update nI, neI, ZI and ZmaxI if using KPRAD
						IF(ikprad.ne.0) then
									ZmaxImp = 1.0 *  real(kprad_z) ! Impurity atomic no.
									Zimp = (Dens1 - &
											1.0*Dens_ion)/Dens_imp ! Impurity Ionization
									Zimp = min(ZmaxImp,Zimp) ! Maximum Impurity Ionization = ZmaxImp
						END IF
						if (x > 0.13) then ! Bound of the Neural Network
								x = 0.13
						endif
						sd = getDreicerHesslow(Dens_ion,Dens_imp,&
												ZmaxImp,ZImp,x,Temp) ! [per cubic m per s]
!~ 							!! *** PRINT FOR DEBGGING ***** !!					
!~                         IF (myrank.eq.1) THEN ! Printing
!~ 							print *, "Dens_ion, Dens_imp [per cubic m] = ", Dens_ion, Dens_imp
!~ 							print *, "ZmaxImp, ZImp = ", ZmaxImp, ZImp
!~ 							print *, "Temp [eV] = ", Temp
!~ 							print *, "E/Ed [-] = ", x
!~ 							print *, "sd [per cubic m per s] = ", sd
!~                         ENDIF
!~ 							!! ******************* !!
					endif
              endif ! End Classical Dreicer 
              ! ***** END OF DREICER GENERATION ****** !
              
              ! ***** ACTIVATED SOURCES ****** !
              ! Tritium Beta Source - by defualt assumes 0.5*ni is Tritium
              if (iTritBeta.eq.1) then ! Tritium Beta Calculation
				sbeta = 0.5 * Dens_ion * beta_source(Ed) ! Tritium beta source [per cubic m per s]
              endif
              ! Compton Source
              if (iCompton.eq.1) then ! Compton Source for SPARC
                                if (ikprad.ne.0) then
                                       scomp = (kprad_z * Dens_imp &
                                       + Dens_ion) * compton_source(Ed) 
                                else 
				       scomp = Dens1 * compton_source(Ed) ! Compton Source [per cubic m per s]
                                endif
				if (Temp < 2e3) then ! Non-Prompt Compton
					scomp = 1e-3 * scomp ! RiD: Reduce compton contribution when temp. < 2 keV 
				endif
			  endif
			  ! ***** END OF ACTIVATED SOURCES ****** !
			  
			  ! ***** AVALANCHING MODEL  ****** !
              ! avalanche growth only when epar and nre have opposite sign
              if (re_epar*nre<0) then
				if (iAvalanche.eq.1) then ! Avalanche Term Calculation (RP Model)
					sa = nra/tau/Clog*sqrt(pi*gamma/3/(Zeff+5))*(Ed-1)* &
                    1/sqrt(1-1/Ed+4*pi*(Zeff+1)**2/3/gamma/(Zeff+5)&
                    /(Ed**2+4/gamma**2-1)) 
                else if (iAvalanche.eq.2) then ! Hesslow 2019 Model 
					Zval = (/1., real(kprad_z) /) ! Max Ionization
					Zimp = (Dens1 - Dens_ion)/Dens_imp ! Avg. impurity ionization
					Z0val = (/1.,  Zimp/)
					nj = (/Dens_ion, Dens_imp /) ! Density [per cubic m]
					ne_free = sum(Z0val * nj) ! Free elec. density [per cubic m]
					lnAc = 14.6_dp + 0.5_dp * log( Temp / (ne_free/1e20_dp)) ! Coul. Log
					Ech = ec**3 * ne_free * lnAc / (me * c**2 * 4.0_dp * pi * eps0**2) ! Connor Hastie Field [V/m]
					Estar = abs(re_epar)/Ech ! Normalized E-field [-]
					! Interpolation !
					Z0high = (/1.,  real(CEILING(Zimp))/)
					Z0low = (/1.,  real(FLOOR(Zimp))/)
					weight = Zimp - real(FLOOR(Zimp))
					
					Ylow = 0.0 ! Avalanching rates
					Yhigh = 0.0
					
					call getHesslowAvalanchingRate(Estar,Zval,Z0low,nj, &
						2,Temp,"Ne",btoroidal*ri,Ylow)
					call getHesslowAvalanchingRate(Estar,Zval,Z0high,nj, &
						2,Temp,"Ne",btoroidal*ri,Yhigh)
					sa = (1.0_dp - weight) * Ylow + weight * Yhigh ! Combine [per second]
!~ 					!! ****** RiD: For Debugging *****!!
!~ 					if (myrank.eq.8) then ! Printing
!~ 						print *, "Yava [s-1] = ", sa
!~ 						print *, "ne, ni, nz = ", Dens1, Dens_ion, Dens_imp
!~ 						print *, "Estar = ", Estar
!~ 						print *, "Temp = ", Temp
!~ 						print *, "Btor = ", btoroidal
!~ 						print *, "ri = ", ri
!~ 						print *, "z = ", z
!~                        print *, "rank = ", myrank
!~                     end if
!~                     !! ********** !!
					sa = nra * sa 
				else if (iAvalanche.eq.3) then ! Simple RP Model
					lnAc = 14.6_dp + 0.5_dp * log( Temp / (Dens1/1e20_dp)) ! Coul. Log
					Ech = ec**3 * Dens1 * lnAc / (me * c**2 * 4.0_dp * pi * eps0**2) ! Connor Hastie Field [V/m]
					Estar = abs(re_epar)/Ech ! Normalized E-field [-]
					ntot = Dens_ion * 1.0 + Dens_imp * kprad_z ! Total e density
					Zimp = (Dens1 - Dens_ion)/Dens_imp ! Avg. impurity ionization
					Zeff_val = (Dens_ion * 1**2 + Dens_imp * Zimp**2)/Dens1
					sa = 1 * (ntot/Dens1) * ec/(me*c*lnAc) * Ech/sqrt(Zeff_val+5) * (Estar-1)
!~ 					if (myrank.eq.1) then
!~ 						print *, "Yava [s-1] = ", sa
!~ 						print *, "ne, ni, nz = ", Dens1, Dens_ion, Dens_imp
!~ 						print *, "Estar = ", Estar
!~ 						print *, "Temp = ", Temp
!~ 					end if
					sa = nra * sa 
                endif  ! Ending avalanching options iAvalanche
              else
                      sa = 0.
              endif ! ending if (re_epar*nre<0)
              
              dndt = ((sd*esign + &
					sbeta*esign + &
					scomp*esign + sa)*cre*ec*va + jre_const) ! Combining all sources
					
              nrel = nre + dndt*dt_si 
              
                 
          else !  the next section executes when the e-field is less than Ecrit
			
              nrel = nre + jre_const * dt_si  ! RiD Adding jre_const [A/m2/s]
              dndt = 0. + jre_const  ! RiD Adding jre_const [A/m2/s]
          end if ! ending if(abs(re_epar).gt.Ecrit);

       else
           nrel = nre
       endif ! ending if (Temp.ge.teval .and. mr.eq.0)
    
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
    use auxiliary_fields
    
    implicit none

    integer, intent(in) :: itri,izone
    vectype, dimension(MAX_PTS) :: epar, te, ne, nre, eta, &
            n_ion, kp_den, kp_z, btoroidal
    vectype, dimension(MAX_PTS) :: dndt, re_79, re_j79, &
                                   re_epar, ecrit, bz, bi  
    vectype, dimension(dofs_per_element) :: dofs
    integer, dimension(MAX_PTS) :: mr, tmp
    


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
    btoroidal = bztx79(:,OP_1) ! Toroidal magnetic field * R
    
    ! RiD: Adding quantities required for new sources
    n_ion = nt79(:,OP_1)*n0_norm*1e6 ! Ion density [per cubic m]
    kp_den = 0.0 * n_ion ! Impurity density [per cubic m]
    IF(ikprad.ne.0) THEN
        call calculate_kprad_totden(itri, kp_den)
        kp_den = kp_den*n0_norm*1e6
	!kp_den = kprad_fz*nt79(:,OP_1)*n0_norm*1e6
	END IF
    
    ! Call RE subroutine
    call runaway_current(nre,epar,te,ne,z_ion,eta,ecrit,&
                         re_79,re_j79,re_epar,dndt,mr,bz,bi,ri_79,z_79,&
                         n_ion,kp_den,btoroidal) ! RiD: added ion density and impurity density

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
	  lnLambda = 14.9d0 - 0.5d0 * log(density / 1.0d20) &
				+ log(temp_e / 1.0d3)
	  ! Calculate tau_ee
	  num = 8D0*sqrt(2.D0)*3.14159265D0*(e*temp_e)**1.5*eps0**2*sqrt(me)
	  den = lnLambda* e**4 * density   
	  outval = num/den
	end function get_tau_ee
	
	
	!!!!!!!  Hesslow + (2019) Impurity Avalanching Model !!!!!!!!
	
	! Effective electric field !
	
	pure function calcDelta(Eceff,nuD0,nuS1,tauRadInv,phib1,phib2,pc0) result(delta)
		implicit none
		integer, parameter :: dp = kind(1.0d0)
		real(dp), intent(in) :: Eceff,nuD0,nuS1,tauRadInv,phib1,phib2,pc0
		real(dp) :: delta
		delta = (nuD0/(nuS1**2)) * (nuD0*tauRadInv/Eceff + phib1 + phib2*log(pc0))
	end function calcDelta

	pure function calcEceff(delta,nuS0,nuS1,nuD1,nuD0,pc0) result(Eceff)
		implicit none
		integer, parameter :: dp = kind(1.0d0)
		real(dp), intent(in) :: delta,nuS0,nuS1,nuD1,nuD0,pc0
		real(dp) :: Eceff
		Eceff = nuS0 + nuS1 * ( (1.0_dp + nuD1/nuD0)*log(pc0) + sqrt(2.0_dp*delta + 1.0_dp) )
	end function calcEceff
	
	pure function getEc_eff(Z,Z0,nj,nSpecies,T,B,element) result(Eceff_out)
	
	!~ # Returns the effective critical electric field in [V/m] with Ne or Ar impurities 
	!~ # INPUTS:
	!~ #       Z   : atomic number of each ion species
	!~ #       Z0  : net charge of plasma for each ion species
	!~ #       nj  : density of each ion species [m^(-3)]
	!~ #       Z, Z0 and nj must be arrays of length nSpecies
	!~ #       T   : temperature [eV]  (for calculation of the Coulomb logarithm)
	!~ #       B   : magnetic field [T]
	!~ #       element: "Ne" or "Ar"
	
		implicit none
		
		! Constants !
		integer, parameter :: dp = kind(1.0d0)
		real(dp), parameter :: e = 1.602176634d-19
		real(dp), parameter :: epsilon0 = 8.8541878128d-12
		real(dp), parameter :: me = 9.1093837015d-31
		real(dp), parameter :: c = 2.99792458d8
		real(dp), parameter :: pi = 3.14159265358979323846d0
		
		integer, intent(in) :: nSpecies ! NUmber of species
		real(dp), intent(in) :: Z(nSpecies), Z0(nSpecies), nj(nSpecies) ! 
		real(dp), intent(in) :: T, B
		character(len=*), intent(in) :: element
		real(dp) :: Eceff_out

		real(dp) :: ln_aBar(18), ln_I(18)
		real(dp) :: ne_free, ntot, Zeff, Zfulleff, lnLStar
		real(dp) :: nuD0, nuD1, nuS0, nuS1
		real(dp) :: Ne(nSpecies)
		real(dp) :: tauRadInv, ALPHA, bremsprefactor, phib1, phib2
		real(dp) :: pc0, Ectot, Eceff0, delta0, delta1, EceffOverEctot
		integer :: i
		
		! pick the element
		if (trim(element) == "Ne") then
			ln_aBar(1:10) = Ne_ln_aBar_const
			ln_I(1:10)    = Ne_ln_I_const
		else if (trim(element) == "Ar") then
			ln_aBar(1:18) = Ar_ln_aBar_const
			ln_I(1:18)    = Ar_ln_I_const
		else
			Eceff_out = 0.0_dp
			return
		end if
		
		! free and total electron densities
		Ne = Z - Z0
		ne_free = sum(Z0*nj)
		ntot    = sum(Z*nj)

		Zeff = sum((Z0**2)*nj) / ne_free
		Zfulleff = sum((Z**2)*nj) / ne_free

		lnLStar = 14.6_dp + 0.5_dp * log( T / (ne_free/1.0e20_dp) ) ! Coul. Log
		
		! compute nu constants
		nuD0 = (1.0_dp + Zeff) - (2.0_dp/3.0_dp) * sum( (Ne**2)*nj )/ne_free / lnLStar
		nuD1 = Zfulleff / lnLStar
		nuS0 = 1.0_dp
		nuS1 = (1.0_dp + 3.0_dp*sum(Ne*nj)/ne_free) / (2.0_dp*lnLStar)

		do i = 2, nSpecies
			nuD0 = nuD0 + (Z(i)**2 - Z0(i)**2) * ln_aBar(int(Z0(i))+1) * (nj(i)/ne_free) / lnLStar
			nuS0 = nuS0 + Ne(i) * (-ln_I(int(Z0(i))+1) - 1.0_dp) * (nj(i)/ne_free) / lnLStar
		end do
		
		! radiation terms
		tauRadInv = B**2 / (ne_free/1e20_dp) / (15.44_dp * lnLStar)

		ALPHA = 1.0_dp/137.036_dp
		bremsprefactor = ALPHA * Zfulleff / lnLStar
		phib1 = bremsprefactor * 0.35_dp
		phib2 = bremsprefactor * 0.20_dp

		pc0 = nuD0/(2.0_dp*nuS1)

		! iterate calculation
		Ectot = sum(Z*nj)/ne_free
		Eceff0 = Ectot
		delta0 = calcDelta(Eceff0,nuD0,nuS1,tauRadInv,phib1,phib2,pc0)
		delta1 = calcDelta(calcEceff(delta0,nuS0,nuS1,nuD1,nuD0,pc0), &
						   nuD0,nuS1,tauRadInv,phib1,phib2,pc0)

		EceffOverEctot = calcEceff(delta1,nuS0,nuS1,nuD1,nuD0,pc0)/Ectot

		! convert Ectot to SI
		Ectot = sum(Z*nj) * e**3 * lnLStar /(4.0_dp*pi*epsilon0**2 * me * c**2)

		Eceff_out = Ectot * EceffOverEctot
	
	end function getEc_eff
	
	! Collison Frequencies and Effective Momentum !
	pure subroutine get_nuD_nuS(p,Z,Z0,nj,nSpecies,T,element,nuD,nuS)
!~ 	    # Outputs the deflection, slowing down freq (Hesslow 2018)
!~ 		# p = normalized momentum
!~ 		# Z = max. ionization
!~ 		# Z0 = ionization level
!~ 		# nj = density [per cubic m]
!~ 		# T = temperature [eV]
!~ 		# element = 'Ne' or 'Ar'

		! Constants !
		integer, parameter :: dp = kind(1.0d0)
		real(dp), parameter :: e = 1.602176634d-19
		real(dp), parameter :: epsilon0 = 8.8541878128d-12
		real(dp), parameter :: me = 9.1093837015d-31
		real(dp), parameter :: c = 2.99792458d8
		real(dp), parameter :: pi = 3.14159265358979323846d0
		
		real(dp), intent(in) :: p, T
		real(dp), intent(in) :: Z(nSpecies), Z0(nSpecies), nj(nSpecies)
		integer, intent(in) :: nSpecies
		character(len=*), intent(in) :: element
		real(dp), intent(out) :: nuD, nuS

		real(dp) :: ln_aBar(18), ln_I(18)
		real(dp) :: aBar_Ne(18), I_Ne(18)
		real(dp) :: ne_free, Zeff
		real(dp) :: Ne(nSpecies)
		real(dp) :: gamma_rel, beta
		real(dp) :: PT, lnAc, lnA0, lnAee, lnAei
		real(dp) :: hj, F1j, pabar1p5
		integer :: j
		real(dp), parameter :: k = 5.0_dp

		if (trim(element) == "Ne") then
			ln_aBar(1:10) = Ne_ln_aBar_const
			ln_I(1:10)    = Ne_ln_I_const
		else if (trim(element) == "Ar") then
			ln_aBar(1:18) = Ar_ln_aBar_const
			ln_I(1:18)    = Ar_ln_I_const
		else
			nuD = 1.0_dp + sum((Z0**2)*nj) / sum(Z0*nj) ! Just use the fully screened frequencies
			nuS = 1.0_dp
			return
		end if

		aBar_Ne = exp(ln_aBar)
		I_Ne    = exp(ln_I)

		Ne = Z - Z0
		ne_free = sum(Z0*nj)
		Zeff = sum((Z0**2)*nj) / ne_free

		gamma_rel = sqrt(1.0_dp + p*p)
		beta = p*p / (1.0_dp + p*p)

		PT = sqrt( 2.0_dp * T / (me * c*c) )

		lnAc = 14.6_dp + 0.5_dp * log( T / (ne_free/1e20_dp) )
		lnA0 = 14.9_dp - 0.5_dp*log(ne_free/1e20_dp) + log(T/1e3_dp)

		lnAee = lnA0 + 0.2_dp * log( 1.0_dp + ( 2.0_dp*(gamma_rel-1.0_dp)/PT**2 )**2.5_dp )
		lnAei = lnA0 + 0.2_dp * log( 1.0_dp + ( 2.0_dp*p/PT )**5 )

		nuS = lnAee
		nuD = lnAee + Zeff * lnAei

		do j = 2, nSpecies      ! matches python: skip index 0
			hj = p * sqrt(gamma_rel-1.0_dp) / I_Ne(int(Z0(j))+1)

			F1j = (1.0_dp/k) * log(1.0_dp + hj**k) - beta**2

			pabar1p5 = (p * aBar_Ne(int(Z0(j))+1))**1.5_dp

			nuS = nuS + (nj(j)/ne_free) * Ne(j) * F1j

			nuD = nuD + (nj(j)/ne_free) * ( &
				 (2.0_dp/3.0_dp)*(Z(j)**2 - Z0(j)**2) * log(pabar1p5 + 1.0_dp)   &
			   - (2.0_dp/3.0_dp)*(Ne(j)**2) * (pabar1p5/(pabar1p5+1.0_dp)) )
		end do

		nuS = nuS / lnAc
		nuD = nuD / lnAc


  end subroutine get_nuD_nuS
  
  pure function fn(pstar, Estar, Z, Z0, nj, nSpecies, T, element) result(val)
!~ 	   # This describes the function to solve critical momentum
!~     # T [eV]
!~     # Pstar = normalized momentum
!~     # Estar = Epar / E_CH (E_CH = Connor Hastie Critical Field)
	integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: pstar, Estar, T
    real(dp), intent(in) :: Z(nSpecies), Z0(nSpecies), nj(nSpecies)
    integer, intent(in) :: nSpecies
    character(len=*), intent(in) :: element
    real(dp) :: val, nuD, nuS

    call get_nuD_nuS(pstar, Z, Z0, nj, nSpecies, T, element, nuD, nuS)

    val = pstar * sqrt(Estar) - (nuS * nuD)**0.25_dp
  end function fn
  
  pure subroutine get_pstar(E_star, Z, Z0, nj, nSpecies, T, element, pstar)
  ! Calculates effective momentum using Newton's method
    implicit none
    ! Inputs
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: E_star, T
    real(dp), intent(in) :: Z(nSpecies), Z0(nSpecies), nj(nSpecies)
    integer, intent(in) :: nSpecies
    character(len=*), intent(in) :: element

    ! Output
    real(dp), intent(out) :: pstar

    ! Parameters
    real(dp), parameter :: x0 = 1.0_dp
    real(dp), parameter :: eps = 1.0e-3_dp
    real(dp), parameter :: thres = 1.0e-4_dp
    integer, parameter  :: nmax = 100

    ! Locals
    integer :: ii
    real(dp) :: x, xnew, f_x, fp_x, res

    x = x0
    res = 1.0e3_dp
    ii  = 0

    do while (ii < nmax .and. res > thres)

        f_x = fn(x, E_star, Z, Z0, nj, nSpecies, T, element)

        fp_x = ( fn(x+eps, E_star, Z, Z0, nj, nSpecies, T, element) &
               - fn(x-eps, E_star, Z, Z0, nj, nSpecies, T, element) ) / (2.0_dp*eps)

        if (fp_x == 0.0_dp) then
            exit
        end if

        xnew = x - f_x / fp_x

        if (x /= 0.0_dp) then
            res = abs((xnew - x)/x)
        else
            res = abs(xnew - x)
        end if

        x = xnew
        ii = ii + 1
    end do

    pstar = x

   end subroutine get_pstar

  pure subroutine getHesslowAvalanchingRate(E_star,Z,Z0,nj,nSpecies,T,element,B,Yav)
!~ 	# Returns the avalanching rate in per second
!~ 	# E_star [-] = Parallel E-field / E_CH (E_CH = Connor Hastie Electric field)
!~ 	# Z = max. ionization; first element of array is deuterium, next is Impurity (Ne or Ar)
!~ 	# Z0 = ionization level; first element of array is 1, next is Neon ionization
!~ 	# nj = density of each species [per cubic m]
!~ 	# T = temperature [eV]
!~ 	# element = either 'Ar' or 'Ne'
!~  # Yav = output rate [per s]

	implicit none
	integer, parameter :: dp = kind(1.0d0)
	real(dp), parameter :: e = 1.602176634d-19
	real(dp), parameter :: epsilon0 = 8.8541878128d-12
	real(dp), parameter :: me = 9.1093837015d-31
	real(dp), parameter :: c = 2.99792458d8
	real(dp), parameter :: pi = 3.14159265358979323846d0
	
    real(dp), intent(in) :: E_star, T, B
    real(dp), intent(in) :: Z(nSpecies), Z0(nSpecies), nj(nSpecies)
    integer, intent(in) :: nSpecies
    character(len=*), intent(in) :: element
    real(dp), intent(out) :: Yav

    real(dp) :: ne_free, ntot, lnAc
    real(dp) :: Ec, Epar, Eceff
    real(dp) :: pstar
    real(dp) :: nuD, nuS, Zeff
    real(dp) :: Ectmp

    ne_free = sum(Z0 * nj) ! Free electron density [per cubic m]
    ntot    = sum(Z  * nj) ! Free + bound electron density [per cubic m]

    lnAc = 14.6_dp + 0.5_dp * log( T / (ne_free/1e20_dp)) ! Coul. Log

    Ec = e**3 * ne_free * lnAc / (me * c**2 * 4.0_dp * pi * epsilon0**2) ! Connor Hastie Field [V/m]
    Epar = Ec * E_star ! Parallel Electric Field [V/m]

    if (abs(ne_free - ntot) < 1.0e-5_dp) then ! Check if fully ionized
        Zeff = sum((Z0**2)*nj) / ne_free ! Use Rosenbluth-Putvinski Model
        Eceff = Ec
        nuD = 1.0_dp + Zeff
        nuS = 1.0_dp
    else
        Eceff = getEc_eff(Z,Z0,nj,nSpecies,T,B,element)
        call get_pstar(E_star, Z, Z0, nj, nSpecies, T, element, pstar)
        call get_nuD_nuS(pstar, Z, Z0, nj, nSpecies, T, element, nuD, nuS)
    end if
	
	if (Epar > Eceff) then
		Yav = ( e / (me*c*lnAc) ) * (ntot/ne_free) * (Epar - Eceff) / sqrt(4.0_dp + nuS*nuD)
	else
		Yav = 0.0
	end if

  end subroutine getHesslowAvalanchingRate

     
end module runaway_mod

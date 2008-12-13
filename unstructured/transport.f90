module transport_coefficients

contains

! Density Sources/Sinks
! ~~~~~~~~~~~~~~~~~~~~~
vectype function sigma_func(i)
  use basic
  use nintegrate_mod
  use diagnostics

  implicit none

  integer, intent(in) :: i
  integer :: j
  vectype :: temp

  temp = 0.

  ! Pellet injection model
  if(ipellet.eq.1) then
     temp79a = ri_79*pellet_rate/(2.*pi*pellet_var**2) & 
          *exp(-((x_79 - pellet_x)**2 + (z_79 - pellet_z)**2) &
          /(2.*pellet_var**2))
     temp = temp + int2(g79(:,OP_1,i),temp79a)
  endif

  ! Ionization model
  if(ionization.eq.1) then
     temp79d = pt79(:,OP_1) / nt79(:,OP_1)
     
     do j=1,npoints
        if(real(temp79d(j)) .gt. ionization_temp) then
           temp79e(j) = exp(-(temp79d(j) - ionization_temp) &
                / ionization_depth)
        else
           temp79e(j) = 1.
        endif
     enddo
     
     temp79a = ionization_rate * temp79e * &
          exp(-ionization_temp / temp79d)
     temp = temp + int2(g79(:,OP_1,i),temp79a)
  endif

  ! Localized sink(s)
  if(isink.ge.1) then
     temp79a = &
          - nt79(:,OP_1)*ri_79*sink1_rate/(2.*pi*sink1_var**2) & 
          *exp(-((x_79 - sink1_x)**2 + (z_79 - sink1_z)**2) &
          /(2.*sink1_var**2))
     temp = temp + int2(g79(:,OP_1,i),temp79a)
  endif
  if(isink.ge.2) then
     temp79a = &
          - nt79(:,OP_1)*ri_79*sink2_rate/(2.*pi*sink2_var**2) & 
          *exp(-((x_79 - sink2_x)**2 + (z_79 - sink2_z)**2) &
          /(2.*sink2_var**2))
     temp = temp + int2(g79(:,OP_1,i),temp79a)
  endif

  sigma_func = temp
  return
end function sigma_func


! Resistivity
! ~~~~~~~~~~~
vectype function resistivity_func(i)
  use basic
  use nintegrate_mod
  use diagnostics

  implicit none

  integer, intent(in) :: i
  vectype :: temp

  temp = 0.

  if(eta0 .ne. 0) then
     select case (iresfunc)
     case(0)  ! resistivity = 1/Te**(3/2) = sqrt((n/pe)**3)
        if(linear.eq.1) then
           temp79a = eta0*sqrt((n079(:,OP_1)/(pefac*pe079(:,OP_1)))**3)
        else
           temp79a = eta0*sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
        endif

     case(1)      ! added 08/05/08 for stability benchmarking
        temp79a = eta0*.5* &
             (1. + tanh((real(ps079(:,OP_1))-(psilim+etaoff*(psilim-psimin)))&
             /(etadelt*(psilim-psimin))))
     end select
     temp = temp + int2(g79(:,OP_1,i),temp79a)
  endif

  resistivity_func = temp
  return
end function resistivity_func


! Viscosity
! ~~~~~~~~~
vectype function viscosity_func(i)
  use basic
  use nintegrate_mod
  use diagnostics

  implicit none

  integer, intent(in) :: i
  integer :: j
  real :: factor
  vectype :: temp

  temp = 0.

  if(amu_edge.ne.0.) then

     ! added 10/18/08  to make viscosity function more like resistivity 
     ! for iresfunc=1  scj
     select case (iresfunc)
     case(0)
        do j=1,npoints
           call mask(x_79(j)-xzero, z_79(j)-zzero, factor)
           temp79a(j) = amu*amu_edge*(1.-factor)
        end do
        temp = temp + int2(g79(:,OP_1,i),temp79a)
     case(1)
        temp79a = amu_edge*.5* &
             (1. + tanh((real(ps079(:,OP_1))-(psilim+etaoff*(psilim-psimin))) &
             /(etadelt*(psilim-psimin))))
        temp = temp + int2(g79(:,OP_1,i),temp79a)
     end select
  endif

  viscosity_func = temp
  return
end function viscosity_func

! Kappa
! ~~~~~
vectype function kappa_func(i)
  use basic
  use nintegrate_mod
  use diagnostics
  
  implicit none
  
  integer, intent(in) :: i
  vectype :: temp

  if(numvar.lt.3) then
     kappa_func = 0.
     return
  endif

  temp = 0.

  ! kappa = p/T**(3/2) = sqrt(n**3/p)
  if(kappa0.ne.0.) then
     temp79a = kappa0*sqrt(nt79(:,OP_1)**3/pt79(:,OP_1))
     temp = temp + int2(g79(:,OP_1,i),temp79a)
  endif

  if(kappah.ne.0.) then
     ! insert h-mode model here!
  end if
  
  kappa_func = temp
  return
end function kappa_func


! Electron viscosity
! ~~~~~~~~~~~~~~~~~~
vectype function electron_viscosity_func(i)
  use basic
  use nintegrate_mod
  use diagnostics

  implicit none

  integer, intent(in) :: i
  vectype :: temp

  temp = 0.

  if(amue.ne.0) then
     temp79f = -amue * r2_79 * &
          (bzt79(:,OP_DZ)*pst79(:,OP_DZ) + bzt79(:,OP_DR)*pst79(:,OP_DR)) &
          / (nt79(:,OP_1)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 + 1e-1)**2)
     temp = temp + int2(g79(:,OP_1,i),temp79a)
  endif

  electron_viscosity_func = temp
  return
end function electron_viscosity_func
end module transport_coefficients



! define_transport_coefficients
! =============================
subroutine define_transport_coefficients()

  use basic
  use arrays
  use nintegrate_mod
  use newvar_mod
  use sparse
  use transport_coefficients

  implicit none

  integer :: i, itri, ibegin, iendplusone, numnodes
  integer :: ione, numelms, def_fields

  logical, save :: first_time = .true.
  logical :: solve_sigma, solve_kappa, solve_visc, solve_resistivity, &
       solve_visc_e

  if((linear.eq.1).and.(.not.first_time)) return
  first_time = .false.

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating auxiliary variables"

  solve_resistivity = .false.
  solve_visc = .false.
  solve_kappa = .false.
  solve_sigma = .false.
  solve_visc_e = .false.

  resistivity = 0.
  kappa = 0.
  sigma = 0.
  visc = 0.
  visc_e = 0.
  tempvar = 0.

  def_fields = FIELD_N + FIELD_PE + FIELD_P + FIELD_PSI + FIELD_I

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' defining...'

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     call define_fields(itri, def_fields, int_pts_aux, 1)

     do i=1,18
        ione = isval1(itri,i)

        resistivity(ione) = resistivity(ione) + resistivity_func(i)
        if(.not.solve_resistivity) solve_resistivity = resistivity(ione).ne.0.

        kappa(ione) = kappa(ione) + kappa_func(i)
        if(.not.solve_kappa) solve_kappa = kappa(ione).ne.0.

        sigma(ione) = sigma(ione) + sigma_func(i)
        if(.not.solve_sigma) solve_sigma = sigma(ione).ne.0.

        visc(ione) = visc(ione) + viscosity_func(i)
        if(.not.solve_visc) solve_visc = visc(ione).ne.0.

        if(ibootstrap.eq.1) then
           visc_e(ione) = visc_e(ione) + electron_viscosity_func(i)
           if(.not.solve_visc_e) solve_visc_e = visc_e(ione).ne.0.
        endif

        tempvar(ione) = tempvar(ione) &
             + int2(g79(:,OP_1,i),sz79(:,OP_1))
     end do
  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' solving...'

  if(solve_resistivity) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  resistivity'
     call solve_newvar(resistivity, NV_NOBOUND, mass_matrix_lhs)
  end if

  if(solve_kappa) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  kappa'
     call solve_newvar(kappa, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(solve_sigma) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  sigma'
     call solve_newvar(sigma, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(solve_visc) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  viscosity'
     call solve_newvar(visc, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(solve_visc_e) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  electron viscosity'
     call solve_newvar(visc_e, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(myrank.eq.0 .and. iprint.ge.1) print *, '  size field'
  call solve_newvar(tempvar, NV_NOBOUND, mass_matrix_lhs)

  visc_c = visc

  ! add in constant components
  call numnod(numnodes)
  do i=1,numnodes
     call entdofs(1,i,0,ibegin,iendplusone)
     resistivity(ibegin) = resistivity(ibegin) + etar
     visc(ibegin) = visc(ibegin) + amu
     visc_c(ibegin) = visc_c(ibegin) + amuc
     kappa(ibegin) = kappa(ibegin) + kappat
  enddo

end subroutine define_transport_coefficients

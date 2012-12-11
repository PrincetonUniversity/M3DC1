 module transport_coefficients
  use spline

  type(spline1d), private :: kappa_spline

contains

! Density Sources/Sinks
! ~~~~~~~~~~~~~~~~~~~~~
vectype function sigma_func(i)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  integer, intent(in) :: i
  vectype :: temp
  integer :: iregion, j, magnetic_region

  temp = 0.

  ! Pellet injection model
  if(ipellet.eq.1) then
     temp79a = ri_79*pellet_rate/(2.*pi*pellet_var**2) & 
          *exp(-((x_79 - pellet_x)**2 + (z_79 - pellet_z)**2) &
          /(2.*pellet_var**2))
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  endif

!......distributed source added 11/23/2011   (scj)
  if(ipellet.eq.2) then
     temp79b = max(pedge,real(pt79(:,OP_1)))
     temp79a = pellet_rate*den0*(temp79b/p0)**expn
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
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
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  endif

  if(ibeam.eq.1) then 
     temp79a = neutral_beam_deposition(x_79,z_79)
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  end if

  ! Localized sink(s)
  if(isink.ge.1) then
     temp79a = &
          - nt79(:,OP_1)*ri_79*sink1_rate/(2.*pi*sink1_var**2) & 
          *exp(-((x_79 - sink1_x)**2 + (z_79 - sink1_z)**2) &
          /(2.*sink1_var**2))
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  endif
  if(isink.ge.2) then
     temp79a = &
          - nt79(:,OP_1)*ri_79*sink2_rate/(2.*pi*sink2_var**2) & 
          *exp(-((x_79 - sink2_x)**2 + (z_79 - sink2_z)**2) &
          /(2.*sink2_var**2))
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  endif

  ! Enforce density floor
  if(idenfloor.ge.1) then

  endif

  do j=1, npoints
     temp79a(j) = 0.
     iregion = magnetic_region(pst79(j,:), x_79(j), z_79(j))
     if(iregion.ge.1) temp79a(j) = alphadenfloor*( den_edge - nt79(j,OP_1))
  end do
  temp = temp + int2(mu79(:,OP_1,i),temp79a)

  sigma_func = temp
  return
end function sigma_func


! Momentum Sources/Sinks
! ~~~~~~~~~~~~~~~~~~~~~~
vectype function force_func(i)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  integer, intent(in) :: i
  integer :: j
  vectype :: temp

  temp = 0.

  ! Beam source
  if(ibeam.eq.1) then
     temp79a = neutral_beam_deposition(x_79,z_79)
     temp = temp + nb_v*int2(mu79(:,OP_1,i),temp79a)
     if(ivform.eq.0) then
        temp = temp - int4(ri_79,mu79(:,OP_1,i),temp79a,vzt79(:,OP_1))
     else
        temp = temp - int4(r_79,mu79(:,OP_1,i),temp79a,vzt79(:,OP_1))
     endif
  endif

  force_func = temp
  return
end function force_func

! Poloidal Momentum Sources/Sinks
! ~~~~~~~~~~~~~~~~~~~~~~
vectype function pforce_func(i)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  integer, intent(in) :: i
  integer :: iregion, j, magnetic_region
  real :: psimaxl, psiminl
  vectype, dimension(MAX_PTS,OP_NUM) :: psi

  select case(ipforce)
  case(1)

  if(linear.eq.1) then
      psi = ps079
  else
      psi = pst79
  end if
  temp79a = (psi(:,OP_1)-psimin)/(psibound - psimin)

  do j=1, npoints
    temp79b(j) = aforce*(1.-temp79a(j))**nforce  &
          * dforce**2/((temp79a(j) - xforce)**2 + dforce**2)
    iregion = magnetic_region(psi(j,:), x_79(j), z_79(j))
    if(iregion.ge.1) temp79b(j) = 0.
  end do


  pforce_func = int2(mu79(:,OP_1,i),temp79b)
  case(2)

  if(linear.eq.1) then
      psi = ps079
  else
      psi = pst79
  end if
  temp79a = (psi(:,OP_1)-psibound)/(psimin - psibound)
  psimaxl = 0.6
  psiminl = 1.e-3
  temp79b = 1. - (temp79a - psiminl)/(psimaxl - psiminl)
  temp79c = atan2(z_79-zmag,x_79-xmag)
  temp79d = sqrt(max(real((pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)*ri2_79),1.e-6))
  temp79e =  aforce*temp79b**2*nt79(:,OP_1)*cos(temp79c/2.)/temp79d 
  do j=1,npoints
    if(real(temp79a(j)) .lt. psiminl .or. real(temp79a(j)) .gt. psimaxl) temp79e(j) = 0.
  enddo
  
   pforce_func = int2(mu79(:,OP_1,i),temp79e)
  end select
  return
end function pforce_func
vectype function pmach_func(i)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  integer, intent(in) :: i
  integer :: j

! calculate the poloidal mach number
!  do j=1,npoints
!    temp79a(j) = max(real((pst79(j,OP_DR)**2 + pst79(j,OP_DZ)**2)*ri2_79(j)),1.e-6)
!    temp79b(j) = temp79a(j) + bzt79(j,OP_1)**2*ri2_79(j)
!    temp79c(j) = max(real(gam*pt79(j,OP_1)*ni79(j,OP_1)*temp79a(j)/temp79b(j)),1.e-6)
!  enddo
  temp79a = max(real((pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)*ri2_79),1.e-6)
  temp79b = temp79a + bzt79(:,OP_1)**2*ri2_79
  temp79c = max(real(gam*pt79(:,OP_1)*ni79(:,OP_1)*temp79a/temp79b),1.e-6)
! note: temp79c can vanish at x-point and magnetic axis

  temp79d = max(real(r2_79*(pht79(:,OP_DR)**2 + pht79(:,OP_DZ)**2)   &
         + ri4_79*(cht79(:,OP_DR)**2 + cht79(:,OP_DZ)**2)  &
         + 2.*ri_79*(cht79(:,OP_DZ)*pht79(:,OP_DR)-cht79(:,OP_DR)*pht79(:,OP_DZ))),1.e-9)
  temp79e = sqrt(temp79d/temp79c)

  pmach_func = int2(mu79(:,OP_1,i),temp79e)


  return
end function pmach_func


! Heat Sources/Sinks
! ~~~~~~~~~~~~~~~~~~
vectype function q_func(i)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  integer, intent(in) :: i
  integer :: j
  vectype :: temp

  temp = 0.

  ! Pellet injection model
  if(igaussian_heat_source.eq.1) then
     temp79a = ri_79*ghs_rate/(2.*pi*ghs_var**2) & 
          *exp(-((x_79 - ghs_x)**2 + (z_79 - ghs_z)**2) &
          /(2.*ghs_var**2))
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  endif

  ! Beam source
  if(ibeam.eq.1) then
     temp79a = 0.5*neutral_beam_deposition(x_79,z_79)
     temp = temp + (nb_v**2 + nb_dv**2)*int2(mu79(:,OP_1,i),temp79a)
     if(ivform.eq.0) then
        temp = temp &
             - 2.*nb_v*int4(ri_79,mu79(:,OP_1,i),temp79a,vzt79(:,OP_1)) &
             + int5(ri_79,mu79(:,OP_1,i),temp79a,vzt79(:,OP_1),vzt79(:,OP_1))
     else
        temp = temp &
             - 2.*nb_v*int4(r_79,mu79(:,OP_1,i),temp79a,vzt79(:,OP_1)) &
             + int5(r_79,mu79(:,OP_1,i),temp79a,vzt79(:,OP_1),vzt79(:,OP_1))
     endif
  endif

  q_func = temp
  return
end function q_func


! Resistivity
! ~~~~~~~~~~~
vectype function resistivity_func(i)
  use basic
  use m3dc1_nint
  use diagnostics

  implicit none

  integer, intent(in) :: i
  integer :: ii

  select case (iresfunc)
  case(0)  ! resistivity = 1/Te**(3/2) = sqrt((n/pe)**3)
     if(eta0.ne.0.) then
        if(linear.eq.1) then
           temp79a = eta0*sqrt((ne079(:,OP_1)/pe079(:,OP_1))**3)
        else
           if(itemp.eq.1) then
              temp79a = eta0*tet79(:,OP_1)**(-1.5)
           else
              temp79b = max(pedge*pefac,real(pet79(:,OP_1)))
              temp79a = eta0*sqrt((net79(:,OP_1)/temp79b)**3)
           endif
        endif
     else
        temp79a = 0.
     end if


  case(1)      ! added 08/05/08 for stability benchmarking
     if(linear.eq.1) then
       temp79a = eta0*.5* &
          (1. + &
          tanh((real(ps079(:,OP_1))-(psilim+etaoff*(psilim-psimin)))&
          /(etadelt*(psilim-psimin))))
     else
       temp79a = eta0*.5* &
          (1. + &
          tanh((real(pst79(:,OP_1))-(psilim+etaoff*(psilim-psimin)))&
          /(etadelt*(psilim-psimin))))
     endif
  case(2)
     if(linear.eq.1) then
       temp79b = (ps079(:,OP_1)-psimin)/(psibound-psimin)
       temp79a = eta0*.5* &
          (1. + tanh((real(temp79b) - etaoff)/etadelt))
     else
       temp79b = (pst79(:,OP_1)-psimin)/(psibound-psimin)
       temp79a = eta0*.5* &
          (1. + tanh((real(temp79b) - etaoff)/etadelt))
     endif
  case(3)
     temp79a = eta79(:,OP_1) - etar

  case(4)
     temp79a = eta79(:,OP_1) - etar

  case default
     temp79a = 0.
  end select

  resistivity_func = int2(mu79(:,OP_1,i),temp79a)
end function resistivity_func


! Viscosity
! ~~~~~~~~~
vectype function viscosity_func(i)
  use basic
  use m3dc1_nint
  use diagnostics

  implicit none

  integer, intent(in) :: i
  vectype :: temp
  integer :: iregion, j
  integer :: magnetic_region
  vectype, dimension(MAX_PTS,OP_NUM) :: psi

  temp = 0.

  if(amu_edge.ne.0.) then

     select case (ivisfunc)
     case(0)
        temp79a = 0.

     case(1)
        if(linear.eq.1) then
          temp79a = amu_edge*.5* &
             (1. + &
             tanh((real(ps079(:,OP_1))-(psibound+amuoff*(psibound-psimin))) &
             /(amudelt*(psibound-psimin))))
        else
          temp79a = amu_edge*.5* &
             (1. + &
             tanh((real(pst79(:,OP_1))-(psibound+amuoff*(psibound-psimin))) &
             /(amudelt*(psibound-psimin))))
        endif

     case(2)
        if(linear.eq.1) then
           psi = ps079
        else
           psi = pst79
        end if
        temp79b = (psi(:,OP_1)-psimin)/(psibound - psimin)

        do j=1, npoints
           iregion = magnetic_region(psi(j,:), x_79(j), z_79(j))
           if(iregion.eq.2) temp79b(j) = 2. - temp79b(j)
        end do

        temp79a = amu_edge*.5* &
             (1. + tanh((real(temp79b) - amuoff)/amudelt))
        if(amuoff2.ne.0. .and. amudelt2.ne.0.) then
           temp79a = temp79a + amu_edge*.5* &
                (1. + tanh((real(temp79b) - amuoff2)/amudelt2))
           temp79a = temp79a / 2.
        endif

     case(3)
        temp79a = vis79(:,OP_1) - amu
     end select
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  endif

  viscosity_func = temp
  return
end function viscosity_func

! Kappa
! ~~~~~
vectype function kappa_func(i)
  use math
  use read_ascii
  use basic
  use m3dc1_nint
  use diagnostics
  
  implicit none
  
  integer, intent(in) :: i
  integer :: nvals, j, ierr
  real :: val, valp, valpp, pso
  real, allocatable :: xvals(:), yvals(:)
  integer :: magnetic_region
  vectype :: temp

  if(numvar.lt.3) then
     kappa_func = 0.
     return
  endif

  temp = 0.

  select case (ikappafunc)
  case(0)
     temp79b = max(pedge,real(pt79(:,OP_1)))
     ! kappa = p/T**(3/2) = sqrt(n**3/p)

     if(kappa0.eq.0) then
        temp79a = 0.
     else
        temp79a = kappa0*sqrt(nt79(:,OP_1)**3/pt79(:,OP_1))
     end if
        
  case(1)
     if(linear.eq.1) then
        temp79a = kappa0*.5* &
             (1. + &
             tanh((real(ps079(:,OP_1))-(psilim+kappaoff*(psilim-psimin)))&
             /(kappadelt*(psilim-psimin))))
     else
        temp79a = kappa0*.5* &
             (1. + &
             tanh((real(pst79(:,OP_1))-(psilim+kappaoff*(psilim-psimin)))&
             /(kappadelt*(psilim-psimin)))) 
     endif
  case(2)
     if(linear.eq.1) then
        temp79b = (ps079(:,OP_1)-psimin)/(psibound-psimin)
        temp79a = kappa0*.5* &
             (1. + tanh((real(temp79b) - kappaoff)/kappadelt))
     else
        temp79b = (pst79(:,OP_1)-psimin)/(psibound-psimin)
        temp79a = kappa0*.5* &
             (1. + tanh((real(temp79b) - kappaoff)/kappadelt))
     endif
     !
     !.....added 11/26/2011     scj
  case(3)
     ! kappa = sqrt(1./ (p*n))
     if(kappa0.eq.0) then
        temp79a = 0.
     else
        temp79a = kappa0*sqrt(1./(nt79(:,OP_1)*pt79(:,OP_1)))      
     end if
     
  case(10,11)
     if(.not.allocated(kappa_spline%x)) then
        ! Read in m^2/s (10) or normalized units (11)
        nvals = 0
        call read_ascii_column('profile_kappa', xvals, nvals, icol=1)
        call read_ascii_column('profile_kappa', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(6)
        if(ikappafunc.eq.10) then
           yvals = yvals / &
                (l0_norm * b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm))
        end if
        call create_spline(kappa_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
     
     do j=0, npoints
        if(magnetic_region(pst79(j,:),x_79(j),z_79(j)).ne.0) &
             then
           pso = 1.
        else
           pso = (real(pst79(j,OP_1)) - psimin)/(psilim-psimin)
        end if
        call evaluate_spline(kappa_spline,pso,val,valp,valpp)
        temp79a(j) = val
     end do
     
  case default
     temp79a = 0.
  end select
  temp = temp + int2(mu79(:,OP_1,i),temp79a)

  if(kappah.ne.0.) then
     temp79b = (pst79(:,OP_1) - psimin)/(psibound - psimin)
     temp79a = kappah*tanh((real(temp79b) - 1.)/.2)**2
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  end if
  
  kappa_func = temp
  return
end function kappa_func


! Electron viscosity
! ~~~~~~~~~~~~~~~~~~
vectype function electron_viscosity_func(i)
  use basic
  use m3dc1_nint
  use diagnostics

  implicit none

  integer, intent(in) :: i
  vectype :: temp

  temp = 0.

  if(amue.ne.0) then
     temp79f = -amue * r2_79 * &
          (bzt79(:,OP_DZ)*pst79(:,OP_DZ) + bzt79(:,OP_DR)*pst79(:,OP_DR)) &
          / (nt79(:,OP_1)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 + 1e-1)**2)
     temp = temp + int2(mu79(:,OP_1,i),temp79a)
  endif

  electron_viscosity_func = temp
  return
end function electron_viscosity_func


! define_transport_coefficients
! =============================
subroutine define_transport_coefficients()

  use basic
  use arrays
  use m3dc1_nint
  use newvar_mod
  use sparse
  use neutral_beam

  implicit none

  include 'mpif.h'

  integer :: i, itri
  integer :: numelms, def_fields,ier

  logical, save :: first_time = .true.
  logical :: solve_sigma, solve_kappa, solve_visc, solve_resistivity, &
       solve_visc_e, solve_q, solve_f, solve_fp

  integer, dimension(8) :: temp, temp2
  vectype, dimension(dofs_per_element) :: dofs

  ! transport coefficients are only calculated once in linear mode
  if((linear.eq.1).and.(.not.first_time)) return
  first_time = .false.

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating auxiliary variables"

  ! which transport coefficients need matrix solve
  solve_resistivity = .false.
  solve_visc = .false.
  solve_kappa = .false.
  solve_sigma = .false.
  solve_visc_e = .false.
  solve_f = .false.
  solve_q = .false.
  solve_fp = .false.

  ! clear variables
  resistivity_field = 0.
  kappa_field = 0.

  visc_field = 0.
  if(density_source) sigma_field = 0.  
  if(momentum_source) Fphi_field = 0.
  if(heat_source) q_field = 0.
  if(ibootstrap.ne.0) visc_e_field = 0.
  if(ipforce.gt.0) pforce_field = 0.
  if(ipforce.gt.0) pmach_field = 0.

  call finalize(field0_vec)
  call finalize(field_vec)

  ! specify which primitive fields are to be evalulated
  def_fields = FIELD_N + FIELD_PE + FIELD_P + FIELD_PSI + FIELD_I + FIELD_B2I
  if(itemp.ge.1) def_fields = def_fields + FIELD_TE
  if(iresfunc.eq.3 .or. iresfunc.eq.4) def_fields = def_fields + FIELD_ETA
  if(ivisfunc.eq.3) def_fields = def_fields + FIELD_MU
  if(ibeam.eq.1) def_fields = def_fields + FIELD_V
  if(ipforce.gt.0) def_fields = def_fields + FIELD_PHI + FIELD_CHI + FIELD_NI

  if(myrank.eq.0 .and. iprint.ge.2) print *, '  defining...'

  ! Calculate RHS
  numelms = local_elements()
  do itri=1,numelms

     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, def_fields, 1, linear)

     do i=1, dofs_per_element
        dofs(i) = resistivity_func(i)
        if(.not.solve_resistivity) solve_resistivity = dofs(i).ne.0.
     end do
     if(solve_resistivity) &
          call vector_insert_block(resistivity_field%vec,itri,1,dofs,VEC_ADD)

     do i=1, dofs_per_element
        dofs(i) = kappa_func(i)
        if(.not.solve_kappa) solve_kappa = dofs(i).ne.0.
     end do
     if(solve_kappa) &
          call vector_insert_block(kappa_field%vec,itri,1,dofs,VEC_ADD)

     if(density_source) then
        do i=1, dofs_per_element
           dofs(i) = sigma_func(i)
           if(.not.solve_sigma) solve_sigma = dofs(i).ne.0.
        end do
        if(solve_sigma) &
             call vector_insert_block(sigma_field%vec,itri,1,dofs,VEC_ADD)
     end if

     do i=1, dofs_per_element
        dofs(i) = viscosity_func(i)
        if(.not.solve_visc) solve_visc = dofs(i).ne.0.
     end do
     if(solve_visc) &
          call vector_insert_block(visc_field%vec,itri,1,dofs,VEC_ADD)

     if(momentum_source) then 
        do i=1, dofs_per_element
           dofs(i) = force_func(i)
           if(.not.solve_f) solve_f = dofs(i).ne.0.
        end do
        if(solve_f) &
             call vector_insert_block(Fphi_field%vec,itri,1,dofs,VEC_ADD)
     end if

     if(ipforce.gt.0) then
        do i=1, dofs_per_element
           dofs(i) = pforce_func(i)
           if(.not.solve_fp) solve_fp = dofs(i).ne.0.
        end do
        if(solve_fp) &
             call vector_insert_block(pforce_field%vec,itri,1,dofs,VEC_ADD)

        do i=1, dofs_per_element
           dofs(i) = pmach_func(i)
        end do
        if(solve_fp) &
             call vector_insert_block(pmach_field%vec,itri,1,dofs,VEC_ADD)

     end if

     if(heat_source) then
        do i=1, dofs_per_element
           dofs(i) = q_func(i)
           if(.not.solve_q) solve_q = dofs(i).ne.0.
        end do
        if(solve_q) &
             call vector_insert_block(q_field%vec,itri,1,dofs,VEC_ADD)
     end if

     if(ibootstrap.ne.0) then
        do i=1, dofs_per_element
           dofs(i) = electron_viscosity_func(i)
           if(.not.solve_visc_e) solve_visc_e = dofs(i).ne.0.
        end do
        if(solve_visc_e) &
             call vector_insert_block(visc_e_field%vec,itri,1,dofs,VEC_ADD)
     end if
  end do

  ! Solve all the variables that have been defined
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! make sure all processes agree on what needs to be solved
  if(maxrank.gt.1) then 
     temp = 0
     if(solve_resistivity) temp(1) = 1
     if(solve_kappa)       temp(2) = 1
     if(solve_sigma)       temp(3) = 1
     if(solve_visc)        temp(4) = 1
     if(solve_visc_e)      temp(5) = 1
     if(solve_f)           temp(6) = 1
     if(solve_q)           temp(7) = 1
     if(solve_fp)          temp(8) = 1
     call mpi_allreduce(temp,temp2,8,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
     solve_resistivity = temp2(1).eq.1
     solve_kappa       = temp2(2).eq.1
     solve_sigma       = temp2(3).eq.1
     solve_visc        = temp2(4).eq.1
     solve_visc_e      = temp2(5).eq.1
     solve_f           = temp2(6).eq.1
     solve_q           = temp2(7).eq.1
     solve_fp          = temp2(8).eq.1
  end if

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' solving...'

  if(solve_resistivity) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  resistivity'
     call newvar_solve(resistivity_field%vec, mass_mat_lhs)
  end if

  if(solve_kappa) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  kappa'
     call newvar_solve(kappa_field%vec, mass_mat_lhs)
  endif

  if(solve_sigma) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  sigma'
     call newvar_solve(sigma_field%vec, mass_mat_lhs)
  endif

  if(solve_visc) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  viscosity'
     call newvar_solve(visc_field%vec, mass_mat_lhs)
  endif

  if(solve_visc_e) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  electron viscosity'
     call newvar_solve(visc_e_field%vec, mass_mat_lhs)
  endif

  if(solve_f) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  fphi'
     call newvar_solve(Fphi_field%vec, mass_mat_lhs)
  endif

  if(solve_q) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Q'
     call newvar_solve(q_field%vec, mass_mat_lhs)
  endif

  if(solve_fp) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  pforce'
     call newvar_solve(pforce_field%vec, mass_mat_lhs)

     if(myrank.eq.0 .and. iprint.ge.1) print *, '  pmach'
     call newvar_solve(pmach_field%vec, mass_mat_lhs)

  endif

  ! the "compressible" viscosity is the same as the "incompressible"
  ! viscosity up to a constant
  visc_c_field = visc_field

  
  ! add in constant components
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  call add(resistivity_field, etar)
  call add(visc_field, amu)
  call add(visc_c_field, amuc)
  call add(kappa_field, kappat)

  if(myrank.eq.0 .and. iprint.ge.2) &
       print *, 'done define_transport_coefficients'

end subroutine define_transport_coefficients

end module transport_coefficients

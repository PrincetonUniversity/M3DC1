module transport_coefficients
  use spline

  type(spline1d), private :: eta_spline
  type(spline1d), private :: kappa_spline
  type(spline1d), private :: denm_spline
  type(spline1d), private :: amu_spline
  type(spline1d), private :: heatsource_spline
  type(spline1d), private :: particlesource_spline
  type(spline1d), private :: prad_nz_spline
  type(spline1d), private :: cd_spline         ! current drive 
  type(spline1d), private :: coef_spline_L31
  type(spline1d), private :: coef_spline_L32
  type(spline1d), private :: coef_spline_L34
  type(spline1d), private :: coef_spline_alpha
  type(spline1d), private :: coef_spline_fluxavg_iBsq
  type(spline1d), private :: coef_spline_dtedpsit
  type(spline1d), private :: coef_spline_fluxavg_G
contains

! bootstrap model coeffcients L31, L32, L34, alpha, 1/Bsq, dte/dpsit 
! ~~~~~~~~~~~
function bootstrapCoeff_func(col_number)
   use basic
   use m3dc1_nint
   use diagnostics
   use math
   use read_ascii
   use resistive_wall
   use bootstrap
 
   implicit none
 
   vectype, dimension(dofs_per_element) :: bootstrapCoeff_func
   integer, intent(in) :: col_number
   real :: tmin
   integer :: nvals, j, mr, bootstrap_coef_0flag
   real, allocatable :: xvals(:), yvals(:)
   real :: val, valp, valpp, pso, psib
   integer :: izone
   integer, dimension(MAX_PTS) :: izarr
   character(len=256) :: Coeff_filename

   bootstrap_coef_0flag = 0
   if(ibootstrap.ne.0) then
         if (ibootstrap_map_te == 0) then
            Coeff_filename='ProfileJBSCoeff_Psi_L31_32_34_alpha_B2_dtedpsit_G'
         elseif (ibootstrap_map_te == 1) then
            Coeff_filename='ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G'
         endif
   else
      print *, 'Error: Incorrect ibootstrap method 2: to use da/dpsit=da/dte dte/dpsit 1: '
      !(B×∇a).∇ϕ=-‖B_p ‖^2  ∂a/∂ψ ---- ∂/∂ψ=-1/‖B_p ‖^2   1/R ((1/R ψ_z+f_ϕR )  ∂/∂z+(1/R ψ_r-f_ϕz )  ∂/∂r)'
      stop
   endif

   if(ibootstrap.ne.0) then
     if(col_number==2)then
         if(.not.allocated(coef_spline_L31%x)) then
            ! Read in L31 or L32 or L34 or alpha or 1/<B^2> or dte/dpsit from col number 2,3,4,5,6,7
            nvals = 0
            call read_ascii_column(Coeff_filename, xvals, nvals,skip=1, icol=1)
            call read_ascii_column(Coeff_filename, yvals, nvals,skip=1, icol=col_number)
            if(nvals.eq.0) call safestop(6)
           ! print*,'Reading L31,32,34,alpha from column number',col_number
           ! print*,'Reading X',xvals
           ! print*,'Reading Y',yvals
            call create_spline(coef_spline_L31, nvals, xvals, yvals)
            deallocate(xvals, yvals)
         end if
      elseif(col_number==3)then
         if(.not.allocated(coef_spline_L32%x)) then
            ! Read in L31 or L32 or L34 or alpha or 1/<B^2> or dte/dpsit from col number 2,3,4,5,6,7
            nvals = 0
            call read_ascii_column(Coeff_filename, xvals, nvals,skip=1, icol=1)
            call read_ascii_column(Coeff_filename, yvals, nvals,skip=1, icol=col_number)
            if(nvals.eq.0) call safestop(6)
            !print*,'Reading L31,32,34,alpha from column number',col_number
            !print*,'Reading X',xvals
            !print*,'Reading Y',yvals
            call create_spline(coef_spline_L32, nvals, xvals, yvals)
            deallocate(xvals, yvals)
         end if
      elseif(col_number==4)then        
        if(.not.allocated(coef_spline_L34%x)) then
            ! Read in L31 or L32 or L34 or alpha or 1/<B^2> or dte/dpsit from col number 2,3,4,5,6,7
            nvals = 0
            call read_ascii_column(Coeff_filename, xvals, nvals,skip=1, icol=1)
            call read_ascii_column(Coeff_filename, yvals, nvals,skip=1, icol=col_number)
            if(nvals.eq.0) call safestop(6)
            !print*,'Reading L31,32,34,alpha from column number',col_number
            !print*,'Reading X',xvals
            !print*,'Reading Y',yvals
            call create_spline(coef_spline_L34, nvals, xvals, yvals)
            deallocate(xvals, yvals)
         end if
      elseif(col_number==5)then    
         if(.not.allocated(coef_spline_alpha%x)) then
            ! Read in L31 or L32 or L34 or alpha or 1/<B^2> or dte/dpsit from col number 2,3,4,5,6,7
            nvals = 0
            call read_ascii_column(Coeff_filename, xvals, nvals,skip=1, icol=1)
            call read_ascii_column(Coeff_filename, yvals, nvals,skip=1, icol=col_number)
            if(nvals.eq.0) call safestop(6)
            !print*,'Reading L31,32,34,alpha from column number',col_number
            !print*,'Reading X',xvals
            !print*,'Reading Y',yvals
            call create_spline(coef_spline_alpha, nvals, xvals, yvals)
            deallocate(xvals, yvals) 
         end if
      elseif(col_number==6)then    
            if(.not.allocated(coef_spline_fluxavg_iBsq%x)) then
               ! Read in L31 or L32 or L34 or alpha or 1/<B^2> or dte/dpsit from col number 2,3,4,5,6,7
               nvals = 0
               call read_ascii_column(Coeff_filename, xvals, nvals,skip=1, icol=1)
               call read_ascii_column(Coeff_filename, yvals, nvals,skip=1, icol=col_number)
               if(nvals.eq.0) call safestop(6)
               !print*,'Reading L31,32,34,alpha,B from column number',col_number
               !print*,'Reading X',xvals
               !print*,'Reading Y',yvals
               call create_spline(coef_spline_fluxavg_iBsq, nvals, xvals, yvals)
               deallocate(xvals, yvals) 
            end if
      elseif(col_number==7)then    
            if(.not.allocated(coef_spline_dtedpsit%x)) then
               ! Read in L31 or L32 or L34 or alpha or 1/<B^2> or dte/dpsit from col number 2,3,4,5,6,7
               nvals = 0
               call read_ascii_column(Coeff_filename, xvals, nvals,skip=1, icol=1)
               call read_ascii_column(Coeff_filename, yvals, nvals,skip=1, icol=col_number)
               if(nvals.eq.0) call safestop(6)
               !print*,'Reading L31,32,34,alpha,1/<B^2> or dte/dpsit from column number',col_number
               !print*,'Reading X',xvals
               !print*,'Reading Y',yvals
               yvals = yvals * 1.6022e-9 * (4.*pi*n0_norm)/ (b0_norm**2)
               call create_spline(coef_spline_dtedpsit, nvals, xvals, yvals)
               deallocate(xvals, yvals) 
            end if
      elseif(col_number==8)then    
            if(.not.allocated(coef_spline_fluxavg_G%x)) then
               ! Read in L31 or L32 or L34 or alpha or 1/<B^2> or dte/dpsit or G from col number 2,3,4,5,6,7,8
               nvals = 0
               call read_ascii_column(Coeff_filename, xvals, nvals,skip=1, icol=1)
               call read_ascii_column(Coeff_filename, yvals, nvals,skip=1, icol=col_number)
               if(nvals.eq.0) call safestop(6)
               !print*,'Reading L31,32,34,alpha,1/<B^2> or dte/dpsit or G from column number',col_number
               !print*,'Reading X',xvals
               !print*,'Reading Y',yvals
               call create_spline(coef_spline_fluxavg_G, nvals, xvals, yvals)
               deallocate(xvals, yvals) 
            end if
     end if
     
         
        
     
      do j=1, npoints
         if (ibootstrap_map_te == 0) then
            call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
                 x_79(j),z_79(j),mr,psib)
            if(mr.eq.REGION_PF) then
               pso = 2.*psib - pso
            end if
            pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)
         else if (ibootstrap_map_te == 1) then
            if(itemp == 1) then
               pso = real(tet79(j,OP_1))
            else
               pso=pet79(j,OP_1)/net79(j,OP_1)
            endif
            !convert to pso to eV
            pso=pso * (p0_norm / n0_norm) / 1.6022e-12
         endif 
         if(bootstrap_coef_0flag==1)then
            temp79a(j) = 0
         else
            if(col_number==2)then  
               call evaluate_spline(coef_spline_L31,pso,val,valp,valpp,extrapolate=1)
            elseif(col_number==3)then  
               call evaluate_spline(coef_spline_L32,pso,val,valp,valpp,extrapolate=1)
            elseif(col_number==4)then  
               call evaluate_spline(coef_spline_L34,pso,val,valp,valpp,extrapolate=1)
            elseif(col_number==5)then 
               call evaluate_spline(coef_spline_alpha,pso,val,valp,valpp,extrapolate=1) 
            elseif(col_number==6)then 
               call evaluate_spline(coef_spline_fluxavg_iBsq,pso,val,valp,valpp,extrapolate=1) 
            elseif(col_number==7)then 
               call evaluate_spline(coef_spline_dtedpsit,pso,val,valp,valpp,extrapolate=1) 
            elseif(col_number==8)then 
               call evaluate_spline(coef_spline_fluxavg_G,pso,val,valp,valpp,extrapolate=1)    
            endif
            temp79a(j) = val
            !if(myrank.eq.0) print *,'pso,val', pso, val
         endif
      end do
                 
            
         
      bootstrapCoeff_func = intx2(mu79(:,:,OP_1),temp79a)
   end if
 
 end function bootstrapCoeff_func

! Density Sources/Sinks
! ~~~~~~~~~~~~~~~~~~~~~
function sigma_func(izone)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use pellet
  use read_ascii

  implicit none

  vectype, dimension(dofs_per_element) :: sigma_func
  integer, intent(in) :: izone
  vectype, dimension(dofs_per_element) :: temp
  integer :: iregion, j
  integer :: nvals
  real :: val, valp, valpp, pso
  real :: rate, psib
  real, allocatable :: xvals(:), yvals(:)
  integer :: ip, mr

  ! Don't allow particle source in wall or vacuum region
  if(izone.ne.ZONE_PLASMA) then
     sigma_func = 0.
     return
  end if

  temp = 0.

  ! Pellet injection model

  if(ipellet.gt.0 .and. (ipellet_z.eq.0 .or. any(pellet_mix.gt.0.))) then

     do ip=1, npellets
        if(ipellet_abl.gt.0. .and. pellet_var(ip).lt.1.e-8) then
           pellet_var(ip) = 0.
           temp79a = 0.
        else
           if(pellet_mix(ip).eq.0.) then
              rate = pellet_rate(ip)
           else
              if(ipellet_abl.eq.0) then
                 pellet_rate_D2(ip) = pellet_rate(ip)* &
                      pellet_mix(ip) / (1. - pellet_mix(ip))
              end if
              rate = pellet_rate_D2(ip)*2.0 ! two deuterium ions per D2 molecule
           end if
           temp79a = rate*pellet_distribution(ip, x_79, phi_79, z_79, real(pt79(:,OP_1)), 1, izone)
        endif
     
        temp = temp + intx2(mu79(:,:,OP_1),temp79a)
     end do
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
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif

  if(ibeam.eq.1 .or. ibeam.eq.2) then
#ifdef USEST
     temp79a = neutral_beam_deposition(xl_79,zl_79)
#else
     temp79a = neutral_beam_deposition(x_79,z_79)
#endif
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  end if

  ! Read numerical particle source profile
  if(iread_particlesource.eq.1) then
     if(.not.allocated(particlesource_spline%x)) then
        nvals = 0
        call read_ascii_column('profile_particlesource', xvals, nvals, icol=1)
        call read_ascii_column('profile_particlesource', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(6)
        call create_spline(particlesource_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if

     do j=1, npoints
#ifdef USEST
        pso = xl_79(j)**2+zl_79(j)**2 
#else
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j),z_79(j),mr,psib)
        if(mr.eq.REGION_PF) then
           pso = 2.*psib - pso
        end if
        pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)
#endif
        call evaluate_spline(particlesource_spline,pso,val,valp,valpp)
        temp79a(j) = val * pellet_rate_scl
     end do

     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif 

  ! Localized sink(s)
  if(isink.ge.1) then
     temp79a = &
          - nt79(:,OP_1)*ri_79*sink1_rate/(2.*pi*sink1_var**2) & 
          *exp(-((x_79 - sink1_x)**2 + (z_79 - sink1_z)**2) &
          /(2.*sink1_var**2))
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif
  if(isink.ge.2) then
     temp79a = &
          - nt79(:,OP_1)*ri_79*sink2_rate/(2.*pi*sink2_var**2) & 
          *exp(-((x_79 - sink2_x)**2 + (z_79 - sink2_z)**2) &
          /(2.*sink2_var**2))
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif


  if(iarc_source.eq.1) then
     ! temp79a = normal current density = -j.grad(wall_dist)
     temp79a = bzt79(:,OP_DZ)*wall79(:,OP_DR)*ri_79 &
          -    bzt79(:,OP_DR)*wall79(:,OP_DZ)*ri_79
#if defined(USE3D) || defined(COMPLEX)
     temp79a = temp79a  &
          + bfpt79(:,OP_DZP)*wall79(:,OP_DR)*ri_79 &
          - bfpt79(:,OP_DRP)*wall79(:,OP_DZ)*ri_79 &
          - pst79(:,OP_DRP)*wall79(:,OP_DR)*ri2_79 &
          - pst79(:,OP_DZP)*wall79(:,OP_DZ)*ri2_79
#endif
     temp79b = arc_source_alpha*temp79a* &
          (wall79(:,OP_1)/arc_source_eta)*exp(-wall79(:,OP_1)/arc_source_eta)
     where(real(temp79b).lt.0.)
        temp79b = 0.
     end where
     temp = temp + intx2(mu79(:,:,OP_1),temp79b)
  end if

  ! Enforce density floor
  if(idenfloor.ge.1) then
     temp79a = 0.
     do j=1, npoints
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j), z_79(j), iregion)
        if(iregion.ne.REGION_PLASMA) &
             temp79a(j) = alphadenfloor*( den_edge - nt79(j,OP_1))
     end do
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif

  sigma_func = temp
end function sigma_func


! Momentum Sources/Sinks
! ~~~~~~~~~~~~~~~~~~~~~~
function force_func(izone)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  vectype, dimension(dofs_per_element) :: force_func
  integer, intent(in) :: izone
  vectype, dimension(dofs_per_element) :: temp

  ! Don't allow momentum source in wall or vacuum region
  if(izone.ne.ZONE_PLASMA) then
     force_func = 0.
     return
  end if

  temp = 0.

  ! Beam source
  if(ibeam.eq.1 .or. ibeam.eq.4 .or. ibeam.eq.5) then
     temp79a = neutral_beam_deposition(x_79,z_79)
     temp = temp + nb_v*beam_fracpar*intx2(mu79(:,:,OP_1),temp79a)
     temp = temp - intx4(mu79(:,:,OP_1),r_79,temp79a,vzt79(:,OP_1))
  endif

  force_func = temp
end function force_func

! Poloidal Momentum Sources/Sinks
! ~~~~~~~~~~~~~~~~~~~~~~
function pforce_func()
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  vectype, dimension(dofs_per_element) :: pforce_func
  integer :: iregion, j
  real :: psimaxl, psiminl

  select case(ipforce)
  case(1)

     temp79a = (pst79(:,OP_1)-psimin)/(psibound - psimin)

     do j=1, npoints
        temp79b(j) = aforce*(1.-temp79a(j))**nforce  &
             * dforce**2/((temp79a(j) - xforce)**2 + dforce**2)
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j), z_79(j),iregion)
        if(iregion.ne.REGION_PLASMA) temp79b(j) = 0.
     end do

     pforce_func = intx2(mu79(:,:,OP_1),temp79b)
  case(2)

     temp79a = (pst79(:,OP_1)-psibound)/(psimin - psibound)
     psimaxl = 0.6
     psiminl = 1.e-3
     temp79b = 1. - (temp79a - psiminl)/(psimaxl - psiminl)
     temp79c = atan2(z_79-zmag,x_79-xmag)
     temp79d = sqrt(max(real((pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)*ri2_79),1.e-6))
     temp79e =  aforce*temp79b**2*nt79(:,OP_1)*cos(temp79c/2.)/temp79d 
     do j=1,npoints
        if(real(temp79a(j)) .lt. psiminl .or. real(temp79a(j)) .gt. psimaxl) temp79e(j) = 0.
     enddo
  
     pforce_func = intx2(mu79(:,:,OP_1),temp79e)

  case default
     pforce_func = 0.

  end select

end function pforce_func


function pmach_func()
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam

  implicit none

  vectype, dimension(dofs_per_element) :: pmach_func
!  integer :: j

! calculate the poloidal mach number
!  do j=1,npoints
!    temp79a(j) = max(real((pst79(j,OP_DR)**2 + pst79(j,OP_DZ)**2)*ri2_79(j)),1.e-6)
!    temp79b(j) = temp79a(j) + bzt79(j,OP_1)**2*ri2_79(j)
!    temp79c(j) = max(real(gam*pt79(j,OP_1)*ni79(j,OP_1)*temp79a(j)/temp79b(j)),1.e-6)
!  enddo
  temp79a = max(real((pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)*ri2_79),1.e-6)
  temp79b = temp79a + bzt79(:,OP_1)**2*ri2_79
  temp79c = max(real(gam*pt79(:,OP_1)*ni79(:,OP_1)*temp79a/temp79b),1.e-7)
! note: temp79c can vanish at x-point and magnetic axis

  temp79d = max(real(r2_79*(pht79(:,OP_DR)**2 + pht79(:,OP_DZ)**2)   &
         + ri4_79*(cht79(:,OP_DR)**2 + cht79(:,OP_DZ)**2)  &
         + 2.*ri_79*(cht79(:,OP_DZ)*pht79(:,OP_DR)-cht79(:,OP_DR)*pht79(:,OP_DZ))),1.e-9)
  temp79e = sqrt(temp79d/temp79c)

  pmach_func = intx2(mu79(:,:,OP_1),temp79e)
end function pmach_func


! ==================================================
! Heat Sources/Sinks
! ~~~~~~~~~~~~~~~~~~
!
! NOTE: When adding heat source and radiation source, make sure that
! heat_source and rad_source are set to .true.
! ==================================================
function q_func(izone)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii
  use radiation
  use basicq

  implicit none

  vectype, dimension(dofs_per_element) :: q_func
  integer, intent(in) :: izone
  vectype, dimension(dofs_per_element) :: temp
  integer :: nvals, j, mr
  real :: val, valp, valpp, pso, rsq, psib
  real, allocatable :: xvals(:), yvals(:)
  real, dimension(MAX_PTS) :: r

  ! Don't allow heating in wall or vacuum region
  if(izone.ne.ZONE_PLASMA) then
     q_func = 0.
     return
  end if

  temp = 0.

  ! Gaussian heat source model
  if(igaussian_heat_source.eq.1) then
#ifdef USEST
     temp79a = ri_79*ghs_rate/(2.*pi*ghs_var**2) & 
          *exp(-((xl_79 - ghs_x)**2 + (zl_79 - ghs_z)**2) &
          /(2.*ghs_var**2))
#else
     temp79a = ri_79*ghs_rate/(2.*pi*ghs_var**2) & 
          *exp(-((x_79 - ghs_x)**2 + (z_79 - ghs_z)**2) &
          /(2.*ghs_var**2))
#endif
#ifdef USE3D
     if(ghs_var_tor .gt. 0) then
        if(itor.eq.1) then 
           temp79a = temp79a*exp( &
                -2.*r_79*ghs_x*(1.-cos(phi_79 - ghs_phi)) &
                /(2.*ghs_var_tor**2)) &
                / (sqrt(2.*pi)*ghs_var_tor)
        else
           temp79b = &
                (exp(-(phi_79-ghs_phi)**2/(2.*ghs_var_tor**2)) &
                +exp(-(phi_79-ghs_phi+toroidal_period)**2/(2.*ghs_var_tor**2)) &
                +exp(-(phi_79-ghs_phi-toroidal_period)**2/(2.*ghs_var_tor**2))) &
                / (sqrt(2.*pi)*ghs_var_tor)
           temp79a = temp79a*temp79b
        end if
     end if
#endif
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif

  ! Beam source
  if(ibeam.ge.1 .and. ibeam.le.4) then
     temp79a = 0.5*neutral_beam_deposition(x_79,z_79)
     temp = temp + (nb_v**2 + nb_dv**2)*intx2(mu79(:,:,OP_1),temp79a)
     temp = temp &
          - 2.*nb_v*intx4(mu79(:,:,OP_1),r_79,temp79a,vzt79(:,OP_1)) &
          + intx5(mu79(:,:,OP_1),r2_79,temp79a,vzt79(:,OP_1),vzt79(:,OP_1))
  endif

  ! Read numerical heat source profile
  if(iread_heatsource.eq.1) then
     if(.not.allocated(heatsource_spline%x)) then
        nvals = 0
        call read_ascii_column('profile_heatsource', xvals, nvals, icol=1)
        call read_ascii_column('profile_heatsource', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(6)
           yvals = yvals * ghs_rate
        call create_spline(heatsource_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if

     do j=1, npoints
#ifdef USEST
        pso = (xl_79(j)**2+zl_79(j)**2) 
#else
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j),z_79(j), mr, psib)
        if(mr.eq.REGION_PF) then
           pso = 2.*psib - pso
        end if
        pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)
#endif
        call evaluate_spline(heatsource_spline,pso,val,valp,valpp)
        temp79a(j) = val
     end do

     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif

  ! Heat sink for use with itaylor=27
  if(iheat_sink.eq.1 .and. itaylor.eq.27) then
     r = sqrt((x_79-xmag)**2 + (z_79-zmag)**2)
     do j=1,npoints
        rsq = r(j)**2
!       temp79a(j) = coolrate*(pfunc(rsq)-pt79(j,OP_1))
        temp79a(j) = pefac*coolrate*(pfunc(rsq)) ! now use new time p in pressure_lin
     end do
     temp79a = temp79a*(1. + tanh((r-libetap)/p1))
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  endif

  q_func = temp
end function q_func

function totrad_func(itri)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii
  use radiation
  use basicq
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: itri

  vectype, dimension(dofs_per_element) :: totrad_func
  vectype, dimension(dofs_per_element) :: temp
  integer :: j, ierr, nvals
  real :: val, pso
  real, allocatable :: xvals(:), yvals(:)

  temp = 0.

  ! Radiation
  if(iprad.eq.1) then

     if(itemp.eq.1) then
        temp79b = tet79(:,OP_1)
     else
        temp79b = pet79(:,OP_1)/net79(:,OP_1)
     end if

     ! convert temperature to keV
     temp79b = temp79b * (p0_norm / n0_norm) / 1.6022e-12 / 1000.

     ! convert density to /m^3
     temp79c = net79(:,OP_1) * n0_norm * 1e6

    
     if(iread_prad.eq.1) then 
        if(.not.allocated(prad_nz_spline%x)) then
           ! Read in units of 10^20 / m^3
           nvals = 0
           call read_ascii_column('profile_nz', xvals, nvals, icol=1)
           call read_ascii_column('profile_nz', yvals, nvals, icol=2)
           if(nvals.eq.0) call safestop(6)
           yvals = yvals * 1e14 / n0_norm
           call create_spline(prad_nz_spline, nvals, xvals, yvals)
           deallocate(xvals, yvals)
        end if

        ! evaluate impurity density
        do j=1, npoints
           pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)
           call evaluate_spline(prad_nz_spline,pso,val)
        
           temp79d(j) = val
        end do

        ! convert impurity density to /m^3
        temp79d = temp79d * n0_norm * 1e6
     else
      
        ! impurity density is prad_fz * electron density
        temp79d = temp79c*prad_fz
     end if

     ierr = 0
     do j=1, npoints
        call get_Prad_simple(temp79a(j), temp79b(j), temp79d(j), &
             prad_z, temp79c(j), ierr)
     end do

     ! convert output to normalized units
     temp79a = temp79a * 10. / (p0_norm / t0_norm)
     
     temp = temp - intx2(mu79(:,:,OP_1),temp79a)
  end if

  if(ikprad.ne.0) then
     call eval_ops(itri, kprad_rad, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
     call eval_ops(itri, kprad_brem, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
     call eval_ops(itri, kprad_ion, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
     call eval_ops(itri, kprad_reck, tm79, rfac) ! only kinetic recombination
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
  end if

  totrad_func = temp
end function totrad_func

function linerad_func(itri)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii
  use radiation
  use basicq
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: itri

  vectype, dimension(dofs_per_element) :: linerad_func
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  if(ikprad.ne.0) then
     call eval_ops(itri, kprad_rad, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
  end if

  linerad_func = temp
end function linerad_func

function bremrad_func(itri)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii
  use radiation
  use basicq
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: itri

  vectype, dimension(dofs_per_element) :: bremrad_func
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  if(ikprad.ne.0) then
     call eval_ops(itri, kprad_brem, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
  end if

  bremrad_func = temp
end function bremrad_func

function ionrad_func(itri)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii
  use radiation
  use basicq
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: itri

  vectype, dimension(dofs_per_element) :: ionrad_func
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  if(ikprad.ne.0) then
     call eval_ops(itri, kprad_ion, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
  end if

  ionrad_func = temp
end function ionrad_func

function reckrad_func(itri)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii
  use radiation
  use basicq
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: itri

  vectype, dimension(dofs_per_element) :: reckrad_func
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  if(ikprad.ne.0) then
     call eval_ops(itri, kprad_reck, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
  end if

  reckrad_func = temp
end function reckrad_func

function recprad_func(itri)
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii
  use radiation
  use basicq
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: itri

  vectype, dimension(dofs_per_element) :: recprad_func
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  if(ikprad.ne.0) then
     call eval_ops(itri, kprad_recp, tm79, rfac)
     temp = temp - intx2(mu79(:,:,OP_1),tm79(:,OP_1))
  end if

  recprad_func = temp
end function recprad_func



! Current Drive sources
! ~~~~~~~~~~~~~~~~~~
function cd_func()
  use math
  use basic
  use m3dc1_nint
  use diagnostics
  use neutral_beam
  use read_ascii

  implicit none

  vectype, dimension(dofs_per_element) :: cd_func
  integer :: iregion, j, nvals
  vectype, dimension(dofs_per_element) :: temp
  real, allocatable :: xvals(:), yvals(:)

  temp = 0.

  ! Gaussian source
  if(icd_source.eq.1) then
     do j=1,npoints
        temp79a(j) = J_0cd * exp( -(x_79(j)-R_0cd)**2/w_cd**2 &
             - (z_79(j)-Z_0cd)**2/w_cd**2 ) - delta_cd
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j),z_79(j),iregion)
        if(iregion.ne.REGION_PLASMA) temp79a(j) = 0.
     enddo
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
#ifdef USEST
  else if(icd_source.eq.2) then
     temp79b = sqrt((xl_79-xcenter)**2 +(zl_79-zcenter)**2)
     temp79a = J_0cd/sqrt(2.*pi*w_cd**2) & 
          *exp(-(temp79b - delta_cd)**2/(2.*w_cd**2))
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  else if(icd_source.eq.3) then
     if(.not.allocated(cd_spline%x)) then
        nvals = 0
        call read_ascii_column('profile_cd', xvals, nvals, icol=1)
        call read_ascii_column('profile_cd', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(6)
        yvals = yvals / 7.96e5 
        call create_spline(cd_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     endif

     temp79b = sqrt((xl_79-xcenter)**2 +(zl_79-zcenter)**2)
     do j=1,npoints
        call evaluate_spline(cd_spline, temp79b(j)**2, temp79a(j))
     end do
     temp79a = temp79a*r_79
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
#endif
  endif

  cd_func = temp
end function cd_func

! Resistivity
! ~~~~~~~~~~~
function resistivity_func(izone_index)
  use basic
  use m3dc1_nint
  use diagnostics
  use math
  use read_ascii
  use resistive_wall

  implicit none

  vectype, dimension(dofs_per_element) :: resistivity_func
  integer, intent(in) :: izone_index
  real :: tmin
  integer :: nvals, j, mr, iz
  real, allocatable :: xvals(:), yvals(:)
  real :: val, valp, valpp, pso, psib
  integer :: izone
  integer, dimension(MAX_PTS) :: izarr

  izone = zone_type(izone_index)

  if(izone.eq.ZONE_PLASMA) then
     select case (iresfunc)
     case(0)  ! resistivity = 1/Te**(3/2) = sqrt((n/pe)**3)
        if(eta0.ne.0.) then
           tmin = (eta_fac*eta0/(eta_max - etar*eta_fac))**(2./3.)

           if(itemp.eq.1) then
              temp79b = tet79(:,OP_1) - eta_te_offset
           else
              temp79b = pet79(:,OP_1)/net79(:,OP_1) - eta_te_offset
           endif
           where(real(temp79b).lt.tmin)
              temp79a = eta_max - etar*eta_fac
           elsewhere
              temp79a = eta_fac*eta0*temp79b**(-1.5)
           end where
           where(real(temp79a).lt.(eta_min-etar*eta_fac))
              temp79a = eta_min - etar*eta_fac
           endwhere

        else
           temp79a = 0.
        end if

     case(1)      ! added 08/05/08 for stability benchmarking
          temp79a = eta_fac*eta0*.5* &
               (1. + &
               tanh((real(pst79(:,OP_1))-(psilim+etaoff*(psilim-psimin)))&
               /(etadelt*(psilim-psimin))))

     case(2)
        temp79a = eta79(:,OP_1) - etar*eta_fac

     case(3)
        temp79a = eta79(:,OP_1) - etar*eta_fac

     case(4)
        temp79a = eta79(:,OP_1) - etar*eta_fac

     case(5)  ! resistivity = 1/Te**(3/2) = sqrt((n/pe)**3)/(1 - 2 sqrt(eps))
              ! neoclassical:  Park, et al NF 30 2413 (1990)
        if(eta0.ne.0.) then
           if(itemp.eq.1) then
              temp79c = eta_fac*eta0*tet79(:,OP_1)**(-1.5)
           else
              temp79c = eta_fac*eta0*sqrt((net79(:,OP_1)/pet79(:,OP_1))**3)
           endif
        else
           temp79c = 0.
        endif
        temp79b = sqrt(((x_79 - xmag)**2 + (z_79 - zmag)**2)/rzero**2)
        temp79a = temp79c/(1. - 1.46*sqrt(temp79b))

     case(10,11)
        if(.not.allocated(eta_spline%x)) then
           ! Read in ohm-m (10) or normalized units (11)
           nvals = 0
           call read_ascii_column('profile_eta', xvals, nvals, icol=1)
           call read_ascii_column('profile_eta', yvals, nvals, icol=2)
           if(nvals.eq.0) call safestop(6)
           if(iresfunc.eq.10) then
              yvals = yvals / (c_light**2 / 1.e11)            ! convert from ohm-m to s
              yvals = yvals / (4.*pi*(v0_norm/c_light)**2 * t0_norm)  ! convert from s
           end if
           call create_spline(eta_spline, nvals, xvals, yvals)
           deallocate(xvals, yvals)
        end if
        
        do j=1, npoints
           call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
                x_79(j),z_79(j),mr,psib)
           if(mr.eq.REGION_PF) then
              pso = 2.*psib - pso
           end if
           pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)

           call evaluate_spline(eta_spline,pso,val,valp,valpp)
           temp79a(j) = val
           if(myrank.eq.0) print *, pso, val
        end do

#ifdef USEST
     case(21)
        if(igeometry.eq.1) then
           temp79b = sqrt((xl_79-xcenter)**2 + (zl_79-zcenter)**2 + regular**2)
           temp79a = eta0* &
                (1. + tanh((temp79b-(1.+etaoff))/etadelt))
        else
           if(myrank.eq.0) print *, 'iresfunc = 21 requires igeometry = 1'
           call safestop(73)
        end if
#endif     

     case default
        if(myrank.eq.0) print *, 'Error: invalid value for iresfunc: ', iresfunc
        call safestop(73)

     end select
  else if(izone.eq.ZONE_CONDUCTOR) then
     izarr = iz
     temp79a = wall_resistivity(x_79,phi_79,z_79,izarr) - etar*eta_fac
  else if(izone.eq.ZONE_VACUUM) then
     temp79a = eta_vac - etar*eta_fac
  end if

  resistivity_func = intx2(mu79(:,:,OP_1),temp79a)
end function resistivity_func


! Viscosity
! ~~~~~~~~~
function viscosity_func()
  use basic
  use m3dc1_nint
  use diagnostics
  use read_ascii
  use basicj

  implicit none

  vectype, dimension(dofs_per_element) :: viscosity_func
  integer :: iregion, j, nvals, mr
  real :: val, valp, valpp, pso, rsq, psib
  real, allocatable :: xvals(:), yvals(:)

  temp79a = 0.

  select case (ivisfunc)
  case(0)
     temp79a = 0.
     
  case(1)
     temp79a = amu_edge*.5* &
          (1. + &
          tanh((real(pst79(:,OP_1))-(psibound+amuoff*(psibound-psimin))) &
          /(amudelt*(psibound-psimin))))
     
  case(2)
     temp79b = (pst79(:,OP_1)-psimin)/(psibound - psimin)
     
     do j=1, npoints
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j), z_79(j), iregion)
        if(iregion.eq.REGION_PF) temp79b(j) = 2. - temp79b(j)
     end do
     
     temp79a = amu_edge*.5* &
          (1. + tanh((real(temp79b) - amuoff)/amudelt))
     if(amuoff2.ne.0. .and. amudelt2.ne.0.) then
        temp79a = temp79a + amu_edge*.5* &
             (1. + tanh((real(temp79b) - amuoff2)/amudelt2))
        temp79a = temp79a / 2.
     endif

  case(3,4)
     temp79a = vis79(:,OP_1) - amu
     
  case(10,11)
     if(.not.allocated(amu_spline%x)) then
        ! Read in m^2/s (10) or normalized units (11)
        nvals = 0
        call read_ascii_column('profile_amu', xvals, nvals, icol=1)
        call read_ascii_column('profile_amu', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(6)
        if(ivisfunc.eq.10) then
           yvals = yvals * 1e4 / (l0_norm * v0_norm)
        end if
        call create_spline(amu_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
     
     do j=1, npoints
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j),z_79(j),mr,psib)
        if(mr.eq.REGION_PF) then
           pso = 2.*psib - pso
        end if
        pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)

        call evaluate_spline(amu_spline,pso,val,valp,valpp)
        temp79a(j) = val
     end do
     if(ivisfunc.eq.10) temp79a = temp79a*nt79(:,OP_1)

  case(12)          !  option to go with itaylor=27, iresfunc=4
     do j=1, npoints
        rsq = (x_79(j)-xmag)**2+(z_79(j)-zmag)**2
        if(itaylor.eq.29) then
           val = amu_edge*basicj_dscale(rsq) - amu
        end if
        temp79a(j) = val
     end do

#ifdef USEST
  case(21)
     if(igeometry.eq.1) then
        temp79b = sqrt((xl_79-xcenter)**2 + (zl_79-zcenter)**2 + regular**2)
        temp79a = amu_edge* &
             (1. + tanh((temp79b-(1.+amuoff))/amudelt))
     else
        if(myrank.eq.0) print *, 'ivisfunc = 21 requires igeometry = 1'
        call safestop(73)
     end if
#endif     

  case default
     if(myrank.eq.0) print *, 'Error: invalid value for ivisfunc: ', ivisfunc
     call safestop(73)
     
  end select

  if(amu_wall.ne.0) then
     temp79a = temp79a + amu_wall* &
          2.*(1. - tanh((wall79(:,OP_1) - amu_wall_off)/amu_wall_delt))
  end if

  viscosity_func = intx2(mu79(:,:,OP_1),temp79a)
end function viscosity_func

! Kappa
! ~~~~~
function kappa_func()
  use math
  use read_ascii
  use basic
  use m3dc1_nint
  use diagnostics
  use basicq
  use basicj

  implicit none

  vectype, dimension(dofs_per_element) :: kappa_func  
  integer :: nvals, j, iregion
  real :: val, valp, valpp, pso, rsq, psib
  real, allocatable :: xvals(:), yvals(:)
  vectype, dimension(dofs_per_element) :: temp

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
     temp79a = kappa0*.5* &
          (1. + &
          tanh((real(pst79(:,OP_1))-(psilim+kappaoff*(psilim-psimin)))&
          /(kappadelt*(psilim-psimin)))) 

  case(2)
     temp79b = (pst79(:,OP_1)-psimin)/(psibound - psimin)
     
     do j=1, npoints
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j), z_79(j), iregion, psib)
        if(iregion.eq.REGION_PF) temp79b(j) = 2.*psib - temp79b(j)
     end do

     temp79a = kappa0*.5* &
          (1. + tanh((real(temp79b) - kappaoff)/kappadelt))

     !
     !.....added 11/26/2011     scj
  case(3)
     ! kappa = sqrt(1./ (p*n))
     if(kappa0.eq.0) then
        temp79a = 0.
     else
        temp79b = max(den_edge*pedge,real(nt79(:,OP_1))*real(pt79(:,OP_1)))
        temp79a = kappa0*sqrt(1./temp79b)      
     end if

  case(4)
     !....added 3/4/2014      scj
     if(itemp.eq.1) then
        temp79a = kappa0*(1. + kappadelt*(tet79(:,OP_DR)*tet79(:,OP_DR) &
             + tet79(:,OP_DZ)*tet79(:,OP_DZ)))
#if defined(USE3D) || defined(USECOMPLEX)
        if(itor.eq.1) temp79a = temp79a + kappa0*kappadelt*tet79(:,OP_DP)*tet79(:,OP_DP)*ri2_79
#endif
     else
        temp79b = pet79(:,OP_DR)**2 + pet79(:,OP_DZ)**2
        temp79c = net79(:,OP_DR)**2 + net79(:,OP_DZ)**2
        temp79d = pet79(:,OP_DR)*net79(:,OP_DR) + pet79(:,OP_DZ)*net79(:,OP_DZ)
#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = temp79b + pet79(:,OP_DP)**2 * ri2_79
        temp79c = temp79c + net79(:,OP_DP)**2 * ri2_79
        temp79d = temp79d + pet79(:,OP_DP)*net79(:,OP_DP) * ri2_79
#endif

        temp79a = kappa0*(1. + kappadelt* &
             (temp79b/net79(:,OP_1)**2 &
             +temp79c*pet79(:,OP_1)**2 / net79(:,OP_1)**4 &
             -2.*temp79d*pet79(:,OP_1) / net79(:,OP_1)**3))
     end if
     
  case(5)
     ! kappa ~ 1/Te (with a maximum)

     temp79a = kap79(:,OP_1) - kappat

  case(10,11)
     if(.not.allocated(kappa_spline%x)) then
        ! Read in m^2/s (10) or normalized units (11)
        nvals = 0
        call read_ascii_column('profile_kappa', xvals, nvals, icol=1)
        call read_ascii_column('profile_kappa', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(6)
        if(ikappafunc.eq.10) then
           yvals = yvals * 1e4 / (l0_norm * v0_norm)
        end if
        call create_spline(kappa_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
     
     do j=1, npoints
        pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)
        call  magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j),z_79(j), iregion, psib)
        if(iregion.eq.REGION_PF) then
           pso = 2.*psib - pso
        end if
        call evaluate_spline(kappa_spline,pso,val,valp,valpp)
        temp79a(j) = val
     end do
     if(ikappafunc.eq.10 .and. itemp.eq.0) temp79a = temp79a*nt79(:,OP_1)

  case(12)          !  option to go with itaylor=27, iresfunc=4
     do j=1, npoints
        rsq = (x_79(j)-xmag)**2+(z_79(j)-zmag)**2
        if(itaylor.eq.27) then 
           val = get_kappa(rsq)
        elseif(itaylor.eq.29) then
           val = kappa0*basicj_dscale(rsq) - kappat
        end if
        temp79a(j) = val
     end do

#ifdef USEST
  case(21)
     if(igeometry.eq.1) then
        temp79b = sqrt((xl_79-xcenter)**2 + (zl_79-zcenter)**2 + regular**2)
        temp79a = kappa0 * &
             (1. + tanh((temp79b-(1.+kappaoff))/kappadelt))
     else
        if(myrank.eq.0) print *, 'ikappafunc = 21 requires igeometry = 1'
        call safestop(73)
     end if
#endif

  case default
     if(myrank.eq.0) print *, 'Error: invalid value for ikappafunc: ', ikappafunc
     call safestop(73)
  end select

  temp79a = temp79a + kappat

  ! BCL 9/30/19: ikappafunc condition added here since defined in m3dc1_nint.f90
  if(kappaf.ge.0. .and. gradp_crit.ne.0 .and. ikappafunc.ne.5) then
     temp79b = pt79(:,OP_DR)**2 + pt79(:,OP_DZ)**2
#ifdef USE3D
     temp79b = temp79b + ri2_79*pt79(:,OP_DP)**2
#endif
     where(real(temp79b).lt.gradp_crit**2)
        temp79a = temp79a * kappaf
     end where
  end if

  temp = temp + intx2(mu79(:,:,OP_1),temp79a)

  ! BCL 9/30/19: ikappafunc condition added here since defined in m3dc1_nint.f90
  if(kappah.ne.0. .and. ikappafunc.ne.5) then
     temp79b = (pst79(:,OP_1) - psimin)/(psibound - psimin)
     temp79a = kappah*tanh((real(temp79b) - 1.)/.2)**2
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
  end if

  kappa_func = temp
end function kappa_func
function kappar_func()
  use math
  use read_ascii
  use basic
  use m3dc1_nint
  use diagnostics
  use basicq
  use basicj

  implicit none

  vectype, dimension(dofs_per_element) :: kappar_func
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  select case (ikapparfunc)
  case(0)
     temp79a = kappar
        
  case(1)
     temp79a = kappar/( (tcrit/tet79(:,OP_1))**2.5 + 1.)
     where(temp79a.ne.temp79a) temp79a = 0.

  case(2)
     temp79a = kar79(:,OP_1)

  case default
     if(myrank.eq.0) print *, 'Error: invalid value for ikapparfunc: ', ikapparfunc
     call safestop(73)
  end select

  temp = temp + intx2(mu79(:,:,OP_1),temp79a)

  kappar_func = temp
end function kappar_func

! denm
! ~~~~
function denm_func()
  use math
  use read_ascii
  use basic
  use m3dc1_nint
  use diagnostics
  use basicq
  use basicj

  implicit none

  vectype, dimension(dofs_per_element) :: denm_func
  integer :: nvals, j, mr
  real :: val, valp, valpp, pso, psib
  real, allocatable :: xvals(:), yvals(:)

  select case (idenmfunc)
  case(0)
       temp79a = denm

  case(1)
     ! denm ~ 1/Te (with minimum & maximum)

     temp79a = denm79(:,OP_1)
     temp79a = min(real(temp79a),denmmax)
     temp79a = max(real(temp79a),denmmin)

  case(10,11)
     if(.not.allocated(denm_spline%x)) then
        ! Read in m^2/s (10) or normalized units (11)
        nvals = 0
        call read_ascii_column('profile_denm', xvals, nvals, icol=1)
        call read_ascii_column('profile_denm', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(6)
        if(idenmfunc.eq.10) then
           yvals = yvals * 1e4 / (l0_norm * v0_norm)
        end if
        call create_spline(denm_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
     
     do j=1, npoints
        pso = (real(pst79(j,OP_1)) - psimin)/(psibound - psimin)
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j),z_79(j), mr, psib)
        if(mr.eq.REGION_PF) then
           pso = 2.*psib - pso
        end if
        call evaluate_spline(denm_spline,pso,val,valp,valpp)
        temp79a(j) = val
     end do

  case default
     if(myrank.eq.0) print *, 'Error: invalid value for idenmfunc: ', idenmfunc
     call safestop(73)
  end select

  denm_func = intx2(mu79(:,:,OP_1),temp79a)
end function denm_func


! Electron viscosity
! ~~~~~~~~~~~~~~~~~~
function electron_viscosity_func()
  use basic
  use m3dc1_nint
  use diagnostics

  implicit none

  vectype, dimension(dofs_per_element) :: electron_viscosity_func
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.

  if(amue.ne.0) then
     temp79f = -amue * r2_79 * &
          (bzt79(:,OP_DZ)*pst79(:,OP_DZ) + bzt79(:,OP_DR)*pst79(:,OP_DR)) &
          / (nt79(:,OP_1)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 + regular)**2)
     temp = temp + intx2(mu79(:,:,OP_1),temp79f)
  endif

  electron_viscosity_func = temp
end function electron_viscosity_func


function be_func()
  use math
  use basic
  use m3dc1_nint
  use diagnostics

  implicit none

  vectype, dimension(dofs_per_element) :: be_func
!
!   need to define this to be p_perp
  if(kinetic.eq.2) then
     be_func = intx2(mu79(:,:,OP_1),p179(:,OP_1))
  else 
     be_func = 0.
  endif

end function be_func


function al_func()
  use math
  use basic
  use m3dc1_nint
  use diagnostics

  implicit none

  vectype, dimension(dofs_per_element) :: al_func
!
!   need to define this as (p_parallel - p_perp)/B**2
  if(kinetic.eq.2) then
     al_func = intx3(mu79(:,:,OP_1),pe179(:,OP_1),b2i79(:,OP_1))   &
          - intx3(mu79(:,:,OP_1), p179(:,OP_1),b2i79(:,OP_1))
  else
     al_func = 0.
  endif

end function al_func


function bs_func()
  use math
  use basic
  use m3dc1_nint
  use diagnostics

  implicit none

  vectype, dimension(dofs_per_element) :: bs_func

  temp79a = ri2_79* &
          (pstx79(:,OP_DR)**2 + pstx79(:,OP_DZ)**2 + bztx79(:,OP_1)**2)

#if defined(USECOMPLEX) || defined(USE3D)
     temp79b = (bfptx79(:,OP_DR)**2 + bfptx79(:,OP_DZ)**2) &
             + 2.*ri_79*(pstx79(:,OP_DZ)*bfptx79(:,OP_DR) - pstx79(:,OP_DR)*bfptx79(:,OP_DZ))
     temp79c  =  (temp79a  + temp79b )
#else
     temp79c  =  temp79a
#endif

  bs_func = intx2(mu79(:,:,OP_1),temp79c)
end function bs_func


! define_transport_coefficients
! =============================
subroutine define_transport_coefficients()

  use basic
  use arrays
  use m3dc1_nint
  use newvar_mod
  use sparse
  use neutral_beam
  use pellet
  use diagnostics
  use kprad_m3dc1

  implicit none

  include 'mpif.h'

  integer :: itri, izone, izone_index
  integer :: numelms, def_fields,ier

  logical, save :: first_time = .true.
  logical :: solve_sigma, solve_kappa, solve_kappar, solve_visc, solve_resistivity, &
       solve_visc_e, solve_q, solve_totrad, solve_linerad, solve_bremrad, &
       solve_ionrad, solve_reckrad, solve_recprad, solve_cd, solve_f, &
       solve_fp, solve_denm, solve_L31,solve_L32,solve_L34,solve_alpha,&
       solve_fluxavg_iBsq,solve_dtedpsit,solve_fluxavg_G

  integer, parameter :: num_scalars = 24
  integer, dimension(num_scalars) :: temp, temp2
  vectype, dimension(dofs_per_element) :: dofs

  ! transport coefficients are only calculated once in linear mode
  if((linear.eq.1).and.(.not.first_time)) return
  first_time = .false.
  
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating transport coefficients"

  ! which transport coefficients need matrix solve
  solve_resistivity = .false.
  solve_visc = .false.
  solve_kappa = .false.
  solve_kappar = .false.
  solve_denm = .false.
  solve_sigma = .false.
  solve_visc_e = .false.
  solve_f = .false.
  solve_q = .false.
  solve_totrad = .false.
  solve_linerad = .false.
  solve_bremrad = .false.
  solve_ionrad = .false.
  solve_reckrad = .false.
  solve_recprad = .false.
  solve_cd = .false.
  solve_fp = .false.
  solve_L31 = .false.
  solve_L32 = .false.
  solve_L34 = .false.
  solve_alpha = .false.
  solve_fluxavg_iBsq = .false.
  solve_dtedpsit = .false.
  solve_fluxavg_G = .false.

  ! clear variables
  resistivity_field = 0.
  kappa_field = 0.
  kappar_field = 0.
  denm_field = 0.

  visc_field = 0.
  if(density_source) sigma_field = 0.  
  if(momentum_source) Fphi_field = 0.
  if(heat_source) Q_field = 0.
  if(rad_source) then
     Totrad_field = 0.
     Linerad_field = 0.
     Bremrad_field = 0.
     Ionrad_field = 0.
     Reckrad_field = 0.
     Recprad_field = 0.
  end if
  if(icd_source .gt. 0) cd_field = 0.
  if(ibootstrap.ne.0) then
   visc_e_field = 0.
   Jbs_L31_field = 0. 
   Jbs_L32_field = 0.
   Jbs_L34_field = 0.
   Jbs_alpha_field = 0.
   Jbs_fluxavg_iBsq_field = 0.
   Jbs_fluxavg_G_field = 0.
   if(ibootstrap.eq.2) then
      Jbs_dtedpsit_field = 0.
   endif
  end if
  if(ipforce.gt.0) pforce_field = 0.
  if(ipforce.gt.0) pmach_field = 0.

  call finalize(field0_vec)
  call finalize(field_vec)

  ! specify which primitive fields are to be evalulated
  def_fields = FIELD_N + FIELD_PE + FIELD_P + FIELD_PSI + FIELD_I + FIELD_B2I
  if(itemp.ge.1) def_fields = def_fields + FIELD_TE
  if(iresfunc.eq.2 .or. iresfunc.eq.3 .or. iresfunc.eq.4) &
       def_fields = def_fields + FIELD_ETA
  if(idenmfunc.eq.1) def_fields = def_fields + FIELD_DENM
  if(ikappafunc.eq.5 .or. ikapparfunc.eq.1 .or. ikapparfunc.eq.2) &
       def_fields = def_fields + FIELD_KAP
  if(ivisfunc.eq.3 .or. ivisfunc.eq.4) def_fields = def_fields + FIELD_MU
  if(ibeam.ge.1) def_fields = def_fields + FIELD_V
  if(ipforce.gt.0) def_fields = def_fields + FIELD_PHI + FIELD_CHI + FIELD_NI

  if(iarc_source.ne.0) def_fields = def_fields + FIELD_WALL

  if(ibootstrap.ne.0) def_fields = def_fields + FIELD_JBS

  if(myrank.eq.0 .and. iprint.ge.2) print *, '  defining...'

  if(ipellet.ne.0) then
     ! make sure normalization for pellet_distribution defined
     call calculate_Lor_vol
  end if
  
  ! Calculate RHS
  numelms = local_elements()
!!$OMP PARALLEL DO &
!!$OMP& PRIVATE(dofs)
  do itri=1,numelms

     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, def_fields, 1, linear)

     call get_zone_index(itri, izone_index)
     izone = zone_type(izone_index)

     dofs = resistivity_func(izone_index)
     if(.not.solve_resistivity) solve_resistivity = any(dofs.ne.0.)

!!$OMP CRITICAL
     if(solve_resistivity) &
          call vector_insert_block(resistivity_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL

     dofs = kappa_func()
     if(.not.solve_kappa) solve_kappa = any(dofs.ne.0.)
!!$OMP CRITICAL
     if(solve_kappa) &
          call vector_insert_block(kappa_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL

     dofs = kappar_func()
     if(.not.solve_kappar) solve_kappar = any(dofs.ne.0.)
!!$OMP CRITICAL
     if(solve_kappar) &
          call vector_insert_block(kappar_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL

     dofs = denm_func()
     if(.not.solve_denm) solve_denm = any(dofs.ne.0.)
!!$OMP CRITICAL
     if(solve_denm) &
          call vector_insert_block(denm_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL


     if(density_source) then
        dofs = sigma_func(izone)
        if(.not.solve_sigma) solve_sigma = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_sigma) &
             call vector_insert_block(sigma_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
     end if

     dofs = viscosity_func()
     if(.not.solve_visc) solve_visc = any(dofs.ne.0.)
!!$OMP CRITICAL
     if(solve_visc) &
          call vector_insert_block(visc_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL

     if(momentum_source) then 
        dofs = force_func(izone)
        if(.not.solve_f) solve_f = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_f) &
             call vector_insert_block(Fphi_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
     end if
     
     if(ipforce.gt.0) then
        dofs = pforce_func()
        if(.not.solve_fp) solve_fp = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_fp) &
             call vector_insert_block(pforce_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL

        dofs = pmach_func()
!!$OMP CRITICAL
        if(solve_fp) &
             call vector_insert_block(pmach_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
     end if

     if(heat_source) then
        dofs = q_func(izone)
        if(.not.solve_q) solve_q = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_q) &
             call vector_insert_block(Q_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
     end if

     if(rad_source) then
        dofs = totrad_func(itri)
        if(.not.solve_totrad) solve_totrad = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_totrad) &
             call vector_insert_block(Totrad_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
        
        dofs = linerad_func(itri)
        if(.not.solve_linerad) solve_linerad = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_linerad) &
             call vector_insert_block(Linerad_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
        
        dofs = bremrad_func(itri)
        if(.not.solve_bremrad) solve_bremrad = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_bremrad) &
             call vector_insert_block(Bremrad_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
        
        dofs = ionrad_func(itri)
        if(.not.solve_ionrad) solve_ionrad = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_ionrad) &
             call vector_insert_block(Ionrad_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
        
        dofs = reckrad_func(itri)
        if(.not.solve_reckrad) solve_reckrad = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_reckrad) &
             call vector_insert_block(Reckrad_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
        
        dofs = recprad_func(itri)
        if(.not.solve_recprad) solve_recprad = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_recprad) &
             call vector_insert_block(Recprad_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
     end if

     if(icd_source .gt. 0) then
        dofs = cd_func()
        if(.not.solve_cd) solve_cd = any(dofs.ne.0)
        if(solve_cd) &
             call vector_insert_block(cd_field%vec,itri,1,dofs,VEC_ADD)
     end if

     if(ibootstrap.ne.0) then
        dofs = electron_viscosity_func()
        if(.not.solve_visc_e) solve_visc_e = any(dofs.ne.0.)
!!$OMP CRITICAL
        if(solve_visc_e) &
             call vector_insert_block(visc_e_field%vec,itri,1,dofs,VEC_ADD)
!!$OMP END CRITICAL
     end if

   !Adding bootstrap coefficients
   if(ibootstrap.ne.0) then
      !Adding L31
      dofs = bootstrapCoeff_func(2)
      if(.not.solve_L31) solve_L31 = .true. !any(dofs.ne.0.)
!$OMP CRITICAL
      if(solve_L31) &
           call vector_insert_block(Jbs_L31_field%vec,itri,1,dofs,VEC_ADD)
!$OMP END CRITICAL
    
      dofs = bootstrapCoeff_func(3)
      if(.not.solve_L32) solve_L32 = .true. !any(dofs.ne.0.)
!$OMP CRITICAL
      if(solve_L32) &
            call vector_insert_block(Jbs_L32_field%vec,itri,1,dofs,VEC_ADD)
!$OMP END CRITICAL

      dofs = bootstrapCoeff_func(4)
      if(.not.solve_L34) solve_L34 = .true. !any(dofs.ne.0.)
!$OMP CRITICAL
      if(solve_L34) &
            call vector_insert_block(Jbs_L34_field%vec,itri,1,dofs,VEC_ADD)
!$OMP END CRITICAL

      dofs = bootstrapCoeff_func(5)
      if(.not.solve_alpha) solve_alpha = .true. !any(dofs.ne.0.)
!$OMP CRITICAL
      if(solve_alpha) &
            call vector_insert_block(Jbs_alpha_field%vec,itri,1,dofs,VEC_ADD)
!$OMP END CRITICAL
     
      dofs = bootstrapCoeff_func(6)
      if(.not.solve_fluxavg_iBsq) solve_fluxavg_iBsq = .true. !any(dofs.ne.0.)
!$OMP CRITICAL
      if(solve_fluxavg_iBsq) &
            call vector_insert_block(Jbs_fluxavg_iBsq_field%vec,itri,1,dofs,VEC_ADD)
!$OMP END CRITICAL  
   endif
   if(ibootstrap.eq.2) then          
      dofs = bootstrapCoeff_func(7)
      if(.not.solve_dtedpsit) solve_dtedpsit = .true. !any(dofs.ne.0.)
!$OMP CRITICAL
      if(solve_dtedpsit) &
            call vector_insert_block(Jbs_dtedpsit_field%vec,itri,1,dofs,VEC_ADD)
!$OMP END CRITICAL 
   end if

   if(ibootstrap.ne.0) then
      dofs = bootstrapCoeff_func(8)
      if(.not.solve_fluxavg_G) solve_fluxavg_G = .true. !any(dofs.ne.0.)
!$OMP CRITICAL
      if(solve_fluxavg_G) &
            call vector_insert_block(Jbs_fluxavg_G_field%vec,itri,1,dofs,VEC_ADD)
!$OMP END CRITICAL  
   endif
   
  end do
!!$OMP END PARALLEL DO


  ! Solve all the variables that have been defined
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! make sure all processes agree on what needs to be solved
  if(maxrank.gt.1) then 
     temp = 0
     temp2 = 0
     if(solve_resistivity) temp(1) = 1
     if(solve_kappa)       temp(2) = 1
     if(solve_sigma)       temp(3) = 1
     if(solve_visc)        temp(4) = 1
     if(solve_visc_e)      temp(5) = 1
     if(solve_f)           temp(6) = 1
     if(solve_q)           temp(7) = 1
     if(solve_fp)          temp(8) = 1
     if(solve_cd)          temp(9) = 1
     if(solve_totrad)      temp(10) = 1
     if(solve_linerad)     temp(11) = 1
     if(solve_bremrad)     temp(12) = 1
     if(solve_ionrad)      temp(13) = 1
     if(solve_reckrad)     temp(14) = 1
     if(solve_recprad)     temp(15) = 1
     if(solve_denm)        temp(16) = 1
     if(solve_kappar)      temp(17) = 1
     if(solve_L31)         temp(18) = 1
     if(solve_L32)         temp(19) = 1
     if(solve_L34)         temp(20) = 1
     if(solve_alpha)       temp(21) = 1
     if(solve_fluxavg_iBsq)temp(22) = 1
     if(solve_dtedpsit)    temp(23) = 1
     if(solve_fluxavg_G)   temp(24) = 1
     

     call mpi_allreduce(temp, temp2, num_scalars, MPI_INTEGER, &
          MPI_MAX, MPI_COMM_WORLD, ier)

     solve_resistivity = temp2(1).eq.1
     solve_kappa       = temp2(2).eq.1
     solve_sigma       = temp2(3).eq.1
     solve_visc        = temp2(4).eq.1
     solve_visc_e      = temp2(5).eq.1
     solve_f           = temp2(6).eq.1
     solve_q           = temp2(7).eq.1
     solve_fp          = temp2(8).eq.1
     solve_cd          = temp2(9).eq.1
     solve_totrad      = temp2(10).eq.1
     solve_linerad     = temp2(11).eq.1
     solve_bremrad     = temp2(12).eq.1
     solve_ionrad      = temp2(13).eq.1
     solve_reckrad     = temp2(14).eq.1
     solve_recprad     = temp2(15).eq.1
     solve_denm        = temp2(16).eq.1
     solve_kappar      = temp2(17).eq.1
     solve_L31         = temp2(18).eq.1
     solve_L32         = temp2(19).eq.1
     solve_L34         = temp2(20).eq.1
     solve_alpha       = temp2(21).eq.1
     solve_fluxavg_iBsq= temp2(22).eq.1
     solve_dtedpsit    = temp2(23).eq.1
     solve_fluxavg_G   = temp2(24).eq.1
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

  if(solve_kappar) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  kappar'
     call newvar_solve(kappar_field%vec, mass_mat_lhs)
  endif

  if(solve_denm) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  denm'
     call newvar_solve(denm_field%vec, mass_mat_lhs)
  endif

  if(solve_sigma) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  sigma'
     call newvar_solve(sigma_field%vec, mass_mat_lhs_dc)
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
     call newvar_solve(Fphi_field%vec, mass_mat_lhs_dc)
  endif

  if(solve_q) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Q'
     call newvar_solve(Q_field%vec, mass_mat_lhs_dc)
  endif

  if(solve_totrad) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Totrad'
     call newvar_solve(Totrad_field%vec, mass_mat_lhs)
  endif

  if(solve_linerad) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Linerad'
     call newvar_solve(Linerad_field%vec, mass_mat_lhs)
  endif

  if(solve_bremrad) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Bremrad'
     call newvar_solve(Bremrad_field%vec, mass_mat_lhs)
  endif

  if(solve_ionrad) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Ionrad'
     call newvar_solve(Ionrad_field%vec, mass_mat_lhs)
  endif

  if(solve_reckrad) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Reckrad'
     call newvar_solve(Reckrad_field%vec, mass_mat_lhs)
  endif

  if(solve_recprad) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  Recprad'
     call newvar_solve(Recprad_field%vec, mass_mat_lhs)
  endif

  if(solve_cd) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' cd'
     call newvar_solve(cd_field%vec, mass_mat_lhs)
  endif

  if(solve_fp) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  pforce'
     call newvar_solve(pforce_field%vec, mass_mat_lhs)

     if(myrank.eq.0 .and. iprint.ge.1) print *, '  pmach'
     call newvar_solve(pmach_field%vec, mass_mat_lhs)

  endif

  if(solve_L31) then
   if(myrank.eq.0 .and. iprint.ge.1) print *, '  Jbs_L31'
   call newvar_solve(Jbs_L31_field%vec, mass_mat_lhs)
  endif

  if(solve_L32) then
   if(myrank.eq.0 .and. iprint.ge.1) print *, '  Jbs_L32'
   call newvar_solve(Jbs_L32_field%vec, mass_mat_lhs)
  endif

  if(solve_L34) then
   if(myrank.eq.0 .and. iprint.ge.1) print *, '  Jbs_L34'
   call newvar_solve(Jbs_L34_field%vec, mass_mat_lhs)
  endif

  if(solve_alpha) then
   if(myrank.eq.0 .and. iprint.ge.1) print *, '  Jbs_alpha'
   call newvar_solve(Jbs_alpha_field%vec, mass_mat_lhs)
  endif
  
  if(solve_fluxavg_iBsq) then
   if(myrank.eq.0 .and. iprint.ge.1) print *, '  Jbs_fluxavg_iBsq'
   call newvar_solve(Jbs_fluxavg_iBsq_field%vec, mass_mat_lhs)
  endif

  if(solve_dtedpsit) then
   if(myrank.eq.0 .and. iprint.ge.1) print *, '  Jbs_dtepdpsit'
   call newvar_solve(Jbs_dtedpsit_field%vec, mass_mat_lhs)
  endif

  if(solve_fluxavg_G) then
   if(myrank.eq.0 .and. iprint.ge.1) print *, '  Jbs_fluxavg_G'
   call newvar_solve(Jbs_fluxavg_G_field%vec, mass_mat_lhs)
  endif
  ! the "compressible" viscosity is the same as the "incompressible"
  ! viscosity up to a constant
  visc_c_field = visc_field

  
  ! add in constant components
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  call add(resistivity_field, etar*eta_fac)
  call add(visc_field, amu)
  call add(visc_c_field, amuc)


  ! Read LP data
  if(iread_lp_source.gt.0) then
     call read_lp_source('cloud.txt', ier)
     if(ier.ne.0) then
        if(myrank.eq.0) print *, 'Error reading LP source ', 'cloud.txt'
        call safestop(8)
     end if
  end if

  if(myrank.eq.0 .and. iprint.ge.2) &
       print *, 'done define_transport_coefficients'
       
end subroutine define_transport_coefficients

end module transport_coefficients

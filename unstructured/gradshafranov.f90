module gradshafranov

  use field
  use spline
  use coils

  implicit none

  integer, parameter :: numvargs = 1

  type(spline1d), private :: omega_spline   ! ion toroidal angular frequency 
  type(spline1d), private :: omega_spline_0 ! profile read from profile_omega
  type(spline1d), private :: alpha_spline   ! 0.5*rzero**2 * omega**2 / Ti 
  type(spline1d), private :: n0_spline      ! Ion density
  type(spline1d), private :: g0_spline      ! 0.5*(F**2 - F1**2)
  type(spline1d), private :: p0_spline      ! Total pressure
  type(spline1d), private :: g2_spline, g3_spline
  type(spline1d), private :: te_spline
#ifdef USEPARTICLES
  type(spline1d) :: rho_spline      ! Ion density
  type(spline1d) :: nf_spline      ! Ion density
  type(spline1d) :: tf_spline      ! Ion density
  type(spline1d) :: nfi_spline      ! Ion density
  type(spline1d) :: tfi_spline      ! Ion density
#endif
  
  type(spline1d), private :: psi_spline     ! Psi as a function of rho
  type(spline1d), private :: ffprime_spline ! F*F'
  type(spline1d), private :: pprime_spline  ! p'

  type(field_type), private :: psi_vec
  type(field_type), private :: fun1_vec, fun2_vec, fun3_vec, fun4_vec

  real, private :: dpsii = 0.
  real, private :: gamma2, gamma3, gamma4  

  logical, private :: constraint = .false.
  logical, private :: do_feedback, do_feedback_x

  real, private :: gnorm, libetapeff, fac2

  integer, private :: int_tor

  integer :: igs_feedfac

  real, dimension(maxfilaments) :: xc_vac, zc_vac
  complex, dimension(maxfilaments) :: ic_vac, ic_out
  integer :: numcoils_vac
  integer, dimension(maxfilaments) :: coil_mask
  integer, dimension(maxfilaments) :: filaments

  real, dimension(maxcoils) :: gs_vertical_feedback
  real, dimension(maxcoils) :: gs_radial_feedback
  real, dimension(maxcoils) :: gs_vertical_feedback_i
  real, dimension(maxcoils) :: gs_radial_feedback_i
  real, dimension(maxcoils) :: gs_vertical_feedback_x
  real, dimension(maxcoils) :: gs_radial_feedback_x
  real, dimension(maxcoils) :: gs_vertical_feedback_x_i
  real, dimension(maxcoils) :: gs_radial_feedback_x_i
  real :: xmag0, zmag0, xmagi, zmagi
  real :: xnull0, znull0, xnulli, znulli

  integer :: igs_start_xpoint_search
  integer :: igs_forcefree_lcfs

  logical :: igs_calculate_pf_fields = .false.
  logical :: igs_calculate_ip_fields = .false.

  real :: gs_pf_psi_width
  real :: tiedge

contains

subroutine gradshafranov_init()

  use basic
  use arrays
  use diagnostics

  implicit none

  real :: tstart, tend

  ! Define initial values of psi
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.gt.0) print *, "before call to vacuum field"
  if(ifixedb.eq.0) then
     igs_calculate_pf_fields = .true.
     igs_calculate_ip_fields = .true.
     call vacuum_field
  end if
  if(myrank.eq.0 .and. iprint.gt.0) print *, "after call to vacuum field"
     
  ! define initial field associated with delta-function or gaussian source
  !     corresponding to current tcuro at location (xmag,zmag)
  if(myrank.eq.0 .and. iprint.gt.0) write(*,1001) xmag,zmag,tcuro,sigma0
1001 format(' in gradshafranov_init',/,'   xmag,zmag,tcuro,sigma0',1p4e12.4)
  if(sigma0 .eq.0) then
     call deltafun(xmag,zmag,tcuro,jphi_field)
  else
     call gaussianfun(xmag,zmag,tcuro,sigma0,jphi_field)
  endif

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(igs.ne.0) call gradshafranov_solve()
  if(myrank.eq.0 .and. itimer.eq.1) then 
     call second(tend)
     t_gs = tend - tstart
  endif

  call gradshafranov_per()

end subroutine gradshafranov_init

!==============================================================
! gradshafranov_per
! ~~~~~~~~~~~~~~~~~
! Imposes perturbations on GS solution
!==============================================================
subroutine gradshafranov_per()
  use init_common

  implicit none

  call init_perturbations
end subroutine gradshafranov_per

subroutine coil_feedback(itnum)
  use basic
  use arrays
  use coils
  use math

  implicit none

  include 'mpif.h'

  integer, intent(in) :: itnum

  integer :: i, ierr
  
  if(myrank.eq.0 .and. iprint.ge.2) then 
     print *, 'Doing feedback', xmag-xmag0, zmag-zmag0,xnull-xnull0,znull-znull0
  end if

  if(do_feedback) then
     xmagi = xmagi + (xmag-xmag0)
     zmagi = zmagi + (zmag-zmag0)
  else
     xmagi = 0.
     zmagi = 0.
  end if

  if((itnum.gt.10).and.(do_feedback_x)) then
     xnulli = xnulli + (xnull-xnull0)
     znulli = znulli + (znull-znull0)
  else
     xnulli = 0.
     znulli = 0.
  end if

  if(myrank.eq.0) then
     do i=1, numcoils_vac
        ic_out(i) = ic_vac(i)

        ! Do magnetic axis control
        if((xmag0.gt.0.).and.(do_feedback)) then
           ic_out(i) = ic_out(i) &
                + (amu0 * 1000. / twopi) / filaments(i) * &
                (gs_vertical_feedback(coil_mask(i))*(zmag-zmag0) &
                +gs_radial_feedback(coil_mask(i))*(xmag-xmag0) &
                +gs_vertical_feedback_i(coil_mask(i))*zmagi &
                +gs_radial_feedback_i(coil_mask(i))*xmagi)
        end if

        ! Do x-point control
        if((xnull0.gt.0. .and. itnum.gt.10).and.(do_feedback_x)) then
           ic_out(i) = ic_out(i) &
                + (amu0 * 1000. / twopi) / filaments(i) * &
                (gs_vertical_feedback_x(coil_mask(i))*(znull-znull0) &
                +gs_radial_feedback_x(coil_mask(i))*(xnull-xnull0) &
                +gs_vertical_feedback_x_i(coil_mask(i))*znulli &
                +gs_radial_feedback_x_i(coil_mask(i))*xnulli)
        end if
     end do
  end if
  call mpi_bcast(ic_out, numcoils_vac, &
       MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)

  ! Field due to coil currents
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating fields due to coils in feedback loop"
  psi_coil_field = 0.
  call field_from_coils(xc_vac,zc_vac,ic_out,numcoils_vac, &
       psi_coil_field,0,ierr)
  if(myrank.eq.0 .and. iprint.ge.2) &
       print *, "Done calculating fields due to coils"
  if(ierr.ne.0) call safestop(5)
 
end subroutine coil_feedback

subroutine write_feedback(filename)
  use basic
  use math

  implicit none

  character(len=*), intent(in) :: filename
  integer :: ic, i, subcoils
  integer, parameter :: ifile = 123
  
  if(myrank.ne.0 .or. numcoils_vac.eq.0) return

  open(unit=ifile,file=filename,action='write',status='replace')

  ic = coil_mask(1)
  subcoils = 0
  do i=1, numcoils_vac
     if(coil_mask(i).ne.ic) then
        write(ifile, '(2F12.4)') &
             real(ic_out(i-1))*twopi/(1000.*amu0)*subcoils, 0.0
        ic = coil_mask(i)
        subcoils = 0
     end if
     
     subcoils = subcoils + 1
  end do
  write(ifile, '(2F12.4)') &
       real(ic_out(numcoils_vac))*twopi/(1000.*amu0)*subcoils, 0.0

  close(ifile)
end subroutine write_feedback

subroutine pf_coil_field(ierr)
  use basic
  use coils
  use math
  use coil_sets
  use arrays

  implicit none

  integer, intent(out) :: ierr

  real, dimension(maxfilaments) :: xc, zc
  complex, dimension(maxfilaments) :: ic
  integer :: ipole, numcoils
  real :: aminor, bv, rnorm, fac

  ierr = 0
  rnorm = 10.
  ! based on filiment with current tcuro
  ! and vertical field of strength bv given by shafranov formula
  ! NOTE:  This formula assumes (li/2 + beta_P) = libetap
  fac  = tcuro/twopi
  fac2 = tcuro / (8.*pi**2*xmag)
  ! minor radius
  aminor = abs(xmag-xlim)
  if(itor.eq.1) then
     bv =  alog(8.*xmag/aminor) - 1.5 + libetap
     libetapeff = libetap
  else
     bv = 0.
  endif
  
  ipole = 0
  select case(idevice)
  case(-1)
     call load_coils(xc,zc,ic,numcoils,'coil.dat','current.dat',coil_mask,&
          filaments)
     xc_vac = xc
     zc_vac = zc
     ic_vac = ic
     numcoils_vac = numcoils
     
  case(0) ! Generic

     if(myrank.eq.0) print *, "Using generic (dipole) configuration"

     numcoils = 1
     xc(1) = 102.
     zc(1) = rnorm
     ipole = 1
     ic = bv*fac2

  case default
     numcoils = 0
     
  end select

  ! Field due to coil currents
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating fields due to coils"
  call field_from_coils(xc,zc,ic,numcoils,psi_coil_field,ipole,ierr)
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Done calculating fields due to coils"
  if(ierr.ne.0) call safestop(5)

  ! Field due to extra divertor currents
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating fields due to divertors"
  if(divertors.ge.1) then
     xc(1:2) = xdiv
     zc(1) = zdiv
     if(divertors.eq.2) zc(2) = -zdiv
     ic(1:2) = fac*divcur
     call field_from_coils(xc,zc,ic,divertors,psi_coil_field,0,ierr)
     if(ierr.ne.0) call safestop(5)
  endif


end subroutine pf_coil_field

!==========================================================
! vacuum_field
! ~~~~~~~~~~~~
! Calculate the field due to external coils
!==========================================================
subroutine vacuum_field()
  use math
  use basic
  use arrays
  use coils

  implicit none

  real, dimension(maxfilaments) :: xp, zp, xc, zc
  complex, dimension(maxfilaments) :: ic
  integer :: ierr
  real, dimension(1,1) :: g1, g2
  real, parameter :: rnorm = 10.
  
  gnorm = 0
  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " calculating vacuum field...."

  !......define feedback parameters needed for normalization
  if(idevice .eq. 0) then
     xc(1) = 102.
     zc(1) = rnorm
     xp = xlim
     zp = zlim
     call gvect0(xp,zp,1,xc,zc,1,g1,1,ierr)
     xp = xlim2
     zp = zlim2
     call gvect0(xp,zp,1,xc,zc,1,g2,1,ierr)
     gnorm = g1(1,1) - g2(1,1)
  endif

  if(igs_calculate_pf_fields) then 
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "Calculating fields due to poloidal field coils"
     call pf_coil_field(ierr)
  end if
 
  if(igs_calculate_ip_fields) then 
     ! Field due to plasma current
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "Calculating fields due to plasma"
     xc(1) = xmag
     zc(1) = zmag
     ic(1) = tcuro/(2.*pi)
     call field_from_coils(xc,zc,ic,1,psi_field(0),0,ierr)
     if(ierr.ne.0) call safestop(5)
  end if

  if(icsubtract.eq.0) call add(psi_field(0),psi_coil_field)

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Done calculating fields"
end subroutine vacuum_field

!======================================================================
! define_profiles
! ~~~~~~~~~~~~~~~
! set up profile splines
!======================================================================
subroutine define_profiles
  use math
  use read_ascii
  use basic
  use iterdb
  implicit none

  real, allocatable :: xvals(:), yvals(:)
  real :: teold, pval, ppval
  integer :: nvals, i, ierr
  type(spline1d) :: fpol_spline, bscale_spline, pscale_spline
  real :: fval, fpval

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining profiles'

  ! define pressure profile
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, ' defining pressure profile'
  ierr = 0
  select case(iread_p)
     
  case(1)
     ! Read in J/m^3
     nvals = 0
     call read_ascii_column('profile_p', xvals, nvals, icol=1)
     call read_ascii_column('profile_p', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals * 10. / (b0_norm**2/(4.*pi))

  case default

  end select

  if(allocated(yvals)) then
     call destroy_spline(p0_spline)
     call create_spline(p0_spline, nvals, xvals, yvals)
     deallocate(xvals, yvals)

     call destroy_spline(pprime_spline)
     call copy_spline(pprime_spline, p0_spline)
     do i=1, pprime_spline%n
        call evaluate_spline(p0_spline, pprime_spline%x(i), pval, ppval)
        pprime_spline%y(i) = ppval
     end do
  end if


  ! define toroidal field profile
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, ' defining R*Bphi profile'
  ierr = 0
  select case(iread_f)
     
  case(1)
     ! Read in T*m
     nvals = 0
     call read_ascii_column('profile_f', xvals, nvals, icol=1)
     call read_ascii_column('profile_f', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals * 1e6 / (b0_norm*l0_norm)

  case default

  end select

  if(allocated(yvals)) then
     bzero = yvals(nvals)/rzero
     yvals = 0.5*(yvals**2 - yvals(nvals)**2)
     call destroy_spline(g0_spline)
     call create_spline(g0_spline, nvals, xvals, yvals)
     deallocate(xvals, yvals)

     call destroy_spline(ffprime_spline)
     call copy_spline(ffprime_spline, g0_spline)
     do i=1, ffprime_spline%n
        call evaluate_spline(g0_spline, ffprime_spline%x(i), pval, ppval)
        ffprime_spline%y(i) = ppval
     end do
  end if


  ! If p' and ff' profiles are not yet defined, define them
  if(.not.allocated(p0_spline%x)) then
     if(inumgs .eq. 1) then
        ! read p' and ff' profiles
        call readpgfiles
     else
        ! use analytic p' and ff' profiles
        call default_profiles
     end if
  else
     constraint = .true.
  end if

  if(iread_te.eq.20 .or. iread_ne.eq.20 .or. iread_omega.eq.20) then
     call load_iterdb('iterdb', ierr)
     if(ierr.ne.0) call safestop(5)
  end if

  ! scale profiles
  if(myrank.eq.0 .and. iprint.ge.2) print *, ' scaling profiles'
  p0_spline%y = p0_spline%y*pscale
  pprime_spline%y = pprime_spline%y*pscale
  g0_spline%y = g0_spline%y*bscale**2
  ffprime_spline%y = ffprime_spline%y*bscale**2
  bzero = bzero*bscale

  if(iread_pscale.eq.1) then
     nvals = 0
     call read_ascii_column('profile_pscale', xvals, nvals, icol=1)
     call read_ascii_column('profile_pscale', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)

     call create_spline(pscale_spline, nvals, xvals, yvals)
     deallocate(xvals, yvals)

     do i=1, pprime_spline%n
        call evaluate_spline(pscale_spline, pprime_spline%x(i), fval, fpval)
        call evaluate_spline(p0_spline, pprime_spline%x(i), pval)
        pprime_spline%y(i) = pprime_spline%y(i)*fval + pval*fpval
     end do
     do i=1, p0_spline%n
        call evaluate_spline(pscale_spline, p0_spline%x(i), fval)
        p0_spline%y(i) = p0_spline%y(i)*fval
     end do

     call destroy_spline(pscale_spline)
  end if

  if(bzero.gt.0) then 
     ffprime_spline%y(:) = ffprime_spline%y(:)* &
          (bpscale**2 + &
          (1.-bpscale)*bpscale*bzero*rzero / &
          sqrt(2.*g0_spline%y(:)+(bzero*rzero)**2))
     g0_spline%y(:) = &
          .5*bpscale*(sqrt(2.*g0_spline%y(:)+(bzero*rzero)**2)-(bzero*rzero)) &
          *(2.*bzero*rzero &
          + (sqrt(2.*g0_spline%y(:)+(bzero*rzero)**2)-(bzero*rzero))*bpscale)
  else
     ffprime_spline%y(:) = ffprime_spline%y(:)* &
          (bpscale**2 - &
          (1.-bpscale)*bpscale*bzero*rzero / &
          sqrt(2.*g0_spline%y(:)+(bzero*rzero)**2))
     g0_spline%y(:) = &
          .5*bpscale*(-sqrt(2.*g0_spline%y(:)+(bzero*rzero)**2)-(bzero*rzero))&
          *(2.*bzero*rzero &
          + (-sqrt(2.*g0_spline%y(:)+(bzero*rzero)**2)-(bzero*rzero))*bpscale)
  endif

  if(iread_bscale.eq.1) then
     nvals = 0
     call read_ascii_column('profile_bscale', xvals, nvals, icol=1)
     call read_ascii_column('profile_bscale', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)

     call create_spline(bscale_spline, nvals, xvals, yvals)
     deallocate(xvals, yvals)

     call copy_spline(fpol_spline, g0_spline)
     fpol_spline%y(:) = sqrt(2.*g0_spline%y(:)+(bzero*rzero)**2)
     if(bzero.lt.0) fpol_spline%y(:) = -fpol_spline%y(:)
        
     do i=1, fpol_spline%n
        call evaluate_spline(bscale_spline, fpol_spline%x(i), fval)
        fpol_spline%y(i) = fpol_spline%y(i)*fval
     end do
     do i=1, g0_spline%n
        call evaluate_spline(fpol_spline, g0_spline%x(i), fval, fpval)
        g0_spline%y(i) = 0.5*(fval**2 - (bzero*rzero)**2)
        ffprime_spline%y(i) = fval*fpval
     end do

     call destroy_spline(fpol_spline)
     call destroy_spline(bscale_spline)
  end if

  ! add pedge to pressure
  if(pedge.gt.0.) p0_spline%y = p0_spline%y - p0_spline%y(p0_spline%n) + pedge

  ! define Te profile
  ! ~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, ' defining temperature profile'
  ierr = 0
  select case(iread_te)
     
  case(1)
     ! Read in keV vs Psi
     nvals = 0
     call read_ascii_column('profile_te', xvals, nvals, icol=1)
     call read_ascii_column('profile_te', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals * 1.6022e-9 / (b0_norm**2/(4.*pi*n0_norm))
  case(2)
     ! Read in eV vs Psi
     nvals = 0
     call read_ascii_column('profile_te', xvals, nvals, icol=1)
     call read_ascii_column('profile_te', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals * 1.6022e-12 / (b0_norm**2/(4.*pi*n0_norm))

  case(4)
     ! Read in keV vs Rho (sqrt Phi)
     nvals = 0
     call read_ascii_column('profile_te_rho_3', xvals, nvals, icol=1)
     call read_ascii_column('profile_te_rho_3', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals * 1.6022e-9 / (b0_norm**2/(4.*pi*n0_norm))
     xvals = xvals / xvals(nvals) ! normalize rho
     call rho_to_psi(nvals, xvals, xvals)

  case(10)
     ! Read from Corsica file (keV vs Psi)
     nvals = 0
     call read_ascii_column('corsica', xvals, nvals, icol=4, &
          read_until='profile data:', skip=2)
     call read_ascii_column('corsica', yvals, nvals, icol=6, &
          read_until='profile data:', skip=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals * 1.6022e-9 / (b0_norm**2/(4.*pi*n0_norm))

  case(20)
     ! Read from iterdb text file (keV vs Psi)
     nvals = idb_nj
     allocate(xvals(nvals), yvals(nvals))
     xvals = idb_psi
     xvals = (xvals(:) - xvals(1)) / (xvals(nvals) - xvals(1))
     yvals = idb_Te
     yvals = yvals * 1.6022e-9 / (b0_norm**2/(4.*pi*n0_norm))

  case default
     
  end select

  if(ierr.ne.0) call safestop(5)

  if(allocated(yvals)) then
     call create_spline(te_spline, nvals, xvals, yvals)
     deallocate(xvals, yvals)
  end if

#ifdef USEPARTICLES
  if ((kinetic.eq.1).and.(kinetic_fast_ion.eq.1)) then
     nvals = 0
     call read_ascii_column('profile_tf_rho', xvals, nvals, icol=1)
     call read_ascii_column('profile_tf_rho', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     xvals = xvals / xvals(nvals) ! normalize rho
     !yvals=yvals*0.25
     if (allocated(psi_spline%y)) then
        call rho_to_psi(nvals, xvals, xvals)
     else
        xvals=xvals**2
     endif
     if(allocated(yvals)) then
        call create_spline(tf_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
  endif
  
     allocate(xvals(101))
     allocate(yvals(101))
     xvals = [(i * 0.01, i = 0, 100)]
     yvals = xvals
     nvals = 101
     if (allocated(psi_spline%y)) then
        call rho_to_psi(nvals, xvals, xvals)
     else
        xvals=xvals**2
     endif
     if(allocated(yvals)) then
        call create_spline(rho_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if

  if ((kinetic.eq.1).and.(kinetic_fast_ion.eq.1)) then
     nvals = 0
     call read_ascii_column('profile_nf_rho', xvals, nvals, icol=1)
     call read_ascii_column('profile_nf_rho', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals / n0_norm !rsae
     xvals = xvals / xvals(nvals) ! normalize rho
     if (allocated(psi_spline%y)) then
        call rho_to_psi(nvals, xvals, xvals)
     else
        xvals=xvals**2
     endif
     if(allocated(yvals)) then
        call create_spline(nf_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
  endif

  if ((kinetic.eq.1).and.(kinetic_thermal_ion.eq.1)) then
     nvals = 0
     call read_ascii_column('profile_ti_rho', xvals, nvals, icol=1)
     call read_ascii_column('profile_ti_rho', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     xvals = xvals / xvals(nvals) ! normalize rho
     ! yvals=yvals*0.5
     if (allocated(psi_spline%y)) then
        call rho_to_psi(nvals, xvals, xvals)
     else
        xvals=xvals**2
     endif
     if(allocated(yvals)) then
        call create_spline(tfi_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
  endif
  
  if ((kinetic.eq.1).and.(kinetic_thermal_ion.eq.1)) then
     nvals = 0
     call read_ascii_column('profile_ni_rho', xvals, nvals, icol=1)
     call read_ascii_column('profile_ni_rho', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals / n0_norm !rsae
     xvals = xvals / xvals(nvals) ! normalize rho
     if (allocated(psi_spline%y)) then
        call rho_to_psi(nvals, xvals, xvals)
     else
        xvals=xvals**2
     endif
     if(allocated(yvals)) then
        call create_spline(nfi_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
  endif
#endif

  ! define density profile
  ! ~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' defining density profile'
  ! If Te is specified but pe equation is not included
  ! then define density based on p and Te
  if(ipres.eq.0 .and. allocated(te_spline%y) .and. eqsubtract.eq.0) then
     if(iread_ne.ne.0) then
        if(myrank.eq.0) &
             print *, 'ERROR: cannot read both ne and Te profiles with&
        & ipres=0 and eqsubtract=0'
        call safestop(17)
     end if

     print *, 'Defining ne = pe/Te'
     call copy_spline(n0_spline, te_spline)
     do i=1, te_spline%n
        call evaluate_spline(p0_spline, n0_spline%x(i), pval)
        n0_spline%y(i) = pefac*pval/(z_ion*te_spline%y(i))
     end do
     call destroy_spline(te_spline)

  else

     ! Otherwise, read ne profile directly
     ierr = 0
     select case(iread_ne)
     
     case(1)
        ! Read in 10^20/m^3 vs Psi
        nvals = 0
        call read_ascii_column('profile_ne', xvals, nvals, icol=1)
        call read_ascii_column('profile_ne', yvals, nvals, icol=2)
        if(nvals.eq.0) then
             if(myrank.eq.0 .and. iprint.ge.0) print *,"nvals=0 when reading profile_ne"
             call safestop(5)
        endif
        yvals = yvals * 1e14 / n0_norm / z_ion
        
     case(2)
        ! Read in 10^19/m^3 vs Psi
        call read_ascii_column('dne.xy', xvals, nvals, skip=3, icol=1)
        call read_ascii_column('dne.xy', yvals, nvals, skip=3, icol=7)
        if(nvals.eq.0) call safestop(5)
        yvals = yvals * 1e13 / n0_norm / z_ion
        
     case(4)
        ! Read in cm^-3 vs Rho (sqrt Phi)
        nvals = 0
        call read_ascii_column('profile_ne_rho_0', xvals, nvals, icol=1)
        call read_ascii_column('profile_ne_rho_0', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(5)
        yvals = yvals / n0_norm / z_ion
        xvals = xvals / xvals(nvals) ! normalize rho
        call rho_to_psi(nvals, xvals, xvals)

     case(10)
        ! Read in corsica (10^20 m^-3 vs Psi)
        nvals = 0
        call read_ascii_column('corsica', xvals, nvals, icol=4, &
          read_until='profile data:', skip=2)
        call read_ascii_column('corsica', yvals, nvals, icol=8, &
          read_until='profile data:', skip=2)
        if(nvals.eq.0) call safestop(5)
        yvals = yvals * 1e14 / n0_norm / z_ion

     case(20)
        ! Read from iterdb text file (m^-3 vs Psi)
        nvals = idb_nj
        allocate(xvals(nvals), yvals(nvals))
        xvals = idb_psi
        xvals = (xvals(:) - xvals(1)) / (xvals(nvals) - xvals(1))
        yvals = idb_ne
        yvals = yvals / 1e6 / n0_norm / z_ion
        
     case default
        call density_profile
        
     end select
     
     if(allocated(yvals)) then
        call create_spline(n0_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
  end if

  if(igs_extend_p.ne.0) call extend_pressure

  ! Set outermost ion temperature to tiedge
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(tiedge.gt.0.) then
     if(pedge.gt.0. .and. myrank.eq.0) then
        print *, 'Warning: tiedge is overriding setting for pedge.'
     end if

     if(allocated(te_spline%y)) then
        ! pedge = ni*(Z*Te + Tiedge)
        pedge = n0_spline%y(n0_spline%n)* &
             (z_ion*te_spline%y(te_spline%n) + tiedge)
     else
        ! (1-pefac) * pedge = ni*Tiedge
        pedge = n0_spline%y(n0_spline%n)*tiedge / (1. - pefac)
     end if

     p0_spline%y = p0_spline%y - p0_spline%y(p0_spline%n) + pedge
  end if

  ! add tedge to temperature
  if(tedge.gt.0.) then
     if(allocated(te_spline%y)) then
        teold = te_spline%y(te_spline%n)
        te_spline%y = te_spline%y - teold + tedge
     else
        teold = pefac*p0_spline%y(p0_spline%n)/n0_spline%y(n0_spline%n)/z_ion
     end if
     ! add difference to pressure profile, so ion temp is not affected.
     if(pedge.le.0.) then
        p0_spline%y = p0_spline%y + n0_spline%n*(tedge - teold)*z_ion
     end if
  end if

  ! Check for negative values in profiles
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do i=1, p0_spline%n
     if(p0_spline%n .le. 0.) then
        if(myrank.eq.0) print *, 'Error: p profile has negative values'
        call safestop(33)
     end if
  end do
  if(allocated(n0_spline%y)) then
     do i=1, n0_spline%n
        if(n0_spline%n .le. 0.) then
           if(myrank.eq.0) print *, 'Error: n profile has negative values'
           call safestop(33)
        end if
     end do
  end if
  if(allocated(te_spline%y)) then 
     do i=1, te_spline%n
        if(te_spline%n .le. 0.) then
           if(myrank.eq.0) print *, 'Error: Te profile has negative values'
           call safestop(33)
        end if
     end do
  end if

  ! define rotation profile
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, ' defining rotation profile'
  if(irot.ne.0) then
     ierr = 0
     select case(iread_omega)

     case(1)
        ! Read in krad/sec
        nvals = 0
        call read_ascii_column('profile_omega', xvals, nvals, icol=1)
        call read_ascii_column('profile_omega', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(5)
        yvals = 1000.* yvals / &
             (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)
        
     case(2)
        ! Read in rad/sec
        call read_ascii_column('dtrot.xy', xvals, nvals, skip=3, icol=1)
        call read_ascii_column('dtrot.xy', yvals, nvals, skip=3, icol=7)
        if(nvals.eq.0) call safestop(5)
        yvals = yvals / &
             (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)

     case(3)
        ! Read in m/sec
        nvals = 0
        call read_ascii_column('profile_vphi', xvals, nvals, icol=1)
        call read_ascii_column('profile_vphi', yvals, nvals, icol=2)
        if(nvals.eq.0) call safestop(5)
        yvals = yvals / &
             (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm) &
             / rzero

     case(4)
        ! Read in rad/sec vs Rho (sqrt Phi)
        nvals = 0
        call read_ascii_column('profile_omega_rho_0',xvals,nvals,icol=1)
        call read_ascii_column('profile_omega_rho_0',yvals,nvals,icol=2)
        if(nvals.eq.0) call safestop(5)
        yvals = yvals / &
             (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)
        xvals = xvals / xvals(nvals) ! normalize rho
        call rho_to_psi(nvals, xvals, xvals)

     case(5)    ! added 1/1/2014 to read J.Menard files   (scj)
        ! Read in rad/sec
        nvals = 0
        call read_ascii_column('profile_omega', xvals, nvals, 1, icol=1)
! Question:  do we need to square this to get normalized poloidal flux?
        call read_ascii_column('profile_omega', yvals, nvals, 1, icol=2)
        if(nvals.eq.0) call safestop(5)
        yvals = yvals / &
             (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)

     case(20)
        ! Read from iterdb text file (rad/s vs Psi)
        nvals = idb_nj
        allocate(xvals(nvals), yvals(nvals))
        xvals = idb_psi
        xvals = (xvals(:) - xvals(1)) / (xvals(nvals) - xvals(1))
        yvals = idb_omega
        yvals = yvals / &
             (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)

     case default
        call default_omega
        call copy_spline(omega_spline_0, omega_spline)
     end select

     if(allocated(yvals)) then
        call create_spline(omega_spline_0, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if

     ! scale rotation
     omega_spline_0%y = omega_spline_0%y*vscale

     ! If we've read in electron rotation,
     ! add in diamagnetic term to get ion rotation
     call calc_omega_profile
  endif

  call unload_iterdb

  ! output profiles
  if(myrank.eq.0 .and. iprint.ge.1) call write_profile

  if(myrank.eq.0) then
     print *, 'pprime at max psi: ', &
          pprime_spline%x(pprime_spline%n), &
          pprime_spline%y(pprime_spline%n)
     print *, 'ffprime at max psi: ', &
          ffprime_spline%x(ffprime_spline%n), &
          ffprime_spline%y(ffprime_spline%n)
  end if

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Done defining profiles'

end subroutine define_profiles

subroutine calc_omega_profile
  use basic
  implicit none

  real :: ppval, nval, np, dia, te0, tep
  integer :: i, iout

  call copy_spline(omega_spline, omega_spline_0)

  if((iread_omega_e.ne.0 .or. iread_omega_ExB.ne.0) .and. db.ne.0.) then 
     if(dpsii.eq.0.) then
        if(myrank.eq.0) then
           print *, 'Error: psi bounds are required when reading omega_e'
        end if
        call safestop(14)
     endif
     
     do i=1, omega_spline%n
        if(omega_spline%x(i).ge.1. .and. igs_extend_diamag.eq.0) cycle
           
        call evaluate_spline(pprime_spline, omega_spline%x(i), ppval, &
             iout=iout)
        if(iout.eq.1) ppval = 0.
        call evaluate_spline(n0_spline, omega_spline%x(i), nval, np, &
             iout=iout)
        if(iout.eq.1) np = 0.
        
        if(iread_omega_e.ne.0) then 
           ! we're reading in electron rotation; 
           ! add full diamagnetic term
           dia = db*ppval/nval
           
        else if(iread_omega_ExB.ne.0) then 
           ! we're reading in ExB rotation; add ion diamagnetic term
           if(allocated(te_spline%y)) then
              call evaluate_spline(te_spline, omega_spline%x(i), &
                   te0,tep,iout=iout)
              if(iout.eq.1) tep = 0.
              dia = db*(ppval/nval - z_ion*(1.+thermal_force_coeff)*tep &
                   - z_ion*te0*np/nval)
           else
              dia = db*(1.-pefac)*ppval/nval
           endif
        endif
        if(iflip_j.eq.1) dia = -dia
        if(iflip_v.eq.1) dia = -dia
        
        omega_spline%y(i) = omega_spline_0%y(i) - dpsii*dia
     end do
  end if
  
  ! calculate alpha from omega
  ! ensure that derivatives at LCFS are zero
  if(igs_forcefree_lcfs.eq.1) then 
     where(omega_spline%x .ge. 1.) omega_spline%y = 0.
  else if(igs_forcefree_lcfs.eq.2) then
     call evaluate_spline(omega_spline_0, 1., te0, iout=iout)
     if(iout.eq.1) te0 = 0.
     where(omega_spline%x .ge. 1.) omega_spline%y = te0
  end if
  
  call calculate_alpha

end subroutine calc_omega_profile


subroutine extend_pressure
  use basic
  use spline
  
  implicit none

  type(spline1d) :: pe_spline, temp_spline
  real :: p0_val, n0_val, ti0, pe0
  real :: pe_val, p_val, pp_val, t_val, n_val, tex, nex, px
  integer :: i

  tex = 0.
  nex = 0.
  px = p0_spline%x(p0_spline%n)
  if(allocated(te_spline%y)) tex = te_spline%x(te_spline%n)
  nex = n0_spline%x(n0_spline%n)
  if(tex.le.px .and. nex.le.px) return

  ! calculate p & ni at last p location
  call evaluate_spline(p0_spline, px, p0_val)
  call evaluate_spline(n0_spline, px, n0_val)

  ! create pe spline
  if(tex.gt.nex) then 
     call copy_spline(pe_spline, te_spline)
  else
     call copy_spline(pe_spline, n0_spline)
  end if

  do i=1, pe_spline%n
        call evaluate_spline(n0_spline, pe_spline%x(i), n_val)
     if(allocated(te_spline%y)) then
        call evaluate_spline(te_spline, pe_spline%x(i), t_val)
        pe_spline%y(i) = t_val*n_val*z_ion
     else
        pe_spline%y(i) = p0_val*(p0-pi0)/p0*(n_val/n0_val)
     end if
  end do

  call evaluate_spline(pe_spline, px, pe0)
  ti0 = (p0_val - pe0) / n0_val
  if(ti0.lt.0.) ti0 = 0.
  if(myrank.eq.0) print *, 'ti0 = ', ti0

  ! create new pressure spline
  call copy_spline(temp_spline, p0_spline)
  call destroy_spline(p0_spline)
  call copy_spline(p0_spline, pe_spline)
  do i=1, p0_spline%n
     if(p0_spline%x(i).le.px) then
        ! if inside px, copy old pressure value
        call evaluate_spline(temp_spline, p0_spline%x(i), p_val)
        p0_spline%y(i) = p_val
     else
        ! otherwise, extend pressure keeping Ti constant
        call evaluate_spline(pe_spline, p0_spline%x(i), pe_val)
        call evaluate_spline(n0_spline, p0_spline%x(i), n_val)
        p0_spline%y(i) = pe_val + n_val*ti0
     end if
  end do
  call destroy_spline(temp_spline)
  call destroy_spline(pe_spline)

  ! re-evaluate pprime_spline
  call destroy_spline(pprime_spline)
  call copy_spline(pprime_spline, p0_spline)
  do i=1, pprime_spline%n
     call evaluate_spline(p0_spline, pprime_spline%x(i), p_val, pp_val)
     pprime_spline%y(i) = pp_val
  end do

!  call safestop(3)

end subroutine extend_pressure

!============================================================
subroutine gradshafranov_solve

  use math
  use mesh_mod
  use basic
  use arrays
  use sparse
  use diagnostics
  use newvar_mod
  use m3dc1_nint
  use matrix_mod
  use boundary_conditions
  use model

  implicit none
  include 'mpif.h'

  type(field_type) :: b1vecini_vec, b2vecini_vec
  type(field_type) :: b3vecini_vec, b4vecini_vec

  type(matrix_type) :: gs_matrix

  integer :: itri, i, ier, itnum, ibound, izone
  integer :: numelms, numnodes
  real :: feedfac

  real :: error, error2, error3 

  vectype :: tf, tf2
  vectype, dimension(dofs_per_element,dofs_per_element) :: temp

  real :: crit_v, max_v
  real :: norm_1, norm_2, dnorm


!!$  integer :: is_edge(3)  ! is inode on boundary
!!$  real :: n(2,3)
!!$  integer :: iedge, idim(3)

  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, "Calculating GS- Equilibrium"

  if(imulti_region.eq.1) then
     int_tor = 5
  else
     int_tor = 0
  end if

  numnodes = owned_nodes()
  numelms = local_elements()

  ! allocate memory for arrays
if (ispradapt .eq. 1) then
  call create_field(b1vecini_vec, "b1vecini_vec")
  call create_field(b2vecini_vec, "b2vecini_vec")
  call create_field(psi_vec, "psi_vec")
  call create_field(fun1_vec, "fun1_vec")
  call create_field(fun2_vec, "fun2_vec")
  call create_field(fun3_vec, "fun3_vec")
  call create_field(fun4_vec, "fun4_vec")
else
  call create_field(b1vecini_vec)
  call create_field(b2vecini_vec)
  call create_field(psi_vec)
  call create_field(fun1_vec)
  call create_field(fun2_vec)
  call create_field(fun3_vec)
  call create_field(fun4_vec)
endif
  psi_vec = psi_field(0)
  b1vecini_vec = jphi_field
  if(iread_eqdsk .eq. 1) psilim = psibound

  ! form the grad-sharfranov matrix
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " forming the GS matrix..."

  call set_matrix_index(gs_matrix, gsmatrix_sm)
  call create_mat(gs_matrix, numvargs, numvargs, icomplex, 1)

  if(int_tor.eq.0) then
     ibound = BOUNDARY_DIRICHLET + BOUNDARY_AXISYMMETRIC
  else
     ibound = BOUNDARY_DIRICHLET
  endif

  ! populate the matrix
  do itri=1,numelms

     call define_element_quadrature(itri,int_pts_main,int_tor)
     call define_fields(itri,0,1,0)

     temp = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_GS),ri_79)
#ifdef USE3D
     temp = temp - eta_gs*intxx3(mu79(:,:,OP_DP),nu79(:,:,OP_DP),ri3_79)
#endif
!!$     temp = &
!!$         -intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),ri_79) &
!!$         -intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),ri_79)
!!$     if(itor.eq.1) then 
!!$         temp = temp &
!!$             -intxx3(mu79(:,:,OP_1),nu79(:,:,OP_DR),ri2_79)
!!$     endif
!!$ 
!!$     ! add surface terms
!!$     call boundary_edge(itri, is_edge, n, idim)
!!$     
!!$     do iedge=1,3
!!$        if(is_edge(iedge).eq.0) cycle
!!$
!!$        call define_boundary_quadrature(itri, iedge, 5, 5, n, idim)
!!$        call define_fields(itri, 0, 1, 0)
!!$
!!$        temp = temp &
!!$             +intxx4(mu79(:,:,OP_1),nu79(:,:,OP_DR),norm79(:,1),ri_79) &
!!$             +intxx4(mu79(:,:,OP_1),nu79(:,:,OP_DZ),norm79(:,2),ri_79)
!!$     end do

     call apply_boundary_mask(itri, ibound, temp, tags=BOUND_DOMAIN)

     call insert_block(gs_matrix, itri, 1, 1, temp, MAT_ADD)
  enddo

  feedfac = 0.

  ! insert boundary conditions
  call flush(gs_matrix)

  call boundary_gs(b2vecini_vec%vec, feedfac, gs_matrix)
  call finalize(gs_matrix)

  call define_profiles

  ! determine if feedback will be used
  do_feedback = .false.
  do_feedback_x = .false.
  if(idevice .eq. -1) then
     do i=1, maxcoils
        if((gs_vertical_feedback(i) .ne. 0) .or. &
             (gs_radial_feedback(i) .ne. 0)) then
           do_feedback = .true.
        end if
        if((gs_vertical_feedback_x(i) .ne. 0) .or. &
             (gs_radial_feedback_x(i) .ne. 0)) then
           do_feedback_x = .true.
        end if
     end do
  end if
  if(do_feedback) then 
     if(xmag0.eq.0) then
        xmag0 = xmag
        zmag0 = zmag
     end if
     xmagi = 0.
     zmagi = 0.
  end if
  if(do_feedback_x) then
     if(xnull0.eq.0) then
        xnull0 = xnull
        znull0 = znull
     end if
     xnulli = 0.
     znulli = 0.
  end if

  if(igs.ne.0) call lcfs(psi_vec, iwall_is_limiter.eq.1, &
       igs_start_xpoint_search.eq.0)

  error2 = 0.
  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
  mainloop: do itnum=1, iabs(igs)

     if(myrank.eq.0) print *, "GS iteration = ", itnum, error2
     
     ! apply boundary conditions
     if(iread_eqdsk.ne.1 .or. itnum.gt.1) then
        feedfac = 0.
        if(itnum.gt.1 .and. gnorm.ne.0 .and. xlim2.ne.0 .and. igs_feedfac.eq.1) then
           feedfac = -0.25*(psilim - psilim2)/gnorm
           !......as a diagnostic, calculate the effective value of libetap (including feedback term)
           libetapeff =  libetapeff + feedfac/fac2
           if(myrank.eq.0 .and. iprint.ge.2) &
                write(*,'(A,2E12.4)') "feedfac, gnorm", feedfac,  gnorm
        endif

        if(myrank.eq.0 .and. iprint.ge.2) print *, '  applying bcs'
        call boundary_gs(b1vecini_vec%vec, feedfac)
        
        ! perform LU backsubstitution to get psi solution
        if(myrank.eq.0 .and. iprint.ge.2) print *, '  solving'

        if(myrank.eq.0 .and. iprint.ge.1) print *, 'before solve'
        call newsolve(gs_matrix,b1vecini_vec%vec,ier)
        if(myrank.eq.0 .and. iprint.ge.1) print *, ' after solve'

        if(ier.ne.0) then
           if(myrank.eq.0) print *, 'Error in GS solve'
           call safestop(10)
        end if
        if(is_nan(b1vecini_vec)) then 
           print *, 'Error: solution is NaN'
           call safestop(11)
        endif

        ! combine solve result with old solution to get new solution
        if(itnum.eq.1) then
           psi_vec = b1vecini_vec
           !call m3dc1_field_display(psi_vec%vec%id)
        else
           ! psi_vec = th_gs*b1vecini_vec + (1.-th_gs)*psi_vec
           b2vecini_vec = b1vecini_vec
           !call m3dc1_field_display(b2vecini_vec%vec%id)
           call mult(b2vecini_vec,th_gs)
           !call m3dc1_field_display(b2vecini_vec%vec%id)
           call mult(psi_vec,1.-th_gs)
           !call m3dc1_field_display(b2vecini_vec%vec%id)
           call add(psi_vec,b2vecini_vec)
           !call m3dc1_field_display(b2vecini_vec%vec%id)
           !print *, "finish"
        endif
     endif

     ! Find new magnetic axis and lcfs
     call lcfs(psi_vec,iwall_is_limiter.eq.1, &
          itnum.ge.igs_start_xpoint_search)
     if(psibound.eq.psimin) then
        if(myrank.eq.0) print *, 'ERROR: psimin = psilim = ', psibound
        call safestop(4)
     end if

     ! define the pressure and toroidal field functions
     if(constraint) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, '  calling fundef2'
        call fundef2(error3)
     else
        if(myrank.eq.0 .and. iprint.ge.1) print *, '  calling fundef'
        call fundef
     end if

     if(itnum.gt.1) then
        ! Calculate error in new solution
        if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating error'
        call calculate_error(error,error2,b1vecini_vec)
        if(constraint) error = error3

        if(myrank.eq.0 .and. iprint.ge.1) then
           write(*,'(A,1p4e12.4)') ' Error in GS solution: ', error, error2, xmag, zmag
        endif
        ! if error is NaN, quit if itnum > 2   (needed to run on KNL)
        if(error.ne.error) then
             if(itnum.le.2) then
                   error = 1.  
             else
                   call safestop(11)
             endif
        endif

        ! if error is sufficiently small, stop iterating
        if(itnum .gt. 1 .and. error2 .lt. tol_gs) exit mainloop
    endif
     ! calculate gammas to constrain current, etc.
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating gammas'
     call calculate_gamma(gamma2,gamma3,gamma4)
     if(myrank.eq.0 .and. iprint.ge.2) print *, "gamma2,gamma3,gamma4", gamma2,gamma3,gamma4
     ! do feedback
     if(do_feedback .or. do_feedback_x) call coil_feedback(itnum)

     ! Define RHS vector
     b2vecini_vec = fun1_vec
     if(gamma2.ne.0.) then
        b1vecini_vec = fun2_vec
        call mult(b1vecini_vec, gamma2)
        call add(b2vecini_vec, b1vecini_vec)
     endif
     if(gamma3.ne.0.) then
        b1vecini_vec = fun3_vec
        call mult(b1vecini_vec, gamma3)
        call add(b2vecini_vec, b1vecini_vec)
     endif
     if(gamma4.ne.0.) then
        b1vecini_vec = fun4_vec
        call mult(b1vecini_vec, gamma4)
        call add(b2vecini_vec, b1vecini_vec)
     endif
     call mult(b2vecini_vec, -1.)
     call matvecmult(mass_mat_rhs%mat, b2vecini_vec%vec, b1vecini_vec%vec)
    
  end do mainloop

  ! do feedback
  if(do_feedback) call write_feedback('current.dat.out')


  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, "Converged GS error =",error2
     if(feedfac.ne.0) print *, "init&fin(effective) libetap", libetap, libetapeff
  endif

  ! if igs is positive, stop after iabs(igs) iterations
  ! continue for igs negative
!  if(itnum.eq.igs) call safestop(3)

  ! recalculate lcfs
  if(igs.ne.0) call lcfs(psi_vec, iwall_is_limiter.eq.1)

  ! Define equilibrium fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' defining equilibrium fields'

! bateman scale for igs_pp_ffp_rescale .ne. 0
  if(igs_pp_ffp_rescale.ne.0 .and. batemanscale.ne.1.) then
     if(myrank.eq.0) bzero = bzero*batemanscale
     call mpi_bcast(bzero,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  endif

  ! solve for p and f fields that best approximate gs solution
  b1vecini_vec = 0.
  b2vecini_vec = 0.

if (ispradapt .eq. 1) then
  call create_field(b3vecini_vec, "b3vecini_vec")
  if(irot.ne.0) call create_field(b4vecini_vec, "b4vecini_vec")
else
  call create_field(b3vecini_vec)
  if(irot.ne.0) call create_field(b4vecini_vec)
endif

  !recalculate omega
!  if(irot.eq.-1) call calc_omega_profile

  if(myrank.eq.0 .and. iprint.ge.2) print *, '  populating'
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_main, int_tor)
     call define_fields(itri, 0, 1, 0)
     
     call eval_ops(itri, psi_vec, ps079)
     if(icsubtract.eq.1) then 
        call eval_ops(itri, psi_coil_field, psc79)
        ps079 = ps079 + psc79
     end if

     call get_zone(itri, izone)

     do i=1, npoints
        call calc_toroidal_field(ps079(i,:),tf,x_79(i),z_79(i),izone)
        temp79b(i) = tf
        call calc_pressure(ps079(i,:),tf,x_79(i),z_79(i),izone)
        temp79a(i) = tf
        call calc_density(ps079(i,:),tf,x_79(i),z_79(i),izone)
        temp79c(i) = tf
        if(irot.ne.0) then
           call calc_rotation(ps079(i,:),tf,x_79(i),z_79(i),izone)
           temp79d(i) = tf
        endif
     end do
     
     temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
     temp(:,2) = intx2(mu79(:,:,OP_1),temp79b)
     temp(:,3) = intx2(mu79(:,:,OP_1),temp79c)
     if(irot.ne.0) temp(:,4) = intx2(mu79(:,:,OP_1),temp79d)
     call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     call vector_insert_block(b2vecini_vec%vec,itri,1,temp(:,2),VEC_ADD)
     call vector_insert_block(b3vecini_vec%vec,itri,1,temp(:,3),VEC_ADD)
     if(irot.ne.0) then
        call vector_insert_block(b4vecini_vec%vec,itri,1,temp(:,4),VEC_ADD)
     endif
  end do

  if(myrank.eq.0 .and. iprint.ge.2) print *, '  solving...'

  call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
  !call mult(b1vecini_vec, 0.5)

  p_field(0) = b1vecini_vec

  call newvar_solve(b2vecini_vec%vec,mass_mat_lhs)
  bz_field(0) = b2vecini_vec

  call newvar_solve(b3vecini_vec%vec,mass_mat_lhs)
  den_field(0) = b3vecini_vec

  if(irot.ne.0) then
     call newvar_solve(b4vecini_vec%vec,mass_mat_lhs)
     vz_field(0) = b4vecini_vec
  endif

  call destroy_field(b3vecini_vec)
  if(irot.ne.0) call destroy_field(b4vecini_vec)
     
  psi_field(0) = psi_vec
  psi_field(1) = 0.

#ifdef USEPARTICLES
  ! Define rho field
  if(allocated(rho_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating rho...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)
        
        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then 
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if
        
        call get_zone(itri, izone)
        
        do i=1, npoints 
           call calc_rho(ps079(i,:),tf,x_79(i),z_79(i),izone)
           temp79a(i) = tf
        end do
        
        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do
     
     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     !call mult(b1vecini_vec,0.25)
     rho_field = b1vecini_vec
  end if
  ! Define nf field
  if(allocated(nf_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating nf...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)
        
        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then 
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if
        
        call get_zone(itri, izone)
        
        do i=1, npoints 
           call calc_fdensity(ps079(i,:),tf,x_79(i),z_79(i),izone)
           temp79a(i) = tf
        end do
        
        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do
     
     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     !call mult(b1vecini_vec,0.25)
     nf_field = b1vecini_vec
  end if
  ! Define tf field
  if(allocated(tf_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating nf...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)
        
        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then 
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if
        
        call get_zone(itri, izone)
        
        do i=1, npoints 
           call calc_ftemp(ps079(i,:),tf,x_79(i),z_79(i),izone)
           temp79a(i) = tf
        end do
        
        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do
     
     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     !call mult(b1vecini_vec,0.25)
     tf_field = b1vecini_vec
  end if
  ! Define pf field
  if(allocated(tf_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating pf...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)

        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if

        call get_zone(itri, izone)

        do i=1, npoints
           call calc_fdensity(ps079(i,:),tf,x_79(i),z_79(i),izone)
           call calc_ftemp(ps079(i,:),tf2,x_79(i),z_79(i),izone)
           if (fast_ion_dist==1) then
              temp79a(i) =tf*tf2* 1.6022e-12 / (b0_norm**2/(4.*pi*n0_norm))!rsae
           else
              crit_v=sqrt(2*tf2*1.6e-19/fast_ion_mass/m_p)
              max_v=sqrt(2*fast_ion_max_energy*1.6e-19/fast_ion_mass/m_p)
              norm_1=1./6.*(-crit_v**2*log(crit_v**2-crit_v*max_v+max_v**2)+2*crit_v**2*log(crit_v+max_v)-2*sqrt(3.)*&
                 crit_v**2*atan((2*max_v-crit_v)/(sqrt(3.)*crit_v))+3*max_v**2)
              norm_2=1./6.*(-crit_v**2*log(crit_v**2)+2*crit_v**2*log(crit_v)-2*sqrt(3.)*crit_v**2*atan(-1./sqrt(3.)))
              dnorm=1./3.*(log(max_v**3+crit_v**3)-log(crit_v**3))
              temp79a(i) =tf*(norm_1-norm_2)/dnorm*fast_ion_mass*m_p/3.0&
                 *1.e7 / (b0_norm**2/(4.*pi*n0_norm))!rsae
           endif
        end do

        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do

     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     !call mult(b1vecini_vec,0.5)
     pf_field = b1vecini_vec
  end if 
  ! Define nfi field
  if(allocated(nfi_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating nfi...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)
        
        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then 
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if
        
        call get_zone(itri, izone)
        
        do i=1, npoints 
           call calc_fidensity(ps079(i,:),tf,x_79(i),z_79(i),izone)
           temp79a(i) = tf
        end do
        
        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do
     
     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     !call mult(b1vecini_vec,0.25)
     nfi_field = b1vecini_vec
   end if
  ! Define tfi field
  if(allocated(tfi_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating tfi...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)
        
        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then 
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if
        
        call get_zone(itri, izone)
        
        do i=1, npoints 
           call calc_fitemp(ps079(i,:),tf,x_79(i),z_79(i),izone)
           temp79a(i) = tf
        end do
        
        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do
     
     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     !call mult(b1vecini_vec,0.25)
     tfi_field = b1vecini_vec
  end if
  ! Define pfi field
  if(allocated(tfi_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating pfi...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)

        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if

        call get_zone(itri, izone)

        do i=1, npoints
           call calc_fidensity(ps079(i,:),tf,x_79(i),z_79(i),izone)
           call calc_fitemp(ps079(i,:),tf2,x_79(i),z_79(i),izone)
           temp79a(i) =tf*tf2* 1.6022e-12 / (b0_norm**2/(4.*pi*n0_norm))!rsae
        end do

        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do

     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     !call mult(b1vecini_vec,0.5)
     pfi_field = b1vecini_vec
  end if 

  if ((kinetic.eq.1).and.(kinetic_fast_ion.eq.1)) then
     call mult(pf_field, -1.)
     call add(p_field(0), pf_field)
     call mult(pf_field, -1.)
  endif

  if ((kinetic.eq.1).and.(kinetic_thermal_ion.eq.1)) then
     call mult(pfi_field, -1.)
     call add(p_field(0), pfi_field)
     call mult(pfi_field, -1.)
  endif
#endif

  ! Define pe field
  if(allocated(te_spline%y)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating pe...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)
        
        call eval_ops(itri, psi_field(0), ps079)
        if(icsubtract.eq.1) then 
           call eval_ops(itri, psi_coil_field, psc79)
           ps079 = ps079 + psc79
        end if
        
        call get_zone(itri, izone)
        
        do i=1, npoints 
           call calc_electron_pressure(ps079(i,:),tf,x_79(i),z_79(i),izone)
           temp79a(i) = tf
        end do
        
        temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do
     
     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     pe_field(0) = b1vecini_vec
  else
     pe_field(0) = p_field(0)
     call mult(pe_field(0), pefac)
  end if

  ! Define Ti and Te field
  if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating Ti and Te...'
  b1vecini_vec = 0.
  b2vecini_vec = 0.
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_main, int_tor)
     call define_fields(itri, 0, 1, 0)
     
     call eval_ops(itri, p_field(0), p079)
     call eval_ops(itri, pe_field(0), pe079)
     call eval_ops(itri, den_field(0), n079)
     
     temp79a = (pe079(:,OP_1)) / (z_ion*n079(:,OP_1))
     temp79b = (p079(:,OP_1) - pe079(:,OP_1)) / n079(:,OP_1)
     
     temp(:,1) = intx2(mu79(:,:,OP_1),temp79a)
     temp(:,2) = intx2(mu79(:,:,OP_1),temp79b)
     call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     call vector_insert_block(b2vecini_vec%vec,itri,1,temp(:,2),VEC_ADD)
  end do

  call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
  te_field(0) = b1vecini_vec  
  call newvar_solve(b2vecini_vec%vec,mass_mat_lhs)
  ti_field(0) = b2vecini_vec

  call finalize(field0_vec)

  ! free memory
  call destroy_field(b1vecini_vec)
  call destroy_field(b2vecini_vec)
  call destroy_field(psi_vec)
  call destroy_field(fun1_vec)
  call destroy_field(fun2_vec)
  call destroy_field(fun3_vec)
  call destroy_field(fun4_vec)

  call destroy_mat(gs_matrix)

  call destroy_spline(p0_spline)
  call destroy_spline(g0_spline)
  call destroy_spline(n0_spline)
  call destroy_spline(ffprime_spline)
  call destroy_spline(pprime_spline)
  if(.not.constraint) then
     call destroy_spline(g2_spline)
     call destroy_spline(g3_spline)
  end if
  if(irot.ne.0) then
     call destroy_spline(alpha_spline)
     call destroy_spline(omega_spline)
  end if
  call destroy_spline(te_spline)

  ! calculate final error
  call calculate_gs_error(error)
  if(myrank.eq.0) print *, 'Final error in GS solution: ', error

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'done gradshafranov_solve.'

end subroutine gradshafranov_solve

subroutine calculate_error(error, error2, psinew)
  use basic
  use field
  use boundary_conditions
  use mesh_mod

  implicit none

  include 'mpif.h'

  real, intent(out) :: error, error2
  type(field_type), intent(in) :: psinew

  integer :: i,inode, numnodes, izone, izonedim, ier
  real :: sum, sum2, norm, norm2, normal(2), curv(3), x, phi, z, lhs, rhs
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: psi0, psi1, f1, f2, f3, f4
  real, dimension(5) :: temp1, temp2

  sum = 0.
  norm = 0.
  sum2 = 0.
  norm2 = 0.
  temp2 = 0

  numnodes = owned_nodes()
  do i=1,numnodes
     inode = nodes_owned(i)
     call boundary_node(inode,is_boundary,izone,izonedim,normal,curv,x,phi,z,&
          BOUND_DOMAIN)
     if(is_boundary) cycle

     call get_node_data(psi_vec, inode, psi0)
     call get_node_data(psinew, inode, psi1)
     call get_node_data(fun1_vec, inode, f1)
     call get_node_data(fun2_vec, inode, f2)
     call get_node_data(fun3_vec, inode, f3)
     call get_node_data(fun4_vec, inode, f4)

     lhs = (psi0(4) + psi0(6) - psi0(2)/x)/x
     rhs =  -(f1(1) + gamma2*f2(1) + gamma3*f3(1) + gamma4*f4(1))

     sum = sum + (lhs-rhs)**2
     norm = norm + lhs**2
     sum2 = sum2 + abs(psi0(1) - psi1(1))
     norm2 = norm2 + abs(psi0(1))
  enddo

  if(maxrank.gt.1) then
     temp1(1) = sum
     temp1(2) = norm
     temp1(3) = sum2
     temp1(4) = norm2
     call mpi_allreduce(temp1, temp2, 4, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ier)
     error = sqrt(temp2(1)/temp2(2))
     error2= temp2(3)/temp2(4)
  else
     error = sqrt(sum/norm)
     error2= sum2/norm2
  endif
end subroutine calculate_error

!============================================================
! calculate_gamma
! ~~~~~~~~~~~~~~~
! calculates the values of gamma2, gamma3, and gamma4 to
! constrain solution to have the specified current, etc.
!============================================================
subroutine calculate_gamma(g2, g3, g4)
  use basic
  use mesh_mod
  use arrays
  use m3dc1_nint
  use math

  implicit none

  include 'mpif.h'

  real, intent(out) :: g2, g3, g4

  integer :: itri, numelms, ier

  real :: gsint1, gsint2, gsint3, gsint4, curr, g0
  real, dimension(5) :: temp1, temp2

  vectype, dimension(MAX_PTS,OP_NUM) :: psi_n, psic_n, &
       fun1_n, fun2_n, fun3_n, fun4_n

  ! start of loop over triangles to compute integrals needed to keep
  !     total current and q_0 constant using gamma4, gamma2, gamma3
  if(nv1equ.eq.1) then
     g2 = 0.
     g3 = 0.
     g4 = 0.
     return
  endif
       
  if(constraint) then
     g2 = 0.
     g3 = 0.
     g4 = 1.
     return
  endif

  g0 = bzero*rzero

  gsint1 = 0.
  gsint4 = 0.
  gsint2 = 0.
  gsint3 = 0.
  curr = 0.

  numelms = local_elements()

  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_main, int_tor)
     call define_fields(itri, 0, 0, 0)

     call eval_ops(itri, psi_vec, psi_n)
     if(icsubtract.eq.1) then
        call eval_ops(itri, psi_coil_field, psic_n)
        psi_n = psi_n + psic_n
     end if
     call eval_ops(itri, fun1_vec, fun1_n)
     call eval_ops(itri, fun2_vec, fun2_n)
     call eval_ops(itri, fun3_vec, fun3_n)
     call eval_ops(itri, fun4_vec, fun4_n)

     curr = curr - int2(ri2_79,psi_n(:,OP_GS))
     gsint1 = gsint1 + int2(ri_79,fun1_n)
     gsint2 = gsint2 + int2(ri_79,fun2_n)
     gsint3 = gsint3 + int2(ri_79,fun3_n)
     gsint4 = gsint4 + int2(ri_79,fun4_n)
  enddo

  if(iprint.ge.2) write(80+myrank,1080) myrank,curr,gsint1,gsint2,gsint3,gsint4
 1080 format(i12,1p5e12.4)
     
  if(maxrank.gt.1) then
     temp1(1) = curr
     temp1(2) = gsint1
     temp1(3) = gsint2
     temp1(4) = gsint3
     temp1(5) = gsint4
     call mpi_allreduce(temp1, temp2, 5, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ier)
     curr   = temp2(1)
     gsint1 = temp2(2)
     gsint2 = temp2(3)
     gsint3 = temp2(4)
     gsint4 = temp2(5)
  end if
!
#ifdef USE3D
     gsint1 = gsint1/twopi
     gsint2 = gsint2/twopi
     gsint3 = gsint3/twopi
     gsint4 = gsint4/twopi
#endif

  ! choose gamma2 to fix q0/qstar.  Note that there is an additional
  ! degree of freedom in gamma3.  Could be used to fix qprime(0)

#ifdef USEST
  if (igeometry.eq.1) then ! use physical coordinates of magnetic axis
     g2 =  -xmagp**2*p0*p1 - 2.*abs(g0)/(xmagp*q0*abs(dpsii))
     g3 = -4.*(abs(g0)/xmagp)*djdpsi/dpsii - xmagp**2*p0*p2
  else
#endif  
  g2 =  -xmag**2*p0*p1 - 2.*abs(g0)/(xmag*q0*abs(dpsii))
  g3 = -4.*(abs(g0)/xmag)*djdpsi/dpsii - xmag**2*p0*p2
#ifdef USEST
  end if
#endif
  if(gsint4.eq.0.) then
     g4 = 0.
  else
     g4 = -(-tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4
  end if

  if(myrank.eq.0 .and. iprint.ge.2) write(79,1079) dpsii,curr,gsint1,gsint2,gsint3,gsint4
!
 1079 format("dpsii,curr,gsint1,gsint2,gsint3,gsint4",1p6e12.4)
  if(myrank.eq.0 .and. iprint.ge.1) write(*,'(A,1p1e12.4)') ' current = ', curr

end subroutine calculate_gamma


! ===========================================================
! deltafun
! ~~~~~~~~
! sets jout_i =  <mu_i | -val*delta(R-x)delta(Z-z)> 
! ===========================================================
subroutine deltafun(x,z,val,jout)

  use mesh_mod
  use basic
  use arrays
  use field
  use m3dc1_nint
  use math

  implicit none

  include 'mpif.h'

  type(element_data) :: d
  real, intent(in) :: x, z, val
  real :: val2
  type(field_type), intent(inout) :: jout

  integer :: itri, i, k,in_domain, in_domains, ier
  real :: x1, z1, si, zi, eta
  vectype, dimension(dofs_per_element) :: temp
  real, dimension(dofs_per_element,coeffs_per_element) :: c

  itri = 0
  call whattri(x, 0., z, itri, x1, z1, IGNORE_PHI)

  val2 = val
  if(maxrank.gt.1) then
     in_domain = 0
     if(itri.gt.0) in_domain = 1
     call mpi_allreduce(in_domain,in_domains,1,MPI_INTEGER, &
          MPI_SUM,MPI_COMM_WORLD,ier)
     call mpi_barrier(MPI_COMM_WORLD, ier)
     if(in_domains.eq.0) then
        if(myrank.eq.0) print *, 'Error: point not found:', x, z
        call safestop(4)
     end if
     val2 = val/in_domains
  end if

  if(itri.gt.0) then
     temp = 0.

     ! calculate local coordinates
     call get_element_data(itri, d)
     call global_to_local(d, x, 0., z, si, zi, eta)

     ! calculate temp_i = -val*mu_i(si,eta)
     if(iprecompute_metric.eq.1) then
        c = ctri(:,:,itri)
     else
        call local_coeff_vector(itri, c)
     endif 
     do i=1,dofs_per_element
        do k=1, coeffs_per_tri
           temp(i) = temp(i) - val2*c(i,k)*si**mi(k)*eta**ni(k)
        end do
        if(equilibrate.ne.0) temp(i) = temp(i)*equil_fac(i, itri)
     end do

#ifdef USE3D
     temp = temp*twopi
#endif
!     call vector_insert_block(jout%vec, itri, jout%index, temp, VEC_SET)
     call vector_insert_block(jout%vec, itri, jout%index, temp, VEC_ADD)
  end if
  call sum_shared(jout%vec)

end subroutine deltafun

! ===========================================================
subroutine gaussianfun(x,z,val,denom,jout)

  use mesh_mod
  use basic
  use arrays
  use field
  use m3dc1_nint
  use math

  implicit none

  include 'mpif.h'

  real, intent(in) :: x, z, val, denom
  type(field_type), intent(inout) :: jout

  integer :: itri, nelms
  real :: befo
  vectype, dimension(MAX_PTS) :: temp
  vectype, dimension(dofs_per_element) :: dofs
#ifdef USEST
  if (igeometry.eq.1) then ! use physical coordinates of magnetic axis
      call physical_geometry(xmagp,zmagp,x,0.,z)
      if(myrank.eq.0 .and. iprint.ge.2) print *, 'Physical magnetic axis:', xmagp, zmagp 
  end if
#endif

  befo =-val/(pi*denom**2)
  
  nelms = local_elements()
  if (myrank.eq.0 .and. iprint .ge.1) write(*,1001) x,z,val,denom,befo,nelms
1001 format('Gaussian Called:',/,'    x,z,val,denom,befo,nelms',1p5e12.4,i5)
  do itri=1,nelms
      call define_element_quadrature(itri,int_pts_diag, int_pts_tor)
      call define_fields(itri,0,1,0)  !  defines x_79,z_79,mu,nu
#ifdef USEST
      if (igeometry.eq.1) then ! use physical coordinates of magnetic axis
          temp = befo*exp(-((x_79-xmagp)**2 + (z_79-zmagp)**2)/denom**2)
      else
#endif 
!     assemble matrix
      temp = befo*exp(-((x_79-x)**2 + (z_79-z)**2)/denom**2)
#ifdef USEST
      end if
#endif 
      dofs = intx2(mu79(:,:,OP_1),temp)
      ! call vector_insert_block(jout%vec, itri, jout%index, dofs, VEC_ADD)
      call vector_insert_block(jout%vec, itri, 1, dofs, VEC_ADD)
  enddo
  call sum_shared(jout%vec)

end subroutine gaussianfun


!============================================================
subroutine fundef

!.....defines the source functions for the GS equation:
  ! fun1 = r*p'
  ! fun4 = G4'/r
  ! fun2 = G2'/r
  ! fun3 = G3'/r

  use arrays
  use basic
  use diagnostics
  use mesh_mod
  
  implicit none 
  integer :: inode,ii, numnodes
  real :: x, phi, z, pso, pspx, pspy, pspxx, pspxy, pspyy
  real :: fbig0, fbig, fbigp, fbigpp
  real :: g2big, g2bigp, g2bigpp
  real :: g3big, g3bigp, g3bigpp
  real :: g4big, g4bigp, g4bigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha
  integer :: izone, izonedim, iout, mr
#ifdef USEST
  real :: xp, zp ! physical coordinates 
#endif

  vectype, dimension(dofs_per_node) :: temp, temp2

  dpsii =  (1./(psibound - psimin))

  numnodes = owned_nodes()
  do inode=1,numnodes
     ii=nodes_owned(inode) 
     call get_node_pos(ii, x, phi, z)
#ifdef USEST
     if (igeometry.eq.1) then ! use physical coordinates 
        call physical_geometry(xp,zp,x,phi,z)
        x = xp
        z = zp
     end if
#endif
    
     call get_node_data(psi_vec, ii, temp)
     if(icsubtract.eq.1) then
        call get_node_data(psi_coil_field, ii, temp2)
        temp = temp + temp2
     end if
     pso =  (temp(1) - psimin)*dpsii

     call m3dc1_ent_getgeomclass(0,ii-1,izonedim,izone) 
     call magnetic_region(temp(1),temp(2),temp(3),x,z,mr)
     if(mr.ne.REGION_PLASMA .or. izone.ne.ZONE_PLASMA) then
        temp = 0.
        call set_node_data(fun1_vec, ii, temp)
        call set_node_data(fun2_vec, ii, temp)
        call set_node_data(fun3_vec, ii, temp)
        call set_node_data(fun4_vec, ii, temp)
     else
        ! derivatives of real flux
        pspx  = temp(2)
        pspy  = temp(3)
        pspxx = temp(4)
        pspxy = temp(5)
        pspyy = temp(6)
     
        call evaluate_spline(pprime_spline,pso,fbig,fbigp,fbigpp,iout=iout)
        if(iout.eq.1) then
           fbig = 0.
           fbigp = 0.
           fbigpp = 0.
        else
           ! convert from normalized to real flux
           fbig = fbig*dpsii
           fbigp = fbigp*dpsii**2
           fbigpp = fbigpp*dpsii**3
        end if

        if(irot.eq.1) then
           ! include toroidal rotation in equilibrium
           ! this section calculates pressure derivatives 
           ! in presence of rotation
           !   scj   11/19/10
           !
           ! assumes pressure of the form 
           ! p(x,psi) = p(psi) exp (alpha(psi)*r1)
           
           r0m = 1./rzero**2
           r1 = (x**2 - rzero**2)/rzero**2
           r1m= x**2/rzero**2
           r2 = (x**2 - rzero**2)**2/rzero**4
           r3 = (x**2 - rzero**2)**3/rzero**6

           call evaluate_spline(p0_spline,pso,fbig0)
           call evaluate_spline(alpha_spline,pso, &
                alphap0,alphap,alphapp,alphappp)
           ! convert from normalized to real flux
           alphap = alphap*dpsii
           alphapp = alphapp*dpsii**2
           alphappp = alphappp*dpsii**3

           ealpha = exp(alphap0*r1)

           temp(1) = ealpha*(x*fbig + fbig0*alphap*x*r1)

           temp(2) = ealpha*                                             &
                (fbig + fbig0*alphap*r1                                  &
                +fbig0*alphap0*alphap*2*r1m*r1                           &
                + (alphap0*fbig + fbig0*alphap) * 2.*r1m                 &
                + pspx*( x*fbigp + (2.*fbig*alphap + fbig0*alphapp)*x*r1 &
                + fbig0 * alphap**2 * x*r2))

           temp(3) = ealpha*                                             &
                pspy*(x*fbigp + (2.*fbig*alphap + fbig0*alphapp)*x*r1   &
                + fbig0 * alphap**2 * x*r2)

           temp(4) = ealpha*                                                &
                ((3*alphap0*fbig + 3*fbig0*alphap)*2*x*r0m                  &
                + 3*fbig0*alphap0*alphap*2*x*r0m*r1                         &
                + fbig0*alphap0**2*alphap*4*r1m*r0m*r1                      &
                + alphap0*(alphap0*fbig+2*fbig0*alphap)*4.*r1m*r0m          &
                + pspx*( 2*fbigp                                            &
                +(2*fbig0*alphapp + 4*alphap*fbig)*r1                       &
                +(4*alphap*fbig + 2*alphap0*fbigp + 2*fbig0*alphapp)*2*r1m  &
                +(4*alphap*fbig + 2*alphap0*fbigp + 2*fbig0*alphapp)*r1m    &
                +2*fbig0*alphap**2*r2                                       &
                +(4*fbig0*alphap**2 + 2*fbig0*alphap0*alphapp               &
                + 4*alphap0*fbig*alphap)*r1m*r1                             &
                + 2*fbig0*alphap0*alphap**2*2*r1m*r2)                       &
                + pspx*pspx*(x*fbigpp                                       &
                +(3*fbigp*alphap + 3*fbig*alphapp + fbig0*alphappp)*x*r1    &
                +3*(fbig0*alphapp + fbig*alphap)*alphap*x*r2                &
                + fbig0*alphap**3*x*r3)                                     &  
                + pspxx*(x*fbigp                                            &
                + (2*fbig*alphap + fbig0*alphapp)*x*r1                      &
                + p0*alphap**2*x*r2))

           temp(5) = ealpha*                                                &
                (pspy*(fbigp                                                &
                + (2.*fbig*alphap + fbig0*alphapp)*r1                       &
                + (2.*fbig*alphap + fbig0*alphapp + fbigp*alphap0)*2*r1m    &
                + (2*fbig*alphap*alphap0 + fbig0*alphapp*alphap0            &
                + 2*fbig0*alphap**2)*2*r1m*r1                               &
                +  fbig0*alphap**2*r2                                       &
                +  fbig0*alphap**2*alphap0*2*r1m*r2)                        &
                +pspx*pspy*(x*fbigpp                                        &
                + (3*fbigp*alphap + 3*fbig*alphapp + fbig0*alphappp)*x*r1   &
                + (3*fbig*alphap + 3*fbig0*alphapp)*alphap*x*r2             &
                + fbig0*alphap**3*x*r3)                                     &
                + pspxy*(x*fbigp                                            &
                +   (2*fbig*alphap + fbig0*alphapp)*x*r1                    &
                + p0*alphap**2*x*r2))

           temp(6) = ealpha*                                                &
                (pspy*pspy*(x*fbigpp                                        &
                + (3*fbigp*alphap + 3*fbig*alphapp + fbig0*alphappp)*x*r1   &
                + (3*fbig*alphap + 3*fbig0*alphapp)*alphap*x*r2             &
                + fbig0*alphap**3*x*r3)                                     &
                + pspyy*(x*fbigp                                            &
                + (2*fbig*alphap + fbig0*alphapp)*x*r1                      &
                + p0*alphap**2*x*r2))
        else        
           ! no toroidal rotation in equilibrium
           temp(1) = x*fbig
           temp(2) = fbig + x*fbigp*pspx
           temp(3) =        x*fbigp*pspy
           temp(4) = 2.*fbigp*pspx + x*(fbigpp*pspx**2+fbigp*pspxx)
           temp(5) = fbigp*pspy + x*(fbigpp*pspx*pspy +fbigp*pspxy)
           temp(6) = x*(fbigpp*pspy**2 + fbigp*pspyy)
        endif   !...end of branch on irot

        call set_node_data(fun1_vec, ii, temp)
        
        call evaluate_spline(ffprime_spline, pso, g4big, g4bigp, g4bigpp, iout=iout)
        if(iout.eq.1) then
           g4big = 0.
           g4bigp = 0.
           g4bigpp = 0.
        else
           ! convert from normalized to real flux
           g4big = g4big*dpsii
           g4bigp = g4bigp*dpsii**2
           g4bigpp = g4bigpp*dpsii**3
        end if
        
        temp(1) = g4big/x
        temp(2) = g4bigp*pspx/x - g4big/x**2
        temp(3) = g4bigp*pspy/x
        temp(4) = (g4bigpp*pspx**2 + g4bigp*pspxx)/x                  &
             - 2.*g4bigp*pspx/x**2 + 2.*g4big/x**3
        temp(5) = (g4bigpp*pspx*pspy+g4bigp*pspxy)/x                  &
             - g4bigp*pspy/x**2
        temp(6) = (g4bigpp*pspy**2 + g4bigp*pspyy)/x
        call set_node_data(fun4_vec, ii, temp)

        if(.not.constraint) then
           call evaluate_spline(g2_spline, pso, g2big, g2bigp, g2bigpp)
           call evaluate_spline(g3_spline, pso, g3big, g3bigp, g3bigpp)
           ! convert from normalized to real flux
           g2big = g2big*dpsii
           g2bigp = g2bigp*dpsii**2
           g2bigpp = g2bigpp*dpsii**3
           g3big = g3big*dpsii
           g3bigp = g3bigp*dpsii**2
           g3bigpp = g3bigpp*dpsii**3
        
           temp(1) = g2big/x
           temp(2) = g2bigp*pspx/x - g2big/x**2
           temp(3) = g2bigp*pspy/x
           temp(4) = (g2bigpp*pspx**2 + g2bigp*pspxx)/x                  &
                - 2.*g2bigp*pspx/x**2 + 2.*g2big/x**3
           temp(5) = (g2bigpp*pspx*pspy+g2bigp*pspxy)/x                  &
                - g2bigp*pspy/x**2
           temp(6) = (g2bigpp*pspy**2 + g2bigp*pspyy)/x
           call set_node_data(fun2_vec, ii, temp)
           
           temp(1) = g3big/x
           temp(2) = g3bigp*pspx/x - g3big/x**2
           temp(3) = g3bigp*pspy/x
           temp(4) = (g3bigpp*pspx**2 + g3bigp*pspxx)/x                  &
                - 2.*g3bigp*pspx/x**2 + 2.*g3big/x**3
           temp(5) = (g3bigpp*pspx*pspy+g3bigp*pspxy)/x                  &
                - g3bigp*pspy/x**2
           temp(6) = (g3bigpp*pspy**2 + g3bigp*pspyy)/x
           call set_node_data(fun3_vec, ii, temp)
        end if
     end if
  end do

  call finalize(fun1_vec%vec)
  call finalize(fun2_vec%vec)
  call finalize(fun3_vec%vec)
  call finalize(fun4_vec%vec)
end subroutine fundef


subroutine fundef2(error)

  use basic
  use mesh_mod
  use arrays
  use sparse
  use newvar_mod
  use m3dc1_nint
  use diagnostics

  implicit none

  include 'mpif.h'

  real, intent(out) :: error

  integer :: i, itri, numelms, ier, iout
  real :: pso, psm, dpsii, norm, temp1(2), temp2(2)
  vectype, dimension(dofs_per_element) :: temp3, temp4
      
  integer :: izone, mr
  real :: pp0, a0, ap, p, ffp0, f
  real :: psinb
!!$  real :: w0, wp, n0, np

  dpsii =  1./(psibound - psimin)

  fun1_vec = 0.
  fun2_vec = 0.
  fun3_vec = 0.
  fun4_vec = 0.
  norm = 0.
  error = 0.

  if(myrank.eq.0 .and. iprint.ge.1) print *, '  In fundef2'

  numelms = local_elements()
  do itri=1,numelms

     call get_zone(itri, izone)

!     if(izone.ne.1) then
     if(izone.ne.ZONE_PLASMA) then
        temp3 = 0.
        temp4 = 0.
        call vector_insert_block(fun1_vec%vec,itri,1,temp3,VEC_ADD)
        call vector_insert_block(fun4_vec%vec,itri,1,temp4,VEC_ADD)
        cycle
     end if

     call define_element_quadrature(itri, int_pts_main, int_tor)
     call define_fields(itri, 0, 1, 0)

     call eval_ops(itri, psi_vec, ps079)
     if(icsubtract.eq.1) then 
        call eval_ops(itri, psi_coil_field, psc79)
        ps079 = ps079 + psc79
     end if
     
     do i=1, npoints
        
        pso = (ps079(i,OP_1)-psimin)*dpsii
        psm = pso
        f = 1.
        ! if we are in private flux region, make sure Psi > 1
        call magnetic_region(ps079(i,OP_1),ps079(i,OP_DR),ps079(i,OP_DZ), &
             x_79(i),z_79(i),mr,psinb)
        
        if(igs_forcefree_lcfs.ge.1 .and. mr.ne.REGION_PLASMA) then
           ffp0 = 0.
           pp0 = 0.
           a0 = 0.
           ap = 0.
           p = 0.
        else
           if(mr.eq.REGION_PF) then
              psm = 1.
              pso = psin_in_pf_region(pso, psinb, f)
           else
              f = 1.
           end if

           if(mr.eq.REGION_PLASMA) then
              call evaluate_spline(ffprime_spline,psm,ffp0, iout=iout)
              if(iout.eq.1) ffp0 = 0.
           else
              ffp0 = 0.
           end if
           
           call evaluate_spline(pprime_spline,pso,pp0, iout=iout)
           if(iout.eq.1) pp0 = 0.
           if(irot.eq.1) then
              call evaluate_spline(alpha_spline,pso,a0,ap,iout=iout)
              if(iout.eq.1) ap = 0.
              call evaluate_spline(p0_spline,pso,p)
              
!!$                 call evaluate_spline(omega_spline,pso,w0,wp,iout=iout)
!!$                 if(iout.eq.1) wp = 0.
!!$                 call evaluate_spline(n0_spline,pso,n0,np,iout=iout)
!!$                 if(iout.eq.1) np = 0.
!!$                 a0 = 0.5*rzero**2*w0*w0*n0/p
!!$                 ap = 0.5*rzero**2*w0*(2.*wp*n0 + w0*np - w0*n0*pp0/p)/p
           end if
        end if
        temp79a(i) = pp0*f
        temp79b(i) = ffp0*f
        if(irot.eq.1) then
           temp79c(i) = a0
           temp79d(i) = ap*f
           temp79e(i) = p
        end if
     end do

     
     ! convert from normalized to real flux
     temp79a = temp79a*dpsii
     temp79b = temp79b*dpsii
     if(irot.eq.1) temp79d = temp79d*dpsii
     
     ! p(R, psi) = p0(psi)*Exp[ alpha(psi) * (R^2-R0^2)/R0^2 ]
     
     ! dp/dpsi = [ p0' + p0*alpha'*(R^2-R0^2)/R0^2 ] * Exp[ ... ]  
     
     if(irot.eq.1) then
        temp79a = exp(temp79c*(x_79**2 - rzero**2)/rzero**2)* &
             (temp79a + temp79e*temp79d*(x_79**2 - rzero**2)/rzero**2)
     endif
     
     temp3 = intx3(mu79(:,:,OP_1),r_79,temp79a)
     temp4 = intx3(mu79(:,:,OP_1),ri_79,temp79b)
     
     call vector_insert_block(fun1_vec%vec,itri,1,temp3,VEC_ADD)
     call vector_insert_block(fun4_vec%vec,itri,1,temp4,VEC_ADD)
     
     ! Del*[psi] + R^2 p' + FF' = 0 
     temp79c = ps079(:,OP_GS) + r2_79*temp79a + temp79b
     
     norm = norm + abs(int2(ri_79,ps079(:,OP_GS)))
     error = error + abs(int2(ri_79,temp79c))
  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, '    solving...'
  call newvar_solve(fun1_vec%vec, mass_mat_lhs)
  call newvar_solve(fun4_vec%vec, mass_mat_lhs)

  if(maxrank.gt.1) then
     temp1(1) = norm
     temp1(2) = error
     call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ier)
     norm     = temp2(1)
     error    = temp2(2)
  end if
  error = error / norm

  if(myrank.eq.0 .and. iprint.ge.1) print *, '   Done fundef2'
end subroutine fundef2

subroutine readpgfiles
  use basic

  implicit none

  integer :: j, n
  real :: dum
!  real :: y, yp, fac
  real :: dpsidpsin, ptot
  real, allocatable :: psinorm(:), pres0(:), g0(:), ffn(:), ppn(:)

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Reading profiles files"

  ! Read p and p' profiles
  open(unit=76,file="profiles-p",status="old")
  read(76,803) n
      if(myrank.eq.0 .and. iprint.ge.1) print *, 'number of pressure points',n    ! DEBUG
  allocate(psinorm(n), pres0(n), ppn(n))
  do j=1,n
     read(76,802) psinorm(j), pres0(j), ppn(j), dum, dum
  enddo
  close(76)

  p0 = pres0(1)
  call create_spline(p0_spline, n, psinorm, pres0)
  call create_spline(pprime_spline, n, psinorm, ppn)

  ! calculate dpsidpsin
  ptot = 0.
  do j=2, n
     ptot = ptot + (psinorm(j)-psinorm(j-1))*(ppn(j)+ppn(j-1))/2.
  end do
  dpsidpsin = (pres0(n) - pres0(1)) / ptot
  dpsii = 1./dpsidpsin
  if(myrank.eq.0 .and. iprint.ge.1) then
     print * , 'PTOT, DPSIDPSIN = ', pres0(n) - pres0(1), &
          dpsidpsin
  end if

  deallocate(psinorm, pres0, ppn)

  ! Read g and FF' profiles
  open(unit=77,file="profiles-g",status="old")
  read(77,804) n
  allocate(psinorm(n), g0(n), ffn(n))
  do j=1,n
    read(77,802) psinorm(j), g0(j), ffn(j), dum, dum
  enddo
  close(77)

  call create_spline(g0_spline, n, psinorm, g0)
  call create_spline(ffprime_spline, n, psinorm, ffn)
  deallocate(psinorm, g0, ffn)

  constraint = .true.
  
  ! Change p' and FF' profiles to use derivatives wrt. psi_norm
  pprime_spline%y = pprime_spline%y * dpsidpsin
  ffprime_spline%y = ffprime_spline%y * dpsidpsin

return
  802 format(5x,5e18.9)
  803 format(i5)
  804 format(i5)
end subroutine readpgfiles


 subroutine default_profiles
   use basic

   implicit none

   integer :: j, npsi
   real :: psii
   real, allocatable :: psinorm(:), pres0(:)
   real, allocatable :: g0(:), g2(:), g3(:), ffn(:), ppn(:)

   if(myrank.eq.0) print *, "Using analytic p' and FF' profiles"

   npsi = 100
   allocate(psinorm(npsi), pres0(npsi))
   allocate(g0(npsi), g2(npsi), g3(npsi), ffn(npsi), ppn(npsi))
  
   do j=1, npsi
      psii = (j-1.)/(npsi-1.)
      psinorm(j) = psii
      pres0(j) = p0*(1.+p1*psii+p2*psii**2 &
           - (20. + 10.*p1 + 4.*p2)*psii**3 &
           + (45. + 20.*p1 + 6.*p2)*psii**4 &
           - (36. + 15.*p1 + 4.*p2)*psii**5 &
           + (10. +  4.*p1 +    p2)*psii**6)
      ppn(j) = p0*(p1+2.*p2*psii &
           - 3.*(20. + 10.*p1 + 4.*p2)*psii**2 &
           + 4.*(45. + 20.*p1 + 6.*p2)*psii**3 &
           - 5.*(36. + 15.*p1 + 4.*p2)*psii**4 &
           + 6.*(10. +  4.*p1 +    p2)*psii**5)

      g0(j)  = 1.      - 20.*psii**3 + 45.*psii**4 - 36.*psii**5 + 10.*psii**6
      g2(j)  = 1.      - 30.*psii**2 + 80.*psii**3 - 75.*psii**4 + 24.*psii**5
      g3(j)  = 2.*psii - 12.*psii**2 + 24.*psii**3 - 20.*psii**4 +  6.*psii**5
      ffn(j) =         - 60.*psii**2 +180.*psii**3 -180.*psii**4 + 60.*psii**5
   end do

   call create_spline(p0_spline, npsi, psinorm, pres0)
   call create_spline(g0_spline, npsi, psinorm, g0)
   call create_spline(g2_spline, npsi, psinorm, g2)
   call create_spline(g3_spline, npsi, psinorm, g3)
   call create_spline(ffprime_spline, npsi, psinorm, ffn)
   call create_spline(pprime_spline, npsi, psinorm, ppn)

   deallocate(psinorm, pres0, g0, g2, g3, ffn, ppn)
 end subroutine default_profiles

 subroutine density_profile
   use basic

   call copy_spline(n0_spline, p0_spline)

   if(expn.eq.0.) then
      n0_spline%y = den0
   else
      if(idenfunc.ne.4 .and. den_edge.gt.0.) then
         n0_spline%y = den0*((p0_spline%y-p0_spline%y(p0_spline%n))/p0)**expn + den_edge
      else
         n0_spline%y = den0*(p0_spline%y/p0)**expn
      end if
   end if
 end subroutine density_profile

 !=======================================================================
 ! default_omega
 ! ~~~~~~~~~~~~~
 ! use analytic omega profile
 !=======================================================================
 subroutine default_omega
   use basic
   implicit none

   integer :: j, npsi
   real :: psii, pres, dens, alpha
   real, allocatable :: omega(:), psinorm(:)

   if(myrank.eq.0 .and. iprint.ge.1) print *, "Using analytic alpha profile"

   if(.not.allocated(p0_spline%y) .or. .not.allocated(n0_spline%y)) then
      print *, 'Error: p and n profiles not allocated'
      call safestop(53)
   end if

   npsi = 100
   allocate(omega(npsi), psinorm(npsi))

   do j=1, npsi
      psii = (j-1.)/(npsi-1.)
      psinorm(j) = psii

      alpha = alpha0 + alpha1*psii + alpha2*psii**2 + alpha3*psii**3 

      call evaluate_spline(p0_spline, psii, pres)
      call evaluate_spline(n0_spline, psii, dens)
      if(iscale_rot_by_p.eq.0) alpha = alpha * dens/pres
      if(iscale_rot_by_p.eq.2) then
        alpha = alpha0 + alpha1*exp(-((psii-alpha2)/alpha3)**2)
        alpha = alpha * dens/pres
      end if   
      omega(j) = sqrt(2./rzero**2 * alpha * pres/dens)
   end do
   if(iscale_rot_by_p.eq.1) omega = omega - omega(npsi)

   call create_spline(omega_spline, npsi, psinorm, omega)

   deallocate(omega, psinorm)

   if(myrank.eq.0 .and. iprint.ge.1) print *, "Done defining alpha"
 end subroutine default_omega

 !=======================================================================
 ! calculate_alpha
 ! ~~~~~~~~~~~~~~~
 ! calculate alpha from omega profile
 !=======================================================================
 subroutine calculate_alpha
   use basic
   implicit none

   integer :: j
   real :: psii, pres, dens, omega

   if(myrank.eq.0) print *, "Calculating alpha"

   ! make sure alpha spline is big enough to capture all variation in n and p
   if((omega_spline%x(omega_spline%n).ge.p0_spline%x(p0_spline%n)) .and. &
        (omega_spline%x(omega_spline%n).ge.n0_spline%x(n0_spline%n))) then 
      call copy_spline(alpha_spline, omega_spline)
   else if (p0_spline%x(p0_spline%n).ge.n0_spline%x(n0_spline%n)) then
      call copy_spline(alpha_spline, p0_spline)
   else 
      call copy_spline(alpha_spline, n0_spline)
   end if

   do j=1, alpha_spline%n
      psii = alpha_spline%x(j)

      call evaluate_spline(p0_spline, psii, pres)
      call evaluate_spline(n0_spline, psii, dens)
      call evaluate_spline(omega_spline, psii, omega)

      if(pres.le.0. .or. dens.le.0.) then
         print *, 'ERROR: negative pressure or density', pres, dens
         call safestop(24)
      end if
      alpha_spline%y(j) = 0.5*rzero**2*omega**2*dens/pres
   end do
   
 end subroutine calculate_alpha


!================================================================
! create_profile
! ~~~~~~~~~~~~~~
!
! Sets up the GS solver to use a specific profile for p and I.
! n = number of radial points in profile
! p = pressure profile
! pp = p' = dp/dpsi (with psi the non-normalized flux)
! f = I = R*B_phi
! ffp =  I*I'
! flux = psi
!================================================================
 subroutine create_profile(n, p, pp, f, ffp, flux)
   use basic
   use math

   implicit none

   integer, intent(in) :: n
   real, dimension(n), intent(in) :: p, pp, f, ffp, flux

   real, allocatable :: pres0(:), g0(:), psinorm(:), ffn(:), ppn(:)

   integer :: i
   real :: temp, scale

   allocate(psinorm(n), pres0(n), g0(n), ffn(n), ppn(n))

   pres0 = p
   g0    = 0.5*(f**2 - f(n)**2)
   ffn   = ffp*(flux(n) - flux(1))  ! convert to derivative wrt normalized psi
   ppn   = pp *(flux(n) - flux(1))  ! convert to derivative wrt normalized psi
   dpsii =  (1./(flux(n) - flux(1)))
   psinorm = (flux - flux(1)) * dpsii

   if(igs_pp_ffp_rescale.eq.1) then
      ! do sanity check on pprime, ffprime
     temp = 0.
     do i=1, n-1
        temp = temp - pp(i)*(flux(i+1)-flux(i))
     end do
     if(myrank.eq.0) print *, "int(p'), p = ", temp, p(1)
     scale = p(1)/temp
     if(myrank.eq.0) print *, "Scaling p' by", scale
     ppn = ppn*scale

     temp = (bzero*rzero)**2/2.
     do i=1, n-1
        temp = temp - ffp(i)*(flux(i+1)-flux(i))
     end do
     if(myrank.eq.0) print *, "int(ff'), ff/2 = ", temp, f(1)**2/2.
     scale = (f(1)**2/2.-f(n)**2/2.) / (temp-f(n)**2/2.)
     if(myrank.eq.0) print *, "Scaling FF' by", scale
     ffn = ffn*scale


   end if

   call create_spline(p0_spline, n, psinorm, pres0)
   call create_spline(g0_spline, n, psinorm, g0)
   call create_spline(ffprime_spline, n, psinorm, ffn)
   call create_spline(pprime_spline, n, psinorm, ppn)

   deallocate(psinorm, pres0, g0, ffn, ppn)

   constraint = .true.
 end subroutine create_profile


 !============================================================
 ! write_profile
 ! ~~~~~~~~~~~~~~
 ! Write profile data to file
 !============================================================
 subroutine write_profile
   implicit none

   integer :: j
   real :: y, yp, yp_p, ypp, yppp

   open(unit=75,file="testout",status="unknown")
   do j=1, 1000
      yppp = (j-1.)/(1000.-1.)
      call evaluate_spline(p0_spline, yppp, y, yp_p)
      call evaluate_spline(pprime_spline, yppp, yp, ypp)
      write(75,'(5e12.5)') yppp, y, yp_p, yp, ypp
   enddo
   close(75)

   open(unit=76,file="profilesdb-p",status="unknown")
   do j=1, p0_spline%n
      call evaluate_spline(p0_spline, p0_spline%x(j), y, yp, ypp, yppp)
      write(76,802) j, p0_spline%x(j), y, yp, ypp, yppp
   enddo
   close(76)
   
   open(unit=77,file="profilesdb-g",status="unknown")
   do j=1, g0_spline%n
      call evaluate_spline(g0_spline, g0_spline%x(j), y, yp, ypp, yppp)
      write(77,802) j, g0_spline%x(j), y, yp, ypp, yppp
   enddo
   close(77)

   open(unit=78,file="profilesdb-pprime",status="unknown")
   do j=1, pprime_spline%n
      call evaluate_spline(pprime_spline, pprime_spline%x(j), y, yp, ypp, yppp)
      write(78,802) j, pprime_spline%x(j), y, yp, ypp, yppp
   enddo
   close(78)
   
   open(unit=79,file="profilesdb-ffprime",status="unknown")
   do j=1, ffprime_spline%n
      call evaluate_spline(ffprime_spline, ffprime_spline%x(j), y, yp, ypp, yppp)
      write(79,802) j, ffprime_spline%x(j), y, yp, ypp, yppp
   enddo
   close(79)

   open(unit=80,file="profilesdb-w",status="unknown")
   do j=1, omega_spline%n
      call evaluate_spline(omega_spline, omega_spline%x(j), y, yp, ypp, yppp)
      write(80,802) j, omega_spline%x(j), y, yp, ypp, yppp
   enddo
   close(80)

802 format(i5,1p6e18.9)   
 end subroutine write_profile

 !=============================================================
 ! calculate_gs_error
 ! ~~~~~~~~~~~~~~~~~~
 ! calculates the fractional error in the GS solution
 !=============================================================
 subroutine calculate_gs_error(error)
   use basic
   use m3dc1_nint
   
   implicit none
   
   include 'mpif.h'
   
   real, intent(out) :: error
   
   real :: norm, temp1(2), temp2(2)
   integer :: itri, nelms, def_fields, ier, ieqs_temp, izone
   
   norm = 0.
   error = 0.
   
   ! must set eqsubtract=1 so that equilibrium field is read
   ieqs_temp = eqsubtract
   eqsubtract = 1

   def_fields = FIELD_PSI + FIELD_P + FIELD_I
   if(irot.eq.1) then
      def_fields = def_fields + FIELD_V + FIELD_N
   endif

   nelms = local_elements()
   do itri=1, nelms
      call define_element_quadrature(itri, int_pts_diag, int_tor)
      call define_fields(itri, def_fields, 0, 1)
      
      call get_zone(itri, izone)

      if(izone.ne.ZONE_PLASMA) cycle

      ! Del*[psi]<psi,psi> = -F<F,psi> - R^2<p,psi> + R^3 rho w^2 <R,psi>

      temp79a = ps079(:,OP_GS)*(ps079(:,OP_DR)**2 + ps079(:,OP_DZ)**2)
      temp79b = -r2_79* &
           (p079(:,OP_DR)*ps079(:,OP_DR) + p079(:,OP_DZ)*ps079(:,OP_DZ))
      temp79c = -bz079(:,OP_1)* &
           (bz079(:,OP_DR)*ps079(:,OP_DR) + bz079(:,OP_DZ)*ps079(:,OP_DZ))
      
      if(irot.eq.1) then
         temp79d = r3_79*n079(:,OP_1)*vz079(:,OP_1)**2*ps079(:,OP_DR)
      else 
         temp79d = 0.
      endif

      temp79e = temp79a - (temp79b + temp79c + temp79d)
         
      norm = norm + abs(int2(ri_79,temp79a))
      error = error + abs(int2(ri_79,temp79e))
   end do

   eqsubtract = ieqs_temp

   if(maxrank.gt.1) then
      temp1(1) = norm
      temp1(2) = error
      call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, &
           MPI_SUM, MPI_COMM_WORLD, ier)
      norm     = temp2(1)
      error    = temp2(2)
   end if
   if(norm.ne.0.) error = error / norm
 end subroutine calculate_gs_error

!======================================================================
! calc_toroidal_field
! ~~~~~~~~~~~~~~~~~~~
!
! calculates the toroidal field (I) as a function of the 
! normalized flux
!======================================================================
pure subroutine calc_toroidal_field(ps0,tf,x,z,izone)
  use basic
  use diagnostics

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: tf    ! toroidal field (I)
  real, intent(in) :: x, z
  integer, intent(in) :: izone
  
  real :: g2, g3, g4
  real :: psii     ! normalized flux
  integer :: mr
 
  call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr)
  if(mr.ne.REGION_PLASMA .or. izone.ne.ZONE_PLASMA) then
     tf = bzero*rzero
  else
     psii = (real(ps0(1)) - psimin)/(psibound - psimin)

     if(.not.constraint) then
        g2 = psii - 10.*psii**3 + 20.*psii**4       &
             - 15.*psii**5 + 4.*psii**6
     
        g3 = psii**2 - 4.*psii**3 + 6.*psii**4      &
             - 4.*psii**5 + psii**6
     else
        g2 = 0.
        g3 = 0.
     end if
     
     call evaluate_spline(g0_spline, psii, g4)
     
!.....convert from gg' = .5(g^2)' to (g^2)'
     g2 = 2.*g2
     g3 = 2.*g3
     g4 = 2.*g4

     tf = sqrt((bzero*rzero)**2 + gamma2*g2 + gamma3*g3 + gamma4*g4)
     if(bzero.lt.0) tf = -tf
  endif
end subroutine calc_toroidal_field


!======================================================================
! calc_pressure
! ~~~~~~~~~~~~~
!
! calculates the pressure as a function of the poloidal flux 
! (and major radius if rotation is present)
!======================================================================
subroutine calc_pressure(ps0, pres, x, z, izone)
  use basic
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: pres     ! pressure
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  real :: p
  real :: alpha
  real :: psii     ! normalized flux
  real :: psinb    ! normalized flux of private flux boundary
!!$  real :: n0, w0

  integer :: mr

  if(izone.ne.ZONE_PLASMA) then 
     pres = p0_spline%y(p0_spline%n)
     return
  end if

  psii = (real(ps0(1)) - psimin)/(psibound - psimin)

  call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)

  if(igs_forcefree_lcfs.ge.1 .and. mr.ne.REGION_PLASMA) then
     call evaluate_spline(p0_spline, 1., p)
     pres = p
     return
  end if

  ! if we are in private flux region, make sure Psi > 1
  if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii, psinb)

  call evaluate_spline(p0_spline, psii, p)

  if(irot.eq.1) then
     call evaluate_spline(alpha_spline, psii, alpha)
     pres = p*exp(alpha*(x**2 - rzero**2)/rzero**2)
!!$     call evaluate_spline(omega_spline,psii,w0)
!!$     call evaluate_spline(p0_spline,psii,p)
!!$     pres = p*exp(0.5*(x**2 - rzero**2)*w0**2*n0/p)
  else
     pres = p
  endif
end subroutine calc_pressure


!======================================================================
! calc_density
! ~~~~~~~~~~~~
!
! calculates the density as a function of the poloidal flux 
! (and major radius if rotation is present)
!======================================================================
subroutine calc_density(ps0, dens, x, z, izone)
  use basic
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: dens     ! density
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  real :: n0
  real :: alpha
  real :: psii     ! normalized flux
  real :: psinb
!!$  real :: p, w0

  integer :: mr

  if(izone.ne.ZONE_PLASMA) then
     dens = n0_spline%y(n0_spline%n)
     return
  end if

  call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)

  psii = (real(ps0(1)) - psimin)/(psibound - psimin)

  ! if we are in private flux region, make sure Psi > 1
  if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)

  call evaluate_spline(n0_spline, psii, n0)

  if(irot.eq.1) then
     call evaluate_spline(alpha_spline,psii,alpha)
     dens = n0*exp(alpha*(x**2 - rzero**2)/rzero**2)
!!$     call evaluate_spline(omega_spline,psii,w0)
!!$     call evaluate_spline(p0_spline,psii,p)
!!$     dens = n0*exp(0.5*(x**2 - rzero**2)*w0**2*n0/p)
  else
     dens = n0
  endif
end subroutine calc_density
!======================================================================
#ifdef USEPARTICLES
subroutine calc_rho(ps0, dens, x, z, izone)
  use basic
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: dens     ! density
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  real :: n0
  real :: alpha
  real :: psii     ! normalized flux
  real :: psinb
!!$  real :: p, w0

  integer :: mr

  if(izone.gt.2) then
     dens = rho_spline%y(rho_spline%n)
     return
  end if

  call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)

  psii = (real(ps0(1)) - psimin)/(psibound - psimin)

  ! if we are in private flux region, make sure Psi > 1
  if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)

  call evaluate_spline(rho_spline, psii, n0)

  dens = n0
end subroutine calc_rho

subroutine calc_fdensity(ps0, dens, x, z, izone)
  use basic
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: dens     ! density
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  real :: n0
  real :: alpha
  real :: psii     ! normalized flux
  real :: psinb
!!$  real :: p, w0

  integer :: mr

  if(izone.gt.2) then
     dens = nf_spline%y(nf_spline%n)
     return
  end if

  call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)

  psii = (real(ps0(1)) - psimin)/(psibound - psimin)

  ! if we are in private flux region, make sure Psi > 1
  if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)

  call evaluate_spline(nf_spline, psii, n0)

  dens = n0
end subroutine calc_fdensity

subroutine calc_ftemp(ps0, tf, x, z, izone)
  use basic
  use math
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: tf
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  vectype :: pres0, n0
  real :: psii          ! normalized flux
  real :: tf0, psinb
  integer :: mr

  if(allocated(tf_spline%y)) then
!     if(izone.ne.1) then 
     if(izone.ne.ZONE_PLASMA) then 
        tf0 = tf_spline%y(tf_spline%n)
     else
        psii = (real(ps0(1)) - psimin)/(psibound - psimin)
        call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)
        if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)
        call evaluate_spline(tf_spline,psii,tf0)
     end if
     tf = tf0
  else
     tf = 0.
  end if
end subroutine calc_ftemp
subroutine calc_fidensity(ps0, dens, x, z, izone)
  use basic
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: dens     ! density
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  real :: n0
  real :: alpha
  real :: psii     ! normalized flux
  real :: psinb
!!$  real :: p, w0

  integer :: mr

  if(izone.gt.2) then
     dens = nf_spline%y(nf_spline%n)
     return
  end if

  call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)

  psii = (real(ps0(1)) - psimin)/(psibound - psimin)

  ! if we are in private flux region, make sure Psi > 1
  if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)

  call evaluate_spline(nfi_spline, psii, n0)

  dens = n0
end subroutine calc_fidensity

subroutine calc_fitemp(ps0, tf, x, z, izone)
  use basic
  use math
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: tf
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  vectype :: pres0, n0
  real :: psii          ! normalized flux
  real :: tf0, psinb
  integer :: mr

  if(allocated(tfi_spline%y)) then
!     if(izone.ne.1) then 
     if(izone.ne.ZONE_PLASMA) then 
        tf0 = tfi_spline%y(tf_spline%n)
     else
        psii = (real(ps0(1)) - psimin)/(psibound - psimin)
        call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)
        if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)
        call evaluate_spline(tfi_spline,psii,tf0)
     end if
     tf = tf0
  else
     tf = 0.
  end if
end subroutine calc_fitemp
#endif
!======================================================================
! calc_electron_pressure
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the electron pressure as a function of the poloidal flux
!======================================================================
subroutine calc_electron_pressure(ps0, pe, x, z, izone)
  use basic
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  vectype, intent(out) :: pe
  real, intent(in) :: x, z
  integer, intent(in) :: izone

  vectype :: pres0, n0
  real :: psii          ! normalized flux
  real :: te0, psinb
  integer :: mr

  if(allocated(te_spline%y)) then
!     if(izone.ne.1) then 
     if(izone.ne.ZONE_PLASMA) then 
        te0 = te_spline%y(te_spline%n)
     else
        psii = (real(ps0(1)) - psimin)/(psibound - psimin)
        call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)
        if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)
        call evaluate_spline(te_spline,psii,te0)
     end if
     call calc_density(ps0, n0, x, z, izone)
     pe = z_ion*n0*te0    ! (ni = n0; ne = ni*Z_ion)
  else
     call calc_pressure(ps0, pres0, x, z, izone)
     pe = pres0*pefac
  end if
end subroutine calc_electron_pressure

!======================================================================
! calc_rotation
! ~~~~~~~~~~~~~
!
! calculates the rotation as a function of the poloidal flux
!======================================================================
subroutine calc_rotation(ps0,omega, x, z, izone)
  use basic
  use diagnostics

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: ps0
  real, intent(in) :: x, z
  vectype, intent(out) :: omega     ! rotation
  integer, intent(in) :: izone

  real :: psii     ! normalized flux
  real :: psinb
  real :: w0
  integer :: mr

  if(irot.eq.0) then
     omega = 0.
     return
  endif
  if(izone.ne.ZONE_PLASMA) then
     omega = omega_spline%y(omega_spline%n)
     return
  end if
 
  call magnetic_region(ps0(1),ps0(2),ps0(3),x,z,mr,psinb)

  if(igs_forcefree_lcfs.eq.1 .and. mr.ne.REGION_PLASMA) then
     omega = 0.
     return
  end if    

  psii = (real(ps0(1)) - psimin)/(psibound - psimin)
  if(mr.eq.REGION_PF) psii = psin_in_pf_region(psii,psinb)

  call evaluate_spline(omega_spline,psii,w0)

  omega = w0
end subroutine calc_rotation

!=======================================================
! boundary_gs
! ~~~~~~~~~~~
!
! sets boundary conditions on psi in the GS solver
!=======================================================
subroutine boundary_gs(rhs, feedfac, mat)
  use basic
  use arrays
  use mesh_mod
  use coils
  use vector_mod
  use matrix_mod
  use boundary_conditions

  implicit none
  
  real, intent(in) :: feedfac
  type(vector_type), intent(inout) :: rhs
  type(matrix_type), intent(inout), optional :: mat
    
  integer :: i, inode, izone, izonedim, index
  integer :: numnodes, ineg
  real :: normal(2), curv(3)
  real :: x, z, phi
  real, dimension(1,1,6) :: g
  real, dimension(1) :: xp, zp, xc, zc
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

#ifdef USE3D
  integer :: j
#endif

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.1) print *, "boundary_gs called"

  numnodes = owned_nodes()
  do i=1, numnodes
     inode = nodes_owned(i)
     index = node_index(rhs, inode, 1)

     call boundary_node(inode,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_DOMAIN)
     if(is_boundary) then

        ! add feedback field
        if(idevice .eq. 0 .and. ifixedb .eq. 0 .and. feedfac.ne.0.) then
           xp(1) = x
           zp(1) = z
           xc(1) = 102.
           zc(1) = 10.
           call gvect(xp,zp,1,xc,zc,1,g,1,ineg)
           temp(1:6) = g(1,1,1:6)*feedfac
           call add(psi_field(0), inode, temp)
        endif
        
        if(ifixedb.ge.1) then
           temp = 0.
        else
           call get_node_data(psi_field(0), inode, temp)
        endif
       
        call set_dirichlet_bc(index,rhs,temp,normal,curv,izonedim,mat)
     endif

#ifdef USE3D
     if(int_tor.eq.0) then
        do j=1, 6
           ! set toroidal derivatives to zero
           if(present(mat)) call identity_row(mat, index+j+5)
           call insert(rhs, index+j+5, 0., VEC_SET)
        end do
     end if
#endif
  end do
end subroutine boundary_gs

subroutine create_rho_from_q(npsi, psi, q)
  use basic
  use math
  implicit none

  integer, intent(in) :: npsi
  real, intent(in), dimension(npsi) :: psi, q

  integer :: i
  real :: x(npsi), y(npsi)

  x(1) = 0.
  y(1) = 0.
  do i=2, npsi
     x(i) = x(i-1) + 2.*pi*q(i) * (psi(i) - psi(i-1))   ! Phi
     y(i) = (psi(i) - psi(1))/(psi(npsi)-psi(1))  ! Psi
  end do

  ! convert from phi to rho
  x = sqrt(abs(x))

  ! normalize rho
  x = x/x(npsi)

  call create_spline(psi_spline, npsi, x, y)
end subroutine create_rho_from_q

subroutine rho_to_psi(n, rho, psi)
  use basic
  implicit none

  integer, intent(in) :: n
  real, intent(in) :: rho(n)
  real, intent(out) :: psi(n)

  integer :: i

  if(.not.allocated(psi_spline%y)) then
     if(myrank.eq.0) print *, 'Error: no rho data found'
     call safestop(16)
  end if

  do i=1, n
     call evaluate_spline(psi_spline, rho(i), psi(i))
  end do
end subroutine rho_to_psi

real function psin_in_pf_region(psin, psin0, dpsout)
  implicit none

  real, intent(in) :: psin                ! psi_norm
  real, intent(in) :: psin0               ! psi_norm of reversal surface
  real, intent(out), optional :: dpsout   ! dpsin_in_pf_region / dpsin

  real :: t, dt

  if(gs_pf_psi_width.eq.0.) then 
     t = 0.
     dt = 0.
  else
     t = tanh((psin - psin0)/gs_pf_psi_width)
     dt = 1. - t**2
  end if

  psin_in_pf_region = 2.*psin0 - psin + 2.*gs_pf_psi_width*t
  if(present(dpsout)) then 
     dpsout = -1. + 2.*dt
  end if
end function psin_in_pf_region

end module gradshafranov

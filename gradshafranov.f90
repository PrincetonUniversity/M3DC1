module gradshafranov

  use field
  use spline

  implicit none

  integer, parameter :: numvargs = 1

  type(spline1d), private :: omega_spline
  type(spline1d), private :: alpha_spline
  type(spline1d), private :: n0_spline
  type(spline1d), private :: g0_spline
  type(spline1d), private :: p0_spline
  type(spline1d), private :: g2_spline, g3_spline
  type(spline1d), private :: te_spline

  type(spline1d), private :: ffprime_spline
  type(spline1d), private :: pprime_spline

  type(field_type), private :: psi_vec
  type(field_type), private :: fun1_vec, fun2_vec, fun3_vec, fun4_vec

  real, private :: dpsii
  real, private :: gamma2, gamma3, gamma4  

  logical, private :: constraint = .false.

  real, private :: gnorm, libetapeff, fac2

  integer, private :: int_tor = 0

  ! if use_norm_psi==1, pprime and ffprime are derivs wrt normalized flux
  integer, private :: use_norm_psi = 1
contains

subroutine gradshafranov_init()

  use basic
  use arrays
  use diagnostics

  implicit none

  real :: tstart, tend

  ! Define initial values of psi
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ifixedb.eq.0) call vacuum_field
     
  ! define initial field associated with delta-function source
  !     corresponding to current tcuro at location (xmag,zmag)
  call deltafun(xmag,zmag,tcuro,jphi_field)

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

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  integer :: i, numnodes
  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask


  if(myrank.eq.0 .and. iprint.ge.1) print *, 'in gradshafranov_per'

  numnodes = owned_nodes()

  do i=1, numnodes

     call get_node_pos(i, x, phi, z)

     call get_local_vals(i)

     vmask = 1.
     vmask(1:6) = p0_l(1:6)/p0
     vmask(1) = vmask(1) - pedge/p0
     
     ! initial parallel rotation
     u1_l = phizero*vmask

     ! allow for initial toroidal rotation
     vz1_l = 0.
     if(vzero.ne.0) call add_angular_velocity(vz1_l, x+xzero, vzero*vmask)

     ! add random perturbations
     if(nonrect.eq.0) then
        vmask(1) = 1.
        vmask(2:6) = 0.
     endif
     call random_per(x,phi,z,vmask)

     call set_local_vals(i)
  enddo

  call finalize(field_vec)

end subroutine gradshafranov_per

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
  use coil_sets

  implicit none

  real, dimension(6) :: g1, g2
  real, dimension(maxcoils) :: xp, zp, xc, zc
  complex, dimension(maxcoils) :: ic
  real :: aminor, bv, rnorm, fac
  integer :: ineg, ipole, numcoils
  

  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " calculating vacuum field..."

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

  rnorm = 10.
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "gradshafranov_solve xmag zmag ", xmag, zmag
  
  !......define feedback parameters needed for normalization
  if(idevice .eq. 0) then
     xc(1) = 102.
     zc(1) = rnorm
     xp = xlim
     zp = zlim
     call gvect(xp,zp,xc,zc,1,g1,1,ineg)
     xp = xlim2
     zp = zlim2
     call gvect(xp,zp,xc,zc,1,g2,1,ineg)
     gnorm = g1(1) - g2(1)
  endif

  ipole = 0
  select case(idevice)
  case(-1)
     call load_coils(xc,zc,ic,numcoils,'coil.dat','current.dat',0)
     
  case(1) ! CDX-U
     if(myrank.eq.0) print *, "Using standard CDX-U configuration"
     numcoils = 4
     xc(1) = 0.846
     zc(1) = 0.360
     xc(2) = 0.846
     zc(2) =-0.360
     xc(3) = 0.381
     zc(3) = 0.802
     xc(4) = 0.381
     zc(4) =-0.802
     ic = -.2*fac
     
  case(2) ! NSTX
     if(myrank.eq.0) print *, "Using standard NSTX configuration"
     numcoils = NSTX_coils
     xc(1:NSTX_coils) = NSTX_r
     zc(1:NSTX_coils) = NSTX_z
     ic(1:NSTX_coils) = fac*NSTX_I

  case(3) ! ITER
     if(myrank.eq.0) print *, "Using standard ITER configuration"
     numcoils = ITER_coils
     xc(1:ITER_coils) = ITER_r
     zc(1:ITER_coils) = ITER_z
     ic(1:ITER_coils) = fac*ITER_I
     
  case(4) ! DIII
     if(myrank.eq.0) print *, "Using standard DIII-D configuration"
     numcoils = DIII_coils
     xc(1:DIII_coils) = DIII_r
     zc(1:DIII_coils) = DIII_z
     ic(1:DIII_coils) = fac*DIII_I
     
  case default ! Generic

     if(myrank.eq.0) print *, "Using generic (dipole) configuration"

     numcoils = 1
     xc(1) = 102.
     zc(1) = rnorm
     ipole = 1
     ic = bv*fac2
  end select


  ! Field due to coil currents
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating fields due to coils"
  call field_from_coils(xc,zc,ic,numcoils,psi_field(0),ipole)
 
  ! Field due to extra divertor currents
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating fields due to divertors"
  if(divertors.ge.1) then
     xc(1:2) = xdiv
     zc(1) = zdiv
     if(divertors.eq.2) zc(2) = -zdiv
     ic(1:2) = fac*divcur
     call field_from_coils(xc,zc,ic,divertors,psi_field(0),0)
  endif

  ! Field due to plasma current
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Calculating fields due to plasma"
  xc(1) = xmag
  zc(1) = zmag
  ic(1) = tcuro/(2.*pi)
  call field_from_coils(xc,zc,ic,1,psi_field(0),0)
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
  implicit none

  real, allocatable :: xvals(:), yvals(:)
  real :: teold
  integer :: nvals

  ! If p' and ff' profiles are not yet defined, define them
  if(.not.allocated(p0_spline%x)) then
     if(inumgs .eq. 1) then
        ! read p' and ff' profiles
        call readpgfiles
     else
        ! use analytic p' and ff' profiles
        call default_profiles
     end if
  end if

  ! ensure that derivatives in SOL are zero
  pprime_spline%y(pprime_spline%n) = 0.
  ffprime_spline%y(ffprime_spline%n) = 0.

  ! scale profiles
  p0_spline%y = p0_spline%y*pscale
  g0_spline%y = g0_spline%y*sqrt(bscale)
  bzero = bzero*bscale

  ! add pedge to pressure
  if(pedge.ge.0.) p0_spline%y = p0_spline%y - p0_spline%y(p0_spline%n) + pedge

  ! define Te profile
  ! ~~~~~~~~~~~~~~~~~
  select case(iread_te)
     
  case(1)
     ! Read in keV
     nvals = 0
     call read_ascii_column('profile_te', xvals, nvals, icol=1)
     call read_ascii_column('profile_te', yvals, nvals, icol=2)
     yvals = yvals * 1.6022e-9 / (b0_norm**2/(4.*pi*n0_norm))
  case(2)

  case default
     
  end select

  if(allocated(yvals)) then
     call create_spline(te_spline, nvals, xvals, yvals)
     deallocate(xvals, yvals)
  end if


  ! define density profile
  ! ~~~~~~~~~~~~~~~~~~~~~~
  select case(iread_ne)
     
  case(1)
     ! Read in 10^20/m^3
     nvals = 0
     call read_ascii_column('profile_ne', xvals, nvals, icol=1)
     call read_ascii_column('profile_ne', yvals, nvals, icol=2)
     yvals = yvals * 1e14 / n0_norm
     
  case(2)
     ! Read in 10^19/m^3
     call read_ascii_column('dne.xy', xvals, nvals, skip=3, icol=1)
     call read_ascii_column('dne.xy', yvals, nvals, skip=3, icol=7)
     yvals = yvals * 1e13 / n0_norm

  case default
     call density_profile
     
  end select

  if(allocated(yvals)) then
     call create_spline(n0_spline, nvals, xvals, yvals)
     deallocate(xvals, yvals)
  end if


  ! add tedge to temperature
  if(tedge.ge.0.) then
     if(allocated(te_spline%y)) then
        teold = te_spline%y(te_spline%n)
        te_spline%y = te_spline%y - teold + tedge
     else
        teold = pefac*p0_spline%y(p0_spline%n)/n0_spline%y(n0_spline%n)
     end if
     ! add difference to pressure profile, so ion temp is not affected.
     if(pedge.lt.0.) then
        p0_spline%y = p0_spline%y + n0_spline%n*(tedge - teold)
     end if
  end if


  ! define rotation profile
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(irot.eq.1) then
     select case(iread_omega)

     case(1)
        ! Read in krad/sec
        nvals = 0
        call read_ascii_column('profile_omega', xvals, nvals, icol=1)
        call read_ascii_column('profile_omega', yvals, nvals, icol=2)
        yvals = 1000.* yvals / (b0_norm/sqrt(4.*pi*1.6726e-24*n0_norm)/l0_norm)
        
     case(2)
        ! Read in rad/sec
        call read_ascii_column('dtrot.xy', xvals, nvals, skip=3, icol=1)
        call read_ascii_column('dtrot.xy', yvals, nvals, skip=3, icol=7)
        yvals = yvals / (b0_norm/sqrt(4.*pi*1.6726e-24*n0_norm)/l0_norm)

     case default
        call default_omega
     end select

     if(allocated(yvals)) then
        call create_spline(omega_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if

    ! calculate alpha from omega
     call calculate_alpha
  end if

  ! output profiles
  if(myrank.eq.0) call write_profile

end subroutine define_profiles

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

  implicit none

  type(element_data) :: d
  type(field_type) :: b1vecini_vec, b2vecini_vec
  type(field_type) :: b3vecini_vec, b4vecini_vec

  type(matrix_type) :: gs_matrix
  type(newvar_matrix) :: dp_mat_lhs

  integer :: itri,i,j,ier, itnum, ibound
  integer :: numelms, numnodes
  real :: feedfac

  real :: x, phi, z, error, error2, error3 

  integer, dimension(dofs_per_element) :: imask

  vectype, dimension(dofs_per_node) :: tf
  vectype, dimension(coeffs_per_element) :: avec
  vectype, dimension(dofs_per_element,dofs_per_element) :: temp

  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, "Calculating Grad-Shafranov Equilibrium"

  numnodes = owned_nodes()
  numelms = local_elements()

  ! allocate memory for arrays
  call create_field(b1vecini_vec)
  call create_field(b2vecini_vec)
  call create_field(psi_vec)
  call create_field(fun1_vec)
  call create_field(fun2_vec)
  call create_field(fun3_vec)
  call create_field(fun4_vec)

  psi_vec = psi_field(0)
  b1vecini_vec = jphi_field
  if(iread_eqdsk .eq. 1) psilim = psibound

  ! form the grad-sharfranov matrix
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " forming the GS matrix..."

  call set_matrix_index(gs_matrix, gsmatrix_sm)
  call create_mat(gs_matrix, numvargs, numvargs, icomplex, .true.)
#ifdef CJ_MATRIX_DUMP
  print *, "create_mat gradshafranov gs_matrix", gs_matrix%imatrix     
#endif 

  if(int_tor.eq.0) then
     ibound = BOUNDARY_DIRICHLET + BOUNDARY_AXISYMMETRIC
  else
     ibound = BOUNDARY_DIRICHLET
  endif

  ! populate the matrix
  do itri=1,numelms

     call define_element_quadrature(itri,int_pts_main,int_tor)
     call define_fields(itri,0,1,0)

     call get_boundary_mask(itri, ibound, imask)
     do i=1,dofs_per_element
        if(imask(i).eq.0) then
           temp(i,:) = 0.
        else
           do j=1,dofs_per_element
              temp(i,j) = int3(ri_79,mu79(:,OP_1,i),nu79(:,OP_GS,j))
           enddo
        endif
     enddo

     call insert_block(gs_matrix, itri, 1, 1, temp, MAT_ADD)
  enddo

  feedfac = 0.

  ! insert boundary conditions
  call flush(gs_matrix)

  call boundary_gs(b2vecini_vec%vec, feedfac, gs_matrix)
  call finalize(gs_matrix)

  call define_profiles

  error2 = 0.
  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
  mainloop: do itnum=1, iabs(igs)

     if(myrank.eq.0) print *, "GS iteration = ", itnum, error2
     
     ! apply boundary conditions
     if(iread_eqdsk.ne.1 .or. itnum.gt.1) then
        feedfac = 0.
        if(itnum.gt.1 .and. gnorm.ne.0 .and. xlim2.ne.0) then
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

#ifdef CJ_MATRIX_DUMP
        if(itnum.le.1) then 
           call write_matrix(gs_matrix,'gs_matrix')
           call write_vector(b1vecini_vec%vec, 'gs_matrix_rhs.out')
        endif
#endif 
        call newsolve(gs_matrix,b1vecini_vec%vec,ier)
        if(myrank.eq.0 .and. iprint.ge.2) print *, '  done solve'


#ifdef CJ_MATRIX_DUMP
        if(itnum.le.1) then 
           call write_vector(b1vecini_vec%vec, 'gs_matrix_sol.out')
        endif
#endif 

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
        else
           ! psi_vec = th_gs*b1vecini_vec + (1.-th_gs)*psi_vec
           b2vecini_vec = b1vecini_vec
           call mult(b2vecini_vec,th_gs)
           call mult(psi_vec,1.-th_gs)
           call add(psi_vec,b2vecini_vec)
        endif
     endif

     ! Find new magnetic axis and lcfs
     call lcfs(psi_vec)

     ! define the pressure and toroidal field functions
     if(constraint .and. igs_method.ne.1) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, '  calling fundef2'
        call fundef2(error3)
     else
        if(myrank.eq.0 .and. iprint.ge.2) print *, '  calling fundef'
        call fundef
     end if

     ! Calculate error in new solution
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating error'
     call calculate_error(error,error2,b1vecini_vec)
     if(constraint .and. igs_method.ne.1) error = error3

     if(myrank.eq.0 .and. iprint.ge.1) then
        write(*,'(A,1p2e12.4)') ' Error in GS solution: ', error, error2
     endif

     ! if error is sufficiently small, stop iterating
     if(itnum .gt. 1 .and. error2 .lt. tol_gs) exit mainloop
    
     ! calculate gammas to constrain current, etc.
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating gammas'
     call calculate_gamma(gamma2,gamma3,gamma4)

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

  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, "Converged GS error =",error2
     print *, "initial and final(effective) libetap", libetap, libetapeff
  endif

  ! if igs is positive, stop after iabs(igs) iterations
  ! continue for igs negative
!  if(itnum.eq.igs) call safestop(3)


  ! Define equilibrium fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' defining equilibrium fields'
  if(igs_method.eq.3) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  creating solution matrix'
     call set_matrix_index(dp_mat_lhs%mat, dp_mat_lhs_index)
     call create_newvar_matrix(dp_mat_lhs, NV_DCBOUND, NV_DP_MATRIX, .true.)
#ifdef CJ_MATRIX_DUMP
     print *, "create_mat gradshafranov dp_mat_lhs", dp_mat_lhs%mat%imatrix     
#endif 
  endif
  if(igs_method.eq.2 .or. igs_method.eq.3) then
     ! solve for p and f fields which best approximate gs solution
     b1vecini_vec = 0.
     b2vecini_vec = 0.

     call create_field(b3vecini_vec)
     if(irot.eq.1) call create_field(b4vecini_vec)

     if(myrank.eq.0 .and. iprint.ge.2) print *, '  populating'
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)

        call get_element_data(itri, d)
        call calcavector(itri, psi_vec, avec)
        call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
             npoints, ps079)

        if(igs_method.eq.2) then 
           do i=1, npoints       
              call calc_toroidal_field(ps079(i,:),tf,x_79(i),z_79(i))
              temp79b(i) = tf(1)
              call calc_pressure(ps079(i,:),tf,x_79(i),z_79(i))
              temp79a(i) = tf(1)
              call calc_density(ps079(i,:),tf,x_79(i),z_79(i))
              temp79c(i) = tf(1)
              if(irot.eq.1) then
                 call calc_rotation(ps079(i,:),tf,x_79(i),z_79(i))
                 temp79d(i) = tf(1)
              endif
           end do

           do i=1,dofs_per_element
              temp(i,1) = int2(mu79(:,OP_1,i),temp79a)
              temp(i,2) = int2(mu79(:,OP_1,i),temp79b)
              temp(i,3) = int2(mu79(:,OP_1,i),temp79c)
              if(irot.eq.1) temp(i,4) = int2(mu79(:,OP_1,i),temp79d)
           end do
           call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
           call vector_insert_block(b2vecini_vec%vec,itri,1,temp(:,2),VEC_ADD)
           call vector_insert_block(b3vecini_vec%vec,itri,1,temp(:,3),VEC_ADD)
           if(irot.eq.1) then
              call vector_insert_block(b4vecini_vec%vec,itri,1,temp(:,4),VEC_ADD)
           endif
        else if(igs_method.eq.3) then

           call calcavector(itri, fun1_vec, avec)
           call eval_ops(avec,xi_79,zi_79,eta_79, d%co, d%sn, ri_79, &
                npoints,ph079)
           call calcavector(itri, fun4_vec, avec)
           call eval_ops(avec,xi_79,zi_79,eta_79, d%co, d%sn, ri_79, &
                npoints,vz079)

           do i=1, dofs_per_element
              temp(i,1) = &
                   int4(ri_79,ph079(:,OP_1),mu79(:,OP_DR,i),ps079(:,OP_DR)) + &
                   int4(ri_79,ph079(:,OP_1),mu79(:,OP_DZ,i),ps079(:,OP_DZ))
              temp(i,2) = 2.*gamma4* &
                   (int4(r_79,vz079(:,OP_1),mu79(:,OP_DR,i),ps079(:,OP_DR)) &
                   +int4(r_79,vz079(:,OP_1),mu79(:,OP_DZ,i),ps079(:,OP_DZ)))
           end do
           call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
           call vector_insert_block(b2vecini_vec%vec,itri,1,temp(:,2),VEC_ADD)
        endif
     end do

     if(myrank.eq.0 .and. iprint.ge.2) print *, '  solving...'
     if(igs_method.eq.2) then
        call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
        p_field(0) = b1vecini_vec

        call newvar_solve(b2vecini_vec%vec,mass_mat_lhs)
        bz_field(0) = b2vecini_vec

        call newvar_solve(b3vecini_vec%vec,mass_mat_lhs)
        den_field(0) = b3vecini_vec

        if(irot.eq.1) then
           call newvar_solve(b4vecini_vec%vec,mass_mat_lhs)
           vz_field(0) = b4vecini_vec
        endif

     else if(igs_method.eq.3) then
        call newvar_solve(b1vecini_vec%vec,dp_mat_lhs)
        call add(b1vecini_vec, pedge)
        p_field(0) = b1vecini_vec

        b1vecini_vec = (bzero*rzero)**2
        call newvar_solve(b2vecini_vec%vec,dp_mat_lhs,b1vecini_vec%vec)

        call pow(b2vecini_vec, 0.5)
        if(bzero*rzero .lt. 0.) call mult(b2vecini_vec, -1.)
        bz_field(0) = b2vecini_vec
     endif

     call destroy_field(b3vecini_vec)
     if(irot.eq.1) call destroy_field(b4vecini_vec)
  endif
  if(igs_method.eq.3) call destroy_mat(dp_mat_lhs%mat)
     
  ! calculate density profile for igs_method.eq.3
  if(igs_method.eq.3) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, 'calculating density...'
     b1vecini_vec = 0.
     do itri=1,numelms
        call define_element_quadrature(itri, int_pts_main, int_tor)
        call define_fields(itri, 0, 1, 0)
        call get_element_data(itri, d)
        
        call calcavector(itri, psi_field(0), avec)
        call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, npoints, &
             ps079)

        do i=1, npoints       
           call calc_density(ps079(i,:),tf,x_79(i),z_79(i))
           temp79c(i) = tf(1)
        end do

        do i=1, dofs_per_element
           temp(i,1) =  int2(mu79(:,OP_1,i),temp79c)
        end do
        call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
     end do

     call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
     den_field(0) = b1vecini_vec
  endif


  psi_field(0) = psi_vec
  psi_field(1) = 0.

  if(igs_method.eq.1) then
     do i=1,numnodes
        call get_node_pos(i, x, phi, z)
        call get_local_vals(i)
        call calc_toroidal_field(psi0_l, bz0_l, x, z)
        call calc_pressure(psi0_l, p0_l, x, z)
        call calc_density(psi0_l,den0_l,x,z)
        call calc_rotation(psi0_l,vz0_l,x,z)
        call calc_electron_pressure(psi0_l, pe0_l, x, z)
        call set_local_vals(i)
     end do
  end if

  ! Define pe field
  if(igs_method.ne.1) then
     if(allocated(te_spline%y)) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating Te...'
        b1vecini_vec = 0.
        do itri=1,numelms
           call define_element_quadrature(itri, int_pts_main, int_tor)
           call define_fields(itri, 0, 1, 0)
           call get_element_data(itri, d)
           
           call calcavector(itri, psi_field(0), avec)
           call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
                npoints, ps079)

           do i=1, npoints 
              call calc_electron_pressure(ps079(i,:),tf,x_79(i),z_79(i))
              temp79a(i) = tf(1)
           end do

           do i=1, dofs_per_element
              temp(i,1) = int2(mu79(:,OP_1,i),temp79a)
           end do
           call vector_insert_block(b1vecini_vec%vec,itri,1,temp(:,1),VEC_ADD)
        end do

        call newvar_solve(b1vecini_vec%vec,mass_mat_lhs)
        pe_field(0) = b1vecini_vec
     else
        pe_field(0) = p_field(0)
        call mult(pe_field(0), pefac)
     end if
  end if

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
  if(irot.eq.1) then
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

  integer :: i, numnodes, izone, izonedim, ier
  real :: sum, sum2, norm, norm2, normal(2), curv, x, phi, z, lhs, rhs
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: psi0, psi1, f1, f2, f3, f4
  real, dimension(5) :: temp1, temp2

  sum = 0.
  norm = 0.
  sum2 = 0.
  norm2 = 0.

  numnodes = owned_nodes()
  do i=1,numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(is_boundary) cycle

     call get_node_pos(i,x,phi,z)

     call get_node_data(psi_vec, i, psi0)
     call get_node_data(psinew, i, psi1)
     call get_node_data(fun1_vec, i, f1)
     call get_node_data(fun2_vec, i, f2)
     call get_node_data(fun3_vec, i, f3)
     call get_node_data(fun4_vec, i, f4)

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

  type(element_data) :: d
  integer :: itri, numelms, ier

  real :: gsint1, gsint2, gsint3, gsint4, curr, g0
  real, dimension(5) :: temp1, temp2

  vectype, dimension(MAX_PTS,OP_NUM) :: psi_n, fun1_n, fun2_n, fun3_n, fun4_n
  vectype, dimension(coeffs_per_element) :: avec

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
     call get_element_data(itri, d)

     call calcavector(itri, psi_vec, avec)
     call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
          npoints, psi_n)
     call calcavector(itri, fun1_vec, avec)
     call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
          npoints, fun1_n)
     call calcavector(itri, fun2_vec, avec)
     call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
          npoints, fun2_n)
     call calcavector(itri, fun3_vec, avec)
     call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
          npoints, fun3_n)
     call calcavector(itri, fun4_vec, avec)
     call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
          npoints, fun4_n)

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
  
  g2 =  -xmag**2*p0*p1 - 2.*abs(g0)/(xmag*q0*abs(dpsii))
  g3 = -4.*(abs(g0)/xmag)*djdpsi/dpsii - xmag**2*p0*p2
  g4 = -(-tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4

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
  type(field_type), intent(out) :: jout

  integer :: itri, i, k,in_domain, in_domains, ier
  real :: x1, z1, si, zi, eta
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,coeffs_per_element) :: c

  itri = 0
  call whattri(x, 0., z, itri, x1, z1, IGNORE_PHI)

  val2 = val
  if(maxrank.gt.1) then
     in_domain = 0
     if(itri.gt.0) in_domain = 1
     call mpi_allreduce(in_domain,in_domains,1,MPI_INTEGER, &
          MPI_SUM,MPI_COMM_WORLD,ier)
     val2 = val/in_domains
  end if

  if(itri.gt.0) then
     temp = 0.

     ! calculate local coordinates
     call get_element_data(itri, d)
     call global_to_local(d, x, 0., z, si, zi, eta)

     ! calculate temp_i = -val*mu_i(si,eta)
     call local_coeff_vector(itri, c, .false.)
     do i=1,dofs_per_element
        do k=1, coeffs_per_tri
           temp(i) = temp(i) - val2*c(i,k)*si**mi(k)*eta**ni(k)
        end do
        if(equilibrate.eq.1) temp(i) = temp(i)*equil_fac(i, itri)
     end do

#ifdef USE3D
     temp = temp*twopi
#endif

     call vector_insert_block(jout%vec, itri, jout%index, temp, VEC_SET)
  end if

  call sum_shared(jout%vec)
end subroutine deltafun


!============================================================
subroutine fundef

!.....defines the source functions for the GS equation:
  ! fun1 = r*p'
  ! fun4 = G4'/r
  ! fun2 = G2'/r
  ! fun3 = G3'/r

  use basic
  use diagnostics
  use mesh_mod
  
  implicit none 
  integer :: inode, numnodes
  real :: x, phi, z, pso, pspx, pspy, pspxx, pspxy, pspyy
  real :: fbig0, fbig, fbigp, fbigpp
  real :: g2big, g2bigp, g2bigpp
  real :: g3big, g3bigp, g3bigpp
  real :: g4big, g4bigp, g4bigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha
  logical :: inside_lcfs

  vectype, dimension(dofs_per_node) :: temp

  dpsii = 1./(psilim - psimin)

  numnodes = owned_nodes()
  do inode=1,numnodes
     
     call get_node_pos(inode, x, phi, z)
     
     call get_node_data(psi_vec, inode, temp)
     pso =  (temp(1) - psimin)*dpsii
     
     if(.not.inside_lcfs(temp,x,z,.true.)) then
        temp = 0.
        call set_node_data(fun1_vec, inode, temp)
        call set_node_data(fun2_vec, inode, temp)
        call set_node_data(fun3_vec, inode, temp)
        call set_node_data(fun4_vec, inode, temp)
     else
        ! derivatives of real flux
        pspx  = temp(2)
        pspy  = temp(3)
        pspxx = temp(4)
        pspxy = temp(5)
        pspyy = temp(6)
     
        call evaluate_spline(pprime_spline,pso,fbig,fbigp,fbigpp)
        ! convert from normalized to real flux
        fbig = fbig*dpsii**(use_norm_psi)
        fbigp = fbigp*dpsii**(use_norm_psi+1)
        fbigpp = fbigpp*dpsii**(use_norm_psi+2)

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

        call set_node_data(fun1_vec, inode, temp)
        
        call evaluate_spline(ffprime_spline, pso, g4big, g4bigp, g4bigpp)
        ! convert from normalized to real flux
        g4big = g4big*dpsii**use_norm_psi
        g4bigp = g4bigp*dpsii**(use_norm_psi+1)
        g4bigpp = g4bigpp*dpsii**(use_norm_psi+2)
        
        temp(1) = g4big/x
        temp(2) = g4bigp*pspx/x - g4big/x**2
        temp(3) = g4bigp*pspy/x
        temp(4) = (g4bigpp*pspx**2 + g4bigp*pspxx)/x                  &
             - 2.*g4bigp*pspx/x**2 + 2.*g4big/x**3
        temp(5) = (g4bigpp*pspx*pspy+g4bigp*pspxy)/x                  &
             - g4bigp*pspy/x**2
        temp(6) = (g4bigpp*pspy**2 + g4bigp*pspyy)/x
        call set_node_data(fun4_vec, inode, temp)

        if(.not.constraint) then
           call evaluate_spline(g2_spline, pso, g2big, g2bigp, g2bigpp)
           call evaluate_spline(g3_spline, pso, g3big, g3bigp, g3bigpp)
           ! convert from normalized to real flux
           g2big = g2big*dpsii**use_norm_psi
           g2bigp = g2bigp*dpsii**(use_norm_psi+1)
           g2bigpp = g2bigpp*dpsii**(use_norm_psi+2)
           g3big = g3big*dpsii**use_norm_psi
           g3bigp = g3bigp*dpsii**(use_norm_psi+1)
           g3bigpp = g3bigpp*dpsii**(use_norm_psi+2)
        
           temp(1) = g2big/x
           temp(2) = g2bigp*pspx/x - g2big/x**2
           temp(3) = g2bigp*pspy/x
           temp(4) = (g2bigpp*pspx**2 + g2bigp*pspxx)/x                  &
                - 2.*g2bigp*pspx/x**2 + 2.*g2big/x**3
           temp(5) = (g2bigpp*pspx*pspy+g2bigp*pspxy)/x                  &
                - g2bigp*pspy/x**2
           temp(6) = (g2bigpp*pspy**2 + g2bigp*pspyy)/x
           call set_node_data(fun2_vec, inode, temp)
           
           temp(1) = g3big/x
           temp(2) = g3bigp*pspx/x - g3big/x**2
           temp(3) = g3bigp*pspy/x
           temp(4) = (g3bigpp*pspx**2 + g3bigp*pspxx)/x                  &
                - 2.*g3bigp*pspx/x**2 + 2.*g3big/x**3
           temp(5) = (g3bigpp*pspx*pspy+g3bigp*pspxy)/x                  &
                - g3bigp*pspy/x**2
           temp(6) = (g3bigpp*pspy**2 + g3bigp*pspyy)/x
           call set_node_data(fun3_vec, inode, temp)
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

  implicit none

  include 'mpif.h'

  real, intent(out) :: error

  type(element_data) :: d
  integer :: i, itri, numelms, ier
  real :: pso, dpsii, norm, temp1(2), temp2(2)
  vectype, dimension(coeffs_per_element) :: avec
  vectype, dimension(dofs_per_element) :: temp3, temp4
      
  logical :: inside_lcfs
  real :: temp(5)

  dpsii = 1./(psilim - psimin)
  fun1_vec = 0.
  fun2_vec = 0.
  fun3_vec = 0.
  fun4_vec = 0.
  norm = 0.
  error = 0.

  if(myrank.eq.0 .and. iprint.ge.2) print *, '  poulating...'
  numelms = local_elements()
  do itri=1,numelms

     call define_element_quadrature(itri, int_pts_main, int_tor)
     call define_fields(itri, 0, 1, 0)
     call get_element_data(itri, d)

     call calcavector(itri, psi_vec, avec)
     call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
          npoints, ps079)

     do i=1, npoints
        
        pso = (ps079(i,OP_1)-psimin)*dpsii
        if(.not.inside_lcfs(ps079(i,:),x_79(i),z_79(i),.true.)) pso = 1.

        call evaluate_spline(pprime_spline,pso,temp(1))
        call evaluate_spline(ffprime_spline,pso,temp(2))
        if(irot.eq.1) then
           call evaluate_spline(alpha_spline, pso, temp(3), temp(4))
           call evaluate_spline(p0_spline,pso,temp(5))
        else
           temp(3) = 0.
           temp(4) = 0.
        end if
        temp79a(i) = temp(1)
        temp79b(i) = temp(2)
        temp79c(i) = temp(3)
        temp79d(i) = temp(4)
        temp79e(i) = temp(5)
     end do
     
     ! convert from normalized to real flux
     temp79a = temp79a*dpsii**use_norm_psi
     temp79b = temp79b*dpsii**use_norm_psi
     temp79d = temp79d*dpsii

     if(irot.eq.1) then
        temp79a = exp(temp79c*(x_79**2 - rzero**2)/rzero**2)* &
             (temp79a + temp79e*temp79d*(x_79**2 - rzero**2)/rzero**2)
     endif

     do i=1,dofs_per_element
        temp3(i) = int3(r_79, mu79(:,OP_1,i),temp79a)
        temp4(i) = int3(ri_79,mu79(:,OP_1,i),temp79b)
     end do
     call vector_insert_block(fun1_vec%vec,itri,1,temp3,VEC_ADD)
     call vector_insert_block(fun4_vec%vec,itri,1,temp4,VEC_ADD)
     
     temp79c = ps079(:,OP_GS) + r2_79*temp79a + temp79b

     norm = norm + abs(int2(ri_79,ps079(:,OP_GS)))
     error = error + abs(int2(ri_79,temp79c))
  end do

  if(myrank.eq.0 .and. iprint.ge.2) print *, '   solving...'
  call newvar_solve(fun1_vec%vec, mass_mat_lhs)
  call newvar_solve(fun4_vec%vec, mass_mat_lhs)

  if(maxrank.gt.1) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, '   allreducing...'
     temp1(1) = norm
     temp1(2) = error
     call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ier)
     norm     = temp2(1)
     error    = temp2(2)
  end if
  error = error / norm
end subroutine fundef2

subroutine readpgfiles
  use basic

  implicit none

  integer :: j, n
  real :: dum
  real, allocatable :: psinorm(:), pres0(:), g0(:), ffn(:), ppn(:)

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Reading profiles files"

  open(unit=76,file="profiles-p",status="old")
  read(76,803) n
  allocate(psinorm(n), pres0(n), ppn(n))
  do j=1,n
     read(76,802) psinorm(j), pres0(j), ppn(j), dum, dum
  enddo
  close(76)

  p0 = pres0(1)

  call create_spline(p0_spline, n, psinorm, pres0)
  call create_spline(pprime_spline, n, psinorm, ppn)
  deallocate(psinorm, pres0, ppn)

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
  
  ! in this case, ffprime and pprime are derivatives wrt actual flux
  use_norm_psi = 0

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
      n0_spline%y = den0*(p0_spline%y/p0)**expn
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

   if(myrank.eq.0) print *, "Using analytic alpha profile"

   npsi = 100
   allocate(omega(npsi), psinorm(npsi))

   do j=1, npsi
      psii = (j-1.)/(npsi-1.)
      psinorm(j) = psii

      alpha = alpha0 + alpha1*psii + alpha2*psii**2 + alpha3*psii**3 

      call evaluate_spline(p0_spline, psii, pres)
      call evaluate_spline(n0_spline, psii, dens)
      if(iscale_rot_by_p.eq.0) alpha = alpha * dens/pres
      
      omega(j) = sqrt(2./rzero**2 * alpha * pres/dens)
   end do

   call create_spline(omega_spline, npsi, psinorm, omega)

   deallocate(omega, psinorm)
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
   real :: psii, pres, dens

   if(myrank.eq.0) print *, "Calculating alpha"

   call copy_spline(alpha_spline, omega_spline)

   do j=1, alpha_spline%n
      psii = alpha_spline%x(j)

      call evaluate_spline(p0_spline, psii, pres)
      call evaluate_spline(n0_spline, psii, dens)

      alpha_spline%y(j) = 0.5*rzero**2*omega_spline%y(j)**2*dens/pres
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

   implicit none

   integer, intent(in) :: n
   real, dimension(n), intent(in) :: p, pp, f, ffp, flux

   real, allocatable :: pres0(:), g0(:), psinorm(:), ffn(:), ppn(:)

   allocate(psinorm(n), pres0(n), g0(n), ffn(n), ppn(n))

   pres0 = p
   g0    = 0.5*(f**2 - f(n)**2)
   ffn   = ffp*(flux(n) - flux(1))  ! convert to derivative wrt normalized psi
   ppn   = pp *(flux(n) - flux(1))  ! convert to derivative wrt normalized psi
   psinorm = (flux - flux(1)) / (flux(n) - flux(1))

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
   real :: y, yp, ypp, yppp

   open(unit=75,file="testout",status="unknown")
   do j=1, 1000
      yppp = (j-1.)/(1000.-1.)
      call evaluate_spline(p0_spline, yppp, y)
      call evaluate_spline(pprime_spline, yppp, yp, ypp)
      write(75,'(4e12.5)') yppp, y, yp, ypp
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
      call evaluate_spline(g0_spline, p0_spline%x(j), y, yp, ypp, yppp)
      write(77,802) j, g0_spline%x(j), y, yp, ypp, yppp
   enddo
   close(77)

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
   integer :: itri, nelms, def_fields, ier, ieqs_temp
   
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
      
      temp79a = ps079(:,OP_GS)*(ps079(:,OP_DR)**2 + ps079(:,OP_DZ)**2)
      temp79b = -r2_79* &
           (p079(:,OP_DR)*ps079(:,OP_DR) + p079(:,OP_DZ)*ps079(:,OP_DZ))
      temp79c = -bz079(:,OP_1)* &
           (bz079(:,OP_DR)*ps079(:,OP_DR) + bz079(:,OP_DZ)*ps079(:,OP_DZ))
      
      if(irot.eq.1) then
         select case(ivform)
         case(0)
            temp79d = ri_79*n079(:,OP_1)*vz079(:,OP_1)**2*ps079(:,OP_DR)
         case(1)
            temp79d = r3_79*n079(:,OP_1)*vz079(:,OP_1)**2*ps079(:,OP_DR)
         end select
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
   if(myrank.eq.0) write(*,'(A,1p2E12.4)') 'Final error, norm: ', error, norm
   if(norm.ne.0.) error = error / norm
 end subroutine calculate_gs_error

!======================================================================
! calc_toroidal_field
! ~~~~~~~~~~~~~~~~~~~
!
! calculates the toroidal field (I) as a function of the 
! normalized flux
!======================================================================
subroutine calc_toroidal_field(psi0,tf,x,z)
  use basic

  vectype, intent(in), dimension(dofs_per_node)  :: psi0
  vectype, intent(out), dimension(dofs_per_node) :: tf    ! toroidal field (I)
  real, intent(in) :: x, z
  
  vectype :: g0
  real, dimension(dofs_per_node) :: g2, g3, g4
  real :: g4big0, g4big, g4bigp
  real :: g2big, g2bigp, g3big, g3bigp
  real, dimension(dofs_per_node)  :: psii     ! normalized flux

  logical :: inside_lcfs
  
  if(.not.inside_lcfs(psi0,x,z,.true.)) then
     g0 = bzero*rzero
     call constant_field(tf,real(g0))
  else
     psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
     psii(2:6) = real(psi0(2:6))/(psilim - psimin)

     if(.not.constraint) then
        g2(1) = psii(1) - 10.*psii(1)**3 + 20.*psii(1)**4       &
             - 15.*psii(1)**5 + 4.*psii(1)**6
        g2big =  (1. - 30.*psii(1)**2 + 80.*psii(1)**3          &
             - 75.*psii(1)**4 + 24.*psii(1)**5)
        g2bigp =  (-60.*psii(1) + 240.*psii(1)**2               &
             - 300.*psii(1)**3 + 120.*psii(1)**4)
        g2(2) = (psii(2))*g2big
        g2(3) = (psii(3))*g2big
        g2(4) = (psii(4)*g2big + psii(2)**2*g2bigp)
        g2(5) = (psii(5)*g2big + psii(2)*psii(3)*g2bigp)
        g2(6) = (psii(6)*g2big + psii(3)**2*g2bigp)
     
        g3(1) = psii(1)**2 - 4.*psii(1)**3 + 6.*psii(1)**4      &
             - 4.*psii(1)**5 + psii(1)**6
        
        g3big =  (2.*psii(1) - 12.*psii(1)**2 + 24.*psii(1)**3  &
             - 20.*psii(1)**4 + 6.*psii(1)**5)
        g3bigp =  (2. - 24.*psii(1) + 72.*psii(1)**2            &
             - 80.*psii(1)**3 + 30.*psii(1)**4)
        g3(2) = (psii(2))*g3big
        g3(3) = (psii(3))*g3big
        g3(4) = (psii(4)*g3big + psii(2)**2*g3bigp)
        g3(5) = (psii(5)*g3big + psii(2)*psii(3)*g3bigp)
        g3(6) = (psii(6)*g3big + psii(3)**2*g3bigp)
     end if
     
     call evaluate_spline(g0_spline, psii(1), g4big0, g4big, g4bigp)
     
     g4(1) = g4big0
     g4(2) = (psii(2))*g4big
     g4(3) = (psii(3))*g4big
     g4(4) = (psii(4)*g4big + psii(2)**2*g4bigp)
     g4(5) = (psii(5)*g4big + psii(2)*psii(3)*g4bigp)
     g4(6) = (psii(6)*g4big + psii(3)**2*g4bigp)
     
!
!.....convert from gg' = .5(g^2)' to (g^2)'
     g2 = 2.*g2
     g3 = 2.*g3
     g4 = 2.*g4

     tf(1) = sqrt((bzero*rzero)**2 + &
          gamma2*g2(1) + gamma3*g3(1) + gamma4*g4(1))
     tf(2) = 0.5*(gamma2*g2(2) + gamma3*g3(2) + gamma4*g4(2)) / tf(1)
     tf(3) = 0.5*(gamma2*g2(3) + gamma3*g3(3) + gamma4*g4(3)) / tf(1)
     tf(4) = 0.5*(gamma2*g2(4) + gamma3*g3(4) + gamma4*g4(4)) / tf(1) &
          - (0.5*(gamma2*g2(2) + gamma3*g3(2) + gamma4*g4(2)))**2 / tf(1)**3
     tf(5) = 0.5*(gamma2*g2(5) + gamma3*g3(5) + gamma4*g4(5)) / tf(1) &
          -  0.5*(gamma2*g2(2) + gamma3*g3(2) + gamma4*g4(2)) &
          *  0.5*(gamma2*g2(3) + gamma3*g3(3) + gamma4*g4(3)) / tf(1)**3
     tf(6) = 0.5*(gamma2*g2(6) + gamma3*g3(6) + gamma4*g4(6)) / tf(1) &
          - (0.5*(gamma2*g2(3) + gamma3*g3(3) + gamma4*g4(3)))**2 / tf(1)**3
  
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
subroutine calc_pressure(psi0, pres, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: psi0
  vectype, intent(out), dimension(dofs_per_node) :: pres     ! pressure
  real, intent(in) :: x, z

  real :: fbig0, fbig, fbigp
  real :: alphap0, alphap, alphapp
  real :: ealpha, r0m, r1, r1m, r2, r3
  real :: pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(dofs_per_node) :: psii     ! normalized flux

  logical :: inside_lcfs

  psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
  psii(2:6) = real(psi0(2:6))/(psilim - psimin)

  if(.not.inside_lcfs(psi0,x,z,.true.)) psii(1) = 1.

  pspx = real(psi0(2))
  pspy = real(psi0(3))
  pspxx= real(psi0(4))
  pspxy= real(psi0(5))
  pspyy= real(psi0(6))

  call evaluate_spline(p0_spline, psii(1), fbig0, fbig, fbigp)
  fbig = fbig*dpsii
  fbigp = fbigp*dpsii**2

  if(irot.eq.1) then
     r0m = 1./rzero**2
     r1 = (x**2 - rzero**2)/rzero**2
     r1m= x**2/rzero**2
     r2 = (x**2 - rzero**2)**2/rzero**4
     r3 = (x**2 - rzero**2)**3/rzero**6

     call evaluate_spline(alpha_spline, psii(1), alphap0, alphap, alphapp)

     !...convert all derivatives to wrt psi, not normalized psi
     alphap = alphap*dpsii
     alphapp = alphapp*dpsii**2

     ealpha = exp(alphap0*r1)

     pres(1) = fbig0

     pres(2) = fbig0*alphap0*2*x*r0m  + (fbig + fbig0*alphap*r1)*pspx

     pres(3) =(fbig + fbig0*alphap*r1)*pspy

     pres(4) = fbig0*alphap0*2*r0m + fbig0*alphap0*alphap0*4*r0m*r1m &
          + (2*fbig0*alphap + 2.*fbig*alphap0           &
          +  2*fbig0*alphap0*alphap*r1)*2*x*r0m*pspx    &
          + (fbigp + (2*fbig*alphap + fbig0*alphapp)*r1 &
          +  fbig0*alphap**2*r2)*pspx*pspx              &
          + (fbig + fbig0*alphap*r1)*pspxx

     pres(5) = (fbig*alphap0 + fbig0*alphap &
          +     fbig0*alphap0*alphap*r1)*2*x*r0m*pspy       &
          +    (fbigp + (2.*fbig*alphap + fbig0*alphapp)*r1 &
          +     fbig0 *alphap**2*r2)*pspx*pspy              &
          +    (fbig + fbig0*alphap*r1)*pspxy

     pres(6) = (fbigp + (2.*fbig*alphap + fbig0*alphapp)*r1 &
          +     fbig0 *alphap**2*r2)*pspy*pspy              &
          +    (fbig + fbig0*alphap*r1)*pspyy

     pres = pres*ealpha
  else
     pres(1) = fbig0
     pres(2) = fbig*pspx
     pres(3) = fbig*pspy
     pres(4) = fbig*pspxx + fbigp*pspx**2
     pres(5) = fbig*pspxy + fbigp*pspx*pspy
     pres(6) = fbig*pspyy + fbigp*pspy**2
  endif     !....end of branch on irot            
end subroutine calc_pressure


!======================================================================
! calc_density
! ~~~~~~~~~~~~
!
! calculates the density as a function of the poloidal flux 
! (and major radius if rotation is present)
!======================================================================
subroutine calc_density(psi0, dens, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: psi0
  real, intent(in) :: x, z
  vectype, intent(out), dimension(dofs_per_node) :: dens     ! density

  real :: rbig0, rbig, rbigp
  real :: alphap0, alphap, alphapp
  real :: r0m, r1, r1m, r2, r3
  real :: ealpha, pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(dofs_per_node) :: psii     ! normalized flux

  logical :: inside_lcfs

  pspx = real(psi0(2))
  pspy = real(psi0(3))
  pspxx= real(psi0(4))
  pspxy= real(psi0(5))
  pspyy= real(psi0(6))

  psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
  psii(2:6) = real(psi0(2:6))/(psilim - psimin)

  if(.not.inside_lcfs(psi0,x,z,.true.)) psii(1) = 1.

  call evaluate_spline(n0_spline, psii(1), rbig0, rbig, rbigp)
  rbig = rbig/(psilim-psimin)
  rbigp = rbigp/(psilim-psimin)**2

  if(irot.eq.1) then     
     call evaluate_spline(alpha_spline,psii(1),alphap0,alphap,alphapp)
     alphap = alphap/(psilim-psimin)
     alphapp = alphapp/(psilim-psimin)**2

     r0m = 1./rzero**2
     r1 = (x**2-rzero**2)/rzero**2
     r1m= x**2/rzero**2
     r2 = (x**2 - rzero**2)**2/rzero**4
     r3 = (x**2 - rzero**2)**3/rzero**6

     ealpha = exp(alphap0*r1)

     dens(1) = rbig0

     dens(2) = rbig0*alphap0*2*x*r0m + (rbig + rbig0*alphap*r1)*pspx

     dens(3) = (rbig + rbig0*alphap*r1)*pspy

     dens(4) =   rbig0*alphap0*2*r0m + rbig0*alphap0*alphap0*4*r0m*r1m      &
          + ((2*rbig0*alphap + 2.*rbig*alphap0) &
          +  2*rbig0*alphap0*alphap*r1)*2*x*r0m*pspx  &
          + (rbigp + (2*rbig*alphap + rbig0*alphapp)*r1 &
          +  rbig0*alphap**2*r2)*pspx*pspx    &
          + (rbig + rbig0*alphap*r1)*pspxx

     dens(5) =  (rbig*alphap0 &
          +      rbig0*alphap +rbig0*alphap0*alphap*r1)*2*x*r0m*pspy  &
          + (rbigp + (2.*rbig*alphap + rbig0*alphapp)*r1 + &
             rbig0 *alphap**2*r2)*pspx*pspy &
          + (rbig + rbig0*alphap*r1)*pspxy

     dens(6) = (rbigp + (2.*rbig*alphap + rbig0*alphapp)*r1 &
          +     rbig0 *alphap**2*r2)*pspy*pspy  &
          + (rbig + rbig0*alphap*r1)*pspyy

     dens = dens*ealpha
  else
     dens(1) = rbig0
     dens(2) = rbig*pspx
     dens(3) = rbig*pspy
     dens(4) = rbig*pspxx + rbigp*pspx**2
     dens(5) = rbig*pspxy + rbigp*pspx*pspy
     dens(6) = rbig*pspyy + rbigp*pspy**2
  endif     !....end of branch on irot
  
  return
end subroutine calc_density


!======================================================================
! calc_electron_pressure
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the electron pressure as a function of the poloidal flux
!======================================================================
subroutine calc_electron_pressure(psi0, pe, x, z)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: psi0
  real, intent(in) :: x, z
  vectype, intent(out), dimension(dofs_per_node) :: pe     ! rotation

  vectype, dimension(dofs_per_node) :: pres0, n0
  real, dimension(dofs_per_node) :: psii          ! normalized flux
  real :: te0,tep,tepp

  logical :: inside_lcfs

  psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
  psii(2:6) = real(psi0(2:6))/(psilim - psimin)

  if(allocated(te_spline%y)) then
     if(.not.inside_lcfs(psi0,x,z,.true.)) psii(1) = 1.

     call evaluate_spline(te_spline, psii(1),te0,tep,tepp)
     call calc_density(psi0, n0, x, z)

     pe(1) = n0(1)*te0
     pe(2) = n0(2)*te0 + n0(1)*tep*psii(2)
     pe(3) = n0(3)*te0 + n0(1)*tep*psii(3)
     pe(4) = n0(4)*tep + 2.*n0(2)*tep*psii(2) &
          + n0(1)*tepp*psii(2)**2 + n0(1)*tep*psii(4)
     pe(5) = n0(5)*tep + n0(2)*tep*psii(3) + n0(3)*tep*psii(2) &
          + n0(1)*tepp*psii(2)*psii(3) + n0(1)*tep*psii(5)
     pe(6) = n0(6)*tep + 2.*n0(3)*tep*psii(3) &
          + n0(1)*tepp*psii(3)**2 + n0(1)*tep*psii(6)
  else
     call calc_pressure(psi0, pres0, x, z)
     pe = pres0*pefac
  end if
end subroutine calc_electron_pressure


!======================================================================
! calc_rotation
! ~~~~~~~~~~~~~
!
! calculates the rotation as a function of the poloidal flux
!======================================================================
subroutine calc_rotation(psi0,omega, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: psi0
  real, intent(in) :: x, z
  vectype, intent(out), dimension(dofs_per_node) :: omega     ! rotation

  real :: w0, wp, wpp
  real, dimension(dofs_per_node) :: psii     ! normalized flux

  logical :: inside_lcfs

  if(irot.eq.0) then
     omega = 0.
     return
  endif

  psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
  psii(2:6) = real(psi0(2:6))/(psilim - psimin)

  if(.not.inside_lcfs(psi0,x,z,.true.)) psii(1) = 1.

  call evaluate_spline(omega_spline, psii(1),w0,wp,wpp)

  omega(1) = w0
  omega(2) = wp*psii(2)
  omega(3) = wp*psii(3)
  omega(4) = wpp*psii(2)**2 + wp*psii(4)
  omega(5) = wpp*psii(2)*psii(3) + wp*psii(5)
  omega(6) = wpp*psii(3)**2 + wp*psii(6)

  if(ivform.eq.0) then
     ! V = omega*r^2
     omega(6) = x**2 * omega(6)
     omega(5) = x**2 * omega(5) + 2.*x*omega(3)
     omega(4) = x**2 * omega(4) + 2.*x*omega(2) + 2.*omega(1)
     omega(3) = x**2 * omega(3)
     omega(2) = x**2 * omega(2) + 2.*x*omega(1)
     omega(1) = x**2 * omega(1)
  endif
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
    
  integer :: i, izone, izonedim, index
  integer :: numnodes, ineg
  real :: normal(2), curv
  real :: x, z
  real, dimension(6) :: g
  real, dimension(1) :: xp, zp, xc, zc
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

#ifdef USE3D
  integer :: iplane, itri, nelms, nvals, j
  integer, dimension(nodes_per_element) :: inode
  integer, dimension(2) :: icol
  vectype, dimension(2) :: val
#endif

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_gs called"

  numnodes = owned_nodes()
  do i=1, numnodes

     index = node_index(rhs, i, 1)

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(is_boundary) then

        ! add feedback field
        if(idevice .eq. 0 .and. ifixedb .eq. 0) then
           xp(1) = x
           zp(1) = z
           xc(1) = 102.
           zc(1) = 10.
           call gvect(xp,zp,xc,zc,1,g,1,ineg)
           temp(1:6) = g*feedfac
           call add(psi_field(0), i, temp)
        endif
        
        if(ifixedb.ge.1) then
           temp = 0.
        else
           call get_node_data(psi_field(0), i, temp)
        endif
        
        call set_dirichlet_bc(index,rhs,temp,normal,curv,izonedim,mat)
     endif
  end do


#ifdef USE3D
  ! enforce axisymmetry

  if(int_tor.eq.0) then
     iplane = local_plane()
     nelms = local_elements()
     nvals = 2
     val(1) = -1.
     val(2) =  1.
     
     do itri=1, nelms
        call get_element_nodes(itri, inode)
        do i=1, nodes_per_element
           index = node_index(rhs, inode(i), 1)
           
           do j=1, 6
              ! if the node is not on the first plane,
              ! set its value to be the same as on the next plane
              If(iplane.ne.1 .and. i.le.pol_nodes_per_element) then
                 icol(1) = index+j-1
                 icol(2) = node_index(rhs,inode(i+pol_nodes_per_element),1)+j-1
                 if(present(mat)) &
                      call set_row_vals(mat, index+j-1, nvals, icol, val)
                 call insert(rhs, index+j-1, 0., VEC_SET)
              endif
              
              ! set toroidal derivatives to zero
              if(present(mat)) call identity_row(mat, index+j+5)
              call insert(rhs, index+j+5, 0., VEC_SET)
           end do
        end do
     end do
  endif
#endif
end subroutine boundary_gs


end module gradshafranov

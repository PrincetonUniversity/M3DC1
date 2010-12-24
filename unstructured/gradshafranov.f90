module gradshafranov

  use field

  implicit none

  integer, parameter :: numvargs = 1

  type(field_type), private :: psi_vec
  type(field_type), private :: fun1_vec, fun2_vec, fun3_vec, fun4_vec

  real, private :: dpsii
  real, private :: gamma2, gamma3, gamma4  

  logical, private :: constraint = .false.

  real, private :: gnorm, libetapeff, fac2

  integer, private :: npsi ! number of radial points in profile

  real, private, allocatable :: psinorm(:)
  real, private, allocatable :: g4big0t(:), g4bigt(:), g4bigpt(:), g4bigppt(:)
  real, private, allocatable :: fbig0t(:), fbigt(:), fbigpt(:), fbigppt(:)
  real, private, allocatable :: alphap0t(:), alphapt(:), alphappt(:), alphapppt(:)

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

     vmask = p0_l/p0
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
     call random_per(x,z,23,vmask)

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


!============================================================
subroutine gradshafranov_solve

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

  integer :: itri,i,j,ier, itnum
  integer :: numelms, numnodes
  real :: feedfac

  real :: x, phi, z, error, error2, error3 
  real :: tstart, tend

  vectype, dimension(dofs_per_node) :: tf, tm
  vectype, dimension(coeffs_per_element) :: avec
  vectype, dimension(dofs_per_element,dofs_per_element) :: temp

  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, "Calculating Grad-Shafranov Equilibrium"

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart) ! t_gs_init

  t_gs_magaxis = 0.
  t_gs_solve = 0.
  t_gs_fundef = 0.

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

  ! populate the matrix
  do itri=1,numelms

     call define_element_quadrature(itri,25,5)
     call define_fields(itri,0,1,0)

     do i=1,dofs_per_element
        do j=1,dofs_per_element
           temp(i,j) = int3(ri_79,mu79(:,OP_1,i),nu79(:,OP_GS,j))
        enddo
     enddo

     call apply_boundary_mask(itri, BOUNDARY_DIRICHLET, temp)
     call insert_block(gs_matrix, itri, 1, 1, temp, MAT_ADD)
  enddo

  feedfac = 0.

  ! insert boundary conditions
  call flush(gs_matrix)

  call boundary_gs(b2vecini_vec%vec, feedfac, gs_matrix)
  call finalize(gs_matrix)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_gs_init = tend - tstart
  endif

  ! read in numerical values for p and g functions for inumgs = 1
  if(inumgs .eq. 1) then
     call readpgfiles
  else
     if(.not.allocated(psinorm)) call default_profiles
  endif

  fbig0t = fbig0t*pscale
  fbigt = fbigt*pscale
  fbigpt = fbigpt*pscale
  fbigppt = fbigppt*pscale
  g4big0t = g4big0t*sqrt(bscale)
  g4bigt = g4bigt*sqrt(bscale)
  g4bigpt = g4bigpt*sqrt(bscale)
  g4bigppt = g4bigppt*sqrt(bscale)
  bzero = bzero*bscale

  if(myrank.eq.0) call write_profile

  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
  mainloop: do itnum=1, iabs(igs)

     if(myrank.eq.0) print *, "GS iteration = ", itnum
     
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
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        if(myrank.eq.0 .and. iprint.ge.2) print *, '  solving'

#ifdef CJ_MATRIX_DUMP
        if(itnum.le.1) then 
           call write_matrix(gs_matrix,'gs_matrix')
           call write_vector(b1vecini_vec%vec, 'gs_matrix_rhs.out')
        endif
#endif 
        call newsolve(gs_matrix,b1vecini_vec%vec,ier)

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
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_gs_solve = t_gs_solve + tend - tstart
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
     if(myrank.eq.0 .and. iprint.ge.2) print *, '  calculating funs'
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(constraint .and. igs_method.ne.1) then
        call fundef2(error3)
     else
        call fundef
     end if
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_gs_fundef = t_gs_fundef + tend - tstart
     endif

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
        call define_element_quadrature(itri, int_pts_aux, int_pts_tor)
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
              tm = tf
              call calc_density(ps079(i,:),tm,tf,x_79(i),z_79(i))
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
        call define_element_quadrature(itri, int_pts_aux, int_pts_tor)
        call define_fields(itri, 0, 1, 0)
        call get_element_data(itri, d)
        
        call calcavector(itri, psi_field(0), avec)
        call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, npoints, &
             ps079)
        call calcavector(itri, p_field(0), avec)
        call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, npoints, &
             p079)

        do i=1, npoints       
           call calc_density(ps079(i,:),p079(i,:),tf,x_79(i),z_79(i))
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
        call calc_density(psi0_l,p0_l,den0_l,x,z)
        call calc_rotation(psi0_l,vz0_l,x,z)
        call set_local_vals(i)
     end do
  end if

  pe_field(0) = p_field(0)
  call mult(pe_field(0), 1.-ipres*pi0/p0)

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
     call define_element_quadrature(itri, 25, 5)
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

  ! choose gamma2 to fix q0/qstar.  Note that there is an additional
  ! degree of freedom in gamma3.  Could be used to fix qprime(0)
  
  g2 =  -xmag**2*p0*p1 - 2.*abs(g0)/(xmag*q0*abs(dpsii))
  g3 = -4.*(abs(g0)/xmag)*djdpsi/dpsii - xmag**2*p0*p2
  g4 = -(-tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4

  if(myrank.eq.0 .and. iprint.ge.1) write(*,'(A,1p1e12.4)') ' current = ', curr
end subroutine calculate_gamma


! ===========================================================
subroutine deltafun(x,z,val,jout)

  use mesh_mod
  use basic
  use arrays
  use field
  use m3dc1_nint
  use math

  implicit none

  type(element_data) :: d
  real, intent(in) :: x, z, val
  type(field_type), intent(out) :: jout

  integer :: itri, i, k
  real :: x1, z1, si, zi, eta
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(dofs_per_element,coeffs_per_element) :: c

  itri = 0
  call whattri(x, 0., z, itri, x1, z1, IGNORE_PHI)

  if(itri.gt.0) then

     temp = 0.

     ! calculate local coordinates
     call get_element_data(itri, d)
     call global_to_local(d, x, 0., z, si, zi, eta)

     ! calculate temp_i = -val*mu_i(si,eta)
     call local_coeff_vector(itri, c, .false.)
     do i=1,dofs_per_element
        do k=1, coeffs_per_tri
           temp(i) = temp(i) - val*c(i,k)*si**mi(k)*eta**ni(k)
        end do
     end do

#ifdef USE3D
     temp = temp*twopi/nplanes
#endif

     call vector_insert_block(jout%vec, itri, jout%index, temp, VEC_ADD)
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
  real :: x, phi, z, pso, psox, psoy, psoxx, psoxy, psoyy, fbig, fbig0
  real :: fbigp, fbigpp, g4big0, g4big, g4bigp, g4bigpp, g2big, g2bigp
  real :: g2bigpp, g3big, g3bigp, g3bigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha,  pspx, pspy, pspxx, pspxy, pspyy
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
        psox = temp(2)*dpsii
        psoy = temp(3)*dpsii
        psoxx= temp(4)*dpsii
        psoxy= temp(5)*dpsii
        psoyy= temp(6)*dpsii
     
        call fget(pso, fbig0, fbig, fbigp, fbigpp)

        if(inumgs.eq.0) then
           fbig = fbig*dpsii
           fbigp = fbigp*dpsii
           fbigpp = fbigpp*dpsii
        endif
        if(irot.eq.1) then
           ! include toroidal rotation in equilibrium
           ! this section calculates pressure derivatives 
           ! in presence of rotation
           !   scj   11/19/10
           !
           ! assumes pressure of the form 
           ! p(x,psi) = p(psi) exp (alpha(psi)*r1)
           
           ! spatial derivatives of psi, not normalized psi
           pspx  = temp(2)
           pspy  = temp(3)
           pspxx = temp(4)
           pspxy = temp(5)
           pspyy = temp(6)

           r0m = 1./rzero**2
           r1 = (x**2-rzero**2)/rzero**2
           r1m= x**2/rzero**2
           r2 = (x**2 - rzero**2)**2/rzero**4
           r3 = (x**2 - rzero**2)**3/rzero**6
           call alphaget(pso,alphap0,alphap,alphapp,alphappp)

           !...convert all derivatives to wrt psi, not normalized psi
           fbigp = fbigp*dpsii
           fbigpp= fbigpp*dpsii**2
           alphap = alphap*dpsii
           alphapp = alphapp*dpsii**2
           alphappp = alphapp*dpsii**3

           ealpha = exp(alphap0*r1)

           temp(1) = ealpha*( x *fbig + fbig0*alphap*x*r1)

           temp(2) = ealpha*                                             &
                (fbig + fbig0*alphap*x*r1                                &
                +fbig0*alphap0*alphap*2*r1m*r1                           &
                + (alphap0*fbig + fbig0*alphap) * 2.*r1m                 &
                + pspx*( x*fbigp + (2.*fbig*alphap + fbig0*alphapp)*x*r1 &
                + fbig0 * alphap**2 * x*r2))

           temp(3) = ealpha*                                             &
                pspx*( x*fbigp + (2.*fbig*alphap + fbig0*alphapp)*x*r1   &
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
           temp(2) = fbig + x*fbigp*psox
           temp(3) =        x*fbigp*psoy
           temp(4) = 2.*fbigp*psox + x*(fbigpp*psox**2+fbigp*psoxx)
           temp(5) = fbigp*psoy + x*(fbigpp*psox*psoy +fbigp*psoxy)
           temp(6) = x*(fbigpp*psoy**2 + fbigp*psoyy)
           
        endif   !...end of branch on irot
        
        call set_node_data(fun1_vec, inode, temp)
        
        call g4get(pso, g4big0, g4big, g4bigp, g4bigpp)
        
        if(inumgs.eq.0) then
           g4big = g4big*dpsii
           g4bigp = g4bigp*dpsii
           g4bigpp = g4bigpp*dpsii
        endif
        
        temp(1) = g4big/x
        temp(2) = g4bigp*psox/x - g4big/x**2
        temp(3) = g4bigp*psoy/x
        temp(4) = (g4bigpp*psox**2 + g4bigp*psoxx)/x                  &
             - 2.*g4bigp*psox/x**2 + 2.*g4big/x**3
        temp(5) = (g4bigpp*psox*psoy+g4bigp*psoxy)/x - g4bigp*psoy/x**2
        temp(6) = (g4bigpp*psoy**2 + g4bigp*psoyy)/x
        call set_node_data(fun4_vec, inode, temp)
        
        g2big =  dpsii*(1. - 30.*pso**2 + 80.*pso**3                     &
             - 75.*pso**4 + 24.*pso**5)
        g2bigp =  dpsii*(-60.*pso + 240.*pso**2                         &
             - 300.*pso**3 + 120.*pso**4)
        g2bigpp =  dpsii*(-60. + 480.*pso                               &
             - 900.*pso**2 + 480.*pso**3)
        temp(1) =  g2big/x
        temp(2) =  g2bigp*psox/x - g2big/x**2
        temp(3) =  g2bigp*psoy/x
        temp(4) =  (g2bigpp*psox**2 + g2bigp*psoxx)/x                 &
             - 2.*g2bigp*psox/x**2 + 2.*g2big/x**3
        temp(5) =(g2bigpp*psox*psoy+g2bigp*psoxy)/x                   &
             - g2bigp*psoy/x**2
        temp(6) = (g2bigpp*psoy**2 + g2bigp*psoyy)/x
        call set_node_data(fun2_vec, inode, temp)
        
        g3big =  dpsii*(2.*pso - 12.*pso**2 + 24.*pso**3                &
             - 20.*pso**4 + 6.*pso**5)
        g3bigp =  dpsii*(2. - 24.*pso + 72.*pso**2                      &
             - 80.*pso**3 + 30.*pso**4)
        g3bigpp =  dpsii*(- 24. + 144.*pso                              &
             - 240.*pso**2 + 120.*pso**3)
        temp(1) = g3big/x
        temp(2) = g3bigp*psox/x - g3big/x**2
        temp(3) = g3bigp*psoy/x
        temp(4) = (g3bigpp*psox**2 + g3bigp*psoxx)/x                  &
             - 2.*g3bigp*psox/x**2 + 2.*g3big/x**3
        temp(5) = (g3bigpp*psox*psoy+g3bigp*psoxy)/x                  &
             - g3bigp*psoy/x**2
        temp(6) =  (g3bigpp*psoy**2 + g3bigp*psoyy)/x
        call set_node_data(fun3_vec, inode, temp)
     endif
  enddo

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

  numelms = local_elements()

  do itri=1,numelms

     call define_element_quadrature(itri, int_pts_aux, int_pts_tor)
     call define_fields(itri, 0, 1, 0)
     call get_element_data(itri, d)

     call calcavector(itri, psi_vec, avec)
     call eval_ops(avec, xi_79, zi_79, eta_79, d%co, d%sn, ri_79, &
          npoints, ps079)

     do i=1, npoints
        
        pso = (ps079(i,OP_1)-psimin)*dpsii
        if(inside_lcfs(ps079(i,:),x_79(i),z_79(i),.true.)) then
           call cubic_interpolation(npsi,psinorm,pso,fbigt,temp(1))
           call cubic_interpolation(npsi,psinorm,pso,g4bigt,temp(2))
           if(irot.eq.1) then
              call cubic_interpolation(npsi,psinorm,pso,alphap0t,temp(3))
              call cubic_interpolation(npsi,psinorm,pso,alphapt,temp(4))
              call cubic_interpolation(npsi,psinorm,pso,fbig0t,temp(5))
           endif
        else
           temp = 0.
        endif
        temp79a(i) = temp(1)
        temp79b(i) = temp(2)
        temp79c(i) = temp(3)
        temp79d(i) = temp(4)
        temp79e(i) = temp(5)
     end do
     
     if(inumgs.eq.0) then
        temp79a = temp79a*dpsii
        temp79b = temp79b*dpsii
        temp79d = temp79d*dpsii
     endif

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

end subroutine fundef2


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
  real :: g4big0, g4big, g4bigp, g4bigpp
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
     
     call g4get(psii(1), g4big0, g4big, g4bigp, g4bigpp)
     
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
     
!
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
! calculates the pressure as a function of the poloidal flux (and major radius if rotation is present)
!======================================================================
subroutine calc_pressure(psi0,pres, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: psi0
  real, intent(in) :: x, z
  vectype, intent(out), dimension(dofs_per_node) :: pres     ! pressure

  real :: fbig0, fbig, fbigp, fbigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha,  pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(dofs_per_node) :: psii     ! normalized flux

  logical :: inside_lcfs

  if(.not.inside_lcfs(psi0,x,z,.true.)) then
     call constant_field(pres, 0.)
  else
     psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
     psii(2:6) = real(psi0(2:6))/(psilim - psimin)
     
     call fget(psii(1), fbig0, fbig, fbigp, fbigpp)

     if(irot.eq.1) then

        pspx = real(psi0(2))
        pspy = real(psi0(3))
        pspxx= real(psi0(4))
        pspxy= real(psi0(5))
        pspyy= real(psi0(6))
        
        r0m = 1./rzero**2
        r1 = (x**2-rzero**2)/rzero**2
        r1m= x**2/rzero**2
        r2 = (x**2 - rzero**2)**2/rzero**4
        r3 = (x**2 - rzero**2)**3/rzero**6
        call alphaget(psii(1),alphap0,alphap,alphapp,alphappp)

        !...convert all derivatives to wrt psi, not normalized psi
        fbig = fbig*dpsii
        fbigp = fbigp*dpsii**2
        alphap = alphap*dpsii
        alphapp = alphapp*dpsii**2
        
        ealpha = exp(alphap0*r1)
        
        pres(1) = ealpha*fbig0
        
        pres(2) = ealpha*(fbig0*alphap0*2*x*r0m                          &
             + (fbig + fbig0*alphap*r1)*pspx)
        
        pres(3) = ealpha*(fbig + fbig0*alphap*r1)*pspy
        
        pres(4) = ealpha*(                                                  &
             fbig0*alphap0*2*r0m + fbig0*alphap0*alphap0*4*r0m*r1m          &
             +((2*fbig0*alphap + 2.*fbig*alphap0) + 2*fbig0*alphap0*alphap*r1)*2*x*r0m*pspx  &
             +(fbigp + (2*fbig*alphap + fbig0*alphapp)*r1 + fbig0*alphap**2*r2)*pspx*pspx    &
             +(fbig + fbig0*alphap*r1)*pspxx)
        
        pres(5) = ealpha*(                                               &
             (fbig*alphap0 + fbig0*alphap +fbig0*alphap0*alphap*r1)*2*x*r0m*pspy       &
             +(fbigp + (2.*fbig*alphap + fbig0*alphapp)*r1 + fbig0 *alphap**2*r2)*pspx*pspy  &
             +(fbig + fbig0*alphap*r1)*pspxy)
        
        pres(6) = ealpha*(                                               &
             +(fbigp + (2.*fbig*alphap + fbig0*alphapp)*r1 + fbig0 *alphap**2*r2)*pspy*pspy  &
             +(fbig + fbig0*alphap*r1)*pspyy)
        
     else
        
        pres(1) = fbig0
        pres(2) = psii(2)*fbig
        pres(3) = psii(3)*fbig
        pres(4) = (psii(4)*fbig + psii(2)**2*fbigp)
        pres(5) = (psii(5)*fbig + psii(2)*psii(3)*fbigp)
        pres(6) = (psii(6)*fbig + psii(3)**2*fbigp)
     endif     !....end of branch on irot
  endif

  pres(1) = pres(1) + pedge
    
end subroutine calc_pressure

subroutine readpgfiles
  use basic

  implicit none

  integer :: j, n
  real :: psii

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Reading profiles files"

  open(unit=76,file="profiles-p",status="old")
  read(76,803) n
  allocate(psinorm(n))
  allocate(fbig0t(n),fbigt(n),fbigpt(n),fbigppt(n))
  do j=1,n
    read(76,802) psinorm(j),fbig0t(j),fbigt(j),fbigpt(j),fbigppt(j)
  enddo
  close(76)

  p0 = fbig0t(1)

  open(unit=77,file="profiles-g",status="old")
  read(77,804) n
  allocate(g4big0t(n),g4bigt(n),g4bigpt(n),g4bigppt(n))
  do j=1,n
    read(77,802) psinorm(j),g4big0t(j),g4bigt(j),g4bigpt(j),g4bigppt(j)
  enddo
  close(77)

  do npsi=1, n
     if (psinorm(npsi) .ge. psiscale) exit
  end do

  psinorm = psinorm/psinorm(npsi)


  if(irot.eq.1) then
     allocate(alphap0t(npsi),alphapt(npsi),alphappt(npsi),alphapppt(npsi))
     do j=1, npsi
        psii = (j-1.)/(npsi-1.)
        alphap0t(j) = alpha0 + alpha1*psii + alpha2*psii**2
        alphapt(j)  =          alpha1   + 2.*alpha2*psii
        alphappt(j) =   2.*alpha2
        alphapppt(j) = 0.
     end do
  endif

  constraint = .true.

return
  802 format(5x,5e18.9)
  803 format(i5)
  804 format(i5)
end subroutine readpgfiles

subroutine g4get(pso, g4big0, g4big, g4bigp, g4bigpp)
  implicit none

  real, intent(in) :: pso
  real, intent(out) :: g4big0,g4big, g4bigp, g4bigpp

  call cubic_interpolation(npsi,psinorm,pso,g4big0t,g4big0)
  call cubic_interpolation(npsi,psinorm,pso,g4bigt,g4big)
  call cubic_interpolation(npsi,psinorm,pso,g4bigpt,g4bigp)
  call cubic_interpolation(npsi,psinorm,pso,g4bigppt,g4bigpp)
  return
end subroutine g4get

subroutine fget(pso, fbig0, fbig, fbigp, fbigpp)
  implicit none

  real, intent(in) :: pso
  real, intent(out) :: fbig0,fbig, fbigp, fbigpp

  call cubic_interpolation(npsi,psinorm,pso,fbig0t,fbig0)
  call cubic_interpolation(npsi,psinorm,pso,fbigt,fbig)
  call cubic_interpolation(npsi,psinorm,pso,fbigpt,fbigp)
  call cubic_interpolation(npsi,psinorm,pso,fbigppt,fbigpp)  
  return
end subroutine fget

subroutine alphaget(pso, alphap0, alphap, alphapp, alphappp)
  implicit none

  real, intent(in) :: pso
  real, intent(out) :: alphap0,alphap, alphapp, alphappp

  call cubic_interpolation(npsi,psinorm,pso,alphap0t,alphap0)
  call cubic_interpolation(npsi,psinorm,pso,alphapt,alphap)
  call cubic_interpolation(npsi,psinorm,pso,alphappt,alphapp)
  call cubic_interpolation(npsi,psinorm,pso,alphapppt,alphappp)
  return
end subroutine alphaget

 subroutine default_profiles
   
   use basic

   implicit none

   integer :: j
   real :: psii

   if(myrank.eq.0) print *, "Using default p, alpha, and I profiles"

   npsi = 500
   allocate(psinorm(npsi))
   allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
   allocate(alphap0t(npsi),alphapt(npsi),alphappt(npsi),alphapppt(npsi))
   allocate(g4big0t(npsi),g4bigt(npsi),g4bigpt(npsi),g4bigppt(npsi))
  
   do j=1, npsi
      psii = (j-1.)/(npsi-1.)
      psinorm(j) = psii
      fbig0t(j) = p0*(1.+p1*psii+p2*psii**2 &
           - (20. + 10.*p1 + 4.*p2)*psii**3 &
           + (45. + 20.*p1 + 6.*p2)*psii**4 &
           - (36. + 15.*p1 + 4.*p2)*psii**5 &
           + (10. +  4.*p1 +    p2)*psii**6)
      fbigt(j) = p0*(p1 + 2.*p2*psii - 3.*(20. + 10.*p1+4.*p2)*psii**2      &
           + 4.*(45.+20.*p1+6.*p2)*psii**3 - 5.*(36.+15.*p1+4.*p2)*psii**4  &
           + 6.*(10.+4.*p1+p2)*psii**5)
      fbigpt(j) = p0*(2.*p2 - 6.*(20. + 10.*p1+4.*p2)*psii                  &
           + 12.*(45.+20.*p1+6*p2)*psii**2 - 20.*(36.+15*p1+4*p2)*psii**3   &
           + 30.*(10.+4.*p1+p2)*psii**4)
      fbigppt(j) = p0*(- 6.*(20. + 10.*p1+4.*p2)                            &
           + 24.*(45.+20.*p1+6*p2)*psii - 60.*(36.+15*p1+4*p2)*psii**2      &
           + 120.*(10.+4.*p1+p2)*psii**3)

      alphap0t(j) = alpha0 + alpha1*psii + alpha2*psii**2
      alphapt(j)  =          alpha1   + 2.*alpha2*psii
      alphappt(j) =   2.*alpha2
      alphapppt(j) = 0.

      g4big0t(j) = 1.- 20.*psii**3+  45.*psii**4-  36.*psii**5+  10.*psii**6
      g4bigt(j)  =   - 60.*psii**2+ 180.*psii**3- 180.*psii**4+  60.*psii**5
      g4bigpt(j) =   -120.*psii   + 540.*psii**2- 720.*psii**3+ 300.*psii**4
      g4bigppt(j)=   -120.        +1080.*psii   -2160.*psii**2+1200.*psii**3
   end do

 end subroutine default_profiles

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

   real, dimension(4) :: a
   real :: dp,dpp
   integer :: j

   do npsi=1, n
      if ((flux(npsi) - flux(1)) / (flux(n) - flux(1)) .ge. psiscale) exit
   end do

   allocate(psinorm(npsi))
   allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
   allocate(alphap0t(npsi),alphapt(npsi),alphappt(npsi),alphapppt(npsi))
   allocate(g4big0t(npsi),g4bigt(npsi),g4bigpt(npsi),g4bigppt(npsi))

   fbig0t(1:n) = p                              ! p
   fbigt(1:n) = pp * (flux(n) - flux(1))        ! p' = dp/dPsi
   g4big0t(1:n) = 0.5*(f**2 - f(n)**2)          ! g
   g4bigt(1:n) = ffp * (flux(n) - flux(1))      ! f f' = f * df/dPsi

   ! calculate normalized flux
   ! and derivatives wrt acutal flux
   do j=1,npsi
      call cubic_interpolation_coeffs(flux,npsi,j,a)
      psinorm(j) = (flux(j) - flux(1)) / (flux(npsi) - flux(1))
      dp = a(2)      ! d psi/dn
      dpp = 2.*a(3)  ! d^2 psi/dn^2

      call cubic_interpolation_coeffs(fbigt,npsi,j,a)
      fbigpt(j) =     a(2)/dp                       ! p''
      fbigppt(j) = 2.*a(3)/dp**2 - a(2)*dpp/dp**3   ! p'''

      call cubic_interpolation_coeffs(g4bigt,npsi,j,a)
      g4bigpt(j)  =    a(2)/dp                      ! (f f')'
      g4bigppt(j) = 2.*a(3)/dp**2 - a(2)*dpp/dp**3  ! (f f')''

      if(irot.eq.1) then
         alphap0t(j) =  alpha0 + alpha1*psinorm(j) + alpha2*psinorm(j)**2
         alphapt(j)  =           alpha1         + 2.*alpha2*psinorm(j)
         alphappt(j) = 2.*alpha2
         alphapppt(j) = 0.
      endif
   end do

   ! change derivatives from actual to normalized flux
   fbigpt = fbigpt*(flux(npsi) - flux(1))
   fbigppt = fbigppt*(flux(npsi) - flux(1))**2
   g4bigpt = g4bigpt*(flux(npsi) - flux(1))
   g4bigppt = g4bigppt*(flux(npsi) - flux(1))**2
   alphapt = alphapt*(flux(npsi) - flux(1))
   alphappt = alphappt*(flux(npsi) - flux(1))**2
  
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

   open(unit=76,file="profilesdb-p",status="unknown")
   do j=1,npsi
      write(76,802) j,psinorm(j),fbig0t(j),fbigt(j),fbigpt(j),fbigppt(j)
   enddo
   close(76)
   
   open(unit=77,file="profilesdb-g",status="unknown")
   do j=1,npsi
      write(77,802) j,psinorm(j),g4big0t(j),g4bigt(j),g4bigpt(j),g4bigppt(j)
   enddo
802 format(i5,1p6e18.9)
   close(77)
   
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
      call define_element_quadrature(itri, int_pts_main, int_pts_tor)
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
! calc_density
! ~~~~~~~~~~~~~
!
! calculates the density as a function of the poloidal flux (and major radius if rotation is present)
!======================================================================
subroutine calc_density(psi0,pres,dens, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: psi0, pres
  real, intent(in) :: x, z
  vectype, intent(out), dimension(dofs_per_node) :: dens     ! density

  real :: fbig0, fbig, fbigp, fbigpp
  real :: rbig0, rbig, rbigp, p0ni
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha,  pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(dofs_per_node) :: psii     ! normalized flux

  logical :: inside_lcfs

  if(irot.eq.1) then
     if(.not.inside_lcfs(psi0,x,z,.true.)) then
        fbig0 = 0.
        fbig = 0.
        fbigp = 0.
        fbigpp = 0.
     else
        psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
        psii(2:6) = real(psi0(2:6))/(psilim - psimin)
        
        call fget(psii(1), fbig0, fbig, fbigp, fbigpp)

        ! for inumgs.ne.0 fbig is derivative wrt total psi, not normalized
!!$        if(inumgs.ne.0) then
!!$           fbig = fbig*(psilim-psimin)
!!$           fbigp = fbigp*(psilim-psimin)
!!$        endif
     endif
     fbig0 = fbig0 + pedge

!.....include toroidal rotation in equilibrium
!...this section calculates density derivatives in presence of rotation
!   scj   11/19/10
!
!....assumes density of the form p(x,psi) = p(psi) exp (alpha(psi)*r1)
!
!...spatial derivatives of psi, not normalized psi
     pspx = real(psi0(2))
     pspy = real(psi0(3))
     pspxx= real(psi0(4))
     pspxy= real(psi0(5))
     pspyy= real(psi0(6))

     r0m = 1./rzero**2
     r1 = (x**2-rzero**2)/rzero**2
     r1m= x**2/rzero**2
     r2 = (x**2 - rzero**2)**2/rzero**4
     r3 = (x**2 - rzero**2)**3/rzero**6
     call alphaget(psii(1),alphap0,alphap,alphapp,alphappp)

!...convert all derivatives to wrt psi, not normalized psi
     fbig = fbig/(psilim-psimin)
     fbigp = fbigp/(psilim-psimin)**2
     alphap = alphap/(psilim-psimin)
     alphapp = alphapp/(psilim-psimin)**2

     p0ni = den0*(1./p0)**expn
     rbig0 = p0ni*fbig0**expn
     rbig = 0.
     rbigp = 0.
     if(fbig0.gt.0.) then
        rbig = p0ni*expn*fbig0**(expn-1)*fbig
        rbigp= p0ni*expn*((expn-1)*fbig0**(expn-2)*fbig**2 + fbig0**(expn-1.)*fbigp)
     endif

     ealpha = exp(alphap0*r1)

     dens(1) = ealpha*rbig0

     dens(2) = ealpha*(rbig0*alphap0*2*x*r0m                                                &
                     + (rbig + rbig0*alphap*r1)*pspx)

     dens(3) = ealpha*(rbig + rbig0*alphap*r1)*pspy

     dens(4) = ealpha*(                                                                     &
             rbig0*alphap0*2*r0m + rbig0*alphap0*alphap0*4*r0m*r1m                          &
            +((2*rbig0*alphap + 2.*rbig*alphap0) + 2*rbig0*alphap0*alphap*r1)*2*x*r0m*pspx  &
            +(rbigp + (2*rbig*alphap + rbig0*alphapp)*r1 + rbig0*alphap**2*r2)*pspx*pspx    &
            +(rbig + rbig0*alphap*r1)*pspxx)

     dens(5) = ealpha*(                                                                     &
             (rbig*alphap0 + rbig0*alphap +rbig0*alphap0*alphap*r1)*2*x*r0m*pspy       &
            +(rbigp + (2.*rbig*alphap + rbig0*alphapp)*r1 + rbig0 *alphap**2*r2)*pspx*pspy  &
            +(rbig + rbig0*alphap*r1)*pspxy)

     dens(6) = ealpha*(                                                                     &
            +(rbigp + (2.*rbig*alphap + rbig0*alphapp)*r1 + rbig0 *alphap**2*r2)*pspy*pspy  &
            +(rbig + rbig0*alphap*r1)*pspyy)

  else
     if(expn .gt. 0.) then
        dens(1) = pres(1)**expn
        dens(2) = pres(1)**(expn-1.)*pres(2)*expn
        dens(3) = pres(1)**(expn-1.)*pres(3)*expn
        dens(4) = pres(1)**(expn-1.)*pres(4)*expn &
             + pres(1)**(expn-2.)*pres(2)**2.*expn*(expn-1.)
        dens(5) = pres(1)**(expn-1.)*pres(5)*expn &
             + pres(1)**(expn-2.)*pres(2)*pres(3) &
             * expn*(expn-1.)
        dens(6) = pres(1)**(expn-1.)*pres(6)*expn &
             + pres(1)**(expn-2.)*pres(3)**2.*expn*(expn-1.)
        dens = den0*dens/p0**expn
     else   ! expn.eq.0
        dens(1) = 1.
        dens(2:6) = 0.
     endif
  endif     !....end of branch on irot
  
  return
end subroutine calc_density


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

  real :: fbig0, fbig, fbigp, fbigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: tp0,tp,tpp,p0n,befoh,befomh,befom3h
  real :: r0m, r1, r1m, r2, r3,  pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(dofs_per_node) :: psii     ! normalized flux

  logical :: inside_lcfs

  if(irot.eq.0) then
     omega = 0.
     return
  endif

  if(.not.inside_lcfs(psi0,x,z,.true.)) then
     fbig0 = 0.
     fbig = 0.
     fbigp = 0.
     fbigpp = 0.
  else
     psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
     psii(2:6) = real(psi0(2:6))/(psilim - psimin)
     
     call fget(psii(1), fbig0, fbig, fbigp, fbigpp)
     
     ! for inumgs.ne.0 fbig is derivative wrt total psi, not normalized
!!$     if(inumgs.ne.0) then
!!$        fbig = fbig*(psilim-psimin)
!!$        fbigp = fbigp*(psilim-psimin)
!!$     endif
  endif
  fbig0 = fbig0 + pedge

!.....include toroidal rotation in equilibrium
!...this section calculates rotation derivatives in presence of rotation
!   scj   11/19/10
!
!....assumes rotation of the form p(x,psi) = p(psi) exp (alpha(psi)*r1)
!
!...spatial derivatives of psi, not normalized psi
  pspx = real(psi0(2))
  pspy = real(psi0(3))
  pspxx= real(psi0(4))
  pspxy= real(psi0(5))
  pspyy= real(psi0(6))
  
  r0m = 1./rzero**2
  r1 = (x**2-rzero**2)/rzero**2
  r1m= x**2/rzero**2
  r2 = (x**2 - rzero**2)**2/rzero**4
  r3 = (x**2 - rzero**2)**3/rzero**6
  call alphaget(psii(1),alphap0,alphap,alphapp,alphappp)

!...convert all derivatives to wrt psi, not normalized psi
  fbig = fbig/(psilim-psimin)
  fbigp = fbigp/(psilim-psimin)**2
  alphap = alphap/(psilim-psimin)
  alphapp = alphapp/(psilim-psimin)**2

!...define temperature and derivatives
  p0n = p0**expn/den0
  tp0 = p0n*fbig0**(1.-expn)
  tp  = p0n*(1.-expn)*fbig0**(-expn)*fbig
  tpp = p0n*(expn*(expn-1.)*fbig0**(-1.-expn)*fbig**2 + (1.-expn)*fbig0**(-expn)*fbigp)
  
  befoh   = (2.*alphap0*tp0/rzero**2)**0.5
  befomh  = (2.*alphap0*tp0/rzero**2)**(-0.5)
  befom3h = (2.*alphap0*tp0/rzero**2)**(-1.5)
  
  if(befomh.gt.0.) then      
     omega(1) = befoh
  
     omega(2) = befomh*r0m*(alphap*tp0 + alphap0*tp)*pspx
  
     omega(3) = befomh*r0m*(alphap*tp0 + alphap0*tp)*pspy
  
     omega(4) = -befom3h*r0m**2*(alphap*tp0 + alphap0*tp)**2*pspx*pspx                   &
              +   befomh*r0m*((alphapp*tp0 + 2.*alphap*tp + alphap0*tpp)*pspx*pspx       &
                             +(alphap*tp0 + alphap0*tp)*pspxx)
  
     omega(5) =  -befom3h*r0m**2*(alphap*tp0 + alphap0*tp)**2*pspx*pspy                  &
              +   befomh*r0m*((alphapp*tp0 + 2.*alphap*tp + alphap0*tpp)*pspx*pspy       &
                                   +(alphap*tp0 + alphap0*tp)*pspxy)

     omega(6) =  -befom3h*r0m**2*(alphap*tp0 + alphap0*tp)**2*pspy*pspy                  &
              +   befomh*r0m*((alphapp*tp0 + 2.*alphap*tp + alphap0*tpp)*pspy*pspy       &
                             +(alphap*tp0 + alphap0*tp)*pspyy)
  else
     omega = 0.
  endif

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

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_gs called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     index = node_index(rhs, i, 1)

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
  end do
end subroutine boundary_gs


end module gradshafranov

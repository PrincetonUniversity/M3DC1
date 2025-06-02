module init_common

  implicit none

contains

subroutine init_planewave(outarr,x,phi,z,kx,kphi,kz,amp,phasex,phasephi,phasez)
  use element

  implicit none

  vectype, intent(out), dimension(MAX_PTS) :: outarr
  real, intent(in), dimension(MAX_PTS)  :: x, phi, z
  real, intent(in) :: kx, kphi, kz, amp, phasex, phasephi, phasez

  real, dimension(MAX_PTS) :: argx,argp,argz
  argx = kx*x     + phasex
  argp = kphi*phi + phasephi
  argz = kz*z     + phasez

  outarr = amp*cos(argp)*sin(argx)*sin(argz)

#ifdef USECOMPLEX
  outarr = outarr + amp*(0,1)*sin(argp)*sin(argx)*sin(argz)
#endif
end subroutine init_planewave


subroutine init_random(x,phi,z,outarr)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in), dimension(MAX_PTS) :: x, phi, z
  vectype, intent(out), dimension(MAX_PTS) :: outarr
  integer, allocatable :: seed(:)
  integer :: i, j, n
  real, dimension(MAX_PTS) :: xx, zz, rsq, r, theta
  vectype, dimension(MAX_PTS) :: temp, phase
  real :: alx, alz, kx, kp, kz, px, pp, pz, random, roundoff

  outarr = 0.
  if(eps.eq.0.) return

  call get_bounding_box_size(alx, alz)

!
! changed to be consistent with fortran95 7/12/2011 ... scj
  call random_seed(SIZE = n)
  allocate(seed(n))
  seed = 23
  call random_seed(PUT = seed)
  deallocate(seed)

  roundoff = 1.e-12

  xx = x - xzero
  zz = z - zzero
  theta = atan2(zz,xx)
#ifdef USECOMPLEX
  if(itor.eq.1) then
    phase = exp((0,1)*(ntor*phi - mpol*theta))
  else  
    phase = exp((0,1)*(ntor*phi/rzero - mpol*theta))
  end if  
#else
  if(itor.eq.1) then
    phase = cos(ntor*phi - mpol*theta)
  else  
    phase = cos(ntor*phi/rzero - mpol*theta)
  end if
#endif

  select case (icsym)

  case (0)   !  original option...no symmetry imposed
     do i=1,maxn
        kx = pi*i/alx
        do j=1, maxn
           kz = j*pi/alz
           kp = j
           call random_number(px)
           call random_number(pp)
           call random_number(pz)
           px = 2.*pi*px
           pp = 2.*pi*pp
           pz = 2.*pi*pz
           call random_number(random)
           call init_planewave(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
                px,pp,pz)
           outarr = outarr + temp
        end do
     end do

  case (1)  !   odd symmetry about midplane
     do i=1,maxn
        kx = pi*i/alx
        do j=1, maxn/2
           kz = 2.*j*pi/alz
           kp = j
           call random_number(random)
           call init_planewave(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
                0.,0.,0.)
           outarr = outarr + temp
        end do
     end do

  case (2)  !   even symmetry about midplane
     do i=1,maxn
        kx = pi*i/alx
        do j=1, maxn/2
           kz = (2.*j-1)*pi/alz
           kp = j
           call random_number(random)
           call init_planewave(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
                0.,0.,0.)
           outarr = outarr + temp
        end do
     end do

  case (3)  !   NOT RANDOM....start in (1,1) eigenfunction
     rsq = xx**2 + zz**2 + roundoff
     r = sqrt(rsq)

     outarr = eps*r*phase
     if(ln.gt.0.) outarr = outarr*exp(-rsq/ln)
             
  end select
end subroutine init_random

subroutine init_perturbations
  use basic
  use arrays
  use field
  use m3dc1_nint
  use newvar_mod
  use diagnostics

  implicit none

  type(field_type) :: psi_vec, phi_vec
  integer :: itri, numelms, i, izone, ifield
  integer, dimension(MAX_PTS) :: imr
  vectype, dimension(dofs_per_element) :: dofs

  call create_field(psi_vec)
  call create_field(phi_vec)

  psi_vec = 0.
  phi_vec = 0.

  ifield = FIELD_PSI + FIELD_P

  numelms = local_elements()

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining initial perturbations'

  do itri=1,numelms
     call get_zone(itri, izone)

     if(izone.ne.1) cycle

     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,ifield,1,0)

     ps179 = 0.
     ph179 = 0.

     ! calculate perturbed fields
#ifdef USEST
     if(igeometry.eq.1) then
        call init_random(xl_79, phi_79, zl_79, ph179(:,OP_1))
     else
        call init_random(x_79-xmag, phi_79, z_79, ph179(:,OP_1))
     endif
#else
     call init_random(x_79-xmag, phi_79, z_79, ph179(:,OP_1))
#endif

     ph179(:,OP_1) = ph179(:,OP_1) + r_79*verzero

     ! apply mask
     if(p0 .gt. pedge) then
        call magnetic_region(pst79(:,OP_1),pst79(:,OP_DR),pst79(:,OP_DZ), &
             x_79(:), z_79(:), imr)

#ifdef USEPARTICLES
        where((imr.eq.REGION_PLASMA).and.(abs(pt79(:,OP_1))>10*pedge))
#else
        where(imr.eq.REGION_PLASMA)
#endif
           temp79a(:) = (pt79(:,OP_1) - pedge)/p0
        elsewhere
           temp79a(:) = 0.
        end where
        ph179(:,OP_1) = ph179(:,OP_1)*temp79a
     end if

     ! populate vectors for solves
        
     ! psi
     dofs = intx2(mu79(:,:,OP_1),ps179(:,OP_1))
     call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)
     
     ! phi
     dofs = intx2(mu79(:,:,OP_1),ph179(:,OP_1))
     call vector_insert_block(phi_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  ! do solves
  call newvar_solve(psi_vec%vec,mass_mat_lhs)
  psi_field(1) = psi_vec

  ! use dirichlet boundary conditions to avoid perturbations on boundary
!  call newvar_solve(phi_vec%vec,mass_mat_lhs_dc)
  call newvar_solve(phi_vec%vec,mass_mat_lhs)
  u_field(1) = phi_vec

  call destroy_field(psi_vec)
  call destroy_field(phi_vec)

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, 'Done defining initial perturbations'
end subroutine init_perturbations

subroutine den_eq
  use basic
  use arrays
  use diagnostics
  use math
  use mesh_mod
  use m3dc1_nint
  use newvar_mod
  use pellet
  use read_ascii
  use spline

  implicit none

  type(field_type) :: den_vec
  integer :: itri, numelms, def_fields
#ifdef USEPARTICLES
  integer :: nvals, j
  real, allocatable :: xvals(:), yvals(:)
  real :: val, val2
  type(spline1d) :: nf_spline      ! Ion density
  type(spline1d) :: tf_spline      ! Ion density
#endif
  real :: rate
  vectype, dimension(dofs_per_element) :: dofs
  real, dimension(MAX_PTS) :: n, p
  integer :: ip, izone
  
  if((idenfunc.eq.0 .or. idenfunc.eq.4) .and. .not.(ipellet.gt.0 .and. linear.eq.1)) return
  if(ipellet.ne.0) then
     if(ipellet_z.ne.0 .and. all(pellet_mix.eq.0.)) return
  end if

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Defining density equilibrium'
  call create_field(den_vec)
  
  def_fields = FIELD_PSI + FIELD_N + FIELD_P + FIELD_TE

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,def_fields,1,0,1)
     call get_zone(itri, izone)

     if(iread_ne.eq.0) then 
        select case(idenfunc)
        case(1)
           n079(:,OP_1) = den0*.5* &
                (1. + &
                tanh((real(ps079(:,OP_1))-(psibound+denoff*(psibound-psimin)))&
                /(dendelt*(psibound-psimin))))
           
        case(2)        
           n079(:,OP_1) = den0
           if(den0.ne.den_edge) then
              temp79a = ((real(ps079(:,OP_1))-psimin)/(psibound-psimin) - denoff)&
                   /dendelt
              n079(:,OP_1) = n079(:,OP_1) &
                   + .5*(den_edge-den0)*(1. + tanh(real(temp79a)))
           end if

        ! idenfunc = 20 in lieu of 0 (unclear why 0 is disabled here). YZ
        case(20) 
           if(expn.eq.0.) then
              n079(:,OP_1) = den0 + den_edge
           else
#ifndef USECOMPLEX
              do ip = 1, MAX_PTS
                 temp79a(ip) = max(p079(ip,OP_1), 0.)
              end do
#endif     
              n079(:,OP_1) = den0*(temp79a/p0)**expn + den_edge
           end if
#ifdef USEST
        case(21)
           if(igeometry.eq.1) then
              temp79b = sqrt((xl_79-xcenter)**2 + (zl_79-zcenter)**2 + regular**2)
              n079(:,OP_1) = den0 + (den_edge-den0)*.5*&
                   (1. + tanh((temp79b-(1.+denoff))/dendelt))
           else
              if(myrank.eq.0) print *, 'idenfunc = 21 requires igeometry = 1'
           end if
#endif     
        end select
     end if

     if(ipellet.gt.0 .and. linear.eq.1) then
        n = 0.
        p = 0.
        do ip=1,npellets
           if(pellet_mix(ip).eq.0.) then
              rate = pellet_rate(ip)
           else
              rate = pellet_rate_D2(ip)*2.0 ! two deuterium ions per D2 molecule
           end if
           n079(:,OP_1) = n079(:,OP_1) + rate*pellet_distribution(ip, x_79, phi_79, z_79, p, 1, izone)
        end do
     end if

     dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
     call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(den_vec%vec,mass_mat_lhs)
  den_field(0) = den_vec

#if defined(USEST) && defined(USEPARTICLES)
  if ((kinetic.eq.1).and.(kinetic_fast_ion.eq.1)) then
     nvals = 0
     call read_ascii_column('nf_profile', xvals, nvals, icol=1)
     call read_ascii_column('nf_profile', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     yvals = yvals / 1e6 / n0_norm !rsae
     xvals = xvals / xvals(nvals) ! normalize rho
     if(allocated(yvals)) then
        call create_spline(nf_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
     den_vec=0.
     do itri=1,numelms
        call define_element_quadrature(itri,int_pts_main,int_pts_tor)
        call define_fields(itri,def_fields,1,0)
        call get_zone(itri, izone)
        temp79b = (xl_79-xcenter)**2 + (zl_79-zcenter)**2
        do j=1, npoints
              call evaluate_spline(nf_spline,temp79b(j),val)
              n079(j,OP_1) = val
        end do
        dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
        call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
     end do
     call newvar_solve(den_vec%vec,mass_mat_lhs)
     nf_field = den_vec

     nvals = 0
     call read_ascii_column('tf_profile', xvals, nvals, icol=1)
     call read_ascii_column('tf_profile', yvals, nvals, icol=2)
     if(nvals.eq.0) call safestop(5)
     !yvals = yvals / n0_norm !rsae
     xvals = xvals / xvals(nvals) ! normalize rho
     if(allocated(yvals)) then
        call create_spline(tf_spline, nvals, xvals, yvals)
        deallocate(xvals, yvals)
     end if
     den_vec=0.
     do itri=1,numelms
        call define_element_quadrature(itri,int_pts_main,int_pts_tor)
        call define_fields(itri,def_fields,1,0)
        call get_zone(itri, izone)
        temp79b = (xl_79-xcenter)**2 + (zl_79-zcenter)**2
        do j=1, npoints
              call evaluate_spline(tf_spline,temp79b(j),val)
              n079(j,OP_1) = val
        end do
        dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
        call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
     end do
     call newvar_solve(den_vec%vec,mass_mat_lhs)
     tf_field = den_vec

     den_vec=0.
     do itri=1,numelms
        call define_element_quadrature(itri,int_pts_main,int_pts_tor)
        call define_fields(itri,def_fields,1,0)
        call get_zone(itri, izone)
        temp79b = (xl_79-xcenter)**2 + (zl_79-zcenter)**2
        do j=1, npoints
          call evaluate_spline(nf_spline,temp79b(j),val)
          call evaluate_spline(tf_spline,temp79b(j),val2)
          if (fast_ion_dist==1) then
            n079(j,OP_1) = val*val2* 1.6022e-12 / (b0_norm**2/(4.*pi*n0_norm))!rsae
          endif
        end do
        dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
        call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
     end do
     call newvar_solve(den_vec%vec,mass_mat_lhs)
     pf_field = den_vec
  endif

  if (kinetic.eq.1) then
     den_vec=0.
     do itri=1,numelms
        call define_element_quadrature(itri,int_pts_main,int_pts_tor)
        call define_fields(itri,def_fields,1,0)
        call get_zone(itri, izone)
        if(igeometry.eq.1) then
           temp79b = sqrt((xl_79-xcenter)**2 + (zl_79-zcenter)**2+1e-8)
           n079(:,OP_1) = temp79b
        else
           if(myrank.eq.0) print *, 'idenfunc = 21 requires igeometry = 1'
        end if
        dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
        call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
     end do

     call newvar_solve(den_vec%vec,mass_mat_lhs)
     rho_field = den_vec
  endif
#endif     

  !call mult(ti_field(0), 0.5)
  !call mult(te_field(0), 0.5)
  !call mult(p_field(0), 0.5)
  !call mult(pe_field(0), 0.5)
  !pfi_field = p_field(0)

  call destroy_field(den_vec)

end subroutine den_eq

subroutine den_per
  use basic
  use arrays
  use pellet
  use diagnostics
  use field
  use m3dc1_nint
  use newvar_mod
  implicit none

  type(field_type) :: den_vec
  integer :: numelms, itri
  real :: rate
  vectype, dimension(dofs_per_element) :: dofs
  real, dimension(MAX_PTS) :: n, p
  integer :: ip, izone

  if(ipellet.ge.0) return
  if(ipellet_z.ne.0 .and. all(pellet_mix.eq.0.)) return

  call create_field(den_vec)
  den_vec = 0.

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,0,1,0)
     call get_zone(itri, izone)

     n179(:,OP_1) = 0.
     ! negative ipellet in initial perturbation
     if(ipellet.lt.0) then
        n = 0.
        p = 0.
        do ip=1,npellets
           if(pellet_mix(ip).eq.0) then
              rate = pellet_rate(ip)
           else
              rate = pellet_rate_D2(ip)*2.0 ! two deuterium ions per D2 molecule
           end if
           n179(:,OP_1) = n179(:,OP_1) + rate*pellet_distribution(ip, x_79, phi_79, z_79, p, 1, izone)
        end do
     end if

     dofs = intx2(mu79(:,:,OP_1),n179(:,OP_1))

     call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(den_vec%vec,mass_mat_lhs)
  den_field(1) = den_vec

  call destroy_field(den_vec)

end subroutine den_per

subroutine nre_eq
  use basic
  use arrays
  use diagnostics
  use math
  use mesh_mod
  use m3dc1_nint
  use newvar_mod
  use pellet

  implicit none

  type(field_type) :: nre_vec
  integer :: itri, numelms, def_fields
  real :: rate, nr_a, nr_b, nr_w
  vectype, dimension(dofs_per_element) :: dofs
  real, dimension(MAX_PTS) :: n, rr, mask, nr_f1, nr_f2

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Defining RE equilibrium'
  call create_field(nre_vec)

  def_fields = FIELD_PSI + FIELD_RE

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,def_fields,1,0,1)

     nre079(:,OP_1) = 1.0*ri_79*ps079(:,OP_GS)/1.000
     if(irunaway == 1) nre079(:,OP_1) = 0.8e-0*nre079(:,OP_1)
     if(irunaway == 2) nre079(:,OP_1) = 0.

     dofs = intx2(mu79(:,:,OP_1),nre079(:,OP_1))
     call vector_insert_block(nre_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(nre_vec%vec,mass_mat_lhs)
  nre_field(0) = nre_vec

  call destroy_field(nre_vec)

end subroutine nre_eq

subroutine nre_per
  use basic
  use arrays
  use pellet
  use diagnostics
  use field
  use m3dc1_nint
  use newvar_mod
  use math
  implicit none

  type(field_type) :: nre_vec
  integer :: numelms, itri, def_fields
  real :: rate, nr_a, nr_b, nr_w
  vectype, dimension(dofs_per_element) :: dofs
  real, dimension(MAX_PTS) :: n, p
  real,dimension(MAX_PTS) :: rr,mask,nr_f1,nr_f2
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Defining RE perturbation'
  call create_field(nre_vec)
  nre_vec = 0.

  def_fields = FIELD_PSI + FIELD_RE

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,def_fields,1,0)
     nre179(:,OP_1) = 0
     dofs = intx2(mu79(:,:,OP_1),nre179(:,OP_1))

     call vector_insert_block(nre_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(nre_vec%vec,mass_mat_lhs)
  nre_field(1) = nre_vec

  call destroy_field(nre_vec)

end subroutine nre_per

end module init_common

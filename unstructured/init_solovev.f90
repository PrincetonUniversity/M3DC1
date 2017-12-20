! This module defines the Solov'ev equilibrium specified in
! Berkery et al. Phys. Plasmas 21, 052505 (2014)

module solovev
  implicit none

  vectype, private :: psi_offset

contains

  subroutine solovev_init()
    use basic
    use arrays
    use m3dc1_nint
    use mesh_mod
    use newvar_mod
    use field
    use fit_magnetics

    implicit none

    type(field_type) :: psi_vec, bz_vec, den_vec, p_vec
    integer :: itri, numelms, i, k
    real, dimension(1) :: r(1), z(1)
    vectype, dimension(1) :: psi
    vectype, dimension(dofs_per_element) :: dofs

    if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining Solov'ev Equilibrium"

    call solovev_vac

    xlim = sqrt(2.*(ln*rzero) + rzero**2)
    zlim = 0.
    r = xlim
    z = zlim
!    call evaluate_bessel_fields_2D(1,r,z,psi)
    call evaluate_multipole_fields_2D(1,r,z,psi)
    psi_offset = psi(1)
    call psi_vac(1,r,z,psi)
    psi_offset = psi_offset + psi(1)
    if(myrank.eq.0) then
       print *, "Limiter = ", xlim
       print *, "psi_lim = ", psi_offset
    end if

    call create_field(psi_vec)
    call create_field(bz_vec)
    call create_field(p_vec)
    call create_field(den_vec)

    numelms = local_elements()

    do k=0,1
       psi_vec = 0.
       bz_vec = 0.
       p_vec = 0.
       den_vec = 0.
       
       do itri=1,numelms
          call define_element_quadrature(itri,int_pts_main,int_pts_tor)
          call define_fields(itri,0,1,0)

          if(k.eq.0) then 
             ! calculate equilibrium fields
             call solovev_equ
          else
             ! calculate perturbed fields
             call solovev_per
          end if

          ! populate vectors for solves

          ! psi
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),ps079(:,OP_1))
          end do
          call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! bz
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),bz079(:,OP_1))
          end do
          call vector_insert_block(bz_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! p
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),p079(:,OP_1))
          end do
          call vector_insert_block(p_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! den
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),n079(:,OP_1))
          end do
          call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
       end do

       ! do solves
       call newvar_solve(psi_vec%vec,mass_mat_lhs)
       psi_field(k) = psi_vec
       
       call newvar_solve(bz_vec%vec,mass_mat_lhs)
       bz_field(k) = bz_vec
       
       call newvar_solve(p_vec%vec,mass_mat_lhs)
       p_field(k) = p_vec
       pe_field(k) = p_vec
       call mult(pe_field(k), pefac)
       
       call newvar_solve(den_vec%vec,mass_mat_lhs)
       den_field(k) = den_vec
    end do

    call destroy_field(psi_vec)
    call destroy_field(bz_vec)
    call destroy_field(den_vec)
    call destroy_field(p_vec)

!    call clear_bessel_fields

  end subroutine solovev_init

  subroutine solovev_equ()
    use basic
    use m3dc1_nint
    use fit_magnetics

    implicit none

    real, dimension(MAX_PTS) :: g
    real, dimension(MAX_PTS) :: psin

    bz079(:,OP_1) = 1.

    g = (x_79*z_79/elongation)**2 + &
         0.25*(x_79**2 - rzero**2)**2 - (ln*rzero)**2

!    ps079(:,OP_1) = elongation/(2.*rzero**3*q0) * g
!    call psi_vac(npoints,x_79(1:npoints),z_79(1:npoints),ps079(1:npoints,OP_1))
!    ps079(:,OP_1) = ps079(:,OP_1) - psi_offset

    ! total field is field from plasma current plus vacuum field
!    call evaluate_bessel_fields_2D(npoints, x_79, z_79, temp79a)
    call evaluate_multipole_fields_2D(npoints, x_79, z_79, temp79a)
    call psi_vac(npoints,x_79,z_79,temp79b)
    ps079(:,OP_1) = temp79a + temp79b - psi_offset

!    where(g .lt. 0)
    where(real(ps079(:,OP_1)).lt.0.)
!       ps079(:,OP_1) = elongation/(2.*rzero**3*q0) * g
       p079(:,OP_1) = pedge - g * (1. + elongation**2) &
            / (2.*(rzero**3*q0)**2)
       p079(:,OP_1) = pedge - ps079(:,OP_1) * (1. + elongation**2) &
            / (elongation*rzero**3*q0)
       psin = 1. + g/(ln*rzero)**2
       n079(:,OP_1) = den0*(1. - 0.7*psin)
    elsewhere
       p079(:,OP_1) = pedge
       n079(:,OP_1) = 0.3*den0
    end where

  end subroutine solovev_equ

  subroutine solovev_per()
    use basic
    use m3dc1_nint
    use init_common

    implicit none

    bz079(:,OP_1) = 0.
    n079(:,OP_1) = 0.
    p079(:,OP_1) = 0.
    call init_random(x_79-rzero, phi_79, z_79, ps079(:,OP_1))
    
  end subroutine solovev_per

  subroutine solovev_vac()
    use basic
    use math
    use fit_magnetics

    implicit none

    integer, parameter :: npts = 20
    integer, parameter :: n = 2*npts         ! BR and BZ for each point
    real, dimension(n) :: r, z, rn, zn, b
    integer :: i, ierr
    real :: theta
    vectype, dimension(npts) :: psi, br, bz

    ! calculate BR, BZ on LCFS; subract field from plasma current
    if(myrank.eq.0 .and. iprint.ge.1)  print *, 'Calculating vacuum field...'

    do i=1, npts
       theta = i*twopi/npts

       ! BR
       r(i) = sqrt(rzero**2 + 2.*ln*rzero*cos(theta))
       z(i) = rzero*ln*elongation*sin(theta) / r(i)
       b(i) = -ln*sin(theta)/(q0*rzero**2)
       rn(i) = 1.
       zn(i) = 0.

       !BZ
       r(i+npts) = r(i)
       z(i+npts) = z(i)
       b(i+npts) = ln*elongation* &
            (2.*rzero*cos(theta) + ln*(3.+cos(2.*theta))) &
            / (2.*q0*rzero**2*(rzero + 2.*ln*cos(theta)))
       rn(i+npts) = 0.
       zn(i+npts) = 1.
    end do

    ! Calculate field from plasma current
    call psi_vac(npts,r,z,psi,br,bz)
    b(1:npts) = b(1:npts) - br
    b(npts+1:n) = b(npts+1:n) - bz

!    call solve_bessel_fields_2D(n, r, z, rn, zn, b, ierr)
    call solve_multipole_fields_2D(n, r, z, rn, zn, b, ierr)
    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Done'
  end subroutine solovev_vac


  subroutine psi_vac(np,r,z,psi,br,bz)
    use basic
    use coils
    use math
    
    integer, intent(in) :: np
    real, intent(in), dimension(np) :: r, z
    vectype, intent(out), dimension(np) :: psi
    vectype, intent(out), dimension(np), optional :: br, bz

    real :: rmin, rmax, zmin, zmax
    real :: r0, z0, dr, dz, j0, jphi
    real, dimension(1) :: rp, zp
    real, dimension(np,1,6) :: g
    integer, parameter :: n = 200
    integer, parameter :: ipole = 0
    integer :: ierr, i, j

    j0 = -(1.+elongation**2)/(q0*rzero**3*elongation)

    rmin = sqrt(rzero**2 - 2.*ln*rzero)
    rmax = sqrt(rzero**2 + 2.*ln*rzero)
    dr = (rmax - rmin)/(n-1)
    
    psi = 0.
    if(present(br)) br = 0.
    if(present(bz)) bz = 0.
    do i=1, n
       r0 = (i-1)*dr + rmin
       rp = r0
       jphi = j0*r0
       zmax = elongation*sqrt((2.*ln*rzero)**2-(r0**2-rzero**2)**2)/(2.*r0)
       zmin = -zmax
       dz = (zmax-zmin)/(n-1)
       
       do j=1, n
          z0 = (j-1)*dz + zmin
          zp = z0

          call gvect(r,z,np,rp,zp,1,g,ipole,ierr)

          psi = psi - g(:,1,1) * jphi*dr*dz / twopi
          if(present(br)) br = br + (g(:,1,3)/r) * jphi*dr*dz / twopi
          if(present(bz)) bz = bz - (g(:,1,2)/r) * jphi*dr*dz / twopi
       end do
    end do
  end subroutine psi_vac

end module solovev

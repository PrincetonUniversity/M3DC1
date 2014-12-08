!==============================================================================
! Resistive wall test
!==============================================================================
module resistive_wall_test
  real, parameter, private :: a = .1
  real, private :: k, Im_ka

!!$  private :: analytic_response_matrix
contains

!!$  subroutine analytic_response_matrix
!!$    use math
!!$    use basic
!!$    use vacuum_interface
!!$
!!$    implicit none
!!$
!!$    real :: ka, bessk, besskp
!!$    real, allocatable :: theta(:)
!!$    integer :: ierr, i, j, m
!!$    complex :: fac, fac1, fac2, expimth
!!$
!!$    print *, 'Using analytic response matrix.'
!!$
!!$    call load_boundary_nodes(ierr)
!!$    if(ierr.ne.0) return
!!$
!!$    allocate(theta(nodes+1))
!!$    do i=1, nodes+1
!!$       theta(i) = atan2(znode(i)-zzero,xnode(i)-xzero)
!!$    end do
!!$
!!$    !...  define matrices with analytic formula
!!$    ka = a*ntor/rzero
!!$
!!$    zgrbth = 0.
!!$    zgrbph = 0.
!!$    zgrbphp = 0.
!!$    m = mpol
!!$    do m=mpol,mpol
!!$       fac =  (1./(nodes))*bessk(m,ka)/besskp(m,ka)
!!$       fac1 =  (0.,1.)*m/ka * fac
!!$       fac2 = +(0.,1.) * fac
!!$!      fac2 = -(0.,1.) * fac
!!$       do i=1,nodes+1
!!$          do j=1,nodes+1
!!$             expimth = cos(m*(theta(i) - theta(j))) + (0.,1.)*sin(m*(theta(i) - theta(j)))
!!$             zgrbth(i,j) = zgrbth(i,j) + fac1*expimth
!!$             zgrbph(i,j) = zgrbph(i,j) + fac2*expimth
!!$          enddo
!!$       enddo
!!$    end do
!!$
!!$    ! calculate derivatives (wrt i)
!!$    zgrbthp(1,:) = (zgrbth(2,:) - zgrbth(nodes+1,:)) &
!!$         /(theta(2) - theta(nodes+1))/a
!!$    zgrbphp(1,:) = (zgrbph(2,:) - zgrbph(nodes+1,:)) &
!!$         /(theta(2) - theta(nodes+1))/a
!!$    do i=2,nodes
!!$       zgrbthp(i,:) = (zgrbth(i+1,:) - zgrbth(i-1,:)) &
!!$            /(theta(i+1) - theta(i-1))/a 
!!$       zgrbphp(i,:) = (zgrbph(i+1,:) - zgrbph(i-1,:)) &
!!$            /(theta(i+1) - theta(i-1))/a
!!$    end do
!!$
!!$    deallocate(theta)
!!$
!!$  end subroutine analytic_response_matrix


subroutine resistive_wall_test_init()
  use basic
  use arrays
!!$  use vacuum_interface

  implicit none

  integer :: i,j, numnodes, icounter_t
  real :: x, phi, z, Km_ka, Imp_ka, Kmp_ka
  real :: bessi, bessk, bessip, besskp

  k = abs(ntor)/rzero
  Im_ka = bessi(mpol, k*a)
  Imp_ka = bessip(mpol, k*a)
  Km_ka = bessk(mpol, k*a)
  Kmp_ka = besskp(mpol, k*a)
  
  print *, 'analytic growth rate = ', &
       (eta_wall/delta_wall)*((mpol**2 + (k*a)**2)/(k*a)) * &
       ((Imp_ka*Km_ka - Im_ka*Kmp_ka)/(Imp_ka*Kmp_ka))/a

!!$  if(itaylor.eq.12 .and. eta_wall .gt. 0) call analytic_response_matrix
!!$
!!$  open(unit=97, file='response_poloidal', status='unknown')
!!$  do i=1, nodes
!!$     write(97,'(1p10e12.4)') (zgrbth(i,j),j=1,nodes+1)
!!$  end do
!!$  close(97)
!!$  open(unit=98, file='response_toroidal', status='unknown')
!!$  do i=1, nodes
!!$     write(98,'(1p10e12.4)') (zgrbph(i,j),j=1,nodes+1)
!!$  end do
!!$  close(98)



  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call get_node_pos(i, x, phi, z)

     call get_local_vals(i)

     x = x - xzero
     z = z - zzero

     call resistive_wall_test_equ(x, z)
     call resistive_wall_test_per(x, z)

     call set_local_vals(i)
  enddo

end subroutine resistive_wall_test_init

subroutine resistive_wall_test_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  psi0_l(1) = x**2 + z**2
  psi0_l(2) = 2.*x
  psi0_l(3) = 2.*z
  psi0_l(4) = 2.
  psi0_l(5) = 0.
  psi0_l(6) = 2.
    
  call constant_field(den0_l, 1.)
  if(numvar.le.1) return

  call constant_field(bz0_l , bzero*rzero)
  call constant_field(p0_l , p0)
  call constant_field(pe0_l, p0-pi0)

end subroutine resistive_wall_test_equ

subroutine resistive_wall_test_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: theta, r, ram, Im_kr, Imp_kr, Impp_kr
  real :: bessi
  complex :: fac, co, sn


  r = sqrt(x**2 + z**2)

  if(r.eq.0.) then
     psi1_l = 0.
     if(numvar.gt.1) bz1_l = 0.
  else
     theta = atan2(z, x)
     co = cos(mpol*theta) + (0.,1.)*sin(mpol*theta)
     sn = sin(mpol*theta) - (0.,1.)*cos(mpol*theta)
     ram = (r/a)**mpol

     psi1_l(1) = co*ram
     psi1_l(2) = mpol * (ram/r**2) * (x*co + z*sn)
     psi1_l(3) = mpol * (ram/r**2) * (z*co - x*sn)
     psi1_l(4) = mpol*(mpol-1) * (ram/r**4) * ((x**2 - z**2)*co + 2.*x*z*sn)
     psi1_l(5) = mpol*(mpol-1) * (ram/r**4) * (2.*x*z*co - (x**2 - z**2)*sn)
     psi1_l(6) = -psi1_l(4)
     psi1_l = psi1_l*eps
     
  if(numvar.le.1) return
     fac = (0.,1.)*k/Im_ka
     if(itor.eq.1) fac = fac*rzero
     Im_kr = bessi(mpol, k*r)
     Imp_kr = (bessi(mpol-1, k*r) + bessi(mpol+1, k*r))/2.
     Impp_kr = (bessi(mpol-2, k*r) + bessi(mpol+2, k*r) + 2.*Im_kr)/4.
   
     bz1_l(1) = -fac*sn*Im_kr
     bz1_l(2) =  (fac/r**2)*(mpol*z*co*Im_kr - k*r*x*sn*Imp_kr)
     bz1_l(3) = -(fac/r**2)*(mpol*x*co*Im_kr + k*r*z*sn*Imp_kr)
     bz1_l(4) = -(fac/r**4)*( mpol*z*(2.*x*co-mpol*z*sn)*Im_kr &
          - k*r*z*(2.*mpol*x*co-z*sn)*Imp_kr + (k*r*x)**2*sn*Impp_kr)
     bz1_l(5) = -(fac/r**4)*(mpol*(-(x**2-z**2)*co + mpol*x*z*sn)*Im_kr &
          + k*r*(mpol*(x**2-z**2)*co-x*z*sn)*Imp_kr + (k*r)**2*x*z*sn*Impp_kr)
     bz1_l(6) = -(fac/r**4)*(-mpol*x*(2.*z*co+mpol*x*sn)*Im_kr &
          + k*r*x*(2.*mpol*z*co+x*sn)*Imp_kr + (k*r*z)**2*sn*Impp_kr)
     
     bz1_l = bz1_l*eps
  end if

end subroutine resistive_wall_test_per

end module resistive_wall_test



!==============================================================================
! Circular shell for resistive wall test
!==============================================================================
module circ_shell_only

  private :: analytic_response_matrix
contains

!!$  subroutine analytic_response_matrix
!!$    use basic
!!$    use vacuum_interface
!!$
!!$    implicit none
!!$
!!$    real, allocatable :: theta(:)
!!$    real :: ka,  bessk, besskp, grate, a
!!$    integer :: ierr, i, j, m
!!$  complex :: fac
!!$
!!$    print *, 'Using analytic response matrix.'
!!$
!!$    call load_boundary_nodes(ierr)
!!$    allocate(theta(nodes+1))
!!$    do i=1, nodes+1
!!$       theta(i) = atan2(znode(i)-zzero,xnode(i)-xzero)
!!$    end do
!!$    if(ierr.ne.0) return
!!$
!!$    !...  define matrices with analytic formula
!!$    a = .1
!!$    ka = a*ntor/rzero
!!$
!!$    zgrbth = 0.
!!$    do m=mpol,mpol
!!$       fac = (0.,1.)*m/(nodes*ka)*bessk(m,ka)/besskp(m,ka)
!!$       do i=1,nodes+1
!!$          do j=1,nodes+1
!!$             zgrbth(i,j) = zgrbth(i,j) + fac*(cos(m*(theta(i) - theta(j)))  &
!!$                                   + (0.,1.)* sin(m*(theta(i) - theta(j))))
!!$          enddo
!!$       enddo
!!$    end do
!!$    ! calculate derivatives (wrt i)
!!$    zgrbthp(1,:) = (zgrbth(2,:) - zgrbth(nodes+1,:)) &
!!$         /(theta(2) - theta(nodes+1))/a
!!$    do i=2,nodes
!!$       zgrbthp(i,:) = (zgrbth(i+1,:) - zgrbth(i-1,:)) &
!!$            /(theta(i+1) - theta(i-1))/a
!!$    end do
!!$    grate = (eta_wall/delta_wall)*(mpol/a)* &
!!$         (1. - (mpol/ka)*bessk(mpol,ka)/besskp(mpol,ka))
!!$    write(*,'(A,1pe12.4)') " Analytic decay rate",  grate
!!$    write(97,'(A,1pe12.4)') " Analytic decay rate",  grate
!!$    deallocate(theta)
!!$
!!$  end subroutine analytic_response_matrix


  subroutine circ_shell_only_init()
    use basic
    use arrays
!!$    use vacuum_interface
    
    implicit none
    
    integer :: i,j, numnodes, icounter_t
    real :: x, phi, z
    
!!$    open(unit=97,file="response_matrix",status="unknown")
!!$    if(itaylor.eq.10 .and. eta_wall .gt. 0) call analytic_response_matrix
!!$    
!!$    do i=1, nodes
!!$       write(97,'(1p10e12.4)') (zgrbth(i,j),j=1,nodes+1)
!!$    end do
!!$    close(97)
    
    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       i = nodes_owned(icounter_t)
       call get_node_pos(i, x, phi, z)
       
       call get_local_vals(i)
       
       x = x - xzero
       z = z - zzero
       
       call circ_shell_only_equ(x, z)
       call circ_shell_only_per(x, z)
       
       call set_local_vals(i)
    enddo
    
  end subroutine circ_shell_only_init
  
  subroutine circ_shell_only_equ(x, z)
    use basic
    use arrays
    
    implicit none
    
    real, intent(in) :: x, z

    
    call constant_field(den0_l, 1.)
    u0_l = 0.
    
    psi0_l(1) = x**2 + z**2
    psi0_l(2) = 2.*x
    psi0_l(3) = 2.*z
    psi0_l(4) = 2.
    psi0_l(5) = 0.
    psi0_l(6) = 2.
    
    if(numvar.le.1) return
    
    vz0_l = 0.
    call constant_field(bz0_l , bzero)
    
    if(numvar.le.2) return
    chi0_l = 0.
    call constant_field(p0_l  , p0)
    call constant_field(pe0_l , p0-pi0)
    
  end subroutine circ_shell_only_equ
  
  
  subroutine circ_shell_only_per(x, z)
    use basic
    use arrays
    
    implicit none
    
    real, intent(in) :: x, z
    real ::  k, r, kr, theta
    complex :: cosm, sinm

    den1_l = 0.
    
    u1_l = 0.
    psi1_l = 0.
    r = sqrt(x**2 + z**2)
    if(r.eq.0) return
    k = ntor/rzero
    kr = k*r
    theta = atan2(z,x)
    cosm = cos(mpol*theta) + (0.,1.)*sin(mpol*theta)
    sinm = - (0.,1.)*cosm

    !..for del_perp_squared
    psi1_l(1) = eps*cosm * r**mpol
    psi1_l(2) = eps*(( z*sinm + x*cosm)*mpol*r**(mpol-2))
    psi1_l(3) = eps*((-x*sinm + z*cosm)*mpol*r**(mpol-2))
    psi1_l(4) = eps*((2*x*z*sinm + (x**2 - z**2)*cosm)*mpol*(mpol-1)*r**(mpol-4))
    psi1_l(5) = eps*((2*x*z*cosm + (z**2 - x**2)*sinm)*mpol*(mpol-1)*r**(mpol-4))
    psi1_l(6) = eps*((-2*x*z*sinm + (z**2 - x**2)*cosm)*mpol*(mpol-1)*r**(mpol-4))    
    
  end subroutine circ_shell_only_per

end module circ_shell_only

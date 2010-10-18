module element

  ! The following give the meaning of each dof at one node
  integer, parameter :: DOF_1   = 1
  integer, parameter :: DOF_DR  = 2
  integer, parameter :: DOF_DZ  = 3
  integer, parameter :: DOF_DRR = 4
  integer, parameter :: DOF_DRZ = 5
  integer, parameter :: DOF_DZZ = 6
#ifdef USE3D
  integer, parameter :: DOF_DP   = 7
  integer, parameter :: DOF_DRP  = 8
  integer, parameter :: DOF_DZP  = 9
  integer, parameter :: DOF_DRRP = 10
  integer, parameter :: DOF_DRZP = 11
  integer, parameter :: DOF_DZZP = 12
#endif

  ! The dofs per element are the concatenation of the dofs per node
  ! for each node in the element

  integer, parameter :: maxpol = 3
  integer, parameter :: pol_dofs_per_node = 6
  integer, parameter :: pol_nodes_per_element = 3
#ifdef USE3D
  integer, parameter :: tor_dofs_per_node = 2
  integer, parameter :: tor_nodes_per_element = 2
  integer, parameter :: maxtor = 3
  integer, parameter :: edges_per_element = 9
  integer, parameter :: coeffs_per_dphi = 4
  integer, parameter :: dofs_per_dphi = 4
#else
  integer, parameter :: tor_dofs_per_node = 1
  integer, parameter :: tor_nodes_per_element = 1
  integer, parameter :: maxtor = 1
  integer, parameter :: edges_per_element = 3
  integer, parameter :: coeffs_per_dphi = 1
  integer, parameter :: dofs_per_dphi = 1
#endif
  integer, parameter :: dofs_per_node = tor_dofs_per_node*pol_dofs_per_node
  integer, parameter :: nodes_per_element = &
       pol_nodes_per_element*tor_nodes_per_element
  integer, parameter :: coeffs_per_tri = 20
  integer, parameter :: dofs_per_tri = 18

  integer, parameter :: nodes_per_edge = 2
  integer, parameter :: dofs_per_element = nodes_per_element*dofs_per_node
  integer, parameter :: coeffs_per_element = coeffs_per_tri*coeffs_per_dphi

  type element_data
     real :: R, Phi, Z, a, b, c, d, co, sn, itri
  end type element_data

  integer :: ni(coeffs_per_tri),mi(coeffs_per_tri)  
  data mi /0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0/
  data ni /0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5/
#ifdef USE3D
  integer :: li(coeffs_per_dphi)
  data li /0,1,2,3/
#endif

  real, allocatable :: gtri(:,:,:), gtri_old(:,:,:)
  real, allocatable :: htri(:,:,:)

contains

  
  !=======================================================
  ! global_to_local
  ! ~~~~~~~~~~~~~~~
  ! transforms from global coordinates 
  ! to local (element) coordinates
  !=======================================================
  elemental subroutine global_to_local(d, R, Phi, Z, xi, zi, eta)
    implicit none

    type(element_data), intent(in) :: d
    real, intent(in) :: R, Phi, Z
    real, intent(out) :: xi, zi, eta

    xi  =  (R-d%R)*d%co + (Z-d%Z)*d%sn - d%b
    eta = -(R-d%R)*d%sn + (Z-d%Z)*d%co
    zi  =  Phi - d%Phi
  end subroutine global_to_local

  !=======================================================
  ! local_to_global
  ! ~~~~~~~~~~~~~~~
  ! transforms from local (element) coordinates
  ! to global coordinates
  !=======================================================   
  elemental subroutine local_to_global(d, xi, zi, eta, R, Phi, Z)
    implicit none

    type(element_data), intent(in) :: d
    real, intent(in) :: xi, zi, eta
    real, intent(out) :: R, Phi, Z

    R = d%R + (d%b+xi)*d%co - eta*d%sn
    Z = d%Z + (d%b+xi)*d%sn + eta*d%co
    Phi = d%Phi + zi
  end subroutine local_to_global

  logical function is_in_element(d, R, phi, z)
    implicit none
    type(element_data), intent(in) :: d
    real, intent(in) :: R, Phi, Z
    real :: f, xi, zi, eta

    call global_to_local(d, R, Phi, Z, xi, zi, eta)

    is_in_element = .false.
    if(eta.lt.0.) return
    if(eta.gt.d%c) return

    f = 1. - eta/d%c
    if(xi.lt.-f*d%b) return
    if(xi.gt. f*d%a) return
    
#ifdef USE3D
    if(zi.le.0.) return
    if(zi.ge.d%d) return
#endif

    is_in_element = .true.
  end function is_in_element


  !======================================================================
  ! rotate_dofs
  ! ~~~~~~~~~~~
  !
  ! Performs coordinate rotation from (R, Z) to (n, t) on invec,
  ! returns result in outvec.
  !======================================================================
  subroutine rotate_dofs(invec, outvec, normal, curv, ic)
    implicit none

    real, intent(in) :: curv, normal(2) 
    integer, intent(in) :: ic
    
    vectype, dimension(dofs_per_node), intent(in) :: invec
    vectype, dimension(dofs_per_node), intent(out) :: outvec

    if(ic.eq.1) then
       outvec(1) = invec(1)
       outvec(2) = normal(1)*invec(2) + normal(2)*invec(3)
       outvec(3) = normal(1)*invec(3) - normal(2)*invec(2)
       outvec(4) = normal(1)**2*invec(4) + normal(2)**2*invec(6) &
            + 2.*normal(1)*normal(2)*invec(5)
       outvec(5) = (normal(1)**2 - normal(2)**2)*invec(5) &
            + normal(1)*normal(2)*(invec(6) - invec(4)) &
            + curv*outvec(3)
       outvec(6) = normal(1)**2*invec(6) + normal(2)**2*invec(4) &
            - 2.*normal(1)*normal(2)*invec(5) &
            - curv*outvec(2)
#ifdef USE3D
       outvec(7) = invec(7)
       outvec(8) = normal(1)*invec(8) + normal(2)*invec(9)
       outvec(9) = normal(1)*invec(9) - normal(2)*invec(8)
       outvec(10) = normal(1)**2*invec(10) + normal(2)**2*invec(12) &
            + 2.*normal(1)*normal(2)*invec(11)
       outvec(11) = (normal(1)**2 - normal(2)**2)*invec(11) &
            + normal(1)*normal(2)*(invec(12) - invec(10)) &
            + curv*outvec(9)
       outvec(12) = normal(1)**2*invec(12) + normal(2)**2*invec(10) &
            - 2.*normal(1)*normal(2)*invec(11) &
            - curv*outvec(8)
#endif
    else if (ic.eq.-1) then
       outvec(1) = invec(1)
       outvec(2) = normal(1)*invec(2) - normal(2)*invec(3) &
            - curv*normal(2)*invec(5) - curv*normal(1)*invec(6)
       outvec(3) = normal(2)*invec(2) + normal(1)*invec(3) &
            + curv*normal(1)*invec(5) - curv*normal(2)*invec(6)
       outvec(4) = normal(1)**2*invec(4) + normal(2)**2*invec(6) &
            - normal(1)*normal(2)*invec(5)
       outvec(5) = 2.*normal(1)*normal(2)*(invec(4) - invec(6)) &
            + (normal(1)**2 - normal(2)**2)*invec(5)
       outvec(6) = normal(2)**2*invec(4) + normal(1)**2*invec(6) &
               + normal(2)*normal(2)*invec(5)
#ifdef USE3D
       outvec(7) = invec(7)
       outvec(8) = normal(1)*invec(8) - normal(2)*invec(9) &
            - curv*normal(2)*invec(11) - curv*normal(1)*invec(12)
       outvec(9) = normal(2)*invec(8) + normal(1)*invec(9) &
            + curv*normal(1)*invec(11) - curv*normal(2)*invec(12)
       outvec(10) = normal(1)**2*invec(10) + normal(2)**2*invec(12) &
            - normal(1)*normal(2)*invec(11)
       outvec(11) = 2.*normal(1)*normal(2)*(invec(10) - invec(12)) &
            + (normal(1)**2 - normal(2)**2)*invec(11)
       outvec(12) = normal(2)**2*invec(10) + normal(1)**2*invec(12) &
               + normal(2)*normal(2)*invec(11)
#endif
    else
       outvec(1) = invec(1)
       outvec(2) = normal(1)*invec(2) + normal(2)*invec(3) &
            + curv*normal(2)**2*invec(4) &
            - curv*normal(1)*normal(2)*invec(5) &
            + curv*normal(1)**2*invec(6)
       outvec(3) = normal(1)*invec(3) - normal(2)*invec(2) &
            + 2.*curv*normal(1)*normal(2)*(invec(4) - invec(6)) &
            - curv*(normal(1)**2 - normal(2)**2)*invec(5)
       outvec(4) = normal(1)**2*invec(4) + normal(2)**2*invec(6) &
            + normal(1)*normal(2)*invec(5)
       outvec(5) = (normal(1)**2 - normal(2)**2)*invec(5) &
            + 2.*normal(1)*normal(2)*(invec(6) - invec(4))
       outvec(6) = normal(1)**2*invec(6) + normal(2)**2*invec(4) &
            - normal(1)*normal(2)*invec(5)
#ifdef USE3D
       outvec(7) = invec(7)
       outvec(8) = normal(1)*invec(8) + normal(2)*invec(9) &
            + curv*normal(2)**2*invec(10) &
            - curv*normal(1)*normal(2)*invec(11) &
            + curv*normal(1)**2*invec(12)
       outvec(9) = normal(1)*invec(9) - normal(2)*invec(8) &
            + 2.*curv*normal(1)*normal(2)*(invec(10) - invec(12)) &
            - curv*(normal(1)**2 - normal(2)**2)*invec(11)
       outvec(10) = normal(1)**2*invec(10) + normal(2)**2*invec(12) &
            + normal(1)*normal(2)*invec(11)
       outvec(11) = (normal(1)**2 - normal(2)**2)*invec(11) &
            + 2.*normal(1)*normal(2)*(invec(12) - invec(10))
       outvec(12) = normal(1)**2*invec(12) + normal(2)**2*invec(10) &
            - normal(1)*normal(2)*invec(11)
#endif
    endif
  end subroutine rotate_dofs

  !============================================================
  ! tmatrix
  ! ~~~~~~~
  ! define the 20 x 20 Transformation Matrix that enforces the condition that
  ! the nomal slope between triangles has only cubic variation..
  !============================================================
  subroutine tmatrix(ti,ndim,a,b,c)
    implicit none
    integer, intent(in) :: ndim
    real, intent(out) :: ti(ndim,*)
    real, intent(in) :: a, b, c
    
    integer :: ifail, info1, info2
    real :: danaly, det, percent, diff
    real :: t(coeffs_per_tri,coeffs_per_tri), wkspce(9400)
    integer :: ipiv(coeffs_per_tri)
    
    integer, parameter :: ierrorchk = 0
    
    ! first initialize to zero
    t = 0
    
    ! Table 1 of Ref. [2]
    t(1,1)   = 1.
    t(1,2)   = -b
    t(1,4)   = b**2
    t(1,7)   = -b**3
    t(1,11)  = b**4
    t(1,16)  = -b**5
    
    t(2,2)   = 1
    t(2,4)   = -2*b
    t(2,7)   = 3*b**2
    t(2,11)  = -4*b**3
    t(2,16)  = 5*b**4
    
    t(3,3)   = 1
    t(3,5)   = -b
    t(3,8)   = b**2
    t(3,12)  = -b**3
    
    t(4,4)   = 2.
    t(4,7)   = -6.*b
    t(4,11)  = 12*b**2
    t(4,16)  = -20*b**3
    
    t(5,5)   = 1.
    t(5,8)   = -2.*b
    t(5,12)  = 3*b**2
    
    t(6,6)   = 2.
    t(6,9)   = -2*b
    t(6,13)  = 2*b**2
    t(6,17)  = -2*b**3
    
    t(7,1)   = 1.
    t(7,2)   = a
    t(7,4)   = a**2
    t(7,7)   = a**3
    t(7,11)  = a**4
    t(7,16)  = a**5
    
    t(8,2)   = 1.
    t(8,4)   = 2*a
    t(8,7)   = 3*a**2
    t(8,11)  = 4*a**3
    t(8,16)  = 5*a**4
    
    t(9,3)   = 1.
    t(9,5)   = a
    t(9,8)   = a**2
    t(9,12)  = a**3
    
    t(10,4)  = 2
    t(10,7)  = 6*a
    t(10,11) = 12*a**2
    t(10,16) = 20*a**3
    
    t(11,5)  = 1.
    t(11,8)  = 2.*a
    t(11,12) = 3*a**2
    
    t(12,6)  = 2.
    t(12,9)  = 2*a
    t(12,13) = 2*a**2
    t(12,17) = 2*a**3
    
    t(13,1)  = 1
    t(13,3)  = c
    t(13,6)  = c**2
    t(13,10) = c**3
    t(13,15) = c**4
    t(13,20) = c**5
    
    t(14,2)  = 1.
    t(14,5)  = c
    t(14,9)  = c**2
    t(14,14) = c**3
    t(14,19) = c**4
    
    t(15,3)  = 1.
    t(15,6)  = 2*c
    t(15,10) = 3*c**2
    t(15,15) = 4*c**3
    t(15,20) = 5*c**4
    
    t(16,4)  = 2.
    t(16,8)  = 2*c
    t(16,13) = 2*c**2
    t(16,18) = 2*c**3
    
    t(17,5)  = 1.
    t(17,9)  = 2*c
    t(17,14) = 3*c**2
    t(17,19) = 4*c**3
    
    t(18,6)  = 2.
    t(18,10) = 6*c
    t(18,15) = 12*c**2
    t(18,20) = 20*c**3
    
    t(19,16) = 5*a**4*c
    t(19,17) = 3*a**2*c**3 - 2*a**4*c
    t(19,18) = -2*a*c**4+3*a**3*c**2
    t(19,19) = c**5-4*a**2*c**3
    t(19,20) = 5*a*c**4
    
    t(20,16) = 5*b**4*c
    t(20,17) = 3*b**2*c**3 - 2*b**4*c
    t(20,18) = 2*b*c**4 - 3*b**3*c**2
    t(20,19) = c**5 - 4*b**2*c**3
    t(20,20) = -5*b*c**4
    
    if(ierrorchk.eq.1) then
       ! analytic formula for determinant
       danaly = -64*(a+b)**17*c**20*(a**2+c**2)*(b**2+c**2)
       
       ! calculate determinant using nag
       ifail = 0
       ti(1:20,1:20) = t
       det = 0.
       !     call f03aaf(ti,20,20,det,wkspce,ifail)
       
       diff = det - danaly
       percent = 100* diff / danaly
       
       if(abs(percent) .gt. 1.e-12) &
            print *, "percent error in determinant: ", percent
    endif
    
    ! calculate the inverse of t using NAG library routines
    info1 = 0
    info2 = 0
    ti(1:20,1:20) = t
    call f07adf(20,20,ti,20,ipiv,info1)
    call f07ajf(20,ti,20,ipiv,wkspce,400,info2)
    if(info1.ne.0.or.info2.ne.0) write(*,'(3I5)') info1,info2
    
  end subroutine tmatrix


  !======================================================================
  ! local_coeff_vector
  ! ~~~~~~~~~~~~~~~~~~
  ! calculates the coefficients of the polynomial expansion of the
  ! field in the element domain
  !======================================================================
  subroutine local_coeff_vector(itri, c, iold)
    implicit none

    integer, intent(in) :: itri
    real, intent(out), dimension(dofs_per_element,coeffs_per_element) :: c
    logical, intent(in) :: iold

    integer :: i, j, k, l, m, n
    integer :: idof, icoeff, ip, it

    c = 0.

    icoeff = 0
    do i=1,coeffs_per_dphi
       do j=1,coeffs_per_tri
          icoeff = icoeff + 1
          idof = 0
          it = 0
          do k=1,tor_nodes_per_element
             do l=1,pol_nodes_per_element
                do m=1,tor_dofs_per_node
                   do n=1,pol_dofs_per_node
                      idof = idof + 1
                      ip = n + (l-1)*pol_dofs_per_node
                      it = m + (k-1)*tor_dofs_per_node
                      if(iold) then
                         c(idof,icoeff) = c(idof,icoeff) &
                              + htri(i,it,itri)*gtri_old(j,ip,itri)
                      else
                         c(idof,icoeff) = c(idof,icoeff) &
                              + htri(i,it,itri)*gtri(j,ip,itri)
                      endif
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine local_coeff_vector


  !======================================================================
  ! local_coeffs
  ! ~~~~~~~~~~~~
  ! calculates the coefficients of the polynomial expansion of the
  ! field in the element domain
  !======================================================================
  subroutine local_coeffs(itri, dof, c)
    implicit none

    integer, intent(in) :: itri
    vectype, intent(in), dimension(dofs_per_element) :: dof
    vectype, intent(out), dimension(coeffs_per_element) :: c

    real, dimension(dofs_per_element,coeffs_per_element) :: cl
    integer :: i, j

    call local_coeff_vector(itri, cl, .true.)
    c = 0.
    do i=1, coeffs_per_element
       do j=1, dofs_per_element
          c(i) = c(i) + cl(j,i)*dof(j)
       end do
    end do
  end subroutine local_coeffs


!!$  !======================================================================
!!$  ! local_value
!!$  ! ~~~~~~~~~~~
!!$  ! calculates the value of a field at n points within an element
!!$  ! v(m,i,j,k) is d^(i-1)_R d^(j-1)_Phi d^(k-1)_Z of the field at point m
!!$  !======================================================================
!!$  subroutine local_value(d, n, xi, zi, eta, dpol, dtor, v)
!!$    type(element_data), intent(in) :: d
!!$    integer, intent(in) :: n
!!$    vectype, intent(in), dimension(dofs_per_element) :: dof
!!$    real, intent(in), dimension(n) :: xi, zi, eta
!!$    vectype, intent(out), dimension(n, maxpol, maxtor, maxpol) :: v
!!$
!!$    vectype, dimension(coeffs_per_element) :: c
!!$
!!$    integer :: p,q,i,j,k
!!$    real :: co, sn, co2, cosn, sn2
!!$    real, dimension(maxpol, maxpol) :: lval
!!$    real, dimension(maxpol, maxtor, maxpol) :: val
!!$
!!$    call local_coeffs(dof, c)
!!$
!!$    ! need inverse rotation to get from local to global coords
!!$    co =  d%co
!!$    sn = -d%sn
!!$
!!$    co2 = co*co
!!$    cosn = co*sn
!!$    sn2 = sn*sn
!!$
!!$    do k=1, n
!!$       do p=1, coeffs_per_tri
!!$          val = 0.
!!$          
!!$          ! calculate values in local coordinates
!!$          lval(0,0) = xi(k)**mi(p) * eta(k)**ni(p)
!!$          
!!$          if(dpol.ge.1) then
!!$             if(mi(p).ge.1) then
!!$                ! d_xi terms
!!$                lval(1,0) = mi(p)*xi(k)**(mi(p)-1) * eta(k)**ni(p)
!!$             endif
!!$             if(ni(p).ge.1) then
!!$                ! d_eta terms
!!$                lval(0,1) = xi(k)**mi(p) * eta(k)**(ni(p)-1)*ni(p)
!!$             endif
!!$          endif
!!$          
!!$          if(dpol.ge.2) then
!!$             if(mi(p).ge.2) then
!!$                ! d_xi^2 terms
!!$                lval(2,0) = xi(k)**(mi(p)-2)*(mi(p)-1)*mi(p) * eta(k)**ni(p)
!!$             endif
!!$             
!!$             if(ni(p).ge.2) then
!!$                ! d_eta^2 terms
!!$                lval(0,2) = xi(k)**mi(p) * eta(k)**(ni(p)-2)*(ni(p)-1)*ni(p)
!!$             endif
!!$             
!!$             if(mi(p).ge.1 .and. ni(p).ge.1) then
!!$                ! d_xi d_eta terms
!!$                lval(1,1) = xi(k)**(mi(p)-1)*mi(p) * eta(k)**(ni(p)-1)*ni(p)
!!$             endif
!!$          endif
!!$
!!$          if(dpol.ge.3) then
!!$             if(mi(p).ge.3) then
!!$                ! d_xi^3 terms
!!$                lval(3,0) = xi(k)**(mi(p)-3)*(mi(p)-2)*(mi(p)-1)*mi(p) &
!!$                     *     eta(k)**ni(p)
!!$             endif
!!$             if(mi(p).ge.2 .and. ni(p).ge.1) then
!!$                ! d_xi^2 d_eta terms
!!$                lval(2,1) = xi(k)**(mi(p)-2)*(mi(p)-1)*mi(p) &
!!$                     *     eta(k)**(ni(p)-1)* ni(p)
!!$             endif
!!$             if(mi(p).ge.1 .and. ni(p).ge.2) then
!!$                ! d_xi d_eta^2 terms
!!$                lval(1,2) = xi(k)**(mi(p)-1)*mi(p) &
!!$                     *     eta(k)**(ni(p)-2)*(ni(p)-1)*ni(p)
!!$             endif
!!$             if(ni(p).ge.3) then
!!$                ! d_eta^3 terms
!!$                lval(0,3) = xi(k)**mi(p) &
!!$                     *     eta(k)**(ni(p)-3)*(ni(p)-2)*(ni(p)-1)*ni(p)
!!$             endif             
!!$          endif
!!$
!!$
!!$          ! rotate values to global coordinates
!!$          val(0,0,0) = lval(0,0)
!!$          val(1,0,0) = co*lval(1,0) + sn*lval(0,1)
!!$          val(0,0,1) = co*lval(0,1) - sn*lval(1,0)
!!$          val(2,0,0) = co2*lval(2,0) + sn2*lval(0,2) + 2.*cosn*lval(1,1)
!!$          val(1,0,1) = (co2 - sn2)*lval(1,1) + cosn*(lval(0,2) - lval(2,0))
!!$          val(0,0,2) = co2*lval(0,2) + sn2*lval(2,0) - 2.*cosn*lval(1,1)
!!$
!!$          ! NEED TO INCLUDE ROTATION OF 3RD DERIVATIVE TERMS HERE
!!$
!!$          ! include toroidal derivatives
!!$          do q=1, coeffs_per_dphi
!!$#ifdef USE3D
!!$             if(dtor.ge.3) then
!!$                if(li(q).ge.3) then
!!$                   val(:,3,:) = val(:,0,:) &
!!$                        *zi**(li(q)-3)*(li(q)-2)*(li(q)-1)*li(q)
!!$                else
!!$                   val(:,3,:) = 0.
!!$                endif
!!$             endif
!!$             if(dtor.ge.2) then
!!$                if(li(q).ge.2) then
!!$                   val(:,2,:) = val(:,0,:)*zi**(li(q)-2)*(li(q)-1)*li(q)
!!$                else
!!$                   val(:,2,:) = 0.
!!$                endif
!!$             endif
!!$             if(dtor.ge.1) then
!!$                if(li(q).ge.1) then
!!$                   val(:,1,:) = val(:,0,:)*zi**(li(q)-1)*li(q)
!!$                else
!!$                   val(:,1,:) = 0.
!!$                endif
!!$             endif
!!$             val(:,0,:) = val(:,0,:)*zi**li(q)
!!$#endif
!!$             v(k,:,:,:) = v(k,:,:,:) + c(j)*val(:,:,:)
!!$
!!$             j = j + 1
!!$          end do
!!$       end do
!!$    end do
!!$    
!!$  end subroutine local_value
  
end module element

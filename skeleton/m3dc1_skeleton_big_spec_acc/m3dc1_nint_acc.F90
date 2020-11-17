!======================================================================
! m3dc1_nint
! ~~~~~~~~~~
! This module defines arrays of the values of fields at numerical
! integration quadrature sampling points, and routines for populating
! these arrays.
!======================================================================
module m3dc1_nint
  use nintegrate

  implicit none

  ! The following give the meaning of the value array returned by local_value
  integer, parameter :: OP_1    = 1
  integer, parameter :: OP_DR   = 2
  integer, parameter :: OP_DZ   = 3
  integer, parameter :: OP_DRR  = 4
  integer, parameter :: OP_DRZ  = 5
  integer, parameter :: OP_DZZ  = 6
  integer, parameter :: OP_LP   = 7
  integer, parameter :: OP_GS   = 8
  integer, parameter :: OP_NUM_POL = 8
#if defined(USE3D) 
  integer, parameter :: OP_DP    = 9
  integer, parameter :: OP_DRP   = 10
  integer, parameter :: OP_DZP   = 11
  integer, parameter :: OP_DRRP  = 12
  integer, parameter :: OP_DRZP  = 13
  integer, parameter :: OP_DZZP  = 14
  integer, parameter :: OP_LPP   = 15
  integer, parameter :: OP_GSP   = 16
  integer, parameter :: OP_DPP   = 17
  integer, parameter :: OP_DRPP  = 18
  integer, parameter :: OP_DZPP  = 19
  integer, parameter :: OP_DRRPP = 20
  integer, parameter :: OP_DRZPP = 21
  integer, parameter :: OP_DZZPP = 22
  integer, parameter :: OP_LPPP  = 23
  integer, parameter :: OP_GSPP  = 24
  integer, parameter :: OP_LPR   = 25
  integer, parameter :: OP_LPZ   = 26
  integer, parameter :: OP_NUM   = 26
#else
  integer, parameter :: OP_LPR  = 9
  integer, parameter :: OP_LPZ  = 10
  integer, parameter :: OP_NUM  = 10
#endif

  integer, parameter :: nfield = 20

  real, dimension(MAX_PTS, OP_NUM, dofs_per_element) :: mu79, nu79
  real, dimension(MAX_PTS, OP_NUM) :: eta79, psi79
  real, dimension(MAX_PTS, OP_NUM, nfield) :: f79
!$OMP THREADPRIVATE(mu79,nu79,eta79,psi79,f79)

  !$acc declare create (nu79, mu79)
  ! precalculated terms
  !real, private :: fterm(MAX_PTS, OP_NUM, coeffs_per_element)
  !acc declare create (fterm)
!$OMP THREADPRIVATE(fterm)
contains

!==================================================
! precalculate_terms
! ~~~~~~~~~~~~~~~~~~
! Precalculates the values of each term in the 
! finite element expansion at each sampling point
!==================================================
subroutine precalculate_terms(xi,zi,eta,co,sn,fterm)
  use element

!$acc routine

  implicit none
  real :: fterm(MAX_PTS, OP_NUM, coeffs_per_element)
  real, intent(in) :: co, sn
  real, dimension(MAX_PTS), intent(in) :: xi, zi, eta

  integer :: i,j,p,op
  real, dimension(MAX_PTS) :: temp
  real :: xpow(MAX_PTS,-3:5), ypow(MAX_PTS,-3:5)
#ifdef USE3D
  real :: zpow(MAX_PTS,-2:3)
#endif
  real :: co2, sn2, cosn

  co2 = co*co
  sn2 = sn*sn
  cosn = co*sn

  ! precalculate powers
  xpow(:,-3:-1) = 0.
  ypow(:,-3:-1) = 0.
  xpow(:,0) = 1.
  ypow(:,0) = 1.
#ifdef USE3D
  zpow(:,-2:-1) = 0.
  zpow(:,0) = 1.
#endif

  do p=1, 5
     xpow(:,p) = xpow(:,p-1)*xi(:)
     ypow(:,p) = ypow(:,p-1)*eta(:)
  end do
#ifdef USE3D
  do p=1, 3
     zpow(:,p) = zpow(:,p-1)*zi(:)
  end do
#endif


  fterm = 0.

  do p=1, coeffs_per_tri

     fterm(:,OP_1,p) = xpow(:,mi(p))*ypow(:,ni(p))
          
     if(mi(p).ge.1) then
        ! d_si terms
        temp = mi(p)*xpow(:,mi(p)-1) * ypow(:,ni(p))
        fterm(:,OP_DR,p) = fterm(:,OP_DR,p) + co*temp
        fterm(:,OP_DZ,p) = fterm(:,OP_DZ,p) + sn*temp           
             
        if(mi(p).ge.2) then
           ! d_si^2 terms
           temp = xpow(:,mi(p)-2)*(mi(p)-1)*mi(p) * ypow(:,ni(p))
           fterm(:,OP_DRR,p) = fterm(:,OP_DRR,p) + co2*temp
           fterm(:,OP_DZZ,p) = fterm(:,OP_DZZ,p) + sn2*temp
           fterm(:,OP_DRZ,p) = fterm(:,OP_DRZ,p) + cosn*temp
           fterm(:,OP_LP,p) = fterm(:,OP_LP,p) + temp
        endif
     endif
     if(ni(p).ge.1) then
        ! d_eta terms
        temp = xpow(:,mi(p)) * ypow(:,ni(p)-1)*ni(p)
        fterm(:,OP_DR,p) = fterm(:,OP_DR,p) - sn*temp
        fterm(:,OP_DZ,p) = fterm(:,OP_DZ,p) + co*temp
             
        if(ni(p).ge.2) then
           ! d_eta^2 terms
           temp = xpow(:,mi(p)) * ypow(:,ni(p)-2)*(ni(p)-1)*ni(p)
           fterm(:,OP_DRR,p) = fterm(:,OP_DRR,p) + sn2*temp
           fterm(:,OP_DZZ,p) = fterm(:,OP_DZZ,p) + co2*temp
           fterm(:,OP_DRZ,p) = fterm(:,OP_DRZ,p) - cosn*temp
           fterm(:,OP_LP,p) = fterm(:,OP_LP,p) + temp
        endif
             
        if(mi(p).ge.1) then
           ! d_eta_si terms
           temp = xpow(:,mi(p)-1)*mi(p) * ypow(:,ni(p)-1)*ni(p)
                
           fterm(:,OP_DRR,p) = fterm(:,OP_DRR,p) - 2.*cosn*temp
           fterm(:,OP_DZZ,p) = fterm(:,OP_DZZ,p) + 2.*cosn*temp
           fterm(:,OP_DRZ,p) = fterm(:,OP_DRZ,p) + (co2-sn2)*temp
        endif
     endif
                   
     ! Grad-Shafranov operator, and
     ! cylindrical correction to Laplacian
     fterm(:,OP_GS,p) = fterm(:,OP_LP,p)

#ifdef USE3D
     do op=1, OP_NUM_POL
        temp = fterm(:,op,p)
        do i=1, coeffs_per_dphi
           j = p + (i-1)*coeffs_per_tri

           fterm(:,op,j) = temp(:)*zpow(:,li(i))

           ! first toroidal derivative
           if(li(i).ge.1) then
              fterm(:,op+OP_NUM_POL,j) = temp(:) &
                   *zpow(:,li(i)-1)*li(i)
           endif
           ! second toroidal derivative
           if(li(i).ge.2) then
              fterm(:,op+2*OP_NUM_POL,j) = temp(:) &
                   *zpow(:,li(i)-2)*li(i)*(li(i)-1)
           endif
        end do
     end do
#endif
  end do
end subroutine precalculate_terms

  subroutine define_basis(itri,fterm,mu79_l,nu79_l)
    implicit none

    !$acc routine
    !$acc routine (local_coeff_vector)
    !$acc routine (dgemm) bind(d_dgemm)

    integer, intent(in) :: itri
    real, dimension(MAX_PTS, OP_NUM, dofs_per_element) :: nu79_l, mu79_l
    real :: fterm(MAX_PTS, OP_NUM, coeffs_per_element)
    real, dimension(dofs_per_element,coeffs_per_element) :: cl
    real, dimension(MAX_PTS*OP_NUM, dofs_per_element) :: temp

    integer :: i, p, op

    call local_coeff_vector(itri, cl, .false.)
    call dgemm('N','T',MAX_PTS*OP_NUM,dofs_per_element,dofs_per_element,&
               1.,fterm,MAX_PTS*OP_NUM,cl,dofs_per_element,0.,mu79_l,MAX_PTS*OP_NUM)

    nu79_l = mu79_l
  end subroutine define_basis

  !===============================================
  ! eval_ops
  ! --------
  !
  ! evaluates linear unitary operators
  !===============================================
  subroutine eval_ops(itri,fin,outarr,nu79_l)
    use field
    implicit none
      
    real, dimension(MAX_PTS, OP_NUM), intent(out) :: outarr
    real, dimension(MAX_PTS, OP_NUM, dofs_per_element) :: nu79_l
    integer :: i, op

    integer, intent(in) :: itri
    type(field_type), intent(in) :: fin

    real, dimension(dofs_per_element) :: dofs

    !$acc routine
    !$acc routine (get_element_dofs)
    !$acc routine (dgemv) bind(d_dgemv)

    call get_element_dofs(fin, itri, dofs)
  
    call dgemv('N',MAX_PTS*OP_NUM,dofs_per_element,1.,nu79_l,MAX_PTS*OP_NUM,dofs,1,0.,outarr,1)
  end subroutine eval_ops

  !===============================================
  ! eval_ops
  ! --------
  !
  ! evaluates linear unitary operators
  !===============================================
  subroutine eval_ops_old(avector,ngauss,outarr,fterm)
    implicit none
      
    integer, intent(in) :: ngauss
    real, dimension(coeffs_per_element), intent(in) :: avector
    real, dimension(ngauss, OP_NUM), intent(out) :: outarr
    real :: fterm(MAX_PTS, OP_NUM, coeffs_per_element)

    integer :: p, op

    outarr = 0.
    do op=1, OP_NUM
       do p=1, coeffs_per_element
          outarr(:,op) = outarr(:,op) + avector(p)*fterm(:,op,p)
       end do
    end do
  end subroutine eval_ops_old


  !=====================================================
  ! define_fields
  !=====================================================
  subroutine define_fields(itri,psi_field, eta_field,eta79_l,psi79_l,f79_l,fterm,nu79_l,mu79_l)
    use m3dc1_data
    use mesh

    implicit none
    type(field_type) :: psi_field, eta_field 
    real, dimension(MAX_PTS, OP_NUM) :: eta79_l, psi79_l
    real, dimension(MAX_PTS, OP_NUM, nfield) :: f79_l
    real :: fterm(MAX_PTS, OP_NUM, coeffs_per_element)
    real, dimension(MAX_PTS, OP_NUM, dofs_per_element) :: nu79_l, mu79_l
!$acc routine 
!$acc routine (define_basis, precalculate_terms)

    integer, intent(in) :: itri

    integer :: i
    type(element_data) :: d

    !!$acc routine (local_to_global)
    !!$acc routine (eval_ops)
    !!$acc routine

    ! calculate the major radius, and useful powers
    call get_element_data(itri, d)
    call local_to_global(d, xi_79, zi_79, eta_79, x_79, phi_79, z_79)
    
    call precalculate_terms(xi_79,zi_79,eta_79,d%co,d%sn,fterm)
    call define_basis(itri,fterm,mu79_l,nu79_l)

    ! PSI
    ! ~~~
    call eval_ops(itri, psi_field, psi79_l,nu79_l)

    ! ETA
    ! ~~~
    call eval_ops(itri, eta_field, eta79_l,nu79_l)    

    ! define a lot of fields
    do i=1, nfield
       call eval_ops(itri, eta_field, f79_l(:,:,i),nu79_l)
    end do
  end subroutine define_fields
  
end module m3dc1_nint

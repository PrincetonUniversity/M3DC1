module matdef

  implicit none

  integer, parameter :: num_fields = 1
  integer, parameter :: psi_g = 1

contains

  !======================================================================
  ! Flux Equation
  !     call flux_lin(mu79(:,:,:),nu79(:,:,j),ss(:,j,:),dd(:,j,:))
  !x    call flux_lin(mu79(:,:,i),nu79(:,:,j),ss(i,j,:),dd(i,j,:))
  !======================================================================
  subroutine flux_lin(trial, lin, ssterm, ddterm, j, eta79_l, f79_l)
    use m3dc1_data
    use m3dc1_nint
    
    implicit none
!$acc routine
    
!   real, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
    real, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trial
    real, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: lin 
!   real, dimension(num_fields), intent(out) :: ssterm, ddterm
    real, dimension(dofs_per_element,dofs_per_element,num_fields), intent(out) :: ssterm, ddterm
    real, dimension(MAX_PTS, OP_NUM) :: eta79_l, psi79_l
    real, dimension(MAX_PTS, OP_NUM, nfield) :: f79_l
!   real :: temp
    real :: temp(dofs_per_element)
    real :: tempg(dofs_per_element)
    
    real :: thimp

    integer :: i, j
    
     thimp = 0.6

    !!$acc data present (f79, trial, lin, ssterm, ddterm) create (tempg)
!!$acc kernels
    ssterm = 0.
    ddterm = 0.
!!$acc end kernels
    
!!$acc kernels
    ! Time Derivatives
    ! ~~~~~~~~~~~~~~~~
!   temp = int2(trial(:,OP_1),lin(:,OP_1))
    tempg = int2(trial(:,:,OP_1),lin(:,OP_1,j))
!!$acc end kernels

!$acc loop independent
do i=1,dofs_per_element
    ssterm(i,psi_g) = ssterm(i,j,psi_g) + tempg(i)
    ddterm(i,psi_g) = ddterm(i,j,psi_g) + tempg(i)
enddo

    
!!$acc update device (ssterm, ddterm)
!!$acc serial 
    
    ! Resistive Terms
    ! ~~~~~~~~~~~~~~~
!   temp = -int3(trial(:,OP_DR),lin(:,OP_DR),eta79(:,OP_1)) &
!        -  int3(trial(:,OP_DZ),lin(:,OP_DZ),eta79(:,OP_1))     
    tempg = -int3(trial(:,:,OP_DR),lin(:,OP_DR,j),eta79_l(:,OP_1)) &
          -  int3(trial(:,:,OP_DZ),lin(:,OP_DZ,j),eta79_l(:,OP_1))     

!!$acc end serial 

!$acc  loop independent
do i=1,dofs_per_element
    ssterm(i,j,psi_g) = ssterm(i,j,psi_g) -     thimp *dt*tempg(i)
    ddterm(i,j,psi_g) = ddterm(i,j,psi_g) + (1.-thimp)*dt*tempg(i)
enddo


    ! Do a lot of metricterms calculations
!$acc  loop independent private(temp)
    do i=1, 100
       j = mod((i-1), nfield) + 1
!      temp = -int3(trial(:,OP_DR),lin(:,OP_DR),f79(:,OP_1,j)) &
!           -  int3(trial(:,OP_DZ),lin(:,OP_DZ),f79(:,OP_1,j)) &
!           +  int3(trial(:,OP_1),lin(:,OP_1),f79(:,OP_GS,j))     
       temp = -int3(trial(:,:,OP_DR),lin(:,OP_DR,j),f79_l(:,OP_1,j)) &
            -  int3(trial(:,:,OP_DZ),lin(:,OP_DZ,j),f79_l(:,OP_1,j)) &
            +  int3(trial(:,:,OP_1),lin(:,OP_1,j),f79_l(:,OP_GS,j))     
       ssterm(:,j,psi_g) = ssterm(:,j,psi_g) -     thimp *dt*temp(:)
       ddterm(:,j,psi_g) = ddterm(:,j,psi_g) + (1.-thimp)*dt*temp(:)
    end do
!!$acc end data
  end subroutine flux_lin


  !======================================================================
  ! matdefall
  ! ---------
  !
  ! Clears, populates, and finalizes all matrices for the implicit 
  ! time advance.  Does not insert boundary conditions, or finalize the
  ! s* matrices.
  !
  !======================================================================
  subroutine matdefall
    use mesh
    use nintegrate
    use m3dc1_data
    use m3dc1_nint

    implicit none
    
    integer :: itri, numelms
    
    logical :: is_edge(3)  ! is inode on boundary
    real :: n(2,3)
    integer :: iedge, idim(3), it
    real, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
    real, dimension(MAX_PTS, OP_NUM) :: eta79_l, psi79_l
    real, dimension(MAX_PTS, OP_NUM, dofs_per_element) :: mu79_l, nu79_l
    real, dimension(MAX_PTS, OP_NUM, nfield) :: f79_l
    real, dimension(MAX_PTS) :: xi_79_l, zi_79_l, eta_79_l, weight_79_l
    real, dimension(MAX_PTS,2) :: norm79_l
    real :: fterm(MAX_PTS, OP_NUM, coeffs_per_element)
    !$acc routine (matdef) vector
    !$acc routine (define_element_quadrature,define_fields,int0,int1) seq

    numelms = local_elements()
    
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP REDUCTION(+:volume) &
!$OMP REDUCTION(+:current)

    ! Loop over elements
    !$acc parallel loop gang present(psi_field, eta_field) &
    !$acc  private(ss,dd,eta79_l,psi79_l,f79_l,mu79_l,nu79_l,weight_79_l,fterm, &
    !$acc          xi_79_l,zi_79_l,norm79_l) &
    !$acc  reduction(+:volume,current)
    do itri=1,numelms
       ! Interior terms
       ! ~~~~~~~~~~~~~~
       it=itri 
       ! calculate the field values and derivatives at the sampling points
       call define_element_quadrature(it, 12, 5,weight_79_l, &
                                      xi_79_l, zi_79_l, eta_79_l, norm79_l)
       call define_fields(it,psi_field, eta_field,eta79_l,psi79_l,f79_l,fterm,&
                          nu79_l,mu79_l)
       
       ! add element's contribution to matrices
       call matdefphi(it,ss,dd,mu79_l,nu79_l,eta_79,f79_l)

       ! calculate some diagnostic quantities
       volume = volume + int0(weight_79_l)
       current = current + int1(psi79_l(:,OP_LP),weight_79_l)     
    end do

!$OMP END PARALLEL DO
  end subroutine matdefall
  
  
  !======================================================================
  ! matdefphi
  ! ---------
  !
  ! populates the matrices for implicit field advance
  !
  ! itri: index of finite element
  !======================================================================
  subroutine matdefphi(itri,ss,dd,mu79_l,nu79_l,eta79_l,f79_l)
    
    use element
    use nintegrate
    use m3dc1_nint

    implicit none
    
    integer, intent(in) :: itri
    
    integer :: i, j
    
    real, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
    real, dimension(MAX_PTS, OP_NUM, dofs_per_element) :: nu79_l, mu79_l
    real, dimension(MAX_PTS, OP_NUM, nfield) :: f79_l
    real, dimension(MAX_PTS, OP_NUM) :: eta79_l
 
    !$acc routine vector

    ss = 0.
    dd = 0.

    if(surface_int) return

    !$acc  loop vector  
       do j=1,dofs_per_element
!         call flux_lin(mu79(:,:,i),nu79(:,:,j),ss(i,j,:),dd(i,j,:))
!         call flux_lin(mu79(:,:,:),nu79(:,:,j),ss(:,j,:),dd(:,j,:))
          call flux_lin(mu79_l(:,:,:),nu79_l(:,:,:),ss(:,:,:),dd(:,:,:),j, &
                eta79_l,f79_l)
       end do       

  end subroutine matdefphi
end module matdef

module matdef

  implicit none

  integer, parameter :: num_fields = 1
  integer, parameter :: psi_g = 1

contains

  !======================================================================
  ! Flux Equation
  !======================================================================
  subroutine flux_lin(trial, lin, ssterm, ddterm)
    use m3dc1_data
    use m3dc1_nint
    
    implicit none
    
    real, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
    real, dimension(num_fields), intent(out) :: ssterm, ddterm
    real :: temp
    
    real :: thimp = 0.6

    integer :: i, j
    
    ssterm = 0.
    ddterm = 0.
    
    ! Time Derivatives
    ! ~~~~~~~~~~~~~~~~
    temp = int2(trial(:,OP_1),lin(:,OP_1))
    ssterm(psi_g) = ssterm(psi_g) + temp
    ddterm(psi_g) = ddterm(psi_g) + temp
    
    
    ! Resistive Terms
    ! ~~~~~~~~~~~~~~~
    temp = -int3(trial(:,OP_DR),lin(:,OP_DR),eta79(:,OP_1)) &
         -  int3(trial(:,OP_DZ),lin(:,OP_DZ),eta79(:,OP_1))     
    ssterm(psi_g) = ssterm(psi_g) -     thimp *dt*temp
    ddterm(psi_g) = ddterm(psi_g) + (1.-thimp)*dt*temp

    ! Do a lot of metricterms calculations
    do i=1, 100
       j = mod((i-1), nfield) + 1
       temp = -int3(trial(:,OP_DR),lin(:,OP_DR),f79(:,OP_1,j)) &
            -  int3(trial(:,OP_DZ),lin(:,OP_DZ),f79(:,OP_1,j)) &
	    +  int3(trial(:,OP_1),lin(:,OP_1),f79(:,OP_GS,j))     
       ssterm(psi_g) = ssterm(psi_g) -     thimp *dt*temp
       ddterm(psi_g) = ddterm(psi_g) + (1.-thimp)*dt*temp
    end do
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
    integer :: iedge, idim(3)
        
    numelms = local_elements()
    
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP REDUCTION(+:volume) &
!$OMP REDUCTION(+:current)

    ! Loop over elements
    do itri=1,numelms
       ! Interior terms
       ! ~~~~~~~~~~~~~~
       
       ! calculate the field values and derivatives at the sampling points
       call define_element_quadrature(itri, 12, 5)
       call define_fields(itri)
       
       ! add element's contribution to matrices
       call matdefphi(itri)

       ! calculate some diagnostic quantities
       volume = volume + int0()
       current = current + int1(psi79(:,OP_LP))     
       
       ! Surface terms
       ! ~~~~~~~~~~~~~
!       call boundary_edge(itri, is_edge, n, idim)
!       do iedge=1,3
!          ! if this edge is not on boundary, skip it
!          if(.not.is_edge(iedge)) cycle
!          
!          ! calculate the field values and derivatives at the sampling points
!          call define_boundary_quadrature(itri, iedge, 5, 5, n, idim)
!          call define_fields(itri)
!          
!          ! add edge's contribution to matrices
!          call matdefphi(itri)
!       end do
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
  subroutine matdefphi(itri)
    
    use element
    use nintegrate
    use m3dc1_nint

    implicit none
    
    integer, intent(in) :: itri
    
    integer :: i, j
    
    real, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
    
    ss = 0.
    dd = 0.

    if(surface_int) return

    do i=1,dofs_per_element
       do j=1,dofs_per_element
          call flux_lin(mu79(:,:,i),nu79(:,:,j),ss(i,j,:),dd(i,j,:))
       end do       
       ! Here is where the matrix elements would be inserted   
       ! call insert_block(bb1,itri,ieq(k),psi_i,ss(:,:,psi_g),MAT_ADD)
       ! call insert_block(bb0,itri,ieq(k),psi_i,dd(:,:,psi_g),MAT_ADD)
    end do
  end subroutine matdefphi
end module matdef

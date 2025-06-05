module biharmonic

  implicit none

  real, parameter :: a = 1
  real, parameter :: a2 = a**2
  real, dimension(MAX_PTS) :: r, r2, phi
  integer :: biharmonic_operator
  integer, parameter :: solution_order = 3

contains

  subroutine analytic_solution(n, x, z, soln)
    use math
    use element

    implicit none

    integer, intent(in) :: n
    real, intent(in) :: x, z
    real, intent(out), dimension(dofs_per_node) :: soln

    integer :: i
    real :: fac, s(dofs_per_node)
    
    soln = 0.

    do i=0, n
       fac = binomial(2*n+1, 2*i)*(-1)**i
       s = 0.
       s(1) = x**(2*(n-i)+1) * z**(2*i  ) &
            + z**(2*(n-i)+1) * x**(2*i  )
       if(2*(n-i)+1.gt.0) then
          s(2) = s(2) + x**(2*(n-i)  ) * z**(2*i  ) * (2*(n-i)+1)
          s(3) = s(3) + z**(2*(n-i)  ) * x**(2*i  ) * (2*(n-i)+1)
          if(i.gt.0) then
             s(5) = x**(2*(n-i)  ) * z**(2*i-1) * (2*(n-i)+1)*(2*i) &
                  + z**(2*(n-i)  ) * x**(2*i-1) * (2*(n-i)+1)*(2*i)
          endif
       endif
       if(2*(n-i).gt.0) then
          s(4) = s(4) + x**(2*(n-i)-1) * z**(2*i  ) * (2*(n-i)+1)*(2*(n-i))
          s(6) = s(6) + z**(2*(n-i)-1) * x**(2*i  ) * (2*(n-i)+1)*(2*(n-i))
       endif
       if(i.gt.0) then
          s(2) = s(2) + z**(2*(n-i)+1) * x**(2*i-1) * (2*i)
          s(3) = s(3) + x**(2*(n-i)+1) * z**(2*i-1) * (2*i)
       endif
       if(2*i-1.gt.0) then
          s(4) = s(4) + z**(2*(n-i)+1) * x**(2*i-2) * (2*i)*(2*i-1)
          s(6) = s(6) + x**(2*(n-i)+1) * z**(2*i-2) * (2*i)*(2*i-1)
       endif
       
       soln = soln + fac*s
    end do

  end subroutine analytic_solution

  subroutine biharmonic_init(bi)

    use basic
    use sparse
    use arrays
    use mesh_mod
    use m3dc1_nint
    use boundary_conditions

    implicit none
    
    integer, intent(in) :: bi
    integer :: itri, i, ii, j, numnodes, numelms, ier, inode_t
    integer :: is_edge(3)  ! is inode on boundary
    real :: n(2,3), sum, sum2, x, z, phi
    integer :: idim(3), izone, izonedim
    logical :: is_boundary
    real :: normal(2), curv(3)

    type(field_type) :: rhs, bcs
    real :: soln(MAX_PTS)

    integer, dimension(nodes_per_element) :: inode
    integer, parameter :: biharmonic_mat_index = 1
    type(matrix_type) :: biharmonic_mat
    integer, dimension(1, dofs_per_element) :: ind

    biharmonic_operator = bi

    call set_matrix_index(biharmonic_mat, biharmonic_mat_index)
    call create_mat(biharmonic_mat, 1, 1, icomplex, 1)

    numelms = local_elements()

    print *, 'populating matrix'
    ! populate the matrix
    do itri=1,numelms
       call get_element_nodes(itri,inode)

       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,0,1,0)

       call get_element_indices(1, itri, ind)
    
       do i=1, dofs_per_element
          do j=1, dofs_per_element
             
             select case(biharmonic_operator)
             case(0)
                temp79a = &
                     -(mu79(i,:,OP_DZ)*nu79(j,:,OP_DZ) &
                      +mu79(i,:,OP_DR)*nu79(j,:,OP_DR))
             case(1)
                temp79a = &
                     (mu79(i,:,OP_LP)*nu79(j,:,OP_LP) &
                     +mu79(i,:,OP_LP)*nu79(j,:,OP_LP))
             end select
             sum = int1(temp79a)
             call insert(biharmonic_mat, sum, ind(1,i), ind(1,j), MAT_ADD)
          enddo
       enddo
       
       if(isurface.eq.0) cycle
       
       ! add surface terms
       call boundary_edge(itri, is_edge, n, idim)
       
       do ii=1,edges_per_element
          if(is_edge(ii).eq.0) cycle
          
          call define_boundary_quadrature(itri, ii, 5, 5, n, idim)
          call define_fields(itri, 0, 1, 0)
          
          do i=1,dofs_per_element
             do j=1,dofs_per_element
                select case(biharmonic_operator)
                case(0)
                   temp79a = &
                        (norm79(:,1)*nu79(j,:,OP_DR) &
                        +norm79(:,2)*nu79(j,:,OP_DZ))
                   sum = int2(mu79(i,:,OP_1),temp79a)
                case(1)
                   temp79a = &
                        (norm79(:,1)*nu79(j,:,OP_LPR)*mu79(i,:,OP_1) &
                        +norm79(:,2)*nu79(j,:,OP_LPZ)*mu79(i,:,OP_1) &
                        -norm79(:,1)*mu79(i,:,OP_DR)*nu79(j,:,OP_LP) &
                        -norm79(:,2)*mu79(i,:,OP_DZ)*nu79(j,:,OP_LP))
                   sum = int1(temp79a)
                end select
                call insert(biharmonic_mat, sum, ind(1,i), ind(1,j), MAT_ADD)
             enddo
          enddo
       enddo
    end do

    print *, 'allocating arrays'
    call create_field(rhs)
    call create_field(bcs)

    print *, 'defining bcs'
    numnodes = owned_nodes()

    ! define boundary conditions
    do i=1, numnodes
       inode_t = nodes_owned(i)
       call boundary_node(inode_t,is_boundary,izone,izonedim,normal,curv, &
            x,phi,z)

       if(.not.is_boundary) cycle

!!$       bcs(ibegin  ) = x**3 * z**3
!!$       bcs(ibegin+1) = 3.*x**2 * z**3
!!$       bcs(ibegin+2) = 3.*x**3 * z**2
!!$       bcs(ibegin+3) = 6.*x    * z**3
!!$       bcs(ibegin+4) = 9.*x**2 * z**2
!!$       bcs(ibegin+5) = 6.*x**3 * z
!!$
!!$       bcs(ibegin  ) =    x**3 - 3.*x**2*z - 3.*x*z**2 +    z**3
!!$       bcs(ibegin+1) = 3.*x**2 - 6.*x   *z - 3.  *z**2
!!$       bcs(ibegin+2) =         - 3.*x**2   - 6.*x*z    + 3.*z**2
!!$       bcs(ibegin+3) = 6.*x    - 6.     *z
!!$       bcs(ibegin+4) =         - 6.*x      - 6.  *z
!!$       bcs(ibegin+5) =                     - 6.*x      + 6.*z
!!$
!!$
       call analytic_solution(solution_order,x,z,soln(1:dofs_per_node))
       call set_node_data(bcs, inode_t, soln(1:dofs_per_node))
    end do
   
    call boundary_biharmonic(rhs, bcs, biharmonic_mat)

    ! solve equation
    print *, 'solving'
    call finalize(biharmonic_mat)

    call newsolve(biharmonic_mat, rhs%vec, ier)
   
    ! calculate error
    sum = 0.
    sum2 = 0.
    do itri=1, numelms
       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,0,0,0)

       call eval_ops(itri, rhs, ps079)

!!$       call biharmonic_solution(x_79,z_79,1.,soln)
!!$       ps179(:,OP_1) = soln
!!$       ps179(:,OP_1) = x_79**3 - 3.*x_79**2*z_79 - 3.*x_79*z_79**2 + z_79**3
       do i=1, npoints
          call analytic_solution(solution_order,x_79(i),z_79(i), &
               soln(1:dofs_per_node))
          ps179(i,OP_1) = soln(1)
       end do

       sum = sum + int1((ps079(:,OP_1) - ps179(:,OP_1))**2)
       sum2 = sum2 + int1(ps179(:,OP_1)**2)
    end do

    print *, 'Error in biharmonic solution = ', sum/sum2, sum2

    ! copy solution to psi0
    print *, 'copying solution bcs'
    psi_field(0) = rhs

    ! free temporary arrays
    print *, 'deallocating arrays'
    call destroy_field(rhs)
    call destroy_field(bcs)
    call destroy_mat(biharmonic_mat)
   
  end subroutine biharmonic_init

  function biharmonic_integrand(eta)
    implicit none

    real, dimension(MAX_PTS) :: biharmonic_integrand
    real, intent(in) :: eta

    real, dimension(MAX_PTS) :: denom, co
    real :: f, g

    f = (a2*cos(eta)*sin(eta))**3
    g = 6.*f/a
    co = cos(eta-phi)
    denom = 1./(r2 + a2 - 2.*a*r*co)

    select case(biharmonic_operator)
    case(0)
       biharmonic_integrand = f*denom
    case(1)
       biharmonic_integrand = (f*(a - r*co)*denom - 0.5*g)*denom
    end select

  end function biharmonic_integrand
 
  subroutine biharmonic_solution(x, z, a, w)
    use math
    
    implicit none

    real, intent(in), dimension(MAX_PTS) :: x, z
    real, intent(in) :: a
    real, intent(out), dimension(MAX_PTS) :: w
    real :: eta, deta
    integer :: i, n

    r2 = x**2 + z**2
    r = sqrt(r2)
    phi = atan2(z,x)

    ! integrate using trapezoid rule
    n = 100000

    w = biharmonic_integrand(0.) + biharmonic_integrand(2.*pi)

    deta = twopi/n
    eta = deta
    do i=1, n-1
       w = w + 2.*biharmonic_integrand(eta)
       eta = eta + deta
    end do
    w = w*twopi/(2.*n)

    select case(biharmonic_operator)
    case(0)
       w = w*(a2 - r2)/(twopi)
    case(1)
       w = w*(a2 - r2)**2/(twopi*a)
    end select
   
  end subroutine biharmonic_solution

  subroutine boundary_biharmonic(rhs, bvec, mat)
    use basic
    use boundary_conditions
    use field
    use matrix_mod

    implicit none
  
    type(field_type), intent(inout) :: rhs
    type(field_type), intent(in) :: bvec
    type(matrix_type), intent(inout), optional :: mat
    
    integer :: i, inode_t, izone, izonedim, index
    integer :: numnodes
    real :: normal(2), curv(3)
    real :: x, phi, z
    logical :: is_boundary
    vectype, dimension(dofs_per_node) :: temp

    if(iper.eq.1 .and. jper.eq.1) return
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_biharmonic called"
    
    numnodes = owned_nodes()
    do i=1, numnodes
       inode_t = nodes_owned(i)
       call boundary_node(inode_t,is_boundary,izone,izonedim,normal,curv, &
            x,phi,z)
       if(.not.is_boundary) cycle
       
       index = node_index(rhs, inode_t)
       
       call get_node_data(bvec, inode_t, temp)
       call set_dirichlet_bc(index,rhs%vec,temp,normal,curv,izonedim,mat)

       if(biharmonic_operator.eq.1) then
          call set_normal_bc(index,rhs%vec,temp,normal,curv,izonedim,mat)
       endif
    end do
    
  end subroutine boundary_biharmonic

end module biharmonic

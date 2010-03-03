module biharmonic

  implicit none

  integer, parameter :: npts_biharmonic = 25
  real, parameter :: a = 1
  real, parameter :: a2 = a**2
  real, dimension(npts_biharmonic) :: r, r2, phi
  integer :: biharmonic_operator
  integer, parameter :: solution_order = 3

contains

  integer function factorial(n)
    implicit none

    integer, intent(in) :: n

    integer :: i

    factorial = 1.
    do i=2,n
       factorial = factorial*i
    end do
  end function factorial

  real function binomial(n, m)
    implicit none

    integer, intent(in) :: n, m

    binomial = factorial(n)/(factorial(m)*factorial(n-m))
  end function binomial

  subroutine analytic_solution(n, x, z, soln)
    implicit none

    integer, intent(in) :: n
    real, intent(in) :: x, z
    real, intent(out), dimension(6) :: soln

    integer :: i
    real :: fac, s(6)
    
    soln(1:6) = 0.

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
    use t_data
    use nintegrate_mod

    implicit none
    
    integer, intent(in) :: bi
    integer :: itri, i, ii, i1, j, j1, numnodes, numelms, ier
    integer :: ibegin, iendplusone
    logical :: is_edge(3)  ! is inode on boundary
    real :: n(2,3), sum, sum2, x, z
    integer :: idim(3), izone, izonedim
    logical :: is_boundary
    real :: normal(2), curv
    vectype, dimension(20) :: avec
  
    vectype, allocatable :: rhs(:), bcs(:)
    real :: soln(MAX_PTS)

    integer, parameter :: biharmonic_sm = 1

    biharmonic_operator = bi

    print *, 'zeroing matrix'
    call zerosuperlumatrix(biharmonic_sm, icomplex, 1)

    call numfac(numelms)

    print *, 'populating matrix'
    ! populate the matrix
    do itri=1,numelms

       call define_triangle_quadrature(itri,25)
       call define_fields(itri,0,1,0)
    
       do i=1,18
          i1 = isval1(itri,i)
          do j=1,18
             j1 = isval1(itri,j)
             select case(biharmonic_operator)
             case(0)
                temp79a = &
                     -(g79(:,OP_DZ,i)*g79(:,OP_DZ,j) &
                      +g79(:,OP_DR,i)*g79(:,OP_DR,j))
             case(1)
                temp79a = &
                     (g79(:,OP_LP,i)*g79(:,OP_LP,j) &
                     +g79(:,OP_LP,i)*g79(:,OP_LP,j))
             end select
             sum = int1(temp79a)
             call insertval(biharmonic_sm, sum, 0, i1,j1,1)
          enddo
       enddo
       
       if(isurface.eq.0) cycle
       
       ! add surface terms
       call boundary_edge(itri, is_edge, n, idim)
       
       do ii=1,3
          if(.not.is_edge(ii)) cycle
          
          call define_edge_quadrature(itri, ii, 5, n, idim)
          call define_fields(itri, 0, 1, 0)
          
          do j=1,18
             j1 = isval1(itri,j)
             do i=1,18
                i1 = isval1(itri,i)
                select case(biharmonic_operator)
                case(0)
                   temp79a = &
                        (norm79(:,1)*g79(:,OP_DR,j) &
                        +norm79(:,2)*g79(:,OP_DZ,j))
                   sum = int2(g79(:,OP_1,i),temp79a)
                case(1)
                   temp79a = &
                        (norm79(:,1)*g79(:,OP_LPR,j)*g79(:,OP_1,i) &
                        +norm79(:,2)*g79(:,OP_LPZ,j)*g79(:,OP_1,i) &
                        -norm79(:,1)*g79(:,OP_DR,i)*g79(:,OP_LP,j) &
                        -norm79(:,2)*g79(:,OP_DZ,i)*g79(:,OP_LP,j))
                   sum = int1(temp79a)
                end select
                call insertval(biharmonic_sm, sum, 0, i1,j1,1)
             enddo
          enddo
       end do
    enddo

    print *, 'allocating arrays'
    call createvec(rhs, 1)
    call createvec(bcs, 1)

    print *, 'defining bcs'
    call numnod(numnodes)

    ! define boundary conditions
    do i=1, numnodes
       call entdofs(1, i, 0, ibegin, iendplusone)

       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)

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
       call analytic_solution(solution_order,x,z,soln(1:6))
       bcs(ibegin:ibegin+5) = soln(1:6)

    end do
   
    call boundary_biharmonic(biharmonic_sm, rhs, bcs)

    ! solve equation
    print *, 'solving'
    call finalizematrix(biharmonic_sm)

    call solve(biharmonic_sm, rhs, ier)
   
    ! calculate error
    sum = 0.
    sum2 = 0.
    do itri=1, numelms
       call define_triangle_quadrature(itri,npts_biharmonic)
       call define_fields(itri,0,0,0)

       call calcavector(itri, rhs, 1, 1, avec)
       call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ps079)

!!$       call biharmonic_solution(x_79,z_79,1.,soln)
!!$       ps179(:,OP_1) = soln
!!$       ps179(:,OP_1) = x_79**3 - 3.*x_79**2*z_79 - 3.*x_79*z_79**2 + z_79**3
       do i=1, npts_biharmonic
          call analytic_solution(solution_order,x_79(i),z_79(i),soln(1:6))
          ps179(i,OP_1) = soln(1)
       end do

       sum = sum + int1((ps079(:,OP_1) - ps179(:,OP_1))**2)
       sum2 = sum2 + int1(ps179(:,OP_1)**2)
    end do

    print *, 'Error in biharmonic solution = ', sum/sum2, sum2

    ! copy solution to psi0
    print *, 'copying solution bcs'
    call copyvec(rhs, 1, 1, field0, psi_g, num_fields)

    ! free temporary arrays
    print *, 'deallocating arrays'
    call deletevec(rhs)
    call deletevec(bcs)
   
  end subroutine biharmonic_init

  function biharmonic_integrand(eta)
    implicit none

    real, dimension(npts_biharmonic) :: biharmonic_integrand
    real, intent(in) :: eta

    real, dimension(npts_biharmonic) :: denom, co
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

    implicit none

    real, intent(in), dimension(npts_biharmonic) :: x, z
    real, intent(in) :: a
    real, intent(out), dimension(npts_biharmonic) :: w
    real :: eta, deta
    real, parameter :: pi = 3.14159265358979323846
    integer :: i, n

    r2 = x**2 + z**2
    r = sqrt(r2)
    phi = atan2(z,x)

    ! integrate using trapezoid rule
    n = 100000

    w = biharmonic_integrand(0.) + biharmonic_integrand(2.*pi)

    deta = 2.*pi/n
    eta = deta
    do i=1, n-1
       w = w + 2.*biharmonic_integrand(eta)
       eta = eta + deta
    end do
    w = w*2.*pi/(2.*n)

    select case(biharmonic_operator)
    case(0)
       w = w*(a2 - r2)/(2.*pi)
    case(1)
       w = w*(a2 - r2)**2/(2.*pi*a)
    end select
   
  end subroutine biharmonic_solution

  subroutine boundary_biharmonic(imatrix, rhs, bvec)
    use basic

    implicit none
  
    integer, intent(in) :: imatrix
    vectype, intent(inout), dimension(*) :: rhs
    vectype, intent(in), dimension(*) :: bvec
    
    integer :: i, izone, izonedim
    integer :: ibegin, iendplusone, numnodes
    real :: normal(2), curv
    real :: x, z
    logical :: is_boundary
    vectype, dimension(6) :: temp

    if(iper.eq.1 .and. jper.eq.1) return
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_biharmonic called"
    
    call numnod(numnodes)
    do i=1, numnodes
       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
       if(.not.is_boundary) cycle
       
       call entdofs(1, i, 0, ibegin, iendplusone)
       call rotate_matrix(imatrix, ibegin, normal, curv, rhs, icurv)
       
       temp = bvec(ibegin:ibegin+5)
       call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim)

       if(biharmonic_operator.eq.1) then
          call set_normal_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim)
       endif
    end do
    
  end subroutine boundary_biharmonic

end module biharmonic

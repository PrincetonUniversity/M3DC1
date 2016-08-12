module fit_magnetics
  implicit none

  integer, private :: nk
  real, private, allocatable :: x(:)
  real, private :: dk

contains

  subroutine solve_bessel_fields_2D(n, r, z, rn, zn, b,ierr)
    use math

    implicit none

    integer, intent(in) :: n                 ! number of measurements
    real, intent(in), dimension(n) :: r, z   ! position of measurements
    real, intent(in), dimension(n) :: rn, zn ! direction of measurements
    real, intent(in), dimension(n) :: b      ! value of measurements
    integer, intent(out) :: ierr
    real, allocatable :: matrix(:,:), work(:), s(:)
    integer :: i,j,m, lwork, rank
    real :: k
    real :: rcond = -1

    ! define the matrix
    nk = 3   ! number of modes to include
    m = 4*nk
    dk = pi  ! spacing of modes

    if(m.gt.n) then
       print *, 'Error: underdetrmined matrix'
       ierr = 1
       return
    end if

    allocate(matrix(n,m))

    do i=1, n
       do j=1, nk
          k = j*dk
          matrix(i,j     ) =  k*exp( k*z(i))* &
               (zn(i)*Bessel_J(0,k*r(i)) - rn(i)*Bessel_J(1,k*r(i)))
          matrix(i,j+  nk) =  k*exp( k*z(i))* &
               (zn(i)*Bessel_Y(0,k*r(i)) - rn(i)*Bessel_Y(1,k*r(i)))
          matrix(i,j+2*nk) = -k*exp(-k*z(i))* &
               (zn(i)*Bessel_J(0,k*r(i)) + rn(i)*Bessel_J(1,k*r(i)))
          matrix(i,j+3*nk) = -k*exp(-k*z(i))* &
               (zn(i)*Bessel_Y(0,k*r(i)) + rn(i)*Bessel_Y(1,k*r(i)))
       end do
    end do

    ! get least-squares solution
    lwork = 3*m+1*(m+n)
    allocate(s(m),work(lwork))
    call dgelss(n,m,1,matrix,n,b,n,s,rcond,rank,work,lwork,ierr)
    print *, ' dgelss returned with info = ', ierr
!    print *, ' Singular values = ', s
!    print *, ' rank = ', rank
    
    deallocate(work,s,matrix)

    allocate(x(m))
    x = 0.
    x = b(1:m)
  end subroutine solve_bessel_fields_2D


  subroutine evaluate_bessel_fields_2D(n, r, z, psi)
    use math

    implicit none

    integer, intent(in) :: n
    real, intent(in), dimension(*) :: r, z
    vectype, intent(out), dimension(*) :: psi

    integer :: i, j
    real :: k, bj, by

    psi(1:n) = 0.
    do j=1, n
       do i=1, nk
          k = i*dk
          bj = Bessel_J(1,k*r(j))
          by = Bessel_Y(1,k*r(j))
          psi(j) = psi(j) &
               + r(j)*exp( k*z(j))*(x(i     )*bj + x(i+  nk)*by) &
               - r(j)*exp(-k*z(j))*(x(i+2*nk)*bj + x(i+3*nk)*by)
       end do
    end do    
  end subroutine evaluate_bessel_fields_2D

  subroutine clear_bessel_fields()
    implicit none

    if(allocated(x)) deallocate(x)
  end subroutine clear_bessel_fields


  subroutine solve_multipole_fields_2D(n, r, z, rn, zn, b,ierr)
    use math
    use coils

    implicit none

    integer, intent(in) :: n                 ! number of measurements
    real, intent(in), dimension(n) :: r, z   ! position of measurements
    real, intent(in), dimension(n) :: rn, zn ! direction of measurements
    real, intent(in), dimension(n) :: b      ! value of measurements
    integer, intent(out) :: ierr
    real, allocatable :: matrix(:,:), work(:), s(:)
    integer :: i,j,m, lwork, rank
    real :: rcond = -1
    real, dimension(n,1,6) :: g
    real, dimension(1) :: xi, zi

    ! define the matrix
    nk = 4   ! number of modes to include
    m = nk   ! columns of matrix

    if(m.gt.n) then
       print *, 'Error: underdetrmined matrix'
       ierr = 1
       return
    end if

    if(nk.gt.8) then
       print *, 'Error: can only use up to 8 multipole fields'
       ierr = 1
       return
    end if

    allocate(matrix(n,m))

    do j=1, nk
       xi = 101.+j
       zi = 10.
       call gvect(r,z,n,xi,zi,1,g,1,ierr)

       matrix(:,j) = -rn*g(:,1,3)/r + zn*g(:,1,2)/r
    end do

    ! get least-squares solution
    lwork = 3*m+1*(m+n)
    allocate(s(m),work(lwork))
    call dgelss(n,m,1,matrix,n,b,n,s,rcond,rank,work,lwork,ierr)
    print *, ' dgelss returned with info = ', ierr
!    print *, ' Singular values = ', s
!    print *, ' rank = ', rank
    
    deallocate(work,s,matrix)

    allocate(x(m))
    x = 0.
    x = b(1:m)
  end subroutine solve_multipole_fields_2D

  subroutine evaluate_multipole_fields_2D(n, r, z, psi)
    use math
    use coils

    implicit none

    integer, intent(in) :: n
    real, intent(in), dimension(n) :: r, z
    vectype, intent(out), dimension(n) :: psi

    integer :: i, j, ierr
    real, dimension(n,1,6) :: g
    real, dimension(n) :: xi, zi

    psi = 0.
    do j=1, nk

       xi = 101.+j
       zi = 10.

       call gvect(r,z,n,xi,zi,1,g,1,ierr)

       psi = psi + x(j)*g(:,1,1)
    end do    
  end subroutine evaluate_multipole_fields_2D

end module fit_magnetics

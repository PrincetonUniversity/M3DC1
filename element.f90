module element

!!$  integer, parameter :: VAL       = 1
!!$  integer, parameter :: VAL_DX    = 2
!!$  integer, parameter :: VAL_DZ    = 3
!!$  integer, parameter :: VAL_DXX   = 4
!!$  integer, parameter :: VAL_DXZ   = 5
!!$  integer, parameter :: VAL_DZZ   = 6
#ifdef USE3D
!!$  integer, parameter :: VAL_DY    = 8
!!$  integer, parameter :: VAL_DXY   = 9
!!$  integer, parameter :: VAL_DYZ   = 10
!!$  integer, parameter :: VAL_DXXY  = 11
!!$  integer, parameter :: VAL_DXYZ  = 12
!!$  integer, parameter :: VAL_DYZZ  = 13
!!$  integer, parameter :: VAL_DYY   = 14
!!$  integer, parameter :: VAL_DXYY  = 15
!!$  integer, parameter :: VAL_DYYZ  = 16
!!$  integer, parameter :: VAL_DXXYY = 17
!!$  integer, parameter :: VAL_DXYYZ = 18
!!$  integer, parameter :: VAL_DYYZZ = 19
  integer, parameter :: dofs_per_node = 19
#else
  integer, parameter :: dofs_per_node = 6
#endif

  integer, parameter :: nodes_per_edge = 2
  integer, parameter :: nodes_per_element = 3
  integer, parameter :: edges_per_element = 3
  integer, parameter :: dofs_per_element = nodes_per_element*dofs_per_node

  type element_data
     real :: x, z, a, b, c, co, sn
  end type element_data

contains

  
  !=======================================================
  ! global_to_local
  ! ~~~~~~~~~~~~~~~
  ! transforms from global coordinates 
  ! to local (element) coordinates
  !=======================================================
  elemental subroutine global_to_local(d, x, z, xi, eta)
    implicit none

    type(element_data), intent(in) :: d
    real, intent(in) :: x, z
    real, intent(out) :: xi, eta

    xi  =  (x-d%x)*d%co + (z-d%z)*d%sn - d%b
    eta = -(x-d%x)*d%sn + (z-d%z)*d%co
  end subroutine global_to_local

  !=======================================================
  ! local_to_global
  ! ~~~~~~~~~~~~~~~
  ! transforms from local (element) coordinates
  ! to global coordinates
  !=======================================================   
  elemental subroutine local_to_global(d, xi, eta, x, z)
    implicit none

    type(element_data), intent(in) :: d
    real, intent(in) :: xi, eta
    real, intent(out) :: x, z

    x = d%x + (d%b+xi)*d%co - eta*d%sn
    z = d%z + (d%b+xi)*d%sn + eta*d%co
  end subroutine local_to_global

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
    real :: t(20,20), wkspce(9400)
    integer :: ipiv(20)
    
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
  
end module element

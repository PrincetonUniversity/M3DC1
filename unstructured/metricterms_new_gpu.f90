module metricterms_new
  use m3dc1_nint

implicit none

  type muarray
     sequence
     integer :: len
     vectype, dimension(MAX_PTS, OP_NUM) :: value79
     integer, dimension(OP_NUM) :: op
  end type muarray

  type nuarray
     sequence
     integer :: len
     vectype, dimension(MAX_PTS, OP_NUM) :: value79
     integer, dimension(OP_NUM) :: op
  end type nuarray

  type prodarray
     sequence
     integer :: len
     vectype, dimension(MAX_PTS, 51) :: value79
     integer, dimension(51) :: op1
     integer, dimension(51) :: op2
  end type prodarray

  interface operator(+)
     module procedure prod_add
     module procedure mu_add
  end interface

  interface operator(-)
     module procedure prod_sub
  end interface

  interface operator(*)
     module procedure mu_mult_array
     module procedure prod_mult_array
     module procedure real_mult_mu
     module procedure real_mult_prod
  end interface

  interface mu
     module procedure mu_new_array
  end interface

  interface prod
     module procedure prod_mu_mu
     module procedure prod_new_array
  end interface

contains

function mu_new_array(x2, op)

  implicit none

  type(muarray) :: mu_new_array
  vectype, dimension(MAX_PTS), intent(in) :: x2
  integer, intent(in) :: op
  integer :: i, ind

  mu_new_array%len=1
  mu_new_array%value79(:,mu_new_array%len)=x2
  mu_new_array%op(mu_new_array%len)=op

end function mu_new_array

function mu_add(v1, v2)

  implicit none

  type(muarray) :: mu_add
  type(muarray), intent(in) :: v1, v2
  integer, dimension(OP_NUM) :: op_index
  integer :: i, ind

  mu_add=v1
  op_index=0
  do i=1,v1%len
     op_index(v1%op(i))=i
  enddo
  do i=1,v2%len
     ind=op_index(v2%op(i))
     if (ind>0) then
        mu_add%value79(:,ind)=v1%value79(:,ind)+v2%value79(:,i)
     else
        mu_add%len=mu_add%len+1
        mu_add%value79(:,mu_add%len)=v2%value79(:,i)
        mu_add%op(mu_add%len)=v2%op(i)
     endif
  enddo
  !mu_add%len=v1%len+v2%len
  !mu_add%op(1:v1%len)=v1%op(1:v1%len)
  !mu_add%value79(:,1:v1%len)=v1%value79(:,1:v1%len)
  !mu_add%op(v1%len+1:v1%len+v2%len)=v2%op(1:v2%len)
  !mu_add%value79(:,v1%len+1:v1%len+v2%len)=v2%value79(:,1:v2%len)

end function mu_add

function prod_new_array(x2, op1, op2)

  implicit none

  type(prodarray) :: prod_new_array
  vectype, dimension(MAX_PTS), intent(in) :: x2
  integer, intent(in) :: op1, op2

  prod_new_array%len=1
  prod_new_array%value79(:,prod_new_array%len)=x2
  prod_new_array%op1(prod_new_array%len)=op1
  prod_new_array%op2(prod_new_array%len)=op2

end function prod_new_array

subroutine prod_add_array(p1, x2, op1, op2)

  implicit none

  type(prodarray), intent(inout) :: p1
  vectype, dimension(MAX_PTS), intent(in) :: x2
  integer, intent(in) :: op1, op2
  integer :: i, ind

  do i=1,p1%len
     if ((p1%op1(i)==op1).and.(p1%op2(i)==op2)) then
        p1%value79(:,i)=p1%value79(:,i)+x2
        exit
     endif
  enddo
  p1%len=p1%len+1
  p1%value79(:,p1%len)=x2
  p1%op1(p1%len)=op1
  p1%op2(p1%len)=op2

end subroutine prod_add_array

function prod_add(p1, p2)

  implicit none

  type(prodarray) :: prod_add
  type(prodarray), intent(in) :: p1, p2
  integer, dimension(OP_NUM,OP_NUM) :: op_index
  integer :: i, ind

  prod_add=p1
  op_index=0
  do i=1,p1%len
     op_index(p1%op1(i),p1%op2(i))=i
  enddo
  do i=1,p2%len
     ind=op_index(p2%op1(i),p2%op2(i))
     if (ind>0) then
        prod_add%value79(:,ind)=p1%value79(:,ind)+p2%value79(:,i)
     else
        prod_add%len=prod_add%len+1
        prod_add%value79(:,prod_add%len)=p2%value79(:,i)
        prod_add%op1(prod_add%len)=p2%op1(i)
        prod_add%op2(prod_add%len)=p2%op2(i)
     endif
  enddo
  !prod_add%len=p1%len+p2%len
  !prod_add%op1(1:p1%len)=p1%op1(1:p1%len)
  !prod_add%op2(1:p1%len)=p1%op2(1:p1%len)
  !prod_add%value79(:,1:p1%len)=p1%value79(:,1:p1%len)
  !prod_add%op1(p1%len+1:p1%len+p2%len)=p2%op1(1:p2%len)
  !prod_add%op2(p1%len+1:p1%len+p2%len)=p2%op2(1:p2%len)
  !prod_add%value79(:,p1%len+1:p1%len+p2%len)=p2%value79(:,1:p2%len)

end function prod_add

function prod_sub(p1, p2)

  implicit none

  type(prodarray) :: prod_sub
  type(prodarray), intent(in) :: p1, p2
  integer, dimension(OP_NUM,OP_NUM) :: op_index
  integer :: i, ind

  prod_sub=p1
  op_index=0
  do i=1,p1%len
     op_index(p1%op1(i),p1%op2(i))=i
  enddo
  do i=1,p2%len
     ind=op_index(p2%op1(i),p2%op2(i))
     if (ind>0) then
        prod_sub%value79(:,ind)=p1%value79(:,ind)-p2%value79(:,i)
     else
        prod_sub%len=prod_sub%len+1
        prod_sub%value79(:,prod_sub%len)=-p2%value79(:,i)
        prod_sub%op1(prod_sub%len)=p2%op1(i)
        prod_sub%op2(prod_sub%len)=p2%op2(i)
     endif
  enddo
end function prod_sub

function mu_mult_array(v1, x2)

  implicit none

  type(muarray) :: mu_mult_array
  type(muarray), intent(in) :: v1
  vectype, dimension(MAX_PTS), intent(in) :: x2
  integer :: i

  mu_mult_array%len=v1%len
  mu_mult_array%op=v1%op
  do i=1,v1%len
     mu_mult_array%value79(:,i)=v1%value79(:,i)*x2
  enddo

end function mu_mult_array

function real_mult_mu(x2, v1)

  implicit none

  type(muarray) :: real_mult_mu
  type(muarray), intent(in) :: v1
  real, intent(in) :: x2
  integer :: i

  real_mult_mu%len=v1%len
  real_mult_mu%op=v1%op
  do i=1,v1%len
     real_mult_mu%value79(:,i)=v1%value79(:,i)*x2
  enddo

end function real_mult_mu

function prod_mult_array(p1, x2)

  implicit none

  type(prodarray) :: prod_mult_array
  type(prodarray), intent(in) :: p1
  vectype, dimension(MAX_PTS), intent(in) :: x2
  integer :: i

  prod_mult_array%len=p1%len
  prod_mult_array%op1=p1%op1
  prod_mult_array%op2=p1%op2
  do i=1,p1%len
     prod_mult_array%value79(:,i)=p1%value79(:,i)*x2
  enddo

end function prod_mult_array

function real_mult_prod(x2, p1)

  implicit none

  type(prodarray) :: real_mult_prod
  type(prodarray), intent(in) :: p1
  real, intent(in) :: x2
  integer :: i

  real_mult_prod%len=p1%len
  real_mult_prod%op1=p1%op1
  real_mult_prod%op2=p1%op2
  do i=1,p1%len
     real_mult_prod%value79(:,i)=p1%value79(:,i)*x2
  enddo

end function real_mult_prod

function prod_mu_mu(v1, v2)

  implicit none

  type(prodarray) :: prod_mu_mu
  type(muarray), intent(in) :: v1, v2
  integer :: i, j

  prod_mu_mu%len=0
  do i=1,v1%len
     do j=1,v2%len
        prod_mu_mu%len=prod_mu_mu%len+1
        prod_mu_mu%op1(prod_mu_mu%len)=v1%op(i)
        prod_mu_mu%op2(prod_mu_mu%len)=v2%op(j)
        prod_mu_mu%value79(:,prod_mu_mu%len)=v1%value79(:,i)*v2%value79(:,j)
     enddo
  enddo

end function prod_mu_mu

!============================================================================
! V1 TERMS
!============================================================================
  

! V1umu 
! =====
! function v1umu(e,f,g,h)

!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1umu
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
!   vectype, dimension(dofs_per_element) :: temp

!      temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
!      if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

!      if(surface_int) then
!         if(inoslip_pol.eq.1) then
!            temp = 0.
!         else
!            temp = intx5(e(:,:,OP_DZ),r2_79,g(:,OP_1),norm79(:,2),temp79b) &
!                 - intx5(e(:,:,OP_DR),r2_79,g(:,OP_1),norm79(:,1),temp79b) &
!                 + 2.* &
!                 (intx5(e(:,:,OP_DZ),r2_79,g(:,OP_1),norm79(:,1),f(:,OP_DRZ)) &
!                 +intx5(e(:,:,OP_DR),r2_79,g(:,OP_1),norm79(:,2),f(:,OP_DRZ)))
           
!            if(itor.eq.1) then
!               temp79a = h(:,OP_1) - g(:,OP_1)
!               temp = temp &
!                    + 2.*intx5(e(:,:,OP_DZ),r_79,g(:,OP_1),norm79(:,1),f(:,OP_DZ))&
!                    + 4.* &
!                    (intx5(e(:,:,OP_DZ),r_79,norm79(:,1),f(:,OP_DZ),temp79a) &
!                    -intx5(e(:,:,OP_DR),r_79,norm79(:,2),f(:,OP_DZ),temp79a)) &
!                    + 4.* &
!                    (intx5(e(:,:,OP_1),r_79,norm79(:,1),f(:,OP_DZZ),h(:,OP_1)) &
!                    -intx5(e(:,:,OP_1),r_79,norm79(:,2),f(:,OP_DRZ),h(:,OP_1)))
!            endif

! #if defined(USE3D) || defined(USECOMPLEX)
!            temp = temp &
!                 - intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DRPP),g(:,OP_1)) &
!                 - intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DZPP),g(:,OP_1))
! #endif
!         endif
!      else 
!         temp79d = 2.*f(:,OP_DRZ)
!         if(itor.eq.1) then
!            temp79d = temp79d + ri_79*f(:,OP_DZ)
!         endif

!         temp = &
!              - intx4(e(:,:,OP_DZZ),r2_79,temp79b,g(:,OP_1)) &
!              + intx4(e(:,:,OP_DRR),r2_79,temp79b,g(:,OP_1)) &
!              - 2.*intx4(e(:,:,OP_DRZ),r2_79,temp79d,g(:,OP_1))
        
!         if(itor.eq.1) then
!            temp = temp &
!                 + intx4(e(:,:,OP_DR),r_79,temp79b,g(:,OP_1)) &
!                 - intx4(e(:,:,OP_DZ),r_79,temp79d,g(:,OP_1)) &
!                 + 5.*intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
!                 - 8.*intx3(e(:,:,OP_DZ),f(:,OP_DZ),h(:,OP_1))
!         endif
     
! #if defined(USE3D) || defined(USECOMPLEX)
!         temp = temp &
! #ifdef USEST
!              - intx3(e(:,:,OP_DZP),f(:,OP_DZP),g(:,OP_1)) &
!              - intx3(e(:,:,OP_DRP),f(:,OP_DRP),g(:,OP_1)) &
!              - intx3(e(:,:,OP_DZ),f(:,OP_DZP),g(:,OP_DP)) &
!              - intx3(e(:,:,OP_DR),f(:,OP_DRP),g(:,OP_DP))
! #else
!              + intx3(e(:,:,OP_DZ),f(:,OP_DZPP),g(:,OP_1)) &
!              + intx3(e(:,:,OP_DR),f(:,OP_DRPP),g(:,OP_1))
! #endif
! #endif
!      end if

!   v1umu = temp
! end function v1umu

function v1umu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1umu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  type(muarray) :: tempb, tempd

     tempb = mu(g(:,OP_1),OP_DZZ) + mu(-g(:,OP_1),OP_DRR)
     if(itor.eq.1) tempb = tempb + mu(-g(:,OP_1)*ri_79,OP_DR)

     if(surface_int) then
           temp%len = 0
     else
        tempd = mu(2.*g(:,OP_1),OP_DRZ)
        if(itor.eq.1) then
           tempd = tempd + mu(g(:,OP_1)*ri_79,OP_DZ)
        endif

        temp = prod(mu(-r2_79,OP_DZZ),tempb) &
             + prod(mu( r2_79,OP_DRR),tempb) &
             + prod(mu(-2.*r2_79,OP_DRZ),tempd)
        
        if(itor.eq.1) then
           temp = temp &
                + prod(mu( r_79,OP_DR),tempb) &
                + prod(mu(-r_79,OP_DZ),tempd) &
                + prod( 5.*g(:,OP_1),OP_DZ,OP_DZ) &
                + prod(-8.*h(:,OP_1),OP_DZ,OP_DZ)
         endif
     
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
#ifdef USEST
             + prod(-g(:,OP_1),OP_DZP,OP_DZP) &
             + prod(-g(:,OP_1),OP_DRP,OP_DRP) &
             + prod(-g(:,OP_DP),OP_DZ,OP_DZP) &
             + prod(-g(:,OP_DP),OP_DR,OP_DRP)
#else
             + prod(g(:,OP_1),OP_DZ,OP_DZPP) &
             + prod(g(:,OP_1),OP_DR,OP_DRPP)
#endif
#endif
     end if

  v1umu = temp
end function v1umu

! V1vmu
! =====
!function v1vmu(e,f,g,h)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1vmu
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  !vectype, dimension(dofs_per_element) :: temp

  !temp = 0.

!#if defined(USE3D) || defined(USECOMPLEX)
     !if(surface_int) then
        !if(inoslip_pol.eq.1) then
           !temp = 0.
        !else
           !temp = intx5(e(:,:,OP_1),r_79,norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
                !- intx5(e(:,:,OP_1),r_79,norm79(:,1),f(:,OP_DZP),g(:,OP_1))
           !if(itor.eq.1) then
              !temp = temp &
                   !+ 2.*intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DP),g(:,OP_1)) &
                   !- 4.*intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DP),h(:,OP_1))
           !endif
        !end if
     !else 
        !temp = intx4(e(:,:,OP_DR),r_79,f(:,OP_DZP),g(:,OP_1)) &
             !- intx4(e(:,:,OP_DZ),r_79,f(:,OP_DRP),g(:,OP_1))
        
        !if(itor.eq.1) then
           !temp = temp &
                !+ 4.*intx3(e(:,:,OP_DZ),f(:,OP_DP),h(:,OP_1)) &
                !- 2.*intx3(e(:,:,OP_DZ),f(:,OP_DP),g(:,OP_1))
        !endif
     !endif

!#endif
  !v1vmu = temp
!end function v1vmu

function v1vmu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vmu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  temp%len = 0

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len = 0
        else
        end if
     else
        temp = prod( r_79*g(:,OP_1),OP_DR,OP_DZP) &
             + prod(-r_79*g(:,OP_1),OP_DZ,OP_DRP)
        if(itor.eq.1) then
           temp = temp &
                + prod( 4.*h(:,OP_1),OP_DZ,OP_DP) &
                + prod(-2.*g(:,OP_1),OP_DZ,OP_DP)
        endif
     endif
#endif
  v1vmu = temp
end function v1vmu



! V1chimu
! =======
! function v1chimu(e,f,g,h)

!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1chimu
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
!   vectype, dimension(dofs_per_element) :: temp

!      if(surface_int) then
!         if(inoslip_pol.eq.1) then
!            temp = 0.
!         else
!            temp = -2.* &
!                 (intx5(e(:,:,OP_DZ),ri_79,norm79(:,1),f(:,OP_DRR),g(:,OP_1)) &
!                 +intx5(e(:,:,OP_DZ),ri_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_1)) &
!                 -intx5(e(:,:,OP_DR),ri_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_1)) &
!                 -intx5(e(:,:,OP_DR),ri_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_1))) &
!                 + 2.* &
!                 (intx5(e(:,:,OP_DZ),ri_79,norm79(:,1),f(:,OP_GS),g(:,OP_1)) &
!                 -intx5(e(:,:,OP_DR),ri_79,norm79(:,2),f(:,OP_GS),g(:,OP_1)))
           
!            if(itor.eq.1) then
!               temp = temp + 2.* &
!                    (intx5(e(:,:,OP_DZ),ri2_79,norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
!                    +intx5(e(:,:,OP_DZ),ri2_79,norm79(:,2),f(:,OP_DZ),g(:,OP_1)) &
!                    +intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_1)) &
!                    +intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_1)) &
!                    +intx5(e(:,:,OP_DZ),ri2_79,norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
!                    -intx5(e(:,:,OP_DR),ri2_79,norm79(:,1),f(:,OP_DZ),g(:,OP_1)) &
!                    -2.*intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_GS),h(:,OP_1)) &
!                    -2.*intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DZ),g(:,OP_1)))
!            end if
           
! #if defined(USE3D) || defined(USECOMPLEX)
!            temp = temp &
!                 + intx5(e(:,:,OP_1),ri3_79,g(:,OP_1),norm79(:,2),f(:,OP_DRPP)) &
!                 - intx5(e(:,:,OP_1),ri3_79,g(:,OP_1),norm79(:,1),f(:,OP_DZPP))
! #endif
!         endif
!      else 
!         temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
!         temp79c = f(:,OP_DRZ)
!         if(itor.eq.1) then
!            temp79b = temp79b + 2.*ri_79*f(:,OP_DR)
!            temp79c = temp79c -    ri_79*f(:,OP_DZ)
!         endif
     
!         temp = -2.* &
!              (intx4(e(:,:,OP_DRZ),ri_79,temp79b,g(:,OP_1)) &
!              -intx4(e(:,:,OP_DZZ),ri_79,temp79c,g(:,OP_1)) &
!              +intx4(e(:,:,OP_DRR),ri_79,temp79c,g(:,OP_1)))
        
!         if(itor.eq.1) then
!            temp = temp &
!                 -    intx4(e(:,:,OP_DZ),ri2_79,temp79b,g(:,OP_1)) &
!                 - 2.*intx4(e(:,:,OP_DR),ri2_79,temp79c,g(:,OP_1)) &
!                 + 4.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),h(:,OP_1)) &
!                 - 3.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),g(:,OP_1))
!         endif
     
! #if defined(USE3D) || defined(USECOMPLEX)
!         temp = temp &
! #ifdef USEST
!              - intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZP),g(:,OP_1)) &
!              + intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DRP),g(:,OP_1)) &
!              - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_DP)) &
!              + intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_DP))
! #else
!              + intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZPP),g(:,OP_1)) &
!              - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRPP),g(:,OP_1))
! #endif
! #endif
!      endif

!   v1chimu = temp
! end function v1chimu

function v1chimu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chimu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  type(muarray) :: tempb, tempc

     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len = 0
        else
           temp =prod(-2.*ri_79*norm79(:,1)*g(:,OP_1),OP_DZ,OP_DRR) &
                +prod(-2.*ri_79*norm79(:,2)*g(:,OP_1),OP_DZ,OP_DRZ) &
                +prod( 2.*ri_79*norm79(:,1)*g(:,OP_1),OP_DR,OP_DRZ) &
                +prod( 2.*ri_79*norm79(:,2)*g(:,OP_1),OP_DR,OP_DZZ) &
                +prod( 2.*ri_79*norm79(:,1)*g(:,OP_1),OP_DZ,OP_GS) &
                +prod(-2.*ri_79*norm79(:,2)*g(:,OP_1),OP_DR,OP_GS)
 
           if(itor.eq.1) then
              temp = temp + &
                   (prod(2.*ri2_79*norm79(:,1)*g(:,OP_1),OP_DZ,OP_DR) &
                   +prod(2.*ri2_79*norm79(:,2)*g(:,OP_1),OP_DZ,OP_DZ) &
                   +prod(2.*ri2_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DRZ) &
                   +prod(2.*ri2_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DZZ) &
                   +prod(2.*ri2_79*norm79(:,1)*g(:,OP_1),OP_DZ,OP_DR) &
                   +prod(-2.*ri2_79*norm79(:,1)*g(:,OP_1),OP_DR,OP_DZ) &
                   +prod(-4.*ri2_79*norm79(:,2)*h(:,OP_1),OP_1,OP_GS) &
                   +prod(-4.*ri3_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DZ))
            end if
           
#if defined(USE3D) || defined(USECOMPLEX)
           temp = temp &
                + prod( ri3_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DRPP) &
                + prod(-ri3_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DZPP)
#endif
        endif
     else 
        temp79b(:) = 1.
        tempb = mu(temp79b,OP_DZZ) + mu(-temp79b,OP_DRR)
        tempc = mu(temp79b,OP_DRZ)
        if(itor.eq.1) then
           tempb = tempb + mu(2.*ri_79,OP_DR)
           tempc = tempc + mu(-ri_79,OP_DZ)
        endif
     
        temp = &
             (prod(mu(-2.*ri_79*g(:,OP_1),OP_DRZ),tempb) &
             +prod(mu( 2.*ri_79*g(:,OP_1),OP_DZZ),tempc) &
             +prod(mu(-2.*ri_79*g(:,OP_1),OP_DRR),tempc))
        
 
        if(itor.eq.1) then
           temp = temp &
                + prod(mu(-   ri2_79*g(:,OP_1),OP_DZ),tempb) &
                + prod(mu(-2.*ri2_79*g(:,OP_1),OP_DR),tempc) &
                + prod( 4.*ri2_79*h(:,OP_1),OP_DZ,OP_GS) &
                + prod(-3.*ri2_79*g(:,OP_1),OP_DZ,OP_GS)
         endif
     
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
#ifdef USEST
             + prod(-ri3_79*g(:,OP_1),OP_DRP,OP_DZP) &
             + prod( ri3_79*g(:,OP_1),OP_DZP,OP_DRP) &
             + prod(-ri3_79*g(:,OP_DP),OP_DR,OP_DZP) &
             + prod( ri3_79*g(:,OP_DP),OP_DZ,OP_DRP)
#else
             + prod( ri3_79*g(:,OP_1),OP_DR,OP_DZPP) &
             + prod(-ri3_79*g(:,OP_1),OP_DZ,OP_DRPP)
#endif
#endif
     endif

  v1chimu = temp
end function v1chimu



! V1un
! ====
!function v1un(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1un
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then 
        !if(inoslip_pol.eq.1) then
           !temp = 0.
        !else
           !temp = &
                !- intx5(e(:,:,OP_1),r2_79,norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
                !- intx5(e(:,:,OP_1),r2_79,norm79(:,2),f(:,OP_DZ),g(:,OP_1))
        !endif
     !else
        !temp = intx4(e(:,:,OP_DR),r2_79,f(:,OP_DR),g(:,OP_1)) &
             !+ intx4(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),g(:,OP_1))
     !end if

  !v1un = temp
!end function v1un

function v1un(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1un
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then 
           temp%len = 0.
     else
        temp = prod(r2_79*g(:,OP_1),OP_DR,OP_DR) &
             + prod(r2_79*g(:,OP_1),OP_DZ,OP_DZ)
     end if

  v1un = temp
end function v1un


! V1chin
! ======
!function v1chin(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1chin
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !if(inoslip_pol.eq.1) then
           !temp = 0.
        !else
           !temp = intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DR),g(:,OP_1)) &
                !- intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DZ),g(:,OP_1))
        !end if
     !else        
        !temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
             !- intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1))
     !endif

  !v1chin = temp
!end function v1chin

function v1chin(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chin
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len = 0
        else
           temp = prod( ri_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DR) &
                + prod(-ri_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DZ)
         end if
     else        
        temp = prod( ri_79*g(:,OP_1),OP_DR,OP_DZ) &
             + prod(-ri_79*g(:,OP_1),OP_DZ,OP_DR)
      endif

  v1chin = temp
end function v1chin


! V1psipsi
! ========
function v1psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1psipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DZ),g(:,OP_GS)) &
                - intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DR),g(:,OP_GS))
        endif
     else
        temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_GS)) &
             - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_GS))
     endif


  v1psipsi = temp
end function v1psipsi

function v1psipsi1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1psipsi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then
           temp%len = 0
     else
        temp = prod(ri_79*g(:,OP_GS),OP_DZ,OP_DR) &
             + prod(-ri_79*g(:,OP_GS),OP_DR,OP_DZ)
     endif

  v1psipsi1 = temp
end function v1psipsi1

function v1psipsi2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1psipsi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

     if(surface_int) then
           temp%len = 0
     else
        temp = prod(ri_79*f(:,OP_DR),OP_DZ,OP_GS) &
             + prod(-ri_79*f(:,OP_DZ),OP_DR,OP_GS)
     endif

  v1psipsi2 = temp
end function v1psipsi2


! V1psib
! ======
function v1psib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1psib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
                - intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
                - intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZP),g(:,OP_1))
        end if
     else
        temp = intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
             + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1))
     endif

  v1psib = temp
#else
  v1psib = 0.
#endif
end function v1psib

function v1psib1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1psib1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
           temp%len = 0
     else
        temp = prod(ri2_79*g(:,OP_1),OP_DZ,OP_DZP) &
             + prod(ri2_79*g(:,OP_1),OP_DR,OP_DRP)
     endif

  v1psib1 = temp
#else
  v1psib1%len = 0
#endif
end function v1psib1

function v1psib2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1psib2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len = 0
        else
        end if
     else
        temp = prod(ri2_79*f(:,OP_DZP),OP_DZ,OP_1) &
             + prod(ri2_79*f(:,OP_DRP),OP_DR,OP_1)
     endif

  v1psib2 = temp
#else
  v1psib2%len = 0.
#endif
end function v1psib2


! V1bb
! ====
function v1bb(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1bb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        ! Seems to work better without this...
        temp = 0.
!!$        temp = 0.25* &
!!$             (int5(ri_79,norm79(:,1),e(:,OP_DZ),f(:,OP_1),g(:,OP_1)) &
!!$             -int5(ri_79,norm79(:,2),e(:,OP_DR),f(:,OP_1),g(:,OP_1))) &
!!$             + 0.5 * &
!!$             (int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1)) &
!!$             -int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_1)))
     else
        temp = 0.
     endif

  v1bb = temp
end function v1bb


! V1uun 
! =====
function v1uun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1uun
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1uun = 0.
     return
  end if

     if(surface_int) then
        temp = 0.
     else
        temp = -intx5(e(:,:,OP_DZ),r3_79,f(:,OP_DR),g(:,OP_LP),h(:,OP_1)) &
             +  intx5(e(:,:,OP_DR),r3_79,f(:,OP_DZ),g(:,OP_LP),h(:,OP_1)) &
             -  intx5(e(:,:,OP_DR),r3_79,f(:,OP_DRZ),g(:,OP_DR),h(:,OP_1)) &
             -  intx5(e(:,:,OP_DR),r3_79,f(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1)) &
             +  intx5(e(:,:,OP_DZ),r3_79,f(:,OP_DRR),g(:,OP_DR),h(:,OP_1)) &
             +  intx5(e(:,:,OP_DZ),r3_79,f(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp &
                + intx5(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
                + intx5(e(:,:,OP_DZ),r2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
        end if
     end if


  v1uun = temp
end function v1uun

function v1uun1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uun1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1uun1%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp =  prod(-r3_79*g(:,OP_LP)*h(:,OP_1),OP_DZ,OP_DR) &
             +  prod( r3_79*g(:,OP_LP)*h(:,OP_1),OP_DR,OP_DZ) &
             +  prod(-r3_79*g(:,OP_DR)*h(:,OP_1),OP_DR,OP_DRZ) &
             +  prod(-r3_79*g(:,OP_DZ)*h(:,OP_1),OP_DR,OP_DZZ) &
             +  prod( r3_79*g(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DRR) &
             +  prod( r3_79*g(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DRZ)


        if(itor.eq.1) then
           temp = temp &
                + prod(r2_79*g(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DZ) &
                + prod(r2_79*g(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DR)
         end if
     end if

  v1uun1 = temp
end function v1uun1

function v1uun2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uun2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1uun2%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp =  prod(-r3_79*f(:,OP_DR)*h(:,OP_1),OP_DZ,OP_LP) &
             +  prod( r3_79*f(:,OP_DZ)*h(:,OP_1),OP_DR,OP_LP) &
             +  prod(-r3_79*f(:,OP_DRZ)*h(:,OP_1),OP_DR,OP_DR) &
             +  prod(-r3_79*f(:,OP_DZZ)*h(:,OP_1),OP_DR,OP_DZ) &
             +  prod( r3_79*f(:,OP_DRR)*h(:,OP_1),OP_DZ,OP_DR) &
             +  prod( r3_79*f(:,OP_DRZ)*h(:,OP_1),OP_DZ,OP_DZ)

        if(itor.eq.1) then
           temp = temp &
                + prod(r2_79*f(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DZ) &
                + prod(r2_79*f(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DR)
      end if
     end if

  v1uun2 = temp
end function v1uun2

function v1uun3(f,g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uun3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1uun3%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp =  prod(-r3_79*f(:,OP_DR)*g(:,OP_LP),OP_DZ,OP_1) &
             +  prod( r3_79*f(:,OP_DZ)*g(:,OP_LP),OP_DR,OP_1) &
             +  prod(-r3_79*f(:,OP_DRZ)*g(:,OP_DR),OP_DR,OP_1) &
             +  prod(-r3_79*f(:,OP_DZZ)*g(:,OP_DZ),OP_DR,OP_1) &
             +  prod( r3_79*f(:,OP_DRR)*g(:,OP_DR),OP_DZ,OP_1) &
             +  prod( r3_79*f(:,OP_DRZ)*g(:,OP_DZ),OP_DZ,OP_1)

        if(itor.eq.1) then
           temp = temp &
                + prod(r2_79*f(:,OP_DZ)*g(:,OP_DZ),OP_DZ,OP_1) &
                + prod(r2_79*f(:,OP_DR)*g(:,OP_DR),OP_DZ,OP_1)
        end if
     end if


  v1uun3 = temp
end function v1uun3

! V1uvn
! =====
function v1uvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1uvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1uvn = 0.
     return
  end if

     if(surface_int) then
        temp = 0.
     else
        temp = &
             - intx5(e(:,:,OP_DZ),r2_79,f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
             - intx5(e(:,:,OP_DR),r2_79,f(:,OP_DRP),g(:,OP_1),h(:,OP_1))
     end if

#else
  temp = 0.
#endif
  v1uvn = temp
end function v1uvn

function v1uvn1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uvn1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1uvn1%len = 0
     return
  end if
     if(surface_int) then
        temp%len = 0
     else
        temp = &
               prod(-r2_79*g(:,OP_1)*h(:,OP_1),OP_DZ,OP_DZP) &
             + prod(-r2_79*g(:,OP_1)*h(:,OP_1),OP_DR,OP_DRP)
      end if
#else
  temp%len = 0
#endif
  v1uvn1 = temp
end function v1uvn1

function v1uvn2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uvn2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1uvn2%len = 0
     return
  end if
     if(surface_int) then
        temp%len = 0
     else
        temp = &
               prod(-r2_79*f(:,OP_DZP)*h(:,OP_1),OP_DZ,OP_1) &
             + prod(-r2_79*f(:,OP_DRP)*h(:,OP_1),OP_DR,OP_1)
      end if
#else
  temp%len = 0
#endif
  v1uvn2 = temp
end function v1uvn2

function v1uvn3(f,g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uvn3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1uvn3%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp = &
               prod(-r2_79*f(:,OP_DZP)*g(:,OP_1),OP_DZ,OP_1) &
             + prod(-r2_79*f(:,OP_DRP)*g(:,OP_1),OP_DR,OP_1)
     end if

#else
  temp%len = 0
#endif
  v1uvn3 = temp
end function v1uvn3

! v1uchin
! =======
function v1uchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1uchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1uchin = 0.
     return
  end if


     if(surface_int) then
        temp = 0.
     else
        temp79a = (f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &
             +     f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ))
        temp79b = (f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &
             +     f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ))
        temp = -intx4(e(:,:,OP_DR),f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) &
             -  intx4(e(:,:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
             -  intx3(e(:,:,OP_DZ),temp79a,h(:,OP_1)) &
             +  intx3(e(:,:,OP_DR),temp79b,h(:,OP_1))

        if(itor.eq.1) then
           temp = temp &
                - 2.*intx5(e(:,:,OP_DR),ri_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))&
                - 2.*intx5(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1))&
                + intx5(e(:,:,OP_DZ),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
                - intx5(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
        end if
     end if

  v1uchin = temp
end function v1uchin

function v1uchin1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uchin1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(inertia.eq.0) then
     v1uchin1%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        tempa = mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &
           +    mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR)
        tempb = mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &
           +    mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR)
        temp =  prod(-g(:,OP_DR)*h(:,OP_1),OP_DR,OP_GS) &
             +  prod(-g(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_GS) &
             +  prod(mu(-h(:,OP_1),OP_DZ),tempa) &
             +  prod(mu(h(:,OP_1),OP_DR),tempb)

        if(itor.eq.1) then
           temp = temp &
                + prod(-2.*ri_79*g(:,OP_DR)*h(:,OP_1),OP_DR,OP_DR)&
                + prod(-2.*ri_79*g(:,OP_DZ)*h(:,OP_1),OP_DR,OP_DZ)&
                + prod( ri_79*g(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DZ) &
                + prod(-ri_79*g(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DR)
        end if
     end if

  v1uchin1 = temp
end function v1uchin1

function v1uchin2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uchin2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(inertia.eq.0) then
     v1uchin2%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        tempa = (mu(f(:,OP_DRZ),OP_DR) + mu(-f(:,OP_DRR),OP_DZ) &
           +     mu(f(:,OP_DZ),OP_DRR) + mu(-f(:,OP_DR),OP_DRZ))
        tempb = (mu(f(:,OP_DZZ),OP_DR) + mu(-f(:,OP_DRZ),OP_DZ) &
           +     mu(f(:,OP_DZ),OP_DRZ) + mu(-f(:,OP_DR),OP_DZZ))
        temp =  prod(-f(:,OP_GS)*h(:,OP_1),OP_DR,OP_DR) &
             +  prod(-f(:,OP_GS)*h(:,OP_1),OP_DZ,OP_DZ) &
             +  prod(mu(-h(:,OP_1),OP_DZ),tempa) &
             +  prod(mu( h(:,OP_1),OP_DR),tempb)

        if(itor.eq.1) then
           temp = temp &
                + prod(-2.*ri_79*f(:,OP_DR)*h(:,OP_1),OP_DR,OP_DR)&
                + prod(-2.*ri_79*f(:,OP_DZ)*h(:,OP_1),OP_DR,OP_DZ)&
                + prod( ri_79*f(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DR) &
                + prod(-ri_79*f(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DZ)
        end if
     end if

  v1uchin2 = temp
end function v1uchin2

function v1uchin3(f,g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uchin3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1uchin3%len = 0
     return
  end if


     if(surface_int) then
        temp%len = 0
     else
        temp79a = (f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &
             +     f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ))
        temp79b = (f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &
             +     f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ))
        temp79c = 1.
        temp =  prod(-f(:,OP_GS)*g(:,OP_DR),OP_DR,OP_1) &
             +  prod(-f(:,OP_GS)*g(:,OP_DZ),OP_DZ,OP_1) &
             +  prod(-temp79c,OP_DZ,OP_1)*temp79a &
             +  prod( temp79c,OP_DR,OP_1)*temp79b

        if(itor.eq.1) then
           temp = temp &
                + prod(-2.*ri_79*f(:,OP_DR)*g(:,OP_DR),OP_DR,OP_1)&
                + prod(-2.*ri_79*f(:,OP_DZ)*g(:,OP_DZ),OP_DR,OP_1)&
                + prod( ri_79*f(:,OP_DZ)*g(:,OP_DR),OP_DZ,OP_1) &
                + prod(-ri_79*f(:,OP_DR)*g(:,OP_DZ),OP_DZ,OP_1)
        end if
     end if

  v1uchin3 = temp
end function v1uchin3

! V1vvn
! =====
function v1vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1vvn = 0.
     return
  end if

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        if(itor.eq.1) then
           temp = -intx5(e(:,:,OP_DZ),r2_79,f(:,OP_1),g(:,OP_1),h(:,OP_1))
        endif
     end if

  v1vvn = temp
end function v1vvn

function v1vvn1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vvn1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1vvn1%len = 0
     return
  end if

  temp%len = 0
     if(surface_int) then
        temp%len = 0
     else
        if(itor.eq.1) then
           temp = prod(-r2_79*g(:,OP_1)*h(:,OP_1),OP_DZ,OP_1)
        endif
     end if

  v1vvn1 = temp
end function v1vvn1

function v1vvn2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vvn2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1vvn2%len = 0
     return
  end if

  temp%len = 0
     if(surface_int) then
        temp%len = 0
     else
        if(itor.eq.1) then
           temp = prod(-r2_79*f(:,OP_1)*h(:,OP_1),OP_DZ,OP_1)
        endif
     end if

  v1vvn2 = temp
end function v1vvn2

function v1vvn3(f,g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vvn3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1vvn3%len = 0
     return
  end if

  temp%len = 0
     if(surface_int) then
        temp%len = 0
     else
        if(itor.eq.1) then
           temp = prod(-r2_79*f(:,OP_1)*g(:,OP_1),OP_DZ,OP_1)
        endif
     end if

  v1vvn3 = temp
end function v1vvn3

! V1vchin
! =======
function v1vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1vchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1vchin = 0.
     return
  end if

     if(surface_int) then
        temp = 0.
     else
        temp = intx5(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1)) &
             - intx5(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1))
     endif

#else
  temp = 0.
#endif
  v1vchin = temp
end function v1vchin

function v1vchin1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vchin1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1vchin1%len = 0
     return
  end if
     if(surface_int) then
        temp%len = 0
     else
        temp = prod( ri_79*g(:,OP_DRP)*h(:,OP_1),OP_DZ,OP_1) &
             + prod(-ri_79*g(:,OP_DZP)*h(:,OP_1),OP_DR,OP_1)
     endif
#else
  temp%len = 0
#endif
  v1vchin1 = temp
end function v1vchin1

function v1vchin2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vchin2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1vchin2%len = 0
     return
  end if
     if(surface_int) then
        temp%len = 0
     else
        temp = prod( ri_79*f(:,OP_1)*h(:,OP_1),OP_DZ,OP_DRP) &
             + prod(-ri_79*f(:,OP_1)*h(:,OP_1),OP_DR,OP_DZP)
     endif
#else
  temp%len = 0
#endif
  v1vchin2 = temp
end function v1vchin2

function v1vchin3(f,g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vchin3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1vchin3%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp = prod( ri_79*f(:,OP_1)*g(:,OP_DRP),OP_DZ,OP_1) &
             + prod(-ri_79*f(:,OP_1)*g(:,OP_DZP),OP_DR,OP_1)
     endif

#else
  temp%len = 0
#endif
  v1vchin3 = temp
end function v1vchin3

! v1chichin
! =========
function v1chichin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chichin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0) then
     v1chichin = 0.
     return
  end if

     if(surface_int) then
        temp = 0.
     else
        temp = -intx5(e(:,:,OP_DR),ri3_79,f(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1)) &
             -  intx5(e(:,:,OP_DR),ri3_79,f(:,OP_DRZ),g(:,OP_DR),h(:,OP_1)) &
             +  intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1)) &
             -  intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_DRR),g(:,OP_DR),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp + 2.*&
                (intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) & 
                +intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)))
        end if
     endif
  
  v1chichin = temp
end function v1chichin

function v1chichin1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chichin1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1chichin1%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp =  prod(-ri3_79*g(:,OP_DZ)*h(:,OP_1),OP_DR,OP_DZZ) &
             +  prod(-ri3_79*g(:,OP_DR)*h(:,OP_1),OP_DR,OP_DRZ) &
             +  prod( ri3_79*g(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DRZ) &
             +  prod(-ri3_79*g(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DRR)

        if(itor.eq.1) then
           temp = temp + &
                (prod(2.*ri4_79*g(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DZ) & 
                +prod(2.*ri4_79*g(:,OP_DZ)*h(:,OP_1),OP_DR,OP_DR))
        end if
     endif
  
  v1chichin1 = temp
end function v1chichin1

function v1chichin2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chichin2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1chichin2%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp =  prod(-ri3_79*f(:,OP_DZZ)*h(:,OP_1),OP_DR,OP_DZ) &
             +  prod(-ri3_79*f(:,OP_DRZ)*h(:,OP_1),OP_DR,OP_DR) &
             +  prod( ri3_79*f(:,OP_DRZ)*h(:,OP_1),OP_DZ,OP_DZ) &
             +  prod(-ri3_79*f(:,OP_DRR)*h(:,OP_1),OP_DZ,OP_DR)

        if(itor.eq.1) then
           temp = temp + &
                (prod(2.*ri4_79*f(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DZ) & 
                +prod(2.*ri4_79*f(:,OP_DR)*h(:,OP_1),OP_DR,OP_DZ))
        end if
     endif
  
  v1chichin2 = temp
end function v1chichin2

function v1chichin3(f,g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chichin3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

  if(inertia.eq.0) then
     v1chichin3%len = 0
     return
  end if

     if(surface_int) then
        temp%len = 0
     else
        temp =  prod(-ri3_79*f(:,OP_DZZ)*g(:,OP_DZ),OP_DR,OP_1) &
             +  prod(-ri3_79*f(:,OP_DRZ)*g(:,OP_DR),OP_DR,OP_1) &
             +  prod( ri3_79*f(:,OP_DRZ)*g(:,OP_DZ),OP_DZ,OP_1) &
             +  prod(-ri3_79*f(:,OP_DRR)*g(:,OP_DR),OP_DZ,OP_1)

        if(itor.eq.1) then
           temp = temp + 2.*&
                (prod(ri4_79*f(:,OP_DZ)*g(:,OP_DZ),OP_DZ,OP_1) & 
                +prod(ri4_79*f(:,OP_DR)*g(:,OP_DZ),OP_DR,OP_1))
        end if
     endif
  
  v1chichin3 = temp
end function v1chichin3

! V1upsipsi
! =========
!function v1upsipsi(e,f,g,h)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1upsipsi
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  !vectype, dimension(dofs_per_element) :: temp
  
  !vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe, tempf
  !integer :: j
  

  !! |u, psi(1)|,r
  !temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ ) &
       !+    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
  !! |u, psi(1)|,z
  !temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ ) &
       !+    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)

  !do j=1, dofs_per_element
     !! |nu, psi(2)|,r
     !tempe(j,:) = e(j,:,OP_DRZ)*h(:,OP_DR ) - e(j,:,OP_DRR)*h(:,OP_DZ ) &
          !+       e(j,:,OP_DZ )*h(:,OP_DRR) - e(j,:,OP_DR )*h(:,OP_DRZ)
     !! |nu, psi(2)|,z
     !tempf(j,:) = e(j,:,OP_DZZ)*h(:,OP_DR ) - e(j,:,OP_DRZ)*h(:,OP_DZ ) &
          !+       e(j,:,OP_DZ )*h(:,OP_DRZ) - e(j,:,OP_DR )*h(:,OP_DZZ)
  !end do

     !if(surface_int) then
        !temp = intx4(e(:,:,OP_DZ),h(:,OP_DR),norm79(:,1),temp79b) &
             !- intx4(e(:,:,OP_DR),h(:,OP_DZ),norm79(:,1),temp79b) &
             !+ intx4(e(:,:,OP_DZ),h(:,OP_DR),norm79(:,2),temp79c) &
             !- intx4(e(:,:,OP_DR),h(:,OP_DZ),norm79(:,2),temp79c)
        !if(itor.eq.1) then
           !temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           !temp = temp &
                !+ intx5(e(:,:,OP_DZ),h(:,OP_DR),ri_79,norm79(:,1),temp79d) &
                !- intx5(e(:,:,OP_DR),h(:,OP_DZ),ri_79,norm79(:,1),temp79d)
        !endif
     !else
        !temp = -intx2(tempe,temp79b)  &
             !- intx2(tempf,temp79c)  &
             !+ intx3(e(:,:,OP_DZ),temp79b,h(:,OP_GS))  &
             !- intx3(e(:,:,OP_DR),temp79c,h(:,OP_GS))

        !if(itor.eq.1) then
           !! |u, psi(1)|
           !temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           
           !do j=1, dofs_per_element
              !! |nu,psi(2)|
              !tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
           !end do
           !temp = temp             &
                !- intx3(tempd,ri_79,temp79b)   &
                !- intx3(tempe,ri_79,temp79a)   &
                !- intx3(tempd,ri2_79,temp79a)   &
                !+ intx4(e(:,:,OP_DZ),ri_79,h(:,OP_GS),temp79a)
        !endif
     !endif

  !v1upsipsi = temp
!end function v1upsipsi

function v1upsipsi(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1upsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  
  type(muarray) :: tempa, tempb, tempc
  type(muarray) :: tempd, tempe, tempf
  integer :: j
  

  ! |u, psi(1)|,r
  tempb = mu(g(:,OP_DR ),OP_DRZ) + mu(-g(:,OP_DZ ),OP_DRR) &
       +  mu(g(:,OP_DRR),OP_DZ)  + mu(-g(:,OP_DRZ),OP_DR)
  ! |u, psi(1)|,z
  tempc = mu(g(:,OP_DR ),OP_DZZ) + mu(-g(:,OP_DZ ),OP_DRZ) &
       +  mu(g(:,OP_DRZ),OP_DZ)  + mu(-g(:,OP_DZZ),OP_DR)

  ! |nu, psi(2)|,r
  tempe = mu(-h(:,OP_DR ),OP_DRZ) + mu(h(:,OP_DZ ),OP_DRR) &
        + mu(-h(:,OP_DRR),OP_DZ)  + mu(h(:,OP_DRZ),OP_DR)
  ! |nu, psi(2)|,z
  tempf = mu(-h(:,OP_DR ),OP_DZZ) + mu(h(:,OP_DZ ),OP_DRZ) &
        + mu(-h(:,OP_DRZ),OP_DZ)  + mu(h(:,OP_DZZ),OP_DR)

     if(surface_int) then
        temp = prod(mu( h(:,OP_DR)*norm79(:,1),OP_DZ),tempb) &
             + prod(mu(-h(:,OP_DZ)*norm79(:,1),OP_DR),tempb) &
             + prod(mu( h(:,OP_DR)*norm79(:,2),OP_DZ),tempc) &
             + prod(mu(-h(:,OP_DZ)*norm79(:,2),OP_DR),tempc)
        if(itor.eq.1) then
           tempd = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
           temp = temp &
                + prod(mu( h(:,OP_DR)*ri_79*norm79(:,1),OP_DZ),tempd) &
                + prod(mu(-h(:,OP_DZ)*ri_79*norm79(:,1),OP_DR),tempd)
        endif
     else
        temp = prod(tempe,tempb)  &
             + prod(tempf,tempc)  &
             + prod(mu( h(:,OP_GS),OP_DZ),tempb)  &
             + prod(mu(-h(:,OP_GS),OP_DR),tempc)
     
        if(itor.eq.1) then
           ! |u, psi(1)|
           tempa = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
           ! |nu,psi(2)|
           tempd = mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR)
           temp = temp             &
                + prod(tempd*(-ri_79),tempb)   &
                + prod(tempe*(ri_79),tempa)   &
                + prod(tempd*(-ri2_79),tempa)   &
                + prod(mu(h(:,OP_GS)*ri_79,OP_DZ),tempa)
        endif
     endif

  v1upsipsi = temp
end function v1upsipsi


! V1upsib
! =======
! function v1upsib(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1upsib
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

! #if defined(USE3D) || defined(USECOMPLEX)
!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempc, tempd
!   integer :: j

!      if(surface_int) then
!         temp79a = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
!              +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)
!         temp79b = ri_79*h(:,OP_DP)
!         temp = intx5(e(:,:,OP_DR),ri_79,norm79(:,1),temp79a,h(:,OP_1)) &
!              + intx5(e(:,:,OP_DZ),ri_79,norm79(:,2),temp79a,h(:,OP_1)) &
!              + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
!              - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
!              + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
!              - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
!              + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
!              - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
!              + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
!              - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ))
!         if(itor.eq.1) then
!            temp79c = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
!            temp = temp &
!                 + intx5(e(:,:,OP_1),ri_79,temp79b,temp79c,norm79(:,1))
!         endif
!      else
!         do j=1, dofs_per_element
!            tempa(j,:) = h(:,OP_1)*e(j,:,OP_GS)  &
!                 + h(:,OP_DZ)*e(j,:,OP_DZ) + h(:,OP_DR)*e(j,:,OP_DR)
!            tempc(j,:) = e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP)
!            tempd(j,:) = h(:,OP_DP)* &
!                 (e(j,:,OP_DZ)*f(:,OP_DR )-e(j,:,OP_DR)*f(:,OP_DZ )) &
!                 +    h(:,OP_1 )* &
!                 (e(j,:,OP_DZ)*f(:,OP_DRP)-e(j,:,OP_DR)*f(:,OP_DZP))
!         end do
!         temp79b = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
!              +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)
! #ifdef USEST
!         temp79e = h(:,OP_1)*f(:,OP_GS) &
!              + h(:,OP_DZ)*f(:,OP_DZ ) + h(:,OP_DR)*f(:,OP_DR ) 
! #else
!         temp79e = h(:,OP_DP)*f(:,OP_GS) + h(:,OP_1)*f(:,OP_GSP) &
!              + h(:,OP_DZP)*f(:,OP_DZ ) + h(:,OP_DRP)*f(:,OP_DR ) &
!              + h(:,OP_DZ )*f(:,OP_DZP) + h(:,OP_DR )*f(:,OP_DRP)
! #endif
     
!         temp = intx3(tempa,ri_79,temp79b) &
!              + intx4(tempc,ri_79,f(:,OP_DR),h(:,OP_DZ)) &
!              - intx4(tempc,ri_79,f(:,OP_DZ),h(:,OP_DR)) &
!              - intx3(tempd,ri_79,g(:,OP_GS)) &
! #ifdef USEST
!              + intx4(e(:,:,OP_DZP),ri_79,g(:,OP_DR),temp79e) &
!              - intx4(e(:,:,OP_DRP),ri_79,g(:,OP_DZ),temp79e) &
!              + intx4(e(:,:,OP_DZ),ri_79,g(:,OP_DRP),temp79e) &
!              - intx4(e(:,:,OP_DR),ri_79,g(:,OP_DZP),temp79e)
! #else
!              - intx4(e(:,:,OP_DZ),ri_79,g(:,OP_DR),temp79e) &
!              + intx4(e(:,:,OP_DR),ri_79,g(:,OP_DZ),temp79e)
! #endif
!         temp = -temp
!      endif

!   v1upsib = temp
! #else
!   v1upsib = 0.
! #endif

! end function v1upsib

function v1upsib(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1upsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempe
  type(prodarray) :: tempd
  integer :: j

     if(surface_int) then
        tempa = mu(g(:,OP_DR ),OP_DZP) + mu(-g(:,OP_DZ ),OP_DRP) &
              + mu(g(:,OP_DRP),OP_DZ)  + mu(-g(:,OP_DZP),OP_DR)
        temp79b = ri_79*h(:,OP_DP)
        temp = prod(mu(ri_79*norm79(:,1)*h(:,OP_1),OP_DR),tempa) &
             + prod(mu(ri_79*norm79(:,2)*h(:,OP_1),OP_DZ),tempa) &
             + prod( temp79b*norm79(:,1)*g(:,OP_DR ),OP_1,OP_DRZ) &
             + prod(-temp79b*norm79(:,1)*g(:,OP_DZ ),OP_1,OP_DRR) &
             + prod( temp79b*norm79(:,1)*g(:,OP_DRR),OP_1,OP_DZ) &
             + prod(-temp79b*norm79(:,1)*g(:,OP_DRZ),OP_1,OP_DR) &
             + prod( temp79b*norm79(:,2)*g(:,OP_DR ),OP_1,OP_DZZ) &
             + prod(-temp79b*norm79(:,2)*g(:,OP_DZ ),OP_1,OP_DRZ) &
             + prod( temp79b*norm79(:,2)*g(:,OP_DRZ),OP_1,OP_DZ) &
             + prod(-temp79b*norm79(:,2)*g(:,OP_DZZ),OP_1,OP_DR)
        if(itor.eq.1) then
           tempc = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
           temp = temp &
                + prod(mu(ri_79*temp79b*norm79(:,1),OP_1),tempc)
         endif
     else
        tempa = mu(h(:,OP_1),OP_GS)  &
              + mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
        tempc = mu(g(:,OP_DZP),OP_DZ) + mu(g(:,OP_DRP),OP_DR)
        tempd = prod(h(:,OP_DP),OP_DZ,OP_DR )+prod(-h(:,OP_DP),OP_DR,OP_DZ ) &
              + prod(h(:,OP_1),OP_DZ,OP_DRP)+prod(-h(:,OP_1),OP_DR,OP_DZP)
        tempb = mu(g(:,OP_DR ),OP_DZP) + mu(-g(:,OP_DZ ),OP_DRP) &
              + mu(g(:,OP_DRP),OP_DZ)  + mu(-g(:,OP_DZP),OP_DR)

#ifdef USEST
        tempe = mu(h(:,OP_1),OP_GS) &
             +  mu(h(:,OP_DZ),OP_DZ ) + mu(h(:,OP_DR),OP_DR ) 
#else
        tempe = mu(h(:,OP_DP),OP_GS)  + mu(h(:,OP_1),OP_GSP) &
             +  mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR ) &
             +  mu(h(:,OP_DZ),OP_DZP) + mu(h(:,OP_DR ),OP_DRP)
#endif
     
        temp = prod(tempa,tempb*(-ri_79)) &
             + prod(tempc,mu(-ri_79*h(:,OP_DZ),OP_DR)) &
             + prod(tempc,mu( ri_79*h(:,OP_DR),OP_DZ)) &
             + tempd*(ri_79*g(:,OP_GS)) &
#ifdef USEST
             + prod(mu(-ri_79*g(:,OP_DR),OP_DZP),tempe) &
             + prod(mu( ri_79*g(:,OP_DZ),OP_DRP),tempe) &
             + prod(mu(-ri_79*g(:,OP_DRP),OP_DZ),tempe) &
             + prod(mu( ri_79*g(:,OP_DZP),OP_DR),tempe)
#else
             + prod(mu( ri_79*g(:,OP_DR),OP_DZ),tempe) &
             + prod(mu(-ri_79*g(:,OP_DZ),OP_DR),tempe)
#endif
     endif

  v1upsib = temp
#else
  v1upsib%len = 0
#endif

end function v1upsib


! V1ubb 
! =====
!function v1ubb(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v1ubb
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
!  integer :: j

!  temp = 0.

!     if(surface_int) then
!        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
!        temp = intx4(e(:,:,OP_1),temp79a,norm79(:,2),h(:,OP_DR)) &
!             - intx4(e(:,:,OP_1),temp79a,norm79(:,1),h(:,OP_DZ))
!#if defined(USE3D) || defined(USECOMPLEX)
!        temp79a = ri2_79*h(:,OP_DP)
!        temp = temp &
!             + intx5(e(:,:,OP_1),temp79a,norm79(:,1),f(:,OP_DR),g(:,OP_DP)) &
!             + intx5(e(:,:,OP_1),temp79a,norm79(:,2),f(:,OP_DZ),g(:,OP_DP)) &
!             + intx5(e(:,:,OP_1),temp79a,norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
!             + intx5(e(:,:,OP_1),temp79a,norm79(:,2),f(:,OP_DZP),g(:,OP_1))
!#endif
!     else
!#if defined(USE3D) || defined(USECOMPLEX)
!#ifdef USEST 
!        do j=1, dofs_per_element
!           tempa(j,:) = &
!                -(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
!                *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) &
!                -(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
!                *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP))
!        end do
!        temp = intx2(tempa,ri2_79)
!        !temp = intx5(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1),ri2_79,h(:,OP_DPP))&
!        !     + intx5(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1),ri2_79,h(:,OP_DPP))
!#else
!        do j=1, dofs_per_element
!           tempa(j,:) = &
!                (e(j,:,OP_DZ)*f(:,OP_DZPP) + e(j,:,OP_DR)*f(:,OP_DRPP)) &
!                *g(:,OP_1) &
!                + 2.*(e(j,:,OP_DZ)*f(:,OP_DZP) + e(j,:,OP_DR)*f(:,OP_DRP)) &
!                *g(:,OP_DP) &
!                +    (e(j,:,OP_DZ)*f(:,OP_DZ) + e(j,:,OP_DR)*f(:,OP_DR)) &
!                *g(:,OP_DPP)
!        end do
!        temp = intx3(tempa,ri2_79,h(:,OP_1))
!#endif
!#endif
!     end if

!  v1ubb = temp
!end function v1ubb

function v1ubb(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1ubb
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp, tempb
  type(muarray) :: tempa
  integer :: j

  temp%len = 0.

     if(surface_int) then
        tempa = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
        temp = prod(mu( norm79(:,2)*h(:,OP_DR),OP_1),tempa) &
             + prod(mu(-norm79(:,1)*h(:,OP_DZ),OP_1),tempa)
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = ri2_79*h(:,OP_DP)
        temp = temp &
             + prod(temp79a*norm79(:,1)*g(:,OP_DP),OP_1,OP_DR) &
             + prod(temp79a*norm79(:,2)*g(:,OP_DP),OP_1,OP_DZ) &
             + prod(temp79a*norm79(:,1)*g(:,OP_1),OP_1,OP_DRP) &
             + prod(temp79a*norm79(:,2)*g(:,OP_1),OP_1,OP_DZP)
#endif
     else
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST 
        tempb = &
                prod(mu(-h(:,OP_DP),OP_DZ) + mu(-h(:,OP_1),OP_DZP), &
                     mu(g(:,OP_1),OP_DZP) + mu(g(:,OP_DP),OP_DZ)) &
              + prod(mu(-h(:,OP_DP),OP_DR) + mu(-h(:,OP_1),OP_DRP), &
                     mu(g(:,OP_1),OP_DRP) + mu(g(:,OP_DP),OP_DR))
        temp = tempb*(ri2_79)
        !temp = intx5(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1),ri2_79,h(:,OP_DPP))&
        !     + intx5(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1),ri2_79,h(:,OP_DPP))
#else
        tempb = prod(g(:,OP_1),OP_DZ,OP_DZPP) + prod(g(:,OP_1),OP_DR,OP_DRPP) &
              + prod(2.*g(:,OP_DP),OP_DZ,OP_DZP) + prod(2.*g(:,OP_DP),OP_DR,OP_DRP) &
              + prod(g(:,OP_DPP),OP_DZ,OP_DZ) + prod(g(:,OP_DPP),OP_DR,OP_DR)
        temp = tempb*(ri2_79*h(:,OP_1))
#endif
#endif
     end if

  v1ubb = temp
end function v1ubb

#ifdef USE3D
! V1upsif
! =====
! function v1upsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1upsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         ! [u, psi]_R*R
!         temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)
!         if(itor.eq.1) then
!            temp79a = temp79a - ri_79* &
!                   (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
!         end if 
!         ! [u, psi]_Z*R
!         temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ) 
!         ! (u, f')_R
!         temp79c = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
!         ! (u, f')_Z
!         temp79d = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
!         if(itor.eq.1) then
!            ! (u, f')
!            temp79e = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ)           
!         endif

!         do j=1, dofs_per_element
!            ! (nu, f')_R
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) &
!               + e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) 
!            ! (nu, f')_Z
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) &
!               + e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) 
!            ! [nu, psi]_R*R
!            tempc(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DRR) - e(j,:,OP_DR)*g(:,OP_DRZ) &
!               + e(j,:,OP_DRZ)*g(:,OP_DR) - e(j,:,OP_DRR)*g(:,OP_DZ)
!            if(itor.eq.1) then
!               tempc(j,:) = tempc(j,:) - ri_79* &
!                    (e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)) 
!            end if 
!            ! [nu, psi]_Z*R
!            tempd(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DRZ) - e(j,:,OP_DR)*g(:,OP_DZZ) &
!               + e(j,:,OP_DZZ)*g(:,OP_DR) - e(j,:,OP_DRZ)*g(:,OP_DZ) 
!            ! [nu, R^2*(U, f')]/R
!            tempe(j,:) = e(j,:,OP_DZ)*temp79c - e(j,:,OP_DR)*temp79d
!            if(itor.eq.1) then
!               tempe(j,:) = tempe(j,:) + 2*e(j,:,OP_DZ)*temp79e*ri_79 
!            endif
!         end do

!         temp = intx3(tempa,temp79a,r_79) &
!              + intx3(tempb,temp79b,r_79) &
!              + intx3(tempc,temp79c,r_79) &
!              + intx3(tempd,temp79d,r_79) &
!              - intx3(tempe,g(:,OP_GS),r_79) 
!      end if

!   v1upsif = temp
! end function v1upsif

function v1upsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1upsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp, tempe
  type(muarray) :: tempa, tempb, tempc, tempd
  type(muarray) :: tempa1, tempb1, tempc1, tempd1, tempe1
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
        ! [u, psi]_R*R
        tempa1 = mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &           
               + mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR)
        if(itor.eq.1) then
           tempa1 = tempa1 + &
                    mu(-ri_79*g(:,OP_DR),OP_DZ) + mu(ri_79*g(:,OP_DZ),OP_DR) 
        end if 
        ! [u, psi]_Z*R
        tempb1 = mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &           
               + mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR)
        ! (u, f')_R
        tempc1 = mu(h(:,OP_DR),OP_DRR) + mu(h(:,OP_DZ),OP_DRZ) &           
               + mu(h(:,OP_DRR),OP_DR) + mu(h(:,OP_DRZ),OP_DZ)            
        ! (u, f')_Z
        tempd1 = mu(h(:,OP_DR),OP_DRZ) + mu(h(:,OP_DZ),OP_DZZ) &           
               + mu(h(:,OP_DRZ),OP_DR) + mu(h(:,OP_DZZ),OP_DZ)            
        if(itor.eq.1) then
           ! (u, f')
           tempe1 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ)           
        endif

           ! (nu, f')_R
           tempa = &
                mu(h(:,OP_DRZ),OP_DZ) + mu(h(:,OP_DRR),OP_DR) &
              + mu(h(:,OP_DZ),OP_DRZ) + mu(h(:,OP_DR),OP_DRR)
           ! (nu, f')_Z
           tempb = &
                mu(h(:,OP_DZZ),OP_DZ) + mu(h(:,OP_DRZ),OP_DR) &
              + mu(h(:,OP_DZ),OP_DZZ) + mu(h(:,OP_DR),OP_DRZ)
           ! [nu, psi]_R*R
           tempc = &
                mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR) &
              + mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR)
           if(itor.eq.1) then
              tempc = tempc + &
                   mu(-ri_79*g(:,OP_DR),OP_DZ) + mu(ri_79*g(:,OP_DZ),OP_DR) 
           end if 
           ! [nu, psi]_Z*R
           tempd = &
                mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR) &
              + mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) 
           ! [nu, R^2*(U, f')]/R
           temp79a=1.
           tempe = prod(mu(temp79a,OP_DZ),tempc1) + prod(mu(-temp79a,OP_DR),tempd1)
           if(itor.eq.1) then
              tempe = tempe + prod(mu(2.*ri_79,OP_DZ),tempe1)
           endif

        temp = prod(tempa*r_79,tempa1) &
             + prod(tempb*r_79,tempb1) &
             + prod(tempc*r_79,tempc1) &
             + prod(tempd*r_79,tempd1) &
             + tempe*(-g(:,OP_GS)*r_79) 
     end if

  v1upsif = temp
end function v1upsif

! V1ubf
! =====
!function v1ubf(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v1ubf
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!  integer :: j


!     if(surface_int) then
!        temp = 0.
!     else
!        do j=1, dofs_per_element
!           ! (nu, f')'
!           tempa(j,:) = &
!                e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) &
!              + e(j,:,OP_DZP)*h(:,OP_DZ) + e(j,:,OP_DRP)*h(:,OP_DR) 
!           ! [nu,f'']*R
!           tempb(j,:) = &
!                e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
!        end do
!        ! F*u_LP + R^2*(u, F/R^2)
!        temp79a = g(:,OP_1)*f(:,OP_LP)  & 
!                + g(:,OP_DR)*f(:,OP_DR) + g(:,OP_DZ)*f(:,OP_DZ)
!        if(itor.eq.1) then
!           temp79a = temp79a - 2*ri_79*g(:,OP_1)*f(:,OP_DR)
!        end if 
!        ! (R^2(u, f'))_R/R^2
!        temp79b = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
!        if(itor.eq.1) then
!           temp79b = temp79b + 2*ri_79* &
!                 (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
!        end if 
!        ! (R^2(u, f'))_Z/R^2
!        temp79c = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
!        ![u,F]*R
!        temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ) 

!        temp = intx2(tempa,temp79a) &
!             + intx3(e(:,:,OP_DRP),temp79b,g(:,OP_1)) &
!             + intx3(e(:,:,OP_DZP),temp79c,g(:,OP_1)) &
!             + intx3(e(:,:,OP_DR),temp79b,g(:,OP_DP)) &
!             + intx3(e(:,:,OP_DZ),temp79c,g(:,OP_DP)) &
!             + intx2(tempb,temp79d)
!     end if

!  v1ubf = temp
!end function v1ubf

function v1ubf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1ubf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2, tempc, tempd
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f')'
           tempa = &
                mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR) &
              + mu(h(:,OP_DZ),OP_DZP) + mu(h(:,OP_DR),OP_DRP) 
           ! [nu,f'']*R
           tempb = &
                mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR) 
        ! F*u_LP + R^2*(u, F/R^2)
        tempa2 = mu(g(:,OP_1),OP_LP)  & 
              + mu(g(:,OP_DR),OP_DR) + mu(g(:,OP_DZ),OP_DZ)
        if(itor.eq.1) then
           tempa2 = tempa2 + mu(-2*ri_79*g(:,OP_1),OP_DR)
        end if 
        ! (R^2(u, f'))_R/R^2
        tempb2 = mu(h(:,OP_DR),OP_DRR) + mu(h(:,OP_DZ),OP_DRZ) &           
              + mu(h(:,OP_DRR),OP_DR) + mu(h(:,OP_DRZ),OP_DZ)            
        if(itor.eq.1) then
           tempb2 = tempb2 + (mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ)) &
                 * (2*ri_79)
        end if 
        ! (R^2(u, f'))_Z/R^2
        tempc = mu(h(:,OP_DR),OP_DRZ) + mu(h(:,OP_DZ),OP_DZZ) &           
              + mu(h(:,OP_DRZ),OP_DR) + mu(h(:,OP_DZZ),OP_DZ)            
        ![u,F]*R
        tempd = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR) 

        temp = prod(tempa,tempa2) &
             + prod(mu(g(:,OP_1),OP_DRP),tempb2) &
             + prod(mu(g(:,OP_1),OP_DZP),tempc) &
             + prod(mu(g(:,OP_DP),OP_DR),tempb2) &
             + prod(mu(g(:,OP_DP),OP_DZ),tempc) &
             + prod(tempb,tempd)
     end if

  v1ubf = temp
end function v1ubf

! V1uff
! =====
!function v1uff(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v1uff
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!  integer :: j


!     if(surface_int) then
!        temp = 0.
!     else
!        do j=1, dofs_per_element
!           !(nu, f')_R
!           tempa(j,:) = &
!                e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) &
!              + e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) 
!           !(nu, f')_Z
!           tempb(j,:) = &
!                e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) &
!              + e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) 
!        end do
!        !(u, f')_R
!        temp79a = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &           
!                + f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ)            
!        !(u, f')_Z
!        temp79b = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &           
!                + f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ)            

!        temp = - intx3(tempa,temp79a,r2_79) &
!               - intx3(tempb,temp79b,r2_79) 
!     end if

!  v1uff = temp
!end function v1uff

function v1uff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1uff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa,tempb, tempa2, tempb2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           !(nu, f')_R
           tempa = &
                mu(h(:,OP_DRZ),OP_DZ) + mu(h(:,OP_DRR),OP_DR) &
              + mu(h(:,OP_DZ),OP_DRZ) + mu(h(:,OP_DR),OP_DRR) 
           !(nu, f')_Z
           tempb = &
                mu(h(:,OP_DZZ),OP_DZ) + mu(h(:,OP_DRZ),OP_DR) &
              + mu(h(:,OP_DZ),OP_DZZ) + mu(h(:,OP_DR),OP_DRZ) 
        !(u, f')_R
        tempa2 = mu(g(:,OP_DR),OP_DRR) + mu(g(:,OP_DZ),OP_DRZ) &           
                + mu(g(:,OP_DRR),OP_DR) + mu(g(:,OP_DRZ),OP_DZ)            
        !(u, f')_Z
        tempb2 = mu(g(:,OP_DR),OP_DRZ) + mu(g(:,OP_DZ),OP_DZZ) &           
               + mu(g(:,OP_DRZ),OP_DR) + mu(g(:,OP_DZZ),OP_DZ)            

        temp = prod(tempa,tempa2*(-r2_79)) &
             + prod(tempb,tempb2*(-r2_79)) 
     end if

  v1uff = temp
end function v1uff
#endif

! V1up
! ====
!function v1up(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1up
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

     !if(itor.eq.0) then
        !temp = 0.
     !else
        !if(surface_int) then
           !temp = 2.* &
                !(intx5(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_DR),norm79(:,2)) &
                !-intx5(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_DZ),norm79(:,2))) &
                !+4.*gam*intx4(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_1),norm79(:,2))
        !else
           !temp = 2.* &
                !(intx4(e(:,:,OP_DZ),r_79,f(:,OP_DR),g(:,OP_DZ)) &
                !-intx4(e(:,:,OP_DZ),r_79,f(:,OP_DZ),g(:,OP_DR))) &
                !-4.*gam*intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
        !endif
     !end if


  !v1up = temp
!end function v1up

function v1up(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1up
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(itor.eq.0) then
        temp%len = 0
     else
        if(surface_int) then
           temp = (prod(2.*r_79*g(:,OP_DR)*norm79(:,2),OP_1,OP_DZ) &
                + prod(-2.*r_79*g(:,OP_DZ)*norm79(:,2),OP_1,OP_DR)) &
                +prod(4.*gam*g(:,OP_1)*norm79(:,2),OP_1,OP_DZ)
         else
           temp = ( prod(2.*r_79*g(:,OP_DZ),OP_DZ,OP_DR) &
                +  prod(-2.*r_79*g(:,OP_DR),OP_DZ,OP_DZ)) &
                +  prod(-4.*gam*g(:,OP_1),OP_DZ,OP_DZ)
        endif
     end if

  v1up = temp
end function v1up



! V1vpsipsi
! =========
! function v1vpsipsi(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1vpsipsi
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

! #if defined(USE3D) || defined(USECOMPLEX)
!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempc
!   integer :: j

!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP)
!            tempc(j,:) = &
!                 f(:,OP_DP)* &
!                 (e(j,:,OP_DZ)*h(:,OP_DR ) - e(j,:,OP_DR)*h(:,OP_DZ )) &
!                 +  f(:,OP_1 )* &
!                 (e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP))
!         end do
! #ifdef USEST
!         temp79b = f(:,OP_1)*g(:,OP_GS) &
!              + f(:,OP_DZ )*g(:,OP_DZ) + f(:,OP_DR )*g(:,OP_DR)
! #else
!         temp79b = f(:,OP_DP)*g(:,OP_GS) + f(:,OP_1)*g(:,OP_GSP) &
!              + f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
!              + f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)
! #endif
!         temp = intx4(tempa,ri_79,g(:,OP_DZ),f(:,OP_DR)) &
!              - intx4(tempa,ri_79,g(:,OP_DR),f(:,OP_DZ)) &
! #ifdef USEST
!              - intx4(e(:,:,OP_DZP),ri_79,h(:,OP_DR),temp79b) &
!              + intx4(e(:,:,OP_DRP),ri_79,h(:,OP_DZ),temp79b) &
!              - intx4(e(:,:,OP_DZ),ri_79,h(:,OP_DRP),temp79b) &
!              + intx4(e(:,:,OP_DR),ri_79,h(:,OP_DZP),temp79b) &
! #else
!              + intx4(e(:,:,OP_DZ),ri_79,h(:,OP_DR),temp79b) &
!              - intx4(e(:,:,OP_DR),ri_79,h(:,OP_DZ),temp79b) &
! #endif
!              + intx3(tempc,ri_79,g(:,OP_GS))
!         temp= -temp
!      end if

!   v1vpsipsi = temp
! #else
!   v1vpsipsi = 0.
! #endif
! end function v1vpsipsi

function v1vpsipsi(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vpsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp, tempc
  type(muarray) :: tempa, tempb

     if(surface_int) then
        temp%len = 0
     else
           tempa = mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR)
           tempc = &
                prod(h(:,OP_DR ),OP_DZ,OP_DP) + prod(-h(:,OP_DZ ),OP_DR,OP_DP) &
                +  &
                prod(h(:,OP_DRP),OP_DZ,OP_1) + prod(-h(:,OP_DZP),OP_DR,OP_1)
#ifdef USEST
        tempb = mu(g(:,OP_GS),OP_1) &
              + mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
#else
        tempb = mu(g(:,OP_GS),OP_DP) + mu(g(:,OP_GSP),OP_1) &
             +  mu(g(:,OP_DZ),OP_DZP) +mu(g(:,OP_DR ),OP_DRP) &
             +  mu(g(:,OP_DZP),OP_DZ) +mu(g(:,OP_DRP),OP_DR)
#endif
        temp = prod(tempa,mu(-ri_79*g(:,OP_DZ),OP_DR)) &
             + prod(tempa,mu(ri_79*g(:,OP_DR),OP_DZ)) &
#ifdef USEST
             + prod(mu( ri_79*h(:,OP_DR),OP_DZP),tempb) &
             + prod(mu(-ri_79*h(:,OP_DZ),OP_DRP),tempb) &
             + prod(mu( ri_79*h(:,OP_DRP),OP_DZ),tempb) &
             + prod(mu(-ri_79*h(:,OP_DZP),OP_DR),tempb) &
#else
             + prod(mu(-ri_79*h(:,OP_DR),OP_DZ),tempb) &
             + prod(mu( ri_79*h(:,OP_DZ),OP_DR),tempb) &
#endif
             + tempc*(-ri_79*g(:,OP_GS))
     end if

  v1vpsipsi = temp
#else
  v1vpsipsi%len = 0
#endif
end function v1vpsipsi


! V1vpsib
! =======
!function v1vpsib(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v1vpsib
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa
!  integer :: j

!  temp = 0.

!     if(surface_int) then
!        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
!        temp = intx4(e(:,:,OP_1),temp79a,norm79(:,2),h(:,OP_DR)) &
!             - intx4(e(:,:,OP_1),temp79a,norm79(:,1),h(:,OP_DZ))
!#if defined(USE3D) || defined(USECOMPLEX)
!        temp79a = ri2_79*h(:,OP_DP)
!        temp = temp &
!             - intx5(e(:,:,OP_1),temp79a,norm79(:,1),g(:,OP_DR),f(:,OP_DP)) &
!             - intx5(e(:,:,OP_1),temp79a,norm79(:,2),g(:,OP_DZ),f(:,OP_DP)) &
!             - intx5(e(:,:,OP_1),temp79a,norm79(:,1),g(:,OP_DRP),f(:,OP_1)) &
!             - intx5(e(:,:,OP_1),temp79a,norm79(:,2),g(:,OP_DZP),f(:,OP_1))
!#endif
!     else
!        temp = 0.

!#if defined(USE3D) || defined(USECOMPLEX)
!#ifdef USEST 
!        !temp = - intx5(e(:,:,OP_DR),f(:,OP_1),g(:,OP_DR),ri2_79,h(:,OP_DPP))&
!        !       - intx5(e(:,:,OP_DZ),f(:,OP_1),g(:,OP_DZ),ri2_79,h(:,OP_DPP))
!        do j=1, dofs_per_element
!           tempa(j,:) = &
!                +(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
!                *(g(:,OP_DZP)*f(:,OP_1) + g(:,OP_DZ)*f(:,OP_DP)) &
!                +(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
!                *(g(:,OP_DRP)*f(:,OP_1) + g(:,OP_DR)*f(:,OP_DP))
!        end do
!        temp = intx2(tempa,ri2_79)
!#else
!        do j=1, dofs_per_element
!           tempa(j,:) = &
!                f(:,OP_DPP)* &
!                (e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR)) &
!                +2.*f(:,OP_DP)* &
!                (e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP)) &
!                +   f(:,OP_1)* &
!                (e(j,:,OP_DZ)*g(:,OP_DZPP) + e(j,:,OP_DR)*g(:,OP_DRPP))
!        end do
!        temp = -intx3(tempa,ri2_79,h(:,OP_1))
!#endif
!#endif
!     end if

!  v1vpsib = temp
!end function v1vpsib

function v1vpsib(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vpsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp, tempb
  type(muarray) :: tempa

  temp%len = 0
     if(surface_int) then
        tempa = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
        temp = prod(mu( norm79(:,2)*h(:,OP_DR),OP_1),tempa) &
             + prod(mu(-norm79(:,1)*h(:,OP_DZ),OP_1),tempa)
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = ri2_79*h(:,OP_DP)
        temp = temp &
             + prod(-temp79a*norm79(:,1)*g(:,OP_DR),OP_1,OP_DP) &
             + prod(-temp79a*norm79(:,2)*g(:,OP_DZ),OP_1,OP_DP) &
             + prod(-temp79a*norm79(:,1)*g(:,OP_DRP),OP_1,OP_1) &
             + prod(-temp79a*norm79(:,2)*g(:,OP_DZP),OP_1,OP_1)
#endif
     else
        temp%len = 0

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST 
        !temp = - intx5(e(:,:,OP_DR),f(:,OP_1),g(:,OP_DR),ri2_79,h(:,OP_DPP))&
        !       - intx5(e(:,:,OP_DZ),f(:,OP_1),g(:,OP_DZ),ri2_79,h(:,OP_DPP))
           tempb = &
                  prod(mu(h(:,OP_DP),OP_DZ) + mu(h(:,OP_1),OP_DZP), &
                  mu(g(:,OP_DZP),OP_1) + mu(g(:,OP_DZ),OP_DP)) &
                 +prod(mu(h(:,OP_DP),OP_DR) + mu(h(:,OP_1),OP_DRP), &
                  mu(g(:,OP_DRP),OP_1) + mu(g(:,OP_DR),OP_DP))
        temp = tempb*(ri2_79)
#else
        tempb = &
                prod(g(:,OP_DZ),OP_DZ,OP_DPP) + prod(g(:,OP_DR),OP_DR,OP_DPP) &
              + prod(2.*g(:,OP_DZP),OP_DZ,OP_DP) + prod(2.*g(:,OP_DRP),OP_DR,OP_DP) &
              + prod(g(:,OP_DZPP),OP_DZ,OP_1) + prod(g(:,OP_DRPP),OP_DR,OP_1)
        temp = tempb*(-ri2_79*h(:,OP_1))
#endif
#endif
     end if

  v1vpsib = temp
end function v1vpsib

#ifdef USE3D
! V1vpsif
! =====
! function v1vpsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1vpsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe, tempf
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! [nu, psi]*R
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
!            ! [nu, f'']*R
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
!            ! (nu, psi')
!            tempc(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP) 
!            ! (nu, f')
!            tempd(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
!            ! (nu, f'')
!            tempe(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
!            ! (nu', f')
!            tempf(j,:) = &
!                 e(j,:,OP_DZP)*h(:,OP_DZ) + e(j,:,OP_DRP)*h(:,OP_DR) 
!         end do

!         ! [v, f']'*R
!         temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
!         ! [v, psi]*R
!         temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)            
!         ! (v, f') + v*f'_LP
!         temp79c = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
!                 + f(:,OP_1)*h(:,OP_LP)
!         ! (v, psi)
!         temp79d = f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ) 

!         temp = - intx2(tempa,temp79a) &
!                + intx2(tempb,temp79b) &
!                - intx2(tempc,temp79c) &
!                - intx2(tempe,temp79d) & 
!                - intx2(tempf,temp79d) & 
!                + intx3(tempd,g(:,OP_GS),f(:,OP_DP)) &
!                - intx3(tempf,g(:,OP_GS),f(:,OP_1)) 
!      end if

!   v1vpsif = temp
! end function v1vpsif

function v1vpsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vpsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd, tempe, tempf
  type(muarray) :: tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! [nu, psi]*R
           tempa = &
                mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR) 
           ! [nu, f'']*R
           tempb = &
                mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR) 
           ! (nu, psi')
           tempc = &
                mu(g(:,OP_DZP),OP_DZ) + mu(g(:,OP_DRP),OP_DR) 
           ! (nu, f')
           tempd = &
                mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR) 
           ! (nu, f'')
           tempe = &
                mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR) 
           ! (nu', f')
           tempf = &
                mu(h(:,OP_DZ),OP_DZP) + mu(h(:,OP_DR),OP_DRP) 

        ! [v, f']'*R
        tempa2 = mu(h(:,OP_DR),OP_DZP) + mu(-h(:,OP_DZ),OP_DRP) &           
               + mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR)            
        ! [v, psi]*R
        tempb2 = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)            
        ! (v, f') + v*f'_LP
        tempc2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) & 
                + mu(h(:,OP_LP),OP_1)
        ! (v, psi)
        tempd2 = mu(g(:,OP_DR),OP_DR) + mu(g(:,OP_DZ),OP_DZ) 

        temp = prod((-1.)*tempa,tempa2) &
             + prod(tempb,tempb2) &
             + prod((-1.)*tempc,tempc2) &
             + prod((-1.)*tempe,tempd2) & 
             + prod((-1.)*tempf,tempd2) & 
             + prod(tempd,mu( g(:,OP_GS),OP_DP)) &
             + prod(tempf,mu(-g(:,OP_GS),OP_1)) 
     end if

  v1vpsif = temp
end function v1vpsif

! V1vbf
! =====
! function v1vbf(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1vbf
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! [nu, f']*R
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ) 
!            ! [nu, f'']*R
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
!            ! [nu', f'']*R
!            tempc(j,:) = &
!                 e(j,:,OP_DZP)*h(:,OP_DRP) - e(j,:,OP_DRP)*h(:,OP_DZP) 
!         end do

!         temp = - intx4(tempa,g(:,OP_1),f(:,OP_DPP),ri_79) &
!                - intx4(tempb,g(:,OP_1),f(:,OP_DP),ri_79) &
!                + intx4(tempc,g(:,OP_1),f(:,OP_1),ri_79) 
!      end if

!   v1vbf = temp
! end function v1vbf

function v1vbf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vbf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! [nu, f']*R
           tempa = &
                mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR) 
            ! [nu, f'']*R
           tempb = &
                mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR) 
           ! [nu', f'']*R
           tempc = &
                mu(h(:,OP_DRP),OP_DZP) + mu(-h(:,OP_DZP),OP_DRP) 

        temp = prod(tempa,mu(-g(:,OP_1)*ri_79,OP_DPP)) &
             + prod(tempb,mu(-g(:,OP_1)*ri_79,OP_DP)) &
             + prod(tempc,mu( g(:,OP_1)*ri_79,OP_1)) 
     end if

  v1vbf = temp
end function v1vbf

! V1vff
! =====
! function v1vff(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1vff
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (nu, f')
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR) 
!            ! [nu, f'']*R
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DRP) - e(j,:,OP_DR)*g(:,OP_DZP) 
!         end do
!         ! [v, f']'*R 
!         temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
!         ! (v, f')
!         temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) 

!         temp = intx3(tempa,temp79a,r_79) &
!              - intx3(tempb,temp79b,r_79) 
!      end if

!   v1vff = temp
! end function v1vff

function v1vff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb
  type(muarray) :: tempa2, tempb2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f')
           tempa = &
                mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR) 
           ! [nu, f'']*R
           tempb = &
                mu(g(:,OP_DRP),OP_DZ) + mu(-g(:,OP_DZP),OP_DR) 
        ! [v, f']'*R 
        tempa2 = mu(h(:,OP_DR),OP_DZP) + mu(-h(:,OP_DZ),OP_DRP) &           
               + mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR)            
        ! (v, f')
        tempb2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) 

        temp = prod(tempa,tempa2*r_79) &
             + prod(tempb,tempb2*(-r_79)) 
     end if

  v1vff = temp
end function v1vff
#endif

! V1vp
! ====
!function v1vp(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1vp
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp


     !temp = 0.
!#if defined(USE3D) || defined(USECOMPLEX)
     !if(itor.eq.1) then
        !if(surface_int) then
           !temp = -2.* &
                !(    intx4(e(:,:,OP_1),f(:,OP_1),g(:,OP_DP),norm79(:,2)) &
                !+gam*intx4(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1),norm79(:,2)))
        !else
           !temp79a =  f(:,OP_1 )*g(:,OP_DP) &
                !+ gam*f(:,OP_DP)*g(:,OP_1 )
           !temp = 2.*intx2(e(:,:,OP_DZ),temp79a)
        !end if
     !end if
!#endif

  !v1vp = temp
!end function v1vp

function v1vp(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1vp
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp
  type(muarray) :: tempa

     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     if(itor.eq.1) then
        if(surface_int) then
           temp = &
                ( prod(-2.*g(:,OP_DP)*norm79(:,2),OP_1,OP_1) &
                + prod(-2.*gam*g(:,OP_1)*norm79(:,2),OP_1,OP_DP))
         else
           tempa = mu(g(:,OP_DP),OP_1) &
                + mu(gam*g(:,OP_1 ),OP_DP)
           temp79a = 1.
           temp = prod(mu(2.*temp79a,OP_DZ),tempa)
         end if
     end if
#endif

  v1vp = temp
end function v1vp


! V1chipsipsi
! ===========
!function v1chipsipsi(e,f,g,h)
  !use basic
  !use arrays
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1chipsipsi
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  !vectype, dimension(dofs_per_element) :: temp
  !vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempd, tempe, tempf
  !integer :: j

     !if(surface_int) then
        !do j=1, dofs_per_element
           !tempa(j,:) = e(j,:,OP_DR)*h(:,OP_DZ) - e(j,:,OP_DZ)*h(:,OP_DR)
        !end do
        !temp = intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
             !+ intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
             !+ intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
             !+ intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
             !+ intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
             !+ intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
             !+ intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
             !+ intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
        !if(itor.eq.1) then
           !temp = temp - 2.* &
                !(intx5(tempa,ri4_79,norm79(:,1),f(:,OP_DZ),g(:,OP_DZ)) &
                !+intx5(tempa,ri4_79,norm79(:,1),f(:,OP_DR),g(:,OP_DR)))
        !endif
     !else
        !temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
             !+    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
        !temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
             !+    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)
        
        !do j=1, dofs_per_element
           !tempe(j,:) = e(j,:,OP_DRZ)*h(:,OP_DR ) - e(j,:,OP_DRR)*h(:,OP_DZ ) &
                !+    e(j,:,OP_DZ )*h(:,OP_DRR) - e(j,:,OP_DR )*h(:,OP_DRZ)
           !tempf(j,:) = e(j,:,OP_DZZ)*h(:,OP_DR ) - e(j,:,OP_DRZ)*h(:,OP_DZ ) &
                !+    e(j,:,OP_DZ )*h(:,OP_DRZ) - e(j,:,OP_DR )*h(:,OP_DZZ)
        !end do

        !temp = intx3(tempe,ri3_79,temp79b) &
             !+ intx3(tempf,ri3_79,temp79c) &
             !+ intx4(e(:,:,OP_DR),ri3_79,temp79c,h(:,OP_GS)) &
             !- intx4(e(:,:,OP_DZ),ri3_79,temp79b,h(:,OP_GS))
     
        !if(itor.eq.1) then
           !temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
           !do j=1, dofs_per_element
              !tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
           !end do
           
           !temp = temp + &
                !(2.*intx4(e(:,:,OP_DZ),ri4_79,temp79a,h(:,OP_GS)) &
                !-2.*intx3(tempe,ri4_79,temp79a) &
                !+   intx3(tempd,ri4_79,temp79b) &
                !-2.*intx3(tempd,ri5_79,temp79a))
        !endif
     !end if


  !v1chipsipsi = temp
!end function v1chipsipsi

function v1chipsipsi(g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chipsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) ::  g, h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd, tempe, tempf

     if(surface_int) then
        tempa = mu(h(:,OP_DZ),OP_DR) + mu(-h(:,OP_DR),OP_DZ)
        temp = prod(tempa,mu(ri3_79*norm79(:,1)*g(:,OP_DZ ),OP_DRZ)) &
             + prod(tempa,mu(ri3_79*norm79(:,1)*g(:,OP_DR ),OP_DRR)) &
             + prod(tempa,mu(ri3_79*norm79(:,1)*g(:,OP_DRZ),OP_DZ)) &
             + prod(tempa,mu(ri3_79*norm79(:,1)*g(:,OP_DRR),OP_DR)) &
             + prod(tempa,mu(ri3_79*norm79(:,2)*g(:,OP_DZ ),OP_DZZ)) &
             + prod(tempa,mu(ri3_79*norm79(:,2)*g(:,OP_DR ),OP_DRZ)) &
             + prod(tempa,mu(ri3_79*norm79(:,2)*g(:,OP_DZZ),OP_DZ)) &
             + prod(tempa,mu(ri3_79*norm79(:,2)*g(:,OP_DRZ),OP_DR))
        if(itor.eq.1) then
           temp = temp + &
                 prod(tempa,mu(-2.*ri4_79*norm79(:,1)*g(:,OP_DZ),OP_DZ)) &
                +prod(tempa,mu(-2.*ri4_79*norm79(:,1)*g(:,OP_DR),OP_DR))
         endif
     else
        tempb = mu(g(:,OP_DZ ),OP_DRZ) + mu(g(:,OP_DR ),OP_DRR) &
             +  mu(g(:,OP_DRZ),OP_DZ)  + mu(g(:,OP_DRR),OP_DR)
        tempc = mu(g(:,OP_DZ ),OP_DZZ) + mu(g(:,OP_DR ),OP_DRZ) &
             +  mu(g(:,OP_DZZ),OP_DZ)  + mu(g(:,OP_DRZ),OP_DR)
         
        tempe = mu(h(:,OP_DR ),OP_DRZ) + mu(-h(:,OP_DZ ),OP_DRR) &
              + mu(h(:,OP_DRR),OP_DZ)  + mu(-h(:,OP_DRZ),OP_DR)
        tempf = mu(h(:,OP_DR ),OP_DZZ) + mu(-h(:,OP_DZ ),OP_DRZ) &
              + mu(h(:,OP_DRZ),OP_DZ)  + mu(-h(:,OP_DZZ),OP_DR)
        
        temp = prod(tempe*ri3_79,tempb) &
             + prod(tempf*ri3_79,tempc) &
             + prod(mu( ri3_79*h(:,OP_GS),OP_DR),tempc) &
             + prod(mu(-ri3_79*h(:,OP_GS),OP_DZ),tempb)
     
        if(itor.eq.1) then
           tempa = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
           tempd = mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR)
           
           temp = temp + &
                (prod(mu(2.*ri4_79*h(:,OP_GS),OP_DZ),tempa) &
                +prod(tempe*(-2.*ri4_79),tempa) &
                +prod(tempd*ri4_79,tempb) &
                +prod(tempd*(-2.*ri5_79),tempa))
         endif
     end if

  v1chipsipsi = temp
end function v1chipsipsi

! V1chipsib
! =========
!function v1chipsib(e,f,g,h) 
  !use basic
  !use arrays
  !use m3dc1_nint

  !implicit none
  
  !vectype, dimension(dofs_per_element) :: v1chipsib
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

!#if defined(USE3D) || defined(USECOMPLEX)
  !vectype, dimension(dofs_per_element) :: temp
  !vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
  !integer :: j

     !if(surface_int) then
        !temp79a = f(:,OP_DRP)*g(:,OP_DR ) + f(:,OP_DZP)*g(:,OP_DZ ) &
             !+    f(:,OP_DR )*g(:,OP_DRP) + f(:,OP_DZ )*g(:,OP_DZP)
        !temp79b = ri4_79*h(:,OP_DP)
        !temp = &
             !- intx5(e(:,:,OP_DR),ri4_79,norm79(:,1),temp79a,h(:,OP_1)) &
             !- intx5(e(:,:,OP_DZ),ri4_79,norm79(:,2),temp79a,h(:,OP_1)) &
             !- intx5(e(:,:,OP_1),norm79(:,1),temp79b,f(:,OP_DRR),g(:,OP_DR )) &
             !- intx5(e(:,:,OP_1),norm79(:,1),temp79b,f(:,OP_DRZ),g(:,OP_DZ )) &
             !- intx5(e(:,:,OP_1),norm79(:,1),temp79b,f(:,OP_DR ),g(:,OP_DRR)) &
             !- intx5(e(:,:,OP_1),norm79(:,1),temp79b,f(:,OP_DZ ),g(:,OP_DRZ)) &
             !- intx5(e(:,:,OP_1),norm79(:,2),temp79b,f(:,OP_DRZ),g(:,OP_DR )) &
             !- intx5(e(:,:,OP_1),norm79(:,2),temp79b,f(:,OP_DZZ),g(:,OP_DZ )) &
             !- intx5(e(:,:,OP_1),norm79(:,2),temp79b,f(:,OP_DR ),g(:,OP_DRZ)) &
             !- intx5(e(:,:,OP_1),norm79(:,2),temp79b,f(:,OP_DZ ),g(:,OP_DZZ))
        !if(itor.eq.1) then
           !temp79b = ri5_79*h(:,OP_DP)

           !temp = temp + 2.* &
                !(intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ),g(:,OP_DZ)) &
                !+intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR),g(:,OP_DR)))
        !endif
     !else
        !temp79a = h(:,OP_DZP)*f(:,OP_DR ) - h(:,OP_DRP)*f(:,OP_DZ ) &
             !+    h(:,OP_DZ )*f(:,OP_DRP) - h(:,OP_DR )*f(:,OP_DZP)
        !do j=1, dofs_per_element
           !tempb(j,:) = (e(j,:,OP_DZ)*f(:,OP_DZ )+e(j,:,OP_DR)*f(:,OP_DR )) &
                !*h(:,OP_DP) &
                !+    (e(j,:,OP_DZ)*f(:,OP_DZP)+e(j,:,OP_DR)*f(:,OP_DRP)) &
                !*h(:,OP_1 )
        !end do
        !temp79c = f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
             !+    f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)
        !temp79d = h(:,OP_1)*f(:,OP_GS) &
             !+ h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
        !if(itor.eq.1) then
           !temp79a = temp79a + 4.*ri_79* &
                !(h(:,OP_DP)*f(:,OP_DZ) + h(:,OP_1)*f(:,OP_DZP))
           !temp79d = temp79d - 2.*ri_79*h(:,OP_1)*f(:,OP_DR)
        !endif
     
        !temp = intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DR),temp79a)  &
             !- intx4(e(:,:,OP_DR),ri4_79,g(:,OP_DZ),temp79a)  &
             !- intx3(tempb,ri4_79,g(:,OP_GS))               &
             !- intx4(e(:,:,OP_GS),ri4_79,h(:,OP_1 ),temp79c)  &
             !- intx4(e(:,:,OP_DZ),ri4_79,h(:,OP_DZ),temp79c)  &
             !- intx4(e(:,:,OP_DR),ri4_79,h(:,OP_DR),temp79c)  &
             !+ intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),temp79d) &
             !+ intx4(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),temp79d)
        !temp = -temp
     !end if

  !v1chipsib = temp
!#else
  !v1chipsib = 0.
!#endif

!end function v1chipsib

function v1chipsib(g,h) 
  use basic
  use arrays
  use m3dc1_nint

  implicit none
  
  type(prodarray) :: v1chipsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp, tempb
  type(muarray) :: tempa, tempc, tempd

     if(surface_int) then
        tempa = mu(g(:,OP_DR ),OP_DRP) + mu(g(:,OP_DZ ),OP_DZP) &
             +    mu(g(:,OP_DRP),OP_DR)  + mu(g(:,OP_DZP),OP_DZ)
        temp79b = ri4_79*h(:,OP_DP)
        temp = &
               prod(mu(-ri4_79*norm79(:,1)*h(:,OP_1),OP_DR),tempa) &
             + prod(mu(-ri4_79*norm79(:,2)*h(:,OP_1),OP_DZ),tempa) &
             + prod(-norm79(:,1)*temp79b*g(:,OP_DR),OP_1,OP_DRR) &
             + prod(-norm79(:,1)*temp79b*g(:,OP_DZ),OP_1,OP_DRZ) &
             + prod(-norm79(:,1)*temp79b*g(:,OP_DRR),OP_1,OP_DR ) &
             + prod(-norm79(:,1)*temp79b*g(:,OP_DRZ),OP_1,OP_DZ ) &
             + prod(-norm79(:,2)*temp79b*g(:,OP_DR),OP_1,OP_DRZ) &
             + prod(-norm79(:,2)*temp79b*g(:,OP_DZ),OP_1,OP_DZZ) &
             + prod(-norm79(:,2)*temp79b*g(:,OP_DRZ),OP_1,OP_DR ) &
             + prod(-norm79(:,2)*temp79b*g(:,OP_DZZ),OP_1,OP_DZ )
         if(itor.eq.1) then
           temp79b = ri5_79*h(:,OP_DP)

           temp = temp + &
                (prod(2.*temp79b*norm79(:,1)*g(:,OP_DZ),OP_1,OP_DZ) &
                +prod(2.*temp79b*norm79(:,1)*g(:,OP_DR),OP_1,OP_DR))
         endif
     else
        tempa = mu(h(:,OP_DZP),OP_DR) + mu(-h(:,OP_DRP),OP_DZ) &
             +    mu(h(:,OP_DZ ),OP_DRP)+ mu(-h(:,OP_DR ),OP_DZP)
          tempb = (prod(h(:,OP_DP),OP_DZ,OP_DZ )+prod(h(:,OP_DP),OP_DR,OP_DR )) &
                     + (prod(h(:,OP_1),OP_DZ,OP_DZP) +prod(h(:,OP_1),OP_DR,OP_DRP))
         tempc = mu(g(:,OP_DZ ),OP_DZP) + mu(g(:,OP_DR ),OP_DRP) &
              +  mu(g(:,OP_DZP),OP_DZ)  + mu(g(:,OP_DRP),OP_DR)
         tempd = mu(h(:,OP_1),OP_GS) &
               + mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
          if(itor.eq.1) then
           tempa = tempa + &
                (mu(h(:,OP_DP),OP_DZ) + mu(h(:,OP_1),OP_DZP))*(4.*ri_79)
            tempd = tempd + mu(-2.*ri_79*h(:,OP_1),OP_DR)
        endif
     
        temp = prod(mu(-ri4_79*g(:,OP_DR),OP_DZ),tempa)  &
             + prod(mu( ri4_79*g(:,OP_DZ),OP_DR),tempa)  &
             + tempb*(ri4_79*g(:,OP_GS))               &
             + prod(mu( ri4_79*h(:,OP_1 ),OP_GS),tempc)  &
             + prod(mu( ri4_79*h(:,OP_DZ),OP_DZ),tempc)  &
             + prod(mu( ri4_79*h(:,OP_DR),OP_DR),tempc)  &
             + prod(mu(-ri4_79*g(:,OP_DZP),OP_DZ),tempd) &
             + prod(mu(-ri4_79*g(:,OP_DRP),OP_DR),tempd)
     end if
  v1chipsib = temp
#else
  v1chipsib%len = 0
#endif

end function v1chipsib

! V1chibb
! =======
!function v1chibb(e,f,g,h)
!  use basic
!  use arrays
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v1chibb
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa
!  integer :: j

!     if(surface_int) then
!        temp79a = norm79(:,1)*h(:,OP_DZ) - norm79(:,2)*h(:,OP_DR)
!        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,f(:,OP_GS),g(:,OP_1 )) &
!             + intx5(e(:,:,OP_1),ri3_79,temp79a,f(:,OP_DZ),g(:,OP_DZ)) &
!             + intx5(e(:,:,OP_1),ri3_79,temp79a,f(:,OP_DR),g(:,OP_DR)) 
!        if(itor.eq.1) then
!           temp = temp &
!                - 2.*intx5(e(:,:,OP_1),ri4_79,temp79a,f(:,OP_DR),g(:,OP_1))
!        endif
!#if defined(USE3D) || defined(USECOMPLEX)
!        temp79b = ri5_79*h(:,OP_DP)
!        temp = temp &
!             + intx5(e(:,:,OP_1),temp79b,g(:,OP_1 ),norm79(:,1),f(:,OP_DZP)) &
!             - intx5(e(:,:,OP_1),temp79b,g(:,OP_1 ),norm79(:,2),f(:,OP_DRP)) &
!             + intx5(e(:,:,OP_1),temp79b,g(:,OP_DP),norm79(:,1),f(:,OP_DZ )) &
!             - intx5(e(:,:,OP_1),temp79b,g(:,OP_DP),norm79(:,2),f(:,OP_DR ))
!#endif
!     else
!        temp = 0.
!#if defined(USE3D) || defined(USECOMPLEX)
!#ifdef USEST 
!        !temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),g(:,OP_1),ri5_79,h(:,OP_DPP))&
!        !     - intx5(e(:,:,OP_DR),f(:,OP_DZ),g(:,OP_1),ri5_79,h(:,OP_DPP))
!        do j=1, dofs_per_element
!        tempa(j,:) = &
!                 (e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
!                *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP)) &
!                -(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
!                *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP))
!        end do
!        temp = intx2(tempa,ri5_79)
!#else
!        do j=1, dofs_per_element
!        tempa(j,:) = &
!             (e(j,:,OP_DZ)*f(:,OP_DR) - e(j,:,OP_DR)*f(:,OP_DZ))*g(:,OP_DPP) &
!        + 2.*(e(j,:,OP_DZ)*f(:,OP_DRP) - e(j,:,OP_DR)*f(:,OP_DZP))*g(:,OP_DP) &
!        +    (e(j,:,OP_DZ)*f(:,OP_DRPP) - e(j,:,OP_DR)*f(:,OP_DZPP))*g(:,OP_1)
!        end do
!        temp = -intx3(tempa,ri5_79,h(:,OP_1))
!#endif
!#endif
!     end if

!  v1chibb = temp
!end function v1chibb

function v1chibb(g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chibb
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h

  type(prodarray) :: temp, tempa

     if(surface_int) then
        temp79a = norm79(:,1)*h(:,OP_DZ) - norm79(:,2)*h(:,OP_DR)
        temp = prod(ri3_79*temp79a*g(:,OP_1),OP_1,OP_GS) &
             + prod(ri3_79*temp79a*g(:,OP_DZ),OP_1,OP_DZ) &
             + prod(ri3_79*temp79a*g(:,OP_DR),OP_1,OP_DR) 
        if(itor.eq.1) then
           temp = temp &
                + prod(-2.*ri4_79*temp79a*g(:,OP_1),OP_1,OP_DR)
        endif
#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = ri5_79*h(:,OP_DP)
        temp = temp &
             + prod( temp79b*g(:,OP_1 )*norm79(:,1),OP_1,OP_DZP) &
             + prod(-temp79b*g(:,OP_1 )*norm79(:,2),OP_1,OP_DRP) &
             + prod( temp79b*g(:,OP_DP)*norm79(:,1),OP_1,OP_DZ ) &
             + prod(-temp79b*g(:,OP_DP)*norm79(:,2),OP_1,OP_DR )
#endif
     else
        temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST 
        !temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),g(:,OP_1),ri5_79,h(:,OP_DPP))&
        !     - intx5(e(:,:,OP_DR),f(:,OP_DZ),g(:,OP_1),ri5_79,h(:,OP_DPP))
        tempa = &
                 prod(mu(h(:,OP_DP),OP_DZ) + mu(h(:,OP_1),OP_DZP), &
                 mu(g(:,OP_1),OP_DRP) + mu(g(:,OP_DP),OP_DR)) &
                +prod(mu(-h(:,OP_DP),OP_DR) + mu(-h(:,OP_1),OP_DRP), &
                 mu(g(:,OP_1),OP_DZP) + mu(g(:,OP_DP),OP_DZ))
        temp = tempa*ri5_79
#else
        tempa = &
             prod(g(:,OP_DPP),OP_DZ,OP_DR) + prod(-g(:,OP_DPP),OP_DR,OP_DZ) &
        +    prod(2.*g(:,OP_DP),OP_DZ,OP_DRP) + prod(-2.*g(:,OP_DP),OP_DR,OP_DZP) &
        +    prod(g(:,OP_1),OP_DZ,OP_DRPP) + prod(-g(:,OP_1),OP_DR,OP_DZPP)
         temp = tempa*(-ri5_79*h(:,OP_1))
#endif
#endif
     end if

  v1chibb = temp
end function v1chibb

#ifdef USE3D
! V1chipsif
! =====
! function v1chipsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1chipsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempc, tempd
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! ([nu, psi]*R^2)_R / R
!            tempa(j,:) = &
!                 e(j,:,OP_DRZ)*g(:,OP_DR) - e(j,:,OP_DRR)*g(:,OP_DZ) &
!               + e(j,:,OP_DZ)*g(:,OP_DRR) - e(j,:,OP_DR)*g(:,OP_DRZ)
!            if(itor.eq.1) then 
!               tempa(j,:) = tempa(j,:) + ri_79* &
!                   (e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)) 
!            end if
!            ! ([nu, psi]*R^2)_Z / R
!            tempb(j,:) = &
!                 e(j,:,OP_DZZ)*g(:,OP_DR) - e(j,:,OP_DRZ)*g(:,OP_DZ) &
!               + e(j,:,OP_DZ)*g(:,OP_DRZ) - e(j,:,OP_DR)*g(:,OP_DZZ)
!            ! (R^2*(nu, f'))_R / R^2
!            tempc(j,:) = &
!                 e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) & 
!               + e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) 
!            if(itor.eq.1) then 
!               tempc(j,:) = tempc(j,:) + 2*ri_79* &
!                   (e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)) 
!            end if
!            ! (R^2*(nu, f'))_Z / R^2
!            tempd(j,:) = &
!                 e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) & 
!               + e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) 
!         end do

!         ! [chi, f']_R*R 
!         temp79a = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
!         if(itor.eq.1) then 
!            temp79a = temp79a - ri_79* &        
!                 (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
!         end if 
!         ! [chi, f']_Z*R 
!         temp79b = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
!         ! ((chi, psi)/R^2)_R * R^2
!         temp79c = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &
!                 + f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ) 
!         if(itor.eq.1) then 
!            temp79c = temp79c - 2*ri_79*  & 
!                  (f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ)) 
!         end if 
!         ! ((chi, psi)/R^2)_Z * R^2
!         temp79d = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &
!                 + f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ) 

!         temp = intx3(tempa,temp79a,ri2_79) &
!              + intx3(tempb,temp79b,ri2_79) &
!              - intx3(tempc,temp79c,ri2_79) &
!              - intx3(tempd,temp79d,ri2_79) &
!              - intx4(e(:,:,OP_DZ),temp79a,ri2_79,g(:,OP_GS)) &
!              + intx4(e(:,:,OP_DR),temp79b,ri2_79,g(:,OP_GS)) 
!      end if

!   v1chipsif = temp
! end function v1chipsif

function v1chipsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chipsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd
  type(muarray) :: tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! ([nu, psi]*R^2)_R / R
           tempa = &
                mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &
              + mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR)
           if(itor.eq.1) then 
              tempa = tempa + &
                  (mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR))*ri_79 
           end if
           ! ([nu, psi]*R^2)_Z / R
           tempb = &
                mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &
              + mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR)
           ! (R^2*(nu, f'))_R / R^2
           tempc = &
                mu(h(:,OP_DZ),OP_DRZ) + mu(h(:,OP_DR),OP_DRR) & 
              + mu(h(:,OP_DRZ),OP_DZ) + mu(h(:,OP_DRR),OP_DR) 
           if(itor.eq.1) then 
              tempc = tempc + &
                  (mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR))*(2*ri_79) 
           end if
           ! (R^2*(nu, f'))_Z / R^2
           tempd = &
                mu(h(:,OP_DZ),OP_DZZ) + mu(h(:,OP_DR),OP_DRZ) & 
              + mu(h(:,OP_DZZ),OP_DZ) + mu(h(:,OP_DRZ),OP_DR) 

        ! [chi, f']_R*R 
        tempa2 = mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR) &           
               + mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then 
           tempa2 = tempa2 + &        
                (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79) 
        end if 
        ! [chi, f']_Z*R 
        tempb2 = mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ) &           
               + mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR)            
        ! ((chi, psi)/R^2)_R * R^2
        tempc2 = mu(g(:,OP_DR),OP_DRR) + mu(g(:,OP_DZ),OP_DRZ) &
               + mu(g(:,OP_DRR),OP_DR) + mu(g(:,OP_DRZ),OP_DZ) 
        if(itor.eq.1) then 
           tempc2 = tempc2 + & 
                 (mu(g(:,OP_DR),OP_DR) + mu(g(:,OP_DZ),OP_DZ))*(-2*ri_79) 
        end if 
        ! ((chi, psi)/R^2)_Z * R^2
        tempd2 = mu(g(:,OP_DR),OP_DRZ) + mu(g(:,OP_DZ),OP_DZZ) &
                + mu(g(:,OP_DRZ),OP_DR) + mu(g(:,OP_DZZ),OP_DZ) 

        temp = prod(tempa,tempa2*ri2_79) &
             + prod(tempb,tempb2*ri2_79) &
             + prod(tempc,tempc2*(-ri2_79)) &
             + prod(tempd,tempd2*(-ri2_79)) &
             + prod(mu(-ri2_79*g(:,OP_GS),OP_DZ),tempa2) &
             + prod(mu( ri2_79*g(:,OP_GS),OP_DR),tempb2) 
     end if

  v1chipsif = temp
end function v1chipsif

! V1chibf
! =====
! function v1chibf(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1chibf
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (nu, f')
!            tempa(j,:) = &
!               + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
!            ! [nu, f'']*R
!            tempb(j,:) = &
!               + e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
!         end do

!         ! [chi, F/R^4]'*R^5
!         temp79a = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP) 
!         if(itor.eq.1) then 
!            temp79a = temp79a - 4*ri_79* &        
!                 (f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) 
!         end if 
!         ! [chi, f']_R*R 
!         temp79b = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
!         if(itor.eq.1) then 
!            temp79b = temp79b - ri_79* &        
!                 (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
!         end if 
!         ! [chi, f']_Z*R 
!         temp79c = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
!         ! (chi, F/R^2)*R^2 + F*chi_GS
!         temp79d = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
!                 + f(:,OP_GS)*g(:,OP_1) 
!         if(itor.eq.1) then 
!            temp79d = temp79d - 2*ri_79*f(:,OP_DR)*g(:,OP_1) 
!         end if 

!         temp = - intx3(tempa,temp79a,ri3_79) &
!                - intx3(tempb,temp79d,ri3_79) &
!                + intx4(e(:,:,OP_DZ),temp79c,ri3_79,g(:,OP_DP)) &
!                + intx4(e(:,:,OP_DR),temp79b,ri3_79,g(:,OP_DP)) &
!                + intx4(e(:,:,OP_DZP),temp79c,ri3_79,g(:,OP_1)) &
!                + intx4(e(:,:,OP_DRP),temp79b,ri3_79,g(:,OP_1)) 
!      end if

!   v1chibf = temp
! end function v1chibf

function v1chibf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chibf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb
  type(muarray) :: tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f')
           tempa = &
              mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR) 
           ! [nu, f'']*R
           tempb = &
              mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR) 

        ! [chi, F/R^4]'*R^5
        tempa2 = mu(g(:,OP_DR),OP_DZP) + mu(-g(:,OP_DZ),OP_DRP) &           
               + mu(g(:,OP_DRP),OP_DZ) + mu(-g(:,OP_DZP),OP_DR) 
        if(itor.eq.1) then 
           tempa2 = tempa2 + &        
                (mu(g(:,OP_1),OP_DZP) + mu(g(:,OP_DP),OP_DZ))*(-4*ri_79) 
        end if 
        ! [chi, f']_R*R 
        tempb2 = mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR) &           
               + mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then 
           tempb2 = tempb2 + &        
                (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79) 
        end if 
        ! [chi, f']_Z*R 
        tempc2 = mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ) &           
               + mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR)            
        ! (chi, F/R^2)*R^2 + F*chi_GS
        tempd2 = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR) &
               + mu(g(:,OP_1),OP_GS) 
        if(itor.eq.1) then 
           tempd2 = tempd2 + mu(-2*ri_79*g(:,OP_1),OP_DR) 
        end if 

        temp = prod(tempa,tempa2*(-ri3_79)) &
             + prod(tempb,tempd2*(-ri3_79)) &
             + prod(mu(ri3_79*g(:,OP_DP),OP_DZ),tempc2) &
             + prod(mu(ri3_79*g(:,OP_DP),OP_DR),tempb2) &
             + prod(mu(ri3_79*g(:,OP_1),OP_DZP),tempc2) &
             + prod(mu(ri3_79*g(:,OP_1),OP_DRP),tempb2) 
     end if

  v1chibf = temp
end function v1chibf

! V1chiff
! =====
! function v1chiff(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v1chiff
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j

!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (R^2*(nu, f'))_R / R^2
!            tempa(j,:) = &
!                 e(j,:,OP_DRZ)*h(:,OP_DZ) + e(j,:,OP_DRR)*h(:,OP_DR) & 
!               + e(j,:,OP_DZ)*h(:,OP_DRZ) + e(j,:,OP_DR)*h(:,OP_DRR) 
!            if(itor.eq.1) then 
!               tempa(j,:) = tempa(j,:) + 2*ri_79* &
!                   (e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)) 
!            end if
!            ! (R^2*(nu, f'))_Z / R^2
!            tempb(j,:) = &
!                 e(j,:,OP_DZZ)*h(:,OP_DZ) + e(j,:,OP_DRZ)*h(:,OP_DR) & 
!               + e(j,:,OP_DZ)*h(:,OP_DZZ) + e(j,:,OP_DR)*h(:,OP_DRZ) 
!         end do

!         ! [chi, f']_R*R
!         temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
!         if(itor.eq.1) then 
!            temp79a = temp79a - ri_79* &        
!                 (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
!         end if 
!         ! [chi, f']_Z*R
!         temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            

!         temp = - intx3(tempa,temp79a,ri_79) &
!                - intx3(tempb,temp79b,ri_79) 
!      end if

!   v1chiff = temp
! end function v1chiff

function v1chiff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chiff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb
  type(muarray) :: tempa2, tempb2
  integer :: j

     if(surface_int) then
        temp%len = 0
     else
           ! (R^2*(nu, f'))_R / R^2
           tempa = &
                mu(h(:,OP_DZ),OP_DRZ) + mu(h(:,OP_DR),OP_DRR) & 
              + mu(h(:,OP_DRZ),OP_DZ) + mu(h(:,OP_DRR),OP_DR) 
           if(itor.eq.1) then 
              tempa = tempa + &
                  (mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR))*(2*ri_79) 
           end if
           ! (R^2*(nu, f'))_Z / R^2
           tempb = &
                mu(h(:,OP_DZ),OP_DZZ) + mu(h(:,OP_DR),OP_DRZ) & 
              + mu(h(:,OP_DZZ),OP_DZ) + mu(h(:,OP_DRZ),OP_DR) 

        ! [chi, f']_R*R
        tempa2 = mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &           
               + mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then 
           tempa2 = tempa2 + &        
                (mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR))*(-ri_79) 
        end if 
        ! [chi, f']_Z*R
        tempb2 = mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &           
               + mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR)            

        temp = prod(tempa,tempa2*(-ri_79)) &
             + prod(tempb,tempb2*(-ri_79)) 
     end if

  v1chiff = temp
end function v1chiff
#endif

! V1chip
! ======
!function v1chip(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v1chip
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

     !if(itor.eq.0) then
        !temp = 0.
     !else
        !temp79a = gam*f(:,OP_GS)*g(:,OP_1) + &
             !f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

        !if(surface_int) then
           !temp = -2.*intx4(e(:,:,OP_1),ri2_79,temp79a,norm79(:,2))
        !else
           !temp = 2.*intx3(e(:,:,OP_DZ),ri2_79,temp79a)
        !end if
     !end if

  !v1chip = temp
!end function v1chip

function v1chip(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chip
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp
  type(muarray) :: tempa

     if(itor.eq.0) then
        temp%len = 0
     else
        tempa = mu(gam*g(:,OP_1),OP_GS) + &
             mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)

        if(surface_int) then
           temp = prod(mu(-2.*ri2_79*norm79(:,2),OP_1),tempa)
        else
           temp = prod(mu(2.*ri2_79,OP_DZ),tempa)
        end if
     end if

  v1chip = temp
end function v1chip


! V1ngrav
! =======
function v1ngrav(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ngrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ngrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = gravz*intx3(e(:,:,OP_1), r_79,f(:,OP_DR)) &
          - gravr*intx3(e(:,:,OP_1),ri_79,f(:,OP_DZ))
  end if

  v1ngrav = temp
end function v1ngrav


! V1ungrav
! ========
function v1ungrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ungrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ungrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)
  
     temp = gravz*intx2(e(:,:,OP_DR),temp79a) &
          - gravr*intx3(e(:,:,OP_DZ),ri2_79,temp79a)
     
     if(itor.eq.1) &
          temp = temp + 2.*gravz*intx3(e(:,:,OP_1),ri_79,temp79a)
  end if

  v1ungrav = temp
end function v1ungrav


! V1chingrav
! ==========
function v1chingrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chingrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1chingrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = r_79*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
          + g(:,OP_1)*f(:,OP_LP))

     temp = gravz*intx2(e(:,:,OP_DR),temp79a) &
          - gravr*intx3(e(:,:,OP_DZ),ri2_79,temp79a)

     if(itor.eq.1) &
          temp = temp + 2.*gravz*intx3(e(:,:,OP_1),ri_79,temp79a)
  endif

  v1chingrav = temp
end function v1chingrav


! V1ndenmgrav
! ===========
function v1ndenmgrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1ndenmgrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  real, intent(in) :: g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ndenmgrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = -g*r_79*f(:,OP_LP)

     temp = gravz*intx2(e(:,:,OP_DR),temp79a) &
          - gravr*intx3(e(:,:,OP_DZ),ri2_79,temp79a)

     if(itor.eq.1) &
          temp = temp + 2.*gravz*intx3(e(:,:,OP_1),ri_79,temp79a)
  end if

  v1ndenmgrav = temp
end function v1ndenmgrav


! V1us
! ====
function v1us(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1us
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp =  &
                + intx5(e(:,:,OP_1),r2_79,norm79(:,1),f(:,OP_DR),temp79a) &
                + intx5(e(:,:,OP_1),r2_79,norm79(:,2),f(:,OP_DZ),temp79a)
        endif
     else
        temp = -intx4(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),temp79a) &
             -  intx4(e(:,:,OP_DR),r2_79,f(:,OP_DR),temp79a)
     end if

  v1us = temp
end function v1us

function v1us1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1us1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1us1%len = 0
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len = 0
        else
           temp =  &
                  prod(r2_79*norm79(:,1)*temp79a,OP_1,OP_DR) &
                + prod(r2_79*norm79(:,2)*temp79a,OP_1,OP_DZ)
        endif
     else
        temp = prod(-r2_79*temp79a,OP_DZ,OP_DZ) &
             + prod(-r2_79*temp79a,OP_DR,OP_DR)
     end if

  v1us1 = temp
end function v1us1

! V1chis
! ======
function v1chis(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1chis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)


     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
                + intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DZ),temp79a) &
                - intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DR),temp79a)
        endif
     else
        temp = intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DR),temp79a) &
             - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZ),temp79a)
     endif

  v1chis = temp
end function v1chis

function v1chis1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1chis1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1chis1%len = 0
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)


     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len = 0
        else
           temp = &
                  prod( ri_79*norm79(:,1)*temp79a,OP_1,OP_DZ) &
                + prod(-ri_79*norm79(:,2)*temp79a,OP_1,OP_DR)
        endif
     else
        temp = prod( ri3_79*temp79a,OP_DZ,OP_DR) &
             + prod(-ri3_79*temp79a,OP_DR,OP_DZ)
     endif

  v1chis1 = temp
end function v1chis1

! V1psif
! ======
function v1psif(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1psif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp
  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp =  &
                + intx4(e(:,:,OP_1),f(:,OP_GS),norm79(:,1),g(:,OP_DR)) &
                + intx4(e(:,:,OP_1),f(:,OP_GS),norm79(:,2),g(:,OP_DZ))
        endif
     else
        temp = &
             - intx3(e(:,:,OP_DZ),f(:,OP_GS),g(:,OP_DZ)) &
             - intx3(e(:,:,OP_DR),f(:,OP_GS),g(:,OP_DR))
     end if
  
#else
  temp = 0.
#endif
  v1psif = temp
end function v1psif

function v1psif1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1psif1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp
  temp%len =0
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len =0
        else
           temp =  &
                  prod(norm79(:,1)*g(:,OP_DR),OP_1,OP_GS) &
                + prod(norm79(:,2)*g(:,OP_DZ),OP_1,OP_GS)
        endif
     else
        temp = &
               prod(-g(:,OP_DZ),OP_DZ,OP_GS) &
             + prod(-g(:,OP_DR),OP_DR,OP_GS)
     end if
  
#else
  temp%len =0
#endif
  v1psif1 = temp
end function v1psif1

function v1psif2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1psif2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp
  temp%len =0
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len =0
        else
           temp =  &
                  prod(f(:,OP_GS)*norm79(:,1),OP_1,OP_DR) &
                + prod(f(:,OP_GS)*norm79(:,2),OP_1,OP_DZ)
        endif
     else
        temp = &
               prod(-f(:,OP_GS),OP_DZ,OP_DZ) &
             + prod(-f(:,OP_GS),OP_DR,OP_DR)
     end if
  
#else
  temp%len =0
#endif
  v1psif2 = temp
end function v1psif2


! V1bf
! ====
function v1bf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1bf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
                + intx5(e(:,:,OP_1),ri_79,f(:,OP_1),norm79(:,1),g(:,OP_DZP)) &
                - intx5(e(:,:,OP_1),ri_79,f(:,OP_1),norm79(:,2),g(:,OP_DRP))
        end if
     else
        temp = &
             + intx4(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DRP)) &
             - intx4(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZP))
     end if
#else
  temp = 0.
#endif
  v1bf = temp
end function v1bf

function v1bf1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1bf1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len =0
        else
           temp = &
                  prod( ri_79*norm79(:,1)*g(:,OP_DZP),OP_1,OP_1) &
                + prod(-ri_79*norm79(:,2)*g(:,OP_DRP),OP_1,OP_1)
        end if
     else
        temp = &
               prod( ri_79*g(:,OP_DRP),OP_DZ,OP_1) &
             + prod(-ri_79*g(:,OP_DZP),OP_DR,OP_1)
     end if
#else
  temp%len =0
#endif
  v1bf1 = temp
end function v1bf1

function v1bf2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1bf2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp%len =0
        else
           temp = &
                  prod( ri_79*f(:,OP_1)*norm79(:,1),OP_1,OP_DZP) &
                + prod(-ri_79*f(:,OP_1)*norm79(:,2),OP_1,OP_DRP)
        end if
     else
        temp = &
               prod( ri_79*f(:,OP_1),OP_DZ,OP_DRP) &
             + prod(-ri_79*f(:,OP_1),OP_DR,OP_DZP)
     end if
#else
  temp%len =0
#endif
  v1bf2 = temp
end function v1bf2


! V1p
! ===
function v1p(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1p
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  vectype, dimension(dofs_per_element) :: temp

  if(itor.eq.0) then
     v1p = 0.
     return
  end if

  temp = 0.
     if(surface_int) then
        if(inoslip_pol.eq.1 .or. iconst_p.ge.1) then
           temp = 0.
        else
           temp = &
                + intx4(e(:,:,OP_1),r_79,norm79(:,1),f(:,OP_DZ)) &
                - intx4(e(:,:,OP_1),r_79,norm79(:,2),f(:,OP_DR))
        endif
     else
        temp = intx3(e(:,:,OP_DZ),r_79,f(:,OP_DR)) &
             - intx3(e(:,:,OP_DR),r_79,f(:,OP_DZ))
     end if

  v1p = temp
end function v1p

function v1p1
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v1p1

  type(prodarray) :: temp

  if(itor.eq.0) then
     v1p1%len = 0
     return
  end if

     if(surface_int) then
        if(inoslip_pol.eq.1 .or. iconst_p.ge.1) then
           temp%len = 0
        else
           temp = prod( r_79*norm79(:,1),OP_1,OP_DZ) &
                + prod(-r_79*norm79(:,2),OP_1,OP_DR)
         endif
     else
        temp = prod( r_79,OP_DZ,OP_DR) &
             + prod(-r_79,OP_DR,OP_DZ)
      end if

  v1p1 = temp
end function v1p1

! V1psiforce
! ===
vectype function v1psiforce(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp = int3(g(:,OP_1),e(:,OP_DR),f(:,OP_DR)) &
             + int3(g(:,OP_1),e(:,OP_DZ),f(:,OP_DZ))
     end if


  v1psiforce = temp
  return
end function v1psiforce

! V1be
! ===

function v1par(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1par
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp = - intx2(e(:,:,OP_DZ),f(:,OP_1))
     end if

  v1par = temp
end function v1par

function v1parb2ipsipsi(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1parb2ipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79a = f(:,OP_1)*g(:,OP_1)*ri_79
        temp =  intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DR),i(:,OP_DRR))   &
             +  intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DZ),i(:,OP_DRZ))   &
             -  intx4(e(:,:,OP_DR),temp79a,h(:,OP_DR),i(:,OP_DRZ))   &
             -  intx4(e(:,:,OP_DR),temp79a,h(:,OP_DZ),i(:,OP_DZZ))   

        temp79b = -f(:,OP_1)*ri_79
        temp = temp                                                         &
             +  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DR))  &
             +  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ))  &
             -  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DR))  &
             -  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

        temp79c = -g(:,OP_1)*ri_79
        temp = temp                                                         &
             +  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DR)) &
             +  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ)) &
             -  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DR)) &
             -  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

        temp79d = -f(:,OP_1)*g(:,OP_1)*h(:,OP_GS)*ri_79
        temp = temp                                   &
             +  intx3(e(:,:,OP_DZ),temp79d,i(:,OP_DR))   &
             -  intx3(e(:,:,OP_DR),temp79d,i(:,OP_DZ))   
     end if

  v1parb2ipsipsi = temp
end function v1parb2ipsipsi

function v1parb2ipsib(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1parb2ipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  if(surface_int) then
     temp = 0.
  else
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = - f(:,OP_DP)*g(:,OP_1)*i(:,OP_1)*ri2_79
     temp = +intx3(e(:,:,OP_DR),temp79a,h(:,OP_DR)) &
          +  intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DZ))
#endif
  end if

  v1parb2ipsib = temp
end function v1parb2ipsib

#ifdef USEPARTICLES
function v1p_2(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1p_2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        if(inoslip_pol.eq.1 .or. iconst_p.eq.1) then
           temp = 0.
        else
           temp = &
                + intx5(e(:,:,OP_1),r_79,norm79(:,1),f(:,OP_DZ),g) &
                - intx5(e(:,:,OP_1),r_79,norm79(:,2),f(:,OP_DR),g)
        endif
     else
        temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_DR),g) &
             - intx4(e(:,:,OP_DR),r_79,f(:,OP_DZ),g)
     end if

  v1p_2 = temp
end function v1p_2

function v1pbb(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1pbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79b = -g*ri_79
        temp =  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DR)) &
              + intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DZ)) &
              - intx5(e(:,:,OP_DR),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DR)) &
              - intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DZ))
  #if defined(USE3D) || defined(USECOMPLEX)
        temp79a = - f(:,OP_DP)*bzt79(:,OP_1)*g*ri2_79
        temp = temp +intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DR)) &
              +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ))
        temp79a = g
        temp = temp +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
        temp79a = -g
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DZ)) 
        temp79a = -f(:,OP_DP)*g*bzt79(:,OP_1)*ri_79
        temp = temp +intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR)) &
              -intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ)) 
        temp79a = g*r_79
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZP),f(:,OP_DR),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZP),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
        #endif
    end if

  v1pbb = temp
end function v1pbb

function v1jxb(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v1jxb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79a = -f*pst79(:,OP_GS)*ri_79
        temp = temp                                   &
             +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR))   &
             -  intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ))   
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = - f*bzt79(:,OP_1)*ri_79
     temp = temp + intx3(e(:,:,OP_DZ),temp79a,(bzt79(:,OP_DR)+bfpt79(:,OP_DRP))) &
          - intx3(e(:,:,OP_DR),temp79a,bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))
     temp79a = -f*bzt79(:,OP_1)*ri2_79
     temp = temp+intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DRP)) &
          +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZP))
     temp79a = f*pst79(:,OP_GS)
     temp = temp+intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR)) &
          +  intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ))
     #endif
    end if

  v1jxb = temp
end function v1jxb

#endif
!============================================================================
! V2 TERMS
!============================================================================


! V2vn
! ====
!function v2vn(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none
  
  !vectype, dimension(dofs_per_element) :: v2vn
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  !if(surface_int) then
     !v2vn = 0.
     !return
  !endif

  !v2vn = intx4(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_1))
!end function v2vn

function v2vn(g)

  use basic
  use m3dc1_nint

  implicit none
  
  type(prodarray) :: v2vn
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  if(surface_int) then
     v2vn%len = 0
     return
  endif

     v2vn = prod(r2_79*g(:,OP_1),OP_1,OP_1)
end function v2vn


! V2umu
! =====
!function v2umu(e,f,g,h)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v2umu
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  !vectype, dimension(dofs_per_element) :: temp

  !temp = 0.
!#if defined(USE3D) || defined(USECOMPLEX)
     !if(surface_int) then
        !temp = intx5(e(:,:,OP_1),r_79,norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
             !- intx5(e(:,:,OP_1),r_79,norm79(:,1),f(:,OP_DZP),g(:,OP_1))
     !else
        !temp = intx4(e(:,:,OP_DR),r_79,f(:,OP_DZP),g(:,OP_1)) &
             !- intx4(e(:,:,OP_DZ),r_79,f(:,OP_DRP),g(:,OP_1))
        !if(itor.eq.1) then
           !temp = temp &
                !+ 2.*intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_1)) &
                !- 4.*intx3(e(:,:,OP_1),f(:,OP_DZP),h(:,OP_1))
        !endif
     !end if

!#endif
  !v2umu = temp
!end function v2umu

function v2umu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2umu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp = prod( r_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DRP) &
             + prod(-r_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DZP)
      else
        temp = prod( r_79*g(:,OP_1),OP_DR,OP_DZP) &
             + prod(-r_79*g(:,OP_1),OP_DZ,OP_DRP)
        if(itor.eq.1) then
           temp = temp &
                + prod( 2.*g(:,OP_1),OP_1,OP_DZP) &
                + prod(-4.*h(:,OP_1),OP_1,OP_DZP)
         endif
     end if
#endif
  v2umu = temp
end function v2umu


! V2vmu
! =====
!function v2vmu(e,f,g,h)

!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v2vmu
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
!  vectype, dimension(dofs_per_element) :: temp

!  if(surface_int) then
!     temp = intx5(e(:,:,OP_1),r2_79,norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
!          + intx5(e(:,:,OP_1),r2_79,norm79(:,2),f(:,OP_DZ),g(:,OP_1))
!  else
!     temp = -intx4(e(:,:,OP_DZ),r2_79,f(:,OP_DZ),g(:,OP_1)) &
!          -  intx4(e(:,:,OP_DR),r2_79,f(:,OP_DR),g(:,OP_1))
     
!#if defined(USE3D) || defined(USECOMPLEX)
!     temp = temp + 2.*intx3(e(:,:,OP_1),f(:,OP_DPP),h(:,OP_1))
!#endif
!  end if

!  v2vmu = temp
!end function v2vmu

function v2vmu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vmu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

     if(surface_int) then
        temp = prod(r2_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DR) &
             + prod(r2_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DZ)
      else
        temp = prod(-r2_79*g(:,OP_1),OP_DZ,OP_DZ) &
             + prod(-r2_79*g(:,OP_1),OP_DR,OP_DR)
 
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp + prod(2.*h(:,OP_1),OP_1,OP_DPP)
#endif
    end if

  v2vmu = temp
end function v2vmu


! V2chimu
! =======
! function v2chimu(e,f,g,h)

!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2chimu
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
!   vectype, dimension(dofs_per_element) :: temp

!   temp = 0.
! #if defined(USE3D) || defined(USECOMPLEX)
!      if(surface_int) then
!         temp = intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
!              + intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZP),g(:,OP_1))
!      else
!         temp79a = h(:,OP_1) - g(:,OP_1)
!         temp = &
!              - intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
!              - intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1)) &
! #ifdef USEST
!              - 2.*intx4(e(:,:,OP_DP),ri2_79,f(:,OP_GS),temp79a) &
!              - 2.*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GS),h(:,OP_DP)) &
!              + 2.*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GS),g(:,OP_DP))
! #else
!              + 2.*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GSP),temp79a)
! #endif
!         if(itor.eq.1) then
!            temp = temp &
!                 +2.*intx4(e(:,:,OP_1),ri3_79,f(:,OP_DRP),g(:,OP_1))
!         endif
!      end if

! #endif
!   v2chimu = temp
! end function v2chimu

function v2chimu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chimu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp = prod(ri2_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DRP) &
             + prod(ri2_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DZP)
      else
        temp79a = h(:,OP_1) - g(:,OP_1)
        temp = &
               prod(-ri2_79*g(:,OP_1),OP_DZ,OP_DZP) &
             + prod(-ri2_79*g(:,OP_1),OP_DR,OP_DRP) &
#ifdef USEST
             + prod(-2.*ri2_79*temp79a,OP_DP,OP_GS) &
             + prod(-2.*ri2_79*h(:,OP_DP),OP_1,OP_GS) &
             + prod( 2.*ri2_79*g(:,OP_DP),OP_1,OP_GS)
#else
              + prod(2.*ri2_79*temp79a,OP_1,OP_GSP)
#endif
         if(itor.eq.1) then
           temp = temp &
                +prod(2.*ri3_79*g(:,OP_1),OP_1,OP_DRP)
        endif
     end if
#endif
  v2chimu = temp
end function v2chimu



! V2vun
! =====
function v2vun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vun
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vun = 0.
     return
  end if


     temp = intx5(e(:,:,OP_1),r3_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
          - intx5(e(:,:,OP_1),r3_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1))

     if(itor.eq.1) then
        temp = temp + &
             2.*intx5(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_1))
     end if

  v2vun = temp
end function v2vun

function v2vun1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vun1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vun1%len = 0
     return
  end if

     temp = prod( r3_79*g(:,OP_DZ)*h(:,OP_1),OP_1,OP_DR) &
          + prod(-r3_79*g(:,OP_DR)*h(:,OP_1),OP_1,OP_DZ)

     if(itor.eq.1) then
        temp = temp + &
             prod(2.*r2_79*g(:,OP_DZ)*h(:,OP_1),OP_1,OP_1)
     end if

  v2vun1 = temp
end function v2vun1

function v2vun2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vun2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vun2%len = 0
     return
  end if

     temp = prod( r3_79*f(:,OP_DR)*h(:,OP_1),OP_1,OP_DZ) &
          + prod(-r3_79*f(:,OP_DZ)*h(:,OP_1),OP_1,OP_DR)

     if(itor.eq.1) then
        temp = temp + &
             prod(2.*r2_79*f(:,OP_1)*h(:,OP_1),OP_1,OP_DZ)
     end if

  v2vun2 = temp
end function v2vun2

! V2vvn
! =====
function v2vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vvn = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp = -intx5(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_DP),h(:,OP_1))
#else
  temp = 0.
#endif
  v2vvn = temp
end function v2vvn

function v2vvn1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vvn1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vvn1%len = 0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
     temp = prod(-r2_79*g(:,OP_DP)*h(:,OP_1),OP_1,OP_1)
#else
  temp%len = 0
#endif
  v2vvn1 = temp
end function v2vvn1

function v2vvn2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vvn2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vvn2%len = 0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
     temp = prod(-r2_79*f(:,OP_1)*h(:,OP_1),OP_1,OP_DP)
#else
  temp%len = 0
#endif
  v2vvn2 = temp
end function v2vvn2


! V2up
! ====
!function v2up(e,f,g)
  !use basic
  !use arrays
  !use m3dc1_nint
  
  !implicit none

  !vectype, dimension(dofs_per_element) :: v2up
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g
  !vectype, dimension(dofs_per_element) :: temp

  !temp = 0.

  !if(surface_int) then
     !v2up = 0.
     !return
  !end if

!#if defined(USE3D) || defined(USECOMPLEX)
     !temp = intx4(e(:,:,OP_1),r_79,f(:,OP_DRP),g(:,OP_DZ)) &
          !- intx4(e(:,:,OP_1),r_79,f(:,OP_DZP),g(:,OP_DR)) &
          !+ intx4(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_DZP)) &
          !- intx4(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_DRP))
     !if(itor.eq.1) then
        !temp = temp - 2.*gam* &
             !(intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_1)) &
             !+intx3(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DP)))
     !endif
!#endif

  !v2up = temp
!end function v2up

function v2up(g)
  use basic
  use arrays
  use m3dc1_nint
  
  implicit none

  type(prodarray) :: v2up
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  temp%len = 0

  if(surface_int) then
     v2up%len = 0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
     temp = prod( r_79*g(:,OP_DZ),OP_1,OP_DRP) &
          + prod(-r_79*g(:,OP_DR),OP_1,OP_DZP) &
          + prod( r_79*g(:,OP_DZP),OP_1,OP_DR) &
          + prod(-r_79*g(:,OP_DRP),OP_1,OP_DZ)
     if(itor.eq.1) then
        temp = temp + &
            (prod(-2.*gam*g(:,OP_1),OP_1,OP_DZP) &
            +prod(-2.*gam*g(:,OP_DP),OP_1,OP_DZ))
     endif
#endif

  v2up = temp
end function v2up

! V2vp
! ====
!function v2vp(e,f,g)
  !use basic
  !use arrays
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v2vp
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g
  !vectype, dimension(dofs_per_element) :: temp

  !if(surface_int) then
     !v2vp = 0.
     !return
  !end if

  !temp = 0.

!#if defined(USE3D) || defined(USECOMPLEX)
     !temp =          intx3(e(:,:,OP_1),g(:,OP_DPP),f(:,OP_1)) &
          !+ (1.+gam)*intx3(e(:,:,OP_1),g(:,OP_DP),f(:,OP_DP)) &
          !+ gam*     intx3(e(:,:,OP_1),g(:,OP_1),f(:,OP_DPP))
!#endif

  !v2vp = temp
!end function v2vp

function v2vp(g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vp
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(surface_int) then
     v2vp%len = 0
     return
  end if

  temp%len = 0

#if defined(USE3D) || defined(USECOMPLEX)
     temp =          prod(g(:,OP_DPP),OP_1,OP_1) &
          +          prod((1.+gam)*g(:,OP_DP),OP_1,OP_DP) &
          +          prod(gam*g(:,OP_1),OP_1,OP_DPP)
#endif

  v2vp = temp
end function v2vp


! V2chip
! ======
! function v2chip(e,f,g)

!   use basic
!   use arrays
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2chip
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g
!   vectype, dimension(dofs_per_element) :: temp

!   if(surface_int) then
!      v2chip = 0.
!      return
!   end if

!   temp = 0.

! #if defined(USE3D) || defined(USECOMPLEX)
!      temp =     intx4(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR))    &
!               + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ))    &
!               + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP))    &
!               + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP))    &
! #ifdef USEST
!           - gam*intx4(e(:,:,OP_DP),ri2_79,f(:,OP_GS),g(:,OP_1 ))    
! #else
!           + gam*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GSP),g(:,OP_1 ))    &
!           + gam*intx4(e(:,:,OP_1),ri2_79,f(:,OP_GS),g(:,OP_DP))
! #endif
! #endif

!   v2chip = temp
! end function v2chip

function v2chip(g)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chip
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(surface_int) then
     v2chip%len = 0
     return
  end if

  temp%len = 0

#if defined(USE3D) || defined(USECOMPLEX)
     temp =     prod(ri2_79*g(:,OP_DR),OP_1,OP_DRP)    &
              + prod(ri2_79*g(:,OP_DZ),OP_1,OP_DZP) &
              + prod(ri2_79*g(:,OP_DRP),OP_1,OP_DR) &
              + prod(ri2_79*g(:,OP_DZP),OP_1,OP_DZ) &
#ifdef USEST
          + prod(-gam*ri2_79*g(:,OP_1 ),OP_DP,OP_GS)    
#else
          + prod(gam*ri2_79*g(:,OP_1),OP_1,OP_GSP) &
          + prod(gam*ri2_79*g(:,OP_DP),OP_1,OP_GS)
#endif
#endif

  v2chip = temp
end function v2chip


! V2p
! ===
function v2p(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2p
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  if(surface_int) then
     v2p = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2p = -intx2(e(:,:,OP_1),f(:,OP_DP))
#else
  v2p = 0.
#endif

end function v2p

function v2p1
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2p1

  if(surface_int) then
     v2p1%len = 0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = 1.
  v2p1 = prod(-temp79a,OP_1,OP_DP)
#else
  v2p1%len = 0
#endif

end function v2p1


! V2psipsi
! ========
function v2psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(surface_int) then
     v2psipsi = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psipsi = - &
       (intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ)) &
       +intx4(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR)))
#else
  v2psipsi = 0.
#endif
end function v2psipsi

function v2psipsi1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psipsi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  if(surface_int) then
     v2psipsi1%len = 0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psipsi1 = &
       (prod(-ri2_79*g(:,OP_DZ),OP_1,OP_DZP) &
       +prod(-ri2_79*g(:,OP_DR),OP_1,OP_DRP))
#else
  v2psipsi1%len = 0
#endif
end function v2psipsi1

function v2psipsi2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psipsi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  if(surface_int) then
     v2psipsi2%len = 0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psipsi2 = &
       (prod(-ri2_79*f(:,OP_DZP),OP_1,OP_DZ) &
       +prod(-ri2_79*f(:,OP_DRP),OP_1,OP_DR))
#else
  v2psipsi2%len = 0
#endif
end function v2psipsi2


! V2psib
! ======
function v2psib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(surface_int) then
     v2psib = 0.
     return
  end if

  v2psib = intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
       -   intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR))
end function v2psib

function v2psib1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psib1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  if(surface_int) then
     v2psib1%len = 0
     return
  end if

  v2psib1 = prod( ri_79*g(:,OP_DZ),OP_1,OP_DR) &
          + prod(-ri_79*g(:,OP_DR),OP_1,OP_DZ)
end function v2psib1

function v2psib2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psib2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  if(surface_int) then
     v2psib2%len = 0
     return
  end if

  v2psib2 = prod( ri_79*f(:,OP_DR),OP_1,OP_DZ) &
         + prod(-ri_79*f(:,OP_DZ),OP_1,OP_DR)
end function v2psib2

! V2vpsipsi
! =========
! function v2vpsipsi(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2vpsipsi
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  
!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
!   integer :: j

!   do j=1, dofs_per_element
!      ! [nu,psi(2)]
!      tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
!   end do

!      if(surface_int) then
!         temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
!         temp = intx4(e(:,:,OP_1),temp79a,norm79(:,2),h(:,OP_DR)) &
!              - intx4(e(:,:,OP_1),temp79a,norm79(:,1),h(:,OP_DZ))
!      else
!         temp = intx3(tempa,f(:,OP_DR),g(:,OP_DZ)) &
!              - intx3(tempa,f(:,OP_DZ),g(:,OP_DR))

! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         temp79b = &
!               f(:,OP_DP)*(g(:,OP_DZ )*h(:,OP_DZ ) + g(:,OP_DR )*h(:,OP_DR )) &
!           + f(:,OP_1 )*(g(:,OP_DZP)*h(:,OP_DZ ) + g(:,OP_DRP)*h(:,OP_DR )) 
!         temp = temp - intx3(e(:,:,OP_DP),ri2_79,temp79b)
! #else
!         temp79b = &
!               f(:,OP_DPP)*(g(:,OP_DZ )*h(:,OP_DZ ) + g(:,OP_DR )*h(:,OP_DR )) &
!           + 2.*f(:,OP_DP)*(g(:,OP_DZP)*h(:,OP_DZ ) + g(:,OP_DRP)*h(:,OP_DR )) &
!           +   f(:,OP_DP )*(g(:,OP_DZ )*h(:,OP_DZP) + g(:,OP_DR )*h(:,OP_DRP)) &
!           + f(:,OP_1 )*(g(:,OP_DZPP)*h(:,OP_DZ ) + g(:,OP_DRPP)*h(:,OP_DR )) &
!           + f(:,OP_1 )*(g(:,OP_DZP )*h(:,OP_DZP) + g(:,OP_DRP )*h(:,OP_DRP))
!         temp = temp + intx3(e(:,:,OP_1),ri2_79,temp79b)
! #endif
! #endif
!      end if

!   v2vpsipsi = temp
! end function v2vpsipsi

function v2vpsipsi(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vpsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  tempa = mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR)

     if(surface_int) then
        tempb = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
        temp = prod(mu( norm79(:,2)*h(:,OP_DR),OP_1),tempb) &
             + prod(mu(-norm79(:,1)*h(:,OP_DZ),OP_1),tempb)
      else
        temp = prod(tempa,mu( g(:,OP_DZ),OP_DR)) &
             + prod(tempa,mu(-g(:,OP_DR),OP_DZ))

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
        tempb = &
              mu(g(:,OP_DZ )*h(:,OP_DZ ) + g(:,OP_DR )*h(:,OP_DR ),OP_DP) &
          +   mu(g(:,OP_DZP)*h(:,OP_DZ ) + g(:,OP_DRP)*h(:,OP_DR ),OP_1) 
        temp = temp + prod(mu(-ri2_79,OP_DP),tempb)
#else
         tempb = &
              mu((g(:,OP_DZ )*h(:,OP_DZ ) + g(:,OP_DR )*h(:,OP_DR )),OP_DPP) &
          +   mu(2.*(g(:,OP_DZP)*h(:,OP_DZ ) + g(:,OP_DRP)*h(:,OP_DR )),OP_DP) &
          +   mu((g(:,OP_DZ )*h(:,OP_DZP) + g(:,OP_DR )*h(:,OP_DRP)),OP_DP) &
          +   mu((g(:,OP_DZPP)*h(:,OP_DZ ) + g(:,OP_DRPP)*h(:,OP_DR )),OP_1) &
          +   mu((g(:,OP_DZP )*h(:,OP_DZP) + g(:,OP_DRP )*h(:,OP_DRP)),OP_1)
        temp = temp + prod(mu(ri2_79,OP_1),tempb)
#endif
#endif
     end if

  v2vpsipsi = temp
end function v2vpsipsi


! V2vpsib
! =======
!function v2vpsib(e,f,g,h)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v2vpsib
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!#if defined(USE3D) || defined(USECOMPLEX)
  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !temp = 0.
     !else
        !temp79a = f(:,OP_DP)*(g(:,OP_DZ )*h(:,OP_DR)-g(:,OP_DR )*h(:,OP_DZ)) &
             !+    f(:,OP_1 )*(g(:,OP_DZP)*h(:,OP_DR)-g(:,OP_DRP)*h(:,OP_DZ))
        !temp = intx3(e(:,:,OP_1),ri_79,temp79a)
     !end if

  !v2vpsib = temp
!#else
  !v2vpsib = 0.
!#endif
!end function v2vpsib

function v2vpsib(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vpsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

     if(surface_int) then
        temp%len = 0
     else
        tempa = mu((g(:,OP_DZ )*h(:,OP_DR)-g(:,OP_DR )*h(:,OP_DZ)),OP_DP) &
             +  mu((g(:,OP_DZP)*h(:,OP_DR)-g(:,OP_DRP)*h(:,OP_DZ)),OP_1)
        temp = prod(mu(ri_79,OP_1),tempa)
     end if
  v2vpsib = temp
#else
  v2vpsib%len = 0
#endif
end function v2vpsib

#ifdef USE3D
! V2vpsif
! =====
! function v2vpsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2vpsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (nu, f') + nu*f'_LP
!            tempa(j,:) = e(j,:,OP_1)*h(:,OP_LP) &
!               + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
!            ! [nu, psi]*R
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
!         end do

!         ! [v, psi]*R
!         temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)            
!         ! (v, f') + v*f'_LP
!         temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
!                 + f(:,OP_1)*h(:,OP_LP)
!         ! [psi,f']*R
!         temp79c = g(:,OP_DZ)*h(:,OP_DR) - g(:,OP_DR)*h(:,OP_DZ)            
!         ! [psi,f']'*R
!         temp79d = g(:,OP_DZP)*h(:,OP_DR) - g(:,OP_DRP)*h(:,OP_DZ) & 
!                 + g(:,OP_DZ)*h(:,OP_DRP) - g(:,OP_DR)*h(:,OP_DZP)        

!         temp = + intx3(tempa,temp79a,r_79) &
!                + intx3(tempb,temp79b,r_79) &
!                - 2*intx4(e(:,:,OP_DP),f(:,OP_DP),temp79c,ri_79) &
!                - intx4(e(:,:,OP_DP),f(:,OP_1),temp79d,ri_79) 
!      end if

!   v2vpsif = temp
! end function v2vpsif

function v2vpsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vpsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f') + nu*f'_LP
           tempa = mu(h(:,OP_LP),OP_1) &
              + mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR) 
           ! [nu, psi]*R
           tempb = &
                mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR) 

        ! [v, psi]*R
        tempa2 = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)            
        ! (v, f') + v*f'_LP
        tempb2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) & 
               + mu(h(:,OP_LP),OP_1)
        ! [psi,f']*R
        temp79c = g(:,OP_DZ)*h(:,OP_DR) - g(:,OP_DR)*h(:,OP_DZ)            
        ! [psi,f']'*R
        temp79d = g(:,OP_DZP)*h(:,OP_DR) - g(:,OP_DRP)*h(:,OP_DZ) & 
                + g(:,OP_DZ)*h(:,OP_DRP) - g(:,OP_DR)*h(:,OP_DZP)        

        temp = prod(tempa,tempa2*r_79) &
             + prod(tempb,tempb2*r_79) &
             + prod(-2*temp79c*ri_79,OP_DP,OP_DP) &
             + prod(-temp79d*ri_79,OP_DP,OP_1) 
     end if

  v2vpsif = temp
end function v2vpsif

! V2vbf
! =====
! function v2vbf(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2vbf
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp

!      if(surface_int) then
!         temp = 0.
!      else
!         ! (F,f')
!         temp79a = g(:,OP_DZ)*h(:,OP_DZ) + g(:,OP_DR)*h(:,OP_DR)
!         ! (F,f'')
!         temp79b = g(:,OP_DZ)*h(:,OP_DZP) + g(:,OP_DR)*h(:,OP_DRP)        

!         temp = intx3(e(:,:,OP_1),f(:,OP_DP),temp79a) &
!              + intx3(e(:,:,OP_1),f(:,OP_1),temp79b) 
!      end if

!   v2vbf = temp
! end function v2vbf

function v2vbf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vbf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp

     if(surface_int) then
        temp%len = 0
     else
        ! (F,f')
        temp79a = g(:,OP_DZ)*h(:,OP_DZ) + g(:,OP_DR)*h(:,OP_DR)
        ! (F,f'')
        temp79b = g(:,OP_DZ)*h(:,OP_DZP) + g(:,OP_DR)*h(:,OP_DRP)        

        temp = prod(temp79a,OP_1,OP_DP) &
             + prod(temp79b,OP_1,OP_1) 
     end if

  v2vbf = temp
end function v2vbf

! V2vff
! =====
! function v2vff(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2vff
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (nu, f') + nu*f'_LP
!            tempa(j,:) = e(j,:,OP_1)*g(:,OP_LP) &
!               + e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR) 
!         end do
!         ! (v, f') + v*f'_LP
!         temp79a = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
!                 + f(:,OP_1)*h(:,OP_LP)
!         ! (f',f')
!         temp79b = g(:,OP_DZ)*h(:,OP_DZ) + g(:,OP_DR)*h(:,OP_DR)            
!         ! (f',f'') 
!         temp79c = g(:,OP_DZ)*h(:,OP_DZP) + g(:,OP_DR)*h(:,OP_DRP) 

!         temp = - intx3(tempa,temp79a,r2_79) &
!                - intx3(e(:,:,OP_DP),f(:,OP_DP),temp79b) &
!                - intx3(e(:,:,OP_DP),f(:,OP_1),temp79c) 
!      end if

!   v2vff = temp
! end function v2vff

function v2vff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempa2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f') + nu*f'_LP
           tempa = mu(g(:,OP_LP),OP_1) &
              + mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR) 
        ! (v, f') + v*f'_LP
        tempa2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) & 
               + mu(h(:,OP_LP),OP_1)
        ! (f',f')
        temp79b = g(:,OP_DZ)*h(:,OP_DZ) + g(:,OP_DR)*h(:,OP_DR)            
        ! (f',f'') 
        temp79c = g(:,OP_DZ)*h(:,OP_DZP) + g(:,OP_DR)*h(:,OP_DRP) 

        temp = prod(tempa,tempa2*(-r2_79)) &
             + prod(-temp79b,OP_DP,OP_DP) &
             + prod(-temp79c,OP_DP,OP_1) 
     end if

  v2vff = temp
end function v2vff
#endif

! V2upsipsi
! =========
! function v2upsipsi(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2upsipsi
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

! #if defined(USE3D) || defined(USECOMPLEX)  
!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element, MAX_PTS) :: tempc, tempd
!   integer :: j

!      temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
!      temp79b = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
!           +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)

!      if(surface_int) then
!         temp = &
!              - intx5(e(:,:,OP_1),ri_79,temp79a,norm79(:,1),h(:,OP_DRP)) &
!              - intx5(e(:,:,OP_1),ri_79,temp79a,norm79(:,2),h(:,OP_DZP)) &
!              - intx5(e(:,:,OP_1),ri_79,temp79b,norm79(:,1),h(:,OP_DR)) &
!              - intx5(e(:,:,OP_1),ri_79,temp79b,norm79(:,2),h(:,OP_DZ))
!      else
!         do j=1, dofs_per_element
! #ifdef USEST
!            tempc(j,:) =  &
! #else
!            tempc(j,:) = e(j,:,OP_1)*h(:,OP_GS) &
! #endif
!                 +    e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)

! #ifdef USEST
!            tempd(j,:) = - e(j,:,OP_DP)*h(:,OP_GS) &
! #else
!            tempd(j,:) = e(j,:,OP_1)*h(:,OP_GSP) &
! #endif
!                 +    e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP)
!         end do
!         temp = intx3(tempd,ri_79,temp79a) &
!              + intx3(tempc,ri_79,temp79b)
!      end if

!   v2upsipsi = temp
! #else
!   v2upsipsi = 0.
! #endif
! end function v2upsipsi

function v2upsipsi(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2upsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)  
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd

     tempa = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
     tempb = mu(g(:,OP_DR ),OP_DZP) + mu(-g(:,OP_DZ ),OP_DRP) &
          +    mu(g(:,OP_DRP),OP_DZ)  + mu(-g(:,OP_DZP),OP_DR)

     if(surface_int) then
        temp = &
               prod(mu(-ri_79*norm79(:,1)*h(:,OP_DRP),OP_1),tempa) &
             + prod(mu(-ri_79*norm79(:,2)*h(:,OP_DZP),OP_1),tempa) &
             + prod(mu(-ri_79*norm79(:,1)*h(:,OP_DR),OP_1),tempb) &
             + prod(mu(-ri_79*norm79(:,2)*h(:,OP_DZ),OP_1),tempb)
      else
#ifdef USEST
           tempc =  &
#else
           tempc = mu(h(:,OP_GS),OP_1) + &
#endif
                   mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
#ifdef USEST
           tempd = mu(-h(:,OP_GS),OP_DP) &
#else
           tempd = mu(h(:,OP_GSP),OP_1) &
#endif
                +  mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR)
        temp = prod(tempd*ri_79,tempa) &
             + prod(tempc*ri_79,tempb)
     end if
  v2upsipsi = temp
#else
  v2upsipsi%len = 0
#endif
end function v2upsipsi



! V2upsib
! =======
! function v2upsib(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2upsib
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element, MAX_PTS) :: tempa
!   integer :: j

!   do j=1, dofs_per_element
!      tempa(j,:) = e(j,:,OP_DZ)*f(:,OP_DR) - e(j,:,OP_DR)*f(:,OP_DZ)
!   end do

!      if(surface_int) then
!         temp79a = h(:,OP_DZ)*g(:,OP_DR) - h(:,OP_DR)*g(:,OP_DZ)
!         temp = intx4(e(:,:,OP_1),temp79a,norm79(:,1),f(:,OP_DZ)) &
!              - intx4(e(:,:,OP_1),temp79a,norm79(:,2),f(:,OP_DR))
!      else
!         temp = (intx3(tempa,g(:,OP_DR),h(:,OP_DZ)) &
!              -  intx3(tempa,g(:,OP_DZ),h(:,OP_DR)))

! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         temp79b = &
!               h(:,OP_DP)*(f(:,OP_DZ )*g(:,OP_DZ ) + f(:,OP_DR )*g(:,OP_DR )) &
!           + h(:,OP_1 )*(f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR )) 
!         temp = temp + intx3(e(:,:,OP_DP),ri2_79,temp79b)
! #else
!         temp79b = &
!           2.*(f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ))*h(:,OP_DP ) &
!           +  (f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP))*h(:,OP_DP ) &
!           +  (f(:,OP_DZ )*g(:,OP_DZ ) + f(:,OP_DR )*g(:,OP_DR ))*h(:,OP_DPP) &
!           +  (f(:,OP_DZPP)*g(:,OP_DZ ) + f(:,OP_DRPP)*g(:,OP_DR ))*h(:,OP_1) &
!           +  (f(:,OP_DZP )*g(:,OP_DZP) + f(:,OP_DRP )*g(:,OP_DRP))*h(:,OP_1)
!         temp = temp - intx3(e(:,:,OP_1),ri2_79,temp79b)
! #endif
! #endif
!      end if

!   v2upsib = temp
! end function v2upsib

function v2upsib(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2upsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp, tempa
  type(muarray) :: tempb

  temp79a(:)=1.
  tempa = prod(temp79a,OP_DZ,OP_DR) + prod(-temp79a,OP_DR,OP_DZ)

     if(surface_int) then
        temp79a = h(:,OP_DZ)*g(:,OP_DR) - h(:,OP_DR)*g(:,OP_DZ)
        temp = prod( temp79a*norm79(:,1),OP_1,OP_DZ) &
             + prod(-temp79a*norm79(:,2),OP_1,OP_DR)
      else
        temp = tempa*( g(:,OP_DR)*h(:,OP_DZ)) &
             + tempa*(-g(:,OP_DZ)*h(:,OP_DR))

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
        tempb = &
              (mu(g(:,OP_DZ ),OP_DZ) + mu(g(:,OP_DR ),OP_DR))*h(:,OP_DP) &
          +   (mu(g(:,OP_DZ ),OP_DZP)+ mu(g(:,OP_DR ),OP_DRP))*h(:,OP_1) 
        temp = temp + prod(mu(ri2_79,OP_DP),tempb)
#else
         tempb = &
             (mu(g(:,OP_DZ ),OP_DZP) +  mu(g(:,OP_DR ),OP_DRP))*(2.*h(:,OP_DP )) &
          +  (mu(g(:,OP_DZP),OP_DZ)  +  mu(g(:,OP_DRP),OP_DR))*h(:,OP_DP ) &
          +  (mu(g(:,OP_DZ ),OP_DZ)  +  mu(g(:,OP_DR ),OP_DR))*h(:,OP_DPP) &
          +  (mu(g(:,OP_DZ ),OP_DZPP) + mu(g(:,OP_DR ),OP_DRPP))*h(:,OP_1) &
          +  (mu(g(:,OP_DZP),OP_DZP) +  mu(g(:,OP_DRP),OP_DRP))*h(:,OP_1)
         temp = temp + prod(mu(-ri2_79,OP_1),tempb)
#endif
#endif
     end if

  v2upsib = temp
end function v2upsib

! V2ubb
! =====
!function v2ubb(e,f,g,h)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v2ubb
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!#if defined(USE3D) || defined(USECOMPLEX)
  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !temp = 0.
     !else
        !temp79a = h(:,OP_DP)*(f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)) &
             !+    h(:,OP_1)*(f(:,OP_DRP)*g(:,OP_DZ) - f(:,OP_DZP)*g(:,OP_DR))
        !temp = intx3(e(:,:,OP_1),ri_79,temp79a)
     !end if

  !v2ubb = temp
!#else
  !v2ubb = 0.
!#endif
!end function v2ubb

function v2ubb(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2ubb
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

     if(surface_int) then
        temp%len = 0
     else
        tempa = mu(h(:,OP_DP)*g(:,OP_DZ),OP_DR) + mu(-h(:,OP_DP)*g(:,OP_DR),OP_DZ) &
             +  mu(h(:,OP_1)*g(:,OP_DZ),OP_DRP) + mu(-h(:,OP_1)*g(:,OP_DR),OP_DZP)
        temp = prod(mu(ri_79,OP_1),tempa)
     end if
  v2ubb = temp
#else
  v2ubb%len = 0
#endif
end function v2ubb

#ifdef USE3D
! V2upsif
! =====
! function v2upsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2upsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp


!      if(surface_int) then
!         temp = 0.
!      else
!         ! ([u, psi]*R^2)_R/R
!         temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
!         if(itor.eq.1) then
!            temp79a = temp79a + ri_79* &
!                   (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
!         end if 
!         ! ([u, psi]*R^2)_Z/R
!         temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            
!         ! ((u, f')*R^2)_R/R^2
!         temp79c = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
!         if(itor.eq.1) then
!            temp79c = temp79c + 2*ri_79* &
!                   (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ)) 
!         end if 
!         ! ((u, f')*R^2)_Z/R^2
!         temp79d = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            

!         temp = - intx3(e(:,:,OP_DP),h(:,OP_DZ),temp79a) &
!                + intx3(e(:,:,OP_DP),h(:,OP_DR),temp79b) &
!                - intx3(e(:,:,OP_DP),g(:,OP_DR),temp79c) &
!                - intx3(e(:,:,OP_DP),g(:,OP_DZ),temp79d) 
!      end if

!   v2upsif = temp
! end function v2upsif

function v2upsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2upsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd


     if(surface_int) then
        temp%len = 0
     else
        ! ([u, psi]*R^2)_R/R
        tempa = mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &           
              + mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then
           tempa = tempa + &
                  (mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR))*ri_79 
        end if 
        ! ([u, psi]*R^2)_Z/R
        tempb = mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &           
              + mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR)            
        ! ((u, f')*R^2)_R/R^2
        tempc = mu(h(:,OP_DR),OP_DRR) + mu(h(:,OP_DZ),OP_DRZ) &           
              + mu(h(:,OP_DRR),OP_DR) + mu(h(:,OP_DRZ),OP_DZ)            
        if(itor.eq.1) then
           tempc = tempc + &
                  (mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ))*(2*ri_79) 
        end if 
        ! ((u, f')*R^2)_Z/R^2
        tempd = mu(h(:,OP_DR),OP_DRZ) + mu(h(:,OP_DZ),OP_DZZ) &           
              + mu(h(:,OP_DRZ),OP_DR) + mu(h(:,OP_DZZ),OP_DZ)            

        temp = prod(mu(-h(:,OP_DZ),OP_DP),tempa) &
             + prod(mu( h(:,OP_DR),OP_DP),tempb) &
             + prod(mu(-g(:,OP_DR),OP_DP),tempc) &
             + prod(mu(-g(:,OP_DZ),OP_DP),tempd) 
     end if

  v2upsif = temp
end function v2upsif

! V2ubf
! =====
!function v2ubf(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v2ubf
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa
!  integer :: j


!     if(surface_int) then
!        temp = 0.
!     else
!        do j=1, dofs_per_element
!           ! (nu, f') + nu*f'_LP
!           tempa(j,:) = e(j,:,OP_1)*h(:,OP_LP) &
!              + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
!        end do
!        ! (F'[u,f'] + F[u',f'])*R 
!        temp79a = g(:,OP_1)* &
!                  (h(:,OP_DR)*f(:,OP_DZP) - h(:,OP_DZ)*f(:,OP_DRP)) &
!                + g(:,OP_DP)* &
!                  (h(:,OP_DR)*f(:,OP_DZ) - h(:,OP_DZ)*f(:,OP_DR))
!        ! (R^2(u, f'))_R/R^2
!        temp79b = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
!        if(itor.eq.1) then
!           temp79b = temp79b + 2*ri_79* &
!                 (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
!        end if 
!        ! (R^2(u, f'))_Z/R^2
!        temp79c = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
!        ![u,F]*R
!        temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ) 

!        temp = intx3(e(:,:,OP_DP),temp79a,ri_79) &
!             + intx4(e(:,:,OP_1),temp79c,g(:,OP_DR),r_79) &
!             - intx4(e(:,:,OP_1),temp79b,g(:,OP_DZ),r_79) &
!             + intx3(tempa,temp79d,r_79)
!     end if

!  v2ubf = temp
!end function v2ubf

function v2ubf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2ubf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f') + nu*f'_LP
           tempa = mu(h(:,OP_LP),OP_1) &
              + mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR) 
        ! (F'[u,f'] + F[u',f'])*R 
        tempa2 = (mu(h(:,OP_DR),OP_DZP) + mu(-h(:,OP_DZ),OP_DRP))*g(:,OP_1) &
               + (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*g(:,OP_DP)
        ! (R^2(u, f'))_R/R^2
        tempb2 = mu(h(:,OP_DR),OP_DRR) + mu(h(:,OP_DZ),OP_DRZ) &           
              + mu(h(:,OP_DRR),OP_DR) + mu(h(:,OP_DRZ),OP_DZ)            
        if(itor.eq.1) then
           tempb2 = tempb2 + &
                 (mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ))*(2*ri_79)
        end if 
        ! (R^2(u, f'))_Z/R^2
        tempc2 = mu(h(:,OP_DR),OP_DRZ) + mu(h(:,OP_DZ),OP_DZZ) &           
              + mu(h(:,OP_DRZ),OP_DR) + mu(h(:,OP_DZZ),OP_DZ)            
        ![u,F]*R
        tempd2 = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR) 

        temp = prod(mu(ri_79,OP_DP),tempa2) &
             + prod(mu( g(:,OP_DR)*r_79,OP_1),tempc2) &
             + prod(mu(-g(:,OP_DZ)*r_79,OP_1),tempb2) &
             + prod(tempa,tempd2)*r_79
     end if

  v2ubf = temp
end function v2ubf

! V2uff
! =====
!function v2uff(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v2uff
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp

!     if(surface_int) then
!        temp = 0.
!     else
!        !((u, f')*R^2)_R/R^2
!        temp79a = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &           
!                + f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ)            
!        if(itor.eq.1) then
!           temp79a = temp79a + 2*ri_79* &
!                 (f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ))
!        end if 
!        !((u, f')*R^2)_Z/R^2
!        temp79b = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &           
!                + f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ)            

!        temp = intx4(e(:,:,OP_DP),h(:,OP_DZ),temp79a,r_79) &
!             - intx4(e(:,:,OP_DP),h(:,OP_DR),temp79b,r_79) 
!     end if

!  v2uff = temp
!end function v2uff

function v2uff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2uff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

     if(surface_int) then
        temp%len = 0
     else
        !((u, f')*R^2)_R/R^2
        tempa = mu(g(:,OP_DR),OP_DRR) + mu(g(:,OP_DZ),OP_DRZ) &           
              + mu(g(:,OP_DRR),OP_DR) + mu(g(:,OP_DRZ),OP_DZ)            
        if(itor.eq.1) then
           tempa = tempa + &
                 (mu(g(:,OP_DR),OP_DR) + mu(g(:,OP_DZ),OP_DZ))*(2*ri_79)
        end if 
        !((u, f')*R^2)_Z/R^2
        tempb = mu(g(:,OP_DR),OP_DRZ) + mu(g(:,OP_DZ),OP_DZZ) &           
              + mu(g(:,OP_DRZ),OP_DR) + mu(g(:,OP_DZZ),OP_DZ)            

        temp = prod(mu( h(:,OP_DZ)*r_79,OP_DP),tempa) &
             + prod(mu(-h(:,OP_DR)*r_79,OP_DP),tempb) 
     end if

  v2uff = temp
end function v2uff
#endif

! v2upsisb2
! ========
function v2upsisb2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2upsisb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  v2upsisb2 = 0.
end function v2upsisb2


! v2ubsb1
! =======
function v2ubsb1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2ubsb1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  v2ubsb1 = 0.
end function v2ubsb1


! v2chipsipsi
! ===========
! function v2chipsipsi(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2chipsipsi
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

! #if defined(USE3D) || defined(USECOMPLEX)
!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempc, tempd
!   integer :: j

!      temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
!      temp79b = f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
!           +    f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)

!      if(surface_int) then
!         temp = intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_DRP)) &
!              + intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,2),h(:,OP_DZP)) &
!              + intx5(e(:,:,OP_1),ri4_79,temp79b,norm79(:,1),h(:,OP_DR )) &
!              + intx5(e(:,:,OP_1),ri4_79,temp79b,norm79(:,2),h(:,OP_DZ ))
!      else
!         do j=1, dofs_per_element
! #ifdef USEST
!            tempc(j,:) =  &
! #else
!            tempc(j,:) = e(j,:,OP_1)*h(:,OP_GS) &
! #endif
!              + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
! #ifdef USEST
!            tempd(j,:) = -e(j,:,OP_DP)*h(:,OP_GS) &
! #else
!            tempd(j,:) = e(j,:,OP_1)*h(:,OP_GSP) &
! #endif
!              + e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP)
!         end do
        
!         temp = -intx3(tempd,ri4_79,temp79a) &
!                -intx3(tempc,ri4_79,temp79b)
!      end if

!   v2chipsipsi = temp
! #else
!   v2chipsipsi = 0.
! #endif
! end function v2chipsipsi

function v2chipsipsi(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chipsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd

     tempa = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
     tempb = mu(g(:,OP_DZ ),OP_DZP) + mu(g(:,OP_DR ),OP_DRP) &
          +  mu(g(:,OP_DZP),OP_DZ) +  mu(g(:,OP_DRP),OP_DR)


     if(surface_int) then
        temp = prod(mu(ri4_79*norm79(:,1)*h(:,OP_DRP),OP_1),tempa) &
             + prod(mu(ri4_79*norm79(:,2)*h(:,OP_DZP),OP_1),tempa) &
             + prod(mu(ri4_79*norm79(:,1)*h(:,OP_DR ),OP_1),tempb) &
             + prod(mu(ri4_79*norm79(:,2)*h(:,OP_DZ ),OP_1),tempb)
       else
#ifdef USEST
           tempc = &
#else
           tempc = mu(h(:,OP_GS),OP_1) + &
#endif
               mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
#ifdef USEST
           tempd = mu(-h(:,OP_GS),OP_DP) &
#else
           tempd = mu(h(:,OP_GSP),OP_1) &
#endif
             + mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR)
         
        temp =  prod(tempd*(-ri4_79),tempa) &
               +prod(tempc*(-ri4_79),tempb)
     end if
  v2chipsipsi = temp
#else
  v2chipsipsi%len = 0
#endif
end function v2chipsipsi


! v2chipsib
! =========
! function v2chipsib(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2chipsib
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempb, tempc
!   integer :: j


!   temp79a = h(:,OP_1 )*f(:,OP_GS) &
!        +    h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)

!   do j=1, dofs_per_element
!      tempb(j,:) = e(j,:,OP_DR)*h(:,OP_DZ) - e(j,:,OP_DZ)*h(:,OP_DR)

!      tempc(j,:) = e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)
!   end do
  

!      if(surface_int) then
!         temp79a = norm79(:,1)*f(:,OP_DR) + norm79(:,2)*f(:,OP_DZ)
!         temp79b = h(:,OP_1)* &
!              (norm79(:,1)*g(:,OP_DZ) - norm79(:,2)*g(:,OP_DR))

!         temp = intx5(e(:,:,OP_1),ri3_79,temp79a,g(:,OP_DZ),h(:,OP_DR)) &
!              - intx5(e(:,:,OP_1),ri3_79,temp79a,g(:,OP_DR),h(:,OP_DZ)) &
!              + intx4(e(:,:,OP_1),ri3_79,temp79b,f(:,OP_GS))
!         if(itor.eq.1) then
!            temp = temp - 2.*intx4(e(:,:,OP_1),ri4_79,temp79b,f(:,OP_DR))
!         endif
!      else
!         temp = intx3(tempc,ri3_79,temp79a) &
!              + intx4(tempb,ri3_79,f(:,OP_DZ),g(:,OP_DZ)) &
!              + intx4(tempb,ri3_79,f(:,OP_DR),g(:,OP_DR))

!         if(itor.eq.1) then
!            temp = temp - &
!                 2.*intx4(tempc,ri4_79,f(:,OP_DR),h(:,OP_1))
!         endif
! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         temp79d = &
!           (f(:,OP_DZ)*g(:,OP_DR ) - f(:,OP_DR)*g(:,OP_DZ ))*h(:,OP_DP ) &
!         + (f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ))*h(:,OP_1 ) 
!         temp = temp + intx3(e(:,:,OP_DP),ri5_79,temp79d)
! #else
!         temp79d = &
!          2.*(f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ))*h(:,OP_DP ) &
!          +  (f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP))*h(:,OP_DP ) &
!          +  (f(:,OP_DZ )*g(:,OP_DR ) - f(:,OP_DR )*g(:,OP_DZ ))*h(:,OP_DPP) &
!          +  (f(:,OP_DZPP)*g(:,OP_DR ) - f(:,OP_DRPP)*g(:,OP_DZ ))*h(:,OP_1) &
!          +  (f(:,OP_DZP )*g(:,OP_DRP) - f(:,OP_DRP )*g(:,OP_DZP))*h(:,OP_1)
!         temp = temp - intx3(e(:,:,OP_1),ri5_79,temp79d)
! #endif
! #endif
!      end if

!   v2chipsib = temp
! end function v2chipsib

function v2chipsib(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chipsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd
  integer :: j


  tempa = mu(h(:,OP_1 ),OP_GS) &
       +  mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)

     tempb = mu(h(:,OP_DZ),OP_DR) + mu(-h(:,OP_DR),OP_DZ)
     tempc = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
   
     if(surface_int) then
        tempa = mu(norm79(:,1),OP_DR) + mu(norm79(:,2),OP_DZ)
        temp79b = h(:,OP_1)* &
             (norm79(:,1)*g(:,OP_DZ) - norm79(:,2)*g(:,OP_DR))

        temp = prod(mu( ri3_79*g(:,OP_DZ)*h(:,OP_DR),OP_1),tempa) &
             + prod(mu(-ri3_79*g(:,OP_DR)*h(:,OP_DZ),OP_1),tempa) &
             + prod(ri3_79*temp79b,OP_1,OP_GS)
        if(itor.eq.1) then
           temp = temp + prod(-2.*ri4_79*temp79b,OP_1,OP_DR)
        endif
     else
        temp = prod(tempc*ri3_79,tempa) &
             + prod(tempb,mu(ri3_79*g(:,OP_DZ),OP_DZ)) &
             + prod(tempb,mu(ri3_79*g(:,OP_DR),OP_DR))

        if(itor.eq.1) then
           temp = temp + &
                prod(tempc,mu(-2.*ri4_79*h(:,OP_1),OP_DR))
        endif
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
        tempd = &
            (mu(g(:,OP_DR ),OP_DZ) + mu(-g(:,OP_DZ ),OP_DR))*h(:,OP_DP ) &
        +   (mu(g(:,OP_DR ),OP_DZP)+ mu(-g(:,OP_DZ ),OP_DRP))*h(:,OP_1 ) 
        temp = temp + prod(mu(ri5_79,OP_DP),tempd)
#else
         tempd = &
            (mu(g(:,OP_DR ),OP_DZP) + mu(-g(:,OP_DZ ),OP_DRP))*(2.*h(:,OP_DP )) &
         +  (mu(g(:,OP_DRP),OP_DZ)  + mu(-g(:,OP_DZP),OP_DR))*h(:,OP_DP ) &
         +  (mu(g(:,OP_DR ),OP_DZ)  + mu(-g(:,OP_DZ ),OP_DR))*h(:,OP_DPP) &
         +  (mu(g(:,OP_DR ),OP_DZPP)+ mu(-g(:,OP_DZ ),OP_DRPP))*h(:,OP_1) &
         +  (mu(g(:,OP_DRP),OP_DZP) + mu(-g(:,OP_DZP),OP_DRP))*h(:,OP_1)
         temp = temp + prod(mu(-ri5_79,OP_1),tempd)
#endif
#endif
     end if

  v2chipsib = temp
end function v2chipsib

! v2chibb
! =======
!function v2chibb(e,f,g,h)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v2chibb
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!#if defined(USE3D) || defined(USECOMPLEX)
  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !temp = 0.
     !else
        !temp = intx5(e(:,:,OP_1),ri4_79,g(:,OP_1),h(:,OP_DR),f(:,OP_DRP)) &
             !+ intx5(e(:,:,OP_1),ri4_79,g(:,OP_1),h(:,OP_DZ),f(:,OP_DZP)) &
             !+ intx5(e(:,:,OP_1),ri4_79,g(:,OP_DP),h(:,OP_DR),f(:,OP_DR)) &
             !+ intx5(e(:,:,OP_1),ri4_79,g(:,OP_DP),h(:,OP_DZ),f(:,OP_DZ))
     !end if

  !v2chibb = temp
!#else
  !v2chibb = 0.
!#endif
!end function v2chibb

function v2chibb(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chibb
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

     if(surface_int) then
        temp%len = 0
     else
        temp = prod(ri4_79*g(:,OP_1)*h(:,OP_DR),OP_1,OP_DRP) &
             + prod(ri4_79*g(:,OP_1)*h(:,OP_DZ),OP_1,OP_DZP) &
             + prod(ri4_79*g(:,OP_DP)*h(:,OP_DR),OP_1,OP_DR) &
             + prod(ri4_79*g(:,OP_DP)*h(:,OP_DZ),OP_1,OP_DZ)
      end if

  v2chibb = temp
#else
  v2chibb%len = 0
#endif
end function v2chibb

#ifdef USE3D
! V2chipsif
! =====
! function v2chipsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2chipsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp

!      if(surface_int) then
!         temp = 0.
!      else
!         ! [chi, f']_R*R 
!         temp79a = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
!         if(itor.eq.1) then 
!            temp79a = temp79a - ri_79* &        
!                 (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
!         end if 
!         ! [chi, f']_Z*R 
!         temp79b = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
!         ! ((chi, psi)/R^2)_R * R^2
!         temp79c = f(:,OP_DRZ)*g(:,OP_DZ) + f(:,OP_DRR)*g(:,OP_DR) &
!                 + f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DR)*g(:,OP_DRR) 
!         if(itor.eq.1) then 
!            temp79c = temp79c - 2*ri_79*  & 
!                  (f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)) 
!         end if 
!         ! ((chi, psi)/R^2)_Z * R^2
!         temp79d = f(:,OP_DZZ)*g(:,OP_DZ) + f(:,OP_DRZ)*g(:,OP_DR) &
!                 + f(:,OP_DZ)*g(:,OP_DZZ) + f(:,OP_DR)*g(:,OP_DRZ) 

!         temp = &
!              - intx4(e(:,:,OP_DP),temp79a,ri3_79,g(:,OP_DR)) &
!              - intx4(e(:,:,OP_DP),temp79b,ri3_79,g(:,OP_DZ)) &
!              - intx4(e(:,:,OP_DP),temp79d,ri3_79,h(:,OP_DR)) &
!              + intx4(e(:,:,OP_DP),temp79c,ri3_79,h(:,OP_DZ)) 
!      end if

!   v2chipsif = temp
! end function v2chipsif

function v2chipsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chipsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd

     if(surface_int) then
        temp%len = 0
     else
        ! [chi, f']_R*R 
        tempa = mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR) &           
              + mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then 
           tempa = tempa + &        
                (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79) 
        end if 
        ! [chi, f']_Z*R 
        tempb = mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ) &           
              + mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR)            
        ! ((chi, psi)/R^2)_R * R^2
        tempc = mu(g(:,OP_DZ),OP_DRZ) + mu(g(:,OP_DR),OP_DRR) &
              + mu(g(:,OP_DRZ),OP_DZ) + mu(g(:,OP_DRR),OP_DR) 
        if(itor.eq.1) then 
           tempc = tempc + & 
                 (mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR))*(-2*ri_79) 
        end if 
        ! ((chi, psi)/R^2)_Z * R^2
        tempd = mu(g(:,OP_DZ),OP_DZZ) + mu(g(:,OP_DR),OP_DRZ) &
              + mu(g(:,OP_DZZ),OP_DZ) + mu(g(:,OP_DRZ),OP_DR) 

        temp = &
               prod(mu(-ri3_79*g(:,OP_DR),OP_DP),tempa) &
             + prod(mu(-ri3_79*g(:,OP_DZ),OP_DP),tempb) &
             + prod(mu(-ri3_79*h(:,OP_DR),OP_DP),tempd) &
             + prod(mu( ri3_79*h(:,OP_DZ),OP_DP),tempc) 
     end if

  v2chipsif = temp
end function v2chipsif

! V2chibf
! =====
! function v2chibf(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2chibf
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (nu, f') + nu*f'_LP
!            tempa(j,:) = e(j,:,OP_1)*h(:,OP_LP) &
!               + e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
!            ! [nu, F]*R
!            tempb(j,:) = &
!               + e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
!         end do

!         ! [chi, f']*R
!         temp79a = f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)  
!         ! (chi, f')
!         temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) 
!         ! (chi', f')
!         temp79c = f(:,OP_DRP)*h(:,OP_DR) + f(:,OP_DZP)*h(:,OP_DZ) 
!         ! (chi, F/R^4)*R^4 + F*chi_LP
!         temp79d = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
!                 + f(:,OP_LP)*g(:,OP_1) 
!         if(itor.eq.1) then 
!            temp79d = temp79d - 4*ri_79*f(:,OP_DR)*g(:,OP_1) 
!         end if 

!         temp = - intx3(tempb,temp79a,ri2_79) &
!                - intx3(tempa,temp79d,ri2_79) &
!                - intx4(e(:,:,OP_DP),temp79b,ri4_79,g(:,OP_DP)) &
!                - intx4(e(:,:,OP_DP),temp79c,ri4_79,g(:,OP_1)) 
!      end if

!   v2chibf = temp
! end function v2chibf

function v2chibf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chibf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f') + nu*f'_LP
           tempa = mu(h(:,OP_LP),OP_1) &
              + mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR) 
           ! [nu, F]*R
           tempb = &
                mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR) 

        ! [chi, f']*R
        tempa2 = mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR)  
        ! (chi, f')
        tempb2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) 
        ! (chi', f')
        tempc2 = mu(h(:,OP_DR),OP_DRP) + mu(h(:,OP_DZ),OP_DZP) 
        ! (chi, F/R^4)*R^4 + F*chi_LP
        tempd2 = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR) &
               + mu(g(:,OP_1),OP_LP) 
        if(itor.eq.1) then 
           tempd2 = tempd2 + mu(-4*ri_79*g(:,OP_1),OP_DR) 
        end if 

        temp = prod(tempb,tempa2*(-ri2_79)) &
             + prod(tempa,tempd2*(-ri2_79)) &
             + prod(mu(-ri4_79*g(:,OP_DP),OP_DP),tempb2) &
             + prod(mu(-ri4_79*g(:,OP_1),OP_DP),tempc2) 
     end if

  v2chibf = temp
end function v2chibf

! V2chiff
! =====
! function v2chiff(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v2chiff
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j

!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! [nu, f']*R
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)  
!            ! [nu, f'']*R
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP)  
!         end do

!         ! [chi, f']*R
!         temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)        
!         ! [chi, f']'*R
!         temp79b = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP)            

!         temp = - intx3(tempa,temp79b,ri2_79) &
!                - intx3(tempb,temp79a,ri2_79) 
!      end if

!   v2chiff = temp
! end function v2chiff

function v2chiff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2chiff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2
  integer :: j

     if(surface_int) then
        temp%len = 0
     else
           ! [nu, f']*R
           tempa = &
                mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR)  
           ! [nu, f'']*R
           tempb = &
                mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR)  

        ! [chi, f']*R
        tempa2 = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)        
        ! [chi, f']'*R
        tempb2 = mu(g(:,OP_DR),OP_DZP) + mu(-g(:,OP_DZ),OP_DRP) &           
                + mu(g(:,OP_DRP),OP_DZ) + mu(-g(:,OP_DZP),OP_DR)            

        temp = prod(tempa,tempb2*(-ri2_79)) &
             + prod(tempb,tempa2*(-ri2_79)) 
     end if

  v2chiff = temp
end function v2chiff
#endif

! v2vchin
! =======
function v2vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vchin = 0.
     return
  end if


     temp =-intx4(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          - intx4(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
     if(itor.eq.1) then
        temp = temp &
             - 2.*intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DR),h(:,OP_1))
     endif

  v2vchin = temp
end function v2vchin

function v2vchin1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vchin1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vchin1%len = 0
     return
  end if

     temp = prod(-g(:,OP_DZ)*h(:,OP_1),OP_1,OP_DZ) &
          + prod(-g(:,OP_DR)*h(:,OP_1),OP_1,OP_DR)
     if(itor.eq.1) then
        temp = temp &
             + prod(-2.*ri_79*g(:,OP_DR)*h(:,OP_1),OP_1,OP_1)
     endif

  v2vchin1 = temp
end function v2vchin1

function v2vchin2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vchin2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vchin2%len = 0
     return
  end if

     temp = prod(-f(:,OP_DZ)*h(:,OP_1),OP_1,OP_DZ) &
          + prod(-f(:,OP_DR)*h(:,OP_1),OP_1,OP_DR)
     if(itor.eq.1) then
        temp = temp &
             + prod(-2.*ri_79*f(:,OP_1)*h(:,OP_1),OP_1,OP_DR)
     endif

  v2vchin2 = temp
end function v2vchin2


! v2chibsb1
! =========
function v2chibsb1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2chibsb1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  v2chibsb1 = 0.
end function v2chibsb1

! v2psisb2
! ========
function v2psisb2(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psisb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
          - intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR))
  endif

  v2psisb2 = temp
end function v2psisb2

! v2bsb1
! ======
function v2bsb1(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2bsb1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR)) &
          - intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ))
  end if

  v2bsb1 = temp
end function v2bsb1


! V2vs
! ====
function v2vs(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2vs
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(idens.eq.0 .or. nosig.eq.1 .or. surface_int) then
     v2vs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  v2vs = -intx4(e(:,:,OP_1),r2_79,f(:,OP_1),temp79a)

end function v2vs

function v2vs1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2vs1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  if(idens.eq.0 .or. nosig.eq.1 .or. surface_int) then
     v2vs1%len = 0
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  v2vs1 = prod(-r2_79*temp79a,OP_1,OP_1)

end function v2vs1

! V2psif1
! =======
function v2psif1(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psif1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(surface_int) then
     v2psif1 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif1 = intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZP)) &
       -    intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DRP))
#else
  v2psif1 = 0.
#endif
end function v2psif1

function v2psif11(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psif11
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  if(surface_int) then
     v2psif11%len =0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif11 = prod( ri_79*g(:,OP_DZP),OP_1,OP_DR) &
           + prod(-ri_79*g(:,OP_DRP),OP_1,OP_DZ)
#else
  v2psif11%len =0
#endif
end function v2psif11

function v2psif12(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psif12
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  if(surface_int) then
     v2psif12%len =0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif12 = prod( ri_79*f(:,OP_DR),OP_1,OP_DZP) &
           + prod(-ri_79*f(:,OP_DZ),OP_1,OP_DRP)
#else
  v2psif12%len =0
#endif
end function v2psif12

! V2psif2
! =======
function v2psif2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2psif2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(surface_int) then
     v2psif2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif2 = intx4(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ)) &
       -    intx4(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR))
#else
  v2psif2 = 0.
#endif
end function v2psif2

function v2psif21(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psif21
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  if(surface_int) then
     v2psif21%len =0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif21 = prod( ri_79*g(:,OP_DZ),OP_1,OP_DRP) &
           + prod(-ri_79*g(:,OP_DR),OP_1,OP_DZP)
#else
  v2psif21%len =0
#endif
end function v2psif21

function v2psif22(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2psif22
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  if(surface_int) then
     v2psif22%len =0
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  v2psif22 = prod( ri_79*f(:,OP_DRP),OP_1,OP_DZ) &
           + prod(-ri_79*f(:,OP_DZP),OP_1,OP_DR)
#else
  v2psif22%len =0
#endif
end function v2psif22



! V2bf
! ====
function v2bf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2bf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  
  if(surface_int) then
     v2bf = 0.
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  v2bf = &
       - intx3(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
       - intx3(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR))
#else
  v2bf = 0.
#endif
end function v2bf

function v2bf1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2bf1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  
  if(surface_int) then
     v2bf1%len = 0
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  v2bf1 = &
          prod(-g(:,OP_DZ),OP_1,OP_DZ) &
        + prod(-g(:,OP_DR),OP_1,OP_DR)
#else
  v2bf1%len = 0
#endif
end function v2bf1

function v2bf2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2bf2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  
  if(surface_int) then
     v2bf2%len = 0
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  v2bf2 = &
          prod(-f(:,OP_DZ),OP_1,OP_DZ) &
        + prod(-f(:,OP_DR),OP_1,OP_DR)
#else
  v2bf2%len = 0
#endif
end function v2bf2

! V2ff
! ====
function v2ff(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2ff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  
  if(surface_int) then
     v2ff = 0.
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  v2ff = &
       - intx3(e(:,:,OP_1),f(:,OP_DZP),g(:,OP_DZ)) &
       - intx3(e(:,:,OP_1),f(:,OP_DRP),g(:,OP_DR))
#else
  v2ff = 0.
#endif
end function v2ff

function v2ff1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2ff1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  
  if(surface_int) then
     v2ff1%len = 0
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  v2ff1 = &
          prod(-g(:,OP_DZ),OP_1,OP_DZP) &
        + prod(-g(:,OP_DR),OP_1,OP_DRP)
#else
  v2ff1%len = 0
#endif
end function v2ff1

function v2ff2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v2ff2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  
  if(surface_int) then
     v2ff2%len = 0
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  v2ff2 = &
          prod(-f(:,OP_DZP),OP_1,OP_DZ) &
        + prod(-f(:,OP_DRP),OP_1,OP_DR)
#else
  v2ff2%len = 0
#endif
end function v2ff2


! V2be
! ===

function v2parpb2ipsipsi(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2parpb2ipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     v2parpb2ipsipsi = 0.
     return
  end if

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = -f(:,OP_DP)*g(:,OP_1)*ri2_79
  temp = intx4(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR))    &
       + intx4(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ)) 
#endif

  v2parpb2ipsipsi = temp
end function v2parpb2ipsipsi

function v2parpb2ibb(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2parpb2ibb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     v2parpb2ibb = 0.
     return
  end if

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = -f(:,OP_DP)*g(:,OP_1)*h(:,OP_1)*i(:,OP_1)*ri2_79
  temp = temp + intx2(e(:,:,OP_1),temp79a)  
#endif

  v2parpb2ibb = temp
end function v2parpb2ibb


function v2parpb2ipsib(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2parpb2ipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i

  if(surface_int) then
     v2parpb2ipsib = 0.
     return
  end if

  temp79a =  f(:,OP_1)*g(:,OP_1)*i(:,OP_1)*ri_79

  v2parpb2ipsib = intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DR))    &
       -          intx3(e(:,:,OP_DR),temp79a,h(:,OP_DZ)) 
end function v2parpb2ipsib

#ifdef USEPARTICLES
function v2p_2(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2p_2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  if(surface_int) then
     v2p_2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  v2p_2 = -intx3(e(:,:,OP_1),f(:,OP_DP),g)
#else
  v2p_2 = 0.
#endif

end function v2p_2

function v2pbb(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2pbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = -f(:,OP_DP)*g*bzt79(:,OP_1)*bzt79(:,OP_1)*ri2_79
  temp = temp + intx2(e(:,:,OP_1),temp79a)  
  temp79a = g*bzt79(:,OP_1)
  temp = temp + intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),bfpt79(:,OP_DR))    &
       +  intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),bfpt79(:,OP_DZ)) 
   #endif
  temp79a = -g*bzt79(:,OP_1)*ri_79
  temp = temp + intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),pst79(:,OP_DR))    &
       -          intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),pst79(:,OP_DZ)) 
    end if

  v2pbb = temp
end function v2pbb

function v2jxb(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v2jxb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = f*ri_79
        temp = temp                                   &
             +  intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DZ),bzt79(:,OP_DR)+bfpt79(:,OP_DRP))   &
             -  intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DR),bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))   
  temp79a = f*ri2_79
  temp = temp + intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DR),pst79(:,OP_DRP))    &
       + intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DZ),pst79(:,OP_DZP)) 
  temp79a = f
  temp = temp + intx4(e(:,:,OP_1),temp79a,bfpt79(:,OP_DR),(bzt79(:,OP_DR)+bfpt79(:,OP_DRP)))    &
       + intx4(e(:,:,OP_1),temp79a,bfpt79(:,OP_DZ),(bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))) 
  temp79a = f*ri_79
  temp = temp + intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DZP),bfpt79(:,OP_DR))    &
       - intx4(e(:,:,OP_1),temp79a,pst79(:,OP_DRP),bfpt79(:,OP_DZ)) 
#endif

    end if

  v2jxb = temp
end function v2jxb

#endif
!==============================================================================
! V3 TERMS
!==============================================================================

! V3chin
! ======
!function v3chin(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3chin
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !if(inonormalflow.eq.1) then
           !temp = 0.
        !else
           !temp = intx5(e(:,:,OP_1),ri4_79,g(:,OP_1),norm79(:,1),f(:,OP_DR)) &
                !+ intx5(e(:,:,OP_1),ri4_79,g(:,OP_1),norm79(:,2),f(:,OP_DZ))
        !end if
     !else
        !temp = - intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_1)) &
               !- intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_1))
     !end if

  !v3chin = temp
!end function v3chin

function v3chin(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chin
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           temp = prod(ri4_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DR) &
                + prod(ri4_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DZ)
         end if
     else
        temp =   prod(-ri4_79*g(:,OP_1),OP_DR,OP_DR) &
               + prod(-ri4_79*g(:,OP_1),OP_DZ,OP_DZ)
     end if

  v3chin = temp
end function v3chin



! V3chimu
! =======
! function v3chimu(e,f,g,h)

!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3chimu
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
!   vectype, dimension(dofs_per_element) :: temp

!      if(surface_int) then
!         temp = 2.* &
!              (intx5(e(:,:,OP_DR),ri4_79,g(:,OP_1),norm79(:,1),f(:,OP_DZZ)) &
!              -intx5(e(:,:,OP_DR),ri4_79,g(:,OP_1),norm79(:,2),f(:,OP_DRZ)) &
!              -intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_1),norm79(:,1),f(:,OP_DRZ)) &
!              +intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_1),norm79(:,2),f(:,OP_DRR)) &
!              -intx5(e(:,:,OP_DR),ri4_79,h(:,OP_1),norm79(:,1),f(:,OP_GS)) &
!              -intx5(e(:,:,OP_DZ),ri4_79,h(:,OP_1),norm79(:,2),f(:,OP_GS)))
!         if(itor.eq.1) then
!            temp = temp + 2.* &
!                 (intx5(e(:,:,OP_DZ),ri5_79,g(:,OP_1),norm79(:,1),f(:,OP_DZ)) &
!                 -intx5(e(:,:,OP_DZ),ri5_79,g(:,OP_1),norm79(:,2),f(:,OP_DR)) &
!                 +intx5(e(:,:,OP_DR),ri5_79,g(:,OP_1),norm79(:,1),f(:,OP_DR)) &
!                 +intx5(e(:,:,OP_DR),ri5_79,g(:,OP_1),norm79(:,2),f(:,OP_DZ)) &
!                 +intx5(e(:,:,OP_1),ri5_79,g(:,OP_1),norm79(:,1),f(:,OP_DZZ)) &
!                 -intx5(e(:,:,OP_1),ri5_79,g(:,OP_1),norm79(:,2),f(:,OP_DRZ)) &
!                 +2.*intx5(e(:,:,OP_1),ri6_79,g(:,OP_1),norm79(:,2),f(:,OP_DZ)))
!         endif

! #if defined(USE3D) || defined(USECOMPLEX)
!         temp = temp &
!              + intx5(e(:,:,OP_1),ri6_79,g(:,OP_1),norm79(:,1),f(:,OP_DRPP)) &
!              + intx5(e(:,:,OP_1),ri6_79,g(:,OP_1),norm79(:,2),f(:,OP_DZPP))
! #endif
!      else
!         temp79b = f(:,OP_DRR)
!         temp79d = f(:,OP_DRZ)
!         if(itor.eq.1) then
!            temp79b = temp79b - 2.*ri_79*f(:,OP_DR)
!            temp79d = temp79d -    ri_79*f(:,OP_DZ)
!         endif
!         temp = 2.* &
!              (intx4(e(:,:,OP_DZZ),ri4_79,f(:,OP_DZZ),g(:,OP_1)) &
!              +intx4(e(:,:,OP_DRR),ri4_79,temp79b,g(:,OP_1)) &
!              +2.*intx4(e(:,:,OP_DRZ),ri4_79,temp79d,g(:,OP_1)) & 
!              +intx4(e(:,:,OP_GS),ri4_79,f(:,OP_GS),h(:,OP_1)) &
!              -intx4(e(:,:,OP_GS),ri4_79,f(:,OP_GS),g(:,OP_1)))
!         if(itor.eq.1) then
!            temp = temp &
!                 + 2.*intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DR),g(:,OP_1)) &
!                 - 4.*intx4(e(:,:,OP_DR),ri5_79,temp79b,g(:,OP_1)) &
!                 - 4.*intx4(e(:,:,OP_DZ),ri5_79,temp79d,g(:,OP_1))
!         endif
! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         temp = temp + &
!              (intx4(e(:,:,OP_DZP),ri6_79,f(:,OP_DZP),g(:,OP_1)) &
!              +intx4(e(:,:,OP_DRP),ri6_79,f(:,OP_DRP),g(:,OP_1))) &
!              +(intx4(e(:,:,OP_DZ),ri6_79,f(:,OP_DZP),g(:,OP_DP)) &
!              +intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DRP),g(:,OP_DP)))
! #else
!         temp = temp - &
!              (intx4(e(:,:,OP_DZ),ri6_79,f(:,OP_DZPP),g(:,OP_1)) &
!              +intx4(e(:,:,OP_DR),ri6_79,f(:,OP_DRPP),g(:,OP_1)))
! #endif
! #endif
!      end if

!   v3chimu = temp
!   return
! end function v3chimu

function v3chimu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chimu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  type(muarray) :: tempb, tempd

     if(surface_int) then
        temp = &
             (prod( 2.*ri4_79*g(:,OP_1)*norm79(:,1),OP_DR,OP_DZZ) &
             +prod(-2.*ri4_79*g(:,OP_1)*norm79(:,2),OP_DR,OP_DRZ) &
             +prod(-2.*ri4_79*g(:,OP_1)*norm79(:,1),OP_DZ,OP_DRZ) &
             +prod( 2.*ri4_79*g(:,OP_1)*norm79(:,2),OP_DZ,OP_DRR) &
             +prod(-2.*ri4_79*h(:,OP_1)*norm79(:,1),OP_DR,OP_GS) &
             +prod(-2.*ri4_79*h(:,OP_1)*norm79(:,2),OP_DZ,OP_GS))
        if(itor.eq.1) then
           temp = temp + &
                (prod( 2.*ri5_79*g(:,OP_1)*norm79(:,1),OP_DZ,OP_DZ) &
                +prod(-2.*ri5_79*g(:,OP_1)*norm79(:,2),OP_DZ,OP_DR) &
                +prod( 2.*ri5_79*g(:,OP_1)*norm79(:,1),OP_DR,OP_DR) &
                +prod( 2.*ri5_79*g(:,OP_1)*norm79(:,2),OP_DR,OP_DZ) &
                +prod( 2.*ri5_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DZZ) &
                +prod(-2.*ri5_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DRZ) &
                +prod( 4.*ri6_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DZ))
         endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + prod(ri6_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DRPP) &
             + prod(ri6_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DZPP)
#endif
     else
        temp79b = 1.
        tempb = mu(temp79b,OP_DRR)
        tempd = mu(temp79b,OP_DRZ)
        if(itor.eq.1) then
           tempb = tempb + mu(-2.*ri_79,OP_DR)
           tempd = tempd + mu(-ri_79,OP_DZ)
        endif
        temp = &
             (prod(2.*ri4_79*g(:,OP_1),OP_DZZ,OP_DZZ) &
             +prod(mu(2.*ri4_79*g(:,OP_1),OP_DRR),tempb) &
             +prod(mu(4.*ri4_79*g(:,OP_1),OP_DRZ),tempd) & 
             +prod( 2.*ri4_79*h(:,OP_1),OP_GS,OP_GS) &
             +prod(-2.*ri4_79*g(:,OP_1),OP_GS,OP_GS))
        if(itor.eq.1) then
           temp = temp &
                + prod(2.*ri6_79*g(:,OP_1),OP_DR,OP_DR) &
                + prod(mu(-4.*ri5_79*g(:,OP_1),OP_DR),tempb) &
                + prod(mu(-4.*ri5_79*g(:,OP_1),OP_DZ),tempd)
        endif
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
        temp = temp + &
              prod(ri6_79*g(:,OP_1),OP_DZP,OP_DZP) &
             +prod(ri6_79*g(:,OP_1),OP_DRP,OP_DRP) &
             +prod(ri6_79*g(:,OP_DP),OP_DZ,OP_DZP) &
             +prod(ri6_79*g(:,OP_DP),OP_DR,OP_DRP)
#else
         temp = temp + &
             (prod(-ri6_79*g(:,OP_1),OP_DZ,OP_DZPP) &
             +prod(-ri6_79*g(:,OP_1),OP_DR,OP_DRPP))
#endif
#endif
     end if

  v3chimu = temp
  return
end function v3chimu


! V3umu
! =====
! function v3umu(e,f,g,h)

!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3umu
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
!   vectype, dimension(dofs_per_element) :: temp

!      temp79c = f(:,OP_DZZ) - f(:,OP_DRR)
!      if(itor.eq.1) temp79c = temp79c - ri_79*f(:,OP_DR)

!      if(surface_int) then
!         temp = intx5(e(:,:,OP_DZ),ri_79,g(:,OP_1),norm79(:,1),temp79c) &
!              + intx5(e(:,:,OP_DR),ri_79,g(:,OP_1),norm79(:,2),temp79c) &
!              + 2.* &
!              (intx5(e(:,:,OP_DR),ri_79,g(:,OP_1),norm79(:,1),f(:,OP_DRZ)) &
!              -intx5(e(:,:,OP_DZ),ri_79,g(:,OP_1),norm79(:,2),f(:,OP_DRZ))) &
!              + intx5(e(:,:,OP_DZ),ri_79,g(:,OP_1),norm79(:,1),f(:,OP_LP)) &
!              - intx5(e(:,:,OP_DR),ri_79,g(:,OP_1),norm79(:,2),f(:,OP_LP))
!         if(itor.eq.1) then
!            temp79a = h(:,OP_1) - g(:,OP_1)
!            temp = temp &
!                 + 2.*intx5(e(:,:,OP_DR),ri2_79,g(:,OP_1),norm79(:,1),f(:,OP_DZ)) &
!                 + 2.*intx5(e(:,:,OP_1),ri2_79,g(:,OP_1),norm79(:,2),f(:,OP_LP)) &
!                 + 4.* &
!                 (intx5(e(:,:,OP_DR),ri2_79,temp79a,norm79(:,1),f(:,OP_DZ)) &
!                 +intx5(e(:,:,OP_DZ),ri2_79,temp79a,norm79(:,2),f(:,OP_DZ))) &
!                 - 4.* &
!                 (intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRZ),h(:,OP_1)) &
!                 +intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZZ),h(:,OP_1)))
!         endif

! #if defined(USE3D) || defined(USECOMPLEX)
!         temp = temp &
!              + intx5(e(:,:,OP_1),ri3_79,norm79(:,2),f(:,OP_DRPP),g(:,OP_1)) &
!              - intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DZPP),g(:,OP_1))
! #endif        

!      else
!         temp = 2.* &
!              (intx4(e(:,:,OP_DZZ),ri_79,f(:,OP_DRZ),g(:,OP_1)) &
!              -intx4(e(:,:,OP_DRR),ri_79,f(:,OP_DRZ),g(:,OP_1)) &
!              -intx4(e(:,:,OP_DRZ),ri_79,temp79c,g(:,OP_1)))
!         if(itor.eq.1) then
!            temp = temp - 2.* &
!                 (intx4(e(:,:,OP_DRR),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
!                 -intx4(e(:,:,OP_DR ),ri3_79,f(:,OP_DZ),g(:,OP_1))) &
!                 + 4.*intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRZ),g(:,OP_1)) &
!                 + 2.*intx4(e(:,:,OP_DZ),ri2_79,temp79c,g(:,OP_1))
           
!            temp79d = g(:,OP_1)-h(:,OP_1)
!            temp = temp + 4.*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DZ),temp79d)
!         endif
        
! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         temp = temp &
!              - intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZP),g(:,OP_1)) &
!              + intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DRP),g(:,OP_1)) &
!              - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_DP)) &
!              + intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_DP))
! #else
!         temp = temp &
!              + intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZPP),g(:,OP_1)) &
!              - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRPP),g(:,OP_1))
! #endif
! #endif
!      end if

!   v3umu = temp
! end function v3umu

function v3umu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3umu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  type(muarray) :: tempc

     temp79c = 1.
     tempc = mu(temp79c,OP_DZZ) + mu(-temp79c,OP_DRR)
     if(itor.eq.1) tempc = tempc + mu(-ri_79,OP_DR)

     if(surface_int) then
        temp = prod(mu(ri_79*g(:,OP_1)*norm79(:,1),OP_DZ),tempc) &
             + prod(mu(ri_79*g(:,OP_1)*norm79(:,2),OP_DR),tempc) &
             +(prod( 2.*ri_79*g(:,OP_1)*norm79(:,1),OP_DR,OP_DRZ) &
             + prod(-2.*ri_79*g(:,OP_1)*norm79(:,2),OP_DZ,OP_DRZ)) &
             + prod( ri_79*g(:,OP_1)*norm79(:,1),OP_DZ,OP_LP) &
             + prod(-ri_79*g(:,OP_1)*norm79(:,2),OP_DR,OP_LP)
        if(itor.eq.1) then
           temp79a = h(:,OP_1) - g(:,OP_1)
           temp = temp &
                + prod(2.*ri2_79*g(:,OP_1)*norm79(:,1),OP_DR,OP_DZ) &
                + prod(2.*ri2_79*g(:,OP_1)*norm79(:,2),OP_1,OP_LP) &
                +(prod(4.*ri2_79*temp79a*norm79(:,1),OP_DR,OP_DZ) &
                +prod(4.*ri2_79*temp79a*norm79(:,2),OP_DZ,OP_DZ)) &
                +(prod(-4.*ri2_79*norm79(:,1)*h(:,OP_1),OP_1,OP_DRZ) &
                +prod(-4.*ri2_79*norm79(:,2)*h(:,OP_1),OP_1,OP_DZZ))
         endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + prod( ri3_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DRPP) &
             + prod(-ri3_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DZPP)
#endif        

     else
        temp = &
             (prod( 2.*ri_79*g(:,OP_1),OP_DZZ,OP_DRZ) &
             +prod(-2.*ri_79*g(:,OP_1),OP_DRR,OP_DRZ) &
             +prod(mu(-2.*ri_79*g(:,OP_1),OP_DRZ),tempc))
         if(itor.eq.1) then
            temp = temp + &
                (prod(-2.*ri2_79*g(:,OP_1),OP_DRR,OP_DZ) &
                +prod( 2.*ri3_79*g(:,OP_1),OP_DR,OP_DZ)) &
                + prod(4.*ri2_79*g(:,OP_1),OP_DR,OP_DRZ) &
                + prod(mu(2.*ri2_79*g(:,OP_1),OP_DZ),tempc)
            
           temp79d = g(:,OP_1)-h(:,OP_1)
           temp = temp + prod(4.*ri2_79*temp79d,OP_GS,OP_DZ)
        endif
        
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
        temp = temp &
             + prod(-ri3_79*g(:,OP_1),OP_DRP,OP_DZP) &
             + prod( ri3_79*g(:,OP_1),OP_DZP,OP_DRP) &
             + prod(-ri3_79*g(:,OP_DP),OP_DR,OP_DZP) &
             + prod( ri3_79*g(:,OP_DP),OP_DZ,OP_DRP)
#else
        temp = temp &
             + prod( ri3_79*g(:,OP_1),OP_DR,OP_DZPP) &
             + prod(-ri3_79*g(:,OP_1),OP_DZ,OP_DRPP)
#endif
#endif
     end if

  v3umu = temp
end function v3umu


! V3vmu
! =====
!function v3vmu(e,f,g,h)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3vmu
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  !vectype, dimension(dofs_per_element) :: temp

  !temp = 0.
!#if defined(USE3D) || defined(USECOMPLEX)
     !if(surface_int) then
        !temp79a = h(:,OP_1) - g(:,OP_1)
        !temp = -2.* &
             !(intx5(e(:,:,OP_DR),ri2_79,norm79(:,1),f(:,OP_DP),temp79a) &
             !+intx5(e(:,:,OP_DZ),ri2_79,norm79(:,2),f(:,OP_DP),temp79a)) &
             !+ 2.* &
             !(intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRP),temp79a) &
             !+intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZP),temp79a)) &
             !+ intx5(e(:,:,OP_1),ri2_79,g(:,OP_1),norm79(:,1),f(:,OP_DRP)) &
             !+ intx5(e(:,:,OP_1),ri2_79,g(:,OP_1),norm79(:,2),f(:,OP_DZP))

        !if(itor.eq.1) then
           !temp = temp &
                !- 2.*intx5(e(:,:,OP_1),ri3_79,g(:,OP_1),norm79(:,1),f(:,OP_DP))
        !endif
     !else
        !temp = &
             !- intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
             !- intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1)) &
             !+ 2.*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),h(:,OP_1)) &
             !- 2.*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),g(:,OP_1))
        !if(itor.eq.1) then
           !temp = temp &
                !+ 2.*intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DP),g(:,OP_1))
        !endif
     !end if

!#endif
  !v3vmu = temp
!end function v3vmu

function v3vmu(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vmu
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp79a = h(:,OP_1) - g(:,OP_1)
        temp = &
             (prod(-2.*ri2_79*norm79(:,1)*temp79a,OP_DR,OP_DP) &
             +prod(-2.*ri2_79*norm79(:,2)*temp79a,OP_DZ,OP_DP)) &
             +(prod(2.*ri2_79*norm79(:,1)*temp79a,OP_1,OP_DRP) &
             + prod(2.*ri2_79*norm79(:,2)*temp79a,OP_1,OP_DZP)) &
             + prod(ri2_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DRP) &
             + prod(ri2_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DZP)

        if(itor.eq.1) then
           temp = temp &
                + prod(-2.*ri3_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DP)
         endif
     else
        temp = &
               prod(-ri2_79*g(:,OP_1),OP_DZ,OP_DZP) &
             + prod(-ri2_79*g(:,OP_1),OP_DR,OP_DRP) &
             + prod( 2.*ri2_79*h(:,OP_1),OP_GS,OP_DP) &
             + prod(-2.*ri2_79*g(:,OP_1),OP_GS,OP_DP)
         if(itor.eq.1) then
           temp = temp &
                + prod(2.*ri3_79*g(:,OP_1),OP_DR,OP_DP)
        endif
     end if
#endif
  v3vmu = temp
end function v3vmu


! V3un
! ====
!function v3un(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3un
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !if(inonormalflow.eq.1) then
           !temp = 0.
        !else
           !temp = intx5(e(:,:,OP_1),ri_79,g(:,OP_1),norm79(:,2),f(:,OP_DR)) &
                !- intx5(e(:,:,OP_1),ri_79,g(:,OP_1),norm79(:,1),f(:,OP_DZ))
        !end if
     !else
        !temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
             !- intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1)) 
     !end if

  !v3un = temp
!end function v3un

function v3un(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3un
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           temp = prod( ri_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DR) &
                + prod(-ri_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DZ)
         end if
     else
        temp = prod( ri_79*g(:,OP_1),OP_DR,OP_DZ) &
             + prod(-ri_79*g(:,OP_1),OP_DZ,OP_DR)
     end if

  v3un = temp
end function v3un


! V3p
! ===
function v3p(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3p
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        temp = 0.
!!$        temp = &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR)) &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ))
     else
        temp = intx3(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ)) &
             + intx3(e(:,:,OP_DR),ri2_79,f(:,OP_DR))
     end if

  v3p = temp
end function v3p

function v3p1

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3p1
  type(prodarray) :: temp

     if(surface_int) then
        temp%len = 0
!!$        temp = &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR)) &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ))
     else
        temp = prod(ri2_79,OP_DZ,OP_DZ) &
             + prod(ri2_79,OP_DR,OP_DR)
      end if

  v3p1 = temp
end function v3p1


! V3up
! ====
!function v3up(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3up
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !temp = intx5(e(:,:,OP_DR),ri_79,norm79(:,1),f(:,OP_DZ),g(:,OP_DR)) &
             !+ intx5(e(:,:,OP_DZ),ri_79,norm79(:,2),f(:,OP_DZ),g(:,OP_DR)) &
             !- intx5(e(:,:,OP_DR),ri_79,norm79(:,1),f(:,OP_DR),g(:,OP_DZ)) &
             !- intx5(e(:,:,OP_DZ),ri_79,norm79(:,2),f(:,OP_DR),g(:,OP_DZ)) &
             !+ intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             !- intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             !+ intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             !- intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             !+ intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             !- intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             !+ intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ)) &
             !- intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ))
        !if(itor.eq.1) then
           !temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           !temp = temp + 2.*gam* &
                !(intx5(e(:,:,OP_DR),norm79(:,1),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
                !+intx5(e(:,:,OP_DZ),norm79(:,2),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
                !-intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_1 )) &
                !-intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DR)) &
                !-intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_1 )) &
                !-intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZ))) &
                !-intx4(e(:,:,OP_1),ri2_79,norm79(:,1),temp79b)
        !endif
     !else
        !temp = intx4(e(:,:,OP_GS),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
             !- intx4(e(:,:,OP_GS),ri_79,f(:,OP_DZ),g(:,OP_DR))
        !if(itor.eq.1) then
           !temp = temp - 2.*gam* &
                !intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DZ),g(:,OP_1))
        !endif
     !end if

  !v3up = temp
!end function v3up

function v3up(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3up
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp
  type(muarray) :: tempb

     if(surface_int) then
        temp = prod( ri_79*norm79(:,1)*g(:,OP_DR),OP_DR,OP_DZ) &
             + prod( ri_79*norm79(:,2)*g(:,OP_DR),OP_DZ,OP_DZ) &
             + prod(-ri_79*norm79(:,1)*g(:,OP_DZ),OP_DR,OP_DR) &
             + prod(-ri_79*norm79(:,2)*g(:,OP_DZ),OP_DZ,OP_DR) &
             + prod( ri_79*norm79(:,1)*g(:,OP_DZ ),OP_1,OP_DRR) &
             + prod(-ri_79*norm79(:,1)*g(:,OP_DR ),OP_1,OP_DRZ) &
             + prod( ri_79*norm79(:,1)*g(:,OP_DRZ),OP_1,OP_DR) &
             + prod(-ri_79*norm79(:,1)*g(:,OP_DRR),OP_1,OP_DZ) &
             + prod( ri_79*norm79(:,2)*g(:,OP_DZ ),OP_1,OP_DRZ) &
             + prod(-ri_79*norm79(:,2)*g(:,OP_DR ),OP_1,OP_DZZ) &
             + prod( ri_79*norm79(:,2)*g(:,OP_DZZ),OP_1,OP_DR) &
             + prod(-ri_79*norm79(:,2)*g(:,OP_DRZ),OP_1,OP_DZ)
        if(itor.eq.1) then
           tempb = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
           temp = temp + &
                (prod( 2.*gam*norm79(:,1)*ri2_79*g(:,OP_1),OP_DR,OP_DZ) &
                +prod( 2.*gam*norm79(:,2)*ri2_79*g(:,OP_1),OP_DZ,OP_DZ) &
                +prod(-2.*gam*ri2_79*norm79(:,1)*g(:,OP_1 ),OP_1,OP_DRZ) &
                +prod(-2.*gam*ri2_79*norm79(:,1)*g(:,OP_DR),OP_1,OP_DZ) &
                +prod(-2.*gam*ri2_79*norm79(:,2)*g(:,OP_1 ),OP_1,OP_DZZ) &
                +prod(-2.*gam*ri2_79*norm79(:,2)*g(:,OP_DZ),OP_1,OP_DZ)) &
                +prod(mu(-ri2_79*norm79(:,1),OP_1),tempb)
        endif
     else
        temp = prod( ri_79*g(:,OP_DZ),OP_GS,OP_DR) &
             + prod(-ri_79*g(:,OP_DR),OP_GS,OP_DZ)
        if(itor.eq.1) then
           temp = temp + &
                prod(-2.*gam*ri2_79*g(:,OP_1),OP_GS,OP_DZ)
        endif
     end if

  v3up = temp
end function v3up


! V3vp
! ====
!function v3vp(e,f,g)

  !use basic
  !use arrays
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3vp
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g

  !vectype, dimension(dofs_per_element) :: temp
  !temp = 0.

!#if defined(USE3D) || defined(USECOMPLEX)
     !if(surface_int) then
        !temp = &
             !- intx5(e(:,:,OP_DR),ri2_79,norm79(:,1),f(:,OP_1 ),g(:,OP_DP)) &
             !- intx5(e(:,:,OP_DZ),ri2_79,norm79(:,2),f(:,OP_1 ),g(:,OP_DP)) &
             !- gam*intx5(e(:,:,OP_DR),ri2_79,norm79(:,1),f(:,OP_DP),g(:,OP_1 )) &
             !- gam*intx5(e(:,:,OP_DZ),ri2_79,norm79(:,2),f(:,OP_DP),g(:,OP_1 )) &
             !+ intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR),g(:,OP_DP )) &
             !+ intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_1 ),g(:,OP_DRP)) &
             !+ intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ),g(:,OP_DP )) &
             !+ intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_1 ),g(:,OP_DZP)) &
             !+ gam * &
             !(intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRP),g(:,OP_1 )) &
             !+intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DP ),g(:,OP_DR)) &
             !+intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZP),g(:,OP_1 )) &
             !+intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DP ),g(:,OP_DZ)))
     !else
        !temp = gam*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),g(:,OP_1)) &
             !+ intx4(e(:,:,OP_GS),ri2_79,f(:,OP_1),g(:,OP_DP))
     !end if
!#endif

  !v3vp = temp
!end function v3vp

function v3vp(g)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vp
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp
  temp%len = 0

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp = &
               prod(-ri2_79*norm79(:,1)*g(:,OP_DP),OP_DR,OP_1) &
             + prod(-ri2_79*norm79(:,2)*g(:,OP_DP),OP_DZ,OP_1) &
             + prod(-gam*ri2_79*norm79(:,1)*g(:,OP_1 ),OP_DR,OP_DP) &
             + prod(-gam*ri2_79*norm79(:,2)*g(:,OP_1 ),OP_DZ,OP_DP) &
             + prod(ri2_79*norm79(:,1)*g(:,OP_DP ),OP_1,OP_DR) &
             + prod(ri2_79*norm79(:,1)*g(:,OP_DRP),OP_1,OP_1) &
             + prod(ri2_79*norm79(:,2)*g(:,OP_DP ),OP_1,OP_DZ) &
             + prod(ri2_79*norm79(:,2)*g(:,OP_DZP),OP_1,OP_1) &
             +(prod(gam*ri2_79*norm79(:,1)*g(:,OP_1 ),OP_1,OP_DRP) &
             + prod(gam*ri2_79*norm79(:,1)*g(:,OP_DR),OP_1,OP_DP) &
             + prod(gam*ri2_79*norm79(:,2)*g(:,OP_1 ),OP_1,OP_DZP) &
             + prod(gam*ri2_79*norm79(:,2)*g(:,OP_DZ),OP_1,OP_DP))
      else
        temp = prod(gam*ri2_79*g(:,OP_1),OP_GS,OP_DP) &
             + prod(ri2_79*g(:,OP_DP),OP_GS,OP_1)
     end if
#endif

  v3vp = temp
end function v3vp



! V3chip
! ======
!function v3chip(e,f,g)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3chip
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

     !if(surface_int) then
        !temp = &
             !- intx5(e(:,:,OP_DR),ri4_79,norm79(:,1),f(:,OP_DZ),g(:,OP_DZ)) &
             !- intx5(e(:,:,OP_DZ),ri4_79,norm79(:,2),f(:,OP_DZ),g(:,OP_DZ)) &
             !- intx5(e(:,:,OP_DR),ri4_79,norm79(:,1),f(:,OP_DR),g(:,OP_DR)) &
             !- intx5(e(:,:,OP_DZ),ri4_79,norm79(:,2),f(:,OP_DR),g(:,OP_DR)) &
             !- gam*intx5(e(:,:,OP_DR),ri4_79,norm79(:,1),f(:,OP_GS),g(:,OP_1))  &
             !- gam*intx5(e(:,:,OP_DZ),ri4_79,norm79(:,2),f(:,OP_GS),g(:,OP_1))  &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
             !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
        !if(itor.eq.1) then
           !temp = temp - 2.* &
                !(intx5(e(:,:,OP_1),ri5_79,norm79(:,1),f(:,OP_DZ),g(:,OP_DZ)) &
                !+intx5(e(:,:,OP_1),ri5_79,norm79(:,1),f(:,OP_DR),g(:,OP_DR)))
        !endif
        !!! missing some terms
     !else
        !temp = intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DZ),g(:,OP_DZ)) &
             !+ intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DR),g(:,OP_DR)) &
             !+ gam*intx4(e(:,:,OP_GS),ri4_79,f(:,OP_GS),g(:,OP_1))
     !end if

  !v3chip = temp
!end function v3chip

function v3chip(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chip
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then
        temp = &
               prod(-ri4_79*norm79(:,1)*g(:,OP_DZ),OP_DR,OP_DZ) &
             + prod(-ri4_79*norm79(:,2)*g(:,OP_DZ),OP_DZ,OP_DZ) &
             + prod(-ri4_79*norm79(:,1)*g(:,OP_DR),OP_DR,OP_DR) &
             + prod(-ri4_79*norm79(:,2)*g(:,OP_DR),OP_DZ,OP_DR) &
             + prod(-gam*ri4_79*norm79(:,1)*g(:,OP_1),OP_DR,OP_GS)  &
             + prod(-gam*ri4_79*norm79(:,2)*g(:,OP_1),OP_DZ,OP_GS)  &
             + prod(ri4_79*norm79(:,1)*g(:,OP_DZ ),OP_1,OP_DRZ) &
             + prod(ri4_79*norm79(:,1)*g(:,OP_DR ),OP_1,OP_DRR) &
             + prod(ri4_79*norm79(:,1)*g(:,OP_DRZ),OP_1,OP_DZ) &
             + prod(ri4_79*norm79(:,1)*g(:,OP_DRR),OP_1,OP_DR) &
             + prod(ri4_79*norm79(:,2)*g(:,OP_DZ ),OP_1,OP_DZZ) &
             + prod(ri4_79*norm79(:,2)*g(:,OP_DR ),OP_1,OP_DRZ) &
             + prod(ri4_79*norm79(:,2)*g(:,OP_DZZ),OP_1,OP_DZ) &
             + prod(ri4_79*norm79(:,2)*g(:,OP_DRZ),OP_1,OP_DR)
        if(itor.eq.1) then
           temp = temp + &
                (prod(-2.*ri5_79*norm79(:,1)*g(:,OP_DZ),OP_1,OP_DZ) &
                +prod(-2.*ri5_79*norm79(:,1)*g(:,OP_DR),OP_1,OP_DR))
        endif
        !! missing some terms
     else
        temp = prod(ri4_79*g(:,OP_DZ),OP_GS,OP_DZ) &
             + prod(ri4_79*g(:,OP_DR),OP_GS,OP_DR) &
             + prod(gam*ri4_79*g(:,OP_1),OP_GS,OP_GS)
     end if

  v3chip = temp
end function v3chip



! V3psipsi
! ========
function v3psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3psipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        temp = &
             - intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DR),g(:,OP_GS)) &
             - intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DZ),g(:,OP_GS))
     else
        temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_GS)) &
             + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_GS))
     end if

  v3psipsi = temp
end function v3psipsi

function v3psipsi1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3psipsi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

     if(surface_int) then
        temp = &
               prod(-ri4_79*norm79(:,1)*g(:,OP_GS),OP_1,OP_DR) &
             + prod(-ri4_79*norm79(:,2)*g(:,OP_GS),OP_1,OP_DZ)
      else
        temp = prod(ri4_79*g(:,OP_GS),OP_DZ,OP_DZ) &
             + prod(ri4_79*g(:,OP_GS),OP_DR,OP_DR)
      end if

  v3psipsi1 = temp
end function v3psipsi1

function v3psipsi2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3psipsi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  type(prodarray) :: temp

     if(surface_int) then
        temp = &
               prod(-ri4_79*norm79(:,1)*f(:,OP_DR),OP_1,OP_GS) &
             + prod(-ri4_79*norm79(:,2)*f(:,OP_DZ),OP_1,OP_GS)
      else
        temp = prod(ri4_79*f(:,OP_DZ),OP_DZ,OP_GS) &
             + prod(ri4_79*f(:,OP_DR),OP_DR,OP_GS)
      end if

  v3psipsi2 = temp
end function v3psipsi2


! V3psib
! ======
function v3psib(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3psib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then 
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),ri5_79,norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
                - intx5(e(:,:,OP_1),ri5_79,norm79(:,1),f(:,OP_DZP),g(:,OP_1))
        end if
     else
        temp = - &
             (intx4(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP),g(:,OP_1)) &
             -intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DZP),g(:,OP_1)))
     end if
#else
  temp = 0.
#endif
  v3psib = temp
end function v3psib

function v3psib1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3psib1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then 
           temp%len = 0
        else
           temp = prod( ri5_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DRP) &
                + prod(-ri5_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DZP)
         end if
     else
        temp = &
             (prod(-ri5_79*g(:,OP_1),OP_DZ,OP_DRP) &
             +prod( ri5_79*g(:,OP_1),OP_DR,OP_DZP))
      end if
#else
  temp%len = 0
#endif
  v3psib1 = temp
end function v3psib1

function v3psib2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3psib2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then 
           temp%len = 0
        else
           temp = prod( ri5_79*norm79(:,2)*f(:,OP_DRP),OP_1,OP_1) &
                + prod(-ri5_79*norm79(:,1)*f(:,OP_DZP),OP_1,OP_1)
        end if
     else
        temp = &
             (prod(-ri5_79*f(:,OP_DRP),OP_DZ,OP_1) &
             +prod( ri5_79*f(:,OP_DZP),OP_DR,OP_1))
      end if
#else
  temp%len = 0
#endif
  v3psib2 = temp
end function v3psib2


! V3psif
! ======
function v3psif(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3psif
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_tor.eq.1 .and. imulti_region.eq.0) then
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),ri3_79,f(:,OP_GS),norm79(:,1),g(:,OP_DZ)) &
                - intx5(e(:,:,OP_1),ri3_79,f(:,OP_GS),norm79(:,2),g(:,OP_DR))
        end if
     else
        temp = intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DR),f(:,OP_GS)) &
             - intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DZ),f(:,OP_GS))
     end if
#else
  temp = 0.
#endif

  v3psif = temp
end function v3psif

function v3psif1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3psif1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_tor.eq.1 .and. imulti_region.eq.0) then
           temp%len = 0
        else
           temp = prod( ri3_79*norm79(:,1)*g(:,OP_DZ),OP_1,OP_GS) &
                + prod(-ri3_79*norm79(:,2)*g(:,OP_DR),OP_1,OP_GS)
        end if
     else
        temp = prod( ri3_79*g(:,OP_DR),OP_DZ,OP_GS) &
             + prod(-ri3_79*g(:,OP_DZ),OP_DR,OP_GS)
     end if
#else
  temp%len = 0
#endif

  v3psif1 = temp
end function v3psif1

function v3psif2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3psif2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_tor.eq.1 .and. imulti_region.eq.0) then
           temp%len = 0
        else
           temp = prod( ri3_79*f(:,OP_GS)*norm79(:,1),OP_1,OP_DZ) &
                + prod(-ri3_79*f(:,OP_GS)*norm79(:,2),OP_1,OP_DR)
        end if
     else
        temp = prod( ri3_79*f(:,OP_GS),OP_DZ,OP_DR) &
             + prod(-ri3_79*f(:,OP_GS),OP_DR,OP_DZ)
     end if
#else
  temp%len = 0
#endif

  v3psif2 = temp
end function v3psif2


! V3bb
! ====
function v3bb(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3bb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
           temp = 0.
        else
           temp = &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DZ),g(:,OP_1))
        end if
     else
        temp = intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),g(:,OP_1)) &
             + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),g(:,OP_1))
     end if

  v3bb = temp
end function v3bb

function v3bb1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3bb1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
           temp%len = 0
        else
           temp = &
                  prod(-ri4_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DR) &
                + prod(-ri4_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DZ)
         end if
     else
        temp = prod(ri4_79*g(:,OP_1),OP_DR,OP_DR) &
             + prod(ri4_79*g(:,OP_1),OP_DZ,OP_DZ)
      end if

  v3bb1 = temp
end function v3bb1

function v3bb2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3bb2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  type(prodarray) :: temp

     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
           temp%len = 0
        else
           temp = &
                  prod(-ri4_79*norm79(:,1)*f(:,OP_DR),OP_1,OP_1) &
                + prod(-ri4_79*norm79(:,2)*f(:,OP_DZ),OP_1,OP_1)
         end if
     else
        temp = prod(ri4_79*f(:,OP_DR),OP_DR,OP_1) &
             + prod(ri4_79*f(:,OP_DZ),OP_DZ,OP_1)
      end if

  v3bb2 = temp
end function v3bb2


! V3bf
! ====
function v3bf(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3bf
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then 
           temp = 0.
        else
           temp = &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,1),g(:,OP_DRP),f(:,OP_1)) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,2),g(:,OP_DZP),f(:,OP_1))
        end if
     else
        temp = intx4(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_1)) &
             + intx4(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_1))
     end if
#else
  temp = 0.
#endif

  v3bf = temp
end function v3bf

function v3bf1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3bf1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then 
           temp%len = 0
        else
           temp = &
                  prod(-ri4_79*norm79(:,1)*g(:,OP_DRP),OP_1,OP_1) &
                + prod(-ri4_79*norm79(:,2)*g(:,OP_DZP),OP_1,OP_1)
        end if
     else
        temp = prod(ri4_79*g(:,OP_DRP),OP_DR,OP_1) &
             + prod(ri4_79*g(:,OP_DZP),OP_DZ,OP_1)
     end if
#else
  temp%len = 0
#endif

  v3bf1 = temp
end function v3bf1

function v3bf2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3bf2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then 
           temp%len = 0
        else
           temp = &
                  prod(-ri4_79*norm79(:,1)*f(:,OP_1),OP_1,OP_DRP) &
                + prod(-ri4_79*norm79(:,2)*f(:,OP_1),OP_1,OP_DZP)
         end if
     else
        temp = prod(ri4_79*f(:,OP_1),OP_DR,OP_DRP) &
             + prod(ri4_79*f(:,OP_1),OP_DZ,OP_DZP)
     end if
#else
  temp%len = 0
#endif

  v3bf2 = temp
end function v3bf2


! V3psisb1
! ========
vectype function v3psisb1(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ)) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR)) &
          + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS)) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS))
  end if

  v3psisb1 = temp

end function v3psisb1


! V3bsb2
! ======
vectype function v3bsb2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 )) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1 )) &
          + int4(ri2_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ)) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_1 ),g(:,OP_DR))
  end if

  v3bsb2 = temp
  return
end function v3bsb2


! V3upsipsi
! =========
!function v3upsipsi(e,f,g,h)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3upsipsi
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  !vectype, dimension(dofs_per_element) :: temp
  !vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempd
  !integer :: j

  !if(surface_int) then
        !do j=1, dofs_per_element
           !tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
        !end do
        !temp = intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             !- intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             !+ intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             !- intx5(tempa,ri3_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             !+ intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             !- intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             !+ intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
             !- intx5(tempa,ri3_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ))
        !if(itor.eq.1) then
           !temp = temp &
                !+ intx5(tempa,ri4_79,norm79(:,1),f(:,OP_DZ),g(:,OP_DR)) &
                !- intx5(tempa,ri4_79,norm79(:,1),f(:,OP_DR),g(:,OP_DZ))
        !endif
  !else
  
  !! [f,g],r
  !temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ) &
       !+    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
  
  !! [f,g],z
  !temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ) &
       !+    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)
  
  !temp = intx4(e(:,:,OP_DZ ),ri3_79,temp79c,h(:,OP_GS )) &
       !+ intx4(e(:,:,OP_DR ),ri3_79,temp79b,h(:,OP_GS )) &
       !- intx4(e(:,:,OP_DZZ),ri3_79,temp79c,h(:,OP_DZ )) &
       !- intx4(e(:,:,OP_DR ),ri3_79,temp79c,h(:,OP_DRZ)) &
       !- intx4(e(:,:,OP_DZ ),ri3_79,temp79c,h(:,OP_DZZ)) &
       !- intx4(e(:,:,OP_DRZ),ri3_79,temp79c,h(:,OP_DR )) &
       !- intx4(e(:,:,OP_DRZ),ri3_79,temp79b,h(:,OP_DZ )) &
       !- intx4(e(:,:,OP_DR ),ri3_79,temp79b,h(:,OP_DRR)) &
       !- intx4(e(:,:,OP_DZ ),ri3_79,temp79b,h(:,OP_DRZ)) &
       !- intx4(e(:,:,OP_DRR),ri3_79,temp79b,h(:,OP_DR ))
    
  !if(itor.eq.1) then
     !! [f,g]
     !temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
  

        !do j=1, dofs_per_element
           !tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
        !end do
        
        !temp = temp &
             !+2.*intx3(tempd,ri4_79,temp79b) &
             !-   intx4(e(:,:,OP_DRZ),ri4_79,temp79a,h(:,OP_DZ )) &
             !-   intx4(e(:,:,OP_DZ ),ri4_79,temp79a,h(:,OP_DRZ)) &
             !-   intx4(e(:,:,OP_DRR),ri4_79,temp79a,h(:,OP_DR )) &
             !-   intx4(e(:,:,OP_DR ),ri4_79,temp79a,h(:,OP_DRR)) &
             !+   intx4(e(:,:,OP_DR),ri4_79,temp79a,h(:,OP_GS)) &
             !+2.*intx3(tempd,ri5_79,temp79a)
  !endif

  !end if

  !v3upsipsi = temp
!end function v3upsipsi

function v3upsipsi(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3upsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd

  if(surface_int) then
        tempa = mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
        temp = prod(tempa,mu( ri3_79*norm79(:,1)*g(:,OP_DR ),OP_DRZ)) &
             + prod(tempa,mu(-ri3_79*norm79(:,1)*g(:,OP_DZ ),OP_DRR)) &
             + prod(tempa,mu( ri3_79*norm79(:,1)*g(:,OP_DRR),OP_DZ)) &
             + prod(tempa,mu(-ri3_79*norm79(:,1)*g(:,OP_DRZ),OP_DR)) &
             + prod(tempa,mu( ri3_79*norm79(:,2)*g(:,OP_DR ),OP_DZZ)) &
             + prod(tempa,mu(-ri3_79*norm79(:,2)*g(:,OP_DZ ),OP_DRZ)) &
             + prod(tempa,mu( ri3_79*norm79(:,2)*g(:,OP_DRZ),OP_DZ)) &
             + prod(tempa,mu(-ri3_79*norm79(:,2)*g(:,OP_DZZ),OP_DR))
         if(itor.eq.1) then
           temp = temp &
                + prod(tempa,mu( ri4_79*norm79(:,1)*g(:,OP_DR),OP_DZ)) &
                + prod(tempa,mu(-ri4_79*norm79(:,1)*g(:,OP_DZ),OP_DR))
         endif
  else
  
  ! [f,g],r
  tempb = mu(g(:,OP_DR ),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &
       +  mu(g(:,OP_DRR),OP_DZ)  + mu(-g(:,OP_DRZ),OP_DR)
 
  ! [f,g],z
  tempc = mu(g(:,OP_DR ),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &
       +    mu(g(:,OP_DRZ),OP_DZ)  + mu(-g(:,OP_DZZ),OP_DR)
 
  temp = prod(mu( ri3_79*h(:,OP_GS),OP_DZ),tempc) &
       + prod(mu( ri3_79*h(:,OP_GS),OP_DR),tempb) &
       + prod(mu(-ri3_79*h(:,OP_DZ),OP_DZZ),tempc) &
       + prod(mu(-ri3_79*h(:,OP_DRZ),OP_DR),tempc) &
       + prod(mu(-ri3_79*h(:,OP_DZZ),OP_DZ),tempc) &
       + prod(mu(-ri3_79*h(:,OP_DR),OP_DRZ),tempc) &
       + prod(mu(-ri3_79*h(:,OP_DZ),OP_DRZ),tempb) &
       + prod(mu(-ri3_79*h(:,OP_DRR),OP_DR),tempb) &
       + prod(mu(-ri3_79*h(:,OP_DRZ),OP_DZ),tempb) &
       + prod(mu(-ri3_79*h(:,OP_DR),OP_DRR),tempb)
     
  if(itor.eq.1) then
     ! [f,g]
     tempa = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)
  
           tempd = mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
        
        temp = temp &
             +   prod(tempd*(2.*ri4_79),tempb) &
             +   prod(mu(-ri4_79*h(:,OP_DZ),OP_DRZ),tempa) &
             +   prod(mu(-ri4_79*h(:,OP_DRZ),OP_DZ),tempa) &
             +   prod(mu(-ri4_79*h(:,OP_DR),OP_DRR),tempa) &
             +   prod(mu(-ri4_79*h(:,OP_DRR),OP_DR),tempa) &
             +   prod( mu(ri4_79*h(:,OP_GS),OP_DR),tempa) &
             +   prod(tempd*(2.*ri5_79),tempa)
  endif

  end if

  v3upsipsi = temp
end function v3upsipsi


! V3upsib
! =======
! function v3upsib(e,f,g,h)
!   use basic
!   use arrays
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3upsib
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

! #if defined(USE3D) || defined(USECOMPLEX) 
!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempd, tempf
!   integer :: j

!      if(surface_int) then
!         do j=1, dofs_per_element
!            tempa(j,:) = h(:,OP_1)* &
!                 (norm79(:,1)*e(j,:,OP_DZ) - norm79(:,2)*e(j,:,OP_DR))
!         end do
!         temp = intx4(tempa,ri4_79,f(:,OP_DRP),g(:,OP_DZ )) &
!              - intx4(tempa,ri4_79,f(:,OP_DZP),g(:,OP_DR )) &
!              + intx4(tempa,ri4_79,f(:,OP_DR ),g(:,OP_DZP)) &
!              - intx4(tempa,ri4_79,f(:,OP_DZ ),g(:,OP_DRP))
!      else
!         do j=1, dofs_per_element
!            tempa(j,:) = e(j,:,OP_DZ)*g(:,OP_DRP) - e(j,:,OP_DR)*g(:,OP_DZP)
!            tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DR)  - e(j,:,OP_DR)*h(:,OP_DZ)
!            tempf(j,:) = h(:,OP_DP)* &
!                 (e(j,:,OP_DZ)*f(:,OP_DZ )+e(j,:,OP_DR)*f(:,OP_DR )) &
!                 +    h(:,OP_1 )* &
!                 (e(j,:,OP_DZ)*f(:,OP_DZP)+e(j,:,OP_DR)*f(:,OP_DRP))
!            if(itor.eq.1) then 
!               tempd(j,:) = tempd(j,:) - 4.*ri_79*e(j,:,OP_DZ)*h(:,OP_1)
!            endif
!         end do

!         temp79b = h(:,OP_DZ)*f(:,OP_DR) - h(:,OP_DR)*f(:,OP_DZ)
!         temp79c = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &
!              +    f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP)
! #ifdef USEST
!         temp79e =  h(:,OP_1)*f(:,OP_GS) &
!              + f(:,OP_DZ )*h(:,OP_DZ) + f(:,OP_DR )*h(:,OP_DR)
! #else
!         temp79e = h(:,OP_DP)*f(:,OP_GS) + h(:,OP_1)*f(:,OP_GSP) &
!              + f(:,OP_DZP)*h(:,OP_DZ ) + f(:,OP_DRP)*h(:,OP_DR ) &
!              + f(:,OP_DZ )*h(:,OP_DZP) + f(:,OP_DR )*h(:,OP_DRP)
! #endif
!         temp = intx3(tempa,ri4_79,temp79b) &
!              + intx3(tempd,ri4_79,temp79c) &
! #ifdef USEST
!              - intx4(e(:,:,OP_DZP),ri4_79,temp79e,g(:,OP_DZ)) &
!              - intx4(e(:,:,OP_DRP),ri4_79,temp79e,g(:,OP_DR)) &
!              - intx4(e(:,:,OP_DZ),ri4_79,temp79e,g(:,OP_DZP)) &
!              - intx4(e(:,:,OP_DR),ri4_79,temp79e,g(:,OP_DRP)) &
! #else
!              + intx4(e(:,:,OP_DZ),ri4_79,temp79e,g(:,OP_DZ)) &
!              + intx4(e(:,:,OP_DR),ri4_79,temp79e,g(:,OP_DR)) &
! #endif
!              + intx3(tempf,ri4_79,g(:,OP_GS))
!      end if

!   v3upsib = temp
! #else
!   v3upsib = 0.
! #endif
! end function v3upsib

function v3upsib(g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3upsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX) 
  type(prodarray) :: temp, tempf
  type(muarray) :: tempa, tempb, tempc, tempd, tempe

     if(surface_int) then
           tempa = &
                (mu(norm79(:,1),OP_DZ) + mu(-norm79(:,2),OP_DR))*h(:,OP_1)
         temp = prod(tempa,mu( ri4_79*g(:,OP_DZ ),OP_DRP)) &
             +  prod(tempa,mu(-ri4_79*g(:,OP_DR ),OP_DZP)) &
             +  prod(tempa,mu( ri4_79*g(:,OP_DZP),OP_DR)) &
             +  prod(tempa,mu(-ri4_79*g(:,OP_DRP),OP_DZ))
      else
        tempa = mu(g(:,OP_DRP),OP_DZ) + mu(-g(:,OP_DZP),OP_DR)
        tempd = mu(h(:,OP_DR),OP_DZ)  + mu(-h(:,OP_DZ),OP_DR)
        tempf = prod(h(:,OP_DP),OP_DZ,OP_DZ )+prod(h(:,OP_DP),OP_DR,OP_DR ) &
           +    prod(h(:,OP_1),OP_DZ,OP_DZP)+prod(h(:,OP_1),OP_DR,OP_DRP)
        if(itor.eq.1) then
           tempd = tempd + mu(-4.*ri_79*h(:,OP_1),OP_DZ)
        endif
            
        tempb = mu(h(:,OP_DZ),OP_DR) + mu(-h(:,OP_DR),OP_DZ)
        tempc = mu(g(:,OP_DR),OP_DZP)+ mu(-g(:,OP_DZ),OP_DRP) &
             +  mu(g(:,OP_DRP),OP_DZ)+ mu(-g(:,OP_DZP),OP_DR)

#ifdef USEST
        tempe = mu(h(:,OP_1),OP_GS) &
             +  mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
#else
        tempe = mu(h(:,OP_DP),OP_GS) + mu(h(:,OP_1),OP_GSP) &
             +  mu(h(:,OP_DZ ),OP_DZP)+mu(h(:,OP_DR ),OP_DRP) &
             +  mu(h(:,OP_DZP),OP_DZ) +mu(h(:,OP_DRP),OP_DR)
#endif
        temp = prod(tempa*ri4_79,tempb) &
             + prod(tempd*ri4_79,tempc) &
#ifdef USEST
             + prod(mu(-ri4_79*g(:,OP_DZ),OP_DZP),tempe) &
             + prod(mu(-ri4_79*g(:,OP_DR),OP_DRP),tempe) &
             + prod(mu(-ri4_79*g(:,OP_DZP),OP_DZ),tempe) &
             + prod(mu(-ri4_79*g(:,OP_DRP),OP_DR),tempe) &
#else
             + prod(mu(ri4_79*g(:,OP_DZ),OP_DZ),tempe) &
             + prod(mu(ri4_79*g(:,OP_DR),OP_DR),tempe) &
#endif
             + tempf*(ri4_79*g(:,OP_GS))
     end if
  v3upsib = temp
#else
  v3upsib%len = 0
#endif
end function v3upsib


! V3ubb
! =====
! function v3ubb(e,f,g,h)
!   use basic
!   use arrays
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3ubb
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb
!   integer :: j

!      if(surface_int) then
!         do j=1, dofs_per_element
!            tempa(j,:) = h(:,OP_1)* &
!                 (norm79(:,1)*e(j,:,OP_DR) + norm79(:,2)*e(j,:,OP_DZ))
!            if(itor.eq.1) then
!               tempa(j,:) = tempa(j,:) &
!                    - 2.*ri_79*norm79(:,1)*e(j,:,OP_1)*h(:,OP_1)
!            end if
!         end do
!         temp = intx4(tempa,ri3_79,f(:,OP_DZ),g(:,OP_DR)) &
!              - intx4(tempa,ri3_79,f(:,OP_DR),g(:,OP_DZ))
!      else
!         temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

!         temp = intx3(e(:,:,OP_GS),ri3_79,temp79a)

! !  scj removed 4/1/2011
!         if(itor.eq.1) then
!            temp = temp - &
!                 2.*intx3(e(:,:,OP_DR),ri4_79,temp79a)
!         endif

! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         do j=1, dofs_per_element 
!            tempb(j,:) = &
!                  (e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
!                 *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP)) &
!                 -(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
!                 *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP))
!         end do
!         temp = temp + intx2(tempb,ri5_79)
! #else
!         do j=1, dofs_per_element
!            tempb(j,:) = &
!               (e(j,:,OP_DZ)*f(:,OP_DR)-e(j,:,OP_DR)*f(:,OP_DZ))*g(:,OP_DPP) &
!          + 2.*(e(j,:,OP_DZ)*f(:,OP_DRP)-e(j,:,OP_DR)*f(:,OP_DZP))*g(:,OP_DP) &
!          +    (e(j,:,OP_DZ)*f(:,OP_DRPP)-e(j,:,OP_DR)*f(:,OP_DZPP))*g(:,OP_1)
!         end do
!         temp = temp - intx3(tempb,ri5_79,h(:,OP_1))
! #endif
! #endif
!      end if

!   v3ubb = temp
! end function v3ubb

function v3ubb(g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3ubb
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp, tempb
  type(muarray) :: tempa

     if(surface_int) then
           tempa = (mu(norm79(:,1),OP_DR) + mu(norm79(:,2),OP_DZ))*h(:,OP_1)
           if(itor.eq.1) then
              tempa = tempa &
                   + mu(-2.*ri_79*norm79(:,1)*h(:,OP_1),OP_1)
           end if
         temp = prod(tempa,mu( ri3_79*g(:,OP_DR),OP_DZ)) &
              + prod(tempa,mu(-ri3_79*g(:,OP_DZ),OP_DR))
      else
        tempa = (mu(g(:,OP_DZ),OP_DR) + mu(-g(:,OP_DR),OP_DZ))*h(:,OP_1)

        temp = prod(mu(ri3_79,OP_GS),tempa)

!  scj removed 4/1/2011
        if(itor.eq.1) then
           temp = temp + &
                  prod(mu(-2.*ri4_79,OP_DR),tempa)
        endif

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
           tempb = &
                 prod(mu(h(:,OP_DP),OP_DZ) + mu(h(:,OP_1),OP_DZP), &
                 mu(g(:,OP_1),OP_DRP) + mu(g(:,OP_DP),OP_DR)) &
                +prod(mu(-h(:,OP_DP),OP_DR) + mu(-h(:,OP_1),OP_DRP), &
                 mu(g(:,OP_1),OP_DZP) + mu(g(:,OP_DP),OP_DZ))
        temp = temp + tempb*ri5_79
#else
            tempb = &
              prod(g(:,OP_DPP),OP_DZ,OP_DR)+prod(-g(:,OP_DPP),OP_DR,OP_DZ) &
         +    prod(2.*g(:,OP_DP),OP_DZ,OP_DRP)+prod(-2.*g(:,OP_DP),OP_DR,OP_DZP) &
         +    prod(g(:,OP_1),OP_DZ,OP_DRPP)+prod(-g(:,OP_1),OP_DR,OP_DZPP)
         temp = temp + tempb*(-ri5_79*h(:,OP_1))
#endif
#endif
     end if

  v3ubb = temp
end function v3ubb
#ifdef USE3D
! V3upsif
! =====
! function v3upsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3upsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempd
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! ((nu, psi)/R^2)_R*R^2
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DRZ) + e(j,:,OP_DR)*g(:,OP_DRR) &
!               + e(j,:,OP_DRZ)*g(:,OP_DZ) + e(j,:,OP_DRR)*g(:,OP_DR) 
!            if(itor.eq.1) then
!               tempa(j,:) = tempa(j,:) - 2*ri_79*&
!                    (e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DZ)) 
!            end if 
!            ! ((nu, psi)/R^2)_Z*R^2
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DZZ) + e(j,:,OP_DR)*g(:,OP_DRZ) &
!               + e(j,:,OP_DZZ)*g(:,OP_DZ) + e(j,:,OP_DRZ)*g(:,OP_DR) 
!            ! [nu, f']_R*R
!            tempc(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DRR) - e(j,:,OP_DR)*h(:,OP_DRZ) &
!               + e(j,:,OP_DRZ)*h(:,OP_DR) - e(j,:,OP_DRR)*h(:,OP_DZ)
!            if(itor.eq.1) then
!               tempc(j,:) = tempc(j,:) - ri_79*&
!                    (e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)) 
!            end if 
!            ! [nu, f']_Z*R
!            tempd(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DRZ) - e(j,:,OP_DR)*h(:,OP_DZZ) &
!               + e(j,:,OP_DZZ)*h(:,OP_DR) - e(j,:,OP_DRZ)*h(:,OP_DZ)
!         end do

!         ! ([u, psi]*R^2)_R/R
!         temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
!         if(itor.eq.1) then
!            temp79a = temp79a + ri_79* &
!                   (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
!         end if 
!         ! ([u, psi]*R^2)_Z/R
!         temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            
!         ! (R^2*(u, f'))_R/R^2
!         temp79c = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
!         if(itor.eq.1) then
!            temp79c = temp79c + 2*ri_79* &           
!                   (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
!         endif
!         ! (R^2*(u, f'))_Z/R^2
!         temp79d = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            

!         temp = intx3(tempa,temp79c,ri2_79) &
!              + intx3(tempb,temp79d,ri2_79) &
!              - intx3(tempc,temp79a,ri2_79) &
!              - intx3(tempd,temp79b,ri2_79) &
!              - intx4(e(:,:,OP_DZ),temp79d,g(:,OP_GS),ri2_79) & 
!              - intx4(e(:,:,OP_DR),temp79c,g(:,OP_GS),ri2_79) 
!      end if

!   v3upsif = temp
! end function v3upsif

function v3upsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3upsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd
  type(muarray) :: tempa2, tempb2, tempc2, tempd2

  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! ((nu, psi)/R^2)_R*R^2
           tempa = &
                mu(g(:,OP_DRZ),OP_DZ) + mu(g(:,OP_DRR),OP_DR) &
              + mu(g(:,OP_DZ),OP_DRZ) + mu(g(:,OP_DR),OP_DRR) 
           if(itor.eq.1) then
              tempa = tempa + &
                   (mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DZ),OP_DR))*(-2*ri_79) 
           end if 
           ! ((nu, psi)/R^2)_Z*R^2
           tempb = &
                mu(g(:,OP_DZZ),OP_DZ) + mu(g(:,OP_DRZ),OP_DR) &
              + mu(g(:,OP_DZ),OP_DZZ) + mu(g(:,OP_DR),OP_DRZ) 
           ! [nu, f']_R*R
           tempc = &
                mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR) &
              + mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR)
           if(itor.eq.1) then
              tempc = tempc + &
                   (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79) 
           end if 
           ! [nu, f']_Z*R
           tempd = &
                mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR) &
              + mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ)

        ! ([u, psi]*R^2)_R/R
        tempa2 = mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &           
               + mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then
           tempa2 = tempa2 + &
                  (mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR))*(ri_79) 
        end if 
        ! ([u, psi]*R^2)_Z/R
        tempb2 = mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &           
               + mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR)            
        ! (R^2*(u, f'))_R/R^2
        tempc2 = mu(h(:,OP_DR),OP_DRR) + mu(h(:,OP_DZ),OP_DRZ) &           
               + mu(h(:,OP_DRR),OP_DR) + mu(h(:,OP_DRZ),OP_DZ)            
        if(itor.eq.1) then
           tempc2 = tempc2 + &           
                  (mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ))*(2*ri_79)
        endif
        ! (R^2*(u, f'))_Z/R^2
        tempd2 = mu(h(:,OP_DR),OP_DRZ) + mu(h(:,OP_DZ),OP_DZZ) &           
               + mu(h(:,OP_DRZ),OP_DR) + mu(h(:,OP_DZZ),OP_DZ)            

        temp = prod(tempa,tempc2*ri2_79) &
             + prod(tempb,tempd2*ri2_79) &
             + prod(tempc,tempa2*(-ri2_79)) &
             + prod(tempd,tempb2*(-ri2_79)) &
             + prod(mu(-g(:,OP_GS)*ri2_79,OP_DZ),tempd2) & 
             + prod(mu(-g(:,OP_GS)*ri2_79,OP_DR),tempc2) 
     end if

  v3upsif = temp
end function v3upsif

! V3ubf
! =====
!function v3ubf(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v3ubf
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!  integer :: j


!     if(surface_int) then
!        temp = 0.
!     else
!        do j=1, dofs_per_element
!           ! [nu, f']'*R
!           tempa(j,:) = &
!                e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) &
!              + e(j,:,OP_DZP)*h(:,OP_DR) - e(j,:,OP_DRP)*h(:,OP_DZ) 
!           ! (nu,f'')
!           tempb(j,:) = &
!                e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
!        end do
!        ! F*u_GS + (u, F)
!        temp79a = g(:,OP_1)*f(:,OP_GS) & 
!                + g(:,OP_DR)*f(:,OP_DR) + g(:,OP_DZ)*f(:,OP_DZ)
!        ! (R^2(u, f'))_R/R^2
!        temp79b = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
!        if(itor.eq.1) then
!           temp79b = temp79b + 2*ri_79* &
!                 (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
!        end if 
!        ! (R^2(u, f'))_Z/R^2
!        temp79c = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)            
!        ![u,F]*R
!        temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ) 

!        temp = - intx3(tempa,temp79a,ri3_79) &
!               + intx3(tempb,temp79d,ri3_79) &
!               - intx4(e(:,:,OP_DZP),temp79b,g(:,OP_1),ri3_79) &
!               + intx4(e(:,:,OP_DRP),temp79c,g(:,OP_1),ri3_79) &
!               - intx4(e(:,:,OP_DZ),temp79b,g(:,OP_DP),ri3_79) &
!               + intx4(e(:,:,OP_DR),temp79c,g(:,OP_DP),ri3_79) 
!     end if

!  v3ubf = temp
!end function v3ubf

function v3ubf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3ubf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! [nu, f']'*R
           tempa = &
                mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR) &
              + mu(h(:,OP_DR),OP_DZP) + mu(-h(:,OP_DZ),OP_DRP) 
           ! (nu,f'')
           tempb = &
                mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR) 
        ! F*u_GS + (u, F)
        tempa2 = mu(g(:,OP_1),OP_GS) & 
               + mu(g(:,OP_DR),OP_DR) + mu(g(:,OP_DZ),OP_DZ)
        ! (R^2(u, f'))_R/R^2
        tempb2 = mu(h(:,OP_DR),OP_DRR) + mu(h(:,OP_DZ),OP_DRZ) &           
                + mu(h(:,OP_DRR),OP_DR) + mu(h(:,OP_DRZ),OP_DZ)            
        if(itor.eq.1) then
           tempb2 = tempb2 + &
                 (mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ))*(2*ri_79)
        end if 
        ! (R^2(u, f'))_Z/R^2
        tempc2 = mu(h(:,OP_DR),OP_DRZ) + mu(h(:,OP_DZ),OP_DZZ) &           
               + mu(h(:,OP_DRZ),OP_DR) + mu(h(:,OP_DZZ),OP_DZ)            
        ![u,F]*R
        tempd2 = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR) 

        temp =   prod(tempa,tempa2*(-ri3_79)) &
               + prod(tempb,tempd2*ri3_79) &
               + prod(mu(-g(:,OP_1)*ri3_79,OP_DZP),tempb2) &
               + prod(mu( g(:,OP_1)*ri3_79,OP_DRP),tempc2) &
               + prod(mu(-g(:,OP_DP)*ri3_79,OP_DZ),tempb2) &
               + prod(mu( g(:,OP_DP)*ri3_79,OP_DR),tempc2) 
     end if

  v3ubf = temp
end function v3ubf

! V3uff
! =====
!function v3uff(e,f,g,h)
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v3uff
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!  integer :: j


!     if(surface_int) then
!        temp = 0.
!     else
!        do j=1, dofs_per_element
!           ![nu, f']_R*R
!           tempa(j,:) = &
!                e(j,:,OP_DZ)*g(:,OP_DRR) - e(j,:,OP_DR)*g(:,OP_DRZ) &
!              + e(j,:,OP_DRZ)*g(:,OP_DR) - e(j,:,OP_DRR)*g(:,OP_DZ)
!           if(itor.eq.1) then
!              tempa(j,:) = tempa(j,:) - ri_79*&
!                   (e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ)) 
!           end if 
!           ![nu, f']_Z*R
!           tempb(j,:) = &
!                e(j,:,OP_DZ)*g(:,OP_DRZ) - e(j,:,OP_DR)*g(:,OP_DZZ) &
!              + e(j,:,OP_DZZ)*g(:,OP_DR) - e(j,:,OP_DRZ)*g(:,OP_DZ)
!        end do
!        ! (R^2*(u, f'))_R/R^2
!        temp79a = f(:,OP_DRR)*h(:,OP_DR) + f(:,OP_DRZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRR) + f(:,OP_DZ)*h(:,OP_DRZ)            
!        if(itor.eq.1) then
!           temp79a = temp79a + 2*ri_79* &
!                  (f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))
!        endif
!        ! (R^2*(u, f'))_Z/R^2
!        temp79b = f(:,OP_DRZ)*h(:,OP_DR) + f(:,OP_DZZ)*h(:,OP_DZ) &           
!                + f(:,OP_DR)*h(:,OP_DRZ) + f(:,OP_DZ)*h(:,OP_DZZ)     
!        temp = intx3(tempa,temp79a,ri_79) &
!             + intx3(tempb,temp79b,ri_79) 
!     end if

!  v3uff = temp
!end function v3uff

function v3uff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3uff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ![nu, f']_R*R
           tempa = &
                mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR) &
              + mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR)
           if(itor.eq.1) then
              tempa = tempa + &
                   (mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR))*(-ri_79) 
           end if 
           ![nu, f']_Z*R
           tempb = &
                mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR) &
              + mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ)
        ! (R^2*(u, f'))_R/R^2
        tempa2 = mu(h(:,OP_DR),OP_DRR) + mu(h(:,OP_DZ),OP_DRZ) &           
               + mu(h(:,OP_DRR),OP_DR) + mu(h(:,OP_DRZ),OP_DZ)            
        if(itor.eq.1) then
           tempa2 = tempa2 + &
                  (mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ))*(2*ri_79)
        endif
        ! (R^2*(u, f'))_Z/R^2
        tempb2 = mu(h(:,OP_DR),OP_DRZ) + mu(h(:,OP_DZ),OP_DZZ) &           
               + mu(h(:,OP_DRZ),OP_DR) + mu(h(:,OP_DZZ),OP_DZ)     
        temp = prod(tempa,tempa2*ri_79) &
             + prod(tempb,tempb2*ri_79) 
     end if

  v3uff = temp
end function v3uff
#endif

! v3vpsipsi
! =========
!function v3vpsipsi(e,f,g,h)
!!
!!  e trial
!!  f lin
!!  g psi
!!  h psi
!  use basic
!  use m3dc1_nint

!  implicit none

!  vectype, dimension(dofs_per_element) :: v3vpsipsi
!  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!#if defined(USE3D) || defined(USECOMPLEX)
!  vectype, dimension(dofs_per_element) :: temp
!  vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempf
!  integer :: j

!     if(surface_int) then
!        temp = 0.
!     else
!        do j=1, dofs_per_element
!           tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP)
!           tempf(j,:) = f(:,OP_DP)* &
!                (e(j,:,OP_DZ)*g(:,OP_DZ )+e(j,:,OP_DR)*g(:,OP_DR )) &
!                +       f(:,OP_1 )* &
!                (e(j,:,OP_DZ)*g(:,OP_DZP)+e(j,:,OP_DR)*g(:,OP_DRP))
!        end do

!        temp79b = g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ)
!#ifdef USEST
!        temp79e =  f(:,OP_1)*g(:,OP_GS) &
!             + g(:,OP_DZ )*f(:,OP_DZ) + g(:,OP_DR )*f(:,OP_DR)
!#else
!        temp79e = f(:,OP_DP)*g(:,OP_GS) + f(:,OP_1)*g(:,OP_GSP) &
!             + g(:,OP_DZP)*f(:,OP_DZ ) + g(:,OP_DRP)*f(:,OP_DR ) &
!             + g(:,OP_DZ )*f(:,OP_DZP) + g(:,OP_DR )*f(:,OP_DRP)
!#endif
!        temp = intx3(tempa,ri4_79,temp79b) &
!#ifdef USEST
!             + intx4(e(:,:,OP_DZ),ri4_79,temp79e,h(:,OP_DZP)) &
!             + intx4(e(:,:,OP_DR),ri4_79,temp79e,h(:,OP_DRP)) &
!             + intx4(e(:,:,OP_DZP),ri4_79,temp79e,h(:,OP_DZ)) &
!             + intx4(e(:,:,OP_DRP),ri4_79,temp79e,h(:,OP_DR)) &
!#else
!             - intx4(e(:,:,OP_DZ),ri4_79,temp79e,h(:,OP_DZ)) &
!             - intx4(e(:,:,OP_DR),ri4_79,temp79e,h(:,OP_DR)) &
!#endif
!             - intx3(tempf,ri4_79,h(:,OP_GS))
!     end if

!  v3vpsipsi = temp
!#else
!  v3vpsipsi = 0.
!#endif
!end function v3vpsipsi

function v3vpsipsi(g,h)
!
!  e trial
!  f lin
!  g psi
!  h psi
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vpsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp, tempf
  type(muarray) :: tempa, tempb, tempe

     if(surface_int) then
        temp%len = 0
     else
           tempa = mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR)
           tempf = prod(g(:,OP_DZ ),OP_DZ,OP_DP)+prod(g(:,OP_DR ),OP_DR,OP_DP) &
                +  prod(g(:,OP_DZP),OP_DZ,OP_1)+prod(g(:,OP_DRP),OP_DR,OP_1)
   
        tempb = mu(g(:,OP_DZ),OP_DR) + mu(-g(:,OP_DR),OP_DZ)
#ifdef USEST
        tempe =  mu(g(:,OP_GS),OP_1) &
             + mu(g(:,OP_DZ ),OP_DZ) + mu(g(:,OP_DR ),OP_DR)
#else
        tempe = mu(g(:,OP_GS),OP_DP) + mu(g(:,OP_GSP),OP_1) &
             + mu(g(:,OP_DZP),OP_DZ ) + mu(g(:,OP_DRP),OP_DR ) &
             + mu(g(:,OP_DZ ),OP_DZP) + mu(g(:,OP_DR ),OP_DRP)
#endif

        temp = prod(tempa*ri4_79,tempb) &
#ifdef USEST
             + prod(mu(ri4_79*h(:,OP_DZP),OP_DZ),tempe) &
             + prod(mu(ri4_79*h(:,OP_DRP),OP_DR),tempe) &
             + prod(mu(ri4_79*h(:,OP_DZ),OP_DZP),tempe) &
             + prod(mu(ri4_79*h(:,OP_DR),OP_DRP),tempe) &
#else
             + prod(mu(-ri4_79*h(:,OP_DZ),OP_DZ),tempe) &
             + prod(mu(-ri4_79*h(:,OP_DR),OP_DR),tempe) &
#endif
             + tempf*(-ri4_79*h(:,OP_GS))
      end if
  v3vpsipsi = temp
#else
  v3vpsipsi%len = 0
#endif
end function v3vpsipsi



! v3vpsib
! =======
! function v3vpsib(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3vpsib
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb
!   integer :: j

!      if(surface_int) then
!         do j=1, dofs_per_element
!            tempa(j,:) = norm79(:,1)*e(j,:,OP_DR) + norm79(:,2)*e(j,:,OP_DZ)
!            if(itor.eq.1) then
!               tempa(j,:) = tempa(j,:) &
!                    - 2.*ri_79*norm79(:,1)*e(j,:,OP_1)
!            end if
!         end do
!         temp = intx5(tempa,ri3_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
!              - intx5(tempa,ri3_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
!      else
!         temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

!         temp = intx3(e(:,:,OP_GS),ri3_79,temp79a)

! !   scj removed 4/1/2011        
!         if(itor.eq.1) then
!            temp = temp - &
!                 2.*intx3(e(:,:,OP_DR),ri4_79,temp79a)
!         endif

! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         do j=1, dofs_per_element 
!            tempb(j,:) = &
!                 -(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
!                 *(g(:,OP_DRP)*f(:,OP_1) + g(:,OP_DR)*f(:,OP_DP)) &
!                 +(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
!                 *(g(:,OP_DZP)*f(:,OP_1) + g(:,OP_DZ)*f(:,OP_DP))
!         end do
!         temp = temp + intx2(tempb,ri5_79)
! #else
!         do j=1, dofs_per_element
!            tempb(j,:) = f(:,OP_DPP)* &
!                 (e(j,:,OP_DZ)*g(:,OP_DR)-e(j,:,OP_DR)*g(:,OP_DZ)) &
!                 + 2.*f(:,OP_DP)* &
!                 (e(j,:,OP_DZ)*g(:,OP_DRP)-e(j,:,OP_DR)*g(:,OP_DZP)) &
!                 +    f(:,OP_1)* &
!                 (e(j,:,OP_DZ)*g(:,OP_DRPP)-e(j,:,OP_DR)*g(:,OP_DZPP))
!         end do
!         temp = temp + intx3(tempb,ri5_79,h(:,OP_1))
! #endif
! #endif
!      end if

!   v3vpsib = temp
! end function v3vpsib

function v3vpsib(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vpsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp, tempb
  type(muarray) :: tempa

     if(surface_int) then
        tempa = mu(norm79(:,1),OP_DR) + mu(norm79(:,2),OP_DZ)
           if(itor.eq.1) then
              tempa = tempa &
                   + mu(-2.*ri_79*norm79(:,1),OP_1)
           end if
        temp = prod(tempa,mu( ri3_79*g(:,OP_DR)*h(:,OP_1),OP_DZ)) &
             + prod(tempa,mu(-ri3_79*g(:,OP_DZ)*h(:,OP_1),OP_DR))
      else
        tempa = (mu(g(:,OP_DZ),OP_DR) + mu(-g(:,OP_DR),OP_DZ))*h(:,OP_1)

        temp = prod(mu(ri3_79,OP_GS),tempa)

!   scj removed 4/1/2011        
        if(itor.eq.1) then
           temp = temp + &
                prod(mu(-2.*ri4_79,OP_DR),tempa)
        endif

#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
           tempb = &
                 prod(mu(-h(:,OP_DP),OP_DZ) + mu(-h(:,OP_1),OP_DZP), &
                 mu(g(:,OP_DRP),OP_1) + mu(g(:,OP_DR),OP_DP)) &
                +prod(mu(h(:,OP_DP),OP_DR) + mu(h(:,OP_1),OP_DRP), &
                 mu(g(:,OP_DZP),OP_1) + mu(g(:,OP_DZ),OP_DP))
        temp = temp + tempb*ri5_79
#else
           tempb = prod(g(:,OP_DR),OP_DZ,OP_DPP)+prod(-g(:,OP_DZ),OP_DR,OP_DPP) &
                 + prod(2.*g(:,OP_DRP),OP_DZ,OP_DP)+prod(-2.*g(:,OP_DZP),OP_DR,OP_DP) &
                 + prod(g(:,OP_DRPP),OP_DZ,OP_1)+prod(-g(:,OP_DZPP),OP_DR,OP_1)
         temp = temp + tempb*(ri5_79*h(:,OP_1))
#endif
#endif
     end if

  v3vpsib = temp
end function v3vpsib


! V3vbb
! =====
function v3vbb(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = 0.
  end if

  v3vbb = temp
#else
  v3vbb= 0.
#endif

end function v3vbb

#ifdef USE3D
! V3vpsif
! =====
! function v3vpsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3vpsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempd, tempe, tempf
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! [nu, psi']*R
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DRP) - e(j,:,OP_DR)*g(:,OP_DZP) 
!            ! (nu, f'')
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
!            ! (nu, psi)
!            tempc(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR) 
!            ! [nu, f']*R
!            tempd(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ) 
!            ! [nu, f'']*R
!            tempe(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DRP) - e(j,:,OP_DR)*h(:,OP_DZP) 
!            ! [nu', f']*R
!            tempf(j,:) = &
!                 e(j,:,OP_DZP)*h(:,OP_DR) - e(j,:,OP_DRP)*h(:,OP_DZ) 
!         end do

!         ! [v, f']'*R 
!         temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
!         ! [v, psi]*R
!         temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)            
!         ! (v, f') + v*f'_LP
!         temp79c = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
!                 + f(:,OP_1)*h(:,OP_LP)
!         ! (v, psi)
!         temp79d = f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ) 

!         temp = - intx3(tempc,temp79a,ri3_79) &
!                + intx3(tempb,temp79b,ri3_79) &
!                + intx3(tempa,temp79c,ri3_79) &
!                + intx3(tempe,temp79d,ri3_79) &
!                + intx3(tempf,temp79d,ri3_79) &
!                - intx4(tempd,g(:,OP_GS),f(:,OP_DP),ri3_79) &
!                + intx4(tempf,g(:,OP_GS),f(:,OP_1),ri3_79) 
!      end if

!   v3vpsif = temp
! end function v3vpsif

function v3vpsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vpsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd, tempe, tempf
  type(muarray) :: tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! [nu, psi']*R
           tempa = &
                mu(g(:,OP_DRP),OP_DZ) + mu(-g(:,OP_DZP),OP_DR) 
           ! (nu, f'')
           tempb = &
                mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR) 
           ! (nu, psi)
           tempc = &
                mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR) 
           ! [nu, f']*R
           tempd = &
                mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR) 
           ! [nu, f'']*R
           tempe = &
                mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR) 
           ! [nu', f']*R
           tempf = &
                mu(h(:,OP_DR),OP_DZP) + mu(-h(:,OP_DZ),OP_DRP) 

        ! [v, f']'*R 
        tempa2 = mu(h(:,OP_DR),OP_DZP) + mu(-h(:,OP_DZ),OP_DRP) &           
               + mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR)            
        ! [v, psi]*R
        tempb2 = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)            
        ! (v, f') + v*f'_LP
        tempc2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) & 
               + mu(h(:,OP_LP),OP_1)
        ! (v, psi)
        tempd2 = mu(g(:,OP_DR),OP_DR) + mu(g(:,OP_DZ),OP_DZ) 

        temp =   prod(tempc,tempa2*(-ri3_79)) &
               + prod(tempb,tempb2*ri3_79) &
               + prod(tempa,tempc2*ri3_79) &
               + prod(tempe,tempd2*ri3_79) &
               + prod(tempf,tempd2*ri3_79) &
               + prod(tempd,mu(-g(:,OP_GS)*ri3_79,OP_DP)) &
               + prod(tempf,mu( g(:,OP_GS)*ri3_79,OP_1)) 
     end if

  v3vpsif = temp
end function v3vpsif

! V3vbf
! =====
! function v3vbf(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3vbf
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb, tempc
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (nu, f')
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR) 
!            ! (nu, f'')
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
!            ! (nu', f'')
!            tempc(j,:) = &
!                 e(j,:,OP_DZP)*h(:,OP_DZP) + e(j,:,OP_DRP)*h(:,OP_DRP) 
!         end do
!         ! (v, f') + v*f'_LP
!         temp79a = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ) & 
!                 + f(:,OP_1)*h(:,OP_LP)

!         temp = - intx4(tempa,g(:,OP_1),f(:,OP_DPP),ri4_79) &
!                - intx4(tempb,g(:,OP_1),f(:,OP_DP),ri4_79) &
!                + intx4(tempc,g(:,OP_1),f(:,OP_1),ri4_79) &
!                + intx4(e(:,:,OP_GS),g(:,OP_1),temp79a,ri2_79)
!         if(itor.eq.1) then
!            temp = temp &
!                - 2*intx4(e(:,:,OP_DR),g(:,OP_1),temp79a,ri3_79)
!         end if 
!      end if

!   v3vbf = temp
! end function v3vbf

function v3vbf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vbf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempa2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f')
           tempa = &
                mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR) 
           ! (nu, f'')
           tempb = &
                mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR) 
           ! (nu', f'')
           tempc = &
                mu(h(:,OP_DZP),OP_DZP) + mu(h(:,OP_DRP),OP_DRP) 
        ! (v, f') + v*f'_LP
        tempa2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) & 
               + mu(h(:,OP_LP),OP_1)

        temp =   prod(tempa,mu(-g(:,OP_1)*ri4_79,OP_DPP)) &
               + prod(tempb,mu(-g(:,OP_1)*ri4_79,OP_DP)) &
               + prod(tempc,mu(g(:,OP_1)*ri4_79,OP_1)) &
               + prod(mu(g(:,OP_1)*ri2_79,OP_GS),tempa2)
        if(itor.eq.1) then
           temp = temp &
                + prod(mu(-2*g(:,OP_1)*ri3_79,OP_DR),tempa2)
        end if 
     end if

  v3vbf = temp
end function v3vbf

! V3vff
! =====
! function v3vff(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3vff
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! (nu, f'')
!            tempa(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DZP) + e(j,:,OP_DR)*g(:,OP_DRP) 
!            ! [nu, f']*R
!            tempb(j,:) = &
!                 e(j,:,OP_DZ)*g(:,OP_DR) - e(j,:,OP_DR)*g(:,OP_DZ) 
!         end do
!         ! [v, f']'*R
!         temp79a = f(:,OP_DZP)*h(:,OP_DR) - f(:,OP_DRP)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRP) - f(:,OP_DR)*h(:,OP_DZP)            
!         ! (v, f') 
!         temp79b = f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ)  

!         temp = - intx3(tempa,temp79b,ri2_79) &
!                - intx3(tempb,temp79a,ri2_79) 
!      end if

!   v3vff = temp
! end function v3vff

function v3vff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! (nu, f'')
           tempa = &
                mu(g(:,OP_DZP),OP_DZ) + mu(g(:,OP_DRP),OP_DR) 
           ! [nu, f']*R
           tempb = &
                mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR) 
        ! [v, f']'*R
        tempa2 = mu(h(:,OP_DR),OP_DZP) + mu(-h(:,OP_DZ),OP_DRP) &           
               + mu(h(:,OP_DRP),OP_DZ) + mu(-h(:,OP_DZP),OP_DR)            
        ! (v, f') 
        tempb2 = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ)  

        temp = prod(tempa,tempb2*(-ri2_79)) &
             + prod(tempb,tempa2*(-ri2_79)) 
     end if

  v3vff = temp
end function v3vff
#endif

! V3chipsipsi
! ===========
!function v3chipsipsi(e,f,g,h)
  !use basic
  !use arrays
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3chipsipsi
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

  !vectype, dimension(dofs_per_element) :: temp
  !vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempd, tempe, tempf
  !integer :: j

  !! <f,g>,r
  !temp79b = f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DRZ)*g(:,OP_DZ ) &
       !+    f(:,OP_DR )*g(:,OP_DRR) + f(:,OP_DZ )*g(:,OP_DRZ)
  
  !! <f,g>,z
  !temp79c = f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZZ)*g(:,OP_DZ ) &
       !+    f(:,OP_DR )*g(:,OP_DRZ) + f(:,OP_DZ )*g(:,OP_DZZ)

  !do j=1, dofs_per_element
     !! <e,h>,r
     !tempe(j,:) = e(j,:,OP_DRR)*h(:,OP_DR ) + e(j,:,OP_DRZ)*h(:,OP_DZ ) &
          !+       e(j,:,OP_DR )*h(:,OP_DRR) + e(j,:,OP_DZ )*h(:,OP_DRZ)
  
     !! <e,h>,z
     !tempf(j,:) = e(j,:,OP_DRZ)*h(:,OP_DR ) + e(j,:,OP_DZZ)*h(:,OP_DZ ) &
          !+       e(j,:,OP_DR )*h(:,OP_DRZ) + e(j,:,OP_DZ )*h(:,OP_DZZ)
  !end do
  
     !if(surface_int) then
        !do j=1, dofs_per_element
           !tempa(j,:) = e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
        !end do
        !temp = &
             !- intx5(tempa,ri6_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
             !- intx5(tempa,ri6_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
             !- intx5(tempa,ri6_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
             !- intx5(tempa,ri6_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
             !- intx5(tempa,ri6_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
             !- intx5(tempa,ri6_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
             !- intx5(tempa,ri6_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
             !- intx5(tempa,ri6_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
        !if(itor.eq.1) then
           !temp = temp + 2.* &
                !(intx5(tempa,ri7_79,f(:,OP_DZ),g(:,OP_DZ),norm79(:,1)) &
                !+intx5(tempa,ri7_79,f(:,OP_DR),g(:,OP_DR),norm79(:,1)))
        !endif
     !else
        !temp = intx3(tempe,ri6_79,temp79b) &
             !+ intx3(tempf,ri6_79,temp79c) &
             !- intx4(e(:,:,OP_DZ),ri6_79,temp79c,h(:,OP_GS)) &
             !- intx4(e(:,:,OP_DR),ri6_79,temp79b,h(:,OP_GS))

        !if(itor.eq.1) then
           !temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

           !do j=1, dofs_per_element
              !tempd(j,:) = e(j,:,OP_DZ)*h(:,OP_DZ) + e(j,:,OP_DR)*h(:,OP_DR)
           !end do
           
           !temp = temp &
                !-2.*intx3(tempd,ri7_79,temp79b) &
                !-2.*intx3(tempe,ri7_79,temp79a) &
                !+2.*intx4(e(:,:,OP_DR),ri7_79,temp79a,h(:,OP_GS)) &
                !+4.*intx3(tempd,ri8_79,temp79a)
        !endif
     !end if

  !v3chipsipsi = temp
!end function v3chipsipsi

function v3chipsipsi(g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chipsipsi
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempc, tempd, tempe, tempf

  ! <f,g>,r
  tempb = mu(g(:,OP_DR ),OP_DRR) + mu(g(:,OP_DZ ),OP_DRZ) &
       +  mu(g(:,OP_DRR),OP_DR)  + mu(g(:,OP_DRZ),OP_DZ)
 
  ! <f,g>,z
  tempc = mu(g(:,OP_DR ),OP_DRZ) + mu(g(:,OP_DZ ),OP_DZZ) &
       +  mu(g(:,OP_DRZ),OP_DR) +  mu(g(:,OP_DZZ),OP_DZ)

     ! <e,h>,r
     tempe = mu(h(:,OP_DR ),OP_DRR) + mu(h(:,OP_DZ ),OP_DRZ) &
          +  mu(h(:,OP_DRR),OP_DR) +  mu(h(:,OP_DRZ),OP_DZ)
  
     ! <e,h>,z
     tempf = mu(h(:,OP_DR ),OP_DRZ) + mu(h(:,OP_DZ ),OP_DZZ) &
          +  mu(h(:,OP_DRZ),OP_DR) +  mu(h(:,OP_DZZ),OP_DZ)
    
     if(surface_int) then
           tempa = mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
        temp = &
               prod(tempa,mu(-ri6_79*norm79(:,1)*g(:,OP_DZ ),OP_DRZ)) &
             + prod(tempa,mu(-ri6_79*norm79(:,1)*g(:,OP_DR ),OP_DRR)) &
             + prod(tempa,mu(-ri6_79*norm79(:,1)*g(:,OP_DRZ),OP_DZ)) &
             + prod(tempa,mu(-ri6_79*norm79(:,1)*g(:,OP_DRR),OP_DR)) &
             + prod(tempa,mu(-ri6_79*norm79(:,2)*g(:,OP_DZ ),OP_DZZ)) &
             + prod(tempa,mu(-ri6_79*norm79(:,2)*g(:,OP_DR ),OP_DRZ)) &
             + prod(tempa,mu(-ri6_79*norm79(:,2)*g(:,OP_DZZ),OP_DZ)) &
             + prod(tempa,mu(-ri6_79*norm79(:,2)*g(:,OP_DRZ),OP_DR))
         if(itor.eq.1) then
           temp = temp + &
                (prod(tempa,mu(2.*ri7_79*g(:,OP_DZ)*norm79(:,1),OP_DZ)) &
                +prod(tempa,mu(2.*ri7_79*g(:,OP_DR)*norm79(:,1),OP_DR)))
         endif
     else
        temp = prod(tempe*ri6_79,tempb) &
             + prod(tempf*ri6_79,tempc) &
             + prod(mu(-ri6_79*h(:,OP_GS),OP_DZ),tempc) &
             + prod(mu(-ri6_79*h(:,OP_GS),OP_DR),tempb)

        if(itor.eq.1) then
           tempa = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)

              tempd = mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
           
           temp = temp &
                + prod(tempd*(-2.*ri7_79),tempb) &
                + prod(tempe*(-2.*ri7_79),tempa) &
                + prod(mu(2.*ri7_79*h(:,OP_GS),OP_DR),tempa) &
                + prod(tempd*( 4.*ri8_79),tempa)
        endif
     end if

  v3chipsipsi = temp
end function v3chipsipsi


! V3chipsib
! =========
!function v3chipsib(e,f,g,h)
  !use basic
  !use arrays
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: v3chipsib
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!#if defined(USE3D) || defined(USECOMPLEX)
  !vectype, dimension(dofs_per_element) :: temp
  !vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb, tempd
  !integer :: j

     !if(surface_int) then
        !do j=1, dofs_per_element
           !tempa(j,:) = h(:,OP_1)* &
                !(norm79(:,1)*e(j,:,OP_DZ) - norm79(:,2)*e(j,:,OP_DR))
        !end do
        !temp = intx4(tempa,ri7_79,f(:,OP_DZP),g(:,OP_DZ )) &
             !+ intx4(tempa,ri7_79,f(:,OP_DRP),g(:,OP_DR )) &
             !+ intx4(tempa,ri7_79,f(:,OP_DZ ),g(:,OP_DZP)) &
             !+ intx4(tempa,ri7_79,f(:,OP_DR ),g(:,OP_DRP))
     !else
        !temp79a = h(:,OP_1)*f(:,OP_GS) &
             !+ h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
        !temp79c = f(:,OP_DZP)*h(:,OP_DR ) - f(:,OP_DRP)*h(:,OP_DZ ) &
             !+    f(:,OP_DZ )*h(:,OP_DRP) - f(:,OP_DR )*h(:,OP_DZP)
        !if(itor.eq.1) then
           !temp79a = temp79a - 2.*ri_79*f(:,OP_DR)*h(:,OP_1)
           !temp79c = temp79c - 4.*ri_79* &
                !(f(:,OP_DZP)*h(:,OP_1) + f(:,OP_DZ)*h(:,OP_DP))
        !endif


        !do j=1, dofs_per_element
           !tempb(j,:) = e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)
           !tempd(j,:) = (f(:,OP_DZ )* &
                !e(j,:,OP_DR)-f(:,OP_DR )*e(j,:,OP_DZ))*h(:,OP_DP) &
                !+    (f(:,OP_DZP)* &
                !e(j,:,OP_DR)-f(:,OP_DRP)*e(j,:,OP_DZ))*h(:,OP_1 )
           !if(itor.eq.1) then
              !tempb(j,:) = tempb(j,:) - 4.*ri_79*e(j,:,OP_DZ)*h(:,OP_1)
           !end if
        !end do

        !temp = intx4(e(:,:,OP_DZ),ri7_79,g(:,OP_DRP),temp79a) &
             !- intx4(e(:,:,OP_DR),ri7_79,g(:,OP_DZP),temp79a) &
             !- intx4(tempb,ri7_79,f(:,OP_DZP),g(:,OP_DZ)) &
             !- intx4(tempb,ri7_79,f(:,OP_DRP),g(:,OP_DR)) &
             !- intx4(tempb,ri7_79,f(:,OP_DZ),g(:,OP_DZP)) &
             !- intx4(tempb,ri7_79,f(:,OP_DR),g(:,OP_DRP)) &
             !+ intx4(e(:,:,OP_DZ),ri7_79,temp79c,g(:,OP_DZ)) &
             !+ intx4(e(:,:,OP_DR),ri7_79,temp79c,g(:,OP_DR)) &
             !+ intx3(tempd,ri7_79,g(:,OP_GS))
     !end if

  !v3chipsib = temp
!#else
  !v3chipsib = 0.
!#endif
!end function v3chipsib

function v3chipsib(g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chipsib
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp, tempd
  type(muarray) :: tempa, tempb, tempc

     if(surface_int) then
        tempa = &
                (mu(norm79(:,1),OP_DZ) + mu(-norm79(:,2),OP_DR))*h(:,OP_1)
        temp = prod(tempa,mu(ri7_79*g(:,OP_DZ ),OP_DZP)) &
             + prod(tempa,mu(ri7_79*g(:,OP_DR ),OP_DRP)) &
             + prod(tempa,mu(ri7_79*g(:,OP_DZP),OP_DZ)) &
             + prod(tempa,mu(ri7_79*g(:,OP_DRP),OP_DR))
      else
        tempa = mu(h(:,OP_1),OP_GS) &
              + mu(h(:,OP_DZ),OP_DZ) + mu(h(:,OP_DR),OP_DR)
        tempc = mu(h(:,OP_DR ),OP_DZP) + mu(-h(:,OP_DZ ),OP_DRP) &
             +  mu(h(:,OP_DRP),OP_DZ)  + mu(-h(:,OP_DZP),OP_DR)
        if(itor.eq.1) then
           tempa = tempa + mu(-2.*ri_79*h(:,OP_1),OP_DR)
           tempc = tempc + &
                (mu(h(:,OP_1),OP_DZP) + mu(h(:,OP_DP),OP_DZ))*(-4.*ri_79)
        endif

           tempb = mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR)
           tempd = prod(h(:,OP_DP),OP_DR,OP_DZ )+prod(-h(:,OP_DP),OP_DZ,OP_DR) &
                +  prod(h(:,OP_1),OP_DR,OP_DZP) +prod(-h(:,OP_1),OP_DZ,OP_DRP)
           if(itor.eq.1) then
              tempb = tempb + mu(-4.*ri_79*h(:,OP_1),OP_DZ)
           end if
          
        temp = prod(mu( ri7_79*g(:,OP_DRP),OP_DZ),tempa) &
             + prod(mu(-ri7_79*g(:,OP_DZP),OP_DR),tempa) &
             + prod(tempb,mu(-ri7_79*g(:,OP_DZ),OP_DZP)) &
             + prod(tempb,mu(-ri7_79*g(:,OP_DR),OP_DRP)) &
             + prod(tempb,mu(-ri7_79*g(:,OP_DZP),OP_DZ)) &
             + prod(tempb,mu(-ri7_79*g(:,OP_DRP),OP_DR)) &
             + prod(mu(ri7_79*g(:,OP_DZ),OP_DZ),tempc) &
             + prod(mu(ri7_79*g(:,OP_DR),OP_DR),tempc) &
             + tempd*(ri7_79*g(:,OP_GS))
       end if
  v3chipsib = temp
#else
  v3chipsib%len = 0
#endif
end function v3chipsib


! V3chibb
! =======
! function v3chibb(e,f,g,h)
!   use basic
!   use arrays
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3chibb
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element, MAX_PTS) :: tempa, tempb
!   integer :: j

!      if(surface_int) then
!         do j=1, dofs_per_element
!            tempa(j,:) = h(:,OP_1)* &
!                 (norm79(:,1)*e(j,:,OP_DR) + norm79(:,2)*e(j,:,OP_DZ))
!            !        if(itor.eq.1) then
!            !           temp79a = temp79a &
!            !                - 2.*ri_79*norm79(:,1)*e(:,OP_1)*h(:,OP_1)
!            !        end if
!         end do
!         temp = &
!              - intx4(tempa,ri6_79,f(:,OP_GS),g(:,OP_1)) &
!              - intx4(tempa,ri6_79,f(:,OP_DZ),g(:,OP_DZ)) &
!              - intx4(tempa,ri6_79,f(:,OP_DR),g(:,OP_DR))
!         if(itor.eq.1) then
!            temp = temp + 2.*intx4(tempa,ri7_79,f(:,OP_DR),g(:,OP_1))
!         endif
!      else
!         temp79a = g(:,OP_1)*f(:,OP_GS) &
!              + g(:,OP_DZ)*f(:,OP_DZ) + g(:,OP_DR)*f(:,OP_DR)
        
!         if(itor.eq.1) then
!            temp79a = temp79a - &
!                 2.*ri_79*f(:,OP_DR)*g(:,OP_1)
!         end if
        
!         temp = intx4(e(:,:,OP_GS),ri6_79,temp79a,h(:,OP_1))

! !    scj removed 4/2/2011        
!         if(itor.eq.1) then
!            temp = temp - &
!                 2.*intx4(e(:,:,OP_DR),ri7_79,temp79a,h(:,OP_1))
!         endif
        
! #if defined(USE3D) || defined(USECOMPLEX)
! #ifdef USEST
!         do j=1, dofs_per_element 
!            tempb(j,:) = &
!                 +(e(j,:,OP_DZ)*h(:,OP_DP) + e(j,:,OP_DZP)*h(:,OP_1)) &
!                 *(f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) &
!                 +(e(j,:,OP_DR)*h(:,OP_DP) + e(j,:,OP_DRP)*h(:,OP_1)) &
!                 *(f(:,OP_DRP)*g(:,OP_1) + f(:,OP_DR)*g(:,OP_DP))
!         end do
!         temp = temp + intx2(tempb,ri8_79)
! #else
!         do j=1, dofs_per_element
!            tempb(j,:) = &
!                 (e(j,:,OP_DZ)*f(:,OP_DZPP) + e(j,:,OP_DR)*f(:,OP_DRPP))*g(:,OP_1  ) &
!                 + 2.*(e(j,:,OP_DZ)*f(:,OP_DZP ) + e(j,:,OP_DR)*f(:,OP_DRP ))*g(:,OP_DP ) &
!                 +    (e(j,:,OP_DZ)*f(:,OP_DZ  ) + e(j,:,OP_DR)*f(:,OP_DR  ))*g(:,OP_DPP)
!         end do
!         temp = temp - intx3(tempb,ri8_79,h(:,OP_1))
! #endif
! #endif
!      end if

!   v3chibb = temp
! end function v3chibb

function v3chibb(g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chibb
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h

  type(prodarray) :: temp, tempb
  type(muarray) :: tempa

     if(surface_int) then
          tempa = &
                (mu(norm79(:,1),OP_DR) + mu(norm79(:,2),OP_DZ))*h(:,OP_1)
        temp = &
               prod(tempa,mu(-ri6_79*g(:,OP_1),OP_GS)) &
             + prod(tempa,mu(-ri6_79*g(:,OP_DZ),OP_DZ)) &
             + prod(tempa,mu(-ri6_79*g(:,OP_DR),OP_DR))
        if(itor.eq.1) then
           temp = temp + prod(tempa,mu(2.*ri7_79*g(:,OP_1),OP_DR))
        endif
     else
        tempa = mu(g(:,OP_1),OP_GS) &
             +  mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
         
        if(itor.eq.1) then
           tempa = tempa + &
                mu(-2.*ri_79*g(:,OP_1),OP_DR)
        end if
        
        temp = prod(mu(ri6_79*h(:,OP_1),OP_GS),tempa)

!    scj removed 4/2/2011        
        if(itor.eq.1) then
           temp = temp + &
                prod(mu(-2.*ri7_79*h(:,OP_1),OP_DR),tempa)
        endif
 
#if defined(USE3D) || defined(USECOMPLEX)
#ifdef USEST
           tempb = &
                 prod(mu(h(:,OP_DP),OP_DZ) + mu(h(:,OP_1),OP_DZP), &
                 mu(g(:,OP_1),OP_DZP) + mu(g(:,OP_DP),OP_DZ)) &
                +prod(mu(h(:,OP_DP),OP_DR) + mu(h(:,OP_1),OP_DRP), &
                 mu(g(:,OP_1),OP_DRP) + mu(g(:,OP_DP),OP_DR))
        temp = temp + tempb*ri8_79
#else
           tempb = &
                prod(g(:,OP_1),OP_DZ,OP_DZPP) + prod(g(:,OP_1),OP_DR,OP_DRPP) &
                + prod(2.*g(:,OP_DP),OP_DZ,OP_DZP ) + prod(2.*g(:,OP_DP),OP_DR,OP_DRP ) &
                + prod(g(:,OP_DPP),OP_DZ,OP_DZ  ) + prod(g(:,OP_DPP),OP_DR,OP_DR  )
        temp = temp + tempb*(-ri8_79*h(:,OP_1))
#endif
#endif
     end if

  v3chibb = temp
end function v3chibb

#ifdef USE3D
! V3chipsif
! =====
! function v3chipsif(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3chipsif
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempc, tempd
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! [nu, f']_R * R
!            tempa(j,:) = &
!                 e(j,:,OP_DRZ)*h(:,OP_DR) - e(j,:,OP_DRR)*h(:,OP_DZ) &
!               + e(j,:,OP_DZ)*h(:,OP_DRR) - e(j,:,OP_DR)*h(:,OP_DRZ)
!            if(itor.eq.1) then 
!               tempa(j,:) = tempa(j,:) - ri_79* &
!                   (e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)) 
!            end if
!            ! [nu, f']_Z * R
!            tempb(j,:) = &
!                 e(j,:,OP_DZZ)*h(:,OP_DR) - e(j,:,OP_DRZ)*h(:,OP_DZ) &
!               + e(j,:,OP_DZ)*h(:,OP_DRZ) - e(j,:,OP_DR)*h(:,OP_DZZ)
!            ! ((nu, psi)/R^2)_R * R^2
!            tempc(j,:) = &
!                 e(j,:,OP_DRZ)*g(:,OP_DZ) + e(j,:,OP_DRR)*g(:,OP_DR) & 
!               + e(j,:,OP_DZ)*g(:,OP_DRZ) + e(j,:,OP_DR)*g(:,OP_DRR) 
!            if(itor.eq.1) then 
!               tempc(j,:) = tempc(j,:) - 2*ri_79* &
!                   (e(j,:,OP_DZ)*g(:,OP_DZ) + e(j,:,OP_DR)*g(:,OP_DR)) 
!            end if
!            ! ((nu, psi)/R^2)_Z * R^2
!            tempd(j,:) = &
!                 e(j,:,OP_DZZ)*g(:,OP_DZ) + e(j,:,OP_DRZ)*g(:,OP_DR) & 
!               + e(j,:,OP_DZ)*g(:,OP_DZZ) + e(j,:,OP_DR)*g(:,OP_DRZ) 
!         end do

!         ! [chi, f']_R*R 
!         temp79a = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
!         if(itor.eq.1) then 
!            temp79a = temp79a - ri_79* &        
!                 (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
!         end if 
!         ! [chi, f']_Z*R 
!         temp79b = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
!         ! ((chi, psi)/R^2)_R * R^2
!         temp79c = f(:,OP_DRR)*g(:,OP_DR) + f(:,OP_DRZ)*g(:,OP_DZ) &
!                 + f(:,OP_DR)*g(:,OP_DRR) + f(:,OP_DZ)*g(:,OP_DRZ) 
!         if(itor.eq.1) then 
!            temp79c = temp79c - 2*ri_79* & 
!                  (f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ)) 
!         end if 
!         ! ((chi, psi)/R^2)_Z * R^2
!         temp79d = f(:,OP_DRZ)*g(:,OP_DR) + f(:,OP_DZZ)*g(:,OP_DZ) &
!                 + f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZ)*g(:,OP_DZZ) 

!         temp = intx3(tempa,temp79c,ri5_79) &
!              + intx3(tempb,temp79d,ri5_79) &
!              + intx3(tempc,temp79a,ri5_79) &
!              + intx3(tempd,temp79b,ri5_79) &
!              - intx4(e(:,:,OP_DZ),temp79b,ri5_79,g(:,OP_GS)) &
!              - intx4(e(:,:,OP_DR),temp79a,ri5_79,g(:,OP_GS)) 
!      end if

!   v3chipsif = temp
! end function v3chipsif

function v3chipsif(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chipsif
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb ,tempc, tempd
  type(muarray) :: tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! [nu, f']_R * R
           tempa = &
                mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR) &
              + mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR)
           if(itor.eq.1) then 
              tempa = tempa + &
                  (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79)
           end if
           ! [nu, f']_Z * R
           tempb = &
                mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ) &
              + mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR)
           ! ((nu, psi)/R^2)_R * R^2
           tempc = &
                mu(g(:,OP_DZ),OP_DRZ) + mu(g(:,OP_DR),OP_DRR) & 
              + mu(g(:,OP_DRZ),OP_DZ) + mu(g(:,OP_DRR),OP_DR) 
           if(itor.eq.1) then 
              tempc = tempc + &
                  (mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR))*(-2*ri_79) 
           end if
           ! ((nu, psi)/R^2)_Z * R^2
           tempd = &
                mu(g(:,OP_DZ),OP_DZZ) + mu(g(:,OP_DR),OP_DRZ) & 
              + mu(g(:,OP_DZZ),OP_DZ) + mu(g(:,OP_DRZ),OP_DR) 

        ! [chi, f']_R*R 
        tempa2 = mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR) &           
               + mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then 
           tempa2 = tempa2 + &        
                (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79) 
        end if 
        ! [chi, f']_Z*R 
        tempb2 = mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ) &           
                + mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR)            
        ! ((chi, psi)/R^2)_R * R^2
        tempc2 = mu(g(:,OP_DR),OP_DRR) + mu(g(:,OP_DZ),OP_DRZ) &
               + mu(g(:,OP_DRR),OP_DR) + mu(g(:,OP_DRZ),OP_DZ) 
        if(itor.eq.1) then 
           tempc2 = tempc2 + & 
                 (mu(g(:,OP_DR),OP_DR) + mu(g(:,OP_DZ),OP_DZ))*(-2*ri_79) 
        end if 
        ! ((chi, psi)/R^2)_Z * R^2
        tempd2 = mu(g(:,OP_DR),OP_DRZ) + mu(g(:,OP_DZ),OP_DZZ) &
                + mu(g(:,OP_DRZ),OP_DR) + mu(g(:,OP_DZZ),OP_DZ) 

        temp = prod(tempa,tempc2*ri5_79) &
             + prod(tempb,tempd2*ri5_79) &
             + prod(tempc,tempa2*ri5_79) &
             + prod(tempd,tempb2*ri5_79) &
             + prod(mu(-ri5_79*g(:,OP_GS),OP_DZ),tempb2) &
             + prod(mu(-ri5_79*g(:,OP_GS),OP_DR),tempa2) 
     end if

  v3chipsif = temp
end function v3chipsif

! V3chibf
! =====
! function v3chibf(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3chibf
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j


!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! [nu, f']*R
!            tempa(j,:) = &
!               + e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ) 
!            ! (nu, f'')
!            tempb(j,:) = &
!               + e(j,:,OP_DZ)*h(:,OP_DZP) + e(j,:,OP_DR)*h(:,OP_DRP) 
!         end do

!         ! [chi, F/R^4]'*R^5
!         temp79a = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP) 
!         if(itor.eq.1) then 
!            temp79a = temp79a - 4*ri_79* &        
!                 (f(:,OP_DZP)*g(:,OP_1) + f(:,OP_DZ)*g(:,OP_DP)) 
!         end if 
!         ! [chi, f']_R*R 
!         temp79b = f(:,OP_DRZ)*h(:,OP_DR) - f(:,OP_DRR)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRR) - f(:,OP_DR)*h(:,OP_DRZ)            
!         if(itor.eq.1) then 
!            temp79b = temp79b - ri_79* &        
!                 (f(:,OP_DZ)*h(:,OP_DR) - f(:,OP_DR)*h(:,OP_DZ)) 
!         end if 
!         ! [chi, f']_Z*R 
!         temp79c = f(:,OP_DZZ)*h(:,OP_DR) - f(:,OP_DRZ)*h(:,OP_DZ) &           
!                 + f(:,OP_DZ)*h(:,OP_DRZ) - f(:,OP_DR)*h(:,OP_DZZ)            
!         ! (chi, F/R^2)*R^2 + F*chi_GS
!         temp79d = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
!                 + f(:,OP_GS)*g(:,OP_1) 
!         if(itor.eq.1) then 
!            temp79d = temp79d - 2*ri_79*f(:,OP_DR)*g(:,OP_1) 
!         end if 

!         temp =   intx3(tempa,temp79a,ri6_79) &
!                - intx3(tempb,temp79d,ri6_79) &
!                - intx4(e(:,:,OP_DZ),temp79b,ri6_79,g(:,OP_DP)) &
!                + intx4(e(:,:,OP_DR),temp79c,ri6_79,g(:,OP_DP)) &
!                - intx4(e(:,:,OP_DZP),temp79b,ri6_79,g(:,OP_1)) &
!                + intx4(e(:,:,OP_DRP),temp79c,ri6_79,g(:,OP_1)) 
!      end if

!   v3chibf = temp
! end function v3chibf

function v3chibf(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chibf
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2, tempc2, tempd2
  integer :: j


     if(surface_int) then
        temp%len = 0
     else
           ! [nu, f']*R
           tempa = &
                mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR) 
           ! (nu, f'')
           tempb = &
                mu(h(:,OP_DZP),OP_DZ) + mu(h(:,OP_DRP),OP_DR) 

        ! [chi, F/R^4]'*R^5
        tempa2 = mu(g(:,OP_DR),OP_DZP) + mu(-g(:,OP_DZ),OP_DRP) &           
               + mu(g(:,OP_DRP),OP_DZ) + mu(-g(:,OP_DZP),OP_DR) 
        if(itor.eq.1) then 
           tempa2 = tempa2 + &        
                (mu(g(:,OP_1),OP_DZP) + mu(g(:,OP_DP),OP_DZ))*(-4*ri_79) 
        end if 
        ! [chi, f']_R*R 
        tempb2 = mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR) &           
               + mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then 
           tempb2 = tempb2 + &        
                (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79) 
        end if 
        ! [chi, f']_Z*R 
        tempc2 = mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ) &           
               + mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR)            
        ! (chi, F/R^2)*R^2 + F*chi_GS
        tempd2 = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR) &
               + mu(g(:,OP_1),OP_GS) 
        if(itor.eq.1) then 
           tempd2 = tempd2 + mu(-2*ri_79*g(:,OP_1),OP_DR) 
        end if 

        temp =   prod(tempa,tempa2*ri6_79) &
               + prod(tempb,tempd2*(-ri6_79)) &
               + prod(mu(-ri6_79*g(:,OP_DP),OP_DZ),tempb2) &
               + prod(mu( ri6_79*g(:,OP_DP),OP_DR),tempc2) &
               + prod(mu(-ri6_79*g(:,OP_1),OP_DZP),tempb2) &
               + prod(mu( ri6_79*g(:,OP_1),OP_DRP),tempc2) 
     end if

  v3chibf = temp
end function v3chibf

! V3chiff
! =====
! function v3chiff(e,f,g,h)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: v3chiff
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

!   vectype, dimension(dofs_per_element) :: temp
!   vectype, dimension(dofs_per_element,MAX_PTS) :: tempa, tempb
!   integer :: j

!      if(surface_int) then
!         temp = 0.
!      else
!         do j=1, dofs_per_element
!            ! [nu, f']_R*R
!            tempa(j,:) = &
!                 e(j,:,OP_DRZ)*h(:,OP_DR) - e(j,:,OP_DRR)*h(:,OP_DZ) & 
!               + e(j,:,OP_DZ)*h(:,OP_DRR) - e(j,:,OP_DR)*h(:,OP_DRZ) 
!            if(itor.eq.1) then 
!               tempa(j,:) = tempa(j,:) - ri_79* &
!                   (e(j,:,OP_DZ)*h(:,OP_DR) - e(j,:,OP_DR)*h(:,OP_DZ)) 
!            end if
!            ! [nu, f']_Z*R
!            tempb(j,:) = &
!                 e(j,:,OP_DZZ)*h(:,OP_DR) - e(j,:,OP_DRZ)*h(:,OP_DZ) & 
!               + e(j,:,OP_DZ)*h(:,OP_DRZ) - e(j,:,OP_DR)*h(:,OP_DZZ) 
!         end do

!         ! [chi, f']_R*R
!         temp79a = f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRR) - f(:,OP_DR)*g(:,OP_DRZ)            
!         if(itor.eq.1) then 
!            temp79a = temp79a - ri_79* &        
!                 (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) 
!         end if 
!         ! [chi, f']_Z*R 
!         temp79b = f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ) &           
!                 + f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ)            

!         temp = intx3(tempa,temp79a,ri4_79) &
!              + intx3(tempb,temp79b,ri4_79) 
!      end if

!   v3chiff = temp
! end function v3chiff

function v3chiff(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chiff
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h

  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempa2, tempb2
  integer :: j

     if(surface_int) then
        temp%len = 0
     else
           ! [nu, f']_R*R
           tempa = &
                mu(h(:,OP_DR),OP_DRZ) + mu(-h(:,OP_DZ),OP_DRR) & 
              + mu(h(:,OP_DRR),OP_DZ) + mu(-h(:,OP_DRZ),OP_DR) 
           if(itor.eq.1) then 
              tempa = tempa + &
                  (mu(h(:,OP_DR),OP_DZ) + mu(-h(:,OP_DZ),OP_DR))*(-ri_79) 
           end if
           ! [nu, f']_Z*R
           tempb = &
                mu(h(:,OP_DR),OP_DZZ) + mu(-h(:,OP_DZ),OP_DRZ) & 
              + mu(h(:,OP_DRZ),OP_DZ) + mu(-h(:,OP_DZZ),OP_DR) 

        ! [chi, f']_R*R
        tempa2 = mu(g(:,OP_DR),OP_DRZ) + mu(-g(:,OP_DZ),OP_DRR) &           
               + mu(g(:,OP_DRR),OP_DZ) + mu(-g(:,OP_DRZ),OP_DR)            
        if(itor.eq.1) then 
           tempa2 = tempa2 + &        
                (mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR))*(-ri_79) 
        end if 
        ! [chi, f']_Z*R 
        tempb2 = mu(g(:,OP_DR),OP_DZZ) + mu(-g(:,OP_DZ),OP_DRZ) &           
                + mu(g(:,OP_DRZ),OP_DZ) + mu(-g(:,OP_DZZ),OP_DR)            

        temp = prod(tempa,tempa2*ri4_79) &
             + prod(tempb,tempb2*ri4_79) 
     end if

  v3chiff = temp
end function v3chiff
#endif

! V3uun
! =====
function v3uun(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3uun
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3uun = 0.
     return
  end if

     temp = intx4(e(:,:,OP_DZ),f(:,OP_DR),g(:,OP_DRZ),h(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_DRR),h(:,OP_1)) &
          + intx4(e(:,:,OP_DR),f(:,OP_DZ),g(:,OP_DRZ),h(:,OP_1)) &
          - intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_DZZ),h(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_LP),h(:,OP_1)) &
          - intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_LP),h(:,OP_1))
     if(itor.eq.1) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
        temp = temp + intx4(e(:,:,OP_DZ),ri_79,temp79a,h(:,OP_1))
     endif

  v3uun = temp
end function v3uun

function v3uun1(g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3uun1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h
  type(prodarray) :: temp
  type(muarray) :: tempa

  if(inertia.eq.0 .or. surface_int) then
     v3uun1%len = 0
     return
  end if

     temp = prod( g(:,OP_DRZ)*h(:,OP_1),OP_DZ,OP_DR) &
          + prod(-g(:,OP_DRR)*h(:,OP_1),OP_DZ,OP_DZ) &
          + prod( g(:,OP_DRZ)*h(:,OP_1),OP_DR,OP_DZ) &
          + prod(-g(:,OP_DZZ)*h(:,OP_1),OP_DR,OP_DR) &
          + prod(-g(:,OP_LP)*h(:,OP_1),OP_DZ,OP_DZ) &
          + prod(-g(:,OP_LP)*h(:,OP_1),OP_DR,OP_DR)
     if(itor.eq.1) then
        tempa = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
        temp = temp + prod(mu(ri_79*h(:,OP_1),OP_DZ),tempa)
     endif

  v3uun1 = temp
end function v3uun1

function v3uun2(f,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3uun2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  type(prodarray) :: temp
  type(muarray) :: tempa

  if(inertia.eq.0 .or. surface_int) then
     v3uun2%len = 0
     return
  end if

     temp = prod( f(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DRZ) &
          + prod(-f(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DRR) &
          + prod( f(:,OP_DZ)*h(:,OP_1),OP_DR,OP_DRZ) &
          + prod(-f(:,OP_DR)*h(:,OP_1),OP_DR,OP_DZZ) &
          + prod(-f(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_LP) &
          + prod(-f(:,OP_DR)*h(:,OP_1),OP_DR,OP_LP)
     if(itor.eq.1) then
        tempa = mu(f(:,OP_DZ),OP_DZ) + mu(f(:,OP_DR),OP_DR)
        temp = temp + prod(mu(ri_79*h(:,OP_1),OP_DZ),tempa)
     endif

  v3uun2 = temp
end function v3uun2


! V3uvn
! =====
function v3uvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3uvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3uvn = 0.
     return
  end if

  temp = intx5(e(:,:,OP_DZ),ri_79,f(:,OP_DRP),g(:,OP_1),h(:,OP_1)) &
       - intx5(e(:,:,OP_DR),ri_79,f(:,OP_DZP),g(:,OP_1),h(:,OP_1))

#else
  temp = 0.
#endif

  v3uvn = temp
end function v3uvn

function v3uvn1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3uvn1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3uvn1%len = 0
     return
  end if

     temp = prod( ri_79*g(:,OP_1)*h(:,OP_1),OP_DZ,OP_DRP) &
          + prod(-ri_79*g(:,OP_1)*h(:,OP_1),OP_DR,OP_DZP)
#else
  temp%len = 0
#endif

  v3uvn1 = temp
end function v3uvn1

function v3uvn2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3uvn2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3uvn2%len = 0
     return
  end if

     temp = prod( ri_79*f(:,OP_DRP)*h(:,OP_1),OP_DZ,OP_1) &
          + prod(-ri_79*f(:,OP_DZP)*h(:,OP_1),OP_DR,OP_1)
#else
  temp%len = 0
#endif

  v3uvn2 = temp
end function v3uvn2



! V3vvn
! =====
function v3vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vvn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3vvn = 0.
     return
  end if

  if(itor.eq.0) then
     temp = 0.
  else
     temp = -intx5(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_1),h(:,OP_1))
  endif

  v3vvn = temp
end function v3vvn

function v3vvn1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vvn1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3vvn1%len = 0
     return
  end if

  if(itor.eq.0) then
     temp%len = 0
  else
        temp = prod(-ri_79*g(:,OP_1)*h(:,OP_1),OP_DR,OP_1)
  endif

  v3vvn1 = temp
end function v3vvn1

function v3vvn2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vvn2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,  h
  type(prodarray) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3vvn2%len = 0
     return
  end if

  if(itor.eq.0) then
     temp%len = 0
  else
        temp = prod(-ri_79*f(:,OP_1)*h(:,OP_1),OP_DR,OP_1)
  endif

  v3vvn2 = temp
end function v3vvn2


! V3uchin
! =======
function v3uchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3uchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3uchin = 0.
     return
  end if


     temp79a = g(:,OP_DZZ)*f(:,OP_DR) - g(:,OP_DRZ)*f(:,OP_DZ) &
          +    g(:,OP_DZ)*f(:,OP_DRZ) - g(:,OP_DR)*f(:,OP_DZZ)
     temp79b = g(:,OP_DRZ)*f(:,OP_DR) - g(:,OP_DRR)*f(:,OP_DZ) &
          +    g(:,OP_DZ)*f(:,OP_DRR) - g(:,OP_DR)*f(:,OP_DRZ)

     temp = intx4(e(:,:,OP_DZ),ri3_79,temp79a,h(:,OP_1)) &
          + intx4(e(:,:,OP_DR),ri3_79,temp79b,h(:,OP_1)) &
          + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) &
          - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1))

     if(itor.eq.1) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
        temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

        temp = temp + &
             2.*intx4(e(:,:,OP_DZ),ri4_79,temp79a,h(:,OP_1)) &
             +  intx4(e(:,:,OP_DR),ri4_79,temp79b,h(:,OP_1))
     end if
  
  v3uchin = temp
end function v3uchin

function v3uchin1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3uchin1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(inertia.eq.0 .or. surface_int) then
     v3uchin1%len = 0
     return
  end if


     tempa = mu(g(:,OP_DZZ),OP_DR) + mu(-g(:,OP_DRZ),OP_DZ) &
          +  mu(g(:,OP_DZ),OP_DRZ) + mu(-g(:,OP_DR),OP_DZZ)
     tempb = mu(g(:,OP_DRZ),OP_DR) + mu(-g(:,OP_DRR),OP_DZ) &
          +  mu(g(:,OP_DZ),OP_DRR) + mu(-g(:,OP_DR),OP_DRZ)

     temp = prod( mu(ri3_79*h(:,OP_1),OP_DZ),tempa) &
          + prod( mu(ri3_79*h(:,OP_1),OP_DR),tempb) &
          + prod( ri3_79*g(:,OP_DR)*h(:,OP_1),OP_DZ,OP_GS) &
          + prod(-ri3_79*g(:,OP_DZ)*h(:,OP_1),OP_DR,OP_GS)

     if(itor.eq.1) then
        tempa = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
        tempb = mu(g(:,OP_DR),OP_DZ) + mu(-g(:,OP_DZ),OP_DR)

        temp = temp + &
                prod(mu(2.*ri4_79*h(:,OP_1),OP_DZ),tempa) &
             +  prod(mu(   ri4_79*h(:,OP_1),OP_DR),tempb)
     end if
  
  v3uchin1 = temp
end function v3uchin1

function v3uchin2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3uchin2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(inertia.eq.0 .or. surface_int) then
     v3uchin2%len = 0
     return
  end if


     tempa = mu(f(:,OP_DR),OP_DZZ) + mu(-f(:,OP_DZ),OP_DRZ) &
          +  mu(f(:,OP_DRZ),OP_DZ) + mu(-f(:,OP_DZZ),OP_DR)
     tempb = mu(f(:,OP_DR),OP_DRZ) + mu(-f(:,OP_DZ),OP_DRR) &
          +  mu(f(:,OP_DRR),OP_DZ) + mu(-f(:,OP_DRZ),OP_DR)

     temp = prod(mu(ri3_79*h(:,OP_1),OP_DZ),tempa) &
          + prod(mu(ri3_79*h(:,OP_1),OP_DR),tempb) &
          + prod( ri3_79*f(:,OP_GS)*h(:,OP_1),OP_DZ,OP_DR) &
          + prod(-ri3_79*f(:,OP_GS)*h(:,OP_1),OP_DR,OP_DZ)

     if(itor.eq.1) then
        tempa = mu(f(:,OP_DZ),OP_DZ) + mu(f(:,OP_DR),OP_DR)
        tempb = mu(f(:,OP_DZ),OP_DR) + mu(-f(:,OP_DR),OP_DZ)

        temp = temp + &
                prod(mu(2.*ri4_79*h(:,OP_1),OP_DZ),tempa) &
             +  prod(mu(   ri4_79*h(:,OP_1),OP_DR),tempb)
     end if
  
  v3uchin2 = temp
end function v3uchin2


! V3vchin
! =======
function v3vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3vchin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3vchin = 0.
     return
  end if

  temp = intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1)) &
       + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1))
#else
  temp = 0.
#endif
  v3vchin = temp
end function v3vchin

function v3vchin1(g,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vchin1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3vchin1%len = 0
     return
  end if

     temp = prod(ri4_79*g(:,OP_DZP)*h(:,OP_1),OP_DZ,OP_1) &
          + prod(ri4_79*g(:,OP_DRP)*h(:,OP_1),OP_DR,OP_1)
#else
  temp%len = 0
#endif
  v3vchin1 = temp
end function v3vchin1

function v3vchin2(f,h)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3vchin2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3vchin2%len = 0
     return
  end if

     temp = prod(ri4_79*f(:,OP_1)*h(:,OP_1),OP_DZ,OP_DZP) &
          + prod(ri4_79*f(:,OP_1)*h(:,OP_1),OP_DR,OP_DRP)
#else
  temp%len = 0
#endif
  v3vchin2 = temp
end function v3vchin2

! V3chichin
! =========
function v3chichin(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chichin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h
  vectype, dimension(dofs_per_element) :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3chichin = 0.
     return
  end if

     temp79a = f(:,OP_DZ)*g(:,OP_DZZ) + f(:,OP_DR)*g(:,OP_DRZ)
     temp79b = f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DR)*g(:,OP_DRR)

     temp = intx4(e(:,:,OP_DZ),ri6_79,temp79a,h(:,OP_1)) &
          + intx4(e(:,:,OP_DR),ri6_79,temp79b,h(:,OP_1))

     if(itor.eq.1) then
        temp = temp &
             - 2.*intx5(e(:,:,OP_DZ),ri7_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             - 2.*intx5(e(:,:,OP_DR),ri7_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
     endif

  v3chichin = temp
end function v3chichin

function v3chichin1(g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chichin1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g, h
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(inertia.eq.0 .or. surface_int) then
     v3chichin1%len = 0
     return
  end if

     tempa = mu(g(:,OP_DZZ),OP_DZ) + mu(g(:,OP_DRZ),OP_DR)
     tempb = mu(g(:,OP_DRZ),OP_DZ) + mu(g(:,OP_DRR),OP_DR)

     temp = prod(mu(ri6_79*h(:,OP_1),OP_DZ),tempa) &
          + prod(mu(ri6_79*h(:,OP_1),OP_DR),tempb)

     if(itor.eq.1) then
        temp = temp &
             + prod(-2.*ri7_79*g(:,OP_DR)*h(:,OP_1),OP_DZ,OP_DZ) &
             + prod(-2.*ri7_79*g(:,OP_DR)*h(:,OP_1),OP_DR,OP_DR)
     endif

  v3chichin1 = temp
end function v3chichin1

function v3chichin2(f,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chichin2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, h
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(inertia.eq.0 .or. surface_int) then
     v3chichin2%len = 0
     return
  end if

     tempa = mu(f(:,OP_DZ),OP_DZZ) + mu(f(:,OP_DR),OP_DRZ)
     tempb = mu(f(:,OP_DZ),OP_DRZ) + mu(f(:,OP_DR),OP_DRR)

     temp = prod(mu(ri6_79*h(:,OP_1),OP_DZ),tempa) &
          + prod(mu(ri6_79*h(:,OP_1),OP_DR),tempb)

     if(itor.eq.1) then
        temp = temp &
             + prod(-2.*ri7_79*f(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_DR) &
             + prod(-2.*ri7_79*f(:,OP_DR)*h(:,OP_1),OP_DR,OP_DR)
     endif

  v3chichin2 = temp
end function v3chichin2

! V3ngrav
! =======
function v3ngrav(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ngrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ngrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = gravz*intx2(e(:,:,OP_DZ),f(:,OP_1)) & 
          + gravr*intx3(e(:,:,OP_DR),ri2_79,f(:,OP_1)) 
  end if

  v3ngrav = temp
end function v3ngrav


! V3ungrav
! ========
function v3ungrav(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ungrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ungrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
     
     temp = gravz*intx3(e(:,:,OP_DZ), ri_79,temp79a) &
          + gravr*intx3(e(:,:,OP_DR),ri3_79,temp79a)
  end if

  v3ungrav = temp
end function v3ungrav


! V3chingrav
! ==========
function v3chingrav(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chingrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3chingrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = -(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
          + f(:,OP_LP)*g(:,OP_1))

     temp = gravz*intx2(e(:,:,OP_DZ),temp79a) &
          + gravr*intx3(e(:,:,OP_DR),ri2_79,temp79a)
  end if

  v3chingrav = temp
end function v3chingrav


! V3ndenmgrav
! ===========
function v3ndenmgrav(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3ndenmgrav
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  real, intent(in) :: g
  vectype, dimension(dofs_per_element) :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ndenmgrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = gravz*intx2(e(:,:,OP_DZ),f(:,OP_LP)) &
          + gravr*intx3(e(:,:,OP_DR),ri2_79,f(:,OP_LP))
  end if

  v3ndenmgrav = g*temp
end function v3ndenmgrav


! V3us
! ====
function v3us(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3us
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  if(surface_int) then
     temp = 0.
  else
     temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),temp79a) &
          - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),temp79a)
  end if

  v3us = temp
end function v3us

function v3us1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3us1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3us1%len = 0
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

  if(surface_int) then
     temp%len = 0
  else
     temp = prod( ri_79*temp79a,OP_DZ,OP_DR) &
          + prod(-ri_79*temp79a,OP_DR,OP_DZ)
  end if

  v3us1 = temp
end function v3us1

! V3chis
! ======
function v3chis(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3chis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

     if(surface_int) then
        temp = 0.
     else
        temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ),temp79a) &
             + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR),temp79a)
     end if

  v3chis = temp
end function v3chis

function v3chis1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: v3chis1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3chis1%len = 0
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) ! + denm*nt79(:,OP_LP)

     if(surface_int) then
        temp%len = 0
     else
        temp = prod(ri4_79*temp79a,OP_DZ,OP_DZ) &
             + prod(ri4_79*temp79a,OP_DR,OP_DR)
     end if

  v3chis1 = temp
end function v3chis1

! V3psiforce
! ===
vectype function v3psiforce(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp = -int4(ri3_79,g(:,OP_1),e(:,OP_DZ),f(:,OP_DR)) &
              + int4(ri3_79,g(:,OP_1),e(:,OP_DR),f(:,OP_DZ))
     end if


  v3psiforce = temp
  return
end function v3psiforce



! V3par
! =====
function v3par(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3par
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        if(itor.eq.1) then
           temp = - intx3(e(:,:,OP_DR),ri3_79,f(:,OP_1))
        endif
     end if

  v3par = temp
end function v3par

function v3parb2ipsipsi(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3parb2ipsipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79a = f(:,OP_1)*g(:,OP_1)*ri4_79
        temp =  temp                                                 &
             +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DR),i(:,OP_DRR))   &
             +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DZ),i(:,OP_DRZ))   &
             +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DRR),i(:,OP_DR))   &
             +  .5*intx4(e(:,:,OP_DR),temp79a,h(:,OP_DRZ),i(:,OP_DZ))   &
             +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DR),i(:,OP_DRZ))   &
             +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DZ),i(:,OP_DZZ))   &
             +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DRZ),i(:,OP_DR))   &
             +  .5*intx4(e(:,:,OP_DZ),temp79a,h(:,OP_DZZ),i(:,OP_DZ))   

        temp79b = f(:,OP_1)*ri4_79
        temp = temp                                                          &
             +  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DR)) &
             -  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ)) &
             -  intx5(e(:,:,OP_DZ),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DR)) &
             +  intx5(e(:,:,OP_DR),temp79b,g(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

        temp79c = g(:,OP_1)*ri4_79
        temp = temp                                                         &
             +  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DR)) &
             -  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DZ),h(:,OP_DR),i(:,OP_DZ)) &
             -  intx5(e(:,:,OP_DZ),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DR)) &
             +  intx5(e(:,:,OP_DR),temp79c,f(:,OP_DR),h(:,OP_DZ),i(:,OP_DZ))

        temp79d = -f(:,OP_1)*g(:,OP_1)*h(:,OP_GS)*ri4_79
        temp = temp                                   &
             +  intx3(e(:,:,OP_DR),temp79d,i(:,OP_DR)) &
             +  intx3(e(:,:,OP_DZ),temp79d,i(:,OP_DZ))   
     end if

  v3parb2ipsipsi = temp
end function v3parb2ipsipsi

function v3parb2ipsib(e,f,g,h,i)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3parb2ipsib
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g, h, i

  vectype,dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp79a = f(:,OP_DP)*g(:,OP_1)*i(:,OP_1)*ri5_79
     temp = intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DR)) &
          - intx3(e(:,:,OP_DR),temp79a,h(:,OP_DZ))
  end if
#endif

  v3parb2ipsib = temp
end function v3parb2ipsib

#ifdef USEPARTICLES
function v3p_2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3p_2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        temp = 0.
!!$        temp = &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR)) &
!!$             - intx4(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ))
     else
        temp = intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g) &
             + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g)
     end if

  v3p_2 = temp
end function v3p_2

function v3pbb(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3pbb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: g

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79b = g*ri4_79
        temp = temp                                                         &
             +  intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DR)) &
             -  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DZ),pst79(:,OP_DR),pst79(:,OP_DZ)) &
             -  intx5(e(:,:,OP_DZ),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DR)) &
             +  intx5(e(:,:,OP_DR),temp79b,f(:,OP_DR),pst79(:,OP_DZ),pst79(:,OP_DZ))

  #if defined(USE3D) || defined(USECOMPLEX)
     temp79a = f(:,OP_DP)*g*bzt79(:,OP_1)*ri5_79
     temp = temp + intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR)) &
          - intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ))
     temp79a = -g*ri3_79
        temp = temp +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
        temp79a = -g*ri3_79
        temp = temp +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DR)) &
              -intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DR),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              -intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DZ)) 
        temp79a = -f(:,OP_DP)*g*bzt79(:,OP_1)*ri4_79
        temp = temp +intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR)) &
              +intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ)) 
        temp79a = g*ri2_79
        temp = temp +intx5(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DR),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DZ),f(:,OP_DR),bfpt79(:,OP_DR)) &
              +intx5(e(:,:,OP_DR),temp79a,pst79(:,OP_DRP),f(:,OP_DZ),bfpt79(:,OP_DZ)) &
              +intx5(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZP),f(:,OP_DZ),bfpt79(:,OP_DZ)) 
      #endif
    end if

  v3pbb = temp
end function v3pbb

function v3jxb(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: v3jxb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS) :: f

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
     if(surface_int) then
        temp = 0.
     else
        temp79a = -f*pst79(:,OP_GS)*ri4_79
        temp = temp                                   &
             +  intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DR)) &
             +  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DZ))   
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = - f*bzt79(:,OP_1)*ri4_79
     temp = temp + intx3(e(:,:,OP_DR),temp79a,bzt79(:,OP_DR)+bfpt79(:,OP_DRP)) &
          - intx3(e(:,:,OP_DZ),temp79a,bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))
     temp79a = -f*bzt79(:,OP_1)*ri5_79
     temp = temp+intx3(e(:,:,OP_DR),temp79a,pst79(:,OP_DZP)) &
          -  intx3(e(:,:,OP_DZ),temp79a,pst79(:,OP_DRP))
     temp79a = -f*pst79(:,OP_GS)*ri3_79
     temp = temp+intx3(e(:,:,OP_DZ),temp79a,bfpt79(:,OP_DR)) &
          -  intx3(e(:,:,OP_DR),temp79a,bfpt79(:,OP_DZ))
     #endif
    end if

  v3jxb = temp
end function v3jxb

#endif
!==============================================================================
! B1 TERMS
!==============================================================================


! B1psi
! =====
!function b1psi(e,f)

  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: b1psi
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  !vectype, dimension(dofs_per_element) :: temp

!!!$  if(jadv.eq.0) then
!!!$     temp = int2(e(:,OP_1),f(:,OP_1))
!!!$  else
!!!$!     temp79a = e(:,OP_DZ)*f(:,OP_DZ) + e(:,OP_DR)*f(:,OP_DR)
!!!$!     temp = -int2(ri2_79,temp79a) 
!!!$!...changed 7/22/08    scj
!!!$     temp = -int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ)) &
!!!$            -int3(ri2_79,e(:,OP_DR),f(:,OP_DR))
!!!$  endif

  !if(surface_int) then
     !if(jadv.eq.0) then
        !temp = 0.
     !else
        !! this term doesn't seem to make a difference
!!        temp = 0.
        !temp = intx4(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR)) &
             !+ intx4(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ)) &
             !- intx4(e(:,:,OP_DR),f(:,OP_1),ri2_79,norm79(:,1)) &
             !- intx4(e(:,:,OP_DZ),f(:,OP_1),ri2_79,norm79(:,2))
     !end if
  !else
     !if(jadv.eq.0) then
        !temp = intx2(e(:,:,OP_1),f(:,OP_1))
     !else
        !temp = intx3(e(:,:,OP_GS),ri2_79,f(:,OP_1))
     !endif
  !end if

  !b1psi = temp
!end function b1psi

function b1psi

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1psi
  type(prodarray) :: temp

  if(surface_int) then
        ! this term doesn't seem to make a difference
        !temp%len = 0
        temp = prod(ri2_79*norm79(:,1),OP_1,OP_DR) &
             + prod(ri2_79*norm79(:,2),OP_1,OP_DZ) &
             + prod(-ri2_79*norm79(:,1),OP_DR,OP_1) &
             + prod(-ri2_79*norm79(:,2),OP_DZ,OP_1)
    else
        temp = prod(ri2_79,OP_GS,OP_1)
  end if

  b1psi = temp
end function b1psi


! B1psiu
! ======
function b1psiu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psiu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     ! surface terms
     if(surface_int) then
        if(jadv.eq.0) then
           temp = 0.
        else
!!$           temp79a = f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)
!!$           temp79b = f(:,OP_DRR)*g(:,OP_DZ ) - f(:,OP_DRZ)*g(:,OP_DR ) &
!!$                +    f(:,OP_DR )*g(:,OP_DRZ) - f(:,OP_DZ )*g(:,OP_DRR)
!!$           temp79c = f(:,OP_DRZ)*g(:,OP_DZ ) - f(:,OP_DZZ)*g(:,OP_DR ) &
!!$                +    f(:,OP_DR )*g(:,OP_DZZ) - f(:,OP_DZ )*g(:,OP_DRZ)
!!$
!!$           temp = int4(ri_79,e(:,OP_1 ),norm79(:,1),temp79b) &
!!$                + int4(ri_79,e(:,OP_1 ),norm79(:,2),temp79c) &
!!$                - int4(ri_79,e(:,OP_DR),norm79(:,1),temp79a) &
!!$                - int4(ri_79,e(:,OP_DZ),norm79(:,2),temp79a)
!!$           if(itor.eq.1) then
!!$              temp = temp + 2.*int4(ri2_79,e(:,OP_1),norm79(:,1),temp79a)
!!$           end if
           temp = 0.
        endif

     ! volume terms
     else
        if(jadv.eq.0) then
           temp = intx4(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_DZ)) &
                - intx4(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_DR))
        else
           temp = intx4(e(:,:,OP_GS),ri_79,f(:,OP_DR),g(:,OP_DZ)) &
                - intx4(e(:,:,OP_GS),ri_79,f(:,OP_DZ),g(:,OP_DR))
        endif
     end if
  
  b1psiu = temp
end function b1psiu

function b1psiu1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1psiu1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     ! surface terms
     if(surface_int) then
           temp%len = 0

     ! volume terms
     else
           temp = prod( ri_79*g(:,OP_DZ),OP_GS,OP_DR) &
                + prod(-ri_79*g(:,OP_DR),OP_GS,OP_DZ)
      end if
  
  b1psiu1 = temp
end function b1psiu1


function b1psiu2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1psiu2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

     ! surface terms
     if(surface_int) then
           temp%len = 0

     ! volume terms
     else
           temp = prod(ri_79*f(:,OP_DR),OP_GS,OP_DZ) &
                + prod(-ri_79*f(:,OP_DZ),OP_GS,OP_DR)
     end if
  b1psiu2 = temp
end function b1psiu2


! B1psiv
! ======
function b1psiv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psiv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
    temp = 0.
  else
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp = 0.
           else
           temp = &
                - intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRP),g(:,OP_1 ))&
                - intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZP),g(:,OP_1 ))&
                - intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DP))&
                - intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DP))
           endif
        else
           temp = intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1 )) &
                + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1 )) &
                + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR ),g(:,OP_DP)) &
                + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ ),g(:,OP_DP))
        endif
  endif
#else
  temp = 0
#endif

  b1psiv = temp
end function b1psiv

function b1psiv1(g)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1psiv1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp%len = 0
           else
           temp = &
                  prod(-ri2_79*norm79(:,1)*g(:,OP_1 ),OP_1,OP_DRP)&
                + prod(-ri2_79*norm79(:,2)*g(:,OP_1 ),OP_1,OP_DZP)&
                + prod(-ri2_79*norm79(:,1)*g(:,OP_DP),OP_1,OP_DR)&
                + prod(-ri2_79*norm79(:,2)*g(:,OP_DP),OP_1,OP_DZ)
           endif
        else
           temp = prod(ri2_79*g(:,OP_1 ),OP_DR,OP_DRP) &
                + prod(ri2_79*g(:,OP_1 ),OP_DZ,OP_DZP) &
                + prod(ri2_79*g(:,OP_DP),OP_DR,OP_DR) &
                + prod(ri2_79*g(:,OP_DP),OP_DZ,OP_DZ)
        endif
#else
  temp%len = 0
#endif

  b1psiv1 = temp
end function b1psiv1


function b1psiv2(f)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1psiv2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
    temp%len = 0.
  else
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp%len = 0
           else
           endif
        else
           temp = prod(ri2_79*f(:,OP_DRP),OP_DR,OP_1 ) &
                + prod(ri2_79*f(:,OP_DZP),OP_DZ,OP_1 ) &
                + prod(ri2_79*f(:,OP_DR ),OP_DR,OP_DP) &
                + prod(ri2_79*f(:,OP_DZ ),OP_DZ,OP_DP)
        endif
  endif
#else
  temp%len = 0
#endif

  b1psiv2 = temp
end function b1psiv2

! B1psid
! ======
!function b1psid(e,f,g)
  !use basic
  !use m3dc1_nint

  !vectype, dimension(dofs_per_element) :: b1psid
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

  !if(mass_ratio.eq.0. .or. db.eq.0.) then
     !b1psid = 0.
     !return
  !endif

  !if(surface_int) then
     !temp = 0.
  !else
     !temp = intx3(e(:,:,OP_1),f(:,OP_GS),g(:,OP_1))
  !endif

  !b1psid = temp*me_mp*mass_ratio*db**2
!end function b1psid

function b1psid(g)
  use basic
  use m3dc1_nint

  type(prodarray) :: b1psid
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(mass_ratio.eq.0. .or. db.eq.0.) then
     b1psid%len = 0
     return
  endif

  if(surface_int) then
     temp%len = 0
  else
     temp = prod(g(:,OP_1),OP_1,OP_GS)
  endif

  temp79a= me_mp*mass_ratio*db**2
  b1psid = temp*temp79a
end function b1psid

! B1bu
! ====
function b1bu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
        if(surface_int) then
           if(inoslip_pol.eq.1) then
              temp = 0.
           else
           temp = intx5(e(:,:,OP_1),ri2_79,norm79(:,1),g(:,OP_DRP),f(:,OP_1 ))&
                + intx5(e(:,:,OP_1),ri2_79,norm79(:,2),g(:,OP_DZP),f(:,OP_1 ))&
                + intx5(e(:,:,OP_1),ri2_79,norm79(:,1),g(:,OP_DR ),f(:,OP_DP))&
                + intx5(e(:,:,OP_1),ri2_79,norm79(:,2),g(:,OP_DZ ),f(:,OP_DP))
           endif
        else
           temp = -(intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZP)) &
                  + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DRP)) &
                  + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DP),g(:,OP_DZ)) &
                  + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DP),g(:,OP_DR)))
        end if
  endif

#else
  temp  = 0.
#endif
  b1bu = temp
end function b1bu

function b1bu1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1bu1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
        if(surface_int) then
           if(inoslip_pol.eq.1) then
              temp%len = 0
           endif
        else
           temp =  (prod(-ri2_79*g(:,OP_DZP),OP_DZ,OP_1) &
                  + prod(-ri2_79*g(:,OP_DRP),OP_DR,OP_1) &
                  + prod(-ri2_79*g(:,OP_DZ),OP_DZ,OP_DP) &
                  + prod(-ri2_79*g(:,OP_DR),OP_DR,OP_DP))
        end if

#else
  temp%len  = 0
#endif
  b1bu1 = temp
end function b1bu1

function b1bu2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1bu2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
        if(surface_int) then
              temp%len = 0
        else
           temp =  (prod(-ri2_79*f(:,OP_1),OP_DZ,OP_DZP) &
                  + prod(-ri2_79*f(:,OP_1),OP_DR,OP_DRP) &
                  + prod(-ri2_79*f(:,OP_DP),OP_DZ,OP_DZ) &
                  + prod(-ri2_79*f(:,OP_DP),OP_DR,OP_DR))
        end if
#endif
  b1bu2 = temp
end function b1bu2

! B1jrepsieta
! ====
function b1jrepsieta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jrepsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp = -(intx7(e(:,:,OP_1),norm79(:,2),ri3_79, &
                f(:,OP_1),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
               +intx7(e(:,:,OP_1),norm79(:,1),ri3_79, &
                f(:,OP_1),g(:,OP_DZP),h(:,OP_1),i(:,OP_1))) &
               -(intx7(e(:,:,OP_1),norm79(:,2),ri3_79,&
                f(:,OP_DP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
               +intx7(e(:,:,OP_1),norm79(:,1),ri3_79, &
                f(:,OP_DP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     else
        temp = -(intx6(e(:,:,OP_DZ),ri3_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
               -intx6(e(:,:,OP_DR),ri3_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)) &
               +intx6(e(:,:,OP_DZ),ri3_79,f(:,OP_DP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
               -intx6(e(:,:,OP_DR),ri3_79,f(:,OP_DP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     end if
     temp = temp*1.000
  endif

#else
  temp  = 0.
#endif
  b1jrepsieta = temp
end function b1jrepsieta

! B1jrebeta
! ====
function b1jrebeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jrebeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if (jadv .eq. 0) then
     temp = -intx5(e(:,:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_1),i(:,OP_1))
  else
     if(surface_int) then
        temp = 0.
     else
        temp = -intx6(e(:,:,OP_GS),ri2_79,f(:,OP_1),g(:,OP_1),h(:,OP_1),i(:,OP_1))
     end if
  endif
  temp = temp*1.000
  b1jrebeta = temp
end function b1jrebeta

! B1jrefeta
! ====
function b1jrefeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jrefeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp = (intx7(e(:,:,OP_1),norm79(:,2),ri2_79, &
                f(:,OP_DP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)) &
               +intx7(e(:,:,OP_1),norm79(:,1),ri2_79, &
                f(:,OP_DP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
               +intx7(e(:,:,OP_1),norm79(:,2),ri2_79, &
                f(:,OP_1),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)) &
               +intx7(e(:,:,OP_1),norm79(:,1),ri2_79, &
                f(:,OP_1),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)))
     else
        temp = (intx6(e(:,:,OP_DZ),ri2_79,f(:,OP_DP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1))&
               +intx6(e(:,:,OP_DR),ri2_79,f(:,OP_DP),g(:,OP_DR),h(:,OP_1),i(:,OP_1))&
               +intx6(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1),i(:,OP_1))&
               +intx6(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)))
     end if
     temp = temp*1.000
  endif

#else
  temp  = 0.
#endif
  b1jrefeta = temp
end function b1jrefeta

! B2jrepsieta
! ====
function b2jrepsieta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2jrepsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp = -(intx7(e(:,:,OP_1),norm79(:,2),ri2_79, &
                f(:,OP_1),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)) &
               +intx7(e(:,:,OP_1),norm79(:,1),ri2_79, &
                f(:,OP_1),g(:,OP_DR),h(:,OP_1),i(:,OP_1)))
     else
        temp = -(intx6(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)) &
                +intx6(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DR),h(:,OP_1),i(:,OP_1)))
     end if
     temp = temp*1.000
  endif

  b2jrepsieta = temp
end function b2jrepsieta

! B2jrefeta
! ====
function b2jrefeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2jrefeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp = -(intx7(e(:,:,OP_1),norm79(:,2),ri_79, &
                f(:,OP_1),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
               -intx7(e(:,:,OP_1),norm79(:,1),ri_79, &
                f(:,OP_1),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     else
        temp = -(intx6(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
               -intx6(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     end if
     temp = temp*1.000
  endif

#else
  temp  = 0.
#endif
  b2jrefeta = temp
end function b2jrefeta

! B1bv
! ====
function b1bv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  b1bv = 0.
end function b1bv

function b1bv1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1bv1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  b1bv1%len = 0
end function b1bv1

function b1bv2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1bv2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  b1bv2%len = 0
end function b1bv2

! B1psichi
! ========
function b1psichi(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b1psichi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = -intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ)) &
                  -intx4(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR))
        end if
     else
        if(surface_int) then
           temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
           temp = intx4(e(:,:,OP_DR),ri4_79,temp79a,norm79(:,1)) &
                + intx4(e(:,:,OP_DZ),ri4_79,temp79a,norm79(:,2)) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
                - intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
           if(itor.eq.1) then
              temp = temp + 2.*intx4(e(:,:,OP_1),ri5_79,temp79a,norm79(:,1))
           endif
        else
           temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZZ),g(:,OP_DZ )) &
                + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_DZZ)) &
                + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DRZ),g(:,OP_DR )) &
                + intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DR ),g(:,OP_DRZ)) &
                + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRZ),g(:,OP_DZ )) &
                + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DZ ),g(:,OP_DRZ)) &
                + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRR),g(:,OP_DR )) &
                + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_DRR))
           if(itor.eq.1) then
              temp = temp - 2.* &
                   (intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DR),g(:,OP_DR)) &
                   +intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DZ),g(:,OP_DZ)))
           endif
        end if
     endif

  b1psichi = temp
end function b1psichi

function b1psichi1(g)
  use basic
  use m3dc1_nint

  type(prodarray) :: b1psichi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp
  type(muarray) :: tempa

        if(surface_int) then
           tempa = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
           temp = prod(mu(ri4_79*norm79(:,1),OP_DR),tempa) &
                + prod(mu(ri4_79*norm79(:,2),OP_DZ),tempa) &
                + prod(-ri4_79*norm79(:,1)*g(:,OP_DZ),OP_1,OP_DRZ) &
                + prod(-ri4_79*norm79(:,1)*g(:,OP_DR),OP_1,OP_DRR) &
                + prod(-ri4_79*norm79(:,1)*g(:,OP_DRZ),OP_1,OP_DZ ) &
                + prod(-ri4_79*norm79(:,1)*g(:,OP_DRR),OP_1,OP_DR ) &
                + prod(-ri4_79*norm79(:,2)*g(:,OP_DZ),OP_1,OP_DZZ) &
                + prod(-ri4_79*norm79(:,2)*g(:,OP_DR),OP_1,OP_DRZ) &
                + prod(-ri4_79*norm79(:,2)*g(:,OP_DZZ),OP_1,OP_DZ ) &
                + prod(-ri4_79*norm79(:,2)*g(:,OP_DRZ),OP_1,OP_DR )
           if(itor.eq.1) then
              temp = temp + prod(mu(2.*ri5_79*norm79(:,1),OP_1),tempa)
           endif
        else
           temp = prod(ri4_79*g(:,OP_DZ ),OP_DZ,OP_DZZ) &
                + prod(ri4_79*g(:,OP_DZZ),OP_DZ,OP_DZ) &
                + prod(ri4_79*g(:,OP_DR ),OP_DZ,OP_DRZ) &
                + prod(ri4_79*g(:,OP_DRZ),OP_DZ,OP_DR) &
                + prod(ri4_79*g(:,OP_DZ ),OP_DR,OP_DRZ) &
                + prod(ri4_79*g(:,OP_DRZ),OP_DR,OP_DZ) &
                + prod(ri4_79*g(:,OP_DR ),OP_DR,OP_DRR) &
                + prod(ri4_79*g(:,OP_DRR),OP_DR,OP_DR)
            if(itor.eq.1) then
              temp = temp + &
                   (prod(-2.*ri5_79*g(:,OP_DR),OP_DR,OP_DR) &
                   +prod(-2.*ri5_79*g(:,OP_DZ),OP_DR,OP_DZ))
            endif
        end if

  b1psichi1 = temp
end function b1psichi1

function b1psichi2(f)
  use basic
  use m3dc1_nint

  type(prodarray) :: b1psichi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp
  type(muarray) :: tempa

        if(surface_int) then
           tempa = mu(f(:,OP_DZ),OP_DZ) + mu(f(:,OP_DR),OP_DR)
           temp = prod(mu(ri4_79*norm79(:,1),OP_DR),tempa) &
                + prod(mu(ri4_79*norm79(:,2),OP_DZ),tempa) &
                + prod(-ri4_79*norm79(:,1)*f(:,OP_DRZ),OP_1,OP_DZ ) &
                + prod(-ri4_79*norm79(:,1)*f(:,OP_DRR),OP_1,OP_DR ) &
                + prod(-ri4_79*norm79(:,1)*f(:,OP_DZ ),OP_1,OP_DRZ) &
                + prod(-ri4_79*norm79(:,1)*f(:,OP_DR ),OP_1,OP_DRR) &
                + prod(-ri4_79*norm79(:,2)*f(:,OP_DZZ),OP_1,OP_DZ ) &
                + prod(-ri4_79*norm79(:,2)*f(:,OP_DRZ),OP_1,OP_DR ) &
                + prod(-ri4_79*norm79(:,2)*f(:,OP_DZ ),OP_1,OP_DZZ) &
                + prod(-ri4_79*norm79(:,2)*f(:,OP_DR ),OP_1,OP_DRZ)
           if(itor.eq.1) then
              temp = temp + prod(mu(2.*ri5_79*norm79(:,1),OP_1),tempa)
           endif
        else
           temp = prod(ri4_79*f(:,OP_DZZ),OP_DZ,OP_DZ ) &
                + prod(ri4_79*f(:,OP_DZ ),OP_DZ,OP_DZZ) &
                + prod(ri4_79*f(:,OP_DRZ),OP_DZ,OP_DR ) &
                + prod(ri4_79*f(:,OP_DR ),OP_DZ,OP_DRZ) &
                + prod(ri4_79*f(:,OP_DRZ),OP_DR,OP_DZ ) &
                + prod(ri4_79*f(:,OP_DZ ),OP_DR,OP_DRZ) &
                + prod(ri4_79*f(:,OP_DRR),OP_DR,OP_DR ) &
                + prod(ri4_79*f(:,OP_DR ),OP_DR,OP_DRR)
            if(itor.eq.1) then
              temp = temp + &
                   (prod(-2.*ri5_79*f(:,OP_DR),OP_DR,OP_DR) &
                   +prod(-2.*ri5_79*f(:,OP_DZ),OP_DR,OP_DZ))
            endif
        end if

  b1psichi2= temp
end function b1psichi2


! B1bchi
! ======
function b1bchi(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b1bchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
        if(surface_int) then
           if(inoslip_pol.eq.1) then
              temp = 0.
           else
           temp = intx5(e(:,:,OP_1),ri5_79,f(:,OP_1 ),norm79(:,1),g(:,OP_DZP))&
                - intx5(e(:,:,OP_1),ri5_79,f(:,OP_1 ),norm79(:,2),g(:,OP_DRP))&
                + intx5(e(:,:,OP_1),ri5_79,f(:,OP_DP),norm79(:,1),g(:,OP_DZ ))&
                - intx5(e(:,:,OP_1),ri5_79,f(:,OP_DP),norm79(:,2),g(:,OP_DR ))
           endif
        else
           temp = intx4(e(:,:,OP_DZ),ri5_79,g(:,OP_DRP),f(:,OP_1))  &
                - intx4(e(:,:,OP_DR),ri5_79,g(:,OP_DZP),f(:,OP_1))  &
                + intx4(e(:,:,OP_DZ),ri5_79,g(:,OP_DR),f(:,OP_DP))  &
                - intx4(e(:,:,OP_DR),ri5_79,g(:,OP_DZ),f(:,OP_DP))
        end if
  endif

#else
  temp = 0.
#endif
  b1bchi = temp
end function b1bchi

function b1bchi1(g)
  use basic
  use m3dc1_nint

  type(prodarray) :: b1bchi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
        if(surface_int) then
              temp%len = 0
        else
           temp = prod( ri5_79*g(:,OP_DRP),OP_DZ,OP_1)  &
                + prod(-ri5_79*g(:,OP_DZP),OP_DR,OP_1)  &
                + prod( ri5_79*g(:,OP_DR),OP_DZ,OP_DP)  &
                + prod(-ri5_79*g(:,OP_DZ),OP_DR,OP_DP)
        end if

#else
  temp%len = 0
#endif
  b1bchi1 = temp
end function b1bchi1


function b1bchi2(f)
  use basic
  use m3dc1_nint

  type(prodarray) :: b1bchi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
        if(surface_int) then
           if(inoslip_pol.eq.1) then
              temp%len = 0
           else
           temp = prod( ri5_79*f(:,OP_1 )*norm79(:,1),OP_1,OP_DZP)&
                + prod(-ri5_79*f(:,OP_1 )*norm79(:,2),OP_1,OP_DRP)&
                + prod( ri5_79*f(:,OP_DP)*norm79(:,1),OP_1,OP_DZ )&
                + prod(-ri5_79*f(:,OP_DP)*norm79(:,2),OP_1,OP_DR )
            endif
        else
           temp = prod( ri5_79*f(:,OP_1),OP_DZ,OP_DRP)  &
                + prod(-ri5_79*f(:,OP_1),OP_DR,OP_DZP)  &
                + prod( ri5_79*f(:,OP_DP),OP_DZ,OP_DR)  &
                + prod(-ri5_79*f(:,OP_DP),OP_DR,OP_DZ)
         end if
#else
  temp%len = 0
#endif
  b1bchi2 = temp
end function b1bchi2


! B1psieta
! ========
!function b1psieta1(e,f,g,h,imod)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: b1psieta1
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  !vectype, dimension(dofs_per_element) :: temp
  !logical, intent(in) :: imod

  !if(jadv.eq.0) then

     !if(surface_int) then
        !temp = 0.
     !else
        !temp = intx3(e(:,:,OP_1),f(:,OP_GS),g(:,OP_1))
!#if defined(USE3D) || defined(USECOMPLEX)
        !if(iupstream.eq.1) then 
          !temp79a = abs(h(:,OP_1))*magus
          !temp = temp + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a)
        !elseif(iupstream.eq.2) then
          !temp79a = abs(h(:,OP_1))*magus
          !temp = temp - intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),temp79a)
        !endif
!#endif

         !if(hypf.ne.0) then
           !if(ihypeta.eq.1) then
              !temp = temp - hypf* &
                   !(intx3(e(:,:,OP_1),g(:,OP_LP),f(:,OP_GS)) &
                   !+intx3(e(:,:,OP_LP),g(:,OP_1),f(:,OP_GS)) &
                   !+2.*intx3(e(:,:,OP_DZ),g(:,OP_DZ),f(:,OP_GS)) &
                   !+2.*intx3(e(:,:,OP_DR),g(:,OP_DR),f(:,OP_GS)))
              !if(itor.eq.1) temp = temp  &
                  !- 2.*hypf*intx4(e(:,:,OP_1 ),ri_79 ,f(:,OP_GS),g(:,OP_DR)) &
                  !- 2.*hypf*intx4(e(:,:,OP_DR),ri_79 ,f(:,OP_GS),g(:,OP_1))
           !else
              !temp = temp - hypf*intx2(e(:,:,OP_LP),f(:,OP_GS))

              !if(itor.eq.1) then
                 !temp = temp - 2.*hypf*intx3(e(:,:,OP_DR),ri_79,f(:,OP_GS))
              !endif
           !end if
        !end if
     !end if
  !else
     !if(surface_int) then
        !if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           !temp = 0.
        !else
           !temp = 0.
        !endif

!#if defined(USE3D) || defined(USECOMPLEX)
        !if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           !temp = temp
        !else
           !if(.not.imod) then
              !temp = temp 
           !end if
        !endif
!#endif
     !else
        !temp = intx4(e(:,:,OP_GS),ri2_79,g(:,OP_1),f(:,OP_GS))

!#if defined(USE3D) || defined(USECOMPLEX)
        !if(.not.imod) then
           !temp = temp 
           !if(iupstream.eq.1) then   
              !temp79a = abs(h(:,OP_1))*magus
              !temp = temp - &
                !(intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZPP),temp79a) &
                !+intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRPP),temp79a))
              !if(itor.eq.1) then
                 !temp = temp + intx4(e(:,:,OP_DR),ri5_79,f(:,OP_DPP),temp79a)
              !endif
           !elseif(iupstream.eq.2) then
              !temp79a = abs(h(:,OP_1))*magus
              !temp = temp + &
                !(intx4(e(:,:,OP_DZPP),ri6_79,f(:,OP_DZPP),temp79a) &
                !+intx4(e(:,:,OP_DRPP),ri6_79,f(:,OP_DRPP),temp79a))
              !if(itor.eq.1) then
                 !temp = temp - intx4(e(:,:,OP_DRPP),ri7_79,f(:,OP_DPP),temp79a)
              !endif

           !endif
        !end if
!#endif
     !end if
  !endif

  !b1psieta1 = temp
!end function b1psieta1
function b1psieta1(g,h,imod)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1psieta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  logical, intent(in) :: imod

     if(surface_int) then
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp%len = 0
        else
           temp%len = 0
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp = temp
        else
           if(.not.imod) then
              temp = temp
           end if
        endif
#endif
     else
        temp = prod(ri2_79*g(:,OP_1),OP_GS,OP_GS)

#if defined(USE3D) || defined(USECOMPLEX)
        if(.not.imod) then
           temp = temp
        end if
#endif
     end if

  b1psieta1 = temp
end function b1psieta1

!function b1psieta2(e,f,g,h,imod)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: b1psieta2
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  !vectype, dimension(dofs_per_element) :: temp
  !logical, intent(in) :: imod

  !if(jadv.eq.0) then

     !if(surface_int) then
        !temp = 0.
     !else
        !temp = 0.
     !endif

  !else
     !if(surface_int) then
        !temp = 0

!#if defined(USE3D) || defined(USECOMPLEX)
        !if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           !temp = temp
        !else
           !if(.not.imod) then
              !temp = temp &
                   !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,1),f(:,OP_DRPP),g(:,OP_1)) &
                   !+ intx5(e(:,:,OP_1),ri4_79,norm79(:,2),f(:,OP_DZPP),g(:,OP_1))
           !end if
        !endif
!#endif
     !else

!#if defined(USECOMPLEX)
        !if(.not.imod) then
           !temp = - &
                !(intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZPP),g(:,OP_1)) &
                !+intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRPP),g(:,OP_1)) &
                !+intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_DP)) &
                !+intx4(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_DP)))
        !end if
!#endif
!#if defined(USE3D)
        !if(.not.imod) then
           !temp =  &
                 !intx4(e(:,:,OP_DZP),ri4_79,f(:,OP_DZP),g(:,OP_1)) &
                !+intx4(e(:,:,OP_DRP),ri4_79,f(:,OP_DRP),g(:,OP_1)) 
        !end if
!#endif

     !end if
  !endif

  !b1psieta2 = temp
!end function b1psieta2

function b1psieta2(g,h,imod)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1psieta2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp
  logical, intent(in) :: imod

  if(jadv.eq.0) then

     if(surface_int) then
        temp%len = 0
     else
        temp%len = 0
     endif

  else
     if(surface_int) then
        temp%len = 0

#if defined(USE3D) || defined(USECOMPLEX)
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp = temp
        else
           if(.not.imod) then
              temp = temp &
                   + prod(ri4_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DRPP) &
                   + prod(ri4_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DZPP)
           end if
        endif
#endif
     else

#if defined(USECOMPLEX)
        if(.not.imod) then
           temp = &
                (prod(-ri4_79*g(:,OP_1),OP_DZ,OP_DZPP) &
                +prod(-ri4_79*g(:,OP_1),OP_DR,OP_DRPP) &
                +prod(-ri4_79*g(:,OP_DP),OP_DZ,OP_DZP) &
                +prod(-ri4_79*g(:,OP_DP),OP_DR,OP_DRP))
        end if
#endif
#if defined(USE3D)
        if(.not.imod) then
           temp =  &
                 prod(ri4_79*g(:,OP_1),OP_DZP,OP_DZP) &
                +prod(ri4_79*g(:,OP_1),OP_DRP,OP_DRP)
        end if
#endif

     end if
  endif

  b1psieta2 = temp
end function b1psieta2

! B1jeta
! ========
function b1jeta(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1jeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp = 0.
        else
           temp = 0.
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp = temp
        else
!           temp = temp &
!                + int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DRPP),g(:,OP_1)) &
!                + int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZPP),g(:,OP_1))
           temp = temp
        endif
#endif
     else

           ! note: f = -delstar(psi)

           if(ihypeta.eq.1) then

              temp = hypf*intx4(e(:,:,OP_GS),ri2_79,f(:,OP_GS),g(:,OP_1))

#if defined(USE3D) || defined(USECOMPLEX)
              temp = temp + hypf* &
                  (intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DP) ,g(:,OP_DP)) &
                  +intx4(e(:,:,OP_GS),ri4_79,f(:,OP_DPP),g(:,OP_1)))
#endif

           else
              temp = hypf*intx3(e(:,:,OP_GS),ri2_79,f(:,OP_GS))

#if defined(USE3D) || defined(USECOMPLEX)
              temp = temp + hypf* &
                   intx3(e(:,:,OP_GS),ri4_79,f(:,OP_DPP))
#endif

           endif
     end if

  b1jeta = temp
end function b1jeta


! B1beta
! ======
function b1beta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1beta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else 
     if(surface_int) then
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),ri3_79,g(:,OP_1),norm79(:,2),f(:,OP_DRP)) &
                - intx5(e(:,:,OP_1),ri3_79,g(:,OP_1),norm79(:,1),f(:,OP_DZP))
        endif
     else
#if defined(USE3D)
        temp = - (intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZ),g(:,OP_1 )) &
                - intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DR),g(:,OP_1 )))
#endif
#if defined(USECOMPLEX)
        temp = intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_1 )) &
             - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_1 )) &
             + intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZ ),g(:,OP_DP)) &
             - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DR ),g(:,OP_DP))
#endif
#ifndef USEST
        if(hypf.gt.0. .and. imp_hyper.le.1) then
           if(ihypeta.eq.0) then
              temp = temp - hypf*intx3(e(:,:,OP_DZP),ri5_79,f(:,OP_DRPP)) &
                          + hypf*intx3(e(:,:,OP_DRP),ri5_79,f(:,OP_DZPP))
              if(itor.eq.1) then
                 temp = temp + 4.*hypf*intx3(e(:,:,OP_DZP),ri6_79,f(:,OP_DPP))
              endif
           else
              temp = temp + hypf* &
                   (-intx4(e(:,:,OP_DZP),ri5_79,g(:,OP_1),f(:,OP_DRPP))   &
                   + intx4(e(:,:,OP_DRP),ri5_79,g(:,OP_1),f(:,OP_DZPP))   &
                   - intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DRP),f(:,OP_GS))   &
                   + intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DZP),f(:,OP_GS))   &
                   - intx4(e(:,:,OP_DZ),ri3_79,g(:,OP_DR),f(:,OP_GSP))   &
                   + intx4(e(:,:,OP_DR),ri3_79,g(:,OP_DZ),f(:,OP_GSP)))
              if(itor.eq.1) then
                 temp = temp + &
                      4.*hypf*intx4(e(:,:,OP_DZP),ri6_79,g(:,OP_1),f(:,OP_DPP))
              endif

           endif
        endif
#endif !USEST
     endif
  endif
#else
  temp = 0.
#endif

  b1beta = temp
end function b1beta

function b1beta1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1beta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp%len = 0
  else 
     if(surface_int) then
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp%len = 0.
        else
           temp = prod( ri3_79*g(:,OP_1)*norm79(:,2),OP_1,OP_DRP) &
                + prod(-ri3_79*g(:,OP_1)*norm79(:,1),OP_1,OP_DZP)
        endif
     else
#if defined(USE3D)
        temp =  (prod(-ri3_79*g(:,OP_1 ),OP_DRP,OP_DZ) &
                + prod(ri3_79*g(:,OP_1 ),OP_DZP,OP_DR))
#endif
#if defined(USECOMPLEX)
        temp = prod( ri3_79*g(:,OP_1 ),OP_DR,OP_DZP) &
             + prod(-ri3_79*g(:,OP_1 ),OP_DZ,OP_DRP) &
             + prod( ri3_79*g(:,OP_DP),OP_DR,OP_DZ) &
             + prod(-ri3_79*g(:,OP_DP),OP_DZ,OP_DR)
#endif
#ifndef USEST
        if(hypf.gt.0 .and. imp_hyper.le.1) then
           if(ihypeta.eq.0) then
              temp = temp + prod(-hypf*ri5_79,OP_DZP,OP_DRPP) &
                          + prod( hypf*ri5_79,OP_DRP,OP_DZPP)
              if(itor.eq.1) then
                 temp = temp + prod(4.*hypf*ri6_79,OP_DZP,OP_DPP)
              endif
           else

           endif
        endif
#endif !USEST
     endif
  endif
#else
  temp%len = 0
#endif

  b1beta1 = temp
end function b1beta1

! B1psij
! ======
function b1psij(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psij
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = hypf*(g(:,OP_DRR)+ri2_79*g(:,OP_DPP) +g(:,OP_DZZ))
     if(itor.eq.1) temp79a = temp79a + hypf*ri_79*g(:,OP_DR)
     if     (ihypeta.eq.1) then
        temp79a = eta79(:,OP_1)*temp79a          &
          + hypf*(eta79(:,OP_DR)*g(:,OP_DR)      &
         + ri2_79*eta79(:,OP_DP)*g(:,OP_DP)      &
                + eta79(:,OP_DZ)*g(:,OP_DZ))
     else if(ihypeta.eq.2) then
        temp79a = pt79(:,OP_1)*temp79a          &
          + hypf*(pt79(:,OP_DR)*g(:,OP_DR)      &
         + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
                + pt79(:,OP_DZ)*g(:,OP_DZ))
     else if(ihypeta.gt.2) then
        temp79a = (bharhypeta)**beta*pt79(:,OP_1)*temp79a          &
          + hypf*(bharhypeta)**beta*(pt79(:,OP_DR)*g(:,OP_DR)      &
         + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
                + pt79(:,OP_DZ)*g(:,OP_DZ))
     endif
     temp = -intx5(e(:,:,OP_DZP),ri3_79,b2i79(:,OP_1),f(:,OP_DR),temp79a)    &
            +intx5(e(:,:,OP_DRP),ri3_79,b2i79(:,OP_1),f(:,OP_DZ),temp79a)
     
#else
  temp = 0.
#endif

  b1psij = temp
end function b1psij

! B1bj
! ====
function b1bj(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bj
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp
  real :: hypfm


     temp79a = hypf*(g(:,OP_DRR)+g(:,OP_DZZ))
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a + hypf*ri2_79*g(:,OP_DPP)
#endif
     if(itor.eq.1) temp79a = temp79a + hypf*ri_79*g(:,OP_DR)
     if     (ihypeta.eq.1) then             
         temp79a = eta79(:,OP_1)*temp79a          &
          + hypf*(eta79(:,OP_DR)*g(:,OP_DR)      &
                 + eta79(:,OP_DZ)*g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp79a = temp79a + hypf*ri2_79*eta79(:,OP_DP)*g(:,OP_DP)
#endif
     else if(ihypeta.eq.2) then             
         temp79a = pt79(:,OP_1)*temp79a          &
          + hypf*(pt79(:,OP_DR)*g(:,OP_DR)      &
                 + pt79(:,OP_DZ)*g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp79a = temp79a + hypf*ri2_79*pt79(:,OP_DP)*g(:,OP_DP)
#endif
     else if(ihypeta.gt.2) then
         hypfm = hypf*(bharhypeta)**beta
         temp79a = pt79(:,OP_1)*temp79a*(bharhypeta)**beta &
          + hypfm*(pt79(:,OP_DR)*g(:,OP_DR)      &
                 + pt79(:,OP_DZ)*g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp79a = temp79a + hypfm*ri2_79*pt79(:,OP_DP)*g(:,OP_DP)
#endif
     endif

     temp = intx5(e(:,:,OP_GS),ri2_79,b2i79(:,OP_1),f(:,OP_1),temp79a)

  b1bj = temp
end function b1bj

! B1fj
! ======
function b1fj(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fj
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = hypf*(g(:,OP_DRR)+ri2_79*g(:,OP_DPP) +g(:,OP_DZZ))
     if(itor.eq.1) temp79a = temp79a + hypf*ri_79*g(:,OP_DR)
     if     (ihypeta.eq.1) then
        temp79a = eta79(:,OP_1)*temp79a          &
                + hypf*(eta79(:,OP_DR)*g(:,OP_DR)      &
               + ri2_79*eta79(:,OP_DP)*g(:,OP_DP)      &
                      + eta79(:,OP_DZ)*g(:,OP_DZ))
     else if(ihypeta.eq.2) then
        temp79a = pt79(:,OP_1)*temp79a          &
                + hypf*(pt79(:,OP_DR)*g(:,OP_DR)      &
               + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
                      + pt79(:,OP_DZ)*g(:,OP_DZ))
     else if(ihypeta.gt.2) then
        temp79a = pt79(:,OP_1)*temp79a*(bharhypeta)**beta  &
                + hypf*(bharhypeta)**beta*(pt79(:,OP_DR)*g(:,OP_DR)      &
                                   + ri2_79*pt79(:,OP_DP)*g(:,OP_DP)      &
                                          + pt79(:,OP_DZ)*g(:,OP_DZ))
     endif
     temp =  intx5(e(:,:,OP_DRP),ri2_79,b2i79(:,OP_1),f(:,OP_DR),temp79a)    &
            +intx5(e(:,:,OP_DZP),ri2_79,b2i79(:,OP_1),f(:,OP_DZ),temp79a)
#else
  temp = 0.
#endif

  b1fj = temp
end function b1fj



! B1psipsid
! =========
function b1psipsid(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psipsid
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psipsid = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
             + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1))
     endif
  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZP) + f(:,OP_DR)*g(:,OP_DRP)
        temp79b = ri4_79*h(:,OP_1)
        temp79c = ri4_79*h(:,OP_DP)
        temp = intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             + intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DRP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRRP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DRP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZP))&
             - intx5(e(:,:,OP_DR),ri4_79,temp79a,norm79(:,1),h(:,OP_1)) &
             - intx5(e(:,:,OP_DZ),ri4_79,temp79a,norm79(:,2),h(:,OP_1)) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GSP),norm79(:,1),g(:,OP_DR )) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GSP),norm79(:,2),g(:,OP_DZ )) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GS ),norm79(:,1),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GS ),norm79(:,2),g(:,OP_DZP)) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_GS ),norm79(:,1),g(:,OP_DR )) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_GS ),norm79(:,2),g(:,OP_DZ ))
        if(itor.eq.1) then
           temp = temp &
                - 2.*intx5(e(:,:,OP_1),ri5_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = intx5(e(:,:,OP_GS),ri4_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
             + intx5(e(:,:,OP_GS),ri4_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
#ifdef USEST
             - intx5(e(:,:,OP_DZP),ri4_79,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_1)) &
             - intx5(e(:,:,OP_DRP),ri4_79,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_1))
#else
             + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_GS ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_GS ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_GSP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_GSP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DP)) &
             + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_DP))
#endif
     end if
  endif
  b1psipsid = temp
#else
  b1psipsid = 0.
#endif
end function b1psipsid


! B1psibd1
! ========
function b1psibd1(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibd1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibd1 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             - intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
     endif

  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
        temp79b = ri3_79*h(:,OP_1)
        temp = &
             + intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             + intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ)) &
             - intx5(e(:,:,OP_DR),ri3_79,temp79a,norm79(:,1),h(:,OP_1)) &
             - intx5(e(:,:,OP_DZ),ri3_79,temp79a,norm79(:,2),h(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                - intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             - intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
     end if
  endif

  b1psibd1 = temp
end function b1psibd1

! B1psibd2
! ========
function b1psibd2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibd2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibd2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp79b = ri5_79*h(:,OP_1)
        temp79c = ri5_79*h(:,OP_DP)
        temp = &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRPP),g(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZPP),g(:,OP_1 )) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRP ),g(:,OP_DP)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZP ),g(:,OP_DP)) &
             + intx5(e(:,:,OP_1),temp79c,norm79(:,2),f(:,OP_DRP ),g(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79c,norm79(:,1),f(:,OP_DZP ),g(:,OP_1 ))
     else
        temp = &
#ifdef USEST
             - intx5(e(:,:,OP_DRP),ri5_79,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_1)) &
             + intx5(e(:,:,OP_DZP),ri5_79,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_1))
#else
             + intx5(e(:,:,OP_DR),ri5_79,f(:,OP_DZPP),g(:,OP_1 ),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DZ),ri5_79,f(:,OP_DRPP),g(:,OP_1 ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),ri5_79,f(:,OP_DZP ),g(:,OP_DP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP ),g(:,OP_DP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),ri5_79,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_DP)) &
             - intx5(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_DP))
#endif
     end if
  endif
#else
  temp = 0.
#endif

  b1psibd2 = temp
end function b1psibd2



! B1psifd1
! ========
function b1psifd1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifd1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifd1 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR ),h(:,OP_1)) &
             - intx5(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ ),h(:,OP_1))
     end if
  else
     if(surface_int) then
        temp79b = ri3_79*h(:,OP_1)
        temp79c = ri3_79*h(:,OP_DP)
        temp = &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DRP),f(:,OP_GS ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DZP),f(:,OP_GS ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DR ),f(:,OP_GSP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DZ ),f(:,OP_GSP))&
             - intx5(e(:,:,OP_1),temp79c,norm79(:,2),g(:,OP_DR ),f(:,OP_GS ))&
             + intx5(e(:,:,OP_1),temp79c,norm79(:,1),g(:,OP_DZ ),f(:,OP_GS ))
     else
        temp = intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DZP),g(:,OP_DR ),h(:,OP_1)) &
             - intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DRP),g(:,OP_DZ ),h(:,OP_1)) &
#ifdef USEST
             - intx5(e(:,:,OP_DZP),ri3_79,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_1)) &
             + intx5(e(:,:,OP_DRP),ri3_79,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_1))
#else
             + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GSP),g(:,OP_DR ),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GSP),g(:,OP_DZ ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS ),g(:,OP_DRP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS ),g(:,OP_DZP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_DP)) &
             - intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_DP))
#endif
     endif
  endif
  b1psifd1 = temp
#else
  b1psifd1 = 0.
#endif
end function b1psifd1


! B1psifd2
! ========
function b1psifd2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifd2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifd2 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ ),g(:,OP_DRP),h(:,OP_1)) &
             - intx5(e(:,:,OP_1),ri_79,f(:,OP_DR ),g(:,OP_DZP),h(:,OP_1))
     end if
  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
        temp79b = ri3_79*h(:,OP_1)
        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             + intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ))&
             - intx5(e(:,:,OP_DR),ri3_79,temp79a,norm79(:,1),h(:,OP_1))     &
             - intx5(e(:,:,OP_DZ),ri3_79,temp79a,norm79(:,2),h(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                + intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DZ ),g(:,OP_DRP),h(:,OP_1)) &
             - intx5(e(:,:,OP_GS),ri3_79,f(:,OP_DR ),g(:,OP_DZP),h(:,OP_1))
     endif
  endif
  b1psifd2 = temp
#else
  b1psifd2 = 0.
#endif
end function b1psifd2




! B1bbd
! =====
function b1bbd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bbd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bbd = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if
  else
     if(surface_int) then
        temp79a = ri4_79*f(:,OP_DP)
        temp79b = ri4_79*f(:,OP_1)
        temp = &
             - intx5(e(:,:,OP_1),temp79a,norm79(:,1),g(:,OP_DR ),h(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79a,norm79(:,2),g(:,OP_DZ ),h(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DRP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DZP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DR ),h(:,OP_DP)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DZ ),h(:,OP_DP))
     else
        temp = intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_1 ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_1 ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_DP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_DP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ ),g(:,OP_1 ),h(:,OP_DP)) &
             + intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DR ),g(:,OP_1 ),h(:,OP_DP))
     endif
  endif
  b1bbd = temp
#else
  b1bbd = 0.
#endif
end function b1bbd


! B1bfd1
! ======
function b1bfd1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfd1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfd1 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = intx4(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             + intx4(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
     end if
  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
        temp79b = ri4_79*h(:,OP_1)
        temp = intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             + intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))&
             - intx5(e(:,:,OP_DR),ri4_79,temp79a,norm79(:,1),h(:,OP_1)) &
             - intx5(e(:,:,OP_DZ),ri4_79,temp79a,norm79(:,2),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp &
                - 2.*intx5(e(:,:,OP_1),ri5_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = intx5(e(:,:,OP_GS),ri2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             + intx5(e(:,:,OP_GS),ri2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
     endif
  endif
  b1bfd1 = temp
#else
  b1bfd1 = 0.
#endif
end function b1bfd1

! B1bfd2
! ======
function b1bfd2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfd2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfd2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp79b = ri4_79*h(:,OP_1)
        temp79c = ri4_79*h(:,OP_DP)
        temp = &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_DP),norm79(:,1),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_DP),norm79(:,2),g(:,OP_DZP)) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_1 ),norm79(:,1),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_1 ),norm79(:,2),g(:,OP_DZP))

     else
#ifdef USECOMPLEX
        temp = intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_DP),h(:,OP_1 ))  &
             + intx5(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_DP),h(:,OP_1 ))  &
             + intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_1 ),h(:,OP_DP))  &
             + intx5(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_1 ),h(:,OP_DP))

        ! f''' term hack
        temp = temp + rfac* &
             (intx5(e(:,:,OP_DZ),ri4_79,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
             +intx5(e(:,:,OP_DR),ri4_79,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#elif defined(USE3D)
        ! here, we can integrate by parts
        temp = - &
             (intx5(e(:,:,OP_DZP),ri4_79,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
             +intx5(e(:,:,OP_DRP),ri4_79,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#endif
     endif
  endif
  b1bfd2 = temp
#else
  b1bfd2 = 0.
#endif
end function b1bfd2



! B1ped
! =====
function b1ped(e,f,g)
  use basic
  use m3dc1_nint

  implicit none
  vectype, dimension(dofs_per_element) :: b1ped
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1ped = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = intx3(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1))
     end if
  else
     if(surface_int) then
        temp = intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),norm79(:,1),g(:,OP_DR)) &
             + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),norm79(:,2),g(:,OP_DZ)) &
             - intx5(e(:,:,OP_1),ri2_79,g(:,OP_DP),norm79(:,1),f(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri2_79,g(:,OP_DP),norm79(:,2),f(:,OP_DZ)) &
             - intx5(e(:,:,OP_DR),ri2_79,g(:,OP_1),f(:,OP_DP),norm79(:,1)) &
             - intx5(e(:,:,OP_DZ),ri2_79,g(:,OP_1),f(:,OP_DP),norm79(:,2))
     else
        temp = intx4(e(:,:,OP_GS),ri2_79,f(:,OP_DP),g(:,OP_1)) &
             + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1 )) &
             + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1 )) &
             + intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ ),g(:,OP_DP)) &
             + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR ),g(:,OP_DP))
     end if
  endif
  b1ped = temp
#else
  b1ped = 0.
#endif
end function b1ped

! B1psipsin
! =========
function b1psipsin(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psipsin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psipsin = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp79a = ri2_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),g(:,OP_DZP)) &
             + intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),g(:,OP_DRP))
     endif
  else
     if(surface_int) then
        temp79a = (f(:,OP_DZ)*g(:,OP_DZP) + f(:,OP_DR)*g(:,OP_DRP)) &
             *ni79(:,OP_1)**2
        temp79b = ri4_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp79c = -ri4_79*h(:,OP_DP)*ni79(:,OP_1)**2
        temp =-intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DRP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRRP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DRP ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZP))&
             - intx5(e(:,:,OP_DR),ri4_79,temp79a,norm79(:,1),h(:,OP_1)) &
             - intx5(e(:,:,OP_DZ),ri4_79,temp79a,norm79(:,2),h(:,OP_1)) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GSP),norm79(:,1),g(:,OP_DR )) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GSP),norm79(:,2),g(:,OP_DZ )) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GS ),norm79(:,1),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_GS ),norm79(:,2),g(:,OP_DZP)) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_GS ),norm79(:,1),g(:,OP_DR )) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_GS ),norm79(:,2),g(:,OP_DZ ))
        if(itor.eq.1) then
           temp = temp &
                - 2.*intx5(e(:,:,OP_1),ri5_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp79a = ri4_79*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
             + intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
#ifdef USEST
             - intx5(e(:,:,OP_DZP),temp79a,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_1)) &
             - intx5(e(:,:,OP_DRP),temp79a,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_1)) &
             - 2*intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DP)) &
             - 2*intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_DP))
#else
             + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_GS ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_GS ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_GSP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_GSP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DP)) &
             - intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_GS ),h(:,OP_DP))
#endif
     end if
  endif
  b1psipsin = temp
#else
  b1psipsin = 0.
#endif
end function b1psipsin



! B1psibn1
! ========
function b1psibn1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibn1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibn1 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp79a = ri_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),g(:,OP_DR)) &
             - intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),g(:,OP_DZ))
     endif

  else
     if(surface_int) then
        temp79a = (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) &
             * ni79(:,OP_1)**2
        temp79b = ri3_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp = &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ)) &
             - intx5(e(:,:,OP_DR),ri3_79,temp79a,norm79(:,1),h(:,OP_1)) &
             - intx5(e(:,:,OP_DZ),ri3_79,temp79a,norm79(:,2),h(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                - intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp79a = ri3_79*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             - intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
     end if
  endif

  b1psibn1 = temp
end function b1psibn1

! B1psibd2
! ========
function b1psibn2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psibn2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psibn2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp79b = ri5_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp79c = -ri5_79*h(:,OP_DP)*ni79(:,OP_1)**2
        temp = &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRPP),g(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZPP),g(:,OP_1 )) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRP ),g(:,OP_DP)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZP ),g(:,OP_DP)) &
             + intx5(e(:,:,OP_1),temp79c,norm79(:,2),f(:,OP_DRP ),g(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79c,norm79(:,1),f(:,OP_DZP ),g(:,OP_1 ))
     else
        temp79a = ri5_79*ni79(:,OP_1)**2
        temp = &
#ifdef USEST
             - intx5(e(:,:,OP_DRP),temp79a,f(:,OP_DZP ),g(:,OP_1),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZP),temp79a,f(:,OP_DRP ),g(:,OP_1),h(:,OP_1 )) &
             - 2*intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_DP)) &
             + 2*intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_DP))
#else
             + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZPP),g(:,OP_1 ),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRPP),g(:,OP_1 ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZP ),g(:,OP_DP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRP ),g(:,OP_DP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_DP)) &
             + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_DP))
#endif
     end if
  endif
  b1psibn2 = temp
#else
  b1psibn2 = 0.
#endif
end function b1psibn2


! B1psifn1
! ========
function b1psifn1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifn1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifn1 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp79a =ri_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZP),g(:,OP_DR )) &
             - intx4(e(:,:,OP_1),temp79a,f(:,OP_DRP),g(:,OP_DZ ))
     end if
  else
     if(surface_int) then
        temp79b = ri3_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp79c = -ri3_79*h(:,OP_DP)*ni79(:,OP_1)**2
        temp = &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DRP),f(:,OP_GS ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DZP),f(:,OP_GS ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DR ),f(:,OP_GSP))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DZ ),f(:,OP_GSP))&
             - intx5(e(:,:,OP_1),temp79c,norm79(:,2),g(:,OP_DR ),f(:,OP_GS ))&
             + intx5(e(:,:,OP_1),temp79c,norm79(:,1),g(:,OP_DZ ),f(:,OP_GS ))
     else
        temp79a = ri3_79*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZP),g(:,OP_DR ),h(:,OP_1)) &
             - intx5(e(:,:,OP_GS),temp79a,f(:,OP_DRP),g(:,OP_DZ ),h(:,OP_1)) &
#ifdef USEST
             - intx5(e(:,:,OP_DZP),temp79a,f(:,OP_GS ),g(:,OP_DR),h(:,OP_1 ))&
             + intx5(e(:,:,OP_DRP),temp79a,f(:,OP_GS ),g(:,OP_DZ),h(:,OP_1 ))&
             - 2*intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_DP))&
             + 2*intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_DP))
#else
             + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GSP),g(:,OP_DR ),h(:,OP_1 ))&
             - intx5(e(:,:,OP_DR),temp79a,f(:,OP_GSP),g(:,OP_DZ ),h(:,OP_1 ))&
             + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS ),g(:,OP_DRP),h(:,OP_1 ))&
             - intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS ),g(:,OP_DZP),h(:,OP_1 ))&
             - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS ),g(:,OP_DR ),h(:,OP_DP))&
             + intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS ),g(:,OP_DZ ),h(:,OP_DP))
#endif
     endif
  endif
  b1psifn1 = temp
#else
  b1psifn1 = 0.
#endif
end function b1psifn1


! B1psifn2
! ========
function b1psifn2(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psifn2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1psifn2 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp79a = ri_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ ),g(:,OP_DRP)) &
             - intx4(e(:,:,OP_1),temp79a,f(:,OP_DR ),g(:,OP_DZP))
     end if
  else
     if(surface_int) then
        temp79a = (f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)) &
             *ni79(:,OP_1)**2
        temp79b = ri3_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp =-intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ ))&
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ))&
             - intx5(e(:,:,OP_DR),ri3_79,temp79a,norm79(:,1),h(:,OP_1))     &
             - intx5(e(:,:,OP_DZ),ri3_79,temp79a,norm79(:,2),h(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                + intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp79a = ri3_79*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ ),g(:,OP_DRP),h(:,OP_1)) &
             - intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR ),g(:,OP_DZP),h(:,OP_1))
     endif
  endif
  b1psifn2 = temp
#else
  b1psifn2 = 0.
#endif
end function b1psifn2




! B1bbn
! =====
function b1bbn(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bbn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bbn = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if
  else
     if(surface_int) then
        temp79a = ri4_79*f(:,OP_DP)*ni79(:,OP_1)**2
        temp79b = ri4_79*f(:,OP_1)*ni79(:,OP_1)**2
        temp = &
             - intx5(e(:,:,OP_1),temp79a,norm79(:,1),g(:,OP_DR ),h(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79a,norm79(:,2),g(:,OP_DZ ),h(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DRP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DZP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),g(:,OP_DR ),h(:,OP_DP)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),g(:,OP_DZ ),h(:,OP_DP))
     else
        temp79a = ri4_79*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_1 ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_1 ),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_DP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_DP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_1 ),h(:,OP_DP)) &
             - intx5(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_1 ),h(:,OP_DP))
     endif
  endif
  b1bbn = temp
#else
  b1bbn = 0.
#endif
end function b1bbn


! B1bfn1
! ======
function b1bfn1(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfn1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfn1 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp79a = h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DZ),g(:,OP_DZ)) &
             + intx4(e(:,:,OP_1),temp79a,f(:,OP_DR),g(:,OP_DR))
     end if
  else
     if(surface_int) then
        temp79a = (f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)) &
             * ni79(:,OP_1)**2
        temp79b = ri4_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp =-intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),h(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,2),h(:,OP_DZ)) &
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ))&
             + intx5(e(:,:,OP_1),temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))&
             - intx5(e(:,:,OP_DR),ri4_79,temp79a,norm79(:,1),h(:,OP_1)) &
             - intx5(e(:,:,OP_DZ),ri4_79,temp79a,norm79(:,2),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp &
                - 2.*intx5(e(:,:,OP_1),ri5_79,temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp79a = ri2_79*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_GS),temp79a,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             + intx5(e(:,:,OP_GS),temp79a,f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
     endif
  endif
  b1bfn1 = temp
#else
  b1bfn1 = 0.
#endif
end function b1bfn1

! B1bfn2
! ======
function b1bfn2(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1bfn2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1bfn2 = 0.
     return
  end if

  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp79b = ri4_79*h(:,OP_1)*ni79(:,OP_1)**2
        temp79c =-ri4_79*h(:,OP_DP)*ni79(:,OP_1)**2
        temp = &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_DP),norm79(:,1),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),temp79b,f(:,OP_DP),norm79(:,2),g(:,OP_DZP)) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_1 ),norm79(:,1),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),temp79c,f(:,OP_1 ),norm79(:,2),g(:,OP_DZP))

     else
#ifdef USECOMPLEX
        temp79a = ri4_79*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_DZ),temp79a,g(:,OP_DZP),f(:,OP_DP),h(:,OP_1 )) &
             + intx5(e(:,:,OP_DR),temp79a,g(:,OP_DRP),f(:,OP_DP),h(:,OP_1 )) &
             - intx5(e(:,:,OP_DZ),temp79a,g(:,OP_DZP),f(:,OP_1 ),h(:,OP_DP)) &
             - intx5(e(:,:,OP_DR),temp79a,g(:,OP_DRP),f(:,OP_1 ),h(:,OP_DP))

        ! f''' term hack
        temp = temp + rfac* &
             (intx5(e(:,:,OP_DZ),temp79a,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
             +intx5(e(:,:,OP_DR),temp79a,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#elif defined(USE3D)
        temp79a = ri4_79*ni79(:,OP_1)**2
        ! here, we can integrate by parts
        temp = - &
             (intx5(e(:,:,OP_DZP),temp79a,g(:,OP_DZP),f(:,OP_1),h(:,OP_1))  &
             +intx5(e(:,:,OP_DRP),temp79a,g(:,OP_DRP),f(:,OP_1),h(:,OP_1)))
#endif
     endif
  endif
  b1bfn2 = temp
#else
  b1bfn2 = 0.
#endif
end function b1bfn2



! B1pen
! =====
function b1pen(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1pen
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b1pen = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp79a = ni79(:,OP_1)**2
        temp = intx4(e(:,:,OP_1),temp79a,f(:,OP_DP),g(:,OP_1))
     end if
  else
     if(surface_int) then
        temp79a = ri2_79*ni79(:,OP_1)**2
        temp =-intx5(e(:,:,OP_1),temp79a,f(:,OP_DP),norm79(:,1),g(:,OP_DR)) &
             - intx5(e(:,:,OP_1),temp79a,f(:,OP_DP),norm79(:,2),g(:,OP_DZ)) &
             + intx5(e(:,:,OP_1),temp79a,g(:,OP_DP),norm79(:,1),f(:,OP_DR)) &
             + intx5(e(:,:,OP_1),temp79a,g(:,OP_DP),norm79(:,2),f(:,OP_DZ)) &
             - intx5(e(:,:,OP_DR),temp79a,g(:,OP_1),f(:,OP_DP),norm79(:,1)) &
             - intx5(e(:,:,OP_DZ),temp79a,g(:,OP_1),f(:,OP_DP),norm79(:,2))
     else
        temp79a = ri2_79*ni79(:,OP_1)**2
        temp = intx4(e(:,:,OP_GS),temp79a,f(:,OP_DP),g(:,OP_1)) &
             + intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_1 )) &
             + intx4(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_1 )) &
             - intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DZ ),g(:,OP_DP)) &
             - intx4(e(:,:,OP_DR),temp79a,f(:,OP_DR ),g(:,OP_DP))
     end if
  endif
  b1pen = temp
#else
  b1pen = 0.
#endif
end function b1pen



! B1feta
! ======
function b1feta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1feta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp = 0.
        else
           temp = -intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DZP),g(:,OP_1)) &
                  +intx5(e(:,:,OP_1),ri3_79,norm79(:,2),f(:,OP_DRP),g(:,OP_1))
        end if
     else
#ifdef USECOMPLEX
        temp = intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZPP),g(:,OP_1)) &
             - intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRPP),g(:,OP_1))
#else
        temp = -intx4(e(:,:,OP_DRP),ri3_79,f(:,OP_DZP),g(:,OP_1)) &
             +  intx4(e(:,:,OP_DZP),ri3_79,f(:,OP_DRP),g(:,OP_1))
#endif
     end if
  end if

  b1feta = temp
#else
  b1feta = 0.
#endif
end function b1feta

function b1feta1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1feta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

  if(jadv.eq.0) then
     temp%len = 0
  else
     if(surface_int) then
        if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
           temp%len = 0
        else
           temp =  prod(-ri3_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DZP) &
                  +prod( ri3_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DRP)
        end if
     else
#ifdef USECOMPLEX
        temp = prod( ri3_79*g(:,OP_1),OP_DR,OP_DZPP) &
             + prod(-ri3_79*g(:,OP_1),OP_DZ,OP_DRPP)
#else
        temp =  prod(-ri3_79*g(:,OP_1),OP_DRP,OP_DZP) &
             +  prod( ri3_79*g(:,OP_1),OP_DZP,OP_DRP)
#endif
     end if
  end if

  b1feta1 = temp
#else
  b1feta1%len = 0
#endif
end function b1feta1


! b1fu
! ====
function b1fu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = &
                - intx4(e(:,:,OP_1),r2_79,f(:,OP_DZ),g(:,OP_DZ)) &
                - intx4(e(:,:,OP_1),r2_79,f(:,OP_DR),g(:,OP_DR))
        endif
     else
        if(surface_int) then
           temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
           temp = intx3(e(:,:,OP_DR),temp79a,norm79(:,1)) &
                + intx3(e(:,:,OP_DZ),temp79a,norm79(:,1)) &
                - intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
                - intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
                - intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
                - intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
                - intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
                - intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
                - intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
                - intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
           if(itor.eq.1) then
              temp = temp &
                   - 2.*intx4(e(:,:,OP_1),ri_79,temp79a,norm79(:,1))
           endif
           
        else
           temp = &
                - intx3(e(:,:,OP_GS),f(:,OP_DZ),g(:,OP_DZ)) &
                - intx3(e(:,:,OP_GS),f(:,OP_DR),g(:,OP_DR))
        endif
     endif

  b1fu = temp
#else
  b1fu = 0.
#endif
end function b1fu

function b1fu1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1fu1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

     if(jadv.eq.0) then
    !    if(surface_int) then
           !temp%len = 0
        !else
           !temp = &
                !- prod(e(:,:,OP_1),r2_79,f(:,OP_DZP),g(:,OP_DZ)) &
                !- prod(e(:,:,OP_1),r2_79,f(:,OP_DRP),g(:,OP_DR))
    !    endif
     else
        if(surface_int) then
           tempa = mu(g(:,OP_DZ),OP_DZ) + mu(g(:,OP_DR),OP_DR)
           temp = prod(mu(norm79(:,1),OP_DR),tempa) &
                + prod(mu(norm79(:,1),OP_DZ),tempa) &
                + prod(-norm79(:,1)*g(:,OP_DZ ),OP_1,OP_DRZ) &
                + prod(-norm79(:,1)*g(:,OP_DR ),OP_1,OP_DRR) &
                + prod(-norm79(:,1)*g(:,OP_DRZ),OP_1,OP_DZ) &
                + prod(-norm79(:,1)*g(:,OP_DRR),OP_1,OP_DR) &
                + prod(-norm79(:,2)*g(:,OP_DZ ),OP_1,OP_DZZ) &
                + prod(-norm79(:,2)*g(:,OP_DR ),OP_1,OP_DRZ) &
                + prod(-norm79(:,2)*g(:,OP_DZZ),OP_1,OP_DZ) &
                + prod(-norm79(:,2)*g(:,OP_DRZ),OP_1,OP_DR)
           if(itor.eq.1) then
              temp = temp &
                   + prod(mu(-2.*ri_79*norm79(:,1),OP_1),tempa)
           endif
           
        else
           temp = &
                  prod(-g(:,OP_DZ),OP_GS,OP_DZ) &
                + prod(-g(:,OP_DR),OP_GS,OP_DR)
        endif
     endif

  b1fu1 = temp
#else
  b1fu1%len = 0
#endif
end function b1fu1

function b1fu2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1fu2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

     if(jadv.eq.0) then
     !   if(surface_int) then
           !temp%len = 0
        !else
           !temp = &
                !- prod(e(:,:,OP_1),r2_79,f(:,OP_DZP),g(:,OP_DZ)) &
                !- prod(e(:,:,OP_1),r2_79,f(:,OP_DRP),g(:,OP_DR))
     !   endif
     else
        if(surface_int) then
           tempa = mu(f(:,OP_DZ),OP_DZ) + mu(f(:,OP_DR),OP_DR)
           temp = prod(mu(norm79(:,1),OP_DR),tempa) &
                + prod(mu(norm79(:,1),OP_DZ),tempa) &
                + prod(-norm79(:,1)*f(:,OP_DRZ),OP_1,OP_DZ) &
                + prod(-norm79(:,1)*f(:,OP_DRR),OP_1,OP_DR) &
                + prod(-norm79(:,1)*f(:,OP_DZ ),OP_1,OP_DRZ) &
                + prod(-norm79(:,1)*f(:,OP_DR ),OP_1,OP_DRR) &
                + prod(-norm79(:,2)*f(:,OP_DZZ),OP_1,OP_DZ) &
                + prod(-norm79(:,2)*f(:,OP_DRZ),OP_1,OP_DR) &
                + prod(-norm79(:,2)*f(:,OP_DZ ),OP_1,OP_DZZ) &
                + prod(-norm79(:,2)*f(:,OP_DR ),OP_1,OP_DRZ)
           if(itor.eq.1) then
              temp = temp &
                   + prod(mu(-2.*ri_79*norm79(:,1),OP_1),tempa)
           endif
           
        else
           temp = &
                  prod(-f(:,OP_DZ),OP_GS,OP_DZ) &
                + prod(-f(:,OP_DR),OP_GS,OP_DR)
        endif
     endif

  b1fu2 = temp
#else
  b1fu2%len = 0
#endif
end function b1fu2

! b1fv
! ====
function b1fv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

     if(jadv.eq.0) then
        temp = 0.
     else
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp = 0.
           else
           temp = intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DZP),g(:,OP_1 ))&
                - intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DRP),g(:,OP_1 ))&
                + intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DP))&
                - intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DP))
           endif
        else
           temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DRP),g(:,OP_1)) &
                - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZP),g(:,OP_1)) &
                + intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR ),g(:,OP_DP)) &
                - intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ ),g(:,OP_DP))
        endif
     endif


  b1fv = temp
#else
  b1fv = 0.
#endif
end function b1fv

function b1fv1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1fv1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

     if(jadv.eq.0) then
        temp%len = 0
     else
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp%len = 0
           else
           temp = prod( ri_79*norm79(:,1)*g(:,OP_1 ),OP_1,OP_DZP)&
                + prod(-ri_79*norm79(:,2)*g(:,OP_1 ),OP_1,OP_DRP)&
                + prod( ri_79*norm79(:,1)*g(:,OP_DP),OP_1,OP_DZ)&
                + prod(-ri_79*norm79(:,2)*g(:,OP_DP),OP_1,OP_DR)
           endif
        else
           temp = prod( ri_79*g(:,OP_1),OP_DZ,OP_DRP) &
                + prod(-ri_79*g(:,OP_1),OP_DR,OP_DZP) &
                + prod( ri_79*g(:,OP_DP),OP_DZ,OP_DR) &
                + prod(-ri_79*g(:,OP_DP),OP_DR,OP_DZ)
        endif
     endif

  b1fv1 = temp
#else
  b1fv1%len = 0
#endif
end function b1fv1

function b1fv2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1fv2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

     if(jadv.eq.0) then
        temp%len = 0
     else
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp%len = 0
           else
           temp = prod( ri_79*norm79(:,1)*f(:,OP_DZP),OP_1,OP_1)&
                + prod(-ri_79*norm79(:,2)*f(:,OP_DRP),OP_1,OP_1)&
                + prod( ri_79*norm79(:,1)*f(:,OP_DZ ),OP_1,OP_DP)&
                + prod(-ri_79*norm79(:,2)*f(:,OP_DR ),OP_1,OP_DP)
           endif
        else
           temp = prod( ri_79*f(:,OP_DRP),OP_DZ,OP_1) &
                + prod(-ri_79*f(:,OP_DZP),OP_DR,OP_1) &
                + prod( ri_79*f(:,OP_DR ),OP_DZ,OP_DP) &
                + prod(-ri_79*f(:,OP_DZ ),OP_DR,OP_DP)
        endif
     endif

  b1fv2 = temp
#else
  b1fv2%len = 0
#endif
end function b1fv2


! b1fchi
! ======
function b1fchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1fchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = intx4(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR)) &
                - intx4(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ))
        endif
     else
        if(surface_int) then
           temp79a = f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)
           temp = intx4(e(:,:,OP_DR),ri3_79,temp79a,norm79(:,1)) &
                + intx4(e(:,:,OP_DZ),ri3_79,temp79a,norm79(:,2)) &
                + intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
                - intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
                + intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
                - intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
                + intx5(e(:,:,OP_1),ri3_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
                - intx5(e(:,:,OP_1),ri3_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
                + intx5(e(:,:,OP_1),ri3_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
                - intx5(e(:,:,OP_1),ri3_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ))
           if(itor.eq.1) then
              temp = temp &
                   + 2.*intx4(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1))
           endif
           
        else
           temp = intx4(e(:,:,OP_GS),ri3_79,f(:,OP_DZ),g(:,OP_DR)) &
                - intx4(e(:,:,OP_GS),ri3_79,f(:,OP_DR),g(:,OP_DZ))
        endif
     endif

  b1fchi = temp
#else
  b1fchi = 0.
#endif
end function b1fchi

function b1fchi1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1fchi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

     if(jadv.eq.0) then
        !if(surface_int) then
           !temp%len = 0
        !else
           !temp = prod(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR)) &
                !- prod(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ))
        !endif
     else
        if(surface_int) then
           tempa = mu(g(:,OP_DZ),OP_DR) + mu(-g(:,OP_DR),OP_DZ)
           temp = prod(mu(ri3_79*norm79(:,1),OP_DR),tempa) &
                + prod(mu(ri3_79*norm79(:,2),OP_DZ),tempa) &
                + prod( ri3_79*norm79(:,1)*g(:,OP_DR ),OP_1,OP_DRZ) &
                + prod(-ri3_79*norm79(:,1)*g(:,OP_DZ ),OP_1,OP_DRR) &
                + prod( ri3_79*norm79(:,1)*g(:,OP_DRR),OP_1,OP_DZ) &
                + prod(-ri3_79*norm79(:,1)*g(:,OP_DRZ),OP_1,OP_DR) &
                + prod( ri3_79*norm79(:,2)*g(:,OP_DR ),OP_1,OP_DZZ) &
                + prod(-ri3_79*norm79(:,2)*g(:,OP_DZ ),OP_1,OP_DRZ) &
                + prod( ri3_79*norm79(:,2)*g(:,OP_DRZ),OP_1,OP_DZ) &
                + prod(-ri3_79*norm79(:,2)*g(:,OP_DZZ),OP_1,OP_DR)
           if(itor.eq.1) then
              temp = temp &
                   + prod(mu(2.*ri4_79*norm79(:,1),OP_1),tempa) &
                   + prod(mu(2.*ri4_79*norm79(:,1),OP_1),tempa)
           endif
           
        else
           temp = prod( ri3_79*g(:,OP_DR),OP_GS,OP_DZ) &
                + prod(-ri3_79*g(:,OP_DZ),OP_GS,OP_DR)
         endif
     endif

  b1fchi1 = temp
#else
  b1fchi1%len = 0
#endif
end function b1fchi1

function b1fchi2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b1fchi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

     if(jadv.eq.0) then
        !if(surface_int) then
           !temp%len = 0
        !else
           !temp = prod(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR)) &
                !- prod(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ))
        !endif
     else
        if(surface_int) then
           tempa = mu(f(:,OP_DR),OP_DZ) + mu(-f(:,OP_DZ),OP_DR)
           temp = prod(mu(ri3_79*norm79(:,1),OP_DR),tempa) &
                + prod(mu(ri3_79*norm79(:,2),OP_DZ),tempa) &
                + prod( ri3_79*norm79(:,1)*f(:,OP_DRZ),OP_1,OP_DR) &
                + prod(-ri3_79*norm79(:,1)*f(:,OP_DRR),OP_1,OP_DZ) &
                + prod( ri3_79*norm79(:,1)*f(:,OP_DZ ),OP_1,OP_DRR) &
                + prod(-ri3_79*norm79(:,1)*f(:,OP_DR ),OP_1,OP_DRZ) &
                + prod( ri3_79*norm79(:,2)*f(:,OP_DZZ),OP_1,OP_DR) &
                + prod(-ri3_79*norm79(:,2)*f(:,OP_DRZ),OP_1,OP_DZ) &
                + prod( ri3_79*norm79(:,2)*f(:,OP_DZ ),OP_1,OP_DRZ) &
                + prod(-ri3_79*norm79(:,2)*f(:,OP_DR ),OP_1,OP_DZZ)
           if(itor.eq.1) then
              temp = temp &
                   + prod(mu(2.*ri4_79*norm79(:,1),OP_1),tempa) &
                   + prod(mu(2.*ri4_79*norm79(:,1),OP_1),tempa)
          endif
           
        else
           temp = prod( ri3_79*f(:,OP_DZ),OP_GS,OP_DR) &
                + prod(-ri3_79*f(:,OP_DR),OP_GS,OP_DZ)
        endif
     endif

  b1fchi2 = temp
#else
  b1fchi2%len = 0
#endif
end function b1fchi2

! B1e
! ===
function b1e(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1e
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(jadv.eq.1) then
     temp = 0.
  else
     if(surface_int) then
        temp = 0.
     else
        temp = -intx2(e(:,:,OP_1),f(:,OP_DP))
     end if
  endif
  b1e = temp
#else
  b1e = 0.
#endif
end function b1e

function b1vzdot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1vzdot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
        temp = 0.
  else
     if(jadv.eq.0) then
        temp = intx3(e(:,:,OP_1),r2_79,f(:,OP_1))
     else
        temp = -intx2(e(:,:,OP_DR),f(:,OP_DR))    &
               -intx2(e(:,:,OP_DZ),f(:,OP_DZ))
        if(itor.eq.1) temp = temp - 2.*intx3(e(:,:,OP_DR),ri_79,f(:,OP_1))
     endif

  end if

  b1vzdot = temp
end function b1vzdot

function b1chidot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1chidot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
        temp = 0.
  else
     if(jadv.eq.0) then
        temp = 0.
     else
        temp = -intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRP))    &
               -intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP))   
     endif

  end if
  b1chidot = temp
#else
  b1chidot = 0.
#endif
end function b1chidot


function b1psi2bfpe(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b1psi2bfpe
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j
  vectype, dimension(dofs_per_element) :: temp
  temp = 0.

  temp79a = ri_79*(j(:,OP_DZ)*f(:,OP_DR) - j(:,OP_DR)*f(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = temp79a - (j(:,OP_DR)*i(:,OP_DR) + j(:,OP_DZ)*i(:,OP_DZ))   &
                    + ri2_79*h(:,OP_1)*j(:,OP_DP)
#endif
  temp79a = temp79a*b2i79(:,OP_1)*ri2_79*ni79(:,OP_1)


  if(surface_int) then
        temp = 0.
  else
     if(jadv.eq.0) then
        temp = 0.
     else
        temp = intx3(e(:,:,OP_GS),temp79a,h(:,OP_1))
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp - intx4(e(:,:,OP_DZP),ri_79,temp79a,g(:,OP_DR))   &
                    + intx4(e(:,:,OP_DRP),ri_79,temp79a,g(:,OP_DZ))   &
                    + intx3(e(:,:,OP_DRP),temp79a,i(:,OP_DR))        &
                    + intx3(e(:,:,OP_DZP),temp79a,i(:,OP_DZ))
#endif
     endif

  end if

  b1psi2bfpe = temp
end function b1psi2bfpe


!==============================================================================
! B2 TERMS
!==============================================================================

! B2b
! ===
!function b2b(e,f)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: b2b
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  !vectype, dimension(dofs_per_element) :: temp

  !if(surface_int) then
     !temp = 0.
  !else
     !temp = intx3(e(:,:,OP_1),ri2_79,f(:,OP_1))
  !end if

  !b2b = temp
!end function b2b

function b2b
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2b
  type(prodarray) :: temp

  if(surface_int) then
     temp%len = 0
  else
     temp = prod(ri2_79,OP_1,OP_1)
  end if

  b2b = temp
end function b2b


! B2psieta
! ========
function b2psieta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp = intx5(e(:,:,OP_1),ri3_79,norm79(:,1),f(:,OP_DZP),g(:,OP_1)) &
             - intx5(e(:,:,OP_1),ri3_79,norm79(:,2),f(:,OP_DRP),g(:,OP_1))
     endif
  else
     temp = intx4(e(:,:,OP_DZ),ri3_79,f(:,OP_DRP),g(:,OP_1)) &
          - intx4(e(:,:,OP_DR),ri3_79,f(:,OP_DZP),g(:,OP_1))

#ifndef USEST
     if(hypi.ne.0 .and. imp_hyper.le.1) then
        if(ihypeta.eq.0) then          
           temp = temp + 2.*hypi* &
                (intx3(e(:,:,OP_DZZ),ri3_79,f(:,OP_DRZP)) &
                -intx3(e(:,:,OP_DRR),ri3_79,f(:,OP_DRZP)) &
                -intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZZP)) &
                +intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DRRP)))
           
           if(itor.eq.1) then
              temp = temp - 2.*hypi* &
                   (   intx3(e(:,:,OP_DZZ),ri4_79,f(:,OP_DZP)) &
                   -   intx3(e(:,:,OP_DRR),ri4_79,f(:,OP_DZP)) &
                   +   intx3(e(:,:,OP_DR),ri5_79,f(:,OP_DZP)) &
                   +2.*intx3(e(:,:,OP_DRZ),ri4_79,f(:,OP_DRP)) &
                   -4.*intx3(e(:,:,OP_DZ),ri5_79,f(:,OP_DRP)) &
                   -   intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRZP)) &
                   +   intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZZP)) &
                   +2.*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DRRP)))
           endif
        endif
     end if
#endif !USEST
  end if
  b2psieta = temp
#else
  b2psieta = 0.
#endif
end function b2psieta

function b2psieta1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2psieta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray):: temp

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp%len = 0
     else
        temp = prod( ri3_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DZP) &
             + prod(-ri3_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DRP)
      endif
  else
     temp = prod( ri3_79*g(:,OP_1),OP_DZ,OP_DRP) &
          + prod(-ri3_79*g(:,OP_1),OP_DR,OP_DZP)

#ifndef USEST
     if(hypi.ne.0 .and. imp_hyper.le.1) then
        if(ihypeta.eq.0) then          
           temp = temp + &
                (prod( 2.*hypi*ri3_79,OP_DZZ,OP_DRZP) &
                +prod(-2.*hypi*ri3_79,OP_DRR,OP_DRZP) &
                +prod(-2.*hypi*ri3_79,OP_DRZ,OP_DZZP) &
                +prod( 2.*hypi*ri3_79,OP_DRZ,OP_DRRP))

           if(itor.eq.1) then
             temp = temp + &
                   (   prod(-2.*hypi*ri4_79,OP_DZZ,OP_DZP) &
                   +   prod( 2.*hypi*ri4_79,OP_DRR,OP_DZP) &
                   +   prod(-2.*hypi*ri5_79,OP_DR,OP_DZP) &
                   +   prod(-4.*hypi*ri4_79,OP_DRZ,OP_DRP) &
                   +   prod(+8.*hypi*ri5_79,OP_DZ,OP_DRP) &
                   +   prod( 2.*hypi*ri4_79,OP_DR,OP_DRZP) &
                   +   prod(-2.*hypi*ri4_79,OP_DZ,OP_DZZP) &
                   +   prod(-4.*hypi*ri4_79,OP_DZ,OP_DRRP))
            endif
        endif
     end if
#endif !USEST
  end if
  b2psieta1 = temp
#else
  b2psieta1%len = 0
#endif
end function b2psieta1


! B2psimue
! ========
vectype function b2psimue(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
  
  if(surface_int) then
     temp = 0.
  else
     temp = &
          - int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
          - int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR)) &
          - int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_1 ))
  end if

  b2psimue = temp
  return
end function b2psimue



! B2beta
! ======
function b2beta(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2beta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        ! better to exclude this term
!!$        temp = 0.
        temp = intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
             + intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZ),g(:,OP_1))
     end if
  else
     temp = &
          - intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
          - intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_1)) 
#if defined(USE3D) || defined(USECOMPLEX)
     if(iupstream.eq.1) then    
        temp79a = abs(h(:,OP_1))*magus
        temp = temp + intx4(e(:,:,OP_1),ri4_79,f(:,OP_DPP),temp79a)
     elseif(iupstream.eq.2) then
        temp79a = abs(h(:,OP_1))*magus
        temp = temp - intx4(e(:,:,OP_DPP),ri6_79,f(:,OP_DPP),temp79a)
     endif
#endif     

     if(imp_hyper.le.1) then
!    the following coding should be checked.  It does not agree with my derivation  scj 4/30/14
        if(hypi.ne.0.) then
           if(ihypeta.eq.1) then
              temp = temp - hypi* &
                   (intx4(e(:,:,OP_GS),ri2_79,f(:,OP_GS),g(:,OP_1)) &
                   +intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),g(:,OP_DZ)) &
                   +intx4(e(:,:,OP_DR),ri2_79,f(:,OP_GS),g(:,OP_DR)))
                   
           else
              temp = temp + hypi* &
                   (-intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DZZ)) &
                   + intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DZZ)) &
                   + intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DRR)) &
                   - intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DRR)) &
                   - 4.*intx3(e(:,:,OP_DRZ),ri2_79,f(:,OP_DRZ)))

              if(itor.eq.1) then
                 temp = temp + hypi*&
                      (-intx3(e(:,:,OP_DZZ),ri3_79,f(:,OP_DR)) &
                      + intx3(e(:,:,OP_DRR),ri3_79,f(:,OP_DR)) &
                      - intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DR)) &
                      - 2.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZ)) &
                      + 4.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZ)) &
                      - 4.*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZ)) &
                      - intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DZZ)) &
                      + intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DRR)) &
                      + 2.*intx3(e(:,:,OP_DZ),ri3_79,f(:,OP_DRZ)))
              endif

#if defined(USE3D) || defined(USECOMPLEX)
              temp = temp &
#ifdef USEST
                   - hypi*intx3(e(:,:,OP_DZP),ri4_79,f(:,OP_DZP)) &
                   - hypi*intx3(e(:,:,OP_DRP),ri4_79,f(:,OP_DRP))
#else
                   + hypi*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZPP)) &
                   + hypi*intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRPP))
#endif
#endif
           end if
        endif
     endif
  end if

  b2beta = temp
end function b2beta

function b2beta1(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2beta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp%len = 0
     else
        ! better to exclude this term
!!$        temp = 0.
        temp = prod(ri2_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DR) &
             + prod(ri2_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DZ)
       end if
  else
     temp = &
            prod(-ri2_79*g(:,OP_1),OP_DZ,OP_DZ) &
          + prod(-ri2_79*g(:,OP_1),OP_DR,OP_DR)


     if(imp_hyper.le.1) then
!    the following coding should be checked.  It does not agree with my derivation  scj 4/30/14
        if(hypi.ne.0.) then
           if(ihypeta.eq.1) then
                   
           else
              temp = temp + &
                   ( prod(-hypi*ri2_79,OP_DZZ,OP_DZZ) &
                   + prod( hypi*ri2_79,OP_DRR,OP_DZZ) &
                   + prod( hypi*ri2_79,OP_DZZ,OP_DRR) &
                   + prod(-hypi*ri2_79,OP_DRR,OP_DRR) &
                   + prod(-4.*hypi*ri2_79,OP_DRZ,OP_DRZ))

              if(itor.eq.1) then
                 temp = temp + &
                      ( prod(-hypi*ri3_79,OP_DZZ,OP_DR) &
                      + prod( hypi*ri3_79,OP_DRR,OP_DR) &
                      + prod(-hypi*ri4_79,OP_DR,OP_DR) &
                      + prod(-2.*hypi*ri3_79,OP_DRZ,OP_DZ) &
                      + prod(4.*hypi*ri3_79,OP_DRZ,OP_DZ) &
                      + prod(-4.*hypi*ri4_79,OP_DZ,OP_DZ) &
                      + prod(-hypi*ri3_79,OP_DR,OP_DZZ) &
                      + prod( hypi*ri3_79,OP_DR,OP_DRR) &
                      + prod(2.*hypi*ri3_79,OP_DZ,OP_DRZ))
               endif

#if defined(USE3D) || defined(USECOMPLEX)
              temp = temp &
#ifdef USEST
                   + prod(-hypi*ri4_79,OP_DZP,OP_DZP) &
                   + prod(-hypi*ri4_79,OP_DRP,OP_DRP)
#else
                   + prod(hypi*ri4_79,OP_DZ,OP_DZPP) &
                   + prod(hypi*ri4_79,OP_DR,OP_DRPP)
#endif
#endif

           end if
        endif
     endif
  end if

  b2beta1 = temp
end function b2beta1


! B2feta
! ======
function b2feta(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2feta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp = intx5(e(:,:,OP_1),ri2_79,norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
             + intx5(e(:,:,OP_1),ri2_79,norm79(:,2),f(:,OP_DZP),g(:,OP_1))
     end if
  else
     temp = &
          - intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZP),g(:,OP_1)) &
          - intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DRP),g(:,OP_1))

#ifndef USEST
     if(imp_hyper.le.1) then

!   the following coding should be checked.  does not agree with my derivation scj 4/30/2014
        if(hypi.ne.0.) then
           if(ihypeta.eq.0) then
              temp = temp + hypi*&
                   (-intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DZZP)) &
                   + intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DZZP)) &
                   + intx3(e(:,:,OP_DZZ),ri2_79,f(:,OP_DRRP)) &
                   - intx3(e(:,:,OP_DRR),ri2_79,f(:,OP_DRRP)) &
                   - 4.*intx3(e(:,:,OP_DRZ),ri2_79,f(:,OP_DRZP)))

              if(itor.eq.1) then
                 temp = temp + hypi* &
                      (-intx3(e(:,:,OP_DZZ),ri3_79,f(:,OP_DRP)) &
                      + intx3(e(:,:,OP_DRR),ri3_79,f(:,OP_DRP)) &
                      - intx3(e(:,:,OP_DR),ri4_79,f(:,OP_DRP)) &
                      - 2.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZP)) &
                      + 4.*intx3(e(:,:,OP_DRZ),ri3_79,f(:,OP_DZP)) &
                      - 4.*intx3(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP)) &
                      - intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DZZP)) &
                      + intx3(e(:,:,OP_DR),ri3_79,f(:,OP_DRRP)) &
                      + 2.*intx3(e(:,:,OP_DZ),ri3_79,f(:,OP_DRZP)))
              endif
           endif
        endif

     end if
#endif !USEST
  end if

  b2feta = temp
#else
  b2feta = 0.
#endif
end function b2feta

function b2feta1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2feta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp%len = 0
     else
        temp = prod(ri2_79*norm79(:,1)*g(:,OP_1),OP_1,OP_DRP) &
             + prod(ri2_79*norm79(:,2)*g(:,OP_1),OP_1,OP_DZP)
     end if
  else
     temp = &
            prod(-ri2_79*g(:,OP_1),OP_DZ,OP_DZP) &
          + prod(-ri2_79*g(:,OP_1),OP_DR,OP_DRP)

#ifndef USEST
     if(imp_hyper.le.1) then

!   the following coding should be checked.  does not agree with my derivation scj 4/30/2014
        if(hypi.ne.0.) then
           if(ihypeta.eq.0) then
              temp = temp + &
                   ( prod(-ri2_79*hypi,OP_DZZ,OP_DZZP) &
                   + prod( ri2_79*hypi,OP_DRR,OP_DZZP) &
                   + prod( ri2_79*hypi,OP_DZZ,OP_DRRP) &
                   + prod(-ri2_79*hypi,OP_DRR,OP_DRRP) &
                   + prod(-4.*ri2_79*hypi,OP_DRZ,OP_DRZP))

              if(itor.eq.1) then
                 temp = temp + &
                      ( prod(-ri3_79*hypi,OP_DZZ,OP_DRP) &
                      + prod( ri3_79*hypi,OP_DRR,OP_DRP) &
                      + prod(-ri4_79*hypi,OP_DR,OP_DRP) &
                      + prod(-2.*ri3_79*hypi,OP_DRZ,OP_DZP) &
                      + prod( 4.*ri3_79*hypi,OP_DRZ,OP_DZP) &
                      + prod(-4.*ri4_79*hypi,OP_DZ,OP_DZP) &
                      + prod(-ri3_79*hypi,OP_DR,OP_DZZP) &
                      + prod( ri3_79*hypi,OP_DR,OP_DRRP) &
                      + prod(2.*ri3_79*hypi,OP_DZ,OP_DRZP))
              endif
           endif
        endif

     end if
#endif !USEST
  end if

  b2feta1 = temp
#else
  b2feta1%len = 0
#endif
end function b2feta1

! B2FJ
! ====
function b2fj(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2fj
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  real :: hypfm

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp = 0
     end if
  else
     temp = 0
     if     (ihypeta.eq.1) then
        temp79a = hypf*eta79(:,OP_1)*g(:,OP_DR)
        temp79b = hypf*eta79(:,OP_1)*g(:,OP_DZ)
        temp79c = hypf*(eta79(:,OP_DP)*g(:,OP_DP) + eta79(:,OP_1)*g(:,OP_DPP))
     else if(ihypeta.eq.2) then
        temp79a = hypf*pt79(:,OP_1)*g(:,OP_DR)
        temp79b = hypf*pt79(:,OP_1)*g(:,OP_DZ)
        temp79c = hypf*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
     else if(ihypeta.gt.2) then
        hypfm = hypf*(bharhypeta)**beta
        temp79a = hypfm*pt79(:,OP_1)*g(:,OP_DR)
        temp79b = hypfm*pt79(:,OP_1)*g(:,OP_DZ)
        temp79c = hypfm*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
     else
        temp79a = hypf*g(:,OP_DR)
        temp79b = hypf*g(:,OP_DZ)
        temp79c = hypf*g(:,OP_DPP)
     endif
     
     temp = -intx5(e(:,:,OP_DZ),b2i79(:,OP_DR),ri_79,f(:,OP_DR ),temp79a) &
          +intx5(e(:,:,OP_DR ),b2i79(:,OP_DR),ri_79,f(:,OP_DZ ),temp79a) &
          -intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri_79,f(:,OP_DR ),temp79a) &
          +intx5(e(:,:,OP_DRR),b2i79(:,OP_1 ),ri_79,f(:,OP_DZ ),temp79a) &
          -intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri_79,f(:,OP_DRR),temp79a) &
          +intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri_79,f(:,OP_DRZ),temp79a) &
          -intx5(e(:,:,OP_DZ ),b2i79(:,OP_DZ),ri_79,f(:,OP_DR ),temp79b) &
          +intx5(e(:,:,OP_DR ),b2i79(:,OP_DZ),ri_79,f(:,OP_DZ ),temp79b) &
          -intx5(e(:,:,OP_DZZ),b2i79(:,OP_1 ),ri_79,f(:,OP_DR ),temp79b) &
          +intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri_79,f(:,OP_DZ ),temp79b) &
          -intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri_79,f(:,OP_DRZ),temp79b) &
          +intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri_79,f(:,OP_DZZ),temp79b) 
     if(itor.eq.1) temp = temp &
          +intx5(e(:,:,OP_DZ ),b2i79(:,OP_1),ri2_79,f(:,OP_DR ),temp79a) &
          -intx5(e(:,:,OP_DR ),b2i79(:,OP_1),ri2_79,f(:,OP_DZ ),temp79a) 
     temp = temp                                                         &
          +intx5(e(:,:,OP_DZ ),b2i79(:,OP_1),ri3_79,f(:,OP_DR ),temp79c) &
          -intx5(e(:,:,OP_DR ),b2i79(:,OP_1),ri3_79,f(:,OP_DZ ),temp79c) 
  endif

  b2fj = temp
#else
  b2fj = 0.
#endif
end function b2fj

! B2PSIJ
! ======
function b2psij(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psij
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp
  real :: hypfm

  if(surface_int) then
     if(inocurrent_pol.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp = 0
     end if
  else
   
     if     (ihypeta.eq.1) then
        temp79a = hypf*eta79(:,OP_1)*g(:,OP_DR)
        temp79b = hypf*eta79(:,OP_1)*g(:,OP_DZ)
     else if(ihypeta.eq.2) then
        temp79a = hypf*pt79(:,OP_1)*g(:,OP_DR)
        temp79b = hypf*pt79(:,OP_1)*g(:,OP_DZ)
     else if(ihypeta.gt.2) then
        hypfm = hypf*(bharhypeta)**beta
        temp79a = hypfm*pt79(:,OP_1)*g(:,OP_DR)
        temp79b = hypfm*pt79(:,OP_1)*g(:,OP_DZ)
     else
        temp79a = hypf*g(:,OP_DR)
        temp79b = hypf*g(:,OP_DZ)
     endif
     
     temp = -intx5(e(:,:,OP_DR ),b2i79(:,OP_DR),ri2_79,f(:,OP_DR ),temp79a) &
          -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_DR),ri2_79,f(:,OP_DZ ),temp79a) &
          -  intx5(e(:,:,OP_DRR),b2i79(:,OP_1 ),ri2_79,f(:,OP_DR ),temp79a) &
          -  intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DZ ),temp79a) &
          -  intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DRR),temp79a) &
          -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DRZ),temp79a) &
          -  intx5(e(:,:,OP_DR ),b2i79(:,OP_DZ),ri2_79,f(:,OP_DR ),temp79b) &
          -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_DZ),ri2_79,f(:,OP_DZ ),temp79b) &
          -  intx5(e(:,:,OP_DRZ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DR ),temp79b) &
          -  intx5(e(:,:,OP_DZZ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DZ ),temp79b) &
          -  intx5(e(:,:,OP_DR ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DRZ),temp79b) &
          -  intx5(e(:,:,OP_DZ ),b2i79(:,OP_1 ),ri2_79,f(:,OP_DZZ),temp79b) 
     if(itor.eq.1) temp = temp &
          +2.*intx5(e(:,:,OP_DR ),b2i79(:,OP_1),ri3_79,f(:,OP_DR ),temp79a) &
          +2.*intx5(e(:,:,OP_DZ ),b2i79(:,OP_1),ri3_79,f(:,OP_DZ ),temp79a) 
#if defined(USE3D) || defined(USECOMPLEX)
     if     (ihypeta.eq.1) then
        temp79c = hypf*(eta79(:,OP_DP)*g(:,OP_DP) + eta79(:,OP_1)*g(:,OP_DPP))
     else if(ihypeta.eq.2) then
        temp79c = hypf*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
    else if(ihypeta.gt.2) then
        temp79c = hypfm*(pt79(:,OP_DP)*g(:,OP_DP) + pt79(:,OP_1)*g(:,OP_DPP))
     else
        temp79c = hypf*g(:,OP_DPP)
     endif
     temp = temp                                                        &
          +intx5(e(:,:,OP_DR),b2i79(:,OP_1),ri4_79,f(:,OP_DR ),temp79c) &
          +intx5(e(:,:,OP_DZ),b2i79(:,OP_1),ri4_79,f(:,OP_DZ ),temp79c)
#endif
  end if

  b2psij = temp
end function b2psij

! B2bu
! ====
function b2bu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2bu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           ! this must be included
           temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
                - intx5(e(:,:,OP_1),ri_79,f(:,OP_1),norm79(:,2),g(:,OP_DR))
        endif
     else
        temp = intx4(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_DR)) &
             - intx4(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_DZ))
     endif

  b2bu = temp
end function b2bu

function b2bu1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2bu1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           ! this must be included
           !temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
           !     - intx5(e(:,:,OP_1),ri_79,f(:,OP_1),norm79(:,2),g(:,OP_DR))
        endif
     else
        temp = prod( ri_79*g(:,OP_DR),OP_DZ,OP_1) &
             + prod(-ri_79*g(:,OP_DZ),OP_DR,OP_1)
      endif

  b2bu1 = temp
end function b2bu1

function b2bu2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2bu2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           ! this must be included
           temp = prod( ri_79*f(:,OP_1)*norm79(:,1),OP_1,OP_DZ) &
                + prod(-ri_79*f(:,OP_1)*norm79(:,2),OP_1,OP_DR)
        endif
     else
        temp = prod( ri_79*f(:,OP_1),OP_DZ,OP_DR) &
             + prod(-ri_79*f(:,OP_1),OP_DR,OP_DZ)
     endif

  b2bu2 = temp
end function b2bu2


! B2bchi
! ======
function b2bchi(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2bchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           ! this must be included
           temp = &
                - intx5(e(:,:,OP_1),ri4_79,f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
                - intx5(e(:,:,OP_1),ri4_79,f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        endif
     else
        temp = intx4(e(:,:,OP_DZ),ri4_79,f(:,OP_1),g(:,OP_DZ)) &
             + intx4(e(:,:,OP_DR),ri4_79,f(:,OP_1),g(:,OP_DR))        
     end if

  b2bchi = temp
end function b2bchi

function b2bchi1(g)
  use basic
  use m3dc1_nint

  type(prodarray) :: b2bchi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           ! this must be included
           !temp = &
                !- intx5(e(:,:,OP_1),ri4_79,f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
           !     - intx5(e(:,:,OP_1),ri4_79,f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        endif
     else
        temp = prod(ri4_79*g(:,OP_DZ),OP_DZ,OP_1) &
             + prod(ri4_79*g(:,OP_DR),OP_DR,OP_1)        
     end if

  b2bchi1 = temp
end function b2bchi1


function b2bchi2(f)
  use basic
  use m3dc1_nint

  type(prodarray) :: b2bchi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           ! this must be included
           temp = &
                  prod(-ri4_79*f(:,OP_1)*norm79(:,1),OP_1,OP_DR) &
                + prod(-ri4_79*f(:,OP_1)*norm79(:,2),OP_1,OP_DZ)
         endif
     else
        temp = prod(ri4_79*f(:,OP_1),OP_DZ,OP_DZ) &
             + prod(ri4_79*f(:,OP_1),OP_DR,OP_DR)        
      end if

  b2bchi2 = temp
end function b2bchi2


! B2bd
! ====
!function b2bd(e,f,g)
  !use basic
  !use m3dc1_nint

  !vectype, dimension(dofs_per_element) :: b2bd
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  !vectype, dimension(dofs_per_element) :: temp

  !if(mass_ratio.eq.0. .or. db.eq.0.) then
     !b2bd = 0.
     !return
  !endif

  !if(surface_int) then
     !temp = 0.
  !else
     !temp = - &
          !(intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_DZ),g(:,OP_1)) &
          !+intx4(e(:,:,OP_DR),ri2_79,f(:,OP_DR),g(:,OP_1)))
  !end if

  !b2bd = temp*me_mp*mass_ratio*db**2
!end function b2bd

function b2bd(g)
  use basic
  use m3dc1_nint

  type(prodarray) :: b2bd
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(mass_ratio.eq.0. .or. db.eq.0.) then
     b2bd%len = 0
     return
  endif

  if(surface_int) then
     temp%len = 0
  else
     temp = &
          (prod(-ri2_79*g(:,OP_1),OP_DZ,OP_DZ) &
          +prod(-ri2_79*g(:,OP_1),OP_DR,OP_DR))
  end if

  temp79a=me_mp*mass_ratio*db**2
  b2bd = temp*temp79a
end function b2bd

! B2psiv
! ======
function b2psiv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psiv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),ri_79,g(:,OP_1),norm79(:,2),f(:,OP_DR)) &
                - intx5(e(:,:,OP_1),ri_79,g(:,OP_1),norm79(:,1),f(:,OP_DZ))
        endif
     else
        temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
             - intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1))
     endif

  b2psiv = temp
end function b2psiv

function b2psiv1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2psiv1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp%len = 0
        else
           !temp = intx5(e(:,:,OP_1),ri_79,g(:,OP_1),norm79(:,2),f(:,OP_DR)) &
           !     - intx5(e(:,:,OP_1),ri_79,g(:,OP_1),norm79(:,1),f(:,OP_DZ))
        endif
     else
        temp = prod( ri_79*g(:,OP_1),OP_DR,OP_DZ) &
             + prod(-ri_79*g(:,OP_1),OP_DZ,OP_DR)
     endif

  b2psiv1 = temp
end function b2psiv1

function b2psiv2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2psiv2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp%len = 0
        endif
     else
        temp = prod( ri_79*f(:,OP_DZ),OP_DR,OP_1) &
             + prod(-ri_79*f(:,OP_DR),OP_DZ,OP_1)
     endif

  b2psiv2 = temp
end function b2psiv2


! B2fv
! ====
function b2fv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2fv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp = 0.
        else
           temp = intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_1)) &
                - intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1))
        endif
     else
        temp = intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
             + intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1))
     endif

#else
  temp = 0.
#endif

  b2fv = temp
end function b2fv

function b2fv1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2fv1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp%len = 0
        else
           temp = prod( norm79(:,2)*g(:,OP_1),OP_1,OP_DR) &
                + prod(-norm79(:,1)*g(:,OP_1),OP_1,OP_DZ)
        endif
     else
        temp = prod(g(:,OP_1),OP_DZ,OP_DZ) &
             + prod(g(:,OP_1),OP_DR,OP_DR)
     endif
#else
  temp%len = 0
#endif

  b2fv1 = temp
end function b2fv1

function b2fv2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b2fv2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp%len = 0
        else
           temp = prod( norm79(:,2)*f(:,OP_DR),OP_1,OP_1) &
                + prod(-norm79(:,1)*f(:,OP_DZ),OP_1,OP_1)
        endif
     else
        temp = prod(f(:,OP_DZ),OP_DZ,OP_1) &
             + prod(f(:,OP_DR),OP_DR,OP_1)
     endif
#else
  temp%len = 0
#endif

  b2fv2 = temp
end function b2fv2


! B2psipsid
! =========
function b2psipsid(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psipsid
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psipsid = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_tor.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = f(:,OP_GS)*h(:,OP_1)
        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),g(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),g(:,OP_DZ))
     end if
  else
     temp = intx5(e(:,:,OP_DR),ri3_79,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
          - intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))
  endif

  b2psipsid = temp
end function b2psipsid


! B2psibd
! =======
function b2psibd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psibd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psibd = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = g(:,OP_1)*h(:,OP_1)
        temp = intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),f(:,OP_DRP)) &
             + intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,2),f(:,OP_DZP))
     endif
  else
     temp = &
          -(intx5(e(:,:,OP_DZ),ri4_79,f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
           +intx5(e(:,:,OP_DR),ri4_79,f(:,OP_DRP),g(:,OP_1),h(:,OP_1)))
  end if

  b2psibd = temp
#else
  b2psibd = 0.
#endif
end function b2psibd


! B2bbd
! =====
function b2bbd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2bbd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bbd = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = g(:,OP_1)*h(:,OP_1)
        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),f(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),f(:,OP_DZ))
     endif
  else
     temp = intx5(e(:,:,OP_DR),ri3_79,f(:,OP_DZ),g(:,OP_1),h(:,OP_1)) &
          - intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_DR),g(:,OP_1),h(:,OP_1))
  endif
  
  b2bbd = temp
end function b2bbd



! B2ped
! =====
function b2ped(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2ped
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2ped = 0.
     return
  end if

  if(surface_int) then
     temp = intx5(e(:,:,OP_1),ri_79,norm79(:,2),f(:,OP_DR),g(:,OP_1)) &
          - intx5(e(:,:,OP_1),ri_79,norm79(:,1),f(:,OP_DZ),g(:,OP_1))
  else
     temp = intx4(e(:,:,OP_DR),ri_79,f(:,OP_DZ),g(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),ri_79,f(:,OP_DR),g(:,OP_1))
  end if

  b2ped = temp
end function b2ped


! B2psifd
! =======
function b2psifd(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2psifd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psifd = 0.
     return
  end if
  
  if(surface_int) then
     if(inocurrent_tor.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = f(:,OP_GS)*h(:,OP_1)
        temp = &
             - intx5(e(:,:,OP_1),ri2_79,temp79a,norm79(:,1),g(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri2_79,temp79a,norm79(:,2),g(:,OP_DZ))
     endif
  else
     temp = intx5(e(:,:,OP_DZ),ri2_79,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
          + intx5(e(:,:,OP_DR),ri2_79,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))
  end if
  b2psifd = temp
#else
  b2psifd = 0.
#endif
end function b2psifd


! B2bfd
! =====
function b2bfd(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2bfd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bfd = 0.
     return
  end if
  
  if(surface_int) then
     if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else 
        temp79a = f(:,OP_1)*h(:,OP_1)
        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),g(:,OP_DZP))
     endif
  else
     temp = - &
          (intx5(e(:,:,OP_DZ),ri3_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_1)) &
          -intx5(e(:,:,OP_DR),ri3_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_1)))
  end if
  
  b2bfd = temp
#else
  b2bfd = 0.
#endif
end function b2bfd

! B2psipsin
! =========
function b2psipsin(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psipsin
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psipsin = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_tor.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = f(:,OP_GS)*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),g(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),g(:,OP_DZ))
     end if
  else
     temp79a = ri3_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
          - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))
  endif

  b2psipsin = temp
end function b2psipsin


! B2psibn
! =======
function b2psibn(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psibn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psibn = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = g(:,OP_1)*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,1),f(:,OP_DRP)) &
             + intx5(e(:,:,OP_1),ri4_79,temp79a,norm79(:,2),f(:,OP_DZP))
     endif
  else
     temp79a = ri4_79*ni79(:,OP_1)**2
     temp = &
          -(intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
           +intx5(e(:,:,OP_DR),temp79a,f(:,OP_DRP),g(:,OP_1),h(:,OP_1)))
  end if

  b2psibn = temp
#else
  b2psibn = 0.
#endif
end function b2psibn


! B2bbn
! =====
function b2bbn(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2bbn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bbn = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = g(:,OP_1)*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),f(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),f(:,OP_DZ))
     endif
  else
     temp79a = ri3_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_DR),temp79a,f(:,OP_DZ),g(:,OP_1),h(:,OP_1)) &
          - intx5(e(:,:,OP_DZ),temp79a,f(:,OP_DR),g(:,OP_1),h(:,OP_1))
  endif
  
  b2bbn = temp
end function b2bbn



! B2pen
! =====
function b2pen(e,f,g)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2pen
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2pen = 0.
     return
  end if

  temp79a = ri_79*ni79(:,OP_1)**2
  if(surface_int) then
     temp = intx5(e(:,:,OP_1),temp79a,norm79(:,2),f(:,OP_DR),g(:,OP_1)) &
          - intx5(e(:,:,OP_1),temp79a,norm79(:,1),f(:,OP_DZ),g(:,OP_1))
  else
     temp = intx4(e(:,:,OP_DR),temp79a,f(:,OP_DZ),g(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DR),g(:,OP_1))
  end if

  b2pen = temp
  return
end function b2pen


! B2psifn
! =======
function b2psifn(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2psifn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2psifn = 0.
     return
  end if
  
  if(surface_int) then
     if(inocurrent_tor.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else
        temp79a = f(:,OP_GS)*h(:,OP_1)*ni79(:,OP_1)**2
        temp = &
             - intx5(e(:,:,OP_1),ri2_79,temp79a,norm79(:,1),g(:,OP_DR)) &
             - intx5(e(:,:,OP_1),ri2_79,temp79a,norm79(:,2),g(:,OP_DZ))
     endif
  else
     temp79a = ri2_79*ni79(:,OP_1)**2
     temp = intx5(e(:,:,OP_DZ),temp79a,f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
          + intx5(e(:,:,OP_DR),temp79a,f(:,OP_GS),g(:,OP_DR),h(:,OP_1))
  end if
  b2psifn = temp
#else
  b2psifn = 0.
#endif
end function b2psifn


! B2bfn
! =====
function b2bfn(e,f,g,h)
  use basic
  use m3dc1_nint

  vectype, dimension(dofs_per_element) :: b2bfn
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0) then
     b2bfn = 0.
     return
  end if
  
  if(surface_int) then
     if(inocurrent_norm.eq.1 .and. imulti_region.eq.0) then
        temp = 0.
     else 
        temp79a = f(:,OP_1)*h(:,OP_1)*ni79(:,OP_1)**2
        temp = intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,2),g(:,OP_DRP)) &
             - intx5(e(:,:,OP_1),ri3_79,temp79a,norm79(:,1),g(:,OP_DZP))
     endif
  else
     temp79a = ri3_79*ni79(:,OP_1)**2
     temp = - &
          (intx5(e(:,:,OP_DZ),temp79a,f(:,OP_1),g(:,OP_DRP),h(:,OP_1)) &
          -intx5(e(:,:,OP_DR),temp79a,f(:,OP_1),g(:,OP_DZP),h(:,OP_1)))
  end if
  
  b2bfn = temp
#else
  b2bfn = 0.
#endif
end function b2bfn

function b2phidot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2phidot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = intx2(e(:,:,OP_DR),f(:,OP_DR))   &
          + intx2(e(:,:,OP_DZ),f(:,OP_DZ))
  end if

  b2phidot = temp
end function b2phidot

function b2chidot(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2chidot
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp
  temp = 0.

  if(surface_int) then
     temp = 0.
  else
     if(itor.eq.1) temp = 2.*intx3(e(:,:,OP_1),ri4_79,f(:,OP_DZ))
  end if

  b2chidot = temp
end function b2chidot

function b2psi2bfpe(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b2psi2bfpe
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j
  vectype, dimension(dofs_per_element) :: temp
  temp = 0.

  temp79a = ri_79*(j(:,OP_DZ)*f(:,OP_DR) - j(:,OP_DR)*f(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = temp79a - (j(:,OP_DR)*i(:,OP_DR) + j(:,OP_DZ)*i(:,OP_DZ))   &
                    + ri2_79*h(:,OP_1)*j(:,OP_DP)
#endif
  temp79a = temp79a*b2i79(:,OP_1)*ni79(:,OP_1)


  if(surface_int) then
        temp = 0.
  else
        temp = intx4(e(:,:,OP_DR),ri2_79,temp79a,g(:,OP_DR))    &
             + intx4(e(:,:,OP_DZ),ri2_79,temp79a,g(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp + intx4(e(:,:,OP_DZ),ri_79,temp79a,i(:,OP_DR))   &
                    - intx4(e(:,:,OP_DR),ri_79,temp79a,i(:,OP_DZ))
#endif

  end if

  b2psi2bfpe = temp
end function b2psi2bfpe
!=============================================================================
! B3 TERMS
!=============================================================================

! B3pe
! ====
!function b3pe(e,f)
  !use basic
  !use m3dc1_nint

  !implicit none

  !vectype, dimension(dofs_per_element) :: b3pe
  !vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  !vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  !vectype, dimension(dofs_per_element) :: temp

  !if(surface_int) then
     !temp = 0.
  !else
     !temp = intx2(e(:,:,OP_1),f(:,OP_1))
  !end if

  !b3pe = temp
!end function b3pe

function b3pe
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3pe
  type(prodarray) :: temp

  if(surface_int) then
     temp%len = 0
  else
     temp79a=1.
     temp=prod(temp79a,OP_1,OP_1)
  end if

  b3pe = temp
end function b3pe

  function kappat_p1
    use basic
    use m3dc1_nint

    implicit none

    type(prodarray) :: kappat_p1
    type(prodarray) :: temp

    if(gam.le.1) then 
       kappat_p1%len = 0
       return
    end if

    if(surface_int) then
       temp%len = 0
    else
       temp = prod(-ni79(:,OP_1)*kap79(:,OP_1),OP_DR,OP_DR) &
            + prod(-ni79(:,OP_1)*kap79(:,OP_1),OP_DZ,OP_DZ)

#if defined(USE3D) || defined(USECOMPLEX)
       temp79a = kap79(:,OP_1)
       !if(iupstream.eq.1) then    
          !temp79a = temp79a + abs(vzt79(:,OP_1))*magus
       !endif
       temp = temp +                       &
            prod(-ri2_79*ni79(:,OP_1)*temp79a,OP_DP,OP_DP)
       !if(iupstream.eq.2) then
          !!temp79a = abs(vzt79(:,OP_1))*magus2
          !temp79a = magus2
          !temp = temp -                    &
               !intx5(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),ni79(:,OP_1),temp79a)
       !endif
#endif
       if(hypp.ne.0.) then
          ! Laplacian[f g]

          temp79b = hypp
          temp = temp + prod(-temp79b,OP_LP,OP_LP)
          !temp = temp + &
                   !( prod( hypp*ri4_79,OP_DZZ,OP_DZZPP) &
                   !+ prod( hypp*ri4_79,OP_DRR,OP_DZZPP) &
                   !+ prod( hypp*ri4_79,OP_DZZ,OP_DRRPP) &
          !         + prod( hypp*ri4_79,OP_DRR,OP_DRRPP))
        endif

    end if
    
    temp79b = gam - 1.
    kappat_p1 = temp*temp79b
  end function kappat_p1

  function kappat_lin_p1
    use basic
    use m3dc1_nint

    implicit none

    type(prodarray) :: kappat_lin_p1
    type(prodarray) :: temp

    if(gam.le.1) then 
       kappat_lin_p1%len = 0
       return
    end if

    if(surface_int) then
       temp%len = 0
    else
       temp79a = ni79(:,OP_1)**2

       temp =  prod(-ni79(:,OP_1)*kap79(:,OP_1),OP_DR,OP_DR) &
            +  prod(-ni79(:,OP_1)*kap79(:,OP_1),OP_DZ,OP_DZ) &
            +  prod(n079(:,OP_DR)*temp79a*kap79(:,OP_1),OP_DR,OP_1) &
            +  prod(n079(:,OP_DZ)*temp79a*kap79(:,OP_1),OP_DZ,OP_1)


#if defined(USECOMPLEX)
       temp = temp + &
            prod(ri2_79*ni79(:,OP_1)*kap79(:,OP_1),OP_1,OP_DPP)
#endif
#if defined(USE3D) || defined(USECOMPLEX)
       !temp79a = kap79(:,OP_1)
       !if(iupstream.eq.1) then    
          !temp79a = temp79a + abs(vzt79(:,OP_1))*magus
       !endif
       !temp = temp +                       &
       !     prod(-ri2_79*ni79(:,OP_1)*temp79a,OP_DP,OP_DP)
       !if(iupstream.eq.2) then
          !!temp79a = abs(vzt79(:,OP_1))*magus2
          !temp79a = magus2
          !temp = temp -                    &
               !intx5(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),ni79(:,OP_1),temp79a)
       !endif
#endif
       if(hypp.ne.0.) then
        ! Laplacian[f g]

          temp79b = hypp
          temp = temp + prod(-temp79b,OP_LP,OP_LP)
       endif

     end if
    
    temp79b = (gam-1.)
    kappat_lin_p1 = temp * temp79b
  end function kappat_lin_p1

function kappat_lin_pn2(f)
    use basic
    use m3dc1_nint

    implicit none

    type(prodarray) :: kappat_lin_pn2
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    type(prodarray) :: temp

    if(gam.le.1) then
       kappat_lin_pn2%len = 0
       return
    end if

    temp79a = ni79(:,OP_1)**2*kap79(:,OP_1)
    temp79b = temp79a*ni79(:,OP_1)

    if(surface_int) then
       temp%len = 0
    else
       temp =  prod(f(:,OP_DR)*temp79a,OP_DR,OP_1) &
            +  prod(f(:,OP_DZ)*temp79a,OP_DZ,OP_1) &
            +  prod(f(:,OP_1)*temp79a,OP_DR,OP_DR) &
            +  prod(f(:,OP_1)*temp79a,OP_DZ,OP_DZ) &
            +  prod(-2.*f(:,OP_1)*n079(:,OP_DR)*temp79b,OP_DR,OP_1) &
            +  prod(-2.*f(:,OP_1)*n079(:,OP_DZ)*temp79b,OP_DZ,OP_1)

#if defined(USECOMPLEX)
       temp = temp + &
            prod(-ri2_79*f(:,OP_1)*temp79a,OP_1,OP_DPP)
#endif
    end if
    
    temp79b = (gam-1.)
    kappat_lin_pn2 = temp * temp79b
  end function kappat_lin_pn2

! B3pe27
! ======
function b3pe27(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pe27
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp
  real, dimension(MAX_PTS) :: r

  if(surface_int) then
     temp = 0.
  else
     r = sqrt((x_79-xmag)**2 + (z_79-zmag)**2)
     temp79b = 1. + tanh((r-libetap)/p1)

     temp = intx3(e(:,:,OP_1),f(:,OP_1),temp79b)
  end if

  b3pe27 = temp
end function b3pe27


! B3q
! ===
function b3q(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3q
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = intx2(e(:,:,OP_1),f(:,OP_1))
  end if

  b3q = temp
end function b3q

! B3psipsieta
! ===========
function b3psipsieta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3psipsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else
     temp = (gam-1.)* &
           intx5(e(:,:,OP_1),ri2_79,f(:,OP_GS), g(:,OP_GS), h(:,OP_1))   
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + (gam-1)*   &
           (intx5(e(:,:,OP_1),ri4_79,f(:,OP_DRP),g(:,OP_DRP),h(:,OP_1))   &
         +  intx5(e(:,:,OP_1),ri4_79,f(:,OP_DZP),g(:,OP_DZP),h(:,OP_1)))
     if(irunaway .gt. 2) then
        temp = temp + 1.*(gam-1.) * &
              (-intx6(e(:,:,OP_1),ri3_79,f(:,OP_DZ),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
              + intx6(e(:,:,OP_1),ri3_79,f(:,OP_DR),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)))
     endif
#endif
  end if

  b3psipsieta = temp
end function b3psipsieta

function b3psipsieta1(g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3psipsieta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else
     temp = prod((gam-1.)*ri2_79*g(:,OP_GS)*h(:,OP_1),OP_1,OP_GS)   
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + &
           (prod((gam-1.)*ri4_79*g(:,OP_DRP)*h(:,OP_1),OP_1,OP_DRP)   &
         +  prod((gam-1.)*ri4_79*g(:,OP_DZP)*h(:,OP_1),OP_1,OP_DZP))
     if(irunaway .gt. 2) then
        !temp = temp + 1.*(gam-1.) * &
              !(-prod(e(:,:,OP_1),ri3_79,f(:,OP_DZ),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
        !      + prod(e(:,:,OP_1),ri3_79,f(:,OP_DR),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)))
     endif
#endif
  end if

  b3psipsieta1 = temp
end function b3psipsieta1

! B3psibeta
! ===========
function b3psibeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3psibeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (intx5(e(:,:,OP_1),ri3_79,f(:,OP_DZP),g(:,OP_DR),h(:,OP_1))  &
          -intx5(e(:,:,OP_1),ri3_79,f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1)))
#endif
  end if

  if(irunaway .gt. 2) then
     temp = temp + 1.*(gam-1.) * &
                   (-intx6(e(:,:,OP_1),ri2_79,f(:,OP_GS),g(:,OP_1),h(:,OP_1),i(:,OP_1)) &
                  + intx6(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
                  + intx6(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
  end if

  b3psibeta = temp
end function b3psibeta

function b3psibeta1(g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3psibeta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else
     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     temp = &
          (prod( 2.*(gam-1.)*ri3_79*g(:,OP_DR)*h(:,OP_1),OP_1,OP_DZP)  &
          +prod(-2.*(gam-1.)*ri3_79*g(:,OP_DZ)*h(:,OP_1),OP_1,OP_DRP))
#endif
  end if

  if(irunaway .gt. 2) then
     !temp = temp + 1.*(gam-1.) * &
                   !(-prod(e(:,:,OP_1),ri2_79,f(:,OP_GS),g(:,OP_1),h(:,OP_1),i(:,OP_1)) &
                  !+ prod(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
     !             + prod(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
  end if

  b3psibeta1 = temp
end function b3psibeta1

function b3psibeta2(f,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3psibeta2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else
     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (prod( ri3_79*f(:,OP_DZP)*h(:,OP_1),OP_1,OP_DR)  &
          +prod(-ri3_79*f(:,OP_DRP)*h(:,OP_1),OP_1,OP_DZ))
#endif
  end if

  ! if(irunaway .gt. 0) then
  !    temp = temp + 1.*(gam-1.) * &
  !                  (prod(-ri2_79*f(:,OP_GS)*h(:,OP_1)*i(:,OP_1),OP_1,OP_1) &
  !                 + prod(ri2_79*f(:,OP_DR)*h(:,OP_1)*i(:,OP_1),OP_1,OP_DR) &
  !                 + prod(ri2_79*f(:,OP_DZ)*h(:,OP_1)*i(:,OP_1),OP_1,OP_DZ))
  ! end if

  b3psibeta2 = temp
end function b3psibeta2

! B3psifeta
! =========
function b3psifeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3psifeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (intx5(e(:,:,OP_1),ri3_79,f(:,OP_DZP),g(:,OP_DRP),h(:,OP_1))  &
          -intx5(e(:,:,OP_1),ri3_79,f(:,OP_DRP),g(:,OP_DZP),h(:,OP_1)))
      if(irunaway .gt. 2) then
         temp = temp + 1.*(gam-1.) * &
                       (intx6(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
                      + intx6(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)) &
                      - intx6(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
                      - intx6(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
      end if 
#endif
  end if

  b3psifeta = temp
end function b3psifeta

function b3psifeta1(g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3psifeta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else
     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (prod( ri3_79*g(:,OP_DRP)*h(:,OP_1),OP_1,OP_DZP)  &
          +prod(-ri3_79*g(:,OP_DZP)*h(:,OP_1),OP_1,OP_DRP))
      ! if(irunaway .gt. 0) then
      !    temp = temp + 1.*(gam-1.) * &
      !                  (prod(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
      !                 + prod(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)) &
      !                 - prod(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
      !                 - prod(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
      ! end if 
#endif
  end if

  b3psifeta1 = temp
end function b3psifeta1

function b3psifeta2(f,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3psifeta2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else
     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (prod( ri3_79*f(:,OP_DZP)*h(:,OP_1),OP_1,OP_DRP)  &
          +prod(-ri3_79*f(:,OP_DRP)*h(:,OP_1),OP_1,OP_DZP))
      ! if(irunaway .gt. 0) then
      !    temp = temp + 1.*(gam-1.) * &
      !                  (prod(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1),i(:,OP_1)) &
      !                 + prod(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1),i(:,OP_1)) &
      !                 - prod(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
      !                 - prod(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
      ! end if 
#endif
  end if

  b3psifeta2 = temp
end function b3psifeta2

! B3bbeta
! =======
function b3bbeta(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3bbeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else 
     temp = (gam-1.)* &
          (intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          +intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))
  end if

  b3bbeta = temp
end function b3bbeta

function b3bbeta1(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3bbeta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else 
     temp = &
          (prod((gam-1.)*ri2_79*g(:,OP_DZ)*h(:,OP_1),OP_1,OP_DZ) &
          +prod((gam-1.)*ri2_79*g(:,OP_DR)*h(:,OP_1),OP_1,OP_DR))
  end if

  b3bbeta1 = temp
end function b3bbeta1

! B3bfeta
! =======
function b3bfeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3bfeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else 
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
          +intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)))
     if(irunaway .gt. 2) then
        temp = temp + 1.*(gam-1.)* &
                      (intx6(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
                     - intx6(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     endif
#endif
  endif

  b3bfeta = temp
end function b3bfeta

function b3bfeta1(g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3bfeta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else 
     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (prod(ri2_79*g(:,OP_DZP)*h(:,OP_1),OP_1,OP_DZ) &
          +prod(ri2_79*g(:,OP_DRP)*h(:,OP_1),OP_1,OP_DR))
     ! if(irunaway .gt. 0) then
     !    temp = temp + 1.*(gam-1.)* &
     !                  (prod(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
     !                 - prod(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     ! endif
#endif
  endif

  b3bfeta1 = temp
end function b3bfeta1

function b3bfeta2(f,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3bfeta2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else 
     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     temp = &
          (prod((gam-1.)*ri2_79*f(:,OP_DZ)*h(:,OP_1),OP_1,OP_DZP) &
          +prod((gam-1.)*ri2_79*f(:,OP_DR)*h(:,OP_1),OP_1,OP_DRP))
     if(irunaway .gt. 0) then
        !temp = temp + 1.*(gam-1.)* &
                      !(prod(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
        !             - prod(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     endif
#endif
  endif

  b3bfeta2 = temp
end function b3bfeta2

! B3ffeta
! =======
function b3ffeta(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3ffeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else 
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZP),g(:,OP_DZP),h(:,OP_1)) &
          +intx5(e(:,:,OP_1),ri2_79,f(:,OP_DRP),g(:,OP_DRP),h(:,OP_1)))
     if(irunaway .gt. 0) then
        temp = temp + 1.*(gam-1.)* &
               (intx6(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
              - intx6(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     endif
#endif
  end if

  b3ffeta = temp
end function b3ffeta

function b3ffeta1(g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3ffeta1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  
  type(prodarray) :: temp

  if(gam.le.1. .or. surface_int) then
     temp%len = 0
  else 
     temp%len = 0
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (prod(ri2_79*g(:,OP_DZP)*h(:,OP_1),OP_1,OP_DZP) &
          +prod(ri2_79*g(:,OP_DRP)*h(:,OP_1),OP_1,OP_DRP))
     ! if(irunaway .gt. 0) then
     !    temp = temp + 1.*(gam-1.)* &
     !           (prod(e(:,:,OP_1),ri_79,f(:,OP_DZP),g(:,OP_DR),h(:,OP_1),i(:,OP_1)) &
     !          - prod(e(:,:,OP_1),ri_79,f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1),i(:,OP_1)))
     ! endif
#endif
  end if

  b3ffeta1 = temp
end function b3ffeta1

! B3pepsid
! ========
function b3pepsid(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pepsid
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0 .or. surface_int) then
     b3pepsid = 0.
     return
  end if

  temp = intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
       + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
       - intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_GS),h(:,OP_1)) &
       + gam* &
       (intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_DZ)) &
       +intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_DR)) &
       -intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_GS),h(:,OP_DP)))

  b3pepsid = temp
#else
  b3pepsid = 0.
#endif
end function b3pepsid


! B3pebd
! ======
function b3pebd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pebd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0 .or. surface_int) then
     b3pebd = 0.
     return
  end if

  temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
        -intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
       + gam* &
       (intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DR),h(:,OP_DZ)) &
       -intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DZ),h(:,OP_DR)))

  b3pebd = temp
end function b3pebd


! B3pefd
! ======
function b3pefd(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pefd
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(itwofluid.eq.0 .or. surface_int) then
     b3pefd = 0.
     return
  end if

  temp = intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DRP),h(:,OP_1)) &
        -intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZP),h(:,OP_1)) &
       + gam* &
       (intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DRP),h(:,OP_DZ)) &
       -intx5(e(:,:,OP_1),ri_79,f(:,OP_1),g(:,OP_DZP),h(:,OP_DR)))

  b3pefd = temp
#else
  b3pefd = 0.
#endif
end function b3pefd



! B3pedkappa
! ==========
function b3pedkappa(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pedkappa
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3pedkappa = 0.
     return
  end if

  if(surface_int) then
     temp = intx5(e(:,:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1 ),h(:,OP_1)) &
          + intx5(e(:,:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1)) &
          + intx5(e(:,:,OP_1),norm79(:,1),f(:,OP_1 ),g(:,OP_DR),h(:,OP_1)) &
          + intx5(e(:,:,OP_1),norm79(:,2),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1))
  else
     temp = &
          - intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1)) &
          - intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1 ),h(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1)) &
          - intx4(e(:,:,OP_DR),f(:,OP_1 ),g(:,OP_DR),h(:,OP_1))
  
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = h(:,OP_1)
     if(iupstream.eq.1) then    
        temp79a = temp79a + abs(i(:,OP_1))*magus
     endif
     temp = temp +                       &
          intx5(e(:,:,OP_1),ri2_79,f(:,OP_DPP),g(:,OP_1),temp79a)
     if(iupstream.eq.2) then
       temp79a = abs(i(:,OP_1))*magus
     temp = temp -                       &
          intx5(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),g(:,OP_1),temp79a)
     endif
#endif
     if(hypp.ne.0.) then
        ! Laplacian[f g]
        temp79a = f(:,OP_LP)*g(:,OP_1) + f(:,OP_1)*g(:,OP_LP) &
             + 2.*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR))

        if(ihypkappa.eq.1) then        
           temp = temp - hypp* &
                (intx3(e(:,:,OP_LP),temp79a,h(:,OP_1 )) &
                +intx3(e(:,:,OP_DZ),temp79a,h(:,OP_DZ)) &
                +intx3(e(:,:,OP_DR),temp79a,h(:,OP_DR)))
        else
           temp = temp - hypp*intx2(e(:,:,OP_LP),temp79a)
        endif
     endif
  end if

  b3pedkappa = (gam-1.)*temp  
end function b3pedkappa


! B3tekappa
! ==========
function b3tekappa(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3tekappa
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3tekappa = 0.
     return
  end if

  if(surface_int) then
     temp = intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
          + intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1))
  else
     temp = &
          - intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
          - intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1))
  
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = g(:,OP_1)
     if(iupstream.eq.1) then    
        temp79a = temp79a + abs(h(:,OP_1))*magus
     endif
     temp = temp +                       &
          intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a)
     if(iupstream.eq.2) then    
        temp79a = abs(h(:,OP_1))*magus
        temp = temp -                    &
          intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),temp79a)
     endif
#endif
     if(hypp.ne.0.) then

        if(ihypkappa.eq.1) then        
           temp = temp - hypp* &
                (intx3(e(:,:,OP_LP),f(:,OP_LP),g(:,OP_1 )) &
                +intx3(e(:,:,OP_DZ),f(:,OP_LP),g(:,OP_DZ)) &
                +intx3(e(:,:,OP_DR),f(:,OP_LP),g(:,OP_DR)))
        else
           temp = temp - hypp*intx2(e(:,:,OP_LP),f(:,OP_LP))
        endif
     endif
  end if

  b3tekappa = (gam-1.)*temp
end function b3tekappa

function b3tekappa1(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: b3tekappa1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(gam.le.1.) then
     b3tekappa1%len = 0
     return
  end if

  if(surface_int) then
     temp = prod(norm79(:,1)*g(:,OP_1),OP_1,OP_DR) &
          + prod(norm79(:,2)*g(:,OP_1),OP_1,OP_DZ)
  else
     temp = &
            prod(-g(:,OP_1),OP_DZ,OP_DZ) &
          + prod(-g(:,OP_1),OP_DR,OP_DR)
  
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = g(:,OP_1)
     if(iupstream.eq.1) then    
        temp79a = temp79a + abs(h(:,OP_1))*magus
     endif
     temp = temp +                       &
          prod(ri2_79*temp79a,OP_1,OP_DPP)
     !if(iupstream.eq.2) then    
        !temp79a = abs(h(:,OP_1))*magus
        !temp = temp -                    &
          !prod(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),temp79a)
     !endif
#endif
     if(hypp.ne.0.) then

        !if(ihypkappa.eq.1) then        
           !temp = temp - hypp* &
                !(prod(e(:,:,OP_LP),f(:,OP_LP),g(:,OP_1 )) &
                !+prod(e(:,:,OP_DZ),f(:,OP_LP),g(:,OP_DZ)) &
                !+prod(e(:,:,OP_DR),f(:,OP_LP),g(:,OP_DR)))
        !else
           !temp = temp - hypp*prod(e(:,:,OP_LP),f(:,OP_LP))
        !endif
     endif
  end if

  b3tekappa1 = (gam-1.)*temp
end function b3tekappa1

! B3pedkappag
! ===========
function b3pppkappag(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pppkappag
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, intent(in), dimension(MAX_PTS) :: i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3pppkappag = 0.
     return
  end if

#ifdef USECOMPLEX
  temp79a = h(:,OP_DR)*conjg(g(:,OP_DR)) + h(:,OP_DZ)*conjg(g(:,OP_DZ)) &
       + h(:,OP_DP)*conjg(g(:,OP_DP))*ri2_79
#else
  temp79a = h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)
#endif

#ifdef USE3D
  temp79a = temp79a + h(:,OP_DP)*g(:,OP_DP)*ri2_79
#endif

  if(surface_int) then
     temp = intx5(e(:,:,OP_1),norm79(:,1),f(:,OP_DR),temp79a,i) &
          + intx5(e(:,:,OP_1),norm79(:,2),f(:,OP_DZ),temp79a,i)
  else
     temp = &
          - intx4(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,i) &
          - intx4(e(:,:,OP_DR),f(:,OP_DR),temp79a,i)
  
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + intx5(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a,i)
#endif
  end if

  b3pppkappag = (gam-1.)*kappag*temp
end function b3pppkappag

! B3pkappag
! =========
function b3pkappag(e,f,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: b3pkappag
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, intent(in), dimension(MAX_PTS) :: i
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     b3pkappag = 0.
     return
  end if

  if(surface_int) then
     temp = intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DR),i) &
          + intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DZ),i)
  else
     temp = &
          - intx3(e(:,:,OP_DZ),f(:,OP_DZ),i) &
          - intx3(e(:,:,OP_DR),f(:,OP_DR),i)
  
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),i)
#endif
  end if

  b3pkappag = -gradp_crit**2*(gam-1.)*kappag*temp
end function b3pkappag



!============================================================================
! N1 TERMS
!============================================================================

! N1n
! ===
function n1n(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element,dofs_per_element) :: n1n
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     n1n = 0.
  else
     n1n = intxx2(e(:,:,OP_1),f(:,:,OP_1))
  end if
end function n1n


! N1ndenm
! =======
function n1ndenm(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1ndenm
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     if(inograd_n.eq.1) then
        temp = 0.
     else
        temp =  &
             (intx4(e(:,:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
             +intx4(e(:,:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1)))
     end if
  else
     temp = -(intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
          +   intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_1)))

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = g(:,OP_1)
     if(iupstream .eq. 1) then   
        temp79a = temp79a+abs(h(:,OP_1))*magus
     endif
     temp = temp + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DPP),temp79a) &
          + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_DP))
     if(iupstream .eq. 2) then   
        temp79a = abs(h(:,OP_1))*magus
        temp = temp - intx4(e(:,:,OP_DPP),ri4_79,f(:,OP_DPP),temp79a)
     endif
#endif

     if(hypp.ne.0.) then
        if(ihypkappa.eq.1) then
           temp = temp - hypp*intx3(e(:,:,OP_LP),f(:,OP_LP),g(:,OP_1))
        else
           temp = temp - hypp*intx2(e(:,:,OP_LP),f(:,OP_LP))
        endif
     endif
  end if

  n1ndenm = temp
end function n1ndenm


! N1nu
! ====
function n1nu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1nu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
                - intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,2),g(:,OP_DR))
        endif
     else
        temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
             - intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ))
     endif

  n1nu = temp
end function n1nu


! N1nv
! ====
function n1nv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1nv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp = 0.
     else
        temp = -intx3(e(:,:,OP_1),f(:,OP_1 ),g(:,OP_DP)) &
             -  intx3(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1 ))
     end if

  n1nv = temp
#else
  n1nv = 0.
#endif
end function n1nv

! N1nchi
! ======
function n1nchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1nchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = &
                - intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
                - intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        endif
     else
        temp = intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZ)) &
             + intx4(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DR))
     end if

  n1nchi = temp
end function n1nchi


! N1s
! ===
function n1s(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: n1s
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  if(surface_int) then
     n1s = 0.
  else
     n1s = intx2(e(:,:,OP_1),f(:,OP_1))
  end if
end function n1s

function t3tndenm(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tndenm
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     if(inograd_n.eq.1) then
        temp = 0.
     else
        temp = - &
             (intx5(e(:,:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DR),h(:,OP_1)) &
             +intx5(e(:,:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DZ),h(:,OP_1)))
     end if
  else
     temp =  &
          (intx4(e(:,:,OP_DZ),f(:,OP_1),g(:,OP_DZ),h(:,OP_1)) &
          +intx4(e(:,:,OP_DR),f(:,OP_1),g(:,OP_DR),h(:,OP_1)) &
          +intx4(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          +intx4(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))

#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp - &
          (   intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DPP),h(:,OP_1 )) &
          +   intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DP ),h(:,OP_DP)))
#endif

     if(hypp.ne.0.) then
        if(ihypkappa.eq.1) then
           temp = temp + hypp* &
                ( intx4(e(:,:,OP_LP),f(:,OP_1),g(:,OP_LP),h(:,OP_1))   &
                + intx4(e(:,:,OP_1),f(:,OP_LP),g(:,OP_LP),h(:,OP_1))   &
                + 2.*intx4(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_LP),h(:,OP_1))&
                + 2.*intx4(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_LP),h(:,OP_1)))
        else
           temp = temp + hypp* &
                (intx3(e(:,:,OP_LP),f(:,OP_1),g(:,OP_LP))   &
                + intx3(e(:,:,OP_1),f(:,OP_LP),g(:,OP_LP))   &
                + 2.*intx3(e(:,:,OP_DR),f(:,OP_DR),g(:,OP_LP)) &
                + 2.*intx3(e(:,:,OP_DZ),f(:,OP_DZ),g(:,OP_LP)))
        endif
     endif
  end if

  t3tndenm = temp
end function t3tndenm

function t3tndenm1(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tndenm1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

  if(surface_int) then
     if(inograd_n.eq.1) then
        temp%len = 0
     else
        temp = &
             (prod(-norm79(:,1)*g(:,OP_DR)*h(:,OP_1),OP_1,OP_1) &
             +prod(-norm79(:,2)*g(:,OP_DZ)*h(:,OP_1),OP_1,OP_1))
     end if
  else
     temp =  &
          (prod(g(:,OP_DZ)*h(:,OP_1),OP_DZ,OP_1) &
          +prod(g(:,OP_DR)*h(:,OP_1),OP_DR,OP_1) &
          +prod(g(:,OP_DZ)*h(:,OP_1),OP_1,OP_DZ) &
          +prod(g(:,OP_DR)*h(:,OP_1),OP_1,OP_DR))

#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + &
          (   prod(-ri2_79*g(:,OP_DPP)*h(:,OP_1 ),OP_1,OP_1) &
          +   prod(-ri2_79*g(:,OP_DP )*h(:,OP_DP),OP_1,OP_1))
#endif

     if(hypp.ne.0.) then
        if(ihypkappa.eq.1) then
           temp = temp + hypp* &
                ( prod(g(:,OP_LP)*h(:,OP_1),OP_LP,OP_1)   &
                + prod(g(:,OP_LP)*h(:,OP_1),OP_1,OP_LP)   &
                + prod(2.*g(:,OP_LP)*h(:,OP_1),OP_DR,OP_DR)&
                + prod(2.*g(:,OP_LP)*h(:,OP_1),OP_DZ,OP_DZ))
        else
           temp = temp + hypp* &
                (prod(g(:,OP_LP),OP_LP,OP_1)   &
                + prod(g(:,OP_LP),OP_1,OP_LP)   &
                + prod(2.*g(:,OP_LP),OP_DR,OP_DR) &
                + prod(2.*g(:,OP_LP),OP_DZ,OP_DZ))
        endif
     endif
  end if

  t3tndenm1 = temp
end function t3tndenm1

function t3ts(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3ts
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = -intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_1))
  end if

  t3ts = temp
end function t3ts

function t3ts1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3ts1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(surface_int) then
     temp%len = 0
  else
     temp = prod(-g(:,OP_1),OP_1,OP_1)
  end if

  t3ts1 = temp
end function t3ts1

!============================================================================
! NRE1 TERMS
!============================================================================

! NRE1nrediff
! =======
function nre1nrediff(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nrediff
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
#if defined(USE3D)
  ! B . Grad(p)
  ! [ p, psi] / R + 1/R^2 F p' - <f', p>
  temp79a = ri_79*(f(:,OP_DZ)*b2i79(:,OP_1)*pstx79(:,OP_DR) - &
            f(:,OP_DR)*b2i79(:,OP_1)*pstx79(:,OP_DZ))
  temp79a = temp79a+0.5*ri_79*(f(:,OP_1)*b2i79(:,OP_DZ)*pstx79(:,OP_DR) - &
            f(:,OP_1)*b2i79(:,OP_DR)*pstx79(:,OP_DZ))
  temp79a = temp79a &
          + ri2_79*bztx79(:,OP_1)*b2i79(:,OP_1)*f(:,OP_DP) &
          - b2i79(:,OP_1)*(bfptx79(:,OP_DZ)*f(:,OP_DZ) & 
          + bfptx79(:,OP_DR)*f(:,OP_DR))
  if(surface_int) then
      temp = 0.
  else
      temp = intx4(e(:,:,OP_DZ),ri_79,-temp79a,pstx79(:,OP_DR)) &
           - intx4(e(:,:,OP_DR),ri_79,-temp79a,pstx79(:,OP_DZ))
      temp = temp &
         + intx4(e(:,:,OP_DP),ri2_79,-temp79a,bztx79(:,OP_1)) &
         - intx3(e(:,:,OP_DZ),-temp79a,bfptx79(:,OP_DZ)) &
         - intx3(e(:,:,OP_DR),-temp79a,bfptx79(:,OP_DR))
  end if
#else
  temp79a = ri_79*(f(:,OP_DZ)*b2i79(:,OP_1)*pstx79(:,OP_DR) - &
            f(:,OP_DR)*b2i79(:,OP_1)*pstx79(:,OP_DZ))
  temp79a = temp79a+0.5*ri_79*(f(:,OP_1)*b2i79(:,OP_DZ)*pstx79(:,OP_DR) - &
            f(:,OP_1)*b2i79(:,OP_DR)*pstx79(:,OP_DZ))

  if(surface_int) then
      temp = 0.
  else
      temp = intx4(e(:,:,OP_DZ),ri_79,-temp79a,pstx79(:,OP_DR)) &
           - intx4(e(:,:,OP_DR),ri_79,-temp79a,pstx79(:,OP_DZ))
  end if
#endif

  nre1nrediff = temp
end function nre1nrediff

! NRE1nrepsi
! ====
function nre1nrepsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nrepsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = intx5(e(:,:,OP_DZ),ri_79,f(:,OP_1),g(:,OP_1),h(:,OP_DR)) &
        - intx5(e(:,:,OP_DR),ri_79,f(:,OP_1),g(:,OP_1),h(:,OP_DZ))
     temp = temp * 1.000 * cre 
  endif

  nre1nrepsi = temp
end function nre1nrepsi

! NRE1nreb
! ====
function nre1nreb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nreb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  temp = 0
#if defined(USE3D)
  temp = -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_1),h(:,OP_1)) - &
          intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_DP),h(:,OP_1)) - &
          intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_1),h(:,OP_DP))
#elif defined(USECOMPLEX)
  temp = -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DP),g(:,OP_1),h(:,OP_1))
#endif
  temp = temp * 1.000 * cre

  nre1nreb = temp
end function nre1nreb

! NRE1nreu
! ====
function nre1nreu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: nre1nreu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
                - intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,2),g(:,OP_DR))
        endif
     else
        temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
             - intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ))
     endif

  nre1nreu = temp
end function nre1nreu

!============================================================================
! P1 TERMS
!============================================================================

! P1pu
! ====
function p1pu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1pu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then 
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
                - intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,2),g(:,OP_DR))
        end if
     else
        temp = intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
             - intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ))

        if(itor.eq.1) then
           temp = temp + &
                2.*(gam-1.)*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DZ))
        endif
     end if

  p1pu = temp
end function p1pu

function p1pu1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1pu1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then 
           temp%len = 0
          else
           !temp = intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
           !     - intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,2),g(:,OP_DR))
        end if
     else
        temp =   prod( r_79*g(:,OP_DR),OP_DZ,OP_1) &
               + prod(-r_79*g(:,OP_DZ),OP_DR,OP_1)

        if(itor.eq.1) then
           temp = temp + &
                prod(2.*(gam-1.)*g(:,OP_DZ),OP_1,OP_1)
        endif
     end if

  p1pu1 = temp
end function p1pu1

function p1pu2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1pu2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  type(prodarray) :: temp

     if(surface_int) then
           temp%len = 0
     else
        temp =   prod( r_79*f(:,OP_1),OP_DZ,OP_DR) &
               + prod(-r_79*f(:,OP_1),OP_DR,OP_DZ)

        if(itor.eq.1) then
           temp = temp + &
                prod(2.*(gam-1.)*f(:,OP_1),OP_1,OP_DZ)
         endif
     end if

  p1pu2 = temp
end function p1pu2


! P1pv
! ====
function p1pv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1pv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp = 0.
     else
        temp = - intx3(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1)) &
             - gam*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DP))
     endif
#else
  temp = 0.
#endif

  p1pv = temp
end function p1pv

function p1pv1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1pv1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp%len = 0
     else
        temp = prod(-g(:,OP_1),OP_1,OP_DP) &
             + prod(-gam*g(:,OP_DP),OP_1,OP_1)
     endif
#else
  temp%len = 0
#endif

  p1pv1 = temp
end function p1pv1

function p1pv2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1pv2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp%len = 0
     else
        temp = prod(-f(:,OP_DP),OP_1,OP_1) &
             + prod(-gam*f(:,OP_1),OP_1,OP_DP)
      endif
#else
  temp%len = 0
#endif

  p1pv2 = temp
end function p1pv2


! P1pchi
! ======
function p1pchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1pchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = &
                - intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
                - intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        end if
     else
        temp = intx4(e(:,:,OP_DR),ri2_79,g(:,OP_DR),f(:,OP_1)) &
             + intx4(e(:,:,OP_DZ),ri2_79,g(:,OP_DZ),f(:,OP_1)) &
             - (gam-1.)*intx4(e(:,:,OP_1),ri2_79,g(:,OP_GS),f(:,OP_1))
     endif

  p1pchi = temp
end function p1pchi

function p1pchi1(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1pchi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g

  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
        !   temp = &
                !- intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
        !        - intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        end if
     else
          temp = prod(ri2_79*g(:,OP_DR),OP_DR,OP_1) &
               + prod(ri2_79*g(:,OP_DZ),OP_DZ,OP_1) &
               + prod(-(gam-1.)*ri2_79*g(:,OP_GS),OP_1,OP_1)
     endif

  p1pchi1 = temp
end function p1pchi1

function p1pchi2(f)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1pchi2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           temp = &
                  prod(-ri2_79*f(:,OP_1)*norm79(:,1),OP_1,OP_DR) &
                + prod(-ri2_79*f(:,OP_1)*norm79(:,2),OP_1,OP_DZ)
         end if
     else
        temp =   prod(ri2_79*f(:,OP_1),OP_DR,OP_DR) &
               + prod(ri2_79*f(:,OP_1),OP_DZ,OP_DZ) &
               + prod(-(gam-1.)*ri2_79*f(:,OP_1),OP_1,OP_GS)
     endif

  p1pchi2 = temp
end function p1pchi2


! P1psipsikappar
! ==============
function p1psipsikappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psipsikappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psipsikappar = 0.
     return
  end if

  if(surface_int) then
!!$     ! does better without this term
!!$     ! justification: assert natural b.c. B.grad(T) = 0
!!$     temp = 0.
     temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))
     temp = intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1 )) &
          - intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1 )) &
          + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_1 ),i(:,OP_DR)) &
          - intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_1 ),i(:,OP_DZ))
  else
     temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)*i(:,OP_1)
     temp79b = ri2_79*k(:,OP_1)*j(:,OP_1)*h(:,OP_1)

     temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DZ),i(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DZ),i(:,OP_DR)) &
          - intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DR),i(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DR),i(:,OP_DZ))
  end if

  p1psipsikappar = (gam - 1.) * temp
end function p1psipsikappar

! P1psipsipnkappar
! ================
function p1psipsipnkappar(e,f,g,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psipsipnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1. .or. fac1.eq.0) then
     p1psipsipnkappar = 0.
     return
  end if

  temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! [T,psi]*n = [n p/n^2,psi]*n
  temp79b = ni79(:,OP_1)* &
       (i(:,OP_1)*(h(:,OP_DZ)*g(:,OP_DR)-h(:,OP_DR)*g(:,OP_DZ)) &
       +h(:,OP_1)*(i(:,OP_DZ)*g(:,OP_DR)-i(:,OP_DR)*g(:,OP_DZ))) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*g(:,OP_DR) - ni79(:,OP_DR)*g(:,OP_DZ))

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
!!$     temp = 0.
     temp = intx5(e(:,:,OP_1),temp79a,temp79b,norm79(:,1),f(:,OP_DZ)) &
          - intx5(e(:,:,OP_1),temp79a,temp79b,norm79(:,2),f(:,OP_DR))
  else
     temp = intx4(e(:,:,OP_DZ),temp79a,temp79b,f(:,OP_DR)) &
          - intx4(e(:,:,OP_DR),temp79a,temp79b,f(:,OP_DZ))
  end if

  p1psipsipnkappar = (gam - 1.) * temp
end function p1psipsipnkappar

function p1psipsipnkappar1(g,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psipsipnkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  integer, intent(in) :: fac1
  type(prodarray) :: temp

  if(gam.le.1. .or. fac1.eq.0) then
     p1psipsipnkappar1%len = 0
     return
  end if

  temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! [T,psi]*n = [n p/n^2,psi]*n
  temp79b = ni79(:,OP_1)* &
       (i(:,OP_1)*(h(:,OP_DZ)*g(:,OP_DR)-h(:,OP_DR)*g(:,OP_DZ)) &
       +h(:,OP_1)*(i(:,OP_DZ)*g(:,OP_DR)-i(:,OP_DR)*g(:,OP_DZ))) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*g(:,OP_DR) - ni79(:,OP_DR)*g(:,OP_DZ))

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
!!$     temp%len = 0
     temp = prod( temp79a*temp79b*norm79(:,1),OP_1,OP_DZ) &
          + prod(-temp79a*temp79b*norm79(:,2),OP_1,OP_DR)
  else
     temp = prod( temp79a*temp79b,OP_DZ,OP_DR) &
          + prod(-temp79a*temp79b,OP_DR,OP_DZ)
  end if

  temp79a=(gam-1.)
  p1psipsipnkappar1 = temp*temp79a
end function p1psipsipnkappar1

function p1psipsipnkappar2(f,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psipsipnkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  integer, intent(in) :: fac1
  type(prodarray) :: temp
  type(muarray) tempb

  if(gam.le.1. .or. fac1.eq.0) then
     p1psipsipnkappar2%len = 0
     return
  end if

  temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! [T,psi]*n = [n p/n^2,psi]*n
  tempb = &
          mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DR),OP_DZ) &
        + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DR),OP_DZ) &
        + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZ),OP_DR)+mu(-2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DR),OP_DZ)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
!!$     temp%len = 0
     temp = prod(mu( temp79a*norm79(:,1)*f(:,OP_DZ),OP_1),tempb) &
          + prod(mu(-temp79a*norm79(:,2)*f(:,OP_DR),OP_1),tempb)
  else
     temp = prod(mu( temp79a*f(:,OP_DR),OP_DZ),tempb) &
          + prod(mu(-temp79a*f(:,OP_DZ),OP_DR),tempb)
  end if

  temp79a=(gam-1.)
  p1psipsipnkappar2 = temp*temp79a
end function p1psipsipnkappar2

function p1psipsipnkappar3(f,g,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psipsipnkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i
  integer, intent(in) :: fac1
  type(prodarray) :: temp
  type(muarray) tempb

  if(gam.le.1. .or. fac1.eq.0) then
     p1psipsipnkappar3%len = 0
     return
  end if

  temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! [T,psi]*n = [n p/n^2,psi]*n
  tempb = &
          mu(ni79(:,OP_1)*i(:,OP_1)*g(:,OP_DR),OP_DZ)+mu(-ni79(:,OP_1)*i(:,OP_1)*g(:,OP_DZ),OP_DR) &
        + mu(ni79(:,OP_1)*g(:,OP_DR)*i(:,OP_DZ),OP_1)+mu(-ni79(:,OP_1)*g(:,OP_DZ)*i(:,OP_DR),OP_1) &
        + mu(2.*g(:,OP_DR)*i(:,OP_1)*ni79(:,OP_DZ),OP_1)+mu(-2.*g(:,OP_DZ)*i(:,OP_1)*ni79(:,OP_DR),OP_1)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
!!$     temp%len = 0
     temp = prod(mu( temp79a*norm79(:,1)*f(:,OP_DZ),OP_1),tempb) &
          + prod(mu(-temp79a*norm79(:,2)*f(:,OP_DR),OP_1),tempb)
  else
     temp = prod(mu( temp79a*f(:,OP_DR),OP_DZ),tempb) &
          + prod(mu(-temp79a*f(:,OP_DZ),OP_DR),tempb)
  end if

  temp79a=(gam-1.)
  p1psipsipnkappar3 = temp*temp79a
end function p1psipsipnkappar3

function p1psipsipnkappar4(f,g,h,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psipsipnkappar4
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  integer, intent(in) :: fac1
  type(prodarray) :: temp
  type(muarray) tempb

  if(gam.le.1. .or. fac1.eq.0) then
     p1psipsipnkappar4%len = 0
     return
  end if

  temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! [T,psi]*n = [n p/n^2,psi]*n
  tempb = &
          mu(ni79(:,OP_1)*h(:,OP_DZ)*g(:,OP_DR),OP_1)+mu(-ni79(:,OP_1)*h(:,OP_DR)*g(:,OP_DZ),OP_1) &
        + mu(ni79(:,OP_1)*g(:,OP_DR)*h(:,OP_1),OP_DZ)+mu(-ni79(:,OP_1)*g(:,OP_DZ)*h(:,OP_1),OP_DR) &
        + mu(2.*g(:,OP_DR)*h(:,OP_1)*ni79(:,OP_DZ),OP_1)+mu(-2.*g(:,OP_DZ)*h(:,OP_1)*ni79(:,OP_DR),OP_1)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
!!$     temp%len = 0
     temp = prod(mu( temp79a*norm79(:,1)*f(:,OP_DZ),OP_1),tempb) &
          + prod(mu(-temp79a*norm79(:,2)*f(:,OP_DR),OP_1),tempb)
  else
     temp = prod(mu( temp79a*f(:,OP_DR),OP_DZ),tempb) &
          + prod(mu(-temp79a*f(:,OP_DZ),OP_DR),tempb)
  end if

  temp79a=(gam-1.)
  p1psipsipnkappar4 = temp*temp79a
end function p1psipsipnkappar4

! P1psibkappar
! ============
function p1psibkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psibkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psibkappar = 0.
     return
  end if

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
!!$     temp = 0.
     temp79a = -ri3_79*k(:,OP_1)*j(:,OP_1)*g(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))

     temp = intx4(e(:,:,OP_1),temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + intx4(e(:,:,OP_1),temp79a,h(:,OP_1 ),i(:,OP_DP))
  else
     temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)*j(:,OP_1)  

     temp79b = f(:,OP_DR)*(h(:,OP_DZ)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DZ)) &
          -    f(:,OP_DZ)*(h(:,OP_DR)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DR))

     temp79c = f(:,OP_DRP)*(h(:,OP_DZ )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DZ )) &
          -    f(:,OP_DZP)*(h(:,OP_DR )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DR )) &
          +    f(:,OP_DR )*(h(:,OP_DZP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DZ )) &
          -    f(:,OP_DZ )*(h(:,OP_DRP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DR )) &
          +    f(:,OP_DR )*(h(:,OP_DZ )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DZP)) &
          -    f(:,OP_DZ )*(h(:,OP_DR )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DRP))

     temp79d = temp79c*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*g(:,OP_DP)*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*g(:,OP_1 )*j(:,OP_DP)*k(:,OP_1 ) &
          +    temp79b*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_DP)

     temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,h(:,OP_DP),i(:,OP_1 )) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
          + intx3(e(:,:,OP_1),ri3_79,temp79d)
  end if
  p1psibkappar = (gam - 1.) * temp
#else
  p1psibkappar = 0.
#endif
end function p1psibkappar

! P1psibpnkappar
! ==============
function p1psibpnkappar(e,f,g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psibpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psibpnkappar = 0.
     return
  end if

  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1)

  ! n*dT/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
!     temp = fac1*int5(ri3_79,temp79a,temp79e,norm79(:,2),f(:,OP_DR)) &
!          - fac1*int5(ri3_79,temp79a,temp79e,norm79(:,1),f(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_DP)*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_DP)*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]*n = [n p/n^2,psi]*n
     temp79c = ni79(:,OP_1)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))

     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZP)*f(:,OP_DR) - ni79(:,OP_DRP)*f(:,OP_DZ)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DRP)-h(:,OP_DR)*f(:,OP_DZP)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DRP)-i(:,OP_DR)*f(:,OP_DZP))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DRP) - ni79(:,OP_DR)*f(:,OP_DZP)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_1 )*(h(:,OP_DZP)*f(:,OP_DR)-h(:,OP_DRP)*f(:,OP_DZ)) &
          +h(:,OP_DP)*(i(:,OP_DZ )*f(:,OP_DR)-i(:,OP_DR )*f(:,OP_DZ))) &
          + 2.*h(:,OP_DP)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_DP)*(h(:,OP_DZ )*f(:,OP_DR)-h(:,OP_DR )*f(:,OP_DZ)) &
          +h(:,OP_1 )*(i(:,OP_DZP)*f(:,OP_DR)-i(:,OP_DRP)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_DP)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))

     temp = fac2*intx4(e(:,:,OP_1),ri3_79,temp79a,temp79d) &
          + fac2*intx4(e(:,:,OP_1),ri3_79,temp79b,temp79c) &
          + fac1*intx5(e(:,:,OP_DR),ri3_79,temp79a,f(:,OP_DZ),temp79e) &
          - fac1*intx5(e(:,:,OP_DZ),ri3_79,temp79a,f(:,OP_DR),temp79e)
  end if
  p1psibpnkappar = (gam - 1.) * temp
#else
  p1psibpnkappar = 0.
#endif
end function p1psibpnkappar

function p1psibpnkappar1(g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psibpnkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc, tempd

  if(gam.le.1.) then
     p1psibpnkappar1%len = 0
     return
  end if

  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1)

  ! n*dT/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
!     temp = fac1*int5(ri3_79,temp79a,temp79e,norm79(:,2),f(:,OP_DR)) &
!          - fac1*int5(ri3_79,temp79a,temp79e,norm79(:,1),f(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_DP)*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_DP)*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]*n = [n p/n^2,psi]*n
     tempc = &
          mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DR),OP_DZ) &
        + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DR),OP_DZ) &
        + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZ),OP_DR)+mu(-2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DR),OP_DZ)

     ! d(temp79c)/dphi
     tempd = &
          mu(ni79(:,OP_DP)*i(:,OP_1)*h(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_DP)*i(:,OP_1)*h(:,OP_DR),OP_DZ) &
        + mu(ni79(:,OP_DP)*h(:,OP_1)*i(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_DP)*h(:,OP_1)*i(:,OP_DR),OP_DZ) &
        + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZP),OP_DR)+mu(-2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DRP),OP_DZ) &
        + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZ),OP_DRP)+mu(-ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DR),OP_DZP) &
        + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZ),OP_DRP)+mu(-ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DR),OP_DZP) &
        + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZ),OP_DRP)+mu(-2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DR),OP_DZP) &
        + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZP),OP_DR)+mu(-ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DRP),OP_DZ) &
        + mu(ni79(:,OP_1)*h(:,OP_DP)*i(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_1)*h(:,OP_DP)*i(:,OP_DR),OP_DZ) &
        + mu(2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DZ),OP_DR)+mu(-2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DR),OP_DZ) &
        + mu(ni79(:,OP_1)*i(:,OP_DP)*h(:,OP_DZ),OP_DR)+mu(-ni79(:,OP_1)*i(:,OP_DP)*h(:,OP_DR),OP_DZ) &
        + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZP),OP_DR)+mu(-ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DRP),OP_DZ) &
        + mu(2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DZ),OP_DR)+mu(-2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DR),OP_DZ)

     temp = prod(mu(fac2*ri3_79*temp79a,OP_1),tempd) &
          + prod(mu(fac2*ri3_79*temp79b,OP_1),tempc) &
          + prod(fac1*ri3_79*temp79a*temp79e,OP_DR,OP_DZ) &
          + prod(-fac1*ri3_79*temp79a*temp79e,OP_DZ,OP_DR)
  end if
  temp79a=(gam-1.)
  p1psibpnkappar1 = temp*temp79a
#else
  p1psibpnkappar1%len = 0
#endif
end function p1psibpnkappar1

function p1psibpnkappar2(f,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psibpnkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(gam.le.1.) then
     p1psibpnkappar2%len = 0
     return
  end if

  tempa = mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1),OP_1)

  ! n*dT/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
!     temp = fac1*int5(ri3_79,temp79a,temp79e,norm79(:,2),f(:,OP_DR)) &
!          - fac1*int5(ri3_79,temp79a,temp79e,norm79(:,1),f(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     tempb = mu(kar79(:,OP_DP)*b2i79(:,OP_1 )*ni79(:,OP_1 ),OP_1 ) &
          +    mu(kar79(:,OP_1 )*b2i79(:,OP_DP)*ni79(:,OP_1 ),OP_1 ) &
          +    mu(kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_DP),OP_1 ) &
          +    mu(kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_1 ),OP_DP)

     ! [T,psi]*n = [n p/n^2,psi]*n
     temp79c = ni79(:,OP_1)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))

     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZP)*f(:,OP_DR) - ni79(:,OP_DRP)*f(:,OP_DZ)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DRP)-h(:,OP_DR)*f(:,OP_DZP)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DRP)-i(:,OP_DR)*f(:,OP_DZP))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DRP) - ni79(:,OP_DR)*f(:,OP_DZP)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_1 )*(h(:,OP_DZP)*f(:,OP_DR)-h(:,OP_DRP)*f(:,OP_DZ)) &
          +h(:,OP_DP)*(i(:,OP_DZ )*f(:,OP_DR)-i(:,OP_DR )*f(:,OP_DZ))) &
          + 2.*h(:,OP_DP)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_DP)*(h(:,OP_DZ )*f(:,OP_DR)-h(:,OP_DR )*f(:,OP_DZ)) &
          +h(:,OP_1 )*(i(:,OP_DZP)*f(:,OP_DR)-i(:,OP_DRP)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_DP)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))

     temp = prod(mu(fac2*ri3_79*temp79d,OP_1),tempa) &
          + prod(mu(fac2*ri3_79*temp79c,OP_1),tempb) &
          + prod(mu(fac1*ri3_79*f(:,OP_DZ)*temp79e,OP_DR),tempa) &
          + prod(mu(-fac1*ri3_79*f(:,OP_DR)*temp79e,OP_DZ),tempa)
  end if
  temp79a=(gam-1.)
  p1psibpnkappar2 = temp*temp79a
#else
  p1psibpnkappar2%len = 0
#endif
end function p1psibpnkappar2

function p1psibpnkappar3(f,g,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psibpnkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc, tempd,tempe

  if(gam.le.1.) then
     p1psibpnkappar3%len = 0
     return
  end if

  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1)

  ! n*dT/dphi
  tempe = mu(ni79(:,OP_1)*i(:,OP_DP),OP_1) + mu(ni79(:,OP_1)*i(:,OP_1),OP_DP) &
          + mu(2.*i(:,OP_1)*ni79(:,OP_DP),OP_1)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
!     temp = fac1*int5(ri3_79,temp79a,temp79e,norm79(:,2),f(:,OP_DR)) &
!          - fac1*int5(ri3_79,temp79a,temp79e,norm79(:,1),f(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_DP)*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_DP)*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]*n = [n p/n^2,psi]*n
     tempc = &
          mu(ni79(:,OP_1)*i(:,OP_1)*f(:,OP_DR),OP_DZ)+mu(-ni79(:,OP_1)*i(:,OP_1)*f(:,OP_DZ),OP_DR) &
        + mu(ni79(:,OP_1)*f(:,OP_DR)*i(:,OP_DZ),OP_1)+mu(-ni79(:,OP_1)*f(:,OP_DZ)*i(:,OP_DR),OP_1) &
        + mu(2.*f(:,OP_DR)*i(:,OP_1)*ni79(:,OP_DZ),OP_1)+mu(-2.*f(:,OP_DZ)*i(:,OP_1)*ni79(:,OP_DR),OP_1)

     ! d(temp79c)/dphi
     tempd = &
          mu(ni79(:,OP_DP)*i(:,OP_1)*f(:,OP_DR),OP_DZ)+mu(-ni79(:,OP_DP)*i(:,OP_1)*f(:,OP_DZ),OP_DR) &
        + mu(ni79(:,OP_DP)*f(:,OP_DR)*i(:,OP_DZ),OP_1)+mu(-ni79(:,OP_DP)*f(:,OP_DZ)*i(:,OP_DR),OP_1) &
        + mu(2.*f(:,OP_DR)*i(:,OP_1)*ni79(:,OP_DZP),OP_1)+mu(-2.*f(:,OP_DZ)*i(:,OP_1)*ni79(:,OP_DRP),OP_1) &
        + mu(ni79(:,OP_1)*i(:,OP_1)*f(:,OP_DRP),OP_DZ)+mu(-ni79(:,OP_1)*i(:,OP_1)*f(:,OP_DZP),OP_DR) &
        + mu(ni79(:,OP_1)*f(:,OP_DRP)*i(:,OP_DZ),OP_1)+mu(-ni79(:,OP_1)*f(:,OP_DZP)*i(:,OP_DR),OP_1) &
        + mu(2.*f(:,OP_DRP)*i(:,OP_1)*ni79(:,OP_DZ),OP_1)+mu(-2.*f(:,OP_DZP)*i(:,OP_1)*ni79(:,OP_DR),OP_1) &
        + mu(ni79(:,OP_1)*i(:,OP_1)*f(:,OP_DR),OP_DZP)+mu(-ni79(:,OP_1)*i(:,OP_1)*f(:,OP_DZ),OP_DRP) &
        + mu(ni79(:,OP_1)*f(:,OP_DR)*i(:,OP_DZ),OP_DP)+mu(-ni79(:,OP_1)*f(:,OP_DZ)*i(:,OP_DR),OP_DP) &
        + mu(2.*f(:,OP_DR)*i(:,OP_1)*ni79(:,OP_DZ),OP_DP)+mu(-2.*f(:,OP_DZ)*i(:,OP_1)*ni79(:,OP_DR),OP_DP) &
        + mu(ni79(:,OP_1)*i(:,OP_DP)*f(:,OP_DR),OP_DZ)+mu(-ni79(:,OP_1)*i(:,OP_DP)*f(:,OP_DZ),OP_DR) &
        + mu(ni79(:,OP_1)*f(:,OP_DR)*i(:,OP_DZP),OP_1)+mu(-ni79(:,OP_1)*f(:,OP_DZ)*i(:,OP_DRP),OP_1) &
        + mu(2.*f(:,OP_DR)*i(:,OP_DP)*ni79(:,OP_DZ),OP_1)+mu(-2.*f(:,OP_DZ)*i(:,OP_DP)*ni79(:,OP_DR),OP_1)

     temp = prod(mu(fac2*ri3_79*temp79a,OP_1),tempd) &
          + prod(mu(fac2*ri3_79*temp79b,OP_1),tempc) &
          + prod(mu(fac1*ri3_79*temp79a*f(:,OP_DZ),OP_DR),tempe) &
          + prod(mu(-fac1*ri3_79*temp79a*f(:,OP_DR),OP_DZ),tempe)
  end if
  temp79a=(gam-1.)
  p1psibpnkappar3 = temp*temp79a
#else
  p1psibpnkappar3%len = 0
#endif
end function p1psibpnkappar3

function p1psibpnkappar4(f,g,h,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psibpnkappar4
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc, tempd,tempe

  if(gam.le.1.) then
     p1psibpnkappar4%len = 0
     return
  end if

  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1)

  ! n*dT/dphi
  tempe = mu(ni79(:,OP_1)*h(:,OP_1),OP_DP) + mu(ni79(:,OP_1)*h(:,OP_DP),OP_1) &
          + mu(2.*h(:,OP_1)*ni79(:,OP_DP),OP_1)

  if(surface_int) then
!!$     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
!     temp = fac1*int5(ri3_79,temp79a,temp79e,norm79(:,2),f(:,OP_DR)) &
!          - fac1*int5(ri3_79,temp79a,temp79e,norm79(:,1),f(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_DP)*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_DP)*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]*n = [n p/n^2,psi]*n
     tempc = &
          mu(ni79(:,OP_1)*h(:,OP_DZ)*f(:,OP_DR),OP_1)+mu(-ni79(:,OP_1)*h(:,OP_DR)*f(:,OP_DZ),OP_1) &
        + mu(ni79(:,OP_1)*f(:,OP_DR)*h(:,OP_1),OP_DZ)+mu(-ni79(:,OP_1)*f(:,OP_DZ)*h(:,OP_1),OP_DR) &
        + mu(2.*f(:,OP_DR)*h(:,OP_1)*ni79(:,OP_DZ),OP_1)+mu(-2.*f(:,OP_DZ)*h(:,OP_1)*ni79(:,OP_DR),OP_1)

     ! d(temp79c)/dphi
     tempd = &
          mu(ni79(:,OP_DP)*h(:,OP_DZ)*f(:,OP_DR),OP_1)+mu(-ni79(:,OP_DP)*h(:,OP_DR)*f(:,OP_DZ),OP_1) &
        + mu(ni79(:,OP_DP)*f(:,OP_DR)*h(:,OP_1),OP_DZ)+mu(-ni79(:,OP_DP)*f(:,OP_DZ)*h(:,OP_1),OP_DR) &
        + mu(2.*f(:,OP_DR)*h(:,OP_1)*ni79(:,OP_DZP),OP_1)+mu(-2.*f(:,OP_DZ)*h(:,OP_1)*ni79(:,OP_DRP),OP_1) &
        + mu(ni79(:,OP_1)*h(:,OP_DZ)*f(:,OP_DRP),OP_1)+mu(-ni79(:,OP_1)*h(:,OP_DR)*f(:,OP_DZP),OP_1) &
        + mu(ni79(:,OP_1)*f(:,OP_DRP)*h(:,OP_1),OP_DZ)+mu(-ni79(:,OP_1)*f(:,OP_DZP)*h(:,OP_1),OP_DR) &
        + mu(2.*f(:,OP_DRP)*h(:,OP_1)*ni79(:,OP_DZ),OP_1)+mu(-2.*f(:,OP_DZP)*h(:,OP_1)*ni79(:,OP_DR),OP_1) &
        + mu(ni79(:,OP_1)*h(:,OP_DZP)*f(:,OP_DR),OP_1)+mu(-ni79(:,OP_1)*h(:,OP_DRP)*f(:,OP_DZ),OP_1) &
        + mu(ni79(:,OP_1)*f(:,OP_DR)*h(:,OP_DP),OP_DZ)+mu(-ni79(:,OP_1)*f(:,OP_DZ)*h(:,OP_DP),OP_DR) &
        + mu(2.*f(:,OP_DR)*h(:,OP_DP)*ni79(:,OP_DZ),OP_1)+mu(-2.*f(:,OP_DZ)*h(:,OP_DP)*ni79(:,OP_DR),OP_1) &
        + mu(ni79(:,OP_1)*h(:,OP_DZ)*f(:,OP_DR),OP_DP)+mu(-ni79(:,OP_1)*h(:,OP_DR)*f(:,OP_DZ),OP_DP) &
        + mu(ni79(:,OP_1)*f(:,OP_DR)*h(:,OP_1),OP_DZP)+mu(-ni79(:,OP_1)*f(:,OP_DZ)*h(:,OP_1),OP_DRP) &
        + mu(2.*f(:,OP_DR)*h(:,OP_1)*ni79(:,OP_DZ),OP_DP)+mu(-2.*f(:,OP_DZ)*h(:,OP_1)*ni79(:,OP_DR),OP_DP)

     temp = prod(mu(fac2*ri3_79*temp79a,OP_1),tempd) &
          + prod(mu(fac2*ri3_79*temp79b,OP_1),tempc) &
          + prod(mu(fac1*ri3_79*temp79a*f(:,OP_DZ),OP_DR),tempe) &
          + prod(mu(-fac1*ri3_79*temp79a*f(:,OP_DR),OP_DZ),tempe)
  end if
  temp79a=(gam-1.)
  p1psibpnkappar4 = temp*temp79a
#else
  p1psibpnkappar4%len = 0
#endif
end function p1psibpnkappar4

! P1bbkappar
! ==========
function p1bbkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bbkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1bbkappar = 0.
     return
  end if

  if(surface_int) then
     temp = 0.
  else
     temp79a = h(:,OP_DP)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DP)
     temp79b = h(:,OP_DPP)*i(:,OP_1  ) &
          + 2.*h(:,OP_DP )*i(:,OP_DP ) &
          +    h(:,OP_1  )*i(:,OP_DPP)

     temp79c = f(:,OP_DP)*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_DP)*temp79a*j(:,OP_1 )*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_1 )*temp79b*j(:,OP_1 )*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_1 )*temp79a*j(:,OP_DP)*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_DP)

     temp = intx3(e(:,:,OP_1),ri4_79,temp79c)
  end if
  p1bbkappar = (gam - 1.) * temp
#else
  p1bbkappar = 0.
#endif
end function p1bbkappar


! P1bbpnkappar
! ===========
function p1bbpnkappar(e,f,g,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bbpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  
  if(gam.le.1. .or. fac1.eq.0) then
     p1bbpnkappar = 0.
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
  else
     temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1)

     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

     ! n*dT/dphi
     temp79c = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)
     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DPP) &
          +    ni79(:,OP_1)*(h(:,OP_DP)*i(:,OP_DP) + h(:,OP_DPP)*i(:,OP_1)) &
          + 2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DP) &
          +    ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DPP) + h(:,OP_DP)*i(:,OP_DP)) &
          + 2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DP)

     temp = intx4(e(:,:,OP_1),ri4_79,temp79a,temp79d) &
          + intx4(e(:,:,OP_1),ri4_79,temp79b,temp79c)
  end if
  p1bbpnkappar = (gam - 1.) * temp
#else
  p1bbpnkappar = 0.
#endif
end function p1bbpnkappar

function p1bbpnkappar1(g,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1bbpnkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i
  integer, intent(in) :: fac1

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa
  type(muarray) :: tempb
  
  if(gam.le.1. .or. fac1.eq.0) then
     p1bbpnkappar1%len = 0
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
  else
     tempa = mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1),OP_1)

     ! d(temp79a)/dphi
     tempb = mu(kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1),OP_1) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*g(:,OP_1),OP_1) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*g(:,OP_1),OP_1) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1),OP_DP) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_DP),OP_1)

     ! n*dT/dphi
     temp79c = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)
     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DPP) &
          +    ni79(:,OP_1)*(h(:,OP_DP)*i(:,OP_DP) + h(:,OP_DPP)*i(:,OP_1)) &
          + 2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DP) &
          +    ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DPP) + h(:,OP_DP)*i(:,OP_DP)) &
          + 2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DP)

     temp = prod(mu(ri4_79*temp79d,OP_1),tempa) &
          + prod(mu(ri4_79*temp79c,OP_1),tempb)
  end if
  temp79a=(gam-1.)
  p1bbpnkappar1 = temp*temp79a
#else
  p1bbpnkappar1%len = 0
#endif
end function p1bbpnkappar1

function p1bbpnkappar2(f,h,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1bbpnkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  integer, intent(in) :: fac1

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa
  type(muarray) :: tempb
  
  if(gam.le.1. .or. fac1.eq.0) then
     p1bbpnkappar2%len = 0
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
  else
     tempa = mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1),OP_1)

     ! d(temp79a)/dphi
     tempb = mu(kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1),OP_1) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1),OP_1) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1),OP_1) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP),OP_1) &
          +    mu(kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1),OP_DP)

     ! n*dT/dphi
     temp79c = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)
     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DPP) &
          +    ni79(:,OP_1)*(h(:,OP_DP)*i(:,OP_DP) + h(:,OP_DPP)*i(:,OP_1)) &
          + 2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DP) &
          +    ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DPP) + h(:,OP_DP)*i(:,OP_DP)) &
          + 2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DP)

     temp = prod(mu(ri4_79*temp79d,OP_1),tempa) &
          + prod(mu(ri4_79*temp79c,OP_1),tempb)
  end if
  temp79a=(gam-1.)
  p1bbpnkappar2 = temp*temp79a
#else
  p1bbpnkappar2%len = 0
#endif
end function p1bbpnkappar2

function p1bbpnkappar3(f,g,i,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1bbpnkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i
  integer, intent(in) :: fac1

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc
  type(muarray) :: tempd
  
  if(gam.le.1. .or. fac1.eq.0) then
     p1bbpnkappar3%len = 0
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
  else
     temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1)

     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

     ! n*dT/dphi
     tempc = mu(ni79(:,OP_1)*i(:,OP_DP),OP_1) + mu(ni79(:,OP_1)*i(:,OP_1),OP_DP) &
           + mu(2.*i(:,OP_1)*ni79(:,OP_DP),OP_1)
     ! d(temp79c)/dphi
     tempd = mu(ni79(:,OP_DP)*i(:,OP_DP),OP_1) + mu(ni79(:,OP_DP)*i(:,OP_1),OP_DP) &
           + mu(2.*i(:,OP_1)*ni79(:,OP_DPP),OP_1) &
           + mu(ni79(:,OP_1)*i(:,OP_DP),OP_DP) + mu(ni79(:,OP_1)*i(:,OP_1),OP_DPP) &
           + mu(2.*i(:,OP_1)*ni79(:,OP_DP),OP_DP) &
           + mu(ni79(:,OP_1)*i(:,OP_DPP),OP_1) + mu(ni79(:,OP_1)*i(:,OP_DP),OP_DP) &
           + mu(2.*i(:,OP_DP)*ni79(:,OP_DP),OP_1)
       
     temp = prod(mu(ri4_79*temp79a,OP_1),tempd) &
          + prod(mu(ri4_79*temp79b,OP_1),tempc)
  end if
  temp79a=(gam-1.)
  p1bbpnkappar3 = temp*temp79a
#else
  p1bbpnkappar3%len = 0
#endif
end function p1bbpnkappar3

function p1bbpnkappar4(f,g,h,fac1)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1bbpnkappar4
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  integer, intent(in) :: fac1

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc
  type(muarray) :: tempd
  
  if(gam.le.1. .or. fac1.eq.0) then
     p1bbpnkappar4%len = 0
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
  else
     temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1)

     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

     ! n*dT/dphi
     tempc = mu(ni79(:,OP_1)*h(:,OP_1),OP_DP) + mu(ni79(:,OP_1)*h(:,OP_DP),OP_1) &
           + mu(2.*h(:,OP_1)*ni79(:,OP_DP),OP_1)
     ! d(temp79c)/dphi
     tempd = mu(ni79(:,OP_DP)*h(:,OP_1),OP_DP) + mu(ni79(:,OP_DP)*h(:,OP_DP),OP_1) &
           + mu(2.*h(:,OP_1)*ni79(:,OP_DPP),OP_1) &
           + mu(ni79(:,OP_1)*h(:,OP_DP),OP_DP) + mu(ni79(:,OP_1)*h(:,OP_DPP),OP_1) &
           + mu(2.*h(:,OP_DP)*ni79(:,OP_DP),OP_1) &
           + mu(ni79(:,OP_1)*h(:,OP_1),OP_DPP) + mu(ni79(:,OP_1)*h(:,OP_DP),OP_DP) &
           + mu(2.*h(:,OP_1)*ni79(:,OP_DP),OP_DP)
       
     temp = prod(mu(ri4_79*temp79a,OP_1),tempd) &
          + prod(mu(ri4_79*temp79b,OP_1),tempc)
  end if
  temp79a=(gam-1.)
  p1bbpnkappar4 = temp*temp79a
#else
  p1bbpnkappar4%len = 0
#endif
end function p1bbpnkappar4


! P1psifkappar
! ============
function p1psifkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psifkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psifkappar = 0.
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
!!$     temp79a = k(:,OP_1)*ri_79*e(:,OP_1)* &
!!$          (norm79(:,2)*f(:,OP_DR) - norm79(:,1)*f(:,OP_DZ))*j(:,OP_1)
!!$     temp79b = -k(:,OP_1)*ri_79*e(:,OP_1)* &
!!$          (norm79(:,2)*g(:,OP_DZ) + norm79(:,1)*g(:,OP_DR))*j(:,OP_1)
!!$
!!$     temp = int4(temp79a,g(:,OP_DZ),h(:,OP_DZ),i(:,OP_1 )) &
!!$          + int4(temp79a,g(:,OP_DR),h(:,OP_DR),i(:,OP_1 )) &
!!$          + int4(temp79a,g(:,OP_DZ),h(:,OP_1 ),i(:,OP_DZ)) &
!!$          + int4(temp79a,g(:,OP_DR),h(:,OP_1 ),i(:,OP_DR)) &
!!$          + int4(temp79b,f(:,OP_DR ),h(:,OP_DZ),i(:,OP_1 )) &
!!$          - int4(temp79b,f(:,OP_DZ ),h(:,OP_DR),i(:,OP_1 )) &
!!$          + int4(temp79b,f(:,OP_DR ),h(:,OP_1 ),i(:,OP_DZ)) &
!!$          - int4(temp79b,f(:,OP_DZ ),h(:,OP_1 ),i(:,OP_DR))
  else
     temp79a = k(:,OP_1)*ri_79*j(:,OP_1)*i(:,OP_1)
     temp79b = k(:,OP_1)*ri_79*j(:,OP_1)*h(:,OP_1)

     temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79b,g(:,OP_DR),i(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79b,g(:,OP_DR),i(:,OP_DR)) &
          + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
          - intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DZ ),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DZ ),h(:,OP_DR)) &
          + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79b,f(:,OP_DR ),i(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79b,f(:,OP_DR ),i(:,OP_DZ)) &
          - intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79b,f(:,OP_DZ ),i(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),g(:,OP_DR),temp79b,f(:,OP_DZ ),i(:,OP_DR))
  end if
  p1psifkappar = (gam - 1.) * temp
#else
  p1psifkappar = 0.
#endif
end function p1psifkappar


! P1psifpnkappar
! ==============
function p1psifpnkappar(e,f,g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1psifpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1psifpnkappar = 0.
     return
  end if

  temp79a = ri_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! n*<T,f'>
  temp79b = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))
  ! n*[T,psi]
  temp79c = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DZ)*f(:,OP_DR) - i(:,OP_DR)*f(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DZ)*f(:,OP_DR) - h(:,OP_DR)*f(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))
  
  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
!!$     temp = fac1*int5(temp79a,temp79b,e(:,OP_1),norm79(:,1),f(:,OP_DZ)) &
!!$          - fac1*int5(temp79a,temp79b,e(:,OP_1),norm79(:,2),f(:,OP_DR)) &
!!$          - fac2*int5(temp79a,temp79c,e(:,OP_1),norm79(:,1),g(:,OP_DR)) &
!!$          - fac2*int5(temp79a,temp79c,e(:,OP_1),norm79(:,2),g(:,OP_DZ))
  else
     temp = fac1*intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DR),temp79b) &
          - fac1*intx4(e(:,:,OP_DR),temp79a,f(:,OP_DZ),temp79b) &
          + fac2*intx4(e(:,:,OP_DR),temp79a,g(:,OP_DR),temp79c) &
          + fac2*intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DZ),temp79c)
  end if
  p1psifpnkappar = (gam - 1.) * temp
#else
  p1psifpnkappar = 0.
#endif
end function p1psifpnkappar

function p1psifpnkappar2(f,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1psifpnkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempb

  if(gam.le.1.) then
     p1psifpnkappar2%len = 0
     return
  end if

  temp79a = ri_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! n*<T,f'>
  tempb = mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DR),OP_DR) + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZ),OP_DZ) &
          + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DR),OP_DR) + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZ),OP_DZ) &
          + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DR),OP_DR) + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZ),OP_DZ)
  ! n*[T,psi]
  temp79c = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DZ)*f(:,OP_DR) - i(:,OP_DR)*f(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DZ)*f(:,OP_DR) - h(:,OP_DR)*f(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))
  
  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
!!$     temp = fac1*int5(temp79a,temp79b,e(:,OP_1),norm79(:,1),f(:,OP_DZ)) &
!!$          - fac1*int5(temp79a,temp79b,e(:,OP_1),norm79(:,2),f(:,OP_DR)) &
!!$          - fac2*int5(temp79a,temp79c,e(:,OP_1),norm79(:,1),g(:,OP_DR)) &
!!$          - fac2*int5(temp79a,temp79c,e(:,OP_1),norm79(:,2),g(:,OP_DZ))
  else
     temp = prod(mu( fac1*temp79a*f(:,OP_DR),OP_DZ),tempb) &
          + prod(mu(-fac1*temp79a*f(:,OP_DZ),OP_DR),tempb) &
          + prod(fac2*temp79a*temp79c,OP_DR,OP_DR) &
          + prod(fac2*temp79a*temp79c,OP_DZ,OP_DZ)
  end if
  temp79a=(gam-1.)
  p1psifpnkappar2 = temp*temp79a
#else
  p1psifpnkappar2%len = 0
#endif
end function p1psifpnkappar2

! P1qpsikappar
! ==============
function p1qpsikappar(e,f,g,i,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1qpsikappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,k
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1qpsikappar = 0.
     return
  end if

  temp79a = ri_79*k(:,OP_1)*i(:,OP_1)
  
  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
  else
     temp = intx4(e(:,:,OP_DR),temp79a,g(:,OP_DZ),f(:,OP_1)) &
          - intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DR),f(:,OP_1)) 
  end if

  p1qpsikappar = (gam - 1.) * temp
  return
end function p1qpsikappar

! P1bfkappar
! ==========
function p1bfkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bfkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1bfkappar = 0.
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
!!$     temp79a = -ri2_79*k(:,OP_1)*j(:,OP_1)*e(:,OP_1)*f(:,OP_1)* &
!!$          (norm79(:,1)*g(:,OP_DR) + norm79(:,2)*g(:,OP_DZ))
!!$
!!$     temp = int3(temp79a,h(:,OP_DP),i(:,OP_1 )) &
!!$          + int3(temp79a,h(:,OP_1 ),i(:,OP_DP))
  else
     temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)*j(:,OP_1)

     temp79b = g(:,OP_DZ)*(h(:,OP_DZ)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DZ)) &
          +    g(:,OP_DR)*(h(:,OP_DR)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DR))

     temp79c = g(:,OP_DZP)*(h(:,OP_DZ )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DZ )) &
          +    g(:,OP_DRP)*(h(:,OP_DR )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DR )) &
          +    g(:,OP_DZ )*(h(:,OP_DZP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DZ )) &
          +    g(:,OP_DR )*(h(:,OP_DRP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DR )) &
          +    g(:,OP_DZ )*(h(:,OP_DZ )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DZP)) &
          +    g(:,OP_DR )*(h(:,OP_DR )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DRP))

     temp79d = temp79c*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*f(:,OP_DP)*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*f(:,OP_1 )*j(:,OP_DP)*k(:,OP_1 ) &
          +    temp79b*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_DP)

     temp = intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
          + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,h(:,OP_1 ),i(:,OP_DP)) &
          - intx3(e(:,:,OP_1),ri2_79,temp79d)
  end if
  p1bfkappar = (gam - 1.) * temp
#else
  p1bfkappar = 0.
#endif
end function p1bfkappar


! P1qbkappar
! ==========
function p1qbkappar(e,f,g,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1qbkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1qbkappar = 0.
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp79a = 0.
  else
     temp79a =  ri2_79*i(:,OP_1)*j(:,OP_1)*g(:,OP_1)


     temp = -intx3(e(:,:,OP_DP),temp79a,f(:,OP_1 )) 
  end if
  p1qbkappar = (gam - 1.) * temp
#else
  p1qbkappar = 0.
#endif

end function p1qbkappar

! ==========
function p1bfpnkappar(e,f,g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1bfpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1bfpnkappar = 0.
     return
  end if

  temp79a = ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)

  ! n*dT/dphi = n*d(n p/n^2)/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
!!$     temp = -fac2*intx5(e(:,:,OP_1),temp79a,temp79e,norm79(:,1),g(:,OP_DR)) &
!!$          -  fac2*intx5(e(:,:,OP_1),temp79a,temp79e,norm79(:,2),g(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = ri2_79 * &
          (kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP))

     ! n*<T, f'> = n*<n p/n^2, f'>
     temp79c = ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))

     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
          + ni79(:,OP_DP)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DRP)*g(:,OP_DR) + ni79(:,OP_DZP)*g(:,OP_DZ)) &
          +    ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DRP) + i(:,OP_DZ)*g(:,OP_DZP)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DRP) + ni79(:,OP_DZ)*g(:,OP_DZP)) &
          +    ni79(:,OP_1)*h(:,OP_DP)* &
          (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DRP)*g(:,OP_DR) + h(:,OP_DZP)*g(:,OP_DZ)) &
          + 2.*h(:,OP_DP)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ)) &
          +    ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DRP)*g(:,OP_DR) + i(:,OP_DZP)*g(:,OP_DZ)) &
          + ni79(:,OP_1)*i(:,OP_DP)* &
          (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
          + 2.*h(:,OP_1)*i(:,OP_DP)* &
          (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))

     temp = fac2*intx4(e(:,:,OP_DR),temp79a,g(:,OP_DR),temp79e) &
          + fac2*intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DZ),temp79e) &
          - fac1*intx3(e(:,:,OP_1),temp79a,temp79d) &
          - fac1*intx3(e(:,:,OP_1),temp79b,temp79c)
  end if
  p1bfpnkappar = (gam - 1.) * temp
#else
  p1bfpnkappar = 0.
#endif
end function p1bfpnkappar

function p1bfpnkappar2(f,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: p1bfpnkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i
  integer, intent(in) :: fac1, fac2

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc
  type(muarray) :: tempd

  if(gam.le.1.) then
     p1bfpnkappar2%len = 0
     return
  end if

  temp79a = ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)

  ! n*dT/dphi = n*d(n p/n^2)/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp%len = 0
!!$     temp = -fac2*prod(e(:,:,OP_1),temp79a,temp79e,norm79(:,1),g(:,OP_DR)) &
!!$          -  fac2*prod(e(:,:,OP_1),temp79a,temp79e,norm79(:,2),g(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = ri2_79 * &
          (kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP))

     ! n*<T, f'> = n*<n p/n^2, f'>
    tempc = mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DR),OP_DR) + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZ),OP_DZ) &
          + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DR),OP_DR) + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZ),OP_DZ) &
          + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DR),OP_DR) + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZ),OP_DZ)
  
     ! d(temp79c)/dphi
    tempd = mu(ni79(:,OP_DP)*h(:,OP_1)*i(:,OP_DR),OP_DR) + mu(ni79(:,OP_DP)*h(:,OP_1)*i(:,OP_DZ),OP_DZ) &
          + mu(ni79(:,OP_DP)*i(:,OP_1)*h(:,OP_DR),OP_DR) + mu(ni79(:,OP_DP)*i(:,OP_1)*h(:,OP_DZ),OP_DZ) &
          + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DRP),OP_DR) + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZP),OP_DZ) &
          + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DR),OP_DRP) + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZ),OP_DZP) &
          + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DR),OP_DRP) + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZ),OP_DZP) &
          + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DR),OP_DRP) + mu(2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DZ),OP_DZP) &
          + mu(ni79(:,OP_1)*h(:,OP_DP)*i(:,OP_DR),OP_DR) + mu(ni79(:,OP_1)*h(:,OP_DP)*i(:,OP_DZ),OP_DZ) &
          + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DRP),OP_DR) + mu(ni79(:,OP_1)*i(:,OP_1)*h(:,OP_DZP),OP_DZ) &
          + mu(2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DR),OP_DR) + mu(2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DZ),OP_DZ) &
          + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DRP),OP_DR) + mu(ni79(:,OP_1)*h(:,OP_1)*i(:,OP_DZP),OP_DZ) &
          + mu(ni79(:,OP_1)*i(:,OP_DP)*h(:,OP_DR),OP_DR) + mu(ni79(:,OP_1)*i(:,OP_DP)*h(:,OP_DZ),OP_DZ) &
          + mu(2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DR),OP_DR) + mu(2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DZ),OP_DZ)


     temp = prod(fac2*temp79a*temp79e,OP_DR,OP_DR) &
          + prod(fac2*temp79a*temp79e,OP_DZ,OP_DZ) &
          + prod(mu(-fac1*temp79a,OP_1),tempd) &
          + prod(mu(-fac1*temp79b,OP_1),tempc)
  end if
  temp79a=(gam-1.)
  p1bfpnkappar2 = temp*temp79a
#else
  p1bfpnkappar2%len = 0
#endif
end function p1bfpnkappar2

! P1ffkappar
! ==========
function p1ffkappar(e,f,g,h,i,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1ffkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1ffkappar = 0.
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
!!$     temp79a =  k(:,OP_1)*e(:,OP_1)* &
!!$          (norm79(:,2)*f(:,OP_DZ) + norm79(:,1)*f(:,OP_DR))*j(:,OP_1)
!!$
!!$     temp = int4(temp79a,g(:,OP_DZ),h(:,OP_DZ),i(:,OP_1 )) &
!!$          + int4(temp79a,g(:,OP_DR),h(:,OP_DR),i(:,OP_1 )) &
!!$          + int4(temp79a,g(:,OP_DZ),h(:,OP_1 ),i(:,OP_DZ)) &
!!$          + int4(temp79a,g(:,OP_DR),h(:,OP_1 ),i(:,OP_DR)) 
  else
     temp79a = -k(:,OP_1)*j(:,OP_1)*i(:,OP_1)
     temp79b = -k(:,OP_1)*j(:,OP_1)*h(:,OP_1)

     temp = intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79b,g(:,OP_DZ),i(:,OP_DZ)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79b,g(:,OP_DR),i(:,OP_DR)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79b,g(:,OP_DR),i(:,OP_DR)) 
  end if
  p1ffkappar = (gam - 1.) * temp
#else
  p1ffkappar = 0.
#endif
end function p1ffkappar


! P1ffpnkappar
! ============
function p1ffpnkappar(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1ffpnkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1ffpnkappar = 0.
     return
  end if

  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
  else
     temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

     ! <T, f'>*n = <n p/n^2, f'>*n
     temp79c = ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DR) + i(:,OP_DZ)*g(:,OP_DZ)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DR) + ni79(:,OP_DZ)*g(:,OP_DZ))

     temp = -intx4(e(:,:,OP_DR),temp79a,f(:,OP_DR),temp79c) &
          -  intx4(e(:,:,OP_DZ),temp79a,f(:,OP_DZ),temp79c)
  end if
  p1ffpnkappar = (gam - 1.) * temp
#else
  p1ffpnkappar = 0.
#endif
end function p1ffpnkappar


! P1qfkappar
! ============
vectype function p1qfkappar(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  if(gam.le.1.) then
     p1qfkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     ! assert natural b.c. B.grad(T) = 0
     temp = 0.
  else
     temp79a = i(:,OP_1)*h(:,OP_1)


     temp = int4(temp79a,e(:,OP_DR),g(:,OP_DR),f(:,OP_1)) &
          + int4(temp79a,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_1))
  end if
#else
  temp = 0.
#endif

  p1qfkappar = (gam - 1.) * temp
  return
end function p1qfkappar

! P1kappax
! ========
function p1kappax(e,f,g,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1kappax
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,j
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     p1kappax = 0.
     return
  end if

  if(surface_int) then
     temp = 0.
  else
     temp79a = ri_79*i(:,OP_1)*g(:,OP_1)*j(:,OP_1)

     temp = intx3(e(:,:,OP_DZ),f(:,OP_DR),temp79a) &
          - intx3(e(:,:,OP_DR),f(:,OP_DZ),temp79a)
  endif

  p1kappax = (gam - 1.) * temp
end function p1kappax



! P1uus
! =====
function p1uus(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1uus
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1uus = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = 0.5*(gam-1.)* &
       (intx5(e(:,:,OP_1),r2_79,f(:,OP_DZ),g(:,OP_DZ),temp79a) &
       +intx5(e(:,:,OP_1),r2_79,f(:,OP_DR),g(:,OP_DR),temp79a))

  p1uus = temp
end function p1uus


! P1vvs
! =====
function p1vvs(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1vvs
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1vvs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = 0.5*(gam-1.)* &
       intx5(e(:,:,OP_1),r2_79,f(:,OP_1),g(:,OP_1),temp79a)

  p1vvs = temp
end function p1vvs


! P1chichis
! =========
function p1chichis(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1chichis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1chichis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = 0.5*(gam-1.)* &
       (intx5(e(:,:,OP_1),ri4_79,f(:,OP_DZ),g(:,OP_DZ),temp79a) &
       +intx5(e(:,:,OP_1),ri4_79,f(:,OP_DR),g(:,OP_DR),temp79a))
     
  p1chichis = temp
end function p1chichis


! P1uchis
! =======
function p1uchis(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: p1uchis
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1uchis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) ! + denm*nt79(:,OP_LP)

  temp = -(gam-1.)* & 
       (intx5(e(:,:,OP_1),ri_79,f(:,OP_DZ),g(:,OP_DR),temp79a) &
       -intx5(e(:,:,OP_1),ri_79,f(:,OP_DR),g(:,OP_DZ),temp79a))

  p1uchis = temp
end function p1uchis

!
! Extra diffusion to model upstream differencing
vectype function p1uspu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp 

     if(surface_int) then
!....needs to be defined
        temp = 0.
     else
        temp79a = abs(g(:,OP_DZ))
        temp79b = abs(g(:,OP_DR))
!
        temp = - int4(r_79,e(:,OP_DR),f(:,OP_DR),temp79a)  &
               - int4(r_79,e(:,OP_DZ),f(:,OP_DZ),temp79b)
     end if

  p1uspu = temp

  return
end function p1uspu

!
! Extra diffusion to model upstream differencing
vectype function p1uspchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp 

     if(surface_int) then
!....needs to be defined
        temp = 0.
     else
        temp79a = abs(g(:,OP_DR))
        temp79b = abs(g(:,OP_DZ))
!
        temp = - int4(ri2_79,e(:,OP_DR),f(:,OP_DR),temp79a)  &
               - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),temp79b)
     end if

  p1uspchi = temp

  return
end function p1uspchi

! Extra diffusion to model upstream differencing
vectype function p1uspv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

#if defined(USE3D) || defined(USECOMPLEX)
  vectype :: temp 

     if(surface_int) then
!....needs to be defined
        temp = 0.
     else
        temp79a = abs(g(:,OP_1))
!
        temp =  int3(e(:,OP_1),f(:,OP_DPP),temp79a)  
     end if

  p1uspv = temp
#else

  p1uspv = 0.
#endif

  return
end function p1uspv


!======================================================================
! Parallel Viscous Terms
!======================================================================

! PVS1
! ====
subroutine PVS1(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

     o = 0.
     if(itor.eq.1) then
        o = o - (1./3.)*i(:,OP_DZ)
     endif

end subroutine PVS1


! PVS1psipsi
! ==========
subroutine PVS1psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o


     o = r_79* &
          (j(:,OP_DR)*k(:,OP_DZ)*(i(:,OP_DZZ) - i(:,OP_DRR)) &
          -(j(:,OP_DZ)*k(:,OP_DZ) - j(:,OP_DR)*k(:,OP_DR))*i(:,OP_DRZ))
     if(itor.eq.1) then
        o = o + j(:,OP_DR)* &
             (i(:,OP_DZ)*k(:,OP_DR) - i(:,OP_DR)*k(:,OP_DZ))
     endif

     o = o * ri2_79*b2i79(:,OP_1)

end subroutine PVS1psipsi

! PVS1psib
! ========
subroutine PVS1psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

     o = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + k(:,OP_1)* &
          (j(:,OP_DZ)*i(:,OP_DZP) + j(:,OP_DR)*i(:,OP_DRP))
     o = o * ri2_79*b2i79(:,OP_1)
#endif     

end subroutine PVS1psib


! PVS1bb
! ======
subroutine PVS1bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k 
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = 0.

end subroutine PVS1bb


! PVS2
! ====
subroutine PVS2(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

#if defined(USE3D) || defined(USECOMPLEX)
  o = -i(:,OP_DP)/3.
#else
  o = 0.
#endif

end subroutine PVS2


! PVS2psipsi
! ==========
subroutine PVS2psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = 0.
end subroutine PVS2psipsi


! PVS2psib
! ========
subroutine PVS2psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = r_79*k(:,OP_1)* &
       (i(:,OP_DZ)*j(:,OP_DR) - i(:,OP_DR)*j(:,OP_DZ))

  o = o * ri2_79*b2i79(:,OP_1)

end subroutine PVS2psib

! PVS2bb
! ======
subroutine PVS2bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k
  vectype, intent(out), dimension(MAX_PTS) :: o

#if defined(USE3D) || defined(USECOMPLEX)
  o = j(:,OP_1)*k(:,OP_1)*i(:,OP_DP)
  o = o * ri2_79*b2i79(:,OP_1)
#else
  o = 0.
#endif

end subroutine PVS2bb


! PVS3
! ====
subroutine PVS3(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = -(1./3.)*ri2_79*i(:,OP_GS)

  if(itor.eq.1) then
     o = o + (1./3.)*ri3_79*i(:,OP_DR)
  end if
 
end subroutine PVS3


! PVS3psipsi
! ==========
subroutine PVS3psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = -ri2_79* &
       (j(:,OP_DZ)*k(:,OP_DZ)*i(:,OP_DZZ) &
       +j(:,OP_DR)*k(:,OP_DR)*i(:,OP_DRR) &
       +2.*j(:,OP_DZ)*k(:,OP_DR)*i(:,OP_DRZ) &
       -i(:,OP_GS)*(j(:,OP_DZ)*k(:,OP_DZ) + j(:,OP_DR)*k(:,OP_DR))) 
  
  if(itor.eq.1) then
     o = o + 2.*ri3_79*j(:,OP_DZ) * &
          (i(:,OP_DZ)*k(:,OP_DR) - i(:,OP_DR)*k(:,OP_DZ))
     
  endif

  o = o * ri2_79*b2i79(:,OP_1)
 
end subroutine PVS3psipsi


! PVS3psib
! ========
subroutine PVS3psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

#if defined(USE3D) || defined(USECOMPLEX)
  o = ri3_79*k(:,OP_1) * &
       (i(:,OP_DZP)*j(:,OP_DR) - i(:,OP_DRP)*j(:,OP_DZ))
  o = o * ri2_79*b2i79(:,OP_1)
#else
  o = 0.
#endif
 
end subroutine PVS3psib


! PVS3bb
! ======
subroutine PVS3bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = 0.
 
end subroutine PVS3bb



! PVV1
! ====
subroutine PVV1(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

     if(surface_int) then
        o = 3.*ri_79*b2i79(:,OP_1)* &
             (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR))* &
             (norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR)) &
             - r_79*(norm79(:,1)*e(:,OP_DZ) - norm79(:,2)*e(:,OP_DR))
!        o = 0.
     else
        o =   (e(:,OP_DZZ) - e(:,OP_DRR))*pst79(:,OP_DR)*pst79(:,OP_DZ) &
             + e(:,OP_DRZ)*(pst79(:,OP_DR)**2 - pst79(:,OP_DZ)**2)
        
        if(itor.eq.1) then
           o = o - ri_79* &
                (pst79(:,OP_DZ)* &
                (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR)) &
                +e(:,OP_DZ)*bzt79(:,OP_1)**2)
        endif
        
        o = 3.*ri_79*b2i79(:,OP_1)*o
        
        if(itor.eq.1) then
           o = o + 2.*e(:,OP_DZ)
        endif
        
#if defined(USECOMPLEX)
     o = o - rfac*3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
          (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR))
#endif
     end if

  o = -o
end subroutine  PVV1

! PVV2
! ====
subroutine PVV2(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

     if(surface_int) then
!!$     o = 3.*ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)*e(:,OP_1)* &
!!$          (norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR))
        o = 0.
     else
        o = ri_79*(e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ))

        o = -3.*b2i79(:,OP_1)*bzt79(:,OP_1)*o

#if defined(USECOMPLEX)
        ! This term is a total derivative in phi
        ! and therefore integrates out of the 3D case

        o = o - rfac*e(:,OP_1) * &
             (1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)
#endif
     end if
end subroutine  PVV2


! PVV3
! ====
subroutine PVV3(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

     if(surface_int) then
        o = 0.
     else
        o =      e(:,OP_DZZ)*pst79(:,OP_DR)**2 &
             +   e(:,OP_DRR)*pst79(:,OP_DZ)**2 &
             -2.*e(:,OP_DRZ)*pst79(:,OP_DZ)*pst79(:,OP_DR)
        
        if(itor.eq.1) then
           o = o + 2.*ri_79*pst79(:,OP_DZ)* &
                (e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ)) &
                + ri_79*e(:,OP_DR)*bzt79(:,OP_1)**2
        endif
        
        o = 3.*ri4_79*b2i79(:,OP_1)*o - ri2_79*e(:,OP_GS)
        
#if defined(USECOMPLEX)
     o = o - rfac*3.*ri5_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
          (e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ))
#endif
     end if
     
end subroutine  PVV3


! P1vip
! =====
vectype function P1vip(e)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e

  if(amupar.eq.0) then
     P1vip = 0.
     return
  endif

  call PVS1      (pht79,temp79b)
  call PVS1psipsi(pht79,pst79,pst79,temp79c)
  call PVS1psib  (pht79,pst79,bzt79,temp79d)
  call PVS1bb    (pht79,bzt79,bzt79,temp79e)
  temp79a = temp79b + temp79c + temp79d + temp79e

  if(numvar.ge.2) then
     call PVS2      (vzt79,temp79b)
     call PVS2psipsi(vzt79,pst79,pst79,temp79c)
     call PVS2psib  (vzt79,pst79,bzt79,temp79d)
     call PVS2bb    (vzt79,bzt79,bzt79,temp79e)
     temp79a = temp79a + temp79b + temp79c + temp79d + temp79e
  endif

  if(numvar.ge.3) then
     call PVS3      (cht79,temp79b)
     call PVS3psipsi(cht79,pst79,pst79,temp79c)
     call PVS3psib  (cht79,pst79,bzt79,temp79d)
     call PVS3bb    (cht79,bzt79,bzt79,temp79e)
     temp79a = temp79a + temp79b + temp79c + temp79d + temp79e
  endif

  P1vip = 3.*int4(e(:,OP_1),vip79(:,OP_1),temp79a,temp79a)
end function P1vip


!======================================================================
! Gyroviscous terms
!======================================================================

! g1u
! ===
vectype function g1u(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g1u = 0.
     return
  end if
  
  g1u = gyro_vor_u(e,f)

end function g1u

! g1v
! ===
vectype function g1v(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g1v = 0.
     return
  end if

  g1v = gyro_vor_v(e,f)

end function g1v

! g1chi
! =====
vectype function g1chi(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g1chi = 0.
     return
  end if

  g1chi = gyro_vor_x(e,f)

end function g1chi


! g2u
! ===
vectype function g2u(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g2u = 0.
     return
  end if

  g2u = gyro_tor_u(e,f)

end function g2u


! g2v
! ===
vectype function g2v(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g2v = 0.
     return
  end if

  g2v = gyro_tor_v(e,f)

end function g2v


! g2chi
! =====
vectype function g2chi(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g2chi = 0.
     return
  end if

  g2chi = gyro_tor_x(e,f)

end function g2chi


! g3u
! ===
vectype function g3u(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g3u = 0.
     return
  end if

  g3u = gyro_com_u(e,f)

end function g3u


! g3v
! ===
vectype function g3v(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g3v = 0.
     return
  end if

  g3v = gyro_com_v(e,f)

end function g3v


! g3chi
! =====
vectype function g3chi(e,f)

  use basic
  use m3dc1_nint
  use gyroviscosity

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g3chi = 0.
     return
  end if

  g3chi = gyro_com_x(e,f)

end function g3chi





! ==============================================================
! Ohmic heating terms
! ==============================================================

! qpsipsieta
! ==========
function qpsipsieta(e)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qpsipsieta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, dimension(dofs_per_element) :: temp

  if(hypf.eq.0. .or. jadv.eq.1 .or. surface_int) then
     qpsipsieta = 0.
     return
  end if

  temp79a = ri2_79* &
       (jt79(:,OP_DZ)**2 + jt79(:,OP_DR)**2)
  if(itor.eq.1) then
     temp79a = temp79a - 2.*jt79(:,OP_1)* &
          (ri3_79*jt79(:,OP_DR) - ri4_79*jt79(:,OP_1))
  endif

  if(ihypeta.eq.1) then
     temp = hypf*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = hypf*intx2(e(:,:,OP_1),temp79a)
  endif


  temp79a = jt79(:,OP_DZ)*ni79(:,OP_DZ) + jt79(:,OP_DR)*ni79(:,OP_DR)
  if(itor.eq.1) then
     temp79a = temp79a - ri_79*jt79(:,OP_1)*ni79(:,OP_DR)
  endif
  temp79a = temp79a * ri2_79*nt79(:,OP_1)*jt79(:,OP_1)
     
  if(ihypeta.eq.1) then
     temp = temp + hypf*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = temp + hypf*intx2(e(:,:,OP_1),temp79a)
  endif

  qpsipsieta = temp
end function qpsipsieta

! qbbeta
! ======
function qbbeta(e)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qbbeta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, dimension(dofs_per_element) :: temp

  if(hypi.eq.0. .or. surface_int) then
     qbbeta = 0.
     return
  end if

  temp79a = ri2_79*(bzt79(:,OP_GS)*bzt79(:,OP_GS) &
       + 2.*(bzt79(:,OP_DRZ)**2 - bzt79(:,OP_DRR)*bzt79(:,OP_DZZ)))
  
  if(itor.eq.1) then 
     temp79a = temp79a + 2.* &
          (ri3_79*bzt79(:,OP_DZZ)*bzt79(:,OP_DR) &
          -ri3_79*bzt79(:,OP_DRZ)*bzt79(:,OP_DZ) &
          +ri4_79*bzt79(:,OP_DZ )*bzt79(:,OP_DZ))
  endif

  if(ihypeta.eq.1) then
     temp = hypi*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = hypi*intx2(e(:,:,OP_1),temp79a)
  endif

  temp79a = &
       bzt79(:,OP_DZ)*bzt79(:,OP_DZZ)*ni79(:,OP_DZ) &
       + bzt79(:,OP_DR)*bzt79(:,OP_DRZ)*ni79(:,OP_DZ) &
       + bzt79(:,OP_DZ)*bzt79(:,OP_DRZ)*ni79(:,OP_DR) &
       + bzt79(:,OP_DR)*bzt79(:,OP_DRR)*ni79(:,OP_DR)
  if(itor.eq.1) then
     temp79a = temp79a - ri_79*ni79(:,OP_DR)* &
          (bzt79(:,OP_DZ)**2 + bzt79(:,OP_DR)**2)
  endif
  temp79a = temp79a * ri2_79*nt79(:,OP_1)
     
  if(ihypeta.eq.1) then
     temp = temp + hypi*intx3(e(:,:,OP_1),eta79(:,OP_1),temp79a)
  else
     temp = temp + hypi*intx2(e(:,:,OP_1),temp79a)
  endif

  qbbeta = temp
end function qbbeta


! ==============================================================
! Viscous heating terms
! ==============================================================

! quumu
! =====
function quumu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: quumu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
        
  ! Not yet implemented

  quumu = 0.
end function quumu


! quchimu
! =======
function quchimu(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: quchimu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i

  quchimu = 0.
end function quchimu


! qvvmu
! =====
function qvvmu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qvvmu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        temp = 0.
     else
        temp = - &
             (intx5(e(:,:,OP_1),r2_79,f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             +intx5(e(:,:,OP_1),r2_79,f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))
        
        if(hypv.ne.0.) then
           temp = temp - &
                hypv*intx5(e(:,:,OP_1),r2_79,f(:,OP_GS),g(:,OP_GS),h(:,OP_1))
        endif
     end if

  qvvmu = temp
end function qvvmu


! qchichimu
! =========
function qchichimu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: qchichimu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h

  ! Not yet implemented

  qchichimu = 0.
end function qchichimu


!======================================================================
! ENERGY
!======================================================================

#ifdef USECOMPLEX
#define CONJUGATE(x) conjg(x)
#else
#define CONJUGATE(x) x
#endif

! Poloidal magnetic
! -----------------
real function energy_mp(mask)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), optional, dimension(MAX_PTS) :: mask

  vectype :: temp

  if(present(mask)) then
     temp79a = mask
  else
     temp79a = 1.
  end if


  if(linear.eq.1) then
     temp = .5* &
          (int4(ri2_79,ps179(:,OP_DZ),CONJUGATE(ps179(:,OP_DZ)),temp79a) &
          +int4(ri2_79,ps179(:,OP_DR),CONJUGATE(ps179(:,OP_DR)),temp79a))
#ifdef USECOMPLEX
     if(numvar.gt.1) then
        temp = temp + .5* &
             (int3(bfp179(:,OP_DZ),CONJUGATE(bfp179(:,OP_DZ)),temp79a) &
             +int3(bfp179(:,OP_DR),CONJUGATE(bfp179(:,OP_DR)),temp79a) &
             +int4(ri_79,ps179(:,OP_DZ),CONJUGATE(bfp179(:,OP_DR)),temp79a) &
             -int4(ri_79,ps179(:,OP_DR),CONJUGATE(bfp179(:,OP_DZ)),temp79a) &
             +int4(ri_79,CONJUGATE(ps179(:,OP_DZ)),bfp179(:,OP_DR),temp79a) &
             -int4(ri_79,CONJUGATE(ps179(:,OP_DR)),bfp179(:,OP_DZ),temp79a))
     endif
#endif
  else
!    nonlinear:  do not subtract off equilibrium piece
     temp = .5* &
          (int4(ri2_79,pst79(:,OP_DZ),pst79(:,OP_DZ),temp79a) &
          +int4(ri2_79,pst79(:,OP_DR),pst79(:,OP_DR),temp79a)) ! &
!          - .5* &
!          (int4(ri2_79,ps079(:,OP_DZ),ps079(:,OP_DZ),temp79a) &
!          +int4(ri2_79,ps079(:,OP_DR),ps079(:,OP_DR),temp79a))
#if defined(USE3D)
     if(numvar.gt.1) then
        temp = temp   &
             + .5* &
             (int3(bfpt79(:,OP_DZ),bfpt79(:,OP_DZ),temp79a) &
             +int3(bfpt79(:,OP_DR),bfpt79(:,OP_DR),temp79a) &
             +2.*int4(ri_79,pst79(:,OP_DZ),bfpt79(:,OP_DR),temp79a) &
             -2.*int4(ri_79,pst79(:,OP_DR),bfpt79(:,OP_DZ),temp79a) )
     endif
#endif

#ifdef USECOMPLEX
     if(numvar.gt.1) then
        temp = temp + .5* &
             (int3(bfpt79(:,OP_DZ),CONJUGATE(bfpt79(:,OP_DZ)),temp79a) &
             +int3(bfpt79(:,OP_DR),CONJUGATE(bfpt79(:,OP_DR)),temp79a) &
             +int4(ri_79,pst79(:,OP_DZ),CONJUGATE(bfpt79(:,OP_DR)),temp79a) &
             -int4(ri_79,pst79(:,OP_DR),CONJUGATE(bfpt79(:,OP_DR)),temp79a) &
             +int4(ri_79,CONJUGATE(pst79(:,OP_DZ)),bfpt79(:,OP_DR),temp79a) &
             -int4(ri_79,CONJUGATE(pst79(:,OP_DR)),bfpt79(:,OP_DZ),temp79a))
     endif
#endif
  endif

  energy_mp = temp
  return
end function energy_mp


! Toroidal magnetic
! -----------------
real function energy_mt()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = .5*int3(ri2_79,bz179(:,OP_1),CONJUGATE(bz179(:,OP_1)))
  else
!....nonlinear:  do not subtract off equilibrium piece
     temp = .5*int3(ri2_79,bzt79(:,OP_1),bzt79(:,OP_1))!   &
!          - .5*int3(ri2_79,bz079(:,OP_1),bz079(:,OP_1))
  endif

  energy_mt = temp
  return
end function energy_mt


! Pressure
! --------
real function energy_p(mask)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), optional, dimension(MAX_PTS) :: mask

  vectype :: temp

  if(present(mask)) then
     temp79a = mask
  else
     temp79a = 1.
  end if

  if(gam.le.1.) then 
     temp = 0.
  else
     if(linear.eq.1) then
        temp = int2(p179,temp79a) / (gam - 1.)
     else
!.......nonlinear: subtract off equilibrium piece
        !temp = (int2(pt79,temp79a) - int2(p079,temp79a))/ (gam - 1.)
        temp = int2(pt79,temp79a) / (gam - 1.)
     endif
  endif

  energy_p = temp
  return
end function energy_p

! Electron Pressure
! -----------------
real function energy_pe()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gam.le.1.) then 
     temp = 0.
  else
     if(linear.eq.1) then
        temp = int1(pe179) / (gam - 1.)
     else
!.......nonlinear: subtract off equilibrium piece
        temp = int1(pet79) / (gam - 1.)
     endif
  endif

  energy_pe = temp
  return
end function energy_pe


! Poloidal kinetic
! ----------------
real function energy_kp()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = .5* &
             (int4(r2_79,ph179(:,OP_DZ),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
             +int4(r2_79,ph179(:,OP_DR),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)))
     else
        temp = .5* &
             (int4(r2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),rho79(:,OP_1)) &
             +int4(r2_79,pht79(:,OP_DR),pht79(:,OP_DR),rho79(:,OP_1)))
     endif

  energy_kp = temp
  return
end function energy_kp


! Toroidal kinetic
! ----------------
real function energy_kt()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = .5*int4(r2_79,vz179(:,OP_1),CONJUGATE(vz179(:,OP_1)),n079(:,OP_1))
     else
        temp = .5*int4(r2_79,vzt79(:,OP_1),vzt79(:,OP_1),rho79(:,OP_1))
     endif

  energy_kt = temp
  return
end function energy_kt


! Compressional kinetic
! ---------------------
real function energy_k3()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = .5* &
          (int4(ri4_79,ch179(:,OP_DZ),CONJUGATE(ch179(:,OP_DZ)),n079(:,OP_1)) &
          +int4(ri4_79,ch179(:,OP_DR),CONJUGATE(ch179(:,OP_DR)),n079(:,OP_1)) &
          +int4(ri_79,ch179(:,OP_DZ),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)) &
          -int4(ri_79,ch179(:,OP_DR),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
          +int4(ri_79,CONJUGATE(ch179(:,OP_DZ)),ph179(:,OP_DR),n079(:,OP_1)) &
          -int4(ri_79,CONJUGATE(ch179(:,OP_DR)),ph179(:,OP_DZ),n079(:,OP_1)))
     else
        temp = .5* &
          (int4(ri4_79,cht79(:,OP_DZ),cht79(:,OP_DZ),rho79(:,OP_1)) &
          +int4(ri4_79,cht79(:,OP_DR),cht79(:,OP_DR),rho79(:,OP_1)) &
          +int4(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),rho79(:,OP_1)) &
          -int4(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),rho79(:,OP_1)) &
          +int4(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),rho79(:,OP_1)) &
          -int4(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),rho79(:,OP_1)))
     endif

  energy_k3 = temp
  return
end function energy_k3


! Poloidal resistive
! ------------------
real function energy_mpd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = -int4(ri2_79,ps179(:,OP_GS),CONJUGATE(ps179(:,OP_GS)),eta79(:,OP_1))
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,ps179(:,OP_DZP),CONJUGATE(ps179(:,OP_DZP)),eta79(:,OP_1)) &
          +int4(ri4_79,ps179(:,OP_DRP),CONJUGATE(ps179(:,OP_DRP)),eta79(:,OP_1)))
#endif
  else
     temp = -int4(ri2_79,pst79(:,OP_GS),CONJUGATE(pst79(:,OP_GS)),eta79(:,OP_1))
  endif

  energy_mpd = temp
  return
end function energy_mpd


! Toroidal resistive
! ------------------
real function energy_mtd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = - &
          (int4(ri2_79,bz179(:,OP_DZ),CONJUGATE(bz179(:,OP_DZ)),eta79(:,OP_1))&
          +int4(ri2_79,bz179(:,OP_DR),CONJUGATE(bz179(:,OP_DR)),eta79(:,OP_1)))
  else
     temp = - &
          (int4(ri2_79,bzt79(:,OP_DZ),CONJUGATE(bzt79(:,OP_DZ)),eta79(:,OP_1))&
          +int4(ri2_79,bzt79(:,OP_DR),CONJUGATE(bzt79(:,OP_DR)),eta79(:,OP_1)))
  end if

  energy_mtd = temp
  return
end function energy_mtd


! Poloidal viscous
! ----------------
real function energy_kpd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = -int4(ri2_79,ph179(:,OP_GS),CONJUGATE(ph179(:,OP_GS)),vis79(:,OP_1))
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,ph179(:,OP_DZP),CONJUGATE(ph179(:,OP_DZP)),vis79(:,OP_1)) &
          +int4(ri4_79,ph179(:,OP_DRP),CONJUGATE(ph179(:,OP_DRP)),vis79(:,OP_1)))
#endif
  else
     temp = -int4(ri2_79,pht79(:,OP_GS),CONJUGATE(pht79(:,OP_GS)),vis79(:,OP_1))
  endif

  energy_kpd = temp
  return
end function energy_kpd


! Toroidal viscous
! ----------------
real function energy_ktd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     if(linear.eq.1) then
        temp = - &
             (int4(r2_79,vz179(:,OP_DZ),CONJUGATE(vz179(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(r2_79,vz179(:,OP_DR),CONJUGATE(vz179(:,OP_DR)),vis79(:,OP_1)))
     else
        temp = - &
             (int4(r2_79,vzt79(:,OP_DZ),CONJUGATE(vzt79(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(r2_79,vzt79(:,OP_DR),CONJUGATE(vzt79(:,OP_DR)),vis79(:,OP_1)))
     endif

  energy_ktd = temp
  return
end function energy_ktd

! Compressional viscous
! ---------------------
real function energy_k3d()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = - 2.*int3(ch179(:,OP_LP),CONJUGATE(ch179(:,OP_LP)),vic79(:,OP_1))
  else
     temp = - 2.*int3(cht79(:,OP_LP),CONJUGATE(cht79(:,OP_LP)),vic79(:,OP_1))
  endif

  energy_k3d = temp
  return
end function energy_k3d


! Poloidal hyper-viscous
! ----------------------
real function energy_kph()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

!!$  if(hypc.ne.0.) then
!!$     temp = - hypc* &
!!$          (int4(ri2_79,vot79(:,OP_DZ),CONJUGATE(vot79(:,OP_DZ)),vis79(:,OP_1)) &
!!$          +int4(ri2_79,vot79(:,OP_DR),CONJUGATE(vot79(:,OP_DR)),vis79(:,OP_1)))
!!$  else
     temp = 0.
!!$  end if
  energy_kph = temp
  return
end function energy_kph


! Toroidal hyper-viscous
! ----------------------
real function energy_kth()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = -hypv*int4(r2_79,vzt79(:,OP_GS),CONJUGATE(vzt79(:,OP_GS)),vis79(:,OP_1))

  energy_kth = temp
  return
end function energy_kth

! Compressional hyper-viscous
! ---------------------------
real function energy_k3h()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

!!$  if(hypc.ne.0.) then
!!$     temp = -2.*hypc* &
!!$          (int3(cot79(:,OP_DZ),CONJUGATE(cot79(:,OP_DZ)),vic79(:,OP_1)) &
!!$          +int3(cot79(:,OP_DR),CONJUGATE(cot79(:,OP_DR)),vic79(:,OP_1)))
!!$  else
     temp = 0.
!!$  end if
  energy_k3h = temp
  return
end function energy_k3h



!======================================================================
! FLUXES
!======================================================================


! Pressure convection
! -------------------
real function flux_pressure()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if((numvar.lt.3 .and. ipres.eq.0) .or. gam.eq.1. .or. gam.eq.0.) then
     flux_pressure = 0.
     return
  endif

     temp = 0.5* &
          (int4(r_79,pt79(:,OP_1),norm79(:,2),pht79(:,OP_DR)) &
          -int4(r_79,pt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
          +int4(ri2_79,pt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
          +int4(ri2_79,pt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ)))

  if(db .ne. 0.) then
     temp = temp + db* &
          (int5(ri_79,pet79(:,OP_1),ni79(:,OP_1),norm79(:,1),bzt79(:,OP_DZ)) &
          -int5(ri_79,pet79(:,OP_1),ni79(:,OP_1),norm79(:,2),bzt79(:,OP_DR)))
  endif

  flux_pressure = gam*real(temp)/(gam-1.)

  return
end function flux_pressure


! Kinetic Energy Convection
! -------------------------
real function flux_ke()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

     temp79a = r2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2) &
          +    r2_79*vzt79(:,OP_1)**2 &
          +   ri4_79*(cht79(:,OP_DZ)**2 + cht79(:,OP_DR)**2) &
          + 2.*ri_79* &
          (cht79(:,OP_DZ)*pht79(:,OP_DR)-cht79(:,OP_DR)*pht79(:,OP_DZ))

     temp = 0.5* &
          (int5(r_79,nt79(:,OP_1),norm79(:,2),pht79(:,OP_DR),temp79a) &
          -int5(r_79,nt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ),temp79a) &
          +int5(ri2_79,nt79(:,OP_1),norm79(:,1),cht79(:,OP_DR),temp79a) &
          +int5(ri2_79,nt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ),temp79a))

  flux_ke = real(temp)
  return
end function flux_ke


! Poynting flux
! -------------
real function flux_poynting()

  use math
  use basic
  use m3dc1_nint

  implicit none
  vectype :: temp

  temp = -vloop/toroidal_period * &
       (int3(ri2_79,norm79(:,1),pst79(:,OP_DR)) &
       +int3(ri2_79,norm79(:,2),pst79(:,OP_DZ)))

  flux_poynting = real(temp)
  return
end function flux_poynting


! Heat flux
! ---------
real function flux_heat()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(numvar.lt.3 .and. ipres.eq.0) then
     flux_heat = 0.
     return
  endif

!!$  temp = int4(kap79(:,OP_1),norm79(:,1),pt79(:,OP_DR),ni79(:,OP_1)) &
!!$       + int4(kap79(:,OP_1),norm79(:,2),pt79(:,OP_DZ),ni79(:,OP_1)) &
!!$       + int4(kap79(:,OP_1),norm79(:,1),pt79(:,OP_1),ni79(:,OP_DR)) &
!!$       + int4(kap79(:,OP_1),norm79(:,2),pt79(:,OP_1),ni79(:,OP_DZ))
!!$
!!$  if(kappar.ne.0.) then
!!$     temp79a = ni79(:,OP_1)* &
!!$          (pt79(:,OP_DZ)*pst79(:,OP_DR) - pt79(:,OP_DR)*pst79(:,OP_DZ)) &
!!$          +    pt79(:,OP_1)* &
!!$          (ni79(:,OP_DZ)*pst79(:,OP_DR) - ni79(:,OP_DR)*pst79(:,OP_DZ))
!!$     temp79b = norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR)
!!$     temp = temp &
!!$          + int5(ri2_79,kar79(:,OP_1),b2i79(:,OP_1),temp79a,temp79b)
!!$  endif

  temp = -int3(kap79(:,OP_1),norm79(:,1),tet79(:,OP_DR)) &
       - int3(kap79(:,OP_1),norm79(:,2),tet79(:,OP_DZ))
  temp = temp &
       - kappai_fac*(  int3(kap79(:,OP_1),norm79(:,1),tit79(:,OP_DR)) &
                     + int3(kap79(:,OP_1),norm79(:,2),tit79(:,OP_DZ)))

  if(kappar.ne.0.) then
     temp79a = (tet79(:,OP_DZ)*pst79(:,OP_DR)-tet79(:,OP_DR)*pst79(:,OP_DZ))
     temp79c = (tit79(:,OP_DZ)*pst79(:,OP_DR)-tit79(:,OP_DR)*pst79(:,OP_DZ))
     temp79b = norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR)
     temp = temp &
          + int5(ri2_79,kar79(:,OP_1),b2i79(:,OP_1),temp79a,temp79b) 
     temp = temp &
          + int5(ri2_79,kar79(:,OP_1),b2i79(:,OP_1),temp79c,temp79b)
  endif

  flux_heat = real(temp)
  return
end function flux_heat


! Grav_pot
! --------
real function grav_pot()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     grav_pot = 0.
     return
  endif

  temp = gravr*int3(ri3_79,pht79(:,OP_DZ),nt79(:,OP_1)) &
       - gravz*int3( ri_79,pht79(:,OP_DR),nt79(:,OP_1))
     
  if(numvar.ge.3) then
     temp = -gravr*int3(ri2_79,cht79(:,OP_DR),nt79(:,OP_1)) &
          -gravz*int2(       cht79(:,OP_DZ),nt79(:,OP_1))
  endif

  grav_pot = real(temp)
  return
end function grav_pot


!======================================================================
! Toroidal (angular) momentum
!======================================================================

! torque_em
! ~~~~~~~~~
vectype function torque_em()
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = int4(ri_79,bzt79(:,OP_1),norm79(:,2),pst79(:,OP_DR)) &
       - int4(ri_79,bzt79(:,OP_1),norm79(:,1),pst79(:,OP_DZ))
  
  torque_em = temp
end function torque_em

! torque_sol
! ~~~~~~~~~~
vectype function torque_sol()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = int5(r3_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
       - int5(r3_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),pht79(:,OP_DR))
  
  torque_sol = temp
end function torque_sol

! torque_com
! ~~~~~~~~~~
vectype function torque_com()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(numvar.lt.3) then
     torque_com = 0.
     return
  end if

  temp = &
       - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
       - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ))     

  torque_com = temp
end function torque_com


! torque_visc
! ~~~~~~~~~~~
vectype function torque_visc()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = int4(r2_79,vis79(:,OP_1),norm79(:,1),vzt79(:,OP_DR)) &
       + int4(r2_79,vis79(:,OP_1),norm79(:,2),vzt79(:,OP_DZ))

  torque_visc = temp
end function torque_visc


! torque_parvisc
! ~~~~~~~~~~~~~~
vectype function torque_parvisc()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp
  
  if(amupar.eq.0.) then
     torque_parvisc = 0.
     return
  endif
       
  call PVS1(pht79,temp79a)

  if(numvar.ge.2) then
     call PVS2(vzt79,temp79b)
     temp79a = temp79a + temp79b
  endif
  
  if(numvar.ge.3) then
     call PVS3(cht79,temp79c)
     temp79a = temp79a + temp79c
  endif

  temp79d = 3.*vip79(:,OP_1)*bzt79(:,OP_1)*b2i79(:,OP_1)*temp79a
  temp = int4(ri_79,temp79d,norm79(:,2),pst79(:,OP_DR)) &
       - int4(ri_79,temp79d,norm79(:,1),pst79(:,OP_DZ))

  torque_parvisc = temp
end function torque_parvisc



! torque_gyro
! ~~~~~~~~~~~
vectype function torque_gyro()
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gyro.eq.0 .or. numvar.lt.2) then
     torque_gyro = 0.
     return
  endif

  temp79a = 0.25*db*b2i79(:,OP_1)
  temp79b = temp79a*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)

     temp79c = norm79(:,1)*pst79(:,OP_DZ) + norm79(:,2)*pst79(:,OP_DR)
     temp79d = norm79(:,1)*pst79(:,OP_DR) - norm79(:,2)*pst79(:,OP_DZ)
     temp79e = 3.*temp79a*b2i79(:,OP_1)* &
          (norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR))
     temp79f = pht79(:,OP_DZZ) - pht79(:,OP_DRR)
     if(itor.eq.1) then
        temp79f = temp79f - ri_79*pht79(:,OP_DR)
     endif

     ! U contribution
     temp = int4(r_79,temp79b,temp79c,temp79f) &
        +2.*int4(r_79,temp79b,temp79d,pht79(:,OP_DRZ)) &
        +   int5(ri_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),temp79f) &
        -   int5(ri_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),temp79f) &
        +4.*int5(ri_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),pht79(:,OP_DRZ))

     if(itor.eq.1) then
        temp = temp &
           + 2.*int4(temp79b,norm79(:,1),pst79(:,OP_DR),pht79(:,OP_DZ)) &
           - 2.*int4(temp79a,pht79(:,OP_DZ),norm79(:,1),pst79(:,OP_DR)) &
           - 2.*int4(temp79a,pht79(:,OP_DZ),norm79(:,2),pst79(:,OP_DZ)) &
           + 2.*int5(ri2_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),pht79(:,OP_DZ))
     endif
     
     ! omega contribution
     temp = temp + 2.* &
          (int5(r_79,temp79b,bzt79(:,OP_1),norm79(:,1),vzt79(:,OP_DZ)) &
          -int5(r_79,temp79b,bzt79(:,OP_1),norm79(:,2),vzt79(:,OP_DR)))

     temp79f = cht79(:,OP_DZZ) - cht79(:,OP_DRR)
     if(itor.eq.1) then
        temp79f = temp79f + 2.*ri_79*cht79(:,OP_DR)
     endif

     ! chi constribution
     temp = temp &
          +    int4(ri2_79,temp79b,temp79d,temp79f) &
          - 2.*int4(ri2_79,temp79b,temp79c,cht79(:,OP_DRZ)) &
          +    int5(ri2_79,temp79b,cht79(:,OP_GS),norm79(:,1),pst79(:,OP_DR)) &
          +    int5(ri2_79,temp79b,cht79(:,OP_GS),norm79(:,2),pst79(:,OP_DZ)) &
          - 2.*int4(ri2_79,cht79(:,OP_GS),norm79(:,1),pst79(:,OP_DR)) &
          - 2.*int4(ri2_79,cht79(:,OP_GS),norm79(:,2),pst79(:,OP_DZ)) &
          + 2.*int5(ri4_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),temp79f) &
          - 2.*int5(ri4_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),cht79(:,OP_DRZ)) &
          + 2.*int5(ri4_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),cht79(:,OP_DRZ))

     if(itor.eq.1) then
        temp = temp &
             + 2.*int4(ri3_79,temp79b,temp79c,cht79(:,OP_DZ)) &
             + 3.*int5(ri3_79,temp79b,cht79(:,OP_DR),norm79(:,1),pst79(:,OP_DR)) &
             + 3.*int5(ri3_79,temp79b,cht79(:,OP_DR),norm79(:,2),pst79(:,OP_DZ)) &
             - 6.*int4(ri3_79,cht79(:,OP_DR),norm79(:,1),pst79(:,OP_DR)) &
             - 6.*int4(ri3_79,cht79(:,OP_DR),norm79(:,2),pst79(:,OP_DZ)) &
             + 2.*int5(ri5_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),cht79(:,OP_DZ)) &
             - 2.*int5(ri5_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),cht79(:,OP_DZ))
     endif

  torque_gyro = 0.
end function torque_gyro


! torque_denm
! ~~~~~~~~~~~
vectype function torque_denm()
  use basic
  use m3dc1_nint

  implicit none

  torque_denm = 0.

end function torque_denm





function tepsipsikappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsipsikappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsipsikappar = 0.
     return
  end if

  if(surface_int) then
     temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))
     temp = intx4(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
          - intx4(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DZ))

  else
     temp79a = k(:,OP_1)*ri2_79*j(:,OP_1)

     temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DZ))
  end if

  tepsipsikappar = (gam - 1.) * temp
end function tepsipsikappar

function tepsipsikappar1(g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsipsikappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,j,k
  type(prodarray) :: temp
  type(muarray) :: tempa

  if(gam.le.1.) then
     tepsipsikappar1%len = 0
     return
  end if

  if(surface_int) then
     tempa = (mu(norm79(:,1),OP_DZ) + mu(-norm79(:,2),OP_DR)) &
           * ri2_79*k(:,OP_1)*j(:,OP_1)
     temp = prod(mu( g(:,OP_DZ)*h(:,OP_DR),OP_1),tempa) &
          + prod(mu(-g(:,OP_DR)*h(:,OP_DZ),OP_1),tempa)

  else
     temp79a = k(:,OP_1)*ri2_79*j(:,OP_1)

     temp = prod( temp79a*g(:,OP_DZ)*h(:,OP_DR),OP_DZ,OP_DR) &
          + prod(-temp79a*g(:,OP_DZ)*h(:,OP_DR),OP_DR,OP_DZ) &
          + prod(-temp79a*g(:,OP_DR)*h(:,OP_DZ),OP_DZ,OP_DR) &
          + prod( temp79a*g(:,OP_DR)*h(:,OP_DZ),OP_DR,OP_DZ)
  end if

  tepsipsikappar1 = (gam - 1.) * temp
end function tepsipsikappar1

function tepsipsikappar2(f,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsipsikappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,j,k
  type(prodarray) :: temp

  if(gam.le.1.) then
     tepsipsikappar2%len = 0
     return
  end if

  if(surface_int) then
     temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))
     temp = prod( temp79a*h(:,OP_DR),OP_1,OP_DZ) &
          + prod(-temp79a*h(:,OP_DZ),OP_1,OP_DR)

  else
     temp79a = k(:,OP_1)*ri2_79*j(:,OP_1)

     temp = prod( f(:,OP_DR)*temp79a*h(:,OP_DR),OP_DZ,OP_DZ) &
          + prod(-f(:,OP_DZ)*temp79a*h(:,OP_DR),OP_DR,OP_DZ) &
          + prod(-f(:,OP_DR)*temp79a*h(:,OP_DZ),OP_DZ,OP_DR) &
          + prod( f(:,OP_DZ)*temp79a*h(:,OP_DZ),OP_DR,OP_DR)
  end if

  tepsipsikappar2 = (gam - 1.) * temp
end function tepsipsikappar2

function tepsipsikappar3(f,g,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsipsikappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,j,k
  type(prodarray) :: temp

  if(gam.le.1.) then
     tepsipsikappar3%len = 0
     return
  end if

  if(surface_int) then
     temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))
     temp = prod( temp79a*g(:,OP_DZ),OP_1,OP_DR) &
          + prod(-temp79a*g(:,OP_DR),OP_1,OP_DZ)

  else
     temp79a = k(:,OP_1)*ri2_79*j(:,OP_1)

     temp = prod( f(:,OP_DR)*temp79a*g(:,OP_DZ),OP_DZ,OP_DR) &
          + prod(-f(:,OP_DZ)*temp79a*g(:,OP_DZ),OP_DR,OP_DR) &
          + prod(-f(:,OP_DR)*temp79a*g(:,OP_DR),OP_DZ,OP_DZ) &
          + prod( f(:,OP_DZ)*temp79a*g(:,OP_DR),OP_DR,OP_DZ)
  end if

  tepsipsikappar3 = (gam - 1.) * temp
end function tepsipsikappar3

function tepsibkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsibkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsibkappar = 0.
     return
  end if

  if(surface_int) then
     temp79a = -ri3_79*k(:,OP_1)*j(:,OP_1)*g(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))

     temp = intx3(e(:,:,OP_1),temp79a,h(:,OP_DP))
  else
     temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)*j(:,OP_1)  

     temp79b = f(:,OP_DR)*(h(:,OP_DZ) ) &
          -    f(:,OP_DZ)*(h(:,OP_DR) )

     temp79d = temp79b*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1)

     temp = intx4(e(:,:,OP_DZ),f(:,OP_DR),temp79a,h(:,OP_DP)) &
          - intx4(e(:,:,OP_DR),f(:,OP_DZ),temp79a,h(:,OP_DP)) &
          - intx3(e(:,:,OP_DP),ri3_79,temp79d)
  end if
  tepsibkappar = (gam - 1.) * temp
#else
  tepsibkappar = 0.
#endif
end function tepsibkappar

function tepsibkappar1(g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsibkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempd

  if(gam.le.1.) then
     tepsibkappar1%len = 0
     return
  end if

  if(surface_int) then
     tempa = (mu(norm79(:,1),OP_DZ) + mu(-norm79(:,2),OP_DR)) &
           * (-ri3_79*k(:,OP_1)*j(:,OP_1)*g(:,OP_1))
     temp = prod(mu(h(:,OP_DP),OP_1),tempa)
  else
     temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)*j(:,OP_1)  

     tempb = mu( h(:,OP_DZ),OP_DR) &
        +    mu(-h(:,OP_DR),OP_DZ)

     tempd = tempb*(g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1))

     temp = prod( temp79a*h(:,OP_DP),OP_DZ,OP_DR) &
          + prod(-temp79a*h(:,OP_DP),OP_DR,OP_DZ) &
          + prod(mu(-ri3_79,OP_DP),tempd)
  end if
  tepsibkappar1 = (gam - 1.) * temp
#else
  tepsibkappar1%len = 0
#endif
end function tepsibkappar1

function tepsibkappar2(f,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsibkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempd

  if(gam.le.1.) then
     tepsibkappar2%len = 0
     return
  end if

  if(surface_int) then
     tempa = mu(-ri3_79*k(:,OP_1)*j(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR)),OP_1)

     temp = prod(mu(h(:,OP_DP),OP_1),tempa)
  else
     tempa = mu(-k(:,OP_1)*ri3_79*j(:,OP_1),OP_1)  

     temp79b = f(:,OP_DR)*(h(:,OP_DZ) ) &
          -    f(:,OP_DZ)*(h(:,OP_DR) )

     tempd = mu(temp79b*j(:,OP_1 )*k(:,OP_1),OP_1)

     temp = prod(mu( f(:,OP_DR)*h(:,OP_DP),OP_DZ),tempa) &
          + prod(mu(-f(:,OP_DZ)*h(:,OP_DP),OP_DR),tempa) &
          + prod(mu(-ri3_79,OP_DP),tempd)
  end if
  tepsibkappar2 = (gam - 1.) * temp
#else
  tepsibkappar2%len = 0
#endif
end function tepsibkappar2

function tepsibkappar3(f,g,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsibkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempd

  if(gam.le.1.) then
     tepsibkappar3%len = 0
     return
  end if

  if(surface_int) then
     temp79a = -ri3_79*k(:,OP_1)*j(:,OP_1)*g(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))

     temp = prod(temp79a,OP_1,OP_DP)
  else
     temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)*j(:,OP_1)  

     tempb = mu( f(:,OP_DR),OP_DZ) &
        +    mu(-f(:,OP_DZ),OP_DR)

     tempd = tempb*(g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1))

     temp = prod( temp79a*f(:,OP_DR),OP_DZ,OP_DP) &
          + prod(-temp79a*f(:,OP_DZ),OP_DR,OP_DP) &
          + prod(mu(-ri3_79,OP_DP),tempd)
  end if
  tepsibkappar3 = (gam - 1.) * temp
#else
  tepsibkappar3%len = 0
#endif
end function tepsibkappar3

function tepsibkapparl(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsibkapparl
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsibkapparl = 0.
     return
  end if

  temp79a = i(:,OP_1)*j(:,OP_1)*g(:,OP_1)

  if(surface_int) then
     ! this can't be right .. no trial function here
!     temp = int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,2),f(:,OP_DR)) &
!          - int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,1),f(:,OP_DZ))
     temp = 0.
  else
     ! d(temp79a)/dphi
     temp79b = i(:,OP_DP)*j(:,OP_1 )*g(:,OP_1 ) &
          +    i(:,OP_1 )*j(:,OP_DP)*g(:,OP_1 ) &
          +    i(:,OP_1 )*j(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]
     temp79c = (h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ))

     ! d(temp79c)/dphi
     temp79d =      &
          (h(:,OP_DZ)*f(:,OP_DRP)-h(:,OP_DR)*f(:,OP_DZP)) &
        + (h(:,OP_DZP)*f(:,OP_DR)-h(:,OP_DRP)*f(:,OP_DZ))

     temp = intx4(e(:,:,OP_1),ri3_79,temp79a,temp79d) &
          + intx4(e(:,:,OP_1),ri3_79,temp79b,temp79c) &
          + intx5(e(:,:,OP_DR),ri3_79,temp79a,f(:,OP_DZ),h(:,OP_DP)) &
          - intx5(e(:,:,OP_DZ),ri3_79,temp79a,f(:,OP_DR),h(:,OP_DP))
  end if
  tepsibkapparl = (gam - 1.) * temp
#else
  tepsibkapparl = 0.
#endif
end function tepsibkapparl

function tepsibkapparl1(g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsibkapparl1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc, tempd

  if(gam.le.1.) then
     tepsibkapparl1%len = 0
     return
  end if

  temp79a = i(:,OP_1)*j(:,OP_1)*g(:,OP_1)

  if(surface_int) then
     ! this can't be right .. no trial function here
!     temp = int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,2),f(:,OP_DR)) &
!          - int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,1),f(:,OP_DZ))
     temp%len = 0
  else
     ! d(temp79a)/dphi
     temp79b = i(:,OP_DP)*j(:,OP_1 )*g(:,OP_1 ) &
          +    i(:,OP_1 )*j(:,OP_DP)*g(:,OP_1 ) &
          +    i(:,OP_1 )*j(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]
     tempc = mu(h(:,OP_DZ),OP_DR)+mu(-h(:,OP_DR),OP_DZ)

     ! d(temp79c)/dphi
     tempd =      &
          mu(h(:,OP_DZ),OP_DRP)+mu(-h(:,OP_DR),OP_DZP) &
        + mu(h(:,OP_DZP),OP_DR)+mu(-h(:,OP_DRP),OP_DZ)


     temp = prod(mu(ri3_79*temp79a,OP_1),tempd) &
          + prod(mu(ri3_79*temp79b,OP_1),tempc) &
          + prod( ri3_79*temp79a*h(:,OP_DP),OP_DR,OP_DZ) &
          + prod(-ri3_79*temp79a*h(:,OP_DP),OP_DZ,OP_DR)
  end if
  tepsibkapparl1 = (gam - 1.) * temp
#else
  tepsibkapparl1%len = 0
#endif
end function tepsibkapparl1

function tepsibkapparl2(f,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsibkapparl2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb

  if(gam.le.1.) then
     tepsibkapparl2%len = 0
     return
  end if

  tempa = mu(i(:,OP_1)*j(:,OP_1),OP_1)

  if(surface_int) then
     ! this can't be right .. no trial function here
!     temp = int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,2),f(:,OP_DR)) &
!          - int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,1),f(:,OP_DZ))
     temp%len = 0
  else
     ! d(temp79a)/dphi
     tempb = mu(i(:,OP_DP)*j(:,OP_1 ),OP_1 ) &
          +    mu(i(:,OP_1 )*j(:,OP_DP),OP_1 ) &
          +    mu(i(:,OP_1 )*j(:,OP_1 ),OP_DP)


     ! [T,psi]
     temp79c = (h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ))

     ! d(temp79c)/dphi
     temp79d =      &
          (h(:,OP_DZ)*f(:,OP_DRP)-h(:,OP_DR)*f(:,OP_DZP)) &
        + (h(:,OP_DZP)*f(:,OP_DR)-h(:,OP_DRP)*f(:,OP_DZ))

     temp = prod(mu(ri3_79*temp79d,OP_1),tempa) &
          + prod(mu(ri3_79*temp79c,OP_1),tempb) &
          + prod(mu( ri3_79*f(:,OP_DZ)*h(:,OP_DP),OP_DR),tempa) &
          + prod(mu(-ri3_79*f(:,OP_DR)*h(:,OP_DP),OP_DZ),tempa)
  end if
  tepsibkapparl2 = (gam - 1.) * temp
#else
  tepsibkapparl2%len = 0
#endif
end function tepsibkapparl2

function tepsibkapparl3(f,g,i,j)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsibkapparl3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc, tempd

  if(gam.le.1.) then
     tepsibkapparl3%len = 0
     return
  end if

  temp79a = i(:,OP_1)*j(:,OP_1)*g(:,OP_1)

  if(surface_int) then
     ! this can't be right .. no trial function here
!     temp = int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,2),f(:,OP_DR)) &
!          - int5(ri3_79,temp79a,h(:,OP_DP),norm79(:,1),f(:,OP_DZ))
     temp%len = 0
  else
     ! d(temp79a)/dphi
     temp79b = i(:,OP_DP)*j(:,OP_1 )*g(:,OP_1 ) &
          +    i(:,OP_1 )*j(:,OP_DP)*g(:,OP_1 ) &
          +    i(:,OP_1 )*j(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]
     tempc = mu(f(:,OP_DR),OP_DZ)+mu(-f(:,OP_DZ),OP_DR)

     ! d(temp79c)/dphi
     tempd =      &
          mu(f(:,OP_DRP),OP_DZ)+mu(-f(:,OP_DZP),OP_DR) &
        + mu(f(:,OP_DR),OP_DZP)+mu(-f(:,OP_DZ),OP_DRP)


     temp = prod(mu(ri3_79*temp79a,OP_1),tempd) &
          + prod(mu(ri3_79*temp79b,OP_1),tempc) &
          + prod( ri3_79*temp79a*f(:,OP_DZ),OP_DR,OP_DP) &
          + prod(-ri3_79*temp79a*f(:,OP_DR),OP_DZ,OP_DP)
  end if
  tepsibkapparl3 = (gam - 1.) * temp
#else
  tepsibkapparl3%len = 0
#endif
end function tepsibkapparl3

function tebbkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebbkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tebbkappar = 0.
     return
  end if

  if(surface_int) then
     temp = 0.
  else
     temp79a = h(:,OP_DP)

     temp79c = f(:,OP_1)*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 )

     temp = -intx3(e(:,:,OP_DP),ri4_79,temp79c)
  end if
  tebbkappar = (gam - 1.) * temp
#else
  tebbkappar = 0.
#endif
end function tebbkappar

function tebbkappar1(g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebbkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc

  if(gam.le.1.) then
     tebbkappar1%len = 0
     return
  end if

  if(surface_int) then
     temp%len = 0
  else
     temp79a = h(:,OP_DP)

     tempc = mu(g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 ),OP_1)

     temp = prod(mu(-ri4_79,OP_DP),tempc)
  end if
  tebbkappar1 = (gam - 1.) * temp
#else
  tebbkappar1%len = 0
#endif
end function tebbkappar1

function tebbkappar2(f,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebbkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc

  if(gam.le.1.) then
     tebbkappar2%len = 0
     return
  end if

  if(surface_int) then
     temp%len = 0
  else
     temp79a = h(:,OP_DP)

     tempc = mu(f(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 ),OP_1)

     temp = prod(mu(-ri4_79,OP_DP),tempc)
  end if
  tebbkappar2 = (gam - 1.) * temp
#else
  tebbkappar2%len = 0
#endif
end function tebbkappar2

function tebbkappar3(f,g,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebbkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempc

  if(gam.le.1.) then
     tebbkappar3%len = 0
     return
  end if

  if(surface_int) then
     temp%len = 0
  else
     temp79a = 1.
     tempa = mu(temp79a,OP_DP)

     tempc = tempa*(f(:,OP_1)*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ))

     temp = prod(mu(-ri4_79,OP_DP),tempc)
  end if
  tebbkappar3 = (gam - 1.) * temp
#else
  tebbkappar3%len = 0
#endif
end function tebbkappar3
!
!...the following function must replace tebbkappar for linear runs.
function tebbkapparl(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebbkapparl
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp
  
  if(gam.le.1.) then
     tebbkapparl = 0.
     return
  end if

  if(surface_int) then
     temp = 0.
  else
     temp79a = i(:,OP_1)*j(:,OP_1)*f(:,OP_1)*g(:,OP_1)

     ! d(temp79a)/dphi
     temp79b = i(:,OP_DP)*j(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    i(:,OP_1)*j(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
          +    i(:,OP_1)*j(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
          +    i(:,OP_1)*j(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

     temp = intx4(e(:,:,OP_1),ri4_79,temp79a,h(:,OP_DPP)) &
          + intx4(e(:,:,OP_1),ri4_79,temp79b,h(:,OP_DP))
  end if
  tebbkapparl = (gam - 1.) * temp
#else
  tebbkapparl = 0.
#endif
end function tebbkapparl

function tebbkapparl1(g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebbkapparl1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb
  
  if(gam.le.1.) then
     tebbkapparl1%len = 0
     return
  end if

  if(surface_int) then
     temp%len = 0
  else
     tempa = mu(i(:,OP_1)*j(:,OP_1)*g(:,OP_1),OP_1)

     ! d(temp79a)/dphi
     tempb = mu(i(:,OP_DP)*j(:,OP_1)*g(:,OP_1),OP_1) &
        +    mu(i(:,OP_1)*j(:,OP_DP)*g(:,OP_1),OP_1) &
        +    mu(i(:,OP_1)*j(:,OP_1)*g(:,OP_1),OP_DP) &
        +    mu(i(:,OP_1)*j(:,OP_1)*g(:,OP_DP),OP_1)


     temp = prod(mu(ri4_79*h(:,OP_DPP),OP_1),tempa) &
          + prod(mu(ri4_79*h(:,OP_DP),OP_1),tempb)
  end if
  tebbkapparl1 = (gam - 1.) * temp
#else
  tebbkapparl1%len = 0
#endif
end function tebbkapparl1

function tebbkapparl2(f,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebbkapparl2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb
  
  if(gam.le.1.) then
     tebbkapparl2%len = 0
     return
  end if

  if(surface_int) then
     temp%len = 0
  else
     tempa = mu(i(:,OP_1)*j(:,OP_1)*f(:,OP_1),OP_1)

     ! d(temp79a)/dphi
     tempb = mu(i(:,OP_DP)*j(:,OP_1)*f(:,OP_1),OP_1) &
        +    mu(i(:,OP_1)*j(:,OP_DP)*f(:,OP_1),OP_1) &
        +    mu(i(:,OP_1)*j(:,OP_1)*f(:,OP_DP),OP_1) &
        +    mu(i(:,OP_1)*j(:,OP_1)*f(:,OP_1),OP_DP)


     temp = prod(mu(ri4_79*h(:,OP_DPP),OP_1),tempa) &
          + prod(mu(ri4_79*h(:,OP_DP),OP_1),tempb)
  end if
  tebbkapparl2 = (gam - 1.) * temp
#else
  tebbkapparl2%len = 0
#endif
end function tebbkapparl2

function tebbkapparl3(f,g,i,j)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebbkapparl3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  
  if(gam.le.1.) then
     tebbkapparl3%len = 0
     return
  end if

  if(surface_int) then
     temp%len = 0
  else
     temp79a = i(:,OP_1)*j(:,OP_1)*f(:,OP_1)*g(:,OP_1)

     ! d(temp79a)/dphi
     temp79b = i(:,OP_DP)*j(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    i(:,OP_1)*j(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
          +    i(:,OP_1)*j(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
          +    i(:,OP_1)*j(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

     temp = prod(ri4_79*temp79a,OP_1,OP_DPP) &
          + prod(ri4_79*temp79b,OP_1,OP_DP)
  end if
  tebbkapparl3 = (gam - 1.) * temp
#else
  tebbkapparl3%len = 0
#endif
end function tebbkapparl3

function tepsifkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tepsifkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tepsifkappar = 0.
     return
  end if

  if(surface_int) then
     temp79a = k(:,OP_1)*ri_79* &
          (norm79(:,2)*f(:,OP_DR) - norm79(:,1)*f(:,OP_DZ))*j(:,OP_1)
     temp79b = -k(:,OP_1)*ri_79* &
          (norm79(:,2)*g(:,OP_DZ) + norm79(:,1)*g(:,OP_DR))*j(:,OP_1)

     temp = intx4(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx4(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          + intx4(e(:,:,OP_1),temp79b,f(:,OP_DR ),h(:,OP_DZ)) &
          - intx4(e(:,:,OP_1),temp79b,f(:,OP_DZ ),h(:,OP_DR))
  else
     temp79a = k(:,OP_1)*ri_79*j(:,OP_1)

     temp = intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          + intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DR ),h(:,OP_DZ)) &
          - intx5(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,f(:,OP_DZ ),h(:,OP_DR)) &
          - intx5(e(:,:,OP_DR),g(:,OP_DR),temp79a,f(:,OP_DZ ),h(:,OP_DR))
  end if
  tepsifkappar = (gam - 1.) * temp
#else
  tepsifkappar = 0.
#endif
end function tepsifkappar

function tepsifkappar1(g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsifkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

  if(gam.le.1.) then
     tepsifkappar1%len = 0
     return
  end if

  if(surface_int) then
     tempa = (mu(norm79(:,2),OP_DR) + mu(-norm79(:,1),OP_DZ)) &
           * j(:,OP_1)*(k(:,OP_1)*ri_79)
     temp79b = -k(:,OP_1)*ri_79* &
          (norm79(:,2)*g(:,OP_DZ) + norm79(:,1)*g(:,OP_DR))*j(:,OP_1)

     temp = prod(mu(g(:,OP_DZ )*h(:,OP_DZ),OP_1),tempa) &
          + prod(mu(g(:,OP_DR )*h(:,OP_DR),OP_1),tempa) &
          + prod( temp79b*h(:,OP_DZ),OP_1,OP_DR) &
          + prod(-temp79b*h(:,OP_DR),OP_1,OP_DZ)
  else
     temp79a = k(:,OP_1)*ri_79*j(:,OP_1)

     temp = prod( g(:,OP_DZ)*temp79a*h(:,OP_DZ),OP_DZ,OP_DR) &
          + prod(-g(:,OP_DZ)*temp79a*h(:,OP_DZ),OP_DR,OP_DZ) &
          + prod( g(:,OP_DR)*temp79a*h(:,OP_DR),OP_DZ,OP_DR) &
          + prod(-g(:,OP_DR)*temp79a*h(:,OP_DR),OP_DR,OP_DZ) &
          + prod( temp79a*g(:,OP_DZ )*h(:,OP_DZ),OP_DZ,OP_DR) &
          + prod( temp79a*g(:,OP_DR )*h(:,OP_DZ),OP_DR,OP_DR) &
          + prod(-temp79a*g(:,OP_DZ )*h(:,OP_DR),OP_DZ,OP_DZ) &
          + prod(-temp79a*g(:,OP_DR )*h(:,OP_DR),OP_DR,OP_DZ)
  end if
  tepsifkappar1 = (gam - 1.) * temp
#else
  tepsifkappar1%len = 0
#endif
end function tepsifkappar1

function tepsifkappar2(f,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsifkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempb

  if(gam.le.1.) then
     tepsifkappar2%len = 0
     return
  end if

  if(surface_int) then
     temp79a = k(:,OP_1)*ri_79* &
          (norm79(:,2)*f(:,OP_DR) - norm79(:,1)*f(:,OP_DZ))*j(:,OP_1)
     tempb = (mu(norm79(:,2),OP_DZ) + mu(norm79(:,1),OP_DR)) &
           * j(:,OP_1)*(-k(:,OP_1)*ri_79)


     temp = prod(temp79a*h(:,OP_DZ),OP_1,OP_DZ) &
          + prod(temp79a*h(:,OP_DR),OP_1,OP_DR) &
          + prod(mu( f(:,OP_DR )*h(:,OP_DZ),OP_1),tempb) &
          + prod(mu(-f(:,OP_DZ )*h(:,OP_DR),OP_1),tempb)
  else
     temp79a = k(:,OP_1)*ri_79*j(:,OP_1)

     temp = prod( f(:,OP_DR)*temp79a*h(:,OP_DZ),OP_DZ,OP_DZ) &
          + prod(-f(:,OP_DZ)*temp79a*h(:,OP_DZ),OP_DR,OP_DZ) &
          + prod( f(:,OP_DR)*temp79a*h(:,OP_DR),OP_DZ,OP_DR) &
          + prod(-f(:,OP_DZ)*temp79a*h(:,OP_DR),OP_DR,OP_DR) &
          + prod( temp79a*f(:,OP_DR )*h(:,OP_DZ),OP_DZ,OP_DZ) &
          + prod( temp79a*f(:,OP_DR )*h(:,OP_DZ),OP_DR,OP_DR) &
          + prod(-temp79a*f(:,OP_DZ )*h(:,OP_DR),OP_DZ,OP_DZ) &
          + prod(-temp79a*f(:,OP_DZ )*h(:,OP_DR),OP_DR,OP_DR)
  end if
  tepsifkappar2 = (gam - 1.) * temp
#else
  tepsifkappar2%len = 0
#endif
end function tepsifkappar2

function tepsifkappar3(f,g,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tepsifkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

  if(gam.le.1.) then
     tepsifkappar3%len = 0
     return
  end if

  if(surface_int) then
     temp79a = k(:,OP_1)*ri_79* &
          (norm79(:,2)*f(:,OP_DR) - norm79(:,1)*f(:,OP_DZ))*j(:,OP_1)
     temp79b = -k(:,OP_1)*ri_79* &
          (norm79(:,2)*g(:,OP_DZ) + norm79(:,1)*g(:,OP_DR))*j(:,OP_1)

     temp = prod( temp79a*g(:,OP_DZ),OP_1,OP_DZ) &
          + prod( temp79a*g(:,OP_DR),OP_1,OP_DR) &
          + prod( temp79b*f(:,OP_DR ),OP_1,OP_DZ) &
          + prod(-temp79b*f(:,OP_DZ ),OP_1,OP_DR)
  else
     temp79a = k(:,OP_1)*ri_79*j(:,OP_1)

     temp = prod( f(:,OP_DR)*temp79a*g(:,OP_DZ),OP_DZ,OP_DZ) &
          + prod(-f(:,OP_DZ)*temp79a*g(:,OP_DZ),OP_DR,OP_DZ) &
          + prod( f(:,OP_DR)*temp79a*g(:,OP_DR),OP_DZ,OP_DR) &
          + prod(-f(:,OP_DZ)*temp79a*g(:,OP_DR),OP_DR,OP_DR) &
          + prod( g(:,OP_DZ)*temp79a*f(:,OP_DR ),OP_DZ,OP_DZ) &
          + prod( g(:,OP_DR)*temp79a*f(:,OP_DR ),OP_DR,OP_DZ) &
          + prod(-g(:,OP_DZ)*temp79a*f(:,OP_DZ ),OP_DZ,OP_DR) &
          + prod(-g(:,OP_DR)*temp79a*f(:,OP_DZ ),OP_DR,OP_DR)
  end if
  tepsifkappar3 = (gam - 1.) * temp
#else
  tepsifkappar3%len = 0
#endif
end function tepsifkappar3

function tebfkappar(e,f,g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebfkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tebfkappar = 0.
     return
  end if

  if(surface_int) then
     temp79a = -ri2_79*k(:,OP_1)*j(:,OP_1)*f(:,OP_1)* &
          (norm79(:,1)*g(:,OP_DR) + norm79(:,2)*g(:,OP_DZ))

     temp = intx3(e(:,:,OP_1),temp79a,h(:,OP_DP))
  else
     temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)*j(:,OP_1)

     temp79b = g(:,OP_DZ)*(h(:,OP_DZ) )&
          +    g(:,OP_DR)*(h(:,OP_DR) )

     temp79d = temp79b*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 )

     temp = intx4(e(:,:,OP_DZ),g(:,OP_DZ),temp79a,h(:,OP_DP)) &
          + intx4(e(:,:,OP_DR),g(:,OP_DR),temp79a,h(:,OP_DP)) &
          + intx3(e(:,:,OP_DP),ri2_79,temp79d)
  end if
  tebfkappar = (gam - 1.) * temp
#else
  tebfkappar = 0.
#endif
end function tebfkappar

function tebfkappar1(g,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebfkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempd

  if(gam.le.1.) then
     tebfkappar1%len = 0
     return
  end if

  if(surface_int) then
     tempa = mu(-ri2_79*k(:,OP_1)*j(:,OP_1),OP_1)* &
          (norm79(:,1)*g(:,OP_DR) + norm79(:,2)*g(:,OP_DZ))

     temp = prod(mu(h(:,OP_DP),OP_1),tempa)
  else
     tempa = mu(k(:,OP_1)*ri2_79*j(:,OP_1),OP_1)

     temp79b = g(:,OP_DZ)*(h(:,OP_DZ) )&
          +    g(:,OP_DR)*(h(:,OP_DR) )

     tempd = mu(temp79b*j(:,OP_1 )*k(:,OP_1 ),OP_1)

     temp = prod(mu(g(:,OP_DZ)*h(:,OP_DP),OP_DZ),tempa) &
          + prod(mu(g(:,OP_DR)*h(:,OP_DP),OP_DR),tempa) &
          + prod(mu(ri2_79,OP_DP),tempd)
  end if
  tebfkappar1 = (gam - 1.) * temp
#else
  tebfkappar1%len = 0
#endif
end function tebfkappar1

function tebfkappar2(f,h,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebfkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa, tempb, tempd

  if(gam.le.1.) then
     tebfkappar2%len = 0
     return
  end if

  if(surface_int) then
     tempa = (mu(norm79(:,1),OP_DR) + mu(norm79(:,2),OP_DZ))* &
          (-ri2_79*k(:,OP_1)*j(:,OP_1)*f(:,OP_1))

     temp = prod(mu(h(:,OP_DP),OP_1),tempa)
  else
     temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)*j(:,OP_1)

     tempb = mu(h(:,OP_DZ),OP_DZ )&
        +    mu(h(:,OP_DR),OP_DR )

     tempd = tempb*(f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ))

     temp = prod(temp79a*h(:,OP_DP),OP_DZ,OP_DZ) &
          + prod(temp79a*h(:,OP_DP),OP_DR,OP_DR) &
          + prod(mu(ri2_79,OP_DP),tempd)
  end if
  tebfkappar2 = (gam - 1.) * temp
#else
  tebfkappar2%len = 0
#endif
end function tebfkappar2

function tebfkappar3(f,g,j,k)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebfkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempb, tempd

  if(gam.le.1.) then
     tebfkappar3%len = 0
     return
  end if

  if(surface_int) then
     temp79a = -ri2_79*k(:,OP_1)*j(:,OP_1)*f(:,OP_1)* &
          (norm79(:,1)*g(:,OP_DR) + norm79(:,2)*g(:,OP_DZ))

     temp = prod(temp79a,OP_1,OP_DP)
  else
     temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)*j(:,OP_1)

     tempb = mu(g(:,OP_DZ),OP_DZ)&
        +    mu(g(:,OP_DR),OP_DR)

     tempd = tempb*(f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ))

     temp = prod(g(:,OP_DZ)*temp79a,OP_DZ,OP_DP) &
          + prod(g(:,OP_DR)*temp79a,OP_DR,OP_DP) &
          + prod(mu(ri2_79,OP_DP),tempd)
  end if
  tebfkappar3 = (gam - 1.) * temp
#else
  tebfkappar3%len = 0
#endif
end function tebfkappar3

function tebfkapparl(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: tebfkapparl
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     tebfkapparl = 0.
     return
  end if

  temp79a = ri2_79*j(:,OP_1)*i(:,OP_1)* f(:,OP_1)


  if(surface_int) then
     temp = -intx5(e(:,:,OP_1),temp79a,h(:,OP_DP),norm79(:,1),g(:,OP_DR)) &
          -  intx5(e(:,:,OP_1),temp79a,h(:,OP_DP),norm79(:,2),g(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = ri2_79 * &
          (j(:,OP_DP)*i(:,OP_1)* f(:,OP_1) &
          +j(:,OP_1)*i(:,OP_DP)* f(:,OP_1) &
          +j(:,OP_1)*i(:,OP_1)* f(:,OP_DP))

     !  <T, f'> 
     temp79c =  (h(:,OP_DR)*g(:,OP_DR) + h(:,OP_DZ)*g(:,OP_DZ)) 
    

     ! d(temp79c)/dphi
     temp79d =  &
         +(h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
         +(h(:,OP_DRP)*g(:,OP_DR) + h(:,OP_DZP)*g(:,OP_DZ)) 

     temp = intx4(e(:,:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DP)) &
          + intx4(e(:,:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DP)) &
          - intx3(e(:,:,OP_1),temp79a,temp79d) &
          - intx3(e(:,:,OP_1),temp79b,temp79c)
  end if
  tebfkapparl = (gam - 1.) * temp
#else
  tebfkapparl = 0.
#endif
  return
end function tebfkapparl

function tebfkapparl2(f,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: tebfkapparl2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,i,j

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempc, tempd

  if(gam.le.1.) then
     tebfkapparl2%len = 0
     return
  end if

  temp79a = ri2_79*j(:,OP_1)*i(:,OP_1)* f(:,OP_1)


  if(surface_int) then
     temp = prod(-temp79a*h(:,OP_DP)*norm79(:,1),OP_1,OP_DR) &
          + prod(-temp79a*h(:,OP_DP)*norm79(:,2),OP_1,OP_DZ)
  else
     ! d(temp79a)/dphi
     temp79b = ri2_79 * &
          (j(:,OP_DP)*i(:,OP_1)* f(:,OP_1) &
          +j(:,OP_1)*i(:,OP_DP)* f(:,OP_1) &
          +j(:,OP_1)*i(:,OP_1)* f(:,OP_DP))

     !  <T, f'> 
     tempc = mu(h(:,OP_DR),OP_DR) + mu(h(:,OP_DZ),OP_DZ) 
    

     ! d(temp79c)/dphi
     tempd =  &
        mu(h(:,OP_DR),OP_DRP) + mu(h(:,OP_DZ),OP_DZP) &
       +mu(h(:,OP_DRP),OP_DR) + mu(h(:,OP_DZP),OP_DZ) 

     temp = prod(temp79a*h(:,OP_DP),OP_DR,OP_DR) &
          + prod(temp79a*h(:,OP_DP),OP_DZ,OP_DZ) &
          + prod(mu(-temp79a,OP_1),tempd) &
          + prod(mu(-temp79b,OP_1),tempc)
  end if
  tebfkapparl2 = (gam - 1.) * temp
#else
  tebfkapparl2%len = 0
#endif
  return
end function tebfkapparl2

function teffkappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: teffkappar
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  vectype, dimension(dofs_per_element) :: temp

  if(gam.le.1.) then
     teffkappar = 0.
     return
  end if

  if(surface_int) then
     temp79a =  k(:,OP_1)* &
          (norm79(:,2)*f(:,OP_DZ) + norm79(:,1)*f(:,OP_DR))*j(:,OP_1)

     temp = intx4(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx4(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DR))
  else
     temp79a = -k(:,OP_1)*j(:,OP_1)

     temp = intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DZ),h(:,OP_DZ)) &
          + intx5(e(:,:,OP_DZ),f(:,OP_DZ),temp79a,g(:,OP_DR),h(:,OP_DR)) &
          + intx5(e(:,:,OP_DR),f(:,OP_DR),temp79a,g(:,OP_DR),h(:,OP_DR))
  end if

  teffkappar = (gam - 1.) * temp
#else
  teffkappar = 0.
#endif
  return
end function teffkappar

function teffkappar1(g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: teffkappar1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp
  type(muarray) :: tempa

  if(gam.le.1.) then
     teffkappar1%len = 0
     return
  end if

  if(surface_int) then
     tempa = (mu(norm79(:,2),OP_DZ) + mu(norm79(:,1),OP_DR))*(k(:,OP_1)*j(:,OP_1))

     temp = prod(mu(g(:,OP_DZ)*h(:,OP_DZ),OP_1),tempa) &
          + prod(mu(g(:,OP_DR)*h(:,OP_DR),OP_1),tempa)
  else
     temp79a = -k(:,OP_1)*j(:,OP_1)

     temp = prod(temp79a*g(:,OP_DZ)*h(:,OP_DZ),OP_DZ,OP_DZ) &
          + prod(temp79a*g(:,OP_DZ)*h(:,OP_DZ),OP_DR,OP_DR) &
          + prod(temp79a*g(:,OP_DR)*h(:,OP_DR),OP_DZ,OP_DZ) &
          + prod(temp79a*g(:,OP_DR)*h(:,OP_DR),OP_DR,OP_DR)
  end if

  teffkappar1 = (gam - 1.) * temp
#else
  teffkappar1%len = 0
#endif
  return
end function teffkappar1

function teffkappar2(f,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: teffkappar2
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,h,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

  if(gam.le.1.) then
     teffkappar2%len = 0
     return
  end if

  if(surface_int) then
     temp79a =  k(:,OP_1)* &
          (norm79(:,2)*f(:,OP_DZ) + norm79(:,1)*f(:,OP_DR))*j(:,OP_1)

     temp = prod(temp79a*h(:,OP_DZ),OP_1,OP_DZ) &
          + prod(temp79a*h(:,OP_DR),OP_1,OP_DR)
  else
     temp79a = -k(:,OP_1)*j(:,OP_1)

     temp = prod(f(:,OP_DZ)*temp79a*h(:,OP_DZ),OP_DZ,OP_DZ) &
          + prod(f(:,OP_DR)*temp79a*h(:,OP_DZ),OP_DR,OP_DZ) &
          + prod(f(:,OP_DZ)*temp79a*h(:,OP_DR),OP_DZ,OP_DR) &
          + prod(f(:,OP_DR)*temp79a*h(:,OP_DR),OP_DR,OP_DR)
  end if

  teffkappar2 = (gam - 1.) * temp
#else
  teffkappar2%len = 0
#endif
  return
end function teffkappar2

function teffkappar3(f,g,j,k)

  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: teffkappar3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,j,k

#if defined(USE3D) || defined(USECOMPLEX)
  type(prodarray) :: temp

  if(gam.le.1.) then
     teffkappar3%len = 0
     return
  end if

  if(surface_int) then
     temp79a =  k(:,OP_1)* &
          (norm79(:,2)*f(:,OP_DZ) + norm79(:,1)*f(:,OP_DR))*j(:,OP_1)

     temp = prod(temp79a*g(:,OP_DZ),OP_1,OP_DZ) &
          + prod(temp79a*g(:,OP_DR),OP_1,OP_DR)
  else
     temp79a = -k(:,OP_1)*j(:,OP_1)

     temp = prod(f(:,OP_DZ)*temp79a*g(:,OP_DZ),OP_DZ,OP_DZ) &
          + prod(f(:,OP_DR)*temp79a*g(:,OP_DZ),OP_DR,OP_DZ) &
          + prod(f(:,OP_DZ)*temp79a*g(:,OP_DR),OP_DZ,OP_DR) &
          + prod(f(:,OP_DR)*temp79a*g(:,OP_DR),OP_DR,OP_DR)
  end if

  teffkappar3 = (gam - 1.) * temp
#else
  teffkappar3%len = 0
#endif
  return
end function teffkappar3

function q_delta(e)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: q_delta
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e

  if(surface_int) then
     q_delta = 0.
  else
     q_delta = intx4(e(:,:,OP_1),net79(:,OP_1),tit79(:,OP_1),qd79) &
          -    intx4(e(:,:,OP_1),net79(:,OP_1),tet79(:,OP_1),qd79)
  end if
end function q_delta

function q_delta1(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: q_delta1
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f

  if(surface_int) then
     q_delta1 = 0.
  else
     q_delta1 = intx4(e(:,:,OP_1),f(:,OP_1),net79(:,OP_1),qd79) 
  end if
end function q_delta1

vectype function q1ppsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


  if(surface_int) then
     temp = 0.
  else
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
  end if

  q1ppsi =   temp
  return
end function q1ppsi

vectype function q1pb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp = int5(ri2_79,e(:,OP_1),f(:,OP_DP),g(:,OP_1),h(:,OP_1))
  end if
#else
  temp = 0.
#endif

  q1pb =   temp
  return
end function q1pb

vectype function q1pf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

#if defined(USE3D) || defined(USECOMPLEX)
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = - int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
            - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  end if

  q1pf = temp
#else
  q1pf = 0.
#endif
  return
end function q1pf

! function t3tn(e,f,g)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: t3tn
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
!   vectype, dimension(dofs_per_element) :: temp

!   if(surface_int) then
!      temp = 0.
!   else
!      temp = intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_1))
!   end if

!   t3tn = temp
! end function t3tn

function t3tn(g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tn
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g
  type(prodarray) :: temp

  if(surface_int) then
     temp%len = 0
  else
     temp = prod(g(:,OP_1),OP_1,OP_1)
  end if

  t3tn = temp
end function t3tn

! function t3t(e,f)
!   use basic
!   use m3dc1_nint

!   implicit none

!   vectype, dimension(dofs_per_element) :: t3t
!   vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
!   vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
!   vectype, dimension(dofs_per_element) :: temp

!   if(surface_int) then
!      temp = 0.
!   else
!      temp = intx2(e(:,:,OP_1),f(:,OP_1))
!   end if

!   t3t = temp
! end function t3t

function t3t
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3t
  type(prodarray) :: temp

  if(surface_int) then
     temp%len = 0
  else
     temp79a = 1.
     temp = prod(temp79a,OP_1,OP_1)
  end if

  t3t = temp
end function t3t

function t3tnu(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tnu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,1),h(:,OP_DZ)) &
                - intx5(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,2),h(:,OP_DR))
        end if
     else
           temp = intx5(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_1),h(:,OP_DZ)) &
                - intx5(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_1),h(:,OP_DR))
        if(itor.eq.1) then
           temp = temp + &
                2.*(gam-1.)*intx4(e(:,:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_DZ))
        endif
     end if

  t3tnu = temp
end function t3tnu

function t3tnu1(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tnu1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           !temp = prod(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,1),h(:,OP_DZ)) &
                !- prod(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,2),h(:,OP_DR))
        end if
     else
           temp = prod( r_79*g(:,OP_1)*h(:,OP_DZ),OP_1,OP_DR) &
                + prod(-r_79*g(:,OP_1)*h(:,OP_DR),OP_1,OP_DZ)
        if(itor.eq.1) then
           temp = temp + &
                prod(2.*(gam-1)*g(:,OP_1)*h(:,OP_DZ),OP_1,OP_1)
        endif
     end if

  t3tnu1 = temp
end function t3tnu1

function t3tnu3(f,g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tnu3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           !temp = prod(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,1),h(:,OP_DZ)) &
                !- prod(e(:,:,OP_1),r_79,f(:,OP_1),norm79(:,2),h(:,OP_DR))
        end if
     else
           temp = prod( r_79*f(:,OP_DR)*g(:,OP_1),OP_1,OP_DZ) &
                + prod(-r_79*f(:,OP_DZ)*g(:,OP_1),OP_1,OP_DR)
        if(itor.eq.1) then
           temp = temp + &
                prod(2.*(gam-1)*f(:,OP_1)*g(:,OP_1),OP_1,OP_DZ)
        endif
     end if

  t3tnu3 = temp
end function t3tnu3

function t3tnv(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tnv
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp = 0.
     else
        temp = - intx4(e(:,:,OP_1),f(:,OP_DP),g(:,OP_1),h(:,OP_1)) &
             - (gam-1.)*intx4(e(:,:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_DP))
     endif
#else
  temp = 0.
#endif

  t3tnv = temp

  return
end function t3tnv

function t3tnv1(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tnv1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp%len = 0
     else
        temp = prod(-g(:,OP_1)*h(:,OP_1),OP_1,OP_DP) &
             + prod(-(gam-1.)*g(:,OP_1)*h(:,OP_DP),OP_1,OP_1)
     endif
#else
  temp%len = 0
#endif

  t3tnv1 = temp

  return
end function t3tnv1

function t3tnv3(f,g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tnv3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

#if defined(USE3D) || defined(USECOMPLEX)
     if(surface_int) then
        temp%len = 0
     else
        temp = prod(-f(:,OP_DP)*g(:,OP_1),OP_1,OP_1) &
             + prod(-(gam-1.)*f(:,OP_1)*g(:,OP_1),OP_1,OP_DP)
     endif
#else
  temp%len = 0
#endif

  t3tnv3 = temp

  return
end function t3tnv3

function t3tnchi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: t3tnchi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = &
                - intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,1),h(:,OP_DR)) &
                - intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,2),h(:,OP_DZ))
        endif
     else
        temp = -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_1),h(:,OP_DR))  &
               -intx5(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_1),h(:,OP_DZ)) &
               -(gam-1.)* &
               intx5(e(:,:,OP_1),ri2_79,f(:,OP_1),g(:,OP_1),h(:,OP_GS))
     endif

  t3tnchi = temp
end function t3tnchi

function t3tnchi1(g,h)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tnchi1
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: g,h
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           !temp = &
                !- prod(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,1),h(:,OP_DR)) &
                !- prod(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,2),h(:,OP_DZ))
        endif
     else
        temp =  prod(-ri2_79*g(:,OP_1)*h(:,OP_DR),OP_1,OP_DR)  &
               +prod(-ri2_79*g(:,OP_1)*h(:,OP_DZ),OP_1,OP_DZ) &
               +prod(-(gam-1.)*ri2_79*g(:,OP_1)*h(:,OP_GS),OP_1,OP_1)
     endif

  t3tnchi1 = temp
end function t3tnchi1

function t3tnchi3(f,g)
  use basic
  use m3dc1_nint

  implicit none

  type(prodarray) :: t3tnchi3
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g
  type(prodarray) :: temp

     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp%len = 0
        else
           !temp = &
                !- prod(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,1),h(:,OP_DR)) &
                !- prod(e(:,:,OP_1),ri2_79,f(:,OP_1),norm79(:,2),h(:,OP_DZ))
        endif
     else
        temp =  prod(-ri2_79*f(:,OP_DR)*g(:,OP_1),OP_1,OP_DR)  &
               +prod(-ri2_79*f(:,OP_DZ)*g(:,OP_1),OP_1,OP_DZ) &
               +prod(-(gam-1.)*ri2_79*f(:,OP_1)*g(:,OP_1),OP_1,OP_GS)
     endif

  t3tnchi3 = temp
end function t3tnchi3

vectype function j1b2ipsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

  temp = - int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_GS),h(:,OP_1))

  j1b2ipsib = temp
  return
end function j1b2ipsib
vectype function j1b2ibpsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

  temp = int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DR))    &
       + int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DZ))

  j1b2ibpsi = temp
  return
end function j1b2ibpsi

vectype function j1b2ipsif(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp
#if defined(USE3D) || defined(USECOMPLEX)
  temp = int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DRP))    &
       + int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DZP))    &
       - int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DRP),h(:,OP_DR))    &
       - int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZP),h(:,OP_DZ))
#else
  temp = 0
#endif

  j1b2ipsif = temp
  return
end function j1b2ipsif

vectype function j1b2ifb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = - int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR))    &
         + int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ))
#else
  temp = 0
#endif


  j1b2ifb = temp
  return
end function j1b2ifb

vectype function j1b2iff(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = - int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DRP))    &
         + int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZP))
#else
  temp = 0
#endif

  j1b2iff = temp
  return
end function j1b2iff

vectype function j1b2ipsipsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  temp = + int5(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZP),h(:,OP_DR))    &
         - int5(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DRP),h(:,OP_DZ))
#else
  temp = 0
#endif

  j1b2ipsipsi = temp
  return
end function j1b2ipsipsi


vectype function pparpu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp = int4(r_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR)) &
             - int4(r_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ))

     end if

  pparpu = temp

  return
end function pparpu

vectype function pparpupsipsib2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
        temp = - 2.*int5(ri_79,temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DR))  &
               + 2.*int5(ri_79,temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DZ))
        temp = temp                                                    &
               + int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DR)) &
               + int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DR)) &
               + int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
               + int5(ri_79,temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DR)) &
               - int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DZ)) &
               - int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DZ)) &
               - int5(ri_79,temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DZ)) &
               - int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
        temp = temp                                                    &
               -2*int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DRZ)) &
               -2*int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
               -2*int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DR),g(:,OP_DZZ)) &
               -2*int5(ri_79,temp79a,h(:,OP_DZZ),i(:,OP_DR),g(:,OP_DZ)) &
               +2*int5(ri_79,temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DRR)) &
               +2*int5(ri_79,temp79a,h(:,OP_DRR),i(:,OP_DZ),g(:,OP_DR)) &
               +2*int5(ri_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DRZ)) &
               +2*int5(ri_79,temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
        if(itor.eq.1) then
           temp = temp                                                    &
                  + 2.*int5(ri2_79,temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DZ)) &
                  + 2.*int5(ri2_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DZ)) 
        endif
     end if

  pparpupsipsib2 = temp

  return
end function pparpupsipsib2

vectype function pparpupsibb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp = 0
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
        temp =  -2.*int5(ri2_79,temp79a,i(:,OP_1),g(:,OP_DRP),h(:,OP_DR))  &
                -2.*int5(ri2_79,temp79a,i(:,OP_1),g(:,OP_DZP),h(:,OP_DZ))
#else
        temp = 0
#endif
     end if

     pparpupsibb2 = temp
     return
end function pparpupsibb2

vectype function pparpubbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp = 0
        if(itor.eq.1) then
        temp79a = 2.*e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
           temp = temp                                                    &
                  + int5(ri2_79,temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DZ)) 
        endif
     end if

     pparpubbb2 = temp
     return
end function pparpubbb2

! ===========
vectype function pparpvpsibb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
        temp = - 2.*int5(ri_79,temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1))  &
               + 2.*int5(ri_79,temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1))
     end if

  pparpvpsibb2 = temp

  return
end function pparpvpsibb2

vectype function pparpvbbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp
     temp = 0.
     if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp = -int3(e(:,OP_1),f(:,OP_1),g(:,OP_DP)) 
        temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
        temp = temp                                                    &
             - 2.*int5(ri2_79,temp79a,g(:,OP_DP),h(:,OP_1),i(:,OP_1))  
#endif
     end if

  pparpvbbb2 = temp

  return
end function pparpvbbb2

vectype function pparpchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR)) &
             + int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ))

     end if

  pparpchi = temp

  return
end function pparpchi
vectype function pparpchipsipsib2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
        temp = - 2.*int5(ri4_79,temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DR))  &
               - 2.*int5(ri4_79,temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DZ))
        temp = temp                                                     &
               + int5(ri4_79,temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DR)) &
               + int5(ri4_79,temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DR)) &
               + int5(ri4_79,temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DR)) &
               + int5(ri4_79,temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DR)) &
               + int5(ri4_79,temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DZ)) &
               + int5(ri4_79,temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DZ)) &
               + int5(ri4_79,temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DZ)) &
               + int5(ri4_79,temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DZ))
        temp = temp                                                       &
               -2.*int5(ri4_79,temp79a,g(:,OP_DZZ),h(:,OP_DR),i(:,OP_DR)) &
               -2.*int5(ri4_79,temp79a,g(:,OP_DZ),h(:,OP_DRZ),i(:,OP_DR)) &
               +2.*int5(ri4_79,temp79a,g(:,OP_DRZ),h(:,OP_DZ),i(:,OP_DR)) &
               +2.*int5(ri4_79,temp79a,g(:,OP_DR),h(:,OP_DZZ),i(:,OP_DR)) &
               +2.*int5(ri4_79,temp79a,g(:,OP_DRZ),h(:,OP_DR),i(:,OP_DZ)) &
               +2.*int5(ri4_79,temp79a,g(:,OP_DZ),h(:,OP_DRR),i(:,OP_DZ)) &
               -2.*int5(ri4_79,temp79a,g(:,OP_DRR),h(:,OP_DZ),i(:,OP_DZ)) &
               -2.*int5(ri4_79,temp79a,g(:,OP_DR),h(:,OP_DRZ),i(:,OP_DZ))
        if(itor.eq.1) then
           temp = temp                                                       &
                  - 2.*int5(ri5_79,temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DR)) &
                  - 2.*int5(ri5_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR))
           temp = temp                                                       &
                  -6.*int5(ri5_79,temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DZ))  &
                  +6.*int5(ri5_79,temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR)) 
        endif
     end if

  pparpchipsipsib2 = temp

  return
end function pparpchipsipsib2

vectype function pparpchipsibb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp
     temp = 0.
     if(surface_int) then
        temp = 0.
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
        temp = temp                                                      &
             - 2.*int5(ri2_79,temp79a,g(:,OP_DZP),h(:,OP_DR),i(:,OP_1))  &
             + 2.*int5(ri2_79,temp79a,g(:,OP_DRP),h(:,OP_DZ),i(:,OP_1))  
#endif
     end if

  pparpchipsibb2 = temp

  return
end function pparpchipsibb2

vectype function pparpchibbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  vectype :: temp

     if(surface_int) then
        temp = 0.
     else
        temp = 0
        if(itor.eq.1) then
        temp79a = 2.*e(:,OP_1)*f(:,OP_1)*j(:,OP_1)
           temp = temp                                                    &
                  - int5(ri5_79,temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DR)) 
        endif
     end if

     pparpchibbb2 = temp
     return
end function pparpchibbb2


function pperpu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpu
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = 2.*intx4(e(:,:,OP_DZ),r_79,f(:,OP_1),g(:,OP_DR)) &
          - 2.*intx4(e(:,:,OP_DR),r_79,f(:,OP_1),g(:,OP_DZ)) &
          + intx4(e(:,:,OP_1),r_79,f(:,OP_DZ),g(:,OP_DR)) &
          - intx4(e(:,:,OP_1),r_79,f(:,OP_DR),g(:,OP_DZ))
  end if

  pperpu = temp
end function pperpu

function pperpupsipsib2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpupsipsib2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp79a = ri_79*f(:,OP_1)*j(:,OP_1)
     temp =   intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DR))  &
          -   intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DZ))
     temp = temp                                                    &
          - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DR)) &
          - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DR)) &
          - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
          - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DR)) &
          + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DZ)) &
          + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DZ)) &
          + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DZ)) &
          + .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
     temp = temp                                                    &
          + intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DRZ)) &
          + intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DR)) &
          + intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DR),g(:,OP_DZZ)) &
          + intx5(e(:,:,OP_1),temp79a,h(:,OP_DZZ),i(:,OP_DR),g(:,OP_DZ)) &
          - intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DRR)) &
          - intx5(e(:,:,OP_1),temp79a,h(:,OP_DRR),i(:,OP_DZ),g(:,OP_DR)) &
          - intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DRZ)) &
          - intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DZ))
     if(itor.eq.1) then
        temp79a = ri2_79*f(:,OP_1)*j(:,OP_1)

        temp = temp                                                    &
             - intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DZ)) &
             - intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DZ)) 
     endif
  end if
  
  pperpupsipsib2 = temp
end function pperpupsipsib2

function pperpupsibb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpupsibb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  if(surface_int) then
     temp = 0.
  else
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = ri2_79*f(:,OP_1)*j(:,OP_1)
     temp = temp                                                      &
          + intx5(e(:,:,OP_1),temp79a,g(:,OP_DRP),h(:,OP_DR),i(:,OP_1))  &
          + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZP),h(:,OP_DZ),i(:,OP_1))  
#endif
  end if

  pperpupsibb2 = temp
end function pperpupsibb2


function pperpubbb2(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpubbb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = 0
     if(itor.eq.1) then
        temp79a =   ri2_79*f(:,OP_1)*j(:,OP_1)
        temp = temp                                                    &
             - intx5(e(:,:,OP_1),temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DZ)) 
     endif
  end if

  pperpubbb2 = temp
end function pperpubbb2

! ===========
function pperpvpsibb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpvpsibb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp79a =   ri_79*f(:,OP_1)*j(:,OP_1)
     temp =   intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1))  &
          - intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1))
  end if

  pperpvpsibb2 = temp
end function pperpvpsibb2

function pperpvbbb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpvbbb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  temp = 0.
  if(surface_int) then
     temp = 0.
  else
#if defined(USE3D) || defined(USECOMPLEX)
     temp = -2.*intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_DP)) 
     temp79a =   ri2_79*f(:,OP_1)*j(:,OP_1)
     temp = temp                                       &
          + intx5(e(:,:,OP_1),temp79a, g(:,OP_DP),h(:,OP_1),i(:,OP_1))
#endif
  end if

  pperpvbbb2 = temp
end function pperpvbbb2

function pperpchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchi

  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else   
     temp = 2.*intx4(e(:,:,OP_DR),ri2_79,f(:,OP_1),g(:,OP_DR)) &
          + 2.*intx4(e(:,:,OP_DZ),ri2_79,f(:,OP_1),g(:,OP_DZ)) &
          + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DR),g(:,OP_DR)) &
          + intx4(e(:,:,OP_1),ri2_79,f(:,OP_DZ),g(:,OP_DZ))
  end if

  pperpchi = temp
end function pperpchi


function pperpchipsipsib2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchipsipsib2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

     if(surface_int) then
        temp = 0.
     else
        temp79a =   ri4_79*f(:,OP_1)*j(:,OP_1)
        temp =   intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DR),g(:,OP_DR))  &
               + intx5(e(:,:,OP_1),temp79a,h(:,OP_GS),i(:,OP_DZ),g(:,OP_DZ))
        temp = temp                                                    &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRR),g(:,OP_DR)) &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DRZ),g(:,OP_DR)) &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRR),i(:,OP_DR),g(:,OP_DR)) &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DZ),g(:,OP_DR)) &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DRZ),g(:,OP_DZ)) &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZZ),g(:,OP_DZ)) &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DRZ),i(:,OP_DR),g(:,OP_DZ)) &
               - .5*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZZ),i(:,OP_DZ),g(:,OP_DZ))

         temp = temp                                                       &
               + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZZ),h(:,OP_DR),i(:,OP_DR)) &
               + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DRZ),i(:,OP_DR)) &
               - intx5(e(:,:,OP_1),temp79a,g(:,OP_DRZ),h(:,OP_DZ),i(:,OP_DR)) &
               - intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DZZ),i(:,OP_DR)) &
               - intx5(e(:,:,OP_1),temp79a,g(:,OP_DRZ),h(:,OP_DR),i(:,OP_DZ)) &
               - intx5(e(:,:,OP_1),temp79a,g(:,OP_DZ),h(:,OP_DRR),i(:,OP_DZ)) &
               + intx5(e(:,:,OP_1),temp79a,g(:,OP_DRR),h(:,OP_DZ),i(:,OP_DZ)) &
               + intx5(e(:,:,OP_1),temp79a,g(:,OP_DR),h(:,OP_DRZ),i(:,OP_DZ))
        if(itor.eq.1) then
           temp79a =   ri5_79*f(:,OP_1)*j(:,OP_1)
           temp = temp                                                    &
                  + intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DR),g(:,OP_DR)) &
                  + intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR))
           temp = temp                                                       &
                  +3.*intx5(e(:,:,OP_1),temp79a,h(:,OP_DR),i(:,OP_DZ),g(:,OP_DZ))  &
                  -3.*intx5(e(:,:,OP_1),temp79a,h(:,OP_DZ),i(:,OP_DZ),g(:,OP_DR))  
        endif
     end if

  pperpchipsipsib2 = temp
end function pperpchipsipsib2


function pperpchipsibb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchipsibb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp
  temp = 0.
  if(surface_int) then
     temp = 0.
  else
#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = ri2_79*f(:,OP_1)*j(:,OP_1)
     temp = temp                                                      &
          + intx5(e(:,:,OP_1),temp79a,g(:,OP_DZP),h(:,OP_DR),i(:,OP_1))  &
          - intx5(e(:,:,OP_1),temp79a,g(:,OP_DRP),h(:,OP_DZ),i(:,OP_1))  
#endif
  end if

  pperpchipsibb2 = temp
end function pperpchipsibb2

function pperpchibbb2(e,f,g,h,i,j)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: pperpchibbb2
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i,j

  vectype, dimension(dofs_per_element) :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = 0
     if(itor.eq.1) then
        temp79a =  ri5_79*f(:,OP_1)*j(:,OP_1)
        temp = temp                                                    &
             + intx5(e(:,:,OP_1),temp79a,h(:,OP_1),i(:,OP_1),g(:,OP_DR)) 
     endif
  end if

  pperpchibbb2 = temp
end function pperpchibbb2


function incvb(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: incvb
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(surface_int) then
     incvb = 0.
  else
     incvb = intx3(e(:,:,OP_1),f(:,OP_1),g(:,OP_1))
  end if
end function incvb

function incupsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: incupsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(surface_int) then
     incupsi = 0.
  else
     incupsi = intx3(e(:,:,OP_1),f(:,OP_DR),g(:,OP_DR)) &
          +    intx3(e(:,:,OP_1),f(:,OP_DZ),g(:,OP_DZ))
  end if
end function incupsi

function incchipsi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: incchipsi
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g

  if(surface_int) then
     incchipsi = 0.
  else
     incchipsi = intx4(e(:,:,OP_1),ri3_79,f(:,OP_DZ),g(:,OP_DR)) &
          -      intx4(e(:,:,OP_1),ri3_79,f(:,OP_DR),g(:,OP_DZ))
     
  end if
end function incchipsi

subroutine JxB_r(o, opol)
  use m3dc1_nint

  implicit none

  vectype, intent(out), dimension(MAX_PTS) :: o, opol
  vectype, dimension(MAX_PTS) :: otor

  otor = - ri2_79*pst79(:,OP_GS)*pstx79(:,OP_DR)
  opol = - ri2_79*bztx79(:,OP_1)*bzt79(:,OP_DR)
#if defined(USE3D) || defined(USECOMPLEX)
  otor = otor + ri_79*pst79(:,OP_GS)*bfptx79(:,OP_DZ)
  opol = opol &
       - ri2_79*bztx79(:,OP_1)*bfpt79(:,OP_DRP) &
       - ri3_79*bztx79(:,OP_1)*pst79(:,OP_DZP)
#endif
  o = opol + otor
end subroutine JxB_r

subroutine JxB_phi(o)
  use m3dc1_nint

  implicit none
  
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = &
       + ri2_79*bzt79(:,OP_DZ)*pstx79(:,OP_DR) &
       - ri2_79*bzt79(:,OP_DR)*pstx79(:,OP_DZ)

#if defined(USE3D) || defined(USECOMPLEX)
  o = o &
       + ri2_79*bfpt79(:,OP_DZP)*pstx79(:,OP_DR) &
       - ri2_79*bfpt79(:,OP_DRP)*pstx79(:,OP_DZ) &
       - ri_79*bzt79(:,OP_DZ)*bfptx79(:,OP_DZ) &
       - ri_79*bzt79(:,OP_DR)*bfptx79(:,OP_DR) &
       - ri_79*bfpt79(:,OP_DZP)*bfptx79(:,OP_DZ) &
       - ri_79*bfpt79(:,OP_DRP)*bfptx79(:,OP_DR) &
       - ri3_79*pst79(:,OP_DZP)*pstx79(:,OP_DZ) &
       - ri3_79*pst79(:,OP_DRP)*pstx79(:,OP_DR) &
       - ri2_79*pst79(:,OP_DZP)*bfptx79(:,OP_DR) &
       + ri2_79*pst79(:,OP_DRP)*bfptx79(:,OP_DZ)
#endif       

end subroutine JxB_phi
       
subroutine JxB_z(o, opol)
  use m3dc1_nint

  implicit none

  vectype, intent(out), dimension(MAX_PTS) :: o, opol
  vectype, dimension(MAX_PTS) :: otor

  otor = - ri2_79*pst79(:,OP_GS)*pstx79(:,OP_DZ)
  opol = - ri2_79*bztx79(:,OP_1)*bzt79(:,OP_DZ)

#if defined(USE3D) || defined(USECOMPLEX)
  otor = otor - ri_79*pst79(:,OP_GS)*bfptx79(:,OP_DR)
  opol = opol  &
       - ri2_79*bztx79(:,OP_1)*bfpt79(:,OP_DZP) &
       + ri3_79*bztx79(:,OP_1)*pst79(:,OP_DRP)
#endif
  o = opol + otor
  
end subroutine JxB_z

end module metricterms_new

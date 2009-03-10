!======================================================================
! boundary_node
!
! determines if node is on boundary, and returns relevant info
! about boundary surface
!======================================================================
subroutine boundary_node(inode,is_boundary,izone,izonedim,normal,curv,x,z)
  use basic

  implicit none

  integer, intent(in) :: inode              ! node index
  integer, intent(out) :: izone,izonedim    ! zone type/dimension
  real, intent(out) :: normal(2), curv
  real, intent(out) :: x,z                  ! coordinates of inode
  logical, intent(out) :: is_boundary       ! is inode on boundary

  integer :: ibottom, iright, ileft, itop
  real :: angler

  curv = 0.

  select case(nonrect)
  case(0)
     call zonenod(inode,izone,izonedim)

     if(izonedim.ge.2) then
        is_boundary = .false.
        return
     end if

     call getmodeltags(ibottom, iright, itop, ileft)

     ! for periodic bc's
     ! skip if on a periodic boundary
     ! and convert corner to an edge when one edge is periodic
     if(iper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) then
           is_boundary = .false.
           return
        end if
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) then
           is_boundary = .false.
           return
        end if
     endif
     
     is_boundary = .true.

     ! assign normal vector
     !....convert tiltangle to radians and define normal vectors
     angler = tiltangled*pi/180.
     if(izone.eq.iright) then
        normal(1) = cos(angler)  !cos
        normal(2) = sin(angler)  !sin
     else if(izone.eq.ileft) then
        normal(1) = cos(pi+angler)  !cos
        normal(2) = sin(pi+angler)  !sin
     else if(izone.eq.itop) then
        normal(1) = cos(pi/2. + angler)  !cos
        normal(2) = sin(pi/2. + angler)  !sin
     else if(izone.eq.ibottom) then
        normal(1) = cos(3.*pi/2. + angler)  !cos
        normal(2) = sin(3.*pi/2. + angler)  !sin
     else
        print *, "Error: unknown zone tag ", izone
        normal = 0.
     endif

     if(izonedim.eq.0) then
        normal(1) = 0.
        normal(2) = 1.
     end if

  case(1)
     call nodNormalVec(inode, normal, is_boundary)

     if(.not.is_boundary) return
     izonedim=1      !cj set izonedim always 1 for the shaped boundary
     izone=1         !cj dummy for shaped boundary

     call nodcurvature2(inode, curv, is_boundary)
     if(.not.is_boundary) return
  end select
  
  call nodcoord(inode,x,z)

end subroutine boundary_node

!======================================================================
! rotate_vector
! ~~~~~~~~~~~~~
!
! Performs coordinate rotation from (R, Z) to (n, t) on invec,
! returns result in outvec.
!======================================================================
subroutine rotate_vector(invec, outvec, normal, curv)
  implicit none
  vectype, intent(in), dimension(6) :: invec
  vectype, intent(out), dimension(6) :: outvec
  real, intent(in) :: curv, normal(2)

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
end subroutine rotate_vector


!======================================================================
! rotate_matrix
! ~~~~~~~~~~~~~
!
! Performs coordinate rotation from (R, Z) to (n, t) on imatrix
!======================================================================
subroutine rotate_matrix(imatrix, ibegin, normal, curv, rhs)
  use basic
  implicit none
  integer, intent(in) :: imatrix, ibegin
  real, intent(in) :: curv, normal(2)
  vectype, dimension(*), intent(inout) :: rhs
  integer :: row(5), i
  real :: fac(2)
  vectype, dimension(6) :: temp
  
  if(imatrix.ne.0) then
     do i=1,5
        row(i) = ibegin+i
     end do
     call applyLinCombinationForMatrix2(imatrix,&
          row(1),row(2),row(3),row(4),row(5),normal,curv,icomplex)
  endif

  call rotate_vector(rhs(ibegin:ibegin+5),temp,normal,curv)
  rhs(ibegin:ibegin+5) = temp

end subroutine




!======================================================================
! set_dirichlet_bc
!======================================================================
subroutine set_dirichlet_bc(imatrix,ibegin,rhs,bv,normal,curv,izonedim)
  use basic
  implicit none

  integer, intent(in) :: imatrix              ! matrix handle
  integer, intent(in) :: ibegin               ! first dof of field
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(6) :: bv     ! boundary values
  real, intent(in) :: normal(2), curv
  integer, intent(in) :: izonedim             ! dimension of boundary
  
  !clamp value
  if(imatrix.ne.0) call setdiribc(imatrix, ibegin)
  rhs(ibegin) = bv(1)
     
  ! clamp tangential derivative
  call set_tangent_bc(imatrix,ibegin,rhs,bv,normal,curv,izonedim)

end subroutine set_dirichlet_bc


!======================================================================
! set_tangent_bc
!======================================================================
subroutine set_tangent_bc(imatrix,ibegin,rhs,bv,normal,curv,izonedim)
  use basic
  use arrays
  use gradshafranov

  implicit none
  
  integer, intent(in) :: imatrix              ! matrix handle
  integer, intent(in) :: ibegin               ! first dof of field
  real, intent(in) :: normal(2), curv
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, dimension(6) :: temp
  vectype, intent(in), dimension(6) :: bv     ! boundary values
  integer, intent(in) :: izonedim             ! dimension of boundary

  integer :: irow, numvals, i
  integer, dimension(5) :: cols
  vectype, dimension(5) :: vals
  vectype, dimension(6) :: bv_rotated
  real :: newnormal(2)

  if(imatrix.ne.0) then
     do i=1,5 
        cols(i) = ibegin + i
     end do
  endif

  call rotate_vector(bv, bv_rotated, normal, curv)
     
  ! t
  irow = ibegin+2
  if(imatrix.ne.0) then
     numvals = 2
     vals(1) = -normal(2)
     vals(2) =  normal(1)
     call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
  endif
  rhs(irow) = bv_rotated(3)
  
  ! tt
  irow = ibegin+5
  if(imatrix.ne.0) then
     numvals = 5
     vals(1) = -curv*normal(1)
     vals(2) = -curv*normal(2)
     vals(3) = normal(2)**2
     vals(4) = -2.*normal(1)*normal(2)
     vals(5) = normal(1)**2
     call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
  endif
  rhs(irow) = bv_rotated(6)
  
  if(izonedim.eq.0) then
     ! n
     irow = ibegin+1
     if(imatrix.ne.0) then     
        numvals = 2
        vals(1) = normal(1)
        vals(2) = normal(2)
        call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     endif
     rhs(irow) = bv_rotated(2)

     ! nn
     irow = ibegin+3
     if(imatrix.ne.0) then
        numvals = 3
        vals(3) = normal(1)**2
        vals(4) = 2.*normal(1)*normal(2)
        vals(5) = normal(2)**2
        call setgeneralbc(imatrix,irow,numvals,cols(3:5),vals(3:5),icomplex)
     endif
     rhs(irow) = bv_rotated(4)
  endif

end subroutine set_tangent_bc


!======================================================================
! set_normal_bc
!======================================================================
subroutine set_normal_bc(imatrix,ibegin,rhs,bv,normal,curv,izonedim)
  use basic

  implicit none
  
  integer, intent(in) :: imatrix              ! matrix handle
  integer, intent(in) :: ibegin               ! first dof of field
  real, intent(in) :: normal(2), curv
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(6) :: bv     ! boundary values
  integer, intent(in) :: izonedim             ! dimension of boundary

  integer :: irow, numvals, i
  integer, dimension(5) :: cols
  vectype, dimension(5) :: vals
  vectype, dimension(6) :: bv_rotated
  real :: newnormal(2)

  if(imatrix.ne.0) then
     do i=1,5 
        cols(i) = ibegin + i
     end do
  endif

  call rotate_vector(bv, bv_rotated, normal, curv)
     
  ! n
  irow = ibegin+1
  if(imatrix.ne.0) then
     numvals = 2
     vals(1) = normal(1)
     vals(2) = normal(2)
     call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
  endif
  rhs(irow) = bv_rotated(2)

  ! nt
  irow = ibegin+4
  if(imatrix.ne.0) then
     numvals = 5
     vals(1) = -curv*normal(2)
     vals(2) = curv*normal(1)
     vals(3) = -normal(1)*normal(2)
     vals(4) = normal(1)**2 - normal(2)**2
     vals(5) = normal(1)*normal(2)
     call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
  endif
  rhs(irow) = bv_rotated(5)

  if(izonedim.eq.0) then
     ! t
     irow = ibegin+2
     if(imatrix.ne.0) then
        numvals = 2
        vals(1) = -normal(2)
        vals(2) =  normal(1)
        call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     endif
     rhs(irow) = bv_rotated(3)
  endif
     
end subroutine set_normal_bc
   
   
!======================================================================
! set_laplacian_bc
!======================================================================
subroutine set_laplacian_bc(imatrix,ibegin,rhs,bv,normal,curv,izonedim,radius)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field
  real, intent(in) :: normal(2), curv
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(6) :: bv     ! boundary values
  integer, intent(in) :: izonedim    ! dimension of boundary
  real, intent(in) :: radius         ! radial coordinate of node
                                     ! for gs operator, let radius -> -radius
  
  integer :: numvals, irow
  integer, dimension(3) :: cols
  vectype, dimension(3) :: vals

  if(itor.eq.1) then
     numvals = 3
     cols(1) = ibegin + 1
     cols(2) = ibegin + 3
     cols(3) = ibegin + 5
     vals(1) =  1./radius
     vals(2) =  1.
     vals(3) =  1.
  else
     numvals = 2
     cols(1) = ibegin + 3
     cols(2) = ibegin + 5
     vals(1) = 1.
     vals(2) = 1.
  endif

  if(izonedim.eq.1) then
     ! edge points
     ! ~~~~~~~~~~~
     irow = ibegin + 3
     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = 0.
  end if
  
end subroutine set_laplacian_bc


!=======================================================
! boundary_vel
! ~~~~~~~~~~~~
!
! sets boundary conditions for velocity fields
!=======================================================
subroutine boundary_vel(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  real :: normal(2), curv
  real :: x, z
  vectype, dimension(6) :: temp

  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return 

  call numnod(numnodes)
  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call entdofs(vecsize_vel, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin+u_off, normal, curv, rhs)
     if(numvar.ge.2) &
          call rotate_matrix(imatrix, ibegin+vz_off, normal, curv, rhs)
     if(numvar.ge.3) &
          call rotate_matrix(imatrix, ibegin+chi_off, normal, curv, rhs)

     call assign_local_pointers(i)

     ! no normal flow
     if(inonormalflow.eq.1) then
        temp = 0.
        call set_dirichlet_bc(imatrix,ibegin+u_off,rhs,temp, &
             normal,curv,izonedim)
        if(numvar.ge.3) then
           call set_normal_bc(imatrix,ibegin+chi_off,rhs,temp, &
                normal,curv,izonedim)
        endif
     end if
     
     ! no poloidal slip
     if(inoslip_pol.eq.1) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+u_off,rhs,temp, &
             normal,curv,izonedim)
        if(numvar.ge.3) then
           call set_dirichlet_bc(imatrix,ibegin+chi_off,rhs,temp, &
                normal,curv,izonedim)
        endif
     end if

     ! toroidal velocity
     if(numvar.ge.2) then
        ! no slip
        if(inoslip_tor.eq. 1) then
           temp = vzs_l
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*vzo_v(ibegin+vz_off:ibegin+vz_off+5)
           endif
           call set_dirichlet_bc(imatrix,ibegin+vz_off,rhs,temp, &
                normal,curv,izonedim)
        end if
        
        ! no toroidal stress
        if(inostress_tor.eq.1) then
           temp = 0.
           call set_normal_bc(imatrix,ibegin+vz_off,rhs,temp, &
                normal,curv,izonedim)
        end if
     endif
       
     ! no vorticity
     temp = 0.
     call set_laplacian_bc(imatrix,ibegin+u_off,rhs,temp,normal, &
          curv,izonedim,-x)

     ! no compression
     if(com_bc.eq.1 .and. numvar.ge.3) then
        call set_laplacian_bc(imatrix,ibegin+chi_off,rhs,temp, &
             normal,curv,izonedim,x)
     endif
  end do

end subroutine boundary_vel


!=======================================================
! boundary_mag
! ~~~~~~~~~~~~
!
! sets boundary conditions for magnetic fields
! and electron pressure 
!=======================================================
subroutine boundary_mag(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  integer :: ibegin1, iendplusone1
  real :: normal(2), curv
  real :: x, z
  vectype, dimension(6) :: temp

  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_mag called"

  call numnod(numnodes)

  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call entdofs(1, i, 0, ibegin1, iendplusone1)
     call entdofs(vecsize_phi, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin+psi_off, normal, curv, rhs)
     if(numvar.ge.2) &
          call rotate_matrix(imatrix, ibegin+bz_off, normal, curv, rhs)
     if(numvar.ge.3) &
          call rotate_matrix(imatrix, ibegin+pe_off, normal, curv, rhs)

     call assign_local_pointers(i)

     ! clamp poloidal field
     temp = psis_l
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp = 1.5*temp + 0.5*psio_v(ibegin+psi_off:ibegin+psi_off+5)
     endif
     call set_dirichlet_bc(imatrix,ibegin+psi_off,rhs,temp,normal,curv,izonedim)

     ! add loop voltage
     if(jadv.eq.0 .and. igauge.eq.0) then
        if(integrator.eq.1 .and. ntime.gt.1) then
           rhs(ibegin+psi_off) = rhs(ibegin+psi_off) + 1.5*fbound
        else
           rhs(ibegin+psi_off) = rhs(ibegin+psi_off) + fbound
        endif
     endif


     ! clamp toroidal field
     if(iconst_bz.eq.1 .and. numvar.ge.2) then
        temp = bzs_l
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*bzo_v(ibegin+bz_off:ibegin+bz_off+5)
        endif
        call set_dirichlet_bc(imatrix,ibegin+bz_off,rhs,temp, &
             normal,curv,izonedim)
     endif

     ! no toroidal current
     if(inocurrent_tor.eq.1) then
        temp = 0.
        if(jadv.eq.1 .and. igauge.eq.0) then
           temp(1) = vloop/(2.*pi*resistivity(ibegin1))
        endif
        call set_laplacian_bc(imatrix,ibegin+psi_off,rhs,temp, &
             normal,curv,izonedim,-x)
     end if

     ! no tangential current
     if(inocurrent_pol.eq.1 .and. numvar.ge.2) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+bz_off,rhs,temp,normal,curv,izonedim)
     end if

     if(numvar.ge.3) then 
        if(inograd_p.eq.1) then
           temp = 0.
           call set_normal_bc(imatrix,ibegin+pe_off,rhs,temp, &
                normal,curv,izonedim)
        end if
        if(iconst_p.eq.1) then
           temp = pes_l
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*peo_v(ibegin+pe_off:ibegin+pe_off+5)
           endif
           call set_dirichlet_bc(imatrix,ibegin+pe_off,rhs,temp, &
                normal,curv,izonedim)
        else if(iconst_t.eq.1) then
           temp = pes_l*den1_l(1)/dens_l(1)
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*peo_v(ibegin+pe_off:ibegin+pe_off+5)
           endif
           call set_dirichlet_bc(imatrix,ibegin+pe_off,rhs,temp, &
                normal,curv,izonedim)
        end if

     endif
  end do

end subroutine boundary_mag


!=======================================================
! boundary_den
! ~~~~~~~~~~~~
!
! sets boundary conditions for density
!=======================================================
subroutine boundary_den(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  real :: normal(2), curv
  real :: x,z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_den called"

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call entdofs(vecsize_n, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin+den_off, normal, curv, rhs)
     call assign_local_pointers(i)

     if(inograd_n.eq.1) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+den_off,rhs,temp,normal,curv,izonedim)
     end if
     if(iconst_n.eq.1) then
        temp = dens_l
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*deno_v(ibegin+den_off:ibegin+den_off+5)
        endif
        call set_dirichlet_bc(imatrix,ibegin+den_off,rhs,temp,normal,curv,izonedim)
     end if

  end do

end subroutine boundary_den


!=======================================================
! boundary_pres
! ~~~~~~~~~~~~~
!
! sets boundary conditions for total pressure
!=======================================================
subroutine boundary_pres(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_pres called"

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call entdofs(vecsize_p, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin+p_off, normal, curv, rhs)
     call assign_local_pointers(i)

     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+p_off,rhs,temp,normal,curv,izonedim)
     end if
     if(iconst_p.eq.1) then
        temp = ps_l
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*po_v(ibegin+p_off:ibegin+p_off+5)
        endif
        call set_dirichlet_bc(imatrix,ibegin+p_off,rhs,temp,normal,curv,izonedim)
     end if
  end do

end subroutine boundary_pres


!=======================================================
! boundary_dc
! ~~~~~~~~~~~
!
! sets homogeneous dirichlet boundary condition
!=======================================================
subroutine boundary_dc(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_dc called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call entdofs(1, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin, normal, curv, rhs)

     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim)
  end do

end subroutine boundary_dc


!=======================================================
! boundary_nm
! ~~~~~~~~~~~
!
! sets homogeneous Neumann boundary condition
!=======================================================
subroutine boundary_nm(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_nm called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call entdofs(1, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin, normal, curv, rhs)

     call set_normal_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim)
  end do

end subroutine boundary_nm


!=======================================================
! boundary_gs
! ~~~~~~~~~~~
!
! sets boundary conditions on psi in the GS solver
!=======================================================
subroutine boundary_gs(imatrix, rhs, feedfac)
  use basic
  use arrays
  use gradshafranov

  implicit none
  
  integer, intent(in) :: imatrix
  real, intent(in) :: feedfac
  vectype, intent(inout), dimension(*) :: rhs
  
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, ineg, ier
  real :: normal(2), curv
  real :: x, z,  alx, alz
  real, dimension(6) :: g
  real, dimension(1) :: xp, zp, xc, zc
  double precision :: coords(3)
  logical :: is_boundary
  vectype, dimension(6) :: temp
  real, dimension(2) :: temp1, temp2
  include 'mpif.h'

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_gs called"

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)

  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call assign_local_pointers(i)
     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin, normal, curv, rhs)

!......add feedback field
     if(idevice .eq. 0 .and. ifixedb .eq. 0) then
        xp(1) = x
        zp(1) = z
        xc(1) = 102.
        zc(1) = 10.
        call gvect(xp,zp,xc,zc,1,g,1,ineg)
        psis_l = psis_l + g*feedfac
     endif

     temp = psis_l
     if(ifixedb.ge.1) temp=0.
     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim)

    ! no toroidal current
    temp = 0.
    call set_laplacian_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim,-x)
  end do

end subroutine boundary_gs


!=======================================================
! boundary_vor
! ~~~~~~~~~~~~
!
! sets boundary conditions on Delta*(phi) 
! in the smoother
!=======================================================
subroutine boundary_vor(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vor called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call assign_local_pointers(i)
     call entdofs(numvarsm, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin, normal, curv, rhs)
     call rotate_matrix(imatrix, ibegin+6, normal, curv, rhs)

     if(inonormalflow.eq.1) then
        call set_dirichlet_bc(imatrix,ibegin+6,rhs,temp,normal,curv,izonedim)
     end if
     
     if(inoslip_pol.eq.1) then
        call set_normal_bc(imatrix,ibegin+6,rhs,temp,normal,curv,izonedim)
     end if

     ! no vorticity
     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim)
     call set_laplacian_bc(imatrix,ibegin+6,rhs,temp,normal,curv,izonedim,-x)
  end do

end subroutine boundary_vor

!=======================================================
! boundary_com
! ~~~~~~~~~~~~
!
! sets boundary conditions on Del^2(chi) in the smoother
!=======================================================
subroutine boundary_com(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  vectype, intent(inout), dimension(*) :: rhs
  
  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_com called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     call assign_local_pointers(i)
     call entdofs(numvarsm, i, 0, ibegin, iendplusone)
     call rotate_matrix(imatrix, ibegin, normal, curv, rhs)
     call rotate_matrix(imatrix, ibegin+6, normal, curv, rhs)

     ! clamp compression
     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,curv,izonedim)

     if(inonormalflow.eq.1) then
        call set_normal_bc(imatrix,ibegin+6,rhs,temp,normal,curv,izonedim)
     end if

     if(inoslip_pol.eq.1) then
        call set_dirichlet_bc(imatrix,ibegin+6,rhs,temp,normal,curv,izonedim)
     end if

     ! no compression
     if(com_bc.eq.1) then
        call set_laplacian_bc(imatrix,ibegin+6,rhs,temp,normal,curv,izonedim,x)
     endif
  end do

end subroutine boundary_com

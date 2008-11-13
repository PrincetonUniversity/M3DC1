!======================================================================
! boundary_node
!
! determines if node is on boundary, and returns relevant info
! about boundary surface
!======================================================================
subroutine boundary_node(inode,is_boundary,izone,izonedim,normal,x,z)
  use basic

  implicit none

  integer, intent(in) :: inode             ! node index
  integer, intent(out) :: izone,izonedim   ! zone type/dimension
!cj  real, intent(out) :: normal              ! outward normal of boundary
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real, intent(out) :: x,z                 ! coordinates of inode
  logical, intent(out) :: is_boundary      ! is inode on boundary

  integer :: ibottom, iright, ileft, itop
  double precision :: coords(3)

!cj 0: rect mesh; 1: curve mesh
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
  if(izone.eq.iright) then
!cj     normal = 0.
     normal(1) = 1.  !cos
     normal(2) = 0.  !sin
  else if(izone.eq.ileft) then
!cj     normal = pi
     normal(1) =-1.  !cos
     normal(2) = 0.  !sin
  else if(izone.eq.itop) then
!cj     normal = pi/2.
     normal(1) = 0.  !cos
     normal(2) = 1.  !sin
  else if(izone.eq.ibottom) then
!cj     normal = -pi/2.
     normal(1) = 0.  !cos
     normal(2) =-1.  !sin
  endif

          case(1)
  call nodNormalVec(inode, normal, is_boundary) 
  if(.not.is_boundary) return
  izonedim=1      !cj set izonedim always 1 for the shaped boundary
  izone=1         !cj dummy for shaped boundary

! if(myrank.eq.0) print *,"You are working with curved mesh."
          end select 

  call xyznod(inode,coords)
          select case(nonrect)
          case(0)
  x = coords(1) + xzero
  z = coords(2) + zzero
          case(1)
  x = coords(1) !cjdebug + xzero
  z = coords(2) !cjdebug + zzero
! if(myrank.eq.0) print *,"You are working with curved mesh."
          end select
end subroutine boundary_node


!======================================================================
! set_dirichlet_bc
!======================================================================
subroutine set_dirichlet_bc(imatrix,ibegin,rhs,bv,normal,izonedim)
  use basic
  implicit none

  integer, intent(in) :: imatrix              ! matrix handle
  integer, intent(in) :: ibegin               ! first dof of field
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(6) :: bv     ! boundary values
!cj  real, intent(in) :: normal                  ! angle of normal from horizontal
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  integer, intent(in) :: izonedim             ! dimension of boundary
  

  if(myrank.eq.0 .and. iprint.ge.2) print *, "set_dirichlet_bc called"

  if(izonedim.eq.1) then
     ! edge points
     ! ~~~~~~~~~~~
     !clamp value
     if(imatrix.ne.0) call setdiribc(imatrix, ibegin)
     rhs(ibegin) = bv(1)
     
     ! clamp tangential derivative
     call set_tangent_bc(imatrix,ibegin,rhs,bv,normal,izonedim)

  else if(izonedim.eq.0) then
     ! corner points
     ! ~~~~~~~~~~~~~
     if(imatrix.ne.0) then 
        call setdiribc(imatrix, ibegin)
        call setdiribc(imatrix, ibegin+1)
        call setdiribc(imatrix, ibegin+2)
        call setdiribc(imatrix, ibegin+3)
        call setdiribc(imatrix, ibegin+5)
     endif
     rhs(ibegin  ) = bv(1)
     rhs(ibegin+1) = bv(2)
     rhs(ibegin+2) = bv(3)
     rhs(ibegin+3) = bv(4)
     rhs(ibegin+5) = bv(6)
  endif
end subroutine set_dirichlet_bc


!======================================================================
! set_tangent_bc
!======================================================================
subroutine set_tangent_bc(imatrix,ibegin,rhs,bv,normal,izonedim)
  use basic

  implicit none
  
  integer, intent(in) :: imatrix              ! matrix handle
  integer, intent(in) :: ibegin               ! first dof of field
!cj  real, intent(in) :: normal                  ! angle of normal from horizontal
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(6) :: bv     ! boundary values
  integer, intent(in) :: izonedim             ! dimension of boundary

  integer :: irow, irow_n, irow_m
  integer :: numvals
  integer, dimension(2) :: cols
  vectype, dimension(2) :: vals


!cjdebug  return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "set_tangent_bc called"

  numvals = 2
!cj  vals(1) = -sin(normal)
!cj  vals(2) =  cos(normal)
  vals(1) = -normal(2)
  vals(2) =  normal(1)

  if(izonedim.eq.1) then
     ! edge points
     ! ~~~~~~~~~~~
     ! clamp tangential 1st-derivative
     cols(1) = ibegin + 1
     cols(2) = ibegin + 2

!cj     if(abs(cos(normal)) .gt. abs(sin(normal))) then
     if(abs(normal(1)) .gt. abs(normal(2))) then
        irow = ibegin + 2
     else
        irow = ibegin + 1
     endif
        irow_n = ibegin + 1
        irow_m = ibegin + 2
     
          select case(nonrect)
          case(0) 
     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = vals(1)*bv(2) + vals(2)*bv(3)
          case(1)
     if(imatrix.ne.0) then
          call applyLinCombinationForMatrix(imatrix, irow_n, irow_m, normal, icomplex) 
          call setgeneralbc(imatrix, irow_m, numvals, cols, vals, icomplex)
          rhs(irow_n) = normal(1)*rhs(irow_n) + normal(2)*bv(irow_m)   !cj normal
     endif
     rhs(irow_m) = -normal(2)*bv(2) + normal(1)*bv(3)    !cj tangent
!    if(myrank.eq.0) print *,"You are working with curved mesh."
          end select
     if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "set_tangent_bc : ibegin irow_n irow_m", &
                     ibegin, irow_n, irow_m

     ! clamp tangential 2nd-derivative
     cols(1) = ibegin + 3
     cols(2) = ibegin + 5
     
!cj     if(abs(cos(normal)) .gt. abs(sin(normal))) then
     if(abs(normal(1)) .gt. abs(normal(2))) then
        irow = ibegin + 5
     else
        irow = ibegin + 3
     endif
        irow_n = ibegin + 3
        irow_m = ibegin + 5

          select case(nonrect)
          case(0) 
     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = vals(1)*bv(4) + vals(2)*bv(6)
          case(1)
     if(imatrix.ne.0) then
          call applyLinCombinationForMatrix(imatrix, irow_n, irow_m, normal, icomplex) 
          call setgeneralbc(imatrix, irow_m, numvals, cols, vals, icomplex)
          rhs(irow_n) = normal(1)*rhs(irow_n) + normal(2)*bv(irow_m)   !cj normal
     endif
     rhs(irow_m) = -normal(2)*bv(4) + normal(1)*bv(6)    !cj tangent
!    if(myrank.eq.0) print *,"You are working with curved mesh."
          end select
     if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "set_tangent_bc : ibegin irow_n irow_m", &
                     ibegin, irow_n, irow_m

  else if(izonedim.eq.0) then
     ! corner points
     ! ~~~~~~~~~~~~~
     ! clamp tangential 1st-derivative
     cols(1) = ibegin + 1
     cols(2) = ibegin + 2

     irow = ibegin + 1    
     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = vals(1)*bv(2) + vals(2)*bv(3)

     irow = ibegin + 2
     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = vals(1)*bv(2) + vals(2)*bv(3)

     ! clamp tangential 2nd-derivative
     cols(1) = ibegin + 3
     cols(2) = ibegin + 5

     irow = ibegin + 3
     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = vals(1)*bv(4) + vals(2)*bv(6)

     irow = ibegin + 5
     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = vals(1)*bv(4) + vals(2)*bv(6)
  endif

end subroutine set_tangent_bc


!======================================================================
! set_normal_bc
!======================================================================
subroutine set_normal_bc(imatrix,ibegin,rhs,bv,normal,izonedim)
  use basic

  implicit none
  
  integer, intent(in) :: imatrix              ! matrix handle
  integer, intent(in) :: ibegin               ! first dof of field
!cj  real, intent(in) :: normal                  ! angle of normal from horizontal
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(6) :: bv     ! boundary values
  integer, intent(in) :: izonedim             ! dimension of boundary

  integer :: irow, irow_n, irow_m
  integer :: numvals
  integer, dimension(2) :: cols
  vectype, dimension(2) :: vals


!cjdebug  return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "set_normal_bc called"

  if(izonedim.eq.1) then
     ! edge points
     ! ~~~~~~~~~~~
     numvals = 2
     cols(1) = ibegin + 1
     cols(2) = ibegin + 2
!cj     vals(1) = cos(normal)
!cj     vals(2) = sin(normal)
     vals(1) = normal(1)
     vals(2) = normal(2)

!cj     if(abs(cos(normal)) .gt. abs(sin(normal))) then
     if(abs(normal(1)) .gt. abs(normal(2))) then
        irow = ibegin + 1
     else
        irow = ibegin + 2
     endif
        irow_n = ibegin + 1
        irow_m = ibegin + 2

     if(imatrix.ne.0) then
          select case(nonrect)
          case(0)
        call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
        call setdiribc(imatrix, ibegin+4)
          case(1)
        call applyLinCombinationForMatrix(imatrix, irow_n, irow_m, vals, icomplex) 
        call setgeneralbc(imatrix, irow_n, numvals, cols, vals, icomplex)
        call setdiribc(imatrix, ibegin+4)
        rhs(irow_m) = -vals(2)*rhs(irow_n) + vals(1)*rhs(irow_m)  !cj tangent
!     if(myrank.eq.0) print *,"You are working with curved mesh."
          end select
     endif
          select case(nonrect)
          case(0)
     rhs(irow) = vals(1)*bv(2) + vals(2)*bv(3)
          case(1)
     rhs(irow_n) = vals(1)*bv(2) + vals(2)*bv(3)   !cj normal
!     if(myrank.eq.0) print *,"You are working with curved mesh."
          end select
     rhs(ibegin+4) = bv(5)
     if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "set_normal_bc : ibegin irow_n irow_m", &
                     ibegin, irow_n, irow_m

  else if(izonedim.eq.0) then
     ! corner points
     ! ~~~~~~~~~~~~~
     if(imatrix.ne.0) then
        call setdiribc(imatrix, ibegin+1)
        call setdiribc(imatrix, ibegin+2)
        call setdiribc(imatrix, ibegin+4)
     end if
     rhs(ibegin+1) = bv(2)
     rhs(ibegin+2) = bv(3)
     rhs(ibegin+4) = bv(5)
  endif

end subroutine set_normal_bc


!======================================================================
! set_laplacian_bc
!======================================================================
subroutine set_laplacian_bc(imatrix,ibegin,rhs,bv,normal,izonedim,radius)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field
!cj  real, intent(in) :: normal         ! angle of normal vector from horizontal
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  vectype, intent(inout), dimension(*) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(6) :: bv     ! boundary values
  integer, intent(in) :: izonedim    ! dimension of boundary
  real, intent(in) :: radius         ! radial coordinate of node
                                     ! for gs operator, let radius -> -radius
  
  integer :: numvals, irow
  integer, dimension(3) :: cols
  vectype, dimension(3) :: vals

!cj 9/17/08
  if(nonrect .eq. 1) return

  if(myrank.eq.0 .and. iprint.ge.2) print *, "set_laplacian_bc called"

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
!cj     if(abs(cos(normal)) .gt. abs(sin(normal))) then
     if(abs(normal(1)) .gt. abs(normal(2))) then
        irow = ibegin + 3
     else
        irow = ibegin + 5
     endif

     if(imatrix.ne.0) &
          call setgeneralbc(imatrix, irow, numvals, cols, vals, icomplex)
     rhs(irow) = bv(1)
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
!cj  real :: normal, x, z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z
  vectype, dimension(6) :: temp

  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return 
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vel called"

  call numnod(numnodes)
  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call entdofs(vecsize_vel, i, 0, ibegin, iendplusone)
     call assign_local_pointers(i)

     ! no normal flow
     if(inonormalflow.eq.1) then
        temp = 0.
        call set_dirichlet_bc(imatrix,ibegin+u_off,rhs,temp,normal,izonedim)
        if(numvar.ge.3) then
           call set_normal_bc(imatrix,ibegin+chi_off,rhs,temp,normal,izonedim)
        endif
     end if
     
     ! no poloidal slip
     if(inoslip_pol.eq.1) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+u_off,rhs,temp,normal,izonedim)
        if(numvar.ge.3) then
           call set_dirichlet_bc(imatrix,ibegin+chi_off,rhs,temp,normal,izonedim)
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
           call set_dirichlet_bc(imatrix,ibegin+vz_off,rhs,temp,normal,izonedim)
        end if
        
        ! no toroidal stress
        if(inostress_tor.eq.1) then
           temp = 0.
           call set_normal_bc(imatrix,ibegin+vz_off,rhs,temp,normal,izonedim)
        end if
     endif
       
     ! no vorticity
     temp = 0.
     call set_laplacian_bc(imatrix,ibegin+u_off,rhs,temp,normal,izonedim,-x)

     ! no compression
     if(com_bc.eq.1 .and. numvar.ge.3) then
        call set_laplacian_bc(imatrix,ibegin+chi_off,rhs,temp,normal,izonedim,x)
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
!cj  real :: normal, x, z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z
  vectype, dimension(6) :: temp

  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_mag called"

  call numnod(numnodes)

  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call entdofs(vecsize_phi, i, 0, ibegin, iendplusone)
     call entdofs(1, i, 0, ibegin1, iendplusone1)
     call assign_local_pointers(i)

     ! clamp poloidal field
     temp = psis_l
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp = 1.5*temp + 0.5*psio_v(ibegin+psi_off:ibegin+psi_off+5)
     endif
     call set_dirichlet_bc(imatrix,ibegin+psi_off,rhs,temp,normal,izonedim)

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
        call set_dirichlet_bc(imatrix,ibegin+bz_off,rhs,temp,normal,izonedim)
     endif

     ! no toroidal current
     if(inocurrent_tor.eq.1) then
        temp = 0.
        if(jadv.eq.1 .and. igauge.eq.0) then
           temp(1) = vloop/(2.*pi*resistivity(ibegin1))
        endif
        call set_laplacian_bc(imatrix,ibegin+psi_off,rhs,temp,normal,izonedim,-x)
     end if

     ! no tangential current
     if(inocurrent_pol.eq.1 .and. numvar.ge.2) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+bz_off,rhs,temp,normal,izonedim)
     end if

     if(numvar.ge.3) then 
        if(inograd_p.eq.1) then
           temp = 0.
           call set_normal_bc(imatrix,ibegin+pe_off,rhs,temp,normal,izonedim)
        end if
        if(iconst_p.eq.1) then
           temp = pes_l
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*peo_v(ibegin+pe_off:ibegin+pe_off+5)
           endif
           call set_dirichlet_bc(imatrix,ibegin+pe_off,rhs,temp,normal,izonedim)
        else if(iconst_t.eq.1) then
           temp = pes_l*den1_l(1)/dens_l(1)
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*peo_v(ibegin+pe_off:ibegin+pe_off+5)
           endif
           call set_dirichlet_bc(imatrix,ibegin+pe_off,rhs,temp,normal,izonedim)
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
!cj  real :: normal,x,z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x,z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_den called"

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call entdofs(vecsize_n, i, 0, ibegin, iendplusone)
     call assign_local_pointers(i)

     if(inograd_n.eq.1) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+den_off,rhs,temp,normal,izonedim)
     end if
     if(iconst_n.eq.1) then
        temp = dens_l
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*deno_v(ibegin+den_off:ibegin+den_off+5)
        endif
        call set_dirichlet_bc(imatrix,ibegin+den_off,rhs,temp,normal,izonedim)
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
!cj  real :: normal, x, z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_pres called"

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call entdofs(vecsize_p, i, 0, ibegin, iendplusone)
     call assign_local_pointers(i)

     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(imatrix,ibegin+p_off,rhs,temp,normal,izonedim)
     end if
     if(iconst_p.eq.1) then
        temp = ps_l
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*po_v(ibegin+p_off:ibegin+p_off+5)
        endif
        call set_dirichlet_bc(imatrix,ibegin+p_off,rhs,temp,normal,izonedim)
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
!cj  real :: normal, x, z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_dc called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call entdofs(1, i, 0, ibegin, iendplusone)

     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,izonedim)
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
!cj  real :: normal, x, z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_nm called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call entdofs(1, i, 0, ibegin, iendplusone)

     call set_normal_bc(imatrix,ibegin,rhs,temp,normal,izonedim)
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
  integer :: ibegin, iendplusone, numnodes, ineg
!cj  real :: normal, x, z, xmin, zmin, alx, alz
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z, xmin, zmin, alx, alz
  real, dimension(6) :: g
  real, dimension(1) :: xp, zp, xc, zc
  double precision :: coords(3)
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_gs called"

  call getboundingboxsize(alx, alz)

  call numnod(numnodes)
  call getmincoord(xmin,zmin)
  if(myrank.eq.0 .and. iprint.ge.2) print *, "xmin  zmin = ", xmin,zmin

!cj added 10/06/08
!cjdebug  call dump_matrix(imatrix)

  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call assign_local_pointers(i)
     call entdofs(numvargs, i, 0, ibegin, iendplusone)

!......add feedback field
     if(idevice .eq. 0) then
        call xyznod(i,coords)
          select case(nonrect)
          case(0)
        xp = coords(1) - xmin + xzero
        zp = coords(2) - zmin + zzero
        xc(1) = 102.
        zc(1) = xzero + alx/2.
          case(1)
        xp = coords(1) !cjdebug - xmin + xzero
        zp = coords(2) !cjdebug - zmin + zzero
        xc(1) = 102.
        zc(1) = rzero + alx/2.
!       if(myrank.eq.0) print *,"You are working with curved mesh."
          end select
        call gvect(xp,zp,xc,zc,1,g,1,ineg)
        psis_l = psis_l + g*feedfac
     endif

     ! clamp magnetic field
     temp = psis_l
     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,izonedim)

     ! no toroidal current
     temp = 0.
     call set_laplacian_bc(imatrix,ibegin,rhs,temp,normal,izonedim,-x)
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
!cj  real :: normal, x, z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vor called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call assign_local_pointers(i)
     call entdofs(numvarsm, i, 0, ibegin, iendplusone)

     if(inonormalflow.eq.1) then
        call set_dirichlet_bc(imatrix,ibegin+6,rhs,temp,normal,izonedim)
     end if
     
     if(inoslip_pol.eq.1) then
        call set_normal_bc(imatrix,ibegin+6,rhs,temp,normal,izonedim)
     end if

     ! no vorticity
     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,izonedim)
     call set_laplacian_bc(imatrix,ibegin+6,rhs,temp,normal,izonedim,-x)    
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
!cj  real :: normal, x, z
  real, dimension(2) :: normal          !cj normal is changed from scalar to array
  real :: x, z
  logical :: is_boundary
  vectype, dimension(6) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_com called"

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,x,z)
     if(.not.is_boundary) cycle

     call assign_local_pointers(i)
     call entdofs(numvarsm, i, 0, ibegin, iendplusone)

     ! clamp compression
     call set_dirichlet_bc(imatrix,ibegin,rhs,temp,normal,izonedim)

     if(inonormalflow.eq.1) then
        call set_normal_bc(imatrix,ibegin+6,rhs,temp,normal,izonedim)
     end if

     if(inoslip_pol.eq.1) then
        call set_dirichlet_bc(imatrix,ibegin+6,rhs,temp,normal,izonedim)
     end if

     ! no compression
     if(com_bc.eq.1) then
        call set_laplacian_bc(imatrix,ibegin+6,rhs,temp,normal,izonedim,x)
     endif
  end do

end subroutine boundary_com

subroutine boundary_vel(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, irow
  real :: normal, x
  real, dimension(6) :: temp
  double precision :: coords(3)

  if(iper.eq.1 .and. jper.eq.1) return

  call getmodeltags(ibottom, iright, itop, ileft)

  temp = 0.

  call numnod(numnodes)
  do i=1, numnodes
     call zonenod(i,izone,izonedim)

     ! skip interior points
     if(izonedim.ge.2) cycle

     ! for periodic bc's
     ! skip if on a periodic boundary
     ! and convert corner to an edge when one edge is periodic
     if(iper.eq.1) then 
        if(izone.eq.ileft .or. izone.eq.iright) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
     endif
     if(jper.eq.1) then 
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
     endif

     call entdofs(numvar, i, 0, ibegin, iendplusone)
     call xyznod(i,coords)
     x = coords(i) + xzero

     ! fix stream function gauge
     call setdiribc(imatrix, ibegin)
     rhs(ibegin) = 0.

     ! edges
     if(izonedim.eq.1) then
        if(izone.eq.ileft) then
           normal = 0.
        else if(izone.eq.iright) then
           normal = pi
        else if(izone.eq.itop) then
           normal = pi/2.
        else if(izone.eq.ibottom) then
           normal = -pi/2.
        endif
  
        ! no normal flow
        call boundary_tangential_deriv(imatrix, ibegin, normal, rhs, temp)
        if(numvar.ge.3) then
           call boundary_normal_deriv(imatrix, ibegin+12, normal, irow)
           rhs(irow) = 0.
        endif
   
        ! no normal stress
!!$        call boundary_laplacian(imatrix, ibegin, normal, -x, irow)
!!$        rhs(irow) = 0.
!!$        if(numvar.ge.2) then
!!$           call boundary_normal_deriv(imatrix, ibegin+6, normal, irow)
!!$           rhs(irow) = 0.
!!$        endif
!!$        if(numvar.ge.3) then
!!$           call boundary_laplacian(imatrix, ibegin+12, normal, x, irow)
!!$           rhs(irow) = 0.
!!$        endif

     ! corners
     else if(izonedim.eq.0) then

        ! no normal flow
        call boundary_tangential_deriv(imatrix, ibegin, 0., rhs, temp)
        call boundary_tangential_deriv(imatrix, ibegin, pi/2., rhs, temp)
        if(numvar.ge.3) then
           call boundary_normal_deriv(imatrix, ibegin+12, 0., irow)
           rhs(irow) = 0.
           call boundary_normal_deriv(imatrix, ibegin+12, pi/2., irow)
           rhs(irow) = 0.
        endif

        ! no normal stress
!!$        call boundary_laplacian(imatrix, ibegin, 0., -x, irow)
!!$        rhs(irow) = 0.
!!$        call boundary_laplacian(imatrix, ibegin, pi/2., -x, irow)
!!$        rhs(irow) = 0.
!!$        if(numvar.ge.2) then
!!$           call boundary_normal_deriv(imatrix, ibegin+6, 0., irow)
!!$           rhs(irow) = 0.
!!$           call boundary_normal_deriv(imatrix, ibegin+6, pi/2., irow)
!!$           rhs(irow) = 0.
!!$        endif
!!$        if(numvar.ge.3) then
!!$           call boundary_laplacian(imatrix, ibegin+12, 0., x, irow)
!!$           rhs(irow) = 0.
!!$           call boundary_laplacian(imatrix, ibegin+12, pi/2., x, irow)
!!$           rhs(irow) = 0.
!!$        endif
     endif
        
  end do

end subroutine boundary_vel


subroutine boundary_mag(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, irow
  real :: normal, x
  real, dimension(6) :: temp
  double precision :: coords(3)

  if(iper.eq.1 .and. jper.eq.1) return

  call getmodeltags(ibottom, iright, itop, ileft)

  call numnod(numnodes)

  do i=1, numnodes
     call zonenod(i,izone,izonedim)

     ! skip interior points
     if(izonedim.ge.2) cycle

     ! for periodic bc's
     ! skip if on a periodic boundary
     ! and convert corner to an edge when one edge is periodic
     if(iper.eq.1) then 
        if(izone.eq.ileft .or. izone.eq.iright) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
     endif
     if(jper.eq.1) then 
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
     endif

     call entdofs(numvar, i, 0, ibegin, iendplusone)
     call xyznod(i,coords)
     x = coords(i) + xzero

     if(izonedim.eq.1) then
        if(izone.eq.ileft) then
           normal = 0.
        else if(izone.eq.iright) then
           normal = pi
        else if(izone.eq.itop) then
           normal = pi/2.
        else if(izone.eq.ibottom) then
           normal = -pi/2.
        endif

        ! clamp magnetic field at boundary
        temp = phis(ibegin:ibegin+5)
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*phiold(ibegin:ibegin+5)
        endif
        call boundary_clamp(imatrix, ibegin, normal, rhs, temp)
        if(numvar.ge.2) then
           temp = phis(ibegin+6:ibegin+11)
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*phiold(ibegin+6:ibegin+11)
           endif
           call boundary_clamp(imatrix, ibegin+6, normal, rhs, temp)
        endif
        if(numvar.ge.3) then
           temp = phis(ibegin+12:ibegin+17)
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*phiold(ibegin+12:ibegin+17)
           endif
           call boundary_clamp(imatrix, ibegin+12, normal, rhs, temp)
        endif

        ! no surface currents
!!$        call boundary_laplacian(imatrix, ibegin, normal, -x, irow)
!!$        rhs(irow) = 0.
!!$        if(numvar.ge.2) then
!!$           call boundary_normal_deriv(imatrix, ibegin+6, normal, irow)
!!$           rhs(irow) = 0.
!!$        endif

        ! no pressure gradient
!!$        if(numvar.ge.3) then
!!$           call boundary_normal_deriv(imatrix, ibegin+12, normal, irow)
!!$           rhs(irow) = 0.
!!$        endif

        ! add loop voltage
        if(integrator.eq.1 .and. ntime.gt.1) then
           rhs(ibegin) = rhs(ibegin) + 1.5*fbound
        else
           rhs(ibegin) = rhs(ibegin) + fbound
        endif
         

     else if(izonedim.eq.0) then

        ! clamp magnetic field
        call boundary_clamp_all(imatrix, ibegin)
        rhs(ibegin:ibegin+5) = phis(ibegin:ibegin+5)
        if(integrator.eq.1 .and. ntime.gt.1) then
           rhs(ibegin:ibegin+5) = 1.5*rhs(ibegin:ibegin+5) + 0.5*phiold(ibegin:ibegin+5)
        endif
        if(numvar.ge.2) then
           call boundary_clamp_all(imatrix, ibegin+6)
           rhs(ibegin+6:ibegin+11) = phis(ibegin+6:ibegin+11)
           if(integrator.eq.1 .and. ntime.gt.1) then
              rhs(ibegin+6:ibegin+11) = 1.5*rhs(ibegin+6:ibegin+11) + &
                   0.5*phiold(ibegin+6:ibegin+11)
           endif
        endif
        if(numvar.ge.3) then
           call boundary_clamp_all(imatrix, ibegin+12)
           rhs(ibegin+12:ibegin+17) = phis(ibegin+12:ibegin+17)
           if(integrator.eq.1 .and. ntime.gt.1) then
              rhs(ibegin+12:ibegin+17) = 1.5*rhs(ibegin+12:ibegin+17) + &
                   0.5*phiold(ibegin+12:ibegin+17)
           endif
        endif

        ! no surface currents
!!$        call boundary_laplacian(imatrix, ibegin, 0., -x, irow)
!!$        rhs(irow) = 0.
!!$        call boundary_laplacian(imatrix, ibegin, pi/2., -x, irow)
!!$        rhs(irow) = 0.
!!$        if(numvar.ge.2) then
!!$           call boundary_normal_deriv(imatrix, ibegin+6, 0., irow)
!!$           rhs(irow) = 0.
!!$           call boundary_normal_deriv(imatrix, ibegin+6, pi/2., irow)
!!$           rhs(irow) = 0.
!!$        endif

        ! no pressure gradient
!!$        if(numvar.ge.3) then
!!$           call boundary_normal_deriv(imatrix, ibegin+12, 0., irow)
!!$           rhs(irow) = 0.
!!$           call boundary_normal_deriv(imatrix, ibegin+12, pi/2., irow)
!!$           rhs(irow) = 0.
!!$        endif


        ! add loop voltage
        if(integrator.eq.1 .and. ntime.gt.1) then
           rhs(ibegin) = rhs(ibegin) + 1.5*fbound
        else
           rhs(ibegin) = rhs(ibegin) + fbound
        endif

     endif
        
  end do

end subroutine boundary_mag


subroutine boundary_gs(imatrix, rhs)
  use basic
  use arrays
  use gradshafranov

  implicit none
  
  integer, intent(in) :: imatrix
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, ibeginn, iendplusonen, numnodes, irow
  real :: normal, x
  real, dimension(6) :: temp
  double precision :: coords(3)

  if(iper.eq.1 .and. jper.eq.1) return

  call getmodeltags(ibottom, iright, itop, ileft)

  call numnod(numnodes)

  do i=1, numnodes
     call zonenod(i,izone,izonedim)

     ! skip interior points
     if(izonedim.ge.2) cycle

     ! for periodic bc's
     ! skip if on a periodic boundary
     ! and convert corner to an edge when one edge is periodic
     if(iper.eq.1) then 
        if(izone.eq.ileft .or. izone.eq.iright) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
     endif
     if(jper.eq.1) then 
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
     endif

     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call entdofs(numvar, i, 0, ibeginn, iendplusonen)
     call xyznod(i,coords)
     x = coords(i) + xzero

!!$     if(imatrix.ne.0) call setdiribc(imatrix, ibegin)
!!$     rhs(ibegin) = phis(ibeginn)

!!$     call boundary_clamp_all(imatrix, ibegin)
!!$     rhs(ibegin:ibegin+5) = phis(ibeginn:ibeginn+5)

     if(izonedim.eq.1) then
        if(izone.eq.ileft) then
           normal = 0.
        else if(izone.eq.iright) then
           normal = pi
        else if(izone.eq.itop) then
           normal = pi/2.
        else if(izone.eq.ibottom) then
           normal = -pi/2.
        endif

        ! clamp magnetic field at boundary
        temp = phis(ibeginn:ibeginn+5)
        call boundary_clamp(imatrix, ibegin, normal, rhs, temp)

        ! no surface currents
!!$        call boundary_laplacian(imatrix, ibegin, normal+pi/2., -x, irow)
!!$        rhs(irow) = 0.
     else if(izonedim.eq.0) then

        ! clamp magnetic field
        call boundary_clamp_all(imatrix, ibegin)
        rhs(ibegin:ibegin+5) = phis(ibeginn:ibeginn+5)

!!$        ! no surface currents
!!$        call boundary_laplacian(imatrix, ibegin, 0., -x, irow)
!!$        rhs(irow) = 0.
!!$        call boundary_laplacian(imatrix, ibegin, pi/2., -x, irow)
!!$        rhs(irow) = 0.
     endif        
  end do

end subroutine boundary_gs


subroutine boundary_sm(imatrix, rhs)
  use basic
  use arrays

  implicit none
  
  integer, intent(in) :: imatrix
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, irow
  real :: normal, x
  double precision :: coords(3)

  if(iper.eq.1 .and. jper.eq.1) return

  call getmodeltags(ibottom, iright, itop, ileft)

  call numnod(numnodes)

  do i=1, numnodes
     call zonenod(i,izone,izonedim)

     ! skip interior points
     if(izonedim.ge.2) cycle

     ! for periodic bc's
     ! skip if on a periodic boundary
     ! and convert corner to an edge when one edge is periodic
     if(iper.eq.1) then 
        if(izone.eq.ileft .or. izone.eq.iright) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
     endif
     if(jper.eq.1) then 
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
     endif

     call entdofs(numvar, i, 0, ibegin, iendplusone)
     call xyznod(i,coords)
     x = coords(i) + xzero

     ! no normal stress
     call setdiribc(imatrix, ibegin)
     rhs(ibegin) = 0.       
  end do

end subroutine boundary_sm


! Clamps boudary values to specified values
subroutine boundary_clamp_all(imatrix, ibegin)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field

  integer :: irow

  if(imatrix.eq.0) return

  ! clamp value
  call setdiribc(imatrix, ibegin)
  call setdiribc(imatrix, ibegin+1)
  call setdiribc(imatrix, ibegin+2)
  call setdiribc(imatrix, ibegin+3)
  call setdiribc(imatrix, ibegin+4)
  call setdiribc(imatrix, ibegin+5)

end subroutine boundary_clamp_all

! Clamps boudary values to specified values
subroutine boundary_clamp(imatrix, ibegin, normal, rhs, bv)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field
  real, intent(in) :: normal         ! angle of normal vector from horizontal
  real, intent(in), dimension(6) :: bv        ! boundary values
  real, intent(inout), dimension(*) :: rhs    ! rhs vector

  ! clamp value
  if(imatrix.ne.0) call setdiribc(imatrix, ibegin)
  rhs(ibegin) = bv(1)

  ! clamp tangential derivative
  call boundary_tangential_deriv(imatrix, ibegin, normal, rhs, bv)

end subroutine boundary_clamp

! Sets tangential derivative boundary condition
subroutine boundary_tangential_deriv(imatrix, ibegin, normal, rhs, bv)

  use basic

  implicit none

  integer, intent(in) :: imatrix             ! matrix handle
  integer, intent(in) :: ibegin              ! first dof of field
  real, intent(in) :: normal                 ! angle of normal vector from horizontal
  real, intent(inout), dimension(*) :: rhs   ! right-hand-side of matrix equation
  real, intent(in), dimension(6) :: bv       ! boundary values

  integer :: irow
  integer :: numvals
  integer, dimension(2) :: cols
  real, dimension(2) :: vals

  ! clamp tangential 1st-derivative
  numvals = 2
  cols(1) = ibegin + 1
  cols(2) = ibegin + 2
  vals(1) = -sin(normal)
  vals(2) =  cos(normal)

  if(abs(cos(normal)) .gt. abs(sin(normal))) then
     irow = ibegin + 2
  else
     irow = ibegin + 1
  endif

  if(imatrix.ne.0) call setgeneralbc(imatrix, irow, numvals, cols, vals)
  rhs(irow) = cos(normal)*bv(3) - sin(normal)*bv(2)

  ! clamp tangential 2nd-derivative
  numvals = 2
  cols(1) = ibegin + 3
  cols(2) = ibegin + 5
  vals(1) = -sin(normal)
  vals(2) =  cos(normal)

  if(abs(cos(normal)) .gt. abs(sin(normal))) then
     irow = ibegin + 5
  else
     irow = ibegin + 3
  endif

  if(imatrix.ne.0) call setgeneralbc(imatrix, irow, numvals, cols, vals)
  rhs(irow) = cos(normal)*bv(6) - sin(normal)*bv(4)

end subroutine boundary_tangential_deriv


! Sets normal derivative boundary condition
subroutine boundary_normal_deriv(imatrix, ibegin, normal, irow)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field
  real, intent(in) :: normal         ! angle of normal vector from horizontal
  integer, intent(out) :: irow       ! matrix row replaced by bc

  integer :: numvals
  integer, dimension(2) :: cols
  real, dimension(2) :: vals

  numvals = 2
  cols(1) = ibegin + 1
  cols(2) = ibegin + 2
  vals(1) = cos(normal)
  vals(2) = sin(normal)

  if(abs(cos(normal)) .gt. abs(sin(normal))) then
     irow = ibegin + 1
  else
     irow = ibegin + 2
  endif

  if(imatrix.ne.0) call setgeneralbc(imatrix, irow, numvals, cols, vals)
  
end subroutine boundary_normal_deriv


! Sets laplacian boundary condition
!  (for grad-shafranov operator, radius -> -radius).
subroutine boundary_laplacian(imatrix, ibegin, normal, radius, irow)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field
  real, intent(in) :: normal         ! angle of normal vector from horizontal
  real, intent(in) :: radius         ! radial coordinate of node
  integer, intent(out) :: irow       ! matrix row replaced by bc
  
  integer :: numvals
  integer, dimension(3) :: cols
  real, dimension(3) :: vals

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

  if(abs(cos(normal)) .gt. abs(sin(normal))) then
     irow = ibegin + 5
  else
     irow = ibegin + 3
  endif

  if(imatrix.ne.0) call setgeneralbc(imatrix, irow, numvals, cols, vals)
  
end subroutine










! Dirichlet Boundary Conditions
! =============================
subroutine boundaryds(ibound,nbc,jsymtype)
  use basic
  implicit none

  ! jsymtype = 1 for even up-down symmetry
  !          = 2 for odd   "        "

  integer, dimension(*), intent(out) :: ibound
  integer, intent(in) :: jsymtype
  integer, intent(out) :: nbc

  integer :: numnodes, i, izone, izonedim
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone

  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)
  nbc = 0

  if(iper.eq.1 .and. jper.eq.1) return

  do i=1,numnodes
     call entdofs(1, i, 0, ibegin, iendplusone)
     call zonenod(i,izone,izonedim)
     if(izonedim .eq. 1) then

        ! left/right boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        if(iper.eq.0 .and. (izone.eq.ileft .or. izone.eq.iright)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+5
           nbc =  nbc+3

           ibound(nbc+1) = ibegin+1
           ibound(nbc+2) = ibegin+3
           ibound(nbc+3) = ibegin+4
           nbc = nbc+3
           
        ! top/bottom boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        else if(jper.eq.0 .and. (izone.eq.itop .or. izone.eq.ibottom)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+1
           ibound(nbc+3) = ibegin+3
           nbc =  nbc+3

           ibound(nbc+1) = ibegin+2
           ibound(nbc+2) = ibegin+4
           ibound(nbc+3) = ibegin+5
           nbc = nbc+3

        endif

        ! corner points
        ! ~~~~~~~~~~~~~
     else if(izonedim .eq. 0) then
        ibound(nbc+1) = ibegin
        nbc =  nbc+1

        ibound(nbc+1) = ibegin+1
        ibound(nbc+2) = ibegin+2
        ibound(nbc+3) = ibegin+3
        ibound(nbc+4) = ibegin+4
        ibound(nbc+5) = ibegin+5
        nbc = nbc+5

!!$        if(iper.eq.0) then
!!$           ibound(nbc+1) = ibegin+2
!!$           ibound(nbc+2) = ibegin+5
!!$           nbc =  nbc+2
!!$        endif
!!$        if(jper.eq.0) then
!!$           ibound(nbc+1) = ibegin+1
!!$           ibound(nbc+2) = ibegin+3
!!$           nbc =  nbc+2
!!$        endif
!!$        if(iper.eq.0 .and. jper.eq.0) then
!!$           ibound(nbc+1) = ibegin+4
!!$           nbc =  nbc+1
!!$        endif
     endif
  enddo
  return
end subroutine boundaryds



! Pressure Boundary Conditions
! ============================
subroutine boundarypres(ibound,nbc)
  use basic
  use arrays
  implicit none
  
  integer :: numnodes, i, izone, ibound(*), nbc, izonedim
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone

  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)

  nbc = 0

  if(iper.eq.1 .and. jper.eq.1) return

  do i=1,numnodes
     call entdofs(1, i, 0, ibegin, iendplusone)
     call zonenod(i,izone,izonedim)
     if(izonedim .eq. 1) then

        ! left/right boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        if(iper.eq.0 .and. (izone.eq.ileft .or. izone.eq.iright)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+5
           nbc =  nbc+3
           
        ! top/bottom boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        else if(jper.eq.0 .and. (izone.eq.itop .or. izone .eq. ibottom)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+1
           ibound(nbc+3) = ibegin+3
           nbc =  nbc+3
        end if

        ! corner points
        ! ~~~~~~~~~~~~~
     else if(izonedim .eq. 0) then
        ibound(nbc+1) = ibegin
        nbc =  nbc+1
        if(iper.eq.0) then
           ibound(nbc+1) = ibegin+2
           ibound(nbc+2) = ibegin+5
           nbc =  nbc+2
        endif
        if(jper.eq.0) then
           ibound(nbc+1) = ibegin+1
           ibound(nbc+2) = ibegin+3
           nbc =  nbc+2
        endif
     endif
  enddo

  return
end subroutine boundarypres


! Smoother boundary conditions
! ============================
subroutine boundarysm(ibound,ibound2,nbc,iplace)
  use basic
  implicit none
  
  integer, dimension(*) :: ibound, ibound2
  integer, intent(in) :: iplace
  integer, intent(out) :: nbc

  integer :: numnodes, i, izone, izonedim
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone, numvarsm

  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)

  numvarsm = 2
  nbc = 0

  if(iper.eq.1 .and. jper.eq.1) return

  do i=1,numnodes
     call entdofs(numvarsm, i, 0, ibegin, iendplusone)
     call zonenod(i,izone,izonedim)
     if(izonedim .eq. 1) then

        ! left/right boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        if(iper.eq.0 .and. (izone.eq.ileft .or. izone.eq.iright)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+5
           nbc =  nbc+3
           if(iplace.lt.3) then
              ibound(nbc+1) = ibegin+6
              ibound(nbc+2) = ibegin+8
              ibound(nbc+3) = ibegin+9
              ibound(nbc+4) = ibegin+11
              nbc =  nbc+4
           else
              ibound(nbc+1) = ibegin+7
              ibound(nbc+2) = ibegin+10
              nbc =  nbc+2
           endif
           
        ! top/bottom boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        else if(jper.eq.0 .and. (izone.eq.itop .or. izone .eq. ibottom)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+1
           ibound(nbc+3) = ibegin+3
           nbc =  nbc+3

           if(iplace.lt.3) then
              ibound(nbc+1) = ibegin+6
              ibound(nbc+2) = ibegin+7
              ibound(nbc+3) = ibegin+9
              ibound(nbc+4) = ibegin+11
              nbc =  nbc+4
           else
              ibound(nbc+1) = ibegin+8
              ibound(nbc+2) = ibegin+9
              ibound2(nbc+2)= ibegin+11
              ibound(nbc+3) = ibegin+10
              nbc =  nbc+3
           endif
        endif

        ! corner points
        ! ~~~~~~~~~~~~~
     else if(izonedim .eq. 0) then
        ibound(nbc+1) = ibegin
        nbc =  nbc+1
        if(iper.eq.0) then
           ibound(nbc+1) = ibegin+2
           ibound(nbc+2) = ibegin+5
           nbc = nbc +2
        endif
        if(jper.eq.0) then
           ibound(nbc+1) = ibegin+1
           ibound(nbc+2) = ibegin+3
           nbc =  nbc+2
        endif

        if(iplace.lt.3) then
           ibound(nbc+1) = ibegin+6
           ibound(nbc+2) = ibegin+9
           ibound(nbc+3) = ibegin+11
           nbc =  nbc+3
           if(iper.eq.0) then
              ibound(nbc+1) = ibegin+8
              nbc =  nbc+1
           endif
           if(jper.eq.0) then
              ibound(nbc+1) = ibegin+7
              nbc =  nbc+1
           endif
        else
           ibound(nbc+1) = ibegin+10
           nbc = nbc+1
           if(iper.eq.0) then
              ibound(nbc+1) = ibegin+7
              nbc =  nbc+1
           endif
           if(jper.eq.0) then
              ibound(nbc+1) = ibegin+8
              ibound(nbc+2) = ibegin+9
              ibound2(nbc+2)= ibegin+11
              nbc =  nbc+2
           endif
        endif
     endif
  enddo

  return

end subroutine boundarysm

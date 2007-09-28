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
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

     call entdofs(vecsize, i, 0, ibegin, iendplusone)
     call xyznod(i,coords)
     x = coords(1) + xzero

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
        temp = 0.
        call boundary_clamp(imatrix, ibegin+phi_off, normal, rhs, temp)
        if(numvar.ge.3) then
           call boundary_normal_deriv(imatrix, ibegin+chi_off, normal, rhs, temp)
!!$           if(itor.eq.0) then
!!$              if(imatrix.ne.0) call setdiribc(imatrix, ibegin+chi_off+4)
!!$              rhs(ibegin+chi_off+4) = 0.
!!$           endif
        endif

!!$        ! no velocity 
!!$        call boundary_normal_deriv(imatrix, ibegin+phi_off, normal, rhs, temp)
!!$        if(imatrix.ne.0) call setdiribc(imatrix, ibegin+phi_off+4)
!!$        rhs(ibegin+phi_off+4) = 0.

        ! clamp toroidal velocity
        if(numvar.ge.2) then
           select case(v_bc)
           case(1)               ! no normal stress
              temp = 0.
              call boundary_normal_deriv(imatrix, ibegin+vz_off, normal, rhs, temp)
              if(imatrix.ne.0) call setdiribc(imatrix, ibegin+vz_off+4)
              rhs(ibegin+vz_off+4) = 0.

           case default          ! no slip
              temp = vzs_v(ibegin+vz_off:ibegin+vz_off+5)
              if(integrator.eq.1 .and. ntime.gt.1) then
                 temp = 1.5*temp + 0.5*vzo_v(ibegin+vz_off:ibegin+vz_off+5)
              endif
              call boundary_clamp(imatrix, ibegin+vz_off, normal, rhs, temp)
           end select
        endif
       
        ! no vorticity
        call boundary_laplacian(imatrix, ibegin+phi_off, normal, -x, irow)
        rhs(irow) = 0.

        if(numvar.ge.3) then
           ! no compression
           if(com_bc.eq.1) then
              call boundary_laplacian(imatrix, ibegin+chi_off, normal, x, irow)
              rhs(irow) = 0.
           endif
        endif
   
     ! corners
     else if(izonedim.eq.0) then

        ! no normal flow
        temp = 0.
        call boundary_clamp_all(imatrix, ibegin+phi_off, rhs, temp)
        if(numvar.ge.3) then
           call boundary_normal_deriv(imatrix, ibegin+chi_off, 0., rhs, temp)
           call boundary_normal_deriv(imatrix, ibegin+chi_off, pi/2., rhs, temp)
           if(imatrix.ne.0) call setdiribc(imatrix, ibegin+chi_off+4)
           rhs(ibegin+chi_off+4) = 0.
        endif

        if(numvar.ge.2) then
           select case (v_bc)
           case(1)               ! no normal stress
              call boundary_normal_deriv(imatrix, ibegin+vz_off, 0., rhs, temp)
              call boundary_normal_deriv(imatrix, ibegin+vz_off, pi/2., rhs, temp)
              if(imatrix.ne.0) call setdiribc(imatrix, ibegin+vz_off+4)
              rhs(ibegin+vz_off+4) = 0.           

           case default          ! no-slip         
              temp = vzs_v(ibegin+vz_off:ibegin+vz_off+5)
              if(integrator.eq.1 .and. ntime.gt.1) then
                 temp = 1.5*temp + 0.5*vzo_v(ibegin+vz_off:ibegin+vz_off+5)
              endif
              call boundary_clamp_all(imatrix, ibegin+vz_off, rhs, temp)
           end select
        end if
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
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, irow
  real :: normal, x, z
  real, dimension(6) :: temp
  double precision :: coords(3)

  if(iper.eq.1 .and. jper.eq.1) return

  call getmodeltags(ibottom, iright, itop, ileft)

  call numnod(numnodes)

  do i=1, numnodes
     call zonenod(i,izone,izonedim)

     call entdofs(vecsize, i, 0, ibegin, iendplusone)
     call xyznod(i,coords)
     x = coords(1) + xzero
     z = coords(2) + zzero

     ! skip interior points
     if(izonedim.ge.2) cycle
!!$     if(izonedim.ge.2) then
!!$        ! limiter
!!$        if(numvar.ge.3 .and. itor.eq.1 .and. &
!!$             ((xlim.lt.xmag .and. x.lt.xlim) .or. &
!!$             (xlim.gt.xmag .and. x.gt.xlim))) then
!!$           temp = phis(ibegin+12:ibegin+17)
!!$           if(integrator.eq.1 .and. ntime.gt.1) then
!!$              temp = 1.5*temp + 0.5*phiold(ibegin+12:ibegin+17)
!!$           endif
!!$           call boundary_clamp(imatrix, ibegin+12, normal, rhs, temp)
!!$        else
!!$           cycle
!!$        endif
!!$     endif

     ! for periodic bc's
     ! skip if on a periodic boundary
     ! and convert corner to an edge when one edge is periodic
     if(iper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

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

        ! clamp poloidal field
        temp = psis_v(ibegin+psi_off:ibegin+psi_off+5)
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*psio_v(ibegin+psi_off:ibegin+psi_off+5)
        endif
        call boundary_clamp(imatrix, ibegin+psi_off, normal, rhs, temp)

        ! clamp toroidal field
        if(numvar.ge.2) then
           temp = bzs_v(ibegin+bz_off:ibegin+bz_off+5)
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*bzo_v(ibegin+bz_off:ibegin+bz_off+5)
           endif
           call boundary_clamp(imatrix, ibegin+bz_off, normal, rhs, temp)
        endif  

        ! no toroidal current
        call boundary_laplacian(imatrix, ibegin+psi_off, normal, -x, irow)
        rhs(irow) = 0.

        ! no tangential current
        if(numvar.ge.2) then
           temp = 0.
           call boundary_normal_deriv(imatrix, ibegin+bz_off, normal, rhs, temp)
           if(imatrix.ne.0) call setdiribc(imatrix, ibegin+bz_off+4)
           rhs(ibegin+bz_off+4) = 0.
        endif

        if(numvar.ge.3) then 
           select case (p_bc)
           case(1)       ! no normal pressure gradient (insulating)
              temp = 0.
              call boundary_normal_deriv(imatrix, ibegin+pe_off, normal, rhs, temp)
              if(imatrix.ne.0) call setdiribc(imatrix, ibegin+pe_off+4)
              rhs(ibegin+pe_off+4) = 0.

           case default  ! clamp pressure
              temp = pes_v(ibegin+pe_off:ibegin+pe_off+5)
              if(integrator.eq.1 .and. ntime.gt.1) then
                 temp = 1.5*temp + 0.5*peo_v(ibegin+pe_off:ibegin+pe_off+5)
              endif
              call boundary_clamp(imatrix, ibegin+pe_off, normal, rhs, temp)
           end select
        endif

        ! add loop voltage
        if(integrator.eq.1 .and. ntime.gt.1) then
           rhs(ibegin+psi_off) = rhs(ibegin+psi_off) + 1.5*fbound
        else
           rhs(ibegin+psi_off) = rhs(ibegin+psi_off) + fbound
        endif
         

     else if(izonedim.eq.0) then

        ! clamp poloidal field
        temp = psis_v(ibegin+psi_off:ibegin+psi_off+5)
        if(integrator.eq.1 .and. ntime.gt.1) then
           temp = 1.5*temp + 0.5*psio_v(ibegin+psi_off:ibegin+psi_off+5)
        endif
        call boundary_clamp_all(imatrix, ibegin+psi_off, rhs, temp)

        ! clamp toroidal field
        if(numvar.ge.2) then
           temp = bzs_v(ibegin+bz_off:ibegin+bz_off+5)
           if(integrator.eq.1 .and. ntime.gt.1) then
              temp = 1.5*temp + 0.5*bzo_v(ibegin+bz_off:ibegin+bz_off+5)
           endif
           ! no tangential current
           temp(ibegin+bz_off+1) = 0.
           temp(ibegin+bz_off+2) = 0.
           if(imatrix.ne.0) call setdiribc(imatrix, ibegin+bz_off+4)
           rhs(ibegin+bz_off+4) = 0.
           ! clamp field
           call boundary_clamp_all(imatrix, ibegin+bz_off, rhs, temp)
        endif


        if(numvar.ge.3) then
           select case (p_bc)
           case(1)      ! no normal pressure gradient (insulating)
              temp = 0.
              call boundary_normal_deriv(imatrix, ibegin+pe_off, 0., rhs, temp)
              call boundary_normal_deriv(imatrix, ibegin+pe_off, pi/2., rhs, temp)
              if(imatrix.ne.0) call setdiribc(imatrix, ibegin+pe_off+4)
              rhs(ibegin+pe_off+4) = 0.
              
           case default ! clamp pressure
              temp = pes_v(ibegin+pe_off:ibegin+pe_off+5)
              if(integrator.eq.1 .and. ntime.gt.1) then
                 temp = 1.5*temp + 0.5*phiold(ibegin+pe_off:ibegin+pe_off+5)
              endif
              call boundary_clamp_all(imatrix, ibegin+pe_off, rhs, temp)
           end select
        endif

        ! add loop voltage
        if(integrator.eq.1 .and. ntime.gt.1) then
           rhs(ibegin+psi_off) = rhs(ibegin+psi_off) + 1.5*fbound
        else
           rhs(ibegin+psi_off) = rhs(ibegin+psi_off) + fbound
        endif

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
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, irow
  real :: normal
  real, dimension(6) :: temp

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
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

     call entdofs(vecsize1, i, 0, ibegin, iendplusone)

     temp = dens_v(ibegin+den_off:ibegin+den_off+5)
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp = 1.5*temp + 0.5*deno_v(ibegin+den_off:ibegin+den_off+5)
     endif

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

        ! clamp density
        call boundary_clamp(imatrix, ibegin+den_off, normal, rhs, temp)
!!$        temp = 0.
!!$        call boundary_normal_deriv(imatrix, ibegin, normal, rhs, temp)
!!$        if(imatrix.ne.0) call setdiribc(imatrix, ibegin+4)
!!$        rhs(ibegin+4) = 0.

        
     ! corners
     else if(izonedim.eq.0) then
        ! clamp density
        call boundary_clamp_all(imatrix, ibegin+den_off, rhs, temp)

!!$        temp = 0.
!!$        call boundary_normal_deriv(imatrix, ibegin, 0., rhs, temp)
!!$        call boundary_normal_deriv(imatrix, ibegin, pi/2., rhs, temp)
!!$        if(imatrix.ne.0) call setdiribc(imatrix, ibegin+4)
!!$        rhs(ibegin+4) = 0.
     endif
        
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
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, irow
  real :: normal
  real, dimension(6) :: temp

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
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

     call entdofs(vecsize1, i, 0, ibegin, iendplusone)

     temp = ps_v(ibegin+p_off:ibegin+p_off+5)
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp = 1.5*temp + 0.5*po_v(ibegin+p_off:ibegin+p_off+5)
     endif

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

        ! clamp pressure
        call boundary_clamp(imatrix, ibegin+p_off, normal, rhs, temp)

     ! corners
     else if(izonedim.eq.0) then
        ! clamp pressure
        call boundary_clamp_all(imatrix, ibegin+p_off, rhs, temp)
     endif
        
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
  real, intent(inout), dimension(*) :: rhs
  
  integer :: ibottom, iright, itop, ileft, i, izone, izonedim
  integer :: ibegin, iendplusone, numnodes, irow
  real :: normal
  real, dimension(6) :: temp

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
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

     call entdofs(1, i, 0, ibegin, iendplusone)

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

        ! homogeneous dirichlet boundary conditions
        call boundary_clamp(imatrix, ibegin, normal, rhs, temp)

     ! corners
     else if(izonedim.eq.0) then
        ! clamp density
        call boundary_clamp_all(imatrix, ibegin, rhs, temp)
     endif
        
  end do

end subroutine boundary_dc

!=======================================================
! boundary_gs
! ~~~~~~~~~~~
!
! sets boundary conditions on psi in the GS solver
!=======================================================
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
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call entdofs(vecsize, i, 0, ibeginn, iendplusonen)
     call xyznod(i,coords)
     x = coords(1) + xzero

     temp = psis_v(ibeginn+psi_off:ibeginn+psi_off+5)

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
        call boundary_clamp(imatrix, ibegin, normal, rhs, temp)

        ! no toroidal current
        call boundary_laplacian(imatrix, ibegin, normal, -x, irow)
        rhs(irow) = 0.

     else if(izonedim.eq.0) then
        
        ! clamp magnetic field
        call boundary_clamp_all(imatrix, ibegin, rhs, temp)
        
     endif        
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
  real, intent(inout), dimension(*) :: rhs
  
  integer, parameter :: numvarsm = 2
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
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

     call entdofs(numvarsm, i, 0, ibegin, iendplusone)
     call xyznod(i,coords)
     x = coords(1) + xzero

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

        ! clamp vorticity
        call boundary_clamp(imatrix, ibegin, normal, rhs, temp)

        ! clamp stream function
        call boundary_clamp(imatrix, ibegin+6, normal, rhs, temp)

        ! no vorticity
        call boundary_laplacian(imatrix, ibegin+6, normal, -x, irow)
        rhs(irow) = 0.

     ! corners
     else if(izonedim.eq.0) then
        ! clamp vorticity
        call boundary_clamp_all(imatrix, ibegin, rhs, temp)

        ! clamp stream function
        call boundary_clamp_all(imatrix, ibegin+6, rhs, temp)
     endif
     
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
  real, intent(inout), dimension(*) :: rhs
  
  integer, parameter :: numvarsm = 2
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
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ibottom
        endif
        if(izone.eq.ileft .or. izone.eq.iright) cycle
     endif
     if(jper.eq.1) then 
        if(izonedim.eq.0) then 
           izonedim = 1
           izone = ileft
        endif
        if(izone.eq.ibottom .or. izone.eq.itop) cycle
     endif

     call entdofs(numvarsm, i, 0, ibegin, iendplusone)
     call xyznod(i,coords)
     x = coords(1) + xzero

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

        ! clamp compression
        call boundary_clamp(imatrix, ibegin, normal, rhs, temp)

        ! no normal flow
        call boundary_normal_deriv(imatrix, ibegin+6, normal, rhs, temp)
        if(imatrix.ne.0) call setdiribc(imatrix, ibegin+10)
        rhs(ibegin+10) = 0.

        ! no compression
        if(com_bc.eq.1) then
           call boundary_laplacian(imatrix, ibegin+6, normal, x, irow)
           rhs(irow) = 0.
        endif

     ! corners
     else if(izonedim.eq.0) then
        ! clamp compression
        call boundary_clamp_all(imatrix, ibegin, rhs, temp)

        ! no normal flow
        call boundary_normal_deriv(imatrix, ibegin+6, 0., rhs, temp)
        call boundary_normal_deriv(imatrix, ibegin+6, pi/2., rhs, temp)
        if(imatrix.ne.0) call setdiribc(imatrix, ibegin+10)
        rhs(ibegin+10) = 0.
     endif
  end do

end subroutine boundary_com


! Clamps boundary values to specified values
subroutine boundary_clamp_all(imatrix, ibegin, rhs, bv)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field
  real, intent(inout), dimension(*) :: rhs
  real, intent(in), dimension(6) :: bv

  integer :: irow

  ! clamp value
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


end subroutine boundary_clamp_all

! Clamps boundary values to specified values
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
subroutine boundary_normal_deriv(imatrix, ibegin, normal, rhs, bv)

  use basic

  implicit none

  integer, intent(in) :: imatrix     ! matrix handle
  integer, intent(in) :: ibegin      ! first dof of field
  real, intent(in) :: normal         ! angle of normal vector from horizontal
  real, intent(inout), dimension(*) :: rhs   ! right-hand-side of matrix equation
  real, intent(in), dimension(6) :: bv       ! matrix row replaced by bc
  
  integer :: numvals, irow
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
  rhs(irow) = cos(normal)*bv(2) + sin(normal)*bv(3)

!!$  if(imatrix.ne.0) call setdiribc(imatrix, ibegin+4)
!!$  rhs(ibegin+4) = bv(5)

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
     irow = ibegin + 3
  else
     irow = ibegin + 5
  endif

  if(imatrix.ne.0) call setgeneralbc(imatrix, irow, numvals, cols, vals)
  
end subroutine

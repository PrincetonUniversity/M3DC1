! Dirichlet Boundary Conditions
! =============================
subroutine boundaryds(ibound,nbc,iplace)
  use basic
  implicit none

  integer :: numnodes, i, izone, ibound(*), nbc, iplace, izonedim, j
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone

  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)
  nbc = 0
  if(iper .eq. 0) then 
     ! iper = 0
     ! ~~~~~~~~
     do i=1,numnodes
        call entdofs(1, i, 0, ibegin, iendplusone)
        call zonenod(i,izone, izonedim)
        if(izonedim .eq. 1) then ! an edge...
           ibound(nbc+1) =ibegin
           if(izone .eq. ibottom .or. izone .eq. itop) then
              ibound(nbc+2) = ibegin+1
              ibound(nbc+3) = ibegin+3
           else
              ibound(nbc+2) = ibegin+2
              ibound(nbc+3) = ibegin+5
           endif
           nbc = nbc+3
        else if(izonedim .eq. 0) then
           do j=0,5
              ibound(nbc+j+1) = ibegin+j
           enddo
           nbc = nbc + 6
        endif
     enddo
  else 
     ! iper = 1
     ! ~~~~~~~~
     do i=1,numnodes
        call zonenod(i,izone, izonedim)
        if(izonedim .le. 1) then
           call entdofs(1, i, 0, ibegin, iendplusone)
           if(izone .eq. itop .or. izone .eq. ibottom) then 
              ibound(nbc+1) =ibegin
              ibound(nbc+2) = ibegin+1
              ibound(nbc+3) = ibegin+3
              nbc = nbc+3
           else if(izonedim .eq. 0) then
              ibound(nbc+1) =ibegin
              ibound(nbc+2) = ibegin+1
              ibound(nbc+3) = ibegin+3                  
              nbc = nbc+3
           endif
        endif
     enddo
  endif
  
end subroutine boundaryds


! Velocity Boundary Conditions (no-slip)
! ======================================
subroutine boundaryv(ibound,nbc)
  use basic
  implicit none
  integer :: numnodes, i, izone, ibound(*), nbc, izonedim
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone

  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)
  nbc = 0

  if(iper.eq.1 .and. jper.eq.1) return

  do i=1,numnodes
     call entdofs(numvar, i, 0, ibegin, iendplusone)
     call zonenod(i,izone,izonedim)
     if(izonedim .eq. 1) then

        ! left/right boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        if(iper.eq.0 .and. (izone.eq.ileft .or. izone.eq.iright)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+3
           ibound(nbc+4) = ibegin+5
           nbc =  nbc+4
           if(numvar.ge.2) then
              ibound(nbc+1) = ibegin+6
              ibound(nbc+2) = ibegin+8
              ibound(nbc+3) = ibegin+11
              nbc =  nbc+3
           endif
           if(numvar.ge.3) then
              ibound(nbc+1) = ibegin+13
              ibound(nbc+2) = ibegin+16
              nbc =  nbc+2
           endif
           
        ! top/bottom boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        else if(jper.eq.0 .and. (izone.eq.itop .or. izone.eq.ibottom)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+1
           ibound(nbc+3) = ibegin+3
           ibound(nbc+4) = ibegin+5
           nbc =  nbc+4
           if(numvar.ge.2) then
              ibound(nbc+1) = ibegin+6
              ibound(nbc+2) = ibegin+7
              ibound(nbc+3) = ibegin+9
              nbc =  nbc+3
           endif
           if(numvar.ge.3) then
              ibound(nbc+1) = ibegin+14
              ibound(nbc+2) = ibegin+16
              nbc =  nbc+2
           endif
        endif

        ! corner points
        ! ~~~~~~~~~~~~~
     else if(izonedim .eq. 0) then
        ibound(nbc+1) = ibegin
        ibound(nbc+2) = ibegin+3
        ibound(nbc+3) = ibegin+5
        nbc =  nbc+3
        if(iper.eq.0) then
           ibound(nbc+1) = ibegin+2
           nbc =  nbc+1
        endif
        if(jper.eq.0) then
           ibound(nbc+1) = ibegin+1
           nbc =  nbc+1
        endif

        if(numvar.ge.2) then
           ibound(nbc+1) = ibegin+6
           nbc =  nbc+1
           if(iper.eq.0) then
              ibound(nbc+1) = ibegin+8
              ibound(nbc+2) = ibegin+11
              nbc =  nbc+2
           endif
           if(jper.eq.0) then
              ibound(nbc+1) = ibegin+7
              ibound(nbc+2) = ibegin+9
              nbc =  nbc+2
           endif
        endif

        if(numvar.ge.3) then
           ibound(nbc+1) = ibegin+16
           nbc =  nbc+1
           if(iper.eq.0) then
              ibound(nbc+1) = ibegin+13
              nbc =  nbc+1
           endif
           if(jper.eq.0) then
              ibound(nbc+1) = ibegin+14
              nbc =  nbc+1
           endif
        endif

     endif
  enddo
  return
end subroutine boundaryv


! Pressure Boundary Conditions
! ============================
subroutine boundaryp(ibound,nbc)
  use basic
  use arrays
  implicit none
  
  integer :: numnodes, i, izone, ibound(*), nbc, izonedim
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone

  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)
  psibounds = 0

  nbc = 0

  if(iper.eq.1 .and. jper.eq.1) return

  do i=1,numnodes
     call entdofs(numvar, i, 0, ibegin, iendplusone)
     call zonenod(i,izone,izonedim)
     if(izonedim .eq. 1) then

        ! left/right boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        if(iper.eq.0 .and. (izone.eq.ileft .or. izone.eq.iright)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+5
           nbc =  nbc+3
           if(numvar.ge.2) then
              psibounds(nbc+1) = gbound
              ibound(nbc+1) = ibegin+6
              ibound(nbc+2) = ibegin+8
              ibound(nbc+3) = ibegin+11
              nbc =  nbc+3
           endif
           if(numvar.ge.3) then
              ibound(nbc+1) = ibegin+12
              ibound(nbc+2) = ibegin+14
              ibound(nbc+3) = ibegin+17
              nbc =  nbc+3
           endif
           
        ! top/bottom boundaries
        ! ~~~~~~~~~~~~~~~~~~~~~
        else if(jper.eq.0 .and. (izone.eq.itop .or. izone .eq. ibottom)) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+1
           ibound(nbc+3) = ibegin+3
           nbc =  nbc+3
           if(numvar.ge.2) then
              psibounds(nbc+1) = gbound
              ibound(nbc+1) = ibegin+6
              ibound(nbc+2) = ibegin+7
              ibound(nbc+3) = ibegin+9
              nbc =  nbc+3
           endif
           if(numvar.ge.3) then
              ibound(nbc+1) = ibegin+12
              ibound(nbc+2) = ibegin+13
              ibound(nbc+3) = ibegin+15
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

        if(numvar.ge.2) then
           psibounds(nbc+1) = gbound
           ibound(nbc+1) = ibegin+6
           nbc = nbc+1
           if(iper.eq.0) then
              ibound(nbc+1) = ibegin+8
              ibound(nbc+2) = ibegin+11
              nbc =  nbc+2
           endif
           if(jper.eq.0) then
              ibound(nbc+1) = ibegin+7
              ibound(nbc+2) = ibegin+9
              nbc =  nbc+2
           endif
        endif

        if(numvar.ge.3) then
           ibound(nbc+1) = ibegin+12
           nbc = nbc+1
           if(iper.eq.0) then
              ibound(nbc+1) = ibegin+14
              ibound(nbc+2) = ibegin+17
              nbc =  nbc+2
           endif
           if(jper.eq.0) then
              ibound(nbc+1) = ibegin+13
              ibound(nbc+2) = ibegin+15
              nbc =  nbc+2
           endif
        endif

     endif
  enddo

  return
end subroutine boundaryp


! Smoother boundary conditions
! ============================
subroutine boundarysm(ibound,nbc,iplace)
  use basic
  implicit none
  
  integer :: numnodes, i, izone, ibound(*), nbc, izonedim, iplace
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone

  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)
  nbc = 0
  if(numvar.eq. 2) then
     do i=1,numnodes
        call entdofs(2, i, 0, ibegin, iendplusone)
        call zonenod(i,izone, izonedim)
        if(izonedim .eq. 1) then ! an edge...
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+6
           nbc = nbc+2
           if(izone .eq. ibottom .or. izone .eq. itop) then
              ibound(nbc+1) = ibegin+1
              ibound(nbc+2) = ibegin+3
              ibound(nbc+3) = ibegin+7
              ibound(nbc+4) = ibegin+9
              ibound(nbc+5) = ibegin+11
              nbc = nbc+5
           else
              ibound(nbc+1) = ibegin+2
              ibound(nbc+2) = ibegin+5
              ibound(nbc+3) = ibegin+8
              ibound(nbc+4) = ibegin+9
              ibound(nbc+5) = ibegin+11
              nbc = nbc+5               
           endif
        else if(izonedim .eq. 0) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+6
           nbc = nbc+2
           ibound(nbc+1) = ibegin+1
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+3
           ibound(nbc+4) = ibegin+5
           nbc = nbc+4
           ibound(nbc+1) = ibegin+7
           ibound(nbc+2) = ibegin+8
           ibound(nbc+3) = ibegin+9
           ibound(nbc+4) = ibegin+11
           nbc = nbc+4
        endif
     enddo
     return
  endif

  !     other bc's not implemented yet
  call safestop(8892) 
  return
end subroutine boundarysm

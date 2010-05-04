program struct2vac

  implicit none

  include 'mpif.h'

  integer :: ifirst, itotal, ilast, inode, ier
  integer, parameter :: ifile = 5
  character (len=*), parameter :: filename = 'ordered.points'

  integer :: next_node
  real, dimension(4) :: coords
  real :: normal(2)
  logical :: is_boundary

  call MPI_Init(ier)

  ! load mesh data
  call loadmesh("struct.dmg", "struct-dmg.sms")

  ! open file for output
  open(unit=ifile,file=filename,action='write',status='replace')

  ! find first boundary node and total number of boundary nodes
  call first_node(ifirst, itotal)
  print *, 'Found ', itotal, ' boundary nodes.  First node = ', ifirst

  ! write total number of nodes
  write(ifile, '(I8)') itotal

  ! write node entry for first node
  call xyznod(ifirst,coords)
  call nodNormalVec(ifirst,normal,is_boundary)

!!$  write(ifile, '(I8,4f12.8)') ifirst, coords(1), coords(2), &
!!$       normal(1), normal(2)
  write(ifile, '(I8,2f12.6)') ifirst, coords(1), coords(2)

  inode = ifirst
  do
     ilast = inode
     inode = next_node(ilast)    ! find next node

     ! if next node is the first node, we're done
     if(inode.eq.ifirst .or. inode.eq.-1) exit

     ! write node entry
     call xyznod(inode,coords)
     call nodNormalVec(inode,normal,is_boundary)

!!$     write(ifile, '(I8,4f12.8)') inode, coords(1), coords(2), &
!!$       normal(1), normal(2)
     write(ifile, '(I8,2f12.6)') inode, coords(1), coords(2)

  end do

  ! close output file
  close(ifile)

  call deletesearchstructure()
  call clearscorecdata()
  call MPI_finalize(ier)

end program struct2vac


logical function is_boundary(inode)
  implicit none

  integer, intent(in) :: inode
  integer :: izone, izonedim
  
  call zonenod(inode, izone, izonedim)

  is_boundary = (izonedim.le.1)
end function is_boundary


subroutine first_node(ifirst, itotal)
  implicit none

  integer, intent(out) :: ifirst, itotal

  integer :: numnodes, i
  logical :: is_boundary
  
  call numnod(numnodes)

  ifirst = -1
  itotal = 0

  do i=1, numnodes
     if(is_boundary(i)) then
        if(ifirst.eq.-1) ifirst = i
        itotal = itotal + 1
     endif
  end do
end subroutine first_node


integer function next_node(inode)

  implicit none

  integer, intent(in) :: inode

  integer :: ntri, i, j, k
  integer, dimension(4) :: id
  logical :: is_boundary

  call numfac(ntri)

  do i=1, ntri
     call nodfac(i,id)
     
     do j=1, 3
        if(id(j).eq.inode) then
           k = mod(j,3)+1
           if(is_boundary(id(k))) then
              next_node = id(k)
              return
           endif
        endif
     end do
  end do

  print *, 'Error: cannot find next node'
  next_node = -1
end function next_node

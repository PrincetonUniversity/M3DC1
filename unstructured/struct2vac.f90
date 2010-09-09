program struct2vac

  use mesh_mod

  implicit none

  include 'mpif.h'

  integer :: inode, ier, num_global_nodes, num_local_nodes
  integer :: i, j, ilast, ifirst, numnodes
  integer, parameter :: ifile = 5, inew = 6, ntor = 1
  character (len=*), parameter :: filename = 'ordered.points'

  integer :: next_node, myrank, maxrank
  real, dimension(4) :: coords
  real :: normal(2)
  logical :: is_boundary, myturn, first_time
  integer :: idum

#ifdef USESCOREC
  integer :: maxdofs1
#endif

  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank, ier)


#ifdef USESCOREC
  ! initialize autopack
  if(myrank.eq.0) print *, 'Initializing Autopack'
  call AP_INIT()
#endif

  ! load mesh data
  print *, 'loading mesh'
  call load_mesh

#ifdef USESCOREC
  call createdofnumbering(1, iper, jper, 6, 0, 0, 0, maxdofs1)
#endif

  numnodes = local_nodes()

  ! find the boundary node with the largest global id
  ifirst = -1
  do i=1, numnodes
     if(is_boundary_node(i)) then
        call entglobalid(i,0,inode)
        if(inode.gt.ifirst) ilast = inode
     endif
  end do

  first_time = .true.
  num_local_nodes = 0
  do
     call mpi_allreduce(ilast, inode, 1, MPI_INTEGER, &
          MPI_MAX, MPI_COMM_WORLD, ier)

     if(myrank.eq.0) then
        print *, "=====Boundary section starting at node", inode, "======="
     endif
 
     if(first_time) then
        ifirst = inode
     else if (ifirst.eq.inode) then
        exit
     end if

     myturn = .false.
     ! Check to see if I own current node
     do i=1, numnodes
        call entglobalid(i,0,j)
        if(j.eq.inode) then
           myturn = .true.
           exit
        end if
     end do
     
     ! Check to see if I own the next node
     if(myturn) then
        j = next_node(i)
        myturn = j.ne.-1
     endif
 
     ! If I do own the next node, I have to write the data
     if(myturn) then
        print *, 'node ', myrank, ' owns node ', inode

        if(first_time) then
           open(unit=ifile, file='scratch',action='write',status='replace')

           ! write toroidal mode number
           write(ifile, '(I8)') ntor

           ! write total number of nodes
           write(ifile, '(I8)') 0
        else
           open(unit=ifile,file='scratch',action='write', &
                status='old',access='append')           
        endif

        do
           ! We own this node, so write it to the list
           call xyznod(i,coords)
           call nodNormalVec(i,normal,is_boundary)
           write(ifile, '(I8,4f12.6)') inode, coords(1), coords(2), &
                normal(1), normal(2)           
           num_local_nodes = num_local_nodes + 1

           ! Move to next node
           i = j

           ! Check if we have arrived at the first global node
           call entglobalid(i,0,inode)
           if(inode.eq.ifirst) exit

           ! Check if we own the subsequent node
           ! If not, it is not our responsibility to write the current node
           j = next_node(i)
           if(j.eq.-1) exit
        end do
       
        close(ifile)

        print *, myrank, ' wrote ', num_local_nodes, ' nodes.'

        ilast = inode
     else
        ilast = -1
     endif
    
     first_time = .false.
  end do

  call mpi_reduce(num_local_nodes, num_global_nodes, 1, MPI_INTEGER, &
       MPI_SUM, 0, MPI_COMM_WORLD, ier)

  if(myrank.eq.0) then
     print *, 'Wrote ', num_global_nodes, ' nodes.'

     ! overwrite header
     ! ~~~~~~~~~~~~~~~~
     open(unit=ifile,file='scratch',action='readwrite',status='old', &
          disp='delete')
     open(unit=inew,file=filename,action='write',status='replace')

     ! write toroidal mode number
     read(ifile, '(I8)') idum
     write(inew, '(I8)') ntor

     ! number of global nodes
     read(ifile, '(I8)') idum
     write(inew, '(I8)') num_global_nodes

     ! re-write node data
     ! ~~~~~~~~~~~~~~~~~~
     do i=1, num_global_nodes
        read(ifile, '(I8,4f12.6)') inode, coords(1), coords(2), &
             normal(1), normal(2)           
        write(inew, '(I8,4f12.6)') inode, coords(1), coords(2), &
             normal(1), normal(2)           
     end do

     close(ifile)
     close(inew)

     ! delete scratch file
     ! ~~~~~~~~~~~~~~~~~~~
     
  endif

  call unload_mesh

  call MPI_finalize(ier)

end program struct2vac

!=======================================================
! next_node
! ~~~~~~~~~
! returns the local id of the boundary node that
! succeeds the boundary node with local id inode
!=======================================================
integer function next_node(inode)

  use element
  use mesh_mod

  implicit none

  integer, intent(in) :: inode

  integer :: ntri, i, j, k
  integer, dimension(nodes_per_element) :: id

  call numfac(ntri)

  do i=1, ntri
     call nodfac(i,id)
     
     do j=1, nodes_per_element
        if(id(j).eq.inode) then
           k = mod(j,nodes_per_element)+1
           if(is_boundary_node(id(k))) then
              next_node = id(k)
              return
           endif
        endif
     end do
  end do

  next_node = -1
end function next_node

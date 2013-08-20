module scorec_mesh_mod
  use element

  implicit none

  logical :: is_rectilinear

  real :: xzero, zzero
  real :: tiltangled
  integer :: iper, jper
  integer :: icurv
  integer :: nplanes

  integer, parameter :: maxi = 20

  logical, private :: has_bounding_box = .false.
  real, private :: bb(4)

  logical, parameter :: IGNORE_PHI = .true.
  logical, private :: initialized = .false.

  integer :: plane_comm

  character(len=256) mesh_model
  character(len=256) mesh_filename
  
  integer :: imulti_region

contains

  subroutine load_mesh()
    use math
    implicit none

#ifdef USE3D
    include 'mpif.h'
    real :: minphi, maxphi
    integer :: i,procs_per_plane, myrank, maxrank, ier, full_group, plane_group
    integer, allocatable :: ranks(:)
#endif

    ! initialize scorec solvers
    call initsolvers

    if(imulti_region.eq.1) then
       call create_tag_list(inner_wall, 2)
       inner_wall%tags(1) = 1
       inner_wall%tags(2) = 2
       call create_tag_list(outer_wall, 2)
       outer_wall%tags(1) = 3
       outer_wall%tags(2) = 4
       call create_tag_list(domain_boundary, 2)
       domain_boundary%tags(1) = 5
       domain_boundary%tags(2) = 6
!!$       call create_tag_list(all_boundaries, 2)
!!$       all_boundaries%tags(1) = 5
!!$       all_boundaries%tags(2) = 6
       call create_tag_list(all_boundaries, 6)
       all_boundaries%tags(1) = 1
       all_boundaries%tags(2) = 2
       all_boundaries%tags(3) = 3
       all_boundaries%tags(4) = 4
       all_boundaries%tags(5) = 5
       all_boundaries%tags(6) = 6
    else
       call create_tag_list(inner_wall, 1)
       inner_wall%tags(1) = 1
       call create_tag_list(outer_wall, 1)
       outer_wall%tags(1) = 1
       call create_tag_list(domain_boundary, 1)
       domain_boundary%tags(1) = 1
       call create_tag_list(all_boundaries, 1)
       all_boundaries%tags(1) = 1
    end if

    ! load mesh
#ifdef USE3D 
    call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ier)

    if(myrank.eq.0) print *, 'setting number of planes = ', nplanes
    call setTotNbPlane(nplanes)

    procs_per_plane = maxrank/nplanes
    if(myrank.eq.0) &
         print *, 'setting number of processes per plane = ', procs_per_plane
    call setNbProcPlane(procs_per_plane)

    minphi = 0.
    maxphi = twopi/nplanes * (nplanes-1)
    if(myrank.eq.0) print *, 'setting phi range', minphi, maxphi
    call setPhiRange(minphi, maxphi)

    if(myrank.eq.0) print *, 'loading partitioned mesh...'
    if(is_rectilinear) then
       if(myrank.eq.0) print *, 'loading rectilinear mesh model...'
       call loadPtnMesh('struct.dmg')
    else
       if(myrank.eq.0) print *, 'loading curved mesh model...'
       call loadPtnMesh('AnalyticModel')
    endif

    if(myrank.eq.0) print *, 'setting up 3D mesh...'
    call threeDMeshSetup(0)

    ! set up communications groups
    allocate(ranks(procs_per_plane))
    call MPI_Comm_group(MPI_COMM_WORLD, full_group, ier)
    do i=1, procs_per_plane
       ranks(i) = local_plane()*procs_per_plane + i-1
    end do
    call MPI_Group_incl(full_group, procs_per_plane, ranks, &
         plane_group, ier)
    call MPI_Comm_create(MPI_COMM_WORLD, plane_group, plane_comm, ier)
    deallocate(ranks)
#else 
    call loadmesh(trim(mesh_model), trim(mesh_filename))
#endif

    initialized = .true.
  end subroutine load_mesh

  subroutine unload_mesh
    if(.not. initialized) return

    call deletesearchstructure
    call clearscorecdata
!    call scorecfinalize
    call finalizesolvers
  end subroutine unload_mesh

  subroutine print_node_data()
    implicit none
    integer :: i,numnodes
    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv, x, z

    numnodes = local_nodes()

    do i=1, numnodes
       call boundary_node(i, is_boundary, izone, izonedim, normal, curv, x, z)
       if(.not.is_boundary) cycle
       write(*,'(5f12.4)') x, z, normal(1), normal(2), curv
    end do
  end subroutine print_node_data

  !======================================================================
  ! local_plane
  ! ~~~~~~~~~~~
  ! returns toroidal plane associated with current process
  !======================================================================
  integer function local_plane()
    implicit none
    call getplaneid(local_plane)
  end function local_plane

  !==============================================================
  ! local_nodes
  ! ~~~~~~~~~~~
  ! returns the number of elements local to this process
  !==============================================================
  integer function local_elements()
    implicit none
#ifdef USE3D
    call numwed(local_elements)
#else
    call numfac(local_elements)
#endif
  end function local_elements

  !==============================================================
  ! owned_nodes
  ! ~~~~~~~~~~~
  ! returns the number of elements owned by to this process
  !==============================================================
  integer function owned_nodes()
    implicit none
    call numnod(owned_nodes)
  end function owned_nodes


  !==============================================================
  ! local_nodes
  ! ~~~~~~~~~~~
  ! returns the number of nodes local to this process
  !==============================================================
  integer function local_nodes()
    implicit none
    call numnod(local_nodes)
  end function local_nodes


  !==============================================================
  ! global_node_id
  ! ~~~~~~~~~~~~~~
  !==============================================================
  integer function global_node_id(inode)
    implicit none
    integer, intent(in) :: inode
    call entglobalid(inode, 0, global_node_id)
  end function global_node_id

  !==============================================================
  ! get_element_nodes
  ! ~~~~~~~~~~~~~~~~~
  ! returns the indices of each element node
  !==============================================================
  subroutine get_element_nodes(ielm,n)
    implicit none
    integer, intent(in) :: ielm
    integer, intent(out), dimension(nodes_per_element) :: n

#ifdef USE3D
    integer :: nelms
    call numwed(nelms)
    if(ielm.gt.nelms) then
       print *, 'error!', ielm, nelms
       return
    end if

    call nodwed(ielm, n)
#else
    call nodfac(ielm, n)
#endif
  end subroutine get_element_nodes


  !==============================================================
  ! get_bounding_box
  ! ~~~~~~~~~~~~~~~~
  ! returns the lower-left and upper-right coordinates of the 
  ! mesh bounding box
  !==============================================================
  subroutine get_bounding_box(x1, z1, x2, z2)
    implicit none
    real, intent(out) :: x1, z1, x2, z2
    
    if(.not.has_bounding_box) then
       call getmincoord2(bb(1), bb(2))
       call getmaxcoord2(bb(3), bb(4))
       if(is_rectilinear) then
          bb(3) = bb(3) - bb(1) + xzero
          bb(4) = bb(4) - bb(2) + zzero
          bb(1) = xzero
          bb(2) = zzero
       endif
       has_bounding_box = .true.
    endif
    x1 = bb(1)
    z1 = bb(2)
    x2 = bb(3)
    z2 = bb(4) 
  end subroutine get_bounding_box


  !==============================================================
  ! get_bounding_box
  ! ~~~~~~~~~~~~~~~~
  ! returns the width and heightof the
  ! mesh bounding box
  !==============================================================
  subroutine get_bounding_box_size(alx, alz)
    implicit none

    real, intent(out) :: alx, alz
    real :: x1, z1, x2, z2
    call get_bounding_box(x1,z1,x2,z2)
    alx = x2 - x1
    alz = z2 - z1
  end subroutine get_bounding_box_size


  !============================================================
  ! get_node_pos
  ! ~~~~~~~~~~~~
  ! get the global coordinates (x,z) of node inode
  !============================================================
  subroutine get_node_pos(inode,x,phi,z)
    implicit none
    
    integer, intent(in) :: inode
    real, intent(out) :: x, phi, z
    
    real :: coords(3)
    real :: x1, z1
    
    call xyznod(inode,coords)
    x = coords(1)
    z = coords(2)

    if(is_rectilinear) then
       call getmincoord2(x1, z1)
       x = x + xzero - x1
       z = z + zzero - z1
    endif
#ifdef USE3D
    phi = coords(3)
#endif
  end subroutine get_node_pos


  !============================================================
  ! get_element_data
  ! ~~~~~~~~~~~~~~~~
  ! calculates element parameters a, b, c, and theta from
  ! node coordinates
  !============================================================
  subroutine get_element_data(itri, d)
    use math
    implicit none
    integer, intent(in) :: itri
    type(element_data), intent(out) :: d
    real :: x2, x3, z2, z3, phi2, phi3, x2p, x3p, z2p, z3p, hi
    integer :: nodeids(nodes_per_element)
    
    d%itri = itri

    call get_element_nodes(itri, nodeids)
    
    call get_node_pos(nodeids(1), d%R, d%phi, d%Z)
    call get_node_pos(nodeids(2), x2, phi2, z2)
    call get_node_pos(nodeids(3), x3, phi3, z3)
    x2 = x2 - d%R
    z2 = z2 - d%Z
    x3 = x3 - d%R
    z3 = z3 - d%Z
    
    hi = 1./sqrt(x2**2 + z2**2)
    d%co = x2*hi
    d%sn = z2*hi

    x2p =  d%co * x2 + d%sn * z2
    z2p = -d%sn * x2 + d%co * z2
    x3p =  d%co * x3 + d%sn * z3
    z3p = -d%sn * x3 + d%co * z3
    if(z2p .gt. 0.0001 .or. z2p .lt. -0.0001) then
       print *, "z2p should be 0.d0 but is ", z2p
    endif
    d%a = x2p-x3p
    d%b = x3p
    d%c = z3p
    if(d%c .le. 0.) then
       print *, 'ERROR: clockwise node ordering for element',itri
    endif

#ifdef USE3D
    call get_node_pos(nodeids(4), x2, phi2, z2)
    d%d = phi2 - d%Phi
    if(d%d .le. 0.) d%d = d%d + twopi
#else
    d%d = 0.
#endif
  end subroutine get_element_data


  !============================================================
  ! whattri
  ! ~~~~~~~
  ! Gets the element itri in which global coordinates (x,z)
  ! fall.  If no element on this process contains (x,z) then
  ! itri = -1.  Otherwise, xref and zref are global coordinates
  ! of the first node of element itri.
  !============================================================
  subroutine whattri(x,phi,z,ielm,xref,zref,nophi)
    implicit none
    
    real, intent(in) :: x, phi, z
    integer, intent(inout) :: ielm
    real, intent(out) :: xref, zref
    logical, intent(in), optional :: nophi
    
    integer :: nelms, i
    type(element_data) :: d

    ! first, try ielm
    if(ielm.gt.0) then
       call get_element_data(ielm,d)
       if(is_in_element(d,x,phi,z,nophi)) then
          xref = d%R
          zref = d%Z
          return
       endif
    endif

    ! try all elements
    ielm = -1
    nelms = local_elements()
    do i=1, nelms
       call get_element_data(i,d)
       if(is_in_element(d,x,phi,z,nophi)) then
          xref = d%R
          zref = d%Z
          ielm = i
          return
       end if
    end do
  end subroutine whattri


  !=========================================
  ! is_boundary_node
  ! ~~~~~~~~~~~~~~~~
  ! returns true if node lies on boundary
  !=========================================
  logical function is_boundary_node(inode, tags)
    implicit none
    integer, intent(in) :: inode
    type(tag_list), intent(in), optional :: tags
    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv, x, z

    call boundary_node(inode,is_boundary,izone,izonedim,normal,curv,x,z,tags)
    is_boundary_node = is_boundary
  end function is_boundary_node

  !======================================================================
  ! boundary_node
  ! ~~~~~~~~~~~~~
  ! determines if node is on boundary, and returns relevant info
  ! about boundary surface
  !======================================================================
  subroutine boundary_node(inode,is_boundary,izone,izonedim,normal,curv,x,z, &
       tags)

    use math
    
    implicit none
    
    integer, intent(in) :: inode              ! node index
    integer, intent(out) :: izone,izonedim    ! zone type/dimension
    real, intent(out) :: normal(2), curv
    real, intent(out) :: x,z                  ! coordinates of inode
    logical, intent(out) :: is_boundary       ! is inode on boundary
    type(tag_list), intent(in), optional :: tags
    
    integer :: ibottom, iright, ileft, itop, ib
    real :: angler, phi
    real, dimension(3) :: norm

    curv = 0.

    call zonenod(inode,izone,izonedim)

    if(is_rectilinear) then
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
          normal(1) = 1.
          normal(2) = 0.
       end if
       
    else 
       call nodNormalVec(inode, norm, ib)
       normal = norm(1:2)
       is_boundary = ib.eq.1
       if(.not.is_boundary) return

       if(icurv.eq.0) then
          curv = 0.
       else
          call nodcurvature2(inode, curv, ib)
          is_boundary = ib.eq.1
          if(.not.is_boundary) return
       end if

       if(present(tags)) then
          is_boundary = in_tag_list(tags, izone)         
       else
          is_boundary = in_tag_list(domain_boundary, izone)
       end if
    end if

    call get_node_pos(inode,x,phi,z)
  end subroutine boundary_node

  subroutine boundary_edge(itrin, is_edge, normal, idim)

    implicit none
    integer, intent(in) :: itrin
    integer, intent(out) :: is_edge(3)
    real, intent(out) :: normal(2,3)
    integer, intent(out) :: idim(3)
    
    integer :: inode(nodes_per_element), izone, i, j, jj
    integer :: in(2)
    real :: x, z, c(3)
    logical :: is_bound(3), found_edge

    integer :: iedge(3), izonedim, itri, ifaczone,ifaczonedim

#ifdef USE3D
    integer :: ifac(5), iplane
    call facregion(itrin,ifac)
    itri = 0
    do i=1, 5
       call locationplane(ifac(i), 2, iplane)
       if(iplane.eq.1) then
          itri = ifac(i)
          if(i.ne.1) print *, 'first face is not on plane.'
          exit
       end if
    end do
    if(itri.eq.0) print *, 'Error: no face found on plane!!!'
#else
    itri = itrin
#endif

    call nodfac(itri,inode)
    call edgfac(itri,iedge)
    call zonfac(itri,ifaczone,ifaczonedim)

    do i=1,3
       call boundary_node(inode(i),is_bound(i),izone,idim(i), &
            normal(:,i),c(i),x,z,all_boundaries)
    end do

    is_edge = 0
    do i=1,3
       call zonedg(iedge(i),izone,izonedim)
       if(izonedim.gt.1) cycle

       call nodedg(iedge(i), in)

       ! check nodes to see which one is the first node of iedge(i)
       found_edge = .false.
       do j=1, 3
          jj = mod(j,3)+1
          if(in(1).eq.inode(j) .and. in(2).eq.inode(jj)) then
             found_edge = .true.
             exit
          else if(in(2).eq.inode(j) .and. in(1).eq.inode(jj)) then
!             print *, 'Error: clockwise edge!'
             found_edge = .true.
             exit
          end if
       end do
       if(found_edge) then
          if(is_edge(j).ne.0) print *, 'Warning: edge counted twice'
          is_edge(j) = izone
!          call boundary_node(inode(j),is_bound(j),izone,idim(j), &
!               normal(:,j),c(j),x,z)
!          write(*, '(6F10.4)') x, z, normal(:,j)
       else
          print *, 'Error: phantom edge!'
       end if

       if(imulti_region.eq.1) then
          if(ifaczone.eq.2 .and. (izone.eq.1 .or. izone.eq.2)) normal=-normal
          if(ifaczone.eq.3 .and. (izone.eq.3 .or. izone.eq.4)) normal=-normal
       end if
    end do
  end subroutine boundary_edge

  subroutine get_zone(itri, izone)
    integer, intent(in) :: itri
    integer, intent(out) :: izone
    integer :: izonedim

    call zonfac(itri, izone, izonedim)
  end subroutine get_zone

end module scorec_mesh_mod

module scorec_mesh_mod
  use element
  use physical_mesh 

  implicit none

  logical :: is_rectilinear

  real :: xzero, zzero
  real :: tiltangled
  integer :: iper, jper
  integer :: icurv
  integer :: nplanes
  integer :: iread_planes

  integer, parameter :: maxi = 20

  logical, private :: has_bounding_box = .false.
  real, private :: bb(4)

  logical, parameter :: IGNORE_PHI = .true.
  logical, private :: initialized = .false.

  integer :: plane_comm

  character(len=256) mesh_model
  character(len=256) mesh_filename
  character(len=256) name_buff
  integer :: ipartitioned
  integer :: imatassemble  
  integer :: imulti_region
  integer :: iadapt_snap

  real :: toroidal_pack_factor
  real :: toroidal_pack_angle
  real :: toroidal_period

  integer, dimension (:), allocatable :: nodes_owned

  integer, parameter :: ZONE_ANY = -1
  integer, parameter :: ZONE_UNKNOWN = 0
  integer, parameter :: ZONE_PLASMA = 1
  integer, parameter :: ZONE_CONDUCTOR = 2
  integer, parameter :: ZONE_VACUUM = 3

  integer, parameter :: BOUND_ANY = -1
  integer, parameter :: BOUND_UNKNOWN = 0
  integer, parameter :: BOUND_FIRSTWALL = 1
  integer, parameter :: BOUND_DOMAIN = 2

  integer, parameter :: max_bounds = 20
  integer, parameter :: max_zones = 20
  integer :: boundary_type(max_bounds)
  integer :: zone_type(max_zones)

  integer :: offset_elms, global_elms, elms_per_plane

  ! adjacency matrix info
  integer, allocatable, dimension(:) :: n_adjacent ! num of adjacent elements
  integer, allocatable, dimension(:,:) :: adjacent ! adjacent elements
#ifdef USE3D
  integer, parameter :: max_adj = 5
#else
  integer, parameter :: max_adj = 3
#endif

contains

  subroutine load_mesh
    use math
    use read_ascii
    implicit none

    integer :: myrank, maxrank, ier
    include 'mpif.h'
#ifdef USE3D
    real :: angle, beta
    integer :: i,procs_per_plane, full_group, plane_group
    integer, allocatable :: ranks(:)
    real, allocatable :: xvals(:)
    integer :: nvals
#endif

    ! load mesh
    call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ier)

    ! Define boundary types and zone types
    if(imulti_region.eq.1) then
       if(boundary_type(1).eq.BOUND_UNKNOWN) then
          if(myrank.eq.0) &
               print *, "Boundary types not specified.  Using defaults."
          boundary_type(1) = BOUND_FIRSTWALL
          boundary_type(2) = BOUND_FIRSTWALL
          boundary_type(5) = BOUND_DOMAIN
          boundary_type(6) = BOUND_DOMAIN
       end if

       if(zone_type(1).eq.ZONE_UNKNOWN) then
          if(myrank.eq.0) &
               print *, "Zone types not specified.  Using defaults."
          zone_type(1) = ZONE_PLASMA
          zone_type(2) = ZONE_CONDUCTOR
          zone_type(3) = ZONE_VACUUM
       end if

    else
       boundary_type(1) = BOUND_ANY

       zone_type(1) = ZONE_PLASMA
    end if

#ifdef USEST
    ilog = -1  ! before reading in geometry
    if (igeometry.gt.0) then ! do nothing when igeometry==0 
       if (iread_vmec.ge.1) then ! read geometry from VMEC file
          call process_vmec(myrank)
          if (mod(nfp, nperiods).ne.0 .and. ifull_torus.eq.0) then
             if(myrank.eq.0) print *, 'Error: nfp must be divisible by nperiods'
             call safestop(7)
          end if 
       else ! read boudary geometry
          call read_boundary_geometry(myrank)
       end if
       ilog = 1       ! use logical basis funtions first
       call physical_mesh_setup(toroidal_period)
    end if
#endif
    if(ifull_torus.eq.0 .and. nperiods.ne.1) then
       toroidal_period = toroidal_period/nperiods 
       if(myrank.eq.0) print *, '1/', nperiods, ' of a torus...'
    else if(ifull_torus.eq.1 .or. nperiods.eq.1) then
       if(myrank.eq.0) print *, 'full torus...'
    end if

    if (iadapt_snap.eq.1) call m3dc1_model_settopo()

#ifdef USE3D   
    if(myrank.eq.0) print *, 'setting number of planes = ', nplanes
    call m3dc1_model_setnumplane(nplanes)

    procs_per_plane = maxrank/nplanes
    if(myrank.eq.0) &
         print *, 'number of processes per plane = ', procs_per_plane

    if(myrank.eq.0) print *, 'loading partitioned mesh...'
    if(is_rectilinear) then
       if(myrank.eq.0) print *, 'rectilinear mesh model...'
    else
       if(myrank.eq.0) print *, 'curved mesh model...'
    endif
    write(name_buff,"(A,A)") mesh_model(1:len_trim(mesh_model)),0
    call m3dc1_model_load(name_buff)
    write(name_buff,"(A,A)")  mesh_filename(1:len_trim(mesh_filename)),0
    call m3dc1_mesh_load (name_buff)

    ! set up toroidal angles
    if(iread_planes.eq.1) then
       nvals = 0
       call read_ascii_column('plane_positions', xvals, nvals, icol=1)
       if(nvals.ne.nplanes) then
          if(myrank.eq.0) &
               print *, 'Error: number of planes in plane_positions != nplanes'
          call safestop(8)
       end if
       do i=1, nvals
          if(xvals(i).lt.0. .or. xvals(i).ge.toroidal_period) then 
             if(myrank.eq.0) print *, 'Error: plane ', i, ' is at angle ', xvals(i), &
             ' which is outside the range [ 0, ', toroidal_period, ').'
             call safestop(8)
          end if
          call m3dc1_plane_setphi(i-1, xvals(i))
          if(myrank.eq.0) print *, 'Plane ', i, 'at angle ', xvals(i)
       end do
       deallocate(xvals)
    else
       do i=0, nplanes-1
          if(toroidal_pack_factor.gt.1. .and. i.gt.0) then
             beta = 2.*sqrt(alog(toroidal_pack_factor))
             angle = toroidal_pack_angle + &
                  (toroidal_period/2.)*(1. + &
                  erf(beta*(real(i)/real(nplanes) - 0.5)) / &
                  erf(beta/2.))
          else
             angle = toroidal_period*real(i)/real(nplanes)
          end if
          call m3dc1_plane_setphi(i, angle)
          if(myrank.eq.0) print *, 'Plane ', i+1, 'at angle ', angle
       end do
    end if


    ! build mesh
    if(myrank.eq.0) print *, 'setting up 3D mesh...'
    call m3dc1_mesh_build3d(0,0,0)

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
    ! there are two ways to use 2D mesh
    ! in serial, it loads the whole mesh
    ! in parallel, it loads N-part distributed mesh.
    ! to get N-part distributed mesh, run split_smb provided as mesh utilities
    write(name_buff,"(A,A)")  mesh_model(1:len_trim(mesh_model)),0
    call m3dc1_model_load(name_buff)
    write(name_buff,"(A,A)")  mesh_filename(1:len_trim(mesh_filename)),0
    call m3dc1_mesh_load (name_buff)
#endif
    call update_nodes_owned

    initialized = .true.
  end subroutine load_mesh

  subroutine update_nodes_owned
    integer :: numnodes, inode2, inode, ier, myrank, own_process
    include 'mpif.h'
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ier)
    if(allocated(nodes_owned)) deallocate(nodes_owned)
    allocate(nodes_owned(owned_nodes()))
    numnodes =  local_nodes()
    inode2=1;
    do inode=1, numnodes
       call m3dc1_ent_getownpartid(0, inode-1, own_process)
       if (own_process .eq. myrank) then
          nodes_owned(inode2) = inode
          inode2 = inode2 +1;
       end if
    end do
    if(inode2-1 .ne. owned_nodes()) call abort()

    call get_global_dims
  end subroutine update_nodes_owned
  subroutine unload_mesh
    if(.not. initialized) return
    deallocate(nodes_owned)
    call m3dc1_domain_finalize
  end subroutine unload_mesh

  subroutine print_node_data()
    implicit none
    integer :: i,numnodes
    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv(3), x, phi, z

    numnodes = local_nodes()

    do i=1, numnodes
       call boundary_node(i, is_boundary, izone, izonedim, normal, curv, &
            x, phi, z)
       if(.not.is_boundary) cycle
       write(*,'(7f12.4)') x, z, normal(1), normal(2), curv
    end do
  end subroutine print_node_data

  subroutine get_global_dims

    implicit none

    include 'mpif.h'

    integer :: nelms, error

    nelms = local_elements()

  ! Calculate offset of current process
    call mpi_scan(nelms, offset_elms, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, error)

    offset_elms = offset_elms - nelms
    ! print *, "[",myrank,"] hdf5_get_local_elms Offset ", offset

!  call numglobalents(global_nodes, gobal_edges, global_elms, global_regions)
!  print *, 'myrank, local_elms, global_elms, offset', &
!       myrank, nelms, global_elms, offset
    call mpi_allreduce(nelms, global_elms, 1, MPI_INTEGER, &
         MPI_SUM, MPI_COMM_WORLD, error)

    elms_per_plane = global_elms / nplanes

  end subroutine get_global_dims


  !======================================================================
  ! local_plane
  ! ~~~~~~~~~~~
  ! returns toroidal plane associated with current process
  !======================================================================
  integer function local_plane()
    implicit none
    call m3dc1_model_getplaneid(local_plane)
  end function local_plane

  !==============================================================
  ! local_nodes
  ! ~~~~~~~~~~~
  ! returns the number of elements local to this process
  !==============================================================
  integer function local_elements()
    implicit none
    integer :: elem_dim
    elem_dim = 2
#ifdef USE3D
    elem_dim = 3
#endif
    call m3dc1_mesh_getnument(elem_dim, local_elements)
  end function local_elements

  !==============================================================
  ! owned_nodes
  ! ~~~~~~~~~~~
  ! returns the number of nodes owned by to this process
  !==============================================================
  integer function owned_nodes()
    implicit none
    call m3dc1_mesh_getnumownent (0, owned_nodes)
  end function owned_nodes


  !==============================================================
  ! local_nodes
  ! ~~~~~~~~~~~
  ! returns the number of nodes local to this process
  !==============================================================
  integer function local_nodes()
    implicit none
    call m3dc1_mesh_getnument(0, local_nodes)
  end function local_nodes


  !==============================================================
  ! get_element_nodes
  ! ~~~~~~~~~~~~~~~~~
  ! returns the indices of each element node
  !==============================================================
  subroutine get_element_nodes(ielm,n)
    implicit none
    integer, intent(in) :: ielm
    integer, intent(out), dimension(nodes_per_element) :: n
    integer :: elem_dim, numelm, nodes_per_element_get
    elem_dim = 2
#ifdef USE3D
    elem_dim = 3
#endif
#ifdef DEBUG
    call m3dc1_mesh_getnument(elem_dim, numelm)
    if (ielm-1 .ge. numelm) then
       print *, "ielm ", ielm, " should be equal to or less than ", numelm
    endif
#endif
    call m3dc1_ent_getadj (elem_dim, ielm-1, 0, n, nodes_per_element, nodes_per_element_get)
    if (nodes_per_element_get .ne. nodes_per_element) then
      print *, "error get_element_nodes: nodes_per_element_get",nodes_per_element_get
      call abort 
    end if
    n(:)=n(:)+1
  end subroutine get_element_nodes

  !==============================================================
  subroutine populate_adjacency_matrix()
    implicit none

    integer :: nelms, gid_sz
    integer :: i, j, k
    integer, allocatable :: gids(:), pids(:), local_id(:)
    integer :: adj_sz

    nelms = local_elements()

    gid_sz = max_adj*nelms

    allocate(n_adjacent(nelms))
    allocate(adjacent(max_adj, nelms))

    allocate(local_id(nelms))
    allocate(pids(gid_sz))
    allocate(gids(gid_sz))

    do i=1, nelms
       local_id(i) = i-1
    end do

    n_adjacent = 0
#ifdef USE3D
    call m3dc1_ent_getglobaladj(3, local_id, nelms, 3, &
         n_adjacent, pids, gids, gid_sz, adj_sz)
#else
    call m3dc1_ent_getglobaladj(2, local_id, nelms, 2, &
         n_adjacent, pids, gids, gid_sz, adj_sz)
#endif

    k = 1
    do i=1, nelms
       do j=1, max_adj
          if(j.le.n_adjacent(i)) then
             adjacent(j,i) = gids(k)
             k = k+1
          else
             adjacent(j,i) = -1
          endif
       end do
!!$#ifdef USE3D
!!$       ! Manually add toroidal adjacencies
!!$       if(n_adjacent(i)+2 .gt. max_adj) then
!!$          print *, 'Error: adding toroidal adjacencies results in n_adjacent > max_adj'
!!$       else
!!$          j = i-1-elms_per_plane
!!$          if(j .lt. 0) j = j + global_elms
!!$          adjacent(n_adjacent(i)+1,i) = j
!!$
!!$          j = i-1+elms_per_plane
!!$          if(j .ge. global_elms) j = j - global_elms
!!$          adjacent(n_adjacent(i)+2,i) = j
!!$       end if
!!$#endif
    end do

    deallocate(local_id, pids, gids)
  end subroutine populate_adjacency_matrix

  subroutine clear_adjacency_matrix
    implicit none

    deallocate(n_adjacent)
    deallocate(adjacent)
  end subroutine clear_adjacency_matrix

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
       call m3dc1_model_getmincoord(bb(1), bb(2))
       call m3dc1_model_getmaxcoord(bb(3), bb(4))
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
    
    call m3dc1_node_getcoord(inode-1,coords)
#ifdef USEST
    if (igeometry.gt.1) then ! Use original mesh when igeometry==0 or 1 
       call physical_geometry(x, z, coords(1), coords(3), coords(2))
    else
#endif
    x = coords(1)
    z = coords(2)
#ifdef USEST
    endif
#endif

    if(is_rectilinear) then
       call m3dc1_model_getmincoord(x1, z1)
       x = x + xzero - x1
       z = z + zzero - z1
    endif

#ifdef USE3D
    phi = coords(3)
#else
    phi = 0.
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
       call abort()
    endif

#ifdef USE3D
    call get_node_pos(nodeids(4), x2, phi2, z2)
    d%d = phi2 - d%Phi
    if(d%d .le. 0.) d%d = d%d + toroidal_period
!    if(itri.eq.1) print *, 'd%phi = ', d%phi
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

!!  !============================================================
!!  ! whattri2
!!  subroutine whattri2(x,phi,z,ielm,xref,zref,nophi)
!!    implicit none
!!    
!!    real, intent(in) :: x, phi, z
!!    integer, intent(inout) :: ielm
!!    real, intent(out) :: xref, zref
!!    logical, intent(in), optional :: nophi
!!    
!!    integer :: nelms, i
!!    type(element_data) :: d
!! 
!!    ! first, try ielm
!!    if(ielm.gt.0) then
!!       call get_element_data(ielm,d)
!!       if(is_in_element_notol(ielm,d,x,phi,z,nophi)) then
!!          xref = d%R
!!          zref = d%Z
!!          return
!!       endif
!!    endif
!! 
!!    ! try all elements
!!    ielm = -1
!!    nelms = local_elements()
!!    do i=1, nelms
!!       call get_element_data(i,d)
!!       if(is_in_element_notol(i,d,x,phi,z,nophi)) then
!!          xref = d%R
!!          zref = d%Z
!!          ielm = i
!!          return
!!       endif
!!    end do
!! 
!!  end subroutine whattri2
!! 
!!  !! This is almost same with is_in_element.
!!  !! However, tol in is_in_element is relavtive tolerance.
!!  !! Therefore, overlapping region cannot be calculated with local information
!!  logical function is_in_element_notol(itri,d, R, phi, z, nophi)
!!     implicit none
!!     integer, intent(in) :: itri
!!     type(element_data), intent(in) :: d
!!     real, intent(in) :: R, phi, Z
!!     logical, intent(in), optional :: nophi
!! 
!!     real :: f, xi, zi, eta
!!     real, parameter :: zero = 0e0
!!     integer :: edge3, edge1, edge2, edge_sum
!!     integer :: ent_dim, loc_indx
!! #ifdef USE3D    
!!     logical :: np
!! #endif    
!! 
!!     call global_to_local(d, R, Phi, Z, xi, zi, eta)
!! 
!!     is_in_element_notol = .false.
!!     edge1 = 0
!!     edge2 = 0
!!     edge3 = 0 
!!     if(eta.lt.zero) then
!!        return
!!     elseif(eta .lt. epsilon(zero)) then
!!        edge3 = 1
!!     endif
!!     if(eta.gt.d%c) return
!! 
!!     f = 1. - eta/d%c
!!     if(xi.lt.-f*d%b) then
!!        return
!!     elseif(xi .lt. -f*d%b+epsilon(zero)) then
!!        edge2 = 1
!!     endif   
!!        
!!     if(xi.gt. f*d%a) then
!!        return
!!     elseif (xi .gt. f*d%a-epsilon(zero)) then
!!        edge1 = 1
!!     endif   
!!     
!! #ifdef USE3D
!!     if(present(nophi)) then
!!        np=nophi 
!!     else 
!!        np=.false.
!!     endif
!! 
!!     !! i.e. [0, d%d) : mine
!!     if(.not.np) then 
!!        if(zi.lt.zero) return
!!        if(zi.ge.d%d) return
!!     end if
!! #endif
!! 
!!     edge_sum = edge1+edge2+edge3
!!     if(edge_sum .ne. 0) then
!!        !! need tie treatment, edge_sum should be 1(edge) or 2(node)
!!        ent_dim = 2-edge_sum
!!        if(edge1 .eq. 1 ) then
!!           if(edge2 .eq. 1 ) then
!!              loc_indx=3  !!node
!!           elseif(edge3 .eq. 1) then
!!              loc_indx=2  !!node
!!           else
!!              loc_indx=1  !!edge
!!           endif
!!        else
!!           if(edge2 .eq. 1) then
!!              if(edge3 .eq. 1) then
!!                 loc_indx=1 !!node
!!              else
!!                 loc_indx=2 !!edge
!!              endif
!!           else
!!              loc_indx=3 !!edge
!!           endif
!!        endif
!!        if(is_shared_ent_others(itri,ent_dim,loc_indx)) return
!!     endif   
!! 
!!     is_in_element_notol = .true.
!!  end function is_in_element_notol
!! 
!!  logical function is_shared_ent_others(itri,ent_dim,loc_indx)
!!     implicit none
!!     integer, intent(in) :: itri, ent_dim, loc_indx
!! 
!!     integer :: nodeids(nodes_per_element)
!!     integer :: ent_indx, mine, num_adj_ent, node1, node2, i
!!     integer :: iedge(edges_per_element)
!!     integer :: inodes(2)
!! 
!!     if(ent_dim .eq. 0) then
!!        call get_element_nodes(itri, nodeids)
!!        ent_indx=nodeids(loc_indx)
!!     else
!!        call m3dc1_ent_getadj (2, itri-1, 1, iedge, edges_per_element, num_adj_ent)
!!        node1 = mod(loc_indx+1,3)-1
!!        node2 = mod(loc_indx+2,3)-1
!!        do i=1,edges_per_element
!!           call m3dc1_ent_getadj (1, itri-1, 0, inodes, nodes_per_element, num_adj_ent)
!!           if(node1 .eq. inodes(1) .and. node2 .eq. inodes(2) .or. &
!!              node2 .eq. inodes(1) .and. node1 .eq. inodes(2) ) then
!!              ent_indx=iedge(i)+1 
!!           endif   
!!        enddo   
!!     endif
!!     !! check if the given entity(dim, indx) is shared
!!     !! if shared(actually, should be shared), check remote ranks with mine
!!     !! if mine, return false
!!     !! if NOT mine, return true
!!     call m3dc1_ent_ismine(ent_dim, ent_indx, mine)
!!     if(mine == 1) then
!!        is_shared_ent_others = .false.
!!     else
!!        is_shared_ent_others = .true.
!!     endif   
!!  end function is_shared_ent_others

  !=========================================
  ! is_boundary_node
  ! ~~~~~~~~~~~~~~~~
  ! returns true if node lies on boundary
  !=========================================
  logical function is_boundary_node(inode, tags)
    implicit none
    integer, intent(in) :: inode
    integer, intent(in), optional :: tags
    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv(3), x, phi, z

    call boundary_node(inode,is_boundary,izone,izonedim,normal,curv,&
         x,phi,z,tags)
    is_boundary_node = is_boundary
  end function is_boundary_node

  !======================================================================
  ! boundary_node
  ! ~~~~~~~~~~~~~
  ! determines if node is on boundary, and returns relevant info
  ! about boundary surface
  !======================================================================
  subroutine boundary_node(inode,is_boundary,izone,izonedim,normal,curv,&
       x,phi,z,tags)

    use math
    
    implicit none
    
    integer, intent(in) :: inode              ! node index
    integer, intent(out) :: izone,izonedim    ! zone type/dimension
    real, intent(out) :: normal(2), curv(3)
    real, intent(out) :: x,phi,z              ! coordinates of inode
    logical, intent(out) :: is_boundary       ! is inode on boundary
    integer, intent(in), optional :: tags
    
    integer :: ibottom, iright, ileft, itop, ib
    real :: angler
    real, dimension(3) :: norm

    !curv = 0.

    call m3dc1_ent_getgeomclass(0,inode-1,izonedim,izone)
    call get_node_pos(inode,x,phi,z)

    if(is_rectilinear) then
       if(izonedim.ge.2) then
          is_boundary = .false.
          return
       end if
       
       !call m3dc1_model_getedge(ileft, iright, ibottom, itop)
       !for rect domain, always use the following order of model edges
       ileft=4
       iright=2
       ibottom=1
       itop=3
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
       call m3dc1_node_isongeombdry(inode-1,ib)
       is_boundary = ib.eq.1
       if(.not.is_boundary) return
#ifdef USEST  
       if (igeometry.gt.0 .and. ilog.ne.1) then ! do nothing when igeometry==0 
          call get_boundary_curv(normal,curv,inode)
       else
#endif
          call m3dc1_node_getnormvec(inode-1, norm)
          normal = norm(1:2)
          if(icurv.eq.0) then
             curv = 0.
          else
             call m3dc1_node_getcurv(inode-1, curv(1))
          end if
#ifdef USEST  
       end if
#endif 
       if(imulti_region.gt.0) then
          is_boundary = .false.
          if(boundary_type(izone).eq.BOUND_ANY) then
             is_boundary = .true.
          else
             if(present(tags)) then
                if(tags.eq.BOUND_ANY) then
                   is_boundary = .true.
                else if(tags.eq.boundary_type(izone)) then
                   is_boundary = .true.
                end if
             else if(boundary_type(izone).eq.BOUND_DOMAIN) then
                is_boundary = .true.
             end if
          end if
       end if
    end if
  end subroutine boundary_node

  subroutine boundary_edge(itrin, is_edge, normal, idim, tags)

    implicit none
    integer, intent(in) :: itrin
    integer, intent(out) :: is_edge(3)
    real, intent(out) :: normal(2,3)
    integer, intent(out) :: idim(3)
    
    integer :: inode(nodes_per_element), izone, i, j, jj
    integer :: inode2(nodes_per_element)
    integer :: in(2)
    real :: x, phi, z, c(3)
    logical :: is_bound(3), found_edge

    integer :: iedge(3), iedge2(3), izonedim, itri, ifaczone,ifaczonedim
    integer :: num_adj_ent, numelm
    integer, intent(in), optional :: tags

#ifdef USE3D
    ! the first face must be on the original plane
    call m3dc1_region_getoriginalface (itrin-1, itri)
    itri = itri+1
#else
    itri = itrin
#endif
    !call nodfac(itri,inode)
    !call edgfac(itri,iedge)
    !call zonfac(itri,ifaczone,ifaczonedim)
    call m3dc1_ent_getadj (2, itri-1, 0, inode, 3, num_adj_ent)
    inode = inode+1
    inode2 = inode
    call m3dc1_ent_getadj (2, itri-1, 1, iedge, 3, num_adj_ent)
    iedge = iedge+1
    iedge2 = iedge
    call m3dc1_ent_getgeomclass(2, itri-1, ifaczonedim, ifaczone)
    do i=1,3
       if(present(tags)) then 
          call boundary_node(inode(i),is_bound(i),izone,idim(i), &
               normal(:,i),c(i),x,phi,z,tags)
       else
          call boundary_node(inode(i),is_bound(i),izone,idim(i), &
               normal(:,i),c(i),x,phi,z,BOUND_ANY)
       end if
    end do

    is_edge = 0
    do i=1,3
       !call zonedg(iedge(i),izone,izonedim)
       call m3dc1_ent_getgeomclass(1, iedge2(i)-1, izonedim, izone)
       if(izonedim.gt.1) cycle

       !call nodedg(iedge(i), in)
       call m3dc1_ent_getadj (1, iedge2(i)-1, 0, in, 2, num_adj_ent)
       in = in + 1
       ! check nodes to see which one is the first node of iedge(i)
       found_edge = .false.
       do j=1, 3
          jj = mod(j,3)+1
          if(in(1).eq.inode2(j) .and. in(2).eq.inode2(jj)) then
             found_edge = .true.
             exit
          else if(in(2).eq.inode2(j) .and. in(1).eq.inode2(jj)) then
!             print *, 'Error: clockwise edge!'
             found_edge = .true.
             exit
          end if
       end do
       if(found_edge) then
          if(is_edge(j).ne.0) print *, 'Warning: edge counted twice'
          is_edge(j) = izone
!          call boundary_node(inode(j),is_bound(j),izone,idim(j), &
!               normal(:,j),c(j),x,phi,z)
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

  subroutine get_zone_index(itri, izone)
    integer, intent(in) :: itri
    integer, intent(out) :: izone
    integer :: izonedim
    integer :: elem_dim
    elem_dim = 2
#ifdef USE3D
    elem_dim = 3
#endif
    call m3dc1_ent_getgeomclass(elem_dim, itri-1,izonedim,izone)
  end subroutine get_zone_index

  subroutine get_zone(itri, izone)
    integer, intent(in) :: itri
    integer, intent(out) :: izone
    integer :: izone_ind

    call get_zone_index(itri, izone_ind)
    izone = zone_type(izone_ind)
  end subroutine get_zone



  !======================================================================
  ! local_dof_vector
  ! ~~~~~~~~~~~~~~~~~~
  ! calculates the coefficients of the polynomial expansion of the
  ! field in the element domain
  !======================================================================
  subroutine local_dof_vector(itri, c)
    implicit none

    integer, intent(in) :: itri
    real, intent(out), dimension(dofs_per_element,coeffs_per_element) :: c

    type(element_data) :: d
    real, dimension(coeffs_per_tri,coeffs_per_tri) :: t
    real, dimension(coeffs_per_dphi,coeffs_per_dphi) :: h
    integer :: i, j, k, l, m, n
    integer :: idof, icoeff, ip, it

    call get_element_data(itri,d)
    call tmatrix(t,d%a,d%b,d%c)
    call hmatrix(h,d%d)

    c = 0.

    icoeff = 0
    do i=1,coeffs_per_dphi
       do j=1,coeffs_per_tri
          icoeff = icoeff + 1
          idof = 0
          do k=1,tor_nodes_per_element
             do l=1,pol_nodes_per_element
                do m=1,tor_dofs_per_node
                   do n=1,pol_dofs_per_node
                      idof = idof + 1
                      
                      ip = n + (l-1)*pol_dofs_per_node
                      it = m + (k-1)*tor_dofs_per_node
                      c(idof,icoeff) = c(idof,icoeff) &
                           + h(it,i)*t(ip,j)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine local_dof_vector

  subroutine local_dofs(itri, dof, c)
    implicit none

    type(element_data) :: d
    integer, intent(in) :: itri
    vectype, intent(out), dimension(dofs_per_element) :: dof
    vectype, intent(in), dimension(coeffs_per_element) :: c

    real, dimension(dofs_per_element,coeffs_per_element) :: cl
    integer :: i, j, k
    real :: norm(2), curv(3)
    vectype, dimension(dofs_per_element) :: temp

#ifdef USEST
    integer :: inode(nodes_per_element)
    real, dimension(dofs_per_element, dofs_per_element) :: l2p_mat 
    real, dimension(dofs_per_node, dofs_per_node) :: p2l_mat 

    integer :: info1, info2
    real :: wkspce(9400)
    integer :: ipiv(dofs_per_node)
#endif

    call local_dof_vector(itri, cl)

    temp = 0
    do j=1, coeffs_per_element
       temp(:) = temp(:) + cl(:,j)*c(j)
    end do

    ! calculate the rotation matrix rot
    call get_element_data(itri,d)
    norm(1) = d%co
    norm(2) = d%sn
    curv = 0.

#ifdef USEST
    if(igeometry.eq.1.and.ilog.eq.0) then
       call get_element_nodes(itri, inode)
    end if
    l2p_mat = 0.
#endif
    do i=1, nodes_per_element
       j = (i-1)*dofs_per_node+1
       k = j + dofs_per_node - 1
       call rotate_dofs(temp(j:k), dof(j:k), norm, curv, -1)
#ifdef USEST
       if(igeometry.eq.1.and.ilog.eq.0) then
          info1 = 0
          info2 = 0
          call p2l_matrix(p2l_mat, inode(i)) 
          call dgetrf(dofs_per_node,dofs_per_node,p2l_mat,dofs_per_node,ipiv,info1)
          call dgetri(dofs_per_node,p2l_mat,dofs_per_node,ipiv,wkspce,400,info2)
          if(info1.ne.0.or.info2.ne.0) write(*,'(3I5)') info1,info2
          l2p_mat(j:k,j:k) = p2l_mat
          !call l2p_matrix(l2p_mat(j:k,j:k),inode(i))
       end if
#endif
    end do
#ifdef USEST
    if(igeometry.eq.1.and.ilog.eq.0) then
       dof = matmul(l2p_mat,dof)
    end if
#endif
  end subroutine local_dofs

end module scorec_mesh_mod

module adapt
  use vector_mod
  real :: ke_previous
  real :: target_error
  real :: error_tol = 4.
  real , dimension(2) :: abs_size, rel_size
  data rel_size /0.3,1.5/
  !type(vector_type), private :: error_vec 
  contains
  subroutine adapt_by_psi
    use basic
    use mesh_mod
    use arrays
    use newvar_mod
    use sparse
    use hdf5_output
    use diagnostics
    use boundary_conditions
    use time_step
    use m3dc1_output
    use auxiliary_fields
    use pellet
    use scorec_mesh_mod

    vectype, dimension(dofs_per_node) :: dat
    integer :: izone, izonedim, inode(nodes_per_element)
    integer :: numelms, itri, num_adj_ent, ier
    character(len=32) :: mesh_file_name

    call create_field(temporary_field)
    if(eqsubtract.eq.1) then
       temporary_field = psi_field(0)
       call add(temporary_field, psi_field(1))
    else
       temporary_field = psi_field(1)
    end if
    if(icsubtract.eq.1) call add(temporary_field, psi_coil_field)
    if(imulti_region.eq.1 .and. adapt_psin_vacuum.ne.0.) then
       dat = 0.
       dat(1) = (psibound - psimin)*adapt_psin_vacuum + psimin
       numelms = local_elements()
       do itri=1, numelms
          !call zonfac(itri,izone,izonedim)
          call m3dc1_ent_getgeomclass(2, itri-1, izonedim, izone)
          if(izone.ge.3) then
             !call nodfac(itri,inode)
             call m3dc1_ent_getadj (2, itri-1, 0, inode, 3, num_adj_ent)
             inode = inode+1
             do i=1,3
                call set_node_data(temporary_field,inode(i),dat)
             end do
         end if
       end do
       call sum_shared(temporary_field%vec)
    end if

    call straighten_fields()
    write(mesh_file_name,"(A7,I0,A)"),'meshOrg', ntime,0
    call m3dc1_mesh_write (mesh_file_name, 0)
    call adapt_by_field(temporary_field%vec%id,psimin,psibound)
    write(mesh_file_name,"(A7,A)"),'adapted', 0
    call m3dc1_mesh_write (mesh_file_name,0)
    call m3dc1_mesh_write (mesh_file_name,1)

    call destroy_field(temporary_field)
    call space(0)
    call update_nodes_owned()
    call tridef
    call unstraighten_fields()

    call create_newvar_matrices
    if(irestart .ne. 0) return
    print *, "re-calculate equlibrium after adapt .."
    call initial_conditions
    ! combine the equilibrium and perturbed fields of linear=0
    ! unless eqsubtract = 1
    if(eqsubtract.eq.0) then
       call add(field_vec, field0_vec)
       field0_vec = 0.
    endif
    call derived_quantities(1)
    ke_previous = ekin
  end subroutine adapt_by_psi
#ifdef USEADAPTBYERROR   
  subroutine adapt_by_error
    use diagnostics
    use basic
    use error_estimate
    use scorec_mesh_mod
    use basic
    use mesh_mod
    use arrays
    use newvar_mod
    use sparse
    use time_step
    use auxiliary_fields
    use scorec_mesh_mod
!#include "mpif.h"
    vectype, allocatable :: edge_error(:,:)
    vectype, allocatable :: elm_error(:,:), elm_error_res(:), elm_error_res_psi(:), elm_error_sum(:)
    real, allocatable :: node_error(:)
    integer :: num_edge, ii, jj,kk, num_get, num_get2, num_elm, num_node
    integer, dimension(3) :: nodes
    integer, dimension(256) :: elms 
    type(element_data) :: d

    integer, parameter :: vecsize=1 
    integer :: ndofs, num_node_total
    vectype, dimension(vecsize*dofs_per_node) :: dofData
    real :: facarea, facarea_tol, max_error, buff
    character(len=32) :: mesh_file_name, file_name1, file_name2, file_name3

    vectype, dimension(dofs_per_node*num_fields) :: max_val, min_val
    vectype :: maxPhi, maxPs

    write(mesh_file_name,"(A5,I0,A)"),'adapt', ntime,0 
    write(file_name1, "(A9,I0,A)"),'errorJump', ntime,0
    write(file_name2, "(A8,I0,A)"),'errorElm', ntime,0
    write(file_name3,"(A8,I0,A)"),'errorSum', ntime,0

    !call m3dc1_field_max(jphi_field%vec%id, max_val, min_val)
    call m3dc1_mesh_getnumglobalent (0, num_node_total)

    !if(myrank .eq. 0) print*, "time", ntime, "current max", max_val(1), min_val(1), "mesh size before adapt", num_node_total
    call m3dc1_field_max(field_vec%id, max_val,min_val)
    maxPhi = max(abs(max_val(1+(u_g-1)*dofs_per_node)),abs(min_val(1+(u_g-1)*dofs_per_node)))
    maxPs =  max(abs(max_val((psi_g-1)*dofs_per_node)+1),abs(min_val((psi_g-1)*dofs_per_node)+1))

    call m3dc1_mesh_getnument(1, num_edge)
    allocate(edge_error(num_edge, NUMTERM))
    call m3dc1_mesh_getnument(2, num_elm)
    allocate(elm_error(num_elm, NUMTERM))
    allocate(elm_error_sum(num_elm))
    allocate(elm_error_res(num_elm))
    allocate(elm_error_res_psi(num_elm))
    elm_error_res = 0
    elm_error_res_psi = 0
    call  m3dc1_mesh_getnument(0, num_node)
    allocate(node_error(num_node))
    node_error=0.
    edge_error=0.
    elm_error=0.
    call jump_discontinuity(edge_error)
    do ii=1, num_edge
       call m3dc1_ent_getadj (1, ii-1, 2, elms, 2, num_get)
       do jj=1, num_get 
          elm_error(elms(jj)+1,:)=elm_error(elms(jj)+1,:)+0.5*edge_error(ii,:);
       end do
    end do

    if(itor .eq. 0 .and. numvar .eq. 1 .and. linear .eq. 0 .and. ivform .eq. 0) call elem_residule (elm_error_res, elm_error_res_psi)

    elm_error_sum (:) = elm_error(:,NUMTERM) + elm_error_res (:) + elm_error_res_psi(:)
    !if(isplitstep .eq. 1) elm_error_sum (:) = elm_error_sum (:) + elm_error(:,TSQ) 
    call output_face_data (NUMTERM, sqrt(real(elm_error)), file_name1);
    call output_face_data (1, sqrt(real(elm_error_res)), file_name2);
    call output_face_data (1, sqrt(real(elm_error_sum)), file_name3);

    call get_node_error_from_elm (real(elm_error_sum), 1, node_error);
       
    node_error = sqrt(node_error)
    !print *, edge_error
    deallocate(edge_error)
    deallocate(elm_error)
    deallocate(elm_error_sum)
    deallocate(elm_error_res)

    max_error= maxval(node_error(:))
    buff = max_error
    call mpi_allreduce (buff, max_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier )
    if (myrank .eq. 0) print *, "max extimated error ", max_error
    if (iadapt .eq. 2 .and. max_error .gt. error_tol * target_error) then
       if (myrank .eq. 0) print *, " error exceeds tolerance, start adapting mesh"
       call straighten_fields()
       abs_size(1) = adapt_hmin
       abs_size(2) = adapt_hmax
       call set_mesh_size_bound (abs_size, rel_size)
       call adapt_by_error_field(node_error, target_error);
       call m3dc1_mesh_write (mesh_file_name,0)
       call space(0)
       call update_nodes_owned()
       call tridef
       call unstraighten_fields()

       call create_newvar_matrices
       call derived_quantities(1)
       meshAdapted =1
    end if
    deallocate(node_error)
  end subroutine adapt_by_error
#endif
  integer function run_adapt ()
    use diagnostics
    use basic
    !print *, "check run_adaptA: ke_previous, ekin, ntime",ke_previous,ekin,ntime
    run_adapt = 0
    !if(mod(ntime,1) .eq. 0) then
     ! run_adapt=1
      !ke_previous = ekin;
    !end if
    if(linear .eq. 0 .and. mod(ntime,1) .eq. 0) run_adapt=1
    !if(linear .eq. 1 .and. mod(ntime,50) .eq. 0) run_adapt=1
  end function

  subroutine straighten_fields ()
    use diagnostics
    use basic
    use error_estimate
    use scorec_mesh_mod
    use basic
    use mesh_mod
    use arrays
    use newvar_mod
    use sparse
    use time_step
    use auxiliary_fields
    use scorec_mesh_mod

    integer :: ifield
    
    do ifield=0,1
       call straighten_field(u_field(ifield))
       call straighten_field(vz_field(ifield))
       call straighten_field(chi_field(ifield))

       call straighten_field(psi_field(ifield))
       call straighten_field(bz_field(ifield))
       call straighten_field(pe_field(ifield))

       call straighten_field(den_field(ifield))
       call straighten_field(p_field(ifield))
       call straighten_field(te_field(ifield))
       call straighten_field(ti_field(ifield))
    end do
  end subroutine straighten_fields 

  subroutine unstraighten_fields ()
    use diagnostics
    use basic
    use error_estimate
    use scorec_mesh_mod
    use basic
    use mesh_mod
    use arrays
    use newvar_mod
    use sparse
    use time_step
    use auxiliary_fields
    use scorec_mesh_mod

    integer :: ifield

    do ifield=0,1
       call unstraighten_field(u_field(ifield))
       call unstraighten_field(vz_field(ifield))
       call unstraighten_field(chi_field(ifield))

       call unstraighten_field(psi_field(ifield))
       call unstraighten_field(bz_field(ifield))
       call unstraighten_field(pe_field(ifield))

       call unstraighten_field(den_field(ifield))
       call unstraighten_field(p_field(ifield))
       call unstraighten_field(te_field(ifield))
       call unstraighten_field(ti_field(ifield))
    end do
  end subroutine unstraighten_fields
end module adapt


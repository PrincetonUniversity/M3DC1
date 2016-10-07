subroutine wrrestart
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet

  implicit none

#ifdef USESCOREC

  integer :: mmnn18, j1, numnodes, numelms
  integer :: ndofs
  integer, save :: ifirstrs = 1
  character (len=30) :: fname
  vectype, allocatable :: data_buff(:)

  call createfilename(fname)
  numnodes = local_nodes()
  numelms = local_elements()
  mmnn18 = 0
  
  ifirstrs = 0

 if (myrank .eq. 0) &
      write(*,*) '[P',myrank,'] write file ',fname

  open(56,file=fname,form='unformatted',status='replace', action='write')
  ! first put in information to check on run information
  write(56) numnodes
  write(56) numelms
  call m3dc1_field_getnumlocaldof(num_fields, ndofs)
  write(56) ndofs
  call m3dc1_field_getnumlocaldof(1,ndofs)
  write(56) ndofs
!  write(56) mmnn18
!  write(56) numvar
  write(56) iper
  write(56) jper 
  write(56) myrank
  write(56) maxrank
  write(56) eqsubtract
  write(56) linear
  write(56) icomplex

  call m3dc1_field_getnumlocaldof(num_fields, ndofs)
  allocate (data_buff(ndofs)) 
  call m3dc1_field_retrieve (field_vec%id, data_buff, ndofs) 
  do j1=1,ndofs 
     write(56) data_buff(j1)
  enddo
  call m3dc1_field_retrieve (field0_vec%id, data_buff, ndofs)
  do j1=1,ndofs 
     write(56) data_buff(j1)
  enddo

  write(56) ntime,time,dt
  write(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  write(56) psimin,psilim,psibound
  write(56) xnull, znull
  write(56) xmag, zmag
  deallocate (data_buff)
  call m3dc1_field_getnumlocaldof(1, ndofs)
  allocate (data_buff(ndofs))
  call m3dc1_field_retrieve (bf_field(1)%vec%id, data_buff, ndofs)
  do j1=1,ndofs 
     write(56) data_buff(j1)
  enddo
  call m3dc1_field_retrieve (bf_field(0)%vec%id, data_buff, ndofs)
  do j1=1,ndofs 
     write(56) data_buff(j1)
  enddo
  write(56) pellet_x, pellet_phi, pellet_z, &
       pellet_velx, pellet_velphi, pellet_velz, pellet_var
  write(56) version
  write(56) icsubtract
  if(icsubtract.eq.1) then
     call m3dc1_field_retrieve (psi_coil_field%vec%id, data_buff, ndofs)
     do j1=1,ndofs 
        write(56) data_buff(j1)
     enddo
  end if
  write(56) extsubtract, use_external_fields
  if(use_external_fields) then
     call m3dc1_field_retrieve(psi_ext%vec%id, data_buff, ndofs)
     do j1=1,ndofs 
        write(56) data_buff(j1)
     enddo
     call m3dc1_field_retrieve(bz_ext%vec%id, data_buff, ndofs)
     do j1=1,ndofs 
        write(56) data_buff(j1)
     enddo
     call m3dc1_field_retrieve(bf_ext%vec%id, data_buff, ndofs)
     do j1=1,ndofs 
        write(56) data_buff(j1)
     enddo
  end if
  write(56) pellet_rate
  write(56) xnull2, znull2
  write(56) r_p, r_p2

  deallocate(data_buff)
  close(56)

#endif

end subroutine wrrestart
!============================================================
subroutine rdrestart
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet

  implicit none

#ifdef USESCOREC

  integer :: j1, numnodes, inumnodes
  integer :: inumelms, immnn18, inumvar, iiper, ijper, imyrank
  integer :: imaxrank, numelms, ieqsubtract, ilinear, icomp
  character (len=30) :: fname
  integer :: ndofs
  integer :: iversion
  real :: vloopsave, pelletratesave
  vectype, allocatable :: data_buff(:)

  ! check if 2D-to-3D case
  open(76,file='C1restart00000',form='unformatted',status='unknown')
    read(76) inumnodes
    read(76) inumelms
    read(76) immnn18
    read(76) inumvar
    read(76) iiper
    read(76) ijper
    read(76) imyrank
    read(76) imaxrank
  close(76)

if (imaxrank .ne. maxrank) then
    if (myrank .eq. 0) &
      print *, '[M3D-C1 INFO] 3D Simulation with restart files: #ranks - 2D ', imaxrank
    call rdrestart_2d23d
else
   if (myrank .eq. 0 .and. nplanes .eq. 1) &
     print *, '[M3D-C1 INFO] 2D Simulation with restart files: #ranks - ',maxrank
   if (myrank .eq. 0 .and. nplanes .ne. 1) &
     print *, '[M3D-C1 INFO] 3D Simulation with restart files: #ranks - ',maxrank

  call createfilename(fname)
  call m3dc1_field_getnumlocaldof(num_fields, ndofs)
  numnodes = local_nodes()
  numelms = local_elements()

  open(56,file=fname,form='unformatted',status='unknown')
  read(56) inumnodes
  read(56) inumelms
  read(56) immnn18
  read(56) inumvar
  read(56) iiper
  read(56) ijper
  read(56) imyrank
  read(56) imaxrank
  read(56) ieqsubtract
  read(56) ilinear
  read(56) icomp

  if(inumnodes .ne. numnodes .or. inumelms .ne. numelms .or. &
       iiper .ne. iper .or. ijper .ne. jper .or. &
       imyrank .ne. myrank .or. imaxrank .ne. maxrank) then
     write(*,*) 'Restart file information does not match!'
     close(56)
     if(inumnodes .ne. numnodes) then
        write(*,*) 'numnodes ',inumnodes, numnodes, myrank
     endif
     if(inumelms .ne. numelms) then
        write(*,*) 'numelms ',inumnodes, numnodes, myrank
     endif
     if(iiper .ne. iper) then
        write(*,*) 'iper',iiper, iper, myrank
     endif
     if(ijper .ne. jper) then
        write(*,*) 'jper',ijper, jper, myrank
     endif
     if(imyrank .ne. myrank) then
        write(*,*) 'myrank',imyrank,myrank
     endif
     if(imaxrank .ne. maxrank) then
        write(*,*) 'maxrank',imaxrank, maxrank, myrank
     endif
     call safestop(2)
  endif

  allocate (data_buff(ndofs))

  do j1=1,ndofs
     read(56) data_buff(j1)
  enddo
  call m3dc1_field_set(field_vec%id, data_buff, ndofs)
  do j1=1,ndofs
     read(56) data_buff(j1)
  enddo
  call m3dc1_field_set(field0_vec%id, data_buff,ndofs)

  ! If we are running a linear simulation, but the restart file was
  ! a nonlinear simulation, make the restart data be the equilibrium
  if(linear.eq.1 .and. ilinear.eq.0) then
     call m3dc1_field_add(field0_vec%id, field_vec%id)
     call m3dc1_field_assign (field_vec%id, 0., 0)
  endif

  vloopsave = vloop
  read(56) ntime,time,dt
  read(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  read(56,END=1199) psimin,psilim,psibound
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag
  if(control_type .eq. -1) vloop = vloopsave  ! use vloop from input if no control on I

  deallocate (data_buff)
  call m3dc1_field_getnumlocaldof(1, ndofs)
  allocate (data_buff(ndofs))
  do j1=1,ndofs
     read(56,END=1199) data_buff(j1)
  enddo
  call m3dc1_field_set(bf_field(1)%vec%id, data_buff,ndofs)
  do j1=1,ndofs
     read(56,END=1199) data_buff(j1)
  enddo
  call m3dc1_field_set(bf_field(0)%vec%id, data_buff, ndofs)

  read(56, END=1199) pellet_x, pellet_phi, pellet_z, &
       pellet_velx, pellet_velphi, pellet_velz, pellet_var

  read(56, END=1199) iversion

  if(iversion.ge.7) then
     read(56, END=1199) icsubtract
     if(icsubtract.eq.1) then
        do j1=1,ndofs
           read(56,END=1199) data_buff(j1)
        enddo
        call m3dc1_field_set(psi_coil_field%vec%id, data_buff, ndofs)
     end if
  end if

  if(iversion.ge.10) then
     read(56, END=1199) extsubtract, use_external_fields
     if(use_external_fields) then
        call create_field(psi_ext)
        call create_field(bz_ext)
        call create_field(bf_ext)
        do j1=1,ndofs
           read(56,END=1199) data_buff(j1)
        enddo
        call m3dc1_field_set(psi_ext%vec%id, data_buff, ndofs)
        do j1=1,ndofs
           read(56,END=1199) data_buff(j1)
        enddo
        call m3dc1_field_set(bz_ext%vec%id, data_buff, ndofs)
        do j1=1,ndofs
           read(56,END=1199) data_buff(j1)
        enddo
        call m3dc1_field_set(bf_ext%vec%id, data_buff, ndofs)
     end if
  end if

  if(iversion.ge.12) then 
     pelletratesave = pellet_rate
     read(56,END=1199) pellet_rate
     ! use pellet_rate from input if no control
     if(n_control%icontrol_type .eq. -1) pellet_rate = pelletratesave  
  end if

  if(iversion.ge.13) then
     read(56,END=1199) xnull2, znull2
  end if
  if(iversion.ge.14) then
     read(56,END=1199) r_p, r_p2
  end if

  deallocate (data_buff)
  goto 1200
1199 if(myrank.eq.0) &
        print *, '[M3D-C1 ERROR] failed reading restart file'
1200 close(56)
endif
#endif
end subroutine rdrestart

!============================================================
subroutine rdrestart_2d23d
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet

  implicit none
  integer :: i, j, numnodes, prev_numnodes, iversion
  integer :: prev_numelms, prev_mmnn18, prev_numvar, prev_iper, prev_jper, prev_myrank
  integer :: prev_maxrank, numelms, prev_eqsubtract, prev_linear, prev_comp
  character (len=30) :: fname
  integer :: prev_ndofs1, ndofs1, prev_ndofs2, ndofs2, group_rank
  integer :: prev_ndofs1_pernode, prev_ndofs2_pernode, cur_ndofs1_pernode, cur_ndofs2_pernode 
  real :: vloopsave, pelletratesave
 real, allocatable :: data_buf(:) 
  real, dimension(num_fields*12*2):: dofs_node ! buffer for dofs per node

  call get2dfilename(fname)
  group_rank = modulo(myrank, maxrank/nplanes)
  numnodes = local_nodes()
  numelms = local_elements()

  open(56,file=fname,form='unformatted',status='unknown')
  read(56) prev_numnodes
  read(56) prev_numelms
  read(56) prev_ndofs1
  read(56) prev_ndofs2
  read(56) prev_iper
  read(56) prev_jper 
  read(56) prev_myrank
  read(56) prev_maxrank
  read(56) prev_eqsubtract
  read(56) prev_linear
  read(56) prev_comp
  
  if ((prev_numnodes*2 .ne. numnodes .and. prev_numnodes .ne. numnodes) .or. &
       prev_iper .ne. iper .or. prev_jper .ne. jper) then
     write(*,*) 'restart file information does not match!'
     close(56)
     if (prev_numnodes*2 .ne. numnodes .and. prev_numnodes .ne. numnodes) then
        write(*,*) 'numnodes: prev ',prev_numnodes, ', cur ', numnodes, myrank 
     endif
     if (prev_iper .ne. iper) then
        write(*,*) 'iper: prev ',prev_iper, ', cur ',iper, myrank
     endif
     if (prev_jper .ne. jper) then
        write(*,*) 'jper: prev ',prev_jper, ', cur ',jper, myrank
     endif
     call safestop(2)
  endif

! Allocate space for the arrays tmp_
  call m3dc1_field_getnumlocaldof(num_fields,ndofs1)
  call m3dc1_field_getnumlocaldof(1,ndofs2)

    prev_ndofs1_pernode = prev_ndofs1/prev_numnodes
    cur_ndofs1_pernode = ndofs1/numnodes
    allocate(data_buf(prev_ndofs1))
    !fill field_vec
    do i=1,prev_ndofs1
      read(56) data_buf(i)
    enddo
    call m3dc1_field_assign(field_vec%id, 0., 0)
    do i=1,prev_numnodes
       dofs_node =0.
       do j=0, prev_ndofs1_pernode/6-1
          dofs_node(1+j*12:j*12+6) &
            = data_buf((i-1)*prev_ndofs1_pernode+1+j*6:(i-1)*prev_ndofs1_pernode+j*6+6)
       end do
       call m3dc1_ent_setdofdata (0, i-1, field_vec%id, cur_ndofs1_pernode, dofs_node(1:cur_ndofs1_pernode))
    enddo
    !call m3dc1_field_printcompnorm(field_vec%id, "field_vec before sync "//char(0))
    call m3dc1_field_sync (field_vec%id)
    !call m3dc1_field_printcompnorm(field_vec%id, "field_vec after sync "//char(0))
    !fill field0_vec 
    do i=1,prev_ndofs1
      read(56) data_buf(i)
    enddo 
    call m3dc1_field_assign(field0_vec%id, 0., 0)
    do i=1,prev_numnodes
       dofs_node =0.
       do j=0, prev_ndofs1_pernode/6-1
           dofs_node(1+j*12:j*12+6)  &
            = data_buf((i-1)*prev_ndofs1_pernode+1+j*6:(i-1)*prev_ndofs1_pernode+j*6+6)
       end do
       call m3dc1_ent_setdofdata (0, i-1, field0_vec%id, cur_ndofs1_pernode, dofs_node(1:cur_ndofs1_pernode))
    enddo
    !call m3dc1_field_printcompnorm(field0_vec%id, "field0_vec before sync "//char(0))
    call m3dc1_field_sync (field0_vec%id)
    !call m3dc1_field_printcompnorm(field0_vec%id, "field0_vec after sync "//char(0)) 

  deallocate(data_buf)

  ! If we are running a linear simulation, but the restart file was
  ! a nonlinear simulation, make the restart data be the equilibrium
  if(linear.eq.1 .and. prev_linear.eq.0) then
     call m3dc1_field_add(field0_vec%id, field_vec%id)
     call m3dc1_field_assign (field_vec%id, 0., 0)
  endif

  vloopsave = vloop
  read(56) ntime,time,dt
  read(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  read(56,END=1199) psimin,psilim,psibound
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag
  if (control_type .eq. -1) &
      vloop = vloopsave  ! use vloop from input if no control on I

  if (myrank .eq. 0) then 
      write(*,*) '[M3DC1 INFO] Setting 3D fields from 2D restart files and initializing timestep (ntime=0)'
      write(*,*) '      #ranks: prev ', prev_maxrank, ', cur ', maxrank
      write(*,*) '      #nodes: prev ', prev_numnodes, ', cur ', numnodes
      write(*,*) '	#dofs1: prev ',prev_ndofs1, ', cur ',ndofs1
      write(*,*) '      #dofs1: prev ',prev_ndofs2, ', cur ',ndofs2
  endif
  ntime=0

    prev_ndofs2_pernode = prev_ndofs2/prev_numnodes
    cur_ndofs2_pernode = ndofs2/numnodes
    allocate(data_buf(prev_ndofs2))
    do i=1,prev_ndofs2
      read(56,END=1199) data_buf(i)
    enddo
    do i=1,prev_numnodes
      dofs_node=0.
      dofs_node(1:prev_ndofs2_pernode)&
         = data_buf((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
      call m3dc1_ent_setdofdata (0, i-1, bf_field(1)%vec%id, cur_ndofs2_pernode, dofs_node(1:cur_ndofs2_pernode))
    enddo
    call m3dc1_field_sync (bf_field(1)%vec%id)
    do i=1,prev_ndofs2
      read(56,END=1199) data_buf(i)
    enddo
    do i=1,prev_numnodes
      dofs_node=0.
      dofs_node(1:prev_ndofs2_pernode)&
         = data_buf((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
      call m3dc1_ent_setdofdata (0, i-1, bf_field(0)%vec%id, cur_ndofs2_pernode, dofs_node(1:cur_ndofs2_pernode))
    enddo
    call m3dc1_field_sync (bf_field(0)%vec%id)

  read(56, END=1199) pellet_x, pellet_phi, pellet_z, &
       pellet_velx, pellet_velphi, pellet_velz, pellet_var

  read(56, END=1199) iversion

  if(iversion.ge.7) then
     read(56, END=1199) icsubtract
     if(icsubtract.eq.1) then
          do i=1,prev_ndofs2
            read(56,END=1199) data_buf(i)
          enddo
          do i=1,prev_numnodes
            dofs_node =0.
            dofs_node(1:prev_ndofs2_pernode) &
              = data_buf((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
            call m3dc1_ent_setdofdata (0, i-1, psi_coil_field%vec%id, cur_ndofs2_pernode, dofs_node(1:cur_ndofs2_pernode))
          end do
          call m3dc1_field_printcompnorm(psi_coil_field%vec%id, "psi_coil_field before sync "//char(0))
          call m3dc1_field_sync (psi_coil_field%vec%id)
          call m3dc1_field_printcompnorm(psi_coil_field%vec%id, "psi_coil_field after sync "//char(0))
     end if
  end if

  if(version.ge.10) then
     read(56, END=1199) extsubtract, use_external_fields
     if(use_external_fields) then
        call create_field(psi_ext)
        call create_field(bz_ext)
        call create_field(bf_ext)
          do i=1,prev_ndofs2
            read(56,END=1199) data_buf(i)
          enddo
          do i=1,prev_numnodes
            dofs_node =0.
            dofs_node(1:prev_ndofs2_pernode) &
                = data_buf((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
            call m3dc1_ent_setdofdata (0, i-1, psi_ext%vec%id, cur_ndofs2_pernode, dofs_node(1:cur_ndofs2_pernode))
          enddo
          call m3dc1_field_sync (psi_ext%vec%id)
          do i=1,prev_ndofs2
            read(56,END=1199) data_buf(i)
          enddo
          do i=1,prev_numnodes
            dofs_node =0.
            dofs_node(1:prev_ndofs2_pernode) &
                = data_buf((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
             call m3dc1_ent_setdofdata (0, i-1, bz_ext%vec%id, cur_ndofs2_pernode, dofs_node(1:cur_ndofs2_pernode))
          enddo
          call m3dc1_field_sync (bz_ext%vec%id)
          do i=1,prev_ndofs2
            read(56,END=1199) data_buf(i)
          enddo
          do i=1,prev_numnodes
            dofs_node =0.
            dofs_node(1:prev_ndofs2_pernode) &
                = data_buf((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
             call m3dc1_ent_setdofdata (0, i-1, bf_ext%vec%id, cur_ndofs2_pernode, dofs_node(1:cur_ndofs2_pernode))
          enddo
          call m3dc1_field_sync (bf_ext%vec%id)
     end if
  end if

  if(iversion.ge.12) then 
     pelletratesave = pellet_rate
     read(56,END=1199) pellet_rate
     ! use pellet_rate from input if no control
     if(n_control%icontrol_type .eq. -1) pellet_rate = pelletratesave  
  end if

  if(iversion.ge.13) then
     read(56,END=1199) xnull2, znull2
  end if
  if(iversion.ge.14) then
     read(56,END=1199) r_p, r_p2
  end if

  deallocate(data_buf)
  goto 1200
1199 if (myrank.eq.0) &
          print *, '[M3D-C1 ERROR] failed reading restart file'
1200 close(56)
end subroutine rdrestart_2d23d

!============================================================
subroutine rdrestart_cplx
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet
  use init_common

  implicit none
  
#ifdef USESCOREC

  integer :: j1, numnodes, inumnodes
  integer :: inumelms, immnn18, inumvar, iiper, ijper, imyrank
  integer :: imaxrank, numelms, ieqsubtract, ilinear, icomp
  character (len=30) :: fname
  integer :: ndofs
  integer :: iversion
  real :: vloopsave, pelletratesave
  vectype, allocatable :: data_buff(:)
  real :: tmprestart


  call createfilename(fname)
  call m3dc1_field_getnumlocaldof(num_fields, ndofs)
  numnodes = local_nodes()
  numelms = local_elements()

  open(56,file=fname,form='unformatted',status='unknown')
  read(56) inumnodes
  read(56) inumelms
  read(56) immnn18
  read(56) inumvar
  read(56) iiper
  read(56) ijper 
  read(56) imyrank
  read(56) imaxrank
  read(56) ieqsubtract
  read(56) ilinear
  read(56) icomp

  if(inumnodes .ne. numnodes .or. inumelms .ne. numelms .or. &
       iiper .ne. iper .or. ijper .ne. jper .or. &
       imyrank .ne. myrank .or. imaxrank .ne. maxrank) then
     write(*,*) 'Restart file information does not match!'
     close(56)
     if(inumnodes .ne. numnodes) then
        write(*,*) 'numnodes ',inumnodes, numnodes, myrank 
     endif
     if(inumelms .ne. numelms) then
        write(*,*) 'numelms ',inumnodes, numnodes, myrank 
     endif
     if(iiper .ne. iper) then
        write(*,*) 'iper',iiper, iper, myrank
     endif
     if(ijper .ne. jper) then
        write(*,*) 'jper',ijper, jper, myrank
     endif
     if(imyrank .ne. myrank) then
        write(*,*) 'myrank',imyrank,myrank
     endif
     if(imaxrank .ne. maxrank) then
        write(*,*) 'maxrank',imaxrank, maxrank, myrank
     endif
     call safestop(2)
  endif

  allocate(data_buff(ndofs))

  do j1=1,ndofs 
     read(56) tmprestart
#ifdef USECOMPLEX
     data_buff(j1)=cmplx(tmprestart, 0.)
#else
     data_buff(j1)=tmprestart
#endif
  enddo
  call m3dc1_field_set(field_vec%id, data_buff, ndofs)
  do j1=1,ndofs 
     read(56) tmprestart
#ifdef USECOMPLEX
     data_buff(j1)=cmplx(tmprestart, 0.)
#else
     data_buff(j1)=tmprestart
#endif
  enddo
  call m3dc1_field_set(field0_vec%id, data_buff,ndofs)

                         
  ! If we are running a linear simulation, but the restart file was
  ! a nonlinear simulation, make the restart data be the equilibrium
  if(linear.eq.1 .and. ilinear.eq.0) then
     call m3dc1_field_add(field0_vec%id, field_vec%id)
     call m3dc1_field_assign (field_vec%id, 0., 0)
  endif

  vloopsave = vloop
  read(56) ntime,time,dt
  read(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  read(56,END=1199) psimin,psilim,psibound
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag
  if(control_type .eq. -1) vloop = vloopsave  ! use vloop from input if no control on I

  deallocate (data_buff)
  call m3dc1_field_getnumlocaldof(1, ndofs)
  allocate (data_buff(ndofs))
  do j1=1,ndofs 
     read(56,END=1199) tmprestart
#ifdef USECOMPLEX
     data_buff(j1)=cmplx(tmprestart, 0.)
#else
     data_buff(j1)=tmprestart
#endif
  enddo
  call m3dc1_field_set(bf_field(1)%vec%id, data_buff,ndofs)
  do j1=1,ndofs 
     read(56,END=1199) tmprestart
#ifdef USECOMPLEX
     data_buff(j1)=cmplx(tmprestart, 0.)
#else
     data_buff(j1)=tmprestart
#endif
  enddo
  call m3dc1_field_set(bf_field(0)%vec%id, data_buff, ndofs)

  read(56, END=1199) pellet_x, pellet_phi, pellet_z, &
       pellet_velx, pellet_velphi, pellet_velz, pellet_var

  read(56, END=1199) iversion

  if(iversion.ge.7) then
     read(56, END=1199) icsubtract
     if(icsubtract.eq.1) then
        do j1=1,ndofs 
           read(56,END=1199) tmprestart
#ifdef USECOMPLEX
           data_buff(j1)=cmplx(tmprestart, 0.)
#else
           data_buff(j1)=tmprestart
#endif
        enddo
        call m3dc1_field_set(psi_coil_field%vec%id, data_buff, ndofs)
     end if
  end if

  if(iversion.ge.10) then
     read(56, END=1199) extsubtract, use_external_fields
     if(use_external_fields) then
        call create_field(psi_ext)
        call create_field(bz_ext)
        call create_field(bf_ext)
        do j1=1,ndofs
!          read(56,END=1199) data_buff(j1)
           read(56,END=1199) tmprestart
#ifdef USECOMPLEX
           data_buff(j1)=cmplx(tmprestart, 0.)
#else
           data_buff(j1)=tmprestart
#endif
        enddo
        call m3dc1_field_set(psi_ext%vec%id, data_buff, ndofs)
        do j1=1,ndofs
!          read(56,END=1199) data_buff(j1)
           read(56,END=1199) tmprestart
#ifdef USECOMPLEX
           data_buff(j1)=cmplx(tmprestart, 0.)
#else
           data_buff(j1)=tmprestart
#endif
        enddo
        call m3dc1_field_set(bz_ext%vec%id, data_buff, ndofs)
        do j1=1,ndofs
!          read(56,END=1199) data_buff(j1)
           read(56,END=1199) tmprestart
#ifdef USECOMPLEX
           data_buff(j1)=cmplx(tmprestart, 0.)
#else
           data_buff(j1)=tmprestart
#endif
        enddo
        call m3dc1_field_set(bf_ext%vec%id, data_buff, ndofs)
     end if
  end if

  if(iversion.ge.12) then 
     pelletratesave = pellet_rate
     read(56,END=1199) pellet_rate
     ! use pellet_rate from input if no control
     if(n_control%icontrol_type .eq. -1) pellet_rate = pelletratesave  
  end if

  if(iversion.ge.13) then
     read(56,END=1199) xnull2, znull2
  end if
  if(iversion.ge.14) then
     read(56,END=1199) r_p, r_p2
  end if


  deallocate(data_buff)
  goto 1200
1199 if(myrank.eq.0) &
print *, '[M3D-C1 ERROR] failed reading restart file'
1200 close(56)
  ntime = 0
  time = 0
  call init_perturbations

#endif
end subroutine rdrestart_cplx

!============================================================
subroutine createfilename(filename)
  implicit none

  include 'mpif.h'
  character (len=30) :: filename
  character (len=5) :: charprocnum
  integer :: myrank, j, ier, i

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_rank:',ier
     call safestop(1)
  endif
                                ! initialize the SUPERLU process grid
  do j=1,5
     i = (myrank / 10**(5-j)) - (myrank / 10**(6-j)) * 10
     select case(i)
     case(0)
        charprocnum(j:j) = '0'
     case(1)
        charprocnum(j:j) = '1'
     case(2)
        charprocnum(j:j) = '2'
     case(3)
        charprocnum(j:j) = '3'
     case(4)
        charprocnum(j:j) = '4'
     case(5)
        charprocnum(j:j) = '5'
     case(6)
        charprocnum(j:j) = '6'
     case(7)
        charprocnum(j:j) = '7'
     case(8)
        charprocnum(j:j) = '8'
     case(9)
        charprocnum(j:j) = '9'
     end select
  enddo

  filename = 'C1restart'//charprocnum
  
  return
end subroutine createfilename

!============================================================
subroutine get2dfilename(filename)
  use basic

  implicit none

  include 'mpif.h'
  character (len=30) :: filename
  character (len=5) :: charprocnum
  integer ::j, ier, i, group_rank

  group_rank = modulo(myrank, maxrank/nplanes)
                                ! initialize the SUPERLU process grid
  do j=1,5
     i = (group_rank / 10**(5-j)) - (group_rank / 10**(6-j)) * 10
     select case(i)
     case(0)
        charprocnum(j:j) = '0'
     case(1)
        charprocnum(j:j) = '1'
     case(2)
        charprocnum(j:j) = '2'
     case(3)
        charprocnum(j:j) = '3'
     case(4)
        charprocnum(j:j) = '4'
     case(5)
        charprocnum(j:j) = '5'
     case(6)
        charprocnum(j:j) = '6'
     case(7)
        charprocnum(j:j) = '7'
     case(8)
        charprocnum(j:j) = '8'
     case(9)
        charprocnum(j:j) = '9'
     end select
  enddo

  filename = 'C1restart'//charprocnum
  !write(*,*) '[P',myrank,'] group_rank: ', group_rank,' get2dfilename: ',filename
  return
end subroutine get2dfilename

!============================================================
subroutine wrrestartglobal
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  implicit none

#ifdef USESCOREC
  include 'mpif.h'
  integer :: ifirstrs2, numdofsglobal, mpitag
  integer :: ierr, numentsglobal(4), i, numnodes, j, globalid, ndofs
  integer :: begindofnumber, enddofnumberplusone
  integer :: recvinfo(1)
  integer :: status(MPI_STATUS_SIZE)
  vectype, allocatable :: data_buff(:) ,data_buff2(:)
  data ifirstrs2/1/

  if(ifirstrs2 .ne. 1 .and. myrank .eq. 0) then
  endif

  ifirstrs2 = 0
  mpitag = 9070
  recvinfo(1) = 1
  ierr = 0

  if(myrank .ne. 0) then
     ! unaligned access warnings on PPPL sgi machines 
     ! is from communication here
     call MPI_RECV(recvinfo, 1, MPI_INTEGER, myrank-1,  &
          mpitag, MPI_COMM_WORLD, status, ierr)
     if(ierr .ne. MPI_SUCCESS) then
        write(*,*) 'MPI recv problem in restart', ierr
        call safestop(828)
     endif

#ifdef _AIX
     print *, "must use wrrestart on ibm instead of wrrestartglobal"
     call safestop(9229)
#else
     open(56,file='C1restart',form='unformatted',status='old',&
          position='append')
#endif
     endfile 56
  else
     ifirstrs2 = 0
     open(56,file='C1restart',form='unformatted',status='replace', &
          action='write')
     ! first put in information to check on run information
     call m3dc1_field_getnumglobaldof(num_fields, numdofsglobal)
     do i=0,3
        call m3dc1_mesh_getnumglobalent (i, numentsglobal(i+1))
     end do
     write(56) numentsglobal(1)
     write(56) numentsglobal(2)
     write(56) numentsglobal(3)
     write(56) numdofsglobal
     write(56) numvar
     write(56) iper
     write(56) jper   
     write(56) ntime,time,dt
     write(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
     write(56) icomplex
  endif

  call m3dc1_field_getnumlocaldof(field_vec%id, ndofs)
  allocate (data_buff(ndofs))
  call m3dc1_field_retrieve (field_vec%id, data_buff, ndofs)

  call m3dc1_field_getnumlocaldof(field0_vec%id, ndofs)
  allocate (data_buff2(ndofs))
  call m3dc1_field_retrieve (field0_vec%id, data_buff2, ndofs)

  numnodes = local_nodes()
  do j=1,numnodes
     call m3dc1_node_getglobalid(0, j, globalid)
     write(56) globalid
     call m3dc1_ent_getlocaldofid(0, j, num_fields, j, begindofnumber, &
          enddofnumberplusone)
     do i=begindofnumber,enddofnumberplusone-1
        write(56) data_buff(i)
     enddo
     do i=begindofnumber,enddofnumberplusone-1
        write(56) data_buff2(i)
     enddo
  enddo
  i = -1
  !  below specifies that this is the last line
  if(myrank .eq. maxrank-1) write(56) i

  if(myrank .ne. maxrank-1) then
     call MPI_SEND(recvinfo, 1, MPI_INTEGER, myrank+1, &
          mpitag, MPI_COMM_WORLD, ierr)
     if(ierr .ne. MPI_SUCCESS) then
        write(*,*) 'MPI send problem in restart'
        call safestop(827)
     endif
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if(myrank.eq.0) then
     write(56) psimin,psilim,psibound,xnull,znull
     write(56) xmag,zmag
  endif
  close(56)
  deallocate (data_buff)
  deallocate (data_buff2) 
1201 format("* * * restart file written at cycle,time = ",i5,1pe12.4,    &
          "  * * *")
#endif
end subroutine wrrestartglobal
!============================================================
subroutine rdrestartglobal
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  
  implicit none

#ifdef USESCOREC

  include 'mpif.h'

  integer :: numdofsglobal, numentsglobal(4), i, numnodes, globalid, ndofs
  integer :: in1, in2, in3
  integer :: icomp
  real :: vloopsave
  vectype, allocatable :: data_buff(:) ,data_buff2(:)
  open(56,file='C1restart',form='unformatted',status='old')
      
  ! first put in information to check on run information
  call m3dc1_field_getnumglobaldof(num_fields, numdofsglobal)
  do i=0,3
     call m3dc1_mesh_getnumglobalent (i, numentsglobal(i+1))
  end do

  read(56) in1!, in2, in3
  read(56) in2
  read(56) in3
  write(*,*) 'in',in1, in2, in3
  if(in1 .ne. numentsglobal(1) .or. in2 .ne. numentsglobal(2) .or. &
       in3 .ne. numentsglobal(3)) then
     write(*,*) 'Mesh does not match for rdrestartglobal', &
          in1, in2, in3, &
          numentsglobal(1), numentsglobal(2), numentsglobal(3)
     call safestop(567)
  endif
  vloopsave = vloop
  read(56) numdofsglobal
  read(56) numvar
  read(56) iper
  read(56) jper   
  read(56) ntime,time,dt
  read(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  read(56) icomp
  if(control_type .eq. -1) vloop = vloopsave  ! vloop from input if no I control
  write(*,*) 'in restartglobal time, dts are ', ntime, time, dt

  call m3dc1_field_getnumlocaldof(field_vec%id, ndofs)
  allocate (data_buff(ndofs))

  call m3dc1_field_getnumlocaldof(field0_vec%id, ndofs)
  allocate (data_buff2(ndofs))

  numnodes = local_nodes()
  read(56) globalid
  do while(globalid .ne. -1)
     print *, "does not support globalidnod"
     call abort();
!     if(nodeid .ne. -1) then ! this node exists on this proc
!     write(56) globalid
!     call m3dc1_ent_getlocaldofid(0, j, num_fields, j, begindofnumber, &
!          enddofnumberplusone)
!
!           call entdofs(num_fields, nodeid, 0, begindofnumber, &
!                enddofnumberplusone)
!           do i=begindofnumber,enddofnumberplusone-1
!              read(56) realtemp
!              field_vec%data(i) = realtemp
!           enddo
!           do i=begindofnumber,enddofnumberplusone-1
!              read(56) realtemp
!              field0_vec%data(i) = realtemp
!           enddo
!        else
!           call entdofs(num_fields, nodeid, 0, begindofnumber, &
!                enddofnumberplusone)
!          do i=begindofnumber,enddofnumberplusone-1
!              read(56) field_vec%data(i)
!           enddo
!           do i=begindofnumber,enddofnumberplusone-1
!              read(56) field0_vec%data(i)
!           enddo
!        endif
!     else ! node is not on this proc
!      if(icomp.eq.0 .and. icomplex.eq.1) then
!          ! the number of dof values at this node
!           do i=1,num_fields*6*3 
!              read(56) realtemp
!           enddo
!        else
           ! the number of dof values at this node
!           do i=1,num_fields*6*3 
!              read(56) tempvariable
!           enddo
!        endif
!     endif
!     read(56) globalid
  enddo
  
  read(56,END=1199) psimin,psilim,psibound
  write(*,*) 'psimin,psilim,psibound',psimin,psilim,psibound
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag
  goto 1200
1199 if(myrank.eq.0) &
print *, '[M3D-C1 ERROR] failed reading restart file'
1200 close(56)
  deallocate (data_buff)
  deallocate (data_buff2)
#endif
end subroutine rdrestartglobal



subroutine wrrestart_adios
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet

  implicit none

   include 'mpif.h'

!cj aug05-2011 #ifdef USESCOREC 
#ifdef USEADIOS
  integer :: mmnn18, j1, numnodes, numelms
  integer :: ndofs_1, ndofs_2, i, j
  integer, save :: ifirstrs = 1
  character (len=30) :: fname

  ! ADIOS variables declarations for matching gwrite_restart_c1.fh 
  integer             :: comm, ierr
  integer                 :: adios_err
  integer*8               :: adios_handle 
  real, allocatable :: tmp_field_vec(:), tmp_field0_vec(:)
  real, allocatable :: tmp_bf_field_1(:), tmp_bf_field_0(:)
  real, allocatable :: tmp_psi_ext(:), tmp_bz_ext(:), tmp_bf_ext(:)
  real, allocatable :: tmp_psi_coil_field(:)

  integer :: useext

  fname="restart.bp"
  numnodes = local_nodes()
  numelms = local_elements()


  !call numdofs(num_fields, ndofs_1)
  call m3dc1_field_getnumlocaldof(num_fields, ndofs_1)
  allocate(tmp_field_vec(ndofs_1))
  tmp_field_vec=0.
  allocate(tmp_field0_vec(ndofs_1))
  tmp_field0_vec=0.
  !tmp_field_vec(1:ndofs_1) = field_vec%data(1:ndofs_1)
  call m3dc1_field_retrieve (field_vec%id, tmp_field_vec, ndofs_1)
  !tmp_field0_vec(1:ndofs_1) = field0_vec%data(1:ndofs_1)
  call m3dc1_field_retrieve (field0_vec%id, tmp_field0_vec, ndofs_1)

  !call numdofs(1, ndofs_2)
  call m3dc1_field_getnumlocaldof(1, ndofs_2)
  allocate(tmp_bf_field_1(ndofs_2))
  tmp_bf_field_1=0.
  allocate(tmp_bf_field_0(ndofs_2))
  tmp_bf_field_0=0.
  !tmp_bf_field_1(1:ndofs_2)= bf_field(1)%vec%data(1:ndofs_2)
  call m3dc1_field_retrieve (bf_field(1)%vec%id, tmp_bf_field_1, ndofs_2)
  !tmp_bf_field_0(1:ndofs_2) = bf_field(0)%vec%data(1:ndofs_2)
  call m3dc1_field_retrieve (bf_field(0)%vec%id, tmp_bf_field_0, ndofs_2)

  allocate(tmp_psi_ext(ndofs_2))
  allocate(tmp_bz_ext(ndofs_2))
  allocate(tmp_bf_ext(ndofs_2))
  tmp_psi_ext=0
  tmp_bz_ext=0
  tmp_bf_ext=0
  if(use_external_fields) then
     !tmp_psi_ext(1:ndofs_2)= psi_ext%vec%data(1:ndofs_2)
     call m3dc1_field_retrieve (psi_ext%vec%id, tmp_psi_ext, ndofs_2)
     !tmp_bz_ext(1:ndofs_2) = bz_ext%vec%data(1:ndofs_2)
     call m3dc1_field_retrieve (bz_ext%vec%id, tmp_bz_ext, ndofs_2)
     !tmp_bf_ext(1:ndofs_2) = bf_ext%vec%data(1:ndofs_2)
     call m3dc1_field_retrieve (bf_ext%vec%id, tmp_bf_ext, ndofs_2)
     useext = 1
  else
     useext = 0
  end if

  allocate(tmp_psi_coil_field(ndofs_2))
  !ndofs_2=0
  if(icsubtract.eq.1) then
     call m3dc1_field_retrieve (psi_coil_field%vec%id,tmp_psi_coil_field, ndofs_2)
     call m3dc1_field_printcompnorm(psi_coil_field%vec%id, "psi_coil_field in wrrestart_adios"//char(0))
  end if

    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr) 
    if( ifirstrs .eq. 1 ) then
#ifdef USECOMPLEX
    call adios_init ("m3dc1_cplx.xml"//char(0), adios_err)
#else
    call adios_init ("m3dc1.xml"//char(0), adios_err)
#endif
    endif
    call adios_open (adios_handle, "restart", fname, "w", comm, adios_err)
#ifdef USECOMPLEX
#include "gwrite_restart_c1_cplx.fh" 
#else
#include "gwrite_restart_c1.fh" 
#endif
    call adios_close (adios_handle, adios_err)
    call MPI_Barrier (comm, ierr)
    if( (ntimemax-(ntime-ntime0)) .lt. ntimepr ) then
    call adios_finalize (myrank, adios_err)
    endif

  if(myrank.eq.0) &
      write(*,*) "OUTPUT: wrrestart_adios groupsize totalsize", &
                  adios_groupsize, adios_totalsize

  deallocate(tmp_field_vec, tmp_field0_vec, tmp_bf_field_1, tmp_bf_field_0, &
       tmp_psi_ext, tmp_bz_ext, tmp_bf_ext, tmp_psi_coil_field)

  ifirstrs = 0
#endif
end subroutine wrrestart_adios

!============================================================
subroutine rdrestart_adios
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet
  use gradshafranov
  implicit none
  
 include 'mpif.h'

#ifdef USEADIOS
  integer :: cur_nnodes, prev_nnodes
  integer :: prev_nelms, prev_mmnn18, prev_numvar, prev_iper, prev_jper, prev_myrank
  integer :: prev_maxrank, cur_nelms, prev_eqsubtract, prev_linear, prev_comp
  character (len=30) :: fname
  integer :: cur_ndofs1, cur_ndofs2, prev_ndofs1, prev_ndofs2
  real :: vloopsave

  ! ADIOS variables declarations for matching gread_restart_c11.fh &  gread_restart_c12.fh
  integer             :: comm, ierr, ier, itmp
  integer                 :: adios_err
  real, allocatable :: tmp_field_vec(:), tmp_field0_vec(:)
  real, allocatable :: tmp_bf_field_1(:), tmp_bf_field_0(:)
  real, allocatable :: tmp_psi_ext(:), tmp_bz_ext(:), tmp_bf_ext(:)
  real, allocatable :: tmp_psi_coil_field(:)

  !!! ADIOS variables for reading 
  integer*8                    :: fh, gh          ! file handler and group handler
  character(10), dimension(1)  :: gnamelist       ! expect 1 group only in restart file
  character(30), dimension(50) :: vnamelist    ! list of all vars 
  character(30), dimension(10) :: anamelist   ! list of all attributes 
  integer                      :: tstart, ntsteps ! timesteps available
  integer*8                    :: start, readsize, read_bytes ! bytes read by adios_read_var()
  integer                      :: i,j, gcnt, vcnt, acnt, ts, lastts
  integer                      :: vartype, ndim, timedim ! adios_inq_var()
  integer*8, dimension(2)      :: dims                   ! adios_inq_var()
  integer                      :: elemsize               ! double complex or double
  integer :: prev_useext, prev_version, group_size, group_rank
  integer :: prev_ndofs1_pernode, prev_ndofs2_pernode, cur_ndofs1_pernode, cur_ndofs2_pernode
  integer :: vec_created
 
  real, dimension(num_fields*12*2):: dofs_node1, dofs_node2, dofs_node3 ! buffer for dofs per node
  group_size = maxrank/nplanes
  group_rank = modulo(myrank, group_size)

  fname="restart.bp"
  cur_nnodes = local_nodes()
  cur_nelms = local_elements() 
  !call numdofs(num_fields, cur_ndofs1)
  call m3dc1_field_getnumlocaldof(num_fields, cur_ndofs1)
  !call numdofs(1, cur_ndofs2)
  call m3dc1_field_getnumlocaldof(1, cur_ndofs2)

  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  !! set read method, default=0 the BP reader, 1=BP subfiles with better performance
  call adios_set_read_method (0, adios_err);
  call adios_fopen (fh,fname,MPI_COMM_WORLD,gcnt,adios_err)
  if (adios_err .ne. 0) call safestop(2)
  call adios_inq_file (fh,vcnt,acnt,tstart,ntsteps,gnamelist,adios_err)

  !! Have to open a group from the file to read anything
  call adios_gopen (fh, gh, gnamelist(1), vcnt, acnt, ierr)
  ! not necessary to inquire if you know what to read
  call adios_inq_group(gh, vnamelist, anamelist, ts, lastts, ierr)

    ! Read in scalars (offset = 0, read size = 1 "in all dimensions")
    start = 0
    readsize = 1 ! number of elements of a type, not bytes!
    vloopsave = vloop
    pelletratesave = pellet_rate


    ! These scalars are varying over processors so we use adios_read_local_var() to 
    !  read the group_rank'th value of the scalar from the file
    call adios_read_local_var (gh, "numnodes",   group_rank, start, readsize, prev_nnodes, read_bytes)
    call adios_read_local_var (gh, "numelms",    group_rank, start, readsize, prev_nelms, read_bytes)
    if (prev_nnodes .eq. cur_nnodes .and. prev_nelms .eq. cur_nelms) then !2d to 2d or 3d to 3d
      group_rank = myrank ! turn the rank back to myrank
    endif
    call adios_read_local_var (gh, "mmnn18",     group_rank, start, readsize, prev_mmnn18, read_bytes)
    call adios_read_local_var (gh, "numvar",     group_rank, start, readsize, prev_numvar, read_bytes)
    call adios_read_local_var (gh, "iper",       group_rank, start, readsize, prev_iper, read_bytes)
    call adios_read_local_var (gh, "jper",       group_rank, start, readsize, prev_jper, read_bytes)
    call adios_read_local_var (gh, "myrank",     group_rank, start, readsize, prev_myrank, read_bytes)
    call adios_read_local_var (gh, "maxrank",    group_rank, start, readsize, prev_maxrank, read_bytes)
    call adios_read_local_var (gh, "eqsubtract", group_rank, start, readsize, prev_eqsubtract, read_bytes)
    call adios_read_local_var (gh, "useext",     group_rank, start, readsize, prev_useext, read_bytes)
    call adios_read_local_var (gh, "linear",     group_rank, start, readsize, prev_linear, read_bytes)
    call adios_read_local_var (gh, "icomplex",   group_rank, start, readsize, prev_comp, read_bytes)
    call adios_read_local_var (gh, "ntime",      group_rank, start, readsize, ntime, read_bytes)
    call adios_read_local_var (gh, "time",       group_rank, start, readsize, time, read_bytes)
    call adios_read_local_var (gh, "dt",         group_rank, start, readsize, dt, read_bytes)
    call adios_read_local_var (gh, "totcur0",    group_rank, start, readsize, totcur0, read_bytes)
    call adios_read_local_var (gh, "tflux0",     group_rank, start, readsize, tflux0, read_bytes)
    call adios_read_local_var (gh, "gbound",     group_rank, start, readsize, gbound, read_bytes)
    call adios_read_local_var (gh, "ptot",       group_rank, start, readsize, ptot, read_bytes)
    call adios_read_local_var (gh, "vloop",      group_rank, start, readsize, vloop, read_bytes)
    call adios_read_local_var (gh, "i_control%err_i",      group_rank, start, readsize, i_control%err_i, read_bytes)
    call adios_read_local_var (gh, "i_control%err_p_old",  group_rank, start, readsize, i_control%err_p_old, read_bytes)
    call adios_read_local_var (gh, "n_control%err_i",      group_rank, start, readsize, n_control%err_i, read_bytes)
    call adios_read_local_var (gh, "n_control%err_p_old",  group_rank, start, readsize, n_control%err_p_old, read_bytes)
    call adios_read_local_var (gh, "psimin",     group_rank, start, readsize, psimin, read_bytes)
    call adios_read_local_var (gh, "psilim",     group_rank, start, readsize, psilim, read_bytes)
    call adios_read_local_var (gh, "psibound",   group_rank, start, readsize, psibound, read_bytes)
    call adios_read_local_var (gh, "xnull",      group_rank, start, readsize, xnull, read_bytes)
    call adios_read_local_var (gh, "znull",      group_rank, start, readsize, znull, read_bytes)
    call adios_read_local_var (gh, "xmag",       group_rank, start, readsize, xmag, read_bytes)
    call adios_read_local_var (gh, "zmag",       group_rank, start, readsize, zmag, read_bytes)
    call adios_read_local_var (gh, "ndofs_1",    group_rank, start, readsize, prev_ndofs1, read_bytes)
    call adios_read_local_var (gh, "ndofs_2",    group_rank, start, readsize, prev_ndofs2, read_bytes)
    call adios_read_local_var (gh, "pellet_x",   group_rank, start, readsize, pellet_x, read_bytes)
    call adios_read_local_var (gh, "pellet_phi", group_rank, start, readsize, pellet_phi, read_bytes)
    call adios_read_local_var (gh, "pellet_z",   group_rank, start, readsize, pellet_z, read_bytes)
    call adios_read_local_var (gh, "pellet_velx",group_rank, start, readsize, pellet_velx, read_bytes)
    call adios_read_local_var (gh, "pellet_velphi",        group_rank, start, readsize, pellet_velphi, read_bytes)
    call adios_read_local_var (gh, "pellet_velz",group_rank, start, readsize, pellet_velz, read_bytes)
    call adios_read_local_var (gh, "pellet_var", group_rank, start, readsize, pellet_var, read_bytes)
    call adios_read_local_var (gh, "version",    group_rank, start, readsize, prev_version, read_bytes)
    if(prev_version.ge.12) then
       call adios_read_local_var (gh, "pellet_rate", group_rank, start, readsize, pellet_rate, read_bytes)
    end if
    if(prev_version.ge.13) then
       call adios_read_local_var (gh, "xnull2",      group_rank, start, readsize, xnull, read_bytes)
       call adios_read_local_var (gh, "znull2",      group_rank, start, readsize, znull, read_bytes)
    end if
    if(prev_version.ge.14) then
       call adios_read_local_var (gh, "r_p",      group_rank, start, readsize, xnull, read_bytes)
       call adios_read_local_var (gh, "r_p2",      group_rank, start, readsize, znull, read_bytes)
    end if


    if(control_type .eq. -1) vloop = vloopsave  !  vloop from input if no I control
    if(n_control%icontrol_type .eq. -1) pellet_rate = pelletratesave  !  vloop from input if no I control

  if (numvar .ne. prev_numvar .or. prev_myrank .ne. group_rank .or. &
      prev_iper .ne. iper .or. prev_jper .ne. jper) then
    write(*,*) 'Restart file information does not match! - program will stop'
    if (prev_numvar .ne. numvar) write(*,*) '(P ', myrank,') numvar: prev-',prev_numvar,', cur-', numvar
    if (prev_iper .ne. iper) write(*,*) '(P ', myrank,') iper: prev-',prev_iper,', cur-', iper
    if (prev_jper .ne. jper) write(*,*) '(P ', myrank,') jper: prev-',prev_jper,', cur-', jper
    if (prev_myrank .ne. group_rank) &
      write(*,*) '(P', myrank, ') rank in group: prev-',prev_myrank,', cur-', group_rank
    ierr = 1
  else
    ierr = 0
  end if

  prev_ndofs1_pernode = prev_ndofs1/prev_nnodes
  prev_ndofs2_pernode = prev_ndofs2/prev_nnodes
  cur_ndofs1_pernode = cur_ndofs1/cur_nnodes
  cur_ndofs2_pernode = cur_ndofs2/cur_nnodes

  if(prev_ndofs1_pernode .ne. cur_ndofs1_pernode) then
    if(myrank .eq. 0) print *, "Restarting 3d from 2d. Setting ntime = 0"
    ntime = 0
  end if
  if (prev_nnodes .ne. cur_nnodes .and. prev_nnodes*2 .ne. cur_nnodes) then
      write(*,*) '[M3DC1 ERROR] rdrestart_adios: #nodes mismatch. prev-',prev_nnodes,', cur-',cur_nnodes
      call safestop(2)
  else
    if (myrank .eq. 0) write(*,*) '[M3DC1 INFO] rdrestart_adios: maxrank: prev-',prev_maxrank, 'cur-', maxrank, &
                                  ', nplanes=', nplanes, ', group_size-', group_size, &
                                  ', #nodes: prev-', prev_nnodes,', cur-', cur_nnodes
  endif

  if (cur_ndofs1 .ne. prev_ndofs1 .and. myrank .eq. 0) &
    write(*,*) '[M3DC1 INFO] rdrestart_adios: ndofs1_pernode: prev-',prev_ndofs1_pernode, ', cur-',cur_ndofs1_pernode
  if (cur_ndofs2 .ne. prev_ndofs2 .and. myrank .eq. 0) &
    write(*,*) '[M3DC1 INFO] rdrestart_adios: ndofs2_pernode: prev-',prev_ndofs2_pernode, ', cur-',cur_ndofs2_pernode

  ! check for errors
  call mpi_allreduce(ierr,itmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  if(itmp.ne.0) call safestop(6)

    ! Allocate space for the arrays tmp_
  allocate(tmp_field_vec(prev_ndofs1))
  allocate(tmp_field0_vec(prev_ndofs1)) 
  allocate(tmp_bf_field_1(prev_ndofs2))
  allocate(tmp_bf_field_0(prev_ndofs2)) 
  allocate(tmp_psi_ext(prev_ndofs2))
  allocate(tmp_bz_ext(prev_ndofs2))
  allocate(tmp_bf_ext(prev_ndofs2))
  allocate(tmp_psi_coil_field(prev_ndofs2))

    ! Check types of the array variables
  call adios_inq_var (gh, "tmp_field_vec", vartype, ndim, dims, timedim, adios_err)
#ifdef USECOMPLEX 
  elemsize = 16 ! sizeof double complex
  if (vartype .ne. 11) then  ! 11 = double complex in adios
#else
  elemsize = 8 ! sizeof double
  if (vartype .ne. 6) then  ! 6 = double in adios
#endif
    write(*,*) '[M3DC1 ERROR] rdrestart_adios: scalar type mismatch (complex vs real). cur vartype-',vartype
    call safestop(2)
  endif

    ! Read in arrays organized as a global array: offset = (0, rank), read size = (*,1)
    start    = 0
    readsize = prev_ndofs1 ! number of elements of a type, not bytes!
    
    call adios_read_local_var (gh, "tmp_field_vec", group_rank, start, readsize, tmp_field_vec, read_bytes)

    if (read_bytes .ne. elemsize*prev_ndofs1) then
         write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading tmp_field_vec! read_bytes-', &
                    read_bytes, ', elemsize*prev_ndofs1-', elemsize*prev_ndofs1
         call safestop(2)
    endif
    
    call adios_read_local_var (gh, "tmp_field0_vec", group_rank, start, readsize, tmp_field0_vec, read_bytes)


    if (read_bytes .ne. elemsize*prev_ndofs1) then 
         write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading tmp_field0_vec! read_bytes-', &
                     read_bytes, ', elemsize*prev_ndofs1-', elemsize*prev_ndofs1
         call safestop(2)
    endif
    
    readsize = prev_ndofs2 ! number of elements of a type, not bytes!
    
    call adios_read_local_var (gh, "tmp_bf_field_1", group_rank, start, readsize, tmp_bf_field_1, read_bytes)

    if (read_bytes .ne. elemsize*prev_ndofs2) then  
         write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading tmp_bf_field_1! read_bytes-', &
                   read_bytes, ', elemsize*prev_ndofs2-', elemsize*prev_ndofs2
         call safestop(2)
    endif
    
    call adios_read_local_var (gh, "tmp_bf_field_0", group_rank, start, readsize, tmp_bf_field_0, read_bytes)

    if (read_bytes .ne. elemsize*prev_ndofs2) then 
         write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading tmp_bf_field_0! read_bytes-', &
                     read_bytes, ', elemsize*prev_ndofs2-', elemsize*prev_ndofs2
         call safestop(2)
    endif

    call adios_read_local_var (gh, "tmp_psi_ext", group_rank, start, readsize, tmp_psi_ext, read_bytes)

    if (read_bytes .ne. elemsize*prev_ndofs2) then 
         write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading tmp_psi_ext! read_bytes-', &
                     read_bytes, ', elemsize*prev_ndofs2-', elemsize*prev_ndofs2
         call safestop(2)
    endif

    call adios_read_local_var (gh, "tmp_bz_ext", group_rank, start, readsize, tmp_bz_ext, read_bytes)

    if (read_bytes .ne. elemsize*prev_ndofs2) then 
         write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading tmp_bz_ext! read_bytes-', &
                     read_bytes, ', elemsize*prev_ndofs2-', elemsize*prev_ndofs2
         call safestop(2)
    endif

    call adios_read_local_var (gh, "tmp_bf_ext", group_rank, start, readsize, tmp_bf_ext, read_bytes)

    if (read_bytes .ne. elemsize*prev_ndofs2) then 
         write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading tmp_bf_ext! read_bytes-', &
                     read_bytes, ', elemsize*prev_ndofs2-', elemsize*prev_ndofs2
         call safestop(2)
    endif
   
    if (icsubtract.eq.1) then
      if(myrank .eq. 0)  print*, "read psi_coil_field data"
      call adios_read_local_var (gh, "tmp_psi_coil_field", group_rank, start, readsize, tmp_psi_coil_field, read_bytes)

      if (read_bytes .ne. elemsize*prev_ndofs2) then
           write(*,*) '[M3DC1 ERROR] rdrestart_adios: Size mismatch at reading psi_coil_field! read_bytes-', &
                       read_bytes, ', elemsize*prev_ndofs2-', elemsize*prev_ndofs2
           call safestop(2)
      endif
    end if
 
    !! Close group and file
    call MPI_Barrier (comm, ierr)
    call adios_gclose(gh, ierr)
    call adios_fclose(fh, ierr)

    if (myrank.eq.0) &
      write(*,*) '[M3DC1 INFO] rdrestart_adios: reading ', fname

  if (prev_ndofs1==cur_ndofs1) then
    !field_vec%data(1:prev_ndofs1) = tmp_field_vec(1:prev_ndofs1)
    call m3dc1_field_set(field_vec%id, tmp_field_vec, prev_ndofs1)
    !field0_vec%data(1:prev_ndofs1) = tmp_field0_vec(1:prev_ndofs1)
    call m3dc1_field_set(field0_vec%id, tmp_field0_vec, prev_ndofs1)
  else
    !field_vec%data=0
    call m3dc1_field_assign(field_vec%id, 0., 0)
    !field0_vec%data=0
    call m3dc1_field_assign(field0_vec%id, 0., 0)
    do i=1,prev_nnodes
       dofs_node1 =0.
       dofs_node2 =0.
       do j=0, prev_ndofs1_pernode/6-1
          !field_vec
          dofs_node1(1+j*12:j*12+6) & 
            = tmp_field_vec((i-1)*prev_ndofs1_pernode+1+j*6:(i-1)*prev_ndofs1_pernode+j*6+6)
          !field0_vec
          dofs_node2(1+j*12:j*12+6)  & 
            = tmp_field0_vec((i-1)*prev_ndofs1_pernode+1+j*6:(i-1)*prev_ndofs1_pernode+j*6+6)
       end do
       call m3dc1_ent_setdofdata (0, i-1, field_vec%id, cur_ndofs1_pernode, dofs_node1(1:cur_ndofs1_pernode))
       call m3dc1_ent_setdofdata (0, i-1, field0_vec%id, cur_ndofs1_pernode, dofs_node2(1:cur_ndofs1_pernode))
    enddo
    call m3dc1_field_printcompnorm(field_vec%id, "field_vec before sync "//char(0))
    call m3dc1_field_sync (field_vec%id)
    call m3dc1_field_printcompnorm(field_vec%id, "field_vec after sync "//char(0))
    call m3dc1_field_printcompnorm(field0_vec%id, "field0_vec before sync "//char(0))
    call m3dc1_field_sync (field0_vec%id)
    call m3dc1_field_printcompnorm(field0_vec%id, "field0_vec after sync "//char(0))
  endif

  ! If we are running a linear simulation, but the restart file was
  ! a nonlinear simulation, make the restart data be the equilibrium
  if(linear.eq.1 .and. prev_linear.eq.0) then
     !field0_vec%data = field0_vec%data + field_vec%data
     call m3dc1_field_add(field0_vec%id, field_vec%id)
     !field_vec%data = 0.
     call m3dc1_field_assign (field_vec%id, 0., 0)
  endif 

  if (prev_ndofs2==cur_ndofs2) then
    !bf_field(1)%vec%data(1:prev_ndofs2) = tmp_bf_field_1(1:prev_ndofs2)
    call m3dc1_field_set(bf_field(1)%vec%id, tmp_bf_field_1, prev_ndofs2)
    !bf_field(0)%vec%data(1:prev_ndofs2) = tmp_bf_field_0(1:prev_ndofs2) 
    call m3dc1_field_set(bf_field(0)%vec%id, tmp_bf_field_0, prev_ndofs2)
  else
    do i=1,prev_nnodes
      dofs_node1 =0.
      dofs_node2 =0.
      !bf_field(1)%vec
      dofs_node1(1:prev_ndofs2_pernode)& 
        = tmp_bf_field_1((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
      !bf_field(0)%vec
      dofs_node2(1:prev_ndofs2_pernode)&                                                                     
        = tmp_bf_field_0((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
      call m3dc1_ent_setdofdata (0, i-1, bf_field(1)%vec%id, cur_ndofs2_pernode, dofs_node1(1:cur_ndofs2_pernode))
      call m3dc1_ent_setdofdata (0, i-1, bf_field(0)%vec%id, cur_ndofs2_pernode, dofs_node2(1:cur_ndofs2_pernode))
    enddo
    call m3dc1_field_sync (bf_field(0)%vec%id)
    call m3dc1_field_sync (bf_field(1)%vec%id)
  endif

  if (prev_useext.eq.1) then 
    use_external_fields = .true.
    call create_field(psi_ext)
    call create_field(bz_ext)
    call create_field(bf_ext)
    if (prev_ndofs2==cur_ndofs2) then
      !psi_ext%vec%data(1:prev_ndofs2) = tmp_psi_ext(1:prev_ndofs2)
      call m3dc1_field_set(psi_ext%vec%id, tmp_psi_ext, prev_ndofs2)
      !bz_ext%vec%data(1:prev_ndofs2) = tmp_bz_ext(1:prev_ndofs2)
      call m3dc1_field_set(bz_ext%vec%id, tmp_bz_ext, prev_ndofs2)
      !bf_ext%vec%data(1:prev_ndofs2) = tmp_bf_ext(1:prev_ndofs2)
      call m3dc1_field_set(bf_ext%vec%id, tmp_bf_ext, prev_ndofs2)
    else !(prev_ndofs2!=cur_ndofs2)
      do i=1,prev_nnodes
        dofs_node1 =0.
        dofs_node2 =0.
        dofs_node3 =0.
        !psi_ext%vec
        dofs_node1(1:prev_ndofs2_pernode) & 
          = tmp_psi_ext((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
        !bz_ext%vec
        dofs_node2(1:prev_ndofs2_pernode) & 
          = tmp_bz_ext((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
        !bf_ext%vec
        dofs_node3(1:prev_ndofs2_pernode) &
          = tmp_bf_ext((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
        call m3dc1_ent_setdofdata (0, i-1, psi_ext%vec%id, cur_ndofs2_pernode, dofs_node1(1:cur_ndofs2_pernode))
        call m3dc1_ent_setdofdata (0, i-1, bz_ext%vec%id, cur_ndofs2_pernode, dofs_node2(1:cur_ndofs2_pernode))
        call m3dc1_ent_setdofdata (0, i-1, bf_ext%vec%id, cur_ndofs2_pernode, dofs_node3(1:cur_ndofs2_pernode))
      enddo
      call m3dc1_field_sync (psi_ext%vec%id)
      call m3dc1_field_sync (bz_ext%vec%id)
      call m3dc1_field_sync (bf_ext%vec%id)
    endif !if (prev_ndofs2==cur_ndofs2) then
  end if !if (prev_useext.eq.1) 

  if (icsubtract.eq.1) then
    if (prev_ndofs2==cur_ndofs2) then
      call m3dc1_field_set(psi_coil_field%vec%id, tmp_psi_coil_field, prev_ndofs2)
    else
      do i=1,prev_nnodes
        dofs_node1 =0.
        dofs_node2 =0.
        dofs_node3 =0.
        !psi_ext%vec
        dofs_node1(1:prev_ndofs2_pernode) &
          = tmp_psi_coil_field((i-1)*prev_ndofs2_pernode+1:(i-1)*prev_ndofs2_pernode+prev_ndofs2_pernode)
        call m3dc1_ent_setdofdata (0, i-1, psi_coil_field%vec%id, cur_ndofs2_pernode, dofs_node1(1:cur_ndofs2_pernode))
      end do
      call m3dc1_field_printcompnorm(psi_coil_field%vec%id, "psi_coil_field before sync "//char(0))
      call m3dc1_field_sync (psi_coil_field%vec%id)
      call m3dc1_field_printcompnorm(psi_coil_field%vec%id, "psi_coil_field after sync "//char(0))
    end if
  end if
  deallocate(tmp_field_vec, tmp_field0_vec, tmp_bf_field_1, tmp_bf_field_0, &
       tmp_psi_ext, tmp_bz_ext, tmp_bf_ext, tmp_psi_coil_field) 
  if (prev_ndofs1 .ne. cur_ndofs1 .and. itaylor .eq. 1 .and. itor .eq. 1) then
     if(myrank .eq. 0) print *, "call gradshafranov_per after restart"
     call gradshafranov_per()
  end if
#endif
end subroutine rdrestart_adios

module newvar_mod

  implicit none

  integer, allocatable :: ibound_newvar(:)
  integer :: nbc_newvar

contains


! create_newvar_matrix
! ====================
subroutine create_newvar_matrix(matrix, ibound)

  use p_data
  use t_data
  use arrays
  use sparse

  implicit none

  integer, intent(in) :: matrix, ibound

  integer :: numelms, numnodes, nbound, itri, i, j, ione, jone, izone, izonedim
  real :: fintl(-6:maxi,-6:maxi), dterm(18,18)

  call numnod(numnodes)
  call numfac(numelms)

  ! calculate number of boundary nodes and allocate array
  if(ibound.eq.1) then
     nbound = 0
     do i=1,numnodes
        call zonenod(i,izone,izonedim)
        if(izonedim .ne. 2) nbound = nbound + 1
     enddo
     nbound = 12*nbound    

     allocate(ibound_newvar(nbound))
  endif

  ! populate matrix
  call zeroarray4solve(matrix, numvar1_numbering)
  do itri=1,numelms
     call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
     call calcdterm(itri, dterm, fintl)
     do j=1,18
        jone = isval1(itri,j)
        do i=1,18
           ione = isval1(itri,i)
           call insertval(matrix, dterm(i,j), ione, jone, 1)
        enddo
     enddo
  enddo

  ! apply boundary conditions
  if(ibound.eq.1) then
     call boundaryds(ibound_newvar,nbc_newvar,0)
     do i=1,nbc_newvar
        call setdiribc(matrix, ibound_newvar(i))
     enddo
  end if

  call finalizearray4solve(matrix)

end subroutine create_newvar_matrix


! define_sources
! ==============
subroutine define_sources()

  use basic
  use t_data
  use arrays
  use nintegrate_mod
  use metricterms_n

  implicit none
  
  integer :: itri, numelms, i, ione, ndof, def_fields
  real :: x, z, xmin, zmin, hyp, hypi, hypv, hypc, hypp, deex, dbf, factor, temp

  double precision :: cogcoords(3)

  real, dimension(20) :: avec

  sb1 = 0
  if(numvar.ge.2) sb2 = 0
  if(numvar.ge.3) sp1 = 0

  ! Determine which fields need to be calculated
  def_fields = FIELD_PSI
  if(numvar.ge.2) def_fields = def_fields + FIELD_I
  if(numvar.ge.3) def_fields = def_fields + FIELD_PHI + FIELD_V + FIELD_CHI + &
          FIELD_PE + FIELD_P + FIELD_B2I + FIELD_J + FIELD_COM + FIELD_VOR
  if(idens.eq.1) def_fields = def_fields + FIELD_NI 

  call getmincoord(xmin,zmin)
  
  ! Define RHS of equation
  call numfac(numelms)
  do itri=1,numelms

     if(imask.eq.1) then
        call cogfac(itri,cogcoords)
        x = cogcoords(1)-xmin
        z = cogcoords(2)-zmin
        call mask(x,z,factor)
     else
        factor = 1.
     endif
     dbf = db*factor

     call getdeex(itri,deex)
     hyp  = hyper *deex**2
     hypi = hyperi*deex**2
     hypv = hyperv*deex**2
     hypc = hyperc*deex**2
     hypp = hyperp*deex**2
     call define_fields_79(itri, def_fields)
    
     do i=1,18
        ione = isval1(itri,i)
        
        sb1(ione) = sb1(ione) + b1psieta(g79(:,:,i),pst79,etar,hyp)

        if(numvar.ge.2) then
           sb1(ione) = sb1(ione) + b1psibd(g79(:,:,i),pst79,bzt79,ni79)*dbf

           sb2(ione) = sb2(ione)  &
                + b2psipsid(g79(:,:,i),pst79,pst79,ni79)*dbf &
                + b2bbd    (g79(:,:,i),bzt79,bzt79,ni79)*dbf &
                + b2beta   (g79(:,:,i),bzt79,etar,hypi)
        endif

        if(numvar.ge.3) then
           sb2(ione) = sb2(ione) + b2ped(g79(:,:,i),pet79,ni79)*dbf*pefac

           sp1(ione) = sp1(ione) &
                + b3psipsieta(g79(:,:,i),pst79,pst79)*etar   &
                + b3bbeta    (g79(:,:,i),bzt79,bzt79)*etar   &
                + b3pedkappa (g79(:,:,i),pt79, ni79,kappat,hypp) &
                + p1kappar   (g79(:,:,i),pst79,pst79,pt79,ni79,b2i79)*kappar &
                + b3pebd(g79(:,:,i),pet79,bzt79,ni79)*dbf*pefac

           ! ohmic heating         
           sp1(ione) = sp1(ione) + (gam-1.)* &
                (qpsipsieta(g79(:,:,i),pst79,pst79,etar,hyp,jt79) &
                +qbbeta    (g79(:,:,i),bzt79,bzt79,etar,hypi))

           ! viscous heating
           sp1(ione) = sp1(ione) - (gam-1.)* &
                (quumu    (g79(:,:,i),pht79,pht79,amu,amuc,hypc) &
                +qvvmu    (g79(:,:,i),vzt79,vzt79,amu,     hypv) &
                +quchimu  (g79(:,:,i),pht79,cht79,amu,amuc,hypc) &
                +0.*qchichimu(g79(:,:,i),cht79,cht79,amu,amuc,hypc))
        endif ! on numvar.ge.3

     end do

  end do 


  ! Solve equations
  call numdofs(1, ndof)

  call solve_newvar(sb1, 1, ndof)
  if(numvar.ge.2) call solve_newvar(sb2, 1, ndof)
  if(numvar.ge.3) call solve_newvar(sp1, 1, ndof)

end subroutine define_sources


! newvar_gs
! =========
subroutine newvar_gs(inarray,outarray,numvard,itype,ibound)

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: numvard, itype, ibound
  real, intent(in) :: inarray(*) ! length using numvard ordering
  real, intent(out) :: outarray(*) ! length using numvar=1 ordering

!!$  real, allocatable :: temp(:)

  integer :: ndof, ier, numelms, itri, i, j, ione, j1
  real :: sum

!!$  call createvec(temp, 1)
!!$
!!$  temp = 0
  call numdofs(1, ndof)
  do i=1, ndof
     outarray(i) = 0.
  end do

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     ! calculate the local sampling points and weights for numerical integration
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
     end do

     do i=1,18
        ione = isval1(itri,i)
        sum = 0.
        do j=1,18
           j1 = isvaln(itri,j)
           sum = sum - inarray(j1 + 6*(itype-1)) * &
                (int2(g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79) &
                +int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79))
           if(itor.eq.1) then
              sum = sum - inarray(j1 + 6*(itype-1)) * &
                   2.*int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),weight_79,79)
           endif
        end do
!!$        temp(ione) = temp(ione) + sum
        outarray(ione) = outarray(ione) + sum
     end do
  enddo

  ! solve linear equation
  call solve_newvar(outarray, ibound, ndof)
  
!!$  call deletevec(temp)

end subroutine newvar_gs


! solve_newvar
! ============
subroutine solve_newvar(rhs, ibound, ndof)

  use sparse

  implicit none

  integer, intent(in) :: ibound, ndof
  real, dimension(*), intent(inout) :: rhs

  integer :: i, ier

  call sumshareddofs(rhs)

  if(ibound.eq.1) then
     do i=1,nbc_newvar
        rhs(ibound_newvar(i)) = 0.
     enddo
     call solve(s6matrix_sm,rhs,ier)
  else
     call solve(s3matrix_sm,rhs,ier)
  endif

end subroutine solve_newvar

end module newvar_mod

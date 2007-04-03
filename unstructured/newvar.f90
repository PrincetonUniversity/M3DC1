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
  use nintegrate_mod

  implicit none

  integer, intent(in) :: matrix, ibound

  integer :: numelms, numnodes, nbound, itri, i, j, ione, jone, izone, izonedim
  real :: temp

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

     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
     end do

     if(ijacobian.eq.1) weight_79 = weight_79 * r_79

     do j=1,18
        jone = isval1(itri,j)
        do i=1,18
           ione = isval1(itri,i)
           temp = int2(g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)
           call insertval(matrix, temp, ione, jone, 1)
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

  include 'mpif.h'
  
  integer :: itri, numelms, i, ione, ndof, def_fields
  real :: x, z, xmin, zmin, hypf, hypi, hypv, hypc, hypp, dbf, factor

  double precision, dimension(3)  :: cogcoords
  double precision, dimension(18) :: temp, temp2

!!$  integer :: izone, izonedim

  sb1 = 0.
  if(numvar.ge.2) sb2 = 0.
  if(numvar.ge.3) sp1 = 0.

!!$  tempvar = 0.

  ekino = ekin
  emago = emag
  ekindo = ekind
  emagdo = emagd
!
  ekinpo = ekinp
  emagpo = emagp
  ekinpdo = ekinpd
  emagpdo = emagpd
!
  ekinto = ekint
  emagto = emagt
  ekintdo = ekintd
  emagtdo = emagtd
!
  ekinpho = ekinph
  ekintho = ekinth
  emagpho = emagph
  emagtho = emagth
!
  ekin3o = ekin3 
  ekin3do = ekin3d
  ekin3ho = ekin3h 
  emag3o = emag3
  emag3do = emag3d
  emag3ho = emag3h
!  
  ekin = 0.
  emag = 0.
  ekind = 0.
  emagd = 0.
  ekinp = 0.
  emagp = 0.
  ekinpd = 0.
  emagpd = 0.      
  ekint = 0.
  emagt = 0.
  ekintd = 0.
  emagtd = 0.
  ekinph = 0.
  ekinth = 0.
  emagph = 0.
  emagth = 0.
  ekin3 = 0.
  ekin3d = 0.
  ekin3h = 0.
  emag3 = 0.
  emag3d = 0.
  emag3h = 0.


  hypf = hyper *deex**2
  hypi = hyperi*deex**2
  hypv = hyperv*deex**2
  hypc = hyperc*deex**2
  hypp = hyperp*deex**2

  ! Specify which fields need to be calculated
  def_fields = FIELD_PSI + FIELD_PHI + FIELD_VOR + FIELD_J + FIELD_ETA
  if(numvar.ge.2) def_fields = def_fields + FIELD_I + FIELD_V
  if(numvar.ge.3) def_fields = def_fields + FIELD_CHI + &
          FIELD_PE + FIELD_P + FIELD_B2I + FIELD_COM
  if(idens.eq.1) def_fields = def_fields + FIELD_N + FIELD_NI 

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

     call define_fields_79(itri, def_fields)

     if(isources.eq.1) then
   
     do i=1,18
        ione = isval1(itri,i)

!!$        tempvar(ione) = tempvar(ione) &
!!$             + v1psipsi(g79(:,:,i), pst79, pst79) &
!!$             + v1bb(g79(:,:,i), bzt79, bzt79)
        
        ! Definition of Source Terms
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        sb1(ione) = sb1(ione) + b1psieta(g79(:,:,i),pst79,eta79,hypf)

        if(numvar.ge.2) then
           sb1(ione) = sb1(ione) + b1psibd(g79(:,:,i),pst79,bzt79,ni79)*dbf

           sb2(ione) = sb2(ione)  &
                + b2psipsid(g79(:,:,i),pst79,pst79,ni79)*dbf &
                + b2bbd    (g79(:,:,i),bzt79,bzt79,ni79)*dbf &
                + b2beta   (g79(:,:,i),bzt79,eta79,hypi)
        endif

        if(numvar.ge.3) then
           sb2(ione) = sb2(ione) + b2ped(g79(:,:,i),pet79,ni79)*dbf*pefac

           sp1(ione) = sp1(ione) &
                + b3psipsieta(g79(:,:,i),pst79,pst79,eta79)   &
                + b3bbeta    (g79(:,:,i),bzt79,bzt79,eta79)   &
!!$                + b3pedkappa (g79(:,:,i),pt79,ni79,kappat,hypp)*(gam-1.) &
                + p1kappar   (g79(:,:,i),pst79,pst79,pet79,ni79,b2i79)*kappar*(gam-1.) &
                + b3pebd(g79(:,:,i),pet79,bzt79,ni79)*dbf*pefac

           ! ohmic heating         
           sp1(ione) = sp1(ione) + (gam-1.)* &
                (qpsipsieta(g79(:,:,i),pst79,pst79,eta79,hypf,jt79) &
                +qbbeta    (g79(:,:,i),bzt79,bzt79,eta79,hypi))

           ! viscous heating
           sp1(ione) = sp1(ione) - (gam-1.)* &
                (quumu    (g79(:,:,i),pht79,pht79,amu,amuc,hypc) &
                +qvvmu    (g79(:,:,i),vzt79,vzt79,amu,     hypv) &
                +quchimu  (g79(:,:,i),pht79,cht79,amu,amuc,hypc) &
                +0.*qchichimu(g79(:,:,i),cht79,cht79,amu,amuc,hypc))
        endif ! on numvar.ge.3
     end do

     endif ! on isources

     if(ijacobian.eq.1) weight_79 = weight_79 * ri_79

     ! Definition of energy
     ! ~~~~~~~~~~~~~~~~~~~~
     if(idens.eq.0) then
        ekinp = ekinp + .5* &
             (int3(r2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),weight_79,79) &
             +int3(r2_79,pht79(:,OP_DR),pht79(:,OP_DR),weight_79,79))
     else
        ekinp = ekinp + .5* &
             (int4(r2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             +int4(r2_79,pht79(:,OP_DR),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79))
     endif

     ekinpd = ekinpd - amu*int3(r2_79,pht79(:,OP_LP ),pht79(:,OP_LP ),weight_79,79)
        
     ekinph = ekinph - hypc*amu* &
          (int2(vot79(:,OP_DZ),vot79(:,OP_DZ),weight_79,79) &
          +int2(vot79(:,OP_DR),vot79(:,OP_DR),weight_79,79))

     emagp = emagp + .5* &
          (int3(ri2_79,pst79(:,OP_DZ),pst79(:,OP_DZ),weight_79,79) &
          +int3(ri2_79,pst79(:,OP_DR),pst79(:,OP_DR),weight_79,79))

     emagpd = emagpd - &
          int4(ri2_79,pst79(:,OP_GS),pst79(:,OP_GS),eta79(:,OP_1),weight_79,79)

     emagph = emagph - hypf* &
          (int3(jt79(:,OP_DZ),jt79(:,OP_DZ),eta79(:,OP_1),weight_79,79) &
          +int3(jt79(:,OP_DR),jt79(:,OP_DR),eta79(:,OP_1),weight_79,79))

     if(numvar.ge.2) then
        if(idens.eq.0) then
           ekint = ekint + .5*int2(vzt79(:,OP_1),vzt79(:,OP_1),weight_79,79)
        else
           ekint = ekint + .5*int3(vzt79(:,OP_1),vzt79(:,OP_1),nt79(:,OP_1),weight_79,79)
        endif

        ekintd = ekintd - amu* &
             (int2(vzt79(:,OP_DZ),vzt79(:,OP_DZ),weight_79,79) &
             +int2(vzt79(:,OP_DR),vzt79(:,OP_DR),weight_79,79))

        ekinth = ekinth - amu*hypv*int2(vzt79(:,OP_LP),vzt79(:,OP_LP),weight_79,79)
        
        emagt = emagt + .5*int3(ri2_79,bzt79(:,OP_1),bzt79(:,OP_1),weight_79,79)
           
        emagtd = emagtd - &
             (int3(bzt79(:,OP_DZ),bzt79(:,OP_DZ),eta79(:,OP_1),weight_79,79) &
             +int3(bzt79(:,OP_DR),bzt79(:,OP_DR),eta79(:,OP_1),weight_79,79))
           
        emagth = emagth - &
             hypi*int3(bzt79(:,OP_LP),bzt79(:,OP_LP),eta79(:,OP_1),weight_79,79)
     endif

     if(numvar.ge.3) then
        if(idens.eq.0) then
           ekin3 = ekin3 + .5* &
                (int2(cht79(:,OP_DZ),cht79(:,OP_DZ),weight_79,79) &
                +int2(cht79(:,OP_DR),cht79(:,OP_DR),weight_79,79))
           if(itor.eq.1) then
              ekin3 = ekin3 + &
                   (int3(r_79,cht79(:,OP_DZ),pht79(:,OP_DR),weight_79,79) &
                   -int3(r_79,cht79(:,OP_DR),pht79(:,OP_DZ),weight_79,79))
           endif
        else
           ekin3 = ekin3 + .5* &
                (int3(cht79(:,OP_DZ),cht79(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
                +int3(cht79(:,OP_DR),cht79(:,OP_DR),nt79(:,OP_1),weight_79,79))
           if(itor.eq.1) then
              ekin3 = ekin3 + &
                   (int4(r_79,cht79(:,OP_DZ),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
                   -int4(r_79,cht79(:,OP_DR),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79))
           endif
        endif
           
        ekin3d = ekin3d - amuc*int2(pht79(:,OP_LP),pht79(:,OP_LP),weight_79,79)
           
        ekin3h = ekin3h - hypc*amuc* &
             (int2(cot79(:,OP_DZ),cot79(:,OP_DZ),weight_79,79) &
             +int2(cot79(:,OP_DR),cot79(:,OP_DR),weight_79,79))
        
        emag3 = emag3 + int1(pt79,weight_79,79) / (gam - 1.)
           
     endif
  end do

  ! Solve source term equations
  call numdofs(1, ndof)

  if(isources.eq.1) then
     call solve_newvar(sb1, 1)
     if(numvar.ge.2) call solve_newvar(sb2, 1)
     if(numvar.ge.3) call solve_newvar(sp1, 1)
  endif
!!$ 
!!$  call solve_newvar(tempvar, 1)

  ! Allreduce energy terms
  if(maxrank .gt. 1) then
     temp(1) = ekinp
     temp(2) = emagp
     temp(3) = ekinpd
     temp(4) = emagpd      
     temp(5) = ekint
     temp(6) = emagt
     temp(7) = ekintd
     temp(8) = emagtd
     temp(9) = ekinph
     temp(10) = ekinth
     temp(11) = emagph
     temp(12) = emagth
     temp(13) = ekin3
     temp(14) = ekin3d
     temp(15) = ekin3h
     temp(16) = emag3
     temp(17) = emag3d
     temp(18) = emag3h
         
     !checked that this should be MPI_DOUBLE_PRECISION
     call mpi_allreduce(temp, temp2, 18, MPI_DOUBLE_PRECISION,  &
          MPI_SUM, MPI_COMM_WORLD, i) 
         
     ekinp = temp2(1)
     emagp = temp2(2)
     ekinpd = temp2(3)
     emagpd = temp2(4)      
     ekint = temp2(5)
     emagt = temp2(6)
     ekintd = temp2(7)
     emagtd = temp2(8)
     ekinph = temp2(9)
     ekinth = temp2(10)
     emagph = temp2(11)
     emagth = temp2(12)
     ekin3 = temp2(13)
     ekin3d = temp2(14)
     ekin3h = temp2(15)
     emag3 = temp2(16)
     emag3d = temp2(17)
     emag3h = temp2(18)
  endif !if maxrank .gt. 1

  ekin = ekinp + ekint + ekin3
  emag = emagp + emagt + emag3
  ekind = ekinpd + ekintd + ekin3d
  emagd = emagpd + emagtd + emag3d

  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, "Energy at ntime = ", ntime
     print *, "ekinp, ekint, ekin3 = ", ekinp, ekint, ekin3
     print *, "ekinpd, ekintd, ekin3d = ", ekinpd, ekintd, ekin3d
     print *, "ekinph, ekinth, ekin3h = ", ekinph, ekinth, ekin3h
     print *, "emagp, emagt, emag3 = ", emagp, emagt, emag3
     print *, "emagpd, emagd, emag3d = ", emagpd, emagtd, emag3d
     print *, "emagph, emagth, emag3h = ", emagph, emagth, emag3h
  endif

end subroutine define_sources


! newvar_gs
! =========
subroutine newvar_gs(inarray,outarray,itype,ibound)

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itype, ibound
  real, intent(in) :: inarray(*) ! length using numvard ordering
  real, intent(out) :: outarray(*) ! length using numvar=1 ordering

  integer :: ndof, numelms, itri, i, j, ione, j1
  real :: sum

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

     if(ijacobian.eq.1) weight_79 = weight_79*r_79

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
        outarray(ione) = outarray(ione) + sum
     end do
  enddo

  ! solve linear equation
  call solve_newvar(outarray, ibound)

end subroutine newvar_gs


! newvar_eta
! ==========
subroutine newvar_eta()

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  implicit none

  real :: minpe
  integer :: i, ione, itri, numelms, def_fields

!!$  minpe = p0/100.
  minpe = 0.

  resistivity = 0.
  def_fields = 0.
  if(idens.eq.1)  def_fields = def_fields + FIELD_N
  if(numvar.ge.3) def_fields = def_fields + FIELD_PE

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

     call define_fields_79(itri, def_fields)

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
     end do

     do i=1,18
        ione = isval1(itri,i)       

        if(numvar.ge.3) then
           temp79b = max(pefac*pet79(:,OP_1),minpe)
                     
           if(idens.eq.0) then
              temp79a = sqrt((1./(pefac*pet79(:,OP_1) + minpe))**3)
!!$              temp79a = sqrt(1./temp79b**3)
           else
              temp79a = sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1) + minpe))**3)
!!$              temp79a = sqrt((nt79(:,OP_1)/temp79b)**3)
           endif
        else
           if(idens.eq.0) then
              temp79a = sqrt((1./(p0-pi0 + minpe))**3)
!!$              temp79a = sqrt((1./max(p0 - pi0,minpe))**3)
           else
              temp79a = sqrt((nt79(:,OP_1)/(p0-pi0 + minpe))**3)
!!$              temp79a = sqrt((nt79(:,OP_1)/max(p0 - pi0,minpe))**3)
           endif
        endif

        resistivity(ione) = resistivity(ione) &
             + etar*int1(g79(:,OP_1,i),weight_79,79) &
             + eta0*int2(g79(:,OP_1,i),temp79a, weight_79,79)
     end do
  enddo

  ! solve linear equation
  call solve_newvar(resistivity, 0)

end subroutine newvar_eta


! solve_newvar
! ============
subroutine solve_newvar(rhs, ibound)

  use sparse

  implicit none

  integer, intent(in) :: ibound
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

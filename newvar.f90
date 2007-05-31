module newvar_mod

  implicit none

  integer, parameter :: NV_NOBOUND = 0
  integer, parameter :: NV_DCBOUND = 1

  integer, parameter :: NV_LP = 0
  integer, parameter :: NV_GS = 1


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

  integer :: numelms, itri, i, j, ione, jone, izone, izonedim
  real :: temp
  real, allocatable :: rhs(:)

  call numfac(numelms)

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
  if(ibound.eq.NV_DCBOUND) then
     call createvec(rhs, 1)
     call boundary_dc(matrix, rhs)
     call deletevec(rhs)
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

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none

  include 'mpif.h'
  
  integer :: itri, numelms, i, ione, ndof, def_fields
  real :: x, z, xmin, zmin, hypf, hypi, hypv, hypc, hypp, dbf, factor

  double precision, dimension(3)  :: cogcoords
  double precision, dimension(23) :: temp, temp2

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
  efluxd = 0.
  efluxp = 0.
  efluxk = 0.
  efluxs = 0.
  efluxt = 0.

  hypf = hyper *deex**2
  hypi = hyperi*deex**2
  hypv = hyperv*deex**2
  hypc = hyperc*deex**2
  hypp = hyperp*deex**2

  ! Specify which fields need to be calculated
  def_fields = FIELD_PSI + FIELD_PHI + FIELD_J + FIELD_ETA
  if(numvar.ge.2) def_fields = def_fields + FIELD_I + FIELD_V
  if(numvar.ge.3) then
     def_fields = def_fields + FIELD_CHI + &
          FIELD_PE + FIELD_P
     if(kappar.ne.0) def_fields = def_fields + FIELD_B2I
  endif
  if(idens.eq.1) def_fields = def_fields + FIELD_N

  if(isources.eq.1) then
     if(idens.eq.1) def_fields = def_fields + FIELD_NI
  endif

  if(hypc.ne.0.) then 
     def_fields = def_fields + FIELD_VOR
     if(numvar.ge.3) def_fields = def_fields + FIELD_COM
  end if

  tm79 = 0.
  tm79(:,OP_1) = 1.

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
                (qpsipsieta(g79(:,:,i),hypf) &
                +qbbeta    (g79(:,:,i),hypi))

           ! viscous heating
           sp1(ione) = sp1(ione) - (gam-1.)* &
                (quumu    (g79(:,:,i),pht79,pht79,amu,amuc,hypc) &
                +qvvmu    (g79(:,:,i),vzt79,vzt79,amu,     hypv) &
                +quchimu  (g79(:,:,i),pht79,cht79,amu,amuc,hypc) &
                +qchichimu(g79(:,:,i),cht79,cht79,amu,amuc,hypc))
        endif ! on numvar.ge.3
     end do

     endif ! on isources

     ! Definition of energy
     ! ~~~~~~~~~~~~~~~~~~~~

#ifdef NEW_VELOCITY
     if(idens.eq.0) then
        ekinp = ekinp + .5* &
             (int3(ri2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),weight_79,79) &
             +int3(ri2_79,pht79(:,OP_DR),pht79(:,OP_DR),weight_79,79))
     else
        ekinp = ekinp + .5* &
             (int4(ri2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             +int4(ri2_79,pht79(:,OP_DR),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79))
     endif

     ekinpd = ekinpd - amu*int3(ri2_79,pht79(:,OP_GS),pht79(:,OP_GS),weight_79,79)
        
     ekinph = ekinph - hypc*amu* &
          (int3(ri2_79,vot79(:,OP_DZ),vot79(:,OP_DZ),weight_79,79) &
          +int3(ri2_79,vot79(:,OP_DR),vot79(:,OP_DR),weight_79,79))
#else
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
          (int3(r2_79,vot79(:,OP_DZ),vot79(:,OP_DZ),weight_79,79) &
          +int3(r2_79,vot79(:,OP_DR),vot79(:,OP_DR),weight_79,79))
#endif

     emagp = emagp + .5* &
          (int3(ri2_79,pst79(:,OP_DZ),pst79(:,OP_DZ),weight_79,79) &
          +int3(ri2_79,pst79(:,OP_DR),pst79(:,OP_DR),weight_79,79))

     emagpd = emagpd - &
          int4(ri2_79,pst79(:,OP_GS),pst79(:,OP_GS),eta79(:,OP_1),weight_79,79)

     emagph = emagph - qpsipsieta(tm79,hypf)

     if(numvar.ge.2) then

#ifdef NEW_VELOCITY
        if(idens.eq.0) then
           ekint = ekint + &
                .5*int3(ri2_79,vzt79(:,OP_1),vzt79(:,OP_1),weight_79,79)
        else
           ekint = ekint + &
                .5*int4(ri2_79,vzt79(:,OP_1),vzt79(:,OP_1),nt79(:,OP_1),weight_79,79)
        endif

        ekintd = ekintd - amu* &
             (int3(ri2_79,vzt79(:,OP_DZ),vzt79(:,OP_DZ),weight_79,79) &
             +int3(ri2_79,vzt79(:,OP_DR),vzt79(:,OP_DR),weight_79,79))

        ekinth = ekinth - amu*hypv*int3(ri2_79,vzt79(:,OP_GS),vzt79(:,OP_GS),weight_79,79)
#else
        if(idens.eq.0) then
           ekint = ekint + .5*int2(vzt79(:,OP_1),vzt79(:,OP_1),weight_79,79)
        else
           ekint = ekint + .5*int3(vzt79(:,OP_1),vzt79(:,OP_1),nt79(:,OP_1),weight_79,79)
        endif

        ekintd = ekintd - amu* &
             (int2(vzt79(:,OP_DZ),vzt79(:,OP_DZ),weight_79,79) &
             +int2(vzt79(:,OP_DR),vzt79(:,OP_DR),weight_79,79))

        ekinth = ekinth - amu*hypv*int2(vzt79(:,OP_GS),vzt79(:,OP_GS),weight_79,79)
#endif
        
        emagt = emagt + .5*int3(ri2_79,bzt79(:,OP_1),bzt79(:,OP_1),weight_79,79)
           
        emagtd = emagtd - &
             (int4(ri2_79,bzt79(:,OP_DZ),bzt79(:,OP_DZ),eta79(:,OP_1),weight_79,79) &
             +int4(ri2_79,bzt79(:,OP_DR),bzt79(:,OP_DR),eta79(:,OP_1),weight_79,79))

        emagth = emagth - qbbeta(tm79,hypi)
     endif

     if(numvar.ge.3) then
#ifdef NEW_VELOCITY
        if(idens.eq.0) then
           ekin3 = ekin3 + .5* &
                (int2(cht79(:,OP_DZ),cht79(:,OP_DZ),weight_79,79) &
                +int2(cht79(:,OP_DR),cht79(:,OP_DR),weight_79,79))
           if(itor.eq.1) then
              ekin3 = ekin3 + &
                   (int3(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),weight_79,79) &
                   -int3(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),weight_79,79))
           endif
        else
           ekin3 = ekin3 + .5* &
                (int3(cht79(:,OP_DZ),cht79(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
                +int3(cht79(:,OP_DR),cht79(:,OP_DR),nt79(:,OP_1),weight_79,79))
           if(itor.eq.1) then
              ekin3 = ekin3 + &
                   (int4(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),nt79(:,OP_1),weight_79,79) &
                   -int4(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),nt79(:,OP_1),weight_79,79))
           endif
        endif

        ekin3d = ekin3d - 2.*amuc*int2(pht79(:,OP_LP),pht79(:,OP_LP),weight_79,79)
           
        ekin3h = ekin3h - 2.*hypc*amuc* &
             (int2(cot79(:,OP_DZ),cot79(:,OP_DZ),weight_79,79) &
             +int2(cot79(:,OP_DR),cot79(:,OP_DR),weight_79,79))
#else
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
#endif
                   
        emag3 = emag3 + int1(pt79,weight_79,79) / (gam - 1.)

     endif

     ! Calculate fluxes through boundary of domain
     efluxd = efluxd + flux_diffusive()
     efluxk = efluxk +flux_ke()
     efluxp = efluxp +flux_pressure(dbf)
     efluxs = efluxs +flux_poynting(dbf)
     efluxt = efluxt +flux_heat()
  end do

  ! Solve source term equations
  call numdofs(1, ndof)

  if(isources.eq.1) then
     call solve_newvar(sb1, NV_DCBOUND)
     if(numvar.ge.2) call solve_newvar(sb2, NV_DCBOUND)
     if(numvar.ge.3) call solve_newvar(sp1, NV_DCBOUND)
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
     temp(19) = efluxd
     temp(20) = efluxp
     temp(21) = efluxk
     temp(22) = efluxs
     temp(23) = efluxt
         
     !checked that this should be MPI_DOUBLE_PRECISION
     call mpi_allreduce(temp, temp2, 23, MPI_DOUBLE_PRECISION,  &
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
     efluxd = temp2(19)
     efluxp = temp2(20)
     efluxk = temp2(21)
     efluxs = temp2(22)
     efluxt = temp2(23)
  endif !if maxrank .gt. 1

  ekin = ekinp + ekint + ekin3
  emag = emagp + emagt + emag3
  ekind = ekinpd + ekintd + ekin3d
  emagd = emagpd + emagtd + emag3d

!!$  if(myrank.eq.0 .and. iprint.ge.1) then
!!$     print *, "Energy at ntime = ", ntime
!!$     print *, "ekinp, ekint, ekin3 = ", ekinp, ekint, ekin3
!!$     print *, "ekinpd, ekintd, ekin3d = ", ekinpd, ekintd, ekin3d
!!$     print *, "ekinph, ekinth, ekin3h = ", ekinph, ekinth, ekin3h
!!$     print *, "emagp, emagt, emag3 = ", emagp, emagt, emag3
!!$     print *, "emagpd, emagd, emag3d = ", emagpd, emagtd, emag3d
!!$     print *, "emagph, emagth, emag3h = ", emagph, emagth, emag3h
!!$  endif

end subroutine define_sources


! newvar_d2
! =========
subroutine newvar_d2(inarray,outarray,itype,ibound,gs)

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itype, ibound
  real, intent(in) :: inarray(*) ! length using numvard ordering
  real, intent(out) :: outarray(*) ! length using numvar=1 ordering
  integer, intent(in) :: gs ! NV_GS for grad-shafranov operator, NV_LP for laplacian

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
           if(itor.eq.1 .and. gs.eq.NV_GS) then
              sum = sum - inarray(j1 + 6*(itype-1)) * &
                   2.*int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),weight_79,79)
           endif
        end do
        outarray(ione) = outarray(ione) + sum
     end do
  enddo

  ! solve linear equation
  call solve_newvar(outarray, ibound)

end subroutine newvar_d2


! newvar_eta
! ==========
subroutine newvar_eta()

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  use gradshafranov

  implicit none

  real :: minpe, ajlim, xguess, zguess
  integer :: i, itri
  integer :: ione, numelms, def_fields
!!$  integer :: numnodes, ibegin, iendplusone, ibegin1, iendplusone1

  resistivity = 0.
  def_fields = 0.
  if(idens.eq.1)  def_fields = def_fields + FIELD_N
  if(numvar.ge.3) def_fields = def_fields + FIELD_PE

  if(itor.eq.1 .and. numvar.lt.3) &
       def_fields = def_fields + FIELD_PSI

  if(itor.eq.1 .and. numvar.lt.3) then
     itri = 0.
     call evaluate(xlim-xzero,zlim-zzero,psilim,ajlim,phi,1,numvar,itri)

     xguess = xmag - xzero
     zguess = zmag - zzero
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call magaxis(xguess,zguess,phi+phi0,numvar)
     else
        call magaxis(xguess,zguess,phi,numvar)
     endif
     xmag = xguess + xzero
     zmag = zguess + zzero
  endif

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     call define_fields_79(itri, def_fields)

     if(eta0.ne.0) then

        ! for the grad-shafranov simulation with numvar < 3,
        ! calculate the pressure assuming that p(psi) = p0(psi)
        if(numvar.lt.3 .and. itor.eq.1) then
           temp79c = (pst79(:,OP_1) - psimin)/(psilim - psimin)
           
           do i=1,79
              if(temp79c(i).lt.0) then
                 pet79(i,OP_1) = p0-pi0*ipres
              else if(temp79c(i).gt.1) then
                 pet79(i,OP_1) = pedge*(p0-pi0*ipres)/p0
              else
                 pet79(i,OP_1) = pedge*(p0-pi0*ipres)/p0 + &
                      (p0-pi0*ipres)* &
                      (1.+p1*temp79c(i)+p2*temp79c(i)**2 &
                      -(20. + 10.*p1 + 4.*p2)*temp79c(i)**3 &
                      +(45. + 20.*p1 + 6.*p2)*temp79c(i)**4 &
                      -(36. + 15.*p1 + 4.*p2)*temp79c(i)**5 &
                      +(10. + 4.*p1 + p2)*temp79c(i)**6)
              endif
           end do
        endif


        if(itor.eq.1) then
           if(idens.eq.0) then
              temp79a = sqrt((1./(pefac*pet79(:,OP_1)))**3)
           else
              temp79a = sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
           endif
        else
           if(numvar.ge.3) then         
              if(idens.eq.0) then
                 temp79a = sqrt((1./(pefac*pet79(:,OP_1)))**3)
              else
                 temp79a = sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
              endif
           else
              if(idens.eq.0) then
                 temp79a = sqrt((1./(p0-pi0))**3)
              else
                 temp79a = sqrt((nt79(:,OP_1)/(p0-pi0))**3)
              endif
           endif
        endif
     else
        temp79a = 0.
     endif

     do i=1,18
        ione = isval1(itri,i)       

        resistivity(ione) = resistivity(ione) &
             + etar*int1(g79(:,OP_1,i),weight_79,79) &
             + eta0*int2(g79(:,OP_1,i),temp79a, weight_79,79)
     end do
  enddo

  ! solve linear equation
  call solve_newvar(resistivity, NV_NOBOUND)

!!$  if(idens.eq.0) dens=0.
!!$
!!$  ! Calculate RHS
!!$  call numnod(numnodes)
!!$  do i=1,numnodes
!!$
!!$     if(eta0.ne.0) then
!!$        call entdofs(1, i, 0, ibegin1, iendplusone1)
!!$        call entdofs(numvar, i, 0, ibegin, iendplusone)
!!$
!!$        if(idens.eq.0) den(ibegin1) = 1.
!!$
!!$        resistivity(ibegin1  ) = eta0*sqrt((den(ibegin1)/(pefac*phi(ibegin+12)))**3)
!!$        resistivity(ibegin1+1) = 1.5*resistivity(ibegin1)* &
!!$             (den(ibegin1+1)/den(ibegin1) - phi(ibegin+13)/phi(ibegin+12))
!!$        resistivity(ibegin1+2) = 1.5*resistivity(ibegin1)* &
!!$             (den(ibegin1+2)/den(ibegin1) - phi(ibegin+14)/phi(ibegin+12))
!!$        resistivity(ibegin1+3) = 1.5* &
!!$             ( resistivity(ibegin1+1)* &
!!$             (den(ibegin1+1)/den(ibegin1) - phi(ibegin+13)/phi(ibegin+12)) &
!!$             + resistivity(ibegin1)* &
!!$             (den(ibegin1+3)/den(ibegin1  ) - (den(ibegin1+1)/den(ibegin1  ))**2 & 
!!$             -phi(ibegin+15)/phi(ibegin+12) + (phi(ibegin+13)/phi(ibegin+12))**2))
!!$        resistivity(ibegin1+4) = &
!!$             resistivity(ibegin1+1)*resistivity(ibegin1+2)/resistivity(ibegin1) &
!!$             + 1.5*resistivity(ibegin1)* &
!!$             (den(ibegin1+4)/den(ibegin1  ) &
!!$             -den(ibegin1+1)*den(ibegin1+2)/(den(ibegin1  )**2) & 
!!$             -phi(ibegin+16)/phi(ibegin+12) &
!!$             +phi(ibegin+13)*phi(ibegin+14)/(phi(ibegin+12)**2))
!!$        resistivity(ibegin1+5) = 1.5* &
!!$             ( resistivity(ibegin1+2)* &
!!$             (den(ibegin1+2)/den(ibegin1) - phi(ibegin+14)/phi(ibegin+12)) &
!!$             + resistivity(ibegin1)* &
!!$             (den(ibegin1+5)/den(ibegin1  ) - (den(ibegin1+2)/den(ibegin1  ))**2 & 
!!$             -phi(ibegin+17)/phi(ibegin+12) + (phi(ibegin+14)/phi(ibegin+12))**2))
!!$     endif
!!$
!!$  enddo

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

  if(ibound.eq.NV_DCBOUND) then
     call boundary_dc(0,rhs)
     call solve(s6matrix_sm,rhs,ier)
  else
     call solve(s3matrix_sm,rhs,ier)
  endif

end subroutine solve_newvar

end module newvar_mod

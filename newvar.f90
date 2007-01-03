module newvar_mod

integer, parameter :: VAR_J   = 1
integer, parameter :: VAR_VOR = 2
integer, parameter :: VAR_COM = 3
integer, parameter :: VAR_SB1 = 4
integer, parameter :: VAR_SB2 = 5
integer, parameter :: VAR_SB3 = 6
integer, parameter :: VAR_SP1 = 7

contains

!============================================================
subroutine newvar(inarray,outarray,numvard,itype,iop,ibound)

  ! define a new variable with no boundary conditions applied
  !
  ! NOTE:  the operations for a given iop are defined by
  !             computed go-to below

  !     For iop=(1,3) this routine takes the Laplacian of the scalar field 
  !     in location iplace of the array inarray (with numvard components)  and
  !     puts the result in the stand-alone array outarray(mmnn6).  
  !     For iop=(4,7) ignores inarray, but defines new array in outarray.

  !     The LU decomposition of the mass-matrix takes place the first-time 
  !     called only.

  use p_data
  use t_data
  use basic
  use arrays
  use sparse

  implicit none
  
  integer, intent(in) :: numvard, itype, iop, ibound
  real, intent(in) :: inarray(*) ! length using numvard ordering
  real, intent(out) :: outarray(*) ! length using numvar=1 ordering

  integer :: i, j, ione, jone, itri, ier, ifirst
  integer :: numnodes, numelms, ndof
  integer :: izone, izonedim, nbound, nbcds
  real :: fintl(-6:maxi,-6:maxi)
  real :: dterm(18,18)
  real(r8), allocatable:: temp(:)

  integer, save, allocatable :: iboundds(:)

  !inarray sum matches first time through
  call createvec(temp, 1)
  select case(ibound)
  case(0)
     ifirst = ifirsts3_lu
  case(1)
     ifirst = ifirsts6_lu
  case default
     print *, 'Error: ibound not defined.', ibound
     return
  end select
  call numnod(numnodes)

  ! compute LU decomposition only once
  if(ifirst.eq.0) then

     call numfac(numelms)

     if(ibound.gt.0) then
        nbound = 0
        do i=1,numnodes
           call zonenod(i,izone,izonedim)
           if(izonedim .ne. 2) nbound = nbound + 1
        enddo
        nbound = 12*nbound    
     endif

     ! allocate matrix
     select case(ibound)
     case(0)
        call zeroarray4solve(s3matrix_sm, numvar1_numbering)
        ifirsts3_lu = 1
     case(1)
        call zeroarray4solve(s6matrix_sm, numvar1_numbering)
        ifirsts6_lu = 1
        allocate(iboundds(nbound))
     end select

     do itri=1,numelms
        call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
        call calcdterm(itri, dterm, fintl)
        do j=1,18
           jone = isval1(itri,j)
           do i=1,18
              ione = isval1(itri,i)
              select case(ibound)
              case(0)
                 call insertval(s3matrix_sm, dterm(i,j), ione, jone, 1)
              case(1)
                 call insertval(s6matrix_sm, dterm(i,j), ione, jone, 1)
              end select
           enddo
        enddo
     enddo

     select case(ibound)
     case(0)
        call finalizearray4solve(s3matrix_sm)
     case(1)
        call boundaryds(iboundds,nbcds,itype)
        do i=1,nbcds
           call setdiribc(s6matrix_sm, iboundds(i))
        enddo
        call finalizearray4solve(s6matrix_sm)
     end select
  endif

  temp = 0.

  ! define RHS temp
  select case(iop)
  case(VAR_J,VAR_VOR,VAR_COM)
     call newvar_gradshafranov(temp,inarray,numvard,itype)
              
  case(VAR_SB1)
     call newvar_SB1(temp)

  case(VAR_SB2)
     call newvar_SB2(temp)

  case(VAR_SP1)
     call newvar_SP1(temp)
              
  end select
  ! solve linear equation

  select case(ibound)
  case(0)
     call solve(s3matrix_sm,temp,ier)
  case(1)
     do i=1,nbcds
        temp(iboundds(i)) = 0.
     enddo
     call solve(s6matrix_sm,temp,ier)
  end select

  if(iprint.ge.1) write(*,*) "newvar: after solve", myrank
  call numdofs(1, ndof)
  do i=1,ndof
     outarray(i) = temp(i)
  enddo
  
  call deletevec(temp)

end subroutine newvar



subroutine newvar_gradshafranov(temp,inarr,numvard,itype)

  use basic
  use t_data
  use arrays
  use nintegrate_mod
  
  implicit none

  integer, intent(in) :: numvard, itype
  real, intent(in) :: inarr(*)
  real, intent(out) :: temp(*)

  integer :: itri, i, j, ione, j1, numelms
  real :: sum

  call numfac(numelms)

  do itri=1,numelms

     ! calculate the local sampling points and weights for numerical integration
     call area_to_local(25,                                            &
          alpha_25,beta_25,gamma_25,area_weight_25,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_25, eta_25, weight_25)
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), &
             ri_79, 79, g79(:,:,i))
     end do
     do i=1,18
        ione = isval1(itri,i)
        sum = 0.
        do j=1,18
           j1 = isvaln(itri,j)
           sum = sum - inarr(j1 + 6*(itype-1)) * &
                (int2(g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79) &
                +int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79))
           if(itor.eq.1) then
              sum = sum - inarr(j1 + 6*(itype-1)) * &
                   2.*int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),weight_79,79)
           endif
        end do
        temp(ione) = temp(ione) + sum
     end do
  enddo
  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumshareddofs so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values
  call sumshareddofs(temp)
end subroutine newvar_gradshafranov

subroutine newvar_SB1(temp)

  use basic
  use t_data
  use arrays
  use nintegrate_mod
  use metricterms_n

  implicit none

  real, intent(out) :: temp(*)

  integer :: itri, i, ione, numelms
  real :: sum, factor, x, z, xmin, zmin, avec(20), dbf, deex, hypf

  double precision :: cogcoords(3)

  call numfac(numelms)
  call getmincoord(xmin,zmin)

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
     hypf = hyper * deex**2
     
     ! calculate the local sampling points and weights for numerical integration
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     ! define r, 1/r
     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79

     ! define fields
     call calcavector(itri, phi, 1, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps179)
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, phi0, 1, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps079)
        pst79 = ps079 + ps179
     else
        pst79 = ps179
     endif

     if(numvar.ge.2) then
        call calcavector(itri, phi, 2, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz179)
     
        if(linear.eq.1 .or. eqsubtract.eq.1) then
           call calcavector(itri, phi0, 2, numvar, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz079)
           bzt79 = bz079 + bz179
        else
           bzt79 = bz179
        endif

        if(idens.eq.1) then
           call calcavector(itri, deni, 1, 1, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ni79)
        endif
     endif

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), &
             ri_79, 79, g79(:,:,i))
     end do

     ! Evaluate matrix elements
     ! ~~~~~~~~~~~~~~~~~~~~~~~~
     do i=1,18
        ione = isval1(itri,i)

        sum = b1psieta(g79(:,:,i),pst79,hypf)*etar
        if(numvar.ge.2) then
           sum = sum + b1psibd(g79(:,:,i),pst79,bzt79,ni79)
        endif

        temp(ione) = temp(ione) + sum
     end do
  enddo
  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumshareddofs so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values
  call sumshareddofs(temp)
end subroutine newvar_SB1

subroutine newvar_SB2(temp)

  use basic
  use t_data
  use arrays
  use nintegrate_mod
  use metricterms_n
  
  implicit none

  real, intent(out) :: temp(*)

  integer :: itri, i, ione, j, j1, j01, numelms
  real :: sum, factor, x, z, xmin, zmin, avec(20), dbf, deex, hypi

  double precision :: cogcoords(3)

  call numfac(numelms)
  call getmincoord(xmin,zmin)

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
     hypi = hyperi* deex**2

     ! calculate the local field values
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79
     ri2_79 = ri_79*ri_79
     ri3_79 = ri2_79*ri_79
     ri4_79 = ri2_79*ri2_79
  
     call calcavector(itri, phi, 1, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps179)

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, phi0, 1, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps079)
        pst79 = ps079 + ps179
     else
        pst79 = ps179
     endif
     
     if(numvar.ge.2) then
        call calcavector(itri, phi, 2, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz179)
        
        if(linear.eq.1 .or. eqsubtract.eq.1) then
           call calcavector(itri, phi0, 2, numvar, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz079)
           bzt79 = bz079 + bz179
        else
           bzt79 = bz179
        endif

        if(numvar.ge.3) then
           call calcavector(itri, phi, 3, numvar, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, pe179)
        
           if(linear.eq.1 .or. eqsubtract.eq.1) then
              call calcavector(itri, phi0, 3, numvar, avec)
              call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, pe079)
              pet79 = pe079 + pe179
           else
              pet79 = pe179
           endif
        endif
     endif
     
     if(idens.eq.1) then
        call calcavector(itri, deni, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ni79)
     endif

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
     end do

     ! Evaluate matrix elements
     ! ~~~~~~~~~~~~~~~~~~~~~~~~
     do i=1,18
        ione = isval1(itri,i)

        sum = b2psipsid(g79(:,:,i),pst79,pst79,ni79) &
             +b2bbd    (g79(:,:,i),bzt79,bzt79,ni79) &
             +b2beta   (g79(:,:,i),bzt79,hypi)*etar

        if(numvar.ge.3) then
           sum = sum + b2ped(g79(:,:,i),pet79,ni79)
        endif

        temp(ione) = temp(ione) + sum
     end do
  enddo
  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumshareddofs so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values
  call sumshareddofs(temp)
end subroutine newvar_SB2



! Pressure source term
! ====================
subroutine newvar_SP1(temp)

  use basic
  use t_data
  use arrays
  use nintegrate_mod
  use metricterms_n
  
  implicit none

  real, intent(out) :: temp(*)

  integer :: itri, i, ione, j, j1, j01, numelms
  real :: sum, factor, x, z, xmin, zmin, avec(20), dbf, deex, hypi

  double precision :: cogcoords(3)

  call numfac(numelms)
  call getmincoord(xmin,zmin)

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
     hypi = hyperi* deex**2

     ! calculate the local field values
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79
     ri2_79 = ri_79*ri_79
     ri3_79 = ri2_79*ri_79
     ri4_79 = ri2_79*ri2_79
  
     call calcavector(itri, phi, 1, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps179)

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, phi0, 1, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ps079)
        pst79 = ps079 + ps179
     else
        pst79 = ps179
     endif
     
     call calcavector(itri, phi, 2, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz179)
        
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, phi0, 2, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, bz079)
        bzt79 = bz079 + bz179
     else
        bzt79 = bz179
     endif

     call calcavector(itri, phi, 3, numvar, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, pe179)
        
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calcavector(itri, phi0, 3, numvar, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, pe079)
        pet79 = pe079 + pe179
     else
        pet79 = pe179
     endif

     if(ipres.eq.1) then
        call calcavector(itri, pres, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, p179)
        
        if(linear.eq.1 .or. eqsubtract.eq.1) then
           call calcavector(itri, pres0, 1, 1, avec)
           call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, p079)
           pt79 = p079 + p179
        else
           pt79 = p179
        endif
     else
        pt79 = pet79        
     endif

     if(idens.eq.1) then
        call calcavector(itri, deni, 1, 1, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, ni79)
     endif

     b2i79(:,OP_1) = 1./(pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2 + bzt79(:,OP_1)**2)

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
     end do

     ! Evaluate matrix elements
     ! ~~~~~~~~~~~~~~~~~~~~~~~~
     do i=1,18
        ione = isval1(itri,i)

        sum = b3pebd     (g79(:,:,i),pet79,bzt79,ni79)   &
             +b3psipsieta(g79(:,:,i),pst79,pst79)*etar   &
             +b3bbeta    (g79(:,:,i),bzt79,bzt79)*etar   &
             +b3pedkappa (g79(:,:,i),pt79, ni79 )*kappat &
             +p1kappar   (g79(:,:,i),pst79,pst79,pt79,ni79,b2i79)*kappar

        temp(ione) = temp(ione) + sum
     end do
  enddo
  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumshareddofs so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values
  call sumshareddofs(temp)
end subroutine newvar_SP1

end module newvar_mod

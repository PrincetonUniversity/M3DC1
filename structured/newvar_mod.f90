module newvar_mod

implicit none

! newvarb
! newvar
! newvarbv


contains

subroutine newvarb(inarray,outarray,mmnn,numvard,iplace,iop)
!....define a new variable with homogeneous boundary conditions applied
!
!      NOTE:  the operations for a given iop are defined by
!             computed go-to below
!
!.....replaces delsquared (05/04/05)...defines a new variable
!     iop=1  defines J = delsquared(psi)
!     iop=2  defines vor = delsquared(phi)
!     iop=3  defines com = delsquared(chi)
!     iop=4  defines sphi1..source term in poloidal flux equation
!     iop=5  defines sphi2..source term in toroidal field equation
!     iop=6  defines sphiE..source term in electron pressure equation
!     iop=7  defines sphiP..source term in pressure equation
!     iop=8-13 different components of the toroidal electric field
!     iop=15 defines tw1
!     iop=16 defines tw2
!     iop=17 defines tw3
!     iop=18 defines tw4 (inarray=pres+pres0) or tw40 (inarray=pres0)
!     iop=19 defines tw5 (inarray=pres+pres0) or tw50 (inarray=pres0)
!
!     For iop=(1,3) this routine takes the Laplacian of the scalar field 
!     in location iplace of the array inarray (with numvard components)  and
!     puts the result in the stand-alone array outarray(mmnn6).  
!     For iop=(4,13) ignores inarray, but defines new array in outarray.

  use p_data
  use t_data
  use basic
  use arrays
  use superlu
#ifdef mpi
  use supralu_dist_mod
  implicit none
  include 'mpif.h'
#else
  implicit none
#endif
  integer :: mmnn, numvard, iplace, iop
  real ::  inarray(mmnn*6*numvard), outarray(mmnn*6)

  integer, parameter :: r8a = selected_real_kind(12,100)
  real(r8a), allocatable :: phin(:),temp(:)
  integer, save, allocatable :: iboundds(:)
  real :: temparr(18,18)
  real :: denf, hypp, termbf
  real :: ssterm, x, z, dbf, sum, factor, gmix, terma, termb
  integer ::  mmnn6, lx, lz, nrads, numvards, msizeds, iodd, filenum
  integer :: ier, jer, itri, irect, jrect, ll
  integer :: ir1, ir2, ir3, nbcds, jsymtype
  integer :: i, ione, i1, i2, i3
  integer :: j, jone, j01, j1, j2, j3
  integer :: k, kone, k01, k1, k2, k3
  integer :: l, lone, l1, l2, l3
  integer :: m_lds, nnz_lds, ifs_rds
  real :: tbeforesolve, taftersolve
  real :: tstart, tend


  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
      if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "newvarb called with iop=",iop

  mmnn6 = 6*m*n
!  allocate(phin(mmnn6),temp(mmnn6),temptemp(mmnn6))
  allocate(phin(mmnn6))

  ! put the scalar in inarray into a numvar=1 storage array
  do lx=1,n
     do lz=1,m
        l = lx + (lz-1)*n
        do i=1,6
           phin((l-1)*6+i) = inarray(6*((l-1)*numvard+(iplace-1))+i)
        enddo
     enddo
  enddo

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tnewvarb = tnewvarb + tend - tstart
  endif

  if(ifirsts6_lu.ne.0) go to 500
  ifirsts6_lu = 1
  allocate(iboundds(12*(m+n-2)))
  
  nrads = n*m*6
  numvards = 1
  msizeds = 9*n*m*(6*numvards)**2
  call jbdecomp1(n*m, 6, m_lds, nnz_lds, ifs_rds)
  base = 0

  ! compute LU decomposition only once
  
  ! form s-matrix matrices
  ss = 0
  itri = 0
  do iodd=1,ioddmx
     do ll=1,nreg
        itri = itri + 1
        do j=1,18
           do i=1,18
              jone = isval1(itri,j)
              ione = isval1(itri,i)
              ssterm = dterm(iodd,i,j)
              call increment(ss,m*n,ione,jone,ssterm,1,8)
           enddo
        enddo
     enddo
  enddo

  ! define indices for boundary arrays
  jsymtype = 1              !  (even up-down symmetry)
  if(iop.eq.2) jsymtype = 2
  call boundaryds(iboundds,nbcds,jsymtype)
  
  ! modify the s-matrix, inserting the boundary conditions
  do l=1,nbcds
     i = iboundds(l)
     do j= 1,nrads
        if(i.le.0 .or. i.ge.m*n*6*numvards+1) then
           write(*,3343) i,l,iboundds(l)
3343       format(" error in iboundds index",3i5)
           call safestop(3343)
        endif
        call assign(ss,n*m,i,j,0.,numvards)
     enddo
  enddo
  do l=1,nbcds
     call assign(ss,n*m,iboundds(l),iboundds(l),1.,numvards)
  enddo
  ! define the superlu matrices "s6matrix_lu" and "samatrix_lu"
#ifdef mpi
      call sparseR8d_init(s6matrix_lu,nrads,nnz_lds,m_lds,ifs_rds,ier)
#else
  call sparseR8_init(s6matrix_lu,nrads,msizeds,base,ier)
  s6matrix_lu %permc_spec = 0
#endif
  call define(s6matrix_lu,ss,m*n,1)

  ! perform LU decomposition of sparse matrix "s6matrix"
  ! store result in opaque object "s6handle"
  jer = 0
#ifdef mpi
  call sparseR8d_new(s6matrix_lu, jer)
  if(jer.ne.0) then
     write(*,*) 'after 3rd sparseR8d_new', jer
     call safestop(46)
  endif
#else
  call dsupralu_new(s6handle, s6matrix_lu%values(1),               &
       s6matrix_lu%irow(1), s6matrix_lu%jcol_ptr(1),               &
       s6matrix_lu%nnz, nrads , jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_new', jer
     call safestop(47)
  endif
  call dsupralu_colperm(s6handle, s6matrix_lu% permc_spec, jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_new', jer
     call safestop(48)
  endif
  call dsupralu_lu(s6handle, jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_lu', jer
     call safestop(49)
  endif
#endif

500 continue

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  allocate(temp(s6matrix_lu%m_loc))

  temp = 0
  filenum = 0

  select case(iop)
! case(7)
!    filenum = 105
! case(11)
!    filenum = 151
! case(12)
!    filenum = 107
! case(13)
!    filenum = 152
! case(14)
!    filenum = 153
!    filenum = 153
  case(15) 
     filenum = 101
  case(16) 
     filenum = 178
  case(17) 
     filenum = 179
  case(18) 
     filenum = 178
  case(19) 
     filenum = 179
  case default 
     filenum = 0
  end select

  ! position files with metric terms
  if(filenum.ne.0) rewind(filenum)
  if(iop.eq.6) then
     if(iread(09).ne.0) rewind(109)
     if(iread(10).ne.0) rewind(110)
  endif
  if(iop.ge.8 .and.iop.le.13) then
     if(iread(01).ne.0) rewind(101)
     if(iread(02).ne.0) rewind(102)
  endif

  ! define RHS temp
  do iodd=1,ioddmx
  do i=1,18
  do j=1,18

     terma = aterm(iodd,i,j)
     termb = bterm(iodd,i,j)
     hypp = hyperp*deex**2*termb

     if(filenum.ne.0) then
        read(filenum) ir1,ir2,ir3,((temparr(k,l),l=1,18),k=1,18)
     endif
     if(iop.eq.6) then
        if(iread(09).eq.1)                                             &
             read(109) ir1,ir2,ir3,((u5terml(k,l),l=1,18),k=1,18)
        if(iread(10).eq.1)                                             &
             read(110) ir1,ir2,ir3,((u6terml(k,l),l=1,18),k=1,18)
     endif
     if(iop.ge.8 .and. iop.le.13) then
        if(iread(01).ne.0) read(101) ir1,ir2,ir3,((u3terml(k,l),l=1,18),k=1,18)
        if(iread(02).ne.0) read(102) ir1,ir2,ir3,((u4terml(k,l),l=1,18),k=1,18)
     endif

     do jrect=1,m-1+jper
     do irect=1,n-1+iper

        ll = irect + (jrect-1)*(n-1+iper)
        itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
        ione = isval1(itri,i)

        if(sparseR8d_is_local_row(s6matrix_lu, ione).eq.0) goto 800

        i1 = isvaln(itri,1)
        i2 = i1 + 6
        i3 = i2 + 6
        jone = isval1(itri,j)
        j1 = isvaln(itri,j)
        j01 = isval0(itri,j)
        j2 = j1 + 6
        j3 = j2 + 6

        x = (irect-iodd/3.)*deex
        z = (jrect-1+iodd/3.)*deez
        factor = 1
        if(imask.eq.1) call mask(x,z,factor)
        dbf = db*factor

        sum = 0. 

        ! add in ohmic and viscous heating terms calculated in ludefphi
        if(iop.eq.6) sum = sum + ohmic(i3)   
        if(iop.eq.7) sum = sum + ohmic(i3) + viscous(i3)
        
        go to (101,101,101,104,105, 106,107,108,109,110, &
               111,112,113,800,115, 115,115,118,118),iop
        ! delsquared
101     sum = sum + aterm(iodd,i,j)*phin(jone)
        go to 200

        ! RHS of poloidal flux equation psi
104     continue
        sum = sum + etar*(aterm(iodd,i,j)                              &
             -hyper*deex**2*bterm(iodd,i,j))*phi(j1)
        if(linear.eq.0) then
           sum = sum + etar*(aterm(iodd,i,j)                           &
                -hyper*deex**2*bterm(iodd,i,j))*phi0(j01)
        endif
        if(numvar.ge.2) then
           do k=1,18
              kone = isval1(itri,k)
              k1= isvaln(itri,k)
              k2= k1 + 6
     
              ! NOTE: this may need to be modified to add back in the
              ! equilibrium part for linear=0
              sum = sum + dbf*k0term(iodd,i,j,k)                       &
                   *(phi(k2)*(0.5*phi(j1)+phi0(j01))                   &
                   + phi(j1)*(0.5*phi(k2)+phi0(k2)))
           enddo               ! on k
        endif                  ! numvar.ge.2
        go to 200
     
        ! RHS of toroidal field equation I
105     continue
        sum = sum + etar*(aterm(iodd,i,j)                              &
             -hyperi*deex**2*bterm(iodd,i,j))*phi(j2)
        if(linear.eq.0) then
           sum = sum + etar*(aterm(iodd,i,j)                            &
                -hyperi*deex**2*bterm(iodd,i,j))*phi0(j2)
        endif
        do k=1,18
           kone = isval1(itri,k)
           k1= isvaln(itri,k)
           k2= k1 + 6
           k01 = isval0(itri,k)
           
           ! NOTE: this may need to be modified to add back in the equilibrium
           !       part for linear=0
           sum = sum + dbf*(g0term(iodd,i,j,k)+g0term(iodd,i,k,j))     &
                *phi(j1)*(0.5*phi(k1)+phi0(k01))
        enddo                  ! on k
        go to 200
        
        ! RHS of electron pressure equation
106     continue
        
        sum = sum   + kappa*(terma-hyperp*deex**2*termb)*phi(j3)
        if(linear.eq.0) then
           sum = sum + kappa*(terma-hyperp*deex**2*termb)*phi0(j3)
        endif
        
        do k=1,18
           kone = isval1(itri,k)
           k1 = isvaln(itri,k)
           k2 = k1 + 6
           k3 = k2 + 6
           gmix = g2term(iodd,i,k,j)+g2term(iodd,i,j,k)
!>> Changed (2./3.) to (gam-1.) NMF
           sum = sum + (gam-1.)*kappat*(gmix                            &
                - hypp*g13erm(iodd,i,k,j))*deni(kone)*phi(j3)
           if(idens.eq.0) then
              sum = sum + dbf*k0term(iodd,i,k,j)                       &
                   *(phi(k3)+phi0(k3))*(phi(j2))
              sum = sum + dbf*k0term(iodd,i,j,k)                       &
                   *(phi(k2)+phi0(k2))*phi(j3)
           endif
           
           
           ! QUADRATIC NONLINEAR TERMS
           do l=1,18
              l1 = isvaln(itri,l)
              l2 = l1 + 6
              l3 = l2 + 6
              lone = isval1(itri,l)
              if(idens.gt.0) then
                 denf = deni(lone)*dbf
                 sum = sum + denf*u6terml(k,l)*(phi(k3)+phi0(k3))*phi(j2)
                 sum = sum + denf*u5terml(k,l)*(phi(k2)+phi0(k2))*phi(j3)
              endif
              
           enddo               ! on l
        enddo                  ! on k
        go to 200
        
        ! RHS of pressure equation
107     continue
        
        do k=1,18
           kone = isval1(itri,k)
           k1= isvaln(itri,k)
           k2= k1 + 6
           do l=1,18
              lone = isval1(itri,l)
              l1= isvaln(itri,l)
              l2= l1 + 6
              l3= l2 + 6
              
!>> Changed -gam to +gam NMF
              sum = sum + temparr(k,l)*(phi(k2)+phi0(k2))              &
                   *(     (phi(j3)+phi0(j3))*(deni(lone))              &
                     +gam*(phi(l3)+phi0(l3))*(deni(jone)))
           enddo               ! on l
        enddo                  ! on k
        go to 200
!
!...108-113 added 06/04/06 to check toroidal electric field as diagnostic (SCJ)
108      continue
!
           terma = aterm(iodd,i,j)
           termb = bterm(iodd,i,j)
           termbf = hyper*deex**2*termb
           dbf = db*factor
!
! LINEAR TERMS
           sum = sum + etar*(terma-termbf)*phi(j1)   
!
! SIMPLE NONLINEAR TERMS
           do k=1,18
             k1= isvaln(itri,k)
             k01 = isval0(itri,k)
             k2= k1 + 6
             k3= k2 + 6
             kone = isval1(itri,k)
!
             sum = sum + (k0term(iodd,i,k,j)*(vel(k1)+vel0(k1)))*(phi(j1)+phi0(j01))
!
             if(numvar.ge.2) then
               if(idens.eq.0) then
!
                 sum = sum + dbf*k0term(iodd,i,j,k)*thimp*(phi(k2)+phi0(k2))      &
                                                         *(phi(j1)+phi0(j01))
!
               endif ! on idens.eq.0
               if(numvar.ge.3) then
!
               sum = sum + g2term(iodd,j,k,i)*(vel(k3)+vel0(k3))*(phi(j1)+phi0(j01))
!
               endif  ! on numvar.ge.3
!
            endif  ! on numvar.ge.2
!
! QUADRATIC NONLINEAR TERMS
            do l=1,18
              l1 = isvaln(itri,l)
              l2 = l1 + 6
              l3 = l2 + 6
              lone = isval1(itri,l)
              if(idens.gt.0) then
              if(numvar.ge.2) then
                denf = deni(lone)*dbf
                sum = sum + u3terml(k,l)*denf*(phi(k2)+phi0(k2))*(phi(j1)+phi0(j01))
              endif ! on numvar.gt.2
              endif ! on idens.gt.0
!
            enddo ! on l
          enddo ! on k
         go to 200
!
!........V x B term (incompressible)
109      continue
!
           do k=1,18
             k1= isvaln(itri,k)
!
             sum = sum + (k0term(iodd,i,k,j)*(vel(k1)+vel0(k1)))               &
                                            *(phi(j1)+phi0(j01))
!
           enddo ! on k
         go to 200
!
!........V x B term (compressible)
110      continue 
!
           do k=1,18
             k1= isvaln(itri,k)
             k01 = isval0(itri,k)
             k2= k1 + 6
             k3= k2 + 6
             if(numvar.ge.2) then
               if(numvar.ge.3) then
!
               sum = sum + g2term(iodd,j,k,i)*(vel(k3)+vel0(k3))               &
                                             *(phi(j1)+phi0(j01))
!
               endif  ! on numvar.ge.3
             endif  ! on numvar.ge.2
           enddo ! on k
         go to 200
!
!........resistivity term
111      continue
           terma = aterm(iodd,i,j)
           sum = sum + etar*(terma)*phi(j1)   
         go to 200
!
!.......J x B term
112      continue
           dbf = db*factor
           do k=1,18
             k1= isvaln(itri,k)
             k01 = isval0(itri,k)
             k2= k1 + 6

            do l=1,18
              l1 = isvaln(itri,l)
              l2 = l1 + 6
              l3 = l2 + 6
              lone = isval1(itri,l)
              if(idens.gt.0) then
              if(numvar.ge.2) then
                denf = deni(lone)*dbf
                sum = sum + u3terml(k,l)*denf*(phi(k2)+phi0(k2))              &
                                             *(phi(j1)+phi0(j01))
              endif ! on numvar.gt.2
              endif ! on idens.gt.0
!
            enddo ! on l
          enddo ! on k
         go to 200
!
!........hyper term
113      continue
           termb = bterm(iodd,i,j)
           termbf = hyper*deex**2*termb
           sum = sum + etar*(-termbf)*phi(j1)   
!
         go to 200

! tw1, tw2, tw3     
115     continue
        
        do k=1,18
           kone = isval1(itri,k)

           if(idens.eq.1) then
              do l=1,18
                 lone = isval1(itri,l)
                 
                 sum = sum + temparr(k,l) &
                      *(phi(j1)+phi0(j01))*deni(kone)*bsqri(lone)
              enddo               ! on l
           else !on idens.eq.1
              select case(iop)
                 ! tw1
              case(15)
                 sum = 0.
                 
                 !tw2
              case(16)
                 sum = sum + y0term(iodd,i,j,k)* &
                      (phi(j1)+phi0(j01))*bsqri(kone)
                 
                 !tw3
              case(17)
                 sum = sum + x0term(iodd,i,j,k)* &
                      (phi(j1)+phi0(j01))*bsqri(kone)
                 
              end select
           endif !on idens.eq.1

        enddo                  ! on k
        go to 200

! tw4, tw5     
118     continue
        
        do k=1,18
           kone = isval1(itri,k)

           if(idens.eq.1) then
              do l=1,18
                 lone = isval1(itri,l)
                 
                 sum = sum + temparr(k,l)*bsqri(lone)* &
                      (phin(jone)*deni(kone) + phin(kone)*deni(jone))
              enddo               ! on l
           else !on idens.eq.1
              select case(iop)
                 !tw4
              case(18)
                 sum = sum + y0term(iodd,i,j,k)*phin(jone)*bsqri(kone)
                 
                 !tw5
              case(19)
                 sum = sum + x0term(iodd,i,j,k)*phin(jone)*bsqri(kone)
                 
              end select
           endif !on idens.eq.1

        enddo                  ! on k
        go to 200
     
200     temp(ione-s6matrix_lu%fst_row) = temp(ione-s6matrix_lu%fst_row) + sum
800     continue
     enddo                     ! on irect
     enddo                     ! on jrect
  enddo                     ! on j
  enddo                     ! on i
  enddo                     ! on iodd

  ! complete the imposition of the boundary conditions,
  do l=1,nbcds
     if(sparseR8d_is_local_row(s6matrix_lu, iboundds(l)).eq.1) then
        temp(iboundds(l)-s6matrix_lu%fst_row) = 0.
     endif
  enddo
  if(myrank.eq.0) call second(tbeforesolve)

  ! perform LU backsubstitution to get outarray solution
  call sparseR8d_solve_part(s6matrix_lu,temp,outarray,ier)

  if(myrank.eq.0) call second(taftersolve)
  tsolve = tsolve + taftersolve-tbeforesolve
  if(myrank.eq.0 .and. iprint.ge.1) write(*,*) ier, "after solve1"

  
  deallocate(phin,temp)


  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tnewvarb = tnewvarb + tend - tstart
  endif

  return
  
end subroutine newvarb


!============================================================
subroutine newvar(inarray,outarray,mmnn,numvard,iplace,iop)
!
!....define a new variable with no boundary conditions applied
!
!      NOTE:  the operations for a given iop are defined by
!             computed go-to below
!
!.....replaces delsquared (05/04/05)...defines a new variable
!     iop=1  defines J = delsquared(psi)
!     iop=2  defines vor = delsquared(phi)
!     iop=3  defines com = delsquared(chi)
!     iop=4  defines sphi1..source term in poloidal flux equation
!     iop=5  defines sphi2..source term in toroidal field equation
!     iop=6  defines sphiE..source term in electron pressure equation
!     iop=7  defines sphiP..source term in pressure equation
!     iop=8  defines bsqr
!     iop=9  defines galph
!     iop=10 defines ggam
!     iop=11 defines gmu
!     iop=12 defines gw1
!     iop=13 defines gw2
!     iop=14 defines gw3
!     iop=15 defines tw1
!     iop=16 defines tw2
!     iop=17 defines tw3
!     iop=18 defines tw4 (inarray=pres+pres0) or tw40 (inarray=pres0)
!     iop=19 defines tw5 (inarray=pres+pres0) or tw50 (inarray=pres0)
!
!     For iop=(1,3) this routine takes the Laplacian of the scalar field 
!     in location iplace of the array inarray (with numvard components)  and
!     puts the result in the stand-alone array outarray(mmnn6).  
!     For iop=(4,7) ignores inarray, but defines new array in outarray.

!     The LU decomposition of the mass-matrix takes place the first-time 
!     called only.
!
!
  use p_data
  use t_data
  use basic
  use arrays
  use superlu
#ifdef mpi
  use supralu_dist_mod
  use superlu_mod
  implicit none
  include 'mpif.h'
#else 
  implicit none
#endif
  integer mmnn, numvard, iplace, iop
  real  inarray(mmnn*6*numvard),outarray(mmnn*6)
  integer, parameter :: r8a = selected_real_kind(12,100)
  real(r8a), allocatable:: phin(:),temp(:)
  integer, save, allocatable:: iboundds(:)
  real temparr(18,18)
  integer filenum
  real ssterm, sum, terma, termb, hypp, denf, gmix, factor, x, z, dbf
  integer mmnn6, lx, lz, nrads, numvards, msizeds, iodd
  integer ier, jer, itri, irect, jrect, ll
  integer ir1, ir2, ir3
  integer i, ione, i1, i2, i3
  integer j, jone, j01, j1, j2, j3
  integer k, kone, k01, k1, k2, k3
  integer l, lone, l1, l2, l3
  integer m_lds, nnz_lds, ifs_rds
  real tbeforesolve, taftersolve, tstart, tend

  mmnn6 = 6*m*n
!  allocate(phin(mmnn6),temp(mmnn6),temptemp(mmnn6))
  allocate(phin(mmnn6))

  ! put the scalar in inarray into a numvar=1 storage array
  do lx=1,n
     do lz=1,m
        l = lx + (lz-1)*n
        do i=1,6
           phin((l-1)*6+i) = inarray(6*((l-1)*numvard+(iplace-1))+i)
        enddo
     enddo
  enddo

  if(ifirsts3_lu.ne.0) go to 500
  ifirsts3_lu = 1
  allocate(iboundds(12*(m+n-2)))
  
  nrads = n*m*6
  numvards = 1
  msizeds = 9*n*m*(6*numvards)**2
  call jbdecomp1(n*m, 6, m_lds, nnz_lds, ifs_rds)
  base = 0
  
  ! compute LU decomposition only once
  
  ! form s-matrix matrices
  ss = 0
  itri = 0
  
  do iodd=1,ioddmx
     do ll=1,nreg
        itri = itri + 1
        do j=1,18
           do i=1,18
              jone = isval1(itri,j)
              ione = isval1(itri,i)
              ssterm = dterm(iodd,i,j)
              call increment(ss,m*n,ione,jone,ssterm,1,8)
              
           enddo
        enddo
     enddo
  enddo

#ifdef mpi
  call sparseR8d_init(s3matrix_lu,nrads,nnz_lds,m_lds,ifs_rds,ier)
#else
  call sparseR8_init(s3matrix_lu,nrads,msizeds,base,ier)
  s3matrix_lu %permc_spec = 0
#endif
  call define(s3matrix_lu,ss,m*n,1)
  
  ! perform LU decomposition of sparse matrix "s3matrix"
  ! store result in opaque object "s3handle"
  jer = 0
#ifdef mpi
  call sparseR8d_new(s3matrix_lu, jer)
  if(jer.ne.0) then
     write(*,*) 'after 3rd sparseR8d_new', jer
     call safestop(46)
  endif
#else
  call dsupralu_new(s3handle, s3matrix_lu%values(1),               &
       s3matrix_lu%irow(1), s3matrix_lu%jcol_ptr(1),               &
       s3matrix_lu%nnz, nrads , jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_new', jer
     call safestop(47)
  endif
  call dsupralu_colperm(s3handle, s3matrix_lu% permc_spec, jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_new', jer
     call safestop(48)
  endif
  call dsupralu_lu(s3handle, jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_lu', jer
     call safestop(49)
  endif
#endif

500 continue
  
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  allocate(temp(s3matrix_lu%m_loc))

  temp = 0

  select case(iop)
  case(7) 
     filenum = 105
  case(11) 
     filenum = 151
  case(12) 
     filenum = 107
  case(13) 
     filenum = 152
  case(14) 
     filenum = 153
  case(15) 
     filenum = 101
  case(16) 
     filenum = 178
  case(17) 
     filenum = 179
  case(18) 
     filenum = 178
  case(19) 
     filenum = 179
  case default 
     filenum = 0
  end select
  
  ! position files with metric terms
  if(filenum.ne.0) rewind(filenum)
  if(iop.eq.6) then
     rewind(109)
     rewind(110)
  endif
  
  do iodd=1,ioddmx
  do i=1,18
  do j=1,18

     terma = aterm(iodd,i,j)
     termb = bterm(iodd,i,j)
     hypp = hyperp*deex**2*termb
          
     if(filenum.ne.0) then
        read(filenum) ir1,ir2,ir3,((temparr(k,l),l=1,18),k=1,18)
     endif
     if(iop.eq.6) then
        read(109) ir1,ir2,ir3,((u5terml(k,l),l=1,18),k=1,18)
        read(110) ir1,ir2,ir3,((u6terml(k,l),l=1,18),k=1,18)
     endif

     do jrect=1,m-1+jper
     do irect=1,n-1+iper
     
        ll = irect + (jrect-1)*(n-1+iper)
        itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
        ione = isval1(itri,i)

        if(sparseR8d_is_local_row(s3matrix_lu, ione).eq.0) goto 800

        i2 = i1 + 6
        i3 = i2 + 6
        jone = isval1(itri,j)
        j01 = isval0(itri,j)
        j1 = isvaln(itri,j)
        j2 = j1 + 6
        j3 = j2 + 6

        x = (irect-iodd/3.)*deex
        z = (jrect-1+iodd/3.)*deez
        factor = 1
        if(imask.eq.1) call mask(x,z,factor)
        dbf = db*factor

        sum = 0.

        go to (111,111,111,114,200, 116,117,140,140,140,                &
               140,140,140,140,125, 125,125,128,128),iop
        write(*,*) "Error: iop not defined for newvar: ", iop

111     sum = sum + aterm(iodd,i,j)*phin(jone)
        go to 200
114     sum = sum + etar*aterm(iodd,i,j)*phi(j1)
        go to 200


        ! iop = 6: RHS of electron pressure equation
        ! ------------------------------------------
116     continue

        sum = sum + ohmic(i3)        
        sum = sum + kappa*(terma-hypp)*phi(j3)
        if(linear.eq.0) then
           sum = sum + kappa*(terma-hypp)*phi0(j3)
        endif
        
        do k=1,18
           kone = isval1(itri,k)
           k1 = isvaln(itri,k)
           k2 = k1 + 6
           k3 = k2 + 6
           gmix = g2term(iodd,i,k,j)+g2term(iodd,i,j,k)
           sum = sum + (gam-1.)*kappat*(gmix                       &
                - hypp*g13erm(iodd,i,k,j))*deni(kone)*phi(j3)
           if(idens.eq.0) then
              sum = sum + dbf*k0term(iodd,i,k,j)*(phi(k3)+phi0(k3))*phi(j2)
              sum = sum + dbf*k0term(iodd,i,j,k)*(phi(k2)+phi0(k2))*phi(j3)
           else
              do l=1,18
                 lone = isval1(itri,l)

                 denf = deni(lone)*dbf
                 sum = sum + denf*u6terml(k,l)*(phi(k3)+phi0(k3))*phi(j2)
                 sum = sum + denf*u5terml(k,l)*(phi(k2)+phi0(k2))*phi(j3)
              enddo               ! on l
           endif ! on idens.eq.0
        enddo                  ! on k
        go to 200

        ! iop = 7: RHS of pressure equation
        ! ---------------------------------
117     continue
        
        sum = sum + ohmic(i3) + viscous(i3)

        do k=1,18
           kone = isval1(itri,k)
           k1= isvaln(itri,k)
           k2= k1 + 6
           do l=1,18
              lone = isval1(itri,l)
              l1= isvaln(itri,l)
              l2= l1 + 6
              l3= l2 + 6
              
!>> NMF changed -gam to +gam 
              sum = sum + temparr(k,l)*(phi(k2)+phi0(k2))              &
                   *(     (phi(j3)+phi0(j3))*(deni(lone))              &
                     +gam*(phi(l3)+phi0(l3))*(deni(jone)))
           enddo               ! on l
        enddo                  ! on k
        go to 200
   
! tw1, tw2, tw3     
125     continue
        
        do k=1,18
           kone = isval1(itri,k)

           if(idens.eq.1) then
              do l=1,18
                 lone = isval1(itri,l)
                 
                 sum = sum + temparr(k,l) &
                      *(phi(j1)+phi0(j01))*deni(kone)*bsqri(lone)
              enddo               ! on l
           else !on idens.eq.1
              select case(iop)
                 ! tw1
              case(15)
                 sum = 0.
                 
                 !tw2
              case(16)
                 sum = sum + y0term(iodd,i,j,k)* &
                      (phi(j1)+phi0(j01))*bsqri(kone)
                 
                 !tw3
              case(17)
                 sum = sum + x0term(iodd,i,j,k)* &
                      (phi(j1)+phi0(j01))*bsqri(kone)
                 
              end select
           endif !on idens.eq.1

        enddo                  ! on k
        go to 200

! tw4, tw5     
128     continue
        
        do k=1,18
           kone = isval1(itri,k)

           if(idens.eq.1) then
              do l=1,18
                 lone = isval1(itri,l)
                 
                 sum = sum + temparr(k,l)*bsqri(lone)* &
                      (phin(jone)*deni(kone) + phin(kone)*deni(jone))
              enddo               ! on l
           else !on idens.eq.1
              select case(iop)
                 !tw4
              case(18)
                 sum = sum + y0term(iodd,i,j,k)*phin(jone)*bsqri(kone)
                 
                 !tw5
              case(19)
                 sum = sum + x0term(iodd,i,j,k)*phin(jone)*bsqri(kone)
                 
              end select
           endif !on idens.eq.1

        enddo                  ! on k
        go to 200


140     do k=1,18
           kone = isval1(itri,k)
           k01 = isval0(itri,k)
           k1 = isvaln(itri,k)
           k2 = k1 + 6
           k3 = k2 + 6
           go to(190,190,190,190,190, 190,160,148,149,150,             &
                 160,160,160,160),iop
           ! bsqr = grad(psi)^2+I^2
148        if(numvar.lt.2) then
              sum = sum - g2term(iodd,k,j,i)                           &
                   *(phi(j1)+phi0(j01))*(phi(k1)+phi0(k01))            &
                   + d2term(iodd,i)*bzero*bzero/(18.*18.)
           else 
              sum = sum - g2term(iodd,k,j,i)                           &
                   *(phi(j1)+phi0(j01))*(phi(k1)+phi0(k01))            &
                   + k1term(iodd,i,j,k)                                &
                   *(phi(j2)+phi0(j2))*(phi(k2)+phi0(k2))
           endif
           go to 190
           ! galph = p / B^2
149        if(numvar.ge.3) then
              if(ipres.eq.1) then
                 sum = sum + k1term(iodd,i,j,k)*bsqri(kone)               &
                      *((pres(jone)+pres0(jone))-(phi(j3)+phi0(j3)))
              else
                 sum = sum + k1term(iodd,i,j,k)*bsqri(kone)               &
                      *(pi0/p0)*(phi(j3)+phi0(j3))
              endif
           else
              sum = sum + dterm(iodd,i,j)*bsqri(jone)*pi0/18.
           endif
           go to 190
           ! ggam = 3*galph/B^2
150        sum = sum + 3.*k1term(iodd,i,j,k)*galph(jone)*bsqri(kone)
           go to 190
                
160        do l=1,18
              lone = isval1(itri,l)
              l1= isvaln(itri,l)
              l2= l1 + 6
              l3= l2 + 6
              go to(180,180,180,180,180, 180,180,180,180,180,          &
                    161,162,162,162),iop
              ! gmu
161           sum = sum + temparr(k,l)                                 &
                   *(phi(j2)+phi0(j2))*(phi(k2)+phi0(k2))              &
                   *ggam(lone)
              go to 180
              ! gw1, gw2, gw3
162           sum = sum + temparr(k,l)                                 &
                   *(phi(j1)+phi0(j01))*(phi(k1)+phi0(k01))            &
                   *ggam(lone)/2.
              go to 180

180        enddo ! on l
190     enddo ! on k
200     temp(ione-s3matrix_lu%fst_row) = temp(ione-s3matrix_lu%fst_row) + sum
800     continue
     enddo ! on irect
     enddo ! on jrect
  enddo ! on j
  enddo ! on i
  enddo ! on iodd

  ! perform LU backsubstitution to get outarray solution
  if(myrank.eq.0) call second(tbeforesolve)

  call sparseR8d_solve_part(s3matrix_lu,temp,outarray,ier)

  if(myrank.eq.0) call second(taftersolve)
  tsolve = tsolve + taftersolve - tbeforesolve
  if(myrank.eq.0 .and. iprint.ge.1) write(*,*) ier, "after solve2"

  deallocate(phin,temp)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tnewvarb = tnewvarb + tend - tstart
  endif

  return
  
end subroutine newvar

!========================================================
subroutine newvarbv(inarray,outarray,mmnn,numvard,iplace,iop)
!
! define a new variable with homogeneous boundary conditions applied
!  note:  this routine is only for iop=2 because it has different bc
!         than the others
!
!      NOTE:  the operations for a given iop are defined by
!             computed go-to below
!
!.....replaces delsquared (05/04/05)...defines a new variable
!     iop=2  defines vor = delsquared(phi)
!
!
  use p_data
  use t_data
  use basic
  use arrays
  use superlu
  
  implicit none

  integer :: mmnn, numvard, iplace, iop
  real  :: inarray(mmnn*6*numvard),outarray(mmnn*6)
  integer, parameter :: r8a = selected_real_kind(12,100)
  real(r8a), allocatable:: phin(:),temp(:)
  integer, save, allocatable:: iboundds(:)
  integer, save :: nbcds
  integer :: nrads, numvards, mmnn6, ier, jer, msizeds
  integer :: irect, jrect, i, ione, j, jone, iodd, itri, ll, l, lx, lz
  integer :: m_lds, nnz_lds, ifs_rds
  real :: x, z, dbf, sum, factor, ssterm
  real :: tbeforesolve, taftersolve

  mmnn6 = 6*m*n
!  allocate(phin(mmnn6),temp(mmnn6))
  allocate(phin(mmnn6))

  ! put the scalar in inarray into a numvar=1 storage array
  do lx=1,n
     do lz=1,m
        l = lx + (lz-1)*n
        do i=1,6
           phin((l-1)*6+i) = inarray(6*((l-1)*numvard+(iplace-1))+i)
        enddo
     enddo
  enddo

  if(ifirstsa_lu.ne.0) go to 500
  ifirstsa_lu = 1
  allocate(iboundds(12*(m+n-2)))

  nrads = n*m*6
  numvards = 1
  msizeds = 9*n*m*(6*numvards)**2
  call jbdecomp1(n*m, 6, m_lds, nnz_lds, ifs_rds)
  base = 0

  ! compute LU decomposition only once

  ! form s-matrix matrices
  ss = 0
  itri = 0
  do iodd=1,ioddmx
     do ll=1,nreg
        itri = itri + 1
        do j=1,18
           do i=1,18
              jone = isval1(itri,j)
              ione = isval1(itri,i)
              ssterm = dterm(iodd,i,j)
              call increment(ss,m*n,ione,jone,ssterm,1,8)
              
           enddo
        enddo
     enddo
  enddo

  ! define indices for boundary arrays appropriate for vorticity symmetry
  call boundaryds(iboundds,nbcds,2)

  ! modify the s-matrix, inserting the boundary conditions
  do l=1,nbcds
     i = iboundds(l)
     do j= 1,nrads
        if(i.le.0 .or. i.ge.m*n*6*numvards+1) then
           write(*,3343) i,l,iboundds(l)
3343       format(" error in iboundds index",3i5)
           call safestop(3343)
        endif
        call assign(ss,n*m,i,j,0.,numvards)
     enddo
  enddo
  do l=1,nbcds
     call assign(ss,n*m,iboundds(l),iboundds(l),1.,numvards)
  enddo
  ! define the superlu matrices "samatrix_lu" and "samatrix_lu"
  ! [needs base=0....C-like]

#ifdef mpi
  call sparseR8d_init(samatrix_lu,nrads,nnz_lds,m_lds,ifs_rds,ier)
#else
  call sparseR8_init(samatrix_lu,nrads,msizeds,base,ier)
  samatrix_lu %permc_spec = 0
#endif
  call define(samatrix_lu,ss,m*n,1)


  ! perform LU decomposition of sparse matrix "samatrix"
  ! store result in opaque object "sahandle"

  jer = 0
#ifdef mpi
  call sparseR8d_new(samatrix_lu, jer)
  if(jer.ne.0) then
     write(*,*) 'after 3rd sparseR8d_new', jer
     call safestop(46)
  endif
#else
  call dsupralu_new(sahandle, samatrix_lu%values(1),               &
       samatrix_lu%irow(1), samatrix_lu%jcol_ptr(1),               &
       samatrix_lu%nnz, nrads , jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_new', jer
     call safestop(47)
  endif
  call dsupralu_colperm(sahandle, samatrix_lu% permc_spec, jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_new', jer
     call safestop(48)
  endif
  call dsupralu_lu(sahandle, jer)
  if(jer.ne.0) then
     write(*,*) 'after dsupralu_lu', jer
     call safestop(49)
  endif
#endif

 500   continue

  allocate(temp(samatrix_lu%m_loc))

  temp = 0

  ! define RHS temp
  do iodd=1,ioddmx
  do i=1,18
  do j=1,18

     do jrect=1,m-1+jper
     do irect=1,n-1+iper
        
        ll = irect + (jrect-1)*(n-1+iper)
        itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
        ione = isval1(itri,i)
        if(sparseR8d_is_local_row(samatrix_lu, ione).eq.0) goto 800

        jone = isval1(itri,j)
        
        x = (irect-iodd/3.)*deex
        z = (jrect-1+iodd/3.)*deez
        factor = 1

        if(imask.eq.1) call mask(x,z,factor)
        dbf = db*factor

        sum = 0.

        ! delsquared
111     sum = sum + aterm(iodd,i,j)*phin(jone)

        temp(ione-samatrix_lu%fst_row) = temp(ione-samatrix_lu%fst_row) + sum
800     continue
     enddo ! on irect
     enddo ! on jrect
  enddo ! on j
  enddo ! on i
  enddo ! on iodd

  ! complete the imposition of the boundary conditions,
  do l=1,nbcds
     if(sparseR8d_is_local_row(samatrix_lu, iboundds(l)).eq.1) then
!     if(iboundds(l).le.0 .or. iboundds(l).gt.mmnn6) then
!        write(*,*) "error in delsquared",l,nbcds,iboundds(l),mmnn6
!        call safestop(98)
!     endif
        temp(iboundds(l)-samatrix_lu%fst_row) = 0.
     endif
  enddo
  if(myrank.eq.0) call second(tbeforesolve)
  
  ! perform LU backsubstitution to get outarray solution
!#ifdef mpi
  call sparseR8d_solve_part(samatrix_lu,temp,outarray,ier)
!#else
!  call dsupralu_solve(sahandle,temp,ier)
!#endif
  if(myrank.eq.0) call second(taftersolve)
  tsolve = tsolve + taftersolve-tbeforesolve
  if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "after solve 5"
!  outarray = temp
  
  deallocate(phin,temp)
  
  return
  
end subroutine newvarbv


end module newvar_mod

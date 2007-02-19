#define REAL64 real

module newvar_mod

implicit none

! newvarb
! newvar


contains

subroutine newvarc(ibc,inarray,outarray,mmnn,numvard,iplace,iop)
!
!....replaces newvar(ibc=0), newvarb(ibc=1), and newvarv    11/24/06 (scj)
!
!      NOTE:  the operations for a given iop are defined by
!             computed go-to below
!
!.....defines a new variable with(ibc=1) or without(ibc=0) boundary condition
!     being applied
!     iop=1  defines J = delsquared(psi)
!     iop=2  defines vor = delsquared(phi)
!     iop=3  defines com = delsquared(chi)
!     iop=4  defines sphi1..source term in poloidal flux equation
!     iop=5  defines sphi2..source term in toroidal field equation
!     iop=6  defines sphik..source term in pressure equation involving heat flux
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
!     iop=20 defines ephi ... toroidal electric field
!     iop=21 defines eph2 ... v x B part (incompressible)
!     iop=22 defines eph3 ... v x B part (compressible)
!     iop=23 defines eph4 ... eta J part
!     iop=24 defines eph5 ... J x B part
!     iop=25 defines eph6 ... hyper-resistivity part
!     iop=26 defines initial ion velocity
!
!
!     The LU decomposition of the mass-matrix takes place the first-time 
!     called for each boundary condition type.
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

#ifndef BIT64
  integer, parameter :: r8a = selected_real_kind(12,100)
  real(r8a), allocatable :: phin(:),temp(:)
  real  inarray(mmnn*6*numvard),outarray(mmnn*6)
#else
  REAL64, allocatable :: phin(:), temp(:)
  REAL64 :: inarray(mmnn*6*numvard),outarray(mmnn*6)
#endif

  integer, save, allocatable:: iboundds(:)
  real temparr(18,18)
  integer filenum
  real ssterm, sum, terma, termb, hypp, denf, gmix, factor, x, z, dbf, termbf
  real hypf, hypi
  integer mmnn6, lx, lz, nrads, numvards, msizeds, iodd
  integer ier, jer, itri, irect, jrect, ll, ibc
  integer ir1, ir2, ir3, nbcds, jsymtype
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
  if(ibc.eq.0) then
    if(ifirsts3_lu.ne.0) go to 500
    ifirsts3_lu = 1
  else
    if(ifirsts6_lu.ne.0) go to 500
    ifirsts6_lu = 1
    allocate(iboundds(12*(m+n-2)))
  endif
  
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
!
if(ibc.eq.0) then
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
else
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
endif
500 continue
  
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

if(ibc.eq.0) then
  allocate(temp(s3matrix_lu%m_loc))
else
  allocate(temp(s6matrix_lu%m_loc))
endif

  temp = 0
  filenum = 0

  select case(iop)
  case(4)
     if(numvar.ge.2) filenum = 105   ! v7term
  case(7) 
     filenum = 105                   ! v7term
  case(11) 
     filenum = 151                   ! gmu
  case(12) 
     filenum = 107                   ! v7term
  case(13) 
     filenum = 152                   ! gw2
  case(14) 
     filenum = 153                   ! gw3
  case(15) 
     filenum = 101                   ! u3
  case(16) 
     filenum = 178                   ! w2
  case(17) 
     filenum = 179                   ! w3
  case(18) 
     filenum = 178                   ! w2
  case(19) 
     filenum = 179                   ! w3
  case default 
     filenum = 0
  end select
  
  ! position files with metric terms
  if(filenum.ne.0) rewind(filenum)

  if(iop.eq.5 .and. idens.eq.1) then
     if(iread(01).eq.1) rewind(101)  ! u3
     if(iread(03).eq.1) rewind(103)  ! u1
  endif

  if(iop.eq.6.or.iop.eq.7) then
     if(iread(09).eq.1) rewind(109)  ! u5
     if(iread(10).eq.1) rewind(110)  ! u6
  endif

  if(iop.ge.20 .and.iop.le.25) then
     if(iread(01).eq.1) rewind(101)  ! u3
     if(iread(02).eq.1) rewind(102)  ! u4
  endif
  
  do iodd=1,ioddmx
  do i=1,18
  do j=1,18

     terma = aterm(iodd,i,j)
     termb = bterm(iodd,i,j)
     hypp = hyperp*deex**2 
          
     if(filenum.ne.0) then
        read(filenum) ir1,ir2,ir3,((temparr(k,l),l=1,18),k=1,18)
     endif

     if(iop.eq.5 .and. idens.eq.1) then
        if(iread(01).eq.1) read(101) ir1,ir2,ir3,((u3terml(k,l),l=1,18),k=1,18)
        if(iread(03).eq.1) read(103) ir1,ir2,ir3,((u1terml(k,l),l=1,18),k=1,18)
     endif

     if(iop.eq.6.or.iop.eq.7) then
        if(iread(09).eq.1) read(109) ir1,ir2,ir3,((u5terml(k,l),l=1,18),k=1,18)
        if(iread(10).eq.1) read(110) ir1,ir2,ir3,((u6terml(k,l),l=1,18),k=1,18)
     endif

     if(iop.ge.20 .and. iop.le.25) then
        if(iread(01).eq.1) read(101) ir1,ir2,ir3,((u3terml(k,l),l=1,18),k=1,18)
        if(iread(02).eq.1) read(102) ir1,ir2,ir3,((u4terml(k,l),l=1,18),k=1,18)
     endif


     do jrect=1,m-1+jper
     do irect=1,n-1+iper
     
        ll = irect + (jrect-1)*(n-1+iper)
        itri = (iodd-1)*(n-1+iper)*(m-1+jper) + ll
        ione = isval1(itri,i)

        if(ibc.eq.0 .and. sparseR8d_is_local_row(s3matrix_lu, ione).eq.0) goto 800
        if(ibc.eq.1 .and. sparseR8d_is_local_row(s6matrix_lu, ione).eq.0) goto 800

        i1 = isvaln(itri,i)
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


        go to (101,101,101,104,105, 116,117,140,140,140,                &
               140,140,140,140,125, 125,125,128,128,108,                &
               109,110,111,112,113, 114),iop
        write(*,*) "Error: iop not defined for newvar: ", iop

101     sum = sum + aterm(iodd,i,j)*phin(jone)
        go to 200
!
!       sphi1 ... source term in poloidal field equation...
!                 should contain only lineval piece for lineval=1
104     continue
!
        sum = sum + etar*(aterm(iodd,i,j)-hyper*deex**2*bterm(iodd,i,j))  &
                        *(phi(j1)+(1-lineval)*phi0(j01))
!
        if(numvar.ge.2) then
           do k=1,18
              kone = isval1(itri,k)
              k1= isvaln(itri,k)
              k2= k1 + 6
     
              if(idens.eq.0) then
                 sum = sum + dbf*k0term(iodd,i,j,k)                         &
                      *( (phi(k2)+phi0(k2)) * (phi(j1)+phi0(j01))           &
                        - phi0(k2)*phi0(j01)*lineval )
              else
                 do l=1,18
                    lone = isval1(itri,l)
                    sum = sum + dbf*temparr(k,l)*deni(lone)                 &
                         *( (phi(k2)+phi0(k2)) * (phi(j1)+phi0(j01))        &
                           - phi0(k2)*phi0(j01)*lineval)
                 end do
              endif
           enddo               ! on k
        endif                  ! numvar.ge.2

        go to 200
!
!       sphi2 ... source term in toroidal field equation...
!                 should contain only lineval piece for lineval=1
105     continue
!
        sum = sum + etar*(aterm(iodd,i,j)-hyperi*deex**2*bterm(iodd,i,j))   &
                        *(phi(j2)+(1-lineval)*phi0(j2))
!
        do k=1,18
           kone = isval1(itri,k)
           k1 = isvaln(itri,k)
           k2 = k1 + 6
           k3 = k2 + 6
           k01 = isval0(itri,k)
           
           if(idens.eq.0) then
              sum = sum + dbf*(g0term(iodd,i,j,k)+g0term(iodd,i,k,j))/2.    &
                   *( (phi(j1)+phi0(j01))*(phi(k1)+phi0(k01))               &
                     - phi0(j01)*phi0(k01)*lineval)
           else
              do l=1,18
                 lone = isval1(itri,l)
                 l1 = isvaln(itri,l)
                 l2 = l1+6

                 sum = sum + dbf*u1terml(k,l)*deni(lone)                    &
                      *((phi(j1)+phi0(j01))*(phi(k1)+phi0(k01))             &
                        -phi0(j01)*phi0(k01)*lineval)
                 sum = sum + dbf*u3terml(k,l)*deni(jone)                    &
                      *((phi(k2)+phi0(k2 ))*(phi(l2)+phi0(l2))              &
                        -phi0(k2)*phi0(l2)*lineval)
              end do ! on l
              if(numvar.ge.3) then
                 sum = sum + pefac*dbf*k0term(iodd,i,j,k)*deni(jone)        &
                                      *(phi(k3)+(1-lineval)*phi0(k3))
              endif
           endif

        enddo                  ! on k
        go to 200
!
        ! iop = 6: sphik--Extra sources for pressure equation.
!                  Note: this is only used in the definition of q4.
!
        ! ------------------------------------------
116     continue        
        sum = sum + kappa*(terma-hypp*termb)*(phi(j3)+(1-linear)*phi0(j3))     
!        
        do k=1,18
           kone = isval1(itri,k)
           k1 = isvaln(itri,k)
           k2 = k1 + 6
           k3 = k2 + 6
           gmix = g2term(iodd,i,k,j)+g2term(iodd,i,j,k)
           sum = sum + (gam-1.)*kappat*(gmix - hypp*g13erm(iodd,i,k,j))   &
                               *deni(kone)*(phi(j3)+(1-linear)*phi0(j3))
        enddo                  ! on k
        go to 200

        ! iop = 7: sphip--RHS of pressure equation for inclusion in momentum 
!                         equation via r40 (also used in the definition of q4)
!                  --> needs to have heat conduction added in for this purpose
        ! ---------------------------------
117     continue
        
      hypf = hyper*deex**2
      hypi = hyperi*deex**2
        do k=1,18
           kone = isval1(itri,k)
           k1 = isvaln(itri,k)
           k2 = k1 + 6
           k3 = k2 + 6
!.......add in ohmic and viscous heating.....
!        
              if(linear.eq.1) then
              sum       = sum       + (gam-1.)*etar*                    &
     &           (k1term(iodd,i,j,k)*((jphi(kone)+jphi0(kone))          &
     &                              *(jphi(jone)+jphi0(jone))           &
     &                               -jphi0(kone)*jphi0(jone))          &
     &        -g2term(iodd,k,j,i)*((phi(k2)+phi0(k2))*(phi(j2)+phi0(j2))&
     &                                -phi0(k2)*phi0(j2))               &
     &        -hypf*g2term(iodd,k,j,i)*((jphi(kone)+jphi0(kone))        &
     &                                 *(jphi(jone)+jphi0(jone))        &
     &                                      -jphi0(kone)*jphi0(jone))   &
     &        +hypi*g12erm(iodd,i,j,k)*((phi(k2)+phi0(k2))              &
     &                                 *(phi(j2)+phi0(j2))              &
     &                                 -phi0(k2)*phi0(j2)))
              else ! on linear
              sum       = sum       + (gam-1.)*etar*(                   &
     &            k1term(iodd,i,j,k)*(jphi(kone)+jphi0(kone))           &
     &                              *(jphi(jone)+jphi0(jone))           &
     &        -g2term(iodd,k,j,i)*(phi(k2)+phi0(k2))*(phi(j2)+phi0(j2)) &
     &        -hypf*g2term(iodd,k,j,i)*(jphi(kone)+jphi0(kone))         &
     &                                 *(jphi(jone)+jphi0(jone))        &
     &        +hypi*g12erm(iodd,i,j,k)*(phi(k2)+phi0(k2))               &
     &                                *(phi(j2)+phi0(j2)))
              endif ! on linear
!
!
              sum         = sum         + (gam-1.)*amu*                 &
     &         (  g14erm(iodd,i,j,k)*vel(j1)*vel(k1)                    &
     &           -g2term(iodd,k,j,i)*vel(k2)*vel(j2))
              sum         = sum         + (gam-1.)*amu*                 &
     &            g16erm(iodd,i,j,k)*vel(j1)*vel(k3)
              sum         =sum         + (gam-1.)*amuc*                 &
     &         (  g15erm(iodd,i,j,k)*vel(j3)*vel(k3))
!
!...should be (after g17erm gets defined):
!               viscous(i3) = viscous(i3) + (gam-1.)*vel(j3)*vel(k3)*   &
!     &      (amu*g15erm(iodd,i,j,k) + 2*(amuc-amu)*g17erm(iodd,i,j,k))
           if(idens.eq.0) then
              sum = sum + dbf*k0term(iodd,i,k,j)                          &
                    *((phi(k3)+phi0(k3))*(phi(j2)+phi0(j2))               &
                      -phi0(k3)*phi0(j2)*linear)
           else
              do l=1,18
                 lone = isval1(itri,l)
                 denf = deni(lone)*dbf
                 sum = sum + denf*u5terml(k,l)*pefac                      &
                    *((phi(j3)+phi0(j3))*(phi(k2)+phi0(k2))               &
                      -phi0(j3)*phi0(k2)*linear)
              enddo               ! on l
           endif ! on idens.eq.0
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
        go to 200

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
114      continue
!....calculate the initial ion current velocity
        do k=1,18
           kone = isval1(itri,k)
           sum = sum - velion*dbf*deni(jone)*jphi0(kone)*k1term(iodd,i,j,k)
        enddo

200     continue
        if(ibc.eq.0) then
          temp(ione-s3matrix_lu%fst_row) = temp(ione-s3matrix_lu%fst_row) + sum
        else
          temp(ione-s6matrix_lu%fst_row) = temp(ione-s6matrix_lu%fst_row) + sum
        endif
800     continue
     enddo ! on irect
     enddo ! on jrect
  enddo ! on j
  enddo ! on i
  enddo ! on iodd

  ! perform LU backsubstitution to get outarray solution
  if(myrank.eq.0) call second(tbeforesolve)

        if(ibc.eq.0) then
          call sparseR8d_solve_part(s3matrix_lu,temp,outarray,ier)
        else

  ! complete the imposition of the boundary conditions,
  do l=1,nbcds
     if(sparseR8d_is_local_row(s6matrix_lu, iboundds(l)).eq.1) then
        temp(iboundds(l)-s6matrix_lu%fst_row) = 0.
     endif
  enddo
          call sparseR8d_solve_part(s6matrix_lu,temp,outarray,ier)
        endif

  if(myrank.eq.0) call second(taftersolve)
  tsolve = tsolve + taftersolve - tbeforesolve
  if(myrank.eq.0 .and. iprint.ge.1) write(*,*) ier, "after solve2"

  deallocate(phin,temp)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tnewvarb = tnewvarb + tend - tstart
  endif

  return
  
end subroutine newvarc


end module newvar_mod

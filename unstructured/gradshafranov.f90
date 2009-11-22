module gradshafranov

  implicit none

  integer, parameter :: numvargs = 1
  
  vectype, private, allocatable :: psi(:)
  vectype, private, allocatable :: fun1(:), fun2(:), fun3(:), fun4(:)
  real, private :: dpsii
  real, private :: gamma2, gamma3, gamma4  

  logical, private :: constraint = .false.

  real, private :: gnorm, libetapeff, fac2

  integer, private :: npsi ! number of radial points in input profile

  real, private, allocatable :: psinorm(:)
  real, private, allocatable :: g4big0t(:), g4bigt(:), g4bigpt(:), g4bigppt(:)
  real, private, allocatable :: fbig0t(:), fbigt(:), fbigpt(:), fbigppt(:)
  real, private, allocatable :: alphap0t(:), alphapt(:), alphappt(:), alphapppt(:)

  integer, private, parameter :: NSTX_coils = 158
  real, private, dimension(NSTX_coils) :: NSTX_r, NSTX_z, NSTX_I

  integer, private, parameter :: ITER_coils = 18
  real, private, dimension(ITER_coils) :: ITER_r, ITER_z, ITER_I

  integer, private, parameter :: DIII_coils = 18
  real, private, dimension(DIII_coils) :: DIII_r, DIII_z, DIII_I


  data NSTX_r &
       / 0.1750, 0.1750, 0.1750, 0.1750, 0.1750, & ! PF1aU
         0.1750, 0.1750, 0.1750, 0.1750, 0.1750, & ! PF1aL
         0.2902, 0.2902, 0.2902, 0.2902,         & ! PF1b
         0.3386, 0.3386, 0.3386, 0.3386,         &
         0.7255, 0.7739, 0.8223, 0.8707,         & ! PF2U
         0.7255, 0.7739, 0.8223, 0.8707,         &
         0.7255, 0.7739, 0.8223, 0.8707,         &
         0.7255, 0.7739, 0.8223, 0.8707,         & ! PF2L
         0.7255, 0.7739, 0.8223, 0.8707,         &
         0.7255, 0.7739, 0.8223, 0.8707,         &
         1.4511, 1.5478, 1.4995,                 & ! PF3U
         1.4511, 1.5478, 1.4995,                 &
         1.4511, 1.5478, 1.4995, 1.4027, 1.4027, 1.4027, &
         1.4511, 1.5478, 1.4995,                 & ! PF3L
         1.4511, 1.5478, 1.4995,                 &
         1.4511, 1.5478, 1.4995, 1.4027, 1.4027, 1.4027, &
         1.9832, 2.0315, 1.9832, 2.0315, 1.9832, 2.0315, & ! PF5U
         1.9832, 2.0315, 1.9832, 2.0315, 1.9832, 2.0315, & ! PF5L /
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &   ! OH
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240 /
         
  data NSTX_z &
       / 1.5000, 1.5500, 1.6000, 1.6500, 1.7000, & ! PF1aU
        -1.5000,-1.5500,-1.6000,-1.6500,-1.7000, & ! PF1aL
        -1.9000,-1.8500,-1.8000,-1.7500,         & ! PF1b
        -1.9000,-1.8500,-1.8000,-1.7500,         &
         1.8500, 1.8500, 1.8500, 1.8500,         & ! PF2U
         1.9000, 1.9000, 1.9000, 1.9000,         &
         1.9500, 1.9500, 1.9500, 1.9500,         &
        -1.8500,-1.8500,-1.8500,-1.8500,         & ! PF2L
        -1.9000,-1.9000,-1.9000,-1.9000,         &
        -1.9500,-1.9500,-1.9500,-1.9500,         &
         1.6500, 1.6500, 1.5500,         & ! PF3U
         1.5500, 1.5500, 1.6500,         &
         1.6000, 1.6000, 1.6000, 1.5500, 1.6000, 1.6500, &
        -1.6500,-1.6500,-1.5500,         & ! PF3L
        -1.5500,-1.5500,-1.6500,         &
        -1.6000,-1.6000,-1.6000,-1.5500,-1.6000,-1.6500, &
         0.6500, 0.6500, 0.6000, 0.6000, 0.5500, 0.5500, & ! PF5U
        -0.6500,-0.6500,-0.6000,-0.6000,-0.5500,-0.5500, & ! PF5L
         0.0000, 0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.3500, & ! OH
         0.4000, 0.4500, 0.5000, 0.5500, 0.6000, 0.6500, 0.7000, 0.7500, &
         0.8000, 0.8500, 0.9000, 0.9500, 1.0000, 1.0500, 1.1000, 1.1500, &
         1.2000, 1.2500, 1.3000, 1.3500, 1.4000, 1.4500, 1.5000, 1.5500, &
         1.6000, 1.6500, 1.7000, 1.7500, 1.8000, 1.8500, 1.9000, 1.9500, &
         2.0000,-0.0500,-0.1000,-0.1500,-0.2000,-0.2500,-0.3000,-0.3500, &
        -0.4000,-0.4500,-0.5000,-0.5500,-0.6000,-0.6500,-0.7000,-0.7500, &
        -0.8000,-0.8500,-0.9000,-0.9500,-1.0000,-1.0500,-1.1000,-1.1500, &
        -1.2000,-1.2500,-1.3000,-1.3500,-1.4000,-1.4500,-1.5000,-1.5500, &
        -1.6000,-1.6500,-1.7000,-1.7500,-1.8000,-1.8500,-1.9000,-1.9500 /

  data NSTX_I &
       / 0.0700, 0.0700, 0.0700, 0.0700, 0.0700, & ! PF1aU
         0.0700, 0.0700, 0.0700, 0.0700, 0.0700, & ! PF1aL
         0.0100, 0.0100, 0.0100, 0.0100,         & ! PF1b
         0.0100, 0.0100, 0.0100, 0.0100,         &
         0.0200, 0.0200, 0.0200, 0.0200,         & ! PF2U
         0.0200, 0.0200, 0.0200, 0.0200,         &
         0.0200, 0.0200, 0.0200, 0.0200,         &
         0.0200, 0.0200, 0.0200, 0.0200,         & ! PF2L
         0.0200, 0.0200, 0.0200, 0.0200,         &
         0.0200, 0.0200, 0.0200, 0.0200,         &
        -0.0150,-0.0150,-0.0150,                 & ! PF3U
        -0.0150,-0.0150,-0.0150,                 & 
        -0.0150,-0.0150,-0.0150,-0.0150,-0.0150,-0.0150, &         
        -0.0150,-0.0150,-0.0150,                 & ! PF3L
        -0.0150,-0.0150,-0.0150,                 &
        -0.0150,-0.0150,-0.0150,-0.0150,-0.0150,-0.0150, &
        -0.0350,-0.0350,-0.0350,-0.0350,-0.0350,-0.0350, & ! PF5U 
        -0.0350,-0.0350,-0.0350,-0.0350,-0.0350,-0.0350, & ! PF5L
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, & ! OH
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082 /

data ITER_r &
     / 3.954, 8.309, 11.935, 11.905, 8.395, 4.287, &
       1.688, 1.688,  1.688,  1.688, 1.688, 1.688, &
       1.688, 1.688,  1.688,  1.688, 1.688, 1.688 /

data ITER_z &
     /  7.557,  6.530,  3.244, -2.263, -6.730, -7.557, &
       -4.500, -5.500, -2.500, -3.500, -0.500, -1.500, &
        0.500,  1.500,  2.500,  3.500,  4.500,  5.500 /

data ITER_I &
     /  0.971,  -2e-5, -0.444, -0.177, -0.882,  1.544, &
        1.041,  1.041, -0.467, -0.467, -0.256, -0.256, &
       -0.256, -0.256,  0.283,  0.283, -0.339, -0.339 /

!!$data DIII_r &
!!$     / 0.8608, 0.8614, 0.8628, 0.8611, &
!!$       1.0041, 2.6124, 2.3733, 1.2518, 1.6890, &
!!$       0.8608, 0.8607, 0.8611, 0.8630, &
!!$       1.0025, 2.6124, 2.3834, 1.2524, 1.6889 /

! Coils 7 and 15 have been moved out 0.2 to avoid domain boundary
data DIII_r &
     / 0.8608, 0.8614, 0.8628, 0.8611, &
       1.0041, 2.6124, 2.4733, 1.2518, 1.6890, &
       0.8608, 0.8607, 0.8611, 0.8630, &
       1.0025, 2.6124, 2.4834, 1.2524, 1.6889 /


data DIII_z &
     / 0.16830, 0.50810, 0.84910, 1.1899, &
       1.5169,  0.4376,  1.1171,  1.6019, 1.5874, &
      -0.17370,-0.51350,-0.85430,-1.1957, &
      -1.5169, -0.4376, -1.1171, -1.6027,-1.5780 /

data DIII_I &
     / -0.129524,  -0.0836822, -0.00297894, 0.0974628, &
        0.0923149, -0.152552,  -0.165889,   0.0772530, -0.0400139, &
       -0.124968,  -0.0823382, -0.00382132, 0.0882870, &
        0.0898185, -0.152320,  -0.153892,   0.0711646, -0.0398154 /

contains

subroutine gradshafranov_init()

  use basic
  use arrays
  use diagnostics

  implicit none

  integer :: i, numnodes
  real :: tstart, tend

  ! Define initial values of psi
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(irestart.ne.2) then

     if(ifixedb.eq.0) call vacuum_field
     
     ! define initial field associated with delta-function source
     !     corresponding to current tcuro at location (xmag,zmag)
     call deltafun(xmag,zmag,tcuro,jphi,1,1)

  else   ! on irestart.ne.2
     psimin = -psimin
     psilim = -psilim

     call numnod(numnodes)
     do i=1,numnodes
       call assign_local_pointers(i)
       psi0_l = -psi0_l
     enddo
  endif   ! on irestart.ne.2

  fieldi = field0

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(igs.ne.0) call gradshafranov_solve()
  if(myrank.eq.0 .and. itimer.eq.1) then 
     call second(tend)
     t_gs = tend - tstart
  endif

  call gradshafranov_per()
end subroutine gradshafranov_init

subroutine gradshafranov_per()

  use basic
  use arrays
  use diagnostics

  implicit none

  integer :: i, numnodes
  real :: x, z
  vectype, dimension(6) :: vmask


  call numnod(numnodes)

  do i=1, numnodes

     call nodcoord(i, x, z)

     call assign_local_pointers(i)

     u0_l = 0.
     vz0_l = 0.
     call calc_rotation(psi0_l,vz0_l,x,z)
     chi0_l = 0.

     vmask = p0_l/p0
     vmask(1) = vmask(1) - pedge/p0
     
     ! initial parallel rotation
     u1_l = phizero*vmask

     ! allow for initial toroidal rotation
     vz1_l = 0.
     if(vzero.ne.0) call add_angular_velocity(vz1_l, x+xzero, vzero*vmask)

     ! add random perturbations
     if(nonrect.eq.0) then
        vmask(1) = 1.
        vmask(2:6) = 0.
     endif
     call random_per(x,z,23,vmask)

  enddo

end subroutine gradshafranov_per


! vacuum_field
subroutine vacuum_field()
  use basic
  use arrays
  use coils

  implicit none

  real, dimension(6) :: g1, g2
  real, dimension(maxcoils) :: xp, zp, xc, zc, ic
  real :: aminor, bv, rnorm, fac
  integer :: ineg, ipole, numcoils
  

  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " calculating vacuum field conditions..."

  ! based on filiment with current tcuro
  ! and vertical field of strength bv given by shafranov formula
  ! NOTE:  This formula assumes (li/2 + beta_P) = libetap
  fac  = tcuro/(2.*pi)
  fac2 = tcuro / (8.*pi**2*xmag)
  ! minor radius
  aminor = abs(xmag-xlim)
  if(itor.eq.1) then
     bv =  alog(8.*xmag/aminor) - 1.5 + libetap
     libetapeff = libetap
  else
     bv = 0.
  endif

  rnorm = 10.
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "gradshafranov_solve xmag zmag xzero zzero= ", &
       xmag, zmag, xzero, zzero
  
  !......define feedback parameters needed for normalization
  if(idevice .eq. 0) then
     xc(1) = 102.
     zc(1) = rnorm
     xp = xlim
     zp = zlim
     call gvect(xp,zp,xc,zc,1,g1,1,ineg)
     xp = xlim2
     zp = zlim2
     call gvect(xp,zp,xc,zc,1,g2,1,ineg)
     gnorm = g1(1) - g2(1)
  endif

  ipole = 0
  select case(idevice)
  case(-1)
     call load_coils(xc,zc,ic,numcoils,'coil.dat','current.dat')
     
  case(1) ! CDX-U
     if(myrank.eq.0) print *, "Using standard CDX-U configuration"
     numcoils = 4
     xc(1) = 0.846
     zc(1) = 0.360
     xc(2) = 0.846
     zc(2) =-0.360
     xc(3) = 0.381
     zc(3) = 0.802
     xc(4) = 0.381
     zc(4) =-0.802
     ic = -.2*fac
     
  case(2) ! NSTX
     if(myrank.eq.0) print *, "Using standard NSTX configuration"
     numcoils = NSTX_coils
     xc(1:NSTX_coils) = NSTX_r
     zc(1:NSTX_coils) = NSTX_z
     ic(1:NSTX_coils) = fac*NSTX_I

  case(3) ! ITER
     if(myrank.eq.0) print *, "Using standard ITER configuration"
     numcoils = ITER_coils
     xc(1:ITER_coils) = ITER_r
     zc(1:ITER_coils) = ITER_z
     ic(1:ITER_coils) = fac*ITER_I
     
  case(4) ! DIII
     if(myrank.eq.0) print *, "Using standard DIII-D configuration"
     numcoils = DIII_coils
     xc(1:DIII_coils) = DIII_r
     zc(1:DIII_coils) = DIII_z
     ic(1:DIII_coils) = fac*DIII_I
     
  case default ! Generic

     if(myrank.eq.0) print *, "Using generic (dipole) configuration"

     numcoils = 1
     xc(1) = 102.
     zc(1) = rnorm
     ipole = 1
     ic = bv*fac2
  end select


  ! Field due to coil currents
  call field_from_coils(xc,zc,ic,numcoils,field0,num_fields,psi_g,ipole)
 
  ! Field due to extra divertor currents
  if(divertors.ge.1) then
     xc(1:2) = xdiv
     zc(1) = zdiv
     if(divertors.eq.2) zc(2) = -zdiv
     ic(1:2) = fac*divcur
     call field_from_coils(xc,zc,ic,divertors,field0,num_fields,psi_g,0)
  endif

  ! Field due to plasma current
  xc(1) = xmag
  zc(1) = zmag
  ic(1) = tcuro/(2.*pi)
  call field_from_coils(xc,zc,ic,1,field0,num_fields,psi_g,0)
  
end subroutine vacuum_field


!============================================================
subroutine gradshafranov_solve

  use t_data
  use basic
  use arrays
  use sparse
  use diagnostics
  use newvar_mod
  use nintegrate_mod

  implicit none

#ifdef _AIX
  include 'mpif.h'
#endif
#include "finclude/petsc.h"
  
  vectype, allocatable :: b1vecini(:)
  vectype, allocatable :: b2vecini(:)

  integer :: itri,i,ii,i1,j,j1,ier, itnum
  integer :: numelms, numnodes
  integer :: ibegin, iendplusone
  real :: dterm(18,18)
  real :: feedfac, fintl(-6:maxi,-6:maxi)

  real :: x, z, error, error2, error3
  real :: sum, g0, sum2
  real :: gamma2a,gamma2b,gamma3a,gamma3b
  
  real :: tstart, tend

  integer ::  idim(3)
  real :: n(2,3)
  logical :: is_edge(3)  ! is inode on boundary
  vectype, dimension(6) :: tf
  vectype, dimension(20) :: avec

  PetscTruth :: flg_petsc, flg_solve2, flg_solve1
  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, "Calculating Grad-Shafranov Equilibrium"

  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier)

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart) ! t_gs_init

  t_gs_magaxis = 0.
  t_gs_solve = 0.
  t_gs_fundef = 0.

  g0 = bzero*rzero

  call numnod(numnodes)
  call numfac(numelms)

  ! allocate memory for arrays
  call createvec(b1vecini, numvargs)
  call createvec(b2vecini, numvargs)
  call createvec(psi, numvargs)
  call createvec(fun1, numvargs)
  call createvec(fun4, numvargs)
  call createvec(fun2, numvargs)
  call createvec(fun3, numvargs)  

  call copyvec(field0, psi_g, num_fields, psi, 1, numvargs)
  if(iread_eqdsk .eq. 1) psilim = psibound
  call copyvec(jphi, 1, 1, b1vecini, 1, numvargs)


  ! form the grad-sharfranov matrix
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " forming the GS matrix..."

  ! default linear solver superlu cj-april-09-2008
  if(flg_petsc.eq.PETSC_TRUE) then
     call zeropetscmatrix(gsmatrix_sm, icomplex, numvargs)
     if(iprint.ge.1) print *, "	gradshafranov_solve zeropetscmatrix", gsmatrix_sm
  else
     call zerosuperlumatrix(gsmatrix_sm, icomplex, numvargs)
     if(iprint.ge.1) print *, "	gradshafranov_solve zerosuperlumatrix", gsmatrix_sm
  endif

  ! populate the matrix
  do itri=1,numelms

     call define_triangle_quadrature(itri,25)
     call define_fields(itri,0,1,0)
    
     do j=1,18
        do i=1,18
           i1 = isval1(itri,i)
           j1 = isval1(itri,j)
           temp79a = -ri2_79* &
                (g79(:,OP_DR,i)*g79(:,OP_DR,j) &
                +g79(:,OP_DZ,i)*g79(:,OP_DZ,j))
           sum = int1(temp79a)
!!$           sum = int3(ri2_79,g79(:,OP_1,i),g79(:,OP_GS,j))
           call insertval(gsmatrix_sm, sum, 0, i1,j1,1)
        enddo
     enddo

     if(isurface.eq.0) cycle

     ! add surface terms
     call boundary_edge(itri, is_edge, n, idim)
     
     do ii=1,3
        if(.not.is_edge(ii)) cycle

        call define_edge_quadrature(itri, ii, 5, n, idim)
        call define_fields(itri, 0, 1,0)

        do j=1,18
           j1 = isval1(itri,j)
           do i=1,18
              i1 = isval1(itri,i)
              temp79a = norm79(:,1)*g79(:,OP_DR,j) &
                      + norm79(:,2)*g79(:,OP_DZ,j)
              sum = int3(ri2_79,g79(:,OP_1,i),temp79a)
              call insertval(gsmatrix_sm, sum, 0, i1,j1,1)
           enddo
        enddo
     end do
  enddo


  feedfac = 0.
  ! insert boundary conditions
!
!.....NOTE:   This first call just modifies the gsmatrix_sm by inserting 1's
!             on the diagonal for boundary points (or vector angles for non-rect)
!
  call boundary_gs(gsmatrix_sm, b2vecini, feedfac)
  call finalizematrix(gsmatrix_sm)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_gs_init = tend - tstart
  endif

  !....read in numerical values for p and g functions for inumgs = 1
  if(inumgs .eq. 1) then
     call readpgfiles
  else
     if(.not.allocated(psinorm)) call default_profiles
  endif


  if(myrank.eq.0) then
  write(*,999) 
999 format("    I    error        error2       xmag         psimin       psilim" &
              ,"       psilim2     xnull       znull")
  endif

  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
  mainloop: do itnum=1, iabs(igs)

     if(myrank.eq.0 .and. iprint.eq.1) print *, "GS: iteration = ", itnum
     
     ! apply boundary conditions
     if((iread_eqdsk .ne. 1 .and. irestart.ne.2) .or. itnum.gt.1) then
        feedfac = 0.
        if(itnum.gt.1 .and. gnorm.ne.0 .and. xlim2.ne.0) then
           feedfac = -0.25*(psilim - psilim2)/gnorm
           !......as a diagnostic, calculate the effective value of libetap (including feedback term)
           libetapeff =  libetapeff + feedfac/fac2
           if(myrank.eq.0 .and. iprint.eq.1) &
                write(*,'(A,4E12.4)') "feedfac, psilim, psilim2,gnorm", &
                feedfac, psilim, psilim2, gnorm
        endif
    
        call boundary_gs(0, b1vecini, feedfac)
        
        ! perform LU backsubstitution to get psi solution
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        b2vecini = b1vecini
        if(flg_petsc.eq.PETSC_TRUE .and. flg_solve1.eq.PETSC_TRUE) then
           call solve1(gsmatrix_sm,b1vecini,ier)
        else
           call solve(gsmatrix_sm,b1vecini,ier)
        endif
        if(ier.ne.0) then
           if(myrank.eq.0) print *, 'Error in GS solve'
           call safestop(10)
        end if
        if(b1vecini(1).ne.b1vecini(1)) then 
           print *, 'Error: solution is NaN'
           call safestop(11)
        endif
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_gs_solve = t_gs_solve + tend - tstart
        endif           

        ! combine solve result with old solution to get new solution
        if(itnum.eq.1) then
           psi = b1vecini
        else
           psi = th_gs*b1vecini + (1.-th_gs)*psi
        endif
     endif
     

     ! Find new magnetic axis and lcfs
     call lcfs(psi,1,1)

        
     ! define the pressure and toroidal field functions
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'calculating funs...'
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(constraint .and. igs_method.eq.2) then
        call fundef2(error3)
     else
        call fundef
     end if
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_gs_fundef = t_gs_fundef + tend - tstart
     endif

     ! Calculate error in new solution
     call calculate_error(error,error2,b1vecini)
     if(constraint .and. igs_method.eq.2) error = error3

     if(myrank.eq.0) then
        write(*,'(i5,1p8e13.5)') &
             itnum,error,error2,xmag,psimin,psilim,psilim2,xnull,znull
     endif

     ! if error is sufficiently small, stop iterating
     if(itnum .gt. 1 .and. error2 .lt. tol_gs) exit mainloop
    
     ! calculate gammas to constrain current, etc.
     call calculate_gamma(gamma2,gamma3,gamma4)
     
     ! Define RHS vector
     b1vecini = 0.
     do itri=1,numelms
        
        call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
        call calcdterm(itri, dterm, fintl)
        
        do i=1,18
           sum = 0.
           sum2 = 0.
           
           do j=1,18
              i1 = isval1(itri,i)
              j1 = isval1(itri,j)
              
              sum = sum - dterm(i,j)* &
                   (       fun1(j1) + gamma4*fun4(j1)               &
                   +gamma2*fun2(j1) + gamma3*fun3(j1))
           enddo
           
           b1vecini(i1) =  b1vecini(i1) + sum
           if(fun1(i1).ne.fun1(i1)) print *, fun1(i1)
        enddo
     enddo
     call sumsharedppplvecvals(b1vecini)

  end do mainloop
!

  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, "Converged GS error =",error2
     print *, "initial and final(effective) libetap", libetap, libetapeff
     gamma2a = -xmag**2*p0*p1
     gamma2b = -2.*(abs(g0)/(xmag*q0*dpsii))
     gamma3a = -4.*(abs(g0)/xmag)*djdpsi/dpsii
     gamma3b = -xmag**2*p0*p2

     write(*,1001) gamma2,gamma3,gamma4,xmag,p0,p1,p2,g0,dpsii,   &
          djdpsi,tcuro,gamma2a,gamma2b,gamma3a,gamma3b
1001 format(" gamma2,gamma3,gamma4       =",1p3e12.4,/,     &
          " xmag, p0,p1,p2,g0          =",1p5e12.4,/,     &
          " dpsii,djdpsi,tcurb         =",1p3e12.4,/,     &
          "gamm2a,gamm2b,gamm3a,gamm3b =",1p4e12.4)
  endif

  ! if igs is positive, stop after iabs(igs) iterations
  ! continue for igs negative
  if(itnum.eq.igs) call safestop(3)


  ! Define equilibrium fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining equilibrium fields...'
  if(igs_method.eq.2) then
     ! solve for p and f fields which best approximate gs solution
     b1vecini = 0.
     b2vecini = 0.
     do itri=1,numelms
        call define_triangle_quadrature(itri, int_pts_aux)
        call define_fields(itri, 0, 1, 0)
        
        call calcavector(itri, psi, 1, numvargs, avec)
        call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ps079)
        
        do i=1, npoints       
           call calc_pressure(ps079(i,:),tf,x_79(i),z_79(i))
           temp79a(i) = tf(1)
           call calc_toroidal_field(ps079(i,:),tf,x_79(i),z_79(i))
           temp79b(i) = tf(1)
        end do
        
        do i=1, 18
           i1 = isval1(itri,i)
           b1vecini(i1) = b1vecini(i1) + int2(g79(:,OP_1,i),temp79a)
           b2vecini(i1) = b2vecini(i1) + int2(g79(:,OP_1,i),temp79b)
        end do
     end do
     call solve_newvar(b1vecini,NV_NOBOUND,mass_matrix_lhs,b1vecini)
     call solve_newvar(b2vecini,NV_NOBOUND,mass_matrix_lhs,b2vecini)
     
     call copyvec(b1vecini, 1, 1, field0, p_g, num_fields)
     call copyvec(b2vecini, 1, 1, field0, bz_g, num_fields)
  end if
     
  do i=1,numnodes
     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call assign_local_pointers(i)

     psi0_l = psi(ibegin:ibegin+5)
     psi1_l = 0.

     if(igs_method.eq.1) then
        call nodcoord(i, x, z)
        call calc_toroidal_field(psi0_l, bz0_l, x, z)
        call calc_pressure(psi0_l, p0_l, x, z)
     end if

     pe0_l = (1. - ipres*pi0/p0)*p0_l

     call calc_density(psi0_l,p0_l,den0_l,x,z)

  end do

  ! free memory
  call deletevec(b1vecini)
  call deletevec(b2vecini)
  call deletevec(psi)
  call deletevec(fun1)
  call deletevec(fun4)
  call deletevec(fun2)
  call deletevec(fun3)

  ! calculate final error
  call calculate_gs_error(error)
  if(myrank.eq.0) print *, 'Final error in GS solution: ', error


  if(myrank.eq.0 .and. iprint.ge.1) print *, 'done gradshafranov_solve.'

end subroutine gradshafranov_solve

subroutine calculate_error(error, error2, psinew)
  use basic

  implicit none

  include 'mpif.h'

  real, intent(out) :: error, error2
  vectype, intent(in) :: psinew(*)

  integer :: i, numnodes, ibegin, iendplusone, izone, izonedim, ier
  real :: sum, sum2, norm, norm2, normal(2), curv, x, z, lhs, rhs
  logical :: is_boundary
  real, dimension(5) :: temp1, temp2

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'calculating error...'

  sum = 0.
  norm = 0.
  sum2 = 0.
  norm2 = 0.

  call numnod(numnodes)
  do i=1,numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(is_boundary) cycle

     call nodcoord(i,x,z)
        
     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     
     lhs = (psi(ibegin+3)-psi(ibegin+1)/x+psi(ibegin+5))/x
     rhs =  -(fun1(ibegin)+gamma2*fun2(ibegin)+                        &
          gamma3*fun3(ibegin)+gamma4*fun4(ibegin))

     sum = sum + (lhs-rhs)**2
     norm = norm + lhs**2
     sum2 = sum2 + abs(psi(ibegin)-psinew(ibegin))
     norm2 = norm2 + abs(psi(ibegin))
  enddo
  
  if(maxrank.gt.1) then
     temp1(1) = sum
     temp1(2) = norm
     temp1(3) = sum2
     temp1(4) = norm2
     call mpi_allreduce(temp1, temp2, 4, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ier)
     error = sqrt(temp2(1)/temp2(2))
     error2= temp2(3)/temp2(4)
  else
     error = sqrt(sum/norm)
     error2= sum2/norm2
  endif
end subroutine calculate_error

!============================================================
! calculate_gamma
! ~~~~~~~~~~~~~~~
! calculates the values of gamma2, gamma3, and gamma4 to
! constrain solution to have the specified current, etc.
!============================================================
subroutine calculate_gamma(g2, g3, g4)
  use basic
  use t_data
  use arrays

  implicit none

  include 'mpif.h'

  real, intent(out) :: g2, g3, g4

  integer :: i, j, itri, numelms, i1, jone, ier

  real :: gsint1, gsint2, gsint3, gsint4, curr, g0
  real :: sterm(18,18), fintl(-6:maxi,-6:maxi), cfac(18)
  real, dimension(5) :: temp1, temp2

  ! start of loop over triangles to compute integrals needed to keep
  !     total current and q_0 constant using gamma4, gamma2, gamma3
  if(myrank.eq.0 .and. iprint.ge.1) print *, 'calculating gammas...'

  if(nv1equ.eq.1) then
     g2 = 0.
     g3 = 0.
     g4 = 0.
     return
  endif
       
  if(constraint) then
     g2 = 0.
     g3 = 0.
     g4 = 1.
     return
  endif

  g0 = bzero*rzero

  gsint1 = 0.
  gsint4 = 0.
  gsint2 = 0.
  gsint3 = 0.
  curr = 0.

  call numfac(numelms)

  do itri=1,numelms
     call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
     call calcsterm(itri, sterm, fintl)
           
     do j=1,18
        cfac(j) = 0.
        do i=1,20
           cfac(j) = cfac(j) + gtri(i,j,itri)*fintl(mi(i),ni(i))
        enddo
     enddo
     
     do i=1,18
        i1 = isval1(itri,i)
              
        gsint1 = gsint1 + cfac(i)*fun1(i1)
        gsint4 = gsint4 + cfac(i)*fun4(i1)
        gsint2 = gsint2 + cfac(i)*fun2(i1)
        gsint3 = gsint3 + cfac(i)*fun3(i1)
              
        do j=1,18
           jone = isval1(itri,j)
           curr = curr + sterm(i,j)*psi(i1)*rinv(jone)
        enddo
     enddo
  enddo
     
  if(maxrank.gt.1) then
     temp1(1) = curr
     temp1(2) = gsint1
     temp1(3) = gsint2
     temp1(4) = gsint3
     temp1(5) = gsint4
     call mpi_allreduce(temp1, temp2, 5, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ier)
     curr   = temp2(1)
     gsint1 = temp2(2)
     gsint2 = temp2(3)
     gsint3 = temp2(4)
     gsint4 = temp2(5)
  end if
        
  if(myrank.eq.0 .and. iprint.ge.1) then 
     write(*,'(A,3E12.4)') "GS: curr, x0, z0 = ", curr, xmag, zmag
  endif
  
  ! choose gamma2 to fix q0/qstar.  Note that there is an additional
  ! degree of freedom in gamma3.  Could be used to fix qprime(0)
  
  g2 =  -xmag**2*p0*p1 - 2.*abs(g0)/(xmag*q0*abs(dpsii))
  g3 = -4.*(abs(g0)/xmag)*djdpsi/dpsii - xmag**2*p0*p2
  g4 = -(-tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4
end subroutine calculate_gamma


! ===========================================================
subroutine deltafun(x,z,val,dum,iplace,isize)

  use t_data
  use basic
  use arrays

  implicit none

  integer, intent(in) :: iplace, isize
  real, intent(in) :: x, z, val
  vectype, intent(out) :: dum(*)

  integer :: itri, i, ii, iii, k, index, ibegin, iendplusone
  real :: x1, z1, b, theta, si, eta, sum
  
  call whattri(x,z,itri,x1,z1)

  if(itri.gt.0) then

     ! calculate local coordinates
     theta = ttri(itri)
     b = btri(itri)
     si  = (x-x1)*cos(theta) + (z-z1)*sin(theta) - b
     eta =-(x-x1)*sin(theta) + (z-z1)*cos(theta)

     ! calculate the contribution to b1vecini
     do iii=1,3
        call entdofs(isize, ist(itri,iii)+1, 0, ibegin, iendplusone)
        do ii=1,6
           i = (iii-1)*6 + ii
           index = ibegin + ii - 1 + 6*(iplace-1)

           sum = 0.
           do k=1,20
              sum = sum + gtri(k,i,itri)*si**mi(k)*eta**ni(k)
           enddo
           dum(index) = dum(index) - sum*val
        enddo
     end do
  end if

  call sumsharedppplvecvals(dum)

end subroutine deltafun
!============================================================
subroutine fundef

!.....defines the source functions for the GS equation:
  ! fun1 = r*p'
  ! fun4 = G4'/r
  ! fun2 = G2'/r
  ! fun3 = G3'/r

  use basic
  use diagnostics
  
  implicit none 
  integer :: l, numnodes, i, ibegin, iendplusone
  real :: x, z, pso, psox, psoy, psoxx, psoxy, psoyy, fbig, fbig0
  real :: fbigp, fbigpp, g4big0, g4big, g4bigp, g4bigpp, g2big, g2bigp
  real :: g2bigpp, g3big, g3bigp, g3bigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha,  pspx, pspy, pspxx, pspxy, pspyy
  logical :: inside_lcfs

  vectype, dimension(6) :: temp

  dpsii = 1./(psilim - psimin)

  call numnod(numnodes)
  do l=1,numnodes

     call nodcoord(l, x, z)

     call entdofs(numvargs, l, 0, ibegin, iendplusone)

     temp = psi(ibegin:ibegin+5)
     pso =  (psi(ibegin)-psimin)*dpsii

     if(.not.inside_lcfs(temp,x,z,.true.)) then
        do i=0,5
           fun1(ibegin+i) = 0.
           fun4(ibegin+i) = 0.
           fun2(ibegin+i) = 0.
           fun3(ibegin+i) = 0.
        enddo
     else
        psox = psi(ibegin+1)*dpsii
        psoy = psi(ibegin+2)*dpsii
        psoxx= psi(ibegin+3)*dpsii
        psoxy= psi(ibegin+4)*dpsii
        psoyy= psi(ibegin+5)*dpsii
     
        call fget(pso, fbig0, fbig, fbigp, fbigpp)

        if(inumgs.eq.0) then
           fbig = fbig*dpsii
           fbigp = fbigp*dpsii
           fbigpp = fbigpp*dpsii
        endif
      if(irot.eq.1) then
!.....include toroidal rotation in equilibrium
!...this section calculates pressure derivatives in presence of rotation
!   scj   11/19/10
!
!....assumes pressure of the form p(x,psi) = p(psi) exp (alpha(psi)*r1)
!
!...spatial derivatives of psi, not normalized psi
     pspx = psi(ibegin+1)
     pspy = psi(ibegin+2)
     pspxx= psi(ibegin+3)
     pspxy= psi(ibegin+4)
     pspyy= psi(ibegin+5)

     r0m = 1./rzero**2
     r1 = (x**2-rzero**2)/rzero**2
     r1m= x**2/rzero**2
     r2 = (x**2 - rzero**2)**2/rzero**4
     r3 = (x**2 - rzero**2)**3/rzero**6
     call alphaget(pso,alphap0,alphap,alphapp,alphappp)

!...convert all derivatives to wrt psi, not normalized psi
     fbigp = fbigp*dpsii
     fbigpp= fbigpp*dpsii**2
     alphap = alphap*dpsii
     alphapp = alphapp*dpsii**2
     alphappp = alphapp*dpsii**3

     ealpha = exp(alphap0*r1)

         fun1(ibegin) = ealpha*( x *fbig + fbig0*alphap*x*r1)

         fun1(ibegin+1) = ealpha*                                           &
                (fbig + fbig0*alphap*x*r1                                   &
                +fbig0*alphap0*alphap*2*r1m*r1                              &
           + (alphap0*fbig + fbig0*alphap) * 2.*r1m                         &
           + pspx*( x*fbigp + (2.*fbig*alphap + fbig0*alphapp)*x*r1         &
            + fbig0 * alphap**2 * x*r2))

         fun1(ibegin+2) = ealpha*                                           &
            pspx*( x*fbigp + (2.*fbig*alphap + fbig0*alphapp)*x*r1          &
            + fbig0 * alphap**2 * x*r2)

         fun1(ibegin+3) = ealpha*                                           &
                 ((3*alphap0*fbig + 3*fbig0*alphap)*2*x*r0m                 &
                 + 3*fbig0*alphap0*alphap*2*x*r0m*r1                        &
                 + fbig0*alphap0**2*alphap*4*r1m*r0m*r1                     &
                 + alphap0*(alphap0*fbig+2*fbig0*alphap)*4.*r1m*r0m         &
                 + pspx*( 2*fbigp                                           &
                +(2*fbig0*alphapp + 4*alphap*fbig)*r1                       &
                +(4*alphap*fbig + 2*alphap0*fbigp + 2*fbig0*alphapp)*2*r1m  &
                +(4*alphap*fbig + 2*alphap0*fbigp + 2*fbig0*alphapp)*r1m    &
                +2*fbig0*alphap**2*r2                                       &
                +(4*fbig0*alphap**2 + 2*fbig0*alphap0*alphapp               &
                                            + 4*alphap0*fbig*alphap)*r1m*r1 &
                + 2*fbig0*alphap0*alphap**2*2*r1m*r2)                       &
                + pspx*pspx*(x*fbigpp                                       &
                +(3*fbigp*alphap + 3*fbig*alphapp + fbig0*alphappp)*x*r1    &
                +3*(fbig0*alphapp + fbig*alphap)*alphap*x*r2                &
                + fbig0*alphap**3*x*r3)                                     &  
                   + pspxx*(x*fbigp                                         &
                + (2*fbig*alphap + fbig0*alphapp)*x*r1                      &
                + p0*alphap**2*x*r2))

         fun1(ibegin+4) = ealpha*                                           &
                (pspy*(fbigp                                                &
              + (2.*fbig*alphap + fbig0*alphapp)*r1                         &
              + (2.*fbig*alphap + fbig0*alphapp + fbigp*alphap0)*2*r1m      &
              + (2*fbig*alphap*alphap0 + fbig0*alphapp*alphap0              &
                                            + 2*fbig0*alphap**2)*2*r1m*r1   &
              +  fbig0*alphap**2*r2                                         &
              +  fbig0*alphap**2*alphap0*2*r1m*r2)                          &
              +pspx*pspy*(x*fbigpp                                          &
              + (3*fbigp*alphap + 3*fbig*alphapp + fbig0*alphappp)*x*r1     &
              + (3*fbig*alphap + 3*fbig0*alphapp)*alphap*x*r2               &
              + fbig0*alphap**3*x*r3)                                       &
              + pspxy*(x*fbigp                                              &
                + (2*fbig*alphap + fbig0*alphapp)*x*r1                      &
                + p0*alphap**2*x*r2))

         fun1(ibegin+5) = ealpha*                                           &
                (pspy*pspy*(x*fbigpp                                        &
              + (3*fbigp*alphap + 3*fbig*alphapp + fbig0*alphappp)*x*r1     &
              + (3*fbig*alphap + 3*fbig0*alphapp)*alphap*x*r2               &
              + fbig0*alphap**3*x*r3)                                       &
              + pspyy*(x*fbigp                                              &
                + (2*fbig*alphap + fbig0*alphapp)*x*r1                      &
                + p0*alphap**2*x*r2))

      else
!
!....no toroidal rotation in equilibrium
        fun1(ibegin)   = x*fbig
        fun1(ibegin+1) = fbig + x*fbigp*psox
        fun1(ibegin+2) =        x*fbigp*psoy
        fun1(ibegin+3) = 2.*fbigp*psox + x*(fbigpp*psox**2+fbigp*psoxx)
        fun1(ibegin+4) = fbigp*psoy + x*(fbigpp*psox*psoy +fbigp*psoxy)
        fun1(ibegin+5) = x*(fbigpp*psoy**2 + fbigp*psoyy)

      endif   !...end of branch on irot

        call g4get(pso, g4big0, g4big, g4bigp, g4bigpp)

        if(inumgs.eq.0) then
           g4big = g4big*dpsii
           g4bigp = g4bigp*dpsii
           g4bigpp = g4bigpp*dpsii
        endif

        fun4(ibegin)  = g4big/x
        fun4(ibegin+1)= g4bigp*psox/x - g4big/x**2
        fun4(ibegin+2)= g4bigp*psoy/x
        fun4(ibegin+3)= (g4bigpp*psox**2 + g4bigp*psoxx)/x                  &
             - 2.*g4bigp*psox/x**2 + 2.*g4big/x**3
        fun4(ibegin+4)= (g4bigpp*psox*psoy+g4bigp*psoxy)/x                  &
             - g4bigp*psoy/x**2
        fun4(ibegin+5)=  (g4bigpp*psoy**2 + g4bigp*psoyy)/x

        g2big =  dpsii*(1. - 30.*pso**2 + 80.*pso**3                     &
             - 75.*pso**4 + 24.*pso**5)
        g2bigp =  dpsii*(-60.*pso + 240.*pso**2                         &
             - 300.*pso**3 + 120.*pso**4)
        g2bigpp =  dpsii*(-60. + 480.*pso                               &
             - 900.*pso**2 + 480.*pso**3)
        fun2(ibegin)  =  g2big/x
        fun2(ibegin+1)=  g2bigp*psox/x - g2big/x**2
        fun2(ibegin+2)=  g2bigp*psoy/x
        fun2(ibegin+3)=  (g2bigpp*psox**2 + g2bigp*psoxx)/x                 &
             - 2.*g2bigp*psox/x**2 + 2.*g2big/x**3
        fun2(ibegin+4)=(g2bigpp*psox*psoy+g2bigp*psoxy)/x                   &
             - g2bigp*psoy/x**2
        fun2(ibegin+5)= (g2bigpp*psoy**2 + g2bigp*psoyy)/x
        
        g3big =  dpsii*(2.*pso - 12.*pso**2 + 24.*pso**3                &
             - 20.*pso**4 + 6.*pso**5)
        g3bigp =  dpsii*(2. - 24.*pso + 72.*pso**2                      &
             - 80.*pso**3 + 30.*pso**4)
        g3bigpp =  dpsii*(- 24. + 144.*pso                              &
             - 240.*pso**2 + 120.*pso**3)
        fun3(ibegin)= g3big/x
        fun3(ibegin+1)= g3bigp*psox/x - g3big/x**2
        fun3(ibegin+2)= g3bigp*psoy/x
        fun3(ibegin+3)= (g3bigpp*psox**2 + g3bigp*psoxx)/x                  &
             - 2.*g3bigp*psox/x**2 + 2.*g3big/x**3
        fun3(ibegin+4)= (g3bigpp*psox*psoy+g3bigp*psoxy)/x                  &
             - g3bigp*psoy/x**2
        fun3(ibegin+5)=  (g3bigpp*psoy**2 + g3bigp*psoyy)/x
     endif
  enddo
  
  return
end subroutine fundef


subroutine fundef2(error)

  use basic
  use t_data
  use arrays
  use sparse
  use newvar_mod
  use nintegrate_mod

  implicit none

  include 'mpif.h'

  real, intent(out) :: error

  integer :: i, ione, itri, numelms, ier
  real :: pso, dpsii, norm, temp1(2), temp2(2)
  vectype, dimension(20) :: avec
      
  logical :: inside_lcfs

  dpsii = 1./(psilim - psimin)
  fun1 = 0.
  fun2 = 0.
  fun3 = 0.
  fun4 = 0.
  norm = 0.
  error = 0.

  call numfac(numelms)

  do itri=1,numelms

     call define_triangle_quadrature(itri, int_pts_aux)
     call define_fields(itri, 0, 1, 0)

     call calcavector(itri, psi, 1, numvargs, avec)
     call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79, npoints, ps079)

     do i=1, npoints
        
        pso = (ps079(i,OP_1)-psimin)*dpsii
        if(inside_lcfs(ps079(i,:),x_79(i),z_79(i),.true.)) then
           call cubic_interpolation(npsi,psinorm,pso,fbigt,temp79a(i))
           call cubic_interpolation(npsi,psinorm,pso,g4bigt,temp79b(i))
        else
           temp79a(i) = 0.
           temp79b(i) = 0.
        endif
     end do
     
     if(inumgs.eq.0) then
        temp79a = temp79a*dpsii
        temp79b = temp79b*dpsii
     endif

     do i=1,18
        ione = isval1(itri,i)
        
        fun1(ione) = fun1(ione) + int3(r_79, g79(:,OP_1,i),temp79a)
        fun4(ione) = fun4(ione) + int3(ri_79,g79(:,OP_1,i),temp79b)
     end do
     
     temp79c = ps079(:,OP_GS) + r2_79*temp79a + temp79b

     norm = norm + abs(int2(ri_79,ps079(:,OP_GS)))
     error = error + abs(int2(ri_79,temp79c))
  end do

  if(myrank.eq.0 .and. iprint.eq.1) print *, 'solving funs...'
  call solve_newvar(fun1, NV_NOBOUND, mass_matrix_lhs, fun1)
  call solve_newvar(fun4, NV_NOBOUND, mass_matrix_lhs, fun4)

  if(maxrank.gt.1) then
     temp1(1) = norm
     temp1(2) = error
     call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ier)
     norm     = temp2(1)
     error    = temp2(2)
  end if
  if(myrank.eq.0 .and. iprint.ge.1) &
       write(*,'(A,1p2E12.4)') 'Error, norm: ', error, norm
  error = error / norm

end subroutine fundef2


!======================================================================
! calc_toroidal_field
! ~~~~~~~~~~~~~~~~~~~
!
! calculates the toroidal field (I) as a function of the 
! normalized flux
!======================================================================
subroutine calc_toroidal_field(psi0,tf,x,z)
  use basic

  vectype, intent(in), dimension(6)  :: psi0
  vectype, intent(out), dimension(6) :: tf    ! toroidal field (I)
  real, intent(in) :: x, z
  
  vectype :: g0
  real, dimension(6) :: g2, g3, g4
  real :: g4big0, g4big, g4bigp, g4bigpp
  real :: g2big, g2bigp, g3big, g3bigp
  real, dimension(6)  :: psii     ! normalized flux

  logical :: inside_lcfs
  
  if(.not.inside_lcfs(psi0,x,z,.true.)) then
     g0 = bzero*rzero
     call constant_field(tf,g0)
  else
     psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
     psii(2:6) = real(psi0(2:6))/(psilim - psimin)

     if(.not.constraint) then
        g2(1) = psii(1) - 10.*psii(1)**3 + 20.*psii(1)**4       &
             - 15.*psii(1)**5 + 4.*psii(1)**6
        g2big =  (1. - 30.*psii(1)**2 + 80.*psii(1)**3          &
             - 75.*psii(1)**4 + 24.*psii(1)**5)
        g2bigp =  (-60.*psii(1) + 240.*psii(1)**2               &
             - 300.*psii(1)**3 + 120.*psii(1)**4)
        g2(2) = (psii(2))*g2big
        g2(3) = (psii(3))*g2big
        g2(4) = (psii(4)*g2big + psii(2)**2*g2bigp)
        g2(5) = (psii(5)*g2big + psii(2)*psii(3)*g2bigp)
        g2(6) = (psii(6)*g2big + psii(3)**2*g2bigp)
     
        g3(1) = psii(1)**2 - 4.*psii(1)**3 + 6.*psii(1)**4      &
             - 4.*psii(1)**5 + psii(1)**6
        
        g3big =  (2.*psii(1) - 12.*psii(1)**2 + 24.*psii(1)**3  &
             - 20.*psii(1)**4 + 6.*psii(1)**5)
        g3bigp =  (2. - 24.*psii(1) + 72.*psii(1)**2            &
             - 80.*psii(1)**3 + 30.*psii(1)**4)
        g3(2) = (psii(2))*g3big
        g3(3) = (psii(3))*g3big
        g3(4) = (psii(4)*g3big + psii(2)**2*g3bigp)
        g3(5) = (psii(5)*g3big + psii(2)*psii(3)*g3bigp)
        g3(6) = (psii(6)*g3big + psii(3)**2*g3bigp)
     end if
     
     call g4get(psii(1), g4big0, g4big, g4bigp, g4bigpp)
     
!.....changed 07/05/09 scj  for inumgs.ne.0 g4big is derivative wrt total psi, not normalized
     if(inumgs.ne.0) then
        g4big = g4big/dpsii
        g4bigp = g4bigp/dpsii
        g4bigpp = g4bigpp/dpsii
     endif
     
     g4(1) = g4big0
     g4(2) = (psii(2))*g4big
     g4(3) = (psii(3))*g4big
     g4(4) = (psii(4)*g4big + psii(2)**2*g4bigp)
     g4(5) = (psii(5)*g4big + psii(2)*psii(3)*g4bigp)
     g4(6) = (psii(6)*g4big + psii(3)**2*g4bigp)
     
!
!.....convert from gg' = .5(g^2)' to (g^2)'
     g2 = 2.*g2
     g3 = 2.*g3
     g4 = 2.*g4
     
!
     tf(1) = sqrt((bzero*rzero)**2 + &
          gamma2*g2(1) + gamma3*g3(1) + gamma4*g4(1))
     tf(2) = 0.5*(gamma2*g2(2) + gamma3*g3(2) + gamma4*g4(2)) / tf(1)
     tf(3) = 0.5*(gamma2*g2(3) + gamma3*g3(3) + gamma4*g4(3)) / tf(1)
     tf(4) = 0.5*(gamma2*g2(4) + gamma3*g3(4) + gamma4*g4(4)) / tf(1) &
          - (0.5*(gamma2*g2(2) + gamma3*g3(2) + gamma4*g4(2)))**2 / tf(1)**3
     tf(5) = 0.5*(gamma2*g2(5) + gamma3*g3(5) + gamma4*g4(5)) / tf(1) &
          -  0.5*(gamma2*g2(2) + gamma3*g3(2) + gamma4*g4(2)) &
          *  0.5*(gamma2*g2(3) + gamma3*g3(3) + gamma4*g4(3)) / tf(1)**3
     tf(6) = 0.5*(gamma2*g2(6) + gamma3*g3(6) + gamma4*g4(6)) / tf(1) &
          - (0.5*(gamma2*g2(3) + gamma3*g3(3) + gamma4*g4(3)))**2 / tf(1)**3
  
     if(bzero.lt.0) tf = -tf
  endif

end subroutine calc_toroidal_field


!======================================================================
! calc_pressure
! ~~~~~~~~~~~~~
!
! calculates the pressure as a function of the poloidal flux (and major radius if rotation is present)
!======================================================================
subroutine calc_pressure(psi0,pres, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(6)  :: psi0
  real, intent(in) :: x, z
  vectype, intent(out), dimension(6) :: pres     ! pressure

  real :: fbig0, fbig, fbigp, fbigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha,  pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(6) :: psii     ! normalized flux

  logical :: inside_lcfs

  if(.not.inside_lcfs(psi0,x,z,.true.)) then
     call constant_field(pres,pedge)
  else
     psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
     psii(2:6) = real(psi0(2:6))/(psilim - psimin)

     call fget(psii(1), fbig0, fbig, fbigp, fbigpp)
!
!.....changed 07/05/09 scj  for inumgs.ne.0 fbig is derivative wrt total psi, not normalized
     if(inumgs.ne.0) then
        fbig = fbig/dpsii
        fbigp = fbigp/dpsii
     endif
!
!
  if(irot.eq.1) then
!.....include toroidal rotation in equilibrium
!...this section calculates pressure derivatives in presence of rotation
!   scj   11/19/10
!
!....assumes pressure of the form p(x,psi) = p(psi) exp (alpha(psi)*r1)
!
!...spatial derivatives of psi, not normalized psi
     pspx = real(psi0(2))
     pspy = real(psi0(3))
     pspxx= real(psi0(4))
     pspxy= real(psi0(5))
     pspyy= real(psi0(6))

     r0m = 1./rzero**2
     r1 = (x**2-rzero**2)/rzero**2
     r1m= x**2/rzero**2
     r2 = (x**2 - rzero**2)**2/rzero**4
     r3 = (x**2 - rzero**2)**3/rzero**6
     call alphaget(psii(1),alphap0,alphap,alphapp,alphappp)

!...convert all derivatives to wrt psi, not normalized psi
     fbig = fbig*dpsii
     fbigp = fbigp*dpsii**2
     alphap = alphap*dpsii
     alphapp = alphapp*dpsii**2

     ealpha = exp(alphap0*r1)

     pres(1) = ealpha*fbig0

     pres(2) = ealpha*(fbig0*alphap0*2*x*r0m                                                &
                     + (fbig + fbig0*alphap*r1)*pspx)

     pres(3) = ealpha*(fbig + fbig0*alphap*r1)*pspy

     pres(4) = ealpha*(                                                                     &
             fbig0*alphap0*2*r0m + fbig0*alphap0*alphap0*4*r0m*r1m                          &
            +((2*fbig0*alphap + 2.*fbig*alphap0) + 2*fbig0*alphap0*alphap*r1)*2*x*r0m*pspx  &
            +(fbigp + (2*fbig*alphap + fbig0*alphapp)*r1 + fbig0*alphap**2*r2)*pspx*pspx    &
            +(fbig + fbig0*alphap*r1)*pspxx)

     pres(5) = ealpha*(                                                                     &
             (fbig*alphap0 + fbig0*alphap +fbig0*alphap0*alphap*r1)*2*x*r0m*pspy       &
            +(fbigp + (2.*fbig*alphap + fbig0*alphapp)*r1 + fbig0 *alphap**2*r2)*pspx*pspy  &
            +(fbig + fbig0*alphap*r1)*pspxy)

     pres(6) = ealpha*(                                                                     &
            +(fbigp + (2.*fbig*alphap + fbig0*alphapp)*r1 + fbig0 *alphap**2*r2)*pspy*pspy  &
            +(fbig + fbig0*alphap*r1)*pspyy)

  else

     pres(1) = fbig0
     pres(2) = psii(2)*fbig
     pres(3) = psii(3)*fbig
     pres(4) = (psii(4)*fbig + psii(2)**2*fbigp)
     pres(5) = (psii(5)*fbig + psii(2)*psii(3)*fbigp)
     pres(6) = (psii(6)*fbig + psii(3)**2*fbigp)
  endif     !....end of branch on irot
  
  endif
  pres(1) = pres(1) + pedge
  
  return
end subroutine calc_pressure

subroutine readpgfiles
  use basic

  implicit none

  integer :: j

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Reading profiles files"

  open(unit=76,file="profiles-p",status="old")
  read(76,803) npsi
  allocate(psinorm(npsi))
  allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
  do j=1,npsi
    read(76,802) psinorm(j),fbig0t(j),fbigt(j),fbigpt(j),fbigppt(j)
  enddo
  close(76)
!
!.....p0 is central pressure.   This is used in calc_density and calc_rotation if irot=1
      p0 = fbig0t(1)

!
!.....note:  removed this on 07/03/09...fbig0t is now total pressure
! fbig0t = p0*fbig0t

  open(unit=77,file="profiles-g",status="old")
  read(77,804) npsi
  allocate(g4big0t(npsi),g4bigt(npsi),g4bigpt(npsi),g4bigppt(npsi))
  do j=1,npsi
    read(77,802) psinorm(j),g4big0t(j),g4bigt(j),g4bigpt(j),g4bigppt(j)
  enddo
  close(77)

  constraint = .true.

return
  802 format(5x,5e18.10)
  803 format(i5)
  804 format(i5)
end subroutine readpgfiles

subroutine g4get(pso, g4big0, g4big, g4bigp, g4bigpp)
  implicit none

  real, intent(in) :: pso
  real, intent(out) :: g4big0,g4big, g4bigp, g4bigpp

  call cubic_interpolation(npsi,psinorm,pso,g4big0t,g4big0)
  call cubic_interpolation(npsi,psinorm,pso,g4bigt,g4big)
  call cubic_interpolation(npsi,psinorm,pso,g4bigpt,g4bigp)
  call cubic_interpolation(npsi,psinorm,pso,g4bigppt,g4bigpp)
  return
end subroutine g4get

subroutine fget(pso, fbig0, fbig, fbigp, fbigpp)
  implicit none

  real, intent(in) :: pso
  real, intent(out) :: fbig0,fbig, fbigp, fbigpp

  call cubic_interpolation(npsi,psinorm,pso,fbig0t,fbig0)
  call cubic_interpolation(npsi,psinorm,pso,fbigt,fbig)
  call cubic_interpolation(npsi,psinorm,pso,fbigpt,fbigp)
  call cubic_interpolation(npsi,psinorm,pso,fbigppt,fbigpp)  
  return
end subroutine fget

subroutine alphaget(pso, alphap0, alphap, alphapp, alphappp)
  implicit none

  real, intent(in) :: pso
  real, intent(out) :: alphap0,alphap, alphapp, alphappp

  call cubic_interpolation(npsi,psinorm,pso,alphap0t,alphap0)
  call cubic_interpolation(npsi,psinorm,pso,alphapt,alphap)
  call cubic_interpolation(npsi,psinorm,pso,alphappt,alphapp)
  call cubic_interpolation(npsi,psinorm,pso,alphapppt,alphappp)
  return
end subroutine alphaget

 subroutine default_profiles
   
   use basic

   implicit none

   integer :: j
   real :: psii

   if(myrank.eq.0) print *, "Using default p, alpha, and I profiles"

   npsi = 500
   allocate(psinorm(npsi))
   allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
   allocate(alphap0t(npsi),alphapt(npsi),alphappt(npsi),alphapppt(npsi))
   allocate(g4big0t(npsi),g4bigt(npsi),g4bigpt(npsi),g4bigppt(npsi))
  
   do j=1, npsi
      psii = (j-1.)/(npsi-1.)
      psinorm(j) = psii
      fbig0t(j) = p0*(1.+p1*psii+p2*psii**2 &
           - (20. + 10.*p1 + 4.*p2)*psii**3 &
           + (45. + 20.*p1 + 6.*p2)*psii**4 &
           - (36. + 15.*p1 + 4.*p2)*psii**5 &
           + (10. +  4.*p1 +    p2)*psii**6)
      fbigt(j) = p0*(p1 + 2.*p2*psii - 3.*(20. + 10.*p1+4.*p2)*psii**2      &
           + 4.*(45.+20.*p1+6.*p2)*psii**3 - 5.*(36.+15.*p1+4.*p2)*psii**4  &
           + 6.*(10.+4.*p1+p2)*psii**5)
      fbigpt(j) = p0*(2.*p2 - 6.*(20. + 10.*p1+4.*p2)*psii                  &
           + 12.*(45.+20.*p1+6*p2)*psii**2 - 20.*(36.+15*p1+4*p2)*psii**3   &
           + 30.*(10.+4.*p1+p2)*psii**4)
      fbigppt(j) = p0*(- 6.*(20. + 10.*p1+4.*p2)                            &
           + 24.*(45.+20.*p1+6*p2)*psii - 60.*(36.+15*p1+4*p2)*psii**2      &
           + 120.*(10.+4.*p1+p2)*psii**3)

      alphap0t(j) = alpha0 + alpha1*psii + alpha2*psii**2
      alphapt(j)  =          alpha1   + 2.*alpha2*psii
      alphappt(j) =   2.*alpha2
      alphapppt(j) = 0.

      g4big0t(j) = 1.- 20.*psii**3+  45.*psii**4-  36.*psii**5+  10.*psii**6
      g4bigt(j)  =   - 60.*psii**2+ 180.*psii**3- 180.*psii**4+  60.*psii**5
      g4bigpt(j) =   -120.*psii   + 540.*psii**2- 720.*psii**3+ 300.*psii**4
      g4bigppt(j)=   -120.        +1080.*psii   -2160.*psii**2+1200.*psii**3
   end do

 end subroutine default_profiles

!================================================================
! create_profile
! ~~~~~~~~~~~~~~
!
! Sets up the GS solver to use a specific profile for p and I.
! n = number of radial points in profile
! p = pressure profile
! pp = p' = dp/dpsi (with psi the non-normalized flux)
! f = I = R*B_phi
! ffp =  I*I'
! flux = psi
!================================================================
 subroutine create_profile(n, p, pp, f, ffp, flux,myrankt)
   implicit none

   integer, intent(in) :: n , myrankt
   real, dimension(n), intent(in) :: p, pp, f, ffp, flux

   real, dimension(4) :: a
   real :: dp,dpp
   integer :: j

   npsi = n

   allocate(psinorm(npsi))
   allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
   allocate(alphap0t(npsi),alphapt(npsi),alphappt(npsi),alphapppt(npsi))
   allocate(g4big0t(npsi),g4bigt(npsi),g4bigpt(npsi),g4bigppt(npsi))

   fbig0t(1:n) = p                                 ! p
   fbigt(1:n) = pp * (flux(n) - flux(1))        ! p' = dp/dPsi
   g4big0t(1:n) = 0.5*(f**2 - f(n)**2)          ! g
   g4bigt(1:n) = ffp * (flux(n) - flux(1))      ! f f' = f * df/dPsi

   do j=1,n
      call cubic_interpolation_coeffs(flux,n,j,a)
      psinorm(j) = (flux(j) - flux(1)) / (flux(n) - flux(1))
      dp = a(2)      ! d psi/dn
      dpp = 2.*a(3)  ! d^2 psi/dn^2

      call cubic_interpolation_coeffs(fbigt,n,j,a)
      fbigpt(j) =     a(2)/dp                       ! p''
      fbigppt(j) = 2.*a(3)/dp**2 - a(2)*dpp/dp**3   ! p'''

      call cubic_interpolation_coeffs(g4bigt,n,j,a)
      g4bigpt(j)  =    a(2)/dp                      ! (f f')'
      g4bigppt(j) = 2.*a(3)/dp**2 - a(2)*dpp/dp**3  ! (f f')''
   end do

   fbigpt = fbigpt*(flux(n) - flux(1))
   fbigppt = fbigppt*(flux(n) - flux(1))**2
   g4bigpt = g4bigpt*(flux(n) - flux(1))
   g4bigppt = g4bigppt*(flux(n) - flux(1))**2

 if( myrankt.eq.0) then
    !
    open(unit=76,file="profilesdb-p",status="unknown")
    do j=1,npsi
       write(76,802) j,psinorm(j),fbig0t(j),fbigt(j),fbigpt(j),fbigppt(j)
    enddo
    close(76)
    !
    open(unit=77,file="profilesdb-g",status="unknown")
    do j=1,npsi
       write(77,802) j,psinorm(j),g4big0t(j),g4bigt(j),g4bigpt(j),g4bigppt(j)
    enddo
802 format(i5,1p6e18.10)
    close(77)
 endif
  !
   constraint = .true.

 end subroutine create_profile

 subroutine calculate_gs_error(error)
   use basic
   use nintegrate_mod
   
   implicit none
   
   include 'mpif.h'
   
   real, intent(out) :: error
   
   real :: norm, temp1(2), temp2(2)
   integer :: itri, nelms, def_fields, ier, ieqs_temp
   
   norm = 0.
   error = 0.
   
   ! must set eqsubtract=1 so that equilibrium field is read
   ieqs_temp = eqsubtract
   eqsubtract = 1

   def_fields = FIELD_PSI + FIELD_P + FIELD_I
   
   call numfac(nelms)
   do itri=1, nelms
      call define_triangle_quadrature(itri, int_pts_main)
      call define_fields(itri, def_fields, 0, 1)
      
      temp79a = ps079(:,OP_GS)*(ps079(:,OP_DR)**2 + ps079(:,OP_DZ)**2)
      temp79b = -r2_79* &
           (p079(:,OP_DR)*ps079(:,OP_DR) + p079(:,OP_DZ)*ps079(:,OP_DZ))
      temp79c = -bz079(:,OP_1)* &
           (bz079(:,OP_DR)*ps079(:,OP_DR) + bz079(:,OP_DZ)*ps079(:,OP_DZ))
      temp79d = temp79a - (temp79b + temp79c)
      
      norm = norm + abs(int2(ri_79,temp79a))
      error = error + abs(int2(ri_79,temp79d))
   end do

   eqsubtract = ieqs_temp

   if(maxrank.gt.1) then
      temp1(1) = norm
      temp1(2) = error
      call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, &
           MPI_SUM, MPI_COMM_WORLD, ier)
      norm     = temp2(1)
      error    = temp2(2)
   end if
   if(myrank.eq.0) write(*,'(A,1p2E12.4)') 'Final error, norm: ', error, norm
   error = error / norm
 end subroutine calculate_gs_error
 

!======================================================================
! calc_density
! ~~~~~~~~~~~~~
!
! calculates the density as a function of the poloidal flux (and major radius if rotation is present)
!======================================================================
subroutine calc_density(psi0,pres,dens, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(6)  :: psi0, pres
  real, intent(in) :: x, z
  vectype, intent(out), dimension(6) :: dens     ! density

  real :: fbig0, fbig, fbigp, fbigpp
  real :: rbig0, rbig, rbigp, p0ni, redge
  real :: alphap0, alphap, alphapp, alphappp
  real :: r0m, r1, r1m, r2, r3, ealpha,  pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(6) :: psii     ! normalized flux

  logical :: inside_lcfs

  if(.not.inside_lcfs(psi0,x,z,.true.)) then
     call constant_field(dens,redge)
      redge = (pedge/p0)**expn
  else
     psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
     psii(2:6) = real(psi0(2:6))/(psilim - psimin)

     call fget(psii(1), fbig0, fbig, fbigp, fbigpp)
      fbig0 = fbig0 + pedge
!
!.....changed 07/05/09 scj  for inumgs.ne.0 fbig is derivative wrt total psi, not normalized
     if(inumgs.ne.0) then
        fbig = fbig*(psilim-psimin)
        fbigp = fbigp*(psilim-psimin)
     endif
!
!
  if(irot.eq.1) then
!.....include toroidal rotation in equilibrium
!...this section calculates density derivatives in presence of rotation
!   scj   11/19/10
!
!....assumes density of the form p(x,psi) = p(psi) exp (alpha(psi)*r1)
!
!...spatial derivatives of psi, not normalized psi
     pspx = real(psi0(2))
     pspy = real(psi0(3))
     pspxx= real(psi0(4))
     pspxy= real(psi0(5))
     pspyy= real(psi0(6))

     r0m = 1./rzero**2
     r1 = (x**2-rzero**2)/rzero**2
     r1m= x**2/rzero**2
     r2 = (x**2 - rzero**2)**2/rzero**4
     r3 = (x**2 - rzero**2)**3/rzero**6
     call alphaget(psii(1),alphap0,alphap,alphapp,alphappp)

!...convert all derivatives to wrt psi, not normalized psi
     fbig = fbig/(psilim-psimin)
     fbigp = fbigp/(psilim-psimin)**2
     alphap = alphap/(psilim-psimin)
     alphapp = alphapp/(psilim-psimin)**2

     p0ni = (1./p0)**expn
     rbig0 = p0ni*fbig0**expn
     rbig = 0.
     rbigp = 0.
     if(fbig0.gt.0.) then
       rbig = p0ni*expn*fbig0**(expn-1)*fbig
       rbigp= p0ni*expn*((expn-1)*fbig0**(expn-2)*fbig**2 + fbig0**(expn-1.)*fbigp)
     endif

     ealpha = exp(alphap0*r1)

     dens(1) = ealpha*rbig0

     dens(2) = ealpha*(rbig0*alphap0*2*x*r0m                                                &
                     + (rbig + rbig0*alphap*r1)*pspx)

     dens(3) = ealpha*(rbig + rbig0*alphap*r1)*pspy

     dens(4) = ealpha*(                                                                     &
             rbig0*alphap0*2*r0m + rbig0*alphap0*alphap0*4*r0m*r1m                          &
            +((2*rbig0*alphap + 2.*rbig*alphap0) + 2*rbig0*alphap0*alphap*r1)*2*x*r0m*pspx  &
            +(rbigp + (2*rbig*alphap + rbig0*alphapp)*r1 + rbig0*alphap**2*r2)*pspx*pspx    &
            +(rbig + rbig0*alphap*r1)*pspxx)

     dens(5) = ealpha*(                                                                     &
             (rbig*alphap0 + rbig0*alphap +rbig0*alphap0*alphap*r1)*2*x*r0m*pspy       &
            +(rbigp + (2.*rbig*alphap + rbig0*alphapp)*r1 + rbig0 *alphap**2*r2)*pspx*pspy  &
            +(rbig + rbig0*alphap*r1)*pspxy)

     dens(6) = ealpha*(                                                                     &
            +(rbigp + (2.*rbig*alphap + rbig0*alphapp)*r1 + rbig0 *alphap**2*r2)*pspy*pspy  &
            +(rbig + rbig0*alphap*r1)*pspyy)

  else

     if(expn .gt. 0.) then
        dens(1) = pres(1)**expn
        dens(2) = pres(1)**(expn-1.)*pres(2)*expn
        dens(3) = pres(1)**(expn-1.)*pres(3)*expn
        dens(4) = pres(1)**(expn-1.)*pres(4)*expn &
             + pres(1)**(expn-2.)*pres(2)**2.*expn*(expn-1.)
        dens(5) = pres(1)**(expn-1.)*pres(5)*expn &
             + pres(1)**(expn-2.)*pres(2)*pres(3) &
             * expn*(expn-1.)
        dens(6) = pres(1)**(expn-1.)*pres(6)*expn &
             + pres(1)**(expn-2.)*pres(3)**2.*expn*(expn-1.)
        dens = dens/p0**expn
     else   ! expn.eq.0
        dens(1) = 1.
        dens(2:6) = 0.
     endif
  endif     !....end of branch on irot
  
  endif
  
  return
end subroutine calc_density


!======================================================================
! calc_rotation
! ~~~~~~~~~~~~~
!
! calculates the rotation as a function of the poloidal flux
!======================================================================
subroutine calc_rotation(psi0,omega, x, z)
  
  use basic

  implicit none

  vectype, intent(in), dimension(6)  :: psi0
  real, intent(in) :: x, z
  vectype, intent(out), dimension(6) :: omega     ! rotation

  real :: fbig0, fbig, fbigp, fbigpp
  real :: alphap0, alphap, alphapp, alphappp
  real :: tp0,tp,tpp,p0n,befoh,befomh,befom3h
  real :: r0m, r1, r1m, r2, r3,  pspx, pspy, pspxx, pspxy, pspyy
  real, dimension(6) :: psii     ! normalized flux

  logical :: inside_lcfs

  if(.not.inside_lcfs(psi0,x,z,.true.)) then
     call constant_field(omega,0.)
  else
     psii(1) = (real(psi0(1)) - psimin)/(psilim - psimin)
     psii(2:6) = real(psi0(2:6))/(psilim - psimin)

     call fget(psii(1), fbig0, fbig, fbigp, fbigpp)
      fbig0 = fbig0 + pedge
!
!.....changed 07/05/09 scj  for inumgs.ne.0 fbig is derivative wrt total psi, not normalized
     if(inumgs.ne.0) then
        fbig = fbig*(psilim-psimin)
        fbigp = fbigp*(psilim-psimin)
     endif
!
!
  if(irot.eq.1) then
!.....include toroidal rotation in equilibrium
!...this section calculates rotation derivatives in presence of rotation
!   scj   11/19/10
!
!....assumes rotation of the form p(x,psi) = p(psi) exp (alpha(psi)*r1)
!
!...spatial derivatives of psi, not normalized psi
     pspx = real(psi0(2))
     pspy = real(psi0(3))
     pspxx= real(psi0(4))
     pspxy= real(psi0(5))
     pspyy= real(psi0(6))

     r0m = 1./rzero**2
     r1 = (x**2-rzero**2)/rzero**2
     r1m= x**2/rzero**2
     r2 = (x**2 - rzero**2)**2/rzero**4
     r3 = (x**2 - rzero**2)**3/rzero**6
     call alphaget(psii(1),alphap0,alphap,alphapp,alphappp)

!...convert all derivatives to wrt psi, not normalized psi
     fbig = fbig/(psilim-psimin)
     fbigp = fbigp/(psilim-psimin)**2
     alphap = alphap/(psilim-psimin)
     alphapp = alphapp/(psilim-psimin)**2
!...define temperature and derivatives
     p0n = 0.
     tp0 = 0.
     tp = 0.
     tpp = 0.
     befoh = 0.
     befomh = 0.
     befom3h = 0.
     omega = 0.
     if(p0 .gt. 0 .and. fbig0 .gt.0) then
        p0n = p0**expn
        tp0 = p0n*fbig0**(1.-expn)
        tp  = p0n*(1.-expn)*fbig0**(-expn)*fbig
        tpp = p0n*(expn*(expn-1.)*fbig0**(-1.-expn)*fbig**2 + (1.-expn)*fbig0**(-expn)*fbigp)

        befoh   = (2.*alphap0*tp0/rzero**2)**0.5
        befomh  = (2.*alphap0*tp0/rzero**2)**(-0.5)
        befom3h = (2.*alphap0*tp0/rzero**2)**(-1.5)

        omega(1) = befoh

        omega(2) = befomh*r0m*(alphap*tp0 + alphap0*tp)*pspx

        omega(3) = befomh*r0m*(alphap*tp0 + alphap0*tp)*pspy

        omega(4) = -befom3h*r0m**2*(alphap*tp0 + alphap0*tp)**2*pspx*pspx                   &
                 +   befomh*r0m*((alphapp*tp0 + 2.*alphap*tp + alphap0*tpp)*pspx*pspx       &
                                +(alphap*tp0 + alphap0*tp)*pspxx)

        omega(5) =  -befom3h*r0m**2*(alphap*tp0 + alphap0*tp)**2*pspx*pspy                  &
                 +   befomh*r0m*((alphapp*tp0 + 2.*alphap*tp + alphap0*tpp)*pspx*pspy       &
                                +(alphap*tp0 + alphap0*tp)*pspxy)

        omega(6) =  -befom3h*r0m**2*(alphap*tp0 + alphap0*tp)**2*pspy*pspy                  &
                 +   befomh*r0m*((alphapp*tp0 + 2.*alphap*tp + alphap0*tpp)*pspy*pspy       &
                                +(alphap*tp0 + alphap0*tp)*pspyy)
     endif
  else

     omega(1) = 0.
     omega(2) = 0.
     omega(3) = 0.
     omega(4) = 0.
     omega(5) = 0.
     omega(6) = 0.
  endif     !....end of branch on irot
  
  endif
  
  return
end subroutine calc_rotation

end module gradshafranov

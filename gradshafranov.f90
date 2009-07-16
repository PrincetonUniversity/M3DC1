module gradshafranov

  implicit none

  integer, parameter :: numvargs = 1
  
  vectype, allocatable :: psi(:)
  real, allocatable :: fun1(:), fun2(:), fun3(:), fun4(:)
  real :: dpsii
  real :: gamma2, gamma3, gamma4  
  integer :: itnum
  real :: separatrix_top, separatrix_bottom

  logical :: constraint = .false.

  integer :: npsi
  real, allocatable :: psinorm(:)
  real, allocatable :: g4big0t(:), g4bigt(:), g4bigpt(:), g4bigppt(:)
  real, allocatable :: fbig0t(:), fbigt(:), fbigpt(:), fbigppt(:)

  real :: fac2, gnorm, libetapeff

  integer, parameter :: NSTX_coils = 158
  real, dimension(NSTX_coils) :: NSTX_r, NSTX_z, NSTX_I

  integer, parameter :: ITER_coils = 18
  real, dimension(ITER_coils) :: ITER_r, ITER_z, ITER_I

  integer, parameter :: DIII_coils = 18
  real, dimension(DIII_coils) :: DIII_r, DIII_z, DIII_I


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
     separatrix_top = 1.25
     separatrix_bottom = -1.25

  case(3) ! ITER
     if(myrank.eq.0) print *, "Using standard ITER configuration"
     numcoils = ITER_coils
     xc(1:ITER_coils) = ITER_r
     zc(1:ITER_coils) = ITER_z
     ic(1:ITER_coils) = fac*ITER_I
     separatrix_top = 4.5
     separatrix_bottom = -3.5
     
  case(4) ! DIII
     if(myrank.eq.0) print *, "Using standard DIII-D configuration"
     numcoils = DIII_coils
     xc(1:DIII_coils) = DIII_r
     zc(1:DIII_coils) = DIII_z
     ic(1:DIII_coils) = fac*DIII_I
     separatrix_top = 1.2
     separatrix_bottom = -1.2
     
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
  use nintegrate_mod

  implicit none

#ifdef _AIX
  include 'mpif.h'
#endif
#include "finclude/petsc.h"
  
  real :: gsint1,gsint4,gsint2,gsint3,lhs,cfac(18)
  real, allocatable :: temp(:)
  vectype, allocatable :: b1vecini(:)
  vectype, allocatable :: b2vecini(:)
  real, dimension (2)  :: normal
  real :: curv
  integer :: izone, izonedim
  logical :: is_boundary

  integer :: itri,i,ii,i1,j,j1,jone, k,ier
  integer :: numelms, numnodes
  integer :: ibegin, iendplusone
  real :: dterm(18,18), sterm(18,18)
  real :: feedfac, fintl(-6:maxi,-6:maxi)

  real :: x, z, error, error2
  real :: sum, rhs, ajlim, curr, norm, g0, sum2, norm2
  real, dimension(5) :: temp1, temp2
  real :: psilim2,gamma2a,gamma2b,gamma3a,gamma3b
  
  real :: tstart, tend

  integer ::  idim(3)
  real :: n(2,3)
  logical :: is_edge(3)  ! is inode on boundary

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

  call numnod(numnodes)
  call numfac(numelms)

  ! allocate memory for arrays
  call createrealvec(temp, numvargs)
  call createvec(b1vecini, numvargs)
  call createvec(b2vecini, numvargs)
  call createvec(psi, numvargs)
  call createrealvec(fun1, numvargs)
  call createrealvec(fun4, numvargs)
  call createrealvec(fun2, numvargs)
  call createrealvec(fun3, numvargs)  

  call copyvec(field0, psi_g, num_fields, psi, 1, numvargs)
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
              ,"       psilim2")
  endif

  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
  mainloop: do itnum=1, iabs(igs)

     if(myrank.eq.0 .and. iprint.eq.1) print *, "GS: iteration = ", itnum
     
     ! apply boundary conditions

     if(irestart.ne.2 .or. itnum.gt.1) then
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
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_gs_solve = t_gs_solve + tend - tstart
        endif
        
        if(itnum.eq.1) then
           psi = b1vecini
        else
           psi = th_gs*b1vecini + (1.-th_gs)*psi
        endif
     endif
     
    
     ! Find new magnetic axis (extremum of psi)
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call magaxis(xmag,zmag,psi,1,numvargs,psimin,0,ier)
      if(ier .gt. 0) call safestop(27)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_gs_magaxis = t_gs_magaxis + tend - tstart
     endif
    
     ! calculate psi at the limiter
     if(ifixedb.eq.1) then
        psilim = 0.
        psilim2 = 0.
     else
        if(xlim.eq.0.) then
           call lcfs(-psi, 1, 1, 0)
           psilim = -psibound
        else
           itri = 0.
           call evaluate(xlim,zlim,psilim,ajlim,psi,1,numvargs,itri)
           
           ! calculate psi at a second limiter point as a diagnostic
           if(xlim2.gt.0) then
              itri = 0.
              call evaluate(xlim2,zlim2,psilim2,ajlim,psi,1,numvargs,itri)
           else
              psilim2 = psilim
           endif
        endif
     endif

     ! define the pressure and toroidal field functions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call fundef
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_gs_fundef = t_gs_fundef + tend - tstart
     endif


     ! Calculate error in new solution
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! (nodes shared between processes are currently overcounted.
     !  also, skip boundary nodes)
     sum = 0.
     norm = 0.
     sum2 = 0.
     norm2 = 0.
     do i=1,numnodes
        
!.....added 06/23/09...scj
       call nodcoord(i,x,z)
        call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
!       if(is_boundary) cycle
        
        call entdofs(numvargs, i, 0, ibegin, iendplusone)

        lhs = (psi(ibegin+3)-psi(ibegin+1)/x+psi(ibegin+5))/x
        rhs =  -(fun1(ibegin)+gamma2*fun2(ibegin)+                        &
             gamma3*fun3(ibegin)+gamma4*fun4(ibegin))
        if(.not. is_boundary) then
           sum = sum + (lhs-rhs)**2
           norm = norm + lhs**2
           sum2 = sum2 + abs(psi(ibegin)-b1vecini(ibegin))
           norm2 = norm2 + abs(psi(ibegin))
        endif
!!$        diff = lhs-rhs
!!$        temp1(1) = sqrt( (x-rzero)**2 + z**2)
!!$        temp1(2) = real(psi(ibegin+3)/x)
!!$        temp1(3) = real(-psi(ibegin+1)/x**2)
!!$        temp1(4) = psi(ibegin+5)/x
!!$        temp1(1) = 0.5*sqrt(psi(ibegin+1)**2 + psi(ibegin+2)**2)
!!$        temp1(2) = normal(1)**2 *psi(ibegin+5) + normal(2)**2*psi(ibegin+3) &
!!$             - 2.*normal(1)*normal(2)*psi(ibegin+4)
!!$!>>>>>debug
!!$      if(is_boundary .and. itnum.eq.iabs(igs)) then
!!$         write(70+myrank,1169) x,z,normal(1),normal(2),(real(psi(ibegin+ii)),ii=0,5), lhs, temp1(1),temp1(2),&
!!$                                (real(b2vecini(ibegin+ii)),ii=0,5)
!!$        write(80+myrank,1179) x,z,lhs,rhs
!!$      endif
!!$      if(.not.is_boundary .and. itnum.eq.iabs(igs)) then
!!$        write(90+myrank,1179) x,z,lhs,rhs
!!$      endif
!!$ 1179 format(1p4e12.4)
!!$ 1169 format(4f6.3,1p9e12.4,/,24x,1p6e12.4)
!!$!>>>>>debug end
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
     
     if(myrank.eq.0) then
        if(iprint.eq.1) then
           write(*,2002) (temp2(ii),ii=1,4)
2002       format(" temp1 array",1p4e12.4)
        endif
        write(*,1002) itnum,error,error2,xmag,psimin,psilim,psilim2
1002    format(i5,1p6e13.5)
     endif
     if(itnum .gt. 1 .and. error2 .lt. tol_gs) exit mainloop
     
     ! start of loop over triangles to compute integrals needed to keep
     !     total current and q_0 constant using gamma4, gamma2, gamma3
     gsint1 = 0.
     gsint4 = 0.
     gsint2 = 0.
     gsint3 = 0.
     curr = 0.
     do itri=1,numelms

        call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
        call calcsterm(itri, sterm, fintl)

        do j=1,18
           cfac(j) = 0.
           do k=1,20
              cfac(j) = cfac(j) + gtri(k,j,itri)*fintl(mi(k),ni(k))
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
     g0 = bzero*rzero

!.....changed 06/04/08 to allow more flexibility
!    if(numvar.eq.1 .or. nv1equ.eq.1) then
     if(nv1equ.eq.1) then
        gamma2 = 0.
        gamma3 = 0.
        gamma4 = 0.
     else
!......see if p and g functions defined numerically.  If so, only enforce total current condition
        if(constraint) then
          gamma2 = 0.
          gamma3 = 0.
          gamma4 = 1.
!!$          gamma4 = -(tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4
       else
          gamma2 =  -xmag**2*p0*p1 - 2.*abs(g0)/(xmag*q0*abs(dpsii))
          gamma3 = -4.*(abs(g0)/xmag)*djdpsi/dpsii - xmag**2*p0*p2
          gamma4 = -(tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4
        endif
     endif
         
     ! start loop over elements to define RHS vector
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
        enddo
     enddo
     call sumsharedppplvecvals(b1vecini)

  end do mainloop

  if(myrank.eq.0 ) then
     print *, "Converged GS: curr =", curr," error =",error2
     print *, "initial and final(effective) libetap", libetap, libetapeff
     gamma2a = -xmag**2*p0*p1
     gamma2b = -2.*(abs(g0)/(xmag*q0*dpsii))
     gamma3a = -4.*(abs(g0)/xmag)*djdpsi/dpsii
     gamma3b = -xmag**2*p0*p2

     write(*,1001) gamma2,gamma3,gamma4,xmag,p0,p1,p2,g0,dpsii,   &
          djdpsi,tcuro,gsint1,gsint2,gsint3,gsint4,         &
          gamma2a,gamma2b,gamma3a,gamma3b
1001 format(" gamma2,gamma3,gamma4       =",1p3e12.4,/,     &
          " xmag, p0,p1,p2,g0          =",1p5e12.4,/,     &
          " dpsii,djdpsi,tcurb         =",1p3e12.4,/,     &
          "gsint1,gint2,gsint3,gsint4  =",1p4e12.4,/,     &
          "gamm2a,gamm2b,gamm3a,gamm3b =",1p4e12.4)
  endif
!
!....if igs is positive, stop after iabs(igs) iterations, continue for igs negative
  if(itnum.eq.igs) then
     call safestop(3)
  endif

  ! populate phi0 array
  ! ~~~~~~~~~~~~~~~~~~~
  do i=1,numnodes
     !.....defines the source functions for the GS equation:
     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call assign_local_pointers(i)

     psi0_l = psi(ibegin:ibegin+5)
     psi1_l = 0.

     ! temp = Psi (normalized flux)
     dpsii = 1./(psilim - psimin)
     temp(ibegin) = (psi(ibegin) - psimin)*dpsii
     temp(ibegin+1:ibegin+5) = psi(ibegin+1:ibegin+5)*dpsii

     call calc_toroidal_field(temp(ibegin:ibegin+5), bz0_l)

     call calc_pressure(temp(ibegin:ibegin+5),p0_l)

     call nodcoord(i, x, z)

     if(linear.eq.0) then
        if(((z.gt.separatrix_top) .or. (z .lt.separatrix_bottom)) .and. &
             (separatrix_top.ne.0. .or. separatrix_bottom.ne.0.)) then
           !        if(myrank.eq.0 .and. iprint.ge.1) &
           !             print *, 'removing points above ', separatrix_top, &
           !             ' and below', separatrix_bottom
           p0_l(1) = pedge
           p0_l(2:6) = 0.
        endif
     end if

     pe0_l = (1. - ipres*pi0/p0)*p0_l

     if(expn .gt. 0.) then
        den0_l(1) = p0_l(1)**expn
        den0_l(2) = p0_l(1)**(expn-1.)*p0_l(2)*expn
        den0_l(3) = p0_l(1)**(expn-1.)*p0_l(3)*expn
        den0_l(4) = p0_l(1)**(expn-1.)*p0_l(4)*expn &
             + p0_l(1)**(expn-2.)*p0_l(2)**2.*expn*(expn-1.)
        den0_l(5) = p0_l(1)**(expn-1.)*p0_l(5)*expn &
             + p0_l(1)**(expn-2.)*p0_l(2)*p0_l(3) &
             * expn*(expn-1.)
        den0_l(6) = p0_l(1)**(expn-1.)*p0_l(6)*expn &
             + p0_l(1)**(expn-2.)*p0_l(3)**2.*expn*(expn-1.)
        den0_l = den0_l/p0**expn
      else   ! expn.eq.0
        den0_l(1) = 1.
        den0_l(2:6) = 0.
      endif

  end do
      psibound = psilim

  ! free memory
  call deleterealvec(temp)
  call deletevec(b1vecini)
  call deletevec(b2vecini)
  call deletevec(psi)
  call deleterealvec(fun1)
  call deleterealvec(fun4)
  call deleterealvec(fun2)
  call deleterealvec(fun3)

  return

end subroutine gradshafranov_solve

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
           dum(index) = dum(index) + sum*val
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

  dpsii = 1./(psilim - psimin)

  call numnod(numnodes)
  do l=1,numnodes

     call nodcoord(l, x, z)

     call entdofs(numvargs, l, 0, ibegin, iendplusone)
     pso =  (psi(ibegin)-psimin)*dpsii

!
!......the following line performs more reliably than: if(pso .lt. 0 .or. pso .gt. 1.) then
     if(pso .gt. 1.) then
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

        fun1(ibegin)   = x*fbig
        fun1(ibegin+1) = fbig + x*fbigp*psox
        fun1(ibegin+2) =        x*fbigp*psoy
        fun1(ibegin+3) = 2.*fbigp*psox + x*(fbigpp*psox**2+fbigp*psoxx)
        fun1(ibegin+4) = fbigp*psoy + x*(fbigpp*psox*psoy +fbigp*psoxy)
        fun1(ibegin+5) = x*(fbigpp*psoy**2 + fbigp*psoyy)

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


!======================================================================
! calc_toroidal_field
! ~~~~~~~~~~~~~~~~~~~
!
! calculates the toroidal field (I) as a function of the 
! normalized flux
!======================================================================
subroutine calc_toroidal_field(psii,tf)
  use basic

  real, intent(in), dimension(6)  :: psii     ! normalized flux
  vectype, intent(out), dimension(6) :: tf    ! toroidal field (I)

  real, dimension(6) :: g2, g3, g4
  real :: g4big0, g4big, g4bigp, g4bigpp
  real :: g2big, g2bigp, g3big, g3bigp
  vectype :: g0

  if(psii(1) .gt. 1.) then
     g0 = bzero*rzero
     call constant_field(tf, g0)
  else
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
! calculates the pressure as a function of the normalized flux
!======================================================================
subroutine calc_pressure(psii,pres)
  
  use basic

  implicit none

  real, intent(in), dimension(6)  :: psii     ! normalized flux
  real :: fbig0, fbig, fbigp, fbigpp
  vectype, intent(out), dimension(6) :: pres     ! pressure

  if(psii(1) .gt. 1.) then
     pres = 0.
  else
     call fget(psii(1), fbig0, fbig, fbigp, fbigpp)
!
!.....changed 07/05/09 scj  for inumgs.ne.0 fbig is derivative wrt total psi, not normalized
      if(inumgs.ne.0) then
        fbig = fbig/dpsii
        fbigp = fbigp/dpsii
      endif
!
     
     pres(1) = fbig0
     pres(2) = psii(2)*fbig
     pres(3) = psii(3)*fbig
     pres(4) = (psii(4)*fbig + psii(2)**2*fbigp)
     pres(5) = (psii(5)*fbig + psii(2)*psii(3)*fbigp)
     pres(6) = (psii(6)*fbig + psii(3)**2*fbigp)
  endif

  pres = pres
  pres(1) = pres(1) + pedge

  return
end subroutine calc_pressure

subroutine readpgfiles
  use basic

  implicit none

  integer :: j

  print *, "Reading profiles files"

  open(unit=76,file="profiles-p",status="old")
  read(76,803) npsi
  allocate(psinorm(npsi))
  allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
  do j=1,npsi
    read(76,802) psinorm(j),fbig0t(j),fbigt(j),fbigpt(j),fbigppt(j)
  enddo
  close(76)

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

 subroutine default_profiles
   
   use basic

   implicit none

   integer :: j
   real :: psii

   if(myrank.eq.0) print *, "Using default p and I profiles"

   npsi = 500
   allocate(psinorm(npsi))
   allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
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
 subroutine create_profile(n, p, pp, f, ffp, flux)
   implicit none

   integer, intent(in) :: n
   real, dimension(n), intent(in) :: p, pp, f, ffp, flux

   real, dimension(4) :: a
   real :: dp,dpp
   integer :: j

   npsi = n

   allocate(psinorm(npsi))
   allocate(fbig0t(npsi),fbigt(npsi),fbigpt(npsi),fbigppt(npsi))
   allocate(g4big0t(npsi),g4bigt(npsi),g4bigpt(npsi),g4bigppt(npsi))

   fbig0t = p                                 ! p
   fbigt = pp * (flux(npsi) - flux(1))        ! p' = dp/dPsi
   g4big0t = 0.5*(f**2 - f(npsi)**2)          ! g
   g4bigt = ffp * (flux(npsi) - flux(1))      ! f f' = f * df/dPsi

   do j=1,npsi 
      call cubic_interpolation_coeffs(flux,npsi,j,a)
      psinorm(j) = (flux(j) - flux(1)) / (flux(npsi) - flux(1))
      dp = a(2)      ! d psi/dn
      dpp = 2.*a(3)  ! d^2 psi/dn^2

      call cubic_interpolation_coeffs(fbigt,npsi,j,a)
      fbigpt(j) =     a(2)/dp                       ! p''
      fbigppt(j) = 2.*a(3)/dp**2 - a(2)*dpp/dp**3   ! p'''

      call cubic_interpolation_coeffs(g4bigt,npsi,j,a)
      g4bigpt(j)  =    a(2)/dp                      ! (f f')'
      g4bigppt(j) = 2.*a(3)/dp**2 - a(2)*dpp/dp**3  ! (f f')''
   end do

   fbigpt = fbigpt*(flux(npsi) - flux(1))
   fbigppt = fbigppt*(flux(npsi) - flux(1))**2
   g4bigpt = g4bigpt*(flux(npsi) - flux(1))
   g4bigppt = g4bigppt*(flux(npsi) - flux(1))**2

   constraint = .true.

 end subroutine create_profile

end module gradshafranov

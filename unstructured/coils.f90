module coils
  implicit none

  integer, parameter :: maxfilaments = 40000
  integer, parameter :: maxcoils = 2000

contains

  !======================================================
  ! load_coils
  ! ~~~~~~~~~~
  !
  ! reads coil and current data from file
  !======================================================
 subroutine load_coils(xc, zc, ic, numcoils, coil_filename, current_filename, &
      coil_mask, filaments, numcoils_base)
   use math
   use read_ascii

   implicit none

   include 'mpif.h'

   real, intent(out), dimension(maxfilaments) :: xc, zc  ! coordinates of each coil
   complex, intent(out), dimension(maxfilaments) :: ic   ! current in each coil
   integer, intent(out) :: numcoils                  ! number of coils read
   character*(*) :: coil_filename, current_filename  ! input files
   integer, intent(out), dimension(maxfilaments), optional :: coil_mask
   integer, intent(out), dimension(maxfilaments), optional :: filaments
   integer, intent(out), optional :: numcoils_base

   integer, parameter :: fcoil = 34

   integer :: i, j, k, s, ier, rank
   integer :: coil_style, ncur, ncoil
   real :: dx, dz, a1, a2

   real, allocatable :: xv(:), zv(:), wv(:), hv(:), a1v(:), a2v(:), &
        cv(:), phasev(:), turnsv(:)
   integer, allocatable :: sxv(:), syv(:)

   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ier)

   xc = 0.
   zc = 0.
   ic = 0.
   numcoils = 0
   ncur = 0
   ncoil = 0

   ! determine formatting of coil file
   if(rank.eq.0) then 
      open(unit=fcoil, file=coil_filename, status='old', action='read')
      read(fcoil, *, err=11) dx
      coil_style=0
      goto 20
11    coil_style=1
20    close(fcoil)
      print *, 'Coil file style = ', coil_style
   end if
   call MPI_bcast(coil_style, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)


   if(coil_style.eq.0) then
      call read_ascii_column(coil_filename, xv,  ncoil, 0, 1)
      call read_ascii_column(coil_filename, zv,  ncoil, 0, 2)
      call read_ascii_column(coil_filename, wv,  ncoil, 0, 3)
      call read_ascii_column(coil_filename, hv,  ncoil, 0, 4)
      call read_ascii_column(coil_filename, a1v, ncoil, 0, 5)
      call read_ascii_column(coil_filename, a2v, ncoil, 0, 6)
      call read_ascii_column(coil_filename, sxv, ncoil, 0, 7)
      call read_ascii_column(coil_filename, syv, ncoil, 0, 8)
   else if(coil_style.eq.1) then
      call read_ascii_column(coil_filename, xv,     ncoil, 0, 4)
      call read_ascii_column(coil_filename, zv,     ncoil, 0, 5)
      call read_ascii_column(coil_filename, wv,     ncoil, 0, 6)
      call read_ascii_column(coil_filename, hv,     ncoil, 0, 7)
      call read_ascii_column(coil_filename, a1v,    ncoil, 0, 8)
      call read_ascii_column(coil_filename, a2v,    ncoil, 0, 9)
      call read_ascii_column(coil_filename, turnsv, ncoil, 0, 10)

      allocate(sxv(ncoil), syv(ncoil))
      do i=1, ncoil
         sxv(i) = nint(wv(i)/0.01)
         syv(i) = nint(hv(i)/0.01)
         if(sxv(i).lt.1) sxv(i) = 1
         if(syv(i).lt.1) syv(i) = 1
      end do
      deallocate(turnsv)
   end if

   if(rank.eq.0) print *, coil_filename, ' ', current_filename
   call read_ascii_column(current_filename, cv,     ncur, 0, 1)
   call read_ascii_column(current_filename, phasev, ncur, 0, 2)
   if(ncur.ne.ncoil) then
      if(rank.eq.0) print *, "Error: number of coils != number of currents"
      call safestop(301)
   end if
   if(ncoil.gt.maxcoils) then
      if(rank.eq.0) print *, "Error: too many coils ", ncoil
      call safestop(301)
   end if

   do i=1, ncoil
      if(sxv(i).lt.1) sxv(i) = 1
      if(syv(i).lt.1) syv(i) = 1
   end do

   if(present(numcoils_base)) numcoils_base = ncoil

   cv = amu0 * 1000. * cv / (sxv*syv) / twopi

   ! divide coils into sub-coils
   s = 0
   do i=1, ncoil
      a1 = a1v(i)*pi/180.
      a2 = a2v(i)*pi/180.
      if(a2.ne.0.) a2 = a2 + pi/2.
      a1 = tan(a1)
      a2 = tan(a2)

      if(rank.eq.0) then 
         print *, 'Coil ', i
         write(*,'(6F12.4,2I5)') xv(i), zv(i), wv(i), hv(i), a1, a2, sxv(i), syv(i)
         write(*,'(A,2F12.4)') "current (kA), phase (deg): ", cv(i), phasev(i)
      end if

      do j=1, sxv(i)
         do k=1, syv(i)
            s = s + 1
            
            if(s.gt.maxfilaments) then
               print *, 'Too many filaments.', s, maxfilaments
               call safestop(301)
            end if
            if(sxv(i).eq.1) then 
               dx = 0.
            else
               dx = wv(i)*(j-1.)/(sxv(i)-1.) - wv(i)/2.
            end if
            if(syv(i).eq.1) then 
               dz = 0.
            else
               dz = hv(i)*(k-1.)/(syv(i)-1.) - hv(i)/2.
            end if

            xc(s) = xv(i) + dx - dz*a2
            zc(s) = zv(i) + dz + dx*a1
            
            ic(s) = cv(i)*(cos(pi*phasev(i)/180.) &
                 - (0.,1.)*sin(pi*phasev(i)/180.))
            
            if(present(coil_mask)) coil_mask(s) = i
            if(present(filaments)) filaments(s) = sxv(i)*syv(i)
         end do
      end do
   end do

100   numcoils = s
   deallocate(xv, zv, wv, hv, a1v, a2v, sxv, syv)

 end subroutine load_coils

 subroutine field_from_coils(xc, zc, ic, nc, f, ipole, ierr, ncs, fs, mask)
   use basic
   use field
   use mesh_mod
   use m3dc1_nint
   use newvar_mod

   implicit none

   include 'mpif.h'

   real, intent(in), dimension(nc) :: xc, zc     ! array of coil positions
   complex, intent(in), dimension(nc) :: ic      ! array of coil currents
   integer, intent(in) :: nc                     ! number of coils
   type(field_type), intent(inout) :: f          ! poloidal flux field
   integer, intent(in) :: ipole                  ! type of field to add
   integer, intent(out) :: ierr
   ! poloidal flux field for each coil
   integer, intent(in), optional :: ncs ! number of real coils
   type(field_type), intent(in), dimension(ncs), optional :: fs(:) ! fields
   integer, intent(in), dimension(nc), optional :: mask  !filaments to coils
   
   integer :: numelms, k, itri
   type(field_type) :: psi_vec
   vectype, dimension(dofs_per_element) :: dofs
   real, allocatable :: g(:,:)
   real :: cur

   integer :: inode, numnodes, icounter_t
   vectype, dimension(dofs_per_node) :: datain, dataout

   if(nc.le.0) return

   ierr = 0

   if(present(ncs)) then
      if(.not.present(fs) .or. .not.present(mask)) then
         if(myrank.eq.0) then
            print *, "ERROR: field_from_coils missing some optional arguments"
         end if
         call safestop(1)
      end if
      numnodes = owned_nodes()
      do icounter_t=1,numnodes
         inode = nodes_owned(icounter_t)
         call get_node_data(f, inode, dataout)
         do k=1, ncs
            cur = sum(ic, mask=(mask.eq.k))
            call get_node_data(fs(k), inode, datain)
            datain = datain*cur
            dataout = dataout + datain
         end do
         call set_node_data(f, inode, dataout)
      enddo
      call finalize(f%vec)

      return
   end if

   allocate(g(MAX_PTS,nc))
   call create_field(psi_vec)
   psi_vec = 0.

   numelms = local_elements()
   do itri=1,numelms
      temp79a = 0.

      ! calculate field due to coil currents
      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
      call define_fields(itri,0,1,0)
      call gvect0(x_79,z_79,npoints,xc,zc,nc,g(1:npoints,:),ipole,ierr)
      do k=1, nc
         temp79a = temp79a - g(:,k)*ic(k)
      end do

      dofs = intx2(mu79(:,:,OP_1),temp79a)

      call vector_insert_block(psi_vec%vec, itri, 1, dofs, MAT_ADD)
   enddo

   if(allocated(g)) deallocate(g)

   call sum_shared(psi_vec%vec)
   call newsolve(mass_mat_lhs%mat,psi_vec%vec,ierr)
   call add(f, psi_vec)
   call destroy_field(psi_vec)

 end subroutine field_from_coils

 subroutine store_field_from_coils(xc, zc, nc, nc_base, fs, mask, ipole, ierr)
   use basic
   use field
   use mesh_mod
   use m3dc1_nint
   use newvar_mod

   implicit none

   include 'mpif.h'

   real, intent(in), dimension(nc) :: xc, zc   ! array of coil positions
   integer, intent(in) :: nc                   ! number of filaments
   integer, intent(in) :: nc_base              ! number of coils
   integer, dimension(maxfilaments), intent(in) :: mask ! relate filaments to coils
   type(field_type), allocatable, intent(inout) :: fs(:)       ! poloidal flux field for each coil
   integer, intent(in) :: ipole                ! type of field to add
   integer, intent(out) :: ierr

   integer :: numelms, k, j, itri, nfil
   vectype, dimension(dofs_per_element) :: dofs
   real, allocatable :: g(:,:)

   if(allocated(fs)) then
      if (myrank.eq.0) print *,"ERROR: store_field_from_coils called twice somehow"
      call safestop(1)
   end if

   if(nc.le.0) return

   ierr = 0

   allocate(fs(nc_base))
   allocate(g(MAX_PTS,nc))

   do j=1,nc_base
      call create_field(fs(j))
   enddo
   
   numelms = local_elements()
   do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       ! Field due to coil currents
       call gvect0(x_79,z_79,npoints,xc,zc,nc,g(1:npoints,:),ipole,ierr)
       do j=1, nc_base
          temp79a = 0.
          nfil = 0
          do k=1, nc
             if (mask(k).eq.j) then
                temp79a = temp79a - g(:,k)
                nfil = nfil + 1
             end if
          end do
          temp79a = temp79a/nfil
          dofs = intx2(mu79(:,:,OP_1), temp79a)
          call vector_insert_block(fs(j)%vec, itri, 1, dofs, MAT_ADD)
       enddo
    enddo

    deallocate(g)

    do j=1,nc_base
       call sum_shared(fs(j)%vec)
       call newsolve(mass_mat_lhs%mat,fs(j)%vec,ierr)
    enddo

  end subroutine store_field_from_coils

!============================================================
subroutine gvect(r,z,npt,xi,zi,n,g,nmult,ierr)
  ! calculates derivatives wrt first argument
  use math

  implicit none

  integer, intent(in) :: n, npt, nmult
  integer, intent(out) :: ierr
  real, dimension(npt), intent(in) :: r, z
  real, dimension(n), intent(in) :: xi, zi
  real, dimension(npt,n,6), intent(out) :: g

  real, parameter :: eps = 1e-6
  
  real :: a0,a1,a2,a3,a4
  real :: b0,b1,b2,b3,b4
  real :: c1,c2,c3,c4
  real :: d1,d2,d3,d4

  real, dimension(npt) :: rpxi, rxi, zmzi, rksq, rk, sqrxi, x 
  real, dimension(npt) :: ce, ck, term1, term2, rz, co

  integer :: i, imult

  data a0,a1,a2,a3,a4/1.38629436112,9.666344259e-2,                 &
       3.590092383e-2,3.742563713e-2,1.451196212e-2/
  data b0,b1,b2,b3,b4/.5,.12498593597,6.880248576e-2,               &
       3.328355346e-2,4.41787012e-3/
  data c1,c2,c3,c4/.44325141463,6.260601220e-2,                     &
       4.757383546e-2,1.736506451e-2/
  data d1,d2,d3,d4/.24998368310,9.200180037e-2,                     &
       4.069697526e-2,5.26449639e-3/

  ierr = 0
  if(nmult.le.0) then
     do i=1,n
        rpxi=r+xi(i)
        rxi=r*xi(i)
        zmzi=z-zi(i)
        rksq=4.*rxi/(rpxi**2+zmzi**2)
        rk=sqrt(rksq)
        sqrxi=sqrt(rxi)
        x=1.-rksq + eps

        ce=1.+x*(c1+x*(c2+x*(c3+x*c4)))+                               &
             x*(d1+x*(d2+x*(d3+x*d4)))*(-alog(x))
        ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                               &
             (b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-alog(x))

        term1=2.*ck-2.*ce-ce*rksq/x
        term2=2.*xi(i)-rksq*rpxi
           
        g(:,i,1) =- sqrxi*(2.*ck-2.*ce-ck*rksq)/rk
        g(:,i,2)=-rk*0.25/sqrxi*(rpxi*term1                            &
             +2.*xi(i)*(ce/x-ck))
        g(:,i,3)=-rk*0.25*zmzi/sqrxi*term1
        g(:,i,5)=0.0625*zmzi*(rk/sqrxi)**3*(rpxi*term1                 &
             +(ce-ck+2.*ce*rksq/x)*                                    &
             (term2)/x)
        g(:,i,6)=-rk*0.25/sqrxi*(term1*                                &
             (1.-rksq*zmzi**2/(4.*rxi))+zmzi**2*rksq**2/(4.*rxi*x)     &
             *(ce-ck+2.*ce*rksq/x))
        g(:,i,4)=-rk*0.25/sqrxi*(-rksq*rpxi/(4.*rxi)*                  &
             (rpxi*term1+2.*xi(i)*(ce/x-ck))+term1-                    &
             rksq*rpxi/(4.*rxi*x)*(ce-ck+2.*ce*rksq/x)*                &
             (term2)+rksq/(2.*r*x)*(2.*ce/x-ck)*term2)
     end do
     return
  endif
     
  ! check for multipolar coils
  do i=1,n
     if(xi(i) .lt. 100.) cycle
     rz = zi(i)
     imult = int(xi(i) - 100.)
     if(imult .lt. 0 .or. imult.gt.10) then
        ! error
        ierr=39
        return
     endif

     select case(imult)
     case(0)
        ! even nullapole
        g(:,i,1) = twopi*rz**2
        g(:,i,2) = 0.
        g(:,i,3) = 0.
        g(:,i,4) = 0.
        g(:,i,5) = 0.
        g(:,i,6) = 0.
        
     case(1)
        ! odd nullapole
        g(:,i,1) = 0.
        g(:,i,2) = 0.
        g(:,i,3) = 0.
        g(:,i,4) = 0.
        g(:,i,5) = 0.
        g(:,i,6) = 0.
        
     case(2)
        ! even dipole
        g(:,i,1) = twopi*(r**2 - rz**2)/2.
        g(:,i,2) = twopi*r
        g(:,i,3) = 0.
        g(:,i,4) = 0.
        g(:,i,5) = 0.
        g(:,i,6) = twopi
        
     case(3)
        ! odd dipole
        co=twopi/rz
        g(:,i,1) = co*(r**2*z)
        g(:,i,2) = co*(2.*r*z)
        g(:,i,3) = co*(r**2)
        g(:,i,5) = co*2.*r
        g(:,i,6) = 0.
        g(:,i,4) = co*2*z
        
     case(4)
        ! even quadrapole
        co=pi/(4.*rz**2)
        g(:,i,1) = co*(r**4-4.*r**2*z**2 - 2.*r**2*rz**2+rz**4)
        g(:,i,2) = co*(4.*r**3-8.*r*z**2-4.*r*rz**2)
        g(:,i,3) = co*(-8.*r**2*z)
        g(:,i,5) = co*(-16.*r*z)
        g(:,i,6) = co*(-8.*r**2)
        g(:,i,4) = co*(12.*r**2 - 8.*z**2 - 4.*rz**2)
        
     case(5)
        ! odd quadrapole
        co=pi/(3.*rz**3)
        g(:,i,1) = co*r**2*z*(3.*r**2-4.*z**2-3.*rz**2)
        g(:,i,2) = co*(12.*r**3*z-8.*r*z**3-6.*r*z*rz**2)
        g(:,i,3) = co*(3.*r**4 - 12.*r**2*z**2 - 3.*r**2*rz**2)
        g(:,i,5) = co*(12.*r**3-24.*r*z**2 - 6.*r*rz**2)
        g(:,i,6) = co*(-24.*r**2*z)
        g(:,i,4) = co*(36.*r**2*z-8.*z**3-6.*z*rz**2)

     case(6)
        ! even hexapole
        co=pi/(12.*rz**4)
        g(:,i,1) = co*(r**6 - 12.*r**4*z**2 - 3.*r**4*rz**2       &
             + 8.*r**2*z**4 + 12.*r**2*z**2*rz**2           &
             + 3.*r**2*rz**4 - rz**6 )
        g(:,i,2)= co*(6.*r**5 - 48.*r**3*z**2                       &
             - 12.*r**3*rz**2 + 16.*r*z**4                     &
             + 24.*r*z**2*rz**2 + 6.*r*rz**4 )
        g(:,i,3)= co*(-24.*r**4*z + 32.*r**2*z**3                &
             + 24.*r**2*z*rz**2)
        g(:,i,5)=co*(-96.*r**3*z+64.*r*z**3                     &
             + 48.*r*z*rz**2 )
        g(:,i,4)=co*(30.*r**4-144.*r**2*z**2-36.*r**2*rz**2     &
             + 16.*z**4 + 24.*z**2*rz**2 + 6.*rz**4)
        g(:,i,6)=co*(-24.*r**4 + 96.*r**2*z**2 + 24.*r**2*rz**2)

     case(7)
        ! odd hexapole
        co=pi/(30.*rz**5)
        g(:,i,1) = co*(15.*r**6*z - 60.*r**4*z**3                 &
             - 30.*r**4*z*rz**2 + 24.*r**2*z**5             &
             + 40.*r**2*z**3*rz**2 + 15.*r**2*z*rz**4)
        g(:,i,2)= co*(90.*r**5*z - 240.*r**3*z**3                &
             - 120.*r**3*z*rz**2 + 48.*r*z**5               &
             + 80.*r*z**3*rz**2 + 30.*r*z*rz**4)
        g(:,i,3)= co*(15.*r**6 - 180.*r**4*z**2                     &
             - 30.*r**4*rz**2 + 120.*r**2*z**4                 &
             +120.*r**2*z**2*rz**2 + 15.*r**2*rz**4)
        g(:,i,5)=co*(90.*r**5 - 720.*r**3*z**2                     &
             - 120.*r**3*rz**2 + 240.*r*z**4                   &
             + 240.*r*z**2*rz**2 + 30.*r*rz**4)
        g(:,i,6)=co*(-360.*r**4*z + 480.*r**2*z**3              &
             + 240.*r**2*z*rz**2)
        g(:,i,4)=co*(450.*r**4*z - 720.*r**2*z**3               &
             - 360.*r**2*z*rz**2 + 48.*z**5                    &
             + 80.*z**3*rz**2 + 30.*z*rz**4)

     case(8)
        ! even octapole
        co=pi/(160.*rz**6)
        g(:,i,1) = co*(5.*r**8 - 120.*r**6*z**2 - 20.*r**6*rz**2  &
             + 240.*r**4*z**4 + 240.*r**4*z**2*rz**2        &
             + 30.*r**4*rz**4 - 64.*r**2*z**6                  &
             - 160.*r**2*z**4*rz**2 - 120.*r**2*z**2*rz**4  &
             - 20.*r**2*rz**6 + 5.*rz**8)
        g(:,i,2)= co*(40.*r**7 - 720.*r**5*z**2 - 120.*r**5*rz**2 &
             + 960.*r**3*z**4 + 960.*r**3*z**2*rz**2        &
             + 120.*r**3*rz**4 - 128.*r*z**6                   &
             - 320.*r*z**4*rz**2 - 240.*r*z**2*rz**4        &
             - 40.*r*rz**6)
        g(:,i,3)= co*(-240.*r**6*z + 960.*r**4*z**3              &
             + 480.*r**4*z*rz**2 - 384.*r**2*z**5           &
             - 640.*r**2*z**3*rz**2 - 240.*r**2*z*rz**4)
        g(:,i,5)=co*(-1440.*r**5*z + 3840.*r**3*z**3            &
             + 1920.*r**3*z*rz**2 - 768.*r*z**5             &
             - 1280.*r*z**3*rz**2 - 480.*r*z*rz**4)
        g(:,i,6)=co*(-240.*r**6 + 2880.*r**4*z**2                  &
             + 480.*r**4*rz**2 - 1920.*r**2*z**4               &
             - 1920.*r**2*z**2*rz**2 - 240.*r**2*rz**4)
        g(:,i,4)=co*(280.*r**6 - 3600.*r**4*z**2                   &
             - 600.*r**4*rz**2 + 2880.*r**2*z**4               &
             + 2880.*r**2*z**2*rz**2                              &
             + 360.*r**2*rz**4 - 128.*z**6                        &
             - 320.*z**4*rz**2 - 240.*z**2*rz**4 - 40.*rz**6)

     case(9)
        ! odd octapole
        co=pi/(140.*rz**7)
        g(:,i,1) = co*r**2*z*(35.*r**6 - 280.*r**4*z**2        &
             - 105.*r**4*rz**2 + 336.*r**2*z**4                &
             + 420.*r**2*z**2*rz**2 + 105.*r**2*rz**4          &
             - 64.*z**6 - 168.*z**4*rz**2                         &
             - 140.*z**2*rz**4 - 35.*rz**6)
        g(:,i,2)= co*(280.*r**7*z - 1680.*r**5*z**3              &
             - 630.*r**5*z*rz**2 + 1344.*r**3*z**5          &
             + 1680.*r**3*z**3*rz**2 + 420.*r**3*z*rz**4    &
             - 128.*r*z**7 - 336.*r*z**5*rz**2              &
             - 280.*r*z**3*rz**4 - 70.*r*z*rz**6)
        g(:,i,3)= co*(35.*r**8-840.*r**6*z**2-105.*r**6*rz**2    &
             + 1680.*r**4*z**4 + 1260.*r**4*z**2*rz**2      &
             + 105.*r**4*rz**4 - 448.*r**2*z**6                &
             - 840.*r**2*z**4*rz**2 - 420.*r**2*z**2*rz**4  &
             - 35.*r**2*rz**6)
        g(:,i,5)=co*(280.*r**7 - 5040.*r**5*z**2                   &
             - 630.*r**5*rz**2 + 6720.*r**3*z**4               &
             + 5040.*r**3*z**2*rz**2 + 420.*r**3*rz**4         &
             - 896.*r*z**6 - 1680.*r*z**4*rz**2             &
             - 840.*r*z**2*rz**4 - 70.*r*rz**6)
        g(:,i,6)=co*(-1680.*r**6*z + 6720.*r**4*z**3            &
             + 2520.*r**4*z*rz**2 - 2688.*r**2*z**5         &
             - 3360.*r**2*z**3*rz**2 - 840.*r**2*z*rz**4)
        g(:,i,4)=co*(1960.*r**6*z - 8400.*r**4*z**3             &
             - 3150.*r**4*z*rz**2 + 4032.*r**2*z**5         &
             + 5040.*r**2*z**3*rz**2 + 1260.*r**2*z*rz**4   &
             - 128.*z**7 - 336.*z**5*rz**2                        &
             - 280.*z**3*rz**4 - 70.*z*rz**6)

     case(10)
        ! even decapole
        co=pi/(560.*rz**8)
        g(:,i,1) = co*(7.*r**10 - 280.*r**8*z**2 - 35.*r**8*rz**2 &
             + 1120.*r**6*z**4 + 840.*r**6*z**2*rz**2       &
             + 70.*r**6*rz**4 - 896.*r**4*z**6                 &
             - 1680.*r**4*z**4*rz**2 - 840.*r**4*z**2*rz**4 &
             - 70.*r**4*rz**6 + 128.*r**2*z**8                 &
             + 448.*r**2*z**6*rz**2 + 560.*r**2*z**4*rz**4  &
             + 280.*r**2*z**2*rz**6 + 35.*r**2*rz**8           &
             - 7.*rz**10)
        g(:,i,2)= co*(70.*r**9 - 2240.*r**7*z**2                    &
             - 280.*r**7*rz**2 + 6720.*r**5*z**4               &
             + 5040.*r**5*z**2*rz**2 + 420.*r**5*rz**4         &
             - 3584.*r**3*z**6 - 6720.*r**3*z**4*rz**2      &
             - 3360.*r**3*z**2*rz**4 - 280.*r**3*rz**6         &
             + 256.*r*z**8 + 896.*r*z**6*rz**2              &
             + 1120.*r*z**4*rz**4 + 560.*r*z**2*rz**6       &
             + 70.*r*rz**8)
        g(:,i,3)= co*(-560.*r**8*z + 4480.*r**6*z**3             &
             + 1680.*r**6*z*rz**2 - 5376.*r**4*z**5         &
             - 6720.*r**4*z**3*rz**2 - 1680.*r**4*z*rz**4   &
             + 1024.*r**2*z**7 + 2688.*r**2*z**5*rz**2      &
             + 2240.*r**2*z**3*rz**4 + 560.*r**2*z*rz**6)
        g(:,i,5)=co*(-4480.*r**7*z + 26880.*r**5*z**3           &
             + 10080.*r**5*z*rz**2 - 21504.*r**3*z**5       &
             - 26880.*r**3*z**3*rz**2 - 6720.*r**3*z*rz**4  &
             + 2048.*r*z**7 + 5376.*r*z**5*rz**2            &
             + 4480.*r*z**3*rz**4 + 1120.*r*z*rz**6)
        g(:,i,6)=co*(-560.*r**8 + 13440.*r**6*z**2                 &
             + 1680.*r**6*rz**2 - 26880.*r**4*z**4             &
             - 20160.*r**4*z**2*rz**2 - 1680.*r**4*rz**4       &
             + 7168.*r**2*z**6 + 13440.*r**2*z**4*rz**2     &
             + 6720.*r**2*z**2*rz**4 + 560.*r**2*rz**6)
        g(:,i,4)=co*(630.*r**8 - 15680.*r**6*z**2                  &
             - 1960.*r**6*rz**2 + 33600*r**4*z**4              &
             + 25200.*r**4*z**2*rz**2 + 2100.*r**4*rz**4       &
             - 10752.*r**2*z**6 - 20160.*r**2*z**4*rz**2    &
             - 10080.*r**2*z**2*rz**4 - 840.*r**2*rz**6        &
             + 256.*z**8 + 896.*z**6*rz**2                        &
             + 1120.*z**4*rz**4 + 560.*z**2*rz**6                 &
             + 70.*rz**8)

     case default
        print *, 'Error: unknown multipole ', imult
     end select
  end do
end subroutine gvect

!============================================================
! gvect0
! ~~~~~~
! This is same as gvect for case when
! only OP_1 element of g is needed
!============================================================
subroutine gvect0(r,z,npt,xi,zi,n,g,nmult,ierr)
  use math

  implicit none

  integer, intent(in) :: n, npt, nmult
  integer, intent(out) :: ierr
  real, dimension(npt), intent(in) :: r, z
  real, dimension(n), intent(in) :: xi, zi
  real, dimension(npt,n), intent(out) :: g
  
  real :: a0,a1,a2,a3,a4
  real :: b0,b1,b2,b3,b4
  real :: c1,c2,c3,c4
  real :: d1,d2,d3,d4

  real, dimension(npt) :: rpxi, rxi, zmzi, rksq, rk, sqrxi, x 
  real, dimension(npt) :: ce, ck, rz, co

  integer :: i, imult

  data a0,a1,a2,a3,a4/1.38629436112,9.666344259e-2,                 &
       3.590092383e-2,3.742563713e-2,1.451196212e-2/
  data b0,b1,b2,b3,b4/.5,.12498593597,6.880248576e-2,               &
       3.328355346e-2,4.41787012e-3/
  data c1,c2,c3,c4/.44325141463,6.260601220e-2,                     &
       4.757383546e-2,1.736506451e-2/
  data d1,d2,d3,d4/.24998368310,9.200180037e-2,                     &
       4.069697526e-2,5.26449639e-3/

  ierr = 0
  if(nmult.le.0) then
     do i=1,n
        rpxi=r+xi(i)
        rxi=r*xi(i)
        zmzi=z-zi(i)
        rksq=4.*rxi/(rpxi**2+zmzi**2)
        rk=sqrt(rksq)
        sqrxi=sqrt(rxi)
        x=1.-rksq

        ce=1.+x*(c1+x*(c2+x*(c3+x*c4)))+                               &
             x*(d1+x*(d2+x*(d3+x*d4)))*(-alog(x))
        ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                               &
             (b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-alog(x))

        g(:,i) = -sqrxi*(2.*ck-2.*ce-ck*rksq)/rk
     end do
     return
  endif
     
  ! check for multipolar coils
  do i=1,n
     if(xi(i) .lt. 100.) cycle
     rz = zi(i)
     imult = int(xi(i) - 100.)
     if(imult .lt. 0 .or. imult.gt.10) then
        ! error
        ierr=39
        return
     endif

     select case(imult)
     case(0)
        ! even nullapole
        g(:,i) = twopi*rz**2
        
     case(1)
        ! odd nullapole
        g(:,i) = 0.
        
     case(2)
        ! even dipole
        g(:,i) = twopi*(r**2 - rz**2)/2.
        
     case(3)
        ! odd dipole
        co=twopi/rz
        g(:,i) = co*(r**2*z)
        
     case(4)
        ! even quadrapole
        co=pi/(4.*rz**2)
        g(:,i) = co*(r**4-4.*r**2*z**2 - 2.*r**2*rz**2+rz**4)
        
     case(5)
        ! odd quadrapole
        co=pi/(3.*rz**3)
        g(:,i) = co*r**2*z*(3.*r**2-4.*z**2-3.*rz**2)

     case(6)
        ! even hexapole
        co=pi/(12.*rz**4)
        g(:,i) = co*(r**6 - 12.*r**4*z**2 - 3.*r**4*rz**2       &
             + 8.*r**2*z**4 + 12.*r**2*z**2*rz**2           &
             + 3.*r**2*rz**4 - rz**6 )

     case(7)
        ! odd hexapole
        co=pi/(30.*rz**5)
        g(:,i) = co*(15.*r**6*z - 60.*r**4*z**3                 &
             - 30.*r**4*z*rz**2 + 24.*r**2*z**5             &
             + 40.*r**2*z**3*rz**2 + 15.*r**2*z*rz**4)

     case(8)
        ! even octapole
        co=pi/(160.*rz**6)
        g(:,i) = co*(5.*r**8 - 120.*r**6*z**2 - 20.*r**6*rz**2  &
             + 240.*r**4*z**4 + 240.*r**4*z**2*rz**2        &
             + 30.*r**4*rz**4 - 64.*r**2*z**6                  &
             - 160.*r**2*z**4*rz**2 - 120.*r**2*z**2*rz**4  &
             - 20.*r**2*rz**6 + 5.*rz**8)

     case(9)
        ! odd octapole
        co=pi/(140.*rz**7)
        g(:,i) = co*r**2*z*(35.*r**6 - 280.*r**4*z**2        &
             - 105.*r**4*rz**2 + 336.*r**2*z**4                &
             + 420.*r**2*z**2*rz**2 + 105.*r**2*rz**4          &
             - 64.*z**6 - 168.*z**4*rz**2                         &
             - 140.*z**2*rz**4 - 35.*rz**6)

     case(10)
        ! even decapole
        co=pi/(560.*rz**8)
        g(:,i) = co*(7.*r**10 - 280.*r**8*z**2 - 35.*r**8*rz**2 &
             + 1120.*r**6*z**4 + 840.*r**6*z**2*rz**2       &
             + 70.*r**6*rz**4 - 896.*r**4*z**6                 &
             - 1680.*r**4*z**4*rz**2 - 840.*r**4*z**2*rz**4 &
             - 70.*r**4*rz**6 + 128.*r**2*z**8                 &
             + 448.*r**2*z**6*rz**2 + 560.*r**2*z**4*rz**4  &
             + 280.*r**2*z**2*rz**6 + 35.*r**2*rz**8           &
             - 7.*rz**10)

     case default
        print *, 'Error: unknown multipole ', imult
     end select
  end do
end subroutine gvect0


! Int(cos(2 ntor x)/(1 + k2*sin(x)^2)^(3/2) dx)/pi     0 < x < pi/2
subroutine integral(npts,k2,ntor,f)
  use math

  implicit none

  real, intent(in) :: k2(npts)
  real, intent(out) :: f(npts)
  integer, intent(in) :: ntor, npts

  real :: dx, dx2, x
  real, dimension(npts) :: f0, f1, f2
  integer :: ix
  integer, parameter :: ixmax = 100

  dx = (pi/2.) / ixmax
  dx2 = dx/2.

  ! Integrate using Simpson's rule
  f = 0.
  f0 = 1.
  x = 0.
  do ix=1, ixmax
     f1 = cos(2.*ntor*(x+dx2))/sqrt((1.+k2*sin(x+dx2)**2)**3)
     f2 = cos(2.*ntor*(x+dx))/sqrt((1.+k2*sin(x+dx)**2)**3)
     f = f + (f0 + 4.*f1 + f2)
     f0 = f2
     x = x + dx
  end do
  f = f * (dx/6.) / pi
end subroutine integral


subroutine coil(curr, r1, z1, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none

  complex, intent(in) :: curr
  real, intent(in) :: r1, z1
  integer, intent(in) :: ntor, npts
  real, intent(in), dimension(npts) :: r0, z0
  complex, intent(inout), dimension(npts) :: fr, fphi, fz

  real, dimension(npts) :: dr2, k2, temp, co
  complex, dimension(npts) :: fac
  real, parameter :: dr2_fac = 1e-4

  dr2 = (r0 - r1)**2 + (z0 - z1)**2 + dr2_fac
  k2 = 4.*r1*r0/dr2

  fac = curr*r1/sqrt(dr2**3)

  ! calculate contribution from coils
  call integral(npts,k2,ntor,temp)
  fz = fz + fac*r1*temp
  
  call integral(npts,k2,ntor+1,co)
  call integral(npts,k2,ntor-1,temp)

  fr = fr + fac*(z0 - z1)*(co + temp)/2.
  fz = fz - fac*r0*(co + temp)/2.
  fphi = fphi - fac*(0.,1.)*(z0 - z1)*(co - temp)/2.
  
end subroutine coil

subroutine surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none

  real, intent(in) :: r, z, dldZ, dldR
  integer, intent(in) :: npts, ntor
  real, intent(in), dimension(npts) :: r0, z0
  complex, intent(out), dimension(npts) :: fr, fphi, fz
  
  real, dimension(npts) :: dr2, fac, k2, temp1, temp2
  real, parameter :: dr2_fac = 1e-4

  dr2 = (r - r0)**2 + (z - z0)**2 + dr2_fac
  k2 = 4.*r*r0/dr2
  fac = -ntor/sqrt(dr2**3)

  call integral(npts, k2, ntor+1, temp1)
  call integral(npts, k2, ntor-1, temp2)

  fr =  fac*(dldZ*r+dldR*(z0-z))*(temp1 - temp2)/2.
  fz = -fac*dldR*r0*(temp1 - temp2)/2.
  fphi =  -fac*(0.,1.)*(dldZ*r + dldR*(z0 - z))*(temp1 + temp2)/2.

  call integral(npts, k2, ntor, temp1)
  fphi = fphi + fac*(0.,1.)*dldZ*r0*temp1
end subroutine surface

!======================================================================
! pane
! ~~~~
! Adds 'window pane' contribution to magnetic field (fr, fphi, fz)
! at observation point r=(r0, z0) from two axisymmetric coils
! at (r1,z1) and (r2,z2) carrying currents curr and -curr
!======================================================================
subroutine pane(curr, r1, r2, z1, z2, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none

  complex, intent(in) :: curr
  real, intent(in) :: r1, r2, z1, z2
  integer, intent(in) :: ntor, npts
  real, intent(in), dimension(npts) :: r0, z0
  complex, intent(inout), dimension(npts) :: fr, fphi, fz

  complex, dimension(npts) :: fr0, fphi0, fz0
  complex, dimension(npts) :: fr1, fphi1, fz1
  complex, dimension(npts) :: fr2, fphi2, fz2
  real :: l, dl, lmax, dldZ, dldR, r, z, dl2
  integer :: il
  integer, parameter :: ilmax = 100

  ! calculate conitrbution from coils
  call coil( curr,r1,z1,npts,r0,z0,ntor,fr,fphi,fz)
  call coil(-curr,r2,z2,npts,r0,z0,ntor,fr,fphi,fz)

  lmax = sqrt((r2-r1)**2 + (z2-z1)**2)
  dl = lmax/ilmax
  dl2 = dl/2.
  
  dldR = (r2 - r1)/lmax
  dldZ = (z2 - z1)/lmax

  ! calculate contribution from surface currents 
  r = r1
  z = z1
  call surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr0, fphi0, fz0)
  l = 0.
  do il=1, ilmax
     
     r = (r2*(l+dl2) + r1*(lmax-(l+dl2)))/lmax
     z = (z2*(l+dl2) + z1*(lmax-(l+dl2)))/lmax
     call surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr1, fphi1, fz1)
       
     r = (r2*(l+dl) + r1*(lmax-(l+dl)))/lmax
     z = (z2*(l+dl) + z1*(lmax-(l+dl)))/lmax
     call surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr2, fphi2, fz2)

     fr = fr + curr*dl*(fr0 + 4.*fr1 + fr2)/6.
     fphi = fphi + curr*dl*(fphi0 + 4.*fphi1 + fphi2)/6.
     fz = fz + curr*dl*(fz0 + 4.*fz1 + fz2)/6.

     fr0 = fr2
     fphi0 = fphi2
     fz0 = fz2
     
     l = l + dl
  end do
  
end subroutine pane

end module coils

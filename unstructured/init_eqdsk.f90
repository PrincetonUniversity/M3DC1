!==============================================================================
! Eqdsk_eq
! ~~~~~~~~
!
! Loads eqdsk equilibrium
!==============================================================================
module eqdsk_eq

contains

subroutine eqdsk_init()
  use math
  use basic
  use arrays
  use eqdsk
  use gradshafranov
  use newvar_mod
  use sparse
  use diagnostics
  use mesh_mod
  use m3dc1_nint

  implicit none

  integer :: l, ll, ierr, itri, k, numelms
  real :: dpsi, ffp2, pp2
  vectype, parameter ::  negone = -1
  vectype, dimension(dofs_per_element) :: dofs
  type(field_type) :: psi_vec, bz_vec, den_vec, p_vec

  real, allocatable :: flux(:), nflux(:)

  if(myrank.eq.0 .and. iprint.gt.0) print *, "before load_eqdsk", iread_eqdsk
  call load_eqdsk(ierr)
  if(ierr.ne.0) call safestop(1)

  press = press*amu0
  pprime = pprime*amu0
  current = current*amu0

  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, 'normalized current ', current
  end if

  ! create spline containing q profile as function of psi_N
  if(qpsi(nw) .eq. qpsi(nw)) then       ! ensure that qpsi isn't NaN
    allocate(nflux(nw))
    do l=1, nw
       nflux(l) = l / (nw-1.)
    end do
    call create_spline(q_spline,nw,nflux,qpsi)
    deallocate(nflux)
  endif

  if(current.lt.0 .and. iflip_j.eq.1) then
     if(myrank.eq.0) then 
        print *, 'WARNING: iflip_j = 1 for case with negative current.'
        print *, 'Default behavior has been changed to take negative current'
        print *, 'into account automatically.  To override this, set'
        print *, 'iflip_j = 2.'
     end if
     call safestop(1)
  end if

  if(iflip_z.eq.1) zmaxis = -zmaxis

  tcuro = current
  xmag = rmaxis
  zmag = zmaxis
  if(xmag0.eq.0.) then
     xmag0 = rmaxis
     zmag0 = zmaxis
  end if
  rzero = rmaxis

  if(ifixedb.eq.0) then 
     if(iread_eqdsk.eq.3 .and. ifixedb.eq.0) then
        igs_calculate_ip_fields = .true.
     end if
     if(iread_eqdsk.eq.3 .or. idevice.eq.-1 .or. imulti_region.eq.1) then
        igs_calculate_pf_fields = .true.
     end if
     call vacuum_field
  end if

  if(iread_eqdsk.eq.3) then 
     
 ! define initial field associated with delta-function or gaussian source
  !     corresponding to current tcuro at location (xmag,zmag)
     if(myrank.eq.0 .and. iprint.gt.0) write(*,1001) xmag,zmag,tcuro,sigma0
1001 format(' in gradshafranov_init',/,'   xmag,zmag,tcuro,sigma0',1p4e12.4)
     if(sigma0 .eq.0) then
        call deltafun(xmag,zmag,tcuro,jphi_field)
     else
        call gaussianfun(xmag,zmag,tcuro,sigma0,jphi_field)
     endif
   
  else
     
    if(myrank.eq.0 .and. iprint.ge.1) &
         print *, "Interpolating geqdsk Equilibrium"

    call create_field(psi_vec)
    call create_field(bz_vec)
    call create_field(p_vec)
    call create_field(den_vec)

    numelms = local_elements()

    do k=0,1
       psi_vec = 0.
       bz_vec = 0.
       p_vec = 0.
       den_vec = 0.
       
       do itri=1,numelms
          call define_element_quadrature(itri,int_pts_main,int_pts_tor)
          call define_fields(itri,0,1,0)

          if(k.eq.0) then 
             ! calculate equilibrium fields
             call eqdsk_equ
          else
             ! calculate perturbed fields
             call eqdsk_per
          end if

          ! populate vectors for solves

          ! psi
          dofs = intx2(mu79(:,:,OP_1),ps079(:,OP_1))
          call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! bz
          dofs = intx2(mu79(:,:,OP_1),bz079(:,OP_1))
          call vector_insert_block(bz_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! p
          dofs = intx2(mu79(:,:,OP_1),p079(:,OP_1))
          call vector_insert_block(p_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! den
          dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
          call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
       end do

       ! do solves
       call newvar_solve(psi_vec%vec,mass_mat_lhs)
       psi_field(k) = psi_vec
       
       call newvar_solve(bz_vec%vec,mass_mat_lhs)
       bz_field(k) = bz_vec
       
       call newvar_solve(p_vec%vec,mass_mat_lhs)
       p_field(k) = p_vec
       pe_field(k) = p_vec
       call mult(pe_field(k), pefac)
       
       call newvar_solve(den_vec%vec,mass_mat_lhs)
       den_field(k) = den_vec
    end do

    call destroy_field(psi_vec)
    call destroy_field(bz_vec)
    call destroy_field(den_vec)
    call destroy_field(p_vec)

  end if
!
! Bateman scaling parameter reintroduced
  if(igs_pp_ffp_rescale.ne.1) fpol(nw) = fpol(nw)*batemanscale
!
  bzero = fpol(nw)/rzero
  if(iprint.ge.1 .and. myrank.eq.0) then 
     write(*,'(A,2F12.4)') 'Setting bzero, rzero = ', bzero, rzero
  end if

  if(igs.ne.0) then
     if(iread_eqdsk.eq.2) then
        call default_profiles
     else
        allocate(flux(nw))
        dpsi = (sibry-simag)/(nw-1.)
        
        do l=1,nw
           flux(l) = (l-1)*dpsi
           ll = nw - l
           if(batemanscale.eq.1.0 .or. igs_pp_ffp_rescale.eq.1) cycle
! ...Apply Bateman scaling --- redefine fpol keeping ffprim fixed
           if(ll.gt.0) fpol(ll) = sign(1.0,fpol(nw)) &
                *sqrt(fpol(ll+1)**2 - dpsi*(ffprim(ll)+ffprim(ll+1)))
        end do
        call create_profile(nw,press,pprime,fpol,ffprim,flux)
        
        if(qpsi(nw) .eq. qpsi(nw)) then   ! ensure that qpsi isn't NaN
           call create_rho_from_q(nw,flux,qpsi)
        else
           if(myrank.eq.0) print *, "toroidal flux not defined"
           flux=0
        endif
        if(myrank.eq.0 .and. iprint.ge.1) then
           open(unit=77,file="debug-out",status="unknown")
           write(77,2010) sibry,simag,tcuro,xmag,zmag
2010       format("sibry,simag,tcuro,xmag,zmag =",1p5e12.4,/,  &
                "  l   press       pprime      fpol        ffprim      ffp2        pp2         flux        qpsi")
           do l=1,nw
              if(l.gt.1 .and. l.lt.nw)  then
                ffp2 = fpol(l)*(fpol(l+1)-fpol(l-1))/(2*dpsi)
                pp2 = (press(l+1)-press(l-1))/(2*dpsi)
              else
                ffp2 = 0
                pp2 = 0
              endif
              write(77,2011) l,press(l),pprime(l),fpol(l),ffprim(l),ffp2,pp2,flux(l),qpsi(l)
           enddo
2011       format(i3,1p8e12.4)
           close(77)
        endif
        deallocate(flux)
     end if

     psibound = sibry
     psimin = simag

     call gradshafranov_solve
     call gradshafranov_per
  else
     psibound = sibry
     psimin = simag
  endif

  call unload_eqdsk

  ! flip psi sign convention
  ! (iread_eqdsk==3 does not use the eqdsk psi)
  if(iread_eqdsk.eq.1 .or. iread_eqdsk.eq.2) then
     call mult(psi_field(0), negone)
     call mult(psi_field(1), negone)
     if(icsubtract.eq.1) call mult(psi_coil_field, negone)
     psibound = -psibound
     psimin = -psimin
     psilim = -psilim
  end if

  if(iprint.ge.1 .and. myrank.eq.0) &
       write(*,2012) sibry,simag,psimin,psilim,psibound
2012 format(" sibry, simag, psimin, psilim,psibound =",1p5e12.4)

end subroutine eqdsk_init

subroutine eqdsk_equ()
  use basic
  use arrays
  use m3dc1_nint

  use eqdsk

  implicit none

  real :: p, q
  integer :: i, j, n, m, k

  real :: dx, dz, temp, dpsi
  real, dimension(4,4) :: a
  real, dimension(4) :: b, c
  
  dx = rdim/(nw - 1.)
  dz = zdim/(nh - 1.)

  ps079(:,OP_1) = 0.
  p079(:,OP_1) = 0.
  n079(:,OP_1) = 0.
  bz079(:,OP_1) = 0.

  do k=1, npoints
     p = (x_79(k)-rleft)/dx + 1.
     q = (z_79(k)-zmid)/dz + nh/2. + .5
     i = p
     j = q

     call bicubic_interpolation_coeffs(psirz,nw,nh,i,j,a)

     do n=1, 4
        do m=1, 4
           temp = a(n,m)
           if(n.gt.1) temp = temp*(p-i)**(n-1)
           if(m.gt.1) temp = temp*(q-j)**(m-1)
           ps079(k,OP_1) = ps079(k,OP_1) + temp
        end do
     end do

     ! calculation of p
     ! ~~~~~~~~~~~~~~~~
     dpsi = (sibry - simag)/(nw - 1.)
     p = (ps079(k,OP_1) - simag)/dpsi + 1.
     i = p

     if(i.gt.nw) then
        p079(k,OP_1) = press(nw)
        bz079(k,OP_1) = fpol(nw)
     else
        ! use press and fpol to calculate values of p and I
        call cubic_interpolation_coeffs(press,nw,i,b)
        call cubic_interpolation_coeffs(fpol,nw,i,c)

        do n=1,4
           temp = b(n)
           if(n.gt.1) temp = temp*(p-i)**(n-1)
           p079(k,OP_1) = p079(k,OP_1) + temp
           
           temp = c(n)
           if(n.gt.1) temp = temp*(p-i)**(n-1)
           bz079(k,OP_1) = bz079(k,OP_1) + temp
        end do
     endif
  end do

  if(pedge.ge.0.) p079(:,OP_1) = p079(:,OP_1) + pedge

  where(real(p079(:,OP_1)).lt.0.) p079(:,OP_1) = 0.

  ! Set density
  if(expn.eq.0.) then
     n079(:,OP_1) = 1.
  else
     n079(:,OP_1) = (p079(:,OP_1)/p0)**expn
  end if
end subroutine eqdsk_equ

subroutine eqdsk_per
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  p079(:,OP_1) = 0.
  ps079(:,OP_1) = 0.
  bz079(:,OP_1) = 0.
  n079(:,OP_1) = 0.

end subroutine eqdsk_per
  
end module eqdsk_eq

      subroutine ludefden
      use p_data
      use t_data
      use basic
      use arrays
      use sparse_matrix
      use sparse
#ifdef mpi
      use supralu_dist_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer numnodes, nrank1, itri, jfirst, jlast, 
     &     ll, i, j, j1, i1, k, numelms, 
     &     icount, l, jer, ier
      real terma, termb, termd, hypp, xmin, zmin, deex
      real ssterm(18,18),ddterm(18,18),rrterm(18,18),qqterm(18,18), 
     &     dofs(20), fintl(-6:maxi,-6:maxi), k0term_vels13(18,18),
     &     k0term_vel013(18,18),k0term_den12(18,18),
     &     k0term_den012(18,18), aterm(18,18), bterm(18,18), 
     &     dterm(18,18)

!
      call numnod(numnodes)
      call numfac(numelms)
      call getmincoord(xmin,zmin)
      write(*,*) 'in ludefden'
      nrank1 = numnodes*6 
!
!.....form field matrices
        call zero_array(s8matrix_sm,spo_numvar1,58)
        call zero_array(d8matrix_sm,spo_numvar1,68)
        if(numvar .eq. 1) then
           call zero_array(q8matrix_sm,spo_numvar1,981)
           call zero_array(r8matrix_sm,spo_numvar1,181)
        else if(numvar .eq. 2) then
           call zero_array(q8matrix_sm,spo_numvar2,982)
           call zero_array(r8matrix_sm,spo_numvar2,182)
        else if(numvar .eq. 3) then
           call zero_array(q8matrix_sm,spo_numvar3,983)
           call zero_array(r8matrix_sm,spo_numvar3,183)
        endif

#ifdef mpi
! acbauer -- do load-balancing here
!
#endif
        do itri=1,numelms
           call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
           call calcaterm(itri, aterm, fintl)
           call calcbterm(itri, bterm, fintl)
           call calcdterm(itri, dterm, fintl)
           call getlocaldofs(itri, numvar, dofs, 1., vels, isvaln,
     &          0., vels, isvaln)
           call calck0term(itri, k0term_vels13, fintl, dofs,1,3,2)
           call getlocaldofs(itri, numvar, dofs, 1., vel0, isvaln,
     &          0., vel0, isvaln)
           call calck0term(itri, k0term_vel013, fintl, dofs,1,3,2)
           call getlocaldofs(itri, 1, dofs, 1., den, isval1,
     &          0., den, isval1)
           call calck0term(itri, k0term_den12, fintl, dofs,1,2,3)
           call getlocaldofs(itri, 1, dofs, 1., den0, isval1,
     &          0., den0, isval1)
           call calck0term(itri, k0term_den012, fintl, dofs,1,2,3)
           
           call getdeex(itri,deex)
           ssterm = 0.
           ddterm = 0.
           rrterm = 0.
           qqterm = 0.
           do j=1,18
              do i=1,18         
                 terma = aterm(i,j)
                 termb = bterm(i,j)
                 termd = dterm(i,j)
                 hypp = hyperp*deex**2*termb
                                !
!...NOTE: s8matrix_lu (formed from ss) multiplies the solution at time NTIME+1
!
! LINEAR TERMS
            ssterm(i,j) =  ssterm(i,j) + termd - thimp*dt*kappa*terma
     &                           + thimp*dt*kappa*hypp
            ddterm(i,j) = ddterm(i,j)+ termd + (1.-thimp)*dt*kappa*terma
     &                           - (1.-thimp)*dt*kappa*hypp
c     simple nonlinear terms
            ssterm(i,j) = ssterm(i,j) - dt*thimp*
     &           (k0term_vels13(i,j) + k0term_vel013(i,j))
            ddterm(i,j) = ddterm(i,j) + dt* 
     &           ((-thimp+.5)*k0term_vels13(i,j) + 
     &           (-thimp+1.)*k0term_vel013(i,j))
            rrterm(i,j) = rrterm(i,j) + thimp*dt*
     &           (k0term_den12(i,j)+k0term_den012(i,j))
            qqterm(i,j) = qqterm(i,j) + dt*
     &           ((-thimp+.5)*k0term_den12(i,j) + 
     &           (-thimp+1.)*k0term_den012(i,j))
            
!
         enddo
      enddo
      do j=1,18
         do i=1,18
            j1 = isvaln(itri,j)
            i1 = isvaln(itri,i)
         call insert_val(s8matrix_sm, ssterm(i,j), i1, j1, 1)
         call insert_val(d8matrix_sm, ddterm(i,j), i1, j1, 1)
         call insert_val(q8matrix_sm, qqterm(i,j), i1, j1, 1)
         call insert_val(r8matrix_sm, rrterm(i,j), i1, j1, 1)
         if(numvar .gt. 1) then
            call safestop(5598) !this is not right below...
            call insert_val(q8matrix_sm,qqterm(1,2),i1,j1+6,1)
            call insert_val(q8matrix_sm,qqterm(2,1),i1+6,j1,1)
            call insert_val(q8matrix_sm,qqterm(2,2),i1+6,j1+6,1)
            call insert_val(r8matrix_sm,rrterm(1,2),i1,j1+6,1)
            call insert_val(r8matrix_sm,rrterm(2,1),i1+6,j1,1)
            call insert_val(r8matrix_sm,rrterm(2,2),i1+6,j1+6,1)
         endif
         if(numvar .eq. 3) then
            call insert_val(q8matrix_sm,qqterm(1,3),i1,j1+12,1)
            call insert_val(q8matrix_sm,qqterm(2,3),i1+6,j+12,1)
            call insert_val(q8matrix_sm,qqterm(3,3),i1+12,j1+12,1)
            call insert_val(q8matrix_sm,qqterm(3,1),i1+12,j1,1)
            call insert_val(q8matrix_sm,qqterm(3,2),i1+12,j1+6,1)
            call insert_val(r8matrix_sm,rrterm(1,3),i1,j1+12,1)
            call insert_val(r8matrix_sm,rrterm(2,3),i1+6,j+12,1)
            call insert_val(r8matrix_sm,rrterm(3,3),i1+12,j1+12,1)
            call insert_val(r8matrix_sm,rrterm(3,1),i1+12,j1,1)
            call insert_val(r8matrix_sm,rrterm(3,2),i1+12,j1+6,1)
         endif
      enddo                     ! on i
      enddo                     ! on j
      enddo                     ! on itri
!
!.....modify the s-matrix, inserting the boundary conditions
!
!.....calculate boundary conditions
      call boundaryds(iboundp,nbcp,1)
      if(nbcp .gt. iboundmax) then
      write(*,4888) nbcp, iboundmax
 4888 format(" ERROR: nbcp > iboundmax", 2i5)
      call safestop(9) 
      endif

      !  bc stuff...
      do i=1,nbcp
         call set_diri_bc(s8matrix_sm, iboundp(i))
      enddo
      call finalize_array(s8matrix_sm)
      call finalize_array(d8matrix_sm)
      call finalize_array(q8matrix_sm)
      call finalize_array(r8matrix_sm)

      return
      end
!=====================================
      subroutine ludefvel
      use p_data
      use t_data
      use basic
      use arrays
      use sparse_matrix
      use sparse
#ifdef mpi
      use supralu_dist_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer itri, l, i, j, j1, i1, numelms, numnodes, j01
      real hypv, tstartludefvel,
     &     hypvf, terma, termb, termd,
     &     tboundludefvel, tluludefvel, tendludefvel, tbc, tmat,tfactor,
     &     tmidludefvel, deex
      real ssterm(18,18),ddterm(18,18),rrterm(18,18),
     &     dofs(20), g2term(18,18), fintl(-6:maxi,-6:maxi), 
     &     c1term2(18,18), v1terml(18,18), v2terml(18,18),
     &     g0term_vel12(18,18), g0term_vel13(18,18),
     &     g0term_vel012(18,18), g0term_vel013(18,18),
     &     g0term_vels12(18,18), g0term_vels13(18,18),
     &     g0term_phi12(18,18), g0term_phi13(18,18),
     &     g0term_phi012(18,18), g0term_phi013(18,18),
     &     g0term_sphi112(18,18), g0term_sphi113(18,18), 
     &     dofs1(20), dofs2(20), aterm(18,18), bterm(18,18),
     &     dterm(18,18)

      call numfac(numelms)
      call numnod(numnodes)
      if(myrank.eq.0 .and. itimer.gt.0) call second(tstartludefvel)
!
!.....form velocity matrices
      if(numvar .eq. 1) then
         call zero_array(s1matrix_sm,spo_numvar1, 511)
         call zero_array(d1matrix_sm,spo_numvar1, 611)
         call zero_array(r1matrix_sm,spo_numvar1, 111)
      else if(numvar .eq. 2) then
         call zero_array(s1matrix_sm,spo_numvar2, 512)
         call zero_array(d1matrix_sm,spo_numvar2, 612)
         call zero_array(r1matrix_sm,spo_numvar2, 112)
      else 
         call zero_array(s1matrix_sm,spo_numvar3, 513)
         call zero_array(d1matrix_sm,spo_numvar3, 613)
         call zero_array(r1matrix_sm,spo_numvar3, 113)
      endif
      r4 = 0.
      vels = vel
      
! acbauer -- load-balance here
!
!
      do itri=1,numelms
         call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
         call calcaterm(itri, aterm, fintl)
         call calcbterm(itri, bterm, fintl)
         call calcdterm(itri, dterm, fintl)
        call getlocaldofs(itri, numvar, dofs1, 1., phi, isvaln,
     &        1., phi0, isvaln) ! second isvaln should use isval0
         call getlocaldofs(itri, numvar, dofs2, 1., phi, isvaln,
     &        1., phi0, isvaln) ! second isvaln should use isval0
         call calcc1term2(itri, c1term2, fintl, dofs2, dofs1, 1,2,3,4)
         
         call getlocaldofs(itri, numvar, dofs1, 1., vel, isvaln, 
     &        0., vel, isvaln)
         call calcg0term(itri, g0term_vel12, fintl, dofs1, 1, 2, 3)
         call calcg0term(itri, g0term_vel13, fintl, dofs1, 1, 3, 2)
         call getlocaldofs(itri, numvar, dofs1, 1., vel0, isvaln,
     &        0., vel0, isvaln)
         call calcg0term(itri, g0term_vel012, fintl, dofs1, 1, 2, 3)
         call calcg0term(itri, g0term_vel013, fintl, dofs1, 1, 3, 2)
         call getlocaldofs(itri, numvar, dofs1, 1., vels, isvaln,
     &        0., vels, isvaln)
         call calcg0term(itri, g0term_vels12, fintl, dofs1, 1, 2, 3)
         call calcg0term(itri, g0term_vels13, fintl, dofs1, 1, 3, 2)
         call getlocaldofs(itri, numvar, dofs1, 1., phi, isvaln,
     &        0., phi, isvaln)
         call calcg0term(itri, g0term_phi12, fintl, dofs1, 1, 2, 3)
         call calcg0term(itri, g0term_phi13, fintl, dofs1, 1, 3, 2)
         call getlocaldofs(itri, numvar, dofs1, 1., phi0, isvaln,
     &        0., phi0, isvaln) ! both isvaln's should use isval0
         call calcg0term(itri, g0term_phi012, fintl, dofs1, 1, 2, 3)
         call calcg0term(itri, g0term_phi013, fintl, dofs1, 1, 3, 2)
         call getlocaldofs(itri, 1, dofs1, 1., sb1, isval1,
     &        0., sb1, isval1)
         call calcg0term(itri, g0term_sphi112, fintl, dofs1, 1, 2, 3)
         call calcg0term(itri, g0term_sphi113, fintl, dofs1, 1, 3, 2)

         if(idens.gt.0) then
            call getlocaldofs(itri, 1, dofs1, 1., den+den0, isval1, 
     &           0., den+den0, isval1)
            call calcg2term(itri, g2term, fintl, dofs1, 1, 2, 3)
            call getlocaldofs(itri, numvar, dofs2, 1., vel, isvaln,
     &           1., vel0, isvaln)
            call calcv1terml(itri, v1terml, fintl, dofs2, dofs1, 
     &           1, 2, 3, 4)
            call calcv2terml(itri, v2terml, fintl, dofs2, dofs1, 
     &           1, 2, 3, 4)
         endif
         
         call getdeex(itri,deex)
         hypv = hyperv*deex**2
         hypvf = hypv
!
         ssterm = 0
         ddterm = 0
         rrterm = 0

         do i=1,18
            do j=1,18
                                
               terma = aterm(i,j)
               termb = bterm(i,j)
               termd = dterm(i,j)
               i1 = isvaln(itri,i)
               j01 = isvaln(itri,j) ! should use isval0
!...NOTE: s1matrix_lu (formed from ss) multiplies the solution at time NTIME+1
!
! LINEAR TERMS
               ssterm(i,j) = terma*(1-idens)      - thimp*dt*amu*termb
               ddterm(i,j) = terma*(1-idens) + (1.-thimp)*dt*amu*termb
!
!
!
!
! SIMPLE NONLINEAR TERMS
               if(idens.gt.0) then
                  ssterm(i,j) = ssterm(i,j)+g2term(i,j)
                  ddterm(i,j) = ddterm(i,j)+g2term(i,j)
                  ssterm(i,j) = ssterm(i,j) + dt*thimp*
     &                 (v1terml(i,j) + v2terml(i,j))
                  ddterm(i,j) = ddterm(i,j) + dt*(thimp-0.5)*
     &                 (v1terml(i,j) + v2terml(i,j))
               endif
               ssterm(i,j) = ssterm(i,j) + dt*(1-idens)*thimp*
     &              (g0term_vel12(i,j) + g0term_vel13(i,j) + 
     &              g0term_vel012(i,j) + g0term_vel013(i,j))
               
               ddterm(i,j) = ddterm(i,j) + dt*(1-idens)*
     &              (thimp*(g0term_vels12(i,j)+g0term_vels13(i,j)+
     &              g0term_vel012(i,j) + g0term_vel013(i,j)) - 
     &              (.5*g0term_vel12(i,j)+.5*g0term_vel13(i,j) +
     &              g0term_vel012(i,j)+g0term_vel013(i,j)))
               
               rrterm(i,j) = rrterm(i,j) + 
     &              dt*((.5*g0term_phi12(i,j)+.5*g0term_phi13(i,j))+
     &              g0term_phi012(i,j)+g0term_phi013(i,j)) +
     &             thimp*dt*dt*(g0term_sphi112(i,j)+g0term_sphi113(i,j))
               
               r4(i1) = r4(i1) + thimp*dt*dt*(g0term_sphi112(i,j)+
     &              g0term_sphi113(i,j))*phi0(j01)
               
               ssterm(i,j) = ssterm(i,j) + 
     &              thimp*thimp*dt*dt*c1term2(i,j)   
               ddterm(i,j) = ddterm(i,j) + 
     &              thimp*(thimp-1.)*dt*dt*c1term2(i,j)
            enddo               ! on j
         enddo                  ! on i
         
         do i=1,18
            i1 = isvaln(itri,i)
            do j=1,18
               j1 = isvaln(itri,j)
               call insert_val(s1matrix_sm,ssterm(i,j),i1,j1,1)
               call insert_val(d1matrix_sm,ddterm(i,j),i1,j1,1)
               call insert_val(r1matrix_sm,rrterm(i,j),i1,j1,1)
            enddo               ! on j
         enddo                  ! on i
         
      enddo                     ! on itri
!
      if(myrank.eq.0 .and. itimer.gt.0) call second(tmidludefvel)
!
!.....modify the s-matrix, inserting the boundary conditions
!
!.....calculate boundary conditions
      velbounds = 0.
      call boundaryv(iboundv,nbcv)
      if(nbcv .gt. iboundmax) then
      write(*,4888) nbcv, iboundmax
 4888 format(" ERROR: nbcv > iboundmax", 2i5)
      call safestop(3)
      endif
      do l=1,nbcv
         call set_diri_bc(s1matrix_sm, iboundv(l))
      enddo
      call finalize_array(s1matrix_sm)
      call finalize_array(d1matrix_sm)
      call finalize_array(r1matrix_sm)
      if(myrank.eq.0 .and. itimer.gt.0) call second(tboundludefvel)

!
      if(myrank.eq.0 .and. itimer.gt.0) call second(tluludefvel)
      if(myrank.eq.0 .and. itimer.gt.0) then
          call second(tendludefvel)
          tsetup = tmidludefvel - tstartludefvel
          tbc    = tboundludefvel - tmidludefvel
          tmat   = tluludefvel - tboundludefvel
          tfactor= tendludefvel - tluludefvel
          write(*,3001) tsetup,tbc,tmat,tfactor
 3001     format("LUDEFVEL:   setup and factor times:",1p4e12.4)
          endif
      return
      end
!============================================================
      subroutine ludefphi
      use p_data
      use t_data
      use basic
      use arrays
      use sparse_matrix
      use sparse
#ifdef mpi
      use supralu_dist_mod
      implicit none
      include 'mpif.h'
#else
      implicit none
#endif
      integer itri, numelms, i, j, l, i1, j1
      real xmin, zmin, fintl(-6:maxi,-6:maxi), dofs(20),g12erm(18,18),
     &     k0term_vels13(18,18), k0term_vel013(18,18), 
     &     k0term_phi12(18,18), k0term_phi012(18,18), factor, x, z,
     &     ssterm(18,18), ddterm(18,18), rrterm(18,18), qqterm(18,18),
     &     terma, termb, termd, hypf, deex, aterm(18,18), bterm(18,18),
     &     dterm(18,18)
      double precision cogcoords(3)

      call getmincoord(xmin,zmin)
      call numfac(numelms)
!
!       implicit parameter for hyperviscosity terms
!
!.....form field matrices
      if(numvar .eq. 1) then
         call zero_array(s2matrix_sm,spo_numvar1,521)
         call zero_array(d2matrix_sm,spo_numvar1,621)
         call zero_array(r2matrix_sm,spo_numvar1,121)
         call zero_array(q2matrix_sm,spo_numvar1,921)
      else if(numvar .eq. 2) then
         call zero_array(s2matrix_sm,spo_numvar2,522)
         call zero_array(d2matrix_sm,spo_numvar2,622)
         call zero_array(r2matrix_sm,spo_numvar2,122)
         call zero_array(q2matrix_sm,spo_numvar2,922)
      else 
         call zero_array(s2matrix_sm,spo_numvar3,523)
         call zero_array(d2matrix_sm,spo_numvar3,623)
         call zero_array(r2matrix_sm,spo_numvar3,123)
         call zero_array(q2matrix_sm,spo_numvar3,923)
      endif
      phis = phi
      
! acbauer -- load-balance here
!.......position files with metric terms
!
        do itri=1,numelms
           call calcfint(fintl, maxi, atri(itri),btri(itri),ctri(itri))
           call calcaterm(itri, aterm, fintl)
           call calcbterm(itri, bterm, fintl)
           call calcdterm(itri, dterm, fintl)
           call getlocaldofs(itri, numvar, dofs, 1.,
     &          deni, isval1, 0., deni, isval1)
           call calcg12erm(itri, g12erm, fintl, dofs, 1, 2, 3)
           
           call getlocaldofs(itri, numvar, dofs, 1., vels, isvaln, 
     &          0., vels, isvaln)
           call calck0term(itri, k0term_vels13, fintl, dofs, 1, 3, 2)
           call getlocaldofs(itri, numvar, dofs, 1., vel0, isvaln, 
     &          0., vel0, isvaln)
           call calck0term(itri, k0term_vel013, fintl, dofs, 1, 3, 2)
           call getlocaldofs(itri, numvar, dofs, 1., phi, isvaln,
     &          0., phi, isvaln)
           call calck0term(itri, k0term_phi12, fintl, dofs, 1, 2, 3)
           call getlocaldofs(itri, numvar, dofs, 1., phi0, isvaln, 
     &          0., phi0, isvaln) ! need to replace both isvaln's with isval0
           call calck0term(itri, k0term_phi012, fintl, dofs, 1, 2, 3)

           call getdeex(itri,deex)

        factor = 1.
        call cogfac(itri,cogcoords)
        x = cogcoords(1)-xmin
        z = cogcoords(2)-zmin
        if(imask.eq.1) call mask(x,z,factor)
        ssterm = 0.
        ddterm = 0.
        rrterm = 0.
        qqterm = 0.
          do i=1,18
          do j=1,18
            terma = aterm(i,j)
            termb = bterm(i,j)
            termd = dterm(i,j)
            hypf = hyper*deex**2


!
!...NOTE: s2matrix_lu (formed from ss) multiplies the solution at time NTIME+1
!
! LINEAR TERMS
            ssterm(i,j) = ssterm(i,j) + termd - thimp  *dt*etar*terma              &
     &                  + (1-idens)*thimp      *dt*etar*hypf*termb
            ddterm(i,j) = ddterm(i,j)+ termd + (1.-thimp) *dt*etar*terma              &
     &                  - (1-idens)*(1.-thimp )*dt*etar*hypf*termb

c SIMPLE NONLINEAR TERMS
            if(idens.ge.1) then
               ssterm(i,j) = ssterm(i,j) + g12erm(i,j)*
     &              thimp      *dt*etar*hypf
               ddterm(i,j) = ddterm(i,j) - g12erm(i,j)*  
     &              (1.-thimp )*dt*etar*hypf
            endif               !  on idens.ge.1
            
            ssterm(i,j) = ssterm(i,j) - dt*thimp*
     &           (k0term_vels13(i,j)+k0term_vel013(i,j))
            ddterm(i,j) = ddterm(i,j) + dt *
     &           ((-thimp+.5)*k0term_vels13(i,j) +
     &           (-thimp+1.)*k0term_vel013(i,j))
            rrterm(i,j) = rrterm(i,j) + thimp*dt *
     &           (k0term_phi12(i,j)+k0term_phi012(i,j))
            qqterm(i,j) = qqterm(i,j) + dt *
     &           ((-thimp+.5)*k0term_phi12(i,j) + 
     &           (-thimp+1.)*k0term_phi012(i,j))
            
         enddo
      enddo
      do j=1,18
         do i=1,18
            j1 = isvaln(itri,j)
            i1 = isvaln(itri,i)
            call insert_val(s2matrix_sm,ssterm(i,j),i1,j1,1)
            call insert_val(d2matrix_sm,ddterm(i,j),i1,j1,1)
            call insert_val(r2matrix_sm,rrterm(i,j),i1,j1,1)
            call insert_val(q2matrix_sm,qqterm(i,j),i1,j1,1)
            if(numvar .gt. 1) then
               call safestop(4429) ! this won't work yet
               call insert_val(s2matrix_sm,ssterm(1,2),i1,j1+6,1)
               call insert_val(s2matrix_sm,ssterm(2,1),i1+6,j1,1)
               call insert_val(s2matrix_sm,ssterm(2,2),i1+6,j1+6,1)
               call insert_val(d2matrix_sm,ddterm(1,2),i1,j1+6,1)
               call insert_val(d2matrix_sm,ddterm(2,1),i1+6,j1,1)
               call insert_val(d2matrix_sm,ddterm(2,2),i1+6,j1+6,1)
               call insert_val(r2matrix_sm,rrterm(1,2),i1,j1+6,1)
               call insert_val(r2matrix_sm,rrterm(2,1),i1+6,j1,1)
               call insert_val(r2matrix_sm,rrterm(2,2),i1+6,j1+6,1)
               call insert_val(q2matrix_sm,qqterm(1,2),i1,j1+6,1)
               call insert_val(q2matrix_sm,qqterm(2,1),i1+6,j1,1)
               call insert_val(q2matrix_sm,qqterm(2,2),i1+6,j1+6,1)
            endif
            if(numvar .eq. 3) then
               call insert_val(s2matrix_sm,ssterm(1,3),i1,j1+12,1)
               call insert_val(s2matrix_sm,ssterm(2,3),i1+6,j1+12,1)
               call insert_val(s2matrix_sm,ssterm(3,3),i1+12,j1+12,1)
               call insert_val(s2matrix_sm,ssterm(3,1),i1+12,j1,1)
               call insert_val(s2matrix_sm,ssterm(3,2),i1+12,j1+6,1)
               call insert_val(d2matrix_sm,ddterm(1,3),i1,j1+12,1)
               call insert_val(d2matrix_sm,ddterm(2,3),i1+6,j1+12,1)
               call insert_val(d2matrix_sm,ddterm(3,3),i1+12,j1+12,1)
               call insert_val(d2matrix_sm,ddterm(3,1),i1+12,j1,1)
               call insert_val(d2matrix_sm,ddterm(3,2),i1+12,j1+6,1)
               call insert_val(r2matrix_sm,rrterm(1,3),i1,j1+12,1)
               call insert_val(r2matrix_sm,rrterm(2,3),i1+6,j1+12,1)
               call insert_val(r2matrix_sm,rrterm(3,3),i1+12,j1+12,1)
               call insert_val(r2matrix_sm,rrterm(3,1),i1+12,j1,1)
               call insert_val(r2matrix_sm,rrterm(3,2),i1+12,j1+6,1)
               call insert_val(q2matrix_sm,qqterm(1,3),i1,j1+12,1)
               call insert_val(q2matrix_sm,qqterm(2,3),i1+6,j1+12,1)
               call insert_val(q2matrix_sm,qqterm(3,3),i1+12,j1+12,1)
               call insert_val(q2matrix_sm,qqterm(3,1),i1+12,j1,1)
               call insert_val(q2matrix_sm,qqterm(3,2),i1+12,j1+6,1)
            endif
          enddo
          enddo
        enddo
!
!.....modify the s-matrix, inserting the boundary conditions
!
!.....calculate boundary conditions
      call boundaryp(iboundp,nbcp)
      if(nbcp .gt. iboundmax) then
      write(*,4888) nbcp, iboundmax
 4888 format(" ERROR: nbcp > iboundmax", 2i5)
      call safestop(9) 
      endif
      do l=1,nbcp
         call set_diri_bc(s2matrix_sm,iboundp(l))
      enddo
      call finalize_array(s2matrix_sm)
      call finalize_array(d2matrix_sm)
      call finalize_array(r2matrix_sm)
      call finalize_array(q2matrix_sm)

      return
      end
!================================

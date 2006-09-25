!=====================================
subroutine ludefvel_t
  use p_data
  use t_data
  use basic
  use arrays
  use sparse_matrix
  use sparse
  use nintegrate_mod

#ifdef mpi
  use supralu_dist_mod
  implicit none
  include 'mpif.h'
#else
  implicit none
#endif
  integer itri, l, i, j, j1, i1, i2, i3, numelms
  real :: ssterm(3,3),ddterm(3,3),rrterm(3,3)
  real :: x, xmin, z, zmin, factor, temp

  double precision :: cogcoords(3)

  real :: tstart, tend, tfield, telm

  real :: v3up, v3chip, v3psipsi, v3upsipsi, v3chipsipsi, v3chibb, v3chichin, &
       v1chin, v1uchin, v1chichin, v2vchin, v2chipsib, v1chipsipsi, v1chibb, &
       v1ubb, v3ubb, v3uun, v3uchin, v3vpsib

  tfield = 0.
  telm = 0.

  call numfac(numelms)
  call getmincoord(xmin,zmin)
  
  ! form velocity matrices
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

  veln = vel
  vels = vel
  r4 = 0.
      
! acbauer -- load-balance here
!
  do itri=1,numelms

     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! calculate the field values and derivatives at the sampling points
     call define_fields_25(itri)
     call define_fields_79(itri)
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        tfield = tfield + tend - tstart
     endif
     
!     call getdeex(itri,deex)   

     if(imask.eq.1) then
        call cogfac(itri,cogcoords)
        x = cogcoords(1)-xmin
        z = cogcoords(2)-zmin
        call mask(x,z,factor)
     else
        factor = 1.
     endif
     dbf = db*factor

     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     do i=1,18

        i1 = isvaln(itri,i)
        i2 = i1 + 6
        i3 = i2 + 6

        do j=1,18
           j1 = isvaln(itri,j)

           ssterm = 0.
           ddterm = 0.
           rrterm = 0.

           ! NUMVAR = 1
           ! ~~~~~~~~~~
           call v1un(itri,i,j,ssterm, ddterm, rrterm)        ! passed: 1,2
           call v1umu(itri,i,j,ssterm, ddterm, rrterm)       ! passed: 1
           call v1uun(itri,i,j,ssterm, ddterm, rrterm)       ! passed: 1
           call v1psipsi(itri,i,j,ssterm, ddterm, rrterm)    ! passed: 1, FAILED: 2
           call v1upsipsi(itri,i,j,ssterm, ddterm, rrterm)   ! passed: 1
           call v1psisb1(itri,i,j,ssterm, ddterm, rrterm)    ! passed: 1,2

           ! NUMVAR = 2
           ! ~~~~~~~~~~
           if(numvar.ge.2) then

              temp = v1ubb(g79(:,:,i),g79(:,:,j),bzt79,bzt79)
              ssterm(1,1) = ssterm(1,1) - thimp*    thimp *dt*dt*temp
              ddterm(1,1) = ddterm(1,1) + thimp*(1.-thimp)*dt*dt*temp

              if(linear.eq.1 .or. eqsubtract.eq.1) then
                 rrterm(1,2) = rrterm(1,2) + thimp*dt*dt* &
                      (v1ubb(g79(:,:,i),ph079,g79(:,:,j),bzs79) &
                      +v1ubb(g79(:,:,i),ph079,bzs79,g79(:,:,j)))

                 ! EQUILIBRIUM TERM
                 r4(i2) = r4(i2) + thimp*dt*dt* &
                      (v1ubb(g79(:,:,i),ph079,bz079,bz079))
              endif

              if(itor.eq.1) then
                 call v1vpsib(itri,i,j,ssterm, ddterm, rrterm)    ! passed: 1
                 call v1bb(itri,i,j,ssterm, ddterm, rrterm)       ! passed: 1
                 call v1bsb2(itri,i,j,ssterm, ddterm, rrterm)     ! passed: 1
              endif
              call v2v(itri,i,j,ssterm, ddterm, rrterm)        ! passed: 1
              call v2vmu(itri,i,j,ssterm, ddterm, rrterm)      ! passed: 1
              call v2psib(itri,i,j,ssterm, ddterm, rrterm)     ! passed: 1
              call v2psisb1(itri,i,j,ssterm, ddterm, rrterm)   ! passed: 1
              call v2psisb2(itri,i,j,ssterm, ddterm, rrterm)   ! passed: 1
              call v2vun(itri,i,j,ssterm, ddterm, rrterm)      ! passed: 1
              call v2vpsipsi(itri,i,j,ssterm, ddterm, rrterm)  ! passed: 1
              call v2upsib(itri,i,j,ssterm, ddterm, rrterm)    ! passed: 1
           endif

           ! NUMVAR = 3
           ! ~~~~~~~~~~
           if(numvar.ge.3) then
              if(itor.eq.1) then
                 call v3umu(itri,i,j,ssterm, ddterm, rrterm)
              endif
              call v3chin(itri,i,j,ssterm, ddterm, rrterm)     ! passed: 1
              call v3chimu(itri,i,j,ssterm, ddterm, rrterm)    ! passed: 1
              call v3p(itri,i,j,ssterm, ddterm, rrterm)        ! passed: 1

              call v3bb(itri,i,j,ssterm, ddterm, rrterm)       ! passed: 1
              call v3psisb1(itri,i,j,ssterm, ddterm, rrterm)   ! passed: 1
              call v3bsb2(itri,i,j,ssterm, ddterm, rrterm)     ! passed: 1

              temp = v1uchin(g79(:,:,i),g79(:,:,j),ch179,nt79)        ! passed: 1
              ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
              ddterm(1,1) = ddterm(1,1) + (.5-thimp)*dt*temp

              temp = v1chin(g79(:,:,i),g79(:,:,j),nt79)               ! passed: 1
              ssterm(1,3) = ssterm(1,3) - temp
              ddterm(1,3) = ddterm(1,3) - temp

              temp = v1uchin  (g79(:,:,i),ph179,g79(:,:,j),nt79) &    ! passed: 1
                   + v1chichin(g79(:,:,i),g79(:,:,j),ch179,nt79) &    ! passed: 1
                   + v1chichin(g79(:,:,i),ch179,g79(:,:,j),nt79)      ! passed: 1
              ssterm(1,3) = ssterm(1,3) -     thimp *dt*temp
              ddterm(1,3) = ddterm(1,3) + (.5-thimp)*dt*temp

              temp = 0*v1chipsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79)  & ! FAILED: 1
                   + v1chibb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)
              ssterm(1,3) = ssterm(1,3) - thimp*    thimp *dt*dt*temp
              ddterm(1,3) = ddterm(1,3) + thimp*(1.-thimp)*dt*dt*temp

              temp = v2vchin(g79(:,:,i),g79(:,:,j),ch179,nt79)        ! passed: 1
              ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) + (.5-thimp)*dt*temp

              temp = v2vchin(g79(:,:,i),vz179,g79(:,:,j),nt79)        ! passed: 1
              ssterm(2,3) = ssterm(2,3) -     thimp *dt*temp
              ddterm(2,3) = ddterm(2,3) + (.5-thimp)*dt*temp

              temp = v2chipsib(g79(:,:,i),g79(:,:,j),pst79,bzt79)     ! passed: 1
              ssterm(2,3) = ssterm(2,3) - thimp*    thimp *dt*dt*temp
              ddterm(2,3) = ddterm(2,3) + thimp*(1.-thimp)*dt*dt*temp

              temp = v3uun  (g79(:,:,i),g79(:,:,j),ph179,nt79) &
                   + v3uun  (g79(:,:,i),ph179,g79(:,:,j),nt79) &
                   + 0.*v3uchin(g79(:,:,i),g79(:,:,j),ch179,nt79)             ! FAILED: 1
              ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
              ddterm(3,1) = ddterm(3,1) + (.5-thimp)*dt*temp

              temp = v3up     (g79(:,:,i),g79(:,:,j),pt79) &               ! passed: 1
!                   + v3upsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79) &       ! FAILED: 1
                   + v3ubb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)          ! passed: 1
              ssterm(3,1) = ssterm(3,1) - thimp*    thimp *dt*dt*temp
              ddterm(3,1) = ddterm(3,1) + thimp*(1.-thimp)*dt*dt*temp

              temp = v3vpsib(g79(:,:,i),g79(:,:,j),pst79,bzt79)             ! passed: 1
              ssterm(3,2) = ssterm(3,2) - thimp*    thimp *dt*dt*temp
              ddterm(3,2) = ddterm(3,2) + thimp*(1.-thimp)*dt*dt*temp

              temp = 0.*v3uchin  (g79(:,:,i),ph179,g79(:,:,j),nt79) &    ! FAILED: 1
                   + v3chichin(g79(:,:,i),g79(:,:,j),ch179,nt79) &    ! passed: 1
                   + v3chichin(g79(:,:,i),ch179,g79(:,:,j),nt79)      ! passed: 1
              ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
              ddterm(3,3) = ddterm(3,3) + (.5-thimp)*dt*temp
              
              temp = v3chip     (g79(:,:,i),g79(:,:,j),pt79) &         ! passed: 1
!                   + v3chipsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79) &  ! FAILED: 1
                   + v3chibb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)    ! passed: 1
              ssterm(3,3) = ssterm(3,3) - thimp*    thimp *dt*dt*temp
              ddterm(3,3) = ddterm(3,3) + thimp*(1.-thimp)*dt*dt*temp

              rrterm(3,1) = rrterm(3,1) + &
                   dt* &
                   (v3psipsi(g79(:,:,i),g79(:,:,j),pss79) &       ! passed: 1
                   +v3psipsi(g79(:,:,i),pss79,g79(:,:,j)))        ! passed: 1

              if(linear.eq.1 .or. eqsubtract.eq.1) then
                 temp = v1uchin(g79(:,:,i),g79(:,:,j),ch079,nt79)
                 ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
                 ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

                 temp = v1uchin(g79(:,:,i),ph079,g79(:,:,j),nt79)
                 ssterm(1,3) = ssterm(1,3) -     thimp *dt*temp
                 ddterm(1,3) = ddterm(1,3) + (1.-thimp)*dt*temp

                 temp = v2vchin(g79(:,:,i),g79(:,:,j),ch079,nt79)
                 ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
                 ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp
                 
                 temp = v2vchin(g79(:,:,i),vz079,g79(:,:,j),nt79)
                 ssterm(2,3) = ssterm(2,3) -     thimp *dt*temp
                 ddterm(2,3) = ddterm(2,3) + (1.-thimp)*dt*temp

                 temp = v3uun  (g79(:,:,i),g79(:,:,j),ph079,nt79) &
                      + v3uun  (g79(:,:,i),ph079,g79(:,:,j),nt79) &
                      + v3uchin(g79(:,:,i),g79(:,:,j),ch079,nt79)
                 ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
                 ddterm(3,1) = ddterm(3,1) + (1.-thimp)*dt*temp

                 temp = v3uchin  (g79(:,:,i),ph079,g79(:,:,j),nt79) &
                      + v3chichin(g79(:,:,i),g79(:,:,j),ch079,nt79) &
                      + v3chichin(g79(:,:,i),ch079,g79(:,:,j),nt79)
                 ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
                 ddterm(3,3) = ddterm(3,3) + (1.-thimp)*dt*temp

                 rrterm(1,1) = rrterm(1,1) + thimp*dt*dt* &
                      (v1chipsipsi(g79(:,:,i),ch079,g79(:,:,j),pss79) &
                      +v1chipsipsi(g79(:,:,i),ch079,pss79,g79(:,:,j)))

                 rrterm(1,2) = rrterm(1,2) + thimp*dt*dt* &
                      (v1chibb(g79(:,:,i),ch079,g79(:,:,j),bzs79) &
                      +v1chibb(g79(:,:,i),ch079,bzs79,g79(:,:,j)))

                 rrterm(2,1) = rrterm(2,1) + thimp*dt*dt* &
                      (v2chipsib(g79(:,:,i),ch079,g79(:,:,j),bzs79))

                 rrterm(2,2) = rrterm(2,2) + thimp*dt*dt* &
                      (v2chipsib(g79(:,:,i),ch079,pss79,g79(:,:,j)))

                 rrterm(3,1) = rrterm(3,1) + thimp*dt*dt* &
                      (v3upsipsi  (g79(:,:,i),ph079,g79(:,:,j),pss79) &
                      +v3upsipsi  (g79(:,:,i),ph079,pss79,g79(:,:,j)) &
                      +v3vpsib    (g79(:,:,i),vz079,g79(:,:,j),bzs79) &
                      +v3chipsipsi(g79(:,:,i),ch079,g79(:,:,j),pss79) &
                      +v3chipsipsi(g79(:,:,i),ch079,pss79,g79(:,:,j)))

                 rrterm(3,2) = rrterm(3,2) + thimp*dt*dt* &
                      (v3ubb  (g79(:,:,i),ph079,g79(:,:,j),bzs79) &
                      +v3ubb  (g79(:,:,i),ph079,bzs79,g79(:,:,j)) &
                      +v3vpsib(g79(:,:,i),vz079,pss79,g79(:,:,j)) &
                      +v3chibb(g79(:,:,i),ch079,g79(:,:,j),bzs79) &
                      +v3chibb(g79(:,:,i),ch079,bzs79,g79(:,:,j)))

                 rrterm(3,3) = rrterm(3,3) + thimp*dt*dt* &
                      (v3up  (g79(:,:,i),ph079,g79(:,:,j)) &
                      +v3chip(g79(:,:,i),ch079,g79(:,:,j)))

                 ! DENSITY TERMS
                 r4(i1) = r4(i1) + dt* &
                      (v1uchin  (g79(:,:,i),ph079,ch079,nt79) &
                      +v1chichin(g79(:,:,i),ch079,ch079,nt79))
                 r4(i3) = r4(i3) + dt* &
                      (v3uun    (g79(:,:,i),ph079,ph079,nt79) &
                      +v3uchin  (g79(:,:,i),ph079,ch079,nt79) &
                      +v3chichin(g79(:,:,i),ch079,ch079,nt79))

                 ! EQUILIBRIUM TERMS
                 r4(i2) = r4(i2) + &
                      thimp*dt*dt* &
                      (v2chipsib  (g79(:,:,i),ch079,ps079,bz079))
                 r4(i3) = r4(i3) + &
                      dt* &
                      (v3psipsi(g79(:,:,i),ps079,ps079)) + &
                      thimp*dt*dt* &
                      (v3up       (g79(:,:,i),ph079,p079) &
                      +v3vpsib    (g79(:,:,i),vz079,ps079,bz079) &
                      +v3chip     (g79(:,:,i),ch079,p079) &
                      +v3chipsipsi(g79(:,:,i),ch079,ps079,ps079) &
                      +v3chibb    (g79(:,:,i),ch079,bz079,bz079))
              endif


           endif

           call insert_val(s1matrix_sm,ssterm(1,1),i1,j1,1)
           call insert_val(d1matrix_sm,ddterm(1,1),i1,j1,1)
           call insert_val(r1matrix_sm,rrterm(1,1),i1,j1,1)
           if(numvar.ge.2) then
               call insert_val(s1matrix_sm,ssterm(1,2),i1  ,j1+6,1)
               call insert_val(s1matrix_sm,ssterm(2,1),i1+6,j1  ,1)
               call insert_val(s1matrix_sm,ssterm(2,2),i1+6,j1+6,1)
               call insert_val(d1matrix_sm,ddterm(1,2),i1  ,j1+6,1)
               call insert_val(d1matrix_sm,ddterm(2,1),i1+6,j1  ,1)
               call insert_val(d1matrix_sm,ddterm(2,2),i1+6,j1+6,1)
               call insert_val(r1matrix_sm,rrterm(1,2),i1  ,j1+6,1)
               call insert_val(r1matrix_sm,rrterm(2,1),i1+6,j1  ,1)
               call insert_val(r1matrix_sm,rrterm(2,2),i1+6,j1+6,1)
           endif
           if(numvar.ge.3) then
              call insert_val(s1matrix_sm,ssterm(1,3),i1,   j1+12,1)
              call insert_val(s1matrix_sm,ssterm(2,3),i1+6, j1+12,1)
              call insert_val(s1matrix_sm,ssterm(3,3),i1+12,j1+12,1)
              call insert_val(s1matrix_sm,ssterm(3,1),i1+12,j1,   1)
              call insert_val(s1matrix_sm,ssterm(3,2),i1+12,j1+6, 1)
              call insert_val(d1matrix_sm,ddterm(1,3),i1,   j1+12,1)
              call insert_val(d1matrix_sm,ddterm(2,3),i1+6, j1+12,1)
              call insert_val(d1matrix_sm,ddterm(3,3),i1+12,j1+12,1)
              call insert_val(d1matrix_sm,ddterm(3,1),i1+12,j1,   1)
              call insert_val(d1matrix_sm,ddterm(3,2),i1+12,j1+6, 1)
              call insert_val(r1matrix_sm,rrterm(1,3),i1,   j1+12,1)
              call insert_val(r1matrix_sm,rrterm(2,3),i1+6, j1+12,1)
              call insert_val(r1matrix_sm,rrterm(3,3),i1+12,j1+12,1)
              call insert_val(r1matrix_sm,rrterm(3,1),i1+12,j1,   1)
              call insert_val(r1matrix_sm,rrterm(3,2),i1+12,j1+6, 1)
           endif
        enddo               ! on j
     enddo                  ! on i

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        telm = telm + tend - tstart
     endif
  enddo                     ! on itri

! modify the s-matrix, inserting the boundary conditions
!
! calculate boundary conditions
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

  if(myrank.eq.0 .and. itimer.eq.1) then
     write(*,*) "  ludefvel_t: Time spent defining fields: ", tfield
     write(*,*) "  ludefvel_t: Time spent calculating matrix: ", telm
  endif

  return
end subroutine ludefvel_t
!============================================================
subroutine ludefphi_t
  use p_data
  use t_data
  use basic
  use arrays
  use sparse_matrix
  use sparse
  use nintegrate_mod
#ifdef mpi
  use supralu_dist_mod
  implicit none
  include 'mpif.h'
#else
  implicit none
#endif
  integer itri, numelms, i, j, l, i1, j1, i2, i3
  real :: x, z, xmin, zmin, factor, pefac
!  real :: deex, hypf
  real :: ssterm(3,3),ddterm(3,3),rrterm(3,3),qqterm(3,3)

  real :: tstart, tend, tfield, telm

  real :: temp
  real :: b3pebd, p1pu, p1pchi

  double precision cogcoords(3)

  tfield = 0.
  telm = 0.

  call numfac(numelms)
  call getmincoord(xmin,zmin)
  
  ! implicit parameter for hyperviscosity terms
  
  ! form field matrices
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

  q4 = 0.

  if(ipres.eq.1) then
     pefac = 1.
  else
     pefac = (p0-pi0)/p0
  endif

  do itri=1,numelms

     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! calculate the field values and derivatives at the sampling points
     call define_fields_79(itri)
     call define_fields_25(itri)

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        tfield = tfield + tend - tstart
     endif
     
!     call getdeex(itri,deex)   

     if(imask.eq.1) then
        call cogfac(itri,cogcoords)
        x = cogcoords(1)-xmin
        z = cogcoords(2)-zmin
        call mask(x,z,factor)
     else
        factor = 1.
     endif
     dbf = db*factor

     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     do i=1,18

        i1 = isvaln(itri,i)
        i2 = i1 + 6
        i3 = i2 + 6

        do j=1,18

           ssterm = 0.
           ddterm = 0.
           rrterm = 0.
           qqterm = 0.

!           hypf = hyper*deex**2

           j1 = isvaln(itri,j)

           call b1psi(itri,i,j,ssterm, ddterm, rrterm, qqterm)         ! passed: 1,2
           call b1psieta(itri,i,j,ssterm, ddterm, rrterm, qqterm)      ! passed: 1
           call b1psiu(itri,i,j,ssterm, ddterm, rrterm, qqterm)        ! passed: 1

           if(numvar.ge.2) then
              call b1psibd(itri,i,j,ssterm, ddterm, rrterm, qqterm)    ! passed: 1

              call b2b(itri,i,j,ssterm, ddterm, rrterm, qqterm)        ! passed: 1
              call b2bu(itri,i,j,ssterm, ddterm, rrterm, qqterm)       ! passed: 1
              call b2beta(itri,i,j,ssterm, ddterm, rrterm, qqterm)     ! passed: 1
              call b2psiv(itri,i,j,ssterm, ddterm, rrterm, qqterm)     ! passed: 1
              call b2psipsid(itri,i,j,ssterm, ddterm, rrterm, qqterm)  ! passed: 1
           endif

           if(numvar.ge.3) then
              call b3pe(itri,i,j,ssterm, ddterm, rrterm, qqterm)       ! passed: 1

              temp = b3pebd(g79(:,:,i),pe179,g79(:,:,j),ni79)*dbf      ! passed: 1

              ssterm(3,2) = ssterm(3,2) - thimp*dt*temp
              ddterm(3,2) = ddterm(3,2) - thimp*dt*temp

              temp = b3pebd(g79(:,:,i),g79(:,:,j),bzt79,ni79)*dbf &    ! passed: 1
                   + p1pu  (g79(:,:,i),g79(:,:,j),pht79) !&            ! passed: 1
!                   + p1pchi(g79(:,:,i),g79(:,:,j),cht79)

              ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
              ddterm(3,3) = ddterm(3,3) + (1.-thimp)*dt*temp

              temp = p1pu(g79(:,:,i),pe179,g79(:,:,j))                 ! passed: 1

              rrterm(3,1) = rrterm(3,1) + thimp*dt*temp
              qqterm(3,1) = qqterm(3,1) - thimp*dt*temp

!!$              temp = p1pchi(g79(:,:,i),pe179,g79(:,:,j))
!!$
!!$              rrterm(3,3) = rrterm(3,3) + thimp*dt*temp
!!$              qqterm(3,3) = qqterm(3,3) - thimp*dt*temp


              if(linear.eq.1 .or. eqsubtract.eq.1) then

                 temp = b3pebd(g79(:,:,i),pe079,g79(:,:,j),ni79)*dbf*pefac
                 
                 ssterm(3,2) = ssterm(3,2) -     thimp *dt*temp
                 ddterm(3,2) = ddterm(3,2) + (1.-thimp)*dt*temp

                 temp = p1pu(g79(:,:,i),pe079,g79(:,:,j))
                 
                 rrterm(3,1) = rrterm(3,1) +     thimp *dt*temp
                 qqterm(3,1) = qqterm(3,1) + (1.-thimp)*dt*temp

                 temp = p1pchi(g79(:,:,i),pe079,g79(:,:,j))
                 
                 rrterm(3,3) = rrterm(3,3) +     thimp *dt*temp
                 qqterm(3,3) = qqterm(3,3) + (1.-thimp)*dt*temp

                 q4(i3) = q4(i3) + dt* &
                      (b3pebd(g79(:,:,i),pe079,bz079,ni79)*dbf*pefac)

                 ! EQUILIBRIUM TERMS
                 q4(i3) = q4(i3) + dt* &
                      (p1pu  (g79(:,:,i),pe079,ph079) &
                      +p1pchi(g79(:,:,i),pe079,ch079))
              endif
           endif

           call insert_val(s2matrix_sm,ssterm(1,1),i1,j1,1)
           call insert_val(d2matrix_sm,ddterm(1,1),i1,j1,1)
           call insert_val(r2matrix_sm,rrterm(1,1),i1,j1,1)
           call insert_val(q2matrix_sm,qqterm(1,1),i1,j1,1)
           if(numvar.ge.2) then
              call insert_val(s2matrix_sm,ssterm(1,2),i1  ,j1+6,1)
              call insert_val(s2matrix_sm,ssterm(2,1),i1+6,j1,1)
              call insert_val(s2matrix_sm,ssterm(2,2),i1+6,j1+6,1)
              call insert_val(d2matrix_sm,ddterm(1,2),i1  ,j1+6,1)
              call insert_val(d2matrix_sm,ddterm(2,1),i1+6,j1,1)
              call insert_val(d2matrix_sm,ddterm(2,2),i1+6,j1+6,1)
              call insert_val(r2matrix_sm,rrterm(1,2),i1  ,j1+6,1)
              call insert_val(r2matrix_sm,rrterm(2,1),i1+6,j1,1)
              call insert_val(r2matrix_sm,rrterm(2,2),i1+6,j1+6,1)
              call insert_val(q2matrix_sm,qqterm(1,2),i1  ,j1+6,1)
              call insert_val(q2matrix_sm,qqterm(2,1),i1+6,j1,1)
              call insert_val(q2matrix_sm,qqterm(2,2),i1+6,j1+6,1)
           endif
           if(numvar .eq. 3) then
              call insert_val(s2matrix_sm,ssterm(1,3),i1,   j1+12,1)
              call insert_val(s2matrix_sm,ssterm(2,3),i1+6, j1+12,1)
              call insert_val(s2matrix_sm,ssterm(3,3),i1+12,j1+12,1)
              call insert_val(s2matrix_sm,ssterm(3,1),i1+12,j1,   1)
              call insert_val(s2matrix_sm,ssterm(3,2),i1+12,j1+6, 1)
              call insert_val(d2matrix_sm,ddterm(1,3),i1,   j1+12,1)
              call insert_val(d2matrix_sm,ddterm(2,3),i1+6, j1+12,1)
              call insert_val(d2matrix_sm,ddterm(3,3),i1+12,j1+12,1)
              call insert_val(d2matrix_sm,ddterm(3,1),i1+12,j1,   1)
              call insert_val(d2matrix_sm,ddterm(3,2),i1+12,j1+6, 1)
              call insert_val(r2matrix_sm,rrterm(1,3),i1,   j1+12,1)
              call insert_val(r2matrix_sm,rrterm(2,3),i1+6, j1+12,1)
              call insert_val(r2matrix_sm,rrterm(3,3),i1+12,j1+12,1)
              call insert_val(r2matrix_sm,rrterm(3,1),i1+12,j1,   1)
              call insert_val(r2matrix_sm,rrterm(3,2),i1+12,j1+6, 1)
              call insert_val(q2matrix_sm,qqterm(1,3),i1,   j1+12,1)
              call insert_val(q2matrix_sm,qqterm(2,3),i1+6, j1+12,1)
              call insert_val(q2matrix_sm,qqterm(3,3),i1+12,j1+12,1)
              call insert_val(q2matrix_sm,qqterm(3,1),i1+12,j1,   1)
              call insert_val(q2matrix_sm,qqterm(3,2),i1+12,j1+6, 1)
           endif
        enddo
     enddo

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        telm = telm + tend - tstart
     endif
  enddo
  
  ! modify the s-matrix, inserting the boundary conditions
  !
  ! calculate boundary conditions
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

  if(myrank.eq.0 .and. itimer.eq.1) then
     write(*,*) "  ludefphi: Time spent defining fields: ", tfield
     write(*,*) "  ludefphi: Time spent calculating matrix: ", telm
  endif
  
  return
end subroutine ludefphi_t
!================================
subroutine ludefden_t
  use p_data
  use t_data
  use basic
  use arrays
  use sparse_matrix
  use sparse
  use nintegrate_mod

#ifdef mpi
  use supralu_dist_mod
  implicit none
  include 'mpif.h'
#else
  implicit none
#endif
  
  integer :: itri
  integer :: i, j, j1, i1, ione, jone, numelms
  real :: xmin, zmin
  real :: ssterm,ddterm,rrterm(3),qqterm(3)
 
  call numfac(numelms)
  call getmincoord(xmin,zmin)

  ! form field matrices
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

#endif

  do itri=1,numelms

     ! calculate the field values and derivatives at the sampling points
     call define_fields_79(itri)
     call define_fields_25(itri)

!!$     call getdeex(itri,deex)

     do i=1,18
        ione = isval1(itri,i)
        i1 = isvaln(itri,i)

        do j=1,18         
           ssterm = 0.
           ddterm = 0.
           rrterm = 0.
           qqterm = 0.

           jone = isval1(itri,j)
           j1 = isvaln(itri,j)

           call n1n(itri,i,j,ssterm,ddterm,rrterm,qqterm)
           call n1nu(itri,i,j,ssterm,ddterm,rrterm,qqterm)

           call insert_val(s8matrix_sm, ssterm, ione, jone, 1)
           call insert_val(d8matrix_sm, ddterm, ione, jone, 1)
           call insert_val(r8matrix_sm, rrterm(1), i1, j1, 1)
           call insert_val(q8matrix_sm, qqterm(1), i1, j1, 1)
           if(numvar.ge.2) then
              call insert_val(r8matrix_sm,rrterm(2), i1,j1+6,1)
              call insert_val(q8matrix_sm,qqterm(2), i1,j1+6,1)
           endif
           if(numvar.ge.3) then
              call insert_val(r8matrix_sm,rrterm(3), i1,j1+12,1)
              call insert_val(q8matrix_sm,qqterm(3), i1,j1+12,1)
           endif

        enddo                     ! on i
     enddo                     ! on j
  enddo                     ! on itri
  !
  ! modify the s-matrix, inserting the boundary conditions
  !
  ! calculate boundary conditions
  call boundaryds(iboundn,nbcn,1)
  if(nbcn .gt. iboundmax) then
     write(*,4888) nbcn, iboundmax
4888 format(" ERROR: nbcn > iboundmax", 2i5)
     call safestop(9) 
  endif
  
  !  bc stuff...
  do i=1,nbcn
     call set_diri_bc(s8matrix_sm, iboundn(i))
  enddo

  call finalize_array(s8matrix_sm)
  call finalize_array(d8matrix_sm)
  call finalize_array(q8matrix_sm)
  call finalize_array(r8matrix_sm)
  
  return
end subroutine ludefden_t

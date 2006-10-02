!=====================================
subroutine ludefvel_t
  use p_data
  use t_data
  use basic
  use arrays
  use sparse_matrix
  use sparse
  use nintegrate_mod
  use metricterms_n

#ifdef mpi
  use supralu_dist_mod
  implicit none
  include 'mpif.h'
#else
  implicit none
#endif
  integer itri, l, i, j, j1, i1, i2, i3, numelms
  real :: ssterm(3,3),ddterm(3,3),rrterm(3,3)
  real :: x, xmin, z, zmin, factor, temp, dbf

  double precision :: cogcoords(3)

  real :: tstart, tend, tfield, telm

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
!     call define_fields_25(itri)
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
           temp = v1un(g79(:,:,i),g79(:,:,j),nt79)                  ! passed: 1,2
           ssterm(1,1) = ssterm(1,1) + temp
           ddterm(1,1) = ddterm(1,1) + temp  

           temp = v1umu(g79(:,:,i),g79(:,:,j))*amu &                ! passed: 1,2
                + thimp*dt* &
                (v1upsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79))      ! passed: 1,2
           ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
           ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

           temp = v1uun(g79(:,:,i),g79(:,:,j),ph179,nt79) &         ! passed: 1,2
                + v1uun(g79(:,:,i),ph179,g79(:,:,j),nt79)           ! passed: 1,2
           ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
           ddterm(1,1) = ddterm(1,1) + (.5-thimp)*dt*temp

           rrterm(1,1) = rrterm(1,1) + dt* &
                (v1psipsi(g79(:,:,i),g79(:,:,j),pss79)  &           ! passed: 1,2
                +v1psipsi(g79(:,:,i),pss79,g79(:,:,j))) &           ! passed: 1,2
                + thimp*dt*dt* &
                (v1psisb1 (g79(:,:,i),g79(:,:,j),sb179))

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = v1uun(g79(:,:,i),g79(:,:,j),ph079,nt79) &      ! passed: 1,2
                   + v1uun(g79(:,:,i),ph079,g79(:,:,j),nt79)        ! passed: 1,2
              ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
              ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

              rrterm(1,1) = rrterm(1,1) + thimp*dt*dt* &
                (v1upsipsi(g79(:,:,i),ph079,g79(:,:,j),pss79) &     ! passed: 1,2
                +v1upsipsi(g79(:,:,i),ph079,pss79,g79(:,:,j)))      ! passed: 1,2
           endif

           ! NUMVAR = 2
           ! ~~~~~~~~~~
           if(numvar.ge.2) then

              temp = v1ubb(g79(:,:,i),g79(:,:,j),bzt79,bzt79)
              ssterm(1,1) = ssterm(1,1) - thimp*    thimp *dt*dt*temp
              ddterm(1,1) = ddterm(1,1) + thimp*(1.-thimp)*dt*dt*temp

              temp = v1vpsib(g79(:,:,i),g79(:,:,j),pst79,bzt79)
              ssterm(1,2) = ssterm(1,2) - thimp*    thimp *dt*dt*temp
              ddterm(1,2) = ddterm(1,2) + thimp*(1.-thimp)*dt*dt*temp

              temp = v2vun(g79(:,:,i),vz179,g79(:,:,j),nt79)
              ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
              ddterm(2,1) = ddterm(2,1) + (.5-thimp)*dt*temp

              temp = v2upsib(g79(:,:,i),g79(:,:,j),pst79,bzt79)
              ssterm(2,1) = ssterm(2,1) - thimp*    thimp *dt*dt*temp
              ddterm(2,1) = ddterm(2,1) + thimp*(1.-thimp)*dt*dt*temp

              temp = v2vn(g79(:,:,i),g79(:,:,j),nt79)
              ssterm(2,2) = ssterm(2,2) + temp
              ddterm(2,2) = ddterm(2,2) + temp

              temp = v2vmu(g79(:,:,i),g79(:,:,j))*amu &
                   + thimp*dt* &
                   (v2vpsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79))
              ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

              temp = v2vun(g79(:,:,i),g79(:,:,j),ph179,nt79)
              ssterm(2,2) = ssterm(2,2) -      thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) + (0.5-thimp)*dt*temp

              rrterm(1,2) = rrterm(1,2) + thimp*dt*dt* &
                   (v1bsb2   (g79(:,:,i),g79(:,:,j),sb279))

              rrterm(2,1) = rrterm(2,1) + dt* &
                   (v2psib(g79(:,:,i),g79(:,:,j),bzs79)) &
                   + thimp*dt*dt* &
                   (v2psib   (g79(:,:,i),g79(:,:,j),sb279))

              rrterm(2,2) = rrterm(2,2) + dt* &
                   (v2psib(g79(:,:,i),pss79,g79(:,:,j))) &
                   + thimp*dt*dt* &
                   (v2psib(g79(:,:,i),sb179,g79(:,:,j)))

              if(linear.eq.1 .or. eqsubtract.eq.1) then
                 
                 temp = v2vun(g79(:,:,i),vz079,g79(:,:,j),nt79)
                 ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
                 ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp

                 temp = v2vun(g79(:,:,i),g79(:,:,j),ph079,nt79)
                 ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
                 ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

                 rrterm(1,1) = rrterm(1,1) + thimp*dt*dt* &
                      (v1vpsib  (g79(:,:,i),vz079,g79(:,:,j),bzs79))

                 rrterm(1,2) = rrterm(1,2) + dt* &
                      (v1bb     (g79(:,:,i),g79(:,:,j),bzs79) &
                      +v1bb     (g79(:,:,i),bzs79,g79(:,:,j))) &
                      + thimp*dt*dt* &
                      (v1ubb    (g79(:,:,i),ph079,g79(:,:,j),bzs79) &
                      +v1ubb    (g79(:,:,i),ph079,bzs79,g79(:,:,j)) &
                      +v1vpsib  (g79(:,:,i),vz079,pss79,g79(:,:,j)))
                 
                 rrterm(2,1) = rrterm(2,1) + thimp*dt*dt* &
                      (v2upsib  (g79(:,:,i),ph079,g79(:,:,j),bzs79) &
                      +v2vpsipsi(g79(:,:,i),vz079,g79(:,:,j),pss79) &
                      +v2vpsipsi(g79(:,:,i),vz079,pss79,g79(:,:,j)))

                 rrterm(2,2) = rrterm(2,2) + thimp*dt*dt* &
                      (v2upsib  (g79(:,:,i),ph079,pss79,g79(:,:,j)))

              endif
           endif

           ! NUMVAR = 3
           ! ~~~~~~~~~~
           if(numvar.ge.3) then

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

              temp = v1chipsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79)  & ! passed: 1
                   + v1chibb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)    ! passed: 1
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
              
              temp = v3un(g79(:,:,i),g79(:,:,j),nt79)
              ssterm(3,1) = ssterm(3,1) + temp
              ddterm(3,1) = ddterm(3,1) + temp

              temp = v3uun  (g79(:,:,i),g79(:,:,j),ph179,nt79) &      ! passed: 1
                   + v3uun  (g79(:,:,i),ph179,g79(:,:,j),nt79) !&      ! passed: 1
!!$                   + v3uchin(g79(:,:,i),g79(:,:,j),ch179,nt79)          ! FAILED: 1
              ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
              ddterm(3,1) = ddterm(3,1) + (.5-thimp)*dt*temp

              temp = v3umu   (g79(:,:,i),g79(:,:,j)) &
                   +thimp*dt* &
                   (v3up     (g79(:,:,i),g79(:,:,j),pt79) &               ! 
!!$                   +v3upsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79) &        !
                   +v3ubb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79))          ! 
              ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
              ddterm(3,1) = ddterm(3,1) + (1.-thimp)*dt*temp

              temp = v3vpsib(g79(:,:,i),g79(:,:,j),pst79,bzt79)        ! passed: 1
              ssterm(3,2) = ssterm(3,2) - thimp*    thimp *dt*dt*temp
              ddterm(3,2) = ddterm(3,2) + thimp*(1.-thimp)*dt*dt*temp

              temp = v3chin(g79(:,:,i),g79(:,:,j),nt79)        ! passed: 1
              ssterm(3,3) = ssterm(3,3) + temp
              ddterm(3,3) = ddterm(3,3) + temp

              temp = v3chimu    (g79(:,:,i),g79(:,:,j))        ! passed: 1
              ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
              ddterm(3,3) = ddterm(3,3) + (1.-thimp)*dt*temp

              temp = 0*v3uchin  (g79(:,:,i),ph179,g79(:,:,j),nt79) &    ! FAILED: 1
                   + v3chichin(g79(:,:,i),g79(:,:,j),ch179,nt79) &    ! passed: 1
                   + v3chichin(g79(:,:,i),ch179,g79(:,:,j),nt79)      ! passed: 1
              ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
              ddterm(3,3) = ddterm(3,3) + (.5-thimp)*dt*temp
                           
              temp =  v3chip     (g79(:,:,i),g79(:,:,j),pt79) &         ! passed: 1
!!$                   + v3chipsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79) &  ! FAILED: 1
                   + v3chibb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)     ! passed: 1
              ssterm(3,3) = ssterm(3,3) - thimp*    thimp *dt*dt*temp
              ddterm(3,3) = ddterm(3,3) + thimp*(1.-thimp)*dt*dt*temp

              rrterm(3,1) = rrterm(3,1) +                    &
                   dt*                                       &
                   (v3psipsi(g79(:,:,i),g79(:,:,j),pss79)    &       ! passed: 1
                   +v3psipsi(g79(:,:,i),pss79,g79(:,:,j))) + &       ! passed: 1
                   thimp*dt*dt*                              &
                   (v3psisb1(g79(:,:,i),g79(:,:,j),sb179)    &       ! passed: 1  
                   +v3psisb1(g79(:,:,i),sb179,g79(:,:,j)))           ! passed: 1

              rrterm(3,2) = rrterm(3,2) + dt*            &
                   (v3bb(g79(:,:,i),g79(:,:,j),bzs79)    &           ! passed: 1
                   +v3bb(g79(:,:,i),bzs79,g79(:,:,j))) + &           ! passed: 1
                   thimp*dt*dt*                          &
                   (v3bsb2(g79(:,:,i),g79(:,:,j),sb279))             ! passed: 1
              
              rrterm(3,3) = rrterm(3,3) + dt* &
                   v3p(g79(:,:,i),g79(:,:,j))                        ! passed: 1

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

        ! Definition of R4
        ! ================
        if(linear.eq.1 .or. eqsubtract.eq.1) then
           r4(i1) = r4(i1) + dt* &
                (v1umu    (g79(:,:,i),ph079)) &                   ! passed: 1, 2
                + thimp*dt*dt* &
                (v1psisb1 (g79(:,:,i),ps079,sb179))               ! passed: 1, 2

           ! DENSITY TERMS
           r4(i1) = r4(i1) + dt* &
                (v1un     (g79(:,:,i),ph079,nt79)         &       ! passed: 1, 2 
                +v1uun    (g79(:,:,i),ph079,ph079,nt79))          ! passed: 1, 2

           ! EQUILIBRIUM TERMS
           r4(i1) = r4(i1) + dt* &
                (v1psipsi (g79(:,:,i),ps079,ps079)) &             ! passed: 1, 2
                + thimp*dt*dt* &
                (v1upsipsi(g79(:,:,i),ph079,ps079,ps079))         ! passed: 1, 2

           if(numvar.ge.2) then
              r4(i1) = r4(i1) + thimp*dt*dt* &
                   (v1bsb2(g79(:,:,i),bz079,sb279))
              r4(i2) = r4(i2) + dt* &
                   (v2vmu(g79(:,:,i),vz079)*amu)
              
              ! DENSITY TERMS
              r4(i2) = r4(i2) + dt* &
                   (v2vun(g79(:,:,i),vz079,ph079,nt79)) &
                   + thimp*dt*dt* &
                   (v2psib(g79(:,:,i),sb179,bz079) &
                   +v2psib(g79(:,:,i),ps079,sb279))
              
              ! EQUILIBRIUM TERM
              r4(i1) = r4(i1) + dt* &
                   (v1bb     (g79(:,:,i),bz079,bz079)) &
                   + thimp*dt*dt* &
                   (v1ubb    (g79(:,:,i),ph079,bz079,bz079) &
                   +v1vpsib  (g79(:,:,i),vz079,ps079,bz079))
              r4(i2) = r4(i2) + dt* &
                   (v2psib   (g79(:,:,i),ps079,bz079)) &
                   + thimp*dt*dt* &
                   (v2upsib  (g79(:,:,i),ph079,ps079,bz079) &
                   +v2vpsipsi(g79(:,:,i),vz079,ps079,ps079))
           endif

           if(numvar.ge.3) then
              
              r4(i3) = r4(i3) + &
                   dt* &
                   (v3umu   (g79(:,:,i),ph079)*amu &
                   +v3chimu (g79(:,:,i),ch079)*amu) + &
                   thimp*dt*dt* &
                   (v3psisb1(g79(:,:,i),ps079,sb179) &
                   +v3psisb1(g79(:,:,i),sb179,ps079) &
                   +v3bsb2  (g79(:,:,i),bz079,sb279))

              ! DENSITY TERMS
              r4(i1) = r4(i1) + dt* &
                   (v1uchin  (g79(:,:,i),ph079,ch079,nt79) &
                   +v1chichin(g79(:,:,i),ch079,ch079,nt79))
              r4(i3) = r4(i3) + dt* &
                   (v3uun    (g79(:,:,i),ph079,ph079,nt79) &
                   +v3uchin  (g79(:,:,i),ph079,ch079,nt79) &
                   +v3chichin(g79(:,:,i),ch079,ch079,nt79))
              
              ! EQUILIBRIUM TERMS
              r4(i1) = r4(i1) + thimp*dt*dt* &
                    (v1chipsipsi(g79(:,:,i),ch079,ps079,ps079) &
                    +v1chibb    (g79(:,:,i),ch079,bz079,bz079))
              r4(i2) = r4(i2) + &
                   thimp*dt*dt* &
                   (v2chipsib  (g79(:,:,i),ch079,ps079,bz079))
              r4(i3) = r4(i3) + &
                   dt* &
                   (v3p        (g79(:,:,i),p079)              &
                   +v3psipsi   (g79(:,:,i),ps079,ps079)       &
                   +v3bb       (g79(:,:,i),bz079,bz079)) +    &
                   thimp*dt*dt* &
                   (v3up       (g79(:,:,i),ph079,p079)        &
                   +v3upsipsi  (g79(:,:,i),ph079,ps079,ps079) &
                   +v3ubb      (g79(:,:,i),ph079,bz079,bz079) &
                   +v3vpsib    (g79(:,:,i),vz079,ps079,bz079) &
                   +v3chip     (g79(:,:,i),ch079,p079)        &
                   +v3chipsipsi(g79(:,:,i),ch079,ps079,ps079) &
                   +v3chibb    (g79(:,:,i),ch079,bz079,bz079))
              
           endif ! on numvar.ge.3
        endif    ! on linear.eq.1 .or. eqsubtract.eq.1

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
  use metricterms_n

#ifdef mpi
  use supralu_dist_mod
  implicit none
  include 'mpif.h'
#else
  implicit none
#endif
  integer itri, numelms, i, j, l, i1, j1, i2, i3
  real :: x, z, xmin, zmin, factor, pefac, dbf
!  real :: deex, hypf
  real :: ssterm(3,3),ddterm(3,3),rrterm(3,3),qqterm(3,3)

  real :: tstart, tend, tfield, telm

  real :: temp
!!$  real :: b3pebd, p1pu, p1pchi, b1psichi, b2ped, b2bchi

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
!     call define_fields_25(itri)

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

           ! NUMVAR = 1
           ! ~~~~~~~~~~
           temp = b1psi(g79(:,:,i),g79(:,:,j))
           ssterm(1,1) = ssterm(1,1) + temp
           ddterm(1,1) = ddterm(1,1) + temp

           temp = b1psieta(g79(:,:,i),g79(:,:,j))*etar &
                + b1psiu  (g79(:,:,i),g79(:,:,j),pht79)
           ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
           ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

           temp = b1psiu(g79(:,:,i),ps179,g79(:,:,j))
           rrterm(1,1) = rrterm(1,1) + thimp*dt*temp
           qqterm(1,1) = qqterm(1,1) - thimp*dt*temp

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = b1psiu(g79(:,:,i),ps079,g79(:,:,j))
              rrterm(1,1) = rrterm(1,1) +     thimp *dt*temp
              qqterm(1,1) = qqterm(1,1) + (1.-thimp)*dt*temp
           endif

           ! NUMVAR = 2
           ! ~~~~~~~~~~
           if(numvar.ge.2) then
              
              temp = b1psibd(g79(:,:,i),g79(:,:,j),bz179,ni79)*dbf              
              ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
              ddterm(1,1) = ddterm(1,1) + (.5-thimp)*dt*temp

              temp = b1psibd(g79(:,:,i),ps179,g79(:,:,j),ni79)*dbf              
              ssterm(1,2) = ssterm(1,2) -     thimp *dt*temp
              ddterm(1,2) = ddterm(1,2) + (.5-thimp)*dt*temp

              temp = b2psiv   (g79(:,:,i),g79(:,:,j),vzt79)
              ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
              ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp

              temp = b2psipsid(g79(:,:,i),g79(:,:,j),ps179,ni79)*dbf &
                   + b2psipsid(g79(:,:,i),ps179,g79(:,:,j),ni79)*dbf
              ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
              ddterm(2,1) = ddterm(2,1) + (.5-thimp)*dt*temp

              temp = b2b(g79(:,:,i),g79(:,:,j))
              ssterm(2,2) = ssterm(2,2) + temp
              ddterm(2,2) = ddterm(2,2) + temp

              temp = b2beta(g79(:,:,i),g79(:,:,j))*etar &
                   + b2bu  (g79(:,:,i),g79(:,:,j),pht79)
              ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

              temp = b2bu  (g79(:,:,i),bz179,g79(:,:,j))
              rrterm(2,1) = rrterm(2,1) + thimp*dt*temp
              qqterm(2,1) = qqterm(2,1) - thimp*dt*temp

              temp = b2psiv(g79(:,:,i),ps179,g79(:,:,j))
              rrterm(2,2) = rrterm(2,2) + thimp*dt*temp
              qqterm(2,2) = qqterm(2,2) - thimp*dt*temp

              if(linear.eq.1 .or. eqsubtract.eq.1) then
                 temp = b1psibd(g79(:,:,i),g79(:,:,j),bz079,ni79)*dbf              
                 ssterm(1,1) = ssterm(1,1) -      thimp *dt*temp
                 ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

                 temp = b1psibd(g79(:,:,i),ps079,g79(:,:,j),ni79)*dbf              
                 ssterm(1,2) = ssterm(1,2) -     thimp *dt*temp
                 ddterm(1,2) = ddterm(1,2) + (1.-thimp)*dt*temp

                 temp = b2psipsid(g79(:,:,i),g79(:,:,j),ps079,ni79)*dbf &
                      + b2psipsid(g79(:,:,i),ps079,g79(:,:,j),ni79)*dbf
                 ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
                 ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp

                 temp = b2bu  (g79(:,:,i),bz079,g79(:,:,j))
                 rrterm(2,1) = rrterm(2,1) +     thimp *dt*temp
                 qqterm(2,1) = qqterm(2,1) + (1.-thimp)*dt*temp

                 temp = b2psiv(g79(:,:,i),ps079,g79(:,:,j))
                 rrterm(2,2) = rrterm(2,2) +     thimp *dt*temp
                 qqterm(2,2) = qqterm(2,2) + (1.-thimp)*dt*temp
              endif
           endif

           ! NUMVAR = 3
           ! ~~~~~~~~~~
           if(numvar.ge.3) then

              temp = b1psichi(g79(:,:,i),g79(:,:,j),cht79)             ! passed: 1
              ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
              ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

              temp = b1psichi(g79(:,:,i),ps179,g79(:,:,j))             ! passed: 1
              rrterm(1,3) = rrterm(1,3) + thimp*dt*temp
              qqterm(1,3) = qqterm(1,3) - thimp*dt*temp

              temp = b2bchi(g79(:,:,i),g79(:,:,j),cht79)               ! passed: 1
              ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

              temp = b2ped(g79(:,:,i),g79(:,:,j),ni79)*dbf             ! passed: 1
              ssterm(2,3) = ssterm(2,3) -     thimp *dt*temp
              ddterm(2,3) = ddterm(2,3) - (1.-thimp)*dt*temp

              temp = b2bchi(g79(:,:,i),bz179,g79(:,:,j))               ! passed: 1
              rrterm(2,3) = ssterm(2,3) + thimp*dt*temp
              qqterm(2,3) = qqterm(2,3) - thimp*dt*temp

              temp = b3pebd(g79(:,:,i),pe179,g79(:,:,j),ni79)*dbf*pefac   ! passed: 1
              ssterm(3,2) = ssterm(3,2) - thimp*dt*temp
              ddterm(3,2) = ddterm(3,2) - thimp*dt*temp

              temp = b3pe(g79(:,:,i),g79(:,:,j))
              ssterm(3,3) = ssterm(3,3) + temp
              ddterm(3,3) = ddterm(3,3) + temp

              temp = b3pebd(g79(:,:,i),g79(:,:,j),bzt79,ni79)*dbf*pefac & ! passed: 1
                   + p1pu  (g79(:,:,i),g79(:,:,j),pht79)          &       ! passed: 1
                   + p1pchi(g79(:,:,i),g79(:,:,j),cht79)                  ! passed: 1
              ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
              ddterm(3,3) = ddterm(3,3) + (1.-thimp)*dt*temp

              temp = p1pu(g79(:,:,i),pe179,g79(:,:,j))                    ! passed: 1
              rrterm(3,1) = rrterm(3,1) + thimp*dt*temp
              qqterm(3,1) = qqterm(3,1) - thimp*dt*temp

              temp = p1pchi(g79(:,:,i),pe179,g79(:,:,j))                  ! passed: 1
              rrterm(3,3) = rrterm(3,3) + thimp*dt*temp
              qqterm(3,3) = qqterm(3,3) - thimp*dt*temp

              if(linear.eq.1 .or. eqsubtract.eq.1) then
                 temp = b1psichi(g79(:,:,i),ps079,g79(:,:,j))
                 rrterm(1,3) = rrterm(1,3) +     thimp *dt*temp
                 qqterm(1,3) = qqterm(1,3) + (1.-thimp)*dt*temp

                 temp = b2bchi(g79(:,:,i),bz079,g79(:,:,j))
                 rrterm(2,3) = rrterm(2,3) +     thimp *dt*temp
                 qqterm(2,3) = qqterm(2,3) + (1.-thimp)*dt*temp

                 temp = b3pebd(g79(:,:,i),pe079,g79(:,:,j),ni79)*dbf*pefac
                 ssterm(3,2) = ssterm(3,2) -     thimp *dt*temp
                 ddterm(3,2) = ddterm(3,2) + (1.-thimp)*dt*temp

                 temp = p1pu(g79(:,:,i),pe079,g79(:,:,j))
                 rrterm(3,1) = rrterm(3,1) +     thimp *dt*temp
                 qqterm(3,1) = qqterm(3,1) + (1.-thimp)*dt*temp

                 temp = p1pchi(g79(:,:,i),pe079,g79(:,:,j))                
                 rrterm(3,3) = rrterm(3,3) +     thimp *dt*temp
                 qqterm(3,3) = qqterm(3,3) + (1.-thimp)*dt*temp
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
        enddo ! on j

        if(linear.eq.1 .or. eqsubtract.eq.1) then

           q4(i1) = q4(i1) + dt* &
                (b1psieta(g79(:,:,i),ps079)*etar)

           ! EQUILIBRIUM TERMS
           q4(i1) = q4(i1) + dt* &
                (b1psiu(g79(:,:,i),ps079,ph079))

           if(numvar.ge.2) then
              q4(i2) = q4(i2) + dt* &
                   b2beta (g79(:,:,i),bz079)*etar
              
              ! DENSITY TERMS
              q4(i1) = q4(i1) + dt* &
                   b1psibd(g79(:,:,i),ps079,bz079,ni79)*dbf
              q4(i2) = q4(i2) + dt* &
                   b2psipsid(g79(:,:,i),ps079,ps079,ni79)*dbf

              ! EQUILIBRIUM TERMS
              q4(i2) = q4(i2) + dt* &
                   (b2psiv(g79(:,:,i),ps079,vz079) &
                   +b2bu  (g79(:,:,i),bz079,ph079))
           endif

           if(numvar.ge.3) then
              ! DENSITY TERMS
              q4(i2) = q4(i2) + dt* &
                   (b2ped (g79(:,:,i),pe079,ni79)*dbf)
              q4(i3) = q4(i3) + dt* &
                   (b3pebd(g79(:,:,i),pe079,bz079,ni79)*dbf*pefac)
              
              ! EQUILIBRIUM TERMS
              q4(i1) = q4(i1) + dt* &
                   (b1psichi(g79(:,:,i),ps079,ch079))
              q4(i2) = q4(i2) + dt* &
                   (b2bchi(g79(:,:,i),bz079,ch079))
              q4(i3) = q4(i3) + dt* &
                   (p1pu  (g79(:,:,i),pe079,ph079) &
                   +p1pchi(g79(:,:,i),pe079,ch079))
           endif
        endif

     enddo ! on i

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
  use metricterms_n

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
  real :: temp
 
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

!!$  en1 = 0

  do itri=1,numelms

     ! calculate the field values and derivatives at the sampling points
     call define_fields_79(itri)
!     call define_fields_25(itri)

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

           ! NUMVAR = 1
           ! ~~~~~~~~~~
           temp = n1n(g79(:,:,i),g79(:,:,j))
           ssterm = ssterm + temp
           ddterm = ddterm + temp

           temp = n1nu  (g79(:,:,i),g79(:,:,j),pht79)
           ssterm = ssterm -     thimp *dt*temp
           ddterm = ddterm + (1.-thimp)*dt*temp

           temp = n1nu(g79(:,:,i),n179,g79(:,:,j))
           rrterm(1) = rrterm(1) + thimp*dt*temp
           qqterm(1) = qqterm(1) - thimp*dt*temp

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = n1nu  (g79(:,:,i),n079,g79(:,:,j))
              rrterm(1) = rrterm(1) +     thimp *dt*temp
              qqterm(1) = qqterm(1) + (1.-thimp)*dt*temp
!!$              ! EQUILIBRIUM TERMS
!!$              en1(ione) = en1(ione) + dt* &
!!$                   (n1nu  (g79(:,:,i),n079,ph079))
           endif

           ! NUMVAR = 3
           ! ~~~~~~~~~~
           if(numvar.ge.3) then
              temp = n1nchi(g79(:,:,i),g79(:,:,j),cht79)
              ssterm = ssterm -     thimp *dt*temp
              ddterm = ddterm + (1.-thimp)*dt*temp

              temp = n1nchi(g79(:,:,i),n179,g79(:,:,j))
              rrterm(3) = rrterm(3) + thimp*dt*temp
              qqterm(3) = qqterm(3) - thimp*dt*temp

              if(linear.eq.1 .or. eqsubtract.eq.1) then
                 temp = n1nchi(g79(:,:,i),n079,g79(:,:,j))
                 rrterm(3) = rrterm(3) +     thimp *dt*temp
                 qqterm(3) = qqterm(3) + (1.-thimp)*dt*temp
                 
!!$              ! EQUILIBRIUM TERMS
!!$              en1(ione) = en1(ione) + dt* &
!!$                   (n1nchi(g79(:,:,i),n079,ch079))
              endif
           endif

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

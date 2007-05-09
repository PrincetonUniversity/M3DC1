subroutine ludefall

  use p_data
  use t_data
  use basic
  use arrays
  use sparse
  use nintegrate_mod
  use diagnostics

  implicit none

  include 'mpif.h'


  integer :: itri, numelms, i
  integer :: def_fields
  real :: x, z, xmin, zmin, factor, dbf

  real :: tstart, tend, tfield, telm

  double precision cogcoords(3)

  tfield = 0.
  telm = 0.

  call numfac(numelms)
  call getmincoord(xmin,zmin)
  
  ! form field matrices
  if(idens.eq.1) then
     call zeroarray4solve(s8matrix_sm,numvar1_numbering)
     call zeroarray4multiply(d8matrix_sm,numvar1_numbering)
     if(integrator.eq.1) &
          call zeroarray4multiply(o8matrix_sm,numvar1_numbering)
  endif
  if(numvar .eq. 1) then
     call zeroarray4solve(s1matrix_sm,numvar1_numbering)
     call zeroarray4multiply(d1matrix_sm,numvar1_numbering)
     call zeroarray4multiply(r1matrix_sm,numvar1_numbering)
     call zeroarray4solve(s2matrix_sm,numvar1_numbering)
     call zeroarray4multiply(d2matrix_sm,numvar1_numbering)
     call zeroarray4multiply(r2matrix_sm,numvar1_numbering)
     call zeroarray4multiply(q2matrix_sm,numvar1_numbering)
     if(integrator.eq.1) then
        call zeroarray4multiply(o1matrix_sm,numvar1_numbering)
        call zeroarray4multiply(o2matrix_sm,numvar1_numbering)
     endif
     if(idens.eq.1) then
        call zeroarray4multiply(q8matrix_sm,numvar1_numbering)
        call zeroarray4multiply(r8matrix_sm,numvar1_numbering)
     endif
  else if(numvar .eq. 2) then
     call zeroarray4solve(s1matrix_sm,numvar2_numbering)
     call zeroarray4multiply(d1matrix_sm,numvar2_numbering)
     call zeroarray4multiply(r1matrix_sm,numvar2_numbering)
     call zeroarray4solve(s2matrix_sm,numvar2_numbering)
     call zeroarray4multiply(d2matrix_sm,numvar2_numbering)
     call zeroarray4multiply(r2matrix_sm,numvar2_numbering)
     call zeroarray4multiply(q2matrix_sm,numvar2_numbering)
     if(integrator.eq.1) then
        call zeroarray4multiply(o1matrix_sm,numvar2_numbering)
        call zeroarray4multiply(o2matrix_sm,numvar2_numbering)
     endif
     if(idens.eq.1) then
        call zeroarray4multiply(q8matrix_sm,numvar2_numbering)
        call zeroarray4multiply(r8matrix_sm,numvar2_numbering)
     endif
  else 
     call zeroarray4solve(s1matrix_sm,numvar3_numbering)
     call zeroarray4multiply(d1matrix_sm,numvar3_numbering)
     call zeroarray4multiply(r1matrix_sm,numvar3_numbering)
     call zeroarray4solve(s2matrix_sm,numvar3_numbering)
     call zeroarray4multiply(d2matrix_sm,numvar3_numbering)
     call zeroarray4multiply(r2matrix_sm,numvar3_numbering)
     call zeroarray4multiply(q2matrix_sm,numvar3_numbering)
     if(integrator.eq.1) then
        call zeroarray4multiply(o1matrix_sm,numvar3_numbering)
        call zeroarray4multiply(o2matrix_sm,numvar3_numbering)
     endif
     if(idens.eq.1) then
        call zeroarray4multiply(q8matrix_sm,numvar3_numbering)
        call zeroarray4multiply(r8matrix_sm,numvar3_numbering)
     endif
  endif
  if(ipres.eq.1) then
     call zeroarray4solve(s9matrix_sm,numvar1_numbering)
     call zeroarray4multiply(d9matrix_sm,numvar1_numbering)
     call zeroarray4multiply(q9matrix_sm,numvar3_numbering)
     call zeroarray4multiply(r9matrix_sm,numvar3_numbering)
     if(integrator.eq.1) &
          call zeroarray4multiply(o9matrix_sm,numvar1_numbering)
  endif
  
  r4 = 0.
  q4 = 0.
  if(idens.eq.1) qn4 = 0.
  if(ipres.eq.1) qp4 = 0.
  
  ! Determine which fields need to be calculated
  def_fields = FIELD_PSI + FIELD_PHI + FIELD_ETA
  if(numvar.ge.2) def_fields = def_fields + FIELD_V + FIELD_I
  if(numvar.ge.3) def_fields = def_fields + &
       FIELD_CHI + FIELD_PE + FIELD_B2I + FIELD_J
  if(idens.eq.1) def_fields = def_fields + FIELD_N + FIELD_NI
  if(ipres.eq.1) def_fields = def_fields + FIELD_P

  if(isources.eq.1) then
     def_fields = def_fields + FIELD_SB1
     if(numvar.ge.2) def_fields = def_fields + FIELD_SB2
     if(numvar.ge.3) def_fields = def_fields + FIELD_SP1
  endif

  ! Loop over local elements
  do itri=1,numelms

     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! calculate the field values and derivatives at the sampling points
     call define_fields_79(itri, def_fields)

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

     call ludefvel_n(itri,dbf)
     call ludefphi_n(itri,dbf)
     if(idens.eq.1) call ludefden_n(itri,dbf)
     if(ipres.eq.1) call ludefpres_n(itri,dbf)

  end do
  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumshareddofs so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values
  call sumshareddofs(r4)
  call sumshareddofs(q4)
  if(idens.eq.1) call sumshareddofs(qn4)

  ! Impose boundary conditions
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ! Velocity boundary conditions
!!$  velbounds = 0.
!!$  call boundaryv(iboundv,iboundv2,nbcv)
!!$  if(nbcv .gt. iboundmax) then
!!$     write(*,4881) nbcv, iboundmax
!!$4881 format(" ERROR: nbcv > iboundmax", 2i5)
!!$     call safestop(3)
!!$  endif
!!$  do i=1,nbcv
!!$     call setdiribc(s1matrix_sm, iboundv(i))
!!$  enddo
!!$  call finalizearray4solve(s1matrix_sm)
  call finalizearray4multiply(d1matrix_sm)
  call finalizearray4multiply(r1matrix_sm)
  if(integrator.eq.1) &
       call finalizearray4multiply(o1matrix_sm)
     

  ! Field boundary conditions
!!$  call boundaryp(iboundp,nbcp)
!!$  if(nbcp .gt. iboundmax) then
!!$     write(*,4882) nbcp, iboundmax
!!$4882 format(" ERROR: nbcp > iboundmax", 2i5)
!!$     call safestop(9) 
!!$  endif
!!$  do i=1,nbcp
!!$     call setdiribc(s2matrix_sm,iboundp(i))
!!$  enddo
!!$  call finalizearray4solve(s2matrix_sm)
  call finalizearray4multiply(d2matrix_sm)
  call finalizearray4multiply(r2matrix_sm)
  call finalizearray4multiply(q2matrix_sm)
  if(integrator.eq.1) &
       call finalizearray4multiply(o2matrix_sm)

  ! Density boundary conditions
  if(idens.eq.1) then
     call boundaryds(iboundn,nbcn,1)
     if(nbcn .gt. iboundmax) then
        write(*,4883) nbcn, iboundmax
4883    format(" ERROR: nbcn > iboundmax", 2i5)
        call safestop(9) 
     endif
     
     do i=1,nbcn
        call setdiribc(s8matrix_sm, iboundn(i))
     enddo

     call finalizearray4solve(s8matrix_sm)
     call finalizearray4multiply(d8matrix_sm)
     call finalizearray4multiply(q8matrix_sm)
     call finalizearray4multiply(r8matrix_sm)
     if(integrator.eq.1) &
          call finalizearray4multiply(o8matrix_sm)
  endif ! on idens.eq.1

  ! Pressure boundary conditions
  if(ipres.eq.1) then
     call boundarypres(iboundpres,nbcpres,1)
     if(nbcpres .gt. iboundmax) then
        write(*,4883) nbcpres, iboundmax
4884    format(" ERROR: nbcpres > iboundmax", 2i5)
        call safestop(9) 
     endif
     
     do i=1,nbcpres
        call setdiribc(s9matrix_sm, iboundpres(i))
     enddo

     call finalizearray4solve(s9matrix_sm)
     call finalizearray4multiply(d9matrix_sm)
     call finalizearray4multiply(q9matrix_sm)
     call finalizearray4multiply(r9matrix_sm)
     if(integrator.eq.1) &
          call finalizearray4multiply(o9matrix_sm)
  endif ! on ipres.eq.1

end subroutine ludefall

subroutine ludefvel_n(itri,dbf)

  use basic
  use nintegrate_mod
  use arrays
  use sparse

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif


  implicit none

  integer, intent(in) :: itri
  real, intent(in) :: dbf

  integer :: i, i1, i2, i3, j, j1
  real, dimension(3,3) :: ssterm, ddterm, rrterm, ooterm
  real :: temp
  real :: hypv, thimpv

  thimpv = 1.
  hypv = hyperv*deex**2

  do i=1,18

     i1 = isvaln(itri,i)
     i2 = i1 + 6
     i3 = i2 + 6


     do j=1,18
        j1 = isvaln(itri,j)
        
        ssterm = 0.
        ddterm = 0.
        rrterm = 0.
        ooterm = 0.

        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = v1un(g79(:,:,i),g79(:,:,j),nt79)
        if(integrator.eq.1 .and. ntime.gt.1) then
           ssterm(1,1) = ssterm(1,1) + 1.5*temp
           ddterm(1,1) = ddterm(1,1) + 2.0*temp
           ooterm(1,1) = ooterm(1,1) - 0.5*temp
        else 
           ssterm(1,1) = ssterm(1,1) + temp
           ddterm(1,1) = ddterm(1,1) + temp
        endif

        temp = v1umu(g79(:,:,i),g79(:,:,j))*amu  &              
             + thimp*dt* &
             (v1upsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79))
        if(grav.ne.0) then          
           temp = temp + thimp*dt*&
                v1ungrav(g79(:,:,i),g79(:,:,j),nt79)*grav
        endif
        ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
        ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
        
        temp = v1uun(g79(:,:,i),g79(:,:,j),ph179,nt79) &
             + v1uun(g79(:,:,i),ph179,g79(:,:,j),nt79)
        ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
        ddterm(1,1) = ddterm(1,1) + (.5-thimp)*dt*temp
        
        rrterm(1,1) = rrterm(1,1)  + dt* &
             (v1psipsi(g79(:,:,i),g79(:,:,j),pss79)  &
             +v1psipsi(g79(:,:,i),pss79,g79(:,:,j)))

        if(isources.eq.1) then
           rrterm(1,1) = rrterm(1,1) + thimp*dt*dt* &
             (v1psisb1 (g79(:,:,i),g79(:,:,j),sb179))
        endif
        
        if(linear.eq.1 .or. eqsubtract.eq.1) then
           temp = v1uun(g79(:,:,i),g79(:,:,j),ph079,nt79) &    
                + v1uun(g79(:,:,i),ph079,g79(:,:,j),nt79)      
           ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
           ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
           
           rrterm(1,1) = rrterm(1,1) + thimp*dt*dt* &
                (v1upsipsi(g79(:,:,i),ph079,g79(:,:,j),pss79) &
                +v1upsipsi(g79(:,:,i),ph079,pss79,g79(:,:,j)))
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
           
           temp = v1vvn(g79(:,:,i),g79(:,:,j),vz179,nt79) &
                + v1vvn(g79(:,:,i),vz179,g79(:,:,j),nt79)
           ssterm(1,2) = ssterm(1,2) -     thimp *dt*temp
           ddterm(1,2) = ddterm(1,2) + (.5-thimp)*dt*temp
           
           temp = v2vun(g79(:,:,i),vz179,g79(:,:,j),nt79)
           ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
           ddterm(2,1) = ddterm(2,1) + (.5-thimp)*dt*temp

           temp = v2upsib(g79(:,:,i),g79(:,:,j),pst79,bzt79)
           ssterm(2,1) = ssterm(2,1) - thimp*    thimp *dt*dt*temp
           ddterm(2,1) = ddterm(2,1) + thimp*(1.-thimp)*dt*dt*temp

           temp = v2vn(g79(:,:,i),g79(:,:,j),nt79)
           if(integrator.eq.1 .and. ntime.gt.1) then
              ssterm(2,2) = ssterm(2,2) + 1.5*temp
              ddterm(2,2) = ddterm(2,2) + 2.0*temp
              ooterm(2,2) = ooterm(2,2) - 0.5*temp
           else
              ssterm(2,2) = ssterm(2,2) + temp
              ddterm(2,2) = ddterm(2,2) + temp
           endif

           temp = v2vmu  (g79(:,:,i),g79(:,:,j),amu) &
                + thimp*dt* &
                (v2vpsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79))
           ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
           ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

           if(hypv.gt.0) then
              temp = v2vhypv(g79(:,:,i),g79(:,:,j),amu,hypv)
              ssterm(2,2) = ssterm(2,2) -     thimpv *dt*temp
              ddterm(2,2) = ddterm(2,2) + (1.-thimpv)*dt*temp
           endif

           temp = v2vun(g79(:,:,i),g79(:,:,j),ph179,nt79)
           ssterm(2,2) = ssterm(2,2) -      thimp *dt*temp
           ddterm(2,2) = ddterm(2,2) + (0.5-thimp)*dt*temp

           rrterm(1,2) = rrterm(1,2) + dt* &
                (v1bb     (g79(:,:,i),g79(:,:,j),bzs79)  &
                +v1bb     (g79(:,:,i),bzs79,g79(:,:,j)))
           
           rrterm(2,1) = rrterm(2,1) + dt* &
                (v2psib(g79(:,:,i),g79(:,:,j),bzs79))

           rrterm(2,2) = rrterm(2,2) + dt* &
                (v2psib(g79(:,:,i),pss79,g79(:,:,j)))

           if(isources.eq.1) then 
              rrterm(1,2) = rrterm(2,1) + thimp*dt*dt* &
                (v1bsb2   (g79(:,:,i),g79(:,:,j),sb279))
              rrterm(2,1) = rrterm(2,1) + thimp*dt*dt* &
                   (v2psib(g79(:,:,i),g79(:,:,j),sb279))
              rrterm(2,2) = rrterm(2,2) + thimp*dt*dt* &
                (v2psib(g79(:,:,i),sb179,g79(:,:,j)))
           endif

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              
              temp = v1vvn(g79(:,:,i),g79(:,:,j),vz079,nt79) &
                   + v1vvn(g79(:,:,i),vz079,g79(:,:,j),nt79)
              ssterm(1,2) = ssterm(1,2) -     thimp *dt*temp
              ddterm(1,2) = ddterm(1,2) + (1.-thimp)*dt*temp
              
              temp = v2vun(g79(:,:,i),vz079,g79(:,:,j),nt79)
              ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
              ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp
              
              temp = v2vun(g79(:,:,i),g79(:,:,j),ph079,nt79)
              ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp
              
              rrterm(1,1) = rrterm(1,1) + thimp*dt*dt* &
                   (v1vpsib  (g79(:,:,i),vz079,g79(:,:,j),bzs79))
              
              rrterm(1,2) = rrterm(1,2) + thimp*dt*dt* &
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

           ! regularize the chi equation
           temp = -regular*int2(g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)
           if(integrator.eq.1  .and. ntime.gt.1) then
              ssterm(3,3) = ssterm(3,3) + 1.5*temp
              ddterm(3,3) = ddterm(3,3) + 2.0*temp
              ooterm(3,3) = ooterm(3,3) - 0.5*temp
           else
              ssterm(3,3) = ssterm(3,3) + temp
              ddterm(3,3) = ddterm(3,3) + temp
           endif
           
           temp = v1uchin(g79(:,:,i),g79(:,:,j),ch179,nt79)       
           ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
           ddterm(1,1) = ddterm(1,1) + (.5-thimp)*dt*temp

           temp = v1chin(g79(:,:,i),g79(:,:,j),nt79)
           if(integrator.eq.1 .and. ntime.gt.1) then
              ssterm(1,3) = ssterm(1,3) + 1.5*temp
              ddterm(1,3) = ddterm(1,3) + 2.0*temp
              ooterm(1,3) = ooterm(1,3) - 0.5*temp
           else
              ssterm(1,3) = ssterm(1,3) + temp
              ddterm(1,3) = ddterm(1,3) + temp
           endif

           temp = v1uchin  (g79(:,:,i),ph179,g79(:,:,j),nt79) &
                + v1chichin(g79(:,:,i),g79(:,:,j),ch179,nt79) &
                + v1chichin(g79(:,:,i),ch179,g79(:,:,j),nt79)
           ssterm(1,3) = ssterm(1,3) -     thimp *dt*temp
           ddterm(1,3) = ddterm(1,3) + (.5-thimp)*dt*temp
           
           temp = v1chipsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79)  &
                + v1chibb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)
           if(grav.ne.0) then          
              temp = temp &
                   + v1chingrav(g79(:,:,i),g79(:,:,j),nt79)*grav
           endif
           ssterm(1,3) = ssterm(1,3) - thimp*    thimp *dt*dt*temp
           ddterm(1,3) = ddterm(1,3) + thimp*(1.-thimp)*dt*dt*temp

           temp = v2vchin(g79(:,:,i),g79(:,:,j),ch179,nt79)
           ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
           ddterm(2,2) = ddterm(2,2) + (.5-thimp)*dt*temp
           
           temp = v2vchin(g79(:,:,i),vz179,g79(:,:,j),nt79)
           ssterm(2,3) = ssterm(2,3) -     thimp *dt*temp
           ddterm(2,3) = ddterm(2,3) + (.5-thimp)*dt*temp

           temp = v2chipsib(g79(:,:,i),g79(:,:,j),pst79,bzt79) 
           ssterm(2,3) = ssterm(2,3) - thimp*    thimp *dt*dt*temp
           ddterm(2,3) = ddterm(2,3) + thimp*(1.-thimp)*dt*dt*temp

           temp = v3un(g79(:,:,i),g79(:,:,j),nt79)
           if(integrator.eq.1 .and. ntime.gt.1) then
              ssterm(3,1) = ssterm(3,1) + 1.5*temp
              ddterm(3,1) = ddterm(3,1) + 2.0*temp
              ooterm(3,1) = ooterm(3,1) - 0.5*temp
           else 
              ssterm(3,1) = ssterm(3,1) + temp
              ddterm(3,1) = ddterm(3,1) + temp
           endif

           temp = v3uun  (g79(:,:,i),g79(:,:,j),ph179,nt79) &
               +  v3uun  (g79(:,:,i),ph179,g79(:,:,j),nt79) &
               +  v3uchin(g79(:,:,i),g79(:,:,j),ch179,nt79)
           ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
           ddterm(3,1) = ddterm(3,1) + (.5-thimp)*dt*temp
           
           temp = v3umu   (g79(:,:,i),g79(:,:,j),amu,amuc) &
                +thimp*dt* &
                (v3up     (g79(:,:,i),g79(:,:,j),pt79) &
                +v3upsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79) &
                +v3ubb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)) 
           if(grav.ne.0) then
              temp = temp + thimp*dt* &
                   v3ungrav(g79(:,:,i),g79(:,:,j),nt79)*grav
           endif
           ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
           ddterm(3,1) = ddterm(3,1) + (1.-thimp)*dt*temp
           
           temp = v3vpsib(g79(:,:,i),g79(:,:,j),pst79,bzt79)   
           ssterm(3,2) = ssterm(3,2) - thimp*    thimp *dt*dt*temp
           ddterm(3,2) = ddterm(3,2) + thimp*(1.-thimp)*dt*dt*temp
           
           temp = v3vvn(g79(:,:,i),g79(:,:,j),vz179,nt79) &
                + v3vvn(g79(:,:,i),vz179,g79(:,:,j),nt79)
           ssterm(3,2) = ssterm(3,2) -     thimp *dt*temp
           ddterm(3,2) = ddterm(3,2) + (.5-thimp)*dt*temp

           temp = v3chin(g79(:,:,i),g79(:,:,j),nt79)
           if(integrator.eq.1 .and. ntime.gt.1) then
              ssterm(3,3) = ssterm(3,3) + 1.5*temp
              ddterm(3,3) = ddterm(3,3) + 2.0*temp
              ooterm(3,3) = ooterm(3,3) - 0.5*temp
           else
              ssterm(3,3) = ssterm(3,3) + temp
              ddterm(3,3) = ddterm(3,3) + temp
           endif

           temp = v3chimu    (g79(:,:,i),g79(:,:,j))*2.*amuc   
           ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
           ddterm(3,3) = ddterm(3,3) + (1.-thimp)*dt*temp

           temp = v3uchin  (g79(:,:,i),ph179,g79(:,:,j),nt79) &
                + v3chichin(g79(:,:,i),g79(:,:,j),ch179,nt79) &
                + v3chichin(g79(:,:,i),ch179,g79(:,:,j),nt79)  
           ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
           ddterm(3,3) = ddterm(3,3) + (.5-thimp)*dt*temp
           
           temp = v3chip     (g79(:,:,i),g79(:,:,j),pt79)        &
                + v3chipsipsi(g79(:,:,i),g79(:,:,j),pst79,pst79) &
                + v3chibb    (g79(:,:,i),g79(:,:,j),bzt79,bzt79)
           if(grav.ne.0) then
              temp = temp + &
                   v3chingrav(g79(:,:,i),g79(:,:,j),nt79)*grav
           endif
           ssterm(3,3) = ssterm(3,3) - thimp*    thimp *dt*dt*temp
           ddterm(3,3) = ddterm(3,3) + thimp*(1.-thimp)*dt*dt*temp

           rrterm(3,1) = rrterm(3,1) +                    &
                dt*                                       &
                (v3psipsi(g79(:,:,i),g79(:,:,j),pss79)    & 
                +v3psipsi(g79(:,:,i),pss79,g79(:,:,j)))
           
           rrterm(3,2) = rrterm(3,2) + dt*            &
                (v3bb(g79(:,:,i),g79(:,:,j),bzs79)    &     
                +v3bb(g79(:,:,i),bzs79,g79(:,:,j)))

           rrterm(3,3) = rrterm(3,3) + dt* &
                v3p(g79(:,:,i),g79(:,:,j))

           if(isources.eq.1) then 
              rrterm(3,1) = rrterm(3,1) + & 
                   thimp*dt*dt*                              &
                   (v3psisb1(g79(:,:,i),g79(:,:,j),sb179)    & 
                   +v3psisb1(g79(:,:,i),sb179,g79(:,:,j)))

              rrterm(3,2) = rrterm(3,2) + &
                   thimp*dt*dt*                          &
                   (v3bsb2(g79(:,:,i),g79(:,:,j),sb279))
           end if


           if(gyro.eq.1) then
              temp = g1ub      (g79(:,:,i),g79(:,:,j),      bzt79,pit79,b2i79) &
                   + g1upsipsib(g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(1,1) = ssterm(1,1) +     thimp *dt*temp
              ddterm(1,1) = ddterm(1,1) - (1.-thimp)*dt*temp

              temp = g1vpsi      (g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g1vpsipsipsi(g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g1vpsibb    (g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(1,2) = ssterm(1,2) +     thimp *dt*temp
              ddterm(1,2) = ddterm(1,2) - (1.-thimp)*dt*temp

              temp = g1chib      (g79(:,:,i),g79(:,:,j),      bzt79,pit79,b2i79) &
                   + g1chipsipsib(g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(1,3) = ssterm(1,3) +     thimp *dt*temp
              ddterm(1,3) = ddterm(1,3) - (1.-thimp)*dt*temp

              temp = g2upsi      (g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g2upsipsipsi(g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g2upsibb    (g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(2,1) = ssterm(2,1) +     thimp *dt*temp
              ddterm(2,1) = ddterm(2,1) - (1.-thimp)*dt*temp

              temp = g2vb(g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(2,2) = ssterm(2,2) +     thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) - (1.-thimp)*dt*temp

              temp = g2chipsi      (g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g2chipsipsipsi(g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g2chipsibb    (g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(2,3) = ssterm(2,3) +     thimp *dt*temp
              ddterm(2,3) = ddterm(2,3) - (1.-thimp)*dt*temp

              temp = g3ub      (g79(:,:,i),g79(:,:,j),bzt79,      pit79,b2i79) &
                   + g3upsipsib(g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(3,1) = ssterm(3,1) +     thimp *dt*temp
              ddterm(3,1) = ddterm(3,1) - (1.-thimp)*dt*temp

              temp = g3vpsi      (g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g3vpsipsipsi(g79(:,:,i),g79(:,:,j),pst79,      pit79,b2i79) &
                   + g3vpsibb    (g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(3,2) = ssterm(3,2) +     thimp *dt*temp
              ddterm(3,2) = ddterm(3,2) - (1.-thimp)*dt*temp

              temp = g3chib      (g79(:,:,i),g79(:,:,j),      bzt79,pit79,b2i79) &
                   + g3chipsipsib(g79(:,:,i),g79(:,:,j),pst79,bzt79,pit79,b2i79)
              ssterm(3,3) = ssterm(3,3) +     thimp *dt*temp
              ddterm(3,3) = ddterm(3,3) - (1.-thimp)*dt*temp
           endif
           
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
              
              temp = v3vvn(g79(:,:,i),g79(:,:,j),vz079,nt79) &
                   + v3vvn(g79(:,:,i),vz079,g79(:,:,j),nt79)
              ssterm(3,2) = ssterm(3,2) -     thimp *dt*temp
              ddterm(3,2) = ddterm(3,2) + (1.-thimp)*dt*temp
              
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
        
        call insertval(s1matrix_sm,ssterm(1,1),i1,j1,1)
        call insertval(d1matrix_sm,ddterm(1,1),i1,j1,1)
        call insertval(r1matrix_sm,rrterm(1,1),i1,j1,1)
        if(integrator.eq.1) then
           call insertval(o1matrix_sm,ooterm(1,1),i1,j1,1)
        endif
        if(numvar.ge.2) then
           call insertval(s1matrix_sm,ssterm(1,2),i1  ,j1+6,1)
           call insertval(s1matrix_sm,ssterm(2,1),i1+6,j1  ,1)
           call insertval(s1matrix_sm,ssterm(2,2),i1+6,j1+6,1)
           call insertval(d1matrix_sm,ddterm(1,2),i1  ,j1+6,1)
           call insertval(d1matrix_sm,ddterm(2,1),i1+6,j1  ,1)
           call insertval(d1matrix_sm,ddterm(2,2),i1+6,j1+6,1)
           call insertval(r1matrix_sm,rrterm(1,2),i1  ,j1+6,1)
           call insertval(r1matrix_sm,rrterm(2,1),i1+6,j1  ,1)
           call insertval(r1matrix_sm,rrterm(2,2),i1+6,j1+6,1)
           if(integrator.eq.1) then
              call insertval(o1matrix_sm,ooterm(1,2),i1  ,j1+6,1)
              call insertval(o1matrix_sm,ooterm(2,1),i1+6,j1  ,1)
              call insertval(o1matrix_sm,ooterm(2,2),i1+6,j1+6,1)
           endif
        endif
        if(numvar.ge.3) then
           call insertval(s1matrix_sm,ssterm(1,3),i1,   j1+12,1)
           call insertval(s1matrix_sm,ssterm(2,3),i1+6, j1+12,1)
           call insertval(s1matrix_sm,ssterm(3,3),i1+12,j1+12,1)
           call insertval(s1matrix_sm,ssterm(3,1),i1+12,j1,   1)
           call insertval(s1matrix_sm,ssterm(3,2),i1+12,j1+6, 1)
           call insertval(d1matrix_sm,ddterm(1,3),i1,   j1+12,1)
           call insertval(d1matrix_sm,ddterm(2,3),i1+6, j1+12,1)
           call insertval(d1matrix_sm,ddterm(3,3),i1+12,j1+12,1)
           call insertval(d1matrix_sm,ddterm(3,1),i1+12,j1,   1)
           call insertval(d1matrix_sm,ddterm(3,2),i1+12,j1+6, 1)
           call insertval(r1matrix_sm,rrterm(1,3),i1,   j1+12,1)
           call insertval(r1matrix_sm,rrterm(2,3),i1+6, j1+12,1)
           call insertval(r1matrix_sm,rrterm(3,3),i1+12,j1+12,1)
           call insertval(r1matrix_sm,rrterm(3,1),i1+12,j1,   1)
           call insertval(r1matrix_sm,rrterm(3,2),i1+12,j1+6, 1)
           if(integrator.eq.1) then
              call insertval(o1matrix_sm,ooterm(1,3),i1,   j1+12,1)
              call insertval(o1matrix_sm,ooterm(2,3),i1+6, j1+12,1)
              call insertval(o1matrix_sm,ooterm(3,3),i1+12,j1+12,1)
              call insertval(o1matrix_sm,ooterm(3,1),i1+12,j1,   1)
              call insertval(o1matrix_sm,ooterm(3,2),i1+12,j1+6, 1)
           endif
        endif
     enddo               ! on j

     ! Definition of R4
     ! ================
     if(numvar.ge.3) then
        r4(i3) = r4(i3) + thimp*dt*dt*v3p(g79(:,:,i),sp179)
     endif

     if(grav.ne.0) then          
        r4(i1) = r4(i1) + dt* &
             v1ngrav(g79(:,:,i),nt79)*grav
        
        if(numvar.ge.3) then
           r4(i3) = r4(i3) + dt* &
                v3ngrav(g79(:,:,i),nt79)*grav
        endif
     endif
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        r4(i1) = r4(i1) + dt* &
             (v1umu    (g79(:,:,i),ph079))
        if(isources.eq.1) then
           r4(i1) = r4(i1) + thimp*dt*dt* &
                (v1psisb1 (g79(:,:,i),ps079,sb179))               ! passed: 1, 2
        endif
        
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
                (v2vmu  (g79(:,:,i),vz079,amu) &
                +v2vhypv(g79(:,:,i),vz079,amu,hypv))
           
           ! DENSITY TERMS
           if(grav.ne.0) then          
              r4(i1) = r4(i1) + thimp*dt*dt*grav* &
                   (v1ungrav   (g79(:,:,i),ph079,nt79) &
                   +v1chingrav (g79(:,:,i),ch079,nt79) &
                   )
!!$                      +v1ndenmgrav(g79(:,:,i),nt79))
              
           endif
           r4(i1) = r4(i1) + dt* &
                (v1vvn(g79(:,:,i),vz079,vz079,nt79))
           r4(i2) = r4(i2) + dt* &
                (v2vun(g79(:,:,i),vz079,ph079,nt79))
           if(isources.eq.1) then
              r4(i2) = r4(i2) + thimp*dt*dt* &
                   (v2psib(g79(:,:,i),sb179,bz079) &
                   +v2psib(g79(:,:,i),ps079,sb279))
           endif
           
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
           
           r4(i3) = r4(i3) + dt* &
                (v3umu   (g79(:,:,i),ph079,amu,amuc) &
                +v3chimu (g79(:,:,i),ch079)*amu)
           if(isources.eq.1) then
              r4(i3) = r4(i3) + thimp*dt*dt* &
                (v3psisb1(g79(:,:,i),ps079,sb179) &
                +v3psisb1(g79(:,:,i),sb179,ps079) &
                +v3bsb2  (g79(:,:,i),bz079,sb279))
           endif
                
           
           ! DENSITY TERMS
           r4(i1) = r4(i1) + dt* &
                (v1uchin  (g79(:,:,i),ph079,ch079,nt79) &
                +v1chichin(g79(:,:,i),ch079,ch079,nt79))
           r4(i3) = r4(i3) + dt* &
                (v3uun    (g79(:,:,i),ph079,ph079,nt79) &
                +v3uchin  (g79(:,:,i),ph079,ch079,nt79) &
                +v3vvn    (g79(:,:,i),vz079,vz079,nt79) &
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
end subroutine ludefvel_n


subroutine ludefphi_n(itri,dbf)
  use basic
  use nintegrate_mod
  use arrays
  use sparse

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif


  implicit none

  integer, intent(in) :: itri
  real, intent(in) :: dbf

  integer :: i, i1, i2, i3, j, j1
  real, dimension(3,3) :: ssterm, ddterm, rrterm, qqterm, ooterm
  real :: temp, hypf, hypi, hypp, hypv, hypc

  hypf = hyper *deex**2
  hypi = hyperi*deex**2
  hypp = hyperp*deex**2
  hypv = hyperv*deex**2
  hypc = hyperc*deex**2

  do i=1,18

     i1 = isvaln(itri,i)
     i2 = i1 + 6
     i3 = i2 + 6    

     do j=1,18
        
        ssterm = 0.
        ddterm = 0.
        rrterm = 0.
        qqterm = 0.
        ooterm = 0.

        j1 = isvaln(itri,j)

        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = b1psi(g79(:,:,i),g79(:,:,j))
        if(integrator.eq.1 .and. ntime.gt.1) then
           ssterm(1,1) = ssterm(1,1) + 1.5*temp
           ddterm(1,1) = ddterm(1,1) + 2.0*temp
           ooterm(1,1) = ooterm(1,1) - 0.5*temp
        else
           ssterm(1,1) = ssterm(1,1) + temp
           ddterm(1,1) = ddterm(1,1) + temp
        endif

        temp = b1psieta(g79(:,:,i),g79(:,:,j),eta79,hypf)  &
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
           if(integrator.eq.1 .and. ntime.gt.1) then
              ssterm(2,2) = ssterm(2,2) + 1.5*temp
              ddterm(2,2) = ddterm(2,2) + 2.0*temp
              ooterm(2,2) = ooterm(2,2) - 0.5*temp
           else
              ssterm(2,2) = ssterm(2,2) + temp
              ddterm(2,2) = ddterm(2,2) + temp
           endif

           temp = b2beta(g79(:,:,i),g79(:,:,j),eta79,hypi) &
                + b2bu  (g79(:,:,i),g79(:,:,j),pht79)
           ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
           ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

           temp = b2bbd(g79(:,:,i),g79(:,:,j),bz179,ni79)*dbf &
                + b2bbd(g79(:,:,i),bz179,g79(:,:,j),ni79)*dbf
           ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
           ddterm(2,2) = ddterm(2,2) + (.5-thimp)*dt*temp

           temp = b2bu  (g79(:,:,i),bz179,g79(:,:,j))
           rrterm(2,1) = rrterm(2,1) + thimp*dt*temp
           qqterm(2,1) = qqterm(2,1) - thimp*dt*temp

           temp = b2psiv(g79(:,:,i),ps179,g79(:,:,j))
           rrterm(2,2) = rrterm(2,2) + thimp*dt*temp
           qqterm(2,2) = qqterm(2,2) - thimp*dt*temp

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = b1psibd(g79(:,:,i),g79(:,:,j),bz079,ni79)*dbf              
              ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
              ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
              
              temp = b1psibd(g79(:,:,i),ps079,g79(:,:,j),ni79)*dbf              
              ssterm(1,2) = ssterm(1,2) -     thimp *dt*temp
              ddterm(1,2) = ddterm(1,2) + (1.-thimp)*dt*temp

              temp = b2psipsid(g79(:,:,i),g79(:,:,j),ps079,ni79)*dbf &
                   + b2psipsid(g79(:,:,i),ps079,g79(:,:,j),ni79)*dbf
              ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
              ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp

              temp = b2bbd (g79(:,:,i),g79(:,:,j),bz079,ni79)*dbf &
                   + b2bbd (g79(:,:,i),bz079,g79(:,:,j),ni79)*dbf
              ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
              ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

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

           temp = b1psichi(g79(:,:,i),g79(:,:,j),cht79)                
           ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
           ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

           temp = b1psichi(g79(:,:,i),ps179,g79(:,:,j))                
           rrterm(1,3) = rrterm(1,3) + thimp*dt*temp
           qqterm(1,3) = qqterm(1,3) - thimp*dt*temp

           temp = b2bchi(g79(:,:,i),g79(:,:,j),cht79)                  
           ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
           ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

           temp = b2ped(g79(:,:,i),g79(:,:,j),ni79)*dbf*pefac
           ssterm(2,3) = ssterm(2,3) -     thimp *dt*temp
           ddterm(2,3) = ddterm(2,3) + (1.-thimp)*dt*temp

           temp = b2bchi(g79(:,:,i),bz179,g79(:,:,j))                  
           rrterm(2,3) = rrterm(2,3) + thimp*dt*temp
           qqterm(2,3) = qqterm(2,3) - thimp*dt*temp

           temp = b3psipsieta(g79(:,:,i),g79(:,:,j),ps179,eta79) &
                + b3psipsieta(g79(:,:,i),ps179,g79(:,:,j),eta79)
           ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
           ddterm(3,1) = ddterm(3,1) + (.5-thimp)*dt*temp

           temp = b3pebd(g79(:,:,i),pe179,g79(:,:,j),ni79)*dbf*pefac
           ssterm(3,2) = ssterm(3,2) - thimp*dt*temp
           ddterm(3,2) = ddterm(3,2) - thimp*dt*temp

           temp = b3bbeta(g79(:,:,i),g79(:,:,j),bz179,eta79) &
                + b3bbeta(g79(:,:,i),bz179,g79(:,:,j),eta79)
           ssterm(3,2) = ssterm(3,2) -     thimp *dt*temp
           ddterm(3,2) = ddterm(3,2) + (.5-thimp)*dt*temp

           temp = b3pe(g79(:,:,i),g79(:,:,j))
           if(integrator.eq.1 .and. ntime.gt.1) then
              ssterm(3,3) = ssterm(3,3) + 1.5*temp
              ddterm(3,3) = ddterm(3,3) + 2.0*temp
              ooterm(3,3) = ooterm(3,3) - 0.5*temp
           else
              ssterm(3,3) = ssterm(3,3) + temp
              ddterm(3,3) = ddterm(3,3) + temp
           endif

           temp = b3pebd(g79(:,:,i),g79(:,:,j),bzt79,ni79)*dbf*pefac &
                + p1pu  (g79(:,:,i),g79(:,:,j),pht79)                & 
                + p1pchi(g79(:,:,i),g79(:,:,j),cht79)                & 
                + b3pedkappa(g79(:,:,i),g79(:,:,j),ni79,kappat,hypp)*(gam-1.)
           ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
           ddterm(3,3) = ddterm(3,3) + (1.-thimp)*dt*temp

           temp = p1pu(g79(:,:,i),pe179,g79(:,:,j))
           rrterm(3,1) = rrterm(3,1) + thimp*dt*temp
           qqterm(3,1) = qqterm(3,1) - thimp*dt*temp

           temp = p1pchi(g79(:,:,i),pe179,g79(:,:,j))
           rrterm(3,3) = rrterm(3,3) + thimp*dt*temp
           qqterm(3,3) = qqterm(3,3) - thimp*dt*temp

           ! Anisotropic Heat Flux
           if(kappar.ne.0) then
              temp = kappar*(gam-1.)* &
                   (p1kappar(g79(:,:,i),g79(:,:,j),pst79,pet79,ni79,b2i79) &
                   +p1kappar(g79(:,:,i),pst79,g79(:,:,j),pet79,ni79,b2i79))
              ssterm(3,1) = ssterm(3,1) - thimp*dt*temp
              ddterm(3,1) = ddterm(3,1) - thimp*dt*temp

              temp = kappar*(gam-1.)* &
                   p1kappar(g79(:,:,i),pst79,pst79,g79(:,:,j),ni79,b2i79)
              ssterm(3,3) = ssterm(3,3) -     thimp *dt*temp
              ddterm(3,3) = ddterm(3,3) + (1.-thimp)*dt*temp
           endif 
              
           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = b1psichi(g79(:,:,i),ps079,g79(:,:,j))
              rrterm(1,3) = rrterm(1,3) +     thimp *dt*temp
              qqterm(1,3) = qqterm(1,3) + (1.-thimp)*dt*temp

              temp = b2bchi(g79(:,:,i),bz079,g79(:,:,j))
              rrterm(2,3) = rrterm(2,3) +     thimp *dt*temp
              qqterm(2,3) = qqterm(2,3) + (1.-thimp)*dt*temp

              temp = b3psipsieta(g79(:,:,i),g79(:,:,j),ps079,eta79) &
                   + b3psipsieta(g79(:,:,i),ps079,g79(:,:,j),eta79)
              ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp
              ddterm(3,1) = ddterm(3,1) + (1.-thimp)*dt*temp

              temp = b3bbeta(g79(:,:,i),g79(:,:,j),bz079,eta79) &
                   + b3bbeta(g79(:,:,i),bz079,g79(:,:,j),eta79)
              ssterm(3,2) = ssterm(3,2) -     thimp *dt*temp
              ddterm(3,2) = ddterm(3,2) + (1.-thimp)*dt*temp

              temp = b3pebd(g79(:,:,i),pe079,g79(:,:,j),ni79)*dbf*pefac
              ssterm(3,2) = ssterm(3,2) -     thimp *dt*temp
              ddterm(3,2) = ddterm(3,2) + (1.-thimp)*dt*temp

              temp = p1pu(g79(:,:,i),pe079,g79(:,:,j))
              rrterm(3,1) = rrterm(3,1) +     thimp *dt*temp
              qqterm(3,1) = qqterm(3,1) + (1.-thimp)*dt*temp

              temp = p1pchi(g79(:,:,i),pe079,g79(:,:,j))                
              rrterm(3,3) = rrterm(3,3) +     thimp *dt*temp
              qqterm(3,3) = qqterm(3,3) + (1.-thimp)*dt*temp

              ! Anisotropic Heat Flux
              if(kappar.ne.0) then
                 temp = kappar*(gam-1.)* &
                      (p1kappar(g79(:,:,i),g79(:,:,j),pss79,pe079,ni79,b2i79) &
                      +p1kappar(g79(:,:,i),pss79,g79(:,:,j),pe079,ni79,b2i79))
                 ddterm(3,1) = ddterm(3,1) + dt*temp
              endif
           endif
        endif

        call insertval(s2matrix_sm,ssterm(1,1),i1,j1,1)
        call insertval(d2matrix_sm,ddterm(1,1),i1,j1,1)
        call insertval(r2matrix_sm,rrterm(1,1),i1,j1,1)
        call insertval(q2matrix_sm,qqterm(1,1),i1,j1,1)
        if(integrator.eq.1) then
           call insertval(o2matrix_sm,ooterm(1,1),i1,j1,1)
        endif
        if(numvar.ge.2) then
           call insertval(s2matrix_sm,ssterm(1,2),i1  ,j1+6,1)
           call insertval(s2matrix_sm,ssterm(2,1),i1+6,j1  ,1)
           call insertval(s2matrix_sm,ssterm(2,2),i1+6,j1+6,1)
           call insertval(d2matrix_sm,ddterm(1,2),i1  ,j1+6,1)
           call insertval(d2matrix_sm,ddterm(2,1),i1+6,j1  ,1)
           call insertval(d2matrix_sm,ddterm(2,2),i1+6,j1+6,1)
           call insertval(r2matrix_sm,rrterm(1,2),i1  ,j1+6,1)
           call insertval(r2matrix_sm,rrterm(2,1),i1+6,j1  ,1)
           call insertval(r2matrix_sm,rrterm(2,2),i1+6,j1+6,1)
           call insertval(q2matrix_sm,qqterm(1,2),i1  ,j1+6,1)
           call insertval(q2matrix_sm,qqterm(2,1),i1+6,j1  ,1)
           call insertval(q2matrix_sm,qqterm(2,2),i1+6,j1+6,1)
           if(integrator.eq.1) then
              call insertval(o2matrix_sm,ooterm(1,2),i1  ,j1+6,1)
              call insertval(o2matrix_sm,ooterm(2,1),i1+6,j1  ,1)
              call insertval(o2matrix_sm,ooterm(2,2),i1+6,j1+6,1)
           endif
        endif
        if(numvar .eq. 3) then
           call insertval(s2matrix_sm,ssterm(1,3),i1,   j1+12,1)
           call insertval(s2matrix_sm,ssterm(2,3),i1+6, j1+12,1)
           call insertval(s2matrix_sm,ssterm(3,3),i1+12,j1+12,1)
           call insertval(s2matrix_sm,ssterm(3,1),i1+12,j1,   1)
           call insertval(s2matrix_sm,ssterm(3,2),i1+12,j1+6, 1)
           call insertval(d2matrix_sm,ddterm(1,3),i1,   j1+12,1)
           call insertval(d2matrix_sm,ddterm(2,3),i1+6, j1+12,1)
           call insertval(d2matrix_sm,ddterm(3,3),i1+12,j1+12,1)
           call insertval(d2matrix_sm,ddterm(3,1),i1+12,j1,   1)
           call insertval(d2matrix_sm,ddterm(3,2),i1+12,j1+6, 1)
           call insertval(r2matrix_sm,rrterm(1,3),i1,   j1+12,1)
           call insertval(r2matrix_sm,rrterm(2,3),i1+6, j1+12,1)
           call insertval(r2matrix_sm,rrterm(3,3),i1+12,j1+12,1)
           call insertval(r2matrix_sm,rrterm(3,1),i1+12,j1,   1)
           call insertval(r2matrix_sm,rrterm(3,2),i1+12,j1+6, 1)
           call insertval(q2matrix_sm,qqterm(1,3),i1,   j1+12,1)
           call insertval(q2matrix_sm,qqterm(2,3),i1+6, j1+12,1)
           call insertval(q2matrix_sm,qqterm(3,3),i1+12,j1+12,1)
           call insertval(q2matrix_sm,qqterm(3,1),i1+12,j1,   1)
           call insertval(q2matrix_sm,qqterm(3,2),i1+12,j1+6, 1)
           if(integrator.eq.1) then
              call insertval(o2matrix_sm,ooterm(1,3),i1,   j1+12,1)
              call insertval(o2matrix_sm,ooterm(2,3),i1+6, j1+12,1)
              call insertval(o2matrix_sm,ooterm(3,3),i1+12,j1+12,1)
              call insertval(o2matrix_sm,ooterm(3,1),i1+12,j1,   1)
              call insertval(o2matrix_sm,ooterm(3,2),i1+12,j1+6, 1)
           endif
        endif
     enddo ! on j

     if(numvar.ge.3) then
        ! ohmic heating
        q4(i3) = q4(i3) + dt*(gam-1.)* &
             (qpsipsieta(g79(:,:,i),pst79,pst79,eta79,hypf,jt79) &
             +qbbeta    (g79(:,:,i),bzt79,bzt79,eta79,hypi))

        ! viscous heating
        if(ipres.eq.0) then
           q4(i3) = q4(i3) - dt*(gam-1.)* &
                (quumu    (g79(:,:,i),pht79,pht79,amu,amuc,hypc) &
                +qvvmu    (g79(:,:,i),vzt79,vzt79,amu,     hypv) &
                +quchimu  (g79(:,:,i),pht79,cht79,amu,amuc,hypc) &
                +qchichimu(g79(:,:,i),cht79,cht79,amu,amuc,hypc))
        endif
     endif
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        
        q4(i1) = q4(i1) + dt* &
             (b1psieta(g79(:,:,i),ps079,eta79,hypf))

        ! EQUILIBRIUM TERMS
        q4(i1) = q4(i1) + dt* &
             (b1psiu(g79(:,:,i),ps079,ph079))

        if(numvar.ge.2) then
           q4(i2) = q4(i2) + dt* &
                (b2beta     (g79(:,:,i),bz079,eta79,hypi))
              
           ! DENSITY TERMS
           q4(i1) = q4(i1) + dt* &
                (b1psibd  (g79(:,:,i),ps079,bz079,ni79)*dbf)
           q4(i2) = q4(i2) + dt* &
                (b2psipsid(g79(:,:,i),ps079,ps079,ni79)*dbf &
                +b2bbd    (g79(:,:,i),bz079,bz079,ni79)*dbf)

           ! EQUILIBRIUM TERMS
           q4(i2) = q4(i2) + dt* &
                (b2psiv(g79(:,:,i),ps079,vz079) &
                +b2bu  (g79(:,:,i),bz079,ph079))
        endif

        if(numvar.ge.3) then
           q4(i3) = q4(i3) + dt* &
                (b3psipsieta(g79(:,:,i),ps079,ps079,eta79) &
                +b3bbeta    (g79(:,:,i),bz079,bz079,eta79) &
                +b3pedkappa (g79(:,:,i),pe079,ni79,kappat,hypp)*(gam-1.))

           ! DENSITY TERMS
           q4(i2) = q4(i2) + dt* &
                (b2ped (g79(:,:,i),pe079,ni79)*dbf*pefac)
           q4(i3) = q4(i3) + dt* &
                (b3pebd(g79(:,:,i),pe079,bz079,ni79)*dbf*pefac &
                +p1kappar(g79(:,:,i),ps079,ps079,pe079,ni79,b2i79)*kappar*(gam-1.))
              
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

end subroutine ludefphi_n


subroutine ludefden_n(itri,dbf)

  use basic
  use nintegrate_mod
  use arrays
  use sparse

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif


  implicit none

  integer, intent(in) :: itri
  real, intent(in) :: dbf

  integer :: i, i1, ione, j, j1, jone
  real :: ssterm, ddterm, ooterm
  real, dimension(3) :: rrterm, qqterm
  real :: temp, hypp

  hypp = hyperp*deex**2

  do i=1,18
     ione = isval1(itri,i)
     i1 = isvaln(itri,i)
     
     do j=1,18         
        ssterm = 0.
        ddterm = 0.
        rrterm = 0.
        qqterm = 0.
        ooterm = 0.

        jone = isval1(itri,j)
        j1 = isvaln(itri,j)

        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = n1n(g79(:,:,i),g79(:,:,j))
        if(integrator.eq.1 .and. ntime.gt.1) then
           ssterm = ssterm + 1.5*temp
           ddterm = ddterm + 2.0*temp
           ooterm = ooterm - 0.5*temp
        else
           ssterm = ssterm + temp
           ddterm = ddterm + temp
        endif
        
        temp = n1ndenm(g79(:,:,i),g79(:,:,j),dt*denm,hypp) &
             + n1nu   (g79(:,:,i),g79(:,:,j),pht79)
        ssterm = ssterm -     thimp *dt*temp
        ddterm = ddterm + (1.-thimp)*dt*temp

        temp = n1nu(g79(:,:,i),n179,g79(:,:,j))
        rrterm(1) = rrterm(1) + thimp*dt*temp
        qqterm(1) = qqterm(1) - thimp*dt*temp

        if(linear.eq.1 .or. eqsubtract.eq.1) then
           temp = n1nu  (g79(:,:,i),n079,g79(:,:,j))
           rrterm(1) = rrterm(1) +     thimp *dt*temp
           qqterm(1) = qqterm(1) + (1.-thimp)*dt*temp
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
           endif
        endif

        call insertval(s8matrix_sm, ssterm, ione, jone, 1)
        call insertval(d8matrix_sm, ddterm, ione, jone, 1)
        call insertval(r8matrix_sm, rrterm(1), i1, j1, 1)
        call insertval(q8matrix_sm, qqterm(1), i1, j1, 1)
        if(integrator.eq.1) then
           call insertval(o8matrix_sm, ooterm, ione, jone, 1)
        endif
        if(numvar.ge.2) then
           call insertval(r8matrix_sm,rrterm(2), i1,j1+6,1)
           call insertval(q8matrix_sm,qqterm(2), i1,j1+6,1)
        endif
        if(numvar.ge.3) then
           call insertval(r8matrix_sm,rrterm(3), i1,j1+12,1)
           call insertval(q8matrix_sm,qqterm(3), i1,j1+12,1)
        endif

     enddo                     ! on j

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        
        qn4(ione) = qn4(ione) + dt* &
             (n1ndenm(g79(:,:,i),n079,dt*denm,hypp))
        
        ! EQUILIBRIUM TERMS
        qn4(ione) = qn4(ione) + dt* &
             (n1nu   (g79(:,:,i),n079,ph079))
        
        if(numvar.ge.3) then
                               
           ! EQUILIBRIUM TERMS
           qn4(ione) = qn4(ione) + dt* &
                (n1nchi(g79(:,:,i),n079,ch079))
        endif
     endif

  enddo                     ! on i
end subroutine ludefden_n


subroutine ludefpres_n(itri,dbf)

  use basic
  use nintegrate_mod
  use arrays
  use sparse

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif


  implicit none

  integer, intent(in) :: itri
  real, intent(in) :: dbf

  integer :: i, i1, ione, j, j1, jone
  real :: ssterm, ddterm, ooterm
  real, dimension(3) :: rrterm, qqterm
  real :: temp, hypp

  hypp = hyperp*deex**2

  do i=1,18
     ione = isval1(itri,i)
     i1 = isvaln(itri,i)
     
     do j=1,18         
        ssterm = 0.
        ddterm = 0.
        rrterm = 0.
        qqterm = 0.
        ooterm = 0.

        jone = isval1(itri,j)
        j1 = isvaln(itri,j)


        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = int2(g79(:,:,i),g79(:,:,j),weight_79,79)
        if(integrator.eq.1 .and. ntime.gt.1) then
           ssterm = ssterm + 1.5*temp
           ddterm = ddterm + 2.0*temp
           ooterm = ooterm - 0.5*temp
        else 
           ssterm = ssterm + temp
           ddterm = ddterm + temp
        endif
        
        temp = p1pu   (g79(:,:,i),g79(:,:,j),pht79)
        ssterm = ssterm -     thimp *dt*temp
        ddterm = ddterm + (1.-thimp)*dt*temp

        temp = p1pu(g79(:,:,i),p179,g79(:,:,j))
        rrterm(1) = rrterm(1) + thimp*dt*temp
        qqterm(1) = qqterm(1) - thimp*dt*temp

        if(linear.eq.1 .or. eqsubtract.eq.1) then
           temp = p1pu  (g79(:,:,i),p079,g79(:,:,j))
           rrterm(1) = rrterm(1) +     thimp *dt*temp
           qqterm(1) = qqterm(1) + (1.-thimp)*dt*temp
        endif

        ! NUMVAR = 3
        ! ~~~~~~~~~~
        if(numvar.ge.3) then
           temp = p1pchi    (g79(:,:,i),g79(:,:,j),cht79) &
                + b3pedkappa(g79(:,:,i),g79(:,:,j),ni79,kappat,hypp)*(gam-1.) &
                + p1kappar  (g79(:,:,i),ps079,ps079,g79(:,:,j),ni79,b2i79)*kappar*(gam-1.)
           ssterm = ssterm -     thimp *dt*temp
           ddterm = ddterm + (1.-thimp)*dt*temp
           
           temp = p1pchi(g79(:,:,i),p179,g79(:,:,j))
           rrterm(3) = rrterm(3) + thimp*dt*temp
           qqterm(3) = qqterm(3) - thimp*dt*temp

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = p1pchi(g79(:,:,i),p079,g79(:,:,j))
              rrterm(3) = rrterm(3) +     thimp *dt*temp
              qqterm(3) = qqterm(3) + (1.-thimp)*dt*temp
           endif
        endif
        
        call insertval(s9matrix_sm, ssterm, ione, jone, 1)
        call insertval(d9matrix_sm, ddterm, ione, jone, 1)
        call insertval(r9matrix_sm, rrterm(1), i1, j1, 1)
        call insertval(q9matrix_sm, qqterm(1), i1, j1, 1)
        if(integrator.eq.1) then
           call insertval(o9matrix_sm, ooterm, ione, jone, 1)
        endif
        if(numvar.ge.2) then
           call insertval(r9matrix_sm,rrterm(2), i1,j1+6,1)
           call insertval(q9matrix_sm,qqterm(2), i1,j1+6,1)
        endif
        if(numvar.ge.3) then
           call insertval(r9matrix_sm,rrterm(3), i1,j1+12,1)
           call insertval(q9matrix_sm,qqterm(3), i1,j1+12,1)
        endif

     enddo                     ! on j

     qp4(ione) = qp4(ione) + dt* &
          (b3psipsieta(g79(:,:,i),pst79,pst79,eta79))

     if(numvar.ge.2) then
        qp4(ione) = qp4(ione) + dt* &
             (b3bbeta(g79(:,:,i),bzt79,bzt79,eta79))
     endif

     if(numvar.ge.3) then
        qp4(ione) = qp4(ione) + dt* &
             (b3pebd(g79(:,:,i),pet79,bzt79,ni79))
     endif

     if(linear.eq.1 .or. eqsubtract.eq.1) then
                
        qp4(ione) = qp4(ione) + dt* &
             (b3pedkappa(g79(:,:,i),p079,ni79,kappat,hypp)*(gam-1.) &
             +p1kappar  (g79(:,:,i),ps079,ps079,p079,ni79,b2i79)*kappar*(gam-1.))

        ! EQUILIBRIUM TERMS
        qp4(ione) = qp4(ione) + dt* &
             (p1pu   (g79(:,:,i),p079,ph079))
       
        if(numvar.ge.3) then
           ! EQUILIBRIUM TERMS
           qp4(ione) = qp4(ione) + dt* &
                (p1pchi(g79(:,:,i),p079,ch079))
        endif
     endif

  enddo                     ! on i
end subroutine ludefpres_n

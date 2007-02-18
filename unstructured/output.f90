module output_files
  integer, parameter :: FILE__C1NEW = 9
  integer, parameter :: FILE__PSI_PER = 11
  integer, parameter :: FILE__PSI_FULL = 12
  integer, parameter :: FILE__J_PER = 13
  integer, parameter :: FILE__J_FULL = 14
  integer, parameter :: FILE__I_PER = 15
  integer, parameter :: FILE__I_FULL = 16
  integer, parameter :: FILE__PE_PER = 17
  integer, parameter :: FILE__PE_FULL = 18
  integer, parameter :: FILE__PHI_FULL = 19
  integer, parameter :: FILE__VOR_FULL = 20
  integer, parameter :: FILE__V_FULL = 21
  integer, parameter :: FILE__CHI_FULL = 22
  integer, parameter :: FILE__DIV_FULL = 23
  integer, parameter :: FILE__DENSITY = 24
  integer, parameter :: FILE__MATRIX = 31
  integer, parameter :: FILE__ERROR = 65

end module output_files

!============================================================
subroutine openf

  ! open output files and initialize NCAR graphics
  use output_files
  use p_data
  use t_data
  use basic
  use arrays
  implicit none
  integer :: jj, maxhdf
  real :: alx, alz

  call getboundingboxsize(alx, alz)

  maxhdf = 0

  if (myrank.eq.0) then !Serialize I/O
     call ncarcgm(1,'C1new.cgm')
     call dders(-1)

     ! define evenly spaced coordinates and create new HDF5 file
     xary = (/ ((alx*jj)/(ires-1), jj=0,ires-1) /)
     yary = (/ ((alz*jj)/(ires-1), jj=0,ires-1) /)
     if(numvar.eq.1) then
        if(linear.eq.1) maxhdf=4
        if(linear.eq.0) maxhdf=6
     else if(numvar.eq.2) then
        if(linear.eq.1) maxhdf=6
        if(linear.eq.0) maxhdf=9
     else if(numvar.eq.3) then
        if(linear.eq.1) maxhdf=9
        if(linear.eq.0) maxhdf=13
     endif
         
!!$     call createHDF5(xary,ires,yary,ires,maxhdf)

     ! initialize minimum and maximum
     maf = -1.e20
     mif =  1.e20

     ! open ascii output files
     open(FILE__C1NEW, file='C1new.out',form='formatted',status='unknown')
     open(FILE__ERROR, file='C1error.out',form='formatted',status='unknown')
     open(FILE__PSI_PER,file='C1psi-per',form='formatted',status='unknown')
     open(FILE__J_PER,file='C1J-per',form='formatted',status='unknown')
     open(FILE__PHI_FULL,file='C1phi-full',form='formatted',status='unknown')
     open(FILE__VOR_FULL,file='C1vor-full',form='formatted',status='unknown')
     if(numvar.ge.2) then
        open(FILE__I_PER,file='C1I-per',form='formatted',status='unknown')
        open(FILE__V_FULL,file='C1v-full',form='formatted',status='unknown')
        if(numvar.ge.3) then
           open(FILE__PE_PER,file='C1pe-per',form='formatted',status='unknown')
           open(FILE__CHI_FULL,file='C1chi-full',form='formatted',status='unknown')
           open(FILE__DIV_FULL,file='C1div-full',form='formatted',status='unknown')
        endif
     endif
     open(FILE__PSI_FULL,file='C1psi-full',form='formatted',status='unknown')
     open(FILE__J_FULL,file='C1J-full',form='formatted',status='unknown')
     if(numvar.ge.2) then
        open(FILE__I_FULL,file='C1I-full',form='formatted',status='unknown')
        if(numvar.ge.3) then
           open(FILE__PE_FULL,file='C1pe-full',form='formatted',status='unknown')
        endif
     endif
     if(idens.eq.1) then
        open(FILE__DENSITY,file='C1density',form='formatted',status='unknown')
     endif
     if(idebug.ge.1) &
          open(FILE__MATRIX,file='matrix.txt',form='formatted', status='unknown')

     write(FILE__C1NEW,1001) datec(1:4),datec(5:6),datec(7:8),                   &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("M3D-C1T VERSION 0.1    DATE: ",a4,1x,a2,1x,a2,3x,               &
          "TIME: ",a2,":",a2,":",a4,/)
  endif

end subroutine openf
!============================================================
subroutine output

  use p_data
  use t_data
  use basic
  use arrays
  use output_files

  implicit none
  integer :: i, j, i1, i3, ivertex, indexmid, numnodes
  real :: x, z, ans, ans2, etoto, etot, ediff, denom, gamma
  real :: etotd, etoth, error, enorm, vmaxsq, vnew, vmax, val1
  real :: dum1, val2, dum2, alx, alz

  call numnod(numnodes)
  call getboundingboxsize(alx, alz)

  ! calculate maximum perturbed current for printout
  
  ajmax = 0.
  do i=2,400
     do j=2,400
        x = (i-1)*alx/400.
        z = (j-1)*alz/400.
        call evaluate(x,z,ans,ans2,phi,2,numvar,0)
        ajmax = max(ans2,ajmax)
     enddo
  enddo

  ! search for the location of the magnetic axis and separatrices
  !      call axis(phi,xsep,zsep,1)
  xsep = 0.
  zsep = 0.
  etoto= ekino+emago
  etot = ekin + emag
  ediff = (etot-etoto)/dt
  write(*,*) etot, ekin, emag, ntime
  denom = dt*(ekin + ekino)
  if(denom.ne.0) gamma = (ekin - ekino)/denom
  etotd = .5*(ekind+emagd+ekindo+emagdo)
  etoth = .5*(ekinph +ekinth +ekin3h +emagph +emagth +emag3h        &
       +ekinpho+ekintho+ekin3ho+emagpho+emagtho+emag3ho)
  error =  ediff - etotd - etoth 
  if(myrank.eq.0) write(FILE__ERROR,2001) ntime,error,                &
       (ekinp-ekinpo)/dt,-ekinpd,-ekinph,                             &
       (emagp-emagpo)/dt,-emagpd,-emagph,                             &
       (ekint-ekinto)/dt,-ekintd,-ekinth,                             &
       (emagt-emagto)/dt,-emagtd,-emagth,                             &
       (ekin3-ekin3o)/dt,-ekin3d,-ekin3h,                             &
       (emag3-emag3o)/dt,-emag3d,-emag3h
  ! special debug write
  if(myrank.eq.0) write(FILE__ERROR,2001) ntime,error,                       &
       (ekinp-ekinpo)/dt,-ekinpdo,-ekinpho,                           &
       (emagp-emagpo)/dt,-emagpdo,-emagpho,                           &
       (ekint-ekinto)/dt,-ekintdo,-ekintho,                           &
       (emagt-emagto)/dt,-emagtdo,-emagtho,                           &
       (ekin3-ekin3o)/dt,-ekin3do,-ekin3ho,                           &
       (emag3-emag3o)/dt,-emag3do,-emag3ho
2001 format(/,i5,1pe12.4,3(/1p6e12.4))
  enorm = max( abs((ekin-ekino)/dt), abs((emag-emago)/dt), abs(etotd), abs(etoth) )

  graphit(ntime,1) = (ekin - ekino)/dt

  ! calculate the maximum poloidal velocity at a grid point
  vmaxsq = 0
  do ivertex=1,numnodes
     if(mod(ivertex,5).eq.1.or.mod(ivertex,5).eq.0) go to 100

     i1 = 6*numvar*(ivertex-1)
     i3 = i1 + 12
     vnew = vel(i1+2)**2 + vel(i1+3)**2
     if(numvar.ge.3) then
        vnew = vnew + vel(i3+2)**2 + vel(i3+3)**2                   &  
             +2.*(vel(i3+2)*vel(i1+3)-vel(i3+3)*vel(i1+2))
     endif
     vmaxsq = max(vmaxsq,vnew)

100  continue
  enddo
  vmax = sqrt(vmaxsq)
  graphit(ntime,2) = (emag - emago)/dt
  graphit(ntime,3) = ediff
  graphit(ntime,4) = gamma
  graphit(ntime,5) = ekind
  graphit(ntime,6) = emagd
  graphit(ntime,7) = etotd
  graphit(ntime,8) = time
  graphit(ntime,9) = xsep(1)
  graphit(ntime,10) = zsep(1)
  graphit(ntime,11) = xsep(2)
  graphit(ntime,12) = zsep(2)
  graphit(ntime,13) = xsep(3) 
  graphit(ntime,14) = zsep(3)
  graphit(ntime,15) = xsep(4)
  graphit(ntime,16) = zsep(4)
  graphit(ntime,17) = -etotd
  graphit(ntime,18) = error
  graphit(ntime,19) = (ekint - ekinto)/dt
  graphit(ntime,20) = (emagt - emagto)/dt
  graphit(ntime,21)= -(etoth)
  graphit(ntime,22) = -ekinth
  graphit(ntime,23) = emagph
  graphit(ntime,24) = emagth
  indexmid = 6*numvar*(numnodes - 1)/2 + 1
  !      graphit(ntime,25) = phi(indexmid)
  ! new definition 10/30/04
  call evaluate(0.0  ,alz/2,val1,dum1,phi,1,numvar,0)
  call evaluate(alx/2,alz/2,val2,dum2,phi,1,numvar,0)
  graphit(ntime,25) = 0.5*(val2-val1)
  if(ntime.gt.1) then
     graphit(ntime,26) = (graphit(ntime,25)-graphit(ntime-1,25))/dt
  endif
  graphit(ntime,27) = (ekin3 - ekin3o)/dt
  graphit(ntime,28) = (emag3 - emag3o)/dt
  graphit(ntime,29) = tflux
  graphit(ntime,30) = chierror
  graphit(ntime,31) = totcur
  graphit(ntime,32) = vmax


  ! check if it is time for plots and to write restart file
  if(mod(ntime,ntimepr).ne.0) then
     if(ntime.eq.0) then
        write(*,4110)
        write(FILE__C1NEW,4110)
     endif
     go to 50
  endif

  ! plot full perturbation
  if(linear.eq.1) then
     if(ntime.eq.0) then
        call plotit(vel0,phi0,1)
        if(idens.eq.1 .and. maxrank .eq. 1) call oneplot(den0,1,1,"n-0",0)
     endif
     call plotit(vel,phi,0)
     if(idens.eq.1 .and. maxrank .eq. 1) call oneplot(den,1,1,"n-p",1)
  else
     if(eqsubtract.eq.1) then
!        call plotit(vel,phi,0)
        call plotit(vel+vel0,phi+phi0,1)
        if(idens.eq.1 .and. maxrank .eq. 1) call oneplot(den+den0,1,1,"n-t",1)
     else
        call plotit(vel,phi,1)
        if(idens.eq.1 .and. maxrank .eq. 1) call oneplot(den,1,1,"n-t",1)
     endif
!!$     if(idens.eq.1) call oneplot(deni,1,1,"n^-1",0)
  endif
  
  ! plot linearized solution
  
  if(maxrank .eq. 1) call oneplot(sb1,1,1,"sb1 ",0)
  if(numvar.ge.2 .and. maxrank .eq. 1) call oneplot(sb2,1,1,"sb2 ",0)
  if(numvar.ge.3 .and. maxrank .eq. 1) then
     call oneplot(sp1,1,1,"sp1 ",0)
  endif
  if(idens.eq.1 .and. maxrank .eq. 1) call oneplot(deni,1,1,"deni",0)

  call wrrestart
  write(*,4110)
  write(FILE__C1NEW,4110)

4110 format(/, "cycle    time",                                        &
          "      ekin      rrate     ediff",                                &
          "     etotd      etoth     error       gamma")
50 continue

  if(ntime.ge.5 .and. facd.eq.0) then
     errori = errori + dt*abs(error)/abs(enorm)
     enormi = enormi + dt 
     if(enormi.ne.0) ratioi = errori/enormi
  endif

  write(*,4111) ntime,time,ekin,graphit(ntime,26),ediff,etotd,      &
       &              etoth,error,gamma,ratioi,vmax 
  write(FILE__C1NEW,4111) ntime,time,ekin,graphit(ntime,26),ediff,etotd,      &
       &              etoth,error,gamma,ratioi,vmax 
4111 format(i5,1p7e10.2,1pe14.6,1pe12.4,1pe14.6)
4411 format(i5,1p4e12.4)
  return

end subroutine output
!===================================================================
subroutine plotenergy(graphit,ntimep,maxts,mstart,numvarp)
  
  implicit none
  integer :: ntimep, maxts, mstart, numvarp, ntime, i
  real :: graphit(0:ntimep,*), tmin, tmax, ymin, ymax, ydiff
  character*80 :: s100

  tmin = graphit(1,8)
  tmax = graphit(maxts,8)
  ntime = maxts+1-mstart

  ymin = 0.
  ymax = 0.
  do i=mstart,maxts
     ymin = min(ymin,graphit(i,1),graphit(i,2),graphit(i,19),        &
          graphit(i,17),graphit(i,18),graphit(i,21),graphit(i,20))
     ymax = max(ymax,graphit(i,1),graphit(i,2),graphit(i,19),        &
          graphit(i,17),graphit(i,18),graphit(i,21),graphit(i,20))
  enddo
  if(numvarp.ge.3) then
     do i=mstart,maxts
        ymin = min(ymin,graphit(i,27),graphit(i,28))
        ymax = max(ymax,graphit(i,27),graphit(i,28))
     enddo
  endif
  ymax = min (ymax,2.)
      
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.200,1.0)
  call tracec(1hK, graphit(mstart,8), graphit(mstart,1),ntime, -1,-1,0.,0.)
  call tracec(1hM, graphit(mstart,8), graphit(mstart,2),ntime, -1,-1,0.,0.)
  call tracec(1hD, graphit(mstart,8), graphit(mstart,17),ntime, -1,-1,0.,0.)
  call tracec(1h*, graphit(mstart,8), graphit(mstart,18),ntime, -1,-1,0.,0.)
  call tracec(1hH,graphit(mstart,8), graphit(mstart,21),ntime, -1,-1,0.,0.)
  if(numvarp.ge.2) then
     call tracec(1hI, graphit(mstart,8), graphit(mstart,20),ntime, -1,-1,0.,0.)
     call tracec(1hV, graphit(mstart,8), graphit(mstart,19),ntime, -1,-1,0.,0.)
     if(numvarp.ge.3) then
        call tracec(1hC, graphit(mstart,8), graphit(mstart,27),ntime, -1,-1,0.,0.)
        call tracec(1hP, graphit(mstart,8), graphit(mstart,28),ntime, -1,-1,0.,0.)
     endif
  endif
  call frame(0)
  ymin = graphit(mstart,25)
  ymax = graphit(mstart,25)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,25))
     ymax = max(ymax,graphit(i,25))
  enddo
  ydiff = ymax - ymin
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.150,.50)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,25),ntime, -1,-1,0.,0.)
  call setch(2.,2.0,1,2,0,-1)
  write(s100,9025) graphit(maxts,25)
9025 format("  psi",1pe12.4)
  call gtext(s100,80,0)
  !
  ymin = graphit(mstart,26)
  ymax = graphit(mstart,26)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,26))
     ymax = max(ymax,graphit(i,26))
  enddo
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.650,1.0)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,26),ntime, -1,-1,0.,0.)
  call setch(2.,19.,1,2,0,-1)
  write(s100,9026) graphit(maxts,26)
9026 format("  psidot",1pe12.4)
  call gtext(s100,80,0)
  call frame(0)
  
  ymin = graphit(mstart,29)
  ymax = graphit(mstart,29)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,29))
     ymax = max(ymax,graphit(i,29))
  enddo
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.150,.50)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,29),ntime, -1,-1,0.,0.)
  call setch(2.,2.0,1,2,0,-1)
  write(s100,9029) graphit(maxts,29)
9029 format("  tflux",1pe12.4)
  call gtext(s100,80,0)
  
  if(numvarp.ge.3) then
     ymin = graphit(mstart,30)
     ymax = graphit(mstart,30)
     do i=mstart+1,maxts
        ymin = min(ymin,graphit(i,30))
        ymax = max(ymax,graphit(i,30))
     enddo
     ydiff = ymax - ymin
     ymax = ymax + .05*ydiff
     ymin = ymin - .05*ydiff
     call mapg(tmin,tmax,ymin,ymax,.142,.800,.650,1.0)
     call tracec(1h*,graphit(mstart,8), graphit(mstart,30),ntime, -1,-1,0.,0.)
     call setch(2.,19.,1,2,0,-1)
     write(s100,9030) graphit(maxts,30)
9030 format("  chierror",1pe12.4)
     call gtext(s100,80,0)
  endif

  call frame(0)

  ymin = graphit(mstart,31)
  ymax = graphit(mstart,31)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,31))
     ymax = max(ymax,graphit(i,31))
  enddo
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.650,1.0)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,31),ntime, -1,-1,0.,0.)
  call setch(2.,19.,1,2,0,-1)
  write(s100,9031) graphit(maxts,31)
9031 format("  totcur",1pe12.4)
  call gtext(s100,80,0)
  
  ymin = graphit(mstart,32)
  ymax = graphit(mstart,32)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,32))
     ymax = max(ymax,graphit(i,32))
  enddo
  ymax = min(ymax,2.)
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.150,0.5)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,32),ntime, -1,-1,0.,0.)
  call setch(2.,2.0,1,2,0,-1)
  write(s100,9032) graphit(maxts,32)
9032 format("   vmax ",1pe12.4)
  call gtext(s100,80,0)
  call frame(0)

  return

  !     plot positions of the critical points
  ymin = 0.
  ymax = 0.
  do i=1,ntime
     ymin = min(ymin,graphit(i,9),graphit(i,11),graphit(i,13),graphit(i,15))
     ymax = max(ymax,graphit(i,9),graphit(i,11),graphit(i,13),graphit(i,15))
  enddo
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,0.1,0.45,0.2,0.8)
  call tracec(1h1, graphit(mstart,8), graphit(mstart,9), ntime, -1,-1,0.,0.)
  call tracec(1h2, graphit(mstart,8), graphit(mstart,11),ntime, -1,-1,0.,0.)
  call tracec(1h3, graphit(mstart,8), graphit(mstart,13), ntime, -1,-1,0.,0.)
  call tracec(1h4, graphit(mstart,8), graphit(mstart,15), ntime, -1,-1,0.,0.)
  ymin = 0.
  ymax = 0.
  do i=1,ntime
     ymin = min(ymin,graphit(i,10),graphit(i,12), graphit(i,14),graphit(i,16))
     ymax = max(ymax,graphit(i,10),graphit(i,12), graphit(i,14),graphit(i,16))
  enddo
  ydiff = ymax - ymin
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff

  call mapg(tmin,tmax,ymin,ymax,.55,0.90,0.2,0.8)
  call tracec(1h1, graphit(mstart,8), graphit(mstart,10), ntime, -1,-1,0.,0.)
  call tracec(1h2, graphit(mstart,8), graphit(mstart,22), ntime, -1,-1,0.,0.)
  call tracec(1h3, graphit(mstart,8), graphit(mstart,14), ntime, -1,-1,0.,0.)
  call tracec(1h4, graphit(mstart,8), graphit(mstart,16), ntime, -1,-1,0.,0.)

  return
end subroutine plotenergy
!============================================================
subroutine plotit(vel,phi,ilin)

  use basic
  use output_files

  implicit none
  integer :: ilin, inum, ix, iz, iresmid, iresqtr, ii, ih
  real ::  vel(*), phi(*), x, z, plotmin, plotmax
  real ::  plotmin2, plotmax2
  real ::  cval(ires),plot(ires,ires),plot2(ires,ires)
  real ::  xval(ires),zval(ires),plotmidx(ires),plotmidz(ires)
  real ::plotqtrx(ires),plotqtrz(ires), alx, alz
  character*80 s100

  call getboundingboxsize(alx, alz)
!  call getmincoord(xmin,zmin)

  ! ilin=0 for perturbed plots
  ! ilin=1 for full solution plots

  do inum=1,numvar

! start of plotting coding
     call setch(1.,32.,1,2,0,-1)
     write(s100,9002) ntime
9002 format(i4)
     call gtext(s100,80,0)
!!$     if(inum.eq.1) then
!!$        if(ilin.eq.0) write(FILE__PSI_PER,1101) ntime,ilin,inum
!!$        if(ilin.eq.0) write(FILE__J_PER,1101) ntime,ilin,inum
!!$     endif
!!$     if(inum.eq.2) then
!!$        if(ilin.eq.0) write(FILE__I_PER,1101) ntime,ilin,inum
!!$     endif
!!$     if(inum.eq.3) then
!!$        if(ilin.eq.0) write(FILE__PE_PER,1101) ntime,ilin,inum
!!$     endif
!!$
!!$     if(inum.eq.1) then
!!$        if(ilin.eq.1) write(FILE__PSI_FULL,1101) ntime,ilin,inum
!!$        if(ilin.eq.1) write(FILE__J_FULL,1101) ntime,ilin,inum
!!$     endif
!!$     if(inum.eq.2) then
!!$        if(ilin.eq.1) write(FILE__I_FULL,1101) ntime,ilin,inum
!!$     endif
!!$     if(inum.eq.3) then
!!$        if(ilin.eq.1) write(FILE__PE_FULL,1101) ntime,ilin,inum
!!$     endif

1101 format("ntime, ilin  inum  =",3i5)
     
     do ix = 1,ires
        x = (ix-1)*alx/(ires-1)
        xval(ix) = x + xzero
        do iz = 1,ires
           z = (iz-1)*alz/(ires-1)
           zval(iz) = z
           call evaluate(x,z,plot(ix,iz),plot2(ix,iz),phi,inum,numvar,whichtri(ix,iz))
        enddo

        if(inum.eq.1) then
           if(ilin.eq.0) &
                write(FILE__PSI_PER,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
           if(ilin.eq.0) &
                write(FILE__J_PER,1102) zval(ix),xval(ix),(plot2(ix,iz),iz=1,ires)
        endif
        if(inum.eq.2) then
           if(ilin.eq.0) &
                write(FILE__I_PER,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
        endif
        if(inum.eq.3) then
           if(ilin.eq.0) &
                write(FILE__PE_PER,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
        endif

        if(inum.eq.1) then
           if(ilin.eq.1) &
                write(FILE__PSI_FULL,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
           if(ilin.eq.1) &
                write(FILE__J_FULL,1102) zval(ix),xval(ix),(plot2(ix,iz),iz=1,ires)
        endif
        if(inum.eq.2) then
           if(ilin.eq.1) &
                write(FILE__I_FULL,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
        endif
        if(inum.eq.3) then
           if(ilin.eq.1) &
                write(FILE__PE_FULL,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
        endif
     enddo

     iresmid = (ires-1)/2 + 1
     iresqtr = (ires-1)/4 + 1
     do ii = 1,ires
        plotmidx(ii) = plot(ii,iresmid)
        plotmidz(ii) = plot(iresmid,ii)
        plotqtrx(ii) = plot(ii,iresqtr)
        plotqtrz(ii) = plot(iresqtr,ii)
     enddo
      
     plotmin = plot(1,1)
     plotmax = plot(1,1)
     do ix=1,ires
        do iz=1,ires
           plotmin = min(plotmin,plot(ix,iz))
           plotmax = max(plotmax,plot(ix,iz))
        enddo
     enddo
!!$     if(plotmin.ge.plotmax-1.e-12) go to 100
     cval(1) = plotmin
     cval(2) = plotmax

     ! lower right
     call map(xzero,alx+xzero,0.,alz,.67,1.0,.170,.500)
     call rcontr(-25,cval,0,plot,ires,xval,1,ires,1,zval,1,ires,1)
     call mapg(xzero,alx+xzero,plotmin,plotmax,.67,1.0,.05,.170)
     call trace(xval,plotmidx,ires,-1,-1,0.,0.)
     call tracec(1hQ,xval,plotqtrx,ires,-1,-1,0.,0.)
     call mapg(plotmin,plotmax,0.,alz,.55,.670,.170,.500)
     call trace(plotmidz,zval,ires,-1,-1,0.,0.)
     call tracec(1hQ,plotqtrz,zval,ires,-1,-1,0.,0.)
     
     call setch(33.,4.0,1,2,0,-1)
     if(inum.eq.1) then
        if(ilin.eq.0)write(s100,9005)
        if(ilin.eq.1)write(s100,8005)
8005    format(" PSI-T")
9005    format(" PSI-P")
     endif
     if(inum.eq.2) then
        if(ilin.eq.0)write(s100,9015)
        if(ilin.eq.1)write(s100,8006)
8006    format(" I-T")
9015    format(" I-P")
     endif
     if(inum.eq.3) then
        if(ilin.eq.0)write(s100,9025)
        if(ilin.eq.1)write(s100,8007)
8007    format(" p-T")
9025    format(" p-P")
     endif

     call gtext(s100,80,0)
     write(s100,9006) plotmin
     call gtext(s100,80,0)
     write(s100,9008) plotmax
     call gtext(s100,80,0)
9006 format(1pe10.2)
9008 format(1pe10.2)

100  continue
     plotmin2 = plot2(1,1)
     plotmax2 = plot2(1,1)
     do ix=1,ires
        do iz=1,ires
           plotmin2 = min(plotmin2,plot2(ix,iz))
           plotmax2 = max(plotmax2,plot2(ix,iz))
        enddo
     enddo
     !      if(plotmin2 .ge. plotmax2) go to 210
     cval(1) = plotmin2
     cval(2) = plotmax2
     do ii = 1,ires
        plotmidx(ii) = plot2(ii,iresmid)
        plotmidz(ii)  = plot2(iresmid,ii)
        plotqtrx(ii) = plot2(ii,iresqtr)
        plotqtrz(ii) = plot2(iresqtr,ii)
     enddo

     ! lower left
     call map(xzero,alx+xzero,0.,alz,0.17,0.5,.17,.500)
     call rcontr(-25,cval,0,plot2,ires,xval,1,ires,1,zval,1,ires,1)
     
     call mapg(xzero,alx+xzero,plotmin2,plotmax2,.17,0.5,.05,.170)
     call trace(xval,plotmidx,ires,-1,-1,0.,0.)
     call tracec(1hQ,xval,plotqtrx,ires,-1,-1,0.,0.)
     call mapg(plotmin2,plotmax2,0.,alz,.05,.170,.170,.500)
     call trace(plotmidz,zval,ires,-1,-1,0.,0.)
     call tracec(1hQ,plotqtrz,zval,ires,-1,-1,0.,0.)
     call setch(1.,4.0,1,2,0,-1)
     if(ilin.eq.0) then
        if(itor.eq.0) then
           if(inum.eq.1) write(s100,9105)
           if(inum.eq.2) write(s100,8205)
           if(inum.eq.3) write(s100,8215)
        endif
        if(itor.eq.1) write(s100,9205)
     endif
     if(ilin.eq.1) then
        if(inum.eq.1)write(s100,8105)
        if(inum.eq.2)write(s100,8205)
        if(inum.eq.3) write(s100,8215)
     endif
9105 format("D^2 PSI")
9205 format("D^* PSI")
8105 format("J-T")
8205 format("D^2 I")
8215 format("D^2 p")
     call gtext(s100,80,0)
     write(s100,9106) plotmin2
9106 format(1pe10.2)
     call gtext(s100,80,0)
     write(s100,9106) plotmax2
     call gtext(s100,80,0)
210  continue
     if(itor.eq.1) go to 300
     
!!$     if(inum.eq.1) then
!!$        if(ilin.eq.1) write(FILE__PHI_FULL,1101) ntime,ilin,inum
!!$        if(ilin.eq.1) write(FILE__VOR_FULL,1101) ntime,ilin,inum
!!$     endif
!!$     if(inum.eq.2) then
!!$        if(ilin.eq.1) write(FILE__V_FULL,1101) ntime,ilin,inum
!!$     endif
!!$     if(inum.eq.3) then
!!$        if(ilin.eq.1) write(FILE__CHI_FULL,1101) ntime,ilin,inum
!!$        if(ilin.eq.1) write(FILE__DIV_FULL,1101) ntime,ilin,inum
!!$     endif
     do ix = 1,ires
        x = (ix-1)*alx/(ires-1)
        xval(ix) = x + xzero
        do iz = 1,ires
           z = (iz-1)*alz/(ires-1)
           zval(iz) = z
           call evaluate(x,z,plot(ix,iz),plot2(ix,iz),vel,inum,numvar,whichtri(ix,iz))
        enddo
        
        if(inum.eq.1) then
           write(FILE__PHI_FULL,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
           write(FILE__VOR_FULL,1102) zval(ix),xval(ix),(plot2(ix,iz),iz=1,ires)
        endif
        if(inum.eq.2) then
           write(FILE__V_FULL,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
        endif
        if(inum.eq.3) then
           write(FILE__CHI_FULL,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
           write(FILE__DIV_FULL,1102) zval(ix),xval(ix),(plot2(ix,iz),iz=1,ires)
        endif
1102    format(1p10e12.4)
     enddo

     plotmin = plot(1,1)
     plotmax = plot(1,1)
     do ix=1,ires
        do iz=1,ires
           plotmin = min(plotmin,plot(ix,iz))
           plotmax = max(plotmax,plot(ix,iz))
        enddo
     enddo
     do ii = 1,ires
        plotmidx(ii) = plot(ii,iresmid)
        plotmidz(ii) = plot(iresmid,ii)
        plotqtrx(ii) = plot(ii,iresqtr)
        plotqtrz(ii) = plot(iresqtr,ii)
     enddo
     
     if(plotmax .gt. plotmin) then
        cval(1) = plotmin
        cval(2) = plotmax
        
        ! upper right
        call map(xzero,alx+xzero,0.,alz,0.67,1.0,.670,1.00)
        call rcontr(-25,cval,0,plot,ires,xval,1,ires,1,zval,1,ires,1)

        call mapg(xzero,alx+xzero,plotmin,plotmax,.67,1.0,.55,.670)
        call trace(xval,plotmidx,ires,-1,-1,0.,0.)
        call tracec(1hQ,xval,plotqtrx,ires,-1,-1,0.,0.)
        call mapg(plotmin,plotmax,0.,alz,.55,.670,.670,1.00)
        call trace(plotmidz,zval,ires,-1,-1,0.,0.)
        call tracec(1hQ,plotqtrz,zval,ires,-1,-1,0.,0.)
        call setch(33.,20.,1,2,0,-1)
        if(inum.eq.1) then
           write(s100,9004)
9004       format(" PHI")
        endif
        if(inum.eq.2) then
           write(s100,9014)
9014       format("Vz")
        endif
        if(inum.eq.3) then
           write(s100,9024)
9024       format("Chi")
        endif
        
        call gtext(s100,80,0)
        write(s100,9007) plotmin
9007    format(1pe10.2)
        call gtext(s100,80,0)
        write(s100,9007) plotmax
        call gtext(s100,80,0)
     endif

     plotmin = plot2(1,1)
     plotmax = plot2(1,1)
     do ix=1,ires
        do iz=1,ires
           plotmin = min(plotmin,plot2(ix,iz))
           plotmax = max(plotmax,plot2(ix,iz))
        enddo
     enddo

     if(plotmax .gt. plotmin) then
        cval(1) = plotmin
        cval(2) = plotmax
        do ii = 1,ires
           plotmidx(ii) = plot2(ii,iresmid)
           plotmidz(ii)  = plot2(iresmid,ii)
           plotqtrx(ii) = plot2(ii,iresqtr)
           plotqtrz(ii) = plot2(iresqtr,ii)
        enddo

        ! upper left
        call map(xzero,alx+xzero,0.,alz,0.17,0.5,.670,1.00)
        call rcontr(-25,cval,0,plot2,ires,xval,1,ires,1,zval,1,ires,1)

        call mapg(xzero,alx+xzero,plotmin,plotmax,.17,.50,.55,.670)
        call trace(xval,plotmidx,ires,-1,-1,0.,0.)
        call tracec(1hQ,xval,plotqtrx,ires,-1,-1,0.,0.)
        call mapg(plotmin,plotmax,0.,alz,.05,.170,.670,1.000)
        call trace(plotmidz,zval,ires,-1,-1,0.,0.)
        call tracec(1hQ,plotqtrz,zval,ires,-1,-1,0.,0.)
        call setch(1.,20.,1,2,0,-1)
        if(inum.eq.1) then
           write(s100,9104)
9104       format("D^2 PHI")
        endif
        if(inum.eq.2) then
           write(s100,9114)
9114       format("D^2 Vz")
        endif
        if(inum.eq.3) then
           write(s100,9124)
9124       format("D^2 Chi")
        endif
        call gtext(s100,80,0)
        write(s100,9107) plotmin
9107    format(1pe10.2)
        call gtext(s100,80,0)
        write(s100,9107) plotmax
        call gtext(s100,80,0)
     endif
     go to 101

300  continue
     call mapg(xzero,alx+xzero,plotmin ,plotmax ,0.60,1.00,.600,1.00)
     ih = (ires-1)/2 + 1
     call trace(xval,plot(1,ih) ,ires,-1,-1,0.,0.)
     call mapg(xzero,alx+xzero,plotmin2,plotmax2,0.05,0.45,.600,1.00)
     call trace(xval,plot2(1,ih),ires,-1,-1,0.,0.)
101  call frame(0)
     ! end of plotting coding

     ! end if inum loop
  end do
  
7001 format(f7.2,5x,11f7.1)

end subroutine plotit
!============================================================
subroutine oneplot(dum,inum,numvare,label,iplot)
  use basic
  use output_files

  implicit none
  character*4 label

  integer, intent(in) :: iplot

  integer :: inum, numvare, ix, iz, iresmid, iresqtr, ii
  real :: x, z, plotmin, plotmax
  real :: alx, alz!, xmin, zmin
  real :: dum(*)
  real :: cval(ires),plot(ires,ires),plot2(ires,ires)
  real :: xval(ires),zval(ires),plotmidx(ires),plotmidz(ires)
  real :: plotqtrx(ires),plotqtrz(ires)
  character*80 :: s100

  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "called oneplot"

  call getboundingboxsize(alx, alz)
  !      call getmincoord(xmin, zmin)
     
  if(myrank.eq.0 .and. iprint.gt.0) then
     write(*,*) "ires, inum, numvare", ires, inum, numvare
  endif

100 continue
  do ix = 1,ires
     x = (ix-1)*alx/(ires-1)
     xval(ix) = x + xzero
     do iz = 1,ires
        z = (iz-1)*alz/(ires-1)
        zval(iz) = z
        call evaluate(x,z,plot(ix,iz),plot2(ix,iz),dum,inum,numvare,whichtri(ix,iz))
     enddo
     if(iplot.eq.1) then
        write(FILE__DENSITY,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,ires)
1102    format(1p10e12.4)
     endif
  enddo

  plotmin = plot(1,1)
  plotmax = plot(1,1)
  do ix=1,ires
     do iz=1,ires
        plotmin = min(plotmin,plot(ix,iz))
        plotmax = max(plotmax,plot(ix,iz))
     enddo
  enddo
  
  cval(1) = plotmin
  cval(2) = plotmax
  call map(xzero,alx+xzero,0.,alz,0.67,1.0,.670,1.00)
  call rcontr(-25,cval,0,plot,ires,xval,1,ires,1,zval,1,ires,1)
  
  iresmid = (ires-1)/2 + 1
  iresqtr = (ires-1)/4 + 1
  do ii = 1,ires
     plotmidx(ii) = plot(ii,iresmid)
     plotmidz(ii) = plot(iresmid,ii)
     plotqtrx(ii) = plot(ii,iresqtr)
     plotqtrz(ii) = plot(iresqtr,ii)
  enddo


  call mapg(xzero,alx+xzero,plotmin,plotmax, .67,1.0,.55,.670)
  call trace(xval,plotmidx,ires,-1,-1,0.,0.)
  call tracec(1hQ,xval,plotqtrx,ires,-1,-1,0.,0.)
  call mapg(plotmin,plotmax,0.,alz,.55,.670,.670,1.00)
  call trace(plotmidz,zval,ires,-1,-1,0.,0.)
  call tracec(1hQ,plotqtrz,zval,ires,-1,-1,0.,0.)
  call setch(33.,20.,1,2,0,-1)
  write(s100,9004) label
9004 format(a4)
  
  call gtext(s100,80,0)
  write(s100,9007) plotmin
9007 format(1pe10.2)
  call gtext(s100,80,0)
  write(s100,9007) plotmax
  call gtext(s100,80,0)

101 call frame(0)
  return

end subroutine oneplot
!===================================

subroutine precalc_whattri()

  use basic
  use arrays

  implicit none

  integer :: ix, iz
  real :: x, z, x1, z1, alx, alz

  call getboundingboxsize(alx, alz)

  do ix = 1,ires
     x = (ix-1)*alx/(ires-1)
     do iz = 1,ires
        z = (iz-1)*alz/(ires-1)
       
        call whattri(x,z,whichtri(ix, iz),x1,z1)
     enddo
  enddo

end subroutine precalc_whattri

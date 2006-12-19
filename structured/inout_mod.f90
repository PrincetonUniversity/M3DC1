#define REAL64 real

module inout_mod

implicit none

contains

subroutine plotit(vel,phi,ilin)
  use basic

  implicit none

#ifndef BIT64
  real, dimension(*) :: vel, phi
#else
  REAL64, dimension(*) :: vel, phi
#endif
  integer, intent(in) :: ilin

  real, dimension(irs,irs) :: plot, plot2
  real, dimension(irs) :: cval, xval, zval, &
       plotmidx, plotmidz, plotqtrx, plotqtrz
  character*80 :: s100
  character*13 :: strs
  character*3 ::  maxst,minst

  integer :: inum, ix, iz, numplots, ii, ih
  integer :: iresmid, iresqtr, jresmid, jresqtr
  real :: x, z, plotmin, plotmin2, plotmax, plotmax2, small

  if(iprint.gt.0) write(*,*) "start of plotit"
  maxst = 'max'
  minst = 'min'

!     ilin=0 for perturbed plots
!     ilin=1 for full solution plots

  do inum=1,numvar
     numplots = 0

     ! start of plotting coding
     call setch(1.,32.,1,2,0,-1)
     write(s100,9002) ntime
9002 format(i4)
     call gtext(s100,80,0)
     if(inum.eq.1) then
        if(ilin.eq.0) write(11,1101) ntime,ilin,inum
        if(ilin.eq.0) write(13,1101) ntime,ilin,inum
     endif
     if(inum.eq.2) then
        if(ilin.eq.0) write(15,1101) ntime,ilin,inum
     endif
     if(inum.eq.3) then
        if(ilin.eq.0) write(17,1101) ntime,ilin,inum
     endif

     if(inum.eq.1) then
        if(ilin.eq.1) write(12,1101) ntime,ilin,inum
        if(ilin.eq.1) write(14,1101) ntime,ilin,inum
     endif
     if(inum.eq.2) then
        if(ilin.eq.1) write(16,1101) ntime,ilin,inum
     endif
     if(inum.eq.3) then
        if(ilin.eq.1) write(18,1101) ntime,ilin,inum
     endif

1101 format("ntime, ilin  inum  =",3i5)

     do ix = 1,irs
        x = (ix-1)*alx*(1./(1.+isym))/(irs-1)
        xval(ix) = x + xzero
        do iz = 1,irs
           z = (iz-1)*alz*(1./(1.+jsym))/(irs-1)
           zval(iz) = z
           call evaluate(x,z,plot(ix,iz),plot2(ix,iz),phi,inum,numvar)
        enddo
        
        if(inum.eq.1) then
           if(ilin.eq.0) write(11,1102) zval(ix),xval(ix),             &
                (plot(ix,iz),iz=1,irs)
           if(ilin.eq.0) write(13,1102) zval(ix),xval(ix),             &
                (plot2(ix,iz),iz=1,irs)
        endif
        if(inum.eq.2) then
           if(ilin.eq.0) write(15,1102) zval(ix),xval(ix),             &
                (plot(ix,iz),iz=1,irs)
      endif
      if(inum.eq.3) then
         if(ilin.eq.0) write(17,1102) zval(ix),xval(ix),               &
              (plot(ix,iz),iz=1,irs)
      endif

      if(inum.eq.1) then
         if(ilin.eq.1) then
           write(12,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,irs)
!          psiarray(ix) = plot(ix,(irs-1)/2+1)
         endif
         if(ilin.eq.1) then
              write(14,1102) zval(ix),xval(ix),(plot2(ix,iz),iz=1,irs)
!             curarray(ix) = plot2(ix,(irs-1)/2+1)
              endif
      endif
      if(inum.eq.2) then
         if(ilin.eq.1) write(16,1102) zval(ix),xval(ix),                   &
              (plot(ix,iz),iz=1,irs)
      endif
      if(inum.eq.3) then
         if(ilin.eq.1) then
           write(18,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,irs)
!          pearray(ix) = plot(ix,(irs-1)/2+1)
         endif
      endif
      
   enddo
   
   ! output waveform for wave propagation test
   if((itaylor.eq.4).and.(inum.eq.1)) then
      if(ntime.eq.0) then
         write(40,*) irs
      endif
      write(40,*) (plot(ix,1),ix=1,irs)
   endif

   ! WRITE HDF5 files
#ifndef nohdf
   if(inum.eq.1) then
      if(ilin.eq.0) then
         strs = 'psi-perturbed'
         call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
         maf(0)=max(maf(0),maxval(plot))
         mif(0)=min(mif(0),minval(plot))
         call writeHDF5scalar(iframe-1,ihdf5,maf(0),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(0),minst)
         ihdf5 =     ihdf5+1
         
         strs = 'J-perturbed'
         call writeHDF5(iframe-1,ihdf5,plot2,irs,irs,strs)
         maf(2)=max(maf(2),maxval(plot2))
         mif(2)=min(mif(2),minval(plot2))
         call writeHDF5scalar(iframe-1,ihdf5,maf(2),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(2),minst)
         ihdf5 =     ihdf5+1
      endif
   endif
   
   if(inum.eq.2) then
      if(ilin.eq.0) then
         strs = 'I-perturbed'
         call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
         maf(4)=max(maf(4),maxval(plot))
         mif(4)=min(mif(4),minval(plot))
         call writeHDF5scalar(iframe-1,ihdf5,maf(4),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(4),minst)
         ihdf5 =     ihdf5+1
      endif
   endif

   if(inum.eq.3) then
      if(ilin.eq.0) then
         strs = 'Pe-perturbed'
         call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
         maf(6)=max(maf(6),maxval(plot))
         mif(6)=min(mif(6),minval(plot))
         call writeHDF5scalar(iframe-1,ihdf5,maf(6),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(6),minst)
         ihdf5 =     ihdf5+1
      endif
   endif

   if(inum.eq.1) then
      if(ilin.eq.1) then
         strs = 'psi-full'
         call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
         maf(1)=max(maf(1),maxval(plot))
         mif(1)=min(mif(1),minval(plot))
         call writeHDF5scalar(iframe-1,ihdf5,maf(1),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(1),minst)
         ihdf5 =     ihdf5+1
         
         strs = 'J-full'
         call writeHDF5(iframe-1,ihdf5,plot2,irs,irs,strs)
         maf(3)=max(maf(3),maxval(plot2))
         mif(3)=min(mif(3),minval(plot2))
         call writeHDF5scalar(iframe-1,ihdf5,maf(3),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(3),minst)
         ihdf5 =     ihdf5+1
      endif
   endif
   if(inum.eq.2) then
      if(ilin.eq.1) then
         strs = 'I-full'
         call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
         maf(5)=max(maf(5),maxval(plot))
         mif(5)=min(mif(5),minval(plot))
         call writeHDF5scalar(iframe-1,ihdf5,maf(5),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(5),minst)
         ihdf5 =     ihdf5+1
      endif
   endif
   if(inum.eq.3) then
      if(ilin.eq.1) then
         strs = 'Pe-full'
         call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
         maf(7)=max(maf(7),maxval(plot))
         mif(7)=min(mif(7),minval(plot))
         call writeHDF5scalar(iframe-1,ihdf5,maf(7),maxst)
         call writeHDF5scalar(iframe-1,ihdf5,mif(7),minst)
         ihdf5 =     ihdf5+1
      endif
   endif
   if(iprint.eq.1) write(*,*) "after HDF5 write"
   ! end of write HDF5 files
#endif


   if(isym.eq.0) then
      iresmid = (irs-1)/2 + 1
      iresqtr = (irs-1)/4 + 1
   else
      iresmid = 1.
      iresqtr = (irs-1)/2 + 1
   endif

   if(jsym.eq.0) then
      jresmid = (irs-1)/2 + 1
      jresqtr = (irs-1)/4 + 1
   else
      jresmid = 1.
      jresqtr = (irs-1)/2 + 1
   endif
   
   do ii = 1,irs
      plotmidx(ii) = plot(ii,jresmid)
      plotmidz(ii) = plot(iresmid,ii)
      plotqtrx(ii) = plot(ii,jresqtr)
      plotqtrz(ii) = plot(iresqtr,ii)
   enddo

   plotmin = plot(1,1)
   plotmax = plot(1,1)
   do ix=1,irs
      do iz=1,irs
         plotmin = min(plotmin,plot(ix,iz))
         plotmax = max(plotmax,plot(ix,iz))
      enddo
   enddo
   small = 1.e-12
    if(plotmin.ge.plotmax-small) go to 100
   cval(1) = plotmin
   cval(2) = plotmax
   numplots = numplots+1

   ! lower right
   if(jsym.eq.0) then
      call map(xzero,alxp+xzero,0.,alzp,.67,1.0,.170,.500)
      call rcontr(-25,cval,0,plot,irs,xval,1,irs,1,zval,1,irs,1+jsym)
      call mapg(xzero,alxp+xzero,plotmin,plotmax,.67,1.0,.05,.170)
      call trace(xval,plotmidx,irs,-1,-1,0.,0.)
      call tracec(1hQ,xval,plotqtrx,irs,-1,-1,0.,0.)
      call mapg(plotmin,plotmax,0.,alzp,.55,.670,.170,.500)
      call trace(plotmidz,zval,irs,-1,-1,0.,0.)
      call tracec(1hQ,plotqtrz,zval,irs,-1,-1,0.,0.)
   else  ! jsym=1 coding follows
      call map(xzero,alxp+xzero,0.,alzp,.67,1.0,.335,.500)
      call rcontr(-25,cval,0,plot,irs,xval,1,irs,1,zval,1,irs,2)
      call map(xzero,alxp+xzero,-alzp,0,.67,1.0,.170,.335)
      call rcontr(-25,cval,0,plot,irs,xval,1,irs,1,-zval,1,irs,2)
      call mapg(xzero,alxp+xzero,plotmin,plotmax,.67,1.0,.05,.170)
      call trace(xval,plotmidx,irs,-1,-1,0.,0.)
      call tracec(1hQ,xval,plotqtrx,irs,-1,-1,0.,0.)
      call mapg(plotmin,plotmax,-alzp,alzp,.55,.670,.170,.500)
      call trace(plotmidz,-zval,irs,-1,-1,0.,0.)
      call tracec(1hQ,plotqtrz,-zval,irs,-1,-1,0.,0.)
      call trace(plotmidz,zval,irs,-1,-1,0.,0.)
      call tracec(1hQ,plotqtrz,zval,irs,-1,-1,0.,0.)
   endif ! on jsym
   
   call setch(33.,4.0,1,2,0,-1)
   if(inum.eq.1) then
      if(ilin.eq.0)write(s100,9005)
      if(ilin.eq.1)write(s100,8005)
8005  format(" PSI-T")
9005  format(" PSI-P")
   endif
   if(inum.eq.2) then
      if(ilin.eq.0)write(s100,9015)
      if(ilin.eq.1)write(s100,8006)
8006  format(" I-T")
9015  format(" I-P")
   endif
   if(inum.eq.3) then
      if(ilin.eq.0)write(s100,9025)
      if(ilin.eq.1)write(s100,8007)
8007  format(" p-T")
9025  format(" p-P")
   endif
   
   
   call gtext(s100,80,0)
   write(s100,9006) plotmin
   call gtext(s100,80,0)
   write(s100,9008) plotmax
   call gtext(s100,80,0)
9006 format(1pe10.2)
9008 format(1pe10.2)
   
100 continue
   plotmin2 = plot2(1,1)
   plotmax2 = plot2(1,1)
   do ix=1,irs
      do iz=1,irs
         plotmin2 = min(plotmin2,plot2(ix,iz))
         plotmax2 = max(plotmax2,plot2(ix,iz))
      enddo
   enddo
    if(plotmin2 .ge. plotmax2) go to 210
   cval(1) = plotmin2
   cval(2) = plotmax2
   numplots = numplots+1
   do ii = 1,irs
      plotmidx(ii) = plot2(ii,jresmid)
      plotmidz(ii)  = plot2(iresmid,ii)
      plotqtrx(ii) = plot2(ii,jresqtr)
      plotqtrz(ii) = plot2(iresqtr,ii)
   enddo
   
   ! lower left
   call map(xzero,alxp+xzero,0.,alzp,0.17,0.5,.17,.500)
   call rcontr(-25,cval,0,plot2,irs,xval,1,irs,1,zval,1,irs,1+jsym)
   
   call mapg(xzero,alxp+xzero,plotmin2,plotmax2,.17,0.5,.05,.170)
   call trace(xval,plotmidx,irs,-1,-1,0.,0.)
   call tracec(1hQ,xval,plotqtrx,irs,-1,-1,0.,0.)
   call mapg(plotmin2,plotmax2,0.,alzp,.05,.170,.170,.500)
   call trace(plotmidz,zval,irs,-1,-1,0.,0.)
   call tracec(1hQ,plotqtrz,zval,irs,-1,-1,0.,0.)
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
210 continue
   if(itor.eq.1) go to 300
   
   if(inum.eq.1) then
      if(ilin.eq.1) write(19,1101) ntime,ilin,inum
      if(ilin.eq.1) write(20,1101) ntime,ilin,inum
   endif
   if(inum.eq.2) then
      if(ilin.eq.1) write(21,1101) ntime,ilin,inum
   endif
   if(inum.eq.3) then
      if(ilin.eq.1) write(22,1101) ntime,ilin,inum
      if(ilin.eq.1) write(23,1101) ntime,ilin,inum
   endif
   do ix = 1,irs
      x = (ix-1)*alx*(1./(1.+isym))/(irs-1)
      xval(ix) = x + xzero
      do iz = 1,irs
         z = (iz-1)*alz*(1./(1.+jsym))/(irs-1)
         zval(iz) = z
         call evaluate(x,z,plot(ix,iz),plot2(ix,iz),vel,inum,numvar)
      enddo
      
      if(inum.eq.1) then
         if(ilin.eq.1) write(19,1102) zval(ix),xval(ix),                   &
              (plot(ix,iz),iz=1,irs)
         if(ilin.eq.1) write(20,1102) zval(ix),xval(ix),                   &
              (plot2(ix,iz),iz=1,irs)
      endif
      if(inum.eq.2) then
         if(ilin.eq.1) write(21,1102) zval(ix),xval(ix),                   &
              (plot(ix,iz),iz=1,irs)
      endif
      if(inum.eq.3) then
        if(ilin.eq.1) write(22,1102) zval(ix),xval(ix),                   &
             (plot(ix,iz),iz=1,irs)
        if(ilin.eq.1) write(23,1102) zval(ix),xval(ix),                   &
             (plot2(ix,iz),iz=1,irs)
     endif
1102 format(1p10e12.4)
  enddo

#ifndef nohdf
  ! write HDF5 files
  if(inum.eq.1) then
     if(ilin.eq.0) then
        strs = 'phi-full'
        call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
        maf(8)=max(maf(8),maxval(plot))
        mif(8)=min(mif(8),minval(plot))
        call writeHDF5scalar(iframe-1,ihdf5,maf(8),maxst)
        call writeHDF5scalar(iframe-1,ihdf5,mif(8),minst)
        ihdf5 =     ihdf5+1
        
        strs = "vor-full"
        call writeHDF5(iframe-1,ihdf5,plot2,irs,irs,strs)
        maf(9)=max(maf(9),maxval(plot2))
        mif(9)=min(mif(9),minval(plot2))
        call writeHDF5scalar(iframe-1,ihdf5,maf(9),maxst)
        call writeHDF5scalar(iframe-1,ihdf5,mif(9),minst)
        ihdf5 =     ihdf5+1
     endif
  endif
  if(inum.eq.2) then
     if(ilin.eq.0) then
        strs = 'V-full'
        call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
        maf(10)=max(maf(10),maxval(plot))
        mif(10)=min(mif(10),minval(plot))
        call writeHDF5scalar(iframe-1,ihdf5,maf(10),maxst)
        call writeHDF5scalar(iframe-1,ihdf5,mif(10),minst)
        ihdf5 =     ihdf5+1
     endif
  endif
  if(inum.eq.3) then
     if(ilin.eq.0) then
        strs = 'chi-full'
        call writeHDF5(iframe-1,ihdf5,plot,irs,irs,strs)
        maf(11)=max(maf(11),maxval(plot))
        mif(11)=min(mif(11),minval(plot))
        call writeHDF5scalar(iframe-1,ihdf5,maf(11),maxst)
        call writeHDF5scalar(iframe-1,ihdf5,mif(11),minst)
        ihdf5 =     ihdf5+1
        
        strs = 'div-full'
        call writeHDF5(iframe-1,ihdf5,plot2,irs,irs,strs)
        maf(12)=max(maf(12),maxval(plot2))
        mif(12)=min(mif(12),minval(plot2))
        call writeHDF5scalar(iframe-1,ihdf5,maf(12),maxst)
        call writeHDF5scalar(iframe-1,ihdf5,mif(12),minst)
        ihdf5 =     ihdf5+1
     endif
  endif
#endif
  
  plotmin = plot(1,1)
  plotmax = plot(1,1)
  do ix=1,irs
     do iz=1,irs
        plotmin = min(plotmin,plot(ix,iz))
        plotmax = max(plotmax,plot(ix,iz))
     enddo
  enddo
  do ii = 1,irs
     plotmidx(ii) = plot(ii,jresmid)
     plotmidz(ii) = plot(iresmid,ii)
     plotqtrx(ii) = plot(ii,jresqtr)
     plotqtrz(ii) = plot(iresqtr,ii)
  enddo
  
  if(plotmax .gt. plotmin) then
     cval(1) = plotmin
     cval(2) = plotmax
     numplots = numplots+1
     
     ! upper right
     call map(xzero,alxp+xzero,0.,alzp,0.67,1.0,.670,1.00)
     call rcontr(-25,cval,0,plot,irs,xval,1,irs,1,zval,1,irs,1+jsym)
     
     call mapg(xzero,alxp+xzero,plotmin,plotmax,.67,1.0,.55,.670)
     call trace(xval,plotmidx,irs,-1,-1,0.,0.)
     call tracec(1hQ,xval,plotqtrx,irs,-1,-1,0.,0.)
     call mapg(plotmin,plotmax,0.,alzp,.55,.670,.670,1.00)
     call trace(plotmidz,zval,irs,-1,-1,0.,0.)
     call tracec(1hQ,plotqtrz,zval,irs,-1,-1,0.,0.)
     call setch(33.,20.,1,2,0,-1)
     if(inum.eq.1) then
        write(s100,9004)
9004    format(" PHI")
     endif
     if(inum.eq.2) then
        write(s100,9014)
9014    format("Vz")
     endif
     if(inum.eq.3) then
        write(s100,9024)
9024    format("Chi")
     endif
     
     call gtext(s100,80,0)
     write(s100,9007) plotmin
9007 format(1pe10.2)
     call gtext(s100,80,0)
     write(s100,9007) plotmax
     call gtext(s100,80,0)
  endif
  
  plotmin = plot2(1,1)
  plotmax = plot2(1,1)
  do ix=1,irs
     do iz=1,irs
        plotmin = min(plotmin,plot2(ix,iz))
        plotmax = max(plotmax,plot2(ix,iz))
     enddo
  enddo

  if(plotmax .gt. plotmin) then
     cval(1) = plotmin
     cval(2) = plotmax
     numplots = numplots+1
     do ii = 1,irs
        plotmidx(ii) = plot2(ii,jresmid)
        plotmidz(ii) = plot2(iresmid,ii)
        plotqtrx(ii) = plot2(ii,jresqtr)
        plotqtrz(ii) = plot2(iresqtr,ii)
     enddo

     ! upper left
     call map(xzero,alxp+xzero,0.,alzp,0.17,0.5,.670,1.00)
     call rcontr(-25,cval,0,plot2,irs,xval,1,irs,1,zval,1,irs,1+jsym)
     
     call mapg(xzero,alxp+xzero,plotmin,plotmax,.17,.50,.55,.670)
     call trace(xval,plotmidx,irs,-1,-1,0.,0.)
     call tracec(1hQ,xval,plotqtrx,irs,-1,-1,0.,0.)
     call mapg(plotmin,plotmax,0.,alzp,.05,.170,.670,1.000)
     call trace(plotmidz,zval,irs,-1,-1,0.,0.)
     call tracec(1hQ,plotqtrz,zval,irs,-1,-1,0.,0.)
     call setch(1.,20.,1,2,0,-1)
     if(inum.eq.1) then
        write(s100,9104)
9104    format("D^2 PHI")
     endif
     if(inum.eq.2) then
        write(s100,9114)
9114    format("D^2 Vz")
     endif
     if(inum.eq.3) then
        write(s100,9124)
9124    format("D^2 Chi")
     endif
     call gtext(s100,80,0)
     write(s100,9107) plotmin
9107 format(1pe10.2)
     call gtext(s100,80,0)
     write(s100,9107) plotmax
     call gtext(s100,80,0)
  endif
  go to 101

  
300 continue
  call mapg(xzero,alxp+xzero,plotmin ,plotmax ,0.60,1.00,.600,1.00)
  ih = (irs-1)/2 + 1
  call trace(xval,plot(1,ih) ,irs,-1,-1,0.,0.)
  call mapg(xzero,alxp+xzero,plotmin2,plotmax2,0.05,0.45,.600,1.00)
  call trace(xval,plot2(1,ih),irs,-1,-1,0.,0.)
101 if(numplots.gt.0) call frame(0)
  ! end of plotting coding
  
901 enddo  ! on inum 
return

7001 format(f7.2,5x,11f7.1)
end subroutine plotit

!============================================================
subroutine oneplot(lu,iplot,dum,inum,numvare,label)
!
!....lu:  if nonzero, writes file on that logical unit
!....iplot:  1-4...which quadrant plot appears, if=1, calls frame
!....dum:    numvar=numvare array, location inum being plotted
!....label:  5 character label
!
  use basic

  integer, intent(in) :: lu, iplot, inum, numvare
#ifndef BIT64
  real, dimension(*) :: dum
#else
  REAL64, dimension(*) :: dum
#endif
  character*5, intent(in) :: label

  real :: cval(irs),plot(irs,irs),plot2(irs,irs),                 &
       xval(irs),zval(irs),plotmidx(irs),plotmidz(irs),                 &
       plotqtrx(irs),plotqtrz(irs)
  character*80 :: s100

  integer :: ix, iz, iresmid, iresqtr, jresmid, jresqtr, ii
  real :: x, z, plotmin, plotmax, x1, x2, z1, z2

      go to (1,2,3,4),iplot
 1    continue
      x1 = 0.67
      x2 = 1.0
      z1 = 0.670
      z2 = 1.00
      go to 5
 2    continue
      x1 = 0.67
      x2 = 1.0
      z1 = 0.170
      z2 = 0.50
      go to 5
 3    continue
      x1 = 0.17
      x2 = 0.50
      z1 = 0.170
      z2 = 0.50
      go to 5
 4    continue
      x1 = 0.17
      x2 = 0.50
      z1 = 0.670
      z2 = 1.00
 5    continue
 
      if(lu.ne.0) then
        write(lu,1101) ntime, inum
      endif

1101 format("ntime,  inum  =",2i5)
  do ix = 1,irs
     x = (ix-1)*alx*(1./(1.+isym))/(irs-1)
     xval(ix) = x + xzero
     do iz = 1,irs
        z = (iz-1)*alz*(1./(1.+jsym))/(irs-1)
        zval(iz) = z
        call evaluate(x,z,plot(ix,iz),plot2(ix,iz),dum,inum,numvare)
     enddo
     if(lu.ne.0) then
        write(lu,1102) zval(ix),xval(ix),(plot(ix,iz),iz=1,irs)
!       if(lu.eq.24) densarray(ix) = plot(ix,(irs-1)/2+1)
!       if(lu.eq.27) etajarray(ix) = plot(ix,(irs-1)/2+1)
!       if(lu.eq.28) ephiarray(ix) = plot(ix,(irs-1)/2+1)
!       if(lu.eq.29)  jcbarray(ix) = plot(ix,(irs-1)/2+1)
!       if(lu.eq.30) hyprarray(ix) = plot(ix,(irs-1)/2+1)
!       if(lu.eq.32)  vcbarray(ix) = plot(ix,(irs-1)/2+1)
1102    format(1p10e12.4)
     endif
  enddo

  plotmin = plot(1,1)
  plotmax = plot(1,1)
  do ix=1,irs
     do iz=1,irs
        plotmin = min(plotmin,plot(ix,iz))
        plotmax = max(plotmax,plot(ix,iz))
     enddo
  enddo
  
  if(plotmax .gt. plotmin + 1.e-10) then
     cval(1) = plotmin
     cval(2) = plotmax
     call map(xzero,alxp+xzero,0.,alzp,x1,x2,z1,z2)
     call rcontr(-25,cval,0,plot,irs,xval,1,irs,1,zval,1,irs,1+jsym)
     
     if(isym.eq.0) then
        iresmid = (irs-1)/2 + 1
        iresqtr = (irs-1)/4 + 1
     else
        iresmid = 1.
        iresqtr = (irs-1)/2 + 1
     endif

     if(jsym.eq.0) then
        jresmid = (irs-1)/2 + 1
        jresqtr = (irs-1)/4 + 1
     else
        jresmid = 1.
        jresqtr = (irs-1)/2 + 1
     endif
     
     do ii = 1,irs
        plotmidx(ii) = plot(ii,jresmid)
        plotmidz(ii) = plot(iresmid,ii)
        plotqtrx(ii) = plot(ii,jresqtr)
        plotqtrz(ii) = plot(iresqtr,ii)
     enddo
     
     call mapg(xzero,alxp+xzero,plotmin,plotmax,x1,x2,z1-0.12,z1)
     call trace(xval,plotmidx,irs,-1,-1,0.,0.)
     call tracec(1hQ,xval,plotqtrx,irs,-1,-1,0.,0.)
     call mapg(plotmin,plotmax,0.,alzp,x1-0.12,x1,z1,z2)
     call trace(plotmidz,zval,irs,-1,-1,0.,0.)
     call tracec(1hQ,plotqtrz,zval,irs,-1,-1,0.,0.)
     if(iplot.le.1) call setch(34.,20.,1,2,0,-1)
     if(iplot.eq.2) call setch(34., 4.,1,2,0,-1)
     if(iplot.eq.3) call setch( 3., 4.,1,2,0,-1)
     if(iplot.eq.4) call setch( 3.,20.,1,2,0,-1)

     write(s100,9004) label
9004 format(a4)
     
     call gtext(s100,80,0)
     write(s100,9007) plotmin
9007 format(1pe10.2)
     call gtext(s100,80,0)
     write(s100,9007) plotmax
     call gtext(s100,80,0)
  endif   ! on plotmax .gt. plotmin + 1.e-10

101  if(iplot.le.1) call frame(0)
  return
  
end subroutine oneplot

!============================================================


subroutine openf

  ! open output files and initialize NCAR graphics
  use p_data
  use t_data
  use basic
  use arrays
  
  implicit none

  integer :: maxhdf, jj

  if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Entering openf."
  endif

  maxhdf = 0

  if (myrank.eq.0) then !Serialize I/O
     call ncarcgm(1,'C1new.cgm')
     print *, "done ncarcgm"
     call dders(-1)
     print *, "done dders"

#ifndef nohdf
     !  define evenly spaced coordinates and create new HDF5 file
     xary = (/ ((alxp*jj)/(irs-1), jj=0,irs-1) /)
     yary = (/ ((alzp*jj)/(irs-1), jj=0,irs-1) /)
     if(numvar.eq.1) then
        if(linear.eq.1) maxhdf=4
        if(linear.eq.0) maxhdf=6
     endif

     if(numvar.eq.2) then
        if(linear.eq.1) maxhdf=6
        if(linear.eq.0) maxhdf=9
     endif

     if(numvar.eq.3) then
        if(linear.eq.1) maxhdf=9
        if(linear.eq.0) maxhdf=13
     endif

     call createHDF5(xary,irs,yary,irs,maxhdf)
     print *, "done createHDF5"
#endif

     ! initialize minimum and maximum
     maf = -1.e20
     mif =  1.e20

     ! open ascii output files
     if(itaylor.eq.4) then
        open(40,file="C1wave.out",form='formatted',status='unknown')
     endif
     open(9, file='C1new.out',form='formatted',status='unknown')
     open(65, file='C1error.out',form='formatted',status='unknown')
     open(66, file='C1ener.out',form='formatted',status='unknown')
     open(11,file='C1psi-per',form='formatted',status='unknown')
     open(13,file='C1J-per',form='formatted',status='unknown')
     open(19,file='C1phi-full',form='formatted',status='unknown')
     open(20,file='C1vor-full',form='formatted',status='unknown')
      if(linear.eq.0 .and. itaylor.eq.3) then
       open(25,file='C1VxBU',form='formatted',status='unknown')
       open(26,file='C1VxBC',form='formatted',status='unknown')
       open(27,file='C1ETAJ',form='formatted',status='unknown')
       open(28,file='C1EPHI',form='formatted',status='unknown')
       open(29,file='C1JxB',form='formatted',status='unknown')
       open(30,file='C1HYPR',form='formatted',status='unknown')
       open(32,file='C1VxBT',form='formatted',status='unknown')
      endif
     if(numvar.ge.2) then
        open(15,file='C1I-per',form='formatted',status='unknown')
        open(21,file='C1v-full',form='formatted',status='unknown')
        if(numvar.ge.3) then
           open(17,file='C1pe-per',form='formatted',status='unknown')
           open(22,file='C1chi-full',form='formatted',status='unknown')
           open(23,file='C1div-full',form='formatted',status='unknown')
        endif
     endif
     open(12,file='C1psi-full',form='formatted',status='unknown')
     open(14,file='C1J-full',form='formatted',status='unknown')
     if(numvar.ge.2) then
        open(16,file='C1I-full',form='formatted',status='unknown')
        if(numvar.ge.3) then
           open(18,file='C1pe-full',form='formatted',status='unknown')
        endif ! on numvar.ge.3
     endif ! on numvar.ge.2
     if(idens.eq.1) then
        open(24,file='C1density',form='formatted',status='unknown')
     endif
     if(idebug.ge.1) open(31,file='matrix.txt',form='formatted',    &
          status='unknown')
     write(9,*) version
     write(9,1001) datec(1:4),datec(5:6),datec(7:8),                   &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("DATE: ",a4,1x,a2,1x,a2,3x,          &
          "TIME: ",a2,":",a2,":",a4,/)
  endif

  if(itor.ne.0 .and. ntridim.lt.ntri) then
     write(*,5222) itor, ntridim, ntri
5222 format(" ERROR:  ntridim must be .eq. ntri for itor.ne.0",      &
          "   itor,ntridim,ntri =", 3i6)
     call safestop(5222)
  endif

  if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Exiting openf."
  endif


  return
end subroutine openf

! ==========================================================
subroutine output
  use p_data
  use t_data
  use basic
  use arrays
  use superlu

  implicit none
  integer :: i, j, ivertex, irect, jrect, i1, i3, indexmid, ix
  real :: x, z, ans, ans2, gamma
  real :: etot, etoto, etotd, etoth, ediff, error, enorm, denom, percerr
  real :: vmaxsq, vnew, vmax, x1, x2, z1, z2, val1, dum1, val2, dum2, superlutime

      if(iprint.ge.1) write(*,*) ntime,  "output called"
  ! calculate maximum perturbed current for printout
  ajmax = 0.
  do i=2,400
     do j=2,400
        x = (i-1)*alx*(1./(1.+isym))/400.
        z = (j-1)*alz*(1./(1.+jsym))/400.
        call evaluate(x,z,ans,ans2,phi,2,numvar)
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
  denom = dt*(ekin + ekino)
  if(denom.ne.0) gamma = (ekin - ekino)/denom
  etotd = .5*(ekind+emagd+ekindo+emagdo)
  etoth = .5*(ekinph +ekinth +ekin3h +emagph +emagth +emag3h        &
       +ekinpho+ekintho+ekin3ho+emagpho+emagtho+emag3ho)
  
  ! NOTE:  changed 1/15/06 when Ohmic & Viscous Heating  added to Pressure
  if(numvar.ne.3) then
     error = ediff - etotd - etoth
  else
     error = ediff - ekinph - ekinth
  endif
  if(myrank.eq.0) write(65,2001) ntime,error,                       &
       (ekinp-ekinpo)/dt,-ekinpd,-ekinph,                             &
       (emagp-emagpo)/dt,-emagpd,-emagph,                             &
       (ekint-ekinto)/dt,-ekintd,-ekinth,                             &
       (emagt-emagto)/dt,-emagtd,-emagth,                             &
       (ekin3-ekin3o)/dt,-ekin3d,-ekin3h,                             &
       (emag3-emag3o)/dt,-emag3d,-emag3h
  enorm = max( abs((ekin-ekino)/dt), abs((emag-emago)/dt),          &
       abs( etotd) , abs(etoth) )

  graphit(ntime,1) = (ekin - ekino)/dt

  ! calculate the maximum poloidal velocity at a grid point
  vmaxsq = 0
  ivertex = 0
  do jrect=1,m
     do irect=1,n
        ivertex = ivertex + 1
        if((jsym.eq.0 .and. jrect.eq.1) .or. jrect.eq.m) go to 100
        i1 = 6*numvar*(ivertex-1)
        i3 = i1 + 12
        vnew = vel(i1+2)**2 + vel(i1+3)**2
        if(numvar.ge.3) then
           vnew = vnew + vel(i3+2)**2 + vel(i3+3)**2                  &  
                +2.*(vel(i3+2)*vel(i1+3)-vel(i3+3)*vel(i1+2))
        endif
        vmaxsq = max(vmaxsq,vnew)
100     continue
     enddo
  enddo
  vmax = sqrt(vmaxsq)
  graphit(ntime,2) = (emag - emago)/dt
  graphit(ntime,3) = ediff
  graphit(ntime,4) = gamma
  graphit(ntime,5) = ekind
  graphit(ntime,6) = emagd
  graphit(ntime,7) = -.5*(emagd+emagdo)
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
  indexmid = 6*numvar*(n*m - 1)/2 + 1
  ! graphit(ntime,25) = phi(indexmid)
  ! new definition 10/30/04
  if(isym.eq.0) then
     x1 = alx
     x2 = alx/2.
  else
     x1 = alx/2.
     x2 = 0.
  endif
  if(jsym.eq.0) then
     z1 = alz/2.
     z2 = alz/2.
  else
     z1 = 0.
     z2 = 0.
  endif
  call evaluate(x1,z1,val1,dum1,phi,1,numvar)
  call evaluate(x2,z2,val2,dum2,phi,1,numvar)
  graphit(ntime,25) = 0.5*(val2-val1)
      graphit(ntime,49) = val2
  if(ntime.gt.1) then
     graphit(ntime,26) = (graphit(ntime,25)-graphit(ntime-1,25))/dt
      graphit(ntime,50) = (graphit(ntime,49)-graphit(ntime-1,49))/dt
  endif
  graphit(ntime,27) = (ekin3 - ekin3o)/dt
  graphit(ntime,28) = (emag3 - emag3o)/dt
  graphit(ntime,29) = tflux
  graphit(ntime,30) = chierror
  graphit(ntime,31) = totcur
  graphit(ntime,32) = vmax
  graphit(ntime,33) = ttotal-tfirst
  graphit(ntime,34) = maxrank
  graphit(ntime,35) = tread
  graphit(ntime,36) = telements + tread
  graphit(ntime,37) = tsolve + telements + tread
  graphit(ntime,38) = tonestep
  graphit(ntime,39) = tsolve + telements + tread + tadotb + tinit    &
       + tmpi
  graphit(ntime,40) = tsolve + telements + tread + tadotb
  graphit(ntime,41) = tsolve + telements + tread + tadotb + tinit
  
  graphit(ntime,42) = ekinp+ekint+ekin3
  graphit(ntime,43) = emagp + emagt
  graphit(ntime,44) = emag3
  graphit(ntime,45) = ekinp+ekint+ekin3+emagp+emagt+emag3
!....added 9/13/06
  graphit(ntime,46) = ntime
  graphit(ntime,47) = dt
  graphit(ntime,48) = thimp
!
!...note.....locations 49 and 50 defined above
!

  ! subtract off initial conditions
  if(ntime.gt.0) then
     graphit(ntime,42) = graphit(ntime,42) - graphit(0,42)
     graphit(ntime,43) = graphit(ntime,43) - graphit(0,43)
     graphit(ntime,44) = graphit(ntime,44) - graphit(0,44)
     graphit(ntime,45) = graphit(ntime,45) - graphit(0,45)
     
     ! special debug write
     etot = graphit(ntime,45)
     enorm=abs(ekinp)+abs(emagp)+abs(ekint)+abs(emagt)                 &
          +abs(ekin3)+abs(emag3)
     percerr = etot/enorm
     if(myrank.eq.0) write(66,2002) ntime,time,                        &
          ekinp,emagp,ekint,emagt,ekin3,emag3,etot,ekin,percerr
      superlutime = graphit(ntime  ,37) - graphit(ntime  ,36)                  &
                  -(graphit(ntime-1,37) - graphit(ntime-1,36))
      write(*,2003) ntime,superlutime
 2003 format("superlutime for timestep", i5, 1pe12.4)
2001 format(/,i5,1pe12.4,3(/1p6e12.4))
2002 format(i5,1p10e12.4)
  endif
!
!
!.....compute some midplane arrays for elvis graphics
     z = alz/2
     do ix = 1,irs
        x = (ix-1)*alx*(1./(1.+isym))/(irs-1)
        call evaluate(x,z,psiarray(ix),dum1,phi,1,numvar)
        call evaluate(x,z,curarray(ix),dum1,jphi,1,1)
        call evaluate(x,z,densarray(ix),dum1,dent,1,1)
        call evaluate(x,z,pearray(ix),dum1,phi,3,numvar)
        if(elvis.gt.0) then
          call evaluate(x,z,ephiarray(ix),dum1,ephi,1,1)
          call evaluate(x,z,jcbarray(ix),dum1,eph5,1,1)
          call evaluate(x,z,vcbarray(ix),dum1,eph7,1,1)
          call evaluate(x,z,etajarray(ix),dum1,eph4,1,1)
          call evaluate(x,z,hyprarray(ix),dum1,eph6,1,1)
        endif
     enddo
!
!.....end of computing midplane arrays for elvis graphics
!
  if(ntime.eq.0) return
  ! check if it is time for plots and to write restart file
  if(mod(ntime,ntimepr).ne.0) then
     if(ntime.eq.0) then
        write(*,4110)
        write(9,4110)
     endif
     go to 50
  endif
  iframe = iframe + 1
  ihdf5 = 0
      if(iprint.gt.0) write(*,*) "before first call to plotit"

  ! plot full perturbation
  phi1 = phi + phi0
  if(linear.eq.0) call plotit(vel,phi1,1)
  if(linear.eq.0 .and. jper.eq.1) call plotit(vel,phi,0)
  !plot equilibrium for linear run
  if(linear.eq.1.and.ntime.eq.0) call plotit(vel0,phi0,1)
  !plot linearized solution
  if(linear.eq.1) call plotit(vel,phi,0)
      if(iprint.gt.0) write(*,*) "before first call to oneplot"
  if(idens.eq.1) then
    call oneplot(0,2,dent,1,1,"dent ")
    call oneplot(0,3,deni,1,1,"deni ")
    call oneplot(0,4,den0,1,1,"den0 ")
    call oneplot(24,1,den,1,1,"den  ")
  endif
!
!....plot source terms
  if(numvar.ge.2) call oneplot(0,2,sphi2,1,1,"sph2 ")
  if(numvar.ge.3) then
     call oneplot(0,3,sphip,1,1,"sphp ")
     call oneplot(0,4,sphik,1,1,"sphk ")
  endif
  call oneplot(0,1,sphi1,1,1,"sph1 ")
!
! debug plotting of ohmic and viscous source terms
  if(numvar.ge.2) call oneplot(0,3,ohmic  ,2,numvar,"ohm2 ")
  if(numvar.ge.3) call oneplot(0,4,ohmic  ,3,numvar,"ohm3 ")
  if(numvar.ge.3) call oneplot(0,2,viscous,3,numvar,"visc ")
                  call oneplot(0,1,ohmic  ,1,numvar,"ohm1 ")
!
!....special diagnostic plot of toroidal electric field added 06/04/06...SCJ
      if(linear.eq.0 .and. itaylor.eq.3) then
!     call oneplot(25,2,eph2,1,1,"VxBU ")
!     call oneplot(26,3,eph3,1,1,"VxBC ")
!     call oneplot(27,4,eph4,1,1,"etaJ ")
!     call oneplot(28,1,ephi,1,1,"ephi ")

!     call oneplot(29,2,eph5,1,1,"JxB  ")
!     call oneplot(30,3,eph6,1,1,"hypr ")
!     call oneplot(32,4,eph7,1,1,"VxBT ")
!     call oneplot(0,1,ephi,1,1,"ephi ")
      endif
!

  if(gyro.eq.1) then
     call oneplot(0,2,bsqr,1,1,"b^2  ")
     call oneplot(0,1,bsqri,1,1,"b^-2 ")
  endif
  ! end of special DEBUG Plotting
   if(ipres.eq.1) then
     call oneplot(0,2,sphip,1,1,"sphp ")
     if(linear.ne.1) then
        call oneplot(0,3,pres,1,1,"pres ")
        call oneplot(0,1,pres+pres0,1,1,'p-t  ')
     else
      call oneplot(0,3,pres0,1,1,'p-0  ')
        call oneplot(0,1,pres,1,1,"p-p  ")
     endif
  endif


4444 format(1p10e12.4)
#ifndef nohdf
  call writeHDF5time(iframe-1,time)
#endif
  call wrrestart
  write(*,4110)
  write(9,4110)
4110 format(/,"cycle    time",                            &
          "      ekin      rrate     ediff",              &
          "     etotd      etoth     error       gamma")
  nlast = ntime
50 continue
  
  if(ntime.ge.5 .and. facd.eq.0) then
     errori = errori + dt*abs(error)/abs(enorm)
     enormi = enormi + dt 
     if(enormi.ne.0) ratioi = errori/enormi
  endif

  write(*,4111) ntime,time,ekin,graphit(ntime,26),ediff,etotd,      &
       etoth,error,gamma,ratioi,vmax 
  write(9,4111) ntime,time,ekin,graphit(ntime,26),ediff,etotd,      &
       etoth,error,gamma,ratioi,vmax 
4111 format(i5,1p7e10.2,1pe14.6,1pe12.4,1pe14.6)
4411 format(i5,1p4e12.4)
  return
end subroutine output
!============================================================================
subroutine plotenergy(graphit,ntimep,maxts,mstart,numvarp)

  implicit none

  real, intent(in), dimension(0:ntimep,*) :: graphit
  integer, intent(in) :: ntimep, maxts, mstart, numvarp
  character*80 :: s100
  integer :: ntime, i
  real :: tmin, tmax
  real :: ymin, ymax, ydiff

! first plot incremental energy change each cycle
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
  ymax = min (ymax,4.)
      
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.200,1.0)
  call tracec(1hK, graphit(mstart,8), graphit(mstart,1) ,ntime, -1,-1,0.,0.)
  call tracec(1hM, graphit(mstart,8), graphit(mstart,2) ,ntime, -1,-1,0.,0.)
  call tracec(1hD, graphit(mstart,8), graphit(mstart,17),ntime, -1,-1,0.,0.)
  call tracec(1hR, graphit(mstart,8), graphit(mstart,7), ntime, -1,-1,0.,0.)
  call tracec(1h*, graphit(mstart,8), graphit(mstart,18),ntime, -1,-1,0.,0.)
  call tracec(1hH,graphit(mstart,8), graphit(mstart,21), ntime, -1,-1,0.,0.)
  if(numvarp.ge.2) then
     call tracec(1hI, graphit(mstart,8), graphit(mstart,20),ntime, -1,-1,0.,0.)
     call tracec(1hV, graphit(mstart,8), graphit(mstart,19),ntime, -1,-1,0.,0.)
     if(numvarp.ge.3) then
        call tracec(1hC, graphit(mstart,8), graphit(mstart,27),ntime, &
             -1,-1,0.,0.)
        call tracec(1hP, graphit(mstart,8), graphit(mstart,28),ntime, &
             -1,-1,0.,0.)
     endif
  endif
  call frame(0)

  ! now plot cumulative energies
  ymin = 0.
  ymax = 0.
  do i=mstart,maxts
     ymin = min(ymin,graphit(i,42),graphit(i,43),graphit(i,44),        &
          graphit(i,45))
     ymax = max(ymax,graphit(i,42),graphit(i,43),graphit(i,44),        &
          graphit(i,45))
  enddo

  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.200,1.0)
  call tracec(1hK, graphit(mstart,8), graphit(mstart,42) ,ntime,     &
       -1,-1,0.,0.)
  call tracec(1hM, graphit(mstart,8), graphit(mstart,43) ,ntime,     &
       -1,-1,0.,0.)
  call tracec(1hP, graphit(mstart,8), graphit(mstart,44),ntime,     &
       -1,-1,0.,0.)
  call tracec(1hE, graphit(mstart,8), graphit(mstart,45),ntime,      &
       -1,-1,0.,0.)
  call frame(0)
  
  ymin = min(graphit(mstart,25),graphit(mstart,49))
  ymax = max(graphit(mstart,25),graphit(mstart,49))
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,25),graphit(i,49))
     ymax = max(ymax,graphit(i,25),graphit(i,49))
  enddo
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.150,.50)
  call tracec(1hD,graphit(mstart,8), graphit(mstart,25),ntime,-1,-1,0.,0.)
  call tracec(1hX,graphit(mstart,8), graphit(mstart,49),ntime,-1,-1,0.,0.)
  call setch(2.,2.0,1,2,0,-1)
  write(s100,9025) graphit(maxts,25),graphit(maxts,49)
9025 format("  psi",1p2e12.4)
  call gtext(s100,80,0)
  
  ymin = min(graphit(mstart,26),graphit(mstart,50))
  ymax = max(graphit(mstart,26),graphit(mstart,50))
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,26),graphit(i,50))
     ymax = max(ymax,graphit(i,26),graphit(i,50))
  enddo
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.650,1.0)
  call tracec(1hD,graphit(mstart,8), graphit(mstart,26),ntime,-1,-1,0.,0.)
  call tracec(1hX,graphit(mstart,8), graphit(mstart,50),ntime,-1,-1,0.,0.)
  call setch(2.,19.,1,2,0,-1)
  write(s100,9026) graphit(maxts,26),graphit(maxts,50)
9026 format("  psidot",1p2e12.4)
  call gtext(s100,80,0)
  call frame(0)
  
  ymin = graphit(mstart,29)
  ymax = graphit(mstart,29)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,29))
     ymax = max(ymax,graphit(i,29))
  enddo
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.150,.50)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,29),ntime,-1,-1,0.,0.)
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
     ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
     ymax = ymax + .05*ydiff
     ymin = ymin - .05*ydiff
     call mapg(tmin,tmax,ymin,ymax,.142,.800,.650,1.0)
     call tracec(1h*,graphit(mstart,8), graphit(mstart,30),ntime,-1,-1,0.,0.)
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
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.650,1.0)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,31),ntime,-1,-1,0.,0.)
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
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.150,0.5)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,32),ntime,-1,-1,0.,0.)
  call setch(2.,2.0,1,2,0,-1)
  write(s100,9032) graphit(maxts,32)
9032 format("   vmax ",1pe12.4)
  call gtext(s100,80,0)
  call frame(0)
  
  ymin = 0.
  ymax = graphit(mstart,33)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,33))
     ymax = max(ymax,graphit(i,33))
  enddo
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(graphit(mstart,46),graphit(ntime,46),ymin,ymax,.142,.800,.650,1.0)
  call tracec(1h*,graphit(mstart,46), graphit(mstart,33),ntime,-1,-1,0.,0.)
  call tracec(1hR,graphit(mstart,46), graphit(mstart,35),ntime,-1,-1,0.,0.)
  call tracec(1hE,graphit(mstart,46), graphit(mstart,36),ntime,-1,-1,0.,0.)
  call tracec(1hS,graphit(mstart,46), graphit(mstart,37),ntime,-1,-1,0.,0.)
  call tracec(1hO,graphit(mstart,46), graphit(mstart,38),ntime,-1,-1,0.,0.)
  call tracec(1hM,graphit(mstart,46), graphit(mstart,39),ntime,-1,-1,0.,0.)
  call tracec(1hI,graphit(mstart,46), graphit(mstart,41),ntime,-1,-1,0.,0.)
  call tracec(1hA,graphit(mstart,46), graphit(mstart,40),ntime,-1,-1,0.,0.)
  call setch(2.,19.,1,2,0,-1)
  write(s100,9035) graphit(maxts,33)
9035 format("  cputim",1pe12.4)
  call gtext(s100,80,0)
  
  ymin = 1.
  ymax = 2.
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,34))
     ymax = max(ymax,graphit(i,34))
  enddo
  call mapg(graphit(mstart,46),graphit(ntime,46),ymin,ymax,.142,.800,.150,0.5)
  call tracec(1h*,graphit(mstart,46), graphit(mstart,34),ntime,-1,-1,0.,0.)
  call setch(2.,2.0,1,2,0,-1)
  write(s100,9034) graphit(maxts,34)
9034 format("no proc.",1pe12.4)
  call gtext(s100,80,0)
  call frame(0)
!...new graph added 9/13/06
  ymin = 0.
  ymax = graphit(mstart,47)
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,47))
     ymax = max(ymax,graphit(i,47))
  enddo
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.650,1.0)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,47),ntime,-1,-1,0.,0.)
  call setch(2.,19.,1,2,0,-1)
  write(s100,9033) graphit(maxts,47)
9033 format("  timestep",1pe12.4)
  call gtext(s100,80,0)
  
  ymin = 0.
  ymax = 1.1
  do i=mstart+1,maxts
     ymin = min(ymin,graphit(i,48))
     ymax = max(ymax,graphit(i,48))
  enddo
  call mapg(tmin,tmax,ymin,ymax,.142,.800,.150,0.5)
  call tracec(1h*,graphit(mstart,8), graphit(mstart,48),ntime,-1,-1,0.,0.)
  call setch(2.,2.0,1,2,0,-1)
  write(s100,9036) graphit(maxts,48)
9036 format("theta-imp",1pe12.4)
  call gtext(s100,80,0)
  call frame(0)
  
  return
  
  ! plot positions of the critical points
  ymin = 0.
  ymax = 0.
  do i=1,ntime
     ymin = min(ymin,graphit(i,9),graphit(i,11),                     &
          graphit(i,13),graphit(i,15))
     ymax = max(ymax,graphit(i,9),graphit(i,11),                     &
          graphit(i,13),graphit(i,15))
  enddo
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff
  call mapg(tmin,tmax,ymin,ymax,0.1,0.45,0.2,0.8)
  call tracec(1h1, graphit(mstart,8), graphit(mstart,9),ntime,-1,-1,0.,0.)
  call tracec(1h2, graphit(mstart,8), graphit(mstart,11),ntime,-1,-1,0.,0.)
  call tracec(1h3, graphit(mstart,8), graphit(mstart,13),ntime,-1,-1,0.,0.)
  call tracec(1h4, graphit(mstart,8), graphit(mstart,15),ntime,-1,-1,0.,0.)
  ymin = 0.
  ymax = 0.
  do i=1,ntime
     ymin = min(ymin,graphit(i,10),graphit(i,12),                    &
          graphit(i,14),graphit(i,16))
     ymax = max(ymax,graphit(i,10),graphit(i,12),                    &
          graphit(i,14),graphit(i,16))
  enddo
  ydiff = ymax - ymin + .01*abs(ymax) + 1.e-8
  ymax = ymax + .05*ydiff
  ymin = ymin - .05*ydiff

  call mapg(tmin,tmax,ymin,ymax,.55,0.90,0.2,0.8)
  call tracec(1h1, graphit(mstart,8), graphit(mstart,10) ,ntime,-1,-1,0.,0.)
  call tracec(1h2, graphit(mstart,8), graphit(mstart,22) ,ntime,-1,-1,0.,0.)
  call tracec(1h3, graphit(mstart,8), graphit(mstart,14),ntime,-1,-1,0.,0.)
  call tracec(1h4, graphit(mstart,8), graphit(mstart,16),ntime,-1,-1,0.,0.)
  
  return
end subroutine plotenergy

!============================================================
subroutine input
  use basic

  implicit none
  !
!.....set itest .ne. 0 for a particular test problem
!     itest = 1 for tilting columns (nonlinear, n=m=21)
!     itest = 2 for Taylor problem  (n=m=31, dt=2  )
!     itest = 3 for tilting columns (linear)
!     itest = 4 for Taylor problem  (Fitzpatrick base case)
!     itest = 5 for tilting columns - two fluid (linear)
!

  if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Entering input."
  endif

  ! switch for second order time advance:
  isecondorder = 0

  ! linear parameter  0-nonlinear,  1-linear
  linear = 0

  ! switch for subtracting out the equilibrium
  eqsubtract=0

  ! density advance parameter: 0 no advance,  1-advance density
  idens = 1

  ! restart parameter 0-no restart, 1 restart
  irestart = 0

  ! associated restart parameter: 0-normal, 1-start time from zero
  istart = 0
  
  ! setup parameter 0-read from disk, 1 compute
  isetup = 1
  
  ! boundary parameter  0-nonperiodic, 1-periodic
  iper = 0
  jper = 0
  
  ! switch for Taylor Problem
  itaylor = 0

  ! switch for semi-implicit method (needs bzerosi)
  isemii = 0
  ! symmetry variables:  0-no symmetry,   1-symmetry
  isym = 0
  jsym = 0
  
  ! implicitness parameter (thimp > 0.5 for stability)
  thimp = 0.5
  
  ! hyper-resistivity coefficients
  hyper = 0.5000
  hyperi= 0.5000
  hyperv= 0.50
  hyperp= 0.50

  ! viscosity, resistivity, and heat conduction
  amu = .00010
  etar = .00010
  kappa = .001
  denm = .0005
  
  ! ratio of specific heats
  gam = 5./3.
  
  ! temperature diffusion
  kappat = 0.
  
  ! toroidal magnetic field
  bzero = 1.
  ! semi-implicit magnetic field
  bzerosi = 0.

  facw = 1.
  ! multiplies the diffusion [test] terms 
  facd = 1.
  ! NOTE:  if facw=1., facd is zeroed after the first cycle

  ! 2-fluid coefficients
  cb = 0.000
  db = 0.00
  ! regularization coefficient
  regular = 1.e-7
  
  ! masking switch (0 no mask,  1 mask hyper terms at boundary
  imask = 0

  ! maximum in taylor expansion loop (5 is maximum)
  maxs = 5
  ! how many cycles to skip before inverting matrix
  nskip = 1

  ! timestep
  dt = 0.1
  ! ntimemax is max time cycles, nimepr is cycles between print cycle
  
  ntimemax = 20
  ntimepr   = 5

  ! dimensions and number of vertices in x and z directions
  alx = 4.
  alz = 4.

  n = 21
  m = 21

  ! actual (physical) coordinates of bottom, left corner of mesh
  xzero = 10.
  zzero = -2.

  ! ratio of fluid pressure to magnetic pressure for test problem
  beta = 0.

  ! gravitational force
  grav = 0.

  ! itor=0 for cylinder, 1 for torus
  itor=0
  
  ! parameters needed for toroidal plasma
  tcuro = 1.
  xmag = 12.1
  zmag = 0.
  xlim = 10.5
  zlim = 0.

  ! pressure = p0*(1 + p1*psi + p2*psi**2 + ...)
  !                where psi is the normalized poloidal flux
  p0 = 0.01      
  p1 = -1.
  p2 = -2.
  ! ion pressure
  pi0 = 0.

  ! derivative of current wrt normalized poloidal flux
  djdpsi = 0.0
  
  ! number of velocity and field variables
  numvar = 3

  ! parameters for taylor problem (itaylor=1)
  eps = .01
  tau = 1.
  
  if(itest.eq.0) then
     open(5,file='C1input',form='formatted',status='old')
     read(5,nml=inputnl)
     if(amuc.eq.0) amuc = amu
     go to 100
  endif

  go to (1,2,3,4,5,6,7,8,9,10), itest

1 continue
  ! itest=1 tilting columns (nonlinear)
  itaylor = 0
  iper = 0
  alx = 4.
  alz = 4.
  dt = 0.05
  facd = 1.
  etar = .0010
  amu = .001
  kappa = .001
  hyper = 0.01
  hyperi= 0.01
  hyperv = 0.01
  hyperc = 0.01
  hyperp= 0.50
  cb = 1.
  db = .0001
  ntimemax = 100
  ntimepr = 1
  thimp = 0.55
  n = 31
  m = 31
  isetup = 0
  imask = 1
  numvar = 2
  bzero = 1.
  go to 100

2 continue
  ! itest=2 Taylor problem, high resistivity
  itaylor = 1
  iper = 1
  alx = 8.
  alz = 2.
  dt = 0.5
  facd = 0.
  etar = .0001
  amu  = .00010
  kappa = .0001
  hyper = 1.00
  hyperi = 1.00
  hyperv = 1.0
  hyperc = 1.0
  hyperp= 0.50
  cb = 1.
  db = 0.0
  ntimemax = 80
  ntimepr   = 10
  n = 61
  m = 61
  thimp = 0.6
  isetup = 1
  numvar = 2
  imask = 0
  bzero = 0.
  go to 100
3 continue
  ! itest=3  tilting columns--linear
  linear = 1
  itaylor = 0
  iper = 0
  alx = 4.
  alz = 4.
  dt = 0.10
  facd = 1.
  etar = .0001
  amu = .0001
  kappa = .0001
  hyper = 0.
  hyperi = 0.
  hyperv = 0.0
  hyperc = 0.0
  cb = 1.0
  db = 0.
  ntimemax = 80
  ntimepr = 5
  thimp = 0.80
  n = 15
  m = 15
  isetup = 0
  imask = 0

  go to 100
4 continue
  ! itest=4 Taylor problem..Fitzpatrick base case with eta=mu=1.e-5
  itaylor = 1
  iper = 1
  alx = 8.
  alz = 2.
  dt = 2.0
  facd = 0.
  etar = .000010
  amu  = .0000100
  kappa = .00001
  hyper = 0.
  hyperi = 0.
  hyperv = 0.0
  hyperc = 0.0
  hyperp= 0.50
  cb = 0.
  db = 0.
  ntimemax = 40
  ntimepr   = 5
  n = 41
  m = 41
  thimp = 0.60
  isetup = 1
  go to 100
5 continue
  ! itest=5 tilting columns--linear-two fluid
  linear = 1
  numvar = 2
  itaylor = 0
  iper = 0
  alx = 4.
  alz = 4.
  dt = 0.05
  facd = 1.
  etar = .001
  amu = .001
  kappa = .001
  hyper = 1.00
  hyperi = 1.00
  hyperv = 4.0
  hyperc = 4.0
  hyperp= 1.0
  cb = 1.0
  db = 0.2
  ntimemax = 81
  ntimepr = 1
  thimp = 0.51
  n = 61
  m = 61
  isetup = 0
  irestart = 0  !0=start from scratch;  1=restart
  imask = 0
  bzero = 1.
  beta = 0.0
  
  go to 100
6 continue
  ! itest=6 Taylor problem..eta=mu=1.e-5, with two-fluid terms
  itaylor = 1
  iper = 1
  alx = 8.
  alz = 2.
  dt = 2.0
  facd = 0.
  etar = 0.00010
  amu  = 0.0000100
  kappa= 0.00001
  hyper = .0010
  hyperi = .0010
  hyperv = .0010
  hyperc = .0010
  hyperp= 0.50
  cb = 1.0
  db = 0.
  ntimemax = 20
  ntimepr   = 10
  n = 41
  m = 41
  thimp = 0.6
  isetup = 0
  imask = 1
  go to 100
  
7 continue
  ! itest=7 tilting columns--non-linear-two fluid
  itaylor = 0
  iper = 0
  alx = 4.
  alz = 4.
  dt = 0.025
  facd = 1.
  etar = .0010
  amu = .0010
  kappa = .001
  hyper =  1.
  hyperi = 1.
  hyperv = 1.
  hyperc = 1.
  hyperp= 0.50
  cb = 1.0
  db = .4
  ntimemax = 5
  ntimepr = 5
  thimp = 0.55
  n = 31
  m = 31
  isetup = 0
  imask = 0
  numvar = 2
  bzero = 0.001
  go to 100
  
8 continue
  ! itest=8 Taylor state
  itaylor = 2
  iper = 0
  alx = 4.
  alz = 4.
  dt = 0.18
  facd = 1.
  etar = .0001
  amu = .0001
  kappa = .0001
  hyper = .00
  hyperi = 0.0
  hyperv = 0.0
  hyperc = 0.0
  hyperp= 0.50
  cb = 0.3
  db = 0.
  ntimemax = 10
  ntimepr = 2
  thimp = 0.6
  n = 41
  m = 41
  isetup = 0
  imask = 0
  go to 100
  
9 continue
  ! itest=9 GEM reconnection
  linear = 0
  itaylor = 3
  iper = 1
  alx = 25.6
  alz = 12.8
  dt = .25
  facd = 0.
  etar = .0010
  amu = .0010
  kappa = .001
  hyper = 1.0  !hyper-resistivity coefficient
  hyperi = 1.0 !hyper-resistivity coefficient for I
  hyperv = 2.0 !hyper-viscosity coefficient for V
  hyperc = 0.0 !hyper-viscosity coefficient for chi
  hyperp= 0.50
  cb = 1.0
  db = 0.0001
  ntimemax = 150
  ntimepr = 50
  thimp = .60
  n = 61
  m = 61
  isetup = 1
  imask = 1
  bzero = 0.
  numvar = 2
  eps = 0.1
  go to 100

10 continue
  ! itest=10 GEM reconnection...(two fluid)
  linear = 0
  itaylor = 3
  iper = 1
  alx = 25.6
  alz = 12.8
  dt = 0.25    !time step size
  facd = 0.
  etar = .0010
  amu = .0010
  kappa = .001
  hyper = 2.00  !hyper-resistivity coefficient
  hyperi = 2.00 !hyper-resistivity coefficient for I
  hyperv = 4.0 !hyper-viscosity coefficient for V
  hyperc = 4.0 !hyper-viscosity coefficient for chi
  hyperp= 1.00
  cb = 1.0
  db = 1.0
  ntimemax = 50  !Number of time steps
  ntimepr = 10
  thimp = 0.6
  n = 61
  m = 61
  isetup = 0  !0=already exists; 1=recalc
  irestart = 1  !0=start from scratch;  1=restart
  imask = 1
  bzero = 0.0
  numvar = 3
  eps = 0.1
  go to 100
100 continue

  if(itaylor.eq.7) then
     grav = 3.*p0/ln
     if(myrank.eq.0) write(*,*) "setting grav =", grav
  endif
  
  if(myrank.eq.0) write(*,nml=inputnl)

    if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Exiting input."
  endif

  return
end subroutine input

!============================================================
subroutine readit
  use p_data
  use t_data
  use basic
  use arrays

  implicit none

#ifdef mpi
  include 'mpif.h'

  integer count, root, ier
  integer nread,mread,ioddr,iperr,jperr
  real alxr,alzr
#endif 

  integer :: iodd, i, j, k, l, itri, itype, ir1, ir2, ir3

  if(myrank.eq.0) then
     open(55,file='setup',form='unformatted',status='old')
     
     read(55) nread,mread,ioddr,iperr,jperr,alxr,alzr
     if(nread.ne.n) then
        write(*,*) "error in setup file, nread,n,itest = ",nread,n,itest
     endif
     if(mread.ne.m .or. ioddr.ne.ioddm .or. iperr.ne.iper.or.          &
          jperr.ne.jper .or.  alxr.ne.alx .or. alzr.ne.alz) then
        write(*,*) "error in setup file",mread,m,ioddr,ioddm,           &
             iperr,iper,jperr,jper,alxr,alx,alzr,alz
     endif
  endif ! on myrank.eq.0

#ifdef mpi
  ! send rsq to the other tasks
  root=0
  count=1
  
  CALL MPI_Bcast(nread,count,MPI_INTEGER,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(mread,count,MPI_INTEGER,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(ioddr,count,MPI_INTEGER,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(iperr,count,MPI_INTEGER,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(jperr,count,MPI_INTEGER,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(alxr, count,MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(alzr, count,MPI_DOUBLE_PRECISION, root,MPI_COMM_WORLD, ier)
#endif
  if(nread.ne.n.or.                                               &
       mread.ne.m .or. ioddr.ne.ioddm .or. iperr.ne.iper.or.        &
       jperr.ne.jper .or.  alxr.ne.alx .or. alzr.ne.alz) then
     call safestop(21)
  endif

  if(myrank.eq.0) then
     read(55) ((d2term(iodd,i),iodd=1,2),i=1,18)
     read(55) (((aterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
     read(55) (((dterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
     read(55) (((bterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
     read(55) (((sterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
     read(55) (((xterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
     read(55) (((yterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
     
     read(55) ((((k0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((k1term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((k2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g4term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g5term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g6term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g7term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g8term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g9term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g10erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g11erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g12erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g13erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g14erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g15erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((g16erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((h3term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((h5term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((x0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((y0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((x1term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((y1term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((x2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
     read(55) ((((y2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)

     if(itor.ne.0) read(55) (((hterm(itri,i,j),itri=1,ntridim),i=1,18),j=1,18)

     close(55)
  endif  ! myrank.eq.0

#ifdef mpi
  ! send rsq to the other tasks
  root=0
  count=2*18
  CALL MPI_Bcast(d2term, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  
  count=2*18*18
  CALL MPI_Bcast(aterm, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(dterm, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(bterm, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(sterm, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(xterm, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(yterm, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  if(iprint.eq.1) then
     call MPI_Barrier(MPI_COMM_WORLD, ier) 
     write(*,112) myrank, (sterm(1,1,i),i=1,4)
     write(9+myrank,112) myrank,  (sterm(1,1,i),i=1,4)
  endif
112 format(i3, "readit-2 ",1p4e12.4 )
  call MPI_Barrier(MPI_COMM_WORLD, ier)

  count=2*18*18*18
  CALL MPI_Bcast(k0term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(k1term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(k2term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g0term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g2term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g4term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g5term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g6term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g7term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g8term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g9term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g10erm,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g11erm,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g12erm,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g13erm,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g14erm,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g15erm,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(g16erm,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(h3term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(h5term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(x0term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(y0term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(x1term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(y1term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(x2term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  CALL MPI_Bcast(y2term,count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  call MPI_Barrier(MPI_COMM_WORLD, ier)
  
  if(itor.ne.0) then
     count=ntridim*18*18
     CALL MPI_Bcast(hterm, count,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD, ier)
  endif !(itor.ne.0)

#endif 
  
  
  allocate (ttermij(18,18))
  do iodd=1,ioddmx
     do itype=1,ntensor
        if(iread(itype).eq.0) go to 400
        if(iodd.eq.1) then
           if(myrank.eq.0) then
              open(itype+200,file=filename(itype), form='unformatted',      &
                   status='old')
              rewind(itype+200)
           endif ! myrank.eq.0

           open(itype+100,form='unformatted',status='scratch')
           rewind(itype+100)
        endif

        ! read tensors from submission directory and save to scratch disk
        do i=1,18
           do j=1,18
              if(myrank.eq.0) then
                 read(itype+200,err=1001,end=1001) ir1,ir2,ir3,              &
                      ((ttermij(k,l),l=1,18),k=1,18)
                 if(ir1.ne.iodd .or. ir2.ne.i .or. ir3.ne.j) go to 1001
              endif !if(myrank.eq.0)

#ifdef mpi 
              root=0

              count=18*18
              CALL MPI_Bcast(ttermij,count,MPI_DOUBLE_PRECISION,            &
                   root,MPI_COMM_WORLD, ier)

#endif

              write(itype+100) iodd,i,j,((ttermij(k,l),l=1,18),k=1,18)
              
           enddo ! on j
        enddo ! on i
400  enddo ! on itype
500 enddo ! on iodd
  do itype=1,ntensor
     if(myrank.eq.0) then
        if(iread(itype).ne.0) then
           close(UNIT=itype+200,STATUS='KEEP')
        endif
     endif !if(myrank.eq.0)
  enddo
  write(*,1002) myrank
1002 format(" metric files read from disk, myrank=",i3)
  deallocate(ttermij)
  return
1001 continue
  write(*,2001) itype,ir1,ir2,ir3,iodd,i,j
2001 format(" error in readit",7i5)
  call safestop(666)
  return
end subroutine readit

!============================================================
subroutine writeit
  use p_data
  use t_data
  use basic
  use arrays

  implicit none

  integer :: iodd, i, j, k, itri

      write(*,*) "before open 55 in writeup"
  open(55,file='setup',form='unformatted',status='unknown')
      write(*,*) "after open 55 in writeup"
  
  write(55) n,m,ioddm,iper,jper,alx,alz
  write(55) ((d2term(iodd,i),iodd=1,2),i=1,18)
  write(55) (((aterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
  write(55) (((dterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
  write(55) (((bterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
  write(55) (((sterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
  write(55) (((xterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
  write(55) (((yterm(iodd,i,j),iodd=1,2),i=1,18),j=1,18)
  
  write(55) ((((k0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((k1term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((k2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g4term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g5term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g6term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g7term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g8term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g9term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g10erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g11erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g12erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g13erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g14erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g15erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((g16erm(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((h3term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((h5term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((x0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((y0term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((x1term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((y1term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((x2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)
  write(55) ((((y2term(iodd,i,j,k),iodd=1,2),i=1,18),j=1,18),k=1,18)

  if(itor.ne.0) write(55)  (((hterm(itri,i,j),itri=1,ntridim),i=1,18),j=1,18)
  
  close(55)

  return
end subroutine writeit

!============================================================
subroutine wrrestart
  use p_data
  use t_data
  use basic
  use arrays

  implicit none

  integer :: ifirstrs, mmnn18, mmnn6, j1, igr, itime

  data ifirstrs/1/
  mmnn18 = m*n*numvar*6
!  if(ifirstrs .ne. 1) call rename('C1restartout','C1restartouto')
  if(ifirstrs .ne. 1) call system("mv C1restartout C1restartouto")

  ifirstrs = 0
  open(56,file='C1restartout',form='unformatted',status='unknown')
  write(56) (vel(j1),j1=1,mmnn18)
  write(56) (vel0(j1),j1=1,mmnn18)
  write(56) (vels(j1),j1=1,mmnn18)
  write(56) (phi(j1),j1=1,mmnn18)
  write(56) (phi0(j1),j1=1,mmnn18)
  write(56) (phis(j1),j1=1,mmnn18)
  write(56) (phiold(j1),j1=1,mmnn18)
  write(56) (veln(j1),j1=1,mmnn18)
  
  mmnn6 = 6*m*n
  write(56) (jphi(j1),j1=1,mmnn6)
  write(56) (vor(j1),j1=1,mmnn6)
  write(56) (com(j1),j1=1,mmnn6)
  write(56) (den(j1),j1=1,mmnn6)
  write(56) (den0(j1),j1=1,mmnn6)
  write(56) (pres(j1),j1=1,mmnn6)
  write(56) (pres0(j1),j1=1,mmnn6)
  write(56) (bsi(j1),j1=1,mmnn6)
  write(56) (vor0(j1),j1=1,mmnn6)
  write(56) (com0(j1),j1=1,mmnn6)
  write(56) (jphi0(j1),j1=1,mmnn6)
  
  write(56) ntime,time
  write(56) ((graphit(itime,igr),igr=1,maxplots),itime=0,ntime)
  write(56) ekin,emag,ekind,emagd,ekinp,emagp,ekinpd,emagpd,ekint,    &
       emagt,ekintd,emagtd,ekinph,ekinth,emagph,emagth,          &
       ekin3,ekin3d,ekin3h,emag3,emag3d,emag3h,tflux,gbound
  close(56)
  write(*,1201) ntime,time
  write(9,1201) ntime,time
1201 format("* * * restart file written at cycle,time = ",i5,1pe12.4,    &
          "  * * *")
  return
end subroutine wrrestart

!============================================================
subroutine rdrestart
  use p_data
  use t_data
  use basic
  use arrays

  implicit none
  integer :: mmnn18, mmnn6, j1, igr, itime

  mmnn18 = m*n*numvar*6
  open(56,file='restart',form='unformatted',status='old')
  read(56) (vel(j1),j1=1,mmnn18)
  read(56) (vel0(j1),j1=1,mmnn18)
  read(56) (vels(j1),j1=1,mmnn18)
  read(56) (phi(j1),j1=1,mmnn18)
  read(56) (phi0(j1),j1=1,mmnn18)
  read(56) (phis(j1),j1=1,mmnn18)
  read(56) (phiold(j1),j1=1,mmnn18)
  read(56) (veln(j1),j1=1,mmnn18)
!
!....DEBUG
!     if(myrank.le.0) then
!     do jj=1,m
!     write(98,1002) jj
!1002 format("jj =",i3)
!     write(98,1001) (vel(j1),j1=jj*18-5,jj*18)
!1001 format(1p6e12.4)
!     enddo
!     endif
!     call safestop(2)
!
!
  mmnn6 = 6*m*n
  read(56) (jphi(j1),j1=1,mmnn6)
  read(56) (vor(j1),j1=1,mmnn6)
  read(56) (com(j1),j1=1,mmnn6)
  read(56) (den(j1),j1=1,mmnn6)
  read(56) (den0(j1),j1=1,mmnn6)
  read(56) (pres(j1),j1=1,mmnn6)
  read(56) (pres0(j1),j1=1,mmnn6)
  read(56) (bsi(j1),j1=1,mmnn6)
  read(56) (vor0(j1),j1=1,mmnn6)
  read(56) (com0(j1),j1=1,mmnn6)
  read(56) (jphi0(j1),j1=1,mmnn6)

  read(56) ntimer,timer
  read(56) ((graphit(itime,igr),igr=1,maxplots),itime=0,ntimer)
  read(56) ekin,emag,ekind,emagd,ekinp,emagp,ekinpd,emagpd,ekint,   &
       emagt,ekintd,emagtd,ekinph,ekinth,emagph,emagth,        &
       ekin3,ekin3d,ekin3h,emag3,emag3d,emag3h,tflux,gbound
  close(56)
  write(*,1201) ntimer,timer
  write(9,1201) ntimer,timer
1201 format("* * * restart file read at cycle,time = ",i5,1pe12.4,     &
          "  * * *")
  return
end subroutine rdrestart
 subroutine netcdfout
 use netcdf

 use basic
 use arrays
!
 implicit none
!
 integer ncid   !  netCDF file ID   (defined by NF90_create)
 integer status ! error status return
 integer ifirst_netcdf ! flag for first call
 data ifirst_netcdf/0/
 integer timedim_id, rdim_id, zdim_id   ! dimension variables
 integer time_id ,r_id    !  independent variables 
 integer delta_ekin_id, delta_emag_id, n_etotd_id, n_emagd_id, delta_error_id,        &
         n_etoth_id, delta_emagt_id, delta_ekint_id, delta_ekin3_id, delta_emag3_id,  &
         ekin_id, emag_id, epress_id, etotal_id, psidif_id, psix_id, psidifdot_id,    &
         psixdot_id, tflux_id, chierror_id, totcur_id, vmax_id, totaltime_id,         &
         tread_id, plus_telements_id, plus_tsolve_id, tonestep_id, plus_tmpi_id,      &
         plus_tinit_id, plus_tadotb_id, maxrank_id, dt_id, thimp_id
 integer current_id, flux_id, ephi_id,                                  &
         jcb_id, vcb_id, etaj_id, hypr_id, dens_id, pe_id          !  1D + time variable id's
 integer jj
 real rarray(irs)
      character*60 string


 if(ifirst_netcdf == 0) then
  ifirst_netcdf = 1
!
! ...create the output file and define the file id ncid

  status = nf90_noerr
  if(irestart.ne.0) then
!
!....try and open existing file, otherwise open new file
     status = NF90_OPEN(netcdffilename,NF90_WRITE,ncid) 
  endif

  if(irestart.eq.0 .or. status.ne.nf90_noerr) then
    status = NF90_CREATE(netcdffilename,NF90_CLOBBER,ncid)
    if(status /= nf90_noerr) call handle_err(status,1)

!
! ...define dimensions that will be used
    status = NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED,timedim_id)
    if(status /= nf90_noerr) call handle_err(status,2)
    status = NF90_DEF_DIM(ncid,'R', irs, rdim_id)
    if(status /= nf90_noerr) call handle_err(status,3)
    status = NF90_DEF_DIM(ncid,'Z', irs, zdim_id)
    if(status /= nf90_noerr) call handle_err(status,4)
!
! ...define independent variables
    status = NF90_DEF_VAR(ncid,'time',NF90_DOUBLE,timedim_id,time_id)
    if(status /= nf90_noerr) call handle_err(status,5)
    status = NF90_DEF_VAR(ncid,'r',NF90_DOUBLE,rdim_id,r_id)
    if(status /= nf90_noerr) call handle_err(status,6)
!
!....define scalar functions of time variable only
    status = NF90_DEF_VAR(ncid,'delta_ekin',NF90_DOUBLE,timedim_id,delta_ekin_id)
    if(status /= nf90_noerr) call handle_err(status,1)
    status = NF90_DEF_VAR(ncid,'delta_emag',NF90_DOUBLE,timedim_id,delta_emag_id)
    if(status /= nf90_noerr) call handle_err(status,2)
    status = NF90_DEF_VAR(ncid,'N_etotd',NF90_DOUBLE,timedim_id,N_etotd_id)
    if(status /= nf90_noerr) call handle_err(status,17)  
    status = NF90_DEF_VAR(ncid,'N_emagd',NF90_DOUBLE,timedim_id,N_emagd_id)
    if(status /= nf90_noerr) call handle_err(status,7)
    status = NF90_DEF_VAR(ncid,'delta_error',NF90_DOUBLE,timedim_id,delta_error_id)
    if(status /= nf90_noerr) call handle_err(status,18)
    status = NF90_DEF_VAR(ncid,'N_etoth',NF90_DOUBLE,timedim_id,N_etoth_id)
    if(status /= nf90_noerr) call handle_err(status,21)
    status = NF90_DEF_VAR(ncid,'delta_emagt',NF90_DOUBLE,timedim_id,delta_emagt_id)
    if(status /= nf90_noerr) call handle_err(status,20)
    status = NF90_DEF_VAR(ncid,'delta_ekint',NF90_DOUBLE,timedim_id,delta_ekint_id)
    if(status /= nf90_noerr) call handle_err(status,19)
    status = NF90_DEF_VAR(ncid,'delta_ekin3',NF90_DOUBLE,timedim_id,delta_ekin3_id)
    if(status /= nf90_noerr) call handle_err(status,27)
    status = NF90_DEF_VAR(ncid,'delta_emag3',NF90_DOUBLE,timedim_id,delta_emag3_id)
    if(status /= nf90_noerr) call handle_err(status,28)
    status = NF90_DEF_VAR(ncid,'ekin',NF90_DOUBLE,timedim_id,ekin_id)
    if(status /= nf90_noerr) call handle_err(status,42)
    status = NF90_DEF_VAR(ncid,'emag',NF90_DOUBLE,timedim_id,emag_id)
    if(status /= nf90_noerr) call handle_err(status,43)
    status = NF90_DEF_VAR(ncid,'epress',NF90_DOUBLE,timedim_id,epress_id)
    if(status /= nf90_noerr) call handle_err(status,44)
    status = NF90_DEF_VAR(ncid,'etotal',NF90_DOUBLE,timedim_id,etotal_id)
    if(status /= nf90_noerr) call handle_err(status,45)  
    status = NF90_DEF_VAR(ncid,'psidif',NF90_DOUBLE,timedim_id,psidif_id)
    if(status /= nf90_noerr) call handle_err(status,25)
    status = NF90_DEF_VAR(ncid,'psix',NF90_DOUBLE,timedim_id,psix_id)
    if(status /= nf90_noerr) call handle_err(status,49)
    status = NF90_DEF_VAR(ncid,'psidifdot',NF90_DOUBLE,timedim_id,psidifdot_id)
    if(status /= nf90_noerr) call handle_err(status,26)
    status = NF90_DEF_VAR(ncid,'psixdot',NF90_DOUBLE,timedim_id,psixdot_id)
    if(status /= nf90_noerr) call handle_err(status,50)
    status = NF90_DEF_VAR(ncid,'tflux',NF90_DOUBLE,timedim_id,tflux_id)
    if(status /= nf90_noerr) call handle_err(status,29)
    status = NF90_DEF_VAR(ncid,'chierror',NF90_DOUBLE,timedim_id,chierror_id)
    if(status /= nf90_noerr) call handle_err(status,30)
    status = NF90_DEF_VAR(ncid,'totcur',NF90_DOUBLE,timedim_id,totcur_id)
    if(status /= nf90_noerr) call handle_err(status,31)
    status = NF90_DEF_VAR(ncid,'vmax',NF90_DOUBLE,timedim_id,vmax_id)
    if(status /= nf90_noerr) call handle_err(status,32)
    status = NF90_DEF_VAR(ncid,'totaltime',NF90_DOUBLE,timedim_id,totaltime_id)
    if(status /= nf90_noerr) call handle_err(status,33)
    status = NF90_DEF_VAR(ncid,'tread',NF90_DOUBLE,timedim_id,tread_id)
    if(status /= nf90_noerr) call handle_err(status,35)
    status = NF90_DEF_VAR(ncid,'plus_telements',NF90_DOUBLE,timedim_id,plus_telements_id)
    if(status /= nf90_noerr) call handle_err(status,36)
    status = NF90_DEF_VAR(ncid,'plus_tsolve',NF90_DOUBLE,timedim_id,plus_tsolve_id)
    if(status /= nf90_noerr) call handle_err(status,37)
    status = NF90_DEF_VAR(ncid,'tonestep',NF90_DOUBLE,timedim_id,tonestep_id)
    if(status /= nf90_noerr) call handle_err(status,38)
    status = NF90_DEF_VAR(ncid,'plus_tmpi',NF90_DOUBLE,timedim_id,plus_tmpi_id)
    if(status /= nf90_noerr) call handle_err(status,39)
    status = NF90_DEF_VAR(ncid,'plus_tinit',NF90_DOUBLE,timedim_id,plus_tinit_id)
    if(status /= nf90_noerr) call handle_err(status,41)
    status = NF90_DEF_VAR(ncid,'plus_tadotb',NF90_DOUBLE,timedim_id,plus_tadotb_id)
    if(status /= nf90_noerr) call handle_err(status,40)
    status = NF90_DEF_VAR(ncid,'maxrank',NF90_DOUBLE,timedim_id,maxrank_id)
    if(status /= nf90_noerr) call handle_err(status,34)
    status = NF90_DEF_VAR(ncid,'dt',NF90_DOUBLE,timedim_id,dt_id)
    if(status /= nf90_noerr) call handle_err(status,47)
    status = NF90_DEF_VAR(ncid,'thimp',NF90_DOUBLE,timedim_id,thimp_id)
    if(status /= nf90_noerr) call handle_err(status,48)
!
!
!
!....define functions of time and of R
    status = NF90_DEF_VAR(ncid,'current',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),current_id)
         if(status /= nf90_noerr) call handle_err(status,9)
    status = NF90_DEF_VAR(ncid,'pressure',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),pe_id)
         if(status /= nf90_noerr) call handle_err(status,10)
    status = NF90_DEF_VAR(ncid,'flux',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),flux_id)
         if(status /= nf90_noerr) call handle_err(status,11)
    status = NF90_DEF_VAR(ncid,'E_phi',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),ephi_id)
         if(status /= nf90_noerr) call handle_err(status,12)
    status = NF90_DEF_VAR(ncid,'JxB',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),jcb_id)
         if(status /= nf90_noerr) call handle_err(status,13)
    status = NF90_DEF_VAR(ncid,'VxB',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),vcb_id)
         if(status /= nf90_noerr) call handle_err(status,14)
    status = NF90_DEF_VAR(ncid,'etaxJ',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),etaj_id)
         if(status /= nf90_noerr) call handle_err(status,15)
    status = NF90_DEF_VAR(ncid,'hypr',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),hypr_id)
         if(status /= nf90_noerr) call handle_err(status,16)
    status = NF90_DEF_VAR(ncid,'density',NF90_DOUBLE,   &
                                     (/rdim_id,timedim_id /),dens_id)
         if(status /= nf90_noerr) call handle_err(status,17)
!
!....define dimensions of time
    status = NF90_PUT_ATT(ncid,time_id,'units','seconds')
         if(status /= nf90_noerr) call handle_err(status,18)
!
! ...define attribute for ElVis to start monitoring netCDF file
    status = NF90_PUT_ATT(ncid,NF90_GLOBAL,'running','true')
         if(status /= nf90_noerr) call handle_err(status,19)
!
!....end define mode
    status = NF90_ENDDEF(ncid)
         if(status /= nf90_noerr) call handle_err(status,20)
!
!
!....input data for r variable
    rarray = (/ ((alxp*(jj-1))/(irs-1), jj=1,irs) /)
         status = NF90_PUT_VAR(ncid,r_id,rarray)
         if(status /= nf90_noerr) call handle_err(status,21)
  else    ! on restart == 0
!
    status = NF90_REDEF(ncid)
         if(status /= nf90_noerr) call handle_err(status,220)
    status = NF90_PUT_ATT(ncid,NF90_GLOBAL,'running','true')
         if(status /= nf90_noerr) call handle_err(status,19)
    status = NF90_ENDDEF(ncid)
         if(status /= nf90_noerr) call handle_err(status,120)
!
!   nf90_inq_varid calls for independent variables
    status = NF90_INQ_VARID(ncid,'time',time_id)
         if(status /= nf90_noerr) call handle_err(status,5)
    status = NF90_INQ_VARID(ncid,'r',r_id)
         if(status /= nf90_noerr) call handle_err(status,6)
!
!...functions of time variable only
    status = NF90_INQ_VARID(ncid,'delta_ekin',delta_ekin_id)
         if(status /= nf90_noerr) call handle_err(status,1)
    status = NF90_INQ_VARID(ncid,'delta_emag',delta_emag_id)
         if(status /= nf90_noerr) call handle_err(status,2)
    status = NF90_INQ_VARID(ncid,'N_etotd',N_etotd_id)
         if(status /= nf90_noerr) call handle_err(status,17)  
    status = NF90_INQ_VARID(ncid,'N_emagd',N_emagd_id)
         if(status /= nf90_noerr) call handle_err(status,7)
    status = NF90_INQ_VARID(ncid,'delta_error',delta_error_id)
         if(status /= nf90_noerr) call handle_err(status,18)
    status = NF90_INQ_VARID(ncid,'N_etoth',N_etoth_id)
         if(status /= nf90_noerr) call handle_err(status,21)
    status = NF90_INQ_VARID(ncid,'delta_emagt',delta_emagt_id)
         if(status /= nf90_noerr) call handle_err(status,20)
    status = NF90_INQ_VARID(ncid,'delta_ekint',delta_ekint_id)
         if(status /= nf90_noerr) call handle_err(status,19)
    status = NF90_INQ_VARID(ncid,'delta_ekin3',delta_ekin3_id)
         if(status /= nf90_noerr) call handle_err(status,27)
    status = NF90_INQ_VARID(ncid,'delta_emag3',delta_emag3_id)
         if(status /= nf90_noerr) call handle_err(status,28)
    status = NF90_INQ_VARID(ncid,'ekin',ekin_id)
         if(status /= nf90_noerr) call handle_err(status,42)
    status = NF90_INQ_VARID(ncid,'emag',emag_id)
         if(status /= nf90_noerr) call handle_err(status,43)
    status = NF90_INQ_VARID(ncid,'epress',epress_id)
         if(status /= nf90_noerr) call handle_err(status,44)
    status = NF90_INQ_VARID(ncid,'etotal',etotal_id)
         if(status /= nf90_noerr) call handle_err(status,45)  
    status = NF90_INQ_VARID(ncid,'psidif',psidif_id)
         if(status /= nf90_noerr) call handle_err(status,25)
    status = NF90_INQ_VARID(ncid,'psix',psix_id)
         if(status /= nf90_noerr) call handle_err(status,49)
    status = NF90_INQ_VARID(ncid,'psidifdot',psidifdot_id)
         if(status /= nf90_noerr) call handle_err(status,26)
    status = NF90_INQ_VARID(ncid,'psixdot',psixdot_id)
         if(status /= nf90_noerr) call handle_err(status,50)
    status = NF90_INQ_VARID(ncid,'tflux',tflux_id)
         if(status /= nf90_noerr) call handle_err(status,29)
    status = NF90_INQ_VARID(ncid,'chierror',chierror_id)
         if(status /= nf90_noerr) call handle_err(status,30)
    status = NF90_INQ_VARID(ncid,'totcur',totcur_id)
         if(status /= nf90_noerr) call handle_err(status,31)
    status = NF90_INQ_VARID(ncid,'vmax',vmax_id)
         if(status /= nf90_noerr) call handle_err(status,32)
    status = NF90_INQ_VARID(ncid,'totaltime',totaltime_id)
         if(status /= nf90_noerr) call handle_err(status,33)
    status = NF90_INQ_VARID(ncid,'tread',tread_id)
         if(status /= nf90_noerr) call handle_err(status,35)
    status = NF90_INQ_VARID(ncid,'plus_telements',plus_telements_id)
         if(status /= nf90_noerr) call handle_err(status,36)
    status = NF90_INQ_VARID(ncid,'plus_tsolve',plus_tsolve_id)
         if(status /= nf90_noerr) call handle_err(status,37)
    status = NF90_INQ_VARID(ncid,'tonestep',tonestep_id)
         if(status /= nf90_noerr) call handle_err(status,38)
    status = NF90_INQ_VARID(ncid,'plus_tmpi',plus_tmpi_id)
         if(status /= nf90_noerr) call handle_err(status,39)
    status = NF90_INQ_VARID(ncid,'plus_tinit',plus_tinit_id)
         if(status /= nf90_noerr) call handle_err(status,41)
    status = NF90_INQ_VARID(ncid,'plus_tadotb',plus_tadotb_id)
         if(status /= nf90_noerr) call handle_err(status,40)
    status = NF90_INQ_VARID(ncid,'maxrank',maxrank_id)
         if(status /= nf90_noerr) call handle_err(status,34)
    status = NF90_INQ_VARID(ncid,'dt',dt_id)
         if(status /= nf90_noerr) call handle_err(status,47)
    status = NF90_INQ_VARID(ncid,'thimp',thimp_id)
         if(status /= nf90_noerr) call handle_err(status,48)
!
!...functions of time and of R
    status = NF90_INQ_VARID(ncid,'current',current_id)
         if(status /= nf90_noerr) call handle_err(status,9)
    status = NF90_INQ_VARID(ncid,'pressure',pe_id)
         if(status /= nf90_noerr) call handle_err(status,10)
    status = NF90_INQ_VARID(ncid,'flux',flux_id)
         if(status /= nf90_noerr) call handle_err(status,11)
    status = NF90_INQ_VARID(ncid,'E_phi',ephi_id)
         if(status /= nf90_noerr) call handle_err(status,12)
    status = NF90_INQ_VARID(ncid,'JxB',jcb_id)
         if(status /= nf90_noerr) call handle_err(status,13)
    status = NF90_INQ_VARID(ncid,'VxB',vcb_id)
         if(status /= nf90_noerr) call handle_err(status,14)
    status = NF90_INQ_VARID(ncid,'etaxJ',etaj_id)
         if(status /= nf90_noerr) call handle_err(status,15)
    status = NF90_INQ_VARID(ncid,'hypr',hypr_id)
         if(status /= nf90_noerr) call handle_err(status,16)
    status = NF90_INQ_VARID(ncid,'density',dens_id)
         if(status /= nf90_noerr) call handle_err(status,17)
!
  endif
 endif    ! on ifirst_netcdf == 0
!
!....file is now in data mode...input data for present timestep
!
!....insert present value of time
  status = NF90_PUT_VAR(ncid,time_id,time,(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,22)
!
!
!....insert scalar functions of time
!
  status = NF90_PUT_VAR(ncid,delta_ekin_id,graphit(ntime,1),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,1)
  status = NF90_PUT_VAR(ncid,delta_emag_id,graphit(ntime,2),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,2)
  status = NF90_PUT_VAR(ncid,n_etotd_id,graphit(ntime,17),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,17)
  status = NF90_PUT_VAR(ncid,n_emagd_id,graphit(ntime,7),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,7)
  status = NF90_PUT_VAR(ncid,delta_error_id,graphit(ntime,18),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,18)
  status = NF90_PUT_VAR(ncid,n_etoth_id,graphit(ntime,21),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,21)
  status = NF90_PUT_VAR(ncid,delta_emagt_id,graphit(ntime,20),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,20)
  status = NF90_PUT_VAR(ncid,delta_ekint_id,graphit(ntime,19),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,19)
  status = NF90_PUT_VAR(ncid,delta_ekin3_id,graphit(ntime,27),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,27)
  status = NF90_PUT_VAR(ncid,delta_emag3_id,graphit(ntime,28),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,28)
  status = NF90_PUT_VAR(ncid,ekin_id,graphit(ntime,42),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,42)
  status = NF90_PUT_VAR(ncid,emag_id,graphit(ntime,43),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,43)
  status = NF90_PUT_VAR(ncid,epress_id,graphit(ntime,44),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,44)
  status = NF90_PUT_VAR(ncid,etotal_id,graphit(ntime,45),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,45)
  status = NF90_PUT_VAR(ncid,psidif_id,graphit(ntime,25),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,25)
  status = NF90_PUT_VAR(ncid,psix_id,graphit(ntime,49),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,49)
  status = NF90_PUT_VAR(ncid,psidifdot_id,graphit(ntime,26),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,26)
  status = NF90_PUT_VAR(ncid,psixdot_id,graphit(ntime,50),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,50)
  status = NF90_PUT_VAR(ncid,tflux_id,graphit(ntime,29),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,29)
  status = NF90_PUT_VAR(ncid,chierror_id,graphit(ntime,30),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,30)
  status = NF90_PUT_VAR(ncid,totcur_id,graphit(ntime,31),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,31)
  status = NF90_PUT_VAR(ncid,vmax_id,graphit(ntime,32),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,32)
  status = NF90_PUT_VAR(ncid,totaltime_id,graphit(ntime,33),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,33)
  status = NF90_PUT_VAR(ncid,tread_id,graphit(ntime,35),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,35)
  status = NF90_PUT_VAR(ncid,plus_telements_id,graphit(ntime,36),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,36)
  status = NF90_PUT_VAR(ncid,plus_tsolve_id,graphit(ntime,37),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,37)
  status = NF90_PUT_VAR(ncid,tonestep_id,graphit(ntime,38),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,38)
  status = NF90_PUT_VAR(ncid,plus_tmpi_id,graphit(ntime,39),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,39)
  status = NF90_PUT_VAR(ncid,plus_tinit_id,graphit(ntime,41),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,41)
  status = NF90_PUT_VAR(ncid,plus_tadotb_id,graphit(ntime,40),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,40)
  status = NF90_PUT_VAR(ncid,maxrank_id,graphit(ntime,34),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,34)
  status = NF90_PUT_VAR(ncid,dt_id,graphit(ntime,47),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,47)
  status = NF90_PUT_VAR(ncid,thimp_id,graphit(ntime,48),(/ ntime /))
         if(status /= nf90_noerr) call handle_err(status,48)

!
!
!
!....insert functions of r and time
!
  status = NF90_PUT_VAR(ncid,current_id,curarray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,25)
  status = NF90_PUT_VAR(ncid,flux_id,psiarray,(/1,ntime/))
         if(status /= nf90_noerr) call handle_err(status,26)
  status = NF90_PUT_VAR(ncid,ephi_id,ephiarray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,27)
  status = NF90_PUT_VAR(ncid,jcb_id,jcbarray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,28)
  status = NF90_PUT_VAR(ncid,vcb_id,vcbarray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,29)
  status = NF90_PUT_VAR(ncid,etaj_id,etajarray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,30)
  status = NF90_PUT_VAR(ncid,hypr_id,hyprarray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,31)
  status = NF90_PUT_VAR(ncid,dens_id,densarray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,32)
  status = NF90_PUT_VAR(ncid,pe_id,pearray,(/1,ntime/),(/irs,1/))
         if(status /= nf90_noerr) call handle_err(status,33)
!
  status = nf90_sync(ncid)
         if(status /= nf90_noerr) call handle_err(status,34)
!
  if(ntime.ne.ntimemax) return
!
!...close file
!
  status = nf90_redef(ncid)
         if(status /= nf90_noerr) call handle_err(status,35)
!
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,'running','false')
         if(status /= nf90_noerr) call handle_err(status,36)

  status = nf90_enddef(ncid)
         if(status /= nf90_noerr) call handle_err(status,37)
!
  status = nf90_close(ncid)
         if(status /= nf90_noerr) call handle_err(status,38)
!
!.....copy netcdf file into local directory for saving
      string(1:3) = "cp "
      string(4:43) = netcdffilename
      string(44:60) = " ./C1monitor.cdf "
      call system(string)
! 
  end subroutine netcdfout
 subroutine handle_err(status,loc)
 use netcdf
  integer, intent (in) :: status,loc
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status)),loc
      call safestop(loc)
  endif
 end subroutine handle_err
end module inout_mod

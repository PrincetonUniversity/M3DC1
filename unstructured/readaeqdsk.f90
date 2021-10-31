module eqdsk_a

  implicit none

  character(len=1) :: dum
  character(len=66) :: lineaftertime
  logical :: mastformat=.false.
  integer, parameter, private :: magpri67=29, magpri322=31, magprirdp=8
!  integer, parameter, private :: magpri=magpri67+magpri322+magprirdp
  integer, parameter, private :: magpri=76
  integer, parameter, private :: nsilop=44, nesum=6, nfcoil=18
  integer, parameter, private :: nxtlim=9, nco2v=3, nco2r=2
  integer, parameter, private :: ntime=400
  integer, private :: nlold=35,nlnew=39

  integer, private :: nsilop0,magpri0,nfcoil0,nesum0
  integer, private :: lflag,mco2v,mco2r
  integer, private :: ishot,ktime1

  real, private :: pasman,betatn,psiq1,betat2,rexpan,fexpan,qqmin, &
       fexpvs,shearc,sepnose,ssi01,znose,rqqmin,rcencm,tavem,fluxx, &
       chipre,cprof,chigamt
  real, private, dimension(2,ntime) :: rseps,zseps

  character(len=10), private :: uday
  character(len=5), dimension(5) :: mfvers
  integer, private, dimension(ntime) :: jflag, jerror
  
  real, private, dimension(ntime) :: time, &
       eout,rout,zout,doutu, &
       doutl,aout,vout,betat,otop, &
       betap,ali,oleft,oright,qsta, &
       rcurrt,zcurrt,qout,olefs, &
       orighs,otops,sibdry,areao, &
       wplasm,elongm,qqmagx,terror, &
       rmagx,zmagx,obott,obots, &
       alpha,rttt,dbpli,delbp,oring, &
       sepexp,shearb, &
       xtch,ytch,qpsib,vertn,aaq1, &
       aaq2,aaq3,btaxp,btaxv, &
       simagx,seplim, &
       wbpol,taumhd,betapd,betatd, &
       alid,wplasmd,taudia,wbpold, &
       qmerci,slantu,slantl,zeff, &
       zeffr,tave,rvsin,zvsin, &
       rvsout,zvsout,wpdot,wbdot, &
       vsurfa,cjor95,pp95,ssep, &
       yyy2,xnnc, &
       wtherm,wfbeam,taujd3,tauthn, &
       qsiwant,cjorsw,cjor0, &
       ssiwant,ssi95,peak,dminux, &
       dminlx,dolubat,dolubafm,diludom, &
       diludomm,ratsol,rvsiu,zvsiu, &
       rvsid,zvsid,rvsou,zvsou, & 
       rvsod,zvsod,condno,psin32, &
       psin21,rq32in,rq21top,chilibt, &
       tsaisq,bcentr,cpasma,pasmat, &
       bpolav,s1,s2,s3,cdflux,psiref,xndnt,vloopt,pbinj, &
       zuperts,cjor99,psurfa,dolubaf, &
       cj1ave,rmidin,rmidout

  real, private, allocatable :: csilop(:,:)
  real, private, allocatable :: cmpr2(:,:)
  real, private, allocatable :: ccbrsp(:,:)
  real, private, allocatable :: eccurt(:,:)

  character(len=4), dimension(ntime) :: limloc

  real, private :: rco2r(nco2r,ntime),rco2v(nco2v,ntime),chordv(nco2v)
  real, private :: dco2r(ntime,nco2r),dco2v(ntime,nco2v)
contains

subroutine load_eqdsk_a(filename)

  implicit none

#ifdef USEMPI
  include 'mpif.h'
#endif

  character(len=*), intent(in) :: filename

  integer, parameter :: neqdsk = 12
  integer :: size, rank, ierr
  integer :: i, j, k, jj, ios
  real :: xdum

  character(len=3) :: qmflag

#ifdef USEMPI
  call MPI_comm_size(MPI_COMM_WORLD,size,ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
#else
  rank = 0
#endif 

  if(rank.eq.0) then
     write(0,*) 'Reading EQDSK a-file: ', trim(filename)
     open(unit=neqdsk,file=trim(filename),status='old')

     read(neqdsk,1055) uday,(mfvers(j),j=1,2)
     read(neqdsk,1053) ishot,ktime1
!    MAST a-files do not contain the following line in which the times are listed.
!    Use content of first line to identify MAST a-files.
     if(index(mfvers(1),'EFIT+')==0)then
        read(neqdsk,1040) (time(j),j=1,1)
     else
        mastformat=.true.
     end if

     write(0,*) uday, mfvers(1), mfvers(2)
     write(0,*) ' SHOT, time = ', ishot, ktime1

     do jj=1, ktime1
!       MAST a-files do not contain nlold,nlnew and thus format 1060 cannot be used
         read (neqdsk,'(A)') lineaftertime
         lineaftertime=trim(lineaftertime)//'     0    0'
         read (lineaftertime,1060) dum,time(jj),jflag(jj),lflag,limloc(jj), &
              mco2v,mco2r,qmflag,nlold,nlnew
        if(mco2v.gt.nco2v) then 
           write(0,*) 'Warning: mco2v > nco2v', mco2v
           mco2v = nco2v
        end if
        if(mco2r.gt.nco2r) then 
           write(0,*) 'Warning: mco2r > nco2r', mco2r
           mco2r = nco2r
        end if

        read (neqdsk,1040) tsaisq(jj),rcencm,bcentr(jj),pasmat(jj)
        read (neqdsk,1040) cpasma(jj),rout(jj),zout(jj),aout(jj)
        read (neqdsk,1040) eout(jj),doutu(jj),doutl(jj),vout(jj)
        read (neqdsk,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
        read (neqdsk,1040) betap(jj),ali(jj),oleft(jj),oright(jj)
        read (neqdsk,1040) otop(jj),obott(jj),qpsib(jj),vertn(jj)
        if(mco2v.gt.0) then
            read (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
            read (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
        end if
        if(mco2r.gt.0) then
            read (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
            read (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)
        end if
        read (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
        read (neqdsk,1040) s3(jj),qout(jj),olefs(jj),orighs(jj)
        read (neqdsk,1040) otops(jj),sibdry(jj),areao(jj),wplasm(jj)
        read (neqdsk,1040) terror(jj),elongm(jj),qqmagx(jj),cdflux(jj)
        read (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
        read (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj),zseps(2,jj)
        read (neqdsk,1040) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)
        read (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)
        read (neqdsk,1040) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
     
!        fluxx=diamag(jj)*1.0e-03
        read (neqdsk,1040) betapd(jj),betatd(jj),wplasmd(jj),fluxx
        read (neqdsk,1040) vloopt(jj),taudia(jj),qmerci(jj),tavem
!       In MAST a-files the following line uses e16.9 float formatting instead of i5
        if(mastformat)then
            read (neqdsk,*) nsilop0,magpri0,nfcoil0,nesum0
        else
            read (neqdsk,1041) nsilop0,magpri0,nfcoil0,nesum0
        end if
        
        write(0,*) 'nfcoil0 = ', nfcoil0

        if(nsilop0.gt.nsilop) then 
           write(0,*) 'Warning: nsilop0 > nsilop', nsilop0
!           nsilop0 = nsilop
        end if
        if(magpri0.gt.magpri) then 
           write(0,*) 'Warning: magpri0 > magpri', magpri0
!           magpri0 = magpri
        end if
        if(nfcoil0.gt.nfcoil) then 
           write(0,*) 'Warning: nfcoil0 > nfcoil', nfcoil0
!           nfcoil0 = nfcoil
        end if
        if(nesum0.gt.nesum) then 
           write(0,*) 'Warning: nesum0 > nesum', nesum0
!           nesum0 = nesum
        end if
        allocate(csilop(nsilop0,ntime))
        allocate(cmpr2(magpri0,ntime))
        allocate(ccbrsp(nfcoil0,ntime))
        allocate(eccurt(nesum0,ntime))

        read (neqdsk,1040,end=98) (csilop(k,jj),k=1,nsilop0), (cmpr2(k,jj),k=1,magpri0)
        read (neqdsk,1040,end=98) (ccbrsp(k,jj),k=1,nfcoil0)
!       Do not read the following in case of MAST
        if(.not. mastformat)then
            read (neqdsk,1040,end=99) (eccurt(jj,k),k=1,nesum0)
            read (neqdsk,1040,end=99) pbinj(jj),rvsin(jj),zvsin(jj),rvsout(jj)
            read (neqdsk,1040,end=99) zvsout(jj),vsurfa(jj),wpdot(jj),wbdot(jj)
            read (neqdsk,1040,end=99) slantu(jj),slantl(jj),zuperts(jj),chipre
            read (neqdsk,1040,end=99) cjor95(jj),pp95(jj),ssep(jj),yyy2(jj)
            read (neqdsk,1040,end=99) xnnc(jj),cprof,oring(jj),cjor0(jj)
            read (neqdsk,1040,end=99) fexpan,qqmin,chigamt,ssi01
            read (neqdsk,1040,end=99) fexpvs,sepnose,ssi95(jj),rqqmin
            read (neqdsk,1040,end=99) cjor99(jj),cj1ave(jj),rmidin(jj),rmidout(jj)
            read (neqdsk,1040,end=99) psurfa(jj), peak(jj),dminux(jj),dminlx(jj)
            
            read (neqdsk,1040,end=99) dolubaf(jj),dolubafm(jj),diludom(jj),diludomm(jj)
            read (neqdsk,1040,end=99) ratsol(jj),rvsiu(jj),zvsiu(jj),rvsid(jj)
            read (neqdsk,1040,end=99) zvsid(jj),rvsou(jj),zvsou(jj),rvsod(jj)
            read (neqdsk,1040,end=99) zvsod(jj),condno(jj),psin32(jj),psin21(jj)
            read (neqdsk,1040,end=99) rq32in(jj),rq21top(jj),chilibt(jj),xdum
        end if
     end do

     close(neqdsk)
  end if

  goto 100
98 write(0,*) 'Error: encountered unexpected EOF (before coil currents)'
  return
99 write(0,*) 'Warning: encountered unexpected EOF (after coil currents)'
100 continue

1040 format (1x,4e16.9)
1041 format (1x,4i5)
1050 format (1x,i5,11x,i5)
1053 format (1x,i6,11x,i5)
1055 format (1x,a10,2a5)
1060 format (1a1,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3,1x,2i5)

  if(nfcoil0.eq.18) then 
     write(0,*) 'Assuming DIII-D'
     if(maxval(abs(ccbrsp(1:18,1))).lt.20000.) then
        write(0,*) 'File data is in Amps'
        write(*,1100) ccbrsp( 4,1)/1000.*58.
        write(*,1100) ccbrsp( 3,1)/1000.*58.
        write(*,1100) ccbrsp( 2,1)/1000.*58.
        write(*,1100) ccbrsp( 1,1)/1000.*58.
        write(*,1100) ccbrsp(10,1)/1000.*58.
        write(*,1100) ccbrsp(11,1)/1000.*58.
        write(*,1100) ccbrsp(12,1)/1000.*58.
        write(*,1100) ccbrsp(13,1)/1000.*55.
        write(*,1100) ccbrsp( 5,1)/1000.*58.
        write(*,1100) ccbrsp(14,1)/1000.*55.
        write(*,1100) ccbrsp( 8,1)/1000.*58.
        write(*,1100) ccbrsp(17,1)/1000.*55.
        write(*,1100) ccbrsp( 9,1)/1000.*58.
        write(*,1100) ccbrsp(18,1)/1000.*55.
        write(*,1100) ccbrsp( 7,1)/1000.*58.
        write(*,1100) ccbrsp(16,1)/1000.*55.
        write(*,1100) ccbrsp( 6,1)/1000.*58.
        write(*,1100) ccbrsp(15,1)/1000.*55.
     else
        write(0,*) 'File data is in Amp-turns'
        write(*,1100) ccbrsp( 4,1)/1000.
        write(*,1100) ccbrsp( 3,1)/1000.
        write(*,1100) ccbrsp( 2,1)/1000.
        write(*,1100) ccbrsp( 1,1)/1000.
        write(*,1100) ccbrsp(10,1)/1000.
        write(*,1100) ccbrsp(11,1)/1000.
        write(*,1100) ccbrsp(12,1)/1000.
        write(*,1100) ccbrsp(13,1)/1000.
        write(*,1100) ccbrsp( 5,1)/1000.
        write(*,1100) ccbrsp(14,1)/1000.
        write(*,1100) ccbrsp( 8,1)/1000.
        write(*,1100) ccbrsp(17,1)/1000.
        write(*,1100) ccbrsp( 9,1)/1000.
        write(*,1100) ccbrsp(18,1)/1000.
        write(*,1100) ccbrsp( 7,1)/1000.
        write(*,1100) ccbrsp(16,1)/1000.
        write(*,1100) ccbrsp( 6,1)/1000.
        write(*,1100) ccbrsp(15,1)/1000.
     end if

  else if(nfcoil0.eq.12) then 
     write(0,*) 'Assuming EAST'
     write(*,1100) ccbrsp( 1,1)/1000.  ! PF1
     write(*,1100) ccbrsp( 7,1)/1000.  ! PF2

     write(*,1100) ccbrsp( 2,1)/1000.  ! PF3
     write(*,1100) ccbrsp( 8,1)/1000.  ! PF4

     write(*,1100) ccbrsp( 3,1)/1000.  ! PF5
     write(*,1100) ccbrsp( 9,1)/1000.  ! PF6

     write(*,1100) ccbrsp( 4,1)*(44./248.)/1000.  ! PF7
     write(*,1100) ccbrsp( 4,1)*(204./248.)/1000.  ! PF9

     write(*,1100) ccbrsp(10,1)*(44./248.)/1000.  ! PF8
     write(*,1100) ccbrsp(10,1)*(204./248.)/1000.  ! PF10

     write(*,1100) ccbrsp( 5,1)/1000.  ! PF11
     write(*,1100) ccbrsp(11,1)/1000.  ! PF12

     write(*,1100) ccbrsp( 6,1)/1000.  ! PF13
     write(*,1100) ccbrsp(12,1)/1000.  ! PF14

  else if(nfcoil0.eq.52) then
     write(0,*) 'Assuming NSTX'
     write(0,*) 'nesum0 = ', nesum0
     write(*,1100) eccurt(1,1)*122.625/1000.  ! OH1U
     write(*,1100) eccurt(1,1)*121.000/1000.  ! OH2U
     write(*,1100) eccurt(1,1)*119.875/1000.  ! OH3U
     write(*,1100) eccurt(1,1)*118.000/1000.  ! OH4U
     write(*,1100) eccurt(1,1)*118.000/1000.  ! OH4L
     write(*,1100) eccurt(1,1)*119.875/1000.  ! OH3L
     write(*,1100) eccurt(1,1)*121.000/1000.  ! OH2L
     write(*,1100) eccurt(1,1)*122.625/1000.  ! OH1L
     write(*,1100) ccbrsp( 1,1)*20/1000.  ! PF1AU
     write(*,1100) ccbrsp( 2,1)*14/1000.  ! PF2U1
     write(*,1100) ccbrsp( 2,1)*14/1000.  ! PF2U2
     write(*,1100) ccbrsp( 3,1)*15/1000.  ! PF3U1
     write(*,1100) ccbrsp( 3,1)*15/1000.  ! PF3U2
     write(*,1100) ccbrsp( 4,1)*9 /1000.  ! PF4U1
     write(*,1100) ccbrsp( 4,1)*8 /1000.  ! PF4U2
     write(*,1100) ccbrsp( 5,1)*12/1000.  ! PF5U1
     write(*,1100) ccbrsp( 5,1)*12/1000.  ! PF5U2
     write(*,1100) ccbrsp( 6,1)*12/1000.  ! PF5L1
     write(*,1100) ccbrsp( 6,1)*12/1000.  ! PF5L2
     write(*,1100) ccbrsp( 7,1)*9 /1000.  ! PF4L1
     write(*,1100) ccbrsp( 7,1)*8 /1000.  ! PF4L2
     write(*,1100) ccbrsp( 8,1)*15/1000.  ! PF3L1
     write(*,1100) ccbrsp( 8,1)*15/1000.  ! PF3L2
     write(*,1100) ccbrsp( 9,1)*14/1000.  ! PF2L1
     write(*,1100) ccbrsp( 9,1)*14/1000.  ! PF2L2
     write(*,1100) ccbrsp(10,1)*20/1000.  ! PF1AL
     write(*,1100) ccbrsp(11,1)*32/1000.  ! PF1B
     write(*,1100) ccbrsp(12,1)* 1/1000.  ! TFUI
     write(*,1100) ccbrsp(13,1)* 1/1000.  ! TFUO
     write(*,1100) ccbrsp(15,1)* 1/1000.  ! TFLO
     write(*,1100) ccbrsp(14,1)* 1/1000.  ! TFLI
     write(*,1100) ccbrsp(16,1)*48/1000.  ! PFAB1
     write(*,1100) ccbrsp(17,1)*48/1000.  ! PFAB2
     write(*,1100) ccbrsp(18,1)* 1/1000.  ! VS1U
     write(*,1100) ccbrsp(19,1)* 1/1000.  ! VS2U
     write(*,1100) ccbrsp(20,1)* 1/1000.  ! VS3U
     write(*,1100) ccbrsp(21,1)*.5/1000.  ! VS4U
     write(*,1100) ccbrsp(21,1)*.5/1000.  ! VS4U
     write(*,1100) ccbrsp(22,1)* 1/1000.  ! VS5U
     write(*,1100) ccbrsp(23,1)*.5/1000.  ! VS6U
     write(*,1100) ccbrsp(23,1)*.5/1000.  ! VS6U
     write(*,1100) ccbrsp(24,1)*.3333/1000.  ! VS7U
     write(*,1100) ccbrsp(24,1)*.3333/1000.  ! VS7U
     write(*,1100) ccbrsp(24,1)*.3333/1000.  ! VS7U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(25,1)*.1250/1000.  ! VS8U
     write(*,1100) ccbrsp(26,1)*.2/1000.  ! VS9U
     write(*,1100) ccbrsp(26,1)*.2/1000.  ! VS9U
     write(*,1100) ccbrsp(26,1)*.2/1000.  ! VS9U
     write(*,1100) ccbrsp(26,1)*.2/1000.  ! VS9U
     write(*,1100) ccbrsp(26,1)*.2/1000.  ! VS9U
     write(*,1100) ccbrsp(27,1)*.5/1000.  ! VS10U
     write(*,1100) ccbrsp(27,1)*.5/1000.  ! VS10U
     write(*,1100) ccbrsp(28,1)*.3333/1000.  ! VS11U
     write(*,1100) ccbrsp(28,1)*.3333/1000.  ! VS11U
     write(*,1100) ccbrsp(28,1)*.3333/1000.  ! VS11U
     write(*,1100) ccbrsp(29,1)* 1/1000.  ! VS12U
     write(*,1100) ccbrsp(30,1)* 1/1000.  ! VS13U
     write(*,1100) ccbrsp(31,1)* 1/1000.  ! VS13L
     write(*,1100) ccbrsp(32,1)* 1/1000.  ! VS12L
     write(*,1100) ccbrsp(33,1)*.3333/1000.  ! VS11L
     write(*,1100) ccbrsp(33,1)*.3333/1000.  ! VS11L
     write(*,1100) ccbrsp(33,1)*.3333/1000.  ! VS11L
     write(*,1100) ccbrsp(34,1)*.5/1000.  ! VS10L
     write(*,1100) ccbrsp(34,1)*.5/1000.  ! VS10L
     write(*,1100) ccbrsp(35,1)*.2/1000.  ! VS9L
     write(*,1100) ccbrsp(35,1)*.2/1000.  ! VS9L
     write(*,1100) ccbrsp(35,1)*.2/1000.  ! VS9L
     write(*,1100) ccbrsp(35,1)*.2/1000.  ! VS9L
     write(*,1100) ccbrsp(35,1)*.2/1000.  ! VS9L
     write(*,1100) ccbrsp(36,1)*.2/1000.  ! VS8L
     write(*,1100) ccbrsp(36,1)*.2/1000.  ! VS8L
     write(*,1100) ccbrsp(36,1)*.2/1000.  ! VS8L
     write(*,1100) ccbrsp(36,1)*.2/1000.  ! VS8L
     write(*,1100) ccbrsp(36,1)*.2/1000.  ! VS8L
     write(*,1100) ccbrsp(37,1)*.5/1000.  ! VS6L
     write(*,1100) ccbrsp(37,1)*.5/1000.  ! VS6L
     write(*,1100) ccbrsp(38,1)* 1/1000.  ! VS5L
     write(*,1100) ccbrsp(39,1)*.5/1000.  ! VS4L
     write(*,1100) ccbrsp(39,1)*.5/1000.  ! VS4L
     write(*,1100) ccbrsp(40,1)* 1/1000.  ! VS3L
     write(*,1100) ccbrsp(41,1)* 1/1000.  ! VS2L
     write(*,1100) ccbrsp(42,1)* 1/1000.  ! VS1L
     write(*,1100) ccbrsp(43,1)*.5/1000.  ! DPU1
     write(*,1100) ccbrsp(43,1)*.5/1000.  ! DPU1
     write(*,1100) ccbrsp(44,1)*.5/1000.  ! DPL1
     write(*,1100) ccbrsp(44,1)*.5/1000.  ! DPL1
     write(*,1100) ccbrsp(45,1)* 1/1000.  ! PPSIUU
     write(*,1100) ccbrsp(46,1)* 1/1000.  ! PPSIUL
     write(*,1100) ccbrsp(47,1)* 1/1000.  ! PPPOUU
     write(*,1100) ccbrsp(48,1)* 1/1000.  ! PPPOUL
     write(*,1100) ccbrsp(49,1)* 1/1000.  ! PPPOLU
     write(*,1100) ccbrsp(50,1)* 1/1000.  ! PPPOLL
     write(*,1100) ccbrsp(51,1)* 1/1000.  ! PPSILU
     write(*,1100) ccbrsp(52,1)* 1/1000.  ! PPSILL

  else if(nfcoil0.eq.54) then
     write(0,*) 'Assuming NSTX-U'
     write(*,1100) 0.                  ! OH1U
     write(*,1100) 0.                  ! OH2U
     write(*,1100) 0.                  ! OH3U
     write(*,1100) 0.                  ! OH4U
     write(*,1100) 0.                  ! OH4L
     write(*,1100) 0.                  ! OH3L
     write(*,1100) 0.                  ! OH2L
     write(*,1100) 0.                  ! OH1L
     write(*,1100) ccbrsp( 1,1)*64/1000.  ! PF1AU
     write(*,1100) ccbrsp( 2,1)*32/1000.  ! PF1BU
     write(*,1100) ccbrsp( 3,1)*20/1000.  ! PF1CU
     write(*,1100) ccbrsp( 4,1)*14/1000.  ! PF2U1
     write(*,1100) ccbrsp( 4,1)*14/1000.  ! PF2U2
     write(*,1100) ccbrsp( 5,1)*15/1000.  ! PF3U1
     write(*,1100) ccbrsp( 5,1)*15/1000.  ! PF3U2
     write(*,1100) ccbrsp( 6,1)*9 /1000.  ! PF4U1
     write(*,1100) ccbrsp( 6,1)*8 /1000.  ! PF4U2
     write(*,1100) ccbrsp( 7,1)*12/1000.  ! PF5U1
     write(*,1100) ccbrsp( 7,1)*12/1000.  ! PF5U2
     write(*,1100) ccbrsp( 8,1)*12/1000.  ! PF5L1
     write(*,1100) ccbrsp( 8,1)*12/1000.  ! PF5L2
     write(*,1100) ccbrsp( 9,1)*9 /1000.  ! PF4L1
     write(*,1100) ccbrsp( 9,1)*8 /1000.  ! PF4L2
     write(*,1100) ccbrsp(10,1)*15/1000.  ! PF3L1
     write(*,1100) ccbrsp(10,1)*15/1000.  ! PF3L2
     write(*,1100) ccbrsp(11,1)*14/1000.  ! PF2L1
     write(*,1100) ccbrsp(11,1)*14/1000.  ! PF2L2
     write(*,1100) ccbrsp(12,1)*20/1000.  ! PF1CL
     write(*,1100) ccbrsp(13,1)*32/1000.  ! PF1BL
     write(*,1100) ccbrsp(14,1)*64/1000.  ! PF1AL
     write(*,1100) ccbrsp(15,1)*.5/1000.  ! VS1U
     write(*,1100) ccbrsp(15,1)*.5/1000.  ! VS1U
     write(*,1100) ccbrsp(16,1)*.5/1000.  ! VS2U
     write(*,1100) ccbrsp(16,1)*.5/1000.  ! VS2U
     write(*,1100) ccbrsp(17,1)* 1/1000.  ! VS3U
     write(*,1100) ccbrsp(18,1)* 1/1000.  ! VS4U
     write(*,1100) ccbrsp(19,1)*.3333/1000.  ! VS5U
     write(*,1100) ccbrsp(19,1)*.3333/1000.  ! VS5U
     write(*,1100) ccbrsp(19,1)*.3333/1000.  ! VS5U
     do i=1, 16
        write(*,1100) ccbrsp(20,1)*.0625/1000.  ! VS6U
     end do
     do i=1, 5
        write(*,1100) ccbrsp(21,1)*.2000/1000.  ! VS7U
     end do
     do i=1, 4
        write(*,1100) ccbrsp(22,1)*.2500/1000.  ! VS8U
     end do
     do i=1, 24
        write(*,1100) ccbrsp(23,1)*.041667/1000.  ! VS9U
     end do
     write(*,1100) ccbrsp(24,1)*.5/1000.  ! VS10U
     write(*,1100) ccbrsp(24,1)*.5/1000.  ! VS10U
     do i=1, 5
        write(*,1100) ccbrsp(25,1)*.2/1000.  ! VS11U
     end do
     write(*,1100) ccbrsp(26,1)*.5/1000.  ! VS12U
     write(*,1100) ccbrsp(26,1)*.5/1000.  ! VS12U
     write(*,1100) ccbrsp(27,1)*.3333/1000.  ! VS13U
     write(*,1100) ccbrsp(27,1)*.3333/1000.  ! VS13U
     write(*,1100) ccbrsp(27,1)*.3333/1000.  ! VS13U
     write(*,1100) ccbrsp(28,1)*1/1000.   ! VS14U
     write(*,1100) ccbrsp(29,1)*1/1000.   ! VS15U
     write(*,1100) ccbrsp(30,1)*1/1000.   ! VS15L
     write(*,1100) ccbrsp(31,1)*1/1000.   ! VS14L
     write(*,1100) ccbrsp(32,1)*.3333/1000.  ! VS13L
     write(*,1100) ccbrsp(32,1)*.3333/1000.  ! VS13L
     write(*,1100) ccbrsp(32,1)*.3333/1000.  ! VS13L
     write(*,1100) ccbrsp(33,1)*.5/1000.  ! VS12L
     write(*,1100) ccbrsp(33,1)*.5/1000.  ! VS12L
     do i=1, 5
        write(*,1100) ccbrsp(34,1)*.2/1000.  ! VS11L
     end do
     write(*,1100) ccbrsp(35,1)*.5/1000.  ! VS10L
     write(*,1100) ccbrsp(35,1)*.5/1000.  ! VS10L
     do i=1, 24
        write(*,1100) ccbrsp(36,1)*.041667/1000.  ! VS9L
     end do
     do i=1, 16
        write(*,1100) ccbrsp(37,1)*.0625/1000.  ! VS8L
     end do
     do i=1, 4
        write(*,1100) ccbrsp(38,1)*.25/1000.  ! VS7L
     end do
     do i=1, 5
        write(*,1100) ccbrsp(39,1)*.2/1000.   ! VS6L
     end do
     write(*,1100) ccbrsp(40,1)*.3333/1000.  ! VS5L
     write(*,1100) ccbrsp(40,1)*.3333/1000.  ! VS5L
     write(*,1100) ccbrsp(40,1)*.3333/1000.  ! VS5L
     write(*,1100) ccbrsp(41,1)* 1/1000.  ! VS4L
     write(*,1100) ccbrsp(42,1)* 1/1000.  ! VS3L
     write(*,1100) ccbrsp(43,1)*.5/1000.  ! VS2L
     write(*,1100) ccbrsp(43,1)*.5/1000.  ! VS2L
     write(*,1100) ccbrsp(44,1)*.5/1000.  ! VS1L
     write(*,1100) ccbrsp(44,1)*.5/1000.  ! VS1L

     write(*,1100) ccbrsp(45,1)*.5/1000.  ! DPU1
     write(*,1100) ccbrsp(45,1)*.5/1000.  ! DPU1
     write(*,1100) ccbrsp(46,1)*.5/1000.  ! DPL1
     write(*,1100) ccbrsp(46,1)*.5/1000.  ! DPL1
     write(*,1100) ccbrsp(47,1)* 1/1000.  ! PPSIUU
     write(*,1100) ccbrsp(48,1)* 1/1000.  ! PPSIUL
     write(*,1100) ccbrsp(49,1)* 1/1000.  ! PPPOUU
     write(*,1100) ccbrsp(50,1)* 1/1000.  ! PPPOUL
     write(*,1100) ccbrsp(51,1)* 1/1000.  ! PPPOLU
     write(*,1100) ccbrsp(52,1)* 1/1000.  ! PPPOLL
     write(*,1100) ccbrsp(53,1)* 1/1000.  ! PPSILU
     write(*,1100) ccbrsp(54,1)* 1/1000.  ! PPSILL
  else if(nfcoil0.eq.31) then
     write(0,*) 'Assuming KSTAR '
!!$     write(*,1100) ccbrsp(27,1)*180./1000.  ! PF1U
!!$     write(*,1100) ccbrsp(28,1)*144./1000.  ! PF2U
!!$     write(*,1100) ccbrsp(19,1)* 72./1000.  ! PF3U
!!$     write(*,1100) ccbrsp(20,1)*108./1000.  ! PF4U
!!$     write(*,1100) ccbrsp(21,1)*208./1000.  ! PF5U
!!$     write(*,1100) ccbrsp(22,1)*128./1000.  ! PF6U
!!$     write(*,1100) ccbrsp(29,1)* 72./1000.  ! PF7U
!!$     write(*,1100) ccbrsp(29,1)* 72./1000.  ! PF7L
!!$     write(*,1100) ccbrsp(26,1)*128./1000.  ! PF6L
!!$     write(*,1100) ccbrsp(25,1)*208./1000.  ! PF5L
!!$     write(*,1100) ccbrsp(24,1)*108./1000.  ! PF4L
!!$     write(*,1100) ccbrsp(23,1)* 72./1000.  ! PF3L
!!$     write(*,1100) ccbrsp(28,1)*144./1000.  ! PF2L
!!$     write(*,1100) ccbrsp(27,1)*180./1000.  ! PF1L
!!$
!!$     write(*,1100) ccbrsp( 1,1)*0.2500/1000.  ! VOCS1U
!!$     write(*,1100) ccbrsp( 1,1)*0.2500/1000.  ! VOCS2U
!!$     write(*,1100) ccbrsp( 2,1)*0.1667/1000.  ! VOCS3U
!!$     write(*,1100) ccbrsp( 2,1)*0.1667/1000.  ! VOD4U
!!$     write(*,1100) ccbrsp( 2,1)*0.1667/1000.  ! VOD5U
!!$     write(*,1100) ccbrsp( 3,1)*0.2500/1000.  ! VOD6U
!!$     write(*,1100) ccbrsp( 3,1)*0.2500/1000.  ! VOD7U
!!$     write(*,1100) ccbrsp( 4,1)*0.2500/1000.  ! VOD8U
!!$     write(*,1100) ccbrsp( 4,1)*0.2500/1000.  ! VOD9U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VOD10U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VOD11U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VOD12U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VOD13U
!!$     write(*,1100) ccbrsp( 6,1)*0.2500/1000.  ! VOD14U
!!$     write(*,1100) ccbrsp( 6,1)*0.2500/1000.  ! VOD15U
!!$     write(*,1100) ccbrsp( 7,1)*0.2500/1000.  ! VOD15L
!!$     write(*,1100) ccbrsp( 7,1)*0.2500/1000.  ! VOD14L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VOD13L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VOD12L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VOD11L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VOD10L
!!$     write(*,1100) ccbrsp( 9,1)*0.2500/1000.  ! VOD9L
!!$     write(*,1100) ccbrsp( 9,1)*0.2500/1000.  ! VOD8L
!!$     write(*,1100) ccbrsp(10,1)*0.2500/1000.  ! VOD7L
!!$     write(*,1100) ccbrsp(10,1)*0.2500/1000.  ! VOD6L
!!$     write(*,1100) ccbrsp(11,1)*0.1667/1000.  ! VOD5L
!!$     write(*,1100) ccbrsp(11,1)*0.1667/1000.  ! VOD4L
!!$     write(*,1100) ccbrsp(11,1)*0.1667/1000.  ! VOD3L
!!$     write(*,1100) ccbrsp(12,1)*0.2500/1000.  ! VOD2L
!!$     write(*,1100) ccbrsp(12,1)*0.2500/1000.  ! VOD1L
!!$
!!$     write(*,1100) ccbrsp( 1,1)*0.2500/1000.  ! VICS1U
!!$     write(*,1100) ccbrsp( 1,1)*0.2500/1000.  ! VICS2U
!!$     write(*,1100) ccbrsp( 2,1)*0.1667/1000.  ! VICS3U
!!$     write(*,1100) ccbrsp( 2,1)*0.1667/1000.  ! VID4U
!!$     write(*,1100) ccbrsp( 2,1)*0.1667/1000.  ! VID5U
!!$     write(*,1100) ccbrsp( 3,1)*0.2500/1000.  ! VID6U
!!$     write(*,1100) ccbrsp( 3,1)*0.2500/1000.  ! VID7U
!!$     write(*,1100) ccbrsp( 4,1)*0.2500/1000.  ! VID8U
!!$     write(*,1100) ccbrsp( 4,1)*0.2500/1000.  ! VID9U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VID10U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VID11U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VID12U
!!$     write(*,1100) ccbrsp( 5,1)*0.1250/1000.  ! VID13U
!!$     write(*,1100) ccbrsp( 6,1)*0.2500/1000.  ! VIW14U
!!$     write(*,1100) ccbrsp( 6,1)*0.2500/1000.  ! VIW15U
!!$     write(*,1100) ccbrsp( 7,1)*0.2500/1000.  ! VIW15L
!!$     write(*,1100) ccbrsp( 7,1)*0.2500/1000.  ! VIW14L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VID13L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VID12L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VID11L
!!$     write(*,1100) ccbrsp( 8,1)*0.1250/1000.  ! VID10L
!!$     write(*,1100) ccbrsp( 9,1)*0.2500/1000.  ! VID9L
!!$     write(*,1100) ccbrsp( 9,1)*0.2500/1000.  ! VID8L
!!$     write(*,1100) ccbrsp(10,1)*0.2500/1000.  ! VID7L
!!$     write(*,1100) ccbrsp(10,1)*0.2500/1000.  ! VID6L
!!$     write(*,1100) ccbrsp(11,1)*0.1667/1000.  ! VID5L
!!$     write(*,1100) ccbrsp(11,1)*0.1667/1000.  ! VID4L
!!$     write(*,1100) ccbrsp(11,1)*0.1667/1000.  ! VICS3L
!!$     write(*,1100) ccbrsp(12,1)*0.2500/1000.  ! VICS2L
!!$     write(*,1100) ccbrsp(12,1)*0.2500/1000.  ! VICS1L

     write(*,1100) ccbrsp(27,1)/1000.  ! PF1U
     write(*,1100) ccbrsp(28,1)/1000.  ! PF2U
     write(*,1100) ccbrsp(19,1)/1000.  ! PF3U
     write(*,1100) ccbrsp(20,1)/1000.  ! PF4U
     write(*,1100) ccbrsp(21,1)/1000.  ! PF5U
     write(*,1100) ccbrsp(22,1)/1000.  ! PF6U
     write(*,1100) ccbrsp(29,1)/1000.  ! PF7U
     write(*,1100) ccbrsp(29,1)/1000.  ! PF7L
     write(*,1100) ccbrsp(26,1)/1000.  ! PF6L
     write(*,1100) ccbrsp(25,1)/1000.  ! PF5L
     write(*,1100) ccbrsp(24,1)/1000.  ! PF4L
     write(*,1100) ccbrsp(23,1)/1000.  ! PF3L
     write(*,1100) ccbrsp(28,1)/1000.  ! PF2L
     write(*,1100) ccbrsp(27,1)/1000.  ! PF1L

     write(*,1100) ccbrsp( 1,1)/1000.  ! VOCS1U
     write(*,1100) ccbrsp( 1,1)/1000.  ! VOCS2U
     write(*,1100) ccbrsp( 2,1)/1000.  ! VOCS3U
     write(*,1100) ccbrsp( 2,1)/1000.  ! VOD4U
     write(*,1100) ccbrsp( 2,1)/1000.  ! VOD5U
     write(*,1100) ccbrsp( 3,1)/1000.  ! VOD6U
     write(*,1100) ccbrsp( 3,1)/1000.  ! VOD7U
     write(*,1100) ccbrsp( 4,1)/1000.  ! VOD8U
     write(*,1100) ccbrsp( 4,1)/1000.  ! VOD9U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VOD10U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VOD11U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VOD12U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VOD13U
     write(*,1100) ccbrsp( 6,1)/1000.  ! VOD14U
     write(*,1100) ccbrsp( 6,1)/1000.  ! VOD15U
     write(*,1100) ccbrsp( 7,1)/1000.  ! VOD15L
     write(*,1100) ccbrsp( 7,1)/1000.  ! VOD14L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VOD13L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VOD12L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VOD11L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VOD10L
     write(*,1100) ccbrsp( 9,1)/1000.  ! VOD9L
     write(*,1100) ccbrsp( 9,1)/1000.  ! VOD8L
     write(*,1100) ccbrsp(10,1)/1000.  ! VOD7L
     write(*,1100) ccbrsp(10,1)/1000.  ! VOD6L
     write(*,1100) ccbrsp(11,1)/1000.  ! VOD5L
     write(*,1100) ccbrsp(11,1)/1000.  ! VOD4L
     write(*,1100) ccbrsp(11,1)/1000.  ! VOD3L
     write(*,1100) ccbrsp(12,1)/1000.  ! VOD2L
     write(*,1100) ccbrsp(12,1)/1000.  ! VOD1L

     write(*,1100) ccbrsp( 1,1)/1000.  ! VICS1U
     write(*,1100) ccbrsp( 1,1)/1000.  ! VICS2U
     write(*,1100) ccbrsp( 2,1)/1000.  ! VICS3U
     write(*,1100) ccbrsp( 2,1)/1000.  ! VID4U
     write(*,1100) ccbrsp( 2,1)/1000.  ! VID5U
     write(*,1100) ccbrsp( 3,1)/1000.  ! VID6U
     write(*,1100) ccbrsp( 3,1)/1000.  ! VID7U
     write(*,1100) ccbrsp( 4,1)/1000.  ! VID8U
     write(*,1100) ccbrsp( 4,1)/1000.  ! VID9U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VID10U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VID11U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VID12U
     write(*,1100) ccbrsp( 5,1)/1000.  ! VID13U
     write(*,1100) ccbrsp( 6,1)/1000.  ! VIW14U
     write(*,1100) ccbrsp( 6,1)/1000.  ! VIW15U
     write(*,1100) ccbrsp( 7,1)/1000.  ! VIW15L
     write(*,1100) ccbrsp( 7,1)/1000.  ! VIW14L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VID13L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VID12L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VID11L
     write(*,1100) ccbrsp( 8,1)/1000.  ! VID10L
     write(*,1100) ccbrsp( 9,1)/1000.  ! VID9L
     write(*,1100) ccbrsp( 9,1)/1000.  ! VID8L
     write(*,1100) ccbrsp(10,1)/1000.  ! VID7L
     write(*,1100) ccbrsp(10,1)/1000.  ! VID6L
     write(*,1100) ccbrsp(11,1)/1000.  ! VID5L
     write(*,1100) ccbrsp(11,1)/1000.  ! VID4L
     write(*,1100) ccbrsp(11,1)/1000.  ! VICS3L
     write(*,1100) ccbrsp(12,1)/1000.  ! VICS2L
     write(*,1100) ccbrsp(12,1)/1000.  ! VICS1L
  else if(nfcoil0.eq.101) then
     write(0,*) 'Assuming MAST'
     write(*,1100) ccbrsp( 1,1)/1000.  ! Each turn in OH receives only 
     write(*,1100) ccbrsp( 2,1)/1000.   ! p2iu
     write(*,1100) ccbrsp( 3,1)/1000.    ! p2ou
     write(*,1100) ccbrsp( 4,1)/1000.   ! p2il
     write(*,1100) ccbrsp( 5,1)/1000.    ! p2ol
     write(*,1100) ccbrsp( 6,1)/1000./3.    ! p3u two left windings
     write(*,1100) ccbrsp( 6,1)/1000./3.    ! p3u four center windings
     write(*,1100) ccbrsp( 6,1)/1000./3.    ! p3u two right windings
     write(*,1100) ccbrsp( 7,1)/1000./3.    ! p3l two left windings
     write(*,1100) ccbrsp( 7,1)/1000./3.    ! p3l four center windings
     write(*,1100) ccbrsp( 7,1)/1000./3.    ! p3l two right windings
     write(*,1100) ccbrsp( 8,1)/1000./3.    ! p4u four slightly shifted windings
     write(*,1100) ccbrsp( 8,1)/1000./3.    ! p4u one winding next to shifted ones
     write(*,1100) ccbrsp( 8,1)/1000./3.   ! p4u 6x3 windings in rectangular arrangement
     write(*,1100) ccbrsp( 9,1)/1000./3.    ! p4l four slightly shifted windings
     write(*,1100) ccbrsp( 9,1)/1000./3.    ! p4l one winding next to shifted ones
     write(*,1100) ccbrsp( 9,1)/1000./3.   ! p4l 6x3 windings in rectangular arrangement
     write(*,1100) ccbrsp(10,1)/1000./3.    ! p5u four slightly shifted windings
     write(*,1100) ccbrsp(10,1)/1000./3.    ! p5u one winding next to shifted ones
     write(*,1100) ccbrsp(10,1)/1000./3.   ! p5u 6x3 windings in rectangular arrangement
     write(*,1100) ccbrsp(11,1)/1000./3.    ! p5l four slightly shifted windings
     write(*,1100) ccbrsp(11,1)/1000./3.    ! p5l one winding next to shifted ones
     write(*,1100) ccbrsp(11,1)/1000./3.   ! p5l 6x3 windings in rectangular arrangement
     write(*,1100) ccbrsp(12,1)/1000./2.    ! p6u two upper windings
     write(*,1100) ccbrsp(12,1)/1000./2.    ! p6u two lower windings
     write(*,1100) ccbrsp(13,1)/1000./2.    ! p6l two upper windings
     write(*,1100) ccbrsp(13,1)/1000./2.    ! p6l two lower windings
     write(*,1100) ccbrsp(14,1)/1000.  ! p2u case part 1
     write(*,1100) ccbrsp(14,1)/1000.  ! p2u case part 2
     write(*,1100) ccbrsp(14,1)/1000.  ! p2u case part 3
     write(*,1100) ccbrsp(14,1)/1000.  ! p2u case part 4
     write(*,1100) ccbrsp(15,1)/1000.  ! p2l case part 1
     write(*,1100) ccbrsp(15,1)/1000.  ! p2l case part 2
     write(*,1100) ccbrsp(15,1)/1000.  ! p2l case part 3
     write(*,1100) ccbrsp(15,1)/1000.  ! p2l case part 4
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 1
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 2
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 3
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 4
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 5
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 6
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 7
     write(*,1100) ccbrsp(16,1)/1000.  ! p3u case part 8
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 1
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 2
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 3
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 4
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 5
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 6
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 7
     write(*,1100) ccbrsp(17,1)/1000.  ! p3l case part 8
     write(*,1100) ccbrsp(18,1)/1000.  ! p4u case part 1
     write(*,1100) ccbrsp(18,1)/1000.  ! p4u case part 2
     write(*,1100) ccbrsp(18,1)/1000.  ! p4u case part 3
     write(*,1100) ccbrsp(18,1)/1000.  ! p4u case part 4
     write(*,1100) ccbrsp(19,1)/1000.  ! p4l case part 1
     write(*,1100) ccbrsp(19,1)/1000.  ! p4l case part 2
     write(*,1100) ccbrsp(19,1)/1000.  ! p4l case part 3
     write(*,1100) ccbrsp(19,1)/1000.  ! p4l case part 4
     write(*,1100) ccbrsp(20,1)/1000.  ! p5l case part 1
     write(*,1100) ccbrsp(20,1)/1000.  ! p5l case part 2
     write(*,1100) ccbrsp(20,1)/1000.  ! p5l case part 3
     write(*,1100) ccbrsp(20,1)/1000.  ! p5l case part 4
     write(*,1100) ccbrsp(21,1)/1000.  ! p6u case part 1
     write(*,1100) ccbrsp(21,1)/1000.  ! p6u case part 2
     write(*,1100) ccbrsp(21,1)/1000.  ! p6u case part 3
     write(*,1100) ccbrsp(21,1)/1000.  ! p6u case part 4
     write(*,1100) ccbrsp(22,1)/1000.  ! p6l case part 1
     write(*,1100) ccbrsp(22,1)/1000.  ! p6l case part 2
     write(*,1100) ccbrsp(22,1)/1000.  ! p6l case part 3
     write(*,1100) ccbrsp(22,1)/1000.  ! p6l case part 4
     write(*,1100) ccbrsp(23,1)/1000.    ! vertw1
     write(*,1100) ccbrsp(24,1)/1000.    ! vertw2
     write(*,1100) ccbrsp(25,1)/1000.    ! vertw3
     write(*,1100) ccbrsp(26,1)/1000.    ! vertw4
     write(*,1100) ccbrsp(27,1)/1000.    ! vertw5
     write(*,1100) ccbrsp(28,1)/1000.    ! vertw6
     write(*,1100) ccbrsp(29,1)/1000.    ! vertw7
     write(*,1100) ccbrsp(30,1)/1000.    ! vertw8
     write(*,1100) ccbrsp(31,1)/1000.    ! uhorw1
     write(*,1100) ccbrsp(32,1)/1000.    ! uhorw2
     write(*,1100) ccbrsp(33,1)/1000.    ! uhorw3
     write(*,1100) ccbrsp(34,1)/1000.    ! uhorw4
     write(*,1100) ccbrsp(35,1)/1000.    ! uhorw5
     write(*,1100) ccbrsp(36,1)/1000.    ! uhorw6
     write(*,1100) ccbrsp(37,1)/1000.    ! lhorw1
     write(*,1100) ccbrsp(38,1)/1000.    ! lhorw2
     write(*,1100) ccbrsp(39,1)/1000.    ! lhorw3
     write(*,1100) ccbrsp(40,1)/1000.    ! lhorw4
     write(*,1100) ccbrsp(41,1)/1000.    ! lhorw5
     write(*,1100) ccbrsp(42,1)/1000.    ! lhorw6
     write(*,1100) ccbrsp(43,1)/1000.    ! p2udivpl1
     write(*,1100) ccbrsp(44,1)/1000.    ! p2udivpl2
     write(*,1100) ccbrsp(45,1)/1000.    ! p2ldivpl1
     write(*,1100) ccbrsp(46,1)/1000.    ! p2ldivpl2
     write(*,1100) ccbrsp(47,1)/1000.    ! p2uarm1
     write(*,1100) ccbrsp(48,1)/1000.    ! p2uarm2
     write(*,1100) ccbrsp(49,1)/1000.    ! p2uarm3
     write(*,1100) ccbrsp(50,1)/1000.    ! p2larm1
     write(*,1100) ccbrsp(51,1)/1000.    ! p2larm2
     write(*,1100) ccbrsp(52,1)/1000.    ! p2larm3
     write(*,1100) ccbrsp(53,1)/1000.    ! topcol
     write(*,1100) ccbrsp(54,1)/1000.    ! incon1
     write(*,1100) ccbrsp(55,1)/1000.    ! incon2
     write(*,1100) ccbrsp(56,1)/1000.    ! incon3
     write(*,1100) ccbrsp(57,1)/1000.    ! incon4
     write(*,1100) ccbrsp(58,1)/1000.    ! incon5
     write(*,1100) ccbrsp(59,1)/1000.    ! incon6
     write(*,1100) ccbrsp(60,1)/1000.    ! incon7
     write(*,1100) ccbrsp(61,1)/1000.    ! incon8
     write(*,1100) ccbrsp(62,1)/1000.    ! incon9
     write(*,1100) ccbrsp(63,1)/1000.    ! incon10
     write(*,1100) ccbrsp(64,1)/1000.    ! botcol
!     write(*,1100) ccbrsp(65,1)/1000.    ! endcrown_u
!     write(*,1100) ccbrsp(66,1)/1000.    ! endcrown_l
!     write(*,1100) ccbrsp(67,1)/1000.    ! ring1
     write(*,1100) ccbrsp(68,1)/1000.    ! ring2
     write(*,1100) ccbrsp(69,1)/1000.    ! ring3
     write(*,1100) ccbrsp(70,1)/1000.    ! ring4
     write(*,1100) ccbrsp(71,1)/1000.    ! ring5
     write(*,1100) ccbrsp(72,1)/1000.    ! ring6
     write(*,1100) ccbrsp(73,1)/1000.    ! ring7
     write(*,1100) ccbrsp(74,1)/1000.    ! ring8
     write(*,1100) ccbrsp(75,1)/1000.    ! ring9
!     write(*,1100) ccbrsp(76,1)/1000.    ! ring10
     write(*,1100) ccbrsp(77,1)/1000.    ! rodgr1
     write(*,1100) ccbrsp(78,1)/1000.    ! rodgr2
     write(*,1100) ccbrsp(79,1)/1000.    ! rodgr3
     write(*,1100) ccbrsp(80,1)/1000.    ! rodgr4
     write(*,1100) ccbrsp(81,1)/1000.    ! rodgr5
     write(*,1100) ccbrsp(82,1)/1000.    ! rodgr6
     write(*,1100) ccbrsp(83,1)/1000.    ! rodgr7
     write(*,1100) ccbrsp(84,1)/1000.    ! rodgr8
     write(*,1100) ccbrsp(85,1)/1000.    ! rodgr9
     write(*,1100) ccbrsp(86,1)/1000.    ! rodgr10
     write(*,1100) ccbrsp(87,1)/1000.    ! rodgr11
     write(*,1100) ccbrsp(88,1)/1000.    ! rodgr12
     write(*,1100) ccbrsp(89,1)/1000.    ! mid1
     write(*,1100) ccbrsp(90,1)/1000.    ! mid2
     write(*,1100) ccbrsp(91,1)/1000.    ! mid3
     write(*,1100) ccbrsp(92,1)/1000.    ! mid4
     write(*,1100) ccbrsp(93,1)/1000.    ! mid5
     write(*,1100) ccbrsp(94,1)/1000.    ! mid6
     write(*,1100) ccbrsp(95,1)/1000.    ! mid7
     write(*,1100) ccbrsp(96,1)/1000.    ! mid8
     write(*,1100) ccbrsp(97,1)/1000.    ! mid9
     write(*,1100) ccbrsp(98,1)/1000.    ! mid10
     write(*,1100) ccbrsp(99,1)/1000.    ! mid11
     write(*,1100) ccbrsp(100,1)/1000.    ! mid12
     write(*,1100) ccbrsp(101,1)/1000.  ! p5u_case_current part 1
     write(*,1100) ccbrsp(101,1)/1000.  ! p5u_case_current part 2
     write(*,1100) ccbrsp(101,1)/1000.  ! p5u_case_current part 3
     write(*,1100) ccbrsp(101,1)/1000.  ! p5u_case_current part 4

  end if

1100 format (f12.4)

end subroutine load_eqdsk_a
end module eqdsk_a

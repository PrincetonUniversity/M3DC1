module eqdsk_a

  implicit none

  character(len=1) :: dum
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
  integer :: j, k, jj, ios
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
     read(neqdsk,1040) (time(j),j=1,1)

     write(0,*) uday, mfvers(1), mfvers(2)
     write(0,*) ' SHOT, time = ', ishot, ktime1

     do jj=1, ktime1
        read (neqdsk,1060) dum,time(jj),jflag(jj),lflag,limloc(jj), &
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
        read (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
        read (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
        read (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
        read (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)
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
        read (neqdsk,1041) nsilop0,magpri0,nfcoil0,nesum0
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
  else if(nfcoil0.eq.52) then
     write(0,*) 'Assuming NSTX'
     write(*,1100) 0.                  ! OH1U
     write(*,1100) 0.                  ! OH2U
     write(*,1100) 0.                  ! OH3U
     write(*,1100) 0.                  ! OH4U
     write(*,1100) 0.                  ! OH4L
     write(*,1100) 0.                  ! OH3L
     write(*,1100) 0.                  ! OH2L
     write(*,1100) 0.                  ! OH1L
     write(*,1100) ccbrsp( 1,1)*20/1000.  ! PF1AU
     write(*,1100) ccbrsp( 2,1)*14/1000.  ! PF2U1
     write(*,1100) ccbrsp( 2,1)*14/1000.  ! PF2U2
     write(*,1100) ccbrsp( 3,1)*15/1000.  ! PF3U1
     write(*,1100) ccbrsp( 3,1)*15/1000.  ! PF3U2
     write(*,1100) ccbrsp( 4,1)*9 /1000.  ! PF4U1
     write(*,1100) ccbrsp( 4,1)*8 /1000.  ! PF4U2
     write(*,1100) ccbrsp( 5,1)*12/1000.  ! PF5U1
     write(*,1100) ccbrsp( 5,1)*12/1000.  ! PF5U2
     write(*,1100) ccbrsp( 6,1)*12/1000.  ! PF5L2
     write(*,1100) ccbrsp( 6,1)*12/1000.  ! PF5L1
     write(*,1100) ccbrsp( 7,1)*8 /1000.  ! PF4L2
     write(*,1100) ccbrsp( 7,1)*9 /1000.  ! PF4L1
     write(*,1100) ccbrsp( 8,1)*15/1000.  ! PF3L2
     write(*,1100) ccbrsp( 8,1)*15/1000.  ! PF3L1
     write(*,1100) ccbrsp( 9,1)*14/1000.  ! PF2L2
     write(*,1100) ccbrsp( 9,1)*14/1000.  ! PF2L1
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
  end if

1100 format (f12.4)

end subroutine load_eqdsk_a
end module eqdsk_a

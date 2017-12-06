c
      include 'jsolp.i'
c
      character*8 name(10)
      common/com1/ bgcm1,
     1  rjtor(z1,z2),   vecr(z1,z2),    vec1(z1,z2),   vec2(z1,z2),
     2  vec3(z1,z2),   vecr1(z1,z2),   vec6(z1,z2),  vecro(z1,z2) ,
     3  vec7(z1,z2),    vec8(z1,z2),   vec9(z1,z2),
     4  hary(z1,z2),     aa1(z2,z1),    vtd(z1,z2), aa2(z1,z2),
     5  vta(z1,z2),      vtb(z1,z2),    vtc(z1,z2),
     6  aj(z1,z2),   vec5(z1,z2), vec5a(z1,z2), vecrs(z1,z2),
     7  vec5b(z1,z2),vra(z2,z1),vrb(z2,z1),vrc(z2,z1),
     8  vgtsc(z1,z2),vrd(z2,z1),     vgtsa(z1,z2),
     9  aj3(z1,z2),    psi(z1,z2),    ajbbs(z2),
     a  x(z1,z2),       z(z1,z2),    ajbtot(z2),
     b  xo(z1,z2),     zo(z1,z2),   ajdbbs(z2),
     c  xtemp(z4),   ytemp(z4),   ztemp(z4),
     d  ytemp2(z4),  ytemp3(z4),  pprime(z2),
     e  btary(z4),      bpary(z4),       prary(z4),
     f  psival(z2), psiform(z2),q(z2),  qpsave(z2),  vstre(z2),
     g  xp(z3),      zp(z3),     xpn(z3),     zpn(z3),  vstri(z2),
     h  theval(z1), ptsave(z1), ftrap(z2)  , ftofc(z2),
     i  ajtary(z2),   absary(z2),   acdary(z2), aa2sv(2*z1*z2),
     j ggp1old(z2),cdens(z2),te3hps(z2),ajavby(z2),btsobs(z2),
     . btsobs2(z2),ctot(z2),cdia(z2),cbs(z2),cps(z2),ccd(z2),
     . bsqary(z2),
     k               fncm1
c
      common/com2/ ibgcm2,
     1  npsi, npsim, npsip, nthe, nthe2, nthe3, 
     2  nthe4, nimax, isym,name, nthep, nisave,jfunc2,
     3  ifunc2,mthin,nsrfin,icur,ifile,it,iczero,
     .  mnthe,mnpsi,ni,ni2,imore,iskip,
     4  numit,itooff,jorgn,icurm,ipsi, ntfcoil,
     5  ipstbin,ishp,nnp(25),ineg,irst1,irst2,iplot,iboot,nspec,
     6  jsaw, ilhcd, idoub, ndoub,irip,npitch,
     f  nsc1,  nout, nin, itime, idate, imach, isuffix, idone,
     j  ndbfs, ndbfp,ilhcd2,ilhcd3,jsao,nsaw, npsit, kmax,
     k              ifncm2
c
      common/com2a/zlam, xjfac,czero,eta0,qmin,qminx,xnrat,scalep,
     1   vloopa,vloopp,det1,det2,ad1,bd1,zmin,zzero,zaxe,
     6  tsf,aasaw,psisaw,aasawo,aasawn,qsaw,ad2,bd2,dzero,ptacur,
     7  sf,phi2,dr,dt,rdr,rdt,at1,bt1,at2,bt2,tet1,tet2,
     8  rdrsq,rdtsq,xzero,ap1,ap2,bp1,bp2,ac1,bc1,ac2,bc2,
     9  p0,psilim,delpsi,psimin,roff,facimp,cet1,cet2,
     a  pi,tpi,psiave,xmin,amu0,xaxe,xli1,xli3,troyonfs,
     b  aguess,betap,betat,betatot,psirnge,twant(z2),rsaw,
     c  btor0,alphaq,qwant,sumdb,zeff,rwant(z2),vloops,dbwrt,
     d  vj(z2),betaq,qwants,p0s, t1, t2, t3,tzero,xmi,zave, omrip,
     g  pet1, pet2, bguess, dzguess, cedge, tsfsav, voltot,esaw,
     h  akapp0,cdens0,qstar,qratio,bpav(z2),rav(z2),tsaw,
     i  betapo,betato,betats,aima,btor,qstar2,qratio2,troyonf,
     j  aicd2,alh2,dlh2,a1lh2,a2lh2,aicd, alh, dlh, a1lh, a2lh, faclh,
     k  aicd3,alh3,dlh3,a1lh3,a2lh3,faclh2,faclh3,
     l  alfatot,aitotx,dlpp,dlp
c
      common/com3/ bgcm3,
     1  times,xaxes,zmags,gzeros,apls,betas,betaps,ali2s,
     2  qsaws,psimins,psilims,prpest2(z5),pppest2(z5),
     3  prpest0(z5),pppest0(z5),
     4  ajpest2(z5),xsv2(z5),xpass(z6),zpass(z6),
     5  xpasr(z6),zpasr(z6),elgu(z2),elgl(z2),triu(z2),tril(z2),
     6  asv(z5,8),bsv(z5,8),csv(z5,8),dsv(z5,8),
     7  xsv(z5), rmaj(z2), rmin(z2), elg(z2),tri(z2),shift(z2),
     8  eps,epsu,epsl,dguess,duguess,dlguess,
     9  f(z2),       g(z2),      gs(z2),  vl(z2),
     a  p(z2),       r(z2),         psisav(z2),  xdl(z2),
     b  qtemp(z2), xnuc(z2,10), xmass(10), xchg(10), tzer(10), dzer(10),
     c  vp(z2),gxmja(z2),gtxmja(z2), xmja(z2),vp2(z2),gxmja2(z2),
     d  xja(z2),xmja2(z2),fc(z2),pc(z2),tgs(z2),gx(z2),bpldl(z2),
     1  gja(z2), gja2(z2),
     e  ppnaa, ppnsa, ppnac, ppnsc, ppnae, ppnse, pppsnc,
     f  cpnaa, cpnsa, cpnac, cpnsc, cpnae, cpnse, cppsnc,
     g  tpnaa, tpnsa, tpnac, tpnsc, tpnae, tpnse, tppsnc,
     h  dpnaa, dpnsa, dpnac, dpnsc, dpnae, dpnse, dppsnc,
     i             fncm3
c
      common/scr1/ vt1(z2),vt2(z2),vt3(z2),vt7(z2),
     1             vt8(z2),vt9(z2),vt0(z2),psinp(z2),
     2             ggpx(z2),ggp1(z2),ggp2(z2)
      common/scr2/ va(z2),vb(z2),vc(z2),vd(z2),ve(z2),vf(z2),
     1             facv(z2),psiav(z2),vd1(z2),vd2(z2),vh(z2),vg(z2)
      common/rippl/ripcon(z1,z2),vseg(z1,z2),bmagfc(z1,z2),
     1delcrit(z2),delmax(z2),aplosts(z2),aplostsr(z2),adens(z2),
     2volume(z2),alphanu(z2),aplost(z1,z2),aplostr(z1,z2),
     3apzone(z1,z2),fraclos,fraclosr, rtfcoil,fdt,
     4tempr(z2),dens(z2),vencl(z2),phi(z2)

      common/new/ te2(z5),ane2(z5),ti2(z5),avez2(z5),zeff2(z5),
     1            te(z5) ,ane (z5),ti (z5),avez (z5),zeffa(z5),chi(z5),
     2          vloop(z5)








ccc      cliche l22com
c
      parameter ( ndimlpnl = 129, ndimlpc = 201 )
c
      logical lfunin, checkd, symvac
c
      common / l222 / lfunin, checkd, symvac, lchecke
c
      common / temstr / vacmat(nfm,nfm), vacmat0(nfm,nfm),
     $     vacmat00(nfm,nfm), vacmt00i(nfm,nfm),
     $     vacmatu(nfm,nfm),vacmtiu(nfm,nfm),
     $     vcmattu0(nfm2,nfm2),vcmattu(nfm2,nfm2),
     $     gatovac(nfm,nfm), venerg, wpp(nfm,nfm), wppi(nfm,nfm),
     $     wpp0(nfm,nfm), wpp0i(nfm,nfm),
     $     oplp(nfm2+nwm4,nfm2+nwm4), clrsp(nths)
c
      common / vacval / aw,awx,awv,bw,cw,dw,rw,tw, abulg, bbulg, tbulg,
     $     epsq, delg, cn0, idgt, idot, nsing, noutv, ishape
c
      common / greeny / xt, zt, xs, zs, xtp, ztp,
     .                  aval, bval, ieps, lmax1
c
      common / elipse / xpl, apl, bpl, dpl, epl, fpl, gpl, hpl,
     $     zpoff, a, b, qain, ldelta, adel
c
      common / diags / xloopnl(ndimlpnl), zloopnl(ndimlpnl), xofsl, 
     $     aloop, bloop, dloop, rloop, deloop, delfac, delx, delz,
     $     lpsub, nphil, nphse, ieig, ismth, mx, mz, nph,
     $     ntloop,nloop,nloopr, nxlpin, nzlpin, epslp,
     $     xlpmin, xlpmax, zlpmin, zlpmax, linterior, lchkinextp,
     $     lclosepts
c
      common / bsrce / xilr(nfm), xili(nfm), xilri(nfm2),
     $     bnlr(nfm), bnli(nfm), bnlri(nfm2),
     $     bnkr(nths), bnki(nths)
c
      common / sprkin / civ,
     $     sp2sgn1, sp2sgn2, sp2sgn3, sp2sgn4, sp2sgn5,
     $     sp3sgn1, sp3sgn2, sp3sgn3, sp3sgn4, sp3sgn5,
     $     ff, fv(nfm), nminus, nplus, mphi, lwrt11, lff
c
c     **************************************************************
c

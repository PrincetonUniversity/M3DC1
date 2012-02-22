c
      COMMON / rshelv / gwin(nths2,nwm2), gwot(nths2,nwm2),
     $     blsp(nwm2,nfm2), xlwm(nwm2)
c
      COMMON / rsheli / lwmn, lwmx, jwal1, jwal2, jwal4, nwcur, gmat,
     $     irotplasma
c
      COMMON / eigens  / neig(nwm2), nwke(nwm), nwko(nwm),
     $     rwval(nwm2), rwvle(nwm), rwvlo(nwm),
     $     rwvec(nths,nwm2), rwvce(nths,nwm), rwvco(nths,nwm),
     $     rwvle0(nwm), rwvlo0(nwm),
     $     rwvce0(nths,nwm), rwvco0 (nths,nwm)
c
      COMMON / sprops / eta_vess, th_vess, tc_vess(nwm2), taugam
c
      COMMON / fcoil / xcw(nths), zcw(nths), ccoil(nths),
     $     xcwp(nths), zcwp(nths), clena_c(nths), clen_c, 
     $     ellc(nths), shcoil(nths),
     $     skcoil(nths), coil00(nfm2), coilp1(nfm2), coilm1(nfm2),
     $     gainc, iopfc, icshp
      COMMON / fcoila / xcwa(nths), zcwa(nths),  ccoila(nths),
     $     xcwap(nths), zcwap(nths), clena_u(nths), clen_u, 
     $     skcoilpa(nths),skcoilwa(nths), coil00a(nfm2), gainca,
     $     iopfca, icshpa
      COMMON / fcoilb / xcwb(nths), zcwb(nths), ccoilb(nths),
     $     xcwbp(nths), zcwbp(nths), clena_l(nths), clen_l, 
     $     skcoilpb(nths),skcoilwb(nths), coil00b(nfm2), gaincb,
     $     iopfcb, icshpb

      COMMON / cprops / drcoil, acoil, bcoil, dcoil, coilang, dzcoil,
     $     taucoil, coilcen, chgt(4), mtcoil, mtcoil1, ncgr,
     $     ctip_c(4), ntip_c
      COMMON / cpropsa / drcoila, acoila, bcoila, dcoila, coilanga,
     $     dzcoila, taucoila, coilcena, chgta(4), mtcoila, mtcoila1,
     $     ctip_u(4), ntip_u
      COMMON / cpropsb / drcoilb, acoilb, bcoilb, dcoilb, coilangb,
     $     dzcoilb, taucoilb, coilcenb, chgtb(4), mtcoilb, mtcoilb1,
     $     ctip_l(4), ntip_l

      COMMON / sensl / sang1, sang2, dtsl, asloop, isloop, nslpt,
     $     lsensplot,lrswcomp,lpolsensor


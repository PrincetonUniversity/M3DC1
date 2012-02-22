c
      logical farwal
      common / map2 / xinf(nths), zinf(nths), xplap(nths), zplap(nths),
     $     dlpdth(nths), xwin(nths), zwin(nths)
      common /wpparam/ pelong, pdness, wmaj,wzma, wlrad, welong, wdness,
     $     pcircm, wcircm, wcentr
c
      common/kplt/xpass(nths),zpass(nths),chiwc(nths,nfm),
     $     chiws(nths,nfm), ell(nths)
      common /kplt1/ mw,jtop,jbot, farwal
      common / axes / zerov(nths), tngrd(nths), wgrd(nths), pgrd(nths),
     $     pigrd(nths),piwgrd(nths), xlfm(nfm), lfm(nfm), lrnge
c
c     *****************************************************************
c


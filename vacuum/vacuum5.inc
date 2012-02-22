c
c..   cplar, cplai are the plasma-plasma response matrices.  
c..   cplae, cplao are the plasma-wall response matrices to even, odd bases.
c..   cwallr, cwalli are the wall-plasma response matrices.  
c..   cwwale, cwwalo are the wall-wall response matrices to even, odd bases.

      common /arys/ xmaj,plrad,
     $     xwal(nths),zwal(nths), xwalp(nths),zwalp(nths), elwal(nths),
     $     xwalq(nths),zwalq(nths), xwalpq(nths),zwalpq(nths),
     $     angwal(nths), angwalq(nths),
     $     cnqd(nths), snqd(nths),
     $     cslth(nths,nfm), snlth(nths,nfm),
     $     coslt(nths,nfm), sinlt(nths,nfm),
     $     cplar(nths,nfm), cplai(nths,nfm),
     $     cplae(nths,nfm), cplao(nths,nfm),
     $     cppr(nths,nfm), cppi(nths,nfm),
     $     cwallr(nths,nfm), cwalli(nths,nfm),
     $     cwwale(nths,nfm), cwwalo(nths,nfm),
     $     cwalr(nths,nfm), cwali(nths,nfm),
     $     cpwr(nths,nfm),  cpwi(nths,nfm),
     $     cwwtr(nths,nfm), cwwti(nths,nfm),
     $     wrktr1(nfm,nths), wrktr2(nfm,nths),
     $     chipr(nths), chipi(nths), bthpr(nths), bthpi(nths),
     $     gpsjp(nths), gpsjw(nths)
c
      common /const2/ mdel,jdel
c
      common /zwork/ zork1(nths),zork2(nths),zork3(nths),zork4(nths)
c
c****************************************************************

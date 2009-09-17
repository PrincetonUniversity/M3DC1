      function alog (x)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      dimension alncs(6), center(4), alncen(5)
      external csevl, inits, r1mach
c
c series for aln        on the interval  0.          to  3.46021d-03
c                                        with weighted error   1.50e-16
c                                         log weighted error  15.82
c                               significant figures required  15.65
c                                    decimal places required  16.21
c
      data aln cs( 1) /   1.3347199877 973882e0 /
      data aln cs( 2) /    .0006937562 83284112e0 /
      data aln cs( 3) /    .0000004293 40390204e0 /
      data aln cs( 4) /    .0000000002 89338477e0 /
      data aln cs( 5) /    .0000000000 00205125e0 /
      data aln cs( 6) /    .0000000000 00000150e0 /
c
      data center(1) / 1.0 /
      data center(2) / 1.25 /
      data center(3) / 1.50 /
      data center(4) / 1.75 /
c
      data alncen(  1) / 0.0e0                                         /
      data alncen(  2) / +.2231435513 14209755 e+0                     /
      data alncen(  3) / +.4054651081 08164381 e+0                     /
      data alncen(  4) / +.5596157879 35422686 e+0                     /
      data alncen(  5) / +.6931471805 59945309 e+0                     /
c
c aln2 = alog(2.0) - 0.625
      data aln2 / 0.0681471805 59945309e0 /
      data nterms / 0 /
c
      if (nterms.eq.0) nterms = inits (alncs, 6, 28.9*r1mach(3))
c
      if (x.le.0.) call seteru (
     1  29halog    x is zero or negative, 29, 1, 2)
c
      call r9upak (x, y, n)
c
      xn = n - 1
      y = 2.0*y
      ntrval = 4.0*y - 2.5
      if (ntrval.eq.5) t = ((y-1.0)-1.0) / (y+2.0)
      if (ntrval.lt.5) t = (y-center(ntrval))/(y+center(ntrval))
      t2 = t*t
c
      alog = 0.625*xn + (aln2*xn + alncen(ntrval) + 2.0*t +
     1  t*t2*csevl(578.0*t2-1.0, alncs, nterms) )
c
      return
      end
      function csevl (x, cs, n)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the n-term chebyshev series cs at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
c and parker, chebyshev polys in numerical analysis, oxford press, p.56.
c
c             input arguments --
c x      value at which the series is to be evaluated.
c cs     array of n terms of a chebyshev series.  in eval-
c        uating cs, only half the first coef is summed.
c n      number of terms in array cs.
c
      dimension cs(1)
c
      if (n.lt.1) call seteru (28hcsevl   number of terms le 0, 28, 2,2)
      if (n.gt.1000) call seteru (31hcsevl   number of terms gt 1000,
     1  31, 3, 2)
      if (x.lt.(-1.1) .or. x.gt.1.1) call seteru (
     1  25hcsevl   x outside (-1,+1), 25, 1, 1)
c
      b1 = 0.
      b0 = 0.
      twox = 2.*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
c
      csevl = 0.5 * (b0-b2)
c
      return
      end
      subroutine d9b0mp (x, ampl, theta)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the modulus and phase for the bessel j0 and y0 functions.
c
      double precision x, ampl, theta, bm0cs(37), bt02cs(39),
     1  bm02cs(40), bth0cs(44), xmax, pi4, z, d1mach, dcsevl,
     2  dsqrt
      external d1mach, dcsevl, dsqrt, initds
c
c series for bm0        on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   4.40e-32
c                                         log weighted error  31.36
c                               significant figures required  30.02
c                                    decimal places required  32.14
c
      data bm0 cs(  1) / +.9211656246 8277427125 7376773018 2 d-1      /
      data bm0 cs(  2) / -.1050590997 2719051024 8071637175 5 d-2      /
      data bm0 cs(  3) / +.1470159840 7687597540 5639285095 2 d-4      /
      data bm0 cs(  4) / -.5058557606 0385542233 4792932770 2 d-6      /
      data bm0 cs(  5) / +.2787254538 6324441766 3035613788 1 d-7      /
      data bm0 cs(  6) / -.2062363611 7809148026 1884101897 3 d-8      /
      data bm0 cs(  7) / +.1870214313 1388796751 3817259626 1 d-9      /
      data bm0 cs(  8) / -.1969330971 1356362002 4173077782 5 d-10     /
      data bm0 cs(  9) / +.2325973793 9992754440 1250881805 2 d-11     /
      data bm0 cs( 10) / -.3009520344 9382502728 5122473448 2 d-12     /
      data bm0 cs( 11) / +.4194521333 8506691814 7120676864 6 d-13     /
      data bm0 cs( 12) / -.6219449312 1884458259 7326742956 4 d-14     /
      data bm0 cs( 13) / +.9718260411 3360684696 0176588526 9 d-15     /
      data bm0 cs( 14) / -.1588478585 7010752073 6663596693 7 d-15     /
      data bm0 cs( 15) / +.2700072193 6713088900 8621732445 8 d-16     /
      data bm0 cs( 16) / -.4750092365 2340089924 7750478677 3 d-17     /
      data bm0 cs( 17) / +.8615128162 6043708731 9170374656 0 d-18     /
      data bm0 cs( 18) / -.1605608686 9561448157 4560270335 9 d-18     /
      data bm0 cs( 19) / +.3066513987 3144829751 8853980159 9 d-19     /
      data bm0 cs( 20) / -.5987764223 1939564306 9650561706 6 d-20     /
      data bm0 cs( 21) / +.1192971253 7482483064 8906984106 6 d-20     /
      data bm0 cs( 22) / -.2420969142 0448054894 8468258133 3 d-21     /
      data bm0 cs( 23) / +.4996751760 5106164533 7100287999 9 d-22     /
      data bm0 cs( 24) / -.1047493639 3511585100 9504051199 9 d-22     /
      data bm0 cs( 25) / +.2227786843 7974681010 4818346666 6 d-23     /
      data bm0 cs( 26) / -.4801813239 3981628623 7054293333 3 d-24     /
      data bm0 cs( 27) / +.1047962723 4709599564 7699626666 6 d-24     /
      data bm0 cs( 28) / -.2313858165 6786153251 0126080000 0 d-25     /
      data bm0 cs( 29) / +.5164823088 4626742116 3519999999 9 d-26     /
      data bm0 cs( 30) / -.1164691191 8500653895 2540159999 9 d-26     /
      data bm0 cs( 31) / +.2651788486 0433192829 5833600000 0 d-27     /
      data bm0 cs( 32) / -.6092559503 8257284976 9130666666 6 d-28     /
      data bm0 cs( 33) / +.1411804686 1442593080 3882666666 6 d-28     /
      data bm0 cs( 34) / -.3298094961 2317372457 5061333333 3 d-29     /
      data bm0 cs( 35) / +.7763931143 0740650317 1413333333 3 d-30     /
      data bm0 cs( 36) / -.1841031343 6614584784 2133333333 3 d-30     /
      data bm0 cs( 37) / +.4395880138 5943107371 0079999999 9 d-31     /
c
c series for bth0       on the interval  0.          to  1.56250e-02
c                                        with weighted error   2.66e-32
c                                         log weighted error  31.57
c                               significant figures required  30.67
c                                    decimal places required  32.40
c
      data bth0cs(  1) / -.2490178086 2128936717 7097937899 67 d+0     /
      data bth0cs(  2) / +.4855029960 9623749241 0486155354 85 d-3     /
      data bth0cs(  3) / -.5451183734 5017204950 6562735635 05 d-5     /
      data bth0cs(  4) / +.1355867305 9405964054 3774459299 03 d-6     /
      data bth0cs(  5) / -.5569139890 2227626227 5832184149 20 d-8     /
      data bth0cs(  6) / +.3260903182 4994335304 0042057194 68 d-9     /
      data bth0cs(  7) / -.2491880786 2461341125 2379038779 93 d-10    /
      data bth0cs(  8) / +.2344937742 0882520554 3524135648 91 d-11    /
      data bth0cs(  9) / -.2609653444 4310387762 1775747661 36 d-12    /
      data bth0cs( 10) / +.3335314042 0097395105 8699550149 23 d-13    /
      data bth0cs( 11) / -.4789000044 0572684646 7507705574 09 d-14    /
      data bth0cs( 12) / +.7595617843 6192215972 6425685452 48 d-15    /
      data bth0cs( 13) / -.1313155601 6891440382 7733974876 33 d-15    /
      data bth0cs( 14) / +.2448361834 5240857495 4268207383 55 d-16    /
      data bth0cs( 15) / -.4880572981 0618777683 2567619183 31 d-17    /
      data bth0cs( 16) / +.1032728502 9786316149 2237563612 04 d-17    /
      data bth0cs( 17) / -.2305763381 5057217157 0047445270 25 d-18    /
      data bth0cs( 18) / +.5404444300 1892693993 0171084837 65 d-19    /
      data bth0cs( 19) / -.1324069519 4366572724 1550328823 85 d-19    /
      data bth0cs( 20) / +.3378079562 1371970203 4247921247 22 d-20    /
      data bth0cs( 21) / -.8945762915 7111779003 0269262922 99 d-21    /
      data bth0cs( 22) / +.2451990688 9219317090 8999086514 05 d-21    /
      data bth0cs( 23) / -.6938842287 6866318680 1399331576 57 d-22    /
      data bth0cs( 24) / +.2022827871 4890138392 9463033377 91 d-22    /
      data bth0cs( 25) / -.6062850000 2335483105 7941953717 64 d-23    /
      data bth0cs( 26) / +.1864974896 4037635381 8237883962 70 d-23    /
      data bth0cs( 27) / -.5878373238 4849894560 2450365308 67 d-24    /
      data bth0cs( 28) / +.1895859144 7999563485 5311795035 13 d-24    /
      data bth0cs( 29) / -.6248197937 2258858959 2916207285 65 d-25    /
      data bth0cs( 30) / +.2101790168 4551024686 6386335290 74 d-25    /
      data bth0cs( 31) / -.7208430093 5209253690 8139339924 46 d-26    /
      data bth0cs( 32) / +.2518136389 2474240867 1564059767 46 d-26    /
      data bth0cs( 33) / -.8951804225 8785778806 1439459536 43 d-27    /
      data bth0cs( 34) / +.3235723747 9762298533 2562358685 87 d-27    /
      data bth0cs( 35) / -.1188301051 9855353657 0471441137 96 d-27    /
      data bth0cs( 36) / +.4430628690 7358104820 5792319417 31 d-28    /
      data bth0cs( 37) / -.1676100964 8834829495 7920101356 81 d-28    /
      data bth0cs( 38) / +.6429294692 1207466972 5323939660 88 d-29    /
      data bth0cs( 39) / -.2499226116 6978652421 2072136827 63 d-29    /
      data bth0cs( 40) / +.9839979429 9521955672 8282603553 18 d-30    /
      data bth0cs( 41) / -.3922037524 2408016397 9891316261 58 d-30    /
      data bth0cs( 42) / +.1581810703 0056522138 5906188456 92 d-30    /
      data bth0cs( 43) / -.6452550614 4890715944 3440983654 26 d-31    /
      data bth0cs( 44) / +.2661111136 9199356137 1770183463 67 d-31    /
c
c series for bm02       on the interval  0.          to  1.56250e-02
c                                        with weighted error   4.72e-32
c                                         log weighted error  31.33
c                               significant figures required  30.00
c                                    decimal places required  32.13
c
      data bm02cs(  1) / +.9500415145 2283813693 3086133556 0 d-1      /
      data bm02cs(  2) / -.3801864682 3656709917 4808156685 1 d-3      /
      data bm02cs(  3) / +.2258339301 0314811929 5182992722 4 d-5      /
      data bm02cs(  4) / -.3895725802 3722287647 3062141260 5 d-7      /
      data bm02cs(  5) / +.1246886416 5120816979 3099052972 5 d-8      /
      data bm02cs(  6) / -.6065949022 1025037798 0383505838 7 d-10     /
      data bm02cs(  7) / +.4008461651 4217469910 1527597104 5 d-11     /
      data bm02cs(  8) / -.3350998183 3980942184 6729879457 4 d-12     /
      data bm02cs(  9) / +.3377119716 5174173670 6326434199 6 d-13     /
      data bm02cs( 10) / -.3964585901 6350127005 6935629582 3 d-14     /
      data bm02cs( 11) / +.5286111503 8838572173 8793974473 5 d-15     /
      data bm02cs( 12) / -.7852519083 4508523136 5464024349 3 d-16     /
      data bm02cs( 13) / +.1280300573 3866822010 1163407344 9 d-16     /
      data bm02cs( 14) / -.2263996296 3914297762 8709924488 4 d-17     /
      data bm02cs( 15) / +.4300496929 6567903886 4641029047 7 d-18     /
      data bm02cs( 16) / -.8705749805 1325870797 4753545145 5 d-19     /
      data bm02cs( 17) / +.1865862713 9620951411 8144277205 0 d-19     /
      data bm02cs( 18) / -.4210482486 0930654573 4508697230 1 d-20     /
      data bm02cs( 19) / +.9956676964 2284009915 8162741784 2 d-21     /
      data bm02cs( 20) / -.2457357442 8053133596 0592147854 7 d-21     /
      data bm02cs( 21) / +.6307692160 7620315680 8735370705 9 d-22     /
      data bm02cs( 22) / -.1678773691 4407401426 9333117238 8 d-22     /
      data bm02cs( 23) / +.4620259064 6739044337 7087813608 7 d-23     /
      data bm02cs( 24) / -.1311782266 8603087322 3769340249 6 d-23     /
      data bm02cs( 25) / +.3834087564 1163028277 4792244027 6 d-24     /
      data bm02cs( 26) / -.1151459324 0777412710 7261329357 6 d-24     /
      data bm02cs( 27) / +.3547210007 5233385230 7697134521 3 d-25     /
      data bm02cs( 28) / -.1119218385 8150046462 6435594217 6 d-25     /
      data bm02cs( 29) / +.3611879427 6298378316 9840499425 7 d-26     /
      data bm02cs( 30) / -.1190687765 9133331500 9264176246 3 d-26     /
      data bm02cs( 31) / +.4005094059 4039681318 0247644953 6 d-27     /
      data bm02cs( 32) / -.1373169422 4522123905 9519391601 7 d-27     /
      data bm02cs( 33) / +.4794199088 7425315859 9649152643 7 d-28     /
      data bm02cs( 34) / -.1702965627 6241095840 0699447645 2 d-28     /
      data bm02cs( 35) / +.6149512428 9363300715 0357516132 4 d-29     /
      data bm02cs( 36) / -.2255766896 5818283499 4430023724 2 d-29     /
      data bm02cs( 37) / +.8399707509 2942994860 6165835320 0 d-30     /
      data bm02cs( 38) / -.3172997595 5626023555 6742393615 2 d-30     /
      data bm02cs( 39) / +.1215205298 8812985545 8333302651 4 d-30     /
      data bm02cs( 40) / -.4715852749 7544386930 1321056804 5 d-31     /
c
c series for bt02       on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   2.99e-32
c                                         log weighted error  31.52
c                               significant figures required  30.61
c                                    decimal places required  32.32
c
      data bt02cs(  1) / -.2454829521 3424597462 0504672493 24 d+0     /
      data bt02cs(  2) / +.1254412103 9084615780 7853317782 99 d-2     /
      data bt02cs(  3) / -.3125395041 4871522854 9734467095 71 d-4     /
      data bt02cs(  4) / +.1470977824 9940831164 4534269693 14 d-5     /
      data bt02cs(  5) / -.9954348893 7950033643 4688503511 58 d-7     /
      data bt02cs(  6) / +.8549316673 3203041247 5787113977 51 d-8     /
      data bt02cs(  7) / -.8698975952 6554334557 9855121791 92 d-9     /
      data bt02cs(  8) / +.1005209953 3559791084 5401010821 53 d-9     /
      data bt02cs(  9) / -.1282823060 1708892903 4836236855 44 d-10    /
      data bt02cs( 10) / +.1773170078 1805131705 6557504510 23 d-11    /
      data bt02cs( 11) / -.2617457456 9485577488 6362841809 25 d-12    /
      data bt02cs( 12) / +.4082835138 9972059621 9664812211 03 d-13    /
      data bt02cs( 13) / -.6675166823 9742720054 6067495542 61 d-14    /
      data bt02cs( 14) / +.1136576139 3071629448 3924695499 51 d-14    /
      data bt02cs( 15) / -.2005118962 0647160250 5592664121 17 d-15    /
      data bt02cs( 16) / +.3649797879 4766269635 7205914641 06 d-16    /
      data bt02cs( 17) / -.6830963756 4582303169 3558437888 00 d-17    /
      data bt02cs( 18) / +.1310758314 5670756620 0571042679 46 d-17    /
      data bt02cs( 19) / -.2572336310 1850607778 7571306495 99 d-18    /
      data bt02cs( 20) / +.5152165744 1863959925 2677809493 33 d-19    /
      data bt02cs( 21) / -.1051301756 3758802637 9407414613 33 d-19    /
      data bt02cs( 22) / +.2182038199 1194813847 3010845013 33 d-20    /
      data bt02cs( 23) / -.4600470121 0362160577 2259054933 33 d-21    /
      data bt02cs( 24) / +.9840700692 5466818520 9536511999 99 d-22    /
      data bt02cs( 25) / -.2133403803 5728375844 7359863466 66 d-22    /
      data bt02cs( 26) / +.4683103642 3973365296 0662869333 33 d-23    /
      data bt02cs( 27) / -.1040021369 1985747236 5133823999 99 d-23    /
      data bt02cs( 28) / +.2334910567 7301510051 7777408000 00 d-24    /
      data bt02cs( 29) / -.5295682532 3318615788 0497493333 33 d-25    /
      data bt02cs( 30) / +.1212634195 2959756829 1962879999 99 d-25    /
      data bt02cs( 31) / -.2801889708 2289428760 2756266666 66 d-26    /
      data bt02cs( 32) / +.6529267898 7012873342 5937066666 66 d-27    /
      data bt02cs( 33) / -.1533798006 1873346427 8357333333 33 d-27    /
      data bt02cs( 34) / +.3630588430 6364536682 3594666666 66 d-28    /
      data bt02cs( 35) / -.8656075571 3629122479 1722666666 66 d-29    /
      data bt02cs( 36) / +.2077990997 2536284571 2383999999 99 d-29    /
      data bt02cs( 37) / -.5021117022 1417221674 3253333333 33 d-30    /
      data bt02cs( 38) / +.1220836027 9441714184 1919999999 99 d-30    /
      data bt02cs( 39) / -.2986005626 7039913454 2506666666 66 d-31    /
c
      data pi4 / 0.7853981633 9744830961 5660845819 876 d0 /
      data nbm0, nbt02, nbm02, nbth0, xmax / 4*0, 0.d0 /
c
      if (nbm0.ne.0) go to 10
      eta = 0.1*sngl(d1mach(3))
      nbm0 = initds (bm0cs, 37, eta)
      nbt02 = initds (bt02cs, 39, eta)
      nbm02 = initds (bm02cs, 40, eta)
      nbth0 = initds (bth0cs, 44, eta)
c
      xmax = 1.0d0/d1mach(4)
c
 10   if (x.lt.4.d0) call seteru (22hd9b0mp  x must be ge 4, 22, 1, 2)
c
      if (x.gt.8.d0) go to 20
      z = (128.d0/(x*x) - 5.d0)/3.d0
      ampl = (.75d0 + dcsevl (z, bm0cs, nbm0))/dsqrt(x)
      theta = x - pi4 + dcsevl (z, bt02cs, nbt02)/x
      return
c
 20   if (x.gt.xmax) call seteru (
     1  37hd9b0mp  no precision because x is big, 37, 2, 2)
c
      z = 128.d0/(x*x) - 1.d0
      ampl = (.75d0 + dcsevl (z, bm02cs, nbm02))/dsqrt(x)
      theta = x - pi4 + dcsevl (z, bth0cs, nbth0)/x
      return
c
      end
      double precision function d9pak (y, n)
c december 1979 edition. w. fullerton, c3, los alamos scientific lab.
c
c pack a base 2 exponent into floating point number x.  this routine is
c almost the inverse of d9upak.  it is not exactly the inverse, because
c dabs(x) need not be between 0.5 and 1.0.  if both d9pak and 2.d0**n
c were known to be in range we could compute
c                d9pak = x * 2.0d0**n
c
      double precision y, aln2b, aln210, d1mach
      external d1mach, i1mach
      data nmin, nmax / 2*0 /
      data aln210 / 3.321928094 8873623478 7031942948 9 d0 /
c
      if (nmin.ne.0) go to 10
      aln2b = 1.0d0
      if (i1mach(10).ne.2) aln2b = d1mach(5)*aln210
      nmin = aln2b*dble(float(i1mach(15)))
      nmax = aln2b*dble(float(i1mach(16)))
c
 10   call d9upak (y, d9pak, ny)
c
      nsum = n + ny
      if (nsum.lt.nmin) go to 40
      if (nsum.gt.nmax) call seteru (
     1  31hd9pak   packed number overflows, 31, 1, 2)
c
      if (nsum.eq.0) return
      if (nsum.gt.0) go to 30
c
 20   d9pak = 0.5d0*d9pak
      nsum = nsum + 1
      if (nsum.ne.0) go to 20
      return
c
 30   d9pak = 2.0d0*d9pak
      nsum = nsum - 1
      if (nsum.ne.0) go to 30
      return
c
 40   call seteru (32hd9pak   packed number underflows, 32, 1, 0)
      d9pak = 0.0d0
      return
c
      end
      subroutine d9upak (x, y, n)
c august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c unpack floating point number x so that x = y * 2.0**n, where
c 0.5 .le. abs(y) .lt. 1.0 .
c
      double precision x, y, absx
c
      absx = dabs(x)
      n = 0
      y = 0.0d0
      if (x.eq.0.0d0) return
c
 10   if (absx.ge.0.5d0) go to 20
      n = n - 1
      absx = absx*2.0d0
      go to 10
c
 20   if (absx.lt.1.0d0) go to 30
      n = n + 1
      absx = absx*0.5d0
      go to 20
c
 30   y = dsign (absx, x)
      return
c
      end
      double precision function dbesj0 (x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
      double precision x, bj0cs(19), ampl, theta, xsml, y,
     1  d1mach, dcsevl, dcos, dsqrt
      external d1mach, dcos, dcsevl, dsqrt, initds
c
c series for bj0        on the interval  0.          to  1.60000e+01
c                                        with weighted error   4.39e-32
c                                         log weighted error  31.36
c                               significant figures required  31.21
c                                    decimal places required  32.00
c
      data bj0 cs(  1) / +.1002541619 6893913701 0731272640 74 d+0     /
      data bj0 cs(  2) / -.6652230077 6440513177 6787578311 24 d+0     /
      data bj0 cs(  3) / +.2489837034 9828131370 4604687266 80 d+0     /
      data bj0 cs(  4) / -.3325272317 0035769653 8843415038 54 d-1     /
      data bj0 cs(  5) / +.2311417930 4694015462 9049241177 29 d-2     /
      data bj0 cs(  6) / -.9911277419 9508092339 0485193365 49 d-4     /
      data bj0 cs(  7) / +.2891670864 3998808884 7339037470 78 d-5     /
      data bj0 cs(  8) / -.6121085866 3032635057 8184074815 16 d-7     /
      data bj0 cs(  9) / +.9838650793 8567841324 7687486364 15 d-9     /
      data bj0 cs( 10) / -.1242355159 7301765145 5158970068 36 d-10    /
      data bj0 cs( 11) / +.1265433630 2559045797 9158272103 63 d-12    /
      data bj0 cs( 12) / -.1061945649 5287244546 9148175129 59 d-14    /
      data bj0 cs( 13) / +.7470621075 8024567437 0989155840 00 d-17    /
      data bj0 cs( 14) / -.4469703227 4412780547 6270079999 99 d-19    /
      data bj0 cs( 15) / +.2302428158 4337436200 5230933333 33 d-21    /
      data bj0 cs( 16) / -.1031914479 4166698148 5226666666 66 d-23    /
      data bj0 cs( 17) / +.4060817827 4873322700 8000000000 00 d-26    /
      data bj0 cs( 18) / -.1414383600 5240913919 9999999999 99 d-28    /
      data bj0 cs( 19) / +.4391090549 6698880000 0000000000 00 d-31    /
c
      data ntj0, xsml / 0, 0.d0 /
c
      if (ntj0.ne.0) go to 10
      ntj0 = initds (bj0cs, 19, 0.1*sngl(d1mach(3)))
      xsml = dsqrt (4.0d0*d1mach(3))
c
 10   y = dabs(x)
      if (y.gt.4.0d0) go to 20
c
      dbesj0 = 1.0d0
      if (y.gt.xsml) dbesj0 = dcsevl (.125d0*y*y-1.d0, bj0cs, ntj0)
      return
c
 20   call d9b0mp (y, ampl, theta)
      dbesj0 = ampl * dcos(theta)
c
      return
      end
      double precision function dcos (x)
c august 1980 edition.  w. fullerton, los alamos scientific lab.
c
c this routine is based on the algorithm of cody and waite in
c argonne tm-321, software manual working note number 1.
c
      double precision x, sincs(15), pihi, pilo, pi2, pirec, pi2rec,
     1  xsml, xwarn, xmax, absx, y, xn, f, dint, dcsevl, d1mach,
     2  dsqrt
      external d1mach, dcsevl, dint, dsqrt, initds
c
c series for sin    on the interval  0.00000e+00 to  2.46740e+00
c                                        with weighted error   2.56e-34
c                                         log weighted error  33.59
c                               significant figures required  33.01
c                                    decimal places required  34.18
c
      data sin cs(  1) / -0.3749911549 5587317583 9919279977 323464d0/
      data sin cs(  2) / -0.1816031552 3725020186 3830316158 004754d0/
      data sin cs(  3) /  0.0058047092 7459863355 9427341722 857921d0/
      data sin cs(  4) / -0.0000869543 1177934075 7113212316 353178d0/
      data sin cs(  5) /  0.0000007543 7014808885 1481006839 927030d0/
      data sin cs(  6) / -0.0000000042 6712966505 5961107126 829906d0/
      data sin cs(  7) /  0.0000000000 1698042294 5488168181 824792d0/
      data sin cs(  8) / -0.0000000000 0005012057 8889961870 929524d0/
      data sin cs(  9) /  0.0000000000 0000011410 1026680010 675628d0/
      data sin cs( 10) / -0.0000000000 0000000020 6437504424 783134d0/
      data sin cs( 11) /  0.0000000000 0000000000 0303969595 918706d0/
      data sin cs( 12) / -0.0000000000 0000000000 0000371357 734157d0/
      data sin cs( 13) /  0.0000000000 0000000000 0000000382 486123d0/
      data sin cs( 14) / -0.0000000000 0000000000 0000000000 336623d0/
      data sin cs( 15) /  0.0000000000 0000000000 0000000000 000256d0/
c
c pihi + pilo = pi.  pihi is exactly representable on all machines
c with at least 8 bits of precision.  whether it is exactly
c represented depends on the compiler.  this routine is more
c accurate if it is exactly represented.
      data pihi / 3.140625d0 /
      data pilo / 9.676535897 9323846264 3383279502 88d-4 /
      data pi2 / 1.5707963267 9489661923 1321691639 75d0 /
      data pirec / 0.3183098861 8379067153 7767526745 03d0 /
      data pi2rec / 0.6366197723 6758134307 5535053490 06d0 /
      data ntsn, xsml, xwarn, xmax / 0, 3*0.0d0 /
c
      if (ntsn.ne.0) go to 10
      ntsn = initds (sincs, 15, 0.1*sngl(d1mach(3)))
c
      xsml = dsqrt (2.0d0*d1mach(3))
      xmax = 1.0d0/d1mach(4)
      xwarn = dsqrt (xmax)
c
 10   absx = dabs (x)
      y = absx + pi2
      if (y.gt.xmax) call seteru (
     1  42hdcos    no precision because abs(x) is big, 42, 2, 2)
      if (y.gt.xwarn) call seteru (
     1  54hdcos    answer lt half precision because abs(x) is big,
     2  54, 1, 1)
c
      dcos = 1.0d0
      if (absx.lt.xsml) return
c
      xn = dint (y*pirec+0.5d0)
      n2 = dmod (xn, 2.0d0) + 0.5d0
      xn = xn - 0.5d0
      f = (absx-xn*pihi) - xn*pilo
c
      dcos = f + f*dcsevl (2.0d0*(f*pi2rec)**2-1.0d0, sincs, ntsn)
      if (n2.ne.0) dcos = -dcos
      if (dabs(dcos).gt.1.0d0) dcos = dsign (1.0d0, dcos)
c
      return
      end
      double precision function dcsevl (x, a, n)
c
c evaluate the n-term chebyshev series a at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
c
c             input arguments --
c x      dble prec value at which the series is to be evaluated.
c a      dble prec array of n terms of a chebyshev series.  in eval-
c        uating a, only half the first coef is summed.
c n      number of terms in array a.
c
      double precision a(n), x, twox, b0, b1, b2
c
      if (n.lt.1) call seteru (28hdcsevl  number of terms le 0, 28, 2,2)
      if (n.gt.1000) call seteru (31hdcsevl  number of terms gt 1000,
     1  31, 3, 2)
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call seteru (
     1  25hdcsevl  x outside (-1,+1), 25, 1, 1)
c
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
c
      dcsevl = 0.5d0 * (b0-b2)
c
      return
      end
      double precision function dint (x)
c december 1983 edition. w. fullerton, c3, los alamos scientific lab.
c
c dint is the double precision equivalent of aint.  this portable
c version is quite efficient when the argument is reasonably small (a
c common case), and so no faster machine-dependent version is needed.
c
      double precision x, xscl, scale, xbig, xmax, part, d1mach,
     1  dlog
      external d1mach, dlog, i1mach, r1mach
      data npart, scale, xbig, xmax / 0, 3*0.0d0 /
c
      if (npart.ne.0) go to 10
      ibase = i1mach(10)
      xmax = 1.0d0/d1mach (4)
      xbig = amin1 (float (i1mach(9)), 1.0/r1mach(4))
      scale = ibase**int(dlog(xbig)/dlog(dble(float(ibase)))-0.5d0)
      npart = dlog(xmax)/dlog(scale) + 1.0d0
c
 10   if (x.lt.(-xbig) .or. x.gt.xbig) go to 20
c
      dint = int(sngl(x))
      return
c
 20   xscl = dabs(x)
      if (xscl.gt.xmax) go to 50
c
      do 30 i=1,npart
        xscl = xscl/scale
 30   continue
c
      dint = 0.0d0
      do 40 i=1,npart
        xscl = xscl*scale
        ipart = xscl
        part = ipart
        xscl = xscl - part
        dint = dint*scale + part
 40   continue
c
      if (x.lt.0.0d0) dint = -dint
      return
c
 50   call seteru (68hdint    dabs(x) may be too big to be represented a
     1s an exact integer, 68, 1, 1)
      dint = x
      return
c
      end
      double precision function dlog (x)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, alncs(11), center(4), alncen(5), aln2, y, t,
     1  t2, xn,  dcsevl, d1mach
      external d1mach, dcsevl, initds
c
c series for aln        on the interval  0.          to  3.46021e-03
c                                        with weighted error   4.15e-32
c                                         log weighted error  31.38
c                               significant figures required  31.21
c                                    decimal places required  31.90
c
      data aln cs(  1) / +.1334719987 7973881561 6893860471 87 d+1     /
      data aln cs(  2) / +.6937562832 8411286281 3724383542 25 d-3     /
      data aln cs(  3) / +.4293403902 0450834506 5592108036 62 d-6     /
      data aln cs(  4) / +.2893384779 5432594580 4664403875 87 d-9     /
      data aln cs(  5) / +.2051251753 0340580901 7418134477 26 d-12    /
      data aln cs(  6) / +.1503971705 5497386574 6151533199 99 d-15    /
      data aln cs(  7) / +.1129454069 5636464284 5216133333 33 d-18    /
      data aln cs(  8) / +.8635578867 1171868881 9466666666 66 d-22    /
      data aln cs(  9) / +.6695299053 4350370613 3333333333 33 d-25    /
      data aln cs( 10) / +.5249155744 8151466666 6666666666 66 d-28    /
      data aln cs( 11) / +.4153054068 0362666666 6666666666 66 d-31    /
c
      data center(1) / 1.0d0 /
      data center(2) / 1.25d0 /
      data center(3) / 1.50d0 /
      data center(4) / 1.75d0 /
c
      data alncen(  1) / 0.0d0                                         /
      data alncen(  2) / +.2231435513 1420975576 6295090309 83 d+0     /
      data alncen(  3) / +.4054651081 0816438197 8013115464 34 d+0     /
      data alncen(  4) / +.5596157879 3542268627 0888500526 82 d+0     /
      data alncen(  5) / +.6931471805 5994530941 7232121458 17 d+0     /
c
c aln2 = alog(2.0) - 0.625
      data aln2 / 0.0681471805 5994530941 7232121458 18d0 /
      data nterms / 0 /
c
      if (nterms.eq.0) nterms = initds (alncs, 11, 28.9*sngl(d1mach(3)))
c
      if (x.le.0.d0) call seteru (
     1  29hdlog    x is zero or negative, 29, 1, 2)
c
      call d9upak (x, y, n)
c
      xn = n - 1
      y = 2.0d0*y
      ntrval = 4.0d0*y - 2.5d0
c
      if (ntrval.eq.5) t = ((y-1.0d0)-1.0d0) / (y+2.0d0)
      if (ntrval.lt.5) t = (y-center(ntrval)) / (y+center(ntrval))
      t2 = t*t
      dlog = 0.625d0*xn + (aln2*xn + alncen(ntrval) + 2.0d0*t +
     1  t*t2*dcsevl(578.d0*t2-1.0d0, alncs, nterms) )
c
      return
      end
      double precision function dsqrt (x)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, sqrt2(3), y,  d9pak, d1mach
      external alog, d1mach, d9pak
      data sqrt2(1) / 0.7071067811 8654752440 0844362104 85 d0 /
      data sqrt2(2) / 1.0 d0 /
      data sqrt2(3) / 1.4142135623 7309504880 1688724209 70 d0 /
c
      data niter / 0 /
c
      if (niter.eq.0) niter = 1.443*alog(-0.104*alog(0.1*sngl(d1mach(3))
     1  )) + 1.0
c
      if (x.le.0.d0) go to 20
c
      call d9upak (x, y, n)
      ixpnt = n/2
      irem = n - 2*ixpnt + 2
c
c the approximation below has accuracy of 4.16 digits.
      z = y
      dsqrt = .261599e0 + z*(1.114292e0 + z*(-.516888e0 + z*.141067e0))
c
      do 10 iter=1,niter
        dsqrt = dsqrt + 0.5d0*(y - dsqrt*dsqrt) / dsqrt
 10   continue
c
      dsqrt = d9pak (sqrt2(irem)*dsqrt, ixpnt)
      return
c
 20   if (x.lt.0.d0) call seteru (21hdsqrt   x is negative, 21, 1, 1)
      dsqrt = 0.0d0
      return
c
      end
      subroutine e9rint(messg,nw,nerr,save)
c
c  this routine stores the current error message or prints the old one,
c  if any, depending on whether or not save = .true. .
c
      integer messg(nw)
      logical save
      external i1mach, i8save
c
c  messgp stores at least the first 72 characters of the previous
c  message. its length is machine dependent and must be at least
c
c       1 + 71/(the number of characters stored per integer word).
c
      integer messgp(36),fmt(14),ccplus
c
c  start with no previous message.
c
      data messgp(1)/1h1/, nwp/0/, nerrp/0/
c
c  set up the format for printing the error message.
c  the format is simply (a1,14x,72axx) where xx=i1mach(6) is the
c  number of characters stored per integer word.
c
      data ccplus  / 1h+ /
c
      data fmt( 1) / 1h( /
      data fmt( 2) / 1ha /
      data fmt( 3) / 1h1 /
      data fmt( 4) / 1h, /
      data fmt( 5) / 1h1 /
      data fmt( 6) / 1h4 /
      data fmt( 7) / 1hx /
      data fmt( 8) / 1h, /
      data fmt( 9) / 1h7 /
      data fmt(10) / 1h2 /
      data fmt(11) / 1ha /
      data fmt(12) / 1hx /
      data fmt(13) / 1hx /
      data fmt(14) / 1h) /
c
      if (.not.save) go to 20
c
c  save the message.
c
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
c
        go to 30
c
 20   if (i8save(1,0,.false.).eq.0) go to 30
c
c  print the message.
c
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(7h error ,i4,4h in )
c
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
c
 30   return
c
      end
      subroutine eprint
c
c  this subroutine prints the last error message, if any.
c
      integer messg(1)
c
      call e9rint(messg,1,1,.false.)
      return
c
      end
      integer function i8save(isw,ivalue,set)
c
c  if (isw = 1) i8save returns the current error number and
c               sets it to ivalue if set = .true. .
c
c  if (isw = 2) i8save returns the current recovery switch and
c               sets it to ivalue if set = .true. .
c
      logical set
c
      integer iparam(2)
c  iparam(1) is the error number and iparam(2) is the recovery switch.
c
c  start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end
      function initds (dos, nos, eta)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c initialize the double precision orthogonal series dos so that initds
c is the number of terms needed to insure the error is no larger than
c eta.  ordinarily eta will be chosen to be one-tenth machine precision.
c
c             input arguments --
c dos    dble prec array of nos coefficients in an orthogonal series.
c nos    number of coefficients in dos.
c eta    requested accuracy of series.
c
      double precision dos(nos)
c
      if (nos.lt.1) call seteru (
     1  35hinitds  number of coefficients lt 1, 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru (28hinitds  eta may be too small, 28,
     1  1, 2)
      initds = i
c
      return
      end
      function inits (os, nos, eta)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c initialize the orthogonal series so that inits is the number of terms
c needed to insure the error is no larger than eta.  ordinarily, eta
c will be chosen to be one-tenth machine precision.
c
c             input arguments --
c os     array of nos coefficients in an orthogonal series.
c nos    number of coefficients in os.
c eta    requested accuracy of series.
c
      dimension os(nos)
c
      if (nos.lt.1) call seteru (
     1  35hinits   number of coefficients lt 1, 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru (28hinits   eta may be too small, 28,
     1  1, 2)
      inits = i
c
      return
      end
      subroutine r9upak (x, y, n)
c august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c unpack floating point number x so that x = y * 2.0**n, where
c 0.5 .le. abs(y) .lt. 1.0 .
c
      absx = abs(x)
      n = 0
      y = 0.0
      if (x.eq.0.0) return
c
 10   if (absx.ge.0.5) go to 20
      n = n - 1
      absx = absx*2.0
      go to 10
c
 20   if (absx.lt.1.0) go to 30
      n = n + 1
      absx = absx*0.5
      go to 20
c
 30   y = sign (absx, x)
      return
c
      end
      subroutine s88fmt( n, w, ifmt )
c
c  s88fmt  replaces ifmt(1), ... , ifmt(n) with
c  the characters corresponding to the n least significant
c  digits of w.
c
      integer n,w,ifmt(n)
c
      integer nt,wt,digits(10)
c
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
c
      end
      subroutine seterr (messg, nmessg, nerr, iopt)
c
c  this version modified by w. fullerton to dump if iopt = 1 and
c  not recovering.
c  seterr sets lerror = nerr, optionally prints the message and dumps
c  according to the following rules...
c
c    if iopt = 1 and recovering      - just remember the error.
c    if iopt = 1 and not recovering  - print, dump and stop.
c    if iopt = 2                     - print, dump and stop.
c
c  input
c
c    messg  - the error message.
c    nmessg - the length of the message, in characters.
c    nerr   - the error number. must have nerr non-zero.
c    iopt   - the option. must have iopt=1 or 2.
c
c  error states -
c
c    1 - message length not positive.
c    2 - cannot have nerr=0.
c    3 - an unrecovered error followed by another error.
c    4 - bad value for iopt.
c
c  only the first 72 characters of the message are printed.
c
c  the error handler calls a subroutine named fdump to produce a
c  symbolic dump. to complete the package, a dummy version of fdump
c  is supplied, but it should be replaced by a locally written version
c  which at least gives a trace-back.
c
      integer messg(1)
      external i1mach, i8save
c
c  the unit for error messages.
c
      iwunit=i1mach(4)
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      if (i8save(2,0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
 60   call fdump
      stop
c
      end
      subroutine seteru (messg, nmessg, nerr, iopt)
      common /cseter/ iunflo
      integer messg(1)
c      data iunflo / 0 /
c
      iunflo = 0
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
c
      if (iunflo.le.0) return
      call seterr (messg, nmessg, nerr, 1)
c
      return
      end



ccc
      parameter ( nc31 = 31, maxc1 = 100, maxa1 = 200,
     $     max2 = 3, nccl3=72 )
      character*(nccl3) dskout, cdfin, cdfout
cc      character*(maxct) date0,time0, datem,timem,timev,datev
cc      common / datim0 / seps, cdfin, cdfout, date0,time0,
cc     $     datem,timem, datev,timev
ccc
c      integer iargc, lena, counta
c      external iargc, getarg
c      character*72 namea, cdfin
ccc
      parameter ( numvar=100 )
      parameter ( ndima=2, nd1=ntsin,nd2=nsf, nd12=nd1*nd2 )
      parameter ( ndim0=5 )
      integer cdfid
      integer vid, vtype, dimsiz(ndim0), vn(numvar), vdims(ndima),
     $     vvdims(numvar,ndima), vnatt, vrtype(numvar), vrnat(numvar)
      integer rcode
      integer attlen, attype
      integer recid, nrecs
      integer vindx(ndima), v1ndx
      integer start(ndima), count(ndima)
      character t1char
      character*(maxc1) ctitle
      character*(maxa1) astrng
ccc   31 in the following statement id parameter maxncnam.
      character*(nc31) vname, vrname(numvar), recnam, attnam
ccc
      dimension vals(nd1,nd2), val0(nd12), nw(5)
      common / cdf1 / cdfid, rcode, recid, nrecs, attlen, attype
      common / cdf2 / vid, vtype, dimsiz, vn, vdims,
     $      vvdims, vnatt
      common / cdf3 / vindex, v1dex, start, count
      common / cdf4 / vals, val0, nw
      common / cdfc / t1char, ctitle, astrng, vname, recnam, attnam

ccc

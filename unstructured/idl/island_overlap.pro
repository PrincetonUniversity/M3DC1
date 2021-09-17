function island_overlap, filename, psin=psin, current=cur, $,
                         netcdf=netcdf, sum_files=sum_files, $
                         plot=pl, ntor=ntor, _EXTRA=extra

   width = island_widths(filename, psin=psin,current=cur,sum_files=sum_files,$
                         q=q, netcdf=netcdf, plot=pl, ntor=ntor, _EXTRA=extra)

   if(width eq 0) then return, 0

   n = n_elements(psin[*,0])
   overlap = fltarr(n)
   

   for j=0, n-1 do begin
       lpl = max(psin[j,*])
       if(psin[j,0] gt psin[j,n_elements(psin[j,*])-1]) then begin
          istart = 0
          iend = n_elements(psin[j,*])-1
          istep = 1
       endif else begin
          istart = n_elements(psin[j,*])-1
          iend = 0
          istep = -1
       end
       for i=istart, iend, istep do begin
           pr = psin[j,i] + width[j,i]/2.
           pl = psin[j,i] - width[j,i]/2.

           if(pr lt lpl) then break
           lpl = pl
       end
       overlap[j] = 1.-lpl

       if(keyword_set(pl)) then begin
          oplot, [1.,1.]*lpl, !y.crange, linestyle=2
       end
   end

   return, overlap
end

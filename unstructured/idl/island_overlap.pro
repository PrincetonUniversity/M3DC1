function island_overlap, filename, current=cur, sum_files=sum_files
   width = island_widths(filename, psin=psin,current=cur,sum_files=sum_files,$
                        q=q)

   n = n_elements(psin[*,0])
   overlap = fltarr(n)

   for j=0, n-1 do begin
       lpl = 1.
       for i=n_elements(psin[j,*])-1, 0, -1 do begin
           pr = psin[j,i] + width[j,i]/2.
           pl = psin[j,i] - width[j,i]/2.

           if(pr lt lpl) then break
           lpl = pl
       end
       overlap[j] = 1.-lpl
   end

   return, overlap
end

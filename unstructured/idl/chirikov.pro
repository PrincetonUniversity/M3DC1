function chirikov, filename, psin=psin, psimid=psimid, current=cur, $
                   sum_files=sum_files
   width = island_widths(filename, psin=psin, current=cur, sum_files=sum_files)

   chi = fltarr(n_elements(width[*,0]), n_elements(psin[0,*])-1)
   psimid = fltarr(n_elements(psin[*,0]), n_elements(psin[0,*])-1)

   for j=0, n_elements(width[0,*])-2 do begin
       chi[*,j] = (width[*,j+1]+width[*,j])/2. $
         / (psin[*,j+1]-psin[*,j])
       psimid[*,j] = (psin[*,j+1]+psin[*,j])/2.
   end

   return, chi
end


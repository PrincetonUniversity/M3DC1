function chirikov, filename, psin=psin, psimid=psimid, current=cur, $
                   sum_files=sum_files, netcdf=netcdf, q=q, ntor=ntor

   width = island_widths(filename, psin=psin, current=cur, ntor=ntor, $
                         sum_files=sum_files, netcdf=netcdf, q=q)

   if(width eq 0) then return, 0

   chi = fltarr(n_elements(width[*,0]), n_elements(psin[0,*])-1)
   psimid = fltarr(n_elements(psin[*,0]), n_elements(psin[0,*])-1)

   for j=0, n_elements(width[0,*])-2 do begin
       chi[*,j] = (width[*,j+1]+width[*,j])/2. $
         / abs(psin[*,j+1]-psin[*,j])
       psimid[*,j] = (psin[*,j+1]+psin[*,j])/2.
   end

   return, chi
end


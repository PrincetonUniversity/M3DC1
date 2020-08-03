function chi95, filename, psi0=psi0, current=cur, sum_files=sum_files, $
                netcdf=netcdf

   if(n_elements(psi0) eq 0) then psi0 = 0.95
   if(n_elements(psi0) lt n_elements(filename)) then $
     psi0 = replicate(psi0, n_elements(filename))

   chi = chirikov(filename, psin=psin, psimid=psimid, current=cur, $
                  sum_files=sum_files, netcdf=netcdf)

   m = n_elements(chi[0,*])
;   if(keyword_set(sum_files)) then n = 1 else n =
;   n_elements(filename)
   n = n_elements(chi[*,0])
   c95 = fltarr(n)
   for i=0, n-1 do begin
       c95[i] = interpol(chi[i,0:m-1], psimid[i,0:m-1], psi0[i])
   end

   return, c95
end

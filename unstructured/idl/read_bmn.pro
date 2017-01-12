function read_bmn, filename, m, bmn, phase, $
                   psin=psin, qval=q, qprime=qprime, area=area, $
                   psiprime=psiprime, sum_files=sum_files, factor=factor, $
                   netcdf=netcdf

   n = n_elements(filename)

   ; if sum_files is set, then sum data from files
   if(keyword_set(sum_files)) then begin
       result = read_bmn(filename, m, bmn0, phase0, $
         psin=psin0, qval=q0, qprime=qprime0, $
         area=area0, psiprime=psiprime0, $
         factor=factor, netcdf=netcdf)
       if(result ne 0) then return, result

       bmn = fltarr(1, n_elements(m))
       phase = fltarr(1, n_elements(m))

       temp = complexarr(n_elements(m))
       for i=0, n_elements(filename)-1 do begin
           temp = temp + bmn0[i,*]*exp(complex(0,1)*phase0[i,*])
       end
       bmn[0,*] = abs(temp)
       phase[0,*] = atan(imaginary(temp),real_part(temp))
       psin = psin0[0,*]
       q = q0[0,*]
       qprime = qprime0[0,*]
       area = area0[0,*]
       psiprime = psiprime0[0,*]
       return, 0
   end

   ; Read NetCDF file
   if(keyword_set(netcdf)) then begin
      for i=0, n-1 do begin
         read_bmncdf, file=filename[i], bmn=bmn0, psi=psi0, q=q0, $
                      ntor=ntor, m=m0, flux_pol=flux0, area=area0, bp=bp0

         ; calculate resonant fields
         if(q0[0] gt 0) then begin
            minm = fix(min(abs(q0))*ntor + 1)
            maxm = fix(max(abs(q0))*ntor)
         endif else begin
            maxm = -fix(min(abs(q0))*ntor + 1)
            minm = -fix(max(abs(q0))*ntor)
         end
         if(maxm gt max(m0)) then maxm = max(m0)
         m1 = findgen(maxm - minm + 1) + minm
         m = m1

         if(i eq 0) then begin
            bmn = fltarr(n, n_elements(m))
            phase = fltarr(n, n_elements(m))
            if(arg_present(psin)) then psin = fltarr(n, n_elements(m))
            if(arg_present(q)) then q = fltarr(n, n_elements(m))
            if(arg_present(qprime)) then qprime = fltarr(n, n_elements(m))
            if(arg_present(area)) then area = fltarr(n, n_elements(m))
            if(arg_present(psiprime)) then psiprime = fltarr(n, n_elements(m))
            if(n_elements(factor) eq 0) then factor = 1.
            if(n_elements(factor) lt n) then factor=replicate(factor[0],n)
         endif

         q1 = m1/ntor
         i0 = interpol(findgen(n_elements(q0)), q0, q1)
         j0 = intarr(n_elements(m1))
         for j=0, n_elements(m1)-1 do begin
            j0[j] = where(m0 eq m1[j], count)
            if(count eq 0) then print, $
               'Error: m = ', m1[j], ' not found!'
         end

         res_bmn = interpolate(bmn0, j0, i0)*factor[i]

         bmn[i,*] = abs(res_bmn)
         phase[i,*] = atan(imaginary(res_bmn), real_part(res_bmn))
         if(arg_present(psin)) then psin[i,*] = interpolate(psi0, i0)
         if(arg_present(q)) then q[i,*] = q1
         if(arg_present(qprime)) then qprime[i,*] = deriv(psin[i,*], q1)
         if(arg_present(area) and n_elements(area0) gt 1) then $
            area[i,*] = interpolate(area0, i0)
         if(arg_present(psiprime)) then $
            psiprime[i,*] = deriv(psin[i,*], interpolate(flux0, i0))
      end

   ; Read ASCII file
   endif else begin
      k = 0
      for i=0, n-1 do begin
         if(file_test(filename[i]) eq 0) then begin
            print, filename[i], ' not found'
            continue
         end
         result = read_ascii(filename[i])
         m = result.field1[0,*]
         k = n_elements(m)
         if(k ne 0) then break
      end
      if(k eq 0) then return, 1
      
      bmn = fltarr(n, n_elements(m))
      phase = fltarr(n, n_elements(m))
      if(arg_present(psin)) then psin = fltarr(n, n_elements(m))
      if(arg_present(q)) then q = fltarr(n, n_elements(m))
      if(arg_present(qprime)) then qprime = fltarr(n, n_elements(m))
      if(arg_present(area)) then area = fltarr(n, n_elements(m))
      if(arg_present(psiprime)) then psiprime = fltarr(n, n_elements(m))
      if(n_elements(factor) eq 0) then factor = 1.
      if(n_elements(factor) lt n) then factor=replicate(factor[0],n)
      
      for i=0, n-1 do begin
         if(file_test(filename[i]) eq 0) then begin
            print, 'Error reading ', filename[i]
            continue
         end
         result = read_ascii(filename[i])
         
         m = result.field1[0,*]
         temp = result.field1[1,*]*exp(complex(0,1)*result.field1[2,*])*factor[i]
         bmn[i,*] = abs(temp)
         phase[i,*] = atan(imaginary(temp),real_part(temp))
         if(arg_present(psin)) then psin[i,*] = result.field1[3,*]
         if(arg_present(q)) then q[i,*] = result.field1[4,*]
         if(arg_present(qprime)) then qprime[i,*] = result.field1[5,*]
         if(arg_present(area)) then area[i,*] = result.field1[6,*]
         if(arg_present(psiprime)) then psiprime[i,*] = result.field1[7,*]
      end
   endelse

   return, 0
end

function select_m, m0, bmn0, phase0, m, bmn, phase
   j = where(m0 eq m, count)
   if(count eq 0) then begin
       return, 1
   endif

   bmn = bmn0[*,j]
   phase = phase0[*,j]
   return, 0
end

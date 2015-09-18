function read_bmn, filename, m, bmn, phase, $
                   psin=psin, qval=q, qprime=qprime, area=area, $
                   psiprime=psiprime, sum_files=sum_files, factor=factor

   n = n_elements(filename)

   if(keyword_set(sum_files)) then begin
       result = read_bmn(filename, m, bmn0, phase0, $
         psin=psin0, qval=q0, qprime=qprime0, $
         area=area0, psiprime=psiprime0, $
         factor=factor)
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

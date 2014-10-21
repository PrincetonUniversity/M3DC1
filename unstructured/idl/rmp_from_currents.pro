function fft_cur, cur, nt=ntor, pts=pts
   if(n_elements(pts) eq 0) then pts=100
   n = pts*n_elements(cur)
   c = fltarr(n)
   ntor = fltarr(n)

   for i=0, n_elements(cur)-1 do begin
       c[i*pts:(i+1)*pts-1] = cur[i]
   end

   result = fft(c)
   result = shift(result, n/2)
   ntor = findgen(n) - n/2

   return, result
end


pro rmp_from_currents, iu=iu_cur, il=il_cur, c=c_cur, ntor=ntor, dir=dir
   max_ntor = 4
   phi = findgen(100)*360./100.
   amp = fltarr(100)
   if(n_elements(ntor) eq 0) then ntor = [1,2,3,4]
   if(n_elements(dir) eq 0) then dir='.'
   ii = complex(0,1)

   for in=0, n_elements(ntor)-1 do begin
       nstr = string(format='(I0)',ntor[in])
       openw, if_curr, dir+'/rmp_current.dat.'+nstr, /get_lun
       openw, if_coil, dir+'/rmp_coil.dat.'+nstr, /get_lun

       plot, phi, [-5,5], /nodata
       if(n_elements(il_cur) gt 0) then begin
           c = fft_cur(il_cur,nt=nt)
           i = where(nt eq ntor[in], count)
           a = 2.*abs(c[i])/1000.
           b = atan(imaginary(c[i]), real_part(c[i]))*180./!pi
           print, 'IL: ', a, b
           printf, if_coil, 2.161, -1.005, 0., 0., format='(4F12.4)'
           printf, if_coil, 2.372, -0.508, 0., 0., format='(4F12.4)'
           printf, if_curr, a, b, format='(2F12.4)'
           printf, if_curr, -a, b, format='(2F12.4)'
           for i=0, 99 do $
             amp[i] = real_part(a*exp(ii*(ntor[in]*phi[i]+b)*!pi/180.))
           oplot, phi, amp
       end
       if(n_elements(iu_cur) gt 0) then begin
           c = fft_cur(iu_cur,nt=nt)
           i = where(nt eq ntor[in], count)
           a = 2.*abs(c[i])/1000.
           b = atan(imaginary(c[i]), real_part(c[i]))*180./!pi
           print, 'IU: ', a, b
           printf, if_coil, 2.372, 0.508, 0., 0., format='(4F12.4)'
           printf, if_coil, 2.161, 1.005, 0., 0., format='(4F12.4)'
           printf, if_curr, a, b, format='(2F12.4)'
           printf, if_curr, -a, b, format='(2F12.4)'
           for i=0, 99 do $
             amp[i] = real_part(a*exp(ii*(ntor[in]*phi[i]+b)*!pi/180.))
           oplot, phi, amp, linestyle=1
       end
       if(n_elements(c_cur) gt 0) then begin
           c = fft_cur(c_cur,nt=nt)
           i = where(nt eq ntor[in], count)
           ; C-coil uses opposite sign convention and has 4 turns
           ; so the magnitude here is multiplied by -4.
           a = 2.*abs(c[i])/1000.  * (-4.)
           b = atan(imaginary(c[i]), real_part(c[i]))*180./!pi - 49.*ntor[in]
           print, 'C:', a, b
           printf, if_coil, 3.23, -0.8, 0., 0., format='(4F12.4)'
           printf, if_coil, 3.23,  0.8, 0., 0., format='(4F12.4)'
           printf, if_curr, a, b, format='(2F12.4)'
           printf, if_curr, -a, b, format='(2F12.4)'
           for i=0, 99 do $
             amp[i] = real_part(a*exp(ii*(ntor[in]*phi[i]+b)*!pi/180.))
           oplot, phi, amp, linestyle=2
       end

       free_lun, if_curr
       free_lun, if_coil
   end

end


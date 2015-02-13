function read_coil_data, directory=dir
   if(n_elements(dir) eq 0) then dir='.'

   coil = read_ascii(dir+'/coil.dat')
   curr = read_ascii(dir+'/current.dat')

   n = n_elements(coil.field1[0,*])
   m = n_elements(size(curr.field1, /dim))
   print, 'n, m = ', n, m

   f = fltarr(10,n)

   f[0,*] = curr.field1[0,*]
   if(m ge 2) then begin
       f[1,*] = curr.field1[1,*]
   endif else f[1,*] = curr.field1[1,*]

   l = n_elements(coil.field1[*,0])
   f[2:l+1,*] = coil.field1[0:l-1,*]
   print, f

   return, f
end

pro plot_coils, filename=file, directory=dir, overplot=overplot, _EXTRA=extra
   ct3
   
   if(n_elements(file) eq 1 and n_elements(dir) eq 0) then begin
       result = str_sep(file, '/')
       dir = ''
       for i=0, n_elements(result)-2 do begin
           dir = dir + result[i] + '/'
       end
       if(strlen(dir) eq 0) then dir = '.'
   end
   if(n_elements(dir) eq 0) then dir = '.'

   dat = read_coil_data(directory=dir)
   n = n_elements(dat[0,*])


   if(not keyword_set(overplot)) then begin
       xrange = [min(dat[2,*]-dat[4,*]), max(dat[2,*]+dat[4,*])]
       zrange = [min(dat[3,*]-dat[5,*]), max(dat[3,*]+dat[5,*])]
       plot, xrange, zrange, /nodata, _EXTRA=extra
   end
   for i=0, n-1 do begin
       m = dat[8,i]*dat[9,i]

       a1 = dat[6,i]*!pi/180.
       a2 = dat[7,i]*!pi/180.
       if(a2 ne 0.) then a2 = a2 + !pi/2.
       a1 = tan(a1)
       a2 = tan(a2)

       xp = [dat[2,i]-(dat[4,i]-dat[5,i]*a2)/2., $
             dat[2,i]+(dat[4,i]+dat[5,i]*a2)/2., $
             dat[2,i]+(dat[4,i]-dat[5,i]*a2)/2., $
             dat[2,i]-(dat[4,i]+dat[5,i]*a2)/2.]

       zp = [dat[3,i]-(dat[5,i]+dat[4,i]*a1)/2., $
             dat[3,i]-(dat[5,i]-dat[4,i]*a1)/2., $
             dat[3,i]+(dat[5,i]+dat[4,i]*a1)/2., $
             dat[3,i]+(dat[5,i]-dat[4,i]*a1)/2.]

       oplot, [xp, xp[0]], [zp, zp[0]], color=color(3)
;       continue
       
       if(m gt 1) then begin
           xc = fltarr(m)
           zc = fltarr(m)

           s = 0
           for j=1, dat[8,i] do begin
               for k=1, dat[9,i] do begin
                   if(dat[8,i] eq 1) then begin
                       dx = 0.
                   endif else begin
                       dx = dat[4,i]*(j-1.)/(dat[8,i]-1.) - dat[4,i]/2.
                   end
                   if(dat[9,i] eq 1) then begin
                       dz = 0.
                   endif else begin
                       dz = dat[5,i]*(k-1.)/(dat[9,i]-1.) - dat[5,i]/2.
                   end

                   xc[s] = dat[2,i] + dx - dz*a2
                   zc[s] = dat[3,i] + dz + dx*a1
                   s = s + 1
               end
           end
           oplot, xc, zc, psym=3, color=color(6)
       end
   end

end

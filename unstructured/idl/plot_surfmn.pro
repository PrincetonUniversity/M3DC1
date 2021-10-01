pro plot_surfmn, file, srnorm=srnorm, reverse_q=reverse_q, fac=fac, $
                 straight=straight, ntor=nt, title=title, _EXTRA=extra

  if(n_elements(file) eq 0) then return

  if(n_elements(fac) eq 0) then begin
     fac=replicate(1.,n_elements(file))
     label ='!7d!8B!Dmn!N!6 (G/kA)!X'
  endif else begin
     label ='!7d!8B!Dmn!N!6 (G)!X'
     if(n_elements(fac) lt n_elements(file)) then begin
        fac = replicate(fac, n_elements(file))
     end
  endelse

  for i=0, n_elements(file)-1 do begin
     read_bmncdf, file=file[i], _EXTRA=extra, bmn=bmn0, psi=psi, m=m0, q=q, $
                  ntor=ntor0
     help, bmn0
     if(i eq 0) then begin
        bmn = bmn0*fac[i]
        m = m0
        ntor = ntor0
     endif else begin
        if(ntor ne ntor0) then begin
           print, 'Error: files have different ntor'
           return
        end
        if(not array_equal(m, m0)) then begin
           print, 'Error: files have different m'
           print, m
           print, m0
           return
        end
        bmn = bmn + bmn0*fac[i]
     endelse
  end

  if(keyword_set(reverse_q)) then q = -q

   if(keyword_set(srnorm)) then begin
      y = sqrt(psi)
      ytitle='!9r!7W!X'
   endif else begin
      y = psi
      ytitle='!7W!X'
   end

   if(1 eq strcmp(!d.name, 'PS', /fold_case)) then begin
       xsize = 1.33
   endif else begin
       xsize = 1.
   endelse

   qmax = max(q,/abs)

   k=0
   if(n_elements(ntor) gt 1) then begin
      k = where(ntor eq nt, ct)
      if(ct ne 1) then $
         print, 'Error: ntor = ', nt, ' not found.'
   endif else begin
      k = 0
      nt = ntor
   endelse

   if(n_elements(title) eq 0) then begin
      title = '!8n!6 = ' + string(format='(I0)',nt) + '!X'
   end
   
   if(keyword_set(straight)) then begin

      xtitle='!8m!6/!8nq!X'

      mu = m/(nt*qmax)
      for i=0, n_elements(q)-1 do begin
         bmn[*,*,i] = interpolate(reform(bmn[*,*,i]), mu*nt*q[i]-m[0]) 
      end
      
      contour_and_legend, abs(reform(bmn[k,*,*])), mu, y,  $
                          table=39, xtitle=xtitle, ytitle=ytitle, $
                          xrange=[-4,4], yrange=[0,1], /lines, c_thick=1, $
                          ccolor=!d.table_size-1, label=label, $
                          _EXTRA=extra, xsize=xsize, title=title
      oplot, [1,1], !y.crange, linestyle=2, color=!d.table_size-1
   endif else begin
      xtitle='!8m!X'
      contour_and_legend, abs(reform(bmn[k,*,*])), m, y,  $
                          table=39, xtitle=xtitle, ytitle=ytitle, $
                          xrange=[-4,4]*abs(qmax), yrange=[0,1], $
                          /lines, c_thick=1, $
                          ccolor=!d.table_size-1, label=label, $
                          _EXTRA=extra, xsize=xsize, title=title
      oplot, nt*q, y, linestyle=2, color=!d.table_size-1
   end

end

pro plot_bmn, filename, vac=vac, names=names, $
              ytitle=ytitle, width=width, current=cur, $
              bmncdf=bmncdf, chirikov=chi, sum_files=sum_files, $
              color=c, overplot=overplot, $
              _EXTRA=extra

   if(n_elements(cur) eq 0) then cur=1.
   if(n_elements(cur) lt n_elements(filename)) then $
     cur = replicate(cur, n_elements(filename))

   xtitle = '!7W!X'

   if(keyword_set(width)) then begin
       bmn = island_widths(filename, psin=psin, cur=cur, sum_files=sum_files)
       if(n_elements(ytitle) eq 0) then $
         ytitle='!6Island Width (!7W!6)!X'
   endif else if(keyword_set(chi)) then begin
       bmn = chirikov(filename,cur=cur,psimid=psin, sum_files=sum_files)
       if(n_elements(ytitle) eq 0) then $
         ytitle='!6Chirikov Parameter!X'
   endif else begin
       if(keyword_set(bmncdf)) then begin
           print, 'reading bmncdf file'
           for i=0, n_elements(filename)-1 do begin
               read_bmncdf, file=filename[i], bmn=bmn0, psi=psi, m=m, q=q, $
                 ntor=ntor
               if(i eq 0) then begin
                   minm = fix(q[0] * ntor) 
                   maxm = fix(q[n_elements(q)-1] * ntor)
                   if(maxm gt max(m)) then maxm = max(m)
                   m0 = findgen(maxm - minm) + minm
                   q0 = m0/ntor
                   print, 'm0 = ', m0
                   print, 'q0 = ', q0
                   i0 = interpol(findgen(n_elements(q)), q, q0)
                   print, 'q = ', interpolate(q, i0)
                   print, 'psi = ', interpolate(psi, i0)
                   j0 = intarr(n_elements(m0))
                   for j=0, n_elements(m0)-1 do begin
                       j0[j] = where(m eq m0[j], count)
                       if(count eq 0) then print, $
                         'Error: m = ', m0[j], ' not found!'
                   end
                   bmn = fltarr(n_elements(filename), n_elements(i0))
                   psin = fltarr(n_elements(filename), n_elements(i0))
               end
               bmn[i,*] = interpolate(abs(bmn0), j0, i0)
               psin[i,*] = interpolate(psi, i0)
           end
           if(keyword_set(sum_files)) then begin
               for i=1, n_elements(bmn[*,0])-1 do begin
                   bmn[0,*] = bmn[0,*] + bmn[i,*]
               end
               bmn = bmn[0,*]
           end
       endif else begin
           result = read_bmn(filename,m,bmn,phase,psin=psin,$
                             sum_files=sum_files,factor=cur)
       endelse
      
       if(n_elements(ytitle) eq 0) then $
         ytitle='!6Total Resonant Field (G/kA)!X'
   endelse

   ct3
   if(n_elements(c) eq 0) then begin
       if(n_elements(filename) gt 1) then begin
           c = shift(get_colors(n_elements(filename)),-1)
       endif else begin
           c = get_colors()
        endelse
   end
   if(not keyword_set(overplot)) then begin
       plot, psin, bmn, /nodata, $
         xtitle=xtitle, ytitle=ytitle, _EXTRA=extra
   end
   for i=0, n_elements(bmn[*,0])-1 do begin
       oplot, psin[i,*], bmn[i,*], color=c[i]
       oplot, psin[i,*], bmn[i,*], color=c[i], psym=4
   end

   if(n_elements(vac) gt 0) then begin
       if(keyword_set(width)) then begin
           bmn_vac = island_widths(vac, psin=psin_vac, cur=cur, $
                                  sum_files=sum_files)
       endif else if(keyword_set(chi)) then begin
           bmn_vac = chirikov(vac,cur=cur,psimid=psin_vac)           
       endif else begin
           result = read_bmn(vac, m_vac, bmn_vac, phase_vac, psin=psin_vac, $
                            sum_files=sum_files)
       endelse
       if(n_elements(bmn_vac[*,0]) eq 1) then begin
           oplot, psin_vac[0,*], bmn_vac[0,*], linestyle=1
       endif else begin
           for i=0, n_elements(bmn_vac[*,0])-1 do begin
               oplot, psin_vac[i,*], bmn_vac[i,*], linestyle=1, color=c[i]
;           oplot, psin_vac[i,*], bmn_vac[i,*], linestyle=1, color=c[i], psym=5
           end
       end
   end

   if(n_elements(filename) gt 1 and not keyword_set(sum_files)) then begin
       if(n_elements(names) eq 0) then names = filename
       plot_legend, names, color=c, psym=replicate(4, n_elements(filename)), $
         _EXTRA=extra, linestyle=replicate(0,n_elements(filename))
   end
end

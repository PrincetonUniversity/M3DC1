pro plot_bmn, filename, vac=vac, names=names, $
              ytitle=ytitle, width=width, current=cur, $
              _EXTRA=extra

   if(n_elements(cur) eq 0) then cur=1.
   if(n_elements(cur) lt n_elements(filename)) then $
     cur = replicate(cur, n_elements(filename))

   xtitle = '!7W!X'

   if(keyword_set(width)) then begin
       bmn = island_widths(filename, psin=psin, cur=cur)
       if(n_elements(ytitle) eq 0) then $
         ytitle='!6Island Width (!7W!6)!X'
   endif else begin
       result = read_bmn(filename, m, bmn, phase, psin=psin)
       if(n_elements(ytitle) eq 0) then $
         ytitle='!6Total Resonant Field (G/kA)!X'
   endelse

   ct3
   if(n_elements(filename) gt 1) then begin
       c = shift(get_colors(),-1)
   endif else begin
       c = get_colors()
   endelse
   plot, psin, bmn, /nodata, $
     xtitle=xtitle, ytitle=ytitle, _EXTRA=extra
   for i=0, n_elements(filename)-1 do begin
       oplot, psin[i,*], bmn[i,*], color=c[i]
       oplot, psin[i,*], bmn[i,*], color=c[i], psym=4
   end

   if(n_elements(vac) gt 0) then begin
       if(keyword_set(width)) then begin
           bmn_vac = island_widths(vac, psin=psin_vac, cur=cur)
       endif else begin
           result = read_bmn(vac, m_vac, bmn_vac, phase_vac, psin=psin_vac)
       endelse
       if(n_elements(vac) eq 1) then begin
           oplot, psin_vac[0,*], bmn_vac[0,*], linestyle=1
       endif else begin
           for i=0, n_elements(vac)-1 do begin
               oplot, psin_vac[i,*], bmn_vac[i,*], linestyle=1, color=c[i]
;           oplot, psin_vac[i,*], bmn_vac[i,*], linestyle=1, color=c[i], psym=5
           end
       end
   end

   if(n_elements(filename) gt 1) then begin
       if(n_elements(names) eq 0) then names = filename
       plot_legend, names, color=c, psym=replicate(4, n_elements(filename)), $
         _EXTRA=extra
   end
end

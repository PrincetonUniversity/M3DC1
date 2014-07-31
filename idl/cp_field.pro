pro plot_cp_field, field, filename=filename, _EXTRA=extra, names=names, $
                   plasma=plasma, vacuum=vacuum, abs=abs, phase=phase

   if(n_elements(field) eq 0) then field = 'bx'
   if(n_elements(filename) eq 0) then filename = 'C1.h5'
   if(n_elements(names) eq 0) then names = filename

   if(not (keyword_set(abs) or keyword_set(phase))) then begin
       filename=filename[0]
   end

; from ntor=1
;   fac = 0.95493

   c = get_colors()
   if(n_elements(filename) gt 1) then c = shift(c, -1)

   for fn=0, n_elements(filename)-1 do begin
       ntor = read_parameter('ntor',filename=filename[fn],_EXTRA=extra)

       xtitle = '!6Toroidal Angle (deg)!X'
       ytitle = '!8Z!6 (m)!X'

       if(not keyword_set(vacuum)) then begin
           bx = read_field(field,x,y,t,_EXTRA=extra,slice=1,/complex,/linear,$
                           filename=filename[fn], symbol=symbol)
       end

       if(keyword_set(plasma) or keyword_set(vacuum)) then begin
           bx_vac = read_field(field,x,y,t,_EXTRA=extra,slice=0,/complex,/linear,$
                               filename=filename[fn],symbol=symbol)
           
           if(keyword_set(plasma)) then begin
               title='!6Field from Plasma!X'
               db = reform(bx - bx_vac)
           end
           if(keyword_set(vacuum)) then begin
               title='!6Vacuum Field!X'
               db = reform(bx_vac)
           end
       endif else begin
           title='!6Total Perturbed Field!X'
           db = reform(bx)
       endelse
       
       label = symbol
       
                                ; convert to Gauss
       db = db*1e4
       
       nz = 100
       nt = 100
       r0 = 0.98
       r1 = 0.98
       z0 = -1.
       z1 = 1
       
       phi = 2.*!pi*findgen(nt)/nt
       z = (z1 - z0)*findgen(nz)/(nz-1.) + z0
       r = (r1 - r0)*findgen(nz)/(nz-1.) + r0
       
       ir = n_elements(x) * (r - min(x)) / (max(x) - min(x))
       iz = n_elements(y) * (z - min(y)) / (max(y) - min(y))
       
       dbi = interpolate(db, ir, iz) ;*fac

       if(keyword_set(abs)) then begin
           f = abs(dbi)
           if(fn eq 0) then begin
               plot, z, f, /nodata, _EXTRA=extra, xtitle=ytitle, $
                 ytitle=label, title=title, subtitle=subtitle
           endif
           oplot, z, f, color=c[fn]
       endif else if(keyword_set(phase)) then begin
           ; minus sign is to convert from M3D-C1 angle to DIII-D angle
           f = -atan(imaginary(dbi), real_part(dbi))*180./!pi
           if(fn eq 0) then begin
               plot, z, f, xtitle=ytitle, yrange=[-180,180], $
                 ytitle='!6Toroidal Phase (deg)!X', $
                 title=title, subtitle=subtitle, _EXTRA=extra
           endif
           oplot, z, f, color=c[fn]
       endif else begin
           db_cp = fltarr(1,nt,nz)
           for i=0, nt-1 do begin
               db_cp[0,i,*] = real_part(dbi*exp(complex(0,1)*ntor*phi[i]))
           end
           
           contour_and_legend, db_cp, phi*180./!pi, z, $
             title=title, xtitle=xtitle, ytitle=ytitle, $
             label=label, subtitle=subtitle, _EXTRA=extra
       endelse
   end

   if(n_elements(names) gt 1) then begin
       plot_legend, names, color=c, _EXTRA=extra
   end
end

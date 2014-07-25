pro plot_field, name, time, x, y, points=p, mesh=plotmesh, $
                mcolor=mc, lcfs=lcfs, title=title, units=units, $
                range=range, rrange=rrange, zrange=zrange, linear=linear, $
                xlim=xlim, cutx=cutx, cutz=cutz, mpeg=mpeg, linfac=linfac, $
                mask_val=mask_val, boundary=boundary, q_contours=q_contours, $
                overplot=overplot, phi=phi0, time=realtime, levels=levels, $
                phase=phase, abs=abs, operation=op, magcoord=magcoord, $
                outfile=outfile, fac=fac, _EXTRA=ex

   ; open mpeg object
   if(n_elements(mpeg) ne 0) then begin
       mpegid = mpeg_open([640,480],quality=100)

       for i=time[0],time[1] do begin
           plot_field, name, i, x, y, points=p, mesh=plotmesh, $
             mcolor=mc, lcfs=lcfs, title=title, units=units, $
             range=range, rrange=rrange, zrange=zrange, linear=linear, $
             xlim=xlim, cutx=cutx, cutz=cutz, linfac=linfac, levels=levels, $
             mask_val=mask_val, boundary=boundary, q_contours=q_contours, $
             overplot=overplot, phi=phi0, time=realtime, $
             phase=phase, abs=abs, operation=op, _EXTRA=ex

           image = tvrd(true=1)
               
           image[0,*,*] = rotate(reform(image[0,*,*]), 7)
           image[1,*,*] = rotate(reform(image[1,*,*]), 7)
           image[2,*,*] = rotate(reform(image[2,*,*]), 7)
               
           mpeg_put, mpegid, image=image, frame=(i-time[0])
       end

       print, 'Writing mpeg...'
       mpeg_save, mpegid, filename=mpeg
       mpeg_close, mpegid
       return
   end


   if(n_elements(time) eq 0) then time = 0
   if(n_elements(p) eq 0) then p = 200
   if(n_elements(title) eq 0) then notitle = 1 else notitle = 0

   if(keyword_set(phase) or keyword_set(abs)) then complex=1

   if(size(name, /type) eq 7) then begin
       field = read_field(name, x, y, t, slices=time, mesh=mesh, $
                          points=p, rrange=rrange, zrange=zrange, $
                          symbol=fieldname, units=u, linear=linear, $
                          mask=mask, phi=phi0, time=realtime, $
                          complex=complex, operation=op, $
                          linfac=linfac, fac=fac,_EXTRA=ex)
       if(n_elements(field) le 1) then return
       if(n_elements(units) eq 0) then units=u       
   endif else begin
       field = name
   endelse

   if(keyword_set(phase)) then begin
       field = atan(imaginary(field), real_part(field))*180./!pi
       units = '!6Degrees!X'
       fieldname = '!6Phase(!X' + fieldname + '!6)!X'
   endif else if(keyword_set(abs)) then begin
       field = abs(field)
   endif


   if(n_elements(realtime) ne 0) then print, 'time = ', realtime

   ; remove NaN's from result
   i = where(not float(finite(field)), count)
   if(count gt 0) then begin
       print, "removing NaN's"
       field[i] = 0.
   endif

   if(n_elements(mask_val) ne 0) then begin
       for k=0, n_elements(field[*,0,0])-1 do begin
           field[k,*,*] = field[k,*,*] - mask*(field[k,*,*] - mask_val)
       end
   endif

   field = real_part(field)

   if(n_elements(range) eq 0) then range = [min(field),max(field)]

   if((notitle eq 1) and (n_elements(t) ne 0)) then begin
       title = fieldname
   end

   if(n_elements(cutx) gt 0) then begin
       dum = min(x-cutx,i,/absolute)
       data = reform(field[0,i,*])
       if(keyword_set(overplot)) then begin
           oplot, y, data, _EXTRA=ex
       endif else plot, y, field[0,i,*], title=title, _EXTRA=ex
       if(n_elements(outfile) eq 1) then begin
           openw, ifile, outfile, /get_lun
           printf, ifile, format='(2E16.6)', transpose([[y], [data]])
           free_lun, ifile
       endif
   endif else if(n_elements(cutz) gt 0) then begin
       dum = min(y-cutz,i,/absolute)
       data = reform(field[0,*,i])
       if(keyword_set(overplot)) then begin
           oplot, x, data, _EXTRA=ex
       endif else plot, x, data, title=title, _EXTRA=ex
       if(n_elements(outfile) eq 1) then begin
           openw, ifile, outfile, /get_lun
           printf, ifile, format='(2E16.6)', transpose([[x], [data]])
           free_lun, ifile
       endif
   endif else if(keyword_set(magcoord)) then begin

       psi = read_field('psi',x,y,t,points=p,/equilibrium,_EXTRA=ex)
       field = flux_coord_field(field, psi, x, y, t, fbins=p, tbins=p, $
                                nflux=nflux, angle=angle, bins=p)
       
       contour_and_legend, transpose(field[0,*,*], [0,2,1]),$
         angle*180./!pi, nflux, title=title, $
         label=units, levels=levels, $
         xtitle='!6Angle (Degrees)!X', $
         ytitle='!7W!X', $
         range=range, _EXTRA=ex

       
   endif else begin          
       contour_and_legend, field[0,*,*], x, y, title=title, $
         label=units, levels=levels, $
         xtitle=make_label('!8R!X', /l0, _EXTRA=ex), $
         ytitle=make_label('!8Z!X', /l0, _EXTRA=ex), $
         range=range, overplot=overplot, _EXTRA=ex

       if(n_elements(q_contours) ne 0) then begin
           fval = flux_at_q(q_contours,points=p,_EXTRA=ex)
           plot_flux_contour, fval, points=p, closed=0, /overplot, $
             thick=!p.thick/2., _EXTRA=ex
       endif

       if(keyword_set(lcfs)) then begin
           print, 'passing slice = ', time[0]
           plot_lcfs, points=p, slice=time[0], $
             last=last, _EXTRA=ex
       endif

       if(keyword_set(boundary)) then plotmesh=1
       if(keyword_set(plotmesh)) then begin
           plot_mesh, mesh=mesh, /oplot, $
             boundary=boundary, _EXTRA=ex
       endif
   endelse
end

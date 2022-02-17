pro plot_field, name, time, x, y, points=p, mesh=plotmesh, $
                mcolor=mc, lcfs=lcfs, title=title, units=units, $
                range=range, rrange=rrange, zrange=zrange, linear=linear, $
                xlim=xlim, cutx=cutx, cutz=cutz, mpeg=mpeg, linfac=linfac, $
                mask_val=mask_val, boundary=boundary, q_contours=q_contours, $
                overplot=overplot, phi=phi0, time=realtime, levels=levels, $
                phase=phase, abs=abs, operation=op, magcoord=magcoord, $
                outfile=outfile, fac=fac, filename=filename, $
                psin=psin, coils=coils, axis=axis, wall_regions=plotwall, $
                _EXTRA=ex

   if(n_elements(time) eq 0) then time = 0
   if(n_elements(p) eq 0) then p = 200
   if(n_elements(title) eq 0) then notitle = 1 else notitle = 0

   if(keyword_set(phase) or keyword_set(abs)) then complex=1

   if(size(name, /type) eq 7) then begin
       field = read_field(name, x, y, t, slices=time, mesh=mesh, $
                          points=p, rrange=rrange, zrange=zrange, $
                          symbol=fieldname, units=u, linear=linear, $
                          mask=mask, phi=phi0, time=realtime, $
                          complex=complex, operation=op,filename=filename, $
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
       data = field_at_point(field, x, y, replicate(cutx, n_elements(x)), y)
       if(keyword_set(psin)) then begin
           psi = read_field('psi_norm',x,y,t,points=p,/equilibrium,$
                            filename=filename,_EXTRA=ex)
           yy = reform(psi[0,i,*])
       endif else begin
          xtitle=make_label('!8Z!X', /l0, _EXTRA=ex)
          yy = y
       end
       ytitle=units
       if(keyword_set(overplot)) then begin
           oplot, yy, data, _EXTRA=ex
       endif else plot, yy, data, title=title, xtitle=xtitle, ytitle=ytitle, _EXTRA=ex
       if(n_elements(outfile) eq 1) then begin
           openw, ifile, outfile, /get_lun
           printf, ifile, format='(2E16.6)', transpose([[yy], [reform(data)]])
           free_lun, ifile
       endif
   endif else if(n_elements(cutz) gt 0) then begin
       data = field_at_point(field, x, y, x, replicate(cutz, n_elements(y)))
       if(keyword_set(psin)) then begin
           psi = read_field('psi_norm',x,y,t,points=p,/equilibrium,$
                            filename=filename,_EXTRA=ex)
           xx = reform(psi[0,*,i])
       endif else begin
          xtitle=make_label('!8R!X', /l0, _EXTRA=ex)
          xx = x
       end
       ytitle=units
       if(keyword_set(overplot)) then begin
           oplot, xx, data, _EXTRA=ex
       endif else plot, xx, data, title=title, xtitle=xtitle, ytitle=ytitle, _EXTRA=ex
       if(n_elements(outfile) eq 1) then begin
           openw, ifile, outfile, /get_lun
           printf, ifile, format='(2E16.6)', transpose([[xx], [reform(data)]])
           free_lun, ifile
       endif
   endif else if(keyword_set(magcoord)) then begin

      psi = read_field('psi',x,y,t,points=p,/equilibrium,_EXTRA=ex)
      field = flux_coord_field(field, psi, x, y, t, fbins=p, tbins=p, $
                               nflux=nflux, angle=angle, bins=p, _EXTRA=ex)
      ;; field[0,*,*] = $
      ;;    flux_coord_field_new(field,x,y,/fast,filename=filename,points=p, $
      ;;                         fc=fc,_EXTRA=extra)
       
      contour_and_legend, transpose(field[0,*,*],[0,2,1]), $
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

       if(keyword_set(boundary)) then plotmesh=1
       if(keyword_set(plotmesh)) then begin
           plot_mesh, mesh=mesh, /oplot, $
             boundary=boundary, filename=filename[0], _EXTRA=ex
       endif
       if(keyword_set(plotwall)) then begin
          plot_wall_regions, /over, color=color(3)
       end
       
       if(n_elements(q_contours) ne 0) then begin
           fval = flux_at_q(q_contours,points=p,_EXTRA=ex,$
                           filename=filename[0])
           if(fval[0] ne 0) then begin
              plot_flux_contour, fval, points=p, closed=0, /overplot, $
                                 thick=!p.thick/2., filename=filename[0], $
                                 _EXTRA=ex
           end
       endif

       if(keyword_set(axis)) then begin
          nulls, axis=ax, xpoints=xpoint, _EXTRA=extra, filename=filename
          dx = (!x.crange[1]-!x.crange[0])/50.
          dy = (!y.crange[1]-!y.crange[0])/50.
          oplot, [ax[0]-dx, ax[0]+dx], $
                 [ax[1]-dy, ax[1]+dy], color=color(6,10)
          oplot, [ax[0]-dx, ax[0]+dx], $
                 [ax[1]+dy, ax[1]-dy], color=color(6,10)
       end

       if(keyword_set(lcfs)) then begin
           print, 'passing slice = ', time[0]
           plot_lcfs, points=p, slice=time[0], $
             last=last, filename=filename[0], _EXTRA=ex
       endif

       if(keyword_set(coils)) then begin
           plot_coils, _EXTRA=extra, filename=filename[0], /overplot
       end
   endelse
end

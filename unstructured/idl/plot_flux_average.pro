;======================================================
; plot_flux_average
; ~~~~~~~~~~~~~~~~~
;
; plots the flux average quantity "name" at a give time
;======================================================
pro plot_flux_average, field, time, filename=filename, complex=complex, $
                       color=colors, names=names, bins=bins, linear=linear, $
                       xlog=xlog, ylog=ylog, overplot=overplot, fac=fac, $
                       lcfs=lcfs, normalized_flux=norm, points=pts, $
                       linestyle=ls, linfac=linfac, sum=sum, $
                       minor_radius=minor_radius, smooth=sm, t=t, rms=rms, $
                       bw=bw, srnorm=srnorm, last=last, mks=mks, cgs=cgs, $
                       q_contours=q_contours, rho=rho, integrate=integrate, $
                       multiply_flux=multiply_flux, abs=abs, phase=phase, $
                       stotal=total, nolegend=nolegend, outfile=outfile, $
                       val_at_q=val_at_q, flux_at_q=qflux, $
                       regularize=regularize, _EXTRA=extra

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(linfac) eq 0) then linfac = 1.
   if(n_elements(fac) eq 0) then fac = 1.
   if(n_elements(ls) eq 0) then ls = 0

   if(n_elements(time) eq 0) then time=0
   if(keyword_set(last)) then $
     time = fix(read_parameter('ntime',filename=filename)-1)

   if(n_elements(multiply_flux) eq 0) then multiply_flux = 0.

   if(n_elements(norm) eq 0) then norm=1
   
   if(n_elements(field) gt 1) then begin
       if(keyword_set(bw)) then begin
           if(n_elements(ls) eq 1) then ls = indgen(nfiles)
           colors = replicate(color(0,1), n_elements(field))
       endif else begin
           if(n_elements(colors) eq 0) then colors = shift(get_colors(),-1)
           ls = replicate(0,n_elements(field))
       endelse
       if(n_elements(linfac) eq 1) then linfac=replicate(linfac, n_elements(field))
       for i=0, n_elements(field)-1 do begin
           if((n_elements(q_contours) ne 0) and (i eq 0)) then begin
               qcon = q_contours
           end
           plot_flux_average, field[i], time, filename=filename, $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=col[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, fac=fac, $
             linear=linear, multiply_flux=multiply_flux, mks=mks, cgs=cgs, $
             integrate=integrate, complex=complex, abs=abs, phase=phase, $
             stotal=total, q_contours=qcon, rho=rho, nolegend=nolegend, $
             linfac=linfac[i], regularize=regularize
       end
       if(n_elements(names) ne 0 and not keyword_set(nolengend)) then begin
           plot_legend, names, colors=col, linestyle=ls, _EXTRA=extra
       end
       return
   end

   nfiles = n_elements(filename)
   if(nfiles gt 1 and not keyword_set(sum)) then begin
       if(n_elements(names) eq 0) then names=filename
       if(keyword_set(bw)) then begin
           if(n_elements(ls) eq 1) then ls = indgen(nfiles)
           colors = replicate(color(0,1), nfiles)
       endif else begin
           if(n_elements(colors) eq 0) then colors = shift(get_colors(),-1)
           if(n_elements(ls) eq 0) then ls = replicate(0,nfiles)
       endelse
       if(n_elements(ls) eq 1) then ls = replicate(ls, nfiles)
       if(n_elements(time) eq 1) then time = replicate(time,nfiles)
       if(n_elements(linfac) eq 1) then linfac=replicate(linfac, nfiles)
       if(n_elements(fac) eq 1) then fac=replicate(fac, nfiles)
       if(n_elements(multiply_flux) eq 1) then $
         multiply_flux = replicate(multiply_flux,nfiles)
       if(n_elements(outfile) eq 0) then outfile = replicate('',nfiles)

       for i=0, nfiles-1 do begin
           newfield = field
           plot_flux_average, newfield, time[i], filename=filename[i], $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, fac=fac[i], $
             linear=linear, multiply_flux=multiply_flux[i], mks=mks, cgs=cgs, $
             integrate=integrate, complex=complex, abs=abs, phase=phase, $
             stotal=total, q_contours=q_contours, rho=rho, nolegend=nolegend, $
             linfac=linfac[i], out=outfile[i], regularize=regularize
       end
       if(n_elements(names) gt 0 and not keyword_set(nolegend)) then begin
           plot_legend, names, color=colors, ylog=ylog, xlog=xlog, $
             linestyle=ls, _EXTRA=extra
       endif
       
       return
   endif

   nt = n_elements(time)
   if(nt gt 1) then begin
       nn=strarr(nt)
       if(keyword_set(bw)) then begin
           ls = indgen(nt)
           colors = replicate(color(0,1), nt)
       endif else begin
           if(n_elements(colors) eq 0) then colors = get_colors()
           if(time[0] gt 0) then colors = shift(colors,-1)
           if(n_elements(ls) eq 1) then begin
              ls = replicate(ls,nt)
           endif
       endelse
       if(n_elements(linfac) eq 1) then linfac=replicate(linfac, nt)
       for i=0, n_elements(time)-1 do begin
           newfield = field
           plot_flux_average, newfield, time[i], filename=filename, $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             t=t, rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, fac=fac,$
             linear=linear, multiply_flux=multiply_flux, mks=mks, cgs=cgs, $
             integrate=integrate, complex=complex, asb=aba, phase=phase, $
             stotal=total, rho=rho, nolegend=nolegend, linfac=linfac[i], $
             q_contours=q_contours, regularize=regularize
           lab = parse_units(dimensions(/t0), cgs=cgs, mks=mks)
           t = get_slice_time(filename=filename, slice=time[i], $
                              cgs=cgs, mks=mks)
           nn[i] = string(format='(%"!8t!6 = %g ",A,"!X")', t, lab)
       end

       if(n_elements(names) eq 0) then names=nn
       if(n_elements(names) gt 0 and not keyword_set(nolegend)) then begin
           plot_legend, names, color=colors, ylog=ylog, xlog=xlog, $
             linestyle=ls, _EXTRA=extra
       endif

       return
   endif

   xtitle='!7w!X'
   title = ''

   fa = flux_average(field,slice=time,flux=flux,points=pts,filename=filename, $
                     name=title, symbol=symbol, units=units, bins=bins, $
                     psi=psi,x=x,z=z,t=t,nflux=nflux,linear=linear, fac=fac, $
                     mks=mks, cgs=cgs, area=area, integrate=integrate, $
                     linfac=linfac, sum=sum, _EXTRA=extra, fc=fc, $
                     complex=complex, abs=abs, phase=phase, stotal=total)

   if(n_elements(fa) le 1) then begin
       print, 'Error in flux_average. returning.'
       return
   endif

   ytitle = symbol
   if(strlen(units) gt 0) then ytitle = ytitle + '!6 ('+units+ '!6)!X'

   if(keyword_set(rms)) then begin
       fa2 = flux_average(field^2,flux=flux,nflux=nflux,points=pts,slice=time,$
                          filename=filename, t=t, linear=linear, fac=fac, $
                         mks=mks, cgs=cgs, complex=complex, $
                          abs=abs, phase=phase)
       fa = sqrt(1. - fa^2/fa2)

       ytitle = '!9S!6(1 - !12<' + symbol + '!12>!U2!n/!12<' + $
         symbol + '!6!U2!N!12>!6)!X'
       title = '!6Poloidal Deviation of ' + title + '!X'
   endif

   if(n_elements(multiply_flux) ne 0) then begin
       if(multiply_flux ne 0) then begin
           print, 'multiplying flux by', multiply_flux
           flux = flux*multiply_flux
           nflux=nflux*multiply_flux
       end
   end

   if(keyword_set(srnorm)) then begin
       flux = sqrt(nflux)
       xtitle = '!9r!7W!X'
       lcfs_psi = 1.
   end else if(keyword_set(minor_radius)) then begin
       flux = flux_average('r',points=pts,file=filename,t=t,linear=linear,$
                    name=xtitle,bins=bins,units=units,slice=time, $
                          mks=mks, cgs=cgs, fac=fac)
       xtitle = '!12<' + xtitle + '!12> !6 ('+units+')!X'
   endif else if(keyword_set(rho)) then begin
       flux = fc.rho
       lcfs_psi = 1.
       xtitle = '!7q!X'
   endif else if(norm ne 0) then begin
       flux = nflux
       xtitle = '!7W!X'
       lcfs_psi = 1.
   end

   if(keyword_set(regularize)) then begin
      print, 'REGULARIZING'
      fa = fa/total(fa)
   endif else begin
      print, 'NOT REGULARIZING'
   end

   if(n_elements(sm) eq 1) then begin
       fa = smooth(fa,sm)
   end

   if(keyword_set(overplot)) then begin
       oplot, flux, fa, color=colors, linestyle=ls, _EXTRA=extra
   endif else begin
       if(n_elements(colors) eq 0) then begin
           plot, flux, fa, xtitle=xtitle, linestyle=ls, $
             ytitle=ytitle, title=title, xlog=xlog, ylog=ylog, $
             _EXTRA=extra
       endif else begin
           plot, flux, fa, xtitle=xtitle, linestyle=ls, $
             ytitle=ytitle, title=title, xlog=xlog, ylog=ylog, /nodata, $
             _EXTRA=extra
;              color=color(0), _EXTRA=extra
           oplot, flux, fa, color=colors, linestyle=ls, _EXTRA=extra
       endelse
   endelse

   if(keyword_set(lcfs)) then begin
       oplot, [lcfs_psi,lcfs_psi], !y.crange, linestyle=2, color=colors
   endif

   if(n_elements(q_contours) ne 0) then begin
       fvals = flux_at_q(q_contours, points=pts, filename=filename, fc=fc, $
                         slice=time, normalized_flux=norm, bins=bins, q=q)
       if(fvals[0] gt 0) then begin

       for k=0, n_elements(fvals)-1 do begin
           oplot, [fvals[k], fvals[k]], !y.crange, linestyle=1, color=colors
           xyouts, fvals[k], (!y.crange[1]-!y.crange[0])*0.8 + !y.crange[0], $
                   '!8q!6=' + string(format='(F0.2)', q[k]) +'!X', orient=90
       end
       result = interpol(fa, flux, fvals)
       print, 'filename = ', filename
       val_at_q = result
       qflux = fvals
       print, "q's =  ", q
       print, "flux at q's =  ", qflux
       print, "values at q's =  ", val_at_q

       end
   endif

   if(n_elements(outfile) eq 1) then begin
      if(strlen(outfile) ge 1) then begin
         openw, ifile, outfile, /get_lun
         if(keyword_set(complex)) then begin
            for i=0, n_elements(flux)-1 do begin
               printf, ifile, format='(3E16.6)', $
                       flux[i], real_part(fa[0,i]), imaginary(fa[0,i])
            end
         endif else begin
            for i=0, n_elements(flux)-1 do begin
               printf, ifile, format='(2E16.6)', flux[i], real_part(fa[i])
            end
         endelse
         free_lun, ifile
      end
   endif
end

pro plot_scalar, scalarname, x, filename=filename, names=names, $
                 overplot=overplot, difference=diff, $
                 ylog=ylog, xlog=xlog, absolute_value=absolute, $
                 power_spectrum=pspec, per_length=per_length, $
                 growth_rate=growth_rate, bw=bw, nolegend=nolegend, $
                 cgs=cgs,mks=mks,linestyle=ls, color=co, outfile=outfile, $
                 smooth=sm, compensate_renorm=comp, integrate=integrate, $
                 xscale=xscale, ipellet=ipellet, factor=fac, versus=versus, $
                 xabs=xabs, _EXTRA=extra

  if(n_elements(filename) eq 0) then filename='C1.h5'
  if(n_elements(xscale) eq 0) then xscale=1.
  if(n_elements(fac) eq 0) then fac=1.

  if(n_elements(names) eq 0) then names=filename

  nfiles = n_elements(filename)
  if(nfiles gt 1) then begin
      if(keyword_set(bw)) then begin
          if(n_elements(ls) eq 0) then ls = indgen(nfiles)
          if(n_elements(co) eq 0) then co = replicate(color(0,1),nfiles)
      endif else begin
          if(n_elements(ls) eq 0) then ls = replicate(0, nfiles)
          if(n_elements(co) eq 0) then co = shift(get_colors(),-1)
      endelse
      if(n_elements(x) eq 1) then x = replicate(x, nfiles)

      for i=0, nfiles-1 do begin
          if(n_elements(x) eq 0) then begin
              plot_scalar, scalarname, filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=co[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=pspec, per_length=per_length, $
                growth_rate=growth_rate, linestyle=ls[i], nolegend=nolegend, $
                absolute_value=absolute,cgs=cgs,mks=mks,difference=diff, $
                comp=comp, integrate=integrate, xscale=xscale, ipellet=ipellet, $
                factor=fac, versus=versus, xabs=xabs
          endif else begin
              plot_scalar, scalarname, x[i], filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=co[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=pspec, per_length=per_length, $
                growth_rate=growth_rate, nolegend=nolegend, $
                absolute_value=absolute,cgs=cgs,mks=mks,difference=diff, $
                comp=comp, integrate=integrate, xscale=xscale, ipellet=ipellet, $
                factor=fac, versus=versus, xabs=xabs
          endelse
      end

      if((n_elements(names) gt 0) and (not keyword_set(nolegend))) then begin
          plot_legend, names, ylog=ylog, xlog=xlog, $
            color=co, linestyle=ls, _EXTRA=extra
      endif    

      return
  endif 

  data = read_scalar(scalarname, filename=filename, time=time, ipellet=ipellet, $
                     title=title, symbol=symbol, units=units, cgs=cgs, mks=mks, integrate=integrate)
  data = data*fac
  if(keyword_set(comp)) then data = compensate_renorm(data)
  if(n_elements(data) le 1) then return

  title = '!6' + title + '!X'

  ytitle = symbol

  if(keyword_set(pspec)) then begin
;      xtitle = '!7x!6 (!7s!D!8A!N!6!U-1!N)!X'
      xtitle = make_label('!6Frequency!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)
      data = power_spectrum(data, frequency=tdata, t=max(time))
      print, 'T = ', max(time)
  endif else begin
      xtitle = make_label('!8t!X', /t0, cgs=cgs, mks=mks, _EXTRA=extra)
      tdata = time
  endelse

  if(n_elements(versus) eq 1) then begin
     tdata = read_scalar(versus, filename=filename, time=time, ipellet=ipellet, $
                         title=vtitle, symbol=vsymbol, units=vunits, cgs=cgs, mks=mks)
     if(keyword_set(xabs)) then begin
        tdata = abs(tdata)
        vsymbol = '!3|' + vsymbol + '!3|!X'
     end

     xtitle = vsymbol + ' !6(' + vunits + ')!X'
  end

  if(keyword_set(per_length)) then begin
      itor = read_parameter('itor', filename=filename)
      if(itor eq 1) then begin
          rzero = read_parameter('rzero', filename=filename)
          data = data / rzero
      endif
  endif

  if(keyword_set(diff)) then begin
      data = data - data[0]
      ytitle = '!7D' + ytitle + '!X'
  end

  if(keyword_set(absolute)) then begin
      data = abs(data)
      ytitle = '!3|' + ytitle + '!3|!X'
  end
  
  if(keyword_set(growth_rate)) then begin
      n = min([n_elements(tdata), n_elements(data)])
      if(n_elements(data) lt n) then print, 'truncating data'
      if(n_elements(tdata) lt n) then print, 'truncating tdata'
      data = deriv(tdata(0:n-1), alog(abs(data(0:n-1))))
;      ytitle = '!7c !6(!7s!D!8A!N!6!U-1!N)!X'
      ytitle = make_label('!7c!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)
   endif else begin
      if(strlen(units) gt 0) then begin
         ytitle = ytitle + '!6 (' + units + ')!X'
      endif
   end

  if(keyword_set(sm)) then data = smooth(data, sm)

  if(n_elements(x) eq 0) then begin
      if ipellet eq -1 then begin
        N = size(data)
        N = N[1]
        c = get_colors(N+3)
        if(not keyword_set(overplot)) then begin
          plot, tdata*xscale, data[0,*], xtitle=xtitle, ytitle=ytitle, $
            title=title, _EXTRA=extra, ylog=ylog, xlog=xlog, $
            /nodata
        end
        for n=0,N-1 do oplot, tdata*xscale, data[n,*], color=c[n+1], linestyle=ls, _EXTRA=extra
      endif else begin
        if(not keyword_set(overplot)) then begin
            plot, tdata*xscale, data, xtitle=xtitle, ytitle=ytitle, $
              title=title, _EXTRA=extra, ylog=ylog, xlog=xlog, $
              /nodata
        end
        oplot, tdata*xscale, data, color=co, linestyle=ls, _EXTRA=extra
      endelse
  endif else begin
      xi = x
      x = fltarr(1)
      z = fltarr(1)
      x[0] = xi
      z[0] = data[n_elements(data)-1]

      if(not keyword_set(overplot)) then begin
          plot, x, z, /nodata, $
            title=title, xtitle=xtitle, ytitle=ytitle, $
            _EXTRA=extra, ylog=ylog, xlog=xlog
      end
      oplot, x, z, color=co, linestyle=ls, _EXTRA=extra
  endelse
      if(n_elements(outfile) eq 1) then begin
         openw,ifile,outfile,/get_lun

         d = size(data,/n_dimensions)
         if(d eq 1) then begin
            if(keyword_set(growth_rate)) then begin
               n = min([n_elements(tdata), n_elements(data)])
               printf,ifile,format='(2E16.6)',transpose([[tdata(1:n-1)],[data(1:n-1)]])
            endif else begin
               printf,ifile,format='(2E16.6)',transpose([[tdata],[data]])
            endelse
         endif else begin
            n = min([n_elements(tdata),n_elements(data[0,*])])
            m = n_elements(data[*,0])
            istart = keyword_set(growth_rate)
            form = string(format='("(",I0,"E16.6)")',m+1)
            for i=istart,n-1 do begin
               printf, ifile, format=form, tdata[i], data[*,i]
            end
         end

         free_lun, ifile
      endif
end

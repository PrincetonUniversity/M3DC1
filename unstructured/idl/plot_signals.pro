pro plot_signals, signal, deriv=der, filename=filename, power_spectrum=pspec, $
                     compensate_renorm=comp, overplot=overplot, $
                  scale=scale, tdata=tdata, data=data, noplot=noplot, $
                  outfile=outfile, _EXTRA=extra

  if(n_elements(filename) eq 0) then filename = 'C1.h5'
  if(hdf5_file_test(filename) eq 0) then return

  if(strcmp('mag_probes', signal, /fold_case) eq 1) then begin
     ifl = read_parameter('imag_probes', filename=filename)
     sigdir = 'mag_probes'
     if(keyword_set(der)) then begin
        dim = dimensions(/b0)
        ytitle='!9p!8!Dt!N!5B!9.!5n!6 (' + $
               make_units(/b0,t0=-1,filename=filename,_EXTRA=extra)+')!X'
     endif else begin
        dim = dimensions(/b0)
        ytitle='!5B!9.!5n!6 (' + $
               make_units(/b0,filename=filename,_EXTRA=extra)+')!X'
     endelse
  endif else if(strcmp('flux_loops', signal, /fold_case) eq 1) then begin
     ifl = read_parameter('iflux_loops', filename=filename)
     sigdir = 'flux_loops'
     if(keyword_set(der)) then begin
        dim = dimensions(/b0, l0=2)
        ytitle='!9p!8!Dt!N!7W!D!8p!N!6 (' + $
               make_units(/pot,filename=filename,_EXTRA=extra)+')!X'
     endif else begin
        dim = dimensions(/b0, l0=2)
        ytitle='!7W!D!8p!N!6 (' + $
               make_units(/b0,l0=2,filename=filename,_EXTRA=extra)+')!X'
     endelse
  end

  if(ifl eq 0) then begin
     print, 'No data for signal: ' + signal
     return
  end

  file_id = h5f_open(filename)
  root_id = h5g_open(file_id, "/")
  fl = h5_parse(root_id, sigdir, /read_data)
  data = fl.value._data
  h5g_close, root_id
  h5f_close, file_id

  convert_units, data, dim, filename=filename, _EXTRA=extra
  t = read_scalar('time',_EXTRA=extra,symbol=tsym,units=tu, filename=filename)

  if(keyword_set(comp)) then begin
     for i=0, ifl-1 do begin
        data[i,*] = compensate_renorm(data[i,*])
     end
  end

  if(keyword_set(scale)) then begin
     gamma = read_gamma(filename=filename, _extra=extra)
     for i=0, ifl-1 do begin
        data[i,*] = data[i,*]/exp(t*gamma[0])/data[i,n_elements(t)-1]
     end
  end

  if(keyword_set(der)) then begin
     for i=0, ifl-1 do begin
        data[i,*] = deriv(t,data[i,*])
     end
  endif

  if(keyword_set(pspec)) then begin
     xtitle = make_label('!6Frequency!X', t0=-1, _EXTRA=extra)
     for i=0, ifl-1 do begin
        data[i,*] = power_spectrum(data[i,*], frequency=tdata, t=max(t))
     end
  endif else begin
     xtitle=tsym + ' !6(' + tu + ')!X'
     tdata = t
  end

  if(n_elements(outfile) eq 1) then begin
     openw, ifile, outfile, /get_lun
     for i=0, n_elements(tdata)-1 do begin
        printf, format='(100G12.4)', ifile, tdata[i], data[*,i]
     end
     free_lun, ifile
  end

  if(keyword_set(noplot)) then return

  if(not keyword_set(overplot)) then begin
     plot, tdata, data, /nodata, $
           xtitle=xtitle, ytitle=ytitle, _EXTRA=extra
  end

  c = shift(get_colors(ifl),-1)
  for i=0, ifl-1 do begin
     oplot, tdata, data[i,*], color=c[i]
  end
end

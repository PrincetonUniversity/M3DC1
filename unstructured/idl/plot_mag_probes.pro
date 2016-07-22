pro plot_mag_probes, deriv=der, filename=filename, power_spectrum=pspec, $
                     compensate_renorm=comp, _EXTRA=extra

  if(n_elements(filename) eq 0) then filename = 'C1.h5'
  if(hdf5_file_test(filename) eq 0) then return

  ifl = read_parameter('imag_probes', filename=filename)
  if(ifl eq 0) then begin
     print, 'No mag probe data'
     return
  end

  file_id = h5f_open(filename)
  root_id = h5g_open(file_id, "/")
  fl = h5_parse(root_id, "mag_probes", /read_data)
  data = fl.value._data
  h5g_close, root_id
  h5f_close, file_id

  convert_units, data, dimensions(/b0), filename=filename, _EXTRA=extra
  t = read_scalar('time',_EXTRA=extra,symbol=tsym,units=tu, filename=filename)

  if(keyword_set(comp)) then begin
     for i=0, ifl-1 do begin
        data[i,*] = compensate_renorm(data[i,*])
     end
  end

  if(keyword_set(der)) then begin
     for i=0, ifl-1 do begin
        data[i,*] = deriv(t,data[i,*])
     end
     ytitle='!9p!8!Dt!N!5B!9.!5n!6 (' + $
            make_units(/b0,t0=-1,filename=filename,_EXTRA=extra) + ')!X'
  endif else begin
     ytitle='!5B!9.!5n!6 (' + $
            make_units(/b0,filename=filename,_EXTRA=extra) + ')!X'
  endelse

  if(keyword_set(pspec)) then begin
     xtitle = make_label('!6Frequency!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)
     for i=0, ifl-1 do begin
        data[i,*] = power_spectrum(data[i,*], frequency=tdata, t=max(t))
     end
     print, 'T = ', max(t)
  endif else begin
     xtitle=tsym + ' !6(' + tu + ')!X'
     tdata = t
  end

  plot, tdata, data, /nodata, $
        xtitle=xtitle, ytitle=ytitle, _EXTRA=extra

  c = shift(get_colors(),-1)
  for i=0, ifl-1 do begin
     oplot, tdata, data[i,*], color=c[i]
  end
end

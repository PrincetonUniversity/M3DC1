pro plot_flux_loops, deriv=der, filename=filename, _EXTRA=extra

  if(n_elements(filename) eq 0) then filename = 'C1.h5'
  if(hdf5_file_test(filename) eq 0) then return

  ifl = read_parameter('iflux_loops', filename=filename)
  if(ifl eq 0) then begin
     print, 'No flux loop data'
     return
  end

  file_id = h5f_open(filename)
  root_id = h5g_open(file_id, "/")
  fl = h5_parse(root_id, "flux_loops", /read_data)
  data = fl.value._data

  convert_units, data, dimensions(/b0), _EXTRA=extra
  t = read_scalar('time',_EXTRA=extra,symbol=tsym,units=tu)

  if(keyword_set(der)) then begin
     for i=0, ifl-1 do begin
        data[i,*] = deriv(t,data[i,*])
     end
     ytitle='!9p!8!Dt!N!5B!9.!5n!6 (' + make_units(/b0,t0=-1, _EXTRA=extra) + ')!X'
  endif else begin
     ytitle='!5B!9.!5n!6 (' + make_units(/b0, _EXTRA=extra) + ')!X'
  endelse

  print, data
  
  plot, t, data, /nodata, $
        xtitle=tsym + ' !6(' + tu + ')!X', $
        ytitle='!5B!9.!5n!6 (' + make_units(/b0, _EXTRA=extra) + ')!X'

  c = shift(get_colors(),-1)
  for i=0, ifl-1 do begin
     oplot, t, data[i,*], color=c[i]
  end
end

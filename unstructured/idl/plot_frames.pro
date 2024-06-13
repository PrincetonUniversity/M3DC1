pro plot_frames, field, filename=filename, slices=s, time=t, dt=dt, $
                 outfile=outfile, min_slice=mins, min_frame=minf, _EXTRA=extra

  if(n_elements(filename) eq 0) then filename='C1.h5'
  if(n_elements(mins) eq 0) then mins=0
  if(n_elements(minf) eq 0) then minf=0
  
  if(n_elements(dt) eq 1) then begin
     ntimes = read_parameter('ntime',filename=filename)
     s0 = indgen(ntimes)
     t0 = get_slice_time(filename=filename, slice=s0)
     tmax = t0[ntimes-1]
     n = tmax/dt
     s = intarr(n)
     t = intarr(n)
     t_req = dt*indgen(n)
     for i=0, n-1 do begin
        diff = min(t0-t_req[i], j, /abs)
        s[i] = j
        t[i] = t0[j]
        print, i, t_req[i], t[i]
     end
  end

  if(n_elements(outfile) eq 0) then begin
     plot, t_req, 0*t_req, yran=[-1,1], psym=4
     oplot, t, 0.*t+0.1, psym=5
     return
  end
  
  n = n_elements(s)
  for i=minf, n-1 do begin
     if(s[i] lt mins) then continue
     setplot, 'ps'
     device, /color, /encapsulated, $
             file=outfile+'_'+string(FORMAT='(I4.4)', i)+'.eps'
     plot_field, field, s[i], filename=filename, _EXTRA=extra
     device, /close
  end

end

function cumtrapz,x,y

  npts = n_elements(x)

  cum = fltarr(npts)

  for i=1,npts-1 do begin

    cum[i] = cum[i-1] + 0.5*(x[i]-x[i-1])*(y[i]+y[i-1])

  endfor

  return, cum

end
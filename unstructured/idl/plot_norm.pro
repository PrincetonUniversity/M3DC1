pro plot_norm, filename, _EXTRA=extra, len=len
  if(n_elements(len) eq 0) then len=0.1
  if(n_elements(filename) eq 0) then filename='normcurv'
  a = read_ascii(filename)

  x = a.field1[0,*]
  z = a.field1[1,*]
  nx = a.field1[2,*]
  nz = a.field1[3,*]

  plot, x, z, psym=4, _EXTRA=extra

  n = n_elements(x)

  for i=0, n-1 do begin
    oplot, [x[i],x[i]+nx[i]*len],[z[i],z[i]+nz[i]*len]
  end
  
end

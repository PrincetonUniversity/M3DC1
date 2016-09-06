pro extract_boundary, shot, time, points=n
  x = readg(shot, time)

  if(n_elements(n) eq 0) then n=100

  rl = reform(x.lim[0,*])
  zl = reform(x.lim[1,*])

  rmag = mean(rl)
  zmag = 0.

  t = atan(zl-zmag, rl-rmag)

  j = sort(t[uniq(t)])
  r0 = rl[j]
  z0 = zl[j]
  theta0 = t[j]

  theta = 2.*!pi*findgen(n)/n - !pi

  m = n_elements(theta0)-1
  if(theta0[0] gt theta[0]) then begin
     f0 = !pi-theta0[m]
     f1 = !pi+theta0[0]
     print, 'appending to start. f0, f1 = ', f0, f1
     r0 = [(r0[0]*f0+r0[m]*f1) / (f0+f1), r0]
     z0 = [(z0[0]*f0+z0[m]*f1) / (f0+f1), z0]
     theta0 = [-!pi, theta0]
  end
  m = n_elements(theta0)-1
  if(theta0[m] lt theta[n-1]) then begin
     f0 = !pi-theta0[m]
     f1 = !pi+theta0[0]
     print, 'appending to end. f0, f1 = ', f0, f1
     r0 = [r0, (r0[0]*f0+r0[m]*f1)/(f0+f1)]
     z0 = [z0, (z0[0]*f0+z0[m]*f1)/(f0+f1)]
     theta0 = [theta0, !pi]
  end
  
  r = interpol(r0, theta0, theta)
  z = interpol(z0, theta0, theta)

  openw, ifile, 'boundary.dat', /get_lun
  printf, ifile, n
  for i=0, n-1 do printf, ifile, r[i], z[i]
  free_lun, ifile

  ct3
  plot, r0, z0, /iso
  oplot, r0, z0, psym=2
  oplot, r, z, color=color(1), linestyle=2
end

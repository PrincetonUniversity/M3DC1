function fit_miller, r, z, plot=plot
  r0 = (max(r) + min(r))/2.
  z0 = (max(z) + min(z))/2.
  a = (max(r) - min(r))/2.
  b = (max(z) - min(z))/2.
  kappa = b/a
  print, a, b, max(z), min(z)

  z_top = max(z, i_top)
  z_low = min(z, i_low)
  
  delta_top = (r0 - r[i_top])/a
  delta_low = (r0 - r[i_low])/a
  delta = (delta_top + delta_low)/2.

  n = 100
  theta = findgen(n)/(n-1.) * 2.*!pi
  r_miller = r0 + a*cos(theta + asin(delta) * sin(theta))
  z_miller = z0 + b*sin(theta)
    
  if(keyword_set(plot)) then begin
     plot, r, z
     oplot, r_miller, z_miller, color=color(1)
  end
  
  geom = { R0:r0, Z0:z0, a:a, kappa:kappa, delta:delta, $
         delta_upper:delta_top, delta_lower:delta_low }
  return, geom

end

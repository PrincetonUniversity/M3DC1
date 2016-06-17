pro find_aug_currents, _EXTRA=extra
  jy = read_field('jy',x,y,t,_EXTRA=extra, /mks)
  r = radius_matrix(x,y,t)
  z = z_matrix(x,y,t)

  umask = (z gt 0.5) and (z lt 0.9) and (z gt -(r-2.1)*2. + 0.5) and (z lt -(r-2.25)*2.+0.6)
  lmask = (z gt -1) and (z lt -0.5) and (z gt (r-2.0)*2. - 0.9) and (z lt (r-1.75)*2. - 0.9)

  jyu = jy*umask
  jyl = jy*lmask
  contour_and_legend, jy - (jyu + jyl), x, z

  dx = x[1]-x[0]
  dy = y[1]-y[0]
  iu = total(jyu*dx*dy)
  il = total(jyl*dx*dy)

  print, 'Upper coil current = ', iu/1e3, ' kA-turns'
  print, 'Lower coil current = ', il/1e3, ' kA-turns'
end

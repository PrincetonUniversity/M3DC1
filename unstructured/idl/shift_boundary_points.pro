pro shift_boundary_points, infile, d, _EXTRA=extra
  in = read_ascii(infile, data=1)

  x = reform(in.field1[0,*])
  z = reform(in.field1[1,*])
  n = n_elements(x)
  dx = deriv([x,x[0]])
  dz = deriv([z,z[0]])
  dl = sqrt(dx^2 + dz^2)

  ct3
  plot, x, z, _EXTRA=extra

  x = x + d*dz/dl
  z = z - d*dx/dl

  oplot, x, z, color=color(1)

  outfile = infile + '.shifted'
  openw, ifile, outfile, /get_lun
  printf, ifile, n
  for i=0, n-1 do printf, ifile, x[i], z[i]
  free_lun, ifile

end

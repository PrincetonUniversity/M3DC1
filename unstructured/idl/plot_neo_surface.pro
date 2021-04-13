pro plot_neo_surface, file, psin_range=psin_range, scale=scale, overplot=op, $
                      phi=phi

  if(n_elements(psin_range) eq 0) then psin_range = [0.,1.]
  if(n_elements(psin_range) eq 1) then psin_range = [psin_range, psin_range]
  if(n_elements(scale) eq 0) then scale = 1.
  if(n_elements(phi) eq 0) then phi = 0.

  a = read_neo_input(file)

  psin = (a.psi - a.psi_0) / (a.psi_1 - a.psi_0)

  surf = where((psin ge psin_range[0]) and (psin le psin_range[1]), count)
  
  if(count eq 0) then begin
     print, 'Error: no surfaces found in psin range ', psin_range
     return
  end

  dp = min(a.Phi - phi, /abs, m)
  print, 'Plotting phi = ', a.Phi[m]
  
  if(not keyword_set(op)) then begin
     plot, [min(a.R[surf,*,*]), max(a.R[surf,*,*])], $
           [min(a.Z[surf,*,*]), max(a.Z[surf,*,*])], /nodata
  end

  c = get_colors(count)

  r_av = fltarr(a.np)
  z_av = fltarr(a.np)

  foreach s, surf, i do begin
     for k=0, a.np-1 do begin
        r_av[k] = mean(a.R[s,*,k])
        z_av[k] = mean(a.Z[s,*,k])
     end

     R = (a.R[s,m,*] - r_av)*scale + r_av
     Z = (a.Z[s,m,*] - z_av)*scale + z_av

     oplot, R, Z, color=c[i]
  end

end

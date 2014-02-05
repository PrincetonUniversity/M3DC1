pro plot_omega, filename=filename, slice=time, points=pts, $
                yrange=yrange, q_val=q_val, $
                mtop=mtop, mslope=mslope, _EXTRA=extra, $
                plot_wstar=plot_wstar, plot_wi=plot_wi, plot_we=plot_we

  if(n_elements(pts) eq 0) then pts=200
  if(n_elements(plot_wi) eq 0) then plot_wi=1
  if(n_elements(plot_wstar) eq 0) then plot_wstar=1
  if(n_elements(plot_we) eq 0) then plot_we=1

  db = read_parameter('db', filename=filename)
  itor = read_parameter('itor', filename=filename)
  print, 'db = ', db

  omega = read_field('omega', x, y,t,slices=time,filename=filename,points=pts)
  p = read_field('p', x, y, t, slices=time, filename=filename, points=pts)
  pe = read_field('pe', x, y, t, slices=time, filename=filename, points=pts)
  u = read_field('phi', x, y, t, slices=time, filename=filename, points=pts)
  chi = read_field('chi', x, y, t, slices=time, filename=filename, points=pts)
  psi = read_field('psi', x, y, t, slices=time, filename=filename, points=pts)
  i = read_field('I', x, y, t, slices=time, filename=filename, points=pts)
  den = read_field('den', x, y, t, slices=time, filename=filename, points=pts)

  if(itor eq 1) then begin
      r = radius_matrix(x,y,t)
  endif else r = 1.
       
  v_omega = omega - i/(r^2*s_bracket(psi,psi,x,y)) * $
    (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r)
  w_star_i = db*s_bracket(p-pe,psi,x,y)/s_bracket(psi,psi,x,y) / den
  w_star_e = -db*s_bracket(pe,psi,x,y)/s_bracket(psi,psi,x,y) / den

  get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, $
    zeff=zeff, ion_mass=ion_mass, filename=filename
  convert_units, v_omega, dimensions(t0=-1), $
    b0_norm, n0_norm, l0_norm, zeff, ion_mass, /mks
  convert_units, w_star_i, dimensions(t0=-1), $
    b0_norm, n0_norm, l0_norm, zeff, ion_mass, /mks
  convert_units, w_star_e, dimensions(t0=-1), $
    b0_norm, n0_norm, l0_norm, zeff, ion_mass, /mks

  omega_ExB = v_omega - w_star_i
  ve_omega = omega_ExB + w_star_e

  v_omega_fa = flux_average_field(v_omega, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, _EXTRA=extra)
  w_star_i_fa = flux_average_field(w_star_i, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, _EXTRA=extra)
  omega_ExB_fa = flux_average_field(omega_ExB, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, _EXTRA=extra)
  ve_omega_fa = flux_average_field(ve_omega, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, _EXTRA=extra)

  xtitle = '!7W!X'
  ytitle = '!6krad/s!X'

  omega_ExB_fa= omega_ExB_fa / 1000.
  v_omega_fa = v_omega_fa / 1000.
  ve_omega_fa = ve_omega_fa / 1000.
  w_star_i_fa = w_star_i_fa / 1000.

  if(n_elements(yrange) eq 0) then begin
      yrange = [min([v_omega_fa, w_star_i_fa, omega_ExB_fa, ve_omega_fa]), $
                max([v_omega_fa, w_star_i_fa, omega_ExB_fa, ve_omega_fa])]   
  end

  ct3
  plot, nflux, omega_ExB_fa, yrange=yrange, xtitle=xtitle, ytitle=ytitle, $
    _EXTRA=extra
  names = '!7x!6!DE!9X!6B!N!X'
  col = color(0)
  if(keyword_set(plot_wi)) then begin
      oplot, nflux, v_omega_fa, color=color(1)
      names = [names, '!7x!X']
      col = [col, color(1)]
  end
  if(keyword_set(plot_we)) then begin
      oplot, nflux, ve_omega_fa, color=color(2)
      names = [names, '!7x!D!8e!N!X']
      col = [col, color(2)]
  end
  if(keyword_set(plot_wstar)) then begin
      oplot, nflux, w_star_i_fa, color=color(3)
      names = [names, '!7x!6!D*!8i!N!X']
      col = [col, color(3)]
  end
  oplot, !x.crange, [0,0], linestyle=2

  plot_legend, names, color=col, _EXTRA=extra

  if(n_elements(q_val) ne 0) then begin
      ntor = read_parameter('ntor', filename=filename)
      m = fix(q_val*ntor+0.001)
      psin = flux_at_q(q_val,slice=slice,filename=filename,$
                       points=pts,bins=bins,/norm,_EXTRA=extra)
      print, psin
      if(n_elements(mtop) eq 0) then mtop = 0.05
      if(n_elements(mslope) eq 0) then mslope = 0.
      top =  $
        (!y.crange[1] - !y.crange[0])*mslope*findgen(n_elements(q_val)) $
        +(!y.crange[1] - !y.crange[0])*(1.-mtop) $
        + !y.crange[0]
      m_str = string(format='(I2)', m)
      m_str[0] = 'm = ' + m_str[0]
      for i=0, n_elements(psin)-1 do begin
          oplot, [psin[i],psin[i]], !y.crange, linestyle=1
          xyouts, psin[i], top[i], m_str[i], charsize=!p.charsize/2.
      end
  end
end

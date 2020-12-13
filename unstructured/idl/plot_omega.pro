pro plot_omega, filename=filename, slice=time, points=pts, $
                yrange=yrange, q_val=q_val, out=out, $
                mtop=mtop, mslope=mslope, _EXTRA=extra

  if(n_elements(pts) eq 0) then pts=200
  if(n_elements(time) eq 0) then time=-1

  db = read_parameter('db', filename=filename)
  itor = read_parameter('itor', filename=filename)
  print, 'db = ', db

  omega = read_field('omega', x, y,t,slices=time,filename=filename,points=pts)
  p_r = read_field('p', x, y, t, slices=time, filename=filename, points=pts, op=2)
  p_z = read_field('p', x, y, t, slices=time, filename=filename, points=pts, op=3)
  pe_r = read_field('pe', x, y, t, slices=time, filename=filename, points=pts, op=2)
  pe_z = read_field('pe', x, y, t, slices=time, filename=filename, points=pts, op=3)
  u_r = read_field('phi', x, y, t, slices=time, filename=filename, points=pts, op=2)
  u_z = read_field('phi', x, y, t, slices=time, filename=filename, points=pts, op=3)
  chi_r = read_field('chi', x, y, t, slices=time, filename=filename, points=pts, op=2)
  chi_z = read_field('chi', x, y, t, slices=time, filename=filename, points=pts, op=3)
  psi = read_field('psi', x, y, t, slices=time, filename=filename, points=pts, /equilibrium)
  psi_r = read_field('psi', x, y, t, slices=time, filename=filename, points=pts, /equilibrium, op=2)
  psi_z = read_field('psi', x, y, t, slices=time, filename=filename, points=pts, /equilibrium, op=3)
  i = read_field('I', x, y, t, slices=time, filename=filename, points=pts, /equilibrium)
  den = read_field('den', x, y, t, slices=time, filename=filename, points=pts, /equilibrium)

  fc = flux_coordinates(_EXTRA=extra, points=pts, psi0=psi, i0=i, x=x, z=y, filename=filename, slice=time)

  if(itor eq 1) then begin
      r = radius_matrix(x,y,t)
  endif else r = 1.
  
  psipsi = psi_r^2 + psi_z^2
  pprime = (p_r*psi_r + p_z*psi_z)/psipsi
  peprime = (pe_r*psi_r + pe_z*psi_z)/psipsi
  piprime = pprime - peprime

  v_omega = omega  $
            - i/(r^2*psipsi) * $
            (r^2*(u_r*psi_r + u_z*psi_z) + $
             (chi_z*psi_r - chi_r*psi_z)/r)
  w_star_i = db*piprime / den
  w_star_e = -db*peprime / den

  get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, $
    ion_mass=ion_mass, filename=filename
  convert_units, v_omega, dimensions(t0=-1), $
    b0_norm, n0_norm, l0_norm, ion_mass, /mks
  convert_units, w_star_i, dimensions(t0=-1), $
    b0_norm, n0_norm, l0_norm, ion_mass, /mks
  convert_units, w_star_e, dimensions(t0=-1), $
    b0_norm, n0_norm, l0_norm, ion_mass, /mks

  omega_ExB = v_omega - w_star_i
  ve_omega = omega_ExB + w_star_e

  v_omega_fa = flux_average_field(v_omega, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, fc=fc, _EXTRA=extra)
  w_star_i_fa = flux_average_field(w_star_i, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, fc=fc, _EXTRA=extra)
  w_star_e_fa = flux_average_field(w_star_e, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, fc=fc, _EXTRA=extra)
  omega_ExB_fa = flux_average_field(omega_ExB, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, fc=fc, _EXTRA=extra)
  ve_omega_fa = flux_average_field(ve_omega, psi, x, y, t, file=filename, $
                           nflux=nflux, bins=pts, fc=fc, _EXTRA=extra)

  xtitle = '!7W!X'
  ytitle = '!6krad/s!X'

  omega_ExB_fa= omega_ExB_fa / 1000.
  v_omega_fa = v_omega_fa / 1000.
  ve_omega_fa = ve_omega_fa / 1000.
  w_star_i_fa = w_star_i_fa / 1000.
  w_star_e_fa = w_star_e_fa / 1000.

  if(n_elements(yrange) eq 0) then begin
      yrange = [min([v_omega_fa, w_star_i_fa, omega_ExB_fa, ve_omega_fa]), $
                max([v_omega_fa, w_star_i_fa, omega_ExB_fa, ve_omega_fa])]
  end

  ct3
  plot, nflux, omega_ExB_fa, yrange=yrange, xtitle=xtitle, ytitle=ytitle, $
    _EXTRA=extra
  names = '!7x!6!DE!9X!6B!N!X'
  col = color(0)
  ls = 0

  ; omega_i
  oplot, nflux, v_omega_fa, color=color(1), _EXTRA=extra
  names = [names, '!7x!X']
  col = [col, color(1)]
  ls = [ls, 0]

  ; omega_e
  oplot, nflux, ve_omega_fa, color=color(2), _EXTRA=extra
  names = [names, '!7x!D!8e!N!X']
  col = [col, color(2)]
  ls = [ls, 0]

  ; omega_*i
  oplot, nflux, w_star_i_fa, color=color(3), _EXTRA=extra, linestyle=1
  names = [names, '!7x!6!D*!8i!N!X']
  col = [col, color(3)]
  ls = [ls, 1]

  ; omega_*e
  oplot, nflux, w_star_e_fa, color=color(4), _EXTRA=extra, linestyle=1
  names = [names, '!7x!6!D*!8e!N!X']
  col = [col, color(4)]
  ls = [ls, 1]

  oplot, !x.crange, [0,0], linestyle=2

  plot_legend, names, color=col, linestyle=ls, _EXTRA=extra

  if(n_elements(out) ne 0) then begin
     openw, ifile, out, /get_lun
     printf, ifile, format='(6A14)', $
             'psi_N', 'omega_ExB', 'omega_i', 'omega_e', 'omega_*i', 'omega_*e'
     for i=0, n_elements(nflux)-1 do begin
        printf, ifile, format='(6E14.5)', $
               nflux[i], omega_ExB_fa[i], $
                v_omega_fa[i], ve_omega_fa[i], w_star_i_fa[i], w_star_e_fa[i]
     end
     free_lun, ifile
  end

  if(n_elements(q_val) ne 0) then begin
      ntor = read_parameter('ntor', filename=filename)
      psin = flux_at_q(q_val,slice=slice,filename=filename,q=q,fc=fc,$
                       points=pts,bins=bins,/norm,_EXTRA=extra)
      m = fix(q*ntor+0.001)
      if(n_elements(mtop) eq 0) then mtop = 0.05
      if(n_elements(mslope) eq 0) then mslope = 0.
      top =  $
        (!y.crange[1] - !y.crange[0])*mslope*findgen(n_elements(q)) $
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

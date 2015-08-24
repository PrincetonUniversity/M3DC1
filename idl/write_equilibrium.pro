pro write_equilibrium, _EXTRA=extra, bmncdf=bmncdf

  if(n_elements(bmncdf) eq 0) then bmncdf='eq.nc'

  fc = flux_coordinates(_EXTRA=extra, /plot)

  F = flux_average('I',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                   bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
  p = flux_average('p',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                   bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
  pe = flux_average('pe',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                    bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
  den_e = flux_average('ne',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                       bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)

  n = fc.n  ; number of flux surfaces
  m = fc.m  ; number of poloidal points

  id = ncdf_create(bmncdf, /clobber)
  n_id = ncdf_dimdef(id, 'npsi', n)
  m_id = ncdf_dimdef(id, 'mpol', m)
  psi_var = ncdf_vardef(id, 'psi', [n_id], /float)
  theta_var = ncdf_vardef(id, 'theta', [n_id], /float)
  psin_var = ncdf_vardef(id, 'psi_norm', [n_id], /float)
  q_var = ncdf_vardef(id, 'q', [n_id], /float)
  p_var = ncdf_vardef(id, 'p', [n_id], /float)
  F_var = ncdf_vardef(id, 'I', [n_id], /float)
  pe_var = ncdf_vardef(id, 'pe', [n_id], /float)
  ne_var = ncdf_vardef(id, 'ne', [n_id], /float)
  r_var = ncdf_vardef(id, 'R', [m_id,n_id], /float)
  z_var = ncdf_vardef(id, 'Z', [m_id,n_id], /float)
  bp_var = ncdf_vardef(id, 'Jacobian', [m_id,n_id], /float)
  ncdf_control, id, /endef
  ncdf_varput, id, 'psi', fc.psi
  ncdf_varput, id, 'psi_norm', fc.psi_norm
  ncdf_varput, id, 'theta', fc.theta
  ncdf_varput, id, 'q', fc.q
  ncdf_varput, id, 'p', reform(p)
  ncdf_varput, id, 'I', reform(F)
  ncdf_varput, id, 'pe', reform(pe)
  ncdf_varput, id, 'ne', reform(den_e)
  ncdf_varput, id, 'R', fc.r
  ncdf_varput, id, 'Z', fc.z
  ncdf_varput, id, 'Jacobian', fc.j
  ncdf_close, id

end

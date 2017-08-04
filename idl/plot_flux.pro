pro plot_flux, slice=slice, _EXTRA=extra, out=out

  db = read_parameter('db', _EXTRA=extra)
  ip = read_scalar('ip', _EXTRA=extra)
  print, 'db = ', db
  print, 'ip = ', ip

  EX1 = read_field('E_R',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  EY1 = read_field('E_Phi',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  EZ1 = read_field('E_Z',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  VX1 = read_field('Vx',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  VY1 = read_field('Vy',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  VZ1 = read_field('Vz',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  BX1 = read_field('Bx',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  BY1 = read_field('By',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  BZ1 = read_field('Bz',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  JX1 = read_field('Jx',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  JY1 = read_field('Jy',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  JZ1 = read_field('Jz',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  p1  = read_field('p', x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  pe1 = read_field('pe',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)
  pi1 = p1 - pe1
  ne1 = read_field('ne',x,z,t,slice=slice,/linear,/complex,_EXTRA=extra)

  BX0 = read_field('Bx',x,z,t,slice=-1,_EXTRA=extra)
  BY0 = read_field('By',x,z,t,slice=-1,_EXTRA=extra)
  BZ0 = read_field('Bz',x,z,t,slice=-1,_EXTRA=extra)
  VY0 = read_field('Vy',x,z,t,slice=-1,_EXTRA=extra)
  JY0 = read_field('Jy',x,z,t,slice=-1,_EXTRA=extra)
  p0  = read_field('p ',x,z,t,slice=-1,_EXTRA=extra)
  pe0 = read_field('pe',x,z,t,slice=-1,_EXTRA=extra)
  pi0 = p0 - pe0
  ne0 = read_field('ne',x,z,t,slice=-1,_EXTRA=extra)
  
  ; Radial direction is B0 x Y.  This is outward when IP > 0.
  ; B0 x Y = Bx Z - Bz X

  B0 = sqrt(BX0^2 + BY0^2 + BZ0^2)

  Br1 = (Bz1*Bx0 - Bx1*Bz0) / B0

  ; V_ExB = (E1 x B0) / B0^2  . (B0 x Y) / B0
  V_ExB = ((Ex1*By0 - Ey1*Bx0) * Bx0 - $
           (Ey1*Bz0 - Ez1*By0) * Bz0) / B0^3
  
  Vi_par0 = Vy0*By0 / B0
  Vi_par1 = (Vx1*Bx0 + Vy1*By0 + Vz1*Bz0) / B0

  Ve_par0 = Vi_par0 - (db/ne0)*Jy0*By0 / B0
  Ve_par1 = Vi_par1 - (db/ne0)*(Jx1*Bx0 + Jy1*By0 + Jz1*Bz0) / B0

  gamma_ne_ExB = real_part(conj(ne1)*V_ExB)/2.
  gamma_pe_ExB = real_part(conj(pe1)*V_ExB)/2.
  gamma_pi_ExB = real_part(conj(pi1)*V_ExB)/2.
  gamma_ne_flutter = real_part((ne1*Ve_par0 + ne0*Ve_par1)*conj(Br1)/B0)/2.
  gamma_pe_flutter = real_part((pe1*Ve_par0 + pe0*Ve_par1)*conj(Br1)/B0)/2.
  gamma_pi_flutter = real_part((pi1*Vi_par0 + pi0*Vi_par1)*conj(Br1)/B0)/2.

  if(ip[0] lt 0.) then begin
     gamma_ne_ExB = -gamma_ne_ExB
     gamma_pe_ExB = -gamma_pe_ExB
     gamma_pi_ExB = -gamma_pi_ExB
     gamma_ne_flutter = -gamma_ne_flutter
     gamma_pe_flutter = -gamma_pe_flutter
     gamma_pi_flutter = -gamma_pi_flutter
  end

  g_ne_ExB = flux_average(gamma_ne_ExB,fc=fc,psi=psi, i0=i0, x=x, z=z, t=t,$
                          _EXTRA=extra)
  g_pe_ExB = flux_average(gamma_pe_ExB,fc=fc,psi=psi, i0=i0, x=x, z=z, t=t,$
                          _EXTRA=extra)
  g_pi_ExB = flux_average(gamma_pi_ExB,fc=fc,psi=psi, i0=i0, x=x, z=z, t=t, $
                          _EXTRA=extra)
  g_ne_flutter = flux_average(gamma_ne_flutter,fc=fc,psi=psi, i0=i0, x=x, z=z,$
                              t=t, _EXTRA=extra)
  g_pe_flutter = flux_average(gamma_pe_flutter,fc=fc,psi=psi, i0=i0, x=x, z=z,$
                              t=t, _EXTRA=extra)
  g_pi_flutter = flux_average(gamma_pi_flutter,fc=fc,psi=psi, i0=i0, x=x, z=z,$
                              t=t, _EXTRA=extra)

  dn = dimensions(/n0, /v0)
  un = parse_units(dn, /mks)
  dp = dimensions(/p0, /v0)
  up = parse_units(dp, /mks)
  convert_units, g_ne_ExB, dn, /mks, _EXTRA=extra
  convert_units, g_ne_flutter, dn, /mks, _EXTRA=extra
  convert_units, g_pe_ExB, dp, /mks, _EXTRA=extra
  convert_units, g_pe_flutter, dp, /mks, _EXTRA=extra
  convert_units, g_pi_ExB, dp, /mks, _EXTRA=extra
  convert_units, g_pi_flutter, dp, /mks, _EXTRA=extra
   
  !p.multi = [0,2,1]

  ct3
  print, [min([g_ne_ExB, g_ne_flutter]), max([g_ne_ExB, g_ne_flutter])]
  plot, [0,1], [min([g_ne_ExB, g_ne_flutter]), max([g_ne_ExB, g_ne_flutter])],$
        /nodata, title='!6Electron Particle Flux Density!X', $
        xtitle='!7W!X', ytitle='!6Outward Particle Flux Density(' + un + ')!X'
  oplot, fc.psi_norm, g_ne_ExB + g_ne_flutter, color=color(0)
  oplot, fc.psi_norm, g_ne_ExB, color=color(1)
  oplot, fc.psi_norm, g_ne_flutter, color=color(2)

  plot_legend, ['Total', 'ExB flux', 'Flutter'], color=get_colors()
  
  plot, [0,1], [min([g_pe_ExB, g_pe_flutter]), max([g_pe_ExB, g_pe_flutter])],$
        /nodata, title='!6Electron Heat Flux Density!X', $
        xtitle='!7W!X', ytitle='!6Outward Energy Flux Density (' + up + ')!X'
  oplot, fc.psi_norm, g_pe_ExB + g_pe_flutter, color=color(0)
  oplot, fc.psi_norm, g_pe_ExB, color=color(1)
  oplot, fc.psi_norm, g_pe_flutter, color=color(2)

  plot_legend, ['Total', 'ExB flux', 'Flutter'], color=get_colors()

  !p.multi = 0

  if(n_elements(out) eq 1) then begin
     openw, ifile, out, /get_lun
     printf, ifile, format='(7A12)', $
             'psi norm', 'ne ExB', 'ne flutter', $
             'pe ExB', 'pe flutter', $
             'pi ExB', 'pi flutter'
     printf, ifile, format='(7A12)', $
             ' ', '[/s m^2]', '[/s m^2]', $
             '[J/s m^2]', '[J/s m^2]', $
             '[J/s m^2]', '[J/s m^2]'
     for i=0, n_elements(fc.psi_norm)-1 do begin
        printf, ifile, format='(7G12.5)', $
                fc.psi_norm[i], g_ne_ExB[i], g_ne_flutter[i], $
                g_pe_ExB[i], g_pe_flutter[i], g_pi_ExB[i], g_pi_flutter[i]
     end
     free_lun, ifile
  end
end

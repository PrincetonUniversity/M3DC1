pro schaffer_plot, field, x,z,t, q=q, bins=bins, q_val=q_val, $
                   psi_val=psi_val, ntor=ntor, psi0=psi0, i0=i0, $
                   m_val=m_val, phase=phase, overplot=overplot, $
                   linestyle=linestyle, outfile=outfile, bmnfile=bmnfile, $
                   bmncdf=bmncdf, rhs=rhs, reverse_q=reverse_q, $
                   sqrtpsin=sqrtpsin, profdata=profdata, $
                   boozer=boozer, pest=pest, hamada=hamada, geo=geo, $
                   symbol=symbol, units=units, $
                   dpsi0_dx=psi0_r, dpsi0_dz=psi0_z, $
                   _EXTRA=extra

   print, 'Drawing schaffer plot'

   if(not keyword_set(boozer) and not keyword_set(hamada) and $
      not keyword_set(geo)) then pest=1

   if(n_elements(sqrtpsin) eq 0) then sqrtpsin=1

   if(n_elements(label) eq 1 and n_elements(units) eq 1) then begin
      label = symbol + '!6 (' + units + '!6)!X'
   endif else begin
      label = ' '
   end

   if(n_elements(psi0) eq 0) then begin
;       psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
      print, 'READING PSI IN SCHAFFER_PLOT'
       psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
   endif
   if(n_elements(i0) eq 0) then begin
;       i0   = read_field('i'  ,x,z,t,/equilibrium,_EXTRA=extra)
       i0   = read_field('i'  ,x,z,t,slice=-1,_EXTRA=extra)
   endif
   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor',_EXTRA=extra)
   endif

   itor = read_parameter('itor', _EXTRA=extra)
   if(itor eq 1) then begin
      r = radius_matrix(x,z,t)
   endif else begin
      r = 1.
   endelse

   ; From Schaffer 2008
   ; Br_mn = [(2*pi)^2/S] * [J Br]_mn
   bp = sqrt(s_bracket(psi0,psi0,x,z))/r

   if(size(field, /type) eq 7) then begin
       field = read_field(field,x,z,t,/complex,_EXTRA=extra)
   endif

   d = field_spectrum(field,x,z,psi0=psi0,i0=i0,fc=fc,tbins=bins,fbins=bins, $
                      m=m,n=n,pest=pest,boozer=boozer,hamada=hamada, $
                      dpsi0_dx=psi0_r, dpsi0_dz=psi0_z, _EXTRA=extra)
   nflux=fc.psi_norm
   q=fc.q

   for i=0, fc.m-1 do begin
      for j=0, n_elements(n)-1 do begin
         d[j,i,*] = d[j,i,*]/fc.area/fc.dpsi_dchi
      end
   end 

   if(n_elements(bmnfile) ne 0) then begin
       openw, ifile, /get_lun, bmnfile

       printf, ifile, format='(3I5)', ntor, n_elements(nflux), n_elements(m)
       printf, ifile, format='(1F13.6)', nflux
       printf, ifile, format='(1I5)', transpose(m)
       for i=0, n_elements(nflux)-1 do begin
           printf, ifile, format='(2F13.6)', $
             real_part(reform(d[0,*,i])), $
             imaginary(reform(d[0,*,i]))
       end

       free_lun, ifile
   end

   mu0 = !pi*4.e-7

   if(n_elements(bmncdf) ne 0) then begin
      if(keyword_set(profdata)) then begin
         db = read_parameter('db', _EXTRA=extra)
         omega_i = flux_average('v_omega',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                                bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra,$
                               fc=fc)
         omega_e = flux_average('ve_omega',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                                bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra,$
                               fc=fc)
         w_star_i = flux_average('omega_*i',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                                bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra,$
                               fc=fc)
         F = flux_average('I',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                          bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
         p = flux_average('p',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                          bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
         pe = flux_average('pe',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                           bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
         den_e = flux_average('ne',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                              bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
         omega_ExB = omega_i - w_star_i
      end
       
      plot, [1.0, 2.3], [-1.5, 1.5], /nodata
      for i=0, n_elements(nflux) - 1, 10 do begin
         oplot, fc.r[i,*], fc.z[i,*]
      end

      bpval = fltarr(n_elements(m), n_elements(nflux))
      bpval = field_at_point(bp, x, z, fc.r, fc.z)
       
      id = ncdf_create(bmncdf, /clobber)
      if(n_elements(n) eq 1) then begin
         ncdf_attput, id, 'ntor', fix(ntor), /short, /global
      endif else begin
         ncdf_attput, id, 'ntor', 0, /short, /global
      endelse
      ncdf_attput, id, 'version', 5, /short, /global
      ; v. 4: changed definition of alpha from ~(2pi)^-4 to ~(2pi)^-2
      print, 'outputting symbol = ', symbol
      print, 'outputting units = ', units
      ncdf_attput, id, 'symbol', string(symbol), /global
      ncdf_attput, id, 'units', string(units), /global
      l_id = ncdf_dimdef(id, 'ns', n_elements(n))
      n_id = ncdf_dimdef(id, 'npsi', n_elements(nflux))
      m_id = ncdf_dimdef(id, 'mpol', n_elements(m))
      psi_norm_var = ncdf_vardef(id, 'psi_norm', [n_id], /float)
      psi_var = ncdf_vardef(id, 'psi', [n_id], /float)
      flux_pol_var = ncdf_vardef(id, 'flux_pol', [n_id], /float)
      q_var = ncdf_vardef(id, 'q', [n_id], /float)
      area_var = ncdf_vardef(id, 'area', [n_id], /float)
      current_var = ncdf_vardef(id, 'current', [n_id], /float)
      if(keyword_set(profdata)) then begin
         p_var = ncdf_vardef(id, 'p', [n_id], /float)
         F_var = ncdf_vardef(id, 'F', [n_id], /float)
         pe_var = ncdf_vardef(id, 'pe', [n_id], /float)
         ne_var = ncdf_vardef(id, 'ne', [n_id], /float)
         omega_i_var = ncdf_vardef(id, 'omega_i', [n_id], /float)
         omega_e_var = ncdf_vardef(id, 'omega_e', [n_id], /float)
         omega_e_var = ncdf_vardef(id, 'omega_ExB', [n_id], /float)
         if(keyword_set(boozer)) then begin
            if(n_elements(n) eq 1) then begin
               alpha_real_var = ncdf_vardef(id, 'alpha_real', [m_id,n_id], /float)
               alpha_imag_var = ncdf_vardef(id, 'alpha_imag', [m_id,n_id], /float)
            endif else begin
               alpha_real_var = ncdf_vardef(id, 'alpha_real', [l_id,m_id,n_id], /float)
               alpha_imag_var = ncdf_vardef(id, 'alpha_imag', [l_id,m_id,n_id], /float)

            endelse
         end
      end
      m_var = ncdf_vardef(id, 'm', [m_id], /short)
      n_var = ncdf_vardef(id, 'n', [l_id], /short)
      if(n_elements(n) eq 1) then begin
         bmn_real_var = ncdf_vardef(id, 'bmn_real', [m_id,n_id], /float)
         bmn_imag_var = ncdf_vardef(id, 'bmn_imag', [m_id,n_id], /float)
      endif else begin
         bmn_real_var = ncdf_vardef(id, 'bmn_real', [l_id,m_id,n_id], /float)
         bmn_imag_var = ncdf_vardef(id, 'bmn_imag', [l_id,m_id,n_id], /float)
      endelse
      jac_var = ncdf_vardef(id, 'jacobian', [m_id,n_id], /float)
      rpath_var = ncdf_vardef(id, 'rpath', [m_id,n_id], /float)
      zpath_var = ncdf_vardef(id, 'zpath', [m_id,n_id], /float)
      bp_var = ncdf_vardef(id, 'Bp', [m_id,n_id], /float)
      ncdf_control, id, /endef
      ncdf_varput, id, 'psi_norm', fc.psi_norm
      ncdf_varput, id, 'psi', fc.psi
      ncdf_varput, id, 'flux_pol', fc.flux_pol
      ncdf_varput, id, 'm', m
      ncdf_varput, id, 'n', n
      ncdf_varput, id, 'q', q
      ncdf_varput, id, 'area', fc.area
      ncdf_varput, id, 'current', fc.current/mu0
      if(keyword_set(profdata)) then begin
         ncdf_varput, id, 'p', reform(p)
         ncdf_varput, id, 'F', reform(F)
         ncdf_varput, id, 'pe', reform(pe)
         ncdf_varput, id, 'ne', reform(den_e)
         ncdf_varput, id, 'omega_e', reform(omega_e)
         ncdf_varput, id, 'omega_i', reform(omega_i)
         ncdf_varput, id, 'omega_ExB', reform(omega_ExB)

         if(keyword_set(boozer)) then begin
            alpha = d[0,*,*]
            for i=0, n_elements(nflux)-1 do begin
               for j=0, n_elements(m)-1 do begin
                  for k=0, n_elements(n)-1 do begin
                     alpha[k,j,i] = complex(0,1)*fc.area[i]*alpha[k,j,i] $
                                    / (m[j]*F[i] + n[k]*fc.current[i]/(2.*!pi)) $
                                    / (2.*!pi)^2
                  end
               end
            end
            if(n_elements(n) eq 1) then begin
               ncdf_varput, id, 'alpha_real', real_part(reform(alpha[0,*,*]))
               ncdf_varput, id, 'alpha_imag', imaginary(reform(alpha[0,*,*]))
            endif else begin
               ncdf_varput, id, 'alpha_real', real_part(alpha)
               ncdf_varput, id, 'alpha_imag', imaginary(alpha)
            endelse
         end
      end
      if(n_elements(n) eq 1) then begin
         ncdf_varput, id, 'bmn_real', real_part(reform(d[0,*,*]))
         ncdf_varput, id, 'bmn_imag', imaginary(reform(d[0,*,*]))
      endif else begin
         ncdf_varput, id, 'bmn_real', real_part(d)
         ncdf_varput, id, 'bmn_imag', imaginary(d)
      endelse         
      ncdf_varput, id, 'jacobian', reform(fc.j)
      ncdf_varput, id, 'rpath', fc.r
      ncdf_varput, id, 'zpath', fc.z
      ncdf_varput, id, 'Bp', reform(bpval)
      ncdf_close, id
   end
   
   if(n_elements(m_val) ne 0) then begin

       q_val = abs(m_val/float(ntor))
       indices = interpol(findgen(n_elements(q)), abs(q), q_val)

       if(n_elements(linestyle) eq 0) then linestyle=0

       for i=0, n_elements(m_val)-1 do begin
           dum = min(m-m_val[i], j, /abs)

           c = get_colors()

           if(keyword_set(phase)) then begin
               data = atan(imaginary(d[0,j,*]), real_part(d[0,j,*]))
               ytitle='!6Phase!X'
           endif else begin
               data = abs(d[0,j,*])
               ytitle=label
           endelse
           if(i eq 0 and not keyword_set(overplot)) then begin
               plot, nflux, data, color=c[0], $
                 /nodata, _EXTRA=extra, xtitle='!7W!X', ytitle=ytitle
           endif

           oplot, nflux, data, color=c[i+1], linestyle=linestyle

           fv = interpolate(nflux,indices[i])
           oplot, [fv,fv], !y.crange, color=c[i+1], linestyle=1
       end

       if(n_elements(m_val) gt 1) then begin
           names = string(format='("!8m!6 = ",I0,"!X")', m_val)
           plot_legend, names, color=c[1:n_elements(m_val)], _EXTRA=extra
       endif

       return
   endif

   if (n_elements(psi_val) ne 0) then begin
       indices = interpol(findgen(n_elements(nflux)), nflux, psi_val)
   endif else if(n_elements(q_val) ne 0) then begin
       indices = interpol(findgen(n_elements(q)), abs(q), q_val)
   endif

   if(n_elements(indices) ne 0) then begin
       print, abs(q[fix(indices)])
       print, abs(q[fix(indices+1)])

       col = fltarr(n_elements(indices))
       c = get_colors()
       for i=0, n_elements(col)-1 do begin
          col[i] = c[i mod n_elements(c)]
       end
       if(n_elements(outfile) ne 0) then begin
          openw, ifile, outfile, /get_lun
       end

       for i=0, n_elements(indices)-1 do begin
           dum = min(m-q_val[i]*ntor, j, /abs)
           dum = min(m+q_val[i]*ntor, k, /abs)

           print, 'q, Psi = ', interpolate(abs(q),indices[i]), $
             interpolate(nflux, indices[i])

           dj = interpolate(reform(d[0,j,*]), indices[i])
           dk = interpolate(reform(d[0,k,*]), indices[i])
           print, 'Resonant field: m (mag, phase) = ', m[j], abs(dj), $
             atan(imaginary(dj),real_part(dj))
           print, 'Resonant field: m (mag, phase) = ', m[k], abs(dk), $
             atan(imaginary(dk),real_part(dk))

;           if(i eq 0) then begin
;               plot, m, abs(d[*,i]), xrange=[-20,20], yrange=[0, max(abs(d))]
;           endif else begin
;               oplot, m, abs(d[*,i]), color=col[i]
;           end
           
           if(n_elements(outfile) ne 0) then begin
               print, 'q[0] = ', q[0]
               if(q[0] gt 0) then begin
                  dd = dj
                  mi = j 
               endif else begin
                  dd = dk
                  mi = k
               end
               print, 'm[mi] = ', m[mi]
               printf, ifile, format='(I5,7F12.6)', m[mi], abs(dd), $
                 atan(imaginary(dd),real_part(dd)), $
                 interpolate(nflux, indices[i]), $
                 interpolate(abs(q), indices[i]), $
                 interpolate(deriv(nflux, abs(q)), indices[i]), $
                 interpolate(fc.area, indices[i]), $
                 interpolate(deriv(nflux, fc.psi), indices[i])
           end
       end
       if(n_elements(outfile) ne 0) then begin
           free_lun, ifile
        end
;       d = dold
   endif

   if(1 eq strcmp(!d.name, 'PS', /fold_case)) then begin
       xsize = 1.33
   endif else begin
       xsize = 1.
   endelse

   xtitle='!8m!X'


   if(keyword_set(sqrtpsin)) then begin
      y = sqrt(nflux)
      ytitle='!9r!7W!X'
   endif else begin
      y = nflux
      ytitle='!7W!X'
   end

   if(n_elements(n) gt 1) then begin
      k = where(n eq ntor, ct)
      if(ct ne 1) then begin
         print, 'schaffer_plot: error: ntor = ', ntor, ' not found.'
         return
      end
      print, 'plotting ntor = ', ntor
   endif else begin
      k = 0
   endelse

   contour_and_legend, abs(d[k,*,*]), m, y,  $
                       table=39, xtitle=xtitle, ytitle=ytitle, $
                       xrange=[-20,20], yrange=[0,1], /lines, c_thick=1, $
                       ccolor=!d.table_size-1, label=label, $
                       _EXTRA=extra, xsize=xsize
   
   oplot, ntor*q, y, linestyle=2, color=!d.table_size-1

end

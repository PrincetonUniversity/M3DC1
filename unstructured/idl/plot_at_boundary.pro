pro plot_at_boundary, name, field=field, angle=ang, smooth=sm, $
                      overplot=overplot, outfile=outfile, slice=slice, $
                      noplot=noplot, xdata=x, xtitle=xtitle, ytitle=u, $
                      _EXTRA=ex

  if(n_elements(slice) gt 1) then begin
     n = slice[1] - slice[0]
     s = findgen(n) + slice[0]
     t = fltarr(n)
     for i=0, n-1 do begin
        plot_at_boundary, name, field=f, slice=s[i], smooth=sm, /noplot, $
                          xdata=x, xtitle=xtitle, ytitle=u, _EXTRA=ex
        if(i eq 0) then begin
           sz = size(f)
           field = fltarr(n, sz[1])
        end
        field[i,*] = f
        t[i] = get_slice_time(slice=s[i], _EXTRA=ex)
        
     end

     if(n_elements(outfile) eq 1) then begin
        openw, ifile, outfile, /get_lun
        for i=0, n-1 do begin
           for j=0, sz[1]-1 do begin
              printf, ifile, t[i], x[j], field[i,j]
           end
           printf, ifile, ''
        end
        free_lun, ifile
     end

     
     contour_and_legend, field, t, x, $
                         ytitle=xtitle, label=u, $
                         xtitle=make_label('t',dimensions(/t0),_EXTRA=ex), $
                         _EXTRA=ex
     return
  end
  
  xy = get_boundary_path(norm=norm, center=center, angle=angle, $
                         length=length, _EXTRA=ex)
  
  nvec = 0

  if(size(name, /type) eq 7) then begin
     if(strcmp('jnorm', name, /fold_case) eq 1) then begin

        field_r = read_field('jx', x, z, t, units=u, slice=slice, _EXTRA=ex)
        field_z = read_field('jz', x, z, t, units=u, slice=slice, _EXTRA=ex)
        
        nvec = 1

     endif else if(strcmp('bnorm', name, /fold_case) eq 1) then begin

        field_r = read_field('bx', x, z, t, units=u, slice=slice, _EXTRA=ex)
        field_z = read_field('bz', x, z, t, units=u, slice=slice, _EXTRA=ex)
        
        nvec = 1

     endif else begin
        field = read_field(name, x, z, t, slices=time, mesh=mesh, $
                           points=p, rrange=rrange, zrange=zrange, $
                           symbol=fieldname, units=u, linear=linear, $
                           mask=mask, phi=phi0, time=realtime, $
                           complex=complex, operation=op,filename=filename, $
                           linfac=linfac, fac=fac, slice=slice, _EXTRA=ex)
        if(n_elements(field) le 1) then return
        if(n_elements(units) eq 0) then units=u
     end
  endif else begin
     field = name
  endelse

  if(nvec eq 1) then begin
     field_r = reform(field_at_point(field_r, x, z, xy[0,*], xy[1,*]))
     field_z = reform(field_at_point(field_z, x, z, xy[0,*], xy[1,*]))
     field = field_r*norm[0,*] + field_z*norm[1,*]
  endif else begin
     field = reform(field_at_point(field, x, z, xy[0,*], xy[1,*]))
  end

  if(keyword_set(ang)) then begin
     x = angle
     xtitle = '!6Poloidal Angle!X'
  endif else begin
     x = length
     xtitle = '!6Length Along Boundary!X'
  endelse

  if(n_elements(sm) eq 1) then begin
     field = smooth(field, sm)
  end

  if(not keyword_set(noplot)) then begin
     if(n_elements(overplot) le 0) then begin
        plot, x, real_part(field), _EXTRA=ex, $
              xtitle=xtitle, $
              ytitle=u
     endif else begin
        oplot, x, real_part(field), _EXTRA=ex
     endelse
  end
  
  if(n_elements(outfile) eq 1) then begin
     openw, ifile, outfile, /get_lun
     for i=0, n_elements(x)-1 do printf, ifile, x[i], field[i]
     free_lun, ifile
  end
end

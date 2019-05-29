pro plot_equation, equation, cutz=cutz, _EXTRA=extra, $
                   func=func, outfile=outfile

  call_procedure, 'eqn_' + equation, nterms=nterms, names=names, term=term, $
                  title=title, x=x, z=z, _EXTRA=extra

  total = term[0,*,*]*0.
  for i=0, nterms-1 do begin
     total = total + term[i,*,*]
  end

  if(n_elements(func) eq 0) then func='real_part'
  if(n_elements(func) eq 1) then begin
     term = call_function(func,term)
     total = call_function(func,total)
  end
  
  if(n_elements(cutr) eq 1) then begin
     x0 = replicate(cutz,n_elements(z))
     z0 = z
     xtitle = 'Z'
     xdat = z
  endif else begin
     if(n_elements(cutz) eq 0) then cutz = 0
     x0 = x
     z0 = replicate(cutz,n_elements(x))
     xtitle = 'R'
     xdat = x
  endelse

  f = fltarr(nterms+1,n_elements(xdat))
  for i=0, nterms-1 do begin
     f[i,*] = field_at_point(term[i,*,*],x,z,x0,z0)
  end
  f[nterms,*] = field_at_point(total,x,z,x0,z0)

  
  window, 1
  plot, [min(xdat),max(xdat)], [min(f), max(f)], /nodata, $
        title=title, xtitle=xtitle, _EXTRA=extra

  ct3
  c = get_colors()
  for i=0, nterms-1 do begin
     oplot, xdat, f[i,*],color=c[i+1]
  end
  oplot, xdat, f[nterms,*],color=c[0], linestyle=2

  names = ['Total',names]
  plot_legend, names, color=c, _EXTRA=extra
  
  window, 0
  contour_and_legend, abs(total), x, z, /zlog
  oplot, !x.crange, [cutz, cutz]

  if(n_elements(outfile)) then begin
     openw, ifile, outfile, /get_lun

     s = strtrim(string(nterms+2))
     printf, format='('+s+'A15)', $
             ifile, [xtitle, names]
     dat = fltarr(nterms+2)
     for k=0, n_elements(xdat)-1 do begin
        dat[0] = xdat[k]
        dat[1] = f[nterms,k]
        for i=0, nterms-1 do begin
           dat[i+2] = f[i,k]
        end
        print, dat
        printf, format='('+s+'G15.5)', $
                ifile, dat
     end
     close, ifile
  end

end

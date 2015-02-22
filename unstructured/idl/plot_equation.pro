pro plot_equation, equation, cutz=cutz, _EXTRA=extra, $
                   func=func

  call_procedure, equation, nterms=nterms, names=names, term=term, $
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

  plot, [min(xdat),max(xdat)], [min(term), max(term)], /nodata, $
        title=title, xtitle=xtitle, _EXTRA=extra

  ct3
  c = get_colors()
  for i=0, nterms-1 do begin
     f = field_at_point(term[i,*,*],x,z,x0,z0)
     oplot, xdat, f,color=c[i+1]
  end
  oplot, xdat, field_at_point(total,x,z,x0,z0), linestyle=2
  names = ['Total',names]
  plot_legend, names, color=c, _EXTRA=extra
  
end

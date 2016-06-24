pro plot_boundary, _EXTRA=extra

  xy = get_boundary_path(norm=norm, center=center, angle=angle, $
                         length=length,_EXTRA=extra)

  plot, xy[0,*], xy[1,*], _EXTRA=extra

  len = max(length)

  logdl = fix(alog10(len)) - 1
  dl = 10.^logdl

  print, 'dl = ', dl

  r = sqrt((xy[0,*] - center[0])^2 + (xy[1,*] - center[1])^2)

  i = 0
  l = 0.
  ticklen = mean(r)*0.05
  while(l lt len) do begin
     x = interpol(xy[0,*], length, l)
     y = interpol(xy[1,*], length, l)
     nx = interpol(norm[0,*], length, l)
     ny = interpol(norm[1,*], length, l)

     if((i mod 10) eq 0) then begin
        ticklen = 0.05
        xyouts, x + nx*ticklen, y + ny*ticklen, string(format='(g0)', l)
     endif else begin
        ticklen = 0.02
     endelse
     oplot, x + nx*ticklen*[-1.,1.], y + ny*ticklen*[-1.,1.]
     l = l + dl
     i = i + 1
  end
end

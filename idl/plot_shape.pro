pro plot_shape, filename, names=names, _EXTRA=extra
   c = get_colors()
   n = n_elements(filename)

   for i=0, n-1 do begin
       psilim = lcfs(file=filename[i], flux0=flux0)
       levels = findgen(12)/10.*(psilim-flux0) + flux0
       if(psilim lt flux0) then levels=reverse(levels)
       psi = read_field('psi',x,z,t,slice=-1,/mks,file=filename[i], $
                       _EXTRA=extra)
       contour, psi, x, z, levels=levels, color=c[i], overplot=(i ne 0), $
          xtitle='!8R!6 (m)!X', ytitle='!8Z!6 (m)!X', _EXTRA=extra
   end
end

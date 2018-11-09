pro plot_analytic_model, r0, a, d, z0, b, overplot=op, _EXTRA=extra
   xy = analytic_model(r0, a, d, z0, b)

   if(keyword_set(op)) then begin
      oplot, xy[0,*], xy[1,*], _EXTRA=extra
   endif else begin
      plot, xy[0,*], xy[1,*], _EXTRA=extra
   endelse
end

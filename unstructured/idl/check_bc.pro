pro check_bc, time, inoslip_pol=inoslip_pol, inonormalflow=inonormalflow, filename=filename, $
               components=components, woffset=woffset, _EXTRA = ex

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(woffset) eq 0) then woffset=0
   
   xy = get_boundary_path(norm=norm, center=center, angle=angle, $
                          length=length, filename=filename, _EXTRA=ex)
   
   if(keyword_set(inoslip_pol)) then begin
     window, woffset+0, title = 'inoslip_pol multi'
     
     ; (1/R^2)*dchi/dt
     chi_x = eval_on_boundary('chi', time, x, y, filename=filename, op=2, _EXTRA=ex)
     chi_y = eval_on_boundary('chi', time, x, y, filename=filename, op=3, _EXTRA=ex)
     chi_t = -norm[1,*]*chi_x + norm[0,*]*chi_y
     Xterm = chi_t/x^2
     
     ; dU/dn
     U_x = eval_on_boundary('phi', time, x, y, filename=filename, op=2, _EXTRA=ex)
     U_y = eval_on_boundary('phi', time, x, y, filename=filename, op=3, _EXTRA=ex)
     U_n = norm[0,*]*U_x + norm[1,*]*U_y
     Uterm = x*U_n
     
     if(keyword_set(components)) then begin
       ymax = max([max(Xterm),max(Uterm)])
       ymin = min([min(Xterm),min(Uterm)])
     endif else begin
       ymax = max(Uterm+Xterm)
       ymin = min(Uterm+Xterm)
     endelse
     
     plot, angle[0:-2], Uterm[0:-2] + Xterm[0:-2], charsize=2, thick=3, yrange=[ymin,ymax], xtitle='Theta'
     if(keyword_set(components)) then begin
       oplot, angle[0:-2], Uterm[0:-2], thick=3, color=color(1)
       oplot, angle[0:-2], Xterm[0:-2], thick=3, color=color(2)
     endif
   endif
   
   if(keyword_set(inonormalflow)) then begin
     window, woffset+1, title = 'inonormalflow multi'

     ; (1/R^2)*dchi/dn
     chi_x = eval_on_boundary('chi', time, x, y, filename=filename, op=2, _EXTRA=ex)
     chi_y = eval_on_boundary('chi', time, x, y, filename=filename, op=3, _EXTRA=ex)
     chi_n = norm[0,*]*chi_x + norm[1,*]*chi_y
     Xterm = chi_n/x^2

     ; -R*dU/dt
     U_x = eval_on_boundary('phi', time, x, y, filename=filename, op=2, _EXTRA=ex)
     U_y = eval_on_boundary('phi', time, x, y, filename=filename, op=3, _EXTRA=ex)
     U_t = -norm[1,*]*U_x + norm[0,*]*U_y
     Uterm = -x*U_t

     if(keyword_set(components)) then begin
       ymax = max([max(Xterm),max(Uterm)])
       ymin = min([min(Xterm),min(Uterm)])
     endif else begin
       ymax = max(Uterm+Xterm)
       ymin = min(Uterm+Xterm)
     endelse


     plot,  angle[0:-2], Uterm[0:-2] + Xterm[0:-2], charsize=2, thick=3, yrange=[ymin,ymax], xtitle='Theta'
     if(keyword_set(components)) then begin
       oplot, angle[0:-2], Uterm[0:-2], thick=3, color=color(1)
       oplot, angle[0:-2], Xterm[0:-2], thick=3, color=color(2)
     endif
   endif

end
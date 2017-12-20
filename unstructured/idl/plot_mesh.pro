;=========================================================
; plot_mesh
; ~~~~~~~~~
;
; Plots the mesh.
; /boundary: only plot the mesh boundary
; /oplot: plots mesh as overlay to previous plot
;=========================================================
pro plot_mesh, mesh=mesh, oplot=oplot, boundary=boundary, _EXTRA=ex

   if(n_elements(mesh) eq 0) then mesh = read_mesh(_EXTRA=ex)
   if(n_tags(mesh) eq 0) then return

   nelms = mesh.nelms._data
   
   if(not keyword_set(oplot)) then begin
       xtitle = make_label('!8R!X',/l0,_EXTRA=ex)
       ytitle = make_label('!8Z!X',/l0,_EXTRA=ex)
       plot, mesh.elements._data[4,*], xtitle=xtitle, ytitle=ytitle, $
         mesh.elements._data[5,*], psym = 3, _EXTRA=ex, /nodata
   endif  

   get_normalizations, b0=b0, n0=n0, l0=l0, zeff=zeff, ion_mass=mi, _EXTRA=ex
   fac = 1.
   convert_units, fac, dimensions(/l0), b0, n0, l0, zeff, mi, _EXTRA=ex

   ct3
   col = color(9,16)
 
   version = read_parameter('version', _EXTRA=ex)
   print, 'Output version = ', version
   if(version eq 0) then begin
       xzero = read_parameter("xzero", _EXTRA=ex)
       zzero = read_parameter("zzero", _EXTRA=ex)
   endif else begin
       xzero = 0.
       zzero = 0.
   endelse

   if(keyword_set(boundary) and version ge 3) then begin
       boundary = 1
   endif else boundary = 0

   maxr = [max(mesh.elements._data[4,*]), max(mesh.elements._data[5,*])]*fac
   minr = [min(mesh.elements._data[4,*]), min(mesh.elements._data[5,*])]*fac

   czone = [1, 2, 2, 3, 3, 4, 4]

   for i=long(0), nelms-1 do begin
       a = mesh.elements._data[0,i]*fac
       b = mesh.elements._data[1,i]*fac
       c = mesh.elements._data[2,i]*fac
       t = mesh.elements._data[3,i]
       x = mesh.elements._data[4,i]*fac
       y = mesh.elements._data[5,i]*fac
       bound = fix(mesh.elements._data[6,i])

       p1 = [x, y]
       p2 = p1 + [(b+a) * cos(t), (b+a) * sin(t)]
       p3 = p1 + [b * cos(t) - c * sin(t), $
                  b * sin(t) + c * cos(t)]
       delta = 0.0
       q1 = (1.-2.*delta)*p1 + delta*p2 + delta*p3
       q2 = delta*p1 + (1.-2.*delta)*p2 + delta*p3
       q3 = delta*p1 + delta*p2 + (1.-2.*delta)*p3
      
       if(boundary) then pp=bound else pp=7
 
       if((pp and 1) eq 1) then begin
           if((bound and 1) eq 1) then begin
               izone = (bound and 120)/2^3 + 1
               c = color(czone[izone mod n_elements(czone)])
               oplot, [q1[0],q2[0]]+xzero, [q1[1],q2[1]]+zzero, $
                 color=c, thick=!p.thick*3
           end else begin
               oplot, [p1[0],p2[0]]+xzero, [p1[1],p2[1]]+zzero, $
                 color=col, thick=!p.thick/2.
           end
       end
       if((pp and 2) eq 2) then begin
           if((bound and 2) eq 2) then begin
               izone = (bound and 1920)/2^7 + 1
               c = color(czone[izone mod n_elements(czone)])
               oplot, [q2[0],q3[0]]+xzero, [q2[1],q3[1]]+zzero, $
                 color=c, thick=!p.thick*3
           end else begin
               oplot, [p2[0],p3[0]]+xzero, [p2[1],p3[1]]+zzero, $
                 color=col, thick=!p.thick/2.
           end

       end
       if((pp and 4) eq 4) then begin
           if((bound and 4) eq 4) then begin
               izone = (bound and 30720)/2^11 + 1
               c = color(czone[izone mod n_elements(czone)])
               oplot, [q3[0],q1[0]]+xzero, [q3[1],q1[1]]+zzero, $
                 color=c, thick=!p.thick*3
           end else begin
               oplot, [p3[0],p1[0]]+xzero, [p3[1],p1[1]]+zzero, $
                 color=col, thick=!p.thick/2.
           end
       end

       maxr[0] = max([maxr[0], p1[0], p2[0], p3[0]])
       minr[0] = min([minr[0], p1[0], p2[0], p3[0]])
       maxr[1] = max([maxr[1], p1[1], p2[1], p3[1]])
       minr[1] = min([minr[1], p1[1], p2[1], p3[1]])
   end

   print, 'Elements: ', nelms
   print, 'sqrt(nodes) (estimated): ', sqrt(nelms/2.)
   print, 'Width: ', maxr[0] - minr[0]
   print, 'Height: ', maxr[1] - minr[1]

end

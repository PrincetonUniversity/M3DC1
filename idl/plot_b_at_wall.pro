function discretize_path, xi, yi, ai, n
   out = fltarr(4, n)
   j = 0
   l0 = 0.
   for i=1, n_elements(xi)-1 do begin
       if(i eq n_elements(xi)-1) then begin
           m = n-j
           k = 1
       endif else begin
           m = n / (n_elements(xi)-1)
           k = 0
       endelse
       out[0,j:j+m-1] = (xi[i] - xi[i-1]) * findgen(m)/(m-k) + xi[i-1]
       out[1,j:j+m-1] = (yi[i] - yi[i-1]) * findgen(m)/(m-k) + yi[i-1]
       out[2,j:j+m-1] = ai[i-1]
       out[3,j:j+m-1] = l0 + sqrt((out[0,j:j+m-1] - xi[i-1])^2 + $
                                  (out[1,j:j+m-1] - yi[i-1])^2)
       j = j + m
       l0 = l0 + sqrt((xi[i]-xi[i-1])^2 + (yi[i] - yi[i-1])^2)
   end
   return, out
end

pro plot_b_at_wall, filename=filename, _EXTRA=extra, $
                    diiid_hfs=diiid_hfs, diiid_lfs=diiid_lfs, $
                    names=names, abs=ab, phase=phase, scale=scale, $
                    scan_phase=scan_phase, sum_fields=sum_fields, $
                    cutz=cutz, overplot=overplot,color=col

   if(n_elements(col) eq 0) then begin
       c = shift(get_colors(),-1)
   endif else begin
       c = col
   end
   if(n_elements(scale) eq 0) then scale = 1.
   if(n_elements(scale) eq 1 and n_elements(filename) gt 1) then $
       scale = replicate(scale, n_elements(filename))

   pts = 100
   if(keyword_set(diiid_hfs)) then begin
       xi = [0.9877, 0.9877]
       yi = [-1., 1.]
       ai = 90.
       ii = 1  ; use Z as X-axis
       xtitle='!8Z!6 (m)!X'
   endif else if(keyword_set(diiid_lfs)) then begin
       xi = [1.81657, 2.14931, 2.40824, 2.40824, 2.14931, 1.66285]
       yi = [-1.42342, -1.02799, -0.39916, 0.40310, 1.03195, 1.42545]
       ai = [-129.5, -112.5, -90., -67.5, -39.0]
       ii = 3                   ; use L as X-axis
       xtitle='!6Distance Along Wall (m)!X'
   endif else begin
       return
   endelse
   ai = ai*!pi/180.
   p = discretize_path(xi, yi, ai, pts)
;   plot, p[0,*], p[1,*]
;   plot, p[ii,*], p[2,*]

   b = complexarr(n_elements(filename), pts)

   for i=0, n_elements(filename)-1 do begin
       b[i,*] = read_b_at_points(p[0,*], p[1,*], p[2,*], filename=filename[i],$
                            _EXTRA=extra, /mks)
       b[i,*] = b[i,*]*1.e4*scale[i]
   end

   if(keyword_set(scan_phase)) then begin
       dphi = 2.*!pi*findgen(pts)/pts

       b0 = complexarr(1, pts, pts)
       for i=0, pts-1 do begin
           b0[0,i,*] = b[0,*] + b[1,*]*exp(complex(0., 1.)*dphi[i])
       end

       if(keyword_set(ab)) then begin
           ytitle = '!6Amplitude (G)!X'
           b0 = abs(b0)
       endif else if(keyword_set(phase)) then begin 
           ytitle = '!6Phase (deg)!X'
           b0 = -atan(imaginary(b0), real_part(b0))*180./!pi
       endif else begin
           ytitle = '!6Field at !9P!6 = 0 (G)!X'
           b0 = real_part(b0)
       end

       if(n_elements(cutz) eq 0) then begin
           contour_and_legend, b0, dphi*180./!pi, p[ii,*], $
             ytitle=xtitle, label=ytitle, xtitle='!7D!7u!6 (deg)!X', $
             _EXTRA=extra
       endif else begin
           iz = interpol(findgen(n_elements(reform(p[1,*]))), $
                          reform(p[1,*]), cutz)
           print, 'interploating at ', iz
           b1 = fltarr(pts)
           for i=0, pts-1 do begin
               b1[i] = interpolate(reform(b0[0,i,*]), iz)
           end
           plot, dphi*180./!pi, b1, $
             xtitle='!7D!7u!6 (deg)!X', ytitle=ytitle, $
             xrange=[0,360.], xstyle=1, _EXTRA=extra
       endelse
   endif else begin
       if(keyword_set(sum_fields)) then begin
           for i=1, n_elements(filename)-1 do begin
               b[0,*] = b[0,*] + b[i,*]
           end
       end

       if(keyword_set(ab)) then begin
           ytitle = '!6Amplitude (G/kA)!X'
           b = abs(b)
       endif else if(keyword_set(phase)) then begin 
           ytitle = '!6Phase (deg)!X'
           b = -atan(imaginary(b), real_part(b))*180./!pi
       endif else begin
;           ytitle = '!6Field at !9P!6 = 0 (G)!X'
;           b = real_part(b)
           bb = fltarr(1,pts,pts)
           ntor = read_parameter('ntor',filename=filename[0])
           phi = findgen(pts)*2.*!pi/pts
           for i=0, pts-1 do begin
               ; minus sign and 30 degree shift here 
               ; is to convert from M3D-C1 to DIII-D angle
               delta = 30.*!pi/180.
               bb[0,i,*] = $
                 real_part(b[0,*]*exp(-complex(0,ntor)*(phi[i]-delta)))
           end
           contour_and_legend, bb, phi*180./!pi, p[ii,*], _EXTRA=extra, $
             xtitle = '!9P!6 (deg)!X', ytitle=xtitle
           return
       end
           
       if(not keyword_set(overplot)) then begin
           plot, p[ii,*], b[0,*], /nodata, _EXTRA=extra, $
             ytitle=ytitle, xtitle=xtitle
       end
       
       for i=0, n_elements(filename)-1 do begin
           oplot, p[ii,*], b[i,*], color=c[i]
           if(keyword_set(sum_fields)) then break
       end
       
       if(not keyword_set(sum_fields)) then begin
           if(n_elements(names) eq 0 and n_elements(filename) gt 1) then $
              names = filename
       end
       if(n_elements(names) gt 0) then begin
           plot_legend, names, color=c, _EXTRA=extra
       end
   end
end

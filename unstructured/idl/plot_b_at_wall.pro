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
                    names=names, abs=ab, phase=phase, scale=scale

   c = shift(get_colors(),-1)
   if(n_elements(scale) eq 0) then scale = 1.

   pts = 100
   if(keyword_set(diiid_hfs)) then begin
       xi = [0.98, 0.98]
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

   for i=0, n_elements(filename)-1 do begin
       b = read_b_at_points(p[0,*], p[1,*], p[2,*], filename=filename[i], $
                            _EXTRA=extra, /mks)

       b = b*1.e4*scale
       if(keyword_set(ab)) then begin
           ytitle = '!6Amplitude (G/kA)!X'
           b = abs(b)
       endif else if(keyword_set(phase)) then begin 
           ytitle = '!6Phase (deg)!X'
           b = -atan(imaginary(b), real_part(b))*180./!pi
       end

       if(i eq 0) then begin
           plot, p[ii,*], b, /nodata, _EXTRA=extra, $
             ytitle=ytitle, xtitle=xtitle
       end
       oplot, p[ii,*], b, color=c[i]
   end

   if(n_elements(names) eq 0) then names = filename
   plot_legend, names, color=c, _EXTRA=extra
end

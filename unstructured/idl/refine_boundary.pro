pro refine_boundary, file, min_dist=min_dist, max_dist=max_dist, out=out, smooth=sm, _EXTRA=extra

  if(n_elements(min_dist) eq 0.) then min_dist = 0.01
  if(n_elements(max_dist) eq 0.) then max_dist = min_dist*10.

  z = read_ascii(file,data_start=1)
  x = reform(z.field1[0,*])
  y = reform(z.field1[1,*])

  plot, x, y, /iso, _EXTRA=extra

  i = 0
  while(1) do begin
     n = n_elements(x)
     print, i, n
     if(i ge n) then break
     
     if(i lt n-1) then begin
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
     endif else begin
        dx = x[0] - x[i]
        dy = y[0] - y[i]
     end

     dl = sqrt(dx^2 + dy^2)

     if(dl lt min_dist) then begin
        print, 'deleting', dl
        ; delete the next point
        if(i ge n-1) then begin
           break
        endif else if(i eq n-2) then begin
           x = x[0:i]
           y = y[0:i]
           break
        endif else begin
           x = [x[0:i], x[i+2:n-1]]
           y = [y[0:i], y[i+2:n-1]]
        endelse

     endif else if(dl gt max_dist) then begin
        print, 'adding', dl
        ; add an intermediate point
        if(i eq n-1) then begin
           x = [x[0:i], (x[0]+x[i]) / 2.]
           y = [y[0:i], (y[0]+y[i]) / 2.]
        endif else begin
           x = [x[0:i], (x[i+1]+x[i]) / 2., x[i+1:n-1]]
           y = [y[0:i], (y[i+1]+y[i]) / 2., y[i+1:n-1]]
        endelse
     endif else begin
        ; no problems; move to next point
        
        i = i+1
     end
  end


  
  if(n_elements(sm) eq 1) then begin
     x = smooth(x, sm, /edge_wrap)
     y = smooth(y, sm, /edge_wrap)
  end

  oplot, x, y, psym=4

  if(n_elements(out) eq 1) then begin
     openw, ifile, out, /get_lun
     printf, ifile, n
     for i=0, n-1 do begin
        printf, ifile, x[i], y[i]
     end
     free_lun, ifile
  end

end

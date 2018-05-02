pro refine_boundary, file, min_dist=min_dist, out=out

  if(n_elements(min_dist) eq 0.) then min_dist = 0.1

  z = read_ascii(file,data_start=1)
  x = reform(z.field1[0,*])
  y = reform(z.field1[1,*])

  plot, x, y, /iso

  while(1) do begin
     n = n_elements(x)

     dx = x
     dy = y

     for i=0, n-2 do begin
        dx[i] = x[i+1] - x[i]
        dy[i] = y[i+1] - y[i]
     end
     dx[n-1] = x[0] - x[n-1]
     dy[n-1] = y[0] - y[n-1]

     dl = sqrt(dx^2 + dy^2)
     j = where(dl gt min_dist, count)

     if(count eq 0) then break
     i = j[0]

     if(i eq n-1) then begin
        x = [x[0:i], (x[0]+x[i]) / 2.]
        y = [y[0:i], (y[0]+y[i]) / 2.]
     endif else begin
        x = [x[0:i], (x[i+1]+x[i]) / 2., x[i+1:n-1]]
        y = [y[0:i], (y[i+1]+y[i]) / 2., y[i+1:n-1]]
     endelse
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

function compensate_renorm, x, i

  y = x
  for k=0, n_elements(x)-2 do begin
     if(x[k]*x[k+1] le 0.) then continue
     if(abs(x[k+1]/x[k]) lt 1e-9) then begin
        print, 'Renormalization found at ', k
        ; Rescale
        if(k lt n_elements(x)-2) then begin
           dx = x[k+2] - x[k+1]
        endif else begin
           dx = 0.
        endelse
        for j=0, k do begin
           y[j] = y[j] * (x[k+1]-dx) / x[k]
        end
     end
  end

  return, y
end

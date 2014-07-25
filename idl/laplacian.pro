function laplacian, a, x, z, toroidal=toroidal
    axx = a
    azz = a

    for k=0, n_elements(a[*,0,0])-1 do begin
        for i=0, n_elements(a[k,0,*])-1 do begin
            axx[k,*,i] = deriv(x, deriv(x, a[k,*,i]))
        endfor

        for i=0, n_elements(a[k,*,0])-1 do begin
            azz[k,i,*] = deriv(z, deriv(z, a[k,i,*]))
        endfor

        if(keyword_set(toroidal)) then begin
            for i=0, n_elements(a[k, 0,*])-1 do begin
                axx[k,*,i] = axx[k,*,i] + deriv(x, a[k,*,i])/x
            endfor
        endif
    end

    return, axx+azz
end

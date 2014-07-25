function dz, a, z
    az = a

    for k=0, n_elements(a[*,0,0])-1 do begin
        for i=0, n_elements(a[k,*,0])-1 do begin
            az[k,i,*] = deriv(z, a[k,i,*])
        endfor
    end

    return, az
end

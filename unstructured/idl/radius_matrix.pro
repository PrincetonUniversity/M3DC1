function radius_matrix, x, z, t
    nx = n_elements(x)
    nz = n_elements(z)
    r = fltarr(n_elements(t), nx, nz)
    for k=0, n_elements(t)-1 do begin
        for j=0, nz-1 do r[k,*,j] = x
    end
    return, r
end

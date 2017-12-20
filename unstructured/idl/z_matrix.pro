function z_matrix, x, z, t
    nx = n_elements(x)
    nz = n_elements(z)
    zz = fltarr(n_elements(t), nx, nz)
    for k=0, n_elements(t)-1 do begin
        for j=0, nx-1 do zz[k,j,*] = z
    end
    return, zz
end

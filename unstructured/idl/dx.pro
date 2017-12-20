function dx, a, x
    ax = a

    for k=0, n_elements(a[*,0,0])-1 do begin
        for i=0, n_elements(a[k,0,*])-1 do begin
            ax[k,*,i] = deriv(x, a[k,*,i])
        endfor
    end

    return, ax
end


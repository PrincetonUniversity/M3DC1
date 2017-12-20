function s_bracket, a, b, x, z
    s1 = size(a, /dim)
    s2 = size(b, /dim)

    if(s1[0] gt 1 and s2[0] eq 1) then begin
        c = fltarr(s1[0], n_elements(x), n_elements(z))
        for i=0, s1[0]-1 do c[i,*,*] = s_bracket(a[i,*,*], b[0,*,*],x,z)
        return, c
    endif else if(s2[0] gt 1 and s1[0] eq 1) then begin
        c = fltarr(s2[0], n_elements(x), n_elements(z))
        for i=0, s2[0]-1 do c[i,*,*] = s_bracket(a[0,*,*], b[i,*,*],x,z)
        return, c
    endif else if(s1[0] eq 1 and s2[0] eq 1) then begin
        return, dx(a,x)*dx(b,x) + dz(a,z)*dz(b,z)
    endif else begin
        print, 'Error: sizes do not conform!'
        return, 0
    endelse
end

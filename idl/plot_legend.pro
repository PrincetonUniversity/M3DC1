pro plot_legend, names, linestyles=ls, colors=cs, left=l, top=t, psyms=p, $
                 ylog=ylog, xlog=xlog, charsize=charsize, box=box
    
    N = n_elements(names)

    if n_elements(ls) eq 0 then begin
        if(n_elements(p) eq 0) then begin
            ls = replicate(0,N)
        endif else begin
            ls = replicate(-1,N)
        endelse
    endif
    if n_elements(cs) eq 0 then begin
        cs = intarr(N)
        cs[*] = color(0,N)
    endif
    if n_elements(l) eq 0 then l=0
    if n_elements(t) eq 0 then t=0

    dx = (!x.crange(1) - !x.crange(0)) / 20.
    dy = (!y.crange(1) - !y.crange(0)) / 15.

    x = !x.crange(0) + dx * (1. + l * 20.)
    y = !y.crange(1) - dy * (1. + t * 15.)

    if(n_elements(box) ne 0) then begin
        bb = [x, y+dy, x+dx*box*20., y-N*dy]
        oplot, [bb[0], bb[0]], [bb[1], bb[3]]
        oplot, [bb[2], bb[2]], [bb[1], bb[3]]
        oplot, [bb[0], bb[2]], [bb[1], bb[1]]
        oplot, [bb[0], bb[2]], [bb[3], bb[3]]
    end

    for i=0, N-1 do begin
        d = [x, x+0.75*dx, x+1.5*dx, x+2.*dx]

        if keyword_set(xlog) then d=10.^d
        if keyword_set(ylog) then z=10.^y else z=y

        if(n_elements(ls) ne 0) then begin
            if(ls[i] ge 0) then begin
                oplot, [d(0), d(2)], [z, z], linestyle=ls(i), color=cs(i)
            endif
        endif

        if (n_elements(p) ne 0) then begin
            if (p[i] gt 0) then begin
                oplot, [d(1)], [z], color=cs(i), psym=p(i)
            endif
        endif

        xyouts, d(3), z, names(i), color=cs(i), charsize=charsize
        y = y - dy
    endfor
end

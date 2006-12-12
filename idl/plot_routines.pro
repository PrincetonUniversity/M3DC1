function colors
    if (1 EQ strcmp(!d.name, 'PS')) then begin
        TEK_COLOR
        col = [0, 2, 4, 6, 8, 10]
    endif else begin
        col = [-1, 1023, 1048575, 523776, 10737418] 
    endelse

    return, col
end

function color, c
    col = colors()
    return, col(c) 
end

pro plot_legend, names, linestyles=ls, colors=cs, left=l, top=t, psyms=p, $
                 ylog=ylog, xlog=xlog
    
    N = n_elements(names)

    if n_elements(ls) eq 0 then ls = intarr(N)
    if n_elements(cs) eq 0 then begin
        cs = intarr(N)
        cs(*) = color(0)
    endif
    if n_elements(l) eq 0 then l=0
    if n_elements(t) eq 0 then t=0

    dx = (!x.crange(1) - !x.crange(0)) / 20.
    dy = (!y.crange(1) - !y.crange(0)) / 15.

    x = !x.crange(0) + dx * (1. + l * 20.)
    y = !y.crange(1) - dy * (1. + t * 20.)

    for i=0, N-1 do begin
        d = [x, x+0.75*dx, x+1.5*dx, x+2.*dx]

        if keyword_set(xlog) then d=10.^d
        if keyword_set(ylog) then z=10.^y else z=y

        if (n_elements(p) ne 0) then begin
            if (p[i] ne 0) then begin
                oplot, [d(1)], [z], color=cs(i), psym=p(i)
            endif else begin
                oplot, [d(0), d(2)], [z, z], linestyle=ls(i), color=cs(i)
            endelse
        endif else begin
            oplot, [d(0), d(2)], [z, z], linestyle=ls(i), color=cs(i)
        endelse

        xyouts, d(3), z, names(i), color=cs(i)
        y = y - dy
    endfor
end

pro contour_and_legend, z, x, y, nlevels=nlevels, label=label, $
                        isotropic=iso, lines=lines, range=range

    z = reform(z)

    if n_elements(lines) eq 0 then lines = 0

    if n_elements(nlevels) eq 0 then begin
        if lines eq 0 then nlevels=100 else nlevels=10
    endif

    if n_elements(x) eq 0 then x = findgen(n_elements(z[*,0]))
    if n_elements(y) eq 0 then y = findgen(n_elements(z[0,*]))

    top = 1.0
    width = 0.8
    if n_elements(iso) then if iso eq 1 then begin
        top = (max(y)-min(y))/(max(x)-min(x))
    endif 
    
    if n_elements(label) eq 0 then label = ''

    if(n_elements(range) lt 2) then begin
        minval = min(z)
        maxval = max(z)
    endif else begin
        minval = range[0]
        maxval = range[1]
    endelse
    levels = findgen(nlevels)/(float(nlevels)-1.)*(maxval-minval) + minval

    !x.style = 1
    !y.style = 1

    xtitle = !x.title
    ytitle = !y.title
    !x.title = ''
    !y.title = ''

    !p.region = [width, 0.0, 1.0, top]

    if strcmp(!d.name, 'PS') then begin
        !p.region = [width, 0.0, 1.0, 1.0]
    endif

    xx = indgen(2)
    yy = levels
    zz = fltarr(2,nlevels)
    zz[0,*] = yy
    zz[1,*] = yy
    !x.range=[xx[0],xx[n_elements(xx)-1]]
    !y.range=[yy[0],yy[n_elements(yy)-1]]

    contour, zz, xx, yy, nlevels=nlevels, /fill, xtitle='', ytitle=label, $
      xticks=1, xtickname=[' ',' '], subtitle='', title=''

    !p.noerase = 1

    !x.range=[x[0],x[n_elements(x)-1]]
    !y.range=[y[0],y[n_elements(y)-1]]

    !p.region = [0.0, 0.0, width, top] 
    
    if strcmp(!d.name, 'PS') then begin
        device, xsize=10, ysize=10*top*width, /inches
        !p.region = [0.0, 0.0, width, 1.0]
    endif

    !x.title = xtitle
    !y.title = ytitle

    contour, z, x, y, /fill, levels=levels

    if(keyword_set(lines)) then begin
        contour, z, x, y, /overplot, levels=levels
    endif

    !p.region=0
    !p.noerase=0
end

pro contour_and_legend_mpeg, filename, z, x, y, nlevels=nlevels, label=label, $
                             isotropic=iso, lines=lines, range=range
    if(n_elements(range) le 2) then begin
        range = [min(z), max(z)]
    endif  

    sz = size(z)
    set_plot, 'z'
    
    mpeg_id = mpeg_open([640,480], quality=100)
  
    print, "writing mpeg..."

    for i=0, sz[1]-1 do begin

        contour_and_legend, z[i,*,*], x, y, nlevels=nlevels, label=label,$
          isotropic=iso, lines=lines, range=range

        image = rotate(tvrd(),7)

;        set_plot, 'x'
;        tv, image

        mpeg_put, mpeg_id, frame=i, image=image

    endfor

    print, "saving mpeg as ", filename

    mpeg_save, mpeg_id, filename=filename
    mpeg_close, mpeg_id

    set_plot, 'x'
end

pro begin_capture, file, xsize=xsize, ysize=ysize
    if strcmp(!d.name, 'PS') then begin 
        if n_elements(xsize) eq 0 then xsize=12.7
        if n_elements(ysize) eq 0 then ysize=10.16
        device, filename=file, /color, /encapsulated, xsize=xsize, ysize=ysize
    endif
end

pro end_capture, wait=w
    if strcmp(!d.name, 'PS') then begin 
        device, /close
    endif else if keyword_set(w) then begin
        print, 'type .cont to continue'
        stop
    endif
end

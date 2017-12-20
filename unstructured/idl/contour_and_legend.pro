; ======================================================================
; contour_and_legend
; ------------------
; 
; draws one or more contour plots with color scales
;
; z: Either a fltarr(nx,ny) of values to plot, or a fltarr(n,nx,ny) of n
;    fields to plot.  The field plots are arranged in a grid.
; x: fltarr(nx) containing the horizontal coordinate values
; y: fltarr(ny) containing vertical coordinate values
; label: a strarr(n) containing the color-bar labels for each field
; title: a strarr(n) containing the plot titles for each field
; range: a fltarr(2,n) containing the plot range for each field
; jpeg:  the name of jpeg image to be output
; _extra: passed to contour_and_legend_single
; ======================================================================

pro contour_and_legend, z, x, y, label=label, range=range, levels=levels, $
                          title=title, _EXTRA=ex, jpeg=jpeg, lines=lines, $
                        nlevels=nlevels, zlog=zlog, xsize=xsize, ysize=ysize, $
                        overplot=overplot

    if(n_elements(z) eq 0) then begin
        print, "contour_and_legend error:  nothing to plot"
        return
    end

    sz = size(z, /dim)
    if(n_elements(sz) gt 2) then begin
        n = sz[0]
    endif else begin
        n = 1
    endelse

    if(n_elements(label) eq 0) then begin
        label = strarr(n)
        label(*) = ''
    end
    if(n_elements(title) eq 0) then begin
        title = strarr(n)
        title(*) = ''
    end
    if(n_elements(zlog) eq 1) then begin
        zlog = replicate(zlog, n)
    endif else if(n_elements(zlog) eq 0) then begin
        zlog = replicate(0, n)
    endif
    if(n_elements(range) eq 0) then begin
        range = [min(z), max(z)]
    endif
    if(n_elements(lines) eq 0) then lines=0
    if(n_elements(lines) lt n) then lines = intarr(n) + lines
    if(n_elements(xsize) eq 0) then xsize = 1.
    if(n_elements(ysize) eq 0) then ysize = 1.

    rows = ceil(sqrt(n))
    cols = rows

    erase

    if(xsize gt ysize) then begin
        ysize = float(ysize)/float(xsize)
        xsize = 1.
    endif else begin
        xsize = float(xsize)/float(ysize)
        ysize = 1.
    endelse

    k = 0
    for i=0, rows-1 do begin
        for j=0, cols-1 do begin
            !p.region = [xsize*float(i)/rows,  ysize*float(j)/cols, $
                         xsize*float(i+1)/rows,ysize*float(j+1)/cols]
            if(n_elements(sz) gt 2) then zz = z[k,*,*] else zz = z

            if(n_elements(nlevels) eq 0) then begin
                contour_and_legend_single, zz, x, y, $
                  label=label[k], title=title[k], range=range, $
                  lines=lines[k], zlog=zlog[k], levels=levels, $
                  overplot=overplot, _EXTRA=ex
            endif else begin
                contour_and_legend_single, zz, x, y, levels=levels, $
                  label=label[k], title=title[k], range=range, $
                  lines=lines[k], nlevels=nlevels[k], zlog=zlog[k], $
                  overplot=overplot, _EXTRA=ex
            endelse
            !p.noerase = 1
            k = k + 1
            if(k ge n) then break
        end
        if(k ge n) then break
    end
    
    if(n_elements(jpeg) ne 0) then begin
        image = tvrd(true=1)
        print, "writing jpeg", jpeg
        write_jpeg, jpeg, image, true=1, quality=100
    end

    !p.noerase=0
    !p.region=0
end

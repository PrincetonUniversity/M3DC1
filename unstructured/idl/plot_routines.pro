function colors, maxcolors

    c = findgen(maxcolors) * !d.table_size / maxcolors

    if (1 EQ strcmp(!d.name, 'PS')) then begin
        c[0] = 0
    endif else begin
        c[0] = !d.table_size-1
    endelse
    
    return, c
end

function color, c, maxcolors
    col = colors(maxcolors)
    return, col(c) 
end

pro plot_legend, names, linestyles=ls, colors=cs, left=l, top=t, psyms=p, $
                 ylog=ylog, xlog=xlog
    
    N = n_elements(names)

    if n_elements(ls) eq 0 then ls = intarr(N)
    if n_elements(cs) eq 0 then begin
        cs = intarr(N)
        cs[*] = color(0,N)
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

pro contour_and_legend, z, x, y, label=label, range=range, $
                          title=title, _EXTRA=ex, jpeg=jpeg, lines=lines, $
                        nlevels=nlevels, zlog=zlog

    sz = size(z, /dim)
    n = sz[0]

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
        range = fltarr(2,n)
        for i=0, n-1 do range[*,i] = [min(z[i,*,*]), max(z[i,*,*])]
    endif
    if(n_elements(lines) eq 0) then lines=0
    if(n_elements(lines) lt n) then lines = intarr(n) + lines

    rows = ceil(sqrt(n))
    cols = rows

    erase

    k = 0
    for i=0, rows-1 do begin
        for j=0, cols-1 do begin
            !p.region = [float(i)/rows, float(j)/cols, $
                         float(i+1)/rows,float(j+1)/cols]

            if(n_elements(nlevels) eq 0) then begin
                contour_and_legend_single, z[k,*,*], x, y, $
                  label=label[k], title=title[k], range=range[*,k], $
                  lines=lines[k], zlog=zlog[k], _EXTRA=ex
            endif else begin
                contour_and_legend_single, z[k,*,*], x, y, $
                  label=label[k], title=title[k], range=range[*,k], $
                  lines=lines[k], nlevels=nlevels[k], zlog=zlog[k], _EXTRA=ex
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

    print, "done"

    !p.noerase=0
    !p.region=0
end


; ======================================================================
; contour_and_legend_single
; -------------------------
; 
; Draws a single contour plots with color scale.
;
; z: fltarr(nx,ny) of values to plot.
; x: fltarr(nx) containing the horizontal coordinate values
; y: fltarr(ny) containing vertical coordinate values
; label: color-bar labels
; title: plot title
; range: fltarr(2) containing the plot range
; nlevels: number of contour levels
; lines: if set, draw contour lines
; color_table: which color table to use
; isotropic: set aspect ratio so that x and y scales are equal
; _extra: passed to contour-line drawing call of contour
; ======================================================================

pro contour_and_legend_single, z, x, y, nlevels=nlevels, label=label, $
                        isotropic=iso, lines=lines, range=range, $
                        color_table=ct, zlog=zlog, _EXTRA = ex

    zed = reform(z)
    
    if(keyword_set(zlog)) then begin
        zed = zed > abs(min(zed,/abs))
    endif

    if n_elements(lines) eq 0 then lines = 0

    if n_elements(nlevels) eq 0 then begin
        if lines eq 0 then nlevels=100 else nlevels=10
    endif

    if n_elements(x) eq 0 then x = findgen(n_elements(zed[*,0]))
    if n_elements(y) eq 0 then y = findgen(n_elements(zed[0,*]))

    region = !p.region
    if(region[2] eq 0.) then region[2]=1.
    if(region[3] eq 0.) then region[3]=1.

    width = 0.8*(region[2]-region[0])
    top = region[3]-region[1]

    if(1 eq strcmp('PS', !d.name)) then begin
        device, xsize=10, ysize=8, /inches
    endif

    if keyword_set(iso) then begin
        if strcmp(!d.name, 'PS') then begin
            screen_size = [4.,3.]
        endif else begin
            device, get_screen_size=screen_size
        endelse
        aspect_ratio = (max(y)-min(y))/(max(x)-min(x)) $
          *screen_size[0]/screen_size[1]
        if(aspect_ratio le 1) then top = width*aspect_ratio $
        else begin
            width = top/aspect_ratio
            region[2] = region[0]+width/0.8
        endelse
    endif 
    
    if(width lt 0.3) then width1 = 0.3 else width1 = width
   
    if n_elements(label) eq 0 then label = ''

    if(n_elements(range) lt 2) then begin
        minval = min(zed)
        maxval = max(zed)
    endif else begin
        minval = range[0]
        maxval = range[1]
    endelse

    if(keyword_set(zlog) and minval le 0) then begin
        minval = abs(min(zed, /absolute))
    endif

    print, "maxval, minval = ", maxval, minval
    if(maxval le minval) then begin
        print, "Error in contour_and_legend: maxval <= minval."
        return
    endif

    if(keyword_set(zlog)) then begin
        levels = 10^(alog10(maxval/minval)*findgen(nlevels+1)/(float(nlevels))$
                     + alog10(minval))
    endif else begin
        levels = (maxval-minval)*findgen(nlevels+1)/(float(nlevels)) + minval
    endelse

    if(n_elements(ct) eq 0) then begin
        if(minval*maxval lt 0.) then loadct, 39 $
        else loadct, 3
    endif else loadct, ct

    charsize = !p.charsize
    !p.charsize = (region[2]-region[0]) + 0.2

    ; plot the color scale
    !p.region = [region[0]+width1, region[1], region[2], top+region[1]]
    
    xx = indgen(2)
    yy = levels
    zz = fltarr(2,nlevels+1)
    zz[0,*] = yy
    zz[1,*] = yy
    xrange=[xx[0],xx[n_elements(xx)-1]]
    yrange=[yy[0],yy[n_elements(yy)-1]]

    contour, zz, xx, yy, nlevels=nlevels, /fill, $
      ytitle=label, xtitle='', $
      xticks=1, xtickname=[' ',' '], levels=levels, title='', $
      xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, ylog=zlog
;    contour, zz, xx, yy, /overplot, nlevels=nlevels, levels=levels

    ; Plot the countours

    !p.noerase = 1

    xrange=[x[0],x[n_elements(x)-1]]
    yrange=[y[0],y[n_elements(y)-1]]

    !p.region = [region[0], region[1], width+region[0], top+region[1]]

    contour, zed, x, y, /fill, levels=levels, nlevels=nlevels, $
      xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, _EXTRA=ex

    if(keyword_set(lines)) then begin
        contour, zed, x, y, /overplot, levels=levels, nlevels=nlevels, $
          xrange=xrange, yrange=yrange, xstyle=1, ystyle=1,_EXTRA=ex
    endif

    !p.noerase = 0
    !p.charsize = charsize
    !p.region = 0
end



; ======================================================================
; contour_and_legend_mpeg
; -----------------------
; 
; draws one or more contour plots with color scales
;
; filename: filename of mpeg to output
; z: Either a fltarr(f,nx,ny) of values to plot, or a fltarr(f*n,nx,ny) of n
;    fields to plot, each having f time frames.  The fields plots are
;    arranged in a grid.  Mutliple fields may be plotted by passing 
;    z = [field1, field2, ..] where field1 etc. are fltarr(f,nx,ny).
; x: fltarr(nx) containing the horizontal coordinate values
; y: fltarr(ny) containing vertical coordinate values
; fields: number of fields (n).  IMPORTANT: if n>1, this must be set
;    for idl to properly interpret the z array.
; title: a strarr(n) containing the plot titles for each field
; range: a fltarr(2,n) containing the plot range for each field
; jpeg: if set, outputs jpegs of each frame, and writes mpeg-2 by
;    combining them.
; _extra: passed to contour_and_legend
; ======================================================================
pro contour_and_legend_mpeg, filename, z, x, y, range=range, jpeg=jpeg, $
                             fields=fields, _EXTRA=ex
    if(n_elements(fields) eq 0) then fields = 1

    sz = size(z)
    frames = sz[1]/fields

    zout = fltarr(fields,sz[2],sz[3])   

    if(n_elements(range) eq 0) then begin
        range = fltarr(2,fields)
        for f=0, fields-1 do begin
            range[*,f] = [min(z[f*frames:(f+1)*frames-1,*,*]), $
                          max(z[f*frames:(f+1)*frames-1,*,*])]
        end
    endif  
 
    if(not keyword_set(jpeg)) then begin
        mpeg_id = mpeg_open([640,480], bitrate=104857200)
    endif

    for i=0, frames-1 do begin

        print, "plotting frame", i

        for f=0, fields-1 do zout[f,*,*] = z[f*frames+i,*,*]

        if(keyword_set(jpeg)) then $
          name =  filename+string(i,format='(I4.4)')+'.jpeg'
        contour_and_legend, zout, x, y, range=range, jpeg=name, _EXTRA=ex
        
        image = tvrd(true=1)
        
        if(not keyword_set(jpeg)) then begin
            image[0,*,*] = rotate(reform(image[0,*,*]), 7)
            image[1,*,*] = rotate(reform(image[1,*,*]), 7)
            image[2,*,*] = rotate(reform(image[2,*,*]), 7)

            mpeg_put, mpeg_id, frame=i, image=image, /color
        endif

    endfor
    
    print, "saving mpeg as ", filename
    if(keyword_set(jpeg)) then begin
        make_mpg, prefix=filename, suffix='jpeg', n_start=0, n_end=sz[1]-1, $
          digits=4, mpeg_file=filename, format=1, frame_rate=1
    endif else begin
        mpeg_save, mpeg_id, filename=filename
        mpeg_close, mpeg_id
    endelse

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

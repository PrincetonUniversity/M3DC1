pro ct1
    x = findgen(64)/63.
    r = bytarr(256)
    g = bytarr(256)
    b = bytarr(256)

    b[0:63] = 255
    b[64:127] = bytscl(-x)
    r[128:191] = bytscl(x)
    r[192:255] = 255

    g[0:63] = bytscl(x)
    g[64:127] = bytscl(-x)
    g[128:191] = bytscl(x)
    g[192:255] = bytscl(-x)

    r[0] = 0
    g[0] = 0
    b[0] = 0
    r[255] = 255
    g[255] = 255
    b[255] = 255

    tvlct, r, g, b
end


pro ct2
    x = indgen(64)
    x[0] = 4
    x[1] = 4
    x[2] = 4
    x[4] = 4
    x[60] = 60
    x[61] = 60
    x[62] = 60
    x[63] = 60
    r = bytarr(256)
    g = bytarr(256)
    b = bytarr(256)

    r[0:63] = 0
    r[64:127] = bytscl(x)
    r[128:191] = 255
    r[192:255] = 255

    g[0:63] = bytscl(-x)
    g[64:127] = bytscl(x)
    g[128:191] = bytscl(-x)
    g[192:255] = bytscl(x)

    b[0:63] = 255
    b[64:127] = 255
    b[128:191] = bytscl(-x)
    b[192:255] = 0

    r[253:255] = [0, 0, 255]
    g[253:255] = [0, 0, 255]
    b[253:255] = [0, 0, 255]

    !p.color=253
    !p.background=255

    tvlct, r, g, b
end


pro ct3
    n = 10
    rgb = bytarr(n,3)

    rgb[0,*] = [  0,   0,   0]
    rgb[1,*] = [255,   0,   0]
    rgb[2,*] = [  0,   0, 255]
    rgb[3,*] = [  0, 128,   0]
    rgb[4,*] = [128,   0, 128]
    rgb[5,*] = [255, 128,   0]
    rgb[6,*] = [  0, 192, 192]
    rgb[7,*] = [192, 192,   0]
    rgb[8,*] = [ 64, 128,  64]
    rgb[9,*] = [255, 255, 255]

    dx = !d.table_size/n
    rgb_big = bytarr(!d.table_size, 3)
    rgb_big[*] = 255.
    for i=0, n-1 do begin
        for j=i*dx, (i+1)*dx-1 do rgb_big[j,*] = rgb[i,*]
    end

    tvlct, rgb_big
end


pro ct33
    x = findgen(85)
    r = bytarr(256)
    g = bytarr(256)
    b = bytarr(256)

    r[0:84] = bytscl(x)
    g[0:84] = bytscl(x)/2
    b[0:84] = bytscl(x)/2
    r[85:169] = bytscl(x)/2
    g[85:169] = bytscl(x)
    b[85:169] = bytscl(x)/2
    r[170:254] = bytscl(x)/2
    g[170:254] = bytscl(x)/2
    b[170:254] = bytscl(x)

    tvlct, r, g, b
end




function colors, maxcolors

   if(n_elements(maxcolors) eq 0) then maxcolors=10
   mc = max([maxcolors,10])

   c = indgen(mc) * !d.table_size / mc

   if (1 EQ strcmp(!d.name, 'PS')) then begin
       c[0] = 0
   endif else begin
       c[0] = !d.table_size-1
   endelse    
   
   return, c[0:maxcolors-1]
end


pro setplot, p
   if(1 eq strcmp(p, 'ps', /fold_case)) then begin
       !p.charthick=5
       !p.thick=5
       !p.charsize=1.5
   endif else if (1 eq strcmp(p, 'x', /fold_case)) then begin
       !p.charthick=1
       !p.thick=1
       !p.charsize=1.5
   endif

   set_plot, p
end

function color, c, maxcolors
    col = colors(maxcolors)
    return, col[c] 
end

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
                        nlevels=nlevels, zlog=zlog, xsize=xsize, ysize=ysize

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
                  lines=lines[k], zlog=zlog[k], levels=levels, _EXTRA=ex
            endif else begin
                contour_and_legend_single, zz, x, y, levels=levels, $
                  label=label[k], title=title[k], range=range, $
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
                               table=ct, zlog=zlog, csym=csym, $
                               nofill=nofill, noautoct=noautoct, $
                               nolegend=nolegend, color=color, levels=levels, $
                               clevels=clevels, ccolor=ccolor, $
                               _EXTRA = ex

    zed = reform(z)
   
    if(keyword_set(zlog)) then begin
        zed = zed > abs(min(zed,/abs))
    endif
    if(keyword_set(nofill)) then fill=0 else fill=1

    if(keyword_set(lines))then begin
        if(n_elements(clevels) eq 0 and n_elements(levels) eq 0) then $
          clevels = 10
    endif else lines = 0

    if n_elements(nlevels) eq 0 then nlevels=100


    if n_elements(x) eq 0 then x = findgen(n_elements(zed[*,0]))
    if n_elements(y) eq 0 then y = findgen(n_elements(zed[0,*]))

    if(1 eq strcmp('PS', !d.name)) then begin
        device, xsize=10, ysize=8, /inches
        screen_aspect = float(!d.y_size)/float(!d.x_size)
    endif else begin
        s = intarr(2)
        device, get_screen_size=s
        screen_aspect = float(s[1])/float(s[0])
    endelse

    region = !p.region
    if(region[2] eq 0.) then region[2]=1.
    if(region[3] eq 0.) then region[3]=1.

    ; width of color bar
    width1 = 0.05*(region[2]-region[0])
    ; width of margins
    lgap = 0.11*(region[2]-region[0])
    cgap = 0.18*(region[2]-region[0])
    tgap = 0.06*(region[3]-region[1])
    bgap = 0.11*(region[2]-region[0])
    rgap = 0.01*(region[2]-region[0])

    ; dimensions of plotting region
    width = region[2]-region[0] - width1 - lgap - cgap - rgap
    top = region[3]-region[1] - bgap - tgap
    charsize = !p.charsize*sqrt((region[2]-region[0])*(region[3]-region[1]))

    if keyword_set(iso) then begin
        aspect_ratio = (max(y)-min(y))/(max(x)-min(x))
        if(aspect_ratio le 1) then top = width*aspect_ratio/screen_aspect $
        else width = screen_aspect*top/aspect_ratio
    endif else aspect_ratio = top/width
      
    if n_elements(label) eq 0 then label = ''

    if(n_elements(range) lt 2) then begin
        minval = min(zed)
        maxval = max(zed)
    endif else begin
        if(range[0] ne range[1]) then begin
            minval = range[0]
            maxval = range[1]
        endif else begin
            minval = min(zed)
            maxval = max(zed)
        endelse
    endelse

    if(keyword_set(zlog) and minval le 0) then begin
        minval = abs(min(zed, /absolute))
    endif

    fracdiff = abs(maxval-minval)
    val = abs(max([maxval,minval],/absolute))
    if(val gt 0.) then fracdiff = fracdiff/val
    print, "maxval, minval, fracdiff = ", maxval, minval, fracdiff
    if(fracdiff le 1e-5) then nolegend = 1

    if(keyword_set(zlog)) then begin
        lev = 10^(alog10(maxval/minval)*findgen(nlevels+1)/(float(nlevels))$
                     + alog10(minval))
    endif else begin
        if(keyword_set(csym)) then begin
            absmax = abs(max([maxval, minval], /absolute))
            print, absmax
            maxval = absmax
            minval = -absmax
        end
        lev = (maxval-minval)*findgen(nlevels+1)/(float(nlevels)) + minval
    endelse

    if(n_elements(ct) eq 0) then begin
        if(not keyword_set(noautoct)) then begin
            if(minval*maxval lt 0.) then loadct, 39 $
            else loadct, 3
        endif
    endif else if(ct eq -1) then begin
        ct2
    endif else loadct, ct
        

    if(not keyword_set(nolegend)) then begin
    ; plot the color scale
    ; ***
    !p.position = [region[0]+width+lgap+cgap,        region[1]+bgap, $
                   region[0]+width+lgap+cgap+width1, region[1]+bgap+top]
    
    xx = indgen(2)
    yy = lev
    zz = fltarr(2,nlevels+1)
    zz[0,*] = yy
    zz[1,*] = yy
    xrange=[xx[0],xx[n_elements(xx)-1]]
    yrange=[yy[0],yy[n_elements(yy)-1]]

    contour, zz, xx, yy, nlevels=nlevels, fill=fill, $
      ytitle=label, xtitle='', charsize=charsize, $
      xticks=1, xtickname=[' ',' '], levels=lev, title='', $
      xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, ylog=zlog, $
      color=color    
;    contour, zz, xx, yy, /overplot, nlevels=nlevels, levels=levels
    ; ***

    !p.noerase = 1
    endif

    ; Plot the countours
    ;***
    xrange=[x[0],x[n_elements(x)-1]]
    yrange=[y[0],y[n_elements(y)-1]]

    !p.position = [region[0]+lgap,       region[1]+bgap, $
                   region[0]+lgap+width, region[1]+bgap+top]

    if(fracdiff le 1e-5) then begin
        plot, x, y, /nodata, isotropic=iso, charsize=charsize, $
          color=color, xrange=xrange, yrange=yrange, $
          xstyle=1, ystyle=1, _EXTRA=ex
    endif else begin
        contour, zed > minval < maxval, x, y, $
          fill=fill, levels=lev, nlevels=nlevels, $
          xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, _EXTRA=ex, $
          isotropic=iso, charsize=charsize, color=color
    endelse

    print, 'lines = ', lines
    if(keyword_set(lines) and fill eq 1) then begin
        if(n_elements(levels) eq 0) then begin 
            if(keyword_set(zlog)) then begin
                levels = 10^(alog10(maxval/minval)* $
                          findgen(clevels+1)/(float(clevels)) $
                          + alog10(minval))
            endif else begin
                levels = (maxval-minval)*findgen(clevels+1)/(float(clevels)) $
                  + minval
            endelse
        end
        print, 'levels = ', levels
        if(n_elements(clevels) ne 0) then print, 'clevels= ', clevels

        if(n_elements(ccolor) ne 0) then begin
            loadct,12
        endif

        contour, zed, x, y, /overplot, levels=levels, nlevels=clevels, $
          xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
          _EXTRA=ex, charsize=charsize, color=ccolor
    endif
    ;***

    !p.noerase = 0
    !p.position = 0
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


pro plot_slice, data, x, y, z, value=value, normal=normal, itor=itor, $
                range=range, reject=reject

   if(n_elements(reject) eq 0) then reject = 1

   if(n_elements(normal) gt 0) then begin
       n = n_elements(value)
       for m=long(0), n-1 do begin
           xdata = data
           for i=0, n_elements(x)-1 do begin
               for j=0, n_elements(y)-1 do begin
                   for k=0, n_elements(z)-1 do begin
                       xdata[i,j,k] = $
                         normal[0,m]*x[i] + normal[1,m]*y[j] + normal[2,m]*z[k]
                   end
               end
           end
           isosurface, xdata, value[m], v1, p1
           
           if(n_elements(v1) le 3) then begin
               print, 'Error: surface does not intersect domain'
               return
           endif
           
           if(m eq 0) then begin
               v = v1
               p = p1
           endif else begin
               for i=long(0), n_elements(p1)-1, 4 do begin
                   p1[i+1] = p1[i+1] + n_elements(v[0,*])
                   p1[i+2] = p1[i+2] + n_elements(v[0,*])
                   p1[i+3] = p1[i+3] + n_elements(v[0,*])
               end
               p = [p, p1]
               v = [[v], [v1]]
           end
       end

       aout = interpolate(data, reform(v[0,*]), reform(v[1,*]), reform(v[2,*]))
       set_shading, reject=0

       if(n_elements(range) eq 2) then begin
           aout = (aout > range[0]) < range[1]
       end
       shades = bytscl(aout)
   endif else begin
       nx = lonarr(n_elements(value))
       for i=0, n_elements(value)-1 do begin
           val = max(data)*value[i]

           isosurface, data, val, v0, p0
           nx[i] = n_elements(v0[0,*])
           if(i eq 0) then begin
               v = v0
               p = p0
           endif else begin
               n = n_elements(v[0,*])
               e = n_elements(p0)/4
               for j=long(0), e-1 do begin
                   p0[j*4+1:j*4+3] = p0[j*4+1:j*4+3] + n
               end
               p = [p, p0]
               v = transpose([transpose(v), transpose(v0)])
           endelse
       end
       set_shading, reject=reject, light=[0.2,0.5,0.]
       aout = 0.
   endelse
       
   if(n_elements(v) le 3) then begin
       print, 'Error: surface does not intersect domain'
       return
   endif

   rr = (max(x) - min(x))*v[0,*]/(n_elements(x)-1.) + min(x)
   yy = (max(y) - min(y))*v[1,*]/(n_elements(y)-1.) + min(y)
   zz = (max(z) - min(z))*v[2,*]/(n_elements(z)-1.) + min(z)

   q = v

   if(keyword_set(itor)) then begin
       q[0,*] = rr*cos(yy*!pi/180.)
       q[1,*] = rr*sin(yy*!pi/180.)
       q[2,*] = zz
   endif else begin
       q[0,*] = rr
       q[1,*] = yy
       q[2,*] = zz
   endelse

   if(n_elements(normal) eq 0) then begin
                                ; reverse polygons
       pr = p
       e = n_elements(p)/4
       for i=long(0), e-1 do begin
           pr[i*4+1:i*4+3] = reverse(p[i*4+1:i*4+3])
       end
       
       if(n_elements(brightness) eq 0) then brightness=0.2
       if(n_elements(contrast) eq 0) then contrast=1.
       if(n_elements(light_brightness) eq 0) then light_brightness=0.8
       if(n_elements(specular) eq 0) then specular=0.1
       
       light_dir = [0., 1., -1.]/sqrt(2.)
       normals = compute_mesh_normals(q, p)
       for i=long(0), n_elements(q[0,*])-1 do begin
           g = transpose(!p.t) # [normals[*,i], 0]
           if(g[2] lt 0) then normals[*,i] = -normals[*,i]
       end
       reflect = -(light_dir[0]*normals[0,*] + $
                   light_dir[1]*normals[1,*] + $
                   light_dir[2]*normals[2,*]) > 0
       
       shades = reform(bytscl(aout) * contrast + brightness*255 $
                       + light_brightness*reflect*255 $
                       + specular*exp((reflect-1.)*5.)*255< 255)
       
       if(n_elements(nx) ne 0) then begin
           shades = shades / n_elements(nx)
           e = n_elements(shades)-1
           j = nx[0]
           for i=1, n_elements(nx)-1 do begin
               shades[j:e] = shades[j:e] + 255/n_elements(nx)
               j = j + nx[i]
           end
       end
       ppp = [pr, p]
   endif else begin
       ppp = p
   endelse

   tv, polyshade(q, ppp, /t3d, shades=reform(shades))
end

pro split_plot, x, y, xrange=xrange, yrange=yrange, xlog=xlog, xtitle=xtitle, $
                linestyle=linestyle,  _EXTRA=extra
   sz = size(x)

   if(n_elements(xtitle) ne 0) then begin
       xtp = xtitle
       xtn = '!6-!X' + xtitle
   endif else begin
       xtp = ''
       xtn = ''
   endelse

   if(n_elements(yrange) eq 0) then begin
       yrange = [min(y), max(y)]
   end

   if(keyword_set(xlog)) then begin
       !p.multi = [0,2,1]
       if(n_elements(xrange) eq 0) then begin
           xrange = [min(abs(x)), max(abs(x))]
       end

       if(sz[0] eq 1) then begin
           j = where(x lt 0, count)
           if(count gt 0) then begin
               plot, -x[j], y[j], xrange=reverse(xrange), yrange=yrange, $
                 /xlog, xtitle=xtn, _EXTRA=extra
           end
           j = where(x gt 0, count)
           if(count gt 0) then begin
               plot, x[j], y[j], xrange=xrange, yrange=yrange, $
                 /xlog, xtitle=xtp, _EXTRA=extra
           end      
       endif else begin
           if(n_elements(linestyle) eq 0) then begin
               linestyle = replicate(0, sz[2])
           end
           c = shift(colors(),-1)
           plot, [0,0], [0,0], /nodata, xrange=reverse(xrange), yrange=yrange, $
             /xlog, xtitle=xtn, _EXTRA=extra
           for i=0, sz[2]-1 do begin
               j = where(x[*,i] lt 0, count)
               print, 'j = ', j
               if(count gt 0) then begin
                   oplot, -x[j,i], y[j,i], color=c[i], linestyle=linestyle[i], $
                     _EXTRA=extra
               end
           end
           plot, [0,0], [0,0], /nodata, xrange=xrange, yrange=yrange, $
             /xlog, xtitle=xtp, _EXTRA=extra
           for i=0, sz[2]-1 do begin
               j = where(x[*,i] gt 0, count)
               if(count gt 0) then begin
                   oplot, x[j,i], y[j,i], color=c[i], linestyle=linestyle[i], $
                     _EXTRA=extra
               end
           end
       endelse
       
   endif else begin

       if(n_elements(xrange) eq 0) then begin
           xrange = [min(x), max(x)]
       end

       if(sz[0] eq 1) then begin
           plot, x, y, xrange=xrange, yrange=yrange, $
             xtitle=xtp, _EXTRA=extra
       endif else begin
           if(n_elements(linestyle) eq 0) then begin
               linestyle = replicate(0, sz[2])
           end
           c = shift(colors(),-1)

           plot, [0,0], [0,0], /nodata, xrange=xrange, yrange=yrange, $
             xtitle=xtp, _EXTRA=extra
           for i=0, sz[2]-1 do begin
               oplot, x[*,i], y[*,i], color=c[i], linestyle=linestyle[i], $
                 _EXTRA=extra
           end
       endelse
   endelse

   !p.multi = 0
end

;=====================================================================
; Functions for finding derivatives of fields
;=====================================================================
function laplacian, a, x, z, toroidal=toroidal
    axx = a
    azz = a

    for k=0, n_elements(a[*,0,0])-1 do begin
        for i=0, n_elements(a[k,0,*])-1 do begin
            axx[k,*,i] = deriv(x, deriv(x, a[k,*,i]))
        endfor

        for i=0, n_elements(a[k,*,0])-1 do begin
            azz[k,i,*] = deriv(z, deriv(z, a[k,i,*]))
        endfor

        if(keyword_set(toroidal)) then begin
            for i=0, n_elements(a[k, 0,*])-1 do begin
                axx[k,*,i] = axx[k,*,i] + deriv(x, a[k,*,i])/x
            endfor
        endif
    end

    return, axx+azz
end

function grad_shafranov, a, x, z, toroidal=toroidal
    axx = a
    azz = a

    for k=0, n_elements(a[*,0,0])-1 do begin
        for i=0, n_elements(a[k,0,*])-1 do begin
            axx[k,*,i] = deriv(x, deriv(x, a[k,*,i]))
        endfor

        for i=0, n_elements(a[k,*,0])-1 do begin
            azz[k,i,*] = deriv(z, deriv(z, a[k,i,*]))
        endfor

        if(keyword_set(toroidal)) then begin
            for i=0, n_elements(a[k, 0,*])-1 do begin
                axx[k,*,i] = axx[k,*,i] - deriv(x, a[k,*,i])/x
            endfor
        endif

    end

    return, axx+azz
end

function dx, a, x
    ax = a

    for k=0, n_elements(a[*,0,0])-1 do begin
        for i=0, n_elements(a[k,0,*])-1 do begin
            ax[k,*,i] = deriv(x, a[k,*,i])
        endfor
    end

    return, ax
end

function dz, a, z
    az = a

    for k=0, n_elements(a[*,0,0])-1 do begin
        for i=0, n_elements(a[k,*,0])-1 do begin
            az[k,i,*] = deriv(z, a[k,i,*])
        endfor
    end

    return, az
end

function a_bracket, a, b, x, z
    return, -dx(a,x)*dz(b,z) + dz(a,z)*dx(b,x)
end

function s_bracket, a, b, x, z
    return, dx(a,x)*dx(b,x) + dz(a,z)*dz(b,z)
end


function aa_bracket, a, b, c, x, z
    ax = dx(a,x)
    az = dz(a,z)
    bxx = dx(dx(b,x),x)
    bxz = dx(dz(b,z),x)
    bzz = dz(dz(b,z),z)
    cx = dx(c,x)
    cz = dz(x,z)

    return, ax*bxz*cz - ax*bzz*cx - az*bxx*cz + az*bxz*cx
end


function as_bracket, a, b, c, x, z
    ax = dx(a,x)
    az = dz(a,z)
    bxx = dx(dx(b,x),x)
    bxz = dx(dz(b,z),x)
    bzz = dz(dz(b,z),z)
    cx = dx(c,x)
    cz = dz(x,z)

    return, ax*bxz*cx + ax*bzz*cz - az*bxx*cx - az*bxz*cz
end

function la_bracket, a, b, x, z
    return, 0.5*(laplacian(a_bracket(b,a,x,z),x,z) $
                 -a_bracket(a, laplacian(b,x,z),x,z) $
                 -a_bracket(b, laplacian(a,x,z),x,z))
end

function ls_bracket, a, b, x, z
    return, 0.5*(laplacian(s_bracket(b,a,x,z),x,z) $
                 -s_bracket(a, laplacian(b,x,z),x,z) $
                 -s_bracket(b, laplacian(a,x,z),x,z)) $
      -laplacian(a,x,z)*laplacian(b,x,z)
end

function radius_matrix, x, z, t
    nx = n_elements(x)
    nz = n_elements(z)
    r = fltarr(n_elements(t), nx, nz)
    for k=0, n_elements(t)-1 do begin
        for j=0, nz-1 do r[k,*,j] = x
    end
    return, r
end

;==================================================================
; Functions for reading the C1 hdf5 output
;==================================================================

;==============================================
; hdf5_file_test
; ~~~~~~~~~~~~~~
; 
; returns 1 if "filename" is a valid hdf5 file;
; otherwise returns 0
; =============================================
function hdf5_file_test, filename
   if(file_test(filename) eq 0) then begin
       print, "Error: ", filename, " is not a valid file."
       return, 0
   endif else if(h5f_is_hdf5(filename) eq 0) then begin
       print, "Error: ", filename, " is not a valid HDF5 file."
       return, 0
   endif else return, 1
end

;================================================
; read_attribute
; ~~~~~~~~~~~~~~
;
; returns the value associated with the attribute
; "name" in "group_id"
;================================================
function read_attribute, group_id, name
   nmembers = h5a_get_num_attrs(group_id)

   for i=0, nmembers-1 do begin
       attr_id = h5a_open_idx(group_id, i)
       if(strcmp(name, h5a_get_name(attr_id), /fold_case)) then begin
           result = h5a_read(attr_id)
           h5a_close, attr_id
           return, result
       endif
       h5a_close, attr_id
   end
   
   print, "Error: cannot find attribute ", name

   return, 0
end

;==================================================
; read_parameter
; ~~~~~~~~~~~~~~
;
; returns the value associated with the attribute
; "name" in the root group of "filename"
;==================================================
function read_parameter, name, filename=filename, print=pr
   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(hdf5_file_test(filename) eq 0) then return, 0

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   attr = read_attribute(root_id, name)
   h5g_close, root_id
   h5f_close, file_id

   if(keyword_set(pr)) then print, name, " = ", attr

   return, attr
end


;=====================================================
; read_scalars
; ~~~~~~~~~~~~
;
; returns a structure populated with the data from
; the "scalars" group of "filename"
;=====================================================
function read_scalars, filename=filename
   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(hdf5_file_test(filename) eq 0) then return, 0

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   scalars = h5_parse(root_id, "scalars", /read_data)
   h5g_close, root_id
   h5f_close, file_id

   return, scalars
end

function time_name, t
   label = string(FORMAT='("time_",I3.3)', t)
   return, label
end


pro plot_mesh, mesh, color=col, linestyle=lin, oplot=oplot, $
               filename=filename, _EXTRA=ex
   nelms = mesh.nelms._data
   
   if(not keyword_set(oplot)) then begin
       plot, mesh.elements._data[4,*], $
         mesh.elements._data[5,*], psym = 3, _EXTRA=ex
   endif  

   if(n_elements(col) ne 0) then loadct, 12
  
   xzero = read_parameter("xzero", filename=filename)
   zzero = read_parameter("zzero", filename=filename)

   for i=0, nelms-1 do begin
       a = mesh.elements._data[0,i]
       b = mesh.elements._data[1,i]
       c = mesh.elements._data[2,i]
       t = mesh.elements._data[3,i]
       x = mesh.elements._data[4,i]
       y = mesh.elements._data[5,i]

       p1 = [x, y]
       p2 = p1 + [(b+a) * cos(t), (b+a) * sin(t)]
       p3 = p1 + [b * cos(t) - c * sin(t), $
                  b * sin(t) + c * cos(t)]
       
       oplot, [p1[0],p2[0]]+xzero, [p1[1],p2[1]]+zzero, $
         color=col, linestyle=lin, thick=.2
       oplot, [p2[0],p3[0]]+xzero, [p2[1],p3[1]]+zzero, $
         color=col, linestyle=lin, thick=.2
       oplot, [p3[0],p1[0]]+xzero, [p3[1],p1[1]]+zzero, $
         color=col, linestyle=lin, thick=.2
   end
end

;======================================================
; is_in_tri
; ~~~~~~~~~
;
; returns 1 if the local coordinates "localp" descibes
; is within the triangle described by a, b, c;
; otherwise, returns 0
;======================================================
function is_in_tri, localp, a, b, c

   small = (a+b+c)*1e-6

   if(localp[1] lt 0. - small) then return, 0
   if(localp[1] gt c + small) then return, 0
   
   x = 1.-localp[1]/c
   if(localp[0] lt -b*x - small) then return, 0
   if(localp[0] gt a*x + small) then return, 0

   return, 1
end

;==============================================================
; eval
; ~~~~
;
; given avec field data "field" for element "elm", the value
; of the field at the local position "localpos"
;==============================================================
function eval, field, localpos, elm, operation=op

   mi = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0]
   ni = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5]
   sum = 0.

   if(n_elements(op) eq 0) then op = 1

   for p=0, 19 do begin
       case op of
           1: sum = sum + field[p,elm]*(localpos[0]^mi[p]*localpos[1]^ni[p])
           7: begin
               if(mi[p] ge 2) then $
                 sum = sum + field[p,elm]* $
                 mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
               if(ni[p] ge 2) then $
                 sum = sum + field[p,elm]* $
                 ni[p]*(ni[p]-1)*localpos[1]^(ni[p]-2)*localpos[0]^mi[p]
           end
       end
   end

   return, sum
end


function eval_field, field, mesh, r=xi, z=yi, points=p, operation=op, $
                     filename=filename, xrange=xrange, yrange=yrange

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nelms = mesh.nelms._data 

   if(n_elements(p) eq 0) then p = 100
   if(n_elements(op) eq 0) then op = 1
   
   minpos = fltarr(2)
   maxpos = fltarr(2)
   pos = fltarr(2)
   localpos = fltarr(2)
   index = intarr(2)

   xzero = read_parameter("xzero", filename=filename)
   zzero = read_parameter("zzero", filename=filename)

   if(n_elements(xrange) lt 2) then $
     xrange = [0.,mesh.width._data] + xzero
   if(n_elements(yrange) lt 2) then $
     yrange = [0.,mesh.height._data] + zzero

   dx = (xrange[1] - xrange[0]) / (p - 1.)
   dy = (yrange[1] - yrange[0]) / (p - 1.)

   result = fltarr(p,p)

   small = 1e-3

   for i=0,nelms-1 do begin
       a = mesh.elements._data[0,i]
       b = mesh.elements._data[1,i]
       c = mesh.elements._data[2,i]
       t = mesh.elements._data[3,i]
       x = mesh.elements._data[4,i]
       y = mesh.elements._data[5,i]
       co = cos(t)
       sn = sin(t)

       p1 = [x, y]
       p2 = p1 + [(b+a)*co, (b+a)*sn]
       p3 = p1 + [b*co - c*sn, b*sn + c*co]

       minpos = [min([p1[0], p2[0], p3[0]]), min([p1[1], p2[1], p3[1]])]
       maxpos = [max([p1[0], p2[0], p3[0]]), max([p1[1], p2[1], p3[1]])]
             
       if(maxpos[0] lt xrange[0]-xzero) then continue
       if(maxpos[1] lt yrange[0]-zzero) then continue

       index[1] = (minpos[1]-yrange[0]+zzero)/dy
       pos[1] = index[1]*dy+yrange[0]-zzero

       while(pos[1] le maxpos[1] + small*dy) do begin
           index[0] = (minpos[0]-xrange[0]+xzero)/dx
           pos[0] = index[0]*dx+xrange[0]-xzero

           if (index[1] ge 0) and (index[1] lt p) then begin

               while(pos[0] le maxpos[0] + small*dx) do begin
                   localpos = [(pos[0]-x)*co + (pos[1]-y)*sn - b, $
                              -(pos[0]-x)*sn + (pos[1]-y)*co]

                   if (index[0] ge 0) and (index[0] lt p) then begin
                       if(is_in_tri(localpos,a,b,c) eq 1) then begin
                           result[index[0], index[1]] = $
                             eval(field, localpos, i, op=op)
                       endif
                   endif
                                  
                   pos[0] = pos[0] + dx
                   index[0] = index[0] + 1
               end
           endif
           
           pos[1] = pos[1] + dy
           index[1] = index[1] + 1
       end
   end

   xi = findgen(p)*(xrange[1]-xrange[0])/(p-1.) + xrange[0]
   yi = findgen(p)*(yrange[1]-yrange[0])/(p-1.) + yrange[0]

   return, result
end

function read_raw_field, name, time, mesh=mesh, filename=filename, time=t
   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(hdf5_file_test(filename) eq 0) then return, 0

   nt = read_parameter("ntime", filename=filename)

   if(time ge nt) then begin
       print, "Error: there are only ", nt-1, " time slices."
       return, 0
   endif

   file_id = h5f_open(filename)

   time_group_id = h5g_open(file_id, time_name(time))
   mesh = h5_parse(time_group_id, 'mesh', /read_data)   

   field_group_id = h5g_open(time_group_id, 'fields')
   field = h5_parse(field_group_id, name, /read_data)

   time_id = h5a_open_name(time_group_id, "time")
   t = h5a_read(time_id)
   h5a_close, time_id

   h5g_close, field_group_id
   h5g_close, time_group_id
   h5f_close, file_id

   return, field._data
end



function translate, name, units=units
   units = ''
   va0 = '!8v!DA!60!N!X'
   b0 = '!8B!D!60!N!X'
   n0 = '!8n!D!60!N!X'
   t0 = '!7s!DA!60!N!X'
   l0 = '!8L!X'
   pi4 = '!64!7p!X'
   sq = '!U!62!N!X'
   cu = '!U!63!N!X'
   c = '!8c!X'

   if(strcmp(name, 'psi', /fold_case) eq 1) then begin
       units = b0+l0+sq
       return, "!7w!X"
   endif else if(strcmp(name, 'phi', /fold_case) eq 1) then begin
       units = va0+l0+sq
       return, "!8U!X"
   endif else if(strcmp(name, 'chi', /fold_case) eq 1) then begin
       units = va0+l0
       return, "!7v!X"
   endif else if(strcmp(name, 'jphi', /fold_case) eq 1) then begin
       units = b0
       return, "!7D!6!U*!N!7w!x"
   endif else if(strcmp(name, 'vor', /fold_case) eq 1) then begin
       units = va0
       return, "!7D!6!U*!N!8U!x"
   endif else if(strcmp(name, 'com', /fold_case) eq 1) then begin
       return, "!9G.!17v!X"
   endif else if(strcmp(name, 'eta', /fold_case) eq 1) then begin
       units = pi4+l0+va0+'/'+c+sq
       return, "!7g!X"
   endif else if(strcmp(name, 'den', /fold_case) eq 1) then begin
       units = n0
       return, "!8n!X"
   endif else if(strcmp(name, 'p', /fold_case) eq 1) then begin
       units = b0+sq+'/'+pi4
       return, "!8p!X"
   endif else if(strcmp(name, 'pe', /fold_case) eq 1) then begin
       units = b0+sq+'/'+pi4
       return, "!8p!De!N!X"
   endif else if(strcmp(name, 'beta', /fold_case) eq 1) then begin
       return, "!7b!X"
   endif else if(strcmp(name, 'thermal velocity', /fold_case) eq 1) or $
     (strcmp(name, 'vt', /fold_case) eq 1) then begin
       units = va0
       return, "!8v!Dt!N!X"
   endif else if(strcmp('sound speed', name, /fold_case) eq 1) or $
     (strcmp('cs', name, /fold_case) eq 1) then begin
       units = va0
       return, "!8c!Ds!N!X"       
   endif else if(strcmp(name, 'toroidal velocity', /fold_case) eq 1) or $
     (strcmp(name, 'vz', /fold_case) eq 1) then begin
       units = va0
       return, "!8v!D!9P!N!X"
   endif else if(strcmp(name, 'temperature', /fold_case) eq 1) or $
     (strcmp(name, 't', /fold_case) eq 1) then begin
       units = b0+sq+'/'+pi4+n0
       return, "!8T!X"
   endif else if(strcmp(name, 'angular momentum', /fold_case) eq 1) or $
     (strcmp(name, 'lz', /fold_case) eq 1) then begin
       units = n0+va0
       return, "!8L!D!9P!N!X"
   endif else if(strcmp(name, 'jz', /fold_case) eq 1) or $
     (strcmp(name, 'lz', /fold_case) eq 1) then begin
       units = c+b0+'/'+pi4
       return, "!8J!D!9P!N!X"
   endif else if(strcmp(name, 'sigma', /fold_case) eq 1) then begin
       units = n0+'/'+t0
       return, "!7r!X"
   endif else if(strcmp(name, 'minor radius', /fold_case) eq 1) or $
     (strcmp(name, 'r', /fold_case) eq 1) then begin
       units = l0
       return, "!8r!X"
   endif else if(strcmp(name, 'kappa', /fold_case) eq 1) then begin
       units = n0 + l0+'!U2!N!X'+'/'+t0
       return, "!7j!X"
   endif else if(strcmp(name, 'energy density', /fold_case) eq 1) then begin
       units = b0+sq+'/'+pi4
       return, "!7e!X"
   endif else if(strcmp(name, 'energy', /fold_case) eq 1) then begin
       units = b0+sq+l0+cu+'/'+pi4
       return, "!8E!X"
   endif  

   return, "!8" + name
end


function read_field, name, x, y, t, slices=time, mesh=mesh, filename=filename,$
                     points=pts, rrange=xrange, zrange=yrange, $
                     h_symmetry=h_symmetry, v_symmetry=v_symmetry, diff=diff

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(keyword_set(diff)) then begin
       if (n_elements(filename) eq 1) then filename = [filename, 'C1.h5']
       diff_filename = filename[1]
       filename = filename[0]
   endif
   if(n_elements(pts) eq 0) then pts = 50

   if(hdf5_file_test(filename) eq 0) then return, 0

   nt = read_parameter("ntime", filename=filename)
   nv = read_parameter("numvar", filename=filename)
   itor = read_parameter("itor", filename=filename)

   if(n_elements(time) eq 0) then begin
       trange = [0,nt-1]
   endif else if(n_elements(time) eq 1) then begin
       trange = [time, time]
   endif

   if((trange[0] ge nt) or (trange[1] ge nt)) then begin
       print, "Error: there are only ", nt-1, " time slices."
       return, 0
   endif

   data = fltarr(trange[1]-trange[0]+1, pts, pts)

   ;==========================================
   ; local_beta = 2*P/B^2
   ;==========================================
   if(strcmp('beta', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(n_elements(psi) le 1) then return, 0

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.

       if(nv ge 2) then begin
           I = read_field('I',x,y,t,slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange)
       endif else begin
           I = read_parameter('bzero',filename=filename) $
             * read_parameter('xmin', filename=filename)
       endelse
      
       if(nv ge 3) then begin
           P = read_field('P',x,y,t,slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange)
       endif else begin
           P = read_parameter('p0',filename=filename)
       endelse

       b2 = (s_bracket(psi,psi,x,y) + i^2)/r^2

       beta = 2.*P/b2

       data = beta

   ;===========================================
   ; toroidal field
   ;===========================================
   endif else if(strcmp('toroidal field', name, /fold_case) eq 1) or $
     (strcmp('bz', name, /fold_case) eq 1) then begin
       
       if(nv lt 2) then begin
           print, "numvar < 2"
           return, 0
       endif

       I = read_field('I',x,y,t,slices=time, mesh=mesh, filename=filename, $
                      points=pts, rrange=xrange, zrange=yrange)
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
   
       data = I/r

   ;===========================================
   ; toroidal velocity
   ;===========================================
   endif else if(strcmp('toroidal velocity', name, /fold_case) eq 1) or $
     (strcmp('vz', name, /fold_case) eq 1) then begin
       
       if(nv lt 2) then begin
           print, "numvar < 2"
           return, 0
       endif

       v = read_field('V',x,y,t,slices=time, mesh=mesh, filename=filename, $
                        points=pts,rrange=xrange,zrange=yrange)
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
   
       data = v/r

   ;===========================================
   ; thermal velocity
   ;===========================================
   endif else if(strcmp('thermal velocity', name, /fold_case) eq 1) or $
     (strcmp('vt', name, /fold_case) eq 1) then begin

       idens = read_parameter('idens', filename=filename)

       if(nv ge 3) then begin
           P = read_field('P',x,y,t,slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange)
       endif else begin
           P = read_parameter('p0', filename=filename)
       endelse

       if(idens eq 1) then begin
           n = read_field('den', x,y,t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       endif else begin
           n = 1.
       endelse
  
       data = sqrt(P/n)

   ;===========================================
   ; sound speed
   ;===========================================
   endif else if(strcmp('sound speed', name, /fold_case) eq 1) or $
     (strcmp('cs', name, /fold_case) eq 1) then begin

       vt = read_field('vt',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       gam = read_parameter('gam', filename=filename)
  
       data = sqrt(gam)*vt

   ;===========================================
   ; temperature
   ;===========================================
   endif else if(strcmp('temperature', name, /fold_case) eq 1) or $
     (strcmp('t', name, /fold_case) eq 1) then begin

       idens = read_parameter('idens', filename=filename)

       if(nv ge 3) then begin
           P = read_field('P', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange)
       endif else begin
           P = read_parameter('p0', filename=filename)
       endelse

       if(idens eq 1) then begin
           n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange)
       endif else begin
           n = 1.
       endelse
  
       data = p/n

   ;===========================================
   ; angular momentum
   ;===========================================
   endif else if(strcmp('angular momentum', name, /fold_case) eq 1) or $
     (strcmp('lz', name, /fold_case) eq 1) then begin

       idens = read_parameter('idens', filename=filename)

       V = read_field('V', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       if(idens eq 1) then begin
           n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange)
       endif else begin
           n = 1.
       endelse
  
       data = n*v

   ;===========================================
   ; toroidal current
   ;===========================================
   endif else if(strcmp('jz', name, /fold_case) eq 1) then begin

       idens = read_parameter('idens', filename=filename)

       jphi = read_field('jphi', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.

       data = -jphi/r

   ;===========================================
   ; minor radius
   ;===========================================
   endif else if(strcmp('minor radius', name, /fold_case) eq 1) or $
     (strcmp('r', name) eq 1) then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       nulls, time, psi=psi, xpoint=xpoint, axis=axis, r=x, z=y, $
         filename=filename, _EXTRA=extra
      
       x0 = x[axis[0,0]]
       z0 = y[axis[0,0]]
       xx = fltarr(n_elements(t),n_elements(x),n_elements(y))
       zz = fltarr(n_elements(t),n_elements(x),n_elements(y))
       for k=0, n_elements(t)-1 do begin
           for i=0, n_elements(y)-1 do xx[k,*,i] = x
           for i=0, n_elements(x)-1 do zz[k,i,*] = y
       end

       data = sqrt((xx-x0)^2 + (zz-z0)^2)

   ;===========================================
   ; toroidal current density
   ;===========================================
    endif else if(strcmp('jphi', name, /fold_case) eq 1) then begin

        psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)

        data = -grad_shafranov(psi,x,y,tor=itor)

   ;===========================================
   ; vorticity
   ;===========================================
;    endif else if(strcmp('vor', name, /fold_case) eq 1) then begin

;        phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
;                         filename=filename, points=pts, $
;                         rrange=xrange, zrange=yrange)

;        data = grad_shafranov(phi,x,y,tor=itor)

   ;===========================================
   ; divergence
   ;===========================================
   endif else if(strcmp('com', name, /fold_case) eq 1) then begin

       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       data = laplacian(chi,x,y,tor=itor)

   ;===========================================
   ; divergence
   ;===========================================
   endif else if(strcmp('iota', name, /fold_case) eq 1) then begin

       minor_r = read_field('r', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, $
                            rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       r = radius_matrix(x,y,t)
       
       bt = sqrt(I^2/r^2)
       bp = sqrt(s_bracket(psi,psi,x,y)/r^2)

       data = 2.*3.14159*(r * bp) / minor_r * bt

   endif else begin

       t = fltarr(trange[1]-trange[0]+1)
       file_id = h5f_open(filename)

       for i=trange[0], trange[1] do begin
           
           print, 'Time ', i

           time_group_id = h5g_open(file_id, time_name(i))
           mesh = h5_parse(time_group_id, 'mesh', /read_data)   

           field_group_id = h5g_open(time_group_id, 'fields')

                                ; check to see if "name" exists
           nmembers = h5g_get_nmembers(time_group_id, 'fields')
           match = 0
           for m=0, nmembers-1 do begin
               thisname = h5g_get_member_name(time_group_id, 'fields', m)
               if(strcmp(thisname, name, /fold_case) eq 1) then begin
                   name = thisname
                   match = 1
                   break
               endif
           end
           if(match eq 0) then begin
               print, "No field named ", name, " at time slice", i
               continue
           end
       
           field = h5_parse(field_group_id, name, /read_data)
           
           time_id = h5a_open_name(time_group_id, "time")
           t[i-trange[0]] = h5a_read(time_id)
           h5a_close, time_id

           h5g_close, field_group_id
           h5g_close, time_group_id

           data[i-trange[0],*,*]=eval_field(field._data, mesh, points=pts, $
                                          r=x, z=y, op=1, filename=filename, $
                                          xrange=xrange, yrange=yrange)
       end

       h5f_close, file_id

   end

   if(n_elements(h_symmetry) eq 1) then begin
       data = (data + h_symmetry*reverse(data, 2)) / 2.
   endif
   if(n_elements(v_symmetry) eq 1) then begin
       print, "v symmetry = ", v_symmetry
       data = (data + v_symmetry*reverse(data, 3)) / 2.
   endif

   if(keyword_set(diff))then begin
       data = data - $
         read_field(name, x, y, t, slices=time, mesh=mesh, $
                    filename=diff_filename, points=pts, $
                    rrange=xrange, zrange=yrange, $
                    h_symmetry=h_symmetry, v_symmetry=v_symmetry)
   endif

   return, data
end




pro plot_field, name, time, x, y, points=p, filename=filename, mesh=plotmesh, $
                mcolor=mc, lcfs=lcfs, title=title, units=units, $
                maskrange=maskrange, maskfield=maskfield, $
                rrange=rrange, zrange=zrange, $
                xrange=xrange, yrange=yrange, $
                cutx=cutx, cutz=cutz, _EXTRA = ex

   if(n_elements(time) eq 0) then time = 0

   if(size(name, /type) eq 7) then begin
       field = read_field(name, x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=p, $
                          rrange=rrange, zrange=zrange, _EXTRA=ex)
       if(n_elements(field) le 1) then return

       fieldname = translate(name, units=u)

       if(n_elements(units) eq 0) then units=u
       if(n_elements(title) eq 0) then begin
           if(t gt 0) then begin
               title = fieldname + $
                 string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!N!6)!X")', t)
           endif else begin
               title = fieldname + $
                 string(FORMAT='("!6(!8t!6 = ",G0,")!X")', t)
           endelse
       endif
   endif else begin
       field = name
       if(n_elements(field) le 1) then return
   endelse

   if(n_elements(maskrange) eq 2) then begin
       if(strcmp(name, maskfield) eq 1) then begin
           psi = field
       endif else begin
           psi = read_field(maskfield, slices=time, mesh=mesh, $
                            filename=filename, points=p, $
                            rrange=rrange, zrange=zrange)
       endelse
       mask = (psi ge maskrange[0]) and (psi le maskrange[1])
       field = mask*field + (1-mask)*(min(field-mask*field,/absolute))
   endif

   if(n_elements(cutx) gt 0) then begin
       dum = min(x-cutx,i,/absolute)
       plot, y, field[0,i,*]
   endif else if(n_elements(cutz) gt 0) then begin
       dum = min(y-cutz,i,/absolute)
       plot, x, field[0,*,i]
   endif else begin
       contour_and_legend, field, x, y, title=title, label=units, $
         xtitle='!8r!X', ytitle='!8z!X', _EXTRA=ex

       if(keyword_set(lcfs) or n_elements(maskrange) ne 0) then begin
           if(n_elements(psi) eq 0) then begin
               plot_lcfs, time, color=130, points=p, filename=filename
           endif else begin
               plot_lcfs, time, color=130, val=maskrange[0], psi=psi, $
                 x=x, y=y, points=p
               plot_lcfs, time, color=130, val=maskrange[1], psi=psi, $
                 x=x, y=y, points=p
           endelse
       endif
       
       if(keyword_set(plotmesh)) then begin
           loadct, 12
           plot_mesh, mesh, color=color(3,5), /oplot, filename=filename
       endif
   endelse
end

pro plot_field_mpeg, fieldname, mpegame=mpegname, range=range, points=pts, $
                     _EXTRA=ex
    if(n_elements(mpegname) eq 0) then mpegname = fieldname + '.mpeg'
    if(n_elements(pts) eq 0) then pts = 50

    nt = get_parameters("ntime")

    data = read_field(fieldname, x, y, mesh=mesh, points=pts)

    contour_and_legend_mpeg, mpegname, data, x, y, _EXTRA=ex
end

pro compare, file1, file2, time, names=names

   if(n_elements(names) eq 0) then begin
       names = ['psi', 'phi']
   endif

   for i=0, n_elements(names)-1 do begin
       print, "====================="
       print, "Comparing ", names[i]
       print, "---------------------"

       vals1 = read_field(names[i], slices=time, filename=file1)
       vals2 = read_field(names[i], slices=time, filename=file2)

       rms = sqrt((vals1 - vals2)^2)
       d = mean(rms)
       v = variance(rms)

       print, "Mean/variance of RMS error in ", names[i], " = ", $
         d, " / ",  v
   end

   if(i eq 1) then begin
       contour_and_legend, vals1-vals2
   end
end

function energy_mag, filename=filename
   if(n_elements(filename) eq 0) then filename='C1.h5'
   
   nv = read_parameter("numvar", filename=filename)
   scalars = read_scalars(filename=filename)

   E_M = scalars.E_MP._data
   if(nv ge 2) then begin
       E_M = E_M + scalars.E_MT._data
   endif

   return, E_M
end


function energy_kin, filename=filename
   if(n_elements(filename) eq 0) then filename='C1.h5'
   
   nv = read_parameter("numvar", filename=filename)
   scalars = read_scalars(filename=filename)

   E_K = scalars.E_KP._data
   if(nv ge 2) then begin
       E_K = E_K + scalars.E_KT._data 
   endif
   if(nv ge 3) then begin
       E_K = E_K + scalars.E_K3._data
   endif

   return, E_K
end


;==============================================================
; particle_flux
; ~~~~~~~~~~~~~
;
; Returns a time series of the total rate of particles lost 
; through the simulation domain boundary.
; 
; Optional outputs:
;  t: the time array
;  names: the name of each flux component
;  components: an array of time series of each flux component
;==============================================================
function particle_flux, filename=filename, components=comp, names=names, t=t
   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   t = s.time._data

   ipellet = read_parameter('ipellet', filename=filename) 

   names = ['- dN/dt', 'Diffusive', 'Convective', 'Sources']

;    if(ipellet eq 1) then names=[names, 'Pellet']

  
   comp = fltarr(n_elements(names), n_elements(t))
   comp[0,*] = -deriv(t,s.particle_number._data)
   comp[1,*] = s.Particle_Flux_diffusive._data
   comp[2,*] = s.Particle_Flux_convective._data
   comp[3,*] = s.Particle_source._data

;    n = 3

;    if(ipellet eq 1) then begin
;        comp[n,*] = read_parameter('pellet_rate', filename=filename) 
;        n = n+1
;    endif

   return, total(comp, 1)
end


;==============================================================
; energy_dissipated
; ~~~~~~~~~~~~~~~~~
;
; Returns a time series of the rate of energy dissipation in a
; non-conservative system
; 
; Optional outputs:
;  t: the time array
;  names: the name of each energy component
;  components: an array of time series of each energy component
;==============================================================
function energy_dissipated, filename=filename, components=comp, names=names, $
                            t=t

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = read_parameter("numvar", filename=filename)
   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   t = s.time._data
   comp = fltarr(nv*2,n_elements(t))

   names = ['Solenoidal Viscous', 'Solenoidal Ohmic']

   comp[0,*] = s.E_KPD._data + s.E_KPH._data
   comp[1,*] = s.E_MPD._data + s.E_MPH._data
   if(nv ge 2) then begin
       names = [names, ['Toroidal Viscous', 'Toroidal Ohmic']]
       comp[2,*] = s.E_KTD._data + s.E_KTH._data
       comp[3,*] = s.E_MTD._data + s.E_MTH._data
   endif
   if(nv ge 3) then begin
       names = [names, ['Compressional Viscous', 'Thermal Dissipation']]
       comp[4,*] = s.E_K3D._data + s.E_K3H._data
       comp[5,*] = s.E_PD._data + s.E_PH._data
   endif

   comp = -comp

   if(nv ge 3) then begin
       return, replicate(0., n_elements(t))
   endif else begin
       return, total(comp,1)
   endelse
end


;==============================================================
; energy_flux
; ~~~~~~~~~~~
;
; Returns a time series of the total rate of energy lost 
; through the simulation domain boundary.
; 
; Optional outputs:
;  t: the time array
;  names: the name of each flux component
;  components: an array of time series of each flux component
;==============================================================
function energy_flux, filename=filename, components=comp, names=names, t=t
   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   t = s.time._data

   names = ['Pressure', 'Kinetic', 'Poynting', 'Thermal']
  
   comp = fltarr(n_elements(names), n_elements(t))
   comp[0,*] = s.Flux_pressure._data
   comp[1,*] = s.Flux_kinetic._data
   comp[2,*] = s.Flux_poynting._data
;   comp[2,*] = -s.toroidal_current._data*s.loop_voltage._data/(2.*3.14159)
   comp[3,*] = s.Flux_thermal._data

   return, total(comp, 1)
end


;==============================================================
; energy
; ~~~~~~~~~~~
;
; Returns a time series of the total energy within the
; simulation domain
; 
; Optional outputs:
;  t: the time array
;  names: the name of each energy component
;  components: an array of time series of each energy component
;==============================================================
function energy, filename=filename, components=comp, names=names, t=t

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = read_parameter("numvar", filename=filename)
   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   t = s.time._data
   comp = fltarr(nv*2,n_elements(t))

   names = ['Solenoidal KE', 'Solenoidal ME']

   comp[0,*] = s.E_KP._data
   comp[1,*] = s.E_MP._data
   if(nv ge 2) then begin
       names = [names, ['Toroidal KE', 'Toroidal ME']]
       comp[2,*] = s.E_KT._data 
       comp[3,*] = s.E_MT._data
   endif
   if(nv ge 3) then begin
       names = [names, ['Compressional KE', 'Thermal Pressure']]
       comp[4,*] = s.E_K3._data
       comp[5,*] = s.E_P._data
   endif

   return, total(comp, 1)
end


;==============================================================
; energy_error
; ~~~~~~~~~~~~
;
; Returns a time series of the total rate of energy lost 
; through the simulation domain boundary.
; 
; Optional outputs:
;  t: the time array
;  names: the name of each flux component
;  components: an array of time series of each flux component
;==============================================================
function energy_error, filename=filename, components=comp, names=names, t=t
   s = read_scalars(filename=filename)

   if(n_tags(s) eq 0) then return, 0

   t = s.time._data

   names = ['dE/dt', 'Div(S)', 'Dissipated']
   
   if(n_elements(t) le 2) then begin
       print, "Not enough data points."
       return, 0
   endif

   comp = fltarr(n_elements(names), n_elements(t))
   comp[0,*] = deriv(t, energy(filename=filename))
   comp[1,*] = -energy_flux(filename=filename)
   comp[2,*] = energy_dissipated(filename=filename)

   return, total(comp, 1)
end


;============================================
; plot_energy
;============================================
pro plot_energy, name, filename=filename, norm=norm, diff=diff, $
                 overplot=overplot, per_length=per_length, $
                 yrange=yrange, _EXTRA=extra

   if(n_elements(name) eq 0) then begin
       print, "Usage: plot_energy, name"
       print, "where name is one of:"
       print, " 'error':      error in energy conservation"
       print, " 'flux':       power fluxes through domain boundaries"
       print, " 'dissipated': energy lost in reduced models from dissipation"
       print, " 'particle flux': error in particle number conservation"
       print, " 'energy':     components of internal energy"
       return
   endif

   xtitle = '!8t !6(!7s!D!8A!N!6)!3'

   power_units = " !6(!8B!D!60!U2!N!8L!U!63!N/4!7ps!D!8A!N!6)!X"
   energy_units = " !6(!8B!D!60!U2!N!8L!U!63!N/4!7p!6)!X"
   number_units = " !6(!8n!D!60!N/!7s!D!8A!N!6)!X"
   rate_units = " !6(!7s!D!8A!6!U-1!N)!X"

   if(strcmp(name, 'error', /fold_case) eq 1) then begin
       tot = energy_error(filename=filename, comp=comp, names=names, t=t)
       norm_tot = energy(filename=filename)
       title  = '!6Error!X'
       ytitle = '!6Power!X'
       units = power_units
       norm_units = rate_units
   endif else if(strcmp(name, 'flux', /fold_case) eq 1) then begin
       tot = energy_flux(filename=filename, comp=comp, names=names, t=t)
       norm_tot = energy(filename=filename)
       title  = '!6Energy Flux!X'
       ytitle = '!6Power!X'
       units = power_units
       norm_units = rate_units
   endif else if(strcmp(name, 'dissipated', /fold_case) eq 1) then begin
       tot = energy_dissipated(filename=filename, comp=comp, names=names, t=t)
       norm_tot = energy(filename=filename)
       title  = '!6Dissipated Power!X'
       ytitle = '!6Power!X'
       units = power_units
       norm_units = rate_units
   endif else if(strcmp(name, 'particle flux', /fold_case) eq 1) then begin
       tot = particle_flux(filename=filename, comp=comp, names=names, t=t)   
       norm_tot = tot
       title  = '!6Particle Flux!X'
       ytitle = '!6Particle Flux!X'
       units = number_units
       norm_units = ""
   endif else begin
       tot = energy(filename=filename, comp=comp, names=names, t=t)   
       norm_tot = tot
       title  = '!6Internal Energy!X'
       ytitle = '!6Energy!X'
       units = energy_units
       norm_units = ""
   endelse

   if(n_elements(tot) le 2) then return

   n = n_elements(names)

   if(keyword_set(per_length)) then begin
       itor = read_parameter('itor', filename=filename)
       if(itor eq 1) then begin
           xzero = read_parameter('xzero', filename=filename)
           tot = tot / xzero
           comp = comp / xzero
       endif
   endif

   if(keyword_set(diff)) then begin
       title = title + "!6 (Difference from Initial)!X"
       tot = tot - tot[0]
       for i=0, n-1 do begin
           comp[i,*] = comp[i,*] - comp[i,0]
       endfor
   endif

   if(keyword_set(norm)) then begin
       ytitle = "!6Fractional !X" + ytitle
       tot = tot / norm_tot
       for i=0, n-1 do begin
           comp[i,*] = comp[i,*] / norm_tot[*]
       endfor
       ytitle = ytitle + norm_units
   endif else begin
       ytitle = ytitle + units
   endelse

   if(n_elements(yrange) eq 0) then begin
       yrange = [min([tot, min(comp)]), max([tot, max(comp)])]*1.2
   endif

   if(keyword_set(overplot)) then begin
       oplot, t, tot, color=color(0,n+1), _EXTRA=extra
   endif else begin
       plot, t, tot, color=color(0,n+1), xtitle=xtitle, ytitle=ytitle, $
         title=title, yrange=yrange, _EXTRA=extra
   endelse

   for i=0, n-1 do begin
       oplot, t, comp[i,*], color=color(i+1,n+1), _EXTRA=extra
   end

   plot_legend, ['Total', names], color=colors(n+1), _EXTRA=extra
end



; ==============================================================
; nulls
; -----
;
;  Finds field nulls.
;  xpoint = fltarr(2,xpoints): indices of x-point locations
;  axis   = fltarr(2,axes)   : indices of axis locations
; ==============================================================
pro nulls, time, axis=axis, xpoints=xpoint, $
              _EXTRA=extra, psi=psi, r=x, z=z

   if(n_elements(time) eq 0) then time = 0

   if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) $
     then begin
       psi = read_field('psi', x, z, slice=time, _EXTRA=extra)
   endif
   
   field = s_bracket(psi,psi,x,z)
   d2 = dz(dz(psi,z),z)*dx(dx(psi,x),x)
   
   nulls = field lt mean(field)/1e2

   sz = size(field)

   xpoints = 0
   xpoint = 0.
   axes = 0
   axis = 0.

   oldxflux = min(psi)
   oldaflux = min(psi)
   
   for i=0, sz[2]-1 do begin
       for j=0, sz[3]-1 do begin
           if(nulls[0,i,j] eq 0) then continue

           currentmin = field[0,i,j]
           currentpos = [i,j]

           ; find local minimum
           for m=i, sz[2]-1 do begin
               for n=j, sz[3]-1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] lt currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                   endif
                   nulls[0,m,n] = 0
               endfor
               for n=j-1, 0, -1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] lt currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                   endif
                   nulls[0,m,n] = 0
               endfor
           endfor

           ; throw out local minima on boundaries
           if (currentpos[0] eq 0) or (currentpos[0] eq sz[2]-1) then continue
           if (currentpos[1] eq 0) or (currentpos[1] eq sz[3]-1) then continue

           ; determine if point is an x-point or an axis and
           ; append the location index to the appropriate array
           if(d2[0,currentpos[0],currentpos[1]] lt 0) then begin
               if(psi[0,currentpos[0],currentpos[1]] gt oldxflux) then begin
                   xpoint = currentpos
                   oldxflux = psi[0,currentpos[0],currentpos[1]]
               endif
           endif else begin
               if(psi[0,currentpos[0],currentpos[1]] gt oldaflux) then begin
                   axis = currentpos
                   oldaflux = psi[0,currentpos[0],currentpos[1]]
               endif
           endelse
       endfor
   endfor 

   if(n_elements(axis) lt 2) then begin
       dum = max(psi, i)
       axis = [i mod n_elements(x), fix(i / n_elements(x))]
   endif

   if(n_elements(xpoint) ge 2) then begin
       print, 'Found X-point at ', x[xpoint[0]], z[xpoint[1]]
   end
   print, 'Found axis at ', x[axis[0]], z[axis[1]]

end


; ========================================================
; lcfs
; ~~~~
;
; returns the flux value of the last closed flux surface
; ========================================================
function lcfs, time, psi=psi, r=x, z=z, axis=axis, xpoint=xpoint, _EXTRA=extra
   if(n_elements(time) eq 0) then time = 0

   if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) $
     then begin
       psi = read_field('psi', x, z, slice=time, _EXTRA=extra)
   endif

   nulls, time, psi=psi, xpoint=xpoint, axis=axis, r=x, z=z, _EXTRA=extra

   ; flux at magnetic axis
   if(n_elements(axis) lt 2) then begin
       print, "Error: no magnetic axis"
       psi0 = max(psi[0,*,*])
   endif else begin
       psi0 = psi[0,axis[0,0],axis[1,0]]
   endelse
   if(n_elements(axis) gt 2) then begin
       print, "Warning: there is more than one magnetic axis"
   endif
   print, "Flux on axis:", psi0

   ; limiting value
   ; Find limiting flux by calculating outward normal derivative of
   ; the normalized flux.  If this derivative is negative, there is a
   ; limiter.
   sz = size(psi)

   psiz = dz(psi,z)
   psix = dx(psi,x)

   normal_mask = psi*0.
   normal_mask[0,      *,      0] = 1.
   normal_mask[0,      *,sz[3]-1] = 1.
   normal_mask[0,      0,      *] = 1.
   normal_mask[0,sz[2]-1,      *] = 1.

   normal_deriv = psi*0.
   normal_deriv[0,      *,      0] = -3.14159625/2.     ; bottom
   normal_deriv[0,      *,sz[3]-1] =  3.14159625/2.     ; top
   normal_deriv[0,      0,      *] =  3.14159625        ; left
   normal_deriv[0,sz[2]-1,      *] =  0.                ; right
   normal_deriv = $
     (psix*cos(normal_deriv) + psiz*sin(normal_deriv))*normal_mask

   normal_deriv = normal_deriv lt 0

   psi_bound = psi*normal_deriv + (1-normal_deriv)*1e10

   psilim = min(psi_bound-psi0, i, /absolute)
   psilim = psi_bound[i]
   print, "Flux at limiter", psilim

   ; flux at separatrix
   sz = size(xpoint)
   if(sz[0] gt 0) then begin
       xfluxes = fltarr(sz[0])
       for i=0, sz[0]-1 do xfluxes[i] = psi[0,xpoint[0,i],xpoint[1,i]]
       psix = min(xfluxes-psi0, i, /absolute)
       psix = xfluxes[i]
       print, "Flux at separatrix:", psix

       if(abs(psix-psi0) gt abs(psilim-psi0)) then begin
           print, "Plasma is limited."
       endif else begin
           print, "Plasma is diverted."
           psilim = psix
       endelse
   endif
   
   return, psilim
end


; ========================================================
; plot_lcfs
; ~~~~~~~~~
;
; plots the last closed flux surface
; ========================================================
pro plot_lcfs, time, color=color, val=psival, psi=psi, x=x, y=y, points=pts, $
               filename=filename

    if(n_elements(psi) eq 0) then begin
        psi = read_field('psi', x, y, slice=time, points=pts, $
                         filename=filename)
    endif

    ; if psival not passed, choose limiter value
    if(n_elements(psival) eq 0) then begin
        psival = lcfs(time, psi=psi, r=x, z=y, points=pts)
    endif

    ; plot contour
    if(n_elements(color) ne 0) then loadct, 12
    contour, psi, x, y, /overplot, nlevels=1, levels=psival, $
      color=color, thick=2
end


pro plot_timings, filename=filename, overplot=overplot, _EXTRA=extra

   if(n_elements(filename) eq 0) then filename = 'C1.h5'

   if(hdf5_file_test(filename) eq 0) then return

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   timings = h5_parse(root_id, "timings", /read_data)
   h5g_close, root_id
   h5f_close, file_id

   t_solve = timings.t_solve_b._data + timings.t_solve_v._data + $
     timings.t_solve_n._data + timings.t_solve_p._data
   t_output = timings.t_output_cgm._data + timings.t_output_hdf5._data + $
     timings.t_output_reset._data

   loadct, 12

   if(keyword_set(overplot)) then begin
       oplot, timings.t_onestep._data
   endif else begin
       plot, timings.t_onestep._data, title='!6Timings!3', $
         xtitle='!6Time Step!3', ytitle='!8t!6 (s)!3', _EXTRA=extra
   endelse
   oplot, timings.t_ludefall._data, linestyle=2, color=30
   oplot, timings.t_sources._data, linestyle=1, color=60
   oplot, timings.t_aux._data, linestyle=1, color=80
   oplot, timings.t_smoother._data, linestyle=1, color=100
   oplot, t_solve, linestyle=2, color=160
   oplot, t_output, linestyle=2, color=200


   plot_legend, ['Onestep', 'ludefall', 'sources', 'aux', $
                 'smoother', 'solve', 'output'], $
     linestyle=[0,2,1,1,1,2,2], color=[-1,30,60,80,100,160,200]

end


function minor_radius, r=x, z=z, _extra=extra
   psi = read_field('psi',x,z,t,_extra=extra)
   sz = size(psi)

   rr = radius_matrix(x,z,t)
   zz = rr
   a = rr

   for k=0,sz[1]-1 do begin
       for j=0,sz[2]-1 do zz[k,j,*] = z
       psimax = max(psi[k,*,*], i)
       xi = i/sz[3]
       zi = i mod sz[3]
       a[k,*,*] = sqrt((rr[k,*,*]-x[xi])^2 + zz[k,*,*]^2)
   endfor

   return, a
end

; =====================================================
; Scalar functions
; =====================================================


; ==============
; beta = 2*p/B^2
; ==============
function beta, filename=filename
   nv = read_parameter("numvar", filename=filename)
   if(nv lt 3) then begin
       print, "Must be numvar = 3 for beta calculation"
       return, 0
   endif

   gamma = read_parameter('gam', filename=filename)
   s = read_scalars(filename=filename)

   return, (gamma - 1.)*s.E_P._data / (s.E_MP._data + s.E_MT._data)
end


; ==========================
; beta_poloidal = 2*P/(Ip^2)
; ==========================
function beta_poloidal, filename=filename
   nv = read_parameter("numvar", filename=filename)
   if(nv lt 3) then begin
       print, "Must be numvar = 3 for beta calculation"
       return, 0
   endif

   gamma = read_parameter('gam', filename=filename)
   s = read_scalars(filename=filename)

   return, 2.*(gamma-1.)*s.E_P._data/s.toroidal_current._data^2
end


; ====================================
; beta_normal = beta_t * (B_T*a / I_p)
; ====================================
function beta_normal, filename=filename
   nv = read_parameter("numvar", filename=filename)
   if(nv lt 3) then begin
       print, "Must be numvar = 3 for beta calculation"
       return, 0
   endif 

   gamma = read_parameter('gam', filename=filename)
   bzero = read_parameter('bzero', filename=filename)
   xmag = read_parameter('xmag', filename=filename)
   xlim = read_parameter('xlim', filename=filename)
   psi = read_field('psi',slice=0, mesh=mesh)
   s = read_scalars(filename=filename)

   area = mesh.width._data*mesh.height._data

   if(xmag eq 0.) then begin
       B_T = sqrt(2.*s.E_MT._data/area)
       B_T = B_T[0]
       a = 1.
   endif else begin
       B_T = bzero/xmag
       a = abs(xmag - xlim)
   endelse

   beta_t = (gamma-1.)*s.E_P._data/(s.E_MT._data)
   
   print, 'B_T =', B_T
   print, 'a = ', a
   print, 'beta_T = ', beta_t
   print, 'area = ', area

   return, beta_t / (s.toroidal_current._data/(a*B_T))
end


function read_scalar, scalarname, filename=filename, title=title, $
                      symbol=symbol, units=units, time=time

   if(n_elements(scalarname) eq 0) then begin
       print, "Error: no scalar name provided"
       return, 0
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'

   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   time = s.time._data

   units = ''
   va0 = '!8v!DA!60!N!X'
   b0 = '!8B!D!60!N!X'
   n0 = '!8n!D!60!N!X'
   t0 = '!7s!DA!60!N!X'
   l0 = '!8L!X'
   pi4 = '!64!7p!X'
   sq = '!U!62!N!X'
   c = '!8c!X'

   if(strcmp("toroidal current", scalarname, /fold_case) eq 1) or $
     (strcmp("it", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current._data
       title = 'Toroidal Current'
       symbol = '!8I!DT!N!X'
       units = c + b0 + '/' + pi4 + l0
   endif else $
     if (strcmp("toroidal flux", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_flux._data
       title = 'Toroidal Flux'
       symbol = 'Flux'
       units = l0 + '!U2!N' + b0
   endif else $
     if (strcmp("reconnected flux", scalarname, /fold_case) eq 1) then begin
       data = abs(s.reconnected_flux._data)
       title = 'Reconnected Flux'
       symbol = translate('psi', units=u)
       units = u
   endif else $
     if (strcmp("loop voltage", scalarname, /fold_case) eq 1) or $
     (strcmp("vl", scalarname, /fold_case) eq 1) then begin
       data = s.loop_voltage._data
       title = 'Loop Voltage'
       symbol = '!8V!DL!N!X'
   endif else $
     if (strcmp("beta", scalarname, /fold_case) eq 1) then begin
       data = beta(filename=filename)
       title = 'Average Beta'
       symbol = '!7b!X'
   endif else if $
     (strcmp("poloidal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bp", scalarname, /fold_case) eq 1) then begin
       data = beta_poloidal(filename=filename)
       title = 'Average Poloidal Beta'
       symbol = '!7b!D!8p!N!X'
   endif else $
     if (strcmp("normal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bn", scalarname, /fold_case) eq 1) then begin
       data = beta_normal(filename=filename)
       title = 'Average Normal Beta'
       symbol = '!7b!D!8N!N!3'
   endif else if $
     (strcmp("kinetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("ke", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_KP._data 
       if(nv ge 2) then data = data + s.E_KT._data 
       if(nv ge 3) then data = data + s.E_K3._data
       title = 'Kinetic Energy'
       symbol = '!8E!X'
       units = b0 + '!U2!N' + l0 + '!U3!N' + '/' + pi4
   endif else if $
     (strcmp("magnetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("me", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_MP._data 
       if(nv ge 2) then data = data + s.E_MT._data 
       title = 'Magnetic Energy'
       symbol = '!8E!X'
       units = b0 + '!U2!N' + l0 + '!U3!N' + '/' + pi4
   endif else if $
     (strcmp("energy", scalarname, /fold_case) eq 1) then begin
       nv = read_parameter("numvar", filename=filename)
       data = energy(filename=filename)
       title = 'Total Energy'
       symbol = '!8E!X'
       units = b0 + '!U2!N' + l0 + '!U3!N' + '/' + pi4
   endif else if $
     (strcmp("particles", scalarname, /fold_case) eq 1) or $
     (strcmp("n", scalarname, /fold_case) eq 1) then begin
       data = s.particle_number._data
       title = 'Particle Number'
       symbol = '!8N!X'
       units = n0 + l0 + '!U3!N!X'
   endif else if $
     (strcmp("angular momentum", scalarname, /fold_case) eq 1) then begin
       data = s.angular_momentum._data
       title = 'Angular Momentum'
       symbol = '!8L!D!9P!N!X'
       units = b0 + '!U2!N' + l0 + '!U3!N' + '/' + pi4 + t0
   endif else if (strcmp("circulation", scalarname, /fold_case) eq 1) or $
     (strcmp("vorticity", scalarname, /fold_case) eq 1) then begin
       data = s.circulation._data
;      data = s.vorticity._data
       title = 'Circulation'
       symbol = '!9I!S!7x!R!A!6_!D!9P!N !8dA !x'
       units = va0 + l0
   endif else begin
       print, 'Scalar ', scalarname, ' not recognized.'
       return, 0
   endelse
   
   itor = read_parameter('itor', filename=filename)
   if(itor eq 0) then units = units + ' / ' + l0

   return, data
end


pro plot_scalar, scalarname, x, filename=filename, names=names, $
                 _EXTRA=extra, overplot=overplot, $
                 ylog=ylog, xlog=xlog, left=left, $
                 power_spectrum=power_spectrum, per_length=per_length

  if(n_elements(filename) eq 0) then filename='C1.h5'

  if(n_elements(names) eq 0) then names=filename

  nfiles = n_elements(filename)
  if(nfiles gt 1) then begin
      if(n_elements(c) eq 0) then begin
          colors = colors(nfiles)
      endif else begin
          colors = fltarr(n_elements(scalarname))
          colors[*] = c
      endelse

      for i=0, nfiles-1 do begin
          if(n_elements(x) eq 0) then begin
              plot_scalar, scalarname, filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=power_spectrum, per_length=per_length
          endif else begin
              plot_scalar, scalarname, x[i], filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=power_spectrum, per_length=per_length
          endelse
      end
      if(n_elements(names) gt 0) then begin
          plot_legend, names, color=colors, ylog=ylog, xlog=xlog, left=left
      endif    

      return
  endif 

  data = read_scalar(scalarname, filename=filename, time=time, $
                     title=title, symbol=symbol, units=units)
  if(n_elements(data) le 1) then return

  title = '!6' + title + '!X'
  ytitle = symbol
  if(strlen(units) gt 0) then begin
      ytitle = ytitle + '!6 (' + units + ')!X'
  endif

  if(keyword_set(power_spectrum)) then begin
      xtitle = '!7x!6 (!7s!D!8A!N!6!U-1!N)!X'
      data = power_spectrum(data, frequency=tdata, t=max(time)) 
  endif else begin
      xtitle = '!8t!6 (!7s!D!8A!N!6)!X'
      tdata = time
  endelse

  if(keyword_set(per_length)) then begin
      itor = read_parameter('itor', filename=filename)
      if(itor eq 1) then begin
          xzero = read_parameter('xzero', filename=filename)
          data = data / xzero
      endif
  endif
  
  if(n_elements(x) eq 0) then begin
      if(keyword_set(overplot)) then begin
          oplot, tdata, data, color=c, _EXTRA=extra
      endif else begin
          plot, tdata, data, xtitle=xtitle, ytitle=ytitle, $
            title=title, _EXTRA=extra, ylog=ylog, xlog=xlog, $
            color=c
      endelse
  endif else begin
      xi = x
      x = fltarr(1)
      z = fltarr(1)
      x[0] = xi
      z[0] = data[n_elements(data)-1]

      if(keyword_set(overplot)) then begin
          oplot, x, z, color=c, _EXTRA=extra
      endif else begin
          plot, x, z, $
            title=title, xtitle=xtitle, ytitle=ytitle, $
            _EXTRA=extra, ylog=ylog, xlog=xlog, color=c
      endelse
  endelse
end


; ==================================================
; plot_pol_velocity
; ~~~~~~~~~~~~~~~~~
;
; makes a vector plot of the poloidal velocity
; ==================================================
pro plot_pol_velocity, time, filename=filename, points=pts, maxval=maxval, $
                       lcfs=lcfs, _EXTRA=extra

  if(n_elements(pts) eq 0) then pts=25

  nv = read_parameter('numvar', filename=filename)
  itor = read_parameter('itor', filename=filename)

  phi = read_field('phi', x, z, t, filename=filename, $
                   slice=time, points=pts)
  if(n_elements(phi) le 1) then return

  itor = read_parameter('itor', filename=filename)
  if(itor eq 1) then r = radius_matrix(x,z,t) else r = 1.

  vx = -dz(phi,z)/r
  vz =  dx(phi,x)/r

  if(nv ge 3) then begin
      chi = read_field('chi', x, z, t, filename=filename, slice=time, $
                       points=pts)
      vx = vx + dx(chi,x)
      vz = vz + dz(chi,z)
  endif

  bigvel = max(vx^2 + vz^2)
  print, "maximum velocity: ", bigvel
  if(n_elements(maxval) ne 0) then begin
      length = bigvel/maxval
  endif else length=1

  if(n_elements(title) eq 0) then begin
       if(t gt 0) then begin
           title = "!6Poloidal Flow " + $
             string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!N!6)!3")', t)
       endif else begin
           title = "!6Poloidal Flow " + $
             string(FORMAT='("!6(!8t!6 = ",G0,")!3")', t)
       endelse
   endif

  maxstr=string(format='("!6max(!8v!Dpol!N!6) = ",G0)',bigvel)

  velovect, reform(vx), reform(vz), x, z, length=length, _EXTRA=extra, $
    xtitle='!8r !6(!8L!6)!3', ytitle='!8z !6(!8L!6)!3', $
    title=title, subtitle=maxstr

   if(keyword_set(lcfs)) then begin
       plot_lcfs, time, color=130, points=pts>100, filename=filename
   endif
end


; ==================================================
; plot_tor_velocity
; ~~~~~~~~~~~~~~~~~
;
; makes a contour plot of the toroidal velocity
; ==================================================
pro plot_tor_velocity, time, filename=filename, points=pts, $
                       lcfs=lcfs, _EXTRA=extra

  nv = read_parameter('numvar', filename=filename)
  itor = read_parameter('itor', filename=filename)

  if(nv lt 2) then begin
      print, "numvar < 2"
      return
  endif

  v = read_field('V', x, z, t, filename=filename, $
                 slice=time, points=pts)

  if(itor eq 1) then r = radius_matrix(x,z,t) else r = 1.

  vz = v/r

  if(n_elements(title) eq 0) then begin
       if(t gt 0) then begin
           title = "!6Toroidal Flow " + $
             string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!N!6)!3")', t)
       endif else begin
           title = "!6Toroidal Flow " + $
             string(FORMAT='("!6(!8t!6 = ",G0,")!3")', t)
       endelse
   endif

  contour_and_legend, vz, x, z, _EXTRA=extra, title=title

   if(keyword_set(lcfs)) then begin
       plot_lcfs, time, color=130, points=pts, filename=filename
   endif
end
  

;==================================================================
; flux_average_field
; ~~~~~~~~~~~~~~~~~~
;
; Computes the flux averages of a field at a given time.
; This function is only for internal use.  Users should instead
; call flux_average
;==================================================================
function flux_average_field, field, psi, x, z, t, bins=bins, flux=flux, $
                             area=area, range=range, _EXTRA=extra

   sz = size(field)

   points = sqrt(sz[2]*sz[3])

   if(n_elements(bins) eq 0) then bins = fix(points/8)

   result = fltarr(sz[1], bins)
   flux = fltarr(sz[1], bins)
   area = fltarr(sz[1], bins)
   sfield = 0

   psival = lcfs(psi=psi, r=x, z=z, axis=axis, xpoint=xpoint, _EXTRA=extra)

   if(n_elements(range) eq 0) then begin
       ; if range not provided, use all flux within lcfs
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = [psival, max(psi[k,*,*])]
   endif else if(n_elements(range) eq 2) then begin
       oldrange = range
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = oldrange
   endif

   ; remove divertor region from consideration
   div_mask = fltarr(n_elements(x), n_elements(z))
   if(n_elements(xpoint) ge 2) then begin
       if(xpoint[1] lt axis[1]) then begin
           div_mask[*,xpoint[1]:n_elements(z)-1] = 1.
       endif else begin
           div_mask[*,0:xpoint[1]] = 1.
       endelse
   endif else begin
       div_mask = 1.
   endelse

   itor = read_parameter('itor', _EXTRA=extra)
   if(itor eq 1) then r = radius_matrix(x,z,t)

   for k=0, sz[1]-1 do begin
       dpsi = float(range[k,1] - range[k,0])/float(bins)

       for p=0, bins-1 do begin
           fval = dpsi*p + range[k,0]

           mask = float((psi[k,*,*] gt fval) and (psi[k,*,*] le fval+dpsi))
           mask = mask*div_mask

           if(itor eq 1) then mask = mask*r[k,*,*]

           if(n_elements(field) gt 1) then $
             sfield = field[k,*,*]*mask[0,*,*]

           tot = total(mask)

           flux[k,p] = fval + dpsi/2.
           if(tot eq 0.) then begin
               print, 'error: no data points in bin', p, fval, fval+dpsi
               result[k, p] = 0
           endif else begin
               result[k, p] = total(sfield)/total(mask)
           endelse
           area[k,p] = total(mask)*mean(deriv(x))*mean(deriv(z))
       endfor
   endfor

   return, result
end


;==================================================================
; flux_average
; ~~~~~~~~~~~~
;
; Computes the flux averages of a field at a given time
;
; input:
;   field: the field to flux-average.  This may be a string (the
;          name of the field) or the field values themselves.
;   psi:   the flux field
;   x:     r-coordinates
;   z:     z-coordinates
;   t:     t-coordinates;
;   bins:  number of bins into which to subdivide the flux
;   range: range of flux values
;
; output: 
;   flux:  flux value of each bin
;   area:  cross-sectional area of each flux bin
;   name:  the formatted name of the field
; symbol:  the formatted symbol of the field
;  units:  the formatted units of the field
;==================================================================
function flux_average, field, time, psi=psi, x=x, z=z, t=t, $
                       flux=flux, area=area, bins=bins, points=points, $
                       name=name, symbol=symbol, units=units, _EXTRA=extra
  
   type = size(field, /type)

   sz = size(field)

   if(n_elements(points) eq 0) then begin
       if(type ne 7) then points = sz[2] else points=50
   endif

   if((n_elements(psi) le 1) $
      or (n_elements(x) eq 0) or (n_elements(z) eq 0) $
      or (n_elements(t) eq 0)) then begin

       psi = read_field('psi', x, z, t, slice=time, points=points, $
                        _EXTRA=extra)

       sz = size(psi)
       if(n_elements(psi) le 1) then return, 0
   endif


   if(type eq 7) then begin ; named field
       if (strcmp(field, 'Safety Factor', /fold_case) eq 1) or $
         (strcmp(field, 'q', /fold_case) eq 1) then begin
           minor_r = read_field('r', x, z, t, points=points, $
                                slice=time, _EXTRA=extra)
           I = read_field('I', x, z, t, points=points, $
                          slice=time, _EXTRA=extra)
           r = radius_matrix(x,z,t)

           bt = sqrt(I^2/r^2)
           bp = sqrt(s_bracket(psi,psi,x,z)/r^2)

           field = minor_r * bt / (r * bp)
           units = ''
           name = '!8q!X'

       endif else if (strcmp(field, 'gam', /fold_case) eq 1) or $
         (strcmp(field, 'q', /fold_case) eq 1) then begin
           minor_r = read_field('r', x, z, t, points=points, $
                                slice=time, _EXTRA=extra)
           I = read_field('I', x, z, t, points=points, $
                          slice=time, _EXTRA=extra)
           P = read_field('P', x, z, t, points=points, $
                          slice=time, _EXTRA=extra)
           den = read_field('den', x, z, t, points=points, $
                            slice=time, _EXTRA=extra)
           gam = read_parameter('gam', _EXTRA=extra)

           r = radius_matrix(x,z,t)
           bt = sqrt(I^2/r^2)
           bp = sqrt(s_bracket(psi,psi,x,z)/r^2)
           q = minor_r * bt / (r * bp)
           
           field = 2*gam*p/(den*r^2)*(1.+2./q^2)
           name = '!4x!D!6GAM!N!X'


       endif else begin
           name = translate(field, units=units)
           field = read_field(field, x, z, t, slice=time, points=points,$
                              _EXTRA=extra)
       endelse
   endif

   symbol = name

   return, flux_average_field(field, psi, x, z, t, $
     flux=flux, area=area, bins=bins, _EXTRA=extra)
end



;======================================================
; plot_flux_average
; ~~~~~~~~~~~~~~~~~
;
; plots the flux average quantity "name" at a give time
;======================================================
pro plot_flux_average, field, time, filename=filename, points=pts, $
                       _EXTRA=extra, color=c, names=names, $
                       xlog=xlog, ylog=ylog, left=left, overplot=overplot, $
                       lcfs=lcfs, normalized_flux=nflux, $
                       minor_radius=minor_radius, smooth=sm, t=t


   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(n_elements(time) eq 0) then time=0

   nfiles = n_elements(filename)
   if(nfiles gt 1) then begin
       if(n_elements(names) eq 0) then names=filename
       if(n_elements(c) eq 0) then colors = colors(nfiles)

       for i=0, nfiles-1 do begin
           newfield = field
           plot_flux_average, newfield, time, filename=filename[i], $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=nflux, minor_radius=minor_radius, smooth=sm
       end
       if(n_elements(names) gt 0) then begin
           plot_legend, names, color=colors, ylog=ylog, xlog=xlog, left=left
       endif    
       
       return
   endif
   if(n_elements(time) gt 1) then begin
       if(n_elements(names) eq 0) then names=strarr(n_elements(time))
       if(n_elements(c) eq 0) then colors=colors(n_elements(time))
       
       for i=0, n_elements(time)-1 do begin
           newfield = field
           plot_flux_average, newfield, time[i], filename=filename, $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=nflux, minor_radius=minor_radius, smooth=sm, t=t
           names[i] = string(format='(%"!8t!6 = %d !7s!D!8A!N!X")', t)
       end

       if(n_elements(names) gt 0) then begin
           plot_legend, names, color=colors, ylog=ylog, xlog=xlog, left=left
       endif

       return
   endif

   xtitle='!7w!3'
   title = ''

   fa = flux_average(field,time,flux=flux,points=pts,filename=filename,t=t, $
                    name=title, symbol=symbol, units=units)
   ytitle = symbol
   if(strlen(units) gt 0) then ytitle = ytitle + '!6 ('+units+ '!6)!X'
   if(n_elements(fa) le 1) then begin
       print, 'Error in flux_average. returning.'
       return
   endif

;   if(t gt 0) then begin
;       title = "!12<" + title + $
;        string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!N!6)!12>!7!Dw!N!X")',t)
;   endif else begin
;       title = "!12<!8" + title + $
;         string(FORMAT='("!6(!8t!6 = ",G0,")!3!12>!7!Dw!N!X")', t)
;   endelse
   title = "!12<" + title + "!12>!7!Dw!N!X"

   if(keyword_set(nflux) or keyword_set(lcfs)) then begin
       lcfs_psi = lcfs(time, filename=filename, psi=psi)
   endif

   if(keyword_set(nflux)) then begin
       flux = (flux - max(psi)) / (lcfs_psi-max(psi))
       xtitle = '!7W!3'
       lcfs_psi = 1.
   endif

   if(keyword_set(minor_radius)) then begin
       flux = flux_average('r',time,points=pts,file=filename,t=t, $
                    name=xtitle, units=units)
       xtitle = '!12<' + xtitle + '!12> !6 ('+units+')!X'
   endif

   if(n_elements(sm) eq 1) then begin
       fa = smooth(fa,sm)
   end

   if(keyword_set(overplot)) then begin
       oplot, flux, fa, color=c
   endif else begin
       plot, flux, fa, xtitle=xtitle, $
         ytitle=ytitle, title=title, xlog=xlog, ylog=ylog, _EXTRA=extra
   endelse

   if(keyword_set(lcfs)) then begin
       oplot, [lcfs_psi,lcfs_psi], !y.crange, linestyle=1, color=c
   endif
end


pro plot_poloidal_rotation, _EXTRA=extra
   phi = read_field('phi', x, z, t, _EXTRA=extra)
   den = read_field('den', x, z, t, _EXTRA=extra)

   phi = den*phi

   ihp = reverse(phi, 3)

   anti = (phi - ihp) / 2.
   sym  = (phi + ihp) / 2.

;   contour_and_legend, [anti[0,*,*],sym[0,*,*]], x, z

   plot, t, total(total(sym,2),2)
   
end

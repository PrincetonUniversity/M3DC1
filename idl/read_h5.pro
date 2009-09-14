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

function z_matrix, x, z, t
    nx = n_elements(x)
    nz = n_elements(z)
    zz = fltarr(n_elements(t), nx, nz)
    for k=0, n_elements(t)-1 do begin
        for j=0, nx-1 do zz[k,j,*] = z
    end
    return, zz
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
   if(n_elements(filename) ne 1) then begin

   endif else if(file_test(filename) eq 0) then begin

   endif else if(h5f_is_hdf5(filename) eq 0) then begin

   endif else return, 1

   print, "Error: ", filename, " is not a valid file."
   return, 0
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

   if(n_elements(filename) gt 1) then begin
       attr = fltarr(n_elements(filename))
       for i=0, n_elements(filename)-1 do begin
           attr[i] = read_parameter(name,filename=filename[i],print=pr)
       end
       return, attr
   endif

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
   if(t lt 0) then begin
       label = "equilibrium"
   endif else begin
       label = string(FORMAT='("time_",I3.3)', t)
   endelse
   return, label
end

function time_string, t
   if(t gt 0) then begin
       str = string(FORMAT='("!8t!6 = ",G0," !7s!D!8A!N!X")', t)
   endif else begin
       str = string(FORMAT='("!8t!6 = ",G0,"!X")', t)
   endelse
   return, str
end

;========================================================
; get_slice_time
; ~~~~~~~~~~~~~~
;
; Returns the physical time associated with time slice
;========================================================
function get_slice_time, filename=filename, slice=slice
   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(slice) eq 0) then slice=0

   file_id = h5f_open(filename)
   time_group_id = h5g_open(file_id, time_name(slice))
   time_id = h5a_open_name(time_group_id, "time")

   t = h5a_read(time_id)

   h5a_close, time_id
   h5g_close, time_group_id
   h5f_close, file_id

   return, t
end


;=========================================================
; read_mesh
; ~~~~~~~~~
;
; Returns the mesh data structure at a given time slice
;=========================================================
function read_mesh, filename=filename, slice=t
   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(t) eq 0) then t = 0

   if(hdf5_file_test(filename) eq 0) then return, 0

   file_id = h5f_open(filename)

   time_group_id = h5g_open(file_id, time_name(t))
   mesh = h5_parse(time_group_id, 'mesh', /read_data)   

   h5g_close, time_group_id
   h5f_close, file_id

   return, mesh
end


pro get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, _EXTRA=extra
   b0_norm = read_parameter('b0_norm', _EXTRA=extra)
   n0_norm = read_parameter('n0_norm', _EXTRA=extra)
   l0_norm = read_parameter('l0_norm', _EXTRA=extra)
end

;===============================================================
; convert_units
; ~~~~~~~~~~~~~
;
; converts x having dimensions d to cgs units
; where b0, n0, and l0 are the normalizations (in cgs units)
;===============================================================
pro convert_units, x, d, b0, n0, l0, cgs=cgs, mks=mks
   if(n_elements(x) eq 0) then return

   if(not (keyword_set(cgs) or keyword_set(mks))) then return

   if(b0 eq 0 or n0 eq 0 or l0 eq 0) then begin
       print, "Warning: unknown conversion factors."
       return
   endif

   val = 1.
   if(keyword_set(cgs)) then begin
       fp = (4.*!pi)
       c0 = 3.e10
       v0 = 2.18e11*b0/sqrt(n0)
       t0 = l0/v0
       temp0 = b0^2/(fp*n0) * 1./(1.6022e-12)
       j0 = c0*b0*l0/fp
       e0 = b0^2*l0^3/fp

       val = fp^d[0] $
         * c0^d[1] $
         * n0^d[2] $
         * v0^d[3] $
         * b0^d[4] $
         * t0^d[8] $
         * l0^d[9] $
         * temp0^d[5] $
         * j0^d[6] $
         * e0^d[7]
       
   endif else if(keyword_set(mks)) then begin
       convert_units, x, d, b0, n0, l0, /cgs

       val = (1.e-2)^d[1] $
         * (1.e6)^d[2] $
         * (1.e-2)^d[3] $
         * (1.e-4)^d[4] $
         * (1.e-2)^d[9] $
         * (3.e9)^(-d[6]) $
         * (1.e-7)^d[7]
   end

   x = x*val
end

;=====================================================================
; dimensions
; ~~~~~~~~~~
;
; returs a vector with specified dimensions
;=====================================================================
function dimensions, energy=ener, eta=eta, j0=j, $
                     p0=pres, temperature=temp, e0=elec, $
                     l0=len, light=c, fourpi=fourpi, n0=den, $
                     v0=vel, t0=time, b0=mag, mu0=visc
  d = intarr(10)

  fp =    [1,0,0,0,0,0,0,0,0,0]
  c0 =    [0,1,0,0,0,0,0,0,0,0]
  n0 =    [0,0,1,0,0,0,0,0,0,0]
  v0 =    [0,0,0,1,0,0,0,0,0,0]
  b0 =    [0,0,0,0,1,0,0,0,0,0]
  temp0 = [0,0,0,0,0,1,0,0,0,0]
  j0 =    [0,0,0,0,0,0,1,0,0,0]
  e0 =    [0,0,0,0,0,0,0,1,0,0]
  t0 =    [0,0,0,0,0,0,0,0,1,0]
  l0 =    [0,0,0,0,0,0,0,0,0,1]
  p0 = e0 - 3*l0
  
  if(keyword_set(fourpi)) then d = d + fp*fourpi
  if(keyword_set(light))  then d = d + c0*light
  if(keyword_set(den))    then d = d + n0*den
  if(keyword_set(vel))    then d = d + v0*vel
  if(keyword_set(mag))    then d = d + b0*mag
  if(keyword_set(time))   then d = d + t0*time
  if(keyword_set(len))    then d = d + l0*len
  if(keyword_set(temp))   then d = d + temp0*temp
  if(keyword_set(j))      then d = d + j0*j - 2*l0
  if(keyword_set(ener))   then d = d + ener*e0

  if(keyword_set(elec)) then d = d + elec*(b0+v0-c0)
  if(keyword_set(eta))  then d = d +  eta*(2*l0-t0+fp-2*c0)
  if(keyword_set(pres)) then d = d + pres*(p0)
  if(keyword_set(visc)) then d = d + visc*(p0+t0)

  return, d
end


;=====================================================================
; parse_units
; ~~~~~~~~~~~
;
; x is a vector containing the dimensions of
; [4pi, c, n0, vA0, B0, T0, j0, e0, tA0, L0]
; [  0, 1,  2,   3,  4,   5,  6,  7,  8,  9]
; 
; output is a string containing units
;=====================================================================
function parse_units, x, cgs=cgs, mks=mks
   if(keyword_set(cgs)) then begin
       x[9] = x[9] - 3*x[2] + x[3] + x[1]
       x[8] = x[8]          - x[3] - x[1]
       x[0] = 0
       x[1] = 0
       x[2] = 0
       x[3] = 0
       u = ['!64!7p', '!8c', '!6cm', '!8v!DA!60!N', $
            '!6G', '!6eV', '!6statamps', '!6erg', '!6s', '!6cm'] + '!X'
   endif else if(keyword_set(mks)) then begin
       x[9] = x[9] - 3*x[2] + x[3] + x[1]
       x[8] = x[8]          - x[3] - x[1]
       x[0] = 0
       x[1] = 0
       x[2] = 0
       x[3] = 0
       u = ['!64!7p', '!8c', '!6cm', '!8v!DA!60!N', $
            '!6T', '!6eV', '!6A', '!6J', '!6s', '!6m'] + '!X'
   endif else begin
       x[0] = x[0]   - x[5] - x[6] -   x[7]
       x[1] = x[1]          + x[6]
       x[4] = x[4] + 2*x[5] + x[6] + 2*x[7]
       x[9] = x[9]          - x[6] - 3*x[7]
       x[5] = 0
       x[6] = 0
       x[7] = 0


       u = ['!64!7p', '!8c', '!8n!D!60!N', '!8v!DA!60!N', $
            '!8B!D!60!N', '!6temp', '!6curr', $
            '!6energy', '!7s!DA!60!N', '!8L!D!60!N'] + '!X'
   endelse
   units = ''

   nu = n_elements(x)

   if(max(x) gt 0) then pos=1 else pos=0
   if(min(x) lt 0) then neg=1 else neg=0

   is = 0
   sscript = '("!U!6",G0,"!N!X")'
   for i=0, nu-1 do begin
       if(x[i] gt 0) then begin 
           if(is eq 1) then units = units + ' '
           units = units + u[i]
           if(x[i] ne 1) then $
             units = units + string(format=sscript,x[i])
           is = 1
       endif
   end

   if(pos eq 1) then begin
       if(neg eq 1) then begin
           units = units + '!6/!X'
           is = 0
       end
       for i=0, nu-1 do begin
           if(x[i] lt  0) then begin
               if(is eq 1) then units = units + ' '
               units = units + u[i]
               if(x[i] ne -1) then $
                   units = units + string(format=sscript,-x[i])
               is = 1
           endif
       end
   endif else begin
       for i=0, nu-1 do begin
           if(x[i] lt 0) then begin
               if(is eq 1) then units = units + ' '
               units = units + u[i]
               units = units + string(format=sscript, x[i])
               is = 1
           endif
       end
   endelse

   return, units
end

function make_label, s, d, _EXTRA=extra
   if(n_elements(d) eq 0) then d = dimensions(_EXTRA=extra)
   if(max(d,/abs) eq 0.) then return, s
   return, s + ' !6(!X' + parse_units(d, _EXTRA=extra) + '!6)!X'
end

;======================================================================
; make_units
; ~~~~~~~~~~
;
; creates a string given the specified unit dimensions
;======================================================================
function make_units, _EXTRA=extra
   return, parse_units(dimensions(_EXTRA=extra), _EXTRA=extra)
end


;======================================================================
; translate
; ~~~~~~~~~
;
; provides the associated symbol and units for fields defined in C1.h5
;======================================================================
function translate, name, units=units, itor=itor
   units = dimensions()

   if(strcmp(name, 'psi', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=1+itor)
       return, "!7w!X"
   endif else if(strcmp(name, 'I', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=itor)
       return, "!8I!X"
   endif else if(strcmp(name, 'phi', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=1+itor)
       return, "!8U!X"
   endif else if(strcmp(name, 'V', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=itor)
       return, "!8V!X"
   endif else if(strcmp(name, 'chi', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=1)
       return, "!7v!X"
   endif else if(strcmp(name, 'eta', /fold_case) eq 1) then begin
       units = dimensions(/eta)
       return, "!7g!X"
   endif else if(strcmp(name, 'den', /fold_case) eq 1) then begin
       units = dimensions(/n0)
       return, "!8n!X"
   endif else if(strcmp(name, 'p', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!8p!X"
   endif else if(strcmp(name, 'pe', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!8p!De!N!X"
   endif else if(strcmp(name, 'sigma', /fold_case) eq 1) then begin
       units = dimensions(/n0,t0=-1)
       return, "!7r!X"
   endif else if(strcmp(name, 'kappa', /fold_case) eq 1) then begin
       units = dimensions(/n0, l0=2, t0=-1)
       return, "!7j!X"
   endif else if((strcmp(name, 'visc', /fold_case) eq 1) or $
     (strcmp(name, 'visc_c', /fold_case) eq 1)) then begin
       return, "!7l!X"
   endif else if(strcmp(name, 'jphi', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=itor-1)
       return, "!7D!6!U*!N!7w!X"
   endif else if(strcmp(name, 'vor', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=itor-1)
       return, "!7D!6!U*!N!8U!X"
   endif else if(strcmp(name, 'com', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=-1)
       return, "!9G.!17v!X"
   endif  

   return, '!8' + name + '!X'
end



;=========================================================
; plot_mesh
; ~~~~~~~~~
;
; Plots the mesh.
; /boundary: only plot the mesh boundary
; /oplot: plots mesh as overlay to previous plot
;=========================================================
pro plot_mesh, mesh=mesh, oplot=oplot, boundary=boundary, _EXTRA=ex

   if(n_elements(mesh) eq 0) then mesh = read_mesh(_EXTRA=ex)
   if(n_tags(mesh) eq 0) then return

   nelms = mesh.nelms._data
   
   if(not keyword_set(oplot)) then begin
       xtitle = make_label('!8R!X',/l0)
       ytitle = make_label('!8Z!X',/l0)
       plot, mesh.elements._data[4,*], xtitle=xtitle, ytitle=ytitle, $
         mesh.elements._data[5,*], psym = 3, _EXTRA=ex, /nodata
   endif  

   loadct, 12
   col = color(3,5)
 
   version = read_parameter('version', _EXTRA=ex)
   if(version eq 0) then begin
       xzero = read_parameter("xzero", _EXTRA=ex)
       zzero = read_parameter("zzero", _EXTRA=ex)
   endif else begin
       xzero = 0.
       zzero = 0.
   endelse

   if(keyword_set(boundary) and version ge 3) then begin
       boundary = 1
   endif else boundary = 0

   maxr = [max(mesh.elements._data[4,*]), max(mesh.elements._data[5,*])]
   minr = [min(mesh.elements._data[4,*]), min(mesh.elements._data[5,*])]

   for i=long(0), nelms-1 do begin
       a = mesh.elements._data[0,i]
       b = mesh.elements._data[1,i]
       c = mesh.elements._data[2,i]
       t = mesh.elements._data[3,i]
       x = mesh.elements._data[4,i]
       y = mesh.elements._data[5,i]
       bound = fix(mesh.elements._data[6,i])

       p1 = [x, y]
       p2 = p1 + [(b+a) * cos(t), (b+a) * sin(t)]
       p3 = p1 + [b * cos(t) - c * sin(t), $
                  b * sin(t) + c * cos(t)]
      
       if(boundary) then pp=bound else pp=7
 
       if((pp and 1) eq 1) then begin
           if((bound and 1) eq 1) then th=!p.thick*3 else th=!p.thick
           oplot, [p1[0],p2[0]]+xzero, [p1[1],p2[1]]+zzero, color=col, thick=th
       end
       if((pp and 2) eq 2) then begin
           if((bound and 2) eq 2) then th=!p.thick*3 else th=!p.thick
           oplot, [p2[0],p3[0]]+xzero, [p2[1],p3[1]]+zzero, color=col, thick=th
       end
       if((pp and 4) eq 4) then begin
           if((bound and 4) eq 4) then th=!p.thick*3 else th=!p.thick
           oplot, [p3[0],p1[0]]+xzero, [p3[1],p1[1]]+zzero, color=col, thick=th
       end

       maxr[0] = max([maxr[0], p1[0], p2[0], p3[0]])
       minr[0] = min([minr[0], p1[0], p2[0], p3[0]])
       maxr[1] = max([maxr[1], p1[1], p2[1], p3[1]])
       minr[1] = min([minr[1], p1[1], p2[1], p3[1]])
   end

   print, 'Elements: ', nelms
   print, 'sqrt(nodes) (estimated): ', sqrt(nelms/2.)
   print, 'Width: ', maxr[0] - minr[0]
   print, 'Height: ', maxr[1] - minr[1]

end


function find_next_boundary_point, xy, mesh=mesh
   tol = 1.e-5
   nelms = mesh.nelms._data

   for i=long(0), nelms-1 do begin
       a = mesh.elements._data[0,i]
       b = mesh.elements._data[1,i]
       c = mesh.elements._data[2,i]
       t = mesh.elements._data[3,i]
       x = mesh.elements._data[4,i]
       y = mesh.elements._data[5,i]
       bound = fix(mesh.elements._data[6,i])

       p1 = [x, y]
       p2 = p1 + [(b+a) * cos(t), (b+a) * sin(t)]
       p3 = p1 + [b * cos(t) - c * sin(t), $
                  b * sin(t) + c * cos(t)]
       
       if((bound and 1) eq 1) then begin
           if(n_elements(xy) eq 0) then return, p1
           if(abs(p1[0] - xy[0]) lt tol and $
              abs(p1[1] - xy[1]) lt tol) then begin
               return, p2
           end
       endif
       if((bound and 2) eq 2) then begin
           if(n_elements(xy) eq 0) then return, p2
           if(abs(p2[0] - xy[0]) lt tol and $
              abs(p2[1] - xy[1]) lt tol) then begin
               return, p3
           end
       endif
       if((bound and 4) eq 4) then begin
           if(n_elements(xy) eq 0) then return, p3
           if(abs(p3[0] - xy[0]) lt tol and $
              abs(p3[1] - xy[1]) lt tol) then begin
               return, p1
           end
       endif
   end

   print, 'Error: no next boundary point found'
   return, xy
end

function get_boundary_path, mesh=mesh, _EXTRA=extra

   if(n_elements(mesh) eq 0) then mesh = read_mesh(_EXTRA=ex)
   if(n_tags(mesh) eq 0) then return, [0,0]

   nelms = mesh.nelms._data

   nbound = 0
   for i=long(0), nelms-1 do begin
       bound = fix(mesh.elements._data[6,i])
       if((bound and 1) eq 1) then nbound = nbound + 1
       if((bound and 2) eq 2) then nbound = nbound + 1
       if((bound and 4) eq 4) then nbound = nbound + 1
   end

   xy = fltarr(2,nbound)

   xy[*,0] = find_next_boundary_point(mesh=mesh)

   for i=1, nbound-1 do begin
       xy[*,i] = find_next_boundary_point(xy[*,i-1],mesh=mesh)
   end
   
   return, xy
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

   small = (a+b+c)*1e-2

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
function eval, field, localpos, theta, elm, operation=op

   mi = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0]
   ni = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5]
   sum = 0.

   if(n_elements(op) eq 0) then op = 1

   co = cos(theta)
   sn = sin(theta)

   for p=0, 19 do begin
        case op of
        1: sum = sum + field[p,elm]*localpos[0]^mi[p]*localpos[1]^ni[p]
        2: sum = sum $
          + co*field[p,elm]*mi[p]*localpos[0]^(mi[p]-1>0)*localpos[1]^ni[p] $
          - sn*field[p,elm]*ni[p]*localpos[0]^mi[p]*localpos[1]^(ni[p]-1>0)
        3: sum = sum $
          + sn*field[p,elm]*mi[p]*localpos[0]^(mi[p]-1>0)*localpos[1]^ni[p] $
          + co*field[p,elm]*ni[p]*localpos[0]^mi[p]*localpos[1]^(ni[p]-1>0)
        4:
        5:
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


function eval_global, field, mesh, pos, elm=i, operation=op

   nelms = mesh.nelms._data

   if(n_elements(i) ne 0) then begin
       if(i ge 0) then begin
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
           
           localpos = [(pos[0]-x)*co + (pos[1]-y)*sn - b, $
                       -(pos[0]-x)*sn + (pos[1]-y)*co]
           
           if(is_in_tri(localpos,a,b,c) eq 1) then begin
               return, eval(field, localpos, t, i, operation=op)
           end
       endif
   endif

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

       localpos = [(pos[0]-x)*co + (pos[1]-y)*sn - b, $
                  -(pos[0]-x)*sn + (pos[1]-y)*co]

       if(is_in_tri(localpos,a,b,c) eq 1) then begin
           return, eval(field, localpos, t, i, operation=op)
       end
   end

   print, 'Error: global position ', pos, ' not found'
   i = -1
   return, 0
end



function eval_field, field, mesh, r=xi, z=yi, points=p, operation=op, $
                     filename=filename, xrange=xrange, yrange=yrange, $
                     mask=mask

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nelms = mesh.nelms._data 

   if(n_elements(p) eq 0) then p = 100
   if(n_elements(op) eq 0) then op = 1
   
   minpos = fltarr(2)
   maxpos = fltarr(2)
   pos = fltarr(2)
   localpos = fltarr(2)
   index = intarr(2)

   version = read_parameter('version', filename=filename)
   if(version eq 0) then begin
       xzero = read_parameter("xzero", filename=filename)
       zzero = read_parameter("zzero", filename=filename)

       nonrect = read_parameter('nonrect', filename=filename)
       if(nonrect eq 0.) then begin
           xmin = min(mesh.elements._data[4,*])
           ymin = min(mesh.elements._data[5,*])
       endif else begin
           xmin = 0.
           ymin = 0.
       endelse

   endif else begin
       xzero = 0.
       zzero = 0.
       xmin = 0.
       ymin = 0.
   endelse

   ; find minimum and maximum node coordinates
   if(n_elements(xrange) lt 2 or n_elements(yrange) lt 2) then begin
       minx = min(mesh.elements._data[4,*])
       maxx = max(mesh.elements._data[4,*])
       miny = min(mesh.elements._data[5,*])
       maxy = max(mesh.elements._data[5,*])
       
       for i=long(0),nelms-1 do begin
           a = mesh.elements._data[0,i]
           b = mesh.elements._data[1,i]
           c = mesh.elements._data[2,i]
           t = mesh.elements._data[3,i]
           x = mesh.elements._data[4,i]
           y = mesh.elements._data[5,i]
           co = cos(t)
           sn = sin(t)

           p1 = [x,y]
           p2 = p1 + [(b+a)*co, (b+a)*sn]
           p3 = p1 + [b*co - c*sn, b*sn + c*co]

           minpos = [min([p2[0], p3[0]]), min([p2[1], p3[1]])]
           maxpos = [max([p2[0], p3[0]]), max([p2[1], p3[1]])]

            if(minpos[0] lt minx) then minx = minpos[0]
            if(maxpos[0] gt maxx) then maxx = maxpos[0]
            if(minpos[1] lt miny) then miny = minpos[1]
            if(maxpos[1] gt maxy) then maxy = maxpos[1]
       endfor

       if(n_elements(xrange) lt 2) then $
         xrange = [minx, maxx] + xzero - xmin
       if(n_elements(yrange) lt 2) then $
         yrange = [miny, maxy] + zzero - ymin
   endif
   
   dx = (xrange[1] - xrange[0]) / (p - 1.)
   dy = (yrange[1] - yrange[0]) / (p - 1.)

   result = fltarr(p,p)
   mask = fltarr(p,p)
   mask[*] = 1.

   small = 1e-3

   ; for each triangle, evaluate points within triangle which fall on
   ; rectilinear output grid
   for i=long(0),nelms-1 do begin
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
                             eval(field, localpos, t, i, op=op)
                           mask[index[0], index[1]] = 0
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

;  determine 'edge' value
   if(max(mask) eq 1) then begin
       num_edge_vals = 0
       edge_val = 0
       for i=0,p-1 do begin
           for j=0,p-1 do begin
               if(mask[i,j] eq 0) then begin 
                   is_edge = 0
                   if(i gt 0) then begin
                       if(mask[i,j] ne mask[i-1,j]) then is_edge = 1
                   endif
                   if(i lt p-1) then begin
                       if(mask[i,j] ne mask[i+1,j]) then is_edge = 1
                   endif
                   if(j gt 0) then begin
                       if(mask[i,j] ne mask[i,j-1]) then is_edge = 1
                   endif
                   if(j lt p-1) then begin
                       if(mask[i,j] ne mask[i,j+1]) then is_edge = 1
                   endif
                   if(is_edge eq 1) then begin
                       edge_val = (edge_val*num_edge_vals $
                                   + result[i,j]) / (num_edge_vals + 1.)
                       num_edge_vals = num_edge_vals + 1.
                   endif
               endif
           endfor
       endfor

;       print, 'Setting value outside boundary = ', edge_val
       print, 'applying mask with edge_val = ', edge_val
       result = result + mask*edge_val
   endif

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


function read_field, name, x, y, t, slices=slices, mesh=mesh, $
                     filename=filename, points=pts, mask=mask, $
                     rrange=xrange, zrange=yrange,equilibrium=equilibrium, $
                     h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                     diff=diff, at_points=at_points,operation=op, $
                     linear=linear, last=last, average=average, $
                     dpsi=dpsi, symbol=symbol, units=units, cgs=cgs, mks=mks

   if(n_elements(slices) ne 0) then time=slices

   if(keyword_set(average)) then begin
       if(n_elements(filename) gt 1) then begin
           n = n_elements(filename)
           if(n_elements(time) eq 1) then time=replicate(time,n)
       endif else if(n_elements(time) gt 1) then begin
           n = n_elements(time)
           if(n_elements(filename) eq 1) then filename=replicate(filename,n)
       endif else n = 1

       data = 0
       for i=0, n-1 do begin
           data = data + $
             read_field(name, x, y, t, slices=time[i], mesh=mesh, $
                        filename=filename[i], points=pts, $
                        rrange=xrange, zrange=yrange, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        diff=diff, at_points=at_points,operation=op, $
                        linear=linear, last=last,symbol=symbol,units=units, $
                       cgs=cgs, mks=mks)
       end
       data = data/n
       return, data
   end
   if(keyword_set(diff)) then begin
       if(n_elements(filename) gt 1) then begin
           n = n_elements(filename)
           if(n_elements(time) eq 1) then time=replicate(time,n)
       endif else if(n_elements(time) gt 1) then begin
           n = n_elements(time)
           if(n_elements(filename) eq 1) then filename=replicate(filename,n)
       endif else n = 1

       data = 0
       for i=0, n-1 do begin
           data = data + $
             read_field(name, x, y, t, slices=time[i], mesh=mesh, $
                        filename=filename[i], points=pts, $
                        rrange=xrange, zrange=yrange, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        at_points=at_points,operation=op, $
                        linear=linear, last=last,symbol=symbol,units=units, $
                       cgs=cgs, mks=mks) $
             *((-1)^i)
       end

       symbol = '!7D!X' + symbol

       return, data
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(filename) gt 1) then filename=filename[0]
   if(n_elements(pts) eq 0) then pts = 50

   if(hdf5_file_test(filename) eq 0) then return, 0

   print, 'reading field ', name, ' from file ', filename, $
     ' linear=', keyword_set(linear), $
     ' points=', pts

   nt = read_parameter("ntime", filename=filename)
   nv = read_parameter("numvar", filename=filename)
   itor = read_parameter("itor", filename=filename)
   version = read_parameter('version', filename=filename)
   ivform = read_parameter('ivform', filename=filename)
   if(version eq 0) then begin
       xzero = read_parameter("xzero", filename=filename)
       zzero = read_parameter("zzero", filename=filename)
   endif else begin
       xzero = 0.
       zzero = 0.
   endelse
   ilin = read_parameter('linear', filename=filename)

   if(keyword_set(last)) then time = [nt-1,nt-1]
   if(ilin eq 1 and keyword_set(equilibrium)) then time=[-1,-1]

   if(n_elements(time) eq 0) then begin
       trange = [0,nt-1]
   endif else if(n_elements(time) eq 1) then begin
       trange = [time, time]
   endif else begin
       trange = time
   endelse

   if((trange[0] ge nt) or (trange[1] ge nt)) then begin
       print, "Error: there are only ", nt-1, " time slices."
       return, 0
   endif

   if(n_elements(at_points) eq 0) then begin
       data = fltarr(trange[1]-trange[0]+1, pts, pts)
       if(ilin eq 1) then base = fltarr(pts,pts)
   endif else begin
       sz = size(at_points, /dimension)
       data = fltarr(trange[1]-trange[0]+1,sz[1])
       if(ilin eq 1) then base = fltarr(sz[1])
   endelse

   d = dimensions()
   symbol=name
 
   ;==========================================
   ; local_beta = 2*P/B^2
   ;==========================================
   if(strcmp('beta', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(n_elements(psi) le 1) then return, 0

       I = read_field('I',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)     
       P = read_field('P',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       b2 = (s_bracket(psi,psi,x,y) + i^2)/r^2

       data = 2.*P/b2
       symbol = '!7b!X'


   ;===========================================
   ; toroidal field
   ;===========================================
   endif else if(strcmp('toroidal field', name, /fold_case) eq 1) or $
     (strcmp('bz', name, /fold_case) eq 1) then begin
       
       I = read_field('I',x,y,t,slices=time, mesh=mesh, filename=filename, $
                      points=pts, rrange=xrange, zrange=yrange, $
                      at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.
   
       data = I/r
       symbol = '!8B!D!9P!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; toroidal velocity
   ;===========================================
   endif else if(strcmp('toroidal velocity', name, /fold_case) eq 1) or $
     (strcmp('vz', name, /fold_case) eq 1) then begin
       
       v = read_field('V',x,y,t,slices=time, mesh=mesh, filename=filename, $
                        points=pts,rrange=xrange,zrange=yrange, $
                      at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.
   
       if(ivform eq 0) then begin
           data = v/r
       endif else if(ivform eq 1) then begin
           data = v*r
       endif
       symbol = '!8u!D!9P!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; thermal velocity
   ;===========================================
   endif else if(strcmp('thermal velocity', name, /fold_case) eq 1) or $
     (strcmp('vt', name, /fold_case) eq 1) then begin

       idens = read_parameter('idens', filename=filename)

       Temp = read_field('T',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
         
       data = sqrt(Temp)
       symbol = '!8v!Dt!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; sound speed
   ;===========================================
   endif else if(strcmp('sound speed', name, /fold_case) eq 1) or $
     (strcmp('cs', name, /fold_case) eq 1) then begin

       vt = read_field('vt',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, at_points=at_points)
       gam = read_parameter('gam', filename=filename)
  
       data = sqrt(gam)*vt
       symbol = '!8c!Ds!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; Mach number
   ;===========================================
   endif else if(strcmp('mach', name, /fold_case) eq 1) or $
     (strcmp('m', name, /fold_case) eq 1) then begin

       cs = read_field('cs',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, at_points=at_points)

       phi = read_field('phi',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, at_points=at_points)
       V = read_field('V',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       v2 = s_bracket(phi,phi,x,y)/r^2 $
         + v^2/r^2 + s_bracket(chi,chi,x,y) $
         + 2.*a_bracket(chi,phi,x,y)/r
  
       data = sqrt(v2)/cs
       symbol = '!8M!X'

   ;===========================================
   ; temperature
   ;===========================================
   endif else if(strcmp('temperature', name, /fold_case) eq 1) or $
     (strcmp('t', name, /fold_case) eq 1) then begin

       P = read_field('P', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
  
       data = p/n
       symbol = '!8T!X'
       d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; electron temperature
   ;===========================================
   endif else if(strcmp('electron temperature', name, /fold_case) eq 1) or $
     (strcmp('te', name, /fold_case) eq 1) then begin

       Pe = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
  
       data = pe/n
       symbol = '!8T!De!N!X'
       d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; ion temperature
   ;===========================================
   endif else if(strcmp('ion temperature', name, /fold_case) eq 1) or $
     (strcmp('ti', name, /fold_case) eq 1) then begin

       P = read_field('P', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       Pe = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
  
       data = (p-pe)/n
       symbol = '!8T!Di!N!X'
       d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; ion pressure
   ;===========================================
   endif else if(strcmp('ion pressure', name, /fold_case) eq 1) or $
     (strcmp('pi', name, /fold_case) eq 1) then begin

       P = read_field('P', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       Pe = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
  
       data = p-pe
       symbol = '!8p!Di!N!X'
       d = dimensions(/p0, _EXTRA=extra)

   ;===========================================
   ; angular momentum
   ;===========================================
   endif else if(strcmp('angular momentum', name, /fold_case) eq 1) or $
     (strcmp('lz', name, /fold_case) eq 1) then begin

       V = read_field('V', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
  
       data = n*v
       symbol = '!8L!D!9P!N!X'
       d = dimensions(/n0, /v0, /l0, _EXTRA=extra)

   ;===========================================
   ; toroidal current
   ;===========================================
   endif else if(strcmp('jz', name, /fold_case) eq 1) then begin

;       jphi = read_field('jphi', x, y, t, slices=time, mesh=mesh, $
;                         filename=filename, points=pts, mask=mask, $
;                         rrange=xrange, zrange=yrange, at_points=at_points)

       lp = read_field('psi', x, y, t, slices=time, mesh=mesh, op=7, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, at_points=at_points)
       psir = read_field('psi', x, y, t, slices=time, mesh=mesh, op=2, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, at_points=at_points)


       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       data = -(lp - psir/r)/r
       symbol = '!8J!D!9P!N!X'
       d = dimensions(/j0, _EXTRA=extra)


   ;===========================================
   ; minor radius
   ;===========================================
   endif else if(strcmp('minor radius', name, /fold_case) eq 1) or $
     (strcmp('r', name) eq 1) then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       nulls, psi, x, y, xpoint=xpoint, axis=axis, $
         filename=filename, _EXTRA=extra
      
       x0 = axis[0]
       z0 = axis[1]
       if(n_elements(at_points) eq 0) then begin
           xx = fltarr(n_elements(t),n_elements(x),n_elements(y))
           zz = fltarr(n_elements(t),n_elements(x),n_elements(y))
           for k=0, n_elements(t)-1 do begin
               for i=0, n_elements(y)-1 do xx[k,*,i] = x
               for i=0, n_elements(x)-1 do zz[k,i,*] = y
           end
           data = sqrt((xx-x0)^2 + (zz-z0)^2)
       endif else begin
           data = sqrt((at_points[0,*]-x0)^2 + (at_points[1,*]-z0)^2)
       endelse

       symbol = '!8r!X'
       d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; polodal angle
   ;===========================================
   endif else if(strcmp('polodal angle', name, /fold_case) eq 1) or $
     (strcmp('theta', name) eq 1) then begin


       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       nulls, psi, x, y, xpoint=xpoint, axis=axis, $
         filename=filename, _EXTRA=extra
      
       x0 = axis[0]
       z0 = axis[1]

       if(n_elements(at_points) eq 0) then begin
           xx = fltarr(n_elements(t),n_elements(x),n_elements(y))
           zz = fltarr(n_elements(t),n_elements(x),n_elements(y))
           for k=0, n_elements(t)-1 do begin
               for i=0, n_elements(y)-1 do xx[k,*,i] = x
               for i=0, n_elements(x)-1 do zz[k,i,*] = y
           end
           data = atan(zz-z0,xx-x0)
       endif else begin
           data = atan(at_points[0,*]-x0,at_points[1,*]-z0)
       endelse

       symbol = '!7h!X'

   ;===========================================
   ; Field strength
   ;===========================================
   endif else if( (strcmp('field strength', name, /fold_case) eq 1) $
                  or (strcmp('b', name, /fold_case) eq 1)) $
     then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       data = sqrt((s_bracket(psi,psi,x,y) + I^2)/r^2)
       symbol = '!3|!5B!3|!X'
       d = dimensions(/b0, _EXTRA=extra)
            
   ;===========================================
   ; Poloidal Field strength
   ;===========================================
   endif else if( (strcmp('poloidal field strength', name, /fold_case) eq 1) $
                  or (strcmp('bp', name, /fold_case) eq 1)) $
     then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       data = sqrt(s_bracket(psi,psi,x,y)/r^2)
       symbol = '!3|!5B!D!8p!N!3|!X'
       d = dimensions(/b0, _EXTRA=extra)


   ;===========================================
   ; Toroidal Field strength
   ;===========================================
   endif else if( (strcmp('toroidal field strength', name, /fold_case) eq 1) $
                  or (strcmp('bt', name, /fold_case) eq 1)) $
     then begin

       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       data = sqrt(I^2/r^2)
       symbol = '!3|!8B!D!9P!N!3|!X'
       d = dimensions(/b0, _EXTRA=extra)


   ;===========================================
   ; Kinetic energy density
   ;===========================================
   endif else if( (strcmp('ke', name, /fold_case) eq 1)) then begin

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       v = read_field('v', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       if(ivform eq 0) then begin
           data = 0.5*n*(s_bracket(u,u,x,y)/r^2+v^2/r^2 + $
                         s_bracket(chi,chi,x,y)+2.*a_bracket(chi,u,x,y)/r)
       endif else if(ivform eq 1) then begin
           data = 0.5*n*(r^2*s_bracket(u,u,x,y)/r^2+r^2*v^2 + $
                         s_bracket(chi,chi,x,y)/r^4+2.*a_bracket(chi,u,x,y)/r)
       endif
       symbol ='!6Kinetic Energy Density!X'
       d = dimensions(/p0, _EXTRA=extra)

   ;===========================================
   ; Alfven velocity
   ;===========================================
   endif else if( (strcmp('va', name, /fold_case) eq 1)) $
     then begin

       b = read_field('b', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, at_points=at_points)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, at_points=at_points)
       data = b/sqrt(den)
       symbol = '!8v!DA!N'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; (minor) radial current density
   ;===========================================
   endif else if(strcmp('jr', name, /fold_case) eq 1) then begin
       
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.
       
       data = a_bracket(i,psi,x,y)/(r*sqrt(s_bracket(psi,psi,x,y)))
       symbol = '!8J!Dr!N!X'
       d = dimensions(/j0, _EXTRA=extra)

   ;===========================================
   ; poloidal current density
   ;===========================================
   endif else if(strcmp('jp', name, /fold_case) eq 1) then begin
       
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.
       
       data = s_bracket(i,psi,x,y)/(r*sqrt(s_bracket(psi,psi,x,y)))
       symbol = '!8J!Dp!N!X'
       d = dimensions(/j0, _EXTRA=extra)

       
   ;===========================================
   ; toroidal current density
   ;===========================================
;      endif else if(strcmp('jphi', name, /fold_case) eq 1) then begin

;          psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
;                           filename=filename, points=pts, $
;                           rrange=xrange, zrange=yrange)

;          data = -grad_shafranov(psi,x,y,tor=itor)
;          symbol = translate('jphi', units=d, itor=itor)

   ;===========================================
   ; vorticity
   ;===========================================
      endif else if(strcmp('vor', name, /fold_case) eq 1) then begin

          phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, at_points=at_points)

          data = grad_shafranov(phi,x,y,tor=itor)
          symbol = translate('vor', units=d, itor=itor)

   ;===========================================
   ; divergence
   ;===========================================
     endif else if(strcmp('com', name, /fold_case) eq 1) then begin

         chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, at_points=at_points)

         data = laplacian(chi,x,y,tor=itor)
         symbol = translate('com', units=d, itor=itor)

   ;===========================================
   ; rotational transform
   ;===========================================
   endif else if(strcmp('iota', name, /fold_case) eq 1) then begin

       minor_r = read_field('r', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, $
                            rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       bt = sqrt(I^2/r^2)
       bp = sqrt(s_bracket(psi,psi,x,y)/r^2)

       data = 2.*!pi*(r * bp) / minor_r * bt
       symbol = '!8i!X'

   ;===========================================
   ; angular velocity
   ;===========================================
   endif else if(strcmp('omega', name, /fold_case) eq 1) then begin

       v = read_field('v', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       if(ivform eq 0) then begin
           data = v/r^2
       endif else if(ivform eq 1) then begin
           data = v
       endif
       symbol = '!7x!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; parallel flow
   ;===========================================
   endif else if(strcmp('vpar', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       v = read_field('v', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, at_points=at_points)
         
       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.


       b2 = (s_bracket(psi,psi,x,y) + I^2)/r^2
       data = ((s_bracket(phi,psi,x,y) + v*I)/r^2 $
               + a_bracket(chi,psi,x,y)/r)/sqrt(b2)
       symbol = '!8u!D!9#!N!X'
       d = dimensions(/v0, _EXTRA=extra)
         
   ;===========================================
   ; radial flow
   ;===========================================
   endif else if(strcmp('vr', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, /equilibrium, slices=time, mesh=mesh, $
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       
       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       if(ivform eq 0) then begin
           data = -(a_bracket(psi,phi,x,y)/r + s_bracket(psi,chi,x,y)) / $
             sqrt(s_bracket(psi,psi,x,y))
       endif else if (ivform eq 1) then begin
           data = -(a_bracket(psi,phi,x,y)*r + s_bracket(psi,chi,x,y)/r^2) / $
             sqrt(s_bracket(psi,psi,x,y))
       endif
       symbol = '!8u!Dr!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   endif else if(strcmp('vx', name, /fold_case) eq 1) then begin

       chi_r = read_field('chi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=2, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       phi_z = read_field('psi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=3, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       
       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       if(ivform eq 0) then begin
           data = -phi_z/r + chi_r
       endif else if (ivform eq 1) then begin
           data = -r*phi_z + chi_r/r^2
       endif
       symbol = '!8u!DR!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   ;===========================================
   ; radial field
   ;===========================================
   endif else if(strcmp('br', name, /fold_case) eq 1) then begin

       psi0 = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, at_points=at_points, $
                       linear=linear)
       
       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       data = a_bracket(psi,psi0,x,y)/(r*sqrt(s_bracket(psi0,psi0,x,y)))
       symbol = '!8B!Dr!N!X'
       d = dimensions(/b0, _EXTRA=extra)


   ;===========================================
   ; poloidal flow
   ;===========================================
   endif else if(strcmp('vp', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.
        
       data = (s_bracket(phi,psi,x,y)/r + a_bracket(chi,psi,x,y)) / $ 
         sqrt(s_bracket(psi,psi,x,y))
       symbol = '!8u!Dp!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   ;===========================================
   ; surface flow
   ;===========================================
   endif else if(strcmp('vs', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       V   = read_field('V',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       psipsi = s_bracket(psi,psi,x,y)
       b2 = psipsi + i^2
       data = (i*s_bracket(phi,psi,x,y) - v*psipsi $
               +i*a_bracket(chi,psi,x,y)*r) / (r^2 * sqrt(psipsi*b2))
       symbol = '!8u!Ds!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   ;===========================================
   ; ideal_k
   ;===========================================
   endif else if(strcmp('k', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       n   = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       
       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.
         
       psipsi = s_bracket(psi,psi,x,y)

       data = n*(s_bracket(phi,psi,x,y) + r*a_bracket(chi,psi,x,y))/psipsi

   ;===========================================
   ; ideal omega
   ;===========================================
   endif else if(strcmp('ideal omega', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       v   = read_field('v',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       
       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       psipsi = s_bracket(psi,psi,x,y)
       
       data = $
         (v - $
          i*(s_bracket(phi,psi,x,y) + r*a_bracket(chi,psi,x,y))/psipsi) $
         /r^2

   ;===========================================
   ; Lundquist number
   ;===========================================
   endif else if(strcmp('S', name, /fold_case) eq 1) then begin

       va = read_field('va', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       eta= read_field('eta',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       
       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       data = va/eta
       symbol = '!8S!X'
       d = dimensions(_EXTRA=extra)
       

   ;===========================================
   ; b.W.b
   ;===========================================
   endif else if(strcmp('wpar', name, /fold_case) eq 1) then begin

       u   = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       v   = read_field('v'  , x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       psipsi = s_bracket(psi,psi,x,y)
       r2b2 = (psipsi + i^2)
       com = laplacian(chi,x,y,tor=itor)

       data = $
         (s_bracket(psi,a_bracket(u,psi,x,y)/r,x,y) $
          -0.5*r*a_bracket(u,psipsi/r^2,x,y) $
          -I*r*a_bracket(psi,v/r^2,x,y) $
          +0.5*r^2*s_bracket(chi,psipsi/r^2,x,y) $
          +com*psipsi $
          -s_bracket(psi,s_bracket(psi,chi,x,y),x,y))/r2b2 $
         - com/3.
           
       if(itor eq 1) then begin
           data = data + I^2*(dx(chi,x)/r-dz(u,y)/r^2)/r2b2
       endif

   ;===========================================
   ; radial electric field
   ;===========================================
   endif else if(strcmp('e_r', name, /fold_case) eq 1) then begin

       u   = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       v   = read_field('v'  , x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       b   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       pe  = read_field('pe',  x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       n   = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)
       eta = read_field('eta', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, at_points=at_points)

       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       jphi = grad_shafranov(psi,x,y,tor=itor)
       db = read_parameter('db', filename=filename)

       data = b*s_bracket(u,psi,x,y)/r^2 $
         + b*a_bracket(chi,psi,x,y)/r $
         - v*s_bracket(psi,psi,x,y)/r^2 $
         - (db/n)* $
         (s_bracket(psi,pe,x,y) $
          + (jphi*s_bracket(psi,psi,x,y) + b*s_bracket(psi,b,x,y))/r^2) $
         + eta*a_bracket(psi,b,x,y)/r

                                ; Normalize field
       data = -data/sqrt(s_bracket(psi,psi,x,y))
       symbol = '!8E!Dr!N!X'
       d = dimensions(/e0, _EXTRA=extra)

       
   endif else begin
       ; Field is not a composite field defined above;
       ; try to read it directly from C1.h5 file


       t = fltarr(trange[1]-trange[0]+1)
       file_id = h5f_open(filename)

       if((max(trange) ge 0) and $
          (ilin eq 1) and (not keyword_set(linear))) then begin
           base = read_field(name, x, y, t, slices=-1, mesh=mesh, $
                             filename=filename, points=pts, $
                             rrange=xrange, zrange=yrange, $
                             h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                             at_points=at_points,operation=op, $
                             last=last,symbol=symbol,units=units, $
                             cgs=cgs, mks=mks)
       endif else base = 0.

           

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

           if(n_elements(at_points) eq 0) then begin
               data[i-trange[0],*,*] = $
                 eval_field(field._data, mesh, points=pts, $
                            r=x, z=y, op=op, filename=filename, $
                            xrange=xrange, yrange=yrange, mask=mask) $
                 + base*(i ne -1)
           endif else begin
               print, sz
               for k=0, sz[1]-1 do begin
                   pos = [at_points[0,k]-xzero, $
                          at_points[1,k]-zzero]
                   data[i-trange[0],k] = $
                     eval_global(field._data,mesh,pos,elm=elm,op=op) $
                     + base*(i ne -1)
               end
           endelse
       end

       h5f_close, file_id

       symbol = translate(name, units=d, itor=itor)
   end

   if(n_elements(h_symmetry) eq 1) then begin
       data = (data + h_symmetry*reverse(data, 2)) / 2.
   endif
   if(n_elements(v_symmetry) eq 1) then begin
       print, "v symmetry = ", v_symmetry
       data = (data + v_symmetry*reverse(data, 3)) / 2.
   endif

   get_normalizations, filename=filename,b0=b0,n0=n0,l0=l0
   convert_units, data, d, b0, n0, l0, cgs=cgs, mks=mks
   convert_units, x, dimensions(/l0), b0, n0, l0, cgs=cgs, mks=mks
   convert_units, y, dimensions(/l0), b0, n0, l0, cgs=cgs, mks=mks
   units = parse_units(d, cgs=cgs, mks=mks)

   return, data
end


; field_at_point
function field_at_point, field, x, z, x0, z0
   i = n_elements(x)*(x0-min(x)) / (max(x)-min(x))
   j = n_elements(z)*(z0-min(z)) / (max(z)-min(z))
   return, interpolate(reform(field),i,j)
end


pro read_nulls, axis=axis, xpoints=xpoint, _EXTRA=extra
   s = read_scalars(_EXTRA=extra)

   t0 = get_slice_time(_EXTRA=extra)

   dum = min(s.time._data - t0, i, /abs)

   xpoint = fltarr(2)
   axis = fltarr(2)
   xpoint[0] = s.xnull._data[i]
   xpoint[1] = s.znull._data[i]
   axis[0] = s.xmag._data[i]
   axis[1] = s.zmag._data[i]
   
   return
end


function path_at_flux, psi,x,z,t,flux,breaks=breaks
   contour, psi[0,*,*], x, z, levels=flux, closed=0, $
     path_xy=xy, path_info=info, /path_data_coords, /overplot

   if(n_elements(xy) eq 0) then begin
       print, 'Error: no points at this flux value', flux
       return, 0
   end

; refine the surface found by contour with a single newton iteration
;    i = n_elements(x)*(xy[0,*]-min(x)) / (max(x)-min(x))
;    j = n_elements(z)*(xy[1,*]-min(z)) / (max(z)-min(z))
;    f = interpolate(reform(psi[0,*,*]),i,j)
;    fx = interpolate(reform(dx(psi[0,*,*],x)),i,j)
;    fz = interpolate(reform(dz(psi[0,*,*],z)),i,j)
;    gf2 = fx^2 + fz^2
;    l = (flux-f)/gf2
;    xy[0,*] = xy[0,*] + l*fx
;    xy[1,*] = xy[1,*] + l*fz
  
   ; find breaks
   
   dx = sqrt(deriv(xy[0,*])^2 + deriv(xy[1,*])^2)
   minforbreak = 10.*median(dx)

   breaks = where(dx gt minforbreak,count)

   return, xy
end


; ==============================================================
; find_nulls
; ----------
;
;  Finds field nulls.
;  xpoint = fltarr(2) : position of active x-point
;  axis   = fltarr(2) : position of magnetic axis
; ==============================================================
pro find_nulls, psi, x, z, axis=axis, xpoints=xpoint, _EXTRA=extra

   if(n_elements(time) eq 0) then time = 0
   
   field = s_bracket(psi,psi,x,z)
   d2 = dz(dz(psi,z),z)*dx(dx(psi,x),x)
   
   nulls = field lt mean(field)/10.

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
           foundlocalmin = 0

           ; find local minimum
           for m=i, sz[2]-1 do begin
               if(nulls[0,m,j] eq 0) then break
               for n=j+1, sz[3]-1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] le currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                       foundlocalmin = 1
                   endif
                   nulls[0,m,n] = 0
               endfor
               for n=j-1, 0, -1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] le currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                       foundlocalmin = 1
                   endif
                   nulls[0,m,n] = 0
               endfor
           endfor

           if (foundlocalmin eq 0) then continue

           ; throw out local minima on boundaries
           if (currentpos[0] eq 0) or (currentpos[0] eq sz[2]-1) then continue
           if (currentpos[1] eq 0) or (currentpos[1] eq sz[3]-1) then continue

           ; determine if point is an x-point or an axis and
           ; append the location index to the appropriate array
           if(d2[0,currentpos[0],currentpos[1]] lt 0) then begin
               if(psi[0,currentpos[0],currentpos[1]] gt oldxflux) then begin
                   xpoint = [x[currentpos[0]], z[currentpos[1]]]
                   oldxflux = psi[0,currentpos[0],currentpos[1]]
               end
           endif else begin
               if(psi[0,currentpos[0],currentpos[1]] gt oldxflux) then begin
                   axis = [x[currentpos[0]], z[currentpos[1]]]
                   oldaflux = psi[0,currentpos[0],currentpos[1]]
               endif
           endelse
       endfor
   endfor 

   if(n_elements(axis) lt 2) then begin
       print, 'Warning: cannot find magnetic axis'
   endif

   fieldr = dx(field,x)
   fieldz = dz(field,z)

   ; do iterative refinement on magnetic axis
   for k=1, 2 do begin
       pt = field_at_point(field, x, z, axis[0], axis[1])
       pt1 = field_at_point(fieldr, x, z, axis[0], axis[1])
       pt2 = field_at_point(fieldz, x, z, axis[0], axis[1])

       denom = pt1^2 + pt2^2           
       dx = pt*pt1/denom
       dz = pt*pt2/denom           

       axis[0] = axis[0] - dx
       axis[1] = axis[1] - dz
   end

   print, 'Found axis at ', axis[0], axis[1]


   if(n_elements(xpoint) ge 2) then begin
       print, 'x-point guess at ', xpoint[0], xpoint[1]
       ; do iterative refinement on x-point
       for k=1, 10 do begin
           pt = field_at_point(field, x, z, xpoint[0], xpoint[1])
           pt1 = field_at_point(fieldr, x, z, xpoint[0], xpoint[1])
           pt2 = field_at_point(fieldz, x, z, xpoint[0], xpoint[1])
           
           denom = pt1^2 + pt2^2           
           dx = pt*pt1/denom
           dz = pt*pt2/denom           
           
           xpoint[0] = xpoint[0] - dx
           xpoint[1] = xpoint[1] - dz

           if(xpoint[0] le min(x) or xpoint[0] ge max(x)) then begin
               xpoint = 0
               break
           endif
           if(xpoint[1] le min(z) or xpoint[1] ge max(z)) then begin
               xpoint = 0
               break
           endif           
       end
   endif
   if(n_elements(xpoint) eq 2) then begin
       print, 'Found X-point at ', xpoint[0], xpoint[1]
   endif else begin
       print, 'No X-point found'
   endelse
end


pro nulls, psi, x, z, axis=axis, xpoints=xpoint, _EXTRA=extra
    version = read_parameter('version', _EXTRA=extra)
    if(version ge 3) then begin
        read_nulls, axis=axis, xpoint=xpoint, _EXTRA=extra
    endif else begin
        find_nulls, psi, x, z, axis=axis, xpoint=xpoint, _EXTRA=extra
    endelse 
end

; ========================================================
; find_lcfs
; ~~~~~~~~~
;
; returns the flux value of the last closed flux surface
; ========================================================
function find_lcfs, psi, x, z, axis=axis, xpoint=xpoint, $
               flux0=flux0, psilim=psilim, xlim=xlim, _EXTRA=extra

   if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) then $
     psi = read_field('psi',x,z,t,_EXTRA=extra,linear=0)

   nulls, psi, x, z, xpoint=xpoint, axis=axis, _EXTRA=extra

   ; flux at magnetic axis
   if(n_elements(axis) lt 2) then begin
       print, "Error: no magnetic axis"
       flux0 = max(psi[0,*,*])
   endif else begin
       flux0 = field_at_point(psi,x,z,axis[0],axis[1])
   endelse
   if(n_elements(axis) gt 2) then begin
       print, "Warning: there is more than one magnetic axis"
   endif
   print, "Flux on axis:", flux0

   ; If xlim is set, then use limiter values
   ; given in C1input file
   if(keyword_set(xlim)) then begin
       limiter = fltarr(2)
       limiter[0] = read_parameter('xlim', _EXTRA=extra)
       limiter[1] = read_parameter('zlim', _EXTRA=extra)
       print, "limiter at: ", limiter
   endif

   ; limiting value
   ; Find limiting flux by calculating outward normal derivative of
   ; the normalized flux.  If this derivative is negative, there is a
   ; limiter.
   if(n_elements(psilim) eq 0) then begin
       if(n_elements(limiter) eq 2) then begin
           xerr = min(x-limiter[0],xi,/absolute)
           zerr = min(z-limiter[1],zi,/absolute)
           psilim = psi[0,xi,zi]
       endif else begin
           print, ' No limiter provided, using wall.'
           sz = size(psi)
           
           psiz = dz(psi,z)
           psix = dx(psi,x)
           
           normal_mask = psi*0.
           normal_mask[0,      *,      0] = 1.
           normal_mask[0,      *,sz[3]-1] = 1.
           normal_mask[0,      0,      *] = 1.
           normal_mask[0,sz[2]-1,      *] = 1.
           
           xx = fltarr(1,sz[2],sz[3])
           zz = fltarr(1,sz[2],sz[3])
           for i=0,sz[3]-1 do xx[0,*,i] = x - axis[0]
           for i=0,sz[2]-1 do zz[0,i,*] = z - axis[1]
           normal_deriv = (psix*xx + psiz*zz)*normal_mask
           
           normal_deriv = normal_deriv lt 0
           
           psi_bound = psi*normal_deriv + (1-normal_deriv)*1e10
           
           psilim = min(psi_bound-flux0, i, /absolute)
           psilim = psi_bound[i]
       endelse
       
       print, "Flux at limiter", psilim
       
       
       ; flux at separatrix
       sz = size(xpoint)
       if(sz[0] gt 0 and (not keyword_set(xlim))) then begin
           psix = field_at_point(psi,x,z,xpoint[0],xpoint[1])
           print, "Flux at separatrix:", psix
           
           if(abs(psix-flux0) gt abs(psilim-flux0)) then begin
               print, "Plasma is limited."
           endif else begin
               print, "Plasma is diverted."
               psilim = psix
           endelse
       endif
   end
   
   return, psilim
end


function read_lcfs, axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra
   s = read_scalars(_EXTRA=extra)

   t0 = get_slice_time(_EXTRA=extra)

   dum = min(s.time._data - t0, i, /abs)

   xpoint = fltarr(2)
   axis = fltarr(2)
   xpoint[0] = s.xnull._data[i]
   xpoint[1] = s.znull._data[i]
   axis[0] = s.xmag._data[i]
   axis[1] = s.zmag._data[i]
   
   flux0 = s.psimin._data[i]

   return, s.psi_lcfs._data[i]
end


; ========================================================
; lcfs
; ~~~~
;
; returns the flux value of the last closed flux surface
; ========================================================
function lcfs, psi, x, z, axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra

    version = read_parameter('version', _EXTRA=extra)
    if(version ge 3) then begin
        psilim = read_lcfs(axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra)
    endif else begin
        psilim = find_lcfs(psi, x, z,axis=axis, xpoint=xpoint, flux0=flux0, $
                           _EXTRA=extra)
    endelse 

    print, 'LCFS: '
    print, ' Magnetic axis found at ', axis
    print, ' Active x-point at ', xpoint
    print, ' psi_0, psi_s = ', flux0, psilim
    
    return, psilim
end




pro plot_flux_contour, fval, _EXTRA=extra
   n = n_elements(fval)

   psi = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)

   contour, psi, x, z, closed=0, levels=fval(sort(fval)), _EXTRA=extra
end


pro compare_fields, file1, file2, time, names=names

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

   v = read_parameter('version', filename=filename) 
   print, 'version = ', v
   if(v lt 2) then begin
       names = ['- dN/dt', 'Diffusive', 'Convective', 'Sources']
       comp = fltarr(n_elements(names), n_elements(t))
       comp[0,*] = -deriv(t,s.particle_number._data)
   endif else begin
       names = ['dN/dt', 'Diffusive', 'Convective', 'Sources']
       comp = fltarr(n_elements(names), n_elements(t))
       comp[0,*] = deriv(t,s.particle_number._data)
   endelse
   comp[1,*] = s.Particle_Flux_diffusive._data
   comp[2,*] = s.Particle_Flux_convective._data
   comp[3,*] = s.Particle_source._data

   return, total(comp, 1)
end


;==============================================================
; momentum flux
; ~~~~~~~~~~~~~
;
; Returns a time series of the total rate of toroidal (angular)
; momenumt lost through the simulation domain boundary.
; 
; Optional outputs:
;  t: the time array
;  names: the name of each flux component
;  components: an array of time series of each flux component
;==============================================================
function momentum_flux, filename=filename, components=comp, names=names, t=t, $
                        norm=norm
   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   t = s.time._data

   v = read_parameter('version', filename=filename) 
   if(v lt 2) then begin
       names = ['- dL/dt', 'Magnetic', 'Solenoidal', $
                'Compressional', 'Viscous', 'Gyroviscous']
       comp = fltarr(n_elements(names), n_elements(t))
       comp[0,*] = -deriv(t, s.angular_momentum._data)
   endif else begin
       names = ['- dL/dt', 'Magnetic', 'Solenoidal', $
                'Compressional', 'Viscous', 'Gyroviscous', 'Parallel Viscous']
       comp = fltarr(n_elements(names), n_elements(t))
       comp[0,*] = deriv(t, s.angular_momentum._data)
       comp[6,*] = s.Torque_parvisc._data
   endelse
  
   comp[1,*] = s.Torque_em._data
   comp[2,*] = s.Torque_sol._data
   comp[3,*] = s.Torque_com._data
   comp[4,*] = s.Torque_visc._data
   comp[5,*] = s.Torque_gyro._data

   if(keyword_set(norm)) then comp = comp/max(abs(s.angular_momentum._data))

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
; energy
; ~~~~~~
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

   en = energy(filename=filename, components=en_comp, names=en_names, t=t)

   v = read_parameter('version', filename=filename) 
   if(v lt 2) then begin
       names = ['- dE/dt', 'Pressure', 'Kinetic', 'Poynting', 'Thermal']
       comp = fltarr(n_elements(names), n_elements(t))
       comp[0,*] = -deriv(t, en)
   endif else begin
       names = ['dE/dt', 'Pressure', 'Kinetic', 'Poynting', 'Thermal']
       comp = fltarr(n_elements(names), n_elements(t))
       comp[0,*] = deriv(t, en)
   endelse

   comp[1,*] = s.Flux_pressure._data
   comp[2,*] = s.Flux_kinetic._data
   comp[3,*] = s.Flux_poynting._data
   comp[4,*] = s.Flux_thermal._data

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
   endif else if(strcmp(name, 'momentum', /fold_case) eq 1) then begin
       tot = momentum_flux(filename=filename, comp=comp, names=names, t=t, $
                           norm=norm)
       norm_tot = 1.
       title  = '!6Momentum Flux!X'
       ytitle = '!6Momentum Flux!X'
       units = ''
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

;    if(keyword_set(per_length)) then begin
;        itor = read_parameter('itor', filename=filename)
;        if(itor eq 1) then begin
;            xzero = read_parameter('xzero', filename=filename)
;            tot = tot / xzero
;            comp = comp / xzero
;        endif
;    endif

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


; ========================================================
; plot_lcfs
; ~~~~~~~~~
;
; plots the last closed flux surface
; ========================================================
pro plot_lcfs, psi, x, z, psival=psival,_EXTRA=extra

    if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) then begin
        print, 'reading psi, plot_lcfs'
        psi = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
    end

    ; if psival not passed, choose limiter value
    if(n_elements(psival) eq 0) then psival = lcfs(psi,x,z,_EXTRA=extra)

    ; plot contour
    loadct, 12
;    contour, psi, x, z, /overplot, nlevels=1, levels=psival, $
;      thick=!p.charthick*2., color=color(3,5)
    xy = path_at_flux(psi, x, z, t, psival, breaks=breaks)

    n = n_elements(breaks)
    nxy = n_elements(xy[0,*])

    if(breaks[0] ne -1) then begin
        oplot, xy[0,0:breaks[0]], xy[1,0:breaks[0]], $
          thick=!p.thick*3, color=color(3,5)
        if(n ge 2) then begin
            for i=1, n-2 do begin
                oplot, xy[0,breaks[i-1]:breaks[i]], $
                  xy[1,breaks[i-1]:breaks[i]], $
                  thick=!p.thick*3, color=color(3,5)
            endfor
        endif
        oplot, xy[0,breaks[n-1]:nxy-1], xy[1,breaks[n-1]:nxy-1], $
          thick=!p.thick*3, color=color(3,5)

    endif else begin
        oplot, xy[0,*], xy[1,*], thick=!p.thick*3, color=color(3,5)
    endelse
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
       plot, timings.t_onestep._data>0, title='!6Timings!3', $
         xtitle='!6Time Step!3', ytitle='!8t!6 (s)!3', _EXTRA=extra
   endelse
   oplot, timings.t_ludefall._data, linestyle=2, color=color(1,8)
   oplot, timings.t_sources._data, linestyle=1, color=color(2,8)
   oplot, timings.t_aux._data, linestyle=1, color=color(3,8)
   oplot, timings.t_smoother._data, linestyle=1, color=color(4,8)
   oplot, timings.t_mvm._data, linestyle=1, color=color(5,8)
   oplot, t_solve, linestyle=2, color=color(6,8)
   oplot, t_output, linestyle=2, color=color(7,8)


   plot_legend, ['Onestep', 'ludefall', 'sources', 'aux', $
                 'smoother', 'mat vec mult', 'solve', 'output'], $
     linestyle=[0,2,1,1,1,1,2,2], color=colors(8)

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

; ======================================
; beta_toroidal = 2*<P>/B_T0^2
; ======================================
function beta_toroidal, filename=filename

print, 'filename = ', filename

   gamma = read_parameter('gam', filename=filename)
   xmag = read_parameter('xmag', filename=filename)
   rzero = read_parameter('rzero', filename=filename)
   if(rzero eq 0) then rzero = xmag
   bzero = read_parameter('bzero', filename=filename)

   bt0 = bzero*(rzero/xmag)
   
   s = read_scalars(filename=filename)

   beta_t = 2.*(gamma-1.)*s.E_P._data/bt0^2
   
   print, 'bt0 =', bt0
   print, 'rzero = ', rzero
   print, 'bzero = ', bzero

   return, beta_t
end


; ====================================
; beta_normal = beta_t * (B_T*a / I_p)
; ====================================
function beta_normal, filename=filename

   a = 1.

   gamma = read_parameter('gam', filename=filename)
   xmag = read_parameter('xmag', filename=filename)
   rzero = read_parameter('rzero', filename=filename)
   if(rzero eq 0) then rzero = xmag
   bzero = read_parameter('bzero', filename=filename)

   bt0 = bzero*(rzero/xmag)
   
   s = read_scalars(filename=filename)
   ip = s.toroidal_current._data

   beta_t = 2.*(gamma-1.)*s.E_P._data/bt0^2
   beta_n = beta_t * abs(bt0*a/ip)

   print, 'ip', ip
   
   return, 100.*beta_n
end

function read_scalar, scalarname, filename=filename, title=title, $
                      symbol=symbol, units=units, time=time, final=final, $
                      _EXTRA=extra

   if(n_elements(scalarname) eq 0) then begin
       print, "Error: no scalar name provided"
       return, 0
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(n_elements(filename) gt 1) then begin
       data = fltarr(n_elements(filename))
       for i=0, n_elements(filename)-1 do begin
           data[i] = read_scalar(scalarname, filename=filename[i], $
                                 title=title, symbol=symbol, units=units, $
                                 time=time, /final)
       end
       return, data
   endif

   s = read_scalars(filename=filename)
   itor = read_parameter('itor', filename=filename)
   if(n_tags(s) eq 0) then return, 0

   time = s.time._data

   d = dimensions()

   if(strcmp("toroidal current", scalarname, /fold_case) eq 1) or $
     (strcmp("it", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current._data
       title = 'Toroidal Current'
       symbol = '!8I!DT!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if (strcmp("toroidal flux", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_flux._data
       title = 'Toroidal Flux'
       symbol = 'Flux'
       d = dimensions(/b0, l0=2, _EXTRA=extra)
   endif else $
     if (strcmp("reconnected flux", scalarname, /fold_case) eq 1) then begin
       data = abs(s.reconnected_flux._data)
       title = 'Reconnected Flux'
       symbol = translate('psi')
       d = dimensions(/b0, l0=1+itor, _EXTRA=extra)
   endif else $
     if (strcmp("loop voltage", scalarname, /fold_case) eq 1) or $
     (strcmp("vl", scalarname, /fold_case) eq 1) then begin
       data = s.loop_voltage._data
       title = 'Loop Voltage'
       symbol = '!8V!DL!N!X'
       d = dimensions(/e0,/l0, _EXTRA=extra)
   endif else $
     if (strcmp("beta", scalarname, /fold_case) eq 1) then begin
       data = beta(filename=filename)
       title = 'Average Beta'
       symbol = '!7b!X'
       d = dimensions(_EXTRA=extra)
   endif else if $
     (strcmp("poloidal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bp", scalarname, /fold_case) eq 1) then begin
       data = beta_poloidal(filename=filename)
       title = 'Poloidal Beta'
       symbol = '!7b!D!8P!N!X'
       d = dimensions(_EXTRA=extra)
   endif else $
     if (strcmp("normal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bn", scalarname, /fold_case) eq 1) then begin
       data = beta_normal(filename=filename)
       title = 'Normal Beta'
       symbol = '!7b!D!8N!N!3'
       d = dimensions(_EXTRA=extra)
   endif else if $
     (strcmp("toroidal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bt", scalarname, /fold_case) eq 1) then begin
       data = beta_toroidal(filename=filename)
       title = 'Toroidal Beta'
       symbol = '!7b!D!8T!N!X'
   endif else if $
     (strcmp("kinetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("ke", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_KP._data 
       if(nv ge 2) then data = data + s.E_KT._data 
       if(nv ge 3) then data = data + s.E_K3._data
       title = 'Kinetic Energy'
       symbol = '!8KE!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("magnetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("me", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_MP._data 
       if(nv ge 2) then data = data + s.E_MT._data 
       title = 'Magnetic Energy'
       symbol = '!8ME!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("thermal energy", scalarname, /fold_case) eq 1) or $
     (strcmp("te", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_P._data 
       title = 'Thermal Energy'
       symbol = '!8TE!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("energy", scalarname, /fold_case) eq 1) then begin
       nv = read_parameter("numvar", filename=filename)
       data = energy(filename=filename)
       title = 'Total Energy'
       symbol = '!8E!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("particles", scalarname, /fold_case) eq 1) or $
     (strcmp("n", scalarname, /fold_case) eq 1) then begin
       data = s.particle_number._data
       title = 'Particle Number'
       symbol = '!8N!X'
       d = dimensions(/n0,l0=3, _EXTRA=extra)
   endif else if $
     (strcmp("angular momentum", scalarname, /fold_case) eq 1) then begin
       data = s.angular_momentum._data
       title = 'Angular Momentum'
       symbol = '!8L!D!9P!N!X'
       d = dimensions(/energy,/t0, _EXTRA=extra)
   endif else if (strcmp("circulation", scalarname, /fold_case) eq 1) or $
     (strcmp("vorticity", scalarname, /fold_case) eq 1) then begin
       data = s.circulation._data
;      data = s.vorticity._data
       title = 'Circulation'
       symbol = '!9I!S!7x!R!A!6_!D!9P!N !8dA !x'
       d = dimensions(/v0,/l0, _EXTRA=extra)
   endif else if (strcmp("parallel viscous heating", scalarname, /fold_case) eq 1) then begin $
       data = s.parallel_viscous_heating._data
       title = 'Parallel Viscous Heating'
;       symbol = '!6(!7l!D!9#!N!5b!9 . !3W!9 . !5b!6/2)!U2!N!X'
       symbol = '-!9I!8dV !7P!D!9#!N!3:!9G!5u!X'
       d = dimensions(/energy,t0=-1, _EXTRA=extra)
   endif else if (strcmp("bwb2", scalarname, /fold_case) eq 1) then begin $
       amupar = read_parameter('amupar', filename=filename)
       data = 4.*s.parallel_viscous_heating._data / (3.*amupar)
       title = ''
;       symbol = '!6(!7l!D!9#!N!5b!9 . !3W!9 . !5b!6/2)!U2!N!X'
       symbol = '!6(!5b!9.!3W!9.!5b!6)!U2!N!X'
       d = dimensions(t0=-2,l0=3, _EXTRA=extra)
   endif else begin
       print, 'Scalar ', scalarname, ' not recognized.'
       return, 0
   endelse
   
   if(keyword_set(final)) then begin
       data = data[n_elements(data)-1]
   endif

   get_normalizations, b0=b0,n0=n0,l0=l0, filename=filename, _EXTRA=extra
   convert_units, data, d, b0, n0, l0, _EXTRA=extra
   convert_units, time, dimensions(/t0), b0, n0, l0, _EXTRA=extra
   units = parse_units(d, _EXTRA=extra)

   return, data
end

pro plot_scalar, scalarname, x, filename=filename, names=names, $
                 _EXTRA=extra, overplot=overplot, $
                 ylog=ylog, xlog=xlog, absolute_value=absolute, $
                 power_spectrum=power_spectrum, per_length=per_length, $
                 growth_rate=growth_rate, bw=bw, nolegend=nolegend, $
                 cgs=cgs,mks=mks,linestyle=ls, color=co

  if(n_elements(filename) eq 0) then filename='C1.h5'

  if(n_elements(names) eq 0) then names=filename

  nfiles = n_elements(filename)
  if(nfiles gt 1) then begin
      if(keyword_set(bw)) then begin
          if(n_elements(ls) eq 0) then ls = indgen(nfiles)
          if(n_elements(co) eq 0) then co = replicate(color(0,1),nfiles)
      endif else begin
          if(n_elements(ls) eq 0) then ls = replicate(0, nfiles)
          if(n_elements(co) eq 0) then co = colors(nfiles)
      endelse

      for i=0, nfiles-1 do begin
          if(n_elements(x) eq 0) then begin
              plot_scalar, scalarname, filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=co[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=power_spectrum, per_length=per_length, $
                growth_rate=growth_rate, linestyle=ls[i], nolegend=nolegend, $
                absolute_value=absolute,cgs=cgs,mks=mks
          endif else begin
              plot_scalar, scalarname, x[i], filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=power_spectrum, per_length=per_length, $
                growth_rate=growth_rate, nolegend=nolegend, $
                absolute_value=absolute,cgs=cgs,mks=mks
          endelse
      end

      if((n_elements(names) gt 0) and (not keyword_set(nolegend))) then begin
          plot_legend, names, ylog=ylog, xlog=xlog, $
            color=co, linestyle=ls, _EXTRA=extra
      endif    

      return
  endif 

  data = read_scalar(scalarname, filename=filename, time=time, $
                     title=title, symbol=symbol, units=units, cgs=cgs, mks=mks)
  if(n_elements(data) le 1) then return

  title = '!6' + title + '!X'

  ytitle = symbol
  if(strlen(units) gt 0) then begin
      ytitle = ytitle + '!6 (' + units + ')!X'
  endif

  if(keyword_set(power_spectrum)) then begin
;      xtitle = '!7x!6 (!7s!D!8A!N!6!U-1!N)!X'
      xtitle = make_label('!7x!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)
      data = power_spectrum(data, frequency=tdata, t=max(time)) 
  endif else begin
      xtitle = make_label('!8t!X', /t0, cgs=cgs, mks=mks, _EXTRA=extra)
      tdata = time
  endelse

  if(keyword_set(per_length)) then begin
      itor = read_parameter('itor', filename=filename)
      if(itor eq 1) then begin
          rzero = read_parameter('rzero', filename=filename)
          data = data / rzero
      endif
  endif
  
  if(keyword_set(growth_rate)) then begin
      data = deriv(tdata, alog(abs(data)))
;      ytitle = '!7c !6(!7s!D!8A!N!6!U-1!N)!X'
      ytitle = make_label('!7c!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)
  endif

  if(keyword_set(absolute)) then data = abs(data)

  if(n_elements(x) eq 0) then begin   
      if(keyword_set(overplot)) then begin
          oplot, tdata, data, color=co, linestyle=ls, _EXTRA=extra
      endif else begin
          plot, tdata, data, xtitle=xtitle, ytitle=ytitle, $
            title=title, _EXTRA=extra, ylog=ylog, xlog=xlog, $
            color=co, linestyle=ls
      endelse
  endif else begin
      xi = x
      x = fltarr(1)
      z = fltarr(1)
      x[0] = xi
      z[0] = data[n_elements(data)-1]

      if(keyword_set(overplot)) then begin
          oplot, x, z, color=co, linestyle=ls, _EXTRA=extra
      endif else begin
          plot, x, z, $
            title=title, xtitle=xtitle, ytitle=ytitle, $
            _EXTRA=extra, ylog=ylog, xlog=xlog, color=co, $
            linstyle=ls
      endelse
  endelse
end


; ==================================================
; plot_pol_velocity
; ~~~~~~~~~~~~~~~~~
;
; makes a vector plot of the poloidal velocity
; ==================================================
pro plot_pol_velocity, time,  maxval=maxval, points=points, $
                       lcfs=lcfs, _EXTRA=extra

  if(n_elements(pts) eq 0) then pts=25

  nv = read_parameter('numvar', _EXTRA=extra)
  itor = read_parameter('itor', _EXTRA=extra)
  if(n_elements(itor) gt 1) then itor=itor[0]

  phi = read_field('phi', x, z, t, points=points, _EXTRA=extra, slice=time)
  if(n_elements(phi) le 1) then return

  if(itor eq 1) then r = radius_matrix(x,z,t) else r = 1.

  vx = -dz(phi,z)/r
  vz =  dx(phi,x)/r

  chi = read_field('chi', x, z, t, points=points, _EXTRA=extra, slice=time)
  vx = vx + dx(chi,x)
  vz = vz + dz(chi,z)

  bigvel = max(sqrt(vx^2 + vz^2))
  print, "maximum velocity: ", bigvel
  if(n_elements(maxval) ne 0) then begin
      length = bigvel/maxval
  endif else length=1

  if(n_elements(title) eq 0) then begin
      title = '!6Poloidal Flow!X'
      if(t gt 0) then begin
          title = title +  $
            string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!60!N)!X")', t)
      endif else begin
          title = title + $
            string(FORMAT='("!6(!8t!6 = ",G0,")!X")', t)
      endelse
  endif
  
  maxstr=string(format='("!6max(!8u!Dpol!N!6) = ",G0.3,"!X")',bigvel) + $
    '!6 ' + make_units(/v0) + '!X'

  velovect, reform(vx), reform(vz), x, z, length=length, _EXTRA=extra, $
    xtitle=make_label('!8R!X', /l0, _EXTRA=extra), $
    ytitle=make_label('!8Z!X', /l0, _EXTRA=extra), $
    title=title, subtitle=maxstr

  if(keyword_set(lcfs)) then plot_lcfs, points=200, _EXTRA=extra
end


; ==============================================
; field_at_flux
;
; evaluates field at points where psi=flux
; ==============================================
function field_at_flux, field, psi, x, z, t, flux, theta=theta, angle=angle, $
                      xp=xp, zp=zp, integrate=integrate, axis=axis, $
                        _EXTRA=extra

   if(n_elements(field) le 1) then return, 0

   xy = path_at_flux(psi,x,z,t,flux)
   if(n_elements(xy) le 1) then return, 0

   i = n_elements(x)*(xy[0,*]-min(x)) / (max(x)-min(x))
   j = n_elements(z)*(xy[1,*]-min(z)) / (max(z)-min(z))

   test = interpolate(reform(field[0,*,*]),i,j)
   if(n_elements(theta) ne 0) then begin
       angle = interpolate(reform(theta[0,*,*]),i,j)
   endif else begin
       if(n_elements(axis) eq 0) then begin
           angle = findgen(n_elements(test))
       endif else begin
           angle = atan(xy[1,*]-axis[1],xy[0,*]-axis[0])
       endelse
   endelse

   ; re-order array
   mint = min(angle,p)
   p = p+1
   test = shift(test, -p)
   angle = shift(angle, -p)
   xp = shift(xy[0,*],-p)
   zp = shift(xy[1,*],-p)

   if(keyword_set(integrate)) then begin
       r = radius_matrix(x,z,t)
       jac = interpolate(r/a_bracket(theta,psi,x,z),i,j)
       jac = shift(jac, -p)
       dtheta = deriv(angle)
       test = total(test*dtheta*jac, /cumulative)
       d = min(angle, i,/abs)
       test = test - test[i]
   endif

   return, test
end


;==================================================================
; flux_coord_field
; ~~~~~~~~~~~~~~~~
;==================================================================
function flux_coord_field, field, psi, x, z, t, slice=slice, area=area, $
                           fbins=fbins,  tbins=tbins, flux=flux, angle=angle, $
                           psirange=frange, nflux=nflux, pest=pest, $
                           _EXTRA=extra

   if(n_elements(psi) eq 0) then begin
       linear = read_parameter('linear',_EXTRA=extra)
       if(linear eq 1) then begin
           psi = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
       endif else begin
           psi = read_field('psi',x,z,t,slice=slice,_EXTRA=extra)
       endelse
   endif

   sz = size(field)

   if(n_elements(fbins) eq 0) then fbins = 10
   if(n_elements(tbins) eq 0) then tbins = 10

   result = fltarr(sz[1], fbins, tbins)
   flux = fltarr(sz[1], fbins)
   angle = fltarr(sz[1], tbins)
   area = fltarr(sz[1], fbins)

   psival = lcfs(psi,x,z,axis=axis,xpoint=xpoint,slice=slice,flux0=flux0, $
                 _EXTRA=extra)
   
   if(n_elements(range) eq 0) then begin
       ; if range not provided, use all flux within lcfs
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = [psival, max(psi[k,*,*])]
   endif else if(n_elements(range) eq 2) then begin
       oldrange = range
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = oldrange
   endif

   r = radius_matrix(x,z,t)

   if(keyword_set(pest)) then begin
       linear = read_parameter('linear',_EXTRA=extra)
       if(linear eq 1) then begin
           i = read_field('i',x,z,t,slice=-1,_EXTRA=extra)
       endif else begin
           i = read_field('i',x,z,t,slice=slice,_EXTRA=extra)
       endelse
       bp = sqrt(s_bracket(psi,psi,x,z)/r^2)
       bt = i/r
       db = bt/(r*bp)
   endif

   print, 'binning with fbins, tbins=', fbins, tbins

   for k=0, sz[1]-1 do begin

       angle[k,*] = 2.*!pi*findgen(tbins)/float(tbins) - !pi
       dpsi = float(range[k,1] - range[k,0])/float(fbins)

       for p=0, fbins-1 do begin
           flux[k,p] = range[k,1] - dpsi*(p+0.5)
       
           f = field_at_flux(field, psi, x, z, t, flux[k,p], $
                             angle=a, xp=xp, zp=zp, axis=axis)

           dx = deriv(xp)
           dz = deriv(zp)
           ds = sqrt(dx^2 + dz^2)
           area[k,p] = total(ds*xp)

           if(keyword_set(pest)) then begin
               g = field_at_flux(db, psi, x, z, t, flux[k,p], $
                                 angle=a, axis=axis)
               dum = min(a, i, /abs)
               
               dt = ds*g
               
               a = total(dt,/cum)
               a = 2.*!pi*a/(max(a)-min(a))
               a = a + dum - a[i]
           endif

           result[k,p,*] = interpol(f,a,angle[k,*])
       end
   endfor

   nflux = (flux-flux0)/(psival-flux0)
       
   return, result
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
                             area=area, dV=dV, psirange=range, $
                             integrate=integrate, r0=r0, $
                             nflux=nflux, _EXTRA=extra

   sz = size(field)

   points = sqrt(sz[2]*sz[3])

   if(n_elements(bins) eq 0) then bins = fix(points/4)

   print, 'flux averaging with ', bins, ' bins'

   result = fltarr(sz[1], bins)
   flux = fltarr(sz[1], bins)
   dV = fltarr(sz[1], bins)
   area = fltarr(sz[1], bins)

   psival = lcfs(psi, x, z, axis=axis, xpoint=xpoint, flux0=flux0,_EXTRA=extra)
   r0 = axis[0]

   if(n_elements(range) eq 0) then begin
       ; if range not provided, use all flux within lcfs
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = [psival, flux0]
   endif else if(n_elements(range) eq 2) then begin
       oldrange = range
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = oldrange
   endif

   ; remove divertor region from consideration
   div_mask = fltarr(n_elements(x), n_elements(z))
   if(n_elements(xpoint) ge 2) then begin
       if(xpoint[0] ne 0. or xpoint[1] ne 0.) then begin
           if(xpoint[1] lt axis[1]) then begin
               for i=0, n_elements(z)-1 do begin
                   if(z[i] ge xpoint[1]) then div_mask[*,i] = 1.
               end
               print, ' removing points below z=', xpoint[1]
           endif else begin
               for i=0, n_elements(z)-1 do begin
                   if(z[i] le xpoint[1]) then div_mask[*,i] = 1.
               end
               print, ' removing points above z=', xpoint[1]
           endelse
       endif else div_mask = 1.
   endif else begin
       div_mask = 1.
   endelse

   bp = sqrt(s_bracket(psi,psi,x,z))

   for k=0, sz[1]-1 do begin
       dpsi = float(range[k,1] - range[k,0])/float(bins)

       for p=0, bins-1 do begin
           fval = range[k,1] - dpsi*(p+1)
           flux[k,p] = fval + dpsi/2.

           faf = field_at_flux(field[k,*,*], psi[k,*,*], $
                               x, z, t, xp=xp, zp=zp, flux[k,p])
           bpf = field_at_flux(bp[k,*,*], psi[k,*,*], $
                               x, z, t, xp=xp, zp=zp, flux[k,p])

           if(n_elements(bpf) lt 3) then continue

           dl = sqrt(deriv(xp)^2 + deriv(zp)^2)

           dv[k,p] = total(dl*xp/bpf)
           area[k,p] = total(dl*xp)
           if(keyword_set(integrate)) then begin
               result[k,p] = total(faf*xp*dl*(-dpsi)/bpf)
           end else begin
;               result[k,p] = total(faf*xp*dl) / total(xp*dl)
               result[k,p] = total(faf*xp*dl/bpf)/total(xp*dl/bpf)
           endelse
       endfor
       if(keyword_set(integrate)) then begin
           result[k,*] = total(result[k,*], /cumulative)
       endif
   endfor

   nflux = (flux0-flux)/(flux0-psival)

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
;     dV:  differential volume of each flux surface
;   area:  surface area of each flux surface
;   name:  the formatted name of the field
; symbol:  the formatted symbol of the field
;  units:  the formatted units of the field
;==================================================================
function flux_average, field, psi=psi, x=x, z=z, t=t, r0=r0, $
                       flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
                       points=points, name=name, symbol=symbol, units=units, $
                       integrate=integrate, _EXTRA=extra

   type = size(field, /type)

   sz = size(field)

   if(n_elements(points) eq 0) then begin
       if(type ne 7) then points = sz[2] else points=50
   endif

   if((n_elements(psi) le 1) $
      or (n_elements(x) eq 0) or (n_elements(z) eq 0) $
      or (n_elements(t) eq 0)) then begin

       psi = read_field('psi', x, z, t, points=points, $
                        mask=mask, /equilibrium, _EXTRA=extra)

       sz = size(psi)
       if(n_elements(psi) le 1) then return, 0
   endif

   if(type eq 7) then begin ; named field
       if (strcmp(field, 'Safety Factor', /fold_case) eq 1) or $
         (strcmp(field, 'q', /fold_case) eq 1) then begin

           flux_t = flux_average('flux_t', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, _EXTRA=extra)
           
           units = ''
           name = '!6Safety Factor!X'
           symbol = '!8q!X'

           return, abs(deriv(flux, flux_t))/(2.*!pi)

       endif else $
         if(strcmp(field, 'q2', /fold_case) eq 1) then begin

           minor_r = read_field('r', x, z, t, points=points, $
                                _EXTRA=extra)
           
           r = radius_matrix(x,z,t)
           
           numvar = read_parameter('numvar', _EXTRA=extra)
           I = read_field('I', x, z, t, points=points, $
                          _EXTRA=extra)
           
           bt = sqrt(I^2/r^2)
           bp = sqrt(s_bracket(psi,psi,x,z)/r^2)
           
           field = minor_r * bt / (r * bp)
           
           units = ''
           name = '!6Safety Factor!X'
           symbol = '!8q!X'

       endif else $
         if(strcmp(field, 'flux_t', /fold_case) eq 1) then begin
           I = read_field('I', x, z, t, points=points, $
                          _EXTRA=extra)

           r = radius_matrix(x,z,t)
           field = I/r^2

           units = ''
           name = '!6Toroidal Flux!X'
           symbol = '!7w!D!8t!N!X'

           integrate = 1

       endif else $
         if(strcmp(field, 'beta_pol', /fold_case) eq 1) then begin
           I = read_field('I', x, z, t, points=points, $
                          _EXTRA=extra)
           
           bzero = read_parameter('bzero',_EXTRA=extra)
           rzero = read_parameter('rzero',_EXTRA=extra)
           izero = bzero*rzero
           print, 'izero = ', izero
       
           r = radius_matrix(x,z,t)

           bpol2 = s_bracket(psi,psi,x,z)/r^2

           ii = flux_average_field(izero^2-i^2,psi,x,z,t, r0=r0, $
             flux=flux, area=area, dV=dV, bins=bins, _EXTRA=extra)
           rr = flux_average_field(r,psi,x,z,t, r0=r0, $
             flux=flux, area=area, dV=dV, bins=bins, _EXTRA=extra)
           bb = flux_average_field(bpol2,psi,x,z,t, r0=r0, $
             flux=flux, area=area, dV=dV, bins=bins, _EXTRA=extra)
         
           symbol = '!7b!6!Dpol!N!X'
           units = ''
           name = '!6Poloidal Beta!X'

           return, 1.+0.5*ii/(bb*rr^2)
           
       endif else if (strcmp(field, 'gam', /fold_case) eq 1) then begin
           minor_r = read_field('r', x, z, t, points=points, $
                                _EXTRA=extra)
           I = read_field('I', x, z, t, points=points, $
                          _EXTRA=extra)
           P = read_field('P', x, z, t, points=points, $
                          _EXTRA=extra)
           den = read_field('den', x, z, t, points=points, $
                            _EXTRA=extra)
           gam = read_parameter('gam', _EXTRA=extra)

           r = radius_matrix(x,z,t)
           bt = sqrt(I^2/r^2)
           bp = sqrt(s_bracket(psi,psi,x,z)/r^2)
           q = minor_r * bt / (r * bp)
           
           field = 2.*gam*p/(den*r^2)*(1.+2./q^2)
           field = sqrt(field)
           name = '!6GAM Frequency!X'
           symbol = '!4x!D!6GAM!N!X'
           units = make_units(t0=-1)

       endif else $
         if(strcmp(field, 'alpha', /fold_case) eq 1) then begin
           q = flux_average('q',psi=psi,x=x,z=z, $
                            flux=flux, bins=bins, $
                            points=points, last=last, _EXTRA=extra)
           beta = flux_average('beta',psi=psi,x=x,z=z, $
                            flux=flux, bins=bins, $
                            points=points, last=last, _EXTRA=extra)
           r = flux_average('r',psi=psi,x=x,z=z, $
                            flux=flux, bins=bins,r0=r0, $
                            points=points, last=last, _EXTRA=extra)

           betap = deriv(r,beta)
         
           symbol = '!7a!X'
           units = ''
           name = '!7a!X'

           return, -q^2*betap*r0

       endif else begin
           field = read_field(field, x, z, t, points=points,$
                              symbol=symbol, units=units, _EXTRA=extra)
           name = symbol
       endelse
   endif else begin
       name = ''
       symbol = ''
       units = make_units()
   endelse

   return, flux_average_field(field, psi, x, z, t, r0=r0, flux=flux, $
                              nflux=nflux, area=area, dV=dV, bins=bins, $
                              integrate=integrate, _EXTRA=extra)
end


function flux_at_q, qval, normalized_flux=norm, $
                    q=q, flux=flux, _EXTRA=extra
   q = flux_average('q', flux=flux, nflux=nflux, /equilibrium, _EXTRA=extra)
   if(keyword_set(norm)) then flux=nflux

   dq_dpsi = deriv(flux,q)

   n = n_elements(qval)
   if(n eq 0) then return, 0

   fval = fltarr(n)
   for k=0, n-1 do begin
       ; make initial guess
       dum = min(q-qval[k],/abs,i)
       fval[k] = flux[i]
       
       ; perform newton iterations to refine result
       for j=0, 5 do begin
           dq = interpol(dq_dpsi,flux,fval[k]) 
           q0 = interpol(q,flux,fval[k])
           dpsi = (qval[k] - q0)/dq_dpsi[i]
           fval[k] = fval[k] + dpsi
           if((fval[k] gt max(flux)) or (fval[k] lt min(flux))) then begin
               print, 'flux_at_q: could find surface with q = ', $
                 qval[k]
               break
           endif
;           print, qval[k], q0, fval[k]
       end
   end

   return, fval
end


pro plot_field, name, time, x, y, points=p, mesh=plotmesh, $
                mcolor=mc, lcfs=lcfs, title=title, units=units, $
                maskrange=maskrange, maskfield=maskfield, range=range, $
                rrange=rrange, zrange=zrange, linear=linear, $
                xlim=xlim, cutx=cutx, cutz=cutz, mpeg=mpeg, $
                mask_val=mask_val, boundary=boundary, q_contours=q_contours, $
                _EXTRA=ex

   if(n_elements(time) eq 0) then time = 0
   if(n_elements(p) eq 0) then p = 50
   if(n_elements(title) eq 0) then notitle = 1 else notitle = 0

   if(size(name, /type) eq 7) then begin
       print, 'points = ', p
       field = read_field(name, x, y, t, slices=time, mesh=mesh, $
                          points=p, rrange=rrange, zrange=zrange, $
                          symbol=fieldname, units=u, linear=linear, $
                          mask=mask, _EXTRA=ex)
       if(n_elements(field) le 1) then return

       if(n_elements(units) eq 0) then units=u
   endif else begin
       field = name
       if(n_elements(field) le 1) then return
   endelse

   ; remove NaN's from result
   i = where(not float(finite(field)), count)
   if(count gt 0) then begin
       print, "removing NaN's"
       field[i] = 0.
   endif

   if(n_elements(maskrange) eq 2) then begin
       if((strcmp(name, maskfield) eq 1) and $
         (not keyword_set(linear))) then begin
           psi = field
       endif else begin
           psi = read_field(maskfield, slices=time, mesh=mesh, $
                            points=p, rrange=rrange, zrange=zrange, $
                           /equilibrium, mask=mask, _EXTRA=ex)
       endelse
       newmask = (psi ge maskrange[0]) and (psi le maskrange[1])
       field = newmask*field + (1-newmask)*(min(field-newmask*field,/absolute))
   endif

   sz = size(field, /dimension)
   nt = sz[0]

   ; open mpeg object
   if(n_elements(mpeg) ne 0) then begin
       mpegid = mpeg_open([640,480],bitrate=104857200, iframe_gap=4)
   end

   if(n_elements(range) eq 0) then range = [min(field),max(field)]

   if(n_elements(mask_val) ne 0) then begin
       for k=0, nt-1 do begin
           field[k,*,*] = field[k,*,*] - mask*(field[k,*,*] - mask_val)
       end
   endif

   if(n_elements(cutx) gt 0) then begin
       dum = min(x-cutx,i,/absolute)
       plot, y, field[0,i,*], _EXTRA=ex
   endif else if(n_elements(cutz) gt 0) then begin
       dum = min(y-cutz,i,/absolute)
       plot, x, field[0,*,i], _EXTRA=ex
   endif else begin
       for k=0, nt-1 do begin
           if((notitle eq 1) and (n_elements(t) ne 0)) then begin
               title = fieldname
           end
           
           contour_and_legend, field[k,*,*], x, y, title=title, $
             label=units, $
             xtitle=make_label('!8R!X', /l0, _EXTRA=ex), $
             ytitle=make_label('!8Z!X', /l0, _EXTRA=ex), $
             range=range, _EXTRA=ex

           if(keyword_set(lcfs)) then $
             plot_lcfs, points=p, slice=time, _EXTRA=ex

           if(n_elements(q_contours) ne 0) then begin
               fval = flux_at_q(q_contours,_EXTRA=ex)
               plot_flux_contour, fval, closed=0, /overplot, _EXTRA=ex
           endif

           if(keyword_set(boundary)) then plotmesh=1
           if(keyword_set(plotmesh)) then begin
               loadct, 12
               plot_mesh, mesh=mesh, color=color(3,5), /oplot, $
                 boundary=boundary, _EXTRA=ex
           endif
           
           if(n_elements(mpeg) ne 0) then begin
               image = tvrd(true=1)
               
               image[0,*,*] = rotate(reform(image[0,*,*]), 7)
               image[1,*,*] = rotate(reform(image[1,*,*]), 7)
               image[2,*,*] = rotate(reform(image[2,*,*]), 7)
               
               mpeg_put, mpegid, image=image, frame=5*k
           end
       end
   endelse

   if(n_elements(mpeg) ne 0) then begin
       print, 'Writing mpeg...'
       mpeg_save, mpegid, filename=mpeg
       mpeg_close, mpegid
   end
end



;======================================================
; plot_flux_average
; ~~~~~~~~~~~~~~~~~
;
; plots the flux average quantity "name" at a give time
;======================================================
pro plot_flux_average, field, time, filename=filename, $
                       color=c, names=names, bins=bins, linear=linear, $
                       xlog=xlog, ylog=ylog, overplot=overplot, $
                       lcfs=lcfs, normalized_flux=norm, points=pts, $
                       minor_radius=minor_radius, smooth=sm, t=t, rms=rms, $
                       bw=bw, srnorm=srnorm, last=last, mks=mks, cgs=cgs, $
                       q_contours=q_contours, _EXTRA=extra

   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(n_elements(time) eq 0) then time=0
   if(keyword_set(last)) then $
     time = fix(read_parameter('ntime',filename=filename)-1)

   nfiles = n_elements(filename)
   if(nfiles gt 1) then begin
       if(n_elements(names) eq 0) then names=filename
       if(keyword_set(bw)) then begin
           ls = indgen(nfiles)
           colors = replicate(color(0,1), nfiles)
       endif else begin
           print, 'hello'
           if(n_elements(c) eq 0) then colors = colors(nfiles)
           ls = replicate(0,nfiles)
       endelse
       if(n_elements(time) eq 1) then time = replicate(time,nfiles)

       for i=0, nfiles-1 do begin
           newfield = field
           plot_flux_average, newfield, time[i], filename=filename[i], $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, $
             linear=linear
       end
       if(n_elements(names) gt 0) then begin
           plot_legend, names, color=colors, ylog=ylog, xlog=xlog, $
             linestyle=ls, _EXTRA=extra
       endif    
       
       return
   endif

   nt = n_elements(time)
   if(nt gt 1) then begin
       if(n_elements(names) eq 0) then names=strarr(nt)
       if(keyword_set(bw)) then begin
           ls = indgen(nt)
           colors = replicate(color(0,1), nt)
       endif else begin
           if(n_elements(c) eq 0) then colors = colors(nt)
           ls = replicate(0,nt)
       endelse
       
       for i=0, n_elements(time)-1 do begin
           newfield = field
           plot_flux_average, newfield, time[i], filename=filename, $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             t=t, rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, $
             linear=linear
           names[i] = string(format='(%"!8t!6 = %d !7s!D!8A!N!X")', t)
       end

       if(n_elements(names) gt 0) then begin
           plot_legend, names, color=colors, ylog=ylog, xlog=xlog, $
             linestyle=ls, _EXTRA=extra
       endif

       return
   endif

   xtitle='!7w!X'
   title = ''

   fa = flux_average(field,slice=time,flux=flux,points=pts,filename=filename, $
                     name=title, symbol=symbol, units=units, bins=bins, $
                     psi=psi,x=x,z=z,t=t,nflux=nflux,linear=linear, $
                     mks=mks, cgs=cgs, EXTRA=extra)

   if(n_elements(fa) le 1) then begin
       print, 'Error in flux_average. returning.'
       return
   endif

   ytitle = symbol
   if(strlen(units) gt 0) then ytitle = ytitle + '!6 ('+units+ '!6)!X'

   if(keyword_set(rms)) then begin
       fa2 = flux_average(field^2,flux=flux,nflux=nflux,points=pts,slice=time,$
                          filename=filename, t=t, linear=linear,_EXTRA=extra, $
                         mks=mks, cgs=cgs)
       fa = sqrt(1. - fa^2/fa2)

       ytitle = '!9S!6(1 - !12<' + symbol + '!12>!U2!n/!12<' + $
         symbol + '!6!U2!N!12>!6)!X'
       title = '!6Poloidal Deviation of ' + title + '!X'
   end


;   if(t gt 0) then begin
;       title = "!12<" + title + $
;        string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!N!6)!12>!7!Dw!N!X")',t)
;   endif else begin
;       title = "!12<!8" + title + $
;         string(FORMAT='("!6(!8t!6 = ",G0,")!3!12>!7!Dw!N!X")', t)
;   endelse
;   title = "!12<" + title + "!12>!7!Dw!N!X"

   if(keyword_set(norm)) then begin
       flux = nflux
       xtitle = '!7W!X'
       lcfs_psi = 1.
   endif
   if(keyword_set(srnorm)) then begin
       flux = sqrt(nflux)
       xtitle = '!9r!7W!X'
       lcfs_psi = 1.
   end       

   if(keyword_set(minor_radius)) then begin
       flux = flux_average('r',points=pts,file=filename,t=t,linear=linear,$
                    name=xtitle,bins=bins,units=units,slice=time,_EXTRA=extra,$
                          mks=mks, cgs=cgs)
       xtitle = '!12<' + xtitle + '!12> !6 ('+units+')!X'
   endif

   if(n_elements(sm) eq 1) then begin
       fa = smooth(fa,sm)
   end

   if(keyword_set(overplot)) then begin
       oplot, flux, fa, color=c, _EXTRA=extra
   endif else begin
       plot, flux, fa, xtitle=xtitle, $
         ytitle=ytitle, title=title, xlog=xlog, ylog=ylog, $
         _EXTRA=extra
   endelse

   if(keyword_set(lcfs)) then begin
       oplot, [lcfs_psi,lcfs_psi], !y.crange, linestyle=2, color=c
   endif

   if(n_elements(q_contours) ne 0) then begin
       fvals = flux_at_q(q_contours, points=pts, filename=filename, $
                         slice=time, normalized_flux=norm, bins=bins)
       for k=0, n_elements(fvals)-1 do begin
           oplot, [fvals[k], fvals[k]], !y.crange, linestyle=1
       end
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


pro animate, name, nslices=nslices, slice=slice, _EXTRA=extra
   ntor = read_parameter('ntor',_EXTRA=extra)
  
   print, 'ntor = ', ntor

   if(n_elements(slice) eq 0) then begin
       field = read_field(name,x,z,t,_EXTRA=extra, /last,$
                          symbol=symbol,units=units)
       field_i = read_field(name+'_i',x,z,t,_EXTRA=extra, /last)
   endif else begin
       field = read_field(name,x,z,t,_EXTRA=extra, slice=slice,$
                          symbol=symbol,units=units)
       field_i = read_field(name+'_i',x,z,t,_EXTRA=extra, slice=slice)
   endelse

   if(n_elements(field) le 1) then return

   if(n_elements(nslices) eq 0) then nslices = 12
      
   f = fltarr(nslices,n_elements(x),n_elements(z))

   for i=0, nslices-1 do begin
       f[i,*,*] = field[0,*,*]*cos(2.*!pi*i/nslices) + $
         field_i[0,*,*]*sin(2.*!pi*i/nslices)
   end

   fstr = '("!9P!6 = ",g0.2,"!7p!X")'
   
   range=[min(f), max(f)]

   for i=0, nslices-1 do begin
       if(i eq 0) then begin
           title = '!9P!6 = 0!X'
       endif else if(2.*i/nslices eq 1) then begin
           title = '!9P!6 = !7p!X'
       endif else begin
           title = string(format=fstr,2.*i/nslices)
       endelse
       title = symbol + ' !6(' + title + '!6)!X'
       print, 'plotting...'
       plot_field, f[i,*,*],-1,x,z, title=title, range=range, $
        _EXTRA=extra
       print, 'reading...'
       image = tvrd(/true)
       
       print, 'setting...'
       if(i eq 0) then begin
           sz = size(image, /dim)
           xinteranimate, set=[sz[1], sz[2],nslices], $
             mpeg_bitrate=104857200
; mpeg_quality=100, mpeg_format=1, mpeg_bitrate=429496729200
       endif
       xinteranimate, frame=i, image=image
   endfor

   xinteranimate, 10
end



pro write_geqdsk, eqfile=eqfile, b0=b0, l0=l0, $
                  psilim=psilim, _EXTRA=extra
  
  if(n_elements(slice) eq 0) then begin
      slice = read_parameter('ntime', _EXTRA=extra) - 1
  end
  if(n_elements(b0) eq 0) then b0 = 10000.
  if(n_elements(l0) eq 0) then l0 = 100.
  if(n_elements(eqfile) eq 0) then eqfile = 'geqdsk.out'

  bzero = read_parameter('bzero', _EXTRA=extra)
  rzero = read_parameter('rzero', _EXTRA=extra)

  ; calculate flux averages
  psi = read_field('psi',x,z,t,mesh=mesh,/equilibrium,_EXTRA=extra)
  jphi= read_field('jphi',x,z,t,/equilibrium,_EXTRA=extra)
  p0 = read_field('p',x,z,t,/equilibrium,_EXTRA=extra)
  I0 = read_field('I',x,z,t,/equilibrium,_EXTRA=extra)
  r = radius_matrix(x,z,t)
  beta = r^2*2.*p0/(s_bracket(psi,psi,x,z) + I0^2)
  beta0 = mean(2.*p0*r^2/(bzero*rzero)^2)
  b2 = (s_bracket(psi,psi,x,z) + I0^2)/r^2
  dx = mean(deriv(x))
  dz = mean(deriv(z))
  tcur = abs(total(jphi*dx*dz/r))
  print, 'current = ', tcur

  ; calculate magnetic axis and xpoint
  lcfs_psi = lcfs(psi,x,z, axis=axis, xpoint=xpoint, $
                  flux0=flux0, _EXTRA=extra)

  ; plot psi
  contour_and_legend, psi, x, z

  ; calculate wall points
  bound_xy = get_boundary_path(mesh=mesh, _EXTRA=extra)
  nwall = n_elements(bound_xy[0,*])
  rwall = fltarr(nwall)
  zwall = fltarr(nwall)
  rwall = bound_xy[0,*]
  zwall = bound_xy[1,*]

  ifixedb = read_parameter('ifixedb', _EXTRA=extra)
  print, 'ifixedb = ', ifixedb

  if(ifixedb eq 1) then begin
      psilim = 0.
      lcfs_psi = 0.
      
      ; for ifixedb eq 1, boundary points are same as wall points
      nlim = nwall
      rlim = rwall
      zlim = zwall

  endif else begin
      ; find points at psi = psilim to use as boundary points
      if(n_elements(psilim) eq 0) then begin
          ; use psilim = lcfs if psilim is not given
          psilim = lcfs_psi
      endif
      print, 'lcfs_psi = ', lcfs_psi
      print, 'psilim = ', psilim
      
      window, 0
      lcfs_xy = path_at_flux(psi,x,z,t,psilim)

      ; count only points on separatrix above the xpoint
      if(n_elements(xpoint) gt 1) then begin
          if(xpoint[0] ne 0 or xpoint[1] ne 0) then begin
              if(xpoint[1] lt axis[1]) then begin
                  lcfs_mask = lcfs_xy[1,*] gt xpoint[1]
              endif else begin
                  lcfs_mask = lcfs_xy[1,*] lt xpoint[1]
              endelse
          endif else begin
              lcfs_mask = fltarr(1, n_elements(lcfs_xy[1,*]))
              lcfs_mask[*] = 1
          endelse
      endif else begin
          lcfs_mask = fltarr(1, n_elements(lcfs_xy[1,*]))
          lcfs_mask[*] = 1
      endelse
      
      oplot, lcfs_xy[0,*], lcfs_xy[1,*]
      
      nlim = fix(total(lcfs_mask))
      rlim = fltarr(nlim)
      zlim = fltarr(nlim)
      j = 0 
      for i=0, n_elements(lcfs_mask)-1 do begin
          if(not lcfs_mask[0, i]) then continue
          rlim[j] = lcfs_xy[0,i]
          zlim[j] = lcfs_xy[1,i]
          j = j+1
      end
      print, 'lim points = ', nlim, j
  endelse

  ; reduce the boundary points if necessary
  while(nlim ge 500) do begin
      print, 'reducing lim points...'
      nlim = nlim / 2
      new_rlim = fltarr(nlim)
      new_zlim = fltarr(nlim)
      for k=0, nlim-1 do begin
          new_rlim[k] = rlim[2*k]
          new_zlim[k] = zlim[2*k]
      endfor
      rlim = new_rlim
      zlim = new_zlim
      print, 'new lim points = ', nlim
  end
      
  oplot, rlim, zlim, psym=4
 

  ; calculate flux averages
  p = flux_average_field(p0,psi,x,z,flux=flux,$
                   limiter=limiter,_EXTRA=extra)
  pp = s_bracket(p0,psi,x,z)/s_bracket(psi,psi,x,z)
  pprime = flux_average_field(pp,psi,x,z,flux=flux,$
                 limiter=limiter,_EXTRA=extra)
  I = flux_average_field(I0,psi,x,z,flux=flux,$
                   limiter=limiter,_EXTRA=extra)
  ffp = I0*s_bracket(I0,psi,x,z)/s_bracket(psi,psi,x,z)
  ffprim = flux_average_field(ffp,psi,x,z,flux=flux,$
                   limiter=limiter,_EXTRA=extra)
  q = flux_average('q',slice=time,psi=psi,x=x,z=z,t=t,flux=flux, $
                   limiter=limiter,_EXTRA=extra)
;  q = smooth(q,5,/edge)
  jb = (s_bracket(I0,psi,x,z) - jphi*I0)/r^2
  jdotb = flux_average_field(jb,psi,x,z,t,flux=flux,$
                   limiter=limiter,_EXTRA=extra)
  r2i = flux_average_field(1./r^2,psi,x,z,flux=flux,$
                   limiter=limiter,_EXTRA=extra)

  betacent = field_at_point(beta[0,*,*], x, z, axis[0], axis[1])

  ; to cgs............. to si
  c = 3e10
  p = p*b0^2/(4.*!pi)           / 10.
  p0 = p0*b0^2/(4.*!pi)         / 10.
  psi = psi*b0*l0^2             / 1e8
  flux = flux*b0*l0^2           / 1e8
  flux0 = flux0*b0*l0^2         / 1e8
  psilim=psilim*b0*l0^2         / 1e8
  pprime = pprime*b0^2/(4.*!pi) / 10.
  pprime = pprime/(b0*l0^2)     * 1e8
  I = I*b0*l0                   / (1e4*100.)
  ffprim = ffprim*(b0*l0)^2     / (1e4*100.)^2
  ffprim = ffprim/(b0*l0^2)     * 1e8
  bzero = bzero*b0              / 1e4
  b2 = b2*b0^2                  / (1e4)^2
  tcur = tcur*b0*c*l0/(4.*!pi)  / 3e9
  r = r*l0                      / 100.
  x = x*l0                      / 100.
  z = z*l0                      / 100.
  axis[0] = axis[0]*l0          / 100.
  axis[1] = axis[1]*l0          / 100.
  rzero = rzero*l0              / 100.
  jdotb = jdotb*b0^2*c/(l0*4.*!pi) / (1e4*3e5)
  r2i = r2i/(l0^2)              * 100.^2

  nr = n_elements(flux)
  print, 'nr = ', nr

  psi = -psi
  I = -I
  psilim = -psilim
  psimin = -flux0
  flux = -flux
  pprime = -pprime
  ffprim = -ffprim

  nflux = (flux - min(flux))/(max(flux)-min(flux))

  name = ['name0001', 'name0002', 'name0003', $
          'name0004', 'name0005', 'name0006']
  idum = 13
  nr = n_elements(flux)
  nz = n_elements(z)
  rdim = max(x) - min(x)
  zdim = max(z) - min(z)
  xplas = 1.
  ccon = min(x)
  zmid = (max(z) + min(z))/2.
  rmag = axis[0]
  zmag = axis[1]
  zip = tcur
  bcentr = bzero*rzero/rmag
  beta0 = beta0
  beta_n = 100.*(bzero*rzero/rmag)*beta0/(zip/1e6)
  xdum = 0.

  print, 'rdim, zdim', rdim, zdim
  print, 'rmag, zmag', rmag, zmag
  print, 'bcentr = ', bcentr
  print, 'beta0 = ', beta0
  print, 'betacent = ', betacent
  print, 'beta_n = ', beta_n
  print, 'zip = ', zip
  print, 'psimin, psilim = ', psimin, psilim
  print, 'min, max (flux) = ', min(flux), max(flux)


  ; set up formatting codes
  f2000 = '(6A8,3I4)'
  f2020 = '(5E16.9)'
  f2022 = '(2i5)' 


  ; output to eqdsk file
  print, 'outputting to eqdsk format...'
  file = 1
  
  openw, file, eqfile

  printf, file, format=f2000, name, idum, nr, nz
  printf, file, format=f2020, rdim, zdim, xplas, ccon, zmid
  printf, file, format=f2020, rmag, zmag, psimin, psilim, bcentr
  printf, file, format=f2020, zip, psimin, beta0, rmag, betacent
  printf, file, format=f2020, zmag, beta_n, psilim, xdum, xdum

  printf, file, format=f2020, I
  printf, file, format=f2020, p
  printf, file, format=f2020, ffprim
  printf, file, format=f2020, pprime
  printf, file, format=f2020, psi[0,*,*]
  printf, file, format=f2020, q
  printf, file, format=f2022, nlim, nwall
  printf, file, format=f2020, transpose([[rlim],[zlim]])
  printf, file, format=f2020, transpose([[rwall],[zwall]])

  close, file

  ; output to jsolver
  print, 'outputting to jsolver format...'

  jsfile = 'jsfile'
  openw, file, jsfile  

  ncycle=  1
  isyms=  0
  ipest=  1
  kmax=nlim-1
  npsit = n_elements(flux)-1 ; remove final point to avoid NaN's
  
  times=  0.1140E-01
  xaxes=  axis[0]
  zmags=  axis[1]
  apls=   0.3656E+05
  betas=  0.5184E-02
  betaps= 0.1472E+00
  ali2s=  0.5458E+00
  qsaws=  0.5000E+00
  psimins=min(flux)
  psilims=max(flux)

  f2201 = '(20x,10a8)'
  f6100 = '(5i10)'
  f6101 = '(5e20.12)'

  ;  gzeros= R * B_T in m-T   (at vacuum)
  gzeros = bzero*rzero

  printf, file, format=f6100, ncycle,isyms,ipest,npsit,kmax
  printf, file, format=f6101, times,xaxes,zmags,gzeros,apls,betas,betaps, $
    ali2s,qsaws,psimins,psilims

  mu0 = (4.*!pi*1.e-7)
  ajpest2 = jdotb/(I*r2i)

  printf, file, format=f6101, mu0*p[0:npsit-1]
  printf, file, format=f6101, mu0*pprime[0:npsit-1]
  printf, file, format=f6101, -mu0*ajpest2[0:npsit-1]
  printf, file, format=f6101, flux[0:npsit-1] - flux[0]
  printf, file, format=f6101, rlim
  printf, file, format=f6101, zlim
  
  close, file

  window, 1
;  contour_and_legend, psi,x,z, /iso
;  loadct,12
;  oplot, rlim, zlim, color=color(1,3), thick=3.0
;  oplot, rwall, zwall, color=color(2,3), thick=3.0
 
;  plot, flux, pprime

  loadct,12
  !p.multi = [0,3,2]
  plot, nflux, mu0*p, title='p'
  plot, nflux, mu0*pprime, title="p'"
  oplot, nflux, mu0*deriv(flux,p), color=color(1,2), linestyle=2
  plot, nflux, I, title='f'
  plot, nflux, ffprim, title="f f'"
  oplot, nflux, I*deriv(flux,I), color=color(1,2), linestyle=2
  plot, nflux, q, title='q', yrange=[0,10]
  plot, nflux, mu0*ajpest2, title='ajpest2'

 
  !p.multi=0

end

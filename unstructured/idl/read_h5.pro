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
    s1 = size(a, /dim)
    s2 = size(b, /dim)

    if(s1[0] gt 1 and s2[0] eq 1) then begin
        c = fltarr(s1[0], n_elements(x), n_elements(z))
        for i=0, s1[0]-1 do c[i,*,*] = a_bracket(a[i,*,*], b[0,*,*],x,z)
        return, c
    endif else if(s2[0] gt 1 and s1[0] eq 1) then begin
        c = fltarr(s2[0], n_elements(x), n_elements(z))
        for i=0, s2[0]-1 do c[i,*,*] = a_bracket(a[0,*,*], b[i,*,*],x,z)
        return, c
    endif else if(s1[0] eq 1 and s2[0] eq 1) then begin
        return, -dx(a,x)*dz(b,z) + dz(a,z)*dx(b,x)
    endif else begin
        print, 'Error: sizes do not conform!'
        return, 0
    endelse
end

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

function clamp_and_shift, vec, shift=n
  ; clamp
  new = vec
  n = 0
  for i=0, n_elements(vec)-1 do begin
      if(new[i] lt -!pi) then new[i] = new[i] + 2.*!pi
      if(new[i] ge !pi) then new[i] = new[i] - 2.*!pi
      if(i gt 0) then begin
          if(abs(last - new[i]) gt !pi) then n = i
      endif
      last = new[i]
  end
  new = shift(new, -n)
  return, new
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
;   on_error, 2
   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(hdf5_file_test(filename) eq 0) then return, 0

;     catch, Error_status
;     if Error_status ne 0 then begin
;         print, 'Error reading: ', filename
;         catch, /cancel
;         return, 0
;     end

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

   n = n_elements(slice)

   t = fltarr(n)

   for i=0, n-1 do begin
       file_id = h5f_open(filename)
       time_group_id = h5g_open(file_id, time_name(slice[i]))
       time_id = h5a_open_name(time_group_id, "time")

       t[i] = h5a_read(time_id)

       h5a_close, time_id
       h5g_close, time_group_id
       h5f_close, file_id
   end

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


pro get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, $
                        zeff=zeff, ion_mass=ion_mass, _EXTRA=extra
   b0_norm = read_parameter('b0_norm', _EXTRA=extra)
   n0_norm = read_parameter('n0_norm', _EXTRA=extra)
   l0_norm = read_parameter('l0_norm', _EXTRA=extra)
   zeff = read_parameter("zeff", _EXTRA=extra)
   ion_mass = read_parameter("ion_mass", _EXTRA=extra)
   i = where(zeff eq 0, count)
   if(count gt 0) then zeff[i] = 1.
   i = where(ion_mass eq 0, count)
   if(count gt 0) then ion_mass[i] = 1.
end

;===============================================================
; convert_units
; ~~~~~~~~~~~~~
;
; converts x having dimensions d to cgs units
; where b0, n0, and l0 are the normalizations (in cgs units)
;===============================================================
pro convert_units, x, d, b0, n0, l0, zeff, mi, cgs=cgs, mks=mks
   if(n_elements(x) eq 0) then return

   if(not (keyword_set(cgs) or keyword_set(mks))) then return

   if(b0 eq 0 or n0 eq 0 or l0 eq 0 or zeff eq 0 or mi eq 0) then begin
       print, "Warning: unknown conversion factors."
       print, "Using l0=100, B0=1e4, n0=1e14."
       l0 = 100.
       b0 = 1.e4
       n0 = 1.e14
       zeff = 1.
       mi = 1.
   endif

   val = 1.
   if(keyword_set(cgs)) then begin
       fp = (4.*!pi)
       c0 = 3.e10
       v0 = 2.18e11*b0/sqrt(mi*n0)
       t0 = l0/v0
       temp0 = b0^2/(fp*n0) * 1./(1.6022e-12)
       i0 = c0*b0*l0/fp
       e0 = b0^2*l0^3/fp

       val = fp^d[0] $
         * c0^d[1] $
         * n0^d[2] $
         * v0^d[3] $
         * b0^d[4] $
         * t0^d[8] $
         * l0^d[9] $
         * temp0^d[5] $
         * i0^d[6] $
         * e0^d[7]
       
   endif else if(keyword_set(mks)) then begin
       convert_units, x, d, b0, n0, l0, zeff, mi, /cgs

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
  i0 =    [0,0,0,0,0,0,1,0,0,0]
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
  if(keyword_set(j))      then d = d + i0*j - 2*l0
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
; [4pi, c, n0, vA0, B0, T0, i0, e0, tA0, L0]
; [  0, 1,  2,   3,  4,  5,  6,  7,   8,  9]
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
       x[9] = x[9]          + x[6] - 3*x[7]
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

   if(strcmp(name, 'psi', /fold_case) eq 1 or $
      strcmp(name, 'psi_i', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=1+itor)
       return, "!7w!X"
   endif else if(strcmp(name, 'I', /fold_case) eq 1 or $
                 strcmp(name, 'I_i', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=itor)
       return, "!8I!X"
   endif else if(strcmp(name, 'phi', /fold_case) eq 1 or $
                 strcmp(name, 'phi_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=1+itor)
       return, "!8U!X"
   endif else if(strcmp(name, 'V', /fold_case) eq 1 or $
                 strcmp(name, 'V_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=itor)
       return, "!8V!X"
   endif else if(strcmp(name, 'chi', /fold_case) eq 1 or $
                 strcmp(name, 'chi_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=1)
       return, "!7v!X"
   endif else if(strcmp(name, 'eta', /fold_case) eq 1 or $
                 strcmp(name, 'eta_i', /fold_case) eq 1) then begin
       units = dimensions(/eta)
       return, "!7g!X"
   endif else if(strcmp(name, 'den', /fold_case) eq 1 or $
                 strcmp(name, 'den_i', /fold_case) eq 1) then begin
       units = dimensions(/n0)
       return, "!8n!Di!N!X"
   endif else if(strcmp(name, 'p', /fold_case) eq 1 or $
                 strcmp(name, 'p_i', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!8p!X"
   endif else if(strcmp(name, 'pe', /fold_case) eq 1 or $
                 strcmp(name, 'pe_i', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!8p!De!N!X"
   endif else if(strcmp(name, 'te', /fold_case) eq 1 or $
                 strcmp(name, 'te_i', /fold_case) eq 1) then begin
       units = dimensions(/temperature)
       return, "!8T!De!N!X"
   endif else if(strcmp(name, 'ti', /fold_case) eq 1 or $
                 strcmp(name, 'ti_i', /fold_case) eq 1) then begin
       units = dimensions(/temperature)
       return, "!8T!Di!N!X"
   endif else if(strcmp(name, 'sigma', /fold_case) eq 1 or $
                 strcmp(name, 'sigma_i', /fold_case) eq 1) then begin
       units = dimensions(/n0,t0=-1)
       return, "!7r!X"
   endif else if(strcmp(name, 'force_phi', /fold_case) eq 1 or $
                 strcmp(name, 'force_phi_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, l0=-1)
       return, "!8F!D!9p!N!X"
   endif else if(strcmp(name, 'pforce', /fold_case) eq 1 or $
                 strcmp(name, 'pforce_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, l0=-1)
       return, "!8F!D!9p!N!X"
   endif else if(strcmp(name, 'heat_source', /fold_case) eq 1 or $
                 strcmp(name, 'heat_source_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, t0=-1)
       return, "!8Q!X"
   endif else if(strcmp(name, 'kappa', /fold_case) eq 1 or $
                 strcmp(name, 'kappa_i', /fold_case) eq 1) then begin
       units = dimensions(/n0, l0=2, t0=-1)
       return, "!7j!X"
   endif else if((strcmp(name, 'visc', /fold_case) eq 1) or $
     (strcmp(name, 'visc_c', /fold_case) eq 1) or $
                 strcmp(name, 'visc_i', /fold_case) eq 1 or $
                 strcmp(name, 'visc_c_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, /t0)
       return, "!7l!X"
   endif else if(strcmp(name, 'jphi', /fold_case) eq 1  or $
                 strcmp(name, 'jphi_i', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=itor-1)
       return, "!7D!6!U*!N!7w!X"
   endif else if(strcmp(name, 'vor', /fold_case) eq 1 or $
                 strcmp(name, 'vor_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=itor-1)
       return, "!7D!6!U*!N!8U!X"
   endif else if(strcmp(name, 'com', /fold_case) eq 1 or $
                 strcmp(name, 'com_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=-1)
       return, "!9G.!17v!X"
   endif else if(strcmp(name, 'torque_em', /fold_case) eq 1  or $
                 strcmp(name, 'torque_em_i', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!7s!D!8EM!N!X"
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
       xtitle = make_label('!8R!X',/l0,_EXTRA=ex)
       ytitle = make_label('!8Z!X',/l0,_EXTRA=ex)
       plot, mesh.elements._data[4,*], xtitle=xtitle, ytitle=ytitle, $
         mesh.elements._data[5,*], psym = 3, _EXTRA=ex, /nodata
   endif  

   get_normalizations, b0=b0, n0=n0, l0=l0, zeff=zeff, ion_mass=mi, _EXTRA=ex
   fac = 1.
   convert_units, fac, dimensions(/l0), b0, n0, l0, zeff, mi, _EXTRA=ex

   ct3
   col = color(1,10)
 
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

   maxr = [max(mesh.elements._data[4,*]), max(mesh.elements._data[5,*])]*fac
   minr = [min(mesh.elements._data[4,*]), min(mesh.elements._data[5,*])]*fac

   for i=long(0), nelms-1 do begin
       a = mesh.elements._data[0,i]*fac
       b = mesh.elements._data[1,i]*fac
       c = mesh.elements._data[2,i]*fac
       t = mesh.elements._data[3,i]
       x = mesh.elements._data[4,i]*fac
       y = mesh.elements._data[5,i]*fac
       bound = fix(mesh.elements._data[6,i])

       p1 = [x, y]
       p2 = p1 + [(b+a) * cos(t), (b+a) * sin(t)]
       p3 = p1 + [b * cos(t) - c * sin(t), $
                  b * sin(t) + c * cos(t)]
       delta = 0.0
       q1 = (1.-2.*delta)*p1 + delta*p2 + delta*p3
       q2 = delta*p1 + (1.-2.*delta)*p2 + delta*p3
       q3 = delta*p1 + delta*p2 + (1.-2.*delta)*p3
      
       if(boundary) then pp=bound else pp=7
 
       if((pp and 1) eq 1) then begin
           if((bound and 1) eq 1) then begin
               izone = (bound and 120)/2^3 + 1
               c = color(izone+1)
               oplot, [q1[0],q2[0]]+xzero, [q1[1],q2[1]]+zzero, $
                 color=c, thick=!p.thick*3
           end else begin
               oplot, [p1[0],p2[0]]+xzero, [p1[1],p2[1]]+zzero, $
                 color=col, thick=!p.thick/2.
           end
       end
       if((pp and 2) eq 2) then begin
           if((bound and 2) eq 2) then begin
               izone = (bound and 1920)/2^7 + 1
               c = color(izone+1)
               oplot, [q2[0],q3[0]]+xzero, [q2[1],q3[1]]+zzero, $
                 color=c, thick=!p.thick*3
           end else begin
               oplot, [p2[0],p3[0]]+xzero, [p2[1],p3[1]]+zzero, $
                 color=col, thick=!p.thick/2.
           end

       end
       if((pp and 4) eq 4) then begin
           if((bound and 4) eq 4) then begin
               izone = (bound and 30720)/2^11 + 1
               c = color(izone+1)
               oplot, [q3[0],q1[0]]+xzero, [q3[1],q1[1]]+zzero, $
                 color=c, thick=!p.thick*3
           end else begin
               oplot, [p3[0],p1[0]]+xzero, [p3[1],p1[1]]+zzero, $
                 color=col, thick=!p.thick/2.
           end
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

   small = (a+b+c)*1e-4

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

   sz = size(field, /dim)
   if(sz[0] eq 80) then begin
       threed = 1
   endif else begin
       threed = 0
   endelse
   
   mi = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0]
   ni = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5]
   sum = 0.

   if(n_elements(op) eq 0) then op = 1

   co = cos(theta)
   sn = sin(theta)

   op2 = (op-1) / 10
   op1 = op - op2*10

   for p=0, 19 do begin
       temp = 0.
        case op1 of
        1: temp = localpos[0]^mi[p]*localpos[1]^ni[p]
        2: begin
            if(mi[p] ge 1) then $
              temp = temp $
              + co*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^ni[p]
            if(ni[p] ge 1) then $
              temp = temp $
              - sn*ni[p]*localpos[0]^mi[p]*localpos[1]^(ni[p]-1)
           end
        3: begin
            if(mi[p] ge 1) then $
              temp = temp $
              + sn*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^ni[p]
            if(ni[p] ge 1) then $
              temp = temp $
              + co*ni[p]*localpos[0]^mi[p]*localpos[1]^(ni[p]-1)
           end
        4: begin
            if(mi[p] ge 2) then $
              temp = temp $
              + co*co*mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
            if(ni[p] ge 2) then $
              temp = temp $
              + sn*sn*ni[p]*(ni[p]-1)*localpos[0]^mi[p]*localpos[1]^(ni[p]-2)
            if(mi[p] ge 1 and ni[p] ge 1) then $
              temp = temp $
              -2.*co*sn*ni[p]*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^(ni[p]-1)
           end
        6: begin
            if(mi[p] ge 2) then $
              temp = temp $
              + sn*sn*mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
            if(ni[p] ge 2) then $
              temp = temp $
              + co*co*ni[p]*(ni[p]-1)*localpos[0]^mi[p]*localpos[1]^(ni[p]-2)
            if(mi[p] ge 1 and ni[p] ge 1) then $
              temp = temp $
              +2.*co*sn*ni[p]*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^(ni[p]-1)
           end
        7: begin
            if(mi[p] ge 2) then $
              temp = temp + $
              mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
            if(ni[p] ge 2) then $
              temp = temp + $
              ni[p]*(ni[p]-1)*localpos[1]^(ni[p]-2)*localpos[0]^mi[p]
           end
        8: begin
            temp = temp + $
             ( co $
              *mi[p]*(mi[p]-1)*(mi[p]-2)*localpos[0]^(mi[p]-3>0) $
              *                          localpos[1]^ ni[p]      $
             - sn $
              *                          localpos[0]^ mi[p]      $
              *ni[p]*(ni[p]-1)*(ni[p]-2)*localpos[1]^(ni[p]-3>0) $
             - sn $
              *mi[p]*(mi[p]-1)*          localpos[0]^(mi[p]-2>0) $
              *ni[p]*                    localpos[1]^(ni[p]-1>0) $
             + co $
              *mi[p]*                    localpos[0]^(mi[p]-1>0) $
              *ni[p]*(ni[p]-1)*          localpos[1]^(ni[p]-2>0))
           end
        9: begin
            temp = temp + $
             ( sn $
              *mi[p]*(mi[p]-1)*(mi[p]-2)*localpos[0]^(mi[p]-3>0) $
              *                          localpos[1]^ ni[p]      $
             + co $
              *                          localpos[0]^ mi[p]      $
              *ni[p]*(ni[p]-1)*(ni[p]-2)*localpos[1]^(ni[p]-3>0) $
             + co $
              *mi[p]*(mi[p]-1)*          localpos[0]^(mi[p]-2>0) $
              *ni[p]*                    localpos[1]^(ni[p]-1>0) $
             + sn $
              *mi[p]*                    localpos[0]^(mi[p]-1>0) $
              *ni[p]*(ni[p]-1)*          localpos[1]^(ni[p]-2>0))
        end
       end


       case op2 of
           0: begin
               sum = sum + field[p,elm]*temp
               if(threed eq 1) then begin
                   sum = sum + temp* $
                     (field[p+20,elm]*localpos[2]   $
                     +field[p+40,elm]*localpos[2]^2 $
                     +field[p+60,elm]*localpos[2]^3)
               endif
           end
           1: begin
               if(threed eq 1) then begin
                   sum = sum + temp* $
                     (field[p+20,elm]   $
                     +field[p+40,elm]*localpos[2]*2. $
                     +field[p+60,elm]*localpos[2]^2*3.)
               endif
           end
           2: begin
               if(threed eq 1) then begin
                   sum = sum + temp* $
                     (field[p+40,elm]*2. $
                     +field[p+60,elm]*localpos[2]*6.)
               endif
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
                     mask=mask, phi=phi0

   if(n_elements(phi0) eq 0) then phi0 = 0.
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

   if(p eq 1) then begin
       dx = 0
       dy = 0
   endif else begin
       dx = (xrange[1] - xrange[0]) / (p - 1.)
       dy = (yrange[1] - yrange[0]) / (p - 1.)
   endelse

   result = fltarr(p,p)
   mask = fltarr(p,p)
   mask[*] = 1.

   small = 1e-3
   localphi = 0.

   sz = size(mesh.elements._data, /dim)
   if(sz[0] gt 7) then begin
       threed = 1
   endif else begin
       threed = 0
   endelse

   ; for each triangle, evaluate points within triangle which fall on
   ; rectilinear output grid
   for i=long(0),nelms-1 do begin
       a = mesh.elements._data[0,i]
       b = mesh.elements._data[1,i]
       c = mesh.elements._data[2,i]
       t = mesh.elements._data[3,i]
       x = mesh.elements._data[4,i]
       y = mesh.elements._data[5,i]
       if(threed eq 1) then begin
           d = mesh.elements._data[7,i]
           phi = mesh.elements._data[8,i]
           localphi = phi0 - phi
           if(localphi lt 0 or localphi gt d) then  continue
       endif
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
                              -(pos[0]-x)*sn + (pos[1]-y)*co, $
                              localphi]

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

function read_lcfs, axis=axis, xpoint=xpoint, flux0=flux0, $
                    filename=filename, slice=time, last=last

   s = read_scalars(filename=filename)
   if(n_elements(time) eq 0) then begin
       slice = 0
   endif else slice = time
   if(keyword_set(last)) then begin
       slice = read_parameter("ntime", filename=filename)-1
   endif

   t0 = get_slice_time(filename=filename, slice=slice)

   tmp = s.time._data[*] - t0[0]
   dum = min(tmp, i, /abs)
   print, 'slice time = ', t0
   print, 'time step time: ', s.time._data[i]
   print, 'time slice: ', i
   
   xpoint = fltarr(2)
   axis = fltarr(2)
   xpoint[0] = s.xnull._data[i]
   xpoint[1] = s.znull._data[i]
   axis[0] = s.xmag._data[i]
   axis[1] = s.zmag._data[i]
   
   flux0 = s.psimin._data[i]

   return, s.psi_lcfs._data[i]
end

function read_field, name, x, y, t, slices=slices, mesh=mesh, $
                     filename=filename, points=pts, mask=mask, $
                     rrange=xrange, zrange=yrange,equilibrium=equilibrium, $
                     h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                     diff=diff, operation=op, complex=complex, fac=fac, $
                     linear=linear, last=last, average=average, linfac=linfac,$
                     dpsi=dpsi, symbol=symbol, units=units, cgs=cgs, mks=mks, $
                     real=real, imaginary=imag, edge_val=edge_val, phi=phi0, $
                     time=realtime, abs=abs, phase=phase, dimensions=d, $
                     flux_average=flux_av, rvector=rvector, zvector=zvector, $
                     taverage=taverage, is_nonlinear=is_nonlinear

   if(n_elements(slices) ne 0) then time=slices else time=0
   is_nonlinear = 0

   if(keyword_set(taverage)) then begin
       data = 0
       if(taverage eq 1) then taverage=16
       phi = 360.*findgen(taverage) / (taverage - 1.)
       for i=0, taverage-1 do begin
           data = data + $
             read_field(name, x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        diff=diff, operation=op, dimensions=d, $
                        linear=linear, last=last,symbol=symbol,units=units, $
                       cgs=cgs, mks=mks, time=realtime, $
                       rvector=rvector, zvector=zvector, phi=phi[i])
       end
       data = data/taverage
       return, data
   end

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
                        rrange=xrange, zrange=yrange, complex=complex, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        diff=diff, operation=op, dimensions=d, $
                        linear=linear, last=last,symbol=symbol,units=units, $
                       cgs=cgs, mks=mks, phi=phi0, time=realtime, $
                       rvector=rvector, zvector=zvector)
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
           if(n_elements(filename) eq 0) then filename='C1.h5'
           if(n_elements(filename) eq 1) then filename=replicate(filename,n)
       endif else n = 1

       data = 0
       for i=0, n-1 do begin
           data = data + $
             read_field(name, x, y, t, slices=time[i], mesh=mesh, $
                        filename=filename[i], points=pts, $
                        rrange=xrange, zrange=yrange, mask=mask, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        operation=op, complex=complex, dimensions=d, $
                        linear=linear, last=last,symbol=symbol,units=units, $
                       cgs=cgs, mks=mks, phi=phi0, time=realtime, $
                       rvector=rvector, zvector=zvector) $
             *((-1)^i)
       end

       symbol = '!7D!X' + symbol

       return, data
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(filename) gt 1) then filename=filename[0]
   if(n_elements(pts) eq 0) then pts = 200
   if(n_elements(op) eq 0) then op = 1

   if(hdf5_file_test(filename) eq 0) then return, 0

   version = read_parameter("version", filename=filename)
   nt = read_parameter("ntime", filename=filename)
   nv = read_parameter("numvar", filename=filename)
   itor = read_parameter("itor", filename=filename)
   ntor = read_parameter("ntor", filename=filename)
   version = read_parameter('version', filename=filename)
   ivform = read_parameter('ivform', filename=filename)
   icomplex = read_parameter('icomplex', filename=filename)
   i3d = read_parameter('3d', filename=filename)
   if(version eq 0) then begin
       xzero = read_parameter("xzero", filename=filename)
       zzero = read_parameter("zzero", filename=filename)
   endif else begin
       xzero = 0.
       zzero = 0.
   endelse
   ilin = read_parameter('linear', filename=filename)
   isubeq = read_parameter('eqsubtract', filename=filename)
   extsubtract = read_parameter('extsubtract', filename=filename)

   if(keyword_set(last)) then time = nt-1
   if(keyword_set(equilibrium)) then begin
       if(ilin eq 1) then time=-1
       if(isubeq eq 1) then linear = 0
   end

   if(time ge nt) then begin
       print, "Error: there are only ", nt-1, " time slices."
       return, 0
   endif

   realtime = get_slice_time(filename=filename, slice=time)

   data = fltarr(1, pts, pts)
   if(isubeq eq 1) then base = fltarr(pts,pts)

   d = dimensions()
   symbol=name
 
   print, 'Reading field ', name, ' at timeslice ', time
   print, 'Eqsubtract? ', isubeq
   print, string(form='(" linear=",I0,"; pts=",I0,";' + $
                 'equilibrium=",I0,"; complex=",I0,"; op=",I0)', $
                 keyword_set(linear), pts, keyword_set(equilibrium), $
                 keyword_set(complex), op)

   if(isubeq eq 1 and (not keyword_set(linear)) and (time ge 0)) $
     then begin
       data1 = read_field(name,x,y,t, slices=time, mesh=mesh, fac=fac, $
                          filename=filename, points=pts, mks=mks, cgs=cgs, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                          diff=diff, operation=op, linfac=linfac, $
                          /linear, last=last,symbol=symbol, $
                          units=units, dimensions=d, phi=phi0, $
                         rvector=rvector, zvector=zvector, is_nonlinear=isnl)
       if(isnl eq 1) then begin
          data = data1
       endif else begin
          t1 = t
          data0 = read_field(name,x,y,t, slices=-1, mesh=mesh, $
                             filename=filename, points=pts, fac=fac, $
                             rrange=xrange, zrange=yrange, complex=0, $
                             h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                             diff=diff, operation=op, mask=mask, $
                             symbol=symbol, mks=mks, cgs=cgs, $
                             units=units, dimensions=d, $
                             rvector=rvector, zvector=zvector)
          data = data0 + data1
          t = t1
       endelse
       return, data
   endif
   
   ; check if this is a primitive field
   file_id = h5f_open(filename)
   time_group_id = h5g_open(file_id, time_name(time))
   mesh = h5_parse(time_group_id, 'mesh', /read_data)
               
   field_group_id = h5g_open(time_group_id, 'fields')             
   nmembers = h5g_get_nmembers(time_group_id, 'fields')
   match = 0
   for m=0, nmembers-1 do begin
       thisname = h5g_get_member_name(time_group_id,'fields',m)
       if(strcmp(thisname, name, /fold_case) eq 1) then begin
           name = thisname
           match = 1
           break
       endif
   end
   
   if(match eq 1) then begin

       if(keyword_set(complex)) then begin
           h5g_close, field_group_id
           h5g_close, time_group_id
           h5f_close, file_id
           print, '  reading complex field.', ntor

           data_r = read_field(name,x,y,t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, complex=0, $
                               h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                               diff=diff, operation=op, mask=mask, $
                               /linear, last=last,symbol=symbol, $
                               units=units, dimensions=d, $
                               equilibrium=equilibrium)
           data_i = read_field(name+'_i',x,y,t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, complex=0, $
                               h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                               diff=diff, operation=op, $
                               /linear, last=last,symbol=symbol, $
                               units=units, dimensions=d, $
                               equilibrium=equilibrium)
           data = complex(data_r, data_i)


           ; evaluate at phi0
           ; it is okay to do this here since complex cases are always linear
           if(n_elements(phi0) ne 0) then begin
               print, 'evaluating at angle ', phi0, ' with ntor = ', ntor
               data = data* $
                 complex(cos(ntor*phi0*!pi/180.), sin(ntor*phi0*!pi/180.))
           end
       endif else begin
           print, '  reading real field'

           field = h5_parse(field_group_id, name, /read_data)
               
           time_id = h5a_open_name(time_group_id, "time")
           t = h5a_read(time_id)
           h5a_close, time_id
           h5g_close, field_group_id
           h5g_close, time_group_id
           h5f_close, file_id
               
           if(n_elements(phi0) eq 0) then phi_rad=0. $
           else phi0_rad = phi0*!pi/180.
  
           data[0,*,*] = $
             eval_field(field._data, mesh, points=pts, $
                        r=x, z=y, op=op, filename=filename, $
                        xrange=xrange, yrange=yrange, mask=mask, $
                        phi=phi0_rad)
           symbol = translate(name, units=d, itor=itor)

           if(version lt 5 and isubeq eq 1 and time ge 0 and $
              ((strcmp('te', name, /fold_case) eq 1) or $
               (strcmp('te_i', name, /fold_case) eq 1))) then begin
               zeff = read_parameter('zeff',filename=filename)
               data = data / zeff
               print, 'Correcting bug in linear Te for version < 4'
           end
       endelse

       
   endif else begin
       h5g_close, field_group_id
       h5g_close, time_group_id
       h5f_close, file_id

       print, '  reading composite field'

   if(strcmp('zero', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       data = 0.*psi
       symbol = ''

   ;==========================================
   ; local_beta = 2*P/B^2
   ;==========================================
   endif else if(strcmp('beta', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(n_elements(psi) le 1) then return, 0

       I = read_field('I',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)     
       P = read_field('P',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                      linear=linear, complex=complex, phi=phi0)

       if(extsubtract eq 1 and version lt 8) then begin
           I = I + read_field('I_ext', x, y, t, mesh=mesh, $
                              filename=filename, points=pts, slices=time, $
                              rrange=xrange, zrange=yrange, complex=complex, $
                              linear=linear, mask=mask, phi=phi0)
       end

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                     linear=linear,complex=complex,phi=phi0)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
   
       if(ivform eq 0) then begin
           data = v/r
       endif else if(ivform eq 1) then begin
           data = v*r
       endif
       symbol = '!8u!D!9P!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; toroidal velocity shear
   ;===========================================
   endif else if(strcmp('toroidal velocity shear', name, /fold_case) eq 1) or $
     (strcmp('vzp', name, /fold_case) eq 1) then begin
       
       v = read_field('V',x,y,t,slices=time, mesh=mesh, filename=filename, $
                        points=pts,rrange=xrange,zrange=yrange, linear=linear)
       psi = read_field('psi',x,y,t,slices=time, mesh=mesh, filename=filename,$
                        points=pts,rrange=xrange,zrange=yrange,/equilibrium)


       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
   
       if(ivform eq 0) then begin
           vz = v/r
       endif else if(ivform eq 1) then begin
           vz = v*r
       endif

       data = s_bracket(vz,psi,x,y)/sqrt(s_bracket(psi,psi,x,y))
       symbol = "!8u!D!9P!N'!X"
       d = dimensions(/t0, _EXTRA=extra)

   ;===========================================
   ; thermal velocity
   ;===========================================
   endif else if(strcmp('vt_i', name, /fold_case) eq 1) or $
     (strcmp('vti', name, /fold_case) eq 1) then begin

       Ti = read_field('Ti',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange)
         
       data = sqrt(2.*Ti)
       symbol = '!8v!Dti!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; sound speed
   ;===========================================
   endif else if(strcmp('sound speed', name, /fold_case) eq 1) or $
     (strcmp('cs', name, /fold_case) eq 1) then begin

       P = read_field('P',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       den = read_field('den',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       gam = read_parameter('gam', filename=filename)
  
       data = sqrt(gam*P/den)
       symbol = '!8c!Ds!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; Mach number
   ;===========================================
   endif else if(strcmp('mach', name, /fold_case) eq 1) or $
     (strcmp('m', name, /fold_case) eq 1) then begin

       cs = read_field('cs',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)

       phi = read_field('phi',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       V = read_field('V',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       chi = read_field('chi',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       if(ivform eq 0) then begin  
           v2 = s_bracket(phi,phi,x,y)/r^2 $
             + v^2/r^2 + s_bracket(chi,chi,x,y) $
             + 2.*a_bracket(chi,phi,x,y)/r
       endif else if(ivform eq 1) then begin
           v2 = r^2*s_bracket(phi,phi,x,y) $
             + r^2*v^2 + s_bracket(chi,chi,x,y)/r^4 $
             + 2.*a_bracket(chi,phi,x,y)/r
       endif
  
       data = sqrt(v2)/cs
       symbol = '!8M!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; electron temperature
   ;===========================================
  endif else if(strcmp('electron temperature', name, /fold_case) eq 1) or $
    (strcmp('te', name, /fold_case) eq 1) then begin

      Pe0 = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, linear=linear, $
                       /equilibrium)

      n0 = read_field('ne', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, linear=linear, $
                      /equilibrium)


      if(keyword_set(isubeq eq 1) and keyword_set(linear) and time ge 0) $
        then begin
          Pe1 = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, $
                           rrange=xrange, zrange=yrange, /linear, $
                           complex=complex, phi=phi0)
          
          n1 = read_field('ne', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange, /linear, $
                          complex=complex, phi=phi0)
          data = pe1/n0 - pe0*n1/n0^2
      endif else begin
          data = pe0/n0
      endelse

;       if(keyword_set(linear) and (isubeq eq 1) and (time ge 0)) then begin
;           Pe0 = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
;                            filename=filename, points=pts, $
;                            rrange=xrange, zrange=yrange, /equilibrium)

;           n0 = read_field('ne', x, y, t, slices=time, mesh=mesh, $
;                           filename=filename, points=pts,  $
;                           rrange=xrange, zrange=yrange, /equilibrium)

;           data = pe1/n0 - pe0*n1/n0^2
;       endif else data = pe1/n1

      symbol = '!8T!De!N!X'
      d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; electron density
   ;===========================================
   endif else if(strcmp('electron density', name, /fold_case) eq 1) or $
     (strcmp('ne', name, /fold_case) eq 1) then begin

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, linear=linear, $
                      complex=complex, phi=phi0)
       zeff = read_parameter("zeff", filename=filename)
       if(zeff eq 0) then zeff = 1.
       data = zeff*n
  
       symbol = '!8n!De!N!X'
       d = dimensions(/n0, _EXTRA=extra)

   ;===========================================
   ; displacement
   ;===========================================
   endif else if(strcmp('displacement', name, /fold_case) eq 1) then begin

       Te1 = read_field('Te', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, linear=linear, $
                       complex=complex, phi=phi0)
       

       Te0 = read_field('Te', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)

       psi0 = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
      
       tprime = s_bracket(Te0,psi0,x,y)

       data = -Te1/tprime*sqrt(s_bracket(psi0,psi0,x,y))

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then begin
           data(where(abs(tprime) lt 1e-6 or psi0 lt psis)) = 0.
       endif else begin
           data(where(abs(tprime) lt 1e-6 or psi0 gt psis)) = 0.
       endelse
  
       symbol = '!7n!N!X'
       d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; overlap
   ;===========================================
   endif else if(strcmp('overlap', name, /fold_case) eq 1) then begin

       xi = read_field('displacement', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, linear=linear, $
                       complex=complex, phi=phi0)

       psi0 = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then jac = -1 else jac = 1

       ; d(xi)/dr = d(xi)/dpsi * |grad(psi)|
       ;          = (<xi,psi>/<psi,psi>)*sqrt(<psi,psi>)
       ;          = <xi,psi>/sqrt(<psi,psi>)
       data = abs(s_bracket(xi,psi0,x,y)/sqrt(s_bracket(psi0,psi0,x,y)))
  
       symbol = '!3|!6d!7n!D!8r!N!6/dr!3|!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; linearity
   ;===========================================
   endif else if(strcmp('linearity', name, /fold_case) eq 1) then begin

       xi = read_field('displacement', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, linear=linear, $
                       complex=complex, phi=phi0)
       psi0 = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
       p0 = read_field('p', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)


       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then jac = -1 else jac = 1

       ; d(xi)/dr = d(xi)/dpsi * |grad(psi)|
       ;          = (<xi,psi>/<psi,psi>)*sqrt(<psi,psi>)
       ;          = <xi,psi>/sqrt(<psi,psi>)
       l = p0/s_bracket(p0,psi0,x,y)*sqrt(s_bracket(psi0,psi0,x,y))
       data = xi/l
  
       symbol = '!3|!7n!D!8r!N!3|!6/!8L!Dp!N!X'
       d = dimensions(_EXTRA=extra)


   ;===========================================
   ; psi_norm
   ;===========================================
   endif else if(strcmp('psi_norm', name, /fold_case) eq 1) then begin

       psi = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time)

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       
       data = (psi-flux0)/(psis-flux0)
 
       symbol = '!7W!D!8N!N!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; grad_psi_norm
   ;===========================================
   endif else if(strcmp('grad_psi_norm', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time, op=2)
       psi_z = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time, op=3)

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       
       data = sqrt(psi_r^2 + psi_z^2)/abs(psis-flux0)
 
       symbol = '!3|!9G!7W!3|!X'
       d = dimensions(l0=-1, _EXTRA=extra)

   ;===========================================
   ; ion temperature
   ;===========================================
;  endif else if(strcmp('ion temperature', name, /fold_case) eq 1) or $
;    (strcmp('ti', name, /fold_case) eq 1) then begin

;      P = read_field('P', x, y, t, slices=time, mesh=mesh, $
;                     filename=filename, points=pts, $
;                     rrange=xrange, zrange=yrange)

;      Pe = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
;                     filename=filename, points=pts, $
;                     rrange=xrange, zrange=yrange)

;      n = read_field('den', x, y, t, slices=time, mesh=mesh, $
;                     filename=filename, points=pts, $
;                     rrange=xrange, zrange=yrange)
;
;      data = (p-pe)/n
;      symbol = '!8T!Di!N!X'
;      d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; new ion temperature (as written by m3dc1 since 1/6/2011)
   ;===========================================
   endif else if(strcmp('ion temperature', name, /fold_case) eq 1) then begin

       ti = read_field('ti', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
  
       data = ti
       symbol = '!8T!Di!N!X'
       d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; ion pressure
   ;===========================================
   endif else if(strcmp('ion pressure', name, /fold_case) eq 1) or $
     (strcmp('pi', name, /fold_case) eq 1) then begin

       P = read_field('P', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       Pe = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
  
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
                      rrange=xrange, zrange=yrange)

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
  
       data = n*v
       symbol = '!8L!D!9P!N!X'
       d = dimensions(/n0, /v0, /l0, _EXTRA=extra)

   ;===========================================
   ; toroidal current
   ;===========================================
   endif else if(strcmp('jy', name, /fold_case) eq 1) then begin

       lp = read_field('psi', x, y, t, slices=time, mesh=mesh, op=7, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                      complex=complex,phi=phi0)
       psir = read_field('psi', x, y, t, slices=time, mesh=mesh, op=2, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                        complex=complex,phi=phi0)


       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
           data = -(lp - psir/r)/r
       endif else data = -lp

       symbol = '!8J!D!9P!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; toroidal current
   ;===========================================
   endif else if(strcmp('jy_plasma', name, /fold_case) eq 1) then begin

       lp = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=7, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                      complex=complex,phi=phi0)
       psir = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=2, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                        complex=complex,phi=phi0)


       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
           data = -(lp - psir/r)/r
       endif else data = -lp

       symbol = '!8J!D!9P!N!X'
       d = dimensions(/j0,_EXTRA=extra)


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

       xx = fltarr(n_elements(t),n_elements(x),n_elements(y))
       zz = fltarr(n_elements(t),n_elements(x),n_elements(y))
       for k=0, n_elements(t)-1 do begin
           for i=0, n_elements(y)-1 do xx[k,*,i] = x
           for i=0, n_elements(x)-1 do zz[k,i,*] = y
       end
       data = sqrt((xx-x0)^2 + (zz-z0)^2)

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

       xx = fltarr(n_elements(t),n_elements(x),n_elements(y))
       zz = fltarr(n_elements(t),n_elements(x),n_elements(y))
       for k=0, n_elements(t)-1 do begin
           for i=0, n_elements(y)-1 do xx[k,*,i] = x
           for i=0, n_elements(x)-1 do zz[k,i,*] = y
       end
       data = atan(zz-z0,xx-x0)

       symbol = '!7h!X'

   ;===========================================
   ; pest angle
   ;===========================================
   endif else if(strcmp('pest angle', name, /fold_case) eq 1) then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange,_EXTRA=extra)

       forward_function find_lcfs
       psilim = find_lcfs(psi, x, y,axis=axis, xpoint=xpoint, flux0=flux0, $
                          _EXTRA=extra, filename=filename)

       rr = radius_matrix(x,y,t)
       zz = z_matrix(x,y,t)

       forward_function flux_coord_field
       rrfc = flux_coord_field(rr, psi, x, y, t, slice=time, $
                                 flux=flux, angle=angle, /pest, $
                                 filename=filename, points=pts, $
                                 _EXTRA=extra)
       zzfc = flux_coord_field(zz, psi, x, y, t, slice=time, $
                                 flux=flux, angle=angle, /pest, $
                                 filename=filename, points=pts, $
                                 _EXTRA=extra)

       psinorm = (psi - flux0) / (psilim - flux0)
       data = psi*0.
       for i=0, n_elements(data)-1 do begin
           if(psinorm[i] gt 1.) then continue

           dist = (rrfc-rr[i])^2 + (zzfc-zz[i])^2
           dum = min(dist, j)
           n_guess = j/n_elements(angle)
           m_guess = j - n_guess*n_elements(angle)
           data[i] = angle(n_guess)
       end

       symbol = '!7h!D!6PEST!N!X'
       d = dimensions()

   ;===========================================
   ; Field strength
   ;===========================================
   endif else if( (strcmp('field strength', name, /fold_case) eq 1) $
                  or (strcmp('b', name, /fold_case) eq 1)) $
     then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, complex=complex)
       if(extsubtract eq 1 and version lt 8) then begin
           psi = psi + read_field('psi_ext', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, complex=complex)
       end

       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex)
       if(extsubtract eq 1 and version lt 8) then begin
           I = I + read_field('I_ext', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex)
       end

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(icomplex eq 1) then begin
           b2 = (s_bracket(psi,conj(psi),x,y) + I*conj(I))/r^2
           f = read_field('f', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex)
           if(extsubtract eq 1 and version lt 8) then begin
               f = f + $
                 read_field('f_ext', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex)
           end
           fp = complex(0., ntor)*f
           b2 = b2 + s_bracket(fp,conj(fp),x,y) $
             - a_bracket(fp, conj(psi),x,y)/r $
             - a_bracket(conj(fp), psi,x,y)/r
           b2 = real_part(b2)
           b2 = b2 / 2. ; this comes from the cos^2 dependence of the field
       endif else begin 
           b2 = s_bracket(psi,psi,x,y) + I^2/r^2
           if(i3d eq 1) then begin
               fp = read_field('f', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange, op=11)
               if(extsubtract eq 1 and version lt 8) then begin
                   fp = fp + $
                     read_field('f_ext', x, y, t, slices=time, mesh=mesh, $
                                filename=filename, points=pts, linear=linear, $
                                rrange=xrange, zrange=yrange, op=11)
               end

               b2 = b2 + s_bracket(fp,fp,x,y) $
                 - 2.*a_bracket(fp,psi,x,y)/r
           end           
       endelse

       data = sqrt(b2)
       symbol = '!3|!5B!3|!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; Field energy
   ;===========================================
   endif else if( (strcmp('field energy', name, /fold_case) eq 1) $
                  or (strcmp('b2', name, /fold_case) eq 1)) $
     then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, complex=complex)

       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       b2 = (s_bracket(psi,conj(psi),x,y) + I*conj(I))/r^2
       if(icomplex eq 1) then begin
           ; if the fields are ~exp(i n phi), then
           ; this is the toroidally-averaged value of |B| !

           f = read_field('f', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex)
           fp = complex(0., ntor)*f
           b2 = b2 + s_bracket(fp,conj(fp),x,y) $
             - a_bracket(fp, conj(psi),x,y)/r $
             - a_bracket(conj(fp), psi,x,y)/r

           b2 = b2 / 2. ; this comes from the cos^2 dependence of the field
           b2 = real_part(b2)
       endif

       data = b2/(8.*!pi)
       symbol = '!3|!5B!3|!U!62!N!X'
       d = dimensions(p0=1, _EXTRA=extra)

   ;===========================================
   ; Poloidal Field strength
   ;===========================================
   endif else if( (strcmp('poloidal field strength', name, /fold_case) eq 1) $
                  or (strcmp('bp', name, /fold_case) eq 1)) $
     then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                        rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       v = read_field('v', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                       rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       data = b/sqrt(den)
       symbol = '!8v!DA!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; (minor) radial current density
   ;===========================================
   endif else if(strcmp('jn', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange, op=2)
       psi0_z = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange, op=3)
       i_r = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, op=2)
       i_z = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, op=3)


       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       psipsi = sqrt(psi0_r^2 + psi0_z^2)
       
       data = (i_z*psi0_r - i_r*psi0_z)/(r*psipsi)

       if(ntor ne 0) then begin

           psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=2)
           psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=3)
           f_r = read_field('f', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=2)
           f_z = read_field('f', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=3)

           data = data $
             - complex(0,ntor)*(psi_r*psi0_r + psi_z*psi0_z)/(r^2*psipsi) $
             - ntor^2*(f_z*psi0_r - f_r*psi0_z)/(r*psipsi)
       endif

       symbol = '!8J!Dr!N!X'
       d = dimensions(/j0, _EXTRA=extra)

   ;===========================================
   ; (major) radial current density
   ;===========================================
   endif else if(strcmp('jr', name, /fold_case) eq 1) then begin
       
       i_z = read_field('i', x, y, t, slices=time, mesh=mesh, op=3, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = -i_z / r

       if(ntor ne 0) then begin
           psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, op=2, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex)

           f_z = read_field('f', x, y, t, slices=time, mesh=mesh, op=3, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex)

           data = data + ntor^2 * f_z / r + complex(0., ntor)*psi_r/r^2
       endif
       
       symbol = '!8J!DR!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; vertical current density
   ;===========================================
   endif else if(strcmp('jz', name, /fold_case) eq 1) then begin
       
       i_r = read_field('i', x, y, t, slices=time, mesh=mesh, op=2, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = i_r / r

       if(ntor ne 0) then begin
           psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, op=3, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex)

           f_r = read_field('f', x, y, t, slices=time, mesh=mesh, op=2, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex)

           data = data - ntor^2 * f_r / r + complex(0., ntor)*psi_z/r^2
       endif
       
       symbol = '!8J!DZ!N!X'
       d = dimensions(/j0,_EXTRA=extra)


   ;===========================================
   ; poloidal current density
   ;===========================================
   endif else if(strcmp('jp', name, /fold_case) eq 1) then begin
       
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear,  $
                      rrange=xrange, zrange=yrange, phi=phi0)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       data = s_bracket(i,psi,x,y)/(r*sqrt(s_bracket(psi,psi,x,y)))
       symbol = '!8J!Dp!N!X'
       d = dimensions(/j0,_EXTRA=extra)

       
   ;===========================================
   ; del*(psi)
   ;===========================================
   endif else if(strcmp('jphi', name, /fold_case) eq 1) then begin

       psi_lp = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, op=7)
       data = psi_lp
       if(itor eq 1) then begin
           psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                              filename=filename, points=pts, $
                              rrange=xrange, zrange=yrange, op=2)
           r = radius_matrix(x,y,t)
           data = data - psi_r / r
       end

       symbol = translate('jphi', units=d, itor=itor)
       d = dimensions(/j0,l0=itor,_EXTRA=extra)

   ;===========================================
   ; vorticity
   ;===========================================
      endif else if(strcmp('vor', name, /fold_case) eq 1) then begin

          phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, $
                           complex=complex,phi=phi0)

          data = grad_shafranov(phi,x,y,tor=itor)
          symbol = translate('vor', units=d, itor=itor)

   ;===========================================
   ; divergence
   ;===========================================
     endif else if(strcmp('com', name, /fold_case) eq 1) then begin

         chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          phi=phi0)

         data = laplacian(chi,x,y,tor=itor)
         symbol = translate('com', units=d, itor=itor)

   ;===========================================
   ; R^2 vorticity
   ;===========================================
      endif else if(strcmp('r2vor', name, /fold_case) eq 1) then begin

          phi_lp = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, $
                           complex=complex,phi=phi0,op=7)
          phi_r = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, $
                           complex=complex,phi=phi0,op=2)
          r = radius_matrix(x,y,t)

          data = r^2*phi_lp + 2.*r*phi_r
          symbol = translate('vor', units=d, itor=itor)

   ;===========================================
   ; helicity
   ;===========================================
      endif else if(strcmp('helicity', name, /fold_case) eq 1) then begin

          psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, phi=phi0, $
                            complex=complex, linear=linear)
          i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, phi=phi0, $
                            complex=complex, linear=linear)
          f = read_field('f', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, phi=phi0, $
                            complex=complex, linear=linear)

          r = radius_matrix(x,y,t)

          data = (psi*conj(i) + conj(psi)*i)/r^4 $
            - s_bracket(f, conj(psi), x, y)/r^2 $
            - s_bracket(conj(f), psi, x, y)/r^2
          d = dimensions(l0=4, b0=2)
          symbol = translate('helicity', units=d, itor=itor)


   ;===========================================
   ; chi_perp
   ;===========================================
     endif else if(strcmp('chi_perp', name, /fold_case) eq 1) then begin

         kappa = read_field('kappa', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, phi=phi0)
         den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, /equilibrium)

         data = kappa/den
         d = dimensions(l0=2, t0=-1, _EXTRA=extra)
         symbol = '!7v!X'

   ;===========================================
   ; mu_perp
   ;===========================================
     endif else if(strcmp('mu_perp', name, /fold_case) eq 1) then begin

         visc = read_field('visc', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, phi=phi0)
         den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, /equilibrium)

         data = visc/den
         d = dimensions(l0=2, t0=-1, _EXTRA=extra)
         symbol = '!7l!D!9x!N!X'

   ;===========================================
   ; rotational transform
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

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       bt = sqrt(I^2/r^2)
       bp = sqrt(s_bracket(psi,psi,x,y)/r^2)

       data = 2.*!pi*(r * bp) / minor_r * bt
       symbol = '!8i!X'

   ;===========================================
   ; angular velocity
   ;===========================================
   endif else if(strcmp('omega', name, /fold_case) eq 1) then begin

       v = read_field('v', x, y, t, slices=time, mesh=mesh, complex=complex, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, op=op, phi=phi0)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           if(n_elements(op) ne 0) then begin
               print, 'Warning: using op on omega'
           end
           data = v/r^2
       endif else if(ivform eq 1) then begin
           data = v
       endif
       symbol = '!7X!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; pprime
   ;===========================================
   endif else if(strcmp('pprime', name, /fold_case) eq 1) then begin

       p = read_field('p', x, y, t, slices=time, mesh=mesh, complex=complex, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, op=op, phi=phi0)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, /equilibrium, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, op=op, phi=phi0)

       data = s_bracket(p,psi,x,y)/s_bracket(psi,psi,x,y)

       symbol = "!8p'!X"
       d = dimensions(p0=1, b0=1, l0=1+itor,_EXTRA=extra)


   ;===========================================
   ; electron angular velocity
   ;===========================================
   endif else if(strcmp('omega_e', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename)
       print, 'db = ', db

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       jy = read_field('jy', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = omega - (db*jy/den)/r
       endif else if(ivform eq 1) then begin
           data = omega - (db*jy/den)
       endif
       symbol = '!7x!D!8e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;==========================================================
   ; v_omega (the omega in v = r^2 omega grad(phi) + (K/n) B
   ;==========================================================
   endif else if(strcmp('v_omega', name, /fold_case) eq 1) then begin

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       data = omega - i/(r^2*s_bracket(psi,psi,x,y)) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r)

       symbol = '!7x!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;==========================================================
   ; v_K (the K in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('v_K', name, /fold_case) eq 1) then begin

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = den/s_bracket(psi,psi,x,y) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r)

       symbol = '!8K!X'
       d = dimensions(/v0, /n0, b0=-1, _EXTRA=extra)

   ;==========================================================
   ; v_K (the K in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('v_K_n', name, /fold_case) eq 1) then begin

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = 1./s_bracket(psi,psi,x,y) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r)

       symbol = '!8K!6/!8n!X'
       d = dimensions(/v0, b0=-1, _EXTRA=extra)

   ;==========================================================
   ; ve_omega (the omega in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('ve_omega', name, /fold_case) eq 1) then begin

       di = read_parameter('db', filename=filename)
       print, 'di = ', di
       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi_lp = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, op=7)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       data = omega - i/(r^2*s_bracket(psi,psi,x,y)) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r $
          - (di/den)*s_bracket(i,psi,x,y)) $
         + (di/den)*(psi_lp - itor*dx(psi,x)/r)/r^2

       symbol = '!7x!D!8e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;==========================================================
   ; ve_K (the K in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('ve_K', name, /fold_case) eq 1) then begin

       di = read_parameter('db', filename=filename)
       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = 1./s_bracket(psi,psi,x,y) * $
         (den*(r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r) $
          - di*s_bracket(i,psi,x,y))

       symbol = '!8K!D!8e!N!X'
       d = dimensions(/v0, /n0, b0=-1, _EXTRA=extra)


   ;===========================================
   ; electron angular velocity
   ;===========================================
   endif else if(strcmp('omega_perp_e', name, /fold_case) eq 1) then begin

       omega = read_field('v_omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       omega_star = read_field('omega_*', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)

       data = omega - omega_star
;       symbol = '!7x!S!D!9x!N!S!U!8e!N!X'
       symbol = '!7x!D!8e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; cyclotron frequency
   ;===========================================
   endif else if(strcmp('omega_ci', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', db
       if(db eq 0.) then begin
           print, 'Warning: Assuming d_i = 1.'
           db = 1.
       endif

       B = read_field('B', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = B/db
       symbol = '!7x!8!Dc!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; diamagnetic frequency
   ;===========================================
   endif else if(strcmp('omega_*', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', filename, db

       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)


       data = db*s_bracket(p,psi,x,y)/s_bracket(psi,psi,x,y) / den

       symbol = '!7x!6!D*!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; diamagnetic frequency
   ;===========================================
   endif else if(strcmp('omega_*i', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', filename, db

       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       pe = read_field('pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       data = db*s_bracket(p-pe,psi,x,y)/s_bracket(psi,psi,x,y) / den

       symbol = '!7x!6!D*i!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; diamagnetic frequency
   ;===========================================
   endif else if(strcmp('omega_*e', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', filename, db

       pe = read_field('pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       data = -db*s_bracket(pe,psi,x,y)/s_bracket(psi,psi,x,y) / den

       symbol = '!7x!6!D*e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; ExB frequency
   ;===========================================
   endif else if(strcmp('omega_ExB', name, /fold_case) eq 1) then begin

       omega = read_field('v_omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       w_star_i =  read_field('omega_*i', x, y, t, slices=time, mesh=mesh, $
                              filename=filename, points=pts, $
                              rrange=xrange, zrange=yrange)
       
       data = omega - w_star_i

       symbol = '!7x!6!DE!9X!6B!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)


   ;===========================================
   ; Larmor radius
   ;===========================================
   endif else if(strcmp('rho_i', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)

       Ti = read_field('Ti', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       B = read_field('B', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)


       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = db*sqrt(2.*Ti)/B
       symbol = '!7q!8!Di!N!X'
       d = dimensions(l0=1, _EXTRA=extra)


   ;===========================================
   ; parallel thermal gradient
   ;===========================================
   endif else if(strcmp('xbdotgradt', name, /fold_case) eq 1) then begin

       Te0 = read_field('Te', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange)
       psi0 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, /equilibrium, $
                         rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ilin eq 1) then begin
           Te1 = read_field('Te', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange)
           psi1 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, linear=linear, $
                           rrange=xrange, zrange=yrange)
       
           data = a_bracket(Te1, psi0, x, y)/r + a_bracket(Te0, psi1, x, y)/r

           if(ntor ne 0) then begin
               i1 = read_field('i', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange)
               i0 = read_field('i', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, /equilibrium, $
                               rrange=xrange, zrange=yrange)
               f1 = read_field('f', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange)

               data = data + complex(0.,ntor)* $
                 (Te0*i1 + Te1*i0 - s_bracket(Te0, f1, x, y))
           end
       end else begin
           data = a_bracket(Te0, psi0, x, y)/r + a_bracket(Te0, psi0, x, y)/r           
       end

;       data = te0

       symbol = '!8B!9.G!8T!X'
       d = dimensions(l0=1, _EXTRA=extra)

   ;===========================================
   ; parallel pressure
   ;===========================================
   endif else if(strcmp('xbdotgradp', name, /fold_case) eq 1) then begin

       if(ilin eq 1) then begin
           p1 = read_field('p', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex)
           p0 = read_field('p', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, /equilibrium, $
                            rrange=xrange, zrange=yrange, complex=complex)
           psi1 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, linear=linear, $
                           rrange=xrange, zrange=yrange, complex=complex)
           psi0 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, /equilibrium, $
                           rrange=xrange, zrange=yrange, complex=complex)
           if(itor eq 1) then begin
               r = radius_matrix(x,y,t)
           endif else r = 1.
       
           data = a_bracket(p1, psi0, x, y)/r + a_bracket(p0, psi1, x, y)/r

           if(ntor ne 0) then begin
               i0 = read_field('i', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, /equilibrium, $
                               rrange=xrange, zrange=yrange, complex=complex)
               f1 = read_field('f', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange, complex=complex)

               data = data + complex(0.,ntor)* $
                 (p1*i0/r^2 - s_bracket(p0, f1, x, y))
           end
       end

;       data = te0

       symbol = '!8B!9.G!8p!X'
       d = dimensions(p0=1, b0=1, l0=-1, _EXTRA=extra)


   ;===========================================
   ; parallel flow
   ;===========================================
   endif else if(strcmp('vpar', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       w = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)
         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       
       b2 = (s_bracket(psi,psi,x,y) + I^2)/r^2
       if(ivform eq 0) then begin
           data = ((s_bracket(phi,psi,x,y) + w*I)/r^2 $
                   + a_bracket(chi,psi,x,y)/r)/sqrt(b2)
       endif else begin
           data = (s_bracket(phi,psi,x,y) + w*I $
                   + a_bracket(chi,psi,x,y)/r^3)/sqrt(b2)
       endelse
       symbol = '!8u!D!9#!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; poloidal flow
   ;===========================================
   endif else if(strcmp('vpol', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       
       if(ivform eq 0) then begin
           data = 0.
       endif else begin
           data = r*(s_bracket(phi,psi,x,y) $
                   + a_bracket(chi,psi,x,y)/r^3) $
             / sqrt(s_bracket(psi,psi,x,y))
       endelse
       symbol = '!8u!D!7h!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('omtest', name, /fold_case) eq 1) then begin

       v_omega = read_field('v_omega', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       v_k_n = read_field('v_k_n', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)

         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = v_omega + i*v_k_n/r^2
       
       symbol = '!7X!D!6test!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)


   ;===========================================
   ; perpendicular flow v.(B x grad(psi))
   ;===========================================
   endif else if(strcmp('vperp', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       w = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)
         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       
       pp = s_bracket(psi,psi,x,y)
       b2 = (pp + I^2)/r^2
       if(ivform eq 0) then begin
           data = 0.
       endif else begin
           data = (-I*s_bracket(phi,psi,x,y) + w*s_bracket(psi,psi,x,y) $
                   -I*a_bracket(chi,psi,x,y)/r^3)/sqrt(b2*pp)
       endelse
       symbol = '!8u!D!9x!N!X'
       d = dimensions(/v0, _EXTRA=extra)

         
   ;===========================================
   ; radial flow
   ;===========================================
   endif else if(strcmp('vr', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange, complex=complex)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange,complex=complex)
       psi = read_field('psi', x, y, t, /equilibrium, slices=time, mesh=mesh, $
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                        rrange=xrange, zrange=yrange)
       phi_z = read_field('phi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=3, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = -phi_z/r + chi_r
       endif else if (ivform eq 1) then begin
           data = -r*phi_z + chi_r/r^2
       endif
       symbol = '!8u!DR!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('vy', name, /fold_case) eq 1) then begin

       chi_z = read_field('chi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=3, $
                        rrange=xrange, zrange=yrange)
       phi_r = read_field('phi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=2, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = phi_r/r + chi_z
       endif else if (ivform eq 1) then begin
           data = r*phi_r + chi_z/r^2
       endif
       symbol = '!8u!DZ!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   ;===========================================
   ; radial field
   ;===========================================
   endif else if(strcmp('bn', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=2)
       psi0_z = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=3)
       psi_r = read_field('psi', x, y, t, mesh=mesh, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=2)
       psi_z = read_field('psi', x, y, t, mesh=mesh, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=3)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = (psi_z*psi0_r - psi_r*psi0_z)/r

        if(ntor ne 0) then begin
            f_r = read_field('f', x, y, t, mesh=mesh, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=2)
            f_z = read_field('f', x, y, t, mesh=mesh, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=3)
            data = data + complex(0.,ntor)*(f_r*psi0_r + f_z*psi0_z)
        endif

       data = data / sqrt(psi0_r^2 + psi0_z^2)
       symbol = '!8B!Dn!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; poloidal field
   ;===========================================
   endif else if(strcmp('bpol', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=2)
       psi0_z = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=3)
       psi_r = read_field('psi', x, y, t, mesh=mesh, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=2)
       psi_z = read_field('psi', x, y, t, mesh=mesh, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=3)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = (psi_z*psi0_z + psi_r*psi0_r)/r

        if(ntor ne 0) then begin
            f_r = read_field('f', x, y, t, mesh=mesh, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=2)
            f_z = read_field('f', x, y, t, mesh=mesh, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=3)
            data = data + complex(0.,ntor)*(f_z*psi0_r - f_r*psi0_z)
        endif

       data = data / sqrt(psi0_r^2 + psi0_z^2)
       symbol = '!8B!D!7h!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; (Major) Radial field
   ;===========================================
   endif else if(strcmp('bx', name, /fold_case) eq 1) then begin

       psi_z = read_field('psi', x, y, t, mesh=mesh, operation=3, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        linear=linear, mask=mask, phi=phi0)

       if(extsubtract eq 1 and version lt 8) then begin
           psi_z = psi_z + read_field('psi_ext', x, y, t, mesh=mesh, operation=3, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, mask=mask, phi=phi0)
       end

       if(i3d eq 1) then begin
           f_rp = read_field('f', x, y, t, mesh=mesh, operation=12, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0)
           if(extsubtract eq 1 and version lt 8) then begin
               f_rp = f_rp + read_field('f_ext', x, y, t, mesh=mesh, operation=12, $
                                        filename=filename, points=pts, slices=time, $
                                        rrange=xrange, zrange=yrange, complex=complex, $
                                        linear=linear, phi=phi0)
           end

       endif else if(icomplex eq 1) then begin
           f_r = read_field('f', x, y, t, mesh=mesh, operation=2, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0)
           if(extsubtract eq 1 and version lt 8) then begin
               f_r = f_r + read_field('f_ext', x, y, t, mesh=mesh, operation=2, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, phi=phi0)
           end

           f_rp = complex(0.,ntor)*f_r
       endif else f_rp = 0.
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = -psi_z/r - f_rp
       symbol = '!8B!DR!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; Vertical field
   ;===========================================
   endif else if(strcmp('by', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi', x, y, t, mesh=mesh, operation=2, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        linear=linear, mask=mask, phi=phi0)

       if(extsubtract eq 1 and version lt 8) then begin
           psi_r = psi_r + read_field('psi_ext', x, y, t, mesh=mesh, operation=2, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, mask=mask, phi=phi0)
       end

       if(i3d eq 1) then begin
           f_zp = read_field('f', x, y, t, mesh=mesh, operation=13, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0)
           if(extsubtract eq 1 and version lt 8) then begin
               f_zp = f_zp + read_field('f_ext', x, y, t, mesh=mesh, operation=13, $
                                 filename=filename, points=pts, slices=time, $
                                 rrange=xrange, zrange=yrange, complex=complex, $
                                 linear=linear, phi=phi0)
           end

       endif else if(icomplex eq 1) then begin
           f_z = read_field('f', x, y, t, mesh=mesh, operation=3, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0)
           if(extsubtract eq 1 and version lt 8) then begin
               f_z = f_z + read_field('f_ext', x, y, t, mesh=mesh, operation=3, $
                                 filename=filename, points=pts, slices=time, $
                                 rrange=xrange, zrange=yrange, complex=complex, $
                                 linear=linear, phi=phi0)
           end
           f_zp = complex(0.,ntor)*f_z
       endif else f_zp = 0.
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = psi_r/r - f_zp
       symbol = '!8B!DZ!N!X'
       d = dimensions(/b0, _EXTRA=extra)


   ;===========================================
   ; poloidal flow
   ;===========================================
   endif else if(strcmp('vp', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       V   = read_field('V',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       psipsi = s_bracket(psi,psi,x,y)
       if(ivform eq 0) then begin
           b2 = psipsi + i^2
           data = -(i*s_bracket(phi,psi,x,y) - v*psipsi $
                   +i*a_bracket(chi,psi,x,y)*r) / (r^2 * sqrt(psipsi*b2))
       endif else begin
           b2 = (psipsi + i^2)/r^2
           data = -(i*s_bracket(phi,psi,x,y) - v*psipsi $
                   + i*a_bracket(chi,psi,x,y)/r^3) / sqrt(b2*psipsi)

       endelse
       symbol = '!8u!Ds!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('v_star', name, /fold_case) eq 1) then begin
       rho_i = read_field('rho_i', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       omega_ci = read_field('omega_ci', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       kr = sqrt(s_bracket(p,p,x,y))/p

       data = rho_i^2 * omega_ci * kr
       symbol = '!8u!6!D*!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('omega_star', name, /fold_case) eq 1) then begin
       v_star = read_field('v_star', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1

       data = v_star/r
       symbol = '!7x!6!D*!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   endif else if(strcmp('bp_over_b', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange)

       bp = s_bracket(psi,psi,x,y)
       data = sqrt(bp)/(sqrt(bp + i^2))
       symbol = '!3|!5B!D!8p!N!3|/|!8B!3|!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; ideal_k
   ;===========================================
   endif else if(strcmp('k', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       n   = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
         
       psipsi = s_bracket(psi,psi,x,y)

       data = n*(s_bracket(phi,psi,x,y) + r*a_bracket(chi,psi,x,y))/psipsi

   ;===========================================
   ; ideal omega
   ;===========================================
   endif else if(strcmp('ideal omega', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       v   = read_field('v',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                        rrange=xrange, zrange=yrange)
       eta= read_field('eta',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
                        rrange=xrange, zrange=yrange)
       v   = read_field('v'  , x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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

   endif else if(strcmp('jpar', name, /fold_case) eq 1) then begin

       psi0 = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       jphi = read_field('jphi', x, y, t, linear=linear, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex)
       i0 = read_field('i', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, /equilibrium, $
                      complex=complex)
       
       r = radius_matrix(x,y,t)
       b0 = s_bracket(psi0,psi0,x,y)/r^2 + i0^2/r^2
       data = (s_bracket(i,psi0,x,y)/r^2 - jphi*i0/r^2)/sqrt(b0)

       symbol = '!8J!D!3||!6!N!X'
       d = dimensions(j0=1,_EXTRA=extra)

   endif else if(strcmp('jpar_B', name, /fold_case) eq 1) then begin

       psi0 = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       jphi = read_field('jphi', x, y, t, linear=linear, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex)
       i0 = read_field('i', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, /equilibrium, $
                      complex=complex)
       
       r = radius_matrix(x,y,t)
       b0 = s_bracket(psi0,psi0,x,y)/r^2 + i0^2/r^2
       data = (s_bracket(i,psi0,x,y)/r^2 - jphi*i0/r^2)/b0

       symbol = '!8J!D!9#!N!3/!8B!X'
       d = dimensions(j0=1,b0=-1,_EXTRA=extra)

   endif else if(strcmp('jb', name, /fold_case) eq 1) then begin

       psi0 = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       jy = read_field('jy', x, y, t, linear=linear, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex)
       i0 = read_field('i', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, /equilibrium, $
                      complex=complex)
       
       r = radius_matrix(x,y,t)
       data = s_bracket(i,psi0,x,y)/r^2 + jy*i0/r

       symbol = '!8J!D!9#!N!8B!X'
       d = dimensions(j0=1,b0=1,_EXTRA=extra)

   ;===========================================
   ; particle flux
   ;===========================================
   endif else if(strcmp('flux_nv', name, /fold_case) eq 1) then begin

       den = read_field('den', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       f_p = read_field('f', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=11)
       
       r = radius_matrix(x,y,t)

       bp2 = s_bracket(psi,psi,x,y)/r^2 + 2.*a_bracket(psi,f_p,x,y)/r + s_bracket(f_p,f_p,x,y)
       bpdotv = (r*s_bracket(phi,psi,x,y) + a_bracket(chi,psi,x,y)/r^3 $
                 + a_bracket(f_p,psi,x,y)/r - s_bracket(chi,f_p,x,y)/r^2)/sqrt(bp2)
       data = den*bpdotv
       d = dimensions(/n0, /v0)
       symbol = '!6Parallel Particle Flux!X'

   ;===========================================
   ; particle flux
   ;===========================================
   endif else if(strcmp('flux_nv2', name, /fold_case) eq 1) then begin

       den=read_field('den', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       chi=read_field('chi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       psi0=read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, /equilibrium)

       
       r = radius_matrix(x,y,t)

       data = -conj(den)*(r*a_bracket(psi0,u,x,y) + s_bracket(psi0,chi,x,y)/r)
       data = data / sqrt(s_bracket(psi0,psi0,x,y))
       d = dimensions(l0=-3,t0=-1)
       symbol = '!6Particle Flux!X'

   ;===========================================
   ; parallel power flux
   ;===========================================
    endif else if((strcmp('parallel heat flux', name, /fold_case) eq 1) or $
       (strcmp('qpar', name, /fold_case) eq 1)) then begin

       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange, phi=phi0)
       den = read_field('den', x, y, t, slices=time, mesh=mesh,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       p_p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0, op=11)
       den_p = read_field('den', x, y, t, slices=time,mesh=mesh,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0, op=11)
;;        te = read_field('te', x, y, t, slices=time, mesh=mesh, $
;;                         filename=filename, points=pts, complex=complex, $
;;                         rrange=xrange, zrange=yrange, phi=phi0)
;;        te_p = read_field('te', x, y, t, slices=time,mesh=mesh,$
;;                         filename=filename, points=pts, complex=complex, $
;;                         rrange=xrange, zrange=yrange, phi=phi0, op=11)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       f_p = read_field('f', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=11)
       kappar = read_parameter(filename=filename, 'kappar')
       zeff = read_parameter(filename=filename, 'zeff')

       r = radius_matrix(x,y,t)

;       den = 1.

;;        if(keyword_set(linear)) then begin
;;           print, "LINEAR IS SET"

;;           p0 = read_field('p', x, y, t, slice=-1, mesh=mesh, $
;;                           filename=filename, points=pts, $
;;                           rrange=xrange, zrange=yrange, phi=phi0)
;;           den0 = read_field('den', x, y, t, slice=-1, mesh=mesh,  $
;;                             filename=filename, points=pts, $
;;                             rrange=xrange, zrange=yrange, phi=phi0)
;;           psi0 = read_field('psi', x, y, t, slices=-1, mesh=mesh, phi=phi0, $
;;                             filename=filename, points=pts, $
;;                             rrange=xrange, zrange=yrange)
;;           i0 = read_field('i', x, y, t, slices=-1, mesh=mesh, phi=phi0, $
;;                           filename=filename, points=pts, $
;;                           rrange=xrange, zrange=yrange)
          
;;           p = p + p0
;;           den = den + den0
;;           psi = psi + psi0
;;           i = i + i0
          
;;           te0 = p0 / den0
;;           b20 = s_bracket(psi0,psi0,x,y)/r^2 + i0^2/r^2
;;           bdotgradte0 = a_bracket(te0, psi0, x, y)/r
;;        end

;       den = den*zeff
;       den_p = den_p*zeff
       
       te = p / den
       te_p = p_p / den ;- p*den_p / den^2
       bp2 = s_bracket(psi,psi,x,y)/r^2 + 2.*a_bracket(psi,f_p,x,y)/r + s_bracket(f_p,f_p,x,y)
       b2 = bp2 + i^2/r^2
       bdotgradte = a_bracket(te, psi, x, y)/r $
                    + i*te_p/r^2 - s_bracket(f_p, te, x, y)

       br = -dz(psi,y)/r - dx(f_p,x)
       bbter = br*bdotgradte/b2
       bz =  dx(psi,x)/r - dz(f_p,y)
       bbtez = bz*bdotgradte/b2
;;        if(keyword_set(linear)) then begin
;;           br0 = -dz(psi0,y)/r
;;           bz0 =  dx(psi0,x)/r
;;           bbter = 0.*bbter - br0*bdotgradte0/b20
;;           bbtez = 0.*bbtez - bz0*bdotgradte0/b20
;;        endif
       
       if(keyword_set(rvector)) then begin
          data = -kappar*bbter
          symbol = '!6q!D!9#!N.G!8R!X'
       endif else if(keyword_set(zvector)) then begin
          data = -kappar*bbtez
          symbol = '!6q!D!9#!N.G!8Z!X'
       endif else begin
;          data = kappar*sqrt(abs(bbtez)^2 + abs(bbter)^2)
;          symbol = '!3|!6q!D!9#!N!3|!X'
;          data = -kappar*bdotgradte/sqrt(b2)
;          data = -kappar*bdotgradte*bp2/b2
          data = -kappar*bdotgradte*sqrt(bp2)/b2
          symbol = '!6q!D!9#!N!X'
       endelse
       
       d = dimensions(/p0,/v0)
       is_nonlinear = 1

   ;===========================================
   ; convective power flux
   ;===========================================
    endif else if((strcmp('convective heat flux', name, /fold_case) eq 1) or $
       (strcmp('qcon', name, /fold_case) eq 1)) then begin

       p=read_field('p', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       phi = read_field('phi', x, y, t, slices=time, mesh=mesh,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       gamma = read_parameter(filename=filename, 'gam')

       r = radius_matrix(x,y,t)

       if(keyword_set(rvector)) then begin
          v = -r*dz(phi,y) + dx(chi,x)/r^2
          symbol = '!6q!D!8V!N!9.G!8R!X'
       endif else if(keyword_set(zvector)) then begin
          v =  r*dx(phi,x) + dz(chi,y)/r^2
          symbol = '!6q!D!8V!N!9.G!8Z!X'
       endif else begin
          omega = read_field('omega', x, y, t, slices=time, mesh=mesh,  $
                             filename=filename, points=pts, complex=complex, $
                             rrange=xrange, zrange=yrange)
          v =  sqrt(r^2*s_bracket(phi,phi,x,y) + r^2*omega^2 + s_bracket(chi,chi,x,y)/r^4 $
            + 2.*a_bracket(chi,phi,x,y)/r)
          symbol = '!3|!6q!D!9#!N!3|!X'
       endelse
       
       data = -gamma/(gamma-1.) * p * v
       d = dimensions(/p0,/v0)


   ;===========================================
   ; dbndt
   ;===========================================
   endif else if(strcmp('dbndt', name, /fold_case) eq 1) then begin

;        psi0 = read_field('psi',x,y,t,slices=time,mesh=mesh, linear=linear,  $
;                          filename=filename, points=pts, /equilibrium, $
;                          rrange=xrange, zrange=yrange)
       psi0_r = read_field('psi',x,y,t,slices=time,mesh=mesh,linear=linear,  $
                           filename=filename, points=pts, /equilibrium, $
                           rrange=xrange, zrange=yrange,op=2)
       psi0_z = read_field('psi',x,y,t,slices=time,mesh=mesh,linear=linear,  $
                           filename=filename, points=pts, /equilibrium, $
                           rrange=xrange, zrange=yrange,op=3)
       i0 = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                       filename=filename, points=pts, /equilibrium, $
                       rrange=xrange, zrange=yrange)
       i0_r = read_field('i',x,y,t, slices=time, mesh=mesh, linear=linear,  $
                         filename=filename, points=pts, /equilibrium, $
                         rrange=xrange, zrange=yrange,op=2)
       i0_z = read_field('i', x,y,t, slices=time, mesh=mesh, linear=linear,  $
                         filename=filename, points=pts, /equilibrium, $
                         rrange=xrange, zrange=yrange,op=3)
       w0 = read_field('omega',x,y,t,slices=time, mesh=mesh, linear=linear,  $
                       filename=filename, points=pts, /equilibrium, $
                       rrange=xrange, zrange=yrange)

       psi = read_field('psi',x,y,t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       f = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x,y,t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
;        w = read_field('omega', x,y,t, slices=time, mesh=mesh, linear=linear,$
;                       filename=filename, points=pts, complex=complex, $
;                       rrange=xrange, zrange=yrange)
;        chi = read_field('chi',x,y,t, slices=time, mesh=mesh, linear=linear,$
;                         filename=filename, points=pts, complex=complex, $
;                         rrange=xrange, zrange=yrange)
       w_r = read_field('omega', x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange,op=2)
       w_z = read_field('omega', x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange,op=3)
       chi_lp = read_field('chi',x,y,t,slices=time, mesh=mesh, linear=linear,$
                           filename=filename, points=pts, complex=complex, $
                           rrange=xrange, zrange=yrange,op=7)
       chi_r = read_field('chi',x,y,t, slices=time, mesh=mesh, linear=linear,$
                          filename=filename, points=pts, complex=complex, $
                          rrange=xrange, zrange=yrange,op=2)
       chi_z = read_field('chi',x,y,t, slices=time, mesh=mesh, linear=linear,$
                          filename=filename, points=pts, complex=complex, $
                          rrange=xrange, zrange=yrange,op=3)

       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)

        data = -(i0*(chi_lp - chi_r/r)/r^2 $
                 + (i0_r*chi_r + i0_z*chi_z)/r^2 - 2.*i0*chi_r/r^3 $
                 + r^2*w0*grad_shafranov(rfac*f,x,y,tor=itor) $
                 + s_bracket(w0*r^2,rfac*f,x,y))/r $
          - a_bracket(i0,u,x,y) $
          + a_bracket(w0,psi,x,y) $
          + (w_z*psi0_r - w_r*psi0_z)

;         data = -(i0*grad_shafranov(chi,x,y,tor=itor)/r^2 $
;                  + s_bracket(i0/r^2,chi,x,y) $
;                  + r^2*w0*grad_shafranov(rfac*f,x,y,tor=itor) $
;                  + s_bracket(w0*r^2,rfac*f,x,y))/r $
;           - a_bracket(i0,u,x,y) $
;           + a_bracket(w0,psi,x,y) $
;           + a_bracket(w,psi0,x,y)


;       data = a_bracket(psi0, r*a_bracket(u, psi0, x, y) $
;                        - s_bracket(chi, psi0, x, y)/r^2, x, y)/r $
;         - w0*a_bracket(psi0, rfac*psi, x, y)/r $
;         + i0*a_bracket(psi0, rfac*u, x, y)/r $
;         + w0*s_bracket(psi0, rfac^2*f, x, y) $
;         + i0*s_bracket(psi0, rfac*chi, x, y)/r^4
;       data = data / sqrt(s_bracket(psi0,psi0,x,y))
       d = dimensions(/b0, t0=-1)
       symbol = '!6Curl[VxB]!X'

   ;===========================================
   ; curletaj
   ;===========================================
   endif else if(strcmp('curletaj', name, /fold_case) eq 1) then begin

       eta = read_field('eta',x,y,t,slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange)

       psi = read_field('psi',x,y,t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       f = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)
       

       data = - s_bracket(eta, i + rfac^2*f, x, y)/r $
         - eta*grad_shafranov(i + rfac^2*f, x, y, tor=itor)/r $
         + a_bracket(eta/r^2, rfac*psi, x, y)

       d = dimensions(/b0, t0=-1)
       symbol = '!6Curl[eta J]!X'


   ;===========================================
   ; torque
   ;===========================================
   endif else if(strcmp('torque', name, /fold_case) eq 1) then begin

       force_phi = read_field('force_phi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                              filename=filename, points=pts, complex=complex, $
                              rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       
       data = force_phi*r

       d = dimensions(/p0)
       symbol = '!6Beam Torque!X'


   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_b2', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=2)
       psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=3)
       i_r = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=2)
       i_z = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=3)
       f_r = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=2)
       f_z = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=3)
       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)

;          data = -(psi_r*conj(rfac*psi_r) + psi_z*conj(rfac*psi_z))/r^2 $
;            -     (conj(psi_r)*(rfac*psi_r) + conj(psi_z)*(rfac*psi_z))/r^2 $
;            + (conj(rfac*f_z)*(rfac*psi_r) - conj(rfac*f_r)*(rfac*psi_z))/r $
;            + ((rfac*f_z)*conj(rfac*psi_r) - (rfac*f_r)*conj(rfac*psi_z))/r $
;            + (conj(i_z+rfac^2*f_z)*psi_r - conj(i_r+rfac^2*f_r)*psi_z)/r $
;            + ((i_z+rfac^2*f_z)*conj(psi_r) - (i_r+rfac^2*f_r)*conj(psi_z))/r $
;            - (conj(i_r+rfac^2*f_r)*(rfac*f_r) + conj(i_z+rfac^2*f_z)*(rfac*f_z)) $
;            - ((i_r+rfac^2*f_r)*conj(rfac*f_r) + $
;              (i_z+rfac^2*f_z)*conj(rfac*f_z))

       data = (conj(i_z)*psi_r - conj(i_r)*psi_z)/r $
         - (conj(i_r)*(rfac*f_r) + conj(i_z)*(rfac*f_z)) $
         + (i_z*conj(psi_r) - i_r*conj(psi_z))/r $
         - (i_r*conj(rfac*f_r) + i_z*conj(rfac*f_z))

       data = data / 2. ; factor of 2 is from toroidal average
       
       d = dimensions(/p0)
       symbol = '!6Magnetic Torque!X'


   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_b1', name, /fold_case) eq 1) then begin

        psi0_r = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=2)
        psi0_z = read_field('psi', x, y, t, /equilibrium, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=3)
        i0_r = read_field('i', x, y, t, /equilibrium, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=2)
        i0_z = read_field('i', x, y, t, /equilibrium, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=3)
        psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=2)
        psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                           filename=filename, points=pts, complex=complex, $
                           rrange=xrange, zrange=yrange, op=3)
        i_r = read_field('i', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=2)
        i_z = read_field('i', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                           filename=filename, points=pts, complex=complex, $
                           rrange=xrange, zrange=yrange, op=3)
        f_r = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=2)
        f_z = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=3)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)
 
       data = -(psi0_r*(rfac*psi_r) + psi0_z*(rfac*psi_z))/r^2 $
         + ((i_z + rfac^2*f_z)*psi0_r - (i_r + rfac^2*f_r)*psi0_z)/r $
         + (i0_z*psi_r - i0_r*psi_z)/r $
         - (i0_r*rfac*f_r + i0_z*rfac*f_z)
       
       d = dimensions(/p0)
       symbol = '!7s!X'

   endif else if(strcmp('torque_p', name, /fold_case) eq 1) then begin

       p = read_field('p', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       
       rfac = complex(0., ntor)

       data = -rfac*p
       
       d = dimensions(/p0)
       symbol = '!7s!X'

   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_mu', name, /fold_case) eq 1) then begin

       mu = read_field('visc', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange,/equilibrium)
       mu_c = read_field('visc_c',x,y,t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange,/equilibrium)

       w = read_field('omega', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.

       data = mu*grad_shafranov(r^2*w,x,y,tor=itor) $
         + r^2*s_bracket(mu,w,x,y)

       if(ntor ne 0) then begin
           u = read_field('phi',x,y,t,slices=time, mesh=mesh, linear=linear, $
                          filename=filename, points=pts, complex=complex, $
                          rrange=xrange, zrange=yrange)
           chi = read_field('chi',x,y,t,slices=time, mesh=mesh, linear=linear,$
                            filename=filename, points=pts, complex=complex, $
                            rrange=xrange, zrange=yrange)
           rfac = complex(0.,1.)*ntor

           data = data + 2.*mu_c*rfac^2*w $
             - 4.*rfac*(mu_c*dz(u, y)) $
             - 2.*rfac*(mu-mu_c)*grad_shafranov(chi,x,y,tor=itor) $
             + mu*rfac*laplacian(chi,x,y,tor=itor)/r^2 $
             + a_bracket(mu, rfac*chi, x, y)/r^2
       endif
       
       d = dimensions(/p0)
       symbol = '!7s!X'

   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_v1', name, /fold_case) eq 1) then begin

       w0 = read_field('omega', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       n0 = read_field('den', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)

       u = read_field('phi',x,y,t,slices=time, mesh=mesh, linear=linear, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       w = read_field('omega',x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0.,1.)*ntor
       
       data = -a_bracket(r^4*n0*w0,u,x,y)/r - 2.*r^2*n0*w0*rfac*w
       d = dimensions(/p0)
       symbol = '!7s!X'

   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_vv2', name, /fold_case) eq 1) then begin

       psi0 = read_field('psi', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange,/equilibrium)
       n0 = read_field('den', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, /equilibrium)


       u = read_field('phi',x,y,t,slices=time, mesh=mesh, linear=linear, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       w = read_field('omega', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi',x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       
       data = -n0*conj(w)* $
         (r^3*a_bracket(psi0,u,x,y) + s_bracket(psi0,chi,x,y)) $
         / sqrt(s_bracket(psi0,psi0,x,y))
       d = dimensions(/p0, l0=1)
       symbol = '!6Angular Momentum Flux!X'

   ;===========================================
   ; Cole NTV
   ;===========================================
   endif else if(strcmp('cole_ntv', name, /fold_case) eq 1) then begin
       psi0 = read_field('psi',x,y,t,filename=filename,points=pts,$
                         /equilibrium,_EXTRA=extra)
       i0 = read_field('I',x,y,t,filename=filename,points=pts,$
                         /equilibrium,_EXTRA=extra)
       w0 = read_field('omega',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       p0 = read_field('p',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       pe0 = read_field('pe',x,y,t,filename=filename,points=pts,$
                        /equilibrium,_EXTRA=extra)
       n0 = read_field('den',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       psi = read_field('psi',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)
       i = read_field('I',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)
       f = read_field('f',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)
       Te1 = read_field('Te',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)


       db = read_parameter('db',filename=filename)
       zeff = read_parameter('zeff',filename=filename)
       ntor = read_parameter('ntor',filename=filename)
       ion_mass = read_parameter('ion_mass',filename=filename)
       n0_norm = read_parameter('n0_norm',filename=filename)
       B0_norm = read_parameter('B0_norm',filename=filename)
       L0_norm = read_parameter('l0_norm',filename=filename)
       v0_norm = B0_norm/sqrt(4.*!pi*n0_norm*ion_mass*1.6726e-24)
       t0_norm = L0_norm/v0_norm
       lambda = 16.              ; coulomb logarithm

       print, 'n0, B0, L0, v0, t0 = ', n0_norm, B0_norm, L0_norm, v0_norm, t0_norm

       r = radius_matrix(x,y,t)
       gradpsi = sqrt(s_bracket(psi0,psi0,x,y))

       pi0 = p0-pe0
       Ti = pi0/n0 > 0.01
       Te = pe0/(zeff*n0) > 0.01

       tprime = abs(s_bracket(Te,psi0,x,y)) 
       xi = -Te1/tprime*gradpsi
       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then begin
           xi(where(abs(tprime) lt 1e-6 or psi0 lt psis)) = 0.
       endif else begin
           xi(where(abs(tprime) lt 1e-6 or psi0 gt psis)) = 0.
       endelse

       forward_function flux_average_field
       fa = flux_average_field(r^2,psi0,x,y,t,r0=r0,flux=flux,$
                               filename=filename,integrate=0,_EXTRA=extra)
       r2_av = interpol(fa, flux, psi0)
       
       epsilon = read_field('minor radius',x,y,t,filename=filename,points=pts,$
                            /equilibrium,_EXTRA=extra)/r0

       w_B = abs(db*Ti * s_bracket(epsilon,psi0,x,y)/gradpsi^2)
       w_E = abs(w0 - db*s_bracket(pi0,psi0,x,y)/gradpsi^2 / n0)
       vti = sqrt(2.*Ti)
       
;        nu_i = (64.*!pi^(5/2)/3.)*(zeff*4.8032e-10*n0_norm)^4*lambda $
;          * l0_norm/(n0_norm*B0_norm^4) * n0/Ti^1.5
       nu_i = 1.988e-35*(n0/Ti^1.5)*lambda*zeff^4*n0_norm^3*l0_norm/B0_norm^4
       print, 'min, max(nu_i)', min(nu_i/t0_norm), max(nu_i/t0_norm)
       print, 'min, max(w_B)', min(w_B/t0_norm), max(w_B/t0_norm)
       print, 'min, max(w_E)', min(w_E/t0_norm), max(w_E/t0_norm)
       nu_eff = nu_i / (ntor*epsilon)
       
       mu_p = 0.21*ntor*vti^2*sqrt(epsilon*nu_eff) / $
         (r2_av*(w_E^1.5 + 0.3*w_B*sqrt(nu_eff) + 0.04*nu_eff^1.5))
       
       b02 = gradpsi^2/r^2 + i0^2/r^2
       
       if(icomplex eq 1) then begin
           b1b0 = (s_bracket(psi,psi0,x,y)/r^2 $
             + i*i0/r^2 $
             - complex(0.,ntor)*a_bracket(f,psi0,x,y)/r) / sqrt(b02)
           b1 = b1b0 + xi*s_bracket(sqrt(b02),psi0,x,y)/gradpsi
           b2 = b1*conj(b1) 
;            b2 = s_bracket(psi,conj(psi),x,y)/r^2 $
;              + i*conj(i)/r^2 $
;              + ntor^2*s_bracket(f,conj(f),x,y) $
;              - 0.5*complex(0.,ntor)*a_bracket(f,conj(psi),x,y)/r $
;              + 0.5*complex(0.,ntor)*a_bracket(conj(f),psi,x,y)/r
       end
       
       fa = flux_average_field(b2/b02,psi0,x,y,t,r0=r0,flux=flux,$
                               filename=filename,integrate=0,_EXTRA=extra)
       b2_av = interpol(fa, flux, psi0)
       
       data = -mu_p*b2_av*n0*r^2*w0
       d = dimensions(/p0)
       symbol = '!6NTV!X'

;       data = mu_p
;       d = dimensions(t0=-1)
       
;        data = nu_i
;        data = w_E
;       data = w_E^1.5 + 0.3*w_B*sqrt(nu_eff) + 0.04*nu_eff^1.5
;       d = dimensions(t0=-1)

   ;===========================================
   ; delta_B / B
   ;===========================================

   endif else if(strcmp('dB_B', name, /fold_case) eq 1) then begin
       psi0 = read_field('psi',x,y,t,filename=filename,points=pts,$
                         /equilibrium,_EXTRA=extra)
       i0 = read_field('I',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       Te0 = read_field('Te',x,y,t,filename=filename,points=pts,$
                        /equilibrium,_EXTRA=extra)
       psi = read_field('psi',x,y,t,filename=filename,points=pts,$
                        slice=time,/linear,complex=icomplex,_EXTRA=extra)
       i = read_field('I',x,y,t,filename=filename,points=pts,$
                      slice=time,/linear,complex=icomplex,_EXTRA=extra)
       f = read_field('f',x,y,t,filename=filename,points=pts,$
                      slice=time,/linear,complex=icomplex,_EXTRA=extra)
       Te1 = read_field('Te',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)

       ntor = read_parameter('ntor',filename=filename)

       r = radius_matrix(x,y,t)
 
       gradpsi = sqrt(s_bracket(psi0,psi0,x,y))

       tprime = abs(s_bracket(Te0,psi0,x,y)) 
       xi = -Te1/tprime*gradpsi
       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then begin
           xi(where(abs(tprime) lt 1e-6 or psi0 lt psis)) = 0.
       endif else begin
           xi(where(abs(tprime) lt 1e-6 or psi0 gt psis)) = 0.
       endelse

       b0 = sqrt(gradpsi^2 + i0^2)/r
       if(icomplex eq 1) then begin
 ;            b1b0 = complex(0.,ntor)*(s_bracket(psi,psi0,x,y)/r^2 $
 ;                    + i*i0/r^2 $
 ;                    - complex(0.,ntor)*a_bracket(f,psi0,x,y)/r) / b0
           b1b0 = (s_bracket(psi0,psi,x,y)/r^2 + i0*i/r^2 $
                   - complex(0.,ntor)*a_bracket(f,psi0,x,y)/r) / b0
           b1 = b1b0 + xi*s_bracket(b0,psi0,x,y)/gradpsi
           b2 = b1*conj(b1) 
       end
       data = b2/b0^2
       
       d = dimensions()
       symbol = '!7d!8B!6!U2!N!3/!8B!6!U2!N!X'


   ;===========================================
   ; radial electric field
   ;===========================================
   endif else if(strcmp('e_r', name, /fold_case) eq 1) then begin

       u   = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       v   = read_field('v'  , x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       b   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       pe  = read_field('pe',  x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       n   = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       eta = read_field('eta', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
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
       print, 'composite field ', name, ' not found'
   endelse

   endelse

   ; scale by linfac
   if((ilin eq 1) and (n_elements(linfac) ne 0) and (time ne -1) and keyword_set(linear)) $
     then begin
       print, 'scaling data by ', linfac
       data = data*linfac
   end
   
   print, 'converting units, mks, cgs=', keyword_set(mks), keyword_set(cgs)

   ; convert to appropriate units
   d0 = d
   get_normalizations, filename=filename,b0=b0,n0=n0,l0=l0,zeff=zeff,ion=mi
   convert_units, data, d0, b0, n0, l0, zeff, mi, cgs=cgs, mks=mks
   convert_units, x, dimensions(/l0), b0, n0, l0, zeff, mi, cgs=cgs, mks=mks
   convert_units, y, dimensions(/l0), b0, n0, l0, zeff, mi, cgs=cgs, mks=mks
   convert_units, realtime, dimensions(/t0), $
     b0, n0, l0, zeff, mi, cgs=cgs, mks=mks
   units = parse_units(d0, cgs=cgs, mks=mks)

   if(n_elements(h_symmetry) eq 1) then begin
       data = (data + h_symmetry*reverse(data, 2)) / 2.
   endif
   if(n_elements(v_symmetry) eq 1) then begin
       print, "v symmetry = ", v_symmetry
       data = (data + v_symmetry*reverse(data, 3)) / 2.
   endif

   ; apply mask
   if(n_elements(mask) ne 0 and n_elements(edge_val) ne 0) then begin
       data = data * (1. - mask) + mask*edge_val
   end

   if(keyword_set(complex)) then begin
       if(keyword_set(abs)) then begin
           print, 'Taking absolute value of data'
           data = abs(data)
       endif else if(keyword_set(phase)) then begin
           print, 'Taking phase of data'
           data = atan(imaginary(data),real_part(data))
       endif
   end

   if(n_elements(fac) ne 0) then begin
       print, 'applying factor = ', fac
       data = data*fac
   end


   ; perform flux-average
   if(keyword_set(flux_av)) then begin
       forward_function flux_average
       fa = flux_average(data, psi=psi, x=x, z=z, t=t, flux=flux, $
                         filename=filename, _EXTRA=extra)
       data = interpol(fa, flux, psi)
   end

   print, 'Done reading field'

   return, data
end

function read_field_3d, name, phi, x, z, t, points=points, tpoints=tpoints, $
                        symbol=symbol, units=u, ntor=ntor, _EXTRA=extra

  if(n_elements(points) eq 0) then points=200
  if(n_elements(tpoints) eq 0) then tpoints=16

  field = fltarr(tpoints, points, points)

  phi = fltarr(tpoints)

  for i=0, tpoints-2 do begin
     phi[i] = 360.*i/(tpoints-1)
     field[i,*,*] = reform(read_field(name,x,z,t,phi=phi[i],points=points,$
                                      symbol=symbol, units=u, _EXTRA=extra))
  end
  phi[tpoints-1] = 360.
  field[tpoints-1,*,*] = field[0,*,*]

  if(n_elements(ntor) eq 0) then begin
     return, field
  endif else begin
     fftfield = complexarr(tpoints,points,points)
     for i=0, points-1 do begin
        for j=0, points-1 do begin
           fftfield[*,i,j] = fft(field[*,i,j])
        end
     end
     lastfield = complexarr(1,points,points)
     lastfield = fftfield[ntor, *, *]
     return, lastfield
  end
end


; field_at_point
function field_at_point, field, x, z, x0, z0
   i = (n_elements(x)-1)*(x0-min(x)) / (max(x)-min(x))
   j = (n_elements(z)-1)*(z0-min(z)) / (max(z)-min(z))
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


function in_plasma, xy, psi, x, z, axis, psilim
   psi0 = field_at_point(psi,x,z,axis[0],axis[1])

   dpsidx = dx(psi,x)
   dpsidz = dz(psi,z)
   psin = field_at_point(psi,x,z,xy[0,*],xy[1,*])
   psinx = field_at_point(dpsidx,x,z,xy[0,*],xy[1,*])
   psinz = field_at_point(dpsidz,x,z,xy[0,*],xy[1,*])

   psin = (psin - psi0)/(psilim - psi0)
   psinx = psinx*(psilim - psi0)
   psinz = psinz*(psilim - psi0)
   dx = xy[0,*] - axis[0]
   dz = xy[1,*] - axis[1]

   return, ((psinx*dx + psinz*dz) gt 0.)
end

function path_at_flux, psi,x,z,t,flux,breaks=breaks,refine=refine,$
                       interval=interval, axis=axis, psilim=psilim, $
                       contiguous=contiguous, path_points=pts

   contour, psi[0,*,*], x, z, levels=flux, closed=0, $
     path_xy=xy, path_info=info, /path_data_coords, /overplot

   if(n_elements(xy) eq 0) then begin
       print, 'Error: no points at this flux value', flux
       return, 0
   end

   ; refine the surface found by contour with a single newton iteration
   if(keyword_set(refine)) then begin
       i = n_elements(x)*(xy[0,*]-min(x)) / (max(x)-min(x))
       j = n_elements(z)*(xy[1,*]-min(z)) / (max(z)-min(z))
       f = interpolate(reform(psi[0,*,*]),i,j)
       fx = interpolate(reform(dx(psi[0,*,*],x)),i,j)
       fz = interpolate(reform(dz(psi[0,*,*],z)),i,j)
       gf2 = fx^2 + fz^2
       l = (flux-f)/gf2
       xy[0,*] = xy[0,*] + l*fx
       xy[1,*] = xy[1,*] + l*fz
   endif

   ; remove points in private flux region
   if(n_elements(axis) gt 0 and n_elements(psilim) gt 0) then begin
       ip = in_plasma(xy,psi,x,z,axis,psilim)
       n_new = fix(total(ip))
       n_old = n_elements(xy[0,*])
       if(n_new gt 0) then begin
           xy_new = fltarr(2,n_new)
           j = 0
           for i=0, n_elements(xy[0,*])-1 do begin
               if(ip[i]) then begin
                   xy_new[*,j] = xy[*,i]
                   j = j + 1
               endif
           end
           xy = xy_new
       endif else begin
           print, 'excluding all points!'
           print, 'axis = ', axis
           print, 'psilim = ', psilim
       endelse
   endif
  
   ; find breaks
   n = n_elements(xy[0,*])
   dx = fltarr(n)
   for i=0, n-2 do begin
       dx[i] = sqrt((xy[0,i+1]-xy[0,i])^2 + (xy[1,i+1]-xy[1,i])^2)
   end
   dx[n-1] = sqrt((xy[0,0]-xy[0,n-1])^2 + (xy[1,0]-xy[1,n-1])^2)
   minforbreak = 5.*median(dx)

   breaks = where(dx gt minforbreak,count)

   if(count gt 0) then begin
       if(breaks[0] ge 0 and keyword_set(contiguous)) then begin
           xy = xy[*,0:breaks[0]-1]
       end
   end

   if(n_elements(interval) ne 0) then begin
       if(interval eq 0) then begin
           spline_p, xy[0,*], xy[1,*], xp_new, zp_new
       endif else begin
           spline_p, xy[0,*], xy[1,*], xp_new, zp_new, interval=interval
       endelse
       xy = transpose([[xp_new],[zp_new]])
   endif

   if(n_elements(pts) ne 0) then begin
       ind = findgen(pts)*n_elements(xy[0,*])/pts
       oldxy = xy
       xy = fltarr(2,pts)
       xy[0,*] = interpolate(oldxy[0,*], ind)
       xy[1,*] = interpolate(oldxy[1,*], ind)
   end

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

    ifixedb = read_parameter('ifixedb', _EXTRA=extra)
    if(ifixedb eq 1) then psilim = 0.

    print, 'LCFS: '
    print, ' Magnetic axis found at ', axis
    print, ' Active x-point at ', xpoint
    print, ' psi_0, psi_s = ', flux0, psilim
    
    return, psilim
end




pro plot_flux_contour, fval, _EXTRA=extra
   n = n_elements(fval)

   psi = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)

   loadct, 12
   contour, psi, x, z, closed=0, levels=fval(sort(fval)), color=color(6,10), $
     _EXTRA=extra
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
pro plot_lcfs, psi, x, z, psival=psival, _EXTRA=extra

    if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) then begin
        print, 'reading psi, plot_lcfs'
        psi = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
        print, min(psi), max(psi)
    end

    ; if psival not passed, choose limiter value
    if(n_elements(psival) eq 0) then $
      psival = lcfs(psi,x,z,_EXTRA=extra)

    ; plot contour
    loadct, 12
;    contour, psi, x, z, /overplot, nlevels=1, levels=psival, $
;      thick=!p.charthick*2., color=color(6,10)
    xy = path_at_flux(psi, x, z, t, psival, breaks=breaks)

    if(n_elements(breaks) eq 0) then begin
        oplot, xy[0,*], xy[1,*], thick=!p.thick, color=color(6,10)
    endif else begin
        breaks = [-1,breaks]

        for i=0, n_elements(breaks)-2 do begin
            oplot, xy[0,breaks[i]+1:breaks[i+1]], xy[1,breaks[i]+1:breaks[i+1]], $
              thick=!p.thick, color=color(6,10)
        end
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

   beta_t = 2.*s.Ave_P._data/bt0^2
   
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

   beta_t = 2.*s.Ave_P._data/bt0^2
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
     if (strcmp("time step", scalarname, /fold_case) eq 1) or $
        (strcmp("dt",scalarname, /fold_case) eq 1) then begin
       data = abs(s.dt._data)
       title = 'Time Step'
       symbol = translate('dt')
       d = dimensions(t0=1, _EXTRA=extra)
   endif else $
     if (strcmp("psimin", scalarname, /fold_case) eq 1) then begin
       data = s.psimin._data
       title = 'Psimin'
       symbol = translate('psim')
       d = dimensions(/b0, l0=2, _EXTRA=extra)
   endif else $
     if (strcmp("psibound", scalarname, /fold_case) eq 1) or $
        (strcmp("psilim",scalarname, /fold_case) eq 1) then begin
       data = s.psi_lcfs._data
       title = 'Psilim'
       symbol = translate('psil')
       d = dimensions(/b0, l0=2, _EXTRA=extra)
   endif else $
     if (strcmp("loop voltage", scalarname, /fold_case) eq 1) or $
     (strcmp("vl", scalarname, /fold_case) eq 1) then begin
       data = s.loop_voltage._data
       title = 'Loop Voltage'
       symbol = '!8V!DL!N!X'
       d = dimensions(/e0,/l0, _EXTRA=extra)
   endif else $
     if (strcmp("pellet rate", scalarname, /fold_case) eq 1) or $
     (strcmp("pelr", scalarname, /fold_case) eq 1) then begin
       data = s.pellet_rate._data
       title = 'Pellet Rate'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, t0=-1, _EXTRA=extra)
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
       data = s.E_KP._data + s.E_KT._data + s.E_K3._data
       title = 'Kinetic Energy'
       symbol = '!8KE!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("magnetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("me", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_MP._data + s.E_MT._data 
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
       s = read_scalars(filename=filename)
       n = tag_names(s)
       match = where(strcmp(n, scalarname, /fold_case) eq 1,count)
       if(count eq 0) then begin
           print, 'Scalar ', scalarname, ' not recognized.'
           return, 0
       endif
       data = s.(match[0])._data
       title = ''
       symbol = scalarname
       d = dimensions()
   endelse
   
   if(keyword_set(final)) then begin
       data = data[n_elements(data)-1]
   endif

   get_normalizations, b0=b0,n0=n0,l0=l0, zeff=zeff, ion_mass=mi, $
     filename=filename, _EXTRA=extra
   convert_units, data, d, b0, n0, l0, zeff, mi, _EXTRA=extra
   convert_units, time, dimensions(/t0), b0, n0, l0, zeff, mi, _EXTRA=extra
   units = parse_units(d, _EXTRA=extra)

   return, data
end

pro plot_scalar, scalarname, x, filename=filename, names=names, $
                 _EXTRA=extra, overplot=overplot, difference=diff, $
                 ylog=ylog, xlog=xlog, absolute_value=absolute, $
                 power_spectrum=power_spectrum, per_length=per_length, $
                 growth_rate=growth_rate, bw=bw, nolegend=nolegend, $
                 cgs=cgs,mks=mks,linestyle=ls, color=co, outfile=outfile

  if(n_elements(filename) eq 0) then filename='C1.h5'

  if(n_elements(names) eq 0) then names=filename

  nfiles = n_elements(filename)
  if(nfiles gt 1) then begin
      if(keyword_set(bw)) then begin
          if(n_elements(ls) eq 0) then ls = indgen(nfiles)
          if(n_elements(co) eq 0) then co = replicate(color(0,1),nfiles)
      endif else begin
          if(n_elements(ls) eq 0) then ls = replicate(0, nfiles)
          if(n_elements(co) eq 0) then co = shift(colors(),-1)
      endelse
      if(n_elements(x) eq 1) then x = replicate(x, nfiles)

      for i=0, nfiles-1 do begin
          if(n_elements(x) eq 0) then begin
              plot_scalar, scalarname, filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=co[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=power_spectrum, per_length=per_length, $
                growth_rate=growth_rate, linestyle=ls[i], nolegend=nolegend, $
                absolute_value=absolute,cgs=cgs,mks=mks,difference=diff
          endif else begin
              plot_scalar, scalarname, x[i], filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=co[i], _EXTRA=extra, ylog=ylog, xlog=xlog, $
                power_spectrum=power_spectrum, per_length=per_length, $
                growth_rate=growth_rate, nolegend=nolegend, $
                absolute_value=absolute,cgs=cgs,mks=mks,difference=diff
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
      n = min([n_elements(tdata), n_elements(data)])
      if(n_elements(data) lt n) then print, 'truncating data'
      if(n_elements(tdata) lt n) then print, 'truncating tdata'
      data = deriv(tdata(0:n-1), alog(abs(data(0:n-1))))
;      ytitle = '!7c !6(!7s!D!8A!N!6!U-1!N)!X'
      ytitle = make_label('!7c!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)
  endif

  if(keyword_set(diff)) then begin
      data = data - data[0]
  end

  if(keyword_set(absolute)) then data = abs(data)

  if(n_elements(x) eq 0) then begin   
      if(not keyword_set(overplot)) then begin
          plot, tdata, data, xtitle=xtitle, ytitle=ytitle, $
            title=title, _EXTRA=extra, ylog=ylog, xlog=xlog, $
            /nodata
      end     
      oplot, tdata, data, color=co, linestyle=ls, _EXTRA=extra
  endif else begin
      xi = x
      x = fltarr(1)
      z = fltarr(1)
      x[0] = xi
      z[0] = data[n_elements(data)-1]

      if(not keyword_set(overplot)) then begin
          plot, x, z, /nodata, $
            title=title, xtitle=xtitle, ytitle=ytitle, $
            _EXTRA=extra, ylog=ylog, xlog=xlog
      end
      oplot, x, z, color=co, linestyle=ls, _EXTRA=extra
  endelse
      if(n_elements(outfile) eq 1) then begin
         openw,ifile,outfile,/get_lun
        if(keyword_set(growth_rate)) then begin
         n = min([n_elements(tdata), n_elements(data)])
         printf,ifile,format='(2E16.6)',transpose([[tdata(1:n-1)],[data(1:n-1)]])
        endif else begin
         printf,ifile,format='(2E16.6)',transpose([[tdata],[data]])
        endelse
         free_lun, ifile
      endif
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

;....modified for new form of velocity
  vx = -dz(phi,z)*r
  vz =  dx(phi,x)*r

  chi = read_field('chi', x, z, t, points=points, _EXTRA=extra, slice=time)
  vx = vx + dx(chi,x)/r^2
  vz = vz + dz(chi,z)/r^2

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

  if(keyword_set(lcfs)) then plot_lcfs, points=200, slice=time, _EXTRA=extra
end


; ==============================================
; field_at_flux
;
; evaluates field at points where psi=flux
; ==============================================
function field_at_flux, field, psi, x, z, t, flux, theta=theta, angle=angle, $
                        xp=xp, zp=zp, integrate=integrate, axis=axis, $
                        refine=refine, interval=interval, $
                        psilim=psilim, contiguous=contiguous, _EXTRA=extra

   if(n_elements(field) le 1) then return, 0

   xy = path_at_flux(psi,x,z,t,flux,refine=refine, interval=interval, $
                    axis=axis,psilim=psilim,contiguous=contiguous)
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
function flux_coord_field, field, psi, x, z, t, slice=slice, area=area, i0=i0,$
                           fbins=fbins,  tbins=tbins, flux=flux, angle=angle, $
                           psirange=frange, nflux=nflux, qval=q, pest=pest, $
                           dV=dV, volume=volume, _EXTRA=extra, qflux=qflux

   if(n_elements(psi) eq 0) then begin
       linear = read_parameter('linear',_EXTRA=extra)
       if(linear eq 1) then begin
           psi = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
       endif else begin
           psi = read_field('psi',x,z,t,slice=slice,_EXTRA=extra)
       endelse
   endif

   sz = size(field)

   if(n_elements(fbins) eq 0) then fbins = sqrt(sz[2]*sz[3])
   if(n_elements(tbins) eq 0) then tbins = sqrt(sz[2]*sz[3])

   type = size(field, /type)
   if((type eq 6) or (type eq 9)) then begin
       result = complexarr(sz[1], fbins, tbins)
   endif else begin
       result = fltarr(sz[1], fbins, tbins)
   endelse
   flux = fltarr(sz[1], fbins)
   angle = fltarr(sz[1], tbins)
   area = fltarr(sz[1], fbins)
   volume = fltarr(sz[1], fbins)
   dV = fltarr(sz[1], fbins)

   psival = lcfs(psi,x,z,axis=axis,xpoint=xpoint,slice=slice,flux0=flux0, $
                 _EXTRA=extra)
   
   if(n_elements(range) eq 0) then begin
       ; if range not provided, use all flux within lcfs
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = [psival, flux0]
   endif else if(n_elements(range) eq 2) then begin
       oldrange = range
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = oldrange
   endif

   r = radius_matrix(x,z,t)

   bp = sqrt(s_bracket(psi,psi,x,z)/r^2)
   dpsidr = dx(psi,x)
   dpsidz = dz(psi,z)

   if(keyword_set(pest)) then begin
       linear = read_parameter('linear',_EXTRA=extra)
       if(n_elements(i0) le 1) then begin
           print, 'DBG: flux_coord_field reading field'
           if(linear eq 1) then begin
               i0 = read_field('i',x,z,t,slice=-1,_EXTRA=extra)
           endif else begin
               i0 = read_field('i',x,z,t,slice=slice,_EXTRA=extra)
           endelse
       endif
       bt = i0/r
       db = bt/(r*bp)
   endif

   print, 'binning with fbins, tbins=', fbins, tbins
   printed = 0
   left_handed = (flux0-psival)*mean(i0) lt 0 
   print, 'left_handed = ', left_handed

   q = fltarr(sz[1], fbins)

   for k=0, sz[1]-1 do begin

       angle[k,*] = 2.*!pi*findgen(tbins)/float(tbins) - !pi
       dpsi = float(range[k,1] - range[k,0])/float(fbins)

       for p=0, fbins-1 do begin
           flux[k,p] = range[k,1] - dpsi*(p+0.5)
       
           f = field_at_flux(field, psi, x, z, t, flux[k,p], $
                             angle=a, xp=xp, zp=zp, axis=axis, $
                             psilim=psival, /contiguous)

           if(n_elements(xp) le 2) then print, 'Too few points!'

           dx = deriv(xp)
           dz = deriv(zp)
           ds = sqrt(dx^2 + dz^2)
           area[k,p] = 2.*!pi*int_tabulated(findgen(n_elements(ds)),ds*xp)
           if(min(x) eq max(x)) then print, 'X ERROR!'
           if(min(z) eq max(z)) then print, 'Z ERROR!'
           ix = n_elements(x)*(xp - min(x))/(max(x) - min(x))
           iz = n_elements(z)*(zp - min(z))/(max(z) - min(z))
           h = interpolate(reform(bp[k,*,*]),ix,iz)
           dV[k,p] = 2.*!pi*int_tabulated(findgen(n_elements(ds)),ds/h)

           thimp = 0.
           if(keyword_set(pest)) then begin
               g = interpolate(reform(db[k,*,*]),ix,iz)

               ; determine position where angle changes sign
               if(p eq 0) then begin
                   ; if this is the first point, 
                   ; minimize the geometric angle
                   func = a
                   dum = min(func, i, /abs)
               endif else begin
                   ; otherwise, minimize the distance to the
                   ; point normal to the last surface
                   dp = [interpolate(reform(dpsidr[k,*,*]),ix,iz), $
                         interpolate(reform(dpsidz[k,*,*]),ix,iz)] $
                     / (xp*h) * (-dpsi*thimp)
                   func = (xp-dp[0]-p1[0])^2 + (zp-dp[1]-p1[1])^2
                   dum = min(func, i, /abs)
                   func = deriv(func)
               endelse
               da = deriv(func)
               if(da[i] eq 0) then print, 'DA ERROR!'
               index = i
               di = func[i]/da[i]
               if(abs(di) lt 0.5) then begin
                   index = index - di
               endif else begin
                   print, 'index correction too large: ', di, func[i], da[i]
               endelse
                   
               if(index lt 0 or index ge n_elements(xp)) then begin
                   print, 'Interpolation error ', index, p
                   index = i
               end

               p0 = [interpolate(xp,index), interpolate(zp,index)]
               ix0 = n_elements(x)*(p0[0] - min(x))/(max(x) - min(x))
               iz0 = n_elements(z)*(p0[1] - min(z))/(max(z) - min(z))

               dp = [interpolate(reform(dpsidr[k,*,*]),ix0,iz0), $
                     interpolate(reform(dpsidz[k,*,*]),ix0,iz0)] $
                 / (p0[0]*interpolate(reform(bp[k,*,*]),ix0,iz0))
               p1 = p0 + (-dpsi*(1.-thimp))*dp
;               print, 'p0, p1', p0, p1
                             
               ; calculate dt
               dt = ds*g
               if(dpsi lt 0) then dt = -dt

               pest_angle = fltarr(n_elements(dt))
               pest_angle[0] = 0.
               for j=1, n_elements(dt)-1 do begin
                   pest_angle[j] = pest_angle[j-1] + (dt[j]+dt[j-1])/2.
               end
             
               ; center pest_angle to change sign where a changes sign
               pest_angle = pest_angle - interpolate(pest_angle,index)

               ; rescale pest_angle
               q[k,p] = (max(pest_angle)-min(pest_angle)) $
                 /(2.*!pi)
               if(left_handed) then q[k,p] = -q[k,p]
               pest_angle = pest_angle/abs(q[k,p])
;               qval = interpol(q, qflux, flux[k,p])
;               pest_angle = pest_angle/qval
               ; constrain pest angle to +/- pi
               pest_angle = clamp_and_shift(pest_angle, shift=count)
               f = shift(f,-count)

               tol = 0.0001
               while(abs(pest_angle[0]-pest_angle[1]) lt tol) do begin
                   pest_angle = pest_angle[1:n_elements(pest_angle)-1]
                   f = f[1:n_elements(f)-1]
               end
               while(abs(pest_angle[n_elements(pest_angle)-2] - $
                      pest_angle[n_elements(pest_angle)-1]) lt tol) do begin
                   pest_angle = pest_angle[0:n_elements(pest_angle)-2]
                   f = f[0:n_elements(f)-2]
               end

               result[k,p,*] = interpol(f,pest_angle,angle[k,*])

;                problems = where(result[k,p,*] eq 1./0., count)
;                if(count gt 0) then begin
;                    print, "result = ", reform(result[k,p,*])
;                    print, "f = ", reform(f)
;                    print, "pest_angle = ", reform(pest_angle)
;                    stop
;                end

           endif else begin
               result[k,p,*] = interpol(f,a,angle[k,*])
           endelse 
       end

       if(dpsi gt 0) then begin
           dV[k,*] = -dV[k,*]
       end

       for i=1, n_elements(flux)-1 do $
         volume[k,i] = volume[k,i-1] $
         + (dV[k,i]+dV[k,i-1])/2. * (flux[k,i] - flux[k,i-1])
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
                             area=area, volume=volume, dV=dV, psirange=range, $
                             integrate=integrate, r0=r0, surface_weight=sw, $
                             nflux=nflux, _EXTRA=extra

   sz = size(field)

   points = sqrt(sz[2]*sz[3])

   if(n_elements(bins) eq 0) then bins = fix(points)

   print, 'flux averaging with ', bins, ' bins'

   type = size(field, /type)
   if(type eq 6) then begin
       result = complexarr(sz[1], bins)
   endif else begin
       result = fltarr(sz[1], bins)
   endelse
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
                               x, z, t, xp=xp, zp=zp, flux[k,p], /contiguous)
           bpf = field_at_flux(bp[k,*,*], psi[k,*,*], $
                               x, z, t, xp=xp, zp=zp, flux[k,p], /contiguous)
           if(dpsi gt 0.) then bpf = -bpf

           if(n_elements(bpf) lt 3) then continue

           dl = sqrt(deriv(xp)^2 + deriv(zp)^2)

           dV[k,p] = 2.*!pi*total(dl*xp/bpf)
           area[k,p] = 2.*!pi*total(dl*xp)
           if(keyword_set(integrate)) then begin
               result[k,p] = 2.*!pi*total(faf*xp*dl*(-dpsi)/bpf)
           end else begin
               if(keyword_set(sw)) then begin
                   ; weight by differential surface area
                   result[k,p] = total(faf*xp*dl)/total(xp*dl)
               endif else begin
                   ; weight by differential volume (standard definition)
                   result[k,p] = total(faf*xp*dl/bpf)/total(xp*dl/bpf)
               end
           endelse
       endfor
       if(keyword_set(integrate)) then begin
           val = reform(result[k,*])
           result[k,0] = val[0]/2.
           for i=1, n_elements(val)-1 do begin
               result[k,i] = result[k,i-1] + (val[i]+val[i-1])/2.
           end
       endif
   endfor

   nflux = (flux0-flux)/(flux0-psival)

   volume = fltarr(n_elements(flux))
   for i=1, n_elements(flux)-1 do $
     volume[i] = volume[i-1] + (dV[i]+dV[i-1])/2. * (flux[i] - flux[i-1])

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
function flux_average, field, psi=psi, i0=i0, x=x, z=z, t=t, r0=r0, $
                       flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
                       points=points, name=name, symbol=symbol, units=units, $
                       integrate=integrate, complex=complex, abs=abs, $
                       phase=phase, stotal=total, fac=fac, $
                       _EXTRA=extra

   type = size(field, /type)

   sz = size(field)

   if(n_elements(points) eq 0) then begin
       if(type ne 7) then points = sz[2] else points=200
   endif

   if((n_elements(psi) le 1) $
      or (n_elements(x) eq 0) or (n_elements(z) eq 0) $
      or (n_elements(t) eq 0)) then begin

       print, 'DBG: flux average reading field'
       psi = read_field('psi', x, z, t, points=points, $
                        mask=mask, /equilibrium, _EXTRA=extra, $
                       linear=0)

       sz = size(psi)
       if(n_elements(psi) le 1) then return, 0
   endif

   if(type eq 7) then begin ; named field
       if (strcmp(field, 'Safety Factor', /fold_case) eq 1) or $
         (strcmp(field, 'q', /fold_case) eq 1) then begin

           flux_t = flux_average('flux_t', psi=psi, i0=i0, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, _EXTRA=extra)
           
           units = ''
           name = '!6Safety Factor!X'
           symbol = '!8q!X'

           return, abs(deriv(flux, flux_t))/(2.*!pi)

       endif else $
         if(strcmp(field, 'alpha', /fold_case) eq 1) then begin

           p = flux_average('p', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, _EXTRA=extra)

           V = fltarr(n_elements(flux))
           for i=1, n_elements(flux)-1 do $
             V[i] = V[i-1] + (dV[i]+dV[i-1])/2. * (flux[i] - flux[i-1])
           pp = deriv(flux, p)

           units = ''
           name = '!8Ballooning Parameter!X'
           symbol = '!7a!X'

           return, -2.*dV/(2.*!pi)^2 * sqrt(abs(V)/(2.*!pi^2*R0)) * pp

       endif else $
         if(strcmp(field, 'shear', /fold_case) eq 1) then begin

           q = flux_average('q', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, _EXTRA=extra)

           V = fltarr(n_elements(flux))
           for i=1, n_elements(flux)-1 do $
             V[i] = V[i-1] + (dV[i]+dV[i-1])/2. * (flux[i] - flux[i-1])
           
           dqdV = deriv(V, q)

           units = ''
           name = '!6Magnetic Shear!X'
           symbol = '!8s!X'

           return, 2.*V*dqdV

       endif else $
         if(strcmp(field, 'dqdrho', /fold_case) eq 1) then begin

           q = flux_average('q', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, _EXTRA=extra)
           rho = flux_average('rho', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, _EXTRA=extra)
           
           units = make_units(l0=-1, _EXTRA=extra)
           name = '!8dq!3/!8d!7q!X'
           symbol = name

           return, deriv(rho, q)

       endif else $
         if(strcmp(field, 'lambda', /fold_case) eq 1) then begin
           ; this is m*lambda from Hegna, Callen Phys. Plasmas 1 (1994) p.2308

           I = read_field('I', x, z, t, points=points, /equilibrium, $
                          units=units, _EXTRA=extra)
           sigma = read_field('jpar_B', x, z, t, points=points, /equilibrium, $
                              units=units, _EXTRA=extra)
           sp = s_bracket(psi,sigma,x,z)/s_bracket(psi,psi,x,z)
           r = radius_matrix(x,z,t)

           q = flux_average('q', psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, /equilibrium, _EXTRA=extra)

           qp = deriv(flux,q)
           qrz = interpol(q, flux, psi)
           qprz = interpol(qp, flux, psi)
           jac = r^2*qrz/I

           dpsi = sqrt(s_bracket(psi,psi,x,z))
           dtheta = I/(r*qrz*dpsi)

           gpsi = flux_average(dpsi/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, _EXTRA=extra)
           gchi = flux_average(dtheta/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, _EXTRA=extra)
           gsp = flux_average(sp/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, _EXTRA=extra)
           gi = flux_average(I, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, _EXTRA=extra)
           
           units = ''
           name = '!8m!7k!X'
           symbol = name          

           return, -gi*q*gsp/(2.*qp)*(1./(gpsi*gchi))

       endif else $
         if(strcmp(field, 'rho', /fold_case) eq 1) then begin

           flux_t = flux_average('flux_t', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, _EXTRA=extra)
           bzero = read_parameter('bzero', _EXTRA=extra)
           print, 'bzero = ', bzero

           units = make_units(/l0, _EXTRA=extra)
           name = '!7q!X'
           symbol = name

           return, sqrt(flux_t/(!pi*bzero))

       endif else $
         if(strcmp(field, 'flux_t', /fold_case) eq 1) then begin
           print, 'DBG: flux_t reading field'

           if(n_elements(i0) le 1) then begin
               i0 = read_field('I', x, z, t, points=points, _EXTRA=extra)
           endif

           r = radius_matrix(x,z,t)
           field = i0/(2.*!pi*r^2)

           units = ''
           name = '!6Toroidal Flux!X'
           symbol = '!7w!D!8t!N!X'

           integrate = 1

       endif else $
         if(strcmp(field, 'beta_pol', /fold_case) eq 1) then begin
           I = read_field('I', x, z, t, points=points, _EXTRA=extra)
           
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
           
       endif else $
         if(strcmp(field, 'alpha', /fold_case) eq 1) then begin
           q = flux_average('q',psi=psi,x=x,z=z,nflux=nflux, $
                            flux=flux, bins=bins, r0=r0, $
                            points=points, last=last, _EXTRA=extra)
 
           beta = read_field('beta',x,z,t,points=points,last=last,_EXTRA=extra)
           r = read_field('r',x,z,t,points=points,last=last,_EXTRA=extra)

           betap = s_bracket(beta,r,x,z)
           alpha = $
             flux_average_field(betap, psi, x, z, t, r0=r0, flux=flux, $
                              nflux=nflux, area=area, dV=dV, bins=bins, $
                              integrate=integrate, _EXTRA=extra)

           symbol = '!7a!X'
           units = ''
           name = '!7a!X'

           return, -q^2*r0*alpha

       endif else $
         if(strcmp(field, 'kappa_implied', /fold_case) eq 1) then begin
           Q =  flux_average('heat_source', psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, _EXTRA=extra, /integrate)
           p = read_field('p',x,z,t,points=points,last=last,_EXTRA=extra)
           n = read_field('den',x,z,t,points=points,last=last,_EXTRA=extra)
           temp = p/n

           pprime = s_bracket(temp,psi,x,z)
           GradP = flux_average_field(pprime, psi, x, z, t, r0=r0, flux=flux, $
                                      nflux=nflux, area=area, dV=dV, $
                                      bins=bins, _EXTRA=extra)

           symbol = '!7j!X'
           d = dimensions(l0=2, t0=-1, n0=1)
           units = parse_units(d,_EXTRA=extra)
           name = '!7j!X'

           return, -Q/(dV*GradP)

       endif else begin
           field = read_field(field, x, z, t, points=points, complex=complex, $
                              symbol=symbol, units=units,dimensions=d,fac=fac,$
                              abs=abs, phase=phase, _EXTRA=extra)
           name = symbol
           if(keyword_set(total)) then begin
               d = d + dimensions(l0=2)
               units = parse_units(d, _EXTRA=extra)
           endif else if(keyword_set(integrate)) then begin
               d = d + dimensions(l0=3)
               units = parse_units(d, _EXTRA=extra)
           endif else begin
               symbol = '!12<!X' + symbol + '!12>!X'
           endelse
       endelse
   endif else begin
       name = ''
       symbol = ''
       units = make_units()
   endelse

   fa = flux_average_field(field, psi, x, z, t, r0=r0, flux=flux, $
                           nflux=nflux, area=area, dV=dV, bins=bins, $
                           integrate=integrate, surface_weight=total, $
                           _EXTRA=extra)

   if(keyword_set(total)) then begin
       fa = fa*area
   end

   return, fa
end


function flux_at_q, qval, normalized_flux=norm, points=pts, $
                    q=q, flux=flux, _EXTRA=extra
   q = flux_average('q', flux=flux, nflux=nflux, /equilibrium, points=pts, $
                    _EXTRA=extra)

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
                range=range, rrange=rrange, zrange=zrange, linear=linear, $
                xlim=xlim, cutx=cutx, cutz=cutz, mpeg=mpeg, linfac=linfac, $
                mask_val=mask_val, boundary=boundary, q_contours=q_contours, $
                overplot=overplot, phi=phi0, time=realtime, levels=levels, $
                phase=phase, abs=abs, operation=op, magcoord=magcoord, $
                outfile=outfile, fac=fac, _EXTRA=ex

   ; open mpeg object
   if(n_elements(mpeg) ne 0) then begin
       mpegid = mpeg_open([640,480],quality=100)

       for i=time[0],time[1] do begin
           plot_field, name, i, x, y, points=p, mesh=plotmesh, $
             mcolor=mc, lcfs=lcfs, title=title, units=units, $
             range=range, rrange=rrange, zrange=zrange, linear=linear, $
             xlim=xlim, cutx=cutx, cutz=cutz, linfac=linfac, levels=levels, $
             mask_val=mask_val, boundary=boundary, q_contours=q_contours, $
             overplot=overplot, phi=phi0, time=realtime, $
             phase=phase, abs=abs, operation=op, _EXTRA=ex

           image = tvrd(true=1)
               
           image[0,*,*] = rotate(reform(image[0,*,*]), 7)
           image[1,*,*] = rotate(reform(image[1,*,*]), 7)
           image[2,*,*] = rotate(reform(image[2,*,*]), 7)
               
           mpeg_put, mpegid, image=image, frame=(i-time[0])
       end

       print, 'Writing mpeg...'
       mpeg_save, mpegid, filename=mpeg
       mpeg_close, mpegid
       return
   end


   if(n_elements(time) eq 0) then time = 0
   if(n_elements(p) eq 0) then p = 200
   if(n_elements(title) eq 0) then notitle = 1 else notitle = 0

   if(keyword_set(phase) or keyword_set(abs)) then complex=1

   if(size(name, /type) eq 7) then begin
       field = read_field(name, x, y, t, slices=time, mesh=mesh, $
                          points=p, rrange=rrange, zrange=zrange, $
                          symbol=fieldname, units=u, linear=linear, $
                          mask=mask, phi=phi0, time=realtime, $
                          complex=complex, operation=op, $
                          linfac=linfac, fac=fac,_EXTRA=ex)
       if(n_elements(field) le 1) then return
       if(n_elements(units) eq 0) then units=u       
   endif else begin
       field = name
   endelse

   if(keyword_set(phase)) then begin
       field = atan(imaginary(field), real_part(field))*180./!pi
       units = '!6Degrees!X'
       fieldname = '!6Phase(!X' + fieldname + '!6)!X'
   endif else if(keyword_set(abs)) then begin
       field = abs(field)
   endif


   if(n_elements(realtime) ne 0) then print, 'time = ', realtime

   ; remove NaN's from result
   i = where(not float(finite(field)), count)
   if(count gt 0) then begin
       print, "removing NaN's"
       field[i] = 0.
   endif

   if(n_elements(mask_val) ne 0) then begin
       for k=0, n_elements(field[*,0,0])-1 do begin
           field[k,*,*] = field[k,*,*] - mask*(field[k,*,*] - mask_val)
       end
   endif

   field = real_part(field)

   if(n_elements(range) eq 0) then range = [min(field),max(field)]

   if((notitle eq 1) and (n_elements(t) ne 0)) then begin
       title = fieldname
   end

   if(n_elements(cutx) gt 0) then begin
       dum = min(x-cutx,i,/absolute)
       data = reform(field[0,i,*])
       if(keyword_set(overplot)) then begin
           oplot, y, data, _EXTRA=ex
       endif else plot, y, field[0,i,*], title=title, _EXTRA=ex
       if(n_elements(outfile) eq 1) then begin
           openw, ifile, outfile, /get_lun
           printf, ifile, format='(2E16.6)', transpose([[y], [data]])
           free_lun, ifile
       endif
   endif else if(n_elements(cutz) gt 0) then begin
       dum = min(y-cutz,i,/absolute)
       data = reform(field[0,*,i])
       if(keyword_set(overplot)) then begin
           oplot, x, data, _EXTRA=ex
       endif else plot, x, data, title=title, _EXTRA=ex
       if(n_elements(outfile) eq 1) then begin
           openw, ifile, outfile, /get_lun
           printf, ifile, format='(2E16.6)', transpose([[x], [data]])
           free_lun, ifile
       endif
   endif else if(keyword_set(magcoord)) then begin

       psi = read_field('psi',x,y,t,points=p,/equilibrium,_EXTRA=ex)
       field = flux_coord_field(field, psi, x, y, t, fbins=p, tbins=p, $
                                nflux=nflux, angle=angle, bins=p)
       
       contour_and_legend, transpose(field[0,*,*], [0,2,1]),$
         angle*180./!pi, nflux, title=title, $
         label=units, levels=levels, $
         xtitle='!6Angle (Degrees)!X', $
         ytitle='!7W!X', $
         range=range, _EXTRA=ex

       
   endif else begin          
       contour_and_legend, field[0,*,*], x, y, title=title, $
         label=units, levels=levels, $
         xtitle=make_label('!8R!X', /l0, _EXTRA=ex), $
         ytitle=make_label('!8Z!X', /l0, _EXTRA=ex), $
         range=range, _EXTRA=ex

       if(n_elements(q_contours) ne 0) then begin
           fval = flux_at_q(q_contours,points=p,_EXTRA=ex)
           plot_flux_contour, fval, points=p, closed=0, /overplot, $
             thick=!p.thick/2., _EXTRA=ex
       endif

       if(keyword_set(lcfs)) then begin
           print, 'passing slice = ', time[0]
           plot_lcfs, points=p, slice=time[0], $
             last=last, _EXTRA=ex
       endif

       if(keyword_set(boundary)) then plotmesh=1
       if(keyword_set(plotmesh)) then begin
           plot_mesh, mesh=mesh, /oplot, $
             boundary=boundary, _EXTRA=ex
       endif
   endelse
end

pro plot_3d, name, xslice=xslice, yslice=yslice, zslice=zslice, $
             _EXTRA=extra

   field = read_field_3d(name,y,x,z,t,_EXTRA=extra)

   data = fltarr(n_elements(x), n_elements(y))
   surface, data, x, y, /nodata, /save, zrange=[min(z), max(z)], _EXTRA=extra


   isosurface, data, max(data)*0.05, v, p

   tv, polyshade(v, p, /t3d, shades=shades)

end

;======================================================
; plot_flux_average
; ~~~~~~~~~~~~~~~~~
;
; plots the flux average quantity "name" at a give time
;======================================================
pro plot_flux_average, field, time, filename=filename, complex=complex, $
                       color=colors, names=names, bins=bins, linear=linear, $
                       xlog=xlog, ylog=ylog, overplot=overplot, fac=fac, $
                       lcfs=lcfs, normalized_flux=norm, points=pts, $
                       linestyle=ls, linfac=linfac, $
                       minor_radius=minor_radius, smooth=sm, t=t, rms=rms, $
                       bw=bw, srnorm=srnorm, last=last, mks=mks, cgs=cgs, $
                       q_contours=q_contours, rho=rho, integrate=integrate, $
                       multiply_flux=multiply_flux, abs=abs, phase=phase, $
                       stotal=total, nolegend=nolegend, outfile=outfile, $
                       val_at_q=val_at_q, flux_at_q=flux_at_q, _EXTRA=extra

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(linfac) eq 0) then linfac = 1.

   if(n_elements(time) eq 0) then time=0
   if(keyword_set(last)) then $
     time = fix(read_parameter('ntime',filename=filename)-1)

   if(n_elements(multiply_flux) eq 0) then multiply_flux = 0.

   if(n_elements(field) gt 1) then begin
       if(keyword_set(bw)) then begin
           ls = indgen(nfiles)
           colors = replicate(color(0,1), n_elements(field))
       endif else begin
           if(n_elements(colors) eq 0) then col = shift(colors(),-1)
           ls = replicate(0,n_elements(field))
       endelse
       if(n_elements(linfac) eq 1) then linfac=replicate(linfac, n_elements(field))
       for i=0, n_elements(field)-1 do begin
           if((n_elements(q_contours) ne 0) and (i eq 0)) then begin
               qcon = q_contours
           end
           plot_flux_average, field[i], time, filename=filename, $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=col[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, fac=fac, $
             linear=linear, multiply_flux=multiply_flux, mks=mks, cgs=cgs, $
             integrate=integrate, complex=complex, abs=abs, phase=phase, $
             stotal=total, q_contours=qcon, rho=rho, nolegend=nolegend, $
             linfac=linfac[i]
       end
       if(n_elements(names) ne 0) then begin
           plot_legend, names, colors=col, linestyle=ls, _EXTRA=extra
       end
       return
   end

   nfiles = n_elements(filename)
   if(nfiles gt 1) then begin
       if(n_elements(names) eq 0) then names=filename
       if(keyword_set(bw)) then begin
           if(n_elements(ls) eq 0) then ls = indgen(nfiles)
           colors = replicate(color(0,1), nfiles)
       endif else begin
           if(n_elements(colors) eq 0) then colors = shift(colors(),-1)
           if(n_elements(ls) eq 0) then ls = replicate(0,nfiles)
       endelse
       if(n_elements(time) eq 1) then time = replicate(time,nfiles)
       if(n_elements(linfac) eq 1) then linfac=replicate(linfac, nfiles)
       if(n_elements(multiply_flux) eq 1) then $
         multiply_flux = replicate(multiply_flux,nfiles)

       for i=0, nfiles-1 do begin
           newfield = field
           plot_flux_average, newfield, time[i], filename=filename[i], $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, fac=fac, $
             linear=linear, multiply_flux=multiply_flux[i], mks=mks, cgs=cgs, $
             integrate=integrate, complex=complex, abs=abs, phase=phase, $
             stotal=total, q_contours=q_contours, rho=rho, nolegend=nolegend, $
             linfac=linfac[i]
       end
       if(n_elements(names) gt 0 and not keyword_set(nolegend)) then begin
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
           if(n_elements(colors) eq 0) then colors = colors()
           if(time[0] gt 0) then colors = shift(colors,-1)
           ls = replicate(0,nt)
       endelse       
       if(n_elements(linfac) eq 1) then linfac=replicate(linfac, nt)
       for i=0, n_elements(time)-1 do begin
           newfield = field
           plot_flux_average, newfield, time[i], filename=filename, $
             overplot=((i gt 0) or keyword_set(overplot)), points=pts, $
             color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog, lcfs=lcfs, $
             normalized_flux=norm, minor_radius=minor_radius, smooth=sm, $
             t=t, rms=rms, linestyle=ls[i], srnorm=srnorm, bins=bins, fac=fac,$
             linear=linear, multiply_flux=multiply_flux, mks=mks, cgs=cgs, $
             integrate=integrate, complex=complex, asb=aba, phase=phase, $
             stotal=total, rho=rho, nolegend=nolegend, linfac=linfac[i], $
             q_contours=q_contours
           lab = parse_units(dimensions(/t0), cgs=cgs, mks=mks)
           get_normalizations, b0=b0, n0=n0, l0=l0, $
                        zeff=zeff, ion_mass=mi, filename=filename
           convert_units, t, dimensions(/t0), b0, n0, l0, zeff, mi, cgs=cgs, mks=mks
           names[i] = string(format='(%"!8t!6 = %g ",A,"!X")', t, lab)
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
                     psi=psi,x=x,z=z,t=t,nflux=nflux,linear=linear, fac=fac, $
                     mks=mks, cgs=cgs, area=area, integrate=integrate, linfac=linfac, $
                    complex=complex, abs=abs, phase=phase, stotal=total)

   if(n_elements(fa) le 1) then begin
       print, 'Error in flux_average. returning.'
       return
   endif

   ytitle = symbol
   if(strlen(units) gt 0) then ytitle = ytitle + '!6 ('+units+ '!6)!X'

   if(keyword_set(rms)) then begin
       fa2 = flux_average(field^2,flux=flux,nflux=nflux,points=pts,slice=time,$
                          filename=filename, t=t, linear=linear, fac=fac, $
                         mks=mks, cgs=cgs, complex=complex, $
                          abs=abs, phase=phase)
       fa = sqrt(1. - fa^2/fa2)

       ytitle = '!9S!6(1 - !12<' + symbol + '!12>!U2!n/!12<' + $
         symbol + '!6!U2!N!12>!6)!X'
       title = '!6Poloidal Deviation of ' + title + '!X'
   endif

   if(n_elements(multiply_flux) ne 0) then begin
       if(multiply_flux ne 0) then begin
           print, 'multiplying flux by', multiply_flux
           flux = flux*multiply_flux
           nflux=nflux*multiply_flux
       end
   end

   if(keyword_set(norm)) then begin
       flux = nflux
       xtitle = '!7W!X'
       lcfs_psi = 1.
   endif else if(keyword_set(srnorm)) then begin
       flux = sqrt(nflux)
       xtitle = '!9r!7W!X'
       lcfs_psi = 1.
   end else if(keyword_set(minor_radius)) then begin
       flux = flux_average('r',points=pts,file=filename,t=t,linear=linear,$
                    name=xtitle,bins=bins,units=units,slice=time, $
                          mks=mks, cgs=cgs, fac=fac)
       xtitle = '!12<' + xtitle + '!12> !6 ('+units+')!X'
   endif else if(keyword_set(rho)) then begin
       flux = flux_average('flux_t',points=pts,file=filename,/equilibrium,$
                           bins=bins,slice=time, fac=fac)
       flux = sqrt((flux - flux[0])/(flux[n_elements(flux)-1] - flux[0]))
       lcfs_psi = 1.
       xtitle = '!7q!X'
   endif

   if(n_elements(sm) eq 1) then begin
       fa = smooth(fa,sm)
   end

   if(keyword_set(overplot)) then begin
       oplot, flux, fa, color=colors, linestyle=ls, _EXTRA=extra
   endif else begin
       if(n_elements(colors) eq 0) then begin
           plot, flux, fa, xtitle=xtitle, linestyle=ls, $
             ytitle=ytitle, title=title, xlog=xlog, ylog=ylog, $
             _EXTRA=extra
       endif else begin
           plot, flux, fa, xtitle=xtitle, linestyle=ls, $
             ytitle=ytitle, title=title, xlog=xlog, ylog=ylog, /nodata, $
             color=color(0), _EXTRA=extra
           oplot, flux, fa, color=colors, linestyle=ls, _EXTRA=extra
       endelse
   endelse

   if(keyword_set(lcfs)) then begin
       oplot, [lcfs_psi,lcfs_psi], !y.crange, linestyle=2, color=colors
   endif

   if(n_elements(q_contours) ne 0) then begin
       fvals = flux_at_q(q_contours, points=pts, filename=filename, $
                         slice=time, normalized_flux=norm, bins=bins)
       for k=0, n_elements(fvals)-1 do begin
           oplot, [fvals[k], fvals[k]], !y.crange, linestyle=1, color=colors
       end
       result = interpol(fa, flux, fvals)
       print, 'filename = ', filename
       val_at_q = result
       flux_at_q = fvals
       print, "flux at q's =  ", flux_at_q
       print, "values at q's =  ", val_at_q
   endif

   if(n_elements(outfile) eq 1) then begin
       openw, ifile, outfile, /get_lun
       if(keyword_set(complex)) then begin
           printf, ifile, format='(3E16.6)', $
             transpose([[reform(flux[0,*])], $
                        [reform(real_part(fa[0,*]))], $
                        [reform(imaginary(fa[0,*]))]])
       endif else begin
           printf, ifile, format='(2E16.6)', $
             transpose([[reform(flux[0,*])], [reform(fa[0,*])]])
       endelse
       free_lun, ifile
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



pro write_geqdsk, eqfile=eqfile, $
                  psilim=psilim, points=pts, _EXTRA=extra
  
  if(n_elements(slice) eq 0) then begin
      slice = read_parameter('ntime', _EXTRA=extra) - 1
  end
  if(n_elements(eqfile) eq 0) then eqfile = 'geqdsk.out'
  get_normalizations, b0=b0, n0=n0, l0=l0, _EXTRA=extra

  bzero = read_parameter('bzero', _EXTRA=extra)
  rzero = read_parameter('rzero', _EXTRA=extra)

  ; calculate flux averages
  psi = read_field('psi',x,z,t,mesh=mesh,/equilibrium,_EXTRA=extra)
  psi_r = read_field('psi',x,z,t,mesh=mesh,/equilibrium,op=2,_EXTRA=extra)
  psi_z = read_field('psi',x,z,t,mesh=mesh,/equilibrium,op=3,_EXTRA=extra)
  psi_lp= read_field('psi',x,z,t,/equilibrium,op=7,_EXTRA=extra)
  p0 = read_field('p',x,z,t,/equilibrium,_EXTRA=extra)
  p0_r = read_field('p',x,z,t,/equilibrium,op=2,_EXTRA=extra)
  p0_z = read_field('p',x,z,t,/equilibrium,op=3,_EXTRA=extra)
  I0 = read_field('I',x,z,t,/equilibrium,_EXTRA=extra)
  I0_r = read_field('I',x,z,t,/equilibrium,op=2,_EXTRA=extra)
  I0_z = read_field('I',x,z,t,/equilibrium,op=3,_EXTRA=extra)
  r = radius_matrix(x,z,t)
  beta = r^2*2.*p0/(s_bracket(psi,psi,x,z) + I0^2)
  beta0 = mean(2.*p0*r^2/(bzero*rzero)^2)
  b2 = (s_bracket(psi,psi,x,z) + I0^2)/r^2
  dx = (max(x)-min(x))/(n_elements(x) - 1.)
  dz = (max(z)-min(z))/(n_elements(z) - 1.)
  jphi = psi_lp - psi_r/r
  tcur = total(jphi*dx*dz/r)
  print, 'current = ', tcur

  ; calculate magnetic axis and xpoint
  lcfs_psi = lcfs(psi,x,z, axis=axis, xpoint=xpoint, $
                  flux0=flux0, _EXTRA=extra)

  ; plot psi
  contour_and_legend, psi, x, z

  ; calculate wall points
  bound_xy = get_boundary_path(mesh=mesh, _EXTRA=extra)
  nwall = n_elements(reform(bound_xy[0,*]))
  rwall = fltarr(nwall)
  zwall = fltarr(nwall)
  rwall[*] = bound_xy[0,*]
  zwall[*] = bound_xy[1,*]

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
  oplot, rwall, zwall, psym=6
 
  r2bp = psi_r^2 + psi_z^2

  ; calculate flux averages
  p = flux_average_field(p0,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)
  pp = (p0_r*psi_r + p0_z*psi_z)/r2bp
  pprime = flux_average_field(pp,psi,x,z,nflux=nflux,flux=flux,_EXTRA=extra)
  I = flux_average_field(I0,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)
  ffp = I0*(I0_r*psi_r + I0_z*psi_z)/r2bp
  ffprim = flux_average_field(ffp,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)
  q = flux_average('q',slice=time,psi=psi,x=x,z=z,t=t,flux=flux,nflux=nflux,_EXTRA=extra)
;  q = smooth(q,5,/edge)
  jb = (I0_z*psi_r - I0_r*psi_z - jphi*I0)/r^2
  jdotb = flux_average_field(jb,psi,x,z,t,flux=flux,nflux=nflux,_EXTRA=extra)
  r2i = flux_average_field(1./r^2,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)

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

  psimin = flux0

  if(psimin gt psilim) then begin
      psi = -psi
      psilim = -psilim
      psimin = -psimin
      flux = -flux
      pprime = -pprime
      ffprim = -ffprim
  end

  name = ['M3DC1', '03/17/', '2013    ', $
          '#000000', '0000', '']
  idum = 3
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
  printf, file, format=f2020, reform(psi[0,*,*])
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

  ; jsolver wants boundary points in clockwise direction
  if(ifixedb eq 1) then begin
      print, 'Reversing boundary points'
      rlim = reverse(rlim)
      zlim = reverse(zlim)
  end

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

  window, 0

end

pro plot_field_3d, fieldname, contrast=contrast, flux=flux, $
                   light_brightness=light_brightness, points=points, $
                   brightness=brightness, specular=specular, $
                   solid=solid, positive_only=positive_only, $
                   power=power, angle=angle, ax=ax, az=az, value=value, $
                   xslice=xslice, zslice=zslice, phislice=phislice, $
                   absolute_value=absolute, range=range, $
                   noaxes=noaxes, reject=reject, tpoints=tpoints, $
                   xrange=xrange, yrange=yrange, zrange=zrange, $
                   _EXTRA=extra


  if(n_elements(points) eq 0) then points = 200
  if(n_elements(tpoints) eq 0) then tpoints = 16
  angles = tpoints

  if(keyword_set(solid)) then begin
      if(n_elements(flux) eq 0) then flux = 1.
      psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,points=points)
      psival = lcfs(psi0,x,z,flux0=flux0,_EXTRA=extra)
      flux1 = (psival-flux0)*flux + flux0
      if(n_elements(angle) eq 0) then angle = .75
  endif else begin
      if(n_elements(angle) eq 0) then angle = 1.
  endelse

  shown_angles = angles*angle

  data = fltarr(points,shown_angles,points)
  if(keyword_set(solid)) then $
    psi = fltarr(points,shown_angles,points)

  threed = read_parameter('3d', _EXTRA=extra)
  print, '3d = ', threed

  phi = fltarr(shown_angles)
  if(threed eq 1) then begin
      for i=0, shown_angles-1 do begin
          phi[i] = 360.*i/(angles-1.)
          if(i eq angles-1) then pp = 0 else pp = phi[i]
          data[*,i,*] = $
            read_field(fieldname,x,z,t,_EXTRA=extra,edge_val=0.,phi=pp, $
                      points=points)
          if(keyword_set(solid)) then psi[*,i,*] = psi0[0,*,*]
      end
  endif else begin
      ntor = read_parameter('ntor',_EXTRA=extra)
      print, 'ntor = ', ntor
      
      field = read_field(fieldname,x,z,t,/complex,_EXTRA=extra,edge_val=0., $
                        points=points, /linear)
      field0 = read_field(fieldname,x,z,t,_EXTRA=extra,edge_val=0., $
                        points=points, /equilibrium)

      for i = 0, shown_angles-1 do begin
          phi[i] = 2.*!pi*i/(angles-1.)
          data[*,i,*] = real_part(field[0,*,*]* $
                                  (cos(ntor*phi[i]) + $
                                   complex(0.,1.)*sin(ntor*phi[i]))) $
            + field0[0,*,*]
          
          if(keyword_set(solid)) then psi[*,i,*] = psi0[0,*,*]
      end
      phi = phi*180./!pi
  endelse

  if(keyword_set(positive_only)) then begin
      data = data > 0.
  endif else if(keyword_set(absolute)) then begin
      data = abs(data)
  endif

  if(n_elements(power) ne 0.) then data = data^power

  itor = read_parameter('itor', _EXTRA=extra)
  if(itor eq 1) then begin
      if(n_elements(xrange) ne 2) then xrange = [-max(x), max(x)]
      if(n_elements(yrange) ne 2) then yrange = xrange
      if(n_elements(zrange) ne 2) then zrange = [min(z), max(z)]
      xtitle='x'
      ytitle='y'
  endif else begin
      if(n_elements(xrange) ne 2) then xrange = [min(x), max(x)]
      if(n_elements(yrange) ne 2) then yrange = [min(phi), max(phi)]
      if(n_elements(zrange) ne 2) then zrange = [min(z), max(z)]
      xtitle='R'
      ytitle='!9P!X'
  endelse

  tmpdat = fltarr(2, 2)
  surface, tmpdat, xrange, yrange, charsize=2.5, $
    /nodata, /save, ax=ax, az=az, xstyle=1, ystyle=1, zstyle=1, $
    xrange=xrange, yrange=yrange, zrange=zrange, $
    xtitle=xtitle, ytitle=ytitle, ztitle='Z'
 

  if(keyword_set(solid)) then begin
      print, 'calling interval_volume'
      interval_volume, psi, flux0, flux1, v, conn, $
        auxdata_in=data, auxdata_out=aout
      print, 'calling tetra_surface'
      p = tetra_surface(v, conn)
      set_shading, reject=1
  
  endif else begin
      if(n_elements(xslice) eq 1) then begin
          if(n_elements(value) eq 0) then value = xslice $
          else value = [value,xslice]
          if(n_elements(normal) eq 0) then normal = [1,0,0] $
          else normal = [[normal], [1,0,0]]
      endif
      if(n_elements(zslice) eq 1) then begin
          if(n_elements(value) eq 0) then value = zslice $
          else value = [value,zslice]
          if(n_elements(normal) eq 0) then normal = [0,0,1] $
          else normal = [[normal], [0,0,1]]
      endif
      for i=0, n_elements(phislice)-1 do begin
         if(n_elements(value) eq 0) then value = phislice[i] $
         else value = [value,phislice[i]]
         if(n_elements(normal) eq 0) then normal = [0,1,0] $
         else normal = [[normal], [0,1,0]]
      end
      print, value
      print, normal
      if(n_elements(value) eq 0) then value=0.9
      plot_slice, data, x, phi, z, value=value, normal=normal, range=range, $
        itor=itor, reject=reject
      if(not keyword_set(noaxes)) then begin
          surface, tmpdat, xrange, yrange, xstyle=1, ystyle=1, zstyle=1, $
            /nodata, /noerase, ax=ax, az=az, charsize=2.5, $
            xrange=xrange, yrange=yrange, zrange=zrange, $
            xtitle=xtitle, ytitle=ytitle, ztitle='Z'
      end
      return
  endelse
   
  ; reduce values at edge
;  result = where(v[1,*] eq min(v[1,*]) or v[1,*] eq max(v[1,*]))
;  aout(result) = aout(result)/2.

  rr = (max(x) - min(x))*v[0,*]/(n_elements(x)-1.) + min(x)
  yy = (max(phi) - min(phi))*v[0,*]/(n_elements(phi)-1.) + min(phi)
  zz = (max(z) - min(z))*v[2,*]/(n_elements(z)-1.) + min(z)

  q = v

  if(itor eq 1) then begin
      q[0,*] = rr*cos(yy)
      q[1,*] = rr*sin(yy)
      q[2,*] = zz
  endif else begin
      q[0,*] = rr
      q[1,*] = yy
      q[2,*] = zz
  endelse

  if(n_elements(az) eq 0) then begin
      if(keyword_set(solid)) then az=-25 else az=0.
  endif
  if(n_elements(ax) eq 0) then begin
      if(keyword_set(solid)) then ax=35 else ax=15.
  endif

  if(n_elements(brightness) eq 0) then brightness=0.15
  if(n_elements(contrast) eq 0) then contrast=1.
  if(n_elements(light_brightness) eq 0) then light_brightness=0.3
  if(n_elements(specular) eq 0) then specular=0.0

  if(keyword_set(solid)) then begin
      light_dir = [0., 1., -1.]/sqrt(2.)
      normals = compute_mesh_normals(q, p)
      reflect = -(light_dir[0]*normals[0,*] + $
                  light_dir[1]*normals[1,*] + $
                  light_dir[2]*normals[2,*]) > 0
      shades = bytscl(aout) * contrast + brightness*255 $
        + light_brightness*reflect*255 $
        + specular*exp((reflect-1.)*5.)*255< 255
  endif else begin
      shades = bytscl(aout)
  endelse


  print, 'plotting'
  tv, polyshade(q, p, /t3d, shades=shades)


  
  surface, tmpdat, xrange, yrange, xstyle=1, ystyle=1, zstyle=1, $
    /nodata, /noerase, ax=ax, az=az, charsize=2.5, $
    xrange=xrange, yrange=yrange, zrange=zrange, $
    xtitle='R', ytitle='!9P!X', ztitle='Z'

end

pro movie_3d, fieldname, mpeg=mpeg, _EXTRA=extra
   if(n_elements(mpeg) eq 0) then mpeg = fieldname + '.mpeg'
   ntime = read_parameter('ntime', _EXTRA=extra)
   print, 'ntime = ', ntime

   mpegid = mpeg_open([640,480],bitrate=104857200, iframe_gap=4)

   for i=0, ntime-1 do begin
       plot_field_3d, fieldname, slice=i, _EXTRA=extra
       image = tvrd(true=1)
               
       image[0,*,*] = rotate(reform(image[0,*,*]), 7)
       image[1,*,*] = rotate(reform(image[1,*,*]), 7)
       image[2,*,*] = rotate(reform(image[2,*,*]), 7)

       mpeg_put, mpegid, image=image, frame=5*i
   end

   print, 'Writing ', mpeg, '...'
   mpeg_save, mpegid, filename=mpeg
   mpeg_close, mpegid
   print, 'Done.'
end

pro test_mesh, filename, nplanes=nplanes, _EXTRA=extra
  mesh = read_mesh(filename=filename, _EXTRA=extra)
  nelms = mesh.nelms._data

  if(n_elements(nplanes) eq 0) then begin
     nplanes = read_parameter('nplanes', filename=filename)
     print, 'read nplanes.', nplanes
  end
  n = nelms/nplanes

  print, 'elms = ', nelms
  print, 'nplanes = ', nplanes
  print, 'elms_per_plane = ', n

  tol = 1e-6

  correct = long(0)
  wrong = lonarr(6)
  for i=long(0), n-1 do begin
     for k=1, nplanes-1 do begin
        for j=0, 5 do begin
           if(abs(mesh.elements._data[j,i]- $
                  mesh.elements._data[j,i+n*k]) gt tol) then begin
              wrong[j] = wrong[j] + 1
              print, 'misalignment of element ', i, ' at plane ', k
           endif else correct = correct +  1
        end
     end
  end
  for j=0, 5 do begin
     print, 'wrong ', j, ' = ', wrong[j]
  end
  print, 'correct = ', correct/6
end

pro plot_perturbed_surface, q, scalefac=scalefac, points=pts, $
                            filename=filename, phi=phi0, _EXTRA=extra
   icomp =read_parameter('icomplex', filename=filename)
   if(n_elements(scalefac) eq 0) then scalefac=1.
   if(n_elements(scalefac) eq 1) then $
     scalefac=replicate(scalefac, n_elements(q))

   psi0 = read_field('psi',x,z,t,filename=filename, /equilibrium, $
                     points=pts, _EXTRA=extra)
   psi0_r = read_field('psi',x,z,t,filename=filename, /equilibrium, $
                       points=pts, _EXTRA=extra, op=2)
   psi0_z = read_field('psi',x,z,t, filename=filename, /equilibrium,$
                       points=pts, _EXTRA=extra, op=3)
   zi = read_field('displacement',x,z,t,filename=filename, /linear, $
                   points=pts, slice=slice, complex=icomp, $
                   phi=phi0, _EXTRA=extra)
   
   if(n_elements(bins) eq 0) then bins = n_elements(x)
   fvals = flux_at_q(q, points=pts, filename=filename, $
                     normalized_flux=norm, bins=bins)
   print, fvals

   xhat = psi0_r/sqrt(psi0_r^2 + psi0_z^2)
   zhat = psi0_z/sqrt(psi0_r^2 + psi0_z^2)

   plot, x, z, /nodata, /iso, _EXTRA=extra, $
     xtitle='!8R!6 (m)!X', ytitle='!8Z!6 (m)!X'
   c = colors()
   c0 = c
   if(n_elements(fvals) eq 1) then begin
       l0 = 0
       c[1] = c0[3]
   endif else begin
       l0 = 1
   endelse
   for k=0, n_elements(fvals)-1 do begin
       xy = path_at_flux(psi0, x, z, t, fvals[k], /contiguous)
       xy_new = xy

       for j=0, n_elements(xy[0,*])-1 do begin
           dx = field_at_point(xhat, x, z, xy[0,j], xy[1,j])
           dz = field_at_point(zhat, x, z, xy[0,j], xy[1,j])
           dr = field_at_point(zi, x, z, xy[0,j], xy[1,j])
           xy_new[0,j] = xy[0,j] + dr*dx*scalefac[k]
           xy_new[1,j] = xy[1,j] + dr*dz*scalefac[k]
       end

       oplot, xy[0,*], xy[1,*], linestyle=l0, color=c0[k+1]
       oplot, [xy[0,n_elements(xy[0,*])-1], xy[0,0]],  $
         [xy[1,n_elements(xy[0,*])-1], xy[1,0]], linestyle=l0, color=c0[k+1]

       oplot, xy_new[0,*], xy_new[1,*], color=c[k+1]
       oplot, [xy_new[0,n_elements(xy[0,*])-1], xy_new[0,0]], $
         [xy_new[1,n_elements(xy[0,*])-1], xy_new[1,0]], color=c[k+1]
   end
end

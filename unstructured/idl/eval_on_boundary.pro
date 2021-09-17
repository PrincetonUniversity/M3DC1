function eval_on_boundary, name, time, r, z, operation=op, filename=filename, phi=phi0, _EXTRA=ex

   if(n_elements(phi0) eq 0) then phi0 = 0.
   if(n_elements(filename) eq 0) then filename='C1.h5'
   file_id = h5f_open(filename)
   time_group_id = h5g_open(file_id, time_name(time))
   field_group_id = h5g_open(time_group_id, 'fields')
   field = h5_parse(field_group_id, name, /read_data)
   field = field._data
   mesh = h5_parse(time_group_id, 'mesh', /read_data)
   
   xy = get_boundary_path(norm=norm, center=center, angle=angle, $
                          length=length, filename=filename, _EXTRA=ex)
   r = xy[0,*]
   z = xy[1,*]
   N = n_elements(r)

   nelms = mesh.nelms._data 

   if(n_elements(op) eq 0) then op = 1
   
   pos = fltarr(2)
   localpos = fltarr(2)
   index = intarr(2)

   version = read_parameter('version', filename=filename)
   print, 'Version = ', version

   elm_data = mesh.elements._data
   sz = size(elm_data, /dimensions)
   if(sz[0] ge 10 and version lt 15) then begin
      print, 'Size = ', sz[0]
      print, 'Warning: mismatch between data size and version type.'
      print, 'Assuming version 16'
      version = 16
   end
   
   if(version eq 0) then begin
       xzero = read_parameter("xzero", filename=filename)
       zzero = read_parameter("zzero", filename=filename)

       nonrect = read_parameter('nonrect', filename=filename)
       if(nonrect eq 0.) then begin
           xmin = min(elm_data[4,*])
           ymin = min(elm_data[5,*])
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

   ; clamp phi0 to period
   if(version ge 11) then begin
      period = mesh.period._data
   endif else begin
      itor = read_parameter('itor', filename=filename)
      if(itor eq 1) then begin
         period = 2.*!pi
      endif else begin
         rzero = read_parameter('rzero', filename=filename)
         period = 2.*!pi*rzero
      end
   end
   phi0 = phi0 - floor(phi0/period)*period
   
   sz = size(elm_data, /dim)
   if(sz[0] gt 8) then begin
     threed = 1
   endif else begin
     threed = 0
   endelse

   result = fltarr(N)
   
   eps = 1e-6
   localphi = 0.

   if(version lt 15) then begin
      ib = 6
      izone = 1
   endif else begin
      ib = 7
   end

   ; for each triangle, evaluate points within triangle which fall on
   ; rectilinear output grid.  This is MUCH faster than looping through
   ; rectilinear grid points and searching for the appropriate triangles
   for i=long(0),nelms-1 do begin
      i_data = elm_data[*,i]
      if(threed eq 1) then begin
         d = i_data[ib+1]
         phi = i_data[ib+2]
         localphi = phi0 - phi
         if(localphi lt 0 or localphi gt d) then begin
            ; this assumes elements are ordered by plane
            i = i + nelms / mesh.nplanes._data - 1
            continue
         end
      endif
      a = i_data[0]
      b = i_data[1]
      c = i_data[2]
      t = i_data[3]
      x = i_data[4]
      y = i_data[5]
      if(version ge 16) then izone = i_data[ib]

      co = cos(t)
      sn = sin(t)

      p1 = [x, y]
      p2 = p1 + [(b+a)*co, (b+a)*sn]
      p3 = p1 + [b*co - c*sn, b*sn + c*co]

      minpos = [min([p1[0], p2[0], p3[0]])-eps, min([p1[1], p2[1], p3[1]])-eps]
      maxpos = [max([p1[0], p2[0], p3[0]])+eps, max([p1[1], p2[1], p3[1]])+eps]
             
      for j=0,N-1 do begin
         pos = xy[*,j]
         
         if ((pos[0] lt minpos[0]) or (pos[0] gt maxpos[0]) or $
             (pos[1] lt minpos[1]) or (pos[1] gt maxpos[1])) then continue
         
         localpos = [ (pos[0]-x)*co + (pos[1]-y)*sn - b, $
                     -(pos[0]-x)*sn + (pos[1]-y)*co, localphi]

         if(is_in_tri(localpos,a,b,c) eq 1) then result[j] = eval(field, localpos, t, i, op=op)
      end
   end

   return, result
end

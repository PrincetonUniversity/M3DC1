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

   ; find minimum and maximum node coordinates
   if(n_elements(xrange) lt 2 or n_elements(yrange) lt 2) then begin
       minx = min(elm_data[4,*])
       maxx = max(elm_data[4,*])
       miny = min(elm_data[5,*])
       maxy = max(elm_data[5,*])
       
       for i=long(0),nelms-1 do begin
           if(threed eq 1 and i ge (nelms/mesh.nplanes._data)) then break
           i_data = elm_data[*,i]
           a = i_data[0]
           b = i_data[1]
           c = i_data[2]
           t = i_data[3]
           x = i_data[4]
           y = i_data[5]
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

   if(version lt 15) then begin
      ib = 6
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
       edge_val = 0
       num_edge_vals = 0
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

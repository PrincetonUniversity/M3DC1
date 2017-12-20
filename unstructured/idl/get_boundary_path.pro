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

function find_next_boundary_point, list, xy, mesh=mesh, index=index
   tol = 1.e-6
   n = n_elements(list)

   if(n_elements(index) eq 0) then index=long(-1)

   i0 = index

   for j=long(0), n-1 do begin
      index = long(j + i0 + 1) mod n
      i = list[index]

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
          izone = (bound and 120)/2^3 + 1
          if(izone eq 1 or izone eq 2) then begin
             if(n_elements(xy) eq 0) then return, p1
             if(abs(p1[0] - xy[0]) lt tol and $
                abs(p1[1] - xy[1]) lt tol) then begin
                return, p2
             end
          end
       endif
       if((bound and 2) eq 2) then begin
          izone = (bound and 1920)/2^7 + 1
          if(izone eq 1 or izone eq 2) then begin
             if(n_elements(xy) eq 0) then return, p2
             if(abs(p2[0] - xy[0]) lt tol and $
                abs(p2[1] - xy[1]) lt tol) then begin
                return, p3
             end
          end
       endif
       if((bound and 4) eq 4) then begin
          izone = (bound and 30720)/2^11 + 1
          if(izone eq 1 or izone eq 2) then begin
             if(n_elements(xy) eq 0) then return, p3
             if(abs(p3[0] - xy[0]) lt tol and $
                abs(p3[1] - xy[1]) lt tol) then begin
                return, p1
             end
          end
       endif
   end

   print, 'Error: no next boundary point found'
   return, xy
end

function get_boundary_path, mesh=mesh, _EXTRA=extra
  tol = 1e-6

   if(n_elements(mesh) eq 0) then mesh = read_mesh(_EXTRA=ex)
   if(n_tags(mesh) eq 0) then return, [0,0]

   nelms = mesh.nelms._data

   nbound = 0
   for i=long(0), nelms-1 do begin
       bound = fix(mesh.elements._data[6,i])
       if((bound and 1) eq 1) then begin
          izone = (bound and 30720)/2^3 + 1 
          if(izone eq 1 or izone eq 2) then nbound = nbound + 1
       end
       if((bound and 2) eq 2) then begin
          izone = (bound and 30720)/2^7 + 1 
          if(izone eq 1 or izone eq 2) then nbound = nbound + 1
       end
       if((bound and 4) eq 4) then begin
          izone = (bound and 30720)/2^11 + 1 
          if(izone eq 1 or izone eq 2) then nbound = nbound + 1
       end
    end
   print, 'Number of boundary points: ', nbound
   list = lonarr(nbound)
   j = 0
   for i=long(0), nelms-1 do begin
       bound = fix(mesh.elements._data[6,i])
       if((bound and 1) eq 1) then begin
          izone = (bound and 30720)/2^3 + 1 
          if(izone eq 1 or izone eq 2) then begin
             list[j] = i
             j = j+1
          end
       end
       if((bound and 2) eq 2) then begin
          izone = (bound and 30720)/2^7 + 1 
          if(izone eq 1 or izone eq 2) then begin
             list[j] = i
             j = j+1
          end
       end
       if((bound and 4) eq 4) then begin
          izone = (bound and 30720)/2^11 + 1 
          if(izone eq 1 or izone eq 2) then begin
             list[j] = i
             j = j+1
          end
       end
    end

   xy = fltarr(2,nbound)

   xy[*,0] = find_next_boundary_point(list,mesh=mesh,index=index)

   j = 1
   for i=1, nbound-1 do begin
      tmp = find_next_boundary_point(list,xy[*,j-1],mesh=mesh,index=index)

      if(j ge 2) then begin
         if(abs(tmp[0] - xy[0,j-2]) lt tol and $
            abs(tmp[1] - xy[1,j-2]) lt tol) then begin
            continue
         end
      end

      xy[*,j] = tmp
      j = j+1
   end
   nbound = j

   print, 'found ', nbound, ' points'
   
   return, xy[*,0:nbound-1]
end

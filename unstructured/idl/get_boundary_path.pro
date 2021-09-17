function find_next_boundary_point, list, xy, mesh=mesh, index=index, $
                                   imultiregion=imulti
   tol = 1.d-6
   n = n_elements(list)

   if(n_elements(index) eq 0) then index=long(-1)

   i0 = index
   elm_data = mesh.elements._data
   for j=long(0), n-1 do begin
      index = long(j + i0 + 1) mod n
      i = list[index]
      i_data = elm_data[*,i]

       a = i_data[0]
       b = i_data[1]
       c = i_data[2]
       t = i_data[3]
       x = i_data[4]
       y = i_data[5]
       bound = fix(i_data[6])

       p1 = [x, y]
       p2 = p1 + [(b+a) * cos(t), (b+a) * sin(t)]
       p3 = p1 + [b * cos(t) - c * sin(t), $
                  b * sin(t) + c * cos(t)]
       
       if((bound and 1) eq 1) then begin
          izone = (bound and 120)/2^3 + 1
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then begin
             if(n_elements(xy) eq 0) then return, p1
             if(abs(p1[0] - xy[0]) lt tol and $
                abs(p1[1] - xy[1]) lt tol) then begin
                return, p2
             end
          end
       endif
       if((bound and 2) eq 2) then begin
          izone = (bound and 1920)/2^7 + 1
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then begin
             if(n_elements(xy) eq 0) then return, p2
             if(abs(p2[0] - xy[0]) lt tol and $
                abs(p2[1] - xy[1]) lt tol) then begin
                return, p3
             end
          end
       endif
       if((bound and 4) eq 4) then begin
          izone = (bound and 30720)/2^11 + 1
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then begin
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

function get_boundary_path, mesh=mesh, imultiregion=imulti, _EXTRA=ex, $
                            normal=norm, center=center, angle=angle, $
                            length=length
  tol = 1d-6
  if(n_elements(imulti) eq 0) then $
     imulti=read_parameter('imulti_region', _EXTRA=ex)
  

   if(n_elements(mesh) eq 0) then mesh = read_mesh(_EXTRA=ex)
   if(n_tags(mesh) eq 0) then return, [0,0]

   nelms = mesh.nelms._data
   elm_data = mesh.elements._data
   nbound = 0
   for i=long(0), nelms-1 do begin
       bound = fix(elm_data[6,i])
       if((bound and 1) eq 1) then begin
          izone = (bound and 30720)/2^3 + 1 
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then nbound = nbound + 1
       end
       if((bound and 2) eq 2) then begin
          izone = (bound and 30720)/2^7 + 1 
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then nbound = nbound + 1
       end
       if((bound and 4) eq 4) then begin
          izone = (bound and 30720)/2^11 + 1 
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then nbound = nbound + 1
       end
    end
   print, 'Number of boundary points: ', nbound
   list = lonarr(nbound)
   j = 0
   for i=long(0), nelms-1 do begin
       bound = fix(elm_data[6,i])
       if((bound and 1) eq 1) then begin
          izone = (bound and 30720)/2^3 + 1 
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then begin
             list[j] = i
             j = j+1
          end
       end
       if((bound and 2) eq 2) then begin
          izone = (bound and 30720)/2^7 + 1 
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then begin
             list[j] = i
             j = j+1
          end
       end
       if((bound and 4) eq 4) then begin
          izone = (bound and 30720)/2^11 + 1 
          if(imulti eq 1) then begin
             ibound = izone eq 1 or izone eq 2
          endif else begin
             ibound = izone ne 0
          end
          if(ibound) then begin
             list[j] = i
             j = j+1
          end
       end
    end

   xy = dblarr(2,nbound)

   xy[*,0] = find_next_boundary_point(list,mesh=mesh,index=index, $
                                      imultiregion=imulti)

   j = 1
   closed = 0
   for i=1, nbound-1 do begin
      tmp = find_next_boundary_point(list,xy[*,j-1],mesh=mesh,index=index, $
                                    imultiregion=imulti)

      if(j ge 2) then begin
         ; if next point is same as previous point, we're backtracking.
         if(abs(tmp[0] - xy[0,j-2]) lt tol and $
            abs(tmp[1] - xy[1,j-2]) lt tol) then begin
            continue
         end

         ; if next point is same as first point, we're done.
         if(abs(tmp[0] - xy[0,0]) lt tol and $
            abs(tmp[1] - xy[1,0]) lt tol) then begin
            closed = 1
            break
         end
      end
         
      xy[*,j] = tmp
      j = j+1
   end
   nbound = j
   print, 'found ', nbound, ' points'
   xy = xy[*,0:nbound-1]

   if(closed eq 0) then begin
      print, 'WARNING: boundary path is not closed'
   end

   center = [(max(xy[0,*]) + min(xy[0,*]))/2.D, $
             (max(xy[1,*]) + min(xy[1,*]))/2.D]

   angle = reform(atan(xy[1,*] - center[1], xy[0,*] - center[0]))

   ; if angle is decreasing, then reverse order
   if(median(deriv(angle)) lt 0.) then begin
      print, 'reversing angle!'
      angle = reverse(angle)
      xy = reverse(xy, 2)
   end

   ; clamp angles to [0, 2pi)
   i = where(angle lt 0., count)
   if(count gt 0) then angle[i] = angle[i] + 2.D*!dpi

   ; shift values to start at minimum angle
   a0 = min(angle, i)
   angle = shift(angle, -(i+1))
   xy[0,*] = shift(xy[0,*], -(i+1))
   xy[1,*] = shift(xy[1,*], -(i+1))

   ; calculate normals
   norm = dblarr(2,nbound)
   for i=0, nbound-1 do begin
      if(i eq nbound-1) then ip = 0 else ip = i+1
      if(i eq 0) then im = nbound-1 else im = i-1

      nm = [xy[1,i] - xy[1,im], -(xy[0,i] - xy[0,im])]
      np = [xy[1,ip] - xy[1,i], -(xy[0,ip] - xy[0,i])]
      nm = nm / sqrt(nm[0]^2 + nm[1]^2)
      np = np / sqrt(np[0]^2 + np[1]^2)

      norm[*,i] = (nm + np)/2.D
      norm[*,i] = norm[*,i] / sqrt(norm[0,i]^2 + norm[1,i]^2)
   end

   length = dblarr(nbound)
   for i=1, nbound-1 do begin
      length[i] = length[i-1] + $
                  sqrt((xy[0,i]-xy[0,i-1])^2 + (xy[1,i]-xy[1,i-1])^2)
   end
   
   return, xy
end

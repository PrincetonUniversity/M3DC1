function path_contains_point, path, pt
  point_in_plane = 0

 ; For each line segment, determine if line running from
 ; (x,y) to (infinity, y) passes through the segment.
 ; If an odd number of segments are crossed, the point is inside the plane
  
  print, 'Searching for path containing point ', pt
  x = pt[0]
  y = pt[1]
  n = n_elements(path[0,*])

  for i=0, n-1 do begin
     x0 = path[0,i]
     y0 = path[1,i]
     if(i lt n-1) then begin
        x1 = path[0,i+1]
        y1 = path[1,i+1]
     endif else begin
        x1 = path[0,0]
        y1 = path[1,0]
     endelse
     if(y gt y0 and y gt y1) then continue

     if(y lt y0 and y lt y1) then continue

     if(y0 eq y1) then begin
        point_in_plane = not point_in_plane
        continue
     endif

     m = (x1 - x0) / (y1 - y0)
     b = x0 - m*y0
     if(x le m*y + b) then point_in_plane = not point_in_plane
  end
  
  return, point_in_plane
end

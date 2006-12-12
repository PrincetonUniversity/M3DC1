function time_slice, data, t
   names = tag_names(data)
   label = string(FORMAT='("TIME_",I3.3)', t)

   ncmp = strcmp(label, names)
   nmax = max(ncmp, i)
   if(nmax eq 0) then begin
       print, label, " not found."
       return, 0
   endif

   return, data.(i)
end

pro plot_mesh, slice, color=col, linestyle=lin, oplot=oplot
   nelms = slice.mesh.nelms._data
   
   if(not keyword_set(oplot)) then begin
       plot, slice.mesh.elements._data[4,*], $
         slice.mesh.elements._data[5,*], psym = 3
   endif  

   for i=0, nelms-1 do begin
       a = slice.mesh.elements._data[0,i]
       b = slice.mesh.elements._data[1,i]
       c = slice.mesh.elements._data[2,i]
       t = slice.mesh.elements._data[3,i]
       x = slice.mesh.elements._data[4,i]
       y = slice.mesh.elements._data[5,i]

       p1 = [x, y]
       p2 = p1 + [(b+a) * cos(t), (b+a) * sin(t)]
       p3 = p1 + [b * cos(t) - c * sin(t), $
                  b * sin(t) + c * cos(t)]
       
        oplot, [p1[0], p2[0]], [p1[1], p2[1]], color=col, linestyle=lin
        oplot, [p2[0], p3[0]], [p2[1], p3[1]], color=col, linestyle=lin
        oplot, [p3[0], p1[0]], [p3[1], p1[1]], color=col, linestyle=lin
   end
end

function is_in_tri, localp, a, b, c

   small = (a+b+c)*1e-6

   if(localp[1] lt 0 - small) then return, 0
   if(localp[1] gt c + small) then return, 0
   
   x = 1.-localp[1]/c
   if(localp[0] lt -b*x - small) then return, 0
   if(localp[0] gt a*x + small) then return, 0

   return, 1
end


function eval, field, localpos, elm

   mi = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0]
   ni = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5]
   sum = 0.

   for p=0, 19 do begin
       sum = sum + field[p,elm]*(localpos[0]^mi[p]*localpos[1]^ni[p])
   end

   return, sum
end

function eval_field, field, mesh, x=xi, y=yi, points=p

   nelms = mesh.nelms._data 
   print, nelms, " elements."

   if(n_elements(p) eq 0) then p = 100
   
   minpos = fltarr(2)
   maxpos = fltarr(2)
   pos = fltarr(2)
   localpos = fltarr(2)
   index = intarr(2)

   xrange = [0.,mesh.width._data]
   yrange = [0.,mesh.height._data]

   dx = (xrange[1] - xrange[0]) / (p - 1)
   dy = (yrange[1] - yrange[0]) / (p - 1)

   result = fltarr(p,p)

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
       p2 = p1 + [(b+a) * co, (b+a) * sn]
       p3 = p1 + [b * co - c * sn, b * sn + c * co]

       minpos = [min([p1[0], p2[0], p3[0]]), min([p1[1], p2[1], p3[1]])]
       maxpos = [max([p1[0], p2[0], p3[0]]), max([p1[1], p2[1], p3[1]])]
             
       index[1] = minpos[1]/dy
       pos[1] = index[1]*dy

       while(pos[1] le maxpos[1]) do begin
           index[0] = minpos[0]/dx
           pos[0] = index[0]*dx

           while(pos[0] le maxpos[0]) do begin
               localpos = [(pos[0]-x)*co + (pos[1]-y)*sn - b, $
                          -(pos[0]-x)*sn + (pos[1]-y)*co]

               if(is_in_tri(localpos,a,b,c) eq 1) then begin
                   result[index[0], index[1]] = eval(field, localpos, i)
               endif
               
               pos[0] = pos[0] + dx
               index[0] = index[0] + 1
           end
           pos[1] = pos[1] + dy
           index[1] = index[1] + 1
       end
   end

   xi = findgen(p)*(xrange[1]-xrange[0])/(p-1) + xrange[0]
   yi = findgen(p)*(yrange[1]-yrange[0])/(p-1) + yrange[0]

   return, result
end


pro plot_field, name, time, lines=lines, nlevels=nlevels, points=p, $
                range=range
  print, "Reading data..."

  result = h5_parse('C1.h5', /read_data)

  slice = time_slice(result, time)

  names = tag_names(slice.fields)
  ncmp = strcmp(name, names, /fold_case)
  nmax = max(ncmp, i)
  if(nmax eq 0) then begin
      print, "Field ", name, " not found."
      print, names
      return
  endif

  field = slice.fields.(i)._data

  print, "Evaluating field..."

  data = eval_field(field, slice.mesh, p=p, x=x, y=y)

  print, "Plotting..."

  contour_and_legend, data, x, y, nlevels=nlevels, lines=lines, range=range
end




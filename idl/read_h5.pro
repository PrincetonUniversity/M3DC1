function to_title, name
   if(strcmp(name, 'psi') eq 1) then return, "!7w"
   if(strcmp(name, 'phi') eq 1) then return, "!8U"
   if(strcmp(name, 'chi') eq 1) then return, "!7v"
   if(strcmp(name, 'vor') eq 1) then return, "!17z!9.GX!17v"
   if(strcmp(name, 'com') eq 1) then return, "!9G.!17v"
   
   return, "!8" + name
end

function ntime, filename

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   ntime_id = h5a_open_name(root_id, "ntime")
   t = h5a_read(ntime_id)
   h5a_close, ntime_id
   h5g_close, root_id
   h5f_close, file_id

   return, t
end

function read_scalars, filename=filename
   if(n_elements(filename) eq 0) then filename='C1.h5'

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   scalars = h5_parse(root_id, "scalars", /read_data)
   h5g_close, root_id
   h5f_close, file_id

   return, scalars
end

function numvar, filename

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   ntime_id = h5a_open_name(root_id, "numvar")
   nv = h5a_read(ntime_id)
   h5a_close, ntime_id
   h5g_close, root_id
   h5f_close, file_id

   return, nv
end

function time_name, t
   label = string(FORMAT='("time_",I3.3)', t)
   return, label
end


pro plot_mesh, mesh, color=col, linestyle=lin, oplot=oplot
   nelms = mesh.nelms._data
   
   if(not keyword_set(oplot)) then begin
       plot, mesh.elements._data[4,*], $
         mesh.elements._data[5,*], psym = 3
   endif  

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
       
       oplot, [p1[0], p2[0]], [p1[1], p2[1]], color=col, linestyle=lin
       oplot, [p2[0], p3[0]], [p2[1], p3[1]], color=col, linestyle=lin
       oplot, [p3[0], p1[0]], [p3[1], p1[1]], color=col, linestyle=lin
   end
end

function is_in_tri, localp, a, b, c

   small = (a+b+c)*1e-6

   if(localp[1] lt 0. - small) then return, 0
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

   if(n_elements(p) eq 0) then p = 100
   
   minpos = fltarr(2)
   maxpos = fltarr(2)
   pos = fltarr(2)
   localpos = fltarr(2)
   index = intarr(2)

   xrange = [0.,mesh.width._data]
   yrange = [0.,mesh.height._data]

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
             
       index[1] = minpos[1]/dy
       pos[1] = index[1]*dy

       while(pos[1] le maxpos[1] + small*dy) do begin
           index[0] = minpos[0]/dx
           pos[0] = index[0]*dx

           while(pos[0] le maxpos[0] + small*dx) do begin
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

   xi = findgen(p)*(xrange[1]-xrange[0])/(p-1.) + xrange[0]
   yi = findgen(p)*(yrange[1]-yrange[0])/(p-1.) + yrange[0]

   return, result
end

function read_raw_field, name, time, mesh=mesh, filename=filename, time=t
   if(n_elements(filename) eq 0) then filename='C1.h5'

   nt = ntime(filename)

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


function read_field, name, slices=time, mesh=mesh, filename=filename, time=t, $
                     x=x, y=y, points=pts
   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(pts) eq 0) then pts = 50

   nt = ntime(filename)

   if(n_elements(time) eq 0) then begin
       time = [0,nt-1]
   endif else if(n_elements(time) eq 1) then begin
       time = [time, time]
   endif

   if((time[0] ge nt) or (time[1] ge nt)) then begin
       print, "Error: there are only ", nt-1, " time slices."
       return, 0
   endif

   data = fltarr(time[1]-time[0]+1, pts, pts)
   t = fltarr(time[1]-time[0]+1)

   file_id = h5f_open(filename)

   for i=time[0], time[1] do begin

       print, 'Time ', i
       print, '  reading...'

       time_group_id = h5g_open(file_id, time_name(i))
       mesh = h5_parse(time_group_id, 'mesh', /read_data)   

       field_group_id = h5g_open(time_group_id, 'fields')
       field = h5_parse(field_group_id, name, /read_data)

       time_id = h5a_open_name(time_group_id, "time")
       t[i-time[0]] = h5a_read(time_id)
       h5a_close, time_id

       h5g_close, field_group_id
       h5g_close, time_group_id

       print, '  evaluating...'

       data[i-time[0],*,*]=eval_field(field._data, mesh, $
                                      points=pts, x=x, y=y)
   end

   h5f_close, file_id

   return, data
end


pro plot_field, name, time, points=p, filename=filename, mesh=plotmesh, $
                _EXTRA = ex

   field = read_field(name, slices=time, mesh=mesh, filename=filename, $,
                      time=t, x=x, y=y, points=p)

   print, "Plotting..."
   !x.title = '!8x/L!3'
   !y.title = '!8y/L!3'
   if(t gt 0) then begin
       title = "!8" + to_title(name) + $
         string(FORMAT='("!6(!8t!6 = ",G0," !7X!D!8i0!N!6!U-1!N)!3")', t)
   endif else begin
       title = "!8" + to_title(name) + $
         string(FORMAT='("!6(!8t!6 = ",G0,")!3")', t)
   endelse

   contour_and_legend, field, x, y, title=title, _EXTRA=ex

   if(keyword_set(plotmesh)) then begin
       plot_mesh, mesh, color=color(2), /oplot
   endif

   print, "Done."
end

pro compare, file1, file2, time, names=names

   if(n_elements(names) eq 0) then begin
       names = ['psi', 'phi']
   endif

   for i=0, n_elements(names)-1 do begin
       print, "====================="
       print, "Comparing ", names[i]
       print, "---------------------"

       vals1 = read_field(names[i], time, filename=file1)
       vals2 = read_field(names[i], time, filename=file2)

       rms = sqrt((vals1 - vals2)^2)
       d = mean(rms)
       v = variance(rms)

       print, "Mean/variance of RMS error in ", names[i], " = ", $
         d, " / ",  v
   end
end

pro plot_energy, filename=filename, diff=diff

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = numvar(filename)

   scalars = read_scalars(filename=filename)

   !x.title = '!8t!3'

   loadct, 12
   dc = !d.table_size / 7

   N = n_elements(scalars.time._data)

   E_K = scalars.E_KP._data 
   E_M = scalars.E_MP._data
   E_D = scalars.E_KPD._data + scalars.E_MPD._data
   E_H = scalars.E_KPH._data + scalars.E_MPH._data
   if(nv ge 2) then begin
       E_K = E_K + scalars.E_KT._data 
       E_M = E_M + scalars.E_MT._data
       E_D = E_D + scalars.E_KTD._data + scalars.E_MTD._data
       E_H = E_H + scalars.E_KTH._data + scalars.E_MTH._data
   endif
   if(nv ge 3) then begin
       E_K = E_K + scalars.E_K3._data
       E_M = E_M + scalars.E_P._data
       E_D = E_D + scalars.E_K3D._data + scalars.E_PD._data
       E_H = E_H + scalars.E_K3H._data + scalars.E_PH._data
   endif

   E = E_K + E_M

;   E_D[1:N-1] = (E_D[1:N-1] + E_D[0:N-2])/2.
;   E_H[1:N-1] = (E_H[1:N-1] + E_H[0:N-2])/2.

   if(nv ge 3) then begin
       Error = deriv(scalars.time._data, E) - E_D - E_H   
   endif else begin
       Error = deriv(scalars.time._data, E) - E_H
   endelse


   if(keyword_set(diff)) then begin
       E = deriv(scalars.time._data,E)
       E_K = deriv(scalars.time._data, E_K)
       E_M = deriv(scalars.time._data, E_M)
       E_D = deriv(scalars.time._data, E_D)
       ylog = 0
       !y.title = '!6d!8E!6/d!8t!3'
   endif else begin
       ylog = 1
       !y.title = '!6 Energy!3'
   endelse

   !x.range = 0
   if(!y.range[0] eq 0 and !y.range[1] eq 0) then begin
       !y.range = [min([E, E_K, E_M, -E_D, -E_H], /nan)/2., $
                   max([E, E_K,  E_M, -E_D, -E_H], /nan)*2.]
       if((ylog eq 1) and (!y.range[0] le 0.)) then !y.range[0] = 1e-5
   endif

   plot, scalars.time._data, E, ylog=ylog
   oplot, scalars.time._data, E_K, color=dc, linestyle = 2
   oplot, scalars.time._data, E_M, color=2*dc, linestyle = 2
   oplot, scalars.time._data, -E_D, color=3*dc, linestyle = 1
   oplot, scalars.time._data, -E_H, color=4*dc, linestyle = 1
   oplot, scalars.time._data, abs(Error), color=5*dc, linestyle = 1

   plot_legend, ['Total', 'Kinetic', 'Magnetic', $
                 'Diffusive', 'Hyper-Diffusive', '|Error|'] , $
     color=[-1,1,2,3,4,5,6,7,8,9]*dc, ylog=ylog

   print, "Error = ", Error

end


pro plot_field_mpeg, fieldname, mpegame=mpegname, range=range, points=pts, $
                     _EXTRA=ex
    if(n_elements(mpegname) eq 0) then mpegname = fieldname + '.mpeg'
    if(n_elements(pts) eq 0) then pts = 50

    nt = ntime('C1.h5')

    data = read_field(fieldname, mesh=mesh, x=x, y=y, points=pts)

    contour_and_legend_mpeg, mpegname, data, x, y, _EXTRA=ex
end

function hdf5_file_test, filename
   if(file_test(filename) eq 0) then begin
       print, "Error: ", filename, " is not a valid file."
       return, 0
   endif else if(h5f_is_hdf5(filename) eq 0) then begin
       print, "Error: ", filename, " is not a valid HDF5 file."
       return, 0
   endif else return, 1
end

function translate, name
   if(strcmp(name, 'psi') eq 1) then return, "!7w"
   if(strcmp(name, 'phi') eq 1) then return, "!8U"
   if(strcmp(name, 'chi') eq 1) then return, "!7v"
   if(strcmp(name, 'vor') eq 1) then return, "!17z!9.GX!17v"
   if(strcmp(name, 'com') eq 1) then return, "!9G.!17v"
   if(strcmp(name, 'eta') eq 1) then return, "!7g!3"
   if(strcmp(name, 'den') eq 1) then return, "!8n!3"
   
   return, "!8" + name
end

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

function read_parameter, name, filename=filename, print=pr
   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(hdf5_file_test(filename) eq 0) then return, 0

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   attr = read_attribute(root_id, name)
   h5g_close, root_id
   h5f_close, file_id

   if(keyword_set(pr)) then print, name, " = ", attr

   return, attr
end

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
   label = string(FORMAT='("time_",I3.3)', t)
   return, label
end


pro plot_mesh, mesh, color=col, linestyle=lin, oplot=oplot, $
               filename=filename, _EXTRA=ex
   nelms = mesh.nelms._data
   
   if(not keyword_set(oplot)) then begin
       plot, mesh.elements._data[4,*], $
         mesh.elements._data[5,*], psym = 3, _EXTRA=ex
   endif  

   if(n_elements(col) ne 0) then loadct, 12
  
   xzero = read_parameter("xzero", filename=filename)
   zzero = read_parameter("zzero", filename=filename)

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
       
       oplot, [p1[0],p2[0]]+xzero, [p1[1],p2[1]]+zzero, $
         color=col, linestyle=lin, thick=.2
       oplot, [p2[0],p3[0]]+xzero, [p2[1],p3[1]]+zzero, $
         color=col, linestyle=lin, thick=.2
       oplot, [p3[0],p1[0]]+xzero, [p3[1],p1[1]]+zzero, $
         color=col, linestyle=lin, thick=.2
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


function eval, field, localpos, elm, operation=op

   mi = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0]
   ni = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5]
   sum = 0.

   if(n_elements(op) eq 0) then op = 1

   for p=0, 19 do begin
       case op of
           1: sum = sum + field[p,elm]*(localpos[0]^mi[p]*localpos[1]^ni[p])
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


function eval_field, field, mesh, r=xi, z=yi, points=p, operation=op, $
                     filename=filename, xrange=xrange, yrange=yrange

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nelms = mesh.nelms._data 

   if(n_elements(p) eq 0) then p = 100
   if(n_elements(op) eq 0) then op = 1
   
   minpos = fltarr(2)
   maxpos = fltarr(2)
   pos = fltarr(2)
   localpos = fltarr(2)
   index = intarr(2)

   xzero = read_parameter("xzero", filename=filename)
   zzero = read_parameter("zzero", filename=filename)

   if(n_elements(xrange) lt 2) then $
     xrange = [0.,mesh.width._data] + xzero
   if(n_elements(yrange) lt 2) then $
     yrange = [0.,mesh.height._data] + zzero

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
                             eval(field, localpos, i, op=op)
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


function read_field, name, slices=time, mesh=mesh, filename=filename, time=t, $
                     r=x, z=y, points=pts, xrange=xrange, yrange=yrange
   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(pts) eq 0) then pts = 50

   if(hdf5_file_test(filename) eq 0) then return, 0

   nt = read_parameter("ntime", filename=filename)

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

   if(1 eq strcmp('j', name)) then begin
       newname = 'psi'
       op = 7
   endif else if (1 eq strcmp('vor', name)) then begin
       newname = 'phi'
       op = 7
   endif else begin
       newname = name
       op = 1
   endelse   

   for i=time[0], time[1] do begin

       print, 'Time ', i
       print, '  reading...'

       time_group_id = h5g_open(file_id, time_name(i))
       mesh = h5_parse(time_group_id, 'mesh', /read_data)   

       field_group_id = h5g_open(time_group_id, 'fields')
       field = h5_parse(field_group_id, newname, /read_data)

       time_id = h5a_open_name(time_group_id, "time")
       t[i-time[0]] = h5a_read(time_id)
       h5a_close, time_id

       h5g_close, field_group_id
       h5g_close, time_group_id

       print, '  evaluating...'

       data[i-time[0],*,*]=eval_field(field._data, mesh, points=pts, $
                                      r=x, z=y, op=op, filename=filename, $
                                      xrange=xrange, yrange=yrange)
   end

   h5f_close, file_id

   if(1 eq strcmp('j', name)) then data = -data

   return, data
end


pro plot_field, name, time, points=p, filename=filename, mesh=plotmesh, $
                mcolor=mc, lcfs=lcfs, title=title, $
                maskrange=maskrange, maskfield=maskfield, $
                xrange=xrange, yrange=yrange, _EXTRA = ex

   field = read_field(name, slices=time, mesh=mesh, filename=filename, $,
                      time=t, r=x, z=y, points=p, xrange=xrange, yrange=yrange)
   if(n_elements(field) le 1) then return

   if(n_elements(maskrange) eq 2) then begin
       if(strcmp(name, maskfield) eq 1) then begin
           psi = field
       endif else begin
           psi = read_field(maskfield, slices=time, mesh=mesh, $
                            filename=filename, points=p, $
                            xrange=xrange, yrange=yrange)
       endelse
       mask = (psi ge maskrange[0]) and (psi le maskrange[1])
       field = mask*field
   endif

   print, "Plotting..."
   !x.title = '!8r!3'
   !y.title = '!8z!3'
   if(n_elements(title) eq 0) then begin
       if(t gt 0) then begin
           title = "!8" + translate(name) + $
             string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!N!6)!3")', t)
       endif else begin
           title = "!8" + translate(name) + $
             string(FORMAT='("!6(!8t!6 = ",G0,")!3")', t)
       endelse
   endif

   contour_and_legend, field, x, y, title=title, _EXTRA=ex

   if(keyword_set(lcfs) or n_elements(maskrange) ne 0) then begin
       if(n_elements(psi) eq 0) then begin
           plot_lcfs, time, color=130, points=p, filename=filename
       endif else begin
           plot_lcfs, time, color=130, val=maskrange[0], psi=psi, x=x, y=y, $
             points=p
           plot_lcfs, time, color=130, val=maskrange[1], psi=psi, x=x, y=y, $
             points=p
       endelse
   endif

   if(keyword_set(plotmesh)) then begin
       plot_mesh, mesh, color=220, /oplot, filename=filename
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

function energy, filename=filename, error=error
   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = read_parameter("numvar", filename=filename)
   s = read_scalars(filename=filename)

   N = n_elements(s.time._data)

   E_K = s.E_KP._data
   E_M = s.E_MP._data
   E_D = s.E_KPD._data + s.E_MPD._data
   E_H = s.E_KPH._data + s.E_MPH._data
   if(nv ge 2) then begin
       E_K = E_K + s.E_KT._data 
       E_M = E_M + s.E_MT._data
       E_D = E_D + s.E_KTD._data + s.E_MTD._data
       E_H = E_H + s.E_KTH._data + s.E_MTH._data
   endif
   if(nv ge 3) then begin
       E_K = E_K + s.E_K3._data
       E_M = E_M + s.E_P._data
       E_D = E_D + s.E_K3D._data + s.E_PD._data
       E_H = E_H + s.E_K3H._data + s.E_PH._data
   endif

   E = E_K + E_M

;  Account for energy terms not included in the physical model
   if(nv le 2) then begin
       dissipated = E_D + E_H
   endif else begin
       dissipated = 0.
   endelse

   ; Account for fluxes across boundary
   dissipated = dissipated    $
     + s.Flux_diffusive._data $
     + s.Flux_pressure._data  $
     + s.Flux_kinetic._data   $
     + 0.*s.Flux_poynting._data  $
     + s.Flux_thermal._data

   eloop = s.loop_voltage._data * s.toroidal_current._data / (2.*3.14159265)
   dissipated = dissipated - eloop

   Error = E - E[0]
   total_lost = fltarr(n_elements(Error))
   total_lost[0] = 0.
   for i=1, n_elements(Error)-1 do begin
       dt = s.time._data[i]-s.time._data[i-1]
       total_lost[i] = total_lost[i-1] + $
         dt*(dissipated[i-1] + dissipated[i])/2.
   endfor

   Error = Error - total_lost

   return, E
end


pro plot_fluxes, filename=filename, ylog=ylog, overplot=overplot, _EXTRA=extra

   s = read_scalars(filename=filename)

   yrange=[min([s.Flux_diffusive._data,s.Flux_pressure._data, $
                s.Flux_kinetic._data,  s.Flux_poynting._data, $
                s.Flux_thermal._data]), $
           max([s.Flux_diffusive._data,s.Flux_pressure._data, $
                s.Flux_kinetic._data,  s.Flux_poynting._data, $
                s.Flux_thermal._data])]*1.2

   xtitle = '!8t !6(!7s!D!8A!N!6)!3'
   ytitle = '!6Energy Flux (!8B!D0!U2!N/4!7p s!D!8A!N L!U2!N)!3'

   if(keyword_set(overplot)) then begin
       oplot,  s.time._data, s.Flux_diffusive._data, color=color(0,5), $
         _EXTRA=extra
   endif else begin
       plot,  s.time._data, s.Flux_diffusive._data, color=color(0,5), $
         ylog=ylog, xlog=xlog, yrange=yrange, xtitle=xtitle, ytitle=ytitle, $
         _EXTRA=extra
   endelse
   oplot, s.time._data, s.Flux_kinetic._data, color=color(1,5)
   oplot, s.time._data, s.Flux_pressure._data, color=color(2,5)
   oplot, s.time._data, s.Flux_poynting._data, color=color(3,5)
   oplot, s.time._data, s.Flux_thermal._data, color=color(4,5)

   plot_legend, ['Diffusive', 'KE Conv', 'Pressure Conv.', $
                 'Poynting', 'Thermal'], color=colors(5), ylog=ylog, xlog=xlog
end


pro plot_energy, filename=filename, diff=diff, norm=norm, ylog=ylog, $
                 _EXTRA=extra

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = read_parameter("numvar", filename=filename)

   s = read_scalars(filename=filename)

   !x.title = '!8t !6(!7s!D!8A!N!6)!3'

   N = n_elements(s.time._data)

   E_K = s.E_KP._data
   E_M = s.E_MP._data
   E_D = s.E_KPD._data + s.E_MPD._data
   E_H = s.E_KPH._data + s.E_MPH._data
   if(nv ge 2) then begin
       E_K = E_K + s.E_KT._data 
       E_M = E_M + s.E_MT._data
       E_D = E_D + s.E_KTD._data + s.E_MTD._data
       E_H = E_H + s.E_KTH._data + s.E_MTH._data
   endif
   if(nv ge 3) then begin
       E_K = E_K + s.E_K3._data
       E_M = E_M + s.E_P._data
       E_D = E_D + s.E_K3D._data + s.E_PD._data
       E_H = E_H + s.E_K3H._data + s.E_PH._data
   endif

   E = E_K + E_M

;  Account for energy terms not included in the physical model
   if(nv le 2) then begin
       dissipated = E_D + E_H
   endif else begin
       dissipated = 0.
   endelse

   ; Account for fluxes across boundary
   dissipated = dissipated    $
     + s.Flux_diffusive._data $
     + s.Flux_pressure._data  $
     + s.Flux_kinetic._data   $
     - 0.*s.Flux_poynting._data  $
     + s.Flux_thermal._data

   ; account for loop voltage
   eloop = s.loop_voltage._data * s.toroidal_current._data / (2.*3.14159265)
   dissipated = dissipated - eloop

   Error = E - E[0]
   total_lost = fltarr(n_elements(Error))
   total_lost[0] = 0.
   for i=1, n_elements(Error)-1 do begin
       dt = s.time._data[i]-s.time._data[i-1]
       total_lost[i] = total_lost[i-1]  $
         + dt*(dissipated[i-1] + dissipated[i])/2.
   endfor

   Error = Error - total_lost

   print, Error / E
   print, deriv(s.time._data, Error)

   if(keyword_set(norm)) then begin
       E = E - E[0]
       E_K = E_K - E_K[0]
       E_M = E_M - E_M[0]
       E_D = E_D - E_D[0]
       E_H = E_H - E_H[0]
       s.E_KP._data = s.E_KP._data - s.E_KP._data[0]
       s.E_KT._data = s.E_KT._data - s.E_KT._data[0]
       s.E_K3._data = s.E_K3._data - s.E_K3._data[0]
       s.E_MP._data = s.E_MP._data - s.E_MP._data[0]
       s.E_MT._data = s.E_MT._data - s.E_MT._data[0]
       s.E_P._data  = s.E_P._data  - s.E_P._data[0]
   endif 

   if(keyword_set(diff)) then begin
       E = deriv(s.time._data,E)
       E_K = deriv(s.time._data, E_K)
       E_M = deriv(s.time._data, E_M)
       !y.title = '!6d!8E!6/d!8t!3'
   endif else begin
       !y.title = '!6 Energy!3'
   endelse

   if(nv le 2) then begin
       !y.range=[min([E,E_K,E_M,-E_D,-E_H,Error]), $
                 max([E,E_K,E_M,-E_D,-E_H,Error])]
   endif else begin
       !y.range=[min([E,s.E_KP._data,s.E_KT._data,s.E_K3._data, $
                      s.E_MP._data,s.E_MT._data,s.E_P._data,Error]), $
                 max([E,s.E_KP._data,s.E_KT._data,s.E_K3._data, $
                      s.E_MP._data,s.E_MT._data,s.E_P._data,Error])]       
   endelse
   if(keyword_set(ylog) and (!y.range[0] lt 0)) then $
     !y.range[0] = !y.range[1]*1e-8
   if(keyword_set(ylog)) then begin
       !y.range[1] = !y.range[1]*10.
       Error = abs(Error)
   endif
   

   if(nv le 2) then begin
       plot, s.time._data, E, ylog=ylog, color=color(0,6), _EXTRA=extra
       oplot, s.time._data, E_K, color=color(1,6), linestyle = 2
       oplot, s.time._data, E_M, color=color(2,6), linestyle = 2
       oplot, s.time._data, -E_D, color=color(3,6), linestyle = 1
       oplot, s.time._data, -E_H, color=color(4,6), linestyle = 1
       oplot, s.time._data, Error, color=color(5,6), linestyle = 3

       plot_legend, ['Total', 'Kinetic', 'Magnetic', $
                     'Diffusive', 'Hyper-Diffusive', '|Error|'] , $
         color=colors(6), ylog=ylog, $
         linestyle = [0,2,2,1,1,3]
   endif else begin
       plot, s.time._data, E, ylog=ylog, color=color(0,8), _EXTRA=extra
       oplot, s.time._data, s.E_KP._data, color=color(1,8), linestyle = 2
       oplot, s.time._data, s.E_KT._data, color=color(2,8), linestyle = 2
       oplot, s.time._data, s.E_K3._data, color=color(3,8), linestyle = 2
       oplot, s.time._data, s.E_MP._data, color=color(4,8), linestyle = 1
       oplot, s.time._data, s.E_MT._data, color=color(5,8), linestyle = 1
       oplot, s.time._data, s.E_P._data,  color=color(6,8), linestyle = 1
       oplot, s.time._data, Error, color=color(7,8), linestyle = 3

       plot_legend, ['Total', $
                     'KE: Solenoidal', 'KE: Toroidal', 'KE: Compressional', $
                     'ME: Poloidal', 'ME: Toroidal', 'Pressure', $
                     '|Error|'] , $
         color=colors(8), ylog=ylog, $
         linestyle = [0,2,2,2,1,1,1,3]
   endelse
end


pro plot_field_mpeg, fieldname, mpegame=mpegname, range=range, points=pts, $
                     _EXTRA=ex
    if(n_elements(mpegname) eq 0) then mpegname = fieldname + '.mpeg'
    if(n_elements(pts) eq 0) then pts = 50

    nt = get_parameters("ntime")

    data = read_field(fieldname, mesh=mesh, r=x, z=y, points=pts)

    contour_and_legend_mpeg, mpegname, data, x, y, _EXTRA=ex
end


pro plot_lcfs, time, color=color, val=psival, psi=psi, x=x, y=y, points=pts, $
               filename=filename

    if(n_elements(psi) eq 0) then begin
        psi = read_field('psi', slice=time, points=pts, r=x, z=y, $
                         filename=filename)
    endif

    ; if psival not passed, choose limiter value
    if(n_elements(psival) eq 0) then begin
        xlim = lcfs(time, psi=psi, r=x, z=y, points=pts)

;        xlim = read_parameter("xlim")
        print, "xlim = ", xlim

        ; find index corresponding to limiter
        d = min(x-xlim, ilim, /absolute)
    
        ; find value of flux at the limiter
        psival = max(psi[0,ilim,*])
        print, "psi(xlim) = ", psival
    endif

    ; plot contour
    if(n_elements(color) ne 0) then loadct, 12
    contour, psi, x, y, /overplot, nlevels=1, levels=psival, $
      color=color, thick=2
end


pro plot_timings, filename=filename

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

   plot, timings.t_onestep._data, title='!6Timings!3', $
     xtitle='!6Time Step!3', ytitle='!8t!6 (s)!3'
   oplot, timings.t_ludefall._data, linestyle=2, color=30
   oplot, timings.t_sources._data, linestyle=1, color=60
   oplot, timings.t_aux._data, linestyle=1, color=80
   oplot, timings.t_smoother._data, linestyle=1, color=100
   oplot, t_solve, linestyle=2, color=160
   oplot, t_output, linestyle=2, color=200


   plot_legend, ['Onestep', 'ludefall', 'sources', 'aux', $
                 'smoother', 'solve', 'output'], $
     linestyle=[0,2,1,1,1,2,2], color=[-1,30,60,80,100,160,200]

end


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

function beta_normal, filename=filename
   nv = read_parameter("numvar", filename=filename)
   if(nv lt 3) then begin
       print, "Must be numvar = 3 for beta calculation"
       return, 0
   endif 

   gamma = read_parameter('gam', filename=filename)
   bzero = read_parameter('bzero', filename=filename)
   xmag = read_parameter('xmag', filename=filename)
   xlim = read_parameter('xlim', filename=filename)
   psi = read_field('psi',slice=0, mesh=mesh)
   s = read_scalars(filename=filename)

   area = mesh.width._data*mesh.height._data

   if(xmag eq 0.) then begin
       B_T = sqrt(2.*s.E_MT._data/area)
       B_T = B_T[0]
       a = 1.
   endif else begin
       B_T = bzero/xmag
       a = abs(xmag - xlim)
   endelse

   beta_t = (gamma-1.)*s.E_P._data/(s.E_MT._data)
   
   print, 'B_T =', B_T
   print, 'a = ', a
   print, 'beta_T = ', beta_t
   print, 'area = ', area

   return, beta_t / (s.toroidal_current._data/(a*B_T))
end


pro plot_scalar, scalarname, x, filename=filename, names=names, $
                 _EXTRA=extra, overplot=overplot, $
                 ylog=ylog, xlog=xlog, left=left

  if(n_elements(filename) eq 0) then filename='C1.h5'

  !x.title = '!8t!6 (!7s!D!8A!N!6)!3'

  nfiles = n_elements(filename)
  if(nfiles gt 1) then begin
      if(n_elements(c) eq 0) then begin
          colors = colors(nfiles)
      endif else begin
          colors = fltarr(n_elements(scalarname))
          colors[*] = c
      endelse

      for i=0, nfiles-1 do begin
          if(n_elements(x) eq 0) then begin
              plot_scalar, scalarname, filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog
          endif else begin
              plot_scalar, scalarname, x[i], filename=filename[i], $
                overplot=((i gt 0) or keyword_set(overplot)), $
                color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog
          endelse
      end
      if(n_elements(names) gt 0) then begin
          plot_legend, names, color=colors, ylog=ylog, xlog=xlog, left=left
      endif    

      return
  endif 

  s = read_scalars(filename=filename)

  if(strcmp("toroidal current", scalarname, /fold_case) eq 1) or $
    (strcmp("it", scalarname, /fold_case) eq 1) then begin
      data = s.toroidal_current._data
      title = '!6Toroidal Current!3'
      ytitle = '!8I!DT!N!6 (!8cB!D0!N/4!7p!8L!6)!3'
  endif else $
    if (strcmp("toroidal flux", scalarname, /fold_case) eq 1) then begin
      data = s.toroidal_flux._data
      title = '!6Toroidal Flux!3'
      ytitle = '!6Flux (!8L!U2!N B!D0!N)!3'
  endif else $
    if (strcmp("reconnected flux", scalarname, /fold_case) eq 1) then begin
      data = s.reconnected_flux._data
      title = '!6Reconnected Flux!3'
  endif else $
    if (strcmp("loop voltage", scalarname, /fold_case) eq 1) or $
    (strcmp("vl", scalarname, /fold_case) eq 1) then begin
      data = s.loop_voltage._data
      title = '!6Loop Voltage!3'
      ytitle = '!8V!DL!N!3'
  endif else $
    if (strcmp("beta", scalarname, /fold_case) eq 1) then begin
      data = beta(filename=filename)
      title = '!7b!3'
      ytitle = '!7b!3'
  endif else if $
    (strcmp("poloidal beta", scalarname, /fold_case) eq 1) or $
    (strcmp("bp", scalarname, /fold_case) eq 1) then begin
      data = beta_poloidal(filename=filename)
      title = '!7b!D!8p!N!3'
      ytitle = '!7b!D!8p!N!3'
  endif else $
    if (strcmp("normal beta", scalarname, /fold_case) eq 1) or $
    (strcmp("bn", scalarname, /fold_case) eq 1) then begin
      data = beta_normal(filename=filename)
      title = '!7b!D!8N!N!3'
      ytitle = '!7b!D!8N!N!3'
  endif else if $
    (strcmp("kinetic energy", scalarname, /fold_case) eq 1) or $
    (strcmp("ke", scalarname, /fold_case) eq 1)then begin
      nv = read_parameter("numvar", filename=filename)
       data = s.E_KP._data 
       if(nv ge 2) then data = data + s.E_KT._data 
       if(nv ge 3) then data = data + s.E_K3._data
       title = '!6Kinetic Energy!3'
       ytitle = '!6KE (!8B!D0!N!U!62!N/4!7p!6)!3'
  endif else if $
    (strcmp("energy", scalarname, /fold_case) eq 1) then begin
      nv = read_parameter("numvar", filename=filename)
       data = energy(filename=filename)
       title = '!6Energy!3'
       ytitle = '!6E (!8B!D0!N!U!62!N/4!7p!6)!3'
  endif else if $
    (strcmp("error", scalarname, /fold_case) eq 1) then begin
      nv = read_parameter("numvar", filename=filename)
      dummy = energy(filename=filename, error=data)
      data = abs(data)
      title = '!6Error!3'
      ytitle = '!7D!6E (!8B!D0!N!U!62!N/4!7p!6)!3'
  endif else if $
    (strcmp("particles", scalarname, /fold_case) eq 1) or $
    (strcmp("n", scalarname, /fold_case) eq 1) then begin
      data = s.electron_number._data
      title = '!6Particle number!3'
      ytitle = '!8N !6(!8n!D0!N / L!U3!N!6)!3'
  endif else if $
    (strcmp("angular momentum", scalarname, /fold_case) eq 1) then begin
      data = s.angular_momentum._data
      title = '!6Angular Momentum!3'
      ytitle = '!8N !6(!8n!D0!N / L!U3!N!6)!3'
  endif else begin
      print, 'Scalar ', scalarname, ' not recognized.'
      return
  endelse
  
  if(n_elements(x) eq 0) then begin
      if(keyword_set(overplot)) then begin
          oplot, s.time._data, data, color=c, _EXTRA=extra
      endif else begin
          plot, s.time._data, data, $
            title=title, ytitle=ytitle, _EXTRA=extra, ylog=ylog, xlog=xlog, $
            color=c
      endelse
  endif else begin
      xi = x
      x = fltarr(1)
      z = fltarr(1)
      x[0] = xi
      z[0] = data[n_elements(data)-1]

      if(keyword_set(overplot)) then begin
          oplot, x, z, color=c, _EXTRA=extra
      endif else begin
          plot, x, z, $
            title=title, ytitle=ytitle, _EXTRA=extra, ylog=ylog, xlog=xlog, $
          color=c
      endelse
  endelse
end


pro plot_pol_velocity, time, filename=filename, points=pts, maxval=maxval, $
                       _EXTRA=extra

  nv = read_parameter('numvar', filename=filename)

  phi = read_field('phi', filename=filename, $
                   slice=time, r=x, z=z, t=t, points=pts)
  if(n_elements(phi) le 1) then return

  r = radius_matrix(x,z,t)

  vx = -dz(phi,z)/r
  vz =  dx(phi,x)/r

  if(nv ge 3) then begin
      chi = read_field('chi', filename=filename, slice=time, $
                       r=x, z=z, t=t, points=pts)
      vx = vx + dx(chi,x)
      vz = vz + dz(chi,z)
  endif

  bigvel = max(vx^2 + vz^2)
  print, "maximum velocity: ", bigvel
  if(n_elements(maxval) ne 0) then begin
      length = bigvel/maxval
  endif else length=1

  maxstr=string(format='("!6max(!8v!Dpol!N!6) = ",G0)',bigvel)

  velovect, reform(vx), reform(vz), x, z, length=length, _EXTRA=extra, $
    xtitle='!8r !6(!8L!6)!3', ytitle='!8z !6(!8L!6)!3', subtitle=maxstr
end


pro plot_tor_velocity, time, filename=filename, points=pts, _EXTRA=extra

  nv = read_parameter('numvar', filename=filename)

  if(nv lt 2) then begin
      print, "numvar < 2"
      return
  endif

  v = read_field('V', filename=filename, $
                 slice=time, r=x, z=z, t=t, points=pts)
  r = radius_matrix(x,z,t)

  vz = v/r

  contour_and_legend, vz, x, z, _EXTRA=extra
end

  
pro plot_beta, time, filename=filename, points=pts, _EXTRA=extra

  nv    = read_parameter('numvar', filename=filename)
  idens = read_parameter('idens', filename=filename)

  psi = read_field('psi', filename=filename, $
                   slice=time, r=x, z=z, t=t, points=pts)

  if(n_elements(psi) le 1) then return

  r = radius_matrix(x,z,t)

  if(nv ge 2) then begin
      I = read_field('I', filename=filename, slice=time, $
                     r=x, z=z, t=t, points=pts)
  endif else begin
      I = read_parameter('bzero')
  endelse
      
  if(nv ge 3) then begin
      P = read_field('P', filename=filename, slice=time, $
                     r=x, z=z, t=t, points=pts)
  endif else begin
      P = 1
  endelse

  beta = 2.*P/(s_bracket(psi,psi,x,z) + i^2)

  contour_and_legend, beta[0,*,*], x, z, title='!6Local !7b!3', $
    xtitle='!8r!3', ytitle='!8z!3', _EXTRA=extra

  print, "Mean = ", mean(beta)

end


pro plot_q, time, filename=filename, points=pts, _EXTRA=extra
  nv    = read_parameter('numvar', filename=filename)
  idens = read_parameter('idens', filename=filename)

  psi = read_field('psi', filename=filename, $
                   slice=time, r=x, z=z, t=t, points=pts)

  if(n_elements(psi) le 1) then return

  r = radius_matrix(x,z,t)

  if(nv ge 2) then begin
      I = read_field('I', filename=filename, slice=time, $
                     r=x, z=z, t=t, points=pts)
  endif else begin
      I = read_parameter('bzero')
  endelse     

  q = sqrt(I^2/s_bracket(psi,psi,x,z)) < 10

  contour_and_legend, q[0,*,*], x, z, title='!6Local !8q!3', $
    xtitle='!8r!3', ytitle='!8z!3', _EXTRA=extra

end


function flux_average, psi, field, points=pts, flux=flux

   sz = size(psi)

   if(n_elements(pts) eq 0) then pts = min(sz[2:3])

   result = fltarr(sz[1], pts)
   flux = fltarr(sz[1], pts)

   for k=0, sz[1]-1 do begin
       psimin = min(psi[k,*,*])
       psimax = max(psi[k,*,*])
       dpsi = (psimax - psimin)/(pts - 1.)

       for p=0, pts-1 do begin
           fval = dpsi*p  + psimin

           mask = ((psi ge fval) and (psi lt fval+dpsi))
           sfield = mask*field
           
           flux[k,p] = fval
           result[k, p] = total(sfield)/total(mask)
       endfor
   endfor

   return, result
end


; ==============================================================
; nulls
; -----
;
;  Finds field nulls.
;  xpoint = fltarr(2,xpoints): indices of x-point locations
;  axis   = fltarr(2,axes)   : indices of axis locations
; ==============================================================
pro nulls, time, axis=axis, xpoints=xpoint, $
              _EXTRA=extra, psi=psi, r=x, z=z

   if(n_elements(time) eq 0) then time = 0

   if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) $
     then begin
       psi = read_field('psi', r=x, z=z, slice=time, _EXTRA=extra)
   endif
   
   field = s_bracket(psi,psi,x,z)
   d2 = dz(dz(psi,z),z)*dx(dx(psi,x),x)
   
   nulls = field lt mean(field)/1e4

   sz = size(field)

   xpoints = 0
   xpoint = 0.
   axes = 0
   axis = 0.
   
   for i=0, sz[2]-1 do begin
       for j=0, sz[3]-1 do begin
           if(nulls[0,i,j] eq 0) then continue

           currentmin = field[0,i,j]
           currentpos = [i,j]

           ; find local minimum
           for m=i, sz[2]-1 do begin
               for n=j, sz[3]-1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] lt currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                   endif
                   nulls[0,m,n] = 0
               endfor
               for n=j-1, 0, -1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] lt currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                   endif
                   nulls[0,m,n] = 0
               endfor
           endfor

           ; throw out local minima on boundaries
           if (currentpos[0] eq 0) or (currentpos[0] eq sz[2]-1) then continue
           if (currentpos[1] eq 0) or (currentpos[1] eq sz[3]-1) then continue

           ; determine if point is an x-point or an axis and
           ; append the location index to the appropriate array
           if(d2[0,currentpos[0],currentpos[1]] lt 0) then begin
               print, "X-point found at", x[currentpos[0]], z[currentpos[1]]
               xpoints = xpoints + 1
               if(xpoints eq 1) then begin
                   xpoint = currentpos
               endif else begin
                   oldxpoint = xpoint
                   xpoint = intarr(2, xpoints)
                   xpoint[*,0:xpoints-2] = oldxpoint
                   xpoint[*,xpoints-1] = currentpos
               endelse
           endif else begin
               print, "Axis found at", x[currentpos[0]], z[currentpos[1]]
               axes = axes + 1
               if(axes eq 1) then begin
                   axis = currentpos
               endif else begin
                   oldaxis = axis
                   axis = intarr(2, axes)
                   axis[*,0:axes-2] = oldaxis
                   axis[*,axes-1] = currentpos
               endelse
           endelse
       endfor
   endfor 

end

; ========================================================
; lcfs
; ~~~~
;
; returns the flux value of the last closed flux surface
; ========================================================
function lcfs, time, psi=psi, r=x, z=z, _EXTRA=extra
   if(n_elements(time) eq 0) then time = 0

   if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) $
     then begin
       psi = read_field('psi', r=x, z=z, slice=time, _EXTRA=extra)
   endif

   nulls, time, psi=psi, xpoint=xpoint, axis=axis, r=x, z=z, _EXTRA=extra

   ; flux at magnetic axis
   if(n_elements(axis) ne 2) then begin
       print, "Error: there is not exactly one magnetic axis"
       return, 0
   endif
   psi0 = psi[0,axis[0,0],axis[1,0]]
   print, "Flux on axis:", psi0

   ; limiting value
   ; Find limiting flux by calculating outward normal derivative of
   ; the normalized flux.  If this derivative is negative, there is a
   ; limiter.
   sz = size(psi)

   psiz = dz(psi,z)
   psix = dx(psi,x)

   normal_mask = psi*0.
   normal_mask[0,      *,      0] = 1.
   normal_mask[0,      *,sz[3]-1] = 1.
   normal_mask[0,      0,      *] = 1.
   normal_mask[0,sz[2]-1,      *] = 1.

   normal_deriv = psi*0.
   normal_deriv[0,      *,      0] = -3.14159625/2.     ; bottom
   normal_deriv[0,      *,sz[3]-1] =  3.14159625/2.     ; top
   normal_deriv[0,      0,      *] =  3.14159625        ; left
   normal_deriv[0,sz[2]-1,      *] =  0.                ; right
   normal_deriv = $
     (psix*cos(normal_deriv) + psiz*sin(normal_deriv))*normal_mask

   normal_deriv = normal_deriv gt 0

   psi_bound = psi*normal_mask + (1-normal_mask)*1e10

   psilim = min(psi_bound-psi0, i, /absolute)
   psilim = psi_bound[i]
   print, "Flux at limiter", psilim

   ; flux at separatrix
   sz = size(xpoint)
   if(sz[0] gt 0) then begin
       xfluxes = fltarr(sz[0])
       for i=0, sz[0]-1 do xfluxes[i] = psi[0,xpoint[0,i],xpoint[1,i]]
       psix = min(xfluxes-psi0, i, /absolute)
       psix = xfluxes[i]
       print, "Flux at separatrix:", psix

       if(abs(psix-psi0) gt (psilim-psi0)) then begin
           print, "Plasma is limited."
       endif else begin
           print, "Plasma is diverted."
           psilim = psix
       endelse
   endif
   
   return, psilim
end


pro plot_flux_average, name, time, filename=filename
   
   if(n_elements(filename) eq 0) then filename='C1.h5'

   psi = read_field('psi', slice=time, filename=filename, p=101)
   if(n_elements(psi) le 1) then return

   field = read_field(name, slice=time, filename=filename, p=101, time=t)

   fa = flux_average(psi, field, flux=flux)

   if(n_elements(title) eq 0) then begin
       if(t gt 0) then begin
           title = "!12<!8" + translate(name) + $
             string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!N!6)!12>!3")', t)
       endif else begin
           title = "!12<!8" + translate(name) + $
             string(FORMAT='("!6(!8t!6 = ",G0,")!3!12>")', t)
       endelse
   endif

   plot, flux[0,*], fa[0,*], xtitle='!7w!3', $
     ytitle="!12<!8" + translate(name) + "!12>!3", title=title
end

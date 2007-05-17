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


function eval_field, field, mesh, x=xi, y=yi, points=p, operation=op, $
                     filename=filename

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nelms = mesh.nelms._data 

   if(n_elements(p) eq 0) then p = 100
   if(n_elements(op) eq 0) then op = 1
   
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
                   result[index[0], index[1]] = eval(field, localpos, i, op=op)
               endif
               
               pos[0] = pos[0] + dx
               index[0] = index[0] + 1
           end
           pos[1] = pos[1] + dy
           index[1] = index[1] + 1
       end
   end

   xi = findgen(p)*(xrange[1]-xrange[0])/(p-1.) + xrange[0] $
     + read_parameter("xzero", filename=filename)
   yi = findgen(p)*(yrange[1]-yrange[0])/(p-1.) + yrange[0] $
     + read_parameter("zzero", filename=filename)

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
                     x=x, y=y, points=pts
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
                                      x=x, y=y, op=op, filename=filename)
   end

   h5f_close, file_id

   if(1 eq strcmp('j', name)) then data = -data

   return, data
end


pro plot_field, name, time, points=p, filename=filename, mesh=plotmesh, $
                mcolor=mc, lcfs=lcfs, title=title, $
                maskrange=maskrange, maskfield=maskfield,_EXTRA = ex

   field = read_field(name, slices=time, mesh=mesh, filename=filename, $,
                      time=t, x=x, y=y, points=p)
   if(n_elements(field) le 1) then return

   if(n_elements(maskrange) eq 2) then begin
       if(strcmp(name, maskfield) eq 1) then begin
           psi = field
       endif else begin
           psi = read_field(maskfield, slices=time, mesh=mesh, $
                            filename=filename, points=p)
       endelse
       mask = (psi ge maskrange[0]) and (psi le maskrange[1])
       field = mask*field
   endif

   print, "Plotting..."
   !x.title = '!8x/L!3'
   !y.title = '!8y/L!3'
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
           plot_lcfs, time, color=130
       endif else begin
           plot_lcfs, time, color=130, val=maskrange[0], psi=psi, x=x, y=y
           plot_lcfs, time, color=130, val=maskrange[1], psi=psi, x=x, y=y
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
   scalars = read_scalars(filename=filename)

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

;  Account for energy terms not included in the physical model
   if(nv le 2) then begin
       dissipated = E_D + E_H
   endif else begin
       dissipated = scalars.E_KPH._data + scalars.E_KTH._data
   endelse

   vloop = read_parameter('vloop', filename=filename)
   eloop = vloop * scalars.toroidal_current._data / (2.*3.14159625)
   dissipated = dissipated - eloop

   Error = E - E[0]
   total_lost = fltarr(n_elements(Error))
   total_lost[0] = 0.
   for i=1, n_elements(Error)-1 do begin
       dt = scalars.time._data[i]-scalars.time._data[i-1]
       total_lost[i] = total_lost[i-1] + $
         dt*(dissipated[i-1] + dissipated[i])/2.
   endfor

   Error = Error - total_lost

   return, E
end

pro plot_energy, filename=filename, diff=diff, norm=norm, ylog=ylog

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = read_parameter("numvar", filename=filename)

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

;  Account for energy terms not included in the physical model
   if(nv le 2) then begin
       dissipated = E_D + E_H
   endif else begin
       dissipated = scalars.E_KPH._data + scalars.E_KTH._data
   endelse

   vloop = read_parameter('vloop', filename=filename)
   eloop = vloop * scalars.toroidal_current._data / (2.*3.14159625)
   dissipated = dissipated - eloop

   Error = E - E[0]
   total_lost = fltarr(n_elements(Error))
   total_lost[0] = 0.
   for i=1, n_elements(Error)-1 do begin
       dt = scalars.time._data[i]-scalars.time._data[i-1]
       total_lost[i] = total_lost[i-1] + $
         dt*(dissipated[i-1] + dissipated[i])/2.
   endfor

   Error = Error - total_lost

   print, Error / E

   if(keyword_set(norm)) then begin
       E = E - E[0]
       E_K = E_K - E_K[0]
       E_M = E_M - E_M[0]
       E_D = E_D - E_D[0]
       E_H = E_H - E_H[0]
   endif 

   if(keyword_set(diff)) then begin
       E = deriv(scalars.time._data,E)
       E_K = deriv(scalars.time._data, E_K)
       E_M = deriv(scalars.time._data, E_M)
       !y.title = '!6d!8E!6/d!8t!3'
   endif else begin
       !y.title = '!6 Energy!3'
   endelse

   !y.range=[min([E,E_K,E_M,-E_D,-E_H,Error]), $
             max([E,E_K,E_M,-E_D,-E_H,Error])]
   if(keyword_set(ylog) and (!y.range[0] lt 0)) then $
     !y.range[0] = !y.range[1]*1e-8
   if(keyword_set(ylog)) then !y.range[1] = !y.range[1]*10.

   plot, scalars.time._data, E, ylog=ylog
   oplot, scalars.time._data, E_K, color=dc, linestyle = 2
   oplot, scalars.time._data, E_M, color=2*dc, linestyle = 2
   oplot, scalars.time._data, -E_D, color=3*dc, linestyle = 1
   oplot, scalars.time._data, -E_H, color=4*dc, linestyle = 1
   oplot, scalars.time._data, Error, color=5*dc, linestyle = 3

   plot_legend, ['Total', 'Kinetic', 'Magnetic', $
                 'Diffusive', 'Hyper-Diffusive', '|Error|'] , $
     color=[-1,1,2,3,4,5,6,7,8,9]*dc, ylog=ylog, $
     linestyle = [0,2,2,1,1,3]
end


pro plot_field_mpeg, fieldname, mpegame=mpegname, range=range, points=pts, $
                     _EXTRA=ex
    if(n_elements(mpegname) eq 0) then mpegname = fieldname + '.mpeg'
    if(n_elements(pts) eq 0) then pts = 50

    nt = get_parameters("ntime")

    data = read_field(fieldname, mesh=mesh, x=x, y=y, points=pts)

    contour_and_legend_mpeg, mpegname, data, x, y, _EXTRA=ex
end


pro plot_lcfs, time, color=color, val=psival, psi=psi, x=x, y=y
    pts = 201

    if(n_elements(psi) eq 0) then begin
        psi = read_field('psi', slice=time, points=pts, x=x, y=y)
    endif

    ; if psival not passed, choose limiter value
    if(n_elements(psival) eq 0) then begin
        xlim = read_parameter("xlim")
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


pro plot_scalar, scalarname, filename=filename, names=names, $
                 _EXTRA=extra, overplot=overplot, $
                 ylog=ylog, xlog=xlog, left=left

  if(n_elements(filename) eq 0) then filename='C1.h5'

  !x.title = '!8t!6 (!7s!D!8A!N!6)!3'

  nfiles = n_elements(filename)
  if(nfiles gt 1) then begin
      colors = colors(nfiles)

      plot_scalar, scalarname, filename=filename[0], $
        color=colors[0], _EXTRA=extra, ylog=ylog, xlog=xlog
      for i=1, nfiles-1 do begin
          plot_scalar, scalarname, filename=filename[i], $
            /overplot, color=colors[i], _EXTRA=extra, ylog=ylog, xlog=xlog
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
    if (strcmp("beta", scalarname, /fold_case) eq 1) then begin
      if(nv lt 3) then begin
          print, "Must be numvar = 3 for beta calculation"
          return
      endif
      gamma = read_parameter('gam', filename=filename[i])
      data = 2.*(gamma-1.)*s.E_P._data/(s.E_MP._data + s.E_MT._data)
      title = '!7b!3'
      ytitle = '!7b!3'
  endif else if $
    (strcmp("poloidal beta", scalarname, /fold_case) eq 1) or $
    (strcmp("bp", scalarname, /fold_case) eq 1) then begin
      nv = read_parameter("numvar", filename=filename[0])
      if(nv lt 3) then begin
          print, "Must be numvar = 3 for beta calculation"
          return
      endif
      gamma = read_parameter('gam', filename=filename[i])
      data = 2.*(gamma-1.)*s.E_P._data/s.toroidal_current._data^2
      title = '!7b!D!8p!N!3'
      ytitle = '!7b!D!8p!N!3'
  endif else if $
    (strcmp("kinetic energy", scalarname, /fold_case) eq 1) or $
    (strcmp("ke", scalarname, /fold_case) eq 1)then begin
      nv = read_parameter("numvar", filename=filename[0])
       data = s.E_KP._data 
       if(nv ge 2) then data = data + s.E_KT._data 
       if(nv ge 3) then data = data + s.E_K3._data
       title = '!6Kinetic Energy!3'
       ytitle = '!6KE (!8B!D0!N!U!62!N/4!7p!6)!3'
  endif else begin
  endelse
  
  if(keyword_set(overplot)) then begin
      oplot, s.time._data, data, _EXTRA=extra
  endif else begin
      plot, s.time._data, data, $
        title=title, ytitle=ytitle, _EXTRA=extra, ylog=ylog, xlog=xlog
  endelse
end

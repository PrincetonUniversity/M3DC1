function read_field, name, x, y, t, slices=slices, mesh=mesh, $
                     filename=filename, points=pts, mask=mask, $
                     rrange=xrange, zrange=yrange,equilibrium=equilibrium, $
                     h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                     diff=diff, operation=op, complex=complex, fac=fac, $
                     linear=linear, last=last, average=average, linfac=linfac,$
                     dpsi=dpsi, symbol=symbol, units=units, cgs=cgs, mks=mks, $
                     real=real, imaginary=imag, edge_val=edge_val, phi=phi0, $
                     time=realtime, abs=abs, phase=phase, dimensions=d, $
                     flux_average=flux_av, rvector=rvector, zvector=zvector, $
                     yvector=yvector, taverage=taverage, sum=sum, $
                     tpoints=nphi, logical=logical, map_r=map_r, map_z=map_z, $
                     is_nonlinear=is_nonlinear, outval=mask_val, wall_mask=wall_mask

  if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(slices) ne 0) then time=slices else time=0
   is_nonlinear = 0

   icomplex = read_parameter('icomplex', filename=filename[0])
   itor = read_parameter('itor', filename=filename[0])
   rzero = read_parameter('rzero', filename=filename[0])
   if(itor eq 1) then begin
      period = 2.*!pi
   endif else begin
      period = 2.*!pi*rzero
   end
   if(n_elements(pts) eq 0) then pts=200.
   if(n_elements(phi0) eq 0) then phi0=0.
   if(n_elements(nphi) eq 0) then begin
      nphi = n_elements(phi0)
   endif else begin
      phi0 = findgen(nphi)/nphi * period
      if(itor eq 1) then phi0 = phi0*180./!pi
   end
   print, 'phi0 = ', phi0

   if(nphi gt 1) then begin
      if(icomplex eq 1) then begin
         ntor = read_parameter('ntor', filename=filename)
         data = complexarr(nphi,pts,pts)
         data[0,*,*] = read_field(name, x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, fac=fac, $
                        rrange=xrange, zrange=yrange, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        diff=diff, operation=op, dimensions=d, $
                        last=last,symbol=symbol,units=units, $
                        cgs=cgs, mks=mks, time=realtime, $
                        rvector=rvector, zvector=zvector, yvector=yvector,$
                        phi=phi0[0], wall_mask=wall_mask, /linear, /complex, $
                        logical=logical)
         for i=1, nphi-1 do begin
            phi_rad = (phi0[i]-phi0[0])*!pi/180.
            data[i,*,*] = data[0,*,*]*exp(complex(0.,ntor)*phi_rad)
         end
         data = real_part(data)

         if(not keyword_set(linear) and time ge 0) then begin
            data0 = read_field(name, x, y, t, slices=-1, $
                        filename=filename, points=pts, fac=fac, $
                        rrange=xrange, zrange=yrange, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        diff=diff, operation=op, dimensions=d, $
                        symbol=symbol,units=units, $
                        cgs=cgs, mks=mks, time=realtime, $
                        rvector=rvector, zvector=zvector, yvector=yvector,$
                        phi=0., wall_mask=wall_mask, logical=logical)
            for i=0, nphi-1 do begin
               data[i,*,*] = data[i,*,*] + data0
            end
         end

         return, real_part(data)
      endif else begin
         data = fltarr(nphi,pts,pts)
         for i=0, nphi-1 do begin
            data[i,*,*] = $
               read_field(name, x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, fac=fac, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                          diff=diff, operation=op, dimensions=d, $
                          linear=linear, last=last,symbol=symbol,units=units, $
                          cgs=cgs, mks=mks, time=realtime, $
                          rvector=rvector, zvector=zvector, yvector=yvector,$
                          phi=phi0[i], wall_mask=wall_mask, logical=logical)
         end
         return, data
      end
   end

   if(keyword_set(taverage)) then begin
       data = 0
       if(taverage eq 1) then taverage=16
       phi = period*findgen(taverage) / taverage
       if(itor eq 1) then phi = phi*180./!pi
       for i=0, taverage-1 do begin
           data = data + $
             read_field(name, x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, fac=fac, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        diff=diff, operation=op, dimensions=d, $
                        linear=linear, last=last,symbol=symbol,units=units, $
                        cgs=cgs, mks=mks, time=realtime, $
                        rvector=rvector, zvector=zvector, yvector=yvector,$
                        phi=phi[i], wall_mask=wall_mask,logical=logical)
       end
       data = data/taverage
       return, data
   end

   if(keyword_set(average) or keyword_set(sum)) then begin
       if(n_elements(filename) gt 1) then begin
           n = n_elements(filename)
           if(n_elements(time) eq 1) then time=replicate(time,n)
       endif else if(n_elements(time) gt 1) then begin
           n = n_elements(time)
           if(n_elements(filename) eq 1) then filename=replicate(filename,n)
       endif else n = 1
       if(n_elements(fac) eq 0) then fac=1.
       if(n_elements(fac) lt n) then fac=replicate(fac,n)
       if(n_elements(linfac) eq 0) then linfac=1.
       if(n_elements(linfac) lt n) then linfac=replicate(linfac,n)

       data = 0
       for i=0, n-1 do begin
           data = data + $
             read_field(name, x, y, t, slices=time[i], mesh=mesh, $
                        filename=filename[i], points=pts, $
                        equilibrium=equilibrium, linfac=linfac[i], $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        operation=op, dimensions=d, fac=fac[i], $
                        linear=linear, last=last,symbol=symbol,units=units, $
                        cgs=cgs, mks=mks, phi=phi0, time=realtime, $
                        rvector=rvector, zvector=zvector, yvector=yvector, $
                       wall_mask=wall_mask,logical=logical)
       end
       if(keyword_set(average)) then data = data/n
       return, data
   end
   if(keyword_set(diff)) then begin
       if(n_elements(filename) gt 1) then begin
           n = n_elements(filename)
           if(n_elements(time) eq 1) then time=replicate(time,n)
       endif else if(n_elements(time) gt 1) then begin
           n = n_elements(time)
           if(n_elements(filename) eq 0) then filename='C1.h5'
           if(n_elements(filename) eq 1) then filename=replicate(filename,n)
       endif else n = 1

       data = 0
       for i=0, n-1 do begin
           data = data + $
             read_field(name, x, y, t, slices=time[i], mesh=mesh, $
                        filename=filename[i], points=pts, fac=fac, $
                        rrange=xrange, zrange=yrange, mask=mask, $
                        h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                        operation=op, complex=complex, dimensions=d, $
                        linear=linear, last=last,symbol=symbol,units=units, $
                       cgs=cgs, mks=mks, phi=phi0, time=realtime, $
                       rvector=rvector, zvector=zvector, yvector=yvector, $
                       wall_mask=wall_mask,logical=logical) $
             *((-1)^i)
       end

       symbol = '!7D!X' + symbol

       return, data
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(filename) gt 1) then filename=filename[0]
   if(n_elements(pts) eq 0) then pts = 200
   if(n_elements(op) eq 0) then op = 1

   if(hdf5_file_test(filename) eq 0) then return, 0

   version = read_parameter("version", filename=filename)
   nt = read_parameter("ntime", filename=filename)
   nv = read_parameter("numvar", filename=filename)
   ntor = read_parameter("ntor", filename=filename)
   version = read_parameter('version', filename=filename)
   ivform = read_parameter('ivform', filename=filename)
   i3d = read_parameter('3d', filename=filename)
   if(version eq 0) then begin
       xzero = read_parameter("xzero", filename=filename)
       zzero = read_parameter("zzero", filename=filename)
   endif else begin
       xzero = 0.
       zzero = 0.
   endelse
   ilin = read_parameter('linear', filename=filename)
   isubeq = read_parameter('eqsubtract', filename=filename)
   extsubtract = read_parameter('extsubtract', filename=filename)

   if(keyword_set(last)) then time = nt-1
   if(keyword_set(equilibrium)) then begin
       if(ilin eq 1) then time=-1
       if(isubeq eq 1) then linear = 0
   end

   if(time ge nt) then begin
       print, "Error: there are only ", nt-1, " time slices."
       return, 0
   endif

   realtime = get_slice_time(filename=filename, slice=time)

   data = fltarr(nphi, pts, pts)
   if(isubeq eq 1) then base = fltarr(pts,pts)

   d = dimensions()
   symbol=name
 
   print, '**********************************************************'
   print, 'Reading ', name, ' at timeslice ', time
   print, 'From file ', filename
   print, 'Eqsubtract? ', isubeq
   print, 'sum = ', keyword_set(sum)
   if(n_elements(fac) ne 0) then print, 'fac = ', fac
   print, string(form='(" linear=",I0,"; pts=",I0,";' + $
                 'equilibrium=",I0,"; complex=",I0,"; op=",I0)', $
                 keyword_set(linear), pts, keyword_set(equilibrium), $
                 keyword_set(complex), op)

   if(isubeq eq 1 and (not keyword_set(linear)) and (time ge 0)) $
     then begin
       data1 = read_field(name,x,y,t, slices=time, mesh=mesh, fac=fac, $
                          filename=filename, points=pts, mks=mks, cgs=cgs, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                          diff=diff, operation=op, linfac=linfac, $
                          /linear, last=last,symbol=symbol, $
                          units=units, dimensions=d, phi=phi0, $
                          rvector=rvector, zvector=zvector, $
                          yvector=yvector, is_nonlinear=isnl, $
                         wall_mask=wall_mask,logical=logical)

       if(n_elements(data1) le 1) then begin
          print, 'Perturbed field not found.'
          isnl = 0
       endif
       if(isnl eq 1) then begin
          data = data1
       endif else begin
          data0 = read_field(name,x,y,t0, slices=-1, $
                             filename=filename, points=pts, fac=fac, $
                             rrange=xrange, zrange=yrange, complex=0, $
                             h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                             diff=diff, operation=op, mask=mask, $
                             symbol=symbol, mks=mks, cgs=cgs, $
                             units=units, dimensions=d, $
                             rvector=rvector, zvector=zvector, $
                             yvector=yvector, wall_mask=wall_mask, $
                             logical=logical)
          if(n_elements(data1) le 1) then begin
             data = data0
          endif else begin
             data = data0 + data1
          endelse
      endelse
      print, '**********************************************************'
      return, data
   endif
   
   ; check if this is a primitive field
   file_id = h5f_open(filename)
   time_group_id = h5g_open(file_id, time_name(time))
   mesh = h5_parse(time_group_id, 'mesh', /read_data)
               
   field_group_id = h5g_open(time_group_id, 'fields')             
   nmembers = h5g_get_nmembers(time_group_id, 'fields')
   match = 0
   for m=0, nmembers-1 do begin
       thisname = h5g_get_member_name(time_group_id,'fields',m)
       if(strcmp(thisname, name, /fold_case) eq 1) then begin
           name = thisname
           match = 1
           break
       endif
   end
   
   if(match eq 1) then begin

       if(keyword_set(complex)) then begin
           h5g_close, field_group_id
           h5g_close, time_group_id
           h5f_close, file_id
           print, '  reading complex field.', ntor

           data_r = read_field(name,x,y,t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, complex=0, $
                               h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                               diff=diff, operation=op, mask=mask, $
                               /linear, last=last,symbol=symbol, $
                               units=units, dimensions=d, $
                               equilibrium=equilibrium, wall_mask=wall_mask)

           data_i = read_field(name+'_i',x,y,t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, complex=0, $
                               h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                               diff=diff, operation=op, $
                               /linear, last=last,symbol=symbol, $
                               units=units, dimensions=d, $
                               equilibrium=equilibrium, wall_mask=wall_mask)
           data = complex(data_r, data_i)

             ; evaluate at phi0
             ; it is okay to do this here since complex cases are always linear
           if(n_elements(phi0) ne 0) then begin
              data = data*exp(complex(0.,ntor)*phi0[0]*!pi/180.)
           end
       endif else if (keyword_set(linear) and isubeq eq 0) then begin
           h5g_close, field_group_id
           h5g_close, time_group_id
           h5f_close, file_id
           print, '  calculating perturbed part of eqsubtract=0 field.', ntor

           data0 = read_field(name,x,y,t, slices=-1, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, $
                               h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                               diff=diff, operation=op, mask=mask, $
                               symbol=symbol, phi=phi0, $
                               units=units, dimensions=d, wall_mask=wall_mask, $
                               logical=logical)
           data1 = read_field(name,x,y,t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, $
                               h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                               diff=diff, operation=op, $
                               last=last,symbol=symbol, phi=phi0, $
                               units=units, dimensions=d, $
                               equilibrium=equilibrium, wall_mask=wall_mask, $
                               logical=logical)
           data = data1 - data0
       endif else begin
           print, '  reading real field'

           field = h5_parse(field_group_id, name, /read_data)
               
           time_id = h5a_open_name(time_group_id, "time")
           t = h5a_read(time_id)
           h5a_close, time_id
           h5g_close, field_group_id
           h5g_close, time_group_id
           h5f_close, file_id
               
           if(n_elements(phi0) eq 0) then begin
              phi0_rad=0.
           endif else begin
              if(itor eq 1) then phi0_rad = phi0*!pi/180. else phi0_rad = phi0
           end
  
           data[0,*,*] = $
             eval_field(field._data, mesh, points=pts, $
                        r=x, z=y, op=op, filename=filename, $
                        xrange=xrange, yrange=yrange, mask=mask, $
                        phi=phi0_rad, wall_mask=wall_mask, edge_val=edge_val)
           symbol = field_data(name, units=d, itor=itor, filename=filename)

           igeometry = read_parameter("igeometry", filename=filename)

           if(igeometry gt 0 and not keyword_set(logical)) then begin
              if(n_elements(map_r) ne n_elements(data) or $
                 n_elements(map_z) ne n_elements(data)) then begin
                 rst = read_field('rst',x,z,t,points=pts, filename=filename, $
                                  phi=phi0, /logical, slice=-1)
                 zst = read_field('zst',x,z,t,points=pts, filename=filename, $
                                  phi=phi0, /logical, slice=-1)
                 create_map, rst, zst, r=x, z=y, ix=map_r, iy=map_z, mask=mask
              end
              data[0,*,*] = map_field(data, map_r, map_z, $
                                      mask=mask, outval=edge_val)
           end

           if(version lt 5 and isubeq eq 1 and time ge 0 and $
              ((strcmp('te', name, /fold_case) eq 1) or $
               (strcmp('te_i', name, /fold_case) eq 1))) then begin
               zeff = read_parameter('zeff',filename=filename)
               data = data / zeff
               print, 'Correcting bug in linear Te for version < 4'
           end
       endelse

       
   endif else begin
       h5g_close, field_group_id
       h5g_close, time_group_id
       h5f_close, file_id

       print, '  reading composite field'

   if(strcmp('zero', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       data = 0.*psi
       symbol = ''


    endif else if((strcmp('fp', name, /fold_case) eq 1) and $
                 (version lt 35)) then begin

       print, 'Calculating fp from f'
       
       if(i3d eq 1) then begin
          data = read_field(name,x,y,t, slices=time, mesh=mesh, fac=fac, $
                            filename=filename, points=pts, mks=mks, cgs=cgs, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                            diff=diff, operation=op+10, linfac=linfac, $
                            last=last,symbol=symbol, $
                            units=units, dimensions=d, phi=phi0, $
                            rvector=rvector, zvector=zvector, $
                            yvector=yvector, is_nonlinear=isnl, $
                            wall_mask=wall_mask,logical=logical)
       endif else begin
          data = read_field('f',x,y,t, slices=time, mesh=mesh, fac=fac, $
                            filename=filename, points=pts, mks=mks, cgs=cgs, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                            diff=diff, operation=op, linfac=linfac, $
                            last=last,symbol=symbol, $
                            units=units, dimensions=d, phi=phi0, $
                            rvector=rvector, zvector=zvector, $
                            yvector=yvector, is_nonlinear=isnl, $
                            wall_mask=wall_mask,logical=logical)
          data = data*complex(0.,ntor)
       end

    endif else if((strcmp('fp_plasma', name, /fold_case) eq 1) and $
                 (version lt 35)) then begin

       print, 'Calculating fp_plasma from f_plasma'
       
       if(i3d eq 1) then begin
          data = read_field('f_plasma',x,y,t, slices=time, mesh=mesh, fac=fac, $
                            filename=filename, points=pts, mks=mks, cgs=cgs, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                            diff=diff, operation=op+10, linfac=linfac, $
                            last=last,symbol=symbol, $
                            units=units, dimensions=d, phi=phi0, $
                            rvector=rvector, zvector=zvector, $
                            yvector=yvector, is_nonlinear=isnl, $
                            wall_mask=wall_mask,logical=logical)
       endif else begin
          data = read_field('f_plasma',x,y,t, slices=time, mesh=mesh, fac=fac, $
                            filename=filename, points=pts, mks=mks, cgs=cgs, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                            diff=diff, operation=op, linfac=linfac, $
                            last=last,symbol=symbol, $
                            units=units, dimensions=d, phi=phi0, $
                            rvector=rvector, zvector=zvector, $
                            yvector=yvector, is_nonlinear=isnl, $
                            wall_mask=wall_mask,logical=logical)
          data = data*complex(0.,ntor)
       end

    endif else if((strcmp('fp_ext', name, /fold_case) eq 1) and $
                 (version lt 35)) then begin

       print, 'Calculating fp_ext from f_ext'
       
       if(i3d eq 1) then begin
          data = read_field('f_ext',x,y,t, slices=time, mesh=mesh, fac=fac, $
                            filename=filename, points=pts, mks=mks, cgs=cgs, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                            diff=diff, operation=op+10, linfac=linfac, $
                            last=last,symbol=symbol, $
                            units=units, dimensions=d, phi=phi0, $
                            rvector=rvector, zvector=zvector, $
                            yvector=yvector, is_nonlinear=isnl, $
                            wall_mask=wall_mask,logical=logical)
       endif else begin
          data = read_field('f_ext',x,y,t, slices=time, mesh=mesh, fac=fac, $
                            filename=filename, points=pts, mks=mks, cgs=cgs, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            h_symmetry=h_symmetry, v_symmetry=v_symmetry, $
                            diff=diff, operation=op, linfac=linfac, $
                            last=last,symbol=symbol, $
                            units=units, dimensions=d, phi=phi0, $
                            rvector=rvector, zvector=zvector, $
                            yvector=yvector, is_nonlinear=isnl, $
                            wall_mask=wall_mask,logical=logical)
          data = data*complex(0.,ntor)
       end
       
   ;==========================================
   ; local_beta = 2*P/B^2
   ;==========================================
   endif else if(strcmp('beta', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(n_elements(psi) le 1) then return, 0

       I = read_field('I',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)     
       P = read_field('P',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       b2 = (s_bracket(psi,psi,x,y) + i^2)/r^2

       data = 2.*P/b2
       symbol = '!7b!X'

   ;===========================================
   ; toroidal field
   ;===========================================
   endif else if(strcmp('toroidal field', name, /fold_case) eq 1 or $
                 strcmp('by', name, /fold_case) eq 1) then begin
       
       
       I = read_field('I',x,y,t,slices=time,mesh=mesh,filename=filename,$
                      points=pts, rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, phi=phi0, $
                      wall_mask=wall_mask,outval=mask_val,mask=mask)

       if(extsubtract eq 1 and version lt 8) then begin
           I = I + read_field('I_ext', x, y, t, mesh=mesh, $
                              filename=filename, points=pts, slices=time, $
                              rrange=xrange, zrange=yrange, complex=complex, $
                              linear=linear, mask=mask, phi=phi0, $
                              wall_mask=wall_mask,outval=mask_val)
       end

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
   
       data = I/r
       symbol = '!8B!D!9P!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; toroidal velocity
   ;===========================================
   endif else if(strcmp('toroidal velocity', name, /fold_case) eq 1) or $
     (strcmp('vy', name, /fold_case) eq 1) then begin
       
       v = read_field('V',x,y,t,slices=time, mesh=mesh, filename=filename, $
                        points=pts,rrange=xrange,zrange=yrange, $
                     linear=linear,complex=complex,phi=phi0)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
   
       if(ivform eq 0) then begin
           data = v/r
       endif else if(ivform eq 1) then begin
           data = v*r
       endif
       symbol = '!8u!D!9P!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; toroidal velocity shear
   ;===========================================
   endif else if(strcmp('toroidal velocity shear', name, /fold_case) eq 1) or $
     (strcmp('vzp', name, /fold_case) eq 1) then begin
       
       v = read_field('V',x,y,t,slices=time, mesh=mesh, filename=filename, $
                        points=pts,rrange=xrange,zrange=yrange, linear=linear)
       psi = read_field('psi',x,y,t,slices=time, filename=filename,$
                        points=pts,rrange=xrange,zrange=yrange,/equilibrium)


       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
   
       if(ivform eq 0) then begin
           vz = v/r
       endif else if(ivform eq 1) then begin
           vz = v*r
       endif

       data = s_bracket(vz,psi,x,y)/sqrt(s_bracket(psi,psi,x,y))
       symbol = "!8u!D!9P!N'!X"
       d = dimensions(/t0, _EXTRA=extra)

   ;===========================================
   ; thermal velocity
   ;===========================================
   endif else if(strcmp('vt_i', name, /fold_case) eq 1) or $
     (strcmp('vti', name, /fold_case) eq 1) then begin

       Ti = read_field('Ti',x,y,t,slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange)
         
       data = sqrt(2.*Ti)
       symbol = '!8v!Dti!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; sound speed
   ;===========================================
   endif else if(strcmp('sound speed', name, /fold_case) eq 1) or $
     (strcmp('cs', name, /fold_case) eq 1) then begin

       P = read_field('P',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       den = read_field('den',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       gam = read_parameter('gam', filename=filename)
  
       data = sqrt(gam*P/den)
       symbol = '!8c!Ds!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; Mach number
   ;===========================================
   endif else if(strcmp('mach', name, /fold_case) eq 1) or $
     (strcmp('m', name, /fold_case) eq 1) then begin

       cs = read_field('cs',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)

       phi = read_field('phi',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       V = read_field('V',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       chi = read_field('chi',x,y,t,slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       if(ivform eq 0) then begin  
           v2 = s_bracket(phi,phi,x,y)/r^2 $
             + v^2/r^2 + s_bracket(chi,chi,x,y) $
             + 2.*a_bracket(chi,phi,x,y)/r
       endif else if(ivform eq 1) then begin
           v2 = r^2*s_bracket(phi,phi,x,y) $
             + r^2*v^2 + s_bracket(chi,chi,x,y)/r^4 $
             + 2.*a_bracket(chi,phi,x,y)/r
       endif
  
       data = sqrt(v2)/cs
       symbol = '!8M!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; ExB drift/velocity
   ;===========================================
   endif else if(strcmp('ExB_x', name, /fold_case) eq 1) then begin

       ;Ex = read_field('E_R', x, y, t, slices=time, mesh=mesh, linear=linear, $
       ;                 filename=filename, points=pts, complex=complex, $
       ;                 rrange=xrange, zrange=yrange, phi=phi0)
       Ey = read_field('E_phi', x, y, t, slices=time, mesh=mesh,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ez = read_field('E_Z', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)

       Bx = read_field('Bx', x, y, t, slices=time,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       By = read_field('By', x, y, t, slices=time,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Bz = read_field('Bz', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       
       B2 = Bx^2 + By^2 + Bz^2

       data = (Ey*Bz - Ez*By)/ B2
       d = dimensions(/v0)
       symbol = '!6V!DExB,R!X'

    endif else if(strcmp('ExB_y', name, /fold_case) eq 1) then begin

       Ex = read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       ;Ey = read_field('E_phi', x, y, t, slices=time, mesh=mesh,$
       ;                 filename=filename, points=pts, $
       ;                 rrange=xrange, zrange=yrange, phi=phi0)
       Ez = read_field('E_Z', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)

       Bx = read_field('Bx', x, y, t, slices=time,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       By = read_field('By', x, y, t, slices=time,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Bz = read_field('Bz', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       
       B2 = Bx^2 + By^2 + Bz^2

       data = (Ez*Bx - Ex*Bz)/ B2
       d = dimensions(/v0)
       symbol = '!6V!DExB,!7u!X'

    endif else if(strcmp('ExB_z', name, /fold_case) eq 1) then begin

       Ex = read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ey = read_field('E_phi', x, y, t, slices=time, mesh=mesh,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       ;Ez = read_field('E_Z', x, y, t, slices=time, mesh=mesh, $
       ;                 filename=filename, points=pts, $
       ;                 rrange=xrange, zrange=yrange, phi=phi0)

       Bx = read_field('Bx', x, y, t, slices=time,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       By = read_field('By', x, y, t, slices=time,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Bz = read_field('Bz', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       
       B2 = Bx^2 + By^2 + Bz^2

       data = (Ex*By - Ey*Bx)/ B2
       d = dimensions(/v0)
       symbol = '!6V!DExB,Z!X'


   ;===========================================
   ; electron density
   ;===========================================
  endif else if(strcmp('electron density', name, /fold_case) eq 1) or $
    (strcmp('ne', name, /fold_case) eq 1) then begin

     print, 'WARNING: calculating ne using ne = z_ion * ni'

     ni = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, linear=linear, $
                   complex=complex, op=op)

     if(version ge 23) then  begin
        z_ion = read_parameter('z_ion', filename=filename)
     endif else begin
        z_ion = read_parameter('zeff', filename=filename)
     end
     data = ni*z_ion

      symbol = '!8n!De!N!X'
      d = dimensions(/n0, _EXTRA=extra)


   ;===========================================
   ; electron temperature
   ;===========================================
  endif else if(strcmp('electron temperature', name, /fold_case) eq 1) or $
    (strcmp('te', name, /fold_case) eq 1) then begin

      Pe0 = read_field('Pe', x, y, t, slices=time, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, linear=linear, $
                       /equilibrium)

      n0 = read_field('ne', x, y, t, slices=time, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, linear=linear, $
                      /equilibrium)


      if(keyword_set(isubeq eq 1) and keyword_set(linear) and time ge 0) $
        then begin
          Pe1 = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, $
                           rrange=xrange, zrange=yrange, /linear, $
                           complex=complex, phi=phi0)
          
          n1 = read_field('ne', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, $
                          rrange=xrange, zrange=yrange, /linear, $
                          complex=complex, phi=phi0)
          data = pe1/n0 - pe0*n1/n0^2
      endif else begin
          data = pe0/n0
      endelse

;       if(keyword_set(linear) and (isubeq eq 1) and (time ge 0)) then begin
;           Pe0 = read_field('Pe', x, y, t, slices=time, $
;                            filename=filename, points=pts, $
;                            rrange=xrange, zrange=yrange, /equilibrium)

;           n0 = read_field('ne', x, y, t, slices=time, $
;                           filename=filename, points=pts,  $
;                           rrange=xrange, zrange=yrange, /equilibrium)

;           data = pe1/n0 - pe0*n1/n0^2
;       endif else data = pe1/n1

      symbol = '!8T!De!N!X'
      d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; delta_W
   ;===========================================
   endif else if(strcmp('delta_W', name, /fold_case) eq 1) then begin

      jr0 = read_field('jx',x,y,t,filename=filename,points=pts,$
                       slice=-1, rrange=xrange, zrange=yrange, mask=mask)
      jphi0 = read_field('jy',x,y,t,filename=filename,points=pts,$
                       slice=-1, rrange=xrange, zrange=yrange, mask=mask)
      jz0 = read_field('jz',x,y,t,filename=filename,points=pts,$
                       slice=-1, rrange=xrange, zrange=yrange, mask=mask)
      jr1 = read_field('jx',x,y,t,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)
      jphi1 = read_field('jy',x,y,t,filename=filename,points=pts,$
                         slice=time, rrange=xrange, zrange=yrange, $
                         complex=complex, linear=linear, phi=phi0)
      jz1 = read_field('jz',x,y,t,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)
      br0 = read_field('bx',x,y,t,filename=filename,points=pts,$
                       slice=-1, rrange=xrange, zrange=yrange)
      bphi0 = read_field('by',x,y,t,filename=filename,points=pts,$
                       slice=-1, rrange=xrange, zrange=yrange)
      bz0 = read_field('bz',x,y,t,filename=filename,points=pts,$
                       slice=-1, rrange=xrange, zrange=yrange)
      br1 = read_field('bx',x,y,t,mesh=mesh,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)
      bphi1 = read_field('by',x,y,t,mesh=mesh,filename=filename,points=pts,$
                         slice=time, rrange=xrange, zrange=yrange, $
                         complex=complex, linear=linear, phi=phi0)
      bz1 = read_field('bz',x,y,t,mesh=mesh,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)
      p1 = read_field('p',x,y,t,mesh=mesh,filename=filename,points=pts,$
                      slice=time, rrange=xrange, zrange=yrange, $
                      complex=complex, linear=linear, phi=phi0)
      pr1 = read_field('p',x,y,t,mesh=mesh,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0,op=2)
      pz1 = read_field('p',x,y,t,mesh=mesh,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0,op=3)

      xir = read_field('xi_x',x,y,t,mesh=mesh,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)
      xiphi = read_field('xi_y',x,y,t,mesh=mesh,filename=filename,points=pts,$
                         slice=time, rrange=xrange, zrange=yrange, $
                         complex=complex, linear=linear, phi=phi0)
      xiz = read_field('xi_z',x,y,t,mesh=mesh,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)

      pphi1 = complex(0,ntor)*p1

      r = radius_matrix(x,y,t)

      fr   = (jphi0*bz1 - jz0*bphi1 + jphi1*bz0 - jz1*bphi0 - pr1)
      fphi = (jz0*br1 - jr0*bz1 + jz1*br0 - jr1*bz0         - pphi1/r)
      fz   = (jr0*bphi1 - jphi0*br1 + jr1*bphi0 - jphi1*br0 - pz1)

      data = -0.5*(conj(xir)*fr + conj(xiphi)*fphi + conj(xiz)*fz)
  
      symbol = '!7d!8w!X'
      d = dimensions(/p0, _EXTRA=extra)

   ;===========================================
   ;  K (kinetic part of delta W)
   ;===========================================
   endif else if(strcmp('K', name, /fold_case) eq 1) then begin

      den0 = read_field('den',x,y,t,filename=filename,points=pts,$
                        slice=-1, rrange=xrange, zrange=yrange)
      xir = read_field('xi_x',x,y,t,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)
      xiphi = read_field('xi_y',x,y,t,mesh=mesh,filename=filename,points=pts,$
                         slice=time, rrange=xrange, zrange=yrange, $
                         complex=complex, linear=linear, phi=phi0)
      xiz = read_field('xi_z',x,y,t,mesh=mesh,filename=filename,points=pts,$
                       slice=time, rrange=xrange, zrange=yrange, $
                       complex=complex, linear=linear, phi=phi0)

      r = radius_matrix(x,y,t)

      data = 0.5*den0*(abs(xir)^2 + abs(xiphi)^2 + abs(xiz)^2)
  
      symbol = '!8K!X'
      d = dimensions(/p0, _EXTRA=extra)

   ;===========================================
   ; xi_n
   ;===========================================
   endif else if(strcmp('xi_n', name, /fold_case) eq 1) then begin
      
      g = read_gamma(filename=filename)
      if(g eq 0.) then begin
         data = read_field('vn', x, y, t, mesh=mesh, phi=phi0, slice=time, $
                           filename=filename, points=pts, complex=complex, $
                           rrange=xrange, zrange=yrange, linear=linear, $
                          phi=phi0)
      endif else begin
         print, 'xi_n = normal displacement using xi = int(v dt) method'
         vn = read_field('vn', x, y, t, mesh=mesh, phi=phi0, slice=time, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                        phi=phi0)
         data = vn / g[0]
      end

      symbol = '!7n!D!8n!N!X'
      d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; xi_x
   ;===========================================
   endif else if(strcmp('xi_x', name, /fold_case) eq 1) then begin

      print, 'xi_x = normal displacement using xi = int(v dt) method'
      vx = read_field('vx', x, y, t, mesh=mesh, phi=phi0, slice=time, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange, linear=linear)
      g = read_gamma(filename=filename)
      data = vx / g[0]     

      symbol = '!7n!D!8R!N!X'
      d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; xi_y
   ;===========================================
   endif else if(strcmp('xi_y', name, /fold_case) eq 1) then begin

      print, 'xi_y = toroidal displacement using xi = int(v dt) method'
      vy = read_field('vy', x, y, t, mesh=mesh, phi=phi0, slice=time, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange, linear=linear)
      g = read_gamma(filename=filename)
      data = vy / g[0]     

      symbol = '!7n!D!9P!N!X'
      d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; xi_z
   ;===========================================
   endif else if(strcmp('xi_z', name, /fold_case) eq 1) then begin

      print, 'xi_Z = verical displacement using xi = int(v dt) method'
      vz = read_field('vz', x, y, t, mesh=mesh, phi=phi0, slice=time, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange, linear=linear)
      g = read_gamma(filename=filename)
      data = vz / g[0]     

      symbol = '!7n!D!8Z!N!X'
      d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; displacement
   ;===========================================
   endif else if(strcmp('displacement', name, /fold_case) eq 1) then begin

       Te1 = read_field('Te', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, linear=linear, $
                       complex=complex, phi=phi0)
       

       Te0_r = read_field('Te', x, y, t, op=2, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
       Te0_z = read_field('Te', x, y, t, op=3, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
       psi0 = read_field('psi', x, y, t, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
       psi0_r = read_field('psi', x, y, t, op=2, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
       psi0_z = read_field('psi', x, y, t, op=3, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
      
       tprime = Te0_r*psi0_r + Te0_z*psi0_z

       data = -Te1/tprime*sqrt(psi0_r^2 + psi0_z^2)

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then begin
           data(where(abs(tprime) lt 1e-6 or psi0 lt psis)) = 0.
       endif else begin
           data(where(abs(tprime) lt 1e-6 or psi0 gt psis)) = 0.
       endelse
  
       symbol = '!7n!N!X'
       d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; overlap
   ;===========================================
   endif else if(strcmp('overlap', name, /fold_case) eq 1) then begin

       xi = read_field('displacement', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, linear=linear, $
                       complex=complex, phi=phi0)

       psi0 = read_field('psi', x, y, t, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then jac = -1 else jac = 1

       ; d(xi)/dr = d(xi)/dpsi * |grad(psi)|
       ;          = (<xi,psi>/<psi,psi>)*sqrt(<psi,psi>)
       ;          = <xi,psi>/sqrt(<psi,psi>)
       data = abs(s_bracket(xi,psi0,x,y)/sqrt(s_bracket(psi0,psi0,x,y)))
  
       symbol = '!3|!6d!7n!D!8r!N!6/dr!3|!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; linearity
   ;===========================================
   endif else if(strcmp('linearity', name, /fold_case) eq 1) then begin

       xi = read_field('displacement', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange, linear=linear, $
                       complex=complex, phi=phi0)
       psi0 = read_field('psi', x, y, t, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)
       p0 = read_field('p', x, y, t, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=-1)


       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then jac = -1 else jac = 1

       ; d(xi)/dr = d(xi)/dpsi * |grad(psi)|
       ;          = (<xi,psi>/<psi,psi>)*sqrt(<psi,psi>)
       ;          = <xi,psi>/sqrt(<psi,psi>)
       l = p0/s_bracket(p0,psi0,x,y)*sqrt(s_bracket(psi0,psi0,x,y))
       data = xi/l
  
       symbol = '!3|!7n!D!8r!N!3|!6/!8L!Dp!N!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; curl(E)
   ;===========================================
   endif else if(strcmp('curl_E', name, /fold_case) eq 1) then begin

       if(not keyword_set(zvector) and not keyword_set(yvector)) then begin
           ephi_z = read_field('E_PHI', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, linear=linear, $
                               complex=complex, phi=phi0, op=3)
           if(icomplex eq 1) then begin 
               ez_phi = complex(0,ntor) * $
                 read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, $
                            rrange=xrange, zrange=yrange, linear=linear, $
                            complex=complex, phi=phi0)
           endif else begin
               ez_phi = read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                                   filename=filename, points=pts, $
                                   rrange=xrange,zrange=yrange,linear=linear, $
                                   complex=complex, phi=phi0, op=11)
           end
       end
       if(not keyword_set(rvector) and not keyword_set(yvector)) then begin
           ephi_r = read_field('E_PHI', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, linear=linear, $
                               complex=complex, phi=phi0, op=2)
           if(icomplex eq 1) then begin 
               er_phi = complex(0,ntor) * $
                 read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, $
                            rrange=xrange, zrange=yrange, linear=linear, $
                            complex=complex, phi=phi0)
           endif else begin
               er_phi = read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                                   filename=filename, points=pts, $
                                   rrange=xrange,zrange=yrange,linear=linear, $
                                   complex=complex, phi=phi0, op=11)           
           end
       end
       if(not keyword_set(rvector) and not keyword_set(zvector)) then begin
           ez_r = read_field('E_Z', x, y, t, slices=time, mesh=mesh, $
                             filename=filename, points=pts, $
                             rrange=xrange, zrange=yrange, linear=linear, $
                             complex=complex, phi=phi0, op=2)
           er_z = read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                             filename=filename, points=pts, $
                             rrange=xrange, zrange=yrange, linear=linear, $
                             complex=complex, phi=phi0,op=3)
       end
       
       symbol = '!9GX!6E!X'
       if(keyword_set(rvector)) then begin
           data = ephi_z - ez_phi
           symbol = symbol + '!9.G!8R!X'
       endif else if(keyword_set(zvector)) then begin
           data = er_phi - ephi_r
           symbol = symbol + '!9.G!8Z!X'
       endif else if(keyword_set(yvector)) then begin
           data = ez_r - er_z
           symbol = symbol + '!9.!6(!8R!9GP!6)!X'
       endif else begin
           data = sqrt((ephi_z - ez_phi)^2 + $
                       (ez_r - er_z)^2 + $
                       (er_phi - ephi_r)^2)
           symbol = '!3|!X'+symbol+'!3|!X'
       endelse
       
       d = dimensions(/potential, l0=-2, _EXTRA=extra)

   ;===========================================
   ; div(E)
   ;===========================================
   endif else if(strcmp('div_E', name, /fold_case) eq 1) then begin

       er_r = read_field('E_R', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange, zrange=yrange, linear=linear, $
                               complex=complex, phi=phi0, op=2)
       if(icomplex eq 1) then begin 
           ephi_phi = complex(0,ntor) * $
             read_field('E_PHI', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, linear=linear, $
                        complex=complex, phi=phi0)
       endif else begin
           ephi_phi = read_field('E_PHI', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, $
                               rrange=xrange,zrange=yrange,linear=linear, $
                               complex=complex, phi=phi0, op=11)           
       end
       ez_z = read_field('E_Z', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                         complex=complex, phi=phi0, op=3)
       
       data = er_r + ephi_phi + ez_z
       symbol = '!9G.!6E!X'
       d = dimensions(/potential, l0=-2, _EXTRA=extra)

   ;===========================================
   ; grad_psi
   ;===========================================
   endif else if(strcmp('grad_psi', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi', x, y, t, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time, op=2)
       psi_z = read_field('psi', x, y, t, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time, op=3)
       
       data = sqrt(psi_r^2 + psi_z^2)
 
       symbol = '!3|!9G!7w!3|!X'
       d = dimensions(_EXTRA=extra,b0=1,l0=1+itor)
       
   ;===========================================
   ; psi_norm
   ;===========================================
   endif else if(strcmp('psi_norm', name, /fold_case) eq 1) then begin

       psi = read_field('psi', x, y, t, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time)

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       
       data = (psi-flux0)/(psis-flux0)
 
       symbol = '!7W!D!8N!N!X'
       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; grad_psi_norm
   ;===========================================
   endif else if(strcmp('grad_psi_norm', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time, op=2)
       psi_z = read_field('psi', x, y, t, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, slice=time, op=3)

       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       
       data = sqrt(psi_r^2 + psi_z^2)/abs(psis-flux0)
 
       symbol = '!3|!9G!7W!3|!X'
       d = dimensions(l0=-1, _EXTRA=extra)

   ;===========================================
   ; ion temperature
   ;===========================================
;  endif else if(strcmp('ion temperature', name, /fold_case) eq 1) or $
;    (strcmp('ti', name, /fold_case) eq 1) then begin

;      P = read_field('P', x, y, t, slices=time, mesh=mesh, $
;                     filename=filename, points=pts, $
;                     rrange=xrange, zrange=yrange)

;      Pe = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
;                     filename=filename, points=pts, $
;                     rrange=xrange, zrange=yrange)

;      n = read_field('den', x, y, t, slices=time, mesh=mesh, $
;                     filename=filename, points=pts, $
;                     rrange=xrange, zrange=yrange)
;
;      data = (p-pe)/n
;      symbol = '!8T!Di!N!X'
;      d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; new ion temperature (as written by m3dc1 since 1/6/2011)
   ;===========================================
   endif else if(strcmp('ion temperature', name, /fold_case) eq 1) then begin

       ti = read_field('ti', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
  
       data = ti
       symbol = '!8T!Di!N!X'
       d = dimensions(/temperature, _EXTRA=extra)

   ;===========================================
   ; ion pressure
   ;===========================================
   endif else if(strcmp('ion pressure', name, /fold_case) eq 1) or $
     (strcmp('pi', name, /fold_case) eq 1) then begin

       P = read_field('P', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       Pe = read_field('Pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
  
       data = p-pe
       symbol = '!8p!Di!N!X'
       d = dimensions(/p0, _EXTRA=extra)

   ;===========================================
   ; angular momentum
   ;===========================================
   endif else if(strcmp('angular momentum', name, /fold_case) eq 1) or $
     (strcmp('lz', name, /fold_case) eq 1) then begin

       V = read_field('V', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
  
       data = n*v
       symbol = '!8L!D!9P!N!X'
       d = dimensions(/n0, /v0, /l0, _EXTRA=extra)

   ;===========================================
   ; toroidal current
   ;===========================================
   endif else if(strcmp('jy', name, /fold_case) eq 1) then begin

       lp = read_field('psi', x, y, t, slices=time, mesh=mesh, op=7, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                      complex=complex,phi=phi0, wall_mask=wall_mask)

       if(itor eq 1) then begin
          psir = read_field('psi', x, y, t, slices=time, mesh=mesh, op=2, $
                            filename=filename, points=pts, mask=mask, $
                            rrange=xrange, zrange=yrange, linear=linear, $
                            complex=complex,phi=phi0, wall_mask=wall_mask)

           r = radius_matrix(x,y,t)
           data = -(lp - psir/r)/r
       endif else data = -lp

       symbol = '!8J!D!9P!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; toroidal current
   ;===========================================
   endif else if(strcmp('jy_ext', name, /fold_case) eq 1) then begin

       lp = read_field('psi_ext', x, y, t, slices=time, mesh=mesh, op=7, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                      complex=complex,phi=phi0, wall_mask=wall_mask)

       if(itor eq 1) then begin
          psir = read_field('psi_ext', x, y, t, slices=time, mesh=mesh, op=2, $
                            filename=filename, points=pts, mask=mask, $
                            rrange=xrange, zrange=yrange, linear=linear, $
                            complex=complex,phi=phi0, wall_mask=wall_mask)

           r = radius_matrix(x,y,t)
           data = -(lp - psir/r)/r
       endif else data = -lp

       symbol = '!8J!D!9P!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; toroidal current
   ;===========================================
   endif else if(strcmp('jy_plasma', name, /fold_case) eq 1) then begin

       lp = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=7, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                      complex=complex,phi=phi0,wall_mask=wall_mask,outval=mask_val)
       psir = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=2, $
                         filename=filename, points=pts, mask=mask, $
                         rrange=xrange, zrange=yrange, linear=linear, $
                        complex=complex,phi=phi0,wall_mask=wall_mask,outval=mask_val)


       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
           data = -(lp - psir/r)/r
       endif else data = -lp

       symbol = '!8J!D!9P!N!X'
       d = dimensions(/j0,_EXTRA=extra)


   ;===========================================
   ; minor radius
   ;===========================================
   endif else if(strcmp('minor radius', name, /fold_case) eq 1) or $
     (strcmp('r', name) eq 1) then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       nulls, psi, x, y, xpoint=xpoint, axis=axis, $
         filename=filename, _EXTRA=extra
      
       x0 = axis[0]
       z0 = axis[1]

       xx = fltarr(n_elements(t),n_elements(x),n_elements(y))
       zz = fltarr(n_elements(t),n_elements(x),n_elements(y))
       for k=0, n_elements(t)-1 do begin
           for i=0, n_elements(y)-1 do xx[k,*,i] = x
           for i=0, n_elements(x)-1 do zz[k,i,*] = y
       end
       data = sqrt((xx-x0)^2 + (zz-z0)^2)

       symbol = '!8r!X'
       d = dimensions(/l0, _EXTRA=extra)

   ;===========================================
   ; qstar (local cylindrical approximation to q)
   ;===========================================
   endif else if(strcmp('q_cyl', name, /fold_case) eq 1) then begin

       rho = read_field('r', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       bp = read_field('bp', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       bt = read_field('bt', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
          r = radius_matrix(x,y,t)
       endif else begin
          r = rzero
       end
       
       data = rho*bt / (r*bp)

       symbol = '!8q!D!3*!N!X'
       d = dimensions(_EXTRA=extra)


   ;===========================================
   ; polodal angle
   ;===========================================
   endif else if(strcmp('polodal angle', name, /fold_case) eq 1) or $
     (strcmp('theta', name) eq 1) then begin


       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       nulls, psi, x, y, xpoint=xpoint, axis=axis, $
         filename=filename, _EXTRA=extra
      
       x0 = axis[0]
       z0 = axis[1]

       xx = fltarr(n_elements(t),n_elements(x),n_elements(y))
       zz = fltarr(n_elements(t),n_elements(x),n_elements(y))
       for k=0, n_elements(t)-1 do begin
           for i=0, n_elements(y)-1 do xx[k,*,i] = x
           for i=0, n_elements(x)-1 do zz[k,i,*] = y
       end
       data = atan(zz-z0,xx-x0)

       symbol = '!7h!X'

   ;===========================================
   ; pest angle
   ;===========================================
   endif else if(strcmp('pest angle', name, /fold_case) eq 1) then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange,_EXTRA=extra)

       forward_function find_lcfs
       psilim = find_lcfs(psi, x, y,axis=axis, xpoint=xpoint, flux0=flux0, $
                          _EXTRA=extra, filename=filename)

       rr = radius_matrix(x,y,t)
       zz = z_matrix(x,y,t)

       forward_function flux_coord_field
       rrfc = flux_coord_field(rr, psi, x, y, t, slice=time, $
                                 flux=flux, angle=angle, /pest, $
                                 filename=filename, points=pts, $
                                 _EXTRA=extra)
       zzfc = flux_coord_field(zz, psi, x, y, t, slice=time, $
                                 flux=flux, angle=angle, /pest, $
                                 filename=filename, points=pts, $
                                 _EXTRA=extra)

       psinorm = (psi - flux0) / (psilim - flux0)
       data = psi*0.
       for i=0, n_elements(data)-1 do begin
           if(psinorm[i] gt 1.) then continue

           dist = (rrfc-rr[i])^2 + (zzfc-zz[i])^2
           dum = min(dist, j)
           n_guess = j/n_elements(angle)
           m_guess = j - n_guess*n_elements(angle)
           data[i] = angle(n_guess)
       end

       symbol = '!7h!D!6PEST!N!X'
       d = dimensions()

   ;===========================================
   ; Field strength
   ;===========================================
   endif else if( (strcmp('field strength', name, /fold_case) eq 1) $
                  or (strcmp('b', name, /fold_case) eq 1)) $
     then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, complex=complex)
       if(extsubtract eq 1 and version lt 8) then begin
           psi = psi + read_field('psi_ext', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, complex=complex)
       end

       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex)
       if(extsubtract eq 1 and version lt 8) then begin
           I = I + read_field('I_ext', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex)
       end

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(icomplex eq 1 and time ne -1) then begin
           b2 = (s_bracket(psi,conj(psi),x,y) + I*conj(I))/r^2
           fp = read_field('fp', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex)
           if(extsubtract eq 1 and version lt 8) then begin
               fp = fp + $
                 read_field('fp_ext', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex)
           end
           b2 = b2 + s_bracket(fp,conj(fp),x,y) $
             - a_bracket(fp, conj(psi),x,y)/r $
             - a_bracket(conj(fp), psi,x,y)/r
           b2 = real_part(b2)
           b2 = b2 / 2. ; this comes from the cos^2 dependence of the field
       endif else begin 
           b2 = s_bracket(psi,psi,x,y)/r^2 + I^2/r^2
           if(i3d eq 1) then begin
               fp = read_field('fp', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange)
               if(extsubtract eq 1 and version lt 8) then begin
                   fp = fp + $
                     read_field('fp_ext', x, y, t, slices=time, mesh=mesh, $
                                filename=filename, points=pts, linear=linear, $
                                rrange=xrange, zrange=yrange)
               end

               b2 = b2 + s_bracket(fp,fp,x,y) $
                 - 2.*a_bracket(fp,psi,x,y)/r
           end           
       endelse

       data = sqrt(b2)
       symbol = '!3|!5B!3|!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; Field energy
   ;===========================================
   endif else if( (strcmp('field energy', name, /fold_case) eq 1) $
                  or (strcmp('b2', name, /fold_case) eq 1)) $
     then begin

       psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                       linfac=linfac,op=2)
       psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                       linfac=linfac,op=3)
       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, $
                     linfac=linfac)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       b2 = (psi_r*conj(psi_r) + psi_z*conj(psi_z) + I*conj(I))/r^2
       if(icomplex eq 1 and time ne -1) then begin
           ; if the fields are ~exp(i n phi), then
           ; this is the toroidally-averaged value of |B| !

           f_r = read_field('f_r', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                         linfac=linfac, op=2)
           f_z = read_field('f_z', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                         linfac=linfac, op=3)

           fpr = complex(0., ntor)*f_r
           fpz = complex(0., ntor)*f_z
           b2 = b2 + fpr*conj(fpr) + fpz*conj(fpz) $
             + (fpr*conj(psi_z) - fpz*conj(psi_r))/r $
             + (conj(fpr)*psi_z - conj(fpz)*psi_r)/r

           b2 = b2 / 2. ; this comes from the cos^2 dependence of the field
           b2 = real_part(b2)
       endif

       nolf = 1
       data = b2/2.  ; W = |B|^2 / 2
       symbol = '!3|!5B!3|!U!62!N!3/!62!X'
       d = dimensions(p0=1, _EXTRA=extra)

   ;===========================================
   ; Poloidal Field strength
   ;===========================================
   endif else if( (strcmp('poloidal field strength', name, /fold_case) eq 1) $
                  or (strcmp('bp', name, /fold_case) eq 1)) $
     then begin

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        complex=complex, phi=phi0, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = s_bracket(psi,psi,x,y)/r^2

       if((i3d eq 1 or icomplex eq 1) and $
          not (time eq -1 and ilin eq 1)) then begin
          fp = read_field('fp', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          complex=complex, phi=phi0, linear=linear, $
                          rrange=xrange, zrange=yrange)

           data = data - 2.*a_bracket(fp,psi,x,y)/r + s_bracket(fp,fp,x,y)
       end

       data = sqrt(data)
       symbol = '!3|!5B!D!8p!N!3|!X'
       d = dimensions(/b0, _EXTRA=extra)


   ;===========================================
   ; Toroidal Field strength
   ;===========================================
   endif else if( (strcmp('toroidal field strength', name, /fold_case) eq 1) $
                  or (strcmp('bt', name, /fold_case) eq 1)) $
     then begin

       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, logical=logical)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = sqrt(I^2/r^2)
       symbol = '!3|!8B!D!9P!N!3|!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; toroidal vector potential
   ;===========================================
   endif else if(strcmp('toroidal vector potential', name, /fold_case) eq 1) or $
     (strcmp('ay', name, /fold_case) eq 1) then begin

       psi = read_field('psi',x,y,t,slices=time, mesh=mesh, filename=filename, $
                        points=pts,rrange=xrange,zrange=yrange, $
                        linear=linear,complex=complex,phi=phi0)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
   
       data = psi/r
       symbol = '!8A!D!9P!N!X'
       d = dimensions(/b0, /l0, _EXTRA=extra)


   ;===========================================
   ; Pitch Angle
   ;===========================================
   endif else if( (strcmp('pitch_angle_z', name, /fold_case) eq 1)) $
     then begin

      if(keyword_set(linear)) then begin
         by0 = read_field('by', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bz0 = read_field('bz', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         by1 = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         bz1 = read_field('bz', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         data = bz1/by0 - by1*bz0/by0
         symbol = '!8B!DZ!N!3/!8B!D!9P!N!x'
      endif else begin
         by = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bz = read_field('bz', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         data = bz/by
         symbol = '!8B!DZ!N!3/!8B!D!9P!N!x'
      endelse 

       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; Pitch Angle
   ;===========================================
   endif else if( (strcmp('pitch_angle_x', name, /fold_case) eq 1)) $
     then begin

      if(keyword_set(linear)) then begin
         by0 = read_field('by', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bx0 = read_field('bx', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         by1 = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         bx1 = read_field('bx', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         data = bx1/by0 - by1*bx0/by0
         symbol = '!8B!DR!N!3/!8B!D!9P!N!x'
      endif else begin
         by = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bx = read_field('bx', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         data = bx/by
         symbol = '!8B!DR!N!3/!8B!D!9P!N!x'
      endelse 

       d = dimensions(_EXTRA=extra)


   ;===========================================
   ; fractional Pitch Angle change
   ;===========================================
   endif else if( (strcmp('dpitch_angle_z', name, /fold_case) eq 1)) $
     then begin

      if(keyword_set(linear)) then begin
         by0 = read_field('by', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bz0 = read_field('bz', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         by1 = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         bz1 = read_field('bz', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         data = bz1/bz0 - by1/by0
         symbol = '!7da!DZ!N!3/!7a!DZ!N!x'
      endif else begin
         by = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bz = read_field('bz', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         data = bz/by
         symbol = '!8B!DZ!N!3/!8B!D!9P!N!x'
      endelse 

       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; fractional pitch angle change
   ;===========================================
   endif else if( (strcmp('dpitch_angle_x', name, /fold_case) eq 1)) $
     then begin

      if(keyword_set(linear)) then begin
         by0 = read_field('by', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bx0 = read_field('bx', x, y, t, slices=-1, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         by1 = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         bx1 = read_field('bx', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, linear=linear)
         data = bx1/bx0 - by1/by0
         symbol = '!7da!DR!N!3/!7a!DR!N!x'
      endif else begin
         by = read_field('by', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         bx = read_field('bx', x, y, t, slices=time, mesh=mesh, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange)
         data = bx/by
         symbol = '!8B!DR!N!3/!8B!D!9P!N!x'
      endelse 

       d = dimensions(_EXTRA=extra)

   ;===========================================
   ; Kinetic energy density
   ;===========================================
   endif else if( (strcmp('ke', name, /fold_case) eq 1)) then begin

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       v = read_field('v', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = 0.5*n*(s_bracket(u,u,x,y)/r^2+v^2/r^2 + $
                         s_bracket(chi,chi,x,y)+2.*a_bracket(chi,u,x,y)/r)
       endif else if(ivform eq 1) then begin
           data = 0.5*n*(r^2*s_bracket(u,u,x,y)/r^2+r^2*v^2 + $
                         s_bracket(chi,chi,x,y)/r^4+2.*a_bracket(chi,u,x,y)/r)
       endif
       symbol ='!6Kinetic Energy Density!X'
       d = dimensions(/p0, _EXTRA=extra)

   ;===========================================
   ; Alfven velocity
   ;===========================================
   endif else if( (strcmp('va', name, /fold_case) eq 1)) $
     then begin

       b = read_field('b', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       data = b/sqrt(den)
       symbol = '!8v!DA!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; (minor) radial current density
   ;===========================================
   endif else if(strcmp('jn', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange, op=2)
       psi0_z = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange, op=3)
       i_r = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, op=2)
       i_z = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, op=3)


       if(itor eq 1) then begin
           if(n_elements(at_points) eq 0) then begin
               r = radius_matrix(x,y,t)
           endif else begin
               r = at_points[0,*]
           endelse
       endif else r = 1.

       psipsi = sqrt(psi0_r^2 + psi0_z^2)
       
       data = (i_z*psi0_r - i_r*psi0_z)/(r*psipsi)

       if(ntor ne 0) then begin

           psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=2)
           psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=3)
           fp_r = read_field('fp', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=2)
           fp_z = read_field('fp', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, op=3)

           data = data $
             - complex(0,ntor)*(psi_r*psi0_r + psi_z*psi0_z)/(r^2*psipsi) $
             + complex(0,ntor)*(fp_z*psi0_r - fp_r*psi0_z)/(r*psipsi)
       endif

       symbol = '!8J!Dr!N!X'
       d = dimensions(/j0, _EXTRA=extra)

   ;===========================================
   ; (major) radial current density
   ;===========================================
   endif else if(strcmp('jx', name, /fold_case) eq 1) then begin
       
       i_z = read_field('i', x, y, t, slices=time, mesh=mesh, op=3, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, $
                       phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = -i_z / r

       if(i3d eq 1) then begin
           psi_rp = read_field('psi', x, y, t, slices=time, mesh=mesh, op=12, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           f_zpp = read_field('fp', x, y, t, slices=time, mesh=mesh, op=13, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           data = data - f_zpp / r + psi_rp/r^2
       endif else if(ntor ne 0) then begin
           psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, op=2, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                              phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)           
           f_pz = read_field('fp', x, y, t, slices=time, mesh=mesh, op=3, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, complex=complex, $
                             phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)
           data = data + complex(0.,ntor) * (-f_pz / r + psi_r/r^2)
       endif
       
       symbol = '!8J!DR!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; (major) radial current density (plasma)
   ;===========================================
   endif else if(strcmp('jx_plasma', name, /fold_case) eq 1) then begin
       
       i_z = read_field('I_plasma', x, y, t, slices=time, mesh=mesh, op=3, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, wall_mask=wall_mask)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = -i_z / r

       if(i3d eq 1) then begin
           psi_rp = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=12, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                               wall_mask=wall_mask)

           f_zpp = read_field('f_plasma', x, y, t, slices=time, mesh=mesh, op=23, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                             wall_mask=wall_mask)

           data = data - f_zpp / r + psi_rp/r^2
       endif else if(ntor ne 0) then begin
           psi_r = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=2, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                             wall_mask=wall_mask)

           f_z = read_field('f_plasma', x, y, t, slices=time, mesh=mesh, op=3, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                            wall_mask=wall_mask)

           data = data + ntor^2 * f_z / r + complex(0., ntor)*psi_r/r^2
       endif
       
       symbol = '!8J!DR!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; (major) radial current density (external)
   ;===========================================
   endif else if(strcmp('jx_ext', name, /fold_case) eq 1) then begin
       
       i_z = read_field('i_ext', x, y, t, slices=time, mesh=mesh, op=3, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, $
                       phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = -i_z / r

       if(i3d eq 1) then begin
           psi_rp = read_field('psi_ext', x, y, t, slices=time, mesh=mesh, op=12, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           f_zpp = read_field('fp_ext', x, y, t, slices=time, mesh=mesh, op=13, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           data = data - f_zpp / r + psi_rp/r^2
       endif else if(ntor ne 0) then begin
           psi_r = read_field('psi_ext', x, y, t, slices=time, mesh=mesh, op=2, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                              phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)           
           f_pz = read_field('fp_ext', x, y, t, slices=time, mesh=mesh, op=3, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, complex=complex, $
                             phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)
           data = data + complex(0.,ntor) * (-f_pz / r + psi_r/r^2)
       endif
       
       symbol = '!8J!DR!N!X'
       d = dimensions(/j0,_EXTRA=extra)
       
   ;===========================================
   ; vertical current density
   ;===========================================
   endif else if(strcmp('jz', name, /fold_case) eq 1) then begin
       
       i_r = read_field('i', x, y, t, slices=time, mesh=mesh, op=2, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, $
                      phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = i_r / r

       if(i3d eq 1) then begin
           psi_zp = read_field('psi', x, y, t, slices=time, mesh=mesh, op=13, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)

           f_rpp = read_field('fp', x, y, t, slices=time, mesh=mesh, op=12, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           data = data + f_rpp / r + psi_zp/r^2
        endif else if(ntor ne 0) then begin
           psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, op=3, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           f_rp = read_field('fp', x, y, t, slices=time, mesh=mesh, op=2, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, complex=complex, $
                             phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)
           data = data + complex(0.,ntor) * (f_rp / r + psi_z/r^2)
        end
       
       symbol = '!8J!DZ!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; vertical current density (plasma)
   ;===========================================
   endif else if(strcmp('jz_plasma', name, /fold_case) eq 1) then begin
       
       i_r = read_field('I_plasma', x, y, t, slices=time, mesh=mesh, op=2, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, $
                        wall_mask=wall_mask)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = i_r / r

       if(i3d eq 1) then begin
           psi_zp = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=13, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                              wall_mask=wall_mask)

           f_rpp = read_field('f_plasma', x, y, t, slices=time, mesh=mesh, op=22, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                             wall_mask=wall_mask)

           data = data + f_rpp / r + psi_zp/r^2
        endif else if(ntor ne 0) then begin
           psi_z = read_field('psi_plasma', x, y, t, slices=time, mesh=mesh, op=3, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                             wall_mask=wall_mask)

           f_r = read_field('f_plasma', x, y, t, slices=time, mesh=mesh, op=2, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                           wall_mask=wall_mask)

           data = data - ntor^2 * f_r / r + complex(0., ntor)*psi_z/r^2
       endif
       
       symbol = '!8J!DZ!N!X'
       d = dimensions(/j0,_EXTRA=extra)

   ;===========================================
   ; vertical current density (external)
   ;===========================================
   endif else if(strcmp('jz_ext', name, /fold_case) eq 1) then begin
       
       i_r = read_field('i_ext', x, y, t, slices=time, mesh=mesh, op=2, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, complex=complex, $
                      phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = i_r / r

       if(i3d eq 1) then begin
           psi_zp = read_field('psi_ext', x, y, t, slices=time, mesh=mesh, op=13, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)

           f_rpp = read_field('fp_ext', x, y, t, slices=time, mesh=mesh, op=12, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           data = data + f_rpp / r + psi_zp/r^2
        endif else if(ntor ne 0) then begin
           psi_z = read_field('psi_ext', x, y, t, slices=time, mesh=mesh, op=3, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            phi=phi0,  wall_mask=wall_mask,outval=mask_val,mask=mask)

           f_rp = read_field('fp_ext', x, y, t, slices=time, mesh=mesh, op=2, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, complex=complex, $
                             phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)
           data = data + complex(0.,ntor) * (f_rp / r + psi_z/r^2)
        end
       
       symbol = '!8J!DZ!N!X'
       d = dimensions(/j0,_EXTRA=extra)
       
   ;===========================================
   ; poloidal current density
   ;===========================================
   endif else if(strcmp('jp', name, /fold_case) eq 1) then begin
       
       psi = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, /equilibrium, $
                        wall_mask=wall_mask)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear,  $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      wall_mask=wall_mask)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       data = s_bracket(i,psi,x,y)/(r*sqrt(s_bracket(psi,psi,x,y)))
       symbol = '!8J!Dp!N!X'
       d = dimensions(/j0,_EXTRA=extra)


   ;===========================================
   ; eta * Jy^2
   ;===========================================
   endif else if(strcmp('etaJ2', name, /fold_case) eq 1) then begin

       jy = read_field('jy', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, $
                       complex=complex, phi=phi0, wall_mask=wall_mask, $
                      logical=logical)
       eta = read_field('eta', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, $
                        complex=complex, phi=phi0, wall_mask=wall_mask, $
                       logical=logical)

       data = eta*jy^2
       symbol = field_data('!7g!8J!6!U2!N', units=d, itor=itor)
       d = dimensions(/p0,t0=-1,_EXTRA=extra)
       
   ;===========================================
   ; del*(psi)
   ;===========================================
   endif else if(strcmp('jphi', name, /fold_case) eq 1) then begin

       psi_lp = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange, op=7, $
                          complex=complex, phi=phi0, wall_mask=wall_mask)
       data = psi_lp
       if(itor eq 1) then begin
           psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                              filename=filename, points=pts, linear=linear, $
                              rrange=xrange, zrange=yrange, op=2, $
                             complex=complex, phi=phi0, wall_mask=wall_mask)
           r = radius_matrix(x,y,t)
           data = data - psi_r / r
       end

       symbol = field_data('jphi', units=d, itor=itor)
       d = dimensions(/j0,l0=itor,_EXTRA=extra)

   ;===========================================
   ; vorticity
   ;===========================================
      endif else if(strcmp('vor', name, /fold_case) eq 1) then begin

          phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, $
                           complex=complex,phi=phi0)

          data = delstar(phi,x,y,tor=itor)
          symbol = field_data('vor', units=d, itor=itor)

   ;===========================================
   ; divergence
   ;===========================================
     endif else if(strcmp('com', name, /fold_case) eq 1) then begin

         chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, complex=complex, $
                          phi=phi0)

         data = laplacian(chi,x,y,tor=itor)
         symbol = field_data('com', units=d, itor=itor)

   ;===========================================
   ; R^2 vorticity
   ;===========================================
      endif else if(strcmp('r2vor', name, /fold_case) eq 1) then begin

          phi_lp = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, $
                           complex=complex,phi=phi0,op=7)
          phi_r = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, $
                           complex=complex,phi=phi0,op=2)
          r = radius_matrix(x,y,t)

          data = r^2*phi_lp + 2.*r*phi_r
          symbol = field_data('vor', units=d, itor=itor)

   ;===========================================
   ; helicity
   ;===========================================
      endif else if(strcmp('helicity', name, /fold_case) eq 1) then begin

          psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, phi=phi0, $
                            complex=complex, linear=linear)
          i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, phi=phi0, $
                            complex=complex, linear=linear)
          f = read_field('f', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, mask=mask, $
                           rrange=xrange, zrange=yrange, phi=phi0, $
                            complex=complex, linear=linear)

          r = radius_matrix(x,y,t)

          data = (psi*conj(i) + conj(psi)*i)/r^4 $
            - s_bracket(f, conj(psi), x, y)/r^2 $
            - s_bracket(conj(f), psi, x, y)/r^2
          d = dimensions(l0=4, b0=2)
          symbol = field_data('helicity', units=d, itor=itor)


   ;===========================================
   ; chi_perp
   ;===========================================
     endif else if(strcmp('chi_perp', name, /fold_case) eq 1) then begin

         kappa = read_field('kappa', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, phi=phi0, $
                           linear=linear, logical=logical)
         den = read_field('den', x, y, t, slices=time, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, /equilibrium,$
                         logical=logical)

         data = kappa/den
         d = dimensions(l0=2, t0=-1, _EXTRA=extra)
         symbol = '!7v!X'

     ;===========================================
     ; chi_parallel
     ;===========================================
     endif else if(strcmp('chi_par', name, /fold_case) eq 1) then begin

         kappar = read_parameter('kappar', filename=filename)
         den = read_field('den', x, y, t, slices=time, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, /equilibrium)

         data = kappar/den
         d = dimensions(l0=2, t0=-1, _EXTRA=extra)
         symbol = '!7v!D!3||!N!X'

   ;===========================================
   ; mu_perp
   ;===========================================
     endif else if(strcmp('mu_perp', name, /fold_case) eq 1) then begin

         visc = read_field('visc', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, phi=phi0, $
                          linear=linear)
         den = read_field('den', x, y, t, slices=time, $
                          filename=filename, points=pts, mask=mask, $
                          rrange=xrange, zrange=yrange, /equilibrium)

         data = visc/den
         d = dimensions(l0=2, t0=-1, _EXTRA=extra)
         symbol = '!7l!D!9x!N!X'

   ;===========================================
   ; rotational transform
   ;===========================================
   endif else if(strcmp('iota', name, /fold_case) eq 1) then begin

       minor_r = read_field('r', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, $
                            rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       I = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       bt = sqrt(I^2/r^2)
       bp = sqrt(s_bracket(psi,psi,x,y)/r^2)

       data = 2.*!pi*(r * bp) / minor_r * bt
       symbol = '!8i!X'

   ;===========================================
   ; angular velocity
   ;===========================================
   endif else if(strcmp('omega', name, /fold_case) eq 1) then begin

       v = read_field('v', x, y, t, slices=time, mesh=mesh, complex=complex, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, op=op, phi=phi0, $
                      logical=logical)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           if(n_elements(op) ne 0) then begin
               print, 'Warning: using op on omega'
           end
           data = v/r^2
       endif else if(ivform eq 1) then begin
           data = v
       endif
       if(itor eq 0) then begin
          data = 2.*!pi*data / period
       end

       symbol = '!7X!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; pprime
   ;===========================================
   endif else if(strcmp('pprime', name, /fold_case) eq 1) then begin

       p = read_field('p', x, y, t, slices=time, mesh=mesh, complex=complex, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, op=op, phi=phi0)
       psi = read_field('psi', x, y, t, slices=time, /equilibrium, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange, op=op, phi=phi0)

       data = s_bracket(p,psi,x,y)/s_bracket(psi,psi,x,y)

       symbol = "!8p'!X"
       d = dimensions(p0=1, b0=1, l0=1+itor,_EXTRA=extra)


   ;===========================================
   ; electron angular velocity
   ;===========================================
   endif else if(strcmp('omega_e', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename)
       print, 'db = ', db

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       jy = read_field('jy', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                       filename=filename, points=pts, $
                       rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = omega - (db*jy/den)/r
       endif else if(ivform eq 1) then begin
           data = omega - (db*jy/den)
       endif
       symbol = '!7x!D!8e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;==========================================================
   ; v_omega (the omega in v = r^2 omega grad(phi) + (K/n) B
   ;==========================================================
   endif else if(strcmp('v_omega', name, /fold_case) eq 1) then begin

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       data = omega - i/(r^2*s_bracket(psi,psi,x,y)) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r)

       symbol = '!7x!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;==========================================================
   ; v_K (the K in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('v_K', name, /fold_case) eq 1) then begin

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = den/s_bracket(psi,psi,x,y) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r)

       symbol = '!8K!X'
       d = dimensions(/v0, /n0, b0=-1, _EXTRA=extra)

   ;==========================================================
   ; v_K (the K in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('v_K_n', name, /fold_case) eq 1) then begin

       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = 1./s_bracket(psi,psi,x,y) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r)

       symbol = '!8K!6/!8n!X'
       d = dimensions(/v0, b0=-1, _EXTRA=extra)

   ;==========================================================
   ; ve_omega (the omega in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('ve_omega', name, /fold_case) eq 1) then begin

       di = read_parameter('db', filename=filename)
       print, 'di = ', di
       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi_lp = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, op=7)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       data = omega - i/(r^2*s_bracket(psi,psi,x,y)) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r $
          - (di/den)*s_bracket(i,psi,x,y)) $
         + (di/den)*(psi_lp - itor*dx(psi,x)/r)/r^2

       symbol = '!7x!D!8e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;==========================================================
   ; ve_K (the K in v = r^2 omega grad(phi) + (K/n) B)
   ;==========================================================
   endif else if(strcmp('ve_K', name, /fold_case) eq 1) then begin

       di = read_parameter('db', filename=filename)
       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = 1./s_bracket(psi,psi,x,y) * $
         (den*(r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r) $
          - di*s_bracket(i,psi,x,y))

       symbol = '!8K!D!8e!N!X'
       d = dimensions(/v0, /n0, b0=-1, _EXTRA=extra)


   ;===========================================
   ; electron angular velocity
   ;===========================================
   endif else if(strcmp('omega_perp_e', name, /fold_case) eq 1) then begin

       di = read_parameter('db', filename=filename)
       print, 'di = ', di
       omega = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi_lp = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, op=7)
       i = read_field('I', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
       
       omega_e = omega - i/(r^2*s_bracket(psi,psi,x,y)) * $
         (r^2*s_bracket(u,psi,x,y) + a_bracket(chi,psi,x,y)/r $
          - (di/den)*s_bracket(i,psi,x,y)) $
         + (di/den)*(psi_lp - itor*dx(psi,x)/r)/r^2

       Bp = sqrt(s_bracket(psi,psi,x,y))/r
       B = sqrt(s_bracket(psi,psi,x,y) + i^2)/r

       data = (Bp/B)*omega_e

;       symbol = '!7x!S!D!9x!N!S!U!8e!N!X'
       symbol = '!7X!D!8e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; cyclotron frequency
   ;===========================================
   endif else if(strcmp('omega_ci', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', db
       if(db eq 0.) then begin
           print, 'Warning: Assuming d_i = 1.'
           db = 1.
       endif

       B = read_field('B', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = B/db
       symbol = '!7x!8!Dc!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; diamagnetic frequency
   ;===========================================
   endif else if(strcmp('omega_*', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', filename, db

       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)


       data = db*s_bracket(p,psi,x,y)/s_bracket(psi,psi,x,y) / den

       symbol = '!7x!6!D*!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; diamagnetic frequency
   ;===========================================
   endif else if(strcmp('omega_*i', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', filename, db

       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       pe = read_field('pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       data = db*s_bracket(p-pe,psi,x,y) / s_bracket(psi,psi,x,y) / den

       symbol = '!7x!6!D*i!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; diamagnetic frequency
   ;===========================================
   endif else if(strcmp('omega_*e', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)
       print, 'db = ', filename, db

       pe = read_field('pe', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       den = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)

       data = -db*s_bracket(pe,psi,x,y)/s_bracket(psi,psi,x,y) / den

       symbol = '!7x!6!D*e!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   ;===========================================
   ; ExB frequency
   ;===========================================
   endif else if(strcmp('omega_ExB', name, /fold_case) eq 1) then begin

       omega = read_field('v_omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       w_star_i =  read_field('omega_*i', x, y, t, slices=time, mesh=mesh, $
                              filename=filename, points=pts, $
                              rrange=xrange, zrange=yrange)
       
       data = omega - w_star_i

       symbol = '!7x!6!DE!9X!6B!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)


   ;===========================================
   ; Larmor radius
   ;===========================================
   endif else if(strcmp('rho_i', name, /fold_case) eq 1) then begin

       db = read_parameter('db', filename=filename, _EXTRA=extra)

       Ti = read_field('Ti', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)
       B = read_field('B', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange)


       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = db*sqrt(2.*Ti)/B
       symbol = '!7q!8!Di!N!X'
       d = dimensions(l0=1, _EXTRA=extra)


   ;===========================================
   ; parallel thermal gradient
   ;===========================================
   endif else if(strcmp('xbdotgradt', name, /fold_case) eq 1) then begin

       Te0 = read_field('Te', x, y, t, slices=time, $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange)
       psi0 = read_field('psi', x, y, t, slices=time, $
                         filename=filename, points=pts, /equilibrium, $
                         rrange=xrange, zrange=yrange)
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ilin eq 1) then begin
           Te1 = read_field('Te', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange)
           psi1 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, linear=linear, $
                           rrange=xrange, zrange=yrange)
       
           data = a_bracket(Te1, psi0, x, y)/r + a_bracket(Te0, psi1, x, y)/r

           if(ntor ne 0) then begin
               i1 = read_field('i', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange)
               i0 = read_field('i', x, y, t, slices=time, $
                               filename=filename, points=pts, /equilibrium, $
                               rrange=xrange, zrange=yrange)
               f1 = read_field('f', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange)

               data = data + complex(0.,ntor)* $
                 (Te0*i1 + Te1*i0 - s_bracket(Te0, f1, x, y))
           end
       end else begin
           data = a_bracket(Te0, psi0, x, y)/r + a_bracket(Te0, psi0, x, y)/r           
       end

;       data = te0

       symbol = '!8B!9.G!8T!X'
       d = dimensions(l0=1, _EXTRA=extra)

   ;===========================================
   ; parallel pressure
   ;===========================================
   endif else if(strcmp('xbdotgradp', name, /fold_case) eq 1) then begin

       if(ilin eq 1) then begin
           p1 = read_field('p', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=complex)
           p0 = read_field('p', x, y, t, slices=time, $
                            filename=filename, points=pts, /equilibrium, $
                            rrange=xrange, zrange=yrange, complex=complex)
           psi1 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           filename=filename, points=pts, linear=linear, $
                           rrange=xrange, zrange=yrange, complex=complex)
           psi0 = read_field('psi', x, y, t, slices=time, $
                           filename=filename, points=pts, /equilibrium, $
                           rrange=xrange, zrange=yrange, complex=complex)
           if(itor eq 1) then begin
               r = radius_matrix(x,y,t)
           endif else r = 1.
       
           data = a_bracket(p1, psi0, x, y)/r + a_bracket(p0, psi1, x, y)/r

           if(ntor ne 0) then begin
               i0 = read_field('i', x, y, t, slices=time, $
                               filename=filename, points=pts, /equilibrium, $
                               rrange=xrange, zrange=yrange, complex=complex)
               f1 = read_field('f', x, y, t, slices=time, mesh=mesh, $
                               filename=filename, points=pts, linear=linear, $
                               rrange=xrange, zrange=yrange, complex=complex)

               data = data + complex(0.,ntor)* $
                 (p1*i0/r^2 - s_bracket(p0, f1, x, y))
           end
       end

;       data = te0

       symbol = '!8B!9.G!8p!X'
       d = dimensions(p0=1, b0=1, l0=-1, _EXTRA=extra)


   ;===========================================
   ; parallel flow
   ;===========================================
   endif else if(strcmp('vpar', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       w = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       I = read_field('I', x, y, t, slices=time, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)
         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       
       b2 = (s_bracket(psi,psi,x,y) + I^2)/r^2
       if(ivform eq 0) then begin
           data = ((s_bracket(phi,psi,x,y) + w*I)/r^2 $
                   + a_bracket(chi,psi,x,y)/r)/sqrt(b2)
       endif else begin
           data = (s_bracket(phi,psi,x,y) + w*I $
                   + a_bracket(chi,psi,x,y)/r^3)/sqrt(b2)
       endelse
       symbol = '!8u!D!9#!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   ;===========================================
   ; poloidal flow
   ;===========================================
   endif else if(strcmp('vpol', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       
       if(ivform eq 0) then begin
           data = 0.
       endif else begin
           data = r*(s_bracket(phi,psi,x,y) $
                   + a_bracket(chi,psi,x,y)/r^3) $
             / sqrt(s_bracket(psi,psi,x,y))
       endelse
       symbol = '!8u!D!7h!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('omtest', name, /fold_case) eq 1) then begin

       v_omega = read_field('v_omega', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       v_k_n = read_field('v_k_n', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)

         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = v_omega + i*v_k_n/r^2
       
       symbol = '!7X!D!6test!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)


   ;===========================================
   ; perpendicular flow v.(B x grad(psi))
   ;===========================================
   endif else if(strcmp('vperp', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       w = read_field('omega', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, linear=linear, $
                      rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, linear=linear, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       I = read_field('I', x, y, t, slices=time, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)
         
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       
       pp = s_bracket(psi,psi,x,y)
       b2 = (pp + I^2)/r^2
       if(ivform eq 0) then begin
           data = 0.
       endif else begin
           data = (-I*s_bracket(phi,psi,x,y) + w*s_bracket(psi,psi,x,y) $
                   -I*a_bracket(chi,psi,x,y)/r^3)/sqrt(b2*pp)
       endelse
       symbol = '!8u!D!9x!N!X'
       d = dimensions(/v0, _EXTRA=extra)

         
   ;===========================================
   ; radial flow
   ;===========================================
   endif else if(strcmp('vn', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        phi=phi0)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange,complex=complex, $
                        phi=phi0)
       psi = read_field('psi', x, y, t, /equilibrium, slices=time, $
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange, $
                        phi=phi0)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = (a_bracket(psi,phi,x,y)/r + s_bracket(psi,chi,x,y)) / $
             sqrt(s_bracket(psi,psi,x,y))
       endif else if (ivform eq 1) then begin
           data = (a_bracket(psi,phi,x,y)*r + s_bracket(psi,chi,x,y)/r^2) / $
             sqrt(s_bracket(psi,psi,x,y))
       endif
       symbol = '!8u!Dn!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   endif else if(strcmp('vx', name, /fold_case) eq 1) then begin

       chi_r = read_field('chi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=2, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                          phi=phi0)
       phi_z = read_field('phi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=3, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                         phi=phi0)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = -phi_z/r + chi_r
       endif else if (ivform eq 1) then begin
           data = -r*phi_z + chi_r/r^2
       endif
       symbol = '!8u!DR!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('vz', name, /fold_case) eq 1) then begin

       chi_z = read_field('chi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=3, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                         phi=phi0)
       phi_r = read_field('phi', x, y, t, slices=time,mesh=mesh,linear=linear,$
                        filename=filename, points=pts, mask=mask, op=2, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                         phi=phi0)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       if(ivform eq 0) then begin
           data = phi_r/r + chi_z
       endif else if (ivform eq 1) then begin
           data = r*phi_r + chi_z/r^2
       endif
       symbol = '!8u!DZ!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   ;===========================================
   ; radial field
   ;===========================================
   endif else if(strcmp('bn', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=2)
       psi0_z = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=3)
       psi_r = read_field('psi', x, y, t, mesh=mesh, phi=phi0, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=2)
       psi_z = read_field('psi', x, y, t, mesh=mesh, phi=phi0, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=3)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = (psi_z*psi0_r - psi_r*psi0_z)/r

        if(ntor ne 0) then begin
            f_r = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=2)
            f_z = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=3)
            data = data + complex(0.,ntor)*(f_r*psi0_r + f_z*psi0_z)
        endif

       data = data / sqrt(psi0_r^2 + psi0_z^2)
       symbol = '!8B!Dn!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; poloidal field
   ;===========================================
   endif else if(strcmp('bpol', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=2)
       psi0_z = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, mask=mask, op=3)
       psi_r = read_field('psi', x, y, t, mesh=mesh, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=2)
       psi_z = read_field('psi', x, y, t, mesh=mesh, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex, op=3)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = (psi_z*psi0_z + psi_r*psi0_r)/r

        if(ntor ne 0) then begin
            fp_r = read_field('fp', x, y, t, mesh=mesh, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=2)
            fp_z = read_field('fp', x, y, t, mesh=mesh, $
                           filename=filename, points=pts, slices=time, $
                           rrange=xrange, zrange=yrange, $
                           linear=linear, complex=complex, op=3)
            data = data + (fp_z*psi0_r - fp_r*psi0_z)
        endif

       data = data / sqrt(psi0_r^2 + psi0_z^2)
       symbol = '!8B!D!7h!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; (Major) Radial field
   ;===========================================
   endif else if(strcmp('bx', name, /fold_case) eq 1) then begin

       psi_z = read_field('psi', x, y, t, mesh=mesh, operation=3, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        linear=linear, mask=mask, phi=phi0, $
                          wall_mask=wall_mask,outval=mask_val)

       if(extsubtract eq 1 and version lt 8) then begin
           psi_z = psi_z + read_field('psi_ext', x, y, t, mesh=mesh, operation=3, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, mask=mask, phi=phi0, $
                                     wall_mask=wall_mask,outval=mask_val)
       end

       if(i3d eq 1) then begin
          f_rp = read_field('fp', x, y, t, mesh=mesh, operation=2, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)
          if(extsubtract eq 1 and version lt 8) then begin
             f_rp = f_rp + read_field('fp_ext', x, y, t, mesh=mesh, operation=2, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)
          end

       endif else if(icomplex eq 1 and time ge 0) then begin
          f_rp = read_field('fp', x, y, t, mesh=mesh, operation=2, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)
          if(extsubtract eq 1 and version lt 8) then begin
             f_rp = f_rp + read_field('fp_ext', x, y, t, mesh=mesh, operation=2, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, phi=phi0, wall_mask=wall_mask,outval=mask_val,mask=mask)
          end
       endif else f_rp = 0.
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = -psi_z/r - f_rp
       symbol = '!8B!DR!N!X'
       d = dimensions(/b0, _EXTRA=extra)

   ;===========================================
   ; Vertical field
   ;===========================================
   endif else if(strcmp('bz', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi', x, y, t, mesh=mesh, operation=2, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        linear=linear, mask=mask, phi=phi0, wall_mask=wall_mask,outval=mask_val)

       if(extsubtract eq 1 and version lt 8) then begin
           psi_r = psi_r + read_field('psi_ext', x, y, t, mesh=mesh, operation=2, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, mask=mask, phi=phi0, $
                                     wall_mask=wall_mask,outval=mask_val)
       end

       if(i3d eq 1) then begin
          f_zp = read_field('fp', x, y, t, mesh=mesh, operation=3, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0, $
                            wall_mask=wall_mask,outval=mask_val,mask=mask)
          if(extsubtract eq 1 and version lt 8) then begin
             f_zp = f_zp + read_field('fp_ext', x, y, t, mesh=mesh, operation=3, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, phi=phi0, $
                                      wall_mask=wall_mask,outval=mask_val,mask=mask)
          end

       endif else if(icomplex eq 1 and time ge 0) then begin
          f_zp = read_field('fp', x, y, t, mesh=mesh, operation=3, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0, $
                            wall_mask=wall_mask,outval=mask_val,mask=mask)
          if(extsubtract eq 1 and version lt 8) then begin
             f_zp = f_zp + read_field('fp_ext', x, y, t, mesh=mesh, operation=3, $
                                      filename=filename, points=pts, slices=time, $
                                      rrange=xrange, zrange=yrange, complex=complex, $
                                      linear=linear, phi=phi0, $
                                      wall_mask=wall_mask,outval=mask_val,mask=mask)
          end
       endif else f_zp = 0.

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = psi_r/r - f_zp
       symbol = '!8B!DZ!N!X'
       d = dimensions(/b0, _EXTRA=extra)


   ;===========================================
   ; Vertical field from plasma response
   ;===========================================
   endif else if(strcmp('bz_plasma', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi_plasma', x, y, t, mesh=mesh, operation=2, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange, complex=complex, $
                        linear=linear, mask=mask, phi=phi0, wall_mask=wall_mask)

       if(i3d eq 1) then begin
           f_zp = read_field('f_plasma', x, y, t, mesh=mesh, operation=13, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0, wall_mask=wall_mask)
       endif else if(icomplex eq 1) then begin
           f_z = read_field('f_plasma', x, y, t, mesh=mesh, operation=3, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, complex=complex, $
                            linear=linear, phi=phi0, wall_mask=wall_mask)
           f_zp = complex(0.,ntor)*f_z
       endif else f_zp = 0.
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = psi_r/r - f_zp
       symbol = '!8B!DZ!N!X'
       d = dimensions(/b0, _EXTRA=extra)
       
    ;===============================
    ; Dreicer field
    ;===============================
    endif else if(strcmp('dreicer', name, /fold_case) eq 1) then begin

       n = read_field('den', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts,/mks, $
                      rrange=xrange, zrange=yrange)
       Te = read_field('Te', x, y, t, slices=time, $
                       filename=filename, points=pts, /equilibrium, /mks,$
               rrange=xrange, zrange=yrange)

       lnL = 14.9 - 0.5*ALOG(n/1.e20) + ALOG(Te/1.e3)

       me_rest=511875    ;mc2
       norm=5.0244e-23   ;q^3/(4*!pi*eps0^2*me*c2)

       Ec = norm*n*lnL
       data = Ec*me_rest/Te 

       ; convert from mks back to normalized units
       d = dimensions(/potential,l0=-1,_EXTRA=extra)
       get_normalizations, filename=filename,b0=b0,n0=n0,l0=l0,ion=mi
       e_norm = 1.
       convert_units, e_norm, d, b0, n0, l0, mi, /mks
       data = data / e_norm

       symbol = '!8E!D!6Dreicer!N!X'

   ;===========================================
   ; poloidal flow
   ;===========================================
   endif else if(strcmp('vp', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
        
       data = (s_bracket(phi,psi,x,y)/r + a_bracket(chi,psi,x,y)) / $ 
         sqrt(s_bracket(psi,psi,x,y))
       symbol = '!8u!Dp!N!X'
       d = dimensions(/v0, _EXTRA=extra)


   ;===========================================
   ; surface flow
   ;===========================================
   endif else if(strcmp('vs', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       V   = read_field('V',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       i   = read_field('i',   x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       psipsi = s_bracket(psi,psi,x,y)
       if(ivform eq 0) then begin
           b2 = psipsi + i^2
           data = -(i*s_bracket(phi,psi,x,y) - v*psipsi $
                   +i*a_bracket(chi,psi,x,y)*r) / (r^2 * sqrt(psipsi*b2))
       endif else begin
           b2 = (psipsi + i^2)/r^2
           data = -(i*s_bracket(phi,psi,x,y) - v*psipsi $
                   + i*a_bracket(chi,psi,x,y)/r^3) / sqrt(b2*psipsi)

       endelse
       symbol = '!8u!Ds!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('v_star', name, /fold_case) eq 1) then begin
       rho_i = read_field('rho_i', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       omega_ci = read_field('omega_ci', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       kr = sqrt(s_bracket(p,p,x,y))/p

       data = rho_i^2 * omega_ci * kr
       symbol = '!8u!6!D*!N!X'
       d = dimensions(/v0, _EXTRA=extra)

   endif else if(strcmp('omega_star', name, /fold_case) eq 1) then begin
       v_star = read_field('v_star', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1

       data = v_star/r
       symbol = '!7x!6!D*!N!X'
       d = dimensions(t0=-1, _EXTRA=extra)

   endif else if(strcmp('bp_over_b', name, /fold_case) eq 1) then begin
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, mask=mask, $
                        rrange=xrange, zrange=yrange)

       bp = s_bracket(psi,psi,x,y)
       data = sqrt(bp)/(sqrt(bp + i^2))
       symbol = '!3|!5B!D!8p!N!3|/|!8B!3|!X'
       d = dimensions(_EXTRA=extra)


   ;===========================================
   ; normal curvature
   ;===========================================
   endif else if( (strcmp('normal curvature', name, /fold_case) eq 1) $
                  or (strcmp('kn', name, /fold_case) eq 1)) $
     then begin

       psi0 = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       I0 = read_field('I', x, y, t, slices=time, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)
       p0 = read_field('p', x, y, t, slices=time, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       B02 = s_bracket(psi0,psi0,x,y)/r^2 + I0*I0/r^2
       p02 = s_bracket(p0,p0,x,y)

       if(not keyword_set(linear)) then begin
          data = (p02 + s_bracket(p0,B02,x,y))/(B02*sqrt(p02))
       endif else begin
          psi1 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=icomplex)
          I1 = read_field('I', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=icomplex)
          p1 = read_field('p', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=icomplex)

          if(icomplex eq 1) then begin
             f1 = read_field('f', x, y, t, slices=time, mesh=mesh, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, complex=icomplex)
             f1p = complex(0,ntor)*f1
          endif else begin
             f1p = read_field('f', x, y, t, slices=time, mesh=mesh, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, op=21)
          endelse
          
          bb = s_bracket(psi1,psi0,x,y)/r^2 + I1*I0/r^2 + $
                 a_bracket(psi0,f1p,x,y)/r
          pp = s_bracket(p0,p1,x,y)

          ;; data = (pp +  $
          ;;         s_bracket(p0,bb,x,y) - $
          ;;         bb/B02*(2.*p02+s_bracket(p0,B02,x,y))) $
          ;;        / (B02*sqrt(p02))

          data = (pp +  $
                  s_bracket(p0,bb,x,y) + $
                  .5*s_bracket(p1,B02,x,y) - $
                  (.5*pp/p02 + bb/B02)*s_bracket(p0,B02,x,y) - $
                  2.*p02*bb/B02) $
                 / (B02*sqrt(p02))

          ;; bp0 = p0 + B02/2.
          ;; data = data - (s_bracket(p0,bb,x,y) - $
          ;;         (.5*pp/p02)*s_bracket(p0,B02,x,y) - $
          ;;         2.*bb*s_bracket(p0,bp0,x,y)/B02 + $
          ;;         s_bracket(p1,bp0,x,y)) $
          ;;        / (B02*sqrt(p02))

       endelse

       data = -data ;*(abs(p0) gt 6e-4)

       symbol = '!7j!D!8n!N!X'
       d = dimensions(l0=-1,_EXTRA=extra)

   ;===========================================
   ; geodesic curvature
   ;===========================================
   endif else if( (strcmp('geodesic curvature', name, /fold_case) eq 1) $
                  or (strcmp('kg', name, /fold_case) eq 1)) $
     then begin

       psi0 = read_field('psi', x, y, t, slices=time, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       I0 = read_field('I', x, y, t, slices=time, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)
       p0 = read_field('p', x, y, t, slices=time, $
                      filename=filename, points=pts, $
                      rrange=xrange, zrange=yrange, /equilibrium)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       B02 = s_bracket(psi0,psi0,x,y)/r^2 + I0*I0/r^2
       p02 = s_bracket(p0,p0,x,y)

       if(not keyword_set(linear)) then begin
          data = I0*a_bracket(B02,p0,x,y)/r
       endif else begin
          psi1 = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                            filename=filename, points=pts, linear=linear, $
                            rrange=xrange, zrange=yrange, complex=icomplex)
          I1 = read_field('I', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, complex=icomplex)
          p1 = read_field('p', x, y, t, slices=time, mesh=mesh, $
                          filename=filename, points=pts, linear=linear, $
                          rrange=xrange, zrange=yrange, /complex)

          if(icomplex eq 1) then begin
             f1 = read_field('f', x, y, t, slices=time, mesh=mesh, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, complex=icomplex)
             f1p = complex(0,ntor)*f1
             p1p = complex(0,ntor)*p1
          endif else begin
             f1p = read_field('f', x, y, t, slices=time, mesh=mesh, $
                             filename=filename, points=pts, linear=linear, $
                             rrange=xrange, zrange=yrange, op=21)
          endelse
          
          bb = s_bracket(psi1,psi0,x,y)/r^2 + I1*I0/r^2 + $
               a_bracket(psi0,f1p,x,y)/r
          pp = s_bracket(p0,p1,x,y)
          bbp = complex(0,ntor)*bb

          
          data = I1*a_bracket(B02,p0,x,y)/r + $
                 2.*(I0*a_bracket(bb,p0,x,y)/r - $
                     bbp*s_bracket(psi0,p0,x,y)/r^2) $
                 + (I0*a_bracket(B02,p1,x,y)/r + $
                    p1p*s_bracket(psi0,p0,x,y)/r^2) $
                 - I0*a_bracket(B02,p0,x,y)/r * (pp/p02 + 3.*bb/B02)
       endelse

       data = data / (2.*B02*sqrt(p02*B02))

       data = data ;*(abs(p0) gt 6e-4)

       symbol = '!7j!D!8g!N!X'
       d = dimensions(l0=-1,_EXTRA=extra)



   ;===========================================
   ; ideal_k
   ;===========================================
   endif else if(strcmp('k', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       n   = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.
         
       psipsi = s_bracket(psi,psi,x,y)

       data = n*(s_bracket(phi,psi,x,y) + r*a_bracket(chi,psi,x,y))/psipsi

   ;===========================================
   ; ideal omega
   ;===========================================
   endif else if(strcmp('ideal omega', name, /fold_case) eq 1) then begin

       phi = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       v   = read_field('v',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       psipsi = s_bracket(psi,psi,x,y)
       
       data = $
         (v - $
          i*(s_bracket(phi,psi,x,y) + r*a_bracket(chi,psi,x,y))/psipsi) $
         /r^2

   ;===========================================
   ; Lundquist number
   ;===========================================
   endif else if(strcmp('S', name, /fold_case) eq 1) then begin

       va = read_field('va', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       eta= read_field('eta',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       data = va/eta
       symbol = '!8S!X'
       d = dimensions(_EXTRA=extra)
       

   ;===========================================
   ; b.W.b
   ;===========================================
   endif else if(strcmp('wpar', name, /fold_case) eq 1) then begin

       u   = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       v   = read_field('v'  , x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       i   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       psipsi = s_bracket(psi,psi,x,y)
       r2b2 = (psipsi + i^2)
       com = laplacian(chi,x,y,tor=itor)

       data = $
         (s_bracket(psi,a_bracket(u,psi,x,y)/r,x,y) $
          -0.5*r*a_bracket(u,psipsi/r^2,x,y) $
          -I*r*a_bracket(psi,v/r^2,x,y) $
          +0.5*r^2*s_bracket(chi,psipsi/r^2,x,y) $
          +com*psipsi $
          -s_bracket(psi,s_bracket(psi,chi,x,y),x,y))/r2b2 $
         - com/3.
           
       if(itor eq 1) then begin
           data = data + I^2*(dx(chi,x)/r-dz(u,y)/r^2)/r2b2
       endif

   endif else if(strcmp('jpar', name, /fold_case) eq 1) then begin

       psi1_r = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=2)       
       psi1_z = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=3)
       psi1_lp = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=7)
       i1_r = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=2)       
       i1_z = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=3)
       itor = read_parameter('itor', filename=filename)
       if(itor eq 1) then begin
          ntor = read_parameter('ntor', filename=filename)
          rfac = complex(0,ntor)
          r = radius_matrix(x,y,t)
       endif else begin
          rzero = read_parameter('rzero', filename=filename)
          rfac = complex(0,ntor)/rzero
          r = 1.
       endelse

       if(ilin eq 1) then begin
          psi0_r = read_field('psi', x, y, t, /equilibrium, $
                              filename=filename, points=pts, slices=time, $
                              rrange=xrange, zrange=yrange,op=2)
          psi0_z = read_field('psi', x, y, t, /equilibrium, $
                              filename=filename, points=pts, slices=time, $
                              rrange=xrange, zrange=yrange,op=3)
          i0 = read_field('i', x, y, t, /equilibrium, $
                          filename=filename, points=pts, slices=time, $
                          rrange=xrange, zrange=yrange)
          f1_r = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, $
                            linear=linear, complex=complex, op=2)       
          f1_z = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, $
                            linear=linear, complex=complex, op=3)
          
          f1_rpp = f1_r*rfac^2
          f1_zpp = f1_z*rfac^2
          psi1_rp = psi1_r*rfac
          psi1_zp = psi1_z*rfac
          
          b0 = sqrt(psi0_r^2 + psi0_z^2 + i0^2)/r
                                ; J1.B0 / |B0|
          data = (-i0*(psi1_lp - itor*psi1_r/r)/r^2 $
                  +(psi0_r*(i1_r+f1_rpp) + psi0_z*(i1_z+f1_zpp))/r^2 $
                  -(psi0_z*psi1_rp - psi0_r*psi1_zp)/r^3) / B0
       endif else begin
          i1 = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                            filename=filename, points=pts, slices=time, $
                            rrange=xrange, zrange=yrange, $
                          linear=linear, complex=complex)
          if(i3d eq 1) then begin
             psi1_rp = read_field('psi', x, y, t, mesh=mesh, phi=phi0, $
                                  filename=filename, points=pts, slices=time, $
                                  rrange=xrange, zrange=yrange, $
                                  linear=linear, complex=complex, op=12)
             psi1_zp = read_field('psi', x, y, t, mesh=mesh, phi=phi0, $
                                  filename=filename, points=pts, slices=time, $
                                  rrange=xrange, zrange=yrange, $
                                  linear=linear, complex=complex, op=13)
             f1_rp = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                                filename=filename, points=pts, slices=time, $
                                rrange=xrange, zrange=yrange, $
                                linear=linear, complex=complex, op=12)
             f1_zp = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                                filename=filename, points=pts, slices=time, $
                                rrange=xrange, zrange=yrange, $
                                linear=linear, complex=complex, op=13)
             f1_rpp = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                                 filename=filename, points=pts, slices=time, $
                                 rrange=xrange, zrange=yrange, $
                                 linear=linear, complex=complex, op=22)
             f1_zpp = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                                 filename=filename, points=pts, slices=time, $
                                 rrange=xrange, zrange=yrange, $
                                 linear=linear, complex=complex, op=23)
          endif else begin
             psi1_rp = 0.
             psi1_zp = 0.
             f1_rp = 0.
             f1_zp = 0.
             f1_rpp = 0.
             f1_zpp = 0.
          end
          
          b0 = sqrt((psi1_r^2 + psi1_z^2 + i1^2)/r^2 $
                    + f1_rp^2 + f1_zp^2 + 2*(psi1_z*f1_rp - psi1_r*f1_zp)/r)
          
          data = (-i1*(psi1_lp - itor*psi1_r/r)/r^2 $
                  +(psi1_r*(i1_r+f1_rpp) + psi1_z*(i1_z+f1_zpp))/r^2 $
                  -(psi1_z*psi1_rp - psi1_r*psi1_zp)/r^3 $
                  +(f1_rp*(i1_z+f1_zpp) - f1_zp*(i1_r+f1_rpp))/r $
                  -(f1_rp*psi1_rp + f1_zp*psi1_zp)/r^2) / B0       
       endelse
       
       symbol = '!8J!D!3||!6!N!X'
       d = dimensions(j0=1,_EXTRA=extra)

   endif else if(strcmp('jpar_B', name, /fold_case) eq 1) then begin

       psi0 = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       jphi = read_field('jphi', x, y, t, linear=linear, mesh=mesh, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex)
       i0 = read_field('i', x, y, t, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, /equilibrium, $
                      complex=complex)
       
       r = radius_matrix(x,y,t)
       b0 = s_bracket(psi0,psi0,x,y)/r^2 + i0^2/r^2
       data = (s_bracket(i,psi0,x,y)/r^2 - jphi*i0/r^2)/b0

       symbol = '!8J!D!9#!N!3/!8B!X'
       d = dimensions(j0=1,b0=-1,_EXTRA=extra)

   endif else if(strcmp('jpar_lin', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange,op=2)
       psi0_z = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange,op=3)
       psi0_lp = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange,op=7)
       i0 = read_field('i', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)
       i0_r = read_field('i', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange,op=2)
       i0_z = read_field('i', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange,op=3)

       psi1_r = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=2)       
       psi1_z = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=3)
       psi1_lp = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=7)
       i1 = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex)       
       i1_r = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=2)       
       i1_z = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=3)
       f1_r = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=2)       
       f1_z = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=3)
       itor = read_parameter('itor', filename=filename)
       if(itor eq 1) then begin
          ntor = read_parameter('ntor', filename=filename)
          rfac = complex(0,ntor)
          r = radius_matrix(x,y,t)
       endif else begin
          rzero = read_parameter('rzero', filename=filename)
          rfac = complex(0,ntor)/rzero
          r = 1.
       endelse
       f1_rp = f1_r*rfac
       f1_zp = f1_z*rfac
       f1_rpp = f1_r*rfac^2
       f1_zpp = f1_z*rfac^2
       psi1_rp = psi1_r*rfac
       psi1_zp = psi1_z*rfac

       b0 = sqrt(psi0_r^2 + psi0_z^2 + i0^2)/r
       data = 0.
       ; J1.B0 / |B0|
       data = (-i0*(psi1_lp - itor*psi1_r/r)/r^2 $
               +(psi0_r*(i1_r+f1_rpp) + psi0_z*(i1_z+f1_zpp))/r^2 $
               -(psi0_z*psi1_rp - psi0_r*psi1_zp)/r^3) / B0

       ; J0.B1 / |B0|
       data = data + $
              (-i1*(psi0_lp - itor*psi0_r/r)/r^2 $
               +(psi1_r*i0_r + psi1_z*i0_z)/r^2 $
               +(i0_z*f1_rp - i0_r*f1_zp)/r)/ B0

       b0b1 = (psi0_r*psi1_r + psi0_z*psi1_z + i0*i1)/r^2 $
              + (psi0_z*f1_rp - psi0_r*f1_zp)/r

       ; - J0.B0 B1.B0 / |B0|^3
       data = data - $
              (-i0*(psi0_lp - itor*psi0_r/r) + $
               (psi0_r*i0_r + psi0_z*i0_z)/r^2) * b0b1 / B0^3

       symbol = '!8J!D!3||!6!N!X'
       d = dimensions(j0=1,_EXTRA=extra)


   endif else if(strcmp('jb', name, /fold_case) eq 1) then begin

       psi0_r = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange,op=2)
       psi0_z = read_field('psi', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange,op=3)
       i0 = read_field('i', x, y, t, /equilibrium, $
                        filename=filename, points=pts, slices=time, $
                        rrange=xrange, zrange=yrange)

       psi1_r = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=2)       
       psi1_z = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=3)
       psi1_lp = read_field('psi', x, y, t, mesh=mesh, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, phi=phi0, $
                      linear=linear, complex=complex, op=7)
       i1_r = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=2)       
       i1_z = read_field('i', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=3)
       f1_r = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=2)       
       f1_z = read_field('f', x, y, t, mesh=mesh, phi=phi0, $
                      filename=filename, points=pts, slices=time, $
                      rrange=xrange, zrange=yrange, $
                      linear=linear, complex=complex, op=3)
       itor = read_parameter('itor', filename=filename)
       if(itor eq 1) then begin
          ntor = read_parameter('ntor', filename=filename)
          rfac = complex(0,ntor)
          r = radius_matrix(x,y,t)
       endif else begin
          rzero = read_parameter('rzero', filename=filename)
          rfac = complex(0,ntor)/rzero
          r = 1.
       endelse
       f1_rpp = f1_r*rfac^2
       f1_zpp = f1_z*rfac^2
       psi1_rp = psi1_r*rfac
       psi1_zp = psi1_z*rfac

       ; J1.B0
       data = (-i0*(psi1_lp - itor*psi1_r/r)/r^2 $
               +(psi0_r*(i1_r+f1_rpp) + psi0_z*(i1_z+f1_zpp))/r^2 $
               -(psi0_z*psi1_rp - psi0_r*psi1_zp)/r^3)

       symbol = '!8J!D!9#!N!8B!X'
       d = dimensions(j0=1,b0=1,_EXTRA=extra)

   ;===========================================
   ; diffusive particle flux
   ;===========================================
   endif else if(strcmp('flux_denm', name, /fold_case) eq 1) then begin

      denm = read_field('denm', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                       filename=filename, points=pts, complex=complex, $
                       rrange=xrange, zrange=yrange, phi=phi0)
      if(max(denm) eq 0.) then begin
         print, 'denm field not found; using denm parameter'
         denm = read_parameter('denm', filename=filename)
      end
      den = read_field('den', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                       filename=filename, points=pts, complex=complex, $
                       rrange=xrange, zrange=yrange, phi=phi0)
      psi = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                       filename=filename, points=pts, complex=complex, $
                       rrange=xrange, zrange=yrange, phi=phi0)

      data = denm*s_bracket(den, psi, x, y) / sqrt(s_bracket(psi, psi, x, y))


       symbol = '!7C!D!8D!9G!8n!N!X'
       d = dimensions(/n0,/l0,t0=-1,_EXTRA=extra)

   ;===========================================
   ; particle flux
   ;===========================================
   endif else if(strcmp('flux_nv', name, /fold_case) eq 1) then begin

       den = read_field('den', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)

       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       f_p = read_field('f', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=11)
       
       r = radius_matrix(x,y,t)

       bp2 = s_bracket(psi,psi,x,y)/r^2 + 2.*a_bracket(psi,f_p,x,y)/r + s_bracket(f_p,f_p,x,y)
       bpdotv = (r*s_bracket(phi,psi,x,y) + a_bracket(chi,psi,x,y)/r^3 $
                 + a_bracket(f_p,psi,x,y)/r - s_bracket(chi,f_p,x,y)/r^2)/sqrt(bp2)
       data = den*bpdotv
       d = dimensions(/n0, /v0)
       symbol = '!6Parallel Particle Flux!X'

   ;===========================================
   ; particle flux
   ;===========================================
   endif else if(strcmp('ExB_flux_n', name, /fold_case) eq 1) then begin

       den = read_field('den', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ex = read_field('E_R', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ey = read_field('E_phi', x, y, t, slices=time, mesh=mesh,linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ez = read_field('E_Z', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)

       Bx0 = read_field('Bx', x, y, t, slices=-1,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       By0 = read_field('By', x, y, t, slices=-1,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Bz0 = read_field('Bz', x, y, t, slices=-1, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       
       B2 = Bx0^2 + By0^2 + Bz0^2

       ; (B x y) is outward for co-IP, inward for counter-IP
       ; (E x B).(B x y) / (B^2 |B x y|)
       ; = [(E.B)(B.y) - (E.y)(B.B)]/[B^2 (B^2 - (B.y)^2)^(1/2)]
       data = conj(den)* $
              ((Ex*Bx0 + Ey*By0 + Ez*Bz0)*By0 - Ey*B2)/ $
              (B2 * sqrt(B2 - By0^2))
       d = dimensions(/n0, /v0)
       symbol = '!6ExB Particle Flux!X'

   ;===========================================
   ; electron energy flux
   ;===========================================
   endif else if(strcmp('ExB_flux_pe', name, /fold_case) eq 1) then begin

       pe = read_field('pe', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ex = read_field('E_R', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ey = read_field('E_phi', x, y, t, slices=time, mesh=mesh,linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ez = read_field('E_Z', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)

       Bx0 = read_field('Bx', x, y, t, slices=-1,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       By0 = read_field('By', x, y, t, slices=-1,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Bz0 = read_field('Bz', x, y, t, slices=-1, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       
       B2 = Bx0^2 + By0^2 + Bz0^2

       ; (B x y) is outward for co-IP, inward for counter-IP
       ; (E x B).(B x y) / (B^2 |B x y|)
       ; = [(E.B)(B.y) - (E.y)(B.B)]/[B^2 (B^2 - (B.y)^2)^(1/2)]
       data = conj(pe)* $
              ((Ex*Bx0 + Ey*By0 + Ez*Bz0)*By0 - Ey*B2)/ $
              (B2 * sqrt(B2 - By0^2))
       d = dimensions(/n0, /v0)
       symbol = '!6ExB Electron Energy Flux!X'

   ;===========================================
   ; Ion energy flux
   ;===========================================
   endif else if(strcmp('ExB_flux_pi', name, /fold_case) eq 1) then begin

       pi = read_field('pi', x, y, t, slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ex = read_field('E_R', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ey = read_field('E_phi', x, y, t, slices=time, mesh=mesh,linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Ez = read_field('E_Z', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)

       Bx0 = read_field('Bx', x, y, t, slices=-1,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       By0 = read_field('By', x, y, t, slices=-1,$
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       Bz0 = read_field('Bz', x, y, t, slices=-1, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       
       B2 = Bx0^2 + By0^2 + Bz0^2

       ; (B x y) is outward for co-IP, inward for counter-IP
       ; (E x B).(B x y) / (B^2 |B x y|)
       ; = [(E.B)(B.y) - (E.y)(B.B)]/[B^2 (B^2 - (B.y)^2)^(1/2)]
       data = conj(pi)* $
              ((Ex*Bx0 + Ey*By0 + Ez*Bz0)*By0 - Ey*B2)/ $
              (B2 * sqrt(B2 - By0^2))
       d = dimensions(/n0, /v0)
       symbol = '!6ExB Ion Energy Flux!X'


   ;===========================================
   ; particle flux
   ;===========================================
   endif else if(strcmp('flux_nv2', name, /fold_case) eq 1) then begin

       den=read_field('den', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       u = read_field('phi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       chi=read_field('chi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       psi0=read_field('psi', x, y, t, slices=time, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, /equilibrium)

       
       r = radius_matrix(x,y,t)

       data = -conj(den)*(r*a_bracket(psi0,u,x,y) + s_bracket(psi0,chi,x,y)/r)
       data = data / sqrt(s_bracket(psi0,psi0,x,y))
       d = dimensions(l0=-3,t0=-1)
       symbol = '!6Particle Flux!X'

   ;===========================================
   ; parallel power flux
   ;===========================================
    endif else if((strcmp('parallel heat flux', name, /fold_case) eq 1) or $
       (strcmp('qpar', name, /fold_case) eq 1)) then begin

       p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange, phi=phi0)
       den = read_field('den', x, y, t, slices=time, mesh=mesh,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       p_p = read_field('p', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0, op=11)
       den_p = read_field('den', x, y, t, slices=time,mesh=mesh,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0, op=11)
;;        te = read_field('te', x, y, t, slices=time, mesh=mesh, $
;;                         filename=filename, points=pts, complex=complex, $
;;                         rrange=xrange, zrange=yrange, phi=phi0)
;;        te_p = read_field('te', x, y, t, slices=time,mesh=mesh,$
;;                         filename=filename, points=pts, complex=complex, $
;;                         rrange=xrange, zrange=yrange, phi=phi0, op=11)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       f_p = read_field('f', x, y, t, slices=time, mesh=mesh, phi=phi0, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=11)
       kappar = read_parameter(filename=filename, 'kappar')
       zeff = read_parameter(filename=filename, 'zeff')

       r = radius_matrix(x,y,t)

;       den = 1.

;;        if(keyword_set(linear)) then begin
;;           print, "LINEAR IS SET"

;;           p0 = read_field('p', x, y, t, slice=-1, $
;;                           filename=filename, points=pts, $
;;                           rrange=xrange, zrange=yrange, phi=phi0)
;;           den0 = read_field('den', x, y, t, slice=-1,  $
;;                             filename=filename, points=pts, $
;;                             rrange=xrange, zrange=yrange, phi=phi0)
;;           psi0 = read_field('psi', x, y, t, slices=-1, phi=phi0, $
;;                             filename=filename, points=pts, $
;;                             rrange=xrange, zrange=yrange)
;;           i0 = read_field('i', x, y, t, slices=-1, phi=phi0, $
;;                           filename=filename, points=pts, $
;;                           rrange=xrange, zrange=yrange)
          
;;           p = p + p0
;;           den = den + den0
;;           psi = psi + psi0
;;           i = i + i0
          
;;           te0 = p0 / den0
;;           b20 = s_bracket(psi0,psi0,x,y)/r^2 + i0^2/r^2
;;           bdotgradte0 = a_bracket(te0, psi0, x, y)/r
;;        end

;       den = den*zeff
;       den_p = den_p*zeff
       
       te = p / den
       te_p = p_p / den ;- p*den_p / den^2
       bp2 = s_bracket(psi,psi,x,y)/r^2 + 2.*a_bracket(psi,f_p,x,y)/r + s_bracket(f_p,f_p,x,y)
       b2 = bp2 + i^2/r^2
       bdotgradte = a_bracket(te, psi, x, y)/r $
                    + i*te_p/r^2 - s_bracket(f_p, te, x, y)

       br = -dz(psi,y)/r - dx(f_p,x)
       bbter = br*bdotgradte/b2
       bz =  dx(psi,x)/r - dz(f_p,y)
       bbtez = bz*bdotgradte/b2
;;        if(keyword_set(linear)) then begin
;;           br0 = -dz(psi0,y)/r
;;           bz0 =  dx(psi0,x)/r
;;           bbter = 0.*bbter - br0*bdotgradte0/b20
;;           bbtez = 0.*bbtez - bz0*bdotgradte0/b20
;;        endif
       
       if(keyword_set(rvector)) then begin
          data = -kappar*bbter
          symbol = '!6q!D!9#!N.G!8R!X'
       endif else if(keyword_set(zvector)) then begin
          data = -kappar*bbtez
          symbol = '!6q!D!9#!N.G!8Z!X'
       endif else begin
;          data = kappar*sqrt(abs(bbtez)^2 + abs(bbter)^2)
;          symbol = '!3|!6q!D!9#!N!3|!X'
;          data = -kappar*bdotgradte/sqrt(b2)
;          data = -kappar*bdotgradte*bp2/b2
          data = -kappar*bdotgradte*sqrt(bp2)/b2
          symbol = '!6q!D!9#!N!X'
       endelse
       
       d = dimensions(/p0,/v0)
       is_nonlinear = 1

   ;===========================================
   ; convective power flux
   ;===========================================
    endif else if((strcmp('convective heat flux', name, /fold_case) eq 1) or $
       (strcmp('qcon', name, /fold_case) eq 1)) then begin

       p=read_field('p', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       phi = read_field('phi', x, y, t, slices=time, mesh=mesh,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, phi=phi0)
       gamma = read_parameter(filename=filename, 'gam')

       r = radius_matrix(x,y,t)

       if(keyword_set(rvector)) then begin
          v = -r*dz(phi,y) + dx(chi,x)/r^2
          symbol = '!6q!D!8V!N!9.G!8R!X'
       endif else if(keyword_set(zvector)) then begin
          v =  r*dx(phi,x) + dz(chi,y)/r^2
          symbol = '!6q!D!8V!N!9.G!8Z!X'
       endif else begin
          omega = read_field('omega', x, y, t, slices=time, mesh=mesh,  $
                             filename=filename, points=pts, complex=complex, $
                             rrange=xrange, zrange=yrange)
          v =  sqrt(r^2*s_bracket(phi,phi,x,y) + r^2*omega^2 + s_bracket(chi,chi,x,y)/r^4 $
            + 2.*a_bracket(chi,phi,x,y)/r)
          symbol = '!3|!6q!D!9#!N!3|!X'
       endelse
       
       data = -gamma/(gamma-1.) * p * v
       d = dimensions(/p0,/v0)


   ;===========================================
   ; dbndt
   ;===========================================
   endif else if(strcmp('dbndt', name, /fold_case) eq 1) then begin

;        psi0 = read_field('psi',x,y,t,slices=time,linear=linear,  $
;                          filename=filename, points=pts, /equilibrium, $
;                          rrange=xrange, zrange=yrange)
       psi0_r = read_field('psi',x,y,t,slices=time,linear=linear,  $
                           filename=filename, points=pts, /equilibrium, $
                           rrange=xrange, zrange=yrange,op=2)
       psi0_z = read_field('psi',x,y,t,slices=time,linear=linear,  $
                           filename=filename, points=pts, /equilibrium, $
                           rrange=xrange, zrange=yrange,op=3)
       i0 = read_field('i', x, y, t, slices=time, linear=linear,  $
                       filename=filename, points=pts, /equilibrium, $
                       rrange=xrange, zrange=yrange)
       i0_r = read_field('i',x,y,t, slices=time, linear=linear,  $
                         filename=filename, points=pts, /equilibrium, $
                         rrange=xrange, zrange=yrange,op=2)
       i0_z = read_field('i', x,y,t, slices=time, linear=linear,  $
                         filename=filename, points=pts, /equilibrium, $
                         rrange=xrange, zrange=yrange,op=3)
       w0 = read_field('omega',x,y,t,slices=time, linear=linear,  $
                       filename=filename, points=pts, /equilibrium, $
                       rrange=xrange, zrange=yrange)

       psi = read_field('psi',x,y,t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       f = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       u = read_field('phi', x,y,t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
;        w = read_field('omega', x,y,t, slices=time, mesh=mesh, linear=linear,$
;                       filename=filename, points=pts, complex=complex, $
;                       rrange=xrange, zrange=yrange)
;        chi = read_field('chi',x,y,t, slices=time, mesh=mesh, linear=linear,$
;                         filename=filename, points=pts, complex=complex, $
;                         rrange=xrange, zrange=yrange)
       w_r = read_field('omega', x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange,op=2)
       w_z = read_field('omega', x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange,op=3)
       chi_lp = read_field('chi',x,y,t,slices=time, mesh=mesh, linear=linear,$
                           filename=filename, points=pts, complex=complex, $
                           rrange=xrange, zrange=yrange,op=7)
       chi_r = read_field('chi',x,y,t, slices=time, mesh=mesh, linear=linear,$
                          filename=filename, points=pts, complex=complex, $
                          rrange=xrange, zrange=yrange,op=2)
       chi_z = read_field('chi',x,y,t, slices=time, mesh=mesh, linear=linear,$
                          filename=filename, points=pts, complex=complex, $
                          rrange=xrange, zrange=yrange,op=3)

       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)

        data = -(i0*(chi_lp - chi_r/r)/r^2 $
                 + (i0_r*chi_r + i0_z*chi_z)/r^2 - 2.*i0*chi_r/r^3 $
                 + r^2*w0*delstar(rfac*f,x,y,tor=itor) $
                 + s_bracket(w0*r^2,rfac*f,x,y))/r $
          - a_bracket(i0,u,x,y) $
          + a_bracket(w0,psi,x,y) $
          + (w_z*psi0_r - w_r*psi0_z)

;         data = -(i0*delstar(chi,x,y,tor=itor)/r^2 $
;                  + s_bracket(i0/r^2,chi,x,y) $
;                  + r^2*w0*delstar(rfac*f,x,y,tor=itor) $
;                  + s_bracket(w0*r^2,rfac*f,x,y))/r $
;           - a_bracket(i0,u,x,y) $
;           + a_bracket(w0,psi,x,y) $
;           + a_bracket(w,psi0,x,y)


;       data = a_bracket(psi0, r*a_bracket(u, psi0, x, y) $
;                        - s_bracket(chi, psi0, x, y)/r^2, x, y)/r $
;         - w0*a_bracket(psi0, rfac*psi, x, y)/r $
;         + i0*a_bracket(psi0, rfac*u, x, y)/r $
;         + w0*s_bracket(psi0, rfac^2*f, x, y) $
;         + i0*s_bracket(psi0, rfac*chi, x, y)/r^4
;       data = data / sqrt(s_bracket(psi0,psi0,x,y))
       d = dimensions(/b0, t0=-1)
       symbol = '!6Curl[VxB]!X'

   ;===========================================
   ; curletaj
   ;===========================================
   endif else if(strcmp('curletaj', name, /fold_case) eq 1) then begin

       eta = read_field('eta',x,y,t,slices=time, linear=linear,  $
                        filename=filename, points=pts, /equilibrium, $
                        rrange=xrange, zrange=yrange)

       psi = read_field('psi',x,y,t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       i = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       f = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)
       

       data = - s_bracket(eta, i + rfac^2*f, x, y)/r $
         - eta*delstar(i + rfac^2*f, x, y, tor=itor)/r $
         + a_bracket(eta/r^2, rfac*psi, x, y)

       d = dimensions(/b0, t0=-1)
       symbol = '!6Curl[eta J]!X'


   ;===========================================
   ; torque
   ;===========================================
   endif else if(strcmp('torque', name, /fold_case) eq 1) then begin

       force_phi = read_field('force_phi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                              filename=filename, points=pts, complex=complex, $
                              rrange=xrange, zrange=yrange)
       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       
       data = force_phi*r

       d = dimensions(/p0)
       symbol = '!6Beam Torque!X'

   ;===========================================
   ; Force density - Radial
   ;===========================================
   endif else if(strcmp('JxB_x', name, /fold_case) eq 1) then begin
      by = read_field('by',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask)
      bz = read_field('bz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask)
      jy = read_field('jy_plasma',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask)
      jz = read_field('jz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask)
      
      data = (jy*bz - jz*by)
      d = dimensions(/p0, l0=-1)
      symbol = '!5J!9X!5B!9.!8R!X'

   endif else if(strcmp('JyBz', name, /fold_case) eq 1) then begin
      bz = read_field('bz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask,outval=mask_val,mask=mask)
      jy = read_field('jy_plasma',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask,outval=mask_val,mask=mask)

      data = (jy*bz)
      d = dimensions(/p0, l0=-1)
      symbol = '!8J!D!9P!N!8B!DZ!X'

   endif else if(strcmp('JzBy', name, /fold_case) eq 1) then begin
      by = read_field('by',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask)
      jz = read_field('jz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask)

      data = -(jz*by)
      d = dimensions(/p0, l0=-1)
      symbol = '-!8J!DZ!N!8B!D!9P!X'

   ;===========================================
   ; Force density - toroidal
   ;===========================================
   endif else if(strcmp('JxB_y', name, /fold_case) eq 1) then begin
      bx = read_field('bx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)
      bz = read_field('bz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)
      jx = read_field('jx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)
      jz = read_field('jz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)
      
      data = (jz*conj(bx) - jx*conj(bz))
      d = dimensions(/p0, l0=-1)
      symbol = '!5J!9X!5B!9.!9P!X'

   endif else if(strcmp('JzBx', name, /fold_case) eq 1) then begin
      bx = read_field('bx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)
      jz = read_field('jz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)

      data = (jz*conj(bx))
      d = dimensions(/p0, l0=-1)
      symbol = '!8J!DZ!N!8B!DR!X'

   endif else if(strcmp('JxBz', name, /fold_case) eq 1) then begin
      bz = read_field('bz',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)
      jx = read_field('jx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,linear=linear,$
                     complex=complex,wall_mask=wall_mask)

      data = -(jx*conj(bz))
      d = dimensions(/p0, l0=-1)
      symbol = '-!8J!DR!N!8B!DZ!X'

   ;===========================================
   ; Total vertical force density
   ;===========================================
   endif else if(strcmp('JxB_z', name, /fold_case) eq 1) then begin
      bx = read_field('bx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)
      by = read_field('by',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)
      jx = read_field('jx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)
      jy = read_field('jy_plasma',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)
      
      data = (jx*by - jy*bx)
      d = dimensions(/p0, l0=-1)
      symbol = '!5J!9X!5B!9.!8Z!X'

   ;===========================================
   ; Vertical force contribution 1
   ;===========================================
      endif else if(strcmp('JxBy', name, /fold_case) eq 1) then begin
      by = read_field('by',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)
      jx = read_field('jx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)
 
      data = (jx*by)
      d = dimensions(/p0, l0=-1)
      symbol = '!8J!DR!N!8B!D!9P!X'

   ;===========================================
   ; Vertical force contribution 2
   ;===========================================
      endif else if(strcmp('JyBx', name, /fold_case) eq 1) then begin
      bx = read_field('bx',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)
      jy = read_field('jy_plasma',x,y,t,filename=filename,slices=time,mesh=mesh,$
                      phi=phi0, rrange=xrange,zrange=yrange,points=pts,wall_mask=wall_mask, mask=mask)

      data = -(jy*bx)
      d = dimensions(/p0, l0=-1)
      symbol = '-!8J!D!9P!N!8B!DR!X'

   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_b2', name, /fold_case) eq 1) then begin

       psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=2)
       psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=3)
       i_r = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=2)
       i_z = read_field('i', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=3)
       f_r = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=2)
       f_z = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear,  $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, op=3)
       
       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)

;          data = -(psi_r*conj(rfac*psi_r) + psi_z*conj(rfac*psi_z))/r^2 $
;            -     (conj(psi_r)*(rfac*psi_r) + conj(psi_z)*(rfac*psi_z))/r^2 $
;            + (conj(rfac*f_z)*(rfac*psi_r) - conj(rfac*f_r)*(rfac*psi_z))/r $
;            + ((rfac*f_z)*conj(rfac*psi_r) - (rfac*f_r)*conj(rfac*psi_z))/r $
;            + (conj(i_z+rfac^2*f_z)*psi_r - conj(i_r+rfac^2*f_r)*psi_z)/r $
;            + ((i_z+rfac^2*f_z)*conj(psi_r) - (i_r+rfac^2*f_r)*conj(psi_z))/r $
;            - (conj(i_r+rfac^2*f_r)*(rfac*f_r) + conj(i_z+rfac^2*f_z)*(rfac*f_z)) $
;            - ((i_r+rfac^2*f_r)*conj(rfac*f_r) + $
;              (i_z+rfac^2*f_z)*conj(rfac*f_z))

       data = (conj(i_z)*psi_r - conj(i_r)*psi_z)/r $
         - (conj(i_r)*(rfac*f_r) + conj(i_z)*(rfac*f_z)) $
         + (i_z*conj(psi_r) - i_r*conj(psi_z))/r $
         - (i_r*conj(rfac*f_r) + i_z*conj(rfac*f_z))

       data = data / 2. ; factor of 2 is from toroidal average
       
       d = dimensions(/p0)
       symbol = '!6Magnetic Torque!X'


   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_b1', name, /fold_case) eq 1) then begin

        psi0_r = read_field('psi', x, y, t, /equilibrium, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=2)
        psi0_z = read_field('psi', x, y, t, /equilibrium, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=3)
        i0_r = read_field('i', x, y, t, /equilibrium, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=2)
        i0_z = read_field('i', x, y, t, /equilibrium, $
                         filename=filename, points=pts, $
                         rrange=xrange, zrange=yrange, op=3)
        psi_r = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=2)
        psi_z = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                           filename=filename, points=pts, complex=complex, $
                           rrange=xrange, zrange=yrange, op=3)
        i_r = read_field('i', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=2)
        i_z = read_field('i', x, y, t, slices=time, mesh=mesh, $
                           linear=linear,  $
                           filename=filename, points=pts, complex=complex, $
                           rrange=xrange, zrange=yrange, op=3)
        f_r = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=2)
        f_z = read_field('f', x, y, t, slices=time, mesh=mesh, linear=linear, $
                         filename=filename, points=pts, complex=complex, $
                         rrange=xrange, zrange=yrange, op=3)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0., ntor)
 
       data = -(psi0_r*(rfac*psi_r) + psi0_z*(rfac*psi_z))/r^2 $
         + ((i_z + rfac^2*f_z)*psi0_r - (i_r + rfac^2*f_r)*psi0_z)/r $
         + (i0_z*psi_r - i0_r*psi_z)/r $
         - (i0_r*rfac*f_r + i0_z*rfac*f_z)
       
       d = dimensions(/p0)
       symbol = '!7s!X'

   endif else if(strcmp('torque_p', name, /fold_case) eq 1) then begin

       p = read_field('p', x, y, t, slices=time, mesh=mesh, linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       
       rfac = complex(0., ntor)

       data = -rfac*p
       
       d = dimensions(/p0)
       symbol = '!7s!X'

   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_mu', name, /fold_case) eq 1) then begin

       mu = read_field('visc', x, y, t, slices=time,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange,/equilibrium)
       mu_c = read_field('visc_c',x,y,t, slices=time,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange,/equilibrium)

       w = read_field('omega', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.

       data = mu*delstar(r^2*w,x,y,tor=itor) $
         + r^2*s_bracket(mu,w,x,y)

       if(ntor ne 0) then begin
           u = read_field('phi',x,y,t,slices=time, mesh=mesh, linear=linear, $
                          filename=filename, points=pts, complex=complex, $
                          rrange=xrange, zrange=yrange)
           chi = read_field('chi',x,y,t,slices=time, mesh=mesh, linear=linear,$
                            filename=filename, points=pts, complex=complex, $
                            rrange=xrange, zrange=yrange)
           rfac = complex(0.,1.)*ntor

           data = data + 2.*mu_c*rfac^2*w $
             - 4.*rfac*(mu_c*dz(u, y)) $
             - 2.*rfac*(mu-mu_c)*delstar(chi,x,y,tor=itor) $
             + mu*rfac*laplacian(chi,x,y,tor=itor)/r^2 $
             + a_bracket(mu, rfac*chi, x, y)/r^2
       endif
       
       d = dimensions(/p0)
       symbol = '!7s!X'

   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_v1', name, /fold_case) eq 1) then begin

       w0 = read_field('omega', x, y, t, slices=time,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)
       n0 = read_field('den', x, y, t, slices=time,linear=linear, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange, /equilibrium)

       u = read_field('phi',x,y,t,slices=time, mesh=mesh, linear=linear, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       w = read_field('omega',x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       rfac = complex(0.,1.)*ntor
       
       data = -a_bracket(r^4*n0*w0,u,x,y)/r - 2.*r^2*n0*w0*rfac*w
       d = dimensions(/p0)
       symbol = '!7s!X'

   ;===========================================
   ; toroidal angular momentum flux
   ;===========================================
   endif else if(strcmp('torque_vv2', name, /fold_case) eq 1) then begin

       psi0 = read_field('psi', x, y, t, slices=time,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange,/equilibrium)
       n0 = read_field('den', x, y, t, slices=time,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange, /equilibrium)


       u = read_field('phi',x,y,t,slices=time, mesh=mesh, linear=linear, $
                      filename=filename, points=pts, complex=complex, $
                      rrange=xrange, zrange=yrange)
       w = read_field('omega', x, y, t, slices=time,mesh=mesh,linear=linear, $
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi',x,y,t,slices=time, mesh=mesh, linear=linear,$
                        filename=filename, points=pts, complex=complex, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then r = radius_matrix(x,y,t) else r = 1.
       
       data = -n0*conj(w)* $
         (r^3*a_bracket(psi0,u,x,y) + s_bracket(psi0,chi,x,y)) $
         / sqrt(s_bracket(psi0,psi0,x,y))
       d = dimensions(/p0, l0=1)
       symbol = '!6Angular Momentum Flux!X'

   ;===========================================
   ; Cole NTV
   ;===========================================
   endif else if(strcmp('cole_ntv', name, /fold_case) eq 1) then begin
       if(time eq -1) then return, 0
       psi0 = read_field('psi',x,y,t,filename=filename,points=pts,$
                         /equilibrium,_EXTRA=extra)
       i0 = read_field('I',x,y,t,filename=filename,points=pts,$
                         /equilibrium,_EXTRA=extra)
       w0 = read_field('omega',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       p0 = read_field('p',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       pe0 = read_field('pe',x,y,t,filename=filename,points=pts,$
                        /equilibrium,_EXTRA=extra)
       n0 = read_field('den',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       ne0 = read_field('ne',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       psi = read_field('psi',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)
       i = read_field('I',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)
       f = read_field('f',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)
       Te1 = read_field('Te',x,y,t,filename=filename,points=pts,$
                        slice=time,/linear,complex=icomplex,_EXTRA=extra)


       db = read_parameter('db',filename=filename)
       ntor = read_parameter('ntor',filename=filename)
       ion_mass = read_parameter('ion_mass',filename=filename)
       n0_norm = read_parameter('n0_norm',filename=filename)
       B0_norm = read_parameter('B0_norm',filename=filename)
       L0_norm = read_parameter('l0_norm',filename=filename)
       v0_norm = B0_norm/sqrt(4.*!pi*n0_norm*ion_mass*1.6726e-24)
       t0_norm = L0_norm/v0_norm
       lambda = 16.              ; coulomb logarithm

       print, 'n0, B0, L0, v0, t0 = ', n0_norm, B0_norm, L0_norm, v0_norm, t0_norm

       r = radius_matrix(x,y,t)
       gradpsi = sqrt(s_bracket(psi0,psi0,x,y))

       pi0 = p0-pe0
       Ti = pi0/n0 > 0.01
       Te = pe0/(ne0) > 0.01

      tprime = abs(s_bracket(Te,psi0,x,y)) 
      xi = -Te1/tprime*gradpsi
       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then begin
           xi(where(abs(tprime) lt 1e-6 or psi0 lt psis)) = 0.
       endif else begin
           xi(where(abs(tprime) lt 1e-6 or psi0 gt psis)) = 0.
       endelse

       forward_function flux_average_field
       fa = flux_average_field(r^2,psi0,x,y,t,r0=r0,flux=flux,$
                               filename=filename,integrate=0,_EXTRA=extra)
       r2_av = interpol(fa, flux, psi0)
       
       epsilon = read_field('minor radius',x,y,t,filename=filename,points=pts,$
                            /equilibrium,_EXTRA=extra)/r0

       w_B = abs(db*Ti * s_bracket(epsilon,psi0,x,y)/gradpsi^2)
       wstar = w0 - db*s_bracket(pi0,psi0,x,y)/gradpsi^2 / n0
       w_E = abs(wstar)
       vti = sqrt(2.*Ti)
       
;        nu_i = (64.*!pi^(5/2)/3.)*(zeff*4.8032e-10*n0_norm)^4*lambda $
;          * l0_norm/(n0_norm*B0_norm^4) * n0/Ti^1.5
       z_ion = read_parameter('z_ion', filename=filename)
       nu_i = 1.988e-35*(n0/Ti^1.5)*lambda*z_ion^4*n0_norm^3*l0_norm/B0_norm^4
       print, 'min, max(nu_i)', min(nu_i/t0_norm), max(nu_i/t0_norm)
       print, 'min, max(w_B)', min(w_B/t0_norm), max(w_B/t0_norm)
       print, 'min, max(w_E)', min(w_E/t0_norm), max(w_E/t0_norm)
       nu_eff = nu_i / (ntor*epsilon)
       
       mu_p = 0.21*ntor*vti^2*sqrt(epsilon*nu_eff) / $
         (r2_av*(w_E^1.5 + 0.3*w_B*sqrt(nu_eff) + 0.04*nu_eff^1.5))
       
       b02 = gradpsi^2/r^2 + i0^2/r^2
       
       if(icomplex eq 1) then begin
          ; this form takes into account the "lagrangian" perturbation
          ; c.f. https://aip.scitation.org/doi/10.1063/1.3122862
          b1b0 = (s_bracket(psi,psi0,x,y)/r^2 $
            + i*i0/r^2 $
            - complex(0.,ntor)*a_bracket(f,psi0,x,y)/r) / sqrt(b02)
          b1 = b1b0 + xi*s_bracket(sqrt(b02),psi0,x,y)/gradpsi
          b2 = b1*conj(b1) 

          ; this form is just the local peturbed field ("eulerian" component)
            ;; b2 = s_bracket(psi,conj(psi),x,y)/r^2 $
            ;;   + i*conj(i)/r^2 $
            ;;   + ntor^2*s_bracket(f,conj(f),x,y) $
            ;;   - 0.5*complex(0.,ntor)*a_bracket(f,conj(psi),x,y)/r $
            ;;   + 0.5*complex(0.,ntor)*a_bracket(conj(f),psi,x,y)/r
       end
       
       fa = flux_average_field(b2/b02,psi0,x,y,t,r0=r0,flux=flux,$
                               filename=filename,integrate=0,_EXTRA=extra)
       b2_av = interpol(fa, flux, psi0)
       
       data = -mu_p*b2_av*n0*r^2*w0
;       data = -mu_p*b2_av*n0*r^2*(w0 - wstar)
;       data = b2_av
       d = dimensions(/p0)
       symbol = '!6NTV!X'

;       data = mu_p
;       d = dimensions(t0=-1)
       
;        data = nu_i
;        data = w_E
;       data = w_E^1.5 + 0.3*w_B*sqrt(nu_eff) + 0.04*nu_eff^1.5
;       d = dimensions(t0=-1)

   ;===========================================
   ; delta_B / B
   ;===========================================

   endif else if(strcmp('dB_B', name, /fold_case) eq 1) then begin
       psi0 = read_field('psi',x,y,t,filename=filename,points=pts,$
                         /equilibrium,_EXTRA=extra)
       i0 = read_field('I',x,y,t,filename=filename,points=pts,$
                       /equilibrium,_EXTRA=extra)
       Te0 = read_field('Te',x,y,t,filename=filename,points=pts,$
                        /equilibrium,_EXTRA=extra)
       psi = read_field('psi',x,y,t,filename=filename,points=pts,$
                        slice=time,/linear,complex=icomplex,_EXTRA=extra)
       i = read_field('I',x,y,t,filename=filename,points=pts,$
                      slice=time,/linear,complex=icomplex,_EXTRA=extra)
       f = read_field('f',x,y,t,filename=filename,points=pts,$
                      slice=time,/linear,complex=icomplex,_EXTRA=extra)
       Te1 = read_field('Te',x,y,t,filename=filename,points=pts,$
                       slice=time,/linear,complex=icomplex,_EXTRA=extra)

       ntor = read_parameter('ntor',filename=filename)

       r = radius_matrix(x,y,t)
 
       gradpsi = sqrt(s_bracket(psi0,psi0,x,y))

       tprime = abs(s_bracket(Te0,psi0,x,y)) 
       xi = -Te1/tprime*gradpsi
       psis = read_lcfs(filename=filename, flux0=flux0, _EXTRA=extra)
       if(psis lt flux0) then begin
           xi(where(abs(tprime) lt 1e-6 or psi0 lt psis)) = 0.
       endif else begin
           xi(where(abs(tprime) lt 1e-6 or psi0 gt psis)) = 0.
       endelse

       b0 = sqrt(gradpsi^2 + i0^2)/r
       if(icomplex eq 1) then begin
 ;            b1b0 = complex(0.,ntor)*(s_bracket(psi,psi0,x,y)/r^2 $
 ;                    + i*i0/r^2 $
 ;                    - complex(0.,ntor)*a_bracket(f,psi0,x,y)/r) / b0
           b1b0 = (s_bracket(psi0,psi,x,y)/r^2 + i0*i/r^2 $
                   - complex(0.,ntor)*a_bracket(f,psi0,x,y)/r) / b0
           b1 = b1b0 + xi*s_bracket(b0,psi0,x,y)/gradpsi
           b2 = b1*conj(b1) 
       end
       data = b2/b0^2
       
       d = dimensions()
       symbol = '!7d!8B!6!U2!N!3/!8B!6!U2!N!X'

   ;===========================================
   ; parallel electric field
   ;===========================================
   endif else if(strcmp('epar', name, /fold_case) eq 1) then begin

       bx   = read_field('bx', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        phi=phi0,rrange=xrange, zrange=yrange)
       by   = read_field('by', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        phi=phi0,rrange=xrange, zrange=yrange)
       bz   = read_field('bz', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        phi=phi0,rrange=xrange, zrange=yrange)
       ex   = read_field('e_r', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        phi=phi0,rrange=xrange, zrange=yrange)
       ey   = read_field('e_phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        phi=phi0,rrange=xrange, zrange=yrange)
       ez   = read_field('e_z', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        phi=phi0,rrange=xrange, zrange=yrange)
       
       data = (ex*bx + ey*by + ez*bz) / sqrt(bx^2 + by^2 + bz^2)
       d = dimensions(/potential,l0=-1)
       symbol = '!8E!D!3||!N!X'

   ;===========================================
   ; radial electric field
   ;===========================================
   endif else if(strcmp('e_r', name, /fold_case) eq 1) then begin

       u   = read_field('phi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       v   = read_field('v'  , x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       chi = read_field('chi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       psi = read_field('psi', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       b   = read_field('i',   x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       pe  = read_field('pe',  x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       n   = read_field('den', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)
       eta = read_field('eta', x, y, t, slices=time, mesh=mesh, $
                        filename=filename, points=pts, $
                        rrange=xrange, zrange=yrange)

       if(itor eq 1) then begin
           r = radius_matrix(x,y,t)
       endif else r = 1.

       jphi = delstar(psi,x,y,tor=itor)
       db = read_parameter('db', filename=filename)

       data = b*s_bracket(u,psi,x,y)/r^2 $
         + b*a_bracket(chi,psi,x,y)/r $
         - v*s_bracket(psi,psi,x,y)/r^2 $
         - (db/n)* $
         (s_bracket(psi,pe,x,y) $
          + (jphi*s_bracket(psi,psi,x,y) + b*s_bracket(psi,b,x,y))/r^2) $
         + eta*a_bracket(psi,b,x,y)/r

                                ; Normalize field
       data = -data/sqrt(s_bracket(psi,psi,x,y))
       symbol = '!8E!Dr!N!X'
       d = dimensions(/pot,l0=-1, _EXTRA=extra)

       
   endif else begin
       print, 'composite field ', name, ' not found'
   endelse

   endelse

   ; scale by linfac
   if((ilin eq 1) and (n_elements(linfac) ne 0) $
      and (time ne -1) and keyword_set(linear) and not keyword_set(nolf)) $
     then begin
       print, 'scaling data by ', linfac
       data = data*linfac
   end
   
   print, 'converting units, mks, cgs=', keyword_set(mks), keyword_set(cgs)

   ; convert to appropriate units
   d0 = d
   get_normalizations, filename=filename,b0=b0,n0=n0,l0=l0,ion=mi
   convert_units, data, d0, b0, n0, l0, mi, cgs=cgs, mks=mks
   convert_units, x, dimensions(/l0), b0, n0, l0, mi, cgs=cgs, mks=mks
   convert_units, y, dimensions(/l0), b0, n0, l0, mi, cgs=cgs, mks=mks
   convert_units, realtime, dimensions(/t0), $
     b0, n0, l0, mi, cgs=cgs, mks=mks
   units = parse_units(d0, cgs=cgs, mks=mks)

   if(n_elements(h_symmetry) eq 1) then begin
       data = (data + h_symmetry*reverse(data, 2)) / 2.
   endif
   if(n_elements(v_symmetry) eq 1) then begin
       print, "v symmetry = ", v_symmetry
       data = (data + v_symmetry*reverse(data, 3)) / 2.
   endif

   ; apply mask
   if(n_elements(mask) ne 0 and n_elements(edge_val) ne 0) then begin
       data = data * (1. - mask) + mask*edge_val
   end

   if(keyword_set(complex)) then begin
       if(keyword_set(abs)) then begin
           print, 'Taking absolute value of data'
           data = abs(data)
       endif else if(keyword_set(phase)) then begin
           print, 'Taking phase of data'
           data = atan(imaginary(data),real_part(data))
       endif
   end

   if(n_elements(fac) ne 0) then begin
       print, 'applying factor = ', fac
       data = data*fac
   end


   ; perform flux-average
   if(keyword_set(flux_av)) then begin
       forward_function flux_average
       fa = flux_average(data, psi=psi, x=x, z=z, t=t, flux=flux, $
                         filename=filename, _EXTRA=extra)
       data = interpol(fa, flux, psi)
   end

   print, 'Done reading field'
   print, '**********************************************************'


   return, data
end

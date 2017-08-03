function read_scalar, scalarname, filename=filename, title=title, $
                      symbol=symbol, units=units, time=time, final=final, $
                      _EXTRA=extra

   if(n_elements(scalarname) eq 0) then begin
       print, "Error: no scalar name provided"
       return, 0
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(n_elements(filename) gt 1) then begin
       data = fltarr(n_elements(filename))
       for i=0, n_elements(filename)-1 do begin
           data[i] = read_scalar(scalarname, filename=filename[i], $
                                 title=title, symbol=symbol, units=units, $
                                 time=time, /final)
       end
       return, data
   endif

   s = read_scalars(filename=filename)
   itor = read_parameter('itor', filename=filename)
   if(n_tags(s) eq 0) then return, 0

   time = s.time._data

   d = dimensions()

   if(strcmp("toroidal current", scalarname, /fold_case) eq 1) or $
     (strcmp("it", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current._data
       title = 'Toroidal Current'
       symbol = '!8I!D!9P!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("plasma current", scalarname, /fold_case) eq 1) or $
       (strcmp("ip", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current_p._data
       title = 'Plasma Current'
       symbol = '!8I!DP!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("wall current", scalarname, /fold_case) eq 1) or $
       (strcmp("iw", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current_w._data
       title = 'Wall Current'
       symbol = '!8I!DW!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("total current", scalarname, /fold_case) eq 1) or $
       (strcmp("itot", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current_w._data + $
              s.toroidal_current._data
       title = 'Total Current'
       symbol = '!8I!D!9P!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("volume", scalarname, /fold_case) eq 1) then begin
       data = s.volume_p._data
       title = 'Plasma Volume'
       symbol = '!8V!X'
       d = dimensions(l0=3, _EXTRA=extra)
   endif else $
     if (strcmp("toroidal flux", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_flux_p._data
       title = 'Toroidal Flux'
       symbol = 'Flux'
       d = dimensions(/b0, l0=2, _EXTRA=extra)
   endif else $
     if (strcmp("reconnected flux", scalarname, /fold_case) eq 1) then begin
       data = abs(s.reconnected_flux._data)
       title = 'Reconnected Flux'
       symbol = translate('psi')
       d = dimensions(/b0, l0=1+itor, _EXTRA=extra)
   endif else $
     if (strcmp("time step", scalarname, /fold_case) eq 1) or $
        (strcmp("dt",scalarname, /fold_case) eq 1) then begin
       data = abs(s.dt._data)
       title = 'Time Step'
       symbol = '!7d!8t!X'
       d = dimensions(t0=1, _EXTRA=extra)
   endif else $
     if (strcmp("psimin", scalarname, /fold_case) eq 1) then begin
       data = s.psimin._data
       title = 'Psimin'
       symbol = '!7w!D!60!N!X'
       d = dimensions(/b0, l0=2, _EXTRA=extra)
   endif else $
     if (strcmp("psibound", scalarname, /fold_case) eq 1) or $
        (strcmp("psilim",scalarname, /fold_case) eq 1) then begin
       data = s.psi_lcfs._data
       title = 'Psilim'
       symbol = '!7w!D!8b!N!X'
       d = dimensions(/b0, l0=2, _EXTRA=extra)
   endif else $
     if (strcmp("loop voltage", scalarname, /fold_case) eq 1) or $
     (strcmp("vl", scalarname, /fold_case) eq 1) then begin
       data = s.loop_voltage._data
       title = 'Loop Voltage'
       symbol = '!8V!DL!N!X'
       d = dimensions(/pot, _EXTRA=extra)
   endif else $
     if (strcmp("pellet rate", scalarname, /fold_case) eq 1) or $
     (strcmp("pelr", scalarname, /fold_case) eq 1) then begin
       data = s.pellet_rate._data
       title = 'Pellet Rate'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, t0=-1, _EXTRA=extra)
   endif else $
 if (strcmp("pellet ablation rate", scalarname, /fold_case) eq 1) or $
     (strcmp("pelablr", scalarname, /fold_case) eq 1) then begin
       data = s.pellet_ablrate._data
       title = 'Pellet Ablation Rate'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, t0=-1, _EXTRA=extra)
    endif else $
     if (strcmp("pellet var", scalarname, /fold_case) eq 1) or $
     (strcmp("pelvar", scalarname, /fold_case) eq 1) then begin
       data = s.pellet_var._data
       title = 'Pellet Var'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, t0=-1, _EXTRA=extra)
   endif else $
     if (strcmp("pellet radius", scalarname, /fold_case) eq 1) or $
     (strcmp("pelrad", scalarname, /fold_case) eq 1) then begin
       data = s.r_p2._data
       title = 'Pellet Radius'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, t0=-1, _EXTRA=extra)
    endif else $
     if (strcmp("pellet R position", scalarname, /fold_case) eq 1) or $
     (strcmp("pelrpos", scalarname, /fold_case) eq 1) then begin
       data = s.pellet_x._data
       title = 'Pellet R position'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, t0=-1, _EXTRA=extra)
    endif else $
     if (strcmp("pellet Z position", scalarname, /fold_case) eq 1) or $
     (strcmp("pelzpos", scalarname, /fold_case) eq 1) then begin
       data = s.pellet_z._data
       title = 'Pellet Z position'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, t0=-1, _EXTRA=extra)
   endif else $
     if (strcmp("beta", scalarname, /fold_case) eq 1) then begin
       data = scalar_beta(filename=filename)
       title = 'Average Beta'
       symbol = '!7b!X'
       d = dimensions(_EXTRA=extra)
   endif else if $
     (strcmp("poloidal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bp", scalarname, /fold_case) eq 1) then begin
       data = scalar_beta_poloidal(filename=filename)
       title = 'Poloidal Beta'
       symbol = '!7b!D!8P!N!X'
       d = dimensions(_EXTRA=extra)
   endif else $
     if (strcmp("normal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bn", scalarname, /fold_case) eq 1) then begin
       data = scalar_beta_normal(filename=filename)
       title = 'Normal Beta'
       symbol = '!7b!D!8N!N!3'
       d = dimensions(_EXTRA=extra)
   endif else if $
     (strcmp("toroidal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bt", scalarname, /fold_case) eq 1) then begin
       data = scalar_beta_toroidal(filename=filename)
       title = 'Toroidal Beta'
       symbol = '!7b!D!8T!N!X'
   endif else if $
     (strcmp("kinetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("ke", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_KP._data + s.E_KT._data + s.E_K3._data
       title = 'Kinetic Energy'
       symbol = '!8KE!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("magnetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("me", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_MP._data + s.E_MT._data 
       title = 'Magnetic Energy'
       symbol = '!8ME!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("thermal energy", scalarname, /fold_case) eq 1) or $
     (strcmp("te", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_P._data 
       title = 'Thermal Energy'
       symbol = '!8TE!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("energy", scalarname, /fold_case) eq 1) then begin
       nv = read_parameter("numvar", filename=filename)
       data = energy(filename=filename)
       title = 'Total Energy'
       symbol = '!8E!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("particles", scalarname, /fold_case) eq 1) or $
     (strcmp("n", scalarname, /fold_case) eq 1) then begin
       data = s.particle_number._data
       title = 'Particle Number'
       symbol = '!8N!X'
       d = dimensions(/n0,l0=3, _EXTRA=extra)
   endif else if $
     (strcmp("angular momentum", scalarname, /fold_case) eq 1) then begin
       data = s.angular_momentum._data
       title = 'Angular Momentum'
       symbol = '!8L!D!9P!N!X'
       d = dimensions(/energy,/t0, _EXTRA=extra)
   endif else if $
     (strcmp("time", scalarname, /fold_case) eq 1) then begin
       data = s.time._data
       title = 'Time'
       symbol = '!8t!X'
       d = dimensions(/t0, _EXTRA=extra)
   endif else if (strcmp("circulation", scalarname, /fold_case) eq 1) or $
     (strcmp("vorticity", scalarname, /fold_case) eq 1) then begin
       data = s.circulation._data
;      data = s.vorticity._data
       title = 'Circulation'
       symbol = '!9I!S!7x!R!A!6_!D!9P!N !8dA !x'
       d = dimensions(/v0,/l0, _EXTRA=extra)
   endif else if (strcmp("parallel viscous heating", scalarname, /fold_case) eq 1) then begin $
       data = s.parallel_viscous_heating._data
       title = 'Parallel Viscous Heating'
;       symbol = '!6(!7l!D!9#!N!5b!9 . !3W!9 . !5b!6/2)!U2!N!X'
       symbol = '-!9I!8dV !7P!D!9#!N!3:!9G!5u!X'
       d = dimensions(/energy,t0=-1, _EXTRA=extra)
   endif else if (strcmp("bwb2", scalarname, /fold_case) eq 1) then begin $
       amupar = read_parameter('amupar', filename=filename)
       data = 4.*s.parallel_viscous_heating._data / (3.*amupar)
       title = ''
;       symbol = '!6(!7l!D!9#!N!5b!9 . !3W!9 . !5b!6/2)!U2!N!X'
       symbol = '!6(!5b!9.!3W!9.!5b!6)!U2!N!X'
       d = dimensions(t0=-2,l0=3, _EXTRA=extra)
   endif else if (strcmp("flux", scalarname, /fold_case) eq 1) then begin
       data = -2.*!pi*(s.psi_lcfs._data - s.psimin._data)
       title = 'Flux'
       symbol = '!7W!X'
       d = dimensions(/b0,l0=2, _EXTRA=extra)
   endif else if (strcmp("li", scalarname, /fold_case) eq 1) then begin
      rzero = read_parameter('rzero', filename=filename)
       data = -4.*!pi*(s.psi_lcfs._data - s.psimin._data) $
              / s.toroidal_current_p._data / rzero
       title = 'Normalized Internal Inductance'
       symbol = '!13l!Di!X'
       d = dimensions(_EXTRA=extra)
   endif else if (strcmp("xmag", scalarname, /fold_case) eq 1) then begin
       data = s.xmag._data
       title = '!8R!6-Coordinate of Magnetic Axis!6'
       symbol = '!8R!D!60!N!X'
       d = dimensions(/l0,_EXTRA=extra)
   endif else if (strcmp("zmag", scalarname, /fold_case) eq 1) then begin
       data = s.zmag._data
       title = '!8Z!6-Coordinate of Magnetic Axis!6'
       symbol = '!8Z!D!60!N!X'
       d = dimensions(/l0,_EXTRA=extra)
   endif else if (strcmp("runaways", scalarname, /fold_case) eq 1) then begin
       data = s.runaways._data
       title = '!6Runaway Electrons!6'
       symbol = '!8N!D!6RE!N!X'
       d = dimensions(/n0,l0=3,_EXTRA=extra)
   endif else if (strcmp("IZ", scalarname, /fold_case) eq 1) then begin
       data = s.M_IZ._data / s.toroidal_current_p._data
       title = '!6Plasma Current Centroid!6'
       symbol = '!8Z!DI!N!X'
       d = dimensions(/l0,_EXTRA=extra)
   endif else begin
       s = read_scalars(filename=filename)
       n = tag_names(s)
       match = where(strcmp(n, scalarname, /fold_case) eq 1,count)
       if(count eq 0) then begin
           print, 'Scalar ', scalarname, ' not recognized.'
           return, 0
       endif
       data = s.(match[0])._data
       title = ''
       symbol = scalarname
       d = dimensions()
   endelse
   
   if(keyword_set(final)) then begin
       data = data[n_elements(data)-1]
   endif

   get_normalizations, b0=b0,n0=n0,l0=l0, zeff=zeff, ion_mass=mi, $
     filename=filename, _EXTRA=extra
   convert_units, data, d, b0, n0, l0, zeff, mi, _EXTRA=extra
   convert_units, time, dimensions(/t0), b0, n0, l0, zeff, mi, _EXTRA=extra
   units = parse_units(d, _EXTRA=extra)

   if(n_elements(data) gt n_elements(time)) then $
      data = data[0:n_elements(time)-1]
   if(n_elements(time) gt n_elements(data)) then $
      time = time[0:n_elements(data)-1]

   return, data
end

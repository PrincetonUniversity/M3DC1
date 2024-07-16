function read_scalar, scalarname, filename=filename, title=title, $
                      symbol=symbol, units=units, time=time, final=final, $
                      integrate=integrate, ipellet=ipellet, _EXTRA=extra

   if(n_elements(scalarname) eq 0) then begin
       print, "Error: no scalar name provided"
       return, 0
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'
   
   if(n_elements(ipellet) eq 0) then ipellet=0

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
   rzero = read_parameter('rzero', filename=filename)
   version = read_parameter('version', filename=filename)
   threed = read_parameter('3d', filename=filename)
   if(version ge 31) then p = read_pellets(filename=filename)

   if(n_tags(s) eq 0) then return, 0

   time = s.time._data
   d = dimensions()

   if(strcmp("toroidal current", scalarname, /fold_case) eq 1) or $
     (strcmp("it", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current._data
       if(itor eq 0 and version lt 36) then begin
          print, 'WARNING: correcting for incorrect plasma current definition with itor=0 and version<36'
          data = data/rzero
       end
       title = 'Toroidal Current'
       symbol = '!8I!D!9P!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("plasma current", scalarname, /fold_case) eq 1) or $
       (strcmp("ip", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current_p._data
       if(itor eq 0 and version lt 36) then begin
          print, 'WARNING: correcting for incorrect plasma current definition with itor=0 and version<36'
          data = data/rzero
       end
       title = 'Plasma Current'
       symbol = '!8I!DP!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("wall current", scalarname, /fold_case) eq 1) or $
       (strcmp("iw", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current_w._data
       if(itor eq 0 and version lt 36) then begin
          print, 'WARNING: correcting for incorrect plasma current definition with itor=0 and version<36'
          data = data/rzero
       end
       title = 'Wall Current'
       symbol = '!8I!DW!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("total current", scalarname, /fold_case) eq 1) or $
       (strcmp("itot", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current_w._data + $
              s.toroidal_current._data
       if(itor eq 0 and version lt 36) then begin
          print, 'WARNING: correcting for incorrect plasma current definition with itor=0 and version<36'
          data = data/rzero
       end
       title = 'Total Current'
       symbol = '!8I!D!9P!N!X'
       d = dimensions(/j0, l0=2, _EXTRA=extra)
   endif else $
     if(strcmp("halo current", scalarname, /fold_case) eq 1) or $
       (strcmp("ih", scalarname, /fold_case) eq 1) then begin
       data = s.toroidal_current._data - $
              s.toroidal_current_p._data
       if(itor eq 0 and version lt 36) then begin
          print, 'WARNING: correcting for incorrect plasma current definition with itor=0 and version<36'
          data = data/rzero
       end
       title = 'Toroidal Halo Current'
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
       if(version lt 31) then begin
          data = s.pellet_rate._data
       endif else if(p eq !NULL) then begin
          print, 'Error: pellet data not present in this file'
          return, 0
       endif else begin
          if ipellet eq -1 then begin
            data = p.pellet_rate._data
          endif else begin
            data = p.pellet_rate._data[ipellet,*]
          endelse
       endelse
       title = 'Pellet Rate'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, l0=3, t0=-1, _EXTRA=extra)
     endif else $
       if (strcmp("pellet rate D2", scalarname, /fold_case) eq 1) or $
       (strcmp("pelrD2", scalarname, /fold_case) eq 1) then begin
       if(version lt 31) then begin
         data = s.pellet_rate_D2._data
       endif else if(p eq !NULL) then begin
         print, 'Error: pellet data not present in this file'
         return, 0
       endif else begin
         if ipellet eq -1 then begin
           data = p.pellet_rate_D2._data
         endif else begin
           data = p.pellet_rate_D2._data[ipellet,*]
         endelse
       endelse
       title = 'Pellet Rate D2'
       symbol = '!8V!DL!N!X'
       d = dimensions(/n0, l0=3, t0=-1, _EXTRA=extra)
   endif else $
   if (strcmp("pellet ablation rate", scalarname, /fold_case) eq 1) or $
      (strcmp("pelablr", scalarname, /fold_case) eq 1) then begin
      if(version lt 26) then begin
         data = s.pellet_ablrate._data
         title = 'Pellet Ablation Rate'
         symbol = '!8V!DL!N!X'
         d = dimensions(/n0, t0=-1, _EXTRA=extra)
      endif else begin
         print, 'Error: this data is not present in this version of M3D-C1.'
         return, 0
      end
    endif else $
     if (strcmp("pellet var", scalarname, /fold_case) eq 1) or $
     (strcmp("pelvar", scalarname, /fold_case) eq 1) then begin
       if(version lt 31) then begin
          data = s.pellet_var._data
        endif else if(p eq !NULL) then begin
          print, 'Error: pellet data not present in this file'
          return, 0
        endif else begin
          if ipellet eq -1 then begin
            data = p.pellet_var._data
          endif else begin
            data = p.pellet_var._data[ipellet,*]
          endelse
       endelse
       title = 'Pellet Var'
       symbol = '!8V!DL!N!X'
       d = dimensions(/l0)
   endif else $
     if (strcmp("pellet radius", scalarname, /fold_case) eq 1) or $
     (strcmp("pelrad", scalarname, /fold_case) eq 1) then begin
      if (version lt 26) then begin
         data = s.r_p2._data
      endif else begin
        if(version lt 31) then begin
          data = s.r_p._data
        endif else if(p eq !NULL) then begin
          print, 'Error: pellet data not present in this file'
          return, 0
        endif else begin
          if ipellet eq -1 then begin
            data = p.r_p._data
          endif else begin
            data = p.r_p._data[ipellet,*]
          endelse
        endelse
      end
       title = 'Pellet Radius'
       symbol = '!8V!DL!N!X'
       d = dimensions(/l0, _EXTRA=extra)
    endif else $
     if (strcmp("pellet R position", scalarname, /fold_case) eq 1) or $
     (strcmp("pelrpos", scalarname, /fold_case) eq 1) then begin
       if(version lt 26) then begin
          data = s.pellet_x._data
       endif else begin
         if(version lt 31) then begin
           data = s.pellet_r._data
         endif else if(p eq !NULL) then begin
           print, 'Error: pellet data not present in this file'
           return, 0
         endif else begin
           if ipellet eq -1 then begin
             data = p.pellet_r._data
           endif else begin
             data = p.pellet_r._data[ipellet,*]
           endelse
         endelse
       end
       title = 'Pellet R position'
       symbol = '!8V!DL!N!X'
       d = dimensions(/l0, _EXTRA=extra)
    endif else $
     if (strcmp("pellet phi position", scalarname, /fold_case) eq 1) or $
     (strcmp("pelphipos", scalarname, /fold_case) eq 1) then begin
       if(version lt 31) then begin
          data = s.pellet_phi._data
       endif else if(p eq !NULL) then begin
          print, 'Error: pellet data not present in this file'
          return, 0
       endif else begin
          if ipellet eq -1 then begin
            data = p.pellet_phi._data
          endif else begin
            data = p.pellet_phi._data[ipellet,*]
          endelse
       endelse
       title = 'Pellet !9P!X position'
       symbol = '!8V!DL!N!X'
       d = dimensions(/l0, _EXTRA=extra)
    endif else $
     if (strcmp("pellet Z position", scalarname, /fold_case) eq 1) or $
     (strcmp("pelzpos", scalarname, /fold_case) eq 1) then begin
       if(version lt 31) then begin
          data = s.pellet_z._data
       endif else if(p eq !NULL) then begin
          print, 'Error: pellet data not present in this file'
          return, 0
       endif else begin
          if ipellet eq -1 then begin
            data = p.pellet_z._data
          endif else begin
            data = p.pellet_z._data[ipellet,*]
          endelse
       endelse
       title = 'Pellet Z position'
       symbol = '!8V!DL!N!X'
       d = dimensions(/l0, _EXTRA=extra)
   endif else $
     if (strcmp("beta", scalarname, /fold_case) eq 1) then begin
      if(version lt 26) then begin
         print, 'Beta = int(P) / int(B^2/2); integration over XMHD region'
         data = scalar_beta(filename=filename)
         title = 'Global Beta'
      endif else begin
         print, 'Beta = int(P) / int(B^2/2); integration over plasma volume'
         gamma = read_parameter('gam', filename=filename)
         data = (gamma - 1.)*s.W_P._data / (s.E_MP._data + s.E_MT._data)
         title = 'Plasma Beta'
      endelse
       symbol = '!7b!X'
       d = dimensions(_EXTRA=extra)
   endif else if $
     (strcmp("poloidal beta", scalarname, /fold_case) eq 1) or $
     (strcmp("bp", scalarname, /fold_case) eq 1) then begin
      if(version lt 26) then begin
         print, 'Beta_poloidal = 2.*int(P) / IP^2; integration over XMHD region'
         data = scalar_beta_poloidal(filename=filename)
         title = 'Global Poloidal Beta'
      endif else begin
         print, 'Beta_poloidal = W_P / W_M; integration over plasma volume'
         gamma = read_parameter('gam',filename=filename)
         data = (gamma-1.)*s.w_p._data/s.w_m._data
         title = 'Plasma Poloidal Beta'
      endelse
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
     (strcmp("poloidal magnetic energy", scalarname, /fold_case) eq 1) or $
     (strcmp("Wm", scalarname, /fold_case) eq 1)then begin
       nv = read_parameter("numvar", filename=filename)
       data = s.E_MP._data
       title = 'Poloidal Magnetic Energy'
       symbol = '!8W!Dm!N!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("thermal energy", scalarname, /fold_case) eq 1) or $
     (strcmp("p", scalarname, /fold_case) eq 1)then begin
       data = s.E_P._data 
       title = 'Thermal Energy'
       symbol = '!8TE!X'
       d = dimensions(/energy, _EXTRA=extra)
   endif else if $
     (strcmp("electron thermal energy", scalarname, /fold_case) eq 1) or $
     (strcmp("pe", scalarname, /fold_case) eq 1)then begin
      if(version ge 20) then begin
         data = s.E_PE._data 
      endif else begin
         print, 'Error, this data is not present in this version of M3D-C1.'
         data = 0.
      end
       title = 'Electron Thermal Energy'
       symbol = '!8W!de!N!X'
       d = dimensions(/energy, _EXTRA=extra)
    endif else if $
       (strcmp("ion thermal energy", scalarname, /fold_case) eq 1) or $
       (strcmp("pi", scalarname, /fold_case) eq 1)then begin
       if(version ge 20) then begin
          data = s.E_P._data - s.E_PE._data 
       endif else begin
          print, 'Error, this data is not present in this version of M3D-C1.'
          data = 0.
       end
       title = 'Ion Thermal Energy'
       symbol = '!8W!di!N!X'
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
     (strcmp("electrons", scalarname, /fold_case) eq 1) or $
     (strcmp("ne", scalarname, /fold_case) eq 1) then begin
      if(version ge 20) then begin
         data = s.electron_number._data
      endif else begin
         zeff = read_parameter('zeff', filename=filename)
         data = s.particle_number._data*zeff
      end
      title = 'Electron Number'
      symbol = '!8N!De!N!X'
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
       if(itor eq 1) then data = data*rzero
       title = 'Flux'
       symbol = '!7W!X'
       d = dimensions(/b0,l0=2, _EXTRA=extra)
   endif else if (strcmp("li", scalarname, /fold_case) eq 1) then begin
      rzero = read_parameter('rzero', filename=filename)
       data = -4.*!pi*(s.psi_lcfs._data - s.psimin._data) $
              / s.toroidal_current_p._data / rzero
       if(itor eq 0 and version lt 36) then begin
          print, 'WARNING: correcting for incorrect plasma current definition with itor=0 and version<36'
          data = data*rzero
       end
       title = 'Normalized Internal Inductance'
       symbol = '!13l!Di!X'
       d = dimensions(_EXTRA=extra)
   endif else if (strcmp("li3", scalarname, /fold_case) eq 1) then begin
      rzero = read_parameter('rzero', filename=filename)
       Wm = s.w_m._data
       data = 4.*Wm / s.toroidal_current_p._data^2 / rzero
       if(itor eq 0 and version lt 36) then begin
          print, 'WARNING: correcting for incorrect plasma current definition with itor=0 and version<36'
          data = data*rzero^2
       end
       title = 'Normalized Internal Inductance'
       symbol = '!13l!Di!N!6(3)!X'
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
       title = '!6Runaway Current!6'
       symbol = '!8N!D!6RE!N!X'
       d = dimensions(/j0,_EXTRA=extra)
   endif else if (strcmp("IZ", scalarname, /fold_case) eq 1 $
                 or strcmp("M_IZ", scalarname, /fold_case) eq 1) then begin
       data = s.M_IZ._data / s.toroidal_current_p._data
       title = '!6Plasma Current Centroid!6'
       symbol = '!8Z!DI!N!X'
       d = dimensions(/l0,_EXTRA=extra)
   endif else if (strcmp("radiation", scalarname, /fold_case) eq 1 $
                  or strcmp("prad", scalarname, /fold_case) eq 1) then begin
       data = -s.radiation._data
       title = '!6Radiated Power!6'
       symbol = '!8P!D!6rad!N!X'
       d = dimensions(/p0,l0=3,t0=-1,_EXTRA=extra)
   endif else if (strcmp("line_rad", scalarname, /fold_case) eq 1 $
       or strcmp("pline", scalarname, /fold_case) eq 1) then begin
       data = -s.line_rad._data
       title = '!6Line Radiation Power!6'
       symbol = '!8P!D!6rad!N!X'
       d = dimensions(/p0,l0=3,t0=-1,_EXTRA=extra)
   endif else if (strcmp("brem_rad", scalarname, /fold_case) eq 1 $
       or strcmp("pbrem", scalarname, /fold_case) eq 1) then begin
       data = -s.brem_rad._data
       title = '!6Bremsstrahlung Radiation Power!6'
       symbol = '!8P!D!6rad!N!X'
       d = dimensions(/p0,l0=3,t0=-1,_EXTRA=extra)
   endif else if (strcmp("ion_loss", scalarname, /fold_case) eq 1 $
       or strcmp("pion", scalarname, /fold_case) eq 1) then begin
       data = -s.ion_loss._data
       title = '!6Ionization Power!6'
       symbol = '!8P!D!6rad!N!X'
       d = dimensions(/p0,l0=3,t0=-1,_EXTRA=extra)
   endif else if (strcmp("reck_rad", scalarname, /fold_case) eq 1 $
       or strcmp("preck", scalarname, /fold_case) eq 1) then begin
       data = -s.reck_rad._data
       title = '!6Recombination Radiation Power (Kinetic)!6'
       symbol = '!8P!D!6rad!N!X'
       d = dimensions(/p0,l0=3,t0=-1,_EXTRA=extra)
   endif else if (strcmp("recp_rad", scalarname, /fold_case) eq 1 $
       or strcmp("precp", scalarname, /fold_case) eq 1) then begin
       data = -s.recp_rad._data
       title = '!6Recombination Radiation Power (Potential)!6'
       symbol = '!8P!D!6rad!N!X'
       d = dimensions(/p0,l0=3,t0=-1,_EXTRA=extra)
   endif else if (strcmp("temax", scalarname, /fold_case) eq 1) then begin
       data = s.temax._data
       title = '!6Maximum Te!X'
       symbol = '!6max[!8T!De!N!6]!X'
       d = dimensions(/temp,_EXTRA=extra)
   endif else if (strcmp("POhm", scalarname, /fold_case) eq 1) then begin
       data = -(s.e_mpd._data + s.e_mtd._data)
       title = '!6Ohmic Heating!X'
       symbol = '!8P!D!6ohm!N!X'
       d = dimensions(/p0,t0=-1,l0=3,_EXTRA=extra)
   endif else if (strcmp("ave_p", scalarname, /fold_case) eq 1) then begin
       data = s.Ave_p._data
       title = '!6Average Pressure!X'
       symbol = '!3<!8p!3>!X'
       d = dimensions(/p0,_EXTRA=extra)
   endif else if (strcmp("helicity", scalarname, /fold_case) eq 1) then begin
       data = s.helicity._data
       title = '!6Magnetic Helicity!X'
       symbol = '!8H!X'
       d = dimensions(b0=2,l0=4,_EXTRA=extra)
   endif else if (strcmp("pinj", scalarname, /fold_case) eq 1) then begin
       data = s.power_injected._data
       title = '!6Power Injected!X'
       symbol = '!8P!Dinj!N!X'
       d = dimensions(p0=1,t0=-1,l0=3,_EXTRA=extra)
   endif else if (strcmp("kprad_n0", scalarname, /fold_case) eq 1) then begin
       data = s.kprad_n0._data
       title = '!6Neutral Impurities!X'
       symbol = '!8n!DI!60!N!X'
       d = dimensions(/n0,l0=3,_EXTRA=extra)
   endif else if (strcmp("kprad_n", scalarname, /fold_case) eq 1) then begin
       data = s.kprad_n._data
       title = '!6Total Impurities!6'
       symbol = '!8n!DI!N!X'
       d = dimensions(/n0,l0=3,_EXTRA=extra)
   endif else if (strcmp("kprad_ion_frac", scalarname, /fold_case) eq 1) then begin
       data = 1. - s.kprad_n0._data / s.kprad_n._data
       title = '!6Impurity Ionization Fraction!X'
       symbol = title
       d = dimensions(_EXTRA=extra)
   endif else if (strcmp("zeff", scalarname, /fold_case) eq 1) then begin
       data = s.electron_number._data / $
             (s.particle_number._data + s.kprad_n._data)
       title = '!6Free electrons per nucleus!X'
       symbol = '!8Z!Deff!N!X'
       d = dimensions(_EXTRA=extra)

   endif else begin
       s = read_scalars(filename=filename)
       n = tag_names(s)
       smatch = where(strcmp(n, scalarname, /fold_case) eq 1,scount)
       if(version ge 31) then begin
          p = read_pellets(filename=filename)
          if (p eq !NULL) then begin
            print, 'Warning: pellet data not present in this file'
            pcount = 0
          endif else begin
            n = tag_names(p)
            pmatch = where(strcmp(n, scalarname, /fold_case) eq 1,pcount)
          endelse
       endif else begin
          pcount = 0
       endelse

       if(scount ne 0) then begin
           data = s.(smatch[0])._data
       endif else if(pcount ne 0) then begin
           if ipellet eq -1 then begin
              data = p.(pmatch[0])._data
           endif else begin
              data = p.(pmatch[0])._data[ipellet,*]
           endelse
       endif else begin
           print, 'Scalar ', scalarname, ' not recognized.'
           return, 0
       endelse

       title = ''
       symbol = scalarname
       d = dimensions()
       
       if(strcmp("wall_force", scalarname, 10, /fold_case) eq 1) then begin
          d = dimensions(/p0,l0=2)
          if((version lt 27) and (threed eq 0)) then begin
             data = data*2.*!pi
          end
       endif
       if(strcmp("flux_thermal", scalarname, 10, /fold_case) eq 1) then begin
          symbol = '!7C!D!8t!N!X'
          d = dimensions(/p0,l0=3,t0=-1)
       endif
   endelse
   
   if(keyword_set(final)) then begin
     if ipellet eq -1 then begin
       data = data[*,n_elements(data)-1]
     endif else begin
       data = data[n_elements(data)-1]
     endelse
   endif

   if(keyword_set(integrate)) then begin
       d = d + dimensions(/t0)
       if ipellet eq -1 then begin
         for ip=0,n_elements(data[*,0])-1 do data[ip,*] = cumtrapz(time,data[ip,*])
       endif else begin
         data = cumtrapz(time,data)
       endelse
       title = 'Integrated '+title
       symbol = '!Mi '+symbol+'!6 dt!6'
   endif

   get_normalizations, b0=b0,n0=n0,l0=l0, ion_mass=mi, $
     filename=filename, _EXTRA=extra
   convert_units, data, d, b0, n0, l0, mi, _EXTRA=extra
   convert_units, time, dimensions(/t0), b0, n0, l0, mi, _EXTRA=extra
   units = parse_units(d, _EXTRA=extra)

   if ipellet ne -1 then begin
     if(n_elements(data) gt n_elements(time)) then $
        data = data[0:n_elements(time)-1]
     if(n_elements(time) gt n_elements(data)) then $
        time = time[0:n_elements(data)-1]
   endif else begin
     N = size(data)
     N = N[2]
     if(N gt n_elements(time)) then $
       data = data[*,0:n_elements(time)-1]
     if(n_elements(time) gt N) then $
       time = time[0:N-1]
   endelse
   return, reform(data)
end

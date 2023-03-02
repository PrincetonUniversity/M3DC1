;==============================================================
; energy
; ~~~~~~
;
; Returns a time series of the total energy within the
; simulation domain
; 
; Optional outputs:
;  t: the time array
;  names: the name of each energy component
;  components: an array of time series of each energy component
;==============================================================
function energy, filename=filename, components=comp, names=names, t=t, $
                 no_ke=no_ke, xtitle=xtitle, ytitle=ytitle, _EXTRA=extra

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = read_parameter("numvar", filename=filename)
   version = read_parameter("version", filename=filename)
   imr = read_parameter("imulti_region", filename=filename)
   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   t = s.time._data

   names = ['Solenoidal KE', 'Solenoidal ME', $
            'Toroidal KE', 'Toroidal ME', $
            'Compressional KE', 'Thermal Pressure']
   ke = [1,0,1,0,1,0]

   if(version ge 43 and imr gt 0) then begin
      comp = fltarr(10,n_elements(t))
      names = [names, ['Solenoidal ME in Conductors', $
                       'Solenoidal ME in Vacuum', $
                       'Toroidal ME in Conductors', $
                       'Toroidal ME in Vacuum']]
      ke = [ke, 0, 0, 0, 0]
   endif else begin
      comp = fltarr(6,n_elements(t))
   end

   comp[0,*] = s.E_KP._data
   comp[1,*] = s.E_MP._data
   comp[2,*] = s.E_KT._data
   comp[3,*] = s.E_MT._data
   comp[4,*] = s.E_K3._data
   comp[5,*] = s.E_P._data
   if(version ge 43 and imr gt 0) then begin
      comp[6,*] = s.E_MPC._data
      comp[7,*] = s.E_MPV._data
      comp[8,*] = s.E_MTC._data
      comp[9,*] = s.E_MTV._data
   end

   ; convert units
   get_normalizations, b0=b0,n0=n0,l0=l0, ion_mass=mi, $
     filename=filename, _EXTRA=extra
   de = dimensions(/p0, l0=3, _EXTRA=extra)
   dt = dimensions(/t0, _EXTRA=extra)
   convert_units, comp, de, b0, n0, l0, mi, _EXTRA=extra
   convert_units, t, dt, b0, n0, l0, mi, _EXTRA=extra
   units = parse_units(de, _EXTRA=extra)
   xtitle='!8t!6 ('+parse_units(dt, _EXTRA=extra)+')!X'
   ytitle='!8E!6 ('+parse_units(de, _EXTRA=extra)+')!X'

   toten = total(comp, 1)

   if(keyword_set(no_ke)) then begin
      i = where(ke eq 0)
      comp = reform(comp[i,*])
      names = names[i]
   end

   return, toten
end

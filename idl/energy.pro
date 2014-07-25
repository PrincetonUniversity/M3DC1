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
function energy, filename=filename, components=comp, names=names, t=t

   if(n_elements(filename) eq 0) then filename='C1.h5'

   nv = read_parameter("numvar", filename=filename)
   s = read_scalars(filename=filename)
   if(n_tags(s) eq 0) then return, 0

   t = s.time._data
   comp = fltarr(nv*2,n_elements(t))

   names = ['Solenoidal KE', 'Solenoidal ME']

   comp[0,*] = s.E_KP._data
   comp[1,*] = s.E_MP._data
   if(nv ge 2) then begin
       names = [names, ['Toroidal KE', 'Toroidal ME']]
       comp[2,*] = s.E_KT._data 
       comp[3,*] = s.E_MT._data
   endif
   if(nv ge 3) then begin
       names = [names, ['Compressional KE', 'Thermal Pressure']]
       comp[4,*] = s.E_K3._data
       comp[5,*] = s.E_P._data
   endif

   return, total(comp, 1)
end

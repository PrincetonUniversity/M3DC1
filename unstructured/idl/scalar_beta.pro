; ==============
; beta = 2*p/B^2
; ==============
function scalar_beta, filename=filename
;   nv = read_parameter("numvar", filename=filename)
;   if(nv lt 3) then begin
;       print, "Must be numvar = 3 for beta calculation"
;       return, 0
;   endif

   gamma = read_parameter('gam', filename=filename)
   s = read_scalars(filename=filename)

   return, (gamma - 1.)*s.E_P._data / (s.E_MP._data + s.E_MT._data)
end

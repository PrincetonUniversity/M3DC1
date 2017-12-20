; ==========================
; beta_poloidal = 2*P/(Ip^2)
; ==========================
function scalar_beta_poloidal, filename=filename
   nv = read_parameter("numvar", filename=filename)
   if(nv lt 3) then begin
       print, "Must be numvar = 3 for beta calculation"
       return, 0
   endif

   gamma = read_parameter('gam', filename=filename)
   s = read_scalars(filename=filename)

   return, 2.*(gamma-1.)*s.E_P._data/s.toroidal_current._data^2
end

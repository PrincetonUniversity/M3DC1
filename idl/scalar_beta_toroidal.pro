; ======================================
; beta_toroidal = 2*<P>/B_T0^2
; ======================================
function scalar_beta_toroidal, filename=filename

   gamma = read_parameter('gam', filename=filename)
   xmag = read_parameter('xmag', filename=filename)
   rzero = read_parameter('rzero', filename=filename)
   if(rzero eq 0) then rzero = xmag
   bzero = read_parameter('bzero', filename=filename)

   bt0 = bzero*(rzero/xmag)
   
   s = read_scalars(filename=filename)

   beta_t = 2.*s.Ave_P._data/bt0^2
   
   print, 'bt0 =', bt0
   print, 'rzero = ', rzero
   print, 'bzero = ', bzero

   return, beta_t
end

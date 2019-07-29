; ======================================
; beta_toroidal = 2*<P>/B_T0^2
; ======================================
function scalar_beta_toroidal, filename=filename

   s = read_scalars(filename=filename)

   bzero = read_parameter('bzero', filename=filename)
   rzero = read_parameter('rzero', filename=filename)

   shape = get_shape(filename=filename, slice=-1)
   a = shape.a
   r0 = shape.r0
   bt0 = abs(bzero*rzero/r0)

   print, 'R0 = ', rzero
   print, 'BT (T) = ', bt0

   beta_t = 2.*s.Ave_P._data/bt0^2
   
   return, beta_t
end

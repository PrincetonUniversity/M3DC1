; ====================================
; beta_normal = beta_t * (B_T*a / I_p)
; ====================================
function scalar_beta_normal, filename=filename
  
   s = read_scalars(filename=filename)
   ip = read_scalar('ip', filename=filename, /mks)/1.e6

   bzero = read_parameter('bzero', filename=filename)
   rzero = read_parameter('rzero', filename=filename)

   shape = get_shape(filename=filename, slice=-1)
   a = shape.a
   r0 = shape.r0
   bt0 = abs(bzero*rzero/r0)

   print, 'R0 = ', rzero
   print, 'IP (MA) = ', ip
   print, 'BT (T) = ', bt0

   beta_t = 2.*s.Ave_P._data/bt0^2
   beta_n = beta_t * abs(bt0*a/ip)
   
   return, 100.*beta_n
end

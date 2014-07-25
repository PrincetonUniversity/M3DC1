; ====================================
; beta_normal = beta_t * (B_T*a / I_p)
; ====================================
function scalar_beta_normal, filename=filename

   a = 1.

   gamma = read_parameter('gam', filename=filename)
   xmag = read_parameter('xmag', filename=filename)
   rzero = read_parameter('rzero', filename=filename)
   if(rzero eq 0) then rzero = xmag
   bzero = read_parameter('bzero', filename=filename)

   bt0 = bzero*(rzero/xmag)
   
   s = read_scalars(filename=filename)
   ip = s.toroidal_current._data

   beta_t = 2.*s.Ave_P._data/bt0^2
   beta_n = beta_t * abs(bt0*a/ip)

   print, 'ip', ip
   
   return, 100.*beta_n
end

function analytic_model, r0, a, d, z0, b, points=n
   if(n_elements(n) eq 0) then n=100

   xy = fltarr(2,n)

   t = 2.*!pi*findgen(n)/n

   xy[0,*] = r0 + a*cos(t + d*sin(t))
   xy[1,*] = z0 + b*sin(t)

   return, xy
end

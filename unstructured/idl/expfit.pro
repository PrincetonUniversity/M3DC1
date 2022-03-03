function exp_func, x, a
  return, [ [a[0]*(1.-exp(-a[1]*x))], $
            [     (1.-exp(-a[1]*x))], $
            [a[0]*x*exp(-a[1]*x) ]]
end

pro expfit, x, y, plot=pl, a=a0, t=t0

  n = n_elements(x)
  y0 = y[n-1]
  
  dxdy = deriv(x, y/y0)
  
  tau = 1./dxdy[0]

  a = [y0, 1./tau]
  if(n_elements(a0) eq 1) then a[0] = a0
  if(n_elements(t0) eq 1) then a[1] = t0
  v = lmfit(x, y, a, func='exp_func')

  print, "Fit = A*[1-Exp(-x/T)]"
  print, "  Guess A = ", y0
  print, "  Guess T = ", tau

  print, "  Fit A = ", a[0]
  print, "  Fit T = ", 1./a[1]
 
; plot, x, -y/alog(1-y/y0)
 
  if(keyword_set(pl)) then begin
     plot, x, y, psym=4
     n = 100
     xx = findgen(n) * (max(x) - min(x))/(n-1) + min(x)
;     yy = y0*(1. - exp(-xx/tau))
     yy = exp_func(xx, a)
     oplot, xx, yy
  end
end

function exp_func, x, a
  return, [ [ a[0]*(1.-exp(-a[1]*(x-a[2])))], $
            [      (1.-exp(-a[1]*(x-a[2])))], $
            [ a[0]*(x-a[2])*exp(-a[1]*(x-a[2])) ], $
            [-a[0]*   a[1] *exp( a[1]*(x-a[2])) ]]
end

pro expfit, xin, yin, plot=pl, a=a0, t=t0, start=start

  if(keyword_set(start)) then begin
     i = where(xin gt start)
     x = xin[i]
     y = yin[i]
  endif else begin
     x = xin
     y = yin
  end
  

  n = n_elements(x)
  y0 = y[n-1]
  
  dxdy = deriv(x, y/y0)
  
  tau = 1./mean(dxdy[0:n/4])

  a = [y0, 1./tau, 0.]
  if(n_elements(a0) eq 1) then a[0] = a0
  if(n_elements(t0) eq 1) then a[1] = t0
  v = lmfit(x, y, a, func='exp_func')

  print, "Fit = A*[1-Exp(-(x-t0)/T)]"
  print, "  Guess A = ", y0
  print, "  Guess T = ", tau
  print, "  Guess t0 = ", a[2]

  print, "  Fit A = ", a[0]
  print, "  Fit T = ", 1./a[1]
  print, "  Fit t0 = ", a[2]
 
; plot, x, -y/alog(1-y/y0)
 
  if(keyword_set(pl)) then begin
     plot, x, y, psym=4, xrange=[xin[0], xin[n-1]]
     n = 100
     xx = findgen(n) * (max(x) - min(x))/(n-1) + min(x)
;     yy = y0*(1. - exp(-xx/tau))
     yy = exp_func(xx, a)
     oplot, xx, yy
  end
end

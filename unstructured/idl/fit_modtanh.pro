pro mod_tanhfit, x, a, f, pder
   t = tanh(a[1]*(x-a[0]))
   s = 0.5*a[2]*(1.+a[3]*(1.-x)+a[4]*(1.-x)^2+a[5]*(1.-x)^3)
   f = s*(1.-t)

   if(n_params() ge 4) then begin
       pder = [[-s*(1.-t^2)*(-a[1])],$
               [-s*(1.-t^2)*(x-a[0])],$
               [0.5*(1.+a[3]*(1.-x)+a[4]*(1.-x)^2+a[5]*(1.-x)^3)*(1.-t)],$
               [0.5*a[2]*(1.-x)*(1.-t)], $
               [0.5*a[2]*(1.-x)^2*(1.-t)], $
               [0.5*a[2]*(1.-x)^3*(1.-t)]]
   end
end

pro mod_tanhfit2, x, a, f, pder
   t = tanh(a[1]*(x-a[0]))
   s = 0.5*a[2]*(1.+a[3]*(1.-x)+a[4]*(1.-x)^2+a[5]*(1.-x)^3)
   f = s*(1.-t) + a[6]

   if(n_params() ge 4) then begin
       pder = [[-s*(1.-t^2)*(-a[1])],$
               [-s*(1.-t^2)*(x-a[0])],$
               [0.5*(1.+a[3]*(1.-x)+a[4]*(1.-x)^2+a[5]*(1.-x)^3)*(1.-t)],$
               [0.5*a[2]*(1.-x)*(1.-t)],$
               [0.5*a[2]*(1.-x)^2*(1.-t)],$
               [0.5*a[2]*(1.-x)^3*(1.-t)],$
               [replicate(1.,n_elements(x))]]
   end
end

pro mod_dtanhfit, x, a, f, pder
  s = 1./cosh(a[1]*(x-a[0]))^2
  t = tanh(a[1]*(x-a[0]))
  p = -0.5*a[1]*a[2]*(1.+a[3]*(1.-x)+   a[4]*(1.-x)^2+   a[5]*(1.-x)^3)
  q = -0.5     *a[2]*   (a[3]       +2.*a[4]*(1.-x)  +3.*a[5]*(1.-x)^2)
  f = p*s + q*(1.-t)
  if(n_params() ge 4) then begin
     pder = [[a[1]*s*(q + 2.*p*t)], $
             [s*(p/a[1]-(q+2.*p*t)*(x-a[0]))], $
             [f/a[2]], $
             [-0.5*a[2]         *(a[1]*(1.-x)*s+   (1.-t))], $
             [-0.5*a[2]*(1.-x)  *(a[1]*(1.-x)*s+2.*(1.-t))], $
             [-0.5*a[2]*(1.-x)^2*(a[1]*(1.-x)*s+3.*(1.-t))]]
  end
end

pro fit_modtanh, x, y, xout, yout, parameters=a, weights=w, minval=minval, $
                 plot=plot, deriv=d, sigma=sig, chisq=chi, itmax=itmax

  n = n_elements(x)

  if(keyword_set(d)) then begin
     if(n_elements(w) eq 0) then w = x * (x lt 1)
;     if(n_elements(w) eq 0) then w = replicate(1., n) * (x lt 1)
     if(n_elements(a) eq 0) then $
        a = [0.99, 50., max(y), 0., 0., 0., 0.]
     a5 = a[0:5]
     yfit = curvefit(x, y, w, a5, sig, chisq=chi, $
                     function_name='mod_dtanhfit', itmax=itmax)
     a[0:5] = a5
     mod_dtanhfit, xout, a, yout
  endif else begin
     if(n_elements(w) eq 0) then w = 1/y
     if(n_elements(minval) eq 0) then begin
        if(n_elements(a) eq 0) then $
           a = [0.98, 50., y[n*0.96], 0., 0., 0., min(y)]
        yfit = curvefit(x, y, w, a, sig, chisq=chi, $
                        function_name='mod_tanhfit2', itmax=itmax)
        mod_tanhfit2, xout, a, yout
     endif else begin
        if(n_elements(a) eq 0) then $
           a = [0.99, 50., y[n*0.96], 0., 0., 0., 0.]
        a5 = a[0:5]
        yfit = curvefit(x, y-minval, w, a5, sig, chisq=chi, $
                        function_name='mod_tanhfit', itmax=itmax)
        a[0:5] = a5
        mod_tanhfit, xout, a, yout
        yout = yout + minval
     endelse
  endelse
   
  if(keyword_set(plot)) then begin
     plot, x, y, xrange=[0,1.1]
     oplot, xout, yout, color=color(1)
  end

end

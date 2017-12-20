pro plot_bmn_vs_m, filename, mrange=mrange, ylog=ylog, factor=factor, $
                 psin=psi0, reverse_q=reverse_q, _EXTRA=extra

  if(n_elements(psi0) eq 0) then psi0 = 0.95
  if(n_elements(factor) eq 0) then factor = 1.

  read_bmncdf, file=filename, _EXTRA=extra, bmn=bmn, psi=psi, m=m, q=q, $
               ntor=ntor

  modbmn = abs(bmn)*factor

  q0 = interpol(q,psi,psi0)
  if(n_elements(mrange) eq 0) then mrange=fix([-3,3]*abs(q0*ntor))

  if(keyword_set(reverse_q)) then q0 = -q0

  n = mrange[1]-mrange[0]+1
  mm = indgen(n) + mrange[0]

  ct3
  f = fltarr(n)
  for j=0,n-1 do begin
     m0 = mm[j]
     i = where(m eq m0, count)
     if(count ne 1) then begin
        print, 'Error evaluating m = ', m0
        continue
     end
     f[j] = interpol(modbmn[i,*],psi,psi0)
  end
  plot, mm, f, _EXTRA=extra, xtitle='!8m!X'
  oplot, [q0,q0]*ntor, !y.crange, linestyle=2

  names = string(format='("!7W!3 = ",G0.3,"!X")', psi0)
  plot_legend, names
end

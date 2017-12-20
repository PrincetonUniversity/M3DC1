pro plot_m_vs_r, filename, mrange=mrange, ylog=ylog, factor=factor, $
                 srnorm=srnorm, rhonorm=rhonorm, _EXTRA=extra, phase=phase
  if(n_elements(factor) eq 0) then factor = 1.

  read_bmncdf, file=filename, _EXTRA=extra, bmn=bmn, psi=psi, m=m, q=q, $
               rho=rho, ntor=ntor, symbol=symbol, units=units

  label = symbol + '!6 (' + units + '!6)!X'

  if(n_elements(mrange) eq 0) then mrange=[-5,5]
  if(mrange[1] lt mrange[0]) then mrange=reverse(mrange)

  if(keyword_set(srnorm)) then begin
     xtitle = '!9r!7W!X'
     psi = sqrt(psi)
  endif else if(keyword_set(rhonorm)) then begin
     xtitle = '!7q!X'
     psi = rho
  endif else xtitle='!7W!X'

  if(keyword_set(phase)) then begin
     data = atan(imaginary(bmn), real_part(bmn))*180./!pi
     yran = [-180.,180.]
     ytitle='!6Phase!X'
  endif else begin
     data = abs(bmn)*factor
     yran = [0, max(data)]
     ytitle=label
  endelse
  plot, [0,1], yran, /nodata, _EXTRA=extra, $
        xtitle=xtitle, ylog=ylog, ytitle=ytitle

  n = mrange[1]-mrange[0]+1
  mm = indgen(n) + mrange[0]
  qq = float(mm)/float(ntor)

  psin = flux_at_q(abs(qq),q=abs(q),flux=psi)

  ct3
  c = intarr(n)
  name = strarr(n)
  for j=0,n-1 do begin
     m0 = mm[j]
     i = where(m eq m0, count)
     if(count ne 1) then begin
        print, 'Error evaluating m = ', m0
        continue
     end
     c[j] = color(abs(m0))
     name[m0-mrange[0]] = string(format='("!8m!3 = ",I0,"!X")',m0) 

     if(keyword_set(ylog)) then begin
        yrange = 10^!y.crange
     endif else yrange = !y.crange
     oplot, psi, data[i,*], color=c[j]
     oplot, [psin[j],psin[j]], yrange, color=c[j], linestyle=2
  end
  plot_legend, name, color=c, ylog=ylog
end

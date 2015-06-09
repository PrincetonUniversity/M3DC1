pro plot_m_vs_r, filename, mrange=mrange, ylog=ylog, factor=factor, $
                 srnorm=srnorm, _EXTRA=extra
  if(n_elements(factor) eq 0) then factor = 1.

  read_bmncdf, file=filename, _EXTRA=extra, bmn=bmn, psi=psi, m=m, q=q, $
               ntor=ntor

  if(n_elements(mrange) eq 0) then mrange=[-5,5]

  if(keyword_set(srnorm)) then begin
     xtitle = '!9r!7W!X'
     psi = sqrt(psi)
  endif else xtitle='!7W!X'

  modbmn = abs(bmn)*factor
  print, max(modbmn)
  plot, [0,1], [0,max(modbmn)], /nodata, _EXTRA=extra, $
        xtitle=xtitle, ylog=ylog, $
        ytitle='!8B!Dmn!N!6 (G/kA)!X'

  n = mrange[1]-mrange[0]+1
  mm = indgen(n) + mrange[0]
  qq = abs(float(mm)/float(ntor))

  psin = flux_at_q(qq,q=q,flux=psi)

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
     oplot, psi, modbmn[i,*], color=c[j]
     oplot, [psin[j],psin[j]], yrange, color=c[j], linestyle=2
  end
  plot_legend, name, color=c, ylog=ylog
end

pro plot_m_vs_r, filename, mrange=mrange, ylog=ylog, factor=factor, $
                 srnorm=srnorm, rhonorm=rhonorm, _EXTRA=extra, phase=phase, $
                 jmn=jmn, out=outfile, overplot=over, nolegend=nolegend, $
                 linestyle=linestyle, current=cur, ntor=ntor0, title=title

  if(n_elements(factor) eq 0) then factor = 1.

  read_bmncdf, file=filename, _EXTRA=extra, bmn=bmn0, psi=psi, m=m, q=q, $
               rho=rho, ntor=ntor, symbol=symbol, units=units, jmn=j, $
               current=cur

  if(n_elements(ntor) gt 1) then begin
     if(n_elements(ntor0) eq 1) then begin
        i = where(ntor eq ntor0, ct)
        if(ct ne 1) then begin
           print, 'Error in plot_m_vs_r: ntor = ', ntor0, ' not found'
           return
        end
     endif else begin
        print, 'Error in plot_m_vs_r: must select ntor'
        return
     endelse
  endif else begin
     ntor0 = ntor
     i = 0
  end

  if(n_elements(title) eq 0) then begin
     title = '!8n!6 = ' + string(format='(I0)',ntor0) + '!X'
  end

  sz = size(bmn0)
  if(sz[0] eq 3) then begin
     bmn = reform(bmn0[i,*,*])
  endif else bmn = bmn0

  if(keyword_set(jmn)) then begin
     bmn = j
     symbol = '!8J!Dmn!N!X'
     units = 'A/m!U2!N'
  end

  label = symbol + '!6 (' + units + '!6)!X'

  if(n_elements(mrange) eq 0) then mrange=[-5,5]
  if(mrange[1] lt mrange[0]) then mrange=reverse(mrange)

  if(keyword_set(srnorm)) then begin
     xtitle = '!9r!7W!X'
     psi = sqrt(psi)
     xlabel = 'Sqrt(Psi_N)'
  endif else if(keyword_set(rhonorm)) then begin
     xtitle = '!7q!X'
     psi = rho
     xlabel = 'Rho_N'
  endif else begin
     xtitle='!7W!X'
     xlabel = 'Psi_N'
  end

  if(keyword_set(phase)) then begin
     data = atan(imaginary(bmn), real_part(bmn))*180./!pi
     yran = [-180.,180.]
     ytitle='!6Phase!X'
  endif else begin
     data = abs(bmn)*factor
     yran = [0, max(data)]
     ytitle=label
  endelse
  if(not keyword_set(over)) then begin
     plot, [0,1], yran, /nodata, _EXTRA=extra, $
           xtitle=xtitle, ylog=ylog, ytitle=ytitle, title=title
  end

  n = mrange[1]-mrange[0]+1
  mm = indgen(n) + mrange[0]
  qq = float(mm)/float(ntor0)

  psin = flux_at_q(abs(qq),abs(q),flux=psi,q=qout)

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
     oplot, psi, data[i,*], color=c[j], linestyle=linestyle
     
     if(psin[0] eq 0) then continue
     for k=0, n_elements(psin)-1 do begin
        if(abs(abs(qout[k])-abs(qq[j])) gt 0.001) then continue
        oplot, [psin[k],psin[k]], yrange, color=c[j], linestyle=2
     end
  end
  if(not keyword_set(nolegend)) then begin
     plot_legend, name, color=c, ylog=ylog, _EXTRA=extra
  end

  if(n_elements(outfile) eq 1) then begin
     openw, ifile, outfile, /get_lun
     printf, ifile, format='(A12,31I12)', xlabel, mm
     k = intarr(n)
     for i=0, n-1 do begin
        k[i] = where(m eq mm[i], count)
     end
     for i=0, n_elements(psi)-1 do begin
        printf, ifile, format='(32F12.4)', psi[i], data[k,i]
     end
     free_lun, ifile
  end

end

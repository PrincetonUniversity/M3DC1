pro plot_bmn, filename, file=fn, vac=vac, names=names, nolegend=nolegend, $
              ytitle=ytitle, width=width, current=cur, $
              netcdf=netcdf, chirikov=chi, sum_files=sum_files, $
              color=c, overplot=overplot, monochrome=bw, $
              linestyle=linestyle, out=out, phase=ph, $
              _EXTRA=extra

   if(n_elements(filename) eq 0 and n_elements(fn) gt 0) then filename=fn

   if(n_elements(cur) eq 0) then cur=1.
   if(n_elements(cur) lt n_elements(filename)) then $
     cur = replicate(cur, n_elements(filename))

   xtitle = '!7W!X'

   phase = 0
   if(keyword_set(width)) then begin
       data = island_widths(filename, psin=psin, cur=cur, sum_files=sum_files, $
                          netcdf=netcdf, q=q)
       if(n_elements(ytitle) eq 0) then $
         ytitle='!6Island Width (!7W!6)!X'
   endif else if(keyword_set(chi)) then begin
       data = chirikov(filename,cur=cur,psimid=psin, sum_files=sum_files, $
                     netcdf=netcdf, qval=q)
       if(n_elements(ytitle) eq 0) then $
         ytitle='!6Chirikov Parameter!X'
   endif else begin
      result = read_bmn(filename,m,bmn,phase,psin=psin,$
                        sum_files=sum_files,factor=cur, $
                        netcdf=netcdf, qval=q)

      if(keyword_set(ph)) then begin
         data = phase*180/!pi
         if(n_elements(ytitle) eq 0) then $
            ytitle='!6Phase (deg)!X'
      endif else begin
         data = bmn
         if(n_elements(ytitle) eq 0) then $
            ytitle='!6Total Resonant Field (G/kA)!X'
      endelse
   endelse

   ct3
   if(n_elements(c) eq 0) then begin
       if(n_elements(filename) gt 1) then begin
          c = shift(get_colors(n_elements(filename)),-1)
       endif else begin
          c = get_colors()
       endelse
    end
   if(n_elements(linestyle) eq 0) then linestyle = 0
   if(n_elements(linestyle) lt n_elements(filename)) then begin
      linestyle = replicate(linestyle, n_elements(filename))
   end
   sym = [4, 5, 6, 7]
   if(not keyword_set(overplot)) then begin
       plot, psin, data, /nodata, $
         xtitle=xtitle, ytitle=ytitle, _EXTRA=extra
   end
   for i=0, n_elements(data[*,0])-1 do begin
       oplot, psin[i,*], data[i,*], color=c[i], linestyle=linestyle[i]
       oplot, psin[i,*], data[i,*], color=c[i], psym=sym[i mod 4]
    end

   if(n_elements(out) gt 0) then begin
      if(n_elements(out) eq 1 and n_elements(data[*,0]) gt 1) then begin
         out = replicate(out, n_elements(data[*,0]))
         out = out + string(indgen(n_elements(data[*,0]))+1,format='(I0)')
      end
      for i=0, n_elements(data[*,0])-1 do begin
         openw, ifile, out[i], /get_lun
         if(n_elements(phase) gt 1) then begin
            for j=1, n_elements(data[0,*])-1 do begin
               printf, ifile, psin[i,j], q[i,j], data[i,j], phase[i,j]*180./!pi
            end
         endif else begin
            for j=1, n_elements(data[0,*])-1 do begin
               printf, ifile, psin[i,j], q[i,j], data[i,j]
            end
         endelse
         free_lun, ifile    
      end
   end

   if(n_elements(vac) gt 0) then begin
       if(keyword_set(width)) then begin
           data_vac = island_widths(vac, psin=psin_vac, cur=cur, $
                                  sum_files=sum_files, netcdf=netcdf)
       endif else if(keyword_set(chi)) then begin
           data_vac = chirikov(vac,cur=cur,psimid=psin_vac, netcdf=netcdf)
       endif else begin
          result = read_bmn(vac, m_vac, bmn_vac, phase_vac, psin=psin_vac, $
                            sum_files=sum_files, factor=cur, netcdf=netcdf)
          if(keyword_set(ph)) then begin
             data_vac = phase_vac*180./!pi
          endif else begin
             data_vac = bmn_vac
          end
       endelse
       if(n_elements(data_vac[*,0]) eq 1) then begin
           oplot, psin_vac[0,*], data_vac[0,*], linestyle=1
       endif else begin
           for i=0, n_elements(data_vac[*,0])-1 do begin
               oplot, psin_vac[i,*], data_vac[i,*], linestyle=1, color=c[i]
;           oplot, psin_vac[i,*], data_vac[i,*], linestyle=1, color=c[i], psym=5
           end
       end
   end

   if(not keyword_set(nolegend)) then begin
      if(n_elements(filename) gt 1 and not keyword_set(sum_files)) then begin
         if(n_elements(names) eq 0) then names = filename
         plot_legend, names, color=c, psym=sym, $
                      _EXTRA=extra, linestyle=replicate(0,n_elements(filename))
      end
   end
end

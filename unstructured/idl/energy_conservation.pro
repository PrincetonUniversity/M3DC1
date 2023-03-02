; WIP: this could be made more featureful and flexible
; Plots the change in various energies over time
; fluxes can be a list of strings corresponding to flux fields that will be integrated and added to lost energy
; For energy conservation, the dashed line should be zero
pro energy_conservation,  filename=filename, poloidal=poloidal, xrange=xrange, yrange=yrange, fluxes=fluxes, comp=comp, overplot=overplot, _EXTRA=ex

  E = energy(filename=filename, components=comp, names=names, t=t, $
             xtitle=xtitle, ytitle=ytitle, _EXTRA=ex)

  Nflux = n_elements(fluxes)
  if (Nflux eq 0) then fluxes = []
  if keyword_set(poloidal) then begin
    fluxes = [fluxes, 'pohm']
    Nflux = Nflux + 1
  endif

  names = [['Total'],names] 

  Norig = n_elements(comp[*,0])
  Ncomp = Norig + Nflux
  Scomp = strarr(Ncomp)
  
  if (Nflux gt 0) then begin
    comp2 = fltarr(Ncomp,n_elements(E))
    comp2[0:Norig-1,*] = comp
    for i=0,Nflux-1 do begin
       comp2[Norig+i,*] = $
          read_scalar(fluxes[i],filename=filename,/int,_EXTRA=ex)
      names = [names, '!Mi '+fluxes[i]+'!X']
    endfor
    names = [names, 'Total + flux']
    comp = comp2
  endif
    
  if (n_elements(xrange) gt 0) then begin
    valid = where((t ge xrange[0]) and(t le xrange[-1]))
    E = E[valid]
    comp2 = fltarr(Ncomp,n_elements(valid))
    for i=0,Ncomp-1 do begin
      comp2[i,*] = comp[i,valid]
    endfor
    comp = comp2
  endif
  
  cols = get_colors(Ncomp+1)
  if (n_elements(cols) gt Ncomp+1) then cols = cols[0:Ncomp]
  
  if (n_elements(yrange) eq 0) then begin
    ylim = max(abs([comp[*,0]-comp[*,-1],E[0]-E[-1]]))
    yrange=[-ylim,ylim]
  endif else begin
    ylim = max(abs(yrange))    
  endelse

  if(not keyword_set(overplot)) then begin
     plot, [t[0],t[-1]], [0,0], thick=1, color=cols[0], $
           xrange=xrange, yrange=yrange, linestyle=1, charsize=2, $
           xtitle=xtitle, ytitle=ytitle, _EXTRA=ex
  end
  if(not keyword_set(poloidal)) then begin
    
    oplot, t, E-E[0], thick=3, color=cols[0]
    for i=0,Ncomp-1 do begin
      oplot, t, comp[i,*]-comp[i,0],thick=3, color=cols[i+1]
    endfor
    if (n_elements(fluxes) ne 0) then begin
      cons = E
      for i=0,Nflux-1 do begin
        cons = cons + comp[Norig+i,*]
      endfor
      oplot, t, cons-cons[0], thick=3, color=cols[0], linestyle=2
    endif
  endif else begin

    oplot, t, comp[1,*]-comp[1,0],thick=3, color=cols[0]
    if max(abs(comp[3,*]-comp[3,0])) gt 1e-3*ylim then begin
      oplot, t, comp[3,*]-comp[3,0], thick=3, color=cols[2]
    endif
    oplot, t, comp[-1,*]-comp[-1,0], thick=3, color=cols[-1]
    cons = comp[1,*]+comp[-1,*]
    oplot, t, cons-cons[0], thick=3, color=cols[0], linestyle=2
    
  endelse
  
  ls = replicate(0,Ncomp+1)
  if (Nflux ne 0) then ls = [ls, 2]
  mid = (Ncomp+1)/2
  plot_legend, names[0:mid], color=cols[0:mid], charsize=2
  plot_legend, names[mid+1:*], color=[cols[mid+1:*], cols[0]], charsize=2,left=0.3,linestyles=ls[mid+1:*]

  return
end

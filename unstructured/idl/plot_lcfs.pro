; ========================================================
; plot_lcfs
; ~~~~~~~~~
;
; plots the last closed flux surface
; ========================================================
pro plot_lcfs, psi, x, z, psival=psival, overplot=over, _EXTRA=extra

  xy = get_lcfs(psi, x, z, psival=psival, _EXTRA=extra)

  loadct, 12
  if(keyword_set(over)) then begin
     oplot, xy[0,*], xy[1,*], thick=!p.thick*1.5, color=color(6,10)
  endif else begin
     plot, xy[0,*], xy[1,*]
  end
end

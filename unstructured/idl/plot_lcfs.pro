; ========================================================
; plot_lcfs
; ~~~~~~~~~
;
; plots the last closed flux surface
; ========================================================
pro plot_lcfs, psi, x, z, psival=psival, _EXTRA=extra

    if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) then begin
        print, 'reading psi, plot_lcfs'
        psi = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
        print, min(psi), max(psi)
    end

    ; if psival not passed, choose limiter value
    if(n_elements(psival) eq 0) then $
      psival = lcfs(psi,x,z,_EXTRA=extra)

    ; plot contour
    loadct, 12
    xy = path_at_flux(psi, x, z, t, psival, /contiguous)

    oplot, xy[0,*], xy[1,*], thick=!p.thick, color=color(6,10)
end

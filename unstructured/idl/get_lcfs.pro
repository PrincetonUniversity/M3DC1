; ========================================================
; get_lcfs
; ~~~~~~~~
;
; plots the last closed flux surface
; ========================================================
function get_lcfs, psi, x, z, psival=psival, axis=axis, _EXTRA=extra

    if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) then begin
        print, 'reading psi, get_lcfs'
        psi = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
    end

    ; if psival not passed, choose limiter value
    if(n_elements(psival) eq 0) then $
      psival = lcfs(psi,x,z,axis=axis,_EXTRA=extra)

    ; plot contour
    xy = path_at_flux(psi, x, z, t, psival, /contiguous, axis=axis)

    return, xy
end

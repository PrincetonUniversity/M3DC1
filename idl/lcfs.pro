; ========================================================
; lcfs
; ~~~~
;
; returns the flux value of the last closed flux surface
; ========================================================
function lcfs, psi, x, z, axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra

    version = read_parameter('version', _EXTRA=extra)
    if(version ge 3) then begin
        psilim = read_lcfs(axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra)
    endif else begin
        psilim = find_lcfs(psi, x, z,axis=axis, xpoint=xpoint, flux0=flux0, $
                           _EXTRA=extra)
    endelse 

    ifixedb = read_parameter('ifixedb', _EXTRA=extra)
    if(ifixedb eq 1) then psilim = 0.

    print, 'LCFS: '
    print, ' Magnetic axis found at ', axis
    print, ' Active x-point at ', xpoint
    print, ' psi_0, psi_s = ', flux0, psilim
    
    return, psilim
end

; ========================================================
; lcfs
; ~~~~
;
; returns the flux value of the last closed flux surface
; ========================================================
function lcfs, psi, x, z, filename=filename, $
               axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra

    if(n_elements(filename) eq 0) then filename='C1.h5'

    version = read_parameter('version', filename=filename[0], _EXTRA=extra)
    if(version[0] ge 3) then begin
        psilim = read_lcfs(axis=axis, xpoint=xpoint, flux0=flux0, $
                           filename=filename[0], _EXTRA=extra)
    endif else begin
        psilim = find_lcfs(psi, x, z,axis=axis, xpoint=xpoint, flux0=flux0, $
                           filename=filename[0], _EXTRA=extra)
    endelse 

    ifixedb = read_parameter('ifixedb', filename=filename[0], _EXTRA=extra)
    if(ifixedb eq 1) then psilim = 0.

    print, 'LCFS: '
    print, ' Magnetic axis found at ', axis
    print, ' Active x-point at ', xpoint
    print, ' psi_0, psi_s = ', flux0, psilim
    
    return, psilim
end

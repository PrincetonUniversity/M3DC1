function save_external_field, filename=filename, time=time, points=p, $
    phi=phi0, rrange=rrange, zrange=zrange, outfile=outfile, $
    r=rgrid, z=zgrid, br=br_ext, bz=bz_ext, psi=psi_ext, _EXTRA=ex

    ; Evaluate BR, BZ, and psi of the external (RMP/REMC) field on an R-Z
    ; grid at a given toroidal location, and save them together with the
    ; R,Z coordinates to an IDL SAVE file.
    ;
    ; Keywords:
    ;   filename - HDF5 file (default 'C1.h5')
    ;   time     - time slice index (default 0)
    ;   points   - grid resolution (default 200)
    ;   phi      - toroidal angle at which to evaluate the field (default 0.)
    ;   rrange   - [Rmin, Rmax] (default: mesh bounding box)
    ;   zrange   - [Zmin, Zmax] (default: mesh bounding box)
    ;   outfile  - name of .sav file to write (default 'external_field.sav')
    ;   r, z, br, bz, psi - output: 2D [points x points] arrays
    ;
    ; Returns 1 on success.

    compile_opt idl2

    if n_elements(filename) eq 0 then filename = 'C1.h5'
    if n_elements(time) eq 0 then time = 0
    if n_elements(p) eq 0 then p = 200
    if n_elements(phi0) eq 0 then phi0 = 0.
    if n_elements(outfile) eq 0 then outfile = 'external_field.sav'

    itor = read_parameter('itor', filename=filename)
    i3d  = read_parameter('i3d', filename=filename)
    if n_elements(i3d) eq 0 then i3d = 0

    if n_elements(rrange) eq 0 or n_elements(zrange) eq 0 then begin
        mesh = read_mesh(filename=filename)
        if n_elements(rrange) eq 0 then rrange = [min(mesh.x), max(mesh.x)]
        if n_elements(zrange) eq 0 then zrange = [min(mesh.z), max(mesh.z)]
    endif

    x = findgen(p)/p*(rrange[1]-rrange[0]) + rrange[0]
    z = findgen(p)/p*(zrange[1]-zrange[0]) + zrange[0]
    t = time

    if itor eq 1 then rad = radius_matrix(x, z, t) else rad = 1.

    ; psi_ext itself (default operation=1), and its R,Z derivatives for BR, BZ
    psi_ext = read_field('psi_ext', x, z, t, slices=time, mesh=mesh, $
                          filename=filename, points=p, $
                          phi=phi0, rrange=rrange, zrange=zrange)

    psi_r = read_field('psi_ext', x, z, t, slices=time, mesh=mesh, $
                        filename=filename, points=p, operation=2, $
                        phi=phi0, rrange=rrange, zrange=zrange)

    psi_z = read_field('psi_ext', x, z, t, slices=time, mesh=mesh, $
                        filename=filename, points=p, operation=3, $
                        phi=phi0, rrange=rrange, zrange=zrange)

    if i3d eq 1 then begin
        fp_r = read_field('fp_ext', x, z, t, slices=time, mesh=mesh, $
                           filename=filename, points=p, operation=2, $
                           phi=phi0, rrange=rrange, zrange=zrange)
        fp_z = read_field('fp_ext', x, z, t, slices=time, mesh=mesh, $
                           filename=filename, points=p, operation=3, $
                           phi=phi0, rrange=rrange, zrange=zrange)
    endif else begin
        fp_r = 0.
        fp_z = 0.
    endelse

    br_ext  = reform(-psi_z/rad - fp_r, p, p)
    bz_ext  = reform( psi_r/rad - fp_z, p, p)
    psi_ext = reform(psi_ext, p, p)

    rgrid = rebin(x, p, p)
    zgrid = rebin(reform(z, 1, p), p, p)

    save, rgrid, zgrid, br_ext, bz_ext, psi_ext, phi0, time, filename=outfile

    print, 'save_external_field: wrote ' + outfile + $
        ' (R, Z, BR, BZ, psi at phi = ' + strtrim(phi0, 2) + ')'

    return, 1

end

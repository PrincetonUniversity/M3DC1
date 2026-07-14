pro calc_external_b_field, filename=filename, time=time, points=p, $
    phi=phi0, rrange=rrange, zrange=zrange, outfile=outfile, $
    component=comp, b_r=b_r, b_phi=b_phi, b_z=b_z, $
    rgrid=rgrid, zgrid=zgrid, csv_outfile=csv_outfile, _EXTRA=ex

    ; Calculate and plot external (RMP) magnetic field components
    ;   bx_ext, by_ext, bz_ext from the *_ext HDF5 datasets.
    ;
    ; Delegates to plot_field for all plotting and file output,
    ; so all of plot_field's keywords work (cutx, cutz, boundary,
    ; lcfs, coils, axis, mesh, q_contours, overplot, _EXTRA, etc.).
    ;
    ; Keywords:
    ;   filename    - HDF5 file (default 'C1.h5')
    ;   time        - time slice index (default 0)
    ;   points      - grid resolution (default 200)
    ;   phi         - toroidal angle
    ;   rrange      - [Rmin, Rmax]
    ;   zrange      - [Zmin, Zmax]
    ;   outfile     - base name for output; appends _bx.dat, _by.dat, _bz.dat
    ;   component   - 'bx', 'by', 'bz', or 'all' (default)
    ;   b_r         - output: 2D radial magnetic field [points x points]
    ;   b_phi       - output: 2D toroidal magnetic field [points x points]
    ;   b_z         - output: 2D vertical magnetic field [points x points]
    ;   rgrid       - output: 2D R-coordinate meshgrid [points x points]
    ;   zgrid       - output: 2D Z-coordinate meshgrid [points x points]
    ;   csv_outfile - base name for CSV output; writes _br.csv, _bphi.csv,
    ;                 _bz.csv, _r.csv, _z.csv
    ;   _EXTRA      - forwarded to plot_field

    compile_opt idl2

    ; --- defaults ---
    if n_elements(filename) eq 0 then filename = 'C1.h5'
    if n_elements(time) eq 0 then time = 0
    if n_elements(p) eq 0 then p = 200
    if n_elements(phi0) eq 0 then phi0 = 0.
    if n_elements(comp) eq 0 then comp = 'all'

    ; --- read geometry parameters ---
    itor = read_parameter('itor', filename=filename)
    i3d  = read_parameter('i3d', filename=filename)
    if n_elements(i3d) eq 0 then i3d = 0

    ; --- generate R-Z grid ---
    if n_elements(rrange) eq 0 or n_elements(zrange) eq 0 then begin
        mesh = read_mesh(filename=filename)
        if n_elements(rrange) eq 0 then begin
            rrange = [min(mesh.x), max(mesh.x)]
        endif
        if n_elements(zrange) eq 0 then begin
            zrange = [min(mesh.z), max(mesh.z)]
        endif
    endif

    nx = p
    ny = p
    x = findgen(nx)/nx*(rrange[1]-rrange[0]) + rrange[0]
    z = findgen(ny)/ny*(zrange[1]-zrange[0]) + zrange[0]

    ; dummy time for read_field positional arg (ignored, slices= used)
    t = time

    ; --- radius for toroidal field ---
    if itor eq 1 then begin
        r = radius_matrix(x, z, t)
    endif else begin
        r = 1.
    endelse

    ; --- read external field datasets ---
    I_ext = read_field('I_ext', x, z, t, slices=time, mesh=mesh, $
                       filename=filename, points=p, phi=phi0, $
                       rrange=rrange, zrange=zrange)

    psi_z = read_field('psi_ext', x, z, t, slices=time, mesh=mesh, $
                       filename=filename, points=p, operation=3, $
                       phi=phi0, rrange=rrange, zrange=zrange)

    psi_r = read_field('psi_ext', x, z, t, slices=time, mesh=mesh, $
                       filename=filename, points=p, operation=2, $
                       phi=phi0, rrange=rrange, zrange=zrange)

    if i3d eq 1 then begin
        fp_r = read_field('fp_ext', x, z, t, slices=time, mesh=mesh, $
                          filename=filename, points=p, operation=2, $
                          phi=phi0, rrange=rrange, zrange=zrange)
        fp_z = read_field('fp_ext', x, z, t, slices=time, mesh=mesh, $
                          filename=filename, points=p, operation=3, $
                          phi=phi0, rrange=rrange, zrange=zrange)
        f_rp = fp_r
        f_zp = fp_z
    endif else begin
        f_rp = 0.
        f_zp = 0.
    endelse

    ; --- compute external field components ---
    bx_ext = -psi_z / r - f_rp
    by_ext =  I_ext / r
    bz_ext =  psi_r / r - f_zp

    ; --- parse component selection ---
    case strlowcase(comp) of
        'bx' : comps = ['bx']
        'by' : comps = ['by']
        'bz' : comps = ['bz']
        'all': comps = ['bx', 'by', 'bz']
        else: begin
            print, "calc_external_b_field: unknown component '" + comp + "'"
            print, "  valid options: 'bx', 'by', 'bz', 'all'"
            return
        end
    endcase

    titles = ['B_R (External)', 'B_!7u!X (External)', 'B_Z (External)']

    ; --- multi-panel for full set ---
    if comp eq 'all' then !p.multi = [0, 1, 3]

    for i = 0, n_elements(comps)-1 do begin

        case comps[i] of
            'bx': data = bx_ext
            'by': data = by_ext
            'bz': data = bz_ext
        endcase

        ; build keyword structure for plot_field
        kw = {title: titles[i], units: '!8T!X', $
              filename: filename, points: p, phi: phi0, $
              rrange: rrange, zrange: zrange}

        if n_elements(outfile) gt 0 then $
            kw = create_struct(kw, 'outfile', outfile + '_' + comps[i] + '.dat')

        if n_elements(ex) gt 0 then $
            kw = create_struct(kw, ex)

        plot_field, data, time, x, z, _EXTRA=kw

    endfor

    !p.multi = 0

    ; --- output keywords ---
    if arg_present(b_r) then b_r = reform(bx_ext, p, p)
    if arg_present(b_phi) then b_phi = reform(by_ext, p, p)
    if arg_present(b_z) then b_z = reform(bz_ext, p, p)
    if arg_present(rgrid) then rgrid = rebin(x, p, p)
    if arg_present(zgrid) then zgrid = rebin(reform(z, 1, p), p, p)

    ; --- CSV output ---
    if n_elements(csv_outfile) gt 0 then begin
        openw, lun, csv_outfile + '_r.csv', /get_lun
        printf, lun, strjoin(strtrim(x, 2), ',')
        free_lun, lun

        openw, lun, csv_outfile + '_z.csv', /get_lun
        printf, lun, strjoin(strtrim(z, 2), ',')
        free_lun, lun

        br_2d = reform(bx_ext, p, p)
        openw, lun, csv_outfile + '_br.csv', /get_lun
        for i = 0, p-1 do $
            printf, lun, strjoin(strtrim(br_2d[i, *], 2), ',')
        free_lun, lun

        bp_2d = reform(by_ext, p, p)
        openw, lun, csv_outfile + '_bphi.csv', /get_lun
        for i = 0, p-1 do $
            printf, lun, strjoin(strtrim(bp_2d[i, *], 2), ',')
        free_lun, lun

        bz_2d = reform(bz_ext, p, p)
        openw, lun, csv_outfile + '_bz.csv', /get_lun
        for i = 0, p-1 do $
            printf, lun, strjoin(strtrim(bz_2d[i, *], 2), ',')
        free_lun, lun
    endif

end

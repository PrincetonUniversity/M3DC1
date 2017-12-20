pro ct2
    x = indgen(64)
    x[0] = 4
    x[1] = 4
    x[2] = 4
    x[4] = 4
    x[60] = 60
    x[61] = 60
    x[62] = 60
    x[63] = 60
    r = bytarr(256)
    g = bytarr(256)
    b = bytarr(256)

    r[0:63] = 0
    r[64:127] = bytscl(x)
    r[128:191] = 255
    r[192:255] = 255

    g[0:63] = bytscl(-x)
    g[64:127] = bytscl(x)
    g[128:191] = bytscl(-x)
    g[192:255] = bytscl(x)

    b[0:63] = 255
    b[64:127] = 255
    b[128:191] = bytscl(-x)
    b[192:255] = 0

    r[253:255] = [0, 0, 255]
    g[253:255] = [0, 0, 255]
    b[253:255] = [0, 0, 255]

    !p.color=253
    !p.background=255

    tvlct, r, g, b
end

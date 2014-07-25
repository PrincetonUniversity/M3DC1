pro ct3
    n = 10
    rgb = bytarr(n,3)

    rgb[0,*] = [  0,   0,   0]
    rgb[1,*] = [255,   0,   0]
    rgb[2,*] = [  0,   0, 255]
    rgb[3,*] = [  0, 128,   0]
    rgb[4,*] = [128,   0, 128]
    rgb[5,*] = [192, 192,   0]
    rgb[6,*] = [255, 128,   0]
    rgb[7,*] = [  0, 192, 192]
    rgb[8,*] = [ 64, 192,  64]
    rgb[9,*] = [192,   0, 192]
;    rgb[10,*]= [128, 128, 128]

    dx = !d.table_size/n
    rgb_big = bytarr(!d.table_size, 3)
    rgb_big[*] = 255.
    for i=0, n-1 do begin
        for j=i*dx, (i+1)*dx-1 do rgb_big[j,*] = rgb[i,*]
    end

    tvlct, rgb_big
end

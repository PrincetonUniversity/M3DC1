pro convert_bialek, nowrite=nowrite

  nphi = 16
  phi = 360.*findgen(nphi)/nphi

  file = string(format='("B_npNSTX.ph",F07.3)',phi)

  if(not keyword_set(nowrite)) then begin
     openw, ifile, 'error_field', /get_lun
     printf, ifile, format='(6A20)', $
             '%phi_tor(deg)', 'R(m)', 'Z(m)', 'B_phi', 'B_R', 'B_Z'
  end

  br0 = fltarr(nphi)
  bphi0 = fltarr(nphi)
  bz0 = fltarr(nphi)

  for i=0, nphi-1 do begin
     print, 'Converting phi = ', phi[i]
     x = read_ascii(file[i], data_start=2)
     n = n_elements(x.field01[0,*])
     s = fix(sqrt(n))
     print, format='("Assuming ",I,"x",I)', s, s
     co = cos(phi[i]*!pi/180.)
     sn = sin(phi[i]*!pi/180.)

     ; calculate fields at center of domain to do test fourier spectrum
     m = fix(n / 2)
     br0[i]   = x.field01[7,m]*co + x.field01[8,m]*sn
     bphi0[i] =-x.field01[7,m]*sn + x.field01[8,m]*co
     bz0[i]   =  x.field01[9,m]     

     if(keyword_set(nowrite)) then continue
          
     ; This is necessary because Bialek's format has Z on innermost loop
     ; whereas probe_g format has R on innermost loop
     for j=0, s-1 do begin
        for k=0, s-1 do begin
           m = s*k + j
           r = sqrt(x.field01[1,m]^2 + x.field01[2,m]^2)
           br   =  x.field01[7,m]*co + x.field01[8,m]*sn
           bphi = -x.field01[7,m]*sn + x.field01[8,m]*co
           printf, ifile, format='(6E20.11)', $
                   phi[i], r, x.field01[3,m], $
                   bphi, br, x.field01[9,m]
        end
     end
  end

  if(not keyword_set(nowrite)) then free_lun, ifile

  ct3
  window, 1
  plot, [0, 360], [-2e-3, 2e-3], /nodata, $
        xtitle = 'Toroidal Angle', ytitle='Tesla'
  oplot, phi, br0, color=color(1)
  oplot, phi, bphi0, color=color(2)
  oplot, phi, bz0, color=color(3)
  plot_legend, ['BR', 'BPHI', 'BZ'], color=[color(1),color(2),color(3)]

  br_fft = power_spectrum(br0, f=f, t=2.*!pi)
  bphi_fft = power_spectrum(bphi0, f=f, t=2.*!pi)
  bz_fft = power_spectrum(bz0, f=f, t=2.*!pi)
  f = f*2.*!pi
  window, 0
  plot, f, abs(br_fft), /ylog, /nodata, xrange=[0, 12], $
        xtitle='Toroidal Mode Number', ytitle='Power'
  oplot, f, abs(br_fft), color=color(1)
  oplot, f, abs(bphi_fft), color=color(2)
  oplot, f, abs(bz_fft), color=color(3)
  plot_legend, ['BR', 'BPHI', 'BZ'], color=[color(1),color(2),color(3)], $
               left = 0.5
end

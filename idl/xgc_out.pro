pro xgc_out, _EXTRA=extra, nplanes=nplanes, fac=fac, slice=slice

  if(n_elements(nplanes) eq 0) then nplanes=32
  if(n_elements(fac) eq 0) then fac=1.

;  psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
;  i0 = read_field('I',x,z,t,slice=-1,_EXTRA=extra)

;  psi1 = read_field('psi',x,z,t,/linear,/complex,_EXTRA=extra)
;  I1 = read_field('I',x,z,t,/linear,/complex,_EXTRA=extra)
;  f1 = read_field('f',x,z,t,/linear,/complex,_EXTRA=extra)

  br0 = read_field('bx',x,z,t,slice=-1,_EXTRA=extra)
  bphi0 = read_field('by',x,z,t,slice=-1,_EXTRA=extra)
  bz0 = read_field('bz',x,z,t,slice=-1,_EXTRA=extra)

  br1 = read_field('bx',x,z,t,/linear,/complex,slice=slice,_EXTRA=extra)
  bphi1 = read_field('by',x,z,t,/linear,/complex,slice=slice,_EXTRA=extra)
  bz1 = read_field('bz',x,z,t,/linear,/complex,slice=slice,_EXTRA=extra)
  
  ntor = read_parameter('ntor', _EXTRA=extra)

;  openw, ifile, 'm3dc1_a.out', /get_lun
  openw, ifile, 'm3dc1_b.out', /get_lun

  printf, ifile, format='("NR = ",I5,"  NZ = ",I5,"  ntor = ",I5)', $
                         n_elements(x), n_elements(z), ntor
;  printf, ifile, format='(10A16)', "R", "Z", "PSI0", "I0", $
;          "Re(PSI1)", "Im(PSI1)", "Re(I1)", "Im(I1)", "Re(f1)", "Im(f1)"
  printf, ifile, format='(11A16)', "R", "Z", "BPHI0", "BR0", "BZ0", $
          "Re(BPHI1)", "Im(BPHI1)", "Re(BR1)", "Im(BR1)", "Re(BZ1)", "Im(BZ1)"

;  for i=0, nplanes-1 do begin
;     phi = 2.*!pi*i/nplanes
;     ex = exp(complex(0,ntor)*phi)
     for j=0, n_elements(z)-1 do begin
        for k=0, n_elements(x)-1 do begin
;           printf, ifile, format='(10E16.8)', $
;                   x[k], z[j], psi0[0,k,j], i0[0,k,j], $
;                   real_part(psi1[0,k,j]), imaginary(psi1[0,k,j]), $
;                   real_part(I1[0,k,j]), imaginary(I1[0,k,j]), $
;                   real_part(f1[0,k,j]), imaginary(f1[0,k,j])
           printf, ifile, format='(11E16.8)', $
                   x[k], z[j], bphi0[0,k,j], br0[0,k,j], bz0[0,k,j], $
                   real_part(bphi1[0,k,j]), imaginary(bphi1[0,k,j]), $
                   real_part(br1[0,k,j]), imaginary(br1[0,k,j]), $
                   real_part(bz1[0,k,j]), imaginary(bz1[0,k,j])
        end
     end
;  end

   free_lun, ifile

;  close, ifile0
;  close, ifile1
   help, x, z, bphi0, br0, bz0

end

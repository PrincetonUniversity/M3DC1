;===============================================================
; convert_units
; ~~~~~~~~~~~~~
;
; converts x having dimensions d to cgs units
; where b0, n0, and l0 are the normalizations (in cgs units)
;===============================================================
pro convert_units, x, d, b0, n0, l0, zeff, mi, cgs=cgs, mks=mks, filename=filename
   if(n_elements(x) eq 0) then return

   if(not (keyword_set(cgs) or keyword_set(mks))) then return

   if(n_params() eq 2) then begin
      print, "Reading dimension from file"
      get_normalizations, filename=filename,b0=b0,n0=n0,l0=l0,zeff=zeff,ion=mi
   end

   if(b0 eq 0 or n0 eq 0 or l0 eq 0 or zeff eq 0 or mi eq 0) then begin
      print, "Warning: unknown conversion factors."
      print, "Using l0=100, B0=1e4, n0=1e14."
      l0 = 100.
      b0 = 1.e4
      n0 = 1.e14
      zeff = 1.
      mi = 1.
   endif

   val = 1.
   if(keyword_set(cgs)) then begin
       fp = (4.*!pi)
       c0 = 3.e10
       v0 = 2.18e11*b0/sqrt(mi*n0)
       t0 = l0/v0
       temp0 = b0^2/(fp*n0) * 1./(1.6022e-12)
       i0 = c0*b0*l0/fp
       e0 = b0^2*l0^3/fp
       pot0 = l0*v0*b0/c0

       val = fp^d[0] $
         * c0^d[1] $
         * n0^d[2] $
         * v0^d[3] $
         * b0^d[4] $
         * t0^d[8] $
         * l0^d[9] $
         * temp0^d[5] $
         * i0^d[6] $
         * e0^d[7] $
         * pot0^d[10]
       
   endif else if(keyword_set(mks)) then begin
       convert_units, x, d, b0, n0, l0, zeff, mi, /cgs

       val = (1.e-2)^d[1] $
         * (1.e6)^d[2] $
         * (1.e-2)^d[3] $
         * (1.e-4)^d[4] $
         * (1.e-2)^d[9] $
         * (3.e9)^(-d[6]) $
         * (1.e-7)^d[7] $
         * (3.e2)^d[10]
   end

   x = x*val
end

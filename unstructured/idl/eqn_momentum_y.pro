pro eqn_momentum_y, _EXTRA=extra, nterms=nterms, terms=term, names=names, $
               title=title, x=x, z=z
  title = 'Momentum Equation: !9P!X'
  nterms = 3
  names = ['dV/dt','-JxB','Grad(p)']
  
  complex = read_parameter('icomplex',_EXTRA=extra)
  itor = read_parameter('itor',_EXTRA=extra)
  ntor = read_parameter('ntor',_EXTRA=extra)
  rzero = read_parameter('rzero',_EXTRA=extra)
  numvar = read_parameter('numvar',_EXTRA=extra)
  ivform = read_parameter('ivform',_EXTRA=extra)

  if(complex eq 1) then begin
     gamma = read_gamma(_EXTRA=extra)

     psi0_r = read_field('psi',x,z,t,_EXTRA=extra,/equilibrium,op=2)
     psi0_z = read_field('psi',x,z,t,_EXTRA=extra,/equilibrium,op=3)
;     psi0_gs = read_field('psi',x,z,t,_EXTRA=extra,/equilibrium,op=7)
;     i0 = read_field('i',x,z,t,_EXTRA=extra,/equilibrium)
     i0_r = read_field('i',x,z,t,_EXTRA=extra,/equilibrium,op=2)
     i0_z = read_field('i',x,z,t,_EXTRA=extra,/equilibrium,op=3)
     den0 = read_field('den',x,z,t,_EXTRA=extra,/equilibrium)

     psi1_r  = read_field('psi',x,z,t,_EXTRA=extra,/linear,/complex,op=2)
     psi1_z  = read_field('psi',x,z,t,_EXTRA=extra,/linear,/complex,op=3)
;     psi1_gs = read_field('psi',x,z,t,_EXTRA=extra,/linear,/complex,op=7)
     if(numvar ge 2) then begin
;        i1 = read_field('I',x,z,t,_EXTRA=extra,/linear,/complex)
        w1 = read_field('omega',x,z,t,_EXTRA=extra,/linear,/complex)
        i1_r = read_field('I',x,z,t,_EXTRA=extra,/linear,/complex,op=2)
        i1_z = read_field('I',x,z,t,_EXTRA=extra,/linear,/complex,op=3)
        f1_r = read_field('f',x,z,t,_EXTRA=extra,/linear,/complex,op=2)
        f1_z = read_field('f',x,z,t,_EXTRA=extra,/linear,/complex,op=3)
     endif else begin
;        i1 = den0*0.
        i1_r = den0*0.
        i1_z = den0*0.
        f1_r = den0*0.
        f1_z = den0*0.
     end
       
     if(numvar eq 3) then begin
        p1 = read_field('p',x,z,t,_EXTRA=extra,/linear,/complex)
     endif else begin
        p1 = den0*0.
     end
     
     if(itor eq 1) then begin
        r = radius_matrix(x,z,t)
        rfac = complex(0,ntor)
;        psi0_gs = psi0_gs - psi0_r/r
;        psi1_gs = psi1_gs - psi1_r/r
     endif else begin
        r = 1.
        rfac = complex(0,ntor)/rzero
     end

     vy1 = w1*r

     bx0 = -psi0_z/r
;     by0 =  i0/r
     bz0 =  psi0_r/r

     jx0 = -i0_z/r
;     jy0 = -psi0_gs/r
     jz0 =  i0_r/r

     bx1 = -psi1_z/r - f1_r*rfac
;     by1 =  i1/r
     bz1 =  psi1_r/r - f1_z*rfac

     jx1 = -(i1_z + f1_z*rfac^2)/r + rfac*psi1_r/r^2
;     jy1 = -psi1_gs/r
     jz1 =  (i1_r + f1_r*rfac^2)/r + rfac*psi1_z/r^2

     d = size(den0,/dim)
     term = complexarr(nterms,d[1],d[2])
     
     ; dV / dt
     print, 'defining dv/dt'
     term[0,*,*] = den0*vy1 * gamma[0]

     ; -JxB
     print, 'defining JxB'
     term[1,*,*] = -(jz1*bx0 - jx1*bz0 + jz0*bx1 - jx0*bz1)

     ; Grad(p)
     print, 'defining grad(p)'
     term[2,*,*] = rfac*p1

     print, 'gamma = ', reform(gamma[0])
  end  

end

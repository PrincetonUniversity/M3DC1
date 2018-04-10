pro get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, $
                        ion_mass=ion_mass, _EXTRA=extra
   b0_norm = read_parameter('b0_norm', _EXTRA=extra)
   n0_norm = read_parameter('n0_norm', _EXTRA=extra)
   l0_norm = read_parameter('l0_norm', _EXTRA=extra)
   ion_mass = read_parameter("ion_mass", _EXTRA=extra)
   i = where(ion_mass eq 0, count)
   if(count gt 0) then ion_mass[i] = 1.
end

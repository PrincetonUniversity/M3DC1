; field_at_point
function field_at_point, field, x, z, x0, z0
   i = (n_elements(x)-1)*(x0-min(x)) / (max(x)-min(x))
   j = (n_elements(z)-1)*(z0-min(z)) / (max(z)-min(z))
   return, interpolate(field,i,j)
end

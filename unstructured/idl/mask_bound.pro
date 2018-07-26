; Given R and Z matrices and boundary points B (from get_boundary_path)
; return a mask matrix for R,Z points inside the boundary

function mask_bound, R, Z, B

  sr = size(R)
  mask = fltarr(sr[1], sr[2], sr[3])+1.
  
  sb = size(B)
  nb = sb[2]
    
  for k=0,nb-1 do begin
    
  
    bx2 = B[0,(k+1) mod nb]-R
    by2 = B[1,(k+1) mod nb]-Z
    bx1 = B[0,k]-R
    by1 = B[1,k]-Z

    mask = mask and ((bx2*by1 - bx1*by2) le 0.)
    
  endfor
  
  return, mask

end
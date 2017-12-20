pro read_sav_profiles, shot, time, post, dir=dir
   if(n_elements(shot) eq 0 or n_elements(time) eq 0) then begin
       print, 'Usage: read_sav_profiles, shot, time, dir=dir'
       return
   end

   if(n_elements(dir) eq 0) then begin
       dir = '.'
   endif else if(n_elements(dir) gt 1) then begin
       for i=0, n_elements(dir)-1 do begin
           read_sav_profiles, shot, time, dir=dir[i]
       end
       return
   endif

   infile = dir+'/dne'+string(format='(I0.6,".",I5.5)',shot,time)
   if(n_elements(post) ne 0) then infile=infile+'_'+post
   outfile = dir+'/profile_ne'
   restore, infile
   openw, ifile, outfile, /get_lun
   printf, ifile, transpose([[ne_str.psi_dens], [ne_str.dens/10.]])
   free_lun, ifile

   infile = dir+'/dte'+string(format='(I0.6,".",I5.5)',shot,time)
   if(n_elements(post) ne 0) then infile=infile+'_'+post
   outfile = dir+'/profile_te'
   restore, infile
   openw, ifile, outfile, /get_lun
   printf, ifile, transpose([[te_str.psi_te], [te_str.te]])
   free_lun, ifile

   infile = dir+'/dtrot'+string(format='(I0.6,".",I5.5)',shot,time)
   if(n_elements(post) ne 0) then infile=infile+'_'+post
   outfile = dir+'/profile_omega'
   restore, infile
   openw, ifile, outfile, /get_lun
   printf, ifile, transpose([[tor_rot_str.psi_tor_rot], $
                    [tor_rot_str.v_tor_rot/tor_rot_str.r_tor_rot/1000.]])
   free_lun, ifile
end

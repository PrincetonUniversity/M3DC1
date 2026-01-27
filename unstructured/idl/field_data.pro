;======================================================================
; field_data
; ~~~~~~~~~~
;
; provides the associated symbol and units for fields defined in C1.h5
;======================================================================
function field_data, name, units=units, itor=itor, filename=filename
   units = dimensions()

   if(strcmp(name, 'psi', /fold_case) eq 1 or $
      strcmp(name, 'psi_i', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=1+itor)
       return, "!7w!X"
   endif else if(strcmp(name, 'I', /fold_case) eq 1 or $
                 strcmp(name, 'I_i', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=itor)
       return, "!8I!X"
   endif else if(strcmp(name, 'phi', /fold_case) eq 1 or $
                 strcmp(name, 'phi_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=1-itor)
       return, "!8U!X"
   endif else if(strcmp(name, 'V', /fold_case) eq 1 or $
                 strcmp(name, 'V_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=-itor)
       return, "!8V!X"
   endif else if(strcmp(name, 'chi', /fold_case) eq 1 or $
                 strcmp(name, 'chi_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=1+2*itor)
       return, "!7v!X"
   endif else if(strcmp(name, 'eta', /fold_case) eq 1 or $
                 strcmp(name, 'eta_i', /fold_case) eq 1) then begin
       units = dimensions(/eta)
       return, "!7g!X"
   endif else if(strcmp(name, 'den', /fold_case) eq 1 or $
                 strcmp(name, 'den_i', /fold_case) eq 1) then begin
       units = dimensions(/n0)
       return, "!8n!Di!N!X"
   endif else if(strcmp(name, 'ne', /fold_case) eq 1 or $
                 strcmp(name, 'ne_i', /fold_case) eq 1) then begin
       units = dimensions(/n0)
       return, "!8n!De!N!X"
   endif else if(strcmp(name, 'n_re', /fold_case) eq 1 or $
                 strcmp(name, 'n_re_i', /fold_case) eq 1) then begin
       units = dimensions(/n0)
       return, "!8n!D!6RE!N!X"
   endif else if(strcmp(name, 'p', /fold_case) eq 1 or $
                 strcmp(name, 'p_i', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!8p!X"
   endif else if(strcmp(name, 'pe', /fold_case) eq 1 or $
                 strcmp(name, 'pe_i', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!8p!De!N!X"
   endif else if(strcmp(name, 'te', /fold_case) eq 1 or $
                 strcmp(name, 'te_i', /fold_case) eq 1) then begin
       units = dimensions(/temperature)
       return, "!8T!De!N!X"
   endif else if(strcmp(name, 'ti', /fold_case) eq 1 or $
                 strcmp(name, 'ti_i', /fold_case) eq 1) then begin
       units = dimensions(/temperature)
       return, "!8T!Di!N!X"
   endif else if(strcmp(name, 'sigma', /fold_case) eq 1 or $
                 strcmp(name, 'sigma_i', /fold_case) eq 1) then begin
       units = dimensions(/n0,t0=-1)
       return, "!7r!X"
   endif else if(strcmp(name, 'force_phi', /fold_case) eq 1 or $
                 strcmp(name, 'force_phi_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, l0=-1)
       return, "!8F!D!9p!N!X"
   endif else if(strcmp(name, 'pforce', /fold_case) eq 1 or $
                 strcmp(name, 'pforce_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, l0=-1)
       return, "!8F!D!9p!N!X"
   endif else if(strcmp(name, 'heat_source', /fold_case) eq 1 or $
                 strcmp(name, 'heat_source_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, t0=-1)
       return, "!8Q!X"
   endif else if(strcmp(name, 'kappa', /fold_case) eq 1) then begin
       units = dimensions(/n0, l0=2, t0=-1)
       return, "!7j!X"
   endif else if(strcmp(name, 'kappar', /fold_case) eq 1) then begin
       units = dimensions(/n0, l0=2, t0=-1)
       return, "!7j!D!9#!N!X"
   endif else if(strcmp(name, 'denm', /fold_case) eq 1) then begin
       units = dimensions(l0=2, t0=-1)
       return, "!8D!Dn!N!X"
   endif else if((strcmp(name, 'visc', /fold_case) eq 1) or $
     (strcmp(name, 'visc_c', /fold_case) eq 1) or $
                 strcmp(name, 'visc_i', /fold_case) eq 1 or $
                 strcmp(name, 'visc_c_i', /fold_case) eq 1) then begin
       units = dimensions(/p0, /t0)
       return, "!7l!X"
   endif else if(strcmp(name, 'jphi', /fold_case) eq 1  or $
                 strcmp(name, 'jphi_i', /fold_case) eq 1) then begin
       units = dimensions(/b0, l0=itor-1)
       return, "!7D!6!U*!N!7w!X"
   endif else if(strcmp(name, 'vor', /fold_case) eq 1 or $
                 strcmp(name, 'vor_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=itor-1)
       return, "!7D!6!U*!N!8U!X"
   endif else if(strcmp(name, 'com', /fold_case) eq 1 or $
                 strcmp(name, 'com_i', /fold_case) eq 1) then begin
       units = dimensions(/v0, l0=-1)
       return, "!9G.!17v!X"
   endif else if(strcmp(name, 'torque_em', /fold_case) eq 1  or $
                 strcmp(name, 'torque_em_i', /fold_case) eq 1) then begin
       units = dimensions(/p0)
       return, "!7s!D!8EM!N!X"
   endif else if(strcmp(name, 'e_r', /fold_case) eq 1  or $
                 strcmp(name, 'e_r_i', /fold_case) eq 1) then begin
       units = dimensions(/pot,l0=-1)
       return, "!8E!DR!N!X"
   endif else if(strcmp(name, 'e_phi', /fold_case) eq 1  or $
                 strcmp(name, 'e_phi_i', /fold_case) eq 1) then begin
       units = dimensions(/pot,l0=-1)
       return, "!8E!D!9P!N!X"
   endif else if(strcmp(name, 'e_z', /fold_case) eq 1  or $
                 strcmp(name, 'e_z_i', /fold_case) eq 1) then begin
       units = dimensions(/pot,l0=-1)
       return, "!8E!DZ!N!X"
   endif else if(strcmp(name, 'e_par', /fold_case) eq 1  or $
                 strcmp(name, 'e_par_i', /fold_case) eq 1) then begin
       units = dimensions(/pot,l0=-1)
       return, "!8E!D!3||!N!X"
   endif else if(strcmp(name, 'potential', /fold_case) eq 1  or $
                 strcmp(name, 'potential_i', /fold_case) eq 1) then begin
      units = dimensions(/pot)
      return, "!7U!X"
   endif else if(strcmp(name, 'frequency', /fold_case) eq 1) then begin
      units = dimensions(t0=-1)
      return, "!7x!X"
   endif else if(strcmp(name, 'kprad_sigma_i', /fold_case) eq 1) then begin
      units = dimensions(n0=1,t0=-1)
      return, "!7r!D!6KPRAD,i!N!X"
   endif else if(strcmp(name, 'kprad_sigma_e', /fold_case) eq 1) then begin
      units = dimensions(n0=1,t0=-1)
      return, "!7r!D!6KPRAD,e!N!X"
   endif else if(strcmp(name, 'kprad_totden', /fold_case) eq 1) then begin
      units = dimensions(n0=1)
      return, "!8n!D!6Z,KPRAD!N!X"
   endif else if(strcmp(name, 'kprad_rad', /fold_case) eq 1) then begin
      units = dimensions(/p0,t0=-1)
      version = read_parameter('version',filename=filename)
      if(version lt 22) then begin
         return, "!8P!D!6KPRAD,total!N!X"
      endif else begin
         return, "!8P!D!6KPRAD,line!N!X"
      end
   endif else if(strcmp(name, 'kprad_brem', /fold_case) eq 1) then begin
      units = dimensions(/p0,t0=-1)
      return, "!8P!D!6KPRAD,brem!N!X"
   endif else if(strcmp(name, 'kprad_ion', /fold_case) eq 1) then begin
      units = dimensions(/p0,t0=-1)
      return, "!8P!D!6KPRAD,ion!N!X"
   endif else if(strcmp(name, 'kprad_reck', /fold_case) eq 1) then begin
      units = dimensions(/p0,t0=-1)
      return, "!8P!D!6KPRAD,reck!N!X"
    endif else if(strcmp(name, 'kprad_recp', /fold_case) eq 1) then begin
      units = dimensions(/p0,t0=-1)
      return, "!8P!D!6KPRAD,recp!N!X"
endif else if(strcmp(name, 'rad_source', /fold_case) eq 1) then begin
      units = dimensions(/p0,t0=-1)
      return, "!8Q!D!6rad!N!X"
   endif else if(strcmp(name, 'zeff', /fold_case) eq 1) then begin
      units = dimensions()
      return, "!8Z!D!8eff!N!X"
   endif else if(strcmp(name, 'kprad_n', 7, /fold_case) eq 1) then begin
      z = fix(read_parameter('kprad_z', filename=filename))
      nz = fix(strmid(name, 8, 2))
      zstr = ['0', 'H', 'He', $
              'Li', 'Be', 'B',  'C',  'N', 'O', 'F',  'Ne', $
              'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
      nzstr = ['0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', $
                'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII', $
                'XVIII', 'XIX', 'XX']
      units = dimensions(/n0)
      return, "!8" + zstr[z] + " " + nzstr[nz] + "!X"
   endif else if(strcmp(name, 'kprad_particle_source', 21, /fold_case) eq 1) then begin
      z = fix(read_parameter('kprad_z', filename=filename))
      nz = fix(strmid(name, 22, 2))
      zstr = ['0', 'H', 'He', $
              'Li', 'Be', 'B',  'C',  'N', 'O', 'F',  'Ne', $
              'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
      nzstr = ['0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', $
                'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII', $
                'XVIII', 'XIX', 'XX']
      units = dimensions(/n0, t0=-1)
      return, "!8" + zstr[z] + " " + nzstr[nz] + " Source!X"
   endif  

   return, '!8' + name + '!X'
end

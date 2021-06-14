pro kprad_out,folder,outfile=outfile

  varnum = 21

  P     = read_scalar('p',filename=folder+'/C1.h5',/mks,time=time)                  ;1
  P_e   = read_scalar('pe',filename=folder+'/C1.h5',/mks)                 ;2
  N_e   = read_scalar('ne',filename=folder+'/C1.h5',/mks)                 ;3
  Ip    = read_scalar('it',filename=folder+'/C1.h5',/mks)                 ;4

  Prad  = read_scalar('prad',filename=folder+'/C1.h5',/mks)               ;5
  Pline = read_scalar('pline',filename=folder+'/C1.h5',/mks)              ;6
  Pbrem = read_scalar('pbrem',filename=folder+'/C1.h5',/mks)              ;7
  Pion  = read_scalar('pion',filename=folder+'/C1.h5',/mks)               ;8
  Preck = read_scalar('preck',filename=folder+'/C1.h5',/mks)              ;9
  Precp = read_scalar('precp',filename=folder+'/C1.h5',/mks)              ;10
  Pohm  = read_scalar('pohm',filename=folder+'/C1.h5',/mks)               ;11

  Erad  = read_scalar('prad',filename=folder+'/C1.h5',/mks,/integrate)    ;12
  Eline = read_scalar('pline',filename=folder+'/C1.h5',/mks,/integrate)   ;13
  Ebrem = read_scalar('pbrem',filename=folder+'/C1.h5',/mks,/integrate)   ;14
  Eion  = read_scalar('pion',filename=folder+'/C1.h5',/mks,/integrate)    ;15
  Ereck = read_scalar('preck',filename=folder+'/C1.h5',/mks,/integrate)   ;16
  Erecp = read_scalar('precp',filename=folder+'/C1.h5',/mks,/integrate)   ;17
  Eohm  = read_scalar('pohm',filename=folder+'/C1.h5',/mks,/integrate)    ;18

  temax = read_scalar('temax',filename=folder+'/C1.h5',/mks)    ;19

  Emag = read_scalar('me',filename=folder+'/C1.h5',/mks)    ;20

  NT = n_elements(time)
  out = fltarr(varnum,NT)
  out[0,*]  = time[0:Nt-1]
  out[1,*]  = P[0:Nt-1]
  out[2,*]  = P_e[0:Nt-1]
  out[3,*]  = N_e[0:Nt-1]
  out[4,*]  = Ip[0:Nt-1]
  out[5,*]  = Prad[0:Nt-1]
  out[6,*]  = Pline[0:Nt-1]
  out[7,*]  = Pbrem[0:Nt-1]
  out[8,*]  = Pion[0:Nt-1]
  out[9,*]  = Preck[0:Nt-1]
  out[10,*] = Precp[0:Nt-1]
  out[11,*] = Pohm[0:Nt-1]
  out[12,*] = Erad[0:Nt-1]
  out[13,*] = Eline[0:Nt-1]
  out[14,*] = Ebrem[0:Nt-1]
  out[15,*] = Eion[0:Nt-1]
  out[16,*] = Ereck[0:Nt-1]
  out[17,*] = Erecp[0:Nt-1]
  out[18,*] = Eohm[0:Nt-1]
  out[19,*] = temax[0:Nt-1]
  out[20,*] = Emag[0:Nt-1]

  if n_elements(outfile) eq 0 then outfile = folder+'/kprad_out.txt'

  openw,1,outfile,width=varnum*20
  printf,1,out,format='(21E20.12)'
  close,1

end

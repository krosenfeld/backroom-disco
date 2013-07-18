function mol_dat

openr,lun,'co.dat',/get_lun
sdum    = ''
specref = ''

readf,lun,sdum
readf,lun,specref
readf,lun,sdum
readf,lun,amass
readf,lun,sdum
readf,lun,nlev
readf,lun,sdum

; - read in energy levels
eterm  = dblarr(nlev)
gstat  = eterm

for i = 0, nlev-1 do begin
	readf,lun,idum,ieterm,igstat,idum,format='(I5,2x,f14.10,f6.1,I6)'
	eterm[i] = ieterm
	gstat[i] = igstat
endfor

; - read in radiative transitions
readf,lun,sdum
readf,lun,nrad
readf,lun,sdum

A21 = dblarr(nrad)
freq = A21
Eum = A21
for i = 0, nrad -1 do begin
	readf,lun,idum,idum,idum,iA21,ifreq,iEum,$
		format='(I5,x,I5,x,I5,e12.3,f16.7,f10.2)'
	A21[i] = iA21
	freq[i] = ifreq
	Eum[i] = iEum
endfor
free_lun,lun

return,{eterm:eterm,gstat:gstat,specref:specref,amass:amass,nlev:nlev,A21:A21}

stop
end

function GASmodel, disk, params, obs 
; - computes Inu(x,y) for a simple gas model.

; - relevant constants in cgs units
	@constants.pro

; - observation settings
	if (n_elements(obs) ne 8) then stop,'ERROR: 7 observation parameters'
	l	    = obs[0]	
	nu0	    = obs[1]*GHz
	nr	    = obs[2]
	nphi    = obs[3]
	nz	    = obs[4]
	Zmax	= obs[5]*AU
	veloc	= obs[6]*kms
	doquick = obs[7]

; - disk parameters
	if (n_elements(params) ne 12) then stop,'ERROR: wrong number of disk params (GASmodel)'
	thet   = params[7]*!DTOR	; - convert inclination into radians
	Mstar  = params[8]*Msun		; - convert mass of star to g
    handed = params[11]         ; - handedness of the disk

	nphi = disk.ngrid[0]
	nr   = disk.ngrid[1]
	nz   = disk.ngrid[2]
	S    = disk.S

; - conversions and derivative data
	nu	= veloc*nu0/c + nu0	; - observed nu 
	moldat  = mol_dat()		; - read in data for CO
	El      = moldat.eterm[l]*h*c	; - energy of lower level
	s0	    = h*nu/(4.*!dPI)*(2.*(l+1.)+1.)/(2.*l+1)*c*c/(2.*h*nu^3)*moldat.A21[l]

; - Define some arrays, sizes, and useful constants
	Inu	= dblarr(nphi,nr)
	kap 	= 0			      ; - no dust
	BBF1	= 2.*h/c^2.	      ; - prefactor for BB function 
	BBF2	= h/kB			  ; - exponent prefactor for BB function
	dVF1	= handed*sin(thet)*sqrt(G*Mstar) ; - velocity shift prefactor
	SignuF1 = s0*c/(nu0*sqrt(!dpi))	         ; - absorbing cross section prefactor 

; - Calculate source function and absorbing coefficient
	;dV	= veloc + dVF1*rebin(disk.X,nphi,nr,nz)/disk.r^(1.5)	        ; - velocity shift
	dV	= veloc + handed*sin(thet)*disk.Omg*rebin(disk.X,nphi,nr,nz)	; - velocity shift w/ pressure grad
	Signu	= SignuF1*exp(-dV^2/disk.dBV^2)/disk.dBV*$	                ; - absorbing cross section
		      (1.-exp(-(BBF2*nu)/disk.T))

	Knu	= disk.nl*Signu + kap*(.01+1.)*disk.rhoG		                ; - absorbing coefficient
	Snu	= BBF1*nu^3/(exp((BBF2*nu)/disk.T)-1.)			                ; - source function

	if (disk.i_notdisk[0] ne -1) then Snu[disk.i_notdisk] = 0
	if (disk.i_notdisk[0] ne -1) then Knu[disk.i_notdisk] = 0

if (doquick) then begin
	ftol = 1d-15
	quick = where(Knu ge max(Knu)*ftol,nquick)

        tphi_index = quick mod nphi
        tr_index   = (quick /long(nphi)) mod nr

        pix_index = dindgen(nphi,nr)
        mypix     = dblarr(nphi,nr) - 1
        mypix[tphi_index,tr_index]  = pix_index[tphi_index,tr_index]
        uniq_pix  = uniq(mypix,sort(mypix))

	phi_index = uniq_pix mod nphi
	r_index   = uniq_pix / long(nphi)
        ntodo = n_elements(uniq_pix)
	print, '# of pixels considered',ntodo,'/',long(nr)*long(nphi)

endif else begin
	phi_index = disk.i_xydisk mod nphi
	r_index   = disk.i_xydisk / nphi
	ntodo     = n_elements(disk.i_xydisk) 
endelse

; Trapezoid method (fast)
if (doquick) then begin
	; only calculate some of the pixels
	knu = (reform(knu,long(nr)*nphi,nz))[uniq_pix,*]
	snu = (reform(snu,long(nr)*nphi,nz))[uniq_pix,*]
	ds  = transpose(rebin((s - shift(s,1))/2.,nz,ntodo))
	arg = ds*(knu + shift(knu,0,1))
	arg[*,0] = 0.
	tau = total(arg,2,/cumulative,/double)
	arg = knu*snu*exp(-tau)
	inu[phi_index,r_index] = total((ds*(arg + shift(arg,0,1)))[*,1:nz-1],2,/double)
endif else begin
	; calculate all the pixels
	ds  = transpose(rebin((s - shift(s,1))/2.,nz,nr,nphi))
	arg = ds*(knu + shift(knu,0,0,1))
	arg[*,*,0] = 0.
	tau = total(arg,3,/cumulative,/double)
	arg = knu*snu*exp(-tau)
	inu = total((ds*(arg + shift(arg,0,0,1)))[*,*,1:nz-1],3,/double)
endelse

	return, Inu
end

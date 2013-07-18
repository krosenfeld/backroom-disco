; Initializes the disk:
;	- sets up the grid 

function initDisk, params, obs
set_plot,'X'
; - relevant constants in cgs units
	@constants.pro

; - stuff for later 
	m0      = mu*mh		; - gas mean molecular density
	Hnuctog = 0.706*mu	; - H nuclei abundance fraction (H nuclei:gas)
    sc      = 1.59d21   ; - Av --> H column density     (C. Qi 08,11)
	H2tog   = 0.8		; - H2 abundance fraction       (H2:gas)

    Tco     = 17.       ; - freeze out
    sigphot = 0.79*sc	; - photo-dissociation column

; - observation settings
	if (n_elements(obs) ne 8) then stop,'ERROR: 7 observation parameters'
	l	    = obs[0]	
	nu0	    = obs[1]*GHz
	nr	    = obs[2]
	nphi    = obs[3]
	nz	    = obs[4]
	zmax	= obs[5]*AU
	veloc	= obs[6]*kms
	doquick = obs[7]

; - disk parameters
	if (n_elements(params) ne 12) then stop, 'ERROR: wrong number of disk parameters'
    T100  = params[0]           ; - temperature at 100 AU [K]
    qq    = params[1]           ; - temperature index
    McoG  = params[2]*Msun      ; - gas mass 
    pp    = params[3]           ; - surface density index
    Rin   = params[4]*AU        ; - inner edge in cm
    Rout  = params[5]*AU        ; - outer edge in cm
    Rc    = params[6]*AU        ; - critical radius in cm 
    thet  = params[7]*!DTOR     ; - convert inclination into radians
    Mstar = params[8]*Msun      ; - convert mass of star to g
    Xco   = params[9]           ; - CO gas fraction 
    vturb = params[10]*kms      ; - CO gas fraction 

; - conversions and derivative data
	moldat  = mol_dat()		; - read in data for CO
	gl	    = 2.*l+1		; - degeneracy factor
	El      = moldat.eterm[l]*h*c	; - energy of lower level
	Te      = 2.*El/(l*(l+1.)*kB)
	costhet = cos(thet)		; - cos(i)
	sinthet = sin(thet)		; - sin(i)

; - define the desired regular cylindrical (r,z) grid
	nrc  = 256			; # of unique r points
	rmin = rin			; minimum r [AU]
	rmax = rout			; maximum r [AU]
	nzc  = 10*nrc		; # of unique z points
	zmin = 1*AU			; minimum z [AU]
	rf   = rmin*(rmax/rmin)^(findgen(nrc)/(nrc-1d0))    ; r vector (1D)
	zf   = zmin*(zmax/zmin)^(findgen(nzc)/(nzc-1d0)) 	; z vector (1D, top)
	rcf  = rf#(1+dblarr(nzc))
	zcf  = (1.+dblarr(nrc))#zf

; - interpolate dust temperature and density onto cylindrical (r,z) grid
	tf    = 0.5*!DPI-atan(zcf/rcf)	; THETA values
	rrf   = sqrt(rcf^2.+zcf^2.)	    ; R values

; - bundle the grid for helper functions
    grid = {nrc:nrc,nzc:nzc,rf:rf,rmax:rmax,zf:zf}

; - defining temperature structure  
;   using Dartois (03) type II temperature structure
   del  = 1.                                    ; shape parameter
   zq   = 20*AU*(rcf/(100.*AU))^1.3             ; height of surface layer
   tmid = (19./(155.*AU)^(-0.5))*rcf^(-0.5)     ; midplane temperature
   tatm = 60.*(rcf/(200*AU))^(-0.5)             ; atmosphere temperature
   tempg = tatm + (tmid - tatm)*cos((!dpi/(2.*zq))*zcf)^(2.*del)
   ii = where(zcf gt zq,nii)
   if (nii gt 0) then $
      tempg[ii] = tatm[ii]

; - calculate vertical density distribution
    Sc     = McoG*(2.-pp)/(2*!dpi*Rc*Rc)
    siggas = Sc*(rf/rc)^(-1*pp)*exp(-1*(rf/rc)^(2-pp))
 	rho0 = calc_hydrostatic(tempg,siggas,rcm,Mstar,grid,m0)

; - calculate radial pressure differential
    Pgas = kB/m0*rho0*tempg
    dPdr = (shift(Pgas,-1,0) - Pgas)/(shift(rcf,-1,0)-rcf)

; - calculate velocity field 
    Omg  = sqrt((dPdr/(rcf*rho0) + G*Mstar/(rcf^2 + zcf^2)^1.5))
    ;Omg  = sqrt(G*Mstar/(rcf^2 + zcf^2)^1.5)
    ;Omg  = sqrt((dPdr/(rcf*rho0) + G*Mstar/(rcf^2)^1.5))
    Omk  = sqrt(G*Mstar/rcf^3)

; - check for NANs
    ii = where(finite(Omg) eq 0, nii)
    if (nii gt 0) then $
      Omg[ii] = Omk[ii]
	ii = where(finite(rho0,/nan) eq 1,nii)
	if (nii gt 0) then begin
	  rho0[ii] = 1d-60
	  print,'Beware: removed NaNs from density (#'+string(nii,format='(I0)')+')'
	endif
	ii = where(finite(tempg,/nan) eq 1,nii)
	if (nii gt 0) then begin
	  tempg[ii] = 2.73
	  print,'Beware: removed NaNs from temperature (#'+string(nii,format='(I0)')+')'
	endif

; - find photodissociation boundary layer from top
; - (this is updated to calculate the column density of H nuclei, not H2)
	nsl  = dblarr(nzc)
	zpht = dblarr(nrc)
	for ir=0,nrc-1 do begin
	   psl = reverse(reform(Hnuctog/m0*rho0[ir,*],nzc))
	   zsl = zmax - reverse(zf)
	   foo = (zsl - shift(zsl,1))*(psl + shift(psl,1))/2.
	   foo[0] = 0
	   nsl = total(foo,/cumulative)
	   pht = where(abs(nsl) ge sigphot[0],npht)
	   if npht eq 0 then zpht[ir] = min(zmax - zsl) $
	   else zpht[ir] = max(zmax - zsl[pht])
	endfor
	zpht  = smooth(zpht,8) 	; smooth it
    szpht = zpht            ; save it

; - plot velocity structure
	;print_veloc_structure_3,rcf*Omk,rcf*Omg,tempg,rho0,grid,$
    ;zpht=szpht,Tco=Tco,zq=zq
    ;save,grid,Omk,Omg,szpht,tco,zq,rho0,tempg,zpht,filename='veloc.sav'


; - Define and initialize the cylindrical grid
	Smin	= 10*AU		        ; - offset from zero for log scale
	Smax	= 2.*zmax/costhet	; - los distance through disk
	Smid	= Smax/2.		    ; - halfway along the los
	ytop	= Smax*sinthet/2.	; - y origin offset for observer xy center

	dphi = 2.*!dPI/nphi				            ; - angular size of phi grid
;	R    =  dindgen(nr)*(Rout-Rin)/(nr-1) + Rin	; - r grid
    R    = 10^(alog10(Rin) + alog10(Rout/Rin)*dindgen(nr)/(nr-1))
	PHI  = dindgen(nphi)*2.*!dPI/(nphi-1)		    ; - phi grid
	foo  = floor(nz/2)
	S    = [Smid+Smin-10.^(alog10(Smid)+alog10(Smin/Smid)*dindgen(foo)/(foo)),$
	    	Smid-Smin+10.^(alog10(Smin)+alog10(Smid/Smin)*dindgen(foo)/(foo))]

; - arrays in [phi, r, s]
	X	= R##cos(PHI)	; - grid X coordinates 
	Y	= R##sin(PHI)	; - grid Y coordinates

; - transform grid
	tdiskZ	= zmax - costhet*S
	tdiskZ  = transpose(rebin(tdiskZ,nz,nr,nphi))
	tdiskY  = + ytop $
		      - transpose(rebin(sinthet*S,nz,nr,nphi)) $
		      + rebin(Y/costhet,nphi,nr,nz)
	tr	= sqrt(rebin(X,nphi,nr,nz)^2 + tdiskY^2)	    ; - cylindrical r
	notdisk = where(tr gt Rout)				            ; - individual grid elements not in disk
	xydisk  = where(tr[*,*,0] le  Rout + Smax*sinthet)	; - tracing outline of disk on observer xy plane

; - interpolate to calculate disk temperatures and densities
	print,'interpolating onto final grid'
	xind  = interpol(indgen(nrc),rf,tr)
	yind  = interpol(indgen(nzc),zf,abs(tdiskZ))
    tT    = interpolate(tempg,xind,yind)			  ; - gas temperature
    Omg   = interpolate(Omg,xind,yind)			  ; - gas temperature
    trhoG = H2tog*Xco/m0*interpolate(rho0,xind,yind)  ; - gas distribution
	zpht  = interpol(zpht,rf,tr)				      ; - photo-dissocation

; - photo-dissociation
    zap     = where(abs(tdiskZ) gt zpht, nzap)
    if (nzap gt 0) then trhoG[zap] = 1d-8*trhoG[zap]

; - freeze out
    zap     = where(tT le Tco, nzap)
	if (nzap gt 0) then trhoG[zap] = 1d-8*trhoG[zap]

; - temperature and turbulence broadening
	tdBV	= sqrt(2.*kB/(Da*mCO)*tT+vturb^2)	

; - approximation for partition function
    parZ    = sqrt(1.+(2./Te)^2*tT^2)			

; - calculate level population
	tnl	    = gl*trhoG*exp(-(El/kB)/tT)/parZ

; - store disk
	Disk = {X:X, Y:Y, Z:tdiskZ, S:S, ngrid:[nphi,nr,nz],$
		r:tr, T:tT, dBV:tdBV,$
		nl:tnl, rhoG:trhoG,Omg:Omg,$
		i_notdisk:notdisk, i_xydisk:xydisk}

; - print structure for inspection
	print_structure_4,rhod,rho0,tempg,grid,$
    zpht=szpht,Tco=Tco,zq=zq,zmod=zmod,incl=thet

stop
    print,'calculating sky image'
	return, Disk
end


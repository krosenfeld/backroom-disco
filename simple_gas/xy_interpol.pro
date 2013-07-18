function xy_interpol, cube, DECo, RAo, xnpix=xnpix,imres=imres,flipme=flipme

if (not keyword_set(flipme)) then flipme = 0

; initialize grids 
nchans = (size(cube))[3]
npix    = [xnpix,xnpix]
nx      = npix[0]
ny      = npix[1]
pPhi    = dblarr(nx,ny)

sqcube  = dblarr(xnpix,xnpix,nchans) 
thing = size(DECo)
nphi = thing[1]
nr   = thing[2]
sY   = transpose(rebin((findgen(npix[1])+0.5)*imres-imres*npix[1]/2.,ny,nx))
sX   =           rebin((findgen(npix[0])+0.5)*imres-imres*npix[0]/2.,nx,ny)

dx  = sX[1,0] - sX[0,0]
dy  = sY[0,1] - sY[0,0]

pR   	= sqrt(sX^2+sY^2)
pPhi	= acos(sX/pR)
pPhi[where(sY le 0)] = 2*!dPI - pPhi[where(sY le 0)]
pPhi	= reform(pPhi,nx*ny)
pR	    = reform(pR,nx*ny)

r	= rebin(transpose(sqrt(DECo^2 + RAo^2)),nr)
phi	= findgen(nphi)*2*!dPI/(nphi-1)

; - do interpolation
iR	= interpol(indgen(nr),r,pR)
iPhi	= interpol(indgen(nphi),phi,pPhi)

if (flipme) then $
    dchans = nchans/2. + 0.5 $
else dchans = nchans

for i = 0, dchans-1 do begin
	im     	    	= cube[*,*,i]
	sqcube[*,*,i]	= reform(interpolate(im,iPhi,iR),nx,ny)
endfor

if (flipme) then begin
sX   = -1*rebin((findgen(npix[0])+0.5)*imres-imres*npix[0]/2.,nx,ny)
dx  = sX[1,0] - sX[0,0]

pR   	= sqrt(sX^2+sY^2)
pPhi	= acos(sX/pR)
pPhi[where(sY le 0)] = 2*!dPI - pPhi[where(sY le 0)]
pPhi	= reform(pPhi,nx*ny)
pR	    = reform(pR,nx*ny)

r	= rebin(transpose(sqrt(DECo^2 + RAo^2)),nr)
phi	= findgen(nphi)*2*!dPI/(nphi-1)

; - do interpolation
iR	= interpol(indgen(nr),r,pR)
iPhi	= interpol(indgen(nphi),phi,pPhi)

for i = dchans, nchans-1 do begin
	im     	    	= cube[*,*,i]
	sqcube[*,*,i]	= reform(interpolate(im,iPhi,iR),nx,ny)
endfor

endif

out = {xycube:sqcube, sX:sX, sY:sY}
return, out
 
end

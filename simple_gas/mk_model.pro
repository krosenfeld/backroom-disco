function mk_model,mfile,params,obs

@constants.pro
@user_inputs.pro

; - preliminaries
xpixscale = imsize/(xnpix-1.)
dd  	 = distance*pc		; - distance in cm
arcsec   = rad/dd		; - angular conversion factor (cm to arcsec)
chans    = chanmin + findgen(nchans)*chanstep
tchans   = string(chans,format='(D0.1)')
print, 'running ', modfile, ' for channels (km/s):'
print, tchans
    
disk = initdisk(params,obs)
    cube = dblarr(disk.ngrid[0],disk.ngrid[1],nchans)
    X = disk.X
    Y = disk.Y

; - do the calculation 
for i=0,(nchans-1) do begin
	obs[6] = chans[i]
	Inu = GASmodel(disk,params,obs)
	cube[*,*,i] = Inu
	print,'FINISHED CHANNEL ' + string(i+1,format='(I2)') $
		+ '/' + string(nchans,format='(I2)')
endfor

; - interpolate onto a square grid
out = xy_interpol(cube, X*arcsec, Y*arcsec, xnpix=xnpix,corner=-1*[imsize,imsize]/2.)
im  = out.xycube

; - make header
hdr = write_h(nchans=nchans,dd=distance,xnpix=xnpix,xpixscale=xpixscale,$
		lstep=chanstep)

; - shift and rotate model
im_s = fltarr(sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2'),sxpar(hdr,'NAXIS3'))
for iv=0,sxpar(hdr,'NAXIS3')-1 do begin

   ; - rotate to the desired position angle
   im_s[*,*,iv] = rot(im[*,*,iv],90.-PA,1.,sxpar(hdr,'CRPIX1')-1,$
                      sxpar(hdr,'CRPIX2')-1,/INTERP)

   ; - integer-pixel shift
   os = [-1.,1.]*offs/(3600.*abs([sxpar(hdr,'CDELT1'),sxpar(hdr,'CDELT2')]))
   if (fix(os[0])-os[0] eq 0. and fix(os[1])-os[1] eq 0.) then $
     im_s[*,*,iv] = shift(im_s[*,*,iv],os[0],os[1])

   ; - sub-pixel registration
   sz = size(im)
   xn = findgen(sz[1])#replicate(1.,sz[2]) & x1 = (xn-os[0])>0<(sz[1]-1.)
   yn = replicate(1.,sz[1])#findgen(sz[2]) & y1 = (yn-os[1])>0<(sz[2]-1.)
   im_s[*,*,iv] = interpolate(im_s[*,*,iv],x1,y1)
endfor

; - write processed model
writefits,mfile,im_s*Jy*(xpixscale/rad)^2,hdr

return,0
end

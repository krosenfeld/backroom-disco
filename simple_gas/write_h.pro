function write_h,nchans=nchans,dd=dd,xnpix=xnpix,xpixscale=xpixscale,lstep=lstep

cen	= [xnpix/2.+.5,xnpix/2.+.5]	; - central pixel location

sxaddpar,hdr,'SIMPLE','T'
sxaddpar,hdr,'NAXIS',3
sxaddpar,hdr,'NAXIS1',xnpix
sxaddpar,hdr,'NAXIS2',xnpix
sxaddpar,hdr,'NAXIS3',nchans
sxaddpar,hdr,'CDELT1',-1.*xpixscale/3600.d
sxaddpar,hdr,'CRPIX1',cen[0]
sxaddpar,hdr,'CRVAL1',0
sxaddpar,hdr,'CTYPE1','RA---SIN'
sxaddpar,hdr,'CDELT2',xpixscale/3600.d
sxaddpar,hdr,'CRPIX2',cen[1]
sxaddpar,hdr,'CRVAL2',0
sxaddpar,hdr,'CTYPE2','DEC--SIN'
sxaddpar,hdr,'CTYPE3','VELO-LSR'
sxaddpar,hdr,'CDELT3',lstep*1000	; - dv in m/s
;sxaddpar,hdr,'CRPIX3',nchans/2		; - 0 velocity channel
sxaddpar,hdr,'CRPIX3',nchans/2+1	; - 0 velocity channel
sxaddpar,hdr,'CRVAL3',0			; - has 0 velocity
sxaddpar,hdr,'OBJECT','model'

return,hdr

end

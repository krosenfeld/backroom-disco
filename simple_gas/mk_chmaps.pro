pro mk_chmaps,offs=offs,PA=PA,incl=incl
if (not keyword_set(offs)) then offs=[0,0]
if (not keyword_set(PA))   then PA = 90.
if (not keyword_set(incl)) then incl=0.
no_axis = replicate(' ',10)

; read in file
cube = readfits('model/debugp.model.fits',hdr)

; make axes
;ra = sxpar(hdr,'CRVAL1')+$
;       sxpar(hdr,'CDELT1')*$
;       (findgen(sxpar(hdr,'NAXIS1'))-sxpar(hdr,'CRPIX1')+1)
;de = sxpar(hdr,'CRVAL1')+$
;       sxpar(hdr,'CDELT1')*$
;       (findgen(sxpar(hdr,'NAXIS1'))-sxpar(hdr,'CRPIX1')+1)
ra =  3600d*sxpar(hdr,'CDELT1')*$
      (findgen(sxpar(hdr,'NAXIS1'))-(sxpar(hdr,'NAXIS1')/2.-0.5))

; crop image
imx = 10.    ; image size in arcseconds
de  = -1*ra
ra  = ra    - offs[0]
de  = de - offs[1]
ira = where(abs(ra) lt imx/2.)
ide = where(abs(de) lt imx/2.)
cube = cube[ira[0]:ira[-1],ide[0]:ide[-1],*]
ra = ra[ira]
de = de[ide]

; make figure
nx  = 5
ny  = 3
idx = 1.
ixm = 0.46*[1,1]
iym = [0.35,0.05]
xs  = nx*idx + total(ixm)
idy = idx
ys  = total(iym) + ny*idy
xm  = ixm/xs
ym  = iym/ys
dx  = idx/xs
dy  = idy/ys

levs = [0.1,0.2,0.3,0.5,0.7,0.9,1.1,1.3,1.5]
set_plot,'PS'
device,filename='chmaps.eps',xs=xs,ys=ys,/inches,/encaps
for iy=0,ny-1 do begin
  for ix=0,nx-1 do begin
    pos = [xm[0]+ix*dx,ym[0]+dy*(ny-1-iy),xm[0]+(ix+1)*dx,ym[0]+dy*(ny-iy)]
    if ix eq 0 and iy eq ny-1 then $
    contour,cube[*,*,ix+iy*nx],ra,de,xra=imx*[0.5,-0.5],yra=imx*[-0.5,0.5],$
        /xsty,/ysty,xthick=3,ythick=3,chars=0.8,charthi=2,c_thick=2,$
        /noera,pos=pos,xtit='!7Da!6 ["]',ytit='!7Dd!6 ["]',lev=levs $
    else $
    contour,cube[*,*,ix+iy*nx],ra,de,xra=imx*[0.5,-0.5],yra=imx*[-0.5,0.5],$
        /xsty,/ysty,xthick=3,ythick=3,c_thick=2,$
        /noera,pos=pos,xtickn=no_axis,ytickn=no_axis,lev=levs
     oplot,0.9*imx*[-1,1]/2.,!DTOR*(90.-[PA,PA]),linesty=0,/polar,color=90
     oplot,0.9*imx*[-1,1]/2.*cos(incl*!DTOR),!DTOR*(180.-[PA,PA]),linesty=0,/polar,color=90
  endfor
endfor
device,/close
set_plot,'X'
end

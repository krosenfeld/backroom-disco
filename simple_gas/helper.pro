; *****************************************************
function calc_hydrostatic,tempg,siggas,rcm,mstar,grid,m0
; *****************************************************

; constants
	G	= 6.67259d-8	; - Gravitational constant (cm^3/g/s^2)
	kB	= 1.3807d-16	; - Boltzmann's constant (erg/K)

; grad grid
	nrc = grid.nrc
	nzc = grid.nzc
	zf  = grid.zf
	rf  = grid.rf

; compute rho structure
	rho0   = dblarr(nrc,nzc)
	;sigint = interpol(siggas,rcm,rf)
	sigint = siggas 
	dlnp   = dblarr(nzc)
	lnp    = dblarr(nzc)
	print,'doing hydrostatic equilibrium'
	for ir=0,nrc-1 do begin
	   ; compute the gravo-thermal constant
	   grvc = G*Mstar*m0/kB

	   ; extract the T(z) profile at a given radius
	   T = reform(tempg[ir,*],nzc)

	   ; differential equation for vertical density profile 
	   dz      = (zf - shift(zf,1))
	   dlnT    = (alog(T) - shift(alog(T),1))/dz
	   dlnp    = -1.*grvc*zf/(T*(rf[ir]^2.+zf^2.)^(1.5)) - dlnT
	   dlnp[0] = -1.*grvc*zf[0]/(T[0]*(rf[ir]^2+zf[0]^2)^(1.5))

	   ; numerical integration to get vertical density profile
	   foo = (zf - shift(zf,1))*(dlnp + shift(dlnp,1))/2.
	   foo[0] = 0d
	   lnp = total(foo,/cumulative)

	   ; normalize the density profile (note: this is just half the sigma value!)
	   dens = 0.5*sigint[ir]*exp(lnp)/int_tabulated(zf,exp(lnp),/DOUBLE)
     
	   ; populate the slice in the 2-D density grid
	   rho0[ir,*] = dens
	endfor

  return,rho0

end

; *****************************************************
pro print_structure_4,rhod,rho0,tempg,grid,$
	zpht=zpht,Tco=Tco,zq=zq,zmod=zmod,incl=incl
; *****************************************************

; - plot the 2-D temperature and density distributions in cylindrical coords
; - also include gas-dust temperature differential

mu = 2.37               ; mean molecular weight of gas
mh = 1.67d-24           ; proton mass
AU = 1.496d13           ; AU

ctindex = 100		; color table index
fcols = 5+8*findgen(35)		
dlevs = -24.*(10./20.)^(findgen(35)/34.)-alog10(mu*mh)
tlevs = [10,20,30,40,50,75,100,150]
tlabs = 1 + intarr(n_elements(tlevs))

idx = 3.0
idy = 2.175
nx  = 1.
ixm = [0.4,0.1]
iym = [0.3,0.55]
iba = 0.15
xs  = total(ixm) + idx*nx
ys  = total(iym) + idy
dx  = idx/xs
dy  = idy/ys
xm  = total(ixm/xs ,/cumulative) + dx*findgen(nx)
ym  = iym/ys
ba  = iba/ys
pos1 = [xm[0],ym[0],xm[0]+dx,ym[0]+dy]
posb = [xm[0]+dx/5,ym[0]+dy+ba/2.,xm[0]+dx,ym[0]+dy+ba*3./2]

set_plot,'PS'
device,filename='model/test.structures.eps',/encapsulated,/inches,/color,bits=8,$
       xsize=xs,ysize=ys

xra  = [0,600]
yra  = [0,150]
tickl = 0.03

; grab grid info
rf = grid.rf
zf = grid.zf
; gas
if (ctindex ne 100) then loadct,ctindex else begin
  readcol,'/home/krosenfe/IDLWorkspace/helpers//mycolors3.tbl',r,g,b,$
          format='(i,i,i)',/SIL
  tvlct,r,g,b
endelse
contour,alog10(rho0/(mu*mh)),rf/AU,zf/AU,/fill,levels=dlevs,c_colors=fcols,$
        xra=xra,xsty=5,xthick=3,ythick=3,charsize=0.7,$
        charthick=2,position=pos1,/NOERASE,yra=yra,ysty=5
mkct
contour,alog10(rho0/(mu*mh)),rf/AU,zf/AU,levels=dlevs,xra=xra,/xsty,$
        yra=yra,/ysty,xthick=3,ythick=3,charsize=0.7,charthick=2,/NODATA,$
        position=pos1,/NOERASE,xtit='!8r!6 [AU]',xminor=5,yticki=50,$
        ytit='!8z!6 [AU]'
;oplot,rf/AU,rf/AU*tan(20*!dtor),thick=6,linesty=1
;if (keyword_set(zmod)) then oplot,rf/AU,zmod/AU,thick=3,linesty=1,col=4
if (keyword_set(Tco))   then contour,tempg,rf/AU,zf/AU,levels=Tco,/over,c_thick=5,$
    c_color=0,c_linesty=2 
if (keyword_set(zpht))  then oplot,rf/AU,zpht/AU,thick=5,linesty=2,color=0
orl = 50.
orx = 20.  + orl*[0., sin(incl)]
ory = 190. - orl*[0., cos(incl)]
;oplot,orx,ory,thick=3
contour,tempg,rf/AU,zf/AU,levels=tlevs,/over,c_thick=2.5,$
        c_labels=tlabs,c_col=120,$
        c_charsize=0.6,c_charthick=2

; - color scalebar
fmin = min(dlevs)
fmax = max(dlevs)
if (ctindex ne 100) then loadct,ctindex else begin
  readcol,'/home/sandrews/IDL_stash/stash/mycolors3.tbl',r,g,b,$
          format='(i,i,i)',/SIL
  tvlct,r,g,b
endelse
tv,bytscl(replicate(1B,10)##bindgen(256),top=254),posb[0]*xs,posb[1]*ys,$
   xsize=(posb[2]-posb[0])*xs,ysize=(posb[3]-posb[1])*ys,/inches
mkct
plot,indgen(11),fmin+(fmax-fmin)*dindgen(11)/10.,yra=[0,10],/ysty,$
     xra=[fmin,fmax],/xsty,position=posb,/NOERASE,/NODATA,$
     xtickname=replicate(' ',11),xthick=3.,ythick=3.,charthick=2,$
     charsize=0.7,yticklen=1d-30,ytickname=replicate(' ',11),xticklen=10*tickl,$
     xminor=1
axis,xaxis=1,xra=[fmin,fmax],/xsty,charsize=0.7,charthick=2,ythick=3.,$
    xticklen=10*tickl,xminor=1
;xyouts,posb[2]+2*ba,0.5*(posb[3]+posb[1]),'!6log !8n!6 [cm!E-3!N]',$
xyouts,0.5*(posb[2]+posb[0]),posb[3]+1.2*ba,'!6log !8n!6 [cm!E-3!N]',$
       /NORMAL,align=0.5,orientation=0,charsize=0.7,charthick=2
device,/close
end

; *****************************************************
pro print_veloc_structure_3,Vk,Vg,tempg,rho0,grid,$
	zpht=zpht,Tco=Tco,zq=zq
; *****************************************************

; - plot the 2-D temperature and density distributions in cylindrical coords
; - also include gas-dust temperature differential

mu = 2.37               ; mean molecular weight of gas
mh = 1.67d-24           ; proton mass
AU = 1.496d13           ; AU

ctindex = 33		; color table index
fcols = 5+8*findgen(35)		
dlevs = -24.*(10./20.)^(findgen(35)/34.)-alog10(mu*mh)
tlevs = [10,20,30,40,50,75,100,120,140]
tlabs = 1 + intarr(n_elements(tlevs))
;vra  = [-0.15,0.1]*1d5
vra  = [-0.25,0.05]*1d5
ncol = 50
vlevs = (vra[0] + (vra[1]- vra[0])*findgen(ncol-1)/(ncol-2))
vlevs = [-100*1d5,vlevs]
fcols = 0 + (255)*findgen(ncol)/(ncol-1)

idx = 3.0
idy = 2.175
nx  = 1.
ixm = [0.4,0.1]
iym = [0.3,0.55]
iba = 0.15
xs  = total(ixm) + idx*nx
ys  = total(iym) + idy
dx  = idx/xs
dy  = idy/ys
xm  = total(ixm/xs ,/cumulative) + dx*findgen(nx)
ym  = iym/ys
ba  = iba/ys
pos1 = [xm[0],ym[0],xm[0]+dx,ym[0]+dy]
posb = [xm[0]+dx/5,ym[0]+dy+ba/2.,xm[0]+dx,ym[0]+dy+ba*3./2]

set_plot,'PS'
device,filename='model/veloc.structures.eps',/encapsulated,/inches,/color,bits=8,$
       xsize=xs,ysize=ys

xra  = [0,600]
;yra  = [0,300]
yra  = [0,150]
tickl = 0.03

; grab grid info
rf = grid.rf
zf = grid.zf

; ########################
; relative error
; #######################

if (ctindex ne 100) then loadct,ctindex else begin
  readcol,'/home/krosenfe/IDLWorkspace/helpers//mycolors3.tbl',r,g,b,$
          format='(i,i,i)',/SIL
  tvlct,r,g,b
endelse

dvlevs = -1*reverse([-0.05,0,0.05,0.1,0.15])
dvlabs = 1 + intarr(n_elements(dvlevs))
vra  = [-.1,0.01]
ncol = 50
vlevs = (vra[0] + (vra[1]- vra[0])*findgen(ncol-1)/(ncol-2))
vlevs = [-100,vlevs]
boo = (Vg - Vk)/Vk
ii = where(rho0/(mu*mh) lt 1d2,nii)
if (nii gt 0) then boo[ii] = vlevs[0]-100
contour,boo,rf/AU,zf/AU,/fill,levels=vlevs,c_col=fcols,$
        xra=xra,xsty=5,xthick=3,ythick=3,charsize=0.7,$
        charthick=2,position=pos1,/NOERASE,yra=yra,ysty=5
mkct
contour,boo,rf/AU,zf/AU,levels=vlevs,xra=xra,/xsty,$
        yra=yra,/ysty,xthick=3,ythick=3,charsize=0.7,charthick=2,/NODATA,$
        position=pos1,/NOERASE,xtickn=no_axis,/closed,yticki=50,$
        ticklen=tickl,xminor=5,ytit='!8z!6 [AU]',xtit='!8r!6 [AU]'
boo = (Vg - Vk)*1d-5
ii = where(rho0/(mu*mh) lt 1d2,nii)
if (nii gt 0) then boo[ii] = dvlevs[0]-100
contour,boo,rf/AU,zf/AU,levels=dvlevs,/overplot,c_thick=3,c_labels=dvlabs,$
    c_charthick=2,c_charsize=0.6,c_col=120
if (keyword_set(Tco))   then contour,tempg,rf/AU,zf/AU,levels=Tco,/over,c_thick=5,$
    c_color=0,c_linesty=2 
if (keyword_set(zpht))  then oplot,rf/AU,zpht/AU,thick=5,linesty=2,color=0

; - color scalebar
if 1 then begin
fmin = 100*vra[0]
fmax = 100*vra[1]
if (ctindex ne 100) then loadct,ctindex else begin
  readcol,'/home/sandrews/IDL_stash/stash/mycolors3.tbl',r,g,b,$
          format='(i,i,i)',/SIL
  tvlct,r,g,b
endelse
tv,bytscl(replicate(1B,10)##bindgen(256),top=254),posb[0]*xs,posb[1]*ys,$
   xsize=(posb[2]-posb[0])*xs,ysize=(posb[3]-posb[1])*ys,/inches
mkct
plot,indgen(11),fmin+(fmax-fmin)*dindgen(11)/10.,yra=[0,10],ysty=9,$
     xra=[fmin,fmax],xsty=1,position=posb,/NOERASE,/NODATA,$
     ytickname=replicate(' ',11),xthick=3.,ythick=3.,charthick=2,$
     charsize=0.7,yticklen=1d-30,xtickname=replicate(' ',11),$
     xticklen=10*tickl,xminor=1
axis,xaxis=1,xra=[fmin,fmax],/xsty,charsize=0.7,charthick=2,xthick=3.,$
   xticklen=10*tickl,xminor=1,xtickf='(I3)'
xyouts,0.5*(posb[2]+posb[0]),posb[3]+1.2*ba,'!7d!8v [%]!6',$
       /NORMAL,align=0.5,orientation=0,charsize=0.7,charthick=2
endif
device,/close
end

;@constants.pro
root	 = 'debug'		; - root name for the model
doquick  = 0			; - speed up the calculation

; - disk info
Xco      = 1d-4         ; - CO/H2 fraction
T100     = -1			; - atmosphere gas temperature at 10 AU (K)
qq       = 0.9		   	; - atmosphere gas temperature index
McoG	 = 0.09			; - Total mass of the disk
Rc       = 150.			; - characteristic radius
pp       = 1. 	        ; - gas surface density index
Mstar    = 2.3			; - mass of the star (msun)
vturb    = 0.01			; - turbulent velocity (km/s)
handed   = -1           ; - 1 (right-handed), -1 (left-handed)
flipme   = 1            ; - reflect about the v axis

; - system info
distance = 122.d		    ; - distance (pc)
PA       = 132.+180 	    ; - position angle (deg)
incl     = 48.			    ; - inclination (deg)
vsys     = -1d       	    ; - LSR of system (km/s)
offs     = [0., 0.]	; - phase center offset ("; + = E, N)

; - observation info
Jlow     = 2			; - lower J quantum number
freq     = 345.79599d	; - transition frequency (GHz)

datfile  = 'data/HD163296.CO32.regridded.cen15' 
nchans   = 15			; - number of velocity channels 
chanmin  = -2.24d		; - minimum channel velocity (km/s)
chanstep = 0.32d		; - channel velocity step (km/s)

imres    = 0.04			; - resolution in arc seconds 
xnpix    = 512L			; - number of pixels along one side

; - grid info
Rin      = 10d			; - inner edge (AU) 
Rout     = 700.d		; - outer edge (AU)
nr       = 50L			; - number of r bins
nphi     = 101L			; - number of phi bins (doing full 2 PI) 
nz	     = 400L		    ; - number of bins in Z direction (even #)
TOPZ	 = 170d 		; - height above 0 to start integration (AU)

; NOTE: Other parameters of interest are contained in constants.pro.
; - bundle the model parameters
obs     = [Jlow, freq, nr, nphi, nz, TOPZ, 0., doquick]                ; - observation params  
params  = [T100,qq,McoG,pp,Rin,Rout,Rc,incl,Mstar,Xco,vturb,handed]    ; - disk params
modfile  = 'model/'+root

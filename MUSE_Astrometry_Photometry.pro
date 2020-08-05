;ToDo:
;better fit of FWHM

@circint_MOD.pro

pro MUSE_Astrometry_Photometry

;===============================================================================================

plsc = 25.	;mas/px

guesspos = [115,102]
psfrad = 20.	;radius to fit the star
planetrad = 10.	;radius to fit the planet

filepca = 'PCA7_FinalCube.fits'

readcol, '/home/amueller/work/IDLlibs/AO/MUSE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
; selp = 1
base = tmp[selp-1]
datadir = base+'Reduced/'
resdir = base+'HRSDI/'

;===============================================================================================

cube = mrdfits(resdir+'FinalCube.fits', 0, /silent)
wl = mrdfits(resdir+'FinalCube.fits', 1, /silent)
restore, resdir+'frame_idx.sav'

cube = cube[*,*,idxgood]
wl = wl[idxgood]

dim = n_elements(cube[*,0,0])
nframes = n_elements(idxgood)

;===============================================================================================

; ;measure PSF as function of lambda
; 
; starfwhm = dblarr(4,nframes)
; starpos = dblarr(4,nframes)
; for i=0,nframes-1 do begin
; 
; 	tmpfwhm = (wl[i]*1.d-10)/8.2*206265./25.d-3
; 	tim = cube[dim/2-psfrad:dim/2+psfrad-1,dim/2-psfrad:dim/2+psfrad-1,i]
; 	estimates = [median(tim), max(tim), tmpfwhm*5., tmpfwhm*5., dim/2., dim/2., 0.01, 1.]
;   xa = dindgen(2.*psfrad)
;   ya = xa
;   pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},8)
;   pi[2].limited[0] = 1
;   pi[2].limits[0] = 0.
;   pi[2].limits[1] = 20.
; 	pi[3].limited[0] = 1
;   pi[3].limits[0] = 0.
;   pi[3].limits[1] = 20.
;   
; 	yfit = mpfit2dpeak(tim, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, parinfo=pi, /quiet)
; 	starfwhm[0:1,i] =  2.*a[2:3]	;2.*sqrt(2.*alog(2.))*[a[2:3]]
; 	starfwhm[2:3,i] =  2.*perror[2:3]	;2.*sqrt(2.*alog(2.))*[perror[2:3]]
; 	starpos[0:1,i] = a[4:5]
; 	starpos[2:3,i] = perror[4:5]
; 	
; 	;window, 0, xs=400, ys=400
; 	;cgimage, tim, stretch=2
; 	;window, 1, xs=400, ys=400
; 	;cgimage, yfit, stretch=2
; 	;hak
; 	
; 	proceeding_text,loop=nframes, i=i, prompt='> Fit PSF Frame   '+string(i+1,form='(I4)')
; 	
; endfor
; 
; ;fit fwhm
; starfwhm_fit = dblarr(2,nframes)
; ;x direction
; c = robust_poly_fit(wl, reform(starfwhm[0,*]), 5., yfit, /double, numit=100)
; tmp = reform(starfwhm[0,*])/yfit
; resistant_mean, tmp, 3., t1, t2, goodvec=good
; c = robust_poly_fit(wl[good], reform(starfwhm[0,good]), 10., yfit, /double, numit=100)
; starfwhm_fit[0,*] = poly(wl, c)
; ;y direction
; c = robust_poly_fit(wl, reform(starfwhm[1,*]), 5., yfit, /double, numit=100)
; tmp = reform(starfwhm[1,*])/yfit
; resistant_mean, tmp, 3., t1, t2, goodvec=good
; c = robust_poly_fit(wl[good], reform(starfwhm[1,good]), 10., yfit, /double, numit=100)
; starfwhm_fit[1,*] = poly(wl, c)
; 
; save, wl, starfwhm, starpos, starfwhm_fit, filename=resdir+'Values.sav'
; 
;===============================================================================================

;PCA cube, extract spectrum and astrometry

restore, resdir+'Values.sav'

cube = 0
cube = mrdfits(resdir+filepca, 0, /silent)
cube = cube[*,*,idxgood]

;find the frame with highest signal

; wlrange = [6563., 6570.]
wlrange = [6555., 6575.]
idx = where(wl ge wlrange[0] and wl le wlrange[1])
;haim = total(cube[*,*,idx], 3)
sum = dblarr(n_elements(idx))
for i=0,n_elements(idx)-1 do sum[i] = total(cube[guesspos[0]-planetrad:guesspos[0]+planetrad-1,guesspos[1]-planetrad:guesspos[1]+planetrad-1,idx[i]])

tmp = max(sum, idxmax)
idxmax = idx[idxmax]

tim = median(cube[guesspos[0]-planetrad:guesspos[0]+planetrad-1,guesspos[1]-planetrad:guesspos[1]+planetrad-1,idxmax-1:idxmax+1], dim=3)
cdim = n_elements(tim[*,0])

estimates = [median(tim), max(tim), starfwhm_fit[0,idxmax], starfwhm_fit[1,idxmax], cdim/2., cdim/2., 0.01, 1.]
xa = dindgen(cdim)
ya = xa
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},8)
; pi[2].limited[0] = 1
; pi[2].limits[0] = 0.
; pi[2].limits[1] = 20.
; pi[3].limited[0] = 1
; pi[3].limits[0] = 0.
; pi[3].limits[1] = 20.
pi[2].fixed = 1
pi[3].fixed = 1

yfit = mpfit2dpeak(tim, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, parinfo=pi, /quiet)

pos_px = dblarr(6)	;x, y, dist, xe, ye, diste
pos_mas = dblarr(6)
posang_deg = dblarr(2)

;x negative because East to the left

x = -1.*((guesspos[0]-dim/2.+(a[4]-cdim/2.))-(starpos[0,idxmax]-psfrad))
y = (guesspos[1]-dim/2.+(a[5]-cdim/2.))-(starpos[1,idxmax]-psfrad)
pos_px[0] = x
pos_px[1] = y
pos_px[2] = sqrt(x^2.+y^2.)
pos_mas[0:2] = pos_px[0:2]*plsc

;compute PA
rarad1 = 0.
rarad2 = (pos_mas[0]/(1.d3*60.*60.))*!dtor
dcrad1 = 0.
dcrad2 = (pos_mas[1]/(1.d3*60.*60.))*!dtor
radif  = rarad2-rarad1
angle  = atan(sin(radif),cos(dcrad1)*tan(dcrad2)-sin(dcrad1)*cos(radif))
posang_deg[0] = angle/!dtor
if (angle/!dtor lt 0.) then posang_deg[0] = 360.+angle/!dtor
if (angle/!dtor gt 360.) then posang_deg[0] = 360.-angle/!dtor

save, wl, starfwhm, starpos, starfwhm_fit, pos_px, pos_mas, posang_deg, filename=resdir+'Values.sav'

;===============================================================================================

;extract spectrum

restore, resdir+'Values.sav'
ndata = n_elements(wl)
flux = dblarr(2,ndata)
fluxsky = dblarr(2,ndata)

for xx=0,ndata-1 do begin


	;empty aperturs distributed next to the target aperture in the same annuli as noise has a radial dependence

	tim = cube[*,*,xx]
	photapr = mean([starfwhm_fit[0,xx],starfwhm_fit[1,xx]])

	dx = dim/2.+(-1.*pos_px[0]);+1
	dy = dim/2.+pos_px[1];+1

	maxaper_full = (2.*!DPI*pos_px[2])/(photapr)
	maxaper = floor(maxaper_full)	;maximum number of apertures (incl. science aperture)
	naper = maxaper
	alpha = 360.d0/maxaper	;angle between two consecutive apertures in deg
	zero = posang_deg[0]+90.d0
	alphas = linspace(zero-((naper-1)/2.)*alpha, zero+((naper-1)/2.)*alpha, naper)

	xap = pos_px[2]*cos(alphas*!dtor)+dim/2.;+1.
	yap = pos_px[2]*sin(alphas*!dtor)+dim/2.;+1.

	;marking science aperture
	idx = closest(alphas, zero)
	flagap = uintarr(naper)
	flagap[idx] = 1	

	;window, 0, xs=700, ys=700
	;cgimage, scale_image_am(tim), /axis
	;for j=0,naper-1 do tvcircle, /data, photapr/2., xap[j], yap[j], color=cgcolor('red')
	;tvcircle, /data, photapr/2., xap[idx], yap[idx], color=cgcolor('yellow')

	tmpflux = dblarr(naper)
	for i=0,naper-1 do begin
	
		circint_MOD, tim, xap[i], yap[i], photapr/2., tot, mtot, meantot, maxpx, sdtot, npx, totdif, npxdif, t8
		tmpflux[i] = tot
		
	endfor
	
	idx = where(flagap ne 1)
	fluxsky[0,xx] = median(tmpflux[idx])
	fluxsky[1,xx] = stddev(tmpflux[idx])
	idx = where(flagap eq 1)
	flux[0,xx] = tmpflux[idx]-fluxsky[0,xx]
	flux[1,xx] = fluxsky[1,xx]
	
	proceeding_text,loop=ndata, i=xx, prompt='> Photometry Frame   '+string(xx+1,form='(I4)')

endfor



stop
end

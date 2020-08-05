;ToDo
;detection of bad pixels

@bsort.pro
@fftshift.pro
@scale_image_am.pro
@shift_sub.pro
@headfits.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@gettok.pro
@valid_num.pro
@mrd_skip.pro
@get_eso_keyword.pro
@sxpar.pro
@sky.pro
@mmm.pro
@asinh.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cggetcolorstate.pro
@cgerase.pro
@cgsetcolorstate.pro
@cgcolor.pro
@cgsnapshot.pro
@cgresizeimage.pro
@cgplot.pro
@setdecomposedstate.pro
@decomposedcolor.pro
@colorsareidentical.pro
@mpfit2dpeak.pro
@mpfit.pro
@mpfit2dfun.pro
@array_indices.pro
@fixpix_mod.pro
@dist_circle.pro
@proceeding_text.pro
@writefits.pro
@check_fits.pro
@fxaddpar.pro
@sxdelpar.pro
@mkhdr.pro
@sxaddpar.pro
@poly_smooth.pro
@get_date.pro
@daycnv.pro
@mrdfits.pro
@fxposit.pro
@dist.pro
@tsum.pro
@readcol.pro
@remchar.pro
@cgcolor24.pro
@cgbitget.pro
@convert_to_type.pro
@cgcheckforsymbols.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@host_to_ieee.pro
@fits_ascii_encode.pro
@savgol.pro
@la_cosmic.pro
@find.pro
@fxread.pro
@fxhread.pro
@ieee_to_host.pro
@stddev.pro

function func_hrsdi_pca, obj, truncate_pca, covMatrix, data

;     oobj = obj
    pca_cube = obj
    nobj = (size(obj))[3]
    dim = (size(obj))[1]
; 
; 		med = median(obj, dim=3, /double, /even)
;     for i=0,nobj-1 do obj[*,*,i] = obj[*,*,i]-med
;     
;     data = transpose(reform(obj,dim*dim,nobj))	;[0:dim-1,0:dim-1,*]
;     covMatrix = matrix_multiply(data, data, /btranspose)

    eigenval = la_eigenql(covMatrix, EIGENVECTORS=eigenvect, range=[nobj-truncate_pca-1,nobj-1], /DOUBLE)
    eigenval = reverse(eigenval)

    eigenvect = reverse(eigenvect,2)
    pc_orig = matrix_multiply(eigenvect,data,/atranspose)

    pc = pc_orig
    for k=0,truncate_pca-1 do pc[k,*] = pc_orig[k,*]/(eigenval[k])

    for i=0,nobj-1 do begin

        data2 = transpose(reform(obj[*,*,i],dim*dim))
        s1 = matrix_multiply(pc_orig, data2, /btranspose)
        sk = reform(matrix_multiply(s1, pc), dim, dim)
        
        pca_cube[*,*,i] = obj[*,*,i]-sk

    endfor

    return, pca_cube
    
end


function MUSE_get_wl, file

  hdr = headfits(file[0], exten=1, /silent)
  naxis3 = double(get_eso_keyword(hdr, 'NAXIS3'))
  crval3 = double(get_eso_keyword(hdr,'CRVAL3'))	;readout coordinate at reference pixel
  crpix3 = double(sxpar(hdr,'CRPIX3'))	;readout coord. at reference pixel
  cd3_3 = double(get_eso_keyword(hdr, 'CD3_3'))
  wl = dblarr(naxis3)
  wl[0] = crval3
  for i=0L,naxis3-2 do wl[i+1] = wl[i] + cd3_3

	return, wl
	
end


pro MUSE_HRSDI, extract=extract


;===============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/MUSE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
; selp = 1
base = tmp[selp-1]
datadir = base+'Reduced/'
resdir = base+'HRSDI/'
spawn, 'mkdir -p '+resdir

perc = 0.01
radius = 50.	;cannot be larger than 100

;===============================================================================================

if keyword_set(extract) then begin

    file = file_search(datadir+'*DATACUBE_FINAL*fits', count=nfiles)

    radius0 = 100
    radius1 = 100

    wl = MUSE_get_wl(file[0])

		hdr = headfits(file[0], exten=1, /silent)
		naxis3 = double(get_eso_keyword(hdr, 'NAXIS3'))
		hdr0 = headfits(file[0], exten=0, /silent)
		ron1 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT1 RON'))
		ron2 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT2 RON'))
		ron3 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT3 RON'))
		ron4 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT4 RON'))
		gain1 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT1 GAIN'))
		gain2 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT2 GAIN'))
		gain3 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT3 GAIN'))
		gain4 = double(get_eso_keyword(hdr0, 'HIERARCH ESO DET OUT4 GAIN'))
		
		ron = mean([ron1, ron2, ron3, ron4])
		gain = mean([gain1, gain2, gain3, gain4])

    allcubes = dblarr(2*radius1,2*radius1,naxis3,nfiles)
    cube = dblarr(2*radius1,2*radius1,naxis3)

    coordinates = dblarr(nfiles,2)
    
    ;read files, cut and recenter and rotate
    for xx=0,nfiles-1 do begin

        im = mrdfits(file[xx], 1, /silent)
        hdr0 = headfits(file[xx], exten=0, /silent)
        
        ;compute stddev for sigma determination in la_cosmic
        tmp = im[120-20:120+19,120-20:120+19,*]
        tmp2 = reform(tmp, 40.*40., n_elements(im[0,0,*]))
        sdev = stddev(tmp2, dim=1, /nan)
        
        
;         ;remove first and last 10 frames
;         nl = n_elements(im[0,0,*])
;         im = im[*,*,10:nl-1-10]
;         wl = wl[10:nl-1-10]
; 
;         sum1 = total(im, 1)
;         sum2 = total(sum1, 1)
;         idxgood = where(sum2 ne 0.)
;         idxbad = where(sum2 eq 0.)
;         im = im[*,*,idxgood]
;         wl = wl[idxgood]
;         save, idxgood, idxbad, filename=resdir+'frame_idx.sav'
        
        aveim = median(im,dim=3,/even)

        sz = (size(im))[1:2]
        nframes = (size(im))[3]
        
        ;window, 0, xs=3.*sz[0], ys=3.*sz[1]
        ;device, cursor_standard=2
        ;cgimage, scale_image_am(aveim), /axis

				;star identification
        ;print, 'Select Star'
        ;cursor, x, y, 3, /data  ;serves as start value
        
        find, aveim, x, y, flux, sharp, roudness, 2.*stddev(aveim,/nan), 8., [-1,1], [0.2,10]
				if (n_elements(flux) gt 1) then begin
					
						tmp = max(flux,idxf)
						x = x[idxf]
						y = y[idxf]
					
				endif
				
				coordinates[xx,*] = [x,y]
				
				window, 0, xs=1000, ys=1000
				cgimage, scale_image_am(aveim), /axis
				tvcircle, /data, 10, x, y, color=cgcolor('green'), thick=3
				
        
        xa = indgen(2*radius0) & ya = xa

        cutim = aveim[x-radius0:x+radius0-1, y-radius0:y+radius0-1]
        estimates = [median(cutim), max(cutim), 2., 2., radius0, radius0, 0];, 1]
        yfit = mpfit2dpeak(cutim, A, xa, ya, /gaussian, estimates=estimates, dof=dof, chisq=chisq, perror=perror, /tilt, parinfo=pi, bestnorm=bestnorm, /quiet);, weights=weights
        
        newcube = dblarr(2.*radius0, 2.*radius0, nframes)
        for i=0,nframes-1 do begin
        
            cutim = im[x-radius0:x+radius0-1, y-radius0:y+radius0-1,i]
            
            tbp = dblarr(2.*radius0, 2.*radius0)
            idx0 = where(finite(cutim) ne 1)

            if (n_elements(idx0) lt radius1*radius1) then begin	;because of empty frames due to the laser wavelength
            
							if (idx0[0] ne -1) then begin
            
								idx = array_indices(cutim, idx0)
								for j=0,n_elements(idx[0,*])-1 do tbp[idx[0,j],idx[1,j]] = 1
								tim = cutim
								fixpix_mod, tim, tbp, outim, npix=14, /weight, /silent
								cutim = outim
                    
							endif
                
							writefits, datadir+'tmp.fits', cutim
							la_cosmic, datadir+'tmp.fits', gain=1./gain, readn=ron, sigclip=sdev[i], niter=3, sigfrac=1.
							
							bp = mrdfits(datadir+'tmp-mask.fits', 0, /silent)
							if total(bp gt 0) then begin

								tim = cutim
								fixpix_mod, tim, bp, outim, npix=14, /weight, /silent
								cutim = outim
								bp = 0
								
							endif
							
							spawn, 'rm '+datadir+'tmp*fits'

            endif

            ;cutim = mrdfits(datadir+'tmp-out.fits', 0, /silent)
            
						stmp = shift_sub(cutim, radius0-A[4], radius0-A[5])
						;stmp = fftshift(cutim, radius0-A[4], radius0-A[5])
						cstmp = stmp[radius0-radius1:radius0+radius1-1,radius0-radius1:radius0+radius1-1]
						allcubes[*,*,i,xx] = cstmp
            
            proceeding_text,loop=nframes, i=i, prompt='> Frame   '+string(i+1,form='(I4)')
            
        endfor
        
    endfor

    if (nfiles gt 1) then cubes = mean(allcubes, dim=4) else cubes = temporary(allcubes)

;     ;remove first and last 10 frames
;     nl = n_elements(cubes[0,0,*])
;     cubes = cubes[*,*,10:nl-1-10]
;     wl = wl[10:nl-1-10]

    sum1 = total(cubes, 1, /nan)
    sum2 = total(sum1, 1, /nan)
    idxgood = where(sum2 ne 0.)
    idxgood = idxgood[10:n_elements(idxgood)-10-1]	;remove first and last 10 frames as well
    idxbad = where(sum2 eq 0.)
    idxbad = [indgen(10),idxbad,indgen(10)+n_elements(cubes[0,0,*])-10]
    idxbad = idxbad[(uniq(idxbad, sort(idxbad)))]
    
    ;cube = temporary(cubes[*,*,idxgood])
    ;wl = wl[idxgood]

    save, idxgood, idxbad, filename=resdir+'frame_idx.sav'
    
    writefits, resdir+'FinalCube.fits', cubes, hdr
    writefits, resdir+'FinalCube.fits', wl, /append

    cube = 0
    
endif
;end extract

;-------------------------------------------------------------------------

print, 'Reading data'

restore, resdir+'frame_idx.sav'

cube = mrdfits(resdir+'FinalCube.fits', 0, hdr, /silent)
wl = mrdfits(resdir+'FinalCube.fits', 1, /silent)
cube = cube[*,*,idxgood]
wl = wl[idxgood]

; idxwl = where(wl ge 6000. and wl le 6800.)
; cube = cube[*,*,idxwl]
; wl = wl[idxwl]

nframes = n_elements(wl)
dim = (size(cube))[1]


for i=0,nframes-1 do cube[dim/2.,dim/2.,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.,dim/2.+1,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.+1,dim/2.,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.,dim/2.-1,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.-1,dim/2.,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.-1,dim/2.+1,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.+1,dim/2.+1,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.+1,dim/2.-1,i] = median(cube[*,*,i])
for i=0,nframes-1 do cube[dim/2.-1,dim/2.-1,i] = median(cube[*,*,i])

;create a median spectrum but center is not consider because spectral form deviates significantly
;Hoeijmakers+2018 https://www.aanda.org/articles/aa/pdf/2018/09/aa32902-18.pdf
;Haffert+2019

; mask_t = shift(dist(dim), dim/2., dim/2.)
; mask = mask_t ge 1.
factor = dblarr(dim,dim)

;print, 'Normalization'

normcube = cube
for xx=0,dim-1 do begin

    for yy=0,dim-1 do begin

				factor[xx,yy] = tsum(wl, cube[xx,yy,*])
        normcube[xx,yy,*] = cube[xx,yy,*]/factor[xx,yy]
        
    endfor

    proceeding_text,loop=dim, i=xx, prompt='> Normalization   '+string(xx+1,form='(I4)')
    
endfor

;take the brightest e.g. 1% of pixels for reference
;print, 'Reference spectrum'
avecube = median(cube, dim=3, /even)
nref = perc*n_elements(avecube[*,0])*n_elements(avecube[0,*])
idx = bsort(avecube, /reverse)
idxs = array_indices(avecube, idx)
idxs = idxs[*,0:nref-1]
tmpspec = dblarr(nframes, nref)
for i=0,nref-1 do tmpspec[*,i] = normcube[idxs[0,i],idxs[1,i],*]
refspec = median(tmpspec, dim=2, /even)

; for i=0,nframes-1 do normcube[*,*,i] = normcube[*,*,i]*mask

;make refspec around PSF instead taking entire image
;med1 = median(normcube, dim=1)
;refspec = median(med1, dim=1)
; refrad = 6.
; tmp = normcube[dim/2.-refrad:dim/2.+refrad-1, dim/2.-refrad:dim/2.+refrad-1,*]
; med1 = median(tmp, dim=1, /even)
; refspec = median(med1, dim=1, /even)


; window, 0
; plot, wl, refspec, /nodata
; for i=0,dim-1 do begin
;     for j=0,dim-1 do begin
;     
;         oplot, wl, normcube[i,j,*]
;     
;     endfor
; endfor
; oplot, wl, refspec, color=rgb(255,0,0)

;-------------------------------------------------------------------------

;divide the refspec and apply smoothing filter

;print, 'differential low order continuum signal'

dlocs = normcube    ;differential low order continuum signal
sdlocs = normcube

savgolfilter = savgol(50,50,0,3,/double)

for i=0,dim-1 do begin
    for j=0,dim-1 do begin
    
        dlocs[i,j,*] = normcube[i,j,*]/refspec
        ;sdlocs[i,j,*] = poly_smooth(reform(dlocs[i,j,*]),101)
        sdlocs[i,j,*] = convol(reform(dlocs[i,j,*]), savgolfilter, /edge_truncate)

    endfor
    
    proceeding_text,loop=dim, i=i, prompt='> Diff. low order   '+string(i+1,form='(I4)')

endfor


;divide each spectrum by its low-pass differential continuum to correct for low order residuals
corr_cube = sdlocs
for i=0,dim-1 do begin
    for j=0,dim-1 do begin
    
        corr_cube[i,j,*] = normcube[i,j,*]/sdlocs[i,j,*]
   
    endfor
    
    proceeding_text,loop=dim, i=i, prompt='> Diff. low order correction   '+string(i+1,form='(I4)')
    
endfor


; window, 1
; plot, wl, corr_cube, /nodata
; for i=0,dim-1 do begin
;     for j=0,dim-1 do begin
;     
;         oplot, wl, corr_cube[i,j,*]
;     
;     endfor
; endfor


;-------------------------------------------------------------------------

final_cube = corr_cube
for i=0,dim-1 do begin
    for j=0,dim-1 do begin
    
        final_cube[i,j,*] = corr_cube[i,j,*]-refspec
   
    endfor
endfor

;-------------------------------------------------------------------------

restore, resdir+'frame_idx.sav'
ntot = n_elements(idxgood)+n_elements(idxbad)
; 
; out_cube = dblarr(n_elements(final_cube[*,0,0]),n_elements(final_cube[0,*,0]),ntot)
; for i=0,n_elements(idxgood)-1 do out_cube[*,*,idxgood[i]] = final_cube[*,*,i]
; 
; writefits, resdir+'cube.fits', out_cube, hdr

;-------------------------------------------------------------------------

;PCA
print, 'PCA'


;truncate_pca = [5,7,10,12,15]
;truncate_pca = [1,2,3,4,5,10,15,20,50,75,100,500,1000,2000]
truncate_pca = [5]
npca = n_elements(truncate_pca)

;------
obj = final_cube[dim/2-radius:dim/2+radius-1, dim/2-radius:dim/2+radius-1, *]
factor = factor[dim/2-radius:dim/2+radius-1, dim/2-radius:dim/2+radius-1]

nobj = (size(obj))[3]
dim = (size(obj))[1]

med = median(obj, dim=3, /double, /even)
for i=0,nobj-1 do obj[*,*,i] = obj[*,*,i]-med
;for i=0,nobj-1 do obj[*,*,i] = obj[*,*,i]-mean(obj[*,*,i], /double)


data = transpose(reform(obj,dim*dim,nobj))	;[0:dim-1,0:dim-1,*]
covMatrix = matrix_multiply(data, data, /btranspose)

;------


for xx=0,npca-1 do begin

;     tmp = final_cube
   
;     pca_cube = func_hrsdi_pca(tmp, truncate_pca[xx])
    pca_cube = func_hrsdi_pca(obj, truncate_pca[xx], covMatrix, data)

    out_cube = dblarr(n_elements(pca_cube[*,0,0]),n_elements(pca_cube[0,*,0]),ntot)
    
    for i=0,n_elements(idxgood)-1 do out_cube[*,*,idxgood[i]] = pca_cube[*,*,i]
    ;out_cube = pca_cube
    
    for i=0,dim-1 do begin
        for j=0,dim-1 do begin
        
            out_cube[i,j,*] = out_cube[i,j,*]*factor[i,j]
        
        endfor
    endfor
    
    fn = resdir+'PCA'+strcompress(fix(truncate_pca[xx]),/rem)+'_FinalCube.fits'
    writefits, fn, out_cube, hdr
;     writefits, fn, wl, /append

;     writefits, resdir+'sumPCA'+strcompress(fix(truncate_pca[xx]),/rem)+'_FinalCube.fits', total(out_cube[*,*,1451-4:1451+4],3)

    proceeding_text,loop=npca, i=xx, prompt='> PCA   '+string(xx+1,form='(I4)')
    
endfor

stop
end



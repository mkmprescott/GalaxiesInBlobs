; remember, to call this, inside IDL type:  @LIONSscript
;.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/Sampler.pro
;.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/DensityGetter.pro
;.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/Plotter.pro 

; INPUTS FOR USER TO ALTER AS DESIRED: 
ap_radius_arcsec = 10.
nApertures = 10000.    
brightestmag = 22.
dimmestmag = 30.
binsize = 0.5
z=2.3
zerr=0.15
band='red'
magband = 'F140W'
 ; filenames 
blobFITSname =  'allblobs_10arcsec'                                  ; will be a fits file 
GOODSSfieldFITSname = 'GOODSSapertures_10arcsec_10000aps'            ; will be a fits file 
GOODSNfieldFITSname = 'GOODSNapertures_10arcsec_10000aps'            ; will be a fits file
densityname = 'densities_10arcsec'                                   ; will be both a fits file and a png, after modifications

; INITIALIZATION (this is really only to make the cut-making easier)
inputzstring = strcompress(string(z, format='(D0.3)'), /remove_all)     ; helps with having an accurate name for the file output by making z cuts  
zstring = strjoin(strsplit(inputzstring, '.', /extract), 'p')           ; helps with having an accurate name for the file output by making z cuts  





; SAMPLER.PRO 
.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/Sampler.pro
   readdata, allblobs, HSTfields, nblob, nfield, ntot            ; this gets m, b, and c for blobs and also the number of blobs 
   zarray = float(allblobs.z)
   plotblobs = WHERE(zarray EQ z)

; MAKING APERTURE FILES (raw) 
;   blobsampling, ap_radius_arcsec, blobFITSname
;   fieldsampling, ap_radius_arcsec, nApertures, GOODSSfieldFITSname, indextosample=0
;   fieldsampling, ap_radius_arcsec, nApertures, GOODSNfieldFITSname, indextosample=1

; MAKING CUTS
;   ; redshift cuts to aperture catalogs made above 
;   cuttype='z'
;   makeacut, blobFITSname, cuttype, list(z, zerr), [0,1,2]      ; make redshift cuts to blob aps, but skip the PRGs since they have no z info 
;   makeacut, GOODSSfieldFITSname, cuttype, list(z, zerr), [-1]  ; GOODS-S, no apertures skipped
;   makeacut, GOODSNfieldFITSname, cuttype, list(z, zerr), [-1]  ; GOODS-N, no apertures skipped
;  ; simple color cuts; format is: makeacut, filename, 'colorline', list(m, b, c, 'useblue(either blue or altblue)', 'nameadd'), skip 
;  cuttype = 'colorline'
;    ; PRG1 
;  makeacut, blobFITSname, cuttype, list(allblobs.m[0], allblobs.b[0], allblobs.cap[0], 'altblue', 'PRG1colors'), indgen(nblob-1)+1
;  makeacut, GOODSSfieldFITSname, cuttype, list(allblobs.m[0], allblobs.b[0], allblobs.cap[0], 'altblue', 'PRG1colors'), [-1]
;  makeacut, GOODSNfieldFITSname, cuttype, list(allblobs.m[0], allblobs.b[0], allblobs.cap[0], 'altblue', 'PRG1colors'), [-1]
;    ; PRG2 
;  makeacut, blobFITSname, cuttype, list(allblobs.m[1], allblobs.b[1], allblobs.cap[1], 'blue', 'PRG2colors'), [0, indgen(nblob-2)+2]
;  makeacut, GOODSSfieldFITSname, cuttype, list(allblobs.m[1], allblobs.b[1], allblobs.cap[1], 'blue', 'PRG2colors'), [-1]
;  makeacut, GOODSNfieldFITSname, cuttype, list(allblobs.m[1], allblobs.b[1], allblobs.cap[1], 'blue', 'PRG2colors'), [-1]
;    ; PRG3 
;  makeacut, blobFITSname, cuttype, list(allblobs.m[2], allblobs.b[2], allblobs.cap[2], 'altblue', 'PRG3colors'), [0, 1, indgen(nblob-3)+3]
;  makeacut, GOODSSfieldFITSname, cuttype, list(allblobs.m[2], allblobs.b[2], allblobs.cap[2], 'altblue', 'PRG3colors'), [-1]
;  makeacut, GOODSNfieldFITSname, cuttype, list(allblobs.m[2], allblobs.b[2], allblobs.cap[2], 'altblue', 'PRG3colors'), [-1]





; DENSITYGETTER.PRO
.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/DensityGetter.pro

; COUNTING GALAXIES
;   ; get number of galaxies in each aperture for blobs and fields 
;     ; raw
;   countgals, blobFITSname, brightestmag, dimmestmag, binsize, 'galcounts_'+blobFITSname
;   countgals, GOODSSfieldFITSname, brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSSfieldFITSname
;   countgals, GOODSNfieldFITSname, brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSNfieldFITSname
;     ; z cuts 
;   countgals, 'z/'+zstring+'/'+blobFITSname+'_z'+zstring, brightestmag, dimmestmag, binsize, 'galcounts_'+blobFITSname+'_z'+zstring
;   countgals, 'z/'+zstring+'/'+GOODSSfieldFITSname+'_z'+zstring, brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSSfieldFITSname+'_z'+zstring
;   countgals, 'z/'+zstring+'/'+GOODSNfieldFITSname+'_z'+zstring, brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSNfieldFITSname+'_z'+zstring
;     ; PRG1 color cuts
;   countgals, 'colorline/PRG1colors/'+blobFITSname+'_PRG1colors', brightestmag, dimmestmag, binsize, 'galcounts_'+blobFITSname+'_PRG1colors'
;   countgals, 'colorline/PRG1colors/'+GOODSSfieldFITSname+'_PRG1colors', brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSSfieldFITSname+'_PRG1colors'
;   countgals, 'colorline/PRG1colors/'+GOODSNfieldFITSname+'_PRG1colors', brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSNfieldFITSname+'_PRG1colors'
;     ; PRG2 color cuts 
;   countgals, 'colorline/PRG2colors/'+blobFITSname+'_PRG2colors', brightestmag, dimmestmag, binsize, 'galcounts_'+blobFITSname+'_PRG2colors'
;   countgals, 'colorline/PRG2colors/'+GOODSSfieldFITSname+'_PRG2colors', brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSSfieldFITSname+'_PRG2colors'
;   countgals, 'colorline/PRG2colors/'+GOODSNfieldFITSname+'_PRG2colors', brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSNfieldFITSname+'_PRG2colors'
;     ; PRG3 color cuts 
;   countgals, 'colorline/PRG3colors/'+blobFITSname+'_PRG3colors', brightestmag, dimmestmag, binsize, 'galcounts_'+blobFITSname+'_PRG3colors'
;   countgals, 'colorline/PRG3colors/'+GOODSSfieldFITSname+'_PRG3colors', brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSSfieldFITSname+'_PRG3colors'
;   countgals, 'colorline/PRG3colors/'+GOODSNfieldFITSname+'_PRG3colors', brightestmag, dimmestmag, binsize, 'galcounts_'+GOODSNfieldFITSname+'_PRG3colors'

; GETTING DENSITIES
;     ; raw
;   density, blobFITSname, [GOODSSfieldFITSname, GOODSNfieldFITSname], ap_radius_arcsec, brightestmag, dimmestmag, binsize, densityname, band=band 
;     ; z cuts 
;   density, 'z/'+zstring+'/'+blobFITSname+'_z'+zstring, ['z/'+zstring+'/'+GOODSSfieldFITSname+'_z'+zstring, 'z/'+zstring+'/'+GOODSNfieldFITSname+'_z'+zstring], ap_radius_arcsec, brightestmag, dimmestmag, binsize, densityname+'_z'+zstring, band=band 
;     ; PRG1 color cuts
;   density, 'colorline/PRG1colors/'+blobFITSname+'_PRG1colors', ['colorline/PRG1colors/'+GOODSSfieldFITSname+'_PRG1colors', 'colorline/PRG1colors/'+GOODSNfieldFITSname+'_PRG1colors'], ap_radius_arcsec, brightestmag, dimmestmag, binsize, densityname+'_PRG1colors', band=band 
;     ; PRG2 color cuts 
;   density, 'colorline/PRG2colors/'+blobFITSname+'_PRG2colors', ['colorline/PRG2colors/'+GOODSSfieldFITSname+'_PRG2colors', 'colorline/PRG2colors/'+GOODSNfieldFITSname+'_PRG2colors'], ap_radius_arcsec, brightestmag, dimmestmag, binsize, densityname+'_PRG2colors', band=band 
;     ; PRG3 color cuts 
;   density, 'colorline/PRG3colors/'+blobFITSname+'_PRG3colors', ['colorline/PRG3colors/'+GOODSSfieldFITSname+'_PRG3colors', 'colorline/PRG3colors/'+GOODSNfieldFITSname+'_PRG3colors'], ap_radius_arcsec, brightestmag, dimmestmag, binsize, densityname+'_PRG3colors', band=band 





; PLOTTER.PRO 
.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/Plotter.pro 

; MAKING HISTOGRAMS
   ; make histograms with countgals output files
;     ; raw 
;   aphist, 'galcounts_'+blobFITSname, 'galcounts_'+GOODSSfieldFITSname, 'GOODS-S apertures', 'ngals_'+GOODSSfieldFITSname, windownumber=1      ;, autosave='y or n'
;   aphist, 'galcounts_'+blobFITSname, 'galcounts_'+GOODSNfieldFITSname, 'GOODS-N apertures', 'ngals_'+GOODSNfieldFITSname, windownumber=2 
;   aphist, 'allblobs_galcounts_10arcsec', 'GOODSSaps_galcounts_10arcsec_10000aps', 'GOODS-S apertures', 'ngals_'+GOODSSfieldFITSname+'_mag25', 1, magbin=25., magband=magband     ; example
;     ; z cuts 
;   aphist, 'galcounts_'+blobFITSname+'_z'+zstring, 'galcounts_'+GOODSSfieldFITSname+'_z'+zstring, 'GOODS-S apertures', 'ngals_'+GOODSSfieldFITSname+'_z'+zstring, windownumber=3, plotblobs=plotblobs  
;   aphist, 'galcounts_'+blobFITSname+'_z'+zstring, 'galcounts_'+GOODSNfieldFITSname+'_z'+zstring, 'GOODS-N apertures', 'ngals_'+GOODSNfieldFITSname+'_z'+zstring, windownumber=4, plotblobs=plotblobs 
;     ; PRG1 color cuts 
;   aphist, 'galcounts_'+blobFITSname+'_PRG1colors', 'galcounts_'+GOODSSfieldFITSname+'_PRG1colors', 'GOODS-S apertures', 'ngals_'+GOODSSfieldFITSname+'_PRG1colors', windownumber=5, plotblobs=[0] 
;   aphist, 'galcounts_'+blobFITSname+'_PRG1colors', 'galcounts_'+GOODSNfieldFITSname+'_PRG1colors', 'GOODS-N apertures', 'ngals_'+GOODSNfieldFITSname+'_PRG1colors', windownumber=6, plotblobs=[0] 
;     ; PRG2 color cuts 
;   aphist, 'galcounts_'+blobFITSname+'_PRG2colors', 'galcounts_'+GOODSSfieldFITSname+'_PRG2colors', 'GOODS-S apertures', 'ngals_'+GOODSSfieldFITSname+'_PRG2colors', windownumber=7, plotblobs=[1]  
;   aphist, 'galcounts_'+blobFITSname+'_PRG2colors', 'galcounts_'+GOODSNfieldFITSname+'_PRG2colors', 'GOODS-N apertures', 'ngals_'+GOODSNfieldFITSname+'_PRG2colors', windownumber=8, plotblobs=[1] 
;     ; PRG3 color cuts 
;   aphist, 'galcounts_'+blobFITSname+'_PRG3colors', 'galcounts_'+GOODSSfieldFITSname+'_PRG3colors', 'GOODS-S apertures', 'ngals_'+GOODSSfieldFITSname+'_PRG3colors', windownumber=9, plotblobs=[2]  
;   aphist, 'galcounts_'+blobFITSname+'_PRG3colors', 'galcounts_'+GOODSNfieldFITSname+'_PRG3colors', 'GOODS-N apertures', 'ngals_'+GOODSNfieldFITSname+'_PRG3colors', windownumber=10, plotblobs=[2]   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;pick up here
; PLOTTING DENSITIES 
   ; plot densities with density output files
;;;   densitynameA = densityname + '_' + band    ; because density procedure alters the filename it's given 
;;;   densityplot, densitynameA, magband, 0                     ;, indicestoplot=[0,12], autosave='y'
;   densitynameB = densityname+'_z'+zstring+ '_' + band    ; because density procedure alters the filename it's given 
;   densityplot, densitynameB, magband, 0                     ;, indicestoplot=[0,12], autosave='y'


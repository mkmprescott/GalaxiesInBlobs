.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/Sampler.pro
.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/DensityGetter.pro
.com /boomerang-data/alhall/GalaxiesInBlobs/LIONS/Plotter.pro 

ap_radius_arcsec = 10.
nApertures = 10000.    
brightestmag = 22.
dimmestmag = 30.
binsize = 0.5
blobFITSname =  'allblobs_mags_z'                                           ; will be a fits file 
GOODSSfieldFITSname = 'GOODSSapertures_10arcsec_10000aps_mags_z'            ; will be a fits file 
GOODSNfieldFITSname = 'GOODSNapertures_10arcsec_10000aps_mags_z'            ; will be a fits file
densityFITSname = 'densities_ZCUT3p157'                               ; will be both a fits file and a png, after modifications
countFITSnames = ['allblobs_galcounts_10arcsec','GOODSSaps_galcounts_10arcsec_10000aps','GOODSNaps_galcounts_10arcsec_10000aps']          ; will be both fits files and pngs 
band='red'
magband = 'F140W'
cuttype = 'z'
z=3.157
zerr=0.15
zstring='3p075'
histnameA = 'GOODSS_10arcsec_10000aps_ngals'
histnameB = 'GOODSN_10arcsec_10000aps_ngals'
histnamebin = 'GOODSS_10arcsec_10000aps_ngals_mag25'

;   ; get galaxy IDs in blobs and random field apertures 
;   blobsampling, ap_radius_arcsec, blobFITSname
;   fieldsampling, ap_radius_arcsec, nApertures, GOODSSfieldFITSname
;   fieldsampling, ap_radius_arcsec, nApertures, GOODSNfieldFITSname 

;   ; get the same as above but only for z=2.3 galaxies 
;   zfilter, blobFITSname, z, zerr, zstring
;   zfilter, GOODSSfieldFITSname, z, zerr, zstring
;   zfilter, GOODSNfieldFITSname, z, zerr, zstring


;   ; get number of galaxies in each aperture for blobs and fields 
;   countgals, blobFITSname, brightestmag, dimmestmag, binsize, countFITSnames[0]
;   countgals, GOODSSfieldFITSname, brightestmag, dimmestmag, binsize, countFITSnames[1]
;   countgals, GOODSNfieldFITSname, brightestmag, dimmestmag, binsize, countFITSnames[2]


   ; calculate densities 
;;;   density, blobFITSname, [GOODSSfieldFITSname, GOODSNfieldFITSname], ap_radius_arcsec, brightestmag, dimmestmag, binsize, densityFITSname, cuttype=cuttype,key=z,band=band
;   density, blobFITSname+'_z'+zstring, [GOODSSfieldFITSname+'_z'+zstring, GOODSNfieldFITSname+'_z'+zstring], ap_radius_arcsec, brightestmag, dimmestmag, binsize, densityFITSname+'_z'+zstring, band=band 

   ; plot densities
;;;   densityFITSnameA = densityFITSname + '_' + band    ; because density procedure alters the filename it's given 
;;;   densityplot, densityFITSnameA, magband, 0                     ;, indicestoplot=[0,12], autosave='y'
;   densityFITSnameB = densityFITSname+'_z'+zstring+ '_' + band    ; because density procedure alters the filename it's given 
;   densityplot, densityFITSnameB, magband, 0                     ;, indicestoplot=[0,12], autosave='y'

   ; make histograms 
   makehist, 'allblobs_galcounts_10arcsec', 'GOODSSaps_galcounts_10arcsec_10000aps', histnameA, 1
   makehist, 'allblobs_galcounts_10arcsec', 'GOODSNaps_galcounts_10arcsec_10000aps', histnameB, 1
   makehist, 'allblobs_galcounts_10arcsec', 'GOODSSaps_galcounts_10arcsec_10000aps', histnamebin, 1, magbin=25.


; remember, to call this, inside IDL type:  @littlescript
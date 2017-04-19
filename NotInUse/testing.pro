modularize:
   dumb histograms 
   random aps!!!!!!!!!!!!!

other stuff: 
   rewrite all the headers
   try radius stuff less than 10 







old thing:    
   ; determine ahead of time whether or not to bother saving any files
   saveanyplots = ''
   READ, saveanyplots, PROMPT='Save any color-color plots, overdensity plots, statistical luminosity functions, or histograms? (y/n)'
   IF ((saveanyplots NE 'y') AND (saveanyplots NE 'n')) THEN BEGIN 
     WHILE ((saveanyplots NE 'y') AND (saveanyplots NE 'n')) DO BEGIN 
       print, 'Invalid response. Choose y or n.' 
       READ, saveanyplots, PROMPT='Save any color-color plots, overdensity plots, statistical luminosity functions, or histograms? (y/n)'
     ENDWHILE
   ENDIF





magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datared, galaxy_maglist, blob_galaxies_binned,blob_ngal_binned, blob_ngalerr_binned, densities, density_errors
magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datamid, galaxy_maglist_mid, blob_galaxies_binned_mid, blob_ngal_binned_mid, blob_ngalerr_binned_mid, densities_mid, density_errors_mid
magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datablue, galaxy_maglist_blue, blob_galaxies_binned_blue, blob_ngal_binned_blue, blob_ngalerr_binned_blue, densities_blue, density_errors_blue
magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], dataredcut, galaxy_maglist_cuts, blob_galaxies_binned_cuts, blob_ngal_binned_cuts, blob_ngalerr_binned_cuts, densities_cuts, density_errors_cuts
magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datamidcut, galaxy_maglist_cuts_mid, blob_galaxies_binned_cuts_mid, blob_ngal_binned_cuts_mid, blob_ngalerr_binned_cuts_mid, densities_cuts_mid, density_errors_cuts_mid
magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], databluecut, galaxy_maglist_cuts_blue, blob_galaxies_binned_cuts_blue, blob_ngal_binned_cuts_blue, blob_ngalerr_binned_cuts_blue, densities_cuts_blue, density_errors_cuts_blue
;
PRO magbinning, brightestmag, dimmestmag, binsize, nbins, blobra, blobdec, ap_radius, data, galaxy_maglist, blob_galaxies_binned, blob_ngal_binned, blob_ngalerr_binned, densities, density_errors
 FORWARD_FUNCTION get_galaxies_binned
 FORWARD_FUNCTION Poisson_error
  galaxy_maglist=list()
  FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
    blob_galaxies_binned = get_galaxies_binned(blobra, blobdec, ap_radius, data, mag, mag+binsize)
    galaxy_maglist.Add, blob_galaxies_binned
  ENDFOR 
  blob_ngal_binned = dblarr(nbins)
  blob_ngalerr_binned = dblarr(nbins) 
  FOR i=0, nbins-1. DO BEGIN
    bin = galaxy_maglist[i]
    IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
    blob_ngal_binned[i] = ngal 
    blob_ngalerr_binned[i] = Poisson_error(ngal) 
  ENDFOR 
  print, total(blob_ngal_binned), " total galaxies"   ; just a test to make sure it's working - this should match the previous no-bin number 
  densities = double(blob_ngal_binned)/(!dPI*ap_radius^2.) 
  density_errors = double(blob_ngalerr_binned)/(!dPI*ap_radius^2.)
END 
;     ; MAG BINNING: 
;       ; start with raw data 
;          galaxy_maglist = list()
;          galaxy_maglist_mid = list()
;          galaxy_maglist_blue = list()
;          FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
;            ; red band 
;            blob_galaxies_binned = get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius[r], datared, mag, mag+binsize)
;            galaxy_maglist.Add, blob_galaxies_binned
;            ; middle band 
;            blob_galaxies_binned_mid = get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius[r], datamid, mag, mag+binsize)
;            galaxy_maglist_mid.Add, blob_galaxies_binned_mid
;            ; blue band 
;            blob_galaxies_binned_blue = get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius[r], datablue, mag, mag+binsize)
;            galaxy_maglist_blue.Add, blob_galaxies_binned_blue
;          ENDFOR      ; Now we have lists for each catalog containing arrays of the IDs of every galaxy in each mag bin, one array per bin
;          ; set up empty arrays containing the number of galaxies in each mag bin and corresponding Poisson errors
;             blob_ngal_binned = dblarr(nbins)
;             blob_ngalerr_binned = dblarr(nbins) 
;             blob_ngal_binned_mid = dblarr(nbins)
;             blob_ngalerr_binned_mid = dblarr(nbins) 
;             blob_ngal_binned_blue = dblarr(nbins)
;             blob_ngalerr_binned_blue = dblarr(nbins) 
;          ; fill in the empty arrays, ensuring that bins with 0 galaxies have error bars of 0 as well
;          FOR i=0, nbins-1. DO BEGIN
;            bin = galaxy_maglist[i]
;            bin_mid = galaxy_maglist_mid[i]
;            bin_blue = galaxy_maglist_blue[i]
;            IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
;            IF (bin_mid[0] EQ -1) THEN ngal_mid = 0. ELSE ngal_mid = n_elements(bin_mid) 
;            IF (bin_blue[0] EQ -1) THEN ngal_blue = 0. ELSE ngal_blue = n_elements(bin_blue) 
;            blob_ngal_binned[i] = ngal 
;            blob_ngalerr_binned[i] = Poisson_error(ngal) 
;            blob_ngal_binned_mid[i] = ngal_mid 
;            blob_ngalerr_binned_mid[i] = Poisson_error(ngal_mid) 
;            blob_ngal_binned_blue[i] = ngal_blue 
;            blob_ngalerr_binned_blue[i] = Poisson_error(ngal_blue) 
;          ENDFOR 
;          print, total(blob_ngal_binned), " total galaxies"   ; just a test to make sure it's working - this should match the previous no-bin number 
;          ; Finally, compute density and errors: 
;             densities = double(blob_ngal_binned)/(!dPI*ap_radius[r]^2.) 
;             density_errors = double(blob_ngalerr_binned)/(!dPI*ap_radius[r]^2.)
;             densities_mid = double(blob_ngal_binned_mid)/(!dPI*ap_radius[r]^2.) 
;             density_errors_mid = double(blob_ngalerr_binned_mid)/(!dPI*ap_radius[r]^2.)
;             densities_blue = double(blob_ngal_binned_blue)/(!dPI*ap_radius[r]^2.) 
;             density_errors_blue = double(blob_ngalerr_binned_blue)/(!dPI*ap_radius[r]^2.)
;       ;
;       ; now do the same thing with the color cuts
;          galaxy_maglist_cuts = list()
;          galaxy_maglist_cuts_mid = list()
;          galaxy_maglist_cuts_blue = list()
;          FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
;            blob_galaxies_binned_cuts=get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius[r], dataredcut, mag, mag+binsize)
;            galaxy_maglist_cuts.Add, blob_galaxies_binned_cuts
;            blob_galaxies_binned_cuts_mid=get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius[r], datamidcut, mag, mag+binsize)
;            galaxy_maglist_cuts_mid.Add, blob_galaxies_binned_cuts_mid
;            blob_galaxies_binned_cuts_blue=get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius[r], databluecut, mag, mag+binsize)
;            galaxy_maglist_cuts_blue.Add, blob_galaxies_binned_cuts_blue
;          ENDFOR      ; Now we have lists containing arrays of the IDs of every galaxy in each mag bin, one array per bin
;          ; set up empty arrays containing the number of galaxies in each mag bin and corresponding Poisson errors
;             blob_ngal_binned_cuts = dblarr(nbins)
;             blob_ngalerr_binned_cuts = dblarr(nbins) 
;             blob_ngal_binned_cuts_mid = dblarr(nbins)
;             blob_ngalerr_binned_cuts_mid = dblarr(nbins) 
;             blob_ngal_binned_cuts_blue = dblarr(nbins)
;             blob_ngalerr_binned_cuts_blue = dblarr(nbins) 
;          ; fill in the empty arrays, ensuring that bins with 0 galaxies have error bars of 0 as well
;          FOR i=0, nbins-1. DO BEGIN
;             ; red band 
;                  bin = galaxy_maglist_cuts[i]
;                  IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
;                  blob_ngal_binned_cuts[i] = ngal 
;                  blob_ngalerr_binned_cuts[i] = Poisson_error(ngal) 
;            ; middle band 
;                  bin_mid = galaxy_maglist_cuts_mid[i]
;                  IF (bin_mid[0] EQ -1) THEN ngal_mid = 0. ELSE ngal_mid = n_elements(bin_mid) 
;                  blob_ngal_binned_cuts_mid[i] = ngal_mid 
;                  blob_ngalerr_binned_cuts_mid[i] = Poisson_error(ngal_mid) 
;            ; blue band 
;                  bin_blue = galaxy_maglist_cuts_blue[i]
;                  IF (bin_blue[0] EQ -1) THEN ngal_blue = 0. ELSE ngal_blue = n_elements(bin_blue) 
;                  blob_ngal_binned_cuts_blue[i] = ngal_blue 
;                  blob_ngalerr_binned_cuts_blue[i] = Poisson_error(ngal_blue) 
;          ENDFOR 
;          print, total(blob_ngal_binned_cuts), " total galaxies with color cut"   ; to show how many galaxies are left after color cuts
;          ; Finally, compute density and errors: 
;             densities_cuts = double(blob_ngal_binned_cuts)/(!dPI*ap_radius[r]^2.) 
;             density_errors_cuts = double(blob_ngalerr_binned_cuts)/(!dPI*ap_radius[r]^2.)
;             densities_cuts_mid = double(blob_ngal_binned_cuts_mid)/(!dPI*ap_radius[r]^2.) 
;             density_errors_cuts_mid = double(blob_ngalerr_binned_cuts_mid)/(!dPI*ap_radius[r]^2.)
;             densities_cuts_blue = double(blob_ngal_binned_cuts_blue)/(!dPI*ap_radius[r]^2.) 
;             density_errors_cuts_blue = double(blob_ngalerr_binned_cuts_blue)/(!dPI*ap_radius[r]^2.)
















; make this stuff into a procedure too

     print, 'Now sampling the field for overdensity analysis.' 

     ; SAMPLING THE FIELD with random apertures 

     ; set up the arrays which will contain all the information about the random apertures placed
     ap_ras = dblarr(nApertures)     ; this will contain RAs of all apertures
     ap_decs = dblarr(nApertures)    ; this will contain declinations of all apertures
     ; red band 
     ap_ngal_binned = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
     ap_ngalerr_binned = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
     ap_ngal_binned_cuts = dblarr(nbins,nApertures)  
     ap_ngalerr_binned_cuts = dblarr(nbins,nApertures)   
     megalist = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 
     megalist_cuts = list() 
     ; middle band 
     ap_ngal_binned_mid = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
     ap_ngalerr_binned_mid = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
     ap_ngal_binned_cuts_mid = dblarr(nbins,nApertures)  
     ap_ngalerr_binned_cuts_mid = dblarr(nbins,nApertures)   
     megalist_mid = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 
     megalist_cuts_mid = list() 
     ; blue band 
     ap_ngal_binned_blue = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
     ap_ngalerr_binned_blue = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
     ap_ngal_binned_cuts_blue = dblarr(nbins,nApertures)  
     ap_ngalerr_binned_cuts_blue = dblarr(nbins,nApertures)   
     megalist_blue = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 
     megalist_cuts_blue = list() 
     ; make an aperture image and blank image into which to insert it 
     dist_ellipse, dim, [ 2.*ap_radius_pix[r], 2.*ap_radius_pix[r] ], ap_radius_pix[r], ap_radius_pix[r], 1.,0.
     blankim = ap_radius_pix[r]+fltarr(npix_x,npix_y)  
     ap_count = 0 ; keep track of how many apertures are made in total  

     ; now place random apertures and count galaxies inside them 
     TIC        ; time this to ensure that it's efficient 
     FOR ap=0, nApertures-1 DO BEGIN




apfinder, npix_x, npix_y, ap_radius_pix[r], dim, blankim, bpm, hdr, blobra[blob], blobdec[blob], ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec 


PRO apfinder, npix_x, npix_y, ap_radius_pix, dim, blankim, bpm, hdr, blobra, blobdec, ap_count     ! ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec 
       repeatflag = 0.       ; this will ensure that apertures don't get wasted on areas with only bad pixels 
       ap_count++
       WHILE (repeatflag EQ 0.) DO BEGIN
         ; make a random aperture:
         ap_x = floor(ap_radius_pix + (double(npix_x)-2.*ap_radius_pix)*randomu(seed))
         ap_y = floor(ap_radius_pix + (double(npix_y)-2.*ap_radius_pix)*randomu(seed)) 

         ; FIND BAD PIXELS:
          ; find corners of box in which to insert aperture image 
          xlo = ap_x - floor(ap_radius_pix) 
          xhi = ap_x + floor(ap_radius_pix)
          ylo = ap_y - floor(ap_radius_pix)  
          yhi = ap_y + floor(ap_radius_pix)
          dimsize = size(dim)
          IF ((xhi-xlo) NE (dimsize[1]-1.)) THEN BEGIN 
            xhi = ap_x + floor(ap_radius_pix) - 1.
          ENDIF
          IF ((yhi-ylo) NE (dimsize[2]-1.)) THEN BEGIN 
            yhi = ap_y + floor(ap_radius_pix) - 1.
          ENDIF
          ; create new image using array arithmetic: blank image with dim inserted into it 
          newap = blankim 
          newap[xlo:xhi,ylo:yhi] = dim 
          ; ignore everything outside the random aperture to make aperture map 
          index = where(newap LT ap_radius_pix) 
          wheretomulti, newap, index, col, row, frame
          newap[*,*] = 0. 
          newap[col,row]  = 1. 
          ; multiply the aperture map by the bpm so only GOOD pixels INSIDE the aperture are left: 
          apbpm = bpm*newap          
          badpix = where(apbpm GT 0, nbadpix)   ; label the bad pixels
          ntotalpix = n_elements(apbpm)         ; get the total number of pixels in the aperture

         ; FIND VICINITY TO BLOB: 
          xyad, hdr,ap_x, ap_y, ap_ra, ap_dec      ; convert pixels into coordnates in sky
          ; find how close the aperture is to the blob:
          vicinity = sqrt(  ( (blobra - ap_ra) * cos(blobdec*!dPI/180.) )^2. + (blobdec - ap_dec)^2. )

         ; if the aperture has at least some good pixels AND isn't near the blob, then keep it; otherwise, throw it out and make a new aperture  
         IF ((nbadpix LT ntotalpix) AND (vicinity GT 2.*ap_radius)) THEN BEGIN 
           repeatflag = 1. 
         ENDIF ELSE BEGIN 
           repeatflag = 0. 
           ap_count++    ; note that another aperture is being made since this one failed to meet the criteria for a good aperture 
         ENDELSE 
       ENDWHILE
       ; once the aperture passes the above WHILE test, it can be used
END









       ; fraction of aperture that is good pixels:
       ap_area_weight = (float(ntotalpix) - float(nbadpix)) / float(ntotalpix) 
 
       ; convert pixels to sky coords and add these coords to arrays of aperture properties
       xyad, hdr,ap_x, ap_y, ap_ra, ap_dec
       ap_ras[ap] = ap_ra
       ap_decs[ap] = ap_dec 

       ; get list of galaxy IDs in each mag bin for this aperture
       ap_galaxy_maglist = list()
       ap_galaxy_maglist_cuts = list()
       ap_galaxy_maglist_mid = list()
       ap_galaxy_maglist_cuts_mid = list()
       ap_galaxy_maglist_blue = list()
       ap_galaxy_maglist_cuts_blue = list()
       FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
         ; red band
         ap_galaxies_binned=get_galaxies_binned(ap_ra, ap_dec, ap_radius[r], datared, mag, mag+binsize)
         ap_galaxies_binned_cuts = get_galaxies_binned(ap_ra, ap_dec, ap_radius[r], dataredcut, mag, mag+binsize)
         ap_galaxy_maglist.Add, ap_galaxies_binned
         ap_galaxy_maglist_cuts.Add, ap_galaxies_binned_cuts
         ; middle band
         ap_galaxies_binned_mid=get_galaxies_binned(ap_ra, ap_dec, ap_radius[r], datamid, mag, mag+binsize)
         ap_galaxies_binned_cuts_mid = get_galaxies_binned(ap_ra, ap_dec, ap_radius[r], datamidcut, mag, mag+binsize)
         ap_galaxy_maglist_mid.Add, ap_galaxies_binned_mid
         ap_galaxy_maglist_cuts_mid.Add, ap_galaxies_binned_cuts_mid
         ; blue band
         ap_galaxies_binned_blue=get_galaxies_binned(ap_ra, ap_dec, ap_radius[r], datablue, mag, mag+binsize)
         ap_galaxies_binned_cuts_blue = get_galaxies_binned(ap_ra, ap_dec, ap_radius[r], databluecut, mag, mag+binsize)
         ap_galaxy_maglist_blue.Add, ap_galaxies_binned_blue
         ap_galaxy_maglist_cuts_blue.Add, ap_galaxies_binned_cuts_blue
       ENDFOR 
       ; put the list of ID numbers per mag bin in this aperture into the mega list for all apertures:
       megalist.Add, ap_galaxy_maglist   
       megalist_cuts.Add, ap_galaxy_maglist_cuts
       megalist_mid.Add, ap_galaxy_maglist_mid   
       megalist_cuts_mid.Add, ap_galaxy_maglist_cuts_mid
       megalist_blue.Add, ap_galaxy_maglist_blue   
       megalist_cuts_blue.Add, ap_galaxy_maglist_cuts_blue

       ; get number of galaxies in each bin and error for this aperture 
       FOR i=0, nbins-1. DO BEGIN
         ; red band 
         bin = ap_galaxy_maglist[i]
         IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
         ap_ngal_binned[i,ap] = ngal 
         ap_ngalerr_binned[i,ap] = Poisson_error(ngal)
         ; middle band 
         bin_mid = ap_galaxy_maglist_mid[i]
         IF (bin_mid[0] EQ -1) THEN ngal_mid = 0. ELSE ngal_mid = n_elements(bin_mid) 
         ap_ngal_binned_mid[i,ap] = ngal_mid 
         ap_ngalerr_binned_mid[i,ap] = Poisson_error(ngal_mid)
         ; blue band 
         bin_blue = ap_galaxy_maglist_blue[i]
         IF (bin_blue[0] EQ -1) THEN ngal_blue = 0. ELSE ngal_blue = n_elements(bin_blue) 
         ap_ngal_binned_blue[i,ap] = ngal_blue 
         ap_ngalerr_binned_blue[i,ap] = Poisson_error(ngal_blue)
       ENDFOR 
       ; same as before but with color cuts 
       FOR i=0, nbins-1. DO BEGIN
         ; red band 
         bin = ap_galaxy_maglist_cuts[i]
         IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
         ap_ngal_binned_cuts[i,ap] = ngal 
         ap_ngalerr_binned_cuts[i,ap] = Poisson_error(ngal)
         ; middle band 
         bin_mid = ap_galaxy_maglist_cuts_mid[i]
         IF (bin_mid[0] EQ -1) THEN ngal_mid = 0. ELSE ngal_mid = n_elements(bin_mid) 
         ap_ngal_binned_cuts_mid[i,ap] = ngal_mid 
         ap_ngalerr_binned_cuts_mid[i,ap] = Poisson_error(ngal_mid)
         ; blue band 
         bin_blue = ap_galaxy_maglist_cuts_blue[i]
         IF (bin_blue[0] EQ -1) THEN ngal_blue = 0. ELSE ngal_blue = n_elements(bin_blue) 
         ap_ngal_binned_cuts_blue[i,ap] = ngal_blue 
         ap_ngalerr_binned_cuts_blue[i,ap] = Poisson_error(ngal_blue)
       ENDFOR 

     ENDFOR    ; all the random apertures 
     print, ap_count ; How many total apertures were generated? 
     TOC             ; How long did this analysis take? 

; end of stuff for new procedure 








































































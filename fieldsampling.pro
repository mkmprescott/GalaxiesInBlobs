PRO fieldsampling
;-------------------------------------------------------------------------------------------------------------------------------------------------;
; fieldsampling procedure
;     Places a user-specified number of random apertures on a field (GOODS-S from 3D-HST) and keeps track of the galaxies in each aperture.
;  <describe output structure here>***********************************************************************************************************************************
;
; Uses the following:
;   FUNCTION: get_galaxies_binned
;   PROCEDURE: struct_replace_field
;   PROCEDURE: apfinder
;
;
;
;-------------------------------------------------------------------------------------------------------------------------------------------------;

 ; INITIAL SETUP 

   ; notify IDL of all functions used in this procedure
      FORWARD_FUNCTION get_galaxies_binned

   ; USER-DEFINED VARIABLES: 
      ap_radius_arcsec = 10.                             ; aperture radius in arcsec
      binsize = 0.5                                      ; size of magnitude bins 
      brightestmag = 22.                                 ; smallest (brightest) magnitude of galaxies we can reasonably see in our images 
      dimmestmag = 30.                                   ; biggest (faintest) magnitude of galaxies we can reasonably see in our images 
      nApertures = 10000.                                ; number of random apertures to use in field density calculation 
      hstpixelscale = 0.06                               ; arcseconds per pixel for the field, GOODS-S (I've also seen 0.128 "/px for this?)
      filename = 'HSTapertures_10arcsec_10000aps.fits'   ; name of the FITS file to be created, containing the results of this analysis


   ; define important quantities based on above user preferences
      ; the radius/radii of the aperture(s) used to calculate density:  
        ap_radius = ap_radius_arcsec/3600.             ; aperture radius in decimal degrees
        hstap_radius_pix = ap_radius_arcsec / hstpixelscale
        n_radii = n_elements(ap_radius)                ; number of aperture sizes being tested
      ; mag binning stuff:
        nbins = (dimmestmag - brightestmag + binsize)/binsize     ; the number of magnitude bins as determined from our range and binsize 

   ; read in 3D-HST field data and rename fields to be useful
      goodss = read_csv('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/goodss.csv')   
      struct_replace_field, goodss, 'FIELD1', goodss.FIELD1, newtag='ID'
      struct_replace_field, goodss, 'FIELD2', goodss.FIELD2, newtag='ra'
      struct_replace_field, goodss, 'FIELD3', goodss.FIELD3, newtag='dec'
      struct_replace_field, goodss, 'FIELD4', goodss.FIELD4, newtag='redmag'
      struct_replace_field, goodss, 'FIELD5', goodss.FIELD5, newtag='midmag'
      struct_replace_field, goodss, 'FIELD6', goodss.FIELD6, newtag='bluemag'   
      struct_replace_field, goodss, 'FIELD7', goodss.FIELD7, newtag='altbluemag'
      struct_replace_field, goodss, 'FIELD8', goodss.FIELD8, newtag='z'

   ; read in HST bad pixel map for random aperture analysis  
TIC 
      print, 'reading in header...'
      hsthdr =  headfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/3DHSTimages/stacks/goodss_3dhst.v4.0.F125W_F140W_F160W_det.fits')
TOC
TIC
      print, 'reading in mask image...'
      hstmask = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/3DHSTimages/masks/GOODS-S_area.fits'
TOC
TIC
      print, 'converting mask image into bpm...'
      hstbpm = -1.*(mrdfits(hstmask,0) - 1.)      ; this is the mask image for weeding out bad pixels in 3DHST field, in the same format as our blob bpms 
TOC 
TIC 
      print, 'determining size of field image (bpm)...'
      hstsizeinfo = size(hstbpm)               ; the dimensions of the mask image in pixels 
      hstnpix_x = hstsizeinfo[1]               ; number of pixels in mask in the x direction
      hstnpix_y = hstsizeinfo[2]               ; number of pixels in mask in the y direction 
TOC 


 ; SAMPLE THE FIELD
 print, 'NOW SAMPLING THE FIELD. GOODS-S FIELD ANALYSIS:'

   ; set up a structure to contain all the random aperture information including IDs of galaxies contained in each aperture
      dummyarray = lonarr(150)       ; a dummy array that provides the framework for storing galaxy IDs, assuming no aperture contains more than 150 galaxies
      dummyarray = dummyarray - 1    ; fill the dummy array with all -1s (so that if an aperture has no galaxies, this is indicated by an array of only -1s)
      dummystruct = {aperture:0L, RA:0.0, dec:0.0, galIDs:dummyarray}     ; a dummy structure for one aperture containing aperture number, center RA & dec, and galaxy ID array 
      HSTstruct = replicate(dummystruct,nApertures)                       ; the main structure containing all random aperture information
     
   ; make an aperture image where everything inside the aperture has a value of 1 and everything outside (just the edge stuff) has a value of 0:
      dist_ellipse, hstdim, [ 2.*hstap_radius_pix[r], 2.*hstap_radius_pix[r] ], hstap_radius_pix[r], hstap_radius_pix[r], 1.,0.
      hstdim[where(hstdim LT ap_radius_pix[r], hstntotalpix)] = 1.      ; note that hstntotalpix is defined in this line as well!
      hstdim[where(hstdim GE ap_radius_pix[r])] = 0.

   ; keep track of how many apertures are made in total
      hstap_count = 0 

TIC
   ; now place random apertures and count galaxies inside them 
      FOR ap=0, nApertures-1 DO BEGIN       ; for the specified number of apertures
        ; find a good aperture to use: 
           apfinder, hstnpix_x, hstnpix_y, hstap_radius_pix[r], ap_radius[r], hstdim, hstbpm, hsthdr, hstap_count, hstap_x, hstap_y, hstnbadpix, hstntotalpix, hstap_ra, hstap_dec 
        ; calculate fraction of aperture that is good pixels:
           hstap_area_weight = (float(hstntotalpix) - float(hstnbadpix)) / float(hstntotalpix) 
        ; convert pixels to sky coords 
           xyad, hsthdr, hstap_x, hstap_y, hstap_ra, hstap_dec    ; for HST field
        ; add aperture number, RA, and dec to structure for this aperture
           HSTstruct(ap).aperture = ap+1
           HSTstruct(ap).RA = hstap_ra
           HSTstruct(ap).dec = hstap_dec
        ; create array of galaxy IDs for each aperture - this will hold all the galaxy IDs
           aparray = fltarr(1)-1.    ; for now, it just contains a single entry, a dummy -1 
        ; go through all magnitude bins and add galaxies to aparray
           FOR mag=brightestmag, dimmestmag, binsize DO BEGIN  
             hst_galaxies_binned_indices = get_galaxies_binned(hstap_ra, hstap_dec, ap_radius[r], goodss, mag, mag+binsize, goodss.(3))    ; get catalog indices of the galaxies in the aperture
             IF (hst_galaxies_binned_indices[0] GT -1.0) THEN BEGIN                 ; if there are any galaxies in this bin in this aperture, do the following
               hst_galaxies_binned = goodss.ID[hst_galaxies_binned_indices]         ; turn catalog indices into actual ID numbers
               aparray = [aparray, hst_galaxies_binned]                             ; add galaxy IDs onto aparray
             ENDIF   ; if there are any galaxies in this bin at all 
           ENDFOR     ; all magnitude bins 
        ; now get rid of the dummy -1 at the beginning of aparray
           ap_ngals = n_elements(aparray)      ; total number of galaxies in the aperture plus one (the dummy -1)
           IF (ap_ngals GT 1) THEN BEGIN       ; if there are any galaxies at all in this aperture, then do the following
             lastindex = ap_ngals - 1                  ; define the last index of aparray
             arrayforstruct = aparray[1:lastindex]     ; take the dummy -1 out of aparray in order to get just galaxy IDs for the structure
             newlastindex = lastindex-1                ; since one element was removed from the array, there's a different last index now
             HSTstruct(ap).galIDs[0:newlastindex] = arrayforstruct       ; put the galaxy IDs into the structure for this aperture
           ENDIF    ; the aperture contained at least one galaxy
      ENDFOR    ; all the random apertures 
      help, HSTstruct   ; a test to make sure the completed main structure looks ok
TOC


   ; last step: write the completed main structure to a FITS file
TIC
      mwrfits, HSTstruct, filename, /create
TOC
TIC
      testHSTstruct = mrdfits(filename,1)     ; read the newly-created FITS file back in as a test 
      help, testHSTstruct                     ; check to make sure the read-in structure looks the same as the output one
TOC


END        ; end of fieldsampling procedure
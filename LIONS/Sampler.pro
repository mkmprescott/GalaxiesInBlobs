;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; SAMPLER.PRO 
;  Contains all modules used for getting the IDs of galaxies located in apertures on blobs and in fields. 
;  CONTENTS: 
;      blobsampling-----------------------procedure   (main module for blobs) 
;      fieldsampling----------------------procedure   (main module for fields) 
;      get_galaxies-----------------------function
;      apfinder---------------------------procedure 
;      readdata---------------------------procedure
;      struct_replace_field---------------procedure
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;





PRO blobsampling, ap_radius_arcsec, filename
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; blobsampling procedure
;     Places apertures on blobs listed in a catalog (blobs.csv) and keeps track of the galaxies in each aperture. Each blob's name, RA, declination, fraction of good pixels,  
;       and the IDs of the galaxies it contains are stored in a structure which is exported as a FITS file. 
;
; INPUT: ap_radius_arcsec - the radius in arcseconds of the apertures to place on each blob  
;        filename - the name of the FITS file in which the results of this analysis will be saved (to be located in LIONScatalogs/apertures in GalaxiesInBlobs) 
; 
; OUTPUT: blobstruct - the main structure in which every blob's name, RA, declination, weight (the fraction of the aperture containing good pixels, which should be 1 for 
;           all blobs), and array full of IDs of galaxies contained within the blob are stored, created by replicating the structure for one blob as many times as there
;           are blobs; this structure is also saved as a FITS file with the name given by filename 
;
; Uses the get_galaxies function and the readdata procedure. 
;
; NOTES: This procedure gets all blob data from the file blobs.csv, located in the LIONScatalogs folder of GalaxiesInBlobs. 
;
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

 ; INITIAL SETUP 

   ; notify IDL of all functions used in this procedure
      FORWARD_FUNCTION get_galaxies

   ; define important quantities based on user input 
      ap_radius = ap_radius_arcsec/3600.             ; aperture radius in decimal degrees - this is what the code actually uses
      path = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'    ; the full path for where the output FITS file should be saved
      fullfilename = path + filename + '.fits'                                    ; add the path onto the filename given by user

   ; read in the file containing all blob information
      readdata, blobs, fields, nblobs, nfields, ntotal  

   ; set up a structure to contain all the blob aperture info including IDs of galaxies contained in each aperture
      dummyarray = lonarr(75)          ; a dummy array that provides the framework for storing galaxy IDs, assuming no blob/aperture contains more than 150 galaxies
      dummyarray = dummyarray - 1      ; fill the dummy array with all -1s (so that -1s mark last indices indicating how many galaxies are in each blob/aperture)
      magdummy = 0.*dummyarray + 99.   ; make another dummy array with all 99s to store magnitudes of the galaxies in dummyarray
      dummyfloats = 0.*dummyarray - 1. ; make another dummy array with all float -1s for galaxy redshifts
      ; create a dummy structure for one blob/aperture containing blob name, RA & dec, fraction of good pixels, galaxy ID array, galaxy mags, and galaxy redshifts
      dummystruct = {blobname:'', ra:0.0, dec:0.0, weight: 1.0, galIDs:dummyarray, galmags:magdummy, galz:dummyfloats}   
      ; replicate the dummy structure to make the full structure 
      blobstruct = replicate(dummystruct,nblobs)                        ; the main structure containing all blob information 


 ; SAMPLE THE BLOBS

   ; fill in the parts of blobstruct that we already know 
      blobstruct.blobname = blobs.blobname                                ; fill in names of blobs in the main structure 
      blobstruct.ra = blobs.blobra                                        ; fill in RAs of blobs in the main structure 
      blobstruct.dec = blobs.blobdec                                      ; fill in declinations of blobs in the main structure 

   ; loop over blobs, getting galaxies in each one 
   FOR blob=0, nblobs-1 DO BEGIN
      blobdata = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+blobs.blobinfo[blob],1)     ; read this blob's catalog into a structure
     ; count the galaxies in this blob/aperture
        blob_galaxies_indices = get_galaxies(blobstruct(blob).ra, blobstruct(blob).dec, ap_radius, blobdata)
        blob_galaxies = blobdata.ID[blob_galaxies_indices]
     ; put the blob galaxies into the main structure 
        ngals = n_elements(blob_galaxies)                                                     ; total number of galaxies in the blob/aperture
        print, blobstruct(blob).blobname, ' has ', strcompress(string(ngals)), ' galaxies'    ; tell the user how many galaxies are in/near this blob 
        lastindex = ngals - 1                                                                 ; define the last index of the blob_galaxies array 
        blobstruct(blob).galIDs[0:lastindex] = blob_galaxies                                  ; put the galaxy IDs into the structure for this aperture
        blobstruct(blob).galmags[0:lastindex] = blobdata.redmag[blob_galaxies_indices]        ; put galaxy mags into the structure for this aperture
        validz = tag_exist(blobdata, 'z')                                                     ; check to see if this blob's galaxies have redshift info
        IF (validz NE 0b) THEN BEGIN                                                          ; if we have redshift info, 
          blobstruct(blob).galz[0:lastindex] = blobdata.z[blob_galaxies_indices]                ; put galaxy redshifts into the structure for this aperture
        ENDIF                                                                                 ; this blob has known galaxy redshifts
   ENDFOR       ; all the blobs 

   ; last step: write the completed main structure to a FITS file
      mwrfits, blobstruct, fullfilename, /create

END        ; end of blobsampling procedure









PRO fieldsampling, ap_radius_arcsec, nApertures, filename
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; fieldsampling procedure
;     Places a user-specified number of random apertures on a field (from 3D-HST, specified by user) and keeps track of the galaxies in each aperture. Each aperture's 
;       number, RA, declination, fraction of good pixels, and the IDs of the galaxies it contains are stored in a structure which is exported as a FITS file. 
;
; INPUT: ap_radius_arcsec - the radius in arcseconds of the apertures to place randomly throughout the field 
;        nApertures - the number of random apertures to place in the field 
;        filename - the name of the FITS file in which the results of this analysis will be saved (to be located in LIONScatalogs/apertures in GalaxiesInBlobs) 
; 
; OUTPUT: HSTstruct - the main structure in which every good aperture's number (1 through nApertures), RA, declination, weight (the fraction of the aperture containing 
;           good pixels), and array full of IDs of galaxies contained within the aperture are stored, created by replicating the structure for one aperture "nApertures" 
;           times; this structure is also saved as a FITS file with the name given by filename 
;
; Uses the get_galaxies function, the readdata procedure, and the apfinder procedure. 
;
; NOTES: This procedure gets all field data from the file fields.csv, located in the LIONScatalogs folder of GalaxiesInBlobs. 
;
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

 ; INITIAL SETUP 

   ; notify IDL of all functions used in this procedure
      FORWARD_FUNCTION get_galaxies 

   readdata, blobs, fields, nblobs, nfields, ntotal      ; read in all field data and rename structure's fields to be useful
   blobcoords = dblarr(2,nblobs)                         ; make an array for coordinates (RA & declination) of regions we want to avoid in our analysis (ie, blobs) 
   blobcoords[0,*] = blobs.blobra                        ; fill in RAs of avoidance regions 
   blobcoords[1,*] = blobs.blobdec                       ; fill in declinations of avoidance regions 

   ; have user pick which field to sample
      lastfieldindex = n_elements(fields.fieldname)-1
      FOR i=0, lastfieldindex DO BEGIN
        print, fix(i), '     ', fields.fieldname(i)        ; print the name of each field and a corresponding integer 
      ENDFOR   ; all the fields 
      sample = 0S 
      read, sample, prompt='Enter number of field you wish to sample from list above: '
      typecode = size(sample, /type)         ; a type code of 2 is an integer, which is what's needed 
      IF (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((sample LT 0) OR (sample GT lastfieldindex)) )   ) THEN BEGIN
        WHILE (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((sample LT 0) OR (sample GT lastfieldindex)) )   ) DO BEGIN 
          print, 'Invalid response. Select from the listed field numbers. Enter response as an integer.' 
          sample = 0S        ; reset sample to be an integer 
          read, sample, prompt='Enter number of field you wish to sample from list above: ' 
          typecode = size(sample, /type)         ; a type code of 2 is an integer, which is what's needed 
        ENDWHILE   ; while the 'sample' response is invalid 
      ENDIF ELSE IF ((typecode EQ 2) AND (sample GE 0) AND (sample LE lastfieldindex)) THEN BEGIN     ; if the user entered 'sample' correctly 
        print, 'selected field: ', fields.fieldname(sample) 
      ENDIF      ; if the 'sample' response is invalid 
 
   ; define important quantities based on user input and field data 
      ; the radius of the apertures used to calculate density:  
         ap_radius = ap_radius_arcsec/3600.             ; aperture radius in decimal degrees
         hstap_radius_pix = ap_radius_arcsec / fields.pixelscale(sample)     ; aperture size in pixels using the pixel scale given in fields.csv
      ; file saving stuff: 
         path = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'    ; the full path for where the output FITS file should be saved
         fullfilename = path + filename + '.fits'                                    ; add the path onto the filename given by user

   ; read in the selected 3D-HST field's catalog into a structure
      fieldfile = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/' + fields.fieldinfo(sample)
      thisfield = mrdfits(fieldfile,1)   
   ; read in HST bad pixel map for random aperture analysis  
      print, 'reading in header...'
        stackfile = fields.stack(sample)
        hsthdr =  headfits(stackfile)
      print, 'reading in mask image and converting into bpm...'
        TIC      ; let's time this part, since it makes the computer freeze sometimes (average time is less than 1 minute; any longer, and terminal needs to be killed) 
          hstmask = fields.mask(sample)
          hstbpm = -1.*(mrdfits(hstmask,0) - 1.)      ; this is the mask image for weeding out bad pixels in 3DHST field, in the same format as our blob bpms 
        TOC      ; how long did it take to read in and convert the mask image? 
      print, 'determining size of field image (bpm)...'
        hstsizeinfo = size(hstbpm)               ; the dimensions of the mask image in pixels 
        hstnpix_x = hstsizeinfo[1]               ; number of pixels in mask in the x direction
        hstnpix_y = hstsizeinfo[2]               ; number of pixels in mask in the y direction 


 ; SAMPLE THE FIELD
   print, 'NOW SAMPLING THE FIELD.'

   ; set up a structure to contain all the random aperture information including IDs of galaxies contained in each aperture
      dummyarray = lonarr(75)          ; a dummy array that provides the framework for storing galaxy IDs, assuming no aperture contains more than 150 galaxies
      dummyarray = dummyarray - 1      ; fill the dummy array with all -1s (so that if an aperture has no galaxies, this is indicated by an array of only -1s)
      magdummy = 0.*dummyarray + 99.   ; make another dummy array with all 99s to store magnitudes of the galaxies in dummyarray
      dummyfloats = 0.*dummyarray - 1. ; make another dummy array with all float -1s for galaxy redshifts
      ; make a dummy structure for one aperture containing aperture number, center RA & dec, % of aperture containing good pixels, galaxy ID array, gal mags, and gal redshifts 
      dummystruct = {aperture:0L, ra:0.0, dec:0.0, weight:0.0, galIDs:dummyarray, galmags:magdummy, galz:dummyfloats}     
      ; replicate the dummy structure to make the full structure 
      HSTstruct = replicate(dummystruct,nApertures)                       ; the main structure containing all random aperture information
     
   ; make an aperture image where everything inside the aperture has a value of 1 and everything outside (just the edge stuff) has a value of 0:
      dist_ellipse, hstdim, [ 2.*hstap_radius_pix, 2.*hstap_radius_pix], hstap_radius_pix, hstap_radius_pix, 1.,0.
      hstdim[where(hstdim LT hstap_radius_pix, hstntotalpix)] = 1.      ; note that hstntotalpix is defined in this line as well!
      hstdim[where(hstdim GE hstap_radius_pix)] = 0.

   ; keep track of how many apertures are made in total
      hstap_count = 0 

   ; now place random apertures and count galaxies inside them 
      TIC   ; let's time this 
      FOR ap=0, nApertures-1 DO BEGIN       ; for the specified number of apertures
        ; find a good aperture to use: 
           apfinder, hstnpix_x, hstnpix_y, hstap_radius_pix, ap_radius, hstdim, hstbpm, hsthdr, hstap_count, hstap_x, hstap_y, hstnbadpix, hstntotalpix, $
             hstap_ra, hstap_dec, blobcoords=blobcoords 
        ; calculate fraction of aperture that is good pixels:
           hstap_area_weight = (float(hstntotalpix) - float(hstnbadpix)) / float(hstntotalpix) 
        ; convert pixels to sky coords 
           xyad, hsthdr, hstap_x, hstap_y, hstap_ra, hstap_dec    ; for HST field
        ; add aperture number, RA, and dec to structure for this aperture
           HSTstruct(ap).aperture = ap+1
           HSTstruct(ap).ra = hstap_ra
           HSTstruct(ap).dec = hstap_dec
           HSTstruct(ap).weight = hstap_area_weight 
        ; find the galaxies in this aperture and put them into add galaxies to aparray
           hst_galaxies_indices = get_galaxies(hstap_ra, hstap_dec, ap_radius, thisfield)    ; get catalog indices of the galaxies in the aperture
           IF (hst_galaxies_indices[0] GT -1.0) THEN BEGIN       ; if there are any galaxies in this aperture, do the following
             hst_galaxies = thisfield.ID[hst_galaxies_indices]        ; turn catalog indices into actual ID numbers
             ap_ngals = n_elements(hst_galaxies)                      ; total number of galaxies in the aperture
             lastindex = ap_ngals - 1                                 ; define the last index of aparray
             HSTstruct(ap).galIDs[0:lastindex] = hst_galaxies         ; put the galaxy IDs into the structure for this aperture
             HSTstruct(ap).galmags[0:lastindex] = thisfield.redmag[hst_galaxies_indices]    ; put galaxy mags into structure for this aperture
             HSTstruct(ap).galz[0:lastindex] = thisfield.z[hst_galaxies_indices]           ; put galaxy redshifts into the structure for this aperture
           ENDIF   ; if the aperture contained at least one galaxy  
      ENDFOR    ; all the random apertures 
      help, HSTstruct   ; a test to make sure the completed main structure looks ok
      print, strcompress(string(hstap_count)), ' TOTAL APERTURES GENERATED IN THIS FIELD.'     ; How many random apertures were generated in total during this analysis? 
      hstthrownout = (float(hstap_count) - float(nApertures))/float(hstap_count)*100.          ; percentage of generated apertures that ended up being unsuitable for analysis 
      print, strcompress(string(hstthrownout)), '% of all random apertures were discarded'     ; as a safety check, print out the above percentage 
      TOC        ; how long did this analysis take?

   ; last step: write the completed main structure to a FITS file
      mwrfits, HSTstruct, fullfilename, /create


END        ; end of fieldsampling procedure









PRO zfilter, filename, z, zerr, zstring
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; zfilter procedure
;  Reads in a FITS file created by the blobsampling or fieldsampling procedures and filters that file to contain only galaxies within a specified redshift range. 
; 
; INPUT: filename - a string giving the name of a FITS file containing galaxy IDs within multiple apertures 
;        z - the target redshift at which we want to keep galaxies 
;        zerr - the maximum difference between a galaxy's redshift and the target redshift for that galaxy to be kept 
;        zstring - the target redshift written as a string in the following format: for z = 2.3, zstring = '2p3'
;
; OUTPUT: Filters galaxies with z not within z +- zerr out of the structure from the original FITS file and writes the filtered structure into a new FITS file.  
;
; Uses no other procedures nor functions. 
;
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  apstruct = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+filename+'.fits',1) ; read the apertures file into a structure
  naps = n_elements(apstruct)                                                                              ; count the number of apertures in the structure
  FOR ap=0, naps-1 DO BEGIN                                                                                ; for each aperture in the structure: 
    badz = WHERE( (apstruct(ap).galz GT (z+zerr)) OR (apstruct(ap).galz LT (z-zerr)) )                       ; find where galaxy redshifts don't fall near target redshift
    apstruct(ap).galIDs[badz] = -1                                                                           ; remove those galaxies' IDs from the structure
    apstruct(ap).galmags[badz] = 99.                                                                         ; remove those galaxies' magnitudes from the structure 
    apstruct(ap).galz[badz] = -1.                                                                            ; remove those galaxies' redshifts from the structure
  ENDFOR                                                                                                   ; all the apertures in the structure have been filtered by z 
  newfile = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+filename+'_z' + zstring + '.fits'  ; make a new filename for the filtered structure
  mwrfits, apstruct, newfile, /create                                                                                  ; write the filtered structure to a fits file 
END       ; of zfilter procedure









  FUNCTION get_galaxies, ra, dec, radius, data                                                                                                    
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; get_galaxies function                                   
  ;  Places an aperture in a specified location within a set of data and finds which galaxies in the data catalog fall within that aperture. 
  ;                                                                                        
  ; INPUTS: ra, dec - the RA and declination of a point at which to place an aperture                                                               
  ;         radius - the radius of the aperture placed at the given coordinates                                                                     
  ;         data - a catalog of sources which contains each source's RA and declination (also magnitudes, but those aren't relevant here)          
  ;                   
  ; OUTPUT: galaxies_in_aperture - the index numbers of the galaxies in the data catalog which are located within the specified aperture     
  ;
  ; Uses no other functions nor procedures.  
  ;      
  ; NOTES: The output is found by determining the distance between the aperture's center and each source, then requiring that distance to be less 
  ;           than the radius of the aperture.
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;
    distance = sqrt(  ( (ra - data.ra) * cos(dec*!dPI/180.) )^2. + (dec - data.dec)^2. )    ; Pythagorean theorem               
    galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)                         ; the sources have to be contained within the aperture 
    RETURN, galaxies_in_aperture      
  END    









  PRO apfinder, npix_x, npix_y, ap_radius_pix, ap_radius, dim, bpm, hdr, ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec, blobcoords=blobcoords       
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;                   
  ; apfinder procedure                                                                                                                                                   
  ;   Generates random RA and declination coordinates for an aperture to be placed in a field, then checks the field's bad pixel map (bpm)
  ;     to determine how much of the aperture is filled with bad pixels. If the aperture contains only bad pixels, it is discarded and a 
  ;     new one is generated in its stead. This process continues until an aperture containing good pixels is found, at which point the aperture's 
  ;     position and ratio of good pixels to total area are saved.    
  ;
  ; INPUTS: npix_x - the total number of pixels in the field image in the x direction (ie, the x-size of the field)                    
  ;         npix_y - the total number of pixels in the field image in the y direction (ie, the y-size of the field) 
  ;         ap_radius_pix - the radius of the aperture to be generated, in pixels 
  ;         ap_radius - the radius of the aperture to be generated, in decimal degrees
  ;         dim - a small image of the aperture, where everything inside the aperture has a value of 1 and everything outside has a value of 0
  ;         bpm - the bad pixel map for the field, where bad pixels have a value of 1 and good pixels have a value of 0
  ;         hdr - the header for the field image, containing information on how pixel position translates to RA and declination 
  ;         ap_count - an integer keeping track of how many apertures are being generated (to know how many are being discarded) 
  ;         ap_x - the x coordinate of the aperture center in pixels (randomly generated by this procedure) 
  ;         ap_y - the y coordinate of the aperture center in pixels (randomly generated by this procedure) 
  ;         nbadpix - the number of bad pixels in the generated aperture (calculated by this procedure) 
  ;         ntotalpix - the total number of pixels in the generated aperture (calculated by this procedure) 
  ;         ap_ra - the RA of the aperture center in decimal degrees (calculated by this procedure using ap_x and hdr) 
  ;         ap_dec - the declination of the aperture center in decimal degrees (calculated by this procedure using ap_y and hdr) 
  ;         blobcoords - an optional keyword giving a 2D array with 2 rows, RA & dec, and columns corresponding to a number of regions to be
  ;                       avoided when generating apertures; if this keyword is set, apertures will be discarded if any of the coordinate pairs
  ;                       in blobcoords fall within ap_ra and ap_dec  
  ;
  ; OUTPUT: Finds an aperture with good pixels that falls outside all regions of avoidance. This aperture's x and y coordinates (ap_x and ap_y), 
  ;           the number of bad pixels it contains (nbadpix), the total number of pixels it contains (ntotalpix), and the number of tries the 
  ;           procedure took to generate a good aperture (ap_count) are saved and returned. 
  ;                 
  ; Uses no other procedures nor functions. 
  ;                      
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;            
    ; SETUP
       ; We only want apertures with good pixels, so need to set up a WHILE loop that generates apertures until a good one comes up        
       repeatflag = 0.       ; this will ensure that apertures don't get wasted on areas with only bad pixels                                                                
       ap_count++            ; start counting the total number of apertures (including bad ones) generated        

    ; Make random apertures until a usable one comes up          
    WHILE (repeatflag EQ 0.) DO BEGIN 
                                                                                                                                   
      ; make a random aperture:                                                                                                                                           
         seed = !NULL      ; reset random seed to prevent corrupting random number sequence
         ap_x = floor(ap_radius_pix + (double(npix_x)-2.*ap_radius_pix)*randomu(seed))          ; generate a random x coordinate for aperture center         
         seed = !NULL      ; reset random seed to prevent corrupting random number sequence
         ap_y = floor(ap_radius_pix + (double(npix_y)-2.*ap_radius_pix)*randomu(seed))          ; generate a random y coordinate for aperture center      

      ; FIND BAD PIXELS:                     
         ; find corners of box into which to insert aperture image ("dim")      
            xlo = ap_x - floor(ap_radius_pix)      ; left edge         
            xhi = ap_x + floor(ap_radius_pix)      ; right edge          
            ylo = ap_y - floor(ap_radius_pix)      ; bottom edge         
            yhi = ap_y + floor(ap_radius_pix)      ; top edge   
         ; make sure the box dimensions match the dimensions of the aperture image (this gets weird because of fractions of pixels       
            dimsize = size(dim)               
            IF ((xhi-xlo) NE (dimsize[1]-1.)) THEN BEGIN    ; if box & dim don't match in x size           
              xhi = ap_x + floor(ap_radius_pix) - 1.           
            ENDIF                                       
            IF ((yhi-ylo) NE (dimsize[2]-1.)) THEN BEGIN    ; if box & dim don't match in y size          
              yhi = ap_y + floor(ap_radius_pix) - 1.                                                                                                                          
            ENDIF                                                                                                                                                             
         ; cut out the piece of the bpm that contains this random aperture                                                                                                  
            minibpm = bpm[xlo:xhi,ylo:yhi]                                                                                                                                    
         ; multiply the aperture map by the bpm so only BAD pixels INSIDE the aperture are left (since we want to count them):                                              
            apbpm = minibpm*dim                                                                                                                                               
            badpix = where(apbpm GT 0, nbadpix)   ; label the bad pixels                                                                                                       
                                                                                                                                                                         
      ; FIND VICINITY TO REGIONS OF AVOIDANCE (ie, blobs in the field):           
         IF (n_elements(blobcoords) NE 0) THEN BEGIN      ; if there are any specified regions to avoid, do the following                                               
           xyad, hdr,ap_x, ap_y, ap_ra, ap_dec            ; convert image coordinates into coordnates in sky                     
           ; find how close the aperture is to any regions we want to avoid:       
              nblobs = n_elements(blobcoords[0,*])     ; number of RAs provided tells us the number of regions (blobs) we want to avoid 
              vicinity = fltarr(nblobs)                ; make an array of distances between aperture center and avoidance regions
              FOR blob=0, nblobs-1 DO BEGIN                                                                                                             
                vicinity[blob] = sqrt(  ( (blobcoords[0,blob] - ap_ra) * cos(blobcoords[1,blob]*!dPI/180.) )^2. + (blobcoords[1,blob] - ap_dec)^2. )       
              ENDFOR    ; all the specified avoidance regions                    
            ; if we're looking in a field that doesn't contain any blobs, it doesn't matter because there's nothing to avoid 
            ; still need a vicinity array though, so:                                        
         ENDIF ELSE vicinity = [3.*ap_radius]   ; make it 1-element and give it a value bigger than what's required for an acceptable aperture      
         ; make a sort of "keyword" to indicate if aperture is within an avoidance region
            tooclose = WHERE(vicinity LE 2.*ap_radius, badlocation)     ; if aperture isn't in any avoidance regions, badlocation will be equal to 0         
        
      ; DETERMINE WHETHER OR NOT THIS APERTURE IS GOOD          
         pixratio = float(nbadpix)/float(ntotalpix)    ; determine what fraction of the aperture is made of bad pixels to report this to user       
         ; Good apertures must meet the following requirements:
         ;  1) have at least some good pixels AND 
         ;  2) not fall into a region to be avoided (ie, near a blob)                              
         IF ((nbadpix LT ntotalpix) AND (badlocation EQ 0)) THEN BEGIN      ; if the aperture meets the requirements, keep it & inform user
           print, 'KEEPING aperture ', strcompress(string(ap_count)), ' with bad/total pixel ratio of ', strcompress(string(pixratio))  
           repeatflag = 1.    ; end the loop because we've gotten our good aperture                           
         ENDIF ELSE BEGIN     ; if this aperture doesn't meet the requirements, throw it out and make a new aperture  
           ; keep a record of why this aperture was thrown out, just to make sure nothing is going amiss:       
              IF (nbadpix GE ntotalpix) THEN BEGIN     ; if this aperture only has bad pixels, that's the main reason for throwing it out   
                print, 'throw out aperture ', strcompress(string(ap_count)), ' with bad/total pixel ratio of ', strcompress(string(pixratio))  
                ; note that pixratio should always be 1 in this case
              ENDIF ELSE BEGIN       ; if the aperture had good pixels but was too close to an avoidance region, indicate that                               
                print, 'throw out aperture ', strcompress(string(ap_count)), ' for being too close to ', strcompress(string(badlocation)), ' blobs'
              ENDELSE
           ; since we threw this aperture out, we're going to need to make a new one                                                                        
              ap_count++         ; note that another aperture is being made since this one failed to meet the criteria for a good aperture                                           
              repeatflag = 0.    ; stay in this WHILE loop
         ENDELSE      ; if the aperture didn't meet all requirements and got thrown out      
                                                                                                                                                 
    ENDWHILE      ; while the repeatflag is 0, indicating that we still don't have a good aperture yet          
    ; Once the aperture passes the above WHILE test, it can be used and this procedure is done.         
                                                                                          
  END       ; end of apfinder procedure 









  PRO readdata, blobs, fields, nblobs, nfields, ntotal  
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; readdata procedure
  ;   Reads files listing blobs and fields into structures and renames the structures' fields to be useful. Used to initialize other procedures. 
  ;
  ; INPUT: none 
  ;
  ; OUTPUT: blobs - a structure containing all relevant information about all blobs in blobs.csv
  ;         fields - a structure containing all relevant information about all fields in fields.csv
  ;         nblobs - the total number of blobs in blobs.csv
  ;         nfields - the total number of fields in fields.csv
  ;         ntotal - the total number of objects in both CSV catalogs (nblobs + nfields) 
  ;
  ; Uses the struct_replace_field procedure. 
  ;
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;
 
    ; read in the file containing all blob information
       blobs = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/blobs.csv', n_table_header=1)
    ; assign names to each field of the "blobs" structure
       struct_replace_field, blobs, 'FIELD01', blobs.FIELD01, newtag='blobname'   ; name of blob 
       struct_replace_field, blobs, 'FIELD02', blobs.FIELD02, newtag='blobra'     ; RA of blob 
       struct_replace_field, blobs, 'FIELD03', blobs.FIELD03, newtag='blobdec'    ; declination of blob 
       struct_replace_field, blobs, 'FIELD04', blobs.FIELD04, newtag='stack'      ; filename of stacked image for each blob
       struct_replace_field, blobs, 'FIELD05', blobs.FIELD05, newtag='mask'       ; filename of mask image for each blob 
       struct_replace_field, blobs, 'FIELD06', blobs.FIELD06, newtag='z'          ; redshift 
       struct_replace_field, blobs, 'FIELD07', blobs.FIELD07, newtag='m'          ; slope of color-cut line 
       struct_replace_field, blobs, 'FIELD08', blobs.FIELD08, newtag='b'          ; y-intercept of color-cut line 
       struct_replace_field, blobs, 'FIELD09', blobs.FIELD09, newtag='cap'        ; mag above which all colors are consistent with redshift 
       struct_replace_field, blobs, 'FIELD10', blobs.FIELD10, newtag='bluename'   ; name of the bluest filter used to observe the blob region 
       struct_replace_field, blobs, 'FIELD11', blobs.FIELD11, newtag='midbname'   ; name of the middle filter used to observe the blob region 
       struct_replace_field, blobs, 'FIELD12', blobs.FIELD12, newtag='redname'    ; name of the reddest filter used to observe the blob region 
       struct_replace_field, blobs, 'FIELD13', blobs.FIELD13, newtag='blobinfo'   ; name of the catalog containing all the magnitude, RA, DEC, etc info of sources* 
       struct_replace_field, blobs, 'FIELD14', blobs.FIELD14, newtag='pixelscale' ; pixel scale in arceseconds per pixel ("/px) for each blob image 
       nblobs = n_elements(blobs.blobname)                                        ; number of blobs in the blobs.csv file 
       ; * note: to use, must be preceded by '/boomerang-data/alhall/GalaxiesInBlobs/LIONSCatalogs/individual/'+ 

    ; read in the file containing all field information
       fields = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/fields.csv', n_table_header=2)
    ; assign names to each field of the "fields" structure
       struct_replace_field, fields, 'FIELD1', fields.FIELD1, newtag='fieldname'  ; name of each field (for example, GOODS-S) 
       struct_replace_field, fields, 'FIELD2', fields.FIELD2, newtag='stack'      ; filename of stacked image for each field
       struct_replace_field, fields, 'FIELD3', fields.FIELD3, newtag='mask'       ; filename of mask image for each field
       struct_replace_field, fields, 'FIELD4', fields.FIELD4, newtag='pixelscale' ; pixel scale in arcseconds per pixel ("/px) for each field image 
       struct_replace_field, fields, 'FIELD5', fields.FIELD5, newtag='fieldinfo'  ; name of the catalog containing all the magnitude, RA, DEC, etc info of field sources* 
       nfields = n_elements(fields.fieldname)                                     ; number of fields in the fields.csv file 
       ; * note: to use, must be preceded by '/boomerang-data/alhall/GalaxiesInBlobs/LIONSCatalogs/individual/'+ 

    ntotal = nblobs + nfields   ; total number of blobs AND fields 

  END       ; end of readdata procedure 









  ; NOTE: THE FOLLOWING PROCEDURE WAS NOT WRITTEN BY AGNAR. HOWEVER, AGNAR DID EDIT IT.  
  ;
  PRO struct_replace_field, struct, tag, data, newtag=newtag
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; struct_replace_field procedure 
  ;   Changes the type, dimensionality, and contents of an existing structure field. The tag name may be changed in the process.
  ;
  ; INPUT: struct - the structure to be modified 
  ;        tag - a string containing the case-insensitive tag name describing structure field to modify; leading and trailing spaces will be ignored; 
  ;               if the field does not exist, the structure is not changed and an error is reported
  ;        data - data that will replace the current contents of the field specified by "tag" 
  ;        newtag - optional keyword; a string containing the new tag name for field being replaced. If not specified, the original tag name will be retained.
  ;
  ; OUTPUT: struct - the input structure, now modified 
  ;
  ; Uses no other procedures nor functions. 
  ;
  ; NOTES: 
  ;    EXAMPLES: Replacing sme.wave with the arbitrary contents of "wave" would look like this:   IDL> struct_replace_field, sme, 'wave', wave
  ;              The tag name for a field can be changed without altering the data:               IDL> struct_replace_field, clients, 'NMAE', clients.nmae, newtag='Name'
  ;    CREDIT: Valenti, 2003-Jul-20, initial coding
  ;
  ;----------------------------------------------------------------------------------------------------------------------------------------------------------------------;
    ; First, make sure user provided the right amount of input
       IF (n_params() LT 3) THEN BEGIN                                             ; if the user gave too little input 
         print, 'syntax: struct_replace_field, struct, tag, data [,newtag=]'          ; tell the user how to call this procedure properly 
         RETURN                                                                       ; end the procedure
       ENDIF                                                                       ; user called this procedure incorrectly 

    ; Check that input is a structure
       IF (size(struct, /tname) NE 'STRUCT') THEN BEGIN                            ; if the user gave input that isn't a structure
         message, 'First argument is not a structure.'                                ; give the user an error message saying so 
       ENDIF                       

    ; Get list of structure tags
       tags = tag_names(struct)         ; an array of strings containing the tag names of the input structure
       ntags = n_elements(tags)         ; the total number of tags in the input structure 

    ; Check that requested field exists in input structure.
       ctag = strupcase(strtrim(tag, 2))             ; canoncial form of tag
       itag = WHERE(tags EQ ctag, nmatch)            ; find the array index of "tags" matching the canonical form of the input tag and the number of matches
       IF (nmatch EQ 0) THEN BEGIN                   ; if the input tag doesn't actually exist in the input structure 
         message, 'Structure does not contain ' + ctag + ' field'       ; give the user an error message saying so
         RETURN                                                         ; end the procedure
       ENDIF                                         ; if the user gave a tag that isn't in the input structure 
       itag = itag[0]                                ; convert "itag" to a scalar

    ; Choose tag name for the output structure
       IF keyword_set(newtag) THEN otag = newtag ELSE otag = ctag       ; if the user gave a new tag name, name the field that, otherwise use the old name

    ; Copy any fields that precede target field, then add target field
       IF (itag EQ 0) THEN BEGIN                         ; if target field occurs first
         new = create_struct(otag, data)                     ; make a new structure where the first tag has the name "otag" and contains the input data
       ENDIF ELSE BEGIN                                  ; if there are other fields before target
         new = create_struct(tags[0], struct.(0))            ; initialize the new structure
         FOR i=1, itag-1 DO BEGIN                            ; for every tag before the target
           new = create_struct(new, tags[i], struct.(i))          ; insert leading fields unchanged
         ENDFOR                                              ; every tag before the target has been put into the structure
         new = create_struct(new, otag, data)                ; now insert new data in the specified tag 
       ENDELSE                                           ; now the new data has been put into the structure 

    ; Replicate remainder of structure after desired tag
       FOR i=itag+1, ntags-1 DO BEGIN                      ; for every tag after the target 
         new = create_struct(new, tags[i], struct.(i))         ; tack this tag onto the new structure
       ENDFOR                                              ; now "new" contains all the data the user wanted 

    ; Replace input structure with new structure
       struct = new

  END     ; of struct_replace_field procedure 

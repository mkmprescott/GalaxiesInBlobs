;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; SAMPLER.PRO 
;  Contains all modules used for getting the IDs of galaxies located in apertures on blobs and in fields. 
;  CONTENTS: 
;      blobsampling-----------------------procedure   (main module for blobs) 
;      fieldsampling----------------------procedure   (main module for fields) 
;      get_galaxies-----------------------function
;      apfinder---------------------------procedure 
;      struct_replace_field---------------procedure
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;





PRO blobsampling
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; blobsampling procedure
;     Places apertures on blobs listed in a catalog (blobs.csv) and keeps track of the galaxies in each aperture. Each blob's name, RA, declination, fraction of good pixels,  
;       and the IDs of the galaxies it contains are stored in a structure which is exported as a FITS file. 
;
; INPUT: ap_radius_arcsec - the radius in arcseconds of the apertures to place on each blob  
;        filename - the name of the FITS file in which the results of this analysis will be saved (to be located in LIONScatalogs/FITS in GalaxiesInBlobs) 
; 
; OUTPUT: blobstruct - the main structure in which every blob's name, RA, declination, weight (the fraction of the aperture containing good pixels, which should be 1 for 
;           all blobs), and array full of IDs of galaxies contained within the blob are stored, created by replicating the structure for one blob as many times as there
;           are blobs; this structure is also saved as a FITS file with the name given by filename 
;
; Uses the get_galaxies function and the struct_replace_field procedure. 
;
; NOTES: This procedure gets all blob data from the file blobs.csv, located in the LIONScatalogs folder of GalaxiesInBlobs. 
;
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

 ; INITIAL SETUP 

   ; notify IDL of all functions used in this procedure
      FORWARD_FUNCTION get_galaxies

   ; USER-DEFINED VARIABLES: 
; NOTE: maybe put these in call so user can type them in then instead of editing the code to change ap size
      ap_radius_arcsec = 10.                             ; aperture radius in arcsec
      filename = 'LABapertures_10arcsec.fits'            ; name of the FITS file to be created, containing the results of this analysis

   ; define important quantities based on above user preferences 
      ap_radius = ap_radius_arcsec/3600.             ; aperture radius in decimal degrees - this is what the code actually uses
      path = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/FITS/'    ; the full path for where the output FITS file should be saved
      filename = path + filename                                             ; add the path onto the filename given by user

   ; read in the file containing all blob information
      blobs = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/blobs.csv', n_table_header=1)
   ; rename relevant fields of the "blobs" structure (only the following fields, out of 17 total)
      struct_replace_field, blobs, 'FIELD01', blobs.FIELD01, newtag='blobname'   ; name of blob 
      struct_replace_field, blobs, 'FIELD02', blobs.FIELD02, newtag='blobra'     ; RA of blob 
      struct_replace_field, blobs, 'FIELD03', blobs.FIELD03, newtag='blobdec'    ; declination of blob 
      struct_replace_field, blobs, 'FIELD16', blobs.FIELD16, newtag='blobinfo'   ; name of catalog containing all the magnitude, RA/dec, etc info of the sources
      nblobs = n_elements(blobs.blobname)                                        ; number of blobs in the blobs.csv file 

   ; set up a structure to contain all the blob aperture info including IDs of galaxies contained in each aperture
      dummyarray = lonarr(150)       ; a dummy array that provides the framework for storing galaxy IDs, assuming no blob/aperture contains more than 150 galaxies
      dummyarray = dummyarray - 1    ; fill the dummy array with all -1s (so that -1s mark last indices indicating how many galaxies are in each blob/aperture)
      ; use the dummy array to create a dummy structure for one blob/aperture containing blob name, RA & dec, fraction of good pixels, and galaxy ID array 
      dummystruct = {blobname:'', RA:0.0, dec:0.0, weight: 1.0, galIDs:dummyarray}   
      ; replicate the dummy structure to make the full structure 
      blobstruct = replicate(dummystruct,nblobs)                        ; the main structure containing all blob information 


 ; SAMPLE THE BLOBS

   ; fill in the parts of blobstruct that we already know 
      blobstruct.blobname = blobs.blobname                                ; fill in names of blobs in the main structure 
      blobstruct.RA = blobs.blobra                                        ; fill in RAs of blobs in the main structure 
      blobstruct.dec = blobs.blobdec                                      ; fill in declinations of blobs in the main structure 

   ; loop over blobs, getting galaxies in each one 
   FOR blob=0, nblobs-1 DO BEGIN
     ; make structure for this blob's catalog and rename the relevant fields to be useful (only the first 3 out of 6 total fields)
        blobdata = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+blobs.blobinfo[blob]) 
        struct_replace_field, blobdata, 'FIELD1', blobdata.FIELD1, newtag='ID'     ; galaxy IDs 
        struct_replace_field, blobdata, 'FIELD2', blobdata.FIELD2, newtag='ra'     ; galaxy RAs
        struct_replace_field, blobdata, 'FIELD3', blobdata.FIELD3, newtag='dec'    ; galaxy declinations
     ; count the galaxies in this blob/aperture
        blob_galaxies_indices = get_galaxies(blobstruct(blob).RA, blobstruct(blob).dec, ap_radius, blobdata)
        blob_galaxies = blobdata.ID[blob_galaxies_indices]
     ; put the blob galaxies into the main structure 
        ngals = n_elements(blob_galaxies)                                                     ; total number of galaxies in the blob/aperture
        print, blobstruct(blob).blobname, ' has ', strcompress(string(ngals)), ' galaxies'    ; tell the user how many galaxies are in/near this blob 
        lastindex = ngals - 1                                      ; define the last index of the blob_galaxies array 
        blobstruct(blob).galIDs[0:lastindex] = blob_galaxies       ; put the galaxy IDs into the structure for this aperture
   ENDFOR       ; all the blobs 

   ; last step: write the completed main structure to a FITS file
      mwrfits, blobstruct, filename, /create

END        ; end of blobsampling procedure









PRO fieldsampling
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; fieldsampling procedure
;     Places a user-specified number of random apertures on a field (from 3D-HST, specified by user) and keeps track of the galaxies in each aperture. Each aperture's 
;       number, RA, declination, fraction of good pixels, and the IDs of the galaxies it contains are stored in a structure which is exported as a FITS file. 
;
; INPUT: ap_radius_arcsec - the radius in arcseconds of the apertures to place randomly throughout the field 
;        nApertures - the number of random apertures to place in the field 
;        filename - the name of the FITS file in which the results of this analysis will be saved (to be located in LIONScatalogs/FITS in GalaxiesInBlobs) 
; 
; OUTPUT: HSTstruct - the main structure in which every good aperture's number (1 through nApertures), RA, declination, weight (the fraction of the aperture containing 
;           good pixels), and array full of IDs of galaxies contained within the aperture are stored, created by replicating the structure for one aperture "nApertures" 
;           times; this structure is also saved as a FITS file with the name given by filename 
;
; Uses the get_galaxies function, the struct_replace_field procedure, and the apfinder procedure. 
;
; NOTES: This procedure gets all field data from the file fields.csv, located in the LIONScatalogs folder of GalaxiesInBlobs. 
;
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

 ; INITIAL SETUP 

   ; notify IDL of all functions used in this procedure
      FORWARD_FUNCTION get_galaxies 

   ; USER-DEFINED VARIABLES: 
      ap_radius_arcsec = 10.                             ; aperture radius in arcsec
      nApertures = 10000.                                ; number of random apertures to use in field density calculation 
      filename = 'HSTapertures_10arcsec_10000aps.fits'   ; name of the FITS file to be created, containing the results of this analysis


   ; read in all field data and rename structure's fields to be useful
      fields = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/fields.csv', n_table_header=3)
      struct_replace_field, fields, 'FIELD1', fields.FIELD1, newtag='fieldname'   ; string containing name of this field (for example, GOODS-S)
      struct_replace_field, fields, 'FIELD2', fields.FIELD2, newtag='stack'       ; string containing filename of stack image
      struct_replace_field, fields, 'FIELD3', fields.FIELD3, newtag='mask'        ; string containing filename of mask image
      struct_replace_field, fields, 'FIELD4', fields.FIELD4, newtag='pixelscale'  ; arcseconds per pixel for images
      struct_replace_field, fields, 'FIELD5', fields.FIELD5, newtag='catalog'     ; catalog containing all the sources and their info (galaxy ID, coordinates, etc) 
      ; note: to use stack and mask, they must be preceded by '/boomerang-data/alhall/GalaxiesInBlobs/Data/3DHST/3DHSTimages/'+
   ; have user pick which field to sample
      lastfieldindex = n_elements(fields)-1
      FOR i=0, lastfieldindex DO BEGIN
        print, fix(i), '     ', fields.fieldname(i)        ; print the name of each field and a corresponding integer 
      ENDFOR   ; all the fields 
      sample = 0S 
      read, sample, prompt='Enter number of field you wish to sample from list above: '
      typecode = size(sample, /type)         ; a type code of 2 is an integer, which is what's needed 
      IF (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((sample LT 0) OR (sample GT lastfieldindex)) )   ) DO BEGIN
        WHILE (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((sample LT 0) OR (sample GT lastfieldindex)) )   ) DO BEGIN 
          print, 'Invalid response. Select from the listed field numbers. Enter response as an integer.' 
          sample = 0S        ; reset sample to be an integer 
          read, sample, prompt='Enter number of field you wish to sample from list above: ' 
          typecode = size(sample, /type)         ; a type code of 2 is an integer, which is what's needed 
        ENDWHILE   ; while the 'sample' response is invalid 
      ENDIF ELSE IF ((typecode EQ 2) AND (sample GE 0) AND (sample LE lastfieldindex)) THEN BEGIN     ; if the user entered 'sample' correctly 
        print, 'selected field: ', fields.fieldname(sample) 
      ENDIF      ; if the 'sample' response is invalid 
 
   ; define important quantities based on user-defined variables and field data 
      ; the radius of the apertures used to calculate density:  
         ap_radius = ap_radius_arcsec/3600.             ; aperture radius in decimal degrees
         hstap_radius_pix = ap_radius_arcsec / fields.pixelscale(sample)     ; aperture size in pixels using the pixel scale given in fields.csv
      ; file saving stuff: 
         path = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/FITS/'    ; the full path for where the output FITS file should be saved
         filename = path + filename                                             ; add the path onto the filename given by user

   ; read in the selected 3D-HST field's data and rename fields to be useful
      fieldfile = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/' + fields.catalog(sample)
      thisfield = read_csv(fieldfile)   
      struct_replace_field, thisfield, 'FIELD1', thisfield.FIELD1, newtag='ID'
      struct_replace_field, thisfield, 'FIELD2', thisfield.FIELD2, newtag='ra'
      struct_replace_field, thisfield, 'FIELD3', thisfield.FIELD3, newtag='dec'
      struct_replace_field, thisfield, 'FIELD4', thisfield.FIELD4, newtag='redmag'
      struct_replace_field, thisfield, 'FIELD5', thisfield.FIELD5, newtag='midmag'
      struct_replace_field, thisfield, 'FIELD6', thisfield.FIELD6, newtag='bluemag'   
      struct_replace_field, thisfield, 'FIELD7', thisfield.FIELD7, newtag='altbluemag'
      struct_replace_field, thisfield, 'FIELD8', thisfield.FIELD8, newtag='z'

   ; read in HST bad pixel map for random aperture analysis  
      print, 'reading in header...'
        stackfile = '/boomerang-data/alhall/GalaxiesInBlobs/Data/3DHST/3DHSTimages/'+fields.stack(sample)
        hsthdr =  headfits(stackfile)
      print, 'reading in mask image and converting into bpm...'
        TIC      ; let's time this part, since it makes the computer freeze sometimes (average time is less than 1 minute; any longer, and terminal needs to be killed) 
          hstmask = '/boomerang-data/alhall/GalaxiesInBlobs/Data/3DHST/3DHSTimages/'+fields.mask(sample)
          hstbpm = -1.*(mrdfits(hstmask,0) - 1.)      ; this is the mask image for weeding out bad pixels in 3DHST field, in the same format as our blob bpms 
        TOC      ; how long did it take to read in and convert the mask image? 
      print, 'determining size of field image (bpm)...'
        hstsizeinfo = size(hstbpm)               ; the dimensions of the mask image in pixels 
        hstnpix_x = hstsizeinfo[1]               ; number of pixels in mask in the x direction
        hstnpix_y = hstsizeinfo[2]               ; number of pixels in mask in the y direction 


 ; SAMPLE THE FIELD
 print, 'NOW SAMPLING THE FIELD. GOODS-S FIELD ANALYSIS:'

   ; set up a structure to contain all the random aperture information including IDs of galaxies contained in each aperture
      dummyarray = lonarr(150)       ; a dummy array that provides the framework for storing galaxy IDs, assuming no aperture contains more than 150 galaxies
      dummyarray = dummyarray - 1    ; fill the dummy array with all -1s (so that if an aperture has no galaxies, this is indicated by an array of only -1s)
      ; use the dummy array to make a dummy structure for one aperture containing aperture number, center RA & dec, % of aperture containing good pixels, and galaxy ID array 
      dummystruct = {aperture:0L, RA:0.0, dec:0.0, weight:0.0, galIDs:dummyarray}     
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
           apfinder, hstnpix_x, hstnpix_y, hstap_radius_pix, ap_radius, hstdim, hstbpm, hsthdr, hstap_count, hstap_x, hstap_y, hstnbadpix, hstntotalpix, hstap_ra, hstap_dec 
        ; calculate fraction of aperture that is good pixels:
           hstap_area_weight = (float(hstntotalpix) - float(hstnbadpix)) / float(hstntotalpix) 
        ; convert pixels to sky coords 
           xyad, hsthdr, hstap_x, hstap_y, hstap_ra, hstap_dec    ; for HST field
        ; add aperture number, RA, and dec to structure for this aperture
           HSTstruct(ap).aperture = ap+1
           HSTstruct(ap).RA = hstap_ra
           HSTstruct(ap).dec = hstap_dec
           HSTstruct(ap).weight = hstap_area_weight 
        ; find the galaxies in this aperture and put them into add galaxies to aparray
           hst_galaxies_indices = get_galaxies(hstap_ra, hstap_dec, ap_radius, thisfield)    ; get catalog indices of the galaxies in the aperture
           IF (hst_galaxies_indices[0] GT -1.0) THEN BEGIN       ; if there are any galaxies in this aperture, do the following
             hst_galaxies = thisfield.ID[hst_galaxies_indices]        ; turn catalog indices into actual ID numbers
             ap_ngals = n_elements(hst_galaxies)                      ; total number of galaxies in the aperture
             lastindex = ap_ngals - 1                                 ; define the last index of aparray
             HSTstruct(ap).galIDs[0:lastindex] = arrayforstruct       ; put the galaxy IDs into the structure for this aperture
           ENDIF   ; if the aperture contained at least one galaxy  
      ENDFOR    ; all the random apertures 
      help, HSTstruct   ; a test to make sure the completed main structure looks ok
      print, strcompress(string(hstap_count)), ' TOTAL APERTURES GENERATED IN GOODS-S FIELD.'  ; How many random apertures were generated in total during this analysis? 
      hstthrownout = (float(hstap_count) - float(nApertures))/float(hstap_count)*100.          ; percentage of generated apertures that ended up being unsuitable for analysis 
      print, strcompress(string(hstthrownout)), '% of all random apertures were discarded'     ; as a safety check, print out the above percentage 
      TOC        ; how long did this analysis take?


   ; last step: write the completed main structure to a FITS file
      mwrfits, HSTstruct, filename, /create
      testHSTstruct = mrdfits(filename,1)     ; read the newly-created FITS file back in as a test 
      help, testHSTstruct                     ; check to make sure the read-in structure looks the same as the output one


END        ; end of fieldsampling procedure









  FUNCTION get_galaxies, ra, dec, radius, data                                                                                                    
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;
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
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;
    distance = sqrt(  ( (ra - data.ra) * cos(dec*!dPI/180.) )^2. + (dec - data.dec)^2. )    ; Pythagorean theorem               
    galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)                         ; the sources have to be contained within the aperture 
    RETURN, galaxies_in_aperture      
  END    









  PRO apfinder, npix_x, npix_y, ap_radius_pix, ap_radius, dim, bpm, hdr, ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec, blobcoords=blobcoords       
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     
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
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;              
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
         print, ap_x, ap_y      ; check to make sure the generated coordinates are reasonable and not repeats

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









; NOTE: THE FOLLOWING PROCEDURE WAS NOT WRITTEN BY AGNAR. 
;
pro struct_replace_field, struct, tag, data, newtag=newtag
;Change the type, dimensionality, and contents of an existing structure
; field. The tag name may be changed in the process.
;
;Inputs:
; tag (string) Case insensitive tag name describing structure field to
;  modify. Leading and trailing spaces will be ignored. If the field does
;  not exist, the structure is not changed and an error is reported.
; data (any) data that will replace current contents  of 
; [newtag=] (string) new tag name for field being replaced. If not
;  specified, the original tag name will be retained.
;
;Input/Output:
; struct (structure) structure to be modified.
;
;Examples:
;
; Replace sme.wave with the arbitrary contents of wave:
;
;   IDL> struct_replace_field, sme, 'wave', wave
;
; The tag name for a field can be changed without altering the data:
;
;   IDL> struct_replace_field, clients, 'NMAE', clients.nmae, newtag='Name'
;
;History:
; 2003-Jul-20 Valenti  Initial coding

if n_params() lt 3 then begin
  print, 'syntax: struct_replace_field, struct, tag, data [,newtag=]'
  return
endif

;Check that input is a structure.
  if size(struct, /tname) ne 'STRUCT' then begin
    message, 'first argument is not a structure'
  endif

;Get list of structure tags.
  tags = tag_names(struct)
  ntags = n_elements(tags)

;Check that requested field exists in input structure.
  ctag = strupcase(strtrim(tag, 2))		;canoncial form of tag
  itag = where(tags eq ctag, nmatch)
  if nmatch eq 0 then begin
    message, 'structure does not contain ' + ctag + ' field'
    return
  endif
  itag = itag[0]				;convert to scalar

;Choose tag name for the output structure.
  if keyword_set(newtag) then otag = newtag else otag = ctag

;Copy any fields that precede target field. Then add target field.
  if itag eq 0 then begin			;target field occurs first
    new = create_struct(otag, data)
  endif else begin				;other fields before target
    new = create_struct(tags[0], struct.(0))	;initialize structure
    for i=1, itag-1 do begin			;insert leading unchange
      new = create_struct(new, tags[i], struct.(i))
    endfor
    new = create_struct(new, otag, data)	;insert new data
  endelse

;Replicate remainder of structure after desired tag.
  for i=itag+1, ntags-1 do begin
    new = create_struct(new, tags[i], struct.(i))
  endfor

;Replace input structure with new structure.
  struct = new

end
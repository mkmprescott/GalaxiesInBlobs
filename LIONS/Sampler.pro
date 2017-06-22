;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; SAMPLER.PRO 
;  Contains all modules used for getting the IDs of galaxies located in apertures on blobs and in fields. 
;  CONTENTS: 
;      blobsampling-----------------------procedure   (main module for blobs) 
;      fieldsampling----------------------procedure   (main module for fields) 
;      makeacut---------------------------procedure   (main module for making cuts) IN PROGRESS
;      zfilter----------------------------procedure 
;      colorline--------------------------procedure 
;      colorgrid--------------------------procedure    IN PROGRESS
;      get_galaxies-----------------------function
;      apfinder---------------------------procedure 
;      readdata---------------------------procedure
;      struct_replace_field---------------procedure
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;





PRO blobsampling, ap_radius_arcsec, filename
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

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
      dummystruct = {blobname:'', ra:0.0, dec:0.0, weight: 1.0, galIDs:dummyarray, galz:dummyfloats, galmags:magdummy, mid:magdummy, blue:magdummy, altblue:magdummy}     
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
      ap_radius_pix = ap_radius_arcsec / blobs.pixelscale[blob]                                                         ; define the aperture radius in pixels for this blob image
     ; get aperture weight          
       ; get bpm and header info - but if we already have it from the previous blob, just keep that 
        stackfile = blobs.stack[blob]                      ; name of this blob's stack file that has header info in it 
        IF (blob NE 0) THEN BEGIN                            ; if this isn't the first blob in the catalog (ie, if there's a blob before this one), do the following:  
          IF (stackfile NE blobs.stack[blob-1]) THEN BEGIN     ; if this blob uses a different stack and mask than the previous one, 
            hdr =  headfits(stackfile)                           ; read in the header (to be used for converting RA & dec into pixel coordinates) 
            mask = blobs.mask[blob]                              ; name of this blob's bpm file 
            bpm = mrdfits(mask,0)                                ; read in the bpm 
            IF (blobs.bpmtype[blob] EQ 1) THEN BEGIN             ; check if bad pixels are marked with 0s
              bpm = -1.*(bpm - 1.)                                 ; if they are, flip them to get them into the right format (1 inside, 0 outside) 
            ENDIF                                                ; this blob's bpm is now in the right format 
          ENDIF                                                ; if this blob doesn't already have a header and bpm 
        ENDIF ELSE IF (blob EQ 0) THEN BEGIN                 ; if this IS the first blob in the catalog, 
          hdr =  headfits(stackfile)                           ; read in the header (to be used for converting RA & dec into pixel coordinates) 
          mask = blobs.mask[blob]                              ; name of this blob's bpm file 
          bpm = mrdfits(mask,0)                                ; read in the bpm 
          IF (blobs.bpmtype[blob] EQ 1) THEN BEGIN             ; check if bad pixels are marked with 0s
            bpm = -1.*(bpm - 1.)                                 ; if they are, flip them to get them into the right format (1 inside, 0 outside) 
          ENDIF                                                ; this blob's bpm is now in the right format 
        ENDIF                                                ; if this blob is the first in the catalog 
        ; make an aperture image where everything inside the aperture has a value of 1 and everything outside (just the edge stuff) has a value of 0:     
        dist_ellipse, dim, [ 2.*ap_radius_pix, 2.*ap_radius_pix], ap_radius_pix, ap_radius_pix, 1.,0.  ; make the aperture image
        dim[where(dim LT ap_radius_pix, ntotalpix)] = 1.                                               ; set stuff inside to 1; note that ntotalpix is defined in this line as well!
        dim[where(dim GE ap_radius_pix)] = 0.                                                          ; set stuff outside to 0        
        ; find corners of box into which to insert aperture image ("dim")  
        adxy, hdr, blobstruct(blob).ra, blobstruct(blob).dec, ap_x, ap_y       ; convert blob coordinates into image coordnates     
        xlo = ap_x - floor(ap_radius_pix)      ; left edge         
        xhi = ap_x + floor(ap_radius_pix)      ; right edge          
        ylo = ap_y - floor(ap_radius_pix)      ; bottom edge         
        yhi = ap_y + floor(ap_radius_pix)      ; top edge   
        ; make sure the box dimensions match the dimensions of the aperture image (this gets weird because of fractions of pixels) 
        dimsize = size(dim)                             ; get the dimensions of the aperture image 
        IF ((xhi-xlo) NE (dimsize[1]-1.)) THEN BEGIN    ; if box & dim don't match in x size           
          xhi = ap_x + floor(ap_radius_pix) - 1.           ; readjust box in x
        ENDIF                                       
        IF ((yhi-ylo) NE (dimsize[2]-1.)) THEN BEGIN    ; if box & dim don't match in y size          
          yhi = ap_y + floor(ap_radius_pix) - 1.          ; readjust box in y                                                                                                                 
        ENDIF   
        ; get the number of bad pixels in the blob to get the weight 
        minibpm = bpm[xlo:xhi,ylo:yhi]                                                    ; cut out the piece of the bpm that contains this blob aperture 
        apbpm = minibpm*dim                                                               ; multiply aperture map by bpm to leave only BAD pixels INSIDE blob (so we can count them) 
        badpix = where(apbpm GT 0, nbadpix)                                               ; label the bad pixels   
        blobstruct(blob).weight = (float(ntotalpix) - float(nbadpix)) / float(ntotalpix)  ; fraction of good pixels to total is aperture weight 
     ; count the galaxies in this blob/aperture
        blob_galaxies_indices = get_galaxies(blobstruct(blob).ra, blobstruct(blob).dec, ap_radius, blobdata)
        blob_galaxies = blobdata.ID[blob_galaxies_indices]
     ; put the blob galaxies into the main structure 
        ngals = n_elements(blob_galaxies)                                                     ; total number of galaxies in the blob/aperture
        print, blobstruct(blob).blobname, ' has ', strcompress(string(ngals)), ' galaxies'    ; tell the user how many galaxies are in/near this blob 
        lastindex = ngals - 1                                                                 ; define the last index of the blob_galaxies array 
        blobstruct(blob).galIDs[0:lastindex] = blob_galaxies                                  ; put the galaxy IDs into the structure for this aperture
        blobstruct(blob).galmags[0:lastindex] = blobdata.redmag[blob_galaxies_indices]        ; put red galaxy mags into the structure for this aperture
        blobstruct(blob).mid[0:lastindex] = blobdata.midmag[blob_galaxies_indices]            ; put mid galaxy mags into the structure for this aperture
        blobstruct(blob).blue[0:lastindex] = blobdata.bluemag[blob_galaxies_indices]          ; put blue galaxy mags into the structure for this aperture
        validalt = tag_exist(blobdata, 'altbluemag')                                          ; check to see if this blob's galaxies have mags in another blue filter
        IF (validalt NE 0b) THEN BEGIN                                                        ; if they do, 
          blobstruct(blob).altblue[0:lastindex] = blobdata.altbluemag[blob_galaxies_indices]   ; put altbluemags into the structure for this aperture
        ENDIF ELSE BEGIN                                                                      ; if not, 
          blobstruct(blob).altblue[0:lastindex] = blobdata.bluemag[blob_galaxies_indices]       ; just use blue mags again for altblue 
        ENDELSE                                                                               ; all the magnitude info has been filled in now 
        validz = tag_exist(blobdata, 'z')                                                     ; check to see if we have galaxy redshift info for this blob 
        IF (validz NE 0b) THEN BEGIN                                                          ; if we have redshift info, 
          blobstruct(blob).galz[0:lastindex] = blobdata.z[blob_galaxies_indices]                ; put galaxy redshifts into the structure for this aperture
        ENDIF                                                                                 ; this blob has known galaxy redshifts
   ENDFOR       ; all the blobs 

   ; last step: write the completed main structure to a FITS file
      mwrfits, blobstruct, fullfilename, /create

END        ; end of blobsampling procedure









PRO fieldsampling, ap_radius_arcsec, nApertures, filename, indextosample=indextosample
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; fieldsampling procedure
;     Places a user-specified number of random apertures on a field (from 3D-HST, specified by user) and keeps track of the galaxies in each aperture. Each aperture's 
;       number, RA, declination, fraction of good pixels, and the IDs of the galaxies it contains are stored in a structure which is exported as a FITS file. 
;
; INPUT: ap_radius_arcsec - the radius in arcseconds of the apertures to place randomly throughout the field 
;        nApertures - the number of random apertures to place in the field 
;        filename - the name of the FITS file in which the results of this analysis will be saved (to be located in LIONScatalogs/apertures in GalaxiesInBlobs) 
;        --
;        indextosample - an optional integer specifying which field in fields.csv to sample with random apertures; if not set, user will be asked which field to use 
; 
; OUTPUT: HSTstruct - the main structure in which every good aperture's number (1 through nApertures), RA, declination, weight (the fraction of the aperture containing 
;           good pixels), and array full of IDs of galaxies contained within the aperture are stored, created by replicating the structure for one aperture "nApertures" 
;           times; this structure is also saved as a FITS file with the name given by filename 
;
; Uses the get_galaxies function, the readdata procedure, and the apfinder procedure. 
;
; NOTES: - This procedure gets all field data from the file fields.csv, located in the LIONScatalogs folder of GalaxiesInBlobs. 
;        - When reading in large BPMs, the average time it takes is less than 1 minute; any longer, and terminal needs to be killed! 
;
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

 ; INITIAL SETUP 

   FORWARD_FUNCTION get_galaxies                         ; notify IDL of all functions used in this procedure
   readdata, blobs, fields, nblobs, nfields, ntotal      ; read in all field data and rename structure's fields to be useful
   blobcoords = dblarr(2,nblobs)                         ; make an array for coordinates (RA & declination) of regions we want to avoid in our analysis (ie, blobs) 
   blobcoords[0,*] = blobs.blobra                        ; fill in RAs of avoidance regions 
   blobcoords[1,*] = blobs.blobdec                       ; fill in declinations of avoidance regions 

   IF (n_elements(indextosample) EQ 1) THEN BEGIN                                                                         ; if indextosample is given and has correct # of elements,
     sample = indextosample                                                                                                 ; then we already know which field we want to sample 
   ENDIF ELSE BEGIN                                                                                                       ; otherwise, have user pick which field to sample
      lastfieldindex = n_elements(fields.fieldname)-1                                                                       ; see how many fields are in the fields.csv catalog
      FOR i=0, lastfieldindex DO BEGIN                                                                                      ; loop over all the fields in the catalog 
        print, fix(i), '     ', fields.fieldname(i)                                                                           ; print name of each field and a corresponding integer 
      ENDFOR                                                                                                                ; all the fields 
      sample = 0S                                                                                                           ; set up an empty integer for index of the selected field 
      read, sample, prompt='Enter number of field you wish to sample from list above: '                                     ; ask user which field to use 
      typecode = size(sample, /type)                                                                                        ; a type code of 2 is an integer, which is what we need 
      IF (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((sample LT 0) OR (sample GT lastfieldindex)) )   ) THEN BEGIN     ; if the user gave bad input, 
        WHILE (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((sample LT 0) OR (sample GT lastfieldindex)) )   ) DO BEGIN    ; as long as that input is still bad, 
          print, 'Invalid response. Select from the listed field numbers. Enter response as an integer.'                        ; tell the user
          sample = 0S                                                                                                           ; reset sample to be an integer 
          read, sample, prompt='Enter number of field you wish to sample from list above: '                                     ; ask the user once again which field to use 
          typecode = size(sample, /type)                                                                                        ; check the type of sample again (need 2) 
        ENDWHILE                                                                                                              ; while the 'sample' response is invalid 
      ENDIF ELSE IF ((typecode EQ 2) AND (sample GE 0) AND (sample LE lastfieldindex)) THEN BEGIN                           ; if the user entered 'sample' correctly 
        print, 'selected field: ', fields.fieldname(sample)                                                                   ; confirm the field selected by the user 
      ENDIF                                                                                                                 ; if the 'sample' response is invalid 
    ENDELSE                                                                                                               ; if user didn't specify a field to sample 
 
   ; define important quantities based on user input and field data 
      ; the radius of the apertures used to calculate density:  
         ap_radius = ap_radius_arcsec/3600.                                          ; aperture radius in decimal degrees
         ap_radius_pix = ap_radius_arcsec / fields.pixelscale(sample)                ; aperture size in pixels using the pixel scale given in fields.csv
      ; file saving stuff: 
         path = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'    ; the full path for where the output FITS file should be saved
         fullfilename = path + filename + '.fits'                                    ; add the path onto the filename given by user

   ; read in the selected 3D-HST field's catalog into a structure
      fieldfile = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/' + fields.fieldinfo(sample)
      thisfield = mrdfits(fieldfile,1)   
   ; read in HST bad pixel map for random aperture analysis  
      print, 'reading in header...' 
        stackfile = fields.stack(sample)                          ; filename of the stack image 
        hdr =  headfits(stackfile)                                ; we only need the header from the stack, not the image itself 
      print, 'reading in mask image and converting into bpm...'
        TIC                                                       ; let's time this part, since it makes the computer freeze sometimes (SEE HEADER NOTES) 
          mask = fields.mask(sample)                              ; filename of the bpm 
          bpm = -1.*(mrdfits(mask,0) - 1.)                        ; this is the mask image for weeding out bad pixels in 3DHST field, in the same format as our blob bpms 
        TOC                                                       ; how long did it take to read in and convert the mask image? 
   ; get the size of the field image's mask (bpm)
      sizeinfo = size(bpm)                                        ; the dimensions of the mask image in pixels 
      npix_x = sizeinfo[1]                                        ; number of pixels in mask in the x direction
      npix_y = sizeinfo[2]                                        ; number of pixels in mask in the y direction 


 ; SAMPLE THE FIELD
   print, 'NOW SAMPLING THE FIELD.'

   ; set up a structure to contain all the random aperture information including IDs of galaxies contained in each aperture
      dummyarray = lonarr(75)          ; a dummy array that provides the framework for storing galaxy IDs, assuming no aperture contains more than 150 galaxies
      dummyarray = dummyarray - 1      ; fill the dummy array with all -1s (so that if an aperture has no galaxies, this is indicated by an array of only -1s)
      magdummy = 0.*dummyarray + 99.   ; make another dummy array with all 99s to store magnitudes of the galaxies in dummyarray
      dummyfloats = 0.*dummyarray - 1. ; make another dummy array with all float -1s for galaxy redshifts
      ; make a dummy structure for one aperture containing aperture number, center RA & dec, % of aperture containing good pixels, galaxy ID array, gal mags, and gal redshifts 
      dummystruct = {aperture:0L, ra:0.0, dec:0.0, weight:0.0, galIDs:dummyarray, galz:dummyfloats, galmags:magdummy, mid:magdummy, blue:magdummy, altblue:magdummy}     
      ; replicate the dummy structure to make the full structure 
      HSTstruct = replicate(dummystruct,nApertures)                       ; the main structure containing all random aperture information
     
   ; make an aperture image where everything inside the aperture has a value of 1 and everything outside (just the edge stuff) has a value of 0:
      dist_ellipse, dim, [ 2.*ap_radius_pix, 2.*ap_radius_pix], ap_radius_pix, ap_radius_pix, 1.,0.
      dim[where(dim LT ap_radius_pix, ntotalpix)] = 1.      ; note that ntotalpix is defined in this line as well!
      dim[where(dim GE ap_radius_pix)] = 0.

   ; now place random apertures and count galaxies inside them 
      ap_count = 0                                                                          ; keep track of how many apertures are made in total
      TIC                                                                                   ; let's time this 
      FOR ap=0, nApertures-1 DO BEGIN                                                       ; for the specified number of apertures
        apfinder, npix_x, npix_y, ap_radius_pix, ap_radius, dim, bpm, hdr, ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec, blobcoords=blobcoords ; find a good aperture to use
        ap_area_weight = (float(ntotalpix) - float(nbadpix)) / float(ntotalpix)               ; calculate fraction of aperture that is good pixels
        xyad, hdr, ap_x, ap_y, ap_ra, ap_dec                                                  ; convert pixels to sky coords 
        HSTstruct(ap).aperture = ap+1                                                         ; add aperture number to structure for this aperture
        HSTstruct(ap).ra = ap_ra                                                              ; add RA to structure for this aperture
        HSTstruct(ap).dec = ap_dec                                                            ; add declination to structure for this aperture
        HSTstruct(ap).weight = ap_area_weight                                                 ; add weight to structure for this aperture
        galaxies_indices = get_galaxies(ap_ra, ap_dec, ap_radius, thisfield)                  ; get catalog indices of the galaxies in the aperture
        IF (galaxies_indices[0] GT -1.0) THEN BEGIN                                           ; if there are any galaxies in this aperture, add them to the structure as follows: 
          galaxies = thisfield.ID[galaxies_indices]                                             ; turn catalog indices into actual ID numbers
          ap_ngals = n_elements(galaxies)                                                       ; total number of galaxies in the aperture
          lastindex = ap_ngals - 1                                                              ; define the last index of aparray
          HSTstruct(ap).galIDs[0:lastindex] = galaxies                                          ; put the galaxy IDs into the structure for this aperture
          HSTstruct(ap).galz[0:lastindex] = thisfield.z[galaxies_indices]                       ; put galaxy redshifts into the structure for this aperture
          HSTstruct(ap).galmags[0:lastindex] = thisfield.redmag[galaxies_indices]               ; put red galaxy mags into structure for this aperture
          HSTstruct(ap).mid[0:lastindex] = thisfield.midmag[galaxies_indices]                   ; put mid galaxy mags into structure for this aperture
          HSTstruct(ap).blue[0:lastindex] = thisfield.bluemag[galaxies_indices]                 ; put blue galaxy mags into structure for this aperture
          HSTstruct(ap).altblue[0:lastindex] = thisfield.altbluemag[galaxies_indices]           ; put blue galaxy mags into structure for this aperture
        ENDIF                                                                                 ; if the aperture contained at least one galaxy  
      ENDFOR                                                                                ; all the random apertures 
      help, HSTstruct                                                                       ; a test to make sure the completed main structure looks ok
      print, strcompress(string(ap_count)), ' TOTAL APERTURES GENERATED IN THIS FIELD.'     ; How many random apertures were generated in total during this analysis? 
      thrownout = (float(ap_count) - float(nApertures))/float(ap_count)*100.                ; percentage of generated apertures that ended up being unsuitable for analysis 
      print, strcompress(string(thrownout)), '% of all random apertures were discarded'     ; as a safety check, print out the above percentage 
      TOC        ; how long did this analysis take?

   ; last step: write the completed main structure to a FITS file
      mwrfits, HSTstruct, fullfilename, /create

END        ; end of fieldsampling procedure









PRO makeacut, catalog, cuttype, parameters, skip
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; makeacut procedure
;   description 
;
; INPUT: catalog - a string giving the name of a FITS file containing galaxy IDs within multiple apertures 
;        cuttype - a string specifying which cutting procedure to use on the catalog - can be 'z' (redshift cut), 'colorline' (simple color cut), or 'colorgrid' (grid color cut) 
;        parameters - a list (meaning the IDL data type) of parameters to be passed into whichever procedure is used to make the cut; the contents of this list depend on cuttype 
;        skip - an integer or array of integers giving the index/indices of a catalog to be skipped when applying cuts; to skip no indices, set skip to -1
;
; OUTPUT: filters the catalog using one of three procedures (whichever is specified by "cuttype," see below) and outputs the filtered catalog as a FITS file 
;
; Uses the zfilter procedure, the colorline procedure, and the colorgrid procedure. 
;
; NOTES: See the zfilter, colorline, and colorgrid procedure for details regarding this procedure's "parameters" argument. 
;          For zfilter (cuttype = 'z'), parameters should be: list(z, zerr) 
;          For colorline (cuttype = 'colorline'), parameters should be: list(m, b, c, useblue, nameadd) 
;          For colorgrid (cuttype = 'colorgrid'), parameters should be: ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
;
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; read in the input structure and make sure all main arguments are valid before doing anything else 
  struct = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+catalog+'.fits',1) ; read in the catalog to see how many apertures are in it
  napertures = n_elements(struct)                                                                       ; count the number of apertures in the catalog 
  integertypes = [2, 3, 12, 13, 14, 15]                                                                 ; make an array of all the different IDL type codes for integers
  skiptype = size(skip, /type)                                                                          ; see what data type the "skip" argument is (should be integer) 
  ok = WHERE(integertypes EQ skiptype, check)                                                           ; use the above two lines to check if "skip" is actually integer(s) 
  badneg = WHERE(skip LT -1, nneg)                                                                      ; check if any elements of "skip" are negative aside from -1 
  badbig = WHERE(skip GE napertures, nbig)                                                              ; check if any elements of "skip" are larger than # of indices in catalog
  skipneg = WHERE(skip EQ -1, nminusone)                                                                ; see if any elements of "skip" are -1 (bad if skip is array, fine otherwise)
  IF (check EQ 0) THEN BEGIN                                                                                  ; if skip isn't integer(s), 
    print, 'Invalid data type for "skip" argument. This variable must be an integer. No cuts will be made.'     ; tell the user so and make no cuts 
  ENDIF ELSE IF ( (nneg NE 0) OR (nbig NE 0) ) THEN BEGIN                                                     ; if some elements of skip are too low or too high, 
    print, 'Invalid indices given in "skip" argument. No cuts will be made.'                                    ; tell the user so and make no cuts 
  ENDIF ELSE IF ( (n_elements(skip) GT 1) AND (nminusone NE 0) ) THEN BEGIN                                   ; if skip is an array but contains a -1, 
    print, 'Invalid indices given in "skip" argument. No cuts will be made.'                                    ; tell the user so and make no cuts
  ; now make sure "cuttype" argument is valid 
  ENDIF ELSE IF ( size(cuttype, /type) NE 7 ) THEN BEGIN                                                      ; if cuttype isn't a string, 
    print, 'Invalid data type for "cuttype" argument. This variable must be a string. No cuts will be made.'    ; tell the user so and make no cuts 
  ENDIF ELSE IF ( (cuttype NE 'z') AND (cuttype NE 'colorline') AND (cuttype NE 'colorgrid') ) THEN BEGIN     ; if an invalid cut type was supplied by user 
    print, 'Invalid cut type given. No cuts will be made.'                                                      ; tell the user so and make no cuts 
  ; now move on to the three different types of cuts 

  ENDIF ELSE IF (cuttype EQ 'z') THEN BEGIN                                                       ; REDSHIFT CUTS 
    IF (n_elements(parameters) NE 2) THEN BEGIN                                                          ; if wrong number of parameters given for this kind of cut, 
      print, 'Incorrect number of parameters. No cuts will be made.'                                       ;  tell the user so and make no cuts 
    ENDIF ELSE IF ( (size(parameters[0], /type) NE 4) AND (size(parameters[0], /type) NE 5) ) THEN BEGIN ; if the user's z isn't a float or double (ie, if it has wrong data type), 
      print, 'Redshift z has incorrect data type. No cuts will be made.'                                   ; tell the user so and make no cuts 
    ENDIF ELSE IF ( (size(parameters[1], /type) NE 4) AND (size(parameters[1], /type) NE 5) ) THEN BEGIN ; if the user's zerr isn't a float or double (ie, if it has wrong data type), 
      print, 'Redshift range zerr has incorrect data type. No cuts will be made.'                          ; tell the user so and make no cuts 
    ENDIF ELSE BEGIN                                                                                     ; if all the right parameters were given, 
      zfilter, catalog, parameters[0], parameters[1], skip                                                 ; run zfilter procedure on the catalog 
    ENDELSE                                                                                              ; redshift cuts are done 

  ENDIF ELSE IF (cuttype EQ 'colorline') THEN BEGIN                                               ; SIMPLE STRAIGHT-LINE COLOR CUTS
    IF (n_elements(parameters) NE 5) THEN BEGIN                                                          ; if wrong number of parameters given for this kind of cut, 
      print, 'Incorrect number of parameters. No cuts will be made.'                                       ;  tell the user so and make no cuts  
    ENDIF ELSE IF ( (size(parameters[0], /type) NE 4) AND (size(parameters[0], /type) NE 5) ) THEN BEGIN ; if the user's m isn't a float or double (ie, if it has wrong data type), 
      print, 'Line slope m has incorrect data type. No cuts will be made.'                                 ; tell the user so and make no cuts 
    ENDIF ELSE IF ( (size(parameters[1], /type) NE 4) AND (size(parameters[1], /type) NE 5) ) THEN BEGIN ; if the user's b isn't a float or double (ie, if it has wrong data type), 
      print, 'Line y-intercept b has incorrect data type. No cuts will be made.'                           ; tell the user so and make no cuts 
    ENDIF ELSE IF ( (size(parameters[2], /type) NE 4) AND (size(parameters[2], /type) NE 5) ) THEN BEGIN ; if the user's c isn't a float or double (ie, if it has wrong data type), 
      print, 'Line cap c has incorrect data type. No cuts will be made.'                                   ; tell the user so and make no cuts 
    ENDIF ELSE IF (size(parameters[3], /type) NE 7) THEN BEGIN                                           ; if useblue isn't a string (ie, if it has wrong data type), 
      print, 'Parameter "useblue" has incorrect data type. No cuts will be made.'                          ; tell the user so and make no cuts 
    ENDIF ELSE IF (size(parameters[4], /type) NE 7) THEN BEGIN                                           ; if nameadd isn't a string (ie, if it has wrong data type), 
      print, 'Parameter "nameadd" has incorrect data type. No cuts will be made.'                          ; tell the user so and make no cuts 
    ENDIF ELSE BEGIN                                                                                     ; if all the right parameters were given, 
      colorline, catalog, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], skip  ; run colorline procedure on the catalog 
    ENDELSE                                                                                              ; simple color cuts are done 

  ENDIF ELSE IF (cuttype EQ 'colorgrid') THEN BEGIN                                               ; PROBABILITY-GRID COLOR CUTS
    print, 'Color gridding has not been implemented yet. Sorry for the inconvenience.' 
  ENDIF                                                                                           ; all cuts have been made 

END       ; of makeacut procedure 









  PRO zfilter, filename, z, zerr, skip 
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; zfilter procedure
  ;  Reads in a FITS file created by the blobsampling or fieldsampling procedures and filters that file to contain only galaxies within a specified redshift range. 
  ; 
  ; INPUT: filename - a string giving the name of a FITS file containing galaxy IDs within multiple apertures 
  ;        z - the target redshift at which we want to keep galaxies 
  ;        zerr - the maximum difference between a galaxy's redshift and the target redshift for that galaxy to be kept 
  ;        skip - an integer or array of integers giving the index/indices of a catalog to be skipped when applying cuts; to skip no indices, set skip to -1 
  ;
  ; OUTPUT: Filters galaxies with z not within z +- zerr out of the structure from the original FITS file and writes the filtered structure into a new FITS file.  
  ;
  ; Uses the remover function. 
  ;
  ; NOTE: This procedure is meant to be run by the makeacut procedure rather than being called on its own. 
  ;
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
    FORWARD_FUNCTION remover                                                                                 ; notify IDL that this function will be used 
    apstruct = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+filename+'.fits',1) ; read the apertures file into a structure
    naps = n_elements(apstruct)                                                                              ; count the number of apertures in the structure
    skiporder = sort(skip)                                                                                   ; get indices of skip such that elements are in ascending order
    skip = skip[skiporder]                                                                                   ; rewrite skip so elements increase to make checks easier 

    ; if we don't have to skip any apertures, do cuts on whole structure  
    IF (skip[0] EQ -1) THEN BEGIN                                                                            ; if no apertures should be skipped over when making the cut, 
      FOR ap=0, naps-1 DO BEGIN                                                                                ; for each aperture in the structure: 
        badz = WHERE( (apstruct(ap).galz GT (z+zerr)) OR (apstruct(ap).galz LT (z-zerr)) )                       ; find where galaxy redshifts don't fall near target redshift
        apstruct = remover(apstruct, ap, badz)                                                                   ; remove those galaxies from the structure 
      ENDFOR                                                                                                   ; all the apertures in the structure have been filtered by z 

    ; if we do have to skip any apertures, check each aperture to see whether or not we're skipping it
    ENDIF ELSE IF (skip[0] GE 0) THEN BEGIN                                                                  ; if at least one aperture should be skipped when making the cut, 
      FOR ap=0, naps-1 DO BEGIN                                                                                ; for each aperture in the structure:   
        skipcheck = WHERE(skip EQ ap, skipthis)                                                                  ; see if this aperture is supposed to be skipped 
        IF (skipthis EQ 0) THEN BEGIN                                                                            ; if this aperture is NOT supposed to be skipped, 
          badz = WHERE( (apstruct(ap).galz GT (z+zerr)) OR (apstruct(ap).galz LT (z-zerr)) )                       ; find where galaxy redshifts don't fall near target redshift
          apstruct = remover(apstruct, ap, badz)                                                                   ; remove those galaxies from the structure 
        ENDIF                                                                                                    ; skipthis = 0 (ie, this aperture wasn't in the skip array)
      ENDFOR                                                                                                   ; all desired apertures have been filtered by z, skips were unaltered

    ; last, make sure that skip wasn't entered wrong and if it was, make note of it 
    ENDIF ELSE BEGIN                                                                                         ; if the first element of skip is less than -1, 
      print, 'Error: Invalid aperture index in "skip" array. Cuts will not be made.'                           ; print an error message and don't alter the input structure
    ENDELSE                                                                                                  ; all possible values of first entry in "skip" have been accounted for 

    ; write the resulting structure to a FITS file 
    inputzstring = strcompress(string(z, format='(D0.3)'), /remove_all)                                                          ; convert input z into a string 
    zstring = strjoin(strsplit(inputzstring, '.', /extract), 'p')                                                                ; replace dot in z with "p" for use in filename
    newfile = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/z/'+zstring+'/'+filename+'_z' + zstring + '.fits'  ; make a new filename for the filtered structure
    mwrfits, apstruct, newfile, /create                                                                                          ; write the filtered structure to a fits file 
  END       ; of zfilter procedure









  PRO colorline, filename, m, b, c, useblue, nameadd, skip
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; colorline procedure
  ;  Reads in a FITS file created by the blobsampling or fieldsampling procedures and filters that file by making simple straight-line color cuts to galaxy members. Galaxies 
  ;    that are kept are the ones whose colors are consistent with some target redshift; the straight line whose parameters are passed into this procedure marks the boundary
  ;    between galaxies with colors consistent with the target redshift and those whose colors are inconsistent with the target redshift. 
  ; 
  ; INPUT: filename - a string giving the name of a FITS file containing galaxy IDs within multiple apertures 
  ;        m - the slope of the color cut boundary line
  ;        b - the y-intercept of the color cut boundary line
  ;        c - a cap value above which all galaxies have colors consistent with the target redshift and are therefore all kept 
  ;        useblue - a string specifying which blue band to use for colors if there are multiple blue bands (should be 'blue' or 'altblue') 
  ;        nameadd - a string to be appended onto the original filename to create the name of this procedure's output file 
  ;        skip - an integer or array of integers giving the index/indices of a catalog to be skipped when applying cuts; to skip no indices, set skip to -1 
  ;
  ; OUTPUT: Filters galaxies with colors below and to the right of the given line out of the input structure and writes the filtered structure into a new FITS file.  
  ;
  ; Uses the remover function. 
  ;
  ; NOTE: This procedure is meant to be run by the makeacut procedure rather than being called on its own. 
  ;
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
    FORWARD_FUNCTION remover                                                                                 ; notify IDL that this function will be used 
    apstruct = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+filename+'.fits',1) ; read the apertures file into a structure
    naps = n_elements(apstruct)                                                                              ; count the number of apertures in the structure
    skiporder = sort(skip)                                                                                   ; get indices of skip such that elements are in ascending order
    skip = skip[skiporder]                                                                                   ; rewrite skip so elements increase to make checks easier 

    IF ( (useblue NE 'blue') AND (useblue NE 'altblue') ) THEN BEGIN                    ; make sure user correctly specified which blue filter to use 
      print, 'Specified blue filter does not exist. No cuts will be made.'                ; if not, tell user so and make no cuts 
    ENDIF ELSE BEGIN                                                                    ; if blue filter was specified correctly, proceed with cuts 

      ; do cuts on whole structure if we don't have to skip any apertures 
      IF (skip[0] EQ -1) THEN BEGIN                                                       ; if no apertures should be skipped over when making the cut, 
        FOR ap=0, naps-1 DO BEGIN                                                           ; for each aperture in the structure: 
          IF (useblue EQ 'blue') THEN blue = apstruct(ap).blue                                ; use the "blue" tag if user specified "blue"
          IF (useblue EQ 'altblue') THEN blue = apstruct(ap).altblue                          ; use the "altblue" tag if user specified "altblue"
          xcolor = blue - apstruct(ap).mid                                                    ; define the x axis of a color diagram (blue - middle filter)
          ycolor = apstruct(ap).mid - apstruct(ap).galmags                                    ; define the y axis of a color diagram (middle - red filter)      
          bad = WHERE ( (ycolor LT (xcolor*m + b)) AND (ycolor GE c), ngals )                 ; get indices of galaxies that fall below the straight-line cut 
          IF (ngals NE 0) THEN BEGIN                                                          ; if there are any galaxies that don't make the cut,
            apstruct = remover(apstruct, ap, bad)                                               ; remove them from the structure 
          ENDIF                                                                               ; this aperture has now been filtered 
        ENDFOR                                                                              ; all apertures have now been filtered 

      ; if we do have to skip any apertures, check each aperture to see whether or not we're skipping it
      ENDIF ELSE IF (skip[0] GE 0) THEN BEGIN                                             ; if at least one aperture should be skipped when making the cut, 
        FOR ap=0, naps-1 DO BEGIN                                                           ; for each aperture in the structure: 
          skipcheck = WHERE(skip EQ ap, skipthis)                                             ; see if this aperture is supposed to be skipped 
          IF (skipthis EQ 0) THEN BEGIN                                                       ; if this aperture is NOT supposed to be skipped, 
            IF (useblue EQ 'blue') THEN blue = apstruct(ap).blue                                ; use the "blue" tag if user specified "blue"
            IF (useblue EQ 'altblue') THEN blue = apstruct(ap).altblue                          ; use the "altblue" tag if user specified "altblue"
            xcolor = blue - apstruct(ap).mid                                                    ; define the x axis of a color diagram (blue - middle filter)
            ycolor = apstruct(ap).mid - apstruct(ap).galmags                                    ; define the y axis of a color diagram (middle - red filter)      
            bad = WHERE ( (ycolor LT (xcolor*m + b)) AND (ycolor GE c), ngals )                 ; get indices of galaxies that fall below the straight-line cut 
            IF (ngals NE 0) THEN BEGIN                                                          ; if there are any galaxies that don't make the cut,
              apstruct = remover(apstruct, ap, bad)                                               ; remove them from the structure 
            ENDIF                                                                               ; this aperture has now been filtered 
          ENDIF                                                                               ; skipthis = 0 (ie, this aperture wasn't in the skip array)
        ENDFOR                                                                              ; all desired apertures have been filtered, skips were unaltered
  
      ; last, make sure that skip wasn't entered wrong and if it was, make note of it 
      ENDIF ELSE BEGIN                                                                    ; if the first element of skip is less than -1, 
        print, 'Error: Invalid aperture index in "skip" array. Cuts will not be made.'      ; print an error message and don't alter the input structure
      ENDELSE                                                                             ; all possible values of first entry in "skip" have been accounted for 

    ; write the resulting structure to a FITS file 
    ENDELSE                                                                             ; the "useblue" parameter was input correctly 
    newfile = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/colorline/'+nameadd+'/'+filename+'_'+nameadd+'.fits' ; make a new filename for the filtered structure
    mwrfits, apstruct, newfile, /create                                                                                            ; write the filtered structure to a fits file 
  END        ; of colorline procedure 









  PRO colorgrid, filename, skip      ; plus other things 
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; colorline procedure
  ;  Reads in a FITS file created by the blobsampling or fieldsampling procedures and filters that file by ?????????????????????????????????????????
  ; 
  ; INPUT: filename - a string giving the name of a FITS file containing galaxy IDs within multiple apertures 
  ;        ?????????????????????????????????????????????????????????????
  ;        skip - an integer or array of integers giving the index/indices of a catalog to be skipped when applying cuts; to skip no indices, set skip to -1 
  ;
  ; OUTPUT: ?????????????????  
  ;
  ; Uses the remover function. 
  ;
  ; NOTE: This procedure is meant to be run by the makeacut procedure rather than being called on its own. 
  ;
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
    FORWARD_FUNCTION remover                                                                                 ; notify IDL that this function will be used 
    skiporder = sort(skip)                                                                                   ; get indices of skip such that elements are in ascending order
    skip = skip[skiporder]                                                                                   ; rewrite skip so elements increase to avoid later loop issues
    print, 'Color gridding has not been implemented yet. Sorry for the inconvenience.' 
  END        ; of colorgrid procedure 









  FUNCTION remover, apstruct, ap, bad 
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; remover function
  ;  Removes user-specified galaxies from an aperture structure (ie, one created by blobsampling or fieldsampling) by replacing their information with default values. 
  ;
  ; INPUT: apstruct - a structure containing galaxy IDs within multiple apertures
  ;        ap - integer giving the aperture in apstruct in which to remove bad galaxies 
  ;        bad - an integer or array of integers corresponding to the index/indices of galaxies in apstruct(ap) to be removed 
  ;
  ; OUTPUT: filters (ie, removes bad galaxies from) apstruct, then returns the filtered structure 
  ; 
  ; Uses no other functions nor procedures. 
  ;
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
    apstruct(ap).galIDs[bad] = -1                                                                           ; remove bad galaxies' IDs from the structure
    apstruct(ap).galz[bad] = -1.                                                                            ; remove bad galaxies' redshifts from the structure
    apstruct(ap).galmags[bad] = 99.                                                                         ; remove bad galaxies' red magnitudes from the structure 
    apstruct(ap).mid[bad] = 99.                                                                             ; remove bad galaxies' middle magnitudes from the structure 
    apstruct(ap).blue[bad] = 99.                                                                            ; remove bad galaxies' blue magnitudes from the structure 
    apstruct(ap).altblue[bad] = 99.                                                                         ; remove bad galaxies' other blue magnitudes from the structure 
    RETURN, apstruct                                                                                        ; return the filtered structure 
  END      ; of remover function 







  FUNCTION get_galaxies, ra, dec, radius, data                                                                                                    
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
    distance = sqrt(  ( (ra - data.ra) * cos(dec*!dPI/180.) )^2. + (dec - data.dec)^2. )    ; Pythagorean theorem               
    galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)                         ; the sources have to be contained within the aperture 
    RETURN, galaxies_in_aperture      
  END    









  PRO apfinder, npix_x, npix_y, ap_radius_pix, ap_radius, dim, bpm, hdr, ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec, blobcoords=blobcoords       
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;                
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
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;           
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
         ; make sure the box dimensions match the dimensions of the aperture image (this gets weird because of fractions of pixels)       
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
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
 
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
       struct_replace_field, blobs, 'FIELD15', blobs.FIELD15, newtag='bpmtype'    ; indicator of whether bpm has bad pixels marked as 0s or 1s 
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
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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
  ;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; DENSITYGETTER.PRO 
;  Contains all modules used for getting the galaxy density within blobs and in fields. 
;  CONTENTS: 
;      density----------------------------procedure   (main module for getting blob densities and average field densities)
;      countgals--------------------------procedure   (main module for getting galaxy counts in apertures)
;      bin_galaxies-----------------------function   
;      Poisson_error----------------------function
;      readdata---------------------------procedure
;      struct_replace_field---------------procedure
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;




PRO density, blobFITSname, fieldFITSnames, ap_radius_arcsec, brightestmag, dimmestmag, binsize, outputname, band=band
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; density procedure
;     Calculates galaxy number density and errors within blobs and fields for a specified magnitude band and type of color cut. 
;
; INPUT: blobFITSname - a string giving the name of a FITS file containing galaxy IDs within all blobs 
;        fieldFITSnames - an array of strings giving the names of FITS files for each field to be analyzed; each FITS file contains the galaxy IDs within a number of random apertures
;        ap_radius_arcsec - the radius in arcseconds of the apertures used in the above FITS files 
;        brightestmag - the lower (brightest) bound on galaxy magnitudes to be binned 
;        dimmestmag - the upper (dimmer) bound on galaxy magnitudes to be binned 
;        binsize - the size of the magnitude bins in which density will be calculated
;        outputname - a string containing the desired name of the file that will be created by this procedure 
;        --
;        band - an optional string keyword specifying which magnitude band to use for binning galaxies; if this isn't set, the red band will be used 
;
; OUTPUT: DensityResults - a structure containing: 
;                            a string for the name of an object (a blob or field), 
;                            the lower bounds of the magnitude bins in which galaxies are sorted (an array with nbins elements, where nbins is the number of magnitude bins), 
;                            the density per mag bin for the object (an array with nbins elements), and 
;                            the error on that density (an array with nbins elements) 
;                          for all blobs and fields; this is saved as a FITS file named as specified by outputname 
;
; Uses the bin_galaxies function, the Poisson_error function, and the readdata procedure. 
;   The latter of these uses the struct_replace_field procedure. 
;
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; SETUP 
     ; notify IDL of functions used in this procedure 
        FORWARD_FUNCTION bin_galaxies 
        FORWARD_FUNCTION Poisson_error 
     ; append paths onto FITS filenames provided by user input  
        fullblobFITSname = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+blobFITSname+'.fits'
        fullfieldFITSnames = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+fieldFITSnames+'.fits'
     readdata, blobs, fields, nblobs, nfields, ntotal           ; read in all the data and get number of blobs and fields 
     ap_radius = ap_radius_arcsec/3600.                         ; put the aperture radius in decimal degrees
     ; mag binning stuff:
         nbins = (dimmestmag - brightestmag + binsize)/binsize     ; the number of magnitude bins as determined from our range and binsize 
         mags = brightestmag + binsize*findgen(nbins)              ; an array of every magnitude being used 
         magrange=[brightestmag,dimmestmag]                        ; a two-element array containing only the endpoints of our magnitude range 
     ; set up a structure to contain all the density info; this will be saved as a FITS file to be plotted separately 
        ; note: mags is included to make plotting easy 
        singleobjectdensity = {name:'', mags:mags, densities:dblarr(nbins), derrs:dblarr(nbins)}   ; a structure with density information for a single blob or field 
        IF (n_elements(band) NE 0) THEN BEGIN   ; check to see if a band was given by the user 
          IF (band EQ 'blue') THEN BEGIN          ; see if the blue band is specified by user
            ; in this case, we need two mini-structures for every field because they each have two blue bands 
            DensityResults = replicate(singleobjectdensity, ntotal+nfields)          ; main structure with density info for all blobs and fields
          ENDIF ELSE IF ( (band NE 'red') AND (band NE 'mid') AND (band NE 'blue') ) THEN BEGIN    ; if a band was specified but doesn't exist, 
            print, 'Specified band does not exist. Using the red band.'                              ; just use the red band and let user know
            band = 'red'
            DensityResults = replicate(singleobjectdensity, ntotal)                                                          ; main structure
          ENDIF ELSE IF ( (band EQ 'red') OR (band EQ 'mid') ) THEN BEGIN                          ; if red or middle bands were specified 
            DensityResults = replicate(singleobjectdensity, ntotal)                                                          ; main structure
          ENDIF    
        ENDIF ELSE BEGIN                                                                           ; if a band wasn't specified 
          print, 'No band specified. Using the red band.'                                            ; just use the red band and let user know
          band = 'red' 
          DensityResults = replicate(singleobjectdensity, ntotal)                                                            ; main structure 
        ENDELSE   ; if a band wasn't given  


  ; BLOB DENSITIES
     blobapertures = mrdfits(fullblobFITSname,1)    ; read in the FITS file containing just the IDs of galaxies within the blobs 

     FOR blob=0, nblobs-1 DO BEGIN              ; loop through all the blobs 

       DensityResults(blob).name = blobs.blobname[blob]  ; first, fill what we know about DensityResults, which is this blob's name 
       blobgalIDs = blobapertures(blob).galIDs           ; get the galaxy IDs from the FITS file into an array - this also contains a bunch of -1s
       good = WHERE(blobgalIDs NE -1, napgals)             ; get the array indices where there are actual galaxy IDs instead of -1s
       IF (napgals LE 0) THEN BEGIN                      ; if there are no galaxies in the aperture, then all bins' densities and their errors are 0
         DensityResults(blob).densities[*] = 0.0
         DensityResults(blob).derrs[*] = 0.0 

       ENDIF ELSE BEGIN                         ; if there are some galaxies in this aperture 

         blobgalIDs = blobgalIDs[good]                     ; trim the -1s out of the array so we're only left with galaxy IDs
         nblobgals = n_elements(blobgalIDs)                ; How many galaxies are in the blob? 

         ; get the relevant information for all the galaxies in/near the blob
            blobcat = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+blobs.blobinfo[blob], 1)   ; read in the blob galaxy info catalog
            ; make a truncated version of blobcat that only contains info about the galaxies in/near the blob 
            truncated_blobcat = {ID:blobgalIDs, ra:fltarr(nblobgals), dec:fltarr(nblobgals), redmag:fltarr(nblobgals), midmag:fltarr(nblobgals), bluemag:fltarr(nblobgals)}
            FOR galaxy=0, nblobgals-1 DO BEGIN      ; fill in the info for each galaxy in truncated_blobcat 
              galindex = WHERE(blobcat.ID EQ truncated_blobcat.ID[galaxy])   ; find the index in blobcat where the ID matches this galaxy's ID 
              ; fill in the other fields of truncated_blobcat with the info for this galaxy 
                 truncated_blobcat.ra[galaxy] = blobcat.ra[galindex]
                 truncated_blobcat.dec[galaxy] = blobcat.dec[galindex]
                 truncated_blobcat.redmag[galaxy] = blobcat.redmag[galindex]
                 truncated_blobcat.midmag[galaxy] = blobcat.midmag[galindex]
                 truncated_blobcat.bluemag[galaxy] = blobcat.bluemag[galindex] 
            ENDFOR     ; all the galaxies in/near the blob (ie, in truncated_blobcat) 
       
         ; now we can bin up those galaxies and calculate density and errors in each bin to put in DensityResults 
            FOR mag=brightestmag, dimmestmag, binsize DO BEGIN        ; loop through the magnitude bins 
              ; before binning galaxies, determine which magnitude band to use for binning 
              IF (band EQ 'red') THEN BEGIN               ; if user specified red band (or something nonexistent) 
                datamags = truncated_blobcat.redmag         ; use the red band for magnitudes 
              ENDIF ELSE IF (band EQ 'mid') THEN BEGIN    ; if user specified mid band 
                datamags = truncated_blobcat.midmag         ; use the middle band for magnitudes 
              ENDIF ELSE IF (band EQ 'blue') THEN BEGIN   ; if user specified blue band 
                datamags = truncated_blobcat.bluemag        ; use the blue band for magnitudes 
              ENDIF      ; we've picked which magnitude band to use for binning galaxies 
              ; now bin galaxies and get densities 
              ngalinbin = bin_galaxies(mag, mag+binsize, datamags)                                                            ; number of galaxies in this bin 
              ngalerrbin = Poisson_error(ngalinbin)                                                                           ; Poisson error on above number 
              thisindex = WHERE(DensityResults(blob).mags EQ mag)                                                             ; the array index for this magnitude bin
              DensityResults(blob).densities[thisindex] = double(ngalinbin)/(!dPI*ap_radius^2.*blobapertures(blob).weight)    ; number density of galaxies in this bin 
              DensityResults(blob).derrs[thisindex] = double(ngalerrbin)/(!dPI*ap_radius^2.*blobapertures(blob).weight)       ; error on above number density 
            ENDFOR   ; all magnitude bins 

       ENDELSE       ; this blob has galaxies in it

     ; now that we have density results for all mag bins for this blob, we can move on to the next blob 
     ENDFOR     ; all the blobs

STOP

  ; FIELD DENSITIES
     FOR field=0, nfields-1 DO BEGIN               ; since each field has its own FITS file, start off with a loop 
       DensityResults(nblobs+field).name = fields.fieldname[field]  ; put the name of this field into DensityResults 
       fieldapertures = mrdfits(fullfieldFITSnames[field], 1)           ; read in the FITS file containing just the IDs of galaxies within each random aperture 
       n_apertures = n_elements(fieldapertures.aperture)            ; get the number of random apertures used to make the FITS file 
       fieldcat = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+fields.fieldinfo[field], 1)  ; read in field galaxy info catalog

       ; to get densities, we need an average over all random apertures for each bin, so go aperture by aperture first to average at the end 
       singleapdensities = {densities:fltarr(nbins), derrs:fltarr(nbins)}       ; make a structure with densities and errors in all mag bins for one aperture 
       ap_densities = replicate(singleapdensities, n_apertures)                 ; replicate the previous structure to include all apertures 
       totalapdensities = dblarr(nbins)                                         ; have an array that just includes all the apertures' densities added up, per mag bin 
       totalaperr = dblarr(nbins)                                               ; same as above but for errors on density 
       IF (band EQ 'blue') THEN BEGIN           ; repeat the above for the second blue band if blue band is specified 
         altap_densities = replicate(singleapdensities, n_apertures)            ; replicate the previous structure to include all apertures 
         alttotalapdensities = dblarr(nbins)                                    ; have an array that just includes all the apertures' densities added up, per mag bin 
         alttotalaperr = dblarr(nbins)                                          ; same as above but for errors on density 
         DensityResults(ntotal+field).name = fields.fieldname[field] + 'alt'    ; put the name of this field into DensityResults for alt band densities 
       ENDIF    ; we're using the blue band(s)  

       ; now place random apertures 
       FOR ap=0, n_apertures-1 DO BEGIN       ; loop over all apertures 

         galIDs = fieldapertures(ap).galIDs           ; get this aperture's galaxy IDs from the FITS file into an array - this also contains a bunch of -1s
         good = WHERE(galIDs NE -1, napgals)          ; get the array indices where there are actual galaxy IDs instead of -1s and the number of galaxies in the aperture
         IF (napgals LE 0) THEN BEGIN                 ; if there are no galaxies in the aperture, then all bins' densities and their errors are 0
           ap_densities(ap).densities[*] = 0.0
           ap_densities(ap).derrs[*] = 0.0 

         ENDIF ELSE BEGIN                         ; if there are some galaxies in this aperture 

           galIDs = galIDs[good]                        ; trim the -1s out of the array so we're only left with galaxy IDs 

           ; make a truncated version of fieldcat that only contains info about the galaxies in this particular aperture  
           apcat = {ID:galIDs, ra:fltarr(napgals), dec:fltarr(napgals), redmag:fltarr(napgals), midmag:fltarr(napgals), $
                    bluemag:fltarr(napgals), altbluemag:fltarr(napgals), z:fltarr(napgals)}
           FOR galaxy=0, napgals-1 DO BEGIN      ; fill in the info for each galaxy in apcat 
             galindex = WHERE(fieldcat.ID EQ apcat.ID[galaxy])   ; find the index in fieldcat where the ID matches this galaxy's ID 
             ; fill in the other fields of apcat with the info for this galaxy 
                apcat.ra[galaxy] = fieldcat.ra[galindex]
                apcat.dec[galaxy] = fieldcat.dec[galindex]
                apcat.redmag[galaxy] = fieldcat.redmag[galindex]
                apcat.midmag[galaxy] = fieldcat.midmag[galindex]
                apcat.bluemag[galaxy] = fieldcat.bluemag[galindex] 
                apcat.altbluemag[galaxy] = fieldcat.altbluemag[galindex] 
                apcat.z[galaxy] = fieldcat.z[galindex] 
           ENDFOR     ; all the galaxies in this aperture 

           ; now we can bin up those galaxies and calculate density and errors in each bin to put in ap_densities 
           FOR mag=brightestmag, dimmestmag, binsize DO BEGIN        ; loop through the magnitude bins 
             ; before binning galaxies, determine which magnitude band to use for binning 
             IF (band EQ 'red') THEN BEGIN               ; if user specified red band (or something nonexistent) 
               datamags = apcat.redmag                     ; use the red band for magnitudes 
             ENDIF ELSE IF (band EQ 'mid') THEN BEGIN    ; if user specified mid band 
               datamags = apcat.midmag                     ; use the middle band for magnitudes 
             ENDIF ELSE IF (band EQ 'blue') THEN BEGIN   ; if user specified blue band 
               datamags = apcat.bluemag                    ; use the blue band for magnitudes  
               altdatamags = apcat.altbluemag              ; use the alt blue band separately 
             ENDIF      ; we've picked which magnitude band to use for binning galaxies 
            ; now bin galaxies and get densities 
             ngalinbin = bin_galaxies(mag, mag+binsize, datamags)                                                     ; number of galaxies in this bin 
             ngalerrbin = Poisson_error(ngalinbin)                                                                    ; Poisson error on above number 
             thisindex = WHERE(DensityResults(nblobs+field).mags EQ mag)                                              ; the array index for this magnitude bin
             ap_densities(ap).densities[thisindex] = double(ngalinbin)/(!dPI*ap_radius^2.*fieldapertures(ap).weight)       ; number density of galaxies in this bin 
             ap_densities(ap).derrs[thisindex] = double(ngalerrbin)/(!dPI*ap_radius^2.*fieldapertures(ap).weight)          ; error on above number density 
             totalapdensities[thisindex] += ap_densities(ap).densities[thisindex]                                          ; add this aperture's bin density onto the total 
             totalaperr[thisindex] += ap_densities(ap).derrs[thisindex]                                                    ; add this aperture's bin error onto the total 
             ; if blue band was specified, we have to do this twice since there are two blue bands in the 3D-HST fields
             IF (band EQ 'blue') THEN BEGIN        ; have to do the alt blue band too
               altngalinbin = bin_galaxies(mag, mag+binsize, altdatamags)                                                      ; number of galaxies in this bin 
               altngalerrbin = Poisson_error(ngalinbin)                                                                        ; Poisson error on above number 
               altap_densities(ap).densities[thisindex] = double(altngalinbin)/(!dPI*ap_radius^2.*fieldapertures(ap).weight)   ; number density of galaxies in this bin 
               altap_densities(ap).derrs[thisindex] = double(altngalerrbin)/(!dPI*ap_radius^2.*fieldapertures(ap).weight)      ; error on above number density 
               alttotalapdensities[thisindex] += altap_densities(ap).densities[thisindex]                                      ; add this aperture's bin density onto the total 
               alttotalaperr[thisindex] += altap_densities(ap).derrs[thisindex]                                                ; add this aperture's bin error onto the total 
             ENDIF   ; we're using the blue band for binning galaxies 
           ENDFOR   ; all magnitude bins 

         ENDELSE    ; if there were any galaxies in this aperture  

       ENDFOR       ; all apertures in this field 

       ; now turn the totals into averages to fill in DensityResults 
       DensityResults(nblobs+field).densities = totalapdensities / double(n_apertures)       ; the average field density per mag bin 
       DensityResults(nblobs+field).derrs = totalaperr / double(n_apertures)                 ; the average error on the above 
       IF (band EQ 'blue') THEN BEGIN     ; add alt blue mag densities as well 
         DensityResults(ntotal+field).densities = alttotalapdensities / double(n_apertures)       ; the average field density per mag bin 
         DensityResults(ntotal+field).derrs = alttotalaperr / double(n_apertures)                 ; the average error on the above 
       ENDIF    ; we're binning galaxies by blue magnitude 

     ENDFOR         ; all the fields 

  ; Finally, output the results in a FITS file to be plotted later 
  bandname = strcompress(string(band))                                                       ; we're going to put the name of the band in the filename 
  outputpath = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/densities/'+outputname+'_'+bandname+'.fits'     ; full path to where the output file will be saved 
  mwrfits, DensityResults, outputpath, /create                                               ; write DensityResults to the output FITS file 

  ; just a test to make sure everything matches up
  help, DensityResults             ; give info on the data written to the FITS file 
  test = mrdfits(outputpath, 1)    ; read in the newly-created FITS file as a test
  help, test                       ; give info on the read-in FITS file to make sure it matches the original data 

END      ; end of density procedure   









PRO countgals, FITSname, brightestmag, dimmestmag, binsize, outputname
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; countgals procedure 
;  Goes through a set of apertures and counts the number of galaxies per magnitude bin in each aperture, then exports this information as a FITS file. 
;
; INPUT: FITSname - a string giving the name of a FITS file containing galaxy IDs within multiple apertures (blobs or random apertures in a field) 
;        brightestmag - the lower (brightest) bound on galaxy magnitudes to be binned 
;        dimmestmag - the upper (dimmer) bound on galaxy magnitudes to be binned 
;        binsize - the size of the magnitude bins in which density will be calculated
;        outputname - a string containing the desired name of the FITS file that will be created by this procedure 
;
; OUTPUT: countstruct - a structure containing each aperture's name and the number of galaxies it has in each magnitude bin 
;
; Uses the bin_galaxies function. 
;
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; read in data and define parameters based on input 
  fullFITSname = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/apertures/'+FITSname+'.fits'  ; tack the path and extension onto the input filename 
  apertures = mrdfits(fullFITSname,1)                                                                ; read in the input FITS file containing the IDs of galaxies within apertures
  nApertures = n_elements(apertures)                                                                 ; number of apertures in the input FITS file 
  nbins = (dimmestmag - brightestmag + binsize)/binsize                                              ; number of magnitude bins as determined from our range and binsize 
  mags = brightestmag + binsize*findgen(nbins)                                                       ; an array of every magnitude being used 
  validblob = tag_exist(apertures, 'blobname')                                                       ; this will be 1b if the first tag is 'blobname' and will be 0b otherwise
  validfield = tag_exist(apertures, 'aperture')                                                      ; this will be 1b if the first tag is 'aperture' and will be 0b otherwise
  IF (validblob EQ 1b) THEN BEGIN                                                                    ; if the first tag of input structure is 'blobname'
    names = apertures(*).blobname                                                                      ; then use 'blobname' input tag as the 'name' output tag for this aperture 
  ENDIF ELSE IF (validfield EQ 1b) THEN BEGIN                                                        ; if the first tag of input structure is 'aperture'
    names = apertures(*).aperture                                                                      ; then use 'aperture' input tag as the 'name' output tag for this aperture 
  ENDIF                                                                                              ; now we have an array containing each aperture's name 
  names = strcompress(string(names), /remove)                                                        ; make sure everything in the names array is a string 

  ; make a structure to export as a FITS file 
  dummystruct = {name:'name', mags:mags, weight:1., galcounts:lonarr(nbins)}  ; dummy structure for one aperture - includes its name, mag bins, weight, and # galaxies in it per mag bin 
  countstruct = replicate(dummystruct, nApertures)                            ; replicate dummy structure so there's one structure per aperture to comprise the main output structure 

  ; fill in countstruct one aperture at a time 
  TIC                                                                                ; time this to ensure it's efficient 
  FOR ap=0, nApertures-1 DO BEGIN                                                    ; loop through all apertures 
    ; first, fill in the 'name' and 'weight' tags of countstruct for this aperture
    countstruct(ap).name = names[ap]                                                    ; get name of this aperture from the names array 
    countstruct(ap).weight = apertures(ap).weight                                       ; get weight of aperture (ie, fraction of good pixels) from input structure 
    ; now fill in the 'galcounts' tag of countstruct for this aperture
    good = WHERE(apertures(ap).galIDs NE -1, napgals)                                   ; get the array indices where there are actual galaxy IDs instead of -1s
    IF (napgals LE 0) THEN BEGIN                                                        ; if there are no galaxies in the aperture, 
      countstruct(ap).galcounts[*] = 0                                                    ; then all bins contain 0 galaxies 
    ENDIF ELSE BEGIN                                                                    ; if this aperture DOES contain at least one galaxy, 
      ; get number of galaxies in each bin to fill in 'galcounts' array 
      FOR mag=brightestmag, dimmestmag, binsize DO BEGIN                                       ; loop through the magnitude bins 
        index = WHERE(mags EQ mag)                                                               ; the array index for this magnitude bin
        countstruct(ap).galcounts[index] = bin_galaxies(mag, mag+binsize, apertures(ap).galmags)   ; put number of galaxies in this bin into corresponding index of galcounts array
      ENDFOR                                                                                   ; all magnitude bins 
    ENDELSE                                                                             ; if this aperture contained at least one galaxy
  ENDFOR                                                                             ; every aperture in the input structure
  TOC                                                                                ; see how long it took to run through this loop 
 
  ; output the results in a FITS file to be plotted later 
  outputpath = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/galcounts/'+outputname+'.fits'     ; full path to where the output file will be saved 
  mwrfits, countstruct, outputpath, /create                                                             ; export countstruct as a FITS file with the above name 
   
END   ; of countgals procedure









FUNCTION bin_galaxies, brightmag, dimmag, datamags 
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------; 
; bin_galaxies function     
;    determines the number of galaxies in a catalog that have magnitudes between two specified magnitudes (ie, galaxies in a specified bin) 
;                                                                                                                
; INPUTS: brightmag - the lower bound (brightest magnitude) of a magnitude range (bin)     
;         dimmag - the upper bound (faintest magnitude) of a magnitude range (bin)  
;         datamags - the magnitudes of galaxies in a catalog            
;                          
; OUTPUT: n_galaxies - the number of galaxies in the catalog which have magnitudes between brightmag and dimmag 
;
; Uses no other functions nor procedures. 
;                                                                                         
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  galaxies_in_bin = WHERE (((datamags GE brightmag) AND (datamags LT dimmag)), n_galaxies)    ; index numbers of catalog galaxies that fall into the specified mag bin 
  RETURN, n_galaxies                                                                          ; return the number of galaxies in the bin 
END  









FUNCTION Poisson_error, ngal                                                                                                                                          
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; Poisson_error function    
;    Calculates the Poisson error on a number of galaxies by taking the square root of that number: S = sqrt(N).   
;                                                                                                                    
; INPUT: ngal - a number (of galaxies) for which we want to determine uncertainty                 
;                                               
; OUTPUT: error - the Poisson error on the input number     
; 
; Uses no other functions nor procedures. 
;                                 
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  error = sqrt(double(ngal)) 
  RETURN, error  
END           









PRO readdata, blobs, fields, nblobs, nfields, ntotal  
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
 
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
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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

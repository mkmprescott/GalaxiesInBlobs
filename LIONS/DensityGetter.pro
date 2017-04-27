;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; DENSITYGETTER.PRO 
;  Contains all modules used for getting the galaxy density within blobs and in fields. 
;  CONTENTS: 
;      density----------------------------procedure   (main module) (WILL BE UPDATED AFTER MAKEACUT IS FINISHED) 
;      bin_galaxies=----------------------function   
;      Poisson_error----------------------function
;      makeacut---------------------------procedure  (IN PROGRESS) 
;      blobcatreader----------------------procedure
;      fieldcatreader---------------------procedure
;      readdata---------------------------procedure
;      struct_replace_field---------------procedure
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;




PRO density, cuts=cuttype, band=band 
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
;        cuts - an optional keyword which, if set, specifies a type of cut to make to the blob and field data 
;        band - an option keyword specifying which magnitude band to use for binning galaxies; if this isn't set, the red band will be used 
;
; OUTPUT: DensityResults - a structure containing: 
;                            a string for the name of an object (a blob or field), 
;                            the lower bounds of the magnitude bins in which galaxies are sorted (an array with nbins elements, where nbins is the number of magnitude bins), 
;                            the density per mag bin for the object (an array with nbins elements), and 
;                            the error on that density (an array with nbins elements) 
;                          for all blobs and fields; this is saved as a FITS file named as specified by outputname 
;
; Uses the bin_galaxies function, the Poisson_error function, the makeacut procedure, the blobcatreader procedure, the fieldcatreader procedure, and the readdata procedure. 
;   The latter three of these use the struct_replace_field procedure. 
;
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; SETUP 
     FORWARD_FUNCTION bin_galaxies 
     FORWARD_FUNCTION Poisson_error 

     ; USER-DEFINED INPUT: 
     blobFITSname = 'LABapertures_10arcsec.fits'                 ; name of blobs FITS file 
     fieldFITSnames = ['HSTapertures_10arcsec_10000aps.fits']    ; array of fields' FITS files (because we can have multiple fields each with its own random apertures)
     ap_radius_arcsec = 10.                                      ; radius of apertures used to count galaxies in blobs and fields, in arcseconds 
     brightestmag = 22.                                          ; smallest (brightest) magnitude of galaxies we can reasonably see in our images 
     dimmestmag = 30.                                            ; biggest (faintest) magnitude of galaxies we can reasonably see in our images 
     binsize = 0.5                                               ; size of magnitude bins
     outputname = 'DENSITIES_10arcsec_10000aps.fits'             ; name of FITS file in which results of this analysis will be stored 

     ; append paths onto FITS filenames provided above  
        blobFITSname = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/FITS/'+blobFITSname
        fieldFITSnames = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/FITS/'+fieldFITSnames 
     readdata, blobs, fields, nblobs, nfields, ntotal           ; read in all the data and get number of blobs and fields 
     ap_radius = ap_radius_arcsec/3600.                         ; put the aperture radius in decimal degrees
     ; mag binning stuff:
         nbins = (dimmestmag - brightestmag + binsize)/binsize     ; the number of magnitude bins as determined from our range and binsize 
         mags = brightestmag + binsize*findgen(nbins)              ; an array of every magnitude being used 
         magrange=[brightestmag,dimmestmag]                        ; a two-element array containing only the endpoints of our magnitude range 
     ; set up a structure to contain all the density info; this will be saved as a FITS file to be plotted separately 
        ; note: mags is included to make plotting easy 
        singleobjectdensity = {name:'', mags:mags, densities:dblarr(nbins), derrs:dblarr(nbins)}   ; a structure with density information for a single blob or field 
        IF (n_elements(band) NE 0) THEN BEGIN   ; check to see if the blue band is specified by user (start with checking if a band was given at all) 
          IF (band EQ 'blue') THEN BEGIN
            ; in this case, we need two mini-structures for every field because they each have two blue bands 
            DensityResults = replicate(singleobjectdensity, ntotal+nfields)          ; main structure with density info for all blobs and fields
          ENDIF      ; blue band was specified 
        ENDIF ELSE BEGIN 
          ; if the band wasn't specified or specified band isn't a real band, use the red band for density analysis and let the user know 
          IF (n_elements(band) EQ 0) THEN BEGIN   ; if no band was specified 
            print, 'No band specified. Using the red band.' 
            band = 'red' 
          ENDIF ELSE IF ( (band NE 'red') AND (band NE 'mid') AND (band NE 'blue') ) THEN BEGIN    ; if the specified band doesn't exist 
            print, 'Specified band does not exist. Using the red band.' 
            band = 'red'
          ENDIF     
          DensityResults = replicate(singleobjectdensity, ntotal)                ; main structure with density info for all blobs and fields
        ENDELSE   ; if blue band isn't specified 


  ; BLOB DENSITIES
     blobapertures = mrdfits(blobFITSname,1)    ; read in the FITS file containing just the IDs of galaxies within the blobs 

     FOR blob=0, nblobs-1 DO BEGIN              ; loop through all the blobs 

       DensityResults(blob).name = blobs.blobname[blob]  ; first, fill what we know about DensityResults, which is this blob's name 
       blobgalIDs = blobapertures(blob).galIDs           ; get the galaxy IDs from the FITS file into an array - this also contains a bunch of -1s
       good = WHERE(blobgalIDs NE -1)                    ; get the array indices where there are actual galaxy IDs instead of -1s
       blobgalIDs = blobgalIDs[good]                     ; trim the -1s out of the array so we're only left with galaxy IDs
       nblobgals = n_elements(blobgalIDs)                ; How many galaxies are in the blob? 

       ; get the relevant information for all the galaxies in/near the blob 
          blobcatreader, blobs.blobinfo[blob], blobcat   ; read in the catalog containing ID, RA, dec, and magnitudes for all galaxies in the blob image 
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

          ; apply cuts to the truncated blobcat if desired 
          IF (n_elements(cuts) NE 0) THEN BEGIN      ; this is where any user-specified cuts would be implemented 
            print, 'The code to apply cuts has not been written yet, so no cuts will be implemented. Sorry.'
            print, 'In the future, this part of the code will alter the truncated_blobcat structure by removing some galaxies.'
          ENDIF  ; if the user wanted any membership cuts applied to the data 
       
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
            ngalinbin = bin_galaxies(mag, mag+binsize, datamags)                                 ; number of galaxies in this bin 
            ngalerrbin = Poisson_error(ngalinbin)                                                ; Poisson error on above number 
            thisindex = WHERE(DensityResults(blob).mags EQ mag)                                 ; the array index for this magnitude bin
            DensityResults(blob).densities[thisindex] = double(ngalinbin)/(!dPI*ap_radius^2.)    ; number density of galaxies in this bin 
            DensityResults(blob).derrs[thisindex] = double(ngalerrbin)/(!dPI*ap_radius^2.)       ; error on above number density 
          ENDFOR   ; all magnitude bins 

     ; now that we have density results for all mag bins for this blob, we can move on to the next blob 
     ENDFOR     ; all the blobs



  ; FIELD DENSITIES
     FOR field=0, nfields-1 DO BEGIN               ; since each field has its own FITS file, start off with a loop 
       DensityResults(nblobs+field).name = fields.fieldname[field]  ; put the name of this field into DensityResults 
       fieldapertures = mrdfits(fieldFITSnames[field], 1)           ; read in the FITS file containing just the IDs of galaxies within each random aperture 
       n_apertures = n_elements(fieldapertures.aperture)            ; get the number of random apertures used to make the FITS file 
       fieldcatreader, fields.fieldinfo[field], fieldcat            ; read in the catalog containing ID, RA, dec, and magnitudes for all galaxies in the entire field 

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
       FOR ap=0, n_apertures-1 DO BEGIN       ; loop over all apertures 
         galIDs = fieldapertures(ap).galIDs           ; get this aperture's galaxy IDs from the FITS file into an array - this also contains a bunch of -1s
         good = WHERE(galIDs NE -1, napgals)          ; get the array indices where there are actual galaxy IDs instead of -1s and the number of galaxies in the aperture
         IF (napgals LE 0) THEN BEGIN             ; if there are no galaxies in the aperture, then all bins' densities and their errors are 0
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

           ; apply cuts to the apcat if desired 
           IF (n_elements(cuts) NE 0) THEN BEGIN      ; this is where any user-specified cuts would be implemented 
             print, 'The code to apply cuts has not been written yet, so no cuts will be implemented. Sorry.'
             print, 'In the future, this part of the code will alter the apcat structure by removing some galaxies.'
           ENDIF  ; if the user wanted any membership cuts applied to the data 

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
  bandname = strcompress(string(band))                                           ; we're going to put the name of the band in the filename 
  outputpath = '/boomerang-data/alhall/GalaxiesInBlobs/'+bandname+outputname     ; full path to where the output file will be saved 
  mwrfits, DensityResults, outputpath, /create                                   ; write DensityResults to a FITS file 

help, DensityResults 
test = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/'+bandname+outputname, 1)
help, test 

END      ; end of density procedure   









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









PRO makeacut, catalog, newcatalog, cuttype
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; makeacut procedure
;    Makes a user-specified cut to the galaxy membership in a specified catalog. 
;
; INPUT: catalog - a catalog containing information about galaxies in some image of a blob or field 
;        newcatalog - the name of the post-cut catalog that will be returned
;        cuttype - the kind of cut to be made to catalog to turn it into newcatalog 
;
; OUTPUT: returns newcatalog with all the same fields as the original input catalog, but with some galaxies removed as per the specified cut
;
; Uses the following: 
;
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

END









PRO blobcatreader, blobcatname, blobcat
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; blobcatreader procedure
;    Reads a csv file containing information about galaxies in a blob and puts that information into a structure whose fields are appropriately named. 
;
; INPUT: blobcatname - string of the filename of a csv file containing information about galaxies in a blob 
;        blobcat - the name of the structure into which the catalog will be read
;
; OUTPUT: blobcat - a structure containing all the information in the csv file organized into appropriately-named fields 
;
; Uses the struct_replace_field procedure. 
;
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  catname = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+ blobcatname   ; the full path to this blob's catalog
  blobcat = read_csv(catname)                                                                 ; read in the blob catalog
  ; next, rename the fields of the catalog structure to be useful 
     struct_replace_field, blobcat, 'FIELD1', blobcat.FIELD1, newtag='ID'           ; galaxy ID
     struct_replace_field, blobcat, 'FIELD2', blobcat.FIELD2, newtag='ra'           ; galaxy RA
     struct_replace_field, blobcat, 'FIELD3', blobcat.FIELD3, newtag='dec'          ; galaxy declination
     struct_replace_field, blobcat, 'FIELD4', blobcat.FIELD4, newtag='redmag'       ; the reddest-band magnitude of each galaxy
     struct_replace_field, blobcat, 'FIELD5', blobcat.FIELD5, newtag='midmag'       ; the middle-band magnitude of each galaxy
     struct_replace_field, blobcat, 'FIELD6', blobcat.FIELD6, newtag='bluemag'      ; the bluest-band magnitude of each galaxy 
END            ; of blobcatreader procedure 









PRO fieldcatreader, fieldcatname, fieldcat
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; fieldcatreader procedure
;    Reads a csv file containing information about galaxies in a field and puts that information into a structure whose fields are appropriately named. 
;
; INPUT: fieldcatname - string of the filename of a csv file containing information about galaxies in a field 
;        fieldcat - the name of the structure into which the catalog will be read
;
; OUTPUT: fieldcat - a structure containing all the information in the csv file organized into relevantly-named fields 
;
; Uses the struct_replace_field procedure. 
;
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  catname = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+ fieldcatname   ; the full path to this field's catalog
  fieldcat = read_csv(catname)                                                                 ; read in the field catalog
  ; next, rename the fields of the catalog structure to be useful 
     struct_replace_field, fieldcat, 'FIELD1', fieldcat.FIELD1, newtag='ID'           ; galaxy ID
     struct_replace_field, fieldcat, 'FIELD2', fieldcat.FIELD2, newtag='ra'           ; galaxy RA
     struct_replace_field, fieldcat, 'FIELD3', fieldcat.FIELD3, newtag='dec'          ; galaxy declination
     struct_replace_field, fieldcat, 'FIELD4', fieldcat.FIELD4, newtag='redmag'       ; the reddest-band (F140W) magnitude of each galaxy
     struct_replace_field, fieldcat, 'FIELD5', fieldcat.FIELD5, newtag='midmag'       ; the middle-band (F814W) magnitude of each galaxy
     struct_replace_field, fieldcat, 'FIELD6', fieldcat.FIELD6, newtag='bluemag'      ; the primary bluest-band (F606W) magnitude of each galaxy   
     struct_replace_field, fieldcat, 'FIELD7', fieldcat.FIELD7, newtag='altbluemag'   ; the secondary bluest-band (F475W) magnitude of each galaxy 
     struct_replace_field, fieldcat, 'FIELD8', fieldcat.FIELD8, newtag='z'            ; galaxy redshift 
END            ; of fieldcatreader procedure









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
     struct_replace_field, blobs, 'FIELD04', blobs.FIELD04, newtag='bluecat'    ; filename of sources catalog with blue-filter mags* 
     struct_replace_field, blobs, 'FIELD05', blobs.FIELD05, newtag='midcat'     ; filename of sources catalog with mid-filter mags* 
     struct_replace_field, blobs, 'FIELD06', blobs.FIELD06, newtag='redcat'     ; filename of sources catalog with red-filter mags*
     struct_replace_field, blobs, 'FIELD07', blobs.FIELD07, newtag='stack'      ; filename of stacked image for each blob**  
     struct_replace_field, blobs, 'FIELD08', blobs.FIELD08, newtag='mask'       ; filename of mask image for each blob** 
     struct_replace_field, blobs, 'FIELD09', blobs.FIELD09, newtag='z'          ; redshift 
     struct_replace_field, blobs, 'FIELD10', blobs.FIELD10, newtag='m'          ; slope of color-cut line 
     struct_replace_field, blobs, 'FIELD11', blobs.FIELD11, newtag='b'          ; y-intercept of color-cut line 
     struct_replace_field, blobs, 'FIELD12', blobs.FIELD12, newtag='cap'        ; mag above which all colors are consistent with redshift 
     struct_replace_field, blobs, 'FIELD13', blobs.FIELD13, newtag='bluename'   ; name of the bluest filter used to observe the blob region 
     struct_replace_field, blobs, 'FIELD14', blobs.FIELD14, newtag='midbname'   ; name of the middle filter used to observe the blob region 
     struct_replace_field, blobs, 'FIELD15', blobs.FIELD15, newtag='redname'    ; name of the reddest filter used to observe the blob region 
     struct_replace_field, blobs, 'FIELD16', blobs.FIELD16, newtag='blobinfo'   ; name of the catalog containing all the magnitude, RA, DEC, etc info of sources 
     struct_replace_field, blobs, 'FIELD17', blobs.FIELD17, newtag='pixelscale' ; pixel scale in arceseconds per pixel ("/px) for each blob image 
     nblobs = n_elements(blobs.blobname)                                        ; number of blobs in the blobs.csv file 
     ; * note: only used by BlobcatMaker.pro; to use, must be preceded by '/boomerang-data/alhall/GalaxiesInBlobs/Data/BlobData/BlobCatalogs/'+ 
     ; ** note: to use, must be preceded by '/boomerang-data/alhall/GalaxiesInBlobs/Data/BlobData/BlobStacks/'+

  ; read in the file containing all field information
     fields = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/fields.csv', n_table_header=3)
  ; assign names to each field of the "fields" structure
     struct_replace_field, fields, 'FIELD1', fields.FIELD1, newtag='fieldname'  ; name of each field (for example, GOODS-S) 
     struct_replace_field, fields, 'FIELD2', fields.FIELD2, newtag='stack'      ; filename of stacked image for each field*** 
     struct_replace_field, fields, 'FIELD3', fields.FIELD3, newtag='mask'       ; filename of mask image for each field*** 
     struct_replace_field, fields, 'FIELD4', fields.FIELD4, newtag='pixelscale' ; pixel scale in arcseconds per pixel ("/px) for each field image 
     struct_replace_field, fields, 'FIELD5', fields.FIELD5, newtag='fieldinfo'  ; name of the catalog containing all the magnitude, RA, DEC, etc info of field sources 
     nfields = n_elements(fields.fieldname)                                     ; number of fields in the fields.csv file 
     ; *** note: to use, must be preceded by '/boomerang-data/alhall/GalaxiesInBlobs/Data/3DHST/3DHSTimages/'+ 

  ntotal = nblobs + nfields   ; total number of blobs AND fields 

END       ; end of readdata procedure 









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
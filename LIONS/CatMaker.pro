;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; CATMAKER.PRO 
;  Contains all modules used for making succinct catalogs for blobs and fields in CSV and/or FITS format (with tag names stored in the file for the latter).  
;  CONTENTS: 
;      PRGcatmaker------------------------procedure   (main module for PRG blobs) 
;      HSTcatmaker------------------------procedure   (main module for fields) 
;      CSVtoFITS--------------------------procedure   (main module for converting CSV files made by the preceding two modules into FITS files) 
;      blobcatreader----------------------procedure
;      fieldcatreader---------------------procedure
;      struct_replace_field---------------procedure 
;      fluxtomag--------------------------function 
;      rsex-------------------------------function 
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;




PRO PRGcatmaker 
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; PRGcatmaker procedure   
;   Creates succinct catalogs containing all relevant information about galaxies in images of PRG blobs. 
;
; INPUTS: blobs.csv - a CSV file containing all information regarding a list of PRG blobs, including the names of SExtractor catalogs in multiple filters 
;
; OUTPUT: creates one csv-format catalog per blob; these catalogs contain each galaxy's ID, RA, declination, and magnitudes in 3 bands (red, middle, blue) 
;                       
; Uses the rsex function.
;
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  FORWARD_FUNCTION rsex     ; this function is used to read in SExtractor catalogs 
  blobs = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/PRGcatalogsinfo.csv', n_table_header=2)   ; read in the file containing all blob information 
  ; assign names to each field of the "blobs" structure
     blobname = blobs.FIELD1       ; name of blob 
     blueband = blobs.FIELD2       ; name of reddest-filter's SExtractor catalog
     midband = blobs.FIELD3        ; name of middle-filter's SExtractor catalog
     redband = blobs.FIELD4        ; name of bluest-filter's SExtractor catalog
     nblobs = n_elements(blobname)  ; number of blobs in the PRGcatalogsinfo.csv file 

  FOR blob=0, nblobs-1 DO BEGIN     ; loop over the blobs to make one catalog per blob 
    ; make structures for each SExtractor catalog (one per each of 3 bands for this blob) 
       datared = rsex(redband[blob])      ; F140W
       datamid = rsex(midband[blob])      ; F814W
       datablue = rsex(blueband[blob])    ; either F475W (PRG1 and PRG3) or F606W (PRG2)
    ; get rid of sources that are stars or have flags 
       good = WHERE((datared.CLASS_STAR LT 0.8) AND (datared.IMAFLAGS_ISO EQ 0.))    ; find the indices of unflagged (good photometry) galaxies (no stars, no flags!) 
       datared = datared[good]    ; keep only good objects in the red catalog
       datamid = datamid[good]    ; keep only good objects in the middle catalog
       datablue = datamid[good]   ; keep only good objects in the blue catalog
    ; now combine the individual catalogs to make one succinct catalog
       catname = blobname[blob] + '.csv'                                  ; name the catalog after the blob 
       header = ['ID', 'RA', 'DEC', 'REDMAG', 'MIDMAG', 'BLUEMAG']        ; write a header labeling each field of the CSV file for easy reading 
       write_csv, catname, datared.NUMBER, datared.ALPHA_J2000, datared.DELTA_J2000, datared.MAG_ISO, datamid.MAG_ISO, datablue.MAG_ISO, header=header  ; make the csv catalog
  ENDFOR    ; all the blobs in blobs.csv

END      ; of PRGcatmaker procedure 









PRO HSTcatmaker, originalcat, zcat, filters, outputname
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; HSTcatmaker procedure
;   Creates a succinct catalog containing all relevant information about galaxies in a 3D-HST field. 
;
; INPUTS: originalcat - a string containing the name (with path starting from 3DHST folder) of the original 3D-HST catalog for the desired field with .fits extension
;         zcat - a string containing the name (with path starting from 3DHST folder) of the original 3D-HST redshift catalog for the desired field with .fits extension
;         filters - a four-element string array of the names of tags in the original catalog for magnitudes in bands we want to keep in order of [red, mid, blue, altblue]
;         outputname - a string containing the name of a csv file that will serve as the new catalog, with .csv extension 
;
; OUTPUT: creates a csv-format catalog containing each field galaxy's ID, RA, declination, magnitudes in 4 bands (1 red, 1 middle, 2 blue), and redshift
;         
; Uses the fluxtomag function.
; 
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  FORWARD_FUNCTION fluxtomag      ; this function will be used to convert the sources' fluxes into magnitudes

  ; READ IN THE DATA 
     field = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/Data/3DHST/'+originalcat, 1)          ; magnitude data for sources in the field 
     Zfield = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/Data/3DHST/'+zcat, 1)                ; redshift data for sources in the field 
     z = Zfield.z_best                                                                        ; define redshifts to use 

  ; FIND STRUCTURE TAGS OF MAGNITUDE BANDS TO KEEP 
     datatags = tag_names(field)
     redtag = WHERE(datatags EQ filters[0], /NULL)
     midtag = WHERE(datatags EQ filters[1], /NULL)
     bluetag = WHERE(datatags EQ filters[2], /NULL)
     altbluetag = WHERE(datatags EQ filters[3], /NULL)

  ; FILTER SOURCES to exclude ones that aren't usable 
   ; filter out everything that's flagged as having bad photometry (including objects that aren't galaxies) 
     use = WHERE(field.use_phot EQ 1.)                                     ; "use" is only indices of galaxies with good photometry
     field = field[use]                                                         ; trim field structure
     z = z[use]                                                                 ; trim z array 
   ; filter out galaxies that don't have measured fluxes in the bands we need 
     use = WHERE( (field.(redtag) GT 0.) AND (field.(midtag) GT 0.) )      ; red and middle bands; "use" is only galaxies with both red- and middle-filter fluxes measured
     field = field[use]                                                         ; trim field structure
     z = z[use]                                                                 ; trim z array 
     use = WHERE( (field.(bluetag) GT 0.) OR (field.(altbluetag) GT 0.) )  ; blue bands; "use" is only galaxies with one of the blue-filter fluxes measured
     field = field[use]                                                         ; trim field structure
     z = z[use]                                                                 ; trim z array 

  ; CONVERT FLUXES to magnitudes since magnitude is the quantity used in later analysis 
     fieldred = fluxtomag(field.(redtag))         ; convert the red band flux into magnitude 
     fieldmid = fluxtomag(field.(midtag))     ; convert the middle band flux into magnitude 
     fieldblue = fluxtomag(field.(bluetag))        ; convert the first blue band flux into magnitude 
     fieldbluealt = fluxtomag(field.(altbluetag))     ; convert the second blue band flux into magnitude 

  ; MAKE THE FINAL CATALOG containing only relevant information 
     red = filters[0]                                            ; name of the reddest-magnitude column in the catalog
     mid = filters[1]                                            ; name of the middle-magnitude column in the catalog
     blu = filters[2]                                            ; name of the bluest-magnitude column in the catalog
     alt = filters[3]                                            ; name of the alternate bluest-magnitude column in the catalog
     header = ['ID', 'ra', 'dec', red, mid, blu, alt, 'z']       ; write a header labeling each field of the CSV file for easy reading 
     filename = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+outputname     ; give the full path to where the output file will be saved 
     write_csv, filename, field.id, field.ra, field.dec, fieldred, fieldmid, fieldblue, fieldbluealt, z, header=header     ; make the csv catalog 

END      ; of HSTcatmaker procedure 









PRO CSVtoFITS
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; CSVtoFITS procedure
;  Converts catalogs of galaxies and their properties from CSV files into FITS files for both blobs and fields. 
;
; INPUT: none
;
; OUTPUT: Each blob and field gets its own FITS file, objectname.fits, containing all the information about the galaxies in its data image. 
;
; Uses the struct_replace_field procedure, the blobcatreader procedure, and the fieldcatreader procedure. 
;
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  files = read_csv('/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/CSVfiles.csv')       ; this file contains the names and object types (blob or field) of all the csv catalogs
  struct_replace_field, files, 'FIELD1', files.FIELD1, newtag='objname'                       ; name of object (PRG1, goodss, etc)
  struct_replace_field, files, 'FIELD2', files.FIELD2, newtag='objtype'                       ; type of object: blob or field 
  nobj = n_elements(files.objname)                                                            ; total number of objects for which there are csv catalogs to be converted 

  FOR obj=0, nobj-1 DO BEGIN                                                                  ; loop through every object with a csv catalog to be converted
    filename = files.objname[obj] + '.csv'                                                    ; get object's csv catalog name from its name 
    IF (files.objtype[obj] EQ 'blob') THEN BEGIN                                              ; if the object is a blob, 
      blobcatreader, filename, cat                                                              ; use blobcatreader to read its csv catalog and rename structure tags
    ENDIF ELSE IF (files.objtype[obj] EQ 'field') THEN BEGIN                                  ; if the object is a field, 
      fieldcatreader, filename, cat                                                              ; use fieldcatreader to read its csv catalog and rename structure tags
    ENDIF                                                                                     ; now we have a structure for the object whether it's a blob or a field 
    fullfilename = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/' + files.objname[obj] + '.fits'  ; make a filename by appending path and extension onto object name
    mwrfits, cat, fullfilename, /create                                                                               ; save the structure catalog as a FITS file
  ENDFOR     ; all the objects in CSVfiles.csv that had a csv catalog to be converted into a FITS file 

END     ; of CSVtoFITS procedure









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
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
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
     struct_replace_field, fieldcat, 'FIELD5', fieldcat.FIELD5, newtag='midmag'       ; the middle-band (F814W or F775W) magnitude of each galaxy
     struct_replace_field, fieldcat, 'FIELD6', fieldcat.FIELD6, newtag='bluemag'      ; the primary bluest-band (F606W) magnitude of each galaxy   
     struct_replace_field, fieldcat, 'FIELD7', fieldcat.FIELD7, newtag='altbluemag'   ; the secondary bluest-band (F435W) magnitude of each galaxy 
     struct_replace_field, fieldcat, 'FIELD8', fieldcat.FIELD8, newtag='z'            ; galaxy redshift 
END            ; of fieldcatreader procedure









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









FUNCTION fluxtomag, flux
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; fluxtomag function
;   Converts a 3D-HST flux into an AB magnitude. 
;
; INPUT: flux - the flux of one or more sources in 3D-HST (can be a float or a float array) 
;
; OUTPUT: magAB - the magnitude of the 3D-HST sources (has the same type and dimensionality as flux)  
;
; Uses no other functions nor procedures.   
;
; NOTES: This function is calibrated to the zeropoint for the 3D-HST survey. 
; 
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  magAB = 25.0 - 2.5*alog10(flux)     ; zeropoint is 25 so magAB = 25.0 - 2.5*log10(flux)
  RETURN, magAB 
END    ; of fluxtomag procedure









; NOTE: THE FOLLOWING FUNCTION WAS NOT WRITTEN BY AGNAR. HOWEVER, AGNAR DID EDIT IT.  
; 
FUNCTION rsex, catalog, use_row=use_row
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; rsex function 
;   Reads in arbitrary SExtractor format catalogs, using native header information & data themselves. Correctly reads longs, strings, and doubles.  
; 
; INPUT: catalog - a string containing the name of a SExtractor-format catalog
;        use_row - optional integer keyword specifying a row of the catalog to use to determine the format of the output data structure 
;                   (zero-indexed, default 0); this is necessary because occasionally a column in the first row is a different format 
;                   than the rest of the rows (for example, first row is LONG while rest of rows are DOUBLE)
;
; OUTPUT: cs - a structure with all catalog entries, using field names for tagnames  
;
; Uses the valid_num procedure, which is already known by IDL. 
;
; NOTES: 
;  SYNTAX: cs = rsex('catalog.cat')
;  CREDIT: L. Moustakas (LAM) & J. Moustakas (JM); see full modification history at https://github.com/moustakas/impro/blob/master/pro/sextractor/rsex.pro
;  ISSUES: - uses str_sep, which according to Harris Geospatial is now obsolete and should be replaced with strsplit 
;          - there are two lines, each marked with 10 asterisks at the start of its comment, where Agnar added a /get_lun keyword that wasn't in the original code  
;
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  ; First, make sure everything has been called correctly 
     ; Check that an argument has been passed:     
        IF (n_params() LE 0) THEN BEGIN                ; if user gave no arguments
          print,'syntax: cat=rsex(catalog)'               ; tell the user what the correct syntax for this function is
          RETURN, -1                                      ; return a -1 as an error and exit this function
        ENDIF                                          
     ; See if the catalog file exists at all:
        IF (file_test(catalog,/regular) EQ 0L) THEN BEGIN                  ; if the user input a nonexistent catalog
          print, 'error - the specified catalog does not seem to exist'       ; tell them so 
          RETURN,-1                                                           ; return a -1 as an error and exit this function
        ENDIF


  ; Count catalog lines
     spawn, 'wc -l < ' + catalog, nlines, /sh      ; tell UNIX to count the number of lines in the input SExtractor catalog
     nlines = long(nlines[0])                      ; convert the UNIX output into a long integer to use here in IDL 


  ; Open the catalog to read
     openr, lun ,catalog, error=err, /get_lun                   ; open the catalog
     IF (err NE 0) THEN BEGIN                                   ; if there was an error in opening the catalog
       message, 'error opening catalog - ' + !error_state.msg      ; give the user an error message saying what the error was 
       free_lun, lun                                               ; release the lun used for this catalog
       RETURN, -1                                                  ; return a -1 as an error and exit this function
     ENDIF                                                      ; the catalog had errors 


  ; Determine which row of catalog to use in formatting the output structure
     IF (n_elements(use_row) EQ 0L) THEN use_row = 0L           ; if the use_row keyword wasn't given, use the first row of the catalog
     ; otherwise, the use_row keyword already says which row to use for formatting 
      

  ; Read in the header
     head = ''           ; create an empty string for names of data columns; this variable will become an array
     numb = 0l           ; empty long integer containing the tag number of each tag in "head"; this variable will become an array
     cstr=''             ; create an empty string that will contain each string line of the catalog as we loop through them   
     nheadcomment = 0L   ; empty long integer containing the number of rows in the catalog that are reserved just for header comments 
     tag = 1             ; this is a keyword for the WHILE loop in the next line, indicating that the line being read is a string and not data (ie, line contains '#')  
     WHILE tag DO BEGIN                                              ; while '#' tag is true
       readf, lun, cstr                                                 ; read a line of the file into cstr
       IF (strmatch(cstr,'##*') EQ 0B) THEN BEGIN                       ; if this line isn't a comment line (comment lines start with '##')
         IF (strpos(cstr,'#') NE -1) THEN BEGIN                            ; if this line of the catalog is a tag name (ie, has a '#' in it) 
           head = [head, strupcase((str_sep(strcompress(cstr),' '))[2])]      ; add the tag name to the "head" array
           numb = [numb, (str_sep(strcompress(cstr),' '))[1]]                 ; add the tag number to the "numb" array
         ENDIF ELSE tag=0                                                  ; if the line doesn't have a '#' in it, it's a line of data and not a string, so exit the WHILE loop
       ENDIF ELSE nheadcomment = nheadcomment + 1L                      ; skip comment lines; increment nheadcomment to reflect that this was a comment line
     ENDWHILE                                                        ; all the string lines at the start of the catalog have been read in 
     nhead = n_elements(head)   ; number of header entries
     head = head[1:(nhead-1)]   ; get rid of the first element of the head array (the empty string) 
     numb = numb[1:(nhead-1)]   ; get rid of the first element of the numb array (the empty long) 
     nhead = n_elements(head)   ; update the number of header entries to reflect that the empty string is gone now 
     free_lun, lun              ; release the lun used for this catalog 


  ; Check header for implicit array entries, and expand relevant tagnames
     tothead=''                      ; make an empty string for the TOTAL number of tag names, including dimensions of arrays 
     ; First, loop through all the fields except the last one to check if each is actually an array 
        FOR i=0l, nhead-2 DO BEGIN      ; go through the "head" array from the first entry to the second-to-last entry 
          tothead = [tothead,head[i]]      ; add this element of "head" onto "tothead" 
          mlen = numb[i+1]-numb[i]-1       ; determine whether the next entry in "numb" is more than 1 bigger than this entry of "numb" 
          IF (mlen NE 0) THEN BEGIN        ; this means that this entry of "head" is actually an array 
            FOR j=1, mlen DO tothead=[tothead,head[i]+strcompress(j,/remove_all)]   ; add all the entries in this array within "head" onto "tothead" 
          ENDIF                            ; mlen wasn't zero 
        ENDFOR                             ; every element in "head" except the last one 
        ntothead = n_elements(tothead)                         ; get the total number of strings in the header
        tothead = [tothead[1:(ntothead-1)], head[(nhead-1)]]   ; add the last element of "head" onto "tothead" and remove the empty string at the beginning of "tothead" 
     ; Now check whether the last header field is actually the first one of an array
       ; Read the first data line; if it has more entries than the total number of header fields so far, we've missed an array at the end 
          junk = strarr(nhead+nheadcomment)  ; create a string array containing the same number of entries as all the header lines 
          openr, lun, catalog, /get_lun      ; open the catalog again
          readf, lun, junk                   ; read the header lines into the "junk" string array
          readf, lun, cstr                   ; read the next line of the catalog into "cstr"
          free_lun, lun                      ; release the lun used for this catalog
          nstr = n_elements(str_sep(strcompress(strtrim(cstr,2)),' '))    ; see how many entries are in cstr (how many columns of data are in the catalog) 
       ; make sure the number of entries in cstr matches ntothead
          IF (nstr GT ntothead) THEN BEGIN       ; if there are more elements in cstr than in tothead, then last entry of tothead is actually an array 
            mlen = nstr - ntothead                                                          ; figure out how many more elements cstr has than tothead
            FOR j=1, mlen DO tothead=[tothead,head[(nhead-1)]+strcompress(j,/remove_all)]   ; add all the entries in this array within "head" onto "tothead"
          ENDIF                                  ; now tothead contains one entry for every column of the catalog 
          ntothead=n_elements(tothead)           ; update ntothead since we've added elements to tothead 


  ; That should do it! Now ntothead contains one element per data column. Onwards, to read the data...
     nbody = nlines - nhead - nheadcomment           ; get the number of lines that are actual data rather than header stuff
     junk=strarr(nhead+nheadcomment)                 ; create a string array containing the same number of entries as all the header lines 
     openr, lun, catalog, /get_lun                   ; **********open the catalog to read it in   
     readf, lun, junk                                ; read all the header stuff into "junk"
     FOR ijunk = 0L, use_row-1L DO readf,lun,junk2   ; read everything before "use_row" into another junk array
     readf, lun, cstr                                ; read the info from the "use_row" row into cstr
     free_lun, lun                                   ; release the lun used for this catalog so we can start from the beginning of the data to actually read it in
     body=strarr(nbody)                              ; make a string array with the same number of elements as there are data lines in the catalog 
     openr, lun, catalog, /get_lun                   ; **********open the catalog again    
     readf, lun, junk                                ; read all the header stuff into "junk" 
     readf, lun, body                                ; read the catalog data into "body" 
     free_lun, lun                                   ; release the lun used for this catalog 
     bodyslam=strarr(ntothead,nbody)                 ; make a 2D string array with same number of rows as the number of tags and same number of columns as lines of data
     FOR i=0l, nbody-1 DO BEGIN                      ; loop through all the data lines 
       bodyslam[*,i]=(str_sep(strcompress(strtrim(body[i],2)),' '))[0L:ntothead-1L]   ; separate this line's string of data into individual entries and fill in the row with it 
     ENDFOR                                            ; every line of data has been separated out into one string per tag 
     carr = str_sep(strcompress(strtrim(cstr,1)),' ')  ; divide the formatting line into individual entries too 


  ; make a mini-structure using the "use_row" row to get the format for every field of the structure 
     tind = valid_num(carr)             ; array that checks to see if carr can be represented as a number; 0=string, 1=number
     tindint = valid_num(carr,/int)     ; array that checks to see if carr can be represented as an integer specifically; 0=non-integer, 1=integer
     ; Now use the above to figure out the data type: string, integer, or double (any float will be a double, any integer will be a long) 
        IF (tind[0] EQ 0) THEN type='' ELSE IF (tindint[0] EQ 1) THEN type=0l ELSE type=0.d0       ; get the data type for the first tag 
        cs=create_struct(idl_validname(tothead[0],/convert_all),type)                          ; make a (mini-)structure containing the first tag name and specify its data type 
        FOR i=1l,ntothead-1 DO BEGIN                                                           ; now loop through the rest of the tags
          IF (tind[i] EQ 0) THEN type='' ELSE IF (tindint[i] EQ 1) THEN type=0l ELSE type=0.d0     ; get the data type for this tag 
          cs=create_struct(cs,idl_validname(tothead[i],/convert_all),type)                         ; add this tag and its data type onto the structure
        ENDFOR                                                                                 ; all the tags in ntothead 


  ; Use the mini-structure created above to make the full structure and fill it in with the data from the catalog
     cs=replicate(cs,nbody)            ; replicate the (mini-)structure to make the main structure containing all the data in the catalog 
     FOR i=0l, ntothead-1 DO BEGIN     ; loop through every tag 
       FOR j=0l, nbody-1 DO BEGIN         ; loop through every line of data 
         cs[j].(i)=bodyslam[i,j]             ; fill in this mini-structure's i-th tag with the data, which was stored in "bodyslam" 
       ENDFOR                             ; every line of data (ie, every mini-structure in the main structure) 
     ENDFOR                            ; every tag of the main structure, cs 
     

  ; That's that - all the data in the catalog has been read into the "cs" structure with appropriate tag names, so we can return the structure and exit the function.  
     RETURN, cs      ; return the structure containing everything that was in the SExtractor catalog

END      ; of rsex function

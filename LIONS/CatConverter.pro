;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; CATCONVERTER.PRO
;  Contains all modules used for converting CSV-format galaxy catalogs in both blobs and fields into FITS-format catalogs with tag names stored in the file. 
;  CONTENTS:
;      CSVtoFITS---------------------------procedure (main module)
;      blobcatreader-----------------------procedure
;      fieldcatreader----------------------procedure
;      struct_replace_field----------------procedure 
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;




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

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; CATMAKER.PRO 
;  Contains all modules used for making succinct catalogs for blobs and fields.  
;  CONTENTS: 
;      blobcatmaker-----------------------procedure   (main module for blobs) 
;      HSTcatmaker------------------------procedure   (main module for fields) 
;      fluxtomag--------------------------function 
;      rsex-------------------------------function 
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;




PRO blobcatmaker 
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; blobcatmaker procedure   
;   Creates succinct catalogs containing all relevant information about galaxies in images of blobs. 
;
; INPUTS: blobs.csv - a CSV file containing all information regarding a list of blobs, including the names of SExtractor catalogs in multiple filters 
;
; OUTPUT: creates one csv-format catalog per blob; these catalogs contain each galaxy's ID, RA, declination, and magnitudes in 3 bands (red, middle, blue) 
;                       
; Uses the rsex function.
;
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  FORWARD_FUNCTION rsex     ; this function is used to read in SExtractor catalogs 
  blobs = read_csv('blobs.csv', n_table_header=1)   ; read in the file containing all blob information 
  ; assign names to each field of the "blobs" structure
     blobname = blobs.FIELD01                                                                                       ; name of blob 
     blueband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobCatalogs/'+blobs.FIELD04   ; name of reddest-filter's SExtractor catalog
     midband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobCatalogs/'+blobs.FIELD05    ; name of middle-filter's SExtractor catalog
     redband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobCatalogs/'+blobs.FIELD06    ; name of bluest-filter's SExtractor catalog
     nblobs = n_elements(blobname)                                                                                  ; number of blobs in the blobs.csv file 

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

END      ; of blobcatmaker procedure 









PRO HSTcatmaker 
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; HSTcatmaker procedure
;   Creates a succinct catalog containing all relevant information about galaxies in a 3D-HST field. 
;
; INPUTS: originalcat - a string containing the name (with path starting from 3DHST folder) of the original 3D-HST catalog for the desired field 
;         zcat - a string containing the name (with path starting from 3DHST folder) of the original 3D-HST redshift catalog for the desired field
;         outputname - a string containing the name of a csv file that will serve as the new catalog
;
; OUTPUT: creates a csv-format catalog containing each field galaxy's ID, RA, declination, magnitudes in 4 bands (1 red, 1 middle, 2 blue), and redshift
;         
; Uses the fluxtomag function.
; 
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
  FORWARD_FUNCTION fluxtomag      ; this function will be used to convert the sources' fluxes into magnitudes

  ; USER-DEFINED VARIABLES 
  originalcat = 'GOODSS/Catalog/goodss_3dhst.v4.1.cat.FITS'
  zcat = 'Z_GOODSS/goodss_3dhst.v4.1.5.zbest.fits'
  outputname = 'goodss.csv'

  ; READ IN THE DATA 
     field = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/3DHST/'+originalcat, 1)          ; magnitude data for sources in the field 
     Zfield = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/3DHST/'+zcat, 1)                ; redshift data for sources in the field 
     z = Zfield.z_best                                                                        ; define redshifts to use 

  ; FILTER SOURCES to exclude ones that aren't usable 
   ; filter out everything that's flagged as having bad photometry (including objects that aren't galaxies) 
     use = WHERE(field.use_phot EQ 1.)                                   ; "use" is only indices of galaxies with good photometry
     field = field[use]                                                     ; trim field structure
     z = z[use]                                                             ; trim z array 
   ; filter out galaxies that don't have measured fluxes in the bands we need 
     use = WHERE( (field.f_F140W GT 0.) AND (field.f_F814Wcand GT 0.) )  ; red and middle bands; "use" is only galaxies with both red- and middle-filter fluxes measured
     field = field[use]                                                     ; trim field structure
     z = z[use]                                                             ; trim z array 
     use = WHERE( (field.f_F606W GT 0.) OR (field.f_F435W GT 0.) )       ; blue bands; "use" is only galaxies with one of the blue-filter fluxes measured
     field = field[use]                                                     ; trim field structure
     z = z[use]                                                             ; trim z array 

  ; CONVERT FLUXES to magnitudes since magnitude is the quantity used in later analysis 
     fieldred = fluxtomag(field.f_F140W)         ; convert the red band flux into magnitude 
     fieldmid = fluxtomag(field.f_F814Wcand)     ; convert the middle band flux into magnitude 
     fieldblue = fluxtomag(field.f_F606W)        ; convert the first blue band flux into magnitude 
     fieldbluealt = fluxtomag(field.f_F435W)     ; convert the second blue band flux into magnitude 

  ; MAKE THE FINAL CATALOG containing only relevant information 
     header = ['ID', 'ra', 'dec', 'F140W','F814W','F606W', 'F435W','z']                           ; write a header labeling each field of the CSV file for easy reading 
     filename = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/individual/'+outputname     ; give the full path to where the output file will be saved 
     write_csv, filename, field.id, field.ra, field.dec, fieldred, fieldmid, fieldblue, fieldbluealt, z, header=header     ; make the csv catalog 

END      ; of HSTcatmaker procedure 









  FUNCTION fluxtomag, flux
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;
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
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;
    magAB = 25.0 - 2.5*alog10(flux)     ; zeropoint is 25 so magAB = 25.0 - 2.5*log10(flux)
    return, magAB 
  END    ; of fluxtomag procedure









  ; NOTE: THE FOLLOWING FUNCTION WAS NOT WRITTEN BY AGNAR. 
  ;
  ;
  ; NAME:
  ;       RSEX()
  ;
  ; PURPOSE:
  ;       Read in arbitrary SExtractor format catalogs, using native
  ;       header information & data themselves.  Correctly reads longs,
  ;       strings, and doubles.  
  ;
  ; INPUTS:
  ;       A SExtractor-format catalog
  ;
  ; OPTIONAL INPUTS:
  ;       use_row - use this row of the catalog to determine the format
  ;                 of the output data structure (zero-indexed; default
  ;                 0); this is necessary because occasionally a column
  ;                 in the first row is a different format (e.g., LONG)
  ;                 than the rest of the rows (e.g., DOUBLE)
  ;
  ; KEYWORD PARAMETERS:
  ;
  ; OUTPUTS:
  ;       Returns a structure with all catalog entries, using field
  ;       names for tagnames. 
  ;
  ; COMMON BLOCKS:
  ;       None.
  ;
  ; RESTRICTIONS:
  ;       None.
  ;
  ; PROCEDURE:
  ;       Use syntax
  ;       cs = rsex('catalog.cat')
  ;
  ; COMMENTS:
  ;
  ; PROCEDURES USED:
  ;       FINDFILE
  ;       CREATE_STRUCT
  ;       VALID_NUM
  ;       VALID_NUM_ARR
  ;
  ; MODIFICATION HISTORY:
  ; LAM = L. Moustakas
  ; JM  = J. Moustakas
  ;       LAM '04may02 - fixed problem case of there being an array of
  ;                      values based on the last header tag position.
  ;                      will now work with that case after special
  ;                      check. 
  ;       LAM '04feb04 - converted from the old lrsex.pro; adapted to
  ;                      detect and read longs, strings, and doubles.  
  ;       JM  '04sep22 - skip comment lines in the header designated
  ;                      with a double-hash (##)
  ;       JM  '07aug18 - added USE_ROW optional input
  ;       JM  '08jul25 - use a variable logical unit number
  ;       LAM+JM '08dec11 - vetted - no changes made
  ;       jm14jul12siena - use IDL_VALIDNAME() to ensure valid structure tags
  ;-
  FUNCTION rsex,catalog, use_row=use_row

  ; Check that an argument has been passed     
      IF n_params() LE 0 THEN BEGIN 
	  print,'cat=rsex(catalog)'
	  return,-1
      ENDIF 

  ; See if the catalog file exists at all

      if (file_test(catalog,/regular) eq 0L) then begin ; jm04sep20uofa
	  print,'error - the specified catalog does not seem to exist'
	  return,-1
      endif

  ;    jnk = findfile(catalog,count=catexist)
  ;    if catexist eq -1 then begin
  ;        message,'error - the specified catalog does not seem to exist'
  ;        return,-1
  ;    endif

  ; Count catalog lines
      spawn,'wc -l < '+catalog, nlines, /sh
      nlines = long(nlines[0]);-1

  ; Open the catalog to read
      openr,lun,catalog,error=err,/get_lun
      if err ne 0 then begin
	  message,'error opening catalog - '+!error_state.msg ; jm04sep20uofa
  ;        message,'error opening catalog - ',!error_state.msg
	  free_lun, lun
	  return,-1
      endif

      if (n_elements(use_row) eq 0L) then use_row = 0L ; jm07aug18nyu
      
  ; Read in the header
      head = ''
      numb = 0l
      cstr=''
      tag = 1
      nheadcomment = 0L ; jm04sep20uofa
      while tag do begin         ; while '#' tag is true
	  readf,lun,cstr
	  if (strmatch(cstr,'##*') eq 0B) then begin ; skip comment lines
	      if strpos(cstr,'#') ne -1 then begin
		head = [head, strupcase((str_sep(strcompress(cstr),' '))[2])]
		numb = [numb, (str_sep(strcompress(cstr),' '))[1]]
	      endif else tag=0
	  endif else nheadcomment = nheadcomment + 1L ; jm04sep20uofa
      endwhile
      nhead = n_elements(head)   ; number of header entries
      head = head[1:(nhead-1)]
      numb = numb[1:(nhead-1)]
      nhead = n_elements(head)
      free_lun, lun

  ; Check header for implicit array entries, and expand relevant tagnames
      tothead=''
      for i=0l,nhead-2 do begin
	  tothead = [tothead,head[i]]
	  mlen = numb[i+1]-numb[i]-1
	  if mlen ne 0 then $
	    for j=1,mlen do $
	    tothead=[tothead,head[i]+strcompress(j,/remove_all)]
      endfor

  ; The header info so far
      ntothead = n_elements(tothead)
      tothead=[tothead[1:(ntothead-1)],head[(nhead-1)]]

  ; Need to check one more thing -- whether the last header field is
  ; actually the first one of an array...

  ; Read the first data line.  If it has more entries than the total
  ; number of header fields so far, we've missed an array at the end. 
      junk=strarr(nhead+nheadcomment) ; jm04sep20uofa
      openr,lun,catalog,/get_lun
      readf,lun,junk
      readf,lun,cstr
      free_lun, lun
      nstr = n_elements(str_sep(strcompress(strtrim(cstr,2)),' '))

      if nstr gt ntothead then begin
	  mlen = nstr - ntothead
	  for j=1,mlen do $
	    tothead=[tothead,head[(nhead-1)]+strcompress(j,/remove_all)]
      endif
      ntothead=n_elements(tothead)

  ; That should do it!  Onwards, to read the data...
      nbody=nlines-nhead-nheadcomment ; jm04sep20uofa
      junk=strarr(nhead+nheadcomment) ; jm04sep20uofa
      openr,lun,catalog
      readf,lun,junk
      for ijunk = 0L, use_row-1L do readf,lun,junk2 ; jm07aug18nyu
      readf,lun,cstr
      free_lun, lun
      body=strarr(nbody)
      openr,lun,catalog
      readf,lun,junk
      readf,lun,body
      free_lun, lun
      bodyslam=strarr(ntothead,nbody)
      for i=0l,nbody-1 do begin
  ;       print, i, body[i]
  ;       bodyslam[*,i]=(str_sep(strcompress(strtrim(body[i],2)),' '))[0L:nhead-1L]
	  bodyslam[*,i]=(str_sep(strcompress(strtrim(body[i],2)),' '))[0L:ntothead-1L] ; jm04nov22uofa
      endfor
      carr=str_sep(strcompress(strtrim(cstr,1)),' ')

      tind    = valid_num(carr) ; 0=string 1=number
      tindint = valid_num(carr,/int) ; 0=non-integer 1=integer
      if tind[0] eq 0 then $
	type='' else $
	if tindint[0] eq 1 then $
	type=0l else $
	type=0.d0
      cs=create_struct(idl_validname(tothead[0],/convert_all),type)

      for i=1l,ntothead-1 do begin
	  if tind[i] eq 0 then $
	    type='' else $
	    if tindint[i] eq 1 then $
	    type=0l else $
	    type=0.d0
	  cs=create_struct(cs,idl_validname(tothead[i],/convert_all),type)
      endfor
      cs=replicate(cs,nbody)

      for i=0l,ntothead-1 do begin
	  for j=0l,nbody-1 do begin
	    cs[j].(i)=bodyslam[i,j]
  ;          print, bodyslam[i,j], i, j
  ;          print, i, j
	  endfor
  ;stop
      endfor
      
      return,cs
  END 
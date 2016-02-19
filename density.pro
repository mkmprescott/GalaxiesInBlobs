PRO density   ;, seed, nApertures

FORWARD_FUNCTION get_galaxies
FORWARD_FUNCTION rsex

;data = rsex("prg1_radectest.cat")     ; this is reading everything in as one giant field...

data = read_csv("prg1_rdtest.cat", n_table_header=16)
; will later modify this to be more automated

; to make mag bins: make "new" catalogs! 
;  make a for loop with i=minimum mag, maximum mag, binsize
;   make a list or file or something for galaxies in this binsize
;   run through original catalog and find galaxies with i <= mag <= i+binsize
;   add those galaxies to the list or file or whatever
;  now you have a number of mini catalogs for mags
;  then you can do the following for each mag catalog

blobra = 218.8020833   ; right ascension of blob PRG1 in decimal degrees
blobdec = 35.18558333  ; declination of blob PRG1 in decimal degrees
; I think I'm going to make a data file that contains the RA and dec of all our blobs 
;   and use that for this part - later

; get blob density: number of galaxies inside aperture centered on blob (radius of say, 10 arcsec or so)
aperture_radius_blob = 10./3600.   ; aperture with radius of 10 arcsec in decimal degrees
blob_galaxies=get_galaxies(blobra, blobdec, aperture_radius_blob, data)
blob_ngal = n_elements(blob_galaxies)
; now get density in number of galaxies per square arcsecond 
blobdensity = float(blob_ngal) / (!dPI*(10./3600.)^2.)   ; the denominator is the area of the aperture around the blob, which has 10" radius  

print, blobdensity, " galaxies per square degree in blob region"  
print, blob_ngal, " galaxies in the blob region"

; !!!!!!! CODE PRINTS:     0.0159155 galaxies per square arcsecond in blob region
;  that seems low...
;  REVISED: code prints
;       206264.81 galaxies per square degree in blob region
;       5 galaxies in the blob region








; later stuff:

;aperture_ra_list = list()     ; this will contain RAs of apertures
;aperture_dec_list = list()    ; this will contain declinations of apertures
;aperture_radius_list = list() ; this will contain aperture radii 
;aperture_ngal_list = list()   ; this will contain # galaxies in each aperture

;FOR i=0, nApertures DO BEGIN
;;  place an aperture randomly and add its properties to the lists of aperture properties
;  aperture_radius = 10.*randomu(seed)
;  aperture_radius_list.Add, aperture_radius
;  aperture_ra = (minimum ra of image) + (max ra - min ra)*randomu(seed)
;  aperture_ra_list.Add, aperture_ra
;  aperture_dec = (minimum dec of image) + (max dec - min dec)*randomu(seed)
;  ; figure out the number of galaxies in the aperture 
;  galaxies = get_galaxies(aperture_ra, aperture_dec, aperture_radius, data)
;  ngal = n_elements(galaxies)
;  ; NOTE: ALSO OUTPUT THE WHOLE LIST OF GALAXIES!!! FOR MAG STUFF TO BE ADDED SOON
;  ; append # of galaxies in the aperture to list of those
;  aperture_ngal_list.Add, ngal
;ENDFOR
;
;; turn all the lists into arrays to work with them now that they're finished being constructed
;aperture_ra_array = aperture_ra_list.ToArray(Type=5)
;aperture_dec_array = aperture_dec_list.ToArray(Type=5)
;aperture_radius_array = aperture_radius_list.ToArray(Type=5)
;aperture_ngal_array = aperture_ngal_list.ToArray(Type=5)
;
;; compute density PER MAG BIN!! FIX THIS!!
;densityarray = aperture_ngal_array / ( !PI*(aperture_radius_array)^2. )    ; array of the densities in each aperture
;field_density = total(densityarray) / double(n_elements(densityarray))     ; this is the field density for the whole image

; make histogram of densities (number of apertures in each density bin on y axis) & oplot blob density
; make plot (galaxy magnitude and density)

; make sure to add IF statement accounting for edges of picture
; think about rejecting things that are on the blob

END





FUNCTION get_galaxies, ra, dec, radius, data
  distance = sqrt(  ( (ra - data.FIELD15) * cos(dec*!dPI/180.) )^2. + (dec - data.FIELD16)^2. )
  galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)
  RETURN, galaxies_in_aperture
END



;+
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
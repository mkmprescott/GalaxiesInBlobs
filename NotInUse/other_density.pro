PRO density
;--------------------------------------------------------------------------------------------------
;  some day this will include a nice thorough header and description... but not today
;--------------------------------------------------------------------------------------------------
FORWARD_FUNCTION get_galaxies
FORWARD_FUNCTION rsex
FORWARD_FUNCTION get_galaxies_binned
FORWARD_FUNCTION Poisson_error
FORWARD_FUNCTION get_galaxies_noblob


datared = rsex("setestF140W.cat")    
datamid = rsex("setestF814W.cat")
datablue = rsex("setestF475W.cat")  
; will later modify this to be more automated

; HARD-CODED PARAMETERS: 
blobra = 218.8020833   ; right ascension of blob PRG1 in decimal degrees
blobdec = 35.18558333  ; declination of blob PRG1 in decimal degrees
; I'm going to make a data file that contains the RA and dec of all blobs and use that for this part later on
ap_radius_arcsec = 10.               ; aperture with radius of 10 arcsec, in arcsec 
ap_radius = ap_radius_arcsec/3600.   ; aperture with radius of 10 arcsec in decimal degrees
; mag binning stuff:
binsize=0.5
brightestmag = 22.
dimmestmag = 35.
nbins = (dimmestmag - brightestmag + binsize)/binsize
mags = brightestmag + binsize*findgen(nbins)
; random aperture stuff:
nApertures = 100.
image = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Images/stackedVJH.fits', 0, hdr)
bpm = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Images/stackedVJH_bpm.fits')
pixelscale = 0.06   ; arceseconds per pixel
ap_radius_pix = ap_radius_arcsec / pixelscale
sizeinfo = size(image)
npix_x = sizeinfo[1]
npix_y = sizeinfo[2]
magrange=[brightestmag,dimmestmag]
plotsym, 0, 1, /fill         ; makes the psym=8 symbol a filled circle for plotting


; First, a test: get raw blob density, ie, number of galaxies inside aperture centered on blob with radius specified above with no mag bins yet
blob_galaxies_all=get_galaxies(blobra, blobdec, ap_radius, datared)
blob_ngalraw = n_elements(blob_galaxies_all)
; now get density in number of galaxies per square arcsecond 
blobdensityraw = float(blob_ngalraw) / (!dPI*ap_radius^2.)    
print, blobdensityraw, " galaxies per square degree in blob region"  
print, blob_ngalraw, " galaxies in the blob region"
;  CODE PRINTS:   
;  2722695.4 galaxies per square degree in blob region
;  66 galaxies in the blob region 

window, 31, retain=2, xsize=1200, ysize=1000
redmags = datared[blob_galaxies_all].MAG_ISO
midmags = datamid[blob_galaxies_all].MAG_ISO
bluemags = datablue[blob_galaxies_all].MAG_ISO
; test to make sure this is working correctly: 
;print, n_elements(redmags), n_elements(midmags), n_elements(bluemags), n_elements(blob_galaxies_all)
; code prints 66 for all four - works fine 
color_x = bluemags - midmags
color_y = midmags - redmags
plot, color_x, color_y, title="PRG1 Color-Color Diagram", xtitle="F475W - F814W", ytitle="F814W - F140W", background=255, color=0, charsize=1.5, psym=8, xrange=[-5,5], /xstyle, yrange=[-5,5], /ystyle
field_galaxies = get_galaxies_noblob(blobra, blobdec, ap_radius, datared)
field_redmags = datared[field_galaxies].MAG_ISO
field_midmags = datamid[field_galaxies].MAG_ISO
field_bluemags = datablue[field_galaxies].MAG_ISO
field_color_x = field_bluemags - field_midmags
field_color_y = field_midmags - field_redmags
oplot, field_color_x, field_color_y, color=0, psym=3
LEGEND, ['blob galaxy','field galaxy'], /left, /bottom, color=0, textcolor=0, psym=[8,3], charsize=1., /box, outline_color=0 
;dist = sqrt(  ( (blobra - datared.ALPHA_J2000) * cos(blobdec*!dPI/180.) )^2. + (blobdec - datared.DELTA_J2000)^2. )
;color_x_all = datablue.MAG_ISO - datamid.MAG_ISO
;weird = WHERE (((dist LT ap_radius) AND (datared.CLASS_STAR LT 0.8) AND (datared.IMAFLAGS_ISO EQ 0.) AND (color_x_all GT 2.)), n_galaxies)
;print, datared[weird].ALPHA_J2000, datared[weird].DELTA_J2000
;STOP
write_png, "PRG1color-color.png", tvrd(/true)




; mag binning: 

galaxy_maglist = list()
FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
  blob_galaxies_binned=get_galaxies_binned(blobra, blobdec, ap_radius, datared, mag, mag+binsize)
  galaxy_maglist.Add, blob_galaxies_binned
ENDFOR      ; Now we have a list containing arrays of the IDs of every galaxy in each mag bin, one array per bin
 
blob_ngal_binned = dblarr(nbins)
blob_ngalerr_binned = dblarr(nbins) 

FOR i=0, nbins-1. DO BEGIN
  bin = galaxy_maglist[i]
  IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
  blob_ngal_binned[i] = ngal 
  blob_ngalerr_binned[i] = Poisson_error(ngal) 
  print, blob_ngal_binned[i], " +- ", blob_ngalerr_binned[i], " galaxies in bin", i+1   ; just a test to make sure it's working
ENDFOR 
print, total(blob_ngal_binned), " total galaxies"   ; just a test to make sure it's working - this should match the previous no-bin number 

; check with user whether to make a plot or not since this plot isn't completely necessary
makeplot = ''
READ, makeplot, PROMPT='Make plot of number of galaxies in blob region? (y/n)'
IF (makeplot EQ 'y') THEN BEGIN 
  ; Now, plot up number in each bin with errors:
  window, 0, retain=2, xsize=1200, ysize=1000
  plot, mags, blob_ngal_binned, title="Galaxies in Blob Region", xtitle="magnitude (F140W)", ytitle="number of galaxies", background=255, color=0, charsize=1.5, psym=10, yrange=[0,15], /ystyle, xrange=magrange, /xstyle 
  errplot, mags, blob_ngal_binned-blob_ngalerr_binned, blob_ngal_binned+blob_ngalerr_binned, color=0
  write_png, "blobngal.png", tvrd(/true)
ENDIF 

; Finally, compute density and errors: 
densities = double(blob_ngal_binned)/(!dPI*ap_radius^2.) 
density_errors = double(blob_ngalerr_binned)/(!dPI*ap_radius^2.)

; check with user whether to make a plot or not since this plot isn't completely necessary
makeplot = ''
READ, makeplot, PROMPT='Make plot of density of galaxies in blob region? (y/n)'
IF (makeplot EQ 'y') THEN BEGIN 
  ; plot density with errors:
  window, 1, retain=2, xsize=1200, ysize=1000
  plot, mags, densities, title="Galaxies in Blob Region", xtitle="magnitude (F140W)", ytitle="N!Igal!N per deg!E2!N per 0.5 mag", background=255, color=0, charsize=1.5, psym=-8, /ylog, yrange=[1d4,1d6], /ystyle, xrange=magrange, /xstyle 
  errplot, mags, densities-density_errors, densities+density_errors, color=0
  write_png, "blobdensity.png", tvrd(/true)
ENDIF  




; SAMPLING THE FIELD

ap_ras = dblarr(nApertures)     ; this will contain RAs of all apertures
ap_decs = dblarr(nApertures)    ; this will contain declinations of all apertures
ap_ngal_binned = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
ap_ngalerr_binned = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
;ap_galaxy_list = list()     ; this will contain the IDs of all galaxies in each aperture 
megalist = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 


FOR ap=0, nApertures-1 DO BEGIN
  ; place an aperture randomly and add its properties to the lists of aperture properties
  ; aperture_radius = 10.*randomu(seed)
  ; aperture_radius_list.Add, aperture_radius
  repeatflag = 0.
  WHILE (repeatflag EQ 0.) DO BEGIN
    ap_x = double(npix_x)*randomu(seed)
    ap_y = double(npix_y)*randomu(seed)

    dist_ellipse, dim, [npix_x, npix_y], ap_x, ap_y, 1.,0.
;    apim = where(dim LT ap_radius_pix)
    newap = dim
    FOR i=0, npix_x-1 DO BEGIN
      FOR j=0, npix_y-1 DO BEGIN
        IF (newap[i,j] LT ap_radius_pix) THEN BEGIN
            newap[i,j] = 1.
        ENDIF ELSE newap[i,j] = 0.
      ENDFOR
    ENDFOR
    apbpm = bpm*newap
    badpix = where(apbpm GT 0, nbadpix)
    ntotalpix = n_elements(apbpm)
    IF (nbadpix LT ntotalpix) THEN BEGIN 
      repeatflag = 1. 
    ENDIF ELSE repeatflag = 0.

    xyad, hdr,ap_x, ap_y, ap_ra, ap_dec
    IF (sqrt(  ( (blobra - ap_ra) * cos(blobdec*!dPI/180.) )^2. + (blobdec - ap_dec)^2. ) GT 2.*ap_radius) THEN BEGIN
      repeatflag = 1.
    ENDIF ELSE repeatflag = 0.
  ENDWHILE

  ap_area_weight = (float(ntotalpix) - float(nbadpix)) / float(ntotalpix) 
  xyad, hdr,ap_x, ap_y, ap_ra, ap_dec
  ap_ras[ap] = ap_ra
  ap_decs[ap] = ap_dec

  ; get list of galaxy IDs in each mag bin for this aperture
  ap_galaxy_maglist = list()
  FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
    ap_galaxies_binned=get_galaxies_binned(ap_ra, ap_dec, ap_radius, datared, mag, mag+binsize)
    ap_galaxy_maglist.Add, ap_galaxies_binned
  ENDFOR 
  ; put the list of ID numbers per mag bin in this aperture into the mega list for all apertures:
  megalist.Add, ap_galaxy_maglist   

  ; get number of galaxies in each bin and error for this aperture 
  FOR i=0, nbins-1. DO BEGIN
    bin = ap_galaxy_maglist[i]
    IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
    ap_ngal_binned[i,ap] = ngal 
    ap_ngalerr_binned[i,ap] = Poisson_error(ngal)
  ENDFOR 

ENDFOR


; calculate densities in each bin in each aperture and errors: 
ap_densities = double(ap_ngal_binned)/(!dPI*ap_radius^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
ap_density_errors = double(ap_ngalerr_binned)/(!dPI*ap_radius^2.*ap_area_weight)   ; same as above except errors 
; average the densities of all the apertures in each bin to get average field density in each bin: 
field_densities = total(ap_densities,2)/double(nApertures)             ; field density in each bin 
field_density_errors = total(ap_density_errors,2)/double(nApertures)   ; errors in field density in each bin 

; check with user whether to make a plot or not since this plot isn't completely necessary
makeplot = ''
READ, makeplot, PROMPT='Make plot of density of galaxies in field? (y/n)'
IF (makeplot EQ 'y') THEN BEGIN 
  ; now, plot up the field density with errors
  window, 2, retain=2, xsize=1200, ysize=1000
  plot, mags, field_densities, title="Galaxies in Field", xtitle="magnitude (F140W)", ytitle="N!Igal!N per deg!E2!N per 0.5 mag", background=255, color=0, charsize=1.5, psym=-8, /ylog, yrange=[1d2,1d6], /ystyle, xrange=magrange, /xstyle 
  errplot, mags, field_densities-field_density_errors, field_densities+field_density_errors, color=0
  write_png, "fielddensity.png", tvrd(/true)
ENDIF 


; set up vertices of polygon for polyfill
ypoints=dblarr(2*nbins)
ypoints[0:nbins-1] = field_densities-field_density_errors
FOR i=0, nbins-1 DO BEGIN
ypoints[i+nbins] = field_densities[nbins-1-i]+field_density_errors[nbins-1-i]
ENDFOR
xpoints = dblarr(2*nbins)
xpoints[0:nbins-1] = mags
FOR i=0, nbins-1 DO BEGIN
 xpoints[i+nbins] = mags[nbins-1-i]
ENDFOR

; make the main plot comparing blob to field with error bars 
window, 3, retain=2, xsize=1200, ysize=1000
plot, mags, densities, title="Galaxy Overdensity in Blob Region", xtitle="magnitude (F140W)", ytitle="N!Igal!N per deg!E2!N per 0.5 mag", background=255, color=0, charsize=2., psym=-8, /ylog, yrange=[1d2,1d6], /ystyle, xrange=magrange, /xstyle 
polyfill, xpoints, ypoints, color=200
errplot, mags, densities-density_errors, densities+density_errors, color=0
oplot, mags, densities, psym=-8, color=0
oplot, mags, field_densities, color=0, linestyle=2
;
;errplot, mags, field_densities-field_density_errors, field_densities+field_density_errors, color=0
LEGEND, ['blob','field'], /right, /top, color=0, textcolor=0, linestyle=[0,2], charsize=2.    
axis, color=0, xaxis=0, xrange=magrange, /xstyle, /data, charsize=0, xtickformat="(A1)"
axis, color=0, yaxis=0, /ylog, yrange=[1d2,1d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)"
;write_png, "overdensity.png", tvrd(/true)
;write_png, "overdensity1.png", tvrd(/true)
write_png, "overdensity2.png", tvrd(/true) 

 


END














FUNCTION get_galaxies, ra, dec, radius, data
  distance = sqrt(  ( (ra - data.ALPHA_J2000) * cos(dec*!dPI/180.) )^2. + (dec - data.DELTA_J2000)^2. )
  galaxies_in_aperture = WHERE (((distance LT radius) AND (data.CLASS_STAR LT 0.8) AND (data.IMAFLAGS_ISO EQ 0.)), n_galaxies)
  RETURN, galaxies_in_aperture
END




FUNCTION get_galaxies_noblob, ra, dec, radius, data
  distance = sqrt(  ( (ra - data.ALPHA_J2000) * cos(dec*!dPI/180.) )^2. + (dec - data.DELTA_J2000)^2. )
  galaxies_outside_aperture = WHERE (((distance GE radius) AND (data.CLASS_STAR LT 0.8) AND (data.IMAFLAGS_ISO EQ 0.)), n_galaxies)
  RETURN, galaxies_outside_aperture
END






FUNCTION get_galaxies_binned, ra, dec, radius, data, brightmag, dimmag
  distance = sqrt(  ( (ra - data.ALPHA_J2000) * cos(dec*!dPI/180.) )^2. + (dec - data.DELTA_J2000)^2. )
  galaxies_in_aperture = WHERE (((distance LT radius) AND (data.CLASS_STAR LT 0.8) AND (data.IMAFLAGS_ISO EQ 0.) AND (data.MAG_ISO GE brightmag) AND (data.MAG_ISO LT dimmag)), n_galaxies)
  RETURN, galaxies_in_aperture
END






FUNCTION Poisson_error, ngal 
  error = sqrt(double(ngal))
  RETURN, error
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

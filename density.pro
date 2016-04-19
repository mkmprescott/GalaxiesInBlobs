PRO density
;--------------------------------------------------------------------------------------------------
;  some day this will include a nice thorough header and description... but not today
;--------------------------------------------------------------------------------------------------
FORWARD_FUNCTION get_galaxies
FORWARD_FUNCTION rsex
FORWARD_FUNCTION get_galaxies_binned
FORWARD_FUNCTION Poisson_error
FORWARD_FUNCTION get_galaxies_noblob
FORWARD_FUNCTION colorsort

; read in the file containing all blob information
blobs = read_csv('blobs.csv', n_table_header=1)
; assign names to each field of the "blobs" structure
blobname = blobs.FIELD01
blobra = blobs.FIELD02
blobdec = blobs.FIELD03
blueband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/catalogs/'+blobs.FIELD04
midband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/catalogs/'+blobs.FIELD05
redband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/catalogs/'+blobs.FIELD06

stack = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/stacks/'+blobs.FIELD07
mask = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/stacks/'+blobs.FIELD08 
z = blobs.FIELD09           ; redshift
m = blobs.FIELD10           ; slope of color-cut line 
b = blobs.FIELD11           ; y-intercept of color-cut line 
cap = blobs.FIELD12         ; mag above which all colors are consistent with redshift
nblobs = n_elements(blobname)        ; number of blobs in the blobs.csv file 
ap_radius_arcsec = 10.               ; aperture with radius of 10 arcsec, in arcsec 
ap_radius = ap_radius_arcsec/3600.   ; aperture with radius of 10 arcsec in decimal degrees
; mag binning stuff:
binsize=0.5
brightestmag = 22.
dimmestmag = 35.
nbins = (dimmestmag - brightestmag + binsize)/binsize
mags = brightestmag + binsize*findgen(nbins)
magrange=[brightestmag,dimmestmag]
; for plots: make the psym=8 symbol a filled circle
plotsym, 0, 1, /fill   
; random aperture stuff:
nApertures = 100.   ; number of random apertures to use in field density calculation
pixelscale = 0.06   ; arceseconds per pixel
ap_radius_pix = ap_radius_arcsec / pixelscale


FOR blob=0, nblobs-2 DO BEGIN           ; eventually it will be nblobs-1 like usual, but PRG3 isn't up and running yet :) 

  ; make structures for each SExtractor catalog (one per each of 3 bands for each blob) 
  datablue = rsex(blueband[blob])    ; either F475W (PRG1 and PRG3) or F606W (PRG2)
  datamid = rsex(midband[blob])      ; F814W
  datared = rsex(redband[blob])      ; F140W

  ; image = mrdfits(stack[blob],0,hdr)
  bpm = mrdfits(mask[blob],0,hdr)    ; this is the mask image for weeding out bad pixels in field 
  sizeinfo = size(bpm)               ; the dimensions of the mask image in pixels 
  npix_x = sizeinfo[1]               ; number of pixels in mask in the x direction
  npix_y = sizeinfo[2]               ; number of pixels in mask in the y direction 



  ; First, a test: get raw blob density, ie, number of galaxies inside aperture 
  ;   centered on blob with radius specified above with no mag bins yet
  blob_galaxies_all=get_galaxies(blobra[blob], blobdec[blob], ap_radius, datared)
  blob_ngalraw = n_elements(blob_galaxies_all)
  ; now get density in number of galaxies per square arcsecond 
  blobdensityraw = float(blob_ngalraw) / (!dPI*ap_radius^2.)    
  print, blobdensityraw, " galaxies per square degree in blob region for "+blobname[blob]  
  print, blob_ngalraw, " galaxies in the blob region for "+blobname[blob] 



  ; make new catalogs with color cuts
  xcolorall = datablue.MAG_ISO - datamid.MAG_ISO
  ycolorall = datamid.MAG_ISO - datared.MAG_ISO
  dataredcut = colorsort(datared, xcolorall, ycolorall, m[blob], b[blob], cap[blob]) 
  datamidcut = colorsort(datamid, xcolorall, ycolorall, m[blob], b[blob], cap[blob]) 
  databluecut = colorsort(datablue, xcolorall, ycolorall, m[blob], b[blob], cap[blob]) 



  ; mag binning: 

  galaxy_maglist = list()
  FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
    blob_galaxies_binned=get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius, datared, mag, mag+binsize)
    galaxy_maglist.Add, blob_galaxies_binned
  ENDFOR      ; Now we have a list containing arrays of the IDs of every galaxy in each mag bin, one array per bin
 
  blob_ngal_binned = dblarr(nbins)
  blob_ngalerr_binned = dblarr(nbins) 

  FOR i=0, nbins-1. DO BEGIN
    bin = galaxy_maglist[i]
    IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
    blob_ngal_binned[i] = ngal 
    blob_ngalerr_binned[i] = Poisson_error(ngal) 
  ENDFOR 
  print, total(blob_ngal_binned), " total galaxies"   ; just a test to make sure it's working - this should match the previous no-bin number 

  ; check with user whether to make a plot or not since this plot isn't completely necessary
  makeplot = ''
  READ, makeplot, PROMPT='Make plot of number of galaxies in blob region? (y/n)'
  IF (makeplot EQ 'y') THEN BEGIN 
    ; Now, plot up number in each bin with errors:
    window, 0, retain=2, xsize=1200, ysize=1000
    plot, mags, blob_ngal_binned, title=("Galaxies in Blob Region for " + blobname[blob]), xtitle="magnitude (F140W)", ytitle="number of galaxies", background=255, color=0, charsize=1.5, psym=10, yrange=[0,15], /ystyle, xrange=magrange, /xstyle 
    errplot, mags, blob_ngal_binned-blob_ngalerr_binned, blob_ngal_binned+blob_ngalerr_binned, color=0
    namestring = string(blobname[blob]) + '_blobngal.png'
    write_png, namestring, tvrd(/true)
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
    plot, mags, densities, title=("Galaxies in Blob Region for " + blobname[blob]), xtitle="magnitude (F140W)", ytitle="N!Igal!N per deg!E2!N per 0.5 mag", background=255, color=0, charsize=1.5, psym=-8, /ylog, yrange=[1d4,1d6], /ystyle, xrange=magrange, /xstyle 
    errplot, mags, densities-density_errors, densities+density_errors, color=0
    namestring = string(blobname[blob]) + '_blobdensity.png'
    write_png, namestring, tvrd(/true)
  ENDIF  



  ; now do the same thing with the color cuts
  galaxy_maglist_cuts = list()
  FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
    blob_galaxies_binned_cuts=get_galaxies_binned(blobra[blob], blobdec[blob], ap_radius, dataredcut, mag, mag+binsize)
    galaxy_maglist_cuts.Add, blob_galaxies_binned_cuts
  ENDFOR      ; Now we have a list containing arrays of the IDs of every galaxy in each mag bin, one array per bin
 
  blob_ngal_binned_cuts = dblarr(nbins)
  blob_ngalerr_binned_cuts = dblarr(nbins) 

  FOR i=0, nbins-1. DO BEGIN
    bin = galaxy_maglist_cuts[i]
    IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
    blob_ngal_binned_cuts[i] = ngal 
    blob_ngalerr_binned_cuts[i] = Poisson_error(ngal) 
  ENDFOR 
  print, total(blob_ngal_binned_cuts), " total galaxies with color cut" 
  ; Finally, compute density and errors: 
  densities_cuts = double(blob_ngal_binned_cuts)/(!dPI*ap_radius^2.) 
  density_errors_cuts = double(blob_ngalerr_binned_cuts)/(!dPI*ap_radius^2.)










  ; color-color plot
  redmags = datared[blob_galaxies_all].MAG_ISO
  midmags = datamid[blob_galaxies_all].MAG_ISO
  bluemags = datablue[blob_galaxies_all].MAG_ISO
  ; make the colors
  color_x = bluemags - midmags
  color_y = midmags - redmags
  ; set up a plot 
  window, 31, retain=2, xsize=1200, ysize=1000
  plot, color_x, color_y, title=(blobname[blob]+" Color-Color Diagram"), xtitle="blue - mid", ytitle="mid - red", background=255, color=0, charsize=1.5, psym=8, xrange=[-5,5], /xstyle, yrange=[-5,5], /ystyle 
  field_galaxies = get_galaxies_noblob(blobra[blob], blobdec[blob], ap_radius, datared)
  field_redmags = datared[field_galaxies].MAG_ISO
  field_midmags = datamid[field_galaxies].MAG_ISO
  field_bluemags = datablue[field_galaxies].MAG_ISO
  field_color_x = field_bluemags - field_midmags
  field_color_y = field_midmags - field_redmags
  oplot, field_color_x, field_color_y, color=0, psym=3
  LEGEND, ['blob galaxy','field galaxy'], /left, /bottom, color=0, textcolor=0, psym=[8,3], charsize=1., /box, outline_color=0 
  xvalues = 0.01*findgen(100000) - 500.
  cut = xvalues*m[blob] + b[blob]
  cut[WHERE(cut GT cap[blob])] = cap[blob]
  oplot, xvalues, cut, color=50, linestyle=2
  saveplot = ''
  READ, saveplot, PROMPT='Save color-color plot? (y/n)'
  IF (saveplot EQ 'y') THEN BEGIN 
    namestring = string(blobname[blob]) + '_color-color_withcut.png'
    write_png, namestring, tvrd(/true)
  ENDIF 



  ; see if user wants to do the overdensity (ie, field sampling) stuff 
  proceed=''
  READ, proceed, PROMPT='Proceed with field density/overdensity analysis? (y/n)'
  IF (proceed EQ 'y') THEN BEGIN 



    ; SAMPLING THE FIELD

    ap_ras = dblarr(nApertures)     ; this will contain RAs of all apertures
    ap_decs = dblarr(nApertures)    ; this will contain declinations of all apertures
    ap_ngal_binned = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
    ap_ngalerr_binned = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
    ap_ngal_binned_cuts = dblarr(nbins,nApertures)  
    ap_ngalerr_binned_cuts = dblarr(nbins,nApertures)   
    ;ap_galaxy_list = list()     ; this will contain the IDs of all galaxies in each aperture 
    megalist = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 
    megalist_cuts = list()


    FOR ap=0, nApertures-1 DO BEGIN

      repeatflag = 0.       ; this will ensure that apertures don't get wasted on areas with only bad pixels 
      WHILE (repeatflag EQ 0.) DO BEGIN
        ; make a random aperture:
        ap_x = double(npix_x)*randomu(seed)
        ap_y = double(npix_y)*randomu(seed)

        ; FIND BAD PIXELS:
        ; make an aperture image for the random aperture generated
        dist_ellipse, dim, [npix_x, npix_y], ap_x, ap_y, 1.,0.
        newap = dim       ; this will be used to overlap with the mask image
        FOR i=0, npix_x-1 DO BEGIN
          FOR j=0, npix_y-1 DO BEGIN
            IF (newap[i,j] LT ap_radius_pix) THEN BEGIN
                newap[i,j] = 1.          ; every pixel in the aperture has a value of 1
            ENDIF ELSE newap[i,j] = 0.   ; pixels outside the aperture are 0
          ENDFOR
        ENDFOR
        ; multiply the aperture map by the bpm so only GOOD pixels INSIDE the aperture are left: 
        apbpm = bpm*newap          
        badpix = where(apbpm GT 0, nbadpix)   ; label the bad pixels
        ntotalpix = n_elements(apbpm)         ; get the total number of pixels in the aperture

        ; FIND VICINITY TO BLOB: 
        xyad, hdr,ap_x, ap_y, ap_ra, ap_dec      ; convert pixels into coordnates in sky
        ; find how close the aperture is to the blob:
        vicinity = sqrt(  ( (blobra[blob] - ap_ra) * cos(blobdec[blob]*!dPI/180.) )^2. + (blobdec[blob] - ap_dec)^2. )

        ; if the aperture has at least some good pixels AND isn't near the blob, then keep it;
        ;    otherwise, throw it out and make a new aperture  
        IF ((nbadpix LT ntotalpix) AND (vicinity GT 2.*ap_radius)) THEN BEGIN 
          repeatflag = 1. 
        ENDIF ELSE repeatflag = 0.
      ENDWHILE
      ; once the aperture passes the above WHILE test, it can be used

      ; fraction of aperture that is good pixels:
      ap_area_weight = (float(ntotalpix) - float(nbadpix)) / float(ntotalpix) 
 
      ; convert pixels to sky coords and add these coords to arrays of aperture properties
      xyad, hdr,ap_x, ap_y, ap_ra, ap_dec
      ap_ras[ap] = ap_ra
      ap_decs[ap] = ap_dec

      ; get list of galaxy IDs in each mag bin for this aperture
      ap_galaxy_maglist = list()
      ap_galaxy_maglist_cuts = list()
      FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
        ap_galaxies_binned=get_galaxies_binned(ap_ra, ap_dec, ap_radius, datared, mag, mag+binsize)
        ap_galaxies_binned_cuts = get_galaxies_binned(ap_ra, ap_dec, ap_radius, dataredcut, mag, mag+binsize)
        ap_galaxy_maglist.Add, ap_galaxies_binned
        ap_galaxy_maglist_cuts.Add, ap_galaxies_binned_cuts
      ENDFOR 
      ; put the list of ID numbers per mag bin in this aperture into the mega list for all apertures:
      megalist.Add, ap_galaxy_maglist   
      megalist_cuts.Add, ap_galaxy_maglist_cuts

      ; get number of galaxies in each bin and error for this aperture 
      FOR i=0, nbins-1. DO BEGIN
        bin = ap_galaxy_maglist[i]
        IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
        ap_ngal_binned[i,ap] = ngal 
        ap_ngalerr_binned[i,ap] = Poisson_error(ngal)
      ENDFOR 

      FOR i=0, nbins-1. DO BEGIN
        bin = ap_galaxy_maglist_cuts[i]
        IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
        ap_ngal_binned_cuts[i,ap] = ngal 
        ap_ngalerr_binned_cuts[i,ap] = Poisson_error(ngal)
      ENDFOR 

    ENDFOR    ; all the random apertures


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
      plot, mags, field_densities, title=("Galaxies in Field for " + blobname[blob]), xtitle="magnitude (F140W)", ytitle="N!Igal!N per deg!E2!N per 0.5 mag", background=255, color=0, charsize=1.5, psym=-8, /ylog, yrange=[1d2,1d6], /ystyle, xrange=magrange, /xstyle 
      errplot, mags, field_densities-field_density_errors, field_densities+field_density_errors, color=0
      namestring = string(blobname[blob]) + '_fielddensity.png'
      write_png, namestring, tvrd(/true)
    ENDIF 


    ; calculate densities in each bin in each aperture and errors WITH COLOR CUTS: 
    ap_densities_cuts = double(ap_ngal_binned_cuts)/(!dPI*ap_radius^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
    ap_density_errors_cuts = double(ap_ngalerr_binned_cuts)/(!dPI*ap_radius^2.*ap_area_weight)   ; same as above except errors 
    ; average the densities of all the apertures in each bin to get average field density in each bin: 
    field_densities_cuts = total(ap_densities_cuts,2)/double(nApertures)             ; field density in each bin 
    field_density_errors_cuts = total(ap_density_errors_cuts,2)/double(nApertures)   ; errors in field density in each bin 


    ; MAIN OVERDENSITY PLOT:

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
    title=("Galaxy Overdensity in Blob Region for " + blobname[blob])
    plot, mags, densities, xtitle="magnitude (F140W)", ytitle="N!Igal!N per deg!E2!N per 0.5 mag", title=title, xthick=4,ythick=4,background=255, color=0, charsize=4., psym=-8, /ylog, yrange=[1d3,1d6], /ystyle, xrange=[22.5,29.5], /xstyle, charthick=4 
    polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,1d6], /data, noclip=0
    errplot, mags, densities-density_errors, densities+density_errors, color=0, thick=4
    oplot, mags, densities, psym=-8, color=0, thick=4     ; have to do this again because polyfill covers it up 
    oplot, mags, field_densities, color=0, linestyle=2, thick=4
    LEGEND, ['blob','field'], /center, /bottom, color=0, textcolor=0, linestyle=[0,2], thick=4., charsize=3., charthick=3. , box=0, number=0.1   
    axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=4, charthick=4, xtickformat="(A1)", xthick=4
    axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=4
    axis, color=0, yaxis=0, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=4, charthick=4, ytickformat="(A1)", ythick=4
    axis, color=0, yaxis=1, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=4
    saveplot = ''
    READ, saveplot, PROMPT='Save overdensity plot? (y/n)'
    IF (saveplot EQ 'y') THEN BEGIN 
       namestring = string(blobname[blob]) + '_overdensity.png'
       write_png, namestring, tvrd(/true) 
    ENDIF 





 ; MAIN OVERDENSITY PLOT WITH COLOR CUTS:

    ; set up vertices of polygon for polyfill
    ypoints=dblarr(2*nbins)
    ypoints[0:nbins-1] = field_densities_cuts-field_density_errors_cuts
    FOR i=0, nbins-1 DO BEGIN
    ypoints[i+nbins] = field_densities_cuts[nbins-1-i]+field_density_errors_cuts[nbins-1-i]
    ENDFOR
    xpoints = dblarr(2*nbins)
    xpoints[0:nbins-1] = mags
    FOR i=0, nbins-1 DO BEGIN
     xpoints[i+nbins] = mags[nbins-1-i]
    ENDFOR

    ; make the main plot comparing blob to field with error bars 
    window, 4, retain=2, xsize=1200, ysize=1000
    title=("Galaxy Overdensity in Blob Region for " + blobname[blob]) + ' (with color cuts)'
    plot, mags, densities_cuts, xtitle="magnitude (F140W)", ytitle="N!Igal!N per deg!E2!N per 0.5 mag", title=title,xthick=2,ythick=2,background=255, color=0, charsize=2., psym=-8, /ylog, yrange=[1d3,1d6], /ystyle, xrange=[22.5,29.5], /xstyle ;, charthick=4 
    polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,1d6], /data, noclip=0
    errplot, mags, densities_cuts-density_errors_cuts, densities_cuts+density_errors_cuts, color=0, thick=4
    oplot, mags, densities_cuts, psym=-8, color=0, thick=2     ; have to do this again because polyfill covers it up 
    oplot, mags, field_densities_cuts, color=0, linestyle=2, thick=2
    LEGEND, ['blob','field'], /center, /bottom, color=0, textcolor=0, linestyle=[0,2], thick=2., charsize=1.5, box=0, number=0.1   
    axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=2, charthick=2, xtickformat="(A1)", xthick=2
    axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=2
    axis, color=0, yaxis=0, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=2, charthick=2, ytickformat="(A1)", ythick=2
    axis, color=0, yaxis=1, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=2
    saveplot = ''
    READ, saveplot, PROMPT='Save overdensity plot? (y/n)'
    IF (saveplot EQ 'y') THEN BEGIN 
       namestring = string(blobname[blob]) + '_overdensity_withcuts.png'
       write_png, namestring, tvrd(/true) 
    ENDIF 



    ; combine everything to make a statistical luminosity function
    xvertices = fltarr(2.0*n_elements(mags) + 2.)
    yvertices = fltarr(2.0*n_elements(mags) + 2.)
    xvertices[0] = brightestmag
    xvertices[-1] = dimmestmag
    yvertices[0] = 0
    yvertices[-1] = 0
    xvertices[1] = brightestmag
    xvertices[-2] = dimmestmag
    FOR i=1, n_elements(mags) -1 DO BEGIN
      xvertices[2*i] = mags[i]-(binsize/2.0)
      xvertices[2*i+1] = mags[i]-(binsize/2.0)
    ENDFOR
    FOR i=1, n_elements(mags) -1 DO BEGIN
      yvertices[2*i+1] = blob_ngal_binned_cuts[i]
      yvertices[2*i+2] = blob_ngal_binned_cuts[i]
    ENDFOR
    ; now make statistical luminosity function
    overdensity_cuts = densities_cuts - field_densities_cuts 
    n_overdensity_cuts = overdensity_cuts * (!dPI*ap_radius^2.)  ; convert to number
    window, 5, retain=2, xsize=1200, ysize=1000
    title=("Statistical Luminosity Function for Galaxies in " + blobname[blob]) 
    plot, mags, n_overdensity_cuts, background=255, color=0, charsize=2, linestyle=2, yrange=[0,10], /ystyle,xrange=[22,30],/xstyle, title=title, xtitle='magnitude (F140W)',ytitle='N!Igal!N'
    polyfill, xvertices,yvertices, color=220, clip=[22,0,30,10], /data, noclip=0
    oplot, mags, blob_ngal_binned, color=0, psym=10
    oplot, mags, blob_ngal_binned_cuts, color=0, psym=10
    oplot, mags, n_overdensity_cuts, color=0, thick=2, linestyle=2
    saveplot = ''
    READ, saveplot, PROMPT='Save statistical luminosity function plot? (y/n)'
    IF (saveplot EQ 'y') THEN BEGIN 
       namestring = string(blobname[blob]) + '_statlum.png'
       write_png, namestring, tvrd(/true) 
    ENDIF 






  ENDIF   ; if the user wanted to do the field/overdensity stuff
 

ENDFOR ;all the blobs



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






FUNCTION colorsort, data, xcolor, ycolor, m, b, c 
  good = data[WHERE ((ycolor GE (xcolor*m + b)) OR (ycolor GE c))]
  RETURN, good
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
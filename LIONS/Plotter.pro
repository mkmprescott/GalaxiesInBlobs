;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; PLOTTER.PRO 
;  Contains all modules used for making plots regarding the galaxy number density in blobs and in fields. 
;  CONTENTS: 
;      densityplot------------------------procedure
;      statlum----------------------------procedure
;      makehist---------------------------procedure
;      zhist------------------------------procedure
;      Rmultirad--------------------------procedure
;      Rmultimag--------------------------procedure
;      makereg----------------------------procedure  ???????????
;      plotsave---------------------------procedure 
;      titlemaker-------------------------function  ?????????????
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;




PRO densityplot, file, magband, windownumber, indicestoplot=indicestoplot, autosave=autosave
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------; 
; densityplot procedure
;   Makes plots of galaxy number density vs magnitude for a user-specified blob and overplots the galaxy density in the field for comparison.
;
; INPUTS: file - a string giving the name of a FITS file containing density data for blobs and fields (must NOT include .FITS extension)
;         magband - a string giving the magnitude band being used to bin the galaxies, for labeling the plot's x-axis 
;         windownumber - an integer specifying a unit number for the plot window (in case user wants to run this code multiple times and see all plots side-by-side) 
;         indicestoplot - an optional 2-element integer array specifying which indices of the data structure to use for the blob and field; the first element should be the 
;                          index of the data structure to use for the blob, and the second should be the index of the data structure to use for the field; if this keyword
;                          isn't set, the code will ask the user to select a blob and field  
;         autosave - an optional string keyword specifying whether or not to automatically save the density plot; must be set to 'y' to save the plot automatically; the plot will
;                          not be saved automatically if this keyword is not set or is set to except 'y' 
;
; OUTPUT: a plot of galaxy number density vs magnitude with error bars for the blob superimposed on a gray swath showing the 1-sigma galaxy density range for the field; 
;          this can be saved as a PNG by the user under the same name as the FITS file specified by "file" (with the name of the blob added on)  
;
; Uses the plotsave procedure. 
;
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------; 

  ; SETUP 
     plotsym, 0, 1, /fill            ; plotting symbol: make the psym=8 symbol a filled circle for our plots 
     loadct, 0                       ; load the grayscale color table   
     !P.MULTI = [0,1,1,0,1]          ; make sure there's only one plot per window 

  ; READ IN DATA 
     datafile = '/boomerang-data/alhall/GalaxiesInBlobs/LIONScatalogs/densities/' + file + '.fits'       ; add the full path to the name of the file given by the user  
     data = mrdfits(datafile, 1)                                                 ; read in the data 
     mags = data(0).mags                                                         ; these are the magnitude bins into which the data are sorted
     nbins = n_elements(mags)                                                    ; determine how many magnitude bins there are in total  

  ; DETERMINE WHICH BLOB AND FIELD TO PLOT VIA USER INPUT
     IF (n_elements(indicestoplot) EQ 2) THEN BEGIN           ; if indicestoplot is given and has the right number of elements
       blob = indicestoplot[0]                                          ; this is the data index to use for the blob
       field = indicestoplot[1]                                        ; this is the data index to use for the field 
     ENDIF ELSE BEGIN                                         ; if indicestoplot is either not specified or has the wrong number of elements, ask user for a blob and field
       print, 'OBJECT LIST: '                   ; we'll list everything in data structure so user can see what's there
       lastdataindex = n_elements(data)-1       ; find the last index of the data structure (the total number of blobs and fields in the file then -1 since array indices start at 0) 
       FOR i=0, lastdataindex DO BEGIN          ; loop over the data structure, printing the names of all objects (blobs and fields) contained in it 
         print, fix(i), '     ', data(i).name        ; assign a number to each object in the data structure 
       ENDFOR                                   ; now there's an object list 
       ; BLOB: get blob from user 
          blob = 0S                                                                                               ; integer, will be "data" index of the blob the user wants to plot 
          read, blob, prompt='Choose a blob to plot from the object list shown: '                                            ; user will enter the index of a blob in "data"
          typecode = size(blob, /type)                                                                                       ; a type code of 2 is an integer, which is what's needed 
          IF (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((blob LT 0) OR (blob GT lastdataindex)) )   ) THEN BEGIN       ; make sure user entered a valid integer
            WHILE (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((blob LT 0) OR (blob GT lastdataindex)) )   ) DO BEGIN    ; while user response is invalid
              print, 'Invalid response. Select from the listed object numbers. Enter response as an integer.'                     ; tell user response was invalid
              blob = 0S                                                                                                           ; reset "blob" to be an integer 
              read, blob, prompt='Choose a blob to plot from the object list shown: '                                             ; get user input again 
              typecode = size(blob, /type)                                                                                        ; check input data type 
            ENDWHILE                                                                                                         ; while the "blob" response is invalid 
          ENDIF ELSE IF ((typecode EQ 2) AND (blob GE 0) AND (blob LE lastdataindex)) THEN BEGIN                             ; if the user entered input correctly 
            print, 'selected blob: ', data(blob).name                                                                             ; show user which object they selected 
          ENDIF                                                                                                              ; the "blob" response is valid 
       ; FIELD: get field from user 
          field = 0S                                                                                              ; integer, will be "data" index of the field the user wants to plot
          read, field, prompt='Choose a field to overplot from the object list shown: '                                      ; user will enter the index of a field in "data"
          typecode = size(field, /type)                                                                                      ; a type code of 2 is an integer, which is what's needed 
          IF (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((field LT 0) OR (field GT lastdataindex)) )   ) THEN BEGIN     ; make sure user entered a valid integer
            WHILE (  (typecode NE 2)   OR   ( (typecode EQ 2) AND ((field LT 0) OR (field GT lastdataindex)) )   ) DO BEGIN  ; while user response is invalid
              print, 'Invalid response. Select from the listed object numbers. Enter response as an integer.'                     ; tell user response was invalid
              field= 0S                                                                                                           ; reset "field" to be an integer 
              read, field, prompt='Choose a field to overplot from the object list shown: '                                       ; get user input again 
              typecode = size(field, /type)                                                                                       ; check input data type 
            ENDWHILE                                                                                                         ; while the "field" response is invalid 
          ENDIF ELSE IF ((typecode EQ 2) AND (field GE 0) AND (field LE lastdataindex)) THEN BEGIN                           ; if the user entered input correctly 
            print, 'selected field: ', data(field).name                                                                           ; show user which object they selected 
          ENDIF                                                                                                              ; the "field" response is valid 
     ENDELSE        ; if the user didn't give the indices to plot for the blob and field as a keyword 
     ; use the "blob" and "field" indices to get the information we want to plot 
        Bdensity = data(blob).densities       ; blob density in each mag bin 
        Berr = data(blob).derrs               ; uncertainty on above 
        Blabel = data(blob).name              ; name of blob to label on plot 
        Fdensity = data(field).densities      ; field density in each mag bin 
        Ferr = data(field).derrs              ; uncertainty on above
        Flabel = data(field).name             ; name of field to label on plot 
             

  ; SET UP FIELD SWATH: vertices of polygon for polyfill
     ; Perimeter of field swath will be drawn starting at lower left (small x, small y) and moving to lower right (big x, small y), then to upper right (big x, big y), 
     ;   then to upper left (small x, big y) and connect back at the original starting point from there. That means that the ordered pairs showing where each point along 
     ;   the perimeter lies must follow this order, so x (independent coordinate) will start small and increase to its max value to form the bottom of the swath while 
     ;   y (dependent) stays constantly small (ie, first half of y coordinates will be lower error bars), and then after the max x value has been reached along the 
     ;   bottom swath edge (ie, when the perimeter reaches the lower right corner), y jumps up to large values (upper error bars) and x decreases back to its minimum
     ;   to create the top edge of the swath. 
        ; start with x coordinates of vertices - x is magnitude
        xpoints = dblarr(2*nbins)                                  ; there are two points for each magnitude because there's a lower error bar and an upper error bar
        xpoints[0:nbins-1] = mags                                  ; fill in first half of xpoints array with increasing x values 
        FOR i=0, nbins-1 DO BEGIN                                  ; fill in second half of xpoints array 
          xpoints[i+nbins] = mags[nbins-1-i]                       ; decreasing x values to track back to where we started 
        ENDFOR                                                     ; now all x coordinates of the swath's vertices are filled in 
        ; now the y coordinates of vertices - y is the actual data (number densities) 
        ypoints=dblarr(2*nbins)                                    ; same number of y points as x points since these are ordered pairs
        ypoints[0:nbins-1] = Fdensity-Ferr                         ; fill in first half of ypoints array with lower bounds of error bars 
        FOR i=0, nbins-1 DO BEGIN                                  ; fill in second half of ypoints array     
          ypoints[i+nbins] = Fdensity[nbins-1-i]+Ferr[nbins-1-i]   ; upper bounds of error bars, running backward along x-axis to match xpoints array
        ENDFOR                                                     ; now all y coordinates of the swath's vertices are filled in  
      
                                                                                        
  ; MAKE THE PLOT     
    ; set up titles for plot and axes                                                                               
       xtitle = 'magnitude (' + magband + ')'                                                                                                                               
       ytitle='N!Igal!N per deg!E2!N per 0.5 mag'   
       nameofblob = strcompress(data(blob).name, /remove)        ; this will be used in the PNG filename as well if the user saves this plot 
       plottitle = nameofblob + ' Density Analysis'
    ; open a window for the plot 
       window, windownumber, retain=2, xsize=700, ysize=550
    ; plot the data                                                                                                                      
       ;     ; data          ; plot colors and symbols         ; text labels with large font                                                ; wide margins                      
       plot, mags, Bdensity, background=255, color=0, psym=-8, title=plottitle, xtitle=xtitle, ytitle=ytitle, charsize = 2.5, charthick=2., xmargin=[8,3], ymargin=[4,2.5], $ 
             xrange=[22.5,29.5], /xstyle, yrange=[1d3,3d6], /ystyle, /ylog               ; how axes are numbered                
       polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,3d6], /data, noclip=0     ; overplot gray swath showing field density uncertainty 
       oplot, mags, Fdensity, color=0, linestyle=2, thick=3                                 ; overplot dashed line showing field density      
       errplot, mags, Bdensity-Berr, Bdensity+Berr, color=0, thick=3                        ; overplot error bars on the blob points
       oplot, mags, Bdensity, psym=-8, color=0, thick=3                                     ; have to plot blob density again because polyfill covers it up   
    ; make a legend showing what represents the blob and what represents the field 
       ;       ; labels         ; properties of legend lines                             ; properties of legend text              ; location and properties of legend box  
       LEGEND, [Blabel,Flabel], color=[0,0], linestyle=[0,2], thick=[3.,3.], number=0.1, textcolor=0, charsize=1.5, charthick=2., /left, /top, /box, outline_color=0. 
    ; redraw the axes to ensure they're visible 
       axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, xtickformat="(A1)", xthick=2, charsize=4, charthick=2.                    ; bottom x axis (has labels) 
       axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, xtickformat="(A1)", xthick=2, charsize=0, charthick=0                     ; top x axis (no labels) 
       axis, color=0, yaxis=0, yrange=[1d3,3d6],   /ystyle, /data, ytickformat="(A1)", ythick=2, charsize=4, charthick=2., /ylog, /ynozero   ; left y axis (has labels) 
       axis, color=0, yaxis=1, yrange=[1d3,3d6],   /ystyle, /data, ytickformat="(A1)", ythick=2, charsize=0, charthick=0,  /ylog, /ynozero   ; right y axis (no labels) 

  ; SAVE THE PLOT 
     savefile = '/boomerang-data/alhall/GalaxiesInBlobs/' + nameofblob + '_' + file + '.png'   ; add the full path onto the filename where the plot will be saved 
     IF (n_elements(autosave) NE 0) THEN BEGIN                            ; check if user has a preference on whether or not to save the plot automatically
       IF (autosave EQ 'y') THEN BEGIN                                       ; autosave has to be set to 'y' in order to save the plot automatically
         write_png, savefile, tvrd(/true)                                       ; write the plot to the PNG file specified by namestring 
         print, 'Plot saved automatically.'                                     ; notify the user that the plot was saved
       ENDIF ELSE BEGIN                                                      ; if autosave is set to something other than 'y',  
         plotsave, savefile                                                     ; ask user whether or not to save plot
       ENDELSE                                                               ; that takes care of anything autosave could've been set to
     ENDIF ELSE BEGIN                                                     ; if the autosave keyword isn't set
       plotsave, savefile                                                    ; ask user whether or not to save plot
     ENDELSE                                                              ; now the plot will have been saved if so desired 

END          ; of densityplot procedure 









PRO makehist, blobsfile, fieldapsfile,              magband, windownumber, bin=bin, indicestoplot=indicestoplot, autosave=autosave             ;ap_totals, blob_ngal,
;-------------------------------------------------------------------------------------------------------------------------------------------------;
; makehist procedure                                                                
;   description                                                                            
; INPUTS: etc                                                                                  
;         etc                                                                     
; NOTES:                                                                            
;-------------------------------------------------------------------------------------------------------------------------------------------------;
  FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                      
  ; set up titles for plot and axes                                                                    
    title = titlemaker('aphist', blobname, nApertures=nApertures, radius=radius, cuts=cuts)         
    xtitle = 'n galaxies in aperture'                                                 
    ytitle='number of apertures'                                                      
  plothist, ap_totals, background=255, color=0, axiscolor=0, title=title, xtitle=xtitle, ytitle=ytitle, bin=1, xrange=[0,blob_ngal+10.],/xstyle, charsize=1.5, $
    thick=2, ymargin=[4,4]                                                             
  oplot, [blob_ngal, blob_ngal], [0., nApertures], linestyle=2, thick=2, color=0.      
  LEGEND, ['blob','field'], /center, /top, color=0, textcolor=0, linestyle=[0,2], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5
; xyouts, blob_ngal+1., nApertures/10., blobgaltext, color=0, charsize=1.5                            
END   









PRO plotsave, namestring                                                
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------; 
; plotsave procedure                                                      
;   Asks user whether or not to save a plot and saves said plot if the user chooses to do so.     
;
; INPUT: namestring - a string giving the name of the PNG file containing the plot if the user chooses to save it (.PNG extension must be included) 
;
; OUTPUT: either writes the current plot window to a PNG or informs the user that the plot won't be saved, depending on the user's response 
;
; Uses no other procedures nor functions. 
;
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;       
  saveplot = ''                                                     ; create a string that will contain user's yes-or-no selection                                                              
  read, saveplot, prompt='Save plot? (y/n)'                         ; ask user whether or not to save the plot and read in user's response    
  IF ((saveplot NE 'y') AND (saveplot NE 'n')) THEN BEGIN           ; make sure user actually answered 'y' or 'n'      
    WHILE ((saveplot NE 'y') AND (saveplot NE 'n')) DO BEGIN             ; as long as the answer is invalid, keep asking            
      print, 'Invalid response. Choose y or n.'                          ; tell the user their response didn't work            
      read, saveplot, prompt='Save plot? (y/n)'                          ; ask for input again and read in user's response 
    ENDWHILE                                                        ; loop will end once user has selected either 'y' or 'n' 
  ENDIF ELSE IF (saveplot EQ 'y') THEN BEGIN                        ; if user wants to save the plot,                            
    write_png, namestring, tvrd(/true)                                   ; write the plot to the PNG file specified by namestring 
  ENDIF ELSE IF (saveplot EQ 'n') THEN BEGIN                        ; if user doesn't want to save the plot,                 
      print, 'Plot not saved.'                                           ; inform the user that the plot wasn't saved              
  ENDIF                                                             ; all possible user responses have been accounted for                                                             
END           ; of plotsave procedure   


















                                                                             








;   
PRO statlum, mags, binsize, blob_ngal_binned, blob_ngal_binned_cuts, densities, field_densities, ap_radius, blobname, nApertures, ap_radius_arcsec, band
;-------------------------------------------------------------------------------------------------------------------------------------------------;
; statlum procedure                                                                      
;   description                                                                       
; INPUTS: etc                                                                        
;         etc                                                                          
; NOTES:                                                                                   
;-------------------------------------------------------------------------------------------------------------------------------------------------;
  FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                      
  ; define the statistical luminosity function                                                         
    overdensity = densities - field_densities                                         
    n_overdensity = overdensity * (!dPI*ap_radius^2.)  ; convert to number                            
  ; set up vertices of shaded area                                                                     
    xvertices = fltarr(2.0*n_elements(mags) + 2.)                                      
    yvertices = fltarr(2.0*n_elements(mags) + 2.)                                      
    xvertices[0] = mags[0]                                                             
    xvertices[-1] = mags[-1]                                                           
    yvertices[0] = 0                                                                   
    yvertices[-1] = 0                                                                  
    xvertices[1] = mags[0]                                                             
    xvertices[-2] = mags[-1]                                                           
    FOR i=1, n_elements(mags) -1 DO BEGIN                                              
      xvertices[2*i] = mags[i]-(binsize/2.0)                                           
      xvertices[2*i+1] = mags[i]-(binsize/2.0)                                         
    ENDFOR                                                                             
    FOR i=1, n_elements(mags) -1 DO BEGIN                                              
      yvertices[2*i+1] = blob_ngal_binned_cuts[i]                                      
      yvertices[2*i+2] = blob_ngal_binned_cuts[i]                                      
    ENDFOR                                                                             
  ; now plot statistical luminosity function                                                           
    ; set up titles for plot and axes                                                                  
      title = titlemaker('statlum', blobname, nApertures=nApertures, radius=ap_radius_arcsec, band=band)        
      xtitle = 'magnitude (' + band + ')'                                             
      ytitle = 'N!Igal!N'                                                             
    plot, mags, n_overdensity, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, linestyle=2, yrange=[0,10],/ystyle, xrange=[22,30],/xstyle, $
      charsize=3, thick=2, xmargin=[10,7], ymargin=[4,4]                               
    polyfill, xvertices,yvertices, color=220, clip=[22,0,30,10], /data, noclip=0       
    oplot, mags, blob_ngal_binned, color=0, psym=10, thick=2                           
    oplot, mags, blob_ngal_binned_cuts, color=0, psym=10, thick=2                      
    oplot, mags, n_overdensity, color=0, thick=2, linestyle=2                          
END                                                                                    
                                                                                









PRO Rmultirad, ap_radius_arcsec, mags, radiusresults, index, blobname, nApertures, band, cuts       
;-------------------------------------------------------------------------------------------------------------------------------------------------;
; Rmultirad procedure                                                                     
;   description                                                                                 
; INPUTS: etc                                                                        
;         etc                                                                       
; NOTES:                                                                                   
;-------------------------------------------------------------------------------------------------------------------------------------------------; 
  FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                      
  ; set up an array of colors to plot each radius in                                                   
    n_radii = n_elements(ap_radius_arcsec)                                             
    indices = findgen(n_radii)                                                         
    colors = indices*254./n_radii                                                      
  ; now make the plot                                                                                  
    ; set up titles for plot and axes                                                                  
      title = titlemaker('Rmultirad', blobname, nApertures=nApertures, band=band, cuts=cuts)        
      xtitle = 'magnitude (' + band + ')'                                             
      ytitle = 'overdensity factor'                                                   
    plot, mags, radiusresults[0,*,index,0], background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, yrange=[0,50],/ystyle, charsize=2, thick=2, $ 
      ymargin=[4,4]                                                                    
    errplot, mags, radiusresults[0,*,index,0]-radiusresults[0,*,index,1], radiusresults[0,*,index,0]+radiusresults[0,*,index,1], color=0, thick=2 
    FOR r=1, n_radii-1. DO BEGIN                                                       
      oplot, mags, radiusresults[r,*,index,0], color = colors[r], thick=2              
      errplot, mags, radiusresults[r,*,index,0]-radiusresults[r,*,index,1], radiusresults[r,*,index,0]+radiusresults[r,*,index,1], color=colors[r], thick=2
    ENDFOR                                                                             
    LEGEND, strcompress(string(fix(ap_radius_arcsec)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., $ 
      number=0.5, charthick=1.5, thick=2                                               
END                                                                                    









PRO Rmultimag, mags, ap_radius_arcsec, radiusresults, index, blobname, nApertures, band, cuts       
;-------------------------------------------------------------------------------------------------------------------------------------------------; 
; Rmultimag procedure      
;   description      
; INPUTS: etc      
;         etc               
; NOTES:                 
;-------------------------------------------------------------------------------------------------------------------------------------------------; 
  FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                      
  ; set up an array of colors to plot each magnitude bin in                                            
    indices = findgen(n_elements(mags))                                                
    colors = indices*254./float(n_elements(mags))                                      
  ; now make the plot                                                                                  
    ; set up titles for plot and axes                                                                  
      title = titlemaker('Rmultimag', blobname, nApertures=nApertures, band=band, cuts=cuts)        
      xtitle='aperture radius (arcseconds)'                                           
      ytitle='overdensity factor'                                                     
    ; set up an empty plot (needs to have no data in case the first thing plotted is NaN):             
    plot, ap_radius_arcsec, radiusresults[*,1,index,0], background=255, color=colors[1], title=title, xtitle=xtitle, ytitle=ytitle, xrange=[0,30],/xstyle, $  
      yrange=[0,100],/ystyle, charsize=2, ymargin=[4,4], /nodata                       
    ; now add data on top of the empty plot:                                                           
    FOR i=0, n_elements(mags)-1 DO BEGIN     ; need to make sure that NaNs are thrown out (because where field density is 0, we divided by 0)   
      y = radiusresults[*,i,index,0]                                                   
      yerr = radiusresults[*,i,index,1]                                                
      keep = WHERE(finite(y) EQ 1, nkeep)                                              
      IF (nkeep GT 0) THEN BEGIN                                                       
	oplot, ap_radius_arcsec(keep), y(keep), color=colors[i], thick=2               
	errplot, ap_radius_arcsec(keep), y(keep)-yerr(keep), y(keep)+yerr(keep), color=colors[i], thick=2        
      ENDIF                                                                            
    ENDFOR                                                                             
    ; add a legend:                                                                                    
    LEGEND, strcompress(string(fix(mags)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., number=0.5, $  
      charthick=1.5, thick=2                                                          
END                                                                                   









PRO makereg, data, subset, blobname, radius, cuts                                     
;-------------------------------------------------------------------------------------------------------------------------------------------------; 
; makereg procedure                                                       
;   description                                                                   
; INPUTS: etc                                                                       
;         etc                                                                   
; NOTE: Uses titlemaker function.                                       
;-------------------------------------------------------------------------------------------------------------------------------------------------; 
  FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                      
  savethis = ''                                                                        
  READ, savethis, PROMPT='Make region file? (y/n)'                                     
  IF (savethis EQ 'y') THEN BEGIN                                                      
    cat = data[subset]                                                                 
    filename = titlemaker('blobgals', blobname, radius=radius, cuts=cuts, filetype='.reg')           
    openw, 1, filename                                                                 
    FOR i=0, n_elements(cat)-1 DO BEGIN                                                
      printf, 1,'J2000; circle ', cat[i].alpha_j2000, cat[i].delta_j2000, ' 10p '    ;#text={', cat[i].number,'}' 
    ENDFOR                                                                             
    close, 1                                                                           
  ENDIF ELSE IF ((savethis NE 'y') AND (savethis NE 'n')) THEN BEGIN                   
    WHILE ((savethis NE 'y') AND (savethis NE 'n')) DO BEGIN                           
      print, 'Invalid response. Choose y or n.'                                        
      READ, savethis, PROMPT='Make region file? (y/n)'                                 
    ENDWHILE                                                                           
  ENDIF ELSE IF (savethis EQ 'n') THEN BEGIN                                           
      print, 'Region file not created.'                                               
  ENDIF                                                                                
END                            







PRO cdfsblobregmaker
  RAs = [53.1504167, 53.3358333, 52.9670833, 53.0245833, 52.9529167, 53.0816667, 53.0133333, 52.9320833, 53.1854167, 53.1558333, 53.1800000, 53.0279167, 53.1362500, 53.1341667, 52.9108333, 53.2741667]
  decs = [-28.0151389, -27.6863333, -27.9151667, -27.6213056, -27.8731389, -27.7918889, -27.7569444, -27.6424722, -27.7195000, -28.0349167, -27.7161944, -27.7486667, -27.6517500, -27.6908889, -27.7186111, -27.9611111]
  number = indgen(16)+1
  stringnumber = strcompress(string(number), /remove)
  openw, 2, 'cdfsblobs.reg'
  FOR i=0, n_elements(RAs)-1 DO BEGIN                                                
    printf, 2,'J2000; circle ', RAs[i], decs[i], ' 10p ', '#text={', stringnumber[i], '}'
  ENDFOR           
END                                                                         
                                                                            

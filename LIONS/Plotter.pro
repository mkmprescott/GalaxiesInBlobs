;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; PLOTTER.PRO 
;  Contains all modules used for making plots regarding the galaxy number density in blobs and in fields. 
;  CONTENTS: 
;      densityplot------------------------procedure
;      statlum----------------------------procedure
;      makehist---------------------------procedure
;      zhist------------------------------procedure
;      Rmultirad--------------------------procedure
;      Rmultimag--------------------------procedure
;      plotsave---------------------------procedure     ?????
;      titlemaker-------------------------function
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;

; NOTE: MAKEREG NEEDS TO GO SOMEWHERE!!!!!!!!!!!!!!!
      ; plotting symbol: make the psym=8 symbol a filled circle for our plots 
         plotsym, 0, 1, /fill   

PRO testplot
test = mrdfits('/boomerang-data/alhall/GalaxiesInBlobs/redDENSITIES_10arcsec_10000aps.fits', 1)
;plot, test(0).mags, test(0).densities, background=255, color=0, title=plottitle, xtitle='magnitude', ytitle='density',/ylog, yrange=[1d3,5d6],/ystyle, $ 
;         xrange=[22.5,29.5],/xstyle, charsize=2., xthick=2, ythick=2, xmargin=[10,7], ymargin=[6,6], linestyle=0    
;plot, test(0).mags, test(0).densities, linestyle=0        ; PRG1 is solid 
;oplot, test(1).mags, test(1).densities, linestyle=5, color=0       ; PRG2 is long dashes
;oplot, test(2).mags, test(2).densities, linestyle=2, color=0       ; PRG3 is dashes
;oplot, test(3).mags, test(3).densities, linestyle=1, color=0       ; GOODS-S is dotted

nbins = n_elements(test(0).mags)
!P.MULTI = [0,1,3,0,1]
window, 1, retain=2, xsize=500, ysize=1200
densityplot, nbins, test(0).mags, test(0).densities, test(0).derrs, test(3).densities, test(3).derrs, 'PRG1 overdensity test' 
print, test(0).name

;window, 2, retain=2, xsize=1200, ysize=1000
densityplot, nbins, test(1).mags, test(1).densities, test(1).derrs, test(3).densities, test(3).derrs, 'PRG2 overdensity test' 
print, test(1).name

;window, 3, retain=2, xsize=1200, ysize=1000
densityplot, nbins, test(2).mags, test(2).densities, test(2).derrs, test(3).densities, test(3).derrs, 'PRG3 overdensity test' 
print, test(2).name 

print, 'field is ', test(3).name


END




  PRO densityplot, nbins, mags, densities, density_errs, field_densities, field_density_errs, plottitle   
  ;-------------------------------------------------------------------------------------------------------------------------------------------------; 
      ; plotting symbol: make the psym=8 symbol a filled circle for our plots 
         plotsym, 0, 1, /fill   
 ;  FORWARD_FUNCTION titlemaker           ; this function is used in this procedure  
     loadct, 0                                                                                                                                                            
   ; set up vertices of polygon for polyfill                      
     xpoints = dblarr(2*nbins)                                                                                                                                            
     xpoints[0:nbins-1] = mags                                                                                                                                            
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                            
       xpoints[i+nbins] = mags[nbins-1-i]                                                                                                                                 
     ENDFOR                                                                                                                                                               
     ; y is specific to each catalog (band)                       
     ypoints=dblarr(2*nbins)                                                                                                                                              
     ypoints[0:nbins-1] = field_densities-field_density_errs                                                                                                              
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                            
       ypoints[i+nbins] = field_densities[nbins-1-i]+field_density_errs[nbins-1-i]                                                                                        
     ENDFOR                                                                                                                                                               
   ; make the main plot comparing blob to field with error bars     
     ; set up titles for plot and axes                                                                               
        xtitle = 'magnitude (F140W)'                                                                                                                               
        ytitle='N!Igal!N per deg!E2!N per 0.5 mag'                                                                                                                        
     plot, mags, densities, background=255, color=0, title=plottitle, xtitle=xtitle, ytitle=ytitle, psym=-8, /ylog, yrange=[1d3,5d6],/ystyle, $                           
         xrange=[22.5,29.5],/xstyle, charsize=2., xthick=2, ythick=2, xmargin=[10,7], ymargin=[6,6]                                                                       
     polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,5d6], /data, noclip=0                                                                                     
     errplot, mags, densities-density_errs, densities+density_errs, color=0, thick=2                                                                                      
     oplot, mags, densities, psym=-8, color=0, thick=2     ; have to do this again because polyfill covers it up   
     oplot, mags, field_densities, color=0, linestyle=2, thick=2                                                                                                                                                                                                                                                                                                                                    
     LEGEND, ['blob','GOODS-S field'], /left, /top, color=[0,0], textcolor=0, linestyle=[0,2], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5 
     axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=2, charthick=2, xtickformat="(A1)", xthick=2                                                    
     axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=2                                                                 
     axis, color=0, yaxis=0, /ylog, yrange=[1d3,5d6], /ystyle,  /ynozero, /data, charsize=2, charthick=2, ytickformat="(A1)", ythick=2                                    
     axis, color=0, yaxis=1, /ylog, yrange=[1d3,5d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=2                                                 
  END          



























;----------------------------------MAKING PLOTS---------------------------------------------------------------------------------------------------------------------------;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO colordiagram, blobgals, fieldgals, datared, datamid, datablue, redname, midname, bluename, m, b, cap, blobname, radius                                             ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; colordiagram procedure                                                                                                                          ;                     ;
  ;   description                                                                                                                                   ;                     ;************
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker                                                                                                                                           ;;
   ; separate out the magnitudes to define colors                                                                                                                         ;
      redmags = datared[blobgals].MAG_ISO                                                                                                                                ;;
      midmags = datamid[blobgals].MAG_ISO                                                                                                                                ;;
      bluemags = datablue[blobgals].MAG_ISO                                                                                                                              ;;
      field_redmags = datared[fieldgals].MAG_ISO                                                                                                                         ;;
      field_midmags = datamid[fieldgals].MAG_ISO                                                                                                                         ;;
      field_bluemags = datablue[fieldgals].MAG_ISO                                                                                                                       ;;
   ; make the colors                                                                                                                                                      ;
      color_x = bluemags - midmags                                                                                                                                       ;;
      color_y = midmags - redmags                                                                                                                                        ;;
      field_color_x = field_bluemags - field_midmags                                                                                                                     ;;
      field_color_y = field_midmags - field_redmags                                                                                                                      ;;
   ; make the plot                                                                                                                                                        ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('colors', blobname, radius=radius)                                                                                                            ;;
        xtitle = bluename + ' - ' + midname                                                                                                                              ;;
        ytitle = midname + ' - ' + redname                                                                                                                               ;;
     plot, color_x, color_y, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, psym=8, xrange=[-5,5],/xstyle, yrange=[-5,5],/ystyle, charsize=1.5, $   ;; 
        thick=2, ymargin=[4,4]                                                                                                                                           ;;
     oplot, field_color_x, field_color_y, color=0, psym=3                                                                                                                ;;
     LEGEND, ['blob galaxy','field galaxy'], /left, /top, color=0, textcolor=0, psym=[8,3], charsize=1, charthick=1, /box, outline_color=0                               ;;
   ; add the cut line                                                                                                                                                     ;
      xvalues = 0.01*findgen(100000) - 500.                                                                                                                              ;;
      cut = xvalues*m + b                                                                                                                                                ;;
      cut[WHERE(cut GT cap)] = cap                                                                                                                                       ;;
      oplot, xvalues, cut, color=50, linestyle=2, thick=2                                                                                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO densityplot, nbins, mags, field_densities, field_density_errs, densities, density_errs, hst_densities, hst_density_errs, blobname, nApertures, radius, band, cuts  ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; densityplot procedure                                                                                                                           ;                     ;******make HST stuff look better
  ;   Makes plots of blob and field densities as a function of magnitude.                                                                           ;                     ;
  ; INPUTS: blobname - name of blob whose density is being plotted (for plot title)                                                                 ;                     ;
  ;         nbins - number of magnitude bins being plotted on the x-axis                                                                            ;                     ;
  ;         mags - array of magnitudes (data for the x-axis)                                                                                        ;                     ;
  ;         field_densities - array of field density measurements (data for the y-axis)                                                             ;                     ;
  ;         field_density_errors - array of errors on field density measurements (to be overplotted as a swath on top of field_densities)           ;                     ;
  ;         densities - array of blob density measurements (also data for the y-axis)                                                               ;                     ;
  ;         density_errors - array of errors on blob density measurements (to be overplotted on top of densities using errplot)                     ;                     ;
  ;         apsnameplot - a string showing the number of random apertures used in the analysis (for plot title)                                     ;                     ;
  ;         xtitle - a string containing the label for the x-axis of the plot                                                                       ;                     ;
  ;         cuts - a string, 'y' or 'n', showing whether the densities used are for raw data or post-color-cut data (for plot title)                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
     loadct, 0                                                                                                                                                           ;;
   ; set up vertices of polygon for polyfill                                                                                                                              ;
     xpoints = dblarr(2*nbins)                                                                                                                                           ;;
     xpoints[0:nbins-1] = mags                                                                                                                                           ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       xpoints[i+nbins] = mags[nbins-1-i]                                                                                                                                ;;
     ENDFOR                                                                                                                                                              ;;
     ; y is specific to each catalog (band)                                                                                                                               ;
     ypoints=dblarr(2*nbins)                                                                                                                                             ;;
     ypoints[0:nbins-1] = field_densities-field_density_errs                                                                                                             ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       ypoints[i+nbins] = field_densities[nbins-1-i]+field_density_errs[nbins-1-i]                                                                                       ;;
     ENDFOR                                                                                                                                                              ;;
   ; make the main plot comparing blob to field with error bars                                                                                                           ;
     ; set up titles for plot and axes                                                                                                                                    ;
        plottitle = titlemaker('overdensity', blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts)                                                      ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle='N!Igal!N per deg!E2!N per 0.5 mag'                                                                                                                       ;;
     plot, mags, densities, background=255, color=0, title=plottitle, xtitle=xtitle, ytitle=ytitle, psym=-8, /ylog, yrange=[1d3,5d6],/ystyle, $                          ;;
         xrange=[22.5,29.5],/xstyle, charsize=2., xthick=2, ythick=2, xmargin=[10,7], ymargin=[6,6]                                                                      ;;
     polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,5d6], /data, noclip=0                                                                                    ;;
     errplot, mags, densities-density_errs, densities+density_errs, color=0, thick=2                                                                                     ;;
     oplot, mags, densities, psym=-8, color=0, thick=2     ; have to do this again because polyfill covers it up                                                          ;
     oplot, mags, field_densities, color=0, linestyle=2, thick=2                                                                                                         ;;
     loadct, 39                                                                                                                                                          ;;
     oplot, mags, hst_densities, color=50, linestyle=5, thick=2                                                                                                          ;;
     errplot, mags, hst_densities-hst_density_errs, hst_densities+hst_density_errs, color=50, thick=2                                                                    ;;
     LEGEND, ['blob','blob field', 'GOODS-S field'], /left, /top, color=[0,0,50], textcolor=0, linestyle=[0,2,5], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5               ;;
     axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=2, charthick=2, xtickformat="(A1)", xthick=2                                                   ;;
     axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=2                                                                ;;
     axis, color=0, yaxis=0, /ylog, yrange=[1d3,5d6], /ystyle,  /ynozero, /data, charsize=2, charthick=2, ytickformat="(A1)", ythick=2                                   ;;
     axis, color=0, yaxis=1, /ylog, yrange=[1d3,5d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=2                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO densityplotsmall, nbins, mags, field_densities, field_density_errs, densities, density_errs, blobname, nApertures, radius, band, cuts                              ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; densityplot procedure                                                                                                                           ;                     ;
  ;   Makes plots of blob and field densities as a function of magnitude.                                                                           ;                     ;
  ; INPUTS: blobname - name of blob whose density is being plotted (for plot title)                                                                 ;                     ;
  ;         nbins - number of magnitude bins being plotted on the x-axis                                                                            ;                     ;
  ;         mags - array of magnitudes (data for the x-axis)                                                                                        ;                     ;
  ;         field_densities - array of field density measurements (data for the y-axis)                                                             ;                     ;
  ;         field_density_errors - array of errors on field density measurements (to be overplotted as a swath on top of field_densities)           ;                     ;
  ;         densities - array of blob density measurements (also data for the y-axis)                                                               ;                     ;
  ;         density_errors - array of errors on blob density measurements (to be overplotted on top of densities using errplot)                     ;                     ;
  ;         apsnameplot - a string showing the number of random apertures used in the analysis (for plot title)                                     ;                     ;
  ;         xtitle - a string containing the label for the x-axis of the plot                                                                       ;                     ;
  ;         cuts - a string, 'y' or 'n', showing whether the densities used are for raw data or post-color-cut data (for plot title)                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
     loadct, 0                                                                                                                                                           ;;
   ; set up vertices of polygon for polyfill                                                                                                                              ;
     xpoints = dblarr(2*nbins)                                                                                                                                           ;;
     xpoints[0:nbins-1] = mags                                                                                                                                           ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       xpoints[i+nbins] = mags[nbins-1-i]                                                                                                                                ;;
     ENDFOR                                                                                                                                                              ;;
     ; y is specific to each catalog (band)                                                                                                                               ;
     ypoints=dblarr(2*nbins)                                                                                                                                             ;;
     ypoints[0:nbins-1] = field_densities-field_density_errs                                                                                                             ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       ypoints[i+nbins] = field_densities[nbins-1-i]+field_density_errs[nbins-1-i]                                                                                       ;;
     ENDFOR                                                                                                                                                              ;;
   ; make the main plot comparing blob to field with error bars                                                                                                           ;
     ; set up titles for plot and axes                                                                                                                                    ;
        plottitle = titlemaker('overdensity', blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts)                                                      ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle='N!Igal!N per deg!E2!N per 0.5 mag'                                                                                                                       ;;
     plot, mags, densities, background=255, color=0, title=plottitle, xtitle=xtitle, ytitle=ytitle, psym=-8, /ylog, yrange=[1d3,1d6],/ystyle, $                          ;;
         xrange=[22.5,29.5],/xstyle, charsize=2.5, xthick=2, ythick=2, xmargin=[10,7], ymargin=[6,6], charthick=2                                                                      ;;
     polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,1d6], /data, noclip=0                                                                                    ;;
     errplot, mags, densities-density_errs, densities+density_errs, color=0, thick=4                                                                                     ;;
     oplot, mags, densities, psym=-8, color=0, thick=4     ; have to do this again because polyfill covers it up                                                          ;
     oplot, mags, field_densities, color=0, linestyle=2, thick=4                                                                                                         ;;                                                                  ;;
     LEGEND, ['blob','field (GOODS-S)'], /right, /top, color=[0,50], textcolor=0, linestyle=[0,2], thick=4., charsize=2.5, /box, outline_color=0., charthick=2, number=0.1            ;;
     axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=4, charthick=2, xtickformat="(A1)", xthick=4                                                   ;;
     axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=4                                                                ;;
     axis, color=0, yaxis=0, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=4, charthick=2, ytickformat="(A1)", ythick=4                                   ;;
     axis, color=0, yaxis=1, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=4                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;   
  PRO statlum, mags, binsize, blob_ngal_binned, blob_ngal_binned_cuts, densities, field_densities, ap_radius, blobname, nApertures, ap_radius_arcsec, band               ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; statlum procedure                                                                                                                               ;                     ;***********
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; define the statistical luminosity function                                                                                                                           ;
      overdensity = densities - field_densities                                                                                                                          ;;
      n_overdensity = overdensity * (!dPI*ap_radius^2.)  ; convert to number                                                                                              ;
   ; set up vertices of shaded area                                                                                                                                       ;
     xvertices = fltarr(2.0*n_elements(mags) + 2.)                                                                                                                       ;;
     yvertices = fltarr(2.0*n_elements(mags) + 2.)                                                                                                                       ;;
     xvertices[0] = mags[0]                                                                                                                                              ;;
     xvertices[-1] = mags[-1]                                                                                                                                            ;;
     yvertices[0] = 0                                                                                                                                                    ;;
     yvertices[-1] = 0                                                                                                                                                   ;;
     xvertices[1] = mags[0]                                                                                                                                              ;;
     xvertices[-2] = mags[-1]                                                                                                                                            ;;
     FOR i=1, n_elements(mags) -1 DO BEGIN                                                                                                                               ;;
       xvertices[2*i] = mags[i]-(binsize/2.0)                                                                                                                            ;;
       xvertices[2*i+1] = mags[i]-(binsize/2.0)                                                                                                                          ;;
     ENDFOR                                                                                                                                                              ;;
     FOR i=1, n_elements(mags) -1 DO BEGIN                                                                                                                               ;;
       yvertices[2*i+1] = blob_ngal_binned_cuts[i]                                                                                                                       ;;
       yvertices[2*i+2] = blob_ngal_binned_cuts[i]                                                                                                                       ;;
     ENDFOR                                                                                                                                                              ;;
   ; now plot statistical luminosity function                                                                                                                             ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('statlum', blobname, nApertures=nApertures, radius=ap_radius_arcsec, band=band)                                                               ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle = 'N!Igal!N'                                                                                                                                              ;;
     plot, mags, n_overdensity, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, linestyle=2, yrange=[0,10],/ystyle, xrange=[22,30],/xstyle, $        ;;
       charsize=3, thick=2, xmargin=[10,7], ymargin=[4,4]                                                                                                                ;;
     polyfill, xvertices,yvertices, color=220, clip=[22,0,30,10], /data, noclip=0                                                                                        ;;
     oplot, mags, blob_ngal_binned, color=0, psym=10, thick=2                                                                                                            ;;
     oplot, mags, blob_ngal_binned_cuts, color=0, psym=10, thick=2                                                                                                       ;;
     oplot, mags, n_overdensity, color=0, thick=2, linestyle=2                                                                                                           ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO makehist, ap_totals, blob_ngal, blobname, nApertures, radius, cuts     ; maybe put optional band keyword in here or mag keyword or smth                             ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; makehist procedure                                                                                                                              ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;***************
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; set up titles for plot and axes                                                                                                                                      ;
      title = titlemaker('aphist', blobname, nApertures=nApertures, radius=radius, cuts=cuts)                                                                            ;;
      xtitle = 'n galaxies in aperture'                                                                                                                                  ;;
      ytitle='number of apertures'                                                                                                                                       ;;
   plothist, ap_totals, background=255, color=0, axiscolor=0, title=title, xtitle=xtitle, ytitle=ytitle, bin=1, xrange=[0,blob_ngal+10.],/xstyle, charsize=1.5, $        ;;
     thick=2, ymargin=[4,4]                                                                                                                                              ;;
   oplot, [blob_ngal, blob_ngal], [0., nApertures], linestyle=2, thick=2, color=0.                                                                                       ;;
   LEGEND, ['blob','field'], /center, /top, color=0, textcolor=0, linestyle=[0,2], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5                ;;
  ; xyouts, blob_ngal+1., nApertures/10., blobgaltext, color=0, charsize=1.5                                                                                              ;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO Rmultirad, ap_radius_arcsec, mags, radiusresults, index, blobname, nApertures, band, cuts                                                                          ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; Rmultirad procedure                                                                                                                             ;                     ;
  ;   description                                                                                                                                   ;                     ;*************
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; set up an array of colors to plot each radius in                                                                                                                     ;
     n_radii = n_elements(ap_radius_arcsec)                                                                                                                              ;;
     indices = findgen(n_radii)                                                                                                                                          ;;
     colors = indices*254./n_radii                                                                                                                                       ;;
   ; now make the plot                                                                                                                                                    ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('Rmultirad', blobname, nApertures=nApertures, band=band, cuts=cuts)                                                                           ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle = 'overdensity factor'                                                                                                                                    ;;
     plot, mags, radiusresults[0,*,index,0], background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, yrange=[0,50],/ystyle, charsize=2, thick=2, $           ;;
       ymargin=[4,4]                                                                                                                                                     ;;
     errplot, mags, radiusresults[0,*,index,0]-radiusresults[0,*,index,1], radiusresults[0,*,index,0]+radiusresults[0,*,index,1], color=0, thick=2                       ;;
     FOR r=1, n_radii-1. DO BEGIN                                                                                                                                        ;;
       oplot, mags, radiusresults[r,*,index,0], color = colors[r], thick=2                                                                                               ;;
       errplot, mags, radiusresults[r,*,index,0]-radiusresults[r,*,index,1], radiusresults[r,*,index,0]+radiusresults[r,*,index,1], color=colors[r], thick=2             ;;
     ENDFOR                                                                                                                                                              ;;
     LEGEND, strcompress(string(fix(ap_radius_arcsec)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., $            ;;
       number=0.5, charthick=1.5, thick=2                                                                                                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO Rmultimag, mags, ap_radius_arcsec, radiusresults, index, blobname, nApertures, band, cuts                                                                          ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; Rmultimag procedure                                                                                                                             ;                     ;
  ;   description                                                                                                                                   ;                     ;***********
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; set up an array of colors to plot each magnitude bin in                                                                                                              ;
     indices = findgen(n_elements(mags))                                                                                                                                 ;;
     colors = indices*254./float(n_elements(mags))                                                                                                                       ;;
   ; now make the plot                                                                                                                                                    ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('Rmultimag', blobname, nApertures=nApertures, band=band, cuts=cuts)                                                                           ;;
        xtitle='aperture radius (arcseconds)'                                                                                                                            ;;
        ytitle='overdensity factor'                                                                                                                                      ;;
     ; set up an empty plot (needs to have no data in case the first thing plotted is NaN):                                                                               ;
     plot, ap_radius_arcsec, radiusresults[*,1,index,0], background=255, color=colors[1], title=title, xtitle=xtitle, ytitle=ytitle, xrange=[0,30],/xstyle, $            ;;
       yrange=[0,100],/ystyle, charsize=2, ymargin=[4,4], /nodata                                                                                                        ;;
     ; now add data on top of the empty plot:                                                                                                                             ;
     FOR i=0, n_elements(mags)-1 DO BEGIN     ; need to make sure that NaNs are thrown out (because where field density is 0, we divided by 0)                            ;
       y = radiusresults[*,i,index,0]                                                                                                                                    ;;
       yerr = radiusresults[*,i,index,1]                                                                                                                                 ;;
       keep = WHERE(finite(y) EQ 1, nkeep)                                                                                                                               ;;
       IF (nkeep GT 0) THEN BEGIN                                                                                                                                        ;;
         oplot, ap_radius_arcsec(keep), y(keep), color=colors[i], thick=2                                                                                                ;;
         errplot, ap_radius_arcsec(keep), y(keep)-yerr(keep), y(keep)+yerr(keep), color=colors[i], thick=2                                                               ;;
       ENDIF                                                                                                                                                             ;;
     ENDFOR                                                                                                                                                              ;;
     ; add a legend:                                                                                                                                                      ;
     LEGEND, strcompress(string(fix(mags)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., number=0.5, $            ;;
        charthick=1.5, thick=2                                                                                                                                           ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO makereg, data, subset, blobname, radius, cuts                                                                                                                      ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; makereg procedure                                                                                                                               ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTE: Uses titlemaker function.                                                                                                                 ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   savethis = ''                                                                                                                                                         ;;
   READ, savethis, PROMPT='Make region file? (y/n)'                                                                                                                      ;;
   IF (savethis EQ 'y') THEN BEGIN                                                                                                                                       ;;
     cat = data[subset]                                                                                                                                                  ;;
     filename = titlemaker('blobgals', blobname, radius=radius, cuts=cuts, filetype='.reg')                                                                              ;;
     openw, 1, filename                                                                                                                                                  ;;
     FOR i=0, n_elements(cat)-1 DO BEGIN                                                                                                                                 ;;
       printf, 1,'J2000; circle ', cat[i].alpha_j2000, cat[i].delta_j2000, ' 10p '    ;#text={', cat[i].number,'}'                                                        ;
     ENDFOR                                                                                                                                                              ;;
     close, 1                                                                                                                                                            ;;
   ENDIF ELSE IF ((savethis NE 'y') AND (savethis NE 'n')) THEN BEGIN                                                                                                    ;;
     WHILE ((savethis NE 'y') AND (savethis NE 'n')) DO BEGIN                                                                                                            ;;
       print, 'Invalid response. Choose y or n.'                                                                                                                         ;;
       READ, savethis, PROMPT='Make region file? (y/n)'                                                                                                                  ;;
     ENDWHILE                                                                                                                                                            ;;
   ENDIF ELSE IF (savethis EQ 'n') THEN BEGIN                                                                                                                            ;;
        print, 'Region file not created.'                                                                                                                                ;;
   ENDIF                                                                                                                                                                 ;;
  END                                                                                                                                                                    ;;                                                                                                                                                                     ; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO plotsave, saveanyplots, namestring                                                                                                                                 ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; plotsave procedure                                                                                                                              ;                     ;
  ;   Asks user whether or not to save a plot and saves said plot if the user chooses to do so.                                                     ;                     ;
  ; INPUTS: saveanyplots - a string, 'y' or 'n', determining if the user should be asked whether or not to save this plot at all (code will only    ;                     ;
  ;                          run if this is set to 'y')                                                                                             ;                     ;
  ;         namestring - a string containing the name of the file containing the plot if the user chooses to save it                                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    IF (saveanyplots EQ 'y') THEN BEGIN                                                                                                                                  ;;
      saveplot = ''                                                                                                                                                      ;;
      READ, saveplot, PROMPT='Save plot? (y/n)'                                                                                                                          ;;
      IF (saveplot EQ 'y') THEN BEGIN                                                                                                                                    ;;
         write_png, namestring, tvrd(/true)                                                                                                                              ;;
      ENDIF ELSE IF ((saveplot NE 'y') AND (saveplot NE 'n')) THEN BEGIN                                                                                                 ;;
        WHILE ((saveplot NE 'y') AND (saveplot NE 'n')) DO BEGIN                                                                                                         ;;
          print, 'Invalid response. Choose y or n.'                                                                                                                      ;;
          READ, saveplot, PROMPT='Save plot? (y/n)'                                                                                                                      ;;
        ENDWHILE                                                                                                                                                         ;;
      ENDIF ELSE IF (saveplot EQ 'n') THEN BEGIN                                                                                                                         ;;
        print, 'Plot not saved.'                                                                                                                                         ;;
      ENDIF                                                                                                                                                              ;;
    ENDIF   ; user wants to save at least some plots                                                                                                                      ;
  END                                                                                                                                                                    ;;                                                                                                                                                                   ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  FUNCTION titlemaker, plotname, blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts, filetype=filetype                                                 ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; titlemaker function                                                                                                                             ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: plotname - string containing the name of the kind of plot being made (ie, "overdensity," "color-color," etc)                            ;                     ;
  ;         blobname - string                                                                                                                       ;                     ;
  ;         nApertures - float                                                                                                                      ;                     ;
  ;         plotname - string containing the name of the kind of plot being made (ie, "overdensity," "color-color," etc)                            ;                     ;
  ;         radius - float                                                                                                                          ;                     ;
  ;         band - string (F140W, F606W, etc)                                                                                                       ;                     ;
  ;         cuts - string/keyword                                                                                                                   ;                     ;
  ;         filetype - specifies that this is going to be a filename if set; string that gives extension of the file                                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    ; first, make checks that see if the optional keywords are provided                                                                                                   ; 
    apcheck = n_elements(nApertures)                                                                                                                                     ;;
    radcheck = n_elements(radius)                                                                                                                                        ;;
    bandcheck = n_elements(band)                                                                                                                                         ;;
    cutscheck = n_elements(cuts)                                                                                                                                         ;;
    filecheck = n_elements(filetype)                                                                                                                                     ;;
                                                                                                                                                                         ;;
    ; next make empty strings set up for each element of a title                                                                                                          ;
    apsName = ''                                                                                                                                                         ;;                           
    radName = ''                                                                                                                                                         ;;
    bandName = ''                                                                                                                                                        ;;
    cutsName = ''                                                                                                                                                        ;;
    titleName = ''                                                                                                                                                       ;;
                                                                                                                                                                         ;;
    ; now make the final string that will eventually be returned                                                                                                          ;
    titlestring = ''                                                                                                                                                     ;;
                                                                                                                                                                         ;;
    ; since the cuts keyword can only be 'y' or 'n', if it's set, make sure it's one of those things                                                                      ;
    IF (cutscheck NE 0) THEN BEGIN                                                                                                                                       ;;
      IF ( (cuts NE 'y') AND (cuts NE 'n') ) THEN BEGIN                                                                                                                  ;;
        WHILE ( (cuts NE 'y') AND (cuts NE 'n') ) DO BEGIN        ; make sure the keyword is valid                                                                        ;
          print, 'Invalid cut keyword provided. Choose "y" or "n" only.'                                                                                                 ;;
          cuts = ''                                                                                                                                                      ;;
          READ, cuts, PROMPT='Color cuts?'                                                                                                                               ;;
        ENDWHILE                                                                                                                                                         ;;
      ENDIF ; cuts was set to something other than y or n                                                                                                                 ;
    ENDIF   ; the cuts keyword was provided at all                                                                                                                        ;
                                                                                                                                                                         ;;
    ; now, determine whether to make a plot title or a file name:                                                                                                         ;
                                                                                                                                                                         ;;
    IF (filecheck NE 0) THEN BEGIN            ; if the user wants a filename                                                                                              ;
      ;                                                                                                                                                                   ;
      titlestring = string(plotname) + '_' + string(blobname)                                                                                                            ;;
      ; now tack on any extra keywords one by one                                                                                                                         ;
      ;                                                                                                                                                                   ;
      ; start with nApertures                                                                                                                                             ;
      IF (apcheck NE 0) THEN BEGIN                                                                                                                                       ;; 
        apsName = strcompress(string(fix(nApertures)), /remove) + 'aps'                                                                                                  ;;
        titlestring = titlestring + '_' + string(apsName)                                                                                                                ;;
      ENDIF                                                                                                                                                              ;; 
      ;                                                                                                                                                                   ;                                                                                                                                                                   ;
      ; radius                                                                                                                                                            ;
      IF (radcheck NE 0) THEN BEGIN                                                                                                                                      ;; 
        radName = strcompress(string(fix(radius)), /remove) + 'arcsec'                                                                                                   ;;
        titlestring = titlestring + '_' + string(radName)                                                                                                                ;;
      ENDIF                                                                                                                                                              ;; 
      ;                                                                                                                                                                   ;
      ; band                                                                                                                                                              ;
      IF (bandcheck NE 0) THEN BEGIN                                                                                                                                     ;; 
        bandName = band                                                                                                                                                  ;;
        titlestring = titlestring + '_' + string(bandName)                                                                                                               ;;
      ENDIF                                                                                                                                                              ;;
      ;                                                                                                                                                                   ;
      ; cuts                                                                                                                                                              ;
      IF (cutscheck NE 0) THEN BEGIN                                                                                                                                     ;;
        IF (cuts EQ 'y') THEN BEGIN        ; make a string for the case that there are color cuts                                                                         ;
          cutsName = 'cuts'                                                                                                                                              ;;
        ENDIF ELSE IF (cuts EQ 'n') THEN BEGIN        ; make a string for data without cuts                                                                               ;
          cutsName = 'raw'                                                                                                                                               ;;
        ENDIF         ; now we've sorted out what string to use to describe the cuts keyword                                                                              ;
       titlestring = titlestring + '_' + string(cutsName)                                                                                                                ;;
      ENDIF                                                                                                                                                              ;;
      ;                                                                                                                                                                   ;
      ; finally, add the file extension to the end of the string                                                                                                          ;
      titlestring = titlestring + string(filetype)                                                                                                                       ;;
      ; now we have a valid file name                                                                                                                                     ;
                                                                                                                                                                         ;;
    ENDIF ELSE IF (filecheck EQ 0) THEN BEGIN         ; if the user is making a title for a plot                                                                          ;
      ;                                                                                                                                                                   ;
      ; generate a basic plot title based on what's being plotted                                                                                                         ;
      IF (plotname EQ 'colors') THEN titleName = 'Color-Color Diagram for '                                                                                              ;;
      IF (plotname EQ 'overdensity') THEN titleName = 'Galaxy Overdensity in Blob Region for '                                                                           ;;
      IF (plotname EQ 'statlum') THEN titleName = 'Statistical Luminosity Function for Galaxies in '                                                                     ;;
      IF (plotname EQ 'aphist') THEN titleName = 'Number of Galaxies in Apertures for '                                                                                  ;;
      IF (plotname EQ 'Rmultirad') THEN titleName = 'Galaxy Overdensity in Different Radii for '                                                                         ;;
      IF (plotname EQ 'Rmultimag') THEN titleName = 'Dependence of Overdensity Factor on Aperture Radius for '                                                           ;;
      titlestring = string(titleName) + string(blobname)                                                                                                                 ;;
      ;                                                                                                                                                                   ;
      ; Now, if there are any keywords, add things to the basic plot title:                                                                                               ;
      ;                                                                                                                                                                   ;
      IF ( (apcheck NE 0) OR (radcheck NE 0) OR (bandcheck NE 0) OR (cutscheck NE 0) ) THEN BEGIN                                                                        ;;
        ; add a comma and a line after the main title, since this is extra stuff                                                                                          ;
        titlestring = titlestring + ',!C'                                                                                                                                ;;
        titlestring = string(titlestring)                                                                                                                                ;;
        ;                                                                                                                                                                 ;
        ; First, see if nApertures and/or radius are specified. Since they go hand-in-hand, treat them together:                                                          ;
        IF (radcheck NE 0) THEN radName = strcompress(string(fix(radius)), /remove) + '"'     ; if radius is provided, here's a string describing it                      ;
        IF (apcheck NE 0) THEN BEGIN            ; if nApertures is provided                                                                                               ;
          apsName = strcompress(string(fix(nApertures)), /remove) + ' apertures'    ; string describing nApertures                                                        ;
          IF (radcheck NE 0) THEN BEGIN        ; if nApertures and radius are both provided                                                                               ;
            titlestring = titlestring + string(apsName) + ' with ' + string(radName) + ' radii'                                                                          ;;
          ENDIF ELSE IF (radcheck EQ 0) THEN BEGIN            ; if nApertures is provided but not radius                                                                  ;
            titlestring = titlestring + string(apsName)                                                                                                                  ;;
          ENDIF                                                                                                                                                          ;;
        ENDIF ELSE IF ( (apcheck EQ 0) AND (radcheck NE 0) ) THEN BEGIN     ; if nApertures is not provided but radius is                                                 ;
          titlestring = titlestring + 'aperture radius of ' + string(radName)                                                                                            ;;
        ENDIF         ; if neither nApertures nor radius are provided, nothing else is added to the title string yet                                                      ;
        ;                                                                                                                                                                 ;
        ; now tack on the name of the band if that has been provided                                                                                                      ;
        IF (bandcheck NE 0) THEN BEGIN                                                                                                                                   ;;
          bandName = string(band)                                                                                                                                        ;;
          IF ( (apcheck NE 0) OR (radcheck NE 0) ) THEN BEGIN     ; if anything else is provided before band, add a comma                                                 ;
            titlestring = titlestring + ', '                                                                                                                             ;;
          ENDIF                                                                                                                                                          ;;
          titlestring = titlestring + string(bandName)                                                                                                                   ;;
        ENDIF    ; the band is provided                                                                                                                                   ;
        ;                                                                                                                                                                 ;
        ; now add info about color cuts if that has been provided                                                                                                         ;
        IF (cutscheck NE 0) THEN BEGIN                                                                                                                                   ;;
          IF (cuts EQ 'y') THEN BEGIN        ; make a string for the case that there are color cuts                                                                       ;
            cutsName = 'with color cuts'                                                                                                                                 ;;
          ENDIF ELSE IF (cuts EQ 'n') THEN BEGIN        ; make a string for data without cuts                                                                             ;
            cutsName = 'raw (no cuts)'                                                                                                                                   ;;
          ENDIF         ; sorts out what string to use to describe the cuts keyword                                                                                       ;
          ; now that the string is set up, add it to the title                                                                                                            ;
          IF ( (apcheck NE 0) OR (radcheck NE 0) OR (bandcheck NE 0) ) THEN BEGIN                                                                                        ;;
            titlestring = titlestring + ', '         ; add a comma if anything comes before the cuts text                                                                 ;
          ENDIF                                                                                                                                                          ;;
          titlestring = string(titlestring) + string(cutsName)       ; tack on the cuts text                                                                              ;
        ENDIF      ; the cuts keyword is provided                                                                                                                         ;
        ;                                                                                                                                                                 ;
      ENDIF     ; there are any keywords specified at all                                                                                                                 ;
      ;                                                                                                                                                                   ;
    ENDIF    ; this is the title for a plot                                                                                                                               ;
                                                                                                                                                                         ;;
    RETURN, string(titlestring)       ; provide the final title for the plot or file                                                                                      ;
                                                                                                                                                                         ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
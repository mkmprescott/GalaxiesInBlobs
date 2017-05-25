;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
; COLORANALYZER.PRO 
;  Contains all modules used for analyzing the colors of sources in blobs and in fields. 
;  CONTENTS: 
;      colordiagram-----------------------procedure
;      colordiagramHST--------------------procedure
;      colorsort--------------------------procedure
;      plotsave---------------------------procedure     ?????
;      titlemaker-------------------------function      ?????
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;




  PRO colordiagram, blobgals, fieldgals, datared, datamid, datablue, redname, midname, bluename, m, b, cap, blobname, radius                                               
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;
  ; colordiagram procedure                                                                                              
  ;   description                                                                 
  ; INPUTS: etc                                                                                          
  ;         etc                                                                                         
  ; NOTES:                                                                              
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;
   FORWARD_FUNCTION titlemaker                                                                                                                                             
   ; separate out the magnitudes to define colors  
      redmags = datared[blobgals].MAG_ISO                                                                                                                                  
      midmags = datamid[blobgals].MAG_ISO                                                                                                                                  
      bluemags = datablue[blobgals].MAG_ISO                                                                                                                                
      field_redmags = datared[fieldgals].MAG_ISO                                                                                                                           
      field_midmags = datamid[fieldgals].MAG_ISO                                                                                                                           
      field_bluemags = datablue[fieldgals].MAG_ISO                                                                                                                         
   ; make the colors  
      color_x = bluemags - midmags                                                                                                                                         
      color_y = midmags - redmags                                                                                                                                          
      field_color_x = field_bluemags - field_midmags                                                                                                                       
      field_color_y = field_midmags - field_redmags                                                                                                                        
   ; make the plot  
     ; set up titles for plot and axes 
        title = titlemaker('colors', blobname, radius=radius)                                                                                                              
        xtitle = bluename + ' - ' + midname                                                                                                                                
        ytitle = midname + ' - ' + redname                                                                                                                                 
     plot, color_x, color_y, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, psym=8, xrange=[-5,5],/xstyle, yrange=[-5,5],/ystyle, charsize=1.5, $      
        thick=2, ymargin=[4,4]                                                                                                                                             
     oplot, field_color_x, field_color_y, color=0, psym=3                                                                                                                  
     LEGEND, ['blob galaxy','field galaxy'], /left, /top, color=0, textcolor=0, psym=[8,3], charsize=1, charthick=1, /box, outline_color=0                                 
   ; add the cut line   
      xvalues = 0.01*findgen(100000) - 500.                                                                                                                                
      cut = xvalues*m + b                                                                                                                                                  
      cut[WHERE(cut GT cap)] = cap                                                                                                                                         
      oplot, xvalues, cut, color=50, linestyle=2, thick=2                                                                                                                  
  END                                                                                                                                                                       





      PRO HSTdensity
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      ; HSTdensity procedure
      ;   description 
      ; INPUTS: etc
      ;         etc   
      ; NOTES: 
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      FORWARD_FUNCTION fluxtomag
      FORWARD_FUNCTION colorsort

      ; 5 fields total - aegis, cosmos, goodsn, goodss, uds
      ;aegis = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/AEGIS/Catalog/aegis_3dhst.v4.1.cat.FITS', 1)
      ;cosmos = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/COSMOS/Catalog/cosmos_3dhst.v4.1.cat.FITS',1)
      ;goodsn = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/GOODSN/Catalog/goodsn_3dhst.v4.1.cat.FITS', )
      goodss = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/GOODSS/Catalog/goodss_3dhst.v4.1.cat.FITS', 1)
      Zgoodss = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/GOODSS/Eazy/goodss_3dhst.v4.1.zout.FITS', 1)
      ;uds = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/UDS/Catalog/uds_3dhst.v4.2.cat.FITS', 1)

      ; bands are F140W, F814W, and F606W or F475W

      ; AEGIS BANDS: F140W, F606W, F814W, F125W, F160W
      ; COSMOS BANDS: F140W, F606W, F814W, F125W, F160W
      ; GOODSN BANDS: F140W, F435W, F606W, F775W, F850LP, F125W, F160W 
      ; GOODSS BANDS: F140W, F435W, F606W, F775W, F814Wcand, F850LP, F850LPcand, F125W, F160W 
      ; UDS BANDS: F140W, F606W, F814W, F125W, F160W

      fluxconverter, goodss.f_F140W, goodss.f_F814Wcand, goodss.f_F606W, goodss.f_F435W, goodssred, goodssmid, goodssblue, goodssbluealt

      z = Zgoodss.z_peak 
      zspec = Zgoodss.z_spec
      spec = WHERE(zspec NE -1.) 
      z[spec] = zspec[spec] 
      ; bad = WHERE(z LT 0) 

      ; PRG1: z=1.67, m=1.8, b=-1.0, c=2.0
      ; PRG2: z=2.27, m=2.5, b=-0.8, c=99.0
      ; PRG3: z=2.14, m=1.8, b=-1.0, c=2.0

      ;colordiagramHST, goodssred, goodssmid, goodssbluealt, 'F140W', 'F814W', 'F435W', z, 1.8, -1.0, 2.0, 1.67, 0.15, 'PRG1'
      ;plotsave, '3DHSTcolordiagram_PRG1.png'
      ;colordiagramHST, goodssred, goodssmid, goodssblue, 'F140W', 'F814W', 'F606W', z, 2.5, -0.8, 99.0, 2.27, 0.15, 'PRG2'
      ;plotsave, '3DHSTcolordiagram_PRG2.png'
      ;colordiagramHST, goodssred, goodssmid, goodssbluealt, 'F140W', 'F814W', 'F435W', z, 1.8, -1.0, 2.0, 2.14, 0.15, 'PRG3'
      ;plotsave, '3DHSTcolordiagram_PRG3.png'

      ; for PRG1 and PRG3
      colorcuts = colorsort(goodssblue-goodssmid, goodssmid-goodssred, 1.8, -1.0, 2.0) 
      zcut = z[colorcuts]
      ; for PRG2
      colorcutsalt = colorsort(goodssbluealt-goodssmid, goodssmid-goodssred, 2.5, -0.8, 99.0) 
      zaltcut = z[colorcutsalt]

      zhist, z, zcut, zaltcut

      use = WHERE(goodss.use_phot EQ 1.) 
      goodss = goodss[use]
      goodssred = goodssred[use]
      goodssmid = goodssmid[use]
      goodssblue = goodssblue[use]
      goodssbluealt = goodssbluealt[use]
      z = z[use]
      ;flags = -1.*(goodss.use_phot - 1.) 
      header = ['ID', 'ra', 'dec', 'F140W','F814W','F606W','F435W','z'] 
      write_csv, 'goodss.csv', goodss.id, goodss.ra, goodss.dec, goodssred, goodssmid, goodssblue, goodssbluealt, z, header=header

      END   ; HSTdensity procedure 



      PRO zhist, z, zcut, zaltcut
      ; histograms
      xtitle='z'
      ytitle='number of galaxies'
      ; PRG1 & PRG3: 
      title= 'Z Distribution for PRG1 and PRG3'
      plothist, z, background=255, color=0, axiscolor=0, title=title, xtitle=xtitle, ytitle=ytitle, bin=0.1, charsize=1.5, thick=2, ymargin=[4,4], xrange=[0,6], /xstyle
      plothist, zcut, bin=0.1, /overplot, /fill, fcolor=150 
      axis, color=0, xaxis=0, /data, charsize=2, xtickformat="(A1)", xticks=6, xminor=10
      axis, color=0, yaxis=0, /data, charsize=2, ytickformat="(A1)"
      plotsave, 'zdist_PRG1PRG3.png'
      ; PRG2: 
      title= 'Z Distribution for PRG2'
      plothist, z, background=255, color=0, axiscolor=0, title=title, xtitle=xtitle, ytitle=ytitle, bin=0.1, charsize=1.5, thick=2, ymargin=[4,4], xrange=[0,6], /xstyle
      plothist, zaltcut, bin=0.1, /overplot, /fill, fcolor=150 
      axis, color=0, xaxis=0, /data, charsize=2, xtickformat="(A1)", xticks=6, xminor=10
      axis, color=0, yaxis=0, /data, charsize=2, ytickformat="(A1)"
      plotsave, 'zdist_PRG2.png'
      END 



      PRO fluxconverter, redflux, midflux, blueflux, altblueflux, redmag, midmag, bluemag, altbluemag 
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      ; fluxconverter procedure
      ;   description 
      ; INPUTS: etc
      ;         etc   
      ; NOTES: 
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      FORWARD_FUNCTION fluxtomag
      redmag = fluxtomag(redflux)
      midmag = fluxtomag(midflux)
      bluemag = fluxtomag(blueflux)
      altbluemag = fluxtomag(altblueflux)
      END 



      FUNCTION fluxtomag, flux
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      ; fluxtomag function
      ;   description 
      ; INPUTS: etc
      ;         etc   
      ; NOTES: 
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      ;zeropoint is 25 so magAB = 25.0 - 2.5*log10(flux)
	magAB = 25.0 - 2.5*alog10(flux) 
	return, magAB
      END 



      PRO colordiagramHST, datared, datamid, datablue, redname, midname, bluename, z, m, b, cap, targetz, dz, blobname  
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      ; colordiagramHST procedure
      ;   description 
      ; INPUTS: etc
      ;         etc   
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      ; plotting symbol: make the psym=8 symbol a filled circle for our plots 
	  plotsym, 0, 1, /fill   
      ; make the colors
	  color_x = datablue - datamid
	  color_y = datamid - datared
      ; mark galaxies at the right redshift 
	  goodz = WHERE((z GE targetz-dz) AND (z LE targetz+dz))
	  color_x_good = datablue[goodz] - datamid[goodz] 
	  color_y_good = datamid[goodz] - datared[goodz] 
      ; make the plot 
	; set up titles for plot and axes
	    title = 'GOODS-S Color-Color Diagram for z of ' + blobname
	    xtitle = bluename + ' - ' + midname
	    ytitle = midname + ' - ' + redname 
	    redshift = strcompress(string(targetz, FORMAT='(F4.2)'), /remove) 
	    dshift = strcompress(string(dz,FORMAT='(F4.2)'), /remove) 
	    consistent = string('z within ' + redshift + ' !9+!3 ' + dshift) 
	plot, color_x, color_y, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, psym=3, xrange=[-5,5],/xstyle, yrange=[-5,5],/ystyle, charsize=1.5, thick=2, ymargin=[4,4] 
	oplot, color_x_good, color_y_good, color=0, psym=8
	LEGEND, [consistent, 'other galaxies'], /left, /bottom, color=0, textcolor=0, psym=[8,3], charsize=1, charthick=1, /box, outline_color=0 
      ; add the cut line 
	  xvalues = 0.01*findgen(100000) - 500.
	  cut = xvalues*m + b
	  cut[WHERE(cut GT cap)] = cap
	  oplot, xvalues, cut, color=50, linestyle=2, thick=2
      END 



      PRO plotsave, namestring 
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
      ; plotsave procedure
      ;   Asks user whether or not to save a plot and saves said plot if the user chooses to do so. 
      ; INPUTS: namestring - a string containing the name of the file containing the plot if the user chooses to save it 
      ;-------------------------------------------------------------------------------------------------------------------------------------------------
	  saveplot = ''
	  READ, saveplot, PROMPT='Save plot? (y/n)'
	  IF (saveplot EQ 'y') THEN BEGIN 
	    write_png, namestring, tvrd(/true) 
	  ENDIF ELSE IF ((saveplot NE 'y') AND (saveplot NE 'n')) THEN BEGIN 
	    WHILE ((saveplot NE 'y') AND (saveplot NE 'n')) DO BEGIN 
	      print, 'Invalid response. Choose y or n.' 
	      READ, saveplot, PROMPT='Save plot? (y/n)'
	    ENDWHILE
	  ENDIF ELSE IF (saveplot EQ 'n') THEN BEGIN
	    print, 'Plot not saved.'
	  ENDIF 
      END 



      FUNCTION colorsort, xcolor, ycolor, m, b, c 
      ;------------------------------------------------------------------------------------------------------------------------------------------------
      ; colorsort function
      ;   Applies a simple linear (y >= mx + b) color cut to a catalog made by SExtractor. 
      ; INPUTS: data - a SExtractor catalog of sources
      ;         xcolor - an array containing the blue-middle band color of every source in the data catalog 
      ;         ycolor - an array containing the middle-red band color of every source in the data catalog 
      ;         m - the slope of a line marking a simple color cut (unitless)  
      ;         b - the y-intercept of the simple color cut line (unitless)    
      ;         c - a "cap" value for color cuts, above which all colors are considered "good" (unitless)
      ; OUTPUT: good - the index numbers of all the sources in the data catalog with colors above the cut
      ; NOTES: The cut made is illustrated below. Asterisks mark "good" sources, while periods mark sources that get thrown out by the cut. 
      ;                               
      ;                               |   *     *    *      *
      ;                               |      *          *      *
      ;                             c-|   *    *   __________________
      ;                               |  *        / .     .    .  .
      ;                      ycolor   |    *     /.    .        .
      ;                               |         /   .      ..      .
      ;                               |     *  /     .  .    .   .
      ;                               | *     /   .   .   .    .    
      ;                               |   *  /       .   . .  .    .
      ;                               |_____/________________________
      ;                     (y = mx + b)----^       xcolor            
      ;------------------------------------------------------------------------------------------------------------------------------------------------
	good = WHERE ((ycolor GE (xcolor*m + b)) OR (ycolor GE c))
	RETURN, good
      END 
PRO HSTdensity
;-------------------------------------------------------------------------------------------------------------------------------------------------
; HSTdensity procedure
;   description 
; INPUTS: etc
;         etc   
; NOTES: 
;-------------------------------------------------------------------------------------------------------------------------------------------------
FORWARD_FUNCTION fluxtomag

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


fluxconverter, goodss.f_F140W, goodss.f_F814Wcand, goodss.f_F606W, goodssred, goodssmid, goodssblue
goodssbluealt = fluxtomag(goodss.f_F435W)

z = Zgoodss.z_peak 
zspec = Zgoodss.z_spec
spec = WHERE(zspec NE -1.) 
z[spec] = zspec[spec] 

; PRG1: z=1.67, m=1.8, b=-1.0, c=2.0
; PRG2: z=2.27, m=2.5, b=-0.8, c=99.0
; PRG3: z=2.14, m=1.8, b=-1.0, c=2.0

colordiagramHST, goodssred, goodssmid, goodssbluealt, 'F140W', 'F814W', 'F435W', z, 1.8, -1.0, 2.0, 1.67, 0.15, 'PRG1'
plotsave, '3DHSTcolordiagram_PRG1.png'
colordiagramHST, goodssred, goodssmid, goodssblue, 'F140W', 'F814W', 'F606W', z, 2.5, -0.8, 99.0, 2.27, 0.15, 'PRG2'
plotsave, '3DHSTcolordiagram_PRG2.png'
colordiagramHST, goodssred, goodssmid, goodssbluealt, 'F140W', 'F814W', 'F435W', z, 1.8, -1.0, 2.0, 2.14, 0.15, 'PRG3'
plotsave, '3DHSTcolordiagram_PRG3.png'



; NEED TO MAKE NEW CATALOGS WITH CUTS
;; make lists of the mags of sources in each catalog
;aegismags = list(aegisred, aegismid, aegisblue, ; red/mid/blue with cuts)
;cosmosmags = list()
;goodsnmags = list()
;goodssmags = list()
;udsmags = list()


END   ; HSTdensity procedure 









PRO fluxconverter, redflux, midflux, blueflux, redmag, midmag, bluemag 
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








; ALTER THIS TO FIT WITH 3DHST DATA
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





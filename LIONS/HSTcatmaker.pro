PRO HSTcatmaker, blueband=bluband
;-------------------------------------------------------------------------------------------------------------------------------------------------
; HSTcatmaker procedure
;   description 
; INPUTS: etc
;         etc   
; NOTES: 
;-------------------------------------------------------------------------------------------------------------------------------------------------
FORWARD_FUNCTION fluxtomag

; read in the data
goodss = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/GOODSS/Catalog/goodss_3dhst.v4.1.cat.FITS', 1)
Zgoodss = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/Z_GOODSS/goodss_3dhst.v4.1.5.zbest.fits', 1) 
;Zgoodss = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/GOODSS/Eazy/goodss_3dhst.v4.1.zout.FITS', 1) 

; define redshifts to use 
z = Zgoodss.z_best 
  ;z = Zgoodss.z_peak 
  ;zspec = Zgoodss.z_spec
  ;spec = WHERE(zspec NE -1.) 
  ;z[spec] = zspec[spec] 

; filter out everything that's flagged as having bad photometry
use = WHERE(goodss.use_phot EQ 1.) 
goodss = goodss[use]
z = z[use] 
; filter out galaxies that don't have the fluxes we need in red and middle bands
use = WHERE( (goodss.f_F140W GT 0.) AND (goodss.f_F814Wcand GT 0.) )
goodss = goodss[use]
z = z[use] 
   ;; filter out galaxies that have neither of the blue fluxes we want 
   ;use = WHERE( (goodss.f_F606W GT 0.) OR (goodss.f_F435W GT 0.) )
   ;goodss = goodss[use]
   ;z = z[use] 

IF (blueband EQ 'F606W') THEN BEGIN
  use = WHERE(goodss.f_F606W GT 0.)
  goodss = goodss[use]
  z = z[use] 
ENDIF ELSE IF (blueband EQ 'F435W') THEN BEGIN
  use = WHERE(goodss.f_F435W GT 0.)
  goodss = goodss[use]
  z = z[use] 
ENDIF ELSE BEGIN
  print, 'Specified blue band is not available. Using all galaxies with either F606W or F435W.'
  use = WHERE( (goodss.f_F606W GT 0.) OR (goodss.f_F435W GT 0.) )
  goodss = goodss[use]
  z = z[use] 
ENDELSE



; now convert fluxes to magnitudes
fluxconverter, goodss.f_F140W, goodss.f_F814Wcand, goodss.f_F606W, goodss.f_F435W, goodssred, goodssmid, goodssblue, goodssbluealt 

;goodss = goodss[use]
;goodssred = goodssred[use]
;goodssmid = goodssmid[use]
;goodssblue = goodssblue[use]
;goodssbluealt = goodssbluealt[use]
;z = z[use]

; finally, make the catalog containing only relevant information
IF (blueband EQ 'F606W') THEN BEGIN
  header = ['ID', 'ra', 'dec', 'F140W','F814W','F606W','z'] 
  write_csv, 'goodssF606W.csv', goodss.id, goodss.ra, goodss.dec, goodssred, goodssmid, goodssblue, z, header=header
ENDIF ELSE IF (blueband EQ 'F435W') THEN BEGIN
  header = ['ID', 'ra', 'dec', 'F140W','F814W','F435W','z'] 
  write_csv, 'goodssF435W.csv', goodss.id, goodss.ra, goodss.dec, goodssred, goodssmid, goodssbluealt, z, header=header
ENDIF ELSE BEGIN
  header = ['ID', 'ra', 'dec', 'F140W','F814W','F606W', 'F435W','z'] 
  write_csv, 'goodss.csv', goodss.id, goodss.ra, goodss.dec, goodssred, goodssmid, goodssblue, goodssbluealt, z, header=header
ENDELSE


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
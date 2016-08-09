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
aegis = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/AEGIS/Catalog/aegis_3dhst.v4.1.cat.FITS',)
cosmos = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/COSMOS/Catalog/cosmos_3dhst.v4.1.cat.FITS',)
goodsn = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/GOODSN/Catalog/goodsn_3dhst.v4.1.cat.FITS',)
goodss = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/GOODSS/Catalog/goodss_3dhst.v4.1.cat.FITS',)
uds = mrdfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/3DHST/UDS/Catalog/uds_3dhst.v4.2.cat.FITS',)

; bands are F140W, F814W, and F606W or F475W

fluxconverter, aegis.f_F140W, aegis.f_F814W, aegis.f_F606W, aegisred, aegismid, aegisblue
fluxconverter, cosmos.f_F140W, cosmos.f_F814W, cosmos.f_F606W, cosmosred, cosmosmid, cosmosblue
fluxconverter, goodsn.f_F140W, goodsn.f_F814W, goodsn.f_F606W, goodsnred, goodsnmid, goodsnblue
fluxconverter, goodss.f_F140W, goodss.f_F814W, goodss.f_F606W, goodssred, goodssmid, goodssblue
fluxconverter, uds.f_F140W, uds.f_F814W, uds.f_F606W, udsred, udsmid, udsblue

; NEED TO MAKE NEW CATALOGS WITH CUTS


; make lists of the mags of sources in each catalog
aegismags = list(aegisred, aegismid, aegisblue, ; red/mid/blue with cuts)
cosmosmags = list()
goodsnmags = list()
goodssmags = list()
udsmags = list()


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
  magAB = 25.0 - 2.5*log10(flux) 
  return, magAB
END 







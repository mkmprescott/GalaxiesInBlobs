PRO experiment 
  lista = list()
  listb = list()
  x = [1,2,3]
  y = [1,2,3,4] 
  lista.Add, x
  listb.Add, y
  array = [lista,listb]
  help, array
  print, array
  print, array[0]
  print, array[1]
END 








PRO galbinning, brightestmag, dimmestmag, binsize, ra, dec, radius, datacube, bigmaglist
;-------------------------------------------------------------------------------------------------------------------------------------------------
; galbinning procedure
;   description 
; INPUTS: etc
;         etc   
; NOTE: Uses get_galaxies_binned function. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION get_galaxies_binned
  bigmaglist = list()
  FOR i = 0, 5 DO BEGIN 
    maglist=list()
    FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
      galaxies_binned = get_galaxies_binned(ra, dec, radius, datacube[i], mag, mag+binsize)
      maglist.Add, galaxies_binned
    ENDFOR 
    bigmaglist.Add, maglist 
  ENDFOR   
END 









PRO get_ngal_binned, nbins, bigmaglist, ngal_binned, ngalerr_binned, ap=ap
;-------------------------------------------------------------------------------------------------------------------------------------------------
; get_ngal_binned procedure
;   description 
; INPUTS: etc
;         etc   
; NOTE: Uses Poisson_error function. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
FORWARD_FUNCTION Poisson_error 
  FOR band= 0, 5 DO BEGIN 
    maglist = bigmaglist[band]
    FOR i=0, nbins-1. DO BEGIN
      bin = maglist[i]
      IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
      IF (n_elements(ap) EQ 0) THEN BEGIN 
        ngal_binned[i,band] = ngal 
        ngalerr_binned[i,band] = Poisson_error(ngal) 
      ENDIF ELSE IF (n_elements(ap) NE 0) THEN BEGIN 
        ngal_binned[i,band,ap] = ngal 
        ngalerr_binned[i,band,ap] = Poisson_error(ngal) 
      ENDIF 
    ENDFOR   ; all the bines 
  ENDFOR   ; all the bands
END 
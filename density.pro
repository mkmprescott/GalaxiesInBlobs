PRO density   ;, seed, nApertures

FORWARD_FUNCTION get_galaxies

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
blobdensity = float(blob_ngal) / (!PI*10.^2.)   ; the denominator is the area of the aperture around the blob, which has 10" radius  

print, blobdensity, " galaxies per square arcsecond in blob region"  

; !!!!!!! CODE PRINTS:     0.0159155 galaxies per square arcsecond in blob region





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
  distance = sqrt(  ( (ra - data.FIELD15) * cos(dec*!PI/180.) )^2. + (dec - data.FIELD16)^2. )
  galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)
  RETURN, galaxies_in_aperture
END

PRO density, seed, nApertures, blobra, blobdec  ; the latter two are just temporary 

FORWARD_FUNCTION get_galaxies

data = read_csv("prg1_bw.cat", n_table_header=14)
; will later modify this to be more automated

; !!!!! include RA & dec in output catalog (currently absent)

; !!!!! ALSO: will have to convert RA to degrees!

blobra = (blob ra)    ; we know these coords already
blobdec = (blob dec)
; I think I'm going to make a data file that contains the RA and dec of all our blobs 
;   and use that for this part - later

; get blob density: number of galaxies inside aperture centered on blob (radius of say, 10 arcsec or so)
blob_galaxies=get_galaxies(blobra, blobdec, 10 arcsec, data)
blob_ngal = n_elements(blob_galaxies)
blobdensity = blob_ngal / (!PI*10.^2.)   ; the denominator is the area of the aperture around the blob, which has 10" radius  

aperture_ra_list = list()     ; this will contain RAs of apertures
aperture_dec_list = list()    ; this will contain declinations of apertures
aperture_radius_list = list() ; this will contain aperture radii 
aperture_ngal_list = list()   ; this will contain # galaxies in each aperture

FOR i=0, nApertures DO BEGIN
;  place an aperture randomly and add its properties to the lists of aperture properties
  aperture_radius = 10.*randomu(seed)
  aperture_radius_list.Add, aperture_radius
  aperture_ra = (minimum ra of image) + (max ra - min ra)*randomu(seed)
  aperture_ra_list.Add, aperture_ra
  aperture_dec = (minimum dec of image) + (max dec - min dec)*randomu(seed)
  ; figure out the number of galaxies in the aperture 
  galaxies = get_galaxies(aperture_ra, aperture_dec, aperture_radius, data)
  ngal = n_elements(galaxies)
  ; NOTE: ALSO OUTPUT THE WHOLE LIST OF GALAXIES!!! FOR MAG STUFF TO BE ADDED SOON
  ; append # of galaxies in the aperture to list of those
  aperture_ngal_list.Add, ngal
ENDFOR

; turn all the lists into arrays to work with them now that they're finished being constructed
aperture_ra_array = aperture_ra_list.ToArray(Type=5)
aperture_dec_array = aperture_dec_list.ToArray(Type=5)
aperture_radius_array = aperture_radius_list.ToArray(Type=5)
aperture_ngal_array = aperture_ngal_list.ToArray(Type=5)

; compute density PER MAG BIN!! FIX THIS!!
densityarray = aperture_ngal_array / ( !PI*(aperture_radius_array)^2. )    ; array of the densities in each aperture
field_density = total(densityarray) / double(n_elements(densityarray))     ; this is the field density for the whole image

; make histogram of densities (number of apertures in each density bin on y axis) & oplot blob density
; make plot (galaxy magnitude and density)

; make sure to add IF statement accounting for edges of picture
; think about rejecting things that are on the blob

END





FUNCTION get_galaxies, ra, dec, radius, data
  distance = sqrt(  ( (ra-data.ra) * cos(dec) )^2. + (dec-data.dec)^2. )
  galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)
  RETURN, galaxies_in_aperture
END

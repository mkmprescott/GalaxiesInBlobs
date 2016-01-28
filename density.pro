PRO density 

FORWARD_FUNCTION get_galaxies

data = read_csv("prg1_bw.cat", n_table_header=14)

; include RA & dec in output catalog (currently absent)

blob = (blob coords)  ; we know these coords already
; blobra and blobdec
; Should I make this something the user can input? Or make it take an input file?
;  since there's more than one blob
; I think I'm going to make a data file that contains the RA and dec of all our blobs 
;   and use that for this part

blobdensity = number of galaxies inside aperture centered on blob (radius of say, 10 arcsec or so)
; use get_galaxies for this, yeah?
blobdensity=get_galaxies(blobra, blobdec, 10 arcsec, data)
; The "data" in this is the same as the read_csv data above, right?  

aperture_list = list(ra, dec, aperture radius, # galaxies in aperture)

FOR i=0, nApertures DO BEGIN
  place an aperture randomly
  append ra, dec, radius to list
  galaxies = get_galaxies(this aperture)
  # = n_elements(galaxies)
  append # to list
ENDFOR

compute density
; does z factor into this? 

make histogram of densities (number of apertures in each density bin on y axis) & oplot blob density
make plot (galaxy magnitude and density)

; make sure to add IF statement accounting for edges of picture
; think about rejecting things that are on the blob

END





FUNCTION get_galaxies, ra, dec, radius, data
  distance = sqrt(  ( (ra-data.ra) * cos(dec) )^2. + (dec-data.dec)^2. )
  galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)
  RETURN, galaxies_in_aperture
END

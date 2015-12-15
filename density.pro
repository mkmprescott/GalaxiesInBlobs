PRO density

data = read_csv("prg1_bw.cat", n_table_header=14)

; include RA & dec in output catalog (currently absent)

blob = (blob coords)

blobdensity = number of galaxies inside aperture centered on blob (radius of say, 10 arcsec or so)

aperture_list = list(ra, dec, aperture radius, # galaxies in aperture)

FOR i=0, number_of_apertures DO BEGIN
  place an aperture randomly
  append ra, dec, radius to list
  galaxies = get_galaxies(this aperture)
  # = n_elements(galaxies)
  append # to list
ENDFOR

compute density

make histogram of densities (number of apertures in each density bin on y axis) & oplot blob density
make plot (galaxy magnitude and density)

; make sure to add IF statement accounting for edges of picture
; think about rejecting things that are on the blob

END





FUNCTION get_galaxies, ra, dec, radius, data
  distance = sqrt(  ( (ra-data.ra) * cos(dec) )^2. + (dec-data.dec)^2. )
  galaxies_in_aperture = WHERE (distance < radius, n_galaxies)
  RETURN, galaxies_in_aperture
END

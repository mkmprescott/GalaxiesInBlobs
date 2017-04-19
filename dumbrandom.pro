PRO dumbrandom, n_apertures
  ap_radius_pix = 10./0.06
  npix_x = 23000.
  npix_y = 20000.

  FOR ap=1, n_apertures DO BEGIN

    print, 'aperture number ', ap
    apfinderdumb, npix_x, npix_y, ap_radius_pix, counter, ap_x, ap_y
    print, ' '
    print, ap_x, ap_y
    print, ' '
    print, ' '
    print, ' '

  ENDFOR

END




PRO apfinderdumb, npix_x, npix_y, ap_radius_pix, counter, ap_x, ap_y

    counter = 0.
    WHILE (counter LT 2) DO BEGIN
       ap_x = floor(ap_radius_pix + (double(npix_x)-2.*ap_radius_pix)*randomu(seed))
       ap_y = floor(ap_radius_pix + (double(npix_y)-2.*ap_radius_pix)*randomu(seed)) 
       print, ap_x, ap_y
       counter++
    ENDWHILE

END
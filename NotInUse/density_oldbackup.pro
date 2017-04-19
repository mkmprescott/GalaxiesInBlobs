;                                                           /   \       
; _                                                 )      ((   ))      (
;(@)                                               /|\      ))_((      /|\                                                 _ 
;|-|                                              / | \   <_`\|/`_>   / | \                                               (@)
;| |---------------------------------------------/--|-voV--\<>|<>/--Vov-|--\----------------------------------------------|-|
;|-|                                                  '^`   (> <)   '^`                                                   | |
;| |                                                        `\Y/'                                                         |-|
;|-| DENSITY.PRO                                                                                                          | |
;| |   This code determines the density of galaxies within Lya blobs and surrounding fields                               |-|
;|-|     to determine whether or not Lya blobs tend to be located in overdense regions.                                   | |
;| |                                                                                                                      |-|
;|-| INPUT:                                                                                                               | |
;| |   blobs.csv - a comma-separated-variable file containing information about each Lya blob as follows:                 |-|
;|-|       blob name (PRG1, PRG2, etc)                                                                                    | |
;| |       RA of blob in decimal degrees                                                                                  |-|
;|-|       declination of blob in decimal degrees                                                                         | |
;| |       name of a SExtractor catalog (.cat file) run on a blue-band image of the blob (F606W or F475W) containing      |-|
;|-|          all sources detected in a stacked image of the blob                                                         | |
;| |       name of a SExtractor catalog (.cat file) run on a mid-band image of the blob (F814W) containing all sources    |-|
;|-|          detected in a stacked image of the blob                                                                     | |
;| |       name of a SExtractor catalog (.cat file) run on a red-band image of the blob (F140W) containing all sources    |-|
;|-|          detected in a stacked image of the blob                                                                     | |
;| |       name of a stack of all the images of the blob (a fits file)                                                    |-|
;|-|       name of a bad pixel mask corresponding to the stack (fits file)                                                | |
;| |       the redshift of the blob                                                                                       |-|
;|-|       the slope of a line marking a simple color cut for determining blob membership (unitless)                      | |
;| |       the y-intercept of the simple color cut line (unitless)                                                        |-|
;|-|       a "cap" value for color cuts, above which all galaxy colors are consistent with the blob redshift (unitless)   | |
;| |   ap_radius_arcsec - the radius in arcseconds of an aperture which defines the size of the Lya blob                  |-|
;|-|   nApertures - a number of random apertures to be placed within the field to determine the average field density     | | 
;| |                                                                                                                      |-|
;|-| OUTPUT:                                                                                                              | |
;| |   First, the SExtractor catalogs are read in using rsex (see online documentation) and bad pixel mask                |-|
;|-|       is read in using mrdfits.                                                                                      | |
;| |   After this initial reading-in, the raw number of galaxies within an aperture of the specified radius placed        |-|
;|-|       at the blob coordinates is determined. A region file (for DS9) is created marking each blob galaxy. The        | |
;| |       number density (number of galaxies per square degree) is also calculated.                                      |-|
;|-|   Next, new catalogs are made which contain only sources with colors consistent with the blob redshift as            | |
;| |       determined by a straight-line cut:                                                                             |-|
;|-|       consistent (middle band - red band) >= b + m*(blue band - middle band)                                         | |
;| |       and the refined number of galaxies within the blob is determined. A second region file containing only         |-|
;|-|       post-color-cut blob galaxies is created.                                                                       | |
;| |     color-color plot is created (middle band - red band vs blue band - middle band) showing the blob and field       |-|
;|-|       galaxies, with the color cut line overlaid to indicate which sources are thrown out by the cuts. The user      | |
;| |       may save this plot if desired.                                                                                 |-|
;|-|   From this point onward, every step is run on both catalogs: before AND after color cuts.                           | |
;| |   The number of galaxies in the blob is then split up into magnitude bins (from 22 to 35 in steps of 0.5) and        |-|
;|-|       the Poisson error on the number of galaxies in each bin is calculated. The total number of galaxies (in        | |
;| |       all bins, added up) is returned as a check to ensure that the mag binning process doesn't change the           |-|
;|-|       total number of galaxies detected in the blob region. The number of galaxies per magnitude bin is then         | |
;| |       converted into a number density of galaxies in number per square degree per 0.5 magnitudes and errors are      |-|
;|-|       likewise converted from number to density.                                                                     | |
;| |   The user is then asked whether or not to proceed with the analysis. This is for two reasons: because the process   |-|
;|-|       of determining the average field density takes longer than other parts of the code, and because the user       | |
;| |       may only want to create color-color plots or find the total number or density of blob galaxies. If the user    |-|
;|-|       opts not to proceed, then the code moves on to the next blob.                                                  | |
;| |   If the user chooses to proceed, then the code samples the field to determine the average field density. This       |-|
;|-|       is done by placing a specified number (nApertures) of random apertures in the field around the blob and        | |
;| |       counting up the number (with error bars) of galaxies in each aperture in each magnitude bin, converting        |-|
;|-|       those numbers into densities (with error bars), and then averaging the densities of all the apertures for      | |
;| |       a final count of average field density in every magnitude bin with errors.                                     |-|
;|-|   When placing random apertures, the blob region is avoided, as are areas that contain only bad pixels. If the       | |
;| |       random aperture lands in one of these forbidden areas, it is replaced until it lands in an area outside        |-|
;|-|       the blob with at least one good pixel.                                                                         | |
;| |   The bad pixel map of each blob is used to determine the actual effective area of each random aperture. Bad         |-|
;|-|       pixels are masked and therefore don't contain data, so the fraction of good pixels in each aperture is         | |
;| |       calculated and multiplied by the total aperture area to give the effective aperture area, which is the         |-|
;|-|       area used in calculating the number density of galaxies within each random aperture.                           | |
;| |   Once the average field density has been determined, a plot of the number density of galaxies vs magnitude          |-|
;|-|       for both the blob and the field is made. The blob density is plotted as a solid line with error bars,          | |
;| |       while the field density is plotted as a dashed line surrounded by a gray swath indicating errors. The          |-| 
;|-|       user can choose whether or not to save this plot.                                                              | |
;| |   A statistical luminosity function is also made to test the accuracy of the color cuts. This is done by making      |-|
;|-|       histograms of the raw number of galaxies in the blob in each magnitude bin and the number of galaxies within   | |
;| |       the blob in each magnitude bin after color cuts, and then overlaying a plot of the expected number of          |-|
;|-|       galaxies in the blob as determined by subtracting the raw field density from the raw blob density and          | |
;| |       converting the result into a number of galaxies in each magnitude bin (by multiplying by the area of the       |-|
;|-|       aperture placed around the blob). The user may save this plot if desired.                                      | |
;| |                                                                                                                      |-|
;|-|                                                                                                                      | |
;| | FUNCTIONS USED:                                                                                                      |-|
;|-|   rsex - reads SExtractor catalogs into structures (see online documentation)                                        | |
;| |   get_galaxies - finds which galaxies in a catalog are located within a specified (location and radius) aperture     |-|
;|-|   get_galaxies_binned - does the same as get_galaxies, but sorts the returned galaxies into magnitude bins           | |
;| |   get_galaxies_noblob - finds which galaxies in a catalog are located OUTSIDE a specified aperture                   |-|
;|-|   Poisson_error - calculates the Poisson uncertainty of a number (noise = square root of signal)                     | |
;| |   colorsort - filters a catalog by color cuts such that only galaxies with colors above the cut are returned         |-|
;|-| Each of these functions is described in more detail in its own header.                                               | |
;| |                                                                                                                      |-|
;|-|                                                                                                                      | |
;|_|______________________________________________________________________________________________________________________|-|
;(@)                                          l   /\ /        ) )       \ /\   l                                        `\|_|
;                                             l /   V        / /         V   \ l                                          (@)
;                                             l/           _( (_              \I
;                                                          `\ /'
;                                                            V
;-----------------------------------------------------------------------------------------------------------------------------
PRO density

 ; INITIAL SETUP 

   ; notify IDL of all functions used in this procedure
      FORWARD_FUNCTION rsex
      FORWARD_FUNCTION get_galaxies
      FORWARD_FUNCTION get_galaxies_binned
      FORWARD_FUNCTION get_galaxies_noblob
      FORWARD_FUNCTION Poisson_error
      FORWARD_FUNCTION colorsort
      FORWARD_FUNCTION titlemaker 

   ; USER-DEFINED VARIABLES: 
      pixelscale = 0.06               ; arceseconds per pixel for our images (may eventually become part of blobs.csv) 
      ap_radius_arcsec = [10.,12.]    ; aperture radius in arcsec
      binsize=0.5                     ; size of magnitude bins 
      brightestmag = 22.              ; smallest (brightest) magnitude of galaxies we can reasonably see in our images 
      dimmestmag = 30.                ; biggest (faintest) magnitude of galaxies we can reasonably see in our images 
      nApertures = 10.                ; number of random apertures to use in field density calculation 
      saveanyplots = 'y'              ; string indicating that user should be asked whether or not to save plots 

   ; read in the file containing all blob information
      blobs = read_csv('blobs.csv', n_table_header=1)
   ; assign names to each field of the "blobs" structure
      blobname = blobs.FIELD01
      blobra = blobs.FIELD02
      blobdec = blobs.FIELD03
      blueband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/catalogs/'+blobs.FIELD04
      midband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/catalogs/'+blobs.FIELD05
      redband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/catalogs/'+blobs.FIELD06
      stack = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/stacks/'+blobs.FIELD07
      mask = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/stacks/'+blobs.FIELD08 
      z = blobs.FIELD09                 ; redshift
      m = blobs.FIELD10                 ; slope of color-cut line 
      b = blobs.FIELD11                 ; y-intercept of color-cut line 
      cap = blobs.FIELD12               ; mag above which all colors are consistent with redshift
      bluename = blobs.FIELD13          ; name of the bluest filter used to observe the blob region
      midname = blobs.FIELD14           ; name of the middle filter used to observe the blob region
      redname = blobs.FIELD15           ; name of the reddest filter used to observe the blob region
      nblobs = n_elements(blobname)     ; number of blobs in the blobs.csv file 

   ; define other important quantities  
      ; the radius/radii of the aperture(s) used to calculate density:  
         ap_radius = ap_radius_arcsec/3600.             ; aperture radius in decimal degrees
         ap_radius_pix = ap_radius_arcsec / pixelscale  ; aperture radius in pixels  
         n_radii = n_elements(ap_radius)                ; number of aperture sizes being tested
      ; mag binning stuff:
         nbins = (dimmestmag - brightestmag + binsize)/binsize     ; the number of magnitude bins as determined from our range and binsize 
         mags = brightestmag + binsize*findgen(nbins)              ; an array of every magnitude being used 
         magrange=[brightestmag,dimmestmag]                        ; a two-element array containing only the endpoints of our magnitude range 
      ; plotting symbol: make the psym=8 symbol a filled circle for our plots 
         plotsym, 0, 1, /fill   



 ; DENSITY ANALYSIS FOR EACH BLOB 

 FOR blob=0, nblobs-1 DO BEGIN         

   ; make structures for each SExtractor catalog (one per each of 3 bands for each blob) 
      datablue = rsex(blueband[blob])    ; either F475W (PRG1 and PRG3) or F606W (PRG2)
      datamid = rsex(midband[blob])      ; F814W
      datared = rsex(redband[blob])      ; F140W
   ; establish info about the images
      bpm = mrdfits(mask[blob],0,hdr)    ; this is the mask image for weeding out bad pixels in field 
      sizeinfo = size(bpm)               ; the dimensions of the mask image in pixels 
      npix_x = sizeinfo[1]               ; number of pixels in mask in the x direction
      npix_y = sizeinfo[2]               ; number of pixels in mask in the y direction 
   ; make arrays of info about aperture radii - will contain overdensity factor & error on overdensity factor at each magnitude
      radiusresults_red_nocuts = dblarr(n_radii, nbins, 2.)
      radiusresults_mid_nocuts = dblarr(n_radii, nbins, 2.)
      radiusresults_blue_nocuts = dblarr(n_radii, nbins, 2.)
      radiusresults_red_cuts = dblarr(n_radii, nbins, 2.)
      radiusresults_mid_cuts = dblarr(n_radii, nbins, 2.)
      radiusresults_blue_cuts = dblarr(n_radii, nbins, 2.) 


   ; now run analysis of density for different aperture radii 
   FOR r=0, n_radii-1. DO BEGIN

     ; First get raw blob density, ie, number of galaxies inside aperture centered on blob with no mag bins yet
     blob_galaxies_all=get_galaxies(blobra[blob], blobdec[blob], ap_radius[r], datared)
     blob_ngalraw = n_elements(blob_galaxies_all)
     ; now get density in number of galaxies per square arcsecond 
     blobdensityraw = float(blob_ngalraw) / (!dPI*ap_radius[r]^2.)    
     print, blobdensityraw, " galaxies per square degree in blob region for "+blobname[blob]  
     print, blob_ngalraw, " galaxies in the blob region for "+blobname[blob] 

     ; make new catalogs with color cuts
     xcolorall = datablue.MAG_ISO - datamid.MAG_ISO
     ycolorall = datamid.MAG_ISO - datared.MAG_ISO
     dataredcut = colorsort(datared, xcolorall, ycolorall, m[blob], b[blob], cap[blob]) 
     datamidcut = colorsort(datamid, xcolorall, ycolorall, m[blob], b[blob], cap[blob]) 
     databluecut = colorsort(datablue, xcolorall, ycolorall, m[blob], b[blob], cap[blob]) 
     ; get the number of galaxies in the aperture AFTER color cuts
     blob_galaxies_all_cuts=get_galaxies(blobra[blob], blobdec[blob], ap_radius[r], dataredcut)
     blob_galaxies_all_cuts_mid=get_galaxies(blobra[blob], blobdec[blob], ap_radius[r], datamidcut)
     blob_galaxies_all_cuts_blue=get_galaxies(blobra[blob], blobdec[blob], ap_radius[r], databluecut)

     ; make region files showing blob galaxies before and after color cuts  
     makereg, datared, blob_galaxies_all, blobname[blob], ap_radius_arcsec[r], 'n'            ; before cuts 
     makereg, dataredcut, blob_galaxies_all_cuts, blobname[blob], ap_radius_arcsec[r], 'y'    ; after cuts 

     ; make color-color diagram: 
     field_galaxies = get_galaxies_noblob(blobra[blob], blobdec[blob], ap_radius[r], datared)   ; get IDs of field galaxies 
       ; set up plot window 
         !P.MULTI = [0,1,1,0,0] 
         loadct, 0
         window, 0, retain=2, xsize=400, ysize=350
     colordiagram, blob_galaxies_all, field_galaxies, datared, datamid, datablue, redname[blob], midname[blob], bluename[blob], m[blob], b[blob], cap[blob], blobname[blob], ap_radius_arcsec[r] 
     namestring = titlemaker('colors', blobname[blob], radius=ap_radius_arcsec[r], filetype='.png')
     plotsave, saveanyplots, namestring

     ; MAG BINNING: order is red raw, mid raw, blue raw, red cuts, mid cuts, blue cuts 
     galbinning, brightestmag, dimmestmag, binsize, blobra[blob], blobdec[blob], ap_radius[r], datared, galaxy_maglist, blob_galaxies_binned 
     galbinning, brightestmag, dimmestmag, binsize, blobra[blob], blobdec[blob], ap_radius[r], datamid, galaxy_maglist_mid, blob_galaxies_binned_mid 
     galbinning, brightestmag, dimmestmag, binsize, blobra[blob], blobdec[blob], ap_radius[r], datablue, galaxy_maglist_blue, blob_galaxies_binned_blue 
     galbinning, brightestmag, dimmestmag, binsize, blobra[blob], blobdec[blob], ap_radius[r], dataredcut, galaxy_maglist_cuts, blob_galaxies_binned_cuts 
     galbinning, brightestmag, dimmestmag, binsize, blobra[blob], blobdec[blob], ap_radius[r], datamidcut, galaxy_maglist_cuts_mid, blob_galaxies_binned_cuts_mid 
     galbinning, brightestmag, dimmestmag, binsize, blobra[blob], blobdec[blob], ap_radius[r], databluecut, galaxy_maglist_cuts_blue, blob_galaxies_binned_cuts_blue 

     blob_ngal_binned = dblarr(nbins)
     blob_ngalerr_binned = dblarr(nbins) 
     blob_ngal_binned_mid = dblarr(nbins)
     blob_ngalerr_binned_mid = dblarr(nbins) 
     blob_ngal_binned_blue = dblarr(nbins)
     blob_ngalerr_binned_blue = dblarr(nbins) 
     blob_ngal_binned_cuts = dblarr(nbins)
     blob_ngalerr_binned_cuts = dblarr(nbins) 
     blob_ngal_binned_cuts_mid = dblarr(nbins)
     blob_ngalerr_binned_cuts_mid = dblarr(nbins) 
     blob_ngal_binned_cuts_blue = dblarr(nbins)
     blob_ngalerr_binned_cuts_blue = dblarr(nbins) 

     get_ngal_binned, nbins, galaxy_maglist, blob_ngal_binned, blob_ngalerr_binned
     get_ngal_binned, nbins, galaxy_maglist_mid, blob_ngal_binned_mid, blob_ngalerr_binned_mid
     get_ngal_binned, nbins, galaxy_maglist_blue, blob_ngal_binned_blue, blob_ngalerr_binned_blue
     get_ngal_binned, nbins, galaxy_maglist_cuts, blob_ngal_binned_cuts, blob_ngalerr_binned_cuts
     get_ngal_binned, nbins, galaxy_maglist_cuts_mid, blob_ngal_binned_cuts_mid, blob_ngalerr_binned_cuts_mid
     get_ngal_binned, nbins, galaxy_maglist_cuts_blue, blob_ngal_binned_cuts_blue, blob_ngalerr_binned_cuts_blue

     print, total(blob_ngal_binned), " total galaxies"   ; just a test to make sure it's working - this should match the previous no-bin number

     densitycalc, blob_ngal_binned, blob_ngalerr_binned, ap_radius[r], densities, density_errors
     densitycalc, blob_ngal_binned_mid, blob_ngalerr_binned_mid, ap_radius[r], densities_mid, density_errors_mid
     densitycalc, blob_ngal_binned_blue, blob_ngalerr_binned_blue, ap_radius[r], densities_blue, density_errors_blue
     densitycalc, blob_ngal_binned_cuts, blob_ngalerr_binned_cuts, ap_radius[r], densities_cuts, density_errors_cuts
     densitycalc, blob_ngal_binned_cuts_mid, blob_ngalerr_binned_cuts_mid, ap_radius[r], densities_cuts_mid, density_errors_cuts_mid
     densitycalc, blob_ngal_binned_cuts_blue, blob_ngalerr_binned_cuts_blue, ap_radius[r], densities_cuts_blue, density_errors_cuts_blue

;     magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datared, galaxy_maglist, blob_galaxies_binned,blob_ngal_binned, blob_ngalerr_binned, densities, density_errors
;     magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datamid, galaxy_maglist_mid, blob_galaxies_binned_mid, blob_ngal_binned_mid, blob_ngalerr_binned_mid, densities_mid, density_errors_mid
;     magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datablue, galaxy_maglist_blue, blob_galaxies_binned_blue, blob_ngal_binned_blue, blob_ngalerr_binned_blue, densities_blue, density_errors_blue
;     magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], dataredcut, galaxy_maglist_cuts, blob_galaxies_binned_cuts, blob_ngal_binned_cuts, blob_ngalerr_binned_cuts, densities_cuts, density_errors_cuts
;     magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], datamidcut, galaxy_maglist_cuts_mid, blob_galaxies_binned_cuts_mid, blob_ngal_binned_cuts_mid, blob_ngalerr_binned_cuts_mid, densities_cuts_mid, density_errors_cuts_mid
;     magbinning, brightestmag, dimmestmag, binsize, nbins, blobra[blob], blobdec[blob], ap_radius[r], databluecut, galaxy_maglist_cuts_blue, blob_galaxies_binned_cuts_blue, blob_ngal_binned_cuts_blue, blob_ngalerr_binned_cuts_blue, densities_cuts_blue, density_errors_cuts_blue





; make this stuff into a procedure too

     print, 'Now sampling the field for overdensity analysis.' 

     ; SAMPLING THE FIELD with random apertures 

     ; set up the arrays which will contain all the information about the random apertures placed
     ap_ras = dblarr(nApertures)     ; this will contain RAs of all apertures
     ap_decs = dblarr(nApertures)    ; this will contain declinations of all apertures
     ; red band 
     ap_ngal_binned = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
     ap_ngalerr_binned = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
     ap_ngal_binned_cuts = dblarr(nbins,nApertures)  
     ap_ngalerr_binned_cuts = dblarr(nbins,nApertures)   
     megalist = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 
     megalist_cuts = list() 
     ; middle band 
     ap_ngal_binned_mid = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
     ap_ngalerr_binned_mid = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
     ap_ngal_binned_cuts_mid = dblarr(nbins,nApertures)  
     ap_ngalerr_binned_cuts_mid = dblarr(nbins,nApertures)   
     megalist_mid = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 
     megalist_cuts_mid = list() 
     ; blue band 
     ap_ngal_binned_blue = dblarr(nbins,nApertures)   ; this will contain # galaxies in each aperture in each mag bin
     ap_ngalerr_binned_blue = dblarr(nbins,nApertures)   ; same as above but errors on # rather than just #
     ap_ngal_binned_cuts_blue = dblarr(nbins,nApertures)  
     ap_ngalerr_binned_cuts_blue = dblarr(nbins,nApertures)   
     megalist_blue = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin 
     megalist_cuts_blue = list() 
     ; make an aperture image and blank image into which to insert it 
     dist_ellipse, dim, [ 2.*ap_radius_pix[r], 2.*ap_radius_pix[r] ], ap_radius_pix[r], ap_radius_pix[r], 1.,0.
     blankim = ap_radius_pix[r]+fltarr(npix_x,npix_y)  
     ap_count = 0 ; keep track of how many apertures are made in total  

     ; now place random apertures and count galaxies inside them 
     TIC        ; time this to ensure that it's efficient 
     FOR ap=0, nApertures-1 DO BEGIN
       apfinder, npix_x, npix_y, ap_radius_pix[r], ap_radius[r], dim, blankim, bpm, hdr, blobra[blob], blobdec[blob], ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec 

       ; fraction of aperture that is good pixels:
       ap_area_weight = (float(ntotalpix) - float(nbadpix)) / float(ntotalpix) 
 
       ; convert pixels to sky coords and add these coords to arrays of aperture properties
       xyad, hdr,ap_x, ap_y, ap_ra, ap_dec
       ap_ras[ap] = ap_ra
       ap_decs[ap] = ap_dec 

       galbinning, brightestmag, dimmestmag, binsize, ap_ra, ap_dec, ap_radius[r], datared, ap_galaxy_maglist, ap_galaxies_binned
       galbinning, brightestmag, dimmestmag, binsize, ap_ra, ap_dec, ap_radius[r], datamid, ap_galaxy_maglist_mid, ap_galaxies_binned_mid
       galbinning, brightestmag, dimmestmag, binsize, ap_ra, ap_dec, ap_radius[r], datablue, ap_galaxy_maglist_blue, ap_galaxies_binned_blue
       galbinning, brightestmag, dimmestmag, binsize, ap_ra, ap_dec, ap_radius[r], dataredcut, ap_galaxy_maglist_cuts, ap_galaxies_binned_cuts
       galbinning, brightestmag, dimmestmag, binsize, ap_ra, ap_dec, ap_radius[r], datamidcut, ap_galaxy_maglist_cuts_mid, ap_galaxies_binned_cuts_mid
       galbinning, brightestmag, dimmestmag, binsize, ap_ra, ap_dec, ap_radius[r], databluecut, ap_galaxy_maglist_cuts_blue, ap_galaxies_binned_cuts_blue

       ; put the list of ID numbers per mag bin in this aperture into the mega list for all apertures:
       megalist.Add, ap_galaxy_maglist   
       megalist_cuts.Add, ap_galaxy_maglist_cuts
       megalist_mid.Add, ap_galaxy_maglist_mid   
       megalist_cuts_mid.Add, ap_galaxy_maglist_cuts_mid
       megalist_blue.Add, ap_galaxy_maglist_blue   
       megalist_cuts_blue.Add, ap_galaxy_maglist_cuts_blue

       get_ngal_binned, nbins, ap_galaxy_maglist, ap_ngal_binned, ap_ngalerr_binned, ap=ap
       get_ngal_binned, nbins, ap_galaxy_maglist_mid, ap_ngal_binned_mid, ap_ngalerr_binned_mid, ap=ap
       get_ngal_binned, nbins, ap_galaxy_maglist_blue, ap_ngal_binned_blue, ap_ngalerr_binned_blue, ap=ap
       get_ngal_binned, nbins, ap_galaxy_maglist_cuts, ap_ngal_binned_cuts, ap_ngalerr_binned_cuts, ap=ap
       get_ngal_binned, nbins, ap_galaxy_maglist_cuts_mid, ap_ngal_binned_cuts_mid, ap_ngalerr_binned_cuts_mid, ap=ap
       get_ngal_binned, nbins, ap_galaxy_maglist_cuts_blue, ap_ngal_binned_cuts_blue, ap_ngalerr_binned_cuts_blue, ap=ap

     ENDFOR    ; all the random apertures 
     print, ap_count ; How many total apertures were generated? 
     TOC             ; How long did this analysis take? 



    ; calculate densities in each bin in each aperture and errors, then average the densities of all the apertures in each bin to get average field density in each bin: 
    ;

    densitycalc, ap_ngal_binned, ap_ngalerr_binned, ap_radius[r], ap_densities, ap_density_errors, weight=ap_area_weight
    densitycalc, ap_ngal_binned_mid, ap_ngalerr_binned_mid, ap_radius[r], ap_densities_mid, ap_density_errors_mid, weight=ap_area_weight
    densitycalc, ap_ngal_binned_blue, ap_ngalerr_binned_blue, ap_radius[r], ap_densities_blue, ap_density_errors_blue, weight=ap_area_weight
    densitycalc, ap_ngal_binned_cuts, ap_ngalerr_binned_cuts, ap_radius[r], ap_densities_cuts, ap_density_errors_cuts, weight=ap_area_weight
    densitycalc, ap_ngal_binned_cuts_mid, ap_ngalerr_binned_cuts_mid, ap_radius[r], ap_densities_cuts_mid, ap_density_errors_cuts_mid, weight=ap_area_weight
    densitycalc, ap_ngal_binned_cuts_blue, ap_ngalerr_binned_cuts_blue, ap_radius[r], ap_densities_cuts_blue, ap_density_errors_cuts_blue, weight=ap_area_weight

    ; RED BAND, NO CUTS 
;    ap_densities = double(ap_ngal_binned)/(!dPI*ap_radius[r]^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
;    ap_density_errors = double(ap_ngalerr_binned)/(!dPI*ap_radius[r]^2.*ap_area_weight)   ; same as above except errors 
    field_densities = total(ap_densities,2)/double(nApertures)             ; field density in each bin 
    field_density_errors = total(ap_density_errors,2)/double(nApertures)   ; errors in field density in each bin 
    radiusresults_red_nocuts[r,*,0] = densities/field_densities            ; overdensity factor
    radiusresults_red_nocuts[r,*,1] = density_errors/field_densities       ; error on overdensity factor
    ;
    ; MIDDLE BAND, NO CUTS 
;    ap_densities_mid = double(ap_ngal_binned_mid)/(!dPI*ap_radius[r]^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
;    ap_density_errors_mid = double(ap_ngalerr_binned_mid)/(!dPI*ap_radius[r]^2.*ap_area_weight)   ; same as above except errors 
    field_densities_mid = total(ap_densities_mid,2)/double(nApertures)             ; field density in each bin 
    field_density_errors_mid = total(ap_density_errors_mid,2)/double(nApertures)   ; errors in field density in each bin 
    radiusresults_mid_nocuts[r,*,0] = densities_mid/field_densities_mid            ; overdensity factor
    radiusresults_mid_nocuts[r,*,1] = density_errors_mid/field_densities_mid       ; error on overdensity factor
    ;
    ; BLUE BAND, NO CUTS
;    ap_densities_blue = double(ap_ngal_binned_blue)/(!dPI*ap_radius[r]^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
;    ap_density_errors_blue = double(ap_ngalerr_binned_blue)/(!dPI*ap_radius[r]^2.*ap_area_weight)   ; same as above except errors 
    field_densities_blue = total(ap_densities_blue,2)/double(nApertures)             ; field density in each bin 
    field_density_errors_blue = total(ap_density_errors_blue,2)/double(nApertures)   ; errors in field density in each bin 
    radiusresults_blue_nocuts[r,*,0] = densities_blue/field_densities_blue           ; overdensity factor
    radiusresults_blue_nocuts[r,*,1] = density_errors_blue/field_densities_blue      ; error on overdensity factor
    ;
    ; RED BAND WITH CUTS
;    ap_densities_cuts = double(ap_ngal_binned_cuts)/(!dPI*ap_radius[r]^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
;    ap_density_errors_cuts = double(ap_ngalerr_binned_cuts)/(!dPI*ap_radius[r]^2.*ap_area_weight)   ; same as above except errors 
    field_densities_cuts = total(ap_densities_cuts,2)/double(nApertures)             ; field density in each bin 
    field_density_errors_cuts = total(ap_density_errors_cuts,2)/double(nApertures)   ; errors in field density in each bin
    radiusresults_red_cuts[r,*,0] = densities_cuts/field_densities_cuts              ; overdensity factor
    radiusresults_red_cuts[r,*,1] = density_errors_cuts/field_densities_cuts         ; error on overdensity factor
    ;
    ; MIDDLE BAND WITH CUTS
;    ap_densities_cuts_mid = double(ap_ngal_binned_cuts_mid)/(!dPI*ap_radius[r]^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
;    ap_density_errors_cuts_mid = double(ap_ngalerr_binned_cuts_mid)/(!dPI*ap_radius[r]^2.*ap_area_weight)   ; same as above except errors 
    field_densities_cuts_mid = total(ap_densities_cuts_mid,2)/double(nApertures)             ; field density in each bin 
    field_density_errors_cuts_mid = total(ap_density_errors_cuts_mid,2)/double(nApertures)   ; errors in field density in each bin 
    radiusresults_mid_cuts[r,*,0] = densities_cuts_mid/field_densities_cuts_mid              ; overdensity factor
    radiusresults_mid_cuts[r,*,1] = density_errors_cuts_mid/field_densities_cuts_mid         ; error on overdensity factor
    ;
    ; BLUE BAND WITH CUTS
;    ap_densities_cuts_blue = double(ap_ngal_binned_cuts_blue)/(!dPI*ap_radius[r]^2.*ap_area_weight)           ; array of the densities in each aperture in each mag bin
;    ap_density_errors_cuts_blue = double(ap_ngalerr_binned_cuts_blue)/(!dPI*ap_radius[r]^2.*ap_area_weight)   ; same as above except errors 
    field_densities_cuts_blue = total(ap_densities_cuts_blue,2)/double(nApertures)             ; field density in each bin 
    field_density_errors_cuts_blue = total(ap_density_errors_cuts_blue,2)/double(nApertures)   ; errors in field density in each bin 
    radiusresults_blue_cuts[r,*,0] = densities_cuts_blue/field_densities_cuts_blue             ; overdensity factor
    radiusresults_blue_cuts[r,*,1] = density_errors_cuts_blue/field_densities_cuts_blue        ; error on overdensity factor

; end of stuff for new procedure 


    ; MAIN OVERDENSITY PLOTS:
    ;
    !P.MULTI = [0,2,3,0,1]
    window, 1, retain=2, xsize=950, ysize=1200
    ;
    ; RED BAND, NO CUTS
    densityplot, nbins, mags, field_densities, field_density_errors, densities, density_errors, blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob], 'n'
      ; namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=redname[blob], cuts='n', filetype='.png')
      ; plotsave, saveanyplots, namestring
      ;
    ; MIDDLE BAND, NO CUTS 
    densityplot, nbins, mags, field_densities_mid, field_density_errors_mid, densities_mid, density_errors_mid, blobname[blob], nApertures, ap_radius_arcsec[r], midname[blob], 'n'
      ; namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=midname[blob], cuts='n', filetype='.png')
      ; plotsave, saveanyplots, namestring
      ;
    ; BLUE BAND, NO CUTS 
    densityplot, nbins, mags, field_densities_blue, field_density_errors_blue, densities_blue, density_errors_blue, blobname[blob], nApertures, ap_radius_arcsec[r], bluename[blob], 'n'
      ; namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=bluename[blob], cuts='n', filetype='.png')
      ; plotsave, saveanyplots, namestring
      ;
    ; RED BAND WITH CUTS
    densityplot, nbins, mags, field_densities_cuts, field_density_errors_cuts, densities_cuts, density_errors_cuts, blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob], 'y'
      ; namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=redname[blob], cuts='y', filetype='.png')
      ; plotsave, saveanyplots, namestring
      ;
    ; MIDDLE BAND WITH CUTS 
    densityplot, nbins, mags, field_densities_cuts_mid, field_density_errors_cuts_mid, densities_cuts_mid, density_errors_cuts_mid, blobname[blob], nApertures, ap_radius_arcsec[r], midname[blob], 'y'
      ; namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=midname[blob], cuts='y', filetype='.png')
      ; plotsave, saveanyplots, namestring
      ;
    ; BLUE BAND WITH CUTS 
    densityplot, nbins, mags, field_densities_cuts_blue, field_density_errors_cuts_blue, densities_cuts_blue, density_errors_cuts_blue, blobname[blob], nApertures, ap_radius_arcsec[r], bluename[blob], 'y'
      ; namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=bluename[blob], cuts='y', filetype='.png')
      ; plotsave, saveanyplots, namestring
      ; 
    namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], filetype='.png')
    plotsave, saveanyplots, namestring



    ; STATISTICAL LUMINOSITY FUNCTIONS 
    ; combine everything from before to make these - see statlum procedure header 
    !P.MULTI = [0,1,3,0,1]
    window, 1, retain=2, xsize=650, ysize=1200
    ;
    ; RED:
    statlum, mags, binsize, blob_ngal_binned, blob_ngal_binned_cuts, densities, field_densities, ap_radius[r], blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob]  
      ; namestring = titlemaker('statlum', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=redname[blob], filetype='.png') 
      ; plotsave, saveanyplots, namestring
      ;
    ; MIDDLE: 
    statlum, mags, binsize, blob_ngal_binned_mid, blob_ngal_binned_cuts_mid, densities_mid, field_densities_mid, ap_radius[r], blobname[blob], nApertures, ap_radius_arcsec[r], midname[blob]  
      ; namestring = titlemaker('statlum', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=midname[blob], filetype='.png') 
      ; plotsave, saveanyplots, namestring
      ;
    ; BLUE: 
    statlum, mags, binsize, blob_ngal_binned_blue, blob_ngal_binned_cuts_blue, densities_blue, field_densities_blue, ap_radius[r], blobname[blob], nApertures, ap_radius_arcsec[r], bluename[blob] 
      ; namestring = titlemaker('statlum', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], band=bluename[blob], filetype='.png') 
      ; plotsave, saveanyplots, namestring
      ;
    namestring = titlemaker('statlum', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], filetype='.png') 
    plotsave, saveanyplots, namestring



   ; ADD EVERYTHING BELOW THIS TO THE HEADER!!!

     ; How likely is it to find overdensities in this field anyway? 
     ap_numbers_corrected = ap_densities*!dPI*ap_radius[r]^2.
     ap_numbers_corrected_cuts = ap_densities_cuts*!dPI*ap_radius[r]^2.
     ap_totals = total(ap_numbers_corrected,1)    ; total number of galaxies in each aperture regardless of magnitude
     ap_totals_cuts = total(ap_numbers_corrected_cuts,1)    ; total number of galaxies in each aperture regardless of magnitude

     !P.MULTI = [0,2,1,0,1]
     window, 1, retain=2, xsize=1200, ysize=400
     ; compare total numbers (regardless of magnitude) to blob, no cuts
     makehist, ap_totals, total(blob_ngal_binned), blobname[blob], nApertures, ap_radius_arcsec[r], 'n'  
     ; namestring = titlemaker('aphist', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], cuts='n', filetype='.png)
     ; plotsave, saveanyplots, namestring
     ;
     ; compare total numbers (regardless of magnitude) to blob WITH cuts
     makehist, ap_totals_cuts, total(blob_ngal_binned_cuts), blobname[blob], nApertures, ap_radius_arcsec[r], 'y'  
     ; namestring = titlemaker('aphist', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], cuts='y', filetype='.png)
     ; plotsave, saveanyplots, namestring
     ; 
     namestring = titlemaker('aphist', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], filetype='.png')
     plotsave, saveanyplots, namestring




; ; fix this shit
;  IF (saveanyplots EQ 'y') THEN BEGIN 
;    histquestion = ''
;    READ, histquestion, PROMPT='Make histograms? (y/n)'
;    IF (histquestion EQ 'y') THEN BEGIN 
;      ; now look within each individual magnitude bin 
;      FOR i=0, nbins-1 DO BEGIN
;        n_in_aps = fltarr(nApertures)
;        n_in_aps[*] = ap_numbers_corrected[i,*]
;        checker = WHERE(n_in_aps NE 0, /null)
;        IF (checker NE !null) THEN BEGIN
;          window, 0, retain=2, xsize=1200, ysize=1000
;          title = 'Galaxy Number Density in Apertures for Galaxies between Magnitudes ' + string(mags[i], format='(F4.1)') + ' and ' + string(mags[i]+binsize, format='(F4.1)') + ' for ' + blobname[blob] + ', no cuts'
;          plothist, n_in_aps, background=255, color=0, xtitle='n galaxies in aperture', ytitle='number of apertures', title=title, axiscolor=0, bin=1, xrange=[0,50], /xstyle
;          blobgaltext = 'blob: ' + string(fix(blob_ngal_binned[i])) + ' galaxies in this bin'
;          xyouts, 30., nApertures/10., blobgaltext, color=0, charsize=1.5
;          saveplot = ''
;          READ, saveplot, PROMPT='Save this plot? (y/n)'
;          IF (saveplot EQ 'y') THEN BEGIN 
;           namestring = string(blobname[blob]) + '_' + apsnamesave + '_' + ap_radius_strings[r] + string(mags[i], format='(F4.1)') + '_aphist.png'
;           write_png, namestring, tvrd(/true) 
;          ENDIF 
;          wdelete, 0
;        ENDIF ELSE print, "no apertures contain any galaxies with magnitudes between ", string(mags[i], format='(F4.1)'), " and ", string(mags[i]+binsize, format='(F4.1)'), ' without cuts'
;        ; now do the same for the cuts 
;        n_in_aps_cuts = fltarr(nApertures)
;        n_in_aps_cuts[*] = ap_numbers_corrected_cuts[i,*]
;        checker = WHERE(n_in_aps_cuts NE 0, /null)
;        IF (checker NE !null) THEN BEGIN
;          window, 0, retain=2, xsize=1200, ysize=1000
;          title = 'Galaxy Number Density in Apertures for Galaxies between Magnitudes ' + string(mags[i], format='(F4.1)') + ' and ' + string(mags[i]+binsize, format='(F4.1)') + ' for ' + blobname[blob] + ', with cuts'
;          plothist, n_in_aps_cuts, background=255, color=0, xtitle='n galaxies in aperture', ytitle='number of apertures', title=title, axiscolor=0, bin=1, xrange=[0,50], /xstyle
;          blobgaltext = 'blob: ' + string(fix(blob_ngal_binned_cuts[i])) + ' galaxies in this bin'
;          xyouts, 30., nApertures/10., blobgaltext, color=0, charsize=1.5
;          saveplot = ''
;          READ, saveplot, PROMPT='Save this plot? (y/n)'
;          IF (saveplot EQ 'y') THEN BEGIN 
;           namestring = string(blobname[blob]) + '_' + apsnamesave + '_' + ap_radius_strings[r] + string(mags[i], format='(F4.1)') + '_aphist_cuts.png'
;           write_png, namestring, tvrd(/true) 
;          ENDIF 
;          wdelete, 0
;        ENDIF ELSE print, "no apertures contain any galaxies with magnitudes between ", string(mags[i], format='(F4.1)'), " and ", string(mags[i]+binsize, format='(F4.1)'), ' with cuts'
;      ENDFOR
;    ENDIF   ; user wanted to make histograms 
;  ENDIF   ; user wants to save at least some plots
;; end of shit that needs fixing



   ENDFOR   ; all aperture sizes


   ; if there were multiple aperture radii being tested, analyze them with plots:  
   IF (n_radii GT 1.) THEN BEGIN 
     ;
     ; now make plots of the results of the aperture size analysis: overdensity factor vs magnitude for each band for each radius size 
     loadct, 39    ; new color table 
     !P.MULTI = [0,2,3,0,1]
     window, 1, retain=2, xsize=950, ysize=1200
     ; RED BAND, NO CUTS
     Rmultirad, ap_radius_arcsec, mags, radiusresults_red_nocuts, blobname[blob], nApertures, redname[blob], 'n'
        ; namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, band=redname[blob], cuts='n', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; MIDDLE BAND, NO CUTS
     Rmultirad, ap_radius_arcsec, mags, radiusresults_mid_nocuts, blobname[blob], nApertures, midname[blob], 'n' 
        ; namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, band=midname[blob], cuts='n', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; BLUE BAND, NO CUTS
     Rmultirad, ap_radius_arcsec, mags, radiusresults_blue_nocuts, blobname[blob], nApertures, bluename[blob], 'n' 
        ; namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, band=bluename[blob], cuts='n', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; RED BAND WITH CUTS
     Rmultirad, ap_radius_arcsec, mags, radiusresults_red_cuts, blobname[blob], nApertures, redname[blob], 'y' 
        ; namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, band=redname[blob], cuts='y', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; MIDDLE BAND WITH CUTS
     Rmultirad, ap_radius_arcsec, mags, radiusresults_mid_cuts, blobname[blob], nApertures, midname[blob], 'y' 
        ; namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, band=midname[blob], cuts='y', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; BLUE BAND WITH CUTS
     Rmultirad, ap_radius_arcsec, mags, radiusresults_blue_cuts, blobname[blob], nApertures, bluename[blob], 'y' 
        ; namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, band=bluename[blob], cuts='y', filetype='.png')
        ; plotsave, saveanyplots, namestring
     namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, filetype='.png')
     plotsave, saveanyplots, namestring
     ;
     ; more plots of the results of the aperture size analysis: overdensity factor vs aperture radius for each band for each mag bin 
     window, 0, retain=2, xsize=950, ysize=1200
     ; RED BAND, NO CUTS
     Rmultimag, mags, ap_radius_arcsec, radiusresults_red_nocuts, blobname[blob], nApertures, redname[blob], 'n' 
        ; namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, band=redname[blob], cuts='n', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; MIDDLE BAND, NO CUTS
     Rmultimag, mags, ap_radius_arcsec, radiusresults_mid_nocuts, blobname[blob], nApertures, midname[blob], 'n' 
        ; namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, band=midname[blob], cuts='n', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; BLUE BAND, NO CUTS
     Rmultimag, mags, ap_radius_arcsec, radiusresults_blue_nocuts, blobname[blob], nApertures, bluename[blob], 'n' 
        ; namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, band=bluename[blob], cuts='n', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; RED BAND WITH CUTS
     Rmultimag, mags, ap_radius_arcsec, radiusresults_red_cuts, blobname[blob], nApertures, redname[blob], 'y' 
        ; namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, band=redname[blob], cuts='y', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; MIDDLE BAND WITH CUTS
     Rmultimag, mags, ap_radius_arcsec, radiusresults_mid_cuts, blobname[blob], nApertures, midname[blob], 'y' 
        ; namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, band=midname[blob], cuts='y', filetype='.png')
        ; plotsave, saveanyplots, namestring
     ; BLUE BAND WITH CUTS
     Rmultimag, mags, ap_radius_arcsec, radiusresults_blue_cuts, blobname[blob], nApertures, bluename[blob], 'y' 
        ; namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, band=bluename[blob], cuts='y', filetype='.png')
        ; plotsave, saveanyplots, namestring
     namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, filetype='.png')
     plotsave, saveanyplots, namestring
     ;
   ENDIF     ; there are multiple aperture radii being tested 


 ENDFOR ;all the blobs


END             ; end of density procedure









PRO galbinning, brightestmag, dimmestmag, binsize, ra, dec, radius, data, maglist, galaxies_binned
;-------------------------------------------------------------------------------------------------------------------------------------------------
; galbinning procedure
;   description 
; INPUTS: etc
;         etc   
; NOTE: Uses get_galaxies_binned function. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION get_galaxies_binned
  maglist=list()
  FOR mag=brightestmag, dimmestmag, binsize DO BEGIN
    galaxies_binned = get_galaxies_binned(ra, dec, radius, data, mag, mag+binsize)
    maglist.Add, galaxies_binned
  ENDFOR 
END 









PRO get_ngal_binned, nbins, galaxy_maglist, ngal_binned, ngalerr_binned, ap=ap
;-------------------------------------------------------------------------------------------------------------------------------------------------
; get_ngal_binned procedure
;   description 
; INPUTS: etc
;         etc   
; NOTE: Uses Poisson_error function. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
FORWARD_FUNCTION Poisson_error
  FOR i=0, nbins-1. DO BEGIN
    bin = galaxy_maglist[i]
    IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin) 
    IF (n_elements(ap) EQ 0) THEN BEGIN 
      ngal_binned[i] = ngal 
      ngalerr_binned[i] = Poisson_error(ngal) 
    ENDIF ELSE IF (n_elements(ap) NE 0) THEN BEGIN 
      ngal_binned[i,ap] = ngal 
      ngalerr_binned[i,ap] = Poisson_error(ngal) 
    ENDIF 
  ENDFOR 
END 









PRO densitycalc, ngal_binned, ngalerr_binned, radius, densities, density_errors, weight=weight 
;-------------------------------------------------------------------------------------------------------------------------------------------------
; densitycalc procedure
;   description 
; INPUTS: etc
;         etc   
; NOTE: Uses Poisson_error function. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
  IF (n_elements(weight) EQ 0) THEN BEGIN 
    densities = double(ngal_binned)/(!dPI*radius^2.) 
    density_errors = double(ngalerr_binned)/(!dPI*radius^2.)
  ENDIF ELSE IF (n_elements(weight) NE 0) THEN BEGIN 
    densities = double(ngal_binned)/(!dPI*radius^2.*weight) 
    density_errors = double(ngalerr_binned)/(!dPI*radius^2.*weight)
  ENDIF
END 









PRO apfinder, npix_x, npix_y, ap_radius_pix, ap_radius, dim, blankim, bpm, hdr, blobra, blobdec, ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec 
       repeatflag = 0.       ; this will ensure that apertures don't get wasted on areas with only bad pixels 
       ap_count++
       WHILE (repeatflag EQ 0.) DO BEGIN
         ; make a random aperture:
         ap_x = floor(ap_radius_pix + (double(npix_x)-2.*ap_radius_pix)*randomu(seed))
         ap_y = floor(ap_radius_pix + (double(npix_y)-2.*ap_radius_pix)*randomu(seed)) 

         ; FIND BAD PIXELS:
          ; find corners of box in which to insert aperture image 
          xlo = ap_x - floor(ap_radius_pix) 
          xhi = ap_x + floor(ap_radius_pix)
          ylo = ap_y - floor(ap_radius_pix)  
          yhi = ap_y + floor(ap_radius_pix)
          dimsize = size(dim)
          IF ((xhi-xlo) NE (dimsize[1]-1.)) THEN BEGIN 
            xhi = ap_x + floor(ap_radius_pix) - 1.
          ENDIF
          IF ((yhi-ylo) NE (dimsize[2]-1.)) THEN BEGIN 
            yhi = ap_y + floor(ap_radius_pix) - 1.
          ENDIF
          ; create new image using array arithmetic: blank image with dim inserted into it 
          newap = blankim 
          newap[xlo:xhi,ylo:yhi] = dim 
          ; ignore everything outside the random aperture to make aperture map 
          index = where(newap LT ap_radius_pix) 
          wheretomulti, newap, index, col, row, frame
          newap[*,*] = 0. 
          newap[col,row]  = 1. 
          ; multiply the aperture map by the bpm so only GOOD pixels INSIDE the aperture are left: 
          apbpm = bpm*newap          
          badpix = where(apbpm GT 0, nbadpix)   ; label the bad pixels
          ntotalpix = n_elements(apbpm)         ; get the total number of pixels in the aperture

         ; FIND VICINITY TO BLOB: 
          xyad, hdr,ap_x, ap_y, ap_ra, ap_dec      ; convert pixels into coordnates in sky
          ; find how close the aperture is to the blob:
          vicinity = sqrt(  ( (blobra - ap_ra) * cos(blobdec*!dPI/180.) )^2. + (blobdec - ap_dec)^2. )

         ; if the aperture has at least some good pixels AND isn't near the blob, then keep it; otherwise, throw it out and make a new aperture  
         IF ((nbadpix LT ntotalpix) AND (vicinity GT 2.*ap_radius)) THEN BEGIN 
           repeatflag = 1. 
         ENDIF ELSE BEGIN 
           repeatflag = 0. 
           ap_count++    ; note that another aperture is being made since this one failed to meet the criteria for a good aperture 
         ENDELSE 
       ENDWHILE
       ; once the aperture passes the above WHILE test, it can be used
END 









PRO colordiagram, blobgals, fieldgals, datared, datamid, datablue, redname, midname, bluename, m, b, cap, blobname, radius
;-------------------------------------------------------------------------------------------------------------------------------------------------
; colordiagram procedure
;   description 
; INPUTS: etc
;         etc   
;         saveanyplots - a string, 'y' or 'n', determining if the user should be asked whether or not to save this plot 
;         namestring - a string containing the name of the file containing the plot if the user chooses to save it 
; NOTE: Uses plotsave procedure. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
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
   plot, color_x, color_y, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, psym=8, xrange=[-5,5],/xstyle, yrange=[-5,5],/ystyle, charsize=1.5, thick=2, ymargin=[4,4] 
   oplot, field_color_x, field_color_y, color=0, psym=3
   LEGEND, ['blob galaxy','field galaxy'], /left, /top, color=0, textcolor=0, psym=[8,3], charsize=1, charthick=1, /box, outline_color=0 
 ; add the cut line 
    xvalues = 0.01*findgen(100000) - 500.
    cut = xvalues*m + b
    cut[WHERE(cut GT cap)] = cap
    oplot, xvalues, cut, color=50, linestyle=2, thick=2
END 









PRO densityplot, nbins, mags, field_densities, field_density_errors, densities, density_errors, blobname, nApertures, radius, band, cuts
;-------------------------------------------------------------------------------------------------------------------------------------------------
; densityplot procedure
;   Makes plots of blob and field densities as a function of magnitude.
; INPUTS: blobname - name of blob whose density is being plotted (for plot title) 
;         nbins - number of magnitude bins being plotted on the x-axis 
;         mags - array of magnitudes (data for the x-axis) 
;         field_densities - array of field density measurements (data for the y-axis) 
;         field_density_errors - array of errors on field density measurements (to be overplotted as a swath on top of field_densities) 
;         densities - array of blob density measurements (also data for the y-axis) 
;         density_errors - array of errors on blob density measurements (to be overplotted on top of densities using errplot) 
;         apsnameplot - a string showing the number of random apertures used in the analysis (for plot title) 
;         xtitle - a string containing the label for the x-axis of the plot 
;         cuts - a string, 'y' or 'n', showing whether the densities used are for raw data or post-color-cut data (for plot title) 
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION titlemaker           ; this function is used in this procedure
 ; set up vertices of polygon for polyfill
   xpoints = dblarr(2*nbins)
   xpoints[0:nbins-1] = mags
   FOR i=0, nbins-1 DO BEGIN
     xpoints[i+nbins] = mags[nbins-1-i]
   ENDFOR 
   ; y is specific to each catalog (band) 
   ypoints=dblarr(2*nbins)
   ypoints[0:nbins-1] = field_densities-field_density_errors
   FOR i=0, nbins-1 DO BEGIN
     ypoints[i+nbins] = field_densities[nbins-1-i]+field_density_errors[nbins-1-i]
   ENDFOR
 ; make the main plot comparing blob to field with error bars 
   ; set up titles for plot and axes
      plottitle = titlemaker('overdensity', blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts)
      xtitle = 'magnitude (' + band + ')'
      ytitle='N!Igal!N per deg!E2!N per 0.5 mag'
   plot, mags, densities, background=255, color=0, title=plottitle, xtitle=xtitle, ytitle=ytitle, psym=-8, /ylog, yrange=[1d3,1d6],/ystyle, xrange=[22.5,29.5],/xstyle, charsize=2., xthick=2, ythick=2, xmargin=[10,7], ymargin=[6,6] 
   polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,1d6], /data, noclip=0
   errplot, mags, densities-density_errors, densities+density_errors, color=0, thick=2
   oplot, mags, densities, psym=-8, color=0, thick=2     ; have to do this again because polyfill covers it up 
   oplot, mags, field_densities, color=0, linestyle=2, thick=2
   LEGEND, ['blob','field'], /right, /top, color=0, textcolor=0, linestyle=[0,2], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5   
   axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=2, charthick=2, xtickformat="(A1)", xthick=2
   axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=2
   axis, color=0, yaxis=0, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=2, charthick=2, ytickformat="(A1)", ythick=2
   axis, color=0, yaxis=1, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=2 
END 








 
PRO statlum, mags, binsize, blob_ngal_binned, blob_ngal_binned_cuts, densities, field_densities, ap_radius, blobname, nApertures, ap_radius_arcsec, band   
;-------------------------------------------------------------------------------------------------------------------------------------------------
; statlum procedure
;   description 
; INPUTS: etc
;         etc   
;         saveanyplots - a string, 'y' or 'n', determining if the user should be asked whether or not to save this plot 
;         namestring - a string containing the name of the file containing the plot if the user chooses to save it 
; NOTE: Uses plotsave procedure. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION titlemaker           ; this function is used in this procedure
 ; define the statistical luminosity function
    overdensity = densities - field_densities
    n_overdensity = overdensity * (!dPI*ap_radius^2.)  ; convert to number 
 ; set up vertices of shaded area 
   xvertices = fltarr(2.0*n_elements(mags) + 2.)
   yvertices = fltarr(2.0*n_elements(mags) + 2.)
   xvertices[0] = mags[0]
   xvertices[-1] = mags[-1]
   yvertices[0] = 0
   yvertices[-1] = 0
   xvertices[1] = mags[0]
   xvertices[-2] = mags[-1]
   FOR i=1, n_elements(mags) -1 DO BEGIN
     xvertices[2*i] = mags[i]-(binsize/2.0)
     xvertices[2*i+1] = mags[i]-(binsize/2.0)
   ENDFOR
   FOR i=1, n_elements(mags) -1 DO BEGIN
     yvertices[2*i+1] = blob_ngal_binned_cuts[i]
     yvertices[2*i+2] = blob_ngal_binned_cuts[i]
   ENDFOR
 ; now plot statistical luminosity function
   ; set up titles for plot and axes
      title = titlemaker('statlum', blobname, nApertures=nApertures, radius=ap_radius_arcsec, band=band) 
      xtitle = 'magnitude (' + band + ')'
      ytitle = 'N!Igal!N'
   plot, mags, n_overdensity, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, linestyle=2, yrange=[0,10],/ystyle, xrange=[22,30],/xstyle, charsize=3, thick=2, xmargin=[10,7], ymargin=[4,4]
   polyfill, xvertices,yvertices, color=220, clip=[22,0,30,10], /data, noclip=0
   oplot, mags, blob_ngal_binned, color=0, psym=10, thick=2
   oplot, mags, blob_ngal_binned_cuts, color=0, psym=10, thick=2
   oplot, mags, n_overdensity, color=0, thick=2, linestyle=2
END 









PRO makehist, ap_totals, blob_ngal, blobname, nApertures, radius, cuts     ; maybe put optional band keyword in here or mag keyword or smth
;-------------------------------------------------------------------------------------------------------------------------------------------------
; makehist procedure
;   description 
; INPUTS: etc
;         etc   
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION titlemaker           ; this function is used in this procedure
 ; set up titles for plot and axes
    title = titlemaker('aphist', blobname, nApertures=nApertures, radius=radius, cuts=cuts) 
    xtitle = 'n galaxies in aperture'
    ytitle='number of apertures'
 plothist, ap_totals, background=255, color=0, axiscolor=0, title=title, xtitle=xtitle, ytitle=ytitle, bin=1, xrange=[0,blob_ngal+10.],/xstyle, charsize=1.5, thick=2, ymargin=[4,4]
 oplot, [blob_ngal, blob_ngal], [0., nApertures], linestyle=2, thick=2, color=0.
 LEGEND, ['blob','field'], /center, /top, color=0, textcolor=0, linestyle=[0,2], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5
; xyouts, blob_ngal+1., nApertures/10., blobgaltext, color=0, charsize=1.5
END









PRO Rmultirad, ap_radius_arcsec, mags, radiusresults, blobname, nApertures, band, cuts 
;-------------------------------------------------------------------------------------------------------------------------------------------------
; Rmultirad procedure
;   description 
; INPUTS: etc
;         etc   
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION titlemaker           ; this function is used in this procedure
 ; set up an array of colors to plot each radius in 
   n_radii = n_elements(ap_radius_arcsec)
   indices = findgen(n_radii)
   colors = indices*254./n_radii 
 ; now make the plot 
   ; set up titles for plot and axes
      title = titlemaker('Rmultirad', blobname, nApertures=nApertures, band=band, cuts=cuts)
      xtitle = 'magnitude (' + band + ')'
      ytitle = 'overdensity factor'
   plot, mags, radiusresults[0,*,0], background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, yrange=[0,50],/ystyle, charsize=2, thick=2, ymargin=[4,4] 
   errplot, mags, radiusresults[0,*,0]-radiusresults[0,*,1], radiusresults[0,*,0]+radiusresults[0,*,1], color=0, thick=2
   FOR r=1, n_radii-1. DO BEGIN
     oplot, mags, radiusresults[r,*,0], color = colors[r], thick=2
     errplot, mags, radiusresults[r,*,0]-radiusresults[r,*,1], radiusresults[r,*,0]+radiusresults[r,*,1], color=colors[r], thick=2 
   ENDFOR
   LEGEND, strcompress(string(fix(ap_radius_arcsec)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., number=0.5, charthick=1.5, thick=2
END 









PRO Rmultimag, mags, ap_radius_arcsec, radiusresults, blobname, nApertures, band, cuts 
;-------------------------------------------------------------------------------------------------------------------------------------------------
; Rmultimag procedure
;   description 
; INPUTS: etc
;         etc   
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION titlemaker           ; this function is used in this procedure
 ; set up an array of colors to plot each magnitude bin in 
   indices = findgen(n_elements(mags))
   colors = indices*254./float(n_elements(mags)) 
 ; now make the plot 
   ; set up titles for plot and axes 
      title = titlemaker('Rmultimag', blobname, nApertures=nApertures, band=band, cuts=cuts) 
      xtitle='aperture radius (arcseconds)'
      ytitle='overdensity factor'
   ; set up an empty plot (needs to have no data in case the first thing plotted is NaN): 
   plot, ap_radius_arcsec, radiusresults[*,1,0], background=255, color=colors[1], title=title, xtitle=xtitle, ytitle=ytitle, xrange=[0,30],/xstyle, yrange=[0,100],/ystyle, charsize=2, ymargin=[4,4], /nodata 
   ; now add data on top of the empty plot:  
   FOR i=0, n_elements(mags)-1 DO BEGIN     ; need to make sure that NaNs are thrown out (because where field density is 0, we divided by 0)
     y = radiusresults[*,i,0] 
     yerr = radiusresults[*,i,1] 
     keep = WHERE(finite(y) EQ 1, nkeep)
     IF (nkeep GT 0) THEN BEGIN
       oplot, ap_radius_arcsec(keep), y(keep), color=colors[i], thick=2
       errplot, ap_radius_arcsec(keep), y(keep)-yerr(keep), y(keep)+yerr(keep), color=colors[i], thick=2
     ENDIF
   ENDFOR
   ; add a legend: 
   LEGEND, strcompress(string(fix(mags)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., number=0.5, charthick=1.5, thick=2
END









PRO makereg, data, subset, blobname, radius, cuts 
;-------------------------------------------------------------------------------------------------------------------------------------------------
; makereg procedure
;   description 
; INPUTS: etc
;         etc   
; NOTE: Uses titlemaker function.
;-------------------------------------------------------------------------------------------------------------------------------------------------
 FORWARD_FUNCTION titlemaker           ; this function is used in this procedure
 cat = data[subset] 
 filename = titlemaker('blobgals', blobname, radius=radius, cuts=cuts, filetype='.reg') 
 openw, 1, filename 
 FOR i=0, n_elements(cat)-1 DO BEGIN
   printf, 1,'J2000; circle ', cat[i].alpha_j2000, cat[i].delta_j2000, ' 10p '    ;#text={', cat[i].number,'}'  
 ENDFOR 
 close, 1
END 









PRO plotsave, saveanyplots, namestring 
;-------------------------------------------------------------------------------------------------------------------------------------------------
; plotsave procedure
;   Asks user whether or not to save a plot and saves said plot if the user chooses to do so. 
; INPUTS: saveanyplots - a string, 'y' or 'n', determining if the user should be asked whether or not to save this plot at all (code will only
;                          run if this is set to 'y')  
;         namestring - a string containing the name of the file containing the plot if the user chooses to save it 
;-------------------------------------------------------------------------------------------------------------------------------------------------
  IF (saveanyplots EQ 'y') THEN BEGIN 
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
  ENDIF   ; user wants to save at least some plots
END 









FUNCTION titlemaker, plotname, blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts, filetype=filetype
;-------------------------------------------------------------------------------------------------------------------------------------------------
; titlemaker function
;   description 
; INPUTS: plotname - string containing the name of the kind of plot being made (ie, "overdensity," "color-color," etc)
;         blobname - string
;         nApertures - float
;         plotname - string containing the name of the kind of plot being made (ie, "overdensity," "color-color," etc)   
;         radius - float
;         band - string (F140W, F606W, etc) 
;         cuts - string/keyword  
;         filetype - specifies that this is going to be a filename if set; string that gives extension of the file 
;-------------------------------------------------------------------------------------------------------------------------------------------------
  ; first, make checks that see if the optional keywords are provided
  apcheck = n_elements(nApertures)       
  radcheck = n_elements(radius)
  bandcheck = n_elements(band)
  cutscheck = n_elements(cuts)
  filecheck = n_elements(filetype) 

  ; next make empty strings set up for each element of a title 
  apsName = ''                           
  radName = ''
  bandName = ''
  cutsName = ''
  titleName = '' 

  ; now make the final string that will eventually be returned
  titlestring = ''                    

  ; since the cuts keyword can only be 'y' or 'n', if it's set, make sure it's one of those things 
  IF (cutscheck NE 0) THEN BEGIN 
    IF ( (cuts NE 'y') AND (cuts NE 'n') ) THEN BEGIN
      WHILE ( (cuts NE 'y') AND (cuts NE 'n') ) DO BEGIN        ; make sure the keyword is valid
        print, 'Invalid cut keyword provided. Choose "y" or "n" only.'
        cuts = ''
        READ, cuts, PROMPT='Color cuts?' 
      ENDWHILE
    ENDIF ; cuts was set to something other than y or n
  ENDIF   ; the cuts keyword was provided at all

  ; now, determine whether to make a plot title or a file name:  

  IF (filecheck NE 0) THEN BEGIN            ; if the user wants a filename 
    ;
    titlestring = string(plotname) + '_' + string(blobname) 
    ; now tack on any extra keywords one by one
    ;
    ; start with nApertures
    IF (apcheck NE 0) THEN BEGIN 
      apsName = strcompress(string(fix(nApertures)), /remove) + 'aps'
      titlestring = titlestring + '_' + string(apsName) 
    ENDIF 
    ;
    ; radius
    IF (radcheck NE 0) THEN BEGIN 
      radName = strcompress(string(fix(radius)), /remove) + 'arcsec'
      titlestring = titlestring + '_' + string(radName)  
    ENDIF 
    ;
    ; band
    IF (bandcheck NE 0) THEN BEGIN 
      bandName = band
      titlestring = titlestring + '_' + string(bandName)  
    ENDIF 
    ;
    ; cuts 
    IF (cutscheck NE 0) THEN BEGIN 
      IF (cuts EQ 'y') THEN BEGIN        ; make a string for the case that there are color cuts 
        cutsName = 'cuts' 
      ENDIF ELSE IF (cuts EQ 'n') THEN BEGIN        ; make a string for data without cuts 
        cutsName = 'raw'
      ENDIF         ; now we've sorted out what string to use to describe the cuts keyword
     titlestring = titlestring + '_' + string(cutsName) 
    ENDIF 
    ;
    ; finally, add the file extension to the end of the string 
    titlestring = titlestring + string(filetype) 
    ; now we have a valid file name 

  ENDIF ELSE IF (filecheck EQ 0) THEN BEGIN         ; if the user is making a title for a plot 
    ;
    ; generate a basic plot title based on what's being plotted 
    IF (plotname EQ 'colors') THEN titleName = 'Color-Color Diagram for ' 
    IF (plotname EQ 'overdensity') THEN titleName = 'Galaxy Overdensity in Blob Region for '
    IF (plotname EQ 'statlum') THEN titleName = 'Statistical Luminosity Function for Galaxies in '
    IF (plotname EQ 'aphist') THEN titleName = 'Number of Galaxies in Apertures for '
    IF (plotname EQ 'Rmultirad') THEN titleName = 'Galaxy Overdensity in Different Radii for '
    IF (plotname EQ 'Rmultimag') THEN titleName = 'Dependence of Overdensity Factor on Aperture Radius for '
    titlestring = string(titleName) + string(blobname) 
    ;
    ; Now, if there are any keywords, add things to the basic plot title: 
    ;
    IF ( (apcheck NE 0) OR (radcheck NE 0) OR (bandcheck NE 0) OR (cutscheck NE 0) ) THEN BEGIN
      ; add a comma and a line after the main title, since this is extra stuff 
      titlestring = titlestring + ',!C'
      titlestring = string(titlestring) 
      ;
      ; First, see if nApertures and/or radius are specified. Since they go hand-in-hand, treat them together: 
      IF (radcheck NE 0) THEN radName = strcompress(string(fix(radius)), /remove) + '"'     ; if radius is provided, here's a string describing it 
      IF (apcheck NE 0) THEN BEGIN            ; if nApertures is provided 
        apsName = strcompress(string(fix(nApertures)), /remove) + ' apertures'    ; string describing nApertures
        IF (radcheck NE 0) THEN BEGIN        ; if nApertures and radius are both provided
          titlestring = titlestring + string(apsName) + ' with ' + string(radName) + ' radii'   
        ENDIF ELSE IF (radcheck EQ 0) THEN BEGIN            ; if nApertures is provided but not radius
          titlestring = titlestring + string(apsName)  
        ENDIF
      ENDIF ELSE IF ( (apcheck EQ 0) AND (radcheck NE 0) ) THEN BEGIN     ; if nApertures is not provided but radius is 
        titlestring = titlestring + 'aperture radius of ' + string(radName)       
      ENDIF         ; if neither nApertures nor radius are provided, nothing else is added to the title string yet 
      ;
      ; now tack on the name of the band if that has been provided
      IF (bandcheck NE 0) THEN BEGIN 
        bandName = string(band) 
        IF ( (apcheck NE 0) OR (radcheck NE 0) ) THEN BEGIN     ; if anything else is provided before band, add a comma
          titlestring = titlestring + ', ' 
        ENDIF
        titlestring = titlestring + string(bandName)  
      ENDIF    ; the band is provided 
      ;
      ; now add info about color cuts if that has been provided 
      IF (cutscheck NE 0) THEN BEGIN     
        IF (cuts EQ 'y') THEN BEGIN        ; make a string for the case that there are color cuts 
          cutsName = 'with color cuts' 
        ENDIF ELSE IF (cuts EQ 'n') THEN BEGIN        ; make a string for data without cuts 
          cutsName = 'raw (no cuts)'
        ENDIF         ; sorts out what string to use to describe the cuts keyword
        ; now that the string is set up, add it to the title 
        IF ( (apcheck NE 0) OR (radcheck NE 0) OR (bandcheck NE 0) ) THEN BEGIN 
          titlestring = titlestring + ', '         ; add a comma if anything comes before the cuts text
        ENDIF 
        titlestring = string(titlestring) + string(cutsName)       ; tack on the cuts text 
      ENDIF      ; the cuts keyword is provided 
      ;
    ENDIF     ; there are any keywords specified at all
    ;
  ENDIF    ; this is the title for a plot 

  
  RETURN, string(titlestring)       ; provide the final title for the plot or file 

END 









FUNCTION get_galaxies, ra, dec, radius, data
;-------------------------------------------------------------------------------------------------------------------------------------------------
; get_galaxies function
; INPUTS: ra, dec - the RA and declination of a point at which to place an aperture 
;         radius - the radius of the aperture placed at the given coordinates
;         data - a SExtractor catalog of sources which contains each source's RA, declination, a flag indicating what type of object it is, 
;                  and a flag indicating whether or not there are any bad pixels in the source in the original image
; OUTPUT: galaxies_in_aperture - the index numbers of the galaxies in the data catalog which are located within the specified aperture 
; NOTES: The output is found by determining the distance between the aperture's center and each source, then requiring that distance to be less
;           than the radius of the aperture, filtering out all objects which are marked as stars in the catalog, and filtering out all objects 
;           with flags that indicate that they lie in areas of bad pixels. This ensures that only real sources that are galaxies are returned. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
  distance = sqrt(  ( (ra - data.ALPHA_J2000) * cos(dec*!dPI/180.) )^2. + (dec - data.DELTA_J2000)^2. )         ; Pythagorean theorem 
  ; the source has to be contained within the aperture, has to be a galaxy (no stars), and has to have good data (no bad pixels) 
  galaxies_in_aperture = WHERE (((distance LT radius) AND (data.CLASS_STAR LT 0.8) AND (data.IMAFLAGS_ISO EQ 0.)), n_galaxies)  
  RETURN, galaxies_in_aperture
END









FUNCTION get_galaxies_noblob, ra, dec, radius, data
;-------------------------------------------------------------------------------------------------------------------------------------------------
; get_galaxies_noblob function
; INPUTS: ra, dec - the RA and declination of a point at which to place an aperture 
;         radius - the radius of the aperture placed at the given coordinates
;         data - a SExtractor catalog of sources which contains each source's RA, declination, a flag indicating what type of object it is, 
;                  and a flag indicating whether or not there are any bad pixels in the source in the original image
; OUTPUT: galaxies_outside_aperture - the index numbers of the galaxies in the data catalog which are located OUTSIDE the specified aperture 
; NOTES: The output is found by determining the distance between the aperture's center and each source, then requiring that distance to be greater
;           than the radius of the aperture, filtering out all objects which are marked as stars in the catalog, and filtering out all objects 
;           with flags that indicate that they lie in areas of bad pixels. This ensures that only real sources that are galaxies are returned. 
;        This function essentially does the exact opposite of what get_galaxies does. 
;-------------------------------------------------------------------------------------------------------------------------------------------------
  distance = sqrt(  ( (ra - data.ALPHA_J2000) * cos(dec*!dPI/180.) )^2. + (dec - data.DELTA_J2000)^2. )         ; Pythagorean theorem 
  ; the source has to be located outside the aperture, has to be a galaxy (no stars), and has to have good data (no bad pixels) 
  galaxies_outside_aperture = WHERE (((distance GE radius) AND (data.CLASS_STAR LT 0.8) AND (data.IMAFLAGS_ISO EQ 0.)), n_galaxies)
  RETURN, galaxies_outside_aperture
END









FUNCTION get_galaxies_binned, ra, dec, radius, data, brightmag, dimmag
;------------------------------------------------------------------------------------------------------------------------------------------------
; get_galaxies function
; INPUTS: ra, dec - the RA and declination of a point at which to place an aperture 
;         radius - the radius of the aperture placed at the given coordinates
;         data - a SExtractor catalog of sources which contains each source's RA, declination, magnitude, a flag indicating what type of object 
;                  it is, and a flag indicating whether or not there are any bad pixels in the source in the original image
;         brightmag - the lower bound (brightest magnitude) of a magnitude range (bin) 
;         dimmag - the upper bound (faintest magnitude) of a magnitude range (bin) 
; OUTPUT: galaxies_in_aperture - the index numbers of the galaxies in the data catalog which are located within the specified aperture and have
;                                  magnitudes between brightmag and dimmag
; NOTES: The output is found by determining the distance between the aperture's center and each source, then requiring that distance to be less
;           than the radius of the aperture, filtering out all objects which are marked as stars in the catalog, filtering out all objects 
;           with flags that indicate that they lie in areas of bad pixels, and filtering out all objects with magnitudes outside the range set
;           by brightmag and dimmag. This ensures that only real sources that are galaxies with specific magnitudes are returned. 
;------------------------------------------------------------------------------------------------------------------------------------------------
  distance = sqrt(  ( (ra - data.ALPHA_J2000) * cos(dec*!dPI/180.) )^2. + (dec - data.DELTA_J2000)^2. )         ; Pythagorean theorem 
  ; the source has to be contained within the aperture, has to be a galaxy (no stars), has to have good data (no bad pixels), and has to have a magnitude between brightmag and dimmag  
  galaxies_in_aperture = WHERE (((distance LT radius) AND (data.CLASS_STAR LT 0.8) AND (data.IMAFLAGS_ISO EQ 0.) AND (data.MAG_ISO GE brightmag) AND (data.MAG_ISO LT dimmag)), n_galaxies)
  RETURN, galaxies_in_aperture
END









FUNCTION Poisson_error, ngal 
;------------------------------------------------------------------------------------------------------------------------------------------------
; Poisson_error function 
; INPUT: ngal - a number (of galaxies) for which we want to determine uncertainty
; OUTPUT: error - the Poisson error on the input number (the square root of the input: N = sqrt(S) ) 
;------------------------------------------------------------------------------------------------------------------------------------------------
  error = sqrt(double(ngal))
  RETURN, error
END 









FUNCTION colorsort, data, xcolor, ycolor, m, b, c 
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
  good = data[WHERE ((ycolor GE (xcolor*m + b)) OR (ycolor GE c))]
  RETURN, good
END 
  








; NOTE: THE FOLLOWING FUNCTION WAS NOT WRITTEN BY AGNAR. 
;+
; NAME:
;       RSEX()
;
; PURPOSE:
;       Read in arbitrary SExtractor format catalogs, using native
;       header information & data themselves.  Correctly reads longs,
;       strings, and doubles.  
;
; INPUTS:
;       A SExtractor-format catalog
;
; OPTIONAL INPUTS:
;       use_row - use this row of the catalog to determine the format
;                 of the output data structure (zero-indexed; default
;                 0); this is necessary because occasionally a column
;                 in the first row is a different format (e.g., LONG)
;                 than the rest of the rows (e.g., DOUBLE)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Returns a structure with all catalog entries, using field
;       names for tagnames. 
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;       Use syntax
;       cs = rsex('catalog.cat')
;
; COMMENTS:
;
; PROCEDURES USED:
;       FINDFILE
;       CREATE_STRUCT
;       VALID_NUM
;       VALID_NUM_ARR
;
; MODIFICATION HISTORY:
; LAM = L. Moustakas
; JM  = J. Moustakas
;       LAM '04may02 - fixed problem case of there being an array of
;                      values based on the last header tag position.
;                      will now work with that case after special
;                      check. 
;       LAM '04feb04 - converted from the old lrsex.pro; adapted to
;                      detect and read longs, strings, and doubles.  
;       JM  '04sep22 - skip comment lines in the header designated
;                      with a double-hash (##)
;       JM  '07aug18 - added USE_ROW optional input
;       JM  '08jul25 - use a variable logical unit number
;       LAM+JM '08dec11 - vetted - no changes made
;       jm14jul12siena - use IDL_VALIDNAME() to ensure valid structure tags
;-

FUNCTION rsex,catalog, use_row=use_row

; Check that an argument has been passed     
     IF n_params() LE 0 THEN BEGIN 
         print,'cat=rsex(catalog)'
         return,-1
     ENDIF 

; See if the catalog file exists at all

     if (file_test(catalog,/regular) eq 0L) then begin ; jm04sep20uofa
         print,'error - the specified catalog does not seem to exist'
         return,-1
     endif

;    jnk = findfile(catalog,count=catexist)
;    if catexist eq -1 then begin
;        message,'error - the specified catalog does not seem to exist'
;        return,-1
;    endif

; Count catalog lines
     spawn,'wc -l < '+catalog, nlines, /sh
     nlines = long(nlines[0]);-1

; Open the catalog to read
     openr,lun,catalog,error=err,/get_lun
     if err ne 0 then begin
         message,'error opening catalog - '+!error_state.msg ; jm04sep20uofa
;        message,'error opening catalog - ',!error_state.msg
         free_lun, lun
         return,-1
     endif

     if (n_elements(use_row) eq 0L) then use_row = 0L ; jm07aug18nyu
     
; Read in the header
     head = ''
     numb = 0l
     cstr=''
     tag = 1
     nheadcomment = 0L ; jm04sep20uofa
     while tag do begin         ; while '#' tag is true
         readf,lun,cstr
         if (strmatch(cstr,'##*') eq 0B) then begin ; skip comment lines
            if strpos(cstr,'#') ne -1 then begin
               head = [head, strupcase((str_sep(strcompress(cstr),' '))[2])]
               numb = [numb, (str_sep(strcompress(cstr),' '))[1]]
            endif else tag=0
         endif else nheadcomment = nheadcomment + 1L ; jm04sep20uofa
     endwhile
     nhead = n_elements(head)   ; number of header entries
     head = head[1:(nhead-1)]
     numb = numb[1:(nhead-1)]
     nhead = n_elements(head)
     free_lun, lun

; Check header for implicit array entries, and expand relevant tagnames
     tothead=''
     for i=0l,nhead-2 do begin
         tothead = [tothead,head[i]]
         mlen = numb[i+1]-numb[i]-1
         if mlen ne 0 then $
           for j=1,mlen do $
           tothead=[tothead,head[i]+strcompress(j,/remove_all)]
     endfor

; The header info so far
     ntothead = n_elements(tothead)
     tothead=[tothead[1:(ntothead-1)],head[(nhead-1)]]

; Need to check one more thing -- whether the last header field is
; actually the first one of an array...

; Read the first data line.  If it has more entries than the total
; number of header fields so far, we've missed an array at the end. 
     junk=strarr(nhead+nheadcomment) ; jm04sep20uofa
     openr,lun,catalog,/get_lun
     readf,lun,junk
     readf,lun,cstr
     free_lun, lun
     nstr = n_elements(str_sep(strcompress(strtrim(cstr,2)),' '))

     if nstr gt ntothead then begin
         mlen = nstr - ntothead
         for j=1,mlen do $
           tothead=[tothead,head[(nhead-1)]+strcompress(j,/remove_all)]
     endif
     ntothead=n_elements(tothead)

; That should do it!  Onwards, to read the data...
     nbody=nlines-nhead-nheadcomment ; jm04sep20uofa
     junk=strarr(nhead+nheadcomment) ; jm04sep20uofa
     openr,lun,catalog
     readf,lun,junk
     for ijunk = 0L, use_row-1L do readf,lun,junk2 ; jm07aug18nyu
     readf,lun,cstr
     free_lun, lun
     body=strarr(nbody)
     openr,lun,catalog
     readf,lun,junk
     readf,lun,body
     free_lun, lun
     bodyslam=strarr(ntothead,nbody)
     for i=0l,nbody-1 do begin
;       print, i, body[i]
;       bodyslam[*,i]=(str_sep(strcompress(strtrim(body[i],2)),' '))[0L:nhead-1L]
        bodyslam[*,i]=(str_sep(strcompress(strtrim(body[i],2)),' '))[0L:ntothead-1L] ; jm04nov22uofa
     endfor
     carr=str_sep(strcompress(strtrim(cstr,1)),' ')

     tind    = valid_num(carr) ; 0=string 1=number
     tindint = valid_num(carr,/int) ; 0=non-integer 1=integer
     if tind[0] eq 0 then $
       type='' else $
       if tindint[0] eq 1 then $
       type=0l else $
       type=0.d0
     cs=create_struct(idl_validname(tothead[0],/convert_all),type)

     for i=1l,ntothead-1 do begin
         if tind[i] eq 0 then $
           type='' else $
           if tindint[i] eq 1 then $
           type=0l else $
           type=0.d0
         cs=create_struct(cs,idl_validname(tothead[i],/convert_all),type)
     endfor
     cs=replicate(cs,nbody)

     for i=0l,ntothead-1 do begin
        for j=0l,nbody-1 do begin
           cs[j].(i)=bodyslam[i,j]
;          print, bodyslam[i,j], i, j
;          print, i, j
        endfor
;stop
     endfor
     
     return,cs
END 





;
;                                   ||`-.___
;                                   ||    _.>
;                                   ||_.-'
;               ==========================================
;                `.:::::::.       `:::::::.       `:::::::.
;                  \:::::::.        :::::::.        :::::::\
;                   L:::::::     AGNAR L HALL        :::::::L
;                   J::::::::        ::::::::        :::::::J
;                    F:::::::        ::::::::        ::::::::L
;                    |:::::::        ::::::::        ::::::::|
;                    |:::::::        ::::::::        ::::::::|    ^..---.__
;                    |:::::::        ::::::::        ::::::::|   |  <)    -`.
;                    |:::::::        ::::::::        ::::::::|   |    /^^^^/
;     __             |:::::::        ::::::::        ::::::::|    \ . \vvv\
;   .'_ \            |:::::::        ::::::::        ::::::::|     \ `----`
;   (( ) |           |:::::::        ::::::::        ::::::::|      \ `.
;    `/ /            |:::::::        ::::::::        ::::::::|       L  \
;    / /             |:::::::        ::::::::        ::::::::|       |   \
;   J J              |:::::::        ::::::::        ::::::::|       |    L
;   | |              |:::::::        ::::::::        ::::::::|       |    |
;   | |              |:::::::        ::::::::        ::::::::|       F    |
;   | J\             F:::::::        ::::::::        ::::::::F      /     |
;   |  L\           J::::::::       .::::::::       .:::::::J      /      F
;   J  J `.     .   F:::::::        ::::::::        ::::::::F    .'      J
;    L  \  `.  //  /:::::::'      .::::::::'      .::::::::/   .'        F
;    J   `.  `//_..---.   .---.   .---.   .---.   .---.   <---<         J
;     L    `-//_=/  _  \=/  _  \=/  _  \=/  _  \=/  _  \=/  _  \       /
;     J     /|  |  (_)  |  (_)  |  (_)  |  (_)  |  (_)  |  (_)  |     /
;      \   / |   \     //\     //\     //\     //\     //\     /    .'
;       \ / /     `---//  `---//  `---//  `---//  `---//  `---'   .'
;________/_/_________//______//______//______//______//_________.'_________
;##########################################################################
;
;
;
;   CODE APPROVED BY DOGE
;░░░░░░░░░▄░░░░░░░░░░░░░░▄░░░░
;░░░░░░░░▌▒█░░░░░░░░░░░▄▀▒▌░░░
;░░░░░░░░▌▒▒█░░░░░░░░▄▀▒▒▒▐░░░
;░░░░░░░▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐░░░
;░░░░░▄▄▀▒░▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐░░░
;░░░▄▀▒▒▒░░░▒▒▒░░░▒▒▒▀██▀▒▌░░░
;░░▐▒▒▒▄▄▒▒▒▒░░░▒▒▒▒▒▒▒▀▄▒▒▌░░
;░░▌░░▌█▀▒▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐░░
;░▐░░░▒▒▒▒▒▒▒▒▌██▀▒▒░░░▒▒▒▀▄▌░
;░▌░▒▄██▄▒▒▒▒▒▒▒▒▒░░░░░░▒▒▒▒▌░
;▀▒▀▐▄█▄█▌▄░▀▒▒░░░░░░░░░░▒▒▒▐░
;▐▒▒▐▀▐▀▒░▄▄▒▄▒▒▒▒▒▒░▒░▒░▒▒▒▒▌
;▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒▒▒░▒░▒░▒▒▐░
;░▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒░▒░▒░▒░▒▒▒▌░
;░▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▒▄▒▒▐░░
;░░▀▄▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▄▒▒▒▒▌░░
;░░░░▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀░░░
;░░░░░░▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀░░░░░
;░░░░░░░░░▒▒▒▒▒▒▒▒▒▒▀▀░░░░░░░░
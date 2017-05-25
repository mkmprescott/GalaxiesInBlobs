;  LIONS: LABs In Overdense Neighborhoods Study 
;                                                           /   \       
; _                                                 )      ((   ))      (
;(@)                                               /|\      ))_((      /|\                                                 _ 
;|-|                                              / | \   <_`\|/`_>   / | \                                               (@)
;| |---------------------------------------------/--|-voV--\<>|<>/--Vov-|--\----------------------------------------------|-|
;|-|                                                  '^`   (> <)   '^`                                                   | |
;| |                                                        `\Y/'                                                         |-|
;|-| LIONS.PRO                                                                                                            | |
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
      FORWARD_FUNCTION get_galaxies
      FORWARD_FUNCTION get_galaxies_binned
      FORWARD_FUNCTION get_galaxies_noblob
      FORWARD_FUNCTION Poisson_error
      FORWARD_FUNCTION colorsort
      FORWARD_FUNCTION titlemaker 

   ; USER-DEFINED VARIABLES: 
      pixelscale = 0.06               ; arceseconds per pixel for our images (may eventually become part of blobs.csv) 
      hstpixelscale = 0.06            ; I've also seen 0.128 "/px for this?
      ap_radius_arcsec = [10.]        ; aperture radius in arcsec
      binsize=0.5                     ; size of magnitude bins 
      brightestmag = 22.              ; smallest (brightest) magnitude of galaxies we can reasonably see in our images 
      dimmestmag = 30.                ; biggest (faintest) magnitude of galaxies we can reasonably see in our images 
      nApertures = 10000.              ; number of random apertures to use in field density calculation 
      saveanyplots = 'y'              ; string indicating that user should be asked whether or not to save plots 
      savereg = 'n'                   ; string indicating that user should be asked whether or not to save region files
TIC 
   ; read in the file containing all blob information
      blobs = read_csv('blobs.csv', n_table_header=1)
   ; assign names to each field of the "blobs" structure
      blobname = blobs.FIELD01
      blobra = blobs.FIELD02
      blobdec = blobs.FIELD03
      ; fields 4-7 are irrelevant information for this analysis 
      ;blueband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobCatalogs/'+blobs.FIELD04
      ;midband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobCatalogs/'+blobs.FIELD05
      ;redband = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobCatalogs/'+blobs.FIELD06
      ;stack = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobStacks/'+blobs.FIELD07
      mask = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/BlobData/BlobStacks/'+blobs.FIELD08 
      z = blobs.FIELD09                 ; redshift
      m = blobs.FIELD10                 ; slope of color-cut line 
      b = blobs.FIELD11                 ; y-intercept of color-cut line 
      cap = blobs.FIELD12               ; mag above which all colors are consistent with redshift
      bluename = blobs.FIELD13          ; name of the bluest filter used to observe the blob region
      midname = blobs.FIELD14           ; name of the middle filter used to observe the blob region
      redname = blobs.FIELD15           ; name of the reddest filter used to observe the blob region
      blobinfo = blobs.FIELD16          ; name of the catalog containing all the magnitude, RA, DEC, etc info of the sources
      nblobs = n_elements(blobname)     ; number of blobs in the blobs.csv file 

   ; read in 3D-HST field data and rename fields to be useful
      goodss = read_csv('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/goodss.csv')   
      struct_replace_field, goodss, 'FIELD1', goodss.FIELD1, newtag='ID'
      struct_replace_field, goodss, 'FIELD2', goodss.FIELD2, newtag='ra'
      struct_replace_field, goodss, 'FIELD3', goodss.FIELD3, newtag='dec'
      struct_replace_field, goodss, 'FIELD4', goodss.FIELD4, newtag='redmag'
      struct_replace_field, goodss, 'FIELD5', goodss.FIELD5, newtag='midmag'
      struct_replace_field, goodss, 'FIELD6', goodss.FIELD6, newtag='bluemag'   
      struct_replace_field, goodss, 'FIELD7', goodss.FIELD7, newtag='altbluemag'
      struct_replace_field, goodss, 'FIELD8', goodss.FIELD8, newtag='z'
   ; separate this into two possible catalogs to use depending on which blue band each blob has 
      use = [WHERE (finite(goodss.bluemag) EQ 1)]
      goodsswithblue = {ID:goodss.ID[use], ra:goodss.ra[use], dec:goodss.dec[use], redmag:goodss.redmag[use], midmag:goodss.midmag[use], bluemag:goodss.bluemag[use], z:goodss.z[use]}
      use = [WHERE (finite(goodss.altbluemag) EQ 1)]
      goodsswithalt = {ID:goodss.ID[use], ra:goodss.ra[use], dec:goodss.dec[use], redmag:goodss.redmag[use], midmag:goodss.midmag[use], bluemag:goodss.altbluemag[use], z:goodss.z[use]}
TOC 
   ; read in HST bad pixel map (for field analysis) 
TIC 
      hsthdr =  headfits('/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/3DHST/3DHSTimages/stacks/goodss_3dhst.v4.0.F125W_F140W_F160W_det.fits')
      hstmask = '/boomerang-data/alhall/LABoverdensity/GalaxiesInBlobs/Data/3DHST/3DHSTimages/masks/GOODS-S_area.fits'
      hstbpm = -1.*(mrdfits(hstmask,0) - 1.)      ; this is the mask image for weeding out bad pixels in 3DHST field, in the same format as our blob bpms 
TOC 
TIC 
      hstsizeinfo = size(hstbpm)               ; the dimensions of the mask image in pixels 
      hstnpix_x = hstsizeinfo[1]               ; number of pixels in mask in the x direction
      hstnpix_y = hstsizeinfo[2]               ; number of pixels in mask in the y direction 

   ; define other important quantities  
      ; the radius/radii of the aperture(s) used to calculate density:  
         ap_radius = ap_radius_arcsec/3600.             ; aperture radius in decimal degrees
         ap_radius_pix = ap_radius_arcsec / pixelscale  ; aperture radius in pixels  
         hstap_radius_pix = ap_radius_arcsec / hstpixelscale
         n_radii = n_elements(ap_radius)                ; number of aperture sizes being tested
      ; mag binning stuff:
         nbins = (dimmestmag - brightestmag + binsize)/binsize     ; the number of magnitude bins as determined from our range and binsize 
         mags = brightestmag + binsize*findgen(nbins)              ; an array of every magnitude being used 
         magrange=[brightestmag,dimmestmag]                        ; a two-element array containing only the endpoints of our magnitude range 
      ; plotting symbol: make the psym=8 symbol a filled circle for our plots 
         plotsym, 0, 1, /fill   
TOC 

;**********************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!????????????????*****************PICK UP HERE*****************************************************************************************************************************************************************

     print, 'Now sampling the field for overdensity analysis.' 
     ;
   FOR r=0, n_radii-1. DO BEGIN 
     ; set up the arrays which will contain all the information about the random apertures placed for HST field: 
     hstap_coords = dblarr(2., nApertures)     ; this will contain RAs and declinations of all apertures
     hstap_ngal_binned = dblarr(nbins,4.,nApertures)   ; this will contain # galaxies in each aperture in each mag bin in each band 
     hstap_ngalerr_binned = dblarr(nbins,4.,nApertures)   ; same as above but errors on # rather than just #
     hstmegalist = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin in each band
dummyarray = lonarr(150)
dummyarray = dummyarray - 1
dummystruct = {aperture:0L, RA:0.0, dec:0.0, galIDs:dummyarray}
help, dummystruct
HSTstruct = replicate(dummystruct,nApertures)
     ;
     ; make an aperture image where everything inside the aperture has a value of 1 and everything outside (just the edge stuff) has a value of 0:
     dist_ellipse, hstdim, [ 2.*hstap_radius_pix[r], 2.*hstap_radius_pix[r] ], hstap_radius_pix[r], hstap_radius_pix[r], 1.,0.
       hstdim[where(hstdim LT ap_radius_pix[r], hstntotalpix)] = 1.           ; note that hstntotalpix is defined in this line as well
       hstdim[where(hstdim GE ap_radius_pix[r])] = 0.
     ; keep track of how many apertures are made in total
     hstap_count = 0     ; for the 3D-HST field
     ;
     ; now place random apertures and count galaxies inside them 
     ; FOR THE 3D-HST FIELD 
     print, 'GOODS-S FIELD ANALYSIS'
       TIC        ; time this to ensure that it's efficient 
;openw, lun, 'teststruct2thetestening.txt', /get_lun
       FOR ap=0, nApertures-1 DO BEGIN
         ; find a good aperture to use: 
            apfinder, hstnpix_x, hstnpix_y, hstap_radius_pix[r], ap_radius[r], hstdim, hstbpm, hsthdr, hstap_count, hstap_x, hstap_y, hstnbadpix, hstntotalpix, hstap_ra, hstap_dec 
         ; calculate fraction of aperture that is good pixels:
            hstap_area_weight = (float(hstntotalpix) - float(hstnbadpix)) / float(hstntotalpix) 
         ; convert pixels to sky coords and add these coords to arrays of aperture properties
            xyad, hsthdr, hstap_x, hstap_y, hstap_ra, hstap_dec    ; for HST field
            hstap_coords[0,ap] = hstap_ra
            hstap_coords[1,ap] = hstap_dec
HSTstruct(ap).aperture = ap+1
HSTstruct(ap).RA = hstap_ra
HSTstruct(ap).dec = hstap_dec
         ; create array of galaxy IDs for each aperture
            aparray = fltarr(1)
            aparray[0] = ap+1     ; first element of array is aperture number! 
              FOR mag=brightestmag, dimmestmag, binsize DO BEGIN  
                hst_galaxies_binned_indices = get_galaxies_binned(hstap_ra, hstap_dec, ap_radius[r], goodss, mag, mag+binsize, goodss.(3)) 
                  IF (hst_galaxies_binned_indices[0] GT -1.0) THEN BEGIN
                    hst_galaxies_binned = goodss.ID[hst_galaxies_binned_indices]
                    aparray = [aparray, hst_galaxies_binned]  ; add galaxy IDs onto array
                  ENDIF   ; if there are any galaxies in this bin at all 
              ENDFOR     ; all magnitude bins 
ap_ngals = n_elements(aparray)          ; note that this includes aperture number
IF (ap_ngals GT 1) THEN BEGIN
lastindex = ap_ngals - 1                ; last index of aparray
arrayforstruct = aparray[1:lastindex]   ; take the aperture number out of aparray in order to get just galaxy IDs for the structure
newlastindex = lastindex-1                 ; since one element was removed from the array, there's a different last index now
HSTstruct(ap).galIDs[0:newlastindex] = arrayforstruct
ENDIF
;printf, lun, aparray
;printf, lun, ' '
;print, hstap_ra, hstap_dec  ; check to make sure there's no glitchy overlapping apertures
;print, ' '
         ; put the array of ID numbers in this aperture into the mega list for all apertures:
;            hstmegalist.Add, aparray  
;         ; use apmaglist to fill in the ngals arrays
;            get_ngal_binned, nbins, hstapmaglist, hstap_ngal_binned, hstap_ngalerr_binned, 4, ap=ap
       ENDFOR    ; all the random apertures 
;printf, lun, ' '
;printf, lun, ' '
;printf, lun, ' '
;printf, lun, 'Structure: '
;printf, lun, HSTstruct
mwrfits, HSTstruct, 'HSTapertures.fits', /create
help, HSTstruct
testHSTstruct = mrdfits('HSTapertures.fits',1)
help, testHSTstruct
;printf, lun, ' '
;printf, lun, ' '
;printf, lun, ' '
;printf, lun, 'TEST: Does writing structure to FITS file work? '
;printf, lun, testHSTstruct

;close, lun
TOC             ; How long did this analysis take? 
STOP
;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       print, hstap_count ; How many total apertures were generated in HST field? 
       hstthrownout = (float(hstap_count) - float(nApertures))/float(hstap_count)*100.     ; percentage of generated apertures that ended up being unsuitable for analysis 
     ; as a safety check, see how many apertures ended up being unsuitable in each analysis 
     print, strcompress(string(hstthrownout)), '% of GOODS-S apertures thrown out' 
     ;
     ; calculate densities in each bin in each aperture and errors: 
        densitycalc, hstap_ngal_binned, hstap_ngalerr_binned, ap_radius[r], hstap_densities, hstap_density_errors, weight=hstap_area_weight
     ; then average the densities of all the apertures in each bin to get average field density in each bin (for each band):
        hstfield_densities = total(hstap_densities,3)/double(nApertures)             ; field density in each bin in each band for HST 
        hstfield_density_errors = total(hstap_density_errors,3)/double(nApertures)   ; errors in field density in each bin in each band for HST 
     ; Now, we want to make a catalog of galaxy IDs in each aperture 
help, hstmegalist 
        apstruct = {aperture:0L, galaxyid:make_array(150,value=-1)}
        fieldstruct = replicate(apstruct, nApertures)
        FOR i=1, nApertures DO BEGIN
          fieldstruct[i].aperture = i
apthing = hstmegalist[i]
bandthing = thing[0]
;  NOTE: ACTUALLY, MAKE GETTING GALAXIES IN APERTURE AND BINNING THEM SEPARATE THINGS 
;       instead of get->bin, count->density, do get, count, bin, density

          fieldstruct[i].galaxyid = hstmegalist[i[0]]
        ENDFOR
ENDFOR ; all radii

help, fieldstruct
STOP

;***************************************************************************************PICK UP HERE*****************************************************************************************************************************************************************

 


 ; DENSITY ANALYSIS FOR EACH BLOB 

 FOR blob=0, nblobs-1 DO BEGIN         

   ; make structure for the blob catalog and rename the fields to be useful  
      blobdata = read_csv(blobinfo[blob]) 
      struct_replace_field, blobdata, 'FIELD1', blobdata.FIELD1, newtag='ID'
      struct_replace_field, blobdata, 'FIELD2', blobdata.FIELD2, newtag='ra'
      struct_replace_field, blobdata, 'FIELD3', blobdata.FIELD3, newtag='dec'
      struct_replace_field, blobdata, 'FIELD4', blobdata.FIELD4, newtag='redmag'
      struct_replace_field, blobdata, 'FIELD5', blobdata.FIELD5, newtag='midmag'
      struct_replace_field, blobdata, 'FIELD6', blobdata.FIELD6, newtag='bluemag'

   ; make new catalogs with color cuts
      xcolorall = blobdata.bluemag - blobdata.midmag   ; x color in a color-color diagram 
      ycolorall = blobdata.midmag - blobdata.redmag    ; y color in a color-color diagram 
      blobcuts = colorsort(xcolorall, ycolorall, m[blob], b[blob], cap[blob])
      blobdatacuts = {ID:blobdata.ID[blobcuts], ra:blobdata.ra[blobcuts], dec:blobdata.dec[blobcuts], $
         redmag:blobdata.redmag[blobcuts], midmag:blobdata.midmag[blobcuts], bluemag:blobdata.bluemag[blobcuts]}

   ; now do the same for 3DHST - but first figure out which blue band to use 
      IF (bluename[blob] EQ 'F606W') THEN goodssdata=goodsswithblue
      IF (bluename[blob] EQ 'F475W') THEN goodssdata=goodsswithalt 
   ; make new 3D-HST catalogs with color cuts
      xcolorhst = goodssdata.bluemag - goodssdata.midmag
      ycolorhst = goodssdata.midmag - goodssdata.redmag
      hstcuts = colorsort(xcolorhst, ycolorhst, m[blob], b[blob], cap[blob])
      goodssdatacuts = {ID:goodssdata.ID[hstcuts], ra:goodssdata.ra[hstcuts], dec:goodssdata.dec[hstcuts], $
         redmag:goodssdata.redmag[hstcuts], midmag:goodssdata.midmag[hstcuts], bluemag:goodssdata.bluemag[hstcuts], z:goodssdata.z[hstcuts]}

   ; establish info about the images
      bpm = mrdfits(mask[blob],0,hdr)    ; this is the mask image for weeding out bad pixels in field 
      sizeinfo = size(bpm)               ; the dimensions of the mask image in pixels 
      npix_x = sizeinfo[1]               ; number of pixels in mask in the x direction
      npix_y = sizeinfo[2]               ; number of pixels in mask in the y direction 

   ; make array of info about aperture radii - will contain overdensity factor & error on overdensity factor at each magnitude for each band
      radiusresults = dblarr(n_radii, nbins, 6., 2.)            ; blob compared to blob field 
      hstradiusresults = dblarr(n_radii, nbins, 6., 2.)         ; blob compared to 3D-HST field 


  ; now run analysis of density for different aperture radii 

   FOR r=0, n_radii-1. DO BEGIN 

    ; The following is split into parts to make things cleaner and easier to understand.



   ;PART 1: COLOR-COLOR STUFF
     ;
     ; First get raw blob density, ie, number of galaxies inside aperture centered on blob with no mag bins yet
     blob_galaxies_all=get_galaxies(blobra[blob], blobdec[blob], ap_radius[r], blobdata)
     ; get the number of galaxies in the aperture AFTER color cuts
     blob_galaxies_all_cuts=get_galaxies(blobra[blob], blobdec[blob], ap_radius[r], blobdatacuts)
     ; make region files showing blob galaxies before and after color cuts  
     IF (savereg EQ 'y') THEN BEGIN
       makereg, blobdata, blob_galaxies_all, blobname[blob], ap_radius_arcsec[r], 'n'            ; before cuts 
       makereg, blobdatacuts, blob_galaxies_all_cuts, blobname[blob], ap_radius_arcsec[r], 'y'    ; after cuts 
     ENDIF   ; user wants to save anything
     ;
     ; make color-color diagram: 
     field_galaxies = get_galaxies_noblob(blobra[blob], blobdec[blob], ap_radius[r], blobdata)   ; get IDs of field galaxies 
       ; set up plot window 
         !P.MULTI = [0,1,1,0,0] 
         loadct, 0
;         window, 1, retain=2, xsize=400, ysize=350
;     colordiagram, blob_galaxies_all, field_galaxies, datared, datamid, datablue, redname[blob], midname[blob], $  ;************************************************************************************fix this when colordiagram procedure is updated*********************************************
;                   bluename[blob], m[blob], b[blob], cap[blob], blobname[blob], ap_radius_arcsec[r] 
;     namestring = titlemaker('colors', blobname[blob], radius=ap_radius_arcsec[r], filetype='.png')
;     plotsave, saveanyplots, namestring



   ;PART 2: BLOB DENSITY 
     ;
     ; now get raw, unbinned blob galaxy count and print as a check
     blob_ngalraw = n_elements(blob_galaxies_all)
     print, blob_ngalraw, " galaxies in the blob region for "+blobname[blob]   
     ; now get raw, unbinned blob density in number of galaxies per square arcsecond (as a check) 
     blobdensityraw = float(blob_ngalraw) / (!dPI*ap_radius[r]^2.)    
     print, blobdensityraw, " galaxies per square degree in blob region for "+blobname[blob]  
     ;
     ; MAG BINNING: 
     ; make arrays containing the total number of blob galaxies in each bin in each band and the error on that number
     blob_ngal_binned = dblarr(nbins, 6.)
     blob_ngalerr_binned = dblarr(nbins, 6.) 
     ; create blobmaglist  - will have six entries, all lists (one for each band; each entry is a list of IDs of galaxies in the blob in each mag bin)
     galbinning, brightestmag, dimmestmag, binsize, blobra[blob], blobdec[blob], ap_radius[r], blobdata, blobdatacuts, blobmaglist
     ; fill in the ngal arrays using the mag list 
     get_ngal_binned, nbins, blobmaglist, blob_ngal_binned, blob_ngalerr_binned, 6 
     ; test to make sure it's working - this should match the previous no-bin number
     print, total(blob_ngal_binned[*,0]), " total galaxies"
     print, total(blob_ngal_binned[*,3]), ' total galaxies after color cuts'
     ; finally, get density and error on density 
     densitycalc, blob_ngal_binned, blob_ngalerr_binned, ap_radius[r], densities, density_errors



   ;PART 3: SAMPLING THE FIELD with random apertures 
     ;
     print, 'Now sampling the field for overdensity analysis.' 
     ;
     ; set up the arrays which will contain all the information about the random apertures placed
     ap_coords = dblarr(2., nApertures)     ; this will contain RAs and declinations of all apertures
     ap_ngal_binned = dblarr(nbins,6.,nApertures)   ; this will contain # galaxies in each aperture in each mag bin in each band 
     ap_ngalerr_binned = dblarr(nbins,6.,nApertures)   ; same as above but errors on # rather than just #
     megalist = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin in each band
     ; same as above but for HST field: 
     hstap_coords = dblarr(2., nApertures)     ; this will contain RAs and declinations of all apertures
     hstap_ngal_binned = dblarr(nbins,6.,nApertures)   ; this will contain # galaxies in each aperture in each mag bin in each band 
     hstap_ngalerr_binned = dblarr(nbins,6.,nApertures)   ; same as above but errors on # rather than just #
     hstmegalist = list()   ; list of lists, each corresponding to one aperture, containing the ID #s of all the galaxies in each mag bin in each band
     ;
     ; make an aperture image where everything inside the aperture has a value of 1 and everything outside (just the edge stuff) has a value of 0 
     dist_ellipse, dim, [ 2.*ap_radius_pix[r], 2.*ap_radius_pix[r] ], ap_radius_pix[r], ap_radius_pix[r], 1.,0.
       dim[where(dim LT ap_radius_pix[r], ntotalpix)] = 1.           ; note that ntotalpix is defined in this line as well
       dim[where(dim GE ap_radius_pix[r])] = 0.
     ; same as above but for HST field: 
     dist_ellipse, hstdim, [ 2.*hstap_radius_pix[r], 2.*hstap_radius_pix[r] ], hstap_radius_pix[r], hstap_radius_pix[r], 1.,0.
       hstdim[where(hstdim LT ap_radius_pix[r], hstntotalpix)] = 1.           ; note that hstntotalpix is defined in this line as well
       hstdim[where(hstdim GE ap_radius_pix[r])] = 0.
     ; keep track of how many apertures are made in total
     ap_count = 0        ; for the blob field  
     hstap_count = 0     ; for the 3D-HST field
     ;
     ; now place random apertures and count galaxies inside them 
     ; FOR THE BLOB FIELD 
     print, 'BLOB FIELD ANALYSIS'
       TIC        ; time this to ensure that it's efficient 
       FOR ap=0, nApertures-1 DO BEGIN
         ; find a good aperture to use: 
            apfinder, npix_x, npix_y, ap_radius_pix[r], ap_radius[r], dim, bpm, hdr, ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec, blobra=blobra[blob], blobdec=blobdec[blob]   
         ; calculate fraction of aperture that is good pixels:
            ap_area_weight = (float(ntotalpix) - float(nbadpix)) / float(ntotalpix) 
         ; convert pixels to sky coords and add these coords to arrays of aperture properties
            xyad, hdr, ap_x, ap_y, ap_ra, ap_dec
            ap_coords[0,ap] = ap_ra
            ap_coords[1,ap] = ap_dec
         ; create apmaglist  - will have six entries, all lists (one for each band; each entry is a list of IDs of galaxies in the blob in each mag bin)
            galbinning, brightestmag, dimmestmag, binsize, ap_ra, ap_dec, ap_radius[r], blobdata, blobdatacuts, apmaglist       
         ; put the list of ID numbers per mag bin in this aperture into the mega list for all apertures:
            megalist.Add, apmaglist   
         ; use apmaglist to fill in the ngals arrays
            get_ngal_binned, nbins, apmaglist, ap_ngal_binned, ap_ngalerr_binned, 6, ap=ap
       ENDFOR    ; all the random apertures 
       print, ap_count ; How many total apertures were generated? 
       thrownout = (float(ap_count) - float(nApertures))/float(ap_count)*100.     ; percentage of generated apertures that ended up being unsuitable for analysis 
       TOC             ; How long did this analysis take? 
     ; FOR THE 3D-HST FIELD 
     print, 'GOODS-S FIELD ANALYSIS'
       TIC        ; time this to ensure that it's efficient 
       FOR ap=0, nApertures-1 DO BEGIN
         ; find a good aperture to use: 
            apfinder, hstnpix_x, hstnpix_y, hstap_radius_pix[r], ap_radius[r], hstdim, hstbpm, hsthdr, hstap_count, hstap_x, hstap_y, hstnbadpix, hstntotalpix, hstap_ra, hstap_dec 
         ; calculate fraction of aperture that is good pixels:
            hstap_area_weight = (float(hstntotalpix) - float(hstnbadpix)) / float(hstntotalpix) 
         ; convert pixels to sky coords and add these coords to arrays of aperture properties
            xyad, hsthdr, hstap_x, hstap_y, hstap_ra, hstap_dec    ; for HST field
            hstap_coords[0,ap] = hstap_ra
            hstap_coords[1,ap] = hstap_dec
         ; create apmaglist  - will have six entries, all lists (one for each band; each entry is a list of IDs of galaxies in the aperture in each mag bin)
            galbinning, brightestmag, dimmestmag, binsize, hstap_ra, hstap_dec, ap_radius[r], goodssdata, goodssdatacuts, hstapmaglist    
         ; put the list of ID numbers per mag bin in this aperture into the mega list for all apertures:
            hstmegalist.Add, hstapmaglist  
         ; use apmaglist to fill in the ngals arrays
            get_ngal_binned, nbins, hstapmaglist, hstap_ngal_binned, hstap_ngalerr_binned, 6, ap=ap
       ENDFOR    ; all the random apertures 
       print, hstap_count ; How many total apertures were generated in HST field? 
       hstthrownout = (float(hstap_count) - float(nApertures))/float(hstap_count)*100.     ; percentage of generated apertures that ended up being unsuitable for analysis 
       TOC             ; How long did this analysis take? 
     ; as a safety check, see how many apertures ended up being unsuitable in each analysis 
     print, strcompress(string(thrownout)), '% of blob-field apertures thrown out'
     print, strcompress(string(hstthrownout)), '% of GOODS-S apertures thrown out' 
     ;
     ; calculate densities in each bin in each aperture and errors: 
        densitycalc, ap_ngal_binned, ap_ngalerr_binned, ap_radius[r], ap_densities, ap_density_errors, weight=ap_area_weight
        densitycalc, hstap_ngal_binned, hstap_ngalerr_binned, ap_radius[r], hstap_densities, hstap_density_errors, weight=hstap_area_weight
     ; then average the densities of all the apertures in each bin to get average field density in each bin (for each band):
        field_densities = total(ap_densities,3)/double(nApertures)             ; field density in each bin in each band 
        field_density_errors = total(ap_density_errors,3)/double(nApertures)   ; errors in field density in each bin in each band 
        hstfield_densities = total(hstap_densities,3)/double(nApertures)             ; for HST 
        hstfield_density_errors = total(hstap_density_errors,3)/double(nApertures)   ; for HST 
     ; calculate the overdensity factor and its error for analysis of how aperture size affects the measured overdensity
        radiusresults[r,*,*,0] = densities/field_densities            ; overdensity factor wrt blob field 
        radiusresults[r,*,*,1] = density_errors/field_densities       ; error on overdensity factor wrt blob field 
        hstradiusresults[r,*,*,0] = densities/hstfield_densities            ; overdensity factor wrt 3D-HST field 
        hstradiusresults[r,*,*,1] = density_errors/hstfield_densities       ; error on overdensity factor wrt 3D-HST field 
 

     !P.MULTI = [0,1,1,0,1]
     window, 1, retain=2, xsize=1200, ysize=1000
     densityplotsmall, nbins, mags, hstfield_densities[*,0], hstfield_density_errors[*,0], densities[*,0], density_errors[*,0], blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob], 'n'
     namestring = strcompress(string(blobname[blob]))+'NMsymposium_overdensity_raw.png'
     plotsave, saveanyplots, namestring
     densityplotsmall, nbins, mags, hstfield_densities[*,3], hstfield_density_errors[*,3], densities[*,3], density_errors[*,3], blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob], 'y'
     namestring = strcompress(string(blobname[blob]))+'NMsymposium_overdensity_cuts.png'
     plotsave, saveanyplots, namestring

proceed = 'n'
IF (proceed NE 'n') THEN BEGIN
   ;PART 4: COMPARING THE BLOB AND THE FIELD 
     ;
     ; MAIN OVERDENSITY PLOTS:
     !P.MULTI = [0,2,3,0,1]
     window, 1, retain=2, xsize=950, ysize=1200
     ; RED BAND, NO CUTS
     densityplot, nbins, mags, field_densities[*,0], field_density_errors[*,0], densities[*,0], density_errors[*,0], hstfield_densities[*,0], hstfield_density_errors[*,0], blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob], 'n'
     ; MIDDLE BAND, NO CUTS 
     densityplot, nbins, mags, field_densities[*,1], field_density_errors[*,1], densities[*,1], density_errors[*,1], hstfield_densities[*,1], hstfield_density_errors[*,1], blobname[blob], nApertures, ap_radius_arcsec[r], midname[blob], 'n'
     ; BLUE BAND, NO CUTS 
     densityplot, nbins, mags, field_densities[*,2], field_density_errors[*,2], densities[*,2], density_errors[*,2], hstfield_densities[*,2], hstfield_density_errors[*,2], blobname[blob], nApertures, ap_radius_arcsec[r], bluename[blob], 'n'
     ; RED BAND WITH CUTS
     densityplot, nbins, mags, field_densities[*,3], field_density_errors[*,3], densities[*,3], density_errors[*,3], hstfield_densities[*,3], hstfield_density_errors[*,3], blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob], 'y'
     ; MIDDLE BAND WITH CUTS 
     densityplot, nbins, mags, field_densities[*,4], field_density_errors[*,4], densities[*,4], density_errors[*,4], hstfield_densities[*,4], hstfield_density_errors[*,4], blobname[blob], nApertures, ap_radius_arcsec[r], midname[blob], 'y'
     ; BLUE BAND WITH CUTS 
     densityplot, nbins, mags, field_densities[*,5], field_density_errors[*,5], densities[*,5], density_errors[*,5], hstfield_densities[*,5], hstfield_density_errors[*,5], blobname[blob], nApertures, ap_radius_arcsec[r], bluename[blob], 'y'
     ; save the window of plots 
     namestring = titlemaker('overdensity', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], filetype='.png')
     plotsave, saveanyplots, namestring
ENDIF           ; user wants to make these plots 


;***************************************************************************************PICK UP HERE*****************************************************************************************************************************************************************


;     ;
;     ; STATISTICAL LUMINOSITY FUNCTIONS 
;     ; combine everything from before to make these - see statlum procedure header 
;     !P.MULTI = [0,1,3,0,1]
;     window, 1, retain=2, xsize=650, ysize=1200
;     ; RED:
;     statlum, mags, binsize, blob_ngal_binned[*,0], blob_ngal_binned[*,3], densities[*,0], field_densities[*,0], ap_radius[r], blobname[blob], nApertures, ap_radius_arcsec[r], redname[blob]  
;     ; MIDDLE: 
;     statlum, mags, binsize, blob_ngal_binned[*,1], blob_ngal_binned[*,4], densities[*,1], field_densities[*,1], ap_radius[r], blobname[blob], nApertures, ap_radius_arcsec[r], midname[blob]  
;     ; BLUE: 
;     statlum, mags, binsize, blob_ngal_binned[*,2], blob_ngal_binned[*,5], densities[*,2], field_densities[*,2], ap_radius[r], blobname[blob], nApertures, ap_radius_arcsec[r], bluename[blob] 
;     ; save the window of plots 
;     namestring = titlemaker('statlum', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], filetype='.png') 
;     plotsave, saveanyplots, namestring
;
;     ;
;     ; HISTOGRAMS: How likely is it to find overdensities in this field anyway? 
;     ap_numbers_corrected = ap_densities*!dPI*ap_radius[r]^2.    ; calculate the number of galaxies in each random aperture corrected for bad pixels 
;     ap_n = ap_numbers_corrected[*,0,*]        ; corrected number of galaxies in each aperture without cuts (using the red band)
;     ap_n_cuts = ap_numbers_corrected[*,3,*]   ; corrected number of galaxies in each aperture with cuts (using the red band) 
;     ap_totals = total(ap_n,1)              ; total number of galaxies in each aperture regardless of magnitude
;     ap_totals_cuts = total(ap_n_cuts,1)    ; total number of galaxies in each aperture regardless of magnitude
     ;
     ; 3DHST HISTOGRAMS: How likely is it to find overdensities in this field anyway? 
         ; NOTE: dimensions of ap_numbers_corrected are [mag bin, band, aperture]!!!!!!!!!!!!!!!!!!!!!!!!
     hstap_numbers_corrected = hstap_densities*!dPI*ap_radius[r]^2.    ; calculate the number of galaxies in each random aperture corrected for bad pixels 
     hstap_n = hstap_numbers_corrected[*,0,*]        ; corrected number of galaxies in each aperture without cuts (using the red band)
     hstap_n_cuts = hstap_numbers_corrected[*,3,*]   ; corrected number of galaxies in each aperture with cuts (using the red band) 
     hstap_totals = total(hstap_n,1)              ; total number of galaxies in each aperture regardless of magnitude
     hstap_totals_cuts = total(hstap_n_cuts,1)    ; total number of galaxies in each aperture regardless of magnitude
     ; set up a new plot window
     !P.MULTI = [0,2,1,0,1]
     window, 1, retain=2, xsize=1200, ysize=400
     ; compare total numbers (regardless of magnitude) to blob, no cuts
     makehist, hstap_totals, total(blob_ngal_binned[*,0]), blobname[blob], nApertures, ap_radius_arcsec[r], 'n'  
     ; compare total numbers (regardless of magnitude) to blob WITH cuts
     makehist, hstap_totals_cuts, total(blob_ngal_binned[*,3]), blobname[blob], nApertures, ap_radius_arcsec[r], 'y'  
     ; save the window of plots 
     namestring = titlemaker('aphist', blobname[blob], nApertures=nApertures, radius=ap_radius_arcsec[r], filetype='.png')
     plotsave, saveanyplots, namestring

   ENDFOR   ; all aperture sizes

;
;   ;PART 5: ANALYSIS OF APERTURE SIZE IMPACT ON MEASURED OVERDENSITY 
;     ;
;   ; if there were multiple aperture radii being tested, analyze them with plots:  
;   IF (n_radii GT 1.) THEN BEGIN 
;     ;
;     ; now make plots of the results of the aperture size analysis: overdensity factor vs magnitude for each band for each radius size 
;     loadct, 39    ; new color table 
;     !P.MULTI = [0,2,3,0,1]
;     window, 1, retain=2, xsize=950, ysize=1200
;     ; RED BAND, NO CUTS
;     Rmultirad, ap_radius_arcsec, mags, radiusresults, 0, blobname[blob], nApertures, redname[blob], 'n'
;     ; MIDDLE BAND, NO CUTS
;     Rmultirad, ap_radius_arcsec, mags, radiusresults, 1, blobname[blob], nApertures, midname[blob], 'n' 
;     ; BLUE BAND, NO CUTS
;     Rmultirad, ap_radius_arcsec, mags, radiusresults, 2, blobname[blob], nApertures, bluename[blob], 'n' 
;     ; RED BAND WITH CUTS
;     Rmultirad, ap_radius_arcsec, mags, radiusresults, 3, blobname[blob], nApertures, redname[blob], 'y' 
;     ; MIDDLE BAND WITH CUTS
;     Rmultirad, ap_radius_arcsec, mags, radiusresults, 4, blobname[blob], nApertures, midname[blob], 'y' 
;     ; BLUE BAND WITH CUTS
;     Rmultirad, ap_radius_arcsec, mags, radiusresults, 5, blobname[blob], nApertures, bluename[blob], 'y' 
;     ; save the window of plots 
;     namestring = titlemaker('Rmultirad', blobname[blob], nApertures=nApertures, filetype='.png')
;     plotsave, saveanyplots, namestring
;     ;
;     ; more plots of the results of the aperture size analysis: overdensity factor vs aperture radius for each band for each mag bin 
;     window, 1, retain=2, xsize=950, ysize=1200
;     ; RED BAND, NO CUTS
;     Rmultimag, mags, ap_radius_arcsec, radiusresults, 0, blobname[blob], nApertures, redname[blob], 'n' 
;     ; MIDDLE BAND, NO CUTS
;     Rmultimag, mags, ap_radius_arcsec, radiusresults, 1, blobname[blob], nApertures, midname[blob], 'n' 
;     ; BLUE BAND, NO CUTS
;     Rmultimag, mags, ap_radius_arcsec, radiusresults, 2, blobname[blob], nApertures, bluename[blob], 'n' 
;     ; RED BAND WITH CUTS
;     Rmultimag, mags, ap_radius_arcsec, radiusresults, 3, blobname[blob], nApertures, redname[blob], 'y' 
;     ; MIDDLE BAND WITH CUTS
;     Rmultimag, mags, ap_radius_arcsec, radiusresults, 4, blobname[blob], nApertures, midname[blob], 'y' 
;     ; BLUE BAND WITH CUTS
;     Rmultimag, mags, ap_radius_arcsec, radiusresults, 5, blobname[blob], nApertures, bluename[blob], 'y' 
;     ; save the window of plots 
;     namestring = titlemaker('Rmultimag', blobname[blob], nApertures=nApertures, filetype='.png')
;     plotsave, saveanyplots, namestring
;     ;
;   ENDIF     ; there are multiple aperture radii being tested 


 ENDFOR ;all the blobs


END             ; end of density procedure









;----------------------------------FINDING DENSITIES----------------------------------------------------------------------------------------------------------------------;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO galbinning, brightestmag, dimmestmag, binsize, ra, dec, radius, data, datacuts, bigmaglist                                                                         ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; galbinning procedure                                                                                                                            ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTE: Uses get_galaxies_binned function.                                                                                                        ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION get_galaxies_binned                                                                                                                                  ;;
    bigmaglist = list()                                                                                                                                                  ;;
    FOR i = 3, 5 DO BEGIN                                                                                                                                                ;;
      maglist=list()                                                                                                                                                     ;;
      FOR mag=brightestmag, dimmestmag, binsize DO BEGIN                                                                                                                 ;;
        galaxies_binned = get_galaxies_binned(ra, dec, radius, data, mag, mag+binsize, data.(i))                                                                         ;;
        maglist.Add, galaxies_binned                                                                                                                                     ;;
      ENDFOR                                                                                                                                                             ;;
      bigmaglist.Add, maglist                                                                                                                                            ;;
    ENDFOR                                                                                                                                                               ;;
    FOR i = 3, 5 DO BEGIN                                                                                                                                                ;;
      maglist=list()                                                                                                                                                     ;;
      FOR mag=brightestmag, dimmestmag, binsize DO BEGIN                                                                                                                 ;;
        galaxies_binned = get_galaxies_binned(ra, dec, radius, datacuts, mag, mag+binsize, datacuts.(i))                                                                 ;;
        maglist.Add, galaxies_binned                                                                                                                                     ;;
      ENDFOR                                                                                                                                                             ;;
      bigmaglist.Add, maglist                                                                                                                                            ;;
    ENDFOR                                                                                                                                                               ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO get_ngal_binned, nbins, bigmaglist, ngal_binned, ngalerr_binned, nbands, ap=ap                                                                                             ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; get_ngal_binned procedure                                                                                                                       ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTE: Uses Poisson_error function.                                                                                                              ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  FORWARD_FUNCTION Poisson_error                                                                                                                                         ;;
    FOR band= 0, nbands-1 DO BEGIN                                                                                                                                              ;;
      maglist = bigmaglist[band]                                                                                                                                         ;;
      FOR i=0, nbins-1. DO BEGIN                                                                                                                                         ;;
        bin = maglist[i]                                                                                                                                                 ;;
        IF (bin[0] EQ -1) THEN ngal = 0. ELSE ngal = n_elements(bin)                                                                                                     ;;
        IF (n_elements(ap) EQ 0) THEN BEGIN                                                                                                                              ;;
          ngal_binned[i,band] = ngal                                                                                                                                     ;;
          ngalerr_binned[i,band] = Poisson_error(ngal)                                                                                                                   ;;
        ENDIF ELSE IF (n_elements(ap) NE 0) THEN BEGIN                                                                                                                   ;;
          ngal_binned[i,band,ap] = ngal                                                                                                                                  ;;
          ngalerr_binned[i,band,ap] = Poisson_error(ngal)                                                                                                                ;;
        ENDIF                                                                                                                                                            ;;
      ENDFOR   ; all the bins                                                                                                                                             ;
    ENDFOR   ; all the bands                                                                                                                                              ;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  FUNCTION get_galaxies, ra, dec, radius, data                                                                                                                           ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; get_galaxies function                                                                                                                           ;                     ;
  ; INPUTS: ra, dec - the RA and declination of a point at which to place an aperture                                                               ;                     ;
  ;         radius - the radius of the aperture placed at the given coordinates                                                                     ;                     ;
  ;         data - a SExtractor catalog of sources which contains each source's RA, declination, a flag indicating what type of object it is,       ;                     ;
  ;                  and a flag indicating whether or not there are any bad pixels in the source in the original image                              ;                     ;
  ; OUTPUT: galaxies_in_aperture - the index numbers of the galaxies in the data catalog which are located within the specified aperture            ;                     ;
  ; NOTES: The output is found by determining the distance between the aperture's center and each source, then requiring that distance to be less   ;                     ;
  ;           than the radius of the aperture, filtering out all objects which are marked as stars in the catalog, and filtering out all objects    ;                     ;
  ;           with flags that indicate that they lie in areas of bad pixels. This ensures that only real sources that are galaxies are returned.    ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    distance = sqrt(  ( (ra - data.ra) * cos(dec*!dPI/180.) )^2. + (dec - data.dec)^2. )         ; Pythagorean theorem                                                    ;
    ; the source has to be contained within the aperture, has to be a galaxy (no stars), and has to have good data (no bad pixels)                                        ;
    galaxies_in_aperture = WHERE ((distance LT radius), n_galaxies)                                                                                                      ;;
    RETURN, galaxies_in_aperture                                                                                                                                         ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  FUNCTION get_galaxies_noblob, ra, dec, radius, data                                                                                                                    ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; get_galaxies_noblob function                                                                                                                    ;                     ;
  ; INPUTS: ra, dec - the RA and declination of a point at which to place an aperture                                                               ;                     ;
  ;         radius - the radius of the aperture placed at the given coordinates                                                                     ;                     ;
  ;         data - a SExtractor catalog of sources which contains each source's RA, declination, a flag indicating what type of object it is,       ;                     ;
  ;                  and a flag indicating whether or not there are any bad pixels in the source in the original image                              ;                     ;
  ; OUTPUT: galaxies_outside_aperture - the index numbers of the galaxies in the data catalog which are located OUTSIDE the specified aperture      ;                     ;
  ; NOTES: The output is found by determining the distance between the aperture's center and each source, then requiring that distance to be greater;                     ;
  ;           than the radius of the aperture, filtering out all objects which are marked as stars in the catalog, and filtering out all objects    ;                     ;
  ;           with flags that indicate that they lie in areas of bad pixels. This ensures that only real sources that are galaxies are returned.    ;                     ;
  ;        This function essentially does the exact opposite of what get_galaxies does.                                                             ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    distance = sqrt(  ( (ra - data.ra) * cos(dec*!dPI/180.) )^2. + (dec - data.dec)^2. )         ; Pythagorean theorem                                                    ;
    ; the source has to be located outside the aperture                                                                                                                   ;
    galaxies_outside_aperture = WHERE ((distance GE radius), n_galaxies)                                                                                                 ;;
    RETURN, galaxies_outside_aperture                                                                                                                                    ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  FUNCTION get_galaxies_binned, ra, dec, radius, data, brightmag, dimmag, datamags                                                                                       ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; get_galaxies_binned function                                                                                                                    ;                     ;
  ; INPUTS: ra, dec - the RA and declination of a point at which to place an aperture                                                               ;                     ;
  ;         radius - the radius of the aperture placed at the given coordinates                                                                     ;                     ;
  ;         data - a SExtractor catalog of sources which contains each source's RA, declination, magnitude, a flag indicating what type of object   ;                     ;
  ;                  it is, and a flag indicating whether or not there are any bad pixels in the source in the original image                       ;                     ;
  ;         brightmag - the lower bound (brightest magnitude) of a magnitude range (bin)                                                            ;                     ;
  ;         dimmag - the upper bound (faintest magnitude) of a magnitude range (bin)                                                                ;                     ;
  ; OUTPUT: galaxies_in_aperture - the index numbers of the galaxies in the data catalog which are located within the specified aperture and have   ;                     ;
  ;                                  magnitudes between brightmag and dimmag                                                                        ;                     ;
  ; NOTES: The output is found by determining the distance between the aperture's center and each source, then requiring that distance to be less   ;                     ;
  ;           than the radius of the aperture, filtering out all objects which are marked as stars in the catalog, filtering out all objects        ;                     ;
  ;           with flags that indicate that they lie in areas of bad pixels, and filtering out all objects with magnitudes outside the range set    ;                     ;
  ;           by brightmag and dimmag. This ensures that only real sources that are galaxies with specific magnitudes are returned.                 ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    distance = sqrt(  ( (ra - data.ra) * cos(dec*!dPI/180.) )^2. + (dec - data.dec)^2. )         ; Pythagorean theorem                                                    ;
    ; the source:                                                                                                                                                         ;
    ;  - has to be contained within the aperture, and                                                                                                                     ;                                                                                                                   ;
    ;  - has to have a magnitude between brightmag and dimmag                                                                                                             ;
    galaxies_in_aperture = WHERE (((distance LT radius) AND (datamags GE brightmag) AND (datamags LT dimmag)), n_galaxies)                                               ;;                                                   ;;
    RETURN, galaxies_in_aperture                                                                                                                                         ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  FUNCTION Poisson_error, ngal                                                                                                                                           ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; Poisson_error function                                                                                                                          ;                     ;
  ; INPUT: ngal - a number (of galaxies) for which we want to determine uncertainty                                                                 ;                     ;
  ; OUTPUT: error - the Poisson error on the input number (the square root of the input: N = sqrt(S) )                                              ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    error = sqrt(double(ngal))                                                                                                                                           ;;
    RETURN, error                                                                                                                                                        ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO densitycalc, ngal_binned, ngalerr_binned, radius, densities, density_errors, weight=weight                                                                         ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; densitycalc procedure                                                                                                                           ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTE: Uses Poisson_error function.                                                                                                              ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    IF (n_elements(weight) EQ 0) THEN BEGIN                                                                                                                              ;;
      densities = double(ngal_binned)/(!dPI*radius^2.)                                                                                                                   ;;
      density_errors = double(ngalerr_binned)/(!dPI*radius^2.)                                                                                                           ;;
    ENDIF ELSE IF (n_elements(weight) NE 0) THEN BEGIN                                                                                                                   ;;
      densities = double(ngal_binned)/(!dPI*radius^2.*weight)                                                                                                            ;;
      density_errors = double(ngalerr_binned)/(!dPI*radius^2.*weight)                                                                                                    ;;
    ENDIF                                                                                                                                                                ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO apfinder, npix_x, npix_y, ap_radius_pix, ap_radius, dim, bpm, hdr, ap_count, ap_x, ap_y, nbadpix, ntotalpix, ap_ra, ap_dec, blobra=blobra,blobdec=blobdec          ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; apfinder procedure                                                                                                                              ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTE:                                                                                                                                           ;                     ; 
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    repeatflag = 0.       ; this will ensure that apertures don't get wasted on areas with only bad pixels                                                                ;
    ap_count++                                                                                                                                                           ;;
    WHILE (repeatflag EQ 0.) DO BEGIN                                                                                                                                    ;;
      ; make a random aperture:                                                                                                                                           ;
seed = !NULL
       ap_x = floor(ap_radius_pix + (double(npix_x)-2.*ap_radius_pix)*randomu(seed))                                                                                     ;;
seed = !NULL
       ap_y = floor(ap_radius_pix + (double(npix_y)-2.*ap_radius_pix)*randomu(seed))                                                                                     ;;
print, ap_x, ap_y
      ; FIND BAD PIXELS:                                                                                                                                                  ;
       ; find corners of box into which to insert aperture image                                                                                                          ; 
       xlo = ap_x - floor(ap_radius_pix)                                                                                                                                 ;; 
       xhi = ap_x + floor(ap_radius_pix)                                                                                                                                 ;;
       ylo = ap_y - floor(ap_radius_pix)                                                                                                                                 ;;  
       yhi = ap_y + floor(ap_radius_pix)                                                                                                                                 ;;
       dimsize = size(dim)                                                                                                                                               ;;
       IF ((xhi-xlo) NE (dimsize[1]-1.)) THEN BEGIN                                                                                                                      ;;
         xhi = ap_x + floor(ap_radius_pix) - 1.                                                                                                                          ;;
       ENDIF                                                                                                                                                             ;;
       IF ((yhi-ylo) NE (dimsize[2]-1.)) THEN BEGIN                                                                                                                      ;;
         yhi = ap_y + floor(ap_radius_pix) - 1.                                                                                                                          ;;
       ENDIF                                                                                                                                                             ;;
       ; cut out the piece of the bpm that contains this random aperture                                                                                                  ;
       minibpm = bpm[xlo:xhi,ylo:yhi]                                                                                                                                    ;;
       ; multiply the aperture map by the bpm so only BAD pixels INSIDE the aperture are left (since we want to count them):                                              ;
       apbpm = minibpm*dim                                                                                                                                               ;;
       badpix = where(apbpm GT 0, nbadpix)   ; label the bad pixels                                                                                                       ;
                                                                                                                                                                         ;;
      ; FIND VICINITY TO BLOB:                                                                                                                                            ;
       IF (n_elements(blobra) NE 0) THEN BEGIN                                                                                                                           ;;
         xyad, hdr,ap_x, ap_y, ap_ra, ap_dec      ; convert pixels into coordnates in sky                                                                                 ;
         ; find how close the aperture is to the blob:                                                                                                                    ;
         vicinity = sqrt(  ( (blobra - ap_ra) * cos(blobdec*!dPI/180.) )^2. + (blobdec - ap_dec)^2. )                                                                    ;;
       ENDIF ELSE vicinity = 3.*ap_radius   ; if we're looking in a field that doesn't contain the blob, it doesn't matter because there's no blob to avoid               ;
                                                                                                                                                                         ;;
      ; if the aperture has at least some good pixels AND isn't near the blob, then keep it; otherwise, throw it out and make a new aperture                              ;
      pixratio = float(nbadpix)/float(ntotalpix)
      IF ((nbadpix LT ntotalpix) AND (vicinity GT 2.*ap_radius)) THEN BEGIN      
        print, 'KEEPING aperture ', strcompress(string(ap_count)), ' with bad/total pixel ratio of ', strcompress(string(pixratio))                                                                                        ;;
        repeatflag = 1.                                                                                                                                                  ;; 
      ENDIF ELSE BEGIN                                                                                                                                                   ;;
        repeatflag = 0.                                                                                                                                                  ;;
        ; keep a record of why this aperture was thrown out, just to make sure nothing is going amiss                                                                     ;
;        pixratio = float(nbadpix)/float(ntotalpix)                                                                                                                       ;;
        print, 'aperture ', strcompress(string(ap_count)), ' was within ', strcompress(string(vicinity/ap_radius)), $                                                    ;;
          ' aperture radii of LAB with bad/total pixel ratio of ', strcompress(string(pixratio))                                                                         ;;
print, ' '
        ap_count++    ; note that another aperture is being made since this one failed to meet the criteria for a good aperture                                           ;
      ENDELSE                                                                                                                                                            ;;
    ENDWHILE                                                                                                                                                             ;;
       ; once the aperture passes the above WHILE test, it can be used                                                                                                    ;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------;








;----------------------------------MAKING PLOTS---------------------------------------------------------------------------------------------------------------------------;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO colordiagram, blobgals, fieldgals, datared, datamid, datablue, redname, midname, bluename, m, b, cap, blobname, radius                                             ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; colordiagram procedure                                                                                                                          ;                     ;
  ;   description                                                                                                                                   ;                     ;************
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker                                                                                                                                           ;;
   ; separate out the magnitudes to define colors                                                                                                                         ;
      redmags = datared[blobgals].MAG_ISO                                                                                                                                ;;
      midmags = datamid[blobgals].MAG_ISO                                                                                                                                ;;
      bluemags = datablue[blobgals].MAG_ISO                                                                                                                              ;;
      field_redmags = datared[fieldgals].MAG_ISO                                                                                                                         ;;
      field_midmags = datamid[fieldgals].MAG_ISO                                                                                                                         ;;
      field_bluemags = datablue[fieldgals].MAG_ISO                                                                                                                       ;;
   ; make the colors                                                                                                                                                      ;
      color_x = bluemags - midmags                                                                                                                                       ;;
      color_y = midmags - redmags                                                                                                                                        ;;
      field_color_x = field_bluemags - field_midmags                                                                                                                     ;;
      field_color_y = field_midmags - field_redmags                                                                                                                      ;;
   ; make the plot                                                                                                                                                        ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('colors', blobname, radius=radius)                                                                                                            ;;
        xtitle = bluename + ' - ' + midname                                                                                                                              ;;
        ytitle = midname + ' - ' + redname                                                                                                                               ;;
     plot, color_x, color_y, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, psym=8, xrange=[-5,5],/xstyle, yrange=[-5,5],/ystyle, charsize=1.5, $   ;; 
        thick=2, ymargin=[4,4]                                                                                                                                           ;;
     oplot, field_color_x, field_color_y, color=0, psym=3                                                                                                                ;;
     LEGEND, ['blob galaxy','field galaxy'], /left, /top, color=0, textcolor=0, psym=[8,3], charsize=1, charthick=1, /box, outline_color=0                               ;;
   ; add the cut line                                                                                                                                                     ;
      xvalues = 0.01*findgen(100000) - 500.                                                                                                                              ;;
      cut = xvalues*m + b                                                                                                                                                ;;
      cut[WHERE(cut GT cap)] = cap                                                                                                                                       ;;
      oplot, xvalues, cut, color=50, linestyle=2, thick=2                                                                                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO densityplot, nbins, mags, field_densities, field_density_errs, densities, density_errs, hst_densities, hst_density_errs, blobname, nApertures, radius, band, cuts  ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; densityplot procedure                                                                                                                           ;                     ;******make HST stuff look better
  ;   Makes plots of blob and field densities as a function of magnitude.                                                                           ;                     ;
  ; INPUTS: blobname - name of blob whose density is being plotted (for plot title)                                                                 ;                     ;
  ;         nbins - number of magnitude bins being plotted on the x-axis                                                                            ;                     ;
  ;         mags - array of magnitudes (data for the x-axis)                                                                                        ;                     ;
  ;         field_densities - array of field density measurements (data for the y-axis)                                                             ;                     ;
  ;         field_density_errors - array of errors on field density measurements (to be overplotted as a swath on top of field_densities)           ;                     ;
  ;         densities - array of blob density measurements (also data for the y-axis)                                                               ;                     ;
  ;         density_errors - array of errors on blob density measurements (to be overplotted on top of densities using errplot)                     ;                     ;
  ;         apsnameplot - a string showing the number of random apertures used in the analysis (for plot title)                                     ;                     ;
  ;         xtitle - a string containing the label for the x-axis of the plot                                                                       ;                     ;
  ;         cuts - a string, 'y' or 'n', showing whether the densities used are for raw data or post-color-cut data (for plot title)                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
     loadct, 0                                                                                                                                                           ;;
   ; set up vertices of polygon for polyfill                                                                                                                              ;
     xpoints = dblarr(2*nbins)                                                                                                                                           ;;
     xpoints[0:nbins-1] = mags                                                                                                                                           ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       xpoints[i+nbins] = mags[nbins-1-i]                                                                                                                                ;;
     ENDFOR                                                                                                                                                              ;;
     ; y is specific to each catalog (band)                                                                                                                               ;
     ypoints=dblarr(2*nbins)                                                                                                                                             ;;
     ypoints[0:nbins-1] = field_densities-field_density_errs                                                                                                             ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       ypoints[i+nbins] = field_densities[nbins-1-i]+field_density_errs[nbins-1-i]                                                                                       ;;
     ENDFOR                                                                                                                                                              ;;
   ; make the main plot comparing blob to field with error bars                                                                                                           ;
     ; set up titles for plot and axes                                                                                                                                    ;
        plottitle = titlemaker('overdensity', blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts)                                                      ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle='N!Igal!N per deg!E2!N per 0.5 mag'                                                                                                                       ;;
     plot, mags, densities, background=255, color=0, title=plottitle, xtitle=xtitle, ytitle=ytitle, psym=-8, /ylog, yrange=[1d3,5d6],/ystyle, $                          ;;
         xrange=[22.5,29.5],/xstyle, charsize=2., xthick=2, ythick=2, xmargin=[10,7], ymargin=[6,6]                                                                      ;;
     polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,5d6], /data, noclip=0                                                                                    ;;
     errplot, mags, densities-density_errs, densities+density_errs, color=0, thick=2                                                                                     ;;
     oplot, mags, densities, psym=-8, color=0, thick=2     ; have to do this again because polyfill covers it up                                                          ;
     oplot, mags, field_densities, color=0, linestyle=2, thick=2                                                                                                         ;;
     loadct, 39                                                                                                                                                          ;;
     oplot, mags, hst_densities, color=50, linestyle=5, thick=2                                                                                                          ;;
     errplot, mags, hst_densities-hst_density_errs, hst_densities+hst_density_errs, color=50, thick=2                                                                    ;;
     LEGEND, ['blob','blob field', 'GOODS-S field'], /left, /top, color=[0,0,50], textcolor=0, linestyle=[0,2,5], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5               ;;
     axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=2, charthick=2, xtickformat="(A1)", xthick=2                                                   ;;
     axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=2                                                                ;;
     axis, color=0, yaxis=0, /ylog, yrange=[1d3,5d6], /ystyle,  /ynozero, /data, charsize=2, charthick=2, ytickformat="(A1)", ythick=2                                   ;;
     axis, color=0, yaxis=1, /ylog, yrange=[1d3,5d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=2                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO densityplotsmall, nbins, mags, field_densities, field_density_errs, densities, density_errs, blobname, nApertures, radius, band, cuts                              ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; densityplot procedure                                                                                                                           ;                     ;
  ;   Makes plots of blob and field densities as a function of magnitude.                                                                           ;                     ;
  ; INPUTS: blobname - name of blob whose density is being plotted (for plot title)                                                                 ;                     ;
  ;         nbins - number of magnitude bins being plotted on the x-axis                                                                            ;                     ;
  ;         mags - array of magnitudes (data for the x-axis)                                                                                        ;                     ;
  ;         field_densities - array of field density measurements (data for the y-axis)                                                             ;                     ;
  ;         field_density_errors - array of errors on field density measurements (to be overplotted as a swath on top of field_densities)           ;                     ;
  ;         densities - array of blob density measurements (also data for the y-axis)                                                               ;                     ;
  ;         density_errors - array of errors on blob density measurements (to be overplotted on top of densities using errplot)                     ;                     ;
  ;         apsnameplot - a string showing the number of random apertures used in the analysis (for plot title)                                     ;                     ;
  ;         xtitle - a string containing the label for the x-axis of the plot                                                                       ;                     ;
  ;         cuts - a string, 'y' or 'n', showing whether the densities used are for raw data or post-color-cut data (for plot title)                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
     loadct, 0                                                                                                                                                           ;;
   ; set up vertices of polygon for polyfill                                                                                                                              ;
     xpoints = dblarr(2*nbins)                                                                                                                                           ;;
     xpoints[0:nbins-1] = mags                                                                                                                                           ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       xpoints[i+nbins] = mags[nbins-1-i]                                                                                                                                ;;
     ENDFOR                                                                                                                                                              ;;
     ; y is specific to each catalog (band)                                                                                                                               ;
     ypoints=dblarr(2*nbins)                                                                                                                                             ;;
     ypoints[0:nbins-1] = field_densities-field_density_errs                                                                                                             ;;
     FOR i=0, nbins-1 DO BEGIN                                                                                                                                           ;;
       ypoints[i+nbins] = field_densities[nbins-1-i]+field_density_errs[nbins-1-i]                                                                                       ;;
     ENDFOR                                                                                                                                                              ;;
   ; make the main plot comparing blob to field with error bars                                                                                                           ;
     ; set up titles for plot and axes                                                                                                                                    ;
        plottitle = titlemaker('overdensity', blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts)                                                      ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle='N!Igal!N per deg!E2!N per 0.5 mag'                                                                                                                       ;;
     plot, mags, densities, background=255, color=0, title=plottitle, xtitle=xtitle, ytitle=ytitle, psym=-8, /ylog, yrange=[1d3,1d6],/ystyle, $                          ;;
         xrange=[22.5,29.5],/xstyle, charsize=2.5, xthick=2, ythick=2, xmargin=[10,7], ymargin=[6,6], charthick=2                                                                      ;;
     polyfill, xpoints, ypoints, color=200, clip=[22.5,1d3,29.5,1d6], /data, noclip=0                                                                                    ;;
     errplot, mags, densities-density_errs, densities+density_errs, color=0, thick=4                                                                                     ;;
     oplot, mags, densities, psym=-8, color=0, thick=4     ; have to do this again because polyfill covers it up                                                          ;
     oplot, mags, field_densities, color=0, linestyle=2, thick=4                                                                                                         ;;                                                                  ;;
     LEGEND, ['blob','field (GOODS-S)'], /right, /top, color=[0,50], textcolor=0, linestyle=[0,2], thick=4., charsize=2.5, /box, outline_color=0., charthick=2, number=0.1            ;;
     axis, color=0, xaxis=0, xrange=[22.5,29.5], /xstyle, /data, charsize=4, charthick=2, xtickformat="(A1)", xthick=4                                                   ;;
     axis, color=0, xaxis=1, xrange=[22.5,29.5], /xstyle, /data, charsize=0, xtickformat="(A1)", xthick=4                                                                ;;
     axis, color=0, yaxis=0, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=4, charthick=2, ytickformat="(A1)", ythick=4                                   ;;
     axis, color=0, yaxis=1, /ylog, yrange=[1d3,1d6], /ystyle,  /ynozero, /data, charsize=0, ytickformat="(A1)", ythick=4                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;   
  PRO statlum, mags, binsize, blob_ngal_binned, blob_ngal_binned_cuts, densities, field_densities, ap_radius, blobname, nApertures, ap_radius_arcsec, band               ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; statlum procedure                                                                                                                               ;                     ;***********
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; define the statistical luminosity function                                                                                                                           ;
      overdensity = densities - field_densities                                                                                                                          ;;
      n_overdensity = overdensity * (!dPI*ap_radius^2.)  ; convert to number                                                                                              ;
   ; set up vertices of shaded area                                                                                                                                       ;
     xvertices = fltarr(2.0*n_elements(mags) + 2.)                                                                                                                       ;;
     yvertices = fltarr(2.0*n_elements(mags) + 2.)                                                                                                                       ;;
     xvertices[0] = mags[0]                                                                                                                                              ;;
     xvertices[-1] = mags[-1]                                                                                                                                            ;;
     yvertices[0] = 0                                                                                                                                                    ;;
     yvertices[-1] = 0                                                                                                                                                   ;;
     xvertices[1] = mags[0]                                                                                                                                              ;;
     xvertices[-2] = mags[-1]                                                                                                                                            ;;
     FOR i=1, n_elements(mags) -1 DO BEGIN                                                                                                                               ;;
       xvertices[2*i] = mags[i]-(binsize/2.0)                                                                                                                            ;;
       xvertices[2*i+1] = mags[i]-(binsize/2.0)                                                                                                                          ;;
     ENDFOR                                                                                                                                                              ;;
     FOR i=1, n_elements(mags) -1 DO BEGIN                                                                                                                               ;;
       yvertices[2*i+1] = blob_ngal_binned_cuts[i]                                                                                                                       ;;
       yvertices[2*i+2] = blob_ngal_binned_cuts[i]                                                                                                                       ;;
     ENDFOR                                                                                                                                                              ;;
   ; now plot statistical luminosity function                                                                                                                             ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('statlum', blobname, nApertures=nApertures, radius=ap_radius_arcsec, band=band)                                                               ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle = 'N!Igal!N'                                                                                                                                              ;;
     plot, mags, n_overdensity, background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, linestyle=2, yrange=[0,10],/ystyle, xrange=[22,30],/xstyle, $        ;;
       charsize=3, thick=2, xmargin=[10,7], ymargin=[4,4]                                                                                                                ;;
     polyfill, xvertices,yvertices, color=220, clip=[22,0,30,10], /data, noclip=0                                                                                        ;;
     oplot, mags, blob_ngal_binned, color=0, psym=10, thick=2                                                                                                            ;;
     oplot, mags, blob_ngal_binned_cuts, color=0, psym=10, thick=2                                                                                                       ;;
     oplot, mags, n_overdensity, color=0, thick=2, linestyle=2                                                                                                           ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO makehist, ap_totals, blob_ngal, blobname, nApertures, radius, cuts     ; maybe put optional band keyword in here or mag keyword or smth                             ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; makehist procedure                                                                                                                              ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;***************
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; set up titles for plot and axes                                                                                                                                      ;
      title = titlemaker('aphist', blobname, nApertures=nApertures, radius=radius, cuts=cuts)                                                                            ;;
      xtitle = 'n galaxies in aperture'                                                                                                                                  ;;
      ytitle='number of apertures'                                                                                                                                       ;;
   plothist, ap_totals, background=255, color=0, axiscolor=0, title=title, xtitle=xtitle, ytitle=ytitle, bin=1, xrange=[0,blob_ngal+10.],/xstyle, charsize=1.5, $        ;;
     thick=2, ymargin=[4,4]                                                                                                                                              ;;
   oplot, [blob_ngal, blob_ngal], [0., nApertures], linestyle=2, thick=2, color=0.                                                                                       ;;
   LEGEND, ['blob','field'], /center, /top, color=0, textcolor=0, linestyle=[0,2], thick=2., charsize=1, /box, outline_color=0.,number=0.1, charthick=1.5                ;;
  ; xyouts, blob_ngal+1., nApertures/10., blobgaltext, color=0, charsize=1.5                                                                                              ;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO Rmultirad, ap_radius_arcsec, mags, radiusresults, index, blobname, nApertures, band, cuts                                                                          ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; Rmultirad procedure                                                                                                                             ;                     ;
  ;   description                                                                                                                                   ;                     ;*************
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; set up an array of colors to plot each radius in                                                                                                                     ;
     n_radii = n_elements(ap_radius_arcsec)                                                                                                                              ;;
     indices = findgen(n_radii)                                                                                                                                          ;;
     colors = indices*254./n_radii                                                                                                                                       ;;
   ; now make the plot                                                                                                                                                    ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('Rmultirad', blobname, nApertures=nApertures, band=band, cuts=cuts)                                                                           ;;
        xtitle = 'magnitude (' + band + ')'                                                                                                                              ;;
        ytitle = 'overdensity factor'                                                                                                                                    ;;
     plot, mags, radiusresults[0,*,index,0], background=255, color=0, title=title, xtitle=xtitle, ytitle=ytitle, yrange=[0,50],/ystyle, charsize=2, thick=2, $           ;;
       ymargin=[4,4]                                                                                                                                                     ;;
     errplot, mags, radiusresults[0,*,index,0]-radiusresults[0,*,index,1], radiusresults[0,*,index,0]+radiusresults[0,*,index,1], color=0, thick=2                       ;;
     FOR r=1, n_radii-1. DO BEGIN                                                                                                                                        ;;
       oplot, mags, radiusresults[r,*,index,0], color = colors[r], thick=2                                                                                               ;;
       errplot, mags, radiusresults[r,*,index,0]-radiusresults[r,*,index,1], radiusresults[r,*,index,0]+radiusresults[r,*,index,1], color=colors[r], thick=2             ;;
     ENDFOR                                                                                                                                                              ;;
     LEGEND, strcompress(string(fix(ap_radius_arcsec)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., $            ;;
       number=0.5, charthick=1.5, thick=2                                                                                                                                ;;
  END                                                                                                                                                                    ;; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO Rmultimag, mags, ap_radius_arcsec, radiusresults, index, blobname, nApertures, band, cuts                                                                          ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; Rmultimag procedure                                                                                                                             ;                     ;
  ;   description                                                                                                                                   ;                     ;***********
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTES:                                                                                                                                          ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   ; set up an array of colors to plot each magnitude bin in                                                                                                              ;
     indices = findgen(n_elements(mags))                                                                                                                                 ;;
     colors = indices*254./float(n_elements(mags))                                                                                                                       ;;
   ; now make the plot                                                                                                                                                    ;
     ; set up titles for plot and axes                                                                                                                                    ;
        title = titlemaker('Rmultimag', blobname, nApertures=nApertures, band=band, cuts=cuts)                                                                           ;;
        xtitle='aperture radius (arcseconds)'                                                                                                                            ;;
        ytitle='overdensity factor'                                                                                                                                      ;;
     ; set up an empty plot (needs to have no data in case the first thing plotted is NaN):                                                                               ;
     plot, ap_radius_arcsec, radiusresults[*,1,index,0], background=255, color=colors[1], title=title, xtitle=xtitle, ytitle=ytitle, xrange=[0,30],/xstyle, $            ;;
       yrange=[0,100],/ystyle, charsize=2, ymargin=[4,4], /nodata                                                                                                        ;;
     ; now add data on top of the empty plot:                                                                                                                             ;
     FOR i=0, n_elements(mags)-1 DO BEGIN     ; need to make sure that NaNs are thrown out (because where field density is 0, we divided by 0)                            ;
       y = radiusresults[*,i,index,0]                                                                                                                                    ;;
       yerr = radiusresults[*,i,index,1]                                                                                                                                 ;;
       keep = WHERE(finite(y) EQ 1, nkeep)                                                                                                                               ;;
       IF (nkeep GT 0) THEN BEGIN                                                                                                                                        ;;
         oplot, ap_radius_arcsec(keep), y(keep), color=colors[i], thick=2                                                                                                ;;
         errplot, ap_radius_arcsec(keep), y(keep)-yerr(keep), y(keep)+yerr(keep), color=colors[i], thick=2                                                               ;;
       ENDIF                                                                                                                                                             ;;
     ENDFOR                                                                                                                                                              ;;
     ; add a legend:                                                                                                                                                      ;
     LEGEND, strcompress(string(fix(mags)), /remove), /right, /top, color=colors, linestyle=0, textcolor=0, charsize=1, /box, outline_color=0., number=0.5, $            ;;
        charthick=1.5, thick=2                                                                                                                                           ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO makereg, data, subset, blobname, radius, cuts                                                                                                                      ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; makereg procedure                                                                                                                               ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: etc                                                                                                                                     ;                     ;
  ;         etc                                                                                                                                     ;                     ;
  ; NOTE: Uses titlemaker function.                                                                                                                 ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
   FORWARD_FUNCTION titlemaker           ; this function is used in this procedure                                                                                        ;
   savethis = ''                                                                                                                                                         ;;
   READ, savethis, PROMPT='Make region file? (y/n)'                                                                                                                      ;;
   IF (savethis EQ 'y') THEN BEGIN                                                                                                                                       ;;
     cat = data[subset]                                                                                                                                                  ;;
     filename = titlemaker('blobgals', blobname, radius=radius, cuts=cuts, filetype='.reg')                                                                              ;;
     openw, 1, filename                                                                                                                                                  ;;
     FOR i=0, n_elements(cat)-1 DO BEGIN                                                                                                                                 ;;
       printf, 1,'J2000; circle ', cat[i].alpha_j2000, cat[i].delta_j2000, ' 10p '    ;#text={', cat[i].number,'}'                                                        ;
     ENDFOR                                                                                                                                                              ;;
     close, 1                                                                                                                                                            ;;
   ENDIF ELSE IF ((savethis NE 'y') AND (savethis NE 'n')) THEN BEGIN                                                                                                    ;;
     WHILE ((savethis NE 'y') AND (savethis NE 'n')) DO BEGIN                                                                                                            ;;
       print, 'Invalid response. Choose y or n.'                                                                                                                         ;;
       READ, savethis, PROMPT='Make region file? (y/n)'                                                                                                                  ;;
     ENDWHILE                                                                                                                                                            ;;
   ENDIF ELSE IF (savethis EQ 'n') THEN BEGIN                                                                                                                            ;;
        print, 'Region file not created.'                                                                                                                                ;;
   ENDIF                                                                                                                                                                 ;;
  END                                                                                                                                                                    ;;                                                                                                                                                                     ; 
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  PRO plotsave, saveanyplots, namestring                                                                                                                                 ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; plotsave procedure                                                                                                                              ;                     ;
  ;   Asks user whether or not to save a plot and saves said plot if the user chooses to do so.                                                     ;                     ;
  ; INPUTS: saveanyplots - a string, 'y' or 'n', determining if the user should be asked whether or not to save this plot at all (code will only    ;                     ;
  ;                          run if this is set to 'y')                                                                                             ;                     ;
  ;         namestring - a string containing the name of the file containing the plot if the user chooses to save it                                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    IF (saveanyplots EQ 'y') THEN BEGIN                                                                                                                                  ;;
      saveplot = ''                                                                                                                                                      ;;
      READ, saveplot, PROMPT='Save plot? (y/n)'                                                                                                                          ;;
      IF (saveplot EQ 'y') THEN BEGIN                                                                                                                                    ;;
         write_png, namestring, tvrd(/true)                                                                                                                              ;;
      ENDIF ELSE IF ((saveplot NE 'y') AND (saveplot NE 'n')) THEN BEGIN                                                                                                 ;;
        WHILE ((saveplot NE 'y') AND (saveplot NE 'n')) DO BEGIN                                                                                                         ;;
          print, 'Invalid response. Choose y or n.'                                                                                                                      ;;
          READ, saveplot, PROMPT='Save plot? (y/n)'                                                                                                                      ;;
        ENDWHILE                                                                                                                                                         ;;
      ENDIF ELSE IF (saveplot EQ 'n') THEN BEGIN                                                                                                                         ;;
        print, 'Plot not saved.'                                                                                                                                         ;;
      ENDIF                                                                                                                                                              ;;
    ENDIF   ; user wants to save at least some plots                                                                                                                      ;
  END                                                                                                                                                                    ;;                                                                                                                                                                   ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
  FUNCTION titlemaker, plotname, blobname, nApertures=nApertures, radius=radius, band=band, cuts=cuts, filetype=filetype                                                 ;;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
  ; titlemaker function                                                                                                                             ;                     ;
  ;   description                                                                                                                                   ;                     ;
  ; INPUTS: plotname - string containing the name of the kind of plot being made (ie, "overdensity," "color-color," etc)                            ;                     ;
  ;         blobname - string                                                                                                                       ;                     ;
  ;         nApertures - float                                                                                                                      ;                     ;
  ;         plotname - string containing the name of the kind of plot being made (ie, "overdensity," "color-color," etc)                            ;                     ;
  ;         radius - float                                                                                                                          ;                     ;
  ;         band - string (F140W, F606W, etc)                                                                                                       ;                     ;
  ;         cuts - string/keyword                                                                                                                   ;                     ;
  ;         filetype - specifies that this is going to be a filename if set; string that gives extension of the file                                ;                     ;
  ;-------------------------------------------------------------------------------------------------------------------------------------------------;                     ;
    ; first, make checks that see if the optional keywords are provided                                                                                                   ; 
    apcheck = n_elements(nApertures)                                                                                                                                     ;;
    radcheck = n_elements(radius)                                                                                                                                        ;;
    bandcheck = n_elements(band)                                                                                                                                         ;;
    cutscheck = n_elements(cuts)                                                                                                                                         ;;
    filecheck = n_elements(filetype)                                                                                                                                     ;;
                                                                                                                                                                         ;;
    ; next make empty strings set up for each element of a title                                                                                                          ;
    apsName = ''                                                                                                                                                         ;;                           
    radName = ''                                                                                                                                                         ;;
    bandName = ''                                                                                                                                                        ;;
    cutsName = ''                                                                                                                                                        ;;
    titleName = ''                                                                                                                                                       ;;
                                                                                                                                                                         ;;
    ; now make the final string that will eventually be returned                                                                                                          ;
    titlestring = ''                                                                                                                                                     ;;
                                                                                                                                                                         ;;
    ; since the cuts keyword can only be 'y' or 'n', if it's set, make sure it's one of those things                                                                      ;
    IF (cutscheck NE 0) THEN BEGIN                                                                                                                                       ;;
      IF ( (cuts NE 'y') AND (cuts NE 'n') ) THEN BEGIN                                                                                                                  ;;
        WHILE ( (cuts NE 'y') AND (cuts NE 'n') ) DO BEGIN        ; make sure the keyword is valid                                                                        ;
          print, 'Invalid cut keyword provided. Choose "y" or "n" only.'                                                                                                 ;;
          cuts = ''                                                                                                                                                      ;;
          READ, cuts, PROMPT='Color cuts?'                                                                                                                               ;;
        ENDWHILE                                                                                                                                                         ;;
      ENDIF ; cuts was set to something other than y or n                                                                                                                 ;
    ENDIF   ; the cuts keyword was provided at all                                                                                                                        ;
                                                                                                                                                                         ;;
    ; now, determine whether to make a plot title or a file name:                                                                                                         ;
                                                                                                                                                                         ;;
    IF (filecheck NE 0) THEN BEGIN            ; if the user wants a filename                                                                                              ;
      ;                                                                                                                                                                   ;
      titlestring = string(plotname) + '_' + string(blobname)                                                                                                            ;;
      ; now tack on any extra keywords one by one                                                                                                                         ;
      ;                                                                                                                                                                   ;
      ; start with nApertures                                                                                                                                             ;
      IF (apcheck NE 0) THEN BEGIN                                                                                                                                       ;; 
        apsName = strcompress(string(fix(nApertures)), /remove) + 'aps'                                                                                                  ;;
        titlestring = titlestring + '_' + string(apsName)                                                                                                                ;;
      ENDIF                                                                                                                                                              ;; 
      ;                                                                                                                                                                   ;                                                                                                                                                                   ;
      ; radius                                                                                                                                                            ;
      IF (radcheck NE 0) THEN BEGIN                                                                                                                                      ;; 
        radName = strcompress(string(fix(radius)), /remove) + 'arcsec'                                                                                                   ;;
        titlestring = titlestring + '_' + string(radName)                                                                                                                ;;
      ENDIF                                                                                                                                                              ;; 
      ;                                                                                                                                                                   ;
      ; band                                                                                                                                                              ;
      IF (bandcheck NE 0) THEN BEGIN                                                                                                                                     ;; 
        bandName = band                                                                                                                                                  ;;
        titlestring = titlestring + '_' + string(bandName)                                                                                                               ;;
      ENDIF                                                                                                                                                              ;;
      ;                                                                                                                                                                   ;
      ; cuts                                                                                                                                                              ;
      IF (cutscheck NE 0) THEN BEGIN                                                                                                                                     ;;
        IF (cuts EQ 'y') THEN BEGIN        ; make a string for the case that there are color cuts                                                                         ;
          cutsName = 'cuts'                                                                                                                                              ;;
        ENDIF ELSE IF (cuts EQ 'n') THEN BEGIN        ; make a string for data without cuts                                                                               ;
          cutsName = 'raw'                                                                                                                                               ;;
        ENDIF         ; now we've sorted out what string to use to describe the cuts keyword                                                                              ;
       titlestring = titlestring + '_' + string(cutsName)                                                                                                                ;;
      ENDIF                                                                                                                                                              ;;
      ;                                                                                                                                                                   ;
      ; finally, add the file extension to the end of the string                                                                                                          ;
      titlestring = titlestring + string(filetype)                                                                                                                       ;;
      ; now we have a valid file name                                                                                                                                     ;
                                                                                                                                                                         ;;
    ENDIF ELSE IF (filecheck EQ 0) THEN BEGIN         ; if the user is making a title for a plot                                                                          ;
      ;                                                                                                                                                                   ;
      ; generate a basic plot title based on what's being plotted                                                                                                         ;
      IF (plotname EQ 'colors') THEN titleName = 'Color-Color Diagram for '                                                                                              ;;
      IF (plotname EQ 'overdensity') THEN titleName = 'Galaxy Overdensity in Blob Region for '                                                                           ;;
      IF (plotname EQ 'statlum') THEN titleName = 'Statistical Luminosity Function for Galaxies in '                                                                     ;;
      IF (plotname EQ 'aphist') THEN titleName = 'Number of Galaxies in Apertures for '                                                                                  ;;
      IF (plotname EQ 'Rmultirad') THEN titleName = 'Galaxy Overdensity in Different Radii for '                                                                         ;;
      IF (plotname EQ 'Rmultimag') THEN titleName = 'Dependence of Overdensity Factor on Aperture Radius for '                                                           ;;
      titlestring = string(titleName) + string(blobname)                                                                                                                 ;;
      ;                                                                                                                                                                   ;
      ; Now, if there are any keywords, add things to the basic plot title:                                                                                               ;
      ;                                                                                                                                                                   ;
      IF ( (apcheck NE 0) OR (radcheck NE 0) OR (bandcheck NE 0) OR (cutscheck NE 0) ) THEN BEGIN                                                                        ;;
        ; add a comma and a line after the main title, since this is extra stuff                                                                                          ;
        titlestring = titlestring + ',!C'                                                                                                                                ;;
        titlestring = string(titlestring)                                                                                                                                ;;
        ;                                                                                                                                                                 ;
        ; First, see if nApertures and/or radius are specified. Since they go hand-in-hand, treat them together:                                                          ;
        IF (radcheck NE 0) THEN radName = strcompress(string(fix(radius)), /remove) + '"'     ; if radius is provided, here's a string describing it                      ;
        IF (apcheck NE 0) THEN BEGIN            ; if nApertures is provided                                                                                               ;
          apsName = strcompress(string(fix(nApertures)), /remove) + ' apertures'    ; string describing nApertures                                                        ;
          IF (radcheck NE 0) THEN BEGIN        ; if nApertures and radius are both provided                                                                               ;
            titlestring = titlestring + string(apsName) + ' with ' + string(radName) + ' radii'                                                                          ;;
          ENDIF ELSE IF (radcheck EQ 0) THEN BEGIN            ; if nApertures is provided but not radius                                                                  ;
            titlestring = titlestring + string(apsName)                                                                                                                  ;;
          ENDIF                                                                                                                                                          ;;
        ENDIF ELSE IF ( (apcheck EQ 0) AND (radcheck NE 0) ) THEN BEGIN     ; if nApertures is not provided but radius is                                                 ;
          titlestring = titlestring + 'aperture radius of ' + string(radName)                                                                                            ;;
        ENDIF         ; if neither nApertures nor radius are provided, nothing else is added to the title string yet                                                      ;
        ;                                                                                                                                                                 ;
        ; now tack on the name of the band if that has been provided                                                                                                      ;
        IF (bandcheck NE 0) THEN BEGIN                                                                                                                                   ;;
          bandName = string(band)                                                                                                                                        ;;
          IF ( (apcheck NE 0) OR (radcheck NE 0) ) THEN BEGIN     ; if anything else is provided before band, add a comma                                                 ;
            titlestring = titlestring + ', '                                                                                                                             ;;
          ENDIF                                                                                                                                                          ;;
          titlestring = titlestring + string(bandName)                                                                                                                   ;;
        ENDIF    ; the band is provided                                                                                                                                   ;
        ;                                                                                                                                                                 ;
        ; now add info about color cuts if that has been provided                                                                                                         ;
        IF (cutscheck NE 0) THEN BEGIN                                                                                                                                   ;;
          IF (cuts EQ 'y') THEN BEGIN        ; make a string for the case that there are color cuts                                                                       ;
            cutsName = 'with color cuts'                                                                                                                                 ;;
          ENDIF ELSE IF (cuts EQ 'n') THEN BEGIN        ; make a string for data without cuts                                                                             ;
            cutsName = 'raw (no cuts)'                                                                                                                                   ;;
          ENDIF         ; sorts out what string to use to describe the cuts keyword                                                                                       ;
          ; now that the string is set up, add it to the title                                                                                                            ;
          IF ( (apcheck NE 0) OR (radcheck NE 0) OR (bandcheck NE 0) ) THEN BEGIN                                                                                        ;;
            titlestring = titlestring + ', '         ; add a comma if anything comes before the cuts text                                                                 ;
          ENDIF                                                                                                                                                          ;;
          titlestring = string(titlestring) + string(cutsName)       ; tack on the cuts text                                                                              ;
        ENDIF      ; the cuts keyword is provided                                                                                                                         ;
        ;                                                                                                                                                                 ;
      ENDIF     ; there are any keywords specified at all                                                                                                                 ;
      ;                                                                                                                                                                   ;
    ENDIF    ; this is the title for a plot                                                                                                                               ;
                                                                                                                                                                         ;;
    RETURN, string(titlestring)       ; provide the final title for the plot or file                                                                                      ;
                                                                                                                                                                         ;;
  END                                                                                                                                                                    ;;
;                                                                                                                                                                         ;
;                                                                                                                                                                         ;
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------;









FUNCTION colorsort, xcolor, ycolor, m, b, c
;------------------------------------------------------------------------------------------------------------------------------------------------;
; colorsort function                                                                                                                             ;
;   Applies a simple linear (y >= mx + b) color cut to a catalog                                                                                 ;
; INPUTS: xcolor - an array containing the blue-middle band color of every source in the data catalog                                            ; 
;         ycolor - an array containing the middle-red band color of every source in the data catalog                                             ; 
;         m - the slope of a line marking a simple color cut (unitless)                                                                          ;
;         b - the y-intercept of the simple color cut line (unitless)                                                                            ;
;         c - a "cap" value for color cuts, above which all colors are considered "good" (unitless)                                              ;
; OUTPUT: good - the index numbers of all the sources in the catalog with colors above the cut                                                   ;
; NOTES: The cut made is illustrated below. Asterisks mark "good" sources, while periods mark sources that get thrown out by the cut.            ;
;                                                                                                                                                ;
;                               |   *     *    *      *                                                                                          ;
;                               |      *          *      *                                                                                       ;
;                             c-|   *    *   __________________                                                                                  ;
;                               |  *        / .     .    .  .                                                                                    ;
;                      ycolor   |    *     /.    .        .                                                                                      ;
;                               |         /   .      ..      .                                                                                   ;
;                               |     *  /     .  .    .   .                                                                                     ;
;                               | *     /   .   .   .    .                                                                                       ;
;                               |   *  /       .   . .  .    .                                                                                   ;
;                               |_____/________________________                                                                                  ;
;                     (y = mx + b)----^       xcolor                                                                                             ;
;------------------------------------------------------------------------------------------------------------------------------------------------;
  good = WHERE ((ycolor GE (xcolor*m + b)) OR (ycolor GE c))
  RETURN, good
END 

  








; NOTE: THE FOLLOWING FUNCTION(S) AND/OR PROCEDURE(S) WAS/WERE NOT WRITTEN BY AGNAR. 
;
;
;
pro struct_replace_field, struct, tag, data, newtag=newtag
;Change the type, dimensionality, and contents of an existing structure
; field. The tag name may be changed in the process.
;
;Inputs:
; tag (string) Case insensitive tag name describing structure field to
;  modify. Leading and trailing spaces will be ignored. If the field does
;  not exist, the structure is not changed and an error is reported.
; data (any) data that will replace current contents  of 
; [newtag=] (string) new tag name for field being replaced. If not
;  specified, the original tag name will be retained.
;
;Input/Output:
; struct (structure) structure to be modified.
;
;Examples:
;
; Replace sme.wave with the arbitrary contents of wave:
;
;   IDL> struct_replace_field, sme, 'wave', wave
;
; The tag name for a field can be changed without altering the data:
;
;   IDL> struct_replace_field, clients, 'NMAE', clients.nmae, newtag='Name'
;
;History:
; 2003-Jul-20 Valenti  Initial coding

if n_params() lt 3 then begin
  print, 'syntax: struct_replace_field, struct, tag, data [,newtag=]'
  return
endif

;Check that input is a structure.
  if size(struct, /tname) ne 'STRUCT' then begin
    message, 'first argument is not a structure'
  endif

;Get list of structure tags.
  tags = tag_names(struct)
  ntags = n_elements(tags)

;Check that requested field exists in input structure.
  ctag = strupcase(strtrim(tag, 2))		;canoncial form of tag
  itag = where(tags eq ctag, nmatch)
  if nmatch eq 0 then begin
    message, 'structure does not contain ' + ctag + ' field'
    return
  endif
  itag = itag[0]				;convert to scalar

;Choose tag name for the output structure.
  if keyword_set(newtag) then otag = newtag else otag = ctag

;Copy any fields that precede target field. Then add target field.
  if itag eq 0 then begin			;target field occurs first
    new = create_struct(otag, data)
  endif else begin				;other fields before target
    new = create_struct(tags[0], struct.(0))	;initialize structure
    for i=1, itag-1 do begin			;insert leading unchange
      new = create_struct(new, tags[i], struct.(i))
    endfor
    new = create_struct(new, otag, data)	;insert new data
  endelse

;Replicate remainder of structure after desired tag.
  for i=itag+1, ntags-1 do begin
    new = create_struct(new, tags[i], struct.(i))
  endfor

;Replace input structure with new structure.
  struct = new

end





;
;                                   ||`-.___
;                                   ||    _.>
;                                   ||_.-'
;               ==========================================
;                `.:::::::.       `:::::::.       `:::::::.
;                  \:::::::.        :::::::.        :::::::\
;                   L:::::::    AGNAR LAERAD HALL    :::::::L
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

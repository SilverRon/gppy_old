#   ZEROPOINT CALCULATION
#   19.03.09    MERRY CHRISTMAS!
#   MADE BY GREGORY S.H. PAEK
#=========================================================================#
def which_obs(obs, path_and_file):
	"""
	=====================================================================
	GIVE CCD INFORMATION
	---------------------------------------------------------------------
	INPUT   :   observatory name, and path+table_name(ascii format)
				i.e.) path_and_file='/home/sonic/Research/table/obs.txt'
	OUTPUT  :   GAIN        []
				Pixel Scale ["/pix]
	---------------------------------------------------------------------
	obs         ccd         gain    RDnoise     dark    pixelscale
	SOAO        FLI         1.43    13.60       0.0     0.44540
	LOAO        E2V         2.68    4.84        0.0     0.794
	LOAO_FLI    FLI         3       15.0        0.0     1.28
	DOAO_sophia PIsophia    0.928   7.63        0.0     0.3855
	DOAO        FLI         1.27    14.0        0.005   0.464
	oldBOAO     UNKNOWN     8.9     20.0        0.0     1.7
	BOAO        UNKNOWN     5.0     2.0         0.0     0.2145
	BOAO        UNKNOWN     5.0     5.0         0.0     0.2145    
	SAO         SBIG        1.35    9.0         0.0     0.31
	30inch      UNKNOWN     1.60    5.87        0.0     1.3553
	CCA250      MLI16803    1.46    10.35       0.0     2.06
	MAIDANAK    SNU4kCAM    1.45    4.7         0.0     0.266
	MAIDANAK    UNKNOWN     5.40    9.0         0.0     0.266
	LSGT        SNUCAMII    1.15    6.0         0.2     0.92
	UKIRT       WFCAM       99.0    0.0         0.0     0.4
	=====================================================================
	"""
	import numpy as np
	from astropy.io import ascii
	obsinfo     = ascii.read(path_and_file)
	indx_obs    = np.where( obs == obsinfo['obs'] )
	gain        = np.copy(obsinfo['gain'][indx_obs])[0]
	pixscale    = np.copy(obsinfo['pixelscale'][indx_obs])[0]
	return gain, pixscale
#-------------------------------------------------------------------------#
def secom(inim, gain, pixscale, det_sigma=5, backsize=str(64), backfiltersize=str(3), dual=False, detect='detection.fits', check=False):
	"""
	SourceEXtractor for search transient
	"""
	import numpy as np
	import os
	from astropy.io import ascii
	#   FILE CHECK
	#	CONFIG FILES (USER BASE PATH)
	configfile      = '/home/sonic/Research/yourpy/config/detect.sex'
	paramfile       = '/home/sonic/Research/yourpy/config/detect.param'
	nnwfile		    = '/home/sonic/Research/yourpy/config/detect.nnw'
	convfile	    = '/home/sonic/Research/yourpy/config/detect.conv'
	try:
		comment = '\nSourceEXtractor START\n' \
				+ 'IMAGE\t\t: '+inim+'\n' \
				+ 'GAIN\t\t: '+str(gain)+'\n' \
				+ 'PIXSCALE\t\t: '+str(pixscale)+'\n' \
				+ 'DETECTION SIGMA\t: '+str(det_sigma)+'\n' \
				+ 'PARAM\t\t: '+paramfile+'\n' \
				+ 'BACKSIZE\t: '+backsize+'\n' \
				+ 'BACKFILTER\t: '+backfiltersize+'\n' \
				+ 'CONFIG\t\t: '+configfile+'\n' \
				+ 'NNW\t\t: '+nnwfile+'\n' \
				+ 'CONVOLVE\t: '+convfile
		print(comment)
	except:
		comment = 'CHECK configfile/paramfile/nnewfile/convfile or others.'
		print(comment)
	#   FILE NAME
	cat     = inim[:-5]+'.cat'
	seg     = inim[:-5]+'.seg.fits'
	bkg     = inim[:-5]+'.bkg.fits'
	sub     = inim[:-5]+'.sub.fits'
	psf     = inim[:-5]+'.psf'
	aper    = inim[:-5]+'.aper.fits'

	#   BASIC INFO.
	det_area        = 5
	det_thresh      = det_sigma/np.sqrt(det_area)
	detecminarea    = str(det_area)
	detectthresh    = str(det_thresh)
	seeing			= 10.0
	#seeing, fwhm_arcsec = psfex(inim, pixscale)
	#pixscale        = pixscalecalc(inim)
	#   OPTION
	aperture        = '%.2f'%(3./pixscale)+','+'%.2f'%(5./pixscale)+','+'%.2f'%(7./pixscale)+','+'%.2f'%(seeing)+','+'%.2f'%(1.2*seeing)+','+'%.2f'%(1.5*seeing)+','+'%.2f'%(1.7*seeing)+','+'%.2f'%(2.0*seeing)

	option1 = ' -CATALOG_NAME '+cat+' -PARAMETERS_NAME '+paramfile
	option2 = ' -DETECT_MINAREA '+detecminarea+' -DETECT_THRESH '+detectthresh \
			+' -FILTER Y '+'-FILTER_NAME '+convfile
	option3 = ' -PHOT_APERTURES '+aperture \
			+' -GAIN '+'%.1f'%(gain)+' -PIXEL_SCALE '+'%.1f'%(pixscale)
	option4 = ' -SEEING_FWHM '+'%.3f'%(seeing)+' -STARNNW_NAME '+nnwfile
	option5 = ' -BACK_SIZE '+ backsize \
			+ ' -BACK_FILTERSIZE '+ backfiltersize+' -BACKPHOTO_TYPE LOCAL'
	option6 = ' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND'
	option7 = ' -CHECKIMAGE_NAME '+seg+','+aper+','+','+bkg+','+sub
	option8 = ' -PSF_NAME '+psf
	#   COMMAND
	#   detect = detection.fits is fine image show target significantly and have good seeing and many stars i
	dualcom ='sex -c '+configfile+' '+detect+' , '+inim+' -CATALOG_NAME dual'+cat+' -PARAMETERS_NAME '+paramfile+ ' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8
	sglcom  ='sex -c '+configfile+' '+inim+' '+option1+' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8
	clearcom='sex -c '+configfile+' '+inim+' '+option1+' '+option2+' '+option3+' '+option4+' '+option5+' '+option8
	if dual == False    :
		if check == False:
			os.system(clearcom)
		else:
			os.system(sglcom)
	else                : os.system(dualcom)
	
	setbl   = ascii.read(cat)
	return setbl, cat#, seeing, fwhm_arcsec
#-------------------------------------------------------------------------#
path_and_file='/home/sonic/Research/table/obs.txt'
gain, pixscale	= which_obs(obs, path_and_file)
setbl, secat	= secom(inim, gain, pixscale)




elong_med	= np.median(setbl['ELONGATION'])
ellip_med	= np.median(setbl['ELLIPTICITY'])

seltbl		= setbl[	(setbl['FLAGS']==0) &
						(setbl['ELONGATION']<elong_med) &
						(setbl['ELLIPTICITY']<ellip_med)]

seltbl		= setbl[	(setbl['FLAGS']==0) &
						(setbl['ELONGATION']<elong_med) &
						(setbl['ELLIPTICITY']<ellip_med) &
						(setbl['MU_MAX']+27<22)]


name, ra, dec = seltbl['NUMBER'], seltbl['ALPHA_J2000'], seltbl['DELTA_J2000']

#-------------------------------------------------------------------------#
setbl = ascii.read('hdCalib-...fits')
setbl.meta = dict(obs=obs, naxis=fits.getheader(inim)['NAXIS1'])
setbl['FROM_CENT_IMAGE'] = np.sqrt(
	(setbl['X_IMAGE']-setbl.meta['naxis']/2)**2+(setbl['Y_IMAGE']-setbl.meta['naxis']/2)**2)


seltbl = setbl[	(setbl['ELONGATION'] <= mean_elong) &
                (setbl['ELLIPTICITY'] <= mean_ellip) &
                (setbl['CLASS_STAR'] <= mean_cs) &
                (setbl['FLAGS'] == 0) &
                (setbl['BACKGROUND'] > 0) &
                (setbl['FWHM_IMAGE'] > mean_fwhm) &
                (setbl['FROM_CENT_IMAGE'] < 0.9*setbl.meta['naxis']/2)]

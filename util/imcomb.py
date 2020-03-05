def epochimcomb(imlist, outim='imcomb.fits', path_save='.'):
	'''
	epochimcomb(imlist, outim='imcomb.fits', path_save='.')

	imlist = glob.glob('Calib*20181229*.fits')
	epochimcomb(imlist)
	'''
	#------------------------------------------------------------
	import os, glob
	import numpy as np
	from astropy.nddata import CCDData, fits_ccddata_reader, fits_ccddata_writer
	from matplotlib import pyplot as plt  
	from ccdproc import Combiner
	from astropy.time import Time
	from astropy.io import fits
	from imsng import phot
	#------------------------------------------------------------
	#	EXTRACT INFO. FROM THE FIRST IMAGE
	#------------------------------------------------------------
	data0 = fits_ccddata_reader(imlist[0])
	meta0 = data0.meta
	wcs0 = data0.wcs
	part = imlist[0].split('-')
	#------------------------------------------------------------
	#	IMAGE COMBINE
	#------------------------------------------------------------
	comlist = []
	dateobslist = []
	explist = []
	print('{} IMAGE COMBINE START\n'.format(len(imlist)))
	for inim in imlist:
		print(inim)
		hdr = fits.getheader(inim)
		dateobslist.append(Time(hdr['DATE-OBS'], format='isot').jd)
		explist.append(hdr['EXPTIME'])
		comlist.append(fits_ccddata_reader(inim))

	dateobs = Time(np.mean(dateobslist), format='jd')
	totexp = np.sum(explist)
	try:
		comin = '{}-{}-{}-{}-{}-{}-{}-com.fits'.format(part[0], part[1], part[2], dateobs.isot[0:10].replace('-', ''), dateobs.isot[11:19].replace(':', ''), part[5], int(totexp))
	except:
		print('IMAGE NAME FORMAT IS NOT Calib-... .fits.)
		comin = outim
	c = Combiner(comlist)
	cdata = c.median_combine()
	cdata.meta = meta0
	cdata.wcs = wcs0
	print('OUTPUT IMAGE :\t{}\n'.format(comin))
	fits_ccddata_writer(cdata, path_save+'/'+comin)
	#------------------------------------------------------------
	phot.puthdr(comim, 'TOTEXP', totexp, hdrcomment='Total exposure time in seconds')
	phot.puthdr(comim, 'JD', dateobs.jd, hdrcomment='Center Julian Date at start of exposure')
	phot.puthdr(comim, 'MJD', dateobs.mjd, hdrcomment='Center Modified Julian Date at start of exposure')
	phot.puthdr(comim, 'DATE-OBS', dateobs.isot, hdrcomment='YYYY-MM-DDThh:mm:ss observation start, UT')
	phot.puthdr(comim, 'NCOMBINE', len(imlist), hdrcomment='THE NUMBER OF COMBINED IMAGES')
	for i, inim in enumerate(imlist):
		phot.puthdr(comin, 'COMBINE{}'.format(i+1), inim, hdrcomment='{} COMBINED IMAGE'.format(i+1))
	print('DONE')
	return comin
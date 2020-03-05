def swarpcom(imkey, listname='obj.list', path_save='.', path_obs = '/home/sonic/Research/table')
	import os, glob
	import numpy as np
	from imsng import tool, phot
	'''
	imkey = 'Calib*-20180129-*0.fits'
	path_save = '.'
	path_obs = '/home/sonic/Research/table'
	listname = 'obj.list'
	'''

	imlist = glob.glob(imkey); imlist.sort()
	f = open(listname, 'w')
	for inim in imlist:
		f.write(inim+'\n')
		# print(inim)
	f.close()

	comim, hdr, dateobs, totexp = tool.combname(imlist)
	part = comim.split('-')

	gain, pixscale, fov = tool.getccdinfo(part[1], path_obs)

	conf = 'default.swarp'
	os.system('swarp -d > {}/{}'.format(path_save, conf))

	com = 'swarp @{} -c {} -IMAGEOUT_NAME {} -COMBINE_TYPE MEDIAN -RESAMPLE N -PIXEL_SCALE {} -GAIN_DEFAULT {} -SUBTRACT_BACK Y'.format(listname, conf, comim, pixscale, gain)

	print(com)
	os.system(com)

	phot.puthdr(comim, 'OBJECT', hdr['OBJECT'], hdrcomment='OBJECT')
	phot.puthdr(comim, 'TOTEXP', totexp, hdrcomment='Total exposure time in seconds')
	phot.puthdr(comim, 'JD', dateobs.jd, hdrcomment='Center Julian Date at start of exposure')
	phot.puthdr(comim, 'MJD', dateobs.mjd, hdrcomment='Center Modified Julian Date at start of exposure')
	phot.puthdr(comim, 'DATE-OBS', dateobs.isot, hdrcomment='YYYY-MM-DDThh:mm:ss observation start, UT')
	phot.puthdr(comim, 'NCOMBINE', len(imlist), hdrcomment='THE NUMBER OF COMBINED IMAGES')
	for i, inim in enumerate(imlist):
		phot.puthdr(comim, 'COMBINE{}'.format(i+1), inim, hdrcomment='{} COMBINED IMAGE'.format(i+1))

	return comim



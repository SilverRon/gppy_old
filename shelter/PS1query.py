#   PS1 QUERY & GREGISTERING & SUBTRACTION
#   2019.11.12 MADE		BY Gregory S.H. Paek
#============================================================#
#	MODULE
#------------------------------------------------------------#
import os, sys, glob
from imsng import query
from astropy.io import fits
import numpy as np
#------------------------------------------------------------#
#	FUNCTION
#------------------------------------------------------------#
def RepeatQuery(param_down):
	try:
		query.downimage_routine(**param_down)
	except:
		try:
			query.downimage_routine(**param_down)
		except:
			try:
				query.downimage_routine(**param_down)		
			except:
				try:
					query.downimage_routine(**param_down)		
				except:
					query.downimage_routine(**param_down)
#------------------------------------------------------------#
# path_base = '/data1/S190425z/SAO'
path_base = '/data1/S190425z/KMTNet/snap/190425'
imlist = glob.glob(path_base+'/aSnap-*.fits')

failist = []
for inim in imlist:
	hdr = fits.getheader(inim)
	name = inim.split('-')[2]
	ra, dec = hdr['CRVAL1'], hdr['CRVAL2']
	param_down	= dict(  name = name,
						ra = ra, dec = dec,
						filters = 'r',
						size = 6000,
						output_size = None,
						format = 'fits',
						save_dir=path_base)
	try:
		RepeatQuery(param_down)
	except:
		failist.append(inim)

#------------------------------------------------------------#
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.nddata.utils import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
#------------------------------------------------------------#
path_cand = '/data1/S190425z/info/Update/S190425z_Update-all_candi.txt'
cantbl = ascii.read(path_cand)
cutsize = 300
#------------------------------------------------------------#
# inim = 'Calib-LOAO-PGC1294635-20190426-101613-R-120_gregister-com.fits'
# ra = np.asscalar( cantbl['ra'][cantbl['name']==obj] )
# dec = np.asscalar( cantbl['dec'][cantbl['name']==obj] )
# imlist = glob.glob('a*.fits')
imlist = glob.glob('*Calib*.fits')
# for inim in imlist[191:]:
for inim in imlist:
	obs = inim.split('-')[1]
	# obs = 'KMTNET'
	hdr = fits.getheader(inim)
	dateobs = hdr['DATE-OBS']
	ymd, hms = dateobs.split('T')[0], dateobs.split('T')[1]
	ymd = ymd.replace('-', '')
	hms = hms.replace(':', '')
	hms = hms[:6]
	w = WCS(inim)

	x1, y1 = 1, 1
	x2, y2 = w.array_shape[0], w.array_shape[1]

	for i, obj in enumerate( cantbl['name'] ):
		ra = np.asscalar( cantbl['ra'][i] )
		dec = np.asscalar( cantbl['dec'][i] )
		radec = np.array([[ra, dec]])

		# if not ((np.isnan(xy[0][0])) | (np.isnan(xy[0][1]))):
		# if ( (np.isnan(xy[0][0])==False) & (np.isnan(xy[0][1])==False) ):
		xy = w.wcs_world2pix(radec, 0)
		try:
			if (not np.isnan(xy[0][0])) & ( (xy[0][0]!=np.inf) & (xy[0][1]!=np.inf) ):
				x, y = np.int(xy[0][0]), np.int(xy[0][1])
				if '[' in obj:
					obj = obj.replace('[', 'p')
					obj = obj.replace(']', 'q')
				outim = 'Snap-{}-{}-{}-{}.fits'.format(obs, obj, ymd, hms)
				if (x1 < x < x2) & (y1 < y < y2):
					print((x, y))

					xmin, xmax = x-cutsize, x+cutsize
					ymin, ymax = y-cutsize, y+cutsize
					if xmin < x1:
						xmin = x1
					if xmax > x2:
						xmax = x2
					if ymin < y1:
						ymin = y1
					if ymax > y2:
						ymax = y2
					imcopycom = 'imcopy {}[{}:{},{}:{}] {}'.format(inim, xmin, xmax, ymin, ymax, outim)
					print(imcopycom)
					os.system(imcopycom)
				# else:
					# print('out of range')
					if (('q' in outim) & ('p' in outim)):
						obj = obj.replace('p', '[')
						obj = obj.replace('q', ']')
						newim = 'Snap-{}-{}-{}-{}.fits'.format(obs, obj, ymd, hms)
						os.system('mv {} {}'.format(outim, newim))
		except:
			failist.append(inim)

for inim in glob.glob('Snap*.fits'):
	#	IMAGES
	w = WCS(inim)
	x1, y1 = 1, 1
	x2, y2 = w.array_shape[0], w.array_shape[1]
	xcent, ycent = np.int(x2/2.), np.int(y2/2.)
	xycent = np.array([[xcent, ycent]])
	imcoord = w.wcs_pix2world(xycent, 0)

	racent, decent = imcoord[0][0], imcoord[0][1]
	cim = SkyCoord(racent, decent, unit='deg')

	#	TABLE TARGET
	obj = inim[12:-21]
	# if ('q' in obj) and ('p' in obj):
		# obj = obj.re
	indx = np.where(cantbl['name'] == obj)
	ra = np.asscalar( cantbl['ra'][indx] )
	dec = np.asscalar( cantbl['dec'][indx] )
	ctarg = SkyCoord(ra, dec, unit='deg')
	sep = cim.separation(ctarg)
	if sep.to_value('arcsec') > 150:
		# print(inim)
		os.system('rm {}'.format(inim))
	# if cim.separation(ctarg)*3600 < 150*u.deg:
		# print(inim)





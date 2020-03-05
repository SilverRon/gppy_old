#============================================================
#	Pan-STARRS DR1 data query + 
#	
#	UPDATE : 20.02.20
#============================================================
from __future__ import print_function
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io import ascii
from astropy import units as u
import numpy as np
import sys
import os, glob
from imsng import query
from astropy.wcs import WCS
#============================================================
#	FUNCTION
#============================================================
def ps1tile_routine(inim):
	# inim = 'Calib-LOAO-UGC10320-20190425-103020-R-60_gregister-com.fits'
	obj = inim.split('-')[2]
	data, hdr = fits.getdata(inim, header=True)
	w = WCS(inim)

	shp = data.shape

	pt0x, pt0y = shp[0]*(1/2), shp[1]*(1/2)

	pt1x, pt1y = shp[0]*(1/4), shp[1]*(1/4)
	pt2x, pt2y = shp[0]*(3/4), shp[1]*(1/4)
	pt3x, pt3y = shp[0]*(3/4), shp[1]*(3/4)
	pt4x, pt4y = shp[0]*(1/4), shp[1]*(3/4)


	pt0 = w.pixel_to_world(pt0x, pt0y)
	pt1 = w.pixel_to_world(pt1x, pt1y)
	pt2 = w.pixel_to_world(pt2x, pt2y)
	pt3 = w.pixel_to_world(pt3x, pt3y)
	pt4 = w.pixel_to_world(pt4x, pt4y)


	for i, pt in enumerate([pt0, pt1, pt2, pt3, pt4]):
		outobj = '{}_{}'.format(obj, i)
		outim = 'Calib-PS1-{}-0-0-r-0.fits'.format(outobj)
		param_down	= dict( outim = outim,
							name = outobj,
							ra = pt.ra.value,
							dec = pt.dec.value,
							size = 1500,
							output_size = None,
							filters = 'i',
							format = 'fits',
							save_dir='.')
		iternumb = 0
		# while ((len(glob.glob(outim)) >= 1) | (iternumb < 5)):
		while iternumb < 5:
			iternumb += 1
			try:
				query.downimage_routine(**param_down)
				break
			except:
				print('FAIL TO DOWNLOAD FOR {} (TRIAL {})'.format(outim, iternumb))
				pass

	# querylist = glob.glob('Calib-PS1-{}_*-*.fits'.format(obj)); querylist.sort()

	comin = 'Calib-PS1-{}-0-0-r-0.fits'.format(obj)
	swarpcom = 'swarp Calib-PS1-{}_*-*.fits -c {} -IMAGEOUT_NAME {}'.format(obj, path_swarp, comin)
	os.system(swarpcom)
#============================================================
path_swarp = '/home/sonic/Research/yourpy/gppy/config/default.swarp'

imlist = glob.glob('aCalib-*.fits'); imlist.sort()
for i, inim in enumerate(imlist):
	print('[{}/{}]'.format(i+1, len(imlist)), inim)
	ps1tile_routine(inim)
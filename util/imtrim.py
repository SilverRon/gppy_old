#============================================================
#	IMAGE TRIM WITH ASTROPY
#	2020.02.13	CREATED BY GREGORY S.H. PAEK
#============================================================
import numpy as np
import os, sys, glob
from astropy.io import fits
# import astropy.coordinates as coord
import astropy.units as u
# import time
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
#============================================================
#	FUNCTION
#============================================================
def trim(inim, position, size, outim='trim.fits'):
	# Load the image and the WCS
	hdu = fits.open(inim)[0]
	wcs = WCS(hdu.header)
	# Make the cutout, including the WCS
	cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
	# Put the cutout image in the FITS HDU
	hdu.data = cutout.data
	# Update the FITS header with the cutout WCS
	hdu.header.update(cutout.wcs.to_header())
	# Write the cutout to a new FITS file
	hdu.writeto(outim, overwrite=True)

#============================================================
#	INPUT
#============================================================
# inim = 'ref.fits'
# outim = 'tr'+inim
# tra, tdec = 44.5438520, -8.9577875				#	GRB 190829A
# tra, tdec = 44.863251, 31.385878					#	G0037111
# tra, tdec = 45.3575160, 31.5859348				#	G0250092
# tra, tdec = 45.3660757, 31.5317386				#	G0250092 FOR DOAO
# tra, tdec = 45.4030528, 31.8194872				#	G0049337
# tra, tdec = 45.4208136, 31.8949936				#	G0232794
tra, tdec = 42.1846250, 12.1372444				#	AT2020yxz

position = SkyCoord(tra, tdec, frame='icrs', unit='deg')
length = 5	#	[']
#============================================================
size = u.Quantity((length, length), u.arcmin)
# trim(inim, position, size, outim)
# plt.imshow(cutout.data, origin='lower')

os.system('ls *.fits')
try:
	#	PYTHON 2.X
	imlist = glob.glob(raw_input('IMAGES TO PROCESS\t:')); imlist.sort()
except:
	#	PYTHON 3.X
	imlist = glob.glob(input('IMAGES TO PROCESS\t:')); imlist.sort()
print(imlist)
for inim in imlist:
	trim(inim, position, size, 'tr'+inim)

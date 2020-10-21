#============================================================
#	IMAGE TRIM - WCSREMAP - SUBTRACTION
#	2020.02.16	CREATED BY GREGORY S.H. PAEK
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
from imsng import tool
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
#------------------------------------------------------------
def hotpants(inim, refim):
	outfile = 'hd'+inim
	convfile = 'hc'+inim
	com = 'hotpants -c t -n i -iu 60000 -tu 6000000000 -tl -100000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
	os.system(com)


#============================================================
# tra, tdec = 45.3660757, 31.5317386				#	G0250092 FOR DOAO
# tra, tdec = 14.7255135, 52.6095117				#	G0010044
# tra, tdec = 18.9285439, 55.9965139				#	G0017086
# tra, tdec = 18.9509304, 56.0405869				#	G0017086 FOR LOAO
# tra, tdec = 19.8631833, 58.5216117				#	G0265297 FOR LOAO
# tra, tdec, length = 44.863251, 31.385878, 5		#	G0037111
# tra, tdec = 45.3575160, 31.5859348				#	G0250092
# tra, tdec = 45.4007437, 31.5461861				#	G0250092 FOR SOAO
# tra, tdec = 45.3660757, 31.5317386				#	G0250092 FOR DOAO
# tra, tdec = 45.4030528, 31.8194872				#	G0049337
# tra, tdec = 45.4208136, 31.8949936				#	G0232794
# tra, tdec = 29.2019712, 36.0878171				#	G0079760 FOR LOAO
# tra, tdec = 19.6274302, 51.9457305				#	G0003842 
# tra, tdec = 15.4419912, 55.9492553				#	G0016240 
# tra, tdec = 15.4787597, 55.8646368				#	G0016240 FOR LOAO
# tra, tdec = 29.2014766, 35.9937460				#	G0079760
# tra, tdec = 22.5399478, 57.8445955				#	G0143225
# tra, tdec = 3.2529610, 30.9158691					#	G0002715
# tra, tdec, length = 2.4604205, 33.3074615, 30		#	G0002650 FOR LOAO
# tra, tdec, length = 15.7359666, 51.3422759, 10	#	G0008283
# tra, tdec, length = 14.6530454, 51.6189949, 10	#	G0189676
# tra, tdec, length = 75.9126992, -23.7941181, 10	#	UNKNOWN FOR GW190425
# tra, tdec, length = 170.083333, +12.990278, 20	#	2020cwh
# tra, tdec, length = 170.1247667, 26.95351111, 10	#	2020dbf
# tra, tdec, length = 44.5440417, -8.9579444, 20	#	GRB 190829A
# tra, tdec, length = 185.7288542, 15.8236250, 5	#	SN 2020oi
# tra, tdec, length = 111.8296542, 80.2328528, 10		#	2020ddy
# tra, tdec, length = 111.8296542, 80.2328528, 10		#	2020ddy
# tra, tdec, length = 185.4547125, 4.3485889, 20		#	2020hxk
# tra, tdec, length = 185.4603292, 4.481705551, 5		#	2020hxk
# tra, tdec, length = 185.4603292, 4.4817056, 5  # 2020hxk
# tra, tdec, length = 185.7338750, 15.8260000, 5
# tra, tdec, length = 258.3414542, -9.9644667, 1  #	ZTF19aarykkb
# tra, tdec, length = 262.7914875, -8.4507222, 1  #	ZTF19aarzaod
# tra, tdec, length = 185.4603292, 4.4817056, 5	#	2020jfo
# tra, tdec, length = 254.034297, -8.1985399, 1.5 # AT2019dnv
# tra, tdec, length = 86.3352843, -26.8477109, 10 # AT2019flz
# tra, tdec, length = 29.799542, 18.981944, 10 # AT2020uex
# tra, tdec, length = 21.0286875, 12.92148055, 10  # AT2020uxz
# tra, tdec, length = 308.7195919, 60.1547330, 5  # AT2020uxz
tra, tdec, length = 271.9568721, 17.6887043, 6.2  # Test for CBNUO (NGC6555)
# tra, tdec, length = 308.7282567, 60.1495581, 13	# Test for CBNUO (NGC6946) - LOAO



#============================================================
position = SkyCoord(tra, tdec, frame='icrs', unit='deg')
#	arcmin
# length = 6.5
# length = 6
size = u.Quantity((length*2, length*2), u.arcmin)
#============================================================
# os.system('rm wr*.fits tr*.fits Ref*.fits h*.fits')
os.system('rm wr*.fits Ref*.fits h*.fits')
os.system('ls *.fits')
# inim = raw_input('INPUT\t: ')
imlist = glob.glob(raw_input('INPUT\t: ')); imlist.sort()
refim = raw_input('REF\t: ')
if refim == '':
	refim = 'ref.fits'

for inim in imlist:
	trim(inim, position, size, 'tr'+inim)
	trim(refim, position, size, 'tr'+refim)

	trinim = 'tr'+inim
	trrefim = 'tr'+refim

	wrim = tool.wcsremap(trrefim, trinim)		#	MATCH REF. TO TARGET IMAGE
	hotpants(trinim, wrim)

	ds9com = 'ds9 {} {} {} -lock frame wcs -tile column &'.format(trinim, wrim, 'hd'+trinim)
	# os.system(ds9com)
	cleancom = 'rm hc*.fits *seg.fits'
	os.system(cleancom)

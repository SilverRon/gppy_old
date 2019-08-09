#	CHANGE HEADER (FUNCTION)
#	MADE BY Gregory S.H. Paek   (2019.07.25)
#============================================================
#	MODULES
#============================================================
import numpy as np
from astropy.io import fits
import os,glob
#------------------------------------------------------------
def changehdr(inim, where, what):
	data, hdr = fits.getdata(inim, header = True)
	hdr[where] = what
	fits.writeto(inim, data, hdr, overwrite=True)
#------------------------------------------------------------
# SQUEAN DATE-OBS CORRECTION
sample_applications='''
imlist = glob.glob('*.fits')
for inim in imlist:
	data, hdr = fits.getdata(inim, header=True)
	if 'T' not in hdr['DATE-OBS']:
		dateobs = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
		changehdr(inim, 'DATE-OBS', dateobs)
		print(inim, dateobs)'''
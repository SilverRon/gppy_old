import astropy.io.fits as fits
import astropy.wcs as wcs
import numpy as np
import os, glob
from astropy.wcs import WCS
os.system('ls *.fits -h')
imlist = glob.glob(input('IMAGES TO PROCESS (Cal*.fits) : '))

def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)	
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'


import os,sys, glob
from pyraf import iraf
#------------------------------------------------------------
def imcopy(inim, ranges):
	outname	= 'tr'+inim
	chinim	= inim+ranges
	iraf.imcopy(chinim,output=outname)
#------------------------------------------------------------
os.system('ls *.fits *.fit')
imlist	= glob.glob(raw_input('IMAGES TO PROCESS\t: '))
xranges	= raw_input('XRANGES (x1:x2)\t: ')
yranges	= raw_input('YRANGES (y1:y2)\t: ')

ranges	= '['+xranges+','+yranges+']'

for inim in imlist:
	imcopy(inim, ranges)


#	UKIRT	w.all_pix2world(2068, 2068, 0)
#	KMTNet	w.all_pix2world(760, 760, 0)
#	SAO		w.all_pix2world(2187, 2140, 0)

for inim in imlist:
	w = WCS(inim)
	print(w)
	ra, dec = w.all_pix2world(760, 760, 0)
	puthdr(inim, 'RA', float(ra))
	puthdr(inim, 'DEC', float(dec))

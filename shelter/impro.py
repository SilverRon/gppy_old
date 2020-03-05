#   IMAGE PROCESS MODUEL
#   2019.11.13 MADE		BY Gregory S.H. Paek
#============================================================
#	MODUEL
#------------------------------------------------------------
import os, glob
import numpy as np
import alipy
from astropy.io import ascii
#------------------------------------------------------------
#	FUNCTION
#------------------------------------------------------------
def wcsremap(inim, tmpim):
	'''
	i.e.)
	tmpim	: PS1 image				(pixel scale ~ 0.25"/pix)
	inim	: LOAO image			(pixel scale ~ 0.80"/pix)
	outim	: Transformed PS1 image	(pixel scale ~ 0.80"/pix)
	'''
	outim = 're'+tmpim
	com = 'wcsremap -template {} -source {} -outIm {}'.format(inim, tmpim, outim)
	os.system(com)
	return outim
#------------------------------------------------------------
def hotpants(inim, refim):
	outim = 'hd'+os.path.basename(inim)
	convim = 'hc'+os.path.basename(inim)
	com = 'hotpants -c t -n i -iu 600000 -tu 6000000 -tl -10000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outim+' -oci '+convim
	os.system(com)
	return outim
#------------------------------------------------------------
#============================================================
# os.system('ls *.fits')
'''
tmpim = input('TEMPLATE IMAGE (i.e., ref.fits) :')
inimkey = input('INPUT IMAGE (i.e., Calib*com.fits) :')
imlist = glob.glob(inimkey)
for inim in imlist:
	reim = wcsremap(inim, tmpim)
	subim = hotpants(inim, reim)
'''
# head = 'aSnap'
head = 'Calib'
reflist = glob.glob('Ref-*.fits')
for refim in reflist:
	obj = refim.split('-')[2]
	if len(glob.glob('{}-*{}*.fits'.format(head, obj))) >= 1:
		tmpim = refim
		imlist = glob.glob('{}-*{}*.fits'.format(head, obj))
		for inim in imlist:
			reim = wcsremap(inim, tmpim)
			subim = hotpants(inim, reim)
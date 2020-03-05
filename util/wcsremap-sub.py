#============================================================
#	IMAGE WCSREMAP - SUBTRACTION
#	2020.02.24	CREATED BY GREGORY S.H. PAEK
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
def hotpants(inim, refim):
	outfile = 'hd'+inim
	convfile = 'hc'+inim
	com = 'hotpants -c t -n i -iu 60000 -tu 6000000000 -tl -100000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
	os.system(com)
#============================================================
os.system('rm tr*.fits Ref*.fits h*.fits')
os.system('ls *.fits')
# inim = raw_input('INPUT\t: ')
imlist = glob.glob(input('INPUT\t: ')); imlist.sort()
# refim = input('REF\t: ')

failist = []

for inim in imlist:
	try:
		refim = glob.glob('Calib-PS1-{}-*.fits'.format(inim.split('-')[2]))[0]
		wrim = tool.wcsremap(refim, inim)		#	MATCH REF. TO TARGET IMAGE
		hotpants(inim, wrim)
	except:
		failist.append(inim)

	# ds9com = 'ds9 {} {} {} -lock frame wcs -tile column &'.format(trinim, 'Ref'+trinim, 'hd'+trinim)
	# os.system(ds9com)

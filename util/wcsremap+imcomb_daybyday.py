import os, glob
from astropy.io import fits
import numpy as np
# import astropy.units as u
# from astropy.time import Time
from imsng import tool
#============================================================
#	FUNCTION
#============================================================

#============================================================
# inim = 'Calib-SAO_STX16803-NGC5350-20190502-112635-V-120.fits'
# refim = 'Ref-SDSS-NGC5350-r-SAO.fits'
# os.system('ls*.fits')
#------------------------------------------------------------
datelist = []
for inim in glob.glob('C*0.fits'):
	part = inim.split('-')
	datelist.append(part[3])
datelist = list(set(datelist))
#------------------------------------------------------------
failist = []
comlists = []
for date in datelist[0:2]:
	comlist = []
	subimlist = glob.glob('C*-{}-*0.fits'.format(date)); subimlist.sort()
	for inim in subimlist:
		comlist.append(tool.wcsremap(subimlist[0], inim))
	comlists.append(tool.epochimcomb(comlist))

'FAIL CODE XXXXX'
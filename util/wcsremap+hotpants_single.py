import os, glob
from imsng import tool

# inim = 'Calib-SAO_STX16803-NGC5350-20190502-112635-V-120.fits'
# refim = 'Ref-SDSS-NGC5350-r-SAO.fits'
os.system('ls*.fits')
inim = input('INPUT IMAGE :\t')
refim = input('REF. IMAGE :\t')

wrim = tool.wcsremap(refim, inim)		#	MATCH REF. TO TARGET IMAGE
subim = tool.hotpants(inim, wrim)		#	SUBTRACTION

ds9com = 'ds9 {} {} {} {} -lock frame wcs &'.format(inim, refim, wrim, subim)
print(ds9com)
# os.system(ds9com)
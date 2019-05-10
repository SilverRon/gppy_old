import glob, os
'''
imlist	= glob.glob('Calib-*2MASS+*-*-*.fits')
for inim in imlist:
	head	= inim[:31]
	tail	= inim[32:]
	new		= head+'+'+tail
	com		= 'cp '+inim+' '+new
	print(com)
	os.system(com)
'''
objlist		= []
for inim in glob.glob('Calib-*0.fits'):
	part	= inim.split('-')
	obj		= part[2]
	objlist.append(obj)
#	part[1]	= 'SQUEAN'
#	new		= '-'.join(part)
#	com		= 'cp '+inim+' '+new
#	print(com)
#	os.system(com)
objlist		= list(set(objlist)); objlist.sort()

for obj in objlist:
	sublist	= glob.glob('Calib-*'+obj+'*0.fits')
	if len(sublist) == 1:
		pass
	else:
		print(obj)


'''
#	20190424
2MASS+17030866+0237020
IC4587
NGC6097
'''

#	ASTROMETRY CHECK
import glob, os
imlist	= glob.glob('aCalib-*.fits')
ds9com	= 'ds9 '
catcom	= '-catalog 2mass '
for inim in imlist:
	ds9com	+= inim+' '+catcom
ds9com	+= '&'
print(ds9com)
os.system(ds9com)


#	IMAGES SUMMARY
from astropy.table import Table
import numpy as np
intbl	= Table()
imlist	= glob.glob('Calib-*0.fits')
alist	= glob.glob('aCalib-*0.fits')
intbl['name']	= np.array(imlist)



#	DRAW ARROW IN LIGHTCURVE
P.arrow(nd['jd'][0]-t0, nd['mag'][0], 0.0, +0.2, fc='k', ec='k', head_width=0.05, head_length=0.1)       

intblR	= intbl[intbl['band']=='R']
cutR = intblR[(intblR['mag']>20.0) & (intblR['seeing']<2.0)]


intblB	= intbl[intbl['band']=='B']
cutB = intblB[(intblB['mag']>19.0) & (intblB['seeing']<3.0)]


intblV	= intbl[intbl['band']=='V']
cutV = intblV[(intblV['mag']>19.5) & (intblV['seeing']<3.0)]


intblI	= intbl[intbl['band']=='I']
cutI	= intblI[(intblI['mag']>19.5) & (intblI['seeing']<3.0)]
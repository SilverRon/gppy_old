import os, glob

'''
i.e.)
Calib-LSGT-16573248-2246391-20190425-182951-r-180.fits
'''

imlist = glob.glob('Calib-*.fits')

#	16573248-2246391 -> 16573248+2246391
for inim in imlist:
	new		= inim[:25]+'+'+inim[26:]
	com		= 'mv '+inim+' '+new
	print(com)
	os.system(com)

#	16573248+224639	-> 2MASS+16573248+224639
for inim in imlist:
	part	= inim.split('-')
	part[2]	= '2MASS+'+part[2]
	new		= '-'.join(part)
	com		= 'mv '+inim+' '+new
	print(com)
	os.system(com)

'''
#	59064 -> PGC59064
imlist	= glob.glob('Calib-*.fits')
for inim in imlist:
	part	= inim.split('-')
	part[2]	= 'PGC'+part[2]
	new		= '-'.join(part)
	com		= 'mv '+inim+' '+new
	print(com)
	os.system(com)
'''

#	KMTNet
#	Calib-KMTNET-17161661-0811018-190425-141933-R-120.fits
#	Calib-KMTNET-2MASS+17161661+0811018-190425-141933-R-120.fits
for inim in imlist:
	part	= inim.split('-')
	part[2]	= '2MASS+'+part[2]
	new		= '-'.join(part)
	new		= new[:27]+'+'+new[28:]
	com		= 'mv '+inim+' '+new
	print(com)
	os.system(com)
#	KMTNet
#	Calib-KMTNET-ESO488-004-190427-172617-R-120.fits
#	Calib-KMTNET-ESO488+004-190427-172617-R-120.fits
for inim in imlist:
	part	= inim.split('-')
	new		= '-'.join(part)
	new		= new[:19]+'+'+new[20:]
	com		= 'mv '+inim+' '+new
	print(com)
	os.system(com)
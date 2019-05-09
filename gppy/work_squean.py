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
#	IRAF IMCOPY
#	2019.05.12 MODIFIED BY	Gregory S.H. Paek
#============================================================
import os, sys, glob
import time
from pyraf import iraf
#------------------------------------------------------------
def imcopy(inim, ranges, outim):
	chinim	= inim+ranges
	iraf.imcopy(chinim, output=outim)
#------------------------------------------------------------
time.sleep(5000)
print('done')
imlist	= glob.glob('a*.*.fits')
for inim in imlist:
	for i in [1, 2, 3, 4, 5, 6, 7, 8]:
		xranges	= '{}:{}'.format(1+1152*(i-1), 1152*i)
		yranges	= '1:9232'
		ranges	= '['+xranges+','+yranges+']'
		outim = '{}-{}.fits'.format(inim[:-5], i)
		imcopy(inim, ranges, outim)

#	IRAF IMCOPY
#	2019.05.12 MODIFIED BY	Gregory S.H. Paek
#============================================================
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

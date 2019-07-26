#============================================================
#	astrometry.net script making program for ONE target
#	2018/08/29 G.Paek revision for skipping astrocrappy
#	JUST RUN AND WRITE INPUT & REF IMAGE
#	2019.05.08	GREGORY S.H. PAEK
#============================================================
import numpy as np
import os, sys, glob
from astropy.io import fits
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
#============================================================
def astrometry(imlist, pixelscale):
	starttime	= time.time()
	scale_low       = str( float(pixelscale) - (float(pixelscale)*0.1) )
	scale_high      = str( float(pixelscale) + (float(pixelscale)*0.1) )
	n	= 1
	for inim in imlist:
		com		= 'solve-field '+inim+' --scale-unit arcsecperpix --scale-low '+scale_low+' --scale-high '+scale_high+ ' --no-plots --new-fits a'+inim+' --overwrite --temp-dir ./ '+ '--cpulimit 90 --use-sextractor\n'
		print('['+str(n)+'/'+str(len(imlist))+']'); n += 1
		print(com)
		os.system(com)
	os.system('rm tmp*')
	os.system('rm *.axy *.corr *.xyls *.match *.rdls *.solved *.wcs *axy *.corr')
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#-------------------------------------------------------------
os.system('ls *.fits *.fit')
imlist			= glob.glob(raw_input('Images to process (Calib*.fits)\t: ')); imlist.sort()
if imlist	== ''   :
	imlist	= 'Calib*.fits'
pixelscale		= raw_input('Pixel scale of Images\t: ')
if pixelscale	== ''	:
	pixelscale	= '0'
#------------------------------------------------------------
#	MULTI PROCESSING
#------------------------------------------------------------
if __name__ == '__main__':
	jobs	= []
	p			= mp.Process(target=astrometry, args=(imlist, pixelscale))
	jobs.append(p)
	p.start()
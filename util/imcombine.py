#============================================================
#	PYTHON SCRIPT FOR MEDIAN COMBINE WITH IRAF
#	Usage : 
#	python imcombine.py
#	OR JUST RUN AND WRITE INPUT & REF IMAGE
#	2019.05.08	GREGORY S.H. PAEK
#============================================================
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from pyraf import iraf
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
#============================================================
def imcombine(imlist):
	refim		= imlist[0]
	starttime	= time.time()
	print('LIST '+str(len(imlist)))
	f	= open('imcombine.list', 'w')
	for inim in imlist:
		print(inim)
		f.write(inim+'\n')
	f.close()
	try:
		jdmean, utmean	= centertime(imlist)
		utpart			= utmean.split('T')[1].split(':')
		utpart[2]		= str(int(float(utpart[2])))
		utformat		= ''.join(utpart)
		compart			= refim.split('-')
		compart[0], compart[4], compart[6]	= 'Calib', utformat, compart[6][:-5]+'-com.fits'
		comin			= '-'.join(compart)
	except:
		jdmean, utmean	= centertime2(imlist)
		utpart			= utmean.split('T')[1].split(':')
		utpart[2]		= str(int(float(utpart[2])))
		utformat		= ''.join(utpart)
		compart			= refim.split('-')
		compart[0], compart[4], compart[6]	= 'Calib', utformat, compart[6][:-5]+'-com.fits'
		comin			= '-'.join(compart)
	iraf.imcombine('@imcombine.list',output=comin,combine="median",project="no",reject="none",scale="none",zero="mode")
	comdata, comhdr	= fits.getdata(comin, header=True)
	comhdr['DATE-OBS']	= utmean
	fits.writeto(comin, comdata, comhdr, clobber=True)
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#-------------------------------------------------------------
def centertime(imlist) :	
	'''
	CALCULATE MEAN DATE-OBS
	RETURN JD AND UTC
	'''
	obsdt	= []
	for inim in imlist:
		header	=fits.getheader(inim)
		hobsdate=header['DATE-OBS']
		obsdt.append(hobsdate)
	tt			= Time(obsdt, format='isot', scale='utc')
	ttjd		= tt.jd
	ttjdmean	= np.mean(ttjd)
	ttjdmeanutc	= Time(ttjdmean,format='jd',scale='utc')
	return ttjdmean, ttjdmeanutc.isot
#-------------------------------------------------------------
def centertime2(imlist) :	
	'''
	CALCULATE MEAN DATE-OBS
	RETURN JD AND UTC
	'''
	obsdt	= []
	for inim in imlist:
		header	=fits.getheader(inim)
		hobsdate=header['DATE-OBS']
		htimeobs=header['TIME-OBS']
		obsdt.append(hobsdate+'T'+htimeobs)
	tt			= Time(obsdt, format='isot', scale='utc')
	ttjd		=tt.jd
	ttjdmean	= np.mean(ttjd)
	ttjdmeanutc	= Time(ttjdmean,format='jd',scale='utc')
	return ttjdmean, ttjdmeanutc.isot
#-------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
os.system('ls *.fits')
imlist	= glob.glob(raw_input('IMAGES TO COMBINE\t: ')); imlist.sort()
imcombine(imlist)

'''
#------------------------------------------------------------
#	MULTI PROCESSING
#------------------------------------------------------------
if __name__ == '__main__':
	jobs	= []
	p			= mp.Process(target=imcombine, args=imlist)
	jobs.append(p)
	p.start()
'''
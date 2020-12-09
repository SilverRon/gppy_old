#============================================================
#	python script for gregister using Alipy code
#	iraf-imcombine.py for combine "gregistered" images
#	Usage :
#	python gregister-script.py str1 ref_image
#	OR JUST RUN AND WRITE INPUT & REF IMAGE
#	2017.12.22	Changsu Choi
#	2019.05.08	GREGORY S.H. PAEK
#============================================================
import os, sys, glob
import alipy
import time
from pyraf import iraf
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
#============================================================
def gregistering(images_to_align, ref_image):
	starttime = time.time()
	if ref_image == '':
		ref_image = images_to_align[0]
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications:  # list of the same length as images_to_align.
		if id.ok == True:  # i.e., if it worked
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
		else:
			print "%20s : no transformation found !" % (id.ukn.name)
	outputshape = alipy.align.shape(ref_image)
	for id in identifications:
		if id.ok == True:
			params_align = dict(	filepath=id.ukn.filepath,
                            uknstarlist=id.uknmatchstars,
                            refstarlist=id.refmatchstars,
                            shape=alipy.align.shape(ref_image),
                            outdir='./',
                            makepng=False)
			alipy.align.irafalign(**params_align)
	deltime = time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#-------------------------------------------------------------
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
	return comin
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
#	MAIN COMMAND
#-------------------------------------------------------------
try:
	imkey = sys.argv[1]
except:
	os.system('ls *.fits')
	imkey = raw_input('IMAGES (Calib*0.fits) :\t')

if imkey == '':
	imkey = 'Calib*0.fits'
else:
	pass

imlist = glob.glob(imkey); imlist.sort()
#	DATE
datelist = []
for inim in imlist:
	part = inim.split('-')
	datelist.append(part[3])
datelist = list(set(datelist))
#	GREGISTERING & IMCOMBINE DAY BY DAY
comlists = []
for date in datelist:
	# sublist = glob.glob('Calib-*-{}-*.fits'.format(date)); sublist.sort()
	# sublist = sorted(glob.glob())
	sublist = []
	for inim in imlist:
		if date in inim:
			sublist.append(inim)
	refim = sublist[0]
	gregistering(sublist, refim)
	comlist = sorted(glob.glob('*-{}-*_gregister.fits'.format(date)))
	comlists.append(imcombine(comlist))
#	CLEAN
os.system('rm imcombine.list Calib-*_gregister.fits')

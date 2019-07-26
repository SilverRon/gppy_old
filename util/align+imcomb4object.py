#============================================================
#	Python CODE FOR GREGISTERING + IMCOMBINE
#	JUST RUN AND WRITE INPUT & REF IMAGE
#	2017.12.22	Changsu Choi
#	2019.05.08	GREGORY S.H. PAEK
#	2019.05.21	GREGORY S.H. PAEK
#============================================================
import os, sys, glob
import alipy
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from pyraf import iraf
import numpy as np
import time
#============================================================
def gregistering(images_to_align, ref_image):
	starttime	= time.time()
	if ref_image == '': ref_image = images_to_align[0]
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications: # list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
		else:
			print "%20s : no transformation found !" % (id.ukn.name)
	outputshape = alipy.align.shape(ref_image)
	for id in identifications:
		if id.ok == True:
			params_align	= dict(	filepath	= id.ukn.filepath,
									uknstarlist	= id.uknmatchstars,
									refstarlist	= id.refmatchstars,
									shape		= alipy.align.shape(ref_image),
									outdir		= './',
									makepng		= False)
			alipy.align.irafalign(**params_align)
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#-------------------------------------------------------------
def imcombine(imlist, timemode):
	refim		= imlist[0]
	starttime	= time.time()
	print('LIST '+str(len(imlist)))
	f	= open('imcombine.list', 'w')
	for inim in imlist:
		print(inim)
		f.write(inim+'\n')
	f.close()
	if timemode == 0:
		jdmean, utmean	= centertime(imlist)
		utpart			= utmean.split('T')[1].split(':')
		utpart[2]		= str(int(float(utpart[2])))
		utformat		= ''.join(utpart)
		compart			= refim.split('-')
		compart[0], compart[4], compart[6]	= 'Calib', utformat, compart[6][:-5]+'-com.fits'
		comin			= '-'.join(compart)
	if timemode == 1:
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
def align_imcomb_seq(imlist, refim, timemode):
	gregistering(imlist, refim)
	alignlist	= []
	for inim in imlist:
		alignlist.append(inim[:-5]+'_gregister.fits')
	imcombine(alignlist, timemode)
#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
os.system('ls -C *.fits')
str1	= raw_input('IMAGES TO PROCESS (Cal*.fits)\t: ')
if str1 == '':
	imlist	= glob.glob('Cal*.fits')
else:
	imlist	= glob.glob(str1)
#------------------------------------------------------------
#	ORGANIZE
#------------------------------------------------------------
objlist	= []
for inim in imlist:
	objlist.append(inim.split('-')[2])
objlist	= list(set(objlist)); objlist.sort()
#------------------------------------------------------------
#	MULTI PROCESSING
#------------------------------------------------------------
singlelist	= []
imfail		= []
for obj in objlist:
	prolist	= glob.glob('*'+obj+'*.fits')
	if len(prolist) > 1:
		try:
			#timemode =0 for others, =1 for SQUEAN
			align_imcomb_seq(prolist, prolist[0], timemode=1)
			#os.system('mv *'+obj+'*.fits multi_images')
		except:
			imfail.append(prolist)
	else:
		singlelist.append(obj)
print('SINGLE IMAGE LIST')
for single in singlelist:
	print(single)
	'''
	if len(prolist) > 1:
		if __name__ == '__main__':
			jobs	= []
			p			= mp.Process(target=align_imcomb_seq, args=(prolist, prolist[0]))
			jobs.append(p)
			p.start()
	'''
#============================================================
#	GREGISTERING+IMCOMBINE DAY BY DAY WITHIN 30 MIN
#	2019.06.24.	CREATED BY Gregory S.H. Paek
#============================================================
from astropy.io import fits
import numpy as np
import os, sys, glob
import alipy
import time
import astropy.units as u
from astropy.time import Time
from pyraf import iraf
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
#============================================================
#	FUNCTION
#------------------------------------------------------------
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
	
	jdmean, utmean	= centertime2(imlist)
	utpart			= utmean.split('T')[1].split(':')
	utpart[2]		= str(int(float(utpart[2])))
	utformat		= ''.join(utpart)
	compart			= refim.split('-')
	compart[0], compart[4], compart[6]	= 'Calib', utformat, compart[6][:-5]+'-com.fits'
	comin			= '-'.join(compart)
	'''
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
#------------------------------------------------------------
def hotpants(imlist, refim):
	starttime	= time.time()
	for inim in imlist:
		outfile = 'hd' + inim
		convfile= 'hc' + inim
		#com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -v 0 -inim '+infile[n]+' -tmplim '+ref_img+' -outim '+outfile+' -oci '+convfile
		#com='hotpants -v 0 -c i -n i -inim '+infile[n]+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
		#com='hotpants -c t -n i -inim '+infile[n]+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
		#com     = 'hotpants -c t -n i -iu 60000 -tl -40yu0 -tu 1000000 -v 0 -inim '+infile[n]+' -tmplim '+ref_img+' -outim '+outfile+' -oci '+convfile
		com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -tl -10000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
		os.system(com)
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#============================================================
#	INPUT
#------------------------------------------------------------
os.system('ls *.fits')
imlist		= glob.glob(raw_input('IMAGES\t: ')); imlist.sort()
refim		= raw_input('REF. IMAGE\t: ')
comlists	= []
#------------------------------------------------------------
#	CLASSIFY EPOCH
#------------------------------------------------------------
for inim in imlist:
	comlist	= []
	jdbase	= fits.getheader(inim)['mjd']
	for inim in imlist:
		jd	= fits.getheader(inim)['mjd']
		deljd	= np.abs(jd - jdbase)
		if deljd < 30/1440.:				# x/1440 : within x min
			#imlist.remove(inim)
			comlist.append(inim)
		else:
			pass
	comlists.append(comlist)
#------------------------------------------------------------
#	GERGISTERING + COMBINE
#------------------------------------------------------------
for comlist in comlists:
	reglist		= []
	gregistering(comlist, refim)
	for inim in comlist:
		reglist.append(inim[:-5]+'_gregister.fits')
	imcombine(reglist)
#------------------------------------------------------------
#	GREGISTERING + HOTPANTS
#------------------------------------------------------------
gregistering(glob.glob('C*com.fits'), refim)
if __name__ == '__main__':
	jobs	= []
	p			= mp.Process(target=hotpants, args=(glob.glob('C*com*ter.fits'), refim))
	jobs.append(p)
	p.start()


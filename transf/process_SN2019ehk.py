#============================================================
#	SN2019ehk/SN2018ein/NGC4321/NGC4321 FIELD
#	REGISTER+COMBINE/REGISTER(/SUBTRACTION/PHOTOMETRY/PLOT)
#	2019.05.06	MADE BY Gregory S.H. Paek
#============================================================
import os, sys, glob
import alipy
from astropy.io import ascii
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, Column
from astropy.time import Time
from pyraf import iraf
import os, sys
import numpy as np
import matplotlib.pyplot as plt
#------------------------------------------------------------
date		= '20190509'
obs			= raw_input('OBS\t: ')
if obs == 'LOAO':
	path_qso	= '/data3/IMSNG/IMSNGgalaxies/NGC4321/LOAO'
	path_gundam	= '/mnt/window/Users/User/Downloads/data/Project/SN2019ehk/LOAO'
if obs == 'DOAO':
	path_qso	= '/data3/IMSNG/IMSNGgalaxies/NGC4321/DOAO'
	path_gundam	= '/mnt/window/Users/User/Downloads/data/Project/SN2019ehk/DOAO'
if obs == 'SOAO':
	path_qso	= '/data3/IMSNG/IMSNGgalaxies/NGC4321/SOAO'
	path_gundam	= '/mnt/window/Users/User/Downloads/data/Project/SN2019ehk/SOAO'

path_ref	= '/mnt/window/Users/User/Downloads/data/Project/SN2019ehk/ref'
#------------------------------------------------------------
def imcombine(group, output):
	iraf.imcombine(group,output=output,combine="median",project="no",reject="none",scale="none",zero="mode")
#-------------------------------------------------------------
def centertime(imlist) :	
	'''
	CALCULATE MEAN 
	DATE-OBS
	RETURN JD AND UTC
	'''
	obsdt	= []
	for inim in imlist:
		header	=fits.getheader(inim)
		hobsdate=header['DATE-OBS']
		obsdt.append(hobsdate)
	tt			= Time(obsdt, format='isot', scale='utc')
	ttjd		=tt.jd
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
#	DOWNLOAD FILES FROM QSO SERVER
#------------------------------------------------------------
# downcom		=
print(downcom);os.system(downcom)
#------------------------------------------------------------
#	ARANGE FILES FOR EACH BANDS
#------------------------------------------------------------
for band in ['B', 'V', 'R', 'I']:
	path	= path_gundam+'/'+band
	imform= 'Calib-'+obs+'*'+date+'*-'+band+'-*0.fits'
	os.system('mkdir '+path+'/')
	mvcom	= 'mv '+path_gundam+'/'+imform+' '+path+'/'
	#mvcom	= 'mv '+path_gundam+'/Calib*'+band+'*'+date+'*0.fits '+path+'/'
	print(mvcom);os.system(mvcom)
	#------------------------------------------------------------
	#	GREGISTERING
	#------------------------------------------------------------
	imlist	= glob.glob(path+'/'+imform)
	for n in range(len(imlist)/3):
		sublist	= imlist[n*3:n*3+3]
		print('LIST '+str(n))
		f	= open(path+'/imcomb.list', 'w')
		for inim in sublist:
			print(inim[:-5]+'_gregister.fits')
			f.write(inim[:-5]+'_gregister.fits'+'\n')
		f.close()
		refim	= sublist[0]
		identifications = alipy.ident.run(refim, sublist, visu=False)
		for id in identifications:	# list of the same length as sublist.
			if id.ok == True:		# i.e.) if it worked
				print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
				params_align	= dict(	filepath	= id.ukn.filepath,
										uknstarlist	= id.uknmatchstars,
										refstarlist	= id.refmatchstars,
										shape		= alipy.align.shape(refim),
										outdir		= path,
										makepng		= False)
				#alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, makepng=False)
				alipy.align.irafalign(**params_align)
			else:
				print "%20s : no transformation found !" % (id.ukn.name)
		#------------------------------------------------------------
		#	MEDIAN COMBINE
		#------------------------------------------------------------
		if obs != 'SOAO':
			jdmean, utmean	= centertime(sublist)
		if obs == 'SOAO':
			jdmean, utmean	= centertime2(sublist)
		utpart			= utmean.split('T')[1].split(':')
		utpart[2]		= str(int(float(utpart[2])))
		utformat		= ''.join(utpart)

		compart			= refim.split('-')
		compart[0], compart[4], compart[6]	= 'Calib', utformat, 'com.fits'
		comin			= '-'.join(compart)

		imcombine('@'+path+'/imcomb.list', path+'/'+comin)

		comdata, comhdr	= fits.getdata(path+'/'+comin, header=True)
		comhdr['DATE-OBS']	= utmean
		fits.writeto(comin, comdata, comhdr, clobber=True)
		os.system('mkdir '+path+'/com/')
		os.system('mv '+path+'/Calib-*com.fits '+path+'/com/')
		os.system('rm '+path+'/imcomb.list')
		os.system('rm '+path+'/Calib-*_gregister.fits')
os.system('rm '+path_gundam+'/Calib-*com.fits')


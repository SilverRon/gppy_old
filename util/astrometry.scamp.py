#	WHAT SOURCES ARE IN FOV?
#	MADE BY Gu Lim					(2019.??.??)
#	MODIFIED BY Gregory S.H. Paek	(2019.07.25)
#============================================================
#	MODULE
#============================================================
import os, glob
#============================================================
#	FUNCTION
#============================================================
def wcsreset(inim,wcschar):
	from pyraf import iraf
	iraf.wcsreset.image = inim
	iraf.wcsreset.wcs   = wcschar
	iraf.wcsreset.mode  = 'h'
	iraf.wcsreset()
#------------------------------------------------------------
def hedit(inim,kword,value):
	from pyraf import iraf
	iraf.hedit.images = inim
	iraf.hedit.fields = kword
	iraf.hedit.value  = value
	iraf.hedit.add    = 'yes'
	iraf.hedit.verify = 'no'
	iraf.hedit.update = 'yes'
	iraf.hedit.mode   = 'h'
	iraf.hedit()
#------------------------------------------------------------
def scampastrom(inim, path_conf, param_scamp, cdmatrix, refcat='SDSS-R9'):
	from pyraf import iraf
	import numpy as np
	import astropy.units as u
	from astropy.io import fits
	from astropy.io import ascii
	import astropy.coordinates as coord
	#------------------------------------------------------------
	#	CONFIGUTRATION
	#------------------------------------------------------------
	sexconfig   = path_conf+'/astrom.sex'
	sexconv     = path_conf+'/astrom.conv'
	sexnnw      = path_conf+'/astrom.nnw'
	sexparam    = path_conf+'/astrom.param'
	psfconfig   = path_conf+'/prepsfex.sex'
	scampconfig = path_conf+'/astrom/astrom.scamp'
	#------------------------------------------------------------
	cd11, cd12, cd21, cd22 = cdmatrix['cd11'], cdmatrix['cd12'], cdmatrix['cd21'], cdmatrix['cd22']
	ra, dec, xpix, ypix = param_scamp['ra'], param_scamp['dec'], param_scamp['xpix'], param_scamp['ypix'] 
	threshold, minarea = param_scamp['threshold'], param_scamp['minarea']
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	rad = coord.Angle(ra,unit=u.hour)	
	radd = rad.degree
	decd = coord.Angle(dec,unit=u.deg)
	decdd = decd.degree
	outcat = inim[:-5]+'.astrosex.cat'
	outascii = outcat[:-4]+'.ascii'
	outhead = outcat[:-4]+'.head'
	#------------------------------------------------------------
	sexcom = 'sex -c '+sexconfig +' '+inim+' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH '+str(threshold) +' -DETECT_MINAREA '+str(minarea)+' -CATALOG_NAME '+outcat+' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+' -CATALOG_TYPE FITS_LDAC'
	# scampcom ='scamp -c '+scampconfig+' astromtest.cat '+ '-ASTREF_CATALOG '+refcat+' -CHECKPLOT_DEV PNG'
	scampcom = 'scamp -c {} {} -ASTREF_CATALOG {} -CHECKPLOT_DEV PNG'.format(scampconfig, outcat, refcat)
	#------------------------------------------------------------
	#	IRAF SETTING
	#------------------------------------------------------------
	wcsreset(inim,'physical')
	wcsreset(inim,'world')
	hedit(inim,'WAT0_001','system=image')
	hedit(inim,'WAT1_001', 'wtype=tan axtype=ra')
	hedit(inim,'WAT2_001', 'wtype=tan axtype=dec')
	hedit(inim,'RADECSYS', 'FK5')
	hedit(inim,'EQUINOX', '2000.')
	hedit(inim,'CTYPE1', 'RA---TAN')
	hedit(inim,'CTYPE2', 'DEC--TAN')
	hedit(inim,'CRVAL1', str(radd))
	hedit(inim,'CRVAL2', str(decdd))
	hedit(inim,'CRPIX1', str(xpix))
	hedit(inim,'CRPIX2', str(ypix))
	hedit(inim,'CD1_1', str(cd11))
	hedit(inim,'CD1_2', str(cd12))
	hedit(inim,'CD2_1', str(cd21))
	hedit(inim,'CD2_2', str(cd22))
	#------------------------------------------------------------
	os.system(sexcom)
	os.system('ldactoasc {} > {}'.format(outcat, outascii))		# BINARY TO ASCII
	sextbl = ascii.read(outascii)
	dolist, nolist = [], []
	os.system('mkdir bad_align')
	if len(sextbl['NUMBER']) < 5 : 
		print('Matched star number = ' + str(len(sextbl['NUMBER'])))
		print('I will not touch this, it has stars less than 5. \n')
		nolist.append(inim+'\n')
		os.system('mv '+inim+' bad_ailgn/')
	else :
		print('Scamp will do this!! \n')
		dolist.append(inim+'\n')
		os.system(scampcom)
		hdr1 = hdr
		hdr1.extend(fits.Header.fromtextfile(outhead), update=True, update_first=True)
		fits.writeto('sc'+inim, fits.getdata(inim),hdr1)
		print(inim,'astrometry work using SCAMP')
#============================================================
#	INITIAL SETTING
#============================================================
os.system('ls *.fits')
imlist = glob.glob(raw_input('IMAGE TO PROCESS (*.fits)\t: ')); imlist.sort()
path_conf = '/home/sonic/Research/yourpy/config/conf_scamp'
refcat = 'SDSS-R9'
# refcat = '2MASS'
# refcat = 'ASTREFCAT_GAIADR2'
"""
# USNO-A1,USNO-A2,USNO-B1,
# GSC-1.3,GSC-2.2,GSC-2.3,
# TYCHO-2, UCAC-1,UCAC-2,UCAC-3,UCAC-4,
# NOMAD-1, PPMX, CMC-14, 2MASS, DENIS-3,
# SDSS-R3,SDSS-R5,SDSS-R6,SDSS-R7,
# SDSS-R8, SDSS-R9
"""
#------------------------------------------------------------
#	FOR SQUEAN
cd11, cd12 = -7.95143538356E-05, -1.51742377918E-06
cd21, cd22 = -1.56550201495E-06, 7.9522048878E-05
cdmatrix = dict(cd11=cd11, cd12=cd12, cd21=cd21, cd22=cd22)
#------------------------------------------------------------
#	SET ra, dec, xpix, ypix WITH BRIGHT STAR (ds9 IMAGE(xpix, ypix), Aladin WEB SITE(ra, dec), ...)
#------------------------------------------------------------
# ra, dec = '14:14:29.903', '+36:23:47.62'
# xpix, ypix = 468.03027, 494.91935
# ra, dec = '14:10:00.843', '+38:43:06.09'
# xpix, ypix = 729.74037, 740.72891
# ra, dec = '13:46:02.824','+55:42:03.36'
# xpix, ypix = 685.94074, 347.76679
ra, dec = '13:44:50.599', '+55:52:20.08'
xpix, ypix = 457.02794, 621.24988

threshold, minarea = 2,4
param_scamp = dict(	ra=ra, dec=dec,
					xpix=xpix, ypix=ypix,
					threshold=threshold, minarea=minarea)
#------------------------------------------------------------
for inim in imlist:
	os.system('delwcs {}'.format(inim))
	if 'sc'+inim in glob.glob('sc'+inim): os.system('rm sc'+inim)
	scampastrom(inim, path_conf, param_scamp, cdmatrix, refcat)
#------------------------------------------------------------

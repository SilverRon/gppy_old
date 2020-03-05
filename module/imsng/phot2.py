#	PHOTOMETRY CODE FOR PYTHON 3.X
#	CREATED	2019.11.25	Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from imsng import phot
import time
#============================================================
#	FUNCTION
#============================================================
def getccdinfo(obs, path_obs):
	'''
	GET CCD INFORMATION (GAIN, PIXEL SCALE, FOV)

	gain, pixscale, fov = getccdinfo(obs, path_obs)
	
	INPUT:
	obs = 'SAO'
	path_obs = '/home/sonic/Research/table'

	OUTPUT:
	gain, pixscale, fov
	'''
	obstbl = ascii.read(path_obs+'/obs.txt')
	indx_obs = np.where(obstbl['obs']==obs)
	gain = obstbl[indx_obs]['gain'][0]
	pixscale = obstbl[indx_obs]['pixelscale'][0]
	fov = obstbl[indx_obs]['fov'][0]
	return gain, pixscale, fov
#------------------------------------------------------------
def secom(inim, param_insex, dualmode=False):
	'''
	
	'''
	param_sex = dict(	CONF_NAME = 'default.sex',
						#------------------------------
						#	CATALOG
						#------------------------------
						CATALOG_NAME = 'test.cat',
						CATALOG_TYPE = 'ASCII_HEAD',
						PARAMETERS_NAME = 'default.param',
						#------------------------------
						#	EXTRACTION
						#------------------------------
						DETECT_TYPE = 'CCD',
						DETECT_MINAREA = '5',
						DETECT_MAXAREA = '0',
						DETECT_THRESH = '1.5',
						# ANALYSIS_THRESH = 'RELATIVE',
						ANALYSIS_THRESH = '1.5',						
						FILTER = 'Y',
						FILTER_NAME = 'default.conv',
						DEBLEND_NTHRESH = '64',
						DEBLEND_MINCONT = '0.0001',
						CLEAN = 'Y',
						CLEAN_PARAM = '1.0',
						MASK_TYPE = 'CORRECT',
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						PHOT_APERTURES = '3',
						PHOT_AUTOPARAMS = '2.5,3.5',
						PHOT_PETROPARAMS = '2.0,3.5',
						SATUR_LEVEL  = '50000.0',
						SATUR_KEY = 'SQTURATE',
						MAG_ZEROPOINT = '0.0',
						MAG_GAMMA = '4.0',
						GAIN = '1.0',
						GAIN_KEY = 'GAIN',   
						PIXEL_SCALE = '1.0',
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = '3.0',
						STARNNW_NAME = 'default.nnw',
						#------------------------------
						#	BACKGROUND
						#------------------------------
						BACK_SIZE = '128',
						BACK_FILTERSIZE = '10',
						BACKPHOTO_TYPE = 'LOCAL',
						#------------------------------
						#	CHECK IMAGE
						#------------------------------
						CHECKIMAGE_TYPE = 'NONE',
						CHECKIMAGE_NAME = 'check.fits',
						#==============================
						#	MEMORY & MISCELLANEOUS
						#==============================
						MEMORY_OBJSTACK = '3000',
						MEMORY_PIXSTACK = '300000',
						MEMORY_BUFSIZE = '1024',
						VERBOSE_TYPE = 'NORMAL',
						HEADER_SUFFIX = '.head',
						WRITE_XML = 'N',
						XML_NAME = 'sex.xml')

	for key in param_insex.keys():
		param_sex[key] = param_insex[key]

	
	secom_normal = 'sex -c {} {} '.format(param_sex['CONF_NAME'], inim)
	secom_dual = 'sex -c {} {} '.format(param_sex['CONF_NAME'], inim)
	for key in param_sex.keys():
		if key != 'CONF_NAME':
			secom_normal += '-{} {} '.format(key, param_sex[key])

	# print(secom_normal)
	os.system(secom_normal)
#------------------------------------------------------------
def querybox(refcatname, obj, racent, decent, path_refcat, radius=0.5):
	'''
	reftbl = querybox(**param_query)
	'''
	#------------------------------------------------------------
	#	REF. CATALOG QUERY
	#------------------------------------------------------------
	refcatlist	= glob.glob(path_refcat+'/*.cat')
	#------------------------------------------------------------
	if refcatname	== 'PS1':
		if path_refcat+'/ps1-'+obj+'.cat' not in refcatlist:
			querytbl = phot.ps1_query(obj, racent, decent, path_refcat, radius=radius)
		else:
			querytbl = ascii.read(path_refcat+'/ps1-'+obj+'.cat')
		reftbl, refcat = phot.ps1_Tonry(querytbl, obj)
	#------------------------------------------------------------
	elif refcatname	== 'SDSS':
		if path_refcat+'/sdss-'+obj+'.cat' not in refcatlist:
			querytbl = phot.sdss_query(obj, racent, decent, path_refcat)
		else:
			querytbl = ascii.read(path_refcat+'/sdss-'+obj+'.cat')
		reftbl, refcat = phot.sdss_Blaton(querytbl, obj)
	#------------------------------------------------------------
	elif refcatname == 'APASS':
		if path_refcat+'/apass-'+obj+'.cat' not in refcatlist:
			querytbl = phot.apass_query(obj, racent, decent, path_refcat)
		else:
			querytbl = ascii.read(path_refcat+'/apass-'+obj+'.cat')
		reftbl, refcat = phot.apass_Blaton(querytbl, obj)
	#------------------------------------------------------------
	elif refcatname	== '2MASS':
		if path_refcat+'/2mass-'+obj+'.cat' not in refcatlist:
			querytbl        = phot.twomass_query(obj, racent, decent, path_refcat, band=refmagkey, radius=radius)
		else:
			querytbl = ascii.read(path_refcat+'/2mass-'+obj+'.cat')
		reftbl, refcat = querytbl, '2mass-'+obj+'.cat'
	return reftbl
#============================================================
#	USER SETTING
#============================================================
path_base = './'
radius = 0.5								#	[DEGREE]
frac = 0.9									#	IMAGE FRACTION [%]
refcatname = '2MASS'						#	REFERENCE CATALOG
inmagkey = 'MAG_APER'
inmagerkey = 'MAGERR_APER'
refmaglower, refmagupper = 10, 16.0			#	REF MAG RANGE [MAG]
refmagerupper = 0.1
inmagerupper = 0.1
# aperture = peeing*2
# aperture = '3'
aperture = str(input('aperture:'))			#	DIAMETER
# tra, tdec = 44.54382542, -8.9577325		#	GRB190829A (MORE ACCURATE)
tra, tdec = 44.62963541, -8.951577778		#	GRB190829A (MORE ACCURATE)
phottype = 'normal'
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_refcat	= '/home/sonic/Research/cat/refcat'
path_obs = '/home/sonic/Research/table'
path_config = '/home/sonic/Research/yourpy/gppy/config'
#============================================================
#	PHOTOMETRY
#============================================================
# inim = '/data1/GRB190829A/ukirt/test2/Calib-UKIRT-GRB190829A-20190830-123028-K-10.fits'
inim = '/data1/Test/Calib-DOAO-NGC5350-20190503-145856-R-60.fits'

head = os.path.splitext(inim)[0]

hdr = fits.getheader(inim)
xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
#------------------------------------------------------------
#	RA, Dec CENTER FOR QUERYING
#------------------------------------------------------------
try:
	w = WCS(inim)
	racent, decent = w.all_pix2world(xcent, ycent, 1)
	racent, decent = np.asscalar(racent), np.asscalar(decent)
except:
	print('BAD WCS INFORMATION?')
	racent, decent = hdr['CRVAL1'], hdr['CRVAL2']
#------------------------------------------------------------
#	DATE-OBS, JD
#------------------------------------------------------------
try:
	date_obs = hdr['date-obs']
	jd = round(Time(date_obs, format='isot', scale='utc').jd, 3)
except:
	date_obs = None
	jd = None
#------------------------------------------------------------
#	NAME INFO
#------------------------------------------------------------
part = inim.split('-')
obs = part[1]
obj = part[2]
refmagkey = part[5]
refmagerkey = refmagkey+'err'
#------------------------------------------------------------
gain, pixscale, fov = getccdinfo(obs, path_obs)
#------------------------------------------------------------
cat = head+'.cat'
seg = head+'.seg.fits'
bkg = head+'.bkg.fits'
sub = head+'.sub.fits'
psf = head+'.psf'
aper = head+'.aper.fits'

param = path_config+'/gregoryphot.param'
conv = path_config+'/gregoryphot.conv'
nnw = path_config+'/gregoryphot.nnw'
conf = path_config+'/gregoryphot.sex'
#------------------------------------------------------------
#	SOURCE EXTRACTOR CONFIGURATION
#------------------------------------------------------------
param_insex = dict(	#------------------------------
					#	CATALOG
					#------------------------------
					CATALOG_NAME = cat,
					#------------------------------
					#	CONFIG FILES
					#------------------------------
					CONF_NAME = conf,
					PARAMETERS_NAME = param,
					FILTER_NAME = conv,    
					STARNNW_NAME = nnw,
					#------------------------------
					#	EXTRACTION
					#------------------------------			
					# PSF_NAME = psf,
					DETECT_MINAREA = '5',
					DETECT_THRESH = '1.5',
					DEBLEND_NTHRESH = '64',
					DEBLEND_MINCONT = '0.0001',
					#------------------------------
					#	PHOTOMETRY
					#------------------------------
					PHOT_APERTURES = aperture,
					SATUR_LEVEL  = '65000.0',
					MAG_ZEROPOINT = '0.0',
					GAIN = str(gain),
					PIXEL_SCALE = str(pixscale),
					#------------------------------
					#	STAR/GALAXY SEPARATION
					#------------------------------
					SEEING_FWHM = '3.0',
					#------------------------------
					#	BACKGROUND
					#------------------------------
					BACK_SIZE = '128',
					BACK_FILTERSIZE = '10',
					BACKPHOTO_TYPE = 'LOCAL',
					#------------------------------
					#	CHECK IMAGE
					#------------------------------
					CHECKIMAGE_TYPE = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND',
					CHECKIMAGE_NAME = '{},{},{},{}'.format(seg, aper, bkg, sub),)

# peeing, seeing	= phot.psfex(inim, pixscale)
secom(inim, param_insex)
setbl = ascii.read(cat)

reftbl = querybox(refcatname, obj, racent, decent, path_refcat, radius=0.5)
#------------------------------------------------------------
#	CENTER POS. & DIST CUT
#------------------------------------------------------------
deldist = phot.sqsum((xcent-setbl['X_IMAGE']), (ycent-setbl['Y_IMAGE']))
indx_dist = np.where(deldist < np.sqrt(frac)*(xcent+ycent)/2.)
intbl = setbl[indx_dist]
intbl.write(cat, format='ascii', overwrite=True)
#------------------------------------------------------------
#	MATCHING
#------------------------------------------------------------
param_match = dict(	intbl=intbl, reftbl=reftbl,
					inra=intbl['ALPHA_J2000'], indec=intbl['DELTA_J2000'],
					refra=reftbl['ra'], refdec=reftbl['dec'],
					sep=2.0)
mtbl = phot.matching(**param_match)
#------------------------------------------------------------
#	ZEROPOINT CALCULATION
#------------------------------------------------------------
param_st4zp	= dict(	intbl=mtbl,
					inmagerkey=inmagkey,
					refmagkey=refmagkey,
					refmagerkey=refmagerkey,
					refmaglower=refmaglower,
					refmagupper=refmagupper,
					refmagerupper=refmagerupper,
					inmagerupper=inmagerupper)
param_zpcal	= dict(	intbl=phot.star4zp(**param_st4zp),
					inmagkey=inmagkey, inmagerkey=inmagerkey,
					refmagkey=refmagkey, refmagerkey=refmagerkey,
					sigma=2.0)
zp, zper, otbl, xtbl = phot.zpcal(**param_zpcal)
#------------------------------------------------------------
#	ZEROPOINT PLOT
#------------------------------------------------------------
phot.zpplot(outname='{}.zpcal.png'.format(head, inmagkey),
			otbl=otbl, xtbl=xtbl,
			inmagkey=inmagkey, inmagerkey=inmagerkey,
			refmagkey=refmagkey, refmagerkey=refmagerkey,
			zp=zp, zper=zper)
param_plot	= dict(	inim		= inim,
					numb_list	= otbl['NUMBER'],
					xim_list	= otbl['X_IMAGE'],
					yim_list	= otbl['Y_IMAGE'],
					add			= True,
					numb_addlist= xtbl['NUMBER'],
					xim_addlist	= xtbl['X_IMAGE'],
					yim_addlist	= xtbl['Y_IMAGE'])
try:
	phot.plotshow(**param_plot)
except:
	print('FAIL TO DRAW ZEROPOINT GRAPH')
	pass
#------------------------------------------------------------
#	TARGET PHOTOMETRY
#------------------------------------------------------------
skymean, skymed, skysig = phot.bkgest_mask(inim)
peeing, seeing	= np.median(intbl['FWHM_IMAGE']), np.median(intbl['FWHM_WORLD'])*3600

aper = 2*peeing
ul_3sig = phot.limitmag(3, zp, aper, skysig)
ul_5sig = phot.limitmag(5, zp, aper, skysig)
#------------------------------------------------------------
#	ADD HEADER INFO
#------------------------------------------------------------
phot.puthdr(inim, 'SEEING',	round(seeing, 3),		hdrcomment='SEEING [arcsec]')
phot.puthdr(inim, 'PEEING',	round(peeing, 3),		hdrcomment='SEEING [pixel]')
phot.puthdr(inim, 'SKYSIG',	round(skysig, 3),		hdrcomment='SKY SIGMA VALUE')
phot.puthdr(inim, 'SKYVAL',	round(skymed, 3),		hdrcomment='SKY MEDIAN VALUE')
phot.puthdr(inim, 'ZP',	round(zp, 3),			hdrcomment='2*SEEING DIAMETER')
phot.puthdr(inim, 'ZPERR',round(zper, 3),		hdrcomment='2*SEEING DIAMETER')
phot.puthdr(inim, 'UL3SIG',	round(ul_3sig, 3),			hdrcomment='2*SEEING 3 sigma limit mag')
phot.puthdr(inim, 'UL5SIG',	round(ul_5sig, 5),			hdrcomment='2*SEEING 5 sigma limit mag')
phot.puthdr(inim, 'STDNUMB',len(otbl),				hdrcomment='# OF STD STARS')
#------------------------------------------------------------
#	NORMAL PHOTOMETRY
#------------------------------------------------------------
if phottype == 'normal':
	intbl['REAL_'+inmagkey] = zp + intbl[inmagkey]
	intbl['REAL_'+inmagerkey] = phot.sqsum(zper, intbl[inmagerkey])
	indx_targ = phot.targetfind(tra, tdec, intbl['ALPHA_J2000'], intbl['DELTA_J2000'], sep=seeing)
	if indx_targ != None:
		mag, mager	= intbl[indx_targ]['REAL_'+inmagkey], intbl[indx_targ]['REAL_'+inmagerkey]
		radeg, dedeg = np.asscalar(intbl['ALPHA_J2000'][indx_targ]), np.asscalar(intbl['DELTA_J2000'][indx_targ])
	else:
		mag, mager	= -99, -99
#------------------------------------------------------------
#	SUBTRACTION PHOTOMETRY
#------------------------------------------------------------
elif phottype == 'subt':
	subim = os.path.dirname(inim)+'/hd'+os.path.basename(inim)
	subpsf = os.path.dirname(inim)+'/hd'+os.path.basename(inim)[:-5]+'.psf'
	subcat = os.path.dirname(inim)+'/hd'+os.path.basename(inim)[:-5]+'.cat'

	# os.system('cp {} {}'.format(psf, subim[:-5]+'.psf'))

	param_subsex= dict(	#------------------------------
						#	CATALOG
						#------------------------------
						CATALOG_NAME = subcat,
						#------------------------------
						#	CONFIG FILES
						#------------------------------
						CONF_NAME = conf,
						PARAMETERS_NAME = param,
						FILTER_NAME = conv,    
						STARNNW_NAME = nnw,
						#------------------------------
						#	EXTRACTION
						#------------------------------			
						# PSF_NAME = psf,
						DETECT_MINAREA = '5',
						DETECT_THRESH = '1.5',
						DEBLEND_NTHRESH = '64',
						DEBLEND_MINCONT = '0.0001',
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						PHOT_APERTURES = '5',
						SATUR_LEVEL  = '65000.0',
						MAG_ZEROPOINT = '0.0',
						GAIN = str(gain),
						PIXEL_SCALE = str(pixscale),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = '3.0',
						#------------------------------
						#	BACKGROUND
						#------------------------------
						BACK_SIZE = '128',
						BACK_FILTERSIZE = '10',
						BACKPHOTO_TYPE = 'LOCAL',
						#------------------------------
						#	CHECK IMAGE
						#------------------------------
						CHECKIMAGE_TYPE = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND',
						CHECKIMAGE_NAME = '{},{},{},{}'.format(seg, aper, bkg, sub),)

	secom(subim, param_subsex)

	subtbl = ascii.read(subcat)
	subtbl['REAL_'+inmagkey] = zp + subtbl[inmagkey]
	subtbl['REAL_'+inmagerkey] = phot.sqsum(zper, subtbl[inmagerkey])
	indx_targ = phot.targetfind(tra, tdec, subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'], sep=seeing)
	if indx_targ != None:
		mag, mager	= subtbl[indx_targ]['REAL_'+inmagkey], subtbl[indx_targ]['REAL_'+inmagerkey]
		radeg, dedeg = np.asscalar(intbl['ALPHA_J2000'][indx_targ]), np.asscalar(intbl['DELTA_J2000'][indx_targ])
	else:
		mag, mager	= -99, -99
#------------------------------------------------------------
#	CALC. DEPTH
#------------------------------------------------------------
elif phottype == 'depth':
	mag, mager	= -99, -99

onetbl	= Table([[inim], [obs], [obj], [round(radeg, 3)], [round(dedeg, 3)], [date_obs], [jd], [refmagkey], [len(otbl)], [round(zp, 3)], [round(zper, 3)], [round(seeing, 3)], [round(skymed, 3)], [round(skysig, 3)], [round(ul_3sig, 3)], [round(ul_5sig, 3)], [mag], [mager], [aperture]],
				names=('image', 'obs', 'obj', 'ra', 'dec', 'date-obs', 'jd', 'filter', 'stdnumb', 'zp', 'zper', 'seeing', 'skyval', 'skysig', 'ul_3sig', 'ul_5sig','mag', 'magerr', 'aperture'))
# return onetbl

# tblist.append(onetbl)

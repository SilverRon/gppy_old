#	PHOTOMETRY CODE FOR PYTHON 3.X
#	CREATED	2019.11.25	Gregory S.H. Paek
#	UPDATE 2020.01.17	Gregory S.H. Paek
#	UPDATE 2020.02.20	Gregory S.H. Paek
#============================================================
import os, glob, subprocess
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
from imsng import query
from imsng import tool
from imsng import gcurve
import time
from datetime import date
#============================================================
#	FUNCTION
#============================================================
def gcurvephot(inim, phottype, tra, tdec, path_base, path_refcat, path_obs, path_config, refmaglower, refmagupper, refmagerupper, inmagerupper, flagcut, apertures, aper_input, radius = 0.5, frac = 0.9, refcatname = 'PS1',):
	head = os.path.splitext(inim)[0]
	hdr = fits.getheader(inim)
	xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
	#------------------------------------------------------------
	#	RA, Dec CENTER FOR QUERYING
	#------------------------------------------------------------
	try:
		w = WCS(inim)
		racent, decent = w.all_pix2world(xcent, ycent, 1)
		racent, decent = racent.item(), decent.item()
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
	# refmagkey = 'B'
	refmagerkey = refmagkey+'err'
	#------------------------------------------------------------
	gain, pixscale, rdnoise = tool.getccdinfo(obs, path_obs)
	#------------------------------------------------------------
	cat = head+'.cat'
	cat_gc = head+'.gcurve.cat'
	seg = head+'.seg.fits'
	bkg = head+'.bkg.fits'
	sub = head+'.sub.fits'
	psf = head+'.psf'
	aper = head+'.aper.fits'
	#	GROWTH CURVE NAMES
	param_gc = path_config+'/growthcurve.param'
	conv_gc = path_config+'/growthcurve.conv'
	nnw_gc = path_config+'/growthcurve.nnw'
	conf_gc = path_config+'/growthcurve.sex'
	#	PHOTOMETRY NAMES
	param = path_config+'/gregoryphot.param'
	conv = path_config+'/gregoryphot.conv'
	nnw = path_config+'/gregoryphot.nnw'
	conf = path_config+'/gregoryphot.sex'
	#------------------------------------------------------------
	#	SOURCE EXTRACTOR CONFIGURATION FOR GROTH CURVE
	#------------------------------------------------------------
	param_gcurve = dict(#------------------------------
						#	CATALOG
						#------------------------------
						CATALOG_NAME = cat_gc,
						#------------------------------
						#	CONFIG FILES
						#------------------------------
						CONF_NAME = conf_gc,
						PARAMETERS_NAME = param_gc,
						FILTER_NAME = conv_gc,    
						STARNNW_NAME = nnw_gc,
						#------------------------------
						#	EXTRACTION
						#------------------------------			
						# PSF_NAME = psf,
						DETECT_MINAREA = '5',
						DETECT_THRESH = '3.0',
						DEBLEND_NTHRESH = '64',
						DEBLEND_MINCONT = '0.0001',
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						PHOT_APERTURES = aper_input,
						SATUR_LEVEL  = '65000.0',
						# MAG_ZEROPOINT = '0.0',
						GAIN = str(gain),
						PIXEL_SCALE = str(pixscale),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = '3.0',
						#------------------------------
						#	BACKGROUND
						#------------------------------
						BACK_SIZE = '32',
						BACK_FILTERSIZE = '10',
						BACKPHOTO_TYPE = 'LOCAL',
						#------------------------------
						#	CHECK IMAGE
						#------------------------------
						# CHECKIMAGE_TYPE = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND',
						# CHECKIMAGE_NAME = '{},{},{},{}'.format(seg, aper, bkg, sub),
						)
	print('1. GROWTH CURVE')
	os.system(phot.sexcom(inim, param_gcurve))
	setbl = ascii.read(cat_gc)
	# setbl.write(cat, format='ascii', overwrite=True)
	reftbl = query.querybox(refcatname, obj, racent, decent, path_refcat, radius=1.0, refmagkey=refmagkey)
	#	CENTER POS. & DIST CUT
	deldist = phot.sqsum((xcent-setbl['X_IMAGE']), (ycent-setbl['Y_IMAGE']))
	# indx_dist = np.where(deldist < np.sqrt(frac)*(xcent+ycent)/2.)
	indx_dist = np.where(deldist < frac*(xcent+ycent)/2.)
	intbl = setbl[indx_dist]
	# intbl.write(cat, format='ascii', overwrite=True)
	#	MATCHING
	param_match = dict(	intbl=intbl, reftbl=reftbl,
						inra=intbl['ALPHA_J2000'], indec=intbl['DELTA_J2000'],
						refra=reftbl['ra'], refdec=reftbl['dec'],
						sep=3.0)
	mtbl = phot.matching(**param_match)
	#	ZEROPOINT CALCULATION
	param_st4zp = dict(	intbl=mtbl,
						inmagerkey=inmagkey,
						refmagkey=refmagkey,
						refmagerkey=refmagerkey,
						refmaglower=refmaglower,
						refmagupper=refmagupper,
						refmagerupper=refmagerupper,
						inmagerupper=inmagerupper,
						flagcut=flagcut,
						)
	param_zpcal = dict(	intbl=phot.star4zp(**param_st4zp),
						inmagkey=inmagkey, inmagerkey=inmagerkey,
						refmagkey=refmagkey, refmagerkey=refmagerkey,
						sigma=2.0,
						# method='weightedmean',
						)
	zp, zper, otbl, xtbl = phot.zpcal(**param_zpcal)

	alltbl = vstack([otbl, xtbl])
	alltbl = alltbl[alltbl['NUMBER'].argsort()]

	seeing_input = str(np.median(alltbl['FWHM_WORLD'])*3600)
	# print('SEEING\t'+seeing_input)

	aperture, optapers = gcurve.gcurveplot(inim, alltbl, apertures)
	#------------------------------------------------------------
	#	SOURCE EXTRACTOR CONFIGURATION FOR PHOTOMETRY
	#------------------------------------------------------------
	# aperture = 50
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
						# DETECT_MINAREA = '3',
						# DETECT_THRESH = '1.0',
						# DETECT_THRESH = '1.25',
						# DETECT_THRESH = '1.3',
						DETECT_THRESH = '1.5',
						# DETECT_THRESH = '3.0',
						DEBLEND_NTHRESH = '64',
						DEBLEND_MINCONT = '0.00001',
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						PHOT_APERTURES = str(aperture),
						SATUR_LEVEL  = '65000.0',
						MAG_ZEROPOINT = '0.0',
						GAIN = str(gain),
						PIXEL_SCALE = str(pixscale),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						# SEEING_FWHM = '3.0',
						SEEING_FWHM = seeing_input,
						#------------------------------
						#	BACKGROUND
						#------------------------------
						BACK_SIZE = '64',
						# BACK_SIZE = '128',
						# BACK_SIZE = '256',
						BACK_FILTERSIZE = '3',
						BACKPHOTO_TYPE = 'LOCAL',
						#------------------------------
						#	CHECK IMAGE
						#------------------------------
						CHECKIMAGE_TYPE = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND',
						CHECKIMAGE_NAME = '{},{},{},{}'.format(seg, aper, bkg, sub),
						)
	print('2. SOURCE EXTRACTOR')
	com = phot.sexcom(inim, param_insex)
	sexout = subprocess.getoutput(com)
	line = [s for s in sexout.split('\n') if 'RMS' in s]
	skymed, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])

	# peeing, seeing	= phot.psfex(inim, pixscale)
	setbl = ascii.read(cat)
	#------------------------------------------------------------
	#	CENTER POS. & DIST CUT
	#------------------------------------------------------------
	deldist = phot.sqsum((xcent-setbl['X_IMAGE']), (ycent-setbl['Y_IMAGE']))
	# indx_dist = np.where(deldist < np.sqrt(frac)*(xcent+ycent)/2.)
	indx_dist = np.where(deldist < frac*(xcent+ycent)/2.)
	intbl = setbl
	# intbl.write(cat, format='ascii', overwrite=True)
	frctbl = setbl[indx_dist]
	#------------------------------------------------------------
	#	MATCHING
	#------------------------------------------------------------
	param_match = dict(	intbl=frctbl, reftbl=reftbl,
						inra=frctbl['ALPHA_J2000'], indec=frctbl['DELTA_J2000'],
						refra=reftbl['ra'], refdec=reftbl['dec'],
						sep=3.0)
						# sep=seeing_input)
	print('3. MATCHING')
	mtbl = phot.matching(**param_match)
	# mtbl.write(cat, format='ascii', overwrite=True)
	#------------------------------------------------------------
	#	ZEROPOINT CALCULATION
	#------------------------------------------------------------
	param_st4zp	= dict(	intbl=mtbl,
						inmagerkey=inmagerkey,
						refmagkey=refmagkey,
						refmagerkey=refmagerkey,
						refmaglower=refmaglower,
						refmagupper=refmagupper,
						refmagerupper=refmagerupper,
						inmagerupper=inmagerupper,
						verbose=True,
						plot=True,
						flagcut=flagcut,
						plotout=inim[:-5]+'.star.png')
	param_zpcal	= dict(	intbl=phot.star4zp(**param_st4zp),
						inmagkey=inmagkey, inmagerkey=inmagerkey,
						refmagkey=refmagkey, refmagerkey=refmagerkey,
						sigma=2.0,
						# method='weightedmean',
						)
	print('4. ZERO POINT CALCULATION')
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
	print('5. PLOT')
	try:
		phot.plotshow(**param_plot)
		plt.close('all')
	except:
		print('FAIL TO DRAW ZEROPOINT GRAPH')
		pass
	#------------------------------------------------------------
	#	TARGET PHOTOMETRY
	#------------------------------------------------------------
	print('6. BACKGROUND ESTIMATION')
	# skymean, skymed, skysig = phot.bkgest_mask(inim)
	peeing, seeing	= np.median(intbl['FWHM_IMAGE']), np.median(intbl['FWHM_WORLD'])*3600
	ellipticity = np.median(intbl['ELLIPTICITY'])

	aper = 2*peeing
	ul_3sig = phot.limitmag(3, zp, aper, skysig)
	ul_5sig = phot.limitmag(5, zp, aper, skysig)
	#------------------------------------------------------------
	#	ADD HEADER INFO
	#------------------------------------------------------------
	print('7. CHANGE HEADER')
	phot.puthdr(inim, 'PHOTIME',date.today().isoformat(),	hdrcomment='PHTOMETRY TIME [KR]')
	phot.puthdr(inim, 'SEEING',	round(seeing, 3),	hdrcomment='SEEING [arcsec]')
	phot.puthdr(inim, 'PEEING',	round(peeing, 3),	hdrcomment='SEEING [pixel]')
	phot.puthdr(inim, 'ELLIPTICITY', round(ellipticity, 3), hdrcomment='ELLIPTICITY [0-1]')
	phot.puthdr(inim, 'APERPIX', aperture,			hdrcomment='APERTURE DIAMETER [pixel]')
	phot.puthdr(inim, 'APER', round(aperture*pixscale, 3), hdrcomment='APERTURE DIAMETER [arcsec]')
	phot.puthdr(inim, 'NAPER', round(aperture/peeing, 3), hdrcomment='N = APERTURE/PEEING')
	phot.puthdr(inim, 'SKYSIG',	round(skysig, 3),	hdrcomment='SKY SIGMA VALUE')
	phot.puthdr(inim, 'SKYVAL',	round(skymed, 3),	hdrcomment='SKY MEDIAN VALUE')
	phot.puthdr(inim, 'REFCAT', refcatname,			hdrcomment='REFERENCE CATALOG NAME')
	phot.puthdr(inim, 'MAGLOW', refmaglower,		hdrcomment='REF MAG RANGE, BRIGHT LIMIT')
	phot.puthdr(inim, 'MAGUP', refmagupper,			hdrcomment='REF MAG RANGE, DIM LIMIT')
	phot.puthdr(inim, 'STDNUMB',len(otbl),			hdrcomment='# OF STD STARS')
	phot.puthdr(inim, 'ZP',	round(zp, 3),			hdrcomment='2*SEEING DIAMETER')
	phot.puthdr(inim, 'ZPERR',round(zper, 3),		hdrcomment='2*SEEING DIAMETER')
	phot.puthdr(inim, 'UL3SIG',	round(ul_3sig, 3),	hdrcomment='2*SEEING 3 sigma limit mag')
	phot.puthdr(inim, 'UL5SIG',	round(ul_5sig, 5),	hdrcomment='2*SEEING 5 sigma limit mag')
	print('8. PHOTOMETRY TABLE')
	#------------------------------------------------------------
	#	NORMAL PHOTOMETRY
	#------------------------------------------------------------
	if phottype == 'normal':
		intbl['REAL_'+inmagkey] = zp + intbl[inmagkey]
		intbl['REAL_'+inmagerkey] = phot.sqsum(zper, intbl[inmagerkey])
		indx_targ = phot.targetfind(tra, tdec, intbl['ALPHA_J2000'], intbl['DELTA_J2000'], sep=seeing)
		if indx_targ != None:
			# print(indx_targ)
			print(intbl['X_IMAGE'][indx_targ], intbl['Y_IMAGE'][indx_targ])
			mag, mager	= intbl[indx_targ]['REAL_'+inmagkey], intbl[indx_targ]['REAL_'+inmagerkey]
			radeg, dedeg = np.asscalar(intbl['ALPHA_J2000'][indx_targ]), np.asscalar(intbl['DELTA_J2000'][indx_targ])
			radeg, dedeg = intbl['ALPHA_J2000'][indx_targ].item(), intbl['DELTA_J2000'][indx_targ].item()		
		else:
			mag, mager	= -99, -99
			radeg, dedeg = -99, -99
	#------------------------------------------------------------
	#	SUBTRACTION PHOTOMETRY
	#------------------------------------------------------------
	elif phottype == 'subt':
		if os.path.dirname(inim) == '':
			interval = './'
		else:
			interval = '/'
		subim = os.path.dirname(inim)+interval+'hd'+os.path.basename(inim)
		subpsf = os.path.dirname(inim)+interval+'hd'+os.path.basename(inim)[:-5]+'.psf'
		subcat = os.path.dirname(inim)+interval+'hd'+os.path.basename(inim)[:-5]+'.cat'

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
							DETECT_THRESH = '3.0',
							DEBLEND_NTHRESH = '64',
							DEBLEND_MINCONT = '0.0000001',
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
							SEEING_FWHM = seeing_input,
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

		os.system(phot.sexcom(subim, param_subsex))

		subtbl = ascii.read(subcat)
		subtbl['REAL_'+inmagkey] = zp + subtbl[inmagkey]
		subtbl['REAL_'+inmagerkey] = phot.sqsum(zper, subtbl[inmagerkey])
		indx_targ = phot.targetfind(tra, tdec, subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'], sep=2*seeing)
		if indx_targ != None:
			mag, mager	= subtbl[indx_targ]['REAL_'+inmagkey], subtbl[indx_targ]['REAL_'+inmagerkey]
			radeg, dedeg = subtbl['ALPHA_J2000'][indx_targ].item(), subtbl['DELTA_J2000'][indx_targ].item()

		else:
			mag, mager	= -99, -99
			# radeg, dedeg = subtbl['ALPHA_J2000'][indx_targ].item(), subtbl['DELTA_J2000'][indx_targ].item()
			radeg, dedeg = tra, tdec
	#------------------------------------------------------------
	#	CALC. DEPTH
	#------------------------------------------------------------
	elif phottype == 'depth':
		mag, mager	= -99, -99
		radeg, dedeg = np.median(intbl['ALPHA_J2000']), np.median(intbl['DELTA_J2000'])
	onetbl	= Table([[inim], [obs], [obj], [radeg], [dedeg], [date_obs], [jd], [refmagkey], [len(otbl)], [round(zp, 3)], [round(zper, 3)], [round(seeing, 3)], [round(ellipticity, 3)], [round(skymed, 3)], [round(skysig, 3)], [round(ul_3sig, 3)], [round(ul_5sig, 3)], [round(mag, 3)], [round(mager, 3)], [round(aperture, 3)]],
					names=('image', 'obs', 'obj', 'ra', 'dec', 'date-obs', 'jd', 'filter', 'stdnumb', 'zp', 'zper', 'seeing', 'ellipticity', 'skyval', 'skysig', 'ul_3sig', 'ul_5sig','mag', 'magerr', 'aper_dia_pix'))

	tbl = ascii.read(cat)
	tbl['REAL_'+inmagkey] = zp + tbl[inmagkey]
	tbl['REAL_'+inmagerkey] = phot.sqsum(zper, tbl[inmagerkey])
	tbl.write(cat, format='ascii.tab', overwrite=True)

	return onetbl
#============================================================
#	USER SETTING
#============================================================
os.system('ls *.fits')
imkey = input('image to process\t: ')
if imkey == '':
	imkey = '*com.fits'
starttime = time.time()
imlist = glob.glob(imkey)
imlist.sort()
for img in imlist: print(img)
print(len(imlist))
#------------------------------------------------------------
# inim = '/data1/Test/Calib-DOAO-NGC5350-20190503-145856-R-60.fits'
path_base = './'
radius = 1.0								#	[DEGREE]
# frac = 1.0									#	IMAGE FRACTION [%]
# frac = 0.9
frac = 0.75
# frac = 0.5
refcatname = 'PS1'							#	REFERENCE CATALOG
# refcatname = 'APASS'
# refcatname = 'SDSS'
# refcatname = '2MASS'
#------------------------------------------------------------
inmagkey = 'MAG_APER'
inmagerkey = 'MAGERR_APER'
# inmagkey = 'MAG_AUTO'
# inmagerkey = 'MAGERR_AUTO'
#------------------------------------------------------------
# refmaglower, refmagupper = 12, 17			#	LOAO
# refmaglower, refmagupper = 12, 17			#	MAO
# refmaglower, refmagupper = 14, 17			#	REF MAG RANGE [MAG]
refmaglower, refmagupper = 12, 20			#	REF MAG RANGE [MAG]
# refmaglower, refmagupper = 12, 18			#	REF MAG RANGE [MAG]
# refmaglower, refmagupper = 10, 19			#	REF MAG RANGE [MAG] deep DOAO
# refmaglower, refmagupper = 0, 16.5			#	REF MAG RANGE [MAG]
# refmaglower, refmagupper = 12, 16.0			#	CBNUO
# refmaglower, refmagupper = 0, 18.5			#	DOAO
# refmaglower, refmagupper = 15, 25			#	DOAO
# refmagerupper = 0.1
refmagerupper = 0.05
inmagerupper = 0.05
#------------------------------------------------------------
# tra, tdec = 208.3714550, 40.2754194
# tra, tdec = 185.7288542, 15.8236250				#	SN2020oi
# tra, tdec = 161.6379008, 13.7418711				#	SN2018kp
# tra, tdec = 44.5438520, -8.9577875				#	GRB 190829A
# tra, tdec = 261.2277917, 31.4283333
# tra, tdec = 260.8090281, 14.3502257  # IceCube201021A
tra, tdec = 42.1846250, 12.1372444 # AT2020yxz
# tra, tdec = 94.1635667, -21.4998250 # AT2020zyy
#------------------------------------------------------------
phottype = 'normal'
# phottype = 'subt'
# phottype = 'depth'
#------------------------------------------------------------
# aperture = str(input('aperture:'))			#	DIAMETER
# apertures = np.arange(1, 16.5, 0.5)			#	LOAO
# apertures = np.arange(1, 17.0, 0.5)			#	LOAO
# apertures = np.arange(1, 20.0, 0.75)
# apertures = np.arange(5, 35.0, 1.0)			#	SOAO
# apertures = np.arange(5, 50, 1.5)			#	SOAO
# apertures = np.arange(1, 36.5, 1.5)			#	deep DOAO
# apertures = np.arange(6.25, 39.25, 1.25)			#	deep DOAO
# apertures = np.arange(30, 70, 1.5)			#	deep DOAO
# apertures = np.arange(10, 62, 2)			#//	deep DOAO
# apertures = np.arange(5, 40, 1.5)			#//	deep DOAO
apertures = np.arange(3, 25, 1.0)			#//	deep DOAO

aper_input = ''
for i in apertures:
	aper_input = aper_input+'{},'.format(i)
aper_input = aper_input[:-1]
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_refcat	= '/home/sonic/Research/cat/refcat'
path_obs = '/home/sonic/Research/table'
path_config = '/home/sonic/Research/yourpy/gppy/config'
#============================================================
#	MAIN COMMAND
#============================================================
tblist = []
plt.rcParams["figure.figsize"] = (10, 10)
# time.sleep(900)
for i, inim in enumerate(imlist):
	print('-'*60)
	print(inim, '[{}/{}]'.format(i+1, len(imlist)))
	print('-'*60)
	param_phot = dict(	
						inim = inim,
						phottype = phottype,
						tra = tra, tdec = tdec,

						radius = radius,
						frac = frac,

						path_base = path_base,
						path_refcat = path_refcat,
						path_obs = path_obs,
						path_config = path_config,

						refcatname = refcatname,
						refmaglower = refmaglower,
						refmagupper = refmagupper,
						refmagerupper = refmagerupper,
						inmagerupper = inmagerupper,
						flagcut = 0,
						# flagcut = 2,

						apertures = apertures,
						aper_input = aper_input,
					)
	'''
	for key in param_phot.keys():
		print(key, param_phot[key])
	'''
	# tblist.append(gcurvephot(**param_phot))
	try:
		tblist.append(gcurvephot(**param_phot))
	except:
		pass
	# os.system('rm *.aper.fits *.bkg.fits *.sub.fits')
	# os.system('rm psf*.fits snap*.fits *.xml *sub.fits *.aper.fits *.bkg.fits *.seg.fits')
#------------------------------------------------------------
#	FINISH
#------------------------------------------------------------
if len(tblist) == 0:
	print('PHOTOMETRY FAILED!')
else:
	photbl = vstack(tblist)
	# if 'phot.dat' in glob.glob(path_base+'/phot.dat'):
		# os.system('mv {} {}'.format(path_base+'/phot.dat', path_base+'/phot.dat.bkg'))
	'''
	#	phot.dat REPLACE
	photlist = glob.glob('phot*.dat')

	if 'phot.dat' in photlist:
		photnumb = 0
		phot_rpl = 'phot.{}.dat'.format(photnumb)
		while phot_rpl not in photlist:
			photnumb += 1
		
		com = 'mv phot.dat {}'.format(phot_rpl)
		print(com)
		os.system(com)
	'''
	photbl.write(path_base+'/phot.dat', format='ascii.tab', overwrite=True)

	print('='*60)
	print('All PROCESS IS DONE.\t('+str(round(time.time() - starttime, 2))+' sec)')
	print('-'*60)
	print(photbl['filter', 'date-obs', 'seeing', 'skyval', 'skysig', 'ul_5sig', 'mag', 'magerr'])

single = '''
inim = inim
phottype = phottype
tra = tra
tdec = tdec

radius = radius
frac = frac

path_base = path_base
path_refcat = path_refcat
path_obs = path_obs
path_config = path_config

refcatname = refcatname
refmaglower = refmaglower
refmagupper = refmagupper
refmagerupper = refmagerupper
inmagerupper = inmagerupper
flagcut = 0

apertures = apertures
aper_input = aper_input
'''

image_quality = '''
plt.close('all')
step = +0.1
med = np.median(photbl['ul_5sig'])
std = np.std(photbl['ul_5sig'])
indx = np.where(
				# (photbl['ul_5sig']>med+2*std) |
				(photbl['ul_5sig']<med-2*std)
				)

print(photbl[indx])

bins = np.arange(np.min(photbl['ul_5sig']), np.max(photbl['ul_5sig'])+step, step)
plt.hist(photbl['ul_5sig'], bins=bins, color='tomato', alpha=0.5)
plt.hist(photbl['ul_5sig'][indx], bins=bins, color='k', alpha=0.5)
plt.axvline(x=med, linestyle='--', color='red')
plt.show()
'''
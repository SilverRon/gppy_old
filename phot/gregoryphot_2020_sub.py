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
from imsng import query
from imsng import tool
from imsng import gcurve
import time
#============================================================
#	FUNCTION
#============================================================
def gcurvephot(inim, phottype, tra, tdec, path_base, path_refcat, path_obs, path_config, refmaglower, refmagupper, refmagerupper, inmagerupper, apertures, aper_input, radius = 0.5, frac = 0.9, refcatname = 'PS1',):
	#------------------------------------------------------------
	#	NAME INFO
	#------------------------------------------------------------
	head = os.path.splitext(inim)[0]
	part = inim.split('-')
	obs = part[1]
	obj = part[2]
	refmagkey = part[5]
	refmagerkey = refmagkey+'err'
	#------------------------------------------------------------
	gain, pixscale, fov = tool.getccdinfo(obs, path_obs)
	#------------------------------------------------------------
	#	SUBTRACTION PHOTOMETRY
	#------------------------------------------------------------
	if phottype == 'subt':
		if os.path.dirname(inim) == '':
			interval = './'
		else:
			interval = '/'
		subim = os.path.dirname(inim)+interval+'hd'+os.path.basename(inim)
		subpsf = os.path.dirname(inim)+interval+'hd'+os.path.basename(inim)[:-5]+'.psf'
		subcat = os.path.dirname(inim)+interval+'hd'+os.path.basename(inim)[:-5]+'.cat'

		hdr = fits.getheader(subim)

		#	PHOTOMETRY NAMES
		param = path_config+'/gregoryphot.param'
		conv = path_config+'/gregoryphot.conv'
		nnw = path_config+'/gregoryphot.nnw'
		conf = path_config+'/gregoryphot.sex'

		seg = head+'.seg.fits'
		bkg = head+'.bkg.fits'
		sub = head+'.sub.fits'
		psf = head+'.psf'
		aper = head+'.aper.fits'

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
							DETECT_MINAREA = '3',
							DETECT_THRESH = '1.0',
							DEBLEND_NTHRESH = '64',
							DEBLEND_MINCONT = '0.0000001',
							#------------------------------
							#	PHOTOMETRY
							#------------------------------
							PHOT_APERTURES = str(hdr['APERPIX']),
							SATUR_LEVEL  = '65000.0',
							MAG_ZEROPOINT = '0.0',
							GAIN = str(gain),
							PIXEL_SCALE = str(pixscale),
							#------------------------------
							#	STAR/GALAXY SEPARATION
							#------------------------------
							SEEING_FWHM = str(hdr['SEEING']),
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

		zp = hdr['ZP']
		zper = hdr['ZPERR']
		seeing = hdr['SEEING']
		aperture = hdr['APERPIX']
		date_obs = hdr['DATE-OBS']
		jd = Time(date_obs, format='isot').jd
		skymed = hdr['SKYVAL']
		skysig = hdr['SKYSIG']
		ul_3sig = hdr['UL3SIG']
		ul_5sig = hdr['UL5SIG']

		os.system(phot.sexcom(subim, param_subsex))
		subtbl = ascii.read(subcat)
		subtbl['REAL_'+inmagkey] = zp + subtbl[inmagkey]
		subtbl['REAL_'+inmagerkey] = phot.sqsum(zper, subtbl[inmagerkey])
		indx_targ = phot.targetfind(tra, tdec, subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'], sep=2*hdr['SEEING'])
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
	onetbl	= Table([[inim], [obs], [obj], [radeg], [dedeg], [date_obs], [jd], [refmagkey], [int(hdr['STDNUMB'])], [round(zp, 3)], [round(zper, 3)], [round(seeing, 3)], [round(skymed, 3)], [round(skysig, 3)], [round(ul_3sig, 3)], [round(ul_5sig, 3)], [round(mag, 3)], [round(mager, 3)], [round(aperture, 3)]],
					names=('image', 'obs', 'obj', 'ra', 'dec', 'date-obs', 'jd', 'filter', 'stdnumb', 'zp', 'zper', 'seeing', 'skyval', 'skysig', 'ul_3sig', 'ul_5sig','mag', 'magerr', 'aper_dia_pix'))

	return onetbl
#============================================================
#	USER SETTING
#============================================================
starttime	= time.time()
os.system('ls *.fits')
imlist		= glob.glob(input('image to process\t: '))
imlist.sort()
for img in imlist: print(img)
print(len(imlist))
#------------------------------------------------------------
# inim = '/data1/Test/Calib-DOAO-NGC5350-20190503-145856-R-60.fits'
path_base = './'
radius = 0.5								#	[DEGREE]
frac = 0.99									#	IMAGE FRACTION [%]
refcatname = 'PS1'							#	REFERENCE CATALOG
# refcatname = 'APASS'
# refcatname = 'SDSS'
# refcatname = '2MASS'
#------------------------------------------------------------
# tra, tdec = 208.3714550, 40.2754194
# tra, tdec = 44.8599883, 31.3863268			#	G0037111
# tra, tdec = 45.4286739, 31.8941469			#	G0232794
# tra, tdec = 75.9126992, -23.7941181	#	UNKNOWN FOR GW190425
# tra, tdec = 185.7288542, 15.8236250			#	SN2020oi
# tra, tdec = 161.6379008, 13.7418711			#	SN2018kp
# tra, tdec = 44.5438520, -8.9577875			#	GRB 190829A
# tra, tdec = 111.8296542, 80.2328528		#	2020ddy
# tra, tdec, length = 185.7288542, 15.8236250, 5	#	SN 2020oi
# tra, tdec, length = 185.4603292, 4.481705551, 5		#	2020jfo
# tra, tdec, length = 185.7338750, 15.8260000, 5
# tra, tdec, length = 185.4603292, 4.4817056, 5	#	2020jfo
# tra, tdec, length = 29.799542, 18.981944, 10 # AT2020uex
tra, tdec, length = 21.0286875, 12.92148055, 10  # AT2020uxz

#------------------------------------------------------------
inmagkey = 'MAG_APER'
inmagerkey = 'MAGERR_APER'
refmaglower, refmagupper = 14, 17			#	REF MAG RANGE [MAG]
refmagerupper = 0.05
inmagerupper = 0.15

# phottype = 'normal'
phottype = 'subt'
# phottype = 'depth'
apertures = np.arange(2, 50.0, 1.5)
# apertures = np.arange(1, 17.0, 0.5)			#	LOAO
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
failist = []
plt.rcParams["figure.figsize"] = (10, 10)
for inim in imlist:
	param_phot = dict(	inim = inim,
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
		failist.append(inim)
		pass
	
	# os.system('rm *.aper.fits *.bkg.fits *.sub.fits')
	# os.system('rm psf*.fits snap*.fits *.xml *.aper.fits *.bkg.fits *.seg.fits')
#------------------------------------------------------------
#	FINISH
#------------------------------------------------------------
if len(tblist) == 0:
	print('PHOTOMETRY FAILED!')
else:
	photbl		= vstack(tblist)
	# if 'phot.dat' in glob.glob(path_base+'/phot.dat'):
		# os.system('mv {} {}'.format(path_base+'/phot.dat', path_base+'/phot.dat.bkg'))
	#	phot.dat REPLACE
	photlist = glob.glob('phot*.dat')

	# if 'phot.dat' in photlist:
	# 	photnumb = 0
	# 	phot_rpl = 'phot.{}.dat'.format(photnumb)
	# 	while phot_rpl in glob.glob('phot*.dat'):
	# 		photnumb += 1
	
	# 	com = 'mv phot.dat {}'.format(phot_rpl)
	# 	print(com)
	# 	os.system(com)

	photbl.write(path_base+'/phot.dat', format='ascii', overwrite=True)

	print('All PROCESS IS DONE.\t('+str(round(time.time() - starttime, 2))+' sec)')

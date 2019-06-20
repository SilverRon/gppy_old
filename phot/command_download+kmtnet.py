#	PHOTOMETRY CODE FOR PYTHON 3.X
#	CREATED	2019.03.09	Gregory S.H. Paek
#	UPDATE	2019.06.01	Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
#from imsng import zpcal
#------------------------------------------------------------
def imcopy(inim, ranges):
	outname	= 'tr'+inim
	chinim	= inim+ranges
	iraf.imcopy(chinim,output=outname)
#============================================================
#	USER SETTING
#============================================================
obsinfo			= ascii.read('/home/sonic/Research/table/obs.txt')
path_catalog	= '/mnt/window/Users/User/Downloads/data/catalog'
path_kmtnet		= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/KMTNet/try1'
path_qso		= '/data3/jkim/GW190425z'
kmtbl			= ascii.read(path_kmtnet+'/kmtn_radec.dat')
downcom			= 'sshpass -prjseka23! scp -ro StrictHostKeyChecking=no paek@qso.snu.ac.kr:'
ra1, de1		= 0, 0
#	REF. CATALOG
refcatname		= 'PS1'					#	PS1/SDSS/APASS/2MASS
detectsig		= 3					#	DETECTION SIGMA
x, y			= 9216, 9232
x0, y0			= x/2, y/2
ranges			= '[{}:{},{}:{}]'.format(int(x0-x/4), int(x0+x/4), int(y0-y/4), int(y0+y/4))
#============================================================
#	MAIN COMMAND
#============================================================
imfail	= []
tblist	= []
for image in kmtbl['image']:
	os.system(downcom+path_qso+'/'+image+' ./')
	query_checklist = glob.glob(path_catalog+'/*.cat')
	os.system('imcopy {}{} {}'.format(image, ranges, 'tr'+image))
	inim	= 'tr'+image

	try:
		tmptbl		= Table()
		hdr			= fits.getheader(inim)
		#------------------------------------------------------------
		try:
			ra, dec		= hdr['ra'], hdr['dec']
			c			= SkyCoord(str(ra)+' '+str(dec), unit=(u.deg, u.deg))
			radeg, dedeg= c.to_string().split(' ')[0], c.to_string().split(' ')[1]
			rahms, dedms= c.to_string('hmsdms').split(' ')[0], c.to_string('hmsdms').split(' ')[1]
		except:
			radeg, dedeg, rahms, dedms	= None, None, None, None
		#------------------------------------------------------------
		date_obs	= hdr['date-obs']
		jd			= round(Time(date_obs, format='isot', scale='utc').jd, 3)
		#------------------------------------------------------------
		part		= inim.split('-')
		obs, obj	= 'KMTNET', inim[2:]
		refmagkey	= 'R'
		refmagerkey = refmagkey+'err'
		gain, pixscale	= obsinfo[obsinfo['obs']==obs]['gain'][0], obsinfo[obsinfo['obs']==obs]['pixelscale'][0]
		#------------------------------------------------------------
		tmptbl['obs'], tmptbl['object'], tmptbl['ra'], tmptbl['dec'], tmptbl['radeg'], tmptbl['decdeg'], tmptbl['date-obs'], tmptbl['jd'], tmptbl['filter'] = \
		[obs], [obj], [rahms], [dedms], [radeg], [dedeg], [date_obs], [jd], [refmagkey]
		#------------------------------------------------------------
		#	SourceEXtractor
		#------------------------------------------------------------
		intbl0, incat, fwhm_pix, fwhm_arcsec	= secom(inim, gain=gain, pixscale=pixscale, det_sigma=detectsig, backsize=str(64))
		#	APPROXIMATE CENTER POS. & DIST CUT
		xim_cent, yim_cent	= np.max(intbl0['X_IMAGE'])/2, np.max(intbl0['Y_IMAGE'])/2
		im_dist		= sqsum((xim_cent-intbl0['X_IMAGE']), (yim_cent-intbl0['Y_IMAGE']))
		indx_dist	= np.where( im_dist < 0.99*(xim_cent+yim_cent)/2. )	# 90% area
		intbl		= intbl0[indx_dist]
		intbl.write(incat, format='ascii', overwrite=True)
		#	NEAR CENTER RA DEC
		radeg       = np.median(intbl['ALPHA_J2000'])
		dedeg       = np.median(intbl['DELTA_J2000'])
		#------------------------------------------------------------
		#	REF. CATALOG QUERY
		#------------------------------------------------------------
		if		refcatname	== 'PS1':
			try:
				if path_catalog+'ps1-'+obj+'.cat' not in query_checklist:
					querytbl        = ps1_query(obj, radeg, dedeg, path_catalog, radius=0.65)
				else:
					querytbl        = ascii.read(path_catalog+'ps1-'+obj+'.cat')
				reftbl, refcat  = ps1_Tonry(querytbl, obj)
			except:
				try:
					querytbl        = sdss_query(obj, radeg, dedeg, path_catalog)
					reftbl, refcat  = sdss_Blaton(querytbl, obj)
				except:
					querytbl        = apass_query(obj, radeg, dedeg, path_catalog)
					reftbl, refcat  = apass_Blaton(querytbl, obj)
	#------------------------------------------------------------
	#	MATCHING
	#------------------------------------------------------------
		mtbl		= matching(incat, './'+refcat)
		colnames    = mtbl.colnames
		maglist     = []
		magerlist   = []
		'''
		#	FOR VARIOUS APERTURE
		for col in colnames:
		if 'MAG_' in col:
			maglist.append(col)
		elif 'MAGERR_' in col:
			magerlist.append(col)
		'''
		for col in colnames:
			if 'MAG_APER_7' in col:
				maglist.append(col)
			elif 'MAGERR_APER_7' in col:
				magerlist.append(col)
		for i in range(0, len(maglist)):
			inmagkey    = maglist[i]
			inmagerkey  = magerlist[i]
			param_st4zp	= dict(	intbl=mtbl,
								inmagerkey=inmagerkey,
								refmagkey=refmagkey,
								refmagerkey=refmagerkey,
								refmaglower=13,
								refmagupper=17,
								refmagerupper=0.05,
								inmagerupper=0.05,
								class_star_cut=0.001)
			stars_zp	= star4zp(**param_st4zp)
			#------------------------------------------------------------
			zp, zper, intbl_alive, intbl_exile	= zpcal(stars_zp, inmagkey, inmagerkey, refmagkey, refmagerkey)
			zp_plot(inim, inmagkey, zp, zper, intbl_alive[inmagkey], intbl_alive[refmagkey], intbl_alive['zp'], intbl_exile[inmagkey], intbl_exile[refmagkey], intbl_exile['zp'])
			intbl['REAL_'+inmagkey]		= zp + intbl[inmagkey]
			intbl['REAL_'+inmagerkey]	= sqsum(zper, intbl[inmagerkey])
	#------------------------------------------------------------
	#	TARGET PHOT
	#------------------------------------------------------------
		ra2, de2	= intbl['ALPHA_J2000'], intbl['DELTA_J2000']
		indx_target		= targetfind(ra1, de1, ra2, de2, sep=10)
		skymean, skymed, skysig		= bkgest_mask(inim)
		#------------------------------------------------------------
		aper			= 2*fwhm_pix
		ul				= limitmag(detectsig, zp, aper, skysig)
		seeing, peeing	= round(fwhm_arcsec, 3), round(fwhm_pix, 3)
		#------------------------------------------------------------
		if len(indx_target[0])!=0:
			mag, magerr		= round(intbl[indx_target]['REAL_MAG_APER_7'][0], 3), round(intbl[indx_target]['REAL_MAGERR_APER_7'][0], 3)
		else:
			mag, magerr		= None, None
		#------------------------------------------------------------
		#	PLOT IMAGE
		#------------------------------------------------------------
		param_plot		= dict(	inim		= inim,
								numb_list	= intbl_alive['NUMBER'],
								xim_list	= intbl_alive['X_IMAGE'],
								yim_list	= intbl_alive['Y_IMAGE'],
								add			= True,
								numb_addlist= intbl_exile['NUMBER'],
								xim_addlist	= intbl_exile['X_IMAGE'],
								yim_addlist	= intbl_exile['Y_IMAGE'])
		plotshow(**param_plot)
		#------------------------------------------------------------
		tmptbl['ul'], tmptbl['zp'], tmptbl['zperr'], tmptbl['mag'], tmptbl['magerr'], tmptbl['seeing'], tmptbl['skyval'], tmptbl['skysig'], tmptbl['stdnumb']	=\
		ul, round(zp, 3), round(zper, 3), mag, magerr, seeing, round(skymed, 3), round(skysig, 3), len(intbl_alive)
		tblist.append(tmptbl)
	except:
		imfail.append(inim)
		pass
	os.system('rm *.fits *.xml *.psf *.cat')
#-------------------------------------------------------------------------#
#	PHOTOMETRY RESULT
#-------------------------------------------------------------------------#
photbl	= vstack(tblist)
phot_checklist = glob.glob('phot.dat')
if 'phot.dat' in phot_checklist:
	os.system('mv phot.dat phot.dat.bkg')
	photbl.write('phot.dat', format='ascii', overwrite=True)
else:
	photbl.write('phot.dat', format='ascii', overwrite=True)
#-------------------------------------------------------------------------#
#	END & CLEAN
#-------------------------------------------------------------------------#
os.system('mkdir zpcal/;mv ./*zpcal.png ./zpcal/')
os.system('rm *aper.fits *xml snap*.fits psf-*.fits')
os.system('mkdir overview/;mv ./*png ./overview/')
comment		= '='*60;print(comment)
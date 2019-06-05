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
#from imsng import zpcal
#============================================================
#	USER SETTING
#============================================================
obsinfo			= ascii.read('/home/sonic/Research/table/obs.txt')
path_catalog	= '/mnt/window/Users/User/Downloads/data/catalog'
#------------------------------------------------------------
#	INITIAL INPUT
#------------------------------------------------------------
#	TARGET COORD.
#ra1, de1	= 54.50492875, -26.94636444		#	GRB 190114C
#ra1, de1	= 173.137, +27.699 				#	GRB 130427A
#ra1, de1	= 196.5942029, +20.35490083
#ra1, de1	= 223.3201092, 34.75006139		#	AT2019ein
#ra1, de1	= 185.733875, 15.826			#	SN2019ehk
ra1, de1	= 161.63775, 13.74194444
#------------------------------------------------------------
#	IMAGES TO CALC.
os.system('ls *.fits')
imlist		= glob.glob(input('image to process\t: '))
imlist.sort()
for img in imlist: print(img)
#	REF. CATALOG
refcatname	= 'APASS'					#PS1/SDSS/APASS/2MASS
detectsig	= 3						#	DETECTION SIGMA
#============================================================
#	MAIN COMMAND
#============================================================
imfail	= []
tblist	= []
for inim in imlist:
	query_checklist = glob.glob(path_catalog+'/*.cat')
	try:
		tmptbl		= Table()
		hdr			= fits.getheader(inim)
		#------------------------------------------------------------
		ra, dec		= hdr['ra'], hdr['dec']
		c			= SkyCoord(str(ra)+' '+str(dec), unit=(u.deg, u.deg))
		radeg, dedeg= c.to_string().split(' ')[0], c.to_string().split(' ')[1]
		rahms, dedms= c.to_string('hmsdms').split(' ')[0], c.to_string('hmsdms').split(' ')[1]
		#------------------------------------------------------------
		date_obs	= hdr['date-obs']
		jd			= round(Time(date_obs, format='isot', scale='utc').jd, 3)
		#------------------------------------------------------------
		part		= inim.split('-')
		obs, obj	= part[1], part[2]
		refmagkey	= part[5]
		refmagerkey = refmagkey+'err'
		gain, pixscale	= obsinfo[obsinfo['obs']==obs]['gain'][0], obsinfo[obsinfo['obs']==obs]['pixelscale'][0]
		#------------------------------------------------------------
		tmptbl['obs'], tmptbl['object'], tmptbl['ra'], tmptbl['dec'], tmptbl['radeg'], tmptbl['decdeg'], tmptbl['date-obs'], tmptbl['jd'], tmptbl['filter'] = \
		[obs], [obj], [rahms], [dedms], [radeg], [dedeg], [date_obs], [jd], [refmagkey]
		#------------------------------------------------------------
		#	SourceEXtractor
		#------------------------------------------------------------
		intbl0, incat, fwhm_pix, fwhm_arcsec	= secom(inim, gain=gain, pixscale=pixscale, det_sigma=3.0, backsize=str(64))
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
			if path_catalog+'ps1-'+obj+'.cat' not in query_checklist:
				querytbl        = ps1_query(obj, radeg, dedeg, path_catalog, radius=0.65)
			else:
				querytbl        = ascii.read(path_catalog+'ps1-'+obj+'.cat')
			reftbl, refcat  = ps1_Tonry(querytbl, obj)
		elif 	refcatname	== 'SDSS':
			if path_catalog+'sdss-'+obj+'.cat' not in query_checklist:
				querytbl        = sdss_query(obj, radeg, dedeg, path_catalog)
			else:
				querytbl        = ascii.read(path_catalog+'sdss-'+obj+'.cat')
			reftbl, refcat  = sdss_Blaton(querytbl, obj)

		elif	refcatname	== 'APASS':
			if path_catalog+'apass-'+obj+'.cat' not in query_checklist:
				querytbl        = apass_query(obj, radeg, dedeg, path_catalog)
			else:
				querytbl        = ascii.read(path_catalog+'apass-'+obj+'.cat')
			reftbl, refcat  = apass_Blaton(querytbl, obj)
		elif	refcatname	== '2MASS':
			if path_catalog+'2mass-'+obj+'.cat' not in query_checklist:
				querytbl        = twomass_query(obj, radeg, dedeg, path_catalog, band=refmagkey, radius=1.0)
			else:
				querytbl        = ascii.read(path_catalog+'2mass-'+obj+'.cat')
			reftbl, refcat  = querytbl, '2mass-'+obj+'.cat'
	#------------------------------------------------------------
	#	MATCHING
	#------------------------------------------------------------
		mtbl		= matching(incat, path_catalog+'/'+refcat)
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
								refmagupper=20,
								refmagerupper=0.05,
								inmagerupper=0.1,
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
		#------------------------------------------------------------
		#	PUT INFO IN HEADER
		#------------------------------------------------------------
		puthdr(inim, 'SEEING',		seeing,					hdrcomment='SEEING [arcsec]')
		puthdr(inim, 'PEEING',		peeing,					hdrcomment='SEEING [pixel]')
		puthdr(inim, 'SKYSIG',		round(skysig, 3),		hdrcomment='SKY SIGMA VALUE')
		puthdr(inim, 'SKYVAL',		round(skymed, 3),		hdrcomment='SKY MEDIAN VALUE')
		puthdr(inim, 'OPTZP',		round(zp, 3),			hdrcomment='2*SEEING DIAMETER')
		puthdr(inim, 'OPTZPERR',	round(zper, 3),			hdrcomment='2*SEEING DIAMETER')
		puthdr(inim, 'STDNUMB',		len(intbl_alive),		hdrcomment='# OF STD STARS FROM {}'.format(refcatname))
		puthdr(inim, 'LIMTMAG',		round(ul, 3),			hdrcomment='{} SIGMA DETECTION'.format(detectsig))
	except:
		imfail.append(inim)
		pass
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
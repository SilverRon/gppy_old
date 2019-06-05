#	PHOTOMETRY CODE (TEST) FOR PYTHON 3.X
#	2019.03.09
#	GREGORY S.H. PAEK
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
#from imsng import zpcal
#============================================================
#	USER SETTING
#============================================================
sharepath		= '/home/sonic/Research/yourpy/config/'

configfile      = sharepath+'targetphot.sex'
paramfile       = sharepath+'targetphot.param'
nnwfile		    = sharepath+'targetphot.nnw'
convfile	    = sharepath+'targetphot.conv'

psfexconf_prese_conf    = sharepath+'prepsfex.sex'
psfexconf_prese_param   = sharepath+'prepsfex.param'
psfexconf_psfex_conf    = sharepath+'default.psfex'
psfexconf_psfex_conv    = sharepath+'default.conv'

obsinfo			= ascii.read('/home/sonic/Research/table/obs.txt')
#------------------------------------------------------------
def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)	
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'
	#print(comment)
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

#	IMAGES TO CALC.
#imlist		= glob.glob('Calib-*com.fits')
os.system('ls *.fits')
imlist		= glob.glob(input('image to process\t: '))
imlist.sort()
for img in imlist: print(img)
#	REF. CATALOG
refcatname	= 'PS1'					#PS1/SDSS/APASS/2MASS
#	RESULT FILE
f		= open('phot.dat', 'w')
colline	= '#obs\tdate-obs\taperture\tseeing\tzp\tzperr\tinstmag\tinstmagerr\tmag\tmagerr\n'
f.write(colline)
#============================================================
#	MAIN COMMAND
#============================================================
imfail	= []
for inim in imlist:
	query_checklist = glob.glob('*.cat')
	try:
		hdr			= fits.getheader(inim)
		part        = inim.split('-')
		obs         = part[1]
		name        = part[2]
		exptime		= part[6]
		refmagkey   = part[5]
		refmagerkey = refmagkey+'err'

		gain		= obsinfo[obsinfo['obs']==obs]['gain'][0]
		pixscale	= obsinfo[obsinfo['obs']==obs]['pixelscale'][0]
		#	SourceEXtractor
		intbl0, incat, fwhm_pix, fwhm_arcsec		= secom(inim, gain=gain, pixscale=pixscale, det_sigma=3.0, backsize=str(64))
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
			if 'ps1-'+name+'.cat' not in query_checklist:
				querytbl        = ps1_query(name, radeg, dedeg, radius=0.65)
			else:
				querytbl        = ascii.read('ps1-'+name+'.cat')
			reftbl, refcat  = ps1_Tonry(querytbl, name)

		elif 	refcatname	== 'SDSS':
			if 'sdss-'+name+'.cat' not in query_checklist:
				querytbl        = sdss_query(name, radeg, dedeg)
			else:
				querytbl        = ascii.read('sdss-'+name+'.cat')
			reftbl, refcat  = sdss_Blaton(querytbl, name)

		elif	refcatname	== 'APASS':
			if 'apass-'+name+'.cat' not in query_checklist:
				querytbl        = apass_query(name, radeg, dedeg)
			else:
				querytbl        = ascii.read('apass-'+name+'.cat')
			reftbl, refcat  = apass_Blaton(querytbl, name)
		elif	refcatname	== '2MASS':
			if '2mass-'+name+'.cat' not in query_checklist:
				querytbl        = twomass_query(name, radeg, dedeg, band=refmagkey, radius=1.0)
			else:
				querytbl        = ascii.read('2mass-'+name+'.cat')
			reftbl, refcat  = querytbl, '2mass-'+name+'.cat'
	#------------------------------------------------------------
	#	MATCHING
	#------------------------------------------------------------
		merge_raw	= matching(incat, refcat)
		colnames    = merge_raw.colnames
		maglist     = []
		magerlist   = []
		for col in colnames:
			if 'MAG_APER_7' in col:
				#print(col)
				maglist.append(col)
			elif 'MAGERR_APER_7' in col:
				#print(col)
				magerlist.append(col)
		#intbl	= ascii.read(incat)
		for i in range(0, len(maglist)):
			mtbl       = merge_raw
			inmagkey    = maglist[i]
			inmagerkey  = magerlist[i]

			param_st4zp	= dict(	intbl=mtbl,
								inmagerkey=inmagerkey,
								refmagkey=refmagkey,
								refmagerkey=refmagerkey,
								refmaglower=13,
								refmagupper=16.5,
								refmagerupper=0.05,
								inmagerupper=0.1,
								class_star_cut=0.001)

			stars_zp	= star4zp(**param_st4zp)
			#stars_zp = star4zp(mtbl, inmagerkey, refmagkey, refmagerkey, refmaglower=14, refmagupper=16.5, refmagerupper=0.05, inmagerupper=0.1, class_star_cut=0.001)
			#stars_zp, stdnumb    = star4zp(mtbl, inmagerkey, refmagkey, refmagerkey, refmaglower=14, refmagupper=18, refmagerupper=0.05, inmagerupper=0.1, class_star_cut=0.01)

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

		if len(indx_target[0]) == 0:
			aper	= 2*fwhm_pix
			ul		= limitmag(3, zp, aper, skysig)
			try:
				comment		= inim+'\t\t'+hdr['date-obs']+'\t\t'+'MAG_APER_7'+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
							+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
							+'\t--\t\t\t'+'--\t' \
							+'\t'+str(round(ul, 3))+'\t'+'0'+'\n'
			except:
				comment		= inim+'\t\t'+'MAG_APER_7'+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
							+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
							+'\t--\t\t\t'+'--\t' \
							+'\t'+str(round(ul, 3))+'\t'+'0'+'\n'
			print(comment)
			f.write(comment)

		else:
			try:
				comment		= inim+'\t\t'+hdr['date-obs']+'\t\t'+'MAG_APER_7'+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
							+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
							+'\t'+str(round(intbl[indx_target]['MAG_APER_7'][0], 3))+'\t\t'+str(round(intbl[indx_target]['MAGERR_APER_7'][0], 3)) \
							+'\t'+str(round(intbl[indx_target]['REAL_MAG_APER_7'][0], 3))+'\t'+str(round(intbl[indx_target]['REAL_MAGERR_APER_7'][0], 3))+'\n'
			except:
				comment		= inim+'\t\t'+'MAG_APER_7'+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
							+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
							+'\t'+str(round(intbl[indx_target]['MAG_APER_7'][0], 3))+'\t\t'+str(round(intbl[indx_target]['MAGERR_APER_7'][0], 3)) \
							+'\t'+str(round(intbl[indx_target]['REAL_MAG_APER_7'][0], 3))+'\t'+str(round(intbl[indx_target]['REAL_MAGERR_APER_7'][0], 3))+'\n'

			print(comment)
			f.write(comment)
		#	PLOT IMAGE
		numb_list	= intbl_alive['NUMBER']
		xim_list	= intbl_alive['X_IMAGE']
		yim_list	= intbl_alive['Y_IMAGE']
		numb_addlist= intbl_exile['NUMBER']
		xim_addlist	= intbl_exile['X_IMAGE']
		yim_addlist	= intbl_exile['Y_IMAGE']
		plotshow(inim, numb_list, xim_list, yim_list, add=True, numb_addlist=numb_addlist, xim_addlist=xim_addlist, yim_addlist=yim_addlist)

		puthdr(inim, 'SEEING',		round(fwhm_arcsec, 3),	hdrcomment='SEEING [arcsec]')
		puthdr(inim, 'PEEING',		round(fwhm_pix, 3),		hdrcomment='SEEING [pixel]')
		puthdr(inim, 'STDNUMB',		len(intbl_alive),		hdrcomment='# OF STD STARS')
		puthdr(inim, 'OPTZP',		round(zp, 3),			hdrcomment='2*SEEING DIAMETER')
		puthdr(inim, 'OPTZPERR',	round(zper, 3),			hdrcomment='2*SEEING DIAMETER')
		puthdr(inim, 'SKYSIG',		round(skysig, 3),		hdrcomment='SKY SIGMA VALUE')
		puthdr(inim, 'SKYVAL',		round(skymed, 3),		hdrcomment='SKY MEDIAN VALUE')
	except:
		imfail.append(inim)
		pass
#-------------------------------------------------------------------------#
f.close()

photbl	= ascii.read('phot.dat')
#photbl[photbl['mag']>20]

comment		= '='*60;print(comment)
os.system('mkdir zpcal/;mv ./*zpcal.png ./zpcal/')
os.system('mkdir zpcal_test/;mv ./*zpcal_test.png ./zpcal_test/')
os.system('rm *aper.fits *xml snap*.fits psf-*.fits')
os.system('mkdir overview/;mv ./*png ./overview/')
os.system('cat phot.dat')

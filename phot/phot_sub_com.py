#	PHOTOMETRY CODE (TEST) FOR PYTHON 3.X
#	2019.03.09
#	GREGORY S.H. PAEK
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table
from astropy.time import Time

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
photbl			= Table()
obslist			= []
bandlist		= []
datelist		= []
photlist		= []
photerlist		= []
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
#ra1, de1	= 54.510, -26.947		#	GRB 190114C
#ra1, de1	= 178.2218750, 44.1239444 		#	SN 2017ein
#ra1, de1	= 94.09275, -21.35991111
#ra1, de1	= 248.5427121, 19.634815
#ra1, de1	= 208.3721617, 40.27532028		#	AT 2019ein
ra1, de1	= 185.733875, 15.826			#	SN2019ehk

#	IMAGES TO CALC.
#imlist		= glob.glob('Calib-*com.fits')
os.system('ls *.fits')
imlist		= glob.glob(input('IMAGE TO PROCESS\t: '))


imlist.sort()
for img in imlist: print(img)
#	REF. CATALOG
refcatname	= 'PS1'					#PS1/SDSS/APASS/2MASS



#	RESULT FILE
f		= open('phot.dat', 'w')
colline	= '#obs\tMJD\taperture\tseeing\tzp\tzperr\tinstmag\tinstmagerr\tmag\tmagerr\n'
f.write(colline)
#============================================================
#	MAIN COMMAND
#============================================================
n	= 1
for inim in imlist:
	print('['+str(n)+'/'+str(len(imlist))+']'); n += 1
	query_checklist = glob.glob('*.cat')
	part        = inim.split('-')
	obs         = part[1]
	name        = part[2]
	exptime		= part[6]
	refmagkey   = part[5]
	refmagerkey = refmagkey+'err'

	gain		= obsinfo[obsinfo['obs']==obs]['gain'][0]
	pixscale	= obsinfo[obsinfo['obs']==obs]['pixelscale'][0]
	#	SourceEXtractor
	intbl0, incat, fwhm_pix, fwhm_arcsec		= secom(inim, gain=gain, pixscale=pixscale, det_sigma=3, backsize=str(64))
	#	APPROXIMATE CENTER POS. & DIST CUT
	xim_cent, yim_cent	= np.max(intbl0['X_IMAGE'])/2, np.max(intbl0['Y_IMAGE'])/2
	im_dist		= sqsum((xim_cent-intbl0['X_IMAGE']), (yim_cent-intbl0['Y_IMAGE']))
	indx_dist	= np.where( im_dist < 0.95*(xim_cent+yim_cent)/2. )	# 90% area
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
			querytbl        = ps1_query(name, radeg, dedeg)
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
		#reftbl, refcat  = apass_Blaton(querytbl, name)
		reftbl, refcat	= apass_Lupton(querytbl, name)
#------------------------------------------------------------
#	MATCHING
#------------------------------------------------------------
	merge_raw       = matching(incat, refcat)

	colnames    = merge_raw.colnames
	maglist     = []
	magerlist   = []
	apname		= 'MAG_APER_7'
	apername	= 'MAGERR_APER_7'
	for col in colnames:
		if apname == col:
			#print(col)
			maglist.append(col)
		elif apername == col:
			#print(col)
			magerlist.append(col)
		# elif 'MAG_PSF' in col:
		# 	maglist.append(col)
		# elif 'MAGERR_PSF' in col:
		# 	magerlist.append(col)

	#intbl	= ascii.read(incat)
	# zp of input image

	skymean, skymed, skysig		= bkgest_mask(inim)
	puthdr(inim, 'SEEING',		round(fwhm_arcsec, 3),	hdrcomment='SEEING [arcsec]')
	puthdr(inim, 'PEEING',		round(fwhm_pix, 3),		hdrcomment='SEEING [pixel]')
	puthdr(inim, 'SKYSIG',		'%.3f'% round(skysig, 3),		hdrcomment='SKY SIGMA VALUE')
	puthdr(inim, 'SKYVAL',		'%.3f'% round(skymed, 3),		hdrcomment='SKY MEDIAN VALUE')

	for i in range(0, len(maglist)):
		mtbl       = merge_raw
		inmagkey    = maglist[i]
		inmagerkey  = magerlist[i]
		
		stars_zp = star4zp(mtbl, inmagerkey, refmagkey, refmagerkey, refmaglower=13, refmagupper=18, refmagerupper=0.05, inmagerupper=0.1, class_star_cut=0.01)
		#stars_zp, stdnumb    = star4zp(mtbl, inmagerkey, refmagkey, refmagerkey, refmaglower=14, refmagupper=18, refmagerupper=0.05, inmagerupper=0.1, class_star_cut=0.01)

		zp, zper, intbl_alive, intbl_exile	= zpcal(stars_zp, inmagkey, inmagerkey, refmagkey, refmagerkey)
		aper	= 2*fwhm_pix
		ul		= limitmag(3, zp, aper, skysig)

		zp_plot(inim, inmagkey, zp, zper, intbl_alive[inmagkey], intbl_alive[refmagkey], intbl_alive['zp'], intbl_exile[inmagkey], intbl_exile[refmagkey], intbl_exile['zp'])
		intbl['REAL_'+inmagkey]		= zp + intbl[inmagkey]
		intbl['REAL_'+inmagerkey]	= sqsum(zper, intbl[inmagerkey])

		
		puthdr(inim, 'OPTZP',		round(zp, 3),			hdrcomment='2*SEEING DIAMETER')
		puthdr(inim, 'OPTZPERR',	round(zper, 3),			hdrcomment='2*SEEING DIAMETER')
		puthdr(inim, 'OPTUL',		round(ul, 3),			hdrcomment='2*SEEING 3 sigma limit mag')
	#puthdr(inim, 'PSFZP',		round(zp, 3),			hdrcomment='2*SEEING DIAMETER')
	#puthdr(inim, 'PSFZPERR',	round(zper, 3),			hdrcomment='2*SEEING DIAMETER')

		puthdr(inim, 'STDNUMB',		len(intbl_alive),		hdrcomment='# OF STD STARS')
			#	PLOT IMAGE
		numb_list	= intbl_alive['NUMBER']
		xim_list	= intbl_alive['X_IMAGE']
		yim_list	= intbl_alive['Y_IMAGE']
		numb_addlist= intbl_exile['NUMBER']
		xim_addlist	= intbl_exile['X_IMAGE']
		yim_addlist	= intbl_exile['Y_IMAGE']
		plotshow(inim, numb_list, xim_list, yim_list, add=True, numb_addlist=numb_addlist, xim_addlist=xim_addlist, yim_addlist=yim_addlist)



	# subtracted image, dual sextracted
	#subim='hd'+inim
	intbl, incat, fwhm_pix, fwhm_arcsec		= sedualcom('hd'+inim, gain, pixscale, det_sigma=1.5, backsize=str(64), backfiltersize=str(3), detect='detection-AT2019eez.fits')
	ra2, de2	= intbl['ALPHA_J2000'], intbl['DELTA_J2000']
	indx_target		= targetfind(ra1, de1, ra2, de2, sep=5)
	intbl['REAL_'+inmagkey]		= zp + intbl[inmagkey]
	intbl['REAL_'+inmagerkey]	= sqsum(zper, intbl[inmagerkey])


	if len(indx_target[0]) == 0:
		aper	= 2*fwhm_pix
		ul		= limitmag(3, zp, aper, skysig)
		comment		= inim+'\t\t'+inmagkey+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
					+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
					+'\t--\t\t\t'+'--\t' \
					+'\t'+str(round(ul, 3))+'\t'+'n/d'+'\n'
		print(comment)
		f.write(comment)

	else:
		'''
		comment		= inim+'\t\t'+ str(round(fits.getheader(inim)['DATE-OBS'],5))+'\t\t'+ 'MAG_APER_7'+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
					+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
					+'\t'+str(round(intbl[indx_target]['MAG_APER_7'][0], 3))+'\t\t'+str(round(intbl[indx_target]['MAGERR_APER_7'][0], 3)) \
					+'\t'+str(round(intbl[indx_target]['REAL_MAG_APER_7'][0], 3))+'\t'+str(round(intbl[indx_target]['REAL_MAGERR_APER_7'][0], 3))+'\n'
		'''
		comment		=  inim+'\t\t'+ str(fits.getheader(inim)['DATE-OBS'])+'\t\t'+ inmagkey+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
					+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
					+'\t'+str(round(intbl[indx_target][inmagkey][0], 3))+'\t\t'+str(round(intbl[indx_target][inmagerkey][0], 3)) \
					+'\t'+str(round(intbl[indx_target]['REAL_'+inmagkey][0], 3))+'\t'+str(round(intbl[indx_target]['REAL_'+inmagerkey][0], 3))+'\n'
		print(comment)
		f.write(comment)
	#-------------------------------------------------------------------------#
	obslist.append(obs)
	bandlist.append(refmagkey)
	datelist.append(fits.getheader(inim)['DATE-OBS'])
	if len(indx_target[0]) != 0:
		photlist.append(round(intbl[indx_target]['REAL_'+inmagkey][0], 3))
		photerlist.append(round(intbl[indx_target]['REAL_'+inmagerkey][0], 3))
	else:
		photlist.append(round(ul, 3))
		photerlist.append(-99)




#-------------------------------------------------------------------------#
f.close()

jdlist	= []
for date_obs in datelist:
	jdlist.append(Time(date_obs, format='isot', scale='utc').jd)
jdlist= np.array(jdlist)

obslist			= np.array(obslist)
bandlist		= np.array(bandlist)
datelist		= np.array(datelist)
jdlist			= np.array(jdlist)
photlist		= np.array(photlist)
photerlist		= np.array(photerlist)



photbl['obs'], photbl['date-obs'], photbl['jd'], photbl['band'], photbl['mag'], photbl['magerr'] = obslist, datelist, jdlist, bandlist, photlist, photerlist

photbl.write('phot_table.dat', format='ascii', overwrite=True)

comment		= '='*60;print(comment)
os.system('mkdir zpcal/;mv ./*zpcal.png ./zpcal/')
os.system('mkdir zpcal_test/;mv ./*zpcal_test.png ./zpcal_test/')
os.system('rm *aper.fits *xml snap*.fits psf-*.fits')
os.system('mkdir overview/;mv ./*png ./overview/')
os.system('cat phot.dat')


"""
# target photometry 
ra2, de2	= intbl['ALPHA_J2000'], intbl['DELTA_J2000']
indx_target		= targetfind(ra1, de1, ra2, de2, sep=3)
	if len(indx_target[0]) == 0:
	aper	= 2*fwhm_pix
	ul		= limitmag(3, zp, aper, skysig)
	comment		= inim+'\t\t'+'MAG_APER_7'+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
				+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
				+'\t--\t\t\t'+'--\t' \
				+'\t'+str(round(ul, 3))+'\t'+'n/d'+'\n'
	print(comment)
	f.write(comment)
	else:
	comment		= inim+'\t\t'+'MAG_APER_7'+'\t\t'+str(round(fwhm_arcsec, 3))+'\t\t' \
				+str(round(zp, 3))+'\t'+str(round(zper, 3)) \
				+'\t'+str(round(intbl[indx_target]['MAG_APER_7'][0], 3))+'\t\t'+str(round(intbl[indx_target]['MAGERR_APER_7'][0], 3)) \
				+'\t'+str(round(intbl[indx_target]['REAL_MAG_APER_7'][0], 3))+'\t'+str(round(intbl[indx_target]['REAL_MAGERR_APER_7'][0], 3))+'\n'
	print(comment)
	f.write(comment)
"""







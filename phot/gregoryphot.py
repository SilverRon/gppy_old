#	PHOTOMETRY CODE FOR PYTHON 3.X
#	CREATED	2019.06.20	Gregory S.H. Paek
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
#from multiprocessing import Process, Pool
#import multiprocessing as mp
from imsng import phot
import time
#============================================================
#	FUNCTION
#============================================================
def phot_routine(inim, refcatname, phottype, tra, tdec, path_base='./', aperture='MAG_APER_7', detsig=3.0, frac=0.95):
	#------------------------------------------------------------
	#	HEADER INFO
	hdr			= fits.getheader(inim)
	try:
		w			= WCS(inim)
		radeg, dedeg= w.all_pix2world(xcent, ycent, 1)
		radeg, dedeg= np.asscalar(radeg), np.asscalar(dedeg)
	except:
		print('BAD WCS INFORMATION?')
		radeg,dedeg	= hdr['CRVAL1'], hdr['CRVAL2']
	#xcent, ycent= w.all_world2pix(radeg, dedeg, 1)
	xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
	#------------------------------------------------------------
	try:
		date_obs	= hdr['date-obs']
		jd			= round(Time(date_obs, format='isot', scale='utc').jd, 3)
	except:
		date_obs	= None
		jd			= None
	#------------------------------------------------------------
	#	NAME INFO
	part			= inim.split('-')
	obs, obj		= part[1], part[2]
	refmagkey		= part[5]
	refmagerkey 	= refmagkey+'err'
	indx_obs		= np.where(obstbl['obs']==obs)
	gain, pixscale	= obstbl[indx_obs]['gain'][0], obstbl[indx_obs]['pixelscale'][0]
	#------------------------------------------------------------
	#	REF. CATALOG QUERY
	#------------------------------------------------------------
	refcatlist	= glob.glob(path_refcat+'/*.cat')
	#------------------------------------------------------------
	if		refcatname	== 'PS1':
		if path_refcat+'/ps1-'+obj+'.cat' not in refcatlist:
			querytbl	= phot.ps1_query(obj, radeg, dedeg, path_refcat, radius=3.0)
		else:
			querytbl	= ascii.read(path_refcat+'/ps1-'+obj+'.cat')
		reftbl, refcat  = phot.ps1_Tonry(querytbl, obj)
	#------------------------------------------------------------
	elif 	refcatname	== 'SDSS':
		if path_refcat+'/sdss-'+obj+'.cat' not in refcatlist:
			querytbl        = phot.sdss_query(obj, radeg, dedeg, path_refcat)
		else:
			querytbl        = ascii.read(path_refcat+'/sdss-'+obj+'.cat')
		reftbl, refcat  = phot.sdss_Blaton(querytbl, obj)
	#------------------------------------------------------------
	elif	refcatname	== 'APASS':
		if path_refcat+'/apass-'+obj+'.cat' not in refcatlist:
			querytbl        = phot.apass_query(obj, radeg, dedeg, path_refcat)
		else:
			querytbl        = ascii.read(path_refcat+'/apass-'+obj+'.cat')
		reftbl, refcat  = phot.apass_Blaton(querytbl, obj)
	#------------------------------------------------------------
	elif	refcatname	== '2MASS':
		if path_refcat+'/2mass-'+obj+'.cat' not in refcatlist:
			querytbl        = phot.twomass_query(obj, radeg, dedeg, path_refcat, band=refmagkey, radius=1.0)
		else:
			querytbl        = ascii.read(path_refcat+'/2mass-'+obj+'.cat')
		reftbl, refcat  = querytbl, '2mass-'+obj+'.cat'
	#------------------------------------------------------------
	#	SourceEXtractor
	#------------------------------------------------------------
	peeing, seeing	= phot.psfex(inim, pixscale)
	param_secom	= dict(	inim=inim,
						gain=gain, pixscale=pixscale, seeing=seeing,
						det_sigma=detsig,
						backsize=str(64), backfiltersize=str(3),
						psf=True, check=False)
	intbl0, incat	= phot.secom(**param_secom)
	#	CENTER POS. & DIST CUT
	deldist		= phot.sqsum((xcent-intbl0['X_IMAGE']), (ycent-intbl0['Y_IMAGE']))
	indx_dist	= np.where(deldist < np.sqrt(frac)*(xcent+ycent)/2.)
	intbl		= intbl0[indx_dist]
	intbl.write(incat, format='ascii', overwrite=True)
	#	MATCHING
	param_match	= dict(	intbl=intbl, reftbl=reftbl,
						inra=intbl['ALPHA_J2000'], indec=intbl['DELTA_J2000'],
						refra=reftbl['ra'], refdec=reftbl['dec'])
	mtbl		= phot.matching(**param_match)
	#------------------------------------------------------------
	#	ZEROPOINT CALCULATION
	#------------------------------------------------------------
	inmagkey	= aperture
	inmagerkey	= '_'.join(['MAGERR', inmagkey.split('_')[1], inmagkey.split('_')[2]]) 
	param_st4zp	= dict(	intbl=mtbl,
						inmagerkey=aperture,
						refmagkey=refmagkey,
						refmagerkey=refmagerkey,
						refmaglower=10,
						refmagupper=19,
						refmagerupper=0.1,
						inmagerupper=0.1)
	param_zpcal	= dict(	intbl=phot.star4zp(**param_st4zp),
						inmagkey=inmagkey, inmagerkey=inmagerkey,
						refmagkey=refmagkey, refmagerkey=refmagerkey,
						sigma=2.0)
	zp, zper, otbl, xtbl	= phot.zpcal(**param_zpcal)

	#------------------------------------------------------------
	#	ZEROPOINT PLOT
	#------------------------------------------------------------
	outname	= path_base+'/{0}.{1}.zpcal.png'.format(inim[:-5], inmagkey)
	phot.zpplot(	outname=outname,
					otbl=otbl,xtbl=xtbl,
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
	skymean, skymed, skysig		= phot.bkgest_mask(inim)
	aper	= 2*peeing
	ul		= phot.limitmag(detsig, zp, aper, skysig)
	#------------------------------------------------------------
	#	ADD HEADER INFO
	#------------------------------------------------------------
	phot.puthdr(inim, 'SEEING',	round(seeing, 3),		hdrcomment='SEEING [arcsec]')
	phot.puthdr(inim, 'PEEING',	round(peeing, 3),		hdrcomment='SEEING [pixel]')
	phot.puthdr(inim, 'SKYSIG',	round(skysig, 3),		hdrcomment='SKY SIGMA VALUE')
	phot.puthdr(inim, 'SKYVAL',	round(skymed, 3),		hdrcomment='SKY MEDIAN VALUE')
	phot.puthdr(inim, 'OPTZP',	round(zp, 3),			hdrcomment='2*SEEING DIAMETER')
	phot.puthdr(inim, 'OPTZPERR',round(zper, 3),		hdrcomment='2*SEEING DIAMETER')
	phot.puthdr(inim, 'OPTUL',	round(ul, 3),			hdrcomment='2*SEEING 3 sigma limit mag')
	phot.puthdr(inim, 'STDNUMB',len(otbl),				hdrcomment='# OF STD STARS')
	#------------------------------------------------------------
	#	NORMAL PHOTOMETRY
	#------------------------------------------------------------
	if phottype == 'normal':
		intbl['REAL_'+inmagkey]		= zp + intbl[inmagkey]
		intbl['REAL_'+inmagerkey]	= phot.sqsum(zper, intbl[inmagerkey])
		indx_targ	= phot.targetfind(tra, tdec, intbl['ALPHA_J2000'], intbl['DELTA_J2000'], sep=seeing)
		if indx_targ != None:
			mag, mager	= intbl[indx_targ]['REAL_'+inmagkey], intbl[indx_targ]['REAL_'+inmagerkey]
		else:
			mag, mager	= -99, -99
	#------------------------------------------------------------
	#	SUBTRACTION PHOTOMETRY
	#------------------------------------------------------------
	elif phottype == 'subt':
		subim	= 'hd'+inim
		phot.puthdr(subim, 'SEEING',	round(seeing, 3),		hdrcomment='SEEING [arcsec]')
		phot.puthdr(subim, 'PEEING',	round(peeing, 3),		hdrcomment='SEEING [pixel]')
		phot.puthdr(subim, 'SKYSIG',	round(skysig, 3),		hdrcomment='SKY SIGMA VALUE')
		phot.puthdr(subim, 'SKYVAL',	round(skymed, 3),		hdrcomment='SKY MEDIAN VALUE')
		phot.puthdr(subim, 'OPTZP',	round(zp, 3),				hdrcomment='2*SEEING DIAMETER')
		phot.puthdr(subim, 'OPTZPERR',round(zper, 3),			hdrcomment='2*SEEING DIAMETER')
		phot.puthdr(subim, 'OPTUL',	round(ul, 3),				hdrcomment='2*SEEING 3 sigma limit mag')
		phot.puthdr(subim, 'STDNUMB',len(otbl),					hdrcomment='# OF STD STARS')
		os.system('cp {} {}'.format(inim[:-5]+'.psf', subim[:-5]+'.psf'))
		param_subcom	= dict(	inim=subim,
								gain=gain, pixscale=pixscale, seeing=seeing,
								det_sigma=3,
								backsize=str(64), backfiltersize=str(3),
								psf=True, check=False)
		subtbl, subcat	= phot.secom(**param_subcom)
		subtbl['REAL_'+inmagkey]	= zp + subtbl[inmagkey]
		subtbl['REAL_'+inmagerkey]	= phot.sqsum(zper, subtbl[inmagerkey])
		indx_targ	= phot.targetfind(tra, tdec, subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'], sep=seeing)
		if indx_targ != None:
			mag, mager	= subtbl[indx_targ]['REAL_'+inmagkey], subtbl[indx_targ]['REAL_'+inmagerkey]
		else:
			mag, mager	= -99, -99
	#------------------------------------------------------------
	#	CALC. DEPTH
	#------------------------------------------------------------
	elif phottype == 'depth':
		mag, mager	= -99, -99

	onetbl	= Table([[inim], [obs], [obj], [round(radeg, 3)], [round(dedeg, 3)], [date_obs], [jd], [refmagkey], [len(otbl)], [round(zp, 3)], [round(zper, 3)], [round(seeing, 3)], [round(skymed, 3)], [round(skysig, 3)], [round(ul, 3)], [mag], [mager]],
					names=('image', 'obs', 'obj', 'ra', 'dec', 'date-obs', 'jd', 'filter', 'stdnumb', 'zp', 'zper', 'seeing', 'skyval', 'skysig', 'ul', 'mag', 'magerr'))
	return onetbl
#============================================================
#	USER SETTING
#============================================================
path_base	= './'
path_obs	= '/home/sonic/Research/table'
path_refcat	= '/home/sonic/Research/cat/refcat'
#------------------------------------------------------------
obstbl		= ascii.read(path_obs+'/obs.txt')
#	TARGET COORD.	[deg]
# tra, tdec = 185.733875, 15.826			#	SN2019ehk
# tra, tdec = 258.3414923, -9.964393723	#	ZTF19aarykkb
tra, tdec = 262.7914654, -8.450713499	#	ZTF19aarzaod
#------------------------------------------------------------
#	IMAGES TO PHOTOMETRY
#	INPUT FORMAT	: Calib-[OBS]-[TARGET]-[DATE]-[TIME]-[BAND]*.fits
#------------------------------------------------------------
os.system('ls *.fits')
imlist		= glob.glob(input('image to process\t: '))
imlist.sort()
for img in imlist: print(img)

photlist	= []
refcatname	= 'PS1'			#	(PS1/APASS/SDSS/2MASS)
# refcatname	= 'SDSS'
# refcatname	= 'APASS'
# refcatname	= '2MASS'
phottype	= 'subt'		#	(normal/subt/depth)
# phottype	= 'depth'
starttime	= time.time()
#============================================================
#	MAIN COMMAND
#============================================================
for inim in imlist:
	try:
		param_phot	= dict(	inim=inim, refcatname=refcatname, phottype=phottype,
							tra=tra, tdec=tdec, path_base='./', aperture='MAG_APER_7',
							detsig=3.0, frac=0.9)
		photlist.append(phot_routine(**param_phot))
		os.system('rm psf*.fits snap*.fits *.xml seg.fits')
	except:
		pass
#------------------------------------------------------------
#	FINISH
#------------------------------------------------------------
if len(photlist) == 0:
	print('PHOTOMETRY FAILED!')
else:
	photbl		= vstack(photlist)
	if 'phot.dat' in glob.glob(path_base+'/phot.dat'):
		os.system('mv {} {}'.format(path_base+'/phot.dat', path_base+'/phot.dat.bkg'))
	photbl.write(path_base+'/phot.dat', format='ascii', overwrite=True)
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')

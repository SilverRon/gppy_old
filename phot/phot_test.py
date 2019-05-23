#   ZEROPOINT CALCULATION
#   19.03.09    MERRY CHRISTMAS!
#   MADE BY GREGORY S.H. PAEK
#=========================================================================#
def which_obs(obs, path_and_file):
	"""
	=====================================================================
	GIVE CCD INFORMATION
	---------------------------------------------------------------------
	INPUT   :   observatory name, and path+table_name(ascii format)
				i.e.) path_and_file='/home/sonic/Research/table/obs.txt'
	OUTPUT  :   GAIN        []
				Pixel Scale ["/pix]
	---------------------------------------------------------------------
	obs         ccd         gain    RDnoise     dark    pixelscale
	SOAO        FLI         1.43    13.60       0.0     0.44540
	LOAO        E2V         2.68    4.84        0.0     0.794
	LOAO_FLI    FLI         3       15.0        0.0     1.28
	DOAO_sophia PIsophia    0.928   7.63        0.0     0.3855
	DOAO        FLI         1.27    14.0        0.005   0.464
	oldBOAO     UNKNOWN     8.9     20.0        0.0     1.7
	BOAO        UNKNOWN     5.0     2.0         0.0     0.2145
	BOAO        UNKNOWN     5.0     5.0         0.0     0.2145    
	SAO         SBIG        1.35    9.0         0.0     0.31
	30inch      UNKNOWN     1.60    5.87        0.0     1.3553
	CCA250      MLI16803    1.46    10.35       0.0     2.06
	MAIDANAK    SNU4kCAM    1.45    4.7         0.0     0.266
	MAIDANAK    UNKNOWN     5.40    9.0         0.0     0.266
	LSGT        SNUCAMII    1.15    6.0         0.2     0.92
	UKIRT       WFCAM       99.0    0.0         0.0     0.4
	=====================================================================
	"""
	import numpy as np
	from astropy.io import ascii
	obsinfo     = ascii.read(path_and_file)
	indx_obs    = np.where( obs == obsinfo['obs'] )
	gain        = obsinfo['gain'][indx_obs]
	pixscale    = obsinfo['pixelscale'][indx_obs]
	return gain, pixscale
#-------------------------------------------------------------------------#
def image_list(imlist):
	"""
	INPUT   :   imlist = glob.glob('Calib-*.fits')
	OUTPUT  :   observatory list
				object list
				fillter list
	"""
	obslist = []
	objlist = []
	fillist = []
	for img in imlist:
		sp  = img.split('-')
		obslist.append(sp[1])
		objlist.append(sp[2])
		fillist.append(sp[5])
	obslist = list(set(obslist))
	objlist = list(set(objlist))
	fillist = list(set(fillist))
	return obslist, objlist, fillist
#-------------------------------------------------------------------------#
def sdss_query(name, radeg, dedeg, radius=1.0):
	"""
	SDSS QUERY
	INPUT   :   NAME, RA [deg], Dec [deg], radius
	OUTPUT  :   QUERY TABLE
				sdss-[NAME].cat
	"""
	from astroquery.vizier import Vizier 
	from astropy.coordinates import Angle
	from astropy.table import Table
	import astropy.units as u
	import astropy.coordinates as coord
	comment = 'NAME'+'\t'+': '+name+'\n' \
			+ 'RA'+'\t'+': '+str(round(radeg, 3))+'\n' \
			+ 'Dec'+'\t'+': '+str(round(dedeg, 3))+'\n' \
			+ 'Radius'+'\t'+': '+str(radius)+' deg'+'\n'*2 \
			+ 'LOADING SDSS Catalog ...'+'\n'
	print(comment)
	outname = 'sdss-'+name+'.cat'
	Vizier.ROW_LIMIT    = -1
	query   = Vizier.query_region(coord.SkyCoord(ra=radeg, dec=dedeg, \
								unit=(u.deg, u.deg), frame='icrs'), \
								width=str(radius*60)+'m', catalog=["SDSS12"])
	querycat= query[query.keys()[0]]
	querycat.write(outname, format='ascii', overwrite=True)
	return querycat
#-------------------------------------------------------------------------#
def apass_query(name, radeg, dedeg, radius=1.0):
	"""
	APASS QUERY
	INPUT   :   NAME, RA [deg], Dec [deg], radius
	OUTPUT  :   QUERY TABLE
				apass-[NAME].cat
	"""
	from astroquery.vizier import Vizier 
	from astropy.coordinates import Angle
	from astropy.table import Table
	import astropy.units as u
	import astropy.coordinates as coord
	comment = 'NAME'+'\t'+': '+name+'\n' \
			+ 'RA'+'\t'+': '+str(round(radeg, 3))+'\n' \
			+ 'Dec'+'\t'+': '+str(round(dedeg, 3))+'\n' \
			+ 'Radius'+'\t'+': '+str(radius)+' deg'+'\n'*2 \
			+ 'LOADING APASS Catalog ...'+'\n'
	print(comment)
	outname = 'apass-'+name+'.cat'
	Vizier.ROW_LIMIT    = -1
	query   = Vizier.query_region(coord.SkyCoord(ra=radeg, dec=dedeg, \
								unit=(u.deg, u.deg), frame='icrs'), \
								width=str(radius*60)+'m', catalog=["APASS9"])
	dum     = query[0]
	colnames= dum.colnames
	for col in colnames:
		indx    = np.where( dum[col].mask == False )
		dum     = dum[indx]
	#   Vega    : B, V
	#   AB      : g, r, i
	#   Vega - AB Magnitude Conversion (Blanton+07)
	#   U       : m_AB - m_Vega = 0.79
	#   B       : m_AB - m_Vega =-0.09
	#   V       : m_AB - m_Vega = 0.02
	#   R       : m_AB - m_Vega = 0.21
	#   I       : m_AB - m_Vega = 0.45
	#   J       : m_AB - m_Vega = 0.91
	#   H       : m_AB - m_Vega = 1.39
	#   K       : m_AB - m_Vega = 1.85
	querycat			= Table()
	querycat['NUMBER']  = dum['recno']
	querycat['RA_ICRS'] = dum['RAJ2000']
	querycat['DE_ICRS'] = dum['DEJ2000']
	querycat['Numb_obs']= dum['nobs']
	querycat['Numb_img']= dum['mobs']
	querycat['B-V']     = dum['B-V']    + (-0.09 - 0.02)
	querycat['e_B-V']   = dum['e_B-V']
	querycat['Bmag']    = dum['Bmag']   - 0.09  # [Vega] to [AB]
	querycat['e_Bmag']  = dum['e_Bmag']
	querycat['Vmag']    = dum['Vmag']   + 0.02  # [Vega] to [AB]
	querycat['e_Vmag']  = dum['e_Vmag']
	querycat['gmag']    = dum['g_mag']
	querycat['e_gmag']  = dum['e_g_mag']
	querycat['rmag']    = dum['r_mag']
	querycat['e_rmag']  = dum['e_r_mag']
	querycat['imag']    = dum['i_mag']
	querycat['e_imag']  = dum['e_i_mag']
	
	querycat.write(outname, format='ascii', overwrite=True)
	return querycat
#-------------------------------------------------------------------------#
def ps1_query(name, radeg, dedeg, radius=1.0):
	"""
	#	SELECT STARS FROM STARS & GALAXIES (iPSF - iKron <= 0.05)
	https://outerspace.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies#
	"""
	from astroquery.vizier import Vizier 
	from astropy.coordinates import Angle
	from astropy.table import Table
	from astroquery.vizier import Vizier
	import astropy.units as u
	import astropy.coordinates as coord
	import numpy as np
	comment = 'NAME'+'\t'+': '+name+'\n' \
			+ 'RA'+'\t'+': '+str(round(radeg, 3))+'\n' \
			+ 'Dec'+'\t'+': '+str(round(dedeg, 3))+'\n' \
			+ 'Radius'+'\t'+': '+str(radius)+' deg'+'\n'*2 \
			+ 'LOADING PS1 Catalog ...'+'\n'
	print(comment)
	outname	= 'ps1-'+name+'.cat'
	#	QUERY PART
	Vizier.ROW_LIMIT    = -1
	query   = Vizier.query_region(coord.SkyCoord(ra=radeg, dec=dedeg, \
								unit=(u.deg, u.deg), frame='icrs'), \
								width=str(radius*60)+'m', catalog=["II/349/ps1"])
	dum0    = query[0]
	colnames= dum0.colnames
	#	REMOVE MASKED VALUE ROW
	for col in colnames:
		indx    = np.where( dum0[col].mask == False )
		dum1    = dum0[indx]
	f_objID_bin		= []
	for i in dum1['f_objID']:
		f_objID_bin.append(bin(i)[2:])
	f_objID_bin		= np.array( f_objID_bin )
	#	SELECT POINT SOURCE & NON-VARIABLE & GOOD QUALITY STARS
	indx_sel		= []
	for j in range(len(f_objID_bin)):
		i	= f_objID_bin[j]
		#	REJECT EXTENDED SOURCES THAT CONFIMED BY PS1 & 2MASS (23, 24)
		#	REJECT QSO, RR Lyra, VARIABLE, TRANSIENT (2, 3, 4, 5, 6, 7, 8)
		#	REJECT POOR-QUALITY STACK OBJECT (30 = 0) -> not applied yet
		try:
			if (i[-23] != '1') and (i[-24] != '1') and (i[-2] != '1') and (i[-3] != '1') and (i[-4] != '1') and (i[-5] != '1') and (i[-6] != '1') and (i[-7] != '1') and (i[-8] != '1'):# and (i[0] != '1'):
				indx_sel.append(j)
		except:
			pass
	dum2	= dum1[indx_sel]
	#	SELECT STARS FROM STARS & GALAXIES (iPSF - iKron <= 0.05)
	indx_stars		= np.where( (dum2['imag'] - dum2['iKmag']) <= 0.05 )
	dum		= dum2[indx_stars]
	#	CHANGE TO GENERTAL COL. NAMES
	querytbl			= Table()
	querytbl['NUMBER']  = dum['objID']
	querytbl['RA_ICRS'] = dum['RAJ2000']
	querytbl['DE_ICRS'] = dum['DEJ2000']
	querytbl['Q']		= dum['Qual']
	querytbl['Numb_obs']= dum['Nd']
	querytbl['Numb_img']= dum['Ns']
	querytbl['gmag']    = dum['gmag']
	querytbl['e_gmag']  = dum['e_gmag']
	querytbl['rmag']    = dum['rmag']
	querytbl['e_rmag']  = dum['e_rmag']
	querytbl['imag']    = dum['imag']
	querytbl['e_imag']  = dum['e_imag']
	querytbl['zmag']    = dum['zmag']
	querytbl['e_zmag']  = dum['e_zmag']
	querytbl['ymag']    = dum['ymag']
	querytbl['e_ymag']  = dum['e_ymag']

	querytbl.write(outname, format='ascii', overwrite=True)
	return querytbl
#-------------------------------------------------------------------------#
def twomass_query(name, radeg, dedeg, band=None, radius=1.0):
	"""
	QUERY Point Source Catalog(PSC) PROVIDED BY 2MASS
	REMOVE COMTAMINATED SOURCE BY EXTENDED SOURCE AND MINOR PLANET
	IF GIVE BAND INPUT, 

	INPUT	:	NAME, RA [deg], DEC [deg], BAND, RADIUS
	OUTPUT	:	TABLE, MAGNITUDE [AB]
	
	"""
	from astroquery.vizier import Vizier 
	from astropy.coordinates import Angle
	from astropy.table import Table
	from astroquery.vizier import Vizier
	import astropy.units as u
	import astropy.coordinates as coord
	import numpy as np
	comment = 'NAME'+'\t'+': '+name+'\n' \
			+ 'RA'+'\t'+': '+str(round(radeg, 3))+'\n' \
			+ 'Dec'+'\t'+': '+str(round(dedeg, 3))+'\n' \
			+ 'Radius'+'\t'+': '+str(radius)+' deg'+'\n'*2 \
			+ 'LOADING 2MASS Catalog ...'+'\n'
	print(comment)
	outname	= '2mass-'+name+'.cat'
	#	QUERY PART
	Vizier.ROW_LIMIT    = -1
	query   = Vizier.query_region(coord.SkyCoord(ra=radeg, dec=dedeg, \
								unit=(u.deg, u.deg), frame='icrs'), \
								width=str(radius*60)+'m', catalog=["II/246"])
	dum0    = query[0]
	colnames= dum0.colnames
	['RAJ2000', 'DEJ2000', '_2MASS',
	 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag',
	 'Qflg', 'Rflg', 'Bflg', 'Cflg', 'Xflg', 'Aflg']
	#	REMOVE MASKED VALUE ROW
	for col in colnames:
		indx    = np.where( dum0[col].mask == False )
		dum1    = dum0[indx]
	#	SELECT POINT SOURCE & NON-VARIABLE & GOOD QUALITY STARS
	if	band	== None:
		indx_flg		= np.where(	(dum1['Aflg'] == 0) &
									(dum1['Xflg'] == 0)	)
	else:
		if band	== 'J':
			order	= 0
		if band	== 'H':
			order	= 1
		if band	== 'K':
			order	= 2
		indx_flg		= []
		for i in range(len(dum1)):
			Qflg	= dum1['Qflg'][i]
			Rflg	= dum1['Rflg'][i]
			Bflg	= dum1['Bflg'][i]
			Cflg	= dum1['Cflg'][i]
			Xflg	= dum1['Xflg'][i]
			Aflg	= dum1['Aflg'][i]
			if	(	(Qflg[order]	== 'A')		|
					(Qflg[order]	== 'B')		|
					(Qflg[order]	== 'C'))	& \
				(	(Bflg[order]	!= '0'))	& \
				(	(Cflg[order]	== '0'))	& \
				(	(Xflg			== 0))		& \
				(	(Aflg			== 0)	):
				indx_flg.append(i)
		indx_flg	= np.array( list(set(indx_flg)) )
	dum				= dum1[indx_flg]
	#	CHANGE TO GENERTAL COL. NAMES
	querytbl			= Table()
	querytbl['name']  	= dum['_2MASS']
	querytbl['ra'] 		= dum['RAJ2000']
	querytbl['dec'] 	= dum['DEJ2000']
	#	AB OFFSET
	querytbl['J']    	= dum['Jmag']	+ 0.91
	querytbl['Jerr']  	= dum['e_Jmag']
	querytbl['H']    	= dum['Hmag']	+ 1.39
	querytbl['Herr']  	= dum['e_Hmag']
	querytbl['K']    	= dum['Kmag']	+ 1.85
	querytbl['Kerr']  	= dum['e_Kmag']
	querytbl['Qflg']	= dum['Qflg']
	querytbl['Rflg']	= dum['Rflg']
	querytbl['Bflg']	= dum['Bflg']
	querytbl['Cflg']	= dum['Cflg']

	querytbl.write(outname, format='ascii', overwrite=True)
	return querytbl
#-------------------------------------------------------------------------#
def secom(inim, gain, pixscale, det_sigma=3, backsize=str(64), backfiltersize=str(3), dual=False, detect='detection.fits', check=False):
	"""
	SourceEXtractor
	APERTURE    3", 5", 7",
				1.0seeing, 1.2seeing ,1.5seeing ,1.7seeing ,2.0seeeing
	INPUT   :   (image).fits
				aperture    []
				seeing_fwhm [pixel]
	OUTPUT  :   no return
				.cat
	"""
	import numpy as np
	import os
	from astropy.io import ascii
	#   FILE CHECK
	#	CONFIG FILES (USER BASE PATH)
	configfile      = '/home/sonic/Research/yourpy/config/targetphot.sex'
	paramfile       = '/home/sonic/Research/yourpy/config/targetphot.param'
	nnwfile		    = '/home/sonic/Research/yourpy/config/targetphot.nnw'
	convfile	    = '/home/sonic/Research/yourpy/config/targetphot.conv'
	try:
		comment = '\nSourceEXtractor START\n' \
				+ 'IMAGE\t\t: '+inim+'\n' \
				+ 'GAIN\t\t: '+str(gain)+'\n' \
				+ 'PIXSCALE\t\t: '+str(pixscale)+'\n' \
				+ 'DETECTION SIGMA\t: '+str(det_sigma)+'\n' \
				+ 'PARAM\t\t: '+paramfile+'\n' \
				+ 'BACKSIZE\t: '+backsize+'\n' \
				+ 'BACKFILTER\t: '+backfiltersize+'\n' \
				+ 'CONFIG\t\t: '+configfile+'\n' \
				+ 'NNW\t\t: '+nnwfile+'\n' \
				+ 'CONVOLVE\t: '+convfile
		print(comment)
	except:
		comment = 'CHECK configfile/paramfile/nnewfile/convfile or others.'
		print(comment)
	#   FILE NAME
	cat     = inim[:-5]+'.cat'
	seg     = inim[:-5]+'.seg.fits'
	bkg     = inim[:-5]+'.bkg.fits'
	sub     = inim[:-5]+'.sub.fits'
	psf     = inim[:-5]+'.psf'
	aper    = inim[:-5]+'.aper.fits'

	#   BASIC INFO.
	det_area        = 5
	det_thresh      = det_sigma/np.sqrt(det_area)
	detecminarea    = str(det_area)
	detectthresh    = str(det_thresh)
	seeing, fwhm_arcsec = psfex(inim, pixscale)
	#pixscale        = pixscalecalc(inim)
	#   OPTION
	aperture        = '%.2f'%(3./pixscale)+','+'%.2f'%(5./pixscale)+','+'%.2f'%(7./pixscale)+','+'%.2f'%(seeing)+','+'%.2f'%(1.2*seeing)+','+'%.2f'%(1.5*seeing)+','+'%.2f'%(1.7*seeing)+','+'%.2f'%(2.0*seeing)

	option1 = ' -CATALOG_NAME '+cat+' -PARAMETERS_NAME '+paramfile
	option2 = ' -DETECT_MINAREA '+detecminarea+' -DETECT_THRESH '+detectthresh \
			+' -FILTER Y '+'-FILTER_NAME '+convfile
	option3 = ' -PHOT_APERTURES '+aperture \
			+' -GAIN '+'%.1f'%(gain)+' -PIXEL_SCALE '+'%.1f'%(pixscale)
	option4 = ' -SEEING_FWHM '+'%.3f'%(seeing)+' -STARNNW_NAME '+nnwfile
	option5 = ' -BACK_SIZE '+ backsize \
			+ ' -BACK_FILTERSIZE '+ backfiltersize+' -BACKPHOTO_TYPE LOCAL'
	option6 = ' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND'
	option7 = ' -CHECKIMAGE_NAME '+seg+','+aper+','+','+bkg+','+sub
	option8 = ' -PSF_NAME '+psf
	#   COMMAND
	#   detect = detection.fits is fine image show target significantly and have good seeing and many stars i
	dualcom ='sex -c '+configfile+' '+detect+' , '+inim+' -CATALOG_NAME dual'+cat+' -PARAMETERS_NAME '+paramfile+ ' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8
	sglcom  ='sex -c '+configfile+' '+inim+' '+option1+' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8
	clearcom='sex -c '+configfile+' '+inim+' '+option1+' '+option2+' '+option3+' '+option4+' '+option5+' '+option8
	if dual == False    :
		if check == False:
			os.system(clearcom)
		else:
			os.system(sglcom)
	else                : os.system(dualcom)
	
	secat   = ascii.read(cat)
	return secat, cat, seeing, fwhm_arcsec
#-------------------------------------------------------------------------#
def psfex(inim, pixscale):
	"""
	PSfextractor
	INPUT   :   (image).fits
	OUTPUT  :   FWHM    [pixel]
				FWHM    [arcsec]
	"""
	import os
   
	#   FILE CHECK
	#	CONFIG FILES (USER BASE PATH)
	psfexconf_prese_conf    = '/home/sonic/Research/yourpy/config/prepsfex.sex'
	psfexconf_prese_param   = '/home/sonic/Research/yourpy/config/prepsfex.param'
	psfexconf_psfex_conf    = '/home/sonic/Research/yourpy/config/default.psfex'
	psfexconf_psfex_conv    = '/home/sonic/Research/yourpy/config/default.conv'
	try:
		comment = '\nPSFex START\n' \
				+ 'IMAGE\t\t: '+inim+'\n' \
				+ 'PRE_CONFIG\t: '+psfexconf_prese_conf+'\n' \
				+ 'PRE_PARAM\t: '+psfexconf_prese_param+'\n' \
				+ 'CONFIG\t\t: '+psfexconf_psfex_conf+'\n' \
				+ 'CONV\t\t: '+psfexconf_psfex_conv
		print(comment)
	except:
		comment = 'CHECK psfexconf_prese/psfexconf_prese_param/psfexconf_psfex_conf/psfexconf_psfex_conv OR OTHERS.'
		print(comment)

	#   FILE NAME
	cat     = inim[:-5]+'.cat'
	xml     = inim[:-5]+'.xml'
	snap    = 'snap_'+inim+'[100:125,100:125]'
	psf     = 'psf-'+inim
	#   OPTION
	presecom1   = psfexconf_prese_conf+" "+inim
	presecom2   = " -CATALOG_NAME "+cat
	presecom3   = " -FILTER_NAME " + psfexconf_psfex_conv + " -PARAMETERS_NAME " + psfexconf_prese_param
	#   COMMAND
	presecom    = "sex -c "+presecom1+presecom2+presecom3
	psfexcom    = "psfex -c "+psfexconf_psfex_conf+" "+cat
	os.system(presecom)
	os.system(psfexcom) 
	os.system('cp psfex.xml '+xml)
	#   SNAP IMAGE
	imcopycom   = 'imcopy '+snap+' '+psf
	print(imcopycom);   os.system(imcopycom)
	#   FWHM [pixel], FWHM [arcsec]
	fwhm_pix    = psfexxml(xml)
	fwhm_arcsec = round(fwhm_pix*pixscale, 3)
	comment     = '\n' \
				+ 'FILE NAME'+'\t'+': '+inim+'\n' \
				+ 'FWHM value'+'\t'+': '+str(fwhm_pix)+'\t'+'[pixel]'+'\n' \
				+ '\t'+'\t'+': '+str(fwhm_arcsec)+'\t'+'[arcsec]'+'\n'
	print(comment)
	return fwhm_pix, fwhm_arcsec
#-------------------------------------------------------------------------#
def sdss_Blaton(intbl, name):
	"""
	=====================================================================
	CONVERSION SDSS FILTER SYSTEM TO JOHNSON FILTER SYSTEM [AB]
	INPUT   :   QUERIED SDSS CATALOG
	OUTPUT  :   1.  ONLY STAR CLASS (=6)
				2.  (2)ACCEPTABLE & (3)GOOD QUALITY
				3.  NAME, RA, Dec, ugriz, BVRI [AB] (Blaton+07), clean
				sdss-conv.cat
	---------------------------------------------------------------------
	Blaton+07
	CONVERSION TABLE to AB
	---------------------------------------------------------------------
	Equation                                        Color Dispersion
	---------------------------------------------------------------------
	U   = u - 0.0682 - 0.0140[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	B   = u - 1.0286 - 0.7981[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	*B   = g + 0.2354 + 0.3915[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	V   = g - 0.3516 - 0.7585[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	R   = r - 0.0576 - 0.3718[ (r - i) - 0.2589 ]   sigma[r - i] = 0.10
	I   = i - 0.0647 - 0.7177[ (i - z) - 0.2083 ]   sigma[i - z] = 0.10
	=====================================================================
	"""
	import numpy as np
	from astropy.table import Table
	from astropy.io import ascii
	
	outfile = 'sdss-Blaton-'+name+'.cat'
	
	clas, Q = intbl['class'],       intbl['Q']
	indx    = np.where( (Q != 1) & (clas == 6) )
	reftbl  = intbl[indx]
	
	name    = reftbl['SDSS12']
	ra, de  = reftbl['RA_ICRS'],    reftbl['DE_ICRS']
	clean   = reftbl['q_mode']

	u, uer  = reftbl['umag'],       reftbl['e_umag']
	g, ger  = reftbl['gmag'],       reftbl['e_gmag']
	r, rer  = reftbl['rmag'],       reftbl['e_rmag']
	i, ier  = reftbl['imag'],       reftbl['e_imag']
	z, zer  = reftbl['zmag'],       reftbl['e_zmag']

	uger, grer, rier, izer	= 0.26, 0.15, 0.10, 0.10
	ug, gr, ri, iz			= u-g, g-r, r-i, i-z
	'''
	U		= u - 0.0682 - 0.0140*(ug-1.2638)
	Uer		= np.sqrt( ((uer)**2.) + ((-0.0140*uger)**2.) )
	'''
	B		= g + 0.2354 + 0.3915*(gr-0.6102)
	Ber		= np.sqrt( ((ger)**2.) + ((0.3915*grer)**2.) )
	V		= g - 0.3516 - 0.7585*(gr-0.6102)
	Ver		= np.sqrt( ((ger)**2.) + ((-0.7585*grer)**2.) )
	R		= r - 0.0576 - 0.3718*(ri-0.2589)
	Rer		= np.sqrt( ((rer)**2.) + ((-0.3718*rier)**2.) )
	I		= i - 0.0647 - 0.7177*(iz-0.2083)
	Ier		= np.sqrt( ((ier)**2.) + ((-0.7177*izer)**2.) )

	outtbl	= Table([name, ra, de, u, uer, g, ger, r, rer, i, ier, z, zer, B, Ber, V, Ver, R, Rer, I, Ier, clean], names=['name', 'ra', 'dec', 'u', 'uerr', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'clean'])
	
	outtbl0	= Table([name, ra, de, u, uer, g, ger, r, rer, i, ier, z, zer, B, Ber, V, Ver, R, Rer, I, Ier, clean], names=['#name', 'ra', 'dec', 'u', 'uerr', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'clean'])
	ascii.write(outtbl0, outfile)#, format='fixed_width', delimiter=' ')
	return outtbl, outfile
#-------------------------------------------------------------------------#
def apass_Blaton(intbl, name):
	"""
	=====================================================================
	CONVERSION SDSS FILTER SYSTEM TO JOHNSON FILTER SYSTEM [AB]
	INPUT   :   QUERIED SDSS CATALOG
	OUTPUT  :   1.  ONLY STAR CLASS (=6)
				2.  (2)ACCEPTABLE & (3)GOOD QUALITY
				3.  NAME, RA, Dec, ugriz, BVRI [AB] (Blaton+07), clean
				sdss-conv.cat
	---------------------------------------------------------------------
	Blaton+07
	CONVERSION TABLE to AB
	---------------------------------------------------------------------
	Equation                                        Color Dispersion
	---------------------------------------------------------------------
	U   = u - 0.0682 - 0.0140[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	B   = u - 1.0286 - 0.7981[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	*B   = g + 0.2354 + 0.3915[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	V   = g - 0.3516 - 0.7585[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	R   = r - 0.0576 - 0.3718[ (r - i) - 0.2589 ]   sigma[r - i] = 0.10
	I   = i - 0.0647 - 0.7177[ (i - z) - 0.2083 ]   sigma[i - z] = 0.10
	=====================================================================
	"""
	import numpy as np
	from astropy.table import Table
	from astropy.io import ascii
	
	outfile = 'apass-Blaton-'+name+'.cat'
	
	reftbl	= intbl
	
	name    = reftbl['NUMBER']
	ra, de  = reftbl['RA_ICRS'],    reftbl['DE_ICRS']
	Numb_obs= reftbl['Numb_obs']
	Numb_img= reftbl['Numb_img']
	B		= reftbl['Bmag']
	Ber		= reftbl['e_Bmag']
	V		= reftbl['Vmag']
	Ver		= reftbl['e_Vmag']	
	BV		= reftbl['B-V']
	e_BV	= reftbl['e_B-V']

	g, ger  = reftbl['gmag'],       reftbl['e_gmag']
	r, rer  = reftbl['rmag'],       reftbl['e_rmag']
	i, ier  = reftbl['imag'],       reftbl['e_imag']

	grer, rier		= 0.15, 0.10
	gr, ri			= g-r, r-i

	R		= r - 0.0576 - 0.3718*(ri-0.2589)
	Rer		= np.sqrt( ((rer)**2.) + ((-0.3718*rier)**2.) )
	'''
	I		= i - 0.0647 - 0.7177*(iz-0.2083)
	Ier		= np.sqrt( ((ier)**2.) + ((-0.7177*izer)**2.) )
	'''
	outtbl	= Table([name, ra, de, Numb_obs, Numb_img, g, ger, r, rer, i, ier, B, Ber, V, Ver, R, Rer], names=['name', 'ra', 'dec', 'numb_obs', 'numb_img', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr'])
	
	outtbl0 = Table([name, ra, de, Numb_obs, Numb_img, g, ger, r, rer, i, ier, B, Ber, V, Ver, R, Rer], names=['#name', 'ra', 'dec', 'numb_obs', 'numb_img', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr'])
	ascii.write(outtbl0, outfile)#, format='fixed_width', delimiter=' ')
	return outtbl, outfile
#-------------------------------------------------------------------------#
def apass_Lupton(intbl, name):
	"""
	=====================================================================
	CONVERSION SDSS FILTER SYSTEM TO JOHNSON FILTER SYSTEM [AB]
	INPUT   :   QUERIED SDSS CATALOG
	OUTPUT  :   1.  ONLY STAR CLASS (=6)
				2.  (2)ACCEPTABLE & (3)GOOD QUALITY
				3.  NAME, RA, Dec, ugriz, BVRI, clean
				sdss-conv.cat
	---------------------------------------------------------------------
	Blaton+07
	CONVERSION TABLE to AB
	---------------------------------------------------------------------
	Equation                                        Color Dispersion
	---------------------------------------------------------------------
	U   = u - 0.0682 - 0.0140[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	B   = u - 1.0286 - 0.7981[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	*B   = g + 0.2354 + 0.3915[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	V   = g - 0.3516 - 0.7585[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	R   = r - 0.0576 - 0.3718[ (r - i) - 0.2589 ]   sigma[r - i] = 0.10
	I   = i - 0.0647 - 0.7177[ (i - z) - 0.2083 ]   sigma[i - z] = 0.10
	=====================================================================
	"""
	import numpy as np
	from astropy.table import Table
	from astropy.io import ascii
	
	outfile = 'apass-Lupton-'+name+'.cat'
	
	reftbl	= intbl
	
	name    = reftbl['NUMBER']
	ra, de  = reftbl['RA_ICRS'],    reftbl['DE_ICRS']
	Numb_obs= reftbl['Numb_obs']
	Numb_img= reftbl['Numb_img']
	B		= reftbl['Bmag']
	Ber		= reftbl['e_Bmag']
	V		= reftbl['Vmag']
	Ver		= reftbl['e_Vmag']	
	BV		= reftbl['B-V']
	e_BV	= reftbl['e_B-V']

	g, ger  = reftbl['gmag'],       reftbl['e_gmag']
	r, rer  = reftbl['rmag'],       reftbl['e_rmag']
	i, ier  = reftbl['imag'],       reftbl['e_imag']

	Rsig, Isig  = 0.0072, 0.0078

	R       = r - 0.2936*(r - i) - 0.1439
	Rer1    = ((1.-0.2936)**2)*(ger**2)+((+0.2936)**2)*(rer**2)
	Rer     = np.sqrt(Rer1**2 + Rsig**2)
	
	I       = r - 1.2444*(r - i) - 0.3820
	Ier1    = ((1.-1.2444)**2)*(ger**2)+((+1.2444)**2)*(rer**2)
	Ier     = np.sqrt(Ier1**2 + Isig**2)

	#	Vega to AB
	R2AB, I2AB  = R + 0.21, I + 0.45

	outtbl	= Table([name, ra, de, Numb_obs, Numb_img, g, ger, r, rer, i, ier, B, Ber, V, Ver, R2AB, Rer, I2AB, Ier], names=['name', 'ra', 'dec', 'numb_obs', 'numb_img', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr'])
	
	outtbl0 = Table([name, ra, de, Numb_obs, Numb_img, g, ger, r, rer, i, ier, B, Ber, V, Ver, R2AB, Rer, I2AB, Ier], names=['#name', 'ra', 'dec', 'numb_obs', 'numb_img', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr'])
	ascii.write(outtbl0, outfile)#, format='fixed_width', delimiter=' ')
	return outtbl, outfile
#-------------------------------------------------------------------------#
def ps1_Tonry(intbl, name):
	'''
	PS1 -> Johnson/COusins [Vega] -> [AB]	(Tonry+12)
	#   Vega - AB Magnitude Conversion (Blanton+07)
	U       : m_AB - m_Vega = 0.79
	B       : m_AB - m_Vega =-0.09
	V       : m_AB - m_Vega = 0.02
	R       : m_AB - m_Vega = 0.21
	I       : m_AB - m_Vega = 0.45
	J       : m_AB - m_Vega = 0.91
	H       : m_AB - m_Vega = 1.39
	K       : m_AB - m_Vega = 1.85
	'''
	import numpy as np
	from astropy.table import Table
	from astropy.io import ascii
	#	REJECT BAD QUALITY
	intbl	= intbl[	(intbl['gmag']>12)	&
						(intbl['rmag']>12)	&
						(intbl['imag']>12)	&
						(intbl['zmag']>12)	&
						(intbl['ymag']>12)]


	outfile	= 'ps1-Tonry-'+name+'.cat'
	Q		= intbl['Q']
	indx    = np.where(Q < 128)

	intbl	= intbl[indx]
	Q		= intbl['Q']
	
	name    = intbl['NUMBER']
	ra, de  = intbl['RA_ICRS'],    intbl['DE_ICRS']

	g, ger  = intbl['gmag'],	intbl['e_gmag']
	r, rer  = intbl['rmag'],	intbl['e_rmag']
	i, ier  = intbl['imag'],	intbl['e_imag']
	z, zer  = intbl['zmag'],	intbl['e_zmag']
	y, yer  = intbl['ymag'],	intbl['e_ymag']
	#	TRANSF. ERROR FOR B CONST. TERMS
	Bsig, Vsig, Rsig, Isig	= 0.034, 0.012, 0.01, 0.016
	#	COLOR TERM
	gr		= intbl['gmag']-intbl['rmag']
	grer	= sqsum(intbl['e_gmag'], intbl['e_rmag'])
	#	CONVERT TO B
	B0		= 0.213
	B1		= 0.587
	B		= B0 + B1*gr + intbl['gmag'] - 0.09
	Ber		= sqsum( Bsig, sqsum(B1*grer, intbl['e_gmag']) )
	#	CONVERT TO V
	B0		= 0.006
	B1		= 0.474
	V		= B0 + B1*gr + intbl['rmag'] + 0.02
	Ver	= sqsum( Bsig, sqsum(B1*grer, intbl['e_rmag']) )
	#	CONVERT TO R
	B0		=-0.138
	B1		=-0.131
	R		= B0 + B1*gr + intbl['rmag'] + 0.21
	Rer		= sqsum( Rsig, sqsum(B1*grer, intbl['e_rmag']) )
	#	CONVERT TO I
	B0		=-0.367
	B1		=-0.149
	I		= B0 + B1*gr + intbl['imag'] + 0.45
	Ier		= sqsum( Isig, sqsum(B1*grer, intbl['e_imag']) )
	outtbl	= Table([name, ra, de, g, ger, r, rer, i, ier, z, zer, y, yer, B, Ber, V, Ver, R, Rer, I, Ier, Q], names=['name', 'ra', 'dec', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'y', 'yerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'Q'])
	outtbl0	= Table([name, ra, de, g, ger, r, rer, i, ier, z, zer, y, yer, B, Ber, V, Ver, R, Rer, I, Ier, Q], names=['#name', 'ra', 'dec', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'y', 'yerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'Q'])
	ascii.write(outtbl0, outfile)#, format='fixed_width', delimiter=' ')
	return outtbl, outfile
#-------------------------------------------------------------------------#
def psfexxml(xmlfile):
	"""
	INPUT   :   .xml
	OUTPUT  :   FWHM    [pixel]
	"""
	from astropy.io.votable import parse
	from astropy.table import Table, Column, MaskedColumn
	votable     = parse(xmlfile)
	table       = votable.get_first_table()
	data        = table.array
	#   EXTRACT FWHM [pixel]
	fwhm        = data['FWHM_Mean'][0]
	fwhm        = round(fwhm, 3)
	return fwhm
#-------------------------------------------------------------------------#
def matching(incat, refcat, sep=2.0):
	"""
	MATCHING TWO CATALOG WITH RA, Dec COORD. WITH python
	INPUT   :   SE catalog, SDSS catalog file name, sepertation [arcsec]
	OUTPUT  :   MATCED CATALOG FILE & TABLE
	"""
	import numpy as np
	import astropy.units as u
	from astropy.table import Table, Column
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii
	
	mergename   = incat[:-4]+'.merge.cat'
	if type(incat) == str: intbl       = ascii.read(incat)
	else: intbl	= incat
	if type(refcat) == str: reftbl      = ascii.read(refcat)
	else: reftbl = refcat
	coo_intbl   = SkyCoord(intbl['ALPHA_J2000'], intbl['DELTA_J2000'], unit=(u.deg, u.deg))
	coo_reftbl  = SkyCoord(reftbl['ra'], reftbl['dec'], unit=(u.deg, u.deg))

	#   INDEX FOR REF.TABLE
	indx, d2d, d3d  = coo_intbl.match_to_catalog_sky(coo_reftbl)
	ref_match       = reftbl[indx]
	ref_match['sep']= d2d
	ref_match_col   = ref_match.colnames

	merge_tbl       = intbl

	for col in ref_match.colnames:
		merge_tbl[col]  = ref_match[col]

	indx_cut        = np.where(merge_tbl['sep']*3600. < sep)
	merge           = merge_tbl[indx_cut]

	merge.write(mergename, format='ascii', overwrite=True)
	comment         = mergename+' is generagted.'
	print(comment)
	return  merge
#-------------------------------------------------------------------------#
def star4zp(intbl, inmagerkey, refmagkey, refmagerkey, 
			refmaglower=14., refmagupper=17., refmagerupper=0.05,
			inmagerupper=0.1, class_star_cut=0.5):
	"""
	SELECT STARS FOR USING ZEROPOINT CALCULATION
	INPUT   :   TABLE, IMAGE MAG.ERR KEYWORD, REF.MAG. KEYWORD, REF.MAG.ERR KEYWORD
	OUTPUT  :   NEW TABLE
	"""
	import numpy as np
	indx    = np.where( (intbl['FLAGS'] == 0) & 
						(intbl['CLASS_STAR'] >= class_star_cut) & 
						(intbl[refmagkey] < refmagupper) & 
						(intbl[refmagkey] > refmaglower) & 
						(intbl[refmagerkey] < refmagerupper) &
						(intbl[inmagerkey] < inmagerupper) 
						)
	indx0   = np.where( (intbl['FLAGS'] == 0) )
	indx1   = np.where( (intbl['CLASS_STAR'] >= class_star_cut) )
	indx2   = np.where( (intbl[refmagkey] < refmagupper) & 
						(intbl[refmagkey] > refmaglower) & 
						(intbl[refmagerkey] < refmagerupper) 
						)
	indx3   = np.where( (intbl[inmagerkey] < inmagerupper) )
	newtbl  = intbl[indx]
	comment = '='*60+'\n' \
			+ 'ALL\t\t\t\t: '+str(len(intbl))+'\n' \
			+ '-'*60+'\n' \
			+ 'FLAG(=0)\t\t\t: '+str(len(indx0[0]))+'\n' \
			+ 'CLASS_STAR > '+str(class_star_cut)+'\t\t: '+str(len(indx1[0]))+'\n' \
			+ refmagkey+' REF. MAGCUT ('+str(refmaglower)+'-'+str(refmagupper)+')'+'\t: '+str(len(indx2[0]))+'\n' \
			+ refmagerkey+' REF. MAGERR CUT < '+str(refmagerupper)+'\n' \
			+ inmagerkey+' OF IMAGE CUT < '+str(inmagerupper)+'\t: '+str(len(indx3[0]))+'\n' \
			+ '-'*60+'\n' \
			+ 'TOTAL #\t\t\t\t: '+str(len(indx[0]))+'\n' \
			+ '-'*60
	print(comment)
	return newtbl
#-------------------------------------------------------------------------#
def zpcal(intbl, inmagkey, inmagerkey, refmagkey, refmagerkey, sigma=2.0):
	"""
	ZERO POINT CALCULATION
	3 SIGMA CLIPPING (MEDIAN)

	"""
	import matplotlib.pyplot as plt
	from numpy import median
	from astropy.stats import sigma_clip
	import numpy as np
	#	REMOVE BLANK ROW (=99)	
	indx_avail      = np.where( (intbl[inmagkey] != 99) & (intbl[refmagkey] != 99) )
	intbl           = intbl[indx_avail]
	zplist          = intbl[refmagkey] - intbl[inmagkey]
	intbl['zp']		= zplist
	#	SIGMA CLIPPING
	zplist_clip     = sigma_clip(zplist, sigma=sigma, maxiters=None, cenfunc=median, copy=False)
	indx_alive      = np.where( zplist_clip.mask == False )
	indx_exile      = np.where( zplist_clip.mask == True )
	#	RE-DEF. ZP LIST AND INDEXING CLIPPED & NON-CLIPPED
	intbl_alive     = intbl[indx_alive]
	intbl_exile     = intbl[indx_exile]
	#	ZP & ZP ERR. CALC.
	zp              = np.median(intbl_alive['zp'])
	zper			= np.std(intbl_alive['zp'])
	return zp, zper, intbl_alive, intbl_exile
#-------------------------------------------------------------------------#
def zp_plot(inim, inmagkey, zp, zper, instmag0, refmag0, zplist0, instmag1, refmag1, zplist1):
	"""
	ZERO POINT PLOT
	INPUT   :	IMAGE, APERTURE KEY, ZP, ZP ERR,
				INST. MAG(NOT CLIPPED), REF. MAG(NOT CLIPPED), ZP LIST(NOT CLIPPED),
				INST. MAG(CLIPPED), REF. MAG(CLIPPED), ZP LIST(CLIPPED),
	OUTPUT  :
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	#   FILE NAME
	outname = inim[:-5]+'.'+inmagkey
	#	POINT
	plt.scatter( refmag0, zplist0, color='dodgerblue', s=15, marker='o', linewidth=1, alpha=0.5, label='NOT CLIPPED ('+str(len(refmag0))+')' )
	plt.scatter( refmag1, zplist1, color='tomato', s=30, marker='x', linewidth=1, alpha=0.5, label='CLIPPED ('+str(len(refmag1))+')' )
	#	LINE
	plt.axhline(zp, linewidth=1, linestyle='--', color='gray', label=str(round(zp, 3)) )
	plt.axhline(zp+zper, linewidth=1, linestyle='-', color='gray', alpha=0.5, label=str(round(zper, 3)) )
	plt.axhline(zp-zper, linewidth=1, linestyle='-', color='gray', alpha=0.5 )#, label=str(round(zp, 3)) )
	plt.fill_between([np.min(refmag0)-0.05, np.max(refmag0)+0.05], zp-zper, zp+zper, color='silver', alpha=0.3)
	#	RANGE
	plt.xlim(np.min(refmag0)-0.05, np.max(refmag0)+0.05)
	plt.ylim(zp-0.5, zp+0.5)
	#	SETTING
	plt.title(outname, {'fontsize': 10})
	plt.gca().invert_yaxis()
	plt.xlabel('REF.MAG.', {'color': 'black', 'fontsize': 20})
	plt.ylabel('ZERO POINT [AB]', {'color': 'black', 'fontsize': 20})
	plt.legend(loc='best', prop={'size': 11})
	plt.tight_layout()
	plt.savefig(outname+'.zpcal.png', dpi = 500)
	plt.close()

	magdif0		= (refmag0) - (instmag0+zp)
	magdif1		= (refmag1) - (instmag1+zp)
	#	POINT
	plt.scatter( refmag0, magdif0, color='dodgerblue', s=15, marker='o', linewidth=1, alpha=0.5, label='NOT CLIPPED' )
	plt.scatter( refmag1, magdif1, color='tomato', s=30, marker='x', linewidth=1, alpha=0.5, label='CLIPPED' )
	#	LINE
	plt.axhline(0, linewidth=1, linestyle='--', color='gray', label=str(round(zp, 3)) )
	plt.axhline(zper, linewidth=1, linestyle='-', color='gray', alpha=0.5, label=str(round(zper, 3)) )
	plt.axhline(-zper, linewidth=1, linestyle='-', color='gray', alpha=0.5)#, label=str(round(zper, 3)) )
	plt.fill_between([np.min(refmag0)-0.05, np.max(refmag0)+0.05], -zper, +zper, color='silver', alpha=0.3)
	#	RANGE
	plt.xlim(np.min(refmag0)-0.05, np.max(refmag0)+0.05)
	plt.ylim(-0.5, +0.5)
	#	SETTING
	plt.title(outname, {'fontsize': 10})
	plt.xlabel('REF.MAG.', {'color': 'black', 'fontsize': 20})
	plt.ylabel('REAL MAG.-(ZP+INST.MAG) [AB]', {'color': 'black', 'fontsize': 10})
	plt.legend(loc='best', prop={'size': 11})
	plt.tight_layout()
	plt.savefig(outname+'.zpcal_test.png', dpi = 500)
	plt.close()
	#	PRINT
	print('MAG TYP\t\t: '+inmagkey)
	print('ZP\t\t: '+str(round(zp, 3)))
	print('ZP ERR\t\t: '+str(round(zper, 3)))
	print('STD.NUMB\t: '+str(int(len(zplist0))))
	print('REJ.NUMB\t: '+str(int(len(zplist1))))
#-------------------------------------------------------------------------#
def bkgest_mask(inim):
	'''
	'''
	import numpy as np
	from astropy.io import fits
	from photutils import make_source_mask
	from numpy import mean,median
	from astropy.stats import sigma_clipped_stats
	
	data    =fits.getdata(inim)
	mask    = make_source_mask(data, snr=3, npixels=5, dilate_size=11)
	mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
	return mean, median, std
#-------------------------------------------------------------------------#
def limitmag(N, zp, aper, skysigma):			# 3? 5?, zp, diameter [pixel], skysigma
	import numpy as np
	R           = float(aper)/2.				# to radius
	braket      = N*skysigma*np.sqrt(np.pi*(R**2))
	upperlimit  = float(zp)-2.5*np.log10(braket)
	return round(upperlimit, 3)
#-------------------------------------------------------------------------#
def puthdr(inim, hdrkey, content):
	test	= 'MAKING'
	return test
#-------------------------------------------------------------------------#
def sqsum(a, b):
	'''
	SQUARE SUM
	USEFUL TO CALC. ERROR
	'''
	import numpy as np
	return np.sqrt(a**2.+b**2.)
#-------------------------------------------------------------------------#
def targetfind(ra1, de1, ra2, de2, sep):
	'''
	FIND TARGET
	'''
	import numpy as np
	dist	= np.sqrt( (ra1-ra2)**2. + (de1-de2)**2. )
	indx	= np.where( (dist == np.min(dist)) &
						(dist < sep/3600.) )
	return indx
#-------------------------------------------------------------------------#
def plotshow(inim, numb_list, xim_list, yim_list, add=None, numb_addlist=None, xim_addlist=None, yim_addlist=None):
	'''
	PLOT IMAGE AND SHOW DESINATED OBJECTS
	'''
	import numpy as np
	import matplotlib
	import matplotlib.pyplot as plt
	from astropy.io import fits
	from matplotlib.colors import LogNorm
	from matplotlib.patches import Circle
	from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
	from astropy.visualization import ZScaleInterval, LinearStretch
	from astropy.wcs import WCS

	outname		= inim[:-5]+'.png'
	data, hdr	= fits.getdata(inim, header=True)
	#fig, ax		= plt.subplots(1)
	wcs			= WCS(hdr)
	norm_zscale	= ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
	fig			= plt.figure() 
	ax			= plt.subplot(projection=wcs)
	im			= ax.imshow(data, cmap='gray', origin='lower', norm=norm_zscale) 

	'''
	ax.set_aspect('equal')
	ax.imshow(data, cmap='gray', origin='lower', norm=LogNorm())
	'''
	for xx, yy in zip(xim_list, yim_list):
		circ = Circle((xx, yy), 15, color='gold', fill=None, linewidth='0.3')
		ax.add_patch(circ)
	for i, txt in enumerate(numb_list):
		xim		= xim_list[i]
		yim		= yim_list[i]
		ax.text(xim+7.5, yim+7.5, str(txt), color='gold', fontsize=5)
	if add != None:
		for xx, yy in zip(xim_addlist, yim_addlist):
			'''
			circ = Circle((xx, yy), 15, color='tomato', fill=None, linewidth='0.5')
			ax.add_patch(circ)
			'''
			ax.scatter(xx, yy, color='tomato', marker='o', alpha=0.1, s=10)
		for i, txt in enumerate(numb_addlist):
			xim		= xim_addlist[i]
			yim		= yim_addlist[i]
			ax.text(xim+7.5, yim+7.5, str(txt), color='tomato', fontsize=5)
	else:
		pass
	fig.savefig(outname, dpi=500, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format=None,
		transparent=False, bbox_inches=None, pad_inches=0.1,
		frameon=None, metadata=None)
	comment		= outname+'\tis generated.';print(comment)
	plt.close()
#-------------------------------------------------------------------------#
def sedualcom(inim, gain, pixscale, det_sigma=1.5, backsize=str(64), backfiltersize=str(3), detect='detection.fits'):
	"""
	SourceEXtractor
	APERTURE    3", 5", 7",
				1.0seeing, 1.2seeing ,1.5seeing ,1.7seeing ,2.0seeeing
	INPUT   :   (image).fits
				aperture    []
				seeing_fwhm [pixel]
	OUTPUT  :   no return
				.cat
	"""
	import numpy as np
	import os
	from astropy.io import ascii
	#   FILE CHECK
	#	CONFIG FILES (USER BASE PATH)
	#sharepath       = '/data7/changsu/code/'
	configfile      = '/home/sonic/Research/yourpy/config/targetphot.sex'
	paramfile       = '/home/sonic/Research/yourpy/config/targetphot.param'
	nnwfile		    = '/home/sonic/Research/yourpy/config/targetphot.nnw'
	convfile	    = '/home/sonic/Research/yourpy/config/targetphot.conv'

	#configfile      = '/d0/code/sex.config/targetphot.sex'
	#paramfile       = '/d0/code/sex.config/targetphot.param'
	#nnwfile		    = '/d0/code/sex.config/targetphot.nnw'
	#convfile	    = '/d0/code/sex.config/targetphot.conv'
	try:
		comment = '\nSourceEXtractor START\n' \
				+ 'IMAGE\t\t: '+inim+'\n' \
				+ 'GAIN\t\t: '+str(gain)+'\n' \
				+ 'PIXSCALE\t\t: '+str(pixscale)+'\n' \
				+ 'DETECTION SIGMA\t: '+str(det_sigma)+'\n' \
				+ 'PARAM\t\t: '+paramfile+'\n' \
				+ 'BACKSIZE\t: '+backsize+'\n' \
				+ 'BACKFILTER\t: '+backfiltersize+'\n' \
				+ 'CONFIG\t\t: '+configfile+'\n' \
				+ 'NNW\t\t: '+nnwfile+'\n' \
				+ 'CONVOLVE\t: '+convfile
		print(comment)
	except:
		comment = 'CHECK configfile/paramfile/nnewfile/convfile or others.'
		print(comment)
	#   FILE NAME
	oriim   = inim[2:]
	cat     = inim[:-5]+'.cat'
	seg     = inim[:-5]+'.seg.fits'
	bkg     = inim[:-5]+'.bkg.fits'
	sub     = inim[:-5]+'.sub.fits'
	psf     = inim[2:-5]+'.psf'
	aper    = inim[:-5]+'.aper.fits'

	#   BASIC INFO.
	det_area        = 5
	#det_thresh      = det_sigma/np.sqrt(det_area)
	det_thresh      = det_sigma
	detecminarea    = str(det_area)
	detectthresh    = str(det_thresh)
	seeing, fwhm_arcsec = psfex(oriim, pixscale)
	#pixscale        = pixscalecalc(inim)
	#   OPTION
	aperture        = '%.2f'%(3./pixscale)+','+'%.2f'%(5./pixscale)+','+'%.2f'%(7./pixscale)+','+'%.2f'%(seeing)+','+'%.2f'%(1.2*seeing)+','+'%.2f'%(1.5*seeing)+','+'%.2f'%(1.7*seeing)+','+'%.2f'%(2.0*seeing)

	option1 = ' -CATALOG_NAME '+cat+' -PARAMETERS_NAME '+paramfile
	option2 = ' -DETECT_MINAREA '+detecminarea+' -DETECT_THRESH '+detectthresh \
			+' -FILTER Y '+'-FILTER_NAME '+convfile
	option3 = ' -PHOT_APERTURES '+aperture \
			+' -GAIN '+'%.1f'%(gain)+' -PIXEL_SCALE '+'%.1f'%(pixscale)
	option4 = ' -SEEING_FWHM '+'%.3f'%(seeing)+' -STARNNW_NAME '+nnwfile
	option5 = ' -BACK_SIZE '+ backsize \
			+ ' -BACK_FILTERSIZE '+ backfiltersize+' -BACKPHOTO_TYPE LOCAL'
	option6 = ' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND'
	option7 = ' -CHECKIMAGE_NAME '+seg+','+aper+','+','+bkg+','+sub
	option8 = ' -PSF_NAME '+psf
	#   COMMAND
	#   detect = detection.fits is fine image show target significantly and have good seeing and many stars i
	dualcom ='sex -c '+configfile+' '+detect+' , '+inim+' -CATALOG_NAME '+cat+' -PARAMETERS_NAME '+paramfile+ ' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8
	sglcom  ='sex -c '+configfile+' '+inim+' '+option1+' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8

	os.system(dualcom)
	secat   = ascii.read(cat)
	return secat, cat, seeing, fwhm_arcsec

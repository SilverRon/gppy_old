#    SELECT GW HOST GALAXY CANDIDATES
#   MADE BY Gregory S.H. Paek   (2019.02.10)
#============================================================#
import os, glob, sys
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import special
#============================================================#
#    FUNCTION
#------------------------------------------------------------#
def sqsum(a, b):
	'''
	SQUARE SUM
	USEFUL TO CALC. ERROR
	'''
	import numpy as np
	return np.sqrt(a**2.+b**2.)
#------------------------------------------------------------#
def calc_MLR(s, colour, z):
	import numpy as np
	MLR        = 10.**(s*colour+z)
	return MLR
#------------------------------------------------------------#
def calc_MLR_error(s, colour, z, band0_err, band1_err, MLR_err):
	import numpy as np
	term0        = np.log(10)*s*(10**(s*colour+z))*band0_err
	term1        = np.log(10)*(-s)*(10**(s*colour+z))*band1_err
	MLR_err0    = sqsum(term0, term1)
	total_err    = sqsum(MLR_err0, MLR_err)
	total_err.name    = 'MLR_err'
	return total_err
#------------------------------------------------------------#
def calc_mass_bol(M_abs_gal, MLR):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    [W]
	M_abs_sun_bol    = 4.74                            [mag]
	mass_gal    = MLR*(10.**(-0.4*(M_abs_gal-M_abs_sun_bol)))
	#mass_gal.name    = 'stellar mass'
	print('RETURN STELLAR MASS [SOLAR MASS]\n'+'1.\t[kg]\n'+'2.\t[solar mass]\n'+'3.\t[log10([solar mass])]')
	return mass_gal, mass_gal/mass_sun, np.log10(mass_gal/mass_sun)
#------------------------------------------------------------#
def calc_mass(M_abs_gal, MLR, L_sun, M_abs_sun):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    #    [W]
	MLR_sun            = L_sun/L_sun_bol
	mass_gal        = MLR*MLR_sun*(10.**(-0.4*(M_abs_gal-M_abs_sun)))
	print('RETURN STELLAR MASS [SOLAR MASS]\n'+'1.\t[solar mass]\n'+'2.\t[log10([solar mass])]')
	return mass_gal, np.log10(mass_gal)
#------------------------------------------------------------#
'''
def calc_mass_error(M_abs_gal, MLR, L_sun, M_abs_sun, M_abs_gal_err, MLR_err):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    [W]
	MLR_sun            = L_sun/L_sun_bol
	term0            = MLR_sun*(10.**(-0.4*(M_abs_gal-M_abs_sun)))
	term1            = sqsum(MLR_err, MLR*(-0.4)*np.log(10)*M_abs_gal_err)
	total_err        = term0*term1
	total_err.name    = 'mass_err'
	return total_err, np.log10(total_err)
'''
#------------------------------------------------------------#
def calc_mass_error(M_abs_gal, MLR, L_sun, M_abs_sun, M_abs_gal_err, MLR_err):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    #    [W]
	MLR_sun            = L_sun/L_sun_bol
	total_err        = MLR_err*MLR_sun*(10.**(-0.4*(M_abs_gal-M_abs_sun)))
	total_err.name    = 'mass_err'
	return total_err, np.log10(total_err)
#------------------------------------------------------------#
def calc_lumi(M_abs):
	'''
	CALCULATE LUMINOSITY FOR DESINATED BAND WITH B BAND.
	'''
	import numpy as np
	L_sun_B        = 1.9*10.**26                    #    [W]
	M_abs_B        = 5.497                            #    [AB]
	L            = L_sun_B*(10.**((M_abs_B-M_abs)/2.5))
	return L
#------------------------------------------------------------#
def calc_mass_K(M_abs_gal, M_abs_sun):
	import numpy as np
	mass_gal    = 10.**((M_abs_sun-M_abs_gal)/2.5)
	return mass_gal
#------------------------------------------------------------#
def prob_fin(map_prob, P_2D, D_gal, D_gal_sig, D_gw, D_gw_sig):
	'''
	P_gal    =    (normalization factor)    *    (Mappeli+18 prob) *
				(P_2D from bayester)    *    (1-error fuction for distance)
	'''
	import numpy as np
	from scipy.special import erf
	term    = (1.-erf((np.abs(D_gal-D_gw))/(sqsum(D_gal_sig, D_gw_sig)**2)))
	P_gal    = map_prob * P_2D * term
	nP_gal    = (P_gal)/np.sum(P_gal)
	return nP_gal
#------------------------------------------------------------#
#   FUNCTION
def replaceNULL(x):
	if x == 'null':
		return -99
	return x
#--------------------------------------------------------------------------------#
def matching(inra, inde, refra, refde):
	"""
	MATCHING TWO CATALOG WITH RA, Dec COORD. WITH python
	INPUT   :   SE catalog, SDSS catalog file name, sepertation [arcsec]
	OUTPUT  :   MATCED CATALOG FILE & TABLE
	"""
	#	MODULE
	import numpy as np
	import astropy.units as u
	from astropy.table import Table, Column
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii
	#comment		= '\nMATCHING START\t'+str(len(inra))+', '+str(len(refra));print(comment)
	#	SKY COORD.
	coo_intbl   = SkyCoord(inra, inde, unit=(u.deg, u.deg))
	coo_reftbl  = SkyCoord(refra, refde, unit=(u.deg, u.deg))
	#   INDEX FOR REF.TABLE
	indx, d2d, d3d  = coo_intbl.match_to_catalog_sky(coo_reftbl, nthneighbor=1)
	comment		= 'MATHCING END';print(comment)
	return indx, d2d
#--------------------------------------------------------------------------------#
def probcut(healpixfits, prob_cut=0.9):
	import healpy as hp
	import numpy as np
	from astropy.table import Table, Column, MaskedColumn
	comment     = '='*80+'\n' \
				+ 'Reading healpix file ...'
	print(comment)
	hpx, hdr    = hp.read_map(healpixfits, h=True, verbose=False)
	hdr         = dict(hdr)
	gwdist      = hdr['DISTMEAN']
	gwdiststd   = hdr['DISTSTD']
	#   h       : GIVE ME HEADER INFO
	#   verbose : DON'T PRINT INFO MESSAGE
	#   AFTER TYPING '%matplotlib' ON THE PYTHON PROMPT.
	#   SHOW MAP
	#   hp.mollview(hpx)
	npix        = hdr['NAXIS2']         # i.e. 50331648
	nside       = hp.npix2nside(npix)   # i.e. 2048
	ordering    = hdr['ORDERING']
	#   SORTING
	#   i.e.
	#   array([  1.72849272e-121,   1.72849272e-121,   1.72849272e-121, ...,
	#           6.41252562e-005,   6.41299517e-005,   6.41368262e-005])
	hpx_sort    = (-1)*np.sort((-1)*hpx)
	indx_sort   = np.argsort((-1)*hpx)
	stack       = np.cumsum(hpx_sort)
	#   MAKE PROB. CUT INDEX (SMALLER THAN PROB.) FOR SORTED DATA
	#   prob_cut    = 0.9
	indx_prob   = np.where(stack <= prob_cut)
	ipix_nparr  = np.linspace(0, npix, npix, dtype=int)
	#   NEED TO CONVERT NUMPY ARRAY TO LIST (WHY?)
	ipix        = ipix_nparr.tolist()
	#   PIXEL COORD. TO SPHERERICAL COORD.
	if ordering == 'NESTED':
		theta, phi  = hp.pix2ang(nside, ipix, nest=False)
	else :
		theta, phi  = hp.pix2ang(nside, ipix, nest=True)
	#   RADIAN TO DEGREE(RA & DEC FORM)
	ra, de      = np.rad2deg(phi), np.rad2deg(0.5*np.pi-theta)
	#   SORTING TO DECREASING ORDER
	rasort      = ra[indx_sort]
	desort      = de[indx_sort]
	#   EXPORT WITHIN PROB.
	racut       = rasort[indx_prob]
	decut       = desort[indx_prob]
	#   MAKE ASCII FILE: ra[deg], dec[deg]
	#   PUT # in front of column name because of running stilts
	outdata		= Table([racut, decut, hpx_sort[indx_prob]], \
						names=['ra', 'dec', 'prob'])	
	#   CALCULATE MEAN DISTANCE BETW. PIXELS
	alldeg		= 41253.         						#	[deg^2]
	mpix_dist   = np.sqrt(alldeg / npix)*3600.			#	["/pix]
	
	return outdata, mpix_dist*2., gwdist, gwdiststd
#--------------------------------------------------------------------------------#
def generate_target(result, dist_lower, dist_upper):
	from astropy.table import Table, Column
	from astropy.io import ascii
	indx_dist   = np.where(((dist_lower < result['dist']) & (dist_upper > result['dist'])))
	num_within  = len(indx_dist[0])
	comment     = 'Cadidates galaxies\t: '+str(num_within)
	print(comment)
	result_in   = result[indx_dist]
	#   CONVERT 'masked' to -99
	result_in   = result_in.filled(-99)
	ra_result   = result_in['ra']
	dec_result  = result_in['dec']
	name        = []
	ra          = []
	dec         = []
	prob        = []
	dist        = []
	disterr     = []
	Bmag        = []
	Bmagerr     = []
	Kmag        = []
	Kmagerr     = []
	W1          = []
	W1err       = []
	for i in range(len(result_in['GWGCname'])):
		if result_in['GWGCname'][i] != '-99':
			name.append(result_in['GWGCname'][i])
		else:
			if result_in['PGC'][i]  != '-99':
				name.append(result_in['PGC'][i])
			else:
				if result_in['HyperLEDAname'][i]    != '-99':
					name.append(result_in['HyperLEDAname'][i])
				else:
					if result_in['2MASSname'][i]    != '-99':
						name.append('2MASS+'+result_in['2MASSname'][i])
					else:
						name.append('SDSS+'+result_in['SDSS-DR12name'][i])
		coord_new   = SkyCoord(ra=ra_result[i]*u.degree, \
								dec=dec_result[i]*u.degree, frame='icrs')
		#   RA [deg] to RA [hh:mm:ss.s]
		ra_h        = int(coord_new.ra.hms[0])
		ra_m        = int(coord_new.ra.hms[1])
		ra_s        = int(coord_new.ra.hms[2])
		if ra_h < 10: ra_h_str      = '0'+str(ra_h)
		else: ra_h_str              = str(ra_h)     
		if ra_m < 10: ra_m_str      = '0'+str(ra_m)
		else: ra_m_str              = str(ra_m)     
		if ra_s < 10: ra_s_str      = '0'+str(ra_s)
		else: ra_s_str              = str(ra_s)     
		ra_new      = ra_h_str+':'+ra_m_str+':'+ra_s_str
		ra.append(ra_new)

		de_d        = int(coord_new.dec.dms[0])
		de_m        = np.abs(int(coord_new.dec.dms[1]))
		de_s        = round(np.abs(coord_new.dec.dms[2]), 1)
		if de_d > 0 :
			if de_d < 10:
				de_d_str            = '+0'+str(de_d)
			else:
				de_d_str            = '+'+str(de_d)
		else:
			de_d_str                = str(de_d)
		if de_m < 10: de_m_str      = '0'+str(de_m)
		else: de_m_str              = str(de_m)     
		if de_s < 10: de_s_str      = '0'+str(de_s)
		else: de_s_str              = str(de_s)     
		de_new      = de_d_str+':'+de_m_str+':'+de_s_str
		dec.append(de_new)

	name        = np.array(name)
	ra          = np.array(ra)
	dec         = np.array(dec)
	data        = Table([name, ra, dec], \
							names=['#name', 'ra', 'dec'])
	outdat      = 'gw_candidates_qv.dat'
	ascii.write(data, outdat, format='fixed_width', delimiter=' ')
	return result_in
#------------------------------------------------------------#
def addwise(result_in, raname='ra_1', dename='dec_1'):
	import numpy as np
	from astropy.table import Table, Column
	from astropy.io import ascii
	ra      = result_in[raname]
	dec     = result_in[dename]
	w1      = []
	w1sig   = []

	for i in range(len(ra)):
		table_  = Irsa.query_region(SkyCoord(ra[i], dec[i], \
									unit=(u.deg,u.deg), frame='icrs'), \
									catalog='allwise_p3as_psd', radius='0d1m0s')
		wra,wdec= table_['ra'], table_['dec']
		dist    = np.sqrt( (ra[i] - wra)**2. + (dec[i] -wdec)**2. )
		indx_min= np.where( dist == np.min(dist) )
		table   = table_[indx_min]
		try:
			print('add wise info ...')
			w1.append(table['w1mpro'][0])
			w1sig.append(table['w1sigmpro'][0])
		except:
			print('exception!')
			w1.append(-99)
			w1sig.append(-99)
	
	w1      = np.array(w1)
	w1sig   = np.array(w1sig)
	result_in['w1']     = w1
	result_in['w1sig']  = w1sig
	new_result          = result_in
	return new_result
#------------------------------------------------------------#
def app2abs(mag_app, dist, appkey='K', distkey='dist'):
	#appmag  = intbl[appkey]
	#dist    = intbl[distkey]
	import numpy as np
	abslist = []
	for i in range(len(mag_app)):
		if (mag_app[i] != -99) & (dist[i] != -99):
			mag_abs  = mag_app[i]-5.*np.log10(dist[i]*1e6)+5
		else:
			mag_abs = -99
		abslist.append(round(mag_abs, 3))
	return np.array(abslist)
#------------------------------------------------------------#
def prob(gwdist, gwdister, gldist, gldister, appmag, prob):
	'''
	Swope Supernova Survey 2017a (SSS17a),
	the Optical Counterpart to a Gravitational Wave Sourve
	'''
	import numpy as np
	from scipy import special
	P_gal   = []
	for i in range(len(prob)):
		if gldister[i] == 'null':
			gldister[i] = float(0)
		else:
			gldister[i]	= float(gldister[i])

		L   = ((gldist[i])**2.)*10.**(-0.4*appmag[i])
		braket  = np.abs(gldist[i]-gwdist)/((gwdister**2.)+(gldister[i]**2.))
		P       = L*prob[i]*(1.-special.erf(braket))

		P_gal.append(P)
	P_gal   = np.array(P_gal)
	nP_gal  = P_gal*(1./np.sum(P_gal))
	return nP_gal
#------------------------------------------------------------#
def generate_final_list(incat, outdat):
	import numpy as np
	from astropy.table import Table, Column
	from astropy.io import ascii
	from astropy.coordinates import SkyCoord
	import astropy.units as u
	num         = len(incat)
	try:
		indx_TheMost= np.where(incat['K_prob'] == np.max(incat['K_prob']))
	except:
		indx_TheMost= np.where(incat['B_prob'] == np.max(incat['B_prob']))
	comment     = '='*80+'\n' \
				+ 'Generating '+outdat+'\n' \
				+ 'Cadidates galaxies\t: '  +str(num)+'\n'+'-'*80+'\n' \
				+ 'The Highest Prob.\t: '   +str(incat['GWGCname'][indx_TheMost][0])
	print(comment)

	ra_result   = incat['ra']
	dec_result  = incat['dec']
	name        = []
	ra          = []
	dec         = []
	dist        = []
	Bmag        = []
	Kmag        = []

	for i in range(len(incat['GWGCname'])):
		if incat['GWGCname'][i] != '-99':
			name.append(incat['GWGCname'][i])
		else:
			if incat['PGC'][i]  != -99:
				name.append(incat['PGC'][i])
			else:
				if incat['HyperLEDAname'][i]    != '-99':
					name.append(incat['HyperLEDAname'][i])
				else:
					if incat['2MASSname'][i]    != '-99':
						name.append('2MASS+'+incat['2MASSname'][i])
					else:
						name.append('SDSS+'+incat['SDSS-DR12name'][i])
		coord_new   = SkyCoord(ra=ra_result[i]*u.degree, \
								dec=dec_result[i]*u.degree, frame='icrs')
		#   RA [deg] to RA [hh:mm:ss.s]
		ra_h        = int(coord_new.ra.hms[0])
		ra_m        = int(coord_new.ra.hms[1])
		ra_s        = int(coord_new.ra.hms[2])
		if ra_h < 10: ra_h_str      = '0'+str(ra_h)
		else: ra_h_str              = str(ra_h)     
		if ra_m < 10: ra_m_str      = '0'+str(ra_m)
		else: ra_m_str              = str(ra_m)     
		if ra_s < 10: ra_s_str      = '0'+str(ra_s)
		else: ra_s_str              = str(ra_s)     
		ra_new      = ra_h_str+':'+ra_m_str+':'+ra_s_str
		ra.append(ra_new)

		de_d        = int(coord_new.dec.dms[0])
		de_m        = np.abs(int(coord_new.dec.dms[1]))
		de_s        = round(np.abs(coord_new.dec.dms[2]), 1)
		if de_d > 0 :
			if de_d < 10:
				de_d_str            = '+0'+str(de_d)
			else:
				de_d_str            = '+'+str(de_d)
		else:
			de_d_str                = str(de_d)
		if de_m < 10: de_m_str      = '0'+str(de_m)
		else: de_m_str              = str(de_m)     
		if de_s < 10: de_s_str      = '0'+str(de_s)
		else: de_s_str              = str(de_s)     
		de_new      = de_d_str+':'+de_m_str+':'+de_s_str
		dec.append(de_new)
		dist.append(round(incat['dist'][i], 3))
		Bmag.append(round(incat['B'][i], 3))
		Kmag.append(round(incat['K'][i], 3))
	name        = np.array(name)
	ra          = np.array(ra)
	dec         = np.array(dec)
	dist        = np.array(dist)
	Bmag        = np.array(Bmag)
	Kmag        = np.array(Kmag)
	data        = Table([name, ra, dec, dist, Bmag, Kmag], names=['name', 'ra', 'dec', 'dist', 'Bmag', 'Kmag'])
	# outdat      = 'gw_candidates_sorted.dat'
	ascii.write(data, outdat, format='fixed_width', delimiter=' ')
	outcat      = ascii.read(outdat)
	return outcat
#------------------------------------------------------------#
def heal2table(healpixfits, confidence=0.9, eventname='GW_signal', save_path='./', hdr=True, view=True):
	import numpy as np
	import healpy as hp
	from astropy.table import Table
	hpx, hdr    = hp.read_map(healpixfits, h=True, verbose=False)
	hdr         = dict(hdr)

	gwdist      = hdr['DISTMEAN']
	gwdiststd   = hdr['DISTSTD']
	#   h       : GIVE ME HEADER INFO
	#   verbose : DON'T PRINT INFO MESSAGE
	#   AFTER TYPING '%matplotlib' ON THE PYTHON PROMPT.
	npix        = hdr['NAXIS2']         # i.e. 50331648
	nside       = hp.npix2nside(npix)   # i.e. 2048
	ordering    = hdr['ORDERING']
	#   SORTING
	#   i.e.
	#   array([  1.72849272e-121,   1.72849272e-121,   1.72849272e-121, ...,
	#           6.41252562e-005,   6.41299517e-005,   6.41368262e-005])
	hpx_sort    = (-1)*np.sort((-1)*hpx)
	indx_sort   = np.argsort((-1)*hpx)
	stack       = np.cumsum(hpx_sort)
	#   MAKE PROB. CUT INDEX (SMALLER THAN PROB.) FOR SORTED DATA
	#   confidence    = 0.9
	indx_prob   = np.where(stack <= confidence)
	ipix_nparr  = np.linspace(0, npix, npix, dtype=int)
	#   NEED TO CONVERT NUMPY ARRAY TO LIST (WHY?)
	ipix        = ipix_nparr.tolist()
	#   PIXEL COORD. TO SPHERERICAL COORD.
	if ordering == 'NESTED':
		theta, phi  = hp.pix2ang(nside, ipix, nest=False)
	else :
		theta, phi  = hp.pix2ang(nside, ipix, nest=True)
	#   RADIAN TO DEGREE(RA & DEC FORM)
	ra, de      = np.rad2deg(phi), np.rad2deg(0.5*np.pi-theta)
	#   SORTING TO DECREASING ORDER
	rasort      = ra[indx_sort]
	desort      = de[indx_sort]
	#   EXPORT WITHIN PROB.
	racut       = rasort[indx_prob]
	decut       = desort[indx_prob]
	#   MAKE ASCII FILE: ra[deg], dec[deg]
	#   PUT # in front of column name because of running stilts
	healtbl		= Table([racut, decut, hpx_sort[indx_prob]], names=['ra', 'dec', 'P_2D'])
	#   SHOW MAP
	if view == True:
		hp.mollview(hpx, title=eventname)
		plt.savefig(save_path+'/'+eventname+'-skymap.png')
	return healtbl, hdr
#------------------------------------------------------------
def kilonova_mag(gwdist, gwdiststd):
	import numpy as np
	m0		= 17.476	#	[AB] in i-band (t0+10h)
	m0err	= 0.018	
	dist0	= 38.4		#	[MPC]	Im et al. 2017
	dist0err= 8.9
	
	m		= m0+5.*np.log10(gwdist/dist0)
	merr	= np.sqrt( (m0err)**2 + ((5.*gwdiststd)/(gwdist*np.log(10)))**2 + ((5.*dist0err)/(dist0*np.log(10)))**2 )
	return m, merr
#------------------------------------------------------------
	
	
	
	

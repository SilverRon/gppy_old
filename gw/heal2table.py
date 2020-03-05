#	HEALPIX PIXEL TO RA DEC PROB. FROM *.fits.gz
#------------------------------------------------------------
#	2020.02.19	CREATED BY Gregory S.H. Paek
#============================================================
def heal2table(healpixfits, confidence=0.9, eventname='GW_signal', save_path='./', hdr=True, view=True):
	import matplotlib.pyplot as plt
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



eventname = 'S200213t-Update'
param_heal = dict(	healpixfits = 'S200213t-Update.fits.gz',
					confidence = 0.1,
					eventname = eventname,
					save_path = '.',
					hdr = True,
					view = True)

healtbl, hdr = heal2table(**param_heal)
healtbl.write(eventname+'_healpix.dat', format='ascii', overwrite=True)

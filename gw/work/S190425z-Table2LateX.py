#	S190425z OBSERVATION TABLE -> LaTeX FORMAT
#   2019.09.01	MADE BY		Gregory S.H. Paek
#============================================================
#	MODULE
#------------------------------------------------------------
import healpy as hp
import numpy as np
import time
import os, glob, sys
from astropy.table import Table, Column, MaskedColumn, vstack, hstack
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time
import matplotlib.pyplot as plt
from imsng import gw
from imsng import tool
import ligo.skymap.plot
from scipy.stats import norm
import scipy.stats
import warnings
#============================================================
path_table = '/data1/S190425z/1.result/Update/phot*sources*.dat'
path_candi = '/data1/S190425z/info/Update/S190425z_Update-all_candi.txt'
path_save = '/data1/S190425z/1.result/table'
reftbl = ascii.read(path_candi)
t0 = Time('2019-04-25T08:18:05.017', format='isot', scale='utc')
#------------------------------------------------------------
tblist = []
for incat in glob.glob(path_table):
	if 'GECKO' in incat:
		pass
	else:
		intbl = ascii.read(incat)
		intbl = intbl[intbl['sources']!='None']
		intbl = intbl[np.argsort(intbl['jd'])]
		for i, objs in enumerate(intbl['sources']):
			obs = intbl['obs'][i]
			ra, dec = intbl['ra'][i], intbl['dec'][i]
			c = SkyCoord(ra, dec, frame='icrs', unit='deg')
			cstr = c.to_string('hmsdms')
			rahms, dedms = cstr.split(' ')[0], cstr.split(' ')[1]
			ul = round(intbl['ul'][i], 2)
			for obj in objs.split(','):
				if (('SDSS' not in obj) & ('PGC' not in obj) & (len(obj)>15)):
					obj = '2MASXJ'+obj[6:]
				tmptbl = Table(	[[obj], [rahms], [dedms], [obs], [ul]],
								names=('Name', 'R.A.', 'Dec.', 'Telescope', 'Limiting Mag'))
				tblist.append(tmptbl)
comtbl = vstack(tblist)
comtbl.write(path_save+'/candi.dat', format='ascii.aastex', overwrite=True)
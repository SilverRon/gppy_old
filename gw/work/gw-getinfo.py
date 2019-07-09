#   GET INFORMATION FROM GW BAYESTAR FITS FILES
#   BASED ON https://emfollow.docs.ligo.org/userguide/tutorial/index.html
#   2019.03.03 MADE		BY Gregory S.H. Paek
#	2019.03.28 UPDATED	BY Gregory S.H. Paek
#	2019.04.18 UPDATED	BY Gregory S.H. Paek
#	2019.05.02 UPDATED	BY Gregory S.H. Paek (ver. 2.0)
#	2019.05.10 UPDATED	BY Gregory S.H. Paek
#	2019.07.08 UPDATED	BY Gregory S.H. Paek
#============================================================#
#	MODULE
#------------------------------------------------------------#
import healpy as hp
import numpy as np
import time
import os, glob
from astropy.table import Table, Column, MaskedColumn
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import matplotlib.pyplot as plt
from imsng import gw
from imsng import tool
from astropy.table import Table, vstack
#------------------------------------------------------------#
#	FUNCTION
#------------------------------------------------------------#
def getgwinfo(healpixfits, confidence, save_path='./'):
	eventname, role	= os.path.basename(healpixfits).split('-')[0], os.path.basename(healpixfits).split('-')[1]
	comment     = 'PROCESSING {}'.format(healpixfits); print(comment)
	param_heal	= dict(	healpixfits=healpixfits,
						confidence=confidence,
						eventname=eventname,
						save_path=save_path,
						hdr=True,
						view=False)
	healtbl, hdr= gw.heal2table(**param_heal)
	gwdist      = hdr['DISTMEAN']
	gwdiststd   = hdr['DISTSTD']
	npix        = hdr['NAXIS2']         # i.e. 50331648
	alldeg		= 41253
	skydeg_cut	= (alldeg/npix)*len(healtbl)
	try:
		detector	= hdr['INSTRUME']
	except:
		detector	= 'UNKNOWN'
	onetbl		= Table(	[[eventname], [role], [gwdist], [gwdiststd], [skydeg_cut], [hdr['DATE-OBS']], [hdr['MJD-OBS']], [detector]],
							names=('event', 'role', 'dist', 'diststd', 'cr', 'date-obs', 'mjd', 'detector'))
	return onetbl
#------------------------------------------------------------#
#	WORK
#------------------------------------------------------------#
confidence	= 0.9
tblist		= []
outname		= 'O3run_summary.dat'
for healpixfits in glob.glob('./*.fits.gz'):
	tblist.append(getgwinfo(healpixfits=healpixfits, confidence=confidence))
fintbl		= vstack(tblist)
if outname in glob.glob(outname): os.system('mv {0} {1}'.format(outname, outname+'.bkg'))
fintbl.write(outname, format='ascii', overwrite=True)
#   UTILIZE GW P_3D LOCALIZATION MAP
#   BASED ON https://github.com/lpsinger/gw-galaxies/blob/master/gw-galaxies.ipynb
#   2019.07.19 MADE		BY Gregory S.H. Paek
#	2019.XX.XX UPDATED	BY Gregory S.H. Paek
#============================================================
#	MODULE
#------------------------------------------------------------
import astropy.utils.data
from matplotlib import pyplot as plt
import numpy as np
import healpy as hp
#------------------------------------------------------------
from astropy.table import Table, vstack, hstack, Column
import astropy.units as u
from astropy.coordinates import SkyCoord
import ligo.skymap.plot
from scipy.stats import norm
import scipy.stats
#============================================================
url		= 'https://emfollow.docs.ligo.org/userguide/_static/bayestar.fits.gz'
filename= astropy.utils.data.download_file(url)
#------------------------------------------------------------
hpx, hdu= hp.read_map(filename, verbose=True, h=True)
# hp.mollview(prob)
npix	= len(hpx)
nside	= hp.npix2nside(npix)
ipix	= np.arange(0, npix, step=1)
ra, dec	= hp.pix2ang(nside, ipix, lonlat=True)
# hp.ang2pix(nside, ra, dec, lonlat=True)		# REVERSE
ipix_max= np.argmax(prob)
ra_max, dec_max		= hp.pix2ang(nside, ipix_max, lonlat=True)
#------------------------------------------------------------
from ligo.skymap.postprocess import find_greedy_credible_levels
credible_levels		= find_greedy_credible_levels(hpx)
# credible_levels[ipix]
# credible_levels[ipix] <= 0.9
for conf in [0.5, 0.9]:
	lcsize	= np.sum(credible_levels <= conf) * hp.nside2pixarea(nside, degrees=True)
	print(lcsize)

hpx50	= hpx[credible_levels <= 0.5]
hpx90	= hpx[credible_levels <= 0.9]


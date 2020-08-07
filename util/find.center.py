#	FIND CENTER RA DEC
#	CREATED	2020.06.11	Gregory S.H. Paek
#============================================================
import os, glob, subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from imsng import phot
from imsng import query
from imsng import tool
from imsng import gcurve
import time
from datetime import date
#============================================================
def center_pixel2world(inim):
	# inim = 'Calib-LSGT-NGC4321-20200104-174639-r-180.fits'
	hdr = fits.getheader(inim)
	w = WCS(inim)
	xim, yim = hdr['NAXIS1'], hdr['NAXIS2']
	xcent, ycent = xim/2, yim/2
	racent, decent = w.all_pix2world(xcent, ycent, 1)
	return racent, decent

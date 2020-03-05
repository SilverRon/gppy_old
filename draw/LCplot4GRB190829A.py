#	GRB 190829A UKIRT PHOTOMETRY PLOT
#	CREATED	2020.02.15	Gregory S.H. Paek
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
from astropy import constants as const
from imsng import phot
from scipy.optimize import curve_fit
import time
#============================================================
#	FUNCTION
#============================================================
def linear(x, a, b):
	return a * np.log10(x) + b
#============================================================
#	USER SETTING
#============================================================
#	TARGET
#------------------------------------------------------------
targetname = 'GRB190829A'
t0 = Time('2019-08-29T19:55:53', format='isot', scale='utc')
z = 0.0785
zerr = 0.005
h = 0.7
H0 = h*100					# [km/s/Mpc]
c = const.c.value*1e-3		# [km/s]
d = c*z/H0
derr = c*zerr/H0
path_base	= '.'
photbl = ascii.read(path_base+'/phot.all.dat')
# reftbl = ascii.read(path_base+'/phot.nosub.dat')
#------------------------------------------------------------
plt.close('all')
fig = plt.figure(figsize=(10, 8))
ax = plt.subplot(111)
'''
for band in ['J', 'H', 'K']:
	subtbl = photbl[(photbl['filter']==band)&(photbl['mag']!=-99)]
	ax.errorbar(subtbl['jd']-t0.jd, subtbl['mag'], yerr=subtbl['magerr'],
				marker='o', ms=15,
				capsize=10, capthick=1,
				# color='grey',
				alpha=0.5,
				linestyle='None', label='UKIRT {}-band'.format(band))
'''
#------------------------------------------------------------
#	COMPARE PREVIOUS PHOTOMETRY RESULT FOR EACH FILTERS
#------------------------------------------------------------
# for band in ['J', 'H', 'K']:
# for band in ['J']:
# for band in ['H']:
for band in ['K']:
	subtbl = photbl[(photbl['filter']==band)&(photbl['mag']!=-99)]
	subtbl0 = reftbl[(reftbl['filter']==band)&(reftbl['mag']!=-99)]
	ax.errorbar(subtbl['jd']-t0.jd, subtbl['mag'], yerr=subtbl['magerr'],
				marker='o', ms=15,
				capsize=10, capthick=1,
				# color='dodgerblue',
				# color='orange',
				color='g',
				alpha=0.5,
				linestyle='None', label='UKIRT {} Subt.'.format(band))
	'''
	ax.errorbar(subtbl0['jd']-t0.jd, subtbl0['mag'], yerr=subtbl0['magerr'],
				marker='o', ms=15,
				capsize=10, capthick=1,
				color='grey',
				alpha=0.5,
				linestyle='None', label='UKIRT {} No Subt.'.format(band))
	'''
#------------------------------------------------------------
#	FITTING
#------------------------------------------------------------
'''
for band, color in zip(['K', 'H', 'J'], ['g', 'orange', 'dodgerblue']):
	# band = 'K'
	fitbl = photbl[(photbl['filter'] == band) & (photbl['mag'] != -99)][:3]
	xdata = fitbl['jd']-t0.jd
	ydata = fitbl['mag']
	popt, pcov = curve_fit(linear, xdata, ydata)

	alpha = round(popt[0], 3)
	xdummy = np.arange(xdata[0], 10, 1)
	ax.plot(xdummy, linear(xdummy, *popt), color=color, linestyle='--', label=r'$\alpha$={}'.format(alpha))
'''
#------------------------------------------------------------

# '''
ax.fill_between(x=np.arange(4.5, 5.5, 0.01),
				y1=-99,
				y2=+99,
				facecolor='tomato',
				alpha=0.5,
				label='SN report')

ax.legend(fontsize=20, loc='upper right', framealpha=1.0)
plt.xlabel(r'$\log(t-t_{0})$ [days]', fontsize=20)
plt.ylabel('AB magnitude', fontsize=20)
#------------------------------------------------------------
ax0 = ax.twinx()
ax.tick_params(which='both', direction='in', labelsize=15)
ax.minorticks_on()
ax.set_ylim(21, 16.5)
ax0.tick_params(which='both', direction='in', labelsize=15)
ax0.set_ylim(21-5*np.log10(d*1e6)+5, 16-5*np.log10(d*1e6)+5)
ax0.minorticks_on()
ax0.set_ylabel('Absolute magnitude', fontsize=20, rotation=270, labelpad=20)
plt.xscale('log')
plt.tight_layout()
#------------------------------------------------------------
# plt.savefig(path_base+'/lc.png', dpi=300, overwrite=True)


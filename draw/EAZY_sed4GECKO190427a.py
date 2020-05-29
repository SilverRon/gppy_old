#	DRAW EAZY OUTPUT AS SED
#	CREATED	2020.05.07	Gregory S.H. Paek
#	UPDATE 
#============================================================
import os, glob, subprocess, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import time
from datetime import date
import bisect
#============================================================
#	FUNCTION
#============================================================

#============================================================
#	SETTING
#------------------------------------------------------------
path_base = '/data1/S190425z/KMTNet/fieldsub/transient/gecko190427a/eazy_result'

obstbl = ascii.read('{}/1.obs_sed'.format(path_base))
tmptbl = ascii.read('{}/1.temp_sed'.format(path_base))
pztbl = ascii.read('{}/photz.zout'.format(path_base))

#	Scaling
scalefactor = 1e4
obstbl['flux_cat'] = obstbl['flux_cat']*scalefactor
tmptbl['tempflux'] = tmptbl['tempflux']*scalefactor

labels = ['PS1 g', 'PS1 r', 'PS1 i', 'PS1 z', 'PS1 y', '2MASS J', '2MASS H', 'W1', 'W2', 'W3']
photz, photz_l68, photz_u68, photz_chi = pztbl['z_a'][0].item(), pztbl['l68'][0].item(), pztbl['u68'][0].item(), pztbl['chi_a'][0].item()
#============================================================
#	PLOT
#------------------------------------------------------------
plt.close('all')
fig		= plt.figure()				# Define "figure" instance
fig.set_size_inches(8,6)			# Physical page size in inches, (lx,ly)
# suptit	= r'SED ($z_{phot}=0.16_{0.091}^{0.329}$, $\chi^2=1.274$)'
suptit	= r'SED (z={}({}-{}), $\chi^2={}$)'.format(photz, photz_l68, photz_u68, photz_chi)
# suptit	= r'SED ($z_{phot}=1.279_{1.281}^{1.299}$, $\chi^2=137.7811$)'


xticks = []
for i in [1e3, 1e4]:
	for j in np.arange(1, 10, 1):
		xticks.append(i*j)
xticks.append(1e5)
#------------------------------------------------------------
plt.subplot(3, 1, (1, 2))
plt.title(suptit, fontsize=18)
plt.xscale('symlog')


'''
#	PLOT ALL
params_plot	= dict(	x=obstbl['lambda'],
					# y=obstbl['lambda']*obstbl['flux_cat'],
					y=obstbl['flux_cat'],
					yerr=obstbl['err_full'],
					marker='o', ms=10, mew=1,
					mec='k', mfc='tomato', fmt='ko', alpha=0.75,
					capsize=5, capthick=2,
					label='input')
plt.errorbar(**params_plot)
'''


colorlist = ['green', 'tomato', 'purple', 'violet', 'pink', 'peru', 'orange', 'cyan', 'dodgerblue', 'royalblue']
for i, c in zip(range(len(obstbl)), colorlist):
	params_plot	= dict(	x=obstbl['lambda'][i],
						# y=obstbl['lambda']*obstbl['flux_cat'],
						y=obstbl['flux_cat'][i],
						yerr=obstbl['err_full'][i],
						marker='o', ms=12, mew=1,
						mec='k', mfc=c, fmt='ko', alpha=0.75,
						capsize=5, capthick=2,
						label=labels[i])
	plt.errorbar(**params_plot)




# plt.plot(tmptbl['lambda'], tmptbl['lambda']*tmptbl['tempflux'])
plt.plot(tmptbl['lambda'], tmptbl['tempflux'])
# plt.xlim(1e3, 1e5)
#------------------------------------------------------------
plt.ylabel(r'$F_{\nu}$ $[10^{-4}Jy]$', fontsize=20)
plt.legend(fontsize=12, loc = 'upper left', framealpha=1.0)
plt.tick_params(which='both', direction='in', labelsize=15)
plt.xticks(xticks)
# plt.yticks(np.arange(0, 10, 2))
# plt.yticks(np.arange(0, 90, 10))
plt.grid(linestyle='--', color='grey', alpha=0.5)
plt.minorticks_on()
#------------------------------------------------------------
#	PLOT 2
#------------------------------------------------------------
plt.subplot(313)
plt.xscale('symlog')

xres, yres = [], []

for i, lamb in enumerate(obstbl['lambda']):
	indx = bisect.bisect(tmptbl['lambda'], lamb)
	xres.append(lamb)
	yres.append(obstbl['flux_cat'][i] - tmptbl['tempflux'][indx])
'''
params_plot	= dict(	x=xres,
					y=yres,
					yerr=obstbl['err_full'],
					marker='o', ms=10, mew=1,
					mec='k', mfc='tomato', fmt='ko', alpha=0.75,
					capsize=5, capthick=2,
					label='input')
plt.errorbar(**params_plot)
'''

for i, c in zip(range(len(obstbl)), colorlist):
	params_plot	= dict(	x=xres[i],
						# y=obstbl['lambda']*obstbl['flux_cat'],
						y=yres[i],
						yerr=obstbl['err_full'][i],
						marker='o', ms=12, mew=1,
						mec='k', mfc=c, fmt='ko', alpha=0.75,
						capsize=5, capthick=2,
						label='_nolegend_')
	plt.errorbar(**params_plot)




plt.axhline(y=0, linestyle='--', color='grey', alpha=0.5)
#------------------------------------------------------------
# plt.xlim(1e3, 1e5)
# plt.ylim(-1e-3*scalefactor, +1e-3*scalefactor)
plt.ylim(-1, +1)
# plt.ylim(-10, +10)
plt.xlabel(r'$\lambda(Angstroms)$', fontsize=20)
# plt.ylabel(r'Residual $[10^{4}Jy]$', fontsize=20)
plt.ylabel('Residual', fontsize=20)
plt.tick_params(which='both', direction='in', labelsize=15)
plt.xticks(xticks)
# plt.yticks(np.arange(-0.5, 1.0, 0.5))
plt.grid(linestyle='--', color='grey', alpha=0.5)
plt.minorticks_on()
#------------------------------------------------------------
# plt.tight_layout()
plt.savefig('sed.png', dpi=500, overwrite=True)
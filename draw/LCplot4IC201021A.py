#	LC for IC201021A
#	CREATED	2020.10.27	Gregory S.H. Paek
#============================================================
import os, glob, subprocess, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
path_save = '/home/sonic/Pictures'
path_phot = '/data1/IceCube-201021A/phot.all.dat'
path_plot = '/data1/S190425z/200625.result/obs.color.table.dat'
#------------------------------------------------------------
photbl = ascii.read(path_phot)
pltbl = ascii.read(path_plot)
obslist = ['LOAO', 'DOAO']
#------------------------------------------------------------
ps1, ps1er = 21.22, 0.04
t0 = Time('2020-10-21T06:37:47.48')
photbl['delt'] = photbl['jd'] - t0.jd
#============================================================
#	PLOT
#------------------------------------------------------------
plt.close('all')
plt.rc('font', family='serif')
fig = plt.figure()
x = 1920 / 2 / fig.dpi
y = 1080 / 2 / fig.dpi
fig.set_figwidth(x)
fig.set_figheight(y)
#------------------------------------------------------------
#	GRID
ncols = 3
nrows = 3
grid = GridSpec(nrows, ncols,
				left=0.1, bottom=0.15, right=0.94, top=0.94, wspace=3, hspace=0.1)
ax1 = fig.add_subplot(grid[0:2, 0:3])
ax2 = fig.add_subplot(grid[2:3, 0:3])
#------------------------------------------------------------
#	AXIS 1
#------------------------------------------------------------
ax1.axhline(y=ps1, color='k', linestyle='--')
ax1.fill_between([0, np.max(photbl['delt'])+1], ps1-ps1er, ps1+ps1er, color='silver', label='PS1')
for obs in obslist:
	obstbl = photbl[photbl['obs']==obs]
	obsptbl = pltbl[pltbl['obs']==obs.lower()]
	o_obstbl = obstbl[obstbl['mag']!=-99]
	v_obstbl = obstbl[obstbl['mag']==-99]

	c1, c2, c3 = obsptbl['c1'].item(), obsptbl['c2'].item(), obsptbl['c3'].item()

	if len(o_obstbl)!=0:
		params_plot	= dict(	
							x=o_obstbl['delt'],
							y=o_obstbl['mag'],
							yerr=o_obstbl['magerr'],
							marker='o', ms=20, mew=2,
							mec='k', mfc=c1, fmt='ko', alpha=0.75,
							capsize=10, capthick=2,
							label=obs,
							)
		ax1.errorbar(**params_plot)
	else:
		pass

	if len(v_obstbl)!=0:
		params_plot	= dict(	
							x=v_obstbl['delt'],
							y=v_obstbl['ul_5sig'],
							yerr=0,
							marker='v', ms=20, mew=2,
							mec='k', mfc=c1, fmt='ko', alpha=0.75,
							capsize=10, capthick=2,
							# label=label,
							)
		ax1.errorbar(**params_plot)
	else:
		pass
#------------------------------------------------------------
#	AXIS 2
#------------------------------------------------------------
ax2.axhline(y=-0.5, color='grey', linestyle='--')
ax2.axhline(y=+0.5, color='grey', linestyle='--')

detbl = photbl[photbl['mag']!=-99]
x = detbl['delt']
y = detbl['mag'] - ps1
yer = np.sqrt((detbl['magerr'])**2 + (ps1er)**2)

params_plot = dict(
					x=x,
					y=y,
					yerr=yer,
					marker='o', ms=10, mew=2,
					mec='k', mfc='grey', fmt='ko', alpha=0.75,
					capsize=10, capthick=2,
					# label=label,
					)
ax2.errorbar(**params_plot)
#============================================================
#	Plot setting
#------------------------------------------------------------
#	AXIS 1
#------------------------------------------------------------
ax1.set_xlim(0, np.max(photbl['delt']))
down, up = ax1.set_ylim()
ax1.set_ylim(ps1+1, ps1-1)
ax1.set_title('IceCube-201021A', fontsize=20)
ax1.set_ylabel('AB mag', fontsize=20)
ax1.legend(fontsize=20, loc='upper right', framealpha=1.0)
ax1.tick_params(which='both', direction='in', labelsize=16)
ax1.grid(linestyle='--', color='grey', alpha=0.5)
ax1.minorticks_on()
ax1.tick_params(
				axis='x',
				which='both',
				bottom=False,
				top=False,
				labelbottom=False
				)
#------------------------------------------------------------
#	AXIS 2
#------------------------------------------------------------
ax2.set_xlim(0, np.max(photbl['delt']))
down, up = ax1.set_ylim()
ax2.set_ylim(-1.0, +1.0)

ax2.set_xlabel('Days after event', fontsize=20)
ax2.set_ylabel(r'$\rm \Delta$m', fontsize=20)
ax2.tick_params(which='both', direction='in', labelsize=16)
ax2.minorticks_on()
ax2.grid(linestyle='--', color='grey', alpha=0.5)
#------------------------------------------------------------
fig.savefig("{}/LC_IC201021A.pdf".format(path_save), overwrite=True)
fig.savefig("{}/LC_IC201021A.png".format(path_save), overwrite=True, dpi=500)


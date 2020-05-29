#	KILONOVA LC PLOT
#	CREATED	2020.04.16	Gregory S.H. Paek
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
from imsng import phot, tool
from scipy.optimize import curve_fit
import time
#------------------------------------------------------------
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
#============================================================
#	FUNCTION
#============================================================

#============================================================
#	USER SETTING
#============================================================
path_base = '/home/sonic/Research/gppy/table'
droutbl = ascii.read(path_base+'/phot_gw170817_Drout.dat')
pirotbl = ascii.read(path_base+'/cocoon_Piro+18.dat')

gwdist0, gwdiststd0 = 38.4, 8.9	# [MPC]

#------------------------------------------------------------
gwdist = 400	# [MPC]
filte = 'r'
#------------------------------------------------------------
#	COCOON MODEL
#------------------------------------------------------------
ptbl = pirotbl[pirotbl['filter']==filte]
param_abs2app = dict(mag=ptbl['absmag'],
					magerr=0.,
					gwdist=gwdist*1e6,	#	[MPC] -> [PC]
					gwdiststd=0.)
ptbl['mag'], ptbl['magerr']	= tool.abs2app(**param_abs2app)
# ptbl['mag']= ptbl['mag']-0.7
#------------------------------------------------------------
#	GW170817-like
#------------------------------------------------------------
dtbl = droutbl[droutbl['filter']==filte]; dtbl['delmjd'].sort()
param_calc_app = dict(mag=dtbl['mag'],
					magerr=dtbl['magerr'],
					gwdist0=gwdist0,
					gwdiststd0=gwdiststd0,
					gwdist1=gwdist,
					gwdiststd1=0.)
mag, magerr= tool.calc_app(**param_calc_app)
#------------------------------------------------------------
#	PLOT 1	: TIME - MAG
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 20})
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(10, 10))
#------------------------------------------------------------
#	GW170817-like
ax0.plot(dtbl['delmjd'], mag, color='red', alpha=0.5, label='GW170817-like')
ax0.fill_between(dtbl['delmjd'], mag-magerr, mag+magerr, color='tomato', alpha=0.15, label='_nolegend_')
#	COCOON MODEL (Piro+2018)
ax0.plot(ptbl['delmjd'], ptbl['mag'], color='dodgerblue', alpha=0.75, label='Shock Cooling')
ax0.fill_between(ptbl['delmjd'], ptbl['mag']-ptbl['magerr'], ptbl['mag']+ptbl['magerr'], color='dodgerblue', alpha=0.3, label='_nolegend_')
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
ax0.set_xlabel(r'Time (Days from merger)', fontsize=20)
ax0.set_ylabel(r'Apparent Magnitude', fontsize=20)

yu, yd = int(np.max(mag))+0.5, int(np.min(ptbl['mag']))-0.5
ax0.set_ylim([yu, yd])
ax0.set_xlim([0,2])
plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM discovery')
ax0.legend(loc='lower right', prop={'size':24})
ax0.minorticks_on()
ax0.tick_params(which='both', direction='in', labelsize=20)
#------------------------------------------------------------
ax1 = ax0.twinx()
ax1.set_ylabel(r'Absolute Magnitude', rotation=270, labelpad=20, fontsize=20)
ax1.set_ylim(yu-5*np.log10(gwdist*1e6)+5, yd-5*np.log10(gwdist*1e6)+5)
ax1.minorticks_on()
ax1.tick_params(which='both', direction='in', labelsize=20)
#------------------------------------------------------------
plt.title('Expected Light Curve in {}-band ({}Mpc)'.format(filte, gwdist))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.minorticks_on()
plt.savefig('{}/kilonovaLC_{}band.{}Mpc.png'.format('.', filte, int(gwdist)), dpi=500, overwrite=True)

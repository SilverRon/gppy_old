#============================================================
#	THIS IS FOR COLLECTING ALL OBSERVATION DATA FOR FOLLOW-UP OBSERVATION TO S190425z
#	2019.06.08	MADE BY Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
from astropy.io import ascii
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy.interpolate import interp1d
#============================================================
#	INPUT
#------------------------------------------------------------
path_base	= '/home/gw/Research'
t0			= Time(hdr['DATE-OBS'], format='isot', scale='utc')
jd0_170817	= 2457983.02852
droutbl		= ascii.read(path_base+'/phot_gw170817_Drout.dat')
pirotbl		= ascii.read(path_base+'/cocoon_Piro+18.dat')
#------------------------------------------------------------
#	COCOON MODEL
prtbl		= pirotbl[pirotbl['filter']=='r']
prtbl['mag'], prtbl['magerr']	=	tool.abs2app(prtbl['absmag'], 0, gwdist*1e6, gwdiststd*1e6)
prtbl['mag']= prtbl['mag']-0.7

#	GW170817-like
rtbl		= droutbl[droutbl['filter']=='r']; rtbl['delmjd'].sort()
rmag, rmagerr= tool.calc_app(rtbl['mag'], rtbl['magerr'], 38.4, 8.9, gwdist, gwdiststd)

#------------------------------------------------------------
#	PLOT 1	: TIME - MAG.	(r-band)
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 16})
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(9, 9))
#------------------------------------------------------------
#	GW170817-like
ax0.plot(rtbl['delmjd'], rmag, color='red', alpha=0.5, label='GW170817-like')
ax0.fill_between(rtbl['delmjd'], rmag-rmagerr, rmag+rmagerr, color='tomato', alpha=0.15, label='_nolegend_')
#	COCOON MODEL (Piro+2018)
ax0.plot(prtbl['delmjd'], prtbl['mag'], color='dodgerblue', alpha=0.5, label='Shock Cooling')
ax0.fill_between(prtbl['delmjd'], prtbl['mag']-prtbl['magerr'], prtbl['mag']+prtbl['magerr'], color='dodgerblue', alpha=0.3, label='_nolegend_')
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
ax0.set(xlabel='Time (Days from merger)', ylabel=r'Magnitude')
ax0.set_ylim([int(np.max(rmag))+1.0, int(np.min(prtbl['mag']))-1.0])
ax0.set_xlim([0,2])
plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM discovery')
ax0.legend(loc='upper right', prop={'size':20})
plt.tight_layout()
plt.minorticks_on()
plt.title('{0} r-band'.format(eventname))
plt.savefig('{0}/{1}_LC_rband.png'.format(save_path, eventname), overwrite=True)
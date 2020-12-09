#	KILONOVA LC PLOT
#	CREATED	2020.11.10	Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from imsng import phot, tool
from matplotlib.gridspec import GridSpec
import time
import bisect
#============================================================
#	FUNCTION
#============================================================

#============================================================
#	USER SETTING
#============================================================
#	PATH
#------------------------------------------------------------
path_imsng = '/data3/IMSNG/IMSNGgalaxies'
#------------------------------------------------------------
filist = dict(
				B='dodgerblue',
				V='lime',
				R='tomato',
				I='purple',
			)
obslist = ['LOAO', 'DOAO', 'SOAO', 'CBNUO']
t_lower = -1e1	#	[days]
t_upper = +1e2	#	[days]
#------------------------------------------------------------
#	Target to investigate
#------------------------------------------------------------
transient = 'SN2020uxz'
targ = 'NGC0514'
t0 = Time('2020-10-05T13:45:30', format='isot', scale='utc')
path_targ = '{}/{}'.format(path_imsng, targ)
#============================================================
#	Main part
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
ax0 = fig.add_subplot(111)
#------------------------------------------------------------
# i = 0
# obs = obslist[i]
for i, obs in enumerate(obslist):
	for j, filte in enumerate(list(filist.keys())):
		# filte = list(filist.keys())[2]
		c = filist[filte]
		offset = j*(1e-1)

		imlist = sorted(glob.glob('{}/{}/{}/Calib*-{}-*-{}-*.fits'.format(path_targ, obs, filte, obs, filte)))
		deltlist = []

		# k = 0
		# inim = imlist[k]
		for k, inim in enumerate(imlist):
			part = inim.split('-')
			utdate, uttime = part[3], part[4]
			dateobs = '{}-{}-{}T{}:{}:{}'.format(utdate[:4], utdate[4:6], utdate[6:], uttime[:2], uttime[2:4], uttime[4:])
			t = Time(dateobs, format='isot', scale='utc')
			delt = t.jd - t0.jd

			if (delt >= t_lower) & (delt <= t_upper):
				deltlist.append(delt)
		ax0.plot(deltlist, [i+offset]*len(deltlist), 'o', color=c, alpha=0.5)
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
ax0.set_xlabel('Days after {}'.format(t0.isot[:19]), fontsize=20)

ax0.axvline(x=0, color='k', linewidth=2, linestyle='--', label='Discovery')
ax0.tick_params(axis='x', which='minor', bottom=False)
ax0.grid(linestyle='--', color='grey', alpha=0.5)

ax0.set_title(targ, fontsize=20)
ax0.set_xticklabels(ax0.get_xticks().astype(dtype=int), fontsize=16)
ax0.set_yticks(np.arange(0, len(obslist)+1, 1))
ax0.set_yticklabels(obslist, fontsize=20)
#------------------------------------------------------------
fig.savefig('{}/{}/{}.{}.status.png'.format(path_imsng, targ, targ, transient), bbox_inches='tight', dpi=500, overwrite=True)
fig.savefig('{}/{}/{}.{}.status.pdf'.format(path_imsng, targ, targ, transient), bbox_inches='tight', overwrite=True)

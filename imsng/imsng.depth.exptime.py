#	This code is for drawing figure in pdf extinction for Paper
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.time import Time
#============================================================
#	Function
#------------------------------------------------------------
def imsng_dist(m):
	d = 50.0 * (10.**((m-19.5)/5.))
	return d
#------------------------------------------------------------
def depth(ul0, N):
	return ul0+2.5*np.log10(np.sqrt(N))
#============================================================
#	User setting
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_save = '/data1/IMSNG/4.plot'
#============================================================
#	PLOT
#------------------------------------------------------------
plt.rc('font', family='serif')
plt.close('all')
fig = plt.figure()
ax0 = fig.add_subplot(111)
#------------------------------------------------------------
x = 1920 / 2 / fig.dpi
y = 1080 / 2 / fig.dpi
fig.set_figwidth(x)
fig.set_figheight(y)
#------------------------------------------------------------
# ax0.spines['top'].set_color('none')
# ax0.spines['bottom'].set_color('none')
# ax0.spines['left'].set_color('none')
# ax0.spines['right'].set_color('none')
# ax0.tick_params(labelcolor='w', top=False,
# 				bottom=False, left=False, right=False)
ax0.set_xlabel('EXP.TIME [sec]', fontsize=20)
ax0.set_ylabel(r'$\rm 5\sigma$ depth [mag]', fontsize=20)
#============================================================
#	DEPTH
#------------------------------------------------------------
ul0 = 18.86
t0 = 300	#	[sec]
N = np.arange(t0, 10*t0+1, 1)/t0
ax0.plot(
			N, depth(ul0, N),
			linewidth='5',
			linestyle='-',
			color='orange',
			alpha=0.5,
			label='SOAO'
		)
ax0.axhline(y=19.5, color='tomato', linestyle='--', label='Required depth')
#------------------------------------------------------------
d, u = ax0.set_ylim()
ax0.legend(fontsize=20, framealpha=1.0)
ax0.set_ylim([u, d])
ax0.grid('both', color='silver', linestyle='--')
ax0.tick_params(axis="x", labelsize=14)
ax0.tick_params(axis="y", labelsize=14)
ax0.minorticks_on()
#============================================================
#	SAVE
#------------------------------------------------------------
fig.savefig("{}/imsng.soao.depth.pdf".format(path_save), bbox_inches='tight', overwrite=True)
fig.savefig("{}/imsng.soao.depth.png".format(path_save), dpi=500, bbox_inches='tight', overwrite=True)
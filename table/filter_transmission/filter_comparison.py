#%%
#	20.12.03	Gregory S.H. Paek
#============================================================
import os
import glob
import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d
from imsng import tool
#============================================================
#	USER SETTING
#------------------------------------------------------------
path_save = '/home/sonic/Research/gppy/table/filter_transmission'
# path_save = '/data1/AT2020yxz/4.plot'
#------------------------------------------------------------
path_filter = '/home/sonic/Research/gppy/table/filter_transmission'
# tels = ['2MASS_2MASS', 'TNG_NICS', 'UKIRT_UKIDSS', 'UKIRT_WFCAM', 'NOT_NOTcam']
tels = ['2MASS_2MASS', 'SMARTS_ANDIcam', 'TNG_NICS', 'UKIRT_UKIDSS', 'UKIRT_WFCAM', 'NOT_NOTcam']
#============================================================
#	Plot
#------------------------------------------------------------
plt.close('all')
plt.rc('font', family='serif')
fig = plt.figure()
figx = 1920 / fig.dpi
figy = 1080 / fig.dpi
fig.set_figwidth(figx)
fig.set_figheight(figy)
#------------------------------------------------------------
#	GRID
#------------------------------------------------------------
ax0 = fig.add_subplot(111)
#============================================================
#	Main
#------------------------------------------------------------
# tel = tels[0]
for tel in tels:
	groups = sorted(glob.glob('{}/{}*.dat'.format(path_filter, tel)))
	# filte = groups[0]
	for filte in groups:
		filtbl = ascii.read(filte)
		#------------------------------------------------------------
		#	Filter color
		if '.J' in filte:
			c = 'indigo'
			label = '{} J'.format(tel)
			alpha=1.0
		elif '.H' in filte:
			c = 'orange'
			label = '{} H'.format(tel)
			alpha=1.0
		elif '.K' in filte:
			c = 'tomato'
			label = '{} K'.format(tel)
			alpha=1.0
		#------------------------------------------------------------
		#	Telescope line style
		if tel == '2MASS_2MASS':
			ls = '-'
			c = 'grey'
			alpha=0.5
		elif tel == 'SMARTS_ANDIcam':
			ls = '-'
		elif tel == 'TNG_NICS':
			ls = '--'
		elif tel == 'UKIRT_UKIDSS':
			ls = '-.'
		elif tel == 'UKIRT_WFCAM':
			ls = ':'
		elif tel == 'NOT_NOTcam':
			ls = '-.'


		x = filtbl['wavelength']
		y = filtbl['response']

		ax0.plot(
					x, y,
					color=c,
					linestyle=ls,
					alpha=alpha,
					label=label,
				)
#------------------------------------------------------------
left, right = ax0.set_xlim()
ax0.set_xlim([left, 3*1e4])

ax0.tick_params(which='both', direction='in', labelsize=16)

ax0.set_xlabel(r'wavelength [$\rm \AA$]', fontsize=20)
ax0.set_ylabel('response', fontsize=20)

ax0.legend(fontsize=20, framealpha=1.0)
# ax0.set_title('GRB {}'.format(name), fontsize=16)
ax0.grid('both', linestyle='--', color='grey', alpha=0.5)


fig.savefig("{}/filter_transmission_JHK.pdf".format(path_save), bbox_inches='tight', overwrite=True)
fig.savefig("{}/filter_transmission_JHK.png".format(path_save), bbox_inches='tight', overwrite=True, dpi=500)

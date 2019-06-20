#   O3 RUN SUMMARY PLOT
#   2019.05.29 MADE		BY Gregory S.H. Paek
#============================================================#
#	MODULE
#------------------------------------------------------------#
import os, glob
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np
#------------------------------------------------------------#
path_table		= '/mnt/window/Users/User/Downloads/data/Project/gw/bayestar/sel/O3run_summary.dat'
path_save	= '/home/sonic/S190425z'

eventbl			= ascii.read(path_table)
eventbl			= eventbl[eventbl['noise']!='o']
#------------------------------------------------------------
#	PLOT 1	: TIME - MAG.
#------------------------------------------------------------
#	https://partrita.github.io/posts/matplotlib-examples/
#	Get the figure and the axes
plt.close('all')
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(10, 10))
plt.rcParams.update({'font.size': 20})
#ax0.fill_between(jdrange, imag-imagerr, imag+imagerr, color='dodgerblue', alpha=0.5, label='i-band')
i=0
for i in range(len(eventbl)):
	x	= eventbl['delmjd'][i]
	y	= eventbl['distmean'][i]
	yerr= eventbl['diststd'][i]
	capsize = 10

	size= np.sqrt(eventbl['region'][i]*2/np.pi)


	if		'Tera'in eventbl['prog'][i]:
		color	= 'grey'
	elif	'BHNS'in eventbl['prog'][i]:
		color	= 'blue'
	elif	'BNS' in eventbl['prog'][i]:
		color	= 'dodgerblue'
	elif	'BBH' in eventbl['prog'][i]:
		color	= 'dimgrey'
	elif	'Tera/BNS' in eventbl['prog'][i]:
		color	= 'violet'
	else:
		color	= 'brown'
	marker		= 'o'
#------------------------------------------------------------
	if		eventbl['noise'][i]=='o':
		marker	= 'x'
		color	= 'black'
		size	= 15
		yerr	= 0
		capsize = 0
#------------------------------------------------------------
	params_plot	= dict(	x=x, y=y, yerr=yerr,
						marker=marker, ms=size,
						c=color, alpha=0.5,
						capthick=1, capsize=capsize,
						label='_nolegend_')
	ax0.errorbar(**params_plot)
#------------------------------------------------------------
#	GW170817
#------------------------------------------------------------
params_plot	= dict(	x=-5, y=44.74, yerr=9,
					marker='o', ms=np.sqrt(30/np.pi),
					c='red', alpha=0.5,
					capthick=1, capsize=10,
					label='_nolegend_')
ax0.errorbar(**params_plot)
#------------------------------------------------------------
#	LEGEND (DUMMY POINTS)
#------------------------------------------------------------
ax0.errorbar(	x=-99, y=-99, yerr=0,
				marker='o',
				c='dodgerblue', alpha=0.4,
				label='BNS')
ax0.errorbar(	x=-99, y=-99, yerr=0,
				marker='o',
				c='blue', alpha=0.4,
				label='BNS/BHNS')
ax0.errorbar(	x=-99, y=-99, yerr=0,
				marker='o',
				c='violet', alpha=0.5,
				label='BNS/False Alarm')
ax0.errorbar(	x=-99, y=-99, yerr=0,
				marker='o',
				c='dimgrey', alpha=0.4,
				label='BBH')
'''
ax0.scatter(	x=-99, y=-99,
				marker='x', s=300,
				color='black', alpha=0.4,
				label='False Alarm')
'''
#	GW170817
params_plot	= dict(	x=-99, y=-99, yerr=0,
					marker='o',
					color='red', alpha=0.4,
					label='GW170817')
ax0.errorbar(**params_plot)
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
#fig.suptitle('S190425z', fontsize=14, fontweight='bold')
#ax0.set(title='LIGO/Virgo O3 run', xlabel=r'$t-t_{0}$ [days]', ylabel='GW Luminosity Distance [Mpc]')#, fontsize=14)
#ax0.set(xlabel=r'$t-t_{0}$ [days]', ylabel='GW Luminosity Distance [Mpc]')#, fontsize=14)
ax0.set(xlabel='Since O3 run start [days]', ylabel='GW Luminosity Distance [Mpc]')#, fontsize=14)
#ax0.set_ylim([24, 18])
ax0.set_xlim([0,np.max(eventbl['delmjd']+5)])
#ax0.set_yscale('log')
param_legend	= dict(	loc='upper left', fontsize=20,
						fancybox=True, framealpha=0.5,
						scatterpoints=1, markerscale=2)
plt.legend(**param_legend)
if 'O3_summary.png' in glob.glob('*png'):
	os.system('mv O3_summary.png O3_summary.png.bkg')
plt.minorticks_on()
plt.tight_layout()
plt.xticks(np.arange(0, 80, 10))
#------------------------------------------------------------
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax0, 20, loc='upper right') # zoom-factor: 2.5, location: upper-left
params_plot	= dict(	x=-5, y=44.74, yerr=9,
					marker='o', ms=np.sqrt(30*2/np.pi),
					c='red', alpha=1.0,
					capthick=1, capsize=10,
					label='_nolegend_')
axins.errorbar(**params_plot)
x1, x2, y1, y2 = -5.25, -4.75, 30, 60 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
axins.minorticks_on()
plt.yticks(visible=True)
plt.xticks(visible=False)
plt.yticks(np.arange(30, 61, 10))

plt.savefig(path_save+'/Figure_X+O3_summary.png')
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
eventbl			= ascii.read(path_table)
#------------------------------------------------------------
#	PLOT 1	: TIME - MAG.
#------------------------------------------------------------
#	https://partrita.github.io/posts/matplotlib-examples/
#	Get the figure and the axes
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(8, 6))
plt.rcParams.update({'font.size': 18})
#ax0.fill_between(jdrange, imag-imagerr, imag+imagerr, color='dodgerblue', alpha=0.5, label='i-band')


i=0
for i in range(len(eventbl)):
	x	= eventbl['delmjd'][i]
	y	= eventbl['distmean'][i]
	yerr= eventbl['diststd'][i]
	capsize = 10
	if		'Tera'in eventbl['prog'][i]:
		color	= 'grey'
		size	= 100
	elif	'BHNS'in eventbl['prog'][i]:
		color	= 'blue'
		size	= eventbl['region'][i]/100
	elif	'BNS' in eventbl['prog'][i]:
		color	= 'dodgerblue'
		size	= eventbl['region'][i]/100
	elif	'BBH' in eventbl['prog'][i]:
		color	= 'dimgrey'
		size	= eventbl['region'][i]/100
	else:
		color	= 'brown'
		size	= eventbl['region'][i]/100
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
						c=color, alpha=0.4,
						capthick=1, capsize=capsize,
						label='_nolegend_')
	ax0.errorbar(**params_plot)
#------------------------------------------------------------
#	GW170817
#------------------------------------------------------------
params_plot	= dict(	x=0, y=44.74, yerr=9,
					marker='*', ms=30,
					c='dodgerblue', alpha=0.4,
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
				c='dimgrey', alpha=0.4,
				label='BBH')
'''
ax0.errorbar(	x=-99, y=-99, yerr=0,
				marker='x',
				c='black', alpha=0.4,
				label='False')
'''
ax0.scatter(	x=-99, y=-99,
				marker='x', s=300,
				color='black', alpha=0.4,
				label='False Alarm')
#	GW170817
params_plot	= dict(	x=-99, y=-99,
					marker='*', s=600,
					color='dodgerblue', alpha=0.4,
					label='GW170817')
ax0.scatter(**params_plot)
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------



#fig.suptitle('S190425z', fontsize=14, fontweight='bold')
#ax0.set(title='LIGO/Virgo O3 run', xlabel=r'$t-t_{0}$ [days]', ylabel='GW Luminosity Distance [Mpc]')#, fontsize=14)
ax0.set(xlabel=r'$t-t_{0}$ [days]', ylabel='GW Luminosity Distance [Mpc]')#, fontsize=14)
#ax0.set_ylim([24, 18])
ax0.set_xlim([-5,np.max(eventbl['delmjd']+5)])
ax0.set_yscale('log')
ax0.legend()
if 'O3_summary.png' in glob.glob('*png'):
	os.system('mv O3_summary.png O3_summary.png.bkg')
plt.savefig('O3_summary.png')

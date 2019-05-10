#============================================================
#	LIGHT CURVE FUNCTION FOR GRB/SN/... etc.
#	2019.05.03. MADE BY GREGORY S.H. PAEK
#============================================================
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy.optimize import curve_fit
#------------------------------------------------------------
def LCplot(ax, t0, t, mag, magerr, offset=0, label='Somewhere', marker='o', color='grey', scale='linear'):
	params_plot	= dict(	x=t-t0,
						y=mag,
						yerr=magerr,
						marker=marker, ms=6, mew=1,
						mec='k', mfc=color, fmt='ko',
						capsize=5, capthick=2,
						label=label)
	ax.errorbar(**params_plot)
#------------------------------------------------------------
phot	= ascii.read('/mnt/window/Users/User/Downloads/data/Project/AT2019ein/result.dat')
target	= 'AT2019ein'
'''
bandlist= []
for filename in phot['obs']:
	if '-B-' in filename: bandlist.append('B')
	if '-V-' in filename: bandlist.append('V')
	if '-R-' in filename: bandlist.append('R')
	if '-I-' in filename: bandlist.append('I')
bandlist= np.array(bandlist)

jdlist	= []
for date_obs in phot['date-obs']:
	jdlist.append(Time(date_obs, format='isot', scale='utc').jd)
jdlist= np.array(jdlist)

obslist	= []
for filename in phot['obs']:
	if '-LOAO-' in filename:	obslist.append('LOAO')
	if '-DOAO' in filename:		obslist.append('DOAO')
	if '-SAO-' in filename:		obslist.append('SAO')
	if '-UKIRT-' in filename:	obslist.append('UKIRT')
	if '-OTHER-' in filename:	obslist.append('OTHER')
obslist= np.array(obslist)

phot['band']	= bandlist
phot['jd']		= jdlist
phot['obs']		= obslist
'''
#------------------------------------------------------------
fig		= plt.figure()				# Define "figure" instance
fig.set_size_inches(8,6)			# Physical page size in inches, (lx,ly)
suptit	= target+' Light curve'
fig.suptitle(suptit, fontsize=15)	# Title for the page
ax		= plt.subplot()

t0		= 2458604.9372453704

photB_LOAO	= phot[	(phot['band']=='B')&
					(phot['obs']=='LOAO')]
photR_LOAO	= phot[	(phot['band']=='R')&
					(phot['obs']=='LOAO')]

photR_DOAO	= phot[	(phot['band']=='R')&
					(phot['obs']=='DOAO')]

photR_SOAO	= phot[	(phot['band']=='R')&
					(phot['obs']=='SOAO')]

photlist	= [	photB_LOAO, photR_LOAO,
				photR_DOAO, photR_SOAO]
#------------------------------------------------------------
'''
params_LCplot_B	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photB_LOAO['jd']),
						mag=np.copy(photB_LOAO['mag']),
						magerr=np.copy(photB_LOAO['magerr']),
						label='LOAO B-band',
						marker='o',
						color='dodgerblue')
LCplot(**params_LCplot_B)
'''
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_LOAO['jd']),
						mag=np.copy(photR_LOAO['mag']),
						magerr=np.copy(photR_LOAO['magerr']),
						label='LOAO R-band',
						marker='o',
						color='tomato')
LCplot(**params_LCplot_R)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_DOAO['jd']),
						mag=np.copy(photR_DOAO['mag']),
						magerr=np.copy(photR_DOAO['magerr']),
						label='DOAO R-band',
						marker='s',
						color='tomato')
LCplot(**params_LCplot_R)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_SOAO['jd']),
						mag=np.copy(photR_SOAO['mag']),
						magerr=np.copy(photR_SOAO['magerr']),
						label='SOAO R-band',
						marker='d',
						color='tomato')
LCplot(**params_LCplot_R)
#------------------------------------------------------------
ax1		= plt.gca().invert_yaxis()
#ax1		= plt.xscale('log')
ax1		= plt.xlabel(r'Time since $T_0$ (days)', {'color': 'black', 'fontsize': 10})
ax1		= plt.ylabel('magnitude [AB]', {'color': 'black', 'fontsize': 10})
ax1		= plt.legend(loc='best')

filename	= target+'_LC'
plt.savefig(filename+'.png', dpi = 500)
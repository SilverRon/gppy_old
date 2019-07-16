#============================================================
#	LIGHT CURVE FUNCTION FOR GRB/SN/... etc.
#	2019.05.03. MADE BY GREGORY S.H. PAEK
#	2019.07.16. MADE BY GREGORY S.H. PAEK
#============================================================
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy.optimize import curve_fit
import pylab as P
#------------------------------------------------------------
def LCplot(ax, t0, t, mag, magerr, offset=0, label='Somewhere', marker='o', color='grey', scale='linear'):
	params_plot	= dict(	x=t-t0,
						y=mag+offset,
						yerr=magerr,
						marker=marker, ms=8, mew=1,
						mec='k', mfc=color, fmt='ko',
						capsize=5, capthick=2,
						label=label)
	ax.errorbar(**params_plot)
#------------------------------------------------------------
path_save	= '/home/sonic/SN2019ehk'
phot0	= ascii.read(path_save+'/phot4SN2019ehk.dat')
#phot0	= ascii.read('/mnt/window/Users/User/Downloads/data/Project/SN2019ehk/result.dat')
phot	= phot0[phot0['magerr']!=-99]
nd		= phot0[phot0['magerr']==-99]
target	= 'SN2019ehk'

#------------------------------------------------------------
plt.close('all')
fig		= plt.figure()				# Define "figure" instance
fig.set_size_inches(8,6)			# Physical page size in inches, (lx,ly)
suptit	= target+' Light curve'
fig.suptitle(suptit, fontsize=15)	# Title for the page
ax		= plt.subplot()

t0		= 2458603.072222			#	LSGT Detection
#off		= -0.1
off		= 0.0
#------------------------------------------------------------
photB_LOAO	= phot[	(phot['filter']=='B')&
					(phot['obs']=='LOAO')]
photR_LOAO	= phot[	(phot['filter']=='R')&
					(phot['obs']=='LOAO')]
#------------------------------------------------------------
photB_DOAO	= phot[	(phot['filter']=='B')&
					(phot['obs']=='DOAO')]
photV_DOAO	= phot[	(phot['filter']=='V')&
					(phot['obs']=='DOAO')]
photR_DOAO	= phot[	(phot['filter']=='R')&
					(phot['obs']=='DOAO')]
#------------------------------------------------------------
photR_SOAO	= phot[	(phot['filter']=='R')&
					(phot['obs']=='SOAO')]
#------------------------------------------------------------
photR_MCD30	= phot[	(phot['filter']=='R')&
					(phot['obs']=='MCD30INCH')]
#============================================================
photlist	= [	photB_LOAO, photR_LOAO,
				photB_DOAO, photV_DOAO, photR_DOAO,
				photR_SOAO]
#------------------------------------------------------------
params_LCplot_B	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photB_LOAO['jd']),
						mag=np.copy(photB_LOAO['mag']),
						magerr=np.copy(photB_LOAO['magerr']),
						label='LOAO B-band',
						marker='o',
						color='dodgerblue')
LCplot(**params_LCplot_B)
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
params_LCplot_B	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photB_DOAO['jd']),
						mag=np.copy(photB_DOAO['mag']),
						magerr=np.copy(photB_DOAO['magerr']),
						offset=off,
						label='DOAO B-band',
						marker='s',
						color='dodgerblue')
LCplot(**params_LCplot_B)
#------------------------------------------------------------
params_LCplot_V	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photV_DOAO['jd']),
						mag=np.copy(photV_DOAO['mag']),
						magerr=np.copy(photV_DOAO['magerr']),
						offset=off,
						label='DOAO V-band',
						marker='s',
						color='green')
LCplot(**params_LCplot_V)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_DOAO['jd']),
						mag=np.copy(photR_DOAO['mag']),
						magerr=np.copy(photR_DOAO['magerr']),
						offset=off,
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
						marker='D',
						color='tomato')
LCplot(**params_LCplot_R)
#------------------------------------------------------------
ndB, ndV, ndR	= nd[nd['filter']=='B'], nd[nd['filter']=='V'], nd[nd['filter']=='R']
ndBx, ndBy	= ndB['jd'][ndB['jd']==np.min(ndB['jd'])].data-t0, ndB['ul'][ndB['jd']==np.min(ndB['jd'])].data
#ndVx, ndVy	= ndV['jd'][ndV['jd']==np.min(ndV['jd'])]-t0, ndV['mag'][ndV['jd']==np.min(ndV['jd'])]
ndRx, ndRy	= ndR['jd'][ndR['jd']==np.min(ndR['jd'])].data-t0, ndR['ul'][ndR['jd']==np.min(ndR['jd'])].data

plt.scatter(ndBx, ndBy, marker='v', c='dodgerblue', s=100)
plt.scatter(ndRx, ndRy, marker='v', c='tomato', s=100)

#P.arrow(ndBx, ndBy, 0.0, +0.2, fc='dodgerblue', ec='k', head_width=0.1, head_length=0.1)       
#P.arrow(ndRx, ndRy, 0.0, +0.2, fc='tomato', ec='k', head_width=0.1, head_length=0.1)       
#------------------------------------------------------------


ax1		= plt.gca().invert_yaxis()
#ax1		= plt.xscale('log')
ax1		= plt.xlabel(r'Time since $T_0$ (days)', {'color': 'black', 'fontsize': 16})
ax1		= plt.ylabel('magnitude [AB]', {'color': 'black', 'fontsize': 16})
ax1		= plt.legend(loc='best')
#ax1		= plt.axhline(y=21, color='gray', linestyle='--')

filename	= target+'_LC'
plt.minorticks_on()
plt.savefig(path_save+'/'+filename+'.png')#, dpi = 500)





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
def LCplot(ax, t0, t, mag, magerr, offset=0, label='Somewhere', marker='o', color='grey', scale='linear', alpha=0.5):
	params_plot	= dict(	x=t-t0,
						y=mag+offset,
						yerr=magerr,
						marker=marker, ms=15, mew=1,
						mec='k', mfc=color, fmt='ko', alpha=0.7,
						capsize=5, capthick=1,
						label=label)
	ax.errorbar(**params_plot)
#------------------------------------------------------------
path_save	= '/home/sonic/SN2020oi'
phot0	= ascii.read(path_save+'/phot4SN2020oi.dat')
#phot0	= ascii.read('/mnt/window/Users/User/Downloads/data/Project/SN2020oi/result.dat')
phot	= phot0[phot0['magerr']!=-99]
nd		= phot0[phot0['magerr']==-99]
target	= 'SN2020oi'

#------------------------------------------------------------
plt.close('all')
fig		= plt.figure()				# Define "figure" instance
fig.set_size_inches(8,6)			# Physical page size in inches, (lx,ly)
suptit	= target+' Light curve'
ax		= plt.subplot()

t0 = 2458856.04229					#	ZTF DISCOVERY
t0_ut = '2020-01-07T13:00:54'
#off		= -0.1
off		= 0.0
ax1 = plt.axvline(x=0, linestyle='--', color='grey')
#------------------------------------------------------------
photB_LOAO	= phot[	(phot['filter']=='B')&
					(phot['obs']=='LOAO')]
photR_LOAO	= phot[	(phot['filter']=='R')&
					(phot['obs']=='LOAO')]
#------------------------------------------------------------
photB_SAO	= phot[	(phot['filter']=='B')&
                  	(phot['obs'] == 'SAO_STX16803')]
photV_SAO	= phot[	(phot['filter']=='V')&
                  	(phot['obs'] == 'SAO_STX16803')]
photR_SAO	= phot[	(phot['filter']=='R')&
                  	(phot['obs'] == 'SAO_STX16803')]
photI_SAO	= phot[	(phot['filter']=='I')&
					(phot['obs'] == 'SAO_STX16803')]
#------------------------------------------------------------
photB_DOAO	= phot[	(phot['filter']=='B')&
					(phot['obs']=='DOAO_FLI')]
photV_DOAO	= phot[	(phot['filter']=='V')&
					(phot['obs']=='DOAO_FLI')]
photR_DOAO	= phot[	(phot['filter']=='R')&
					(phot['obs']=='DOAO_FLI')]
photI_DOAO	= phot[	(phot['filter']=='I')&
					(phot['obs']=='DOAO_FLI')]
#------------------------------------------------------------
photR_SOAO	= phot[	(phot['filter']=='R')&
					(phot['obs']=='SOAO')]
#------------------------------------------------------------
photR_KHAO	= phot[	(phot['filter']=='R')&
					(phot['obs']=='KHAO_MDFTS')]
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
						label='LOAO B',
						marker='o',
						color='dodgerblue',
						alpha=0.5,
						)
# LCplot(**params_LCplot_B)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_LOAO['jd']),
						mag=np.copy(photR_LOAO['mag']),
						magerr=np.copy(photR_LOAO['magerr']),
						label='LOAO R',
						marker='o',
						color='tomato',
						alpha=0.5,
						)
# LCplot(**params_LCplot_R)
#------------------------------------------------------------
params_LCplot_B	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photB_SAO['jd']),
						mag=np.copy(photB_SAO['mag']),
						magerr=np.copy(photB_SAO['magerr']),
						label='SAO_STX16803 B',
						marker='*',
						color='dodgerblue',
						alpha=0.5,
						)
LCplot(**params_LCplot_B)
#------------------------------------------------------------
params_LCplot_V	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photV_SAO['jd']),
						mag=np.copy(photV_SAO['mag']),
						magerr=np.copy(photV_SAO['magerr']),
						label='SAO_STX16803 V',
						marker='*',
						color='green',
						alpha=0.5,
						)
LCplot(**params_LCplot_V)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_SAO['jd']),
						mag=np.copy(photR_SAO['mag']),
						magerr=np.copy(photR_SAO['magerr']),
						label='SAO_STX16803 R',
						marker='*',
						color='tomato',
						alpha=0.5,
						)
LCplot(**params_LCplot_R)
#------------------------------------------------------------
params_LCplot_I	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photI_SAO['jd']),
						mag=np.copy(photI_SAO['mag']),
						magerr=np.copy(photI_SAO['magerr']),
						label='SAO_STX16803 I',
						marker='*',
						color='purple',
						alpha=0.5,
						)
LCplot(**params_LCplot_I)
#------------------------------------------------------------
params_LCplot_B	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photB_DOAO['jd']),
						mag=np.copy(photB_DOAO['mag']),
						magerr=np.copy(photB_DOAO['magerr']),
						offset=off,
						label='DOAO_FLI B',
						marker='s',
						color='dodgerblue',
						alpha=0.5,
						)
# LCplot(**params_LCplot_B)
#------------------------------------------------------------
params_LCplot_V	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photV_DOAO['jd']),
						mag=np.copy(photV_DOAO['mag']),
						magerr=np.copy(photV_DOAO['magerr']),
						offset=off,
						label='DOAO_FLI V',
						marker='s',
						color='green',
						alpha=0.5,
						)
# LCplot(**params_LCplot_V)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_DOAO['jd']),
						mag=np.copy(photR_DOAO['mag']),
						magerr=np.copy(photR_DOAO['magerr']),
						offset=off,
						label='DOAO_FLI R',
						marker='s',
						color='tomato',
						alpha=0.5,
						)
LCplot(**params_LCplot_R)
#------------------------------------------------------------
params_LCplot_I	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photI_DOAO['jd']),
						mag=np.copy(photI_DOAO['mag']),
						magerr=np.copy(photI_DOAO['magerr']),
						offset=off,
						label='DOAO_FLI I',
						marker='s',
						color='purple',
						alpha=0.5,
						)
LCplot(**params_LCplot_I)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_SOAO['jd']),
						mag=np.copy(photR_SOAO['mag']),
						magerr=np.copy(photR_SOAO['magerr']),
						label='SOAO R',
						marker='D',
						color='tomato',
						alpha=0.5,
						)
LCplot(**params_LCplot_R)
#------------------------------------------------------------
params_LCplot_R	= dict(	ax=ax,
						t0=t0,
						t=np.copy(photR_KHAO['jd']),
						mag=np.copy(photR_KHAO['mag']),
						magerr=np.copy(photR_KHAO['magerr']),
						label='KHAO_MDFTS R',
						marker='.',
						color='tomato',
						alpha=0.5,
						)
LCplot(**params_LCplot_R)
#------------------------------------------------------------
ndB, ndV, ndR	= nd[nd['filter']=='B'], nd[nd['filter']=='V'], nd[nd['filter']=='R']
'''
ndBx, ndBy	= ndB['jd'][ndB['jd']==np.min(ndB['jd'])].data-t0, ndB['ul_3sig'][ndB['jd']==np.min(ndB['jd'])].data
ndVx, ndVy	= ndV['jd'][ndV['jd']==np.min(ndV['jd'])]-t0, ndV['mag'][ndV['jd']==np.min(ndV['jd'])]
ndRx, ndRy	= ndR['jd'][ndR['jd']==np.min(ndR['jd'])].data-t0, ndR['ul_3sig'][ndR['jd']==np.min(ndR['jd'])].data
'''
# ndBx, ndBy	= ndB['jd'][ndB['jd']==np.min(ndB['jd'])].data-t0, ndB['ul_3sig'][ndB['jd']==np.min(ndB['jd'])].data
# ndVx, ndVy	= ndV['jd'][ndV['jd']==np.min(ndV['jd'])]-t0, ndV['mag'][ndV['jd']==np.min(ndV['jd'])]
ndRx, ndRy	= ndR['jd'].data-t0, ndR['ul_3sig'].data

# plt.scatter(ndBx, ndBy, marker='v', c='dodgerblue', s=100)
plt.scatter(ndRx, ndRy, marker='v', c='tomato', s=200)

#P.arrow(ndBx, ndBy, 0.0, +0.2, fc='dodgerblue', ec='k', head_width=0.1, head_length=0.1)       
#P.arrow(ndRx, ndRy, 0.0, +0.2, fc='tomato', ec='k', head_width=0.1, head_length=0.1)       
#------------------------------------------------------------


ax1		= plt.gca().invert_yaxis()
#ax1		= plt.xscale('log')
ax1		= plt.xlabel(r'Time since $T_0$ (days)', {'color': 'black', 'fontsize': 20})
ax1		= plt.ylabel('magnitude [AB]', {'color': 'black', 'fontsize': 20})
ax1		= plt.legend(loc='best', fontsize=15)
#ax1		= plt.axhline(y=21, color='gray', linestyle='--')
filename	= target+'_LC'

left, right = plt.xlim()
up, down = plt.ylim()
plt.xlim(0, right)
# plt.xlim(2, 2.4)
# plt.ylim(16.5, 15)
plt.title(suptit, fontsize=20)
plt.tick_params(which='both', direction='in', labelsize=15)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True, linestyle='--', color='grey', alpha=0.75)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(path_save+'/'+filename+'.png', dpi = 300)





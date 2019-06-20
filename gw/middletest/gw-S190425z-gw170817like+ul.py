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
#	FUNCTION
#------------------------------------------------------------
def GW170817_like(gwdist, gwdiststd):
	import numpy as np
	m0		= 17.476	#	[AB] in i-band (t0+10h)
	m0err	= 0.018	
	dist0	= 38.4		#	[MPC]	Im et al. 2017
	dist0err= 8.9
	m		= m0+5.*np.log10(gwdist/dist0)
	merr	= np.sqrt( (m0err)**2 + ((5.*gwdiststd)/(gwdist*np.log(10)))**2 + ((5.*dist0err)/(dist0*np.log(10)))**2 )
	return m, merr
#------------------------------------------------------------
def func_linear(a, x, scaling=[0, 0]):
	xpt, ypt= scaling[0], scaling[1]
	ydel	= ypt - (-1*a*xpt)
	return -1*a*x + ydel
#------------------------------------------------------------
def calc_app(mag, magerr, gwdist0, gwdiststd0, gwdist1, gwdiststd1):
	app		= mag+5*np.log10(gwdist1/gwdist0)
	apperr	= np.sqrt( (magerr)**2 + ((5*gwdiststd1)/(np.log(5)*gwdist1))**2 + ((5*gwdiststd0)/(np.log(5)*gwdist0))**2 )
	return app, apperr
#------------------------------------------------------------
def abs2app(mag, magerr, gwdist, gwdiststd):
	app		= 5*np.log10(gwdist)-5+mag
	apperr	= 5*gwdiststd/(gwdist*np.log(10))
	return app, apperr
#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/20190606'
path_save	= '/home/sonic/S190425z'
gwdist		= 156.1439917314773
gwdiststd	= 41.3733570000998
t0			= Time('2019-04-25T08:18:05.017553', format='isot', scale='utc')
jd0_170817	= 2457983.02852
droutbl		= ascii.read('/mnt/window/Users/User/Downloads/data/Project/gw/phot_gw170817_Drout.dat')
pirotbl		= ascii.read(path_save+'/table/cocoon_Piro+18.dat')
#------------------------------------------------------------
#	COCOON MODEL
prtbl, pitbl, pjtbl		=\
pirotbl[pirotbl['filter']=='r'], pirotbl[pirotbl['filter']=='i'], pirotbl[pirotbl['filter']=='J']

prtbl['mag'], prtbl['magerr']	=\
abs2app(prtbl['absmag'], 0, gwdist*1e6, gwdiststd*1e6)
pitbl['mag'], pitbl['magerr']	=\
abs2app(pitbl['absmag'], 0, gwdist*1e6, gwdiststd*1e6)
pjtbl['mag'], pjtbl['magerr']	=\
abs2app(pjtbl['absmag'], 0, gwdist*1e6, gwdiststd*1e6)

prtbl['mag']	= prtbl['mag']-0.7
pitbl['mag']	= pitbl['mag']+0.0
pjtbl['mag']	= pjtbl['mag']+1.0



#	GW170817-like
rtbl, itbl, jtbl, ktbl	=\
droutbl[droutbl['filter']=='r'], droutbl[droutbl['filter']=='i'], droutbl[droutbl['filter']=='J'], droutbl[droutbl['filter']=='K']
rtbl['delmjd'].sort(), itbl['delmjd'].sort(), jtbl['delmjd'].sort(), ktbl['delmjd'].sort()
#photbl		= ascii.read(path_base+'/phot_update.dat')
photbl		= ascii.read('/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/20190609/phot_update_allnew.dat')
photbl['delmjd']	= photbl['jd']-t0.jd

params_obs	= {	'KMTNET':'purple-v-150-R',
				'LOAO'	:'chocolate-v-150-R',
				'LSGT'	:'dimgrey-v-150-r',
				'SAO'	:'green-v-150-R',
				'SQUEAN':'violet-v-200-i',
				'UKIRT'	:'orange-v-200-J'}
#------------------------------------------------------------
#	PLOT 1	: TIME - MAG.	(r-band)
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 16})
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(12, 9))
#------------------------------------------------------------
#	GW170817-like
rmag, rmagerr	= calc_app(rtbl['mag'], rtbl['magerr'], 38.4, 8.9, gwdist, gwdiststd)
ax0.plot(rtbl['delmjd'], rmag, color='red', alpha=0.5, label='GW170817-like (r-band)')
ax0.fill_between(rtbl['delmjd'], rmag-rmagerr, rmag+rmagerr, color='tomato', alpha=0.15, label='_nolegend_')
#	COCOON MODEL (Piro+2018)
ax0.plot(prtbl['delmjd'], prtbl['mag'], color='dodgerblue', alpha=0.5, label='Shock Cooling (r-band)')
ax0.fill_between(prtbl['delmjd'], prtbl['mag']-prtbl['magerr'], prtbl['mag']+prtbl['magerr'], color='dodgerblue', alpha=0.3, label='_nolegend_')
#------------------------------------------------------------
#	NEW PHOT TABLE FOR PPT IMAGE

#------------------------------------------------------------
for obs in ['KMTNET', 'LOAO', 'LSGT', 'SAO']:
	obstbl	= photbl[photbl['obs']==obs]
	sel	= params_obs[obs].split('-')
#	color, marker, size, band	= sel[0], sel[1], int(sel[2]), sel[3]
	color, marker, size, band	= sel[0], sel[1], int(sel[2])*2, sel[3]
	label	= '{}({}-band)'.format(obs, band)
	params_plot	= dict(	x=obstbl['jd']-t0.jd, y=obstbl['ul'],
						c=color, s=size, marker=marker,
						label=label, alpha=0.5)
	ax0.scatter(**params_plot)
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
ax0.set(xlabel='Time (Days from merger)', ylabel=r'Limit Magnitude ($3 \sigma$)')
ax0.set_ylim([24, 18])
ax0.set_xlim([0,3])
plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM counterpart discovery')
ax0.legend(loc='upper right', prop={'size':20})
plt.tight_layout()
plt.minorticks_on()
plt.savefig(path_save+'/Figure_X+S190425z_Rr-band.png', overwrite=True)
#------------------------------------------------------------
#	PLOT 2	: TIME - MAG.	(i-band)
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 16})
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(12, 9))
#------------------------------------------------------------
#	GW170817-like
imag, imagerr	= calc_app(itbl['mag'], itbl['magerr'], 38.4, 8.9, gwdist, gwdiststd)
ax0.plot(itbl['delmjd'], imag, color='violet', alpha=0.5, label='GW170817-like (i-band)')
ax0.fill_between(itbl['delmjd'], imag-imagerr, imag+imagerr, color='violet', alpha=0.15, label='_nolegend_')
#	COCOON MODEL (Piro+2018)
ax0.plot(pitbl['delmjd'], pitbl['mag'], color='dodgerblue', alpha=0.5, label='Shock Cooling (i-band)')
ax0.fill_between(pitbl['delmjd'], pitbl['mag']-pitbl['magerr'], pitbl['mag']+pitbl['magerr'], color='dodgerblue', alpha=0.3, label='_nolegend_')
#------------------------------------------------------------
for obs in ['SQUEAN']:
	obstbl	= photbl[photbl['obs']==obs]
	sel	= params_obs[obs].split('-')
	color, marker, size, band	= sel[0], sel[1], int(sel[2]), sel[3]
	label	= '{}({}-band)'.format(obs, band)
	params_plot	= dict(	x=obstbl['jd']-t0.jd, y=obstbl['ul'],
						c=color, s=size, marker=marker,
						label=label, alpha=0.5)
	ax0.scatter(**params_plot)
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
ax0.set(xlabel='Time (Days from merger)', ylabel=r'Limit Magnitude ($3 \sigma$)')
ax0.set_ylim([24, 18])
ax0.set_xlim([0,3.5])
plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM counterpart discovery')
ax0.legend(loc='upper right', prop={'size':20})
plt.tight_layout()
plt.minorticks_on()
plt.savefig(path_save+'/Figure_X+S190425z_i-band.png', overwrite=True)
#------------------------------------------------------------
#	PLOT 3	: TIME - MAG.	(J,K-band)
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 16})
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(12, 9))
#------------------------------------------------------------
#	GW170817-like
jmag, jmagerr	= calc_app(jtbl['mag'], jtbl['magerr'], 38.4, 8.9, gwdist, gwdiststd)
ax0.plot(jtbl['delmjd'], jmag, color='darkorange', alpha=0.5, label='GW170817-like (J-band)')
ax0.fill_between(jtbl['delmjd'], jmag-jmagerr, jmag+jmagerr, color='darkorange', alpha=0.15, label='_nolegend_')
#	COCOON MODEL (Piro+2018)
ax0.plot(pjtbl['delmjd'], pjtbl['mag'], color='dodgerblue', alpha=0.5, label='Shock Cooling (J-band)')
ax0.fill_between(pjtbl['delmjd'], pjtbl['mag']-pjtbl['magerr'], pjtbl['mag']+pjtbl['magerr'], color='dodgerblue', alpha=0.3, label='_nolegend_')
for obs in ['UKIRT']:
	obstbl	= photbl[ (photbl['obs']==obs) & (photbl['filter']=='J')]
	sel	= params_obs[obs].split('-')
	color, marker, size, band	= sel[0], sel[1], int(sel[2]), sel[3]
	label	= '{}({}-band)'.format(obs, band)
	params_plot	= dict(	x=obstbl['jd']-t0.jd, y=obstbl['ul'],
						c=color, s=size, marker=marker,
						label=label, alpha=0.5)
	ax0.scatter(**params_plot)
#------------------------------------------------------------
kmag, kmagerr	= calc_app(ktbl['mag'], ktbl['magerr'], 38.4, 8.9, gwdist, gwdiststd)
ax0.plot(ktbl['delmjd'], kmag, color='g', alpha=0.5, label='GW170817-like (K-band)')
ax0.fill_between(ktbl['delmjd'], kmag-kmagerr, kmag+kmagerr, color='g', alpha=0.15, label='_nolegend_')
for obs in ['UKIRT']:
	obstbl	= photbl[ (photbl['obs']==obs) & (photbl['filter']=='K')]
	sel	= params_obs[obs].split('-')
	color, marker, size, band	= 'green', sel[1], int(sel[2]), 'K'
	label	= '{}({}-band)'.format(obs, band)
	params_plot	= dict(	x=obstbl['jd']-t0.jd, y=obstbl['ul'],
						c=color, s=size, marker=marker,
						label=label, alpha=0.5)
	ax0.scatter(**params_plot)
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
ax0.set(xlabel='Time (Days from merger)', ylabel=r'Limit Magnitude ($3 \sigma$)')
ax0.set_ylim([24, 18])
ax0.set_xlim([0,7.5])
plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM counterpart discovery')
ax0.legend(loc='upper right', prop={'size':20})
plt.tight_layout()
plt.minorticks_on()
plt.savefig(path_save+'/Figure_X+S190425z_JK-band.png', overwrite=True)
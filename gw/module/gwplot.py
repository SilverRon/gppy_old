#============================================================
#	GW PLOT MODULE
#	2019.08.11	CREATED BY	Gregory S.H. Paek
#	2019.10.09	UPDATED BY	Gregory S.H. Paek
#============================================================
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
import ligo.skymap.plot
from astropy import units as u
from astropy.io import fits, ascii
import os, glob
import numpy as np
from astropy.io import ascii
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy.interpolate import interp1d
# import ligo.skymap
#============================================================
#	FUNCTION
#------------------------------------------------------------
def skymap(bayestarfits, photbl):
	fig = plt.figure(figsize=(12, 6), dpi=100)
	ax = plt.axes(
		[0.05, 0.05, 0.9, 0.9],
		projection='astro hours mollweide',)
	ax.imshow_hpx(bayestarfits, cmap='cylon')
	for obs in setting_obsplot.keys():
		indx = np.where(obs.upper() == photbl['obs'])
		tc = SkyCoord(photbl['ra'][indx], photbl['dec'][indx], frame='icrs', unit='deg')
		ax.plot(tc.ra.deg, tc.dec.deg,
				transform=ax.get_transform('world'),
				**setting_obsplot[obs],
				 markeredgewidth=2)
	ax.grid()
	ax.minorticks_on()
	ax.legend(fontsize=15,
			# loc='lower left',
			loc='upper right',
			markerscale=2,
			framealpha=1)
	ax.tick_params(axis='both', which='major', labelsize=20)
#------------------------------------------------------------
def cumscore(obslist, photbl, timestep, t0, setting_obsplot):
	fig = plt.figure(figsize=(8, 6), dpi=100)
	for obs in obslist:
		indx = np.where(obs.upper() == photbl['obs'])
		onephotbl = photbl[indx]
		tcut = 0.0
		x = []
		y = []
		yerr = []
		for i in range(10):
			tindx = np.where(	(onephotbl['jd']-t0.jd < tcut+timestep)&
								(onephotbl['jd']-t0.jd > tcut))
			tcut += timestep
			if len(tindx[0]) != 0:
				x.append(np.median(onephotbl['jd'][tindx])-t0.jd)
				y.append(np.median(onephotbl['ul'][tindx]))
				yerr.append(np.std(onephotbl['ul'][tindx]))
		if len(x) != 0:
			plt.errorbar(x, y, yerr=yerr,
						c='k', capsize=10,
						ls='', marker='v', markersize=20,
						mfc=setting_obsplot[obs]['color'], mec='k',
						label=setting_obsplot[obs]['label'])
	plt.xlabel(r'Days after $t_0$', fontsize=20)
	plt.ylabel('AB magnitude', fontsize=20)
	plt.legend(fontsize=20)
	plt.minorticks_on()
	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.gca().invert_yaxis()
	plt.tight_layout()
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
#	PATH
path_initial = '/data1/S190425z/info/Initial/S190425z_Initial-bayestar.fits.gz'
path_update = '/data1/S190425z/info/Update/S190425z_Update-bayestar.fits.gz'
path_phot_ini = '/data1/S190425z/1.result/Initial/phot+score4initial.dat'
path_phot_upd = '/data1/S190425z/1.result/Update/phot+score4update.dat'
path_cand_ini = '/data1/S190425z/info/Initial/S190425z_Initial-all_candi.txt'
path_cand_upd = '/data1/S190425z/info/Update/S190425z_Update-all_candi.txt'
#	+KMTNet WITH OBSERVATORIES
path_phot_all = '/data1/S190425z/1.result/phot-S190425z.dat'
#------------------------------------------------------------
path_save = '/data1/S190425z/1.result/figure'
#------------------------------------------------------------
t0 = Time('2019-04-25T08:18:05.017', format='isot', scale='utc')
t0_upd = Time('2019-04-26T15:32:37', format='isot', scale='utc')
iphotbl = ascii.read(path_phot_ini)
uphotbl = ascii.read(path_phot_upd)
icantbl = ascii.read(path_cand_ini)
ucantbl = ascii.read(path_cand_upd)
tc_ini = SkyCoord(iphotbl['ra'], iphotbl['dec'], frame='icrs', unit='deg')
tc_upd = SkyCoord(uphotbl['ra'], uphotbl['dec'], frame='icrs', unit='deg')

photbl = ascii.read(path_phot_all)

#------------------------------------------------------------
#	PLOT SETTING FOR EACH OBSERVATORIES
#------------------------------------------------------------
'''
setting_SAO = dict(marker='+',
					markersize=10,
					color='hotpink',
					ls='',
					label='SAO')
setting_LOAO = dict(marker='+',
					markersize=10,
					color='yellowgreen',
					ls='',
					label='LOAO')
setting_LSGT = dict(marker='+',
					markersize=10,
					color='gold',
					ls='',
					label='LSGT')
setting_SQUEAN = dict(marker='+',
					markersize=10,
					color='purple',
					ls='',
					label='SQUEAN')
setting_KMTNET = dict(marker='+',
					markersize=10,
					color='dodgerblue',
					ls='',
					label='KMTNet')
setting_UKIRT = dict(marker='+',
					markersize=10,
					color='slategrey',
					ls='',
					label='UKIRT')
setting_obsplot = dict(kmtnet=setting_KMTNET,
					squean=setting_SQUEAN,
					ukirt=setting_UKIRT,
					loao=setting_LOAO,
					sao=setting_SAO,
					lsgt=setting_LSGT,)
'''

setting_SAO = dict(marker='+',
					markersize=10,
					color='hotpink',
					# alpha=0.5,
					ls='',
					label='SAO')
setting_LOAO = dict(marker='+',
					markersize=10,
					color='gold',
					# alpha=0.5,
					ls='',
					label='LOAO')
setting_LSGT = dict(marker='+',
					markersize=10,
					color='yellowgreen',
					# alpha=0.5,
					ls='',
					label='LSGT')
setting_SQUEAN = dict(marker='+',
					markersize=10,
					color='purple',
					# alpha=0.5,
					ls='',
					label='SQUEAN')
#------------------------------------------------------------
#	KMTNet
setting_sso = dict(marker='+',
					markersize=10,
					color='teal',
					alpha=0.5,
					ls='',
					label='KMTNet-SSO')
setting_saao = dict(marker='+',
					markersize=10,
					color='cyan',
					alpha=0.5,
					ls='',
					label='KMTNet-SAAO')
setting_ctio = dict(marker='+',
					markersize=10,
					color='deepskyblue',
					alpha=0.5,
					ls='',
					label='KMTNet-CTIO')
#------------------------------------------------------------
setting_UKIRT = dict(marker='+',
					markersize=10,
					color='slategrey',
					# alpha=0.5,
					ls='',
					label='UKIRT')
setting_obsplot = dict(sso=setting_sso,
					saao=setting_saao,
					ctio=setting_ctio,
					squean=setting_SQUEAN,
					ukirt=setting_UKIRT,
					loao=setting_LOAO,
					sao=setting_SAO,
					lsgt=setting_LSGT,)
#------------------------------------------------------------
#	SKYMAP
#------------------------------------------------------------
'''
plt.close('all')
skymap(bayestarfits=path_initial, photbl=iphotbl[iphotbl['jd']<t0_upd.jd])
plt.savefig(path_save+'/S190425z_Initial_skymap.png', overwrite=True)
'''
plt.close('all')
skymap(bayestarfits=path_update, photbl=photbl)
plt.savefig(path_save+'/S190425z_Update_skymap.png', overwrite=True)
'''
#------------------------------------------------------------
#	CUMULATIVE SCORE
#------------------------------------------------------------
plt.close('all')
cumscore(obslist=['kmtnet', 'loao', 'lsgt', 'sao'],
		setting_obsplot=setting_obsplot,
		photbl=uphotbl, timestep=1.1, t0=t0)
plt.savefig(path_save+'/S190425z_ul_Rband_LC.png', overwrite=True)
#------------------------------------------------------------
plt.close('all')
cumscore(obslist=['squean'],
		setting_obsplot=setting_obsplot,
		photbl=uphotbl, timestep=1.1, t0=t0)
plt.savefig(path_save+'/S190425z_ul_iband_LC.png', overwrite=True)
#------------------------------------------------------------
plt.close('all')
cumscore(obslist=['ukirt'],
         setting_obsplot=setting_obsplot,
         photbl=uphotbl, timestep=1.1, t0=t0)
plt.savefig(path_save+'/S190425z_ul_JHKband_LC.png', overwrite=True)
'''

'''
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
'''

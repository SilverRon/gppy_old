#	S190425z MODEL LC
#	Kasen+17
#	Tanaka+??
#	GW170817 LC
#------------------------------------------------------------
#	2019.08.21	CREATED BY Gregory S.H. Paek
#============================================================
from matplotlib import pyplot as plt
import ligo.skymap.plot
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy.io import ascii
from astropy.time import Time
#============================================================
#	FUNCTION
#============================================================
def asciiread(path_ascii):
	intbl = ascii.read(path_ascii)
	intbl.meta['path'] = path_ascii
	intbl.meta['name'] = os.path.basename(path_ascii)
	return intbl
#------------------------------------------------------------
def splinetable(gwtbl, filterlist):
	tblist = []
	for filte in filterlist:
		tmptbl = gwtbl[(gwtbl['Band']==filte) & (gwtbl['l_Mag']!='>')]
		tblist.append(tmptbl)
	return tblist
#------------------------------------------------------------
def kasen2table(kasentbl, filte, dist0, dist0er):
	kasentbl = kasentbl[kasentbl['day']>=0.0]
	appmag, appmagerr = abs2app(kasentbl[filte], np.zeros(len(kasentbl)), dist0, dist0er)
	filtes = np.array([filte for i in range(len(kasentbl))])
	newtbl = Table([np.copy(kasentbl['day']), filtes, appmag, appmagerr],
					names=('Phase', 'Band','Mag0', 'e_Mag0'))
	return newtbl
#------------------------------------------------------------
def kasen2table2(kasentbl, filte, dist, disterr, dist0, dist0err):
	kasentbl = kasentbl[kasentbl['day']>=0.0]
	# appmag, appmagerr = appscale(kasentbl[filte], np.zeros(len(kasentbl)), dist0, dist0er, dist)
	appmag, appmagerr = app2distscaling(kasentbl[filte], np.zeros(len(kasentbl)), dist, disterr, dist0, dist0err)
	filtes = np.array([filte for i in range(len(kasentbl))])
	newtbl = Table([np.copy(kasentbl['day']), filtes, appmag, appmagerr],
					names=('Phase', 'Band','Mag0', 'e_Mag0'))
	return newtbl
#------------------------------------------------------------
def abs2app(mag, magerr, dist, disterr):
	appmag = mag+5*np.log10(dist)-5
	appmagerr = np.sqrt(magerr**2 + ((5*disterr)/(dist*np.log(10)))**2)
	return appmag, appmagerr
#------------------------------------------------------------
def appscale(mag, magerr, dist0, dist0err, dist):
	'''
	dist0	: GW ?
	dist0err: GW ?
	dist	: GW 170817
	'''
	appmag = mag+5*np.log10(dist0/dist)
	appmagerr = np.sqrt(magerr**2 + ((5*dist0err)/(dist*np.log(10)))**2)
	return appmag, appmagerr
#------------------------------------------------------------
def app2distscaling(mag, magerr, dist, disterr, dist0, dist0err):
	'''
	mag		: APPARENT MAGNITUE
	magerr	: ERROR
	dist	: DISTANCE
	disterr	: DISTANCE ERROR
	dist0	: DISTANCE TO SCALE
	dist0err: DISTANCE ERROR TO SCALE
	'''
	mag0 = mag - 5*np.log10(dist/dist0)
	mag0err = np.sqrt((magerr)**2 + ((5*disterr)/(dist*np.log(10)))
					** 2 + ((5*dist0err)/(dist0*np.log(10)))**2)
	return mag0, mag0err
#------------------------------------------------------------
def drawLC(obstbl, obs, param_obs):
	plotbl = obstbl[obstbl['obs']==obs.upper()]
	param_plot = dict(	capsize=10, capthick=1,
						marker='v', ms=20, mec='k', ecolor='k',
						linestyle='', alpha=1.0)
	if 'kmtnet-' in obs:
		observat = obs.split('-')[1]
		plt.errorbar(plotbl['Phase'], plotbl['ul'], yerr=plotbl['ulstd'], **param_obs[observat.lower()], **param_plot)
	else:
		plt.errorbar(plotbl['Phase'], plotbl['ul'], yerr=plotbl['ulstd'], **param_obs[obs.lower()], **param_plot)
#------------------------------------------------------------
def drawLC2(obstbl, obs, param_obs):
	plotbl = obstbl[obstbl['obs']==obs.upper()]
	'''
	param_plot = dict(	capsize=10, capthick=1,
						marker='v', ms=20, 
						mec='k', ecolor='k',
						# mec=param_obs[obs]['mfc'], ecolor=param_obs[obs]['mfc'],
						linestyle='', alpha=0.75)
	'''
	if 'kmtnet-' in obs:
		observat = obs.split('-')[1]
		param_plot = dict(	capsize=10, capthick=1,
							marker='v', ms=20, 
							# mec='k', ecolor='k',
							mec=param_obs[observat]['mfc'], ecolor=param_obs[observat]['mfc'],
							linestyle='', alpha=0.75)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					yerr=[plotbl['yer1'], plotbl['yer2']], xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[observat.lower()], **param_plot)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[observat.lower()], **param_plot)
	else:
		param_plot = dict(	capsize=10, capthick=1,
							marker='v', ms=20, 
							# mec='k', ecolor='k',
							mec=param_obs[obs]['mfc'], ecolor=param_obs[obs]['mfc'],
							linestyle='', alpha=0.75)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					yerr=[plotbl['yer1'], plotbl['yer2']], xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[obs.lower()], **param_plot)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[obs.lower()], **param_plot)



#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_save = '.'
path_base = '/data1/S190425z/1.result/table'
path_gw170817 = path_base+'/GW170817_Villar+17_LC.dat'
# path_obs = '/data1/S190425z/1.result/Update/phot_upd_sources_GECKO.dat'
# path_obs = '/data1/S190425z/1.result/table/obs_all.dat'
# path_obs = '/data1/S190425z/1.result/table/obs_all_kmtnetsplit.dat'
path_obs = '/data1/S190425z/1.result/table/obs_complete_20191009.dat'
#------------------------------------------------------------
path_kasen_blue = path_base+'/knova_d1_n10_m0.025_vk0.30_Xlan1e-4.0.h5_z0.009787_mag.dat'
path_kasen_red = path_base+'/knova_d1_n10_m0.040_vk0.15_Xlan1e-1.5.h5_z0.009787_mag.dat'
#------------------------------------------------------------
path_piro = path_base+'/cocoon_Piro+18.dat'
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
dist = 37.7*1e6								#	[Mpc]
dister = 8.7*1e6							#	[Mpc]
# dist0 = 141*1e6								#	[Mpc]
# dist0er = 56*1e6							#	[Mpc]

dist0 = 200*1e6								#	[Mpc]
dist0er = 80*1e6							#	[Mpc]


t0 = Time('2020-02-13T04:10:40.330219', format='isot', scale='utc')
#------------------------------------------------------------
gwtbl = asciiread(path_gw170817)
#------------------------------------------------------------
kbltbl = asciiread(path_kasen_blue)
krdtbl = asciiread(path_kasen_red); krdtbl['day'] = krdtbl['t']
#------------------------------------------------------------
#	DIFFERENT EOS (TANAKA)
#------------------------------------------------------------

#------------------------------------------------------------
#	COCOON MODEL (Piro) - DIFFERENT FORMAT WITH OTHERS (VERTICAL)
#------------------------------------------------------------
ctbl = asciiread(path_piro)
#------------------------------------------------------------
#	GW170817 LC TABLE -> FOR EACH FILTERS(BANDS)
tblist = splinetable(gwtbl, ['r', 'R', 'i', 'J', 'Ks'])
comtbl = vstack(tblist)
comtbl.write(path_save+'/gw170817_LC.dat', format='ascii', overwrite=True)
comtbl = ascii.read(path_save+'/gw170817_LC.dat')
comtbl['Mag0'], comtbl['e_Mag0'] = app2distscaling(comtbl['Mag'], comtbl['e_Mag'], dist, dister, dist0, dist0er) 
comtbl.meta['path'] = path_gw170817
comtbl.meta['name'] = 'GW170817'
#------------------------------------------------------------
#	Kasen MODEL (BLUE & RED KN)
bltblist = []
rdtblist = []
for filte in ['r', 'i', 'J', 'Ks']:
	bltblist.append(kasen2table(kbltbl, filte, dist0, dist0er))
	rdtblist.append(kasen2table2(krdtbl, filte, dist, dister, dist0, dist0er))



bluetbl = vstack(bltblist)
redtbl = vstack(rdtblist)

# bluetbl.meta = kbltbl.meta
# redtbl.meta = krdtbl.meta
bluetbl.meta['name'] = 'Blue KN'
redtbl.meta['name'] = 'Red KN'
#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
band = 'J'
plt.close('all')
fig = plt.figure(figsize=(10, 8))
ax = plt.subplot(111)
#------------------------------------------------------------
#	KN MODEL PLOT
#------------------------------------------------------------
for intbl in [bluetbl, redtbl]:
	indx = np.where((intbl['Band']==band)&(intbl['Phase']<10.0))
	if len(indx[0]) == 0:
		indx = np.where((intbl['Band']==band+'s')&(intbl['Phase']<10.0))
	pltbl = intbl[indx]
	if intbl.meta['name'] == 'GW170817':
		ax.errorbar(intbl[indx]['Phase'], intbl[indx]['Mag0'], yerr=intbl[indx]['e_Mag0'],
					marker='o', color='grey', alpha=0.5,
					capsize=5, capthick=1,
					linestyle='None', label=intbl.meta['name'])
	if 'KN' in intbl.meta['name']:
		if 'Blue' in intbl.meta['name']:
			c, fc = 'blue', 'dodgerblue'
		elif 'Red' in intbl.meta['name']:
			c, fc = 'red', 'tomato'

		ax.plot(pltbl['Phase'], pltbl['Mag0']+pltbl['e_Mag0'],
				# y2=pltbl['Mag0']-pltbl['e_Mag0'],
				color=fc,
				alpha=0.5,
				)
		ax.plot(pltbl['Phase'], pltbl['Mag0']-pltbl['e_Mag0'],
				color=fc,
				alpha=0.5,
				)
		ax.fill_between(pltbl['Phase'],y1=pltbl['Mag0']+pltbl['e_Mag0'],
										y2=pltbl['Mag0']-pltbl['e_Mag0'],
										facecolor=fc,
										alpha=0.3,
										interpolate=True,
										label=intbl.meta['name'])
#------------------------------------------------------------
#	COCOON MODEL PLOT
#------------------------------------------------------------
ctbl0 = ctbl[(ctbl['filter']==band)&(ctbl['delmjd']<0.5)]
if band == 'r':
	ctbl0 = ctbl[(ctbl['filter']=='J')&(ctbl['delmjd']<3.0)]
cphase = ctbl0['delmjd']
cmag, cmager = abs2app(ctbl0['absmag'], 0, dist0, dist0er)
ax.plot(cphase, cmag+cmager,
		color='crimson',
		alpha=0.5,
		)
ax.plot(cphase, cmag-cmager,
		color='crimson',
		alpha=0.5,
		)
ax.fill_between(cphase, y1=cmag+cmager,
						y2=cmag-cmager,
						facecolor='crimson',
						alpha=0.3,
						interpolate=True,
						label='Cocoon')
# ax.plot(cphase, cmag, color='blue', alpha=0.5)

#------------------------------------------------------------
#	OBSERVATION
#------------------------------------------------------------
#	GW170817 DISCOVERY TIME
ax.axvline(x=0.48, color='grey', alpha=0.5, linestyle='--')
#------------------------------------------------------------
'''
chartBox = ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
ax.legend(loc='upper center', bbox_to_anchor=(1.15, 1.00), shadow=False, ncol=1, fontsize=15, frameon=False)
'''
#------------------------------------------------------------
plt.gca().invert_yaxis()
plt.xlim([-0.5, 8])
# plt.ylim([26.5, 18.5])
plt.xlabel('Phase [days]', fontsize=20)
plt.ylabel('Apparent magnitude', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(np.arange(15, 30, 1), fontsize=20)
# plt.legend(fontsize=15, loc='upper right')
plt.legend(fontsize=20, loc='lower right', framealpha=1.0)
plt.minorticks_on()
plt.title('{}-band'.format(band), fontsize=20)
#------------------------------------------------------------
ax0 = ax.twinx()
ax.tick_params(which='both', direction='in', labelsize=15)
ax.minorticks_on()
ax0.tick_params(which='both', direction='in', labelsize=15)
ax0.set_ylim(29-5*np.log10(dist0)+5, 18-5*np.log10(dist0)+5)
ax0.minorticks_on()
ax0.set_ylabel('Absolute magnitude', fontsize=20, rotation=270, labelpad=20)
#------------------------------------------------------------
plt.tight_layout()
plt.savefig('{}/S200213t_model_in{}band.png'.format(path_save, band), dpi=500, overwrite=True)
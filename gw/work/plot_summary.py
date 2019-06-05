#============================================================
#	THIS IS FOR COLLECTING ALL OBSERVATION DATA FOR FOLLOW-UP OBSERVATION TO S190425z
#	2019.05.23	MADE BY Gregory S.H. Paek
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
#	INPUT
#------------------------------------------------------------
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/result'
gwdist		= 156.1439917314773
gwdiststd	= 41.3733570000998
t0			= Time('2019-04-25T08:18:05.017553', format='isot', scale='utc')
jd0_170817	= 2457983.02852
#------------------------------------------------------------
#	LC REFERENCE
#------------------------------------------------------------
santos17	= ascii.read(path_base+'/LC_GW170817_Soares-Santos+17.dat')
santos17['jd'] = Time(santos17['mjd'], format='mjd', scale='utc').jd
santos17['deljd'] = santos17['jd']-jd0_170817
santos_r	= santos17[santos17['band']=='r'] 
santos_i	= santos17[santos17['band']=='i']
#------------------------------------------------------------
obsul_tbl	= ascii.read(path_base+'/obs_ul.dat')
kmtnet_tbl	= ascii.read(path_base+'/phot_kmtnet.dat')
loao_tbl	= ascii.read(path_base+'/phot_loao.dat')
lsgt_tbl	= ascii.read(path_base+'/phot_lsgt.dat')
sao_tbl		= ascii.read(path_base+'/phot_sao.dat')
squean_tbl	= ascii.read(path_base+'/phot_squean.dat')
#ukirt_tbl	= ascii.read('phot_ukirt.dat)
tblist		= [kmtnet_tbl, loao_tbl, lsgt_tbl, sao_tbl, squean_tbl]
#------------------------------------------------------------
tmptbl		= vstack(tblist)
obslist, objlist, bandlist	= [], [], []
for inim in tmptbl['obs']:
	part	= inim.split('-')
	obslist.append(part[1])
	objlist.append(part[2])
	bandlist.append(part[5])
obslist			= np.array(obslist)
objlist			= np.array(objlist)
bandlist		= np.array(bandlist)
#------------------------------------------------------------
#	MAKE COMPLETE TABLE
#------------------------------------------------------------
comtbl			= Table()
comtbl['object'], comtbl['obs'], comtbl['band']	= objlist, obslist, bandlist
comtbl['date-obs']			= tmptbl['date-obs']
t = Time(comtbl['date-obs'], format='isot', scale='utc')
comtbl['jd']				= t.jd
comtbl['seeing'], comtbl['mag']					= tmptbl['seeing'], tmptbl['mag']
comtbl.write(path_base+'phot_complete.dat', format='ascii', overwrite=True)
#------------------------------------------------------------
#	PLOT 1	: TIME - MAG.
#------------------------------------------------------------
#	https://partrita.github.io/posts/matplotlib-examples/
#	Get the figure and the axes
fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(8, 6))

jdrange		= np.arange(0, 10, 0.01)
'''
nsantos_r= func_linear(-0.7, jdrange, scaling=[santos_r['deljd'][0], santos_r['mag'][0]])
#func_santos_r	= interp1d(santos_r['deljd'][santos_r['deljd']<4], santos_r['mag'][santos_r['deljd']<4])
rmag, rmagerr = calc_app(nsantos_r, np.zeros(len(nsantos_r)), 38.4, 8.9, gwdist, gwdiststd)
ax0.plot(jdrange, rmag, color='blue', alpha=0.5)
ax0.fill_between(jdrange, rmag-rmagerr, rmag+rmagerr, color='dodgerblue', alpha=0.5)
'''
#------------------------------------------------------------

nsantos_i= func_linear(-0.5, jdrange, scaling=[santos_i['deljd'][0], santos_i['mag'][0]])
#func_santos_r	= interp1d(santos_r['deljd'][santos_r['deljd']<4], santos_r['mag'][santos_r['deljd']<4])
imag, imagerr = calc_app(nsantos_i, np.zeros(len(nsantos_i)), 38.4, 8.9, gwdist, gwdiststd)
ax0.plot(jdrange, imag, color='blue', alpha=0.5)
ax0.fill_between(jdrange, imag-imagerr, imag+imagerr, color='dodgerblue', alpha=0.5, label='i-band')



'''
#	THE FIRST GRAPH
ax0.plot(comtbl['jd']-t0.jd, comtbl['mag'], marker='o')
#ax0.set_ylim([2,20])
ax0.set(title='S190425z', xlabel='Time (Days from merger)', ylabel='Apparent Magnitude')
ax0.set_ylim([24, 17])
'''
#	KMTNet, LOAO, LSGT, SAO, SQUEAN(, UKIRT)
colorlist	= ['tomato', 'tomato', 'tomato', 'tomato', 'violet']
markerlist	= ['D', 'o', 's', '*', 'v']
size		= [10, 10, 10, 100, 10]
params_obs	= {		'KMTNET':'tomato-D-50',
					'LOAO'	:'tomato-o-50',
					'LSGT'	:'tomato-s-50',
					'SAO'	:'tomato-*-50',
					'SQUEAN':'violet-v-50'}
for obs in ['KMTNET', 'LOAO', 'LSGT', 'SAO', 'SQUEAN']:
	subtbl	= comtbl[comtbl['obs']==obs]
	sel		= params_obs[obs].split('-')
	color, marker, size	= sel[0], sel[1], int(sel[2])
	params_plot	= dict(	x=subtbl['jd']-t0.jd, y=subtbl['mag'],
						c=color, s=size, marker=marker,
						label=obs)#,
						#alpha=0.5)
	ax0.scatter(**params_plot)
#------------------------------------------------------------

#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
#fig.suptitle('S190425z', fontsize=14, fontweight='bold')
ax0.set(title='S190425z', xlabel='Time (Days from merger)', ylabel='Apparent Magnitude')#, fontsize=14)
ax0.set_ylim([24, 18])
ax0.set_xlim([0,3.5])
ax0.legend()
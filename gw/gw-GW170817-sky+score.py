#============================================================
#	PLOT GW170817 LOCALIZATION MAP ON THE SKY
#	REFERENCE	:	http://learn.astropy.org/plot-catalog.html
#	2019.06.03	MADE BY Gregory S.H. Paek
#============================================================
from astropy.io import ascii 
import numpy as np 
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table, vstack, Column
from matplotlib import cm
#============================================================
#	FUNCTION
#------------------------------------------------------------
def skyplot(ax, fontsize=20):
	param_healpix	= dict(	x=ra,y=dec,
							s=10,
							c=healtbl['P_2D'],cmap='OrRd')#,alpha=0.5)
	ax.scatter(**param_healpix)
	ax.locator_params(nbins=3)
	ax.set_xlabel('RA [deg]', fontsize=fontsize)
	ax.set_ylabel('Dec. [deg]', fontsize=fontsize)
	#ax.set_title('Title', fontsize=fontsize)
	ax.set_xlim(ax.get_xlim()[::-1])
	return 	ax.scatter(**param_healpix)
#------------------------------------------------------------
def candiplot(ax, fontsize=20):
	param_candi		= dict(	x=tra,y=tdec,
							s=50,
							alpha=0.75,
							c=cantbl['score'],cmap='winter')
	ax.scatter(**param_candi)
	ax.set_xlabel('RA [deg]', fontsize=fontsize)
	ax.set_ylabel('Dec. [deg]', fontsize=fontsize)
	ax.set_xlim(ax.get_xlim()[::-1])
	return  ax.scatter(**param_candi)
def allcandiplot(ax, fontsize=20):
	param_candi		= dict(	x=tra0,y=tdec0,
							s=10,
							alpha=0.25,
							c='grey')
	ax.scatter(**param_candi)
#------------------------------------------------------------
def scorebarplot(ax, fontsize=20):
	#xticklist	= cantbl['name']
	#xticklist	= ['NGC4970', 'NGC4993', 'IC4197']
	xticklist	= []
	param_score		= dict(	x=cantbl['name'],height=cantbl['score'],
							align='center',
							color='dodgerblue',alpha=0.75)
	ax.bar(**param_score)
	ax.set_ylabel('Score', fontsize=fontsize)
	ax.set_xticklabels(xticklist, rotation=45)
	#ax.set_ylim([0,0.25])
#------------------------------------------------------------
def cumplot(ax, fontsize=20):
	'''
	param_cum		= dict(	x=np.arange(0, len(cantbl), 1),y=np.cumsum(cantbl['score']),)
							marker='+',ms=10,mfc='dodgerblue',
							linestyle='--',color='tomato')
	'''
	xticklist	= []
	for i in range(len(cantbl)): xticklist.append('')
	param_cum		= dict(	marker='+',ms=20,mew=2,mec='dodgerblue',
							linestyle='--',linewidth=3,color='tomato')
	ax.plot( cantbl['name'], np.cumsum(cantbl['score']), **param_cum)
	ax.set_ylabel('Score', fontsize=fontsize)
	ax.set_xticklabels(xticklist)

#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
path_save	= '/home/sonic/S190425z'
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/gw170817/GW170817/K90'
healtbl = ascii.read(path_base+'/GW170817_healpix.dat')
#	OVERPLOT WITH CANDIDATES
cantbl	= ascii.read(path_base+'/GW170817-newname_candi.txt')
outbl	= ascii.read(path_base+'/candi_inconf90.dat')
#------------------------------------------------------------
#	COORDINATE
#------------------------------------------------------------
#	HEALPIX
ra = coord.Angle(healtbl['ra']*u.degree)
#ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(healtbl['dec']*u.degree)
#	CANDIDATES
tra = coord.Angle(cantbl['ra']*u.degree)
#tra = tra.wrap_at(180*u.degree)
tdec = coord.Angle(cantbl['dec']*u.degree)
#	ALL CANDIDATES
tra0	= coord.Angle(outbl['ra']*u.degree)
#tra0	= tra0.wrap_at(180*u.degree)
tdec0	= coord.Angle(outbl['dec']*u.degree)
#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 20})

fig = plt.figure(figsize=(16, 12))

ax1 = plt.subplot(221)
plt.minorticks_on()
#plt.xticks(rotation=90)
ax2 = plt.subplot(223)
plt.minorticks_on()
plt.xticks(rotation=90)
ax3 = plt.subplot(122)
plt.minorticks_on()


cb1 = skyplot(ax3)

allcandiplot(ax3)

cb2	= candiplot(ax3)
scorebarplot(ax1)
cumplot(ax2)


fig.colorbar(cb2, ax=ax3)
#fig.colorbar(cb1, ax=ax3)
fig.tight_layout()
plt.minorticks_on()
fig.savefig(path_save+'/Figure_X+GW170817_sky+score.png', overwrite=True)

'''
#------------------------------------------------------------
#	PLOT FOR 
#------------------------------------------------------------
plt.close('all')
fig = plt.figure(figsize=(10,10))
bins = np.arange(0, 120, 5)

param_bkg	= dict(	color='grey', alpha=0.5,
					label='All candidates ({}) within 90% confidecne region'.format(len(outbl)))


param_sel	= dict(	color='dodgerblue', alpha=0.75,
					label='Selected candidates ({}) within distance range'.format(len(cantbl)))

plt.hist(outbl['dist'], bins=bins, **param_bkg)
plt.hist(cantbl['dist'], bins=bins, **param_sel)   
plt.legend(prop={'size':18}, loc='upper left')
plt.xlim([0, 100])
plt.ylim([0, 25])
#------------------------------------------------------------
plt.minorticks_on()
plt.tick_params(labelsize='20', length=7.5)
plt.xlabel('Distance [Mpc]', size=20)
fig.tight_layout()
plt.savefig(path_save+'/Figure_X+GW170817_distance_range.png', overwrite=True)

'''
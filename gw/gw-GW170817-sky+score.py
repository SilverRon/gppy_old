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
def skyplot(ax, fontsize=16):
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
def candiplot(ax, fontsize=16):
	param_candi		= dict(	x=tra,y=tdec,
							s=50,
							alpha=0.75,
							c=cantbl['score'],cmap='winter')
	ax.scatter(**param_candi)
	ax.set_xlabel('RA [deg]', fontsize=fontsize)
	ax.set_ylabel('Dec. [deg]', fontsize=fontsize)
	ax.set_xlim(ax.get_xlim()[::-1])
	return  ax.scatter(**param_candi)
#------------------------------------------------------------
def scorebarplot(ax, fontsize=16):
	xticklist	= []
	for i in range(len(cantbl)): xticklist.append('')
	param_score		= dict(	x=cantbl['name'],height=cantbl['score'],
							align='center',
							color='dodgerblue',alpha=0.75)
	ax.bar(**param_score)
	ax.set_ylabel('Score', fontsize=fontsize)
	ax.set_xticklabels(xticklist)
	#ax.set_ylim([0,0.25])
#------------------------------------------------------------
def cumplot(ax, fontsize=16):
	'''
	param_cum		= dict(	x=np.arange(0, len(cantbl), 1),y=np.cumsum(cantbl['score']),)
							marker='+',ms=10,mfc='dodgerblue',
							linestyle='--',color='tomato')
	'''
	param_cum		= dict(	marker='+',ms=20,mew=2,mec='dodgerblue',
							linestyle='--',linewidth=3,color='tomato')
	ax.plot( cantbl['name'], np.cumsum(cantbl['score']), **param_cum)
	ax.set_ylabel('Score', fontsize=fontsize)
	ax.set_xticklabels(cantbl['name'])
#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/gw170817/GW170817'
healtbl = ascii.read(path_base+'/GW170817_healpix.dat')
#	OVERPLOT WITH CANDIDATES
cantbl	= ascii.read(path_base+'/GW170817-newname_candi.txt')
#------------------------------------------------------------
#	COORDINATE
#------------------------------------------------------------
#	HEALPIX
ra = coord.Angle(healtbl['ra']*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(healtbl['dec']*u.degree)
#	CANDIDATES
tra = coord.Angle(cantbl['ra']*u.degree)
tra = tra.wrap_at(180*u.degree)
tdec = coord.Angle(cantbl['dec']*u.degree)
#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
plt.close('all')
fig = plt.figure(figsize=(16,12))

ax1 = plt.subplot(221)
#plt.minorticks_on()
#plt.xticks(rotation=90)
ax2 = plt.subplot(223)
#plt.minorticks_on()
plt.xticks(rotation=90)
ax3 = plt.subplot(122)
plt.minorticks_on()

cb1 = skyplot(ax3)
cb2	= candiplot(ax3)
scorebarplot(ax1)
cumplot(ax2)


fig.colorbar(cb2, ax=ax3)
fig.colorbar(cb1, ax=ax3)
fig.tight_layout()

#============================================================
import numpy as np
from astropy.io import ascii
from astropy.table import vstack
from astropy import units as u
import matplotlib
import matplotlib.pyplot as plt
#============================================================
#	FUNCTION
#------------------------------------------------------------
def tilemaker(cantbl, ra, dec, fov):
	decfactor	= np.cos(dec*u.deg).value
	nfov		= fov/decfactor
	ramin, ramax= ra-nfov/2, ra+nfov/2
	demin, demax= dec-fov/2, dec+fov/2
	indx		= np.where(	(cantbl['ra']>ra-nfov/2) & (cantbl['ra']<ra+nfov/2) &
							(cantbl['dec']>dec-fov/2) & (cantbl['dec']<dec+fov/2))
	return indx, ramin, ramax, demin, demax
#------------------------------------------------------------
def skyplot(ax, healtbl, fontsize=16):
	param_healpix	= dict(	x=healtbl['ra'],y=healtbl['dec'],
							s=10,
							c=healtbl['P_2D'],cmap='OrRd',
							label='_nolegend_')#,alpha=0.5)
	ax.scatter(**param_healpix)
	ax.locator_params(nbins=3)
	ax.set_xlim(ax.get_xlim()[::-1])
	return 	ax.scatter(**param_healpix)
#------------------------------------------------------------
def candiplot(ax, cantbl):
	tra, tdec	= cantbl['ra'], cantbl['dec']
	param_candi		= dict(	x=tra, y=tdec,
							s=5,
							alpha=1,c=cantbl['score'],
							label='Candidates ({})'.format(len(cantbl)))#,cmap='binary')
	ax.scatter(**param_candi)
#------------------------------------------------------------
def cumplot(ax, cantbl, fontsize=16):
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
'''
def findchangept(cantbl):
	delx		= 0.01
	x0, y0		= np.arange(len(cantbl)), np.cumsum(np.copy(cantbl['score']))
	x1, y1		= [], []
	for i in range(len(cantbl)-1):
		slope	= (y0[i+1]-y0[i])/(x0[i+1]-x0[i])
		x1.append(i)
		y1.append(slope)
	x2, y2		= [], []
	for i in range(len(x1)-1):
		slope	= (y1[i+1]-y1[i])/(x1[i+1]-x1[i])
		x2.append(i)
		y2.append(slope)
'''
#------------------------------------------------------------
def highpt4half(cantbl):
	x0, y0		= np.arange(len(cantbl)), np.cumsum(np.copy(cantbl['score']))
	return x0[y0<=0.9]
#------------------------------------------------------------
def eachobjs(objs):
	objlist	= []
	for obj in objs:
		objlist.append(obj)
	objlist	= list(set(objlist))
	return objlist
#------------------------------------------------------------
def counthigh(intbl, cantbl):
	hightbl	= cantbl[highpt4half(cantbl)]
	objlist		= []
	scorelist	= []
	for obj in eachobjs(hightbl['name']):
		indx	= np.where(obj == intbl['name'])
		if len(indx[0])!=0:
			score	= intbl['score'][indx]
			objlist.append(obj)
			scorelist.append(np.copy(score)[0])
	highnumb	= len(objlist)
	highscore	= np.sum(scorelist)
	highobjs	= objlist
	return highnumb, highscore, highobjs
#============================================================
#	INPUT
#============================================================
path_save		= '/home/sonic/S190425z'
#path_base		= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/20190609/jkim'
path_base = '/data0/Users/User/Downloads/data/Project/gw/S190425z/20190609/jkim'
initbl, updtbl	= ascii.read(path_base+'/gwfollowup_kmtnet_target_G330561_initial.txt'), ascii.read(path_base+'/gwfollowup_kmtnet_target_G330561_update.txt')
icantbl, ucantbl	= ascii.read('S190425z_Initial-all_candidates.txt'), ascii.read('S190425z_Update-all_candidates.txt')
ihealtbl, uhealtbl	= ascii.read(path_base+'/S190425z_Initial_healpix.dat'), ascii.read(path_base+'/S190425z_Update_healpix.dat')
#------------------------------------------------------------
itblist	= []
for i in range(len(initbl)): 
	try: 
		n = 2*i+1 
		itblist.append(initbl[n]) 
	except: 
		pass 
icomtbl = vstack(itblist)

utblist	= []
for i in range(len(updtbl)): 
	try: 
		n = 2*i+1 
		utblist.append(updtbl[n]) 
	except: 
		pass 
ucomtbl = vstack(utblist)
#------------------------------------------------------------
#	PLOT FOR INITIAL
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 16})
fig, ax	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(12, 6))

skyplot(ax, ihealtbl)
icanlist= []
itile	= []
for i in range(len(icomtbl[:200])):
	ra, dec	= icomtbl['radeg'][i], icomtbl['decdeg'][i]
	indx, ramin, ramax, demin, demax	= tilemaker(icantbl, ra, dec, 2)
	rect	= matplotlib.patches.Rectangle(	(ramin, demin), ramax-ramin, demax-demin,
											color='grey', alpha=0.3, label='_nolegend_')
	ax.add_patch(rect)
	itile.append((ramax-ramin)*(demax-demin))
	if len(indx[0])!=0:	icanlist.append(icantbl[indx])
intbl_init	= vstack(icanlist)
intbl_init	= intbl_init[intbl_init['score'].argsort()]

candiplot(ax, intbl_init)
#------------------------------------------------------------
#	LABEL
rect	= matplotlib.patches.Rectangle(	(-999, -999), 1, 1,
										color='grey', alpha=0.5, label='KMTNet Tile')
ax.add_patch(rect)



ax.set_xlabel('RA [deg]', fontsize=16)
ax.set_ylabel('Dec. [deg]', fontsize=16)
#ax.set_title('Title', fontsize=fontsize)
ax.set_ylim([-40, +40])
ax.set_xlim([300,160])
#plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM counterpart discovery')
param_legend	= dict(	loc='upper right', fontsize=20,
						fancybox=True, framealpha=1.0,
						scatterpoints=1, markerscale=5)
ax.legend(**param_legend)
plt.tight_layout()
plt.minorticks_on()
plt.savefig(path_save+'/Figure_X+KMTNet_Tile_Initial.png', overwrite=True)
#------------------------------------------------------------
ihighnumb, ihighscore, ihighobjs = counthigh(intbl_init, icantbl)
print('='*55+'\n\t\tInitial\n'+'='*55)
print('COVERED REGION BY KMTNet TILE ({})\t{} deg2'.format(len(itile), round(np.sum(itile), 3)))
print('NUMBER OF HIGH PRIORITY CANDIDATES\t{}/{}'.format(ihighnumb, len(intbl_init)))
print('HIGH PRIORITY SCORE PER TOTAL SCORE\t{}/{}'.format(round(ihighscore, 3), round(np.sum(intbl_init['score']), 3)))
#------------------------------------------------------------
#	PLOT FOR UPDATE
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 16})
fig, ax	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(12, 6))

skyplot(ax, uhealtbl)
ucanlist= []
utile	= []
for i in range(len(ucomtbl[:200])):
	ra, dec	= ucomtbl['radeg'][i], ucomtbl['decdeg'][i]
	indx, ramin, ramax, demin, demax	= tilemaker(ucantbl, ra, dec, 2)
	rect	= matplotlib.patches.Rectangle(	(ramin, demin), ramax-ramin, demax-demin,
											color='grey', alpha=0.3, label='_nolegend_')
	ax.add_patch(rect)
	utile.append((ramax-ramin)*(demax-demin))
	if len(indx[0])!=0:	ucanlist.append(ucantbl[indx])
intbl_upd	= vstack(ucanlist)
intbl_upd	= intbl_upd[intbl_upd['score'].argsort()]

candiplot(ax, intbl_upd)
#------------------------------------------------------------
#	LABEL
rect	= matplotlib.patches.Rectangle(	(-999, -999), 1, 1,
										color='grey', alpha=0.5, label='KMTNet Tile')
ax.add_patch(rect)



ax.set_xlabel('RA [deg]', fontsize=16)
ax.set_ylabel('Dec. [deg]', fontsize=16)
#ax.set_title('Title', fontsize=fontsize)
#ax0.set_ylim([24, 18])
#ax0.set_xlim([0,3])
#plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM counterpart discovery')
param_legend	= dict(	loc='upper right', fontsize=20,
						fancybox=True, framealpha=1.0,
						scatterpoints=1, markerscale=5)
ax.legend(**param_legend)
plt.tight_layout()
plt.minorticks_on()
plt.savefig(path_save+'/Figure_X+KMTNet_Tile_Update.png', overwrite=True)
#------------------------------------------------------------
uhighnumb, uhighscore, uhighobjs = counthigh(intbl_upd, ucantbl)
print('='*55+'\n\t\tUpdate\n'+'='*55)
print('COVERED REGION BY KMTNet TILE ({})\t{} deg2'.format(len(utile), round(np.sum(utile), 3)))
print('NUMBER OF HIGH PRIORITY CANDIDATES\t{}/{}'.format(uhighnumb, len(intbl_upd)))
print('HIGH PRIORITY SCORE PER TOTAL SCORE\t{}/{}'.format(round(uhighscore, 3), round(np.sum(intbl_upd['score']), 3)))

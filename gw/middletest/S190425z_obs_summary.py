#============================================================
#	
#	REFERENCE	:	http://learn.astropy.org/plot-catalog.html
#	2019.05.30	MADE BY Gregory S.H. Paek
#============================================================
from astropy.io import ascii 
import numpy as np 
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table, vstack
#============================================================
#	FUNCTION
#------------------------------------------------------------
def matching(intbl, reftbl, inrakey='ra',indeckey='dec',refrakey='ra',refdeckey='dec',sep=2.0):
	"""
	MATCHING TWO CATALOG WITH RA, Dec COORD. WITH python
	INPUT   :   SE catalog, SDSS catalog file name, sepertation [arcsec]
	OUTPUT  :   MATCED CATALOG FILE & TABLE
	"""
	import numpy as np
	import astropy.units as u
	from astropy.table import Table, Column
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii

	coo_intbl   = SkyCoord(intbl[inrakey], intbl[indeckey], unit=(u.deg, u.deg))
	coo_reftbl  = SkyCoord(reftbl[inrakey], reftbl[refdeckey], unit=(u.deg, u.deg))

	#   INDEX FOR REF.TABLE
	indx, d2d, d3d  = coo_intbl.match_to_catalog_sky(coo_reftbl)
	ref_match       = reftbl[indx]
	ref_match['sep']= d2d
	ref_match_col   = ref_match.colnames
	merge_tbl       = intbl
	for col in ref_match.colnames:
		merge_tbl[col]  = ref_match[col]
	indx_cut        = np.where(merge_tbl['sep']*3600. < sep)
	merge           = merge_tbl[indx_cut]
	return  merge
#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/info'
#healtbl = ascii.read(path_base+'/S190425z_Initial_healpix.dat')
healtbl = ascii.read(path_base+'/S190425z_Update_healpix.dat')
#cantbl	= ascii.read(path_base+'/S190425z_Initial-all_candidates_rank.txt')
cantbl	= ascii.read(path_base+'/S190425z_Update-all_candidates.txt')
#	REVERSE
#cantbl	= cantbl[(-1*cantbl['rank']).argsort()]#[-500:]
'''
obstbl	= ascii.read(path_base+'/resultphot_complete.dat')
#------------------------------------------------------------
#	TABLE MATCHING
#------------------------------------------------------------
#	CHANGE '-' -> '+' IN NAME
newname	= []
for obj in cantbl['name']:
	newname.append(obj.replace('-','+'))
cantbl['name']	= np.array(newname)
#
ralist, delist, Plist	= [], [], []
outlist					= []
for obj in obstbl['object']:
	indx	= np.where(obj == cantbl['name'])
	if len(indx[0]) != 0:
		ralist.append(cantbl['ra'][indx][0])
		delist.append(cantbl['dec'][indx][0])
		Plist.append(cantbl['score'][indx][0])
	if len(indx[0]) == 0:
		ralist.append(-99)
		delist.append(-99)
		Plist.append(-99)
obstbl['ra'], obstbl['dec'], obstbl['score']	= np.array(ralist), np.array(delist), np.array(Plist)
'''
#	intbl : MATCHED TABLE, outbl : NO MATCHED TABLE
obstbl	= ascii.read(path_base+'/S190425z_observed_complete.dat')
obstbl	= obstbl[(obstbl['score']).argsort()]
intbl	= obstbl[obstbl['score']!=-99]
outbl	= obstbl[obstbl['score']==-99]
#	CHANGE '-' -> '+' IN NAME
Plist	= []
outlist	= []
for obj in obstbl['object']:
	indx	= np.where(obj == cantbl['name'])
	if len(indx[0]) != 0:
		Plist.append(cantbl['score'][indx][0])
	if len(indx[0]) == 0:
		Plist.append(-99)
		outlist.append(obj)
obstbl['score']	= np.array(Plist)
#------------------------------------------------------------
#	COORDINATE
#------------------------------------------------------------
#	HEALPIX
ra = coord.Angle(healtbl['ra']*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(healtbl['dec']*u.degree)
#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
plt.rcParams.update({'font.size': 24})
#	HEALPIX RA&DEC, P_2D
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide")
param_cand	= dict(	x=ra.radian,y=dec.radian,
					s=10,
					c=healtbl['P_2D'],cmap='OrRd')#,alpha=0.5)
ax.scatter(**param_cand)
ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
ax.grid(True)
#	HOST GALAXY CANDIDATES (all)
'''
tra		= coord.Angle(cantbl['ra']*u.degree)
tra		= tra.wrap_at(180*u.degree)
tdec	= coord.Angle(cantbl['dec']*u.degree)
param_cand	= dict(	x=tra.radian,y=tdec.radian,
					marker='+',linewidths=0.5,s=50,
					c=cantbl['score'],cmap='winter',alpha=1)
ax.scatter(**param_cand)
'''
#	HOST GALAXY CANDIDATES (OBSERVATION)
'''
tra		= coord.Angle(intbl['ra']*u.degree)
tra		= tra.wrap_at(180*u.degree)
tdec	= coord.Angle(intbl['dec']*u.degree)
param_cand	= dict(	x=tra.radian,y=tdec.radian,
					marker='+',linewidths=0.5,s=100,
					c=intbl['score'],cmap='winter',alpha=1)
ax.scatter(**param_cand)

tra		= coord.Angle(obstbl['ra'][obstbl['ra']!=-99]*u.degree)
tra		= tra.wrap_at(180*u.degree)
tdec	= coord.Angle(obstbl['dec'][obstbl['ra']!=-99]*u.degree)
param_cand	= dict(	x=tra.radian,y=tdec.radian,
					marker='+',linewidths=0.5,s=100,
					c='',alpha=0.5)
ax.scatter(**param_cand)
'''
#------------------------------------------------------------
cutbl	= obstbl[obstbl['ra']!=-99]
#for obs in ['LOAO', 'SAO', 'SQUEAN', 'LSGT', 'UKIRT', 'KMTNET']:
for obs in ['LOAO', 'UKIRT', 'KMTNET']:
	if obs == 'LOAO':
		color='gold'
		label='LOAO'
	if obs == 'SAO':
		color='green'
		label='SAO'
	if obs == 'SQUEAN':
		color='purple'
		label='SQUEAN'
	if obs == 'LSGT':
		color='black'
		label='LSGT'
	if obs == 'UKIRT':
		color='tomato'
		label='UKIRT'
	if obs == 'KMTNET':
		color='dodgerblue'
		label='KMTNet'
	subtbl	= cutbl[cutbl['obs']==obs]
	tra		= coord.Angle(subtbl['ra'][subtbl['ra']!=-99]*u.degree)
	tra		= tra.wrap_at(180*u.degree)
	tdec	= coord.Angle(subtbl['dec'][subtbl['ra']!=-99]*u.degree)
	param_cand	= dict(	x=tra.radian,y=tdec.radian,
						marker='+',linewidths=1,s=200,
						color=color,alpha=0.75,
						label='{} ({})'.format(label, len(subtbl)))
	ax.scatter(**param_cand)
#------------------------------------------------------------
newtbl		= ascii.read('/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/result/retry/phot_all.dat')
for obs in ['SAO', 'SQUEAN', 'LSGT']:
	if obs == 'LOAO':
		color='gold'
		label='LOAO'
	if obs == 'SAO':
		color='green'
		label='SAO'
	if obs == 'SQUEAN':
		color='purple'
		label='SQUEAN'
	if obs == 'LSGT':
		color='black'
		label='LSGT'
	if obs == 'UKIRT':
		color='tomato'
		label='UKIRT'
	if obs == 'KMTNET':
		color='dodgerblue'
		label='KMTNet'
	subtbl	= newtbl[newtbl['obs']==obs]
	tra		= coord.Angle(subtbl['radeg'][subtbl['radeg']!=-99]*u.degree)
	tra		= tra.wrap_at(180*u.degree)
	tdec	= coord.Angle(subtbl['decdeg'][subtbl['radeg']!=-99]*u.degree)
	param_cand	= dict(	x=tra.radian,y=tdec.radian,
						marker='+',linewidths=1,s=200,
						color=color,alpha=0.75,
						label='{} ({})'.format(label, len(subtbl)))
	ax.scatter(**param_cand)
#------------------------------------------------------------
ax.set(title='S190425z')
ax.legend(loc='upper right')
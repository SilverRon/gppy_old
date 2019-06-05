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
from astropy.table import Table, vstack, Column
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
	coo_reftbl  = SkyCoord(reftbl[refrakey], reftbl[refdeckey], unit=(u.deg, u.deg))

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
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/result/retry'
healtbl = ascii.read(path_base+'/S190425z_Update_healpix.dat')
cantbl	= ascii.read(path_base+'/S190425z_Update-all_candidates.txt')
obstbl	= ascii.read(path_base+'/phot_all.dat')
obstbl['P_2D'], obstbl['score'], obstbl['sort']	= None, None, None
#------------------------------------------------------------
#	NAME MATCHING
tblist	= []
nolist	= []
for i in range(len(obstbl)):
	obj	= obstbl['object'][i]
	onetbl	= obstbl[i]
	if obj in cantbl['name']:
		indx	= np.where(obj == cantbl['name'])
		onetbl['P_2D'], onetbl['score'], onetbl['sort']	= \
		np.copy(cantbl['P_2D'][indx])[0], np.copy(cantbl['score'][indx])[0], np.copy(cantbl['sort'][indx])[0]
		tblist.append(onetbl)
	else:
		nolist.append(onetbl)
namtbl	= vstack(tblist)
namtbl['P_2D'], namtbl['score'], namtbl['sort']	= \
Column(np.copy(namtbl['P_2D']), dtype='float64'), Column(np.copy(namtbl['score']), dtype='float64'), Column(np.copy(namtbl['sort']), dtype='str')
notbl	= vstack(nolist)
#------------------------------------------------------------
#	RA, Dec MATCHING (sep ["])
param_matching	= dict(	intbl=notbl, reftbl=cantbl,
						inrakey='radeg', indeckey='decdeg',
						refrakey='ra', refdeckey='dec',
						sep=300.0)
cootbl	= matching(**param_matching)
cootbl['sort']	= Column(np.copy(cootbl['sort']), dtype='str')
rmkeys	= ['name',
 'dist',
 'dist_err',
 'z',
 'B',
 'B_err',
 'B_Abs',
 'J',
 'J_err',
 'H',
 'H_err',
 'K',
 'K_err',
 'flag1',
 'flag2',
 'flag3',
 'flag4',
 'sep',
 'ra',
 'dec']
for key in rmkeys:
	if key in cootbl.keys():
		del cootbl[key]
	if key in namtbl.keys():
		del namtbl[key]
#------------------------------------------------------------
comtbl	= vstack([namtbl, cootbl])
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
plt.rcParams.update({'font.size': 30})
#	HEALPIX RA&DEC, P_2D
try:
	plt.close('all')
except:
	pass
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide")
param_cand	= dict(	x=ra.radian,y=dec.radian,
					s=10,
					c=healtbl['P_2D'],cmap='OrRd')#,alpha=0.5)
ax.scatter(**param_cand)
ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
ax.grid(True)
#------------------------------------------------------------
#	HOST GALAXY CANDIDATES (all)
#for obs in ['LOAO', 'SAO', 'SQUEAN', 'LSGT', 'UKIRT', 'KMTNET']:
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
	subtbl	= obstbl[obstbl['obs']==obs]
	tra		= coord.Angle(subtbl['radeg']*u.degree)
	tra		= tra.wrap_at(180*u.degree)
	tdec	= coord.Angle(subtbl['decdeg']*u.degree)
	param_cand	= dict(	x=tra.radian,y=tdec.radian,
						marker='+',linewidths=1,s=200,
						color=color,alpha=0.75,
						label=label)
	ax.scatter(**param_cand)
#------------------------------------------------------------
ax.set(title='S190425z')
ax.legend()
plt.tight_layout()
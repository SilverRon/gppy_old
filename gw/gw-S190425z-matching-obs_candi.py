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
from astropy.coordinates import SkyCoord
#============================================================
def matching(intbl, reftbl, inrakey='radeg', indeckey='decdeg', refrakey='ra', refdeckey='dec', sep=180.0):
	mtbl		= intbl

	coo_intbl   = SkyCoord(intbl[inrakey], intbl[indeckey], unit=(u.deg, u.deg))
	coo_reftbl  = SkyCoord(reftbl[refrakey], reftbl[refdeckey], unit=(u.deg, u.deg))
	#   INDEX FOR REF.TABLE
	indx, d2d, d3d  = coo_intbl.match_to_catalog_sky(coo_reftbl)
	mreftbl			= reftbl[indx]
	mreftbl['sep']	= d2d
	'''
	for col	in mreftbl.keys():
		if col not in mtbl.keys():
			mtbl[col]	= mreftbl[col]
	'''

	for col	in mreftbl.keys():
		mtbl[col]	= mreftbl[col]

	indx_cut        = np.where(mtbl['sep']*3600. < sep)
	mtbl			= mtbl[indx_cut]
	return mtbl
#------------------------------------------------------------
def namematching(intbl, reftbl, innamekey='object', refnamekey='name'):
	tblist	= []
	nolist	= []
	for i in range(len(intbl)):
		obj = intbl[innamekey][i]
		check	= False
		for robj in reftbl[refnamekey]:
			if obj == robj.replace('-','+'):
				onetbl	= intbl[i]
				tmptbl	= Table()
				for col in intbl.keys():
					tmptbl[col]	= [onetbl[col]]
				for col in reftbl.keys():
					tmptbl[col] = reftbl[col][reftbl[refnamekey]==robj]
				tblist.append(tmptbl)
				check = True
				pass
		if check == False:
			nolist.append(intbl[i])
	notbl	= vstack(nolist)
	comtbl	= vstack(tblist)
	return comtbl, notbl
#------------------------------------------------------------
def candiplot(ax, ra, dec, score, fontsize=16):
	tra = coord.Angle(ra*u.degree)
	tra = tra.wrap_at(180*u.degree)
	tdec = coord.Angle(dec*u.degree)

	param_candi		= dict(	x=tra,y=tdec,
							s=50,
							alpha=0.75,
							c=score,cmap='winter')
	ax.scatter(**param_candi)
	ax.set_xlabel('RA [deg]', fontsize=fontsize)
	ax.set_ylabel('Dec. [deg]', fontsize=fontsize)
	ax.set_xlim(ax.get_xlim()[::-1])
	return  ax.scatter(**param_candi)


def mollplot(ax, ra, dec, score):
	tra = coord.Angle(ra*u.degree)
	tra = tra.wrap_at(180*u.degree)
	tdec = coord.Angle(dec*u.degree)
	param_cand	= dict(	x=tra.radian,y=tdec.radian,
						s=100, marker='+',
						c=score,cmap='winter')#,alpha=0.5)
	ax.scatter(**param_cand)
	ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
	ax.grid(True)


#------------------------------------------------------------
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/20190606'
photbl		= ascii.read(path_base+'/phot_all.dat')
#kmtbl		= ascii.read(path_base+'/phot_kmtnet.dat')
icantbl		= ascii.read(path_base+'/S190425z_Initial-all_candidates.txt')
ucantbl		= ascii.read(path_base+'/S190425z_Update-all_candidates.txt')
#------------------------------------------------------------
#	INITIAL SIGNAL
inametbl, inotbl= namematching(photbl, icantbl)
icootbl			= matching(inotbl, icantbl)
iphotbl			= vstack([inametbl, icootbl])
#------------------------------------------------------------
#	UPDATE SIGNAL
unametbl, unotbl= namematching(photbl, ucantbl)
ucootbl			= matching(unotbl, ucantbl)
uphotbl			= vstack([unametbl, ucootbl])
#------------------------------------------------------------
#	OVERLAPPED
#------------------------------------------------------------
#	INITIAL
objlist	= []
for obj in iphotbl['name']:
	objlist.append(obj)
objlist	= list(set(objlist))
objlist.sort()

tblist	= []
for obj in objlist:
	for obs in ['UKIRT', 'SQUEAN', 'LOAO', 'LSGT', 'SAO']:
		obstbl	= iphotbl[iphotbl['obs']==obs]
		indx	= np.where(obj == obstbl['name'])
		if		len(indx[0]) == 0:
			pass
		elif	len(indx[0]) == 1:
			tblist.append(obstbl[indx])
			break
		elif	len(indx[0]) > 2:
			tblist.append(obstbl[indx][0])
			break
iovertbl	= vstack(tblist)
#------------------------------------------------------------
#	UPDATE
objlist	= []
for obj in uphotbl['name']:
	objlist.append(obj)
objlist	= list(set(objlist))
objlist.sort()

tblist	= []
for obj in objlist:
	for obs in ['UKIRT', 'SQUEAN', 'LOAO', 'LSGT', 'SAO']:
		obstbl	= uphotbl[uphotbl['obs']==obs]
		indx	= np.where(obj == obstbl['name'])
		if		len(indx[0]) == 0:
			pass
		elif	len(indx[0]) == 1:
			tblist.append(obstbl[indx])
			break
		elif	len(indx[0]) > 2:
			tblist.append(obstbl[indx][0])
			break
uovertbl	= vstack(tblist)








#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
plt.close('all')
plt.rcParams.update({'font.size': 24})
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111, projection="mollweide")
mollplot(ax, iphotbl['ra'], iphotbl['dec'], iphotbl['score'])

mollplot(ax, uphotbl['ra'], uphotbl['dec'], uphotbl['score'])

#mollplot(ax, icantbl['ra'], icantbl['dec'], icantbl['score'])






#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
#	INTITIALIZE
plt.close('all')
fig = plt.figure(figsize=(16,12))
ax1 = plt.subplot(221)
ax2 = plt.subplot(223)
plt.xticks(rotation=90)
ax3 = plt.subplot(122)
plt.minorticks_on()
#------------------------------------------------------------
cb2	= candiplot(ax3, iphotbl['radeg'], iphotbl['decdeg'], iphotbl['score'])


fig.tight_layout()
fig.savefig('S190425z_sky.png', overwrite=True)
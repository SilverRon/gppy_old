#============================================================
#	PLOT GW LOCALIZATION MAP ON THE SKY
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

#	NONE

#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
path_base	= '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/info'
healtbl = ascii.read(path_base+'/S190425z_Update_healpix.dat')
'''
#	OVERPLOT WITH CANDIDATES
cantbl	= ascii.read('/mnt/window/Users/User/Downloads/data/Project/gw/GW170817/GW170817-all_candidates.txt')

tra = coord.Angle(cantbl['ra']*u.degree)
tra = tra.wrap_at(180*u.degree)
tdec = coord.Angle(cantbl['dec']*u.degree)
'''
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
#ax.set(title='S190425z')
#ax.legend()
plt.tight_layout()
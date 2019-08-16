#	KILONOVA MODEL TO K-CORRECTED MAGNITUDE BASED ON Masaomi Tanaka
#	REFERENCE	: 
#------------------------------------------------------------
#	2019.08.15	CREATED BY Gregory S.H.Paek
#============================================================
import numpy as np
import bisect
import matplotlib.pyplot as plt 
import os,sys,glob
from scipy.integrate import simps
from astropy.io import ascii
from astropy import units as u
from astropy import constants as const
from speclite import filters
from astropy.table import Table, Row, Column, vstack
#------------------------------------------------------------
#	FUNCTION
#------------------------------------------------------------
def extractinfo(intbl, z):
	zs = intbl['redshift']
	it = bisect.bisect(zs, z)
	if it == 0:
		it = 1
	newtbl = intbl[	(intbl['redshift']==intbl['redshift'][it-1])]
	return newtbl
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_table = '/home/sonic/Research/yourpy/gppy/table/Masaomi_Tanaka_Kilonova_Models'
path_save = '/data1/S190425z/1.result/table'
#------------------------------------------------------------
z = 0.009787
models = glob.glob(path_table+'/knova*obs.dat')
#------------------------------------------------------------
plt.close('all')
for mdl in models:
	mdltbl = extractinfo(ascii.read(mdl), z)
	# for filte in ['u','g','r','i','z','J','H','K']:
	for filte in ['r']:
		plt.plot(mdltbl['day'], mdltbl[filte], label='{} ({})'.format(filte, os.path.basename(mdl)))

# plt.xlim([0, 10])
# plt.ylim([25, 15])

plt.minorticks_on()
plt.tight_layout()
plt.legend()

'''
droutbl = ascii.read(path_save+'/lc_gw170817_Drout.dat')
rdroutbl = droutbl[droutbl['filter']=='r']
plt.scatter(rdroutbl['delmjd'], rdroutbl['mag'], marker='o', c='grey', label='Drout+17')
plt.xlim([-1, 11])
plt.ylim([25, 10])

plt.minorticks_on()
plt.tight_layout()
'''
#	IMSNG TARGET TIME DISTRIBUTION GENERATOR
#	CREATED	2019.11.06	Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import time
#============================================================
def rtsdayplot(rtsname):
	timestep = [17, 18, 19, 20, 21, 22, 23, 24, 0, 1, 2, 3, 4, 5, 6, 7]
	rtstbl = ascii.read(rtsname, header_start=9)
	trlist = []
	for i in range(len(rtstbl)):
		trlist.append( np.asscalar(rtstbl['transit(LT)'][i])[0:2] )
	subtbl = Table()
	subtbl['name'] = rtstbl['name']
	subtbl['tr'] = np.array(trlist)

	plt.close('all')
	plt.hist(subtbl['tr'])
	plt.xlabel('TIME (LT)', fontsize=20)
	plt.ylabel('Number of Targets', fontsize=20)
	plt.tight_layout()
	plt.savefig(rtsname[:-4]+'.png', overwrite=True)
#------------------------------------------------------------
def rtsmonthplot(rtslist):
	timestep = [17, 18, 19, 20, 21, 22, 23, 24, 0, 1, 2, 3, 4, 5, 6, 7]
	trlist = []
	for rtsname in rtslist:
		rtstbl = ascii.read(rtsname, header_start=9)
		for i in range(len(rtstbl)):
			trlist.append( np.asscalar(rtstbl['transit(LT)'][i])[0:2] )
	subtbl = Table()
	subtbl['tr'] = np.array(trlist)

	part = rtsname.split('_')
	part[4] = part[4][0:6]
	newname = '_'.join(part)

	plt.close('all')
	plt.hist(subtbl['tr'])
	plt.xlabel('TIME (LT)', fontsize=20)
	plt.ylabel('Number of Targets', fontsize=20)
	plt.tight_layout()
	plt.savefig(newname[:-4]+'.png', overwrite=True)
#============================================================
path_base = '/home/sonic/Research/CEOU/IMSNG/targetlist/rts_vis_2019/LOAO'
#------------------------------------------------------------
# rtsname = 'rts_vis_20190101_loao.txt'
rtslist = glob.glob(path_base+'/rts_vis_2019*.txt')	; rtslist.sort()
plt.ioff()
for rtsname in rtslist:
	print(rtsname)
	rtsdayplot(rtsname)
'''
for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
	plt.ioff()
	rtslist = glob.glob(path_base+'/rts_vis_2019'+month+'*.txt')
	rtsmonthplot(rtslist)
'''
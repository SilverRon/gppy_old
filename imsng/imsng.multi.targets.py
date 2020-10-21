#	IMSNG Multi Targets
#	2020.10.10	CREATED BY GREGORY S.H. PAEK
#============================================================
import numpy as np
import os, glob
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, vstack

#============================================================
#	FUNCTION
#============================================================

#============================================================
#	PATH
#------------------------------------------------------------
path_alltarget = '/data1/IMSNG/multi_target/alltarget.dat'
path_save = '/data1/IMSNG/multi_target'
#------------------------------------------------------------
fov = 0.55 * u.deg
#------------------------------------------------------------
alltbl = ascii.read(path_alltarget); alltbl = alltbl[alltbl['priority']<10.0]
targets = []
tblist = []
# i = 157
for i, obj in enumerate(alltbl['obj']):
	if obj not in targets:
		c_obj = SkyCoord(alltbl['ra'][i], alltbl['dec'][i], unit=(u.hourangle, u.deg))
		c = SkyCoord(alltbl['ra'], alltbl['dec'], unit=(u.hourangle, u.deg))

		sep = c_obj.separation(c)
		indx_inthefov = np.where(sep < fov)

		if len(indx_inthefov[0]) > 1:
			# print(obj)
			objlist = []
			for k in alltbl['obj'][indx_inthefov]:
				objlist.append(k)
			ralist = []
			delist = []

			c_inthefov = c[indx_inthefov]
			for indx in indx_inthefov[0]:
				ralist.append(c[indx].ra.value)
				delist.append(c[indx].dec.value)

			ra = np.mean(ralist)
			dec = np.mean(delist)
			c_center = SkyCoord(ra, dec, unit='deg')

			rahms = c_center.to_string('hmsdms').split(' ')[0].replace('h', ':').replace('m', ':').replace('s', '')
			dedms = c_center.to_string('hmsdms').split(' ')[1].replace('d', ':').replace('m', ':').replace('s', '')

			tmptbl = Table()
			tmptbl['obj'] = [obj]
			tmptbl['ra'] = rahms
			tmptbl['dec'] = dedms
			tmptbl['inthefov'] = ','.join(objlist)

			tblist.append(tmptbl)
			print(obj, ','.join(objlist), len(indx_inthefov[0]))
		else:
			objlist.append(obj)
		for l in objlist: targets.append(l)
	else:
		pass

comtbl = vstack(tblist)

comtbl.write('{}/imsng.multi_target.dat'.format(path_save), format='ascii.tab', overwrite=True)

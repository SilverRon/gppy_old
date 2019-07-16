#	WHAT SOURCES ARE IN FOV?
#	MADE BY Gregory S.H. Paek   (2019.07.19)
#============================================================#
import os, glob, sys
import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.table import Table, vstack
from astropy.io import ascii, fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import special
from astropy.wcs import WCS
#============================================================#
#    FUNCTION
#------------------------------------------------------------#
def infov(inim, tname, tra, tdec, namekey='name', draw=False):
	data = fits.getdata(inim)
	imx, imy = data.shape
	w = WCS(inim)
	targx, targy = w.wcs_world2pix(tra, tdec, 0)
	indx = np.where((targx>0)&
					(targx<imx)&
					(targy>0)&
					(targy<imy))
	tnames = tname[indx]
	tx, ty = targx[indx], targy[indx]
	if len(tnames)!=0:
		objs = ''
		for obj in tnames:
			objs += obj
			objs += ','
		objs = objs[:-1]
	else:
		objs = 'None'
	onetbl = Table(	[[inim], [objs]],
					names=['image', 'sources'])
	if draw!=False:
		plotshow(inim, tnames, tx, ty)
	return onetbl
#------------------------------------------------------------#
def plotshow(inim, tnames, tx, ty):
	'''
	PLOT IMAGE AND SHOW DESINATED OBJECTS
	'''
	import numpy as np
	import matplotlib
	import matplotlib.pyplot as plt
	from astropy.io import fits
	from matplotlib.colors import LogNorm
	from matplotlib.patches import Circle
	from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
	from astropy.visualization import ZScaleInterval, LinearStretch
	from astropy.wcs import WCS
	#------------------------------------------------------------#
	outname		= inim[:-5]+'.png'
	data, hdr	= fits.getdata(inim, header=True)
	wcs			= WCS(hdr)
	norm_zscale	= ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
	#------------------------------------------------------------#
	plt.close()
	fig			= plt.figure() 
	ax			= plt.subplot(projection=wcs)
	im			= ax.imshow(data, cmap='gray', origin='lower', norm=norm_zscale) 

	for xx, yy in zip(tx, ty):
		circ = Circle((xx, yy), 15, color='gold', fill=None, linewidth='1.0')
		ax.add_patch(circ)
	for i, txt in enumerate(tnames):
		xim		= tx[i]
		yim		= ty[i]
		ax.text(xim, yim+15, str(txt), color='gold', fontsize=15)
	plt.title(outname)
	plt.minorticks_on()
	fig.savefig(outname, dpi=500, facecolor='w', edgecolor='w',
				orientation='portrait', papertype=None, format=None,
				transparent=False, bbox_inches=None, pad_inches=0.1,
				frameon=None, metadata=None)
#------------------------------------------------------------#
path_table = '/mnt/window/Users/User/Downloads/data/Project/gw/S190425z/info/S190425z_Update-all_candidates.txt'
path_obs = '/home/sonic/Research/test/test/telescope.info'
targtbl = ascii.read(path_table)
obstbl = ascii.read(path_obs)
imlist = glob.glob('*.fits')
namekey = 'name'
#------------------------------------------------------------#
tblist = []
for inim in imlist:
	param_infov = dict(	inim=inim,
						tname=targtbl[namekey],
						tra=targtbl['ra'],
						tdec=targtbl['dec'],
						namekey=namekey,
						draw=True)
	tblist.append(infov(**param))
comtbl = vstack(tblist)
if 'inthefov.dat' in glob.glob('inthefov.dat'):
	os.system('mv inthefov.dat inthefov.dat.bkg')
comtbl.write('inthefov.dat', format='ascii', overwrite=True)
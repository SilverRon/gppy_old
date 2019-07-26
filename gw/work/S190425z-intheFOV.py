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
	#	TABLE TARGET RA DEC -> PIXEL X, Y
	targx, targy = [], []
	for i in range(len(tname)):
		x, y =  w.wcs_world2pix(tra[i], tdec[i], 0)
		targx.append(x)
		targy.append(y)
	targx, targy = np.array(targx), np.array(targy)
	np.warnings.filterwarnings('ignore')

	indx = np.where((targx>imx*0.05)&
					(targx<imx*0.95)&
					(targy>imy*0.05)&
					(targy<imy*0.95))
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
	onetbl = Table(	[[os.path.basename(inim)], [objs]],
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
	# outname		= inim[:-5]+'-sources.png'
	outname = './'+os.path.basename(inim)[:-5]+'-sources.png'
	data, hdr	= fits.getdata(inim, header=True)
	wcs			= WCS(hdr)
	norm_zscale	= ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
	#------------------------------------------------------------#
	plt.close()
	fig			= plt.figure() 
	ax			= plt.subplot(projection=wcs)
	im			= ax.imshow(data, cmap='gray', origin='lower', norm=norm_zscale) 

	for xx, yy in zip(tx, ty):
		circ = Circle((xx, yy), 15, color='gold', fill=None, linewidth='0.5')
		ax.add_patch(circ)
	for i, txt in enumerate(tnames):
		xim		= tx[i]
		yim		= ty[i]
		ax.text(xim, yim+15, str(txt), color='gold', fontsize=7)
	plt.title(outname)
	plt.minorticks_on()
	fig.savefig(outname, dpi=500, facecolor='w', edgecolor='w',
				orientation='portrait', papertype=None, format=None,
				transparent=False, bbox_inches=None, pad_inches=0.1,
				frameon=None, metadata=None)

#------------------------------------------------------------#
#------------------------------------------------------------#
imlist = glob.glob('*Calib*.fits')
# imlist = glob.glob(input('IMAGE TO INVESTIGATE\t: '))
# path_obs = '/home/sonic/Research/table/telescope.info'
# obstbl = ascii.read(path_obs)
namekey = 'name'
#------------------------------------------------------------#
path_table = '/data1/S190425z/info/Initial/S190425z_Initial-all_candi.txt'
targtbl = ascii.read(path_table)
tblist = []
i=0
for inim in imlist:
	i+=1
	print('PROCESS [{}/{}]\t: {}'.format(i, len(imlist), inim))
	param_infov = dict(	inim=inim,
						tname=targtbl[namekey],
						tra=targtbl['ra'],
						tdec=targtbl['dec'],
						namekey=namekey,
						draw=True)
	tblist.append(infov(**param_infov))
comtbl = vstack(tblist)
if 'inthefov_initial.dat' in glob.glob('inthefov_initial.dat'):
	os.system('mv inthefov_initial.dat inthefov_initial.dat.bkg')
comtbl.write('inthefov_initial.dat', format='ascii', overwrite=True)
#------------------------------------------------------------#
path_table = '/data1/S190425z/info/Update/S190425z_Update-all_candi.txt'
targtbl = ascii.read(path_table)
tblist = []
failist = []
i=0
for inim in imlist:
	i+=1
	print('PROCESS [{}/{}]\t: {}'.format(i, len(imlist), inim))
	param_infov = dict(	inim=inim,
						tname=targtbl[namekey],
						tra=targtbl['ra'],
						tdec=targtbl['dec'],
						namekey=namekey,
						draw=True)
	if 'scCalib' not in inim:
		tblist.append(infov(**param_infov))

comtbl = vstack(tblist)
if 'inthefov_update.dat' in glob.glob('inthefov_update.dat'):
	os.system('mv inthefov_update.dat inthefov_update.dat.bkg')
comtbl.write('inthefov_update.dat', format='ascii', overwrite=True)


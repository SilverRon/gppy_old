#	WHAT SOURCES ARE IN FOV?
#	MADE BY Gregory S.H. Paek   	(2019.07.19)
#	MODIFIED BY Gregory S.H. Paek	(2019.07.29)
#	MODIFIED BY Gregory S.H. Paek	(2019.08.05)
#============================================================
import os, glob, sys
import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.table import Table, vstack
from astropy.io import ascii, fits
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from scipy import special
from astropy.wcs import WCS
from astropy import units as u
#============================================================
#    FUNCTION
#------------------------------------------------------------
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
#------------------------------------------------------------
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
	#------------------------------------------------------------
	# outname		= inim[:-5]+'-sources.png'
	outname = './'+os.path.basename(inim)[:-5]+'-sources.png'
	data, hdr	= fits.getdata(inim, header=True)
	wcs			= WCS(hdr)
	norm_zscale	= ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
	#------------------------------------------------------------
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
#------------------------------------------------------------
def centfind(inim):
	hdr = fits.getheader(inim)
	xcent, ycent = hdr['NAXIS1']/2, hdr['NAXIS2']/2
	w = WCS(inim)
	racent0, decent0 = w.wcs_pix2world(xcent, ycent, 0)
	racent, decent = np.asscalar(racent0), np.asscalar(decent0)
	return racent, decent
#------------------------------------------------------------
fov = 1.0*u.deg	# [deg]
#------------------------------------------------------------
#	INITIAL
#------------------------------------------------------------
path_table = '/data1/S190425z/info/Initial/S190425z_Initial-all_candi.txt'
imlist = glob.glob('*.fits')
targtbl = ascii.read(path_table)
# obstbl = ascii.read(path_obs)
#------------------------------------------------------------
tblist = []
tname = np.copy(targtbl['name'])
tra = np.copy(targtbl['ra'])
tdec = np.copy(targtbl['dec'])
tcoord = SkyCoord(tra, tdec, unit='deg')
#------------------------------------------------------------
i=0
for inim in imlist:
	i+=1
	print('PROCESS [{}/{}]\t: {}'.format(i, len(imlist), inim))
	racent, decent = centfind(inim)
	imcoord = SkyCoord(racent, decent, unit='deg')
	sep = imcoord.separation(tcoord)
	indx = np.where(sep < fov/2)
	if len(indx[0]) != 0:
		sources = ''
		for n in range(len(indx[0])):
			sources = sources+tname[indx[0][n]]+','
		sources = sources[:-1]
	else:
		sources = 'None'
	onetbl = Table([[inim], [sources]], names=('image', 'sources'))
	tblist.append(onetbl)
comtbl = vstack(tblist)
#------------------------------------------------------------
if 'inthefov_initial.dat' in glob.glob('inthefov_initial.dat'):
	os.system('mv inthefov_initial.dat inthefov_initial.dat.bkg')
comtbl.write('inthefov_initial.dat', format='ascii', overwrite=True)











#------------------------------------------------------------
#	UPDATE
#------------------------------------------------------------
path_table = '/data1/S190425z/info/Update/S190425z_Update-all_candi.txt'
imlist = glob.glob('*.fits')
targtbl = ascii.read(path_table)
#------------------------------------------------------------
tblist = []
tname = np.copy(targtbl['name'])
tra = np.copy(targtbl['ra'])
tdec = np.copy(targtbl['dec'])
tcoord = SkyCoord(tra, tdec, unit='deg')
#------------------------------------------------------------
i=0
for inim in imlist:
	i+=1
	print('PROCESS [{}/{}]\t: {}'.format(i, len(imlist), inim))
	racent, decent = centfind(inim)
	imcoord = SkyCoord(racent, decent, unit='deg')
	sep = imcoord.separation(tcoord)
	indx = np.where(sep < fov/2)
	if len(indx[0]) != 0:
		sources = ''
		for n in range(len(indx[0])):
			sources = sources+tname[indx[0][n]]+','
		sources = sources[:-1]
	else:
		sources = 'None'
	onetbl = Table([[inim], [sources]], names=('image', 'sources'))
	tblist.append(onetbl)
comtbl = vstack(tblist)
if 'inthefov_update.dat' in glob.glob('inthefov_update.dat'):
	os.system('mv inthefov_update.dat inthefov_update.dat.bkg')
comtbl.write('inthefov_update.dat', format='ascii', overwrite=True)


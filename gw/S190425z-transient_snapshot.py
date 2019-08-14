#	FITS SNAPSHOT
#	2019.08.14	CREATED BY	Gregory S.H. Paek
#============================================================
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
from astropy.visualization import ZScaleInterval, LinearStretch
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
#------------------------------------------------------------
path_image = '/data1/S190425z/SQUEAN/transient'
path_save = '/data1/S190425z/1.result/figure'
# inim = 'trCalib-SQUEAN-ZTF19aarykkb-20190426-084629-i-60_gregister-com.fits'
imlist = glob.glob(path_image+'/*tr*ZTF19aarykkb*ter*.fits')
imlist.remove(path_image+'/hctrCalib-SQUEAN-ZTF19aarykkb-20190426-084629-i-60_gregister-com.fits')
imlist.remove(path_image+'/atrRef-PS1-ZTF19aarykkb-00000000-000000-i_gregister.fits')
#------------------------------------------------------------
for inim in imlist:
	if '/trCalib-' in inim:
		title = 'Science'
		outname = 'ZTF19aarykkb-sci.png'
	elif '/hdtrCalib-' in inim:
		title = 'Difference'
		outname = 'ZTF19aarykkb-dif.png'
	elif '/trRef-' in inim:
		title = 'Reference'
		outname = 'ZTF19aarykkb-ref.png'
	data, hdr	= fits.getdata(inim, header=True)
	# data = -1*data
	wcs			= WCS(hdr)
	norm_zscale	= ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
	#------------------------------------------------------------
	param_imshow = dict(X=data,
						cmap='bone', origin='lower',
						norm=norm_zscale)
	plt.close('all')
	fig			= plt.figure(figsize=[6, 6])
	ax			= plt.subplot(projection=wcs)
	im			= ax.imshow(**param_imshow)
	plt.grid(color='white', ls='solid')
	plt.title(title, fontsize=20)
	plt.xlabel('R.A.', fontsize=20)
	plt.ylabel('Dec.', fontsize=20)
	plt.minorticks_on()
	# plt.tight_layout()
	txt = 'ZTF19aarykkb'
	xx, yy =  279.50034, 309.25818
	circ = Circle((xx, yy), 15, color='gold', fill=None, linewidth='1.0')
	ax.add_patch(circ)
	fig.savefig(path_save+'/'+outname, overwrite=True)
'''
fig.savefig(outname, dpi=500, facecolor='w', edgecolor='w',
	orientation='portrait', papertype=None, format=None,
	transparent=False, bbox_inches=None, pad_inches=0.1,
	frameon=None, metadata=None)
'''

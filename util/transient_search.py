#	TRANSIENT SEARCH ALGORITHMS
#	2020.02.28	CREATED BY Gregory S.H. Paek
#	2020.03.02	UPDATED BY Gregory S.H. Paek, SKYBOT FUNCTION
#============================================================
import os, glob, subprocess
from astropy.io import ascii, fits
from astropy.table import Table, vstack, Row
from astropy import units as u
from astroquery.imcce import Skybot
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.time import Time
import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
from astropy.visualization import ZScaleInterval, LinearStretch
#============================================================
def sexcom(inim, conf_sex, conf_param, conf_conv, conf_nnw, peeing):
	outcat = inim[:-5]+'.cat'
	sexcom = 'sex {} -c {} -CATALOG_NAME {} -PARAMETERS_NAME {} -FILTER_NAME {} -STARNNW_NAME {} -PHOT_APERTURES {}'.format(inim, conf_sex, outcat, conf_param, conf_conv, conf_nnw, 2*peeing)
	# print(sexcom)
	# os.system(sexcom)
	return sexcom
#------------------------------------------------------------
def inverse(inim, outim):
	data, hdr = fits.getdata(inim, header=True)
	# w = WCS(inim)
	invdata = data*(-1)
	fits.writeto(outim, invdata, header=hdr, overwrite=True)
#------------------------------------------------------------
def gethead(inim, key):
	hdr = fits.getheader(inim)
	return hdr[key]
#------------------------------------------------------------
def ds9reg(starname, ra, dec, outname='ds9.reg', size=5, color='yellow'):
	radius	= """ {0}" """.format(size)
	color	= "yellow"
	os.system('rm '+outname)
	f		= open(outname,'w')
	head1	= "# Region file format: DS9 version 4.1\n"; f.write(head1)
	head2	= """global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"""; f.write(head2)
	head3	= "fk5\n"; f.write(head3)
	for n in range(len(starname)):
		body="circle("+str(ra[n])+","+str(dec[n])+","+radius+") # color="+color+" text={"+str(starname[n])+"}\n"	
		f.write(body)
	f.close()								
#------------------------------------------------------------
def snapshot(inim, refim, subim, intbl, i):
	'''
	PLOT IMAGE AND SHOW DESINATED OBJECTS
	'''

	#------------------------------------------------------------
	indata, inhdr = fits.getdata(inim, header=True)
	refdata, refhdr = fits.getdata(refim, header=True)
	subdata, subhdr = fits.getdata(subim, header=True)
	#------------------------------------------------------------
	# i = 0
	number = intbl[i]['NUMBER']
	x, y = int(intbl[i]['X_IMAGE']), int(intbl[i]['Y_IMAGE'])

	if np.sum(refdata[y-1:y+1, x-1:x+1]) != 0.0:

		xim, yim = indata.shape
		ra, dec = round(intbl[i]['ALPHA_J2000'], 3), round(intbl[i]['DELTA_J2000'], 3) 
		zp, zper = inhdr['OPTZP'], inhdr['OPTZPERR']
		mag, mager = intbl[i]['MAG_AUTO']+zp, np.sqrt( (intbl[i]['MAGERR_AUTO'])**2 + (zper)**2 )
		x0, x1, y0, y1 = x-int(xim/r), x+int(xim/r), y-int(yim/r), y+int(yim/r)
		#------------------------------------------------------------
		if intbl['SkyBot_name'][i] == 'None':
			tale = ''
			txt = ''
		else:
			tale = '.skybot'
			txt = intbl['SkyBot_name'][i]
		outname		= '{}.{}_{}_{}_{}_{}{}.png'.format(inim[:-5], number, ra, dec, round(mag, 3), round(mager, 3), tale)
		#------------------------------------------------------------
		plt.close('all')
		fig = plt.figure(figsize=(15, 5))
		#------------------------------------------------------------
		ax0 = plt.subplot(1, 3, 1)
		# inwcs = WCS(inhdr)
		norm_zscale	= ImageNormalize(-1*indata[y0:y1, x0:x1], interval=ZScaleInterval(), stretch=LinearStretch())
		im0 = ax0.imshow(-1*indata[y0:y1, x0:x1], cmap='gray', origin='lower', norm=norm_zscale)
		circ0 = Circle((xim/r, yim/r), r/2, color='gold', fill=None, linewidth='2.0')
		ax0.add_patch(circ0)
		ax0.set_title('Number : {}'.format(intbl[i]['NUMBER']))
		ax0.text(xim/r, yim/r, str(txt), color='gold', fontsize=10)
		#------------------------------------------------------------
		ax1 = plt.subplot(1, 3, 2)
		# refwcs = WCS(refhdr)
		norm_zscale	= ImageNormalize(-1*refdata[y0:y1, x0:x1], interval=ZScaleInterval(), stretch=LinearStretch())
		im1 = ax1.imshow(-1*refdata[y0:y1, x0:x1], cmap='gray', origin='lower', norm=norm_zscale)
		circ1 = Circle((xim/r, yim/r), r/2, color='gold', fill=None, linewidth='2.0')
		ax1.add_patch(circ1)
		ax1.set_title('RA, Dec = {}, {}'.format(ra, dec))
		#------------------------------------------------------------
		ax2 = plt.subplot(1, 3, 3)
		# subwcs = WCS(subhdr)
		norm_zscale	= ImageNormalize(-1*subdata[y0:y1, x0:x1], interval=ZScaleInterval(), stretch=LinearStretch())
		im2 = ax2.imshow(-1*subdata[y0:y1, x0:x1], cmap='gray', origin='lower', norm=norm_zscale) 
		circ2 = Circle((xim/r, yim/r), r/2, color='gold', fill=None, linewidth='2.0')
		ax2.add_patch(circ2)
		ax2.set_title('MAG = {}+/-{}'.format(round(mag, 3), round(mager, 3)))
		#------------------------------------------------------------
		fig.suptitle(inim, fontsize=20)
		# plt.tight_layout()
		# ax.text(xim/r, yim/r, str(txt), color='gold', fontsize=5)
		fig.savefig(outname, dpi=300, facecolor='w', edgecolor='w',
					orientation='portrait', papertype=None, format=None,
					transparent=False, bbox_inches=None, pad_inches=0.1,
					metadata=None)
	else:
		pass
#------------------------------------------------------------
def skybot(intbl, dateobs, sep=10):

	strlist = []
	for i in range(len(intbl)): strlist.append('                    ')

	intbl['SkyBot_name'] = strlist
	intbl['SkyBot_Vmag'] = np.zeros(len(intbl))
	intbl['SkyBot_ra'] = np.zeros(len(intbl))
	intbl['SkyBot_dec'] = np.zeros(len(intbl))
	intbl['SkyBot_epoch'] = np.zeros(len(intbl))
	intbl['SkyBot_sep'] = np.zeros(len(intbl))

	# sep = 10		#	["]
	for i in range(len(intbl)):
		print('[{}/{}]'.format(i+1, len(intbl)))
		field = SkyCoord(intbl['ALPHA_J2000'][i]*u.deg, intbl['DELTA_J2000'][i]*u.deg)
		epoch = Time(dateobs, format='isot')

		try:
			sbtbl = Skybot.cone_search(field, sep*u.arcsec, epoch)
			c_sb = SkyCoord(sbtbl['RA'], sbtbl['DEC'])
			indx, d2d, d3d = field.match_to_catalog_sky(c_sb)

			sbonetbl = sbtbl[d2d == np.min(d2d)]
			print(sbonetbl)
			if len(sbonetbl) == 1:
				intbl['SkyBot_name'][i] = sbonetbl['Name'][0]
				intbl['SkyBot_Vmag'][i] = sbonetbl['V'][0].to_value()
				intbl['SkyBot_ra'][i] = sbonetbl['RA'][0].to_value()
				intbl['SkyBot_dec'][i] = sbonetbl['DEC'][0].to_value()
				intbl['SkyBot_epoch'][i] = sbonetbl['epoch'][0].to_value()
				intbl['SkyBot_sep'][i] = d2d.item().arcsec
			# elif len(sbone) == 0:
		except:
			intbl['SkyBot_name'][i] = 'None'
	return intbl
#------------------------------------------------------------
def routine(inim):
	refim = 'wr{}'.format(inim)
	subim = 'hd{}'.format(inim)
	invim = 'inv{}'.format(subim)
	outcat = inim[:-5]+'.transient.cat'
	outreg = subim[:-5]+'.reg'
	hdr = fits.getheader(subim)
	skyval, skysig = hdr['SKYVAL'], hdr['SKYSIG']
	optzp, optzper = hdr['OPTZP'], hdr['OPTZPERR']
	ul_3sig = hdr['OPTUL']
	peeing = hdr['PEEING']
	#------------------------------------------------------------
	inverse(subim, invim)
	os.system(sexcom(refim, conf_sex, conf_param, conf_conv, conf_nnw, peeing))
	os.system(sexcom(invim, conf_sex, conf_param, conf_conv, conf_nnw, peeing))
	#	SE FOR SUBT. IMAGE & BKG. ESTIMATION
	os.system(sexcom(subim, conf_sex, conf_param, conf_conv, conf_nnw, peeing))
	#	IMAGE X, Y
	xim, yim = fits.getdata(subim).shape
	xcent, ycent = xim/2, yim/2
	#------------------------------------------------------------
	reftbl = ascii.read(refim[:-5]+'.cat')
	subtbl = ascii.read(subim[:-5]+'.cat')
	invtbl = ascii.read(invim[:-5]+'.cat')
	#------------------------------------------------------------
	refc = SkyCoord(reftbl['ALPHA_J2000'], reftbl['DELTA_J2000'], frame='icrs', unit='deg')
	invc = SkyCoord(invtbl['ALPHA_J2000'], invtbl['DELTA_J2000'], frame='icrs', unit='deg')

	#	COSMIC RAY & IMAGE CENTER
	indx_cr = np.where(	
						# (subtbl['A_IMAGE'] < 1.2) &
						(subtbl['ELONGATION'] < np.mean(subtbl['ELONGATION']) + 1*np.std(subtbl['ELONGATION'])) &
						(subtbl['ELLIPTICITY'] < np.mean(subtbl['ELLIPTICITY']) + 1*np.std(subtbl['ELLIPTICITY'])) &
						(np.sqrt((subtbl['X_IMAGE']-xcent)**2 + (subtbl['Y_IMAGE']-ycent)**2) < frac*(xcent+ycent)/2.0 ) &
						(subtbl['MAG_AUTO']+hdr['OPTZP'] < hdr['OPTUL']+1) &
						(subtbl['FLAGS'] <= 1)
					)
	subtbl1 = subtbl[indx_cr]
	subc1 = SkyCoord(subtbl1['ALPHA_J2000'], subtbl1['DELTA_J2000'], frame='icrs', unit='deg')
	#	INV. MATCHING
	indx_inv, sep_inv, _ = subc1.match_to_catalog_sky(invc)
	subtbl2 = subtbl1[sep_inv.arcsec>hdr['SEEING']]
	subc2 = SkyCoord(subtbl2['ALPHA_J2000'], subtbl2['DELTA_J2000'], frame='icrs', unit='deg')
	#	REF. MATCHING
	indx_ref, sep_ref, _ = subc2.match_to_catalog_sky(refc)
	subtbl3 = subtbl2[sep_ref.arcsec>hdr['SEEING']]
	subtbl3['inim'] = inim
	subtbl3['refim'] = refim
	subtbl3['subim'] = subim

	subtbl4 = skybot(subtbl3, hdr['DATE-OBS'], sep=10)

	subtbl4.write(outcat, format='ascii', overwrite=True)

	#	MAKING PLOT
	print('MAKING PLOT')
	for i in range(len(subtbl3)):
		snapshot(inim, refim, subim, subtbl4, i)

	ds9reg(subtbl4['NUMBER'], subtbl4['ALPHA_J2000'], subtbl4['DELTA_J2000'], outname=outreg)

	ds9com = 'ds9 {} {} {} -region {} -frame lock wcs -tile column &'.format(inim, refim, subim, outreg)
	return ds9com
#============================================================
path_config = '/home/sonic/Research/yourpy/gppy/config'
conv, nnw, param, sex = 'transient.conv', 'transient.nnw', 'transient.param', 'transient.sex'
conf_sex = '{}/{}'.format(path_config, sex)
conf_param = '{}/{}'.format(path_config, param)
conf_nnw = '{}/{}'.format(path_config, nnw)
conf_conv = '{}/{}'.format(path_config, conv)
#------------------------------------------------------------
#	SEP. ["]
sep = 1
#	FRACTION
frac = 0.9
#	IMAGE SNAPSHOT (IMAGE SIZE/r)
r = 25
#------------------------------------------------------------
imlist = glob.glob('caCalib-*.fits'); imlist.sort()
failist = []
for inim in imlist:
	try:
		# os.system(routine(inim))
		ds9com = routine(inim)
	except:
		failist.append(inim)
		pass

'''
indata, refdata, subdata = fits.getdata(inim), fits.getdata(refim), fits.getdata(subim)
for i in range(len(subtbl3)):

	x, y = int(subtbl3[i]['X_IMAGE']), int(subtbl3[i]['Y_IMAGE'])
	x0, x1, y0, y1 = x-int(xim/r), x+int(xim/r), y-int(yim/r), y+int(yim/r)

	plt.close()
	plt.subplots(1, 3, figsize=(15, 5))
	plt.subplot(1, 3, 1)
	plt.imshow(indata[x0:x1, y0:y1], cmap='gray', vmin=skyval, vmax=np.max(indata[x0:x1, y0:y1])*1e0, norm=colors.LogNorm())
	plt.subplot(1, 3, 2)
	plt.imshow(-1*refdata[x0:x1, y0:y1], cmap='gray')
	plt.subplot(1, 3, 3)
	plt.imshow(subdata[x0:x1, y0:y1], cmap='gray', vmin=0)
	plt.tight_layout()
	plt.savefig('test_{}.png'.format(i), overwrite=True)
'''
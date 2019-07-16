#   CONSIDER MULTIPLE GALAXIES IN AN IMAGE
#   2019.07.09 MADE		BY Gregory S.H. Paek
#	2019.XX.XX UPDATED	BY Gregory S.H. Paek
#============================================================
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
import time
import os
from astropy.table import Table, Column, MaskedColumn
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import matplotlib.pyplot as plt
from imsng import gw
from imsng import tool
from astropy.table import Table, vstack
from tqdm import tqdm, trange
#============================================================
#	FUNCTION
#------------------------------------------------------------
def func_linear(x, x0, y0, x1, y1):
	a		= (y1-y0)/(x1-x0)
	b		= y0-a*x0
	return a*x+b
#------------------------------------------------------------
rtstbl		= ascii.read('rts_sample.txt')
# rtstbl = ascii.read('/home/sonic/Downloads/S190425z_Update-20190426-rts_vis-LOAO.txt')
targlist	= np.copy(rtstbl['name'])	;targlist	= targlist.tolist()
fov			= 27.3*u.arcmin
foveff		= fov*0.8
#------------------------------------------------------------
tblist		= []
# pbar		= tqdm(altarglist)
pbar		= tqdm(targlist)
for targ in pbar:
	pbar.set_description('Processing %s' % targ)
# for targ in targlist:
	# print(targ)
	i		= np.where(rtstbl['name']==targ)
	coord_cent	= SkyCoord(	rtstbl['ra'][i], rtstbl['dec'][i],
							unit=(u.hourangle, u.deg))
	coord_other	= SkyCoord( rtstbl['ra'], rtstbl['dec'],
							unit=(u.hourangle, u.deg))
	demin, demax= coord_cent.dec-foveff/2., coord_cent.dec+foveff/2.

	raupmin, raupmax	= coord_cent.ra-foveff/2./np.cos(demax), coord_cent.ra+foveff/2./np.cos(demax)
	ralomin, ralomax	= coord_cent.ra-foveff/2./np.cos(demin), coord_cent.ra+foveff/2./np.cos(demin)
	#------------------------------------------------------------
	if (raupmax-raupmin)>(ralomax-ralomin):
		indx_within	= np.where(	(coord_other.dec<demax)&
								(coord_other.dec>demin)&
								(coord_other.dec>func_linear(coord_other.ra, ralomin, demin, raupmin, demax))&
								(coord_other.dec>func_linear(coord_other.ra, ralomax, demin, raupmax, demax)))
	else:
		indx_within	= np.where(	(coord_other.dec<demax)&
								(coord_other.dec>demin)&
								(coord_other.dec<func_linear(coord_other.ra, ralomin, demin, raupmin, demax))&
								(coord_other.dec<func_linear(coord_other.ra, ralomax, demin, raupmax, demax)))
	#------------------------------------------------------------
	onetbl		= Table()
	for key in rtstbl.keys()[:-2]: onetbl[key]	= rtstbl[indx_within[0]][key]
	if len(indx_within[0]) == 1:
		multi_cand	= 'None'
	else:
		multi_cand	= ''
		for gal in rtstbl['name'][indx_within]:
			try:
				targlist.remove(gal)
				gal		+= ','
				multi_cand	+= gal
			except:
				pass
		multi_cand	= multi_cand[:-1]
		onetbl['name'][0]= 'F'+onetbl['name'][0][1:]	# 'F' MEANS 'FIELD'
	onetbl['note']		= [multi_cand]
	tblist.append(onetbl[0])
#------------------------------------------------------------
fintbl		= vstack(tblist)



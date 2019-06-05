#   GET GW INFORMATION
#   BASED ON https://emfollow.docs.ligo.org/userguide/tutorial/index.html
#   2019.03.03 MADE		BY Gregory S.H. Paek
#	2019.03.28 UPDATED	BY Gregory S.H. Paek
#	2019.04.18 UPDATED	BY Gregory S.H. Paek
#	2019.05.02 UPDATED	BY Gregory S.H. Paek (ver. 2.0)
#	2019.05.10 UPDATED	BY Gregory S.H. Paek
#	2019.05.29 UPDATED	BY Gregory S.H. Paek
#============================================================#
#	MODULE
#------------------------------------------------------------#
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
import time
import os, glob
from astropy.table import Table, Column, MaskedColumn
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import matplotlib.pyplot as plt
from imsng import gw
from imsng import tool
from astropy.table import Table, vstack
from astropy.time import Time 
#------------------------------------------------------------#
#	ROUTINE FUNCTION
#------------------------------------------------------------#
def bayestar_routine(healpixfits, save_path, confidence=0.9):
	eventname	= healpixfits.split('_')[0]
	eventtype	= healpixfits.split('_')[1].split('-')[0]
	param_heal	= dict(	healpixfits=healpixfits,
						confidence=confidence,
						eventname=eventname,
						save_path=save_path,
						hdr=True,
						view=False)
	healtbl, hdr= gw.heal2table(**param_heal)
	#healtbl.write(save_path+'/'+eventname+'_healpix.dat', format='ascii', overwrite=True)
	#	COVERED SKY REGION
	npix        = hdr['NAXIS2']         # i.e. 50331648
	alldeg		= (4*np.pi)*((180./np.pi)**2)         	#	[deg^2]
	skydeg_cut	= (alldeg/npix)*len(healtbl)
	#skydeg_cut_N= (alldeg/npix)*len(healtbl[healtbl['dec']>0])
	#skydeg_cut_S= (alldeg/npix)*len(healtbl[healtbl['dec']<0])
	'''
	result		= {	'event'		: eventname,
					'instrume'	: hdr['INSTRUME'],
					'date-obs'	: hdr['DATE-OBS'],
					'mjd-obs'	: hdr['MJD-OBS'],
					'dist'		: hdr['DISTMEAN'],
					'diststd'	: hdr['DISTSTD'],
					'region'	: skydeg_cut}
	'''
	outbl			= Table()
	outbl['event']	= [eventname]
	outbl['type']	= [eventtype]
	keylist			= ['instrume', 'date-obs', 'mjd-obs', 'distmean', 'diststd']
	keylist0		= ['date-obs', 'mjd-obs', 'distmean', 'diststd']
	try:
		for key in keylist:
			outbl[key]	= [hdr[key.upper()]]
	except:
		for key in keylist0:
			outbl[key]	= [hdr[key.upper()]]
	outbl['region']	= [skydeg_cut]
	#outbl		= Table([eventname, hdr['INSTRUME'], hdr['DATE-OBS'], hdr['MJD-OBS'], hdr['DISTMEAN'], hdr['DISTSTD'], skydeg_cut], names=['event', 'instrume', 'date-obs', 'mjd-obs', 'dist', 'diststd', 'region'])
	return outbl
#------------------------------------------------------------#
#	INPUT
#------------------------------------------------------------#
os.system('ls *.fits *fits.gz')
healpixlist	= glob.glob(input('BAYESTARS TO PROCESS\t: '))
#healpixfits	= input('BAYSTAR PATH\t: ')
#confidence	= float(input('CONFIDENCE (0.5,0.9)\t: '))
#if confidence	== '': confidence = 0.5
confidence	= 0.9
#------------------------------------------------------------#
save_healfix	= './'
#save_path		= save_healfix+eventname
save_path		= save_healfix
#os.system('mkdir '+save_path)
#------------------------------------------------------------#
tblist		= []
failist		= []
for healpixfits in healpixlist:
	print(healpixfits)
	try:
		outbl	= bayestar_routine(healpixfits, save_path)
		tblist.append(outbl)
	except:
		failist.append(healpixfits)
comtbl	= vstack(tblist)

#	O3 run START
t0 = Time('2019-04-01T00:00:00', format='isot', scale='utc')
comtbl['delmjd']	= comtbl['mjd-obs']-t0.mjd

if 'O3run_summary_from_bayestar.dat' in glob.glob(save_path+'*.dat'):
	os.system('mv {}O3run_summary_from_bayestar.dat {}O3run_summary_from_bayestar.dat.bkg').format(save_path, save_path)
comtbl.write('O3run_summary_from_bayestar.dat', format='ascii')
'''
#------------------------------------------------------------#
#	WRITE SUMMARY
#------------------------------------------------------------#
f	= open(save_path+'/'+eventname+'-summary.txt', 'w')

f.write('EXPECTED KILONOVA APP.MAG IN i-BAND [AB] : '+str(round(mag, 2))+'+/-'+str(round(magerr, 2))+'\n')
f.write(str(confidence*100)+'% confidence area\t: '+str(skydeg_cut)+'\t[deg^2]\n')
f.write(str(confidence*100)+'% confidence area for Northern sphere\t: '+str(skydeg_cut_N)+'\t[deg^2]\n')
f.write(str(confidence*100)+'% confidence area for Southern sphere\t: '+str(skydeg_cut_S)+'\t[deg^2]\n')
f.write('CENTER_RA\t: '+str(healtbl[ healtbl['P_2D'] == np.max(healtbl['P_2D']) ]['ra'][0])+'\n')
f.write('CENTER_Dec\t: '+str(healtbl[ healtbl['P_2D'] == np.max(healtbl['P_2D']) ]['dec'][0])+'\n')
f.write('MAX_RA\t: '+str(np.max(healtbl['ra']))+'\n')
f.write('MAX_Dec\t: '+str(np.max(healtbl['dec']))+'\n')
f.write('MIN_RA\t: '+str(np.min(healtbl['ra']))+'\n')
f.write('MIN_Dec\t: '+str(np.min(healtbl['dec']))+'\n')
f.close()
os.system('chmod 777 '+save_path+'/*')
'''
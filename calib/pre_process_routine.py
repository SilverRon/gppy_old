#============================================================
#   CALIBTRATION ROUTINE FOR IMSNG TELESCOPES
#	LOAO, DOAO, SOAO, ...
#   18.12.04	UPDATED BY Gregory S.H. Paek
#	19.08.08	UPDATED BY Gregory S.H. Paek	
#============================================================
import os, sys, glob
import numpy as np
import astropy.io.ascii as ascii
from astropy.nddata import CCDData
from imsng import calib
from astropy.io import fits
import warnings
# from __future__ import print_function
warnings.filterwarnings(action='ignore')
#============================================================
#	SETTING
#============================================================
#	PATH
#------------------------------------------------------------
# obs = 'DOAO'
try:
	obs = (sys.argv[1]).upper()
	print(obs)
except:
	obs = (raw_input('OBSERVATORY? (LOAO/DOAO/SOAO)\t: ')).upper()
	print(obs)
path_base = '/data3/paek/factory'
path_gal = '/data3/IMSNG/IMSNGgalaxies'
#------------------------------------------------------------
path_mframe = path_base+'/master_frames'
path_calib = path_base+'/calib'
path_ref = path_base+'/ref_frames/'+obs
path_log = '/home/paek/log/'+obs.lower()+'.log'
if obs == 'DOAO':
	path_raw = path_base+'/'+obs
else:
	path_raw = '/data3/grb/'+obs
path_factory = '/data3/paek/factory/'+obs.lower()
path_save = '/data3/IMSNG/'+obs
#------------------------------------------------------------
path_obsinfo = '/home/paek/table/obs.dat'
path_changehdr = '/home/paek/table/changehdr.dat'
#------------------------------------------------------------
logtbl = ascii.read(path_log)
datalist = np.copy(logtbl['date'])
obstbl = ascii.read(path_obsinfo)
hdrtbl = ascii.read(path_changehdr)
#============================================================
#	MAIN BODY
#============================================================
newlist = calib.folder_check(datalist, path_raw, path_factory)
# path_data=newlist[0]
for path_data in newlist:
	#------------------------------------------------------------
	#	CCD TYPE
	#------------------------------------------------------------
	if fits.getheader(glob.glob(path_data+'/*.fits')[0])['INSTRUME'] == 'Apogee USB/Net':
		obs = 'DOAO_APOGEE'
	if fits.getheader(glob.glob(path_data+'/*.fits')[0])['INSTRUME'] == '':
		obs = 'DOAO_FLI'
	if fits.getheader(glob.glob(path_data+'/*.fits')[0])['INSTRUME'] == 'FLI':
		obs = 'DOAO_FLI'
	obsinfo = calib.getobsinfo(obs, obstbl)
	print(obs)
	objfilterlist, objexptimelist, flatfilterlist, darkexptimelist, obstime = calib.correcthdr_routine(path_data, hdrtbl)
	images = calib.call_images(path_data)
	#------------------------------------------------------------
	#	BIAS
	#------------------------------------------------------------
	try:
		mzero = calib.master_zero(images, fig=False)
	except:
		#	IF THERE IS NO FLAT FRAMES, BORROW FROM CLOSEST OTHER DATE
		print('\nNO BIAS FRAMES\n')
		pastzero = np.array(glob.glob('{}/{}/zero/*zero.fits'.format(path_mframe, obs)))
		#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
		deltime = []
		for date in pastzero:
			zeromjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
			deltime.append(np.abs(obstime.mjd-zeromjd))
		indx_closet = np.where(deltime == np.min(deltime))
		tmpzero = path_data+'/'+os.path.basename(np.asscalar(pastzero[indx_closet]))
		cpcom = 'cp {} {}'.format(np.asscalar(pastzero[indx_closet]), tmpzero)
		print(cpcom)	; os.system(cpcom)
		mzero = CCDData.read(tmpzero, hdu=0)#, unit='adu')
	#------------------------------------------------------------
	#	DARK (ITERATION FOR EACH EXPOSURE TIMES)
	#------------------------------------------------------------
	if len(darkexptimelist) >= 1:
		dark_process = True
		for i, exptime in enumerate(darkexptimelist):
			print('PRE PROCESS FOR DARK ({} sec)\t[{}/{}]'.format(exptime, i+1, len(darkexptimelist)))
			mdark = calib.master_dark(images, mzero=mzero, exptime=exptime, fig=False)
	else:
		print('THERE IS NO DARK FRAME. DARK PROCESS WILL BE IGNORED.')
		dark_process = False
	#------------------------------------------------------------
	#	FLAT (ITERATION FOR EACH FILTERS)
	#------------------------------------------------------------
	if dark_process == False:
		for i, filte in enumerate(objfilterlist):
			print('PRE PROCESS FOR {} FILTER OBJECT\t[{}/{}]'.format(filte, i+1, len(objfilterlist)))
			if filte in flatfilterlist:
				mflat = calib.master_flat(images, mzero, filte, mdark=False, fig=True)
			else:
				#	IF THERE IS NO FLAT FRAMES, BORROW FROM CLOSEST OTHER DATE
				print('\nNO {} FLAT FRAMES\n'.format(filte))
				pastflat = np.array(glob.glob('{}/{}/flat/*n{}.fits'.format(path_mframe, obs, filte)))
				#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
				deltime = []
				for date in pastflat:
					flatmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
					deltime.append(np.abs(obstime.mjd-flatmjd))
				# mframe_zero = np.array( glob.glob(path_mframe+obs+'/zero/*zero.fits') )
				indx_closet = np.where(deltime == np.min(deltime))
				tmpflat = path_data+'/'+os.path.basename(np.asscalar(pastflat[indx_closet]))
				cpcom = 'cp {} {}'.format(np.asscalar(pastflat[indx_closet]), tmpflat)
				print(cpcom)	; os.system(cpcom)
				mflat = CCDData.read(tmpflat, hdu=0)#, unit='adu')
			calib.calibration(images, mzero, mflat, filte)
	#	DARK PROCESS
	elif dark_process == True:
		for i, filte in enumerate(objfilterlist):
			print('PRE PROCESS FOR {} FILTER OBJECT\t[{}/{}]'.format(filte, i+1, len(objfilterlist)))
			if filte in flatfilterlist:
				mflat = calib.master_flat(images, mzero, filte, mdark=mdark, fig=True)
			else:
				#	IF THERE IS NO FLAT FRAMES, BORROW FROM CLOSEST OTHER DATE
				print('\nNO {} FLAT FRAMES\n'.format(filte))
				pastflat = np.array(glob.glob('{}/{}/flat/*n{}.fits'.format(path_mframe, obs, filte)))
				#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
				deltime = []
				for date in pastflat:
					flatmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
					deltime.append(np.abs(obstime.mjd-flatmjd))
				# mframe_zero = np.array( glob.glob(path_mframe+obs+'/zero/*zero.fits') )
				indx_closet = np.where(deltime == np.min(deltime))
				tmpflat = path_data+'/'+os.path.basename(np.asscalar(pastflat[indx_closet]))
				cpcom = 'cp {} {}'.format(np.asscalar(pastflat[indx_closet]), tmpflat)
				print(cpcom)	; os.system(cpcom)
				mflat = CCDData.read(tmpflat, hdu=0)#, unit='adu')
			calib.calibration(images, mzero, mflat, filte, mdark=mdark)
	#------------------------------------------------------------
	#	ASTROMETRY
	#------------------------------------------------------------
	for i, inim in enumerate(glob.glob(path_data+'/fz*.fits')):
		print('ASTROMETRY\t[{}/{}]'.format(i+1,len(glob.glob(path_data+'/fz*.fits'))))
		calib.astrometry(inim, obsinfo['pixscale'])
	os.system('rm '+path_data+'/*.axy '+path_data+'/*.corr '+path_data+'/*.xyls '+path_data+'/*.match '+path_data+'/*.rdls '+path_data+'/*.solved '+path_data+'/*.wcs ')
	print('ASOTROMETRY COMPLETE\n'+'='*60)
	#------------------------------------------------------------
	#	FILE NAME CHANGE
	#------------------------------------------------------------
	for inim in glob.glob(path_data+'/afz*.fits'): calib.fnamechange(inim, obs)
	os.system('chmod 777 {}'.format(path_data))
	os.system('chmod 777 {}/*'.format(path_data))

	#	Calib-*.fits TO SAVE PATH
	caliblist = glob.glob(path_data+'/Calib*.fits'); caliblist.sort()
	for inim in caliblist: calib.movecalib(inim, path_gal)
	f = open(path_data+'/object.txt', 'a')
	f.write('obs obj dateobs filter exptime\n')
	for inim in caliblist:
		img = os.path.basename(inim)
		part = img.split('-')
		line = '{} {} {} {} {}\n'.format(part[1], part[2], part[3]+'T'+part[4], part[5], part[6])
		print(line)
		f.write(line)
	f.close()
	#	DATA FOLDER TO SAVE PATH
	os.system('rm {}/afz*.fits {}/fz*.fits'.format(path_data, path_data))
	os.system('rm -rf {}/{}'.format(path_save, os.path.basename(path_data)))
	os.system('mv {} {}'.format(path_data, path_save))
	#	WRITE LOG
	f = open(path_log, 'a')
	f.write(path_raw+'/'+os.path.basename(path_data)+'\n')
	f.close()

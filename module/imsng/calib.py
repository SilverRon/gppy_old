#============================================================
#   CALIBTRATION ROUTINE FUNCTION FOR IMSNG TELESCOPES
#	LOAO, DOAO, SOAO, ...
#   18.12.04	UPDATED BY Gregory S.H. Paek
#	19.08.08	CREATED BY Gregory S.H. Paek	(FUNCTIONALIZATION)
#============================================================
#	FUNCTION
#============================================================
import os, glob
import numpy as np
from astropy import units as u

def folder_check(datalist, path_raw, path_factory):
	'''
	CHECK NEW OBSERVATION DATA FOLDER
	'''
	newlist = []
	#	IS THERE NEW FILE?
	for data in glob.glob(path_raw+'/*'):
		if data not in datalist:
			print('NEW DATA\t: {}'.format(data))
			newlist.append(data)
	#	MOVE NEW FILE TO FACTORY FOLDER
	for data in newlist:
		#	PREVENT OVERWRITE
		if path_factory+'/'+os.path.basename(data) in glob.glob(path_factory+'/*'):
			rmcom = 'rm -rf '+path_factory+'/'+os.path.basename(data)
			print(rmcom) ; os.system(rmcom)
		cpcom = 'cp -r '+data+' '+path_factory
		print(cpcom) ; os.system(cpcom)
	#	CHANGE TO WORKING PATH
	for i in range(len(newlist)):
		newlist[i] = path_factory+'/'+os.path.basename(newlist[i])
	return newlist
#------------------------------------------------------------
def changehdr(inim, where, what):
	from astropy.io import fits
	data, hdr = fits.getdata(inim, header = True)
	hdr[where] = what
	fits.writeto(inim, data, hdr, clobber=True)
#------------------------------------------------------------
def correcthdr_routine(path_data, hdrtbl):
	'''
	1.	NGC337		-> NGC0337
		ngc0337		-> NGC0337
	'''
	from astropy.io import fits
	from astropy.time import Time
	comment = '-'*60+'\n' \
			+ 'CHANGE TO UNIVERSAL HEADER SERIES ...\n' \
			+ '-'*60+'\n'
	print(comment)
	objfilterlist = []
	objexptimelist = []
	flatfilterlist = []
	darkexptimelist = []
	for inim in glob.glob(path_data+'/*.fits'):
		data, hdr   = fits.getdata(inim, header=True)
		#	CHECK IMAGETYP HDR
		if hdr['IMAGETYP'] == 'Light': hdr['IMAGETYP'] = 'object'
		if hdr['IMAGETYP'] == 'zero' : hdr['IMAGETYP'] = 'Bias'
		if hdr['IMAGETYP'] == 'Dark' : hdr['IMAGETYP'] = 'dark'
		#	CHECK FILTER HDR
		if ((hdr['FILTER'] == 'U101') | (hdr['FILTER'] == 1)): hdr['FILTER'] = 'U'
		if ((hdr['FILTER'] == 'B102') | (hdr['FILTER'] == 2)): hdr['FILTER'] = 'B'
		if ((hdr['FILTER'] == 'V103') | (hdr['FILTER'] == 3)): hdr['FILTER'] = 'V'
		if ((hdr['FILTER'] == 'R104') | (hdr['FILTER'] == 4)): hdr['FILTER'] = 'R'
		if ((hdr['FILTER'] == 'I105') | (hdr['FILTER'] == 5)): hdr['FILTER'] = 'I'
		if (hdr['FILTER'] == '0') | (hdr['FILTER'] == 0):
			hdr['FILTER'] = 'R'
			print('PLEASE CHECK SOAO FILTER HDR (FILTER:0?)')
		#	CHECK DATE-OBS
		if 'T' not in hdr['DATE-OBS']: hdr['DATE-OBS'] = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
		t = Time(hdr['DATE-OBS'], format='isot')
		hdr['JD'] = t.jd
		hdr['MJD'] = t.mjd
		#	CHECK OBJECT HDR
		if 'ngc' in hdr['OBJECT']:
			while len(hdr['OBJECT'])<7:
				head = hdr['OBJECT'][0:3]
				tail = hdr['OBJECT'][3: ]
				tail = '0'+tail
				hdr['OBJECT'] = head+tail
		hdr['OBJECT'] = hdr['OBJECT'].upper()
		#	CHECK USER SETTING HDR
		for i in range(len(hdrtbl)):
			key = hdrtbl['key'][i]
			val = hdrtbl['val'][i]
			newval = hdrtbl['newval'][i]
			if hdr[key] == val:
				# print(hdr[key], key, newval)
				hdr[key] = newval
		fits.writeto(inim, data, hdr, overwrite=True)
		if hdr['IMAGETYP'] == 'object':
			objfilterlist.append(hdr['FILTER'])
			objexptimelist.append(hdr['EXPTIME'])
		if hdr['IMAGETYP'].upper() == 'FLAT':
			flatfilterlist.append(hdr['FILTER'])
		if hdr['IMAGETYP'].upper() == 'DARK':
			darkexptimelist.append(hdr['EXPTIME'])
	#	OBJECT FILTER, FLAT FILTER
	objfilterlist = list(set(objfilterlist));		objfilterlist.sort()
	objexptimelist = list(set(objexptimelist));		objexptimelist.sort()
	flatfilterlist = list(set(flatfilterlist));		flatfilterlist.sort()
	darkexptimelist = list(set(darkexptimelist));	darkexptimelist.sort()
	return objfilterlist, objexptimelist, flatfilterlist, darkexptimelist, t
#------------------------------------------------------------
def isot_to_mjd(time):      # 20181026 to 2018-10-26T00:00:00:000 to MJD form
	from astropy.time import Time
	yr  = time[0:4]         # year
	mo  = time[4:6]         # month
	da  = time[6:8]         # day
	isot = yr+'-'+mo+'-'+da+'T00:00:00.000'         #	ignore hour:min:sec
	t = Time(isot, format='isot', scale='utc')		#	transform to MJD
	return t.mjd
#------------------------------------------------------------
def call_images(path_data):
	import ccdproc
	from astropy.io import ascii
	images = ccdproc.ImageFileCollection(path_data, keywords='*')
	#	DATA SUMMARY
	ascii.write(images.summary, path_data+'/summary.txt', format='fixed_width_two_line', overwrite=True)
	return images
#------------------------------------------------------------
def getobsinfo(obs, obstbl):
	import numpy as np
	from astropy import units as u
	indx_obs = np.where(obs == obstbl['obs'])
	#	CCD INFORMATION
	gain = np.asscalar(obstbl['gain'][indx_obs]) * u.electron / u.adu
	rdnoise = np.asscalar(obstbl['RDnoise'][indx_obs]) * u.electron
	pixscale = np.asscalar(obstbl['pixelscale'][indx_obs])
	# dark        = float(obstbl['dark'][indx_obs][0])
	obsinfo = dict(	gain=gain,
					rdnoise=rdnoise,
					pixscale=pixscale)
	return obsinfo
#------------------------------------------------------------
#	CCD PROCESSING
#------------------------------------------------------------
def master_zero(images, fig=False):
	import ccdproc
	import os
	from astropy.nddata import CCDData
	import matplotlib.pyplot as plt
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER ZERO\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	#	HEADER FOR MASTER FRAME
	zerolist   = []
	for hdu, fname in images.hdus(imagetyp='Bias', return_fname=True):
		meta = hdu.header
		zerolist.append(ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu"))
	#	HEADER FOR MASTER FRAME
	n = 0
	for hdu, fname in images.hdus(imagetyp='Bias', return_fname=True):	
		n += 1
		newmeta = meta
		newmeta['FILENAME'] = 'zero.fits'
		newmeta['COMB{}'.format(n)] = fname
	print('{} ZERO IMAGES WILL BE COMBINED.'.format(len(zerolist)))
	zeros = ccdproc.Combiner(zerolist)
	mzero = zeros.median_combine()
	mzero.header  = newmeta
	#	SAVE MASTER FRAME
	if '{}/zero.fits'.format(path_data) in glob.glob(path_data+'/*'):
		os.system('rm {}/zero.fits'.format(path_data))
	mzero.write('{}/zero.fits'.format(path_data))
	if fig != False:
		imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
		zero_min, zero_max, zero_mean, zero_std = imstats(np.asarray(mzero))
		plt.figure(figsize=(15, 15))
		plt.imshow(mzero, vmax=zero_mean + 4*zero_std, vmin=zero_mean - 4*zero_std)
		plt.savefig(path_data+'/zero.png')
	return mzero
#------------------------------------------------------------
def master_dark(images, mzero, exptime, fig=False):
	import ccdproc
	import os
	from astropy.nddata import CCDData
	import matplotlib.pyplot as plt
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER DARK\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	mdark_name = 'dark-{}.fits'.format(int(exptime))
	#	HEADER FOR MASTER FRAME
	zdarklist   = []
	for hdu, fname in images.hdus(imagetyp='dark', exptime=exptime, return_fname=True):
		meta = hdu.header
		dark = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
		zdark = ccdproc.subtract_bias(dark, mzero)
		zdark.meta['SUBBIAS'] = mzero.meta['FILENAME']
		zdarklist.append(ccdproc.CCDData(data=zdark.data, meta=meta, unit="adu"))
	#	HEADER FOR MASTER FRAME
	n = 0
	for hdu, fname in images.hdus(imagetyp='dark', return_fname=True):	
		n += 1
		newmeta = meta
		newmeta['FILENAME'] = mdark_name
		newmeta['COMB{}'.format(n)] = fname
	# print('{} DARK({} sec) IMAGES WILL BE COMBINED.'.format(len(zdarklist)), exptime)
	darks = ccdproc.Combiner(zdarklist)
	mdark = darks.median_combine()
	mdark.header  = newmeta
	#	SAVE MASTER FRAME
	if '{}/{}'.format(path_data, mdark_name) in glob.glob(path_data+'/*'):
		os.system('rm {}/{}'.format(path_data, mdark_name))
	mdark.write('{}/{}'.format(path_data, mdark_name))
	if fig != False:
		imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
		zero_min, zero_max, zero_mean, zero_std = imstats(np.asarray(mzero))
		plt.figure(figsize=(15, 15))
		plt.imshow(mzero, vmax=zero_mean + 4*zero_std, vmin=zero_mean - 4*zero_std)
		plt.savefig(path_data+'/zero.png')
	return mdark
#------------------------------------------------------------
def master_flat(images, mzero, filte, mdark=False, fig=False):
	import os
	import ccdproc
	from astropy.nddata import CCDData
	import matplotlib.pyplot as plt
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER FLAT\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	if mdark == False:
		#	FLAT - MZERO = zFLAT
		print('SUBTRACT MASTER ZERO FROM {} FLAT'.format(filte))
		outname = 'n'+filte+'.fits'
		zflatlist = []
		for hdu, fname in images.hdus(imagetyp='flat', filter=filte, return_fname=True):
			meta = hdu.header
			flat = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
			zflat = ccdproc.subtract_bias(flat, mzero)
			zflat.meta['SUBBIAS'] = mzero.meta['FILENAME']
			zflatlist.append(ccdproc.CCDData(data=zflat.data, meta=meta, unit="adu"))
	else:
		#	FLAT - MZERO = zFLAT
		print('SUBTRACT MASTER ZERO & DARK FROM {} FLAT'.format(filte))
		outname = 'n'+filte+'.fits'
		zflatlist = []
		for hdu, fname in images.hdus(imagetyp='flat', filter=filte, return_fname=True):
			meta = hdu.header
			flat = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
			zflat = ccdproc.subtract_bias(flat, mzero)
			dzflat = ccdproc.subtract_dark(	ccd=zflat, master=mdark,
											# data_exposure=zflat.meta['EXPTIME'],
											exposure_time='EXPTIME',
											exposure_unit=u.second,
											# dark_exposure=mdark.meta['EXPTIME'],
											scale=True)
			meta['SUBBIAS'] = mzero.meta['FILENAME']
			meta['SUBDARK'] = mdark.meta['FILENAME']
			zflatlist.append(ccdproc.CCDData(data=dzflat.data, meta=meta, unit="adu"))
	#	HEADER FOR MASTER FRAME
	newmeta = meta
	newmeta['FILENAME'] = outname
	n = 0
	for hdu, fname in images.hdus(imagetyp='flat', filter=filte, return_fname=True):
		n += 1
		newmeta['COMB{}'.format(n)] = fname
	#	zFLATs -> MASTER FRAME
	print('{} {} FLAT IMAGES WILL BE COMBINED.'.format(len(zflatlist), filte))
	flat_combiner = ccdproc.Combiner(zflatlist)
	flat_combiner.minmax_clipping()
	def scaling_func(arr): return 1/np.ma.median(arr)
	flat_combiner.scaling = scaling_func
	mflat = flat_combiner.median_combine()
	#	SAVE MASTER FRAME
	mflat.header = newmeta
	if '{}/{}'.format(path_data, outname) in glob.glob(path_data+'/*'):
		os.system('rm {}/{}'.format(path_data, outname))
	mflat.write(path_data+'/'+outname)
	# mflat_electron = ccdproc.gain_correct(mflat, gain=gain)
	if fig != False:
		imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
		f_min, f_max, f_mean, f_std = imstats(np.asarray(mflat))
		plt.figure(figsize=(15, 15))
		plt.imshow(mflat, vmin=f_mean-5*f_std, vmax=f_mean+5*f_std)
		plt.savefig(path_data+'/'+outname[:-5]+'.png')

	return mflat
#------------------------------------------------------------
def calibration(images, mzero, mflat, filte, mdark=False):
	import ccdproc
	comment     = '-'*60+'\n' \
				+ 'OBJECT CORRECTION\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	objlist = []
	for hdu, fname in images.hdus(imagetyp='object', filter=filte, return_fname=True):
		meta = hdu.header
		objlist.append(meta['object'])
		obj = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
		if mdark == False:
			#	ZERO CORRECTION
			zobj = ccdproc.subtract_bias(obj, mzero)
			meta['SUBBIAS'] = mzero.meta['FILENAME']
		else:
			#	ZERO & DARK CORRECTION
			zobj = ccdproc.subtract_bias(obj, mzero)
			meta['SUBBIAS'] = mzero.meta['FILENAME']
			zobj = ccdproc.subtract_dark(	ccd=zobj, master=mdark,
											# dark_exposure=mdark.meta['EXPTIME'],
											# data_exposure=dzobj.meta['EXPTIME'],
											exposure_time='EXPTIME',
											exposure_unit=u.second,
											scale=True)
			meta['SUBDARK'] = mdark.meta['FILENAME']
		#	FLAT CORRECTION
		fzobj = ccdproc.flat_correct(zobj, mflat)
		meta['DIVFLAT'] = mflat.meta['FILENAME']
		fzobj.header = meta
		# fzobj.write(path_data+'/fz'+fname)
		fzobj.write('{}/fz{}'.format(path_data, fname), overwrite=True)
#------------------------------------------------------------
'''
def astrometry(path_data, pixscale):
	import os
	order = 0
	imlist = glob.glob(path_data+'/fz*.fits')
	upscl = str(pixscale + pixscale*0.05)
	loscl = str(pixscale - pixscale*0.05)
	comment = '='*60+'\n' \
			+ 'ASTROMETRY PROCESS START\t'+str(len(imlist))+' IMAGES'
	print(comment)
	for inim in imlist:
		order += 1
		inim = inim.split('/')[-1]
		outname = path_data+'/a'+inim
		com     ='solve-field '+inim \
				+' --scale-unit arcsecperpix --scale-low '+loscl+' --scale-high '+upscl \
				+' --no-plots --new-fits '+outname+' --overwrite --use-sextractor --cpulimit 1800\n'
		os.system(com)
		comment = '-'*60+'\n' \
				+ 'ASTROMETRY PROCESS\t'+str(round(order*100/len(imlist)))+'%\t'\
				+ '['+str(order)+'/'+str(len(imlist))+']'+'\n'
		print(comment)
	os.system('rm '+path_data+'/*.axy '+path_data+'/*.corr '+path_data+'/*.xyls '+path_data+'/*.match '+path_data+'/*.rdls '+path_data+'/*.solved '+path_data+'/*.wcs ')
	print('COMPLETE\n'+'='*60)
'''
#------------------------------------------------------------

def astrometry(inim, pixscale):
	import os
	upscl = str(pixscale + pixscale*0.10)
	loscl = str(pixscale - pixscale*0.10)
	outname = os.path.dirname(inim)+'/a'+os.path.basename(inim)
	com     ='solve-field '+inim \
			+' --scale-unit arcsecperpix --scale-low '+loscl+' --scale-high '+upscl \
			+' --no-plots --new-fits '+outname+' --overwrite --use-sextractor --cpulimit 600\n'
	print(com); os.system(com)

#------------------------------------------------------------
def fnamechange(inim, obs):
	import os
	from astropy.io import fits
	hdr = fits.getheader(inim)
	dateobs = hdr['DATE-OBS']
	datestr = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
	timestr = dateobs[11:13]+dateobs[14:16]+dateobs[17:19]
	objname = hdr['OBJECT']
	objname = objname.upper()	
	filname = hdr['FILTER']  # str(hdr['FILTER'])
	exptime = str(int(hdr['EXPTIME']))
	newname = 'Calib-{}-{}-{}-{}-{}-{}.fits'.format(obs, objname, datestr, timestr, filname, exptime)
	cpcom = 'cp {} {}/{}'.format(inim, os.path.dirname(inim), newname)
	print(cpcom); os.system(cpcom)
#------------------------------------------------------------
def movecalib(inim, path_gal):
	import os, glob
	img = os.path.basename(inim)
	part = img.split('-')
	obs = part[1]
	obj = part[2]
	filte = part[5]
	#	MAKE SAVE PATH
	path_goal = path_gal
	for folder in [obj, obs, filte]:
		path_goal = path_goal+'/'+folder
		mkcom = 'mkdir {}'.format(path_goal)
		print(mkcom); os.system(mkcom)
	mvcom = 'mv {} {}'.format(inim, path_goal)
	print(mvcom); os.system(mvcom)

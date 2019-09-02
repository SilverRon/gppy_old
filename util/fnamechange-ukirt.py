def fnamechange_ukirt_raw(inim, obj):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= hdr['TELESCOP']
	# obj		= hdr['MSBTITLE'].split(':')[0]
	obj = 'GRB190829A'
	date	= hdr['UTDATE'][2:]
	ut_pre	= hdr['DATE-OBS'].split('T')[1]
	ut		= ut_pre[0:2]+ut_pre[3:5]+ut_pre[6:8]
	filte	= hdr['FILTER']
	exptime	= str(int(hdr['EXP_TIME']))
	parts	= ['Calib', obs, obj, date, ut, filte, exptime+'.fits']
	newname	= '-'.join(parts)
	#	COPY RAW -> NEW NAME
	com		= 'cp '+inim+' '+newname
	print(com)
	os.system(com)

import os, glob
obj = input('OBJECT NAME :\t')
for inim in glob.glob('*.fit'):
	fnamechange_ukirt_raw(inim, obj)
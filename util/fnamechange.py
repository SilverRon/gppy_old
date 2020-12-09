#-------------------------------------------------------------------------#
def fnamechange_ukirt_raw(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= hdr['TELESCOP']
	obj		= hdr['MSBTITLE'].split(':')[0]
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
#-------------------------------------------------------------------------#
def fnamechange_ukirt(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= hdr['TELESCOP']
	# obj		= hdr['MSBTITLE'].split(':')[0]
	obj		= hdr['OBJECT']
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
#-------------------------------------------------------------------------#
def fnamechange_boao(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= hdr['OBSERVAT']
	obj		= hdr['OBJECT']
	date	= ''.join(hdr['DATE-OBS'][2:].split('-'))
	ut		= ''.join(hdr['TIME-OBS'].split(':'))[:6]
	filte	= hdr['FILTER']
	exptime	= str(int(hdr['EXPTIME']))
	parts	= ['Calib', obs, obj, date, ut, filte, exptime+'.fits']
	newname	= '-'.join(parts)
	#	COPY RAW -> NEW NAME
	com		= 'cp '+inim+' '+newname
	print(com)
	os.system(com)
#-------------------------------------------------------------------------#
def fnamechange_cca250(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= 'CCA250'
	obj		= hdr['OBJECT']

	date_pre= (hdr['DATE-OBS'].split('T')[0]).split('-')
	date_pre[0]	= date_pre[0][2:]
	date	= ''.join(date_pre)
	ut		= ''.join(((hdr['DATE-OBS'].split('T')[1]).split('-')[0]).split(':'))

	filte	= hdr['FILTER']
	exptime	= str(int(hdr['EXPTIME']))
	parts	= ['Calib', obs, obj, date, ut, filte, exptime+'.fits']
	newname	= '-'.join(parts)
	#	COPY RAW -> NEW NAME
	com		= 'cp '+inim+' '+newname
	print(com)
	os.system(com)
#-------------------------------------------------------------------------#
def fnamechange_kmtnet(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= 'KMTNET'
	date_pre= (hdr['DATE-OBS'].split('T')[0]).split('-')
	date_pre[0]	= date_pre[0][2:]
	date	= ''.join(date_pre)
	ut		= ''.join(((hdr['DATE-OBS'].split('T')[1]).split('-')[0]).split(':'))

	filte	= hdr['FILTER']
	exptime	= str(int(hdr['EXPTIME']))
	obj		= hdr['OBJECT']
	parts	= ['Calib', obs, obj, '20'+date, ut, filte, exptime+'.fits']
	newname	= '-'.join(parts)
	#	COPY RAW -> NEW NAME
	com		= 'cp '+inim+' '+newname
	print(com)
	os.system(com)
#-------------------------------------------------------------------------#
def fnamechange_soao(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= 'SOAO'
	date_pre= (hdr['DATE-OBS'].split('T')[0]).split('-')
	date_pre[0]	= date_pre[0][2:]
	date	= ''.join(date_pre)
	ut		= ''.join(((hdr['DATE-OBS'].split('T')[1]).split('-')[0]).split(':'))

	filte	= hdr['FILTER']
	exptime	= str(int(hdr['EXPTIME']))
	obj		= hdr['OBJECT']
	parts	= ['Calib', obs, obj, date, ut, filte, exptime+'.fits']
	newname	= '-'.join(parts)
	#	COPY RAW -> NEW NAME
	com		= 'cp '+inim+' '+newname
	print(com)
	os.system(com)
#-------------------------------------------------------------------------#
def fnamechange_doaofli(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= 'DOAO_FLI'
	date_pre= (hdr['DATE-OBS'].split('T')[0]).split('-')
	date_pre[0]	= date_pre[0][2:]
	date	= ''.join(date_pre)
	ut		= ''.join(((hdr['DATE-OBS'].split('T')[1]).split('-')[0]).split(':'))

	filte	= hdr['FILTER']
	exptime	= str(int(hdr['EXPTIME']))
	obj		= hdr['OBJECT']
	parts	= ['Calib', obs, obj, date, ut, filte, exptime+'.fits']
	newname	= '-'.join(parts)
	#	COPY RAW -> NEW NAME
	com		= 'cp '+inim+' '+newname
	print(com)
	os.system(com)
'''
#	UKRIT 'w20XX*.fit' DATA
import os, glob

imlist		= glob.glob('w2*st.fit')
for inim in imlist	: fnamechange_ukirt(inim)



'''

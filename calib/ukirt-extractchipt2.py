#	UKIRT RAW IMAGE (CHIP 2) TO Calib-*.fits FORM
#	2019.08.23	CREATED BY Gregory S.H. Paek
#============================================================
import os, glob
#============================================================
#	FUNCTION
#============================================================
def extractcam(inim, number):
	'''
	Filename: Calib-UKIRT-GRB190114C-190116-050540-K-10.fits
	No.    Name      Ver    Type      Cards   Dimensions   Format
	0  PRIMARY       1 PrimaryHDU     216   ()      
	1                1 CompImageHDU    142   (4167, 4167)   int32   
	2                1 CompImageHDU    138   (4167, 4167)   float64   
	3                1 CompImageHDU    138   (4167, 4167)   int32   
	4                1 CompImageHDU    138   (4167, 4167)   int32  
	'''
	import os
	from astropy.io import fits
	# newname	= 'Split-'+inim
	# newname = 'Chip{}-{}'.format(number, inim)
	newname = fnamechange_ukirt(inim)
	hdu		= fits.open(inim)
	hdu[0].header.remove('CC_PRES')
	cam		= hdu[number]
	fits.writeto(newname, cam.data, hdu[0].header+hdu[2].header)
#------------------------------------------------------------
def fnamechange_ukirt(inim):
	from astropy.io import fits
	import os
	hdr		= fits.getheader(inim)
	#	GATHER PARTS
	obs		= hdr['TELESCOP']
	#obj		= hdr['MSBTITLE'].split(':')[0]
	obj		= hdr['OBJECT']
	date	= hdr['UTDATE'][2:]
	ut_pre	= hdr['DATE-OBS'].split('T')[1]
	ut		= ut_pre[0:2]+ut_pre[3:5]+ut_pre[6:8]
	filte	= hdr['FILTER']
	exptime	= str(int(hdr['EXP_TIME']))
	parts	= ['Calib', obs, obj, date, ut, filte, exptime+'.fits']
	newname	= '-'.join(parts)
	#	COPY RAW -> NEW NAME
	# com		= 'cp '+inim+' '+newname
	# print(com)
	# os.system(com)
	return newname
#============================================================
imlist	= glob.glob('w*st.fit'); imlist.sort()
for inim in imlist: extractcam(inim, 2)
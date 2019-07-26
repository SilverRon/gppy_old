#	EXTRACT SPECIFIC CHIP DATA FROM MULTI-HEADER UKIRT IMAGE
#	CREATED	2019.07.24	Gregory S.H. Paek
#============================================================
import os, glob
from astropy.io import fits
#------------------------------------------------------------
def extract_chip(inim, numb):
	newim = 'Chip{0}-{1}'.format(numb, inim) 
	hdul = fits.open(inim)
	'''
	hdul.info()
	Filename: w20190425_02578_sf_st.fit
	No.    Name      Ver    Type      Cards   Dimensions   Format
	0  PRIMARY       1 PrimaryHDU     216   ()      
	1                1 CompImageHDU    118   (4135, 4135)   int32   
	2                1 CompImageHDU    118   (4135, 4135)   int32   
	3                1 CompImageHDU    120   (4135, 4135)   int32   
	4                1 CompImageHDU    118   (4135, 4135)   int32   
	'''
	data = hdul[numb].data
	# hdr = hdul['Primary'].header
	# hdr.set('CC_PRES', 0)
	keys = ['TELESCOP', 'INSTRUME', 'UTDATE', 'MSBTITLE', 'OBJECT', 'OBSTYPE', 'DATE-OBS', 'DATE-END', 'MJD-OBS', 'EXP_TIME', 'FILTER']
	hdr = hdul[numb].header
	for key in keys:
		hdr.set(key, hdul['Primary'].header[key])
	fits.writeto(newim, data, hdr)
	#------------------------------------------------------------
	utdate = hdr['DATE-OBS']
	partm = utdate.split('T')
	partm0, partm1 = partm[0], partm[1]
	partmb = partm0.split('-')
	partmb[0] = partmb[0][2:]
	partmc = partm1.split(':')
	for i in range(len(partmc)):
		partmc[i] = partmc[i][:2]
	utdate_name = ''.join(partmb)+'-'+''.join(partmc)
	#------------------------------------------------------------
	newimpart = ['Calib', hdr['TELESCOP'], hdr['OBJECT'], utdate_name, hdr['FILTER'], str(int(hdr['EXP_TIME']))+'.fits']
	newname = '-'.join(newimpart)
	os.system('mv {0} {1}'.format(newim, newname))
#------------------------------------------------------------
numb = 2
# inim = 'w20190425_02578_sf_st.fit'
for inim in glob.glob('w*.fit'): extract_chip(inim, numb)

'''
for inim in glob.glob('Calib-UKIRT-2MASS+1*-*-19*.fits'): 
	newim = inim[:26]+'_'+inim[27:] 
	os.system('cp {0} {1}'.format(inim, newim)) 
'''
#	EXTRACT ALL CHIP DATA FROM MULTI-HEADER UKIRT IMAGE
#	i.e.) S190814bv TILING OBSERVATION
#	2019.07.24	CREATED	BY	Gregory S.H. Paek
#	2019.09.30	MODIFIED BY	Gregory S.H. Paek
#============================================================
import os, glob
from astropy.io import fits
#------------------------------------------------------------
def extract_chip(inim, numb, objkey=False):
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
	partmb[0] = partmb[0]
	partmc = partm1.split(':')
	for i in range(len(partmc)):
		partmc[i] = partmc[i][:2]
	utdate_name = ''.join(partmb)+'-'+''.join(partmc)
	#------------------------------------------------------------
	if objkey == False:
		obj = hdr['OBJECT']
	newimpart = ['Calib', hdr['TELESCOP'], objkey, utdate_name, hdr['FILTER'], str(int(hdr['EXP_TIME']))+'.fits']
	newname = '-'.join(newimpart)
	os.system('mv {0} {1}'.format(newim, newname))
#------------------------------------------------------------
# numb = 2
# objkey = input('OBJECT \t:')
namedict = dict(a=1, b=2, c=3, d=4)
for inim in glob.glob('w*.fit'):
	hdul = fits.open(inim)
	hdr = hdul[0].header
	objkey = hdr['OBJECT']
	objkey = objkey.replace('-', '_')
	for key in namedict:
		numb = namedict[key]
		nobjkey = objkey+key
		print(inim, numb, nobjkey)
		extract_chip(inim, numb, objkey=nobjkey)

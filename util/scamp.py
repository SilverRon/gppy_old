#	2016.06.22	Created by	Changsu Choi
#	2020.05.12	Modified by	Gregory S.H. Paek
#------------------------------------------------------------
#	Former Description
#------------------------------------------------------------
## astrometry with scamp, sextractror, cdsclient, scamp, astropy.io.fits
## using ra,dec from header
## sextractor - scamp - header editing
## reference catalog query from CDSCLIENT, then reference catalog in SCAMP  
## if you need reference catalog from reference image's astrometry with scamp : option "save reference = yes" 
##
## USE : python scamprun-radec.py inim.fits ra dec
## python ~/Desktop/9codes/automatic_1d_maidanak/scamprun-radec.py inim 20:24:28.104 -43:39:12.71
#============================================================
#	Module
#------------------------------------------------------------
from astropy.io import ascii
import numpy as np
import os,sys, glob
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
# from pyraf import iraf
#============================================================
inim = 'a045289.kk-1.fits'
incat = 'a045289.kk-1.cat'

data, hdr = fits.getdata('a045289.kk-1.fits', header=True)

cd11, cd12, cd21, cd22 = hdr['CD1_1'], hdr['CD1_2'], hdr['CD2_1'], hdr['CD2_2']

os.system('ls *.fits')
inim	= input('INPUT IMAGE\t: ')
hdr		= fits.getheader(inim)
ra		= str(hdr['RABASE'])
dec		= str(hdr['DECBASE'])

# configuration files
astromsex   = '/home/sonic/Research/yourpy/config/targetphot.sex'
astromscamp = './default.scamp'
astromparam = '/home/sonic/Research/yourpy/config/targetphot.param'

# option setup

	# position define
	# (xpix,ypix),  (ra,dec)
	# cd matrix cd11,cd12,cd21,cd22
xpix		= hdr['NAXIS1']/2.
ypix		= hdr['NAXIS2']/2.

#ra='12 12 12.0'
#dec='+15 15 15.0'

#ra=hdr['RA']
#dec=hdr['DEC']
rad			= coord.Angle(ra, unit=u.hour)	
radd		= rad.degree
decd		= coord.Angle(dec, unit=u.deg)
decdd		= decd.degree
	
# CD matrix for Maidanak SNUCAM, if you don't know cd matrix then try nova.astrometry.net
'''
cd11= -0.000127886720156 
cd12= -2.63489094687E-06
cd21= 2.60446707435E-06
cd22= -0.000128805470953 
'''
cd11   =   -7.43364634984E-05 #/ Transformation matrix
cd12   =   -1.98606556866E-06 #/ no comment
cd21   =   -2.00823634932E-06 #/ no comment
cd22   =    7.43498410235E-05 #/ no comment

'''
# sample from astrometry.net result
CRVAL1  =        100.716380035 / RA  of reference point
CRVAL2  =       -74.2443939374 / DEC of reference point
CRPIX1  =        782.111596065 / X reference pixel
CRPIX2  =        497.279253685 / Y reference pixel

CD1_1   =   -0.000176942459758 / Transformation matrix
CD1_2   =   -1.46521919518E-07 / no comment
CD2_1   =    3.30848303174E-08 / no comment
CD2_2   =   -0.000176863297194 / no comment
'''
def wcsreset(inim,wcschar):
	#wcschar='physical' or 'world'
	iraf.wcsreset.image=inim
	iraf.wcsreset.wcs=wcschar
	iraf.wcsreset.mode='h'
	iraf.wcsreset()

def hedit(inim,kword,value):
	iraf.hedit.images=inim
	iraf.hedit.fields=kword
	iraf.hedit.value=value
	iraf.hedit.add='yes'
	iraf.hedit.verify='no'
	iraf.hedit.update='yes'
	iraf.hedit.mode='h'
	iraf.hedit()

# sextractor
'''
CATALOG_NAME     test.cat
PARAMETERS_NAME  astrom.param
DETECT_MINAREA   3              # minimum number of pixels above threshold
DETECT_THRESH    3            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
FILTER_NAME      default.conv   # name of the file containing the filter
STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename
'''
threshold,minarea=10,5
sexcom='sex -c '+astromsex +' '+inim+' -PARAMETERS_NAME '+astromparam+' -DETECT_THRESH '+str(threshold) +' -DETECT_MINAREA '+str(minarea)+' -CATALOG_NAME astromtest.cat '+'-FILTER_NAME /home/lim9/Desktop/9codes/astrom.config/default.conv -STARNNW_NAME /home/lim9/Desktop/9codes/astrom.config/default.nnw'

# scamp
'''
ASTREF_CATALOG         FILE #2MASS         # NONE, FILE,USNO-A1, USNO-A2, USNO-B1,
                                       # GSC-1.3, GSC-2.2, UCAC-1, UCAC-2,
                                       # 2MASS, SDSS-R3, or SDSS-R5
ASTREFCAT_NAME         astrefcat.cat   # Local astrometric reference catalogs
MERGEDOUTCAT_NAME      scamp.cat       # Merged output catalog filename
CHECKPLOT_DEV          PNG             # NULL, XWIN, TK, PS, PSC, XFIG, PNG,
                                       # or JPEG
CHECKPLOT_TYPE         FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR
CHECKPLOT_NAME         fgroups,distort,astr_interror2d,astr_interror1d,astr_referror2d,astr_referror1d,astr_chi2,psphot_error 
'''

scampcom='scamp -c '+astromscamp+' astromtest.cat '+ '-ASTREF_CATALOG UCAC-3 -CHECKPLOT_DEV PNG'

# iraf setting

wcsreset(inim,'physical')
wcsreset(inim,'world')
hedit(inim,'WAT0_001','system=image')
hedit(inim,'WAT1_001', 'wtype=tan axtype=ra')
hedit(inim,'WAT2_001', 'wtype=tan axtype=dec')
hedit(inim,'RADECSYS', 'FK5')
hedit(inim,'EQUINOX', '2000.')
hedit(inim,'CTYPE1', 'RA---TAN')
hedit(inim,'CTYPE2', 'DEC--TAN')
hedit(inim,'CRVAL1', str(radd))
hedit(inim,'CRVAL2', str(decdd))
hedit(inim,'CRPIX1', str(xpix))
hedit(inim,'CRPIX2', str(ypix))
hedit(inim,'CD1_1', str(cd11))
hedit(inim,'CD1_2', str(cd12))
hedit(inim,'CD2_1', str(cd21))
hedit(inim,'CD2_2', str(cd22))

os.system(sexcom)
os.system(scampcom)

# new fits file as result 
# header editing

hdr1=hdr
#hdr1.fromTxtFile('astromtest.head')
hdr1.extend(fits.Header.fromtextfile('astromtest.head'), update=True, update_first=True)
#hdr1.fromtextfile('astromtest.head',update=True,update_first=True)
fits.writeto('a'+inim,fits.getdata(inim),hdr1)



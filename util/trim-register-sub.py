#============================================================
#	TRIM - GREGISTERING - HOTPANTS (SAME OBSERVATORY SCI & REF IMAGES)
#	Usage : 
#	OR JUST RUN AND WRITE INPUT & REF IMAGE
#	2020.02.14	GREGORY S.H. PAEK
#============================================================
import numpy as np
import os, sys, glob
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import alipy
#============================================================
#	FUNCTION
#============================================================
def trim(inim, position, size, outim='trim.fits'):
	# Load the image and the WCS
	hdu = fits.open(inim)[0]
	wcs = WCS(hdu.header)
	# Make the cutout, including the WCS
	cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
	# Put the cutout image in the FITS HDU
	hdu.data = cutout.data
	# Update the FITS header with the cutout WCS
	hdu.header.update(cutout.wcs.to_header())
	# Write the cutout to a new FITS file
	hdu.writeto(outim, overwrite=True)
#-------------------------------------------------------------
def gregistering(imlist, tmpim):
	identifications = alipy.ident.run(tmpim, imlist, visu=False)
	for id in identifications:	#	list of the same length as imlist.
		if id.ok == True:		#	i.e., if it worked
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
		else:
			print "%20s : no transformation found !" % (id.ukn.name)
	outputshape = alipy.align.shape(tmpim)
	for id in identifications:
		if id.ok == True:
			params_align	= dict(	filepath	= id.ukn.filepath,
									uknstarlist	= id.uknmatchstars,
									refstarlist	= id.refmatchstars,
									shape		= alipy.align.shape(tmpim),
									# outdir		= './{}'.format(outim),
									outdir		= './',
									makepng		= False)
			alipy.align.irafalign(**params_align)
#-------------------------------------------------------------
def hotpants(imlist, refim):
	for inim in imlist:
		# outfile = os.path.dirname(inim)+'/hd'+os.path.basename(inim)
		# convfile= os.path.dirname(inim)+'/hc'+os.path.basename(inim)
		outfile = 'hd'+inim
		convfile = 'hc'+inim
		#com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -v 0 -inim '+infile[n]+' -tmplim '+ref_img+' -outim '+outfile+' -oci '+convfile
		#com='hotpants -c t -n i -inim '+infile[n]+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
		com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -tl -100000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
		os.system(com)
#============================================================
#	INPUT
#============================================================
os.system('ls *.fits')
imlist = glob.glob(raw_input('Science Image\t: ')); imlist.sort()
tmpim = raw_input('Template Image\t: ')
length = float(raw_input('Trim Size [arcmin]\t: '))

tra, tdec = 44.5438520, -8.9577875				#	GRB 190829A
# tra, tdec = 44.863251, 31.385878				#	G0037111
#============================================================
position = SkyCoord(tra, tdec, frame='icrs', unit='deg')
# length = 5.0	#	[']
size = u.Quantity((length, length), u.arcmin)

for inim in imlist:
	inhdr = fits.getheader(inim)
	tmphdr = fits.getheader(tmpim)

	# try:
		# inseeing = inhdr['SEEING']
		# tmpseeing = tmphdr['SEEING']
	# except:
		# inseeing = 3.0
		# tmpseeing = 1.0
	
	inseeing = 3.0
	tmpseeing = 1.0

	trinim = 'tr'+inim
	trtmpim = 'tr'+tmpim
	#	STEP 1 TRIM
	trim(inim, position, size, trinim)
	trim(tmpim, position, size, trtmpim)

	#	STEP 2 REGISTER
	if inseeing >= tmpseeing:
		outim = '{}_ref2sci.fits'.format(trinim[:-5])
		gregistering([trtmpim], trinim)
		greim = trtmpim[:-5]+'_gregister.fits'
	elif inseeing < tmpseeing:
		outim = '{}_sci2ref.fits'.format(trinim[:-5])
		gregistering([trinim], trtmpim)
		greim = trinim[:-5]+'_gregister.fits'

	os.system('mv {} {}'.format(greim, outim))

	#	STEP 3 SUBTRACTION
	hotpants([trinim], outim)